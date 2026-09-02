// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// felitronics::analysis::MultiResSpectrumPane — property tests against the physics the pane claims
// (docs/ANALYZER-MULTIRES.md), not against golden numbers. Oracles: a direct DFT of the exact suffix a
// tier transforms; the analytic band power of white noise (4σ²B/fs); a full-scale sine's lobe; the Hann
// window's own energy distribution for the time-resolution tests. The seam tests are statistical with a
// fixed PRNG and ensemble averaging, and are checked at three different hops so the reading cannot
// depend on the UI rate. Adversarial cases come from the design consilium (codex + deepseek + Fable).

#include <felitronics_test.h>
#include <felitronics/analysis/MultiResSpectrumPane.h>

#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <memory>
#include <vector>

static std::atomic<long> g_allocs { 0 };
void* operator new      (std::size_t s) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s ? s : 1); }
void* operator new[]    (std::size_t s) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s ? s : 1); }
void  operator delete   (void* p) noexcept { std::free (p); }
void  operator delete[] (void* p) noexcept { std::free (p); }
void  operator delete   (void* p, std::size_t) noexcept { std::free (p); }
void  operator delete[] (void* p, std::size_t) noexcept { std::free (p); }

using namespace felitronics;
using felitronics::test::ok;
using felitronics::test::approx;
using felitronics::test::group;
using Pane = analysis::MultiResSpectrumPane;

namespace
{
constexpr double kPi = 3.14159265358979323846;
constexpr double fs  = 48000.0;
constexpr int    N   = Pane::kMaxSize;          // 16384, the default frame
constexpr double kBandFactor = 0.028882136;     // 2^(1/48) − 2^(−1/48): B(f)/f for a 1/24-oct band

struct Rng                                       // xorshift64*, fixed seed → deterministic ensembles
{
    std::uint64_t s = 0x9E3779B97F4A7C15ull;
    double uni() noexcept                        // uniform in [−1, 1): variance 1/3
    {
        s ^= s >> 12; s ^= s << 25; s ^= s >> 27;
        const std::uint64_t r = s * 0x2545F4914F6CDD1Dull;
        return (double) (r >> 11) * (2.0 / 9007199254740992.0) - 1.0;
    }
};

struct Pink                                      // Kellet's 7-pole pink filter (±0.05 dB above ~10 Hz at 44.1 k)
{
    double b0 = 0, b1 = 0, b2 = 0, b3 = 0, b4 = 0, b5 = 0, b6 = 0;
    double operator() (double w) noexcept
    {
        b0 = 0.99886 * b0 + w * 0.0555179; b1 = 0.99332 * b1 + w * 0.0750759; b2 = 0.96900 * b2 + w * 0.1538520;
        b3 = 0.86650 * b3 + w * 0.3104856; b4 = 0.55000 * b4 + w * 0.5329522; b5 = -0.7616 * b5 - w * 0.0168980;
        const double p = b0 + b1 + b2 + b3 + b4 + b5 + b6 + w * 0.5362;
        b6 = w * 0.115926;
        return p * 0.11;                         // gain into a sane range; a level test compares slopes only
    }
};

void fillSine (float* f, int n, double hz, double amp, double phase = 0.0)
{
    for (int i = 0; i < n; ++i) f[i] = (float) (amp * std::sin (2.0 * kPi * hz * (double) i / fs + phase));
}
void fillConst (float* f, int n, float v) { for (int i = 0; i < n; ++i) f[i] = v; }

// A full-scale sine at hz, raw (smoothCoeff 1 → the reading IS the frame). Returns readDb (hz).
double rawSineDb (Pane& p, double hz)
{
    p.reset(); p.smoothCoeff = 1.0f;
    fillSine (p.frameInput(), N, hz, 1.0);
    p.ingest (14);
    return p.readDb (hz, fs);
}

double binHzOf (const Pane& p, int k) { return fs / (double) (1 << p.tierOrder (k)); }
double binsPerBandAt (const Pane& p, int k, double f) { return f * kBandFactor / binHzOf (p, k); }

// Slope in dB/oct of a set of (f, dB) points by least squares on log2 f.
double slopeDbPerOct (const std::vector<std::pair<double, double>>& pts)
{
    double sx = 0, sy = 0, sxx = 0, sxy = 0; const double n = (double) pts.size();
    for (auto& [f, v] : pts) { const double x = std::log2 (f); sx += x; sy += v; sxx += x * x; sxy += x * v; }
    return (n * sxy - sx * sy) / (n * sxx - sx * sx);
}

// Stream independent noise frames through the pane: each tick a fresh 16384-sample frame.
template <class Gen>
void streamFrames (Pane& p, Gen&& gen, int ticks)
{
    for (int t = 0; t < ticks; ++t)
    {
        float* in = p.frameInput();
        for (int i = 0; i < N; ++i) in[i] = (float) gen();
        p.ingest (14);
    }
}
}

int main()
{
    std::printf ("felitronics::analysis::MultiResSpectrumPane tests\n");
    auto pane = std::make_unique<Pane>();
    Pane& p = *pane;

    //==========================================================================================
    group ("tiers, seams, and the tier map");
    {
        ok (p.tierCount() == 3 && p.tierOrder (0) == 14 && p.tierOrder (1) == 12 && p.tierOrder (2) == 10, "default tiers 16384 / 4096 / 1024");
        ok (p.frameOrder() == 14 && p.frameSize() == 16384, "the frame to request is the longest tier");
        // seam_k = binsPerBand · binHz_k / (2^(o/2) − 2^(−o/2)), binsPerBand 2, 1/24 oct
        approx (p.seamHz (1, fs), 2.0 * (fs / 4096.0) / kBandFactor, 0.05, "seam 4096 = 811.5 Hz at 48 k");
        approx (p.seamHz (2, fs), 2.0 * (fs / 1024.0) / kBandFactor, 0.2,  "seam 1024 = 3246 Hz at 48 k");
        approx (p.seamHz (2, 96000.0), 2.0 * p.seamHz (2, fs), 1e-6, "seams double at 2× fs");
        ok (p.seamHz (0, fs) == 0.0, "the longest tier has no seam");
        int prev = 0; bool mono = true;
        for (double f = 10.0; f < 24000.0; f *= 1.01) { const int k = p.tierAt (f, fs); if (k < prev) mono = false; prev = k; }
        ok (mono, "tierAt is non-decreasing in f");
        ok (p.tierAt (100.0, fs) == 0 && p.tierAt (1000.0, fs) == 1 && p.tierAt (5000.0, fs) == 2, "tier regions at 48 k: <811 → 16384, <3246 → 4096, above → 1024");
        approx (p.enbwBins (0), 1.5 * 16384.0 / 16383.0, 1e-9, "ENBW of the symmetric Hann = 1.5·N/(N−1)");

        auto q = std::make_unique<Pane>();
        const int weird[] = { 3, 14, 14, 99, 12 };
        ok (q->setTiers (weird, 5) == 3 && q->tierOrder (0) == 14 && q->tierOrder (1) == 12 && q->tierOrder (2) == 8, "setTiers: clamp, dedup, sort longest-first");
        ok (q->setTiers (weird, 0) == 1 && q->tierOrder (0) == 14, "empty list → the single longest window");
        ok (q->setTiers (nullptr, 3) == 1 && q->tierOrder (0) == 14, "null list → the single longest window");
        const int six[] = { 9, 10, 11, 12, 13, 14 };
        ok (q->setTiers (six, 6) == 4 && q->tierOrder (3) == 11, "more orders than MaxTiers keeps the longest");
    }

    //==========================================================================================
    group ("oracle: each tier's bins equal a direct DFT of the exact suffix (offset, window, calibration)");
    {
        p.reset(); p.coverSamples = 0; p.smoothCoeff = 1.0f;
        Rng r; float* in = p.frameInput();
        for (int i = 0; i < N; ++i) in[i] = (float) r.uni();
        std::vector<float> frame (in, in + N);
        p.ingest (14);
        for (int k = 0; k < 3; ++k)
        {
            const int n = 1 << p.tierOrder (k);
            std::vector<double> w ((std::size_t) n); double sw2 = 0.0;
            for (int i = 0; i < n; ++i) { w[(std::size_t) i] = 0.5 - 0.5 * std::cos (2.0 * kPi * (double) i / (double) (n - 1)); sw2 += w[(std::size_t) i] * w[(std::size_t) i]; }
            const int probe[] = { 0, 3, 17, n / 8, n / 2 - 1, n / 2 };
            double worst = 0.0;
            for (int b : probe)
            {
                double re = 0.0, im = 0.0;
                for (int i = 0; i < n; ++i)
                {
                    const double x = (double) frame[(std::size_t) (N - n + i)] * w[(std::size_t) i];   // the SUFFIX
                    re += x * std::cos (2.0 * kPi * (double) b * (double) i / (double) n);
                    im -= x * std::sin (2.0 * kPi * (double) b * (double) i / (double) n);
                }
                const double P    = 4.0 * (re * re + im * im) / ((double) n * sw2);   // 4|X|²/(N·Σw²)
                const double want = 10.0 * std::log10 (P);
                worst = std::max (worst, std::fabs (want - (double) p.tierBinDb (k, b)));
            }
            ok (worst < 0.02, "tier " + std::to_string (n) + ": bins match the direct DFT within 0.02 dB (worst " + std::to_string (worst) + ")");
        }
    }

    //==========================================================================================
    group ("oracle: a full-scale sine — lobe inside the band → ≈ 0 dBFS; never below the half-bin scallop");
    {
        p.coverSamples = 0;
        double worstLo = 0.0, worstHi = -100.0;
        for (double f = 30.0; f < 20000.0; f *= 1.037)
        {
            const double v = rawSineDb (p, f);
            worstLo = std::min (worstLo, v); worstHi = std::max (worstHi, v);
            const int k = p.tierAt (f, fs);
            const bool blend = k > 0 && f < p.seamHz (k, fs) * std::exp2 (1.0 / 3.0);
            if (! blend && binsPerBandAt (p, k, f) >= 4.0)
                ok (std::fabs (v) <= 0.3, "band holds the lobe at " + std::to_string ((int) f) + " Hz → |reading| ≤ 0.3 dB (got " + std::to_string (v) + ")");
        }
        ok (worstLo >= -3.5, "no sine reads below −3.5 dB anywhere (worst " + std::to_string (worstLo) + ")");
        ok (worstHi <= 0.2,  "a full-scale sine never reads above +0.2 dB (worst " + std::to_string (worstHi) + ")");
        // bin-limited lows of the longest tier: the peak bin alone / ENBW
        const double binHz = binHzOf (p, 0);
        approx (rawSineDb (p, 20.0 * binHz), -10.0 * std::log10 (p.enbwBins (0)), 0.05, "bin-centred 58.6 Hz reads −1.76 dB (peak bin / ENBW)");
        const double half = rawSineDb (p, 20.5 * binHz);
        ok (half > -3.4 && half < -2.9, "half-bin 60 Hz reads the Hann scallop, ≈ −3.2 dB (got " + std::to_string (half) + ")");
        // the pane is a display: sweeping a sine through a seam and its blend is smooth (no step > 0.5 dB per 1/60 oct)
        for (int k = 1; k < 3; ++k)
        {
            double prev = rawSineDb (p, p.seamHz (k, fs) * std::exp2 (-0.4)); double maxStep = 0.0;
            for (double o = -0.4 + 1.0 / 60.0; o <= 0.75; o += 1.0 / 60.0)
            {
                const double v = rawSineDb (p, p.seamHz (k, fs) * std::exp2 (o));
                maxStep = std::max (maxStep, std::fabs (v - prev)); prev = v;
            }
            ok (maxStep <= 0.5, "sine sweep through seam " + std::to_string (k) + " + blend: max step per 1/60 oct " + std::to_string (maxStep) + " dB ≤ 0.5");
        }
    }

    //==========================================================================================
    group ("oracle: DC is a half cell (one-sided weight ½)");
    {
        p.reset(); p.smoothCoeff = 1.0f;
        fillConst (p.frameInput(), N, 1.0f);
        p.ingest (14);
        approx ((double) p.tierBinDb (0, 0), 10.0 * std::log10 (4.0 / p.enbwBins (0)), 0.02, "raw DC bin of a unit constant = 4/ENBW (+4.26 dB)");
        approx (p.readDb (1.0, fs), 10.0 * std::log10 (0.5 * 4.0 / p.enbwBins (0)), 0.02, "a bin-limited read at DC = half the cell (+1.25 dB)");
    }

    //==========================================================================================
    group ("Welch cover: sub-windows over the hop, and a click in the hop gap is seen");
    {
        p.coverSamples = 1600;
        ok (p.subWindows (2, 14) == 4 && p.subWindows (1, 14) == 1 && p.subWindows (0, 14) == 1, "hop 1600: the 1024 tier runs 4 half-overlapped windows, longer tiers one");
        p.coverSamples = 0;   ok (p.subWindows (2, 14) == 1, "cover 0 → one suffix window");
        p.coverSamples = 1 << 20; ok (p.subWindows (2, 14) == 1 + (N - 1024) / 512, "a huge cover is clamped to what the frame holds");
        p.coverSamples = 1600; ok (p.subWindows (2, 12) == 4, "a 4096 frame at hop 1600 still runs the 4 the hop needs");
        p.coverSamples = 1 << 20; ok (p.subWindows (2, 12) == 1 + (4096 - 1024) / 512, "a shorter frame holds fewer when the hop wants more");
        ok (p.subWindows (0, 12) == 0, "a tier longer than the frame runs none");

        auto click = [&] (int cover)
        {
            p.reset(); p.coverSamples = cover; p.smoothCoeff = 1.0f;
            float* in = p.frameInput(); fillConst (in, N, 0.0f);
            for (int i = N - 1500; i < N - 1200; ++i) in[i] = (float) std::sin (2.0 * kPi * 8000.0 * (double) i / fs);   // inside the hop, outside the last 1024
            p.ingest (14);
            return p.tierBandDb (2, 8000.0, fs);
        };
        const double blind = click (0), seen = click (1600);
        ok (blind <= -100.0, "cover 0: the 1024 tier is blind to a click in the hop gap (" + std::to_string (blind) + " dB)");
        ok (seen >= -20.0,   "cover 1600: the click registers in the 1024 tier (" + std::to_string (seen) + " dB)");
        p.coverSamples = 1600;
    }

    //==========================================================================================
    group ("time resolution: the tiers read different windows (Hann-weighted, end-aligned)");
    {
        const double f = 171.0 * fs / 1024.0;          // bin-centred in every tier (171·16 in 16384)
        p.reset(); p.coverSamples = 0; p.smoothCoeff = 1.0f;
        float* in = p.frameInput(); fillConst (in, N, 0.0f);
        fillSine (in + (N - 512), 512, f, 1.0);         // burst in the last 512 samples
        p.ingest (14);
        const double s2 = p.tierBandDb (2, f, fs), s0 = p.tierBandDb (0, f, fs);
        ok (s2 > -4.5 && s2 < -2.0, "1024 tier: the trailing half of its Hann → ≈ −3 dB (got " + std::to_string (s2) + ")");
        ok (s0 < -40.0,             "16384 tier: the Hann tail → below −40 dB (got " + std::to_string (s0) + ")");
        p.reset();
        in = p.frameInput(); fillConst (in, N, 0.0f);
        fillSine (in + (N - 768), 512, f, 1.0);         // burst centred in the last 1024
        p.ingest (14);
        const double c2 = p.tierBandDb (2, f, fs), c0 = p.tierBandDb (0, f, fs);
        ok (c2 > -1.5 && c2 < 0.3, "1024 tier: a burst under its Hann centre reads ≈ 0 dB (got " + std::to_string (c2) + ")");
        ok (c0 < -25.0,            "16384 tier: the same burst is 25+ dB down (got " + std::to_string (c0) + ")");
        p.coverSamples = 1600;
    }

    //==========================================================================================
    group ("noise ensemble: seam invariance, +3 dB/oct, analytic calibration, hop independence, pink flat");
    {
        auto seamResidue = [&] (int k, bool peak)      // mean over the blend zone of (longer tier − shorter tier)
        {
            double acc = 0.0; int n = 0;
            const double f0 = p.seamHz (k, fs);
            for (double o = 0.0; o <= 1.0 / 3.0; o += 1.0 / 120.0)
            {
                const double f = f0 * std::exp2 (o);
                acc += peak ? (p.tierBandPeakDb (k - 1, f, fs) - p.tierBandPeakDb (k, f, fs))
                            : (p.tierBandDb     (k - 1, f, fs) - p.tierBandDb     (k, f, fs));
                ++n;
            }
            return acc / (double) n;
        };
        auto runWhite = [&] (int cover, int ticks)
        {
            p.reset(); p.coverSamples = cover; p.smoothCoeff = 0.05f;
            Rng r; streamFrames (p, [&] { return r.uni(); }, ticks);
        };

        runWhite (1600, 800);
        for (int k = 1; k < 3; ++k)
        {
            const double d = seamResidue (k, false);
            ok (std::fabs (d) <= 0.6, "white noise, seam " + std::to_string (k) + ": longer vs shorter tier agree within 0.6 dB (Δ " + std::to_string (d) + ")");
            const double dp = seamResidue (k, true);
            ok (std::fabs (dp) <= 1.5, "  … and the peak trace within 1.5 dB (Δ " + std::to_string (dp) + ")");
        }
        std::vector<std::pair<double, double>> pts;
        for (double f = 300.0; f <= 9600.0; f *= 2.0) pts.push_back ({ f, p.readDb (f, fs) });
        const double sl = slopeDbPerOct (pts);
        ok (sl > 2.5 && sl < 3.5, "white noise reads +3 dB/oct across every tier (" + std::to_string (sl) + " dB/oct)");
        {
            const double f = 1200.0, B = f * kBandFactor, sigma2 = 1.0 / 3.0;   // uniform [−1,1): σ² = 1/3
            const double want = 10.0 * std::log10 (4.0 * sigma2 * B / fs);
            approx (p.readDb (f, fs), want, 0.5, "calibration: white noise reads 10·log10(4σ²B/fs) at 1.2 kHz (" + std::to_string (want) + " dB)");
        }
        for (int cover : { 800, 3200 })
        {
            runWhite (cover, 600);
            const double d = seamResidue (2, false);
            ok (std::fabs (d) <= 0.6, "hop " + std::to_string (cover) + ": the 1024 seam residue stays within 0.6 dB (Δ " + std::to_string (d) + ")");
        }
        {
            p.reset(); p.coverSamples = 1600; p.smoothCoeff = 0.05f;
            Rng r; Pink pk;
            for (int i = 0; i < 4096; ++i) (void) pk (r.uni());   // settle the filter
            streamFrames (p, [&] { return pk (r.uni()); }, 800);
            pts.clear();
            for (double f = 300.0; f <= 9600.0; f *= 1.4142) pts.push_back ({ f, p.readDb (f, fs) });
            const double ps = slopeDbPerOct (pts);
            ok (std::fabs (ps) <= 0.7, "pink noise reads flat across every seam (" + std::to_string (ps) + " dB/oct)");
        }
    }

    //==========================================================================================
    group ("Parseval on a seeded frame: contiguous bands tile the axis without double counting");
    {
        p.reset(); p.coverSamples = 0; p.smoothCoeff = 1.0f;
        Rng r; streamFrames (p, [&] { return r.uni(); }, 1);
        const int    k = 0; const int bins = p.tierBins (k); const double binHz = binHzOf (p, k);
        const double h = std::exp2 (1.0 / 48.0);
        const double fa = 400.0, fb = 20000.0;                                  // bands ≥ 4 bins wide here
        double bandSum = 0.0; double e = fa; double lastEdge = fa;
        while (e * h * h <= fb)
        {
            bandSum += std::pow (10.0, p.tierBandDb (k, e * h, fs) / 10.0);   // band centred at e·h spans [e, e·h²]
            e *= h * h; lastEdge = e;
        }
        double binSum = 0.0;                                                     // the same interval, integrated over bins by hand
        const double uA = fa / binHz + 0.5, uB = lastEdge / binHz + 0.5;
        for (int i = 0; i < bins; ++i)
        {
            const double lo = std::max (uA, (double) i), hi = std::min (uB, (double) i + 1.0);
            if (hi > lo) binSum += (hi - lo) * std::pow (10.0, (double) p.tierBinDb (k, i) / 10.0);
        }
        approx (10.0 * std::log10 (bandSum), 10.0 * std::log10 (binSum), 0.02, "Σ bands over [400 Hz, ~20 kHz] == Σ bin power over the same interval");
    }

    //==========================================================================================
    group ("fractional integration: random bands vs an O(N) reference (incl. the bin-limited rule)");
    {
        p.reset(); p.coverSamples = 0; p.smoothCoeff = 1.0f;
        Rng r; streamFrames (p, [&] { return r.uni(); }, 1);
        double worst = 0.0;
        for (int trial = 0; trial < 60; ++trial)
        {
            const int    k     = trial % 3;
            const double octs  = 0.005 + 0.5 * (r.uni() + 1.0) * 0.25;           // 1/200 .. 1/4 oct
            const double f     = 20.0 * std::pow (600.0, 0.5 * (r.uni() + 1.0));  // 20 Hz .. 12 kHz
            p.bandOctaves = octs;
            const double got   = p.tierBandDb (k, f, fs);
            const double bh    = binHzOf (p, k); const int bins = p.tierBins (k);
            const double hh    = std::exp2 (0.5 * octs);
            const double uLo   = std::clamp ((f / hh) / bh + 0.5, 0.5, bins - 0.5), uHi = std::clamp ((f * hh) / bh + 0.5, 0.5, bins - 0.5);
            double acc = 0.0;
            for (int i = 0; i < bins; ++i)
            {
                const double lo = std::max (uLo, (double) i), hi = std::min (uHi, (double) i + 1.0);
                if (hi > lo) acc += (hi - lo) * std::pow (10.0, (double) p.tierBinDb (k, i) / 10.0);
            }
            const double span = uHi - uLo;
            if (span < 1.0) { const int c = (int) (f / bh + 0.5); acc = acc / span * ((c <= 0 || c >= bins - 1) ? 0.5 : 1.0); }
            const double want = 10.0 * std::log10 (acc);
            worst = std::max (worst, std::fabs (want - got));
        }
        p.bandOctaves = 1.0 / 24.0;
        ok (worst < 0.01, "60 random bands across 3 tiers match the reference within 0.01 dB (worst " + std::to_string (worst) + ")");
    }

    //==========================================================================================
    group ("smoothing + peak-hold law (one tier, bin-centred sine)");
    {
        auto q = std::make_unique<Pane>();
        const int one[] = { 11 }; q->setTiers (one, 1); q->coverSamples = 0; q->smoothCoeff = 0.25f; q->peakFallDb = 0.8f;
        const double f = 100.0 * fs / 2048.0;
        auto feed = [&] (double amp) { fillSine (q->frameInput(), N, f, amp); q->ingest (14); };
        feed (1.0);
        const double p0 = std::pow (10.0, (double) q->tierBinDb (0, 100) / 10.0);
        approx (10.0 * std::log10 (p0), -10.0 * std::log10 (q->enbwBins (0)), 0.02, "seed: the first frame lands directly (peak bin / ENBW)");
        feed (0.1);                                                              // −20 dB target → power 0.01·p0
        const double p1 = p0 + 0.25 * (0.01 * p0 - p0);
        approx ((double) q->tierBinDb (0, 100), 10.0 * std::log10 (p1), 0.03, "tick 1: power moved 25 % of the way (one-pole on POWER)");
        approx ((double) q->tierBinPeak (0, 100), std::max (10.0 * std::log10 (p0) - 0.8, 10.0 * std::log10 (p1)), 0.03, "tick 1: peak fell 0.8 dB");
        feed (0.1);
        const double p2 = p1 + 0.25 * (0.01 * p0 - p1);
        approx ((double) q->tierBinDb (0, 100), 10.0 * std::log10 (p2), 0.03, "tick 2: one-pole again");
        approx ((double) q->tierBinPeak (0, 100), 10.0 * std::log10 (p0) - 1.6, 0.03, "tick 2: peak −1.6 dB");
        const float held = q->tierBinDb (0, 100);
        for (int t = 0; t < 15; ++t) q->starve();
        ok (q->tierBinDb (0, 100) == held, "starve ×15: the fill holds");
        approx ((double) q->tierBinPeak (0, 100), 10.0 * std::log10 (p0) - 1.6 - 15 * 0.8, 0.03, "starve ×15: the peak keeps falling");
        q->starve();
        approx ((double) q->tierBinDb (0, 100), (double) held + 10.0 * std::log10 (0.9), 0.03, "starve #16: the fill fades (−0.46 dB/tick)");
        q->reset();
        ok (q->readDb (f, fs) == (double) Pane::kFloorDb && q->readPeakDb (f, fs) == (double) Pane::kFloorDb, "reset → silence");
    }

    //==========================================================================================
    group ("silence reads the floor whatever the band width; a lone tone fabricates nothing far away");
    {
        p.reset(); p.coverSamples = 1600; p.smoothCoeff = 1.0f;
        bool floorAll = true;
        for (double f = 20.0; f < 24000.0; f *= 1.1) if (p.readDb (f, fs) != (double) Pane::kFloorDb || p.readPeakDb (f, fs) != (double) Pane::kFloorDb) floorAll = false;
        ok (floorAll, "after reset every read is −120 dB");
        fillConst (p.frameInput(), N, 0.0f); p.ingest (14);
        floorAll = true;
        for (double f = 20.0; f < 24000.0; f *= 1.1) if (p.readDb (f, fs) != (double) Pane::kFloorDb) floorAll = false;
        ok (floorAll, "a frame of zeros reads −120 dB in every band, wide or narrow");
        fillSine (p.frameInput(), N, 1000.0, 1.0e-5); p.ingest (14);          // −100 dB tone
        ok (p.readDb (10000.0, fs) == (double) Pane::kFloorDb, "a −100 dB tone at 1 kHz leaves 10 kHz at the floor");
        ok (p.readDb (1000.0, fs) > -103.0 && p.readDb (1000.0, fs) < -99.0, "…and reads ≈ −100 dB at 1 kHz");
    }

    //==========================================================================================
    group ("adversarial inputs");
    {
        const double nan = std::numeric_limits<double>::quiet_NaN(), inf = std::numeric_limits<double>::infinity();
        p.reset(); p.smoothCoeff = 0.25f;
        float* in = p.frameInput(); fillSine (in, N, 1000.0, 0.5);
        in[100] = (float) nan; in[N - 7] = (float) inf; in[N - 3000] = (float) -inf;
        p.ingest (14);
        bool finite = true;
        for (double f = 20.0; f < 24000.0; f *= 1.05) if (! std::isfinite (p.readDb (f, fs)) || ! std::isfinite (p.readPeakDb (f, fs))) finite = false;
        ok (finite, "NaN/Inf samples: every read stays finite");
        p.reset(); fillSine (p.frameInput(), N, 1000.0, 0.5); p.ingest (14);
        ok (p.readDb (0.0, fs) == (double) Pane::kFloorDb && p.readDb (-5.0, fs) == (double) Pane::kFloorDb, "f ≤ 0 → floor");
        ok (p.readDb (1000.0, 0.0) == (double) Pane::kFloorDb && p.readDb (1000.0, nan) == (double) Pane::kFloorDb, "fs ≤ 0 / NaN → floor");
        ok (p.readDb (nan, fs) == (double) Pane::kFloorDb && p.readDb (inf, fs) == (double) Pane::kFloorDb, "f NaN/Inf → floor");
        ok (std::fabs (p.readDb (30000.0, fs) - p.readDb (0.5 * fs / std::exp2 (1.0 / 48.0), fs)) < 1e-9, "above Nyquist repeats the last band that fits");
        // a frame shorter than the longest tier: the long tier holds, the shorter ones update
        p.reset();
        fillSine (p.frameInput(), N, 5000.0, 1.0); p.ingest (12);            // 4096-sample frame
        ok (p.tierBandDb (0, 5000.0, fs) == (double) Pane::kFloorDb, "order-12 frame: the 16384 tier holds (still silent)");
        ok (p.tierBandDb (2, 5000.0, fs) > -1.0, "…while the 1024 tier read it");
        fillSine (p.frameInput(), N, 5000.0, 1.0); p.ingest (99);            // clamps to 14
        ok (p.tierBandDb (0, 5000.0, fs) > -1.0, "order 99 clamps to the frame size and serves every tier");
        p.blendOctaves = 0.0;
        ok (std::isfinite (p.readDb (p.seamHz (1, fs), fs)), "blend width 0 = a hard seam, still finite");
        p.blendOctaves = 1.0 / 3.0;
        p.bandOctaves = -1.0;
        ok (std::isfinite (p.readDb (1000.0, fs)), "a non-positive band width is clamped, not a division by zero");
        p.bandOctaves = 1.0 / 24.0;
    }

    //==========================================================================================
    group ("column geometry + tilt, and no allocation on the hot path");
    {
        p.reset(); p.coverSamples = 1600; p.smoothCoeff = 0.25f;
        Rng r; streamFrames (p, [&] { return r.uni(); }, 3);
        analysis::PlotMap pm; pm.width = 640.0f; pm.height = 96.0f; pm.plotBottom = 96.0f; pm.specTop = 6.0; pm.specBottom = -90.0;
        int count = 0; float lastX = -1.0f; bool ascending = true, firstAtZero = false, consistent = true;
        p.buildColumns (pm, fs, 4.5, 1000.0, [&] (int i, float x, float yF, float yP)
        {
            if (i == 0) firstAtZero = (x == 0.0f);
            if (x < lastX) ascending = false; lastX = x;
            const double f = pm.xToFreq (x), tilt = 4.5 * std::log2 (f / 1000.0);
            if (std::fabs (yF - pm.specDbToY (p.readDb (f, fs) + tilt)) > 1e-4f || std::fabs (yP - pm.specDbToY (p.readPeakDb (f, fs) + tilt)) > 1e-4f) consistent = false;
            ++count;
        });
        ok (count == 641 && ascending && firstAtZero, "N+1 = 641 ascending points, i = 0 at x = 0");
        ok (consistent, "every column is specDbToY (readDb + tilt) — one geometry source");
        analysis::PlotMap z; z.width = 0.0f; z.height = 10.0f;
        count = 0; p.buildColumns (z, fs, 0.0, 1000.0, [&] (int, float, float, float) { ++count; });
        ok (count == 257, "a zero-width map still emits the minimum 257 columns");

        const long before = g_allocs.load();
        streamFrames (p, [&] { return r.uni(); }, 3);
        p.starve();
        p.buildColumns (pm, fs, 4.5, 1000.0, [] (int, float, float, float) {});
        (void) p.readDb (1000.0, fs); (void) p.tierAt (1000.0, fs);
        test::okNoAlloc (g_allocs.load() == before, "ingest / starve / buildColumns / reads allocate nothing");
    }

    return test::report();
}

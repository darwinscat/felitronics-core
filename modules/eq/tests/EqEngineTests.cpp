// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// Verifies the runtime engine by actually pushing audio through it: measured steady-state gain
// of a sine must agree with the analytic magnitudeDb() readout, smoothing converges, the swept
// SVF band matches the matched band, a 24 dB/oct HP rolls off ~24 dB/oct, and process() is clean
// (silence stays silent, no NaN). JUCE-free.

#include "TestUtil.h"

#include <teq/EqBand.h>
#include <teq/EqEngine.h>
#include <teq/Smoother.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <string>
#include <vector>

// RT-safety: count every heap allocation so the audio path can be asserted alloc-free (as in the
// other module suites). Only ONE TU in this executable may replace global operator new.
static std::atomic<long> g_allocs { 0 };
void* operator new      (std::size_t s) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s ? s : 1); }
void* operator new[]    (std::size_t s) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s ? s : 1); }
void  operator delete   (void* p) noexcept { std::free (p); }
void  operator delete[] (void* p) noexcept { std::free (p); }
void  operator delete   (void* p, std::size_t) noexcept { std::free (p); }
void  operator delete[] (void* p, std::size_t) noexcept { std::free (p); }

using namespace teq;

namespace
{
    // Steady-state gain (dB) of a mono in-place processing callback at frequency f.
    template <class Apply>
    double sineGainDb (Apply&& apply, double f, double fs, int warm = 600)
    {
        const int block = 256;                 // warm=600 -> ~3.2 s >> 30 ms smoothing + filter ring
        const double amp = 0.25, dp = 2.0 * kPi * f / fs;
        std::vector<float> buf ((size_t) block);
        double phase = 0.0;
        auto fill = [&] { for (int n = 0; n < block; ++n) { buf[(size_t) n] = (float) (amp * std::sin (phase));
                                                            phase += dp; if (phase > 2.0 * kPi) phase -= 2.0 * kPi; } };

        for (int b = 0; b < warm; ++b) { fill(); float* p = buf.data(); apply (&p, 1, block); }

        fill();
        std::vector<float> in (buf.begin(), buf.end());
        { float* p = buf.data(); apply (&p, 1, block); }

        double ir = 0.0, orr = 0.0;
        for (int n = 0; n < block; ++n)
        {
            ir  += (double) in[(size_t) n]  * in[(size_t) n];
            orr += (double) buf[(size_t) n] * buf[(size_t) n];
        }
        return 20.0 * std::log10 (std::sqrt (orr / std::max (1e-30, ir)));
    }

    // Per-channel steady-state gain (dB) for an N-channel in-place callback. Every channel gets the
    // SAME sine but a DISTINCT amplitude, so any cross-channel state bleed shows up as a wrong gain.
    template <class Apply>
    std::vector<double> multiSineGainDb (Apply&& apply, int C, double f, double fs, int warm = 600)
    {
        const int block = 256;
        const double dp = 2.0 * kPi * f / fs;
        std::vector<std::vector<float>> bufs ((size_t) C, std::vector<float> ((size_t) block));
        std::vector<float*> ptr ((size_t) C);
        for (int c = 0; c < C; ++c) ptr[(size_t) c] = bufs[(size_t) c].data();

        double phase = 0.0;
        auto fill = [&]
        {
            for (int n = 0; n < block; ++n)
            {
                const double s = std::sin (phase);
                phase += dp; if (phase > 2.0 * kPi) phase -= 2.0 * kPi;
                for (int c = 0; c < C; ++c) bufs[(size_t) c][(size_t) n] = (float) ((0.1 + 0.04 * c) * s);
            }
        };

        for (int b = 0; b < warm; ++b) { fill(); apply (ptr.data(), C, block); }

        fill();
        const std::vector<std::vector<float>> in = bufs;     // snapshot the inputs before in-place process
        apply (ptr.data(), C, block);

        std::vector<double> g ((size_t) C);
        for (int c = 0; c < C; ++c)
        {
            double ir = 0.0, orr = 0.0;
            for (int n = 0; n < block; ++n)
            {
                ir  += (double) in[(size_t) c][(size_t) n]   * in[(size_t) c][(size_t) n];
                orr += (double) bufs[(size_t) c][(size_t) n] * bufs[(size_t) c][(size_t) n];
            }
            g[(size_t) c] = 20.0 * std::log10 (std::sqrt (orr / std::max (1e-30, ir)));
        }
        return g;
    }

    bool anyNaN (const float* d, int n) { for (int i = 0; i < n; ++i) if (! std::isfinite (d[i])) return true; return false; }
}

// Frozen v0.1.x process semantics for the T1 reference-NULL: the Stereo per-channel loop and the M/S
// delta-fold, verbatim from the pre-lanes engine, but driven by the SAME designBand() coefficients the
// v2 engine uses (fed through laneView, so the coefficients are shared and only the ROUTING is being
// nulled). Snap-on-set (no ramp), so a static band is bit-exact against EqBand.
namespace refv01
{
    struct Band
    {
        double fs = 48000.0;
        int    ch = 2;
        int    nST = 0, nM = 0, nS = 0;
        bool   stRun = false, midRun = false, sideRun = false;
        BiquadCoeffs cST[EqBand::kMaxSections], cM[EqBand::kMaxSections], cS[EqBand::kMaxSections];
        Biquad bST[EqBand::kMaxSections][EqBand::kMaxChannels];
        Biquad bM[EqBand::kMaxSections], bS[EqBand::kMaxSections];

        void prepare (double sr, int nc) noexcept
        {
            fs = sr; ch = nc;
            for (int s = 0; s < EqBand::kMaxSections; ++s)
            {
                for (int c = 0; c < EqBand::kMaxChannels; ++c) bST[s][c].reset();
                bM[s].reset(); bS[s].reset();
            }
        }

        void set (const BandParams& b) noexcept
        {
            const BandParams vST = laneView (b, Lane::Stereo);
            const BandParams vM  = laneView (b, Lane::Mid);
            const BandParams vS  = laneView (b, Lane::Side);
            stRun   = vST.on && ! vST.bypass;
            midRun  = vM.on  && ! vM.bypass;
            sideRun = vS.on  && ! vS.bypass;
            nST = nM = nS = 0;
            if (stRun)   { const auto d = designBand (vST, fs); nST = d.n; for (int s = 0; s < d.n; ++s) { cST[s] = d.sec[s]; for (int c = 0; c < ch; ++c) bST[s][c].setCoeffs (cST[s]); } }
            if (midRun)  { const auto d = designBand (vM,  fs); nM  = d.n; for (int s = 0; s < d.n; ++s) { cM[s]  = d.sec[s]; bM[s].setCoeffs (cM[s]); } }
            if (sideRun) { const auto d = designBand (vS,  fs); nS  = d.n; for (int s = 0; s < d.n; ++s) { cS[s]  = d.sec[s]; bS[s].setCoeffs (cS[s]); } }
        }

        // {st}-only: the frozen v0.1.x per-channel loop.
        void processStereo (float* const* c, int nc, int n) noexcept
        {
            const int C = nc < ch ? nc : ch;
            for (int ci = 0; ci < C; ++ci)
                for (int i = 0; i < n; ++i)
                { float x = c[ci][i]; for (int s = 0; s < nST; ++s) x = bST[s][ci].processSample (x); c[ci][i] = x; }
        }

        // {m,s}: the frozen v0.1.x delta-fold (2-channel). midRun / sideRun gate each lane.
        void processMS (float* const* c, int n) noexcept
        {
            float* L = c[0]; float* R = c[1];
            for (int i = 0; i < n; ++i)
            {
                const float m = 0.5f * (L[i] + R[i]);
                const float s = 0.5f * (L[i] - R[i]);
                float dM = 0.0f, dS = 0.0f;
                if (midRun)  { float x = m; for (int k = 0; k < nM; ++k) x = bM[k].processSample (x); dM = x - m; }
                if (sideRun) { float y = s; for (int k = 0; k < nS; ++k) y = bS[k].processSample (y); dS = y - s; }
                L[i] += dM + dS;
                R[i] += dM - dS;
            }
        }
    };
}

void runEqEngineTests()
{
    using namespace teqtest;
    const double fs = 48000.0;

    group ("Smoother: converges to target, snap settles");
    {
        Smoother s; s.prepare (fs, 10.0); s.snap (0.0); s.setTarget (1.0);
        s.advance (48000);
        expectNear (s.value(), 1.0, 1e-3, "advance(1s) converges");
        s.snap (5.0);
        expectTrue (s.settled(), "snap settles");
    }

    group ("EqBand bell: measured gain == response() at centre, ~unity far off");
    {
        EqBand band; band.prepare (fs, 1);
        BandParams p; p.on = true; p.type = FilterType::Bell; p.lane (Lane::Stereo).freq = 1000.0; p.lane (Lane::Stereo).Q = 2.0; p.lane (Lane::Stereo).gainDb = 6.0;
        band.setParams (p);
        const double meas = sineGainDb ([&] (float* const* ch, int nc, int n) { band.processBlock (ch, nc, n); }, 1000.0, fs);
        const double resp = 20.0 * std::log10 (std::abs (band.response (2.0 * kPi * 1000.0 / fs)));
        std::printf ("      bell @1k: measured=%.3f dB  response=%.3f dB\n", meas, resp);
        expectNear (meas, 6.0, 0.4,  "measured +6 dB");
        expectNear (resp, 6.0, 0.05, "response +6 dB");
        expectNear (meas, resp, 0.4, "measured ~ response");

        EqBand far; far.prepare (fs, 1);
        BandParams q; q.on = true; q.type = FilterType::Bell; q.lane (Lane::Stereo).freq = 1000.0; q.lane (Lane::Stereo).Q = 4.0; q.lane (Lane::Stereo).gainDb = 9.0;
        far.setParams (q);
        const double low = sineGainDb ([&] (float* const* ch, int nc, int n) { far.processBlock (ch, nc, n); }, 60.0, fs);
        expectNear (low, 0.0, 0.3, "unity an octave+ away");
    }

    group ("EqBand swept (SVF) bell tracks the same +6 dB");
    {
        EqBand band; band.prepare (fs, 1);
        BandParams p; p.on = true; p.type = FilterType::Bell; p.lane (Lane::Stereo).freq = 2000.0; p.lane (Lane::Stereo).Q = 3.0; p.lane (Lane::Stereo).gainDb = 6.0; p.swept = true;
        band.setParams (p);
        const double meas = sineGainDb ([&] (float* const* ch, int nc, int n) { band.processBlock (ch, nc, n); }, 2000.0, fs);
        std::printf ("      swept bell @2k: measured=%.3f dB\n", meas);
        expectNear (meas, 6.0, 0.5, "SVF bell +6 dB");
    }

    group ("EqBand HP 24 dB/oct rolls off ~24 dB/octave in the stopband");
    {
        EqBand band; band.prepare (fs, 1);
        BandParams p; p.on = true; p.type = FilterType::HighPass; p.lane (Lane::Stereo).freq = 1000.0; p.lane (Lane::Stereo).Q = 0.707; p.lane (Lane::Stereo).slope = 24;
        band.setParams (p);
        auto ap = [&] (float* const* ch, int nc, int n) { band.processBlock (ch, nc, n); };
        const double g125 = sineGainDb (ap, 125.0, fs);
        const double g250 = sineGainDb (ap, 250.0, fs);
        std::printf ("      HP24: g(125)=%.2f  g(250)=%.2f  slope=%.1f dB/oct\n", g125, g250, g250 - g125);
        expectTrue ((g250 - g125) > 20.0 && (g250 - g125) < 28.0, "~24 dB/oct");
    }

    group ("EqEngine: two bells cascade, measured == magnitudeDb");
    {
        EqEngine eng; eng.prepare (fs, 512, 1);
        BandParams a; a.on = true; a.type = FilterType::Bell; a.lane (Lane::Stereo).freq = 300.0;  a.lane (Lane::Stereo).Q = 1.0; a.lane (Lane::Stereo).gainDb =  4.0;
        BandParams b; b.on = true; b.type = FilterType::Bell; b.lane (Lane::Stereo).freq = 3000.0; b.lane (Lane::Stereo).Q = 1.0; b.lane (Lane::Stereo).gainDb = -3.0;
        eng.setBand (0, a); eng.setBand (1, b);
        const double m = sineGainDb ([&] (float* const* ch, int nc, int n) { eng.process (ch, nc, n); }, 300.0, fs);
        const double r = eng.magnitudeDb (300.0);
        std::printf ("      @300: measured=%.3f dB  magnitudeDb=%.3f dB\n", m, r);
        expectNear (m, r, 0.4, "measured ~ magnitudeDb");
        expectNear (r, 4.0, 0.6, "~+4 dB at band-0 centre");
    }

    group ("swept SVF BandPass has unity gain at centre");
    {
        EqBand band; band.prepare (fs, 1);
        BandParams p; p.on = true; p.type = FilterType::BandPass; p.lane (Lane::Stereo).freq = 1000.0; p.lane (Lane::Stereo).Q = 4.0; p.swept = true;
        band.setParams (p);
        const double meas = sineGainDb ([&] (float* const* ch, int nc, int n) { band.processBlock (ch, nc, n); }, 1000.0, fs);
        std::printf ("      swept BP @1k Q=4: measured=%.3f dB\n", meas);
        expectNear (meas, 0.0, 0.5, "unity (0 dB) at centre");
    }

    group ("swept Notch is a REAL notch (SVF low+high): deep cut at fc, unity off-band — not a pass-through");
    {
        // Regression: FilterType::Notch with swept=true used to be a bit-exact audio pass-through
        // (the SVF mapped it to m0=1,m1=0,m2=0) while the display showed a full notch curve. The SVF
        // now realises notch = low + high = v0 − k·v1, with an exact null at the prewarped fc.
        EqBand band; band.prepare (fs, 1);
        BandParams p; p.on = true; p.type = FilterType::Notch; p.lane (Lane::Stereo).freq = 1000.0; p.lane (Lane::Stereo).Q = 2.0; p.swept = true;
        band.setParams (p);
        auto ap = [&] (float* const* ch, int nc, int n) { band.processBlock (ch, nc, n); };
        const double atFc  = sineGainDb (ap, 1000.0, fs);
        const double below = sineGainDb (ap, 125.0,  fs);
        const double above = sineGainDb (ap, 8000.0, fs);
        std::printf ("      swept notch: fc=%.1f  125Hz=%.2f  8k=%.2f dB\n", atFc, below, above);
        expectTrue (atFc < -30.0,      "deep cut at fc (audio was a bit-exact pass-through)");
        expectNear (below, 0.0, 0.3,   "unity well below fc");
        expectNear (above, 0.0, 0.3,   "unity well above fc");
    }

    group ("swept Tilt runs the MATCHED two-shelf tilt (non-sweepable): audio == displayed curve");
    {
        // Regression: FilterType::Tilt with swept=true used to be a bit-exact audio pass-through while
        // response()/bandResponse displayed the full ±gain tilt (up to ±30 dB of GUI lying). A one-SVF
        // tilt cannot hold a unity pivot beyond ~6 dB, so Tilt now ignores the swept flag entirely.
        EqBand band; band.prepare (fs, 1);
        BandParams p; p.on = true; p.type = FilterType::Tilt; p.lane (Lane::Stereo).freq = 1000.0; p.lane (Lane::Stereo).Q = 1.0; p.lane (Lane::Stereo).gainDb = 12.0; p.swept = true;
        band.setParams (p);
        auto ap = [&] (float* const* ch, int nc, int n) { band.processBlock (ch, nc, n); };
        const double lows  = sineGainDb (ap, 60.0,    fs);
        const double highs = sineGainDb (ap, 12000.0, fs);
        const double rLow  = 20.0 * std::log10 (std::abs (band.response (2.0 * kPi * 60.0    / fs)));
        const double rHigh = 20.0 * std::log10 (std::abs (band.response (2.0 * kPi * 12000.0 / fs)));
        std::printf ("      swept tilt: lows=%.2f (disp %.2f)  highs=%.2f (disp %.2f) dB\n", lows, rLow, highs, rHigh);
        expectNear (lows,  -12.0, 1.0, "lows ~ -12 dB (tilt is audible, not silent)");
        expectNear (highs,  12.0, 1.0, "highs ~ +12 dB");
        expectNear (lows,  rLow,  0.4, "audio == displayed curve (lows)");
        expectNear (highs, rHigh, 0.4, "audio == displayed curve (highs)");
    }

    group ("swept SVF high shelf is Butterworth (plateau == gain, Q-independent)");
    {
        auto plateau = [&] (double Q)
        {
            EqBand band; band.prepare (fs, 1);
            BandParams p; p.on = true; p.type = FilterType::HighShelf; p.lane (Lane::Stereo).freq = 3000.0; p.lane (Lane::Stereo).Q = Q; p.lane (Lane::Stereo).gainDb = 12.0; p.swept = true;
            band.setParams (p);
            return sineGainDb ([&] (float* const* ch, int nc, int n) { band.processBlock (ch, nc, n); }, 14000.0, fs);
        };
        const double g1 = plateau (1.0), g4 = plateau (4.0);
        std::printf ("      swept HS plateau @14k: Q=1 -> %.3f dB, Q=4 -> %.3f dB\n", g1, g4);
        expectNear (g1, 12.0, 0.6, "plateau ~ +12 dB");
        expectNear (g1, g4,   0.2, "Q ignored for shelves");
    }

    group ("parameter safety: bad Q / huge gain / NaN / Inf do not produce NaN");
    {
        const double inf = std::numeric_limits<double>::infinity();
        EqEngine eng; eng.prepare (fs, 256, 1);
        BandParams p; p.on = true; p.type = FilterType::Bell;
        p.lane (Lane::Stereo).freq = std::nan (""); p.lane (Lane::Stereo).Q = 0.0; p.lane (Lane::Stereo).gainDb = inf;       // NaN freq, Q=0, +Inf gain
        eng.setBand (0, p);
        const int n = 256;
        std::vector<float> buf ((size_t) n);
        for (int i = 0; i < n; ++i) buf[(size_t) i] = 0.2f * (float) std::sin (2.0 * kPi * 1000.0 * i / fs);
        float* ch[1] = { buf.data() };
        for (int b = 0; b < 20; ++b) eng.process (ch, 1, n);
        expectTrue (! anyNaN (buf.data(), n),                                       "sanitised params -> finite output");
        expectTrue (std::isfinite (eng.magnitudeDb (1000.0)),                       "magnitudeDb finite");
        expectTrue (std::isfinite (EqEngine::magnitudeDbFor (&p, 1, 1000.0, fs)),   "magnitudeDbFor finite on bad params");
    }

    group ("denormal flush: tail decays to exact zero after silence");
    {
        EqBand band; band.prepare (fs, 1);
        BandParams p; p.on = true; p.type = FilterType::Bell; p.lane (Lane::Stereo).freq = 1000.0; p.lane (Lane::Stereo).Q = 2.0; p.lane (Lane::Stereo).gainDb = 6.0;
        band.setParams (p);
        const int n = 64;
        std::vector<float> buf ((size_t) n, 0.0f);
        buf[0] = 1.0f;                                           // impulse, then silence
        float* ch[1] = { buf.data() };
        band.processBlock (ch, 1, n);
        for (int b = 0; b < 200; ++b) { for (int i = 0; i < n; ++i) buf[(size_t) i] = 0.0f; band.processBlock (ch, 1, n); }
        double tail = 0.0; for (int i = 0; i < n; ++i) tail += std::fabs (buf[(size_t) i]);
        expectTrue (tail == 0.0, "tail flushed to exact zero (no denormal residue)");
    }

    group ("first setParams snaps (no ramp from defaults on load)");
    {
        EqBand band; band.prepare (fs, 1);
        BandParams p; p.on = true; p.type = FilterType::Bell; p.lane (Lane::Stereo).freq = 1000.0; p.lane (Lane::Stereo).Q = 2.0; p.lane (Lane::Stereo).gainDb = 6.0;
        band.setParams (p);
        // measured after only 2 blocks (< 30 ms smoothing): snapped -> already +6 dB, not mid-ramp.
        const double meas = sineGainDb ([&] (float* const* ch, int nc, int n) { band.processBlock (ch, nc, n); }, 1000.0, fs, 2);
        std::printf ("      gain after 2 blocks: %.3f dB\n", meas);
        expectNear (meas, 6.0, 0.5, "target reached immediately (no ramp from defaults)");
    }

    group ("stateless magnitudeDbFor matches the live magnitudeDb");
    {
        EqEngine eng; eng.prepare (fs, 512, 1);
        BandParams arr[EqEngine::kMaxBands] {};
        arr[0].on = true; arr[0].type = FilterType::Bell;      arr[0].lane (Lane::Stereo).freq = 300.0;  arr[0].lane (Lane::Stereo).Q = 1.0; arr[0].lane (Lane::Stereo).gainDb =  4.0;
        arr[1].on = true; arr[1].type = FilterType::HighShelf; arr[1].lane (Lane::Stereo).freq = 8000.0;                 arr[1].lane (Lane::Stereo).gainDb = -5.0;
        eng.setBand (0, arr[0]); eng.setBand (1, arr[1]);
        const int n = 64; std::vector<float> z ((size_t) n, 0.0f); float* ch[1] = { z.data() };
        for (int b = 0; b < 4; ++b) eng.process (ch, 1, n);     // let the live coeffs update
        for (double f : { 100.0, 300.0, 1000.0, 8000.0, 15000.0 })
        {
            const double live      = eng.magnitudeDb (f);
            const double stateless = EqEngine::magnitudeDbFor (arr, EqEngine::kMaxBands, f, fs);
            expectNear (stateless, live, 1e-6, "stateless == live @" + std::to_string ((int) f));
        }
    }

    group ("swept HP slope=24 is single-stage: audio slope == response slope (~12 dB/oct)");
    {
        EqBand band; band.prepare (fs, 1);
        BandParams p; p.on = true; p.type = FilterType::HighPass; p.lane (Lane::Stereo).freq = 1000.0; p.lane (Lane::Stereo).Q = 0.707; p.lane (Lane::Stereo).slope = 24; p.swept = true;
        band.setParams (p);
        auto ap = [&] (float* const* ch, int nc, int n) { band.processBlock (ch, nc, n); };
        const double a125 = sineGainDb (ap, 125.0, fs), a250 = sineGainDb (ap, 250.0, fs);
        const double audioSlope = a250 - a125;
        const double r125 = 20.0 * std::log10 (std::abs (band.response (2.0 * kPi * 125.0 / fs)));
        const double r250 = 20.0 * std::log10 (std::abs (band.response (2.0 * kPi * 250.0 / fs)));
        const double respSlope = r250 - r125;
        std::printf ("      swept HP24: audioSlope=%.1f  responseSlope=%.1f dB/oct\n", audioSlope, respSlope);
        expectTrue (audioSlope > 9.0 && audioSlope < 15.0, "audio ~12 dB/oct (single SVF, not 24)");
        expectNear (respSlope, audioSlope, 2.5,            "response agrees with audio (no 24-vs-12 mismatch)");
    }

    group ("variable HP slopes: analytic rolloff ~ N dB/oct");
    {
        for (int slope : { 6, 12, 24, 36, 48, 72, 96 })
        {
            BandParams p; p.on = true; p.type = FilterType::HighPass; p.lane (Lane::Stereo).freq = 1000.0; p.lane (Lane::Stereo).slope = slope;
            const double r1 = 20.0 * std::log10 (std::abs (bandResponse (p, fs, 2.0 * kPi * 250.0 / fs)));
            const double r2 = 20.0 * std::log10 (std::abs (bandResponse (p, fs, 2.0 * kPi * 500.0 / fs)));
            const double oct = r2 - r1;   // one octave (250->500) in the stopband
            std::printf ("      HP %2d dB/oct: octave rolloff = %.1f dB\n", slope, oct);
            expectTrue (oct > slope * 0.85 && oct < slope * 1.18, "HP slope ~ " + std::to_string (slope) + " dB/oct");
        }
    }

    group ("HP/LP every order incl. ODD: -3.01 dB at fc (Butterworth invariant) + flat passband, no ripple");
    {
        // Regression for the odd-order pole-Q bug: the cos enumeration made a 3rd-order LP read -7.8 dB
        // at fc (corner ~42% low). A Butterworth of ANY order is exactly -3.0103 dB at fc.
        for (int slope : { 6, 12, 18, 24, 30, 36, 42, 48, 54, 66, 78, 90, 96 })
            for (int hp = 0; hp < 2; ++hp)
            {
                const double fc = 1000.0;
                BandParams p; p.on = true; p.type = hp ? FilterType::HighPass : FilterType::LowPass; p.lane (Lane::Stereo).freq = fc; p.lane (Lane::Stereo).slope = slope;
                const std::string at = " (slope " + std::to_string (slope) + (hp ? " HP)" : " LP)");
                auto dbAt = [&] (double f) { return 20.0 * std::log10 (std::abs (bandResponse (p, fs, 2.0 * kPi * f / fs))); };
                expectNear (dbAt (fc), -3.0103, 0.15, "|H(fc)| == -3.01 dB" + at);                 // THE odd-order discriminator
                expectNear (dbAt (hp ? fc * 8.0 : fc / 8.0), 0.0, 0.2, "flat passband" + at);       // deep passband ~ unity
                // monotone through the corner (no ripple/overshoot): passband -> -3 dB -> stopband
                const double near = dbAt (hp ? fc * 1.3 : fc / 1.3), far = dbAt (hp ? fc / 1.3 : fc * 1.3);
                expectTrue (near > -3.0103 && far < -3.0103, "monotone through fc, no ripple" + at);
            }
    }

    group ("Notch variable order: designBand slope->sections mapping; default slope 12 == legacy single notch");
    {
        BandParams p; p.on = true; p.type = FilterType::Notch; p.lane (Lane::Stereo).freq = 1000.0; p.lane (Lane::Stereo).Q = 4.0;
        const int expect[][2] = { {6,1}, {12,1}, {24,2}, {36,3}, {48,4}, {72,6}, {96,8} };   // order=slope/6; sections=ceil(order/2)
        for (auto& e : expect) { p.lane (Lane::Stereo).slope = e[0]; const auto d = designBand (p, fs);
            expectTrue (d.n == e[1], "slope " + std::to_string (e[0]) + " -> " + std::to_string (e[1]) + " sections"); }

        p.lane (Lane::Stereo).slope = 12; const auto d12 = designBand (p, fs);                 // the default must be the pre-existing notch, bit-for-bit
        const auto ref = matched::notch (1000.0, fs, 4.0);
        expectTrue (d12.n == 1 && d12.sec[0].b0 == ref.b0 && d12.sec[0].b1 == ref.b1 && d12.sec[0].b2 == ref.b2
                 && d12.sec[0].a1 == ref.a1 && d12.sec[0].a2 == ref.a2, "slope 12 == matched::notch bit-for-bit (no session drift)");

        p.lane (Lane::Stereo).slope = 96; p.swept = true; const auto ds = designBand (p, fs);  // swept notch stays the single fallback (SVF = one stage)
        expectTrue (ds.n == 1, "swept notch ignores slope (single fallback)");
    }

    group ("EqBand Notch variable slope: audio == analytic; in-band attenuation grows with slope");
    {
        const double f0 = 1000.0, fdet = f0 * std::pow (2.0, 1.0 / 6.0);   // +1/6 octave: inside the -3 dB band
        auto meas = [&] (int slope)
        {
            EqBand band; band.prepare (fs, 1);
            BandParams p; p.on = true; p.type = FilterType::Notch; p.lane (Lane::Stereo).freq = f0; p.lane (Lane::Stereo).Q = 3.0; p.lane (Lane::Stereo).slope = slope;
            band.setParams (p);
            const double a = sineGainDb ([&] (float* const* ch, int nc, int n) { band.processBlock (ch, nc, n); }, fdet, fs);
            const double r = 20.0 * std::log10 (std::abs (band.response (2.0 * kPi * fdet / fs)));
            return std::pair<double, double> { a, r };
        };
        double prev = 1e9;
        for (int slope : { 12, 24, 48, 96 })
        {
            const auto ar = meas (slope);
            std::printf ("      notch slope=%2d @+1/6oct: audio=%.2f  analytic=%.2f dB\n", slope, ar.first, ar.second);
            expectNear (ar.first, ar.second, 0.6, "audio == analytic response (slope " + std::to_string (slope) + ")");
            expectTrue (ar.first < prev - 2.0, "steeper slope attenuates more in-band (slope " + std::to_string (slope) + ")");
            prev = ar.first;
        }
        EqBand band; band.prepare (fs, 1);                                 // far from f0 the high-order notch is transparent
        BandParams p; p.on = true; p.type = FilterType::Notch; p.lane (Lane::Stereo).freq = f0; p.lane (Lane::Stereo).Q = 3.0; p.lane (Lane::Stereo).slope = 96;
        band.setParams (p);
        const double far = sineGainDb ([&] (float* const* ch, int nc, int n) { band.processBlock (ch, nc, n); }, 250.0, fs);
        expectNear (far, 0.0, 0.4, "slope-96 notch: unity two octaves below f0");
    }

    group ("RT-safety: a moving 8-section (slope-96) Notch processes with ZERO heap allocations");
    {
        EqBand band; band.prepare (fs, 2);
        BandParams p; p.on = true; p.type = FilterType::Notch; p.lane (Lane::Stereo).freq = 1000.0; p.lane (Lane::Stereo).Q = 3.0; p.lane (Lane::Stereo).slope = 96;
        band.setParams (p);
        const int n = 128; std::vector<float> L ((size_t) n), R ((size_t) n);
        float* ch[2] = { L.data(), R.data() };
        auto fill = [&] (int b) { for (int i = 0; i < n; ++i) {
            L[(size_t) i] = (float) (0.10 * std::sin (2.0 * kPi *  900.0 * (b * n + i) / fs));
            R[(size_t) i] = (float) (0.08 * std::sin (2.0 * kPi * 1100.0 * (b * n + i) / fs)); } };
        for (int b = 0; b < 4; ++b) { fill (b); band.processBlock (ch, 2, n); }   // warm up (allocs before the snapshot are fine)

        const long before = g_allocs.load (std::memory_order_relaxed);
        for (int b = 0; b < 300; ++b)
        {
            p.lane (Lane::Stereo).freq = 700.0 + 0.7 * b;                                      // keep the freq smoother moving -> designBand()/notchCascade run every block
            band.setParams (p);
            fill (b);
            band.processBlock (ch, 2, n);
        }
        const long after = g_allocs.load (std::memory_order_relaxed);
        std::printf ("      heap allocations during 300 moving 8-section blocks: %ld\n", after - before);
        expectTrue (after == before, "moving high-order notch: no heap allocation in the audio path");
        expectTrue (! anyNaN (L.data(), n) && ! anyNaN (R.data(), n), "finite output");
    }

    group ("Notch through a LIVE broadband signal: per-band energy is cut at f0, preserved elsewhere");
    {
        // Deterministic, portable white noise (fixed 64-bit LCG — identical bytes on every platform,
        // no std::random locale/impl variance) so the measurement is reproducible on Windows / any CPU.
        std::uint64_t st = 0x9E3779B97F4A7C15ULL;
        auto rnd = [&] { st = st * 6364136223846793005ULL + 1442695040888963407ULL;
                         return (double) (st >> 11) * (1.0 / 9007199254740992.0) * 2.0 - 1.0; };
        const int N = 1 << 15;
        std::vector<float> in ((size_t) N);
        for (int i = 0; i < N; ++i) in[(size_t) i] = (float) (0.2 * rnd());

        // Integrated band energy over [fa,fb] (sum of exact-bin Goertzel powers — averages out the
        // single-realization leakage a per-bin estimate suffers on a steep skirt; stable across builds).
        auto binPow = [&] (const std::vector<float>& x, double f) {
            const double c = 2.0 * std::cos (2.0 * kPi * f / fs); double s1 = 0.0, s2 = 0.0;
            for (int i = 0; i < N; ++i) { const double s0 = (double) x[(size_t) i] + c * s1 - s2; s2 = s1; s1 = s0; }
            return s1 * s1 + s2 * s2 - c * s1 * s2;
        };
        auto bandDb = [&] (const std::vector<float>& out, double fa, double fb) {
            const long ka = std::llround (fa * N / fs), kb = std::llround (fb * N / fs);
            double ei = 0.0, eo = 0.0;
            for (long k = ka; k <= kb; ++k) { const double f = (double) k * fs / N; ei += binPow (in, f); eo += binPow (out, f); }
            return 10.0 * std::log10 (eo / std::max (1e-300, ei));
        };
        auto runNotch = [&] (int slope) {
            EqBand band; band.prepare (fs, 1);
            BandParams p; p.on = true; p.type = FilterType::Notch; p.lane (Lane::Stereo).freq = 1000.0; p.lane (Lane::Stereo).Q = 3.0; p.lane (Lane::Stereo).slope = slope;
            band.setParams (p);
            std::vector<float> out = in;
            for (int o = 0; o < N; o += 256) { float* ch[1] = { out.data() + o }; band.processBlock (ch, 1, std::min (256, N - o)); }
            return out;
        };

        const auto o96 = runNotch (96), o24 = runNotch (24);
        const double stop96 = bandDb (o96, 917.0, 1091.0), stop24 = bandDb (o24, 917.0, 1091.0);   // ±1/6 oct around f0
        std::printf ("      live noise: stopband energy slope24=%.1f  slope96=%.1f dB\n", stop24, stop96);
        expectTrue (stop96 < -25.0,          "slope-96 strips >25 dB of live energy in the ±1/6-oct stop band");
        expectTrue (stop24 < -10.0,          "slope-24 also cuts the stop band");
        expectTrue (stop96 < stop24 - 8.0,   "higher order removes materially MORE in-band energy");
        for (const auto& band : { std::pair<double,double> { 280.0, 360.0 }, { 500.0, 700.0 }, { 3800.0, 4200.0 } })
            expectNear (bandDb (o96, band.first, band.second), 0.0, 0.2,
                        "passband energy preserved [" + std::to_string ((int) band.first) + "," + std::to_string ((int) band.second) + "]");
    }

    group ("Notch denormal flush: high-order tail decays to exact zero (no subnormal spikes on any CPU)");
    {
        EqBand band; band.prepare (fs, 1);
        BandParams p; p.on = true; p.type = FilterType::Notch; p.lane (Lane::Stereo).freq = 1000.0; p.lane (Lane::Stereo).Q = 3.0; p.lane (Lane::Stereo).slope = 96;   // 8 sections
        band.setParams (p);
        const int n = 64; std::vector<float> buf ((size_t) n, 0.0f); buf[0] = 1.0f; float* ch[1] = { buf.data() };
        band.processBlock (ch, 1, n);
        for (int b = 0; b < 300; ++b) { for (int i = 0; i < n; ++i) buf[(size_t) i] = 0.0f; band.processBlock (ch, 1, n); }
        double tail = 0.0; for (int i = 0; i < n; ++i) tail += std::fabs (buf[(size_t) i]);
        expectTrue (tail == 0.0, "8-section notch tail flushed to exact zero (no denormal residue)");
    }

    group ("Notch topology toggle (slope 96<->12 changes section count) clears stale state");
    {
        const int n = 64; std::vector<float> buf ((size_t) n, 0.0f); float* ch[1] = { buf.data() };
        EqBand band; band.prepare (fs, 1);
        BandParams hi; hi.on = true; hi.type = FilterType::Notch; hi.lane (Lane::Stereo).freq = 1000.0; hi.lane (Lane::Stereo).Q = 3.0; hi.lane (Lane::Stereo).slope = 96;
        band.setParams (hi);
        for (int i = 0; i < n; ++i) buf[(size_t) i] = 0.0f; buf[0] = 1.0f; band.processBlock (ch, 1, n);       // build 8-section state
        BandParams lo = hi; lo.lane (Lane::Stereo).slope = 12; band.setParams (lo);
        for (int i = 0; i < n; ++i) buf[(size_t) i] = 0.0f; band.processBlock (ch, 1, n);                       // -> single (resets)
        band.setParams (hi);                                                                                    // -> 8 sections again
        for (int i = 0; i < n; ++i) buf[(size_t) i] = 0.0f; band.processBlock (ch, 1, n);
        double tail = 0.0; for (int i = 0; i < n; ++i) tail += std::fabs (buf[(size_t) i]);
        expectTrue (tail == 0.0, "no stale tail after the notch section count changes");
    }

    group ("parameter safety: a high-order Notch with NaN / Inf / degenerate params stays finite");
    {
        const double inf = std::numeric_limits<double>::infinity();
        EqEngine eng; eng.prepare (fs, 256, 1);
        BandParams p; p.on = true; p.type = FilterType::Notch; p.lane (Lane::Stereo).freq = std::nan (""); p.lane (Lane::Stereo).Q = 0.0; p.lane (Lane::Stereo).gainDb = inf; p.lane (Lane::Stereo).slope = 96;
        eng.setBand (0, p);
        const int n = 256; std::vector<float> buf ((size_t) n);
        for (int i = 0; i < n; ++i) buf[(size_t) i] = 0.3f * (float) std::sin (2.0 * kPi * 1000.0 * i / fs);
        float* ch[1] = { buf.data() };
        for (int b = 0; b < 20; ++b) eng.process (ch, 1, n);
        expectTrue (! anyNaN (buf.data(), n),                  "sanitised high-order notch -> finite output");
        expectTrue (std::isfinite (eng.magnitudeDb (1000.0)),  "magnitudeDb finite for a bad-param high-order notch");
    }

    group ("Notch high order in M/S: an 8-section Mid-lane notch leaves the Side axis untouched");
    {
        EqBand b; b.prepare (fs, 2);
        BandParams p; p.on = true; p.type = FilterType::Notch;                      // {m}-only: Mid notch, Side idle -> Side (L−R) bit-exact
        p.lane (Lane::Stereo).on = false;
        p.lane (Lane::Mid).on = true; p.lane (Lane::Mid).freq = 1000.0; p.lane (Lane::Mid).Q = 3.0; p.lane (Lane::Mid).slope = 96;
        b.setParams (p);
        const int n = 200; std::vector<float> L ((size_t) n), R ((size_t) n), S0 ((size_t) n);
        for (int i = 0; i < n; ++i) { L[(size_t) i] = (float) (0.20 * std::sin (2.0 * kPi * 1122.0 * i / fs));
                                      R[(size_t) i] = (float) (0.15 * std::sin (2.0 * kPi * 1122.0 * i / fs + 0.7));
                                      S0[(size_t) i] = 0.5f * (L[(size_t) i] - R[(size_t) i]); }
        float* ch[2] = { L.data(), R.data() }; b.processBlock (ch, 2, n);
        double sErr = 0.0; for (int i = 0; i < n; ++i) sErr = std::max (sErr, (double) std::fabs (0.5f * (L[(size_t) i] - R[(size_t) i]) - S0[(size_t) i]));
        expectNear (sErr, 0.0, 1e-5, "8-section Mid notch: Side axis preserved to float precision");
    }

    // ===== per-type crash / robustness battery — EVERY FilterType, not just Bell/Notch =====
    struct TypeCase { FilterType type; const char* name; };
    const TypeCase ALL_TYPES[] = {
        { FilterType::Bell, "Bell" }, { FilterType::LowShelf, "LowShelf" }, { FilterType::HighShelf, "HighShelf" },
        { FilterType::HighPass, "HighPass" }, { FilterType::LowPass, "LowPass" }, { FilterType::BandPass, "BandPass" },
        { FilterType::Notch, "Notch" }, { FilterType::AllPass, "AllPass" }, { FilterType::Tilt, "Tilt" },
    };

    group ("every FilterType: designBand coeffs finite + every section stable across the param grid");
    {
        for (const auto& tc : ALL_TYPES)
        {
            long designs = 0, bad = 0;
            for (double sr : { 44100.0, 48000.0, 96000.0, 192000.0 })
                for (double f0 = 10.0; f0 <= 0.49 * sr; f0 *= 1.5)
                    for (double Q : { 0.05, 0.707, 5.0, 40.0 })
                        for (double g : { -30.0, 0.0, 30.0 })
                            for (int slope : { 6, 12, 24, 96 })
                                for (int sw = 0; sw < 2; ++sw)
                                {
                                    BandParams p; p.on = true; p.type = tc.type; p.lane (Lane::Stereo).freq = f0; p.lane (Lane::Stereo).Q = Q; p.lane (Lane::Stereo).gainDb = g; p.lane (Lane::Stereo).slope = slope; p.swept = (sw == 1);
                                    const BandDesign d = designBand (p, sr); ++designs;
                                    for (int s = 0; s < d.n; ++s)
                                    {
                                        const auto& c = d.sec[s];
                                        const bool fin = std::isfinite (c.b0) && std::isfinite (c.b1) && std::isfinite (c.b2)
                                                      && std::isfinite (c.a1) && std::isfinite (c.a2);
                                        if (! fin || ! c.isStable()) ++bad;
                                    }
                                }
            expectTrue (bad == 0, std::string (tc.name) + ": coeffs finite + stable across the grid (" + std::to_string (designs) + " designs)");
        }
    }

    group ("every FilterType: adversarial NaN / Inf / degenerate params -> finite output, no NaN");
    {
        const double inf = std::numeric_limits<double>::infinity(), qnan = std::nan ("");
        struct Bad { double f, q, g; };
        const Bad hostile[] = { { qnan, qnan, qnan }, { inf, inf, inf }, { -inf, -inf, -inf }, { 0.0, 0.0, 0.0 },
                                { -1000.0, -5.0, -1e9 }, { 1e12, 1e12, 1e12 }, { qnan, 0.0, inf }, { 60000.0, 40.0, 30.0 },
                                { 1e12, 0.0, 1e9 }, { 1e12, -1.0, 1e9 } };   // -> f≈Nyquist, Q≈min, +30 dB: the resonant-shelf balloon corner
        for (const auto& tc : ALL_TYPES)
        {
            int nanOut = 0, magBad = 0;
            for (const auto& h : hostile)
                for (int sw = 0; sw < 2; ++sw)
                {
                    EqEngine eng; eng.prepare (fs, 128, 2);
                    BandParams p; p.on = true; p.type = tc.type; p.lane (Lane::Stereo).freq = h.f; p.lane (Lane::Stereo).Q = h.q; p.lane (Lane::Stereo).gainDb = h.g; p.lane (Lane::Stereo).slope = 96; p.swept = (sw == 1);
                    eng.setBand (0, p);
                    const int n = 128; std::vector<float> L ((size_t) n), R ((size_t) n);
                    for (int i = 0; i < n; ++i) { L[(size_t) i] = 0.95f * (float) std::sin (2.0 * kPi * 1000.0 * i / fs);
                                                  R[(size_t) i] = 0.80f * (float) std::sin (2.0 * kPi * 3000.0 * i / fs); }
                    float* ch[2] = { L.data(), R.data() };
                    for (int blk = 0; blk < 8; ++blk) eng.process (ch, 2, n);
                    if (anyNaN (L.data(), n) || anyNaN (R.data(), n)) ++nanOut;
                    if (! std::isfinite (eng.magnitudeDb (1000.0))) ++magBad;
                }
            expectTrue (nanOut == 0, std::string (tc.name) + ": sanitised params -> finite output (no NaN)");
            expectTrue (magBad == 0, std::string (tc.name) + ": magnitudeDb finite under hostile params");
        }
    }

    group ("every FilterType: denormal tail flushes to exact zero (matched + swept)");
    {
        for (const auto& tc : ALL_TYPES)
            for (int sw = 0; sw < 2; ++sw)
            {
                EqBand band; band.prepare (fs, 2);
                BandParams p; p.on = true; p.type = tc.type; p.lane (Lane::Stereo).freq = 1000.0; p.lane (Lane::Stereo).Q = 4.0; p.lane (Lane::Stereo).gainDb = 12.0; p.lane (Lane::Stereo).slope = 96; p.swept = (sw == 1);
                band.setParams (p);
                const int n = 64; std::vector<float> L ((size_t) n, 0.0f), R ((size_t) n, 0.0f); L[0] = 1.0f; R[0] = 1.0f;
                float* ch[2] = { L.data(), R.data() };
                band.processBlock (ch, 2, n);
                for (int b = 0; b < 400; ++b) { for (int i = 0; i < n; ++i) { L[(size_t) i] = 0.0f; R[(size_t) i] = 0.0f; } band.processBlock (ch, 2, n); }
                double tail = 0.0; for (int i = 0; i < n; ++i) tail += std::fabs (L[(size_t) i]) + std::fabs (R[(size_t) i]);
                expectTrue (tail == 0.0, std::string (tc.name) + (sw ? " (swept)" : " (matched)") + ": denormal tail flushed to exact zero");
            }
    }

    group ("bypass: a bypassed band is transparent");
    {
        EqBand band; band.prepare (fs, 1);
        BandParams p; p.on = true; p.bypass = true; p.type = FilterType::Bell; p.lane (Lane::Stereo).freq = 1000.0; p.lane (Lane::Stereo).Q = 2.0; p.lane (Lane::Stereo).gainDb = 12.0;
        band.setParams (p);
        const double g = sineGainDb ([&] (float* const* ch, int nc, int n) { band.processBlock (ch, nc, n); }, 1000.0, fs);
        expectNear (g, 0.0, 0.05, "bypassed band passes audio at 0 dB");
    }

    group ("turning a band off resets state (no stale tail on re-enable)");
    {
        EqBand band; band.prepare (fs, 1);
        BandParams on; on.on = true; on.type = FilterType::Bell; on.lane (Lane::Stereo).freq = 1000.0; on.lane (Lane::Stereo).Q = 2.0; on.lane (Lane::Stereo).gainDb = 6.0;
        band.setParams (on);
        const int n = 64; std::vector<float> buf ((size_t) n, 0.0f); buf[0] = 1.0f; float* ch[1] = { buf.data() };
        band.processBlock (ch, 1, n);                              // build state from an impulse
        BandParams off = on; off.on = false; band.setParams (off);
        for (int i = 0; i < n; ++i) buf[(size_t) i] = 0.0f;
        band.processBlock (ch, 1, n);                              // off block -> resets state
        band.setParams (on);                                      // re-enable
        for (int i = 0; i < n; ++i) buf[(size_t) i] = 0.0f;
        band.processBlock (ch, 1, n);                              // silent in -> must be silent out
        double tail = 0.0; for (int i = 0; i < n; ++i) tail += std::fabs (buf[(size_t) i]);
        expectTrue (tail == 0.0, "re-enabled band starts from a clean state");
    }

    group ("topology toggle (swept<->static, 24<->12) clears stale state");
    {
        const int n = 64; std::vector<float> buf ((size_t) n, 0.0f); float* ch[1] = { buf.data() };

        EqBand band; band.prepare (fs, 1);
        BandParams sw; sw.on = true; sw.type = FilterType::Bell; sw.lane (Lane::Stereo).freq = 1000.0; sw.lane (Lane::Stereo).Q = 4.0; sw.lane (Lane::Stereo).gainDb = 6.0; sw.swept = true;
        band.setParams (sw);
        buf[0] = 1.0f; band.processBlock (ch, 1, n);                                          // SVF builds state
        BandParams st = sw; st.swept = false; band.setParams (st);
        for (int i = 0; i < n; ++i) buf[(size_t) i] = 0.0f; band.processBlock (ch, 1, n);     // -> static (resets)
        band.setParams (sw);                                                                 // -> swept again
        for (int i = 0; i < n; ++i) buf[(size_t) i] = 0.0f; band.processBlock (ch, 1, n);
        double t1 = 0.0; for (int i = 0; i < n; ++i) t1 += std::fabs (buf[(size_t) i]);
        expectTrue (t1 == 0.0, "no stale SVF tail after swept<->static");

        EqBand hp; hp.prepare (fs, 1);
        BandParams p24; p24.on = true; p24.type = FilterType::HighPass; p24.lane (Lane::Stereo).freq = 1000.0; p24.lane (Lane::Stereo).Q = 0.707; p24.lane (Lane::Stereo).slope = 24;
        hp.setParams (p24);
        for (int i = 0; i < n; ++i) buf[(size_t) i] = 0.0f; buf[0] = 1.0f; hp.processBlock (ch, 1, n);   // bq0+bq1 state
        BandParams p12 = p24; p12.lane (Lane::Stereo).slope = 12; hp.setParams (p12);
        for (int i = 0; i < n; ++i) buf[(size_t) i] = 0.0f; hp.processBlock (ch, 1, n);                  // -> single (resets)
        hp.setParams (p24);                                                                             // -> 24 again
        for (int i = 0; i < n; ++i) buf[(size_t) i] = 0.0f; hp.processBlock (ch, 1, n);
        double t2 = 0.0; for (int i = 0; i < n; ++i) t2 += std::fabs (buf[(size_t) i]);
        expectTrue (t2 == 0.0, "no stale bq1 tail after 24<->12 dB/oct");
    }

    group ("EqEngine: silence stays silent, no NaN (stereo)");
    {
        EqEngine eng; eng.prepare (fs, 256, 2);
        BandParams p; p.on = true; p.type = FilterType::HighShelf; p.lane (Lane::Stereo).freq = 8000.0; p.lane (Lane::Stereo).gainDb = 12.0;
        eng.setBand (0, p);
        const int n = 256;
        std::vector<float> L ((size_t) n, 0.0f), R ((size_t) n, 0.0f);
        float* ch[2] = { L.data(), R.data() };
        for (int b = 0; b < 10; ++b) eng.process (ch, 2, n);
        double energy = 0.0;
        for (int i = 0; i < n; ++i) energy += std::fabs (L[(size_t) i]) + std::fabs (R[(size_t) i]);
        expectTrue (energy == 0.0, "silence stays silent");
        expectTrue (! anyNaN (L.data(), n) && ! anyNaN (R.data(), n), "no NaN");
    }

    group ("EqEngine: 6-channel (5.1) static bell — every channel +6 dB, no cross-talk");
    {
        const int C = 6;
        EqEngine eng; eng.prepare (fs, 256, C);
        BandParams p; p.on = true; p.type = FilterType::Bell; p.lane (Lane::Stereo).freq = 1000.0; p.lane (Lane::Stereo).Q = 2.0; p.lane (Lane::Stereo).gainDb = 6.0;
        eng.setBand (0, p);
        const auto g = multiSineGainDb ([&] (float* const* ch, int nc, int n) { eng.process (ch, nc, n); }, C, 1000.0, fs);
        for (int c = 0; c < C; ++c)
            expectNear (g[(size_t) c], 6.0, 0.4, "ch " + std::to_string (c) + " static bell +6 dB");
    }

    group ("EqEngine: 8-channel (7.1) swept SVF bell — per-channel SVF state, every channel +6 dB");
    {
        const int C = 8;                                   // exercises Svf ic1[c]/ic2[c] for c = 0..7
        EqEngine eng; eng.prepare (fs, 256, C);
        BandParams p; p.on = true; p.type = FilterType::Bell; p.lane (Lane::Stereo).freq = 2000.0; p.lane (Lane::Stereo).Q = 3.0; p.lane (Lane::Stereo).gainDb = 6.0; p.swept = true;
        eng.setBand (0, p);
        const auto g = multiSineGainDb ([&] (float* const* ch, int nc, int n) { eng.process (ch, nc, n); }, C, 2000.0, fs);
        for (int c = 0; c < C; ++c)
            expectNear (g[(size_t) c], 6.0, 0.5, "ch " + std::to_string (c) + " swept bell +6 dB");
    }

    group ("EqEngine: 16-channel (max cap) silence stays silent, no NaN");
    {
        const int C = teq::kMaxChannels;                   // 16: top of the supported range
        EqEngine eng; eng.prepare (fs, 256, C);
        BandParams p; p.on = true; p.type = FilterType::HighShelf; p.lane (Lane::Stereo).freq = 6000.0; p.lane (Lane::Stereo).gainDb = 9.0;
        eng.setBand (0, p);
        const int n = 256;
        std::vector<std::vector<float>> bufs ((size_t) C, std::vector<float> ((size_t) n, 0.0f));
        std::vector<float*> ptr ((size_t) C);
        for (int c = 0; c < C; ++c) ptr[(size_t) c] = bufs[(size_t) c].data();
        for (int b = 0; b < 10; ++b) eng.process (ptr.data(), C, n);
        double energy = 0.0; bool nan = false;
        for (int c = 0; c < C; ++c) for (int i = 0; i < n; ++i)
        { energy += std::fabs (bufs[(size_t) c][(size_t) i]); nan = nan || ! std::isfinite (bufs[(size_t) c][(size_t) i]); }
        expectTrue (energy == 0.0, "16ch silence stays silent");
        expectTrue (! nan,         "16ch no NaN");
    }

    group ("EqEngine: 16-channel independence — every channel bit-identical to its own mono engine");
    {
        const int C = teq::kMaxChannels;     // 16, each fed a DISTINCT (uncorrelated) frequency
        const int n = 256;
        auto signalAt = [&] (int c, int i) { return (float) (0.2 * std::sin (2.0 * kPi * (110.0 + 37.0 * c) * i / fs)); };

        // Run C independent single-channel engines, then the one 16-channel engine on the same inputs.
        // Per-channel state must be fully separate -> bit-for-bit equal. Any cross-talk / wrong index
        // / channel-count-dependent math (incl. on the top channel 15) makes maxErr != 0.
        auto checkIndependence = [&] (bool swept, const std::string& label)
        {
            BandParams p; p.on = true; p.type = FilterType::Bell; p.lane (Lane::Stereo).freq = 1000.0; p.lane (Lane::Stereo).Q = 3.0; p.lane (Lane::Stereo).gainDb = 6.0; p.swept = swept;

            std::vector<std::vector<float>> ref ((size_t) C, std::vector<float> ((size_t) n));
            for (int c = 0; c < C; ++c)
            {
                EqEngine mono; mono.prepare (fs, n, 1); mono.setBand (0, p);
                for (int i = 0; i < n; ++i) ref[(size_t) c][(size_t) i] = signalAt (c, i);
                float* m[1] = { ref[(size_t) c].data() };
                mono.process (m, 1, n);
            }

            EqEngine eng; eng.prepare (fs, n, C); eng.setBand (0, p);
            std::vector<std::vector<float>> buf ((size_t) C, std::vector<float> ((size_t) n));
            std::vector<float*> ptr ((size_t) C);
            for (int c = 0; c < C; ++c)
            {
                for (int i = 0; i < n; ++i) buf[(size_t) c][(size_t) i] = signalAt (c, i);
                ptr[(size_t) c] = buf[(size_t) c].data();
            }
            eng.process (ptr.data(), C, n);

            double maxErr = 0.0;
            for (int c = 0; c < C; ++c) for (int i = 0; i < n; ++i)
                maxErr = std::max (maxErr, (double) std::fabs (buf[(size_t) c][(size_t) i] - ref[(size_t) c][(size_t) i]));
            expectTrue (maxErr == 0.0, label + ": all 16 channels bit-identical to an independent engine");
        };

        checkIndependence (false, "static biquad");
        checkIndependence (true,  "swept SVF");
    }

    group ("EqBand M/S dual-mode: independent Mid & Side lanes; the idle axis is preserved");
    {
        const int block = 256;
        // steady-state gain (dB) in the Mid (side=false) or Side (side=true) domain at frequency f.
        // Drive a pure-Mid (L=R) or pure-Side (L=-R) sine so the other domain is exactly zero.
        auto msProbe = [&] (EqBand& band, bool side, double f)
        {
            const double amp = 0.25, dp = 2.0 * kPi * f / fs;
            std::vector<float> L ((size_t) block), R ((size_t) block);
            double phase = 0.0;
            auto fill = [&] { for (int n = 0; n < block; ++n) {
                                  const float v = (float) (amp * std::sin (phase));
                                  phase += dp; if (phase > 2.0 * kPi) phase -= 2.0 * kPi;
                                  L[(size_t) n] = v; R[(size_t) n] = side ? -v : v; } };
            for (int b = 0; b < 400; ++b) { fill(); float* ch[2] = { L.data(), R.data() }; band.processBlock (ch, 2, block); }
            fill();
            std::vector<double> inDom ((size_t) block);
            for (int n = 0; n < block; ++n) inDom[(size_t) n] = side ? 0.5 * (L[(size_t) n] - R[(size_t) n]) : 0.5 * (L[(size_t) n] + R[(size_t) n]);
            { float* ch[2] = { L.data(), R.data() }; band.processBlock (ch, 2, block); }
            double ir = 0.0, orr = 0.0;
            for (int n = 0; n < block; ++n) { const double dm = side ? 0.5 * (L[(size_t) n] - R[(size_t) n]) : 0.5 * (L[(size_t) n] + R[(size_t) n]);
                                              ir += inDom[(size_t) n] * inDom[(size_t) n]; orr += dm * dm; }
            return 10.0 * std::log10 (orr / std::max (1e-30, ir));
        };

        { EqBand b; b.prepare (fs, 2);                                  // Mid +12 @1k, Side flat
          BandParams p; p.on = true; p.type = FilterType::Bell; p.lane (Lane::Stereo).on = false;
          p.lane (Lane::Mid).on  = true; p.lane (Lane::Mid).freq  = 1000.0; p.lane (Lane::Mid).Q  = 2.0; p.lane (Lane::Mid).gainDb  = 12.0;
          p.lane (Lane::Side).on = true; p.lane (Lane::Side).freq = 1000.0; p.lane (Lane::Side).Q = 2.0; p.lane (Lane::Side).gainDb =  0.0;
          b.setParams (p);
          expectNear (msProbe (b, false, 1000.0), 12.0, 0.5, "Mid lane +12 dB at 1k");
          expectNear (msProbe (b, true,  1000.0),  0.0, 0.1, "Side lane flat at 1k"); }

        { EqBand b; b.prepare (fs, 2);                                  // Side -9 @3k, Mid flat
          BandParams p; p.on = true; p.type = FilterType::Bell; p.lane (Lane::Stereo).on = false;
          p.lane (Lane::Mid).on  = true; p.lane (Lane::Mid).freq  = 3000.0; p.lane (Lane::Mid).Q  = 2.0; p.lane (Lane::Mid).gainDb  =  0.0;
          p.lane (Lane::Side).on = true; p.lane (Lane::Side).freq = 3000.0; p.lane (Lane::Side).Q = 2.0; p.lane (Lane::Side).gainDb = -9.0;
          b.setParams (p);
          expectNear (msProbe (b, true,  3000.0), -9.0, 0.5, "Side lane -9 dB at 3k");
          expectNear (msProbe (b, false, 3000.0),  0.0, 0.1, "Mid lane flat at 3k"); }

        { EqBand b; b.prepare (fs, 2);                                  // independent freq: Mid +6 @500, Side -6 @5k
          BandParams p; p.on = true; p.type = FilterType::Bell; p.lane (Lane::Stereo).on = false;
          p.lane (Lane::Mid).on  = true; p.lane (Lane::Mid).freq  = 500.0;  p.lane (Lane::Mid).Q  = 2.0; p.lane (Lane::Mid).gainDb  =  6.0;
          p.lane (Lane::Side).on = true; p.lane (Lane::Side).freq = 5000.0; p.lane (Lane::Side).Q = 2.0; p.lane (Lane::Side).gainDb = -6.0;
          b.setParams (p);
          expectNear (msProbe (b, false, 500.0),   6.0, 0.5, "Mid +6 at its own freq");
          expectNear (msProbe (b, true,  5000.0), -6.0, 0.5, "Side -6 at its own freq");
          expectNear (msProbe (b, false, 5000.0),  0.0, 0.5, "Mid flat at the Side freq");
          expectNear (msProbe (b, true,  500.0),   0.0, 0.5, "Side flat at the Mid freq"); }

        { EqBand b; b.prepare (fs, 2);                                  // {m}-only: Side idle -> Side axis bit-exact
          BandParams p; p.on = true; p.type = FilterType::Bell; p.lane (Lane::Stereo).on = false;
          p.lane (Lane::Mid).on = true; p.lane (Lane::Mid).freq = 1000.0; p.lane (Lane::Mid).Q = 2.0; p.lane (Lane::Mid).gainDb = 12.0;
          b.setParams (p);
          const int n = 200; std::vector<float> L ((size_t) n), R ((size_t) n), S0 ((size_t) n);
          for (int i = 0; i < n; ++i) { L[(size_t) i] = (float) (0.20 * std::sin (2.0 * kPi * 700.0 * i / fs));
                                        R[(size_t) i] = (float) (0.15 * std::sin (2.0 * kPi * 700.0 * i / fs + 0.7));
                                        S0[(size_t) i] = 0.5f * (L[(size_t) i] - R[(size_t) i]); }
          float* ch[2] = { L.data(), R.data() }; b.processBlock (ch, 2, n);
          double sErr = 0.0; for (int i = 0; i < n; ++i) sErr = std::max (sErr, (double) std::fabs (0.5f * (L[(size_t) i] - R[(size_t) i]) - S0[(size_t) i]));
          expectNear (sErr, 0.0, 1e-5, "Side idle: Side axis (L-R) preserved to float precision"); }

        { EqBand b; b.prepare (fs, 2);                                  // Mid bypassed (ghost) + Side active
          BandParams p; p.on = true; p.type = FilterType::Bell; p.lane (Lane::Stereo).on = false;
          p.lane (Lane::Mid).on  = true; p.lane (Lane::Mid).bypass = true;   // Mid lane kept but muted
          p.lane (Lane::Mid).freq  = 2000.0; p.lane (Lane::Mid).Q  = 2.0; p.lane (Lane::Mid).gainDb  = 12.0;
          p.lane (Lane::Side).on = true; p.lane (Lane::Side).freq = 2000.0; p.lane (Lane::Side).Q = 2.0; p.lane (Lane::Side).gainDb = 12.0;
          b.setParams (p);
          expectNear (msProbe (b, true,  2000.0), 12.0, 0.6, "Side processes while the Mid lane is bypassed");
          expectNear (msProbe (b, false, 2000.0),  0.0, 0.1, "Mid axis preserved while its lane is bypassed"); }

        { EqBand b; b.prepare (fs, 1);                                  // mono + {m,s}: v2 runs the ST lane ONLY -> transparent
          BandParams p; p.on = true; p.type = FilterType::Bell; p.lane (Lane::Stereo).on = false;
          p.lane (Lane::Mid).on  = true; p.lane (Lane::Mid).freq  = 1000.0; p.lane (Lane::Mid).Q  = 2.0; p.lane (Lane::Mid).gainDb  = 6.0;
          p.lane (Lane::Side).on = true; p.lane (Lane::Side).freq = 1000.0; p.lane (Lane::Side).Q = 2.0; p.lane (Lane::Side).gainDb = 6.0;
          b.setParams (p);
          const double g = sineGainDb ([&] (float* const* ch, int nc, int n) { b.processBlock (ch, nc, n); }, 1000.0, fs);
          expectNear (g, 0.0, 0.05, "mono + {m,s}: M/S lanes inert on a non-stereo bus (ST-only rule)"); }
    }

    group ("Lanes on a surround bus (nc>=3): only the ST lane runs — L/R/M/S are silently inert");
    {
        // v2 contract (LANES.md): on a non-stereo bus (mono / surround / ambisonics) ONLY the Stereo
        // lane runs; the L/R/M/S lanes need a 2-channel bus and are silently inactive. So a {m,s} point
        // is transparent on surround, while a plain {st} point applies to every channel.
        const int C = 6;
        { EqEngine eng; eng.prepare (fs, 256, C);                       // {m,s} split point: Mid +12 vs Side -12, same freq
          BandParams p; p.on = true; p.type = FilterType::Bell; p.lane (Lane::Stereo).on = false;
          p.lane (Lane::Mid).on  = true; p.lane (Lane::Mid).freq  = 1000.0; p.lane (Lane::Mid).Q  = 1.0; p.lane (Lane::Mid).gainDb  =  12.0;
          p.lane (Lane::Side).on = true; p.lane (Lane::Side).freq = 1000.0; p.lane (Lane::Side).Q = 1.0; p.lane (Lane::Side).gainDb = -12.0;
          eng.setBand (0, p);
          const auto g = multiSineGainDb ([&] (float* const* ch, int nc, int n) { eng.process (ch, nc, n); }, C, 1000.0, fs);
          for (int c = 0; c < C; ++c)
              expectNear (g[(size_t) c], 0.0, 0.02, "surround ch " + std::to_string (c) + " == transparent (M/S lanes inert)");
        }
        { EqEngine eng; eng.prepare (fs, 256, C);                       // plain {st} point: applies to every channel
          BandParams p; p.on = true; p.type = FilterType::Bell;
          p.lane (Lane::Stereo).freq = 1000.0; p.lane (Lane::Stereo).Q = 1.0; p.lane (Lane::Stereo).gainDb = 6.0;
          eng.setBand (0, p);
          const auto g = multiSineGainDb ([&] (float* const* ch, int nc, int n) { eng.process (ch, nc, n); }, C, 1000.0, fs);
          for (int c = 0; c < C; ++c)
              expectNear (g[(size_t) c], 6.0, 0.4, "surround ch " + std::to_string (c) + " == ST lane (+6) on every channel");
        }
    }

    // ============================ PR-A lane tests (T1 .. T4) ============================
    // Portable, seeded PRNG in [-1,1] (fixed 64-bit LCG — identical bytes on every platform).
    auto lcg = [] (std::uint64_t& st) noexcept
    {
        st = st * 6364136223846793005ULL + 1442695040888963407ULL;
        return (float) ((double) (st >> 11) * (1.0 / 9007199254740992.0) * 2.0 - 1.0);
    };

    group ("T1 legacy equivalence (reference-null): v2 lanes reproduce the v0.1.x paths bit-exact");
    {
        const int n = 128, blocks = 12;
        auto nullTest = [&] (const BandParams& band, bool ms, const std::string& label)
        {
            EqBand eng; eng.prepare (fs, 2); eng.setParams (band);
            refv01::Band ref; ref.prepare (fs, 2); ref.set (band);
            std::uint64_t s = 0xABCDEF12u;
            std::vector<float> eL ((size_t) n), eR ((size_t) n), rL ((size_t) n), rR ((size_t) n);
            double maxErr = 0.0;
            for (int b = 0; b < blocks; ++b)
            {
                for (int i = 0; i < n; ++i) { const float l = 0.30f * lcg (s), r = 0.25f * lcg (s);
                                              eL[(size_t) i] = rL[(size_t) i] = l; eR[(size_t) i] = rR[(size_t) i] = r; }
                float* ec[2] = { eL.data(), eR.data() }; eng.processBlock (ec, 2, n);
                float* rc[2] = { rL.data(), rR.data() }; if (ms) ref.processMS (rc, n); else ref.processStereo (rc, 2, n);
                for (int i = 0; i < n; ++i) { maxErr = std::max (maxErr, (double) std::fabs (eL[(size_t) i] - rL[(size_t) i]));
                                             maxErr = std::max (maxErr, (double) std::fabs (eR[(size_t) i] - rR[(size_t) i])); }
            }
            expectTrue (maxErr == 0.0, label);
        };

        { BandParams p; p.on = true; p.type = FilterType::Bell;                                    // (a) {st} == old Stereo path
          p.lane (Lane::Stereo).freq = 900.0; p.lane (Lane::Stereo).Q = 2.5; p.lane (Lane::Stereo).gainDb = 7.0;
          nullTest (p, false, "{st} bell == v0.1.x Stereo per-channel loop (bit-exact)"); }

        { BandParams p; p.on = true; p.type = FilterType::Bell; p.lane (Lane::Stereo).on = false;    // (b) {m,s} == old M/S fold
          p.lane (Lane::Mid).on  = true; p.lane (Lane::Mid).freq  = 1200.0; p.lane (Lane::Mid).Q  = 2.0; p.lane (Lane::Mid).gainDb  =  6.0;
          p.lane (Lane::Side).on = true; p.lane (Lane::Side).freq = 3000.0; p.lane (Lane::Side).Q = 1.5; p.lane (Lane::Side).gainDb = -8.0;
          nullTest (p, true, "{m,s} bells == v0.1.x M/S delta-fold (bit-exact)"); }

        { BandParams p; p.on = true; p.type = FilterType::Bell; p.lane (Lane::Stereo).on = false;    // (c) mid-bypassed + side-active
          p.lane (Lane::Mid).on  = true; p.lane (Lane::Mid).bypass = true;
          p.lane (Lane::Mid).freq  = 1000.0; p.lane (Lane::Mid).Q  = 2.0; p.lane (Lane::Mid).gainDb  = 10.0;
          p.lane (Lane::Side).on = true; p.lane (Lane::Side).freq = 2500.0; p.lane (Lane::Side).Q = 2.0; p.lane (Lane::Side).gainDb = -6.0;
          nullTest (p, true, "mid-bypassed + side-active == v0.1.x bypass semantics (bit-exact)"); }
    }

    group ("T2 axis purity: each lane touches only its own axis (untouched axis exact)");
    {
        const int n = 200;
        // L/R-only: the opposite channel is NEVER written -> bit-exact. M-only injects the SAME delta to
        // L and R (leaving the Side signal L-R), S-only injects +/- (leaving the Mid signal L+R) -> those
        // hold to float precision (the test's own L±R recombination rounds; the engine adds zero net to
        // that axis). See report note: 'bit-exact' in the spec is exact for L/R, float-exact for M/S.
        auto axisTest = [&] (Lane lane, double gainDb)
        {
            BandParams p; p.on = true; p.type = FilterType::Bell; p.lane (Lane::Stereo).on = false;
            p.lane (lane).on = true; p.lane (lane).freq = 1000.0; p.lane (lane).Q = 2.0; p.lane (lane).gainDb = gainDb;
            EqBand b; b.prepare (fs, 2); b.setParams (p);
            std::uint64_t s = 0x2468ACEu;
            std::vector<float> L ((size_t) n), R ((size_t) n), L0 ((size_t) n), R0 ((size_t) n);
            for (int i = 0; i < n; ++i) { L[(size_t) i] = L0[(size_t) i] = 0.30f * lcg (s); R[(size_t) i] = R0[(size_t) i] = 0.27f * lcg (s); }
            float* c[2] = { L.data(), R.data() }; b.processBlock (c, 2, n);
            std::array<double,4> d { 0, 0, 0, 0 };   // {maxΔL, maxΔR, maxΔ(L-R), maxΔ(L+R)}
            for (int i = 0; i < n; ++i)
            {
                d[0] = std::max (d[0], (double) std::fabs (L[(size_t) i] - L0[(size_t) i]));
                d[1] = std::max (d[1], (double) std::fabs (R[(size_t) i] - R0[(size_t) i]));
                d[2] = std::max (d[2], (double) std::fabs ((L[(size_t) i] - R[(size_t) i]) - (L0[(size_t) i] - R0[(size_t) i])));
                d[3] = std::max (d[3], (double) std::fabs ((L[(size_t) i] + R[(size_t) i]) - (L0[(size_t) i] + R0[(size_t) i])));
            }
            return d;
        };

        { const auto d = axisTest (Lane::Left, 9.0);  expectTrue (d[1] == 0.0, "L-only: R channel bit-exact"); expectTrue (d[0] > 0.0, "L-only: L channel changed"); }
        { const auto d = axisTest (Lane::Right, 9.0); expectTrue (d[0] == 0.0, "R-only: L channel bit-exact"); expectTrue (d[1] > 0.0, "R-only: R channel changed"); }
        { const auto d = axisTest (Lane::Mid, 12.0);  expectNear (d[2], 0.0, 1e-6, "M-only: Side axis (L-R) preserved"); expectTrue (d[3] > 0.0, "M-only: Mid axis (L+R) changed"); }
        { const auto d = axisTest (Lane::Side, 12.0); expectNear (d[3], 0.0, 1e-6, "S-only: Mid axis (L+R) preserved"); expectTrue (d[2] > 0.0, "S-only: Side axis (L-R) changed"); }
    }

    group ("T3 matrix truth: measured 2x2 stereo response == analytic matrixResponse (complex)");
    {
        const int W = 4096;                                  // analysis window; probe freqs are bin-aligned (leakage-free)
        auto column = [&] (EqEngine& eng, int inCh, double f) -> std::array<std::complex<double>, 2>
        {
            const double w = 2.0 * kPi * f / fs;
            std::vector<float> L ((size_t) W), R ((size_t) W);
            long nAbs = 0;
            auto fill = [&] { for (int i = 0; i < W; ++i) { const double v = 0.25 * std::sin (w * (double) (nAbs + i));
                                  L[(size_t) i] = (inCh == 0) ? (float) v : 0.0f; R[(size_t) i] = (inCh == 1) ? (float) v : 0.0f; } };
            for (int b = 0; b < 40; ++b) { fill(); float* c[2] = { L.data(), R.data() }; eng.process (c, 2, W); nAbs += W; }   // reach steady state
            fill();
            std::complex<double> X { 0, 0 };
            for (int i = 0; i < W; ++i) { const double s = (inCh == 0) ? (double) L[(size_t) i] : (double) R[(size_t) i]; X += s * std::polar (1.0, -w * (double) (nAbs + i)); }
            float* c[2] = { L.data(), R.data() }; eng.process (c, 2, W);
            std::complex<double> YL { 0, 0 }, YR { 0, 0 };
            for (int i = 0; i < W; ++i) { YL += (double) L[(size_t) i] * std::polar (1.0, -w * (double) (nAbs + i)); YR += (double) R[(size_t) i] * std::polar (1.0, -w * (double) (nAbs + i)); }
            return { YL / X, YR / X };
        };

        auto probe = [&] (const BandParams* bands, int nb, const std::string& label)
        {
            for (int ki : { 85, 171, 341 })                  // 996 / 2003 / 3996 Hz (= ki*fs/W -> integer cycles)
            {
                EqEngine eng; eng.prepare (fs, W, 2);
                for (int i = 0; i < nb; ++i) eng.setBand (i, bands[i]);
                const double f = (double) ki * fs / (double) W, w = 2.0 * kPi * f / fs;
                const auto cL = column (eng, 0, f);           // inject L -> {hLL, hRL}
                const auto cR = column (eng, 1, f);           // inject R -> {hLR, hRR}
                const ResponseMatrix H = matrixResponse (bands, nb, fs, w);
                const std::string at = " " + label + " @k" + std::to_string (ki);
                expectNear (std::abs (cL[0] - H.hLL), 0.0, 0.02, "hLL" + at);   // COMPLEX compare (proves polarity + order)
                expectNear (std::abs (cR[0] - H.hLR), 0.0, 0.02, "hLR" + at);
                expectNear (std::abs (cL[1] - H.hRL), 0.0, 0.02, "hRL" + at);
                expectNear (std::abs (cR[1] - H.hRR), 0.0, 0.02, "hRR" + at);
            }
        };

        std::uint64_t rng = 0xC0FFEE123ULL;
        auto uni = [&] (double lo, double hi) { rng = rng * 6364136223846793005ULL + 1442695040888963407ULL;
            return lo + (hi - lo) * ((double) (rng >> 11) * (1.0 / 9007199254740992.0)); };
        auto mk = [&] (LaneParams& lp, bool on) { lp.on = on; lp.freq = uni (300.0, 5000.0); lp.Q = uni (0.6, 3.5); lp.gainDb = uni (-8.0, 8.0); };

        { BandParams p; p.on = true; p.type = FilterType::Bell;                    // (a) one band, ALL FIVE lanes active
          mk (p.lane (Lane::Stereo), true); mk (p.lane (Lane::Left), true); mk (p.lane (Lane::Right), true);
          mk (p.lane (Lane::Mid), true); mk (p.lane (Lane::Side), true);
          probe (&p, 1, "all-five"); }

        { BandParams b[2];                                                          // (b) {L,R} then {M,S} — cross terms + order
          b[0].on = true; b[0].type = FilterType::Bell; b[0].lane (Lane::Stereo).on = false;
          mk (b[0].lane (Lane::Left), true); mk (b[0].lane (Lane::Right), true);
          b[1].on = true; b[1].type = FilterType::Bell; b[1].lane (Lane::Stereo).on = false;
          mk (b[1].lane (Lane::Mid), true); mk (b[1].lane (Lane::Side), true);
          probe (b, 2, "LR;MS"); }

        { BandParams b[3];                                                          // (c) three bands, seeded-random mixed subsets
          for (int j = 0; j < 3; ++j)
          {
              b[j].on = true; b[j].type = FilterType::Bell;
              bool any = false;
              for (Lane l : { Lane::Stereo, Lane::Left, Lane::Right, Lane::Mid, Lane::Side })
              { const bool on = uni (0.0, 1.0) < 0.5; mk (b[j].lane (l), on); any = any || on; }
              if (! any) b[j].lane (Lane::Stereo).on = true;   // guarantee at least one lane
          }
          probe (b, 3, "random"); }
    }

    group ("T4 topology / reset / cost: hard-step enable, type reset, ST-only buses, swept gating, cost-zero");
    {
        const int n = 128;
        auto fillNoise = [&] (std::uint64_t& s, std::vector<float>& L, std::vector<float>& R)
        { for (int i = 0; i < n; ++i) { L[(size_t) i] = 0.30f * lcg (s); R[(size_t) i] = 0.25f * lcg (s); } };

        // (a) enabling a lane mid-stream: HARD STEP, FRESH state, no NaN — the newly engaged lane settles
        //     to the SAME steady state as an engine that always had it (proves no stale history).
        {
            BandParams full; full.on = true; full.type = FilterType::Bell;
            full.lane (Lane::Stereo).freq = 1000.0; full.lane (Lane::Stereo).Q = 2.0; full.lane (Lane::Stereo).gainDb = 6.0;
            full.lane (Lane::Side).on = true; full.lane (Lane::Side).freq = 3000.0; full.lane (Lane::Side).Q = 2.0; full.lane (Lane::Side).gainDb = 9.0;
            BandParams stOnly = full; stOnly.lane (Lane::Side).on = false;

            EqBand A; A.prepare (fs, 2); A.setParams (full);
            EqBand B; B.prepare (fs, 2); B.setParams (stOnly);
            std::uint64_t sa = 0x9u, sb = 0x9u;                              // identical streams
            std::vector<float> aL ((size_t) n), aR ((size_t) n), bL ((size_t) n), bR ((size_t) n);
            bool nan = false;
            for (int blk = 0; blk < 40; ++blk)
            {
                fillNoise (sa, aL, aR); fillNoise (sb, bL, bR);
                if (blk == 6) B.setParams (full);                           // enable the Side lane mid-stream (hard step)
                float* ac[2] = { aL.data(), aR.data() }; A.processBlock (ac, 2, n);
                float* bc[2] = { bL.data(), bR.data() }; B.processBlock (bc, 2, n);
                nan = nan || anyNaN (bL.data(), n) || anyNaN (bR.data(), n);
            }
            double conv = 0.0; for (int i = 0; i < n; ++i) { conv = std::max (conv, (double) std::fabs (aL[(size_t) i] - bL[(size_t) i]));
                                                             conv = std::max (conv, (double) std::fabs (aR[(size_t) i] - bR[(size_t) i])); }
            expectTrue (! nan, "no NaN across a mid-stream lane enable");
            expectNear (conv, 0.0, 1e-5, "enabled lane settles to the always-on steady state (fresh + hard step)");
        }

        // (b) a type change with the SAME section count still resets the columns — no stale tail.
        {
            EqBand b; b.prepare (fs, 1);
            BandParams bell; bell.on = true; bell.type = FilterType::Bell;
            bell.lane (Lane::Stereo).freq = 1000.0; bell.lane (Lane::Stereo).Q = 4.0; bell.lane (Lane::Stereo).gainDb = 12.0;
            b.setParams (bell);
            std::vector<float> buf ((size_t) n, 0.0f); buf[0] = 1.0f; float* c[1] = { buf.data() };
            b.processBlock (c, 1, n);                                        // build Bell state
            BandParams bp = bell; bp.type = FilterType::BandPass;            // Bell -> BandPass: both 1 section
            b.setParams (bp);
            for (int i = 0; i < n; ++i) buf[(size_t) i] = 0.0f; b.processBlock (c, 1, n);   // silence in -> resets -> silent out
            double tail = 0.0; for (int i = 0; i < n; ++i) tail += std::fabs (buf[(size_t) i]);
            expectTrue (tail == 0.0, "type change (Bell->BandPass, same section count) resets state");
        }

        // (c) a non-stereo bus runs the ST lane ONLY — L/R/M/S are provably inert (transparent).
        {
            EqBand b; b.prepare (fs, 1);
            BandParams p; p.on = true; p.type = FilterType::Bell; p.lane (Lane::Stereo).on = false;
            p.lane (Lane::Left).on = true; p.lane (Lane::Left).gainDb = 12.0;
            p.lane (Lane::Mid).on  = true; p.lane (Lane::Mid).gainDb  = 12.0;
            b.setParams (p);
            const double g = sineGainDb ([&] (float* const* ch, int nc, int nn) { b.processBlock (ch, nc, nn); }, 1000.0, fs);
            expectNear (g, 0.0, 0.05, "true-mono: non-ST lanes inert (transparent)");
        }

        // (d) swept is honored ONLY in the single-ST config: split the point (add a flat Mid lane) and the
        //     ST lane must fall back to the matched 24 dB/oct cascade (not the single 12 dB/oct SVF stage).
        {
            EqBand b; b.prepare (fs, 2);
            BandParams p; p.on = true; p.type = FilterType::HighPass; p.swept = true;
            p.lane (Lane::Stereo).freq = 1000.0; p.lane (Lane::Stereo).Q = 0.707; p.lane (Lane::Stereo).slope = 24;
            p.lane (Lane::Mid).on = true; p.lane (Lane::Mid).freq = 1000.0; p.lane (Lane::Mid).gainDb = 0.0;   // flat Mid -> splits the point
            b.setParams (p);
            std::vector<float> L ((size_t) n, 0.0f), R ((size_t) n, 0.0f); float* c[2] = { L.data(), R.data() }; b.processBlock (c, 2, n);
            const double r250 = 20.0 * std::log10 (std::abs (b.response (2.0 * kPi * 250.0 / fs, Axis::Mid)));
            const double r500 = 20.0 * std::log10 (std::abs (b.response (2.0 * kPi * 500.0 / fs, Axis::Mid)));
            expectTrue ((r500 - r250) > 18.0, "split swept HP24: ST uses the matched 24 dB/oct cascade (swept gated off)");
        }

        // (e) COST-ZERO disabled-lane writes: churning a DISABLED lane's params every block leaves the
        //     output BIT-IDENTICAL to never touching it (no design, snapped smoother, no recompute).
        //     Chosen assertion: output-invariance (the strongest honest observable; the material-change
        //     gate additionally suppresses the recompute).
        {
            auto run = [&] (bool churn)
            {
                EqBand b; b.prepare (fs, 2);
                BandParams p; p.on = true; p.type = FilterType::Bell;
                p.lane (Lane::Stereo).freq = 1000.0; p.lane (Lane::Stereo).Q = 2.0; p.lane (Lane::Stereo).gainDb = 6.0;
                p.lane (Lane::Side).on = false; p.lane (Lane::Side).freq = 2000.0;   // present but DISABLED
                b.setParams (p);
                std::uint64_t s = 0x5EEDu; std::vector<float> L ((size_t) n), R ((size_t) n), out;
                for (int blk = 0; blk < 16; ++blk)
                {
                    fillNoise (s, L, R);
                    if (churn) { p.lane (Lane::Side).freq = 500.0 + 13.0 * blk; b.setParams (p); }   // automate the DISABLED lane
                    float* c[2] = { L.data(), R.data() }; b.processBlock (c, 2, n);
                    for (int i = 0; i < n; ++i) { out.push_back (L[(size_t) i]); out.push_back (R[(size_t) i]); }
                }
                return out;
            };
            const auto a = run (false), bb = run (true);
            double d = 0.0; for (size_t i = 0; i < a.size(); ++i) d = std::max (d, (double) std::fabs (a[i] - bb[i]));
            expectTrue (d == 0.0, "disabled-lane automation is inert: output bit-identical (cost-zero)");
        }
    }

    group ("RT-safety (lanes): a moving all-five-lane band processes with ZERO heap allocations");
    {
        EqBand band; band.prepare (fs, 2);
        BandParams p; p.on = true; p.type = FilterType::Bell;
        for (const Lane l : { Lane::Stereo, Lane::Left, Lane::Right, Lane::Mid, Lane::Side })
        { p.lane (l).on = true; p.lane (l).freq = 800.0; p.lane (l).Q = 2.0; p.lane (l).gainDb = 4.0; }
        band.setParams (p);
        const int nn = 128; std::vector<float> L ((size_t) nn), R ((size_t) nn); float* c[2] = { L.data(), R.data() };
        auto fill = [&] (int b) { for (int i = 0; i < nn; ++i) {
            L[(size_t) i] = (float) (0.10 * std::sin (2.0 * kPi *  700.0 * (b * nn + i) / fs));
            R[(size_t) i] = (float) (0.08 * std::sin (2.0 * kPi * 1100.0 * (b * nn + i) / fs)); } };
        for (int b = 0; b < 4; ++b) { fill (b); band.processBlock (c, 2, nn); }   // warm up (allocs before the snapshot are fine)

        const long before = g_allocs.load (std::memory_order_relaxed);
        for (int b = 0; b < 300; ++b)
        {
            for (const Lane l : { Lane::Stereo, Lane::Left, Lane::Right, Lane::Mid, Lane::Side })
                p.lane (l).freq = 600.0 + 0.5 * b + 5.0 * (double) (int) l;    // keep every lane's smoother moving
            band.setParams (p);
            fill (b);
            band.processBlock (c, 2, nn);
        }
        const long after = g_allocs.load (std::memory_order_relaxed);
        std::printf ("      heap allocations during 300 moving all-five-lane blocks: %ld\n", after - before);
        expectTrue (after == before, "moving 5-lane band: no heap allocation in the audio path");
        expectTrue (! anyNaN (L.data(), nn) && ! anyNaN (R.data(), nn), "finite output");
    }
}

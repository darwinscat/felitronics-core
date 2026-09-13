// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// JUCE-free self-tests for the offline Kaiser IR resampler: unity DC gain, output length, passband
// amplitude preservation, and stopband rejection (anti-aliasing on downsample) >= 55 dB. And what lies
// OUTSIDE the input (P67): zeros in front of or behind an IR change nothing but a whole-sample shift, the
// edge taps are the specification recomputed independently, a cabinet's top octave keeps its own level,
// the length never rounds to zero, and a rate or a length the arithmetic cannot carry is refused.

#include <felitronics_test.h>
#include <felitronics/convolution/IrResampler.h>
#include <felitronics/core/Math.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdio>
#include <limits>
#include <string>
#include <utility>
#include <vector>

using namespace felitronics;

static double rmsRange (const std::vector<float>& v, int from, int to)
{
    double s = 0.0; int c = 0;
    for (int i = from; i < to && i < (int) v.size(); ++i) { s += (double) v[i] * v[i]; ++c; }
    return c ? std::sqrt (s / c) : 0.0;
}

namespace
{
struct Biquad { double b0, b1, b2, a1, a2; };

enum class Shape { highPass, lowPass, peak };

Biquad rbj (Shape shape, double f0, double q, double fs, double gainDb = 0.0)
{
    const double w0 = 2.0 * core::kPi * f0 / fs, cw = std::cos (w0), alpha = std::sin (w0) / (2.0 * q);
    const double A = std::pow (10.0, gainDb / 40.0);
    double b0 = 0, b1 = 0, b2 = 0, a0 = 0, a1 = -2.0 * cw, a2 = 0;
    switch (shape)
    {
        case Shape::highPass: b0 = 0.5 * (1.0 + cw); b1 = -(1.0 + cw); b2 = b0; a0 = 1.0 + alpha; a2 = 1.0 - alpha; break;
        case Shape::lowPass:  b0 = 0.5 * (1.0 - cw); b1 = 1.0 - cw;    b2 = b0; a0 = 1.0 + alpha; a2 = 1.0 - alpha; break;
        case Shape::peak:     b0 = 1.0 + alpha * A;  b1 = -2.0 * cw;   b2 = 1.0 - alpha * A;
                              a0 = 1.0 + alpha / A;  a2 = 1.0 - alpha / A; break;
    }
    return { b0 / a0, b1 / a0, b2 / a0, a1 / a0, a2 / a0 };
}

void filterInPlace (std::vector<double>& x, const Biquad& q)
{
    double z1 = 0.0, z2 = 0.0;
    for (double& v : x)
    {
        const double y = q.b0 * v + z1;
        z1 = q.b1 * v - q.a1 * y + z2;
        z2 = q.b2 * v - q.a2 * y;
        v = y;
    }
}

// A cabinet-like IR at 48 kHz, in the shape the loader meets: its leading silence trimmed, so it starts AT
// sample 0 and rises to its peak a few samples in. An impulse through a 70 Hz high-pass, a low bump, a
// presence peak and an 8th-order Butterworth low-pass at 6 kHz, plus a gentle 14 kHz branch 40 dB down —
// the floor a real cabinet keeps over the top octave, which is where the edge defect showed. NOT tapered at
// the start: the onset IS the feature under test. The end decays past 1e-12 on its own.
std::vector<float> cabinetLikeIr (int length)
{
    constexpr double fs = 48000.0;
    std::vector<double> steep ((std::size_t) length, 0.0), floorBranch ((std::size_t) length, 0.0);
    steep[0] = floorBranch[0] = 1.0;
    for (const Biquad& q : { rbj (Shape::highPass, 70.0, 0.707, fs), rbj (Shape::peak, 110.0, 1.0, fs, 3.0),
                             rbj (Shape::peak, 2500.0, 1.2, fs, 4.0), rbj (Shape::lowPass, 6000.0, 0.5098, fs),
                             rbj (Shape::lowPass, 6000.0, 0.6013, fs), rbj (Shape::lowPass, 6000.0, 0.9000, fs),
                             rbj (Shape::lowPass, 6000.0, 2.5629, fs) })
        filterInPlace (steep, q);
    for (const Biquad& q : { rbj (Shape::highPass, 70.0, 0.707, fs), rbj (Shape::lowPass, 14000.0, 0.707, fs) })
        filterInPlace (floorBranch, q);
    std::vector<double> sum ((std::size_t) length);
    double peak = 0.0;
    for (std::size_t i = 0; i < sum.size(); ++i)
    {
        sum[i] = steep[i] + std::pow (10.0, -40.0 / 20.0) * floorBranch[i];
        peak = std::max (peak, std::fabs (sum[i]));
    }
    std::vector<float> ir ((std::size_t) length);
    for (std::size_t i = 0; i < ir.size(); ++i) ir[i] = (float) (sum[i] / peak);
    return ir;
}

// Noise that does NOT decay: a cabinet's tail is silent, so only a signal still loud at its last sample shows
// what the far edge does.
std::vector<float> noiseIr (int length, unsigned long long seed)
{
    std::vector<float> x ((std::size_t) length);
    for (float& v : x)
    {
        seed = seed * 6364136223846793005ULL + 1442695040888963407ULL;
        v = (float) ((seed >> 40) & 0xFFFFFFULL) / 8388608.0f - 1.0f;
    }
    x[0] = 1.0f;                                                           // and loudest at the onset, like the cabinet
    return x;
}

// Band power in dB of `ir` played at `fs`: the mean of |H(f)|^2 at `points` frequencies spread evenly across
// [lo, hi], each summed directly — a DTFT at those frequencies, so there is no FFT bin grid to land on.
double bandPowerDb (const std::vector<float>& ir, double fs, double lo, double hi, int points = 128)
{
    double power = 0.0;
    for (int m = 0; m < points; ++m)
    {
        const double f = lo + (hi - lo) * ((double) m + 0.5) / (double) points;
        const auto step = std::exp (std::complex<double> (0.0, -2.0 * core::kPi * f / fs));
        std::complex<double> phase (1.0, 0.0), response (0.0, 0.0);
        for (float tap : ir) { response += (double) tap * phase; phase *= step; }
        power += std::norm (response);
    }
    return 10.0 * std::log10 (power / (double) points);
}

// resampleIr recomputed from its SPECIFICATION at one output index, sharing no code with the header: long
// double; the position written as (n + 1/2)*inSr/outSr - 1/2 instead of divided by a ratio; its own Bessel
// series; sinc written as 2fc*sin(u)/u; and an explicit zero for every tap outside the input, with the window
// sum over all 2R of them. A change to the grid, the kernel or the edge rule in the header disagrees with it.
double referenceTap (const std::vector<float>& x, double inSr, double outSr, int n,
                     int halfTaps = 32, double beta = 8.0, double cutoffScale = 0.95)
{
    constexpr long double pi = 3.141592653589793238462643383279502884L;
    const auto i0 = [] (long double v)
    {
        long double sum = 1.0L, term = 1.0L;
        const long double y = v * v / 4.0L;
        for (int k = 1; k < 400; ++k)
        {
            term *= y / ((long double) k * (long double) k);
            sum += term;
            if (term < 1.0e-22L * sum) break;
        }
        return sum;
    };
    const long double t  = ((long double) n + 0.5L) * (long double) inSr / (long double) outSr - 0.5L;
    const long double fc = 0.5L * std::min (1.0L, (long double) outSr / (long double) inSr) * (long double) cutoffScale;
    const long double c  = std::floor (t);
    long double num = 0.0L, den = 0.0L;
    for (int j = 1 - halfTaps; j <= halfTaps; ++j)
    {
        const long double k  = c + (long double) j;
        const long double u  = t - k;
        const long double arg = 2.0L * pi * fc * u;
        const long double sinc = std::fabs (u) < 1.0e-12L ? 2.0L * fc : 2.0L * fc * std::sin (arg) / arg;
        const long double r  = u / (long double) halfTaps;
        const long double w  = (r <= -1.0L || r >= 1.0L) ? 0.0L
                             : sinc * i0 ((long double) beta * std::sqrt (1.0L - r * r)) / i0 ((long double) beta);
        den += w;
        const bool inside = k >= 0.0L && k < (long double) x.size();
        num += (inside ? (long double) x[(std::size_t) k] : 0.0L) * w;
    }
    return (double) (num / den);
}
} // namespace

int main()
{
    std::printf ("felitronics::convolution resampler tests\n");

    // --- unity DC gain: a constant resamples to the same constant ---
    test::group ("Resampler unity DC gain");
    {
        std::vector<float> in (1000, 0.5f);
        auto out = convolution::resampleIr (in, 48000.0, 44100.0);
        test::ok (! out.empty(), "produced output");
        double maxErr = 0.0;
        for (int i = 50; i < (int) out.size() - 50; ++i) maxErr = std::max (maxErr, (double) std::fabs (out[(std::size_t) i] - 0.5f));
        test::ok (maxErr < 5e-3, "constant preserved (DC gain == 1)");
    }

    // --- output length ~ inLen * ratio ---
    test::group ("Resampler output length");
    {
        std::vector<float> in (4800, 0.0f);
        auto out = convolution::resampleIr (in, 48000.0, 44100.0);
        const int expect = (int) std::llround (4800.0 * 44100.0 / 48000.0);
        test::ok (std::abs ((int) out.size() - expect) <= 1, "outLen ~ inLen*ratio");
    }

    // --- passband amplitude preserved (1 kHz through 48k -> 44.1k) ---
    test::group ("Resampler passband amplitude");
    {
        const int n = 4800; const double f = 1000.0, sr = 48000.0;
        std::vector<float> in (n);
        for (int i = 0; i < n; ++i) in[(std::size_t) i] = (float) std::sin (2.0 * core::kPi * f * i / sr);
        auto out = convolution::resampleIr (in, 48000.0, 44100.0);
        const double amp = rmsRange (out, 100, (int) out.size() - 100) * std::sqrt (2.0);
        test::approx (amp, 1.0, 0.05, "1 kHz amplitude preserved across the SR change");
    }

    // --- anti-alias: a 20 kHz tone is rejected when downsampling 48k -> 24k (Nyquist 12 kHz) ---
    test::group ("Resampler anti-alias stopband (>= 55 dB)");
    {
        const int n = 9600; const double sr = 48000.0, fhi = 20000.0;
        std::vector<float> in (n);
        for (int i = 0; i < n; ++i) in[(std::size_t) i] = (float) std::sin (2.0 * core::kPi * fhi * i / sr);
        auto out = convolution::resampleIr (in, 48000.0, 24000.0);
        const double inAmp = 1.0 / std::sqrt (2.0);                       // input RMS
        const double outR  = rmsRange (out, 100, (int) out.size() - 100);
        const double rejDb = 20.0 * std::log10 ((outR > 1e-12 ? outR : 1e-12) / inAmp);
        test::ok (rejDb < -55.0, "20 kHz tone rejected >= 55 dB on 48k->24k");
    }

    // --- FALSIFICATION: a degenerate config (halfTaps=0) must not silently produce an all-zero IR ---
    test::group ("Resampler sanitizes a degenerate config");
    {
        std::vector<float> in (1000, 0.5f);
        convolution::IrResampleConfig cfg; cfg.halfTaps = 0;               // window radius 0 → empty tap loop if unguarded
        auto out = convolution::resampleIr (in, 48000.0, 44100.0, cfg);
        test::ok (! out.empty(), "produced output");
        double maxErr = 1.0;                                               // fail-closed if empty
        if (! out.empty())
        {
            maxErr = 0.0;
            for (int i = 50; i < (int) out.size() - 50; ++i) maxErr = std::max (maxErr, (double) std::fabs (out[(std::size_t) i] - 0.5f));
        }
        test::ok (maxErr < 5e-2, "halfTaps=0 clamped to a usable window (not an all-zero IR)");
    }

    // --- P67: what lies outside the input ---------------------------------------------------------------
    // The law: an IR is silence on both sides of its samples. So resampling [zeros(P), x] must equal
    // resampling x shifted by D = P*outSr/inSr output samples (P chosen to make D whole), and resampling
    // [x, zeros] must equal it exactly — same positions, same weights, and the zeros add nothing. The old
    // edge rule divided by only the part of a window that landed on the input and fails both.
    test::group ("P67 — zeros around an IR change nothing but a whole-sample shift (in front: to 1e-6, behind: exactly)");
    {
        struct Case { double inSr, outSr; int zerosInFront, halfTaps; };
        const Case cases[] { { 48000.0, 44100.0, 160, 32 }, { 44100.0, 48000.0, 147, 32 }, { 48000.0, 96000.0, 16, 32 },
                             { 96000.0, 44100.0, 320, 32 }, { 48000.0, 88200.0, 80, 32 }, { 192000.0, 48000.0, 256, 32 },
                             { 44100.0, 192000.0, 441, 32 }, { 48000.0, 44100.0, 160, 8 }, { 48000.0, 44100.0, 160, 1 },
                             { 96000.0, 48000.0, 2, 64 } };
        const std::vector<float> inputs[] { cabinetLikeIr (4096), noiseIr (700, 7), noiseIr (20, 11),
                                            std::vector<float> { 0.8f }, std::vector<float> { 1.0f, -0.6f, 0.3f } };
        const char* const names[] { "cabinet", "noise 700", "noise 20", "1 tap", "3 taps" };
        int shapes = 0, frontMisses = 0, behindMisses = 0, compared = 0;
        double worstFront[5] {};
        for (const Case& cs : cases)
            for (std::size_t input = 0; input < 5; ++input)
            {
                const auto& x = inputs[input];
                const long long scaled = (long long) cs.zerosInFront * (long long) cs.outSr;
                const int shift = (int) (scaled / (long long) cs.inSr);
                const bool whole = scaled % (long long) cs.inSr == 0;
                const int behind = 2 * cs.halfTaps + 8;                    // enough zeros that every shifted index exists
                convolution::IrResampleConfig cfg; cfg.halfTaps = cs.halfTaps;

                std::vector<float> front ((std::size_t) cs.zerosInFront + x.size() + (std::size_t) behind, 0.0f);
                std::copy (x.begin(), x.end(), front.begin() + cs.zerosInFront);
                std::vector<float> back (x.size() + (std::size_t) behind, 0.0f);
                std::copy (x.begin(), x.end(), back.begin());

                const auto plain   = convolution::resampleIr (x, cs.inSr, cs.outSr, cfg);
                const auto shifted = convolution::resampleIr (front, cs.inSr, cs.outSr, cfg);
                const auto padded  = convolution::resampleIr (back, cs.inSr, cs.outSr, cfg);
                if (! whole || plain.empty() || shifted.size() < plain.size() + (std::size_t) shift
                    || padded.size() < plain.size())
                {
                    ++shapes;
                    continue;
                }
                for (std::size_t n = 0; n < plain.size(); ++n)
                {
                    const double d = std::fabs ((double) plain[n] - (double) shifted[n + (std::size_t) shift]);
                    if (! (d <= 1.0e-9 + 1.0e-6 * std::fabs ((double) plain[n]))) ++frontMisses;
                    if (std::isfinite (d)) worstFront[input] = std::max (worstFront[input], d);
                    if (! (plain[n] == padded[n])) ++behindMisses;
                    ++compared;
                }
            }
        std::printf ("    %d samples compared, %d in-front misses, %d behind misses; worst in-front difference:",
                     compared, frontMisses, behindMisses);
        for (std::size_t input = 0; input < 5; ++input) std::printf (" %s %.2e%s", names[input], worstFront[input], input < 4 ? "," : "\n");
        test::ok (shapes == 0 && compared > 0, "every case has a whole shift and output lengths that cover it");
        test::ok (frontMisses == 0, "zeros IN FRONT of an IR shift its resample by whole samples and change nothing else");
        test::ok (behindMisses == 0, "zeros BEHIND an IR change nothing at all, to the bit");
    }

    // The shift law pins the edges against the interior, but a wrong kernel or grid passes it (a shift of a
    // wrong answer is still a shift). So the edge taps — every output whose window runs off either end — are
    // held against the specification recomputed in referenceTap, which shares no code with the header.
    test::group ("P67 — edge taps equal the specification recomputed independently (grid, kernel, zeros outside)");
    {
        const std::vector<float> inputs[] { cabinetLikeIr (1024), noiseIr (300, 5) };
        constexpr double pairs[][2] { { 48000.0, 44100.0 }, { 48000.0, 96000.0 }, { 96000.0, 44100.0 }, { 44100.0, 48000.0 } };
        int misses = 0, checked = 0;
        for (const auto& x : inputs)
            for (const auto& pr : pairs)
            {
                const auto out = convolution::resampleIr (x, pr[0], pr[1]);
                const int edge = (int) std::ceil ((2.0 * 32.0 + 4.0) * std::max (1.0, pr[1] / pr[0]));
                for (int n = 0; n < (int) out.size(); ++n)
                {
                    if (n >= edge && n < (int) out.size() - edge) continue;
                    const double ref = referenceTap (x, pr[0], pr[1], n);
                    if (! (std::fabs ((double) out[(std::size_t) n] - ref) <= 1.0e-6 * std::fabs (ref) + 1.0e-9)) ++misses;
                    ++checked;
                }
            }
        std::printf ("    %d edge taps checked, %d misses\n", checked, misses);
        test::ok (checked > 0 && misses == 0, "every edge tap matches the independently recomputed specification to 1e-6");
    }

    // And one answer worked by hand, so the denominator is pinned by arithmetic rather than by a second program.
    // One tap, ratio 2, halfTaps 2, a rectangular window (beta 0), the full band (cutoffScale 1, so fc = 1/2):
    // output 0 sits at t = -1/4, its window is k = -2..1 at distances 7/4, 3/4, -1/4, -5/4, and
    // sin(pi*u)/(pi*u) there is (2*sqrt2/pi) * (-1/7, 1/3, 1, -1/5), which sums to (2*sqrt2/pi) * 104/105.
    // Only k = 0 is in the input, so the tap is 105/104; output 1 mirrors it. The old edge rule gave exactly 1.
    test::group ("P67 — a hand-worked answer: one tap at ratio 2, halfTaps 2, rectangular, full band -> 105/104 twice");
    {
        convolution::IrResampleConfig cfg;
        cfg.halfTaps = 2; cfg.beta = 0.0; cfg.cutoffScale = 1.0;
        const auto out = convolution::resampleIr (std::vector<float> { 1.0f }, 48000.0, 96000.0, cfg);
        const double want = 105.0 / 104.0;
        test::ok (out.size() == 2 && std::fabs ((double) out[0] - want) <= 1.0e-6 && std::fabs ((double) out[1] - want) <= 1.0e-6,
                  "both outputs are 105/104 = 1.0096154 (the whole window's weight in the denominator)");
    }

    // What the edge rule is FOR, on the shape it bit: a cabinet at 48 kHz played at 44.1 and at 96 kHz. Each
    // band's power, less the ratio the value-preserving resample scales a response by, is held against the
    // cabinet's OWN response — and against the UNTRUNCATED resample (zeros on both sides, every output kept),
    // which separates the kernel from the edge. The old rule lifted this fixture's top octave by up to +16 dB
    // at 44.1 and +19 dB at 96 kHz. What remains is the pre-ringing before output sample 0, which an IR that
    // starts at sample 0 cannot keep: on this fixture at most +1.2 dB, at 20 kHz, 96 kHz — hence the band bound.
    test::group ("P67 — a cabinet's top octave keeps its own level through 48 -> 44.1 and 48 -> 96 kHz");
    {
        const auto ir = cabinetLikeIr (4096);
        const auto bandOf = [] (double centre) { return std::pair { centre * std::pow (2.0, -1.0 / 12.0), centre * std::pow (2.0, 1.0 / 12.0) }; };
        const double passband = bandPowerDb (ir, 48000.0, 1000.0, 4000.0);
        const auto [lo16, hi16] = bandOf (16000.0);
        const auto [lo18, hi18] = bandOf (18000.0);
        const double own16 = bandPowerDb (ir, 48000.0, lo16, hi16) - passband;
        const double own18 = bandPowerDb (ir, 48000.0, lo18, hi18) - passband;
        std::printf ("    fixture: 16 kHz %.1f dB, 18 kHz %.1f dB under its 1-4 kHz level\n", own16, own18);
        test::ok (own16 < -40.0 && own16 > -60.0 && own18 < -40.0 && own18 > -60.0,
                  "the fixture is cabinet-like: its 16 and 18 kHz bands sit 40-60 dB under 1-4 kHz");

        struct Rate { double outSr, topCentre; };
        for (const Rate& rt : { Rate { 44100.0, 18000.0 }, Rate { 96000.0, 20000.0 } })
        {
            constexpr int zeros = 160;                                     // 160*44100/48000 = 147, 160*2 = 320: whole shifts
            std::vector<float> padded ((std::size_t) zeros + ir.size() + (std::size_t) zeros, 0.0f);
            std::copy (ir.begin(), ir.end(), padded.begin() + zeros);
            const auto played      = convolution::resampleIr (ir, 48000.0, rt.outSr);
            const auto untruncated = convolution::resampleIr (padded, 48000.0, rt.outSr);
            const double ratioDb   = 20.0 * std::log10 (rt.outSr / 48000.0);

            bool low = true, top = true, kernel = true;
            std::string line;
            for (double centre : { 1000.0, 4000.0, 8000.0, 12000.0, 14000.0, 16000.0, 18000.0, 20000.0 })
            {
                if (centre > rt.topCentre) break;
                const auto [lo, hi] = bandOf (centre);
                const double own  = bandPowerDb (ir, 48000.0, lo, hi);
                const double dev  = bandPowerDb (played, rt.outSr, lo, hi) - ratioDb - own;
                const double full = bandPowerDb (untruncated, rt.outSr, lo, hi) - ratioDb - own;
                char cell[64];
                std::snprintf (cell, sizeof cell, " %.0fk %+.2f/%+.2f", centre / 1000.0, dev, full);
                line += cell;
                if (centre <= 8000.0) low = low && std::fabs (dev) <= 0.05;
                else                  top = top && dev >= -1.0 && dev <= 1.5;
                kernel = kernel && std::fabs (full) <= 0.05;
            }
            std::printf ("    48 -> %.1f kHz, dB vs the cabinet's own (played/untruncated):%s\n", rt.outSr / 1000.0, line.c_str());
            const std::string rate = rt.outSr > 50000.0 ? "96 kHz" : "44.1 kHz";
            test::ok (low, "at " + rate + ": 1-8 kHz within 0.05 dB of the cabinet's own response");
            test::ok (top, "at " + rate + ": the top octave is not lifted — within -1/+1.5 dB of the cabinet's own response");
            test::ok (kernel, "at " + rate + ": the untruncated resample is within 0.05 dB in every band — the kernel is clean");
        }
    }

    // JUCE's resampleImpulseResponse floored the length at one sample; this rounded a short IR to NOTHING (a
    // one-tap IR at 96 -> 44.1 kHz is 0.46 of a sample), and the loader then had nothing to publish.
    test::group ("P67 — a load is never resampled to nothing: a short IR keeps one sample at every real ratio");
    {
        struct Case { std::vector<float> x; double inSr, outSr; };
        const Case cases[] { { { 0.8f }, 96000.0, 44100.0 }, { { 0.8f }, 192000.0, 44100.0 }, { { 0.8f }, 384000.0, 8000.0 },
                             { { 0.8f }, 48000.0, 22050.0 }, { { 0.8f }, 176400.0, 48000.0 },
                             { { 0.5f, -0.25f }, 384000.0, 8000.0 }, { { 1.0f, 0.5f, 0.25f }, 192000.0, 8000.0 } };
        for (const Case& cs : cases)
        {
            const auto out = convolution::resampleIr (cs.x, cs.inSr, cs.outSr);
            const double ref = referenceTap (cs.x, cs.inSr, cs.outSr, 0);
            const bool one = out.size() == 1;
            char msg[160];
            std::snprintf (msg, sizeof msg, "%zu tap(s) at %.0f -> %.0f Hz (%.3f of a sample) resample to one sample, the specification's",
                           cs.x.size(), cs.inSr, cs.outSr, (double) cs.x.size() * cs.outSr / cs.inSr);
            test::ok (one && std::isfinite (out[0]) && std::fabs ((double) out[0]) > 0.0
                          && std::fabs ((double) out[0] - ref) <= 1.0e-6 * std::fabs (ref) + 1.0e-9, msg);
        }
    }

    // The floor of one sample is what makes a vanishing ratio reachable: output 0 of a 1e-12 ratio sits 5e11
    // input samples in, and `(int) floor(t)` there is undefined. So those, rates that are not positive finite
    // numbers, and a length past INT_MAX are REFUSED — and refused before the cast, which is the order that
    // matters: 2^32 + 1 samples cast to int first is exactly ONE, a plausible-looking IR made of nothing.
    test::group ("P67 — refused, not wrapped: non-finite or vanishing rates, lengths and positions int cannot address");
    {
        const std::vector<float> one { 0.8f };
        const double inf = std::numeric_limits<double>::infinity(), nan = std::numeric_limits<double>::quiet_NaN();
        test::ok (convolution::resampleIr (one, inf, 48000.0).empty() && convolution::resampleIr (one, 48000.0, inf).empty(),
                  "an infinite rate is refused, on either side");
        test::ok (convolution::resampleIr (one, nan, 48000.0).empty() && convolution::resampleIr (one, 48000.0, nan).empty(),
                  "a NaN rate is refused, on either side");
        test::ok (convolution::resampleIr (one, 1.0e300, 1.0e-300).empty(), "a ratio that underflows to zero is refused");
        test::ok (convolution::resampleIr (one, 1.0, 4294967297.0).empty(), "2^32 + 1 output samples is refused, not wrapped to one");
        test::ok (convolution::resampleIr (one, 1.0e12, 1.0).empty(), "an output position past INT_MAX is refused, not converted");
    }

    return test::report();
}

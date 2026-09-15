// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// JUCE-free self-tests for the offline Kaiser IR resampler: unity DC gain, output length, passband
// amplitude preservation, and stopband rejection (anti-aliasing on downsample) >= 55 dB. And (P68) the
// factor that turns this function's amplitude-preserving output into a convolution kernel of the same
// gain, `convolutionRateGain` — exact where the rates make it exact, and 1 for every pair it cannot be
// made from. And what lies
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
// the start: the onset IS the feature under test. At 4096 taps its end has decayed past 1e-12 on its own (the
// band test); cut to 1024 it ends at about 7e-5, which the edge-tap check wants: signal at BOTH edges.
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

// resampleIr recomputed from its SPECIFICATION at one output index, sharing no code with the header: the
// position written as (n + 1/2)*inSr/outSr - 1/2 instead of divided by a ratio; its own Bessel series; sinc
// written as 2fc*sin(u)/u; an explicit zero for every tap outside the input, with the window sum over all 2R
// of them; and both sums COMPENSATED (Neumaier), so the oracle's own rounding sits far under the 1e-6 it is
// held to. In double, not a wider type — law 9: `long double` is 53 bits on the dev machine, 64 on x86-64
// Linux and software binary128 on wasm32, so it would have bought nothing here and a libcall there. A change
// to the grid, the kernel or the edge rule in the header disagrees with it.
struct CompensatedSum
{
    double sum = 0.0, compensation = 0.0;
    void add (double v) noexcept
    {
        const double next = sum + v;
        compensation += std::fabs (sum) >= std::fabs (v) ? (sum - next) + v : (v - next) + sum;
        sum = next;
    }
    double value() const noexcept { return sum + compensation; }
};

double referenceTap (const std::vector<float>& x, double inSr, double outSr, int n,
                     int halfTaps = 32, double beta = 8.0, double cutoffScale = 0.95)
{
    const auto i0 = [] (double v)
    {
        const double y = v * v / 4.0;
        CompensatedSum series;
        series.add (1.0);
        double term = 1.0;
        for (int k = 1; k < 400; ++k)
        {
            term *= y / ((double) k * (double) k);
            series.add (term);
            if (term < 1.0e-18 * series.value()) break;
        }
        return series.value();
    };
    const double t  = ((double) n + 0.5) * inSr / outSr - 0.5;
    const double fc = 0.5 * std::min (1.0, outSr / inSr) * cutoffScale;
    const double c  = std::floor (t);
    CompensatedSum num, den;
    for (int j = 1 - halfTaps; j <= halfTaps; ++j)
    {
        const double k    = c + (double) j;
        const double u    = t - k;
        const double arg  = 2.0 * core::kPi * fc * u;
        const double sinc = std::fabs (u) < 1.0e-12 ? 2.0 * fc : 2.0 * fc * std::sin (arg) / arg;
        const double r    = u / (double) halfTaps;
        const double w    = (r <= -1.0 || r >= 1.0) ? 0.0 : sinc * i0 (beta * std::sqrt (1.0 - r * r)) / i0 (beta);
        den.add (w);
        const bool inside = k >= 0.0 && k < (double) x.size();
        num.add ((inside ? (double) x[(std::size_t) k] : 0.0) * w);
    }
    return num.value() / den.value();
}
// P70's gate: the SAME kernel with no phase memo in it — the loop exactly as it stood before the memo, so
// the two can be compared to the last bit. This is a DIFFERENTIAL check, not a specification: what an output
// SHOULD be is `referenceTap`'s business (and the hand-worked 105/104 case's), and this one answers a
// different question — whether reusing a window that was computed for an earlier output changes anything.
// The answer has to be "not one bit", because a changed tap is a changed cabinet.
std::vector<float> resampleUnmemoized (const std::vector<float>& in, double inSr, double outSr,
                                       convolution::IrResampleConfig cfg = {})
{
    std::vector<float> out;
    const int inLen = (int) in.size();
    if (inLen <= 0) return out;
    const double ratio = outSr / inSr;
    const double want = (double) inLen * ratio;
    if (! (want < (double) std::numeric_limits<int>::max())) return out;
    const int outLen = std::max (1, (int) std::llround (want));
    if (outLen > convolution::kMaxResampleSamples) return out;
    if (cfg.halfTaps > convolution::IrResampleConfig::kMaxHalfTaps) return out;
    if (cfg.beta     > convolution::IrResampleConfig::kMaxBeta)     return out;
    const int    R      = cfg.halfTaps < 1 ? 1 : cfg.halfTaps;
    const double beta   = cfg.beta >= 0.0 ? cfg.beta : 8.0;
    const double cScale = (std::isfinite (cfg.cutoffScale) && cfg.cutoffScale > 0.0 && cfg.cutoffScale <= 1.0)
                        ? cfg.cutoffScale : 0.95;
    const double fc     = 0.5 * std::min (1.0, ratio) * cScale;
    const double i0beta = convolution::detail::besselI0 (beta);
    const double tLast  = ((double) outLen - 0.5) / ratio - 0.5;
    if (! (tLast + (double) R + 2.0 < (double) std::numeric_limits<int>::max())) return out;
    out.assign ((std::size_t) outLen, 0.0f);
    for (int n = 0; n < outLen; ++n)
    {
        const double t = ((double) n + 0.5) / ratio - 0.5;
        const int    c = (int) std::floor (t);
        double acc = 0.0, wsum = 0.0;
        for (int k = c - R + 1; k <= c + R; ++k)
        {
            const double xx   = t - (double) k;
            const double sinc = (std::fabs (xx) < 1e-12) ? (2.0 * fc)
                                                         : std::sin (2.0 * core::kPi * fc * xx) / (core::kPi * xx);
            const double r    = xx / (double) R;
            const double win  = (r <= -1.0 || r >= 1.0) ? 0.0
                              : convolution::detail::besselI0 (beta * std::sqrt (1.0 - r * r)) / i0beta;
            const double w    = sinc * win;
            wsum += w;
            if (k >= 0 && k < inLen) acc += (double) in[(std::size_t) k] * w;
        }
        out[(std::size_t) n] = (float) (! core::exactlyEqual (wsum, 0.0) ? acc / wsum : 0.0);
    }
    return out;
}

// Bits, not `==`: `0.0f == -0.0f` is true and a NaN is equal to nothing, so both would be read wrongly here.
int firstDifferentBit (const std::vector<float>& a, const std::vector<float>& b)
{
    if (a.size() != b.size()) return -1;
    for (std::size_t i = 0; i < a.size(); ++i)
        if (! core::sameBits (a[i], b[i])) return (int) i;
    return -2;                                                             // -2 = identical
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
        double worst = 0.0;                                                // |difference| / (|ref| + 1e-3), for the record
        for (const auto& x : inputs)
            for (const auto& pr : pairs)
            {
                const auto out = convolution::resampleIr (x, pr[0], pr[1]);
                const int edge = (int) std::ceil ((2.0 * 32.0 + 4.0) * std::max (1.0, pr[1] / pr[0]));
                for (int n = 0; n < (int) out.size(); ++n)
                {
                    if (n >= edge && n < (int) out.size() - edge) continue;
                    const double ref = referenceTap (x, pr[0], pr[1], n);
                    const double d = std::fabs ((double) out[(std::size_t) n] - ref);
                    if (! (d <= 1.0e-6 * std::fabs (ref) + 1.0e-9)) ++misses;
                    if (std::isfinite (d)) worst = std::max (worst, d / (std::fabs (ref) + 1.0e-3));
                    ++checked;
                }
            }
        std::printf ("    %d edge taps checked, %d misses; worst |difference| / (|reference| + 1e-3) = %.2e\n", checked, misses, worst);
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

    // Rounded, not truncated, above the floor: 3 taps at 48 -> 44.1 kHz are 2.756 samples, 600 at 44.1 -> 48 are 653.06.
    test::group ("P67 — the length is inLen*ratio rounded to the nearest sample");
    {
        test::ok (convolution::resampleIr (std::vector<float> { 1.0f, 0.5f, 0.25f }, 48000.0, 44100.0).size() == 3,
                  "3 taps at 48 -> 44.1 kHz (2.756) give 3 samples");
        test::ok (convolution::resampleIr (std::vector<float> (600, 0.1f), 44100.0, 48000.0).size() == 653,
                  "600 taps at 44.1 -> 48 kHz (653.06) give 653 samples");
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
        test::ok (convolution::resampleIr (one, 0.0, 48000.0).empty() && convolution::resampleIr (one, 48000.0, 0.0).empty(),
                  "a zero rate is refused, on either side");
        test::ok (convolution::resampleIr (one, -48000.0, -44100.0).empty() && convolution::resampleIr (one, -48000.0, 44100.0).empty(),
                  "a negative rate is refused — even when the two signs make a positive ratio");
        test::ok (convolution::resampleIr (one, 1.0e300, 1.0e-300).empty(), "a ratio that underflows to zero is refused");
        test::ok (convolution::resampleIr (one, 1.0, 4294967297.0).empty(), "2^32 + 1 output samples is refused, not wrapped to one");
        test::ok (convolution::resampleIr (one, 1.0e12, 1.0).empty(), "an output position past INT_MAX is refused, not converted");
    }

    // P70 — the phase memo. A window depends on the output position only through xx = t - k, and where
    // c = floor(t) >= 2R every one of those subtractions is EXACT (Sterbenz), so xx is frac - j with no
    // rounding and two outputs with the same 64-bit frac have the same 2R weights and the same sum. The
    // header reuses them; this says it changed nothing. Audio rates are small rationals, so the reuse is
    // real (98% for 48 -> 44.1 kHz); a ratio with no period (48000 -> 44101) reuses nothing and must still
    // agree. Both paths run here.
    test::group ("P70 — the phase memo moves no bit: every output float equals the unmemoized kernel's");
    {
        struct Case { double inSr, outSr; int len; convolution::IrResampleConfig cfg; const char* what; };
        convolution::IrResampleConfig d, r1, wide, b50, narrowBand, rect, huge;
        huge.halfTaps = convolution::IrResampleConfig::kMaxHalfTaps;
        r1.halfTaps = 1;
        wide.halfTaps = 300;
        b50.beta = 50.0;
        narrowBand.cutoffScale = 1.0e-6;
        rect.halfTaps = 2; rect.beta = 0.0; rect.cutoffScale = 1.0;
        const Case cases[] {
            { 48000.0,  44100.0, 4000, d,          "48 -> 44.1 kHz (913 phases, 98% reused)" },
            { 44100.0,  48000.0, 4000, d,          "44.1 -> 48 kHz" },
            { 48000.0,  96000.0, 2000, d,          "48 -> 96 kHz (two phases)" },
            { 96000.0,  48000.0, 4000, d,          "96 -> 48 kHz (one phase)" },
            { 192000.0, 44100.0, 8000, d,          "192 -> 44.1 kHz" },
            { 44100.0, 192000.0, 3000, d,          "44.1 -> 192 kHz" },
            { 8000.0,   48000.0, 2000, d,          "8 -> 48 kHz" },
            // An ODD integer ratio puts t at M*n + (M-1)/2, so its phase is exactly +0.0 — the one key whose
            // bit pattern a wrong empty-marker would read as "this bucket is free". Nothing else here has it.
            { 16000.0,  48000.0, 700,  d,          "16 -> 48 kHz — phase exactly +0.0 (ratio 3)" },
            { 48000.0,  16000.0, 2100, d,          "48 -> 16 kHz — the same phase the other way" },
            { 9600.0,   48000.0, 400,  d,          "9.6 -> 48 kHz (ratio 5)" },
            { 48000.0,  44101.0, 3000, d,          "48000 -> 44101 Hz — a ratio with no period at all" },
            { 48000.0,  47999.5, 3000, d,          "48000 -> 47999.5 Hz — nor this one" },
            { 48000.0,  44100.0,    1, d,          "one tap" },
            { 48000.0,  44100.0,   70, d,          "70 taps — every output is an edge output" },
            { 48000.0,  44100.0, 2000, r1,         "halfTaps 1" },
            { 48000.0,  44100.0, 2000, wide,       "halfTaps 300" },
            { 48000.0,  44100.0, 2000, b50,        "beta 50 (the series still converges to 53.04)" },
            { 48000.0,  44100.0, 2000, narrowBand, "cutoffScale 1e-6" },
            { 96000.0,  48000.0,  600, rect,       "a rectangular window, full band" },
            // The two cache branches the cases above never reach. An aperiodic ratio long enough to store
            // 4096 phases without one hit trips the GIVE-UP; a 4096-radius kernel makes each row 8193
            // doubles, so the memo's budget stops it at 127 rows and every later phase meets a FULL store —
            // periodic (it keeps hitting the 127) and aperiodic (it never hits at all) both run here.
            { 48000.0,  44101.0, 6000, d,          "48000 -> 44101 Hz, 5512 outputs — the give-up fires" },
            { 48000.0,  44101.0,  300, huge,       "halfTaps 4096, aperiodic — a full store that never hits" },
            { 48000.0,  44100.0,  300, huge,       "halfTaps 4096, periodic — a full store that does" },
        };
        int checked = 0, differ = 0;
        long long floats = 0;
        for (const Case& cs : cases)
        {
            std::vector<float> x = cabinetLikeIr (cs.len);
            if (cs.len >= 3) { x[(std::size_t) (cs.len / 3)] = 0.9f; x[(std::size_t) (cs.len - 1)] = -0.4f; }
            const auto memoed = convolution::resampleIr (x, cs.inSr, cs.outSr, cs.cfg);
            const auto plain  = resampleUnmemoized (x, cs.inSr, cs.outSr, cs.cfg);
            const int at = firstDifferentBit (memoed, plain);
            ++checked; floats += (long long) plain.size();
            // A REFUSAL IS NOT A MATCH. Two empty vectors are trivially "identical", so a case that both
            // paths refuse would pass this loop while testing nothing — which is exactly what happened to
            // the two halfTaps-4096 cases when a mutant moved the ceiling to 4095.
            if (memoed.empty()) { ++differ; std::printf ("    %s: produced nothing\n", cs.what); }
            else if (at != -2) { ++differ; std::printf ("    %s: differs at output %d\n", cs.what, at); }
        }
        std::printf ("    %d rate/config cases, %lld output floats compared as bits\n", checked, floats);
        test::ok (checked == (int) (sizeof cases / sizeof cases[0]) && floats > 0 && differ == 0,
                  "every output float is bit-identical with the memo and without it");
    }

    // A NaN or an infinity in the SAMPLES is not the resampler's business to repair (P67: broken metadata
    // plays as is, and so do broken samples) — but it must travel the same road either way, and a memo that
    // reused a window across such a sample would be visible here.
    test::group ("P70 — the memo is bit-exact on inputs that are not ordinary numbers");
    {
        const float poison[] { std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::infinity(),
                               -std::numeric_limits<float>::infinity(), 1.4e-45f, 3.4e38f };
        int differ = 0, checked = 0;
        for (float bad : poison)
            for (double outSr : { 44100.0, 96000.0 })
            {
                std::vector<float> x = cabinetLikeIr (1500);
                x[700] = bad;
                const int at = firstDifferentBit (convolution::resampleIr (x, 48000.0, outSr),
                                                  resampleUnmemoized (x, 48000.0, outSr));
                ++checked;
                if (at != -2) ++differ;
            }
        test::ok (checked == 10 && differ == 0, "NaN, both infinities, a subnormal and FLT_MAX in the input change nothing");
    }

    // A FIXTURE, NOT A SWEEP — the one place the memo's precondition can be SEEN. The gate is c >= 0,
    // because that is where `frac = t - c` stops being exact: at 8 -> 48 kHz, outputs 0 (c = -1) and 6
    // (c = 0) carry the same frac to the bit, yet five of their 64 window weights differ in the last
    // place. On almost every input those last places vanish in the float32 result, which is why a memo
    // that ignored the gate passed a sixteen-case differential sweep and 32 million output samples. This
    // input is one of the ones where they do NOT vanish: it was found by trying 18.6 million random
    // windows, and with the gate removed output 6 reads 0x1.97ee8cp-27 where it should read 0x1.97ee8ep-27.
    // (What it pins everywhere is the equality below. Whether it also KILLS an ungated memo depends on the
    // platform's sin, since that decides which side of a float32 boundary the difference lands on — so it
    // is a witness on this arithmetic, and a plain differential case on any other.)
    test::group ("P70 — the witness: an input where reusing a c < 0 phase would move a float");
    {
        const std::vector<float> x {
        0x1.a5f78ap-1f, -0x1.97b31cp-1f, -0x1.6ed8fep-2f, 0x1.1035cep-1f, -0x1.e1657ep-1f, -0x1.5f58e8p-1f,
        0x1.8ddb72p-1f, -0x1.648e56p-2f, -0x1.7c821p-3f, 0x1.de1262p-2f, 0x1.9339fep-1f, -0x1.63d3f8p-2f,
        -0x1.7bcb82p-1f, 0x1.059744p-1f, -0x1.940482p-1f, -0x1.27ec28p-1f, 0x1.e8621cp-4f, -0x1.c3b676p-3f,
        0x1.a33044p-2f, -0x1.3ba08ep-4f, 0x1.236a72p-2f, 0x1.5de8b2p-1f, 0x1.3c628cp-1f, 0x1.a1c6fcp-1f,
        -0x1.42e4ap-1f, -0x1.162cd2p-2f, 0x1.76d35p-1f, 0x1.34f444p-3f, 0x1.73355p-3f, -0x1.081f24p-1f,
        0x1.1d9ef4p-5f, 0x1.eb1a14p-1f, 0x1.5fbfap-1f
        };
        const auto memoed = convolution::resampleIr (x, 8000.0, 48000.0);
        const auto plain  = resampleUnmemoized (x, 8000.0, 48000.0);
        test::ok (memoed.size() > 6 && firstDifferentBit (memoed, plain) == -2,
                  "output 6 of an 8 -> 48 kHz resample is its own window's answer, not output 0's");
    }

    // The differential group above proves the memo changes no ANSWER — and would go on proving it if the
    // memo never ran at all (a `probe()` that always misses passes every one of those cases). So the memo's
    // own contract is pinned here, directly: it hits on a repeated key, it stops storing at its capacity,
    // and it gives up after kGiveUpPhases MISSES without a hit — including when the store filled long
    // before, which is where counting insertions instead of misses left the give-up unreachable.
    test::group ("P70 — the memo's own contract: it hits, it fills, and it gives up on a ratio with no period");
    {
        const int stride = 65;                                             // 2R + 1 at the default radius
        convolution::detail::PhaseMemo memo (stride, 10000);
        std::vector<double> row ((std::size_t) stride);
        for (int i = 0; i < stride; ++i) row[(std::size_t) i] = 1.0 + (double) i;

        const double phase = 0.3125;                                       // exact, so the key is stable
        auto first = memo.probe (phase);
        test::ok (first.row == nullptr && first.storable, "an unseen phase misses, and may be stored");
        memo.keep (first, row.data());
        auto again = memo.probe (phase);
        test::ok (again.row != nullptr && again.row[0] == 1.0 && again.row[stride - 1] == (double) stride,
                  "the same phase comes back with the same row");
        test::ok (memo.probe (0.375).row == nullptr, "a different phase does not");

        convolution::detail::PhaseMemo small (stride, 4);                  // capacity four rows
        int roomy = 0;
        for (int i = 0; i < 4; ++i)
        {
            auto p = small.probe (0.125 + (double) i);
            if (p.storable) { small.keep (p, row.data()); ++roomy; }
        }
        test::ok (roomy == 4, "four rows fit");
        test::ok (! small.probe (9.5).storable, "a fifth phase finds the store full");
        test::ok (small.probe (0.125).row != nullptr, "and the four it kept still hit");

        // The give-up, with a capacity far BELOW the threshold — the case that used to be unreachable,
        // because the only place that counted was the insertion that no longer happens.
        convolution::detail::PhaseMemo futile (stride, 4);
        int stored = 0;
        bool liveBefore = true;
        for (int i = 0; i < convolution::detail::PhaseMemo::kGiveUpPhases + 10; ++i)
        {
            if (i == convolution::detail::PhaseMemo::kGiveUpPhases - 2) liveBefore = futile.live;
            auto p = futile.probe (1.0 + (double) i * 0.5);                // every phase distinct: never a hit
            if (p.storable) { futile.keep (p, row.data()); ++stored; }
        }
        test::ok (stored == 4, "it stored only its four rows");
        test::ok (liveBefore && ! futile.live,
                  "and gave up after kGiveUpPhases misses, though it had stopped storing at the fourth");
        auto after = futile.probe (1.0);
        test::ok (after.row == nullptr && ! after.storable, "a memo that has given up neither hits nor stores");

        // A memo that IS hitting never gives up, however many misses it also takes.
        convolution::detail::PhaseMemo mixed (stride, 8);
        auto seed = mixed.probe (0.5);
        mixed.keep (seed, row.data());
        for (int i = 0; i < convolution::detail::PhaseMemo::kGiveUpPhases + 10; ++i)
        {
            (void) mixed.probe (0.5);                                      // a hit
            (void) mixed.probe (100.0 + (double) i);                       // and a miss
        }
        test::ok (mixed.live, "a memo that hits keeps working, however many misses come with the hits");
    }

    // P69 — the allocation bound. The length gate only ever kept the arithmetic addressable: INT_MAX samples
    // is 8.6 GB asked of the heap in one call, on the message thread, and a file rate does not have to be
    // absurd to ask for it. (The ACCEPTING side of this boundary is not tested here on purpose: an output of
    // exactly kMaxResampleSamples is a 64 MiB allocation and a billion window evaluations. What keeps the
    // bound from being set too LOW is every other group in this file.)
    test::group ("P69 — an output past kMaxResampleSamples is refused, not allocated");
    {
        const std::vector<float> one { 0.8f };
        test::ok (convolution::kMaxResampleSamples == (1 << 24), "the bound is 16.7M samples — 64 MiB of float per channel");
        test::ok (convolution::resampleIr (one, 1.0, (double) convolution::kMaxResampleSamples + 1.0).empty(),
                  "one sample more than the bound is refused");
        // And the bound itself is SERVED — the assertion that keeps it from being set too low. Two taps and a
        // rectangular window, so the 16.7M outputs cost a fraction of a second rather than a billion Bessel
        // series; the 64 MiB it allocates is the point of the number.
        convolution::IrResampleConfig cheap;
        cheap.halfTaps = 1; cheap.beta = 0.0;
        const auto atBound = convolution::resampleIr (one, 1.0, (double) convolution::kMaxResampleSamples, cheap);
        test::ok ((int) atBound.size() == convolution::kMaxResampleSamples, "the bound itself is served, to the sample");
        test::ok (! convolution::resampleIr (one, 1.0, 1000.0).empty(), "and an ordinary ratio is not");
        test::ok (convolution::resampleIr (std::vector<float> (48000, 0.1f), 1.2, 48000.0).empty(),
                  "a 1.2 Hz file rate against a 48 kHz host (7.7 GB at 1.92e9 samples) is refused");
        test::ok (convolution::resampleIr (std::vector<float> (1000, 0.1f), 0.001, 48000.0).empty(),
                  "so is a 1 mHz one (4.8e10 samples)");
    }

    // P69 — the config. A value that is not a request is repaired; a request this implementation will not
    // serve is refused, because answering it with a smaller kernel would be a different filter than the
    // caller asked for. beta's ceiling is where the 64-term Bessel series stops meeting its OWN convergence
    // test (53.038057): past it the series returns a number that is silently wrong — 1.2e-12 relative at 64,
    // 9.3e-5 at 90 — and only overflows to inf, and the window to NaN, at about 1.36e4.
    test::group ("P69 — the config's ceilings: a kernel past them is refused, a non-request is repaired");
    {
        const std::vector<float> x (2000, 0.25f);
        const double inf = std::numeric_limits<double>::infinity(), nan = std::numeric_limits<double>::quiet_NaN();
        convolution::IrResampleConfig cfg;

        test::ok (convolution::IrResampleConfig::kMaxHalfTaps == 4096 && convolution::IrResampleConfig::kMaxBeta == 53.0,
                  "the ceilings are 4096 taps of radius and beta 53 — the numbers the header argues for");
        cfg = {}; cfg.halfTaps = convolution::IrResampleConfig::kMaxHalfTaps;
        const auto atRadiusCeiling = convolution::resampleIr (std::vector<float> (200, 0.2f), 48000.0, 44100.0, cfg);
        test::ok (! atRadiusCeiling.empty(), "the radius ceiling itself is SERVED, not refused");
        cfg = {}; cfg.halfTaps = convolution::IrResampleConfig::kMaxHalfTaps + 1;
        test::ok (convolution::resampleIr (x, 48000.0, 44100.0, cfg).empty(), "halfTaps one past the ceiling is refused");
        cfg = {}; cfg.halfTaps = 1000000;
        test::ok (convolution::resampleIr (x, 48000.0, 44100.0, cfg).empty(), "and a millionfold one is not a billion-cycle loop");
        cfg = {}; cfg.beta = 53.04;
        test::ok (convolution::resampleIr (x, 48000.0, 44100.0, cfg).empty(), "a beta past where the series converges is refused");
        cfg = {}; cfg.beta = 1.0e300;
        test::ok (convolution::resampleIr (x, 48000.0, 44100.0, cfg).empty(), "beta 1e300 — finite, and what the isfinite gate let through");
        cfg = {}; cfg.beta = inf;
        test::ok (convolution::resampleIr (x, 48000.0, 44100.0, cfg).empty(), "an infinite beta is refused by the same range");

        cfg = {}; cfg.beta = 53.0;
        auto atCeiling = convolution::resampleIr (x, 48000.0, 44100.0, cfg);
        bool allFinite = ! atCeiling.empty();
        for (float v : atCeiling) allFinite = allFinite && std::isfinite (v);
        test::ok (allFinite, "beta exactly 53 is served, and every tap of it is a number");

        cfg = {}; cfg.beta = nan;
        const auto nanBeta = convolution::resampleIr (x, 48000.0, 44100.0, cfg);
        cfg = {}; cfg.beta = -1.0;
        const auto negBeta = convolution::resampleIr (x, 48000.0, 44100.0, cfg);
        const auto plain   = convolution::resampleIr (x, 48000.0, 44100.0);
        test::ok (firstDifferentBit (nanBeta, plain) == -2 && firstDifferentBit (negBeta, plain) == -2,
                  "a NaN or negative beta is not a request, and is repaired to the default — exactly it");
    }

    // NAMED, MEASURED, AND DELIBERATELY NOT REPAIRED. `cutoffScale`'s range (0, 1] is ratified, and a
    // subnormal inside it still underflows the product: fc = 0.5 * min(1, ratio) * 5e-324 rounds to zero,
    // every weight with it, and the sum-is-zero guard then emits an all-zero IR. That is a hole — an IR of
    // silence is an answer where there is none — but repairing it means refusing a value the ratified range
    // accepts, which is its own decision and not this task's. What is pinned here is what the code does
    // TODAY, so that the hole cannot change shape unnoticed: zeros, and not NaNs.
    test::group ("P69 — a cutoffScale so small that fc underflows: all-zero, not NaN (a named hole, pinned)");
    {
        convolution::IrResampleConfig cfg;
        cfg.cutoffScale = std::numeric_limits<double>::denorm_min();       // 5e-324, inside the ratified range
        const auto out = convolution::resampleIr (std::vector<float> (300, 0.5f), 48000.0, 48000.0, cfg);
        bool allZero = ! out.empty();
        for (float v : out) allZero = allZero && core::sameBits (v, 0.0f);
        test::ok (allZero, "every tap is +0.0 — the weights vanished, and the guard did not divide by them");
    }

    // P68 — THE FACTOR ITSELF, at the door of the function that owns it. `resampleIr` above keeps a tap's
    // AMPLITUDE; a convolution's gain is a sum over taps, so it keeps neither unless the density change is
    // taken back out. That factor is one line of arithmetic, and the reason it is a named function rather
    // than a literal in the loader is that it was ALREADY written out by hand in a second product
    // (orbit-amp's CabinetIr.h) — the class of defect where two copies of the same arithmetic must agree
    // forever. These are the values a caller may rely on.
    test::group ("P68 — convolutionRateGain: the density factor, exact where it can be and total everywhere");
    {
        using convolution::convolutionRateGain;
        const double inf = std::numeric_limits<double>::infinity();
        const double nan = std::numeric_limits<double>::quiet_NaN();

        // EXACT, and spelled as literals rather than recomputed — a test that divides the same two numbers
        // the same way agrees with any formula that happens to be there, including a wrong one.
        test::ok (convolutionRateGain (48000.0, 96000.0)  == 0.5
               && convolutionRateGain (96000.0, 48000.0)  == 2.0
               && convolutionRateGain (48000.0, 192000.0) == 0.25
               && convolutionRateGain (192000.0, 48000.0) == 4.0,
                  "the power-of-two rate pairs are exact: 0.5, 2, 0.25, 4");
        test::ok (convolutionRateGain (48000.0, 48000.0) == 1.0 && convolutionRateGain (44100.0, 44100.0) == 1.0,
                  "a rate against itself is exactly 1 — the un-resampled path can never be scaled by an epsilon");
        // This next one IS the same division on both sides, so it cannot tell a wrong formula from a right
        // one — it is here to pin the SPELLING (that the factor is inSr/outSr and not, say, the tap-count
        // ratio, which differs from it by up to 1.3e-3 dB). The literals above are what carry the truth.
        test::approx (convolutionRateGain (48000.0, 44100.0), 48000.0 / 44100.0, 0.0,
                      "and 48 -> 44.1 kHz is the ratio itself, to the bit");

        // TOTAL: every double a caller can type has an answer, and the answer for a factor that is not a
        // factor is 1 — no compensation. Returning the garbage instead would silence an IR (0) or blast it
        // (inf), and refusing would drop a load whose SAMPLES are fine; P67 settled that broken metadata
        // plays as is, and this is the same rule one function further in.
        test::ok (convolutionRateGain (nan, 48000.0) == 1.0 && convolutionRateGain (48000.0, nan) == 1.0
               && convolutionRateGain (nan, nan) == 1.0,
                  "a NaN on either side (or both) is 1, not a NaN that would poison every tap");
        test::ok (convolutionRateGain (0.0, 48000.0) == 1.0 && convolutionRateGain (48000.0, 0.0) == 1.0,
                  "a zero rate is 1, not 0 (a silenced IR) and not an infinity (a blasted one)");
        test::ok (convolutionRateGain (-48000.0, 96000.0) == 1.0 && convolutionRateGain (48000.0, -96000.0) == 1.0
               && convolutionRateGain (-48000.0, -96000.0) == 1.0,
                  "a negative rate is 1 — even when the two signs would have made a positive ratio");
        test::ok (convolutionRateGain (inf, 48000.0) == 1.0 && convolutionRateGain (48000.0, inf) == 1.0
               && convolutionRateGain (-inf, 48000.0) == 1.0 && convolutionRateGain (48000.0, -inf) == 1.0,
                  "both infinities, on both sides, are 1");
        test::ok (convolutionRateGain (1.0e300, 1.0e-300) == 1.0 && convolutionRateGain (1.0e-300, 1.0e300) == 1.0,
                  "two finite rates whose ratio overflows to inf or underflows to zero are 1 as well");
    }

    return test::report();
}

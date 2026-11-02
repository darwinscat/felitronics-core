// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// codecgrid::Mp3Hybrid against the formulas of ISO/IEC 11172-3 written out directly, in double, with nothing
// shared with the class but the window table: the polyphase analysis as one 512-term sum per subband sample,
// the hybrid MDCT as the 36 x 18 cosine matrix, the alias butterflies from the c_i of table 3-B.9. The class
// folds, tabulates and runs in single precision; this file does none of that.

#include <felitronics/codecgrid/Mp3Hybrid.h>

#include "felitronics_test.h"

#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

using felitronics::codecgrid::Mp3Hybrid;
using felitronics::codecgrid::kMp3SynthesisWindow;
using felitronics::test::group;
using felitronics::test::ok;

namespace
{
constexpr double kPi = 3.14159265358979323846;

// A deterministic programme: a few partials and a noise, so that every subband carries something.
std::vector<float> programme (int n, std::uint32_t seed)
{
    std::vector<float> x ((std::size_t) n);
    std::uint32_t s = seed;
    for (int i = 0; i < n; ++i)
    {
        s = s * 1664525u + 1013904223u;
        const double noise = ((double) (s >> 8) / 8388608.0 - 1.0) * 0.2;
        x[(std::size_t) i] = (float) (0.3 * std::sin (2.0 * kPi * 0.01 * i) + 0.2 * std::sin (2.0 * kPi * 0.137 * i + 0.4)
                                      + 0.1 * std::sin (2.0 * kPi * 0.31 * i + 1.1) + noise);
    }
    return x;
}

// The standard's analysis, one subband sample: S[sb] = sum_i C[i] X[i] cos ((2 sb + 1) (i mod 64 - 16) pi / 64),
// X[i] the sample i steps back from the newest one, the newest at x[first + 511].
double directSubband (const std::vector<float>& x, int first, int sb)
{
    double acc = 0.0;
    for (int i = 0; i < 512; ++i)
    {
        const double c = kMp3SynthesisWindow[(std::size_t) i] / 32.0;
        const double sample = (double) x[(std::size_t) (first + 511 - i)];
        acc += c * sample * std::cos ((double) ((2 * sb + 1) * ((i % 64) - 16)) * kPi / 64.0);
    }
    return acc;
}

// The standard's hybrid for one granule, from 36 x 32 subband samples (double).
void directGranule (const std::vector<double>& sub /* [36][32] */, std::vector<double>& lines /* 576 */)
{
    for (int sb = 0; sb < 32; ++sb)
        for (int k = 0; k < 18; ++k)
        {
            double acc = 0.0;
            for (int n = 0; n < 36; ++n)
            {
                double v = sub[(std::size_t) (n * 32 + sb)];
                if ((sb & 1) != 0 && ((n % 18) & 1) != 0) v = -v;                    // frequency inversion
                acc += v * std::sin (kPi / 36.0 * (n + 0.5)) * std::cos (kPi / 72.0 * (2.0 * n + 1.0 + 18.0) * (2.0 * k + 1.0));
            }
            lines[(std::size_t) (sb * 18 + k)] = acc;
        }
    const double ci[8] { -0.6, -0.535, -0.33, -0.185, -0.095, -0.041, -0.0142, -0.0037 };
    for (int sb = 1; sb < 32; ++sb)
        for (int i = 0; i < 8; ++i)
        {
            const double cs = 1.0 / std::sqrt (1.0 + ci[i] * ci[i]), ca = ci[i] / std::sqrt (1.0 + ci[i] * ci[i]);
            const double lo = lines[(std::size_t) (18 * sb - 1 - i)], hi = lines[(std::size_t) (18 * sb + i)];
            lines[(std::size_t) (18 * sb - 1 - i)] = lo * cs + hi * ca;
            lines[(std::size_t) (18 * sb + i)] = hi * cs - lo * ca;
        }
}
}

int main()
{
    std::printf ("felitronics::codecgrid MP3 hybrid filterbank tests\n");
    const Mp3Hybrid bank;

    group ("the counts: how many subband samples and granules a stretch holds");
    {
        ok (Mp3Hybrid::subbandSamples (511, 0) == 0, "511 samples hold no subband sample");
        ok (Mp3Hybrid::subbandSamples (512, 0) == 1, "512 samples hold exactly one");
        ok (Mp3Hybrid::subbandSamples (512, 1) == 0, "...and none one sample later");
        ok (Mp3Hybrid::subbandSamples (543, 0) == 1 && Mp3Hybrid::subbandSamples (544, 0) == 2, "one more every 32 samples");
        ok (Mp3Hybrid::subbandSamples (544, 31) == 1, "phase 31 of 544 samples: the last full window");
        ok (Mp3Hybrid::subbandSamples (0, 0) == 0 && Mp3Hybrid::subbandSamples (10, 31) == 0, "short input: none, not a negative count");
        ok (Mp3Hybrid::granules (35, 0) == 0, "35 subband samples: no granule (it needs 36)");
        ok (Mp3Hybrid::granules (36, 0) == 1, "36: one");
        ok (Mp3Hybrid::granules (36, 1) == 0, "36 from start 1: none");
        ok (Mp3Hybrid::granules (54, 0) == 2 && Mp3Hybrid::granules (53, 0) == 1, "one more every 18");
        ok (Mp3Hybrid::granules (0, 17) == 0 && Mp3Hybrid::granules (5, 17) == 0, "fewer samples than the start: none");
    }

    group ("the polyphase analysis nulls against the 512-term sum of the standard, at every phase");
    {
        const int n = 512 + 32 * 40 + 31;
        const auto x = programme (n, 12345u);
        double worst = 0.0, peak = 0.0;
        int compared = 0;
        for (int phase = 0; phase < 32; ++phase)
        {
            const int count = Mp3Hybrid::subbandSamples (n, phase);
            std::vector<float> got ((std::size_t) count * 32);
            bank.subbands (x.data(), n, phase, got.data());
            for (int t = 0; t < count; t += 7)
                for (int sb = 0; sb < 32; ++sb)
                {
                    const double want = directSubband (x, phase + 32 * t, sb);
                    worst = std::fmax (worst, std::fabs ((double) got[(std::size_t) (t * 32 + sb)] - want));
                    peak = std::fmax (peak, std::fabs (want));
                    ++compared;
                }
        }
        std::printf ("    subbands: %d values compared, peak %.4f, worst difference %.3e (%.1f dB under the peak)\n",
                     compared, peak, worst, 20.0 * std::log10 (peak / std::fmax (worst, 1e-300)));
        ok (compared > 5000 && peak > 0.05, "the comparison ran on a signal that reaches the subbands");
        ok (std::isfinite (worst) && worst < 2.0e-6 * peak, "single-precision agreement with the direct sum (120 dB under the peak, less 6)");
    }

    group ("the hybrid nulls against the 36 x 18 matrix and the butterflies of the standard");
    {
        const int n = 512 + 32 * (36 + 18 * 3 + 17);
        const auto x = programme (n, 777u);
        double worst = 0.0, peak = 0.0;
        int compared = 0;
        for (int phase : { 0, 5, 31 })
        {
            const int count = Mp3Hybrid::subbandSamples (n, phase);
            std::vector<float> sub ((std::size_t) count * 32);
            bank.subbands (x.data(), n, phase, sub.data());
            for (int g = 0; g < 18; g += 5)
                for (int q = 0; q < Mp3Hybrid::granules (count, g); ++q)
                {
                    const int first = g + 18 * q;
                    float lines[576];
                    bank.granule (sub.data() + (std::size_t) first * 32, lines);
                    std::vector<double> in ((std::size_t) 36 * 32), want (576);
                    for (int t = 0; t < 36; ++t)
                        for (int sb = 0; sb < 32; ++sb) in[(std::size_t) (t * 32 + sb)] = directSubband (x, phase + 32 * (first + t), sb);
                    directGranule (in, want);
                    for (int i = 0; i < 576; ++i)
                    {
                        worst = std::fmax (worst, std::fabs ((double) lines[i] - want[(std::size_t) i]));
                        peak = std::fmax (peak, std::fabs (want[(std::size_t) i]));
                        ++compared;
                    }
                }
        }
        std::printf ("    lines: %d values compared, peak %.4f, worst difference %.3e (%.1f dB under the peak)\n",
                     compared, peak, worst, 20.0 * std::log10 (peak / std::fmax (worst, 1e-300)));
        ok (compared >= 576 * 12 && peak > 0.05, "the comparison ran over several phases and granule starts");
        ok (std::isfinite (worst) && worst < 4.0e-6 * peak, "single-precision agreement with the direct formulas");
    }

    group ("what the filterbank is: a tone lands in its subband, and silence stays silence");
    {
        // 1 kHz at 44.1 kHz sits in subband 1 (689..1378 Hz): its 18 lines carry the energy.
        const int n = 512 + 32 * 36;
        std::vector<float> x ((std::size_t) n);
        for (int i = 0; i < n; ++i) x[(std::size_t) i] = (float) std::sin (2.0 * kPi * 1000.0 / 44100.0 * i);
        std::vector<float> sub ((std::size_t) Mp3Hybrid::subbandSamples (n, 0) * 32);
        bank.subbands (x.data(), n, 0, sub.data());
        float lines[576];
        bank.granule (sub.data(), lines);
        double in = 0.0, total = 0.0;
        for (int i = 0; i < 576; ++i) { total += (double) lines[i] * lines[i]; if (i >= 18 && i < 36) in += (double) lines[i] * lines[i]; }
        std::printf ("    a 1 kHz tone: %.2f %% of the granule's energy in subband 1\n", 100.0 * in / total);
        ok (total > 1.0 && in / total > 0.99, "at least 99 % of a 1 kHz tone is in the 18 lines of subband 1");

        // silence in, silence out — exactly
        std::vector<float> zeros ((std::size_t) n, 0.0f), zsub (sub.size(), 1.0f);
        bank.subbands (zeros.data(), n, 0, zsub.data());
        bank.granule (zsub.data(), lines);
        bool allZero = true;
        for (float v : lines) allZero = allZero && std::fpclassify (v) == FP_ZERO;
        ok (allZero, "digital silence gives lines that are exactly zero");
    }

    return felitronics::test::report();
}

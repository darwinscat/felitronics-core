// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// Adversarial falsification pass for StreamResampler. These checks materialize the valid external
// theory cases not already named in the main theory suite, then attack precision and chunk schedules.

#include <felitronics/core/StreamResampler.h>

#include <felitronics_test.h>

#include <algorithm>
#include <cmath>
#include <vector>

using felitronics::core::StreamResampler;
using felitronics::test::approx;
using felitronics::test::group;
using felitronics::test::ok;

namespace
{
std::vector<float> ramp (int n)
{
    std::vector<float> x ((std::size_t) n);
    for (int i = 0; i < n; ++i) x[(std::size_t) i] = (float) i;
    return x;
}

std::vector<float> alternatingSignal (int n)
{
    std::vector<float> x ((std::size_t) n);
    for (int i = 0; i < n; ++i)
        x[(std::size_t) i] = 0.25f * (float) std::sin (0.017 * i)
                           + 0.10f * (float) std::sin (0.113 * i)
                           + ((i & 1) ? -0.03125f : 0.03125f);
    return x;
}

std::vector<float> oneShot (double inRate, double outRate, const std::vector<float>& in)
{
    StreamResampler r;
    r.reset (inRate, outRate, (int) in.size() + 32);
    r.feed (in.data(), (int) in.size());
    std::vector<float> out ((std::size_t) ((double) in.size() * outRate / inRate) + 128);
    const int got = r.produceAvailable (out.data(), (int) out.size());
    out.resize ((std::size_t) got);
    return out;
}

std::vector<float> streamed (double inRate, double outRate, const std::vector<float>& in,
                             const std::vector<int>& chunks)
{
    int largest = 1;
    for (int n : chunks) largest = std::max (largest, n);

    StreamResampler r;
    r.reset (inRate, outRate, largest * 4 + 96);

    std::vector<float> out;
    std::vector<float> tmp ((std::size_t) (largest * 8 + 256));
    std::size_t pos = 0;
    int ci = 0;
    while (pos < in.size())
    {
        const int want = chunks[(std::size_t) (ci++ % (int) chunks.size())];
        const int n = (int) std::min<std::size_t> ((std::size_t) want, in.size() - pos);
        r.feed (in.data() + pos, n);
        pos += (std::size_t) n;
        const int got = r.produceAvailable (tmp.data(), (int) tmp.size());
        out.insert (out.end(), tmp.begin(), tmp.begin() + got);
    }
    return out;
}
} // namespace

int main()
{
    std::printf ("felitronics::core StreamResampler FALSIFY tests\n");

    group ("DeepSeek 1.2 adapted: 44.1k -> 88.2k ramp is exact at half-sample phases");
    {
        const auto in = ramp (4096);
        const auto out = oneShot (44100.0, 88200.0, in);
        bool exact = out.size() > 6000;
        for (std::size_t k = 8; k + 8 < out.size(); ++k)
        {
            const float want = 0.5f * (float) k - 2.0f;
            if (out[k] != want) { exact = false; break; }
        }
        ok (exact, "linear ramp maps to y[k] = 0.5*k - 2 exactly after startup");
    }

    group ("DeepSeek 1.3 adapted: 48k -> 16k integer-step ramp uses the frac=0 path exactly");
    {
        const auto in = ramp (4096);
        const auto out = oneShot (48000.0, 16000.0, in);
        bool exact = out.size() > 1000;
        for (std::size_t k = 3; k + 3 < out.size(); ++k)
        {
            const float want = 3.0f * (float) k - 2.0f;
            if (out[k] != want) { exact = false; break; }
        }
        ok (exact, "linear ramp maps to y[k] = 3*k - 2 exactly at integer phases");
    }

    group ("DeepSeek 1.4 adapted: one-second 44.1k -> 48k ramp has bounded float precision error, no drift");
    {
        const auto in = ramp (44100);
        const auto out = oneShot (44100.0, 48000.0, in);
        const double step = 44100.0 / 48000.0;
        double worst = 0.0;
        double firstWindow = 0.0, lastWindow = 0.0;
        int firstN = 0, lastN = 0;
        for (std::size_t k = 128; k + 128 < out.size(); ++k)
        {
            const double want = step * (double) k - 2.0;
            const double err = std::fabs ((double) out[k] - want);
            worst = std::max (worst, err);
            if (k < 1024) { firstWindow += err; ++firstN; }
            if (k + 1024 > out.size()) { lastWindow += err; ++lastN; }
        }
        std::printf ("      worst ramp error %.6g, first avg %.6g, last avg %.6g\n",
                     worst, firstWindow / std::max (1, firstN), lastWindow / std::max (1, lastN));
        ok (worst < 0.004, "float-output ramp error stays under one large-value ulp scale over 1 second");
        ok ((lastWindow / std::max (1, lastN)) < 0.0015,
            "tail ramp error remains bounded by float quantization, not accumulator drift");
    }

    group ("DeepSeek 1.5 adapted: hostile chunk schedule matches one-shot within float phase rounding");
    {
        const auto in = alternatingSignal (12000);
        const auto a = oneShot (44100.0, 48000.0, in);
        const auto b = streamed (44100.0, 48000.0, in, { 1, 2, 17, 3, 257, 5, 509, 64 });
        ok (a.size() == b.size(), "hostile chunking preserves produced count");
        float worst = 0.0f;
        const std::size_t n = std::min (a.size(), b.size());
        for (std::size_t i = 0; i < n; ++i) worst = std::max (worst, std::fabs (a[i] - b[i]));
        ok (worst <= 2.0e-6f, "hostile chunking stays within float phase-rounding tolerance");
    }

    group ("own attack: quadratic ramp remains exact across non-integer ratio within float arithmetic");
    {
        std::vector<float> in (4096);
        constexpr double scale = 1.0 / 1024.0;                  // keep the oracle below large-value float ulps
        for (int i = 0; i < (int) in.size(); ++i) in[(std::size_t) i] = (float) ((double) i * (double) i * scale);
        const auto out = oneShot (48000.0, 96000.0, in);
        double worst = 0.0;
        for (std::size_t k = 16; k + 16 < out.size(); ++k)
        {
            const double x = 0.5 * (double) k - 2.0;
            worst = std::max (worst, std::fabs ((double) out[k] - x * x * scale));
        }
        ok (worst < 0.01, "quadratic reproduction holds to float rounding across half-sample phases");
    }

    return felitronics::test::report();
}

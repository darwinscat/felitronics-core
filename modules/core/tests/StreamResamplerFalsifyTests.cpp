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

    // 🔴 WHAT CHANGED HERE, because it is a contract change and not a slackened bound. The cubic was an
    // INTERPOLATING kernel — h(0) = 1, h(n != 0) = 0 — so it passed through the input samples and
    // reproduced a linear ramp BIT-EXACTLY at every phase. A Kaiser-windowed sinc with a 0.99 cutoff is
    // an APPROXIMATING kernel: h(0) = 0.99·window != 1, so it does not pass through the samples and
    // exact polynomial reproduction is gone for good.
    //
    // What replaces it is not weaker in the way it looks. A symmetric, unit-sum kernel has a vanishing
    // FIRST MOMENT up to its own truncation, so a ramp still comes back correct to float precision:
    // measured worst error 7.3e-4 on values spanning 1000 (7e-7 relative, ~12 ulps) and 2.9e-3 on values
    // spanning 9000 (3e-7 relative). The bounds below are those measurements rounded up one decade, and
    // they are RELATIVE to the ramp's own magnitude — an absolute bound would silently loosen as the
    // fixture grew. A kernel with a broken first moment (an asymmetric tap set, an off-by-one centre)
    // fails them by orders, which is what these groups are for.
    group ("DeepSeek 1.2 adapted: 44.1k -> 88.2k ramp survives to float precision at half-sample phases");
    {
        const double D = felitronics::core::StreamResampler::delayInputSamples();
        const auto in = ramp (4096);
        const auto out = oneShot (44100.0, 88200.0, in);
        double worst = 0.0, span = 0.0;
        for (std::size_t k = 128; k + 128 < out.size(); ++k)
        {
            const double want = 0.5 * ((double) k * 0.5 - D) * 2.0;   // ramp() is x[i] = i; ipo = 0.5
            worst = std::max (worst, std::fabs ((double) out[k] - want));
            span  = std::max (span, std::fabs (want));
        }
        std::printf ("      ipo=0.5 ramp: worst |err| %.3e over a span of %.0f (relative %.2e)\n",
                     worst, span, worst / std::max (1.0, span));
        ok (out.size() > 6000 && worst / std::max (1.0, span) < 1.0e-5,
            "linear ramp maps to y[k] = k·ipo - kHalf to float precision — the first moment vanishes");
    }

    group ("DeepSeek 1.3 adapted: 48k -> 16k integer-step ramp, where every phase is exactly zero");
    {
        // ipo = 3 exactly, so the phase never leaves row 0 of the table. This is the case that would
        // still be bit-exact for an interpolating kernel and is NOT for this one — a useful reminder
        // that "integer ratio" no longer means "sample picking". It also means the whole result rides on
        // ONE table row, so a single mis-normalised row shows up here and nowhere else.
        const double D = felitronics::core::StreamResampler::delayInputSamples();
        const auto in = ramp (4096);
        const auto out = oneShot (48000.0, 16000.0, in);
        double worst = 0.0, span = 0.0;
        for (std::size_t k = 64; k + 64 < out.size(); ++k)
        {
            const double want = (double) k * 3.0 - D;
            worst = std::max (worst, std::fabs ((double) out[k] - want));
            span  = std::max (span, std::fabs (want));
        }
        std::printf ("      ipo=3 ramp: worst |err| %.3e over a span of %.0f (relative %.2e)\n",
                     worst, span, worst / std::max (1.0, span));
        ok (out.size() > 1000 && worst / std::max (1.0, span) < 1.0e-5,
            "…and the same at integer phases, where a single table row carries the whole result");
    }

    group ("DeepSeek 1.4 adapted: one-second 44.1k -> 48k ramp has bounded float precision error, no drift");
    {
        const auto in = ramp (44100);
        const auto out = oneShot (44100.0, 48000.0, in);
        const double step = 44100.0 / 48000.0;
        const double D = felitronics::core::StreamResampler::delayInputSamples();
        double worst = 0.0;
        double firstWindow = 0.0, lastWindow = 0.0;
        int firstN = 0, lastN = 0;
        for (std::size_t k = 128; k + 128 < out.size(); ++k)
        {
            const double want = step * (double) k - D;
            const double err = std::fabs ((double) out[k] - want);
            worst = std::max (worst, err / std::max (1.0, std::fabs (want)));    // RELATIVE, see below
            if (k < 1024) { firstWindow += err / std::max (1.0, std::fabs (want)); ++firstN; }
            if (k + 1024 > out.size()) { lastWindow += err / std::max (1.0, std::fabs (want)); ++lastN; }
        }
        const double head = firstWindow / std::max (1, firstN), tail = lastWindow / std::max (1, lastN);
        std::printf ("      worst RELATIVE ramp error %.3e, head avg %.3e, tail avg %.3e\n", worst, head, tail);
        // 🔴 THE MEASURE IS RELATIVE, and the reason is the whole point of the group. A ramp's VALUE
        // grows linearly, so its float quantization error grows with it: an absolute bound would read
        // "the tail is 20x worse than the head" on a perfectly drift-free resampler, which is what an
        // earlier version of this line did. Accumulator drift is what would show up as a growing
        // RELATIVE error; float quantization is flat in relative terms. Absolute worst here is 0.023 on
        // values near 40 000 — six ulps of a 64-term float dot product, where the cubic's four-term one
        // cost about 1.5 (sqrt(64/4) = 4x, and that is the whole difference).
        ok (worst < 5.0e-6, "float-output ramp error stays at the float ulp scale over 1 second");
        ok (tail < 4.0 * std::max (head, 1.0e-8),
            "…and the tail is no worse RELATIVE to the signal than the head: float quantization, not "
            "accumulator drift, which would grow without bound");
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
        const double D = felitronics::core::StreamResampler::delayInputSamples();
        double worst = 0.0, span = 0.0;
        for (std::size_t k = 128; k + 128 < out.size(); ++k)
        {
            const double x = 0.5 * (double) k - D;
            worst = std::max (worst, std::fabs ((double) out[k] - x * x * scale));
            span  = std::max (span, x * x * scale);
        }
        // A QUADRATIC is where the two kernels genuinely part company: the cubic reproduced it exactly
        // (its second moment vanishes by construction), a truncated sinc does not. The residual is the
        // kernel's own second moment against a signal that is not band-limited at all, so this bound is
        // a MEASUREMENT with a decade of headroom, not a derivation — and it is stated as such.
        std::printf ("      quadratic: worst |err| %.3e over a span of %.0f (relative %.2e)\n",
                     worst, span, worst / std::max (1.0, span));
        ok (worst / std::max (1.0, span) < 1.0e-3,
            "quadratic reproduction is no longer exact — an approximating kernel does not pass through "
            "its samples — but it holds to 1e-3 relative across half-sample phases");
    }

    return felitronics::test::report();
}

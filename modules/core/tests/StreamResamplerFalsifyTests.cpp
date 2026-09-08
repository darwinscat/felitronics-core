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
        const double D = felitronics::core::StreamResampler::delayInputSamples (44100.0, 88200.0);
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
        const double D = felitronics::core::StreamResampler::delayInputSamples (48000.0, 16000.0);
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
        const double D = felitronics::core::StreamResampler::delayInputSamples (44100.0, 48000.0);
        double worst = 0.0, worstAbs = 0.0, maxWant = 0.0;
        double firstWindow = 0.0, lastWindow = 0.0;
        int firstN = 0, lastN = 0;
        for (std::size_t k = 128; k + 128 < out.size(); ++k)
        {
            const double want = step * (double) k - D;
            const double err = std::fabs ((double) out[k] - want);
            worstAbs = std::max (worstAbs, err);
            maxWant  = std::max (maxWant, std::fabs (want));
            worst = std::max (worst, err / std::max (1.0, std::fabs (want)));    // RELATIVE, see below
            if (k < 1024) { firstWindow += err / std::max (1.0, std::fabs (want)); ++firstN; }
            if (k + 1024 > out.size()) { lastWindow += err / std::max (1.0, std::fabs (want)); ++lastN; }
        }
        const double head = firstWindow / std::max (1, firstN), tail = lastWindow / std::max (1, lastN);
        std::printf ("      worst RELATIVE ramp error %.3e, head avg %.3e, tail avg %.3e (worst ABS %.3e)\n",
                     worst, head, tail, worstAbs);
        // 🔴 TWO CORRECTIONS A CREW ROUND FORCED HERE, and the second is the interesting one.
        //
        // (1) The error is NOT float quantization. An earlier version of this comment said "six ulps of
        //     a 64-term float dot product", which sounded plausible and was wrong. The residual is the
        //     kernel's own FIRST MOMENT: the 64 taps cover [-31-t, 32-t], which is symmetric about zero
        //     only at t = 0 and t = 0.5, so for every other phase the edge taps are unequal
        //     (kernelAt(31.7) = -1.72e-5 against kernelAt(-31.3) = +1.09e-6) and a ramp comes out
        //     shifted by up to 3.6e-4 of a sample. Times the ramp's slope, that is a fixed ABSOLUTE
        //     error, identical on clang, gcc and MSVC — deterministic, not noise.
        //
        // (2) Because it is absolute and fixed, a RELATIVE bound is silently fitted to WHERE THE LOOP
        //     STARTS: the ramp's value grows, so the same error divided by a bigger number looks
        //     smaller. Starting at k = 64 instead of 128 makes the identical kernel read 1.3e-5 and
        //     fail a 5e-6 line. So the bound on the SIZE is absolute now, derived from the first moment
        //     (3.6e-4 samples x the 0.919 slope = 3.3e-4; asserted at 1e-3), and the RELATIVE pair is
        //     kept for the claim it can actually make: that the error does not GROW along the ramp,
        //     which is what accumulator drift would do and float quantization would not.
        // The bound has TWO terms because the error does, and single-term versions of this line failed
        // the correct kernel twice while I was writing it — once too tight by 80x, once by 1.2x.
        //   Term 1, the kernel's FIRST MOMENT: a fixed 3.3e-4 absolute (see above), dominant while the
        //     ramp is still small. Rounded up to 1e-3.
        //   Term 2, float SUMMATION of large values: the dot product adds 64 products of magnitude up
        //     to |x|·max|w|, and the coefficient vector's ABSOLUTE sum is 2.77 (measured on the shipped
        //     table — a windowed sinc's negative lobes carry far more than its unit sum suggests). The
        //     realistic accumulation error is sqrt(kTaps)·u·Σ|w|·|x| with u = 2^-24; the adversarial
        //     one is kTaps·u·Σ|w|·|x|, which is 8x larger and would make this line meaningless.
        // Measured 0.0266 against the sqrt form's 0.058 — a factor of 2, on a fixture whose values run
        // to 44 000. NEITHER term grows with k, which is the claim the group exists to make.
        constexpr double kAbsCoeffSum = 2.78;                  // measured max Σ|w| over all phase rows
        const double bound = 1.0e-3 + std::sqrt ((double) felitronics::core::StreamResampler::kTaps)
                                      * std::ldexp (1.0, -24) * kAbsCoeffSum * maxWant;
        std::printf ("      derived bound %.3e (first moment 1e-3 + sqrt(64)·2^-24·2.78·%.0f)\n", bound, maxWant);
        ok (worstAbs < bound,
            "the ramp error is the kernel's first moment plus float storage of a large value ("
            + std::to_string (worstAbs) + " against a derived " + std::to_string (bound)
            + "), not accumulator drift");
        ok (tail < 4.0 * std::max (head, 1.0e-8),
            "…and RELATIVE to the signal the tail is no worse than the head, which is the drift claim: "
            "a drifting phase accumulator grows without bound, a fixed first moment does not");
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
        const double D = felitronics::core::StreamResampler::delayInputSamples (48000.0, 96000.0);
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
        // 🔴 1e-3 RELATIVE WAS TOO LOOSE AND A CREW ROUND SAID SO WITH A NUMBER: on a fixture whose
        // values reach 15621, it admitted 15.6 of absolute error, and a harmful mutation walked through.
        // Measured is 3.9e-7 relative; 1e-5 keeps a factor of 25 and closes the door.
        ok (worst / std::max (1.0, span) < 1.0e-5,
            "quadratic reproduction is no longer exact — an approximating kernel does not pass through "
            "its samples — but it holds to 1e-5 relative across half-sample phases");
    }

    return felitronics::test::report();
}

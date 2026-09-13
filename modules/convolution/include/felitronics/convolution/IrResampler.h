// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/core/Math.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

namespace felitronics::convolution
{

namespace detail
{
    // Modified Bessel I0 (series) — for the Kaiser window. Converges fast for moderate beta.
    inline double besselI0 (double x) noexcept
    {
        double sum = 1.0, term = 1.0;
        const double y = x * x * 0.25;
        for (int k = 1; k < 64; ++k)
        {
            term *= y / ((double) k * (double) k);
            sum  += term;
            if (term < 1e-17 * sum) break;
        }
        return sum;
    }
}

struct IrResampleConfig
{
    int    halfTaps    = 32;     // window radius in input samples (64-tap filter) — more = sharper transition
    double beta        = 8.0;    // Kaiser beta (~80 dB stopband); 5.65 ≈ 60 dB
    double cutoffScale = 0.95;   // fraction of the (lower) Nyquist used as the passband edge
};

//==============================================================================
// Offline windowed-sinc (Kaiser) IR resampler. MESSAGE-THREAD ONLY (double math, allocates) — for
// rate-converting an impulse response to the host SR on load. DC gain is normalized to 1.
//
// OUTSIDE THE INPUT IS SILENCE. Every output sample divides by the weight of its WHOLE window, and a tap
// that falls past either end of the input adds its weight but no signal: the IR is resampled as if it were
// padded with zeros, which is what it is — a cabinet trimmed at its onset had silence before it, and the
// convolver plays silence after its last tap. Those taps used to be skipped BEFORE their weight was added,
// so an edge sample was divided by only the part of its window that landed on the input, and the samples
// that were not there counted as the weighted mean of the ones that were. A cabinet's onset is its loudest
// edge; that invented a broadband floor over the top octave, which is exactly where a cabinet is quietest.
// What zeros cannot give back is the kernel's pre-ringing that would fall BEFORE output sample 0: an IR that
// starts at sample 0 keeps it only by adding delay, and it costs more the more abruptly the IR starts.
// `felitronics_convolution_resampler_tests` owns the numbers — shift invariance, and a cabinet-like IR's
// band response against its own response and against the untruncated resample.
//
// A LOAD IS NEVER RESAMPLED TO NOTHING. The length is inLen*ratio rounded but at least one sample, as JUCE's
// resampleImpulseResponse had it: a one-tap IR at 96 -> 44.1 kHz rounded to zero taps, and the loader had
// nothing to publish. The result is empty only for no input (a null pointer or a non-positive length), a rate
// that is not a positive finite number, or an output `int` cannot address — longer than INT_MAX samples, or a
// ratio so small that output sample 0 alone sits past INT_MAX (the floor of one sample is what makes that
// reachable: `(int) floor(t)` would be undefined).
//
// FAMILY SPLIT vs core::StreamResampler — restated, because the other half of it changed under this
// comment. That one used to be a Catmull-Rom cubic and "too low-SNR for IRs" was the whole argument.
// Since P34 it is a 64-tap polyphase windowed sinc, i.e. the SAME family as this one, so the split is
// no longer about quality. It is about BUDGET and THREAD: this is an offline one-shot that may
// allocate, work in double and size its kernel to the job; that is a streaming rate-match on the audio
// thread with a fixed table built in reset() and a per-block cost that has to stay inside a couple of
// percent of a neural stage. Still not interchangeable — for the opposite reason to the one that used
// to be written here.
inline std::vector<float> resampleIr (const float* in, int inLen, double inSr, double outSr,
                                      IrResampleConfig cfg = {})
{
    std::vector<float> out;
    if (in == nullptr || inLen <= 0 || ! (inSr > 0.0) || ! (outSr > 0.0) || ! std::isfinite (inSr)
        || ! std::isfinite (outSr)) return out;

    const double ratio = outSr / inSr;
    if (! (ratio > 0.0) || ! std::isfinite (ratio)) return out;            // finite rates can still under/overflow here
    const double want = (double) inLen * ratio;
    if (! (want < (double) std::numeric_limits<int>::max())) return out;   // before the cast: a wrapped length is garbage
    const int outLen = std::max (1, (int) std::llround (want));            // never zero taps (see above)

    // Sanitize the config — halfTaps < 1 makes the tap loop empty (an all-zero "IR"), a non-finite
    // beta/cutoffScale poisons every tap.
    const int    R      = cfg.halfTaps < 1 ? 1 : cfg.halfTaps;
    const double beta   = (std::isfinite (cfg.beta) && cfg.beta >= 0.0) ? cfg.beta : 8.0;
    const double cScale = (std::isfinite (cfg.cutoffScale) && cfg.cutoffScale > 0.0 && cfg.cutoffScale <= 1.0)
                        ? cfg.cutoffScale : 0.95;
    const double fc     = 0.5 * std::min (1.0, ratio) * cScale;            // cycles per INPUT sample
    const double i0beta = detail::besselI0 (beta);

    // Every tap index is an int: the last output's input position plus the window radius has to fit.
    const double tLast = ((double) outLen - 0.5) / ratio - 0.5;
    if (! (tLast + (double) R + 2.0 < (double) std::numeric_limits<int>::max())) return out;
    out.assign ((std::size_t) outLen, 0.0f);

    for (int n = 0; n < outLen; ++n)
    {
        const double t = ((double) n + 0.5) / ratio - 0.5;                 // output n → input position (sample-centred)
        const int    c = (int) std::floor (t);
        double acc = 0.0, wsum = 0.0;
        for (int k = c - R + 1; k <= c + R; ++k)
        {
            const double xx   = t - (double) k;
            const double sinc = (std::fabs (xx) < 1e-12) ? (2.0 * fc)
                                                         : std::sin (2.0 * core::kPi * fc * xx) / (core::kPi * xx);
            const double r    = xx / (double) R;                           // window argument in [-1, 1]
            const double win  = (r <= -1.0 || r >= 1.0) ? 0.0
                              : detail::besselI0 (beta * std::sqrt (1.0 - r * r)) / i0beta;
            const double w    = sinc * win;
            wsum += w;                                                     // the WHOLE window: OUTSIDE THE INPUT IS SILENCE
            if (k >= 0 && k < inLen) acc += (double) in[k] * w;
        }
        out[(std::size_t) n] = (float) (! core::exactlyEqual (wsum, 0.0) ? acc / wsum : 0.0);    // normalize → unity DC (intentional exact ==)
    }
    return out;
}

inline std::vector<float> resampleIr (const std::vector<float>& in, double inSr, double outSr,
                                      IrResampleConfig cfg = {})
{
    return resampleIr (in.data(), (int) in.size(), inSr, outSr, cfg);
}

} // namespace felitronics::convolution

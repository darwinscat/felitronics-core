// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/core/Math.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
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
// rate-converting an impulse response to the host SR on load. This is the >=60 dB-class resampler the
// convolution path needs; the Catmull-Rom cab::StreamResampler (for NAM rate-match) is too low-SNR for
// IRs (see dsp-shared-dsp-review §4). DC gain is normalized to 1.
inline std::vector<float> resampleIr (const float* in, int inLen, double inSr, double outSr,
                                      IrResampleConfig cfg = {})
{
    std::vector<float> out;
    if (in == nullptr || inLen <= 0 || inSr <= 0.0 || outSr <= 0.0) return out;

    const double ratio  = outSr / inSr;
    const int    outLen = (int) std::llround ((double) inLen * ratio);
    if (outLen <= 0) return out;
    out.assign ((std::size_t) outLen, 0.0f);

    const double fc     = 0.5 * std::min (1.0, ratio) * cfg.cutoffScale;   // cycles per INPUT sample
    const int    R      = cfg.halfTaps;
    const double i0beta = detail::besselI0 (cfg.beta);

    for (int n = 0; n < outLen; ++n)
    {
        const double t = ((double) n + 0.5) / ratio - 0.5;                 // output n → input position (sample-centred)
        const int    c = (int) std::floor (t);
        double acc = 0.0, wsum = 0.0;
        for (int k = c - R + 1; k <= c + R; ++k)
        {
            if (k < 0 || k >= inLen) continue;
            const double xx   = t - (double) k;
            const double sinc = (std::fabs (xx) < 1e-12) ? (2.0 * fc)
                                                         : std::sin (2.0 * core::kPi * fc * xx) / (core::kPi * xx);
            const double r    = xx / (double) R;                           // window argument in [-1, 1]
            const double win  = (r <= -1.0 || r >= 1.0) ? 0.0
                              : detail::besselI0 (cfg.beta * std::sqrt (1.0 - r * r)) / i0beta;
            const double w    = sinc * win;
            acc  += (double) in[k] * w;
            wsum += w;
        }
        out[(std::size_t) n] = (float) (wsum != 0.0 ? acc / wsum : 0.0);    // normalize → unity DC
    }
    return out;
}

inline std::vector<float> resampleIr (const std::vector<float>& in, double inSr, double outSr,
                                      IrResampleConfig cfg = {})
{
    return resampleIr (in.data(), (int) in.size(), inSr, outSr, cfg);
}

} // namespace felitronics::convolution

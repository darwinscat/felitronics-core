// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/core/DetMath.h>
#include <felitronics/core/Math.h>

#include <cmath>

//==============================================================================
// felitronics::codecgrid::window — the analysis windows of the transform codecs, designed in double through
// core::det so that they are the same numbers on every platform: AAC's sine and Kaiser-Bessel-derived windows,
// CELT's low-overlap window. Each satisfies Princen-Bradley (w[n]^2 + w[n + M]^2 = 1), which is what makes an
// MDCT with it invertible, and a codec's decoder uses the same one — so does an analysis that wants to land on
// the encoder's coefficients.
//==============================================================================
namespace felitronics::codecgrid
{

//==============================================================================
// The windows. Each fills w[0 .. 2M - 1] for a frame of 2M samples.
//==============================================================================
namespace window
{
    // sin (pi / (2M) (n + 1/2)) — AAC's sine window, and MP3's at 2M = 36.
    inline void sine (double* w, int frame) noexcept
    {
        for (int n = 0; n < frame; ++n) w[n] = core::det::sin (core::kPi / (double) frame * ((double) n + 0.5));
    }

    namespace detail
    {
        // The modified Bessel function I0 by its series; 60 terms settle it far past double for the arguments used.
        inline double bessel0 (double x) noexcept
        {
            double sum = 1.0, term = 1.0;
            for (int k = 1; k < 60; ++k)
            {
                const double h = x / (2.0 * (double) k);
                term *= h * h;
                sum += term;
            }
            return sum;
        }
    }

    // The Kaiser-Bessel-derived window: the square root of the running sum of a Kaiser window of M + 1 points
    // with beta = pi alpha, mirrored. AAC's long window uses alpha = 4. `scratch` takes M + 1 doubles.
    inline void kaiserBesselDerived (double* w, int frame, double alpha, double* scratch) noexcept
    {
        const int half = frame / 2;
        const double beta = core::kPi * alpha, denom = detail::bessel0 (beta);
        double acc = 0.0;
        for (int i = 0; i <= half; ++i)
        {
            const double r = 2.0 * (double) i / (double) half - 1.0;
            acc += detail::bessel0 (beta * std::sqrt (1.0 - r * r)) / denom;
            scratch[i] = acc;
        }
        for (int i = 0; i < half; ++i)
        {
            w[i] = std::sqrt (scratch[i] / acc);
            w[frame - 1 - i] = w[i];
        }
    }

    // CELT's long-block window (RFC 6716, 4.3.7): of the 2M samples only the middle carries anything — zeros,
    // a rise of `overlap` samples, ones, the fall, zeros. The rise is sin (pi/2 sin^2 (pi/2 (i + 1/2) / overlap)),
    // which is power-complementary with its own mirror.
    inline void celt (double* w, int frame, int overlap) noexcept
    {
        const int hop = frame / 2, edge = (hop - overlap) / 2;
        for (int n = 0; n < frame; ++n) w[n] = 0.0;
        for (int i = 0; i < overlap; ++i)
        {
            const double s = core::det::sin (0.5 * core::kPi * ((double) i + 0.5) / (double) overlap);
            const double v = core::det::sin (0.5 * core::kPi * s * s);
            w[edge + i] = v;
            w[frame - 1 - edge - i] = v;
        }
        for (int n = edge + overlap; n < frame - edge - overlap; ++n) w[n] = 1.0;
    }
}

} // namespace felitronics::codecgrid

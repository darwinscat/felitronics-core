// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/core/Config.h>
#include <felitronics/core/Math.h>

#include <cmath>

namespace felitronics::analysis
{

//==============================================================================
// felitronics::analysis::KWeightingFilter — the ITU-R BS.1770 K-weighting pre-filter: stage 1 (a high
// shelf, ~+4 dB) then stage 2 (the RLB high-pass, ~38 Hz). Coefficients are RECOMPUTED for the actual
// sample rate (analog prototype + bilinear) — correct at 44.1/48/88.2/96 kHz, not a hardcoded 48 k table
// (the usual porting bug). Per-channel double biquad state for loudness accuracy.
class KWeightingFilter
{
public:
    struct Coeffs { double b0, b1, b2, a1, a2; };

    void prepare (double sampleRate, int numChannels) noexcept
    {
        fs = sampleRate;
        ch = numChannels < 1 ? 1 : (numChannels > kMaxChannels ? kMaxChannels : numChannels);
        computeCoeffs (fs);
        reset();
    }

    void reset() noexcept
    {
        for (int c = 0; c < kMaxChannels; ++c) { z1a[c] = z2a[c] = 0.0; z1b[c] = z2b[c] = 0.0; }
    }

    // K-weighted sample (double) for channel c. Two TDF-II biquads in series.
    inline double process (int c, double x) noexcept
    {
        const double y1 = sa.b0 * x  + z1a[c]; z1a[c] = sa.b1 * x  - sa.a1 * y1 + z2a[c]; z2a[c] = sa.b2 * x  - sa.a2 * y1;
        const double y2 = sb.b0 * y1 + z1b[c]; z1b[c] = sb.b1 * y1 - sb.a1 * y2 + z2b[c]; z2b[c] = sb.b2 * y1 - sb.a2 * y2;
        return y2;
    }

    const Coeffs& shelfCoeffs() const noexcept { return sa; }   // for tests (vs the canonical 48 k values)
    const Coeffs& highpassCoeffs() const noexcept { return sb; }

private:
    void computeCoeffs (double sampleRate) noexcept
    {
        // Stage 1 — high shelf (the published BS.1770 design constants).
        {
            const double f0 = 1681.974450955533, G = 3.999843853973347, Q = 0.7071752369554196;
            const double K = std::tan (core::kPi * f0 / sampleRate);
            const double Vh = std::pow (10.0, G / 20.0);
            const double Vb = std::pow (Vh, 0.4996667741545416);
            const double D = 1.0 + K / Q + K * K;
            sa = { (Vh + Vb * K / Q + K * K) / D, 2.0 * (K * K - Vh) / D, (Vh - Vb * K / Q + K * K) / D,
                   2.0 * (K * K - 1.0) / D, (1.0 - K / Q + K * K) / D };
        }
        // Stage 2 — RLB high-pass.
        {
            const double f0 = 38.13547087602444, Q = 0.5003270373238773;
            const double K = std::tan (core::kPi * f0 / sampleRate);
            const double D = 1.0 + K / Q + K * K;
            sb = { 1.0, -2.0, 1.0, 2.0 * (K * K - 1.0) / D, (1.0 - K / Q + K * K) / D };
        }
    }

    static constexpr int kMaxChannels = core::kMaxChannels;
    double fs = 48000.0; int ch = 2;
    Coeffs sa {}, sb {};
    double z1a[kMaxChannels] {}, z2a[kMaxChannels] {}, z1b[kMaxChannels] {}, z2b[kMaxChannels] {};
};

} // namespace felitronics::analysis

// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/core/FlushToZero.h>

#include <algorithm>
#include <cmath>

namespace felitronics::analysis
{

//==============================================================================
// felitronics::analysis::CorrelationMeter — the inter-channel correlation coefficient ρ ∈ [-1, 1]
// (one-pole windowed cross / auto). +1 = identical (mono), -1 = inverted (anti-phase → mono cancels),
// ~0 = uncorrelated/wide. RT-safe streaming. ρ = E[L·R] / sqrt(E[L²]·E[R²]).
class CorrelationMeter
{
public:
    void prepare (double sampleRate, double windowMs = 300.0) noexcept { fs = sampleRate; setWindow (windowMs); reset(); }
    void reset() noexcept { sLL = sRR = sLR = 0.0; }

    void setWindow (double windowMs) noexcept
    {
        const double t = windowMs * 0.001;
        alpha = (t <= 0.0 || fs <= 0.0) ? 1.0 : (1.0 - std::exp (-1.0 / (t * fs)));
    }

    inline void process (float l, float r) noexcept
    {
        sLL += alpha * ((double) l * l - sLL);
        sRR += alpha * ((double) r * r - sRR);
        sLR += alpha * ((double) l * r - sLR);
    }

    double correlation() const noexcept
    {
        const double d = std::sqrt (sLL * sRR);
        return d > 1e-12 ? std::clamp (sLR / d, -1.0, 1.0) : 1.0;
    }

    void flushDenormals() noexcept { core::flushDenormal (sLL); core::flushDenormal (sRR); core::flushDenormal (sLR); }

private:
    double fs = 48000.0, alpha = 0.0, sLL = 0.0, sRR = 0.0, sLR = 0.0;
};

} // namespace felitronics::analysis

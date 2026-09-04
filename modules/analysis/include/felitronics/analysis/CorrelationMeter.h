// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

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
        alpha = (! (t > 0.0) || fs <= 0.0) ? 1.0 : (1.0 - std::exp (-1.0 / (t * fs)));   // !(t>0) also catches NaN
    }

    inline void process (float l, float r) noexcept
    {
        // A single NaN/inf sample would latch all three one-pole accumulators at NaN FOREVER (the state
        // never heals) and correlation() would read a false +1.0 from then on. Treat it as 0.
        if (! std::isfinite (l)) l = 0.0f;
        if (! std::isfinite (r)) r = 0.0f;
        sLL += alpha * ((double) l * l - sLL);
        sRR += alpha * ((double) r * r - sRR);
        sLR += alpha * ((double) l * r - sLR);
        flushDenormals();   // law 8, HERE and not left to the caller — see the note on flushDenormals()
    }

    double correlation() const noexcept
    {
        const double d = std::sqrt (sLL * sRR);
        return d > 1e-12 ? std::clamp (sLR / d, -1.0, 1.0) : 1.0;
    }

    // Law 8. `process()` already calls this every sample, so a consumer never has to — the method stays
    // public only because it was, and because a host may want to zap the state at a transport jump.
    // NB the stall here is SLOW, not absent: at the default 300 ms window silence needs ~220 s to walk the
    // state down into the subnormal band. That is a long tail on a mix bus, not an exemption.
    //
    // Why per sample rather than per block, and why not left to the owner: this meter has no owner. It is
    // a leaf primitive with no in-repo caller — the adapter that drives it is the only thing that could
    // have called a flush, and for the whole life of the class nothing did, which is the F1
    // `MultibandSplitter` failure exactly (the method existed; no line invoked it). On silence all three
    // one-poles decay geometrically and STICK: for the default 300 ms window at 48 kHz every subnormal
    // k·u with k ≤ ½/alpha ≈ 7200 maps to itself, so the state never reaches zero and every later sample
    // pays the subnormal penalty on three dependent double FMAs. The cost of preventing that is three
    // fabs+compares against a three-FMA dependency chain the loop is latency-bound on anyway.
    void flushDenormals() noexcept { core::flushDenormal (sLL); core::flushDenormal (sRR); core::flushDenormal (sLR); }

private:
    double fs = 48000.0, alpha = 0.0, sLL = 0.0, sRR = 0.0, sLR = 0.0;
};

} // namespace felitronics::analysis

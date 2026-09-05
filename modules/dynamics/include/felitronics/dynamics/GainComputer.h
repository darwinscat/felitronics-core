// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <cmath>

namespace felitronics::dynamics
{

// down-compress  = tame when loud (the common dynamic-EQ move).
// up-compress    = lift when quiet.
// down-expand    = duck when quiet (gate-ish).
enum class Mode { DownCompress, UpCompress, DownExpand };

//==============================================================================
// felitronics::dynamics::GainComputer — the dB-domain static curve: a level in dB → a signed gain
// DELTA in dB, clamped to [-range, +range]. Pure and STATELESS (no signal memory): the product feeds
// it 20·log10(envelope) and applies the returned delta (e.g. through an SVF gain-delta, per the
// matched-static × SVF-delta composition — which is product glue, not part of this primitive).
//
// The soft knee is SYMMETRIC (centred on the threshold, smearing ±knee/2) — the textbook shape: just as
// down-compress eases in below threshold, up-compress / down-expand intentionally apply a small tapering
// amount on the inactive side within the knee. knee == 0 → a hard corner.
class GainComputer
{
public:
    // Finite-guarded: a NaN/inf param can't leak a NaN delta into the gain path (the consumer tests
    // assert finiteness). ±inf input levels are handled by the range clamp in deltaDb().
    void setThresholdDb (double dB) noexcept { thr   = std::isfinite (dB) ? dB : 0.0; }
    void setRatio  (double r) noexcept { ratio = (std::isfinite (r) && r > 1.0) ? r : 1.0; }
    void setKneeDb (double k) noexcept { knee  = (std::isfinite (k) && k > 0.0) ? k : 0.0; }
    // BOUNDED, not merely finite. The delta this caps is cast to float by every consumer, and a finite
    // but enormous range let a finite threshold produce a delta of ~1e300 that became -Inf in float; the
    // GR follower then evaluated `-Inf + c*(0 - -Inf)` = NaN and never recovered, so ordinary parameters
    // restored afterwards produced NaN audio for the rest of the stream. 400 dB is a gain of 1e20 —
    // six times any real range and far inside float — so the cap is unreachable in use.
    void setRangeDb (double r) noexcept
    { range = (std::isfinite (r) && r > 0.0) ? (r < kMaxRangeDb ? r : kMaxRangeDb) : 0.0; }   // ± delta cap
    void setMode (Mode m) noexcept { mode = m; }

    double thresholdDb() const noexcept { return thr; }
    Mode   modeValue()   const noexcept { return mode; }

    // level (dB) → signed gain delta (dB), clamped to [-range, +range]. delta == 0 ⇒ transparent.
    double deltaDb (double levelDb) const noexcept
    {
        const double slope = 1.0 - 1.0 / ratio;      // compressor slope (0 at 1:1, →1 at ∞:1)
        double delta = 0.0;
        switch (mode)
        {
            case Mode::DownCompress: delta = -kneeOver (levelDb - thr) * slope;          break;  // cut as it gets loud
            case Mode::UpCompress:   delta = +kneeOver (thr - levelDb) * slope;          break;  // lift as it gets quiet
            case Mode::DownExpand:   delta = -kneeOver (thr - levelDb) * (ratio - 1.0);  break;  // duck below threshold
        }
        if (delta >  range) delta =  range;
        if (delta < -range) delta = -range;
        // The clamp above is TWO ORDERED COMPARISONS, and both are false for a NaN — so it does not
        // make the header's finiteness promise true on its own. One reachable NaN exists with entirely
        // valid params: ratio == 1.0 (the honest "1:1, off" setting, and also the fallback setRatio()
        // applies to any invalid ratio) makes slope exactly 0, and an infinite level then computes
        // `kneeOver(+Inf) * 0` = `Inf * 0` = NaN, which reaches the gain as dbToGain(NaN) and turns the
        // whole stream into NaN permanently. Substituting 0 (= transparent) is the only answer that
        // carries no invented information. Bit-transparent for every finite and infinite input: those
        // are already resolved by the clamps above.
        if (std::isnan (delta)) delta = 0.0;
        return delta;
    }

private:
    // Soft-knee "amount past threshold": 0 below the knee, a C1-continuous quadratic through it, and
    // linear (== x) above. knee == 0 → hard knee (max(0, x)). Symmetric half-knee = `knee/2` either
    // side of the corner. Value+slope are continuous at both knee edges.
    double kneeOver (double x) const noexcept
    {
        if (knee <= 0.0) return x > 0.0 ? x : 0.0;
        const double h = knee * 0.5;
        if (x <= -h) return 0.0;
        if (x >=  h) return x;
        const double t = x + h;                       // 0 .. knee
        return (t * t) / (2.0 * knee);
    }

    static constexpr double kMaxRangeDb = 400.0;

    double thr = -18.0, ratio = 2.0, knee = 0.0, range = 24.0;
    Mode   mode = Mode::DownCompress;
};

} // namespace felitronics::dynamics

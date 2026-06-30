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
    void setRangeDb (double r) noexcept { range = (std::isfinite (r) && r > 0.0) ? r : 0.0; }   // ± delta cap
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

    double thr = -18.0, ratio = 2.0, knee = 0.0, range = 24.0;
    Mode   mode = Mode::DownCompress;
};

} // namespace felitronics::dynamics

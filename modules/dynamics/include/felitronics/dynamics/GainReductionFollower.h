// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/core/FlushToZero.h>
#include <felitronics/core/Math.h>

#include <cmath>

namespace felitronics::dynamics
{

//==============================================================================
// felitronics::dynamics::GainReductionFollower — attack/release ballistics applied to the GAIN
// REDUCTION (a signed gain delta in dB), NOT to the detector level. "Attack" = the magnitude growing
// (more processing — fast); "release" = shrinking back toward 0 dB (slow). This is the clean-compressor
// topology: the static curve (threshold/ratio/knee, in the stateless GainComputer) and the timing
// (here) are decoupled, so the knee can't warp the attack/release the way detector-level smoothing does.
class GainReductionFollower
{
public:
    void prepare (double sampleRate) noexcept { fs = sampleRate; updateCoeffs(); reset(); }
    void reset() noexcept { currentDb = 0.0f; }

    void setTimes (double attackMs, double releaseMs) noexcept { atkMs = attackMs; relMs = releaseMs; updateCoeffs(); }

    // One target gain delta (dB) in → smoothed gain delta (dB) out. RT-safe.
    inline float process (float targetDb) noexcept
    {
        const float c = (std::fabs (targetDb) > std::fabs (currentDb)) ? atkCoeff : relCoeff;
        currentDb = targetDb + c * (currentDb - targetDb);
        return currentDb;
    }

    float valueDb() const noexcept { return currentDb; }

    // LAW 11c — `n` steps toward a target that is KNOWN CONSTANT for all of them. This is the cheap half
    // of "a pause is silence": once the detector level has reached `core::kGainToDbFloor` the dB
    // conversion returns the same bits for every smaller level, so the static curve behind it returns the
    // same delta, and the whole per-sample path collapses to the one multiply-add below.
    //
    // IT CALLS `process()` RATHER THAN INLINING ITS ARITHMETIC, and that is the point: the coefficient
    // choice is `|target| > |current|`, which a constant target does NOT freeze. It can flip more than
    // once — target +1, current -2, release 0.5, attack 0 goes -2 -> -0.5 -> +1, i.e. release, attack,
    // release — so any "pick the coefficient once" shortcut is wrong, and the follower's own branch is
    // the only thing that gets it right. Nothing is saved by inlining it anyway; what was saved is
    // upstream, where a log10 and a knee no longer run.
    void advanceConstant (float targetDb, int n) noexcept
    {
        for (int i = 0; i < n; ++i)
        {
            const float before = currentDb;
            (void) process (targetDb);
            if (core::sameBits (currentDb, before)) return;             // a fixed point stays fixed
            if (! std::isfinite (currentDb) && ! std::isfinite (before)) return;
        }
    }

    // Law 8: flush the follower state once it decays below the subnormal-risk threshold (release → 0),
    // and clear it if it is ever non-finite. The second half is unreachable through `Compressor` — the
    // detector is gated and `GainComputer::deltaDb` substitutes 0 for a NaN delta — so this is the
    // primitive keeping its own promise for whoever drives it directly, not a fix for a live defect.
    void flushDenormals() noexcept { core::flushPoison (currentDb); }

private:
    void  updateCoeffs() noexcept { atkCoeff = coeff (atkMs); relCoeff = coeff (relMs); }
    float coeff (double ms) const noexcept
    {
        const double t = ms * 0.001;
        // Same guards as EnvelopeFollower::coeff, and for the same reasons: !(t>0) catches NaN, the
        // isfinite pair catches +Inf (exp(-1/Inf) == 1.0 = a FROZEN gain, i.e. a stuck gain, which is
        // worse than the "instant" every other invalid time falls back to), and the two bounds below
        // keep the coefficient strictly inside (0, 1) so the follower is never dead in either direction.
        if (! (t > 0.0) || ! std::isfinite (t) || ! (fs > 0.0) || ! std::isfinite (fs)) return 0.0f;
        const float c = (float) std::exp (-1.0 / (t * fs));
        if (c < 1.0e-15f) return 0.0f;
        return c < 1.0f ? c : 0x1.fffffep-1f;
    }

    double fs = 48000.0, atkMs = 10.0, relMs = 100.0;
    float  atkCoeff = 0.0f, relCoeff = 0.0f, currentDb = 0.0f;
};

} // namespace felitronics::dynamics

// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/core/FlushToZero.h>

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

    // Law 8: flush the follower state once it decays below the subnormal-risk threshold (release → 0).
    void flushDenormals() noexcept { core::flushDenormal (currentDb); }

private:
    void  updateCoeffs() noexcept { atkCoeff = coeff (atkMs); relCoeff = coeff (relMs); }
    float coeff (double ms) const noexcept
    {
        const double t = ms * 0.001;
        if (! (t > 0.0) || fs <= 0.0) return 0.0f;   // !(t>0) also catches a NaN time → instant, never a NaN coeff
        const float c = (float) std::exp (-1.0 / (t * fs));
        return c < 1.0e-15f ? 0.0f : c;     // never leave a subnormal coeff (the EnvelopeFollower fix)
    }

    double fs = 48000.0, atkMs = 10.0, relMs = 100.0;
    float  atkCoeff = 0.0f, relCoeff = 0.0f, currentDb = 0.0f;
};

} // namespace felitronics::dynamics

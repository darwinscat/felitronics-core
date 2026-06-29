// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/core/FlushToZero.h>

#include <cmath>

namespace felitronics::dynamics
{

// Peak tracks |x|; RMS tracks a one-pole average of x^2 (returned as sqrt → amplitude). RMS is the
// default for musical EQ dynamics; Peak for fast limiting-style detection.
//
// NB on time constants: Peak's attack/release act on |x| directly (the envelope reaches 1-1/e of a
// step per TC). RMS's act on the POWER (x^2), so the returned amplitude reaches sqrt(1-1/e)≈0.79 per
// attack TC, and an asymmetric attack≠release biases a steady tone's tracked level — use attack≈release
// for a true RMS meter.
enum class Detector { Peak, Rms };

//==============================================================================
// felitronics::dynamics::EnvelopeFollower — a one-pole attack/release follower on a mono sidechain
// probe. Pure signal-in → envelope-out; no EQ/param/GUI knowledge (the dynamic-EQ composition lives in
// the product). This is exactly the "unguarded feedback kernel" the 3rd review flagged for Law 8 — so
// it `flushDenormals()` its state every block (works on every tier, no hardware FTZ needed).
//
// RT-safe: process() does no alloc/lock/IO/throw. One instance per band+lane in the EQ use case.
class EnvelopeFollower
{
public:
    void prepare (double sampleRate) noexcept { fs = sampleRate; updateCoeffs(); reset(); }
    void reset() noexcept { env = 0.0f; }

    void setTimes (double attackMs, double releaseMs) noexcept
    {
        atkMs = attackMs; relMs = releaseMs;
        updateCoeffs();
    }

    void     setDetector (Detector d) noexcept { det = d; }
    Detector detector() const noexcept { return det; }

    // One sidechain sample in → current envelope OUT (linear amplitude). Attack coeff while rising,
    // release coeff while falling (the standard branch). RT-safe.
    inline float process (float x) noexcept
    {
        const float in = (det == Detector::Rms) ? x * x : std::fabs (x);
        const float c  = (in > env) ? atkCoeff : relCoeff;
        env = in + c * (env - in);
        return (det == Detector::Rms) ? std::sqrt (env) : env;
    }

    // Current envelope as linear amplitude (sqrt of the mean-square in RMS mode).
    float envelope() const noexcept { return (det == Detector::Rms) ? std::sqrt (env) : env; }

    // Law 8: zap the follower state to exact zero once it decays below the subnormal-risk threshold,
    // so a long silence can't sustain subnormals (CPU spike). Call once per block.
    void flushDenormals() noexcept { core::flushDenormal (env); }

private:
    void  updateCoeffs() noexcept { atkCoeff = coeff (atkMs); relCoeff = coeff (relMs); }
    float coeff (double ms) const noexcept
    {
        const double t = ms * 0.001;
        if (t <= 0.0 || fs <= 0.0) return 0.0f;
        const float c = (float) std::exp (-1.0 / (t * fs));
        return c < 1.0e-15f ? 0.0f : c;     // absurdly short time → instant; never leave a subnormal coeff
    }

    double   fs = 48000.0, atkMs = 10.0, relMs = 100.0;
    float    atkCoeff = 0.0f, relCoeff = 0.0f;
    float    env = 0.0f;            // peak: |x| envelope · rms: mean-square accumulator
    Detector det = Detector::Rms;
};

} // namespace felitronics::dynamics

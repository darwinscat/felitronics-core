// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

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
// the product). This is exactly the "unguarded feedback kernel" that is a known risk for Law 8 — so
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

    // `env` is NOT the same quantity in the two modes — it holds |x| in Peak and x^2 in Rms — so
    // flipping the enum alone REINTERPRETS the state: a Peak envelope of 0.01 read as a mean-square
    // reports sqrt(0.01) = 0.1, a silent +20 dB step (measured). Converting keeps the reported
    // AMPLITUDE continuous across the switch, which is what a mode change should sound like: nothing.
    // A no-op when the mode does not change, so the common setParams()-every-block path is untouched.
    void setDetector (Detector d) noexcept
    {
        if (d == det) return;
        if (d == Detector::Rms) { const float sq = env * env; env = std::isfinite (sq) ? sq : 0.0f; }
        else                    { env = env > 0.0f ? std::sqrt (env) : 0.0f; }
        det = d;
    }

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
    //
    // THE THRESHOLD IS SQUARED IN RMS MODE, and that is the difference between a denormal guard and an
    // audible one. `env` holds a POWER there, so the house 1e-15 zaps an AMPLITUDE of 3.2e-8, i.e.
    // -150 dBFS — and because the flush fires once per process() call, whether it fired at all then
    // depended on how the caller happened to cut the stream into blocks. Measured on a -160 dBFS input:
    // -7.5 dB of gain reduction in one 10000-sample call against 0.00 dB in 10000 one-sample calls, the
    // same stream. 1e-30 of power is an amplitude of 1e-15, i.e. -300 dBFS — the level the house
    // constant was chosen to mean — and still 8 orders above the float subnormal range.
    // BE PRECISE ABOUT WHAT THAT BUYS. A threshold flush is clocked by the caller's blocks by
    // construction, so it cannot make the state partition-independent; it can only put the floor where
    // nothing downstream can see it. `core::gainToDb` clamps its argument at 1e-12 (-240 dBFS), so a
    // level below the new floor cannot reach a gain at all — what is left is visible only to a caller
    // reading the raw envelope, which is why the compressor's own output no longer depends on the
    // partition and the meter, below -280 dBFS, still can.
    void flushDenormals() noexcept
    {
        if (det == Detector::Rms) { if (env < 1.0e-30f) env = 0.0f; }
        else core::flushDenormal (env);
    }

private:
    void  updateCoeffs() noexcept { atkCoeff = coeff (atkMs); relCoeff = coeff (relMs); }
    float coeff (double ms) const noexcept
    {
        const double t = ms * 0.001;
        // !(t>0) catches a NaN time; the isfinite pair catches the two that slip past it and are NOT
        // harmless: +Inf ms gives exp(-1/Inf) == 1.0, a coefficient that FREEZES the envelope forever
        // (`env = in + 1*(env-in)` never moves), and a non-finite rate leaves a NaN coeff. Both become
        // "instant", which is the same fallback every other non-finite time already took.
        if (! (t > 0.0) || ! std::isfinite (t) || ! (fs > 0.0) || ! std::isfinite (fs)) return 0.0f;
        const float c = (float) std::exp (-1.0 / (t * fs));
        if (c < 1.0e-15f) return 0.0f;      // absurdly short time → instant; never leave a subnormal coeff
        // ...and never leave a coefficient of exactly 1, which is not "very slow" but FROZEN: `env = in +
        // 1*(env-in)` is `env = env`, so the envelope stops tracking for good. A time constant only has
        // to reach ~1e6 ms at 48 kHz for exp() to round to 1.0f. The largest float below 1 keeps it a
        // (very) slow one-pole instead of a dead one, and is bit-transparent for every c that is < 1.
        return c < 1.0f ? c : 0x1.fffffep-1f;
    }

    double   fs = 48000.0, atkMs = 10.0, relMs = 100.0;
    float    atkCoeff = 0.0f, relCoeff = 0.0f;
    float    env = 0.0f;            // peak: |x| envelope · rms: mean-square accumulator
    Detector det = Detector::Rms;
};

} // namespace felitronics::dynamics

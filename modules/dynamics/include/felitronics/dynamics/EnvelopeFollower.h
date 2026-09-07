// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/core/FlushToZero.h>
#include <felitronics/core/Math.h>

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

    // THE STORED WORD, not the reported amplitude — the two are different quantities in Rms mode, where
    // this holds a POWER. Law 11c's silence advancement decides it has reached a fixed point by comparing
    // this across a step, and the getter above cannot do that job: a power moving from 0x3f800001 to
    // 0x3f800000 has both roots at exactly 1.0f, so a getter-based check would declare a state settled
    // while it is still travelling. Exposed for that one purpose; nothing in a signal path should read it.
    float stateWord() const noexcept { return env; }

    // LAW 11c — advance through `n` samples of DIGITAL SILENCE, sample for sample, the way `process(0)`
    // would. This is not a closed form and deliberately not one: `pow(c, n)` is a DIFFERENT number from
    // `n` rounded multiplications, so it would buy speed with the bit-exactness the law is stated in.
    // What makes it cheap instead is that the silent recurrence is AUTONOMOUS, so it reaches a bitwise
    // fixed point and everything past that point is free.
    //
    // TWO BRANCHES COLLAPSE HERE, and both are proven rather than assumed. `in` is `|0|` in Peak mode and
    // `0*0` in Rms mode, i.e. exactly +0.0f either way; `in > env` is therefore false for every reachable
    // `env`, because `env` is non-negative by induction (it starts at 0, and `in + c*(env-in)` with
    // `in >= 0`, `env >= 0`, `c` in [0,1] cannot go below 0). So the release coefficient is chosen every
    // sample — and in Rms mode LinkedDetector sets attack == release anyway, so the choice does not even
    // matter there.
    //
    // WHAT IT DOES NOT REACH IS ZERO, and the earlier draft of this comment claimed it did. Repeated
    // multiplication by `c < 1` stalls on a SUBNORMAL: measured, at `c = 0x1.fffffep-1` the states
    // 2^-149, 3*2^-149, 2^-126 - 2^-149 and 2^-126 (the smallest NORMAL, where the decrement is exactly
    // half a subnormal ULP and ties-to-even rounds back up) are all fixed points, and at realistic
    // coefficients the resting state is a subnormal too — 0x00000078 for a 5 ms release, 0x00005d9f for
    // 1 s. What turns those into a real zero is `flushDenormals()`, once per call, on both sides of the
    // comparison. And the horizon is a property of the TIME CONSTANT, not of the pause: measured 23 609
    // samples to settle at 5 ms, 457 808 at 100 ms, 4 461 677 at 1 s, and NOT settled after 200 000 000
    // at the coefficient cap. "A long gap costs the same as a short one" is true only past that horizon.
    void advanceSilence (int n) noexcept
    {
        for (int i = 0; i < n; ++i)
        {
            const float before = env;
            (void) process (0.0f);
            if (core::sameBits (env, before)) return;                  // a fixed point stays fixed
            // ...and a state that is not finite cannot come back: c*NaN is NaN, and c*Inf is Inf or (at
            // c == 0) NaN. Two non-finite states in a row is that fact, and it stops the loop from
            // spinning for the whole pause on something that will never move again.
            if (! std::isfinite (env) && ! std::isfinite (before)) return;
        }
    }

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
    // reading the raw envelope, and the meter below -280 dBFS still shows it.
    // THIS PARAGRAPH USED TO END "...which is why the compressor's own output no longer depends on the
    // partition." IT IS FALSE, and the counterexample is two samples long. The flush zeroes the state,
    // but the trajectory OUT of zero and out of a surviving 8.45e-31 then climbs back ABOVE the 1e-12
    // floor, where the clamp no longer hides the difference. With the window set so the coefficient is
    // exactly 0.5f (rmsWindowMs = 1000/(fs·ln2)) and a key of [1.3e-15f, 1.5e-12f], the amplitude after
    // the second sample is 1.060660367e-12 carried against 1.060660168e-12 flushed; through a threshold
    // of -240 dB at 4:1 on a programme of 0.5f that is 0.478396237 in one 2-sample call against
    // 0.478396297 in two 1-sample calls — half an LSB of 24-bit, from a floor that was supposed to make
    // it unreachable. What the floor buys is that ordinary material does not reach the window between
    // 1e-15 and 1e-12; it is not a proof that nothing can. A caller that needs the output to be
    // partition-independent has to stop the caller's blocks reaching the stage at all — which is what
    // `felitronics::mastering` does with a fixed internal quantum, and why it has one.
    // POISON IS CLEARED TOO, and in Rms mode NOT by swapping in `core::flushPoison`. That would look
    // like the obvious edit and would silently undo the paragraph above: `flushPoison` carries the
    // house 1e-15 threshold, which in the POWER domain means an amplitude of 3.2e-8, i.e. -150 dBFS —
    // the very bug this comment exists to record. The non-finite half is what has to be added; the
    // 1e-30 threshold is what has to survive.
    void flushDenormals() noexcept
    {
        if (det == Detector::Rms) { if (! std::isfinite (env) || env < 1.0e-30f) env = 0.0f; }
        else core::flushPoison (env);
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

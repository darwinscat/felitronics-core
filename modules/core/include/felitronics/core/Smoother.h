// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <algorithm>
#include <cmath>
#include <limits>

namespace felitronics::core
{

//==============================================================================
// One-pole exponential parameter smoother (RT-safe, JUCE-free). `next()` advances one sample;
// `advance(n)` jumps n samples in closed form so a consumer can recompute coefficients once per block
// while the underlying parameter still moves in real time (no zipper). Migrated from teq::Smoother —
// the de-JUCE replacement for juce::SmoothedValue / juce::LinearSmoothedValue (note: this glide is
// EXPONENTIAL, not linear → swapping a JUCE linear smoother for it is a versioned-behaviour change).
class Smoother
{
public:
    void prepare (double sampleRate, double timeMs = 30.0) noexcept { fs = sampleRate; setTimeMs (timeMs); }

    void setTimeMs (double timeMs) noexcept
    {
        const double tau = timeMs * 0.001;
        coeff = (! (tau > 0.0) || ! (fs > 0.0)) ? 0.0 : std::exp (-1.0 / (tau * fs));   // !(x>0) also catches NaN → instant, never a NaN coeff
    }

    void   setTarget (double t) noexcept { target = t; }
    void   snap (double v) noexcept { target = current = v; }
    double value() const noexcept { return current; }
    double targetValue() const noexcept { return target; }
    bool   settled (double eps = 1e-7) const noexcept { return std::fabs (current - target) <= eps; }

    inline double next() noexcept { current = target + coeff * (current - target); return current; }

    inline double advance (int n) noexcept
    {
        // Settled early-out: a snapped/converged smoother (current == target) skips the pow — so an
        // idle parameter costs one compare per block, not a transcendental (the eq lanes' cost-zero
        // contract rides on this).
        if (n > 0 && std::fabs (current - target) > 0.0)
            current = target + std::pow (coeff, (double) n) * (current - target);
        return current;
    }

private:
    double fs = 0.0, coeff = 0.0, current = 0.0, target = 0.0;
};

//==============================================================================
// One-pole? No — a fixed-increment LINEAR ramp. JUCE-free, RT-safe, bit-for-bit drop-in for
// juce::SmoothedValue<float, ValueSmoothingTypes::Linear>: same method names, same float numerics, the
// same `floor()` ramp-length, and the same decrement-then-step order — so a consumer swaps the TYPE with
// no call-site changes and the host-facing feel (mute/gain/mix ramps) is identical. A new target restarts
// a full `stepsToTarget`-sample ramp from the CURRENT value. Use the exponential Smoother above for
// coefficient glides; use this where a JUCE linear smoother is being replaced.
class LinearSmoother
{
public:
    LinearSmoother() noexcept = default;
    LinearSmoother (float initialValue) noexcept : currentValue (initialValue), target (initialValue) {}   // non-explicit: juce::SmoothedValue is too, so `LinearSmoother x { 1.0f }` array-init stays a drop-in

    // Ramp length from seconds. NOTE: floor(), matching JUCE exactly (not round()).
    void reset (double sampleRate, double rampLengthInSeconds) noexcept
    {
        reset ((int) std::floor (rampLengthInSeconds * sampleRate));
    }

    // Ramp length directly in samples. Snaps current→target and stops any ramp in progress.
    void reset (int numSteps) noexcept
    {
        stepsToTarget = numSteps;
        setCurrentAndTargetValue (target);
    }

    // Matches juce::approximatelyEqual EXACTLY (abs tol = smallest normal, rel tol = machine epsilon; non-finite
    // → exact ==), so a ~1-ULP automation-jitter target is a no-op just like JUCE — no spurious ramp restart /
    // zipper the original never had. A raw `==` would restart a full ramp on a difference JUCE treats as equal.
    static bool approxEqual (float a, float b) noexcept
    {
        if (! (std::isfinite (a) && std::isfinite (b))) return a == b;
        const float diff = std::fabs (a - b);
        return diff <= std::numeric_limits<float>::min()
            || diff <= std::numeric_limits<float>::epsilon() * std::max (std::fabs (a), std::fabs (b));
    }

    void setTargetValue (float newValue) noexcept
    {
        if (approxEqual (newValue, target)) return;
        if (stepsToTarget <= 0) { setCurrentAndTargetValue (newValue); return; }
        target    = newValue;
        countdown = stepsToTarget;
        step      = (target - currentValue) / (float) countdown;
    }

    void setCurrentAndTargetValue (float newValue) noexcept
    {
        target = currentValue = newValue;
        countdown = 0;
    }

    // One sample. Decrement-then-step, snapping to target on the final step (no float drift past target).
    float getNextValue() noexcept
    {
        if (countdown <= 0) return target;
        --countdown;
        if (countdown > 0) currentValue += step;
        else               currentValue = target;
        return currentValue;
    }

    // Closed-form jump of n samples (== n getNextValue() calls), for once-per-block recompute.
    float skip (int numSamples) noexcept
    {
        if (numSamples >= countdown) { setCurrentAndTargetValue (target); return target; }
        currentValue += step * (float) numSamples;
        countdown    -= numSamples;
        return currentValue;
    }

    void applyGain (float* samples, int numSamples) noexcept
    {
        if (countdown > 0) for (int i = 0; i < numSamples; ++i) samples[i] *= getNextValue();
        else               for (int i = 0; i < numSamples; ++i) samples[i] *= target;
    }

    float getCurrentValue() const noexcept { return currentValue; }
    float getTargetValue()  const noexcept { return target; }
    bool  isSmoothing()     const noexcept { return countdown > 0; }

private:
    float currentValue = 0.0f, target = 0.0f, step = 0.0f;
    int   countdown = 0, stepsToTarget = 0;
};

} // namespace felitronics::core

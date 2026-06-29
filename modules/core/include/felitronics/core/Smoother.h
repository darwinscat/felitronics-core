// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa. Part of felitronics-core — see LICENSE.

#pragma once

#include <cmath>

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
        coeff = (tau <= 0.0 || fs <= 0.0) ? 0.0 : std::exp (-1.0 / (tau * fs));
    }

    void   setTarget (double t) noexcept { target = t; }
    void   snap (double v) noexcept { target = current = v; }
    double value() const noexcept { return current; }
    double targetValue() const noexcept { return target; }
    bool   settled (double eps = 1e-7) const noexcept { return std::fabs (current - target) <= eps; }

    inline double next() noexcept { current = target + coeff * (current - target); return current; }

    inline double advance (int n) noexcept
    {
        if (n > 0) current = target + std::pow (coeff, (double) n) * (current - target);
        return current;
    }

private:
    double fs = 0.0, coeff = 0.0, current = 0.0, target = 0.0;
};

} // namespace felitronics::core

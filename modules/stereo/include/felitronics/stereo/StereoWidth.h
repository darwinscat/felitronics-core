// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/stereo/MidSide.h>
#include <felitronics/core/Smoother.h>

#include <algorithm>
#include <cmath>

namespace felitronics::stereo
{

//==============================================================================
// felitronics::stereo::StereoWidth — broadband stereo-image width via a Mid/Side side-gain.
//
//   M = ½(L+R),  S = ½(L-R)
//   L' = g·(M + width·S),   R' = g·(M − width·S)        (g = outputGain, default 1)
//
//   width = 0 → mono (S removed),  1 → neutral,  2 → wide (+6 dB side). Clamped to [0, 2].
//
// THE MONO-FOLD INVARIANT (the reason no "mono-safe" limiter is needed): the tool
// touches ONLY the Side, so at outputGain=1 the Mid is mathematically untouched — the mono sum ½(L'+R') = M
// is preserved (to float rounding; the two halves round separately). Widening therefore can NEVER weaken the
// mono fold below the source; it only changes how much of the (already-present) stereo difference survives
// the fold. So we do NOT second-guess the engineer with a correlation-driven side limiter (that pumps and is
// content-dependent — rejected for a mastering tool); we hard-clamp and let a host CorrelationMeter do the
// metering. A correlation-limited "mono-safe" mode was considered and declined for this reason.
//
// width/outputGain are LINEAR-smoothed (core::LinearSmoother, ~20 ms) so live automation doesn't click — a
// width step is otherwise an instant Δwidth·S discontinuity. Linear (not exponential) so a settled value
// reaches its target EXACTLY: at width=1, gain=1 the round-trip is skipped → BIT-EXACT neutral passthrough,
// zero latency. During a width ramp at gain=1 the Mid is still untouched per sample, so the mono-fold
// invariant holds throughout the move. Mono in (L==R → S=0) stays mono for any width. Setters reject
// non-finite values (a stray host NaN can't poison the stream). Aliased io[0]==io[1] applies gain once.
//
// Frequency-dependent width (low-narrow / high-wide) is deliberately left to composition —
// multiband::MultibandProcessor<StereoWidth> over the allpass-summed LR4 split reconstructs cleanly (the
// bands sum in-phase at each crossover, so a per-band side-gain is smooth, not combed) — built later;
// mono-below-fc stays its own primitive (MonoBass). RT-safe: no alloc/lock/throw; the only state is the two
// smoothers (no feedback → nothing to denormal-flush).
class StereoWidth
{
public:
    static constexpr float  kMaxWidth     = 2.0f;    // +6 dB side — the mastering-grade ceiling
    static constexpr double kSmoothingMs  = 20.0;    // click-free width/gain automation

    void prepare (double sampleRate, int /*maxBlock*/ = 0, int /*maxChannels*/ = 2) noexcept
    {
        fs_ = sampleRate > 0.0 ? sampleRate : 48000.0;
        widthSm_.reset (fs_, kSmoothingMs * 0.001); gainSm_.reset (fs_, kSmoothingMs * 0.001);
        reset();
    }

    void reset() noexcept                            // snap the smoothers to their targets (settled, no glide)
    {
        widthSm_.setCurrentAndTargetValue (width_);
        gainSm_.setCurrentAndTargetValue (gain_);
    }

    void setEnabled (bool e) noexcept { enabled_ = e; }
    void setWidth (float w) noexcept
    {
        if (! std::isfinite (w)) return;             // reject NaN/inf — keep the last good value
        width_ = std::clamp (w, 0.0f, kMaxWidth);
        widthSm_.setTargetValue (width_);
    }
    void setOutputGain (float g) noexcept
    {
        if (! std::isfinite (g)) return;
        gain_ = std::clamp (g, 0.0f, 4.0f);          // manual loudness trim (no auto-comp)
        gainSm_.setTargetValue (gain_);
    }

    bool  isEnabled()  const noexcept { return enabled_; }
    float width()      const noexcept { return width_; }
    float outputGain() const noexcept { return gain_; }
    static constexpr int latencySamples() noexcept { return 0; }

    // Stereo, in place: io[0]=L, io[1]=R. Bypass (untouched) when disabled, < 2 channels, or settled-neutral.
    void process (float* const* io, int numChannels, int n) noexcept
    {
        if (numChannels < 2 || ! enabled_) return;
        if (! widthSm_.isSmoothing() && ! gainSm_.isSmoothing()
            && widthSm_.getCurrentValue() == 1.0f && gainSm_.getCurrentValue() == 1.0f)
            return;                                  // settled neutral → bit-exact, skip the round-trip
        float* L = io[0];
        float* R = io[1];
        for (int i = 0; i < n; ++i)
        {
            const float w = widthSm_.getNextValue();
            const float g = gainSm_.getNextValue();
            float m, s; MidSide::encode (L[i], R[i], m, s);
            MidSide::decode (g * m, g * (w * s), L[i], R[i]);   // gain folded in → M untouched at g=1; alias-safe
        }
    }

private:
    double fs_ = 48000.0;
    float  width_ = 1.0f, gain_ = 1.0f;
    bool   enabled_ = true;
    core::LinearSmoother widthSm_ { 1.0f }, gainSm_ { 1.0f };
};

} // namespace felitronics::stereo

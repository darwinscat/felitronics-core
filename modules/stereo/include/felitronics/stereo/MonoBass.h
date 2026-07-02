// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/stereo/MidSide.h>
#include <felitronics/eq/Crossover2.h>
#include <felitronics/core/Smoother.h>

#include <algorithm>
#include <cmath>

namespace felitronics::stereo
{

//==============================================================================
// felitronics::stereo::MonoBass — "elliptical EQ" / bass mono-maker: collapse the low end to mono below a
// crossover (vinyl/LP lacquer cutting — out-of-phase lows cause vertical stylus excursion — and general
// mastering tightening). Zero added latency, minimum-phase.
//
// TOPOLOGY (verified by math + the self-tests). Operate ONLY on the Side of a Mid/Side split, so the Mid
// (the mono content you keep) passes UNFILTERED — no transient smear on the kept bass; only the
// low-frequency Side (stereo) energy is removed. A full L/R crossover would needlessly filter the mono
// bass too. Below fc the output is L == R BY CONSTRUCTION (the side goes to zero), so the tool can only
// IMPROVE mono compatibility, never harm it.
//
//   M = ½(L+R),  S = ½(L-R)
//   S_wet = lowWidth·LP4(S) + HP4(S)      // LP4/HP4 = eq::Crossover2, a 4th-order Linkwitz-Riley split
//   S'    = xf·S + (1-xf)·S_wet           // xf = the full-wide bypass crossfade (see PARAMS below)
//   L' = M + S',  R' = M − S'
//
// WHY LINKWITZ-RILEY, fixed at LR4 / 24 dB/oct — consilium-settled; don't re-litigate without new math:
//  * The LR4 identity LP4+HP4 = allpass — (1+s⁴)/(s²+√2s+1)², |·| = (1+ω⁴)/(1+ω⁴) = 1, the branches
//    exactly IN-PHASE — is what makes the `lowWidth` taper bump-free: the blended side magnitude is
//    (lowWidth + r⁴)/(1 + r⁴) with r = tan(πf/fs)/tan(πfc/fs) (the BLT warp the TPT SVF realises exactly).
//  * matched::highpass/lowpass (the Nyquist-matched EQ filters) were REJECTED: they are magnitude fits,
//    not a complementary pair — their sum is not that allpass, so the blend would ripple; and Nyquist
//    accuracy buys nothing at a ~120 Hz crossover.
//  * An allpass-subtraction variant was REJECTED: AP4−LP4 has magnitude 1.5 at fc (a boost, not a reject).
//  * A crossfade blend w·S + (1−w)·HP4(S) was REJECTED: S and HP4(S) are anti-phase at fc (HP4(jω₀) = −½),
//    so that blend NULLS the side at fc when w = ⅓.
//  * LR2 / 12 dB/oct was REJECTED: LP2+HP2 = (1+s²)/(s+1)² has a NULL at fc — its flat sum needs a
//    polarity flip (LP−HP), which breaks S' = w·LP + HP. LR8 / 48 dB/oct would be admissible (the in-phase
//    allpass sum holds) but is YAGNI until a product asks for it.
//
// PARAMS. `lowWidth` ∈ [0,1] — 0 = fully mono below fc (default), 1 = full stereo — is LINEAR-smoothed
// (~20 ms) so live automation doesn't click. Full-wide is special: the wet path at lowWidth=1 is the LR4
// allpass — magnitude-flat but phase-rotated (−1 at fc) — so a hard switch to raw bypass there would
// click. Instead lowWidth=1 also ramps the dry/wet crossfade `xf`→1 (same 20 ms), and only once xf has
// SETTLED at 1.0 does process() take the bit-exact early-return bypass. Entering ANY bypass (settled
// full-wide / disabled / non-stereo) resets the crossover state, so re-entry never replays stale tails —
// it restarts the filters from zero, a bounded click-free settle (tested). `setEnabled` stays a hard
// toggle (host-level bypasses ramp; sibling StereoWidth parity). Setters reject non-finite values and
// clamp (a stray host NaN can't poison the SVF state — std::clamp would pass it through). reset() snaps
// the smoothers to their targets (snap-on-load: no ramp on session recall). `frequency()` reads back the
// clamped value ([20, 0.45·fs], re-clamped if prepare() lowers fs).
//
// STRICTLY STEREO: process() touches the buffer ONLY when numChannels == 2 — a mono or surround bus
// passes through whole (treating a stereo pair inside a wider layout is a routing decision the host
// owns, not this class). RT-safe: no alloc/lock/throw in process(); state = one eq::Crossover2 (the Side
// is one channel) + two LinearSmoothers.
class MonoBass
{
public:
    static constexpr double kSmoothingMs = 20.0;    // click-free lowWidth automation + the bypass crossfade
    static constexpr float  kMinFreq     = 20.0f;

    void prepare (double sampleRate, int /*maxBlock*/ = 0, int /*maxChannels*/ = 2) noexcept
    {
        fs_ = sampleRate > 0.0 ? sampleRate : 48000.0;
        xo_.prepare (fs_, 1);
        widthSm_.reset (fs_, kSmoothingMs * 0.001);
        xfSm_.reset (fs_, kSmoothingMs * 0.001);
        applyFrequency();
        reset();
    }

    void reset() noexcept                           // snap smoothers to targets (settled, no glide) + clear filter state
    {
        widthSm_.setCurrentAndTargetValue (lowWidth_);
        xfSm_.setCurrentAndTargetValue (lowWidth_ >= 1.0f ? 1.0f : 0.0f);
        xo_.reset();
        bypassed_ = false;
    }

    void setEnabled (bool e) noexcept { enabled_ = e; }

    void setFrequency (float hz) noexcept
    {
        if (! std::isfinite (hz)) return;           // reject NaN/inf — keep the last good value
        freq_ = hz;
        applyFrequency();
    }

    void setLowWidth (float w) noexcept             // 0 = mono below fc, 1 = full stereo (settles into bypass)
    {
        if (! std::isfinite (w)) return;
        lowWidth_ = std::clamp (w, 0.0f, 1.0f);
        widthSm_.setTargetValue (lowWidth_);
        xfSm_.setTargetValue (lowWidth_ >= 1.0f ? 1.0f : 0.0f);
    }

    bool  isEnabled() const noexcept { return enabled_; }
    float frequency() const noexcept { return freq_; }
    float lowWidth()  const noexcept { return lowWidth_; }
    static constexpr int latencySamples() noexcept { return 0; }

    // Stereo, in place: io[0]=L, io[1]=R. Bypass (buffer untouched) when disabled, numChannels != 2, or
    // settled full-wide. Entering bypass clears the crossover so re-entry starts from silence, not stale tails.
    void process (float* const* io, int numChannels, int n) noexcept
    {
        if (numChannels != 2 || ! enabled_ || (! xfSm_.isSmoothing() && xfSm_.getCurrentValue() == 1.0f))
        {
            if (! bypassed_) { bypassed_ = true; xo_.reset(); }
            return;
        }
        bypassed_ = false;
        float* L = io[0];
        float* R = io[1];
        for (int i = 0; i < n; ++i)
        {
            const float w  = widthSm_.getNextValue();
            const float xf = xfSm_.getNextValue();
            float m, s; MidSide::encode (L[i], R[i], m, s);
            float lp, hp; xo_.processSample (0, s, lp, hp);
            const float wet  = w * lp + hp;                  // side magnitude (w + r⁴)/(1+r⁴) — bump-free (LR4 in-phase)
            const float sOut = xf * s + (1.0f - xf) * wet;   // ≠ wet only while fading into/out of full-wide
            MidSide::decode (m, sOut, L[i], R[i]);
        }
        xo_.flushDenormals();
    }

private:
    void applyFrequency() noexcept                  // coefficients only — no state reset, so a freq move doesn't click
    {
        freq_ = std::clamp (freq_, kMinFreq, (float) (0.45 * fs_));
        xo_.setFrequency (freq_);
    }

    double fs_ = 48000.0;
    float  freq_ = 120.0f, lowWidth_ = 0.0f;
    bool   enabled_ = true, bypassed_ = false;
    eq::Crossover2 xo_;                             // the Side-channel LR4 split (the primitive extracted from here, reused back)
    core::LinearSmoother widthSm_ { 0.0f }, xfSm_ { 0.0f };
};

} // namespace felitronics::stereo

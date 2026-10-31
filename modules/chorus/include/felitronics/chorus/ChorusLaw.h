// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/chorus/BassChorus.h>
#include <felitronics/core/DetMath.h>

#include <algorithm>
#include <cmath>

namespace felitronics::chorus
{

//==============================================================================
// THE ONE-SLIDER LAW — what the bass amp's CHORUS slider does: `amount` 0 … 1 → the stage's knobs.
//
// The slider owns TWO fields, `depthMs` and `mix`, and drives them together, linearly: depth 0 … 3 ms of head
// swing and the voice's level 0 … 1 (1 = the dry's level, the classic 50/50). The wobble (∝ depth) and the
// voice's level both read evenly along the travel. Amount 0 is mix 0, which idles the stage: the dry bit for
// bit, nothing computed.
//
// THE MENU owns the other three — `rateHz`, `crossoverHz`, `delayMs` — and the slider never moves them: a menu
// RATE then means the same thing at every slider position, and the wobble scales with it as a pedal's does
// (faster is more pitch). The menu must not expose depth or mix; they are the slider's.
//
// A non-finite amount is IGNORED — on one condition the caller keeps: it passes the stage's CURRENT knobs,
// `stage.params()`, as `menu`. Then `menu` comes back untouched and `setParams (paramsForAmount (NaN,
// stage.params()))` changes nothing. A menu kept apart from the stage, or the default `{}`, would come back with
// ITS depth and mix, which is not "ignored"; so the slider's call is always the one above. A finite amount is
// clamped to 0 … 1.
//==============================================================================

inline constexpr float kDepthMsAtFull = 3.0f;   ///< the swing at amount 1: 23 cents peak at the 0.7 Hz default
inline constexpr float kMixAtFull     = 1.0f;   ///< the voice at the dry's level at amount 1

static_assert (kDepthMsAtFull + BassChorus::kHeadroomMs + BassChorus::kCapKneeMs <= BassChorus::kMinDelayMs,
               "the centre's floor must keep the slider's full depth out of the smooth minimum's knee");

/** The slider's law: `amount` sets `depthMs` and `mix`; every other field is `menu`'s. */
inline BassChorusParams paramsForAmount (float amount, BassChorusParams menu = {}) noexcept
{
    if (! std::isfinite (amount))
        return menu;
    const float a = detail::clampf (0.0f, 1.0f, amount);
    menu.depthMs = kDepthMsAtFull * a;
    menu.mix     = kMixAtFull * a;
    return menu;
}

/** The voice's peak pitch deviation for these knobs, in cents — the LARGER of the two excursions, which is the
    flat one: the head's speed swings the frequency ratio by ±x (x = 2π·rate·swing), and 1200·log2(1 - x) is
    farther from 0 than 1200·log2(1 + x) (at 10 Hz x 10 ms: -1713 against +844). The swing is the depth as the
    stage holds it, eased under the centre minus the headroom. Knobs are bounded as the stage bounds them; a
    non-finite one reads as no wobble. A pure function: any thread may ask it. */
inline float pitchDeviationCents (const BassChorusParams& p) noexcept
{
    if (! (std::isfinite (p.depthMs) && std::isfinite (p.rateHz) && std::isfinite (p.delayMs)))
        return 0.0f;
    const double depth = (double) detail::clampf (0.0f, BassChorus::kMaxDepthMs, p.depthMs);
    const double rate  = (double) detail::clampf (BassChorus::kMinRateHz, BassChorus::kMaxRateHz, p.rateHz);
    const double delay = (double) detail::clampf (BassChorus::kMinDelayMs, BassChorus::kMaxDelayMs, p.delayMs);
    const double swing = detail::smoothMin (depth, delay - (double) BassChorus::kHeadroomMs, (double) BassChorus::kCapKneeMs);
    const double x = std::min (detail::kTwoPi * rate * swing * 0.001, 0.999);
    return (float) (-1200.0 * core::det::log2 (1.0 - x));
}

} // namespace felitronics::chorus

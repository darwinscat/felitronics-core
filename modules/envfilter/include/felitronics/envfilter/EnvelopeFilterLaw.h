// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/envfilter/EnvelopeFilter.h>

#include <cmath>

namespace felitronics::envfilter
{

//==============================================================================
// THE ONE-SLIDER LAW — what a face's single FILTER slider means. `amount` 0..1 is "how much wah the playing
// gets", and it moves exactly two knobs:
//
//   · MIX rises with the slider's lower half and is all filter from kAmountFullWet on. Mix is what brings the
//     effect IN: at a quarter the filter is heard beside the dry note rather than instead of it. It is safe to
//     blend because the stage's dry and its wet agree wherever the filter does nothing (the clean lows are a
//     complement of the input, not an allpassed band) — a half mix does not notch the bass.
//   · SENSITIVITY climbs the whole travel, kSensAtZeroDb → kSensAtOneDb: how far a given hit opens the sweep.
//     Low on the slider only the hard notes speak (a -12 dBFS note opens 60 % of the way at a quarter); at the
//     top every fingered note opens it all the way and closes it as it decays.
//
// Everything else — the mode, the resonance, the direction, the range, the follower's times and the clean
// lows — is the VOICE, and lives in the menu: it passes through from `menu` untouched. Range is deliberately
// not on the slider: a sweep that also shrank would make the low half of the slider sound like a different
// filter rather than less of this one.
//
// amount 0 gives mix 0, which the stage answers with the buffer untouched and nothing computed (bit-exact
// dry); a non-finite amount counts as 0.
//==============================================================================

inline constexpr float kAmountFullWet = 0.5f;    ///< from here up the slider is all filter
inline constexpr float kSensAtZeroDb  = -6.0f;   ///< sensitivity at the bottom of the travel
inline constexpr float kSensAtOneDb   = 18.0f;   ///< ...and at the top

/** The slider's law: `amount` (clamped to 0..1) sets `mix` and `sensitivityDb`; every other field is `menu`'s. */
inline EnvelopeFilterParams paramsForAmount (float amount, EnvelopeFilterParams menu = {}) noexcept
{
    const float a = std::isfinite (amount) ? (amount < 0.0f ? 0.0f : (amount > 1.0f ? 1.0f : amount)) : 0.0f;
    menu.mix           = a >= kAmountFullWet ? 1.0f : a / kAmountFullWet;
    menu.sensitivityDb = kSensAtZeroDb + (kSensAtOneDb - kSensAtZeroDb) * a;
    return menu;
}

} // namespace felitronics::envfilter

// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/octave/Octaver.h>

#include <cmath>
#include <limits>

namespace felitronics::octave
{

//==============================================================================
// THE ONE-SLIDER LAW — what the bass amp's OCTAVE slider does: `amount` 0 … 1 → the stage's knobs.
//
// The slider is a MASTER over the octaves the menu has voiced. It OWNS `sub1` and `sub2`: each is the menu's own
// value of that field — the level it reaches at the slider's top — times one taper, c(a) = a/2 + a²/2. The dry
// never moves (the stage has no dry level). Amount 0 is c = 0, both levels an exact zero whatever the menu says,
// and the stage goes idle: the dry bit for bit, nothing computed (Octaver.h, IDLE).
//
// The taper with the menu's defaults (−1 at 1.25, −2 at 0): a quarter of the slider is a hint of weight (−14.2 dB re
// the band's peak), half is plainly there (−6.6 dB), the top goes past unity (+1.9 dB), because an octave down at an
// equal peak SOUNDS weaker than the note: the ear's sensitivity falls steeply below 50 Hz, where the sub of every
// open string lives.
//
// THE −2 OCTAVE IS OFF IN THE MENU'S DEFAULT. On a bass it lands at 7.7 Hz (low B) to 24.5 Hz (G2): below hearing
// or at its edge. What is heard of it is its square's upper partials, while its fundamental spends the headroom and
// the speaker's excursion. A menu that sets `sub2` gets it under the same master.
//
// THE MENU owns `toneHz`, `subHpHz` and `gateDb`; they pass through untouched. The amp's menu shows `sub1` and `sub2`
// as "the −1 / −2 level at the top of the slider" — never as the stage's current level.
//
// A non-finite amount is IGNORED: the two owned fields come back as NaN, which Octaver's setters ignore, so the
// stage's levels stand where they were; the menu's fields still pass. A finite amount is clamped to 0 … 1.
//==============================================================================

/** The slider's taper, 0 … 1 → 0 … 1: half linear, half square. */
inline float sliderTaper (float a) noexcept { return 0.5f * a + 0.5f * a * a; }

/** The slider's law: `amount` scales `menu.sub1` and `menu.sub2`; every other field is `menu`'s. */
inline OctaverParams paramsForAmount (float amount, OctaverParams menu = {}) noexcept
{
    if (! std::isfinite (amount))
    {
        menu.sub1 = menu.sub2 = std::numeric_limits<float>::quiet_NaN();
        return menu;
    }

    const float c = sliderTaper (detail::clampf (0.0f, 1.0f, amount));
    const auto  owned = [c] (float top) { return std::isfinite (top) ? detail::clampf (0.0f, Octaver::kMaxSub, top) * c + 0.0f : 0.0f; };
    menu.sub1 = owned (menu.sub1);
    menu.sub2 = owned (menu.sub2);
    return menu;
}

} // namespace felitronics::octave

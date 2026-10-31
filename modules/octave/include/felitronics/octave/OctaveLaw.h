// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/octave/Octaver.h>

#include <cmath>

namespace felitronics::octave
{

//==============================================================================
// THE ONE-SLIDER LAW — what the bass amp's OCTAVE slider does: `amount` 0 … 1 → the stage's knobs.
//
// The −1 octave rises under the dry, which stays at unity; nothing else moves. 0 is an exact zero, so the stage is
// at rest and the dry passes bit for bit (Octaver.h, AT REST). The taper is half linear, half square — a quarter of
// the slider is a hint of weight (−14 dB re the band's peak), half is plainly there (−6.6 dB), and the top goes past
// unity (+1.9 dB), because an octave down at an equal peak SOUNDS weaker than the note: the ear's sensitivity falls
// steeply below 50 Hz, and that is where the sub of every open string lives.
//
// THE −2 OCTAVE STAYS OUT OF THE SLIDER. On a bass it lands at 7.7 Hz (low B) to 24.5 Hz (G2): below hearing or
// at its edge. What is heard of it is its square's upper partials, while its fundamental spends the headroom and
// the speaker's excursion — on a one-slider face that is mud with a cost. It is there in the menu (`sub2`), and the
// renders let the owner hear it alone.
//==============================================================================

/** The top of the slider's −1 octave level. */
inline constexpr float kSliderTopSub1 = 1.25f;

/** The OCTAVE slider: `amount` 0 … 1 (clamped; a non-finite one reads as 0) → the knobs. Tone and gate keep their
    defaults — they belong to the menu. */
inline OctaverParams paramsForAmount (float amount) noexcept
{
    const float a = std::isfinite (amount) ? detail::clampf (0.0f, 1.0f, amount) : 0.0f;

    OctaverParams p;
    p.sub1 = kSliderTopSub1 * (0.5f * a + 0.5f * a * a);
    p.sub2 = 0.0f;
    p.dry  = 1.0f;
    return p;
}

} // namespace felitronics::octave

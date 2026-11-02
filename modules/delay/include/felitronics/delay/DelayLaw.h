// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/core/DetMath.h>

#include <cmath>

namespace felitronics::delay
{

//==============================================================================
// The delay's knob law — the numbers the amps' faces turn into what `TapeDelay` takes, so a plugin and the
// browser demo land on the same values: the sync ladder, the sync policy, the DARK knob, the knob ranges.
// Deterministic (core::det), so every row computes the same time and the same corner from the same knob.
//==============================================================================

/** The sync ladder: the plain values from the longest (1/1) to the shortest (1/16), each with its dotted
    value just before it and its triplet just after. Not monotonic — 1/2T (1.33 beats) sits before 1/4. (1.5). */
inline constexpr int kNumDivisions = 13;
inline constexpr const char* kDivisionNames[kNumDivisions] = { "1/1", "1/2.", "1/2", "1/2T", "1/4.", "1/4",
                                                               "1/4T", "1/8.", "1/8", "1/8T", "1/16.", "1/16",
                                                               "1/16T" };
/** The same ladder in quarter notes: dotted = 1.5 × plain, triplet = ⅔ × plain; ms = beats × 60000 / BPM. */
inline constexpr float kDivisionBeats[kNumDivisions] = { 4.0f, 3.0f, 2.0f, 4.0f / 3.0f, 1.5f, 1.0f,
                                                         2.0f / 3.0f, 0.75f, 0.5f, 1.0f / 3.0f, 0.375f,
                                                         0.25f, 1.0f / 6.0f };
inline constexpr int   kDivisionDefault = 5;                       // 1/4
inline constexpr float kBpmMin = 40.0f, kBpmMax = 240.0f, kBpmDefault = 120.0f;   // the BPM field's range
inline constexpr float kTimeMinMs = 20.0f, kTimeMaxMs = 2000.0f;   // the free TIME knob's range
inline constexpr float kOffsetMaxMs = 30.0f;                       // the OFFSET knob's reach, ± ms
inline constexpr float kDarkLoHz = 800.0f, kDarkHiHz = 16000.0f;   // DARK's floor and ceiling
inline constexpr float kDarkDefaultPct = 45.0f;                    // ~4 kHz — the tape-ish middle

/** A synced division's time at `bpm`, in ms. The division is clamped onto the ladder, and the tempo to no
    slower than `kBpmMin` (a non-finite tempo counts as `kBpmMin` too): the slowest value is then 1/1 at 40,
    6000 ms, inside the 6.1 s line — so the knob can never name a time the line cannot play. */
inline float divisionMs (int division, double bpm) noexcept
{
    const int    d = division < 0 ? 0 : (division >= kNumDivisions ? kNumDivisions - 1 : division);
    const double b = std::isfinite (bpm) && bpm > (double) kBpmMin ? bpm : (double) kBpmMin;
    return (float) ((double) kDivisionBeats[d] * 60000.0 / b);
}

/** The sync policy: a synced delay plays the division at the tempo the caller chose (the host's, or its own
    BPM field — that choice is the caller's), a free one plays its TIME knob. */
inline float timeMsFor (bool sync, int division, double bpm, float freeMs) noexcept
{
    return sync ? divisionMs (division, bpm) : freeMs;
}

/** DARK is a percent of darkness, so the knob's direction tells the truth: clockwise darkens. The corner falls
    log-evenly from the bright ceiling to the dark floor, every degree darkening by the same ear-step. */
inline float darkHzForPercent (float percent) noexcept
{
    return kDarkHiHz * (float) core::det::pow ((double) (kDarkLoHz / kDarkHiHz), (double) (percent * 0.01f));
}

} // namespace felitronics::delay

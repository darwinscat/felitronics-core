// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

//==============================================================================
// felitronics::blend — mix parameter model (JUCE-free, headless). THE single source of truth + defaults
// for a mixer strip / the Master bus, so the blend engine (IrBlend.h) is a PURE FUNCTION OF DATA, not
// widgets. Promoted VERBATIM from OrbitCapture's app-local model/MixModel.h (de-monolith step 3 →
// promoted step 4) — this is now THE home of the mixer defaults.
//==============================================================================

namespace felitronics::blend
{

// Filter-cutoff drag bounds (were CaptureComponent statics kHpfLo…kLpfHi). "Parked at the extreme" == off.
constexpr double kHpfLo = 20.0, kHpfHi = 2000.0, kLpfLo = 500.0, kLpfHi = 20000.0, kShiftMsMax = 2.0;

struct Filter { bool on = false; double hz = 0.0; int slopeDb = 0; };

// THE one home of the mixer defaults (were hardcoded in ≥3 places; kLpfDef=6500 used to diverge).
struct StripParams {
    double gainDb = 0.0, phaseDeg = 0.0, shiftMs = 0.0;
    bool   solo = false, mute = false;
    Filter hpf { false, 80.0, 24 };
    Filter lpf { false, 8000.0, 12 };
};
struct MasterParams {
    double gainDb = 0.0;
    Filter hpf { false, 80.0, 24 };
    Filter lpf { false, 8000.0, 12 };
};

// The slope a filter is actually running at — 0 when off or parked at its extreme (that side bypassed).
inline int hpActiveSlope (const Filter& f) noexcept { return (f.on && f.hz > kHpfLo + 0.5) ? f.slopeDb : 0; }
inline int lpActiveSlope (const Filter& f) noexcept { return (f.on && f.hz < kLpfHi - 0.5) ? f.slopeDb : 0; }

} // namespace felitronics::blend

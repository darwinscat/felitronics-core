// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <bit>
#include <cstddef>
#include <cstdint>

// felitronics::eq base types (FilterType, Lane, LaneParams, BandParams) + the module's internal ALIAS
// HUB: it pulls the cross-module names the eq engine uses — kMaxChannels + kPi + Smoother
// (felitronics::core) and SpectrumTap (felitronics::analysis) — into felitronics::eq, so the
// MatchedBiquad/Svf/EqBand/EqEngine bodies migrated from TabbyEQ's teq:: sources stay byte-identical
// (they reference these unqualified). The SSOT kMaxChannels lives in felitronics::core (see Config.h).
#include <felitronics/core/Config.h>
#include <felitronics/core/Math.h>
#include <felitronics/core/Smoother.h>
#include <felitronics/analysis/SpectrumTap.h>

namespace felitronics::eq
{

// --- internal alias hub (see the header note) ---
using felitronics::core::kMaxChannels;
using felitronics::core::kPi;
using felitronics::core::Smoother;
using felitronics::analysis::SpectrumTap;

enum class FilterType
{
    Bell,        // peaking, gainDb + Q
    LowShelf,    // gainDb (Butterworth slope, Q ignored)
    HighShelf,   // gainDb (Butterworth slope, Q ignored)
    HighPass,    // Q at 12 dB/oct; 24 dB/oct is Butterworth (Q ignored)
    LowPass,     // Q at 12 dB/oct; 24 dB/oct is Butterworth (Q ignored)
    BandPass,    // unity gain at centre; Q = overall −3 dB bandwidth (order-INVARIANT), skirt
                 // steepness by slope (order = slope/6; slope 6/12 = the frozen single matched
                 // band-pass, bit-for-bit). Near Nyquist CENTRE UNITY WINS: the true peak may read
                 // up to ~+0.29 dB just off f0 — matched::bandpassCascade's documented trade
    Notch,       // band-stop: deep null at f0, unity at DC; width by Q, steepness by slope. Toward
                 // Nyquist it tracks the analog Butterworth band-stop residual (a wide/high notch
                 // legitimately reads below 0 dB at fs/2; slope 6/12 = the frozen single notch,
                 // which under-shoots that residual — see matched::notch)
    AllPass,     // flat magnitude, 360° phase rotation through f0 (Q = sharpness)
    Tilt         // spectral tilt about f0: lows -gainDb, highs +gainDb
};

// Placement lanes (Pro-Q-style, but multi-select — a point may live in several at once). The fixed
// enum order IS the processing order (ST first), which the engine, the analytic 2×2 matrix and the FIR
// builder all compose in — see EqBand::process / matrixResponse. `kNumLanes` is the SSOT count.
enum class Lane : std::uint8_t { Stereo, Left, Right, Mid, Side };
inline constexpr int kNumLanes = 5;

// One placement lane's full parameter set. The point's `type`/`swept` are SHARED across its lanes
// (decision #2); everything that can differ per domain lives here.
struct LaneParams
{
    bool   on     = false;    // lane enabled (the menu checkbox)
    double freq   = 1000.0;   // Hz (engine clamps to [10, 0.49*fs])
    double Q      = 1.0;
    double gainDb = 0.0;      // bells & shelves
    int    slope  = 12;       // HP/LP: 6..96 dB/oct Butterworth. Notch/BandPass: skirt steepness
                              // (order = slope/6; Q stays the −3 dB bandwidth, order-invariant)
    bool   bypass = false;    // lane kept but muted (ghost node) — distinct from on=false

    // Doubles compared by bit pattern (not `==`) so the engine's recompute-skip stays exact without
    // tripping -Wfloat-equal in strict-warning builds.
    bool operator== (const LaneParams& o) const noexcept
    {
        auto bits = [] (double d) noexcept { return std::bit_cast<std::uint64_t> (d); };
        return on == o.on && slope == o.slope && bypass == o.bypass
            && bits (freq) == bits (o.freq) && bits (Q) == bits (o.Q) && bits (gainDb) == bits (o.gainDb);
    }
};

// One band ("point"): a filter of one shared `type`, split across a set of placement lanes. The plugin
// adapter maps APVTS params into this; the engine owns no parameter system of its own (framework-
// agnostic). A fresh band reproduces the pre-lanes default exactly: point off, one Stereo lane enabled
// at 1 kHz / Q 1 / 0 dB / slope 12, all other lanes off.
struct BandParams
{
    bool       on     = false;               // the point exists
    FilterType type   = FilterType::Bell;    // SHARED by all lanes (decision #2)
    bool       swept  = false;               // search band; honored only in the single-ST configuration
    bool       bypass = false;               // whole-point bypass (strip power button)
    LaneParams lanes[kNumLanes] { { true }, {}, {}, {}, {} };   // lanes[Stereo].on = true, the rest off

    // Index a lane by its enum (the array is public, but this reads without the static_cast noise).
    LaneParams&       lane (Lane l)       noexcept { return lanes[static_cast<std::size_t> (l)]; }
    const LaneParams& lane (Lane l) const noexcept { return lanes[static_cast<std::size_t> (l)]; }

    // Exact change-detection across the shared fields + every lane (each LaneParams compares its own
    // doubles by bit pattern), so the engine's recompute-skip stays exact.
    bool operator== (const BandParams& o) const noexcept
    {
        if (! (on == o.on && type == o.type && swept == o.swept && bypass == o.bypass)) return false;
        for (int i = 0; i < kNumLanes; ++i) if (! (lanes[i] == o.lanes[i])) return false;
        return true;
    }
};

} // namespace felitronics::eq

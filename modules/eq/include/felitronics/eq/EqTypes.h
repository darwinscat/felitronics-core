// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <bit>
#include <cstdint>

// felitronics::eq base types (FilterType, BandParams) + the module's internal ALIAS HUB: it pulls the
// cross-module names the eq engine uses — kMaxChannels + kPi + Smoother (felitronics::core) and
// SpectrumTap (felitronics::analysis) — into felitronics::eq, so the MatchedBiquad/Svf/EqBand/EqEngine
// bodies migrated from TabbyEQ's teq:: sources stay byte-identical (they reference these unqualified).
// The SSOT kMaxChannels lives in felitronics::core (see Config.h).
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
    BandPass,    // Q, unity gain at centre
    Notch,       // band-stop: deep null at f0, unity at DC; width by Q, steepness by slope. Toward
                 // Nyquist it tracks the analog Butterworth band-stop residual (a wide/high notch
                 // legitimately reads below 0 dB at fs/2; slope 6/12 = the frozen single notch,
                 // which under-shoots that residual — see matched::notch)
    AllPass,     // flat magnitude, 360° phase rotation through f0 (Q = sharpness)
    Tilt         // spectral tilt about f0: lows -gainDb, highs +gainDb
};

// One band's full specification. The plugin adapter maps APVTS params into this; the engine owns no
// parameter system of its own (stays framework-agnostic).
struct BandParams
{
    bool       on     = false;
    FilterType type   = FilterType::Bell;
    double     freq   = 1000.0;   // Hz (engine clamps to [10, 0.49*fs])
    double     Q      = 1.0;
    double     gainDb = 0.0;      // bells & shelves
    int        slope  = 12;       // HP/LP: 6..96 dB/oct Butterworth. Notch: steepness (order=slope/6),
                                  // sections=ceil(order/2); 6/12→single notch (Q=width), 24→2 … 96→8.
    bool       swept  = false;    // true → zero-delay SVF (smooth fast fc sweeps for search mode).
                                  // Notch sweeps as a real SVF notch; Tilt has no one-SVF realisation
                                  // and always runs the matched two-shelf design (flag ignored).
    bool       bypass = false;    // band kept but muted (ghost) — distinct from on=false (removed)

    // M/S dual-mode: the flat fields above are the Mid/main lane; these are the independent Side lane.
    // Only used when ms==true on a 2-channel signal (mono/surround ignore it -> Mid lane only).
    bool       ms      = false;          // false = Stereo (one lane on L/R), true = Mid/Side (two lanes)
    bool       sOn     = true;           // Side lane enabled
    FilterType sType   = FilterType::Bell;
    double     sFreq   = 1000.0;
    double     sQ      = 1.0;
    double     sGainDb = 0.0;
    int        sSlope  = 12;
    bool       sBypass = false;

    // Exact change-detection. The doubles are compared by bit pattern (not `==`) so the engine's
    // recompute-skip stays exact without tripping -Wfloat-equal in strict-warning builds.
    bool operator== (const BandParams& o) const noexcept
    {
        auto bits = [] (double d) noexcept { return std::bit_cast<std::uint64_t> (d); };
        return on == o.on && type == o.type && slope == o.slope && swept == o.swept && bypass == o.bypass
            && ms == o.ms && sOn == o.sOn && sType == o.sType && sSlope == o.sSlope && sBypass == o.sBypass
            && bits (freq) == bits (o.freq) && bits (Q) == bits (o.Q) && bits (gainDb) == bits (o.gainDb)
            && bits (sFreq) == bits (o.sFreq) && bits (sQ) == bits (o.sQ) && bits (sGainDb) == bits (o.sGainDb);
    }
};

} // namespace felitronics::eq

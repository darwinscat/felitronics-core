// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <cmath>
#include <cstdint>

//==============================================================================
// felitronics::measurement::levelprobe — the ONE deterministic reference stimulus that defines what
// "unity" means across the Felitronics family: LCG white noise through a one-pole low-pass at
// kShapeHz, level-calibrated analytically to kRefRmsDb RMS.
//
// THE CONTRACT (this header + ReferenceUnity.h — read both): a stage or an IR is "at unity" when
// THIS stimulus passes through it with 0 dB RMS gain. The two sides of the contract:
//   • the PROBE (here) — feed the stimulus through a stage at a fixed operating point and read the
//     RMS ratio ("reference gain"). Plugins run it OFF-THREAD: e.g. OrbitCab measures a just-loaded
//     capture model before the atomic swap so the tube stage can be level-matched deterministically
//     (the baked tube reference constant was measured with the SAME stimulus);
//   • the NORMALIZER (ReferenceUnity.h) — the load-time IR gain that makes the IR convolve this
//     stimulus at unity RMS, computed in the frequency domain. Its spectral weight
//     w(f) = 1 / (1 + (f/kShapeHz)^2) IS this stimulus's power spectrum — the constants live HERE,
//     the normalizer includes this header, so the two sides cannot drift apart.
//
// LP-shaped (not white) because guitar into an amp carries little energy above a few kHz — a white
// probe would overweight HF the stages never see and skew the measured gain. −18 dBFS RMS is the
// seam's typical playing level (the level capture-model stages normalize toward). Integer-
// deterministic (LCG on the sample index): identical on every platform / run / rate.
//
// FLOAT ON PURPOSE (the measurement-module double convention does not apply here — documented
// exception): the stimulus fills the float buffers the audio stages consume, and the shipped
// numbers are GOLDEN — the family's unity anchor must stay bit-identical to what products already
// measured against (OrbitCab's tube constant, its factory-IR loudness). Pure std, header-only,
// no FFT — safe for engine cores and offline tools alike.
//==============================================================================
namespace felitronics::measurement::levelprobe
{

constexpr double kRefRmsDb   = -18.0;   // stimulus RMS (dBFS) — the family reference level
constexpr double kShapeHz    = 2000.0;  // one-pole LP corner — "guitar-band" shaping; ReferenceUnity weights with THIS
constexpr double kSettleSec  = 0.10;    // discard: stage + shaping-filter settle
constexpr double kMeasureSec = 0.20;    // measured window after the settle

// n-indexed white noise in [-1, 1] — same LCG family as the test benches. Integer + exact-double
// arithmetic only (no libm): bit-identical on every platform, so its values are pinned as goldens.
inline float white (std::int64_t n) noexcept
{
    std::uint32_t s = (std::uint32_t) (n * 2654435761u + 1013904223u);
    s ^= s >> 15; s *= 2246822519u; s ^= s >> 13;
    return ((float) ((double) s / 4294967296.0) - 0.5f) * 2.0f;
}

// Fill `dst` with the shaped, level-calibrated stimulus. The LP state starts at 0 and the gain
// normalizes the SHAPED noise back to kRefRmsDb: for a one-pole LP with coefficient a on white
// noise of RMS r, the output RMS is r·sqrt(a / (2 − a)) — undo that analytically so the level
// is exact without a measurement pass.
inline void fill (float* dst, int numSamples, double sampleRate) noexcept
{
    const double a = 1.0 - std::exp (-6.283185307179586 * kShapeHz / (sampleRate > 0.0 ? sampleRate : 48000.0));
    const double whiteRms = 1.0 / std::sqrt (3.0);                        // uniform [-1,1]
    const double shapedRms = whiteRms * std::sqrt (a / (2.0 - a));        // one-pole LP on white
    const float g = (float) (std::pow (10.0, kRefRmsDb / 20.0) / shapedRms);
    float lp = 0.0f;
    for (int n = 0; n < numSamples; ++n)
    {
        lp += (float) a * (white (n) - lp);
        dst[n] = lp * g;
    }
}

} // namespace felitronics::measurement::levelprobe

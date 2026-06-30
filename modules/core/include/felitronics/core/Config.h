// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <cstddef>

//==============================================================================
// felitronics::core — the shared base every module derives from: the sample type and the size
// configuration (Law 3 + Law 5). JUCE-free, zero deps.
//==============================================================================
namespace felitronics::core
{

// Law 3: float in the hot path; `double` only for offline coefficient design and meter / LUFS /
// true-peak accumulators. `Sample` keeps raw `float*` out of public signatures so full sample-type
// templating later (the `embedded-fpu` tier, still float) is a flag-flip, not a fork-rewrite.
// NB: `bare-mcu` fixed-point is a SEPARATE codebase (different topology / scaling / quantization) —
// this alias does NOT turn a float kernel into a fixed-point one.
using Sample = float;

// Law 5: configurable sizes. These are the DESKTOP-tier defaults and the single source of truth that
// every fixed-size, zero-allocation per-channel state array derives from (so `process()` stays
// RT-safe). A per-tier build overrides them via a CMake preset / policy later; for now this is the
// one place to change.
//
// kMaxChannels = 16 covers every commercial layout — mono, stereo, quad, 5.1, 7.1, 7.1.4 / 9.1.6
// Atmos (<= 16) and 1st–3rd-order ambisonics. The only cost of the headroom is memory.
constexpr int kMaxChannels = 16;

// A conservative upper bound for block-scoped scratch on the funded tiers. Hosts that exceed it must
// be clamped by the adapter (the engines clamp defensively too).
constexpr int kMaxBlockSize = 8192;

} // namespace felitronics::core

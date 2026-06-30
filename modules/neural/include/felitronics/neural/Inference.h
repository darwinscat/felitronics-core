// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/core/Config.h>

#include <concepts>

namespace felitronics::neural
{

//==============================================================================
// The PROCESS-ONLY inference seam (ADR §3/§8). A neural backend (NAM/Eigen, RTNeural, a tiny MCU net…)
// satisfies this compile-time concept; the core never sees the backend's headers, model format, file
// I/O, or exceptions — those live in the adapter. The seam is multichannel + in place (matches the rest
// of felitronics-core and orbitcab's cab::AmpStage). Reached as a TEMPLATE (no vtable in the hot path).
//
//   prepare(sampleRate, maxBlock, maxChannels)  — off the audio thread; may allocate / prewarm.
//   process(io, numChannels, numSamples)         — RT-safe: no alloc/lock/IO/throw; in place.
//   reset()                                      — RT-safe state clear.
//   latencySamples()                             — host-rate latency (incl. any backend resampling).
//
// NOT in the interface: model loading / parsing / paths. The adapter builds a prepared instance and
// hands it to NeuralStage (below) for a swap-safe handoff.
template <class T>
concept Inference =
    requires (T t, const T ct, double sr, int maxBlock, int maxChannels, float* const* io, int nc, int n)
{
    { t.prepare (sr, maxBlock, maxChannels) } noexcept -> std::same_as<void>;
    { t.process (io, nc, n) }                noexcept -> std::same_as<void>;
    { t.reset() }                            noexcept -> std::same_as<void>;
    { ct.latencySamples() }                  noexcept -> std::same_as<int>;
};

} // namespace felitronics::neural

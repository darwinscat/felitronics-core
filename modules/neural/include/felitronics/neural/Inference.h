// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/core/Config.h>

#include <concepts>

namespace felitronics::neural
{

//==============================================================================
// The PROCESS-ONLY inference seam (ADR §3/§8). A neural backend (a NAM runner, RTNeural, a tiny MCU net…)
// satisfies this compile-time concept; the core never sees the backend's headers, model format, file
// I/O, or exceptions — those live in the adapter. The seam is multichannel + in place, like the rest of
// felitronics-core. Reached as a TEMPLATE (no vtable in the hot path).
//
//   prepare(sampleRate, maxBlock, maxChannels)  — off the audio thread; may allocate / prewarm.
//   process(io, numChannels, numSamples) -> bool — RT-safe: no alloc/lock/IO/throw; in place. Law 11:
//                                                 false = the call was REFUSED and nothing changed.
//   reset()                                      — the STREAM RESTART: audio thread, and what it clears
//                                                 is the audio the caller fed, so the next stream starts
//                                                 where a freshly prepared instance would — as far as
//                                                 the backend's own memory is finite and reachable. A
//                                                 backend that cannot reach all of it (a recurrent cell,
//                                                 a third party's clock) says so with a number rather
//                                                 than promising, in its own header. RT-safe in the sense law 2 means — no alloc,
//                                                 lock, IO or throw beyond whatever process() already
//                                                 costs — but NOT necessarily O(the block): a backend
//                                                 with a receptive field has to spend it, and it is the
//                                                 backend's job to publish what that costs. The NAM one
//                                                 does: 3.77 ms per lane for a real WaveNet at a
//                                                 64-sample block, 282 % of that callback.
//   latencySamples()                             — host-rate latency (incl. any backend resampling).
//
// NOT in the interface: model loading / parsing / paths. The adapter builds a prepared instance and
// hands it to NeuralStage (below) for a swap-safe handoff.
template <class T>
concept Inference =
    requires (T t, const T ct, double sr, int maxBlock, int maxChannels, float* const* io, int nc, int n)
{
    { t.prepare (sr, maxBlock, maxChannels) } noexcept -> std::same_as<void>;
    { t.process (io, nc, n) }                noexcept -> std::same_as<bool>;   // law 11: the verdict is RETURNED
    { t.reset() }                            noexcept -> std::same_as<void>;
    { ct.latencySamples() }                  noexcept -> std::same_as<int>;
};

} // namespace felitronics::neural

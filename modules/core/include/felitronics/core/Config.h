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

// THE LOWEST SAMPLE RATE THE CORE TAKES AS AUDIO (P51). NOT a size and NOT a tier knob like the two above: a
// validity floor, and one number for every entry that takes a rate a caller could have got wrong — the loudness
// search (mastering::TargetLoudnessSolver), the mastering chain and the delivery resampler under the C ABI's
// `fc_master_create`, and the probe ABI's measurers (fcore::Probe, fcore::ShapeProbe and the analysis classes
// behind it). Each keeps its own name for it (`X::kMinSampleRate`), and every one of them is THIS constant, because
// the floor used to be three answers to one question: any rate > 0 (the search), whatever the chain's stages
// happened to refuse (above 50 Hz with the limiter on, 20.4 Hz with only the EQ, nothing without either), and
// 1000 Hz (the probe and its analyzers, in seven copies). The CEILINGS stay per class, with their own reasons
// (768 kHz for the offline measurers, 3 MHz for the real-time stages).
//
// Why 8000, and the value is the pair of these, not either alone:
//   * the BS.1770 K-weighting shelf is designed at 1681.97 Hz (analysis::KWeightingFilter::kShelfHz), so below
//     twice that, 3363.95 Hz, it is past Nyquist: the bilinear tan() has wrapped, the shelf is aliased everywhere
//     down there, and in the bands where tan() is negative — (1682, 3364) Hz, (841, 1121) Hz, … — the filter is
//     UNSTABLE. Measured on a 0 dBFS 400 Hz sine: -3.72 LUFS at 48 kHz, -3.16 at 3364 Hz, +3043 LUFS at 3300 Hz.
//     3364 is an edge, not a floor — the shelf's pole radius is 0.99997 there and 0.43 at 8000. (Stable is not exact: 8000 Hz still reads that sine 0.21 LU
//     away from 48 kHz, because the bilinear design warps.)
//   * 8000 Hz is the lowest standard audio rate, so no real programme sits below it, and what does arrive below
//     it is a mistake to refuse rather than measure: a rate passed in KILOHERTZ (44.1, 48, 88.2, 96, 192 — every
//     one below 3364) and a corrupt file header. Before P51 a search handed 88.2 reported Solved at -14 LUFS,
//     without a single error, over a chain running a thousand times too slow.
// Compared as `sampleRate >= kMinSampleRate`: 8000 itself is audio, and a NaN fails the comparison.
constexpr double kMinSampleRate = 8000.0;

} // namespace felitronics::core

// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

//==============================================================================
// felitronics::measurement — Farina deconvolution. Convolve the recording with the sweep's inverse
// filter: the linear IR lands at (sweepLen-1 + τ) for an unknown round-trip τ; harmonic images sit
// at EARLIER indices. So we search [sweepLen-1, end) for the linear IR's onset — this BOTH excludes
// the harmonics AND absorbs latency (a fixed lag = sweepLen-1 silently returns floor noise when τ
// pushes the IR out of a snug window). OFFLINE, MESSAGE-THREAD-ONLY (double, allocates).
//==============================================================================

#include <felitronics/measurement/Convolve.h>
#include <felitronics/measurement/IrPost.h>
#include <felitronics/measurement/Sweep.h>

#include <cstddef>
#include <span>
#include <vector>

namespace felitronics::measurement
{

struct DeconvResult
{
    std::vector<double> full;             // recording ⊛ inverse
    std::size_t         linearLag = 0;    // index of the linear IR's onset (latency-corrected)
    std::size_t         latencySamples = 0; // measured round-trip τ, samples
    double              sampleRate = 0.0;
};

// Deconvolve a recording against its sweep. `onsetThresholdDb`/`onsetPreRoll` tune the onset search
// (promoted from hardcoded to defaulted args — calibration without a fork).
inline DeconvResult deconvolve (std::span<const double> recording, const Sweep& sw,
                                double onsetThresholdDb = -40.0, std::size_t onsetPreRoll = 8)
{
    DeconvResult r;
    r.full       = convolve (recording, sw.inverse);
    r.sampleRate = sw.spec.sampleRate;
    const std::size_t searchStart = (sw.sweepLen > 0) ? sw.sweepLen - 1 : 0;
    const OnsetResult on = detectOnset (r.full, onsetThresholdDb, onsetPreRoll, searchStart);
    r.linearLag      = on.onset;
    r.latencySamples = (on.onset >= searchStart) ? on.onset - searchStart : 0;
    return r;
}

// Copy the linear IR window [linearLag, linearLag+len) out of the deconvolution (zero-padded past end).
inline std::vector<double> extractIr (const DeconvResult& d, std::size_t len)
{
    std::vector<double> ir (len, 0.0);
    for (std::size_t k = 0; k < len; ++k)
    {
        const std::size_t idx = d.linearLag + k;
        if (idx < d.full.size()) ir[k] = d.full[idx];
    }
    return ir;
}

} // namespace felitronics::measurement

// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

//==============================================================================
// felitronics::measurement — fine time/polarity alignment of IRs by normalized cross-correlation.
// OFFLINE, MESSAGE-THREAD-ONLY (double, allocates).
//
// Promoted from OrbitCapture's mixer "Auto" button (core-reuse audit N1): given a reference IR and
// a channel IR that SHOULD be coherent with it (two mics on one cabinet, two IRs of one rig), find
// the lag within ±maxLag that maximizes |normalized cross-correlation| over a window straddling
// the onsets, and the polarity (a negative correlation peak = the channel is phase-flipped).
// Consumers map the result onto their own controls: OrbitCapture sets a strip's time-shift knob
// (positive shiftSamples = DELAY the channel) + a 180° phase; OrbitCab's auto-polarity uses only
// `invert` (its slots have no time-shift).
//
// The onset window uses IrPost::detectOnset (peak-backwards, floor-immune) — NOT a naive
// first-sample-over-threshold scan.
//==============================================================================

#include "IrPost.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <span>

namespace felitronics::measurement
{

struct XcorrAlignment
{
    double shiftSamples = 0.0;   // apply to the CHANNEL: positive = delay it (it led the reference)
    bool   invert = false;       // polarity flip (the correlation peak was negative)
    double corr = 0.0;           // |peak| of the normalized cross-correlation, 0..1 (confidence)
};

// Align `ch` against `ref`: best lag within ±maxLagSamples by normalized cross-correlation over a
// `window`-sample segment starting just before the earlier onset. corr stays 0 when there is
// nothing to correlate (silence, degenerate sizes) — treat that as "no confident suggestion".
inline XcorrAlignment xcorrAlign (std::span<const double> ref, std::span<const double> ch,
                                  int maxLagSamples, std::size_t window = 8192) noexcept
{
    XcorrAlignment out;
    if (ref.empty() || ch.empty() || maxLagSamples <= 0) return out;

    const std::size_t onR = detectOnset (ref).onset;
    const std::size_t onC = detectOnset (ch).onset;
    const std::size_t onMin = std::min (onR, onC);
    const std::size_t start = onMin > (std::size_t) maxLagSamples ? onMin - (std::size_t) maxLagSamples : 0;

    const std::size_t common = std::min (ref.size(), ch.size());
    if (common <= start + (std::size_t) maxLagSamples + 1) return out;
    const std::size_t win = std::min (window, common - start - (std::size_t) maxLagSamples - 1);
    if (win <= 16) return out;

    double eR = 0.0, eC = 0.0;
    for (std::size_t i = 0; i < win; ++i)
    {
        eR += ref[start + i] * ref[start + i];
        eC += ch[start + i]  * ch[start + i];
    }
    if (eR <= 0.0 || eC <= 0.0) return out;

    double best = 0.0; int bestLag = 0;
    for (int lag = -maxLagSamples; lag <= maxLagSamples; ++lag)
    {
        double s = 0.0;
        for (std::size_t i = 0; i < win; ++i)
        {
            const long long j = (long long) (start + i) + lag;
            if (j >= 0 && j < (long long) ch.size()) s += ref[start + i] * ch[(std::size_t) j];
        }
        if (std::fabs (s) > std::fabs (best)) { best = s; bestLag = lag; }
    }

    out.corr = std::fabs (best) / std::sqrt (eR * eC);
    out.invert = best < 0.0;
    // Peak at lag=+L means ch[i+L] ≈ ref[i]: the channel arrives L samples LATE → advance it
    // (negative shift). Peak at lag=-L: the channel LEADS by L → delay it (positive shift).
    out.shiftSamples = -(double) bestLag;
    return out;
}

// Whole set against irs[ref]; entry [ref] stays identity (it IS the time reference).
inline std::vector<XcorrAlignment> xcorrAlignSet (std::span<const std::vector<double>> irs,
                                                  std::size_t ref, int maxLagSamples,
                                                  std::size_t window = 8192)
{
    std::vector<XcorrAlignment> out (irs.size());
    if (ref >= irs.size()) return out;
    for (std::size_t m = 0; m < irs.size(); ++m)
        if (m != ref)
            out[m] = xcorrAlign (irs[ref], irs[m], maxLagSamples, window);
    return out;
}

} // namespace felitronics::measurement

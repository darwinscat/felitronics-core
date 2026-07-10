// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

//==============================================================================
// felitronics::measurement — post-processing of a derived IR: onset detection, onset-trim, peak gain.
// OFFLINE, MESSAGE-THREAD-ONLY (double, allocates). Non-destructive by intent — the caller decides
// whether/where to apply a gain (the raw recording keeps its true level).
//==============================================================================

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <span>
#include <vector>

namespace felitronics::measurement
{

struct OnsetResult
{
    std::size_t onset     = 0;   // trim point (attack preserved by preRoll)
    std::size_t peakIndex = 0;
    double      peak      = 0.0;
};

// Onset detection robust to leading floor / latency-gap / ring: find the PEAK in [searchStart, end),
// then walk BACK to where the signal drops below a peak-relative threshold (the attack), then back off
// preRoll. Searching back from the peak (not forward from the start) makes it immune to any pre-onset
// floor above the threshold (deconv skirt, harmonic-tail residue, mic self-noise in the latency gap).
inline OnsetResult detectOnset (std::span<const double> ir,
                                double thresholdDb = -40.0,
                                std::size_t preRoll = 8,
                                std::size_t searchStart = 0) noexcept
{
    OnsetResult r;
    if (searchStart >= ir.size()) { r.onset = ir.empty() ? 0 : ir.size() - 1; return r; }
    for (std::size_t i = searchStart; i < ir.size(); ++i)
    {
        const double a = std::fabs (ir[i]);
        if (a > r.peak) { r.peak = a; r.peakIndex = i; }
    }
    if (r.peak <= 0.0) { r.onset = searchStart; return r; }
    const double thr = r.peak * std::pow (10.0, thresholdDb / 20.0);
    std::size_t on = r.peakIndex;
    while (on > searchStart && std::fabs (ir[on - 1]) >= thr) --on;
    r.onset = (on > searchStart + preRoll) ? on - preRoll : searchStart;
    return r;
}

// Copy the IR from the onset forward. len == 0 → to end; else exactly len (zero-padded past the source).
inline std::vector<double> trimToOnset (std::span<const double> ir, const OnsetResult& on, std::size_t len = 0)
{
    const std::size_t start = on.onset;
    const std::size_t n = (len > 0) ? len : (ir.size() > start ? ir.size() - start : 0);
    std::vector<double> out (n, 0.0);
    for (std::size_t k = 0; k < n; ++k)
    {
        const std::size_t idx = start + k;
        if (idx < ir.size()) out[k] = ir[idx];
    }
    return out;
}

// Gain that brings the signal's peak to targetPeak (1.0 if silent). The caller applies it where it wants.
inline double peakGain (std::span<const double> x, double targetPeak = 0.98) noexcept
{
    double pk = 0.0;
    for (double v : x) pk = std::max (pk, std::fabs (v));
    return (pk > 0.0) ? targetPeak / pk : 1.0;
}

inline void applyGain (std::span<double> x, double g) noexcept { for (double& v : x) v *= g; }

} // namespace felitronics::measurement

// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

//==============================================================================
// felitronics::measurement — multi-mic IR set alignment. One sweep, N mics: cut EVERY mic's IR at the
// SAME absolute index (the earliest onset across the set) so the inter-mic delays — the sound of a mic
// set — survive. A per-mic re-align would destroy the comb structure that IS the reason to use two mics.
// OFFLINE, MESSAGE-THREAD-ONLY (double, allocates).
//==============================================================================

#include <felitronics/measurement/Deconvolve.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <vector>

namespace felitronics::measurement
{

struct AlignedIrSet
{
    std::vector<std::vector<double>> irs;          // one per mic, ALL cut at the same absolute index
    std::vector<std::size_t>         delaySamples; // linearLag[m] - refLag — the inter-mic delays (preserved)
    std::size_t                      refLag = 0;   // the ONE common time reference = earliest onset
};

// Cut a deconvolved mic set to `irLen` samples from the common reference lag (= the earliest mic onset).
inline AlignedIrSet alignToCommonOnset (std::span<const DeconvResult> mics, std::size_t irLen)
{
    AlignedIrSet out;
    if (mics.empty()) return out;

    std::size_t refLag = std::numeric_limits<std::size_t>::max();
    for (const auto& m : mics) refLag = std::min (refLag, m.linearLag);
    out.refLag = refLag;

    out.irs.reserve (mics.size());
    out.delaySamples.reserve (mics.size());
    for (const auto& m : mics)
    {
        std::vector<double> ir (irLen, 0.0);
        for (std::size_t k = 0; k < irLen; ++k)
        {
            const std::size_t idx = refLag + k;
            if (idx < m.full.size()) ir[k] = m.full[idx];
        }
        out.irs.push_back (std::move (ir));
        out.delaySamples.push_back (m.linearLag - refLag);
    }
    return out;
}

// ONE common gain across the whole set (relative balance between mics preserved), bringing the set's
// peak to targetPeak. Apply per-IR with applyGain.
inline double peakGain (std::span<const std::vector<double>> set, double targetPeak = 0.98) noexcept
{
    double pk = 0.0;
    for (const auto& ir : set)
        for (double v : ir) if (std::isfinite (v)) pk = std::max (pk, std::fabs (v));   // ignore poisoned samples
    return (pk > 0.0) ? targetPeak / pk : 1.0;
}

} // namespace felitronics::measurement

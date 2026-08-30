// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.
#pragma once
// felitronics::rigplayer — HOW FAR ONE MODEL LAGS ANOTHER. Two NAM models of the same device, run on the same
// probe, come out shifted by a few samples (different trained latencies); the player must delay the
// faster one so a crossfade between them does not comb-filter. The measurement is a normalized
// cross-correlation over a small search window — pure math, extracted here from the player's alignment
// pass (the NAM inference and caching stay in the orchestrator). std-only, ctest (ModelAlignmentTest).

#include <cmath>
#include <cstddef>
#include <span>

namespace felitronics::rigplayer {

inline constexpr int kSearch = 96;   // samples either way; real drift between two models is single digits

struct Best { int lag = 0; double corr = -1.0; };

// The lag (in [-search, search]) at which `cur` best aligns to `ref` by NORMALIZED cross-correlation,
// evaluated over the interior [search, N-search) so every shifted index is in range. The result is the
// DELAY of `cur` relative to `ref` (cur[i] approx ref[i-lag]); read `cur` at index k+lag to align it.
// Returns {0, -1} when the two differ in length, are too short to search, or carry no energy at all —
// "no reading", which is a different answer from "aligned at zero".
inline Best bestLag (std::span<const float> ref, std::span<const float> cur, int search = kSearch) {
    Best out;
    if (ref.size() != cur.size() || ref.size() <= (std::size_t) (2 * search)) return out;
    const std::size_t n = ref.size();
    for (int lag = -search; lag <= search; ++lag) {
        double num = 0.0, xx = 0.0, yy = 0.0;
        for (std::size_t k = (std::size_t) search; k + (std::size_t) search < n; ++k) {
            const double u = ref[k], v = cur[(std::size_t) ((long) k + lag)];
            num += u * v; xx += u * u; yy += v * v;
        }
        // A lag with no energy on either side is not a WORSE alignment, it is no reading at all —
        // scoring it 0 made the first lag tried beat the {0, -1} that means "nothing to say", so two
        // silent models came back aligned at -96 samples and the player delayed one against the other.
        if (xx <= 0.0 || yy <= 0.0) continue;
        const double c = num / std::sqrt (xx * yy);
        if (c > out.corr) { out.corr = c; out.lag = lag; }
    }
    return out;
}

} // namespace felitronics::rigplayer
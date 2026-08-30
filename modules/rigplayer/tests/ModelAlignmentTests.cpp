// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.
#include <felitronics/rigplayer/ModelAlignment.h>
#include <felitronics_test.h>
#include <cmath>
#include <vector>
using namespace felitronics::rigplayer;
namespace t = felitronics::test;

// a deterministic "signal": a sum of a few sinusoids, not periodic within the window
static std::vector<float> sig (std::size_t n) {
    std::vector<float> s (n);
    for (std::size_t i = 0; i < n; ++i)
        s[i] = (float) (std::sin (i * 0.017) + 0.6 * std::sin (i * 0.041 + 1.0) + 0.3 * std::sin (i * 0.089 + 2.0));
    return s;
}
// cur[i] = ref[i - d]: cur is ref DELAYED by d; bestLag returns +d
static std::vector<float> delayed (const std::vector<float>& r, int d) {
    std::vector<float> c (r.size(), 0.0f);
    for (std::size_t i = 0; i < r.size(); ++i) { const long j = (long) i - d; if (j >= 0 && j < (long) r.size()) c[i] = r[(std::size_t) j]; }
    return c;
}

int main() {
    const auto ref = sig (4096);
    t::group ("best lag recovers a known shift");
    t::ok (bestLag (ref, ref).lag == 0, "a signal aligns to itself at lag 0");
    t::ok (bestLag (ref, ref).corr > 0.999, "…with correlation ~1");
    for (int d : { 1, 5, 12, -7, -30, 48 }) {
        const auto cur = delayed (ref, d);
        const auto b = bestLag (ref, cur);
        t::ok (b.lag == d, "delay " + std::to_string (d) + " -> lag " + std::to_string (b.lag));
    }
    t::group ("degenerate inputs");
    t::ok (bestLag (ref, sig (2048)).corr < 0.0, "different lengths -> {0, -1}");
    std::vector<float> zero (4096, 0.0f);
    // NOT `corr == 0`: a no-energy lag is not a worse alignment, it is no reading at all. Scoring it
    // zero let the first lag tried beat the {0, -1} sentinel, so two silent models came back "aligned"
    // at -96 samples — a delay the player would then put between them.
    { const auto b = bestLag (ref, zero);
      t::ok (b.corr < 0.0 && b.lag == 0, "silence -> no reading, at lag 0"); }
    { const auto b = bestLag (zero, zero);
      t::ok (b.corr < 0.0 && b.lag == 0, "…and two silences agree about nothing, not about -96"); }
    std::vector<float> tiny (32, 1.0f);
    t::ok (bestLag (tiny, tiny).corr < 0.0, "shorter than the search window -> {0, -1}");
    t::group ("beyond the search window is not found");
    { const auto cur = delayed (ref, 200); const auto b = bestLag (ref, cur, 96);
      t::ok (std::abs (b.lag) <= 96, "a shift past the window clamps inside it (not the true 200)"); }
    return t::report();
}

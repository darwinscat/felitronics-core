// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// felitronics::analysis — PERFORMANCE of the spectrum panes, on the scalar reference FFT: what one UI
// tick costs (ingest + 900 log-frequency columns) for the classic pane at every order and for the
// multi-resolution pane at a normal hop and after a 100 ms stall. The numbers are printed for the record
// (a machine-dependent table, compare across commits on one machine); the CHECKS are loose ceilings that
// only a regression of an order of magnitude would trip — a tick must stay far below the 33 ms frame,
// and the multi-res pane must cost less than the sum of its tiers' classic panes plus the stitching.
// The pffft counterpart (same panes on the SIMD backend, and the speed-up) lives in the pffft suite.

#include <felitronics_test.h>
#include <felitronics/analysis/SpectrumPane.h>
#include <felitronics/analysis/MultiResSpectrumPane.h>

#include <chrono>
#include <cstdint>
#include <cstdio>
#include <memory>
#include <string>

using namespace felitronics;
using felitronics::test::ok;
using felitronics::test::group;

namespace
{
std::uint64_t g_seed = 0x9E3779B97F4A7C15ull;
float uni() noexcept { g_seed ^= g_seed >> 12; g_seed ^= g_seed << 25; g_seed ^= g_seed >> 27; return (float) ((g_seed * 0x2545F4914F6CDD1Dull) >> 40) / 8388608.0f - 1.0f; }

// One UI tick: fill the frame, ingest, build 900 columns. Returns microseconds per tick (median of 3 runs).
template <class Pane>
double tickMicros (Pane& p, int frameSamples, int order, int ticks)
{
    analysis::PlotMap pm; pm.width = 900.0f; pm.height = 300.0f; pm.plotBottom = 300.0f;
    volatile float sink = 0.0f;
    double best[3];
    for (int run = 0; run < 3; ++run)
    {
        const auto t0 = std::chrono::steady_clock::now();
        for (int t = 0; t < ticks; ++t)
        {
            float* f = p.frameInput(); for (int i = 0; i < frameSamples; ++i) f[i] = uni();
            p.ingest (order);
            p.buildColumns (pm, 48000.0, 1.5, 1000.0, [&] (int, float, float y, float yp) { sink = sink + y + yp; });
        }
        best[run] = std::chrono::duration<double, std::micro> (std::chrono::steady_clock::now() - t0).count() / ticks;
    }
    if (best[0] > best[1]) std::swap (best[0], best[1]);
    if (best[1] > best[2]) std::swap (best[1], best[2]);
    if (best[0] > best[1]) std::swap (best[0], best[1]);
    return best[1];
}
}

int main()
{
    std::printf ("felitronics::analysis spectrum-pane performance (scalar FFT; one tick = ingest + 900 columns)\n");
    constexpr double kTickBudgetUs = 33333.0;   // one 30 fps frame
    const int ticks = 120;

    group ("classic pane, every order");
    double classicUs[5] = { 0, 0, 0, 0, 0 };
    for (int o = 10; o <= 14; ++o)
    {
        auto p = std::make_unique<analysis::SpectrumPane>();
        const double us = tickMicros (*p, 1 << o, o, ticks);
        classicUs[o - 10] = us;
        std::printf ("    classic %-6d %8.0f us/tick  (%5.2f %% of a 30 fps frame)\n", 1 << o, us, 100.0 * us / kTickBudgetUs);
        ok (us < 0.15 * kTickBudgetUs, "classic " + std::to_string (1 << o) + ": a tick stays under 15 % of the frame (" + std::to_string ((int) us) + " us)");
    }
    ok (classicUs[4] > classicUs[0], "16384 costs more than 1024 (the timer measures something)");

    group ("multi-resolution pane, 16384 / 4096 / 1024");
    for (int cover : { 1600, 4800 })
    {
        auto m = std::make_unique<analysis::MultiResSpectrumPane>();
        m->coverSamples = cover;
        const double us = tickMicros (*m, 16384, 14, ticks);
        std::printf ("    multi, hop %-5d %8.0f us/tick  (%5.2f %% of a 30 fps frame)\n", cover, us, 100.0 * us / kTickBudgetUs);
        ok (us < 0.20 * kTickBudgetUs, "multi at hop " + std::to_string (cover) + ": a tick stays under 20 % of the frame (" + std::to_string ((int) us) + " us)");
        if (cover == 1600)
            ok (us < 2.0 * (classicUs[4] + classicUs[2] + 3.0 * classicUs[0]),
                "multi at hop 1600 costs less than twice its tiers' classic panes (16384 + 4096 + 3×1024) — the stitching is not the cost");
    }

    return test::report();
}

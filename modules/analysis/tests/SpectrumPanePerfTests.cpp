// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// felitronics::analysis — PERFORMANCE of the spectrum panes, on the scalar reference FFT: what one UI
// tick costs (ingest + 900 log-frequency columns) for the classic pane at every order and for the
// multi-resolution pane at a normal hop and after a 100 ms stall. The numbers are printed for the record
// (a machine-dependent table, compare across commits on one machine); the CHECKS are loose ceilings that
// only a regression of an order of magnitude would trip — a tick must stay far below the 33 ms frame,
// and the multi-res pane must cost less than the sum of its tiers' classic panes plus the stitching.
// The last group pairs the multi-res pane with its fast sibling; on this backend the FFT dominates, so
// the ratio is modest here and the pffft suite is where it shows.
// The pffft counterpart (same panes on the SIMD backend, and the speed-up) lives in the pffft suite.

#include <felitronics_test.h>
#include <felitronics/analysis/SpectrumPane.h>
#include <felitronics/analysis/MultiResSpectrumPane.h>
#include <felitronics/analysis/MultiResSpectrumPaneFast.h>

#include <felitronics/core/Fft.h>

#include <array>
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

// A counting FFT backend: the scalar reference verbatim, plus a tally of transforms by size. It exists so
// that one claim in this file can be checked WITHOUT a clock. `Tag` only gives each instantiation its own
// tally — the panes hold their backends by value, so the counter has to be static.
template <int Tag>
class CountingFft
{
public:
    using Inner = core::fft::ScalarRadix2Real;
    static constexpr bool kPackedHermitianSpectrum = Inner::kPackedHermitianSpectrum;
    static constexpr int  spectrumFloats (int n) noexcept { return Inner::spectrumFloats (n); }

    static constexpr int kMaxLog2 = 24;          // wider than any order a pane template permits
    static inline std::array<long long, kMaxLog2> tally {};      // [log2(n)] -> transforms of that size
    static void resetTally() noexcept { tally.fill (0); }

    bool prepare (int n) noexcept { n_ = n; return inner_.prepare (n); }
    void forward (const float* r, float* w) noexcept { inner_.forward (r, w); bump(); }    // count COMPLETED
    void inverse (const float* sp, float* w) noexcept { inner_.inverse (sp, w); bump(); }   // transforms
    void spectralMultiplyAdd (const float* a, const float* b, float* acc) noexcept { inner_.spectralMultiplyAdd (a, b, acc); }

private:
    // Bucket k is the smallest power of two >= n_, so it is the exact size for the powers of two a radix-2
    // backend accepts. An unprepared instance (n_ == 0) lands in bucket 0, which the schedule assertions
    // then report as transforms missing from every expected size rather than silently matching.
    void bump() noexcept { int k = 0; while ((1LL << k) < n_ && k < kMaxLog2 - 1) ++k; ++tally[(std::size_t) k]; }
    Inner inner_;
    int   n_ = 0;
};

// One tick through a pane, for the tally rather than the clock.
template <class Pane>
void oneTick (Pane& p, int frameSamples, int order)
{
    analysis::PlotMap pm; pm.width = 900.0f; pm.height = 300.0f; pm.plotBottom = 300.0f;
    volatile float sink = 0.0f;
    float* f = p.frameInput(); for (int i = 0; i < frameSamples; ++i) f[i] = uni();
    p.ingest (order);
    p.buildColumns (pm, 48000.0, 1.5, 1000.0, [&] (int, float, float y, float yp) { sink = sink + y + yp; });
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
        // The one check left in this file that compares SEPARATELY TIMED phases. Margin on a dev Mac is
        // about 2x (384 us against a 794 us threshold) — two orders wider than the 1.05x ratio that used to
        // flake below, but not unconditional: contention landing on this measurement and not on the earlier
        // classic ones is how it would go red. Recorded so the next reader has the number rather than a
        // hunch; deliberately NOT widened, because it has never actually been observed to fail.
        if (cover == 1600)
            ok (us < 2.0 * (classicUs[4] + classicUs[2] + 3.0 * classicUs[0]),
                "multi at hop 1600 costs less than twice its tiers' classic panes (16384 + 4096 + 3×1024) — the stitching is not the cost");
    }

    //==========================================================================================
    // The fast sibling, same ladder, same tick.
    //
    // WHAT IS ASSERTED HERE AND WHAT IS NOT. The saving this sibling exists for — a log10 and an exp per
    // bin, and the column geometry derived once instead of twice — is scalar work AROUND the transform.
    // On this backend the transform dominates, so the difference lands inside the noise of a shared CI
    // machine, and an assertion on the RATIO of two wall-clock timings measures the machine, not the code:
    // `b < 1.05 * a` used to live here and failed roughly one CI run in ten with nothing changed
    // (measured on main: 0.918889x on 2026-09-04, 0.906340x on 2026-09-02, and 0.890892x on a branch whose
    // sibling job on the same commit passed). A check that only ever fires falsely is worse than no check:
    // it teaches everyone that red means nothing.
    //
    // So the COST claim lives where the effect is large enough to measure — the pffft suite, where the
    // sibling is 1.8-2.3x cheaper and the bar is `uf < 0.85 * up` with margin to spare. BE CLEAR ABOUT WHAT
    // THAT COSTS: that suite is built only with FELITRONICS_WITH_PFFFT=ON, which is OFF by default, so in a
    // default build nothing here guards the sibling's defining cost property. It is guarded at the level of
    // the REPOSITORY, by the required pffft=ON jobs on Linux and macOS in .github/workflows/ci.yml — not by
    // this file. The group below is deliberately a WEAKER claim, not a substitute: it pins the FFT schedule,
    // which is what makes any timing comparison meaningful in the first place and catches a sibling that
    // quietly starts transforming more than the pane it copies. It says nothing about the scalar work that
    // IS the saving — a sibling that reintroduced the per-bin log10/exp would pass it untouched. The
    // timings are still printed, for the record.
    group ("multi-resolution: the fast sibling against the pane it copies");
    {
        auto cur  = std::make_unique<analysis::MultiResSpectrumPane>();
        auto fast = std::make_unique<analysis::MultiResSpectrumPaneFast>();
        cur->coverSamples = 1600; fast->coverSamples = 1600;
        const double a = tickMicros (*cur, 16384, 14, ticks);
        const double b = tickMicros (*fast, 16384, 14, ticks);
        std::printf ("    multi current %8.0f us/tick   multi FAST %8.0f us/tick   (%.2fx)\n", a, b, a / b);
        ok (b < 0.20 * kTickBudgetUs, "a fast tick stays under 20 % of a 30 fps frame (" + std::to_string ((int) b) + " us)");
    }

    // Equality between the two panes is NOT enough on its own: anything that moved both of them the same
    // way — a different tier ladder, a different hop cadence, a prepare() that silently stopped working —
    // would keep them equal and slip through, and `total > 0` would not notice. So the EXPECTED schedule is
    // asserted for each pane separately, and equality follows from that. The numbers are the ladder this
    // file already names one assertion earlier: per tick, one 16384 + one 4096 + three 1024 at hop 1600.
    // Nothing golden about them — they are the documented design, written down where a change to it has to
    // be acknowledged rather than absorbed.
    group ("the panes' transform schedule is exactly the documented ladder (no clock involved)");
    {
        constexpr int kOrder = analysis::RollingSpectrumTap::kMaxOrder;
        constexpr int kTicks = 8;
        using CurCount  = analysis::MultiResSpectrumPaneT     <kOrder, 4, CountingFft<0>>;
        using FastCount = analysis::MultiResSpectrumPaneFastT <kOrder, 4, CountingFft<1>>;

        struct Expect { int size, perTick; };
        constexpr Expect kLadder[] { { 1024, 3 }, { 4096, 1 }, { 16384, 1 } };

        CountingFft<0>::resetTally();
        CountingFft<1>::resetTally();
        {
            auto c = std::make_unique<CurCount>();  c->coverSamples = 1600;
            auto f = std::make_unique<FastCount>(); f->coverSamples = 1600;
            for (int t = 0; t < kTicks; ++t) { oneTick (*c, 16384, 14); oneTick (*f, 16384, 14); }
        }

        long long totalCur = 0, totalFast = 0, expected = 0;
        std::string shape;
        for (std::size_t k = 0; k < (std::size_t) CountingFft<0>::kMaxLog2; ++k)
        {
            const long long x = CountingFft<0>::tally[k], y = CountingFft<1>::tally[k];
            totalCur += x; totalFast += y;
            if (x != 0 || y != 0)
                shape += " " + std::to_string (1LL << k) + ":" + std::to_string (x) + "/" + std::to_string (y);
        }
        std::printf ("    transforms over %d ticks (size:current/fast) —%s\n", kTicks, shape.c_str());

        for (const Expect& e : kLadder)
        {
            int k = 0; while ((1 << k) < e.size) ++k;
            const long long want = (long long) e.perTick * kTicks;
            expected += want;
            ok (CountingFft<0>::tally[(std::size_t) k] == want,
                "current pane: " + std::to_string (e.perTick) + " x " + std::to_string (e.size) + " per tick");
            ok (CountingFft<1>::tally[(std::size_t) k] == want,
                "fast sibling: " + std::to_string (e.perTick) + " x " + std::to_string (e.size) + " per tick");
        }
        ok (totalCur == expected && totalFast == expected,
            "and neither pane transforms anything else, at any other size (" + std::to_string (totalCur)
                + " / " + std::to_string (totalFast) + " of " + std::to_string (expected) + ")");
    }

    return test::report();
}

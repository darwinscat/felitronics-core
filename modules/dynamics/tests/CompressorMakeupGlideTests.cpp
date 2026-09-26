// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// The Compressor's MAKEUP GLIDE (Compressor::kMakeupGlideMs; "A CHANGE OF MAKEUP GLIDES" in the header). Each group
// pins one claim, exactly where it can be exact:
//
//   * THE FIRST WRITE SNAPS: prepare -> setParams -> process renders what setParams -> prepare -> process renders.
//   * A MID-STREAM MAKEUP MOVE GLIDES linearly in dB, one step per processed sample, landing on the target: pinned
//     sample for sample on material below the threshold, where the curve adds 0 dB and the output IS the delayed
//     input times dbToGain(makeup) — the replayed ramp (double accumulator) must match bit for bit.
//   * autoMakeup AND THE CURVE MOVES THAT CHANGE IT glide the same way.
//   * CLOCKED BY SAMPLES (law 8a): a timeline of makeup writes with a clock-only pause renders the same bits under any
//     cut, and the pause spends the glide as the silence it stands for would.
//   * AN UNCHANGED WRITE RESTARTS NOTHING, and nothing is allocated.

#include <felitronics_test.h>
#include <alloc_counter.h>   // installs the allocation counter: EVERY form of `new`, over-aligned included
#include <felitronics/dynamics/Compressor.h>

#include <bit>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

using namespace felitronics;
using felitronics::test::ok;
using felitronics::test::group;
using dynamics::Compressor;
using dynamics::CompressorParams;

namespace
{
constexpr double kFs = 48000.0;
using Buf = std::vector<std::vector<float>>;
std::uint32_t bits (float f) noexcept { return std::bit_cast<std::uint32_t> (f); }
long long diffs (const Buf& a, const Buf& b, int from = 0)
{
    long long d = 0;
    for (std::size_t c = 0; c < a.size(); ++c)
        for (std::size_t i = (std::size_t) from; i < a[c].size(); ++i) d += bits (a[c][i]) != bits (b[c][i]);
    return d;
}

CompressorParams P (double makeup, bool autoMk = false, double thr = -18.0)
{
    CompressorParams p;
    p.thresholdDb = thr; p.ratio = 4.0; p.attackMs = 5.0; p.releaseMs = 80.0; p.makeupDb = makeup; p.autoMakeup = autoMk;
    p.lookaheadMs = 1.0;
    return p;
}

Buf tone (int n, double amp)
{
    Buf x (2, std::vector<float> ((std::size_t) n));
    for (int i = 0; i < n; ++i)
        for (int c = 0; c < 2; ++c) x[(std::size_t) c][(std::size_t) i] = (float) (amp * std::sin (0.031 * i + 0.4 * c));
    return x;
}

struct Write { int at; CompressorParams p; };
struct Gap { int from, to; };

template <class Cut>
Buf render (const Buf& x, const CompressorParams& p0, const std::vector<Write>& ws, Cut cut, const std::vector<Gap>& gaps = {},
            bool setBeforePrepare = true)
{
    Buf y = x;
    Compressor c;
    if (setBeforePrepare) c.setParams (p0);
    if (! c.prepare (kFs, 512, 2, 50.0)) return {};
    if (! setBeforePrepare) c.setParams (p0);
    const int n = (int) x[0].size();
    std::size_t w = 0;
    for (int off = 0, k = 0; off < n; ++k)
    {
        while (w < ws.size() && ws[w].at <= off) c.setParams (ws[w++].p);
        int len = std::min (cut (k), n - off);
        if (w < ws.size()) len = std::min (len, ws[w].at - off);
        bool paused = false;
        for (const Gap& g : gaps)
        {
            if (off >= g.from && off < g.to) { paused = true; len = std::min (len, g.to - off); }
            else if (off < g.from) len = std::min (len, g.from - off);
        }
        float* io[2] { y[0].data() + off, y[1].data() + off };
        if (len == 0) { felitronics::test::run (c.process (io, 2, 0)); continue; }
        if (paused) { felitronics::test::run (c.process (nullptr, 0, len)); for (auto& ch : y) std::fill (ch.begin() + off, ch.begin() + off + len, 0.0f); }
        else felitronics::test::run (c.process (io, 2, len));
        off += len;
    }
    return y;
}
auto fixed (int b) { return [b] (int) { return b; }; }
auto ragged (std::uint32_t seed)
{
    return [seed] (int k)
    {
        std::uint32_t s = seed + 2654435761u * (std::uint32_t) (k + 1);
        s ^= s >> 13; s *= 0x5bd1e995u; s ^= s >> 15;
        return (s % 7u == 0u) ? 0 : 1 + (int) (s % 900u);
    };
}
} // namespace

static void testFirstWriteSnaps()
{
    group ("the first write snaps — prepare -> set equals set -> prepare");
    const Buf x = tone (20000, 0.5);
    for (const CompressorParams& p : { P (6.0), P (0.0, true), P (-3.0, true, -30.0) })
    {
        const Buf a = render (x, p, {}, fixed (512), {}, true);
        const Buf b = render (x, p, {}, fixed (512), {}, false);
        ok (! a.empty() && diffs (a, b) == 0, "the two orders render the same bits (makeup " + std::to_string ((int) p.makeupDb)
                                               + (p.autoMakeup ? ", auto" : "") + ")");
    }
}

static void testMakeupGlideExact()
{
    group ("a mid-stream makeup move glides per sample and lands — replayed bit for bit below the threshold");
    // -40 dBFS against a -18 dB threshold with a 6 dB knee: the curve adds exactly 0 dB, so out = x[i - look] * g(i).
    const int n = 12000, at = 3001;
    const Buf x = tone (n, 0.01);
    const Buf y = render (x, P (0.0), { { at, P (6.0) } }, ragged (9));
    Compressor probe;
    ok (probe.prepare (kFs, 512, 2, 50.0), "PRECONDITION: prepare");
    probe.setParams (P (0.0));
    const int look = probe.latencySamples(), L = (int) std::lround (Compressor::kMakeupGlideMs * 0.001 * kFs);
    ok (look == 48 && L == 1440, "PRECONDITION: 48 samples of lookahead, a 1440-sample glide");
    long long bad = 0;
    double cur = 0.0, step = 0.0; int left = 0;
    for (int i = 0; i < n; ++i)
    {
        if (i == at) { left = L; step = (6.0 - cur) / (double) L; }
        if (left > 0) { --left; cur = left > 0 ? cur + step : 6.0; }
        const float g = (float) core::dbToGain (cur);
        for (int c = 0; c < 2; ++c)
        {
            const float want = i < look ? 0.0f : x[(std::size_t) c][(std::size_t) (i - look)] * g;
            bad += bits (want) != bits (y[(std::size_t) c][(std::size_t) i]);
        }
    }
    ok (bad == 0, "every sample is the delayed input times the replayed makeup, the landing included (" + std::to_string (bad) + " differ)");
}

static void testAutoMakeupGlides()
{
    group ("autoMakeup, and a curve move that changes it, glide the same way");
    const int n = 20000, at = 8000 + 13;
    const Buf x = tone (n, 0.25);
    auto maxD2 = [] (const std::vector<float>& v, int a, int b)
    {
        double m = 0.0;
        for (int i = std::max (a, 2); i < b; ++i) m = std::max (m, std::fabs ((double) v[(std::size_t) i] - 2.0 * v[(std::size_t) i - 1] + v[(std::size_t) i - 2]));
        return m;
    };
    // A THRESHOLD move also moves the gain reduction, through the attack ballistics — a kink of its own that is not the
    // makeup's (5 ms of attack on a 9 dB move reads ~3x the steady Δ² here with autoMakeup OFF). So that case is held to
    // the edge the same move makes without autoMakeup: the makeup's own glide must add next to nothing to it.
    struct Case { const char* name; CompressorParams a, b; CompressorParams refA, refB; };
    const Case cases[] = { { "makeup 0 -> +3", P (0.0), P (3.0), {}, {} },
                           { "autoMakeup off -> on", P (0.0), P (0.0, true), {}, {} },
                           { "threshold -18 -> -30 with autoMakeup", P (0.0, true), P (0.0, true, -30.0), P (0.0), P (0.0, false, -30.0) } };
    for (const Case& c : cases)
    {
        const Buf y = render (x, c.a, { { at, c.b } }, fixed (128));
        const double steady = std::max (maxD2 (y[0], 2000, 7000), maxD2 (y[0], at + 12000 - 4000, n));
        const double edge = maxD2 (y[0], at, at + 4000);
        double bound = 2.0 * steady;
        if (c.refA.thresholdDb != c.refB.thresholdDb)
        {
            // …scaled by the LEVEL: autoMakeup puts 6.75 dB on the output before the move, and a kink's Δ² scales
            // with the amplitude it bends.
            const Buf r = render (x, c.refA, { { at, c.refB } }, fixed (128));
            const double level = maxD2 (y[0], 2000, 7000) / maxD2 (r[0], 2000, 7000);
            bound = std::max (bound, 1.2 * level * maxD2 (r[0], at, at + 4000));
        }
        std::printf ("    %-38s steady %.2e  edge %.2e  bound %.2e\n", c.name, steady, edge, bound);
        ok (edge < bound, std::string ("no step from the makeup — ") + c.name);
    }
}

static void testClockedBySamples()
{
    group ("law 8a — makeup writes with a clock-only pause render the same bits under any cut");
    const int n = 16000;
    const Buf x = tone (n, 0.3);
    const std::vector<Write> ws { { 1234, P (5.0) }, { 1700, P (-2.0, true) }, { 7777, P (1.0, true, -26.0) } };
    const std::vector<Gap> gaps { { 1800, 2600 } };                  // a pause while the retargeted glide runs
    const Buf ref = render (x, P (0.0), ws, fixed (n), gaps);
    for (const auto& cut : { std::function<int (int)> (fixed (1)), std::function<int (int)> (fixed (64)),
                             std::function<int (int)> (ragged (1)), std::function<int (int)> (ragged (2)) })
    {
        const Buf got = render (x, P (0.0), ws, cut, gaps);
        ok (! got.empty() && diffs (ref, got) == 0, "bit-identical (" + std::to_string (got.empty() ? -1 : diffs (ref, got)) + " differ)");
    }
    // The pause spends the glide: a makeup move written right before a pause longer than the glide has landed when
    // the audio returns — the same output as if the move had been written before the stream started.
    const Buf mv = render (x, P (0.0), { { 4000, P (6.0) } }, fixed (256), { { 4000, 6000 } });
    const Buf snap = render (x, P (6.0), {}, fixed (256), { { 4000, 6000 } });
    ok (diffs (mv, snap, 6000 + 48) == 0, "a glide spanned by a pause has landed when the audio returns");
}

// ABSURD BUT FINITE MAKEUPS stay finite through a glide: -1e308 -> +1e308 -> 0 overflowed the ramp's difference to inf
// and the retarget to NaN, which the ±400 dB sum clamp passes (the code-review round).
static void testAbsurdMakeupsStayFinite()
{
    group ("absurd but finite makeups glide without leaving the finite numbers");
    const int n = 12000;
    const Buf x = tone (n, 0.3);
    const Buf y = render (x, P (-1.0e308), { { 1000, P (1.0e308) }, { 1500, P (0.0) }, { 5000, P (-1.0e308) }, { 5100, P (3.0) } },
                          fixed (100));
    bool finite = true;
    for (const auto& c : y) for (float v : c) finite = finite && std::isfinite (v);
    ok (! y.empty() && finite, "every output sample is finite");
    const Buf ref = render (x, P (3.0), {}, fixed (100));
    ok (diffs (y, ref, 5100 + 1440 + 48) == 0, "…and past the last glide the makeup is exactly the target's");
}

static void testUnchangedWritesAndAllocation()
{
    group ("an unchanged write restarts nothing, and nothing is allocated");
    const int n = 8000;
    const Buf x = tone (n, 0.3);
    std::vector<Write> again;
    for (int k = 1; k < 40; ++k) again.push_back ({ 64 * k, P (4.0) });
    const Buf a = render (x, P (0.0), { { 64, P (4.0) } }, fixed (64));
    const Buf b = render (x, P (0.0), again, fixed (64));
    ok (diffs (a, b) == 0, "re-sending the same parameters every block leaves the glide exactly as one write did");

    Compressor c;
    c.setParams (P (0.0));
    ok (c.prepare (kFs, 128, 2, 50.0), "PRECONDITION: prepare");
    Buf y = tone (128, 0.3);
    float* io[2] { y[0].data(), y[1].data() };
    felitronics::test::run (c.process (io, 2, 128));
    const long long before = alloc::count.load();
    for (int k = 0; k < 300; ++k)
    {
        if (k % 11 == 0) c.setParams (P ((double) (k % 7), (k % 22) == 0));
        felitronics::test::run (c.process (io, 2, 128));
        if (k == 100) felitronics::test::run (c.process (nullptr, 0, 3000));
    }
    felitronics::test::okNoAlloc (alloc::count.load() == before, "no allocation across 300 gliding blocks and a pause");
}

int main()
{
    std::printf ("felitronics::dynamics — the Compressor's makeup glide\n");
    testFirstWriteSnaps();
    testMakeupGlideExact();
    testAutoMakeupGlides();
    testClockedBySamples();
    testAbsurdMakeupsStayFinite();
    testUnchangedWritesAndAllocation();
    return felitronics::test::report();
}

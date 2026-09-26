// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// MonoBass's live moves (see LIVE MOVES GLIDE in MonoBass.h): the crossover corner and the air corner glide on the
// 64-sample grid, `enabled` rides the xf crossfade, the air's `enabled` rides its plateau — and the first write of a
// stream snaps. Each group pins one of those claims:
//
//   * THE FIRST WRITE SNAPS, and the two orders a caller can configure in are ONE behaviour: prepare -> set equals
//     set -> prepare, bit for bit, for every parameter including lowWidth and the air plateau (which used to glide
//     from the prepared value in the first order).
//   * EVERY GLIDE IS ON AUDIO TIME (law 8a): a timeline of writes — both corners, both enables, lowWidth, the plateau —
//     renders the same bits however the stream is cut, with zero-length calls, a mono stretch and a clock-only pause
//     in it; and the corners keep walking through a bypassed stretch exactly as through audio.
//   * `enabled` IS A FADE: switching off leaves a bit-exact passthrough once the fade has settled, and neither edge
//     has a second difference above the steady tone's.
//   * THE setAir DEFECT: a live corner move restarts the width-measurement interval (it used to keep summing across
//     both bands — 2000 judged samples where the new band had judged 1000).
//   * NOTHING IS ALLOCATED.

#include <felitronics_test.h>
#include <alloc_counter.h>   // installs the allocation counter: EVERY form of `new`, over-aligned included
#include <felitronics/stereo/MonoBass.h>

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
using stereo::MonoBass;
using stereo::MonoBassParams;
using stereo::StereoAirParams;

namespace
{
constexpr double kFs = 48000.0, kPi = 3.14159265358979323846;
using Buf = std::vector<std::vector<float>>;

std::uint32_t bits (float f) noexcept { return std::bit_cast<std::uint32_t> (f); }
long long diffs (const Buf& a, const Buf& b, int from = 0, int to = -1)
{
    long long d = 0;
    for (std::size_t c = 0; c < a.size(); ++c)
        for (std::size_t i = (std::size_t) from; i < (to < 0 ? a[c].size() : (std::size_t) to); ++i) d += bits (a[c][i]) != bits (b[c][i]);
    return d;
}

// Wide low end (L and R 90 degrees apart at 103.7 Hz: a large Side the crossover acts on) plus wide top (5 kHz) the air
// shelf acts on, plus a little noise.
Buf programme (int n, unsigned seed)
{
    Buf x (2, std::vector<float> ((std::size_t) n));
    std::uint32_t s = seed;
    for (int i = 0; i < n; ++i)
    {
        s = s * 1664525u + 1013904223u;
        const float nz = ((float) (s >> 8) / 16777216.0f - 0.5f) * 0.05f;
        const double t = i / kFs;
        x[0][(std::size_t) i] = (float) (0.3 * std::sin (2 * kPi * 103.7 * t) + 0.1 * std::sin (2 * kPi * 5003.1 * t)) + nz;
        x[1][(std::size_t) i] = (float) (0.3 * std::sin (2 * kPi * 103.7 * t + kPi / 2) + 0.1 * std::sin (2 * kPi * 5003.1 * t + kPi / 2)) - nz;
    }
    return x;
}

struct Set { MonoBassParams b; StereoAirParams a; };
struct Write { int at; Set s; };
struct Span { int from, to, width; };            // width 0: a clock-only pause, 1: a mono stretch

void apply (MonoBass& m, const Set& s) { m.setParams (s.b); m.setAir (s.a); }

// A fresh stage, `s0` written before the first sample, the writes at their positions, the stream cut by `cut(k)`.
template <class Cut>
Buf render (const Buf& x, const Set& s0, const std::vector<Write>& ws, Cut cut, const std::vector<Span>& spans = {},
            MonoBass* keep = nullptr, bool setBeforePrepare = true)
{
    Buf y = x;
    MonoBass local;
    MonoBass& m = keep != nullptr ? *keep : local;
    if (setBeforePrepare) apply (m, s0);
    if (! m.prepare (kFs, 4096, 2)) return {};
    if (! setBeforePrepare) apply (m, s0);
    const int n = (int) x[0].size();
    std::size_t w = 0;
    for (int off = 0, k = 0; off < n; ++k)
    {
        while (w < ws.size() && ws[w].at <= off) apply (m, ws[w++].s);
        int len = std::min (cut (k), n - off);
        if (w < ws.size()) len = std::min (len, ws[w].at - off);
        int width = 2;
        for (const Span& sp : spans)
        {
            if (off >= sp.from && off < sp.to) { width = sp.width; len = std::min (len, sp.to - off); }
            else if (off < sp.from) len = std::min (len, sp.from - off);
        }
        float* io[2] { y[0].data() + off, y[1].data() + off };
        if (len == 0) { felitronics::test::run (m.process (io, 2, 0)); continue; }
        felitronics::test::run (m.process (width == 0 ? nullptr : io, width, len));
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

// The same wide tones WITHOUT the noise: a second difference measures a click only where nothing else moves it (the
// code-review round found the noisy programme's own Δ² hiding a hard step entirely).
Buf tones (int n)
{
    Buf x (2, std::vector<float> ((std::size_t) n));
    for (int i = 0; i < n; ++i)
    {
        const double t = i / kFs;
        x[0][(std::size_t) i] = (float) (0.3 * std::sin (2 * kPi * 103.7 * t));
        x[1][(std::size_t) i] = (float) (0.3 * std::sin (2 * kPi * 103.7 * t + kPi / 2));
    }
    return x;
}

Set base()
{
    Set s;
    s.b = { true, 120.0f, 0.0f };
    s.a.enabled = true; s.a.frequencyHz = 6000.0f; s.a.gainDb = 3.0f;
    return s;
}

double maxD2 (const std::vector<float>& y, int a, int b)
{
    double m = 0.0;
    for (int i = std::max (a, 2); i < b; ++i)
        m = std::max (m, std::fabs ((double) y[(std::size_t) i] - 2.0 * y[(std::size_t) i - 1] + y[(std::size_t) i - 2]));
    return m;
}
} // namespace

//==============================================================================
static void testFirstWriteSnapsAndTheOrdersAgree()
{
    group ("the first write snaps — prepare -> set equals set -> prepare, bit for bit, every parameter");
    const Buf x = programme (24000, 3u);
    struct Case { const char* name; std::function<void (Set&)> mut; };
    const Case cases[] = {
        { "lowWidth 0 -> 0.5",           [] (Set& s) { s.b.lowWidth = 0.5f; } },
        { "lowWidth 0 -> 1 (full-wide)", [] (Set& s) { s.b.lowWidth = 1.0f; } },
        { "frequency 120 -> 240",        [] (Set& s) { s.b.frequencyHz = 240.0f; } },
        { "enabled off",                 [] (Set& s) { s.b.enabled = false; } },
        { "air plateau 3 -> 5",          [] (Set& s) { s.a.gainDb = 5.0f; } },
        { "air corner 6000 -> 9000",     [] (Set& s) { s.a.frequencyHz = 9000.0f; } },
        { "air off",                     [] (Set& s) { s.a.enabled = false; } },
    };
    for (const Case& c : cases)
    {
        Set s = base(); c.mut (s);
        const Buf a = render (x, s, {}, fixed (512), {}, nullptr, true);
        const Buf b = render (x, s, {}, fixed (512), {}, nullptr, false);    // prepare(defaults) -> set -> process
        ok (! a.empty() && diffs (a, b) == 0, std::string ("set -> prepare and prepare -> set render the same bits — ") + c.name);
    }
    // …and after a reset() that interrupts glides of every kind, the next write snaps.
    MonoBass m;
    Set s0 = base();
    const Buf ref = render (x, s0, {}, fixed (512));
    apply (m, base());
    ok (m.prepare (kFs, 4096, 2), "PRECONDITION: prepare");
    Buf junk = programme (3100, 99u);                        // some other programme streamed before the restart
    float* jo[2] { junk[0].data(), junk[1].data() };
    felitronics::test::run (m.process (jo, 2, 3000));
    Set moved = base(); moved.b = { false, 300.0f, 0.7f }; moved.a.frequencyHz = 10000.0f; moved.a.gainDb = 6.0f;
    apply (m, moved);
    float* jo2[2] { junk[0].data() + 3000, junk[1].data() + 3000 };
    felitronics::test::run (m.process (jo2, 2, 100));
    m.reset();
    apply (m, s0);
    Buf y = x;
    float* io[2] { y[0].data(), y[1].data() };
    felitronics::test::run (m.process (io, 2, 24000));
    ok (diffs (ref, y) == 0, "after a reset() mid-glide the next write snaps: identical to a fresh stage");
}

static void testGlidesAreOnAudioTime()
{
    group ("law 8a — a timeline of writes renders the same bits under every cut, a mono stretch and a pause included");
    const Buf x = programme (40000, 7u);
    Set s = base();
    std::vector<Write> ws;
    s.b.frequencyHz = 250.0f;                       ws.push_back ({ 3001, s });
    s.b.frequencyHz = 70.0f;                        ws.push_back ({ 3500, s });   // retarget mid-glide
    s.b.enabled = false;                            ws.push_back ({ 7777, s });
    s.b.enabled = true; s.b.lowWidth = 0.4f;        ws.push_back ({ 11111, s });
    s.a.frequencyHz = 3000.0f;                      ws.push_back ({ 15000, s });
    s.a.enabled = false;                            ws.push_back ({ 19003, s });
    s.a.enabled = true; s.a.gainDb = 6.0f;          ws.push_back ({ 23456, s });
    s.b.frequencyHz = 180.0f; s.a.frequencyHz = 11000.0f; ws.push_back ({ 27100, s });   // both corners, then a pause
    const std::vector<Span> spans { { 28000, 30500, 0 }, { 33000, 34000, 1 } };
    const Buf ref = render (x, base(), ws, fixed (40000), spans);
    const Buf still = render (x, base(), {}, fixed (40000), spans);
    ok (! ref.empty() && diffs (ref, still, 3100) > 10000, "PRECONDITION: the writes are audible");
    struct Row { const char* name; std::function<int (int)> cut; };
    const Row rows[] = { { "one sample at a time", fixed (1) }, { "blocks of 64", fixed (64) }, { "blocks of 100", fixed (100) },
                         { "ragged with zero-length calls #1", ragged (1) }, { "ragged #2", ragged (2) }, { "blocks of 4096", fixed (4096) } };
    for (const Row& r : rows)
    {
        const Buf got = render (x, base(), ws, r.cut, spans);
        const long long d = got.empty() ? -1 : diffs (ref, got);
        ok (d == 0, std::string ("bit-identical: ") + r.name + " (" + std::to_string (d) + " differ)");
    }
    // A corner move written while the island is IDLE (the bass settled full-wide, the air off) walks through the
    // skipped stretch on the same boundaries the audio path would have met, so the render after the island wakes up is
    // the same bits under any cut.
    {
        Set idle = base(); idle.b.lowWidth = 1.0f; idle.a.enabled = false;
        Set moving = idle; moving.b.frequencyHz = 300.0f;
        Set back = moving; back.b.lowWidth = 0.0f;
        const std::vector<Write> w2 { { 2000, moving }, { 2500, back } };
        const Buf a = render (x, idle, w2, fixed (256));
        const Buf b = render (x, idle, w2, fixed (777));
        const Buf c = render (x, idle, w2, ragged (3));
        ok (diffs (a, b) == 0 && diffs (a, c) == 0, "a corner glide started while the island is idle lands on the same boundaries under any cut");
    }
}

static void testEnabledIsAFade()
{
    group ("`enabled` is a fade — off settles into a bit-exact passthrough, and neither edge steps");
    const int n = 30000, at = 12000 + 17;
    const Buf x = tones (n);
    Set on = base(); on.a.enabled = false; on.b.frequencyHz = 250.0f;
    Set off = on; off.b.enabled = false;
    const Buf goOff = render (x, on, { { at, off } }, fixed (128));
    // Past the 20 ms fade (960 samples) the bass has retired and the island is skipped: the input, bit for bit.
    ok (diffs (goOff, x, at + 64 + 960 + 64) == 0, "switched off: a bit-exact passthrough once the fade has settled");
    ok (diffs (goOff, x, at, at + 960) > 0, "PRECONDITION: and not before — it faded");
    const Buf goOn = render (x, off, { { at, on } }, fixed (128));
    ok (diffs (goOn, x, 0, at) == 0, "PRECONDITION: disabled before the write is a passthrough");
    // THE CRITERION IS THE HARNESS'S, in absolute terms: max|Δ²y| under -60 dBFS is clean (a -10.5 dBFS 104 Hz tone's own
    // Δ² is -85 dBFS, and a LINEAR fade's kink necessarily sits above that — amplitude times one step, the level the
    // shipped lowWidth fade has always had), and a hard switch here reads above -40.
    constexpr double kClean = 1.0e-3, kClick = 1.0e-2;
    for (const auto* y : { &goOff, &goOn })
    {
        const double edge = maxD2 ((*y)[0], at - 64, at + 4000);
        ok (edge < kClean, std::string (y == &goOff ? "off" : "on") + ": the edge stays under -60 dBFS of Δ² ("
                           + std::to_string (20.0 * std::log10 (edge)) + " dBFS)");
    }
    // And the corner: a 60 -> 250 Hz move glides — the worst Δ² around it stays at the steady tone's.
    Set lo = base(); lo.a.enabled = false; lo.b.frequencyHz = 60.0f;
    Set hi = lo; hi.b.frequencyHz = 250.0f;
    const Buf mv = render (x, lo, { { at, hi } }, fixed (128));
    ok (maxD2 (mv[0], at - 64, at + 6000) < kClean, "the crossover corner glides 60 -> 250 Hz under -60 dBFS of Δ² ("
                                                    + std::to_string (20.0 * std::log10 (maxD2 (mv[0], at - 64, at + 6000))) + " dBFS)");
    // …and the check can see a step: the same move written into a stage whose stream has not started (it SNAPS)
    // and spliced — the pre-write half from `lo`, the post-write half from a stage that had `hi` all along — is a
    // hard switch between two settled renders, and its Δ² is many times the steady tone's.
    const Buf always = render (x, hi, {}, fixed (128));
    std::vector<float> spliced (mv[0]);
    const Buf still = render (x, lo, {}, fixed (128));
    for (int i = 0; i < n; ++i) spliced[(std::size_t) i] = i < at ? still[0][(std::size_t) i] : always[0][(std::size_t) i];
    ok (maxD2 (spliced, at - 64, at + 64) > kClick, "PRECONDITION: a hard switch between the two corners IS visible to this check ("
                                                    + std::to_string (20.0 * std::log10 (maxD2 (spliced, at - 64, at + 64))) + " dBFS)");
}

// A GLIDE IN FLIGHT KEEPS MOVING THROUGH A PAUSE, a mono stretch or an idle island — the same samples of audio time,
// the same place: the designed corners and the crossfade after N samples of any of them equal those after N samples
// of stereo audio. (With the skip's corner ticks removed, the code-review round saw every earlier check still pass.)
static void testGlidesSpendSkippedTime()
{
    group ("a glide in flight spends a pause, a mono stretch and an idle island exactly as it spends audio");
    const Buf x = programme (8000, 21u);
    Set s0 = base();
    Set s1 = s0; s1.b.frequencyHz = 300.0f; s1.b.enabled = false; s1.a.frequencyHz = 11000.0f; s1.a.gainDb = 5.0f;
    struct Way { const char* name; int width; };
    float refXo = 0.0f, refAir = 0.0f, refXf = 0.0f;
    for (const Way w : { Way { "stereo audio", 2 }, Way { "a clock-only pause", 0 }, Way { "a mono stretch", 1 } })
    {
        MonoBass m;
        apply (m, s0);
        ok (m.prepare (kFs, 4096, 2), "PRECONDITION: prepare");
        Buf y = x;
        float* io[2] { y[0].data(), y[1].data() };
        felitronics::test::run (m.process (io, 2, 1000));
        apply (m, s1);
        felitronics::test::run (m.process (io, 2, 100));                     // the glides are in flight
        float* io2[2] { y[0].data() + 1100, y[1].data() + 1100 };
        felitronics::test::run (m.process (w.width == 0 ? nullptr : io2, w.width, 700));
        if (w.width == 2) { refXo = m.crossoverDesignHz(); refAir = m.airDesignHz(); refXf = m.crossfade(); }
        ok (bits (m.crossoverDesignHz()) == bits (refXo) && bits (m.airDesignHz()) == bits (refAir) && bits (m.crossfade()) == bits (refXf),
            std::string ("after 700 samples of ") + w.name + ", the corners and the crossfade stand where audio would have left them ("
            + std::to_string (m.crossoverDesignHz()) + " Hz, " + std::to_string (m.airDesignHz()) + " Hz, xf " + std::to_string (m.crossfade()) + ")");
        ok (m.crossoverDesignHz() > 150.0f && m.crossoverDesignHz() < 300.0f && m.crossfade() > 0.5f && m.crossfade() < 1.0f,
            std::string ("PRECONDITION: mid-glide, not before it and not landed — ") + w.name + " (" + std::to_string (m.crossoverDesignHz())
            + " Hz, xf " + std::to_string (m.crossfade()) + ")");
    }
    // A width written while the bass is disabled arrives while the air keeps the island running.
    MonoBass m;
    Set d = base(); d.b.enabled = false;
    apply (m, d);
    ok (m.prepare (kFs, 4096, 2), "PRECONDITION: prepare");
    Buf y = x;
    float* io[2] { y[0].data(), y[1].data() };
    felitronics::test::run (m.process (io, 2, 500));
    d.b.lowWidth = 0.7f;
    apply (m, d);
    felitronics::test::run (m.process (io, 2, 1500));
    ok (bits (m.lowWidthNow()) == bits (0.7f), "a width written while the bass is disabled has arrived 1500 samples later ("
                                              + std::to_string (m.lowWidthNow()) + ")");
}

// THE AIR RETIRES THE ISLAND ON THE SAMPLE CLOCK: with the bass idle, an air fade-out that lands mid-call must leave the
// rest of the call untouched — whole, per-sample and ragged renders the same bits — and a shelf that faded out with the
// bass in the same sample holds no tail for the next enable.
static void testAirRetiresTheIsland()
{
    group ("the air retires the island on the sample clock, and leaves no tail behind");
    const Buf x = programme (12000, 31u);
    Set s = base(); s.b.enabled = false;
    Set off = s; off.a.enabled = false;
    const std::vector<Write> ws { { 3001, off } };
    const Buf ref = render (x, s, ws, fixed (12000));
    for (const auto& cut : { std::function<int (int)> (fixed (1)), std::function<int (int)> (ragged (5)), std::function<int (int)> (fixed (777)) })
    {
        const Buf got = render (x, s, ws, cut);
        ok (diffs (ref, got) == 0, "bass idle, air fading out mid-call: the same bits under a different cut (" + std::to_string (diffs (ref, got)) + " differ)");
    }
    ok (diffs (ref, x, 3001 + 960 + 1) == 0, "past the fade the island is retired: the input, bit for bit");
    // Both fades landing on the same sample, then silence, then the air back on: nothing but silence comes out of it.
    Buf z = programme (12000, 33u);
    for (int i = 5000; i < 12000; ++i) z[0][(std::size_t) i] = z[1][(std::size_t) i] = 0.0f;
    Set both = base();
    Set bothOff = both; bothOff.b.enabled = false; bothOff.a.enabled = false;
    Set airBack = bothOff; airBack.a.enabled = true;
    for (const auto& cut : { std::function<int (int)> (fixed (4096)), std::function<int (int)> (fixed (1)) })
    {
        const Buf y = render (z, both, { { 3000, bothOff }, { 8000, airBack } }, cut);
        double tail = 0.0;
        for (int i = 8000; i < 12000; ++i) tail = std::max (tail, (double) std::fabs (y[0][(std::size_t) i]) + std::fabs (y[1][(std::size_t) i]));
        ok (tail == 0.0, "re-enabled out of silence after both fades landed together: exactly zero (" + std::to_string (tail) + ")");
    }
}

static void testSetAirRetunesOnALiveMove()
{
    group ("setAir — a live corner move restarts the width-measurement interval");
    const Buf x = programme (4000, 17u);
    MonoBass m;
    Set s = base();
    apply (m, s);
    ok (m.prepare (kFs, 4096, 2), "PRECONDITION: prepare");
    Buf y = x;
    float* io[2] { y[0].data(), y[1].data() };
    felitronics::test::run (m.process (io, 2, 1000));
    ok (m.airJudgedSamples() == 1000, "PRECONDITION: 1000 samples judged at 6 kHz");
    s.a.frequencyHz = 8000.0f;
    apply (m, s);
    ok (m.airJudgedSamples() == 0, "the move to 8 kHz restarts the interval at once");
    float* io2[2] { y[0].data() + 1000, y[1].data() + 1000 };
    felitronics::test::run (m.process (io2, 2, 1000));
    ok (m.airJudgedSamples() == 1000, "1000 judged in the new band, not 2000 across both (" + std::to_string (m.airJudgedSamples()) + ")");
    ok (bits (m.air().frequencyHz) == bits (8000.0f), "air() reads back the written corner at once, while the shelf glides to it");
    // A write that does not move the corner keeps the interval.
    s.a.gainDb = 4.0f;
    apply (m, s);
    ok (m.airJudgedSamples() == 1000, "a plateau-only write keeps the interval");
}

static void testNoAllocation()
{
    group ("RT — nothing is allocated while the corners and the fades move");
    MonoBass m;
    apply (m, base());
    ok (m.prepare (kFs, 256, 2), "PRECONDITION: prepare");
    Buf y = programme (256, 5u);
    float* io[2] { y[0].data(), y[1].data() };
    felitronics::test::run (m.process (io, 2, 256));
    const long long before = alloc::count.load();
    Set s = base();
    for (int k = 0; k < 300; ++k)
    {
        if (k % 13 == 0) { s.b.frequencyHz = 60.0f + (float) (k % 7) * 40.0f; s.b.enabled = (k % 26) != 0; s.a.frequencyHz = 3000.0f + (float) (k % 5) * 2000.0f; apply (m, s); }
        felitronics::test::run (m.process (io, 2, 256));
        if (k == 150) felitronics::test::run (m.process (nullptr, 0, 5000));
    }
    felitronics::test::okNoAlloc (alloc::count.load() == before, "no allocation across 300 moving blocks and a pause");
}

int main()
{
    std::printf ("felitronics::stereo — MonoBass live moves\n");
    testFirstWriteSnapsAndTheOrdersAgree();
    testGlidesAreOnAudioTime();
    testEnabledIsAFade();
    testGlidesSpendSkippedTime();
    testAirRetiresTheIsland();
    testSetAirRetunesOnALiveMove();
    testNoAllocation();
    return felitronics::test::report();
}

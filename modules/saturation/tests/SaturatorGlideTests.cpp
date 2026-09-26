// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// The Saturator's PARAMETER GLIDE (Saturator::kGlideMs). What it must be, and what each group below pins:
//
//   * THE FIRST WRITE OF A STREAM SNAPS — every write between prepare()/reset() and the first sample lands at once,
//     so a stage configured before its first sample renders exactly what it rendered before the glide existed. The
//     offline contract of the mastering chain stands on this.
//   * A WRITE MID-STREAM GLIDES, and the glide is CLOCKED BY AUDIO SAMPLES (law 8a): a timeline of writes at fixed
//     stream positions renders the same bits however the stream is cut — per sample, ragged, around maxBlock, with
//     zero-length and clock-only calls — and the same bits at any prepared maxBlock.
//   * IT LANDS: past its last period plus the downsampler's memory, a symmetric curve is bit-identical to a stage that
//     always had the target.
//   * IT DOES NOT CLICK: on a sine, the worst second difference inside a drive glide stays at the steady tone's own,
//     where the same move as a step is many times larger — a retarget mid-glide included.
//   * IT COSTS NOTHING WHEN NOTHING MOVES: a caller that re-sends the same parameters every block renders the same
//     bits as one that wrote them once, and never starts a glide. And glide 0 is the pre-glide stage (the frozen-engine
//     NULL in SaturationTests pins that bit for bit).
//   * NOTHING IS ALLOCATED while a glide runs.

#include <felitronics_test.h>
#include <alloc_counter.h>   // installs the allocation counter: EVERY form of `new`, over-aligned included
#include <felitronics/saturation/Saturator.h>

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
using Sat = saturation::Saturator;
using Shape = saturation::WaveShaper::Shape;

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

Buf programme (int nch, int n, unsigned seed)
{
    Buf x ((std::size_t) nch, std::vector<float> ((std::size_t) n));
    std::uint32_t s = seed;
    for (int i = 0; i < n; ++i)
        for (int c = 0; c < nch; ++c)
        {
            s = s * 1664525u + 1013904223u;
            const float noise = ((float) (s >> 8) / 16777216.0f - 0.5f) * 0.3f;
            x[(std::size_t) c][(std::size_t) i] = (float) (0.6 * std::sin (0.021 * i + c)) + noise;
        }
    return x;
}

Sat::Params P (Shape shape, float drive, float bias = 0.0f, float mix = 1.0f, float out = 0.0f, float comp = 0.5f)
{
    Sat::Params p; p.shape = shape; p.driveDb = drive; p.bias = bias; p.mix = mix; p.outputDb = out; p.autoComp = comp;
    return p;
}

// A write at a stream position. `nch == 0` for a stretch makes it clock-only (a pause) in the cut runner below.
struct Write { int at; Sat::Params p; };
struct Pause { int from, to; };

// Renders `x` through a fresh stage prepared with (os, maxBlock), `p0` written before the first sample, the writes
// delivered at their positions, and the stream cut by `cut(i)` — the length of the i-th call (0 = a zero-length call).
template <class Cut>
Buf render (const Buf& x, int os, int maxBlock, const Sat::Params& p0, const std::vector<Write>& ws, Cut cut,
            const std::vector<Pause>& pauses = {}, Sat* keep = nullptr)
{
    Buf y = x;
    Sat local;
    Sat& s = keep != nullptr ? *keep : local;
    s.setParams (p0);
    if (! s.prepare (kFs, maxBlock, (int) x.size(), os)) return {};
    const int n = (int) x[0].size(), nch = (int) x.size();
    std::size_t w = 0;
    int k = 0;
    for (int off = 0; off < n; ++k)
    {
        while (w < ws.size() && ws[w].at <= off) s.setParams (ws[w++].p);
        int m = std::min (cut (k), n - off);
        if (w < ws.size()) m = std::min (m, ws[w].at - off);
        bool paused = false;
        for (const Pause& q : pauses)
        {
            if (off >= q.from && off < q.to) { paused = true; m = std::min (m, q.to - off); }
            else if (off < q.from) m = std::min (m, q.from - off);
        }
        float* io[core::kMaxChannels] {};
        for (int c = 0; c < nch; ++c) io[c] = y[(std::size_t) c].data() + off;
        if (m == 0) { felitronics::test::run (s.process (io, nch, 0)); continue; }
        if (paused)
        {
            felitronics::test::run (s.process (nullptr, 0, m));
            for (int c = 0; c < nch; ++c) std::fill (io[c], io[c] + m, 0.0f);      // nothing was emitted
        }
        else felitronics::test::run (s.process (io, nch, m));
        off += m;
    }
    return y;
}

auto fixed (int b) { return [b] (int) { return b; }; }
auto ragged (std::uint32_t seed)
{
    return [seed] (int k) mutable
    {
        std::uint32_t s = seed + 2654435761u * (std::uint32_t) (k + 1);
        s ^= s >> 13; s *= 0x5bd1e995u; s ^= s >> 15;
        return (s % 7u == 0u) ? 0 : 1 + (int) (s % 700u);                      // zero-length calls included
    };
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
static void testFirstWriteSnaps()
{
    group ("the first write of a stream snaps — before the first sample, after prepare() and after reset()");
    const Buf x = programme (2, 6000, 11u);
    for (int os : { 1, 4 })
    {
        const std::string tag = " (os " + std::to_string (os) + ")";
        // b keeps a's SHAPE: a shape change snaps on its own, and would prove nothing about the first-write rule
        // (the code-review round found the first version of this test passing with that rule removed).
        const Sat::Params a = P (Shape::Tanh, 3.0f), b = P (Shape::Tanh, 12.0f, 0.2f, 0.7f, -2.0f, 0.9f);
        // The reference: b written before prepare(), the only order that never involved a glide.
        const Buf ref = render (x, os, 512, b, {}, fixed (512));
        // prepare(a) -> setParams(b) -> process: the write lands at once.
        Buf y = x;
        Sat s; s.setParams (a);
        ok (s.prepare (kFs, 512, 2, os), "PRECONDITION: prepare" + tag);
        s.setParams (P (Shape::Tanh, 20.0f));                  // several writes before the first sample: all snap
        ok (! s.isGliding(), "a write before the first sample starts no glide" + tag);
        s.setParams (b);
        ok (! s.isGliding(), "no glide is pending after writes before the first sample" + tag);
        float* io[2] { y[0].data(), y[1].data() };
        felitronics::test::run (s.process (io, 2, 6000));
        ok (diffs (ref, y) == 0, "prepare -> setParams -> process renders the written set from sample 0" + tag);

        // A reset() that interrupts a glide: nothing of it resumes, and the next write snaps.
        Sat r; r.setParams (a);
        ok (r.prepare (kFs, 512, 2, os), "PRECONDITION: prepare" + tag);
        Buf z = x;
        float* iz[2] { z[0].data(), z[1].data() };
        felitronics::test::run (r.process (iz, 2, 700));
        r.setParams (P (Shape::Tanh, 18.0f, 0.0f, 0.4f));
        felitronics::test::run (r.process (iz, 2, 300));
        ok (r.isGliding(), "PRECONDITION: a glide is in flight at the reset" + tag);
        r.reset();
        ok (! r.isGliding(), "reset() lands the glide" + tag);
        r.setParams (b);
        Buf w = x;
        float* iw[2] { w[0].data(), w[1].data() };
        felitronics::test::run (r.process (iw, 2, 6000));
        ok (diffs (ref, w) == 0, "after reset() the next write snaps: identical to a fresh stage prepared with it" + tag);
    }
}

static void testGlideIsClockedBySamples()
{
    group ("law 8a — a timeline of writes renders the same bits under every cut, clock-only pauses included");
    const Buf x = programme (2, 20000, 23u);
    const std::vector<Write> ws {
        { 1111, P (Shape::Tanh, 9.0f, 0.0f, 0.8f, -1.5f, 0.3f) },
        { 1500, P (Shape::Tanh, 4.0f, 0.0f, 0.5f, 1.0f, 0.9f) },            // a retarget mid-glide
        { 5003, P (Shape::Asym, 10.0f, 0.25f, 0.9f, 0.0f, 0.5f) },          // a shape switch: snaps
        { 7777, P (Shape::Asym, 14.0f, -0.3f, 0.6f, -3.0f, 1.0f) },         // a glide on the Asym curve, bias moving
        { 12345, P (Shape::Cubic, 6.0f) },
        { 13001, P (Shape::Cubic, 20.0f, 0.0f, 0.3f) },
    };
    const std::vector<Pause> gap { { 8100, 9300 } };                          // a pause in the middle of a glide
    for (int os : { 1, 4 })
    {
        const std::string tag = " (os " + std::to_string (os) + ")";
        const Sat::Params p0 = P (Shape::Tanh, 3.0f);
        const Buf ref = render (x, os, 20000, p0, ws, fixed (20000), gap);
        ok (! ref.empty() && diffs (ref, x) > 10000, "PRECONDITION: the reference renders and the stage works" + tag);
        const Buf still = render (x, os, 20000, p0, {}, fixed (20000), gap);
        ok (diffs (ref, still, 1200) > 1000, "PRECONDITION: the writes are audible in the reference" + tag);
        struct Row { const char* name; int maxBlock; std::function<int (int)> cut; };
        const Row rows[] = {
            { "one sample at a time", 20000, fixed (1) }, { "blocks of 64", 20000, fixed (64) },
            { "blocks of 100", 20000, fixed (100) },     { "ragged with zero-length calls #1", 20000, ragged (1) },
            { "ragged #2, maxBlock 37 (the internal chunker)", 37, ragged (2) },
            { "blocks of 4096, maxBlock 512", 512, fixed (4096) },
        };
        for (const Row& r : rows)
        {
            const Buf got = render (x, os, r.maxBlock, p0, ws, r.cut, gap);
            const long long d = got.empty() ? -1 : diffs (ref, got);
            ok (d == 0, std::string ("bit-identical: ") + r.name + tag + " (" + std::to_string (d) + " differ)");
        }
    }
}

static void testGlideLands()
{
    group ("it lands — past its last period and the downsampler's memory, the stage that always had the target");
    const Buf x = programme (2, 16000, 5u);
    for (int os : { 1, 4 })
    {
        const std::string tag = " (os " + std::to_string (os) + ")";
        Sat probe;
        ok (probe.prepare (kFs, 512, 2, os), "PRECONDITION: prepare");
        const int ticks = probe.glideTicks();
        ok (ticks == 23, "30 ms at 48 kHz is 23 grid periods (" + std::to_string (ticks) + ")");
        const Sat::Params a = P (Shape::Tanh, 3.0f, 0.0f, 1.0f), b = P (Shape::Tanh, 11.0f, 0.0f, 0.6f, -2.0f, 0.8f);
        const int at = 3000;                                        // 3000 = 46.875 periods: mid-period
        const Buf moved  = render (x, os, 512, a, { { at, b } }, fixed (256));
        const Buf always = render (x, os, 512, b, {}, fixed (256));
        const Buf before = render (x, os, 512, a, {}, fixed (256));
        const int start  = (at / 64 + 1) * 64;                      // the next grid boundary
        const int landed = start + ticks * 64 + probe.latencySamples() + 1;
        ok (diffs (moved, before, 0) == diffs (moved, before, start), "nothing moves before the next grid boundary" + tag);
        ok (diffs (moved, always, landed) == 0, "from the landing plus the round trip on, bit-identical to the target stage" + tag
                                                 + " (" + std::to_string (diffs (moved, always, landed)) + " differ)");
        ok (diffs (moved, always, start + ticks * 64 - 64) > 0, "PRECONDITION: and not before it — the glide took its time" + tag);
    }
}

static void testNoClick()
{
    group ("it does not click — the worst second difference inside a glide stays at the steady tone's own");
    const int n = 24000;
    Buf x (1, std::vector<float> ((std::size_t) n));
    for (int i = 0; i < n; ++i) x[0][(std::size_t) i] = (float) (0.5 * std::sin (2.0 * 3.141592653589793 * 227.3 * i / kFs));
    struct Case { const char* name; Sat::Params a, b; };
    const Case cases[] = {
        { "drive 3 -> 9",       P (Shape::Tanh, 3.0f),                  P (Shape::Tanh, 9.0f) },
        { "drive 3 -> 3.5",     P (Shape::Tanh, 3.0f),                  P (Shape::Tanh, 3.5f) },
        { "mix 1 -> 0.5",       P (Shape::Tanh, 9.0f),                  P (Shape::Tanh, 9.0f, 0.0f, 0.5f) },
        { "trim 0 -> -3 dB",    P (Shape::Tanh, 3.0f),                  P (Shape::Tanh, 3.0f, 0.0f, 1.0f, -3.0f) },
        { "autoComp 0.5 -> 1",  P (Shape::Tanh, 9.0f),                  P (Shape::Tanh, 9.0f, 0.0f, 1.0f, 0.0f, 1.0f) },
        { "Asym bias 0 -> 0.3", P (Shape::Asym, 9.0f),                  P (Shape::Asym, 9.0f, 0.3f) },
    };
    for (const Case& c : cases)
    {
        const int at = 12000 + 17;
        const Buf glide = render (x, 4, 512, c.a, { { at, c.b } }, fixed (128));
        Sat stepped;
        stepped.setGlideMs (0.0);
        const Buf step = render (x, 4, 512, c.a, { { at, c.b } }, fixed (128), {}, &stepped);
        // The steady tone's own Δ² on EITHER side of the move: a harder curve has a sharper waveform, so the
        // new setting's steady value is part of the honest bound.
        const double steady = std::max (maxD2 (glide[0], 4000, 10000), maxD2 (glide[0], at + 6000, at + 11000));
        const double g = maxD2 (glide[0], at, at + 4000), s = maxD2 (step[0], at, at + 4000);
        std::printf ("    %-22s steady %.2e  glide %.2e  step %.2e\n", c.name, steady, g, s);
        ok (g < 1.6 * steady, std::string ("the glide stays within 1.6x the steady tone's own Δ² — ") + c.name);
        ok (s > 1.2 * g, std::string ("PRECONDITION: the same move as a step is larger — ") + c.name);
    }
    // A retarget mid-glide, twice, back and forth: still no edge.
    const Buf re = render (x, 4, 512, P (Shape::Tanh, 3.0f),
                           { { 12017, P (Shape::Tanh, 12.0f) }, { 12400, P (Shape::Tanh, 1.0f, 0.0f, 0.6f) }, { 12900, P (Shape::Tanh, 8.0f) } },
                           fixed (128));
    ok (maxD2 (re[0], 12017, 16000) < 1.6 * std::max (maxD2 (re[0], 4000, 10000), maxD2 (re[0], 18000, 23000)),
        "retargets mid-glide do not step either");
}

// NO SPIKE ON THE WAY: a constant input through a drive glide stays inside the range the SETTLED stage covers along
// the same path. The first version interpolated the curve's peak normaliser directly — it goes like 1/k at small
// drives — and a 0 -> 3 dB glide on a constant 0.2 peaked at 7.665 (found by the code-review round). The bound is the
// envelope of 65 settled designs between the two ends, because the path itself can rise above both ends (the auto-
// compensated level of a constant is not monotonic in drive). Asym runs with its DC blocker off here (dcBlockHz 0),
// so the curve is measured and not the blocker's answer to a moving DC; its drive stops at 12 dB, past which a 0.2
// bias flattens the curve at zero and auto-compensation divides by that flatness in the settled stage too.
static void testNoSpikeFromZeroDrive()
{
    group ("no spike on the way — a constant through a drive glide stays inside the settled path's own range");
    const int n = 9000, at = 4000;
    for (Shape sh : { Shape::Tanh, Shape::Atan, Shape::Cubic, Shape::Asym })
        for (float to : { 3.0f, 12.0f, 24.0f })
            for (int os : { 1, 4 })
            {
                if (sh == Shape::Asym && to > 12.0f) continue;
                const float bias = sh == Shape::Asym ? 0.2f : 0.0f;
                auto prm = [&] (float db) { Sat::Params q = P (sh, db, bias); q.dcBlockHz = 0.0f; return q; };
                Buf x (1, std::vector<float> ((std::size_t) n, 0.2f));
                float lo = 1e9f, hi = -1e9f;
                for (int k = 0; k <= 64; ++k)
                {
                    const Buf ys = render (x, os, 512, prm (to * (float) k / 64.0f), {}, fixed (512));
                    const float v = ys[0][(std::size_t) (n - 1)];
                    lo = std::min (lo, v); hi = std::max (hi, v);
                }
                const Buf y = render (x, os, 512, prm (0.0f), { { at, prm (to) } }, fixed (256));
                float gl = 1e9f, gh = -1e9f;
                for (int i = at; i < n; ++i) { gl = std::min (gl, y[0][(std::size_t) i]); gh = std::max (gh, y[0][(std::size_t) i]); }
                const float slack = 0.02f * std::max (std::fabs (lo), std::fabs (hi));
                ok (gl >= lo - slack && gh <= hi + slack,
                    "shape " + std::to_string ((int) sh) + ", 0 -> " + std::to_string ((int) to) + " dB, os " + std::to_string (os)
                    + ": glide [" + std::to_string (gl) + ", " + std::to_string (gh) + "] inside the settled path's ["
                    + std::to_string (lo) + ", " + std::to_string (hi) + "]");
            }
}

static void testRepeatedWritesCostNothing()
{
    group ("a write that changes nothing is free — no glide, and the same bits as writing once");
    const Buf x = programme (2, 8000, 99u);
    const Sat::Params a = P (Shape::Asym, 7.0f, 0.1f, 0.8f, -1.0f, 0.6f);
    std::vector<Write> every;
    for (int k = 1; k < 60; ++k) every.push_back ({ 128 * k, a });
    Sat s;
    const Buf many = render (x, 4, 512, a, every, fixed (128), {}, &s);
    const Buf once = render (x, 4, 512, a, {}, fixed (128));
    ok (diffs (many, once) == 0 && ! s.isGliding(), "re-sending the same parameters every block renders the same bits, no glide");
    // Asked RIGHT AFTER each write, not thousands of samples later (a glide lasts 1472): a write that changes
    // nothing continuous — the same set, or only dcBlockHz, a coefficient that lands at once — starts no glide.
    Sat d;
    d.setParams (a);
    ok (d.prepare (kFs, 512, 2, 4), "PRECONDITION: prepare");
    Buf y = x;
    float* io[2] { y[0].data(), y[1].data() };
    felitronics::test::run (d.process (io, 2, 1000));
    bool never = true;
    d.setParams (a);                        never = never && ! d.isGliding();
    Sat::Params b = a; b.dcBlockHz = 25.0f;
    d.setParams (b);                        never = never && ! d.isGliding();
    felitronics::test::run (d.process (io, 2, 7));
    d.setParams (b);                        never = never && ! d.isGliding();
    ok (never, "neither an unchanged write nor a dcBlockHz-only write starts a glide, asked at once");
    Sat::Params c = b; c.mix = 0.5f;
    d.setParams (c);
    ok (d.isGliding(), "PRECONDITION: a mix write does start one — the probe can see a glide");
}

static void testNoAllocation()
{
    group ("RT — nothing is allocated while a glide runs, writes and retargets included");
    const int N = 128;
    Sat s;
    s.setParams (P (Shape::Asym, 3.0f, 0.1f));
    ok (s.prepare (kFs, N, 2, 4), "PRECONDITION: prepare");
    std::vector<float> a ((std::size_t) N, 0.2f), b ((std::size_t) N, -0.1f);
    float* io[2] { a.data(), b.data() };
    felitronics::test::run (s.process (io, 2, N));
    const long long before = alloc::count.load();
    for (int k = 0; k < 200; ++k)
    {
        if (k % 17 == 0) s.setParams (P (Shape::Asym, 3.0f + (float) (k % 5) * 3.0f, 0.05f * (float) (k % 3), 0.5f + 0.1f * (float) (k % 4)));
        felitronics::test::run (s.process (io, 2, N));
        if (k == 100) felitronics::test::run (s.process (nullptr, 0, 3000));
    }
    felitronics::test::okNoAlloc (alloc::count.load() == before, "no allocation across 200 gliding blocks and a pause");
}

int main()
{
    std::printf ("felitronics::saturation — the parameter glide\n");
    testFirstWriteSnaps();
    testGlideIsClockedBySamples();
    testGlideLands();
    testNoClick();
    testNoSpikeFromZeroDrive();
    testRepeatedWritesCostNothing();
    testNoAllocation();
    return felitronics::test::report();
}

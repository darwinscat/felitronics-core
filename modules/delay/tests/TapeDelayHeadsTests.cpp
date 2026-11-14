// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// delay::TapeDelay's playback heads, held to their design. Wherever an oracle exists it is computed OUTSIDE the class:
//   · the LINE never reads a head, so a head at a whole-sample distance d plays the one-head stage's wet (the frozen
//     class, TapeDelayFrozen.h) shifted by D - d, byte for byte — at REPEATS 0 and in a dark ringing loop;
//   · a figure of heads is the sum of those shifted copies, each scaled, rounded once and flushed, in index order;
//   · the pan law's bits, the gains' counted ramps, the positions' double glide (followed outside), the bend a moving
//     head makes, a silent head standing on its target, the mono rule as a ramp, the snap after a restart;
//   · re-slicing at absolute event samples and the held layout re-sent every block; a pause moving the heads as
//     silence; back at the default layout, the frozen class's bytes; law 8; the tail; the setter's bounds; a head
//     under one sample; the picture's taps; no allocation;
//   · the wake after a restart, the default layout's fast path decided per call and on every side, the order of
//     the sum, the exact guard, a retarget from where the gain stands, side 1's own ramp, glide and mono hold, and
//     the ramp's pinned double formula.
// Every group says in a comment which defects its gates are there to catch; a check with no such comment is the
// fixture's own sanity (the comparison is not of silence, the reference reads where the oracle assumes). The suite is
// written against delay::kMaxHeads, not a count: a check that needs more heads than a build carries says so and is
// skipped.
#include <felitronics_test.h>
#include <alloc_counter.h>   // installs the allocation counter

#include "TapeDelayFrozen.h"

#include <felitronics/delay/TapeDelay.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <initializer_list>
#include <limits>
#include <string>
#include <vector>

using namespace felitronics;
using felitronics::test::ok;
using felitronics::test::group;
namespace frz = felitronics::delay::frozen;

namespace
{
using Plane  = std::vector<float>;
using Planes = std::vector<Plane>;
using Head   = delay::TapeDelayHead;
using Heads  = delay::TapeDelayHeads;

constexpr int kHeads = delay::kMaxHeads;

float  get  (const Plane& v, int i) { return v[(std::size_t) i]; }
float& cell (Plane& v, int i)       { return v[(std::size_t) i]; }

template <class T>
bool same (T a, T b)
{
    return std::memcmp (&a, &b, sizeof (T)) == 0;
}

std::uint32_t bitsOf (float v)
{
    std::uint32_t u = 0;
    std::memcpy (&u, &v, sizeof u);
    return u;
}

bool isSub (float v) { return std::fpclassify (v) == FP_SUBNORMAL; }

/** `n` samples of `a` from `fromA` against `b` from `fromB`, as bytes. An empty or out-of-range span is a broken
    fixture, never a pass. */
bool sameSpan (const Plane& a, int fromA, const Plane& b, int fromB, int n)
{
    if (n <= 0 || fromA < 0 || fromB < 0 || fromA + n > (int) a.size() || fromB + n > (int) b.size())
        return false;
    return std::memcmp (a.data() + fromA, b.data() + fromB, (std::size_t) n * sizeof (float)) == 0;
}

bool samePlanes (const Planes& a, const Planes& b)
{
    if (a.size() != b.size())
        return false;
    for (std::size_t ch = 0; ch < a.size(); ++ch)
        if (a[ch].size() != b[ch].size() || ! sameSpan (a[ch], 0, b[ch], 0, (int) a[ch].size()))
            return false;
    return true;
}

/** A 64-bit LCG read from its HIGH bits (the low bits of a power-of-two LCG cycle with short periods). */
struct Rng
{
    std::uint64_t s;
    explicit Rng (std::uint64_t seed) : s (seed * 0x9E3779B97F4A7C15ull + 0xD1B54A32D192ED03ull) { u32(); u32(); }
    std::uint32_t u32() { s = s * 6364136223846793005ull + 1442695040888963407ull; return (std::uint32_t) (s >> 32); }
    int  below (int n) { return (int) (((std::uint64_t) u32() * (std::uint64_t) n) >> 32); }
    int  range (int lo, int hi) { return lo + below (hi - lo + 1); }
    bool chance (int percent) { return below (100) < percent; }
    float unit() { return (float) (u32() >> 8) * (1.0f / 16777216.0f); }   // [0, 1)
    float signedUnit() { return 2.0f * unit() - 1.0f; }
};

/** A head layout from a list: the heads past the list are off, list entries past kMaxHeads are dropped. */
Heads layout (std::initializer_list<Head> list)
{
    Heads out {};
    int k = 0;
    for (const Head& h : list)
    {
        if (k < kHeads)
            out[(std::size_t) k] = h;
        ++k;
    }
    return out;
}

Heads defaultHeads() { return delay::TapeDelayParams {}.heads; }

/** False, with a printed note, when this build carries fewer heads than a check needs. A run-time answer on purpose:
    the code behind it stays compiled and reachable for every kMaxHeads. */
bool haveHeads (int n, const char* what)
{
    if (kHeads >= n)
        return true;
    std::printf ("      [skipped: %s needs %d heads, this build carries %d]\n", what, n, kHeads);
    return false;
}

delay::TapeDelayParams knobs (float timeMs, float repeats, float darkHz, float offsetMs, float mix, const Heads& heads)
{
    delay::TapeDelayParams p;
    p.timeMs = timeMs; p.repeats = repeats; p.darkHz = darkHz; p.offsetMs = offsetMs; p.mix = mix; p.heads = heads;
    return p;
}

/** The knobs at the default layout: one head at 1x, level 1, centre. */
delay::TapeDelayParams knobs (float timeMs, float repeats, float darkHz, float offsetMs, float mix)
{
    return knobs (timeMs, repeats, darkHz, offsetMs, mix, defaultHeads());
}

frz::TapeDelayParams frozenKnobs (float timeMs, float repeats, float darkHz, float offsetMs, float mix)
{
    frz::TapeDelayParams p;
    p.timeMs = timeMs; p.repeats = repeats; p.darkHz = darkHz; p.offsetMs = offsetMs; p.mix = mix;
    return p;
}

/** The planes through a stage in blocks of `block`, every call as wide as the planes. */
template <class Stage>
bool renderBlocks (Stage& d, Planes& io, int block, bool feed = true)
{
    const int n = (int) io[0].size(), nch = (int) io.size();
    bool accepted = true;
    for (int i = 0; i < n; i += block)
    {
        const int len = std::min (block, n - i);
        float* p[2] = { io[0].data() + i, nch > 1 ? io[1].data() + i : nullptr };
        accepted = d.process (p, nch, len, feed) && accepted;
    }
    return accepted;
}

/** Noise for `len` samples, then silence; `dualMono` puts the same noise on both lines. */
Planes burst (int n, int len, std::uint64_t seed, float amp, bool dualMono)
{
    Rng r (seed);
    Planes io (2, Plane ((std::size_t) n, 0.0f));
    for (int i = 0; i < len && i < n; ++i)
    {
        cell (io[0], i) = amp * r.signedUnit();
        cell (io[1], i) = dualMono ? get (io[0], i) : amp * r.signedUnit();
    }
    return io;
}

int lineLength (double sr)       { return (int) std::ceil (delay::TapeDelay::kMaxDelaySeconds * sr) + 8; }
int offsetLineLength (double sr) { return (int) std::ceil (delay::TapeDelay::kMaxOffsetSeconds * sr) + 8; }
int rampLength (double sr)       { return (int) std::lround (delay::TapeDelay::kHeadRampSeconds * sr); }

/** Where a head reads, as the design says: position x the glided time IN DOUBLE, rounded to float once, then the
    motor's floor of one sample (and the line's reach). */
float readDistance (double position, double delaySamples, double sr)
{
    const float d  = (float) (position * delaySamples);
    const float hi = (float) (lineLength (sr) - 4);
    return d < 1.0f ? 1.0f : (hi < d ? hi : d);
}

/** One contribution as the design forms it: the tap itself at a gain of exactly 1, else the product rounded once (the
    double product of two floats is exact, so its float is the float product) and flushed under 1e-15. */
float scaled (float tap, float g)
{
    if (same (g, 1.0f))
        return tap;
    float c = (float) ((double) tap * (double) g);
    if (std::fabs (c) < 1e-15f)
        c = 0.0f;
    return c;
}

/** The pan law, outside the class: gL = level (1 - p)/sqrt(1 + p^2), gR = level (1 + p)/sqrt(1 + p^2), in double. */
std::array<float, 2> panGains (float level, float pan)
{
    const double p = (double) pan, lv = (double) level;
    const double r = std::sqrt (1.0 + p * p);
    return { (float) (lv * ((1.0 - p) / r)), (float) (lv * ((1.0 + p) / r)) };
}

const float kSqrt2f = (float) std::sqrt (2.0);   // 0x3fb504f3, the hard-pan gain
const float kNaN    = std::numeric_limits<float>::quiet_NaN();

/** Every head's witnesses: position, both side gains, the shown tap; and the motor and the tail. */
struct Witness
{
    double pos[kHeads] {};
    float  gain[kHeads][2] {};
    float  shown[kHeads] {};
    double motor = 0.0;
    float  tail = 0.0f;
};

Witness witness (const delay::TapeDelay& d)
{
    Witness w;
    for (int k = 0; k < kHeads; ++k)
    {
        w.pos[k]     = d.headPosition (k);
        w.gain[k][0] = d.headGain (k, 0);
        w.gain[k][1] = d.headGain (k, 1);
        w.shown[k]   = d.shownHeadMs (k);
    }
    w.motor = d.headSamples();
    w.tail  = d.tailSeconds();
    return w;
}

bool sameWitness (const Witness& a, const Witness& b)
{
    bool s = same (a.motor, b.motor) && same (a.tail, b.tail);
    for (int k = 0; k < kHeads; ++k)
        s = s && same (a.pos[k], b.pos[k]) && same (a.gain[k][0], b.gain[k][0]) && same (a.gain[k][1], b.gain[k][1])
              && same (a.shown[k], b.shown[k]);
    return s;
}

/** A stage processing one sample at a time, so a test can read the witnesses after every sample. */
struct OneByOne
{
    delay::TapeDelay& d;
    float l = 0.0f, r = 0.0f;
    bool step (float inL, float inR, int nch = 2, bool feed = true)
    {
        l = inL; r = inR;
        float* p[2] = { &l, &r };
        return d.process (nch == 0 ? nullptr : p, nch, 1, feed);
    }
};

//==============================================================================
// 1 · PLACEMENT — the line is independent of the heads, so a head at d is the one-head wet moved by D - d
//==============================================================================
void placement()
{
    group ("placement: a head at a whole-sample d plays the frozen one-head wet shifted by D - d, at REPEATS 0 and ringing");
    constexpr double sr = 48000.0;
    constexpr int D = 14400, L = 2000, N = 4 * D + 1200;
    const float positions[] = { 0.5f, 0.25f, 1.0f / 3.0f, 2.0f / 3.0f };

    bool whole = true;
    for (float p : positions)
    {
        const double d = (double) readDistance ((double) p, (double) D, sr);
        whole = whole && same (d, std::floor (d)) && d > 1.0;
    }
    ok (whole, "1/2, 1/4, 1/3 and 2/3 of 14400 samples are whole-sample distances (float (p) x 14400 rounds onto an integer)");

    for (int pass = 0; pass < 2; ++pass)
    {
        // Catches (pass 1): a loop fed from the heads' sum instead of the 1x tap — the shift no longer holds once the
        // repeats recirculate what a head at d heard.
        const float repeats = pass == 0 ? 0.0f : 0.6f;
        const float dark    = pass == 0 ? 4000.0f : 1500.0f;
        const Planes in = burst (N, L, 11u + (std::uint64_t) pass, 0.8f, false);

        frz::TapeDelay then;
        test::run (then.prepare (sr, 512, 2));
        then.setParams (frozenKnobs (300.0f, repeats, dark, 0.0f, 1.0f));
        Planes wetD = in;
        test::run (renderBlocks (then, wetD, 512));
        ok (same (then.headSamples(), (double) D), "the reference reads at exactly 14400 samples (300 ms at 48 kHz)");

        double ringing = 0.0;
        for (int i = 2 * D + L; i < N; ++i)
            ringing += (double) get (wetD[0], i) * (double) get (wetD[0], i);
        if (pass == 1)
            ok (ringing > 1.0e-3, "the reference loop rings past its second pass (the compare is not of silence)");

        bool all = true;
        for (float p : positions)
        {
            delay::TapeDelay now;
            test::run (now.prepare (sr, 512, 2));
            now.setParams (knobs (300.0f, repeats, dark, 0.0f, 1.0f, layout ({ { p, 1.0f, 0.0f } })));
            Planes out = in;
            test::run (renderBlocks (now, out, 512));
            const int shift = D - (int) readDistance ((double) p, (double) D, sr);
            bool here = true;
            for (int ch = 0; ch < 2; ++ch)
                here = here && sameSpan (out[(std::size_t) ch], L, wetD[(std::size_t) ch], L + shift, N - L - shift);
            if (! here)
                std::printf ("      [position %.9g, repeats %.2f: the shifted copy differs]\n", (double) p, (double) repeats);
            all = all && here;
        }
        // Catches: a head distance off by a sample; a pan law that does not leave the centre at unity.
        ok (all, std::string ("heads at 1/2, 1/4, 1/3, 2/3: out(t) == frozen wet(t + D - d), both lines, byte for byte, REPEATS ")
                     + (pass == 0 ? "0" : "0.6 with DARK at 1.5 kHz"));
    }
}

//==============================================================================
// 2a · STACKED HEADS ADD — two, then every head, at 1x, level 1, centre: each plays the loop's own tap, and the
//      heads play their sum in index order. The default layout's fast path (one 1x head at unity IS the tap) must
//      not take two such heads for one.
//==============================================================================
void stackedHeads()
{
    group ("stacked heads add: two, then every head, at 1x, level 1, centre — the frozen one-head wet summed, byte for byte");
    if (! haveHeads (2, "stacked heads"))
        return;

    constexpr double sr = 48000.0;
    constexpr int D = 14400, L = 2000, N = 2 * D + 1200;
    const Planes in = burst (N, L, 31u, 0.8f, false);

    frz::TapeDelay then;
    test::run (then.prepare (sr, 512, 2));
    then.setParams (frozenKnobs (300.0f, 0.6f, 1500.0f, 0.0f, 1.0f));
    Planes wetD = in;
    test::run (renderBlocks (then, wetD, 512));

    for (int count : { 2, kHeads })
    {
        Heads h {};
        for (int k = 0; k < count; ++k)
            h[(std::size_t) k] = { 1.0f, 1.0f, 0.0f };

        delay::TapeDelay now;
        test::run (now.prepare (sr, 512, 2));
        now.setParams (knobs (300.0f, 0.6f, 1500.0f, 0.0f, 1.0f, h));
        Planes out = in;
        test::run (renderBlocks (now, out, 512));

        // After the burst the dry is silent and OFFSET is 0, so the output IS what the heads played: each head's
        // contribution is the loop's tap (a gain of exactly 1), summed in index order — computed here from the
        // frozen class's own output.
        bool all = true, heard = false;
        for (int ch = 0; ch < 2; ++ch)
            for (int i = L; i < N; ++i)
            {
                const float t = get (wetD[(std::size_t) ch], i);
                float want = t;
                for (int k = 1; k < count; ++k)
                    want = want + t;
                all = all && same (get (out[(std::size_t) ch], i), want);
                heard = heard || std::fabs (t) > 1.0e-3f;
            }
        // Catches: the fast path taking stacked unity heads for the one default head (they would play one tap).
        ok (all && heard, std::to_string (count) + " heads at 1x, level 1, centre: the frozen wet summed "
                          + std::to_string (count) + " times in index order, both lines, byte for byte");
    }
}

//==============================================================================
// 2 · THE TRIPLET FIGURE — three shifted copies, each scaled, rounded once, flushed, summed in index order
//==============================================================================
void tripletFigure()
{
    group ("the triplet figure: heads at 1/3, 2/3 and 1, levels 0.3, 0.6, 1 — and two of them ramping");
    if (! haveHeads (3, "the triplet figure"))
        return;

    constexpr double sr = 48000.0;
    constexpr int D = 14400, L = 2000, N = 3 * D + 1000;
    const float pos[3] = { 1.0f / 3.0f, 2.0f / 3.0f, 1.0f };
    const float lvl[3] = { 0.3f, 0.6f, 1.0f };
    const Planes in = burst (N, L, 23u, 0.8f, false);

    frz::TapeDelay then;
    test::run (then.prepare (sr, 512, 2));
    then.setParams (frozenKnobs (300.0f, 0.5f, 2500.0f, 0.0f, 1.0f));
    Planes wetD = in;
    test::run (renderBlocks (then, wetD, 512));

    int dist[3];
    for (int k = 0; k < 3; ++k)
        dist[k] = (int) readDistance ((double) pos[k], (double) D, sr);

    // The oracle at sample t: head k's tap is the frozen wet at t + D - d_k (the line is the same line).
    const auto expected = [&] (int ch, int t, const float* g)
    {
        float s = 0.0f;
        bool any = false;
        for (int k = 0; k < 3; ++k)
        {
            if (! (g[k] > 0.0f))
                continue;
            const float c = scaled (get (wetD[(std::size_t) ch], t + D - dist[k]), g[k]);
            s = any ? s + c : c;
            any = true;
        }
        return 0.0f + s;   // through the offset line at OFFSET 0 and MIX 1 onto a silent dry
    };

    {
        delay::TapeDelay now;
        test::run (now.prepare (sr, 512, 2));
        now.setParams (knobs (300.0f, 0.5f, 2500.0f, 0.0f, 1.0f,
                              layout ({ { pos[0], lvl[0], 0.0f }, { pos[1], lvl[1], 0.0f }, { pos[2], lvl[2], 0.0f } })));
        Planes out = in;
        test::run (renderBlocks (now, out, 512));

        bool gains = true;
        for (int k = 0; k < 3; ++k)
            gains = gains && same (now.headGain (k, 0), lvl[k]) && same (now.headGain (k, 1), lvl[k]);
        ok (gains, "at the centre each side gain is the level itself, exactly");

        int bad = 0;
        for (int ch = 0; ch < 2; ++ch)
            for (int t = L; t < N - (D - dist[0]); ++t)
                bad += same (get (out[(std::size_t) ch], t), expected (ch, t, lvl)) ? 0 : 1;
        // Catches: contributions not flushed / not rounded once; the exact-unity path for the 1x head; a head read at
        // the wrong place. (The burst is shorter than the heads' spacing, so at most one copy is nonzero at a sample:
        // the ORDER of the sum is orderedSum()'s to catch.)
        ok (bad == 0, "wet = flush(tap(1/3) x 0.3) + flush(tap(2/3) x 0.6) + tap(1), in index order, byte for byte ("
                      + std::to_string (bad) + " samples differ)");

        double third = 0.0, first = 0.0;   // the figure as designed: the third echo is the loudest
        for (int i = 0; i < L; ++i)
        {
            first = std::max (first, (double) std::fabs (get (out[0], dist[0] + i)));
            third = std::max (third, (double) std::fabs (get (out[0], dist[2] + i)));
        }
        ok (third > first, "the third echo is louder than the first");
    }

    {
        // Ramping gains on the audio: heads 0 and 1 retargeted mid-figure, the gains read after every sample. Catches a
        // fading head whose read is skipped or mislocated while its gain is still on its way down.
        delay::TapeDelay now;
        test::run (now.prepare (sr, 512, 2));
        now.setParams (knobs (300.0f, 0.5f, 2500.0f, 0.0f, 1.0f,
                              layout ({ { pos[0], lvl[0], 0.0f }, { pos[1], lvl[1], 0.0f }, { pos[2], lvl[2], 0.0f } })));
        OneByOne s { now };
        int bad = 0, ramped = 0;
        const int tF = D + 500;
        for (int t = 0; t < N - (D - dist[0]); ++t)
        {
            if (t == tF)
            {
                now.setHead (0, { pos[0], 0.0f, 0.0f });
                now.setHead (1, { pos[1], 0.9f, 0.0f });
            }
            test::run (s.step (get (in[0], t), get (in[1], t)));
            const float gl[3] = { now.headGain (0, 0), now.headGain (1, 0), now.headGain (2, 0) };
            const float gr[3] = { now.headGain (0, 1), now.headGain (1, 1), now.headGain (2, 1) };
            if (t > tF && t < tF + rampLength (sr))
                ++ramped;
            if (t >= L)
                bad += (same (s.l, expected (0, t, gl)) && same (s.r, expected (1, t, gr))) ? 0 : 1;
        }
        ok (ramped > 900 && bad == 0, "with two gains ramping, every sample is the figure at the gains it heard ("
                                       + std::to_string (bad) + " differ)");
        ok (same (now.headGain (0, 0), 0.0f) && same (now.headGain (1, 1), 0.9f), "...and the ramps landed: 0 and 0.9 exactly");
    }
}

//==============================================================================
// 3 · PING-PONG ON DUAL MONO
//==============================================================================
void pingPong()
{
    group ("ping-pong on dual mono: 1/2 hard left, 1 hard right — each side plays its own head at sqrt 2, the far side nothing");
    if (! haveHeads (2, "the ping-pong"))
        return;

    constexpr double sr = 48000.0;
    constexpr int D = 14400, L = 2000, N = 3 * D;
    const Planes in = burst (N, L, 31u, 0.7f, true);

    frz::TapeDelay then;
    test::run (then.prepare (sr, 512, 2));
    then.setParams (frozenKnobs (300.0f, 0.4f, 4000.0f, 0.0f, 1.0f));
    Planes wetD = in;
    test::run (renderBlocks (then, wetD, 512));
    ok (samePlanes ({ wetD[0] }, { wetD[1] }), "dual mono: the reference's two lines are the same bytes");

    delay::TapeDelay now;
    test::run (now.prepare (sr, 512, 2));
    now.setParams (knobs (300.0f, 0.4f, 4000.0f, 0.0f, 1.0f, layout ({ { 0.5f, 1.0f, -1.0f }, { 1.0f, 1.0f, 1.0f } })));
    Planes out = in;
    test::run (renderBlocks (now, out, 512));

    ok (bitsOf (kSqrt2f) == 0x3fb504f3u, "the hard-pan gain is float sqrt 2, 0x3fb504f3");
    ok (bitsOf (now.headGain (0, 0)) == 0x3fb504f3u && bitsOf (now.headGain (0, 1)) == 0u
        && bitsOf (now.headGain (1, 0)) == 0u && bitsOf (now.headGain (1, 1)) == 0x3fb504f3u,
        "head 0 is (sqrt 2, +0), head 1 is (+0, sqrt 2), bit for bit");

    int bad = 0;
    for (int t = L; t < N - D / 2; ++t)
    {
        bad += same (get (out[0], t), 0.0f + scaled (get (wetD[0], t + D / 2), kSqrt2f)) ? 0 : 1;   // left: D/2 only
        bad += same (get (out[1], t), 0.0f + scaled (get (wetD[0], t), kSqrt2f)) ? 0 : 1;           // right: D only
    }
    // Catches: a far side that is not exactly silent; the pan law's ends; a mono-style read of side 0 on a stereo call.
    ok (bad == 0, "left = sqrt2 x line(t - D/2), right = sqrt2 x line(t - D), byte for byte (" + std::to_string (bad) + " differ)");

    // Mirrored: the 1/2 head HARD RIGHT, so its side 0 is silent and only side 1 sounds. Catches a read geometry
    // worked out only for a head whose side 0 sounds (the right line would read a stale or zero place).
    delay::TapeDelay mirror;
    test::run (mirror.prepare (sr, 512, 2));
    mirror.setParams (knobs (300.0f, 0.4f, 4000.0f, 0.0f, 1.0f, layout ({ { 0.5f, 1.0f, 1.0f }, { 1.0f, 1.0f, -1.0f } })));
    Planes outM = in;
    test::run (renderBlocks (mirror, outM, 512));
    int badM = 0;
    for (int t = L; t < N - D / 2; ++t)
    {
        badM += same (get (outM[1], t), 0.0f + scaled (get (wetD[0], t + D / 2), kSqrt2f)) ? 0 : 1;   // right: D/2 only
        badM += same (get (outM[0], t), 0.0f + scaled (get (wetD[0], t), kSqrt2f)) ? 0 : 1;           // left: D only
    }
    ok (badM == 0, "mirrored, 1/2 hard right and 1 hard left: right = sqrt2 x line(t - D/2), left = sqrt2 x line(t - D) ("
                       + std::to_string (badM) + " differ)");
}

//==============================================================================
// 4 · THE PAN LAW'S BITS
//==============================================================================
void panLaw()
{
    group ("the pan law: centre exactly unity, hard pan (sqrt 2, 0), constant power, mirror-exact, the rational law's bits");
    delay::TapeDelay d;
    test::run (d.prepare (1000.0, 64, 2));   // the slowest rate: a restart per pan clears the smallest line
    OneByOne s { d };
    const auto landed = [&] (float level, float pan)
    {
        d.setHead (0, { 1.0f, level, pan });
        d.reset();
        test::run (s.step (0.0f, 0.0f));   // the first carried call after a restart puts the gains on their targets
        return std::array<float, 2> { d.headGain (0, 0), d.headGain (0, 1) };
    };

    const auto centre = landed (1.0f, 0.0f), negCentre = landed (1.0f, -0.0f);
    ok (same (centre[0], 1.0f) && same (centre[1], 1.0f) && same (negCentre[0], 1.0f) && same (negCentre[1], 1.0f),
        "centre (and -0) is exactly 1.0f on both sides");
    const auto left = landed (1.0f, -1.0f), right = landed (1.0f, 1.0f);
    ok (bitsOf (left[0]) == 0x3fb504f3u && bitsOf (left[1]) == 0u && bitsOf (right[0]) == 0u && bitsOf (right[1]) == 0x3fb504f3u,
        "pan -1 is (0x3fb504f3, +0), +1 its mirror");

    std::vector<float> pans;
    for (int i = 0; i <= 4000; ++i)
        pans.push_back ((float) (i - 2000) / 2000.0f);
    for (float extra : { std::nextafter (1.0f, 0.0f), std::nextafter (-1.0f, 0.0f), 1.0e-30f, -1.0e-30f,
                         std::numeric_limits<float>::min(), -std::numeric_limits<float>::min(), 0.7071068f, -0.7071068f })
        pans.push_back (extra);

    double worstPower = 0.0;
    bool mirror = true, rational = true, monotone = true;
    float lastL = 3.0f, lastR = -1.0f;
    for (std::size_t i = 0; i < pans.size(); ++i)
    {
        const float p = pans[i];
        const auto g = landed (1.0f, p), m = landed (1.0f, -p), want = panGains (1.0f, p);
        worstPower = std::max (worstPower, std::fabs ((double) g[0] * g[0] + (double) g[1] * g[1] - 2.0));
        mirror   = mirror && same (g[0], m[1]) && same (g[1], m[0]);
        rational = rational && same (g[0], want[0]) && same (g[1], want[1]);
        if (i <= 4000)
        {
            monotone = monotone && ! (lastL < g[0]) && ! (g[1] < lastR);
            lastL = g[0]; lastR = g[1];
        }
    }
    // Catches: a -3 dB centre (gains / sqrt 2); a sine law or a linear law in place of the rational one.
    ok (worstPower <= 2.0e-7, "gL^2 + gR^2 within 2e-7 of 2 across 4001 pans and the edges (worst " + std::to_string (worstPower) + ")");
    ok (mirror, "mirror-exact: gains(-p) are gains(p) swapped, bit for bit, at every pan swept");
    ok (rational, "every pan's gains are the rational law's, computed outside in double, bit for bit");
    ok (monotone, "the left gain never rises and the right never falls as the pan moves right");

    bool levels = true;
    for (float lv : { 0.37f, 0.001f, 1.0e-12f })
        for (float p : { -0.8f, -0.1f, 0.0f, 0.45f, 1.0f })
        {
            const auto g = landed (lv, p), want = panGains (lv, p);
            levels = levels && same (g[0], want[0]) && same (g[1], want[1]);
        }
    ok (levels, "level x the pan law: the double product, rounded once, at levels 0.37, 0.001, 1e-12");

    // And the same bits when the gains get there by a ramp rather than a snap, at 48 kHz.
    delay::TapeDelay r48;
    test::run (r48.prepare (48000.0, 64, 2));
    OneByOne s48 { r48 };
    test::run (s48.step (0.0f, 0.0f));
    r48.setHead (0, { 1.0f, 1.0f, -1.0f });
    for (int i = 0; i < rampLength (48000.0); ++i)
        test::run (s48.step (0.0f, 0.0f));
    ok (bitsOf (r48.headGain (0, 0)) == 0x3fb504f3u && bitsOf (r48.headGain (0, 1)) == 0u,
        "ramped from the centre to hard left over 960 samples: (0x3fb504f3, +0) exactly");
}

//==============================================================================
// 4b · PAN-LAW ANCHORS — off-centre pans at levels other than 1, against bits computed OUTSIDE the class: the
//      rational law evaluated to 60 digits and rounded to float once (Python's decimal; not the class's
//      double arithmetic). The table is the oracle; a law rounded at a different step would miss it.
//==============================================================================
void panAnchors()
{
    group ("pan-law anchors: five off-centre pans at five levels, the gains' bits against a 60-digit evaluation");
    struct Anchor { float pan, level; std::uint32_t left, right; };
    const Anchor table[] = {
        {  0.5f,  0.5f,  0x3e64f92eu, 0x3f2bbae3u },
        { -0.25f, 0.8f,  0x3f785b43u, 0x3f15038eu },
        {  0.75f, 1.0f,  0x3e4ccccdu, 0x3fb33333u },
        { -0.9f,  0.3f,  0x3ed8ec48u, 0x3cb6ac0au },
        {  0.1f,  0.65f, 0x3f15044bu, 0x3f3621b1u },
    };
    bool all = true;
    for (const Anchor& a : table)
    {
        delay::TapeDelay d;
        test::run (d.prepare (48000.0, 64, 2));
        d.setParams (knobs (300.0f, 0.3f, 4000.0f, 0.0f, 1.0f, layout ({ { 1.0f, a.level, a.pan } })));
        OneByOne s { d };
        test::run (s.step (0.0f, 0.0f));   // the first carried call puts the gains on target
        const bool here = bitsOf (d.headGain (0, 0)) == a.left && bitsOf (d.headGain (0, 1)) == a.right;
        if (! here)
            std::printf ("      [pan %g level %g: got 0x%08x 0x%08x]\n", (double) a.pan, (double) a.level,
                         (unsigned) bitsOf (d.headGain (0, 0)), (unsigned) bitsOf (d.headGain (0, 1)));
        all = all && here;
    }
    // Catches: the pan law rounded to float before the level multiplies it, a different law with the same centre
    // and ends, the level applied in float.
    ok (all, "every anchor's two gains, bit for bit");
}

//==============================================================================
// 5 · RAMPS
//==============================================================================
void ramps()
{
    group ("ramps: a level lands in exactly round(0.02 sr) samples on its exact target; re-sends and unmoved sides keep it");
    for (double sr : { 48000.0, 44100.0, 1000.0 })
    {
        const int len = rampLength (sr);
        delay::TapeDelay d;
        test::run (d.prepare (sr, 64, 2));
        const auto p = knobs (300.0f, 0.3f, 4000.0f, 0.0f, 1.0f, layout ({ { 1.0f, 0.2f, 0.0f } }));
        d.setParams (p);
        OneByOne s { d };
        for (int i = 0; i < 10; ++i)
            test::run (s.step (0.0f, 0.0f));
        d.setHead (0, { 1.0f, 0.9f, 0.0f });
        const auto held = knobs (300.0f, 0.3f, 4000.0f, 0.0f, 1.0f, layout ({ { 1.0f, 0.9f, 0.0f } }));

        int landedAt = -1;
        bool rising = true, linear = true, stays = true;
        float last = 0.2f;
        for (int m = 1; m <= len + 40; ++m)
        {
            if (m == len / 2)
            {
                d.setHead (0, { 1.0f, 0.9f, 0.0f });   // re-sent mid-ramp: must not restart it
                d.setParams (held);
                d.setHeads (held.heads);
            }
            test::run (s.step (0.0f, 0.0f));
            const float g = d.headGain (0, 0);
            if (landedAt < 0 && same (g, 0.9f))
                landedAt = m;
            else if (landedAt >= 0 && ! same (g, 0.9f))
                stays = false;
            if (landedAt < 0)
            {
                rising = rising && last < g;
                const double want = (double) 0.2f + ((double) 0.9f - (double) 0.2f) * (double) m / (double) len;
                linear = linear && std::fabs ((double) g - want) <= 1.0e-6;
            }
            last = g;
        }
        const std::string at = " at " + std::to_string ((int) sr) + " Hz";
        // Catches: a ramp of the wrong length, one that lands a sample early or late, or parks short of its target; a
        // re-send that restarts it.
        ok (landedAt == len && stays, "0.2 -> 0.9 lands on 0.9f exactly at sample " + std::to_string (len) + ", not one before (landed at "
                                      + std::to_string (landedAt) + ")" + at);
        ok (rising && linear, "...rising every sample, on the linear ramp to 1e-6" + at);
    }

    {
        // A retarget that moves one side must leave the other's ramp alone: hard left, then a level change mid-ramp.
        // Catches a retarget without its exact-value guard (the unmoved side would restart and land late).
        constexpr double sr = 48000.0;
        const int len = rampLength (sr);
        delay::TapeDelay d;
        test::run (d.prepare (sr, 64, 2));
        d.setParams (knobs (300.0f, 0.3f, 4000.0f, 0.0f, 1.0f));
        OneByOne s { d };
        test::run (s.step (0.0f, 0.0f));
        d.setHead (0, { 1.0f, 1.0f, -1.0f });
        int rightAt = -1, leftAt = -1;
        const float leftTo = panGains (0.8f, -1.0f)[0];
        for (int m = 1; m <= 2 * len; ++m)
        {
            if (m == 400)
                d.setHead (0, { 1.0f, 0.8f, -1.0f });   // the left gain moves; the right's target stays 0
            test::run (s.step (0.0f, 0.0f));
            if (rightAt < 0 && same (d.headGain (0, 1), 0.0f)) rightAt = m;
            if (leftAt < 0 && same (d.headGain (0, 0), leftTo)) leftAt = m;
        }
        ok (rightAt == len && leftAt == 399 + len,
            "a level change mid-pan: the right side still lands on 0 at 960, the left on its new target at 400 + 959 (right "
                + std::to_string (rightAt) + ", left " + std::to_string (leftAt) + ")");
    }

    {
        // In a mono call side 0's target is the level: a pan move must not restart its ramp. Same defect, other door.
        constexpr double sr = 48000.0;
        const int len = rampLength (sr);
        delay::TapeDelay d;
        test::run (d.prepare (sr, 64, 2));
        d.setParams (knobs (300.0f, 0.3f, 4000.0f, 0.0f, 1.0f, layout ({ { 1.0f, 0.5f, 0.3f } })));
        OneByOne s { d };
        for (int i = 0; i < 2 * len; ++i)
            test::run (s.step (0.0f, 0.0f, 1));
        ok (same (d.headGain (0, 0), 0.5f), "a mono stage plays side 0 at the head's level");
        d.setHead (0, { 1.0f, 0.7f, 0.3f });
        int at = -1;
        for (int m = 1; m <= 2 * len; ++m)
        {
            if (m == 300)
                d.setHead (0, { 1.0f, 0.7f, -0.4f });   // pan only: side 0's mono target is unchanged
            test::run (s.step (0.0f, 0.0f, 1));
            if (at < 0 && same (d.headGain (0, 0), 0.7f)) at = m;
        }
        ok (at == len, "mono: a pan move mid-ramp leaves the level ramp landing at 960 (landed at " + std::to_string (at) + ")");
    }

    {
        // A tiny level just above the flush, a pan one float short of hard: every gain normal at every sample, both
        // sides above zero, landing on the rational law's bits; and back down to exactly 0 with nothing subnormal.
        constexpr double sr = 8000.0;
        const int len = rampLength (sr);
        const float tiny = std::nextafter (1.0e-15f, 1.0f);
        delay::TapeDelay d;
        test::run (d.prepare (sr, 64, 2));
        d.setParams (knobs (100.0f, 0.5f, 4000.0f, 3.0f, 1.0f, layout ({ { 1.0f, 0.5f, 0.2f } })));
        OneByOne s { d };
        Rng r (5u);
        for (int i = 0; i < 2000; ++i)
            test::run (s.step (0.01f * r.signedUnit(), 0.01f * r.signedUnit()));
        bool normal = true, landed = true, clean = true;
        for (float pan : { std::nextafter (1.0f, 0.0f), std::nextafter (-1.0f, 0.0f) })
        {
            d.setHead (0, { 1.0f, tiny, pan });
            for (int m = 0; m < len + 4; ++m)
            {
                test::run (s.step (0.01f * r.signedUnit(), 0.01f * r.signedUnit()));
                normal = normal && ! isSub (d.headGain (0, 0)) && ! isSub (d.headGain (0, 1));
                clean = clean && ! d.anySubnormalState();
            }
            const auto want = panGains (tiny, pan);
            landed = landed && same (d.headGain (0, 0), want[0]) && same (d.headGain (0, 1), want[1])
                            && want[0] > 0.0f && want[1] > 0.0f;
        }
        d.setHead (0, { 1.0f, 0.0f, 0.0f });
        for (int m = 0; m < len + 4; ++m)
        {
            test::run (s.step (0.01f * r.signedUnit(), 0.01f * r.signedUnit()));
            normal = normal && ! isSub (d.headGain (0, 0)) && ! isSub (d.headGain (0, 1));
            clean = clean && ! d.anySubnormalState();
        }
        ok (normal && clean, "a level of 1e-15 (+1 ulp) one float short of hard pan, and back to 0: no subnormal gain or state at any sample");
        ok (landed, "...the tiny gains land on the rational law's bits, both sides above zero");
        ok (same (d.headGain (0, 0), 0.0f) && same (d.headGain (0, 1), 0.0f), "...and level 0 lands on +0 exactly");
    }
}

//==============================================================================
// 6 · THE POSITION GLIDE
//==============================================================================
void positionGlide()
{
    group ("the position glide: TIME's law in double, followed outside to 1e-12, landing exactly on its target");
    constexpr double sr = 48000.0;
    delay::TapeDelay d;
    test::run (d.prepare (sr, 64, 2));
    d.setParams (knobs (300.0f, 0.3f, 4000.0f, 0.0f, 1.0f, layout ({ { 1.0f, 1.0f, 0.0f }, { 0.25f, 0.5f, 0.6f } })));
    OneByOne s { d };
    for (int i = 0; i < 64; ++i)
        test::run (s.step (0.0f, 0.0f));
    const int moving = kHeads >= 2 ? 2 : 1;
    const float tgt[2] = { 0.5f, 0.75f };
    d.setHead (0, { tgt[0], 1.0f, 0.0f });
    if (moving > 1)
        d.setHead (1, { tgt[1], 0.5f, 0.6f });

    // The design's coefficient for a 150 ms glide, from the system exp: the class uses its own det exp, so the two
    // trajectories may part by a rounding, never by a float's quantisation.
    const double coeff = 1.0 - std::exp (-1.0 / (0.15 * sr));
    double model[2] = { 1.0, 0.25 };
    int landedAt[2] = { -1, -1 };
    double worst = 0.0;
    bool stays = true, between = true;
    for (int m = 1; m <= 300000; ++m)
    {
        test::run (s.step (0.0f, 0.0f));
        for (int k = 0; k < moving; ++k)
        {
            const double t = (double) tgt[k];
            const double next = model[k] + (t - model[k]) * coeff;
            model[k] = std::fabs (next - model[k]) > 0.0 ? next : t;
            const double hp = d.headPosition (k);
            worst = std::max (worst, std::fabs (hp - model[k]));
            if (landedAt[k] < 0 && same (hp, t))
                landedAt[k] = m;
            else if (landedAt[k] >= 0 && ! same (hp, t))
                stays = false;
            if (landedAt[k] < 0)
                between = between && (k == 0 ? (hp < 1.0 && hp > t) : (hp > 0.25 && hp < t));
        }
    }
    // Catches: a position that glides in float (it parks ~1e-4 short and jumps), a glide with another law or rate.
    ok (worst <= 1.0e-12, "every sample within 1e-12 of next = cur + (tgt - cur)(1 - exp(-1/(0.15 sr))) (worst "
                              + std::to_string (worst) + ")");
    bool landed = landedAt[0] > 0 && (moving < 2 || landedAt[1] > 0);
    ok (landed && stays, "the glide lands on (double) target exactly and stays there (head 0 at sample "
                             + std::to_string (landedAt[0]) + ")");
    ok (between, "...from the side it started on, never overshooting");
    ok (landedAt[0] > 150000, "...and not early: a 150 ms law takes seconds to close a double's last bits");
}

//==============================================================================
// 7 · THE BEND
//==============================================================================
void bend()
{
    group ("the bend: a sounding head moved from 1/2 to 1/4 reads closer and closer — its repeat rises in pitch, then settles");
    constexpr double sr = 48000.0;
    constexpr int N = 100000, t0 = 30000, Df = 288000;   // the reference reads at 6 s: it holds the whole line
    Plane x ((std::size_t) N, 0.0f);
    for (int i = 0; i < N; ++i)
        cell (x, i) = 0.5f * (float) std::sin (2.0 * 3.14159265358979323846 * (double) i / 48.0);

    // The line, from outside: at REPEATS 0 the recorded line does not depend on TIME, so the frozen class at a 6 s
    // delay plays it back intact — line(s) is its output at s + 288000.
    Plane line;
    {
        frz::TapeDelay then;
        test::run (then.prepare (sr, 4096, 1));
        then.setParams (frozenKnobs (6000.0f, 0.0f, 4000.0f, 0.0f, 1.0f));
        Planes io { x };
        io[0].resize ((std::size_t) (N + Df), 0.0f);
        test::run (renderBlocks (then, io, 4096));
        line.assign (io[0].begin() + Df, io[0].begin() + Df + N);
    }

    delay::TapeDelay d;
    test::run (d.prepare (sr, 4096, 2));
    d.setParams (knobs (300.0f, 0.0f, 4000.0f, 0.0f, 1.0f, layout ({ { 0.5f, 1.0f, 0.0f } })));
    OneByOne s { d };
    Plane wet ((std::size_t) N, 0.0f);
    double worst = 0.0, motorMid = 0.0;
    for (int t = 0; t < N; ++t)
    {
        if (t == t0)
            d.setHead (0, { 0.25f, 1.0f, 0.0f });
        if (t == t0 + 5000)
            d.setTimeMs (200.0f);   // the motor glides 14400 -> 9600 under the gliding head
        test::run (s.step (get (x, t), get (x, t)));
        cell (wet, t) = (float) ((double) s.l - (double) get (x, t));
        if (t == t0 + 20000 - 1)
            motorMid = d.headSamples();

        // The read, outside: the class's own position and motor after this sample, the design's distance, and a
        // linear read of the outside line. Catches the product taken after rounding the position to float, and a
        // head placed at a fraction of the TARGET time instead of the glided one (TIME moves mid-window).
        if (t >= t0 - 100 && t < t0 + 20000)
        {
            const float hd = readDistance (d.headPosition (0), d.headSamples(), sr);
            const float hpos = (float) t - hd;
            const int   i0 = (int) hpos;
            const double f = (double) (hpos - (float) i0);
            const double want = (double) get (line, i0) * (1.0 - f) + (double) get (line, i0 + 1) * f;
            worst = std::max (worst, std::fabs (((double) s.l - (double) get (x, t)) - want));
        }
    }
    ok (motorMid > 9600.5 && motorMid < 14399.5, "the motor was still gliding at the window's end (the head rode a moving time)");
    ok (worst <= 1.0e-6, "while it and TIME glide the head reads the outside line at (float)(position x glided time), to 1e-6 (worst "
                             + std::to_string (worst) + ")");

    const auto halfPeriod = [&] (int from, int to)
    {
        int crossings = 0;
        for (int t = from + 1; t < to; ++t)
            crossings += (get (wet, t - 1) < 0.0f) != (get (wet, t) < 0.0f) ? 1 : 0;
        return crossings > 0 ? (double) (to - from) / (double) crossings : 0.0;
    };
    const double before = halfPeriod (t0 - 6000, t0), during = halfPeriod (t0, t0 + 2400), after = halfPeriod (N - 6000, N);
    // Catches: a sounding head that snaps (no bend) or one that never arrives.
    ok (before > 23.5 && before < 24.5, "before the move the repeat has the tone's period (half " + std::to_string (before) + ")");
    ok (during > 14.0 && during < 19.0, "while the head closes in, the half period shortens to ~2/3 (" + std::to_string (during) + ")");
    ok (after > 23.5 && after < 24.5, "settled at 1/4, the period is the tone's again (half " + std::to_string (after) + ")");
}

//==============================================================================
// 8 · A SILENT HEAD STANDS ON ITS TARGET; A SOUNDING HEAD GLIDES
//==============================================================================
void silentSnapsSoundingGlides()
{
    group ("a silent head stands on its target, a sounding one glides — decided by the gains it has NOW, not the ones it is sent");
    constexpr double sr = 48000.0;
    const float third = 1.0f / 3.0f;
    {
        delay::TapeDelay d;
        test::run (d.prepare (sr, 64, 2));
        d.setParams (knobs (300.0f, 0.3f, 4000.0f, 0.0f, 1.0f, Heads {}));
        OneByOne s { d };
        for (int i = 0; i < 64; ++i)
            test::run (s.step (0.1f, -0.1f));
        d.setHead (0, { third, 1.0f, 0.0f });   // off at 1, then on at 1/3 at once
        test::run (s.step (0.1f, -0.1f));
        // Catches: a silent head that glides; a snap decided by the TARGET gain (here 1, so it would glide).
        ok (same (d.headPosition (0), (double) third), "an off head turned on at 1/3: on 1/3 exactly after one sample");
        ok (d.headGain (0, 0) > 0.0f && d.headGain (0, 0) < 0.01f, "...and its gain starts its ramp from 0");
    }
    {
        delay::TapeDelay d;
        test::run (d.prepare (sr, 64, 2));
        d.setParams (knobs (300.0f, 0.3f, 4000.0f, 0.0f, 1.0f));
        OneByOne s { d };
        for (int i = 0; i < 64; ++i)
            test::run (s.step (0.1f, -0.1f));
        d.setHead (0, { third, 1.0f, 0.0f });
        test::run (s.step (0.1f, -0.1f));
        const double p = d.headPosition (0);
        ok (p < 1.0 && p > 0.99, "a sounding head sent to 1/3 has taken one glide step after one sample, no more");

        // Off and moved at once: its gains are still on their way down, so it still glides; once they reach 0 it
        // stands on its target at the next sample. Catches a snap decided by the target gain (0 here: it would jump).
        d.setHead (0, { 0.75f, 0.0f, 0.0f });
        const int len = rampLength (sr);
        bool gliding = true;
        for (int m = 1; m <= len; ++m)
        {
            test::run (s.step (0.1f, -0.1f));
            gliding = gliding && ! same (d.headPosition (0), 0.75);
        }
        ok (gliding && same (d.headGain (0, 0), 0.0f) && same (d.headGain (0, 1), 0.0f),
            "turned off and moved: it glides while its ramp runs down and the ramp lands on 0");
        test::run (s.step (0.1f, -0.1f));
        ok (same (d.headPosition (0), 0.75), "...and the sample after, silent, it stands on 0.75 exactly");
    }
}

//==============================================================================
// 9 · THE MONO RULE IS A RAMP
//==============================================================================
void monoIsARamp()
{
    group ("the mono rule is a ramp: a hard-right head on a stereo stage, the call width 2 -> 1 -> 2, DC in");
    constexpr double sr = 48000.0;
    const int len = rampLength (sr);
    delay::TapeDelay d;
    test::run (d.prepare (sr, 4096, 2));
    d.setParams (knobs (300.0f, 0.0f, 4000.0f, 0.0f, 1.0f, layout ({ { 1.0f, 0.45f, 1.0f } })));
    Planes dc (2, Plane (30000, 0.5f));
    test::run (renderBlocks (d, dc, 1000));
    ok (same (get (dc[0], 29999), 0.5f), "stereo: the left side of a hard-right head is silent (the dry, exactly)");

    OneByOne s { d };
    const auto edge = [&] (int nch, float target, float prev, double& worstStep, int& landedAt, float& last)
    {
        worstStep = 0.0; landedAt = -1;
        for (int m = 1; m <= 2 * len; ++m)
        {
            test::run (s.step (0.5f, 0.5f, nch));
            worstStep = std::max (worstStep, std::fabs ((double) s.l - (double) prev));
            prev = s.l;
            if (landedAt < 0 && same (d.headGain (0, 0), target))
                landedAt = m;
        }
        last = s.l;
    };

    double upStep = 0.0, downStep = 0.0;
    int upAt = -1, downAt = -1;
    float monoOut = 0.0f, stereoOut = 0.0f;
    edge (1, 0.45f, get (dc[0], 29999), upStep, upAt, monoOut);
    const double lineDc = ((double) monoOut - 0.5) / 0.45;   // the DC the line holds, read from the landed mono wet
    edge (2, 0.0f, monoOut, downStep, downAt, stereoOut);
    const double bound = 1.05 * lineDc * 0.45 / (double) len + 1.0e-6;
    // Catches: mono as a switch (a 0.45 x line step on the edge, ~0.2 here).
    ok (lineDc > 0.3 && upStep <= bound, "2 -> 1: the left wet's largest step is one ramp step (" + std::to_string (upStep)
                                             + " <= " + std::to_string (bound) + ")");
    ok (upAt == len, "...and side 0 lands on the level, 0.45 exactly, at sample 960 (" + std::to_string (upAt) + ")");
    ok (downStep <= bound && downAt == len && same (stereoOut, 0.5f),
        "1 -> 2: back down by ramp steps, landing on +0 at 960, the left silent again (" + std::to_string (downAt) + ")");
}

//==============================================================================
// 10 · THE GAINS SNAP ON THE FIRST CARRIED CALL AFTER A RESTART, NOT ON A PAUSE
//==============================================================================
void snapAfterRestart()
{
    group ("after reset the gains snap on the first call that carries a channel, not on a pause, whatever the past was");
    constexpr double sr = 48000.0;
    const Heads lay = layout ({ { 1.0f, 0.8f, 0.6f }, { 0.5f, 0.7f, -0.4f } });
    const Planes noise = burst (20000, 20000, 41u, 0.4f, false);

    const auto run = [&] (int past)   // 0: never ran; 1: a mono past; 2: a stereo past
    {
        delay::TapeDelay d;
        test::run (d.prepare (sr, 512, 2));
        d.setParams (knobs (300.0f, 0.5f, 3000.0f, 4.0f, 1.0f, lay));
        if (past > 0)
        {
            Planes io = noise;
            for (int i = 0; i < 3000; i += 256)
            {
                float* p[2] = { io[0].data() + i, io[1].data() + i };
                test::run (d.process (p, past, std::min (256, 3000 - i)));
            }
            d.reset();
        }
        test::run (d.process (nullptr, 0, 500));
        Planes io = noise;
        test::run (renderBlocks (d, io, 512));
        return io;
    };
    const Planes never = run (0), mono = run (1), stereo = run (2);
    // Catches: a reset that leaves the gains aimed (a mono past would ramp side 0 back to its balance after it).
    ok (samePlanes (mono, stereo) && samePlanes (mono, never),
        "reset -> pause -> stereo: the same bytes after a mono past, a stereo past and no past at all");

    // The witness: the first carried call after reset -> pause puts the gains ON target — a target set after the pause
    // included — and a mono first call puts side 0 on the level. Catches gains primed on the pause.
    delay::TapeDelay d;
    test::run (d.prepare (sr, 512, 2));
    d.setParams (knobs (300.0f, 0.5f, 3000.0f, 4.0f, 1.0f, lay));
    OneByOne s { d };
    for (int i = 0; i < 100; ++i)
        test::run (s.step (0.1f, 0.2f));
    d.reset();
    test::run (d.process (nullptr, 0, 300));
    d.setHead (0, { 1.0f, 0.3f, 0.5f });
    test::run (s.step (0.1f, 0.2f));
    const auto want = panGains (0.3f, 0.5f);
    ok (same (d.headGain (0, 0), want[0]) && same (d.headGain (0, 1), want[1]),
        "reset -> pause -> new pan -> one stereo sample: the gains are on the new target, no ramp");

    test::run (s.step (0.1f, 0.2f));
    d.reset();
    test::run (d.process (nullptr, 0, 300));
    test::run (s.step (0.1f, 0.2f, 1));
    ok (same (d.headGain (0, 0), 0.3f), "reset -> pause -> one mono sample: side 0 is on the level, no ramp");
}

//==============================================================================
// 11 · RE-SLICING AT ABSOLUTE EVENT SAMPLES, AND THE HELD STATE RE-SENT EVERY BLOCK
//==============================================================================
enum class Kind { position, level, pan, time, offset, feed, reset };
constexpr int kKinds = 6;        // the kinds drawn in rotation; the restarts are placed by hand
constexpr int kAllKinds = 7;

struct Event
{
    int   at;
    Kind  kind;
    int   head;
    float value;
    bool  partial;   // send only the field that moves, the others not-a-number (ignored)
};

struct Timeline
{
    int n = 96000;
    int monoFrom = 40000, monoTo = 46000, pauseFrom = 60000, pauseTo = 63000;
    std::vector<Event> events;
    Planes input;
    delay::TapeDelayParams start;
};

Timeline makeTimeline()
{
    Timeline tl;
    tl.start = knobs (300.0f, 0.55f, 3000.0f, 7.0f, 0.8f,
                      layout ({ { 1.0f, 1.0f, 0.0f }, { 0.5f, 0.6f, -0.5f }, { 0.25f, 0.4f, 0.7f }, { 0.75f, 0.3f, 0.2f } }));
    Rng r (0x5EEDu);
    for (int i = 0; i < 72; ++i)
    {
        Event e { r.range (1, tl.n - 1), static_cast<Kind> (i % kKinds), r.below (kHeads), 0.0f, r.chance (50) };
        switch (e.kind)
        {
            case Kind::position: e.value = r.chance (20) ? 1.0f : 0.05f + 0.95f * r.unit(); break;
            case Kind::level:    e.value = (i / kKinds) % 3 == 0 ? 0.0f : 0.1f + 0.9f * r.unit(); break;
            case Kind::pan:      e.value = r.signedUnit(); break;
            case Kind::time:     e.value = 80.0f + 500.0f * r.unit(); break;
            case Kind::offset:   e.value = 25.0f * r.signedUnit(); break;
            case Kind::feed:     break;
            case Kind::reset:    break;   // never drawn: i % kKinds stops short of it; the restarts are placed below
        }
        tl.events.push_back (e);
    }
    // Two events on one sample, events on the stretches' edges and inside them.
    tl.events.push_back ({ tl.events[3].at, Kind::pan, 0, -0.9f, false });
    tl.events.push_back ({ tl.monoFrom, Kind::position, 0, 0.6f, true });
    tl.events.push_back ({ tl.monoFrom + 1000, Kind::level, kHeads - 1, 0.9f, false });
    tl.events.push_back ({ tl.pauseFrom + 500, Kind::position, kHeads - 1, 0.3f, false });
    tl.events.push_back ({ tl.pauseFrom + 700, Kind::time, 0, 450.0f, false });
    tl.events.push_back ({ tl.pauseTo, Kind::level, 0, 0.0f, true });
    // Restarts at absolute samples: one in the open, one inside the mono stretch, one inside the pause (the gains
    // then wait for the first call that carries a channel).
    tl.events.push_back ({ 20011, Kind::reset, 0, 0.0f, false });
    tl.events.push_back ({ tl.monoFrom + 2500, Kind::reset, 0, 0.0f, false });
    tl.events.push_back ({ tl.pauseFrom + 1200, Kind::reset, 0, 0.0f, false });
    std::stable_sort (tl.events.begin(), tl.events.end(), [] (const Event& a, const Event& b) { return a.at < b.at; });

    Rng nz (99u);
    tl.input.assign (2, Plane ((std::size_t) tl.n, 0.0f));
    for (int i = 0; i < tl.n; ++i)
    {
        const float tone = 0.3f * (float) std::sin (0.0731 * (double) i);
        cell (tl.input[0], i) = tone + 0.1f * nz.signedUnit();
        cell (tl.input[1], i) = -tone + 0.1f * nz.signedUnit();
    }
    return tl;
}

struct Census
{
    int kinds[kAllKinds] {};
    int inMono = 0, inPause = 0, turnedOn = 0, turnedOff = 0;
    long long monoCalls = 0, pauseCalls = 0, noOpCalls = 0;
    bool noOpsAsSaid = true;   // every zero-length call accepted, every refused one refused
};

void apply (delay::TapeDelay& d, delay::TapeDelayParams& held, bool& feed, const Event& e, Census& c)
{
    Head& h = held.heads[(std::size_t) e.head];
    const float was = h.level;
    switch (e.kind)
    {
        case Kind::position: h.position = e.value; d.setHead (e.head, e.partial ? Head { e.value, kNaN, kNaN } : h); break;
        case Kind::level:    h.level    = e.value; d.setHead (e.head, e.partial ? Head { kNaN, e.value, kNaN } : h); break;
        case Kind::pan:      h.pan      = e.value; d.setHead (e.head, e.partial ? Head { kNaN, kNaN, e.value } : h); break;
        case Kind::time:     held.timeMs = e.value; d.setTimeMs (e.value); break;
        case Kind::offset:   held.offsetMs = e.value; d.setOffsetMs (e.value); break;
        case Kind::feed:     feed = ! feed; break;
        case Kind::reset:    d.reset(); break;
    }
    ++c.kinds[(int) e.kind];
    if (e.kind == Kind::level && ! (was > 0.0f) && e.value > 0.0f) ++c.turnedOn;
    if (e.kind == Kind::level && was > 0.0f && ! (e.value > 0.0f)) ++c.turnedOff;
}

/** Renders the timeline. Every block is cut at its slice edge, at the next event's sample and at the stretches' edges,
    so an event lands on its absolute sample in any slicing. `slice` 0 draws the slices from a seeded generator. */
Planes renderTimeline (const Timeline& tl, int slice, bool heldEveryBlock, Census& c, Witness& w, bool noOps = false)
{
    delay::TapeDelay d;
    test::run (d.prepare (48000.0, 4096, 2));
    delay::TapeDelayParams held = tl.start;
    bool feed = true;
    d.setParams (held);
    Planes io = tl.input;
    Rng r (0xB10Cu);
    std::size_t next = 0;
    int pos = 0, sliceEdge = 0;
    while (pos < tl.n)
    {
        while (next < tl.events.size() && tl.events[next].at == pos)
            apply (d, held, feed, tl.events[next++], c);
        if (heldEveryBlock)
        {
            d.setParams (held);
            d.setHeads (held.heads);
        }
        if (pos >= sliceEdge)
            sliceEdge = pos + (slice > 0 ? slice : r.range (1, 2000));
        int end = std::min (sliceEdge, tl.n);
        if (next < tl.events.size())
            end = std::min (end, tl.events[next].at);
        for (int e : { tl.monoFrom, tl.monoTo, tl.pauseFrom, tl.pauseTo })
            if (e > pos)
                end = std::min (end, e);
        const bool mono  = pos >= tl.monoFrom && pos < tl.monoTo;
        const bool pause = pos >= tl.pauseFrom && pos < tl.pauseTo;
        const int  nch   = pause ? 0 : (mono ? 1 : 2);
        float* p[2] = { io[0].data() + pos, io[1].data() + pos };
        if (noOps)
        {
            // Calls that must move NOTHING, of every width: zero-length (accepted) and refused ones. A width
            // aimed, a ramp stepped or a head snapped by any of them would show in the bytes.
            const bool zeroes = d.process (p, 1, 0, ! feed) && d.process (nullptr, 0, 0, ! feed) && d.process (p, 2, 0, feed);
            // The refused width reaches process() as a run-time number, as a host's would: a literal 3 lets gcc clone
            // process() for three channels and warn that the clone's channel loops overrun (-Waggressive-loop-
            // optimizations) — a clone the refusal at its top never lets run.
            volatile int tooWide = 3;
            const bool refused = ! d.process (p, tooWide, 64, feed) && ! d.process (p, 1, -1, feed) && ! d.process (nullptr, 1, 64, feed);
            c.noOpsAsSaid = c.noOpsAsSaid && zeroes && refused;
            c.noOpCalls += 6;
        }
        test::run (d.process (nch == 0 ? nullptr : p, nch, end - pos, feed));
        c.monoCalls  += mono ? 1 : 0;
        c.pauseCalls += pause ? 1 : 0;
        pos = end;
    }
    w = witness (d);
    return io;
}

void reslicing()
{
    group ("re-slicing at absolute event samples: heads, TIME, OFFSET, feed, a mono and a pause stretch — any slicing, the same bytes");
    const Timeline tl = makeTimeline();
    Census first;
    Witness w0;
    const Planes ref = renderTimeline (tl, 4096, false, first, w0);

    bool sameAll = true;
    for (int slice : { 1, 64, 333, 517, 0 })
    {
        Census c;
        Witness w;
        const Planes out = renderTimeline (tl, slice, false, c, w);
        const bool here = samePlanes (out, ref) && sameWitness (w, w0);
        if (! here)
            std::printf ("      [slicing %d differs from 4096]\n", slice);
        sameAll = sameAll && here;
    }
    // Catches: anything that steps per block instead of per sample — a ramp, a glide, the mono edge, the pause —
    // from a layout off the plain path (the fast path's own per-call decision is plainPathReslicing's).
    ok (sameAll, "slicings 1, 64, 333, 517 and a seeded random one: the 4096 render's bytes and witnesses");

    bool heldSame = true;
    for (int slice : { 64, 0 })
    {
        Census c;
        Witness w;
        const Planes out = renderTimeline (tl, slice, true, c, w);
        heldSame = heldSame && samePlanes (out, ref) && sameWitness (w, w0);
    }
    // Catches: a re-send that restarts a ramp or a glide (a guard that is merely loose shows in exactGuard).
    ok (heldSame, "the held knobs and layout re-sent at the top of EVERY block (setParams and setHeads): the same bytes");

    Census nc;
    Witness nw;
    const Planes noOp = renderTimeline (tl, 333, false, nc, nw, true);
    // Catches: the width aimed (or a head stepped) by a zero-length or a refused call — e.g. aimWidth() ahead of the
    // n == 0 return: a zero-length MONO call would send side 0 toward the mono level.
    ok (nc.noOpsAsSaid && nc.noOpCalls > 1000 && samePlanes (noOp, ref) && sameWitness (nw, w0),
        "zero-length calls of width 0, 1 and 2 and refused calls before EVERY block (" + std::to_string (nc.noOpCalls)
            + " of them): the same bytes and witnesses");

    for (const Event& e : tl.events)
    {
        if (e.at >= tl.monoFrom && e.at < tl.monoTo) ++first.inMono;
        if (e.at >= tl.pauseFrom && e.at < tl.pauseTo) ++first.inPause;
    }
    std::printf ("      [timeline: %d position, %d level, %d pan, %d time, %d offset, %d feed, %d restart events; %d heads "
                 "turned on, %d off; %d in the mono stretch, %d in the pause; %lld mono and %lld pause calls]\n",
                 first.kinds[0], first.kinds[1], first.kinds[2], first.kinds[3], first.kinds[4], first.kinds[5],
                 first.kinds[6], first.turnedOn, first.turnedOff, first.inMono, first.inPause, first.monoCalls,
                 first.pauseCalls);
    bool everyKind = true;
    for (int k = 0; k < kKinds; ++k)
        everyKind = everyKind && first.kinds[k] >= 10;
    ok (everyKind && first.turnedOn >= 2 && first.turnedOff >= 2 && first.inMono >= 2 && first.inPause >= 2
            && first.monoCalls >= 2 && first.pauseCalls >= 2,
        "the census: every kind of event at least 10 times, heads turned on and off, events inside both stretches");
    ok (first.kinds[(int) Kind::reset] == 3, "the census: the three restarts were applied (in the open, in mono, in the pause)");
}

//==============================================================================
// 12 · A PAUSE MOVES THE HEADS AS SILENCE
//==============================================================================
void pauseIsSilence()
{
    group ("a pause moves the heads as silence: positions, gains, motor and taps bit-identical to n samples of zeros");
    constexpr double sr = 48000.0;
    const Heads lay = layout ({ { 1.0f, 1.0f, 0.0f }, { 0.5f, 0.6f, 0.4f }, { 0.3f, 0.5f, -0.2f } });
    delay::TapeDelay a, b;
    for (delay::TapeDelay* d : { &a, &b })
    {
        test::run (d->prepare (sr, 1024, 2));
        d->setParams (knobs (300.0f, 0.5f, 3000.0f, 6.0f, 1.0f, lay));
        Planes io = burst (3000, 3000, 51u, 0.3f, false);
        test::run (renderBlocks (*d, io, 256));
        d->setHead (0, { 0.4f, 0.8f, 0.3f });
        d->setHead (1, { 0.9f, 0.2f, -0.6f });
        d->setTimeMs (200.0f);
    }
    bool sameAll = true, mid = false;
    for (int round = 0; round < 2; ++round)
    {
        test::run (a.process (nullptr, 0, 700));
        Planes zeros (2, Plane (700, 0.0f));
        test::run (renderBlocks (b, zeros, 700));
        sameAll = sameAll && sameWitness (witness (a), witness (b));
        if (round == 0)
            mid = ! same (a.headPosition (0), (double) 0.4f) && ! same (a.headGain (0, 0), panGains (0.8f, 0.3f)[0]);
    }
    ok (mid, "the first pause ends mid-glide and mid-ramp (the comparison is of moving heads)");
    // Catches: a pause that does not step the heads, or steps them other than per sample.
    ok (sameAll, "after 700 and 1400 samples, paused vs fed zeros: every head's position, gains and tap, the motor, bit for bit");
}

//==============================================================================
// 13 · NO SOUNDING HEAD, THEN THE DEFAULT HEAD BACK
//==============================================================================
void zeroSoundingHeads()
{
    group ("every head silent: the output IS the input while the loop rings; the default head back: the frozen class's bytes");
    constexpr double sr = 48000.0;
    const int len = rampLength (sr), offN = offsetLineLength (sr);
    delay::TapeDelay now;
    frz::TapeDelay then;
    test::run (now.prepare (sr, 64, 2));
    test::run (then.prepare (sr, 64, 2));
    now.setParams (knobs (250.0f, 0.7f, 3000.0f, 5.0f, 1.0f));
    then.setParams (frozenKnobs (250.0f, 0.7f, 3000.0f, 5.0f, 1.0f));

    constexpr int block = 64, tOff = 6016, tOn = 30016, N = 60032;
    const int quietFrom = tOff + len + offN + 128, matchFrom = tOn + len + offN + 128;
    Rng r (61u);
    bool dryOnly = true, match = true;
    double ringing = 0.0;
    for (int t = 0; t < N; t += block)
    {
        if (t == tOff)
            now.setHeads (Heads {});
        if (t == tOn)
            now.setHead (0, { 1.0f, 1.0f, 0.0f });
        const bool feed = t < tOn;
        Plane l ((std::size_t) block), rr ((std::size_t) block);
        for (int i = 0; i < block; ++i)
        {
            cell (l, i) = 0.3f * r.signedUnit();
            cell (rr, i) = 0.3f * r.signedUnit();
        }
        Plane l2 = l, r2 = rr, l0 = l, r0 = rr;
        float* pn[2] = { l.data(), rr.data() };
        float* pt[2] = { l2.data(), r2.data() };
        test::run (now.process (pn, 2, block, feed));
        test::run (then.process (pt, 2, block, feed));
        if (t >= quietFrom && t < tOn)
            dryOnly = dryOnly && sameSpan (l, 0, l0, 0, block) && sameSpan (rr, 0, r0, 0, block);
        if (t >= matchFrom)
        {
            match = match && sameSpan (l, 0, l2, 0, block) && sameSpan (rr, 0, r2, 0, block);
            for (int i = 0; i < block; ++i)
                ringing += std::fabs ((double) get (l2, i) - (double) get (l0, i));
        }
    }
    // Catches: a silent head that still adds (a +0 that is not one, a skipped head that is not skipped).
    ok (dryOnly, "all levels 0, MIX 1: output == input, byte for byte");
    ok (ringing > 1.0, "...while the loop rang on: the stand-by tail the default head finds there is not silence");
    // Catches: a loop that reads a head; an offset line that does not drain; a landed head that is not the 1x tap.
    ok (match, "head 0 back on at 1x, feed off: after its ramp and the offset line's drain, the frozen class's bytes");
}

//==============================================================================
// 14 · AWAY AND BACK
//==============================================================================
void awayAndBack()
{
    group ("away and back: every head wandered, the default re-sent; once landed and drained, the frozen class's bytes");
    constexpr double sr = 48000.0;
    const int offN = offsetLineLength (sr);
    delay::TapeDelay now;
    frz::TapeDelay then;
    test::run (now.prepare (sr, 64, 2));
    test::run (then.prepare (sr, 64, 2));
    now.setParams (knobs (280.0f, 0.5f, 3000.0f, -7.0f, 0.6f));
    then.setParams (frozenKnobs (280.0f, 0.5f, 3000.0f, -7.0f, 0.6f));

    constexpr int block = 64, tAway = 3008, tTime = 8000, tBack = 15040, maxN = 420000;
    Rng r (71u);
    int landedAt = -1, compareFrom = -1;
    long long compared = 0, differedAway = 0;
    bool match = true;
    for (int t = 0; t < maxN; t += block)
    {
        if (t == tAway)
            now.setHeads (layout ({ { 0.6f, 0.8f, 0.5f }, { 0.3f, 0.9f, -1.0f }, { 0.8f, 0.5f, 0.2f }, { 0.45f, 0.7f, -0.3f } }));
        if (t == tTime)
        {
            now.setTimeMs (200.0f);
            then.setTimeMs (200.0f);
        }
        if (t == tBack)
        {
            now.setParams (knobs (200.0f, 0.5f, 3000.0f, -7.0f, 0.6f));
            then.setParams (frozenKnobs (200.0f, 0.5f, 3000.0f, -7.0f, 0.6f));
        }
        Plane l ((std::size_t) block), rr ((std::size_t) block);
        for (int i = 0; i < block; ++i)
        {
            const int at = t + i;
            const bool loud = at < 30000 || at % 10000 < 500;
            const float tone = loud ? 0.3f * (float) std::sin (0.05 * (double) at) : 0.0f;
            cell (l, i)  = loud ? tone + 0.2f * r.signedUnit() : 0.0f;
            cell (rr, i) = loud ? -tone + 0.2f * r.signedUnit() : 0.0f;
        }
        Plane l2 = l, r2 = rr;
        float* pn[2] = { l.data(), rr.data() };
        float* pt[2] = { l2.data(), r2.data() };
        test::run (now.process (pn, 2, block));
        test::run (then.process (pt, 2, block));
        const bool sameBlock = sameSpan (l, 0, l2, 0, block) && sameSpan (rr, 0, r2, 0, block);
        if (t >= tAway && t < tBack && ! sameBlock)
            ++differedAway;

        if (landedAt < 0 && t > tBack)
        {
            bool home = same (now.headPosition (0), 1.0) && same (now.headGain (0, 0), 1.0f) && same (now.headGain (0, 1), 1.0f);
            for (int k = 1; k < kHeads; ++k)
                home = home && same (now.headGain (k, 0), 0.0f) && same (now.headGain (k, 1), 0.0f);
            if (home)
            {
                landedAt = t + block;
                compareFrom = landedAt + offN + 2 * block;
            }
        }
        if (compareFrom >= 0 && t >= compareFrom)
        {
            match = match && sameBlock;
            compared += block;
            if (compared >= 40000)
                break;
        }
    }
    ok (differedAway > 10, "while the heads wandered the sound was not the one-head stage's");
    ok (landedAt > 0, "head 0 glides home to 1x exactly at unity, the others silent (by sample " + std::to_string (landedAt) + ")");
    // Catches: a loop that read the wandering heads (the line would remember it); a head that parks short of 1.
    ok (compared >= 40000 && match, "from the landing + the offset line's drain on: the frozen class's bytes for 40000 samples");
}

//==============================================================================
// 15 · LAW 8
//==============================================================================
void law8()
{
    group ("law 8: levels fading to 0 under a ringing loop with OFFSET, tiny levels ramping, subnormal input — no subnormal anywhere");
    constexpr double sr = 8000.0;   // the witness scans every line: the slowest useful rate keeps it cheap
    constexpr int N = 120000, block = 100;
    const int tiny = std::min (1, kHeads - 1), fading = std::min (2, kHeads - 1);
    delay::TapeDelay d;
    test::run (d.prepare (sr, block, 2));
    d.setParams (knobs (100.0f, 0.6f, 2000.0f, 13.0f, 1.0f,
                        layout ({ { 1.0f, 1.0f, 0.0f }, { 0.37f, 0.0f, 0.0f }, { 0.71f, 0.5f, -0.3f } })));
    Rng r (81u);
    Planes io (2, Plane ((std::size_t) N, 0.0f));
    for (int i = 0; i < N; ++i)
        for (int ch = 0; ch < 2; ++ch)
        {
            float v = 0.0f;
            if (i < 1500)
                v = 0.5f * r.signedUnit();
            else if (i < 5000)
                v = (r.chance (50) ? 1.0f : -1.0f) * std::numeric_limits<float>::denorm_min() * (float) r.range (1, 1000);
            cell (io[(std::size_t) ch], i) = v;
        }

    int subnormalBlocks = 0, toggles = 0;
    bool tinyOn = false;
    double rang = 0.0;
    for (int t = 0; t < N; t += block)
    {
        if (t == 6000)
            d.setHead (fading, { 0.71f, 0.0f, -0.3f });   // fades under the ringing loop
        if (t == 12000)
            d.setHead (0, { 1.0f, 0.0f, 0.0f });
        if (t >= 12000 && t < 100000 && t % 1100 == 0)
        {
            tinyOn = ! tinyOn;
            const float pan = (toggles / 2) % 2 == 0 ? std::nextafter (1.0f, 0.0f) : std::nextafter (-1.0f, 0.0f);
            d.setHead (tiny, { 0.37f, tinyOn ? std::nextafter (1.0e-15f, 1.0f) : 0.0f, pan });
            ++toggles;
        }
        if (t == 100000)
        {
            d.setHead (tiny, { 0.37f, 0.0f, 0.0f });
            d.setHead (0, { 1.0f, 1.0f, 0.0f });
        }
        float* p[2] = { io[0].data() + t, io[1].data() + t };
        test::run (d.process (p, 2, block));
        bool sub = d.anySubnormalState();
        for (int k = 0; k < kHeads; ++k)
            sub = sub || isSub (d.headGain (k, 0)) || isSub (d.headGain (k, 1));
        subnormalBlocks += sub ? 1 : 0;
        if (t >= 5000 && t < 12000)
            for (int i = 0; i < block; ++i)
                rang += std::fabs ((double) get (io[0], t + i));
    }
    // Catches: contributions that are not flushed (a tiny gain's product is subnormal), a level that keeps 1e-38.
    ok (subnormalBlocks == 0, "no subnormal state or gain at any of " + std::to_string (N / block) + " block ends ("
                                  + std::to_string (subnormalBlocks) + " had one)");
    ok (toggles > 50 && rang > 1.0e-3, "the loop rang after the input and the tiny head ramped in and out " + std::to_string (toggles) + " times");
    bool silent = true;
    for (int ch = 0; ch < 2; ++ch)
        for (int i = 110000; i < N; ++i)
            silent = silent && bitsOf (get (io[(std::size_t) ch], i)) == 0u;
    ok (silent, "the tail reaches exact digital silence: the last 10000 samples are +0 on both lines, the 1x head at unity");
}

//==============================================================================
// 16 · THE TAIL
//==============================================================================
double tailFor (float timeMs, float repeats, float offsetMs, double a)
{
    const double t = (double) timeMs * 0.001, g = (double) repeats, off = std::fabs ((double) offsetMs) * 0.001;
    if (g <= 0.0) return off + t;
    if (g >= 1.0) return off + 1.0e6;
    return a > 0.001 ? off + t * (1.0 + std::log (0.001 / a) / std::log (g)) : off + t;
}

void tail()
{
    group ("the tail: the heads' summed first-echo level A in the count of passes; the default layout's bits unchanged");
    const auto rel = [] (double got, double want) { return std::fabs (got - want) <= 1.0e-6 * std::fabs (want); };
    {
        delay::TapeDelay d;
        Heads all {};
        for (auto& h : all)
            h = { 1.0f, 1.0f, 0.0f };
        d.setParams (knobs (400.0f, 0.5f, 4000.0f, 12.0f, 0.3f, all));
        // Catches a tail that ignores the heads.
        ok (rel ((double) d.tailSeconds(), tailFor (400.0f, 0.5f, 12.0f, (double) kHeads)),
            "every head at unity, centre: A = kMaxHeads (" + std::to_string (d.tailSeconds()) + " s)");
        d.setHeads (layout ({ { 0.5f, 1.0f, -1.0f } }));
        ok (rel ((double) d.tailSeconds(), tailFor (400.0f, 0.5f, 12.0f, (double) kSqrt2f)), "one head hard left: A = float sqrt 2");
        d.setHeads (layout ({ { 0.5f, 0.4f, 0.5f } }));
        ok (rel ((double) d.tailSeconds(), tailFor (400.0f, 0.5f, 12.0f, (double) panGains (0.4f, 0.5f)[1])),
            "level 0.4 at pan 0.5: A = the louder side's gain");
        const float before = d.tailSeconds();
        d.setHead (kHeads - 1, { 0.3f, 0.0f, 0.9f });
        if (kHeads > 1)
            ok (same (d.tailSeconds(), before), "a silent head's pan and position move no estimate");
        d.setHeads (Heads {});
        ok (same (d.tailSeconds(), (float) (std::fabs ((double) 12.0f) * 0.001 + (double) 400.0f * 0.001)),
            "every head off: offset + one TIME, exactly");
        d.setHeads (layout ({ { 1.0f, 0.0009f, 0.0f } }));
        ok (same (d.tailSeconds(), (float) (std::fabs ((double) 12.0f) * 0.001 + (double) 400.0f * 0.001)),
            "A under a thousandth: offset + one TIME, exactly");
        d.setHeads (all);
        d.setRepeats (1.0f);
        ok (same (d.tailSeconds(), (float) (std::fabs ((double) 12.0f) * 0.001 + 1.0e6)), "REPEATS 1 with every head on: offset + 1e6");
        d.setRepeats (0.0f);
        ok (same (d.tailSeconds(), (float) (std::fabs ((double) 12.0f) * 0.001 + (double) 400.0f * 0.001)),
            "REPEATS 0 with every head on: offset + one TIME");
    }
    {
        delay::TapeDelay now;
        frz::TapeDelay then;
        bool sameTail = same (now.tailSeconds(), then.tailSeconds());
        const float sets[][3] = { { 350.0f, 0.35f, 0.0f }, { 23.0f, 0.97f, -11.0f }, { 6100.0f, 0.999f, 31.0f },
                                  { 1.0f, 1.0e-6f, 0.5f }, { 500.0f, 0.0f, -3.0f }, { 777.0f, 1.0f, 30.0f } };
        for (const auto& s : sets)
        {
            now.setParams (knobs (s[0], s[1], 4000.0f, s[2], 0.25f));
            then.setParams (frozenKnobs (s[0], s[1], 4000.0f, s[2], 0.25f));
            sameTail = sameTail && same (now.tailSeconds(), then.tailSeconds());
            now.setHeads (layout ({ { 0.3f, 1.0f, 0.5f }, { 0.7f, 0.6f, -1.0f } }));
            now.setHeads (defaultHeads());
            sameTail = sameTail && same (now.tailSeconds(), then.tailSeconds());
        }
        ok (sameTail, "the default layout's tail is the frozen class's bits: at construction and at six knob sets, after a detour");
    }
}

//==============================================================================
// 17 · THE SETTER'S BOUNDS
//==============================================================================
void validation()
{
    group ("setHead's bounds: an index off the array, a field that is not a number, the clamps, a level of 1e-38, -0 pan");
    constexpr double sr = 8000.0;
    const int len = rampLength (sr);
    const Heads lay = layout ({ { 0.6f, 1.0f, 0.3f }, { 0.4f, 0.5f, -0.5f } });

    // Twins: the same input, one of them sent `a` and the other `b`; the bytes and the witnesses must agree throughout.
    const auto twins = [&] (auto sendA, auto sendB)
    {
        delay::TapeDelay a, b;
        bool s = true;
        for (delay::TapeDelay* d : { &a, &b })
        {
            test::run (d->prepare (sr, 256, 2));
            d->setParams (knobs (200.0f, 0.5f, 3000.0f, 4.0f, 1.0f, lay));
        }
        Planes io = burst (4000, 4000, 92u, 0.3f, false);
        Planes ia = io, ib = io;
        float* pa[2] = { ia[0].data(), ia[1].data() };
        float* pb[2] = { ib[0].data(), ib[1].data() };
        test::run (a.process (pa, 2, 500));
        test::run (b.process (pb, 2, 500));
        sendA (a);
        sendB (b);
        s = s && sameWitness (witness (a), witness (b));
        for (int i = 500; i < 4000; i += 500)
        {
            float* qa[2] = { ia[0].data() + i, ia[1].data() + i };
            float* qb[2] = { ib[0].data() + i, ib[1].data() + i };
            test::run (a.process (qa, 2, 500));
            test::run (b.process (qb, 2, 500));
            s = s && sameWitness (witness (a), witness (b));
        }
        return s && samePlanes (ia, ib);
    };
    const auto nothing = [] (delay::TapeDelay&) {};

    // Catches: an index guard that lets k = kMaxHeads or a negative one through.
    ok (twins ([] (delay::TapeDelay& d) { d.setHead (-1, { 0.5f, 1.0f, -1.0f }); d.setHead (kHeads, { 0.2f, 1.0f, 1.0f });
                                           d.setHead (1000, { 0.1f, 0.5f, 0.0f }); }, nothing),
        "k = -1, kMaxHeads and 1000 change nothing: the same bytes and witnesses as no call");
    ok (twins ([] (delay::TapeDelay& d) { d.setHead (0, { kNaN, 0.5f, kNaN }); },
               [] (delay::TapeDelay& d) { d.setHead (0, { 0.6f, 0.5f, 0.3f }); }),
        "{NaN, 0.5, NaN} changes only the level");
    const float inf = std::numeric_limits<float>::infinity();
    ok (twins ([=] (delay::TapeDelay& d) { d.setHead (0, { inf, -inf, kNaN }); }, nothing),
        "{inf, -inf, NaN} changes nothing");
    const float fmin = std::numeric_limits<float>::min();
    ok (twins ([] (delay::TapeDelay& d) { d.setHead (0, { 0.0f, 0.5f, 2.0f }); },
               [=] (delay::TapeDelay& d) { d.setHead (0, { fmin, 0.5f, 1.0f }); }),
        "position 0 is FLT_MIN, pan 2 is 1");
    ok (twins ([] (delay::TapeDelay& d) { d.setHead (0, { -3.0f, 2.0f, -2.0f }); },
               [=] (delay::TapeDelay& d) { d.setHead (0, { fmin, 1.0f, -1.0f }); }),
        "position -3 is FLT_MIN, level 2 is 1, pan -2 is -1");
    ok (twins ([] (delay::TapeDelay& d) { d.setHead (0, { 2.0f, 0.3f, -0.0f }); },
               [] (delay::TapeDelay& d) { d.setHead (0, { 1.0f, 0.3f, 0.0f }); }),
        "position 2 is 1, pan -0 is the centre");

    {
        // The clamp, witnessed: a silent head sent below the floor stands on FLT_MIN after one sample.
        delay::TapeDelay d;
        test::run (d.prepare (sr, 64, 2));
        d.setParams (knobs (200.0f, 0.5f, 3000.0f, 0.0f, 1.0f, Heads {}));
        OneByOne s { d };
        test::run (s.step (0.0f, 0.0f));
        d.setHead (0, { 0.0f, 0.0f, 0.0f });
        test::run (s.step (0.0f, 0.0f));
        const bool atFloor = same (d.headPosition (0), (double) fmin);
        d.setHead (0, { -3.0f, 0.0f, 0.0f });
        test::run (s.step (0.0f, 0.0f));
        const bool stillFloor = same (d.headPosition (0), (double) fmin);
        d.setHead (0, { 2.0f, 0.0f, 0.0f });
        test::run (s.step (0.0f, 0.0f));
        ok (atFloor && stillFloor && same (d.headPosition (0), 1.0), "a silent head sent to 0 or -3 stands on FLT_MIN, sent to 2 on 1");
    }
    {
        // Catches a level of 1e-38 kept: a head that never turns off, its gain subnormal.
        delay::TapeDelay d;
        test::run (d.prepare (sr, 64, 2));
        d.setParams (knobs (200.0f, 0.5f, 3000.0f, 4.0f, 1.0f, layout ({ { 1.0f, 0.5f, 0.0f } })));
        OneByOne s { d };
        Rng r (93u);
        for (int i = 0; i < 400; ++i)
            test::run (s.step (0.2f * r.signedUnit(), 0.2f * r.signedUnit()));
        d.setHead (0, { 1.0f, 1.0e-38f, 0.0f });
        bool normal = true;
        for (int m = 0; m < len + 2; ++m)
        {
            test::run (s.step (0.2f * r.signedUnit(), 0.2f * r.signedUnit()));
            normal = normal && ! isSub (d.headGain (0, 0)) && ! isSub (d.headGain (0, 1)) && ! d.anySubnormalState();
        }
        const bool off = bitsOf (d.headGain (0, 0)) == 0u && bitsOf (d.headGain (0, 1)) == 0u;
        const Witness w = witness (d);
        d.setHead (0, { 1.0f, 0.0f, 0.0f });
        ok (normal && off && sameWitness (witness (d), w),
            "level 1e-38 is 0: the gains land on +0 with nothing subnormal, and a 0 sent after it moves nothing");
    }
    {
        delay::TapeDelay d;
        test::run (d.prepare (sr, 64, 2));
        d.setParams (knobs (200.0f, 0.5f, 3000.0f, 4.0f, 1.0f, layout ({ { 1.0f, 0.7f, -0.0f } })));
        OneByOne s { d };
        test::run (s.step (0.0f, 0.0f));
        ok (same (d.headGain (0, 0), 0.7f) && same (d.headGain (0, 1), 0.7f), "pan -0: both sides at the level exactly, the centre");
    }
}

//==============================================================================
// 18 · UNDER ONE SAMPLE
//==============================================================================
void underOneSample()
{
    group ("under one sample: a head at 1e-6 of 1 ms plays an impulse at exactly one sample — the motor's floor");
    constexpr double sr = 48000.0;
    constexpr int N = 200, D = 48;
    Planes in (2, Plane ((std::size_t) N, 0.0f));
    cell (in[0], 0) = 1.0f;
    cell (in[1], 0) = -0.5f;

    frz::TapeDelay then;
    test::run (then.prepare (sr, 64, 2));
    then.setParams (frozenKnobs (1.0f, 0.0f, 4000.0f, 0.0f, 1.0f));
    Planes wetD = in;
    test::run (renderBlocks (then, wetD, 64));

    const auto play = [&] (float position)
    {
        delay::TapeDelay d;
        test::run (d.prepare (sr, 64, 2));
        d.setParams (knobs (1.0f, 0.0f, 4000.0f, 0.0f, 1.0f, layout ({ { position, 1.0f, 0.0f } })));
        Planes out = in;
        test::run (renderBlocks (d, out, 64));
        return out;
    };
    const Planes a = play (1.0e-6f);
    // Catches: a head distance off by one; a floor that is not the motor's one sample.
    ok (same (get (a[0], 0), 1.0f) && same (get (a[1], 0), -0.5f), "index 0: the dry alone, the head has nothing yet");
    ok (get (a[0], 1) > 0.1f && sameSpan (a[0], 1, wetD[0], D, 150) && sameSpan (a[1], 1, wetD[1], D, 150),
        "from index 1 on: the line one sample back, byte for byte (the frozen 48-sample wet moved 47 earlier)");
    ok (same (readDistance ((double) (1.0f / 48.0f), 48.0, sr), 1.0f), "position 1/48 of 48 samples is a distance of exactly 1");
    ok (samePlanes (play (1.0f / 48.0f), a) && samePlanes (play (0.0f), a),
        "a distance of exactly one sample, and position 0 (the floor), play the same bytes");
}

//==============================================================================
// 19 · THE PICTURE'S TAPS
//==============================================================================
void shownHeads()
{
    group ("shownHeadMs: the target x TIME before the stage runs, the glided read after each call, 0 off the array");
    // Catches: a picture that shows the target while the head glides, one that misses TIME before the first call,
    // a tap without the motor's floor, an index off the array that reads anything but 0.
    delay::TapeDelay d;
    const float third = 1.0f / 3.0f;
    ok (same (d.shownHeadMs (0), 350.0f), "constructed: the default head at 1 x 350 ms");
    d.setHead (0, { third, 1.0f, 0.0f });
    ok (same (d.shownHeadMs (0), (float) ((double) third * 350.0)), "before prepare: position x TIME");
    d.setTimeMs (500.0f);
    ok (same (d.shownHeadMs (0), (float) ((double) third * 500.0)), "...and it follows setTimeMs");
    if (kHeads > 1)
    {
        d.setHead (kHeads - 1, { 0.25f, 0.0f, 0.0f });
        ok (same (d.shownHeadMs (kHeads - 1), 125.0f), "a silent head is shown where it stands too");
    }
    constexpr double sr = 48000.0;
    test::run (d.prepare (sr, 256, 2));
    d.setTimeMs (400.0f);
    ok (same (d.shownHeadMs (0), (float) ((double) third * 400.0)), "prepared, before the first call: still position x TIME");

    const auto glided = [&] (int k) { return (float) ((double) readDistance (d.headPosition (k), d.headSamples(), sr) / sr * 1000.0); };
    Planes io (2, Plane (256, 0.0f));
    test::run (renderBlocks (d, io, 256));
    bool all = true;
    for (int k = 0; k < kHeads; ++k)
        all = all && same (d.shownHeadMs (k), glided (k));
    ok (all, "after a call: every head's (float)(floored distance / sr x 1000), bit for bit");

    // Catches a running stage that publishes a head's target from the setter (the picture would jump to where the
    // head is SENT, then back to where it reads after the next call).
    const float shownBefore = d.shownHeadMs (0);
    d.setHead (0, { 0.8f, 1.0f, 0.0f });
    ok (same (d.shownHeadMs (0), shownBefore), "a running stage: setHead moves no tap until the next call");
    test::run (renderBlocks (d, io, 256));
    ok (same (d.shownHeadMs (0), glided (0)) && d.shownHeadMs (0) < 0.8f * 400.0f,
        "while it glides: where it reads, not where it is sent");
    const float shownNow = d.shownHeadMs (0);
    d.setTimeMs (100.0f);
    ok (same (d.shownHeadMs (0), shownNow), "a running stage publishes after its calls: a knob between them moves no tap");

    d.setHead (0, { 0.0f, 1.0f, 0.0f });
    Planes more (2, Plane (60000, 0.0f));
    test::run (renderBlocks (d, more, 4096));
    test::run (renderBlocks (d, more, 4096));
    test::run (renderBlocks (d, more, 4096));
    test::run (renderBlocks (d, more, 4096));
    test::run (renderBlocks (d, more, 4096));
    ok (same (d.shownHeadMs (0), (float) (1.0 / sr * 1000.0)), "a head under one sample is shown at the floor, one sample");
    ok (same (d.shownHeadMs (-1), 0.0f) && same (d.shownHeadMs (kHeads), 0.0f) && same (d.shownHeadMs (1000), 0.0f)
            && same (d.headPosition (-1), 0.0) && same (d.headGain (kHeads, 0), 0.0f) && same (d.headGain (0, 2), 0.0f),
        "an index off the array reads 0");
}

//==============================================================================
// 15b · LAW 8 FOR A RAMP RESTARTED FROM WHERE IT STOOD — a level sent up and back down before every sample
//==============================================================================
void restartedRampFloor()
{
    group ("law 8: a level sent up and back down before every sample — each ramp restarted from where the last stood — "
           "never parks a subnormal gain");
    for (double sr : { 1000.0, 8000.0 })
    {
        delay::TapeDelay d;
        test::run (d.prepare (sr, 1, 2));
        OneByOne s { d };
        test::run (s.step (0.0f, 0.0f));
        const int n = (int) (2.0 * 87.4 / (1.0 / (double) rampLength (sr)));   // twice the samples to shrink 1 to 2^-126
        bool anySub = false;
        int zeroAt = -1;
        for (int i = 0; i < n; ++i)
        {
            d.setHead (0, { 1.0f, 1.0f, 0.0f });
            d.setHead (0, { 1.0f, 0.0f, 0.0f });
            test::run (s.step (0.0f, 0.0f));
            anySub = anySub || d.anySubnormalState() || isSub (d.headGain (0, 0)) || isSub (d.headGain (0, 1));
            if (zeroAt < 0 && same (d.headGain (0, 0), 0.0f) && same (d.headGain (0, 1), 0.0f))
                zeroAt = i;
        }
        // Catches: a ramp value stored without the floor at the subnormal edge (the gain then parks at ~2^-146).
        ok (! anySub && zeroAt > 0, "at " + std::to_string ((int) sr) + " Hz: no subnormal in any state over "
                                     + std::to_string (n) + " samples, and the gain reaches exactly 0 (at sample "
                                     + std::to_string (zeroAt) + ")");

        // ...and nothing small is left behind in the ramp's ends: send level 1, let it land.
        d.setHead (0, { 1.0f, 1.0f, 0.0f });
        for (int i = 0; i < rampLength (sr) + 2; ++i)
            test::run (s.step (0.0f, 0.0f));
        ok (same (d.headGain (0, 0), 1.0f) && same (d.headGain (0, 1), 1.0f) && ! d.anySubnormalState(),
            "...then level 1 lands on exactly 1 and no state (ramp ends included) is subnormal");
    }
}

//==============================================================================
// 15c · A REAL POISON RESTART — a finite overflow, not a rejected NaN — with the heads moving
//==============================================================================
void poisonRestart()
{
    group ("poison: a finite overflow restarts the line at the call's end, with heads mid-ramp — what follows is a fresh "
           "restart's bytes");
    constexpr double sr = 48000.0;
    const float big = std::numeric_limits<float>::max();
    const delay::TapeDelayParams p = knobs (20.0f, 1.0f, 100000.0f, 9.0f, 1.0f,
                                            layout ({ { 1.0f, 0.7f, -0.4f }, { 0.5f, 1.0f, 1.0f }, { 0.25f, 0.5f, 0.0f } }));
    delay::TapeDelay a, b;
    test::run (a.prepare (sr, 512, 2));
    test::run (b.prepare (sr, 512, 2));
    a.setParams (p);
    b.setParams (p);

    Planes warm = burst (4000, 4000, 7u, 0.5f, false);
    test::run (renderBlocks (a, warm, 512));
    a.setHead (0, { 0.6f, 0.2f, 0.9f });   // a head mid-glide and mid-ramp when the overflow lands

    // Finite samples that overflow the loop filter: +FLT_MAX then -FLT_MAX, DARK wide open, REPEATS 1.
    Planes boom (2, Plane (256, 0.0f));
    for (int i = 0; i < 256; ++i)
        cell (boom[0], i) = cell (boom[1], i) = (i % 2 == 0 ? big : -big);
    float* q[2] = { boom[0].data(), boom[1].data() };
    test::run (a.process (q, 2, 256));

    // The twin restarts at the same call boundary by hand, with the same targets.
    b.setHead (0, { 0.6f, 0.2f, 0.9f });
    b.reset();

    Planes after = burst (30000, 3000, 9u, 0.5f, false);
    Planes ia = after, ib = after;
    test::run (renderBlocks (a, ia, 300));
    test::run (renderBlocks (b, ib, 300));
    bool finite = true;
    for (const auto& pl : ia)
        for (float v : pl)
            finite = finite && std::isfinite (v);
    // Catches: no restart at the end of a call whose loop state went non-finite (the line would replay inf / NaN).
    ok (finite && samePlanes (ia, ib) && sameWitness (witness (a), witness (b)),
        "after the overflow call: every sample finite, and the bytes and witnesses of a twin reset at that boundary");
}

//==============================================================================
// 21 · WAKING AFTER A RESTART — a head sent between the pause and the first carried call stands on its target
//==============================================================================
void wakeSnap()
{
    group ("waking: after prepare or reset, a pause, then a head turned on at 1/2 — it stands on 1/2, the sound as with no pause");
    constexpr double sr = 48000.0;
    const int k = kHeads - 1;                                      // off at 1 in the layout below
    const Heads lay = kHeads > 1 ? defaultHeads() : Heads {};
    constexpr int N = 18000, L = 6000;   // the head at 1/2 of 300 ms plays its first echo at 7200: N runs well past it
    const Planes noise = burst (N, L, 101u, 0.4f, false);

    // restart 0: a fresh prepare; 1: a reset after a run. TIME, the head and OFFSET 0 put every read on a whole sample,
    // so the clocks a pause moves on cannot move a bit.
    const auto run = [&] (int restart, bool pause, double& afterOne)
    {
        delay::TapeDelay d;
        test::run (d.prepare (sr, 512, 2));
        d.setParams (knobs (300.0f, 0.5f, 3000.0f, 0.0f, 1.0f, lay));
        if (restart == 1)
        {
            Planes past = burst (3000, 3000, 102u, 0.4f, false);
            test::run (renderBlocks (d, past, 256));
            d.reset();
        }
        if (pause)
            test::run (d.process (nullptr, 0, 300));
        d.setHead (k, { 0.5f, 1.0f, 0.0f });
        Planes io = noise;
        float* p[2] = { io[0].data(), io[1].data() };
        test::run (d.process (p, 2, 1));
        afterOne = d.headPosition (k);
        for (int i = 1; i < N; i += 512)
        {
            float* q[2] = { io[0].data() + i, io[1].data() + i };
            test::run (d.process (q, 2, std::min (512, N - i)));
        }
        return io;
    };
    for (int restart = 0; restart < 2; ++restart)
    {
        double withPause = 0.0, without = 0.0;
        const Planes a = run (restart, true, withPause), b = run (restart, false, without);
        const std::string what = restart == 0 ? "prepare" : "reset";
        // Catches a wake that aims the gains but leaves the head where the pause left it: it would glide from 1 to 1/2
        // at full level from the first sample.
        ok (same (withPause, (double) 0.5f), what + " -> pause -> head on at 1/2 -> one stereo sample: on 1/2 exactly");
        bool heard = false;   // past the burst the dry is silent: what sounds is the heads' wet
        for (int i = L; i < N; ++i)
            heard = heard || std::fabs (get (a[0], i)) > 1.0e-3f;
        ok (heard && samePlanes (a, b), what + " -> pause -> head on: the same bytes as with no pause, echoes included");
    }
}

//==============================================================================
// 22 · THE FAST PATH IS DECIDED PER CALL — its own re-slicing case, from the default layout
//==============================================================================
void plainPathReslicing()
{
    group ("the default layout's fast path, per call: level, pan, position and a second head sent at absolute samples — slicings agree");
    constexpr double sr = 48000.0;
    constexpr int N = 40000;
    struct Send { int at; int head; Head h; };
    const Send sends[] = {
        { 3000,  0, { 1.0f, 0.7f, 0.0f } },    // a level: a ramp that starts from exactly 1
        { 7000,  0, { 1.0f, 1.0f, 0.0f } },    // back to the default (plain again once landed)
        { 11000, 0, { 1.0f, 1.0f, 0.3f } },    // a pan
        { 15000, 0, { 1.0f, 1.0f, 0.0f } },
        { 19000, 1, { 0.5f, 0.8f, -0.4f } },   // a second head turned on
        { 24000, 1, { 0.5f, 0.0f, -0.4f } },
        { 30000, 0, { 0.6f, 1.0f, 0.0f } },    // a position: the 1x head leaves the loop's tap
    };
    constexpr int nSends = (int) (sizeof (sends) / sizeof (sends[0]));
    Rng r (111u);
    Planes in (2, Plane ((std::size_t) N, 0.0f));
    for (int i = 0; i < N; ++i)
    {
        const float tone = 0.3f * (float) std::sin (0.0613 * (double) i);
        cell (in[0], i) = tone + 0.1f * r.signedUnit();
        cell (in[1], i) = -tone + 0.1f * r.signedUnit();
    }
    // slice 0: one call from each send to the next (every call 4000 samples or more)
    const auto render = [&] (int slice)
    {
        delay::TapeDelay d;
        test::run (d.prepare (sr, 4096, 2));
        d.setParams (knobs (250.0f, 0.5f, 3000.0f, 4.0f, 0.8f));
        Planes io = in;
        int pos = 0, next = 0;
        while (pos < N)
        {
            while (next < nSends && sends[next].at == pos)
            {
                d.setHead (sends[next].head, sends[next].h);
                ++next;
            }
            int end = slice > 0 ? std::min (pos + slice, N) : N;
            if (next < nSends)
                end = std::min (end, sends[next].at);
            float* p[2] = { io[0].data() + pos, io[1].data() + pos };
            test::run (d.process (p, 2, end - pos));
            pos = end;
        }
        return io;
    };
    const Planes whole = render (0);
    // Catches a fast path chosen from the gains alone: a call that starts at a gain of exactly 1 with a ramp, a glide
    // or a second head on its way would play the loop's tap for the whole call.
    ok (samePlanes (render (1), whole) && samePlanes (render (64), whole),
        "slicings 1 and 64 against one call per gap (4000+ samples): the same bytes");
}

//==============================================================================
// 23 · THE ORDER OF THE SUM — three heads sounding at once
//==============================================================================
void orderedSum()
{
    group ("three heads sounding at once (1/3, 2/3, 1 at 0.3, 0.6, 0.9): scaled, flushed, summed in index order — byte for byte");
    if (! haveHeads (3, "the ordered sum"))
        return;

    constexpr double sr = 48000.0;
    constexpr int D = 14400, L = 12000, N = 3 * D;   // the burst outlasts twice the heads' spacing: all three overlap
    const float pos[3] = { 1.0f / 3.0f, 2.0f / 3.0f, 1.0f };
    const float lvl[3] = { 0.3f, 0.6f, 0.9f };
    const Planes in = burst (N, L, 121u, 0.6f, false);

    frz::TapeDelay then;
    test::run (then.prepare (sr, 512, 2));
    then.setParams (frozenKnobs (300.0f, 0.5f, 2500.0f, 0.0f, 1.0f));
    Planes wetD = in;
    test::run (renderBlocks (then, wetD, 512));

    delay::TapeDelay now;
    test::run (now.prepare (sr, 512, 2));
    now.setParams (knobs (300.0f, 0.5f, 2500.0f, 0.0f, 1.0f,
                          layout ({ { pos[0], lvl[0], 0.0f }, { pos[1], lvl[1], 0.0f }, { pos[2], lvl[2], 0.0f } })));
    Planes out = in;
    test::run (renderBlocks (now, out, 512));

    int dist[3];
    for (int k = 0; k < 3; ++k)
        dist[k] = (int) readDistance ((double) pos[k], (double) D, sr);
    int bad = 0, allThree = 0, orderShows = 0;
    for (int ch = 0; ch < 2; ++ch)
        for (int t = L; t < N - (D - dist[0]); ++t)
        {
            float c[3];
            bool every = true;
            for (int k = 0; k < 3; ++k)
            {
                c[k] = scaled (get (wetD[(std::size_t) ch], t + D - dist[k]), lvl[k]);
                every = every && ! same (c[k], 0.0f);
            }
            const float forward = (c[0] + c[1]) + c[2], backward = (c[2] + c[1]) + c[0];
            bad += same (get (out[(std::size_t) ch], t), 0.0f + forward) ? 0 : 1;
            allThree += every ? 1 : 0;
            orderShows += same (forward, backward) ? 0 : 1;
        }
    ok (allThree > 4000 && orderShows > 100, "the fixture: " + std::to_string (allThree) + " samples with three copies at once, "
                                                 + std::to_string (orderShows) + " where the order of the sum changes the bits");
    // Catches a sum taken in any other order than the heads' index order.
    ok (bad == 0, "each sample is ((c0 + c1) + c2) of the scaled, flushed copies, byte for byte (" + std::to_string (bad) + " differ)");
}

//==============================================================================
// 24 · THE EXACT GUARD — a level one float away is a new level
//==============================================================================
void exactGuard()
{
    group ("the setter's guard is exact: a level one float above the one that landed lands on that float");
    constexpr double sr = 48000.0;
    const int len = rampLength (sr);
    delay::TapeDelay d;
    test::run (d.prepare (sr, 64, 2));
    d.setParams (knobs (300.0f, 0.3f, 4000.0f, 0.0f, 1.0f, layout ({ { 1.0f, 0.5f, 0.0f } })));
    OneByOne s { d };
    for (int i = 0; i < 4; ++i)
        test::run (s.step (0.0f, 0.0f));
    const bool landed = same (d.headGain (0, 0), 0.5f) && same (d.headGain (0, 1), 0.5f);
    bool exact = true;
    for (float next : { std::nextafter (0.5f, 1.0f), 0.5f, std::nextafter (0.5f, 0.0f) })
    {
        d.setHead (0, { 1.0f, next, 0.0f });
        for (int m = 0; m < len; ++m)
            test::run (s.step (0.0f, 0.0f));
        exact = exact && same (d.headGain (0, 0), next) && same (d.headGain (0, 1), next);
    }
    // Catches a guard within a float epsilon: the one-ulp request would be taken for the level it already holds.
    ok (landed && exact, "0.5 lands; then 0.5 + 1 ulp, 0.5 and 0.5 - 1 ulp each land on exactly that float, both sides");
}

//==============================================================================
// 25 · A RETARGET MID-RAMP STARTS FROM WHERE THE GAIN STANDS
//==============================================================================
void retargetFromWhereItStands()
{
    group ("a retarget mid-ramp starts from where the gain stands: 0 -> 1, then 0.2 at sample 300 — no jump, landing N later");
    constexpr double sr = 48000.0;
    const int len = rampLength (sr);
    delay::TapeDelay d;
    test::run (d.prepare (sr, 64, 2));
    d.setParams (knobs (300.0f, 0.3f, 4000.0f, 0.0f, 1.0f, Heads {}));
    OneByOne s { d };
    test::run (s.step (0.0f, 0.0f));
    d.setHead (0, { 1.0f, 1.0f, 0.0f });
    double worstStep = 0.0;
    float last = d.headGain (0, 0);
    for (int m = 1; m <= 300; ++m)
    {
        test::run (s.step (0.0f, 0.0f));
        worstStep = std::max (worstStep, std::fabs ((double) d.headGain (0, 0) - (double) last));
        last = d.headGain (0, 0);
    }
    const bool midway = last > 0.25f && last < 0.35f;
    d.setHead (0, { 1.0f, 0.2f, 0.0f });
    int landedAt = -1;
    for (int m = 1; m <= 2 * len; ++m)
    {
        test::run (s.step (0.0f, 0.0f));
        const float g = d.headGain (0, 0);
        worstStep = std::max (worstStep, std::fabs ((double) g - (double) last));
        last = g;
        if (landedAt < 0 && same (g, 0.2f))
            landedAt = m;
    }
    // Catches a ramp restarted from the old TARGET (the gain would jump from ~0.31 to 1 and ramp down from there).
    ok (midway && worstStep <= 1.0 / (double) len + 1.0e-6,
        "no step between consecutive gains larger than one full-range ramp step 1/N (worst " + std::to_string (worstStep) + ")");
    ok (landedAt == len, "...and it lands on 0.2 exactly N samples after the second send (" + std::to_string (landedAt) + ")");
}

//==============================================================================
// 26 · SIDE 1 — the right side's own ramps, glides and mono hold
//==============================================================================
void sideOne()
{
    group ("side 1: a hard-right head ramps, glides as a sounding head, and a mono call leaves its balance gain in place");
    constexpr double sr = 48000.0;
    const int len = rampLength (sr);
    {
        // Catches an at-rest test that looks only at side 0's ramp: side 1's ramp would never be stepped.
        delay::TapeDelay d;
        test::run (d.prepare (sr, 64, 2));
        d.setParams (knobs (300.0f, 0.3f, 4000.0f, 0.0f, 1.0f, layout ({ { 1.0f, 0.5f, 1.0f } })));
        OneByOne s { d };
        for (int i = 0; i < 4; ++i)
            test::run (s.step (0.0f, 0.0f));
        const float from = panGains (0.5f, 1.0f)[1], to = panGains (0.9f, 1.0f)[1];
        const bool start = same (d.headGain (0, 0), 0.0f) && same (d.headGain (0, 1), from);
        d.setHead (0, { 1.0f, 0.9f, 1.0f });
        int landedAt = -1;
        bool moving = true;
        for (int m = 1; m <= len + 10; ++m)
        {
            test::run (s.step (0.0f, 0.0f));
            const float g = d.headGain (0, 1);
            if (landedAt < 0 && same (g, to))
                landedAt = m;
            if (landedAt < 0)
                moving = moving && g > from && g < to;
        }
        ok (start && moving && landedAt == len, "a hard-right head's level change ramps on side 1 and lands at N ("
                                                    + std::to_string (landedAt) + ")");
    }
    {
        // Catches a silent test that takes "either side at 0" for silence: a hard-right head would snap, not glide.
        delay::TapeDelay d;
        test::run (d.prepare (sr, 64, 2));
        d.setParams (knobs (300.0f, 0.3f, 4000.0f, 0.0f, 1.0f, layout ({ { 0.5f, 1.0f, 1.0f } })));
        OneByOne s { d };
        for (int i = 0; i < 4; ++i)
            test::run (s.step (0.0f, 0.0f));
        d.setHead (0, { 0.25f, 1.0f, 1.0f });
        test::run (s.step (0.0f, 0.0f));
        const double p = d.headPosition (0);
        ok (p < 0.5 && p > 0.49, "a sounding hard-right head sent from 1/2 to 1/4 takes one glide step, not a snap");
    }
    {
        // A mono first call snaps side 0 to the level and side 1 to its BALANCE gain; the stereo call that follows
        // ramps side 0 only and plays side 1 at that balance from its first sample. Catches side 1 aimed at the level
        // in a mono call (it would hold the level and never be sent back). Oracle: the frozen class's right wet.
        constexpr int D = 14400, L = 2000, N = 2 * D;
        const auto g = panGains (0.6f, 0.5f);
        delay::TapeDelay now;
        frz::TapeDelay then;
        test::run (now.prepare (sr, 512, 2));
        test::run (then.prepare (sr, 512, 2));
        now.setParams (knobs (300.0f, 0.4f, 4000.0f, 0.0f, 1.0f, layout ({ { 1.0f, 0.6f, 0.5f } })));
        then.setParams (frozenKnobs (300.0f, 0.4f, 4000.0f, 0.0f, 1.0f));
        Planes monoA = burst (200, 200, 131u, 0.3f, false), monoB = monoA;
        float* ma[2] = { monoA[0].data(), nullptr };
        float* mb[2] = { monoB[0].data(), nullptr };
        test::run (now.process (ma, 1, 200));
        test::run (then.process (mb, 1, 200));
        const bool held = same (now.headGain (0, 0), 0.6f) && same (now.headGain (0, 1), g[1]);
        Planes a = burst (N, L, 132u, 0.5f, false), b = a;
        test::run (renderBlocks (now, a, 512));
        test::run (renderBlocks (then, b, 512));
        int bad = 0;
        for (int t = L; t < N; ++t)
            bad += same (get (a[1], t), 0.0f + scaled (get (b[1], t), g[1])) ? 0 : 1;
        ok (held, "after a mono first call: side 0 at the level, side 1 at its balance gain, exactly");
        ok (bad == 0 && same (now.headGain (0, 1), g[1]),
            "the stereo call after it plays the right line at that balance gain, byte for byte (" + std::to_string (bad) + " differ)");
    }
}

//==============================================================================
// 27 · THE RAMP IS THE PINNED DOUBLE FORMULA
//==============================================================================
void rampFormula()
{
    group ("the ramp is the pinned double formula: (float)(from + [(to - from) x (m x (1/N))]), the product rounded on its own");
    constexpr double sr = 48000.0;
    const int len = rampLength (sr);
    const double step = 1.0 / (double) len;
    const auto pinned = [&] (float from, float to, int m)
    {
        volatile double prod = ((double) to - (double) from) * ((double) m * step);
        const double y = prod + (double) from;
        return (float) y;
    };
    // The first pair was found by a search: at m = 560 of 960 a FUSED multiply-add lands one float away from the
    // pinned formula. The others spread the ramp over directions and sizes.
    const float pairs[][2] = { { 0x1.e66104p-1f, 0x1.0ef32p-5f }, { 0.2f, 0.9f }, { 0.9f, 0.2f }, { 0.123f, 0.987f },
                               { 1.0f, 0.333f }, { 0.0517f, 0.0518f }, { 0.7f, 0.0f } };
    {
        const double x = 560.0 * step, f = (double) pairs[0][0], dd = (double) pairs[0][1] - (double) pairs[0][0];
        ok (! same (pinned (pairs[0][0], pairs[0][1], 560), (float) std::fma (dd, x, f)),
            "the fixture: at m = 560 the pinned and the fused formula land on different floats");
    }
    int bad = 0, checked = 0;
    for (const auto& pr : pairs)
    {
        delay::TapeDelay d;
        test::run (d.prepare (sr, 64, 2));
        d.setParams (knobs (300.0f, 0.3f, 4000.0f, 0.0f, 1.0f, layout ({ { 1.0f, pr[0], 0.0f } })));
        OneByOne s { d };
        test::run (s.step (0.0f, 0.0f));
        bad += same (d.headGain (0, 0), pr[0]) ? 0 : 1;
        d.setHead (0, { 1.0f, pr[1], 0.0f });
        for (int m = 1; m <= len; ++m)
        {
            test::run (s.step (0.0f, 0.0f));
            const float want = m < len ? pinned (pr[0], pr[1], m) : pr[1];
            bad += same (d.headGain (0, 0), want) && same (d.headGain (0, 1), want) ? 0 : 1;
            ++checked;
        }
    }
    // Catches a ramp in float, a ramp whose multiply-add is fused (where the row contracts), another step law.
    ok (bad == 0, "every intermediate gain of " + std::to_string (checked) + " is the pinned formula's float, both sides ("
                      + std::to_string (bad) + " differ)");
}

//==============================================================================
// 28 · THE FAST PATH NEEDS UNITY ON EVERY SIDE — gL exactly 1 at a small pan, gR not
//==============================================================================
void plainBothSides()
{
    group ("the fast path needs unity on every side the call reads: a level at pan -0.001 where gL is exactly 1 and gR is not");
    constexpr double sr = 48000.0;
    constexpr int D = 14400, L = 2000, N = 2 * D + 500;
    const float pan = -0.001f;
    const double ratio = (1.0 - (double) pan) / std::sqrt (1.0 + (double) pan * (double) pan);
    float level = 0.0f;
    bool found = false;
    float c = (float) (1.0 / ratio);
    for (int i = 0; i < 64; ++i)
        c = std::nextafter (c, 0.0f);
    for (int i = 0; i < 128 && ! found; ++i, c = std::nextafter (c, 2.0f))
    {
        const auto g = panGains (c, pan);
        if (same (g[0], 1.0f) && ! same (g[1], 1.0f) && c < 1.0f)
        {
            level = c;
            found = true;
        }
    }
    const auto g = panGains (level, pan);
    ok (found, "found: level " + std::to_string (level) + " at pan -0.001 gives gL == 1.0f exactly, gR " + std::to_string (g[1]));

    const Planes in = burst (N, L, 141u, 0.7f, false);
    frz::TapeDelay then;
    test::run (then.prepare (sr, 512, 2));
    then.setParams (frozenKnobs (300.0f, 0.5f, 3000.0f, 0.0f, 1.0f));
    Planes wetD = in;
    test::run (renderBlocks (then, wetD, 512));

    delay::TapeDelay now;
    test::run (now.prepare (sr, 512, 2));
    now.setParams (knobs (300.0f, 0.5f, 3000.0f, 0.0f, 1.0f, layout ({ { 1.0f, level, pan } })));
    Planes out = in;
    test::run (renderBlocks (now, out, 512));
    int bad = 0;
    for (int t = L; t < N; ++t)
    {
        bad += same (get (out[0], t), get (wetD[0], t)) ? 0 : 1;                       // left: the tap at unity
        bad += same (get (out[1], t), 0.0f + scaled (get (wetD[1], t), g[1])) ? 0 : 1;   // right: gR x the tap
    }
    // Catches a fast path that reads side 0 only: a stereo call would play the right line at unity too.
    ok (same (now.headGain (0, 0), 1.0f) && same (now.headGain (0, 1), g[1]) && bad == 0,
        "the left plays the loop's tap, the right gR x the tap, byte for byte (" + std::to_string (bad) + " differ)");
}

//==============================================================================
// 20 · NO ALLOCATION
//==============================================================================
void noAllocation()
{
    group ("RT: the head setters, process with every head moving (stereo, mono, pause, stand-by) and reset never allocate");
    delay::TapeDelay d;
    test::run (d.prepare (96000.0, 1024, 2));
    std::vector<float> l (1024, 0.25f), r (1024, -0.25f);
    float* lr[2] = { l.data(), r.data() };
    Heads moving {};
    for (int k = 0; k < kHeads; ++k)
        moving[(std::size_t) k] = { 0.2f + 0.2f * (float) k, 0.9f - 0.1f * (float) k, 0.3f * (float) k - 0.4f };

    delay::TapeDelay plainStage;   // the default layout: every call at rest takes the loop's-tap fast path
    test::run (plainStage.prepare (96000.0, 1024, 2));
    plainStage.setParams (knobs (700.0f, 0.6f, 3000.0f, -20.0f, 0.7f));

    const long long before = felitronics::test::alloc::count.load();
    test::run (plainStage.process (lr, 2, 1024));
    test::run (plainStage.process (lr, 2, 1024));
    test::run (plainStage.process (lr, 1, 1024));
    d.setHead (0, { 0.7f, 0.8f, -0.2f });
    d.setHeads (moving);
    d.setParams (knobs (700.0f, 0.6f, 3000.0f, -20.0f, 0.7f, moving));
    test::run (d.process (lr, 2, 1024));
    for (auto& h : moving)
        h.position *= 0.5f;
    d.setHeads (moving);
    test::run (d.process (lr, 1, 1024));
    d.setHead (0, { 0.1f, 0.0f, 1.0f });
    test::run (d.process (nullptr, 0, 1024));
    test::run (d.process (lr, 2, 1024, false));
    d.reset();
    d.reset();
    test::run (d.process (lr, 2, 1024));
    const long long delta = felitronics::test::alloc::count.load() - before;
    // Catches an allocation on any of these paths — the default layout's fast path included.
    ok (delta == 0, "no allocation: the default layout's fast path (stereo, mono), setHead, setHeads, setParams with heads, "
                    "process (stereo, mono, pause, stand-by), reset (" + std::to_string (delta) + ")");
}
//==============================================================================
// 22 · A HAIR SHORT OF THE END, IN A HEAD'S READ — no sanitizer needed
//==============================================================================
void hairShortOfTheEnd()
{
    group ("a head's read a hair short of the line's end is the line's first sample, and its last slot's neighbour is the "
           "first — against the frozen class's line, read outside");
    constexpr double sr = 48000.0;
    constexpr int D = 14400;
    const int lineN = lineLength (sr);

    // The line from outside: the frozen class at 300 ms, REPEATS 0, an impulse at sample 0 — its output at D + s is
    // what the line holds at s.
    Plane line (400, 0.0f);
    {
        frz::TapeDelay then;
        test::run (then.prepare (sr, 4096, 1));
        then.setParams (frozenKnobs (300.0f, 0.0f, 16000.0f, 0.0f, 1.0f));
        Plane io ((std::size_t) D + 400, 0.0f);
        cell (io, 0) = 0.8f;
        float* p[1] = { io.data() };
        test::run (then.process (p, 1, D + 400));
        for (int t = 0; t < 400; ++t)
            cell (line, t) = get (io, D + t);
    }

    for (int which = 0; which < 2; ++which)
    {
        // A position whose distance sits a hair past a whole sample (which 0: the read before the write clock
        // reaches it rounds onto the line's length) or about half one (which 1: the read sits in the last slot).
        float pos = 0.5f, hd = 0.0f;
        bool found = false;
        for (int step = 0; step < 1000000 && ! found; ++step, pos = std::nextafter (pos, 1.0f))
        {
            hd = (float) ((double) pos * (double) D);
            const float frac = hd - std::floor (hd);
            found = which == 0 ? (frac > 0.004f && frac < 0.012f) : (frac > 0.45f && frac < 0.55f);
            if (found)
                break;
        }
        ok (found, std::string ("the fixture: a position whose distance sits ") + (which == 0 ? "a hair" : "half a sample")
                       + " past a whole sample was found");
        if (! found)
            continue;
        const int I = (int) std::floor (hd);

        delay::TapeDelay d;
        test::run (d.prepare (sr, 4096, 1));
        d.setParams (knobs (300.0f, 0.0f, 16000.0f, 0.0f, 1.0f, layout ({ { pos, 1.0f, 0.0f } })));
        Plane io ((std::size_t) I + 8, 0.0f);
        cell (io, 0) = 0.8f;
        float* q[1] = { io.data() };
        test::run (d.process (q, 1, I + 8));

        // The read at each sample around the head, outside: the class's float index arithmetic over the line as
        // written up to that sample (nothing at or after it yet).
        double worst = 0.0, reach = 0.0;
        bool reached = false;
        for (int t = I - 2; t < I + 6; ++t)
        {
            float hp = (float) t - hd;
            if (hp < 0.0f)
                hp += (float) lineN;
            const int   whole = (int) hp;
            const float f = hp - (float) whole;
            const int   i0 = whole == lineN ? 0 : whole, i1 = i0 + 1 == lineN ? 0 : i0 + 1;
            reached = reached || (which == 0 ? whole == lineN : whole == lineN - 1);
            const auto at = [&] (int i) { return i < t && i < 400 ? (double) get (line, i) : 0.0; };
            const double want = at (i0) * (1.0 - (double) f) + at (i1) * (double) f;
            worst = std::max (worst, std::fabs ((double) get (io, t) - want));
            reach = std::max (reach, std::fabs (want));
        }
        // Catches: the head's read without its hair-short rule (index lineN, one past the end) or without its second
        // index's wrap (lineN again) — each reads past the line; the CI sanitizer row saw them, no other check did.
        ok (reached && worst <= 1.0e-6 && reach > 0.05,
            std::string (which == 0 ? "a read a hair short of the end" : "a read in the line's last slot")
                + ": the line's own samples, eight around the head (worst " + std::to_string (worst) + ")");
    }
}

} // namespace

int main()
{
    std::printf ("felitronics::delay::TapeDelay — the playback heads (%d)\n", kHeads);
    placement();
    stackedHeads();
    tripletFigure();
    pingPong();
    panLaw();
    panAnchors();
    ramps();
    positionGlide();
    bend();
    silentSnapsSoundingGlides();
    monoIsARamp();
    snapAfterRestart();
    reslicing();
    pauseIsSilence();
    zeroSoundingHeads();
    awayAndBack();
    law8();
    restartedRampFloor();
    poisonRestart();
    tail();
    validation();
    underOneSample();
    shownHeads();
    wakeSnap();
    hairShortOfTheEnd();
    plainPathReslicing();
    orderedSum();
    exactGuard();
    retargetFromWhereItStands();
    sideOne();
    rampFormula();
    plainBothSides();
    noAllocation();
    return felitronics::test::report();
}

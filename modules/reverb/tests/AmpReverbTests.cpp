// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// reverb::Freeverb and reverb::AmpReverb — regression goldens (the core's own renders; how far they sit from
// the app's original class is the out-of-tree parity harness's report), Freeverb's published constants and
// glide, and the stage contract: prepare binds and a refused one disarms, a fresh room starts at its
// parameters, the tail is right from construction, reset is a restart, absurd knobs stay bounded, chunking,
// a stopped channel, poison, law 8, no allocation. The pause (law 11c) is in modules/laws.
#include <felitronics_test.h>
#include <alloc_counter.h>   // installs the allocation counter

#include "reverb_cases.h"
#include "reverb_golden.h"

#include <cmath>
#include <limits>
#include <string>
#include <vector>

using namespace felitronics;
using felitronics::test::ok;
using felitronics::test::group;
namespace rc = felitronics::test::reverbcases;

namespace
{
double noSin (double) { return 0.0; }
double noExp (double) { return 0.0; }

// The goldens are bit-exact on the row that wrote them. Elsewhere only FMA contraction can move a bit (law 10:
// arm64 fuses a*b + c, x86-64 and wasm cannot) — every transcendental is core::det. One path AMPLIFIES that
// last bit: the MODULATED chorus reads its line at a fractional position of ~576 samples built from fusable
// multiply-adds, so one ulp of position (6e-5 of a sample) times the tail's slope reaches the output —
// measured 1.0e-5 between arm64 and wasm with core::det::sin in the LFO (8e-6 with the system sin: the sine
// was never the cause). A scenario that ever runs it is held to 5e-5; making it bit-identical across rows
// would take pinning every multiply-add in the chain against contraction, which law 10 reserves to one
// header today.
constexpr double kSampleTol = 1.0e-6;
constexpr double kChorusTol = 5.0e-5;
constexpr double kSumRelTol = 1.0e-6;

bool runsChorus (const rc::Case& c)
{
    bool any = c.knobs.character == (int) reverb::Character::modulated;
    for (const auto& e : c.events) any = any || e.knobs.character == (int) reverb::Character::modulated;
    return any;
}

void goldens()
{
    group ("regression goldens: the core's own renders");
    const auto cases = rc::goldenCases();
    ok (std::size (test::reverbgolden::kGolden) == cases.size(), "one golden per scenario");

    for (std::size_t k = 0; k < cases.size(); ++k)
    {
        const auto& c = cases[k];
        const auto& g = test::reverbgolden::kGolden[k];
        ok (std::string (g.name) == c.name, std::string ("golden order: ") + c.name);

        auto io = rc::makeInput (c, noSin, noExp);
        rc::CoreAdapter ad;
        rc::render (c, ad, io);
        ok (ad.ok, std::string (c.name) + ": every call accepted");

        double worst = 0.0;
        for (const auto& w : g.w)
            for (int ch = 0; ch < c.nch; ++ch)
                for (int i = 0; i < 16; ++i)
                    worst = std::max (worst, std::abs ((double) io[(std::size_t) ch][(std::size_t) (w.start + i)] - (double) w.x[ch][i]));
        const double tol = runsChorus (c) ? kChorusTol : kSampleTol;
        ok (worst <= tol, std::string (c.name) + ": 3 windows within " + (runsChorus (c) ? "5e-5" : "1e-6")
                          + " (worst " + std::to_string (worst) + ")");

        for (int ch = 0; ch < c.nch; ++ch)
        {
            double sa = 0.0, sq = 0.0;
            for (float v : io[(std::size_t) ch]) { sa += std::abs ((double) v); sq += (double) v * v; }
            ok (std::abs (sa - g.sumAbs[ch]) <= kSumRelTol * g.sumAbs[ch] + 1e-9
                && std::abs (sq - g.sumSq[ch]) <= kSumRelTol * g.sumSq[ch] + 1e-9,
                std::string (c.name) + ": channel " + std::to_string (ch) + " sums of |x| and x^2 match");
        }
    }
}

void freeverbConstants()
{
    group ("Freeverb: the published tunings, scaled in integer arithmetic, the 10 ms glide, refusals");
    reverb::Freeverb f;
    float x1[4] {}, x2[4] {};
    ok (f.longestCombSamples() == 0, "nothing allocated before prepare");
    ok (! f.processStereo (x1, x2, 4) && ! f.processMono (x1, 4), "an unprepared Freeverb refuses, and says so");
    ok (! f.prepare (0.0) && ! f.prepare (999.0) && ! f.prepare (768001.0) && ! f.prepare (std::numeric_limits<double>::infinity()),
        "prepare refuses a rate outside 1 kHz … 768 kHz");
    test::run (f.prepare (44100.0));
    ok (f.longestCombSamples() == 1617 + 23, "44.1 kHz: the right channel's longest comb is 1617 + 23");
    ok (! f.processMono (x1, -1) && ! f.processStereo (x1, nullptr, 4), "a negative length or a null plane is refused");
    test::run (f.prepare (96000.0));
    ok (f.longestCombSamples() == 1640 * 96000 / 44100, "96 kHz: (int) rate * tuning / 44100 = 3570");
    ok (! f.prepare (-1.0) && ! f.processMono (x1, 4), "a refused re-prepare disarms");

    // The dry glide, observed: dry 0 -> 1 (gain 2) over floor(0.01 * 48000) = 480 samples, wet off, freeze on
    // so nothing but the dry reaches the output.
    reverb::Freeverb g;
    reverb::FreeverbParams p; p.wetLevel = 0.0f; p.dryLevel = 0.0f; p.freezeMode = 1.0f;
    g.setParams (p);
    test::run (g.prepare (48000.0));
    p.dryLevel = 1.0f;
    g.setParams (p);
    std::vector<float> x (600, 1.0f);
    test::run (g.processMono (x.data(), 600));
    ok (x[478] < 2.0f && x[479] == 2.0f && x[599] == 2.0f, "the dry lands on 2 x dryLevel at exactly the 480th sample");
    test::approx (x[239], 1.0, 1e-5, "...linearly: half-way at the 240th");

    group ("Freeverb: freeze holds the tail and takes no input");
    reverb::Freeverb z;
    test::run (z.prepare (48000.0));
    reverb::FreeverbParams q; q.wetLevel = 1.0f / 3.0f; q.dryLevel = 0.0f; q.roomSize = 0.5f;
    z.setParams (q);
    std::vector<float> l (4800, 0.0f), r (4800, 0.0f);
    rc::Lcg n { 5u };
    for (int i = 0; i < 2400; ++i) { l[(std::size_t) i] = 0.5f * n.next(); r[(std::size_t) i] = 0.5f * n.next(); }
    test::run (z.processStereo (l.data(), r.data(), 4800));
    q.freezeMode = 1.0f;
    z.setParams (q);
    auto energy = [] (const std::vector<float>& v) { double e = 0; for (float s : v) e += (double) s * s; return e; };
    std::vector<float> a (48000, 0.3f), b (48000, -0.3f);   // loud input: frozen, it must not get in
    test::run (z.processStereo (a.data(), b.data(), 48000));
    std::vector<float> a2 (a.end() - 4800, a.end()), a1 (a.begin() + 4800, a.begin() + 9600);
    const double e1 = energy (a1), e2 = energy (a2);
    // Unfrozen, this room would lose ~99.8 % of its energy in that second; a leak of the loud input would
    // multiply it. A frozen tail only beats (eight comb periods against a window): within 20 %.
    ok (e1 > 0.0 && e2 / e1 > 0.8 && e2 / e1 < 1.25, "a frozen tail neither grows nor decays over a second ("
                                                     + std::to_string (e1) + " -> " + std::to_string (e2) + ")");
}

void contract()
{
    group ("contract: prepare is binding, a refused one disarms, refusals move nothing, n == 0 is a no-op");
    reverb::AmpReverb r;
    float a[64] {}, b[64] {};
    float* io[2] = { a, b };
    ok (! r.process (io, 1, 64), "an unprepared stage refuses");
    ok (! r.prepare (48000.0, 64, 3) && ! r.prepare (48000.0, 64, 0) && ! r.prepare (0.0, 64, 2) && ! r.prepare (999.0, 64, 2)
        && ! r.prepare (768001.0, 64, 2) && ! r.prepare (48000.0, 0, 2),
        "prepare refuses a width outside 1..2, a rate outside 1 kHz … 768 kHz and an empty block");
    test::run (r.prepare (48000.0, 64, 1));
    a[0] = 0.5f;
    ok (! r.process (io, 2, 64) && a[0] == 0.5f, "a call wider than prepared is refused before anything moves");
    ok (r.process (io, 1, 0) && a[0] == 0.5f, "n == 0 is accepted and changes nothing");
    ok (! r.prepare (48000.0, 64, 4), "a re-prepare at a width it cannot honour is refused");
    a[0] = 0.25f;
    ok (! r.isPrepared() && ! r.process (io, 1, 64) && a[0] == 0.25f,
        "...and DISARMS: the old preparation no longer stands, process refuses and touches nothing");
    ok (r.prepare (48000.0, 1 << 30, 2) && r.maxBlock() == reverb::AmpReverb::kMaxChunk,
        "a whole-file maxBlock is accepted and costs only the capped scratch");

    group ("law 11a: a call past maxBlock is processed whole, as the caller chunking it would be");
    rc::Knobs k; k.character = 5; k.mix = 0.7f; k.predelayMs = 20.0f;
    auto make = [] { std::vector<std::vector<float>> v (2, std::vector<float> (5000)); rc::Lcg g { 9u };
                     for (auto& ch : v) for (auto& s : ch) s = 0.5f * g.next(); return v; };
    auto one = make(), many = make();
    rc::CoreAdapter x, y;
    x.prepare (48000.0, 64, 2); y.prepare (48000.0, 64, 2);
    x.apply (k); y.apply (k);
    float* p1[2] = { one[0].data(), one[1].data() };
    x.process (p1, 2, 5000, true);
    for (int pos = 0; pos < 5000; pos += 64)
    {
        float* p2[2] = { many[0].data() + pos, many[1].data() + pos };
        y.process (p2, 2, std::min (64, 5000 - pos), true);
    }
    ok (x.ok && y.ok && one == many, "one 5000-sample call == the caller's own 64-sample calls, bit for bit");
}

void startsAtParams()
{
    group ("a fresh room starts AT its parameters: nothing of the dry leaks into the tail");
    for (int c = 0; c < reverb::kNumCharacters; ++c)
    {
        reverb::AmpReverb r;
        reverb::AmpReverbParams p; p.character = (reverb::Character) c; p.mix = 1.0f;
        r.setParams (p);
        test::run (r.prepare (48000.0, 1000, 2));
        std::vector<float> l (1000), rr (1000);
        rc::Lcg g { 6u };
        for (int i = 0; i < 1000; ++i) { l[(std::size_t) i] = 0.5f * g.next(); rr[(std::size_t) i] = 0.5f * g.next(); }
        const auto l0 = l, r0 = rr;
        float* lr[2] = { l.data(), rr.data() };
        test::run (r.process (lr, 2, 1000));
        // The first comb returns after 1214 samples at 48 kHz: until then the room adds exactly nothing.
        ok (l == l0 && rr == r0, std::string (reverb::kCharacterNames[c]) + ": the first 1000 samples are the dry, bit for bit");
    }
}

double tailFormula (float roomSize, double predelayMs)
{
    const double g = (double) roomSize * 0.28 + 0.7;
    return predelayMs * 0.001 + 1640.0 / 44100.0 * (1.0 + std::log (0.001) / std::log (g));
}

void tail()
{
    group ("tailSeconds: right from construction, the longest comb's thousandth plus the predelay");
    reverb::AmpReverb fresh;
    test::approx (fresh.tailSeconds(), tailFormula (0.35f, 0.0), 1e-5, "a fresh ROOM already reports the formula (1.176 s)");
    test::run (fresh.prepare (48000.0, 64, 2));
    test::approx (fresh.tailSeconds(), tailFormula (0.35f, 0.0), 1e-5, "...and after a prepare at its default knobs");

    reverb::AmpReverb r;
    reverb::AmpReverbParams p; p.character = reverb::Character::hall; p.decay = 2.0f; p.predelayMs = 50.0f;
    r.setParams (p);
    test::approx (r.tailSeconds(), tailFormula (0.98f, 50.0), 1e-4, "HALL at DECAY x2: the room clamped to 0.98");
    ok (r.tailSeconds() > 9.0f && r.tailSeconds() < 10.5f, "...a shade under ten seconds: " + std::to_string (r.tailSeconds()));
}

void resetIsRestart()
{
    group ("reset() is a stream restart: what follows does not depend on what came before");
    std::vector<std::vector<float>> first;
    bool same = true;
    for (int warm : { 0, 1, 333, 4801, 96001 })
    {
        reverb::AmpReverb r;
        reverb::AmpReverbParams p; p.character = reverb::Character::modulated; p.mix = 1.0f; p.predelayMs = 17.0f;
        r.setParams (p);
        test::run (r.prepare (48000.0, 512, 2));
        std::vector<float> l (512), rr (512);
        float* lr[2] = { l.data(), rr.data() };
        rc::Lcg g { 31u };
        for (int done = 0; done < warm; done += 512)
        {
            for (auto& v : l) v = 0.3f * g.next();
            for (auto& v : rr) v = 0.3f * g.next();
            test::run (r.process (lr, 2, std::min (512, warm - done)));
        }
        p.decay = 1.7f;                  // a glide set just before the restart: still in flight at reset()
        r.setParams (p);
        if (warm > 0) test::run (r.process (lr, 2, 7));
        r.reset();
        std::vector<std::vector<float>> out;
        rc::Lcg h { 77u };
        for (int b = 0; b < 40; ++b)
        {
            for (auto& v : l) v = 0.3f * h.next();
            for (auto& v : rr) v = 0.3f * h.next();
            test::run (r.process (lr, 2, 512));
            if (out.empty()) out.assign (2, {});
            out[0].insert (out[0].end(), l.begin(), l.end());
            out[1].insert (out[1].end(), rr.begin(), rr.end());
        }
        if (first.empty()) first = out; else same = same && out == first;
    }
    ok (same, "after warm-ups of 0, 1, 333, 4801 and 96001 samples (a glide in flight), the post-reset output is bit-identical");
}

void extremes()
{
    group ("absurd but finite knobs stay bounded");
    bool finite = true;
    const float F = std::numeric_limits<float>::max();
    for (double sr : { 1000.0, 48000.0, 768000.0 })
    {
        reverb::AmpReverb r;
        test::run (r.prepare (sr, 256, 2));
        std::vector<float> l (256), rr (256);
        float* lr[2] = { l.data(), rr.data() };
        rc::Lcg g { 2u };
        const float knobs[][4] = { { F, F, F, F }, { -F, -F, -F, -F }, { 1e-30f, 1e-30f, 1e-30f, 1e-30f }, { 0.5f, 1.0f, 30.0f, 120.0f } };
        for (int c = 0; c < reverb::kNumCharacters; ++c)
            for (const auto& k : knobs)
            {
                r.setParams ({ (reverb::Character) c, k[0], k[1], k[2], k[3] });
                for (int b = 0; b < 8; ++b)
                {
                    for (auto& v : l) v = 0.5f * g.next();
                    for (auto& v : rr) v = 0.5f * g.next();
                    test::run (r.process (lr, 2, 256));
                    for (int i = 0; i < 256; ++i) finite = finite && std::isfinite (l[(std::size_t) i]) && std::isfinite (rr[(std::size_t) i]);
                }
            }
        r.setCharacter ((reverb::Character) 99);
        finite = finite && std::isfinite (r.tailSeconds());
    }
    ok (finite, "±FLT_MAX and 1e-30 on every knob of every character, at 1, 48 and 768 kHz: every sample finite");
}

void stoppedChannel()
{
    group ("law 11a: a channel that stopped comes back to silence, not to its old room");
    reverb::AmpReverb r;
    test::run (r.prepare (48000.0, 256, 2));
    reverb::AmpReverbParams p; p.character = reverb::Character::hall; p.mix = 1.0f; p.decay = 2.0f;   // width 1: no cross-feed
    r.setParams (p);
    std::vector<float> l (256), rr (256);
    float* lr[2] = { l.data(), rr.data() };
    rc::Lcg g { 4u };
    for (int b = 0; b < 20; ++b)
    {
        for (int i = 0; i < 256; ++i) { l[(std::size_t) i] = 0.3f * g.next(); rr[(std::size_t) i] = 0.3f * g.next(); }
        test::run (r.process (lr, 2, 256));
    }
    std::fill (l.begin(), l.end(), 0.0f);
    test::run (r.process (lr, 1, 256));   // the right channel STOPS
    float worstR = 0.0f, lastL = 0.0f;
    for (int b = 0; b < 40; ++b)
    {
        std::fill (l.begin(), l.end(), 0.0f); std::fill (rr.begin(), rr.end(), 0.0f);
        test::run (r.process (lr, 2, 256));
        for (float v : rr) worstR = std::max (worstR, std::abs (v));
        lastL = std::abs (l[255]);
    }
    ok (worstR == 0.0f, "the right channel returns to exact silence");
    ok (lastL > 0.0f, "...while the left, which never stopped, still rings");
}

void poisonAndLaw8()
{
    group ("poison: a bad input sample stays the dry's; the room never hears it");
    reverb::AmpReverb r;
    test::run (r.prepare (48000.0, 128, 2));
    reverb::AmpReverbParams p; p.mix = 1.0f; p.character = reverb::Character::plate;
    r.setParams (p);
    std::vector<float> l (128, 0.1f), rr (128, 0.1f);
    l[10] = std::numeric_limits<float>::quiet_NaN();
    rr[20] = -std::numeric_limits<float>::infinity();
    float* lr[2] = { l.data(), rr.data() };
    test::run (r.process (lr, 2, 128));
    bool finiteElsewhere = true;
    for (int i = 0; i < 128; ++i)
    {
        if (i != 10) finiteElsewhere = finiteElsewhere && std::isfinite (l[(std::size_t) i]);
        if (i != 20) finiteElsewhere = finiteElsewhere && std::isfinite (rr[(std::size_t) i]);
    }
    ok (std::isnan (l[10]) && std::isinf (rr[20]) && finiteElsewhere, "the NaN and the inf pass as dry; every other sample is finite");

    group ("law 8: after a long HALL, no subnormal anywhere, and the silence, once reached, holds");
    reverb::AmpReverb s;
    test::run (s.prepare (48000.0, 480, 2));
    reverb::AmpReverbParams q; q.character = reverb::Character::hall; q.mix = 1.0f; q.decay = 2.0f; q.hpfHz = 40.0f;
    s.setParams (q);
    std::vector<float> a (480), b (480);
    float* ab[2] = { a.data(), b.data() };
    rc::Lcg g { 8u };
    for (int i = 0; i < 480; ++i) { a[(std::size_t) i] = 0.5f * g.next(); b[(std::size_t) i] = 0.5f * g.next(); }
    test::run (s.process (ab, 2, 480));
    bool anySub = false;
    int lastLoud = -1;
    for (int blk = 0; blk < 18000; ++blk)   // three minutes
    {
        std::fill (a.begin(), a.end(), 0.0f); std::fill (b.begin(), b.end(), 0.0f);
        test::run (s.process (ab, 2, 480));
        anySub = anySub || s.anySubnormalState();
        const bool silent = std::all_of (a.begin(), a.end(), [] (float v) { return v == 0.0f; })
                         && std::all_of (b.begin(), b.end(), [] (float v) { return v == 0.0f; });
        if (! silent) lastLoud = blk;
    }
    ok (! anySub, "no subnormal in the room, the predelay or the chorus line, at any block");
    ok (lastLoud > 100 && lastLoud < 15000, "the tail rang (" + std::to_string (lastLoud / 100.0)
                                            + " s), then reached exact zero and STAYED there to the end");
}

void noAllocation()
{
    group ("RT: process, reset and the setters never allocate");
    reverb::AmpReverb r;
    test::run (r.prepare (96000.0, 256, 2));
    std::vector<float> l (1024, 0.25f), rr (1024, -0.25f);
    float* lr[2] = { l.data(), rr.data() };
    const long long before = alloc::count.load();
    reverb::AmpReverbParams p; p.character = reverb::Character::modulated; p.predelayMs = 60.0f; p.decay = 1.7f;
    r.setParams (p);
    test::run (r.process (lr, 2, 1024));
    test::run (r.process (lr, 1, 1024));
    test::run (r.process (nullptr, 0, 1024));
    r.reset(); r.reset();
    test::run (r.process (lr, 2, 1024, false));
    felitronics::test::okNoAlloc (alloc::count.load() == before, "no allocation: setParams, process (stereo, mono, pause, stand-by, past maxBlock), reset");
}
} // namespace

int main()
{
    std::printf ("felitronics::reverb::Freeverb + AmpReverb\n");
    goldens();
    freeverbConstants();
    contract();
    startsAtParams();
    tail();
    resetIsRestart();
    extremes();
    stoppedChannel();
    poisonAndLaw8();
    noAllocation();
    return felitronics::test::report();
}

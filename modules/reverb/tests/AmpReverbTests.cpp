// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// reverb::Freeverb and reverb::AmpReverb — the goldens the app's ORIGINAL stage rendered (the parity harness
// that wrote them drives the app's class, on the framework's Freeverb, itself), Freeverb's published constants
// and glide, and the stage contract: refusals, chunking, a stopped channel, poison, law 8, no allocation.
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

// Tolerances that hold on every row: bit-exact on the row that wrote the goldens; elsewhere std::tan (the HPF
// design) and FMA contraction move the last bits. The MODULATED character reads its chorus line at a delay
// steered by std::sin per sample, and one ulp of that sine is ~1e-5 of a sample of delay: measured 8e-6 between
// Apple's libm and emscripten's, so a scenario that ever runs it is held to 5e-5 instead.
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
    group ("goldens: the app's original ReverbStage, rendered by the parity harness");
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
                          + " of the app (worst " + std::to_string (worst) + ")");

        for (int ch = 0; ch < c.nch; ++ch)
        {
            double sa = 0.0, sq = 0.0;
            for (float v : io[(std::size_t) ch]) { sa += std::abs ((double) v); sq += (double) v * v; }
            ok (std::abs (sa - g.sumAbs[ch]) <= kSumRelTol * g.sumAbs[ch] + 1e-9
                && std::abs (sq - g.sumSq[ch]) <= kSumRelTol * g.sumSq[ch] + 1e-9,
                std::string (c.name) + ": channel " + std::to_string (ch) + " sums of |x| and x^2 match the app");
        }
    }
}

void freeverbConstants()
{
    group ("Freeverb: the published tunings, scaled in integer arithmetic, and the 10 ms glide");
    reverb::Freeverb f;
    ok (f.longestCombSamples() == 0, "nothing allocated before prepare");
    ok (! f.prepare (0.0) && ! f.prepare (-1.0) && ! f.prepare (std::numeric_limits<double>::infinity()), "prepare refuses a bad rate");
    test::run (f.prepare (44100.0));
    ok (f.longestCombSamples() == 1617 + 23, "44.1 kHz: the right channel's longest comb is 1617 + 23");
    test::run (f.prepare (96000.0));
    ok (f.longestCombSamples() == 1640 * 96000 / 44100, "96 kHz: (int) rate * tuning / 44100 = 3570");

    // The dry glide, observed: dry 0 -> 1 (gain 2) over floor(0.01 * 48000) = 480 samples, wet off, freeze on
    // so nothing but the dry reaches the output.
    reverb::Freeverb g;
    reverb::FreeverbParams p; p.wetLevel = 0.0f; p.dryLevel = 0.0f; p.freezeMode = 1.0f;
    g.setParams (p);
    test::run (g.prepare (48000.0));
    p.dryLevel = 1.0f;
    g.setParams (p);
    std::vector<float> x (600, 1.0f);
    g.processMono (x.data(), 600);
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
    z.processStereo (l.data(), r.data(), 4800);
    q.freezeMode = 1.0f;
    z.setParams (q);
    auto energy = [] (const std::vector<float>& v) { double e = 0; for (float s : v) e += (double) s * s; return e; };
    std::vector<float> a (48000, 0.3f), b (48000, -0.3f);   // loud input: frozen, it must not get in
    z.processStereo (a.data(), b.data(), 48000);
    std::vector<float> a2 (a.end() - 4800, a.end()), a1 (a.begin() + 4800, a.begin() + 9600);
    const double e1 = energy (a1), e2 = energy (a2);
    // Unfrozen, this room would lose ~99.8 % of its energy in that second; a leak of the loud input would
    // multiply it. A frozen tail only beats (eight comb periods against a window): within 20 %.
    ok (e1 > 0.0 && e2 / e1 > 0.8 && e2 / e1 < 1.25, "a frozen tail neither grows nor decays over a second ("
                                                     + std::to_string (e1) + " -> " + std::to_string (e2) + ")");
}

void contract()
{
    group ("contract: prepare is binding, refusals move nothing, n == 0 is a no-op");
    reverb::AmpReverb r;
    float a[64] {}, b[64] {};
    float* io[2] = { a, b };
    ok (! r.process (io, 1, 64), "an unprepared stage refuses");
    ok (! r.prepare (48000.0, 64, 3) && ! r.prepare (48000.0, 64, 0) && ! r.prepare (0.0, 64, 2)
        && ! r.prepare (48000.0, 0, 2), "prepare refuses a width outside 1..2, a bad rate and an empty block");
    test::run (r.prepare (48000.0, 64, 1));
    a[0] = 0.5f;
    ok (! r.process (io, 2, 64) && a[0] == 0.5f, "a call wider than prepared is refused before anything moves");
    ok (r.process (io, 1, 0) && a[0] == 0.5f, "n == 0 is accepted and changes nothing");

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

    group ("law 8: a long silence after a long HALL ends in exact digital silence");
    reverb::AmpReverb s;
    test::run (s.prepare (48000.0, 480, 2));
    reverb::AmpReverbParams q; q.character = reverb::Character::hall; q.mix = 1.0f; q.decay = 2.0f; q.hpfHz = 40.0f;
    s.setParams (q);
    std::vector<float> a (480), b (480);
    float* ab[2] = { a.data(), b.data() };
    rc::Lcg g { 8u };
    for (int i = 0; i < 480; ++i) { a[(std::size_t) i] = 0.5f * g.next(); b[(std::size_t) i] = 0.5f * g.next(); }
    test::run (s.process (ab, 2, 480));
    int blocksToSilence = -1;
    for (int blk = 0; blk < 12000 && blocksToSilence < 0; ++blk)   // up to two minutes
    {
        std::fill (a.begin(), a.end(), 0.0f); std::fill (b.begin(), b.end(), 0.0f);
        test::run (s.process (ab, 2, 480));
        const bool silent = std::all_of (a.begin(), a.end(), [] (float v) { return v == 0.0f; })
                         && std::all_of (b.begin(), b.end(), [] (float v) { return v == 0.0f; });
        if (silent) blocksToSilence = blk;
    }
    ok (blocksToSilence > 0, "the tail reaches exact zero (" + std::to_string (blocksToSilence / 100.0) + " s), never parking in the subnormal band");
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

void tail()
{
    group ("tailSeconds: the longest comb's thousandth, plus the predelay");
    reverb::AmpReverb r;
    reverb::AmpReverbParams p; p.character = reverb::Character::hall; p.decay = 2.0f; p.predelayMs = 50.0f;
    r.setParams (p);
    const double g = (double) 0.98f * 0.28 + 0.7;
    test::approx (r.tailSeconds(), 0.05 + 1640.0 / 44100.0 * (1.0 + std::log (0.001) / std::log (g)), 1e-4,
                  "HALL at DECAY x2: the room clamped to 0.98 — a shade under ten seconds");
    ok (r.tailSeconds() > 9.0f && r.tailSeconds() < 10.5f, "...in seconds: " + std::to_string (r.tailSeconds()));
}
} // namespace

int main()
{
    std::printf ("felitronics::reverb::Freeverb + AmpReverb\n");
    goldens();
    freeverbConstants();
    contract();
    stoppedChannel();
    poisonAndLaw8();
    noAllocation();
    tail();
    return felitronics::test::report();
}

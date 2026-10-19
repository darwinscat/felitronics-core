// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// delay::TapeDelay and its knob law — regression goldens (the core's own renders; how far they sit from the
// app's original class is the out-of-tree parity harness's report), the ladder and the ranges, and the stage
// contract: prepare binds and a refused one disarms, the tail is right from construction, reset is a restart,
// absurd knobs stay bounded, re-slicing, a stopped channel, a pause, poison, law 8, the glides land, no
// allocation.
#include <felitronics_test.h>
#include <alloc_counter.h>   // installs the allocation counter

#include "delay_cases.h"
#include "delay_golden.h"

#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

using namespace felitronics;
using felitronics::test::ok;
using felitronics::test::group;
namespace dc = felitronics::test::delaycases;

namespace
{
double noSin (double) { return 0.0; }   // the golden inputs never reach the sweep
double noExp (double) { return 0.0; }

std::vector<std::vector<float>> renderCore (const dc::Case& c, bool* accepted = nullptr)
{
    auto io = dc::makeInput (c, noSin, noExp);
    dc::CoreAdapter ad;
    dc::render (c, ad, io);
    if (accepted) *accepted = ad.ok;
    return io;
}

// The goldens are bit-exact on the row that wrote them. Elsewhere only FMA contraction can move a bit (law 10:
// arm64 fuses a*b + c, x86-64 and wasm cannot) — every transcendental is core::det.
constexpr double kSampleTol = 1.0e-6;
constexpr double kSumRelTol = 1.0e-6;

void goldens()
{
    group ("regression goldens: the core's own renders");
    const auto cases = dc::goldenCases();
    ok (std::size (test::delaygolden::kGolden) == cases.size(), "one golden per scenario");

    for (std::size_t k = 0; k < cases.size(); ++k)
    {
        const auto& c = cases[k];
        const auto& g = test::delaygolden::kGolden[k];
        ok (std::string (g.name) == c.name, std::string ("golden order: ") + c.name);

        bool accepted = false;
        const auto io = renderCore (c, &accepted);
        ok (accepted, std::string (c.name) + ": every call accepted");

        double worst = 0.0;
        for (const auto& w : g.w)
            for (int ch = 0; ch < c.nch; ++ch)
                for (int i = 0; i < 16; ++i)
                    worst = std::max (worst, std::abs ((double) io[(std::size_t) ch][(std::size_t) (w.start + i)] - (double) w.x[ch][i]));
        ok (worst <= kSampleTol, std::string (c.name) + ": 3 windows within 1e-6 (worst " + std::to_string (worst) + ")");

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

void parameterLaw()
{
    group ("the knob law: the whole ladder, the ranges, the sync policy, DARK");
    ok (delay::kNumDivisions == 13 && std::size (delay::kDivisionNames) == 13 && std::size (delay::kDivisionBeats) == 13,
        "13 divisions, 13 names, 13 beat values");
    const char* plain[] = { "1/1", "1/2", "1/4", "1/8", "1/16" };
    const float plainBeats[] = { 4.0f, 2.0f, 1.0f, 0.5f, 0.25f };
    const int   plainAt[] = { 0, 2, 5, 8, 11 };
    bool ladder = true;
    for (int k = 0; k < 5; ++k)
    {
        const int i = plainAt[k];
        ladder = ladder && std::strcmp (delay::kDivisionNames[i], plain[k]) == 0 && delay::kDivisionBeats[i] == plainBeats[k];
        if (i > 0)   // 1/1 has neither: the ladder opens on it
            ladder = ladder && std::strcmp (delay::kDivisionNames[i + 1], (std::string (plain[k]) + "T").c_str()) == 0
                            && delay::kDivisionBeats[i + 1] == plainBeats[k] * (2.0f / 3.0f);     // triplet: ⅔ after it
        if (i > 0)
            ladder = ladder && std::strcmp (delay::kDivisionNames[i - 1], (std::string (plain[k]) + ".").c_str()) == 0
                            && delay::kDivisionBeats[i - 1] == plainBeats[k] * 1.5f;          // dotted: 1.5 before it
    }
    ok (ladder, "after 1/1, every plain value has its dotted (×1.5) just before it and its triplet (×⅔) just after, named so");
    ok (delay::kDivisionDefault == 5 && std::strcmp (delay::kDivisionNames[delay::kDivisionDefault], "1/4") == 0, "the default is 1/4");
    ok (delay::kBpmMin == 40.0f && delay::kBpmMax == 240.0f && delay::kBpmDefault == 120.0f, "BPM field 40 … 240, default 120");
    ok (delay::kTimeMinMs == 20.0f && delay::kTimeMaxMs == 2000.0f, "free TIME 20 … 2000 ms");
    ok (delay::kOffsetMaxMs == 30.0f, "OFFSET reaches ±30 ms");
    ok (delay::kDarkLoHz == 800.0f && delay::kDarkHiHz == 16000.0f && delay::kDarkDefaultPct == 45.0f, "DARK 16 kHz … 800 Hz, default 45 %");

    ok (delay::divisionMs (delay::kDivisionDefault, 120.0) == 500.0f, "1/4 at 120 BPM is 500 ms");
    ok (delay::divisionMs (0, 40.0) == 6000.0f, "1/1 at 40 BPM is 6000 ms — the line's design maximum");
    ok (delay::divisionMs (0, 20.0) == 6000.0f && delay::divisionMs (0, 0.0) == 6000.0f
        && delay::divisionMs (0, std::numeric_limits<double>::quiet_NaN()) == 6000.0f,
        "a tempo below 40 BPM (or none) counts as 40: the knob never names more than the line plays");
    ok (delay::divisionMs (99, 120.0) == delay::divisionMs (12, 120.0) && delay::divisionMs (-3, 120.0) == delay::divisionMs (0, 120.0),
        "an index off the ladder is clamped onto it");
    ok (delay::divisionMs (0, 40.0) <= delay::TapeDelay::kLineMaxMs, "the slowest division fits the line");
    ok (delay::timeMsFor (true, 8, 120.0, 777.0f) == 250.0f && delay::timeMsFor (false, 8, 120.0, 777.0f) == 777.0f,
        "the sync policy: synced plays the division, free plays the knob");
    ok (delay::darkHzForPercent (0.0f) == 16000.0f, "DARK 0 % is the bright ceiling, 16 kHz");
    test::approx (delay::darkHzForPercent (100.0f), 800.0, 0.01, "DARK 100 % is the dark floor, 800 Hz");
    ok (delay::darkHzForPercent (45.0f) > 3500.0f && delay::darkHzForPercent (45.0f) < 4500.0f, "the default 45 % is ~4 kHz");
}

void contract()
{
    group ("contract: prepare is binding, a refused one disarms, refusals move nothing, n == 0 is a no-op");
    delay::TapeDelay d;
    float a[64] {}, b[64] {};
    float* io[2] = { a, b };
    ok (! d.process (io, 1, 64), "an unprepared stage refuses");
    ok (! d.prepare (48000.0, 64, 3) && ! d.prepare (48000.0, 64, 0) && ! d.prepare (0.0, 64, 2)
        && ! d.prepare (999.0, 64, 2) && ! d.prepare (768001.0, 64, 2) && ! d.prepare (1.0e9, 64, 2)
        && ! d.prepare (std::numeric_limits<double>::quiet_NaN(), 64, 2) && ! d.prepare (48000.0, 0, 2),
        "prepare refuses a width outside 1..2, a rate outside 1 kHz … 768 kHz, and an empty block");
    ok (! d.isPrepared(), "...and a refused prepare leaves it unprepared");
    test::run (d.prepare (48000.0, 64, 1));
    a[0] = 0.5f;
    ok (! d.process (io, 2, 64) && a[0] == 0.5f, "a call wider than prepared is refused before anything moves");
    ok (! d.process (io, -1, 64) && ! d.process (io, 1, -1), "a negative width or length is refused");
    ok (d.process (io, 1, 0) && a[0] == 0.5f, "n == 0 is accepted and changes nothing");

    ok (! d.prepare (48000.0, 64, 5), "a re-prepare at a width it cannot honour is refused");
    a[0] = 0.25f;
    ok (! d.isPrepared() && ! d.process (io, 1, 64) && a[0] == 0.25f,
        "...and DISARMS: the old preparation no longer stands, process refuses and touches nothing");
    ok (d.prepare (1000.0, 64, 2) && d.prepare (768000.0, 64, 2), "the rate's edges, 1 kHz and 768 kHz, are accepted");

    group ("mix 0: the dry is bit-exact, whatever the loop holds");
    delay::TapeDelay m;
    test::run (m.prepare (44100.0, 512, 2));
    delay::TapeDelayParams p; p.mix = 0.0f; p.repeats = 0.9f; p.timeMs = 3.0f; p.offsetMs = 7.0f;
    m.setParams (p);
    std::vector<float> l (4096), r (4096);
    dc::Lcg g { 7u };
    for (auto& v : l) v = g.next();
    for (auto& v : r) v = g.next();
    auto l0 = l, r0 = r;
    float* lr[2] = { l.data(), r.data() };
    test::run (m.process (lr, 2, 4096));
    ok (l == l0 && r == r0, "the dry passes untouched at mix 0");
}

double tailFormula (double timeMs, double repeats, double offsetMs)
{
    const double t = timeMs * 0.001;
    return std::abs (offsetMs) * 0.001 + (repeats <= 0.0 ? t : repeats >= 1.0 ? 1.0e6 : t * (1.0 + std::log (0.001) / std::log (repeats)));
}

void tail()
{
    group ("tailSeconds: right from construction, and the discrete-echo estimate");
    delay::TapeDelay fresh;
    test::approx (fresh.tailSeconds(), tailFormula (350.0, (double) 0.35f, 0.0), 1e-5,
                  "a fresh instance at its default knobs already reports the formula (2.653 s), not a placeholder");
    fresh.setParams (delay::TapeDelayParams {});
    test::run (fresh.prepare (48000.0, 64, 2));
    test::approx (fresh.tailSeconds(), tailFormula (350.0, (double) 0.35f, 0.0), 1e-5, "...and still after the default knobs and a prepare");

    delay::TapeDelay d;
    delay::TapeDelayParams p; p.timeMs = 500.0f; p.repeats = 0.0f; p.offsetMs = 10.0f;
    d.setParams (p);
    test::approx (d.tailSeconds(), 0.51, 1e-6, "no feedback: one echo plus the offset");
    p.repeats = 0.5f; d.setParams (p);
    test::approx (d.tailSeconds(), tailFormula (500.0, 0.5, 10.0), 1e-5, "0.5: the first echo under a thousandth");
    p.repeats = 1.0f; d.setParams (p);
    ok (d.tailSeconds() > 1.0e5f, "full feedback never decays — the caller caps it");
}

void reslicing()
{
    group ("law 8a: per-sample state, so any slicing of one stream is bit-identical");
    dc::Knobs k; k.timeMs = 23.0f; k.repeats = 0.97f; k.darkPct = 60.0f; k.offsetMs = -11.0f; k.mix = 0.8f;
    std::vector<std::vector<float>> ref;
    bool same = true;
    for (int block : { 4096, 1, 64, 333, 517 })
    {
        dc::Case c { "slice", 48000.0, block, 2, dc::kNoiseBurst, 0.02, 0.25, k, {} };
        const auto out = renderCore (c);
        if (ref.empty()) ref = out; else same = same && out == ref;
    }
    ok (same, "blocks 1, 64, 333, 517 and 4096 give the same bits");
}

void resetIsRestart()
{
    group ("reset() is a stream restart: what follows does not depend on what came before");
    delay::TapeDelayParams p; p.timeMs = 3.1f; p.repeats = 0.7f; p.offsetMs = 4.3f; p.mix = 1.0f;
    std::vector<std::vector<float>> first;
    bool same = true;
    for (int warm : { 0, 1, 777, 48000, 262147 })   // 262147: the write clock past 2^18, where a float's ulp is 1/32
    {
        delay::TapeDelay d;
        test::run (d.prepare (48000.0, 4096, 2));
        d.setParams (p);
        std::vector<float> l (4096), r (4096);
        float* lr[2] = { l.data(), r.data() };
        dc::Lcg g { 21u };
        for (int done = 0; done < warm; done += 4096)
        {
            for (auto& v : l) v = 0.3f * g.next();
            for (auto& v : r) v = 0.3f * g.next();
            test::run (d.process (lr, 2, std::min (4096, warm - done)));
        }
        d.reset();
        dc::Lcg h { 99u };
        for (auto& v : l) v = 0.3f * h.next();
        for (auto& v : r) v = 0.3f * h.next();
        test::run (d.process (lr, 2, 4096));
        std::vector<std::vector<float>> out { l, r };
        if (first.empty()) first = out; else same = same && out == first;
    }
    ok (same, "after warm-ups of 0, 1, 777, 48000 and 262147 samples, the post-reset output is bit-identical");
}

void extremes()
{
    group ("absurd but finite knobs stay bounded: no NaN, no index off the line");
    bool finite = true;
    for (double sr : { 1000.0, 48000.0, 768000.0 })
    {
        delay::TapeDelay d;
        test::run (d.prepare (sr, 256, 2));
        std::vector<float> l (256), r (256);
        float* lr[2] = { l.data(), r.data() };
        dc::Lcg g { 5u };
        const float F = std::numeric_limits<float>::max();
        const float knobs[][5] = { { F, F, F, F, F }, { -F, -F, -F, -F, -F }, { F, 0.9f, 1e-30f, -F, 1.0f },
                                   { 350.0f, 1.0f, F, F, 1.0f }, { 1e-30f, 1.0f, 1.0f, 1e-30f, 1.0f },
                                   { 350.0f, 0.5f, 4000.0f, 0.0f, 0.5f } };
        for (const auto& k : knobs)
        {
            d.setParams ({ k[0], k[1], k[2], k[3], k[4] });
            for (int b = 0; b < 40; ++b)
            {
                for (auto& v : l) v = 0.5f * g.next();
                for (auto& v : r) v = 0.5f * g.next();
                test::run (d.process (lr, 2, 256));
                for (int i = 0; i < 256; ++i) finite = finite && std::isfinite (l[(std::size_t) i]) && std::isfinite (r[(std::size_t) i]);
                finite = finite && std::isfinite (d.shownTimeMs()) && std::isfinite (d.shownOffsetMs());
            }
        }
    }
    ok (finite, "±FLT_MAX and 1e-30 on every knob, at 1, 48 and 768 kHz, then back to sane: every sample finite");
}

void stoppedChannelAndPause()
{
    group ("law 11a/11c: a stopped channel drops its echoes; a pause glides the motor on");
    delay::TapeDelay d;
    test::run (d.prepare (48000.0, 256, 2));
    delay::TapeDelayParams p; p.timeMs = 50.0f; p.repeats = 0.9f; p.mix = 1.0f;
    d.setParams (p);

    std::vector<float> l (256), r (256);
    float* lr[2] = { l.data(), r.data() };
    dc::Lcg g { 3u };
    for (int b = 0; b < 8; ++b)
    {
        for (int i = 0; i < 256; ++i) { l[(std::size_t) i] = 0.3f * g.next(); r[(std::size_t) i] = 0.3f * g.next(); }
        test::run (d.process (lr, 2, 256));
    }
    for (int b = 0; b < 4; ++b)          // mono for a while: the right channel has STOPPED
    {
        std::fill (l.begin(), l.end(), 0.0f);
        test::run (d.process (lr, 1, 256));
    }
    float worstR = 0.0f;
    for (int b = 0; b < 40; ++b)         // back in stereo, digital silence on the right
    {
        std::fill (l.begin(), l.end(), 0.0f); std::fill (r.begin(), r.end(), 0.0f);
        test::run (d.process (lr, 2, 256));
        for (float v : r) worstR = std::max (worstR, std::abs (v));
    }
    ok (worstR == 0.0f, "the right channel comes back to silence, not to the echoes it held when it stopped");

    delay::TapeDelay x, y;
    test::run (x.prepare (48000.0, 64, 2)); test::run (y.prepare (48000.0, 64, 2));
    x.setParams (p); y.setParams (p);
    std::vector<float> z (64, 0.0f), z2 (64, 0.0f);
    float* zz[2] = { z.data(), z2.data() };
    test::run (x.process (zz, 2, 64)); test::run (y.process (zz, 2, 64));
    p.timeMs = 400.0f; x.setParams (p); y.setParams (p);
    test::run (x.process (nullptr, 0, 3000));
    for (int done = 0; done < 3000; done += 64)
    {
        std::fill (z.begin(), z.end(), 0.0f); std::fill (z2.begin(), z2.end(), 0.0f);
        test::run (y.process (zz, 2, std::min (64, 3000 - done)));
    }
    ok (x.shownTimeMs() == y.shownTimeMs() && x.shownTimeMs() > 50.0f && x.shownTimeMs() < 400.0f,
        "after a 3000-sample pause the head stands where 3000 samples of silence left it, mid-glide");
}

void glidesLand()
{
    group ("the heads LAND on their targets (LAW8-AUDIT.md item 6)");
    delay::TapeDelay d;
    test::run (d.prepare (48000.0, 480, 2));
    delay::TapeDelayParams p; p.timeMs = 350.0f; p.offsetMs = 20.0f;
    d.setParams (p);
    std::vector<float> l (480, 0.0f), r (480, 0.0f);
    float* lr[2] = { l.data(), r.data() };
    test::run (d.process (lr, 2, 480));
    p.timeMs = 6000.0f; p.offsetMs = 7.0f;
    d.setParams (p);
    for (int b = 0; b < 1000; ++b)   // ten seconds: sixty-six of the motor's time constants
        test::run (d.process (lr, 2, 480));
    // The targets, spelt as the stage spells them: the motor's in double from TIME, the offset head's in float.
    const double motorTarget  = (double) (float) ((double) 6000.0f * 0.001 * 48000.0);
    const double offsetTarget = (double) (7.0f * 0.001f * (float) 48000.0);
    ok (d.headSamples() == motorTarget, "the motor stands ON its target, bit for bit (a float motor parks 112 samples "
                                        "short, a double one without the landing ~2e-7 short): "
                                        + std::to_string (d.headSamples() - motorTarget));
    ok (d.offsetHeadSamples (1) == offsetTarget && d.offsetHeadSamples (0) == 0.0,
        "the offset heads stand ON theirs: the right on 7 ms of samples, the left on zero");
}

std::int64_t ulpsApart (float a, float b)
{
    std::int32_t ia, ib;
    std::memcpy (&ia, &a, 4); std::memcpy (&ib, &b, 4);
    if (ia < 0) ia = std::numeric_limits<std::int32_t>::min() - ia;   // a monotone integer line through ±0
    if (ib < 0) ib = std::numeric_limits<std::int32_t>::min() - ib;
    return std::abs ((std::int64_t) ia - (std::int64_t) ib);
}

void detTanhAccuracy()
{
    group ("the loop's det tanh: within 1 ulp of the correctly rounded tanh, both branches and the seam");
    std::int64_t worst = 0;
    int points = 0;
    auto probe = [&] (float x)
    {
        const float want = (float) std::tanh ((double) x);   // double tanh, rounded once: the reference
        worst = std::max (worst, ulpsApart (delay::detail::detTanh (x), want));
        ++points;
    };
    for (float sign : { 1.0f, -1.0f })
    {
        for (double x = 1.0e-30; x < 30.0; x *= 1.0001)       // every decade from the flush floor to saturation
            probe (sign * (float) x);
        float s = 0.0625f;                                      // the seam between the series and the exponential:
        for (int k = 0; k < 4096; ++k) s = std::nextafter (s, 0.0f);   // 4096 floats either side of it, one by one
        for (int k = 0; k < 8193; ++k, s = std::nextafter (s, 1.0f))
            probe (sign * s);
    }
    probe (0.0f);
    ok (worst <= 1, std::to_string (points) + " points, worst " + std::to_string (worst) + " ulp");
}

void poisonAndLaw8()
{
    group ("poison: a bad input sample stays the dry's; the loop never records it");
    delay::TapeDelay d;
    test::run (d.prepare (48000.0, 128, 1));
    delay::TapeDelayParams p; p.timeMs = 1.0f; p.repeats = 1.0f; p.mix = 1.0f;
    d.setParams (p);
    std::vector<float> x (128, 0.1f);
    x[10] = std::numeric_limits<float>::quiet_NaN();
    x[20] = std::numeric_limits<float>::infinity();
    float* io[1] = { x.data() };
    test::run (d.process (io, 1, 128));
    bool finiteElsewhere = true;
    for (int i = 0; i < 128; ++i)
        if (i != 10 && i != 20) finiteElsewhere = finiteElsewhere && std::isfinite (x[(std::size_t) i]);
    ok (std::isnan (x[10]) && std::isinf (x[20]) && finiteElsewhere, "the NaN and the inf pass as dry; every other sample is finite");
    std::fill (x.begin(), x.end(), 0.0f);
    test::run (d.process (io, 1, 128));
    bool clean = true;
    for (float v : x) clean = clean && std::isfinite (v);
    ok (clean, "the next block is finite: nothing poisoned the line");

    group ("law 8: a minute of silence after a bloom leaves no subnormal anywhere, and the silence holds");
    delay::TapeDelay s;
    test::run (s.prepare (48000.0, 480, 2));
    delay::TapeDelayParams q; q.timeMs = 37.0f; q.repeats = 0.6f; q.offsetMs = 20.0f; q.mix = 1.0f;
    s.setParams (q);
    std::vector<float> l (480), r (480);
    float* lr[2] = { l.data(), r.data() };
    dc::Lcg g { 11u };
    for (int i = 0; i < 480; ++i) { l[(std::size_t) i] = 0.5f * g.next(); r[(std::size_t) i] = 0.5f * g.next(); }
    test::run (s.process (lr, 2, 480));
    q.offsetMs = 0.0f; s.setParams (q);   // the offset glide now heads for ZERO — the case that parks subnormals
    bool anySub = false;
    int lastLoud = -1;
    for (int b = 0; b < 6000; ++b)
    {
        std::fill (l.begin(), l.end(), 0.0f); std::fill (r.begin(), r.end(), 0.0f);
        test::run (s.process (lr, 2, 480));
        anySub = anySub || s.anySubnormalState();
        const bool silent = std::all_of (l.begin(), l.end(), [] (float v) { return v == 0.0f; })
                         && std::all_of (r.begin(), r.end(), [] (float v) { return v == 0.0f; });
        if (! silent) lastLoud = b;
    }
    ok (! anySub, "no subnormal in the loop filter, the offset glide or the line, at any block");
    ok (lastLoud >= 0 && lastLoud < 5000, "the tail ends in exact digital silence (last sound at block "
                                          + std::to_string (lastLoud) + ") and STAYS silent to the minute's end");
}

void noAllocation()
{
    group ("RT: process, reset and the setters never allocate");
    delay::TapeDelay d;
    test::run (d.prepare (96000.0, 256, 2));
    std::vector<float> l (1024, 0.25f), r (1024, -0.25f);
    float* lr[2] = { l.data(), r.data() };
    const long long before = alloc::count.load();
    delay::TapeDelayParams p; p.timeMs = 700.0f; p.offsetMs = -20.0f;
    d.setParams (p);
    test::run (d.process (lr, 2, 1024));
    test::run (d.process (lr, 1, 1024));
    test::run (d.process (nullptr, 0, 1024));
    d.reset(); d.reset();
    test::run (d.process (lr, 2, 1024, false));
    felitronics::test::okNoAlloc (alloc::count.load() == before, "no allocation: setParams, process (stereo, mono, pause, stand-by), reset");
}
} // namespace

int main()
{
    std::printf ("felitronics::delay::TapeDelay\n");
    goldens();
    parameterLaw();
    contract();
    tail();
    reslicing();
    resetIsRestart();
    extremes();
    stoppedChannelAndPause();
    glidesLand();
    detTanhAccuracy();
    poisonAndLaw8();
    noAllocation();
    return felitronics::test::report();
}

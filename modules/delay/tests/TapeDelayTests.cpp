// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// delay::TapeDelay — the goldens the app's ORIGINAL stage rendered (the parity harness that wrote them runs the
// app's class itself; here core CI holds the line without it), the parameter law, and the stage contract:
// re-slicing, a stopped channel, a pause, poison, law 8, no allocation.
#include <felitronics_test.h>
#include <alloc_counter.h>   // installs the allocation counter

#include "delay_cases.h"
#include "delay_golden.h"

#include <cmath>
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

// Tolerances that hold on every row: the goldens are bit-exact on the row that wrote them, but std::tanh and
// std::exp are the system libm (by design — see the header) and FMA contraction differs between ISAs.
constexpr double kSampleTol = 1.0e-6;
constexpr double kSumRelTol = 1.0e-6;

void goldens()
{
    group ("goldens: the app's original DelayStage, rendered by the parity harness");
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
        ok (worst <= kSampleTol, std::string (c.name) + ": 3 windows within 1e-6 of the app (worst " + std::to_string (worst) + ")");

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

void parameterLaw()
{
    group ("the parameter law: sync ladder and DARK");
    ok (delay::divisionMs (delay::kDivisionDefault, 120.0) == 500.0f, "1/4 at 120 BPM is 500 ms");
    ok (delay::divisionMs (0, 40.0) == 6000.0f, "1/1 at 40 BPM is 6000 ms — the line's design maximum");
    ok (delay::divisionMs (99, 120.0) == delay::divisionMs (12, 120.0) && delay::divisionMs (-3, 120.0) == delay::divisionMs (0, 120.0),
        "an index off the ladder is clamped onto it");
    ok (delay::divisionMs (0, 40.0) <= (float) (delay::TapeDelay::kMaxDelaySeconds * 1000.0), "the slowest division fits the line");
    ok (delay::darkHzForPercent (0.0f) == 16000.0f, "DARK 0 % is the bright ceiling, 16 kHz");
    test::approx (delay::darkHzForPercent (100.0f), 800.0, 0.01, "DARK 100 % is the dark floor, 800 Hz");
    ok (delay::darkHzForPercent (45.0f) > 3500.0f && delay::darkHzForPercent (45.0f) < 4500.0f, "the default 45 % is ~4 kHz");
}

void contract()
{
    group ("contract: prepare is binding, refusals move nothing, n == 0 is a no-op");
    delay::TapeDelay d;
    float a[64] {}, b[64] {};
    float* io[2] = { a, b };
    ok (! d.process (io, 1, 64), "an unprepared stage refuses");
    ok (! d.prepare (48000.0, 64, 3) && ! d.prepare (48000.0, 64, 0) && ! d.prepare (0.0, 64, 2)
        && ! d.prepare (std::numeric_limits<double>::quiet_NaN(), 64, 2) && ! d.prepare (48000.0, 0, 2),
        "prepare refuses a width outside 1..2, a bad rate and an empty block");
    ok (! d.isPrepared(), "...and a refused prepare leaves it unprepared");
    test::run (d.prepare (48000.0, 64, 1));
    a[0] = 0.5f;
    ok (! d.process (io, 2, 64) && a[0] == 0.5f, "a call wider than prepared is refused before anything moves");
    ok (! d.process (io, -1, 64) && ! d.process (io, 1, -1), "a negative width or length is refused");
    ok (d.process (io, 1, 0) && a[0] == 0.5f, "n == 0 is accepted and changes nothing");

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

    // A pause moves the shared motor exactly as silence would.
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

    group ("law 8: a minute of silence after a full-feedback bloom leaves no subnormal anywhere");
    delay::TapeDelay s;
    test::run (s.prepare (48000.0, 480, 2));
    delay::TapeDelayParams q2; q2.timeMs = 37.0f; q2.repeats = 0.6f; q2.offsetMs = 20.0f; q2.mix = 1.0f;
    s.setParams (q2);
    std::vector<float> l (480), r (480);
    float* lr[2] = { l.data(), r.data() };
    dc::Lcg g { 11u };
    for (int i = 0; i < 480; ++i) { l[(std::size_t) i] = 0.5f * g.next(); r[(std::size_t) i] = 0.5f * g.next(); }
    test::run (s.process (lr, 2, 480));
    q2.offsetMs = 0.0f; s.setParams (q2);   // the offset glide now heads for ZERO — the case that parks subnormals
    bool anySub = false, silentAtEnd = false;
    for (int b = 0; b < 6000; ++b)
    {
        std::fill (l.begin(), l.end(), 0.0f); std::fill (r.begin(), r.end(), 0.0f);
        test::run (s.process (lr, 2, 480));
        anySub = anySub || s.anySubnormalState();
        if (b == 5999) silentAtEnd = std::all_of (l.begin(), l.end(), [] (float v) { return v == 0.0f; });
    }
    ok (! anySub, "no subnormal in the loop filter, the offset glide or the line, at any block");
    ok (silentAtEnd, "...and the tail ends in exact digital silence");
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

void tail()
{
    group ("tailSeconds: the discrete-echo estimate");
    delay::TapeDelay d;
    delay::TapeDelayParams p; p.timeMs = 500.0f; p.repeats = 0.0f; p.offsetMs = 10.0f;
    d.setParams (p);
    test::approx (d.tailSeconds(), 0.51, 1e-6, "no feedback: one echo plus the offset");
    p.repeats = 0.5f; d.setParams (p);
    test::approx (d.tailSeconds(), 0.01 + 0.5 * (1.0 + std::log (0.001) / std::log (0.5)), 1e-5, "0.5: the first echo under a thousandth");
    p.repeats = 1.0f; d.setParams (p);
    ok (d.tailSeconds() > 1.0e5f, "full feedback never decays — the caller caps it");
}
} // namespace

int main()
{
    std::printf ("felitronics::delay::TapeDelay\n");
    goldens();
    parameterLaw();
    contract();
    reslicing();
    stoppedChannelAndPause();
    poisonAndLaw8();
    noAllocation();
    tail();
    return felitronics::test::report();
}

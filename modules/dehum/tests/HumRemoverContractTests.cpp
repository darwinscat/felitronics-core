// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.
//
// dehum::HumRemover as a STAGE: re-slicing bit identity (law 8a), a pause is silence (11c), a stopped
// channel keeps its model (11a), prepare refuses and disarms (11b), no allocation, digital black in is
// bit-exact black out, a minute of silence leaves no state in the subnormal band (law 8), NaN/Inf/huge
// input recovers, and a mode change ramps the correction out over 20 ms before re-arming. The law-11
// property sweep runs in the census (modules/laws/tests/CallContractTests.cpp).
#include <felitronics_test.h>
#include <alloc_counter.h>   // installs the allocation counter: EVERY form of `new`, over-aligned included
#include "dehum_fixtures.h"

#include <cstdio>
#include <cstring>
#include <limits>
#include <string>

using namespace dehumtest;
namespace dh = felitronics::dehum;
using felitronics::test::ok;
using felitronics::test::group;
using HR = dh::HumRemover<>;

static constexpr double kFs = 48000.0;

// Stereo programme: hum (different phases per channel) + a tone + low noise.
static void programme (std::vector<float>& l, std::vector<float>& r, std::int64_t start, int n, std::uint32_t seed)
{
    Hum h0; h0.seed (1u); Hum h1; h1.seed (2u);
    Lcg rng { seed };
    l.resize ((std::size_t) n); r.resize ((std::size_t) n);
    for (int i = 0; i < n; ++i)
    {
        const std::int64_t t = start + i;
        l[(std::size_t) i] = (float) (h0.at (t) + 0.1 * sinAt (329.6276, t, kFs, 0.0)) + 0.001f * rng.bipolar();
        r[(std::size_t) i] = (float) (h1.at (t) + 0.1 * sinAt (246.9417, t, kFs, 1.0)) + 0.001f * rng.bipolar();
    }
}

static bool run2 (HR& h, std::vector<float>& l, std::vector<float>& r, int from, int n)
{
    float* io[2] { l.data() + from, r.data() + from };
    return h.process (io, 2, n);
}

// The whole of two planes, in 512-sample calls (the last one short).
static void runAll2 (HR& h, std::vector<float>& l, std::vector<float>& r)
{
    for (std::size_t o = 0; o < l.size(); o += 512)
        felitronics::test::run (run2 (h, l, r, (int) o, (int) std::min<std::size_t> (512, l.size() - o)));
}

static bool sameWeights (const HR& a, const HR& b)
{
    for (int c = 0; c < 2; ++c)
        for (int k = 1; k <= HR::kMaxHarmonics; ++k)
        {
            const auto x = a.weight (c, k), y = b.weight (c, k);
            if (std::memcmp (&x.re, &y.re, sizeof (double)) != 0 || std::memcmp (&x.im, &y.im, sizeof (double)) != 0) return false;
        }
    const double fa = a.mainsHz(), fb = b.mainsHz();   // and the tracker's estimate
    return std::memcmp (&fa, &fb, sizeof (double)) == 0 && a.quiet() == b.quiet() && a.detected() == b.detected();
}

static bool bitsEqual (const std::vector<float>& a, const std::vector<float>& b)
{
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i) if (! core::sameBits (a[i], b[i])) return false;
    return true;
}

// Null (dB) of the hum's harmonics k = 1..24 in `y` relative to `x`, over [from, from + len).
static double worstNull (const std::vector<float>& y, const std::vector<float>& x, std::int64_t from, std::int64_t len, double f0 = 50.0)
{
    std::vector<double> fr;
    for (int k = 1; k <= 24; ++k) fr.push_back (f0 * k);
    const auto ay = fitAmps (y, from, len, fr, kFs), ax = fitAmps (x, from, len, fr, kFs);
    double w = -400.0;
    for (std::size_t k = 0; k < fr.size(); ++k) w = std::max (w, dbOf (ay[k] / ax[k]));
    return w;
}

//------------------------------------------------------------------------------------------------------
// LAW 8a: one stream, three slicings — one call per event stretch, one sample per call, random lengths —
// bit-identical output and weights. The stream carries a mode change (50 -> 60 at a fixed sample, i.e. the
// caller's own event timeline), a run of exact zeros (digital black) and a NaN.
static void testReslicing()
{
    group ("law 8a: bit-identical under re-slicing (whole / 1-sample / random), through a mode change, black and NaN");
    const int N = 96000, SWITCH = 48017;
    std::vector<float> l0, r0;
    programme (l0, r0, 0, N, 9u);
    for (int i = 30000; i < 30200; ++i) { l0[(std::size_t) i] = 0.0f; r0[(std::size_t) i] = 0.0f; }
    l0[40000] = std::numeric_limits<float>::quiet_NaN();

    auto render = [&] (int mode, std::vector<float>& l, std::vector<float>& r, HR& h)
    {
        l = l0; r = r0;
        h.setParams (params (dh::Mains::Hz50, 0.5, 3.0));   // params first: prepare applies them at once
        felitronics::test::run (h.prepare (kFs, 512, 2));
        Lcg rng { 4242u };
        auto span = [&] (int from, int to)
        {
            for (int o = from; o < to;)
            {
                const int len = mode == 0 ? to - o : (mode == 1 ? 1 : 1 + (int) (rng.uniform() * 700.0));
                const int m = std::min (len, to - o);
                felitronics::test::run (run2 (h, l, r, o, m));
                o += m;
            }
        };
        span (0, SWITCH);
        h.setParams (params (dh::Mains::Hz60, 0.5, 3.0));
        span (SWITCH, N);
    };
    std::vector<float> la, ra, lb, rb, lc, rc;
    HR a, b, c;
    render (0, la, ra, a);
    render (1, lb, rb, b);
    render (2, lc, rc, c);
    ok (bitsEqual (la, lb) && bitsEqual (ra, rb), "one sample per call == one call: output bit-identical");
    ok (bitsEqual (la, lc) && bitsEqual (ra, rc), "random slices == one call: output bit-identical");
    ok (sameWeights (a, b) && sameWeights (a, c), "and the weights are bit-identical");
    ok (a.running() && a.mainsHz() > 59.9, "precondition: the stream did switch to 60 Hz");
    double moved = 0.0;
    for (std::size_t i = 0; i < la.size(); ++i) if (std::isfinite (la[i])) moved = std::max (moved, (double) std::fabs (la[i] - l0[i]));
    ok (moved > 0.01, "precondition: the stage changed the stream (max |out - in| " + std::to_string (moved) + ")");
}

//------------------------------------------------------------------------------------------------------
// LAW 11c: a clock-only gap equals the same number of zeros at a live width — output and weights, bit for
// bit — including a gap that lands inside a 20 ms ramp and the grid re-arm after it. The programme after
// the gap opens with exact zeros, so a zero-run counter the gap failed to advance would show.
static void testPauseIsSilence()
{
    group ("law 11c: a pause is silence (gaps 1..48000, with and without a ramp in flight)");
    const int gaps[] = { 1, 2, 7, 63, 64, 65, 127, 480, 1000, 4801, 48000 };
    int cases = 0, same = 0;
    for (bool rampInFlight : { false, true })
        for (int G : gaps)
        {
            HR a, b;
            felitronics::test::run (a.prepare (kFs, 512, 2));
            felitronics::test::run (b.prepare (kFs, 512, 2));
            a.setParams (params (dh::Mains::Hz50, 0.5, 3.0)); b.setParams (params (dh::Mains::Hz50, 0.5, 3.0));
            std::vector<float> l, r, l2, r2;
            programme (l, r, 0, 24000 + 13, 3u); l2 = l; r2 = r;
            felitronics::test::run (run2 (a, l, r, 0, (int) l.size()));
            felitronics::test::run (run2 (b, l2, r2, 0, (int) l2.size()));
            if (rampInFlight) { a.setParams (params (dh::Mains::Off, 3.0)); b.setParams (params (dh::Mains::Off, 3.0)); }

            float* none[2] { nullptr, nullptr };
            felitronics::test::run (a.process (none, 0, G));                                  // THE GAP
            std::vector<float> zl ((std::size_t) G, 0.0f), zr ((std::size_t) G, 0.0f);
            felitronics::test::run (run2 (b, zl, zr, 0, G));                                 // ...vs silence

            std::vector<float> pl, pr;
            programme (pl, pr, 24013 + G, 9600, 5u);
            for (int i = 0; i < 10; ++i) { pl[(std::size_t) i] = 0.0f; pr[(std::size_t) i] = 0.0f; }
            std::vector<float> ql = pl, qr = pr;
            felitronics::test::run (run2 (a, pl, pr, 0, 9600));
            felitronics::test::run (run2 (b, ql, qr, 0, 9600));
            ++cases;
            if (bitsEqual (pl, ql) && bitsEqual (pr, qr) && sameWeights (a, b) && a.running() == b.running()) ++same;
            else std::printf ("      MISMATCH: gap %d, ramp in flight %d\n", G, (int) rampInFlight);
        }
    ok (same == cases, "pause == silence, output + weights bit-identical in " + std::to_string (same) + " of " + std::to_string (cases) + " cases");
}

//------------------------------------------------------------------------------------------------------
// LAW 11a: a channel that stops keeps its weights (a model of the hum on that input) and cancels at once
// when it returns.
static void testStoppedChannelKeepsModel()
{
    group ("law 11a: a stopped channel's weights are kept and cancel at once on return");
    HR h;
    h.setParams (params (dh::Mains::Hz50, 3.0));   // params first: prepare applies them at once
    felitronics::test::run (h.prepare (kFs, 512, 2));
    Hum h1; h1.seed (2u);
    std::vector<float> l ((std::size_t) 96000), r ((std::size_t) 96000);
    for (int i = 0; i < 96000; ++i) { l[(std::size_t) i] = (float) h1.at (i); r[(std::size_t) i] = (float) h1.at (i); }
    runAll2 (h, l, r);
    HR::Weight before[25];
    for (int k = 1; k <= 24; ++k) before[k] = h.weight (1, k);

    std::vector<float> m ((std::size_t) 48000);
    for (int i = 0; i < 48000; ++i) m[(std::size_t) i] = (float) h1.at (96000 + i);
    for (int o = 0; o < 48000; o += 512) { float* io[1] { m.data() + o }; felitronics::test::run (h.process (io, 1, std::min (512, 48000 - o))); }
    bool kept = true;
    for (int k = 1; k <= 24; ++k)
    {
        const auto w = h.weight (1, k);
        kept = kept && std::memcmp (&w, &before[k], sizeof w) == 0;
    }
    ok (kept, "channel 1's weights are bit-identical after a second of mono calls");
    ok (std::fabs (before[1].re) + std::fabs (before[1].im) > 1e-3, "precondition: there was a model to keep");

    const int R = 4800;   // the first 100 ms back
    std::vector<float> rl ((std::size_t) R), rr ((std::size_t) R);
    for (int i = 0; i < R; ++i) { rl[(std::size_t) i] = (float) h1.at (144000 + i); rr[(std::size_t) i] = rl[(std::size_t) i]; }
    const std::vector<float> rin = rr;
    felitronics::test::run (run2 (h, rl, rr, 0, R));
    const double nul = worstNull (rr, rin, 0, R);
    ok (nul <= -50.0, "and on return it cancels from the first sample: worst harmonic " + std::to_string (nul) + " dB over the first 100 ms");
}

//------------------------------------------------------------------------------------------------------
static void testPrepareAndAlloc()
{
    group ("law 11b: prepare refuses and disarms; no allocation in prepare / setParams / process / reset");
    {
        HR h;
        ok (! h.prepare (7999.0, 512, 2), "fs below kMinSampleRate is refused");
        ok (! h.prepare (std::numeric_limits<double>::quiet_NaN(), 512, 2), "a NaN rate is refused");
        ok (! h.prepare (std::numeric_limits<double>::infinity(), 512, 2), "an infinite rate is refused");
        ok (! h.prepare (48000.0, 0, 2), "maxBlock 0 is refused");
        ok (h.prepare (8000.0, 512, 2), "8000 Hz is audio");
        h.setParams (params (dh::Mains::Hz50, 0.5));
        float z[4] {}; float* io[2] { z, z + 2 };
        ok (h.process (io, 2, 2), "precondition: prepared, it runs");
        ok (! h.prepare (48000.0, 512, 0), "zero channels are refused");
        ok (! h.process (io, 2, 2) && ! h.isPrepared(), "and the refused prepare DISARMED it (process refuses)");
    }
    {
        HR h;
        felitronics::test::run (h.prepare (8000.0, 64, 1));
        dh::HumRemoverParams p = params (dh::Mains::Hz50, 0.5);
        p.harmonics = 64; p.maxHz = 5000.0;
        h.setParams (p); felitronics::test::run (h.prepare (8000.0, 64, 1));
        ok (h.activeHarmonics() == 64, "8 kHz, 64 asked, maxHz 5000: 0.9·Nyquist allows 72, the cap is 64 ("
                                        + std::to_string (h.activeHarmonics()) + ")");
        p.maxHz = 2000.0; h.setParams (p); felitronics::test::run (h.prepare (8000.0, 64, 1));
        ok (h.activeHarmonics() == 40, "maxHz 2000 at 50 Hz: 40 harmonics (" + std::to_string (h.activeHarmonics()) + ")");
        p.harmonics = 24; p.mains = dh::Mains::Hz60; h.setParams (p); felitronics::test::run (h.prepare (8000.0, 64, 1));
        ok (h.activeHarmonics() == 24 && std::fabs (h.mainsHz() - 60.0) < 1e-5, "default 24 at 60 Hz");
    }
    std::vector<float> l, r;
    programme (l, r, 0, 9600, 1u);
    const long long before = alloc::count.load();
    {
        HR h;
        h.setParams (params (dh::Mains::Hz50, 0.5));   // params first: prepare applies them at once
        felitronics::test::run (h.prepare (kFs, 512, 2));
        felitronics::test::run (run2 (h, l, r, 0, 4800));
        h.setParams (params (dh::Mains::Hz60, 0.5));
        felitronics::test::run (run2 (h, l, r, 4800, 4800));
        float* none[2] { nullptr, nullptr };
        felitronics::test::run (h.process (none, 0, 1000));
        h.reset();
    }
    felitronics::test::okNoAlloc (alloc::count.load() == before, "no allocation: construct, prepare, setParams, process, gap, reset");
    std::printf ("      footprint: sizeof(HumRemover<64>) = %zu bytes, sizeof(HumRemover<24>) = %zu bytes\n",
                 sizeof (dh::HumRemover<64>), sizeof (dh::HumRemover<24>));
}

//------------------------------------------------------------------------------------------------------
// DIGITAL BLACK: 64 exact zeros mute the correction (bit-exact zeros out, -0.0 kept) and hold the weights;
// the hum is cancelled again from its first sample back.
static void testDigitalBlack()
{
    group ("digital black: bit-exact black out after 64 zeros, weights held, cancels at once on return");
    HR h;
    h.setParams (params (dh::Mains::Hz50, 0.5));   // params first: prepare applies them at once
    felitronics::test::run (h.prepare (kFs, 512, 2));
    // hum only: a tone in the programme would leave its own notch ringing (its beat lives in the nearest
    // harmonic's weight and decays with tau once the tone stops) and the return check would read that
    Hum h0; h0.seed (1u); Hum h1; h1.seed (2u);
    std::vector<float> l ((std::size_t) (6 * 48000)), r ((std::size_t) (6 * 48000));
    for (int i = 0; i < 6 * 48000; ++i) { l[(std::size_t) i] = (float) h0.at (i); r[(std::size_t) i] = (float) h1.at (i); }
    runAll2 (h, l, r);
    HR snap = h;

    const int Z = 4800 + 7;
    std::vector<float> zl ((std::size_t) Z, 0.0f), zr ((std::size_t) Z, -0.0f);
    const std::vector<float> zl0 = zl, zr0 = zr;
    felitronics::test::run (run2 (h, zl, zr, 0, 100));
    felitronics::test::run (run2 (h, zl, zr, 100, Z - 100));
    bool black = true, early = false;
    for (int i = 0; i < Z; ++i)
    {
        if (i >= HR::kBlackRun - 1) black = black && core::sameBits (zl[(std::size_t) i], zl0[(std::size_t) i]) && core::sameBits (zr[(std::size_t) i], zr0[(std::size_t) i]);
        else early = early || zl[(std::size_t) i] != 0.0f;
    }
    ok (black, "from the 64th zero on, the output is the input's bits (+0.0 and -0.0)");
    ok (early, "(the first 63 zeros still carry the correction, as documented)");
    ok (sameWeights (h, snap), "the weights held through the black stretch");

    // back: the same hum, in phase with the stream
    std::vector<float> bl ((std::size_t) 4800), br ((std::size_t) 4800);
    for (int i = 0; i < 4800; ++i) { bl[(std::size_t) i] = (float) h0.at (6 * 48000 + Z + i); br[(std::size_t) i] = (float) h1.at (6 * 48000 + Z + i); }
    const std::vector<float> bin = bl;
    felitronics::test::run (run2 (h, bl, br, 0, 4800));
    const double nul = worstNull (bl, bin, 0, 4800);
    ok (nul <= -50.0, "the hum is cancelled from its first sample back: worst harmonic " + std::to_string (nul) + " dB over 100 ms");
}

//------------------------------------------------------------------------------------------------------
// LAW 8: a minute of silence. Digital black holds normal weights and emits exact zeros. Quiet that is not
// black (±1e-20 noise) lets the weights decay; the grid flush keeps every weight at 0 or >= 1e-15 — never
// the band below, where an unflushed decay would sit and eventually reach the subnormals.
static void testMinuteOfSilence()
{
    group ("law 8: 60 s of black and 60 s of quiet leave no state in the subnormal band");
    for (bool black : { true, false })
    {
        HR h;
        h.setParams (params (dh::Mains::Hz50, 0.5));   // params first: prepare applies them at once
        felitronics::test::run (h.prepare (kFs, 512, 2));
        std::vector<float> l, r;
        programme (l, r, 0, 2 * 48000, 7u);
        runAll2 (h, l, r);
        Lcg rng { 99u };
        bool band = true, subOut = false, zeroOut = true;
        std::vector<float> a (512), b (512);
        for (int blk = 0; blk < 60 * 48000 / 512; ++blk)
        {
            for (int i = 0; i < 512; ++i) { a[(std::size_t) i] = black ? 0.0f : 1e-20f * rng.bipolar(); b[(std::size_t) i] = black ? 0.0f : 1e-20f * rng.bipolar(); }
            felitronics::test::run (run2 (h, a, b, 0, 512));
            for (int i = 0; i < 512; ++i)
            {
                subOut = subOut || std::fpclassify (a[(std::size_t) i]) == FP_SUBNORMAL || std::fpclassify (b[(std::size_t) i]) == FP_SUBNORMAL;
                if (blk > 0) zeroOut = zeroOut && a[(std::size_t) i] == 0.0f && b[(std::size_t) i] == 0.0f;
            }
            for (int c = 0; c < 2; ++c)
                for (int k = 1; k <= 24; ++k)
                {
                    const auto w = h.weight (c, k);
                    for (double v : { w.re, w.im })
                        band = band && (v == 0.0 || std::fabs (v) >= 1e-15) && std::fpclassify (v) != FP_SUBNORMAL;
                }
        }
        const std::string what = black ? "digital black" : "quiet (1e-20 noise)";
        ok (band, what + ": after every 64-aligned call each weight is 0 or >= 1e-15 (none below, none subnormal)");
        ok (! subOut, what + ": no output sample is subnormal");
        if (black) ok (zeroOut, what + ": the output is exact zero");
    }
}

//------------------------------------------------------------------------------------------------------
// POISON: NaN / +Inf / -Inf samples pass through and adapt nothing — the model survives, so the null holds
// right after; a burst near float max cannot leave a model beyond the ceiling, and the null comes back.
static void testPoison()
{
    group ("NaN / Inf / huge input: passes through, the model survives or recovers");
    HR h;
    h.setParams (params (dh::Mains::Hz50, 3.0));   // params first: prepare applies them at once
    felitronics::test::run (h.prepare (kFs, 512, 1));
    Hum hum; hum.seed (31u);
    const int N = 7 * 48000;
    std::vector<float> x ((std::size_t) N);
    for (int i = 0; i < N; ++i) x[(std::size_t) i] = (float) hum.at (i);
    const int P0 = 2 * 48000, B0 = 4 * 48000;
    x[(std::size_t) P0]       = std::numeric_limits<float>::quiet_NaN();
    x[(std::size_t) P0 + 50]  = std::numeric_limits<float>::infinity();
    x[(std::size_t) P0 + 51]  = -std::numeric_limits<float>::infinity();
    for (int i = 0; i < 100; ++i) x[(std::size_t) (P0 + 200 + i)] = std::numeric_limits<float>::quiet_NaN();
    for (int i = 0; i < 256; ++i) x[(std::size_t) (B0 + i)] = (i & 1) ? -3.0e38f : 3.0e38f;
    const std::vector<float> in = x;
    runMono (h, x);

    bool finiteElsewhere = true;
    for (int i = P0; i < P0 + 48000; ++i)
    {
        const bool bad = ! std::isfinite (in[(std::size_t) i]);
        finiteElsewhere = finiteElsewhere && (bad || std::isfinite (x[(std::size_t) i]));
    }
    ok (finiteElsewhere, "only the injected samples are non-finite in the output");
    const double nulAfter = worstNull (x, in, P0 + 4800, 48000);
    ok (nulAfter <= -50.0, "the model survived: worst harmonic " + std::to_string (nulAfter) + " dB in the second after the NaN/Inf");

    bool finiteAfterBurst = true;
    for (int i = B0 + 256; i < N; ++i) finiteAfterBurst = finiteAfterBurst && std::isfinite (x[(std::size_t) i]);
    ok (finiteAfterBurst, "after a +-3e38 burst the output is finite");
    bool bounded = true;
    for (int k = 1; k <= 24; ++k) { const auto w = h.weight (0, k); bounded = bounded && std::fabs (w.re) <= 1.0 && std::fabs (w.im) <= 1.0; }
    ok (bounded, "and the model is a hum model again (|w| <= 1)");
    const double nulBurst = worstNull (x, in, B0 + 2 * 48000, 48000);
    ok (nulBurst <= -50.0, "the null comes back within 2 s of the burst: worst harmonic " + std::to_string (nulBurst) + " dB");
}

//------------------------------------------------------------------------------------------------------
// MODES. 50 -> Off at an arbitrary sample: the correction (in - out) is the unswitched twin's correction
// times a 20 ms linear ramp, then dry bit-exact. 50 -> 60: ramp, re-arm with zero weights, and the 60 Hz
// null forms. Auto abstains (dry) on an input with no hum. Off -> 50 starts at the next grid boundary.
static void testModes()
{
    group ("modes: 20 ms ramp of the correction, re-arm, Auto abstains on no hum, Off -> 50");
    std::vector<float> l, r;
    programme (l, r, 0, 4 * 48000, 21u);
    HR a;
    a.setParams (params (dh::Mains::Hz50, 3.0));   // params first: prepare applies them at once
    felitronics::test::run (a.prepare (kFs, 512, 2));
    const int SW = 2 * 48000 + 17;
    felitronics::test::run (run2 (a, l, r, 0, SW));
    HR twin = a;
    std::vector<float> l2 (l.begin() + SW, l.end()), r2 (r.begin() + SW, r.end());
    std::vector<float> lt = l2, rt = r2;
    const std::vector<float> in = l2;
    a.setParams (params (dh::Mains::Off, 3.0));
    felitronics::test::run (run2 (a, l2, r2, 0, (int) l2.size()));
    felitronics::test::run (run2 (twin, lt, rt, 0, (int) lt.size()));

    core::LinearSmoother g { 1.0f };
    g.reset (kFs, HR::kRampSeconds);
    g.setTargetValue (0.0f);
    const int ramp = (int) std::floor (HR::kRampSeconds * kFs);
    double worst = 0.0, corrPeak = 0.0;
    for (int i = 0; i < ramp; ++i)
    {
        const double gi = (double) g.getNextValue();
        const double c = (double) in[(std::size_t) i] - (double) l2[(std::size_t) i];
        const double ct = (double) in[(std::size_t) i] - (double) lt[(std::size_t) i];
        worst = std::max (worst, std::fabs (c - gi * ct));
        corrPeak = std::max (corrPeak, std::fabs (ct));
    }
    ok (corrPeak > 1e-3, "precondition: there was a correction to ramp (" + std::to_string (corrPeak) + ")");
    ok (worst <= 1e-6, "50 -> Off: the correction is the twin's times the 20 ms linear ramp (max dev " + std::to_string (worst) + ")");
    bool dry = true;
    for (int i = ramp + 64; i < (int) in.size(); ++i) dry = dry && core::sameBits (l2[(std::size_t) i], in[(std::size_t) i]);
    ok (dry && ! a.running(), "after the ramp and the next grid boundary: dry, bit-exact");

    // 50 -> 60 on a 60 Hz hum
    {
        HR h;
        h.setParams (params (dh::Mains::Hz50, 3.0));   // params first: prepare applies them at once
        felitronics::test::run (h.prepare (kFs, 512, 1));
        Hum hum; hum.f0 = 60.0; hum.seed (8u);
        std::vector<float> x ((std::size_t) (5 * 48000));
        for (std::size_t i = 0; i < x.size(); ++i) x[i] = (float) hum.at ((std::int64_t) i);
        const std::vector<float> xin = x;
        { float* io[1] { x.data() }; felitronics::test::run (h.process (io, 1, 48000 + 5)); }
        h.setParams (params (dh::Mains::Hz60, 3.0));
        { float* io[1] { x.data() + 48005 }; felitronics::test::run (h.process (io, 1, ramp + 64)); }
        ok (h.running() && std::fabs (h.mainsHz() - 60.0) < 0.02, "50 -> 60: re-armed at 60 Hz after the ramp");
        { float* io[1] { x.data() + 48005 + ramp + 64 }; felitronics::test::run (h.process (io, 1, (int) x.size() - 48005 - ramp - 64)); }
        const double nul = worstNull (x, xin, 4 * 48000, 48000, 60.0);
        ok (nul <= -50.0, "and the 60 Hz null forms (worst harmonic " + std::to_string (nul) + " dB at 3 s)");
    }
    // Auto abstains on an input that carries no hum: dry, bit-exact
    {
        HR h;
        h.setParams (params (dh::Mains::Auto, 0.5, 3.0));   // params first: prepare applies them at once
        felitronics::test::run (h.prepare (kFs, 512, 2));
        Lcg rng { 5u };
        std::vector<float> al ((std::size_t) (4 * 48000)), ar ((std::size_t) (4 * 48000));
        for (int i = 0; i < 4 * 48000; ++i)
        {
            al[(std::size_t) i] = (float) (0.1 * sinAt (329.6276, i, kFs, 0.0)) + 0.001f * rng.bipolar();
            ar[(std::size_t) i] = (float) (0.1 * sinAt (246.9417, i, kFs, 1.0)) + 0.001f * rng.bipolar();
        }
        const std::vector<float> l0 = al, r0 = ar;
        felitronics::test::run (run2 (h, al, ar, 0, (int) al.size()));
        ok (bitsEqual (al, l0) && bitsEqual (ar, r0) && ! h.running() && h.detected() == dh::Mains::Off,
            "Auto on a hum-free input abstains and passes it untouched");
    }
    // Off -> 50 at an arbitrary sample: dry until the next grid boundary, live from it
    {
        HR h;
        h.setParams (params (dh::Mains::Off, 3.0));   // params first: prepare applies them at once
        felitronics::test::run (h.prepare (kFs, 512, 1));
        std::vector<float> x (l.begin(), l.begin() + 4800);
        const std::vector<float> xin = x;
        { float* io[1] { x.data() }; felitronics::test::run (h.process (io, 1, 1000)); }
        h.setParams (params (dh::Mains::Hz50, 3.0));
        { float* io[1] { x.data() + 1000 }; felitronics::test::run (h.process (io, 1, 3800)); }
        bool dryToBoundary = true, liveAfter = false;
        for (int i = 0; i < 1024; ++i) dryToBoundary = dryToBoundary && core::sameBits (x[(std::size_t) i], xin[(std::size_t) i]);
        for (int i = 1024; i < 4800; ++i) liveAfter = liveAfter || ! core::sameBits (x[(std::size_t) i], xin[(std::size_t) i]);
        ok (dryToBoundary && liveAfter, "Off -> 50: dry up to the grid boundary at 1024, live from it");
    }
}

int main()
{
    std::printf ("felitronics::dehum — HumRemover as a stage\n");
    testReslicing();
    testPauseIsSilence();
    testStoppedChannelKeepsModel();
    testPrepareAndAlloc();
    testDigitalBlack();
    testMinuteOfSilence();
    testPoison();
    testModes();
    return felitronics::test::report();
}

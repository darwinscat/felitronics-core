// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// THE AUDIO-TIME MAINTENANCE GRID, and what it buys.
//
// Law 8's denormal flush used to run "once per process() call", so the caller's block boundary decided
// where a numerical event landed. Three consequences, all measured on this module before the fix and all
// asserted here:
//   1. the output was a function of the SLICING — 38522 of 40000 tail samples differed between a
//      whole-file call and one-sample calls of the same stream;
//   2. a whole-file call never flushed INSIDE itself, so the tail sat subnormal — 37678 of 38000 samples,
//      the 10-100x stall law 8 exists to prevent;
//   3. one non-finite input sample was never healed either: 479900 of the next 480000 outputs non-finite
//      on a whole-file call, against 1 when the same stream was fed one sample at a time.
//
// EVERY FIXTURE HERE ASSERTS ITS OWN PRECONDITION, because each of these is a test that would pass on a
// dead signal: a band that filters nothing is trivially slicing-invariant, a tail that never reaches the
// flush threshold trivially has no subnormals, and a poison fixture whose +Inf never reached the state
// trivially recovers. The precondition says the fixture is live; the assertion says the code is right.

#include <felitronics_test.h>
#include <felitronics/core/StateGrid.h>
#include <felitronics/core/Smoother.h>
#include <felitronics/eq/EqBand.h>
#include <felitronics/eq/MatchedBiquad.h>
#include <complex>
#include <felitronics/eq/EqEngine.h>
#include <felitronics/stereo/MonoBass.h>

#include <atomic>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <vector>

static std::atomic<int> g_allocs { 0 };
void* operator new      (std::size_t s) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s ? s : 1); }
void* operator new[]    (std::size_t s) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s ? s : 1); }
void  operator delete   (void* p) noexcept { std::free (p); }
void  operator delete[] (void* p) noexcept { std::free (p); }
void  operator delete   (void* p, std::size_t) noexcept { std::free (p); }
void  operator delete[] (void* p, std::size_t) noexcept { std::free (p); }

using namespace felitronics;
using felitronics::test::ok;
using felitronics::test::approx;
using felitronics::test::group;

static constexpr double kFs = 48000.0;
static constexpr int    kK  = core::StateGrid::kPeriod;
static constexpr float  kSmallestNormal = 1.17549435e-38f;

//==============================================================================
// The fixture: a tone that stops dead into digital silence. The abrupt stop is deliberate and is NOT
// the trap the true-peak work recorded — nothing here is compared against a passband oracle, both sides
// of every comparison see the same edge, and the decaying tail after it is the only place the flush has
// anything to do.
static void programme (std::vector<float>& L, std::vector<float>& R, int n, int toneSamples = 1000)
{
    L.assign ((std::size_t) n, 0.0f); R.assign ((std::size_t) n, 0.0f);
    for (int i = 0; i < toneSamples && i < n; ++i)
    {
        const float v = (float) (0.5 * std::sin (2.0 * M_PI * 1000.0 * (double) i / kFs));
        L[(std::size_t) i] = v; R[(std::size_t) i] = 0.9f * v;
    }
}

static eq::BandParams bell (double freq, double gainDb, double q = 1.0)
{
    eq::BandParams p; p.on = true; p.type = eq::FilterType::Bell;
    for (auto& l : p.lanes) l.on = false;
    auto& st = p.lanes[(std::size_t) eq::Lane::Stereo];
    st.on = true; st.freq = freq; st.Q = q; st.gainDb = gainDb;
    return p;
}

// Render `n` samples through a fresh engine in blocks of `blk` (0 = one whole-file call).
static void renderEq (const eq::BandParams& p, const eq::BandParams* second,
                      std::vector<float>& L, std::vector<float>& R, int n, int blk)
{
    programme (L, R, n);
    eq::EqEngine e;
    ok (e.prepare (kFs, n, 2), "renderEq: engine prepared");
    e.setBand (0, p);
    if (second != nullptr) e.setBand (0, *second);
    const int step = blk > 0 ? blk : n;
    for (int off = 0; off < n; off += step)
    {
        const int m = (step < n - off) ? step : n - off;
        float* ch[2] = { L.data() + off, R.data() + off };
        e.process (ch, 2, m);
    }
}

//==============================================================================
static void testGridArithmetic()
{
    group ("core::StateGrid — the phase arithmetic itself");
    core::StateGrid g;
    ok (g.phase() == 0, "a fresh grid is at phase 0");
    ok (g.segment (1000000) == kK, "segment() never hands out more than one period");
    ok (g.segment (5) == 5, "segment() never hands out more than the caller has left");

    // Walk a whole period in ragged pieces: the boundary must land on the kK-th sample and nowhere else.
    int taken = 0, boundaries = 0;
    const int pieces[] = { 7, 1, 30, 100, 4, 64, 13 };
    for (int p : pieces)
    {
        int left = p;
        while (left > 0)
        {
            const int seg = g.segment (left);
            ok (seg >= 1 && seg <= kK, "segment() is in [1, kPeriod]");
            taken += seg; left -= seg;
            if (g.advance (seg)) { ++boundaries; ok (taken % kK == 0, "a boundary falls on a multiple of kPeriod"); }
            ok (g.phase() >= 0 && g.phase() < kK, "phase stays inside [0, kPeriod)");
        }
    }
    ok (boundaries == taken / kK, "exactly one boundary per period of audio, however it was cut");

    core::StateGrid h;
    h.advance (h.segment (10));
    ok (h.phase() == 10, "phase tracks the audio consumed");
    h.reset();
    ok (h.phase() == 0, "reset() re-anchors the grid");
    h.skip (kK * 1000 + 5);
    ok (h.phase() == 5, "skip() is modular and survives a large jump");
    h.skip (0);
    ok (h.phase() == 5, "skip(0) moves nothing");
}

//==============================================================================
static void testEqSlicingInvariance()
{
    group ("eq::EqEngine — the same stream, any slicing, bit-identical output");
    const int n = 40000;
    const int blocks[] = { 0, 4096, 512, 100, 64, 7, 1 };

    struct Case { const char* name; eq::BandParams p; };
    // The 0 dB bell is here on purpose: P18 F18 measured that a matched bell at 0 dB is NOT bit-exactly
    // unity — it holds state and rings at ~1.4e-10 — so its state sits AT the flush threshold and the
    // divergence used to land in the middle of the TONE, not only in the tail.
    const Case cases[] = { { "Bell +6 dB", bell (1000.0, 6.0) },
                           { "Bell  0 dB", bell (1000.0, 0.0) },
                           { "Bell -20 dB Q8", bell (4000.0, -20.0, 8.0) } };

    for (const Case& c : cases)
    {
        std::vector<float> ref[2];
        renderEq (c.p, nullptr, ref[0], ref[1], n, blocks[0]);

        // PRECONDITION 1 — the band is LIVE: it changed the signal.
        std::vector<float> dryL, dryR; programme (dryL, dryR, n);
        bool changed = false;
        for (int i = 0; i < 1000; ++i) if (ref[0][(std::size_t) i] != dryL[(std::size_t) i]) { changed = true; break; }
        ok (changed, std::string (c.name) + ": PRECONDITION the band actually filters the tone");

        // PRECONDITION 2 — the FLUSH is live: the tail decays past the threshold and lands on exact zero,
        // which only happens because something zeroed it. Without a flush the tail is a subnormal creep.
        long long zeros = 0, subs = 0;
        for (int i = 2000; i < n; ++i)
        {
            const float v = ref[0][(std::size_t) i];
            if (v == 0.0f) ++zeros;
            else if (std::fabs (v) < kSmallestNormal) ++subs;
        }
        ok (zeros > n / 2, std::string (c.name) + ": PRECONDITION the tail reaches EXACT zero (flush fired)");
        ok (subs == 0, std::string (c.name) + ": no subnormal sample survives a WHOLE-FILE call");

        for (int blk : blocks)
        {
            std::vector<float> L, R;
            renderEq (c.p, nullptr, L, R, n, blk);
            long long diff = 0;
            for (int i = 0; i < n; ++i)
                if (L[(std::size_t) i] != ref[0][(std::size_t) i] || R[(std::size_t) i] != ref[1][(std::size_t) i]) ++diff;
            ok (diff == 0, std::string (c.name) + ": block " + std::to_string (blk) + " is bit-identical to the whole-file call");
        }
    }
}

//==============================================================================
static void testEqSlicingInvarianceUnderARamp()
{
    group ("eq::EqEngine — invariance HOLDS while a parameter ramp is in flight");
    const int n = 8000;
    const int blocks[] = { 0, 4096, 512, 100, 64, 7, 1 };
    const eq::BandParams from = bell (500.0, 0.0);
    const eq::BandParams to   = bell (5000.0, 9.0);

    std::vector<float> ref[2];
    renderEq (from, &to, ref[0], ref[1], n, blocks[0]);

    // PRECONDITION — the ramp is REALLY in flight during the render, not settled before the first sample.
    // Read the band's own response at the target frequency at three instants: it must still be moving.
    {
        eq::EqBand b;
        ok (b.prepare (kFs, 2), "ramp fixture: band prepared");
        b.setParams (from);
        b.setParams (to);
        std::vector<float> buf[2]; buf[0].assign (64, 0.0f); buf[1].assign (64, 0.0f);
        float* ch[2] = { buf[0].data(), buf[1].data() };
        const double w = 2.0 * M_PI * 5000.0 / kFs;
        b.processBlock (ch, 2, 64);
        const double a1 = std::abs (b.response (w));
        for (int k = 0; k < 8; ++k) b.processBlock (ch, 2, 64);
        const double a2 = std::abs (b.response (w));
        for (int k = 0; k < 40; ++k) b.processBlock (ch, 2, 64);
        const double a3 = std::abs (b.response (w));
        ok (std::fabs (a2 - a1) > 1e-3 && std::fabs (a3 - a2) > 1e-3,
            "PRECONDITION the design is still MOVING across the render (|H| " + std::to_string (a1) + " -> "
            + std::to_string (a2) + " -> " + std::to_string (a3) + ")");
    }

    for (int blk : blocks)
    {
        std::vector<float> L, R;
        renderEq (from, &to, L, R, n, blk);
        long long diff = 0; double worst = 0.0;
        for (int i = 0; i < n; ++i)
        {
            if (L[(std::size_t) i] != ref[0][(std::size_t) i]) { ++diff; const double d = std::fabs ((double) L[(std::size_t) i] - (double) ref[0][(std::size_t) i]); if (d > worst) worst = d; }
            if (R[(std::size_t) i] != ref[1][(std::size_t) i]) ++diff;
        }
        ok (diff == 0, "ramping: block " + std::to_string (blk) + " bit-identical (worst " + std::to_string (worst) + ")");
    }
}

//==============================================================================
static void testPoisonWindowIsBounded()
{
    group ("eq::EqEngine — one non-finite sample is healed within ONE grid period, at any block size");
    const int n = 20000;
    for (int blk : { 0, 8192, 4096, 512, 100, 64, 7, 1 })
    {
        eq::EqEngine e;
        ok (e.prepare (kFs, n, 2), "poison fixture: engine prepared");
        e.setBand (0, bell (1000.0, 6.0, 2.0));
        std::vector<float> L ((std::size_t) n), R ((std::size_t) n);
        for (int i = 0; i < n; ++i)
        {
            const float v = (float) (0.25 * std::sin (2.0 * M_PI * 1000.0 * (double) i / kFs));
            L[(std::size_t) i] = v; R[(std::size_t) i] = v;
        }
        L[100] = std::numeric_limits<float>::infinity();
        const int step = blk > 0 ? blk : n;
        for (int off = 0; off < n; off += step)
        {
            const int m = (step < n - off) ? step : n - off;
            float* ch[2] = { L.data() + off, R.data() + off };
            e.process (ch, 2, m);
        }
        long long bad = 0; int last = -1;
        for (int i = 0; i < n; ++i) if (! std::isfinite (L[(std::size_t) i])) { ++bad; last = i; }
        // PRECONDITION — the poison REACHED the state. A fixture whose +Inf was swallowed before the
        // recursion would report a perfect recovery while testing nothing.
        ok (bad >= 1, "block " + std::to_string (blk) + ": PRECONDITION the +Inf actually poisoned the filter");
        ok (bad <= kK, "block " + std::to_string (blk) + ": recovery within one grid period (" + std::to_string (bad) + " non-finite)");
        ok (last < 100 + kK, "block " + std::to_string (blk) + ": last bad sample inside the period that carried the poison");
    }
    // ... and a host with a block SHORTER than a period keeps the immediate recovery it had before: the
    // poison half still runs at the end of every call.
    {
        eq::EqEngine e;
        ok (e.prepare (kFs, 4000, 2), "poison fixture (short blocks): prepared");
        e.setBand (0, bell (1000.0, 6.0, 2.0));
        std::vector<float> L (4000, 0.0f), R (4000, 0.0f);
        for (int i = 0; i < 4000; ++i) { const float v = (float) (0.25 * std::sin (2.0 * M_PI * 1000.0 * (double) i / kFs)); L[(std::size_t) i] = v; R[(std::size_t) i] = v; }
        L[100] = std::numeric_limits<float>::infinity();
        for (int off = 0; off < 4000; ++off) { float* ch[2] = { L.data() + off, R.data() + off }; e.process (ch, 2, 1); }
        long long bad = 0;
        for (int i = 0; i < 4000; ++i) if (! std::isfinite (L[(std::size_t) i])) ++bad;
        ok (bad == 1, "one-sample calls still heal in ONE sample, not one period (" + std::to_string (bad) + ")");
    }
}

//==============================================================================
static void testResetIsAStreamRestart()
{
    group ("eq::EqBand::reset() — a stream restart, not a filter-state wipe (P6 F13)");
    const int n = 4000;
    const eq::BandParams a = bell (900.0, 0.0);
    const eq::BandParams b = bell (900.0, 9.0);

    auto run = [&] (bool viaReset, const eq::BandParams* writeAfter, std::vector<float>& out)
    {
        eq::EqBand band;
        ok (band.prepare (kFs, 2), "reset fixture: prepared");
        // Tone for a third, then digital silence: the ramp is audible in the first part and the FLUSH
        // has something to do in the second. With a full-length tone (the first version of this
        // fixture) the grid phase was unobservable and a mutation that stopped reset() re-anchoring it
        // passed — the tail is what makes the phase visible.
        std::vector<float> L, R; programme (L, R, n, n / 3);
        if (viaReset)
        {
            band.setParams (a);
            band.setParams (b);                                // start a 30 ms glide
            std::vector<float> wL, wR; programme (wL, wR, 480, 480);
            float* wc[2] = { wL.data(), wR.data() };
            band.processBlock (wc, 2, 480);                    // render 10 ms of it
            band.reset();
            if (writeAfter != nullptr) band.setParams (*writeAfter);
        }
        else
        {
            band.setParams (writeAfter != nullptr ? *writeAfter : b);   // a freshly prepared band, one write
        }
        float* ch[2] = { L.data(), R.data() };
        band.processBlock (ch, 2, n);
        out = L;
    };

    // PRECONDITION — the two parameter points are audibly different, so "same output" is a real claim.
    {
        std::vector<float> x, y;
        run (false, &a, x); run (false, &b, y);
        double worst = 0.0;
        for (int i = 0; i < n; ++i) worst = std::fmax (worst, std::fabs ((double) x[(std::size_t) i] - (double) y[(std::size_t) i]));
        ok (worst > 0.1, "PRECONDITION 0 dB and +9 dB differ by " + std::to_string (worst) + " full scale");
    }

    std::vector<float> viaReset, fresh;
    run (true, nullptr, viaReset);
    run (false, nullptr, fresh);
    long long diff = 0; double worst = 0.0;
    for (int i = 0; i < n; ++i)
        if (viaReset[(std::size_t) i] != fresh[(std::size_t) i])
        { ++diff; worst = std::fmax (worst, std::fabs ((double) viaReset[(std::size_t) i] - (double) fresh[(std::size_t) i])); }
    ok (diff == 0, "a band reset mid-ramp renders bit-identically to a freshly prepared one (" + std::to_string (diff)
                   + " differ, worst " + std::to_string (worst) + ")");

    // The half that is easy to miss: a write ARRIVING AFTER the reset must snap, not ramp.
    std::vector<float> afterReset, freshOther;
    const eq::BandParams other = bell (3000.0, -12.0);
    run (true, &other, afterReset);
    run (false, &other, freshOther);
    diff = 0; worst = 0.0;
    for (int i = 0; i < n; ++i)
        if (afterReset[(std::size_t) i] != freshOther[(std::size_t) i])
        { ++diff; worst = std::fmax (worst, std::fabs ((double) afterReset[(std::size_t) i] - (double) freshOther[(std::size_t) i])); }
    ok (diff == 0, "a DIFFERENT write after reset() snaps, exactly as the first write after prepare() does ("
                   + std::to_string (diff) + " differ, worst " + std::to_string (worst) + ")");
}

//==============================================================================
static void testClearAudioStateIsAStop()
{
    group ("eq::EqBand::clearAudioState() — a STOP: clears history, leaves the ramp alone");
    const int n = 4000;
    const eq::BandParams a = bell (900.0, 0.0);
    const eq::BandParams b = bell (900.0, 9.0);

    eq::EqBand stopped, running;
    ok (stopped.prepare (kFs, 2) && running.prepare (kFs, 2), "stop fixture: prepared");
    for (eq::EqBand* band : { &stopped, &running })
    {
        band->setParams (a);
        band->setParams (b);
        std::vector<float> wL, wR; programme (wL, wR, 480, 480);
        float* wc[2] = { wL.data(), wR.data() };
        band->processBlock (wc, 2, 480);
    }
    // PRECONDITION — the glide is genuinely unfinished at the stop, or "keeps ramping" means nothing.
    const double w = 2.0 * M_PI * 900.0 / kFs;
    const double midway = std::abs (running.response (w));
    ok (midway > 1.05 && midway < 2.7, "PRECONDITION the glide is mid-flight (|H| = " + std::to_string (midway) + ")");

    stopped.clearAudioState();
    ok (std::fabs (std::abs (stopped.response (w)) - midway) < 1e-12,
        "clearAudioState() leaves the design exactly where it was");

    std::vector<float> L1, R1, L2, R2; programme (L1, R1, n, n / 3); L2 = L1; R2 = R1;
    float* c1[2] = { L1.data(), R1.data() };
    float* c2[2] = { L2.data(), R2.data() };
    stopped.processBlock (c1, 2, n);
    running.processBlock (c2, 2, n);
    // The two differ only where the cleared filter memory shows (the first samples); the RAMP must be the
    // same, so the two converge and the far end is bit-identical.
    long long tailDiff = 0;
    for (int i = n / 2; i < n; ++i) if (L1[(std::size_t) i] != L2[(std::size_t) i]) ++tailDiff;
    ok (tailDiff == 0, "after the cleared history has decayed the two runs agree bit-for-bit — the glide was not snapped");
}

//==============================================================================
static void testMonoBassSlicingInvariance()
{
    group ("stereo::MonoBass — the same stream, any slicing, bit-identical output");
    const int n = 40000;
    const int blocks[] = { 0, 4096, 512, 100, 64, 7, 1 };
    for (float width : { 0.0f, 0.5f, 1.0f })
    {
        std::vector<float> ref[2];
        auto render = [&] (std::vector<float>& L, std::vector<float>& R, int blk)
        {
            programme (L, R, n);
            stereo::MonoBass m; m.prepare (kFs, n, 2);
            m.setParams ({ true, 150.0f, width });
            const int step = blk > 0 ? blk : n;
            for (int off = 0; off < n; off += step)
            {
                const int k = (step < n - off) ? step : n - off;
                float* ch[2] = { L.data() + off, R.data() + off };
                m.process (ch, 2, k);
            }
        };
        render (ref[0], ref[1], 0);

        // PRECONDITION — the stage is live at this width. At width 1.0 it settles into bypass, and that
        // is the case whose settle-to-bypass edge used to be call-clocked, so it must still be swept.
        std::vector<float> dryL, dryR; programme (dryL, dryR, n);
        long long touched = 0;
        for (int i = 0; i < 1000; ++i) if (ref[0][(std::size_t) i] != dryL[(std::size_t) i]) ++touched;
        ok (width >= 1.0f || touched > 100,
            "width " + std::to_string (width) + ": PRECONDITION the stage moves the signal (" + std::to_string (touched) + ")");

        for (int blk : blocks)
        {
            std::vector<float> L, R; render (L, R, blk);
            long long diff = 0;
            for (int i = 0; i < n; ++i)
                if (L[(std::size_t) i] != ref[0][(std::size_t) i] || R[(std::size_t) i] != ref[1][(std::size_t) i]) ++diff;
            ok (diff == 0, "width " + std::to_string (width) + ", block " + std::to_string (blk) + ": bit-identical");
        }
    }
}

//==============================================================================
// The ramp's ABSOLUTE rate, against an oracle built out of the primitives rather than out of the band.
// Every other check here compares one EqBand render with another, so all of them pass on a band whose
// glide runs at the wrong speed — measured: a mutation advancing by kPeriod-1 instead of kPeriod
// survived the whole suite. This is the only test that can see it.
static void testRampFollowsItsOwnDesign()
{
    group ("eq::EqBand — the glide advances by EXACTLY one grid period per tick (independent oracle)");
    const double smoothMs = 30.0;
    const double f0 = 500.0, f1 = 5000.0, q = 1.4;
    // +12 dB, NOT 0 dB. A matched bell at 0 dB is very nearly unity at every frequency, so |H| moves by
    // ~1e-14 while the design frequency travels a decade — an oracle built on it agrees with a band
    // gliding at the WRONG SPEED, and that is exactly how the first version of this test let a mutation
    // advancing by kPeriod-1 through. The precondition below now asserts that the ORACLE'S OWN READOUT
    // travels, not merely that the smoother behind it does.
    const double gain = 12.0;

    eq::EqBand b;
    ok (b.prepare (kFs, 2, smoothMs), "ramp oracle: band prepared");
    eq::BandParams a = bell (f0, gain, q); b.setParams (a);           // first write snaps
    eq::BandParams t = bell (f1, gain, q); b.setParams (t);           // now glide

    // The oracle: the same primitive, driven by hand, one tick of kPeriod samples at a time.
    core::Smoother ref; ref.prepare (kFs, smoothMs); ref.snap (f0); ref.setTarget (f1);

    std::vector<float> buf[2];
    buf[0].assign ((std::size_t) kK, 0.0f); buf[1].assign ((std::size_t) kK, 0.0f);
    float* ch[2] = { buf[0].data(), buf[1].data() };

    const double w = 2.0 * M_PI * 2000.0 / kFs;
    bool agreed = true;
    double worst = 0.0, minH = 1e300, maxH = -1e300; int worstTick = -1; long long moving = 0;
    for (int tickNo = 1; tickNo <= 2000; ++tickNo)
    {
        b.processBlock (ch, 2, kK);
        ref.advance (kK);

        eq::BandParams oracleP = bell (ref.value(), gain, q);
        const eq::BandDesign d = eq::designBand (oracleP, kFs);
        std::complex<double> h { 1.0, 0.0 };
        for (int i = 0; i < d.n; ++i) h *= eq::evalCoeffs (d.sec[i], w);
        const double want = std::abs (h);
        const double got  = std::abs (b.response (w));
        minH = std::fmin (minH, want); maxH = std::fmax (maxH, want);
        // Compare only WHILE THE GLIDE IS MOVING. `Smoother::settled()` carries an epsilon of 1e-7 in
        // the parameter's own unit, and the band deliberately stops redesigning once it fires — so past
        // that point the oracle keeps creeping the last 1e-7 Hz and the band, correctly, does not. That
        // is the band's contract, not a disagreement, and comparing there would only be measuring the
        // epsilon. The speed is observable exactly where it matters: before the smoother lands.
        if (! ref.settled())
        {
            const double err = std::fabs (got - want);
            if (err > worst) { worst = err; worstTick = tickNo; }
            if (err > 1e-12) agreed = false;
            ++moving;
        }
    }
    // PRECONDITION — the ORACLE'S READOUT travels, by a lot. Not "the smoother moved": the smoother can
    // travel a decade while the quantity being compared sits still, and then the comparison proves
    // nothing about the speed of anything.
    ok (maxH - minH > 0.5, "PRECONDITION the oracle's own |H| travels across the render ("
                           + std::to_string (minH) + " -> " + std::to_string (maxH) + ")");
    ok (ref.settled(), "PRECONDITION the oracle ARRIVES at the target (" + std::to_string (ref.value()) + " Hz)");
    ok (moving > 50, "PRECONDITION the comparison ran over " + std::to_string (moving) + " MOVING ticks");
    ok (agreed, "the band's design tracks a hand-driven Smoother tick for tick (worst |dH| = "
                + std::to_string (worst) + " at tick " + std::to_string (worstTick) + ", over "
                + std::to_string (moving) + " ticks)");

    // ... and the SETTLED value is designed, not the one-tick-early value. This is the half a suite
    // without it misses: the tick on which a smoother lands reads `settled()`, so the last redesign has
    // to be carried over one more tick.
    eq::BandParams finalP = bell (f1, gain, q);
    const eq::BandDesign fd = eq::designBand (finalP, kFs);
    std::complex<double> hf { 1.0, 0.0 };
    for (int i = 0; i < fd.n; ++i) hf *= eq::evalCoeffs (fd.sec[i], w);
    // The settled design, to a tolerance the settle epsilon justifies rather than one picked by hand:
    // the band's last redesign happened within `Smoother::settled()`'s 1e-7 of the target frequency.
    ok (std::fabs (std::abs (b.response (w)) - std::abs (hf)) < 1e-6,
        "the FINAL settled design is the target's, within the smoother's own settle epsilon");

    // THE SAME ORACLE ON A MONO LANE. The ST columns and the L/R/M/S lanes advance through different
    // code, so pinning one says nothing about the other.
    {
        eq::EqBand sb;
        ok (sb.prepare (kFs, 2, smoothMs), "ramp oracle (Side lane): band prepared");
        auto sideBand = [&] (double f) {
            eq::BandParams p; p.on = true; p.type = eq::FilterType::Bell;
            for (auto& l : p.lanes) l.on = false;
            auto& sd = p.lanes[(std::size_t) eq::Lane::Side];
            sd.on = true; sd.freq = f; sd.Q = q; sd.gainDb = gain;
            return p;
        };
        sb.setParams (sideBand (f0));
        sb.setParams (sideBand (f1));
        core::Smoother sref; sref.prepare (kFs, smoothMs); sref.snap (f0); sref.setTarget (f1);
        bool sAgreed = true; double sWorst = 0.0, sMin = 1e300, sMax = -1e300; long long sMoving = 0;
        for (int tickNo = 1; tickNo <= 2000; ++tickNo)
        {
            sb.processBlock (ch, 2, kK);
            sref.advance (kK);
            eq::BandParams op = sideBand (sref.value());
            op.lanes[(std::size_t) eq::Lane::Stereo] = op.lanes[(std::size_t) eq::Lane::Side];   // designBand reads the primary slot
            const eq::BandDesign d = eq::designBand (op, kFs);
            std::complex<double> h { 1.0, 0.0 };
            for (int i = 0; i < d.n; ++i) h *= eq::evalCoeffs (d.sec[i], w);
            const double want = std::abs (h);
            const double got  = std::abs (sb.response (w, eq::Axis::Side));
            sMin = std::fmin (sMin, want); sMax = std::fmax (sMax, want);
            if (! sref.settled())
            {
                const double err = std::fabs (got - want);
                if (err > sWorst) sWorst = err;
                if (err > 1e-12) sAgreed = false;
                ++sMoving;
            }
        }
        ok (sMax - sMin > 0.5, "PRECONDITION the Side-lane oracle's |H| travels (" + std::to_string (sMin) + " -> " + std::to_string (sMax) + ")");
        ok (sMoving > 50, "PRECONDITION the Side-lane comparison ran over " + std::to_string (sMoving) + " MOVING ticks");
        ok (sAgreed, "the SIDE lane's design tracks the same hand-driven Smoother (worst |dH| = " + std::to_string (sWorst)
                     + ", over " + std::to_string (sMoving) + " ticks)");
    }
}

//==============================================================================
// The grid is re-anchored by reset() and only by reset(). Without this a stream restart inherits the
// previous stream's phase, so the flush lands somewhere else and two renders of the same programme
// differ — invisible to any fixture whose tail never reaches the flush threshold.
static void testGridReAnchorsOnReset()
{
    group ("the maintenance grid is re-anchored by reset(), in both modules");
    const int n = 20000;
    const int warm = 37;                     // deliberately NOT a multiple of the period

    {   // eq::EqBand
        eq::EqBand fresh, used;
        ok (fresh.prepare (kFs, 2) && used.prepare (kFs, 2), "grid anchor fixture: prepared");
        fresh.setParams (bell (1000.0, 6.0)); used.setParams (bell (1000.0, 6.0));
        std::vector<float> wL, wR; programme (wL, wR, warm, warm);
        float* wc[2] = { wL.data(), wR.data() };
        used.processBlock (wc, 2, warm);     // leave the phase at 37
        used.reset(); fresh.reset();
        std::vector<float> aL, aR, bL, bR; programme (aL, aR, n); bL = aL; bR = aR;
        float* ac[2] = { aL.data(), aR.data() };
        float* bc[2] = { bL.data(), bR.data() };
        fresh.processBlock (ac, 2, n);
        used.processBlock (bc, 2, n);
        long long diff = 0;
        for (int i = 0; i < n; ++i) if (aL[(std::size_t) i] != bL[(std::size_t) i] || aR[(std::size_t) i] != bR[(std::size_t) i]) ++diff;
        // PRECONDITION — the tail actually reaches the flush, or the phase is unobservable.
        long long zeros = 0; for (int i = 2000; i < n; ++i) if (aL[(std::size_t) i] == 0.0f) ++zeros;
        ok (zeros > n / 2, "PRECONDITION the reference tail reaches exact zero (" + std::to_string (zeros) + ")");
        ok (diff == 0, "EqBand: a band used for " + std::to_string (warm) + " samples then reset renders identically to a fresh one");
    }
    {   // stereo::MonoBass
        stereo::MonoBass fresh, used;
        fresh.prepare (kFs, n, 2); used.prepare (kFs, n, 2);
        fresh.setParams ({ true, 150.0f, 0.0f }); used.setParams ({ true, 150.0f, 0.0f });
        std::vector<float> wL, wR; programme (wL, wR, warm, warm);
        float* wc[2] = { wL.data(), wR.data() };
        used.process (wc, 2, warm);
        used.reset(); fresh.reset();
        std::vector<float> aL, aR, bL, bR; programme (aL, aR, n); bL = aL; bR = aR;
        float* ac[2] = { aL.data(), aR.data() };
        float* bc[2] = { bL.data(), bR.data() };
        fresh.process (ac, 2, n);
        used.process (bc, 2, n);
        long long diff = 0;
        for (int i = 0; i < n; ++i) if (aL[(std::size_t) i] != bL[(std::size_t) i] || aR[(std::size_t) i] != bR[(std::size_t) i]) ++diff;
        long long zeros = 0; for (int i = 2000; i < n; ++i) if (aL[(std::size_t) i] == 0.0f) ++zeros;
        ok (zeros > n / 2, "PRECONDITION the MonoBass tail reaches exact zero (" + std::to_string (zeros) + ")");
        ok (diff == 0, "MonoBass: same, after " + std::to_string (warm) + " samples and a reset");
    }
}

//==============================================================================
// A HALF-poisoned state, which is the case the atomic contract exists for and which an +Inf fixture
// never produces: a finite input can drive one word non-finite while the other is still a large finite
// number, and a per-word heal then leaves the filter to re-poison itself from the half it kept.
static void testPoisonHealIsAtomic()
{
    group ("healPoison() is ATOMIC per channel — a half-healed filter re-poisons itself");
    {
        // SEARCH for the half-poisoned state rather than assume one fixture reaches it — Svf.h records
        // that it exists on this filter (a finite 1e37 sine at Q 40 leaves z1 at +Inf while z2 is a
        // finite -3.27e38) and the whole point of the atomic contract is that a per-word heal then
        // resurrects the filter from the word it kept.
        eq::Biquad b; bool half = false; double foundAmp = 0.0, foundQ = 0.0, foundF = 0.0;
        for (double amp : { 1.0e35, 1.0e36, 1.0e37, 1.0e38, 3.0e38 })
          for (double qq : { 0.7071, 2.0, 8.0, 40.0 })
            for (double ff : { 200.0, 2000.0, 7000.0, 18000.0 })
              for (int typeI = 0; typeI < 3 && ! half; ++typeI)
              {
                  if (half) break;
                  eq::Biquad t;
                  t.setCoeffs (typeI == 0 ? eq::matched::peakingDb (ff, kFs, qq, 24.0)
                             : typeI == 1 ? eq::matched::bandpass (ff, kFs, qq)
                                          : eq::matched::lowpass (ff, kFs, qq));
                  t.reset();
                  for (int i = 0; i < 512; ++i)
                  {
                      t.processSample ((float) (amp * std::sin (2.0 * M_PI * ff * (double) i / kFs)));
                      if (std::isfinite (t.z1) != std::isfinite (t.z2))
                      { half = true; b = t; foundAmp = amp; foundQ = qq; foundF = ff; break; }
                  }
              }
        ok (half, "PRECONDITION a finite input reached a HALF non-finite state (amp " + std::to_string (foundAmp)
                  + ", f " + std::to_string (foundF) + ", Q " + std::to_string (foundQ) + ")");
        if (half)
        {
            const float kept = std::isfinite (b.z1) ? b.z1 : b.z2;
            b.healPoison();
            ok (b.z1 == 0.0f && b.z2 == 0.0f,
                "healPoison() takes BOTH words, not only the bad one (the surviving word was " + std::to_string (kept) + ")");
            bool clean = true;
            for (int i = 0; i < 64; ++i) if (! std::isfinite (b.processSample (0.0f))) clean = false;
            ok (clean, "the healed filter emits finite samples into silence instead of re-poisoning itself");
        }
    }
    {   // ... and the same through MonoBass's public surface, which is where the per-call heal lives.
        stereo::MonoBass m; m.prepare (kFs, 512, 2);
        m.setParams ({ true, 150.0f, 0.0f });
        std::vector<float> L (512, 0.0f), R (512, 0.0f);
        for (int i = 0; i < 512; ++i) { L[(std::size_t) i] = 0.3f; R[(std::size_t) i] = -0.3f; }
        L[10] = std::numeric_limits<float>::infinity();
        float* ch[2] = { L.data(), R.data() };
        m.process (ch, 2, 512);
        long long bad = 0; for (float v : L) if (! std::isfinite (v)) ++bad;
        ok (bad >= 1, "PRECONDITION the +Inf poisoned MonoBass's crossover (" + std::to_string (bad) + " non-finite)");
        std::vector<float> L2 (512, 0.05f), R2 (512, -0.05f);
        float* ch2[2] = { L2.data(), R2.data() };
        m.process (ch2, 2, 512);
        long long bad2 = 0; for (float v : L2) if (! std::isfinite (v)) ++bad2;
        ok (bad2 == 0, "the NEXT call is clean — the per-call poison heal ran (" + std::to_string (bad2) + " non-finite)");
    }
}

//==============================================================================
static void testSmootherStrideIsExact()
{
    group ("core::Smoother — the memoised coeffⁿ is the same number, not a nearby one");
    // The cache is keyed on n, so a fixture that always asks for the SAME n never tests the key —
    // a mutation that filled the cache once and then ignored n survived exactly that way. Alternate.
    for (double ms : { 1.0, 30.0, 250.0 })
    {
        const int ns[] = { 64, 1, 64, 7, 64, 100, 64, 64, 3 };
        core::Smoother cached;
        cached.prepare (kFs, ms); cached.snap (100.0); cached.setTarget (4000.0);
        // Independent oracle: the same glide stepped ONE SAMPLE AT A TIME through next(), which never
        // touches the cache at all. Exponential decay composes, so after the same number of samples the
        // two must agree to rounding — and, more sharply, an n-blind cache diverges by whole periods.
        core::Smoother perSample;
        perSample.prepare (kFs, ms); perSample.snap (100.0); perSample.setTarget (4000.0);
        double worst = 0.0; long long total = 0;
        for (int rep = 0; rep < 12; ++rep)
            for (int n : ns)
            {
                cached.advance (n);
                for (int i = 0; i < n; ++i) perSample.next();
                total += n;
                worst = std::fmax (worst, std::fabs (cached.value() - perSample.value()));
            }
        ok (worst < 1e-6, std::to_string (ms) + " ms: the memoised stride tracks a per-sample glide over "
                          + std::to_string (total) + " samples (worst " + std::to_string (worst) + ")");
        // ... and a coefficient change must drop the cache, or the next advance uses the old glide
        core::Smoother c; c.prepare (kFs, ms); c.snap (100.0); c.setTarget (4000.0);
        c.advance (64);
        c.setTimeMs (ms * 4.0);
        core::Smoother d; d.prepare (kFs, ms * 4.0); d.snap (c.value()); d.setTarget (4000.0);
        c.advance (64); d.advance (64);
        ok (c.value() == d.value(), "setTimeMs() invalidates the memoised power (" + std::to_string (ms) + " ms)");
    }
}

//==============================================================================
int main()
{
    std::printf ("felitronics::core::StateGrid + eq/stereo audio-time maintenance tests\n");
    testGridArithmetic();
    testEqSlicingInvariance();
    testEqSlicingInvarianceUnderARamp();
    testRampFollowsItsOwnDesign();
    testGridReAnchorsOnReset();
    testPoisonWindowIsBounded();
    testPoisonHealIsAtomic();
    testResetIsAStreamRestart();
    testClearAudioStateIsAStop();
    testMonoBassSlicingInvariance();
    testSmootherStrideIsExact();

    group ("RT-safety");
    {
        eq::EqBand b; b.prepare (kFs, 2);
        b.setParams (bell (1000.0, 6.0));
        std::vector<float> v[2]; float* ch[2] {};
        for (int c = 0; c < 2; ++c) { v[c].assign (4096, 0.1f); ch[c] = v[c].data(); }
        b.processBlock (ch, 2, 4096);
        const int before = g_allocs.load();
        for (int k = 0; k < 20; ++k) b.processBlock (ch, 2, 4096);
        b.reset();
        b.clearAudioState();
        for (int k = 0; k < 20; ++k) b.processBlock (ch, 2, 1);
        felitronics::test::okNoAlloc (g_allocs.load() == before, "no allocation across the segment loop, reset() and clearAudioState()");
    }

    return felitronics::test::report();
}

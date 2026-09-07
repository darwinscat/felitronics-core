// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// The same defect one floor up from the filters, where fixing eq::EqBand alone would have left it live.
//
// Both objects here own PER-CHANNEL detector state that advances only for c < nc, while the envelope and
// the gain-reduction follower it feeds are SHARED and keep integrating. So the lane never stops — the
// Stereo lane runs at any channel count — and the per-lane drop cannot see the problem: a channel that
// leaves and returns hands the linked max() a probe column frozen from before the gap, and the shared
// envelope takes it as present energy. Measured before this suite, on DIGITAL SILENCE and against a
// never-left reference reading 0.000: 11.97 dB of unearned gain reduction in LaneDynamics, and in
// DynamicEqBand an 0.388 audio tail eleven samples in with 8.85 dB of delta behind it.
//
// The shared half is deliberately NOT touched. A detector linked across channels is supposed to follow
// whichever channels are actually there; only the frozen per-channel columns are a lie.

#include <felitronics_test.h>
#include <felitronics/dynamiceq/LaneDynamics.h>
#include <felitronics/dynamiceq/DynamicEqBand.h>
#include <felitronics/core/Math.h>

#include <atomic>
#include <cmath>
#include <cstdlib>
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

static double peakOf (const std::vector<float>& v)
{
    double m = 0.0; for (float x : v) m = std::fmax (m, (double) std::fabs (x)); return m;
}

static bool bitEqual (const std::vector<float>& a, const std::vector<float>& b)
{
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i) if (! (a[i] == b[i])) return false;
    return true;
}

static eq::BandParams dynPoint (double freq, double Q, double thrDb = -60.0)
{
    eq::BandParams p;
    p.on = true; p.type = eq::FilterType::Bell;
    eq::LaneParams& st = p.lane (eq::Lane::Stereo);
    st.on = true; st.freq = freq; st.Q = Q; st.gainDb = 0.0;
    p.dyn.on = true; p.dyn.rangeDb = -12.0;
    p.dyn.thrAuto = false; p.dyn.thrDb = thrDb;   // absolute, so the fixture is not at the mercy of the estimator
    return p;
}

//==================================================================================================
static void testLaneDynamicsProbe()
{
    group ("LaneDynamics: a returned channel earns no gain reduction out of digital silence");

    struct Case { double f, Q; const char* name; };
    const Case cases[3] { { 1000.0, 2.0, "1 kHz / Q2" }, { 100.0, 10.0, "100 Hz / Q10" }, { 5000.0, 1.0, "5 kHz / Q1" } };

    for (const Case& cs : cases)
    {
        const int block = 64;
        const eq::BandParams p = dynPoint (cs.f, cs.Q);

        eq::EqBand band; band.prepare (kFs, 2); band.setParams (p);
        dynamiceq::LaneDynamics dyn; felitronics::test::run (dyn.prepare (kFs, 2)); dyn.setParams (p);

        std::vector<float> L ((std::size_t) block), R ((std::size_t) block), sl ((std::size_t) block), sr ((std::size_t) block);
        float* aud[2] { L.data(), R.data() };
        const float* sc[2] { sl.data(), sr.data() };

        // Charge: channel 0 silent, channel 1 at full scale on the band centre. Only ch1's probe column
        // ends up holding energy, which is exactly the column the excursion will freeze.
        long ph = 0;
        double deepest = 0.0;
        for (int done = 0; done < (int) kFs; done += block)
        {
            for (int i = 0; i < block; ++i, ++ph)
            {
                const float v = (float) std::sin (2.0 * core::kPi * cs.f * (double) ph / kFs);
                L[(std::size_t) i] = sl[(std::size_t) i] = 0.0f;
                R[(std::size_t) i] = sr[(std::size_t) i] = v;
            }
            felitronics::test::run (dyn.processBand (aud, sc, 2, block, band));
            deepest = std::fmin (deepest, dyn.deltaDb (eq::Lane::Stereo));
        }
        ok (deepest < -6.0, std::string ("precondition ") + cs.name + ": the tone really earned gain reduction");

        // The excursion: mono, on silence, long enough for the delta to release all the way back.
        for (int done = 0; done < (int) kFs; done += block)
        {
            std::fill (L.begin(), L.end(), 0.0f); std::fill (sl.begin(), sl.end(), 0.0f);
            std::fill (R.begin(), R.end(), 0.0f); std::fill (sr.begin(), sr.end(), 0.0f);
            felitronics::test::run (dyn.processBand (aud, sc, 1, block, band));
        }
        const double before = dyn.deltaDb (eq::Lane::Stereo);
        approx (before, 0.0, 0.05, std::string ("precondition ") + cs.name + ": the delta released before the return");

        // The return: stereo, still DIGITAL SILENCE. Nothing may be EARNED from nothing — which is not the
        // same as "the delta is zero". The shared follower is still finishing its release when the channel
        // comes back (measured 2.6e-05 dB on the way down), and that residual is legitimate: it belongs to
        // the excursion, not to the return. The defect showed as the delta getting DEEPER, all the way to
        // -11.97 dB, so that is what is asserted — a bound on the release residual would have measured the
        // follower and called it the fix.
        double deepestAfter = 0.0, audio = 0.0;
        for (int done = 0; done < block * 40; done += block)
        {
            std::fill (L.begin(), L.end(), 0.0f); std::fill (sl.begin(), sl.end(), 0.0f);
            std::fill (R.begin(), R.end(), 0.0f); std::fill (sr.begin(), sr.end(), 0.0f);
            felitronics::test::run (dyn.processBand (aud, sc, 2, block, band));
            deepestAfter = std::fmin (deepestAfter, dyn.deltaDb (eq::Lane::Stereo));
            audio        = std::fmax (audio, std::fmax (peakOf (L), peakOf (R)));
        }
        ok (deepestAfter >= before - 1.0e-4, std::string (cs.name) + ": the return earns no gain reduction at all (was -11.97 dB)");
        ok (audio == 0.0, std::string (cs.name) + ": and puts no audio into digital silence");
    }
}

static void testLaneDynamicsIsolation()
{
    group ("LaneDynamics: the shared detector still follows the channels that are actually present");

    // The linked envelope is SHARED, so it legitimately changes when a channel leaves — this asserts the
    // opposite of an isolation claim, and it is what stops the fix from over-reaching into the shared half.
    const int block = 64;
    const eq::BandParams p = dynPoint (1000.0, 2.0);
    eq::EqBand band; band.prepare (kFs, 2); band.setParams (p);
    dynamiceq::LaneDynamics dyn; felitronics::test::run (dyn.prepare (kFs, 2)); dyn.setParams (p);

    std::vector<float> L ((std::size_t) block), R ((std::size_t) block), sl ((std::size_t) block), sr ((std::size_t) block);
    float* aud[2] { L.data(), R.data() };
    const float* sc[2] { sl.data(), sr.data() };

    long ph = 0;
    for (int done = 0; done < (int) kFs; done += block)     // only ch1 is loud
    {
        for (int i = 0; i < block; ++i, ++ph)
        {
            const float v = (float) std::sin (2.0 * core::kPi * 1000.0 * (double) ph / kFs);
            L[(std::size_t) i] = sl[(std::size_t) i] = 0.0f;
            R[(std::size_t) i] = sr[(std::size_t) i] = v;
        }
        felitronics::test::run (dyn.processBand (aud, sc, 2, block, band));
    }
    const double withBoth = dyn.deltaDb (eq::Lane::Stereo);
    ok (withBoth < -6.0, "precondition: the loud channel is driving the linked detector");

    for (int done = 0; done < (int) (kFs / 2); done += block)   // drop to mono: the loud channel is gone
    {
        for (int i = 0; i < block; ++i, ++ph) { L[(std::size_t) i] = sl[(std::size_t) i] = 0.0f; }
        felitronics::test::run (dyn.processBand (aud, sc, 1, block, band));
    }
    ok (dyn.deltaDb (eq::Lane::Stereo) > withBoth + 5.0,
        "losing the loud channel really releases the shared gain — the shared half was not frozen by the fix");
}

static void testLaneDynamicsStayingColumn()
{
    group ("LaneDynamics: losing a SILENT channel must not disturb the loud one that stayed");

    // The mirror image of the probe test, and the one that pins the clear to a single column. Here the
    // channel that LEAVES is silent and the channel that STAYS carries the tone, so a wholesale
    // probe.reset() on the edge would restart the staying column's band-pass, the linked envelope would
    // collapse for a few milliseconds and the delta would visibly let go of a signal that never stopped.
    // The threshold is deliberately NOT deep: at -60 dBFS the reduction sits pinned against rangeDb, and a
    // pinned delta cannot register a detector dip at all — the first version of this test measured 0.000
    // for the bug and for the fix alike. At -18 dBFS the delta sits in the linear region, where losing the
    // staying channel's probe column moves it 2.26 dB against the 0.04 dB the honest per-column clear costs.
    const int block = 64;
    const eq::BandParams p = dynPoint (100.0, 10.0, -18.0);
    eq::EqBand band; band.prepare (kFs, 2); band.setParams (p);
    dynamiceq::LaneDynamics dyn; felitronics::test::run (dyn.prepare (kFs, 2)); dyn.setParams (p);

    std::vector<float> L ((std::size_t) block), R ((std::size_t) block), sl ((std::size_t) block), sr ((std::size_t) block);
    float* aud[2] { L.data(), R.data() };
    const float* sc[2] { sl.data(), sr.data() };

    long ph = 0;
    auto feed = [&] (int nc)
    {
        for (int i = 0; i < block; ++i, ++ph)
        {
            const float v = (float) std::sin (2.0 * core::kPi * 100.0 * (double) ph / kFs);
            L[(std::size_t) i] = sl[(std::size_t) i] = v;      // channel 0 is the LOUD one, and it stays
            R[(std::size_t) i] = sr[(std::size_t) i] = 0.0f;   // channel 1 is silent, and it leaves
        }
        felitronics::test::run (dyn.processBand (aud, sc, nc, block, band));
    };

    for (int done = 0; done < (int) kFs; done += block) feed (2);
    const double before = dyn.deltaDb (eq::Lane::Stereo);
    ok (before < -2.0 && before > -11.5, "precondition: the delta is in its LINEAR region, not pinned at the range");

    double worst = 0.0;
    for (int k = 0; k < 40; ++k) { feed (1); worst = std::fmax (worst, std::fabs (dyn.deltaDb (eq::Lane::Stereo) - before)); }
    ok (worst < 0.3, "the delta does not lurch when the silent channel leaves — only its own column was cleared");
}

static void testDynamicEqBandChannelGate()
{
    group ("DynamicEqBand: a returning channel brings back neither audio nor unearned gain");

    struct Case { double f, Q; const char* name; };
    const Case cases[2] { { 1000.0, 2.0, "1 kHz / Q2" }, { 100.0, 10.0, "100 Hz / Q10" } };

    for (const Case& cs : cases)
    {
        const int N = 64;
        dynamiceq::DynamicEqBand b;
        dynamiceq::DynamicEqBandParams p;
        p.freq = cs.f; p.Q = cs.Q; p.thresholdDb = -60.0; p.rangeDb = 12.0; p.staticGainDb = 0.0;
        felitronics::test::run (b.prepare (kFs, 2)); b.setParams (p);

        std::vector<float> L ((std::size_t) N, 0.0f), R ((std::size_t) N, 0.0f);
        float* io[2] { L.data(), R.data() };

        long ph = 0;
        double drove = 0.0;
        for (int done = 0; done < (int) kFs; done += N)      // ch1 loud at the band centre, ch0 silent
        {
            for (int i = 0; i < N; ++i, ++ph)
            {
                L[(std::size_t) i] = 0.0f;
                R[(std::size_t) i] = (float) std::sin (2.0 * core::kPi * cs.f * (double) ph / kFs);
            }
            felitronics::test::run (b.process (io, 2, N));
            drove = std::fmax (drove, std::fabs (b.dynamicDeltaDb()));
        }
        ok (drove > 3.0, std::string ("precondition ") + cs.name + ": the tone really moved the band");

        for (int done = 0; done < (int) kFs; done += N)      // the excursion: mono, silence
        {
            std::fill (L.begin(), L.end(), 0.0f);
            felitronics::test::run (b.process (io, 1, N));
        }

        std::fill (L.begin(), L.end(), 0.0f); std::fill (R.begin(), R.end(), 0.0f);
        double worst = 0.0;
        for (int k = 0; k < 20; ++k)
        {
            std::fill (L.begin(), L.end(), 0.0f); std::fill (R.begin(), R.end(), 0.0f);
            felitronics::test::run (b.process (io, 2, N));
            worst = std::fmax (worst, std::fmax (peakOf (L), peakOf (R)));
        }
        ok (worst == 0.0, std::string (cs.name) + ": exact zero out of silence on the returned channel (was 0.388)");
        approx (b.dynamicDeltaDb(), 0.0, 0.05, std::string (cs.name) + ": and no unearned gain behind it (was 8.85 dB)");
    }
}

static void testDynamicEqBandIsolation()
{
    group ("DynamicEqBand: the channel that stayed is bit-exact; the one that returned is not owed that");

    const int N = 64;
    dynamiceq::DynamicEqBand dut, ref;
    dynamiceq::DynamicEqBandParams p;
    p.freq = 800.0; p.Q = 1.5; p.thresholdDb = -70.0; p.rangeDb = 9.0;
    felitronics::test::run (dut.prepare (kFs, 2)); felitronics::test::run (ref.prepare (kFs, 2));
    dut.setParams (p);   ref.setParams (p);

    std::vector<float> d0 ((std::size_t) N), d1 ((std::size_t) N), r0 ((std::size_t) N), r1 ((std::size_t) N);
    float* dio[2] { d0.data(), d1.data() };
    float* rio[2] { r0.data(), r1.data() };

    // Channel 0 identical on both, channel 1 identical too — the DUT simply stops being told about it for
    // a while. The detector is shared, so ch0's audio may legitimately move; what may NOT happen is a
    // divergence AFTER the return that outlives the shared detector's own recovery. Compare the settled
    // tail rather than the transient, which is the honest form of this claim for a linked detector.
    long ph = 0;
    for (int k = 0; k < 400; ++k)
    {
        for (int i = 0; i < N; ++i, ++ph)
        {
            const float v = (float) (0.4 * std::sin (2.0 * core::kPi * 800.0 * (double) ph / kFs));
            d0[(std::size_t) i] = r0[(std::size_t) i] = v;
            d1[(std::size_t) i] = r1[(std::size_t) i] = v;
        }
        felitronics::test::run (dut.process (dio, (k >= 100 && k < 200) ? 1 : 2, N));
        felitronics::test::run (ref.process (rio, 2, N));
    }
    ok (peakOf (d0) > 0.05, "precondition: the stream is still carrying signal at the end");

    // Channel 0 never left: it is an S cell, and S is where bit-equality is owed and delivered. Note this
    // also proves the SHARED detector reconverged exactly — had it not, ch0's coefficients would differ.
    ok (bitEqual (d0, r0), "the channel that never left is bit-equal to the never-left run");

    // Channel 1 returned: it is an R cell, and R is owed "equal to a fresh column at the return", never
    // "equal to a column that kept running". Its restarted recursion converges back toward the reference
    // and settles within a couple of float ulp of it — measured 2.98e-08 against a 0.14 signal, i.e. two
    // ulp at that magnitude, which is where two recursions fed identical input from different states end
    // up and stay. Demanding bit-equality here would be demanding that the gap never happened.
    double worst = 0.0;
    for (std::size_t i = 0; i < d1.size(); ++i) worst = std::fmax (worst, std::fabs ((double) d1[i] - (double) r1[i]));
    ok (worst < 1.0e-6, "the channel that returned reconverges to within a few ulp of it");
}

//==================================================================================================
// What a PARK does to the programme estimate — the policy, asserted by its contract.
//==================================================================================================
// A parked lane cannot see the programme move. Keeping its picture and discarding it are both wrong in
// opposite cases, and from the inside the lane cannot tell them apart: on return it only sees "a signal
// N dB from my last known norm". The one thing it does know is HOW LONG it was blind, and that is what
// the policy keys on, with thresholds taken from the estimator's own averaging constant.
//
// This replaces a test that could not see the question at all: it inserted the park only into the parked
// run, so both rigs had heard exactly the same programme and "parked == stayed" was true by construction.
// Its recovery check was also unsound — it latched the FIRST near-zero crossing, so a trajectory that
// never recovered still reported a settle time. Both are fixed here: the reference HEARS the change, and
// recovery is the LAST crossing plus a requirement that it stays recovered.
static void testParkPolicy()
{
    group ("auto threshold: how long the lane was blind decides what happens to its programme estimate");

    const int B = 64;
    const double f = 1000.0;

    // Returns the largest |DUT - REF| over the window after the return, and the last time it exceeded
    // 0.5 dB. `heard` selects the reference: one that ran through the stretch, or one that also missed it.
    struct Trace { double worst, lastAboveMs, at10, at100; };
    // `viaDyn` selects HOW the lane is parked: by the channel count, or by switching dynamics off — a
    // different code path (processBand's disengaged branch) that reaches the SAME policy. A version of
    // this test that only ever parked by channel count left that path unmeasured.
    auto compare = [&] (double gapSec, bool refHeard, bool viaDyn = false) -> Trace
    {
        auto build = [&] (dynamiceq::LaneDynamics& dyn, eq::EqBand& band)
        {
            eq::BandParams p;
            p.on = true; p.type = eq::FilterType::Bell;
            p.lane (eq::Lane::Stereo).on = false;                 // Side only: a mono stretch parks it
            eq::LaneParams& sd = p.lane (eq::Lane::Side);
            sd.on = true; sd.freq = f; sd.Q = 2.0; sd.gainDb = 0.0;
            p.dyn.on = true; p.dyn.rangeDb = -24.0;               // inside the rail, so nothing pins
            p.dyn.thrAuto = true;
            band.prepare (kFs, 2); band.setParams (p);
            felitronics::test::run (dyn.prepare (kFs, 2));  dyn.setParams (p);
        };
        eq::EqBand bd, br; dynamiceq::LaneDynamics dd, dr;
        build (dd, bd); build (dr, br);
        auto setDyn = [&] (dynamiceq::LaneDynamics& dyn, eq::EqBand& band, bool on)
        {
            eq::BandParams q = band.params(); q.dyn.on = on; band.setParams (q); dyn.setParams (q);
        };

        std::vector<float> L ((std::size_t) B), R ((std::size_t) B), sl ((std::size_t) B), sr ((std::size_t) B);
        float* aud[2] { L.data(), R.data() };
        const float* sc[2] { sl.data(), sr.data() };
        long ph = 0;
        auto feed = [&] (dynamiceq::LaneDynamics& dyn, eq::EqBand& band, double db, int nc, int blocks)
        {
            const double amp = core::dbToGain (db);
            for (int k = 0; k < blocks; ++k)
            {
                for (int i = 0; i < B; ++i)
                {
                    const float v = (float) (amp * std::sin (2.0 * core::kPi * f * (double) (ph + i) / kFs));
                    L[(std::size_t) i] = sl[(std::size_t) i] =  v;
                    R[(std::size_t) i] = sr[(std::size_t) i] = -v;
                }
                felitronics::test::run (dyn.processBand (aud, sc, nc, B, band));
                ph += B;
            }
        };
        // Both learn the same quiet programme, then it becomes loud DURING the stretch. The DUT is parked
        // through it; the reference either hears it (nc = 2) or is parked too.
        const long save = ph;
        feed (dd, bd, -40.0, 2, (int) (kFs * 4.0 / B));
        ph = save; feed (dr, br, -40.0, 2, (int) (kFs * 4.0 / B));
        const long afterLearn = ph;
        if (viaDyn) setDyn (dd, bd, false);                                      // park by switching dynamics off
        feed (dd, bd, -10.0, viaDyn ? 2 : 1, (int) (kFs * gapSec / B));           // DUT: blind either way
        if (viaDyn) setDyn (dd, bd, true);
        ph = afterLearn; feed (dr, br, -10.0, refHeard ? 2 : 1, (int) (kFs * gapSec / B));
        const long afterGap = ph;

        Trace t { 0.0, -1.0, 0.0, 0.0 };
        const int n = (int) (kFs * 6.0 / B);
        for (int k = 0; k < n; ++k)
        {
            const long here = ph;
            feed (dd, bd, -10.0, 2, 1);
            ph = here; feed (dr, br, -10.0, 2, 1);
            const double d = std::fabs (dd.deltaDb (eq::Lane::Side) - dr.deltaDb (eq::Lane::Side));
            const double ms = (double) k * B / kFs * 1000.0;
            t.worst = std::fmax (t.worst, d);
            if (d > 0.5) t.lastAboveMs = ms;                       // LAST crossing, not the first
            if (ms <= 10.0)  t.at10  = d;
            if (ms <= 100.0) t.at100 = d;
        }
        (void) afterGap;
        return t;
    };

    // PRECONDITION. The two references must actually disagree, or every number below is about nothing.
    const Trace blind = compare (10.0, false);
    approx (blind.worst, 0.0, 1.0, "precondition: against a reference that was ALSO blind, any policy looks fine");

    // SHORT park — below a quarter of the estimator's averaging constant. Too short for the programme to
    // have moved far, so the picture is kept and the lane behaves as it always did.
    const Trace shortPark = compare (0.3, true);
    ok (shortPark.worst > 0.0, "precondition: the short-park fixture is live");
    ok (shortPark.lastAboveMs < 4000.0, "a short park does not leave the lane wrong for the whole settle");

    // LONG park — beyond four times the constant. The programme certainly moved; the returning lane must
    // agree with the one that heard it. This is the case a KEEP policy fails by 17.87 dB for 3.36 s.
    const Trace longPark = compare (10.0, true);
    approx (longPark.worst, 0.0, 1.0,
            "after a long park the lane agrees with the reference that heard the change (keep: 17.87 dB)");
    ok (longPark.lastAboveMs < 0.0 || longPark.lastAboveMs < 200.0,
        "and it agrees from the start, not after a settle");

    // MID park — between the two. The picture is kept but the estimator is told to correct fast, so the
    // divergence must SHRINK: a policy that merely keeps holds it flat for seconds.
    const Trace midPark = compare (2.0, true);
    ok (midPark.at10 > 1.0, "precondition: the mid-park fixture really starts out wrong");
    ok (midPark.at100 < 0.75 * midPark.at10, "a mid-length park corrects fast rather than holding its error");
    ok (midPark.lastAboveMs > 0.0 && midPark.lastAboveMs < 1500.0,
        "and it is done inside the fast window, not after the full averaging constant");

    // THE OTHER PARK PATH. Switching dynamics off leaves through processBand's disengaged branch, never
    // reaching the per-lane code the cases above exercise. It must accrue blind time the same way, or a
    // plain A/B bypass across a programme change comes back with the full unearned reduction.
    const Trace longDyn = compare (10.0, true, true);
    approx (longDyn.worst, 0.0, 1.0, "a long park via dyn.on obeys the same policy as one via the channel count");
    const Trace shortDyn = compare (0.3, true, true);
    ok (shortDyn.lastAboveMs < 4000.0, "and a short one via dyn.on does too");
}

//==================================================================================================
int main()
{
    std::printf ("felitronics::dynamiceq — execution-gate state\n");

    testLaneDynamicsProbe();
    testLaneDynamicsIsolation();
    testLaneDynamicsStayingColumn();
    testDynamicEqBandChannelGate();
    testDynamicEqBandIsolation();
    testParkPolicy();

    group ("RT-safety");
    {
        const int N = 64;
        const eq::BandParams p = dynPoint (1000.0, 2.0);
        eq::EqBand band; band.prepare (kFs, 4); band.setParams (p);
        dynamiceq::LaneDynamics dyn; felitronics::test::run (dyn.prepare (kFs, 4)); dyn.setParams (p);
        std::vector<float> v[4], s[4];
        float* aud[4] {}; const float* sc[4] {};
        for (int c = 0; c < 4; ++c)
        {
            v[c].assign ((std::size_t) N, 0.2f); s[c].assign ((std::size_t) N, 0.2f);
            aud[c] = v[c].data(); sc[c] = s[c].data();
        }
        felitronics::test::run (dyn.processBand (aud, sc, 2, N, band));
        const int before = g_allocs.load();
        for (int k = 0; k < 40; ++k) felitronics::test::run (dyn.processBand (aud, sc, (k % 3) + 1, N, band));
        felitronics::test::okNoAlloc (g_allocs.load() == before, "no allocation across 40 blocks of changing width");
    }

    return felitronics::test::report();
}

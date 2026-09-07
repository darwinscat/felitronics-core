// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// What a cell of this band does when it STOPS being executed and is executed again later.
//
// A filter that is not called does not decay — it freezes. Every gate below can stop a real recursion
// while the band keeps running around it: the channel count (L/R/M/S are stereo-only, the ST columns are
// per channel), a lane switched off, the swept/matched branch, and dyn.on. Before this suite a Side lane
// that sat out a mono stretch came back and replayed 2.547 (+8.12 dBFS) into DIGITAL SILENCE, and it did
// so on the very first sample.
//
// The four acceptance invariants are scoped the only way they can both hold: partition the state at a
// transition into S (cells that kept executing) and R (cells that stopped). PREVENTION is a statement
// about S — a surviving column is bit-equal to the run where nothing ever left. RECOVERY is a statement
// about R — at the return a stopped cell equals a fresh one, in state AND in the coefficients it will
// filter with. Asking one cell for both at once is unsatisfiable, and asking it was the mistake this
// scoping replaces.
//
// Every fixture asserts its own PRECONDITION first. A test that silently measures a lane which never
// engaged reports success about nothing, which is how two checks in an earlier task passed while the
// band under them was switched off.

#include <felitronics_test.h>
#include <felitronics/eq/EqBand.h>
#include <felitronics/eq/EqEngine.h>
#include <felitronics/eq/Svf.h>

#include <atomic>
#include <cmath>
#include <complex>
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
using namespace felitronics::eq;
using felitronics::test::ok;
using felitronics::test::approx;
using felitronics::test::group;

static constexpr double kFs = 48000.0;

// ST lane active and FLAT — the companion that keeps `anyRun` true and hides the stopped lane. That is
// the whole point of the fixture: with the ST lane off the band self-cleans and the defect never appears.
static BandParams stPlusSide (double sideHz, double sideGainDb, double sideQ = 1.0)
{
    BandParams p;
    p.on = true; p.type = FilterType::Bell;
    LaneParams& st = p.lane (Lane::Stereo);
    st.on = true; st.freq = 1000.0; st.Q = 1.0; st.gainDb = 0.0;
    LaneParams& sd = p.lane (Lane::Side);
    sd.on = true; sd.freq = sideHz; sd.Q = sideQ; sd.gainDb = sideGainDb;
    return p;
}

static BandParams stOnly (double hz, double Q, double gainDb)
{
    BandParams p;
    p.on = true; p.type = FilterType::Bell;
    LaneParams& st = p.lane (Lane::Stereo);
    st.on = true; st.freq = hz; st.Q = Q; st.gainDb = gainDb;
    return p;
}

static double peakOf (const std::vector<float>& v)
{
    double m = 0.0; for (float x : v) m = std::fmax (m, (double) std::fabs (x)); return m;
}

static bool bitEqual (const std::vector<float>& a, const std::vector<float>& b)
{
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i) if (! (a[i] == b[i])) return false;   // exact, not approx
    return true;
}

static void fillSine (std::vector<float>& v, double hz, double amp, double& phase)
{
    const double w = 2.0 * core::kPi * hz / kFs;
    for (auto& x : v) { x = (float) (amp * std::sin (phase)); phase += w; }
}

static void fillNoise (std::vector<float>& v, unsigned seed)
{
    unsigned s = seed;
    for (auto& x : v) { s = s * 1664525u + 1013904223u; x = (float) ((int) (s >> 8) % 2000 - 1000) * 0.0005f; }
}

//==================================================================================================
// 1. The defect itself, on both channel-count paths. Exact zero, not "small".
//==================================================================================================
static void testChannelGateSilence()
{
    group ("channel gate: a returning lane emits nothing from digital silence");

    for (int mode = 0; mode < 2; ++mode)   // 0: 2->1->2   1: 2->3->2 with four channels prepared
    {
        const int prepared = (mode == 0) ? 2 : 4;
        const int gapNc    = (mode == 0) ? 1 : 3;
        const char* name   = (mode == 0) ? "2->1->2" : "2->3->2";
        const int N = 1024;

        EqBand b; b.prepare (kFs, prepared);
        b.setParams (stPlusSide (200.0, 12.0));

        // PRECONDITION: the Side lane is genuinely engaged at the fixture's frequency.
        const double w200 = 2.0 * core::kPi * 200.0 / kFs;
        std::vector<float> bufs[4];
        float* ch[4] {};
        for (int c = 0; c < prepared; ++c) { bufs[c].assign ((std::size_t) N, 0.0f); ch[c] = bufs[c].data(); }
        felitronics::test::run (b.processBlock (ch, 2, N));                                    // one block so the design is applied
        approx (core::gainToDb (std::abs (b.response (w200, Axis::Side))), 12.0, 0.5,
                std::string ("precondition ") + name + ": Side lane really is +12 dB at 200 Hz");

        // Charge: antiphase 200 Hz is pure Side, so the lane under test is the one that gets the energy.
        double ph = 0.0;
        for (int k = 0; k < 20; ++k)
        {
            fillSine (bufs[0], 200.0, 1.0, ph);
            for (int i = 0; i < N; ++i) bufs[1][(std::size_t) i] = -bufs[0][(std::size_t) i];
            for (int c = 2; c < prepared; ++c) std::fill (bufs[c].begin(), bufs[c].end(), 0.0f);
            felitronics::test::run (b.processBlock (ch, 2, N));
        }
        ok (peakOf (bufs[0]) > 1.5, std::string ("precondition ") + name + ": the charge really charged the lane");

        // The gap: the Side lane is not executed at all (nc != 2), while the ST lane keeps the band alive.
        for (int k = 0; k < 48; ++k)
        {
            for (int c = 0; c < prepared; ++c) std::fill (bufs[c].begin(), bufs[c].end(), 0.0f);
            felitronics::test::run (b.processBlock (ch, gapNc, N));
        }

        // Return, on DIGITAL SILENCE.
        for (int c = 0; c < prepared; ++c) std::fill (bufs[c].begin(), bufs[c].end(), 0.0f);
        felitronics::test::run (b.processBlock (ch, 2, N));
        ok (peakOf (bufs[0]) == 0.0 && peakOf (bufs[1]) == 0.0,
            std::string (name) + ": silence in, exact zero out (was 2.547 = +8.12 dBFS at i=0)");
    }
}

//==================================================================================================
// 2. PREVENTION, on S: a channel that never left is bit-equal to the run where nothing left.
//    ST-only by construction — with an M/S lane running, ch0's output legitimately depends on ch1.
//==================================================================================================
static void testPreventionOnSurvivingChannel()
{
    group ("prevention: the surviving channel is bit-exact against the never-left run");

    for (int mode = 0; mode < 2; ++mode)
    {
        const int prepared = (mode == 0) ? 2 : 4;
        const int gapNc    = (mode == 0) ? 1 : 3;
        const char* name   = (mode == 0) ? "2->1->2" : "2->3->2";
        const int N = 256, charge = 12, gap = 30, tail = 12;

        EqBand dut, ref;
        dut.prepare (kFs, prepared); ref.prepare (kFs, prepared);
        BandParams p = stOnly (900.0, 1.2, 6.0);
        p.dyn.on = true;                                   // the delta bell is per channel too — include it
        dut.setParams (p); ref.setParams (p);
        dut.setLaneDeltaDb (Lane::Stereo, -7.0);
        ref.setLaneDeltaDb (Lane::Stereo, -7.0);

        std::vector<float> d[4], r[4];
        float* dch[4] {}; float* rch[4] {};
        for (int c = 0; c < prepared; ++c)
        {
            d[c].assign ((std::size_t) N, 0.0f); r[c].assign ((std::size_t) N, 0.0f);
            dch[c] = d[c].data(); rch[c] = r[c].data();
        }

        bool allEqual = true;
        double energy = 0.0;
        for (int k = 0; k < charge + gap + tail; ++k)
        {
            const bool inGap = (k >= charge && k < charge + gap);
            for (int c = 0; c < prepared; ++c)
            {
                fillNoise (d[c], (unsigned) (1234u + 77u * (unsigned) k + 5u * (unsigned) c));
                r[c] = d[c];                               // identical input on every channel, every block
            }
            felitronics::test::run (dut.processBlock (dch, inGap ? gapNc : 2, N));
            felitronics::test::run (ref.processBlock (rch, 2, N));                  // the reference never changes its channel count
            const int surviving = (mode == 0) ? 1 : 2;     // channels that are present in EVERY segment
            for (int c = 0; c < surviving; ++c)
            {
                allEqual = allEqual && bitEqual (d[c], r[c]);
                energy += peakOf (d[c]);
            }
        }
        ok (energy > 1.0, std::string ("precondition ") + name + ": the surviving channel carried signal");
        ok (allEqual, std::string (name) + ": every surviving channel bit-equal to the never-left run");
    }
}

//==================================================================================================
// 3. RECOVERY, on R. Three assertions, none of which compares R against a never-left stream.
//==================================================================================================
static void testRecoveryCoefficients()
{
    group ("recovery: a lane returns with CURRENT coefficients, not the ones in force when it left");

    // A lane whose smoothers settle while it is parked used to come back designed at the value of the
    // first block of the gap: measured 798 Hz for a 500 -> 4000 Hz edit, +0.64 dB where +12 was asked.
    // Straddle the self-healing boundary (tau*ln(|delta|/1e-7) ~ 0.73 s here): the two short stretches
    // recovered even before the fix, the two long ones did not.
    const int N = 64;
    const int stretches[4] { 8, 75, 750, 1500 };
    const double w4k = 2.0 * core::kPi * 4000.0 / kFs;

    for (int si = 0; si < 4; ++si)
    {
        EqBand dut, ref;
        dut.prepare (kFs, 2); ref.prepare (kFs, 2);
        BandParams p = stPlusSide (500.0, 12.0, 2.0);
        dut.setParams (p); ref.setParams (p);

        std::vector<float> d0 (N, 0.0f), d1 (N, 0.0f), r0 (N, 0.0f), r1 (N, 0.0f);
        float* dch[2] { d0.data(), d1.data() };
        float* rch[2] { r0.data(), r1.data() };
        for (int k = 0; k < 20; ++k) { felitronics::test::run (dut.processBlock (dch, 2, N)); felitronics::test::run (ref.processBlock (rch, 2, N)); }

        p.lane (Lane::Side).freq = 4000.0;                 // the edit, issued one block before the gap
        dut.setParams (p); ref.setParams (p);
        felitronics::test::run (dut.processBlock (dch, 2, N)); felitronics::test::run (ref.processBlock (rch, 2, N));

        for (int k = 0; k < stretches[si]; ++k) { felitronics::test::run (dut.processBlock (dch, 1, N)); felitronics::test::run (ref.processBlock (rch, 2, N)); }
        felitronics::test::run (dut.processBlock (dch, 2, N)); felitronics::test::run (ref.processBlock (rch, 2, N));

        const std::complex<double> hd = dut.response (w4k, Axis::Side);
        const std::complex<double> hr = ref.response (w4k, Axis::Side);
        char msg[160];
        std::snprintf (msg, sizeof msg, "stretch %d blocks: returned coefficients bit-equal to never-left", stretches[si]);
        ok (hd.real() == hr.real() && hd.imag() == hr.imag(), msg);

        if (si == 3)   // PRECONDITION on the longest case: the reference really did reach the new design
            approx (core::gainToDb (std::abs (hr)), 12.0, 0.2,
                    "precondition: the never-left reference really settled at +12 dB / 4 kHz");
    }
}

static void testRecoveryStrongNull()
{
    group ("recovery: above t_forget the returned lane nulls a never-left one that heard silence");

    // The one regime where "never left" and "reset at return" coincide bit-for-bit: let the reference's
    // own tail decay to EXACT zero through the per-block flush, and the two definitions agree, so the
    // comparison can be sample-by-sample on real signal rather than on silence.
    const int N = 256;
    EqBand dut, ref;
    dut.prepare (kFs, 2); ref.prepare (kFs, 2);
    const BandParams p = stPlusSide (1000.0, 12.0, 1.0);   // tau = Q/(pi*f) ~ 0.3 ms: forgets in a few hundred samples
    dut.setParams (p); ref.setParams (p);

    std::vector<float> d0 (N), d1 (N), r0 (N), r1 (N);
    float* dch[2] { d0.data(), d1.data() };
    float* rch[2] { r0.data(), r1.data() };

    double ph = 0.0, charged = 0.0;
    for (int k = 0; k < 12; ++k)                            // charge both identically, in stereo
    {
        fillSine (d0, 1000.0, 1.0, ph);
        for (int i = 0; i < N; ++i) d1[(std::size_t) i] = -d0[(std::size_t) i];
        r0 = d0; r1 = d1;
        felitronics::test::run (dut.processBlock (dch, 2, N)); felitronics::test::run (ref.processBlock (rch, 2, N));
        charged = peakOf (d0);
    }
    ok (charged > 1.5, "precondition: both instances really charged the Side lane");

    for (int k = 0; k < 40; ++k)                            // the gap: DUT parks the lane, REF hears silence
    {
        std::fill (d0.begin(), d0.end(), 0.0f); std::fill (d1.begin(), d1.end(), 0.0f);
        std::fill (r0.begin(), r0.end(), 0.0f); std::fill (r1.begin(), r1.end(), 0.0f);
        felitronics::test::run (dut.processBlock (dch, 1, N));
        felitronics::test::run (ref.processBlock (rch, 2, N));
    }
    ok (peakOf (r0) == 0.0 && peakOf (r1) == 0.0, "precondition: the reference's own tail reached EXACT zero");

    bool equal = true; double e = 0.0;
    for (int k = 0; k < 10; ++k)                            // return: identical real signal into both
    {
        fillNoise (d0, 900u + 3u * (unsigned) k); fillNoise (d1, 4300u + 7u * (unsigned) k);
        r0 = d0; r1 = d1;
        felitronics::test::run (dut.processBlock (dch, 2, N)); felitronics::test::run (ref.processBlock (rch, 2, N));
        equal = equal && bitEqual (d0, r0) && bitEqual (d1, r1);
        e += peakOf (d0);
    }
    ok (e > 1.0, "precondition: the post-return segment carried signal");
    ok (equal, "post-return output is bit-equal, sample by sample, to the silence-fed never-left run");
}

//==================================================================================================
// 4. The gates that are NOT the channel count. Each one stops a real recursion at a constant nc.
//==================================================================================================
static void testNonChannelGates()
{
    group ("the other execution gates: swept/matched, dyn.on, a lane switched off");

    const int N = 256;

    // (a) matched -> swept -> matched. The swept branch never touches deltaST_, so its integrators used
    //     to sit out the excursion and resume: measured 0.690 (-3.2 dBFS) out of silence at 100 Hz / Q10.
    {
        EqBand b; b.prepare (kFs, 2);
        BandParams p = stOnly (100.0, 10.0, 0.0);
        p.dyn.on = true;
        b.setParams (p);
        b.setLaneDeltaDb (Lane::Stereo, -12.0);

        std::vector<float> a (N), c (N);
        float* ch[2] { a.data(), c.data() };
        double ph = 0.0, charged = 0.0;
        for (int k = 0; k < 40; ++k) { fillSine (a, 100.0, 1.0, ph); c = a; felitronics::test::run (b.processBlock (ch, 2, N)); charged = peakOf (a); }
        ok (charged > 0.1, "precondition: the ST delta bell really carried signal");

        p.swept = true; b.setParams (p);                    // ST-only point, so this really is the swept path
        for (int k = 0; k < 40; ++k) { std::fill (a.begin(), a.end(), 0.0f); c = a; felitronics::test::run (b.processBlock (ch, 2, N)); }
        p.swept = false; b.setParams (p);
        std::fill (a.begin(), a.end(), 0.0f); std::fill (c.begin(), c.end(), 0.0f);
        felitronics::test::run (b.processBlock (ch, 2, N));
        ok (peakOf (a) == 0.0 && peakOf (c) == 0.0, "swept -> matched: exact zero (was 0.690 = -3.2 dBFS)");
    }

    // (b) dyn.on off and on again. Every delta SVF stops; nothing in setParams treats that as topology.
    //     "Exact zero out of silence" would be the WRONG oracle here and the fixture says why: a matched
    //     Bell at 0 dB is not bit-exactly an identity, so the static column keeps a real, continuous tail
    //     of its own (measured 1.4e-10, decaying) that has nothing to do with the gate. The delta bell is
    //     isolated instead — against a reference that shares the static history exactly (the delta sits
    //     AFTER the static sections, so running it or not cannot change them) and whose delta bell never
    //     held state. If the stopped delta resumes, the two diverge; if it restarts clean, they cannot.
    {
        EqBand dut, ref;
        dut.prepare (kFs, 2); ref.prepare (kFs, 2);
        BandParams p = stOnly (100.0, 10.0, 0.0);
        BandParams q = p;
        p.dyn.on = true;                                    // DUT: the delta bell runs and gets charged
        q.dyn.on = false;                                   // REF: identical static path, delta never runs
        dut.setParams (p); ref.setParams (q);
        dut.setLaneDeltaDb (Lane::Stereo, -12.0);

        std::vector<float> d0 (N), d1 (N), r0 (N), r1 (N);
        float* dch[2] { d0.data(), d1.data() };
        float* rch[2] { r0.data(), r1.data() };
        double ph = 0.0, deltaSeen = 0.0;
        for (int k = 0; k < 40; ++k)
        {
            fillSine (d0, 100.0, 1.0, ph); d1 = d0; r0 = d0; r1 = d0;
            felitronics::test::run (dut.processBlock (dch, 2, N)); felitronics::test::run (ref.processBlock (rch, 2, N));
            deltaSeen = std::fabs (peakOf (d0) - peakOf (r0));
        }
        ok (deltaSeen > 0.1, "precondition: the delta bell really was in the DUT's path and not the REF's");

        p.dyn.on = false; dut.setParams (p);                // the gate closes: the delta cells stop
        for (int k = 0; k < 4; ++k)
        {
            std::fill (d0.begin(), d0.end(), 0.0f); d1 = d0; r0 = d0; r1 = d0;
            felitronics::test::run (dut.processBlock (dch, 2, N)); felitronics::test::run (ref.processBlock (rch, 2, N));
        }
        p.dyn.on = true; dut.setParams (p);                 // and opens again, on a fresh command
        q.dyn.on = true; ref.setParams (q);
        dut.setLaneDeltaDb (Lane::Stereo, -12.0);
        ref.setLaneDeltaDb (Lane::Stereo, -12.0);

        bool equal = true;
        for (int k = 0; k < 4; ++k)
        {
            std::fill (d0.begin(), d0.end(), 0.0f); d1 = d0; r0 = d0; r1 = d0;
            felitronics::test::run (dut.processBlock (dch, 2, N)); felitronics::test::run (ref.processBlock (rch, 2, N));
            equal = equal && bitEqual (d0, r0) && bitEqual (d1, r1);
        }
        ok (equal, "dyn.on off -> on: the delta bell resumes from a clean state, not from before the gap");
    }

    // (c) a lane switched off and on. Its static columns were already cleared by the topology path; its
    //     DELTA bell was not, and that is the half this closes.
    {
        EqBand b; b.prepare (kFs, 2);
        BandParams p = stPlusSide (100.0, 0.0, 10.0);
        p.dyn.on = true; b.setParams (p);
        b.setLaneDeltaDb (Lane::Side, -12.0);

        std::vector<float> a (N), c (N);
        float* ch[2] { a.data(), c.data() };
        double ph = 0.0, charged = 0.0;
        for (int k = 0; k < 40; ++k)
        {
            fillSine (a, 100.0, 1.0, ph);
            for (int i = 0; i < N; ++i) c[(std::size_t) i] = -a[(std::size_t) i];
            felitronics::test::run (b.processBlock (ch, 2, N)); charged = peakOf (a);
        }
        ok (charged > 0.1, "precondition: the Side delta bell really carried signal");

        p.lane (Lane::Side).on = false; b.setParams (p);
        for (int k = 0; k < 4; ++k) { std::fill (a.begin(), a.end(), 0.0f); c = a; felitronics::test::run (b.processBlock (ch, 2, N)); }
        p.lane (Lane::Side).on = true; b.setParams (p);
        b.setLaneDeltaDb (Lane::Side, -12.0);
        std::fill (a.begin(), a.end(), 0.0f); std::fill (c.begin(), c.end(), 0.0f);
        felitronics::test::run (b.processBlock (ch, 2, N));
        ok (peakOf (a) == 0.0 && peakOf (c) == 0.0, "lane off -> on: exact zero out of silence");
    }
}

//==================================================================================================
// 4b. A parameter ramp runs on the caller's clock, and an INERT lane must not decide when it lands.
//==================================================================================================
static void testRampIsWallClock()
{
    group ("a ramp keeps running while the band is idle, with or without a companion lane");

    // The smoothers have always advanced for lanes that were not running; the fully-idle band was the one
    // case that fell out of that rule, because the early return sat above them. The visible cost was that
    // the SAME edit landed at two different times depending on whether an unrelated lane happened to be on
    // — with a flat 0 dB ST companion the design tracked through the gap, without one it froze and finished
    // about 200 ms AFTER the stream came back (measured: parked at 651.5 Hz, then 1548 / 3536 / 3984 Hz at
    // 10 / 50 / 200 ms). An inert lane deciding another lane's behaviour is the shape this file closed once
    // already for state; this is the same shape for design.
    const int N = 64;
    const double target = 4000.0;

    auto peakHz = [] (EqBand& b)
    {
        double best = 0.0, bf = 0.0;
        for (double f = 50.0; f < 12000.0; f *= 1.01)
        {
            const double m = std::abs (b.response (2.0 * core::kPi * f / kFs, Axis::Side));
            if (m > best) { best = m; bf = f; }
        }
        return bf;
    };

    double parked[2] {}, onReturn[2] {};
    for (int companion = 0; companion < 2; ++companion)
    {
        BandParams p;
        p.on = true; p.type = FilterType::Bell;
        p.lane (Lane::Stereo).on = (companion != 0);          // the inert companion, flat at 0 dB
        p.lane (Lane::Stereo).freq = 1000.0; p.lane (Lane::Stereo).gainDb = 0.0;
        LaneParams& sd = p.lane (Lane::Side);
        sd.on = true; sd.freq = 500.0; sd.Q = 2.0; sd.gainDb = 12.0;

        EqBand b; b.prepare (kFs, 2); b.setParams (p);
        std::vector<float> L ((std::size_t) N, 0.01f), R ((std::size_t) N, -0.01f);
        float* ch[2] { L.data(), R.data() };
        auto run = [&] (int nc, int blocks) { for (int k = 0; k < blocks; ++k) felitronics::test::run (b.processBlock (ch, nc, N)); };

        run (2, 20);
        sd.freq = target; p.lane (Lane::Side) = sd; b.setParams (p);   // the edit, one block before the gap
        run (2, 1);
        run (1, (int) (kFs * 1.0 / N));                                 // a second of mono: no Side audio
        parked[companion] = peakHz (b);
        run (2, 1);
        onReturn[companion] = peakHz (b);
    }

    ok (std::fabs (parked[1] - target) < 100.0, "precondition: WITH a companion the design tracked through the gap");
    ok (std::fabs (parked[0] - target) < 100.0, "without one it tracks too — the ramp did not stop with the audio");
    ok (std::fabs (onReturn[0] - onReturn[1]) < 1.0,
        "the edit lands at the same time either way — an inert lane decides nothing");
}

//==================================================================================================
// 5. What must NOT be cleared, and what must NOT count as an edge.
//==================================================================================================
static void testWhatSurvives()
{
    group ("scope: control survives, and a call with no samples is not a transition");

    // (a) The seam is a COMMAND, not memory of past audio. Going idle clears the band's signal history;
    //     it must not silently cancel the gain reduction the producer last asked for.
    {
        EqBand b; b.prepare (kFs, 2);
        BandParams p = stOnly (1000.0, 1.0, 0.0);
        p.dyn.on = true; b.setParams (p);
        b.setLaneDeltaDb (Lane::Stereo, -9.0);

        std::vector<float> a (256), c (256);
        float* ch[2] { a.data(), c.data() };
        double ph = 0.0;
        for (int k = 0; k < 4; ++k) { fillSine (a, 1000.0, 0.5, ph); c = a; felitronics::test::run (b.processBlock (ch, 2, 256)); }

        p.on = false; b.setParams (p);                       // the whole band goes idle: every cell stops
        std::fill (a.begin(), a.end(), 0.0f); c = a;
        felitronics::test::run (b.processBlock (ch, 2, 256));
        approx (b.laneDeltaDb (Lane::Stereo), -9.0, 0.0, "a whole-band idle does NOT cancel the commanded delta");
    }

    // (b) A call carrying no samples ran nothing, so it stopped nothing. Otherwise a host that probes the
    //     graph with an empty block at a different width would silently wipe the stream it is about to run.
    {
        EqBand dut, ref;
        dut.prepare (kFs, 2); ref.prepare (kFs, 2);
        const BandParams p = stPlusSide (200.0, 12.0);
        dut.setParams (p); ref.setParams (p);

        std::vector<float> d0 (128), d1 (128), r0 (128), r1 (128);
        float* dch[2] { d0.data(), d1.data() };
        float* rch[2] { r0.data(), r1.data() };
        bool equal = true; double e = 0.0;
        for (int k = 0; k < 16; ++k)
        {
            fillNoise (d0, 11u + (unsigned) k); fillNoise (d1, 91u + (unsigned) k);
            r0 = d0; r1 = d1;
            if (k == 8) felitronics::test::run (dut.processBlock (dch, 1, 0));        // the empty probe, at a narrower width
            felitronics::test::run (dut.processBlock (dch, 2, 128));
            felitronics::test::run (ref.processBlock (rch, 2, 128));
            equal = equal && bitEqual (d0, r0) && bitEqual (d1, r1);
            e += peakOf (d0);
        }
        ok (e > 0.5, "precondition: the stream carried signal");
        ok (equal, "an empty call at a narrower width moves no edge — the stream is bit-identical");
    }

    // (c) A call that processes no CHANNEL is the same kind of non-event, and the width is caller-supplied
    //     and unclamped from below: a negative one would make the drop's half-open ranges start at a
    //     negative column and write bqST_[s][-1] — inside the object, where a sanitizer sees nothing.
    {
        EqBand dut, ref;
        dut.prepare (kFs, 2); ref.prepare (kFs, 2);
        BandParams p = stPlusSide (200.0, 12.0);
        p.dyn.on = true;
        dut.setParams (p); ref.setParams (p);
        dut.setLaneDeltaDb (Lane::Side, -6.0); ref.setLaneDeltaDb (Lane::Side, -6.0);

        std::vector<float> d0 (128), d1 (128), r0 (128), r1 (128);
        float* dch[2] { d0.data(), d1.data() };
        float* rch[2] { r0.data(), r1.data() };

        // (a) A NEGATIVE width is MALFORMED: refused, and nothing moved — so the stream is bit-identical
        //     to one that never saw the call. This is the half of the old `nc > 0` guard that survives.
        {
            bool equal = true; double e = 0.0;
            for (int k = 0; k < 16; ++k)
            {
                fillNoise (d0, 21u + (unsigned) k); fillNoise (d1, 71u + (unsigned) k);
                r0 = d0; r1 = d1;
                if (k == 6) ok (! dut.processBlock (dch, -1, 128), "a NEGATIVE width is refused (law 11)");
                felitronics::test::run (dut.processBlock (dch, 2, 128));
                felitronics::test::run (ref.processBlock (rch, 2, 128));
                equal = equal && bitEqual (d0, r0) && bitEqual (d1, r1);
                e += peakOf (d0);
            }
            ok (e > 0.5, "precondition: the stream carried signal");
            ok (equal, "a refused (negative-width) call leaves the stream bit-identical, and indexes nothing");
        }

        // (b) A ZERO width with samples is NOT a no-op, and this is the P20 change. It spends audio time
        //     (the grid advances) and every cell stopped for those samples, so the falling edge fires —
        //     law 11(a). Under the old rule ("a call that processes no channel stopped nothing") the
        //     stream stayed bit-identical, and that is exactly how a gap replays: the sibling case,
        //     measured on Compressor, emitted -11.05 dBFS out of digital silence after such a stretch.
        {
            EqBand dut2, ref2;
            (void) dut2.prepare (kFs, 2); (void) ref2.prepare (kFs, 2);
            dut2.setParams (p); ref2.setParams (p);
            dut2.setLaneDeltaDb (Lane::Side, -6.0); ref2.setLaneDeltaDb (Lane::Side, -6.0);
            double worst = 0.0, e = 0.0;
            int    differing = 0, total = 0;
            for (int k = 0; k < 16; ++k)
            {
                fillNoise (d0, 21u + (unsigned) k); fillNoise (d1, 71u + (unsigned) k);
                r0 = d0; r1 = d1;
                if (k == 6) felitronics::test::run (dut2.processBlock (dch, 0, 128));   // the GAP
                felitronics::test::run (dut2.processBlock (dch, 2, 128));
                felitronics::test::run (ref2.processBlock (rch, 2, 128));
                for (int i = 0; i < 128; ++i)
                {
                    worst = std::max (worst, (double) std::fabs (d0[(std::size_t) i] - r0[(std::size_t) i]));
                    worst = std::max (worst, (double) std::fabs (d1[(std::size_t) i] - r1[(std::size_t) i]));
                    if (d0[(std::size_t) i] != r0[(std::size_t) i] || d1[(std::size_t) i] != r1[(std::size_t) i]) ++differing;
                    ++total;
                }
                e += peakOf (d0);
            }
            std::printf ("      [zero-width gap] %d of %d samples differ, worst |delta| = %.6g\n", differing, total, worst);
            ok (e > 0.5, "precondition: the stream carried signal");
            ok (differing > 0, "a ZERO-width call with samples is a GAP: it stops every cell and moves the stream");
            ok (worst > 1e-6, "...and the difference is audio-sized, not a rounding artefact");
        }
    }
}

//==================================================================================================
// 6. The primitive the isolation invariant rests on.
//==================================================================================================
static void testSvfResetChannel()
{
    group ("Svf::resetChannel clears one column and leaves the others bit-exact");

    Svf a, b;
    a.prepare (kFs, 4); b.prepare (kFs, 4);
    a.setParams (FilterType::Bell, 800.0, 2.0, 9.0);
    b.setParams (FilterType::Bell, 800.0, 2.0, 9.0);

    double last[4] {};
    for (int n = 0; n < 500; ++n)
        for (int c = 0; c < 4; ++c)
        {
            const float x = (float) std::sin (0.05 * n + 0.7 * c);
            last[c] = a.processSample (c, x);
            (void) b.processSample (c, x);
        }
    ok (std::fabs (last[1]) > 1e-3, "precondition: the columns really carry state");

    a.resetChannel (2);
    bool others = true;
    for (int n = 0; n < 64; ++n)
        for (int c = 0; c < 4; ++c)
        {
            const float x = (float) std::sin (0.05 * (n + 500) + 0.7 * c);
            const float ya = a.processSample (c, x);
            const float yb = b.processSample (c, x);
            if (c != 2) others = others && (ya == yb);
        }
    ok (others, "columns 0, 1 and 3 are bit-identical to an instance that was never touched");

    Svf fresh; fresh.prepare (kFs, 4);
    fresh.setParams (FilterType::Bell, 800.0, 2.0, 9.0);
    Svf one; one.prepare (kFs, 4);
    one.setParams (FilterType::Bell, 800.0, 2.0, 9.0);
    for (int n = 0; n < 200; ++n) (void) one.processSample (0, (float) std::sin (0.03 * n));
    one.resetChannel (0);
    bool same = true;
    for (int n = 0; n < 64; ++n)
    {
        const float x = (float) std::cos (0.02 * n);
        same = same && (one.processSample (0, x) == fresh.processSample (0, x));
    }
    ok (same, "a reset column behaves exactly like a never-used one");

    a.resetChannel (-1); a.resetChannel (99);               // out of range is a no-op, not an overrun
    ok (true, "an out-of-range channel index is refused without touching memory");
}

//==================================================================================================
// 7. The section-input capture has a WIDTH, and it used to record only a length.
//==================================================================================================
static void testCaptureWidth()
{
    group ("EqEngine::captureSectionInput reports its width and refuses the columns outside it");

    const int N = 64;
    EqEngine e; felitronics::test::run (e.prepare (kFs, N, 4));
    std::vector<float> v[4];
    float* ch[4] {};
    const float* cch[4] {};
    for (int c = 0; c < 4; ++c)
    {
        v[c].assign ((std::size_t) N, 0.25f + 0.1f * (float) c);
        ch[c] = v[c].data(); cch[c] = v[c].data();
    }

    const float* const* wide = e.captureSectionInput (cch, 4, N);
    ok (wide != nullptr && e.sectionInputChannels() == 4, "precondition: a four-channel capture reports width 4");
    ok (wide[3] != nullptr && wide[3][0] == v[3][0], "precondition: the widest column really was captured");

    const float* const* narrow = e.captureSectionInput (cch, 2, N);
    ok (narrow != nullptr && e.sectionInputChannels() == 2, "a narrower capture reports its own width");
    ok (narrow[2] == nullptr && narrow[3] == nullptr,
        "the columns outside it are refused, not left pointing at the PREVIOUS block's audio");
    ok (narrow[0] != nullptr && narrow[1] != nullptr, "and the captured columns are still there");

    const float* const* again = e.captureSectionInput (cch, 4, N);
    ok (again != nullptr && again[3] != nullptr && e.sectionInputChannels() == 4,
        "a wider capture afterwards restores every column");
}

//==================================================================================================
int main()
{
    std::printf ("felitronics::eq — execution-gate state\n");

    testChannelGateSilence();
    testPreventionOnSurvivingChannel();
    testRecoveryCoefficients();
    testRecoveryStrongNull();
    testNonChannelGates();
    testRampIsWallClock();
    testWhatSurvives();
    testSvfResetChannel();
    testCaptureWidth();

    group ("RT-safety");
    {
        EqBand b; b.prepare (kFs, 4);
        BandParams p = stPlusSide (200.0, 12.0);
        p.dyn.on = true; b.setParams (p);
        std::vector<float> v[4];
        float* ch[4] {};
        for (int c = 0; c < 4; ++c) { v[c].assign (512, 0.1f); ch[c] = v[c].data(); }
        felitronics::test::run (b.processBlock (ch, 2, 512));                          // warm: designs, applies, settles
        const int before = g_allocs.load();
        for (int k = 0; k < 40; ++k) felitronics::test::run (b.processBlock (ch, (k % 3 == 0) ? 1 : ((k % 3 == 1) ? 3 : 2), 512));
        felitronics::test::okNoAlloc (g_allocs.load() == before, "no allocation across 40 blocks of changing width");
    }

    return felitronics::test::report();
}

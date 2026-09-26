// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// LaneDynamics' RELEASE ON DISENGAGE (opt-in; see the note above the class). Each group pins one claim:
//
//   * OFF BY DEFAULT: a point switched off mid-duck still snaps its delta to 0 on the edge, exactly as before.
//   * OPTED IN, the edge RELEASES: the delta decays monotonically through the lane's own release follower, the
//     band's dynamic seam is held open for exactly as long (band.params().dyn.on), and when every delta has arrived
//     the point disengages as it would have and the caller's `dyn.on` is back on the band — from then on a band
//     whose ducked lane is the only one it runs renders bit for bit what a band that never had dynamics renders (a
//     lane downstream of the ducked one keeps a recursive filter's memory of the duck; see the class note).
//   * THE DETECTOR KEEPS LISTENING through a release, lanes keep their participation rules, a refused call moves
//     nothing, and reset() hands a held band back.
//   * EVERY DISENGAGE EDGE releases: dyn.on off, rangeDb 0, and a missing sidechain (the mastering chain's case when
//     no other point is armed); a clock-only pause spends the release; switching back on mid-release carries on from
//     where it stands.
//   * THE EDGE NO LONGER STEPS: its second difference is a fraction of the snap's.
//   * NOTHING IS ALLOCATED.

#include <felitronics_test.h>
#include <alloc_counter.h>   // installs the allocation counter: EVERY form of `new`, over-aligned included
#include <felitronics/dynamiceq/LaneDynamics.h>

#include <bit>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

using namespace felitronics;
using felitronics::test::ok;
using felitronics::test::group;

namespace
{
constexpr double kFs = 48000.0;

eq::BandParams bell (bool dynOn, double rangeDb = -9.0)
{
    eq::BandParams p;
    p.on = true; p.type = eq::FilterType::Bell;
    auto& st = p.lane (eq::Lane::Stereo);
    st.on = true; st.freq = 220.0; st.Q = 1.0; st.gainDb = 0.0;
    p.dyn.on = dynOn; p.dyn.rangeDb = rangeDb; p.dyn.thrAuto = false; p.dyn.thrDb = -40.0;
    return p;
}

struct Rig
{
    eq::EqBand band;
    dynamiceq::LaneDynamics dyn;
    std::vector<float> L, R;
    long long pos = 0;
    bool init (const eq::BandParams& p, bool release)
    {
        if (! band.prepare (kFs, 2) || ! dyn.prepare (kFs, 2)) return false;
        band.setParams (p); dyn.setParams (p);
        dyn.setReleaseOnDisengage (release);
        return true;
    }
    void write (const eq::BandParams& p) { band.setParams (p); dyn.setParams (p); }
    // `n` samples of the -12 dBFS 227 Hz tone, in blocks of `blk`; `sidechain` false hands a null key.
    bool run (int n, int blk = 128, bool sidechain = true)
    {
        bool ok2 = true;
        for (int done = 0; done < n; done += blk)
        {
            const int m = std::min (blk, n - done);
            std::vector<float> a ((std::size_t) m), b ((std::size_t) m);
            for (int i = 0; i < m; ++i) a[(std::size_t) i] = b[(std::size_t) i] = (float) (0.25 * std::sin (2.0 * core::kPi * 227.3 * (double) (pos + i) / kFs));
            std::vector<float> sa = a, sb = b;
            float* io[2] { a.data(), b.data() };
            const float* sc[2] { sa.data(), sb.data() };
            ok2 = ok2 && dyn.processBand (io, sidechain ? sc : nullptr, 2, m, band);
            L.insert (L.end(), a.begin(), a.end()); R.insert (R.end(), b.begin(), b.end());
            pos += m;
        }
        return ok2;
    }
};

double maxD2 (const std::vector<float>& y, int a, int b)
{
    double m = 0.0;
    for (int i = std::max (a, 2); i < b; ++i)
        m = std::max (m, std::fabs ((double) y[(std::size_t) i] - 2.0 * y[(std::size_t) i - 1] + y[(std::size_t) i - 2]));
    return m;
}
} // namespace

static void testDefaultSnaps()
{
    group ("off by default — a point switched off mid-duck still snaps, exactly as before");
    Rig r;
    ok (r.init (bell (true), false), "PRECONDITION: prepare");
    ok (r.run (48000), "PRECONDITION: a second of ducking");
    ok (r.dyn.deltaDb (eq::Lane::Stereo) < -6.0, "PRECONDITION: it is ducking (" + std::to_string (r.dyn.deltaDb (eq::Lane::Stereo)) + " dB)");
    r.write (bell (false));
    ok (r.run (128), "PRECONDITION: a block after the edge");
    ok (! r.dyn.isReleasing() && r.dyn.deltaDb (eq::Lane::Stereo) == 0.0, "the delta is 0 at once and nothing releases");
}

static void testReleaseOnDisengage()
{
    group ("opted in — the edge releases through the follower, holds the seam open, and hands the band back");
    Rig r, still;
    ok (r.init (bell (true), true) && still.init (bell (false), true), "PRECONDITION: prepare");
    ok (r.run (48000) && still.run (48000), "PRECONDITION: a second");
    const double d0 = r.dyn.deltaDb (eq::Lane::Stereo);
    ok (d0 < -6.0, "PRECONDITION: ducking (" + std::to_string (d0) + " dB)");
    r.write (bell (false));
    bool monotone = true, heldWhileReleasing = true;
    double prev = d0;
    int finishedAt = -1;
    for (int k = 0; k < 2000 && finishedAt < 0; ++k)
    {
        felitronics::test::run (r.run (64) && still.run (64));
        const double d = r.dyn.deltaDb (eq::Lane::Stereo);
        monotone = monotone && d >= prev && d <= 0.0;
        prev = d;
        if (r.dyn.isReleasing()) heldWhileReleasing = heldWhileReleasing && r.band.params().dyn.on;
        else finishedAt = (int) r.pos;
    }
    ok (monotone, "the delta rises monotonically toward 0 dB through the release");
    ok (heldWhileReleasing, "the band's dynamic seam is held open for as long as the release runs");
    ok (finishedAt > 48000 + 64, "the release took time (finished at " + std::to_string (finishedAt) + ")");
    ok (finishedAt > 0 && finishedAt < 48000 + 48000, "…and finished within a second");
    ok (! r.band.params().dyn.on && r.dyn.deltaDb (eq::Lane::Stereo) == 0.0, "then the caller's dyn.on is back on the band and the delta is 0");
    ok (r.run (4800) && still.run (4800), "PRECONDITION: more audio");
    long long d = 0;
    for (std::size_t i = (std::size_t) finishedAt; i < r.L.size(); ++i)
        d += std::bit_cast<std::uint32_t> (r.L[i]) != std::bit_cast<std::uint32_t> (still.L[i])
          || std::bit_cast<std::uint32_t> (r.R[i]) != std::bit_cast<std::uint32_t> (still.R[i]);
    ok (d == 0, "from the moment it finished, this single-lane band renders bit for bit what a band that never had dynamics "
                "renders (" + std::to_string (d) + " differ)");
}

static void testEveryEdgeReleases()
{
    group ("every disengage edge releases — rangeDb 0, a missing sidechain, a pause; and back on carries on");
    {
        Rig r;
        ok (r.init (bell (true), true) && r.run (48000), "PRECONDITION: ducking");
        r.write (bell (true, 0.0));
        ok (r.run (64) && r.dyn.isReleasing() && r.dyn.deltaDb (eq::Lane::Stereo) < -5.0, "rangeDb -> 0 releases instead of snapping");
    }
    {
        Rig r;
        ok (r.init (bell (true), true) && r.run (48000), "PRECONDITION: ducking");
        ok (r.run (64, 64, false) && r.dyn.isReleasing() && r.dyn.deltaDb (eq::Lane::Stereo) < -5.0, "a missing sidechain releases instead of snapping");
    }
    {
        Rig r, s;
        ok (r.init (bell (true), true) && s.init (bell (true), true) && r.run (48000) && s.run (48000), "PRECONDITION: ducking");
        r.write (bell (false)); s.write (bell (false));
        ok (r.run (64) && s.run (64), "PRECONDITION: the release has started");
        // The same 4800 samples: one rig as a clock-only pause, the other as silence at width 2 with a silent key —
        // a pause IS silence, so the release must stand at the same place.
        felitronics::test::run (r.dyn.processBand (nullptr, nullptr, 0, 4800, r.band));
        std::vector<float> z (4800, 0.0f), z2 (4800, 0.0f), k1 (4800, 0.0f), k2 (4800, 0.0f);
        float* io[2] { z.data(), z2.data() };
        const float* sc[2] { k1.data(), k2.data() };
        for (int o = 0; o < 4800; o += 128)
        {
            float* p[2] { io[0] + o, io[1] + o };
            const float* q[2] { sc[0] + o, sc[1] + o };
            felitronics::test::run (s.dyn.processBand (p, q, 2, std::min (128, 4800 - o), s.band));
        }
        ok (r.dyn.deltaDb (eq::Lane::Stereo) > -5.0 && std::bit_cast<std::uint64_t> (r.dyn.deltaDb (eq::Lane::Stereo))
                                                     == std::bit_cast<std::uint64_t> (s.dyn.deltaDb (eq::Lane::Stereo)),
            "a clock-only pause spends the release exactly as silence does (" + std::to_string (r.dyn.deltaDb (eq::Lane::Stereo)) + " dB)");
    }
    {
        Rig r;
        ok (r.init (bell (true), true) && r.run (48000), "PRECONDITION: ducking");
        r.write (bell (false));
        ok (r.run (1024), "PRECONDITION: releasing");
        const double mid = r.dyn.deltaDb (eq::Lane::Stereo);
        r.write (bell (true));
        ok (r.run (64) && ! r.dyn.isReleasing(), "switched back on mid-release: engaged again");
        ok (std::fabs (r.dyn.deltaDb (eq::Lane::Stereo) - mid) < 1.0, "…from where the release stood, not from 0 or from the old duck ("
                                                                     + std::to_string (mid) + " -> " + std::to_string (r.dyn.deltaDb (eq::Lane::Stereo)) + ")");
    }
}


// FOUND BY THE CODE-REVIEW ROUND, each pinned where it was found.
static void testReviewFindings()
{
    group ("the detector keeps listening, lanes keep their rules, a refused call moves nothing, reset() hands the band back");
    // reset() mid-release: the held band gets its deltas and the caller's dyn.on back.
    {
        Rig r;
        ok (r.init (bell (true), true) && r.run (48000), "PRECONDITION: ducking");
        r.write (bell (false));
        ok (r.run (256) && r.dyn.isReleasing() && r.band.params().dyn.on, "PRECONDITION: releasing, the seam held open");
        r.dyn.reset();
        ok (! r.band.params().dyn.on && r.band.laneDeltaDb (eq::Lane::Stereo) == 0.0 && ! r.dyn.isReleasing(),
            "reset() mid-release restores the caller's dyn.on and zeroes the band's delta");
    }
    // The detector listens through the release: the programme goes silent while releasing, the point is switched back on
    // 100 ms later — it must not re-duck from an envelope frozen at the loud passage.
    {
        Rig r;
        ok (r.init (bell (true), true) && r.run (48000), "PRECONDITION: ducking");
        r.write (bell (false));
        ok (r.run (1024), "PRECONDITION: releasing");
        // silence, key included
        std::vector<float> a (128, 0.0f), b (128, 0.0f), sa (128, 0.0f), sb (128, 0.0f);
        float* io[2] { a.data(), b.data() };
        const float* sc[2] { sa.data(), sb.data() };
        for (int k = 0; k < 37; ++k) felitronics::test::run (r.dyn.processBand (io, sc, 2, 128, r.band));   // ~100 ms
        const double before = r.dyn.deltaDb (eq::Lane::Stereo);
        r.write (bell (true));
        for (int k = 0; k < 8; ++k) felitronics::test::run (r.dyn.processBand (io, sc, 2, 128, r.band));
        const double after = r.dyn.deltaDb (eq::Lane::Stereo);
        ok (after >= before - 1e-6, "switched back on over silence: no re-duck (" + std::to_string (before) + " -> " + std::to_string (after) + " dB)");
    }
    // A lane switched off mid-release drops its delta at once, as it would while engaged.
    {
        Rig r;
        ok (r.init (bell (true), true) && r.run (48000), "PRECONDITION: ducking");
        eq::BandParams off = bell (false);
        r.write (off);
        ok (r.run (256) && r.dyn.deltaDb (eq::Lane::Stereo) < -3.0, "PRECONDITION: releasing");
        off.lane (eq::Lane::Stereo).on = false;
        r.write (off);
        ok (r.run (64) && r.dyn.deltaDb (eq::Lane::Stereo) == 0.0, "a lane switched off mid-release drops its delta, as while engaged");
    }
    // A call the band refuses moves nothing: a band prepared for ONE channel, a producer for two, a stereo call on the edge.
    {
        eq::EqBand band; dynamiceq::LaneDynamics dyn;
        ok (band.prepare (kFs, 1) && dyn.prepare (kFs, 2), "PRECONDITION: prepare");
        band.setParams (bell (true)); dyn.setParams (bell (true)); dyn.setReleaseOnDisengage (true);
        std::vector<float> m (128), s2 (128);
        for (int k = 0; k < 400; ++k)
        {
            for (int i = 0; i < 128; ++i) m[(std::size_t) i] = (float) (0.25 * std::sin (2.0 * core::kPi * 227.3 * (double) (k * 128 + i) / kFs));
            std::vector<float> key = m;
            float* io[1] { m.data() }; const float* sc[1] { key.data() };
            felitronics::test::run (dyn.processBand (io, sc, 1, 128, band));
        }
        ok (dyn.deltaDb (eq::Lane::Stereo) < -3.0, "PRECONDITION: the mono point ducks");
        band.setParams (bell (false)); dyn.setParams (bell (false));
        float* io2[2] { m.data(), s2.data() };
        const float* sc2[2] { m.data(), s2.data() };
        ok (! dyn.processBand (io2, sc2, 2, 128, band), "PRECONDITION: the band refuses a stereo call");
        ok (! dyn.isReleasing() && ! band.params().dyn.on, "the refused call started no release and did not open the seam");
    }
}

static void testTheEdgeNoLongerSteps()
{
    group ("the edge no longer steps — a fraction of the snap's second difference");
    Rig snap, rel;
    ok (snap.init (bell (true), false) && rel.init (bell (true), true), "PRECONDITION: prepare");
    ok (snap.run (48000) && rel.run (48000), "PRECONDITION: ducking");
    snap.write (bell (false)); rel.write (bell (false));
    ok (snap.run (9600) && rel.run (9600), "PRECONDITION: past the edge");
    const double s = maxD2 (snap.L, 48000, 48000 + 2000), r = maxD2 (rel.L, 48000, 48000 + 2000), steady = maxD2 (rel.L, 30000, 46000);
    std::printf ("    snap %.2e (%.1f dBFS)  release %.2e (%.1f dBFS)  steady %.2e\n", s, 20.0 * std::log10 (s), r, 20.0 * std::log10 (r), steady);
    ok (r < 0.05 * s, "the release's worst Δ² is under 5% of the snap's");
}

static void testNoAllocation()
{
    group ("RT — nothing is allocated while a release runs");
    Rig r;
    ok (r.init (bell (true), true), "PRECONDITION: prepare");
    r.L.reserve (1 << 20); r.R.reserve (1 << 20);
    ok (r.run (48000), "PRECONDITION: ducking");
    r.write (bell (false));
    // The rig's own per-block vectors allocate; count only what LaneDynamics and the band do, on fixed buffers.
    std::vector<float> a (128, 0.1f), b (128, 0.1f), sa (128, 0.1f), sb (128, 0.1f);
    float* io[2] { a.data(), b.data() };
    const float* sc[2] { sa.data(), sb.data() };
    const long long before = alloc::count.load();
    for (int k = 0; k < 400; ++k)
    {
        felitronics::test::run (r.dyn.processBand (io, sc, 2, 128, r.band));
        if (k == 100) felitronics::test::run (r.dyn.processBand (nullptr, nullptr, 0, 3000, r.band));
    }
    felitronics::test::okNoAlloc (alloc::count.load() == before, "no allocation across 400 releasing blocks and a pause");
}

int main()
{
    std::printf ("felitronics::dynamiceq — release on disengage\n");
    testDefaultSnaps();
    testReleaseOnDisengage();
    testEveryEdgeReleases();
    testReviewFindings();
    testTheEdgeNoLongerSteps();
    testNoAllocation();
    return felitronics::test::report();
}

// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// What a channel of this stage does when it STOPS being fed and is fed again later — plus the stage's
// own second gate, the DC blocker, which opens and closes on the SHAPE at a constant channel count.
//
// Three recursions hide behind one `c < nc`: the oversampler's polyphase histories (with their ring
// POSITIONS, which are state too), the DC blocker's x1/y1, and the dry delay line. None of them decays
// while the channel is absent — they freeze, and the channel replays them on return. Measured before
// this suite: 0.9337 out of DIGITAL SILENCE, 29 samples after a stereo -> mono -> stereo excursion.
//
// Isolation is owed here in full: this stage has no gain shared between channels (comp_, mix_ and
// outGain_ come from parameters, never from the signal), so a channel that never left must keep its
// history bit-exact and a wholesale reset would be the wrong remedy.

#include <felitronics_test.h>
#include <felitronics/saturation/Saturator.h>
#include <felitronics/oversampling/PolyphaseOversampler.h>

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

static void fillNoise (std::vector<float>& v, unsigned seed)
{
    unsigned s = seed;
    for (auto& x : v) { s = s * 1664525u + 1013904223u; x = (float) ((int) (s >> 8) % 2000 - 1000) * 0.0004f; }
}

//==================================================================================================
static void testChannelGate()
{
    group ("channel gate: a returning channel emits nothing from digital silence");

    const int N = 64;
    saturation::Saturator s;
    saturation::Saturator::Params p;
    p.driveDb = 12.0f;
    p.mix     = 0.5f;   // mix < 1 puts the DRY DELAY LINE in the returned channel's path as well; at full
                        // wet it contributes nothing and a frozen delay line would go unmeasured.
    s.setParams (p);
    ok (s.prepare (kFs, N, 2), "precondition: the stage prepared");

    std::vector<float> L (N, 0.0f), R (N, 0.9f);
    float* io[2] { L.data(), R.data() };
    felitronics::test::run (s.process (io, 2, N));
    ok (peakOf (R) > 0.3, "precondition: the right channel really was driven");

    for (int k = 0; k < 20; ++k) { std::fill (L.begin(), L.end(), 0.0f); felitronics::test::run (s.process (io, 1, N)); }

    std::fill (L.begin(), L.end(), 0.0f); std::fill (R.begin(), R.end(), 0.0f);
    felitronics::test::run (s.process (io, 2, N));
    ok (peakOf (R) == 0.0, "silence in, exact zero out on the returned channel (was 0.9337 at i=29)");
    ok (peakOf (L) == 0.0, "and on the channel that stayed, which has heard nothing but silence");
}

static void testIsolation()
{
    group ("isolation: the channel that never left keeps its history bit-exact");

    const int N = 64;
    saturation::Saturator dut, ref;
    saturation::Saturator::Params p;
    p.driveDb = 9.0f; p.mix = 0.6f;                          // mix < 1 puts the dry delay line in play too
    dut.setParams (p); ref.setParams (p);
    ok (dut.prepare (kFs, N, 2) && ref.prepare (kFs, N, 2), "precondition: both stages prepared");

    std::vector<float> d0 (N), d1 (N), r0 (N), r1 (N);
    float* dio[2] { d0.data(), d1.data() };
    float* rio[2] { r0.data(), r1.data() };

    bool equal = true; double energy = 0.0;
    for (int k = 0; k < 60; ++k)
    {
        const bool gap = (k >= 20 && k < 40);
        fillNoise (d0, 31u + (unsigned) k); fillNoise (d1, 707u + (unsigned) k);
        r0 = d0; r1 = d1;
        felitronics::test::run (dut.process (dio, gap ? 1 : 2, N));
        felitronics::test::run (ref.process (rio, 2, N));
        equal = equal && bitEqual (d0, r0);
        energy += peakOf (d0);
    }
    ok (energy > 1.0, "precondition: the surviving channel carried signal");
    ok (equal, "channel 0 is bit-equal to the run where nothing ever left");
}

static void testShapeGate()
{
    group ("the stage's own gate: the DC blocker stops on a symmetric curve and restarts clean");

    // dcEnabled_ is (dcBlockHz > 0) AND shape == Asym, so Asym -> Tanh -> Asym stops and restarts a real
    // IIR recursion at a CONSTANT channel count. Everything else in the chain is an FIR or a delay line,
    // so a long enough silent stretch drains them exactly — which is what makes "exact zero" the right
    // oracle here rather than a tolerance.
    const int N = 64;
    saturation::Saturator s;
    saturation::Saturator::Params p;
    p.shape = saturation::WaveShaper::Shape::Asym;
    p.bias = 0.3f; p.driveDb = 12.0f;
    s.setParams (p);
    ok (s.prepare (kFs, N, 2), "precondition: the stage prepared");

    std::vector<float> a (N), b (N);
    float* io[2] { a.data(), b.data() };
    double charged = 0.0;
    for (int k = 0; k < 20; ++k)
    {
        for (int i = 0; i < N; ++i) a[(std::size_t) i] = b[(std::size_t) i] = (float) (0.8 * std::sin (0.06 * (k * N + i)));
        felitronics::test::run (s.process (io, 2, N));
        charged = peakOf (a);
    }
    ok (charged > 0.1, "precondition: the asymmetric curve really was driving the blocker");

    p.shape = saturation::WaveShaper::Shape::Tanh;           // the gate closes: x1/y1 stop advancing
    s.setParams (p);
    for (int k = 0; k < 40; ++k)
    {
        std::fill (a.begin(), a.end(), 0.0f); std::fill (b.begin(), b.end(), 0.0f);
        felitronics::test::run (s.process (io, 2, N));
    }
    p.shape = saturation::WaveShaper::Shape::Asym;           // and opens again
    s.setParams (p);
    std::fill (a.begin(), a.end(), 0.0f); std::fill (b.begin(), b.end(), 0.0f);
    felitronics::test::run (s.process (io, 2, N));
    ok (peakOf (a) == 0.0 && peakOf (b) == 0.0, "Asym -> Tanh -> Asym: exact zero out of silence");
}

static void testNoEdgeWithoutSamples()
{
    group ("a call that carries no samples is not a transition");

    const int N = 64;
    saturation::Saturator dut, ref;
    saturation::Saturator::Params p; p.driveDb = 6.0f;
    dut.setParams (p); ref.setParams (p);
    ok (dut.prepare (kFs, N, 2) && ref.prepare (kFs, N, 2), "precondition: both stages prepared");

    std::vector<float> d0 (N), d1 (N), r0 (N), r1 (N);
    float* dio[2] { d0.data(), d1.data() };
    float* rio[2] { r0.data(), r1.data() };
    bool equal = true; double e = 0.0;
    for (int k = 0; k < 24; ++k)
    {
        fillNoise (d0, 5u + (unsigned) k); fillNoise (d1, 55u + (unsigned) k);
        r0 = d0; r1 = d1;
        if (k == 12) felitronics::test::run (dut.process (dio, 1, 0));                // an empty probe at a narrower width
        felitronics::test::run (dut.process (dio, 2, N));
        felitronics::test::run (ref.process (rio, 2, N));
        equal = equal && bitEqual (d0, r0) && bitEqual (d1, r1);
        e += peakOf (d0);
    }
    ok (e > 0.5, "precondition: the stream carried signal");
    ok (equal, "the stream is bit-identical across the empty probe");
}

static void testRePrepareNarrower()
{
    group ("a narrower re-prepare must not leave a ledger pointing past the new buffers");

    // prepare() REALLOCATES every per-channel vector and does NOT call reset(), so a run ledger left over
    // from a wider previous life indexes storage that no longer exists — an AddressSanitizer
    // container-overflow on dryDelay_, reached by nothing more exotic than a host changing its bus width.
    const int N = 64;
    saturation::Saturator reused, fresh;
    saturation::Saturator::Params p; p.driveDb = 9.0f; p.mix = 0.5f;
    reused.setParams (p); fresh.setParams (p);

    ok (reused.prepare (kFs, N, 2), "precondition: the stage prepared wide");
    std::vector<float> a (N), b (N);
    float* io2[2] { a.data(), b.data() };
    for (int k = 0; k < 8; ++k) { fillNoise (a, 3u + (unsigned) k); fillNoise (b, 33u + (unsigned) k); felitronics::test::run (reused.process (io2, 2, N)); }
    ok (peakOf (a) > 0.0, "precondition: the wide life really ran");

    ok (reused.prepare (kFs, N, 1), "the stage re-prepares narrower");
    ok (fresh.prepare  (kFs, N, 1), "and a fresh instance prepares the same way");

    // Everything prepare() touches is reallocated zeroed, so a correctly re-prepared instance is
    // indistinguishable from a new one — which is both the safety check and the behaviour check.
    bool equal = true;
    std::vector<float> u (N), v (N);
    float* iou[1] { u.data() }; float* iov[1] { v.data() };
    for (int k = 0; k < 12; ++k)
    {
        fillNoise (u, 500u + (unsigned) k); v = u;
        felitronics::test::run (reused.process (iou, 1, N));
        felitronics::test::run (fresh.process  (iov, 1, N));
        equal = equal && bitEqual (u, v);
    }
    ok (equal, "the re-prepared instance is bit-identical to a fresh one");
}

static void testOversamplerResetChannel()
{
    group ("PolyphaseOversampler::resetChannel clears the ring AND its position");

    const int N = 32, L = 4;
    oversampling::PolyphaseOversampler a, fresh;
    ok (a.prepare (L, 3, 32) && fresh.prepare (L, 3, 32), "precondition: both oversamplers prepared");

    std::vector<float> in0 (N), in1 (N), in2 (N);
    std::vector<float> up0 ((std::size_t) N * L), up1 ((std::size_t) N * L), up2 ((std::size_t) N * L);
    float* cin[3] { in0.data(), in1.data(), in2.data() };
    float* cup[3] { up0.data(), up1.data(), up2.data() };

    // The block length must NOT divide the ring: 8 blocks of 32 through a 32-slot ring land the cursor
    // back on zero, and a test run from a zeroed cursor cannot tell a cleared cursor from an uncleared
    // one no matter what it asserts. 30 is coprime enough to leave it elsewhere. The DOWNsampler is
    // charged too — its ring is N = L*tpp long and has its own cursor.
    const int M = 30;
    std::vector<float> dn0 (M), dn1 (M), dn2 (M);
    float* cdn[3] { dn0.data(), dn1.data(), dn2.data() };
    for (int k = 0; k < 8; ++k)                              // charge every column with something different
    {
        fillNoise (in0, 1u + (unsigned) k); fillNoise (in1, 9001u + (unsigned) k); fillNoise (in2, 4242u + (unsigned) k);
        a.upsample (cin, 3, M, cup);
        const float* down[3] { up0.data(), up1.data(), up2.data() };
        a.downsample (down, 3, M, cdn);
    }
    ok (peakOf (up1) > 0.0, "precondition: the columns really carry history");

    a.resetChannel (1);

    // A reset column must behave exactly like one from a never-used instance. Note what this does NOT
    // prove: leaving the cursor where it was is UNOBSERVABLE here and everywhere else through this API,
    // because every read is relative to the post-write cursor, so a uniformly zero ring is the same ring
    // at any rotation. resetChannel() clears the cursor anyway — a column whose ring is zeroed while its
    // cursor sits elsewhere is not the fresh column it claims to be — but no test can fail on it, and
    // saying otherwise would be a test that certifies its own comment.
    std::vector<float> t (N, 0.0f); t[0] = 1.0f;             // an impulse reads the whole history out
    std::vector<float> zero (N, 0.0f);
    std::vector<float> aOut ((std::size_t) N * L), fOut ((std::size_t) N * L);
    float* aIn[3] { t.data(), t.data(), t.data() };
    float* aUp[3] { up0.data(), aOut.data(), up2.data() };
    float* fIn[3] { t.data(), t.data(), t.data() };
    float* fUp[3] { up0.data(), fOut.data(), up2.data() };
    a.upsample (aIn, 3, N, aUp);
    fresh.upsample (fIn, 3, N, fUp);
    ok (bitEqual (aOut, fOut), "the reset column is bit-identical to a never-used one, cursor included");

    a.resetChannel (-1); a.resetChannel (999);
    ok (true, "an out-of-range channel index is refused without touching memory");
}

//==================================================================================================
int main()
{
    std::printf ("felitronics::saturation — execution-gate state\n");

    testChannelGate();
    testRePrepareNarrower();
    testIsolation();
    testShapeGate();
    testNoEdgeWithoutSamples();
    testOversamplerResetChannel();

    group ("RT-safety");
    {
        const int N = 128;
        saturation::Saturator s;
        saturation::Saturator::Params p; p.shape = saturation::WaveShaper::Shape::Asym; p.mix = 0.5f;
        s.setParams (p);
        ok (s.prepare (kFs, N, 4), "precondition: the stage prepared");
        std::vector<float> v[4];
        float* io[4] {};
        for (int c = 0; c < 4; ++c) { v[c].assign ((std::size_t) N, 0.05f); io[c] = v[c].data(); }
        felitronics::test::run (s.process (io, 4, N));
        const int before = g_allocs.load();
        for (int k = 0; k < 40; ++k) felitronics::test::run (s.process (io, (k % 3) + 2, N));
        felitronics::test::okNoAlloc (g_allocs.load() == before, "no allocation across 40 blocks of changing width");
    }

    return felitronics::test::report();
}

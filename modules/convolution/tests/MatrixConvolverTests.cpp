// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// JUCE-free self-tests for the MatrixConvolver: the 2×2 OPERATOR convolver on one canonical raw-L/R
// history. Reference-NULL against the proven PartitionedConvolver: each routing topology is nulled against
// an INDEPENDENT time-domain computation of the same math, so the re-plumb (raw-L/R FDL + on-the-fly M/S
// views + cross-input sums) is proven, not assumed. Plus: an atomic topology SWITCH mid-playback proven
// click-free AND warm (its incoming operator carries the input's past, unlike a cold-started instance),
// and no allocation in process() even across a swap.

#include <felitronics_test.h>
#include <felitronics/convolution/MatrixConvolver.h>
#include <felitronics/convolution/PartitionedConvolver.h>   // reference convolver for the null tests

#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

static std::atomic<long> g_allocs { 0 };
void* operator new      (std::size_t s) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s ? s : 1); }
void* operator new[]    (std::size_t s) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s ? s : 1); }
void  operator delete   (void* p) noexcept { std::free (p); }
void  operator delete[] (void* p) noexcept { std::free (p); }
void  operator delete   (void* p, std::size_t) noexcept { std::free (p); }
void  operator delete[] (void* p, std::size_t) noexcept { std::free (p); }

using namespace felitronics;
using MC = convolution::MatrixConvolver<>;

struct Lcg { unsigned long long s; float next() { s = s * 6364136223846793005ULL + 1442695040888963407ULL; return (float) ((s >> 40) & 0xffff) / 32768.0f - 1.0f; } };

// One mono reference convolution (in ∗ ir) via the proven PartitionedConvolver, from zero history.
static std::vector<float> convRef (const std::vector<float>& ir, const std::vector<float>& in, int P, int irMax)
{
    convolution::PartitionedConvolver<> pc; pc.prepare (P, irMax); pc.setIr (ir.data(), (int) ir.size());
    std::vector<float> out (in.size(), 0.0f);
    pc.process (in.data(), out.data(), (int) in.size());
    return out;
}

static double maxDiff (const std::vector<float>& a, const std::vector<float>& b, int from, int to)
{
    double m = 0.0; for (int i = from; i < to; ++i) m = std::max (m, (double) std::fabs (a[(std::size_t) i] - b[(std::size_t) i])); return m;
}
static double maxDeriv (const std::vector<float>& y, int from, int to)
{
    double m = 0.0; for (int i = from + 1; i < to; ++i) m = std::max (m, (double) std::fabs (y[(std::size_t) i] - y[(std::size_t) (i - 1)])); return m;
}

int main()
{
    std::printf ("felitronics::convolution MatrixConvolver tests\n");
    const int P = 64, irMax = 400, len = 200, xfade = 128, n = 4000;
    const int settled = 900;                              // past coldXfade (≈384) + tail fill + margin

    Lcg r { 2718 };
    auto mkIr = [&] { std::vector<float> v ((std::size_t) len); for (auto& x : v) x = 0.15f * r.next(); return v; };
    std::vector<float> irM = mkIr(), irS = mkIr(), irL = mkIr(), irR = mkIr();
    std::vector<float> irLL = mkIr(), irLR = mkIr(), irRL = mkIr(), irRR = mkIr();

    std::vector<float> xL (n), xR (n);
    for (int i = 0; i < n; ++i) { xL[(std::size_t) i] = (float) (0.5 * std::sin (2.0 * core::kPi * 500.0 * i / 48000.0)); xR[(std::size_t) i] = 0.4f * r.next(); }

    auto runStereo = [&] (MC& mc, std::vector<float>& oL, std::vector<float>& oR)
    {
        oL.assign ((std::size_t) n, 0.0f); oR.assign ((std::size_t) n, 0.0f);
        for (int i = 0; i < n; ++i) { oL[(std::size_t) i] = xL[(std::size_t) i]; oR[(std::size_t) i] = xR[(std::size_t) i]; }
        for (int o = 0; o < n; o += 512) { float* io[2] { oL.data() + o, oR.data() + o }; mc.process (io, io, 2, std::min (512, n - o)); }
    };

    // --- LRDiag == two independent mono convolutions (direct L/R routing) ---
    test::group ("MatrixConvolver LRDiag == two independent mono convolutions");
    {
        MC mc; test::ok (mc.prepare (P, irMax, xfade, 2), "prepare stereo");
        const float* banks[2] { irL.data(), irR.data() };
        test::ok (mc.setOperator (MC::Topology::LRDiag, banks, 2, len), "setOperator LRDiag");
        std::vector<float> oL, oR; runStereo (mc, oL, oR);
        const std::vector<float> refL = convRef (irL, xL, P, irMax), refR = convRef (irR, xR, P, irMax);
        test::ok (maxDiff (oL, refL, settled, n) < 1e-4, "yL == irL ∗ xL (independent of xR)");
        test::ok (maxDiff (oR, refR, settled, n) < 1e-4, "yR == irR ∗ xR (independent of xL)");
    }

    // --- MSDiag == time-domain encode→conv→decode (the spectral-view re-plumb proof) ---
    test::group ("MatrixConvolver MSDiag == encode→conv→decode reference");
    {
        MC mc; mc.prepare (P, irMax, xfade, 2);
        const float* banks[2] { irM.data(), irS.data() };
        test::ok (mc.setOperator (MC::Topology::MSDiag, banks, 2, len), "setOperator MSDiag");
        std::vector<float> oL, oR; runStereo (mc, oL, oR);
        std::vector<float> m (n), s (n);
        for (int i = 0; i < n; ++i) { m[(std::size_t) i] = 0.5f * (xL[(std::size_t) i] + xR[(std::size_t) i]); s[(std::size_t) i] = 0.5f * (xL[(std::size_t) i] - xR[(std::size_t) i]); }
        const std::vector<float> yM = convRef (irM, m, P, irMax), yS = convRef (irS, s, P, irMax);
        std::vector<float> refL (n), refR (n);
        for (int i = 0; i < n; ++i) { refL[(std::size_t) i] = yM[(std::size_t) i] + yS[(std::size_t) i]; refR[(std::size_t) i] = yM[(std::size_t) i] - yS[(std::size_t) i]; }
        const double eL = maxDiff (oL, refL, settled, n), eR = maxDiff (oR, refR, settled, n);
        std::printf ("      MSDiag view-vs-time-encode max err L=%.2e R=%.2e\n", eL, eR);
        test::ok (eL < 1e-5, "yL == (yM+yS) — M/S spectral view == time-domain encode (float-exact)");
        test::ok (eR < 1e-5, "yR == (yM−yS) — the on-the-fly view re-plumb is exact");
    }

    // --- Full == direct 4-conv cross sums (cross-input routing + polarity) ---
    test::group ("MatrixConvolver Full == 4-conv cross-sum reference");
    {
        MC mc; mc.prepare (P, irMax, xfade, 2);
        const float* banks[4] { irLL.data(), irLR.data(), irRL.data(), irRR.data() };
        test::ok (mc.setOperator (MC::Topology::Full, banks, 4, len), "setOperator Full (4 banks)");
        std::vector<float> oL, oR; runStereo (mc, oL, oR);
        const std::vector<float> cLL = convRef (irLL, xL, P, irMax), cLR = convRef (irLR, xR, P, irMax);
        const std::vector<float> cRL = convRef (irRL, xL, P, irMax), cRR = convRef (irRR, xR, P, irMax);
        std::vector<float> refL (n), refR (n);
        for (int i = 0; i < n; ++i) { refL[(std::size_t) i] = cLL[(std::size_t) i] + cLR[(std::size_t) i]; refR[(std::size_t) i] = cRL[(std::size_t) i] + cRR[(std::size_t) i]; }
        test::ok (maxDiff (oL, refL, settled, n) < 5e-4, "yL == LL∗xL + LR∗xR");
        test::ok (maxDiff (oR, refR, settled, n) < 5e-4, "yR == RL∗xL + RR∗xR");
        // off-diagonal actually routes: dropping xR must change yL (LR∗xR term is real)
        double lrEnergy = 0.0; for (int i = settled; i < n; ++i) lrEnergy += std::fabs (cLR[(std::size_t) i]);
        test::ok (lrEnergy > 1e-2, "the LR cross term is non-trivial (Full genuinely mixes channels)");
    }

    // --- topology SWITCH mid-playback (2-bank MSDiag → 4-bank Full) : click-free + WARM history ---
    test::group ("MatrixConvolver topology switch MSDiag→Full: click-free + warm history");
    {
        const int T = 2000;                              // switch point (well past the cold prime → warm)
        // A: warm switch. Run MSDiag to T, then switch to Full mid-stream, keep streaming.
        MC A; A.prepare (P, irMax, xfade, 2);
        { const float* b[2] { irM.data(), irS.data() }; A.setOperator (MC::Topology::MSDiag, b, 2, len); }
        std::vector<float> aL (n), aR (n);
        for (int i = 0; i < n; ++i) { aL[(std::size_t) i] = xL[(std::size_t) i]; aR[(std::size_t) i] = xR[(std::size_t) i]; }
        bool switched = false;
        for (int o = 0; o < n; o += 256)
        {
            if (! switched && o >= T)
            { const float* b[4] { irLL.data(), irLR.data(), irRL.data(), irRR.data() };
              if (A.setOperator (MC::Topology::Full, b, 4, len)) switched = true; }
            float* io[2] { aL.data() + o, aR.data() + o }; A.process (io, io, 2, std::min (256, n - o));
        }
        test::ok (switched, "topology switch accepted mid-stream");

        // click-free: the derivative across the switch stays near the steady-state derivative.
        const double steady = maxDeriv (aL, T - 400, T - 50);
        const double across = maxDeriv (aL, T, T + xfade + 300);
        test::ok (across < 4.0 * steady + 1e-6, "no click across the topology change");

        // correctness: after the fade, A == the Full operator fed the WHOLE stream (warm history honoured).
        const std::vector<float> cLL = convRef (irLL, xL, P, irMax), cLR = convRef (irLR, xR, P, irMax);
        const std::vector<float> cRL = convRef (irRL, xL, P, irMax), cRR = convRef (irRR, xR, P, irMax);
        std::vector<float> fullL (n), fullR (n);
        for (int i = 0; i < n; ++i) { fullL[(std::size_t) i] = cLL[(std::size_t) i] + cLR[(std::size_t) i]; fullR[(std::size_t) i] = cRL[(std::size_t) i] + cRR[(std::size_t) i]; }
        const int after = T + xfade + P + 8;
        test::ok (maxDiff (aL, fullL, after, n) < 5e-4, "A converges to the true Full response (switched correctly)");

        // WARM proof: a COLD instance fed ONLY the post-switch samples differs from A right after the switch —
        // A's incoming Full operator reads the shared warm FDL (the input's past), the cold one cannot.
        MC B; B.prepare (P, irMax, xfade, 2);
        { const float* b[4] { irLL.data(), irLR.data(), irRL.data(), irRR.data() }; B.setOperator (MC::Topology::Full, b, 4, len); }
        const int tail = n - T;
        std::vector<float> bL (tail), bR (tail);
        for (int i = 0; i < tail; ++i) { bL[(std::size_t) i] = xL[(std::size_t) (T + i)]; bR[(std::size_t) i] = xR[(std::size_t) (T + i)]; }
        for (int o = 0; o < tail; o += 256) { float* io[2] { bL.data() + o, bR.data() + o }; B.process (io, io, 2, std::min (256, tail - o)); }
        // compare A[T + after..] to B[after..] over a window where the pre-switch tail still matters
        double warmDiff = 0.0;
        for (int i = xfade + P; i < xfade + P + 400; ++i)
            warmDiff = std::max (warmDiff, (double) std::fabs (aL[(std::size_t) (T + i)] - bL[(std::size_t) i]));
        test::ok (warmDiff > 1e-2, "A (warm) ≠ a cold instance fed only post-switch samples → shared warm FDL");
    }

    // --- MONO degenerate == a single PartitionedConvolver ---
    test::group ("MatrixConvolver mono == one PartitionedConvolver");
    {
        MC mc; test::ok (mc.prepare (P, irMax, xfade, 1), "prepare mono"); test::ok (mc.numChannels() == 1, "1 channel");
        test::ok (mc.setIr (irM.data(), len), "setIr (mono convenience)");
        std::vector<float> y (n); for (int i = 0; i < n; ++i) y[(std::size_t) i] = xL[(std::size_t) i];
        for (int o = 0; o < n; o += 512) { float* io[1] { y.data() + o }; mc.process (io, io, 1, std::min (512, n - o)); }
        const std::vector<float> ref = convRef (irM, xL, P, irMax);
        test::ok (maxDiff (y, ref, settled, n) < 1e-4, "mono output == irM ∗ xL");
    }

    // --- no allocation in process(), even across a topology-changing swap ---
    test::group ("MatrixConvolver no-alloc in process() (incl. across a swap)");
    {
        MC mc; mc.prepare (P, irMax, xfade, 2);
        { const float* b[2] { irM.data(), irS.data() }; mc.setOperator (MC::Topology::MSDiag, b, 2, len); }
        std::vector<float> l (512, 0.2f), rr (512, -0.1f); float* io[2] { l.data(), rr.data() };
        mc.process (io, io, 2, 512);                      // consume the initial fade-in
        { const float* b[4] { irLL.data(), irLR.data(), irRL.data(), irRR.data() }; mc.setOperator (MC::Topology::Full, b, 4, len); }
        const long before = g_allocs.load();
        mc.process (io, io, 2, 512);                      // crosses the crossfade (MSDiag→Full)
        mc.process (io, io, 2, 512);
        test::okNoAlloc (g_allocs.load() == before, "process() performed zero heap allocations across a topology swap");
    }

    // --- swap coalescing + unprepared guards ---
    test::group ("MatrixConvolver swap coalescing + unprepared guards");
    {
        MC fresh;
        const float* b[2] { irM.data(), irS.data() };
        test::ok (! fresh.setOperator (MC::Topology::MSDiag, b, 2, len), "setOperator before prepare() rejected");
        float l[16] {}, rr[16] {}; float* io[2] { l, rr }; fresh.process (io, io, 2, 16);   // no-op, no crash
        test::ok (fresh.prepare (P, irMax, xfade, 2), "prepare after the rejected call");
        test::ok (fresh.setOperator (MC::Topology::MSDiag, b, 2, len), "first operator accepted (idle)");
        std::vector<float> l2 (512, 0.1f), r2 (512, 0.1f); float* io2[2] { l2.data(), r2.data() };
        fresh.process (io2, io2, 2, 8);                  // begin the (cold) fade → busy
        test::ok (fresh.isBusy(), "busy during the cold prime");
        const float* b2[4] { irLL.data(), irLR.data(), irRL.data(), irRR.data() };
        test::ok (! fresh.setOperator (MC::Topology::Full, b2, 4, len), "second operator rejected while busy (host coalesces)");
    }

    return test::report();
}

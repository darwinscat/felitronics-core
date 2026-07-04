// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// JUCE-free self-tests for MatrixConvolverNupc — the non-uniform (Gardner) sibling of MatrixConvolver, Phase 2
// (mono + LRDiag, instant swap). Reference-NULL against the proven PartitionedConvolver FROM SAMPLE 0 (unlike
// MatrixConvolver's own tests, which skip the cold crossfade — this convolver is true-conv from sample 0), plus
// a direct time-domain check on short IRs. Proves: mono == one convolution; LRDiag == two INDEPENDENT
// convolutions (yL free of xR, yR free of xL); channel isolation; true zero latency (impulse in → IR out, no
// shift); an instant operator swap that honours the warm history; MSDiag/Full correctly rejected (Phase 3);
// no allocation in process(); the MatrixConvolver-compatible API surface + prepare guards.

#include <felitronics_test.h>
#include <felitronics/convolution/MatrixConvolverNupc.h>
#include <felitronics/convolution/PartitionedConvolver.h>

#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <new>
#if defined(_WIN32)
 #include <malloc.h>
#endif
#include <vector>

static std::atomic<long> g_allocs { 0 };
void* operator new      (std::size_t s) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s ? s : 1); }
void* operator new[]    (std::size_t s) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s ? s : 1); }
void  operator delete   (void* p) noexcept { std::free (p); }
void  operator delete[] (void* p) noexcept { std::free (p); }
void  operator delete   (void* p, std::size_t) noexcept { std::free (p); }
void  operator delete[] (void* p, std::size_t) noexcept { std::free (p); }
static inline void* countedAlignedNew (std::size_t s, std::align_val_t a)
{
    g_allocs.fetch_add (1, std::memory_order_relaxed);
    const std::size_t al = (std::size_t) a < sizeof (void*) ? sizeof (void*) : (std::size_t) a;
   #if defined(_WIN32)
    void* p = _aligned_malloc (s ? s : 1, al);
   #else
    void* p = nullptr; if (::posix_memalign (&p, al, s ? s : 1) != 0) p = nullptr;
   #endif
    if (p == nullptr) throw std::bad_alloc();
    return p;
}
static inline void countedAlignedFree (void* p) noexcept
{
   #if defined(_WIN32)
    _aligned_free (p);
   #else
    std::free (p);
   #endif
}
void* operator new      (std::size_t s, std::align_val_t a) { return countedAlignedNew (s, a); }
void* operator new[]    (std::size_t s, std::align_val_t a) { return countedAlignedNew (s, a); }
void  operator delete   (void* p, std::align_val_t) noexcept { countedAlignedFree (p); }
void  operator delete[] (void* p, std::align_val_t) noexcept { countedAlignedFree (p); }
void  operator delete   (void* p, std::size_t, std::align_val_t) noexcept { countedAlignedFree (p); }
void  operator delete[] (void* p, std::size_t, std::align_val_t) noexcept { countedAlignedFree (p); }

using namespace felitronics;
using MCN = convolution::MatrixConvolverNupc<>;

struct Lcg { unsigned long long s; float next() { s = s * 6364136223846793005ULL + 1442695040888963407ULL; return (float) ((s >> 40) & 0xffff) / 32768.0f - 1.0f; } };

// exact reference via the proven uniform PartitionedConvolver (any P computes the same linear convolution).
static std::vector<float> convRef (const std::vector<float>& x, const std::vector<float>& h, int maxIr)
{
    convolution::PartitionedConvolver<> pc; pc.prepare (128, maxIr); pc.setIr (h.data(), (int) h.size());
    std::vector<float> y (x.size(), 0.0f); pc.process (x.data(), y.data(), (int) x.size());
    return y;
}
static double maxDiff (const std::vector<float>& a, const std::vector<float>& b, int from, int to)
{ double m = 0.0; for (int i = from; i < to; ++i) m = std::max (m, (double) std::fabs (a[(std::size_t) i] - b[(std::size_t) i])); return m; }
static double peak (const std::vector<float>& a) { double m = 0.0; for (float v : a) m = std::max (m, (double) std::fabs (v)); return m; }

static void runMono (MCN& mc, const std::vector<float>& x, std::vector<float>& y, int block)
{
    y.assign (x.size(), 0.0f);
    for (int i = 0; i < (int) x.size(); i += block) { const int m = std::min (block, (int) x.size() - i); const float* in[1] { &x[(std::size_t) i] }; float* out[1] { &y[(std::size_t) i] }; mc.process (in, out, 1, m); }
}
static void runStereo (MCN& mc, const std::vector<float>& xL, const std::vector<float>& xR, std::vector<float>& yL, std::vector<float>& yR, int block)
{
    yL.assign (xL.size(), 0.0f); yR.assign (xR.size(), 0.0f);
    for (int i = 0; i < (int) xL.size(); i += block) { const int m = std::min (block, (int) xL.size() - i); const float* in[2] { &xL[(std::size_t) i], &xR[(std::size_t) i] }; float* out[2] { &yL[(std::size_t) i], &yR[(std::size_t) i] }; mc.process (in, out, 2, m); }
}

int main()
{
    std::printf ("felitronics::convolution MatrixConvolverNupc tests (Phase 2: mono + LRDiag)\n");
    const int maxIr = 131072;

    // --- mono == a single PartitionedConvolver, from sample 0, across blocks ---
    test::group ("mono == one convolution (NULL vs PartitionedConvolver, from sample 0)");
    {
        for (int L : { 60, 400, 5000, 131072 })
            for (int block : { 1, 64, 257, 4096 })
            {
                Lcg r { 11u + (unsigned) (L * 7 + block) };
                std::vector<float> h ((std::size_t) L); for (auto& v : h) v = 0.05f * r.next();
                const int N = L + 3000; std::vector<float> x ((std::size_t) N); for (auto& v : x) v = 0.3f * r.next();
                for (int p : { 0, 127, 128, 511, 512, L - 1 }) if (p >= 0 && p < N) x[(std::size_t) p] += 1.0f;
                MCN mc; test::ok (mc.prepare (128, maxIr, 128, 1), "prepare mono");
                mc.setIr (h.data(), L);
                std::vector<float> y; runMono (mc, x, y, block);
                const double rel = maxDiff (y, convRef (x, h, maxIr), 0, N) / (peak (convRef (x, h, maxIr)) + 1e-30);
                if (! (rel < 2e-4)) std::printf ("      mono L=%d block=%d rel=%.3e\n", L, block, rel);
                test::ok (rel < 2e-4, "mono NUPC == convolution L=" + std::to_string (L) + " block=" + std::to_string (block));
            }
    }

    // --- LRDiag == two INDEPENDENT convolutions (yL free of xR, yR free of xL) ---
    test::group ("LRDiag == two independent convolutions");
    {
        const int L = 4000, N = L + 3000;
        Lcg r { 4242 };
        std::vector<float> hL ((std::size_t) L), hR ((std::size_t) L); for (auto& v : hL) v = 0.05f * r.next(); for (auto& v : hR) v = 0.05f * r.next();
        std::vector<float> xL ((std::size_t) N), xR ((std::size_t) N); for (auto& v : xL) v = 0.3f * r.next(); for (auto& v : xR) v = 0.3f * r.next();
        MCN mc; test::ok (mc.prepare (128, maxIr, 128, 2) && mc.numChannels() == 2, "prepare stereo");
        const float* banks[2] { hL.data(), hR.data() };
        test::ok (mc.setOperator (MCN::Topology::LRDiag, banks, 2, L), "setOperator LRDiag");
        std::vector<float> yL, yR; runStereo (mc, xL, xR, yL, yR, 512);
        const std::vector<float> refL = convRef (xL, hL, maxIr), refR = convRef (xR, hR, maxIr);
        test::ok (maxDiff (yL, refL, 0, N) / (peak (refL) + 1e-30) < 2e-4, "yL == hL ∗ xL (independent of xR)");
        test::ok (maxDiff (yR, refR, 0, N) / (peak (refR) + 1e-30) < 2e-4, "yR == hR ∗ xR (independent of xL)");

        // yL must not change when xR changes (true channel independence)
        std::vector<float> xR2 ((std::size_t) N); for (auto& v : xR2) v = 0.7f * r.next();
        MCN mc2; mc2.prepare (128, maxIr, 128, 2); mc2.setOperator (MCN::Topology::LRDiag, banks, 2, L);
        std::vector<float> yL2, yR2; runStereo (mc2, xL, xR2, yL2, yR2, 512);
        test::ok (maxDiff (yL, yL2, 0, N) < 1e-6, "yL unchanged when xR changes (no cross-channel bleed)");
    }

    // --- channel isolation: an impulse in L only leaves R silent ---
    test::group ("channel isolation: impulse in L → R silent");
    {
        const int L = 2000, N = 5000;
        Lcg r { 5 };
        std::vector<float> hL ((std::size_t) L), hR ((std::size_t) L); for (auto& v : hL) v = 0.05f * r.next(); for (auto& v : hR) v = 0.05f * r.next();
        std::vector<float> xL ((std::size_t) N, 0.0f), xR ((std::size_t) N, 0.0f); xL[100] = 1.0f;
        MCN mc; mc.prepare (128, maxIr, 128, 2);
        const float* banks[2] { hL.data(), hR.data() }; mc.setOperator (MCN::Topology::LRDiag, banks, 2, L);
        std::vector<float> yL, yR; runStereo (mc, xL, xR, yL, yR, 128);
        test::ok (peak (yL) > 1e-3, "yL responds to the L impulse");
        test::ok (peak (yR) < 1e-9, "yR is exactly silent (R input all zero → no bleed)");
        // and yL[100+i] == hL[i]  (zero-latency, correctly placed)
        double e = 0.0; for (int i = 0; i < L; ++i) e = std::max (e, (double) std::fabs (yL[(std::size_t) (100 + i)] - hL[(std::size_t) i]));
        test::ok (e < 1e-4, "yL impulse response == hL, no shift (true zero latency)");
    }

    // --- zero-latency: impulse in → IR out with no shift (mono) ---
    test::group ("TRUE sample-zero-latency (mono impulse == IR)");
    {
        const int L = 3000;
        Lcg r { 909 };
        std::vector<float> h ((std::size_t) L); for (auto& v : h) v = 0.05f * r.next();
        std::vector<float> imp ((std::size_t) (L + 200), 0.0f); imp[0] = 1.0f;
        MCN mc; mc.prepare (128, maxIr, 128, 1); mc.setIr (h.data(), L);
        std::vector<float> y; runMono (mc, imp, y, 100);
        double e = 0.0; for (int i = 0; i < L; ++i) e = std::max (e, (double) std::fabs (y[(std::size_t) i] - h[(std::size_t) i]));
        std::printf ("      impulse@0 max|y-h| = %.3e (y[0]=%.6f h[0]=%.6f)\n", e, y[0], h[0]);
        test::ok (e < 1e-4, "y == h, y[0]==h[0] → zero added latency");
    }

    // --- instant operator swap honours the warm history (mono) ---
    test::group ("instant operator swap (warm history)");
    {
        const int L = 3000, N = 9000, T = 4000;
        Lcg r { 3131 };
        std::vector<float> h1 ((std::size_t) L), h2 ((std::size_t) L); for (auto& v : h1) v = 0.05f * r.next(); for (auto& v : h2) v = 0.05f * r.next();
        std::vector<float> x ((std::size_t) N); for (auto& v : x) v = 0.3f * r.next();
        MCN mc; mc.prepare (128, maxIr, 128, 1); mc.setIr (h1.data(), L);
        std::vector<float> y ((std::size_t) N, 0.0f);
        bool swapped = false;
        for (int i = 0; i < N; i += 256)
        {
            if (! swapped && i >= T) { if (mc.setIr (h2.data(), L)) swapped = true; }
            const int m = std::min (256, N - i); const float* in[1] { &x[(std::size_t) i] }; float* out[1] { &y[(std::size_t) i] }; mc.process (in, out, 1, m);
        }
        test::ok (swapped, "swap accepted mid-stream");
        // after the swap settles (past one big-stage fill), output == h2 convolved over the WHOLE warm stream
        const std::vector<float> ref2 = convRef (x, h2, maxIr);
        const int settled = T + 4096 + 256;
        test::ok (maxDiff (y, ref2, settled, N) / (peak (ref2) + 1e-30) < 2e-4, "post-swap output == h2 ∗ x with WARM history (not cold-restarted)");
    }

    // --- instant swap at an ARBITRARY (non-boundary) sample: every stage's tail must be primed correctly for
    //     its in-progress chunk, even when stages are mid-chunk at different phases. An instant warm swap makes
    //     the output == h2 ∗ x EXACTLY from the swap sample (no transient), so we NULL the WHOLE post-swap region
    //     — a mid-chunk prime bug can't hide in a settle-skip. L>4096 so the 2048 tail stage has >1 active
    //     partition (exercises the multi-partition FDL wrap through the swap). ---
    test::group ("instant swap mid-chunk (mixed phases, multi-partition tail) — exact from the swap sample");
    {
        const int L = 6000, N = 12000, T = 4001;
        Lcg r { 2024 };
        std::vector<float> h1 ((std::size_t) L), h2 ((std::size_t) L); for (auto& v : h1) v = 0.05f * r.next(); for (auto& v : h2) v = 0.05f * r.next();
        std::vector<float> x ((std::size_t) N); for (auto& v : x) v = 0.3f * r.next();
        MCN mc; mc.prepare (128, maxIr, 128, 1); mc.setIr (h1.data(), L);
        std::vector<float> y ((std::size_t) N, 0.0f);
        int swapAt = -1;
        for (int i = 0; i < N; i += 137)   // 137 coprime with the stage sizes → the swap lands off every boundary
        {
            if (swapAt < 0 && i >= T && mc.setIr (h2.data(), L)) swapAt = i;
            const int m = std::min (137, N - i); const float* in[1] { &x[(std::size_t) i] }; float* out[1] { &y[(std::size_t) i] }; mc.process (in, out, 1, m);
        }
        test::ok (swapAt > 0 && (swapAt % 2048) != 0 && (swapAt % 1024) != 0, "swap landed mid-chunk (off the big-stage boundaries)");
        const std::vector<float> ref2 = convRef (x, h2, maxIr);
        test::ok (maxDiff (y, ref2, swapAt, N) / (peak (ref2) + 1e-30) < 2e-4, "post-swap output == h2 ∗ x EXACTLY from the swap sample (mid-chunk tails primed right)");
    }

    // --- stereo LRDiag→LRDiag warm swap with a host block LARGER than B_max ---
    test::group ("stereo LRDiag warm swap, host block > B_max");
    {
        const int L = 6000, N = 24576, T = 8192;
        Lcg r { 71 };
        std::vector<float> aL ((std::size_t) L), aR ((std::size_t) L), bL ((std::size_t) L), bR ((std::size_t) L);
        for (auto& v : aL) v = 0.05f * r.next(); for (auto& v : aR) v = 0.05f * r.next(); for (auto& v : bL) v = 0.05f * r.next(); for (auto& v : bR) v = 0.05f * r.next();
        std::vector<float> xL ((std::size_t) N), xR ((std::size_t) N); for (auto& v : xL) v = 0.3f * r.next(); for (auto& v : xR) v = 0.3f * r.next();
        MCN mc; mc.prepare (128, maxIr, 128, 2);
        const float* op1[2] { aL.data(), aR.data() }; mc.setOperator (MCN::Topology::LRDiag, op1, 2, L);
        std::vector<float> yL ((std::size_t) N, 0.0f), yR ((std::size_t) N, 0.0f);
        int swapAt = -1;
        for (int i = 0; i < N; i += 8192)   // block 8192 > B_max=2048 (the per-sample loop absorbs any block)
        {
            if (swapAt < 0 && i >= T) { const float* op2[2] { bL.data(), bR.data() }; if (mc.setOperator (MCN::Topology::LRDiag, op2, 2, L)) swapAt = i; }
            const int m = std::min (8192, N - i); const float* in[2] { &xL[(std::size_t) i], &xR[(std::size_t) i] }; float* out[2] { &yL[(std::size_t) i], &yR[(std::size_t) i] }; mc.process (in, out, 2, m);
        }
        test::ok (swapAt > 0, "stereo swap accepted with block 8192 > B_max");
        const std::vector<float> rL = convRef (xL, bL, maxIr), rR = convRef (xR, bR, maxIr);
        test::ok (maxDiff (yL, rL, swapAt, N) / (peak (rL) + 1e-30) < 2e-4, "post-swap yL == bL ∗ xL (warm, block>B_max)");
        test::ok (maxDiff (yR, rR, swapAt, N) / (peak (rR) + 1e-30) < 2e-4, "post-swap yR == bR ∗ xR (warm, block>B_max)");
    }

    // --- MSDiag / Full are rejected in Phase 2 ---
    test::group ("MSDiag / Full rejected in Phase 2 (Phase-3 topologies)");
    {
        const int L = 200; std::vector<float> a ((std::size_t) L, 0.01f), b ((std::size_t) L, 0.02f), c ((std::size_t) L, 0.03f), d ((std::size_t) L, 0.04f);
        MCN mc; mc.prepare (128, maxIr, 128, 2);
        const float* two[2] { a.data(), b.data() }; const float* four[4] { a.data(), b.data(), c.data(), d.data() };
        test::ok (! mc.setOperator (MCN::Topology::MSDiag, two, 2, L), "MSDiag rejected (Phase 3)");
        test::ok (! mc.setOperator (MCN::Topology::Full, four, 4, L), "Full rejected (Phase 3)");
        test::ok (  mc.setOperator (MCN::Topology::LRDiag, two, 2, L), "LRDiag still accepted");
    }

    // --- no allocation in process(), including across an instant swap ---
    test::group ("no allocation in process() (incl. across a swap)");
    {
        const int L = 20000;
        std::vector<float> h1 ((std::size_t) L, 0.001f), h2 ((std::size_t) L, -0.001f);
        MCN mc; mc.prepare (128, maxIr, 128, 1); mc.setIr (h1.data(), L);
        std::vector<float> x (2048, 0.2f), y (2048, 0.0f);
        const float* in[1] { x.data() }; float* out[1] { y.data() };
        mc.process (in, out, 1, 2048);
        mc.setIr (h2.data(), L);                                     // stage a swap
        const long before = g_allocs.load();
        mc.process (in, out, 1, 2048);                              // consumes the instant swap
        mc.process (in, out, 1, 2048);
        test::okNoAlloc (g_allocs.load() == before, "process() performed zero heap allocations across a swap");
    }

    // --- in-place aliasing (in == out) ---
    test::group ("in-place aliasing (in == out) matches out-of-place (LRDiag)");
    {
        const int L = 3000, N = L + 2000;
        Lcg r { 88 };
        std::vector<float> hL ((std::size_t) L), hR ((std::size_t) L); for (auto& v : hL) v = 0.05f * r.next(); for (auto& v : hR) v = 0.05f * r.next();
        std::vector<float> xL ((std::size_t) N), xR ((std::size_t) N); for (auto& v : xL) v = 0.3f * r.next(); for (auto& v : xR) v = 0.3f * r.next();
        MCN a; a.prepare (128, maxIr, 128, 2); const float* bk[2] { hL.data(), hR.data() }; a.setOperator (MCN::Topology::LRDiag, bk, 2, L);
        std::vector<float> yL, yR; runStereo (a, xL, xR, yL, yR, 128);
        MCN b; b.prepare (128, maxIr, 128, 2); b.setOperator (MCN::Topology::LRDiag, bk, 2, L);
        std::vector<float> iL = xL, iR = xR;
        for (int i = 0; i < N; i += 128) { const int m = std::min (128, N - i); const float* in[2] { &iL[(std::size_t) i], &iR[(std::size_t) i] }; float* out[2] { &iL[(std::size_t) i], &iR[(std::size_t) i] }; b.process (in, out, 2, m); }
        test::ok (maxDiff (iL, yL, 0, N) < 1e-6 && maxDiff (iR, yR, 0, N) < 1e-6, "in-place bit-matches out-of-place");
    }

    // --- API surface + prepare guards ---
    test::group ("API parity + prepare guards");
    {
        test::ok (MCN::latencySamples() == 0, "latencySamples()==0");
        test::ok (MCN::numBanksFor (MCN::Topology::Full) == 4 && MCN::numBanksFor (MCN::Topology::LRDiag) == 2, "numBanksFor parity");
        MCN mc;
        const float* b[2] { nullptr, nullptr };
        test::ok (! mc.setOperator (MCN::Topology::LRDiag, b, 2, 100), "setOperator before prepare rejected");
        test::ok (! mc.prepare (127, maxIr, 128, 2), "non-pow2 head rejected");
        test::ok (! mc.prepare (128, maxIr, 128, 3), "numChannels=3 rejected");
        test::ok (  mc.prepare (128, maxIr, 128, 2), "prepare after rejected calls");
        test::ok (! mc.isBusy(), "idle after prepare");
        float in0[16] {}, in1[16] {}; const float* in[2] { in0, in1 }; float* out[2] { in0, in1 };
        mc.process (in, out, 2, 16);   // no operator yet → no crash (outputs from empty slot)
        test::ok (true, "process without an operator did not crash");
    }

    return test::report();
}

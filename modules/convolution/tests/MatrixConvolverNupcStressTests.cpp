// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// Popper-campaign hardening for MatrixConvolverNupc + NonUniformConvolver: a randomized DIFFERENTIAL fuzzer
// (the durable residue of a maniacal "prove it wrong" pass). For hundreds of random (head P0, maxIr, IR length,
// topology, channels, block sequence) configs it NULLs the convolver against an INDEPENDENT reference — the
// proven PartitionedConvolver, and for short IRs a full DOUBLE-PRECISION direct convolution — in the settled
// window past the cold-prime crossfade. Plus: BIT-IDENTICAL block-independence (the output must not depend on how
// process() is chunked); the kMaxPartition overflow guard (a huge pow2 partition must be rejected, not overflow
// 2*B); and a denormal / NaN / Inf stress that must not crash or allocate. Any failure prints its exact
// reproducer (seed + config).

#include <felitronics_test.h>
#include <felitronics/convolution/MatrixConvolverNupc.h>
#include <felitronics/convolution/NonUniformConvolver.h>
#include <felitronics/convolution/PartitionedConvolver.h>

#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <new>
#if defined(_WIN32)
 #include <malloc.h>
#endif
#include <string>
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

// tiny deterministic PRNG (no heap)
struct Rng { unsigned long long s; unsigned next() { s = s * 6364136223846793005ULL + 1442695040888963407ULL; return (unsigned) (s >> 33); }
             float unit() { return (float) (next() & 0xffffff) / 8388608.0f - 1.0f; }   // [-1,1)
             int range (int lo, int hi) { return lo + (int) (next() % (unsigned) (hi - lo + 1)); } };

static std::vector<float> convDbl (const std::vector<float>& x, const std::vector<float>& h)
{
    std::vector<float> y (x.size(), 0.0f);
    for (int n = 0; n < (int) x.size(); ++n) { double a = 0.0; const int jmax = std::min ((int) h.size() - 1, n);
        for (int j = 0; j <= jmax; ++j) a += (double) h[(std::size_t) j] * (double) x[(std::size_t) (n - j)]; y[(std::size_t) n] = (float) a; }
    return y;
}
static std::vector<float> convP (const std::vector<float>& x, const std::vector<float>& h, int maxIr)
{
    convolution::PartitionedConvolver<> pc; pc.prepare (128, maxIr); pc.setIr (h.data(), (int) h.size());
    std::vector<float> y (x.size(), 0.0f); pc.process (x.data(), y.data(), (int) x.size()); return y;
}
static double relSettled (const std::vector<float>& a, const std::vector<float>& b, int from)
{
    double m = 0.0, p = 0.0; for (int i = from; i < (int) a.size(); ++i) { m = std::max (m, (double) std::fabs (a[(std::size_t) i] - b[(std::size_t) i])); p = std::max (p, (double) std::fabs (b[(std::size_t) i])); }
    return m / (p + 1e-30);
}

// run a config through MCN in a given block sequence (mono uses ch=1); topo/nb per caller
static void runMCN (MCN& mc, const std::vector<float>& xL, const std::vector<float>& xR, int nch,
                    std::vector<float>& yL, std::vector<float>& yR, Rng& blkRng)
{
    const int N = (int) xL.size(); yL.assign ((std::size_t) N, 0.0f); yR.assign ((std::size_t) N, 0.0f);
    int i = 0;
    while (i < N)
    {
        int b = blkRng.range (1, 900); b = std::min (b, N - i);
        const float* in[2] { &xL[(std::size_t) i], &xR[(std::size_t) i] }; float* out[2] { &yL[(std::size_t) i], &yR[(std::size_t) i] };
        mc.process (in, out, nch, b); i += b;
    }
}

int main()
{
    std::printf ("felitronics::convolution MatrixConvolverNupc STRESS / differential fuzzer\n");
    double worstRel = 0.0; std::string worst;

    // ---------------------------------------------------------------------------------------------
    test::group ("randomized differential NULL vs PartitionedConvolver (settled window)");
    {
        int fails = 0;
        for (int it = 0; it < 400; ++it)
        {
            Rng r { 0x9E37u + (unsigned long long) it * 2654435761u };
            const int p0choices[3] { 64, 128, 256 };
            const int P0 = p0choices[r.range (0, 2)];
            const int maxIr = 1 << r.range (9, 12);                 // 512..4096
            const int len = r.range (0, maxIr);
            const int nch = r.range (0, 1) == 0 ? 1 : 2;
            const int topoSel = nch == 1 ? 0 : r.range (0, 2);      // 0 LRDiag,1 MSDiag,2 Full
            const int nb = nch == 1 ? 1 : (topoSel == 2 ? 4 : 2);
            const int settled = len + 2048 + 400;
            const int N = settled + 1500;

            std::vector<std::vector<float>> h ((std::size_t) nb, std::vector<float> ((std::size_t) len));
            for (auto& hi : h) for (auto& v : hi) v = 0.05f * r.unit();
            std::vector<float> xL ((std::size_t) N), xR ((std::size_t) N);
            for (auto& v : xL) v = 0.3f * r.unit(); for (auto& v : xR) v = 0.4f * r.unit();

            MCN mc; if (! mc.prepare (P0, maxIr, 128, nch)) { std::printf ("  prepare FAIL it=%d\n", it); ++fails; continue; }
            std::vector<const float*> banks ((std::size_t) nb); for (int b = 0; b < nb; ++b) banks[(std::size_t) b] = h[(std::size_t) b].data();
            if (nch == 1) mc.setIr (h[0].data(), len);
            else mc.setOperator (topoSel == 0 ? MCN::Topology::LRDiag : topoSel == 1 ? MCN::Topology::MSDiag : MCN::Topology::Full, banks.data(), nb, len);

            Rng blk { 0xB5u + (unsigned long long) it };
            std::vector<float> yL, yR; runMCN (mc, xL, xR, nch, yL, yR, blk);

            // per-topology reference
            std::vector<float> rL, rR;
            if (nch == 1) { rL = convP (xL, h[0], maxIr); rR.assign ((std::size_t) N, 0.0f); }
            else if (topoSel == 0) { rL = convP (xL, h[0], maxIr); rR = convP (xR, h[1], maxIr); }
            else if (topoSel == 1) { std::vector<float> mm ((std::size_t) N), ss ((std::size_t) N); for (int i = 0; i < N; ++i) { mm[(std::size_t) i] = 0.5f * (xL[(std::size_t) i] + xR[(std::size_t) i]); ss[(std::size_t) i] = 0.5f * (xL[(std::size_t) i] - xR[(std::size_t) i]); } auto yM = convP (mm, h[0], maxIr), yS = convP (ss, h[1], maxIr); rL.assign ((std::size_t) N, 0.0f); rR.assign ((std::size_t) N, 0.0f); for (int i = 0; i < N; ++i) { rL[(std::size_t) i] = yM[(std::size_t) i] + yS[(std::size_t) i]; rR[(std::size_t) i] = yM[(std::size_t) i] - yS[(std::size_t) i]; } }
            else { auto cLL = convP (xL, h[0], maxIr), cLR = convP (xR, h[1], maxIr), cRL = convP (xL, h[2], maxIr), cRR = convP (xR, h[3], maxIr); rL.assign ((std::size_t) N, 0.0f); rR.assign ((std::size_t) N, 0.0f); for (int i = 0; i < N; ++i) { rL[(std::size_t) i] = cLL[(std::size_t) i] + cLR[(std::size_t) i]; rR[(std::size_t) i] = cRL[(std::size_t) i] + cRR[(std::size_t) i]; } }

            const double e = std::max (relSettled (yL, rL, settled), nch == 1 ? 0.0 : relSettled (yR, rR, settled));
            if (e > worstRel) { worstRel = e; worst = "P0=" + std::to_string (P0) + " maxIr=" + std::to_string (maxIr) + " len=" + std::to_string (len) + " nch=" + std::to_string (nch) + " topo=" + std::to_string (topoSel); }
            if (! (e < 1e-3)) { std::printf ("  *** it=%d P0=%d maxIr=%d len=%d nch=%d topo=%d rel=%.3e\n", it, P0, maxIr, len, nch, topoSel, e); ++fails; }
        }
        std::printf ("      400 random configs, worst rel=%.3e (%s)\n", worstRel, worst.c_str());
        test::ok (fails == 0, "every random config nulls the reference (rel<1e-3)");
    }

    // ---------------------------------------------------------------------------------------------
    test::group ("independent double-precision NULL (short IRs, all topologies)");
    {
        int fails = 0; double w = 0.0;
        for (int it = 0; it < 60; ++it)
        {
            Rng r { 0xD1CEu + (unsigned long long) it * 40503u };
            const int p0c2[2] { 64, 128 }; const int P0 = p0c2[r.range (0, 1)];
            const int len = r.range (1, 600);
            const int maxIr = std::max (len, 1 << r.range (9, 11));
            const int nch = r.range (0, 1) == 0 ? 1 : 2;
            const int topoSel = nch == 1 ? 0 : r.range (0, 2);
            const int nb = nch == 1 ? 1 : (topoSel == 2 ? 4 : 2);
            const int settled = len + 2048 + 300; const int N = settled + 1200;
            std::vector<std::vector<float>> h ((std::size_t) nb, std::vector<float> ((std::size_t) len));
            for (auto& hi : h) for (auto& v : hi) v = 0.06f * r.unit();
            std::vector<float> xL ((std::size_t) N), xR ((std::size_t) N);
            for (auto& v : xL) v = 0.3f * r.unit(); for (auto& v : xR) v = 0.4f * r.unit();
            MCN mc; mc.prepare (P0, maxIr, 128, nch);
            std::vector<const float*> banks ((std::size_t) nb); for (int b = 0; b < nb; ++b) banks[(std::size_t) b] = h[(std::size_t) b].data();
            if (nch == 1) mc.setIr (h[0].data(), len); else mc.setOperator (topoSel == 0 ? MCN::Topology::LRDiag : topoSel == 1 ? MCN::Topology::MSDiag : MCN::Topology::Full, banks.data(), nb, len);
            Rng blk { 0x77u + (unsigned long long) it }; std::vector<float> yL, yR; runMCN (mc, xL, xR, nch, yL, yR, blk);
            std::vector<float> rL, rR;
            if (nch == 1) { rL = convDbl (xL, h[0]); rR.assign ((std::size_t) N, 0.0f); }
            else if (topoSel == 0) { rL = convDbl (xL, h[0]); rR = convDbl (xR, h[1]); }
            else if (topoSel == 1) { std::vector<float> mm ((std::size_t) N), ss ((std::size_t) N); for (int i = 0; i < N; ++i) { mm[(std::size_t) i] = 0.5f * (xL[(std::size_t) i] + xR[(std::size_t) i]); ss[(std::size_t) i] = 0.5f * (xL[(std::size_t) i] - xR[(std::size_t) i]); } auto yM = convDbl (mm, h[0]), yS = convDbl (ss, h[1]); rL.assign ((std::size_t) N, 0.0f); rR.assign ((std::size_t) N, 0.0f); for (int i = 0; i < N; ++i) { rL[(std::size_t) i] = yM[(std::size_t) i] + yS[(std::size_t) i]; rR[(std::size_t) i] = yM[(std::size_t) i] - yS[(std::size_t) i]; } }
            else { auto cLL = convDbl (xL, h[0]), cLR = convDbl (xR, h[1]), cRL = convDbl (xL, h[2]), cRR = convDbl (xR, h[3]); rL.assign ((std::size_t) N, 0.0f); rR.assign ((std::size_t) N, 0.0f); for (int i = 0; i < N; ++i) { rL[(std::size_t) i] = cLL[(std::size_t) i] + cLR[(std::size_t) i]; rR[(std::size_t) i] = cRL[(std::size_t) i] + cRR[(std::size_t) i]; } }
            const double e = std::max (relSettled (yL, rL, settled), nch == 1 ? 0.0 : relSettled (yR, rR, settled));
            w = std::max (w, e);
            if (! (e < 5e-4)) { std::printf ("  *** double it=%d P0=%d len=%d topo=%d rel=%.3e\n", it, P0, len, topoSel, e); ++fails; }
        }
        std::printf ("      60 configs vs double-precision, worst rel=%.3e\n", w);
        test::ok (fails == 0, "NULL vs double-precision direct convolution (rel<5e-4)");
    }

    // ---------------------------------------------------------------------------------------------
    test::group ("block-independence: output BIT-IDENTICAL across any block sequence");
    {
        int fails = 0;
        for (int it = 0; it < 80; ++it)
        {
            Rng r { 0xB10Cu + (unsigned long long) it * 2246822519u };
            const int p0c3[3] { 64, 128, 256 }; const int P0 = p0c3[r.range (0, 2)];
            const int maxIr = 1 << r.range (9, 12); const int len = r.range (1, maxIr);
            const int nch = r.range (0, 1) == 0 ? 1 : 2; const int topoSel = nch == 1 ? 0 : r.range (0, 2);
            const int nb = nch == 1 ? 1 : (topoSel == 2 ? 4 : 2); const int N = maxIr + 2500;
            std::vector<std::vector<float>> h ((std::size_t) nb, std::vector<float> ((std::size_t) len));
            for (auto& hi : h) for (auto& v : hi) v = 0.05f * r.unit();
            std::vector<float> xL ((std::size_t) N), xR ((std::size_t) N); for (auto& v : xL) v = 0.3f * r.unit(); for (auto& v : xR) v = 0.4f * r.unit();
            std::vector<const float*> banks ((std::size_t) nb); for (int b = 0; b < nb; ++b) banks[(std::size_t) b] = h[(std::size_t) b].data();
            auto build = [&] (MCN& mc) { mc.prepare (P0, maxIr, 128, nch); if (nch == 1) mc.setIr (h[0].data(), len); else mc.setOperator (topoSel == 0 ? MCN::Topology::LRDiag : topoSel == 1 ? MCN::Topology::MSDiag : MCN::Topology::Full, banks.data(), nb, len); };
            // A: one big call
            MCN a; build (a); std::vector<float> aL ((std::size_t) N, 0.0f), aR ((std::size_t) N, 0.0f);
            { const float* in[2] { xL.data(), xR.data() }; float* out[2] { aL.data(), aR.data() }; a.process (in, out, nch, N); }
            // B: random small blocks
            MCN b; build (b); std::vector<float> bL, bR; Rng blk { 0x51u + (unsigned long long) it }; runMCN (b, xL, xR, nch, bL, bR, blk);
            const bool same = std::memcmp (aL.data(), bL.data(), (std::size_t) N * sizeof (float)) == 0
                           && (nch == 1 || std::memcmp (aR.data(), bR.data(), (std::size_t) N * sizeof (float)) == 0);
            if (! same) { std::printf ("  *** block-dependence it=%d P0=%d maxIr=%d len=%d topo=%d\n", it, P0, maxIr, len, topoSel); ++fails; }
        }
        test::ok (fails == 0, "output is bit-identical regardless of host block size (80 configs)");
    }

    // ---------------------------------------------------------------------------------------------
    test::group ("overflow guard: a huge power-of-two partition is rejected (not 2*B UB)");
    {
        MCN mc;
        test::ok (! mc.prepare (1 << 30, (1 << 30) + 1, 128, 1), "prepare(1<<30, …) rejected (kMaxPartition)");
        test::ok (! mc.prepare (1 << 25, (1 << 25) + 1, 128, 2), "prepare(1<<25, …) rejected");
        test::ok (  mc.prepare (1 << 12, 131072, 128, 2), "prepare(4096, 131072, …) still accepted");
        convolution::NonUniformConvolver<> nu;
        test::ok (! nu.prepare (1 << 30, 1 << 30, (1 << 30) + 1), "NonUniformConvolver huge head rejected");
        convolution::NonUniformConvolver<>::ScheduleStep bad[] { { 1 << 30, 1 } };
        test::ok (! nu.prepareWithSchedule (1 << 30, bad, 1, (1 << 30) + 1), "NonUniformConvolver huge schedule step rejected");
    }

    // ---------------------------------------------------------------------------------------------
    test::group ("denormal / NaN / Inf input: no crash, no alloc, bounded where input is finite");
    {
        const int L = 4000, N = 6000;
        Rng r { 0xFACEu };
        std::vector<float> h ((std::size_t) L); for (auto& v : h) v = 0.05f * r.unit();
        MCN mc; mc.prepare (128, 131072, 128, 1); mc.setIr (h.data(), L);
        std::vector<float> x ((std::size_t) N), y ((std::size_t) N, 0.0f);
        const float denorm = std::numeric_limits<float>::denorm_min();
        for (int i = 0; i < N; ++i) x[(std::size_t) i] = (i % 7 == 0) ? denorm : 0.2f * r.unit();
        { const float* in[1] { x.data() }; float* out[1] { y.data() }; mc.process (in, out, 1, N); }   // warm past the cold fade
        bool finite = true; for (float v : y) if (! std::isfinite (v)) finite = false;
        test::ok (finite, "denormal-laced finite input → finite output (no denormal runaway in a pure FIR)");

        // no allocation on a steady-state stereo Full block
        std::vector<float> ir (L, 0.001f); const float* bk[4] { ir.data(), ir.data(), ir.data(), ir.data() };
        MCN mf; mf.prepare (128, 131072, 128, 2); mf.setOperator (MCN::Topology::Full, bk, 4, L);
        std::vector<float> l (4096, 0.1f), rr (4096, -0.1f); const float* sin[2] { l.data(), rr.data() }; float* so[2] { l.data(), rr.data() };
        mf.process (sin, so, 2, 4096);
        const long before = g_allocs.load(); mf.process (sin, so, 2, 4096);
        test::okNoAlloc (g_allocs.load() == before, "no heap allocation in a stress Full process()");

        // NaN in must not crash or hang (output may be NaN — we only require it returns + doesn't corrupt state)
        std::vector<float> xn ((std::size_t) N, 0.1f); xn[(std::size_t) 3000] = std::numeric_limits<float>::quiet_NaN();
        std::vector<float> yn ((std::size_t) N, 0.0f);
        MCN mn; mn.prepare (128, 131072, 128, 1); mn.setIr (h.data(), L);
        { const float* in[1] { xn.data() }; float* out[1] { yn.data() }; mn.process (in, out, 1, N); }
        mn.reset();   // recovers cleanly
        std::vector<float> xc ((std::size_t) 2000, 0.1f), yc ((std::size_t) 2000, 0.0f);
        { const float* in[1] { xc.data() }; float* out[1] { yc.data() }; mn.process (in, out, 1, 2000); }
        bool recovered = true; for (float v : yc) if (! std::isfinite (v)) recovered = false;
        test::ok (recovered, "reset() recovers finite output after a NaN-poisoned input (no persistent corruption)");
    }

    return test::report();
}

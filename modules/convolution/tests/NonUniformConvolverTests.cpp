// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// JUCE-free self-tests for the NonUniformConvolver: the mono, true-sample-zero-latency, non-uniform (Gardner)
// partitioned convolver. The whole point of Phase 1 is a PROOF the non-uniform schedule computes the same
// convolution as the proven uniform primitives — so this is reference-NULL, hard:
//   • vs a DIRECT double-precision time-domain convolution (the independent analytic truth) for short IRs;
//   • vs the proven PartitionedConvolver (fast, same FFT family) for the full 131072-tap IR;
//   • TRUE zero latency proven by an impulse in → the IR out with NO shift (y[0]=h[0]);
//   • impulses placed on EVERY stage boundary (an off-by-one-block schedule bug shifts the whole tail there);
//   • multiple schedules (capped-doubling / a tuned {block,count} list / pure octave doubling / P0=64);
//   • arbitrary host blocks incl. block > maxBlock, block = 1, and a non-power-of-2 (257) block;
//   • in-place (in==out) aliasing; a short IR loaded into a large schedule (inactive/partial stages);
//   • and no allocation in process().

#include <felitronics_test.h>
#include <felitronics/convolution/NonUniformConvolver.h>
#include <felitronics/convolution/PartitionedConvolver.h>   // reference convolver for the long-IR null

#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <new>
#if defined(_WIN32)
 #include <malloc.h>   // _aligned_malloc / _aligned_free (MSVC has no posix_memalign)
#endif
#include <vector>

static std::atomic<long> g_allocs { 0 };
void* operator new      (std::size_t s) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s ? s : 1); }
void* operator new[]    (std::size_t s) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s ? s : 1); }
void  operator delete   (void* p) noexcept { std::free (p); }
void  operator delete[] (void* p) noexcept { std::free (p); }
void  operator delete   (void* p, std::size_t) noexcept { std::free (p); }
void  operator delete[] (void* p, std::size_t) noexcept { std::free (p); }
// Aligned overloads too — SeamAllocator<64> uses ::operator new(size, align_val_t); without these the counter
// goes blind to every SIMD-aligned seam buffer. Portable: _aligned_malloc on MSVC, posix_memalign elsewhere.
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
using NU   = convolution::NonUniformConvolver<>;
using Step = NU::ScheduleStep;

struct Lcg { unsigned long long s; float next() { s = s * 6364136223846793005ULL + 1442695040888963407ULL; return (float) ((s >> 40) & 0xffff) / 32768.0f - 1.0f; } };

// The independent truth: full double-precision causal time-domain convolution. O(N·L) — short IRs only.
static std::vector<float> convDirect (const std::vector<float>& x, const std::vector<float>& h)
{
    std::vector<float> y (x.size(), 0.0f);
    for (int n = 0; n < (int) x.size(); ++n)
    {
        double acc = 0.0;
        const int jmax = std::min ((int) h.size() - 1, n);
        for (int j = 0; j <= jmax; ++j) acc += (double) h[(std::size_t) j] * (double) x[(std::size_t) (n - j)];
        y[(std::size_t) n] = (float) acc;
    }
    return y;
}

// Reference via the proven uniform PartitionedConvolver (fast; used for the full-length IR null).
static std::vector<float> convPartitioned (const std::vector<float>& x, const std::vector<float>& h, int P, int maxIr)
{
    convolution::PartitionedConvolver<> pc; pc.prepare (P, maxIr); pc.setIr (h.data(), (int) h.size());
    std::vector<float> y (x.size(), 0.0f);
    pc.process (x.data(), y.data(), (int) x.size());
    return y;
}

// Run the NUPC over x in host blocks of `block` (block may exceed any stage size or be 1 / non-pow2).
static std::vector<float> runNU (NU& nu, const std::vector<float>& x, int block)
{
    std::vector<float> y (x.size(), 0.0f);
    for (int i = 0; i < (int) x.size(); i += block)
    {
        const int m = std::min (block, (int) x.size() - i);
        nu.process (&x[(std::size_t) i], &y[(std::size_t) i], m);
    }
    return y;
}

static double maxDiff (const std::vector<float>& a, const std::vector<float>& b, int from, int to)
{
    double m = 0.0; for (int i = from; i < to; ++i) m = std::max (m, (double) std::fabs (a[(std::size_t) i] - b[(std::size_t) i])); return m;
}
static double peak (const std::vector<float>& a) { double m = 0.0; for (float v : a) m = std::max (m, (double) std::fabs (v)); return m; }

// A random IR + a random input with unit impulses added on the stage boundaries (0, P0±1, 2P0±1, Bmax±1, L-1).
static void makeSignals (int L, int N, int P0, int Bmax, unsigned seed, std::vector<float>& h, std::vector<float>& x)
{
    Lcg r { seed };
    h.assign ((std::size_t) L, 0.0f); for (auto& v : h) v = 0.05f * r.next();
    x.assign ((std::size_t) N, 0.0f); for (auto& v : x) v = 0.3f  * r.next();
    for (int p : { 0, P0 - 1, P0, 2*P0 - 1, 2*P0, Bmax - 1, Bmax, 2*Bmax, L - 1 })
        if (p >= 0 && p < N) x[(std::size_t) p] += 1.0f;
}

int main()
{
    std::printf ("felitronics::convolution NonUniformConvolver tests\n");

    // ---------------------------------------------------------------------------------------------
    test::group ("NULL vs direct time-domain convolution (short IRs, capped-doubling default)");
    {
        for (int L : { 40, 130, 400, 1500 })
            for (int block : { 1, 64, 200, 512 })
            {
                std::vector<float> h, x; makeSignals (L, L + 3000, 128, 4096, 11u + (unsigned) (L * 31 + block), h, x);
                NU nu; test::ok (nu.prepare (128, 4096, L), "prepare capped");
                nu.setIr (h.data(), L);
                const std::vector<float> y = runNU (nu, x, block);
                const std::vector<float> ref = convDirect (x, h);
                const double rel = maxDiff (y, ref, 0, (int) x.size()) / (peak (ref) + 1e-30);
                if (! (rel < 2e-4)) std::printf ("      L=%d block=%d rel=%.3e\n", L, block, rel);
                test::ok (rel < 2e-4, "NUPC == direct convolution (rel<2e-4) at L=" + std::to_string (L) + " block=" + std::to_string (block));
            }
    }

    // ---------------------------------------------------------------------------------------------
    test::group ("TRUE sample-zero-latency: impulse in → IR out with NO shift");
    {
        const int L = 1700;
        std::vector<float> h, x; makeSignals (L, L + 500, 128, 4096, 4242u, h, x);
        // impulse at t=0: the whole output must equal h (y[0]==h[0]) — any latency k would zero y[0..k).
        std::vector<float> imp (L + 500, 0.0f); imp[0] = 1.0f;
        NU nu; nu.prepare (128, 4096, L); nu.setIr (h.data(), L);
        const std::vector<float> y = runNU (nu, imp, 64);
        double e = 0.0; for (int i = 0; i < L; ++i) e = std::max (e, (double) std::fabs (y[(std::size_t) i] - h[(std::size_t) i]));
        std::printf ("      impulse@0: max|y-h| over the whole IR = %.3e (y[0]=%.6f, h[0]=%.6f)\n", e, y[0], h[0]);
        test::ok (e < 1e-4, "y == h exactly, y[0]==h[0] → ZERO added latency (not shifted by any stage's block)");

        // impulse at t0 (mid-stream): y[t0+i] == h[i], again no shift. Buffer must span t0 + L so the whole
        // impulse response plays out.
        const int t0 = 1000;
        std::vector<float> imp2 (t0 + L + 200, 0.0f); imp2[(std::size_t) t0] = 1.0f;
        NU nu2; nu2.prepare (128, 4096, L); nu2.setIr (h.data(), L);
        const std::vector<float> y2 = runNU (nu2, imp2, 100);
        double e2 = 0.0; for (int i = 0; i < L; ++i) e2 = std::max (e2, (double) std::fabs (y2[(std::size_t) (t0 + i)] - h[(std::size_t) i]));
        test::ok (e2 < 1e-4, "mid-stream impulse response == IR, correctly time-aligned");
    }

    // ---------------------------------------------------------------------------------------------
    test::group ("NULL vs PartitionedConvolver (full 131072-tap IR; capped / tuned / octave schedules)");
    {
        const int L = 131072;
        const int N = L + 8192 + 512;
        std::vector<float> h, x; makeSignals (L, N, 128, 4096, 20260704u, h, x);
        const std::vector<float> ref = convPartitioned (x, h, 128, L);   // fast exact reference (same value for every NU below)
        const double refPk = peak (ref) + 1e-30;

        auto nullAgainstRef = [&] (const char* name, NU& nu, int block)
        {
            nu.setIr (h.data(), L);
            const std::vector<float> y = runNU (nu, x, block);
            const double rel = maxDiff (y, ref, 0, N) / refPk;
            if (! (rel < 5e-4)) std::printf ("      %s block=%d rel=%.3e\n", name, block, rel);
            test::ok (rel < 5e-4, std::string (name) + " == PartitionedConvolver (rel<5e-4) block=" + std::to_string (block));
        };

        { NU nu; test::ok (nu.prepare (128, 4096, L), "prepare capped (Bmax=4096)"); nullAgainstRef ("capped", nu, 512); }
        { NU nu; test::ok (nu.prepare (128, 4096, L), "prepare capped, block=8192 > Bmax"); nullAgainstRef ("capped", nu, 8192); }
        { NU nu; test::ok (nu.prepare (128, 4096, L), "prepare capped, non-pow2 block=257"); nullAgainstRef ("capped", nu, 257); }

        Step tuned[] { {128,3},{512,7},{4096,31} };      // head128 + 128×3 + 512×7 + 4096×31 = 131072
        { NU nu; test::ok (nu.prepareWithSchedule (128, tuned, 3, L), "prepare tuned {128x3,512x7,4096x31}"); nullAgainstRef ("tuned", nu, 2048); }

        std::vector<Step> oct; for (int b = 128; b <= 65536; b <<= 1) oct.push_back ({ b, 1 });   // pure octave to 65536
        { NU nu; test::ok (nu.prepareWithSchedule (128, oct.data(), (int) oct.size(), L), "prepare pure-octave {128..65536}"); nullAgainstRef ("octave", nu, 64); }

        { NU nu; test::ok (nu.prepare (64, 2048, L), "prepare P0=64 Bmax=2048"); nullAgainstRef ("P0=64", nu, 300); }
    }

    // ---------------------------------------------------------------------------------------------
    test::group ("in-place aliasing (in == out) matches out-of-place");
    {
        const int L = 3000, N = L + 2000;
        std::vector<float> h, x; makeSignals (L, N, 128, 4096, 77u, h, x);
        NU a; a.prepare (128, 4096, L); a.setIr (h.data(), L);
        const std::vector<float> outOfPlace = runNU (a, x, 128);
        NU b; b.prepare (128, 4096, L); b.setIr (h.data(), L);
        std::vector<float> io = x;                                     // process in place
        for (int i = 0; i < N; i += 128) { const int m = std::min (128, N - i); b.process (&io[(std::size_t) i], &io[(std::size_t) i], m); }
        test::ok (maxDiff (io, outOfPlace, 0, N) < 1e-6, "in-place output bit-matches out-of-place");
    }

    // ---------------------------------------------------------------------------------------------
    test::group ("short IR in a large schedule (inactive/partial stages) + IR swap");
    {
        const int maxIr = 131072, N = 5000;
        // A schedule sized for 131072, but a SHORT IR → only the first stages have any active partitions.
        std::vector<float> hShort, x; makeSignals (300, N, 128, 4096, 5u, hShort, x);
        NU nu; test::ok (nu.prepare (128, 4096, maxIr), "prepare for maxIr=131072");
        nu.setIr (hShort.data(), 300);
        const std::vector<float> y = runNU (nu, x, 200);
        const std::vector<float> ref = convDirect (x, hShort);
        test::ok (maxDiff (y, ref, 0, N) / (peak (ref) + 1e-30) < 2e-4, "short IR (300 taps) correct in a 131072-sized schedule");

        // swap to a different, longer IR on the same prepared engine (setIr resets state) → still correct.
        std::vector<float> hLong; Lcg r { 909 }; hLong.assign (9000, 0.0f); for (auto& v : hLong) v = 0.04f * r.next();
        nu.setIr (hLong.data(), 9000);
        const std::vector<float> y2 = runNU (nu, x, 200);
        const std::vector<float> ref2 = convDirect (x, hLong);
        test::ok (maxDiff (y2, ref2, 0, N) / (peak (ref2) + 1e-30) < 2e-4, "IR swap (300→9000 taps) correct");
    }

    // ---------------------------------------------------------------------------------------------
    test::group ("no allocation in process()");
    {
        const int L = 20000;
        std::vector<float> h, x; makeSignals (L, 4096, 128, 4096, 33u, h, x);
        NU nu; nu.prepare (128, 4096, L); nu.setIr (h.data(), L);
        std::vector<float> y (4096, 0.0f);
        nu.process (x.data(), y.data(), 2048);                        // warm: cross a few chunk boundaries
        nu.process (&x[2048], &y[2048], 2048);
        const long before = g_allocs.load();
        nu.process (x.data(), y.data(), 2048);                        // steady-state incl. big-stage FFT firings
        nu.process (&x[2048], &y[2048], 2048);
        test::okNoAlloc (g_allocs.load() == before, "process() performed zero heap allocations");
    }

    // ---------------------------------------------------------------------------------------------
    test::group ("prepare validation + unprepared guards");
    {
        NU nu;
        test::ok (! nu.prepare (127, 4096, 1000), "non-pow2 head rejected");
        test::ok (! nu.prepare (128, 3000, 1000), "non-pow2 maxBlock rejected");
        // an unprepared engine is a safe no-op
        std::vector<float> in (64, 0.1f), out (64, -9.0f);
        nu.process (in.data(), out.data(), 64);
        test::ok (out[0] == -9.0f, "process() before prepare() is a no-op (no write, no crash)");
        float dummy = 1.0f; nu.setIr (&dummy, 1);                     // setIr before prepare — no-op, no crash
        test::ok (true, "setIr() before prepare() did not crash");

        // schedule that violates the zero-latency recurrence B_{s+1} = (count_s+1)*B_s is rejected
        Step bad[] { {128,1},{512,1} };                              // 512 != (1+1)*128 = 256
        test::ok (! nu.prepareWithSchedule (128, bad, 2, 2048), "schedule violating offset==B_s rejected");
        Step good[] { {128,1},{256,3},{1024,1} };                   // 256=2*128, 1024=(3+1)*256 — valid; covers 2048 incl head
        test::ok (nu.prepareWithSchedule (128, good, 3, 2048), "valid non-trivial recurrence accepted (covers maxIr exactly)");
        // a valid recurrence that does NOT cover maxIrSamples is rejected (else setIr would silently truncate the IR)
        test::ok (! nu.prepareWithSchedule (128, good, 3, 4096), "valid recurrence but under-covering maxIr rejected");
        // first stage block must equal the head (B_0 == P0)
        Step badFirst[] { {256,1} };
        test::ok (! nu.prepareWithSchedule (128, badFirst, 1, 256), "first stage block != head size rejected");
        // null steps pointer with numSteps>0 is rejected (no deref)
        test::ok (! nu.prepareWithSchedule (128, nullptr, 1, 256), "null schedule with numSteps>0 rejected");
    }

    // ---------------------------------------------------------------------------------------------
    test::group ("C>1 stage with partial activity (1 < activeParts < C), later partitions nulled");
    {
        // head128 + {128×3}=[128,512) + {512×7}=[512,4096). An IR of 1300 taps makes the 512-stage's partitions
        // at offsets 512 and 1024 active (2 of 7) — exercises the multi-partition FDL delay with a partial fill,
        // and N=3200 > 512 + 2·512 so BOTH active tail partitions' outputs are nulled (not just partition 0).
        Step sch[] { {128,3},{512,7} };
        const int L = 1300, N = 3200;
        std::vector<float> h, x; makeSignals (L, N, 128, 512, 1234u, h, x);
        NU nu; test::ok (nu.prepareWithSchedule (128, sch, 2, 4096), "prepare {128x3,512x7}");
        nu.setIr (h.data(), L);
        const std::vector<float> y = runNU (nu, x, 128);
        const std::vector<float> ref = convDirect (x, h);
        test::ok (maxDiff (y, ref, 0, N) / (peak (ref) + 1e-30) < 2e-4, "partial multi-partition (512) stage correct — both active partitions nulled");
    }

    return test::report();
}

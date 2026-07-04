// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// Cross-backend proof for the optional pffft SIMD backend. The SIMD audio path is validated by NULLING
// MatrixConvolver<PffftRealFft> against MatrixConvolver<ScalarRadix2Real> — identical prepare / operator IR /
// input -> identical output to FFT float tolerance, across EVERY routing topology (LRDiag, MSDiag, Full,
// mono) plus a mid-stream topology swap and an odd host-block size. MSDiag matters most: its on-the-fly M/S
// spectral view (viewSpec = ½(X_L ± X_R), elementwise) must stay valid in pffft's opaque z-order layout — it
// does, because the layout is a fixed linear map of the transform. Plus: a direct forward→inverse identity, a
// prepare() admissibility floor, an alignment PROBE (the engine must hand the backend >=16B-aligned spectrum
// rows — pffft's zconvolve precondition, an x86 movaps would fault otherwise), no heap alloc in process(),
// and the C1 rejection (pffft is a RealFftBackend but NOT PackedHermitianSpectrum, so the design-FFT gate
// refuses it at compile time — asserted here at runtime for visibility).

#include <felitronics_test.h>
#include <felitronics/convolution/MatrixConvolver.h>
#include <felitronics/fftpffft/PffftRealFft.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <new>
#if defined(_WIN32)
 #include <malloc.h>   // _aligned_malloc / _aligned_free (MSVC has no posix_memalign)
#endif
#include <string>
#include <vector>

// --- allocation counter (global operator new/delete; aligned overloads too, so SeamAllocator's aligned
//     new is not invisible). Windows-portable per the house pattern. ---
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
using Scalar = core::fft::ScalarRadix2Real;
using Pf     = fftpffft::PffftRealFft;
using McS    = convolution::MatrixConvolver<Scalar>;
using McP    = convolution::MatrixConvolver<Pf>;

struct Lcg { unsigned long long s; float next() { s = s * 6364136223846793005ULL + 1442695040888963407ULL; return (float) ((s >> 40) & 0xffff) / 32768.0f - 1.0f; } };

static const int P = 64, irMax = 400, len = 200, xfade = 128, N = 6000;

template <class MC>
static void runStereo (MC& mc, const std::vector<float>& xL, const std::vector<float>& xR,
                       std::vector<float>& oL, std::vector<float>& oR, int block = 512)
{
    oL = xL; oR = xR;
    for (int o = 0; o < N; o += block) { float* io[2] { oL.data() + o, oR.data() + o }; mc.process (io, io, 2, std::min (block, N - o)); }
}

static double maxDiff (const std::vector<float>& a, const std::vector<float>& b)
{ double m = 0.0; for (std::size_t i = 0; i < a.size(); ++i) m = std::max (m, (double) std::fabs (a[i] - b[i])); return m; }
static double peak (const std::vector<float>& a)
{ double m = 0.0; for (float x : a) m = std::max (m, (double) std::fabs (x)); return m; }

// alignment probe backend: a scalar transform that flags any SPECTRUM-side pointer pffft would need
// >=16B-aligned (forward/inverse `spec`, zconvolve's a/b/acc). The real (time-domain) side is bounced by
// PffftRealFft and so is intentionally not checked — this probes exactly pffft's zero-copy precondition.
static std::atomic<int> g_misaligned { 0 };
static inline bool ali16 (const void* p) { return (reinterpret_cast<std::uintptr_t> (p) & 15u) == 0; }
struct ProbeFft
{
    Scalar inner;
    static constexpr int spectrumFloats (int n) noexcept { return Scalar::spectrumFloats (n); }
    bool prepare (int n) noexcept { return inner.prepare (n); }
    void forward (const float* r, float* spec) noexcept { if (! ali16 (spec)) g_misaligned.fetch_add (1); inner.forward (r, spec); }
    void inverse (const float* spec, float* r) noexcept { if (! ali16 (spec)) g_misaligned.fetch_add (1); inner.inverse (spec, r); }
    void spectralMultiplyAdd (const float* a, const float* b, float* acc) const noexcept
    { if (! ali16 (a) || ! ali16 (b) || ! ali16 (acc)) g_misaligned.fetch_add (1); inner.spectralMultiplyAdd (a, b, acc); }
    int size() const noexcept { return inner.size(); }
};
static_assert (core::fft::RealFftBackend<ProbeFft>, "ProbeFft must satisfy the seam");

int main()
{
    std::printf ("felitronics::fftpffft PffftRealFft cross-backend tests (pffft simd_size=%d)\n", Pf::simdWidth());

    test::group ("pffft SIMD kernel is active (not a scalar fallback)");
    test::ok (Pf::simdWidth() == 4, "pffft_simd_size()==4 — the SSE/NEON kernel is compiled in (the point of the module)");

    Lcg r { 90210 };
    auto mkIr = [&] { std::vector<float> v ((std::size_t) len); for (auto& x : v) x = 0.15f * r.next(); return v; };
    const std::vector<float> irM = mkIr(), irS = mkIr(), irL = mkIr(), irR = mkIr();
    const std::vector<float> irLL = mkIr(), irLR = mkIr(), irRL = mkIr(), irRR = mkIr();

    std::vector<float> xL (N), xR (N);
    for (int i = 0; i < N; ++i) { xL[(std::size_t) i] = (float) (0.5 * std::sin (2.0 * core::kPi * 500.0 * i / 48000.0)); xR[(std::size_t) i] = 0.4f * r.next(); }

    auto report2 = [&] (const std::string& topo, const std::vector<float>& sL, const std::vector<float>& sR,
                        const std::vector<float>& pL, const std::vector<float>& pR)
    {
        const double eL = maxDiff (sL, pL), eR = maxDiff (sR, pR), pk = std::max (peak (sL), peak (sR));
        std::printf ("      %-18s scalar-vs-pffft err L=%.2e R=%.2e (peak %.2e, rel %.2e)\n",
                     topo.c_str(), eL, eR, pk, std::max (eL, eR) / (pk + 1e-30));
        const double tol = 1e-4 * pk + 1e-6;   // float SSE/NEON pffft vs double-twiddle scalar radix-2; observed ~2e-7,
                                               // so ~500x margin over float noise yet far below any structural (layout/
                                               // normalization/MAC) bug, which would null at >=1e-2
        test::ok (eL < tol, topo + ": pffft nulls scalar (L)");
        test::ok (eR < tol, topo + ": pffft nulls scalar (R)");
    };

    // --- cross-backend NULL over every topology ---
    test::group ("cross-backend NULL: MatrixConvolver<Pffft> == MatrixConvolver<Scalar>");
    {
        std::vector<float> sL, sR, pL, pR;
        { McS s; McP p; s.prepare (P, irMax, xfade, 2); p.prepare (P, irMax, xfade, 2);
          const float* b[2] { irL.data(), irR.data() };
          s.setOperator (McS::Topology::LRDiag, b, 2, len); p.setOperator (McP::Topology::LRDiag, b, 2, len);
          runStereo (s, xL, xR, sL, sR); runStereo (p, xL, xR, pL, pR); report2 ("LRDiag", sL, sR, pL, pR); }

        { McS s; McP p; s.prepare (P, irMax, xfade, 2); p.prepare (P, irMax, xfade, 2);
          const float* b[2] { irM.data(), irS.data() };
          s.setOperator (McS::Topology::MSDiag, b, 2, len); p.setOperator (McP::Topology::MSDiag, b, 2, len);
          runStereo (s, xL, xR, sL, sR); runStereo (p, xL, xR, pL, pR); report2 ("MSDiag", sL, sR, pL, pR); }

        { McS s; McP p; s.prepare (P, irMax, xfade, 2); p.prepare (P, irMax, xfade, 2);
          const float* b[4] { irLL.data(), irLR.data(), irRL.data(), irRR.data() };
          s.setOperator (McS::Topology::Full, b, 4, len); p.setOperator (McP::Topology::Full, b, 4, len);
          runStereo (s, xL, xR, sL, sR); runStereo (p, xL, xR, pL, pR); report2 ("Full", sL, sR, pL, pR); }

        // odd host block (non-power-of-two, straddles internal partition boundaries) — parity must hold
        { McS s; McP p; s.prepare (P, irMax, xfade, 2); p.prepare (P, irMax, xfade, 2);
          const float* b[2] { irM.data(), irS.data() };
          s.setOperator (McS::Topology::MSDiag, b, 2, len); p.setOperator (McP::Topology::MSDiag, b, 2, len);
          runStereo (s, xL, xR, sL, sR, 257); runStereo (p, xL, xR, pL, pR, 257); report2 ("MSDiag@block257", sL, sR, pL, pR); }
    }

    // --- mono ---
    test::group ("cross-backend NULL: mono");
    {
        McS s; McP p; s.prepare (P, irMax, xfade, 1); p.prepare (P, irMax, xfade, 1);
        s.setIr (irM.data(), len); p.setIr (irM.data(), len);
        std::vector<float> ys = xL, yp = xL;
        for (int o = 0; o < N; o += 512) { float* a[1] { ys.data() + o }; float* b[1] { yp.data() + o }; s.process (a, a, 1, std::min (512, N - o)); p.process (b, b, 1, std::min (512, N - o)); }
        const double e = maxDiff (ys, yp), pk = peak (ys);
        std::printf ("      mono err=%.2e (peak %.2e)\n", e, pk);
        test::ok (e < 1e-4 * pk + 1e-6, "mono: pffft nulls scalar");   // same gate as the stereo topologies
    }

    // --- mid-stream topology swap: both backends run the same MSDiag->Full swap; parity incl. the fade ---
    test::group ("cross-backend NULL: mid-stream MSDiag->Full swap (fade + warm FDL)");
    {
        McS s; McP p; s.prepare (P, irMax, xfade, 2); p.prepare (P, irMax, xfade, 2);
        { const float* b[2] { irM.data(), irS.data() }; s.setOperator (McS::Topology::MSDiag, b, 2, len); p.setOperator (McP::Topology::MSDiag, b, 2, len); }
        std::vector<float> sL = xL, sR = xR, pL = xL, pR = xR;
        bool sw = false;
        for (int o = 0; o < N; o += 256)
        {
            if (! sw && o >= 2000)
            { const float* b[4] { irLL.data(), irLR.data(), irRL.data(), irRR.data() };
              s.setOperator (McS::Topology::Full, b, 4, len); p.setOperator (McP::Topology::Full, b, 4, len); sw = true; }
            float* a[2] { sL.data() + o, sR.data() + o }; float* b[2] { pL.data() + o, pR.data() + o };
            const int m = std::min (256, N - o); s.process (a, a, 2, m); p.process (b, b, 2, m);
        }
        report2 ("swap MSDiag>Full", sL, sR, pL, pR);
    }

    // --- direct adapter: forward -> inverse == identity (spectrum buffers must be aligned) ---
    test::group ("PffftRealFft forward->inverse == identity");
    {
        Pf f; test::ok (f.prepare (256), "prepare(256)"); test::ok (f.size() == 256, "size()==256");
        core::fft::AlignedVector<float> spec (256, 0.0f);
        std::vector<float> x (256), y (256);
        for (auto& v : x) v = r.next();
        f.forward (x.data(), spec.data());
        f.inverse (spec.data(), y.data());
        double e = 0.0; for (int i = 0; i < 256; ++i) e = std::max (e, (double) std::fabs (y[(std::size_t) i] - x[(std::size_t) i]));
        std::printf ("      round-trip max err=%.2e\n", e);
        test::ok (e < 1e-4, "forward then inverse recovers the input (1/N normalized)");
    }

    // --- prepare() admissibility floor ---
    test::group ("PffftRealFft prepare admissibility");
    {
        Pf f;
        test::ok (! f.prepare (16),  "N=16 rejected (pffft real needs N>=32)");
        test::ok (! f.prepare (100), "N=100 (non-pow2) rejected");
        test::ok (! f.prepare (0),   "N=0 rejected");
        test::ok (! f.prepare (-256), "N<0 rejected");
        test::ok (! f.prepare (1 << 27), "N=2^27 rejected (>2^26 pre-guard; upstream pffft would assert)");
        test::ok (  f.prepare (32),  "N=32 accepted (pffft real minimum)");
        test::ok (  f.prepare (512), "N=512 accepted"); test::ok (f.size() == 512, "size()==512");
        test::ok (! f.prepare (100), "a FAILED re-prepare is rejected");
        test::ok (  f.size() == 0,   "size()==0 after a failed re-prepare (old plan released, left unprepared)");
    }

    // --- alignment probe: the engine hands the backend >=16B-aligned spectrum rows on every path ---
    test::group ("alignment probe: engine supplies SIMD-aligned spectrum rows (zconvolve precondition)");
    {
        g_misaligned.store (0);
        convolution::MatrixConvolver<ProbeFft> a; a.prepare (P, irMax, xfade, 2);
        { const float* b[2] { irM.data(), irS.data() }; a.setOperator (convolution::MatrixConvolver<ProbeFft>::Topology::MSDiag, b, 2, len); }
        std::vector<float> oL, oR; runStereo (a, xL, xR, oL, oR);
        { const float* b[4] { irLL.data(), irLR.data(), irRL.data(), irRR.data() }; a.setOperator (convolution::MatrixConvolver<ProbeFft>::Topology::Full, b, 4, len); }
        runStereo (a, xL, xR, oL, oR);
        convolution::MatrixConvolver<ProbeFft> c; c.prepare (P, irMax, xfade, 2);
        { const float* b[2] { irL.data(), irR.data() }; c.setOperator (convolution::MatrixConvolver<ProbeFft>::Topology::LRDiag, b, 2, len); }
        runStereo (c, xL, xR, oL, oR);
        std::printf ("      spectrum-side misaligned pointers seen: %d\n", g_misaligned.load());
        test::ok (g_misaligned.load() == 0, "every forward/inverse/MAC spectrum pointer is >=16B aligned");
    }

    // --- RT-safety: no heap allocation in process() ---
    test::group ("no heap allocation in pffft process()");
    {
        McP p; p.prepare (P, irMax, xfade, 2);
        { const float* b[2] { irM.data(), irS.data() }; p.setOperator (McP::Topology::MSDiag, b, 2, len); }
        std::vector<float> l (512, 0.2f), rr (512, -0.1f); float* io[2] { l.data(), rr.data() };
        p.process (io, io, 2, 512);
        const long before = g_allocs.load();
        p.process (io, io, 2, 512);
        p.process (io, io, 2, 512);
        test::okNoAlloc (g_allocs.load() == before, "pffft process() performed zero heap allocations");
    }

    // --- C1: the design/audio safety split (compile-enforced; asserted here for visibility) ---
    test::group ("C1: pffft z-order is rejected from design-time FFTs");
    test::ok (  core::fft::RealFftBackend<Pf>,          "PffftRealFft IS a RealFftBackend (valid audio backend)");
    test::ok (! core::fft::PackedHermitianSpectrum<Pf>, "PffftRealFft is NOT PackedHermitianSpectrum (design-FFT gate rejects it)");

    return felitronics::test::report();
}

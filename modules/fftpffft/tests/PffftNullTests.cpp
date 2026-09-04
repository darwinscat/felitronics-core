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
#include <felitronics/convolution/NonUniformConvolver.h>
#include <felitronics/convolution/MatrixConvolverNupc.h>
#include <felitronics/fftpffft/PffftRealFft.h>
#include <felitronics/fftpffft/PffftOrderedRealFft.h>
#include <felitronics/analysis/SpectrumPane.h>
#include <felitronics/analysis/MultiResSpectrumPane.h>
#include <felitronics/analysis/MultiResSpectrumPaneFast.h>
#include <memory>

#include <algorithm>
#include <atomic>
#include <chrono>
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
   #if defined(__cpp_exceptions) || defined(_CPPUNWIND)
    if (p == nullptr) throw std::bad_alloc();
   #else
    if (p == nullptr) std::abort();   // the wasm-audio tier compiles -fno-exceptions, where `throw` is a PARSE error
   #endif
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

    // --- NUPC cross-backend NULL: the SHIPPING pffft path through the NON-UNIFORM stage array. This is the
    //     only null of NonUniformConvolver<Pffft> — its per-stage distinct FFT sizes + memcpy between FDL rows
    //     are exactly where a z-order / row-alignment bug would hide (the scalar tests can't see it). ---
    test::group ("cross-backend NULL: NonUniformConvolver<Pffft> == NonUniformConvolver<Scalar>");
    {
        using NuS = convolution::NonUniformConvolver<Scalar>;
        using NuP = convolution::NonUniformConvolver<Pf>;
        const int nuLen = 5000, nuMax = 8192, nuN = 9000;
        Lcg rn { 24680 };
        std::vector<float> h ((std::size_t) nuLen); for (auto& v : h) v = 0.06f * rn.next();
        std::vector<float> xin ((std::size_t) nuN);  for (auto& v : xin) v = 0.3f * rn.next();
        for (int p : { 0, 127, 128, 4095, 4096, nuLen - 1 }) if (p >= 0 && p < nuN) xin[(std::size_t) p] += 1.0f;  // stage-boundary impulses

        for (int block : { 128, 512, 257 })   // small + typical + non-power-of-two
        {
            NuS s; NuP pf;
            test::ok (s.prepare (128, nuMax, nuLen) && pf.prepare (128, nuMax, nuLen), "NUPC prepare (both backends)");
            s.setIr (h.data(), nuLen); pf.setIr (h.data(), nuLen);
            std::vector<float> ys (nuN, 0.0f), yp (nuN, 0.0f);
            for (int o = 0; o < nuN; o += block)
            { const int m = std::min (block, nuN - o); s.process (&xin[(std::size_t) o], &ys[(std::size_t) o], m); pf.process (&xin[(std::size_t) o], &yp[(std::size_t) o], m); }
            const double e = maxDiff (ys, yp), pk = peak (ys);
            std::printf ("      NUPC capped block=%3d  err=%.2e (peak %.2e, rel %.2e)\n", block, e, pk, e / (pk + 1e-30));
            test::ok (e < 1e-4 * pk + 1e-6, "NUPC pffft nulls scalar (capped, block=" + std::to_string (block) + ")");
        }
    }

    // --- MatrixConvolverNupc (Phase-2 matrix facade) cross-backend NULL: pffft z-order through the per-stage
    //     channel-indexed FDL + LRDiag routing (mono + LRDiag) ---
    test::group ("cross-backend NULL: MatrixConvolverNupc<Pffft> == <Scalar> (mono / LRDiag / MSDiag / Full)");
    {
        using McnS = convolution::MatrixConvolverNupc<Scalar>;
        using McnP = convolution::MatrixConvolverNupc<Pf>;
        const int nuLen = 5000, nuMax = 8192, nuN = 9000;
        Lcg rn { 13579 };
        std::vector<float> hL ((std::size_t) nuLen), hR ((std::size_t) nuLen);
        for (auto& v : hL) v = 0.06f * rn.next(); for (auto& v : hR) v = 0.06f * rn.next();
        std::vector<float> xL ((std::size_t) nuN), xR ((std::size_t) nuN);
        for (auto& v : xL) v = 0.3f * rn.next(); for (auto& v : xR) v = 0.3f * rn.next();

        // mono
        { McnS s; McnP p; s.prepare (128, nuMax, 128, 1); p.prepare (128, nuMax, 128, 1);
          s.setIr (hL.data(), nuLen); p.setIr (hL.data(), nuLen);
          std::vector<float> ys = xL, yp = xL;
          for (int o = 0; o < nuN; o += 512) { const int m = std::min (512, nuN - o); const float* a[1] { ys.data() + o }; const float* b[1] { yp.data() + o }; float* ao[1] { ys.data() + o }; float* bo[1] { yp.data() + o }; s.process (a, ao, 1, m); p.process (b, bo, 1, m); }
          const double e = maxDiff (ys, yp), pk = peak (ys);
          std::printf ("      MatrixConvolverNupc mono err=%.2e (peak %.2e, rel %.2e)\n", e, pk, e / (pk + 1e-30));
          test::ok (e < 1e-4 * pk + 1e-6, "MatrixConvolverNupc mono: pffft nulls scalar"); }

        // LRDiag, non-power-of-two block
        { McnS s; McnP p; s.prepare (128, nuMax, 128, 2); p.prepare (128, nuMax, 128, 2);
          const float* bk[2] { hL.data(), hR.data() };
          s.setOperator (McnS::Topology::LRDiag, bk, 2, nuLen); p.setOperator (McnP::Topology::LRDiag, bk, 2, nuLen);
          std::vector<float> sL = xL, sR = xR, pL = xL, pR = xR;
          for (int o = 0; o < nuN; o += 257) { const int m = std::min (257, nuN - o);
            const float* a[2] { sL.data() + o, sR.data() + o }; float* ao[2] { sL.data() + o, sR.data() + o };
            const float* b[2] { pL.data() + o, pR.data() + o }; float* bo[2] { pL.data() + o, pR.data() + o };
            s.process (a, ao, 2, m); p.process (b, bo, 2, m); }
          const double e = std::max (maxDiff (sL, pL), maxDiff (sR, pR)), pk = std::max (peak (sL), peak (sR));
          std::printf ("      MatrixConvolverNupc LRDiag@257 err=%.2e (peak %.2e, rel %.2e)\n", e, pk, e / (pk + 1e-30));
          test::ok (e < 1e-4 * pk + 1e-6, "MatrixConvolverNupc LRDiag: pffft nulls scalar"); }

        // MSDiag — the per-stage ½(X_L±X_R) view must stay valid in pffft's z-order layout
        { McnS s; McnP p; s.prepare (128, nuMax, 128, 2); p.prepare (128, nuMax, 128, 2);
          std::vector<float> hS ((std::size_t) nuLen); for (auto& v : hS) v = 0.05f * rn.next();
          const float* bk[2] { hL.data(), hS.data() };
          s.setOperator (McnS::Topology::MSDiag, bk, 2, nuLen); p.setOperator (McnP::Topology::MSDiag, bk, 2, nuLen);
          std::vector<float> sL = xL, sR = xR, pL = xL, pR = xR;
          for (int o = 0; o < nuN; o += 512) { const int m = std::min (512, nuN - o);
            const float* a[2] { sL.data() + o, sR.data() + o }; float* ao[2] { sL.data() + o, sR.data() + o };
            const float* b[2] { pL.data() + o, pR.data() + o }; float* bo[2] { pL.data() + o, pR.data() + o };
            s.process (a, ao, 2, m); p.process (b, bo, 2, m); }
          const double e = std::max (maxDiff (sL, pL), maxDiff (sR, pR)), pk = std::max (peak (sL), peak (sR));
          std::printf ("      MatrixConvolverNupc MSDiag err=%.2e (rel %.2e)\n", e, e / (pk + 1e-30));
          test::ok (e < 1e-4 * pk + 1e-6, "MatrixConvolverNupc MSDiag: pffft nulls scalar"); }

        // Full — 4-bank cross routing on the per-stage z-order FDLs
        { McnS s; McnP p; s.prepare (128, nuMax, 128, 2); p.prepare (128, nuMax, 128, 2);
          std::vector<float> h3 ((std::size_t) nuLen), h4 ((std::size_t) nuLen); for (auto& v : h3) v = 0.04f * rn.next(); for (auto& v : h4) v = 0.04f * rn.next();
          const float* bk[4] { hL.data(), hR.data(), h3.data(), h4.data() };
          s.setOperator (McnS::Topology::Full, bk, 4, nuLen); p.setOperator (McnP::Topology::Full, bk, 4, nuLen);
          std::vector<float> sL = xL, sR = xR, pL = xL, pR = xR;
          for (int o = 0; o < nuN; o += 257) { const int m = std::min (257, nuN - o);
            const float* a[2] { sL.data() + o, sR.data() + o }; float* ao[2] { sL.data() + o, sR.data() + o };
            const float* b[2] { pL.data() + o, pR.data() + o }; float* bo[2] { pL.data() + o, pR.data() + o };
            s.process (a, ao, 2, m); p.process (b, bo, 2, m); }
          const double e = std::max (maxDiff (sL, pL), maxDiff (sR, pR)), pk = std::max (peak (sL), peak (sR));
          std::printf ("      MatrixConvolverNupc Full@257 err=%.2e (rel %.2e)\n", e, e / (pk + 1e-30));
          test::ok (e < 1e-4 * pk + 1e-6, "MatrixConvolverNupc Full: pffft nulls scalar"); }
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

    //==========================================================================================
    // The ORDERED sibling — the packed-Hermitian layout, bin for bin, so the spectrum panes can ride SIMD.
    using Po = fftpffft::PffftOrderedRealFft;
    using Sc = core::fft::ScalarRadix2Real;

    test::group ("PffftOrderedRealFft: the ordered transform IS the packed-Hermitian layout (vs the scalar reference, re AND im)");
    {
        std::uint64_t seed = 0x1234567ull;
        auto uni = [&] { seed ^= seed >> 12; seed ^= seed << 25; seed ^= seed >> 27; return (float) ((seed * 0x2545F4914F6CDD1Dull) >> 40) / 8388608.0f - 1.0f; };
        for (int n : { 32, 64, 256, 1024, 16384 })
        {
            Po po; Sc sc;
            test::ok (po.prepare (n) && sc.prepare (n), "prepare " + std::to_string (n));
            std::vector<float> x ((std::size_t) n), a ((std::size_t) n), b ((std::size_t) n);
            for (auto& v : x) v = uni();
            po.forward (x.data(), a.data()); sc.forward (x.data(), b.data());
            float maxAbs = 0.0f, worst = 0.0f;
            for (int i = 0; i < n; ++i) maxAbs = std::max (maxAbs, std::fabs (b[(std::size_t) i]));
            for (int i = 0; i < n; ++i) worst = std::max (worst, std::fabs (a[(std::size_t) i] - b[(std::size_t) i]));
            test::ok (worst <= 2.0e-6f * maxAbs, "N=" + std::to_string (n) + ": every packed float (DC, Nyq, re/im) matches the scalar within 2e-6 of full scale (worst " + std::to_string (worst / maxAbs) + ")");
            // basis vectors pin each slot on its own: a constant → s[0] = N; (−1)^n → s[1] = N and NOT s[N−1]
            // (FFTPACK's own end slot); an impulse at n = 3 → Re/Im k = cos/−sin (6πk/N), the e^{−jωn} sign
            for (int i = 0; i < n; ++i) x[(std::size_t) i] = 1.0f;
            po.forward (x.data(), a.data());
            test::approx (a[0], (double) n, 1e-3 * n, "N=" + std::to_string (n) + ": constant → s[0] = N (DC slot)");
            test::approx (a[1], 0.0, 1e-3 * n, "  … and s[1] = 0");
            for (int i = 0; i < n; ++i) x[(std::size_t) i] = (i & 1) ? -1.0f : 1.0f;
            po.forward (x.data(), a.data());
            test::approx (a[1], (double) n, 1e-3 * n, "  (−1)^n → s[1] = N (the Nyquist slot)");
            test::approx (a[(std::size_t) (n - 1)], 0.0, 1e-3 * n, "  … and s[N−1] = 0 (not FFTPACK's end slot)");
            for (int i = 0; i < n; ++i) x[(std::size_t) i] = (i == 3) ? 1.0f : 0.0f;
            po.forward (x.data(), a.data());
            for (int k : { 1, 2, n / 4 })
            {
                const double th = 2.0 * 3.14159265358979323846 * 3.0 * (double) k / (double) n;
                test::approx (a[(std::size_t) (2 * k)],     std::cos (th), 1e-5, "  impulse at 3 → Re[" + std::to_string (k) + "] = cos");
                test::approx (a[(std::size_t) (2 * k + 1)], -std::sin (th), 1e-5, "  impulse at 3 → Im[" + std::to_string (k) + "] = −sin (e^{−jωn})");
            }
            // a bin-centred cosine: its Re lands in bin 5 with the scalar's sign convention (e^{-jωn}), Im ≈ 0
            for (int i = 0; i < n; ++i) x[(std::size_t) i] = std::cos (2.0f * 3.14159265f * 5.0f * (float) i / (float) n);
            po.forward (x.data(), a.data());
            test::approx (a[10], (double) n / 2.0, 1e-3 * n, "N=" + std::to_string (n) + ": cos at bin 5 → Re[5] = N/2 at s[10]");
            test::approx (a[11], 0.0, 1e-3 * n, "  … Im[5] = 0 at s[11]");
            for (int i = 0; i < n; ++i) x[(std::size_t) i] = std::sin (2.0f * 3.14159265f * 5.0f * (float) i / (float) n);
            po.forward (x.data(), a.data());
            test::approx (a[11], -(double) n / 2.0, 1e-3 * n, "  … sin at bin 5 → Im[5] = −N/2 (the e^{−jωn} convention the designers assume)");
        }
    }

    test::group ("PffftOrderedRealFft forward->inverse == identity; the scalar's spectrum inverts on pffft; misaligned buffers; admissibility; packed MAC");
    {
        Po po; test::ok (po.prepare (4096), "prepare 4096");
        std::vector<float> x (4096), s (4096), y (4096);
        for (std::size_t i = 0; i < x.size(); ++i) x[i] = std::sin (0.013f * (float) i) * 0.5f + (float) (i % 7) * 0.01f;
        po.forward (x.data(), s.data()); po.inverse (s.data(), y.data());
        float worst = 0.0f; for (std::size_t i = 0; i < x.size(); ++i) worst = std::max (worst, std::fabs (x[i] - y[i]));
        test::ok (worst < 1.0e-5f, "round trip within 1e-5 (worst " + std::to_string (worst) + ")");
        // the inverse must read the layout the SCALAR wrote — not merely undo its own forward
        Sc sc; sc.prepare (4096); sc.forward (x.data(), s.data()); po.inverse (s.data(), y.data());
        worst = 0.0f; for (std::size_t i = 0; i < x.size(); ++i) worst = std::max (worst, std::fabs (x[i] - y[i]));
        test::ok (worst < 1.0e-5f, "scalar forward → pffft inverse == identity (the packed layout is shared both ways)");
        // unaligned caller buffers on every side: pffft's alignment need is met by the bounces, not by the caller
        std::vector<float> xu (4097 + 1), su (4096 + 1), yu (4096 + 1);
        for (std::size_t i = 0; i < 4096; ++i) xu[i + 1] = x[i];
        po.forward (xu.data() + 1, su.data() + 1); po.inverse (su.data() + 1, yu.data() + 1);
        worst = 0.0f; for (std::size_t i = 0; i < 4096; ++i) worst = std::max ({ worst, std::fabs (su[i + 1] - s[i]) * 1e-3f, std::fabs (yu[i + 1] - x[i]) });
        test::ok (worst < 1.0e-5f, "real in / spectrum out / spectrum in / real out all at data()+1: same spectrum, same round trip");
        Po bad;
        test::ok (! bad.prepare (16), "N=16 refused (pffft real needs ≥ 32)");
        test::ok (! bad.prepare (48), "N=48 refused (not a power of two)");
        test::ok (  bad.prepare (32), "N=32 admitted");
        Po m; Sc ms; m.prepare (64); ms.prepare (64);
        std::vector<float> a (64), b (64), acc1 (64, 0.25f), acc2 (64, 0.25f);
        for (int i = 0; i < 64; ++i) { a[(std::size_t) i] = 0.01f * (float) i; b[(std::size_t) i] = 1.0f - 0.02f * (float) i; }
        m.spectralMultiplyAdd (a.data(), b.data(), acc1.data()); ms.spectralMultiplyAdd (a.data(), b.data(), acc2.data());
        worst = 0.0f; for (int i = 0; i < 64; ++i) worst = std::max (worst, std::fabs (acc1[(std::size_t) i] - acc2[(std::size_t) i]));
        test::ok (worst < 1.0e-6f, "spectralMultiplyAdd == the scalar packed MAC");
        test::ok (core::fft::PackedHermitianSpectrum<Po>, "PffftOrderedRealFft IS PackedHermitianSpectrum (admitted where bins are read)");
    }

    test::group ("cross-backend NULL: SpectrumPaneT<PffftOrdered> == SpectrumPaneT<Scalar> (columns, both orders)");
    {
        auto pp = std::make_unique<analysis::SpectrumPaneT<Po>>();
        auto ps = std::make_unique<analysis::SpectrumPaneT<Sc>>();
        std::uint64_t seed = 0xBEEFull;
        auto uni = [&] { seed ^= seed >> 12; seed ^= seed << 25; seed ^= seed >> 27; return (float) ((seed * 0x2545F4914F6CDD1Dull) >> 40) / 8388608.0f - 1.0f; };
        analysis::PlotMap pm; pm.width = 900.0f; pm.height = 96.0f; pm.plotBottom = 96.0f; pm.specTop = 6.0; pm.specBottom = -90.0;   // 1 px per dB
        for (int order : { 11, 14, 10 })
        {
            const int n = 1 << order;
            for (int t = 0; t < 3; ++t)
            {
                std::vector<float> fr ((std::size_t) n); for (auto& v : fr) v = uni() * (1.0f + 0.5f * std::sin (0.01f * (float) t));
                std::copy (fr.begin(), fr.end(), pp->frameInput()); pp->ingest (order);
                std::copy (fr.begin(), fr.end(), ps->frameInput()); ps->ingest (order);
            }
            std::vector<float> yp, ys, pkp, pks;
            pp->buildColumns (pm, 48000.0, 4.5, 1000.0, [&] (int, float, float y, float pk) { yp.push_back (y); pkp.push_back (pk); });
            ps->buildColumns (pm, 48000.0, 4.5, 1000.0, [&] (int, float, float y, float pk) { ys.push_back (y); pks.push_back (pk); });
            float worst = 0.0f;
            for (std::size_t i = 0; i < ys.size(); ++i) worst = std::max ({ worst, std::fabs (yp[i] - ys[i]), std::fabs (pkp[i] - pks[i]) });
            test::ok (ys.size() == 901 && worst <= 0.02f, "order " + std::to_string (order) + ": 901 columns, fill + peak within 0.02 dB (worst " + std::to_string (worst) + ")");
        }
        pp->starve(); ps->starve();
        std::vector<float> fr (16384); for (auto& v : fr) v = uni();   // the test's own buffer, before the count starts
        const long before = g_allocs.load();
        std::copy (fr.begin(), fr.end(), pp->frameInput()); pp->ingest (14);
        pp->buildColumns (pm, 48000.0, 4.5, 1000.0, [] (int, float, float, float) {});
        test::okNoAlloc (g_allocs.load() == before, "the pffft pane's ingest + buildColumns allocate nothing");
    }

    test::group ("cross-backend NULL: MultiResSpectrumPaneT<…, PffftOrdered> == scalar (bins, stitched reads)");
    {
        using Mp = analysis::MultiResSpectrumPaneT<14, 4, Po>;
        using Ms = analysis::MultiResSpectrumPaneT<14, 4, Sc>;
        auto mp = std::make_unique<Mp>(); auto ms = std::make_unique<Ms>();
        mp->coverSamples = 1600; ms->coverSamples = 1600;
        std::uint64_t seed = 0xC0FFEEull;
        auto uni = [&] { seed ^= seed >> 12; seed ^= seed << 25; seed ^= seed >> 27; return (float) ((seed * 0x2545F4914F6CDD1Dull) >> 40) / 8388608.0f - 1.0f; };
        for (int t = 0; t < 4; ++t)
        {
            std::vector<float> fr (16384);
            for (std::size_t i = 0; i < fr.size(); ++i) fr[i] = 0.3f * uni() + 0.5f * std::sin (2.0f * 3.14159265f * 1234.5f * (float) i / 48000.0f);
            std::copy (fr.begin(), fr.end(), mp->frameInput()); mp->ingest (14);
            std::copy (fr.begin(), fr.end(), ms->frameInput()); ms->ingest (14);
        }
        float worstBin = 0.0f; int compared = 0;
        for (int k = 0; k < 3; ++k)
            for (int i = 0; i < ms->tierBins (k); ++i)
                if (ms->tierBinDb (k, i) > -90.0f) { worstBin = std::max (worstBin, std::fabs (mp->tierBinDb (k, i) - ms->tierBinDb (k, i))); ++compared; }
        test::ok (compared > 5000 && worstBin <= 0.02f, "every bin above −90 dB in every tier within 0.02 dB (" + std::to_string (compared) + " bins, worst " + std::to_string (worstBin) + ")");
        double worstRead = 0.0;
        for (double f = 30.0; f < 20000.0; f *= 1.02)
            worstRead = std::max ({ worstRead, std::fabs (mp->readDb (f, 48000.0) - ms->readDb (f, 48000.0)), std::fabs (mp->readPeakDb (f, 48000.0) - ms->readPeakDb (f, 48000.0)) });
        test::ok (worstRead <= 0.02, "the stitched fill + peak reads 30 Hz–20 kHz within 0.02 dB (worst " + std::to_string (worstRead) + ")");
        test::ok (mp->subWindows (2, 14) == 3, "the pffft pane runs the same three sub-windows");
    }

    //==========================================================================================
    // The fast sibling has its own paired NULL against MultiResSpectrumPane on the scalar reference
    // (felitronics_multires_fast_pane_tests). What only this suite can check is that it behaves the same
    // on the SIMD backend: that it still NULLs against the pane it copies WITH pffft underneath, and that
    // it is itself backend-independent to the same 0.02 dB the sibling is held to.
    test::group ("NULL: MultiResSpectrumPaneFast on pffft — vs the pane it copies, and vs itself on scalar");
    {
        auto fp = std::make_unique<analysis::MultiResSpectrumPaneFastT<14, 4, Po>>();
        auto fsc = std::make_unique<analysis::MultiResSpectrumPaneFastT<14, 4, Sc>>();
        auto cp = std::make_unique<analysis::MultiResSpectrumPaneT<14, 4, Po>>();
        fp->coverSamples = 1600; fsc->coverSamples = 1600; cp->coverSamples = 1600;
        std::uint64_t seed = 0xBEEF1234ull;
        auto uni = [&] { seed ^= seed >> 12; seed ^= seed << 25; seed ^= seed >> 27; return (float) ((seed * 0x2545F4914F6CDD1Dull) >> 40) / 8388608.0f - 1.0f; };
        for (int t = 0; t < 6; ++t)
        {
            std::vector<float> fr (16384);
            for (std::size_t i = 0; i < fr.size(); ++i) fr[i] = 0.3f * uni() + 0.5f * std::sin (2.0f * 3.14159265f * 1234.5f * (float) i / 48000.0f);
            std::copy (fr.begin(), fr.end(), fp ->frameInput()); fp ->ingest (14);
            std::copy (fr.begin(), fr.end(), fsc->frameInput()); fsc->ingest (14);
            std::copy (fr.begin(), fr.end(), cp ->frameInput()); cp ->ingest (14);
        }
        double sib = 0.0, backend = 0.0;
        for (double f = 30.0; f < 20000.0; f *= 1.01)
            for (double tilt : { 0.0, 4.5 })
            {
                sib = std::max ({ sib, std::fabs (fp->readDb (f, 48000.0, tilt, 1000.0) - cp->readDb (f, 48000.0, tilt, 1000.0)),
                                       std::fabs (fp->readPeakDb (f, 48000.0, tilt, 1000.0) - cp->readPeakDb (f, 48000.0, tilt, 1000.0)) });
                backend = std::max ({ backend, std::fabs (fp->readDb (f, 48000.0, tilt, 1000.0) - fsc->readDb (f, 48000.0, tilt, 1000.0)),
                                              std::fabs (fp->readPeakDb (f, 48000.0, tilt, 1000.0) - fsc->readPeakDb (f, 48000.0, tilt, 1000.0)) });
            }
        test::ok (sib <= 0.02, "on pffft the fast pane still matches the pane it copies within 0.02 dB (worst " + std::to_string (sib) + ")");
        test::ok (backend <= 0.02, "and the fast pane is backend-independent to the same 0.02 dB (worst " + std::to_string (backend) + ")");
        // The columns too: the plan must survive a change of FFT backend unchanged.
        analysis::PlotMap pm; pm.width = 900.0f; pm.height = 300.0f; pm.plotBottom = 300.0f;
        std::vector<float> a, b; double worstPx = 0.0;
        fp ->buildColumns (pm, 48000.0, 4.5, 1000.0, [&] (int, float, float y, float yp) { a.push_back (y); a.push_back (yp); });
        cp ->buildColumns (pm, 48000.0, 4.5, 1000.0, [&] (int, float, float y, float yp) { b.push_back (y); b.push_back (yp); });
        for (std::size_t i = 0; i < a.size() && i < b.size(); ++i) worstPx = std::max (worstPx, (double) std::fabs (a[i] - b[i]));
        test::ok (a.size() == b.size() && worstPx < 0.05, "and buildColumns agrees to well under a pixel on pffft (worst " + std::to_string (worstPx) + " px)");
    }

    //==========================================================================================
    // PERFORMANCE — the reason the ordered backend exists. One UI tick (ingest + 900 columns) per pane on
    // both backends; printed for the record, and the check is the one that matters: with the SIMD kernel
    // compiled in, the pffft pane must be cheaper than the scalar one. (The scalar table on its own is
    // felitronics_spectrum_pane_perf_tests.)
    test::group ("performance: the panes on pffft vs the scalar reference (one tick = ingest + 900 columns)");
    {
        std::uint64_t seed = 0x51EDull;
        auto uni = [&] { seed ^= seed >> 12; seed ^= seed << 25; seed ^= seed >> 27; return (float) ((seed * 0x2545F4914F6CDD1Dull) >> 40) / 8388608.0f - 1.0f; };
        // MINIMUM OF 300 SINGLE TICKS, with the harness's own PRNG fill OUTSIDE the clock. Both halves of
        // that sentence were bought with a CI failure, so they are worth the paragraph.
        //
        // It used to be the best of three windows of a HUNDRED ticks. Same 300 ticks of work, but the wrong
        // shape: a min over windows needs ONE clean window, and a 100-tick window here is ~15 ms, which
        // necessarily spans several scheduler quanta. On a co-tenanted 3-vCPU CI VM there is no clean 15 ms,
        // so best-of-3 — or best-of-30 — buys nothing, and the two panes are timed in SEPARATE windows, so
        // whatever steal each one caught lands in the ratio. That is exactly what happened: the runner read
        // 308 us for the fast pane against ~95 here, but only 323 against ~156 for the pane it copies —
        // 3.24x versus 2.08x. A slower machine scales both alike; only additive contamination of separately
        // timed windows skews them, and the fast pane, timed last, wore it.
        //
        // A single tick is ~60-150 us and fits inside one quantum, so each of the 300 draws is its own
        // chance at a clean sample and the minimum finds one with overwhelming probability. Measured on an
        // 18-core Mac, this pane pair on pffft: the old estimator reads 1.73/1.90/1.81 quiet and
        // 1.83/1.81/1.84/1.87 under twenty busy loops — a 10 % spread; this one reads 2.361/2.373/2.367
        // quiet and 2.389/2.369/2.394/2.375 loaded — 1.4 %, and the loaded and quiet ranges OVERLAP.
        // Contention finer than a tick (memory bandwidth, DVFS) is multiplicative and cancels in a ratio.
        //
        // The fill moves out because it is the harness, not the pane: 16384 PRNG samples cost ~25 us per
        // tick and were being charged to both panes equally, which drags a ratio toward 1. That single
        // change is most of the gap between the old numbers and these (and reconciles the suite with
        // docs/PERF-ANALYZER-MULTIRES.md, whose harness never timed it). NB a min-of-N ratio sits ABOVE a
        // median-of-N ratio — per-tick jitter is proportionally larger on the cheaper pane — so 2.37 here
        // and 1.82 in the doc's median-of-nine table are the same panes under two estimators, not a
        // disagreement.
        auto tickMicros = [&] (auto& pane, int frameSamples, int order)
        {
            analysis::PlotMap pm; pm.width = 900.0f; pm.height = 300.0f; pm.plotBottom = 300.0f;
            volatile float sink = 0.0f; double best = 1e30;
            for (int t = 0; t < 300; ++t)
            {
                float* f = pane.frameInput(); for (int i = 0; i < frameSamples; ++i) f[i] = uni();   // NOT timed
                const auto t0 = std::chrono::steady_clock::now();
                pane.ingest (order);
                pane.buildColumns (pm, 48000.0, 1.5, 1000.0, [&] (int, float, float y, float yp) { sink = sink + y + yp; });
                best = std::min (best, std::chrono::duration<double, std::micro> (std::chrono::steady_clock::now() - t0).count());
            }
            return best;
        };
        for (int o : { 11, 14 })
        {
            auto sp = std::make_unique<analysis::SpectrumPaneT<Sc>>(); auto pp = std::make_unique<analysis::SpectrumPaneT<Po>>();
            const double us = tickMicros (*sp, 1 << o, o), up = tickMicros (*pp, 1 << o, o);
            std::printf ("      classic %-6d scalar %7.0f us   pffft %7.0f us   (%.1fx)\n", 1 << o, us, up, us / up);
            if (Pf::simdWidth() == 4) test::ok (up < us, "classic " + std::to_string (1 << o) + " on pffft is cheaper than on the scalar reference");
        }
        {
            auto ms = std::make_unique<analysis::MultiResSpectrumPaneT<14, 4, Sc>>(); auto mp = std::make_unique<analysis::MultiResSpectrumPaneT<14, 4, Po>>();
            ms->coverSamples = 1600; mp->coverSamples = 1600;
            const double us = tickMicros (*ms, 16384, 14), up = tickMicros (*mp, 16384, 14);
            std::printf ("      multi, hop 1600  scalar %7.0f us   pffft %7.0f us   (%.1fx)\n", us, up, us / up);
            if (Pf::simdWidth() == 4) test::ok (up < us, "the multi-res pane on pffft is cheaper than on the scalar reference");
            test::ok (up < 33333.0 * 0.2, "a multi-res tick on pffft stays under 20 % of a 30 fps frame");

            // The fast sibling on the same backend. With the transform this cheap, what it removes —
            // a log10 and an exp per bin, and the column geometry the fill and the peak each derived
            // for themselves — is most of what is left, so the gap is at its widest here.
            //
            // The bar stayed at 0.85 through the CI failure that prompted the estimator above, and that is
            // the point: the claim was never what broke. With the intrinsic tick measured rather than a
            // window of them plus the harness's PRNG, this pair reads ~2.37x on pffft — twice the margin
            // the bar asks for, and it did not move under twenty busy loops. It cannot fall below 1 by
            // construction either: the fast pane's work is a strict SUBSET of the pane it copies — same
            // ladder (pinned by #107's transform tally), same magnitude loop, minus two libm calls per bin
            // and minus the column geometry the fill and the peak each used to derive for themselves.
            //
            // What could still make it red on a machine nobody here has: a part whose libm is far cheaper
            // relative to its FFT than Apple's. The printed number is the diagnosis — a red at 2.3 means
            // the code changed, a red near 1.0 means look at the machine before the code.
            auto mf = std::make_unique<analysis::MultiResSpectrumPaneFastT<14, 4, Po>>();
            mf->coverSamples = 1600;
            const double uf = tickMicros (*mf, 16384, 14);
            std::printf ("      multi FAST, hop 1600           pffft %7.0f us   (%.2fx the current pane)\n", uf, up / uf);
            if (Pf::simdWidth() == 4) test::ok (uf < 0.85 * up, "the fast pane is materially cheaper than the pane it copies on pffft");
        }
    }

    return felitronics::test::report();
}

// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa. Part of felitronics-core — see LICENSE.

// JUCE-free self-tests for per-band stereo width (MultibandProcessor<StereoWidth>). The two contracts that
// matter (DSP council): (1) all bands neutral → the splitter's allpass reconstruction (sample-wise null);
// (2) the GLOBAL mono-fold invariant survives ANY per-band widths — ½(L+R) == allpass(½(Lin+Rin)) because
// each band's Mid is untouched. Plus: a per-band width is frequency-selective (mono-ing one band removes
// only that band's side); mono is a true passthrough; no alloc.

#include <felitronics_test.h>
#include <felitronics/core/Math.h>
#include <felitronics/multiband/MultibandWidth.h>
#include <felitronics/eq/MultibandSplitter.h>

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

static double rmsTail (const std::vector<float>& v, int from)
{
    double e = 0.0; int n = 0;
    for (int i = from; i < (int) v.size(); ++i) { e += (double) v[i] * v[i]; ++n; }
    return n ? std::sqrt (e / (double) n) : 0.0;
}

int main()
{
    std::printf ("felitronics::multiband MultibandWidth tests\n");
    const double sr = 48000.0, pi = core::kPi;
    const float xover[2] { 200.0f, 2000.0f };
    auto rng = [s = 1ull] () mutable { s = s * 6364136223846793005ull + 1442695040888963407ull; return (float) ((s >> 40) & 0xffff) / 32768.0f - 1.0f; };

    // --- (1) all bands width=1 → sample-wise null against eq::MultibandSplitter allpass reconstruction ---
    test::group ("MultibandWidth unity bands == allpass reconstruction");
    {
        const int N = 8000; std::vector<float> l (N), r (N), l0 (N), r0 (N);
        for (int i = 0; i < N; ++i) { l[i] = l0[i] = 0.4f * rng(); r[i] = r0[i] = 0.4f * rng(); }
        multiband::MultibandWidth<4> mw; mw.prepare (sr, 1024, 2); mw.setNumBands (3); mw.setCrossovers (xover, 2);
        for (int o = 0; o < N; o += 1024) { float* io[2] { l.data() + o, r.data() + o }; mw.process (io, 2, std::min (1024, N - o)); }
        eq::MultibandSplitter<4> refL, refR; refL.prepare (sr, 1); refR.prepare (sr, 1);
        refL.setNumBands (3); refL.setCrossovers (xover, 2); refR.setNumBands (3); refR.setCrossovers (xover, 2);
        double md = 0; float band[4];
        for (int i = 0; i < N; ++i)
        {
            refL.splitSample (0, l0[i], band); const float rl = refL.sumSample (band);
            refR.splitSample (0, r0[i], band); const float rr = refR.sumSample (band);
            if (i >= 1000) { md = std::max (md, (double) std::fabs (l[i] - rl)); md = std::max (md, (double) std::fabs (r[i] - rr)); }
        }
        test::ok (md < 1e-3, "all widths=1 → L/R == the splitter's allpass reconstruction (sample-wise)");
    }

    // --- (2) THE mono-fold invariant: ½(L+R) == allpass(½(Lin+Rin)) for ARBITRARY per-band widths ---
    test::group ("MultibandWidth global mono-fold invariant");
    {
        const int N = 8000; std::vector<float> l (N), r (N), m0 (N);
        for (int i = 0; i < N; ++i) { l[i] = 0.4f * rng(); r[i] = 0.4f * rng(); m0[i] = 0.5f * (l[i] + r[i]); }
        multiband::MultibandWidth<4> mw; mw.prepare (sr, 1024, 2); mw.setNumBands (3); mw.setCrossovers (xover, 2);
        mw.setBandWidth (0, 0.5f); mw.setBandWidth (1, 1.3f); mw.setBandWidth (2, 1.8f); mw.reset();   // arbitrary, settled
        for (int o = 0; o < N; o += 1024) { float* io[2] { l.data() + o, r.data() + o }; mw.process (io, 2, std::min (1024, N - o)); }
        eq::MultibandSplitter<4> refM; refM.prepare (sr, 1); refM.setNumBands (3); refM.setCrossovers (xover, 2);
        double md = 0; float band[4];
        for (int i = 0; i < N; ++i) { refM.splitSample (0, m0[i], band); const float rm = refM.sumSample (band); if (i >= 1000) md = std::max (md, (double) std::fabs (0.5f * (l[i] + r[i]) - rm)); }
        test::ok (md < 1e-3, "½(L+R) folds to the allpass-reconstructed input Mid (Mid untouched at any per-band width)");
    }

    // --- (3) per-band width is FREQUENCY-SELECTIVE: mono-ing the high band removes only the high-freq side ---
    test::group ("MultibandWidth per-band width is frequency-selective");
    {
        const int N = 12000; std::vector<float> l (N), r (N);
        for (int i = 0; i < N; ++i)                                                       // M = 100 Hz, S = 6 kHz
        {
            const float mid = 0.3f * (float) std::sin (2.0 * pi * 100.0  * i / sr);
            const float sde = 0.3f * (float) std::sin (2.0 * pi * 6000.0 * i / sr);
            l[i] = mid + sde; r[i] = mid - sde;                                            // side lives entirely at 6 kHz
        }
        std::vector<float> sideIn (N); for (int i = 0; i < N; ++i) sideIn[i] = 0.5f * (l[i] - r[i]);
        multiband::MultibandWidth<4> mw; mw.prepare (sr, 1024, 2); mw.setNumBands (3); mw.setCrossovers (xover, 2);
        mw.setBandWidth (2, 0.0f); mw.reset();                                             // mono the > 2 kHz band
        for (int o = 0; o < N; o += 1024) { float* io[2] { l.data() + o, r.data() + o }; mw.process (io, 2, std::min (1024, N - o)); }
        std::vector<float> sideOut (N); for (int i = 0; i < N; ++i) sideOut[i] = 0.5f * (l[i] - r[i]);
        test::ok (rmsTail (sideOut, N / 2) < 0.1 * rmsTail (sideIn, N / 2), "width=0 on the high band kills the 6 kHz side (low-band side would survive)");
    }

    // --- (4) mono (< 2 channels) is a TRUE passthrough (no allpass phase) ---
    test::group ("MultibandWidth mono passthrough");
    {
        const int N = 1024; std::vector<float> x (N), x0 (N);
        for (int i = 0; i < N; ++i) x[i] = x0[i] = 0.5f * rng();
        multiband::MultibandWidth<4> mw; mw.prepare (sr, 1024, 2); mw.setNumBands (3); mw.setCrossovers (xover, 2);
        mw.setBandWidth (0, 0.0f); mw.setBandWidth (2, 2.0f);
        float* io[1] { x.data() }; mw.process (io, 1, N);
        double md = 0; for (int i = 0; i < N; ++i) md = std::max (md, (double) std::fabs (x[i] - x0[i]));
        test::ok (md == 0.0, "numChannels < 2 → bit-exact passthrough");
    }

    // --- (5) no allocation in process() ---
    test::group ("MultibandWidth no-alloc");
    {
        multiband::MultibandWidth<4> mw; mw.prepare (sr, 512, 2); mw.setNumBands (3); mw.setCrossovers (xover, 2);
        mw.setBandWidth (0, 0.7f); mw.setBandWidth (2, 1.6f);
        std::vector<float> l (512, 0.3f), r (512, -0.2f); float* io[2] { l.data(), r.data() };
        mw.process (io, 2, 512);
        const long before = g_allocs.load();
        mw.process (io, 2, 512); mw.process (io, 2, 512);
        test::ok (g_allocs.load() == before, "process() did not allocate");
    }

    // --- (6) width=2 on an isolated band doubles THAT band's side (+6 dB), per the M/S side law ---
    test::group ("MultibandWidth width=2 doubles a band's side");
    {
        const int N = 10000; std::vector<float> l (N), r (N), sideIn (N);
        for (int i = 0; i < N; ++i) { const float s = 0.3f * (float) std::sin (2.0 * pi * 600.0 * i / sr); l[i] = s; r[i] = -s; sideIn[i] = 0.5f * (l[i] - r[i]); }   // M=0, S at 600 Hz (mid-band 1)
        multiband::MultibandWidth<4> mw; mw.prepare (sr, 1024, 2); mw.setNumBands (3); mw.setCrossovers (xover, 2);
        mw.setBandWidth (1, 2.0f); mw.reset();
        for (int o = 0; o < N; o += 1024) { float* io[2] { l.data() + o, r.data() + o }; mw.process (io, 2, std::min (1024, N - o)); }
        std::vector<float> sideOut (N); for (int i = 0; i < N; ++i) sideOut[i] = 0.5f * (l[i] - r[i]);
        const double ratio = rmsTail (sideOut, N / 2) / rmsTail (sideIn, N / 2);
        test::ok (ratio > 1.8 && ratio < 2.15, "width=2 on the 600 Hz band → its side ≈ doubled");
    }

    // --- (7) an out-of-range band index is a safe no-op (no UB on the per-band array) ---
    test::group ("MultibandWidth out-of-range band index is safe");
    {
        multiband::MultibandWidth<4> mw; mw.prepare (sr, 1024, 2); mw.setNumBands (3); mw.setCrossovers (xover, 2);
        mw.setBandWidth (1, 1.5f);
        mw.setBandWidth (-1, 0.0f); mw.setBandWidth (4, 0.0f); mw.setBandWidth (99, 0.0f);   // all ignored
        mw.setBandBypass (-1, true); mw.setBandSolo (7, true);
        const int N = 1024; std::vector<float> l (N), r (N); for (int i = 0; i < N; ++i) { l[i] = 0.3f * rng(); r[i] = 0.3f * rng(); }
        float* io[2] { l.data(), r.data() }; mw.process (io, 2, N);
        bool fin = true; for (int i = 0; i < N; ++i) if (! std::isfinite (l[i]) || ! std::isfinite (r[i])) fin = false;
        test::ok (fin && mw.bandWidth (1) == 1.5f && mw.bandWidth (99) == 1.0f, "OOB indices ignored (no crash; the valid band untouched)");
    }

    // --- (8) neutral holds even with a parallel dry blend: widths=1, mix=0.5 → still allpass recon ---
    test::group ("MultibandWidth neutral holds at mix < 1");
    {
        const int N = 8000; std::vector<float> l (N), r (N), l0 (N), r0 (N);
        for (int i = 0; i < N; ++i) { l[i] = l0[i] = 0.4f * rng(); r[i] = r0[i] = 0.4f * rng(); }
        multiband::MultibandWidth<4> mw; mw.prepare (sr, 1024, 2); mw.setNumBands (3); mw.setCrossovers (xover, 2);
        mw.setMix (0.5f);                                                                  // widths all 1 → wet == dry
        for (int o = 0; o < N; o += 1024) { float* io[2] { l.data() + o, r.data() + o }; mw.process (io, 2, std::min (1024, N - o)); }
        eq::MultibandSplitter<4> refL, refR; refL.prepare (sr, 1); refR.prepare (sr, 1);
        refL.setNumBands (3); refL.setCrossovers (xover, 2); refR.setNumBands (3); refR.setCrossovers (xover, 2);
        double md = 0; float band[4];
        for (int i = 0; i < N; ++i) { refL.splitSample (0, l0[i], band); const float rl = refL.sumSample (band); refR.splitSample (0, r0[i], band); const float rr = refR.sumSample (band); if (i >= 1000) { md = std::max (md, (double) std::fabs (l[i] - rl)); md = std::max (md, (double) std::fabs (r[i] - rr)); } }
        test::ok (md < 1e-3, "all widths=1, mix=0.5 → still the allpass reconstruction (neutral wet == dry)");
    }

    // --- (9) the invariants hold at 4 bands too (not a 3-band accident) ---
    test::group ("MultibandWidth 4-band mono-fold invariant");
    {
        const int N = 8000; std::vector<float> l (N), r (N), m0 (N);
        for (int i = 0; i < N; ++i) { l[i] = 0.4f * rng(); r[i] = 0.4f * rng(); m0[i] = 0.5f * (l[i] + r[i]); }
        const float xf4[3] { 120.0f, 800.0f, 5000.0f };
        multiband::MultibandWidth<4> mw; mw.prepare (sr, 1024, 2); mw.setNumBands (4); mw.setCrossovers (xf4, 3);
        mw.setBandWidth (0, 0.5f); mw.setBandWidth (1, 1.2f); mw.setBandWidth (2, 0.8f); mw.setBandWidth (3, 1.7f); mw.reset();
        for (int o = 0; o < N; o += 1024) { float* io[2] { l.data() + o, r.data() + o }; mw.process (io, 2, std::min (1024, N - o)); }
        eq::MultibandSplitter<4> refM; refM.prepare (sr, 1); refM.setNumBands (4); refM.setCrossovers (xf4, 3);
        double md = 0; float band[4];
        for (int i = 0; i < N; ++i) { refM.splitSample (0, m0[i], band); const float rm = refM.sumSample (band); if (i >= 1000) md = std::max (md, (double) std::fabs (0.5f * (l[i] + r[i]) - rm)); }
        test::ok (md < 1e-3, "4 bands, arbitrary widths → ½(L+R) still folds to allpass(M)");
    }

    return test::report();
}

// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// JUCE-free self-tests for the reusable crossover primitives: the Svf AllPass mode, the LR4 `Crossover2`,
// and the `MultibandSplitter` with allpass-compensated FLAT reconstruction. The headline check proves the
// claim the band-compensation rests on — that the LR4 allpass (LP4+HP4) collapses to a SINGLE 2nd-order
// allpass — by nulling one `Svf` AllPass against a full `Crossover2` low+high.

#include <felitronics_test.h>
#include <felitronics/core/Math.h>
#include <felitronics/eq/Svf.h>
#include <felitronics/eq/Crossover2.h>
#include <felitronics/eq/MultibandSplitter.h>

#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdlib>

static std::atomic<long> g_allocs { 0 };
void* operator new      (std::size_t s) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s ? s : 1); }
void* operator new[]    (std::size_t s) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s ? s : 1); }
void  operator delete   (void* p) noexcept { std::free (p); }
void  operator delete[] (void* p) noexcept { std::free (p); }
void  operator delete   (void* p, std::size_t) noexcept { std::free (p); }
void  operator delete[] (void* p, std::size_t) noexcept { std::free (p); }

using namespace felitronics;
static constexpr double kQ = 0.7071067811865476;

// RMS(sum-of-bands)/RMS(input) at frequency f (≈1 ⇒ flat reconstruction), and per-band RMS out.
static double analyze (int numBands, const float* freqs, double f, double sr, double* bandRmsOut)
{
    eq::MultibandSplitter<4> mb; mb.prepare (sr, 1); mb.setNumBands (numBands); mb.setCrossovers (freqs, numBands - 1);
    const int N = 8000; double inSq = 0, outSq = 0, bSq[4] = { 0, 0, 0, 0 }; float band[4];
    for (int i = 0; i < N; ++i)
    {
        const float x = (float) std::sin (2.0 * core::kPi * f * i / sr);
        mb.splitSample (0, x, band);
        const float y = mb.sumSample (band);
        if (i >= N / 2) { inSq += (double) x * x; outSq += (double) y * y; for (int b = 0; b < numBands; ++b) bSq[b] += (double) band[b] * band[b]; }
    }
    for (int b = 0; b < numBands; ++b) bandRmsOut[b] = std::sqrt (bSq[b] / (0.5 * (double) N));
    return std::sqrt (outSq) / std::sqrt (inSq);
}

int main()
{
    std::printf ("felitronics::eq crossover tests\n");
    const double sr = 48000.0, pi = core::kPi;

    // --- Svf AllPass is genuinely flat (was a passthrough stub before) ---
    test::group ("Svf AllPass flat magnitude");
    {
        eq::Svf ap; ap.prepare (sr, 1); ap.setParams (eq::FilterType::AllPass, 1000.0, kQ, 0.0);
        int N = 8000; double inSq = 0, outSq = 0;
        for (int i = 0; i < N; ++i) { const float x = (float) std::sin (2.0 * pi * 1000.0 * i / sr); const float y = ap.processSample (0, x); if (i >= N / 2) { inSq += (double) x * x; outSq += (double) y * y; } }
        test::approx (std::sqrt (outSq / inSq), 1.0, 0.02, "|H|≈1 at fc");
    }

    // --- THE collapse: one 2nd-order AllPass == the 4th-order LR4 allpass (LP4+HP4) ---
    test::group ("LR4 allpass collapses to 2nd-order (Svf AllPass == Crossover2 low+high)");
    {
        eq::Svf ap; ap.prepare (sr, 1); ap.setParams (eq::FilterType::AllPass, 1234.0, kQ, 0.0);
        eq::Crossover2 xo; xo.prepare (sr, 1); xo.setFrequency (1234.0f);
        unsigned long s = 22699; auto rng = [&]() { s = s * 6364136223846793005ULL + 1442695040888963407ULL; return (float) ((s >> 40) & 0xffff) / 32768.0f - 1.0f; };
        int N = 6000; double md = 0;
        for (int i = 0; i < N; ++i)
        {
            const float x = rng();
            const float a = ap.processSample (0, x);
            float lo, hi; xo.processSample (0, x, lo, hi);
            if (i >= 1000) md = std::max (md, (double) std::fabs (a - (lo + hi)));
        }
        test::ok (md < 1e-3, "one AllPass SVF nulls against LP4+HP4 → the LR4 allpass is 2nd-order");
    }

    // --- Crossover2 low+high is flat (allpass) across the spectrum ---
    test::group ("Crossover2 low+high flat (allpass)");
    {
        double worst = 0.0;
        for (double f : { 50.0, 250.0, 1000.0, 4000.0, 12000.0 })
        {
            eq::Crossover2 xo; xo.prepare (sr, 1); xo.setFrequency (1000.0f);
            int N = 8000; double inSq = 0, outSq = 0;
            for (int i = 0; i < N; ++i) { const float x = (float) std::sin (2.0 * pi * f * i / sr); float lo, hi; xo.processSample (0, x, lo, hi); const float y = lo + hi; if (i >= N / 2) { inSq += (double) x * x; outSq += (double) y * y; } }
            worst = std::max (worst, std::fabs (std::sqrt (outSq / inSq) - 1.0));
        }
        test::ok (worst < 0.03, "|low+high|≈|in| at every test freq");
    }

    // --- MultibandSplitter: FLAT reconstruction (the allpass band-compensation) ---
    test::group ("MultibandSplitter flat reconstruction");
    {
        const double testF[] = { 40.0, 150.0, 500.0, 2000.0, 8000.0, 16000.0 };
        double bands[4];

        { const float xf[1] = { 1000.0f }; double w = 0; for (double f : testF) w = std::max (w, std::fabs (analyze (2, xf, f, sr, bands) - 1.0)); test::ok (w < 0.04, "2-band sum flat across spectrum"); }
        { const float xf[2] = { 250.0f, 2500.0f }; double w = 0; for (double f : testF) w = std::max (w, std::fabs (analyze (3, xf, f, sr, bands) - 1.0)); test::ok (w < 0.04, "3-band sum flat (allpass compensation works at f1 and f2)"); }
        { const float xf[3] = { 150.0f, 800.0f, 5000.0f }; double w = 0; for (double f : testF) w = std::max (w, std::fabs (analyze (4, xf, f, sr, bands) - 1.0)); test::ok (w < 0.04, "4-band sum flat"); }
    }

    // --- bands are actually band-limited ---
    test::group ("MultibandSplitter band-limiting");
    {
        const float xf[2] = { 250.0f, 2500.0f }; double bands[4];
        analyze (3, xf, 40.0, sr, bands);                          // 40 Hz: low band dominates
        test::ok (bands[0] > 5.0 * bands[2], "40 Hz → low band >> high band");
        analyze (3, xf, 10000.0, sr, bands);                       // 10 kHz: high band dominates
        test::ok (bands[2] > 5.0 * bands[0], "10 kHz → high band >> low band");
    }

    // --- no allocation in the per-sample split ---
    test::group ("MultibandSplitter no-alloc");
    {
        eq::MultibandSplitter<4> mb; mb.prepare (sr, 1); mb.setNumBands (4);
        float band[4]; mb.splitSample (0, 0.5f, band);
        const long before = g_allocs.load();
        for (int i = 0; i < 512; ++i) mb.splitSample (0, 0.1f, band);
        test::okNoAlloc (g_allocs.load() == before, "splitSample() did not allocate");
    }

    // --- Svf AllPass |H|≈1 across Q / freq / rate (confirms m1=-2k is correct off-Butterworth) ---
    test::group ("Svf AllPass |H|=1 across Q / freq / rate");
    {
        double worst = 0.0;
        for (double rate : { 44100.0, 48000.0, 96000.0 })
            for (double Q : { 0.5, 0.7071067811865476, 2.0, 5.0 })
                for (double fc : { 60.0, 600.0, 6000.0 })
                {
                    eq::Svf ap; ap.prepare (rate, 1); ap.setParams (eq::FilterType::AllPass, fc, Q, 0.0);
                    const double pf = fc * 1.7;
                    int N = 8000; double inSq = 0, outSq = 0;
                    for (int i = 0; i < N; ++i) { const float x = (float) std::sin (2.0 * pi * pf * i / rate); const float y = ap.processSample (0, x); if (i >= N / 2) { inSq += (double) x * x; outSq += (double) y * y; } }
                    worst = std::max (worst, std::fabs (std::sqrt (outSq / inSq) - 1.0));
                }
        test::ok (worst < 0.03, "allpass magnitude unity for ALL Q (m1=-2k correct beyond Butterworth)");
    }

    // --- STRONGEST recon check: splitter sum == an allpass cascade of the crossovers, SAMPLE-WISE ---
    test::group ("MultibandSplitter sum == allpass cascade (sample-wise null)");
    {
        const float xf[3] = { 200.0f, 1100.0f, 6000.0f };
        for (int bands = 3; bands <= 4; ++bands)
        {
            eq::MultibandSplitter<4> mb; mb.prepare (sr, 1); mb.setNumBands (bands); mb.setCrossovers (xf, bands - 1);
            eq::Svf ref[3]; for (int j = 0; j < bands - 1; ++j) { ref[j].prepare (sr, 1); ref[j].setParams (eq::FilterType::AllPass, xf[j], kQ, 0.0); }
            unsigned long s = 999; auto rng = [&]() { s = s * 6364136223846793005ULL + 1442695040888963407ULL; return (float) ((s >> 40) & 0xffff) / 32768.0f - 1.0f; };
            int N = 6000; double md = 0; float band[4];
            for (int i = 0; i < N; ++i)
            {
                const float x = rng();
                mb.splitSample (0, x, band); const float y = mb.sumSample (band);
                float r = x; for (int j = 0; j < bands - 1; ++j) r = ref[j].processSample (0, r);
                if (i >= 1500) md = std::max (md, (double) std::fabs (y - r));
            }
            test::ok (md < 2e-3, bands == 3 ? "3-band sum nulls against AP(f1)·AP(f2)" : "4-band sum nulls against AP(f1)·AP(f2)·AP(f3)");
        }
    }

    // --- Crossover2 is −6 dB at the cutoff (LR4 signature) ---
    test::group ("Crossover2 −6 dB at cutoff");
    {
        eq::Crossover2 xo; xo.prepare (sr, 1); xo.setFrequency (1000.0f);
        int N = 8000; double inSq = 0, loSq = 0, hiSq = 0;
        for (int i = 0; i < N; ++i) { const float x = (float) std::sin (2.0 * pi * 1000.0 * i / sr); float lo, hi; xo.processSample (0, x, lo, hi); if (i >= N / 2) { inSq += (double) x * x; loSq += (double) lo * lo; hiSq += (double) hi * hi; } }
        test::approx (20.0 * std::log10 (std::sqrt (loSq / inSq)), -6.02, 0.5, "low band −6 dB at fc");
        test::approx (20.0 * std::log10 (std::sqrt (hiSq / inSq)), -6.02, 0.5, "high band −6 dB at fc");
    }

    // --- cross-rate flat reconstruction (44.1k + 96k, not just 48k) ---
    test::group ("MultibandSplitter flat at 44.1k and 96k");
    {
        const double testF[] = { 60.0, 400.0, 3000.0, 12000.0 };
        const float xf[2] = { 250.0f, 2500.0f }; double bands[4];
        double w441 = 0, w96 = 0;
        for (double f : testF) { w441 = std::max (w441, std::fabs (analyze (3, xf, f, 44100.0, bands) - 1.0)); w96 = std::max (w96, std::fabs (analyze (3, xf, f, 96000.0, bands) - 1.0)); }
        test::ok (w441 < 0.04, "3-band flat @ 44.1 kHz");
        test::ok (w96  < 0.04, "3-band flat @ 96 kHz");
    }

    // --- per-channel state independence (stereo: no cross-talk) ---
    test::group ("MultibandSplitter channel independence");
    {
        eq::MultibandSplitter<4> mb; mb.prepare (sr, 2); mb.setNumBands (3); const float xf[2] = { 250.0f, 2500.0f }; mb.setCrossovers (xf, 2);
        float b0[4], b1[4]; double ch1max = 0;
        for (int i = 0; i < 2000; ++i) { mb.splitSample (0, (i % 64 == 0) ? 1.0f : 0.0f, b0); mb.splitSample (1, 0.0f, b1); ch1max = std::max (ch1max, (double) std::fabs (mb.sumSample (b1))); }
        test::ok (ch1max == 0.0, "channel 1 stays silent while channel 0 is driven (independent state)");
    }

    // --- setCrossovers robustness: hostile inputs clamp/repair, no UB, still flat ---
    test::group ("MultibandSplitter setCrossovers robustness");
    {
        eq::MultibandSplitter<4> mb; mb.prepare (sr, 1); mb.setNumBands (3);
        const float hostile[2] = { 30000.0f, 100.0f };   // descending + above Nyquist (would be lo>hi clamp UB)
        mb.setCrossovers (hostile, 2);
        float band[4]; double inSq = 0, outSq = 0;
        for (int i = 0; i < 8000; ++i) { const float x = (float) std::sin (2.0 * pi * 1000.0 * i / sr); mb.splitSample (0, x, band); const float y = mb.sumSample (band); if (i >= 4000) { inSq += (double) x * x; outSq += (double) y * y; } }
        test::ok (std::isfinite (outSq) && std::fabs (std::sqrt (outSq / inSq) - 1.0) < 0.05, "hostile crossovers clamped → still flat, no NaN/UB");
    }

    return test::report();
}

// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// JUCE-free self-tests for the broadband Compressor: static gain reduction == the analytic compressor
// law, channel-LINK (identical gain on all channels → no image shift), lookahead latency/alignment,
// and no-allocation-in-process().

#include <felitronics_test.h>
#include <felitronics/dynamics/Compressor.h>
#include <felitronics/core/Math.h>

#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

// global allocation counter (no-alloc-in-process proof)
static std::atomic<long> g_allocs { 0 };
void* operator new      (std::size_t s) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s ? s : 1); }
void* operator new[]    (std::size_t s) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s ? s : 1); }
void  operator delete   (void* p) noexcept { std::free (p); }
void  operator delete[] (void* p) noexcept { std::free (p); }
void  operator delete   (void* p, std::size_t) noexcept { std::free (p); }
void  operator delete[] (void* p, std::size_t) noexcept { std::free (p); }

using namespace felitronics;

static double rmsTail (const std::vector<float>& v, double frac)
{
    const int from = (int) (v.size() * (1.0 - frac));
    double s = 0.0; int c = 0;
    for (int i = from; i < (int) v.size(); ++i) { s += (double) v[i] * v[i]; ++c; }
    return c ? std::sqrt (s / c) : 0.0;
}

int main()
{
    std::printf ("felitronics::dynamics Compressor tests\n");
    const double fs = 48000.0;

    // --- static GR == analytic law: RMS detector, thr -20, ratio 4, input RMS -8 dBFS → GR -9 dB ---
    test::group ("Compressor static GR == analytic (4:1)");
    {
        const int n = (int) (1.0 * fs);
        const double f = 300.0;
        const double A = std::sqrt (2.0) * core::dbToGain (-8.0);   // sine whose RMS is -8 dBFS
        std::vector<float> x (n);
        for (int i = 0; i < n; ++i) x[i] = (float) (A * std::sin (2.0 * core::kPi * f * i / fs));
        std::vector<float> y = x;
        float* ch[1] { y.data() };

        dynamics::Compressor comp;
        comp.prepare (fs, n, 1);
        dynamics::CompressorParams p;
        p.detector = dynamics::Detector::Rms; p.rmsWindowMs = 30.0;
        p.mode = dynamics::Mode::DownCompress; p.thresholdDb = -20.0; p.ratio = 4.0; p.kneeDb = 0.0;
        p.attackMs = 5.0; p.releaseMs = 80.0; p.makeupDb = 0.0;
        comp.setParams (p);
        comp.process (ch, 1, n);

        const double appliedDb = 20.0 * std::log10 (rmsTail (y, 0.2) / rmsTail (x, 0.2));
        test::approx (appliedDb, -9.0, 0.6, "12 dB over @4:1 → ~-9 dB applied");
        test::ok (comp.gainReductionDb() < -7.5 && comp.gainReductionDb() > -10.5, "GR meter ~ -9 dB");
    }

    // --- channel link: L loud, R quiet → BOTH get the SAME gain (stereo image preserved) ---
    test::group ("Compressor stereo link (identical gain both channels)");
    {
        const int n = (int) (1.0 * fs);
        const double f = 200.0;
        std::vector<float> xl (n), xr (n);
        for (int i = 0; i < n; ++i)
        {
            xl[i] = (float) (0.5 * std::sin (2.0 * core::kPi * f * i / fs));
            xr[i] = (float) (0.1 * std::sin (2.0 * core::kPi * f * i / fs));
        }
        std::vector<float> yl = xl, yr = xr;
        float* ch[2] { yl.data(), yr.data() };

        dynamics::Compressor comp;
        comp.prepare (fs, n, 2);
        dynamics::CompressorParams p;
        p.detector = dynamics::Detector::Rms; p.rmsWindowMs = 30.0; p.link = dynamics::LinkMode::Max;
        p.thresholdDb = -24.0; p.ratio = 4.0; p.kneeDb = 0.0; p.attackMs = 5.0; p.releaseMs = 80.0;
        comp.setParams (p);
        comp.process (ch, 2, n);

        const double gL = 20.0 * std::log10 (rmsTail (yl, 0.2) / rmsTail (xl, 0.2));
        const double gR = 20.0 * std::log10 (rmsTail (yr, 0.2) / rmsTail (xr, 0.2));
        test::approx (gL, gR, 0.1, "L and R receive the same gain (no image shift)");
        test::ok (gL < -1.0, "compression actually happened");
    }

    // --- lookahead: with no compression (high threshold) the compressor is a pure delay of `lookahead` ---
    test::group ("Compressor lookahead latency + alignment");
    {
        const int look = 64, n = 400;
        std::vector<float> x (n, 0.0f); x[0] = 1.0f;
        std::vector<float> y = x;
        float* ch[1] { y.data() };

        dynamics::Compressor comp;
        comp.prepare (fs, n, 1, 5.0);            // maxLookahead 5 ms ⊇ 64 samples
        dynamics::CompressorParams p;
        p.thresholdDb = 24.0;                    // above any signal → gain == 1 → pure delay
        p.lookaheadMs = (double) look / fs * 1000.0;
        comp.setParams (p);
        test::ok (comp.latencySamples() == look, "reports lookahead as latency");

        comp.process (ch, 1, n);
        test::approx (y[(std::size_t) look], 1.0, 1e-4, "impulse appears at exactly `lookahead`");
        double leak = 0.0; for (int i = 0; i < look; ++i) leak = std::max (leak, (double) std::fabs (y[(std::size_t) i]));
        test::ok (leak < 1e-6, "nothing before the lookahead delay (no off-by-one)");
    }

    // --- no allocation during process() ---
    test::group ("Compressor no-alloc in process()");
    {
        const int n = 512;
        std::vector<float> a (n, 0.2f), b (n, 0.2f);
        float* ch[2] { a.data(), b.data() };
        dynamics::Compressor comp;
        comp.prepare (fs, n, 2, 10.0);
        dynamics::CompressorParams p; p.thresholdDb = -30.0; p.ratio = 4.0; p.lookaheadMs = 2.0;
        comp.setParams (p);                      // allocations allowed up to here
        const long before = g_allocs.load();
        comp.process (ch, 2, n);
        comp.process (ch, 2, n);
        const long after = g_allocs.load();
        test::ok (after == before, "process() performed zero heap allocations");
    }

    return test::report();
}

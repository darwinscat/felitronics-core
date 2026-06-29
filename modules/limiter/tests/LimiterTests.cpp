// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa. Part of felitronics-core — see LICENSE.

// JUCE-free self-tests for the true-peak limiter: the GUARANTEE (output true-peak, measured by an
// INDEPENDENT oversampler, stays at/below the ceiling even when the input's inter-sample peaks exceed
// it), transparency below the ceiling, latency, and no-allocation-in-process().

#include <felitronics_test.h>
#include <felitronics/limiter/TruePeakLimiter.h>
#include <felitronics/oversampling/PolyphaseOversampler.h>
#include <felitronics/core/Math.h>

#include <atomic>
#include <algorithm>
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

// Independent true-peak (dBTP) of a baseband buffer, via a FRESH 8x oversampler (not the limiter's).
static double measureTruePeakDb (const std::vector<float>& x)
{
    oversampling::PolyphaseOversampler m; m.prepare (8, 1, 32);
    std::vector<float> osb (x.size() * 8);
    const float* xi[1] { x.data() }; float* oo[1] { osb.data() };
    m.upsample (xi, 1, (int) x.size(), oo);
    double mx = 0.0;
    for (std::size_t i = x.size(); i + x.size() < osb.size(); ++i) mx = std::max (mx, (double) std::fabs (osb[i]));
    return 20.0 * std::log10 (mx > 1e-9 ? mx : 1e-9);
}

static double rmsTail (const std::vector<float>& v, double frac)
{
    const int from = (int) (v.size() * (1.0 - frac));
    double s = 0.0; int c = 0;
    for (int i = from; i < (int) v.size(); ++i) { s += (double) v[i] * v[i]; ++c; }
    return c ? std::sqrt (s / c) : 0.0;
}

int main()
{
    std::printf ("felitronics::limiter tests\n");
    const double sr = 48000.0;

    // --- THE GUARANTEE: a signal whose inter-sample peaks exceed the ceiling → output true-peak ≤ ceiling ---
    test::group ("Limiter true-peak guarantee (ISP)");
    {
        const int n = 8192; const double f = sr * 0.25, A = 1.0;     // sample peaks ~0.765, TRUE peak ~1.0 (0 dBTP)
        std::vector<float> x (n);
        for (int i = 0; i < n; ++i) x[i] = (float) (A * std::sin (2.0 * core::kPi * f * i / sr + 0.7));
        const double inTp = measureTruePeakDb (x);
        test::ok (inTp > -1.0 + 0.3, "input true-peak is genuinely above the ceiling (real ISP case)");

        std::vector<float> y = x; float* ch[1] { y.data() };
        limiter::TruePeakLimiter lim; lim.prepare (sr, n, 1);
        limiter::TruePeakLimiterParams p; p.ceilingDbTp = -1.0; p.releaseMs = 50.0; p.lookaheadMs = 1.0; p.oversampleFactor = 4;
        lim.setParams (p);
        lim.process (ch, 1, n);

        const double outTp = measureTruePeakDb (y);
        test::ok (outTp <= -1.0 + 0.5, "output true-peak ≤ ceiling (+0.5 dB downsample-ripple margin)");
        test::ok (lim.gainReductionDb() < -0.2, "limiter actually engaged");
    }

    // --- transparency: a signal below the ceiling passes ~unchanged (just the latency) ---
    test::group ("Limiter transparent below ceiling");
    {
        const int n = 4096; const double f = 1000.0, A = 0.3;        // ~-10 dBFS, well below -1 dBTP
        std::vector<float> x (n);
        for (int i = 0; i < n; ++i) x[i] = (float) (A * std::sin (2.0 * core::kPi * f * i / sr));
        std::vector<float> y = x; float* ch[1] { y.data() };
        limiter::TruePeakLimiter lim; lim.prepare (sr, n, 1);
        limiter::TruePeakLimiterParams p; p.ceilingDbTp = -1.0; p.lookaheadMs = 1.0;
        lim.setParams (p);
        lim.process (ch, 1, n);
        test::approx (rmsTail (y, 0.3) / rmsTail (x, 0.3), 1.0, 0.05, "amplitude preserved (transparent)");
        test::ok (lim.gainReductionDb() > -0.3, "no meaningful gain reduction below ceiling");
    }

    // --- latency = oversampler round-trip + lookahead (baseband) ---
    test::group ("Limiter latency");
    {
        limiter::TruePeakLimiter lim; lim.prepare (sr, 512, 2);
        limiter::TruePeakLimiterParams p; p.lookaheadMs = 1.0; p.oversampleFactor = 4;
        lim.setParams (p);
        const int look = (int) std::lround (1.0 * 0.001 * sr);       // 48 baseband
        test::ok (lim.latencySamples() == (32 - 1) + look, "latency = (tpp-1) + lookahead");
    }

    // --- no allocation during process() ---
    test::group ("Limiter no-alloc in process()");
    {
        const int n = 512;
        std::vector<float> a (n, 0.6f), b (n, 0.6f);
        float* ch[2] { a.data(), b.data() };
        limiter::TruePeakLimiter lim; lim.prepare (sr, n, 2);
        limiter::TruePeakLimiterParams p; p.ceilingDbTp = -1.0; p.lookaheadMs = 1.0;
        lim.setParams (p);
        const long before = g_allocs.load();
        lim.process (ch, 2, n);
        lim.process (ch, 2, n);
        const long after = g_allocs.load();
        test::ok (after == before, "process() performed zero heap allocations");
    }

    return test::report();
}

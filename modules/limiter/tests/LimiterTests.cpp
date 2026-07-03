// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

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
#include <limits>
#include <string>
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
        test::okNoAlloc (after == before, "process() performed zero heap allocations");
    }

    // --- lifecycle/misuse: process before/after a FAILED prepare + oversized blocks must not OOB osBuf ---
    test::group ("TruePeakLimiter: reject process before / after failed prepare, and oversized blocks");
    {
        limiter::TruePeakLimiter lim;                                // NOT prepared (maxCh == 0)
        float a[64] {}, b[64] {}; float* io[2] { a, b };
        lim.process (io, 2, 16);                                     // maxCh==0 → no-op
        lim.prepare (48000.0, 16, 2, 2);                            // tapsPerPhase=2 < 4 → oversampler prepare fails → stays unprepared
        lim.process (io, 2, 16);                                     // must no-op, not run an unprepared oversampler
        lim.prepare (48000.0, 16, 2, 32);                           // valid
        lim.process (io, 2, 16);                                     // works
        lim.process (io, 2, 64);                                     // numSamples=64 > maxBlock=16 → no-op, must not overrun osBuf
        test::ok (true, "no OOB across failed-prepare / oversized-block process (ASan/UBSan is the check)");
    }

    // --- FALSIFICATION: the SlidingMax deque must survive a strictly-decreasing run that fills it ---
    // (insert-before-expire with capacity == W overwrites the head — the current max — on the W+1-th push)
    test::group ("SlidingMax survives a full deque (4,3,2,1,0 @ W=4)");
    {
        limiter::detail::SlidingMax sm; sm.prepare (4);
        const float in[5]   = { 4.0f, 3.0f, 2.0f, 1.0f, 0.0f };
        const float want[5] = { 4.0f, 4.0f, 4.0f, 4.0f, 3.0f };   // window of 4 → last covers {3,2,1,0}
        for (int i = 0; i < 5; ++i)
            test::approx ((double) sm.push (in[i]), (double) want[i], 0.0, "sliding max after push #" + std::to_string (i));
    }

    // --- SlidingMax == brute-force window max over a hostile sequence (long decreasing ramps + noise) ---
    test::group ("SlidingMax == brute-force window max (hostile sequence)");
    {
        const int W = 5;
        limiter::detail::SlidingMax sm; sm.prepare (W);
        unsigned long long s = 42;
        auto rnd = [&]() { s = s * 6364136223846793005ULL + 1442695040888963407ULL; return (float) ((s >> 40) & 0xffff) / 65536.0f; };
        std::vector<float> xs;
        for (int i = 0; i < 40; ++i)
        {
            if (i % 4 == 0) { float v = 1.0f + rnd(); for (int k = 0; k < 50; ++k) { xs.push_back (v); v *= 0.98f; } }  // strictly-decreasing run ≫ W
            else            { for (int k = 0; k < 30; ++k) xs.push_back (rnd()); }
        }
        bool allOk = true;
        for (int i = 0; i < (int) xs.size() && allOk; ++i)
        {
            const float got = sm.push (xs[i]);
            float want = 0.0f;
            for (int j = std::max (0, i - (W - 1)); j <= i; ++j) want = std::max (want, xs[j]);
            allOk = (got == want);
        }
        test::ok (allOk, "sliding max exact for every push (W=5, decreasing runs fill the deque)");
    }

    // --- FALSIFICATION: ceiling guarantee on a >20 ms decaying ramp (deque-overflow end-to-end) ---
    test::group ("Limiter ceiling guarantee on a decaying ramp");
    {
        const int n = 8192;
        std::vector<float> x (n, 0.0f);
        const int rampLen = (int) (0.025 * sr);                      // 25 ms strictly-decreasing > the 20 ms window
        for (int i = 0; i < rampLen; ++i) x[i] = 1.0f - 0.7f * (float) i / (float) rampLen;   // 1.0 → 0.3
        std::vector<float> y = x; float* ch[1] { y.data() };
        limiter::TruePeakLimiter lim; lim.prepare (sr, n, 1);
        limiter::TruePeakLimiterParams p; p.ceilingDbTp = 20.0 * std::log10 (0.25); p.releaseMs = 0.5; p.lookaheadMs = 1.0;
        lim.setParams (p);
        lim.process (ch, 1, n);
        double mx = 0.0; for (float v : y) mx = std::max (mx, (double) std::fabs (v));
        test::ok (mx <= 0.25 * 1.03, "every output sample ≤ ceiling on the smooth ramp (got max " + std::to_string (mx) + ")");
    }

    // --- FALSIFICATION: release must start after the LOOKAHEAD window, not the fixed 20 ms max window ---
    test::group ("Limiter release follows the lookahead window (no 20 ms hold)");
    {
        const int n = 4096;
        std::vector<float> x (n, 0.0f);
        for (int i = 0; i < 96; ++i)  x[i] = (float) std::sin (2.0 * core::kPi * 0.25 * i + 0.7);            // 2 ms 0 dBFS fs/4 burst
        for (int i = 96; i < n; ++i)  x[i] = (float) (0.25 * std::sin (2.0 * core::kPi * 1000.0 * i / sr));  // then a −12 dB tone
        std::vector<float> y = x; float* ch[1] { y.data() };
        limiter::TruePeakLimiter lim; lim.prepare (sr, n, 1);
        limiter::TruePeakLimiterParams p; p.ceilingDbTp = -6.0; p.releaseMs = 1.0; p.lookaheadMs = 1.0;
        lim.setParams (p);
        lim.process (ch, 1, n);
        auto rmsWin = [] (const std::vector<float>& v, int a, int b) {
            double s2 = 0.0; for (int i = a; i < b; ++i) s2 += (double) v[i] * v[i]; return std::sqrt (s2 / std::max (1, b - a)); };
        const double got  = rmsWin (y, (int) (0.010 * sr), (int) (0.020 * sr));   // 10–20 ms: burst long gone at 1 ms lookahead
        const double want = rmsWin (x, (int) (0.010 * sr), (int) (0.020 * sr));
        test::approx (got / want, 1.0, 0.06, "tone recovered a few ms after the burst (gain not held for 20 ms)");
    }

    // --- non-finite params must not poison the stream (house rule: clamp + reject non-finite) ---
    test::group ("Limiter non-finite params rejected");
    {
        const int n = 1024;
        std::vector<float> y (n);
        for (int i = 0; i < n; ++i) y[i] = (float) (0.8 * std::sin (2.0 * core::kPi * 997.0 * i / sr));
        float* ch[1] { y.data() };
        limiter::TruePeakLimiter lim; lim.prepare (sr, n, 1);
        limiter::TruePeakLimiterParams p;
        p.ceilingDbTp = std::numeric_limits<double>::quiet_NaN();
        p.releaseMs   = std::numeric_limits<double>::quiet_NaN();
        p.lookaheadMs = std::numeric_limits<double>::quiet_NaN();
        lim.setParams (p);
        lim.process (ch, 1, n);
        bool finite = true; for (float v : y) finite &= (bool) std::isfinite (v);
        test::ok (finite, "NaN ceiling/release/lookahead → finite output");
        test::ok (lim.latencySamples() >= 0, "latency stays non-negative under NaN lookahead");
    }

    // --- latency query before prepare() must not report the unprepared oversampler's tpp-1 == -1 ---
    test::group ("Limiter latencySamples before prepare == 0");
    {
        limiter::TruePeakLimiter lim;
        test::ok (lim.latencySamples() == 0, "unprepared limiter reports 0 latency, not -1");
    }

    return test::report();
}

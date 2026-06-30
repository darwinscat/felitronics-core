// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// JUCE-free self-tests for the Dynamic EQ band: a 0 dB band below threshold is a passthrough; cut-when-loud
// (de-ess) attenuates a loud in-band tone; boost-when-quiet lifts a quiet in-band tone; an off-band tone is
// untouched (frequency-selective detection); ratio=1 is a plain static bell; stereo is linked; no alloc.

#include <felitronics_test.h>
#include <felitronics/core/Math.h>
#include <felitronics/dynamiceq/DynamicEqBand.h>

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
using DB   = dynamiceq::DynamicEqBand;
using DP   = dynamiceq::DynamicEqBandParams;
using Mode = dynamiceq::DynamicEqMode;

static double rmsTail (const std::vector<float>& v, int from)
{
    double e = 0.0; int n = 0;
    for (int i = from; i < (int) v.size(); ++i) { e += (double) v[i] * v[i]; ++n; }
    return n ? std::sqrt (e / (double) n) : 0.0;
}

static std::vector<float> runTone (const DP& p, double f, double amp, int N, double sr)
{
    DB d; d.prepare (sr, 1); d.setParams (p);
    std::vector<float> y (N);
    for (int i = 0; i < N; ++i) y[i] = (float) (amp * std::sin (2.0 * core::kPi * f * i / sr));
    for (int o = 0; o < N; o += 512) { float* io[1] { y.data() + o }; d.process (io, 1, std::min (512, N - o)); }
    return y;
}

int main()
{
    std::printf ("felitronics::dynamiceq tests\n");
    const double sr = 48000.0, pi = core::kPi;

    // --- below threshold + 0 dB static → passthrough ---
    test::group ("DynamicEqBand passthrough below threshold");
    {
        DP p; p.type = eq::FilterType::Bell; p.freq = 5000.0; p.Q = 2.0; p.mode = Mode::CutWhenLoud; p.thresholdDb = 12.0; p.ratio = 4.0;
        auto y = runTone (p, 5000.0, 0.3, 8000, sr);                        // −10 dBFS < 12 dB threshold → no cut
        std::vector<float> x (8000); for (int i = 0; i < 8000; ++i) x[i] = (float) (0.3 * std::sin (2.0 * pi * 5000.0 * i / sr));
        test::ok (std::fabs (rmsTail (y, 4000) / rmsTail (x, 4000) - 1.0) < 0.05, "0 dB bell below threshold → passthrough");
    }

    // --- cut when loud (de-ess) ---
    test::group ("DynamicEqBand cut-when-loud (de-ess)");
    {
        DP cut; cut.freq = 5000.0; cut.Q = 2.0; cut.mode = Mode::CutWhenLoud; cut.thresholdDb = -30.0; cut.ratio = 6.0; cut.rangeDb = 18.0;
        DP no = cut; no.thresholdDb = 12.0;                                 // same band, but never triggers
        const double rCut = rmsTail (runTone (cut, 5000.0, 0.7, 12000, sr), 6000);
        const double rNo  = rmsTail (runTone (no,  5000.0, 0.7, 12000, sr), 6000);
        test::ok (rCut < 0.6 * rNo, "a loud in-band tone is attenuated when over threshold");
    }

    // --- boost when quiet ---
    test::group ("DynamicEqBand boost-when-quiet");
    {
        DP boost; boost.freq = 1000.0; boost.Q = 2.0; boost.mode = Mode::BoostWhenQuiet; boost.thresholdDb = -12.0; boost.ratio = 4.0; boost.rangeDb = 12.0;
        DP flat = boost; flat.rangeDb = 0.0;                                // range 0 → no dynamic move
        const double rB = rmsTail (runTone (boost, 1000.0, 0.05, 12000, sr), 6000);   // −26 dBFS quiet → boosted
        const double rF = rmsTail (runTone (flat,  1000.0, 0.05, 12000, sr), 6000);
        test::ok (rB > 1.3 * rF, "a quiet in-band tone is lifted below threshold");
    }

    // --- frequency-selective: an off-band tone is untouched ---
    test::group ("DynamicEqBand off-band unaffected");
    {
        DP p; p.freq = 5000.0; p.Q = 2.0; p.mode = Mode::CutWhenLoud; p.thresholdDb = -30.0; p.ratio = 6.0; p.rangeDb = 18.0;
        auto y = runTone (p, 100.0, 0.7, 8000, sr);                         // loud 100 Hz — the 5 kHz sidechain BP doesn't see it
        std::vector<float> x (8000); for (int i = 0; i < 8000; ++i) x[i] = (float) (0.7 * std::sin (2.0 * pi * 100.0 * i / sr));
        test::ok (std::fabs (rmsTail (y, 4000) / rmsTail (x, 4000) - 1.0) < 0.05, "loud 100 Hz unaffected by a 5 kHz de-ess band");
    }

    // --- ratio=1 / range=0 → a plain static bell ---
    test::group ("DynamicEqBand static bell when not dynamic");
    {
        DP p; p.freq = 1000.0; p.Q = 2.0; p.staticGainDb = 6.0; p.ratio = 1.0; p.rangeDb = 0.0;
        auto y = runTone (p, 1000.0, 0.3, 8000, sr);
        std::vector<float> x (8000); for (int i = 0; i < 8000; ++i) x[i] = (float) (0.3 * std::sin (2.0 * pi * 1000.0 * i / sr));
        test::approx (20.0 * std::log10 (rmsTail (y, 4000) / rmsTail (x, 4000)), 6.0, 0.5, "static +6 dB bell boosts the band ~6 dB");
    }

    // --- no allocation in process() ---
    test::group ("DynamicEqBand no-alloc");
    {
        DB d; d.prepare (sr, 2); DP p; p.freq = 3000.0; p.mode = Mode::CutWhenLoud; p.thresholdDb = -20.0; d.setParams (p);
        std::vector<float> l (512, 0.3f), r (512, -0.2f); float* io[2] { l.data(), r.data() };
        d.process (io, 2, 512);
        const long before = g_allocs.load();
        d.process (io, 2, 512); d.process (io, 2, 512);
        test::okNoAlloc (g_allocs.load() == before, "process() did not allocate");
    }

    // --- stereo linked: one gain, image preserved ---
    test::group ("DynamicEqBand stereo identical");
    {
        DB d; d.prepare (sr, 2); DP p; p.freq = 5000.0; p.Q = 2.0; p.mode = Mode::CutWhenLoud; p.thresholdDb = -30.0; p.ratio = 6.0; d.setParams (p);
        const int N = 8000; std::vector<float> l (N), r (N);
        for (int i = 0; i < N; ++i) { const float v = (float) (0.7 * std::sin (2.0 * pi * 5000.0 * i / sr)); l[i] = v; r[i] = v; }
        for (int o = 0; o < N; o += 512) { float* io[2] { l.data() + o, r.data() + o }; d.process (io, 2, std::min (512, N - o)); }
        double md = 0; for (int i = 0; i < N; ++i) md = std::max (md, (double) std::fabs (l[i] - r[i]));
        test::ok (md == 0.0, "identical L/R in → identical L/R out (linked detector)");
    }

    // --- control-rate gain update is artifact-free (the top concern) ---
    test::group ("DynamicEqBand control-rate K artifact check");
    {
        auto run = [&] (int K)
        {
            DP p; p.freq = 5000.0; p.Q = 2.0; p.mode = Mode::CutWhenLoud; p.thresholdDb = -30.0; p.ratio = 6.0; p.rangeDb = 18.0; p.coeffUpdatePeriod = K;
            return runTone (p, 5000.0, 0.7, 12000, sr);
        };
        const auto k1 = run (1), k16 = run (16);
        double settled = 0; for (int i = 8000; i < 12000; ++i) settled = std::max (settled, (double) std::fabs (k16[i] - k1[i]));
        test::ok (settled < 1.0e-3, "K=16 converges to per-sample K=1 in steady state (no zipper)");
        auto maxStep = [] (const std::vector<float>& v) { double m = 0; for (int i = 1; i < (int) v.size(); ++i) m = std::max (m, (double) std::fabs (v[i] - v[i - 1])); return m; };
        test::ok (maxStep (k16) < 1.5 * maxStep (k1) + 1.0e-4, "no coeff-step click (K=16 max slope ≈ K=1 max slope)");
    }

    // --- ratio=1 / range=0 nulls SAMPLE-WISE against a plain eq::Svf Bell (same shape, no drift) ---
    test::group ("DynamicEqBand static null vs plain Svf Bell");
    {
        DP p; p.freq = 1200.0; p.Q = 1.5; p.staticGainDb = 5.0; p.ratio = 1.0; p.rangeDb = 0.0;
        DB d; d.prepare (sr, 1); d.setParams (p);
        eq::Svf ref; ref.prepare (sr, 1); ref.setParams (eq::FilterType::Bell, 1200.0, 1.5, 5.0);
        unsigned long s = 11; auto rng = [&]() { s = s * 6364136223846793005ULL + 1442695040888963407ULL; return (float) ((s >> 40) & 0xffff) / 32768.0f - 1.0f; };
        const int N = 6000; std::vector<float> x (N), y (N);
        for (int i = 0; i < N; ++i) { x[i] = 0.3f * rng(); y[i] = x[i]; }
        for (int o = 0; o < N; o += 512) { float* io[1] { y.data() + o }; d.process (io, 1, std::min (512, N - o)); }
        double md = 0; for (int i = 0; i < N; ++i) { const float r = ref.processSample (0, x[i]); if (i >= 1000) md = std::max (md, (double) std::fabs (y[i] - r)); }
        test::ok (md < 1.0e-3, "ratio=1 dynamic band == a plain Svf Bell (identical shape)");
    }

    return test::report();
}

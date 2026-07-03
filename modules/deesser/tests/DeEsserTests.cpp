// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// JUCE-free self-tests for the De-esser: both topologies attenuate a loud sibilant-band tone; an off-band
// (low) tone is untouched; the split-band sum is transparent (allpass reconstruction) when idle; `listen`
// solos the sidechain (passes sibilance, rejects lows); no alloc; stereo is linked.

#include <felitronics_test.h>
#include <felitronics/core/Math.h>
#include <felitronics/deesser/DeEsser.h>

#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <vector>

static std::atomic<long> g_allocs { 0 };
void* operator new      (std::size_t s) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s ? s : 1); }
void* operator new[]    (std::size_t s) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s ? s : 1); }
void  operator delete   (void* p) noexcept { std::free (p); }
void  operator delete[] (void* p) noexcept { std::free (p); }
void  operator delete   (void* p, std::size_t) noexcept { std::free (p); }
void  operator delete[] (void* p, std::size_t) noexcept { std::free (p); }

using namespace felitronics;
using DE   = deesser::DeEsser;
using DP   = deesser::DeEsserParams;
using Mode = deesser::DeEsserMode;

static double rmsTail (const std::vector<float>& v, int from)
{
    double e = 0.0; int n = 0;
    for (int i = from; i < (int) v.size(); ++i) { e += (double) v[i] * v[i]; ++n; }
    return n ? std::sqrt (e / (double) n) : 0.0;
}

static std::vector<float> runTone (const DP& p, double f, double amp, int N, double sr)
{
    DE d; d.prepare (sr, 512, 1); d.setParams (p);
    std::vector<float> y (N);
    for (int i = 0; i < N; ++i) y[i] = (float) (amp * std::sin (2.0 * core::kPi * f * i / sr));
    for (int o = 0; o < N; o += 512) { float* io[1] { y.data() + o }; d.process (io, 1, std::min (512, N - o)); }
    return y;
}

int main()
{
    std::printf ("felitronics::deesser tests\n");
    const double sr = 48000.0, pi = core::kPi;

    // --- split-band: a loud sibilant-band tone is ducked ---
    test::group ("DeEsser split-band de-ess");
    {
        DP cut; cut.mode = Mode::SplitBand; cut.fc = 6000.0; cut.scQ = 1.0; cut.thresholdDb = -30.0; cut.ratio = 4.0; cut.rangeDb = 18.0;
        DP no = cut; no.thresholdDb = 12.0;
        const double rCut = rmsTail (runTone (cut, 8000.0, 0.7, 12000, sr), 6000);
        const double rNo  = rmsTail (runTone (no,  8000.0, 0.7, 12000, sr), 6000);
        test::ok (rCut < 0.7 * rNo, "loud 8 kHz (above the 6 kHz split, in the detector band) is ducked");
    }

    // --- a low (off-band) tone is untouched ---
    test::group ("DeEsser off-band unaffected");
    {
        DP p; p.mode = Mode::SplitBand; p.fc = 6000.0; p.scQ = 1.0; p.thresholdDb = -30.0; p.ratio = 4.0; p.rangeDb = 18.0;
        auto y = runTone (p, 200.0, 0.7, 8000, sr);
        std::vector<float> x (8000); for (int i = 0; i < 8000; ++i) x[i] = (float) (0.7 * std::sin (2.0 * pi * 200.0 * i / sr));
        test::ok (std::fabs (rmsTail (y, 4000) / rmsTail (x, 4000) - 1.0) < 0.05, "low tone untouched (not detected; passes the low band)");
    }

    // --- the surgical dynamic-EQ mode also cuts the loud sibilant band ---
    test::group ("DeEsser dynamic-eq mode");
    {
        DP cut; cut.mode = Mode::DynamicEq; cut.fc = 7000.0; cut.scQ = 2.0; cut.thresholdDb = -30.0; cut.ratio = 4.0; cut.rangeDb = 18.0;
        DP no = cut; no.thresholdDb = 12.0;
        const double rCut = rmsTail (runTone (cut, 7000.0, 0.7, 12000, sr), 6000);
        const double rNo  = rmsTail (runTone (no,  7000.0, 0.7, 12000, sr), 6000);
        test::ok (rCut < 0.7 * rNo, "dynamic-eq mode cuts the loud sibilant band");
    }

    // --- split-band is transparent (allpass reconstruction) when idle ---
    test::group ("DeEsser idle transparent");
    {
        DP p; p.mode = Mode::SplitBand; p.fc = 6000.0; p.scQ = 1.0; p.thresholdDb = 12.0; p.ratio = 4.0;   // never triggers
        DE d; d.prepare (sr, 512, 1); d.setParams (p);
        unsigned long long s = 4; auto rng = [&]() { s = s * 6364136223846793005ULL + 1442695040888963407ULL; return (float) ((s >> 40) & 0xffff) / 32768.0f - 1.0f; };
        const int N = 8000; std::vector<float> x (N), y (N);
        for (int i = 0; i < N; ++i) { x[i] = 0.2f * rng(); y[i] = x[i]; }
        for (int o = 0; o < N; o += 512) { float* io[1] { y.data() + o }; d.process (io, 1, std::min (512, N - o)); }
        test::ok (std::fabs (rmsTail (y, 4000) / rmsTail (x, 4000) - 1.0) < 0.05, "idle split-band = allpass reconstruction (RMS flat)");
    }

    // --- listen solos the sidechain (passes sibilance, rejects lows) ---
    test::group ("DeEsser listen solos the sidechain");
    {
        DP p; p.mode = Mode::SplitBand; p.fc = 6000.0; p.scQ = 1.0; p.listen = true;
        const double r8   = rmsTail (runTone (p, 8000.0, 0.5, 4000, sr), 2000);
        const double r200 = rmsTail (runTone (p, 200.0,  0.5, 4000, sr), 2000);
        test::ok (r200 < 0.1 * r8, "listen: sidechain passes the sibilant band, rejects the lows");
    }

    // --- no allocation in process() (both modes) ---
    test::group ("DeEsser no-alloc");
    {
        DE d; d.prepare (sr, 512, 2); DP p; p.mode = Mode::SplitBand; p.fc = 6500.0; p.thresholdDb = -20.0; d.setParams (p);
        std::vector<float> l (512, 0.3f), r (512, -0.2f); float* io[2] { l.data(), r.data() };
        d.process (io, 2, 512);
        const long before = g_allocs.load();
        d.process (io, 2, 512);
        DP dq = p; dq.mode = Mode::DynamicEq; d.setParams (dq); d.process (io, 2, 512);
        test::okNoAlloc (g_allocs.load() == before, "process() did not allocate (split-band + dynamic-eq)");
    }

    // --- stereo linked: one gain, image preserved ---
    test::group ("DeEsser stereo identical");
    {
        DE d; d.prepare (sr, 512, 2); DP p; p.mode = Mode::SplitBand; p.fc = 6000.0; p.scQ = 1.0; p.thresholdDb = -30.0; p.ratio = 4.0; d.setParams (p);
        const int N = 8000; std::vector<float> l (N), r (N);
        for (int i = 0; i < N; ++i) { const float v = (float) (0.7 * std::sin (2.0 * pi * 8000.0 * i / sr)); l[i] = v; r[i] = v; }
        for (int o = 0; o < N; o += 512) { float* io[2] { l.data() + o, r.data() + o }; d.process (io, 2, std::min (512, N - o)); }
        double md = 0; for (int i = 0; i < N; ++i) md = std::max (md, (double) std::fabs (l[i] - r[i]));
        test::ok (md == 0.0, "identical L/R in → identical L/R out (linked detector)");
    }

    // --- toggling mode/listen mid-stream stays clean (the top worry here: stale-path state) ---
    test::group ("DeEsser mode/listen toggle continuity");
    {
        DE d; d.prepare (sr, 512, 1); DP p; p.mode = Mode::SplitBand; p.fc = 6000.0; p.scQ = 1.0; p.thresholdDb = -30.0; p.ratio = 4.0; d.setParams (p);
        const int N = 8000; std::vector<float> y (N); for (int i = 0; i < N; ++i) y[i] = (float) (0.5 * std::sin (2.0 * pi * 8000.0 * i / sr));
        for (int o = 0; o < N; o += 512)
        {
            if (o == 2048) { DP q = p; q.mode = Mode::DynamicEq; d.setParams (q); }                 // toggle topology
            if (o == 4608) { DP q = p; q.mode = Mode::DynamicEq; q.listen = true; d.setParams (q); } // toggle listen
            float* io[1] { y.data() + o }; d.process (io, 1, std::min (512, N - o));   // last block < 512 (N not a multiple)
        }
        bool bad = false; double mx = 0; for (int i = 0; i < N; ++i) { if (! std::isfinite (y[i])) bad = true; mx = std::max (mx, (double) std::fabs (y[i])); }
        test::ok (! bad && mx < 2.0, "mode/listen toggled mid-stream → no NaN, no full-scale spike");
    }

    // --- idle split-band nulls SAMPLE-WISE against a bare eq::Crossover2 reconstruction ---
    test::group ("DeEsser idle == Crossover2 reconstruction");
    {
        DP p; p.mode = Mode::SplitBand; p.fc = 6000.0; p.scQ = 1.0; p.thresholdDb = 12.0; p.ratio = 4.0;   // never triggers
        DE d; d.prepare (sr, 512, 1); d.setParams (p);
        eq::Crossover2 ref; ref.prepare (sr, 1); ref.setFrequency (6000.0f);
        unsigned long long s = 6; auto rng = [&]() { s = s * 6364136223846793005ULL + 1442695040888963407ULL; return (float) ((s >> 40) & 0xffff) / 32768.0f - 1.0f; };
        const int N = 6000; std::vector<float> x (N), y (N);
        for (int i = 0; i < N; ++i) { x[i] = 0.2f * rng(); y[i] = x[i]; }
        for (int o = 0; o < N; o += 512) { float* io[1] { y.data() + o }; d.process (io, 1, std::min (512, N - o)); }
        double md = 0; for (int i = 0; i < N; ++i) { float lo, hi; ref.processSample (0, x[i], lo, hi); if (i >= 1000) md = std::max (md, (double) std::fabs (y[i] - (lo + hi))); }
        test::ok (md < 1e-3, "idle split-band == Crossover2 low+high (allpass reconstruction)");
    }

    // --- the DEFAULT (dynamic-EQ) mode is FREQUENCY-FOCUSED: the bell cuts at fc, far air is spared.
    //     (SplitBand intentionally ducks the WHOLE high band — the documented "can dull air" trade-off.) ---
    test::group ("DeEsser dynamic-eq de-ess is frequency-focused");
    {
        DP p; p.mode = Mode::DynamicEq; p.fc = 6000.0; p.scQ = 2.0; p.thresholdDb = -25.0; p.ratio = 4.0; p.rangeDb = 18.0;
        auto reduction = [&] (double f) -> double
        {
            const auto y = runTone (p, f, 0.6, 10000, sr);
            std::vector<float> x (10000); for (int i = 0; i < 10000; ++i) x[i] = (float) (0.6 * std::sin (2.0 * pi * f * i / sr));
            return rmsTail (y, 5000) / rmsTail (x, 5000);                  // < 1 → cut
        };
        test::ok (reduction (6000.0) < 0.9 * reduction (13000.0), "the 6 kHz bell cuts its band; 13 kHz air (off the bell) is spared");
    }

    // --- non-finite params must not poison the stream (NaN fc reaches the Svf/Crossover otherwise) ---
    test::group ("DeEsser non-finite params rejected");
    {
        const double qnan = std::numeric_limits<double>::quiet_NaN();
        for (auto mode : { Mode::DynamicEq, Mode::SplitBand })
        {
            DP p; p.mode = mode; p.fc = qnan; p.scQ = qnan; p.thresholdDb = qnan;
            p.ratio = qnan; p.kneeDb = qnan; p.rangeDb = qnan; p.attackMs = qnan; p.releaseMs = qnan;
            DE d; d.prepare (sr, 512, 1); d.setParams (p);
            const int N = 4000; std::vector<float> y (N);
            for (int i = 0; i < N; ++i) y[i] = (float) (0.5 * std::sin (2.0 * pi * 8000.0 * i / sr));
            for (int o = 0; o < N; o += 512) { float* io[1] { y.data() + o }; d.process (io, 1, std::min (512, N - o)); }
            bool fin = true; for (float v : y) fin &= (bool) std::isfinite (v);
            test::ok (fin, mode == Mode::DynamicEq ? "all-NaN params → finite output (DynamicEq mode)"
                                                   : "all-NaN params → finite output (SplitBand mode)");
        }
    }

    return test::report();
}

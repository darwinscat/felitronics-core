// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// JUCE-free self-tests for the TransientShaper: a STEADY tone is left alone (the fast/slow detectors share a
// type so norm≈0 in steady state — different crest factors would falsely fire); attack shaping boosts/cuts
// the onset monotonically with attackDb; sustain shaping changes the decay tail; mix=0 is exact dry; latency
// follows the lookahead; process() never allocates; a linked gain preserves the stereo image.

#include <felitronics_test.h>
#include <felitronics/core/Math.h>
#include <felitronics/dynamics/TransientShaper.h>

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
using TS = dynamics::TransientShaper;
using TP = dynamics::TransientShaperParams;

static double rmsRange (const std::vector<float>& v, int a, int b)
{
    double e = 0.0; int n = 0;
    for (int i = a; i < b && i < (int) v.size(); ++i) { e += (double) v[i] * v[i]; ++n; }
    return n ? std::sqrt (e / (double) n) : 0.0;
}

int main()
{
    std::printf ("felitronics::dynamics TransientShaper tests\n");
    const double sr = 48000.0, pi = core::kPi;

    // --- a STEADY tone is not shaped (the both-Peak detector reads norm ≈ 0 in steady state) ---
    test::group ("TransientShaper steady tone unaffected");
    {
        TS ts; ts.prepare (sr, 512, 1);
        TP p; p.attackDb = 12.0; p.sustainDb = 12.0; ts.setParams (p);
        const int N = 12000; std::vector<float> x (N), y (N);
        for (int i = 0; i < N; ++i) { x[i] = 0.5f * (float) std::sin (2.0 * pi * 1000.0 * i / sr); y[i] = x[i]; }
        float* io[1] { y.data() }; ts.process (io, 1, N);
        test::ok (std::fabs (rmsRange (y, 6000, 12000) / rmsRange (x, 6000, 12000) - 1.0) < 0.05,
                  "steady sine: output ≈ input (no false transient shaping)");
    }

    // --- attack shaping is monotonic in attackDb ---
    test::group ("TransientShaper attack shaping");
    {
        auto run = [&] (double atkDb) -> double
        {
            TS ts; ts.prepare (sr, 512, 1); TP p; p.attackDb = atkDb; ts.setParams (p);
            const int N = 12000; std::vector<float> y (N, 0.0f);
            const int onset = 48;                                            // 1 ms of silence, then a tone onset
            for (int i = onset; i < N; ++i) y[i] = 0.5f * (float) std::sin (2.0 * pi * 1000.0 * (i - onset) / sr);
            float* io[1] { y.data() }; ts.process (io, 1, N);
            return rmsRange (y, onset, onset + 480) / rmsRange (y, 9000, 12000);   // onset region / steady region
        };
        const double up = run (12.0), flat = run (0.0), down = run (-12.0);
        test::ok (up > 1.1,   "attackDb +12 lifts the onset above the steady level");
        test::ok (down < 0.9, "attackDb -12 softens the onset");
        test::ok (up > flat && flat > down, "attack shaping monotonic in attackDb");
    }

    // --- sustain shaping changes the decay tail ---
    test::group ("TransientShaper sustain shaping");
    {
        auto tail = [&] (double susDb) -> double
        {
            TS ts; ts.prepare (sr, 512, 1); TP p; p.sustainDb = susDb; ts.setParams (p);
            const int N = 12000; std::vector<float> y (N);
            for (int i = 0; i < N; ++i) { const double t = (double) i / sr; y[i] = (float) (std::exp (-t / 0.03) * std::sin (2.0 * pi * 1000.0 * i / sr)); }
            float* io[1] { y.data() }; ts.process (io, 1, N);
            return rmsRange (y, 4800, 9600);                                 // the 100–200 ms decay tail
        };
        test::ok (tail (12.0) > tail (-12.0) * 1.2, "sustainDb +12 fattens the tail vs -12");
    }

    // --- mix = 0 is exact dry (no lookahead) ---
    test::group ("TransientShaper mix");
    {
        TS ts; ts.prepare (sr, 512, 1); TP p; p.attackDb = 12.0; p.mix = 0.0; ts.setParams (p);
        const int N = 4000; std::vector<float> x (N, 0.0f), y (N, 0.0f);
        const int onset = 48; for (int i = onset; i < N; ++i) { x[i] = 0.5f * (float) std::sin (2.0 * pi * 1000.0 * (i - onset) / sr); y[i] = x[i]; }
        float* io[1] { y.data() }; ts.process (io, 1, N);
        double md = 0; for (int i = 0; i < N; ++i) md = std::max (md, (double) std::fabs (y[i] - x[i]));
        test::ok (md < 1e-6, "mix=0 → exact dry passthrough");
    }

    // --- zero latency (gain aligned with audio; no anticipatory lookahead in v1) ---
    test::group ("TransientShaper zero latency");
    {
        test::ok (TS::latencySamples() == 0, "zero latency");
    }

    // --- no allocation in process() ---
    test::group ("TransientShaper no-alloc");
    {
        TS ts; ts.prepare (sr, 512, 2); TP p; p.attackDb = 6.0; p.sustainDb = -6.0; ts.setParams (p);
        std::vector<float> l (512, 0.2f), r (512, -0.1f); float* io[2] { l.data(), r.data() };
        ts.process (io, 2, 512);
        const long before = g_allocs.load();
        ts.process (io, 2, 512); ts.process (io, 2, 512);
        test::ok (g_allocs.load() == before, "process() did not allocate");
    }

    // --- linked gain preserves the stereo image ---
    test::group ("TransientShaper stereo image preserved");
    {
        TS ts; ts.prepare (sr, 512, 2); TP p; p.attackDb = 12.0; ts.setParams (p);
        const int N = 4000; std::vector<float> l (N), r (N);
        const int onset = 48; for (int i = 0; i < N; ++i) { const float v = (i < onset) ? 0.0f : 0.5f * (float) std::sin (2.0 * pi * 1000.0 * (i - onset) / sr); l[i] = v; r[i] = v; }
        float* io[2] { l.data(), r.data() }; ts.process (io, 2, N);
        double md = 0; for (int i = 0; i < N; ++i) md = std::max (md, (double) std::fabs (l[i] - r[i]));
        test::ok (md == 0.0, "identical L/R in → identical L/R out (one linked gain)");
    }

    // --- steady tones across the spectrum stay unshaped (the deadzone holds at low freq too) ---
    test::group ("TransientShaper steady unshaped across freq");
    {
        double worst = 0.0;
        for (double f : { 50.0, 500.0, 5000.0 })
        {
            TS ts; ts.prepare (sr, 512, 1); TP p; p.attackDb = 24.0; p.sustainDb = 24.0; ts.setParams (p);
            const int N = 12000; std::vector<float> x (N), y (N);
            for (int i = 0; i < N; ++i) { x[i] = 0.5f * (float) std::sin (2.0 * pi * f * i / sr); y[i] = x[i]; }
            float* io[1] { y.data() }; ts.process (io, 1, N);
            worst = std::max (worst, std::fabs (rmsRange (y, 6000, 12000) / rmsRange (x, 6000, 12000) - 1.0));
        }
        test::ok (worst < 0.05, "50 / 500 / 5k Hz steady tones unshaped (deadzone holds)");
    }

    // --- attack and sustain controls are INDEPENDENT (attack→onset, sustain→tail) ---
    test::group ("TransientShaper attack/sustain independence");
    {
        auto run = [&] (double atk, double sus, int a, int b) -> double
        {
            TS ts; ts.prepare (sr, 512, 1); TP p; p.attackDb = atk; p.sustainDb = sus; ts.setParams (p);
            const int N = 12000; std::vector<float> y (N);
            for (int i = 0; i < N; ++i) { const double t = (double) i / sr; y[i] = (float) (std::exp (-t / 0.04) * std::sin (2.0 * pi * 1000.0 * i / sr)); }
            float* io[1] { y.data() }; ts.process (io, 1, N);
            return rmsRange (y, a, b);
        };
        const double onsetFlat = run (0, 0, 0, 480), tailFlat = run (0, 0, 6000, 11000);
        test::ok (run (12, 0, 0, 480)   > onsetFlat * 1.1, "attackDb lifts the onset");
        test::ok (run (0, 12, 6000, 11000) > tailFlat * 1.1, "sustainDb lifts the tail");
        test::ok (std::fabs (run (12, 0, 6000, 11000) / tailFlat  - 1.0) < 0.15, "attackDb leaves the tail ~unchanged");
        test::ok (std::fabs (run (0, 12, 0, 480)       / onsetFlat - 1.0) < 0.15, "sustainDb leaves the onset ~unchanged");
    }

    // --- threshold: higher deadzone → less steady-state shaping ---
    test::group ("TransientShaper threshold");
    {
        auto ratio = [&] (double thr) -> double
        {
            TS ts; ts.prepare (sr, 512, 1); TP p; p.attackDb = 12.0; p.sustainDb = 12.0; p.threshold = thr; ts.setParams (p);
            const int N = 12000; std::vector<float> x (N), y (N);
            for (int i = 0; i < N; ++i) { x[i] = 0.5f * (float) std::sin (2.0 * pi * 200.0 * i / sr); y[i] = x[i]; }
            float* io[1] { y.data() }; ts.process (io, 1, N);
            return rmsRange (y, 6000, 12000) / rmsRange (x, 6000, 12000);
        };
        test::ok (std::fabs (ratio (0.9) - 1.0) <= std::fabs (ratio (0.0) - 1.0) + 1e-6, "threshold 0.9 leaves the steady tone at least as unshaped as 0.0");
    }

    // --- sample-wise null against a hand-recomputed deadzone+smooth+mix pipeline (the key check) ---
    test::group ("TransientShaper null vs hand-computed gain");
    {
        TS ts; ts.prepare (sr, 512, 1); TP p; p.attackDb = 9.0; p.sustainDb = -6.0; p.threshold = 0.2; p.gainSmoothMs = 1.0; p.mix = 0.8; ts.setParams (p);
        dynamics::EnvelopeFollower rf, rs;
        rf.prepare (sr); rf.setDetector (dynamics::Detector::Peak); rf.setTimes (0.3, 20.0);
        rs.prepare (sr); rs.setDetector (dynamics::Detector::Peak); rs.setTimes (15.0, 150.0);
        const float thr = 0.2f, attDb = 9.0f, susDb = -6.0f, mix = 0.8f, sc = (float) std::exp (-1.0 / (0.001 * sr));
        float gsm = 1.0f;
        const int N = 4000; std::vector<float> x (N), y (N);
        for (int i = 0; i < N; ++i) { const float v = (i < 48) ? 0.0f : 0.5f * (float) std::sin (2.0 * pi * 1000.0 * (i - 48) / sr); x[i] = v; y[i] = v; }
        float* io[1] { y.data() }; ts.process (io, 1, N);
        double md = 0;
        for (int i = 0; i < N; ++i)
        {
            const float lk = std::fabs (x[i]); const float fe = rf.process (lk), se = rs.process (lk);
            const float norm = (fe - se) / (fe + se + 1.0e-9f); const float a = std::fabs (norm);
            const float t = a <= thr ? 0.0f : (a - thr) / (1.0f - thr); const float sh = norm < 0.0f ? -t : t;
            float gdb = std::clamp (attDb * std::max (sh, 0.0f) + susDb * std::max (-sh, 0.0f), -24.0f, 24.0f);
            const float g = (float) core::dbToGain ((double) gdb); gsm = g + sc * (gsm - g);
            md = std::max (md, (double) std::fabs (y[i] - x[i] * (1.0f + mix * (gsm - 1.0f))));
        }
        test::ok (md < 1e-5, "output == hand-computed deadzone→smooth→mix pipeline (sample-wise null)");
    }

    return test::report();
}

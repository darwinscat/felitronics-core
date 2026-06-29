// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa. Part of felitronics-core — see LICENSE.

// JUCE-free self-tests for the stereo module: the Mid/Side matrix round-trips, and MonoBass collapses the
// low end to mono while (a) passing IN-PHASE bass bit-exact (the kept mono content is never filtered),
// (b) removing low anti-phase side energy, (c) preserving high anti-phase side energy, with clean bypass.

#include <felitronics_test.h>
#include <felitronics/core/Math.h>
#include <felitronics/stereo/MidSide.h>
#include <felitronics/stereo/MonoBass.h>

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
    double e = 0.0; int c = 0;
    for (int i = from; i < (int) v.size(); ++i) { e += (double) v[i] * v[i]; ++c; }
    return c ? std::sqrt (e / (double) c) : 0.0;
}

int main()
{
    std::printf ("felitronics::stereo tests\n");
    const double sr = 48000.0, pi = core::kPi;

    // --- Mid/Side round-trip ---
    test::group ("MidSide identity");
    {
        float m, s, l, r;
        stereo::MidSide::encode (0.3f, -0.7f, m, s);
        stereo::MidSide::decode (m, s, l, r);
        test::approx (l,  0.3, 1e-6, "decode(encode) recovers L");
        test::approx (r, -0.7, 1e-6, "decode(encode) recovers R");
        test::approx (m, -0.2, 1e-6, "M = 1/2(L+R)");
        test::approx (s,  0.5, 1e-6, "S = 1/2(L-R)");
    }

    // --- the headline property: IN-PHASE bass is NOT filtered (vs a full crossover, which would) ---
    test::group ("MonoBass: mono content passes UNFILTERED");
    {
        stereo::MonoBass mb; mb.prepare (sr); mb.setFrequency (150.0f); mb.setLowWidth (0.0f);
        const int N = 4000; std::vector<float> L (N), R (N), ref (N);
        for (int i = 0; i < N; ++i) { const float v = 0.5f * (float) std::sin (2.0 * pi * 50.0 * i / sr); L[i] = v; R[i] = v; ref[i] = v; }
        float* io[2] { L.data(), R.data() }; mb.process (io, 2, N);
        double md = 0; for (int i = 0; i < N; ++i) { md = std::max (md, (double) std::fabs (L[i] - ref[i])); md = std::max (md, (double) std::fabs (R[i] - ref[i])); }
        test::ok (md < 1e-6, "in-phase 50 Hz (L==R) passes bit-exact — kept bass never crossover-filtered");
    }

    // --- low anti-phase side is removed (monoed) ---
    test::group ("MonoBass: low anti-phase side removed");
    {
        stereo::MonoBass mb; mb.prepare (sr); mb.setFrequency (150.0f); mb.setLowWidth (0.0f);
        const int N = 8000; std::vector<float> L (N), R (N); double inSq = 0.0;
        for (int i = 0; i < N; ++i) { const float v = 0.5f * (float) std::sin (2.0 * pi * 50.0 * i / sr); L[i] = v; R[i] = -v; inSq += (double) v * v; }
        const double inRms = std::sqrt (inSq / N);
        float* io[2] { L.data(), R.data() }; mb.process (io, 2, N);
        test::ok (rmsTail (L, N / 2) < 0.1 * inRms, "50 Hz anti-phase side strongly attenuated (collapsed to mono)");
    }

    // --- high anti-phase side is preserved (stereo above the cutoff) ---
    test::group ("MonoBass: high anti-phase side preserved");
    {
        stereo::MonoBass mb; mb.prepare (sr); mb.setFrequency (150.0f); mb.setLowWidth (0.0f);
        const int N = 8000; std::vector<float> L (N), R (N); double inSq = 0.0;
        for (int i = 0; i < N; ++i) { const float v = 0.5f * (float) std::sin (2.0 * pi * 2000.0 * i / sr); L[i] = v; R[i] = -v; inSq += (double) v * v; }
        const double inRms = std::sqrt (inSq / N);
        float* io[2] { L.data(), R.data() }; mb.process (io, 2, N);
        test::ok (rmsTail (L, N / 2) > 0.9 * inRms, "2 kHz anti-phase side preserved (above the cutoff)");
    }

    // --- bypass paths are bit-exact ---
    test::group ("MonoBass: bypass");
    {
        stereo::MonoBass mb; mb.prepare (sr); mb.setFrequency (150.0f);
        const int N = 256; std::vector<float> L (N), R (N), L0 (N), R0 (N);
        for (int i = 0; i < N; ++i) { L[i] = 0.3f * (float) std::sin (2.0 * pi * 100.0 * i / sr); R[i] = -0.4f * (float) std::sin (2.0 * pi * 300.0 * i / sr); L0[i] = L[i]; R0[i] = R[i]; }
        float* io[2] { L.data(), R.data() };

        mb.setLowWidth (1.0f); mb.process (io, 2, N);             // fully wide -> bypass
        double md = 0; for (int i = 0; i < N; ++i) { md = std::max (md, (double) std::fabs (L[i] - L0[i])); md = std::max (md, (double) std::fabs (R[i] - R0[i])); }
        test::ok (md == 0.0, "lowWidth=1 -> exact passthrough");

        mb.setLowWidth (0.0f); mb.setEnabled (false); mb.process (io, 2, N);
        md = 0; for (int i = 0; i < N; ++i) { md = std::max (md, (double) std::fabs (L[i] - L0[i])); md = std::max (md, (double) std::fabs (R[i] - R0[i])); }
        test::ok (md == 0.0, "disabled -> exact passthrough");
    }

    // --- no allocation in process() ---
    test::group ("MonoBass: no-alloc");
    {
        stereo::MonoBass mb; mb.prepare (sr); mb.setFrequency (150.0f); mb.setLowWidth (0.0f);
        std::vector<float> L (512, 0.2f), R (512, -0.2f); float* io[2] { L.data(), R.data() };
        mb.process (io, 2, 512);
        const long before = g_allocs.load();
        mb.process (io, 2, 512); mb.process (io, 2, 512);
        test::ok (g_allocs.load() == before, "process() did not allocate");
    }

    return test::report();
}

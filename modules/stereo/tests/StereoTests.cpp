// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// JUCE-free self-tests for the stereo module: the Mid/Side matrix round-trips, and MonoBass collapses the
// low end to mono while (a) passing IN-PHASE bass bit-exact (the kept mono content is never filtered),
// (b) removing low anti-phase side energy, (c) preserving high anti-phase side energy, with clean bypass.

#include <felitronics_test.h>
#include <felitronics/core/Math.h>
#include <felitronics/stereo/MidSide.h>
#include <felitronics/stereo/MonoBass.h>
#include <felitronics/stereo/StereoWidth.h>

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
        test::okNoAlloc (g_allocs.load() == before, "process() did not allocate");
    }

    //==========================================================================
    // StereoWidth: a Mid/Side side-gain. The headline safety property is the mono-fold invariant.
    auto rngPair = [] (int N, unsigned long seed, std::vector<float>& L, std::vector<float>& R)
    {
        unsigned long s = seed; auto u = [&] { s = s * 6364136223846793005ULL + 1442695040888963407ULL; return (float) ((s >> 40) & 0xffff) / 32768.0f - 1.0f; };
        L.resize (N); R.resize (N); for (int i = 0; i < N; ++i) { L[i] = 0.6f * u(); R[i] = 0.6f * u(); }
    };

    // --- THE mono-fold invariant: ½(L'+R') == ½(L+R) for ANY width (we touch only the Side) ---
    test::group ("StereoWidth: mono-fold invariant");
    {
        const int N = 4000; std::vector<float> L0, R0; rngPair (N, 7, L0, R0);
        double worst = 0.0;
        for (float w : { 0.0f, 0.5f, 1.3f, 2.0f })
        {
            std::vector<float> L = L0, R = R0; stereo::StereoWidth sw; sw.prepare (sr); sw.setWidth (w);
            float* io[2] { L.data(), R.data() }; sw.process (io, 2, N);
            for (int i = 0; i < N; ++i) worst = std::max (worst, (double) std::fabs (0.5f * (L[i] + R[i]) - 0.5f * (L0[i] + R0[i])));
        }
        test::ok (worst < 1e-6, "mono sum unchanged at width ∈ {0, .5, 1.3, 2} (widening can't weaken the fold)");
    }

    // --- neutral (width=1, gain=1) is BIT-EXACT passthrough (the early-return bypass) ---
    test::group ("StereoWidth: neutral is bit-exact");
    {
        const int N = 512; std::vector<float> L0, R0; rngPair (N, 3, L0, R0); auto L = L0, R = R0;
        stereo::StereoWidth sw; sw.prepare (sr);                              // defaults: width 1, gain 1
        float* io[2] { L.data(), R.data() }; sw.process (io, 2, N);
        double md = 0; for (int i = 0; i < N; ++i) { md = std::max (md, (double) std::fabs (L[i] - L0[i])); md = std::max (md, (double) std::fabs (R[i] - R0[i])); }
        test::ok (md == 0.0, "width=1, gain=1 → exact passthrough");
    }

    // --- width=0 collapses to the mono sum (L'==R'==M), exactly ---
    test::group ("StereoWidth: width=0 → mono");
    {
        const int N = 1000; std::vector<float> L0, R0; rngPair (N, 9, L0, R0); auto L = L0, R = R0;
        stereo::StereoWidth sw; sw.prepare (sr); sw.setWidth (0.0f); sw.reset();   // snap past the smoothing ramp
        float* io[2] { L.data(), R.data() }; sw.process (io, 2, N);
        double diff = 0, toM = 0; for (int i = 0; i < N; ++i) { diff = std::max (diff, (double) std::fabs (L[i] - R[i])); toM = std::max (toM, (double) std::fabs (L[i] - 0.5f * (L0[i] + R0[i]))); }
        test::ok (diff == 0.0 && toM < 1e-6, "width=0 → L==R==½(L+R) (pure mono)");
    }

    // --- width=2 doubles the side: (L'-R') == 2·(L-R) ---
    test::group ("StereoWidth: width=2 doubles the side");
    {
        const int N = 1000; std::vector<float> L0, R0; rngPair (N, 21, L0, R0); auto L = L0, R = R0;
        stereo::StereoWidth sw; sw.prepare (sr); sw.setWidth (2.0f); sw.reset();   // snap past the smoothing ramp
        float* io[2] { L.data(), R.data() }; sw.process (io, 2, N);
        double md = 0; for (int i = 0; i < N; ++i) md = std::max (md, (double) std::fabs ((L[i] - R[i]) - 2.0f * (L0[i] - R0[i])));
        test::ok (md < 1e-5, "side difference scaled ×2 (the S axis tracks width)");
    }

    // --- a MONO source stays bit-exact for any width (you cannot synthesise stereo from S=0) ---
    test::group ("StereoWidth: mono source is invariant");
    {
        const int N = 600; std::vector<float> in; std::vector<float> dummy; rngPair (N, 5, in, dummy);
        auto L = in, R = in;                                                  // L==R → S=0
        stereo::StereoWidth sw; sw.prepare (sr); sw.setWidth (2.0f);
        float* io[2] { L.data(), R.data() }; sw.process (io, 2, N);
        double md = 0; for (int i = 0; i < N; ++i) { md = std::max (md, (double) std::fabs (L[i] - in[i])); md = std::max (md, (double) std::fabs (R[i] - in[i])); }
        test::ok (md == 0.0, "L==R in → unchanged out at width=2 (width scales zero side)");
    }

    // --- outputGain is a clean linear trim on the whole signal ---
    test::group ("StereoWidth: outputGain trims level");
    {
        const int N = 600; std::vector<float> L0, R0; rngPair (N, 13, L0, R0); auto L = L0, R = R0;
        stereo::StereoWidth sw; sw.prepare (sr); sw.setWidth (1.0f); sw.setOutputGain (2.0f); sw.reset();
        float* io[2] { L.data(), R.data() }; sw.process (io, 2, N);
        double md = 0; for (int i = 0; i < N; ++i) { md = std::max (md, (double) std::fabs (L[i] - 2.0f * L0[i])); md = std::max (md, (double) std::fabs (R[i] - 2.0f * R0[i])); }
        test::ok (md < 1e-6, "width=1, gain=2 → output = 2× input (M and S both scale)");
    }

    // --- side energy tracks width monotonically (and ≈ proportional) ---
    test::group ("StereoWidth: side energy ∝ width");
    {
        const int N = 4000; std::vector<float> L0, R0; rngPair (N, 31, L0, R0);
        auto sideRms = [&] (float w) {
            std::vector<float> L = L0, R = R0; stereo::StereoWidth sw; sw.prepare (sr); sw.setWidth (w); sw.reset();
            float* io[2] { L.data(), R.data() }; sw.process (io, 2, N);
            std::vector<float> sde (N); for (int i = 0; i < N; ++i) sde[i] = 0.5f * (L[i] - R[i]);
            return rmsTail (sde, 0);
        };
        const double s05 = sideRms (0.5f), s10 = sideRms (1.0f), s15 = sideRms (1.5f);
        test::ok (s05 < s10 && s10 < s15 && std::fabs (s15 / s10 - 1.5) < 0.02, "narrower→less / wider→more side, ratio ≈ width");
    }

    // --- bypass paths: < 2 channels and disabled are untouched ---
    test::group ("StereoWidth: bypass");
    {
        const int N = 256; std::vector<float> L0, R0; rngPair (N, 99, L0, R0); auto L = L0, R = R0;
        stereo::StereoWidth sw; sw.prepare (sr); sw.setWidth (1.8f);
        float* mono[1] { L.data() }; sw.process (mono, 1, N);                 // < 2 ch → passthrough
        double md = 0; for (int i = 0; i < N; ++i) md = std::max (md, (double) std::fabs (L[i] - L0[i]));
        test::ok (md == 0.0, "mono call (numChannels<2) → untouched");
        sw.setEnabled (false); float* io[2] { L.data(), R.data() }; sw.process (io, 2, N);
        md = 0; for (int i = 0; i < N; ++i) { md = std::max (md, (double) std::fabs (L[i] - L0[i])); md = std::max (md, (double) std::fabs (R[i] - R0[i])); }
        test::ok (md == 0.0, "disabled → untouched");
    }

    // --- no allocation in process() ---
    test::group ("StereoWidth: no-alloc");
    {
        stereo::StereoWidth sw; sw.prepare (sr); sw.setWidth (1.7f); sw.setOutputGain (1.1f);
        std::vector<float> L (512, 0.3f), R (512, -0.2f); float* io[2] { L.data(), R.data() };
        sw.process (io, 2, 512);
        const long before = g_allocs.load();
        sw.process (io, 2, 512); sw.process (io, 2, 512);
        test::okNoAlloc (g_allocs.load() == before, "process() did not allocate");
    }

    // --- a live width step is SMOOTHED (no click) — and the mono fold stays invariant through the ramp ---
    test::group ("StereoWidth: width step is click-free");
    {
        const int N = 4000; std::vector<float> L0 (N), R0 (N);
        for (int i = 0; i < N; ++i) { L0[i] = 0.5f * (float) std::sin (2.0 * pi * 500.0 * i / sr); R0[i] = 0.5f * (float) std::sin (2.0 * pi * 500.0 * i / sr + 0.6); }
        double inSlew = 0; for (int i = 1; i < N; ++i) inSlew = std::max (inSlew, (double) std::fabs (0.5f * (L0[i] - R0[i]) - 0.5f * (L0[i - 1] - R0[i - 1])));
        auto L = L0, R = R0; stereo::StereoWidth sw; sw.prepare (sr); sw.reset();        // settled at width 1
        float* a[2] { L.data(), R.data() }; sw.process (a, 2, N / 2);                     // width 1 (near-passthrough)
        sw.setWidth (2.0f);                                                              // STEP — must ramp, not jump
        float* b[2] { L.data() + N / 2, R.data() + N / 2 }; sw.process (b, 2, N / 2);
        double worstStep = 0, worstMono = 0;
        for (int i = 1; i < N; ++i)
        {
            worstStep = std::max (worstStep, (double) std::fabs (0.5f * (L[i] - R[i]) - 0.5f * (L[i - 1] - R[i - 1])));
            worstMono = std::max (worstMono, (double) std::fabs (0.5f * (L[i] + R[i]) - 0.5f * (L0[i] + R0[i])));
        }
        test::ok (worstStep < 4.0 * inSlew, "width 1→2 ramps smoothly (side step ≪ the instant 2× jump an unsmoothed step gives)");
        test::ok (worstMono < 1e-6, "mono fold stays invariant THROUGH the width ramp (Mid untouched per sample)");
    }

    // --- a stray non-finite parameter is rejected, never poisons the stream ---
    test::group ("StereoWidth: non-finite params rejected");
    {
        const int N = 256; std::vector<float> L0, R0; rngPair (N, 55, L0, R0); auto L = L0, R = R0;
        stereo::StereoWidth sw; sw.prepare (sr); sw.setWidth (1.5f); sw.setOutputGain (1.2f); sw.reset();
        sw.setWidth (std::nanf ("")); sw.setOutputGain (INFINITY); sw.setWidth (-INFINITY);   // all must be ignored
        test::ok (sw.width() == 1.5f && sw.outputGain() == 1.2f, "setWidth(NaN)/setOutputGain(±inf) ignored — last good value kept");
        float* io[2] { L.data(), R.data() }; sw.process (io, 2, N);
        bool fin = true; for (int i = 0; i < N; ++i) if (! std::isfinite (L[i]) || ! std::isfinite (R[i])) fin = false;
        test::ok (fin, "output stays finite after a NaN/inf parameter poke");
    }

    // --- aliased channels (io[0] == io[1], dual-mono) apply gain ONCE, not squared ---
    test::group ("StereoWidth: aliased L==R applies gain once");
    {
        const int N = 200; std::vector<float> in, dummy; rngPair (N, 77, in, dummy); auto buf = in;
        stereo::StereoWidth sw; sw.prepare (sr); sw.setWidth (1.7f); sw.setOutputGain (2.0f); sw.reset();
        float* io[2] { buf.data(), buf.data() };                                  // L and R alias one buffer
        sw.process (io, 2, N);
        double md = 0; for (int i = 0; i < N; ++i) md = std::max (md, (double) std::fabs (buf[i] - 2.0f * in[i]));
        test::ok (md < 1e-5, "io[0]==io[1] → output = gain·in (gain applied once; width inert on S=0)");
    }

    return test::report();
}

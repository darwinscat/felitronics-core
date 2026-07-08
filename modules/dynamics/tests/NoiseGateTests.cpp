// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// JUCE-free self-tests for felitronics::dynamics::NoiseGate — a Schmitt-trigger noise gate (distinct from
// the continuous Compressor). This battery is ADVERSARIAL: several cases encode bugs an earlier design had
// (a permanent fail-OPEN on a finite ±FLT_MAX overflow, a dB-linear open that swallowed the pick attack, a
// GR reading that lied while disabled, low-note chatter) plus the felitronics discipline: block-split NULL,
// SR-invariance, and no-allocation-in-process(). Never loosen these to go green.

#include <felitronics_test.h>
#include <felitronics/dynamics/NoiseGate.h>

#include <atomic>
#include <cfloat>
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
using dynamics::NoiseGate;

static constexpr double kSr = 48000.0;
static constexpr double kTwoPi = 6.283185307179586;

static std::vector<float> sine (double f, double sr, int n, double peak)
{
    std::vector<float> b ((std::size_t) n);
    for (int i = 0; i < n; ++i) b[(std::size_t) i] = (float) (peak * std::sin (kTwoPi * f * i / sr));
    return b;
}

// Run `blocks` of a looped mono `src` through analyse(); return the core VCA gain at the last block.
static float runMono (NoiseGate& g, const std::vector<float>& src, int blocks, int n, bool on, float thr)
{
    std::vector<float> blk ((std::size_t) n);
    for (int b = 0; b < blocks; ++b)
    {
        for (int i = 0; i < n; ++i) blk[(std::size_t) i] = src[(std::size_t) ((b * n + i) % (int) src.size())];
        const float* in[1] { blk.data() };
        g.analyse (in, 1, n, on, thr);
    }
    return g.currentCoreGain();
}

int main()
{
    std::printf ("felitronics::dynamics NoiseGate tests\n");
    const int N = 512;

    test::group ("starts CLOSED — silence + on holds the floor (no startup leak)");
    {
        NoiseGate g; g.prepare (kSr, N, 2); g.seedEnabled (true);   // a session restored gate-ON: no fade-in leak
        std::vector<float> z ((std::size_t) N, 0.0f); const float* in[1] { z.data() };
        g.analyse (in, 1, N, true, -60.0f);
        test::ok (g.currentCoreGain() < 1.0e-3f, "coreGain at floor on silence");
        test::ok (g.currentGain()     < 1.0e-3f, "effective gain at floor (enabled)");
    }

    test::group ("gate OFF = bit-exact unity passthrough");
    {
        NoiseGate g; g.prepare (kSr, N, 2);
        auto s = sine (220.0, kSr, N, 0.3);
        std::vector<float> probe ((std::size_t) N, 1.0f);
        for (int b = 0; b < 8; ++b) { const float* in[1] { s.data() }; g.analyse (in, 1, N, false, -60.0f); }
        float* pb[1] { probe.data() }; g.applyGain (pb, 1, N);
        bool unity = true; for (float v : probe) unity = unity && (v == 1.0f);
        test::ok (unity, "applyGain leaves the buffer exactly unchanged when off");
        test::approx (g.currentGain(), 1.0, 1.0e-6, "effective gain is exactly unity when off");
    }

    test::group ("loud tone above threshold OPENS");
    {
        NoiseGate g; g.prepare (kSr, N, 2);
        test::ok (runMono (g, sine (220.0, kSr, N, 0.1), 8, N, true, -60.0f) > 0.99f, "-20 dBFS > -60 thr opens fully");
    }

    // A FINITE ±FLT_MAX overflows the sidechain IIR to ±Inf; a naive flush passes Inf and sticks the gate open.
    test::group ("NaN/Inf key does not brick the gate (recovers to closed)");
    {
        NoiseGate g; g.prepare (kSr, N, 2);
        std::vector<float> poison ((std::size_t) N, 0.0f); poison[0] = FLT_MAX; poison[1] = -FLT_MAX;
        { const float* in[1] { poison.data() }; g.analyse (in, 1, N, true, -40.0f); }
        std::vector<float> z ((std::size_t) N, 0.0f); float last = 1.0f;
        for (int b = 0; b < 40; ++b) { const float* in[1] { z.data() }; g.analyse (in, 1, N, true, -40.0f); last = g.currentCoreGain(); }
        test::ok (std::isfinite (last), "coreGain finite after +/-FLT_MAX");
        test::ok (last < 1.0e-3f, "gate RECOVERS to closed (not stuck open)");
    }

    test::group ("GR reading is effective — off reads ~unity, not floor");
    {
        NoiseGate g; g.prepare (kSr, N, 2);
        std::vector<float> z ((std::size_t) N, 0.0f);
        for (int b = 0; b < 4; ++b) { const float* in[1] { z.data() }; g.analyse (in, 1, N, false, -60.0f); }
        test::ok (g.currentGain() > 0.99f, "effective gain ~1.0 while disabled");
        test::ok (g.currentCoreGain() < 1.0e-3f, "raw core gain still tracks closed");
    }

    test::group ("block-split NULL: whole vs [1,7,64,513,15] chunks byte-identical");
    {
        auto src = sine (130.0, kSr, 600, 0.2);
        auto render = [&] (std::vector<int> chunks)
        {
            NoiseGate g; g.prepare (kSr, 600, 2);
            std::vector<float> out = src; int off = 0;
            for (int c : chunks)
            {
                int m = std::min (c, (int) src.size() - off); if (m <= 0) break;
                const float* in[1] { src.data() + off }; float* ob[1] { out.data() + off };
                g.analyse (in, 1, m, true, -70.0f); g.applyGain (ob, 1, m); off += m;
            }
            return out;
        };
        auto whole = render ({ 600 });
        auto split = render ({ 1, 7, 64, 513, 15 });
        float maxd = 0.0f; for (std::size_t i = 0; i < whole.size(); ++i) maxd = std::max (maxd, std::fabs (whole[i] - split[i]));
        test::ok (maxd == 0.0f, "split output is byte-identical to the whole-block render");
    }

    test::group ("sustained low-E does not chatter / false-close");
    {
        NoiseGate g; g.prepare (kSr, N, 2);
        auto e = sine (82.41, kSr, N, 0.06);   // ~ -24 dBFS, well above -60
        float minCore = 1.0f;
        for (int b = 0; b < 90; ++b)
        {
            std::vector<float> blk ((std::size_t) N);
            for (int i = 0; i < N; ++i) blk[(std::size_t) i] = e[(std::size_t) ((b * N + i) % (int) e.size())];
            const float* in[1] { blk.data() }; g.analyse (in, 1, N, true, -60.0f);
            if (b > 4) minCore = std::min (minCore, g.currentCoreGain());
        }
        test::ok (minCore > 0.99f, "gate stays fully open across the whole sustained low note");
    }

    test::group ("linear-fast open preserves the transient (unity within ~2 ms)");
    {
        NoiseGate g; g.prepare (kSr, N, 2);
        std::vector<float> z ((std::size_t) N, 0.0f);
        for (int b = 0; b < 20; ++b) { const float* in[1] { z.data() }; g.analyse (in, 1, N, true, -60.0f); }
        std::vector<float> step ((std::size_t) N, 0.5f), probe ((std::size_t) N, 1.0f);
        const float* in[1] { step.data() }; float* pb[1] { probe.data() };
        g.analyse (in, 1, N, true, -60.0f); g.applyGain (pb, 1, N);
        int idx = -1; for (int i = 0; i < N; ++i) if (probe[(std::size_t) i] > 0.5f) { idx = i; break; }
        test::ok (idx >= 0 && idx < 96, "gain reaches ~unity within ~2 ms of the onset");
    }

    test::group ("linked stereo detector (max of L/R)");
    {
        NoiseGate g; g.prepare (kSr, N, 2);
        auto loud = sine (300.0, kSr, N, 0.2); std::vector<float> silent ((std::size_t) N, 0.0f);
        for (int b = 0; b < 8; ++b) { const float* in[2] { silent.data(), loud.data() }; g.analyse (in, 2, N, true, -60.0f); }
        test::ok (g.currentCoreGain() > 0.99f, "loud RIGHT opens the gate even with a silent LEFT");
    }

    test::group ("sample-rate invariance of the open time");
    {
        auto openMs = [] (double sr)
        {
            const int n = (int) sr;
            NoiseGate g; g.prepare (sr, n, 2);
            std::vector<float> z ((std::size_t) n, 0.0f);
            { const float* in[1] { z.data() }; g.analyse (in, 1, n, true, -60.0f); }
            std::vector<float> step ((std::size_t) n, 0.5f), probe ((std::size_t) n, 1.0f);
            const float* in[1] { step.data() }; float* pb[1] { probe.data() };
            g.analyse (in, 1, n, true, -60.0f); g.applyGain (pb, 1, n);
            int idx = 0; for (int i = 0; i < n; ++i) if (probe[(std::size_t) i] > 0.5f) { idx = i; break; }
            return idx * 1000.0 / sr;
        };
        const double m48 = openMs (48000.0), m96 = openMs (96000.0), m441 = openMs (44100.0);
        test::approx (m96,  m48, 0.15, "96 kHz open time matches 48 kHz (ms)");
        test::approx (m441, m48, 0.15, "44.1 kHz open time matches 48 kHz (ms)");
    }

    test::group ("degenerate numChannels = 0 is a safe no-op");
    {
        NoiseGate g; g.prepare (kSr, N, 2);
        const float* none[1] { nullptr }; g.analyse (none, 0, N, true, -60.0f);
        float* nb[1] { nullptr }; g.applyGain (nb, 0, N);
        test::ok (true, "no crash / no OOB on a zero-lane block");
    }

    test::group ("NaN threshold falls back, does not freeze");
    {
        NoiseGate g; g.prepare (kSr, N, 2);
        test::ok (runMono (g, sine (220.0, kSr, N, 0.1), 8, N, true, std::nanf ("")) > 0.99f,
                  "NaN threshold -> fallback, still opens on a loud tone");
    }

    // RT-safety: analyse()/applyGain()/process() must not allocate (all storage sized in prepare()).
    test::group ("no allocation in analyse/applyGain/process");
    {
        NoiseGate g; g.prepare (kSr, N, 2);
        auto s = sine (220.0, kSr, N, 0.2);
        std::vector<float> buf ((std::size_t) N);
        const long before = g_allocs.load();
        for (int b = 0; b < 50; ++b)
        {
            for (int i = 0; i < N; ++i) buf[(std::size_t) i] = s[(std::size_t) i];
            const float* in[1] { buf.data() }; float* io[1] { buf.data() };
            g.analyse (in, 1, N, true, -40.0f);
            g.applyGain (io, 1, N);
            g.process  (io, 1, N, true, -40.0f);
        }
        test::okNoAlloc (g_allocs.load() == before, "analyse/applyGain/process allocated nothing");
    }

    return test::report();
}

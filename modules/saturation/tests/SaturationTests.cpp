// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// JUCE-free self-tests for the saturation module: peak-normalised WaveShaper curves (odd vs even
// harmonics), and the oversampled Saturator (clean dry/wet, ~linear at zero drive, peak-safe under
// full-scale drive, no-alloc in process(), and the DC blocker neutralising the asymmetric curve's offset).

#include <felitronics_test.h>
#include <felitronics/saturation/WaveShaper.h>
#include <felitronics/saturation/Saturator.h>

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
using Shape = saturation::WaveShaper::Shape;

int main()
{
    std::printf ("felitronics::saturation tests\n");
    const double pi = core::kPi;

    // --- WaveShaper transfer curves ---
    test::group ("WaveShaper peak-normalised curves");
    {
        saturation::WaveShaper ws; ws.setShape (Shape::Tanh); ws.setDrive (2.0f);
        test::approx (ws.processSample ( 1.0f),  1.0, 1e-5, "tanh: peak preserved at +1");
        test::approx (ws.processSample (-1.0f), -1.0, 1e-5, "tanh: peak preserved at -1");
        test::approx (ws.processSample ( 0.0f),  0.0, 1e-9, "tanh: through origin");
        test::ok (ws.processSample (0.5f) > 0.5f, "tanh: small-signal boost (norm slope > 1)");
        test::approx (ws.processSample (0.3f) + ws.processSample (-0.3f), 0.0, 1e-6, "tanh: odd (no even harmonics)");
        test::ok (ws.slopeAtZero() > 1.0f, "tanh: slopeAtZero > 1 for drive > 0");

        saturation::WaveShaper a; a.setShape (Shape::Asym); a.setDrive (2.0f); a.setBias (0.4f);
        test::approx (a.processSample (0.0f), 0.0, 1e-6, "asym: through origin");
        test::ok (std::fabs (a.processSample (0.3f) + a.processSample (-0.3f)) > 1e-3, "asym: NOT odd -> even harmonics");
        test::ok (std::fabs (a.processSample ( 1.0f)) <= 1.0001f, "asym: peak bounded (+)");
        test::ok (std::fabs (a.processSample (-1.0f)) <= 1.0001f, "asym: peak bounded (-)");

        saturation::WaveShaper c; c.setShape (Shape::Cubic); c.setDrive (1.5f);
        test::approx (c.processSample (1.0f), 1.0, 1e-5, "cubic: peak preserved at +1");
        test::ok (c.processSample (5.0f) <= 1.0001f, "cubic: hard-bounded beyond the knee");
    }

    // --- Saturator: clean dry/wet + ~linear at zero drive (os = 1) ---
    test::group ("Saturator dry/wet + linearity (os=1)");
    {
        saturation::Saturator s; test::ok (s.prepare (48000.0, 512, 2, 1), "prepare os=1");
        test::ok (s.latencySamples() == 0, "os=1 -> zero latency");

        saturation::Saturator::Params p; p.driveDb = 18.0f; p.mix = 0.0f; p.outputDb = 0.0f; s.setParams (p);
        float a[6] { 0.1f, -0.2f, 0.3f, -0.4f, 0.5f, -0.6f }, ref[6];
        for (int i = 0; i < 6; ++i) ref[i] = a[i];
        float b6[6] {}; float* io[2] { a, b6 };
        s.process (io, 2, 6);
        double md = 0; for (int i = 0; i < 6; ++i) md = std::max (md, (double) std::fabs (a[i] - ref[i]));
        test::ok (md < 1e-6, "mix=0 -> exact dry passthrough");

        saturation::Saturator s2; s2.prepare (48000.0, 512, 1, 1);
        saturation::Saturator::Params lp; lp.driveDb = 0.0f; lp.mix = 1.0f; lp.autoComp = 0.0f; s2.setParams (lp);
        std::vector<float> x (512), y (512);
        for (int i = 0; i < 512; ++i) { x[i] = 0.5f * (float) std::sin (2.0 * pi * 1000.0 * i / 48000.0); y[i] = x[i]; }
        float* io2[1] { y.data() };
        s2.process (io2, 1, 512);
        double dl = 0; for (int i = 64; i < 512; ++i) dl = std::max (dl, (double) std::fabs (y[i] - x[i]));
        test::ok (dl < 5e-3, "drive 0 -> ~linear passthrough (os=1)");
    }

    // --- Saturator oversampled (os = 4): latency, peak-safety, no-alloc ---
    test::group ("Saturator oversampled (os=4)");
    {
        saturation::Saturator s; s.prepare (48000.0, 512, 2, 4);              // tpp = 32 default
        test::ok (s.latencySamples() == 31, "os=4 round-trip latency == tpp-1 == 31");

        saturation::Saturator::Params p; p.shape = Shape::Tanh; p.driveDb = 18.0f; p.mix = 1.0f; p.autoComp = 0.0f; s.setParams (p);
        std::vector<float> L (512), R (512); float* io[2] { L.data(), R.data() };
        double peak = 0;
        for (int blk = 0; blk < 8; ++blk)
        {
            for (int i = 0; i < 512; ++i) { const int n = blk * 512 + i; const float v = (float) std::sin (2.0 * pi * 1000.0 * n / 48000.0); L[i] = v; R[i] = v; }
            s.process (io, 2, 512);
            if (blk >= 4) for (int i = 0; i < 512; ++i) peak = std::max (peak, (double) std::fabs (L[i]));
        }
        test::ok (peak < 1.10, "peak-safe: full-scale sine stays ~bounded (peak-normalised curve)");

        for (int i = 0; i < 512; ++i) { L[i] = 0.3f; R[i] = -0.3f; }
        const long before = g_allocs.load();
        s.process (io, 2, 512); s.process (io, 2, 512);
        test::okNoAlloc (g_allocs.load() == before, "process() did not allocate (os=4)");
    }

    // --- Asymmetric curve -> the DC blocker removes the offset (with a contrast where it's off) ---
    test::group ("Saturator Asym -> DC-blocked");
    {
        const double f = 187.5;   // exactly 4 periods / 1024-sample block @48k -> a block mean = DC only
        saturation::Saturator::Params p; p.shape = Shape::Asym; p.driveDb = 12.0f; p.bias = 0.4f; p.mix = 1.0f; p.autoComp = 0.0f;

        saturation::Saturator on; on.prepare (48000.0, 1024, 1, 4); on.setParams (p);
        saturation::Saturator off; off.prepare (48000.0, 1024, 1, 4);
        saturation::Saturator::Params p0 = p; p0.dcBlockHz = 0.0f; off.setParams (p0);

        std::vector<float> xo (1024), xf (1024); double meanOn = 0, meanOff = 0;
        for (int blk = 0; blk < 8; ++blk)
        {
            for (int i = 0; i < 1024; ++i) { const int n = blk * 1024 + i; const float v = 0.7f * (float) std::sin (2.0 * pi * f * n / 48000.0); xo[i] = v; xf[i] = v; }
            float* ioOn[1] { xo.data() }; float* ioOff[1] { xf.data() };
            on.process (ioOn, 1, 1024); off.process (ioOff, 1, 1024);
            if (blk == 7)
            {
                double sOn = 0, sOff = 0; for (int i = 0; i < 1024; ++i) { sOn += xo[i]; sOff += xf[i]; }
                meanOn = sOn / 1024.0; meanOff = sOff / 1024.0;
            }
        }
        test::ok (std::fabs (meanOn)  < 1e-3, "DC blocker ON  -> output block-mean ~ 0");
        test::ok (std::fabs (meanOff) > 1e-3, "DC blocker OFF -> asym curve leaves a measurable DC offset");
    }

    return test::report();
}

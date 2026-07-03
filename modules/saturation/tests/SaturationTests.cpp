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

    // --- lifecycle/misuse: process() after a FAILED prepare (partial init) or with n>maxBlock must not OOB ---
    // prepare() sets channels_ (clamped >=1) BEFORE ovs_.prepare() can fail — a failed prepare left channels_>=1
    // with empty osBuf_/osPtrs_, and process() (nc>=1) indexed them; an oversized block overran osBuf_. Guarded now.
    test::group ("Saturator: reject process before / after failed prepare, and oversized blocks");
    {
        saturation::Saturator sat;                                   // NOT prepared (channels_ == 0)
        float a[32] {}, b[32] {}; float* io[2] { a, b };
        sat.process (io, 2, 16);                                     // channels_==0 → safe no-op
        test::ok (! sat.prepare (48000.0, 16, 2, 4, 2), "prepare(tapsPerPhase=2) fails (oversampler rejects <4)");
        sat.process (io, 2, 16);                                     // FAILED prepare → no-op, must not index empty osBuf_
        test::ok (sat.prepare (48000.0, 16, 2, 4, 32), "prepare valid");
        sat.process (io, 2, 16);                                     // works
        sat.process (io, 2, 32);                                     // n=32 > maxBlock=16 → no-op, must not overrun osBuf_
        test::ok (true, "no OOB across failed-prepare / oversized-block process (ASan/UBSan is the real check)");
    }

    // --- FALSIFICATION: at os>1 the dry/wet mix must be latency-aligned (undelayed dry combs the wet) ---
    test::group ("Saturator os=4 dry/wet mix is latency-aligned (no comb)");
    {
        const int n = 4096;
        std::vector<float> x (n);
        for (int i = 0; i < n; ++i) x[i] = (float) (0.5 * std::sin (2.0 * pi * 1000.0 * i / 48000.0));
        std::vector<float> y = x; float* ch[1] { y.data() };
        saturation::Saturator sat; sat.prepare (48000.0, n, 1, 4);
        saturation::Saturator::Params p;
        p.driveDb = 0.0f; p.autoComp = 0.0f; p.mix = 0.5f; p.outputDb = 0.0f;   // ≈linear curve → wet ≈ delayed dry
        sat.setParams (p);
        sat.process (ch, 1, n);
        auto rmsHalf = [] (const std::vector<float>& v) {
            double s2 = 0.0; const int from = (int) v.size() / 2;
            for (int i = from; i < (int) v.size(); ++i) s2 += (double) v[i] * v[i];
            return std::sqrt (s2 / ((int) v.size() - from)); };
        // aligned: 0.5·x + 0.5·x == x; unaligned by 31 samples @1 kHz: |0.5 + 0.5·e^{-jθ}| ≈ 0.44
        test::approx (rmsHalf (y) / rmsHalf (x), 1.0, 0.05, "mix=0.5 at ~zero drive is level-neutral (dry delayed to match wet)");
    }

    // --- autoComp=1 → unity small-signal gain regardless of drive ---
    test::group ("Saturator autoComp=1 small-signal unity");
    {
        const int n = 4096;
        std::vector<float> x (n);
        for (int i = 0; i < n; ++i) x[i] = (float) (0.005 * std::sin (2.0 * pi * 500.0 * i / 48000.0));
        std::vector<float> y = x; float* ch[1] { y.data() };
        saturation::Saturator sat; sat.prepare (48000.0, n, 1, 4);
        saturation::Saturator::Params p; p.driveDb = 12.0f; p.autoComp = 1.0f; p.mix = 1.0f;
        sat.setParams (p);
        sat.process (ch, 1, n);
        double sx = 0.0, sy = 0.0;
        for (int i = n / 2; i < n; ++i) { sx += (double) x[i] * x[i]; sy += (double) y[i] * y[i]; }
        test::approx (std::sqrt (sy / sx), 1.0, 0.02, "12 dB drive + autoComp 1 → tiny signal passes at unity");
    }

    // --- non-finite params must not poison the stream (house rule: clamp + reject non-finite) ---
    test::group ("Saturator non-finite params rejected");
    {
        const int n = 512;
        std::vector<float> y (n);
        for (int i = 0; i < n; ++i) y[i] = (float) (0.5 * std::sin (2.0 * pi * 997.0 * i / 48000.0));
        float* ch[1] { y.data() };
        saturation::Saturator sat; sat.prepare (48000.0, n, 1, 4);
        saturation::Saturator::Params p;
        p.shape = Shape::Asym;
        const float qnan = std::numeric_limits<float>::quiet_NaN();
        p.bias = qnan; p.mix = qnan; p.outputDb = qnan; p.autoComp = qnan; p.dcBlockHz = qnan; p.driveDb = qnan;
        sat.setParams (p);
        sat.process (ch, 1, n);
        bool finite = true; for (float v : y) finite &= (bool) std::isfinite (v);
        test::ok (finite, "all-NaN params → finite output (fallbacks applied)");
    }

    return test::report();
}

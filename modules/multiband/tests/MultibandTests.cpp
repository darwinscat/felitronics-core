// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// JUCE-free self-tests for the multiband engine: with no compression the bands recombine to the splitter's
// flat (allpass) reconstruction; a per-band compressor lowers ONLY its band's level + reports gain
// reduction; solo isolates a band; latency follows the max per-band lookahead; process() never allocates.

#include <felitronics_test.h>
#include <felitronics/core/Math.h>
#include <felitronics/multiband/MultibandProcessor.h>
#include <felitronics/multiband/MultibandCompressor.h>

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

static double rms (const std::vector<float>& v, int from)
{
    double e = 0.0; int c = 0;
    for (int i = from; i < (int) v.size(); ++i) { e += (double) v[i] * v[i]; ++c; }
    return c ? std::sqrt (e / (double) c) : 0.0;
}

static dynamics::CompressorParams comp (double thr, double ratio, double lookMs = 0.0)
{
    dynamics::CompressorParams p;
    p.thresholdDb = thr; p.ratio = ratio; p.attackMs = 5.0; p.releaseMs = 50.0; p.lookaheadMs = lookMs;
    return p;
}

static void runBlocks (multiband::MultibandCompressor<4>& mc, std::vector<float>& y, int block = 512)
{
    for (int o = 0; o < (int) y.size(); o += block)
    {
        float* io[1] { y.data() + o };
        mc.process (io, 1, std::min (block, (int) y.size() - o));
    }
}

int main()
{
    std::printf ("felitronics::multiband tests\n");
    const double sr = 48000.0, pi = core::kPi;

    // --- no compression → bands recombine to the flat (allpass) reconstruction ---
    test::group ("MultibandCompressor flat when not compressing");
    {
        multiband::MultibandCompressor<4> mc; mc.prepare (sr, 512, 1); mc.setNumBands (3);
        const float xf[2] = { 250.0f, 2500.0f }; mc.setCrossovers (xf, 2);
        for (int b = 0; b < 3; ++b) mc.setBandParams (b, comp (12.0, 2.0));     // threshold above signal → no GR
        const int N = 8000; std::vector<float> x (N), y (N);
        unsigned long long s = 7; auto rng = [&]() { s = s * 6364136223846793005ULL + 1442695040888963407ULL; return (float) ((s >> 40) & 0xffff) / 32768.0f - 1.0f; };
        for (int i = 0; i < N; ++i) { x[i] = 0.3f * rng(); y[i] = x[i]; }
        runBlocks (mc, y);
        test::ok (std::fabs (rms (y, N / 2) / rms (x, N / 2) - 1.0) < 0.05, "no-compression output == allpass reconstruction (RMS flat)");
    }

    // --- a per-band compressor lowers ONLY its band's level + reports GR ---
    test::group ("MultibandCompressor low band compresses");
    {
        multiband::MultibandCompressor<4> mc; mc.prepare (sr, 512, 1); mc.setNumBands (3);
        const float xf[2] = { 250.0f, 2500.0f }; mc.setCrossovers (xf, 2);
        auto run = [&] (double thr) -> double
        {
            mc.reset();
            mc.setBandParams (0, comp (thr, 4.0)); mc.setBandParams (1, comp (12.0, 2.0)); mc.setBandParams (2, comp (12.0, 2.0));
            const int N = 12000; std::vector<float> y (N);
            for (int i = 0; i < N; ++i) y[i] = 0.7f * (float) std::sin (2.0 * pi * 80.0 * i / sr);   // 80 Hz → band 0
            runBlocks (mc, y);
            return rms (y, N / 2);
        };
        const double on  = run (-30.0);
        const double gr  = mc.bandGainReductionDb (0);
        const double off = run (12.0);
        test::ok (on < 0.6 * off, "compressing the low band lowers its level");
        test::ok (gr < -3.0, "low band reports gain reduction");
    }

    // --- solo isolates a band ---
    test::group ("MultibandCompressor solo isolates a band");
    {
        multiband::MultibandCompressor<4> mc; mc.prepare (sr, 512, 1); mc.setNumBands (3);
        const float xf[2] = { 250.0f, 2500.0f }; mc.setCrossovers (xf, 2);
        for (int b = 0; b < 3; ++b) mc.setBandParams (b, comp (12.0, 2.0));
        mc.setBandSolo (0, true);                                              // only the low band
        const int N = 8000; std::vector<float> y (N);
        for (int i = 0; i < N; ++i) y[i] = 0.5f * (float) std::sin (2.0 * pi * 5000.0 * i / sr);   // 5 kHz → band 2
        runBlocks (mc, y);
        test::ok (rms (y, N / 2) < 0.1 * (0.5 / std::sqrt (2.0)), "solo band 0 mutes the 5 kHz (band 2) content");
    }

    // --- latency follows the max per-band lookahead ---
    test::group ("MultibandCompressor latency");
    {
        multiband::MultibandCompressor<4> mc; mc.prepare (sr, 512, 2, 5.0); mc.setNumBands (3);
        for (int b = 0; b < 3; ++b) mc.setBandParams (b, comp (12.0, 2.0));
        test::ok (mc.latencySamples() == 0, "no lookahead → 0 latency");
        mc.setBandParams (1, comp (12.0, 2.0, 2.0));                          // band 1: 2 ms lookahead
        test::ok (mc.latencySamples() == (int) std::lround (2.0 * 0.001 * sr), "latency == the max band lookahead");
    }

    // --- no allocation in process() ---
    test::group ("MultibandCompressor no-alloc");
    {
        multiband::MultibandCompressor<4> mc; mc.prepare (sr, 512, 2); mc.setNumBands (4);
        for (int b = 0; b < 4; ++b) mc.setBandParams (b, comp (-10.0, 4.0, 1.0));
        std::vector<float> l (512, 0.2f), r (512, -0.2f); float* io[2] { l.data(), r.data() };
        mc.process (io, 2, 512);
        const long before = g_allocs.load();
        mc.process (io, 2, 512); mc.process (io, 2, 512);
        test::okNoAlloc (g_allocs.load() == before, "process() did not allocate (4-band stereo, lookahead)");
    }

    // --- solo overrides the parallel dry/wet (the bug it guards: at mix<1 the dry leaked non-solo bands) ---
    test::group ("MultibandCompressor solo overrides dry/wet mix");
    {
        multiband::MultibandCompressor<4> mc; mc.prepare (sr, 512, 1); mc.setNumBands (3);
        const float xf[2] = { 250.0f, 2500.0f }; mc.setCrossovers (xf, 2);
        for (int b = 0; b < 3; ++b) mc.setBandParams (b, comp (12.0, 2.0));
        mc.setBandSolo (0, true); mc.setMix (0.5f);                          // solo low + 50% parallel mix
        const int N = 8000; std::vector<float> y (N);
        for (int i = 0; i < N; ++i) y[i] = 0.5f * (float) std::sin (2.0 * pi * 5000.0 * i / sr);   // 5 kHz → band 2
        runBlocks (mc, y);
        test::ok (rms (y, N / 2) < 0.1 * (0.5 / std::sqrt (2.0)), "solo+mix=0.5: non-solo band does NOT leak through the dry");
    }

    // --- unity bands ⇒ output nulls SAMPLE-WISE against the splitter's own allpass reconstruction ---
    test::group ("MultibandProcessor perfect reconstruction (unity bands)");
    {
        multiband::MultibandCompressor<4> mc; mc.prepare (sr, 512, 1); mc.setNumBands (3);
        const float xf[2] = { 250.0f, 2500.0f }; mc.setCrossovers (xf, 2);
        for (int b = 0; b < 3; ++b) mc.setBandParams (b, comp (24.0, 2.0));  // never compresses → unity
        eq::MultibandSplitter<4> ref; ref.prepare (sr, 1); ref.setNumBands (3); ref.setCrossovers (xf, 2);
        unsigned long long s = 3; auto rng = [&]() { s = s * 6364136223846793005ULL + 1442695040888963407ULL; return (float) ((s >> 40) & 0xffff) / 32768.0f - 1.0f; };
        const int N = 6000; std::vector<float> x (N), y (N);
        for (int i = 0; i < N; ++i) { x[i] = 0.3f * rng(); y[i] = x[i]; }
        runBlocks (mc, y);
        float band[4]; double md = 0;
        for (int i = 0; i < N; ++i) { ref.splitSample (0, x[i], band); const float r = ref.sumSample (band); if (i >= 1000) md = std::max (md, (double) std::fabs (y[i] - r)); }
        test::ok (md < 1e-3, "unity bands → output == the allpass reconstruction (sample-wise)");
    }

    // --- a lookahead band stays TIME-ALIGNED: output == reconstruction delayed by the latency (no comb) ---
    test::group ("MultibandProcessor latency alignment");
    {
        const double lookMs = 1.0; const int L = (int) std::lround (lookMs * 0.001 * sr);
        multiband::MultibandCompressor<4> mc; mc.prepare (sr, 512, 1, 5.0); mc.setNumBands (3);
        const float xf[2] = { 250.0f, 2500.0f }; mc.setCrossovers (xf, 2);
        mc.setBandParams (0, comp (24.0, 2.0)); mc.setBandParams (1, comp (24.0, 2.0, lookMs)); mc.setBandParams (2, comp (24.0, 2.0));
        eq::MultibandSplitter<4> ref; ref.prepare (sr, 1); ref.setNumBands (3); ref.setCrossovers (xf, 2);
        core::DelayLine refDelay; refDelay.prepare (L + 1); refDelay.setDelay (L);
        unsigned long long s = 5; auto rng = [&]() { s = s * 6364136223846793005ULL + 1442695040888963407ULL; return (float) ((s >> 40) & 0xffff) / 32768.0f - 1.0f; };
        const int N = 6000; std::vector<float> x (N), y (N);
        for (int i = 0; i < N; ++i) { x[i] = 0.3f * rng(); y[i] = x[i]; }
        runBlocks (mc, y);
        float band[4]; double md = 0;
        for (int i = 0; i < N; ++i) { ref.splitSample (0, x[i], band); const float r = refDelay.process (ref.sumSample (band)); if (i >= L + 1000) md = std::max (md, (double) std::fabs (y[i] - r)); }
        test::ok (mc.latencySamples() == L && md < 1e-3, "all bands aligned to the lookahead latency (no comb in the sum)");
    }

    // --- parallel mix law: dry (mix=0) louder than the compressed wet (mix=1), 0.5 between ---
    test::group ("MultibandCompressor parallel mix law");
    {
        auto runMix = [&] (float mix) -> double
        {
            multiband::MultibandCompressor<4> mc; mc.prepare (sr, 512, 1); mc.setNumBands (3);
            const float xf[2] = { 250.0f, 2500.0f }; mc.setCrossovers (xf, 2);
            mc.setBandParams (0, comp (-30.0, 4.0)); mc.setBandParams (1, comp (24.0, 2.0)); mc.setBandParams (2, comp (24.0, 2.0));
            mc.setMix (mix);
            const int N = 12000; std::vector<float> y (N);
            for (int i = 0; i < N; ++i) y[i] = 0.7f * (float) std::sin (2.0 * pi * 80.0 * i / sr);
            runBlocks (mc, y);
            return rms (y, N / 2);
        };
        const double m0 = runMix (0.0f), m05 = runMix (0.5f), m1 = runMix (1.0f);
        test::ok (m0 > m1 * 1.3, "mix=0 (dry) louder than mix=1 (compressed wet)");
        test::ok (m05 > m1 && m05 < m0, "mix=0.5 sits between dry and wet");
    }

    // --- stereo channel independence ---
    test::group ("MultibandCompressor stereo independence");
    {
        multiband::MultibandCompressor<4> mc; mc.prepare (sr, 512, 2); mc.setNumBands (3);
        const float xf[2] = { 250.0f, 2500.0f }; mc.setCrossovers (xf, 2);
        for (int b = 0; b < 3; ++b) mc.setBandParams (b, comp (12.0, 2.0));
        const int N = 6000; std::vector<float> l (N), r (N, 0.0f);
        unsigned long long s = 9; auto rng = [&]() { s = s * 6364136223846793005ULL + 1442695040888963407ULL; return (float) ((s >> 40) & 0xffff) / 32768.0f - 1.0f; };
        for (int i = 0; i < N; ++i) l[i] = 0.3f * rng();
        for (int o = 0; o < N; o += 512) { float* io[2] { l.data() + o, r.data() + o }; mc.process (io, 2, std::min (512, N - o)); }
        test::ok (rms (r, N / 2) == 0.0, "right channel stays silent while left is driven (independent state)");
    }

    return test::report();
}

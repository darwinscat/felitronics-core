// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// JUCE-free self-tests for the Dither. The properties that matter for a mastering dither:
//   (1) TPDF dither LINEARISES — a sub-LSB DC survives as the time-average of the quantized output (without
//       dither it would round to a flat zero); (2) the output really is on the bit-depth grid; (3) the dither
//       level is the textbook TPDF ½LSB+½LSB RMS = lsb/√6; (4) noise shaping is HIGH-PASS (NTF = 1−H, the
//       error's lag-1 autocorrelation is NEGATIVE) and so REDUCES in-band (LF) noise — this also pins the
//       feedback SIGN; (5) auto-blanking yields true digital black on a silent tail; (6) L/R dither is
//       decorrelated (identical input → different output); (7) bypass at 32-bit; (8) no alloc; finite guard.

#include <felitronics_test.h>
#include <felitronics/core/Math.h>
#include <felitronics/dither/Dither.h>

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
using D  = dither::Dither;
using DP = dither::DitherParams;
using NS = dither::NoiseShaping;

// run a constant-input block through a fresh Dither, return the output
static std::vector<float> runConst (const DP& p, double dc, int N)
{
    D d; d.prepare (48000.0, 1024, 1); d.setParams (p);
    std::vector<float> y (N, (float) dc);
    for (int o = 0; o < N; o += 1024) { float* io[1] { y.data() + o }; d.process (io, 1, std::min (1024, N - o)); }
    return y;
}

// the error signal (out − in) for a 1 kHz tone at `amp`, fresh Dither
static std::vector<double> toneError (const DP& p, double amp, int N, double sr)
{
    D d; d.prepare (sr, 1024, 1); d.setParams (p);
    std::vector<float> y (N); std::vector<double> in (N);
    for (int i = 0; i < N; ++i) { in[i] = amp * std::sin (2.0 * core::kPi * 1000.0 * i / sr); y[i] = (float) in[i]; }
    for (int o = 0; o < N; o += 1024) { float* io[1] { y.data() + o }; d.process (io, 1, std::min (1024, N - o)); }
    std::vector<double> e (N); for (int i = 0; i < N; ++i) e[i] = (double) y[i] - in[i];
    return e;
}

int main()
{
    std::printf ("felitronics::dither tests\n");
    const double sr = 48000.0;
    const double lsb16 = 1.0 / 32768.0;                      // 16-bit LSB

    // --- (1) TPDF dither linearises: a sub-LSB DC survives as the mean of the quantized output ---
    test::group ("Dither linearises sub-LSB DC");
    {
        DP p; p.bits = 16; p.shaping = NS::None;
        const double dc = 0.3 * lsb16;                        // well below one LSB → would round to 0 undithered
        const auto y = runConst (p, dc, 200000);
        double mean = 0; for (float v : y) mean += v; mean /= (double) y.size();
        test::ok (std::fabs (mean - dc) < 0.05 * lsb16, "mean(dithered) tracks a 0.3-LSB DC (un-dithered would be flat 0)");
    }

    // --- (2) output lands on the quantization grid, within the signed-PCM code range ---
    test::group ("Dither output is on the bit grid");
    {
        DP p; p.bits = 16; p.shaping = NS::Weighted;
        const auto y = runConst (p, 0.123456, 8000);
        double worst = 0; bool inRange = true;
        for (float v : y) { const double code = (double) v * 32768.0; worst = std::max (worst, std::fabs (code - std::round (code))); if (v < -1.0f || v >= 1.0f) inRange = false; }
        test::ok (worst < 1e-3 && inRange, "every output sample is an integer/32768 within [-1, 1)");
    }

    // --- (3) the total added noise is textbook TPDF: RMS = lsb/2 (= +4.77 dB over undithered lsb/√12) ---
    test::group ("Dither level is TPDF (lsb/2)");
    {
        const auto e = toneError ({ .bits = 16, .shaping = NS::None }, 0.2, 100000, sr);
        double sq = 0; for (double v : e) sq += v * v; const double rms = std::sqrt (sq / (double) e.size());
        const double expect = lsb16 * 0.5;                    // var = lsb²/4 for 2-LSB-pp TPDF (dither + quant residue)
        test::approx (rms, expect, 0.05 * expect, "total error RMS ≈ lsb/2 (the well-known +4.77 dB TPDF penalty)");
    }

    // --- (4) noise shaping is HIGH-PASS (NTF = 1−H): negative lag-1 autocorrelation, reduced LF noise ---
    test::group ("Dither shaping is high-pass (NTF sign)");
    {
        auto ac1 = [] (const std::vector<double>& e) { double n = 0, d = 0; for (std::size_t i = 1; i < e.size(); ++i) { n += e[i] * e[i - 1]; d += e[i] * e[i]; } return d > 0 ? n / d : 0.0; };
        // Welch PSD estimate at frequency f: averaged single-bin Goertzel over K segments (kills the χ²-variance)
        auto welchPow = [&] (const std::vector<double>& e, double f, int seg) { double tot = 0; int K = 0; for (std::size_t base = 0; base + (std::size_t) seg <= e.size(); base += (std::size_t) seg) { double re = 0, im = 0; for (int j = 0; j < seg; ++j) { const double a = 2.0 * core::kPi * f * j / sr; re += e[base + j] * std::cos (a); im += e[base + j] * std::sin (a); } tot += re * re + im * im; ++K; } return K > 0 ? tot / K : 0.0; };
        const auto eNone = toneError ({ .bits = 16, .shaping = NS::None },     0.2, 100000, sr);
        const auto eWtd  = toneError ({ .bits = 16, .shaping = NS::Weighted }, 0.2, 100000, sr);
        const auto ePsy  = toneError ({ .bits = 16, .shaping = NS::Psychoacoustic }, 0.2, 100000, sr);
        test::ok (std::fabs (ac1 (eNone)) < 0.05,           "None: white error (lag-1 autocorr ≈ 0)");
        test::ok (ac1 (eWtd) < -0.2,                        "Weighted: HF-emphasised error (lag-1 autocorr < 0 → SUBTRACT sign is right)");
        test::ok (ac1 (ePsy) < 0.0,                         "Psychoacoustic: net high-pass (lag-1 autocorr < 0)");
        test::ok (welchPow (eWtd, 200.0, 1000) < 0.1 * welchPow (eNone, 200.0, 1000), "Weighted slashes 200 Hz (in-band) noise vs flat TPDF — the audible win");
    }

    // --- (4b) the 9th-order F-weighted curve is STABLE under hot input (bounded feedback) ---
    test::group ("Dither psychoacoustic is stable");
    {
        DP p; p.bits = 16; p.shaping = NS::Psychoacoustic;
        D d; d.prepare (sr, 1024, 1); d.setParams (p);
        unsigned long s = 1; auto rng = [&]() { s = s * 6364136223846793005ULL + 1442695040888963407ULL; return (float) ((s >> 40) & 0xffff) / 32768.0f - 1.0f; };
        const int N = 200000; std::vector<float> y (N); for (int i = 0; i < N; ++i) y[i] = 0.95f * rng();
        double mx = 0; for (int o = 0; o < N; o += 1024) { float* io[1] { y.data() + o }; d.process (io, 1, std::min (1024, N - o)); }
        for (float v : y) mx = std::max (mx, (double) std::fabs (v));
        test::ok (mx < 1.001, "9th-order feedback stays bounded on full-scale noise (no runaway)");
    }

    // --- (5) auto-blanking → true digital black on a silent tail; off → dither keeps running ---
    test::group ("Dither auto-blanking");
    {
        const int N = 8000;
        DP on;  on.bits = 16;  on.shaping = NS::Weighted; on.autoBlank = true;  on.autoBlankSamples = 1000;
        DP off = on; off.autoBlank = false;
        auto tailEnergy = [&] (const DP& p) { D d; d.prepare (sr, 1024, 1); d.setParams (p); std::vector<float> y (N, 0.0f); for (int o = 0; o < N; o += 1024) { float* io[1] { y.data() + o }; d.process (io, 1, std::min (1024, N - o)); } double s = 0; for (int i = N - 2000; i < N; ++i) s += (double) y[i] * y[i]; return s; };
        test::ok (tailEnergy (on)  == 0.0, "autoBlank on: a sustained-silence tail is exactly zero (digital black)");
        test::ok (tailEnergy (off) >  0.0, "autoBlank off: silence still gets dither (a steady noise floor)");
    }

    // --- (6) stereo dither is DECORRELATED: identical input at a code boundary → different output per channel ---
    test::group ("Dither decorrelates L/R");
    {
        DP p; p.bits = 16; p.shaping = NS::None;
        D d; d.prepare (sr, 1024, 2); d.setParams (p);
        const int N = 4000; std::vector<float> l (N, (float) (0.5 * lsb16)), r (N, (float) (0.5 * lsb16));   // sit on the 0↔1-LSB boundary
        for (int o = 0; o < N; o += 1024) { float* io[2] { l.data() + o, r.data() + o }; d.process (io, 2, std::min (1024, N - o)); }
        int diff = 0; for (int i = 0; i < N; ++i) if (l[i] != r[i]) ++diff;
        test::ok (diff > N / 4, "identical boundary input → independent per-channel dither flips rounding differently");
    }

    // --- (7) bits ≥ 32 → bit-exact bypass (float export) ---
    test::group ("Dither 32-bit bypass");
    {
        DP p; p.bits = 32;
        D d; d.prepare (sr, 1024, 1); d.setParams (p);
        const int N = 1000; std::vector<float> y (N), y0 (N);
        unsigned long s = 7; auto rng = [&]() { s = s * 6364136223846793005ULL + 1442695040888963407ULL; return (float) ((s >> 40) & 0xffff) / 32768.0f - 1.0f; };
        for (int i = 0; i < N; ++i) { y[i] = 0.5f * rng(); y0[i] = y[i]; }
        for (int o = 0; o < N; o += 1024) { float* io[1] { y.data() + o }; d.process (io, 1, std::min (1024, N - o)); }
        double md = 0; for (int i = 0; i < N; ++i) md = std::max (md, (double) std::fabs (y[i] - y0[i]));
        test::ok (md == 0.0, "bits=32 → exact passthrough (no quantization, no dither)");
    }

    // --- (8) no alloc in process(); a non-finite input cannot poison the output ---
    test::group ("Dither no-alloc + finite guard");
    {
        DP p; p.bits = 24; p.shaping = NS::Psychoacoustic;
        D d; d.prepare (sr, 1024, 2); d.setParams (p);
        std::vector<float> l (512, 0.3f), r (512, -0.2f); float* io[2] { l.data(), r.data() };
        d.process (io, 2, 512);
        const long before = g_allocs.load();
        d.process (io, 2, 512); d.process (io, 2, 512);
        const bool noAlloc = (g_allocs.load() == before);
        l[10] = std::nanf (""); r[20] = INFINITY;
        d.process (io, 2, 512);
        bool fin = true; for (int i = 0; i < 512; ++i) if (! std::isfinite (l[i]) || ! std::isfinite (r[i])) fin = false;
        test::ok (noAlloc, "process() did not allocate");
        test::ok (fin,     "NaN/inf input → finite output (guarded)");
    }

    // --- (9) NO NOISE MODULATION: error variance is independent of signal level — the property the
    //     "shape the dither too" topology could have broken (the top adversarial concern) ---
    test::group ("Dither has no noise modulation");
    {
        auto errVar = [&] (NS shaping, double dc) {
            DP p; p.bits = 16; p.shaping = shaping; p.autoBlank = false;   // off so DC=0 doesn't blank
            const auto y = runConst (p, dc, 120000);
            double m = 0; for (float v : y) m += (double) v - dc; m /= (double) y.size();
            double s = 0; for (float v : y) { const double e = (double) v - dc - m; s += e * e; } return s / (double) y.size();
        };
        for (NS sh : { NS::None, NS::Weighted })
        {
            double lo = 1e30, hi = 0; for (double dc : { -0.4, -0.1, 0.0, 0.1, 0.4 }) { const double v = errVar (sh, dc); lo = std::min (lo, v); hi = std::max (hi, v); }
            test::ok (hi < 1.05 * lo, sh == NS::None ? "None: error variance flat across a DC sweep (TPDF, not RPDF)"
                                                     : "Weighted: variance still flat with the dither shaped (independence survives)");
        }
    }

    // --- (10) higher bit depths quantise on-grid and the dither linearises there too ---
    test::group ("Dither at 20- and 24-bit");
    {
        for (int bits : { 20, 24 })
        {
            const double lsb = 1.0 / (double) ((std::int64_t) 1 << (bits - 1));
            DP p; p.bits = bits; p.shaping = NS::None;
            const auto y = runConst (p, 0.31 * lsb, 300000);
            double mean = 0, worst = 0; for (float v : y) { mean += v; const double code = (double) v / lsb; worst = std::max (worst, std::fabs (code - std::round (code))); }
            mean /= (double) y.size();
            test::ok (worst < 1e-2 && std::fabs (mean - 0.31 * lsb) < 0.05 * lsb, bits == 20 ? "20-bit: on grid + sub-LSB DC linearised" : "24-bit: on grid + sub-LSB DC linearised");
        }
    }

    // --- (11) total noise power INCREASES under shaping (LF is traded for more HF — the documented cost) ---
    test::group ("Dither shaping raises total noise power");
    {
        const auto eNone = toneError ({ .bits = 16, .shaping = NS::None },     0.2, 60000, sr);
        const auto eWtd  = toneError ({ .bits = 16, .shaping = NS::Weighted }, 0.2, 60000, sr);
        auto pow = [] (const std::vector<double>& e) { double s = 0; for (double v : e) s += v * v; return s; };
        test::ok (pow (eWtd) > 1.2 * pow (eNone), "Weighted total power > flat TPDF (the LF cut is paid for in HF)");
    }

    // --- (12) reset is deterministic; a different seed gives a different stream ---
    test::group ("Dither reset determinism + seed");
    {
        DP p; p.bits = 16; p.shaping = NS::Weighted; p.seed = 12345;
        auto run = [&] (const DP& q) { D d; d.prepare (sr, 1024, 1); d.setParams (q); std::vector<float> y (4000); unsigned long s = 3; auto rng = [&]() { s = s * 6364136223846793005ULL + 1442695040888963407ULL; return (float) ((s >> 40) & 0xffff) / 32768.0f - 1.0f; }; for (int i = 0; i < 4000; ++i) y[i] = 0.2f * rng(); for (int o = 0; o < 4000; o += 1024) { float* io[1] { y.data() + o }; d.process (io, 1, std::min (1024, 4000 - o)); } return y; };
        const auto a = run (p), b = run (p);
        double md = 0; for (int i = 0; i < 4000; ++i) md = std::max (md, (double) std::fabs (a[i] - b[i]));
        test::ok (md == 0.0, "same seed + fresh prepare/reset → identical output (deterministic dither)");
        DP p2 = p; p2.seed = 999;
        const auto c = run (p2);
        int diff = 0; for (int i = 0; i < 4000; ++i) if (a[i] != c[i]) ++diff;
        test::ok (diff > 1000, "a different seed → a different dither stream");
    }

    return test::report();
}

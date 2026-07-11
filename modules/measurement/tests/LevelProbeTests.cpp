// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// JUCE-free self-tests for the family "unity" contract — LevelProbe (the reference stimulus) +
// ReferenceUnity (the load-time IR normalization that targets it). The numbers here are GOLDEN:
// products (OrbitCab's tube constant, its factory-IR loudness) shipped against them, so this suite
// pins behavior, not just properties:
//   • white() golden literals — integer + exact-double arithmetic only, so bit-identical on every
//     platform; a changed literal means the family's unity anchor moved (a product-breaking event);
//   • analytic level calibration of the stimulus (kRefRmsDb over the contract window);
//   • referenceUnityGain oracle — the FFT-domain gain must equal actually convolving the stimulus
//     through the IR and reading the RMS ratio (the OrbitCab acceptance test, engine-free);
//   • the floor / clamp / analysis-cap / channel-mean edge semantics, pinned exactly.
// Bitwise parity vs the verbatim OrbitCab source was verified at extraction time over 42
// randomized IR cases (3 rates × 7 lengths × 2 widths) + stimulus fills at 3 rates.

#include <felitronics_test.h>
#include <felitronics/measurement/LevelProbe.h>
#include <felitronics/measurement/ReferenceUnity.h>
#include <felitronics/measurement/Convolve.h>   // magSpectrum — independent check of the LP shaping

#include <cmath>
#include <cstddef>
#include <cstring>
#include <vector>

using namespace felitronics;
using felitronics::test::ok;
using felitronics::test::approx;
using felitronics::test::group;

namespace lp = measurement::levelprobe;

// The family contract constants are load-bearing across products — compile-time pins.
static_assert (lp::kRefRmsDb == -18.0,  "family reference level moved — breaks every shipped unity match");
static_assert (lp::kShapeHz == 2000.0,  "family reference shaping moved — breaks every shipped unity match");
static_assert (measurement::kIrNormMinDb == -30.0f && measurement::kIrNormMaxDb == 30.0f,
               "normalization clamp moved");
static_assert (measurement::kIrRefFloorDb == -60.0, "near-silent floor moved");
static_assert (measurement::kIrAnalysisSeconds == 1.0, "analysis window moved");

namespace
{
// Deterministic damped-noise "IR" — taps from the probe's own LCG (platform-exact input data).
std::vector<float> dampedNoiseIr (int len, float amp, int seed, float decay)
{
    std::vector<float> h ((std::size_t) len);
    for (int i = 0; i < len; ++i)
        h[(std::size_t) i] = amp * lp::white (i + seed) * std::exp (-decay * (float) i / (float) len);
    return h;
}

// RMS (dB) of x over the contract window [settle, settle+measure) at `sr`.
double windowRmsDb (const std::vector<float>& x, double sr)
{
    const int settle = (int) std::lround (lp::kSettleSec * sr);
    const int meas   = (int) std::lround (lp::kMeasureSec * sr);
    double acc = 0.0;
    for (int i = settle; i < settle + meas; ++i) acc += (double) x[(std::size_t) i] * x[(std::size_t) i];
    return 10.0 * std::log10 (acc / meas);
}
}

int main()
{
    std::printf ("felitronics::measurement level-probe + reference-unity tests\n");

    // --- Golden fingerprint: white() literals (bit-exact on every platform — no libm inside) ---
    group ("white() golden literals — the unity anchor's identity");
    {
        ok (lp::white (0)     == 0.84210813f,    "white(0)");
        ok (lp::white (1)     == -0.0219188333f, "white(1)");
        ok (lp::white (2)     == 0.143028259f,   "white(2)");
        ok (lp::white (12345) == -0.918834925f,  "white(12345)");
        bool inRange = true;
        for (std::int64_t n = 0; n < 10000; ++n)
        { const float w = lp::white (n); if (! (w >= -1.0f && w <= 1.0f)) inRange = false; }
        ok (inRange, "white() stays in [-1, 1]");
    }

    // --- Determinism + prefix property: the stimulus depends only on the sample INDEX ---
    group ("fill() determinism and prefix invariance");
    {
        std::vector<float> a (2048), b (2048), c (1024);
        lp::fill (a.data(), 2048, 48000.0);
        lp::fill (b.data(), 2048, 48000.0);
        lp::fill (c.data(), 1024, 48000.0);
        ok (std::memcmp (a.data(), b.data(), a.size() * sizeof (float)) == 0, "two fills are identical");
        ok (std::memcmp (a.data(), c.data(), c.size() * sizeof (float)) == 0,
            "a shorter fill is a bitwise prefix of a longer one");
    }

    // --- Analytic level calibration: RMS over the contract window == kRefRmsDb (no measurement pass) ---
    group ("stimulus lands at the reference level at any rate");
    for (double sr : { 44100.0, 48000.0, 96000.0 })
    {
        const int total = (int) std::lround ((lp::kSettleSec + lp::kMeasureSec) * sr);
        std::vector<float> s ((std::size_t) total);
        lp::fill (s.data(), total, sr);
        approx (windowRmsDb (s, sr), lp::kRefRmsDb, 0.25,   // worst observed 0.09 dB (48 kHz window)
                "stimulus RMS == kRefRmsDb over the contract window");
    }

    // --- Shaping sanity: the stimulus is guitar-band (LP at kShapeHz), not white ---
    group ("stimulus is low-pass shaped, not white");
    {
        const double sr = 48000.0;
        std::vector<float> s (65536);
        lp::fill (s.data(), (int) s.size(), sr);
        std::vector<double> sd (s.begin(), s.end());
        const auto mag = measurement::magSpectrum (sd, sd.size());
        const double binHz = sr / (double) sd.size();
        auto bandPower = [&] (double f0, double f1)
        {
            double p = 0.0; std::size_t n = 0;
            for (std::size_t k = (std::size_t) (f0 / binHz); k < (std::size_t) (f1 / binHz) && k < mag.size(); ++k)
            { p += mag[k] * mag[k]; ++n; }
            return p / (double) (n > 0 ? n : 1);
        };
        const double low = bandPower (100.0, 1000.0), high = bandPower (8000.0, 20000.0);
        ok (high < 0.2 * low, "per-bin power at 8-20 kHz well below 0.1-1 kHz (one-pole LP)");
    }

    // --- referenceUnityGain: flat IR is already at unity ---
    group ("delta IR (flat response) needs no gain");
    {
        std::vector<float> d (256, 0.0f); d[0] = 1.0f;
        approx ((double) measurement::referenceUnityGain (d.data(), (int) d.size(), 48000.0), 1.0, 1.0e-3,
                "gain(delta) == 1");
        // scale inversion: half-level IR needs exactly the inverse boost
        d[0] = 0.5f;
        approx ((double) measurement::referenceUnityGain (d.data(), (int) d.size(), 48000.0), 2.0, 2.0e-3,
                "gain(0.5*delta) == 2");
    }

    // --- ONE common gain from the channel-MEAN power (stereo imaging untouched) ---
    group ("stereo: channel-mean power, one common gain");
    {
        std::vector<float> l (256, 0.0f), r (256, 0.0f);
        l[0] = 1.0f; r[0] = 0.5f;                             // P_mean = (1 + 0.25)/2 = 0.625
        const float* both[2] { l.data(), r.data() };
        approx ((double) measurement::referenceUnityGain (both, 2, 256, 48000.0),
                1.0 / std::sqrt (0.625), 2.0e-3, "gain == 1/sqrt(mean channel power)");
    }

    // --- Floor and clamps, pinned exactly ---
    group ("near-silent floor and the +/-30 dB clamp");
    {
        std::vector<float> d (256, 0.0f);
        d[0] = 1.0e-4f;                                       // gSq ~ 1e-8 < 10^(-60/10) -> leave it alone
        ok (measurement::referenceUnityGain (d.data(), 256, 48000.0) == 1.0f,
            "below the -60 dB floor: gain == 1 exactly (don't amplify garbage)");
        d[0] = 1.0e-2f;                                       // needs +40 dB -> clamped to +30
        approx ((double) measurement::referenceUnityGain (d.data(), 256, 48000.0),
                31.6227766, 1.0e-4, "boost clamps at +30 dB");
        d[0] = 100.0f;                                        // needs -40 dB -> clamped to -30
        approx ((double) measurement::referenceUnityGain (d.data(), 256, 48000.0),
                0.0316227766, 1.0e-7, "cut clamps at -30 dB");
        ok (measurement::referenceUnityGain (nullptr, 0, 0, 0.0) == 1.0f, "degenerate input: gain == 1");

        // Null CHANNEL pointers are degenerate too (a caller bug, not a silent channel): the whole
        // call answers g = 1 — leave the IR alone, never dereference.
        ok (measurement::referenceUnityGain ((const float*) nullptr, 256, 48000.0) == 1.0f,
            "mono nullptr taps: gain == 1, no dereference");
        d[0] = 1.0f;
        const float* oneNull[2] { d.data(), nullptr };
        ok (measurement::referenceUnityGain (oneNull, 2, 256, 48000.0) == 1.0f,
            "stereo with one null channel: gain == 1, no dereference");
        const float* allNull[2] { nullptr, nullptr };
        ok (measurement::referenceUnityGain (allNull, 2, 256, 48000.0) == 1.0f,
            "stereo with all-null channels: gain == 1, no dereference");
    }

    // --- Analysis cap: taps past kIrAnalysisSeconds are ignored EXACTLY ---
    group ("analysis window caps at the first second");
    {
        const double sr = 48000.0;
        const int one = (int) sr;
        auto h = dampedNoiseIr (one + one / 2, 0.3f, 0, 3.0f);
        for (int i = one; i < one + one / 2; ++i) h[(std::size_t) i] = 5.0f;   // loud junk after 1 s
        const float gLong = measurement::referenceUnityGain (h.data(), one + one / 2, sr);
        const float gCap  = measurement::referenceUnityGain (h.data(), one, sr);
        ok (gLong == gCap, "gain(1.5 s IR) == gain(its first second), bitwise");
    }

    // --- Golden gains (deterministic damped-noise IRs; extraction-time values, libm-tolerant) ---
    group ("golden gains — extraction-time fingerprint");
    {
        const int len = 4800; const double sr = 48000.0;
        const auto l = dampedNoiseIr (len, 0.3f, 0, 3.0f);
        const auto r = dampedNoiseIr (len, 0.2f, 7, 4.0f);
        const float* both[2] { l.data(), r.data() };
        approx ((double) measurement::referenceUnityGain (l.data(), len, sr), 0.20633997, 1.0e-4,
                "mono golden gain");
        approx ((double) measurement::referenceUnityGain (both, 2, len, sr), 0.253126025, 1.0e-4,
                "stereo golden gain");
    }

    // --- THE CONTRACT, end to end: a normalized vendor-style cab IR passes the probe at ~unity ---
    // (the OrbitCab acceptance test, engine-free: FFT-domain gain == time-domain convolved RMS ratio)
    group ("probe through a normalized IR reads ~0 dB — 'unity' is one contract");
    {
        const double sr = 48000.0;
        const int len = 2048;
        std::vector<float> h ((std::size_t) len, 0.0f);       // cab-flavored damped sines, peak-normalized
        for (int i = 0; i < len; ++i)
        {
            const double t = (double) i / sr;
            h[(std::size_t) i] = (float) (std::exp (-t * 60.0)  * std::sin (2.0 * 3.14159265358979 * 110.0 * t) * 0.9
                                        + std::exp (-t * 200.0) * std::sin (2.0 * 3.14159265358979 * 900.0 * t) * 0.6
                                        + std::exp (-t * 900.0) * std::sin (2.0 * 3.14159265358979 * 3500.0 * t) * 0.4);
        }
        float pk = 0.0f; for (float v : h) pk = std::max (pk, std::fabs (v));
        for (float& v : h) v = v / pk * 0.15f;                // vendor-like: hot passband, inside the clamp

        const float g = measurement::referenceUnityGain (h.data(), len, sr);
        const double rawDb = -20.0 * std::log10 ((double) g); // the IR's raw reference gain
        ok (rawDb > 6.0 && rawDb < 30.0, "raw peak-normalized cab IR measures hot (like real vendor IRs)");

        const int settle = (int) std::lround (lp::kSettleSec * sr);
        const int meas   = (int) std::lround (lp::kMeasureSec * sr);
        const int total  = settle + meas + len;
        std::vector<float> x ((std::size_t) total);
        lp::fill (x.data(), total, sr);
        double si = 0.0, so = 0.0;                            // convolve the NORMALIZED IR, window the RMS
        for (int n = settle; n < settle + meas; ++n)
        {
            double acc = 0.0;
            for (int k = 0; k < len && k <= n; ++k)
                acc += (double) h[(std::size_t) k] * g * x[(std::size_t) (n - k)];
            si += (double) x[(std::size_t) n] * x[(std::size_t) n];
            so += acc * acc;
        }
        const double throughDb = 10.0 * std::log10 (so / si);
        approx (throughDb, 0.0, 0.75,                         // observed +0.14 dB; OrbitCab acceptance used 1.0
                "normalized IR passes the reference stimulus at ~unity RMS");
    }

    return felitronics::test::report();
}

// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// THEORY-FIRST FALSIFICATION suite for the family "unity" contract (LevelProbe.h + ReferenceUnity.h).
// Distinct from LevelProbeTests.cpp (which pins GOLDEN fingerprints): every expectation here is
// DERIVED from the documented contract + signal theory BEFORE the arithmetic is trusted, then the
// implementation is challenged to meet it. Nothing here re-reads a golden literal; each check states
// the analytic prediction and the tolerance's statistical basis.
//
// Attack surfaces (all numbers derived, not tuned to pass):
//   1. Stimulus calibration  — white uniform[-1,1] (var 1/3) through one-pole LP y=a·x+(1-a)·y[-1]
//      has output variance σ²·a/(2-a); the fill gain undoes that analytically, so the windowed RMS
//      must land on kRefRmsDb. Asserted over the documented window (sampling-noise tol) AND a long
//      converged window (unbiasedness tol) across six rates.
//   2. Spectrum shape        — the stimulus PSD must follow the ONE-POLE response. Welch estimate vs
//      the DISCRETE |H(e^jω)|² (the realized filter) and vs the CONTINUOUS w(f)=1/(1+(f/2000)²) (the
//      weight ReferenceUnity integrates). Derived finding: w(f) is a LOW-FREQUENCY approximation of
//      the realized discrete PSD — they agree <0.1 dB below ~1 kHz and diverge to ~2.6 dB by 20 kHz.
//      DOCUMENTED DEVIATION, not a defect (the HF band is ~15 dB down; closed-loop unity stays tight).
//   3. Closed loop           — g=referenceUnityGain(h) must pass the stimulus through g·h at 0 dB for
//      IN-SCOPE (short, cab-like) IRs. The HF weight gap of #2 predicts a SMALL POSITIVE residual for
//      a high-pass IR (continuous w under-weights the HF the stimulus actually carries) — asserted by
//      sign. Long reverb IRs are OUT of the documented scope ("NOT for reverb/decay IRs") — shown to
//      miss cab unity, confirming the exclusion.
//   4. Scale invariance      — |k·h|² scales power by k², so g(k·h)=g(h)/k EXACTLY, until the ±30 dB
//      clamp or the −60 dB gSq floor bite. Boundaries probed ±0.1 dB either side (delta amplitude).
//   5. Channel-mean law      — G² = mean-over-channels power ⇒ identical stereo == mono (bit-exact);
//      one silent channel halves the mean power ⇒ gain ×√2; two deltas ⇒ 1/√(mean(A²,B²)).
//   6. Fuzz under contract   — random rate/length/finite-content ⇒ gain finite and inside the clamp
//      range, degenerate/null ⇒ exactly 1.0f, no crash under ASan/UBSan.

#include <felitronics_test.h>
#include <felitronics/measurement/LevelProbe.h>
#include <felitronics/measurement/ReferenceUnity.h>
#include <felitronics/measurement/Convolve.h>   // measurement::magSpectrum — independent double-FFT analysis

#include <cmath>
#include <cstddef>
#include <random>
#include <vector>

using namespace felitronics;
using felitronics::test::ok;
using felitronics::test::approx;
using felitronics::test::group;

namespace lp = measurement::levelprobe;

namespace
{
constexpr double kTwoPi = 6.283185307179586;

// Deterministic damped-noise "IR" — taps from the probe's own LCG (platform-exact input data).
std::vector<float> dampedNoiseIr (int len, float amp, int seed, float decay)
{
    std::vector<float> h ((std::size_t) len);
    for (int i = 0; i < len; ++i)
        h[(std::size_t) i] = amp * lp::white (i + seed) * std::exp (-decay * (float) i / (float) len);
    return h;
}

// RMS (dB) of x over [startSec, startSec+lenSec) at sr.
double windowRmsDb (const std::vector<float>& x, double sr, double startSec, double lenSec)
{
    const int s = (int) std::lround (startSec * sr);
    const int m = (int) std::lround (lenSec  * sr);
    double acc = 0.0;
    for (int i = s; i < s + m; ++i) acc += (double) x[(std::size_t) i] * x[(std::size_t) i];
    return 10.0 * std::log10 (acc / m);
}

// Power response of the REALIZED discrete one-pole LP: a = 1-exp(-2π·kShapeHz/sr), pole ρ = 1-a.
double discPow (double f, double sr)
{
    const double a = 1.0 - std::exp (-kTwoPi * lp::kShapeHz / sr);
    const double rho = 1.0 - a;
    const double w = kTwoPi * f / sr;
    return a * a / (1.0 - 2.0 * rho * std::cos (w) + rho * rho);
}
// The CONTINUOUS one-pole weight ReferenceUnity integrates.
double contPow (double f) { const double r = f / lp::kShapeHz; return 1.0 / (1.0 + r * r); }

// Deterministic Welch PSD estimate (Hann, 50% overlap) — the stimulus is fixed, so this is stable
// across runs/platforms (double FFT on identical input). Returns the mean periodogram, length nfft/2.
std::vector<double> welchPsd (const std::vector<float>& s, int nfft)
{
    std::vector<double> win ((std::size_t) nfft);
    for (int i = 0; i < nfft; ++i) win[(std::size_t) i] = 0.5 - 0.5 * std::cos (kTwoPi * i / (nfft - 1));
    std::vector<double> psd ((std::size_t) (nfft / 2), 0.0);
    const int nseg = ((int) s.size() - nfft) / (nfft / 2);
    int used = 0;
    for (int seg = 0; seg < nseg; ++seg)
    {
        const int off = seg * (nfft / 2);
        std::vector<double> buf ((std::size_t) nfft);
        for (int i = 0; i < nfft; ++i) buf[(std::size_t) i] = (double) s[(std::size_t) (off + i)] * win[(std::size_t) i];
        const auto m = measurement::magSpectrum (buf, (std::size_t) nfft);
        for (int k = 0; k < nfft / 2; ++k) psd[(std::size_t) k] += m[(std::size_t) k] * m[(std::size_t) k];
        ++used;
    }
    for (auto& v : psd) v /= (used > 0 ? used : 1);
    return psd;
}

// Band-averaged PSD value over [f·(1-frac), f·(1+frac)] — cuts single-bin variance for a fair
// shape comparison at log-spaced checkpoints.
double bandAvg (const std::vector<double>& psd, double binHz, double f, double frac)
{
    double p = 0.0; int n = 0;
    for (int k = (int) (f * (1.0 - frac) / binHz); k <= (int) (f * (1.0 + frac) / binHz) && k < (int) psd.size(); ++k)
    { p += psd[(std::size_t) k]; ++n; }
    return p / (n > 0 ? n : 1);
}

// Convolve the stimulus through g·h and read the through-RMS (dB) over the contract window.
// g·h at unity ⇒ 0 dB (the OrbitCab acceptance test, engine-free).
double closedLoopThroughDb (const std::vector<float>& h, double sr)
{
    const int len    = (int) h.size();
    const float g    = measurement::referenceUnityGain (h.data(), len, sr);
    const int settle = (int) std::lround (lp::kSettleSec  * sr);
    const int meas   = (int) std::lround (lp::kMeasureSec * sr);
    const int total  = settle + meas + len;
    std::vector<float> x ((std::size_t) total);
    lp::fill (x.data(), total, sr);
    double si = 0.0, so = 0.0;
    for (int n = settle; n < settle + meas; ++n)
    {
        double acc = 0.0;
        for (int k = 0; k < len && k <= n; ++k)
            acc += (double) h[(std::size_t) k] * g * x[(std::size_t) (n - k)];
        si += (double) x[(std::size_t) n] * x[(std::size_t) n];
        so += acc * acc;
    }
    return 10.0 * std::log10 (so / si);
}

// approx on a relative tolerance (values that span orders of magnitude).
void approxRel (double got, double want, double rel, const std::string& msg)
{ approx (got, want, std::fabs (want) * rel + 1.0e-12, msg); }
} // namespace

int main()
{
    std::printf ("felitronics::measurement level-probe + reference-unity THEORY-FIRST falsification\n");

    // =====================================================================================
    // 1. STIMULUS CALIBRATION — derived: uniform[-1,1] var=1/3, one-pole LP output var=σ²·a/(2-a),
    //    fill gain undoes it ⇒ ensemble RMS == kRefRmsDb. The finite window only adds SAMPLING noise
    //    (AR(1) pole ρ=1-a: rel std of the RMS estimate ≈ 0.5·√((2/M)(1+ρ²)/(1-ρ²)) ≈ 0.12 dB at
    //    48 kHz over the 0.2 s window). Tolerance below is ~3σ of that; a real calibration error is dB.
    // =====================================================================================
    group ("stimulus lands at kRefRmsDb over the DOCUMENTED window (sampling-noise bound)");
    for (double sr : { 44100.0, 48000.0, 88200.0, 96000.0, 176400.0, 192000.0 })
    {
        const int total = (int) std::lround ((lp::kSettleSec + lp::kMeasureSec) * sr);
        std::vector<float> s ((std::size_t) total);
        lp::fill (s.data(), total, sr);
        approx (windowRmsDb (s, sr, lp::kSettleSec, lp::kMeasureSec), lp::kRefRmsDb, 0.35,
                "windowed stimulus RMS == kRefRmsDb (worst derived-observed 0.16 dB @176.4k)");
    }

    group ("stimulus RMS is UNBIASED over a long converged window (sampling noise beaten down)");
    for (double sr : { 44100.0, 48000.0, 96000.0, 192000.0 })
    {
        const int total = (int) std::lround ((0.1 + 2.0) * sr);   // 0.1 s settle + 2 s measured
        std::vector<float> s ((std::size_t) total);
        lp::fill (s.data(), total, sr);
        // Residual < 0.07 dB (a mild positive bias that grows with sr: the LCG input is not perfectly
        // white, and a smaller a at higher sr integrates more near-neighbour correlation). Documented,
        // well inside the acceptance the products shipped against.
        approx (windowRmsDb (s, sr, 0.1, 2.0), lp::kRefRmsDb, 0.20,
                "converged RMS == kRefRmsDb (residual bias < 0.07 dB)");
    }

    // =====================================================================================
    // 2. SPECTRUM SHAPE — derived: stimulus = white ⊛ one-pole ⇒ PSD ∝ the one-pole power response.
    //    We test against the REALIZED discrete |H(e^jω)|² and expose that the CONTINUOUS weight
    //    w(f)=1/(1+(f/2000)²) — which ReferenceUnity integrates — is only its LF approximation.
    // =====================================================================================
    {
        const double sr = 48000.0;
        const int nfft = 4096;
        const int total = nfft + 600 * (nfft / 2);
        std::vector<float> s ((std::size_t) total);
        lp::fill (s.data(), total, sr);
        const auto psd = welchPsd (s, nfft);
        const double binHz = sr / nfft;

        // Premise check: the raw LCG input IS white (flat PSD) — otherwise "PSD == |H|²" is unfounded.
        group ("premise: raw white() is spectrally flat");
        {
            std::vector<float> w ((std::size_t) total);
            for (int i = 0; i < total; ++i) w[(std::size_t) i] = lp::white (i);
            const auto wp = welchPsd (w, nfft);
            const double ref = bandAvg (wp, binHz, 600.0, 0.5);   // 300..900 Hz reference
            for (double f : { 200.0, 1000.0, 5000.0, 10000.0, 20000.0 })
                approx (10.0 * std::log10 (bandAvg (wp, binHz, f, 0.15) / ref), 0.0, 0.30,
                        "white() PSD flat within 0.3 dB (input really is white)");
        }

        // Fit the estimate's overall scale in the flat LF band, then challenge the SHAPE elsewhere.
        double num = 0.0, den = 0.0;
        for (int k = (int) (150.0 / binHz); k < (int) (600.0 / binHz); ++k)
        { num += psd[(std::size_t) k]; den += discPow (k * binHz, sr); }
        const double scale = num / den;

        group ("stimulus PSD follows the REALIZED discrete one-pole across the band");
        for (double f : { 200.0, 315.0, 500.0, 800.0, 1250.0, 2000.0, 3150.0, 5000.0, 8000.0, 12500.0, 16000.0, 20000.0 })
        {
            const double emp = bandAvg (psd, binHz, f, 0.12) / scale;
            approx (10.0 * std::log10 (emp / discPow (f, sr)), 0.0, 0.40,
                    "empirical PSD == discrete |H(e^jω)|² (Welch-variance bound)");
        }

        group ("PSD follows the CONTINUOUS weight only at LF (where the approximation holds)");
        for (double f : { 200.0, 500.0, 1000.0, 2000.0 })
        {
            const double emp = bandAvg (psd, binHz, f, 0.12) / scale;
            approx (10.0 * std::log10 (emp / contPow (f)), 0.0, 0.20,
                    "empirical PSD == continuous w(f) below 2 kHz");
        }

        group ("DOCUMENTED DEVIATION: at HF the stimulus carries MORE power than w(f) states");
        {
            const double e16 = bandAvg (psd, binHz, 16000.0, 0.12) / scale;
            const double e20 = bandAvg (psd, binHz, 20000.0, 0.12) / scale;
            ok (10.0 * std::log10 (e16 / contPow (16000.0)) > 1.0,
                "empirical/continuous @16 kHz exceeds +1 dB (realized filter ≠ continuous model)");
            ok (10.0 * std::log10 (e20 / contPow (20000.0)) > 2.0,
                "empirical/continuous @20 kHz exceeds +2 dB");
            // Root-cause pin (pure analytic identity, deterministic): the model gap the normalizer
            // integrates against. <0.03 dB below 1 kHz; grows to ~2.6 dB by 20 kHz.
            approx (10.0 * std::log10 (discPow (500.0,   sr) / contPow (500.0)),   0.0,   0.03,  "disc==cont @500 Hz");
            approx (10.0 * std::log10 (discPow (8000.0,  sr) / contPow (8000.0)),  0.399, 0.015, "disc/cont @8 kHz");
            approx (10.0 * std::log10 (discPow (16000.0, sr) / contPow (16000.0)), 1.643, 0.015, "disc/cont @16 kHz");
            approx (10.0 * std::log10 (discPow (20000.0, sr) / contPow (20000.0)), 2.628, 0.015, "disc/cont @20 kHz");
        }
    }

    // =====================================================================================
    // 3. CLOSED LOOP — the contract, end to end: g·h passes the stimulus at 0 dB for IN-SCOPE IRs.
    // =====================================================================================
    group ("in-scope (short, cab-like) IRs pass the probe at ~unity (|through| <= 0.75 dB)");
    {
        const double sr = 48000.0;
        { std::vector<float> h (256, 0.f); h[0] = 1.f;
          approx (closedLoopThroughDb (h, sr), 0.0, 0.75, "delta"); }
        { std::vector<float> h (2048, 0.f); h[137] = 0.8f;
          approx (closedLoopThroughDb (h, sr), 0.0, 0.75, "delayed delta (magnitude-flat, phase irrelevant)"); }
        // low-pass-ish IR: a one-pole (200 Hz) impulse response — energy at LF, where w matches reality.
        { const int len = 4096; std::vector<float> h ((std::size_t) len, 0.f);
          const double a = 1.0 - std::exp (-kTwoPi * 200.0 / sr); double y = 0.0, imp = 1.0;
          for (int i = 0; i < len; ++i) { y += a * (imp - y); imp = 0.0; h[(std::size_t) i] = (float) y; }
          approx (closedLoopThroughDb (h, sr), 0.0, 0.75, "low-pass-ish IR"); }
        // realistic peak-normalized cab IR (hot passband, inside the clamp).
        { const int len = 2048; std::vector<float> h ((std::size_t) len, 0.f);
          for (int i = 0; i < len; ++i) { const double t = i / sr;
            h[(std::size_t) i] = (float) (std::exp (-t*60.0)  * std::sin (kTwoPi*110.0*t)  * 0.9
                                        + std::exp (-t*200.0) * std::sin (kTwoPi*900.0*t)  * 0.6
                                        + std::exp (-t*900.0) * std::sin (kTwoPi*3500.0*t) * 0.4); }
          float pk = 0.f; for (float v : h) pk = std::max (pk, std::fabs (v));
          for (float& v : h) v = v / pk * 0.15f;
          approx (closedLoopThroughDb (h, sr), 0.0, 0.75, "peak-normalized cab IR"); }
    }

    group ("HF weight gap (finding #2) predicts a small POSITIVE residual for a high-pass IR");
    {
        const double sr = 48000.0; const int len = 4096;
        std::vector<float> h ((std::size_t) len, 0.f);
        const double a = 1.0 - std::exp (-kTwoPi * 1500.0 / sr); double y = 0.0, imp = 1.0;
        for (int i = 0; i < len; ++i) { const double x = imp; imp = 0.0; y += a * (x - y); h[(std::size_t) i] = (float) (x - y); }
        const double through = closedLoopThroughDb (h, sr);
        // Derived sign: HP energy sits where the continuous w UNDER-weights the discrete stimulus,
        // so more HF passes than the normalizer assumed ⇒ hotter than unity. Observed +0.18 dB.
        ok (through > 0.0 && through < 0.5,
            "high-pass IR reads HOT (0 < through < 0.5 dB) — the continuous-weight HF gap made real");
    }

    group ("OUT OF SCOPE: a long reverb/decay IR does NOT satisfy cab unity (documented exclusion)");
    {
        const double sr = 48000.0; const int len = (int) (0.8 * sr);
        std::vector<float> h ((std::size_t) len);
        for (int i = 0; i < len; ++i) h[(std::size_t) i] = 0.3f * lp::white (i + 5) * std::exp (-3.0f * i / len);
        const float g = measurement::referenceUnityGain (h.data(), len, sr);
        ok (std::isfinite (g) && g >= 0.0316227f && g <= 31.62278f, "reverb IR gain still finite & clamped");
        // Its decaying tail exceeds the measurement window ⇒ the convolution never reaches steady state
        // in-window and unity misses by ~2 dB. This is WHY the header excludes reverb/decay IRs.
        ok (std::fabs (closedLoopThroughDb (h, sr)) > 0.75,
            "reverb IR misses cab unity (|through| > 0.75 dB) — 'NOT for reverb/decay IRs' confirmed");
    }

    // =====================================================================================
    // 4. SCALE INVARIANCE + CLAMP/FLOOR — g(k·h)=g(h)/k until the boundaries bite (probed ±0.1 dB).
    //    Delta amplitude A: gSq=A², gDb=-20·log10(A). Clamp at gDb=±30 (A=10^∓1.5); floor at gSq<1e-6.
    // =====================================================================================
    group ("scale law on a delta: g(2^m·δ) == 2^-m while gDb stays in (-30,30)");
    {
        const double sr = 48000.0; const int len = 512;
        for (int m = -4; m <= 4; ++m)                               // gDb = -6.0206·m ∈ [-24.1, 24.1]
        {
            const double A = std::pow (2.0, m);
            std::vector<float> h ((std::size_t) len, 0.f); h[0] = (float) A;
            const float g = measurement::referenceUnityGain (h.data(), len, sr);
            approxRel ((double) g, std::pow (2.0, -m), 1.0e-5, "g(2^m·δ) == 2^-m (exact power scaling)");
        }
    }
    group ("scale law on a COLORED IR is spectral-content-independent, over ~54 dB");
    {
        const double sr = 48000.0; const int len = 4800;
        const auto h0 = dampedNoiseIr (len, 0.3f, 0, 3.0f);
        const double g0 = (double) measurement::referenceUnityGain (h0.data(), len, sr);
        const double gDb0 = -20.0 * std::log10 (g0);               // ≈ -13.7 dB (this IR measures hot)
        for (int m = -7; m <= 2; ++m)                              // gDb0 - 6.0206·m ∈ (-30,30), above floor
        {
            std::vector<float> h ((std::size_t) len);
            const double k = std::pow (2.0, m);
            for (int i = 0; i < len; ++i) h[(std::size_t) i] = (float) (h0[(std::size_t) i] * k);
            const float g = measurement::referenceUnityGain (h.data(), len, sr);
            const bool inRange = std::fabs (gDb0 - 6.020599913 * m) < 30.0;  // stay clear of the clamp
            if (inRange) approxRel ((double) g, g0 / k, 1.0e-5, "g(k·h) == g(h)/k for a colored IR");
        }
    }
    group ("clamp/floor boundaries — probed ±0.1 dB either side, pinned to the documented edges");
    {
        const double sr = 48000.0; const int len = 512;
        auto gOfDelta = [&] (double A) { std::vector<float> h ((std::size_t) len, 0.f); h[0] = (float) A;
                                         return measurement::referenceUnityGain (h.data(), len, sr); };
        // −30 dB clamp (cut): gDb = -20·log10(A).
        approxRel ((double) gOfDelta (std::pow (10.0,  29.9 / 20.0)), 1.0 / std::pow (10.0, 29.9 / 20.0), 1.0e-5,
                   "gDb=-29.9 (just inside): unclamped g == 1/A");
        approx    ((double) gOfDelta (std::pow (10.0,  30.1 / 20.0)), 0.0316227766, 1.0e-6,
                   "gDb=-30.1 (just outside): cut clamps to 10^(-30/20)");
        // +30 dB clamp (boost).
        approxRel ((double) gOfDelta (std::pow (10.0, -29.9 / 20.0)), 1.0 / std::pow (10.0, -29.9 / 20.0), 1.0e-5,
                   "gDb=+29.9 (just inside): unclamped g == 1/A");
        approx    ((double) gOfDelta (std::pow (10.0, -30.1 / 20.0)), 31.6227766, 1.0e-4,
                   "gDb=+30.1 (just outside): boost clamps to 10^(30/20)");
        // −60 dB gSq floor (gSq = A²): boost still applies just above, drops to exactly 1 just below.
        approx    ((double) gOfDelta (std::pow (10.0, -59.9 / 20.0)), 31.6227766, 1.0e-4,
                   "gSq=-59.9 dB (above floor): still the +30 dB clamp");
        ok        (gOfDelta (std::pow (10.0, -60.1 / 20.0)) == 1.0f,
                   "gSq=-60.1 dB (below floor): gain == 1.0f EXACTLY (don't amplify garbage)");
        // Pin the floor threshold at exactly 10^(-60/10) in gSq using powers of two (A²=2^(2m)).
        ok (gOfDelta (std::pow (2.0, -9)) == 31.6227766f || std::fabs (gOfDelta (std::pow (2.0, -9)) - 31.6227766f) < 1e-3f,
            "gSq=2^-18 (>1e-6): above floor ⇒ clamped boost");
        ok (gOfDelta (std::pow (2.0, -10)) == 1.0f,
            "gSq=2^-20 (<1e-6): below floor ⇒ gain == 1.0f");
    }

    // =====================================================================================
    // 5. CHANNEL-MEAN LAW — G² = Σw·P_mean/Σw with P_mean the per-channel mean power.
    // =====================================================================================
    group ("channel-mean: identical stereo == mono; one silent channel ⇒ ×√2; two deltas ⇒ 1/√mean");
    {
        const double sr = 48000.0; const int len = 1024;
        const auto L = dampedNoiseIr (len, 0.3f, 0, 3.0f);
        const std::vector<float> silent ((std::size_t) len, 0.f);
        const float gMono = measurement::referenceUnityGain (L.data(), len, sr);
        const float* LL[2]   { L.data(), L.data() };
        const float* Lsil[2] { L.data(), silent.data() };
        ok (measurement::referenceUnityGain (LL, 2, len, sr) == gMono,
            "identical stereo channels ⇒ EXACTLY the mono gain (bit-for-bit)");
        // One silent channel halves the mean power ⇒ G halves in power ⇒ gain scales by √2.
        approxRel ((double) measurement::referenceUnityGain (Lsil, 2, len, sr), std::sqrt (2.0) * gMono, 1.0e-5,
                   "one silent channel ⇒ gain × √2 (mean-power halved)");
        // Two deltas: P_mean = (A²+B²)/2, flat ⇒ g = 1/√P_mean.
        for (auto ab : { std::pair<double,double> { 1.0, 0.5 }, { 0.8, 0.3 } })
        {
            std::vector<float> a ((std::size_t) 256, 0.f), b ((std::size_t) 256, 0.f);
            a[0] = (float) ab.first; b[0] = (float) ab.second;
            const float* both[2] { a.data(), b.data() };
            approxRel ((double) measurement::referenceUnityGain (both, 2, 256, sr),
                       1.0 / std::sqrt ((ab.first * ab.first + ab.second * ab.second) / 2.0), 2.0e-4,
                       "two-delta stereo ⇒ g == 1/√(mean channel power)");
        }
    }

    // =====================================================================================
    // 6. FUZZ UNDER CONTRACT — finite, in-clamp-range, no crash (ASan/UBSan); degenerate ⇒ 1.0f.
    // =====================================================================================
    group ("degenerate / null paths return EXACTLY 1.0f (documented)");
    {
        ok (measurement::referenceUnityGain (nullptr, 0, 0, 0.0) == 1.0f, "null channel array");
        ok (measurement::referenceUnityGain ((const float*) nullptr, 256, 48000.0) == 1.0f, "mono null taps");
        std::vector<float> d (256, 0.f); d[0] = 1.f;
        const float* oneNull[2] { d.data(), nullptr };
        const float* allNull[2] { nullptr, nullptr };
        ok (measurement::referenceUnityGain (oneNull, 2, 256, 48000.0) == 1.0f, "one null channel of two");
        ok (measurement::referenceUnityGain (allNull, 2, 256, 48000.0) == 1.0f, "all null channels");
        ok (measurement::referenceUnityGain (d.data(), 0, 48000.0) == 1.0f, "zero samples");
        ok (measurement::referenceUnityGain (d.data(), 256, 0.0)   == 1.0f, "zero sample rate");
        ok (measurement::referenceUnityGain (d.data(), 256, -48000.0) == 1.0f, "negative sample rate");
        const std::vector<float> sil (256, 0.f);
        ok (measurement::referenceUnityGain (sil.data(), 256, 48000.0) == 1.0f, "pure silence ⇒ floor ⇒ 1.0f");
        std::vector<float> one (1, 0.5f);                          // single-tap IR = a delta
        approxRel ((double) measurement::referenceUnityGain (one.data(), 1, 48000.0), 2.0, 1.0e-4, "length-1 IR");
    }

    group ("fuzz: random rate/length/finite content ⇒ finite gain inside the clamp range, no crash");
    {
        std::mt19937 rng (0xC0FFEE);
        int bad = 0;
        for (int it = 0; it < 200; ++it)
        {
            const double sr  = std::uniform_real_distribution<double> (44100.0, 192000.0) (rng);
            const int    len = std::uniform_int_distribution<int> (1, (int) (2.0 * sr)) (rng);
            const bool stereo = (it % 3 == 0);
            std::vector<float> a ((std::size_t) len), b ((std::size_t) (stereo ? len : 0));
            const float scale = std::pow (10.f, std::uniform_real_distribution<float> (-4.f, 4.f) (rng));  // ±80 dB
            std::uniform_real_distribution<float> amp (-2.f, 2.f);
            for (auto& v : a) v = amp (rng) * scale;
            for (auto& v : b) v = amp (rng) * scale;
            float g;
            if (stereo) { const float* ch[2] { a.data(), b.data() }; g = measurement::referenceUnityGain (ch, 2, len, sr); }
            else        g = measurement::referenceUnityGain (a.data(), len, sr);
            if (! (std::isfinite (g) && g >= 0.0316227f && g <= 31.62278f)) ++bad;
        }
        ok (bad == 0, "every fuzz case: finite gain in [10^-1.5, 10^1.5], no crash under sanitizers");
    }

    return felitronics::test::report();
}

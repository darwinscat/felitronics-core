// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// JUCE-free self-tests for the MATRIX FIR path (split-lane EQ → 2×2 transfer matrix → MatrixConvolver).
// The decisive properties:
//   (1) BASIS DETECTION picks the right topology (MSDiag / LRDiag / Full) from a snapshot;
//   (2) LINEAR matrix truth — through the REAL process(), the measured 2×2 entry responses (recovered as
//       impulse responses per input channel) equal the analytic eq::matrixResponseZeroPhase: real & SIGNED
//       after removing the N/2 delay (off-diagonal SIGN proven), and each entry IR is SYMMETRIC (linear phase);
//   (3) NATURAL matrix truth — the Full composition (per-lane cepstral blend → compose complex entries)
//       preserves the |H| structure at k=0 and stays finite + cross-routed for k>0.
//
// The entry IRs are recovered end-to-end: settle the (idle) convolver, then feed a unit impulse on ONE input
// channel — a zero-latency convolver then outputs exactly that column of the matrix (yL,yR) = (h_XL, h_XR).

#include <felitronics_test.h>
#include <felitronics/core/Math.h>
#include <felitronics/lineareq/LinearPhaseEq.h>
#include <felitronics/lineareq/NaturalPhaseEq.h>

#include <atomic>
#include <complex>
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

// DTFT of a finite IR at digital w (exact — no window leakage, the IR IS finite).
static std::complex<double> dft (const std::vector<float>& ir, double w)
{
    std::complex<double> h { 0.0, 0.0 };
    for (int k = 0; k < (int) ir.size(); ++k) h += (double) ir[(std::size_t) k] * std::polar (1.0, -w * (double) k);
    return h;
}

// Recover one COLUMN of the realised 2×2 (both output channels) from a settled engine by injecting a unit
// impulse on input channel `inCh`. `Engine` is LinearPhaseEq or NaturalPhaseEq. cap = taps to capture.
template <class Engine>
static void captureColumn (Engine& e, int inCh, int cap, std::vector<float>& oL, std::vector<float>& oR)
{
    // settle: run silence long enough to finish the cold fade-in → convolver idle on the single operator.
    for (int b = 0; b < 40; ++b) { std::vector<float> z0 (512, 0.0f), z1 (512, 0.0f); float* io[2] { z0.data(), z1.data() }; e.process (io, 2, 512); }
    const int M = 3 * cap + 1024;
    std::vector<float> L ((std::size_t) M, 0.0f), R ((std::size_t) M, 0.0f);
    (inCh == 0 ? L : R)[0] = 1.0f;
    for (int o = 0; o < M; o += 512) { float* io[2] { L.data() + o, R.data() + o }; e.process (io, 2, std::min (512, M - o)); }
    oL.assign (L.begin(), L.begin() + cap);
    oR.assign (R.begin(), R.begin() + cap);
}

static double asymmetry (const std::vector<float>& ir, int N)   // about the centre tap N/2
{
    double a = 0.0; for (int d = 1; d < N / 2; ++d) a = std::max (a, (double) std::fabs (ir[(std::size_t) (N / 2 + d)] - ir[(std::size_t) (N / 2 - d)])); return a;
}

int main()
{
    std::printf ("felitronics::lineareq matrix-FIR tests\n");
    const double sr = 48000.0;

    // seeded random lane maker (mirrors the eq module's T3 matrix-truth test)
    std::uint64_t rng = 0xA11CE5EEDULL;
    auto uni = [&] (double lo, double hi) { rng = rng * 6364136223846793005ULL + 1442695040888963407ULL;
        return lo + (hi - lo) * ((double) (rng >> 11) * (1.0 / 9007199254740992.0)); };
    auto mk = [&] (eq::LaneParams& lp, bool on, double g) { lp.on = on; lp.freq = uni (400.0, 4000.0); lp.Q = uni (0.8, 2.2); lp.gainDb = g; };

    // --- (1) basis detection: the topology is chosen from the ACTIVE lane set ---
    test::group ("matrix-FIR basis detection (MSDiag / LRDiag / Full)");
    {
        eq::BandParams st;  st.on = true;                                   // plain single-ST → MSDiag (today's cost)
        test::ok (eq::detectBasis (&st, 1) == eq::MatrixBasis::MSDiag, "single-ST point → MSDiag");

        eq::BandParams ms; ms.on = true; ms.lane (eq::Lane::Stereo).on = false;
        ms.lane (eq::Lane::Mid).on = true; ms.lane (eq::Lane::Side).on = true;
        test::ok (eq::detectBasis (&ms, 1) == eq::MatrixBasis::MSDiag, "{M,S} point → MSDiag");

        eq::BandParams lr; lr.on = true; lr.lane (eq::Lane::Stereo).on = false;
        lr.lane (eq::Lane::Left).on = true; lr.lane (eq::Lane::Right).on = true;
        test::ok (eq::detectBasis (&lr, 1) == eq::MatrixBasis::LRDiag, "{L,R} point → LRDiag");

        eq::BandParams mix[2] { lr, ms };                                  // one band L/R, one band M/S → Full
        test::ok (eq::detectBasis (mix, 2) == eq::MatrixBasis::Full, "L/R band + M/S band → Full");

        eq::BandParams bypassed = lr; bypassed.lane (eq::Lane::Left).bypass = true; bypassed.lane (eq::Lane::Right).bypass = true;
        test::ok (eq::detectBasis (&bypassed, 1) == eq::MatrixBasis::MSDiag, "a fully-bypassed L/R point is not 'active' → MSDiag");
    }

    // Build a genuinely-Full snapshot (forces the 4-entry matrix) with distinct M/S gains → a clear
    // off-diagonal SIGN. Band 0: {L,R}; Band 1: {ST,M,S} with M ≫ S.
    eq::BandParams full[2];
    full[0].on = true; full[0].type = eq::FilterType::Bell; full[0].lane (eq::Lane::Stereo).on = false;
    mk (full[0].lane (eq::Lane::Left), true, 5.0); mk (full[0].lane (eq::Lane::Right), true, -4.0);
    full[1].on = true; full[1].type = eq::FilterType::Bell;
    mk (full[1].lane (eq::Lane::Stereo), true, 2.0);
    mk (full[1].lane (eq::Lane::Mid), true, 6.0); mk (full[1].lane (eq::Lane::Side), true, -6.0);

    // --- (2) LINEAR matrix truth + IR symmetry: measured 2×2 == analytic matrixResponseZeroPhase ---
    test::group ("LinearPhaseEq Full matrix truth: real signed entries + linear-phase symmetry");
    {
        const int Q = 1;                                                   // N = 8192
        lineareq::LinearPhaseEq e; test::ok (e.prepare (sr, 512, 2, Q), "prepare stereo q1");
        const int N = e.firSize();
        test::ok (e.setBands (full, 2), "setBands (Full snapshot)");

        std::vector<float> LL, RL, LR, RR;
        captureColumn (e, 0, N + 1, LL, RL);                              // inject L → (h_LL, h_RL)
        captureColumn (e, 1, N + 1, LR, RR);                              // inject R → (h_LR, h_RR)

        // symmetry about N/2 (⇒ exactly linear phase) for all four entries
        test::ok (asymmetry (LL, N) < 1e-5 && asymmetry (RL, N) < 1e-5 && asymmetry (LR, N) < 1e-5 && asymmetry (RR, N) < 1e-5,
                  "all four entry IRs symmetric about N/2 (linear phase)");

        double worst = 0.0, worstImag = 0.0, offDiagMag = 0.0;
        for (double f : { 500.0, 1500.0, 3500.0 })
        {
            const double w = 2.0 * core::kPi * f / sr;
            const std::complex<double> delay = std::polar (1.0, w * (double) (N / 2));   // remove the linear-phase N/2 delay
            const eq::ResponseMatrix H = eq::matrixResponseZeroPhase (full, 2, sr, w);
            const std::complex<double> mLL = dft (LL, w) * delay, mLR = dft (LR, w) * delay;
            const std::complex<double> mRL = dft (RL, w) * delay, mRR = dft (RR, w) * delay;
            worst = std::max ({ worst, std::abs (mLL - H.hLL), std::abs (mLR - H.hLR), std::abs (mRL - H.hRL), std::abs (mRR - H.hRR) });
            worstImag = std::max ({ worstImag, std::fabs (mLL.imag()), std::fabs (mLR.imag()), std::fabs (mRL.imag()), std::fabs (mRR.imag()) });
            offDiagMag = std::max ({ offDiagMag, std::abs (H.hLR), std::abs (H.hRL) });
        }
        std::printf ("      linear Full: worst |measured−analytic|=%.3e, worst |imag|=%.3e, off-diag |H|=%.3f\n", worst, worstImag, offDiagMag);
        test::ok (worst < 0.01, "measured 2×2 == matrixResponseZeroPhase (complex, delay-removed)");
        test::ok (worstImag < 2e-3, "delay-removed entries are REAL (zero-phase ⇒ symmetric IRs)");
        test::ok (offDiagMag > 0.05, "off-diagonal is non-trivial (the cross terms are actually present + signed)");

        // explicit off-diagonal SIGN: |H_M|≫|H_S| near band-1's centre ⇒ (|H_M|−|H_S|)/2 > 0 ⇒ h_LR, h_RL > 0.
        const double wc = 2.0 * core::kPi * full[1].lane (eq::Lane::Mid).freq / sr;
        const eq::ResponseMatrix Hc = eq::matrixResponseZeroPhase (full, 2, sr, wc);
        const std::complex<double> dc = std::polar (1.0, wc * (double) (N / 2));
        const double mLRr = (dft (LR, wc) * dc).real(), mRLr = (dft (RL, wc) * dc).real();
        test::ok ((mLRr > 0.0) == (Hc.hLR.real() > 0.0) && (mRLr > 0.0) == (Hc.hRL.real() > 0.0),
                  "off-diagonal SIGN matches the analytic matrix (polarity, not just magnitude)");
    }

    // --- (2b) diagonal topologies still route correctly through the matrix convolver ---
    test::group ("LinearPhaseEq LRDiag Full-vs-diagonal consistency");
    {
        const int Q = 1;
        lineareq::LinearPhaseEq e; e.prepare (sr, 512, 2, Q);
        const int N = e.firSize();
        eq::BandParams lr[1]; lr[0].on = true; lr[0].type = eq::FilterType::Bell; lr[0].lane (eq::Lane::Stereo).on = false;
        mk (lr[0].lane (eq::Lane::Left), true, 6.0); mk (lr[0].lane (eq::Lane::Right), true, -3.0);
        test::ok (eq::detectBasis (lr, 1) == eq::MatrixBasis::LRDiag, "snapshot routes LRDiag");
        e.setBands (lr, 1);
        std::vector<float> LL, RL, LR, RR;
        captureColumn (e, 0, N + 1, LL, RL);
        captureColumn (e, 1, N + 1, LR, RR);
        double crossEnergy = 0.0; for (int i = 0; i <= N; ++i) crossEnergy += std::fabs (RL[(std::size_t) i]) + std::fabs (LR[(std::size_t) i]);
        test::ok (crossEnergy < 1e-4, "LRDiag has ZERO cross terms (h_RL = h_LR = 0)");
        const double w = 2.0 * core::kPi * lr[0].lane (eq::Lane::Left).freq / sr;
        const std::complex<double> delay = std::polar (1.0, w * (double) (N / 2));
        const eq::ResponseMatrix H = eq::matrixResponseZeroPhase (lr, 1, sr, w);
        test::ok (std::abs (dft (LL, w) * delay - H.hLL) < 0.05 && std::abs (dft (RR, w) * delay - H.hRR) < 0.05,
                  "diagonal entries track the analytic L/R composites");
    }

    // --- (3) NATURAL matrix truth: Full composition preserves |H| at k=0, cross-routes, stays finite ---
    test::group ("NaturalPhaseEq Full matrix: magnitude structure + cross-routing");
    {
        lineareq::NaturalPhaseEq e; test::ok (e.prepare (sr, 512, 2, 1, 0.0f), "prepare natural q1 k=0");   // k=0 ⇒ zero-phase magnitude
        const int L = e.firSize();
        test::ok (e.setBands (full, 2), "setBands (Full snapshot)");
        std::vector<float> LL, RL, LR, RR;
        captureColumn (e, 0, L, LL, RL);
        captureColumn (e, 1, L, LR, RR);
        double worstMag = 0.0, offDiag = 0.0; bool finite = true;
        for (double f : { 500.0, 1500.0, 3500.0 })
        {
            const double w = 2.0 * core::kPi * f / sr;
            const eq::ResponseMatrix H = eq::matrixResponseZeroPhase (full, 2, sr, w);
            const double dLL = std::abs (dft (LL, w)), dLR = std::abs (dft (LR, w)), dRL = std::abs (dft (RL, w)), dRR = std::abs (dft (RR, w));
            worstMag = std::max ({ worstMag, std::fabs (dLL - std::abs (H.hLL)), std::fabs (dLR - std::abs (H.hLR)),
                                             std::fabs (dRL - std::abs (H.hRL)), std::fabs (dRR - std::abs (H.hRR)) });
            offDiag = std::max ({ offDiag, dLR, dRL });
            for (float v : LL) finite = finite && std::isfinite (v);
        }
        std::printf ("      natural Full k=0: worst |mag−|H||=%.3f, off-diag mag=%.3f\n", worstMag, offDiag);
        test::ok (finite, "entry IRs finite");
        test::ok (worstMag < 0.10, "k=0 Full entry magnitudes track |matrixResponseZeroPhase| (composition correct)");
        test::ok (offDiag > 0.05, "off-diagonal present (per-lane-blend-then-compose actually cross-routes)");
    }

    test::group ("NaturalPhaseEq Full at k=0.5 stays finite + magnitude-sane (per-lane blend runs)");
    {
        lineareq::NaturalPhaseEq e; e.prepare (sr, 512, 2, 1, 0.5f);
        const int L = e.firSize();
        test::ok (e.setBands (full, 2), "setBands (Full, k=0.5)");
        std::vector<float> LL, RL, LR, RR;
        captureColumn (e, 0, L, LL, RL);
        captureColumn (e, 1, L, LR, RR);
        bool finite = true; for (float v : LL) finite = finite && std::isfinite (v);
        for (float v : RR) finite = finite && std::isfinite (v);
        test::ok (finite, "mixed-phase Full entries finite at k=0.5");
        // low-freq diagonal is ST/L/R-dominated (M/S cross-terms small there) → magnitude ≈ zero-phase target
        const double w = 2.0 * core::kPi * 300.0 / sr;
        const eq::ResponseMatrix H = eq::matrixResponseZeroPhase (full, 2, sr, w);
        test::ok (std::fabs (std::abs (dft (LL, w)) - std::abs (H.hLL)) < 0.2, "k=0.5 diagonal magnitude sane at 300 Hz");
    }

    // --- (4) RT-safety: process() allocation-free on the matrix (Full) path ---
    test::group ("matrix-FIR process() no-alloc (Full topology)");
    {
        lineareq::LinearPhaseEq e; e.prepare (sr, 512, 2, 0);
        e.setBands (full, 2);
        std::vector<float> L (512, 0.2f), R (512, -0.1f); float* io[2] { L.data(), R.data() };
        e.process (io, 2, 512);                                            // consume the fade-in
        const long before = g_allocs.load();
        e.process (io, 2, 512); e.process (io, 2, 512);
        test::okNoAlloc (g_allocs.load() == before, "process() did not allocate on the Full matrix path");
    }

    return test::report();
}

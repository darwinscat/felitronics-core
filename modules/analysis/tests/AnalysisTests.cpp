// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa. Part of felitronics-core — see LICENSE.

// JUCE-free self-tests for felitronics::analysis::SpectrumTap (the lock-free SPSC handshake).

#include <felitronics_test.h>
#include <felitronics/analysis/SpectrumTap.h>
#include <felitronics/analysis/KWeightingFilter.h>
#include <felitronics/analysis/LoudnessMeter.h>
#include <felitronics/analysis/CorrelationMeter.h>
#include <felitronics/core/Math.h>

#include <cmath>
#include <cstdio>
#include <vector>

using namespace felitronics;

// Different orders are distinct, correctly-sized types (Law 5: configurable size without a fork).
static_assert (analysis::SpectrumTapT<10>::kSize == 1024, "order 10 -> 1024");
static_assert (analysis::SpectrumTap::kSize == 2048,       "default -> 2048");

int main()
{
    std::printf ("felitronics::analysis tests\n");
    const int N = analysis::kSpectrumFftSize;
    static float dst[analysis::kSpectrumFftSize];   // static: keep the 8 KB off the stack

    // --- SPSC handshake: a frame becomes ready only on the wrap push, and tryPull re-arms ---
    test::group ("SpectrumTap SPSC handshake");
    {
        analysis::SpectrumTap tap;
        test::ok (! tap.tryPull (dst), "no frame before a full window");

        for (int i = 0; i < N; ++i) tap.push ((float) i);   // fills fifo[0..N-1], idx = N
        test::ok (! tap.tryPull (dst), "not ready until the wrap push");

        tap.push (-1.0f);                                    // wrap → snapshots [0..N-1]
        test::ok (tap.tryPull (dst), "frame ready after the wrap push");
        test::approx (dst[0],     0.0,            1e-9, "data[0] == 0");
        test::approx (dst[N - 1], (double) (N - 1), 1e-9, "data[N-1] == N-1");
        test::ok (! tap.tryPull (dst), "re-armed (consumed) after pull");
    }

    // --- no overwrite while a frame is unread: the reader sees the OLD frame, not the new fill ---
    test::group ("SpectrumTap preserves an unread frame");
    {
        analysis::SpectrumTap t2;
        for (int i = 0; i < N; ++i) t2.push (1.0f);
        t2.push (1.0f);                                      // ready, data == all 1.0
        for (int i = 0; i < N; ++i) t2.push (2.0f);          // another window, but unread → snapshot skipped
        test::ok (t2.tryPull (dst), "frame still ready");
        test::approx (dst[0], 1.0, 1e-9, "unread frame preserved (1.0, not overwritten by 2.0)");
    }

    // --- reset() clears the cursor + ready flag ---
    test::group ("SpectrumTap reset");
    {
        analysis::SpectrumTap t3;
        for (int i = 0; i < N + 1; ++i) t3.push (3.0f);      // make a frame ready
        t3.reset();
        test::ok (! t3.tryPull (dst), "reset clears ready");
    }

    // --- BS.1770 K-weighting: canonical 48 kHz coefficients + recompute at 44.1 k ---
    test::group ("KWeightingFilter canonical 48 kHz coefficients");
    {
        analysis::KWeightingFilter kw; kw.prepare (48000.0, 2);
        const auto s = kw.shelfCoeffs();
        test::approx (s.b0,  1.5351248595869702, 1e-9, "shelf b0");
        test::approx (s.b1, -2.6916961894063807, 1e-9, "shelf b1");
        test::approx (s.b2,  1.1983928108528501, 1e-9, "shelf b2");
        test::approx (s.a1, -1.6906592931824103, 1e-9, "shelf a1");
        test::approx (s.a2,  0.7324807742158501, 1e-9, "shelf a2");
        const auto h = kw.highpassCoeffs();
        test::approx (h.a1, -1.9900474548339797, 1e-9, "hp a1");
        test::approx (h.a2,  0.9900722503662099, 1e-9, "hp a2");
        analysis::KWeightingFilter kw2; kw2.prepare (44100.0, 2);
        test::ok (std::fabs (kw2.shelfCoeffs().b0 - s.b0) > 1e-4, "coeffs recomputed at 44.1 k (not hardcoded)");
    }

    // --- correlation: identical → +1, inverted → -1, quadrature → ~0 ---
    test::group ("CorrelationMeter");
    {
        const double sr = 48000.0, f = 500.0;
        analysis::CorrelationMeter cm; cm.prepare (sr, 50.0);
        for (int i = 0; i < (int) sr; ++i) { const float v = (float) std::sin (2.0 * core::kPi * f * i / sr); cm.process (v, v); }
        test::approx (cm.correlation(), 1.0, 0.01, "identical → +1");
        cm.reset();
        for (int i = 0; i < (int) sr; ++i) { const float v = (float) std::sin (2.0 * core::kPi * f * i / sr); cm.process (v, -v); }
        test::approx (cm.correlation(), -1.0, 0.01, "inverted → -1");
        cm.reset();
        for (int i = 0; i < (int) sr; ++i) { const float a = (float) std::sin (2.0 * core::kPi * f * i / sr), b = (float) std::cos (2.0 * core::kPi * f * i / sr); cm.process (a, b); }
        test::ok (std::fabs (cm.correlation()) < 0.1, "quadrature → ~0");
    }

    // --- loudness: a steady tone has momentary ≈ short-term ≈ integrated (gating doesn't skew it) ---
    test::group ("LoudnessMeter steady-tone consistency");
    {
        const double sr = 48000.0; const int n = (int) (5.0 * sr);
        std::vector<float> L (n), R (n);
        for (int i = 0; i < n; ++i) { const float v = (float) (0.5 * std::sin (2.0 * core::kPi * 1000.0 * i / sr)); L[i] = v; R[i] = v; }
        analysis::LoudnessMeter lm; lm.prepare (sr, 2, 10.0);
        const float* ch[2] { L.data(), R.data() };
        lm.process (ch, 2, n);
        const double m = lm.momentaryLufs(), st = lm.shortTermLufs(), ig = lm.integratedLufs();
        test::ok (std::isfinite (m) && m > -60.0, "momentary finite + sensible");
        test::approx (ig, m,  0.3, "integrated ≈ momentary for a steady tone");
        test::approx (st, m,  0.3, "short-term ≈ momentary for a steady tone");
    }

    return test::report();
}

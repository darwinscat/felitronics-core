// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// JUCE-free self-tests for MixedPhaseFir — the cepstral phase-blend FIR builder behind "Natural phase".
// Reference-NULL: render a magnitude, transform the result back, and check the invariants exactly:
//   • magnitude is preserved EXACTLY (it's set by construction): |DFT(h)| == target for every k.
//   • phase is LINEAR in k: φ(k) == k·φ_min  (so φ(0)=0 zero-phase, φ(0.5)=½·φ(1), φ(1)=φ_min).
//   • k=0 → zero-phase ⇒ the impulse is even (symmetric); k=1 → minimum-phase ⇒ causal (no pre-ring).
//   • pre-ring (anti-causal energy) is monotonically killed as k: 0 → 1.

#include <felitronics_test.h>
#include <felitronics/lineareq/MixedPhaseFir.h>
#include <felitronics/core/Fft.h>

#include <cmath>
#include <cstdio>
#include <vector>

using namespace felitronics;

int main()
{
    std::printf ("felitronics::lineareq MixedPhaseFir tests\n");

    const int D = 4096, H = D / 2;
    core::fft::DefaultRealFft fft; fft.prepare (D);
    const int specF = core::fft::DefaultRealFft::spectrumFloats (D);

    // A smooth magnitude bump (a +6 dB-ish resonance) + a gentle shelf tilt — a realistic, non-trivial |H|.
    std::vector<float> mag ((std::size_t) (H + 1), 1.0f);
    const double b0 = 0.18 * H, w = 0.05 * H;
    for (int b = 0; b <= H; ++b)
    {
        const double bump = 6.0 * std::exp (-0.5 * ((b - b0) / w) * ((b - b0) / w));   // dB
        const double tilt = -2.0 * (double) b / (double) H;                            // gentle HF cut, dB
        mag[(std::size_t) b] = (float) std::pow (10.0, (bump + tilt) / 20.0);
    }

    lineareq::MixedPhaseFir<> mp; test::ok (mp.prepare (D), "prepare");

    // Forward-transform an impulse and read its magnitude + phase per bin (skip DC/Nyquist where φ≡0).
    std::vector<float> spec ((std::size_t) specF, 0.0f);
    auto analyse = [&] (const float* h, std::vector<float>& outMag, std::vector<float>& outPhase)
    {
        std::vector<float> hh (h, h + D);
        fft.forward (hh.data(), spec.data());
        outMag.assign ((std::size_t) (H + 1), 0.0f);
        outPhase.assign ((std::size_t) (H + 1), 0.0f);
        outMag[0] = std::fabs (spec[0]); outMag[(std::size_t) H] = std::fabs (spec[1]);
        for (int b = 1; b < H; ++b)
        {
            const float re = spec[(std::size_t) (2 * b)], im = spec[(std::size_t) (2 * b + 1)];
            outMag[(std::size_t) b]   = std::sqrt (re * re + im * im);
            outPhase[(std::size_t) b] = std::atan2 (im, re);
        }
    };
    auto antiCausalRatio = [&] (const float* h)   // Σ|h[n]| over the anti-causal half ÷ total
    {
        double anti = 0.0, all = 0.0;
        for (int n = 0; n < D; ++n) { all += std::fabs (h[n]); if (n > H) anti += std::fabs (h[n]); }
        return all > 0.0 ? anti / all : 0.0;
    };

    std::vector<float> m0, p0, m5, p5, m1, p1;
    // build each k into its own copy (build() returns the shared internal buffer, so snapshot before reuse)
    std::vector<float> H0, H5, H1;
    { const float* h = mp.build (mag.data(), 0.0f); H0.assign (h, h + D); }
    { const float* h = mp.build (mag.data(), 0.5f); H5.assign (h, h + D); }
    { const float* h = mp.build (mag.data(), 1.0f); H1.assign (h, h + D); }
    analyse (H0.data(), m0, p0);
    analyse (H5.data(), m5, p5);
    analyse (H1.data(), m1, p1);

    // --- magnitude preserved exactly for every k (it's set by construction) ---
    test::group ("MixedPhaseFir preserves the target magnitude for all k");
    {
        double e0 = 0.0, e5 = 0.0, e1 = 0.0;
        for (int b = 0; b <= H; ++b)
        {
            e0 = std::max (e0, (double) std::fabs (m0[(std::size_t) b] - mag[(std::size_t) b]));
            e5 = std::max (e5, (double) std::fabs (m5[(std::size_t) b] - mag[(std::size_t) b]));
            e1 = std::max (e1, (double) std::fabs (m1[(std::size_t) b] - mag[(std::size_t) b]));
        }
        test::ok (e0 < 1e-4, "k=0 magnitude == target");
        test::ok (e5 < 1e-4, "k=0.5 magnitude == target");
        test::ok (e1 < 1e-4, "k=1 magnitude == target");
    }

    // --- phase is LINEAR in k: φ(0)=0, φ(0.5)=½·φ_min, φ(1)=φ_min ---
    test::group ("MixedPhaseFir phase blends linearly: phi = k·phi_min");
    {
        double zero = 0.0, half = 0.0, phiMax = 0.0;
        for (int b = 1; b < H; ++b)
        {
            zero   = std::max (zero,   (double) std::fabs (p0[(std::size_t) b]));
            half   = std::max (half,   (double) std::fabs (p5[(std::size_t) b] - 0.5f * p1[(std::size_t) b]));
            phiMax = std::max (phiMax, (double) std::fabs (p1[(std::size_t) b]));
        }
        test::ok (zero < 1e-4, "k=0 → zero phase");
        test::ok (half < 1e-4, "k=0.5 phase == half of the minimum phase");
        test::ok (phiMax > 0.1, "k=1 has real (minimum) phase to blend against");
    }

    // --- impulse symmetry / causality at the endpoints ---
    test::group ("MixedPhaseFir endpoints: k=0 symmetric, k=1 causal");
    {
        double asym = 0.0;
        for (int n = 1; n < H; ++n) asym = std::max (asym, (double) std::fabs (H0[(std::size_t) n] - H0[(std::size_t) (D - n)]));
        test::ok (asym < 1e-5, "k=0 impulse is even (zero-phase ⇒ symmetric)");
        test::ok (antiCausalRatio (H1.data()) < 0.02, "k=1 impulse is causal (minimum-phase ⇒ ~no pre-ring)");
    }

    // --- pre-ringing is monotonically reduced as the phase blends toward minimum ---
    test::group ("MixedPhaseFir pre-ringing shrinks as k grows");
    {
        const double a0 = antiCausalRatio (H0.data());
        const double a5 = antiCausalRatio (H5.data());
        const double a1 = antiCausalRatio (H1.data());
        test::ok (a0 > a5 + 1e-4, "k=0 pre-ring > k=0.5 pre-ring");
        test::ok (a5 > a1 + 1e-4, "k=0.5 pre-ring > k=1 pre-ring");
    }

    return test::report();
}

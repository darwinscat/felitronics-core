// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// JUCE-free self-tests for the FFT seam (felitronics::core::fft::ScalarRadix2Real reference backend):
// round-trip identity, packed bins vs a brute-force DFT, known signals, and the spectral multiply-add
// (== circular convolution) the partitioned convolver is built on.

#include <felitronics_test.h>
#include <felitronics/core/Fft.h>
#include <felitronics/core/Math.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <vector>

using namespace felitronics;
using core::fft::ScalarRadix2Real;

// Brute-force DFT bin k (double) — the analytic reference.
static std::complex<double> dftBin (const std::vector<float>& x, int k)
{
    const int n = (int) x.size();
    std::complex<double> acc (0.0, 0.0);
    for (int i = 0; i < n; ++i)
    {
        const double ang = -2.0 * core::kPi * k * i / n;
        acc += (double) x[i] * std::complex<double> (std::cos (ang), std::sin (ang));
    }
    return acc;
}

int main()
{
    std::printf ("felitronics::core::fft tests\n");

    // --- round-trip identity across sizes ---
    test::group ("FFT round-trip identity (inverse(forward(x)) == x)");
    for (int order = 4; order <= 12; ++order)
    {
        const int n = 1 << order;
        ScalarRadix2Real fft;
        test::ok (fft.prepare (n), "prepare pow2");
        std::vector<float> x (n), y (n), spec (n);
        for (int i = 0; i < n; ++i) x[i] = std::sin (0.1f * i) + 0.3f * std::cos (0.027f * i * i);
        fft.forward (x.data(), spec.data());
        fft.inverse (spec.data(), y.data());
        double maxErr = 0.0;
        for (int i = 0; i < n; ++i) maxErr = std::max (maxErr, (double) std::fabs (y[i] - x[i]));
        test::ok (maxErr < 1e-3, "round-trip within tol");
    }

    // --- prepare guards ---
    test::group ("FFT prepare guards");
    {
        ScalarRadix2Real f;
        test::ok (! f.prepare (100), "rejects non-pow2");
        test::ok (! f.prepare (2),   "rejects < 4");
        test::ok (  f.prepare (16),  "accepts 16");
    }

    // --- packed bins vs brute DFT ---
    test::group ("FFT packed bins == brute DFT");
    {
        const int n = 64;
        ScalarRadix2Real fft; fft.prepare (n);
        std::vector<float> x (n), spec (n);
        for (int i = 0; i < n; ++i)
            x[i] = std::sin (2 * core::kPi * 5 * i / n) + 0.5f * std::cos (2 * core::kPi * 11 * i / n);
        fft.forward (x.data(), spec.data());
        test::approx (spec[0], dftBin (x, 0).real(),     1e-2, "DC.re");
        test::approx (spec[1], dftBin (x, n / 2).real(), 1e-2, "Nyquist.re");
        for (int k = 1; k < n / 2; ++k)
        {
            const auto X = dftBin (x, k);
            test::approx (spec[2 * k],     X.real(), 1e-2, "Re[k]");
            test::approx (spec[2 * k + 1], X.imag(), 1e-2, "Im[k]");
        }
    }

    // --- known signals: DC + Nyquist concentrate in one bin ---
    test::group ("FFT known signals (DC / Nyquist)");
    {
        const int n = 32;
        ScalarRadix2Real fft; fft.prepare (n);
        std::vector<float> x (n, 1.0f), spec (n);
        fft.forward (x.data(), spec.data());
        test::approx (spec[0], (double) n, 1e-3, "DC signal -> spec[0]==N");
        test::approx (spec[1], 0.0,        1e-3, "DC signal -> Nyquist 0");
        for (int i = 0; i < n; ++i) x[i] = (i & 1) ? -1.0f : 1.0f;
        fft.forward (x.data(), spec.data());
        test::approx (spec[1], (double) n, 1e-3, "Nyquist signal -> spec[1]==N");
        test::approx (spec[0], 0.0,        1e-3, "Nyquist signal -> DC 0");
    }

    // --- spectralMultiplyAdd == circular convolution (the convolver's core op) ---
    test::group ("FFT spectralMultiplyAdd == circular convolution");
    {
        const int n = 64;
        ScalarRadix2Real fft; fft.prepare (n);
        std::vector<float> a (n, 0.0f), b (n, 0.0f), fa (n), fb (n), acc (n, 0.0f), y (n);
        for (int i = 0; i < n; ++i) { a[i] = std::sin (0.3f * i); b[i] = (i < 8) ? (1.0f - i * 0.1f) : 0.0f; }
        fft.forward (a.data(), fa.data());
        fft.forward (b.data(), fb.data());
        fft.spectralMultiplyAdd (fa.data(), fb.data(), acc.data());   // acc starts at 0
        fft.inverse (acc.data(), y.data());
        double maxErr = 0.0;
        for (int m = 0; m < n; ++m)
        {
            double s = 0.0;
            for (int i = 0; i < n; ++i) s += (double) a[i] * b[((m - i) % n + n) % n];
            maxErr = std::max (maxErr, std::fabs (s - y[m]));
        }
        test::ok (maxErr < 1e-3, "FFT product, inverted, equals circular convolution");
    }

    return test::report();
}

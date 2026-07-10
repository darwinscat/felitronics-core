// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

//==============================================================================
// felitronics::blend — pure blend DSP kernels (JUCE-free, headless-testable). Promoted VERBATIM from
// OrbitCapture's app-local core/BlendKernels.h (de-monolith step 2 → promoted step 4; the only change
// there was juce::MathConstants<double>::pi → oc::kPi, a bit-identical π; here oc:: → felitronics::core::).
// No widget or file access — these are pure functions of std::vector<float> + scalars.
// One RBJ biquad + one-pole HP/LP, a TRUE-Butterworth slope cascade, a Hilbert all-pass phase rotation,
// and a fractional-sample time shift — the guts of the mic-blend engine.
//==============================================================================

#include <felitronics/core/Math.h>        // core::kPi
#include <felitronics/core/OfflineFft.h>  // core::offline::fftInplace, core::offline::nextPow2

#include <complex>
#include <vector>
#include <cmath>
#include <cstddef>

namespace felitronics::blend
{

// One RBJ 2nd-order (12 dB/oct) section with an explicit Q, applied in place to a mic's IR. Causal —
// adds the real filter phase near the cutoff. Q defaults to Butterworth (1/sqrt2).
inline void biquadInplace(std::vector<float>& x, double sr, double fc, bool highpass, double Q = 0.7071067811865476) {
    if (x.empty() || fc <= 0.0 || fc >= sr * 0.5) return;
    const double w0 = 2.0 * core::kPi * fc / sr;
    const double c = std::cos(w0), s = std::sin(w0), alpha = s / (2.0 * Q);
    double b0, b1, b2;
    if (highpass) { b0 = (1.0 + c) * 0.5; b1 = -(1.0 + c); b2 = (1.0 + c) * 0.5; }
    else          { b0 = (1.0 - c) * 0.5; b1 =  (1.0 - c); b2 = (1.0 - c) * 0.5; }
    const double a0 = 1.0 + alpha, a1 = -2.0 * c, a2 = 1.0 - alpha;
    b0 /= a0; b1 /= a0; b2 /= a0;
    const double na1 = a1 / a0, na2 = a2 / a0;
    double z1 = 0.0, z2 = 0.0;                                        // transposed Direct Form II
    for (auto& v : x) {
        const double in = v, out = b0 * in + z1;
        z1 = b1 * in - na1 * out + z2;
        z2 = b2 * in - na2 * out;
        v = (float)out;
    }
}

// 1st-order (6 dB/oct) one-pole HP/LP, applied in place.
inline void onePoleInplace(std::vector<float>& x, double sr, double fc, bool highpass) {
    if (x.empty() || fc <= 0.0 || fc >= sr * 0.5) return;
    const double dt = 1.0 / sr, rc = 1.0 / (2.0 * core::kPi * fc);
    if (highpass) {
        const double a = rc / (rc + dt); double yp = 0.0, xp = 0.0;
        for (auto& v : x) { const double in = v, y = a * (yp + in - xp); yp = y; xp = in; v = (float)y; }
    } else {
        const double a = dt / (rc + dt); double yp = 0.0;
        for (auto& v : x) { const double in = v; yp += a * (in - yp); v = (float)yp; }
    }
}

// A TRUE Butterworth HP/LP at the chosen slope. 6 dB/oct = one pole; else order N = slope/6 factored
// into N/2 biquad sections with the PROPER per-section Q's (Q_k = 1/(2·cos((2k+1)π/2N))) so the
// response is −3 dB at fc, matching the drawn EQ curve.
inline void applyBlendSlope(std::vector<float>& v, double sr, double fc, bool hp, int slopeDbOct) {
    if (slopeDbOct <= 0) return;
    if (slopeDbOct == 6) { onePoleInplace(v, sr, fc, hp); return; }
    if ((slopeDbOct / 6) % 2 != 0) return;   // only even-order cascades {12,24,36,48,72,96}; an odd order (e.g. 18) can't be factored into biquads
    const int order = slopeDbOct / 6;                                 // poles: 12→2, 24→4, 36→6, 48→8, 72→12, 96→16
    for (int k = 0; k < order / 2; ++k) {                             // all slopes ≥12 are even order
        const double theta = (2.0 * k + 1.0) * core::kPi / (2.0 * order);
        biquadInplace(v, sr, fc, hp, 1.0 / (2.0 * std::cos(theta)));   // Butterworth section Q
    }
}

// Broadband phase rotation by `degrees` (all-pass): y = x·cosθ − H(x)·sinθ, H = Hilbert transform.
// θ=±180° == polarity flip; θ=0 == no-op. Frequency-INDEPENDENT phase (not a delay).
inline void rotatePhase(std::vector<float>& x, double degrees) {
    if (x.empty() || std::abs(degrees) < 0.5) return;
    const double th = degrees * core::kPi / 180.0, c = std::cos(th), s = std::sin(th);
    const std::size_t n = core::offline::nextPow2(x.size());
    std::vector<std::complex<double>> X(n, std::complex<double>(0.0, 0.0));
    for (std::size_t i = 0; i < x.size(); ++i) X[i] = x[(size_t)i];
    core::offline::fftInplace(X, -1);
    X[0] = std::complex<double>(0.0, 0.0);                         // Hilbert: +f *(-j), -f *(+j), DC/Nyquist -> 0
    const std::size_t half = n / 2;
    for (std::size_t k = 1; k < half; ++k) { X[k] *= std::complex<double>(0.0, -1.0); X[n - k] *= std::complex<double>(0.0, 1.0); }
    X[half] = std::complex<double>(0.0, 0.0);
    core::offline::fftInplace(X, +1);
    const double inv = 1.0 / (double)n;
    for (std::size_t i = 0; i < x.size(); ++i) x[(size_t)i] = (float)(c * (double)x[(size_t)i] - s * X[i].real() * inv);
}

// Fractional-sample time-shift (linear interpolation); +d delays later, −d pulls earlier. Pure delay:
// does NOT change a single mic's magnitude curve — only how it combs with others.
inline std::vector<float> shiftFrac(const std::vector<float>& x, double d) {
    if (std::abs(d) < 1e-4 || x.empty()) return x;
    const int n = (int) x.size();
    std::vector<float> y((size_t) n, 0.0f);
    for (int i = 0; i < n; ++i) {
        const double src = (double) i - d;
        const int i0 = (int) std::floor(src);
        const double f = src - (double) i0;
        const double a = (i0 >= 0 && i0 < n) ? (double) x[(size_t) i0] : 0.0;
        const double b = (i0 + 1 >= 0 && i0 + 1 < n) ? (double) x[(size_t) (i0 + 1)] : 0.0;
        y[(size_t) i] = (float) (a + f * (b - a));
    }
    return y;
}

} // namespace felitronics::blend

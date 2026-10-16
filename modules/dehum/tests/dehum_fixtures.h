// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.
//
// Fixtures shared by the dehum suites: the hum generator, the measuring oracles and the exact transfer.
//
// THE ORACLES. Every level here is a LEAST-SQUARES FIT of sinusoids at KNOWN frequencies over a window —
// exact for a stationary sum of those sinusoids, whatever the window length, so no window, leakage or
// bin grid stands between the signal and the number (the private fixture rules 2 and 3 are about peak and
// zero-padding grids; a fit at the known frequency has neither). Ratios are taken against the INPUT fitted
// the same way, and the notch-shape test cross-checks the fit with an RMS ratio, an oracle of a different
// construction. The reference for a tone's attenuation is the canceller's exact multi-notch transfer,
// computed independently below (exactGainDb), not a value read off the code.
#pragma once

#include <felitronics/core/Math.h>
#include <felitronics/dehum/HumRemover.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <vector>

namespace dehumtest
{
namespace core = felitronics::core;
constexpr double kPi = core::kPi;

struct Lcg
{
    std::uint32_t s;
    double uniform() { s = s * 1664525u + 1013904223u; return (double) s / 4294967296.0; }   // [0, 1)
    float  bipolar() { return (float) (2.0 * uniform() - 1.0); }
};

inline double dbOf (double ratio) { return 20.0 * std::log10 (std::max (ratio, 1e-300)); }

// sin(2π·f·n/fs + ph) with the phase reduced exactly where f·n is an integer (integer f) — long fixtures
// keep full precision.
inline double sinAt (double f, std::int64_t n, double fs, double ph)
{
    const double cyc = std::fmod (f * (double) n, fs) / fs;
    return std::sin (2.0 * kPi * cyc + ph);
}

// Mains hum: K harmonics of f0, amplitude a1/k, deterministic random phases.
struct Hum
{
    double f0 = 50.0, a1 = 0.0316227766, fs = 48000.0;
    int    k = 24;
    double ph[64] {};
    void seed (std::uint32_t s) { Lcg r { s }; for (double& p : ph) p = 2.0 * kPi * r.uniform(); }
    double at (std::int64_t n) const
    {
        double v = 0.0;
        for (int h = 1; h <= k; ++h) v += a1 / h * sinAt (f0 * h, n, fs, ph[h - 1]);
        return v;
    }
};

// 10 ms raised-cosine taper at both ends (fixture rule 1).
inline void taper (std::vector<float>& x, double fs, double ms = 10.0)
{
    const int m = std::min ((int) (ms * 1e-3 * fs), (int) x.size() / 2);
    for (int i = 0; i < m; ++i)
    {
        const float w = (float) (0.5 - 0.5 * std::cos (kPi * (i + 0.5) / m));
        x[(std::size_t) i] *= w;
        x[x.size() - 1 - (std::size_t) i] *= w;
    }
}

// Solve A·x = b in place (Gaussian elimination, partial pivoting). n <= 32.
inline void solve (std::vector<double>& A, std::vector<double>& b, int n)
{
    for (int c = 0; c < n; ++c)
    {
        int piv = c;
        for (int r = c + 1; r < n; ++r) if (std::fabs (A[(std::size_t) (r * n + c)]) > std::fabs (A[(std::size_t) (piv * n + c)])) piv = r;
        if (piv != c) { for (int j = 0; j < n; ++j) std::swap (A[(std::size_t) (c * n + j)], A[(std::size_t) (piv * n + j)]); std::swap (b[(std::size_t) c], b[(std::size_t) piv]); }
        for (int r = c + 1; r < n; ++r)
        {
            const double f = A[(std::size_t) (r * n + c)] / A[(std::size_t) (c * n + c)];
            for (int j = c; j < n; ++j) A[(std::size_t) (r * n + j)] -= f * A[(std::size_t) (c * n + j)];
            b[(std::size_t) r] -= f * b[(std::size_t) c];
        }
    }
    for (int r = n - 1; r >= 0; --r)
    {
        double s = b[(std::size_t) r];
        for (int j = r + 1; j < n; ++j) s -= A[(std::size_t) (r * n + j)] * b[(std::size_t) j];
        b[(std::size_t) r] = s / A[(std::size_t) (r * n + r)];
    }
}

// Joint least-squares amplitudes of sinusoids at `freqs` over x[from, from+len).
inline std::vector<double> fitAmps (const std::vector<float>& x, std::int64_t from, std::int64_t len,
                                    const std::vector<double>& freqs, double fs)
{
    const int P = (int) freqs.size(), n = 2 * P;
    std::vector<double> A ((std::size_t) (n * n), 0.0), b ((std::size_t) n, 0.0), basis ((std::size_t) n);
    std::vector<std::complex<double>> rot ((std::size_t) P), cur ((std::size_t) P);
    for (int p = 0; p < P; ++p)
    {
        const double w = 2.0 * kPi * freqs[(std::size_t) p] / fs;
        rot[(std::size_t) p] = std::polar (1.0, w);
        cur[(std::size_t) p] = std::polar (1.0, std::fmod (w * (double) from, 2.0 * kPi));
    }
    for (std::int64_t i = 0; i < len; ++i)
    {
        if ((i & 4095) == 0)   // re-seed the recursive phasors so their drift never reaches the fit
            for (int p = 0; p < P; ++p)
                cur[(std::size_t) p] = std::polar (1.0, std::fmod (2.0 * kPi * freqs[(std::size_t) p] / fs * (double) (from + i), 2.0 * kPi));
        for (int p = 0; p < P; ++p)
        {
            basis[(std::size_t) (2 * p)]     = cur[(std::size_t) p].real();
            basis[(std::size_t) (2 * p + 1)] = cur[(std::size_t) p].imag();
            cur[(std::size_t) p] *= rot[(std::size_t) p];
        }
        const double v = (double) x[(std::size_t) (from + i)];
        for (int r = 0; r < n; ++r)
        {
            b[(std::size_t) r] += basis[(std::size_t) r] * v;
            for (int c = r; c < n; ++c) A[(std::size_t) (r * n + c)] += basis[(std::size_t) r] * basis[(std::size_t) c];
        }
    }
    for (int r = 0; r < n; ++r) for (int c = 0; c < r; ++c) A[(std::size_t) (r * n + c)] = A[(std::size_t) (c * n + r)];
    solve (A, b, n);
    std::vector<double> amps ((std::size_t) P);
    for (int p = 0; p < P; ++p) amps[(std::size_t) p] = std::hypot (b[(std::size_t) (2 * p)], b[(std::size_t) (2 * p + 1)]);
    return amps;
}

inline double fitAmp (const std::vector<float>& x, std::int64_t from, std::int64_t len, double f, double fs)
{
    return fitAmps (x, from, len, { f }, fs)[0];
}

inline double rms (const std::vector<float>& x, std::int64_t from, std::int64_t len)
{
    double e = 0.0;
    for (std::int64_t i = 0; i < len; ++i) e += (double) x[(std::size_t) (from + i)] * (double) x[(std::size_t) (from + i)];
    return std::sqrt (e / (double) len);
}

// The canceller's exact transfer at frequency f: H = 1 / (1 + Σ_k G_k),
// G_k(z) = μ·(z·cos ω_k − 1) / (z² − 2z·cos ω_k + 1), ω_k = 2π·k·fOsc/fs. Derived from the LMS loop
// (v = w·r obeys v[n+1] = e^{jω}(v[n] + μ·e[n]); the correction is Re v), independently of the code.
inline double exactGainDb (double f, double fOsc, int nActive, double bwHz, double fs)
{
    const double mu = 2.0 * kPi * bwHz / fs;
    const std::complex<double> z = std::polar (1.0, 2.0 * kPi * f / fs);
    std::complex<double> G = 0.0;
    for (int k = 1; k <= nActive; ++k)
    {
        const double c = std::cos (2.0 * kPi * k * fOsc / fs);
        G += mu * (z * c - 1.0) / (z * z - 2.0 * z * c + 1.0);
    }
    return dbOf (std::abs (1.0 / (1.0 + G)));
}

// The single-notch textbook shape: attenuation at offset d from the centre for a −3 dB width bw.
inline double singleNotchDb (double d, double bw) { return 10.0 * std::log10 (d * d / (d * d + 0.25 * bw * bw)); }

// Run a mono buffer through `h` in fixed blocks.
template <class H>
inline void runMono (H& h, std::vector<float>& x, int block = 512)
{
    for (std::size_t o = 0; o < x.size(); o += (std::size_t) block)
    {
        float* io[1] { x.data() + o };
        felitronics::test::run (h.process (io, 1, (int) std::min<std::size_t> ((std::size_t) block, x.size() - o)));
    }
}

inline felitronics::dehum::HumRemoverParams params (felitronics::dehum::Mains m, double bw)
{
    felitronics::dehum::HumRemoverParams p;
    p.mains = m;
    p.bandwidthHz = bw;
    return p;
}

} // namespace dehumtest

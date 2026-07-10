// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

//==============================================================================
// felitronics::measurement — the OFFLINE, MESSAGE-THREAD-ONLY double-precision math floor for IR
// measurement (allocates; NOT RT-safe by design). This is deliberately NOT the RT FFT: real-time
// consumers use the float `felitronics::core::fft` SEAM. Measurement deconvolves whole ~6 s captures
// (~2^20-point transforms) where the numerical noise floor must stay far below the analog chain's —
// a float 2^20 FFT's ~-120 dB roundoff, amplified ~+30 dB by the Farina inverse's +3 dB/oct HF boost,
// would land within ~20-40 dB of a real capture's tail and stop being transparent. Double keeps the
// round-trip floor near -208 dBFS at 2^20 (measured; the `w *= wlen` twiddle recurrence drifts ~O(len)·eps
// but that is still ~-178 dB after the Farina boost — far below a ~-120 dB analog floor, so the recurrence
// is kept for speed). PROMOTION TRIGGER: the day a SECOND offline-double-FFT consumer appears
// (analysis / lineareq), move `detail::fftInplace` + `convolve` to `felitronics/core/OfflineFft.h`
// (a mechanical header move). Correctness is oracle-anchored (sweep⊛inverse≈δ + known-answer tests).
//==============================================================================

#include <felitronics/core/Math.h>   // felitronics::core::kPi

#include <complex>
#include <cstddef>
#include <span>
#include <vector>

namespace felitronics::measurement
{

// Smallest power of two >= n (>= 1). Returns 0 if n exceeds the largest representable power of two
// (n > 2^63 on 64-bit) — the `p != 0` guard stops the shift from wrapping to 0 and looping forever.
inline std::size_t nextPow2 (std::size_t n) noexcept
{
    std::size_t p = 1;
    while (p < n && p != 0) p <<= 1;
    return p;
}

namespace detail
{
// In-place iterative radix-2 Cooley–Tukey. sign = -1 forward, +1 inverse (inverse NOT normalized —
// callers divide by n). a.size() must be a power of two. All twiddles/butterflies in double: this is
// the offline reference path, so accuracy beats speed.
inline void fftInplace (std::vector<std::complex<double>>& a, int sign) noexcept
{
    const std::size_t n = a.size();
    if (n < 2) return;
    for (std::size_t i = 1, j = 0; i < n; ++i)   // bit-reversal permutation
    {
        std::size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap (a[i], a[j]);
    }
    for (std::size_t len = 2; len <= n; len <<= 1)
    {
        const double ang = (double) sign * 2.0 * core::kPi / (double) len;
        const std::complex<double> wlen (std::cos (ang), std::sin (ang));
        for (std::size_t i = 0; i < n; i += len)
        {
            std::complex<double> w (1.0, 0.0);            // twiddle recurrence — floor ~-208 dB at 2^20 (see header)
            for (std::size_t k = 0; k < len / 2; ++k)
            {
                const std::complex<double> u = a[i + k];
                const std::complex<double> v = a[i + k + len / 2] * w;
                a[i + k]           = u + v;
                a[i + k + len / 2] = u - v;
                w *= wlen;
            }
        }
    }
}
} // namespace detail

// Real linear convolution via double FFT. Returns length x.size()+h.size()-1 (empty if either input
// is empty). OFFLINE — allocates. The measurement pipeline's one MAC (recording ⊛ inverse, sweep ⊛ inverse).
inline std::vector<double> convolve (std::span<const double> x, std::span<const double> h)
{
    if (x.empty() || h.empty()) return {};
    const std::size_t out = x.size() + h.size() - 1;
    const std::size_t n   = nextPow2 (out);
    if (n < 2) return {};                          // out exceeds the largest representable power of two
    std::vector<std::complex<double>> X (n, std::complex<double> {}), H (n, std::complex<double> {});
    for (std::size_t i = 0; i < x.size(); ++i) X[i] = x[i];
    for (std::size_t i = 0; i < h.size(); ++i) H[i] = h[i];
    detail::fftInplace (X, -1);
    detail::fftInplace (H, -1);
    for (std::size_t i = 0; i < n; ++i) X[i] *= H[i];
    detail::fftInplace (X, +1);
    const double invn = 1.0 / (double) n;
    std::vector<double> y (out);
    for (std::size_t i = 0; i < out; ++i) y[i] = X[i].real() * invn;
    return y;
}

// Linear magnitude spectrum |X[k]| of a real signal, zero-padded to nfft (pow2). Length nfft/2.
inline std::vector<double> magSpectrum (std::span<const double> x, std::size_t nfft)
{
    if (nfft < 2) return {};
    nfft = nextPow2 (nfft);                        // radix-2 requires a power of two — round up (heal, don't UB)
    if (nfft < 2) return {};                        // nextPow2 overflowed
    std::vector<std::complex<double>> X (nfft, std::complex<double> {});
    for (std::size_t i = 0; i < x.size() && i < nfft; ++i) X[i] = x[i];
    detail::fftInplace (X, -1);
    std::vector<double> m (nfft / 2);
    for (std::size_t i = 0; i < nfft / 2; ++i) m[i] = std::abs (X[i]);
    return m;
}

} // namespace felitronics::measurement

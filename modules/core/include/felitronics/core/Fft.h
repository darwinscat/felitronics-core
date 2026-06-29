// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/core/Math.h>

#include <cmath>
#include <complex>
#include <concepts>
#include <vector>

//==============================================================================
// felitronics::core::fft — the FFT SEAM (the architecture keystone). A real-input FFT reached through a
// COMPILE-TIME backend (template/concept, never a vtable in the convolver hot path). Each backend owns
// its OWN spectrum layout (so a SIMD backend — pffft/vDSP — never pays an O(N) repack to a forced
// lowest-common-denominator); a consumer (the partitioned convolver) only ever calls forward / inverse /
// spectralMultiplyAdd and never indexes spectrum bins. A plan is built in prepare() (all allocation
// there); forward()/inverse()/spectralMultiplyAdd() are RT-safe (no alloc/lock/throw).
//
// Ships a scalar radix-2 reference backend so the core self-tests JUCE-free with no heavy dep. Real
// backends (pffft/kissfft = BSD on desktop+wasm; juce::dsp::FFT in the JUCE adapter; CMSIS embedded)
// plug in later as compiled targets satisfying RealFftBackend — pffft's native zconvolve_accumulate maps
// straight onto spectralMultiplyAdd.
//==============================================================================
namespace felitronics::core::fft
{

inline constexpr bool isPow2 (int n) noexcept { return n > 0 && (n & (n - 1)) == 0; }

// A real-FFT backend for size N (pow2). `spectrumFloats(N)` is the backend's spectrum buffer length in
// floats (layout is the backend's business). `inverse` is 1/N-normalized so a round-trip is identity.
// `spectralMultiplyAdd(a,b,acc)` does acc += a (.*) b in the backend's own layout (the convolver's MAC).
template <class B>
concept RealFftBackend = requires (B b, const float* r, float* w, const float* s, float* acc, int n) {
    { B::spectrumFloats (n) } noexcept -> std::convertible_to<int>;
    { b.prepare (n) }         noexcept -> std::same_as<bool>;
    { b.forward (r, w) }      noexcept -> std::same_as<void>;   // real[N]      -> spectrum[spectrumFloats(N)]
    { b.inverse (s, w) }      noexcept -> std::same_as<void>;   // spectrum     -> real[N]  (1/N normalized)
    { b.spectralMultiplyAdd (s, s, acc) } noexcept -> std::same_as<void>;
};

//==============================================================================
// ScalarRadix2Real — the JUCE-free reference + spike default. Correctness-first (a full complex radix-2
// transform of the real input; a real backend is ~2× faster but this is the analytic reference the tests
// trust). Packed spectrum (N floats): s[0]=Re[0] (DC), s[1]=Re[N/2] (Nyquist), then s[2k]=Re[k],
// s[2k+1]=Im[k] for 1<=k<N/2.
class ScalarRadix2Real
{
public:
    static constexpr int spectrumFloats (int n) noexcept { return n; }

    bool prepare (int n) noexcept
    {
        if (! isPow2 (n) || n < 4) return false;
        n_ = n;
        scratch_.assign ((std::size_t) n, std::complex<float> {});   // ALLOC here (prepare) only
        return true;
    }

    int size() const noexcept { return n_; }

    void forward (const float* real, float* spec) noexcept
    {
        for (int i = 0; i < n_; ++i) scratch_[(std::size_t) i] = std::complex<float> (real[i], 0.0f);
        transform (scratch_.data(), n_, false);
        spec[0] = scratch_[0].real();                       // DC (real)
        spec[1] = scratch_[(std::size_t) (n_ / 2)].real();  // Nyquist (real)
        for (int k = 1; k < n_ / 2; ++k)
        {
            spec[2 * k]     = scratch_[(std::size_t) k].real();
            spec[2 * k + 1] = scratch_[(std::size_t) k].imag();
        }
    }

    void inverse (const float* spec, float* real) noexcept
    {
        scratch_[0]                            = std::complex<float> (spec[0], 0.0f);
        scratch_[(std::size_t) (n_ / 2)]       = std::complex<float> (spec[1], 0.0f);
        for (int k = 1; k < n_ / 2; ++k)       // Hermitian symmetry rebuilds the upper half
        {
            const std::complex<float> c (spec[2 * k], spec[2 * k + 1]);
            scratch_[(std::size_t) k]          = c;
            scratch_[(std::size_t) (n_ - k)]   = std::conj (c);
        }
        transform (scratch_.data(), n_, true);
        const float inv = 1.0f / (float) n_;
        for (int i = 0; i < n_; ++i) real[i] = scratch_[(std::size_t) i].real() * inv;
    }

    // acc[] += a[] (complex .*) b[] in the packed layout. DC + Nyquist are real-only.
    void spectralMultiplyAdd (const float* a, const float* b, float* acc) const noexcept
    {
        acc[0] += a[0] * b[0];
        acc[1] += a[1] * b[1];
        for (int k = 1; k < n_ / 2; ++k)
        {
            const float ar = a[2 * k], ai = a[2 * k + 1], br = b[2 * k], bi = b[2 * k + 1];
            acc[2 * k]     += ar * br - ai * bi;
            acc[2 * k + 1] += ar * bi + ai * br;
        }
    }

private:
    // Iterative Cooley-Tukey radix-2, in place. Twiddles accumulated in double for reference accuracy;
    // butterflies in float (the hot-path type). Caller normalizes the inverse.
    static void transform (std::complex<float>* a, int n, bool inverse) noexcept
    {
        for (int i = 1, j = 0; i < n; ++i)              // bit-reversal permutation
        {
            int bit = n >> 1;
            for (; j & bit; bit >>= 1) j ^= bit;
            j ^= bit;
            if (i < j) std::swap (a[i], a[j]);
        }
        for (int len = 2; len <= n; len <<= 1)
        {
            const double ang = 2.0 * kPi / len * (inverse ? 1.0 : -1.0);
            const std::complex<double> wlen (std::cos (ang), std::sin (ang));
            for (int i = 0; i < n; i += len)
            {
                std::complex<double> w (1.0, 0.0);
                for (int k = 0; k < len / 2; ++k)
                {
                    const std::complex<float> u = a[i + k];
                    const std::complex<float> v = a[i + k + len / 2]
                                                * std::complex<float> ((float) w.real(), (float) w.imag());
                    a[i + k]            = u + v;
                    a[i + k + len / 2]  = u - v;
                    w *= wlen;
                }
            }
        }
    }

    int n_ = 0;
    std::vector<std::complex<float>> scratch_;
};

static_assert (RealFftBackend<ScalarRadix2Real>, "ScalarRadix2Real must satisfy the seam");

// The default real-FFT backend (swap per tier later via the template seam).
using DefaultRealFft = ScalarRadix2Real;

} // namespace felitronics::core::fft

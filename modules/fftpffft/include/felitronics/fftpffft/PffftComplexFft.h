// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// felitronics::fftpffft::PffftComplexFft — pffft's COMPLEX transform, in canonical order: out[k] is bin k,
// real then imaginary, unnormalised. It is not a core::fft::RealFftBackend and is not meant to be: that seam is
// the convolvers' (a real transform of a power-of-two size with a spectral multiply-add). This is the other
// thing pffft does — a complex transform of any N = 2^a 3^b 5^c with a >= 4 — for a consumer that needs exactly
// that. felitronics::codecgrid is the first: an MDCT of 2M samples is one complex FFT of M/2 points, and CELT's
// frame makes that 480.
//
// Same shape as its siblings: the vendored <pffft.h> is private to the .cpp; the plan and the aligned scratch are
// built in prepare(); forward() is memcpy + pffft only (no C++ allocation, no lock, no throw). Input and output
// may be unaligned — pffft needs SIMD alignment on its buffers, so both bounce through an aligned vector. The
// bounce is why forward() is const over mutable scratch: one object serves one thread at a time, as every
// backend here does.

#pragma once

#include <felitronics/core/Fft.h>

#include <complex>

namespace felitronics::fftpffft
{

class PffftComplexFft
{
public:
    PffftComplexFft() noexcept = default;
    ~PffftComplexFft();

    PffftComplexFft (const PffftComplexFft&)            = delete;   // owns the pffft plan + aligned scratch
    PffftComplexFft& operator= (const PffftComplexFft&) = delete;
    PffftComplexFft (PffftComplexFft&&)                 = delete;
    PffftComplexFft& operator= (PffftComplexFft&&)      = delete;

    static constexpr int kMaxSize = 1 << 20;

    // N = 2^a 3^b 5^c with a >= 4 (pffft's complex transform works in blocks of 16), up to kMaxSize.
    static constexpr bool supported (int n) noexcept
    {
        if (n < 16 || n > kMaxSize || n % 16 != 0) return false;
        for (int p : { 2, 3, 5 }) while (n % p == 0) n /= p;
        return n == 1;
    }

    static int simdWidth() noexcept;                                    // pffft_simd_size(): 4 = SSE/NEON kernel, 1 = scalar fallback

    [[nodiscard]] bool prepare (int n);                                 // message thread — plan + scratch
    int size() const noexcept { return n_; }

    // out[k] = sum_n in[n] exp (-2 pi i k n / N). `in` and `out` are distinct arrays of N.
    void forward (const std::complex<float>* in, std::complex<float>* out) const noexcept;

private:
    void release() noexcept;

    void* setup_ = nullptr;                                             // opaque pffft plan
    mutable core::fft::AlignedVector<float> in_, out_, work_;           // aligned bounces + pffft work (2N floats each)
    int n_ = 0;
};

} // namespace felitronics::fftpffft

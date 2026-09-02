// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// felitronics::fftpffft::PffftOrderedRealFft — the pffft SIMD real transform in CANONICAL order: the
// packed-Hermitian layout the scalar reference writes (s[0]=DC, s[1]=Nyquist, s[2k]/s[2k+1]=Re/Im of bin k).
// pffft's rfftf assembles exactly that when asked for an ordered transform (F(0) and F(N/2), both real,
// share the first complex slot), so this backend advertises kPackedHermitianSpectrum and is admissible
// wherever bins are READ: the spectrum panes (SpectrumPaneT / MultiResSpectrumPaneT), and — by the same
// gate — the design-time FFTs, should anyone choose to run them on it.
//
// It sits BESIDE PffftRealFft, not instead of it. The z-order backend is the audio-path backend: the
// convolvers never index a bin and pffft's zconvolve_accumulate vectorises their MAC, which the ordered
// layout cannot. Ordering costs a reorder pass per transform (pffft_transform_ordered = transform +
// zreorder), so this one is for consumers that would otherwise pay the scalar radix-2 for every bin
// they look at — an analyzer running 16384 + 4096 + 3×1024 points thirty times a second.
//
// Same shape as PffftRealFft: the vendored <pffft.h> is private to the .cpp; the plan and the aligned
// scratch are built in prepare(); forward()/inverse() are memcpy + pffft only (no alloc, no lock, no
// throw). Both the real input and the packed spectrum may be unaligned — pffft needs SIMD alignment on
// its buffers, so each side bounces through an aligned vector (N floats, negligible next to the FFT).
// spectralMultiplyAdd is the scalar packed complex MAC (there is nothing to vectorise in the ordered
// layout without knowing the SIMD stride) — correct, and not this backend's job.

#pragma once

#include <felitronics/core/Fft.h>

namespace felitronics::fftpffft
{

class PffftOrderedRealFft
{
public:
    PffftOrderedRealFft() noexcept = default;
    ~PffftOrderedRealFft();

    PffftOrderedRealFft (const PffftOrderedRealFft&)            = delete;   // owns the pffft plan + aligned scratch
    PffftOrderedRealFft& operator= (const PffftOrderedRealFft&) = delete;
    PffftOrderedRealFft (PffftOrderedRealFft&&)                 = delete;
    PffftOrderedRealFft& operator= (PffftOrderedRealFft&&)      = delete;

    static constexpr bool kPackedHermitianSpectrum = true;              // s[0]=DC, s[1]=Nyq, s[2k]/s[2k+1]=Re/Im
    static constexpr int  spectrumFloats (int n) noexcept { return n; }
    static int simdWidth() noexcept;                                    // pffft_simd_size(): 4 = SSE/NEON kernel, 1 = scalar fallback

    bool prepare (int n) noexcept;                                      // message thread — plan + scratch; pow2, 32 ≤ N ≤ 2^26
    void forward (const float* real, float* spec) noexcept;            // real[N] -> packed spectrum[N]
    void inverse (const float* spec, float* real) noexcept;            // packed spectrum[N] -> real[N], 1/N normalized
    void spectralMultiplyAdd (const float* a, const float* b, float* acc) const noexcept;   // acc += a (.*) b, packed layout

    int size() const noexcept { return n_; }

private:
    void release() noexcept;

    void* setup_ = nullptr;                                             // opaque pffft plan
    core::fft::AlignedVector<float> in_, out_, work_;                   // aligned bounces + pffft work (N floats)
    int n_ = 0;
};

static_assert (core::fft::RealFftBackend<PffftOrderedRealFft>,
               "PffftOrderedRealFft must satisfy the FFT seam (RealFftBackend)");
static_assert (core::fft::PackedHermitianSpectrum<PffftOrderedRealFft>,
               "PffftOrderedRealFft writes the packed-Hermitian layout — that is its reason to exist");

} // namespace felitronics::fftpffft

// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// felitronics::fftpffft::PffftRealFft — a RealFftBackend backed by pffft ("Pretty Fast FFT"): the UNORDERED
// (z-order) real transform + pffft_zconvolve_accumulate for a fully VECTORIZED spectral multiply-accumulate.
// This is the SIMD AUDIO-path backend for the partitioned convolvers (MatrixConvolver / ConvolutionEngine /
// PartitionedConvolver). It vectorizes BOTH the transform (~17-24x the scalar reference on long transforms)
// AND the per-partition complex MAC that dominates a long IR's tail — the point of the whole exercise.
//
// LAYOUT — pffft's private z-order spectrum; spectrumFloats(N)=N. The convolver only forwards / inverts /
// MACs and never indexes a bin, so the layout stays fully opaque. It is NOT the packed-Hermitian layout the
// lineareq FIR designers hand-write, so this backend deliberately does NOT advertise kPackedHermitianSpectrum
// — and PR2's PackedHermitianSpectrum concept REJECTS it from any design-time FFT at COMPILE time (the
// static_assert below is the runtime-safety proof). Use it ONLY as the audio backend, e.g.
// MatrixConvolver<PffftRealFft> or BasicLinearPhaseEq<PffftRealFft> (whose design FFTs stay pinned scalar).
//
// The pffft C library is a PRIVATE implementation detail — included only in PffftRealFft.cpp. A consumer sees
// just this header and links felitronics::fftpffft; the vendored <pffft.h> never reaches its include path (no
// pollution, no ODR clash with a consumer's own pffft). Hence the decl-only shape: the hot-path methods are
// per-BLOCK (not per-sample), so being out-of-line costs nothing measurable.
//
// SEAM CONTRACT — the z-order spectrum is a pure function of (N, SIMD width), independent of the instance, so
// the two backend instances a convolver holds (message-thread build + audio) are interchangeable, as the seam
// requires. RT-safety — the pffft plan + aligned scratch are built in prepare(); the hot path is
// memcpy / pffft_transform / pffft_zconvolve_accumulate only (no alloc: 'work' is pre-sized so pffft never
// falls back to alloca).
#pragma once

#include <felitronics/core/Fft.h>

namespace felitronics::fftpffft
{

class PffftRealFft
{
public:
    PffftRealFft() noexcept = default;
    ~PffftRealFft();

    PffftRealFft (const PffftRealFft&)            = delete;   // owns the pffft plan + aligned scratch
    PffftRealFft& operator= (const PffftRealFft&) = delete;
    PffftRealFft (PffftRealFft&&)                 = delete;
    PffftRealFft& operator= (PffftRealFft&&)      = delete;

    static constexpr int spectrumFloats (int n) noexcept { return n; }   // pffft z-order: N floats
    static int simdWidth() noexcept;   // pffft_simd_size(): 4 = SSE/NEON kernel compiled in, 1 = scalar fallback (no speedup)

    bool prepare (int n) noexcept;                                        // message thread — builds plan + scratch
    void forward (const float* real, float* spec) noexcept;              // real[N]  -> z-order spectrum[N]
    void inverse (const float* spec, float* real) noexcept;             // z-order spectrum[N] -> real[N], 1/N normalized
    void spectralMultiplyAdd (const float* a, const float* b, float* acc) const noexcept;   // acc += a (.*) b, VECTORIZED

    int size() const noexcept { return n_; }

private:
    void release() noexcept;                                             // destroy the pffft plan (scratch is RAII)

    void* setup_ = nullptr;   // opaque pffft plan (a PFFFT_Setup*) — the pffft type stays entirely private to the .cpp
    core::fft::AlignedVector<float> in_, out_, work_;   // aligned real-side bounce (in/out) + pffft work buffer (N floats)
    int n_ = 0;
};

// The whole reason this module exists behind the seam — and the runtime-safety proof for PR2's design/audio
// split: pffft satisfies the backend concept, and its z-order layout is (correctly) NOT packed-Hermitian, so
// instantiating a design-time FFT / MixedPhaseFir on it fails to COMPILE rather than silently mis-designing.
static_assert (core::fft::RealFftBackend<PffftRealFft>,
               "PffftRealFft must satisfy the FFT seam (RealFftBackend)");
static_assert (! core::fft::PackedHermitianSpectrum<PffftRealFft>,
               "pffft's z-order spectrum is NOT the hand-packed layout — it must be rejected by the design-FFT gate (C1)");

} // namespace felitronics::fftpffft

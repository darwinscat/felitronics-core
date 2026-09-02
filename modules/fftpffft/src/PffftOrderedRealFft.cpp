// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// PffftOrderedRealFft out-of-line bodies — with PffftRealFft.cpp, the only translation units that see the
// vendored pffft C API (see the header for why this backend exists beside the z-order one).

#include <felitronics/fftpffft/PffftOrderedRealFft.h>

#include <pffft.h>

#include <cstddef>
#include <cstring>

namespace felitronics::fftpffft
{

namespace { inline PFFFT_Setup* asPlan (void* p) noexcept { return static_cast<PFFFT_Setup*> (p); } }

PffftOrderedRealFft::~PffftOrderedRealFft() { release(); }

bool PffftOrderedRealFft::prepare (int n) noexcept
{
    release();
    // The same admissible set and the same primary guard as PffftRealFft: pffft_new_setup assert()s on a
    // bad N rather than returning NULL for every one of them, so the range check must come first.
    if (! core::fft::isPow2 (n) || n < 32 || n > (1 << 26)) return false;

    setup_ = pffft_new_setup (n, PFFFT_REAL);
    if (setup_ == nullptr) return false;

    const auto len = static_cast<std::size_t> (n);
    in_.assign   (len, 0.0f);   // aligned bounce for the (possibly unaligned) real input / packed spectrum
    out_.assign  (len, 0.0f);   // aligned transform output, copied to the caller's (possibly unaligned) buffer
    work_.assign (len, 0.0f);   // N floats for a REAL transform — non-NULL so pffft never falls back to alloca
    n_ = n;
    return true;
}

// real[N] -> packed [DC, Nyq, re1, im1, …]. pffft's ordered real forward is FFTPACK's rfftf with the
// canonical interleaving, F(0) and F(N/2) sharing slot 0/1 — the scalar reference's layout exactly.
void PffftOrderedRealFft::forward (const float* real, float* spec) noexcept
{
    if (n_ <= 0) return;
    const auto bytes = static_cast<std::size_t> (n_) * sizeof (float);
    std::memcpy (in_.data(), real, bytes);
    pffft_transform_ordered (asPlan (setup_), in_.data(), out_.data(), work_.data(), PFFFT_FORWARD);
    std::memcpy (spec, out_.data(), bytes);
}

// packed spectrum[N] -> real[N], 1/N normalized (pffft's backward is unscaled).
void PffftOrderedRealFft::inverse (const float* spec, float* real) noexcept
{
    if (n_ <= 0) return;
    std::memcpy (in_.data(), spec, static_cast<std::size_t> (n_) * sizeof (float));
    pffft_transform_ordered (asPlan (setup_), in_.data(), out_.data(), work_.data(), PFFFT_BACKWARD);
    const float invN = 1.0f / static_cast<float> (n_);
    for (int i = 0; i < n_; ++i) real[i] = out_[(std::size_t) i] * invN;
}

// acc += a (.*) b in the packed layout — the scalar reference's loop. DC and Nyquist are real-only.
void PffftOrderedRealFft::spectralMultiplyAdd (const float* a, const float* b, float* acc) const noexcept
{
    if (n_ <= 0) return;
    acc[0] += a[0] * b[0];
    acc[1] += a[1] * b[1];
    for (int k = 1; k < n_ / 2; ++k)
    {
        const float ar = a[2 * k], ai = a[2 * k + 1], br = b[2 * k], bi = b[2 * k + 1];
        acc[2 * k]     += ar * br - ai * bi;
        acc[2 * k + 1] += ar * bi + ai * br;
    }
}

int PffftOrderedRealFft::simdWidth() noexcept { return pffft_simd_size(); }

void PffftOrderedRealFft::release() noexcept
{
    if (setup_ != nullptr) { pffft_destroy_setup (asPlan (setup_)); setup_ = nullptr; }
    n_ = 0;
}

} // namespace felitronics::fftpffft

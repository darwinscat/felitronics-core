// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// PffftComplexFft out-of-line bodies — with its two siblings, the only translation units that see the vendored
// pffft C API.

#include <felitronics/fftpffft/PffftComplexFft.h>

#include <pffft.h>

#include <cstddef>
#include <cstring>

namespace felitronics::fftpffft
{

namespace { inline PFFFT_Setup* asPlan (void* p) noexcept { return static_cast<PFFFT_Setup*> (p); } }

PffftComplexFft::~PffftComplexFft() { release(); }

bool PffftComplexFft::prepare (int n)
{
    release();
    // pffft_new_setup assert()s on a size it cannot factor rather than returning NULL, so the check comes first.
    if (! supported (n)) return false;

    setup_ = pffft_new_setup (n, PFFFT_COMPLEX);
    if (setup_ == nullptr) return false;

    const auto len = 2u * static_cast<std::size_t> (n);
    in_.assign   (len, 0.0f);
    out_.assign  (len, 0.0f);
    work_.assign (len, 0.0f);   // 2N floats for a COMPLEX transform — non-NULL, so pffft never takes its stack fallback
    n_ = n;
    return true;
}

// std::complex<float> is two floats, real then imaginary — the interleaving pffft's ordered complex transform
// reads and writes.
void PffftComplexFft::forward (const std::complex<float>* in, std::complex<float>* out) const noexcept
{
    if (n_ <= 0) return;
    const auto bytes = 2u * static_cast<std::size_t> (n_) * sizeof (float);
    std::memcpy (in_.data(), in, bytes);
    pffft_transform_ordered (asPlan (setup_), in_.data(), out_.data(), work_.data(), PFFFT_FORWARD);
    std::memcpy (out, out_.data(), bytes);
}

int PffftComplexFft::simdWidth() noexcept { return pffft_simd_size(); }

void PffftComplexFft::release() noexcept
{
    if (setup_ != nullptr) { pffft_destroy_setup (asPlan (setup_)); setup_ = nullptr; }
    n_ = 0;
}

} // namespace felitronics::fftpffft

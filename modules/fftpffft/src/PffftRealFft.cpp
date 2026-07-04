// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// PffftRealFft out-of-line bodies — the ONLY translation unit that sees the vendored pffft C API, keeping
// <pffft.h> private to felitronics::fftpffft (see the header for the design rationale).

#include <felitronics/fftpffft/PffftRealFft.h>

#include <pffft.h>

#include <cstddef>
#include <cstring>

namespace felitronics::fftpffft
{

namespace { inline PFFFT_Setup* asPlan (void* p) noexcept { return static_cast<PFFFT_Setup*> (p); } }

PffftRealFft::~PffftRealFft() { release(); }

bool PffftRealFft::prepare (int n) noexcept
{
    release();

    // pffft's real transform admits N = (2^a)(3^b)(5^c) with a>=5; we only ever use pow2 sizes, so
    // "pow2 and 32 <= N <= 2^26" is exactly the admissible set. This pre-guard is the PRIMARY defense:
    // the vendored pffft_new_setup assert()s on an out-of-range N (it does NOT cleanly return NULL for
    // every bad input), so an unguarded call would trap in Debug / be UB under NDEBUG. Keep it ahead of
    // the setup call; the NULL-check below then only covers a genuine factorization/allocation failure.
    if (! core::fft::isPow2 (n) || n < 32 || n > (1 << 26)) return false;

    setup_ = pffft_new_setup (n, PFFFT_REAL);
    if (setup_ == nullptr) return false;

    // 64-byte AlignedVector (a superset of pffft's 16-byte need) — RAII scratch, no manual aligned malloc/
    // free. assign() throwing on OOM inside this noexcept prepare() terminates, matching ScalarRadix2Real.
    const auto len = static_cast<std::size_t> (n);
    in_.assign   (len, 0.0f);   // aligned bounce for a (possibly unaligned) real input
    out_.assign  (len, 0.0f);   // aligned scratch for the real IFFT output
    work_.assign (len, 0.0f);   // N floats (REAL transform; a COMPLEX transform would need 2N) — non-NULL so
                                //   pffft_transform never falls back to alloca (RT-safety)
    n_ = n;
    return true;
}

// real[N] -> unordered spectrum[N]. `real` is a (possibly unaligned) engine time-domain slice -> bounce
// through the aligned in_; `spec` is an aligned engine spectrum row -> pffft writes straight into it.
void PffftRealFft::forward (const float* real, float* spec) noexcept
{
    if (n_ <= 0) return;
    std::memcpy (in_.data(), real, static_cast<std::size_t> (n_) * sizeof (float));
    pffft_transform (asPlan (setup_), in_.data(), spec, work_.data(), PFFFT_FORWARD);
}

// unordered spectrum[N] -> real[N], 1/N normalized (pffft's backward is unscaled: BACKWARD(FORWARD(x))=N*x).
// `spec` is aligned; transform into the aligned out_ scratch, then scale into the caller's `real` (which may
// be unaligned and need not be touched by pffft) — simpler and, measured, marginally faster than in-place.
void PffftRealFft::inverse (const float* spec, float* real) noexcept
{
    if (n_ <= 0) return;
    pffft_transform (asPlan (setup_), spec, out_.data(), work_.data(), PFFFT_BACKWARD);
    const float invN = 1.0f / static_cast<float> (n_);
    for (int i = 0; i < n_; ++i) real[i] = out_[(std::size_t) i] * invN;
}

// acc += a (.*) b in the z-order layout, VECTORIZED. a/b/acc are aligned engine spectrum rows (the convolver
// holds them in AlignedVector), which pffft_zconvolve_accumulate requires; scaling 1.0 => pure accumulate.
void PffftRealFft::spectralMultiplyAdd (const float* a, const float* b, float* acc) const noexcept
{
    if (n_ <= 0) return;
    pffft_zconvolve_accumulate (asPlan (setup_), a, b, acc, 1.0f);
}

// pffft's compiled SIMD width: 4 on SSE/NEON, 1 if pffft.c fell back to scalar. A test asserting ==4 catches a
// silent scalar build that would pass every null yet deliver none of the SIMD speedup this module exists for.
int PffftRealFft::simdWidth() noexcept { return pffft_simd_size(); }

void PffftRealFft::release() noexcept
{
    if (setup_ != nullptr) { pffft_destroy_setup (asPlan (setup_)); setup_ = nullptr; }
    n_ = 0;   // in_/out_/work_ are AlignedVector — freed by their destructors; a re-prepare reassigns them
}

} // namespace felitronics::fftpffft

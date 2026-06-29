// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa. Part of felitronics-core — see LICENSE.

#pragma once

#include <cmath>

// Law 8 (denormals). Every feedback kernel flushes denormal state in SOFTWARE — `flushDenormal()`
// below — which works on EVERY tier (no FP-control register needed), so `wasm-audio` (and many
// embedded ARMs that expose no MXCSR/FPCR) stay safe. Hardware FTZ/DAZ (`ScopedFlushToZero`) is a
// DESKTOP optimization an adapter MAY add — never a correctness crutch, never set globally by the
// core. Do NOT compile with -ffast-math (it breaks the NaN/inf semantics the tests assert).

#if defined(__SSE2__) || defined(_M_X64) || defined(_M_AMD64) || defined(__x86_64__) || defined(__i386__) || defined(_M_IX86)
  #include <xmmintrin.h>
  #include <pmmintrin.h>
  #define FELITRONICS_FTZ_X86 1
#elif defined(__aarch64__)
  #define FELITRONICS_FTZ_AARCH64 1
#endif

namespace felitronics::core
{

// Per-block (or per-sample) software denormal flush: zap tiny state to exact zero so decaying tails
// don't sustain subnormals (the 10–100× CPU spike). 1e-15 is far above the subnormal range
// (~1.2e-38 float / ~2.2e-308 double) yet inaudible (< ~-300 dB); a block is too short to re-traverse
// the gap, so state never reaches subnormal across blocks.
inline void flushDenormal (float&  x) noexcept { if (std::fabs (x) < 1e-15f) x = 0.0f; }
inline void flushDenormal (double& x) noexcept { if (std::fabs (x) < 1e-15)  x = 0.0;  }

//==============================================================================
// ScopedFlushToZero — OPTIONAL desktop hardware FTZ/DAZ, RAII (set on construct, restore on destroy).
// A belt-and-suspenders bonus over the software flush, NOT a correctness requirement. No-op where the
// FP-control register is unavailable (wasm-audio / many embedded) — which is exactly why the core
// never relies on it for correctness.
class ScopedFlushToZero
{
public:
    ScopedFlushToZero (const ScopedFlushToZero&) = delete;
    ScopedFlushToZero& operator= (const ScopedFlushToZero&) = delete;

#if defined(FELITRONICS_FTZ_X86)
    ScopedFlushToZero() noexcept : saved_ (_mm_getcsr()) { _mm_setcsr (saved_ | 0x8040u); }   // FTZ(0x8000) | DAZ(0x0040)
    ~ScopedFlushToZero() noexcept { _mm_setcsr (saved_); }
private:
    unsigned int saved_;
#elif defined(FELITRONICS_FTZ_AARCH64)
    ScopedFlushToZero() noexcept
    {
        __asm__ volatile ("mrs %0, fpcr" : "=r" (saved_));
        const unsigned long long v = saved_ | (1ull << 24);                                   // FPCR.FZ (bit 24)
        __asm__ volatile ("msr fpcr, %0" : : "r" (v) : "memory");
    }
    ~ScopedFlushToZero() noexcept { __asm__ volatile ("msr fpcr, %0" : : "r" (saved_) : "memory"); }
private:
    unsigned long long saved_ {};
#else
    ScopedFlushToZero() noexcept = default;   // no FP-control register (wasm-audio / embedded) → no-op
#endif
};

} // namespace felitronics::core

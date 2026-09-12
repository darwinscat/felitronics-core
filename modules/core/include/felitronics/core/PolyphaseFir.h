// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <cstddef>

//==============================================================================
// felitronics::core::firDot — THE polyphase-FIR inner product. One function, one summation order,
// four hand-written kernels (scalar / SSE2 / NEON / wasm-SIMD128) that all compute the SAME BITS.
//
// It replaced three separately written scalar loops — `PolyphaseOversampler::upsample`, the same
// class's `downsample`, and `TruePeakMeter::process` — which between them were ~83 % of a mastering
// render (measured: `downsample` 50 %, `upsample` 33 % of a 10 s sample of `fcore_master render`).
// The three were NOT the same loop: two gathered their coefficients with a stride of `L`
// (`proto[k*L+p]`) and one read them contiguously, and all three walked a modulo ring backwards with
// a wrap test inside the inner loop. What makes them one loop is the REPACKING at the call sites,
// not a cleverer kernel:
//   * coefficients are stored PHASE-MAJOR (`protoPhase[p][k]` contiguous in k), so the gather goes;
//   * sample history lives in a DOUBLE-LENGTH ring written NEWEST-FIRST — every sample stored twice,
//     at `pos` and at `pos + len` — so the window [pos, pos+len) is always contiguous and already in
//     the order the coefficients want, and the wrap test leaves the inner loop;
//   * `len` is padded to a multiple of four with +0.0f coefficients in `prepare()`, so there is no
//     tail and not one branch inside the loop. (+0.0f is exact; `acc + (0.0f * x)` returns `acc`
//     unchanged for every finite `x`, so padding is a numeric no-op — see THE ONE CAVEAT below.)
//
// ------------------------------------------------------------------------------------------------
// THE ORDER, which is the whole point. Every kernel computes exactly this and nothing else:
//
//     s0 s1 s2 s3 = +0.0f
//     for i in 0, 4, 8, ... len-4:
//         for j in 0..3:  p = a[i+j] * b[i+j]     (one rounding)
//                         sj = sj + p             (a SECOND rounding — never fused into the first)
//     return (s0 + s1) + (s2 + s3)
//
// FOUR accumulators, on every row, including the scalar one. That is a 128-bit vector's shape, and
// the consequence is accepted deliberately: A 256-BIT KERNEL WOULD BUY NOTHING HERE, because eight
// accumulators are a different summation and therefore different bits. One answer everywhere is
// worth more than the top of an AVX row's throughput, so there is no 256-bit path and adding one
// would break the only property this file exists to provide.
//
// ------------------------------------------------------------------------------------------------
// WHY IT OVERRIDES LAW 10 LOCALLY. The tree states `-ffp-contract=on` (DSP-ARCHITECTURE.md §2 law 10)
// because contraction FIXES two cancellation claims elsewhere. Here it does the opposite: `acc += a*b`
// is precisely the contractible form, arm64 has an `fmadd` and fuses it, and baseline x86-64 has no
// FMA instruction and cannot — so these three loops have been quietly computing different numbers on
// different rows all along. A fused multiply-add is a BETTER number and a DIFFERENT one; this
// primitive wants the same one everywhere, so it turns contraction off for itself.
//
// A header cannot do that with a compiler flag — a consumer compiles these lines with its own flags —
// so it takes a pragma, and the three compilers need three different ones. All six statements below
// were MEASURED on the actual toolchains by compiling with contraction forced ON and grepping the
// emitted asm for a fused op, not read out of a manual:
//
//   gcc 14.2 (x86-64) and 14.4 (aarch64):
//       `#pragma STDC FP_CONTRACT OFF`                       IGNORED in C++ — still emitted fmadd.
//       separate statements (`p = a*b;` then `acc = acc+p;`) STILL FUSED at -ffp-contract=fast,
//                                                            which is gcc's own default in gnu++ modes.
//       `#pragma GCC optimize("fp-contract=off")`             WORKS, and keeps working under -ffast-math.
//   clang 21 (Apple; emscripten is the same front end):
//       `#pragma clang fp contract(off)`                      WORKS at -ffp-contract=on/off (the flag
//                                                            this tree states, and clang's default).
//       `#pragma float_control(precise, on, push)`            ⚠ TURNS CONTRACTION BACK ON. It fused even
//                                                            under -ffp-contract=off. Do NOT add it here
//                                                            "for safety" — it is the opposite of safe.
//       at -ffp-contract=fast / -ffast-math                   NO pragma helps: the flag sets the backend's
//                                                            AllowFPOpFusion::Fast, which fuses below the
//                                                            level any source-level pragma reaches.
//   MSVC 19.44 (x64, /arch:AVX2 /fp:fast):
//       `#pragma fp_contract(off)`                            WORKS; and `float_control(push)` / `(pop)`
//                                                            does save and restore it, so the pragma
//                                                            here cannot leak into the including TU.
//
// So the defence is complete on gcc and MSVC and complete on clang for every flag this tree or its
// consumers actually use — but a clang build with `-ffast-math` defeats it, and nothing in a header
// can stop that. That is why the claim is GATED rather than asserted: `CoreTests.cpp` feeds `firDot`
// two operands whose fused and unfused results differ by one ulp and checks the unfused bits, so a
// row that lost the pragma goes red instead of going quietly wrong.
//
// ------------------------------------------------------------------------------------------------
// WHAT "THE SAME BITS" DOES AND DOES NOT COVER. The order and the contraction are the two things this
// file controls. Three more differences are properties of the RUNTIME, not of the code, and they are
// named here rather than pretended away — the crew round that preceded this file produced all three
// as working counterexamples:
//
//  1. DENORMAL FLUSHING (FTZ/DAZ). x86 flushes via MXCSR, aarch64 via FPCR.FZ, and WASM CANNOT FLUSH
//     AT ALL — it is strictly IEEE and has no control register. So a host (or this tree's own
//     `core::ScopedFlushToZero`) that turns flushing on makes the native rows disagree with wasm. And
//     this does NOT need denormal INPUTS: with a = 2^-126·(1+2^-23) and b = 0.5 every product is
//     2^-127 + 2^-150, an inexact subnormal; four of them sum to 2^-125, a perfectly NORMAL result
//     that reads 0x00000000 with flushing on and 0x01000000 with it off, from four normal operands.
//     The bit-identity gate therefore runs on ordinary audio magnitudes, where no intermediate can go
//     subnormal, and that is the honest scope of the claim.
//  2. ROUNDING MODE. MXCSR.RC / FPCR.RMode change every add and multiply below; wasm has no rounding
//     mode and is always nearest-even. `1.0f + 2^-24` is an exact tie: nearest-even gives 1.0f,
//     round-toward-+inf gives 1.0f+2^-23.
//  3. NaN SIGN AND PAYLOAD. x86 `MULPS(inf, 0)` yields the negative indefinite 0xffc00000; NEON `FMUL`
//     yields 0x7fc00000; and the wasm spec leaves NaN propagation nondeterministic on purpose. The
//     claim is for FINITE inputs.
//
// None of the three is reachable by writing the loop differently, and all three are outside what this
// primitive is for. Inside its scope — finite ordinary-magnitude audio, default rounding, the same
// denormal mode — five rows return identical bits.
//
// ------------------------------------------------------------------------------------------------
// THE ONE CAVEAT the zero padding carries. `0.0f * x` is +0.0f for every finite x and adding it is
// exact, so padding changes no finite answer. It is not free for a NON-finite one: the padded taps
// read history slots that have already left the true window, so a NaN or Inf that entered the stream
// poisons the output for up to three samples longer than it used to. Both shipped configurations
// (`tapsPerPhase` 64 and the meter's 12) are already multiples of four and pad by nothing at all, so
// this is reachable only with a hand-picked odd tap count AND a non-finite sample.
//==============================================================================
// ISA selection. Exactly one kernel is compiled; they are all the same arithmetic in the same order,
// so which one a row picks is a throughput decision and never a numeric one.
#if defined(__wasm_simd128__)
  #define FELITRONICS_FIR_WASM 1
#elif defined(__SSE2__) || defined(_M_X64) || defined(_M_AMD64) || (defined(_M_IX86_FP) && _M_IX86_FP >= 2)
  #define FELITRONICS_FIR_SSE2 1
#elif defined(__aarch64__) || defined(_M_ARM64) || defined(__ARM_NEON) || defined(__ARM_NEON__)
  #define FELITRONICS_FIR_NEON 1
#endif

#if defined(FELITRONICS_FIR_WASM)
  #include <wasm_simd128.h>
#elif defined(FELITRONICS_FIR_SSE2)
  #include <emmintrin.h>
#elif defined(FELITRONICS_FIR_NEON)
  #include <arm_neon.h>
#endif

//==============================================================================
// The contraction pragmas, per the measurements in the header comment. `__clang__` is tested FIRST
// because clang also defines `__GNUC__` and would otherwise take gcc's branch, where
// `#pragma GCC optimize` is silently accepted and does nothing.
#if defined(__clang__)
  #define FELITRONICS_FIR_EXACT_BEGIN
  #define FELITRONICS_FIR_EXACT_END
  #define FELITRONICS_FIR_NO_FMA        _Pragma ("clang fp contract(off)")
#elif defined(__GNUC__)
  #define FELITRONICS_FIR_EXACT_BEGIN   _Pragma ("GCC push_options") \
                                        _Pragma ("GCC optimize(\"fp-contract=off\")")
  #define FELITRONICS_FIR_EXACT_END     _Pragma ("GCC pop_options")
  #define FELITRONICS_FIR_NO_FMA
#elif defined(_MSC_VER)
  #define FELITRONICS_FIR_EXACT_BEGIN   __pragma (float_control (push)) __pragma (fp_contract (off))
  #define FELITRONICS_FIR_EXACT_END     __pragma (float_control (pop))
  #define FELITRONICS_FIR_NO_FMA
#else
  #define FELITRONICS_FIR_EXACT_BEGIN
  #define FELITRONICS_FIR_EXACT_END
  #define FELITRONICS_FIR_NO_FMA
#endif

FELITRONICS_FIR_EXACT_BEGIN

namespace felitronics::core
{

// Round a tap count up to the multiple of four `firDot` requires. Call it in `prepare()`; the
// coefficient array is zero-filled to this length and the history ring is sized from it.
constexpr int firPadLen (int len) noexcept { return (len + 3) & ~3; }

//==============================================================================
// THE inner product. `a` and `b` are `len` contiguous floats each; `len` must be > 0 and a MULTIPLE
// OF FOUR (the caller's `prepare()` guarantees it — see `firPadLen`). Loads may be unaligned.
//
// Returns the order stated at the top of this file, on every row. There is no second implementation
// of this sum in the tree: the naive sequential loop survives only as the ORACLE in
// `modules/core/tests/CoreTests.cpp`, which nulls this function against it in double precision.
inline float firDot (const float* a, const float* b, int len) noexcept
{
    FELITRONICS_FIR_NO_FMA

#if defined(FELITRONICS_FIR_WASM)
    v128_t acc0 = wasm_f32x4_const_splat (0.0f);
    for (int i = 0; i < len; i += 4)
    {
        const v128_t p = wasm_f32x4_mul (wasm_v128_load (a + i), wasm_v128_load (b + i));
        acc0 = wasm_f32x4_add (acc0, p);
    }
    float s[4];
    wasm_v128_store (s, acc0);
    return (s[0] + s[1]) + (s[2] + s[3]);

#elif defined(FELITRONICS_FIR_SSE2)
    __m128 acc0 = _mm_setzero_ps();
    for (int i = 0; i < len; i += 4)
    {
        const __m128 p = _mm_mul_ps (_mm_loadu_ps (a + i), _mm_loadu_ps (b + i));
        acc0 = _mm_add_ps (acc0, p);
    }
    float s[4];
    _mm_storeu_ps (s, acc0);
    return (s[0] + s[1]) + (s[2] + s[3]);

#elif defined(FELITRONICS_FIR_NEON)
    float32x4_t acc0 = vdupq_n_f32 (0.0f);
    for (int i = 0; i < len; i += 4)
    {
        const float32x4_t p = vmulq_f32 (vld1q_f32 (a + i), vld1q_f32 (b + i));
        acc0 = vaddq_f32 (acc0, p);
    }
    float s[4];
    vst1q_f32 (s, acc0);
    return (s[0] + s[1]) + (s[2] + s[3]);

#else
    // The scalar kernel is one of the four, not a fallback with its own arithmetic: four
    // accumulators and the same reduction tree, so a row with no SIMD returns the same bits as one
    // with it. The products go to NAMED variables before the adds — under `-ffp-contract=on` that
    // alone already forbids fusion, and the pragma above covers `fast`.
    float s0 = 0.0f, s1 = 0.0f, s2 = 0.0f, s3 = 0.0f;
    for (int i = 0; i < len; i += 4)
    {
        const float p0 = a[i + 0] * b[i + 0];
        const float p1 = a[i + 1] * b[i + 1];
        const float p2 = a[i + 2] * b[i + 2];
        const float p3 = a[i + 3] * b[i + 3];
        s0 = s0 + p0;
        s1 = s1 + p1;
        s2 = s2 + p2;
        s3 = s3 + p3;
    }
    return (s0 + s1) + (s2 + s3);
#endif
}

} // namespace felitronics::core

FELITRONICS_FIR_EXACT_END

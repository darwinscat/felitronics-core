// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// The wasm-audio tier built the way a consumer that wants Eigen's SSE path builds it: -msimd128 -msse4.1.
// Emscripten then defines __SSE__/__SSE2__ and lowers the SSE intrinsics onto SIMD128, but WebAssembly has
// no MXCSR, so its <xmmintrin.h> has no _mm_setcsr. FlushToZero.h used to key its x86 path on __SSE2__ and
// this file did not compile ("use of undeclared identifier '_mm_setcsr'"). Building it is half the gate.
// Running it is the other half: under these flags ScopedFlushToZero is the same no-op it is in a plain wasm
// build, wasm keeps its subnormals (there is no FTZ mode to switch on), and law 8's software flush is what
// zaps them (DSP-ARCHITECTURE.md §2: the wasm-audio tier "relies on software denormal flushing").

#if ! defined(__EMSCRIPTEN__) || ! defined(__wasm_simd128__) || ! defined(__SSE2__) || ! defined(__SSE4_1__)
  #error "this suite means something only in an emscripten build with -msimd128 -msse4.1"
#endif

#include <smmintrin.h>   // what Eigen includes for its SSE4.1 path, ahead of any core header

#include <felitronics_test.h>
#include <felitronics/core/FlushToZero.h>
#include <felitronics/core/PolyphaseFir.h>   // the other core header with a branch keyed on __SSE2__

#include <bit>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>

#if defined(FELITRONICS_FTZ_X86) || defined(FELITRONICS_FTZ_AARCH64)
  #error "wasm has no FP-control register: FlushToZero.h must not pick a hardware FTZ path here"
#endif
#if ! defined(FELITRONICS_FIR_WASM) || defined(FELITRONICS_FIR_SSE2)
  #error "firDot must stay on its SIMD128 kernel, the one the plain wasm build runs, not the SSE2 emulation"
#endif

using namespace felitronics;

int main()
{
    std::printf ("felitronics::core wasm + SSE-compat tests (-msimd128 -msse4.1)\n");

    test::group ("emscripten's MXCSR stand-in reports no FTZ");
    {
        test::ok ((static_cast<unsigned> (_mm_getcsr()) & _MM_FLUSH_ZERO_MASK) == 0u,
                  "_mm_getcsr() has FTZ off: there is no register behind it to turn on");
    }

    test::group ("ScopedFlushToZero is a no-op: subnormals survive inside its scope");
    {
        core::ScopedFlushToZero ftz;
        volatile float smallest = std::numeric_limits<float>::min();
        const float half = smallest * 0.5f;                                    // FTZ would make this +0
        test::ok (std::fpclassify (half) == FP_SUBNORMAL, "FLT_MIN * 0.5 stays subnormal (no FTZ)");
        volatile float denorm = std::numeric_limits<float>::denorm_min();
        const float twice = denorm * 2.0f;                                     // DAZ would read the input as 0
        test::ok (twice == 2.0f * std::numeric_limits<float>::denorm_min() && twice != 0.0f,
                  "a subnormal operand is read as itself (no DAZ)");
        volatile double smallestD = std::numeric_limits<double>::min();
        const double halfD = smallestD * 0.5;
        test::ok (std::fpclassify (halfD) == FP_SUBNORMAL, "DBL_MIN * 0.5 stays subnormal (no FTZ)");
    }

    test::group ("law 8's software flush covers what the hardware does not");
    {
        volatile float smallest = std::numeric_limits<float>::min();
        float x = smallest * 0.5f;
        core::flushDenormal (x);
        test::ok (std::bit_cast<std::uint32_t> (x) == 0u, "flushDenormal: subnormal float -> +0 exactly");
        float p = smallest * 0.5f;
        core::flushPoison (p);
        test::ok (std::bit_cast<std::uint32_t> (p) == 0u, "flushPoison: subnormal float -> +0 exactly");
        volatile double smallestD = std::numeric_limits<double>::min();
        double d = smallestD * 0.5;
        core::flushDenormal (d);
        test::ok (std::bit_cast<std::uint64_t> (d) == 0u, "flushDenormal: subnormal double -> +0 exactly");
        float keep = 0.25f;
        core::flushDenormal (keep);
        test::ok (keep == 0.25f, "a normal value passes through");
    }

    return test::report();
}

// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// felitronics::core::firDot — P56. The suite that makes the claim in PolyphaseFir.h checkable instead
// of asserted. Four things are gated here and each one fails differently:
//
//   * THE HASH pins the exact bits of 4096 inner products against a constant measured on five rows.
//     Any change of kernel, order, accumulator count or contraction setting moves it, on any row.
//   * THE CONTRACTION GATE fails on a row where the no-FMA pragma did not take — it is the only thing
//     standing between "we turned contraction off" and "we believe we turned contraction off".
//   * THE ORDER ORACLE re-spells the stated summation independently, so a kernel that is merely
//     ACCURATE but differently associated is caught. A double-precision null would pass that.
//   * THE ACCURACY ORACLE is the naive sequential loop the three call sites used to write by hand. It
//     lives ONLY here — the product contains one implementation of this sum and no second one.

#include <felitronics/core/PolyphaseFir.h>

#include "felitronics_test.h"

#include <cstdint>
#include <cstring>
#include <vector>

namespace test = felitronics::test;
using felitronics::core::firDot;
using felitronics::core::firPadLen;

static std::uint32_t bits (float f) noexcept { std::uint32_t u; std::memcpy (&u, &f, 4); return u; }

//==============================================================================
// A libm-free, bit-reproducible corpus: an LCG, a 24-bit integer, an exact power-of-two scale. Nothing
// in here can differ between rows, so a hash mismatch can only be the kernel. Magnitudes are ordinary
// audio, so no product and no partial sum can reach the subnormal range — which is what keeps the hash
// independent of FTZ, the one runtime state that would otherwise split wasm from every native row.
struct Lcg
{
    std::uint32_t s;
    explicit Lcg (std::uint32_t seed) : s (seed) {}
    std::uint32_t next() noexcept { s = s * 1664525u + 1013904223u; return s; }
    float unit()  noexcept { return (float) (std::int32_t) ((next() >> 8) ^ 0x800000u) * (1.0f / 8388608.0f) - 1.0f; }
    float coeff() noexcept { return unit() * (1.0f / 64.0f); }
};

// THE ACCURACY ORACLE — the sequential loop, in double.
static double oracleDot (const float* a, const float* b, int len, double& absSum) noexcept
{
    double acc = 0.0; absSum = 0.0;
    for (int i = 0; i < len; ++i)
    {
        const double t = (double) a[i] * (double) b[i];
        acc += t; absSum += (t < 0 ? -t : t);
    }
    return acc;
}

// THE ORDER ORACLE — the stated order re-spelled in float, independently of the kernel, so that a
// kernel which is merely as ACCURATE but differently associated is caught. This TU carries NO
// contraction pragma: it compiles under the tree's stated `-ffp-contract=on`, where fusing across two
// statements is not allowed on any of the three compilers (measured), so this oracle does not fuse
// either and the comparison is exact. Under a consumer's `-ffp-contract=fast` it WOULD fuse and this
// check would go red while firDot stayed right — which is a true report about the build, not a false
// one about the kernel.
static float orderOracle (const float* a, const float* b, int len) noexcept
{
    float s0 = 0.0f, s1 = 0.0f, s2 = 0.0f, s3 = 0.0f;
    for (int i = 0; i < len; i += 4)
    {
        const float p0 = a[i + 0] * b[i + 0];
        const float p1 = a[i + 1] * b[i + 1];
        const float p2 = a[i + 2] * b[i + 2];
        const float p3 = a[i + 3] * b[i + 3];
        s0 = s0 + p0; s1 = s1 + p1; s2 = s2 + p2; s3 = s3 + p3;
    }
    return (s0 + s1) + (s2 + s3);
}

// A SEQUENTIAL float sum — the shape the three call sites used before P56. It exists to prove the
// order gate can FAIL: if the two oracles agreed on the corpus, the order check would be vacuous.
static float sequentialOracle (const float* a, const float* b, int len) noexcept
{
    float acc = 0.0f;
    for (int i = 0; i < len; ++i) { const float p = a[i] * b[i]; acc = acc + p; }
    return acc;
}

//==============================================================================
int main()
{
    std::printf ("felitronics::core::firDot — P56 (the one polyphase FIR kernel)\n");

    std::vector<float> a (4096), b (4096);
    { Lcg rng (0x9E3779B9u); for (auto& v : a) v = rng.coeff(); for (auto& v : b) v = rng.unit(); }

    // --- the five-row hash ------------------------------------------------------------------------
    {
        test::group ("firDot: the bits, pinned against five rows");
        std::uint64_t h = 1469598103934665603ull;
        double worstRel = 0.0;
        int orderMatches = 0, sequentialDiffers = 0, n = 0;
        for (int len = 4; len <= 256; len += 4)
            for (int off = 0; off < 64; ++off)
            {
                const float* pa = &a[(std::size_t) off];
                const float* pb = &b[(std::size_t) off * 2];
                const float got = firDot (pa, pb, len);
                const std::uint32_t u = bits (got);
                for (int k = 0; k < 4; ++k) { h ^= (u >> (8 * k)) & 0xffu; h *= 1099511628211ull; }

                double absSum = 0.0;
                const double want = oracleDot (pa, pb, len, absSum);
                const double err  = (double) got - want < 0 ? want - (double) got : (double) got - want;
                if (absSum > 0.0 && err / absSum > worstRel) worstRel = err / absSum;

                if (bits (orderOracle (pa, pb, len)) == u) ++orderMatches;
                if (bits (sequentialOracle (pa, pb, len)) != u) ++sequentialDiffers;
                ++n;
            }

        // THE constant. Measured 2026-09-12 on all five rows of the acceptance: win (MSVC 19.44, SSE2),
        // deb (gcc 14.2 x86-64, SSE2), mac (Apple clang 14.0.3 x86-64, SSE2), docker linux/arm64
        // (gcc 14.4, NEON) and wasm (emscripten 6.0.9 in node — the scalar kernel AND, separately,
        // -msimd128). If this line ever needs changing, the summation order changed with it and every
        // consumer's output moved: say so out loud rather than re-measuring the constant.
        test::ok (h == 0xc197e85e7d8e9897ull,
                  "the corpus hash is the cross-row constant c197e85e7d8e9897 (got "
                      + std::to_string (h) + " decimal)");
        test::ok (worstRel < 4.0e-7,
                  "and it is ACCURATE: worst error over sum|a*b| is " + std::to_string (worstRel));
        test::ok (orderMatches == n,
                  "every result equals the stated order re-spelled independently ("
                      + std::to_string (orderMatches) + "/" + std::to_string (n) + ")");
        // Without this the order check above could pass on a kernel that simply summed sequentially.
        test::ok (sequentialDiffers > n / 2,
                  "and the corpus really distinguishes the orders: a sequential sum differs on "
                      + std::to_string (sequentialDiffers) + "/" + std::to_string (n));
    }

    // --- the contraction gate ---------------------------------------------------------------------
    {
        test::group ("firDot: no fused multiply-add, on this row");
        // a = b = 1 + 2^-12, so a*b is exactly 1 + 2^-11 + 2^-24. Rounded to float (ties-to-even, and
        // 2^-24 is exactly half an ulp of 1.0 with an even neighbour) that is 1 + 2^-11; adding
        // c = -(1 + 2^-11) gives EXACTLY zero. Fused, the product keeps its 2^-24 and the answer is
        // 2^-24 = 0x33800000 instead. The operands live in vectors so no row can fold this at compile
        // time and answer from its own constant evaluator instead of from the kernel.
        std::vector<float> ga (8, 0.0f), gb (8, 0.0f);
        ga[0] = -(1.0f + 1.0f / 2048.0f); gb[0] = 1.0f;                 // seeds lane 0 with c
        ga[4] =  (1.0f + 1.0f / 4096.0f); gb[4] = 1.0f + 1.0f / 4096.0f; // then adds a*b into it
        const float gate = firDot (ga.data(), gb.data(), 8);
        test::ok (bits (gate) == 0u,
                  "a*b + c is two roundings, not one — this row did NOT fuse (got "
                      + std::to_string (bits (gate)) + ")");
        // The witness that the fixture can see a fusion at all: the same operands one rounding apart.
        const float fused = (float) ((double) ga[4] * (double) gb[4] + (double) ga[0]);
        test::ok (bits (fused) != 0u, "and the fixture is not blind: the FUSED answer is not zero");
    }

    // --- the zero padding is a numeric no-op ------------------------------------------------------
    {
        test::group ("firDot: +0.0f padding changes no finite answer");
        std::vector<float> pa (a.begin(), a.begin() + 36), pb (b.begin(), b.begin() + 36);
        for (int i = 32; i < 36; ++i) pa[(std::size_t) i] = 0.0f;       // the pad prepare() writes
        const float withPad = firDot (pa.data(), pb.data(), 36);
        const float without = firDot (pa.data(), pb.data(), 32);
        test::ok (bits (withPad) == bits (without),
                  "four +0.0f coefficients are bit-for-bit invisible");
        test::ok (firPadLen (12) == 12 && firPadLen (13) == 16 && firPadLen (64) == 64
                    && firPadLen (1) == 4 && firPadLen (0) == 0,
                  "firPadLen rounds up to a multiple of four and leaves one alone");
    }

    // --- degenerate and edge lengths ---------------------------------------------------------------
    {
        test::group ("firDot: edges");
        std::vector<float> z (64, 0.0f);
        test::ok (bits (firDot (z.data(), z.data(), 64)) == 0u,
                  "all-zero in, +0.0f out (not -0.0f: the accumulators start at +0)");
        const float one4 = firDot (a.data(), b.data(), 4);
        double absSum = 0.0;
        const double want4 = oracleDot (a.data(), b.data(), 4, absSum);
        test::approx ((double) one4, want4, 1.0e-7, "the shortest legal length, 4, is exact enough");
        // 1024 is past every tap count the tree asks for (64 taps x 8x = 512 in the decimator).
        test::ok (std::isfinite (firDot (a.data(), b.data(), 1024)), "and a long one stays finite");
    }

    // --- which kernel this row built ---------------------------------------------------------------
    {
#if defined(FELITRONICS_FIR_WASM)
        const char* k = "wasm-simd128";
#elif defined(FELITRONICS_FIR_SSE2)
        const char* k = "sse2";
#elif defined(FELITRONICS_FIR_NEON)
        const char* k = "neon";
#else
        const char* k = "scalar";
#endif
        std::printf ("       kernel on this row: %s\n", k);
    }

    return test::report();
}

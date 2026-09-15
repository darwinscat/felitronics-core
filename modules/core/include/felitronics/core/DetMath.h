// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <cmath>
#include <cstdint>
#include <limits>

namespace felitronics::core::det
{

//==============================================================================
// felitronics::core::det — transcendental functions that give THE SAME BITS on every row we build on.
// Not more accurate than the system libm: THE SAME as itself, everywhere. That is a different property
// and it is the one the offline analyzers need, because their outputs are compared byte for byte
// (native CLI against the wasm module) and because a threshold that moves by one ulp does not move a
// reported number by one ulp — it flips a decision.
//
// WHY THE SYSTEM libm CANNOT DO THIS — measured on this tree, three rows (Apple clang/arm64,
// gcc 14 + glibc, emcc 6.0.9 + musl), over the ranges these analyzers actually use:
//
//     function                     points    Apple/musl   Apple/glibc   glibc/musl
//     log10 (dB of a ratio)        200000           236          4325         4278
//     tan   (filter prewarp)        20918          7259          7274          656
//     pow10 (dB thresholds)        132001         54356         54342           73
//
// pow10 is the worst at 41 %, and it is what every dB THRESHOLD is built from. glibc and musl agree
// almost everywhere; the outlier is Apple — so the developer's Mac computes something the CI row and
// the browser do not, and nothing today would notice. The Hann window built with std::cos differs in
// 502 of 16384 coefficients (order 14) and 4032 of 131072 (order 17) between Apple and musl, and every
// spectral bin is multiplied by those.
//
// THE PIN, AND WHY IT IS `volatile` AND NOT A PRAGMA. This library builds with -ffp-contract=on on
// purpose (law 10, top-level CMakeLists.txt). Under contraction the compiler fuses `a*b + c` into one
// FMA, which rounds once instead of twice — so a hand-written polynomial stops being deterministic:
// measured 525 differences of 16384 between a contract=on and a contract=off build OF THE SAME CODE ON
// ONE MACHINE. Worse, wasm32 has no scalar FMA at all, so it can never fuse: the native row and the
// wasm row diverge BY CONSTRUCTION. Every multiply-add here therefore goes through mulAdd(), which
// pins the rounding with a volatile temporary — the technique already in analysis::StereoSums::add for
// the same reason. Measured on the canonical FMA detector `a*b + (-round(a*b))`:
//
//     row                  canary at -ffp-contract=on   #pragma STDC FP_CONTRACT OFF   volatile
//     clang Apple arm64    fuses                        holds                          holds
//     gcc 14 x86-64        fuses                        IGNORED — fuses                holds
//     emcc wasm32          cannot fuse                  holds                          holds
//
// gcc ignores the STDC pragma in C++ (DSP-ARCHITECTURE.md §2 law 10 says so; this measured it), so the
// pragma is not an option for a header that must hold under all three. DetMathTests.cpp builds itself
// twice on one machine, contraction on and off, and REQUIRES an unpinned canary to differ — a run where
// it does not is reported as inconclusive rather than green, because that is the shape in which this
// gate would otherwise quietly stop being a gate.
//
// COST, measured: ~5.5 ns per det::cos against ~1.9 ns for the system one, i.e. x2.9 per call — but
// every call site in this library is a COEFFICIENT site (a window built once in prepare(), log2(N)
// twiddle seeds per transform, a handful of thresholds per setParams, one log10 per reported value),
// never a per-sample loop. On one order-17 window that is 0.71 ms against 0.24 ms, once per file,
// inside an analyzer run of 20-30 ms. It must never be put on a per-sample path; that is what the
// native path is for.
//
// DOMAIN. Argument reduction here is Cody-Waite with a two-part pi/2, which is exact only while
// |k * pio2Hi| is exact — that holds comfortably to |x| < 2^24. Outside that this returns NaN rather
// than a plausible number: a value that reads as an answer where none was computed is the one thing
// this library refuses everywhere else too.

//==============================================================================
// a*b + c that is NEVER fused. The volatile store is the pin; see the table above for why a pragma is
// not enough. Everything in this header routes its arithmetic through here.
inline double mulAdd (double a, double b, double c) noexcept
{
    volatile double t = a * b;
    return t + c;
}

// a*b that is never absorbed into a following add by contraction.
inline double mul (double a, double b) noexcept
{
    volatile double t = a * b;
    return t;
}

namespace detail
{
    // pi/2 split so that k * kPio2Hi is EXACT for every |k| < 2^24: kPio2Hi carries 33 significant bits,
    // leaving 20 for k. The residual lives in Lo (and Lo2 for the last bits of pi/2).
    // pi/2 as THREE doubles that sum to it to 123.6 bits, generated by tools/det-oracle/gen.py (mpmath,
    // 120 digits) and NOT copied from another library's table: p1 and p2 each carry 33 significant bits,
    // so k*p1 and k*p2 are EXACT for every |k| < 2^20, and p3 holds the rest.
    // ⚠ The first draft of this header took p1 from fdlibm's 33-bit split but p2 from its 53-bit one and
    // then added p3 (the residual after a 33-bit p2) — counting the tail twice. It cost ~1.6e11 ulp on
    // sin(pi_d), which is the window's own argument at i = N/4 and 3N/4, where the result decides whether
    // the coefficient is exactly 0.5 or one ulp below. Constants that pair must be generated together.
    inline constexpr double kPio2Hi  = 1.57079632673412561417e+00;   // 0x3ff921fb54400000, 33 bits
    inline constexpr double kPio2Lo  = 6.07710050630396597660e-11;   // 0x3dd0b4611a600000, 33 bits
    inline constexpr double kPio2Lo2 = 2.02226624879595063154e-21;   // 0x3ba3198a2e037073, the rest
    inline constexpr double k2OverPi = 6.36619772367581382433e-01;
    inline constexpr double kMaxArg  = 16777216.0;                   // 2^24

    // sin(r) on |r| <= pi/4, Taylor to r^15 (truncation below 5e-17 there).
    inline double polySin (double r) noexcept
    {
        const double r2 = mul (r, r);
        double p = -7.6471637318198164e-13;                // -1/15!
        p = mulAdd (p, r2,  1.6059043836821613e-10);       // +1/13!
        p = mulAdd (p, r2, -2.5052108385441720e-08);       // -1/11!
        p = mulAdd (p, r2,  2.7557319223985893e-06);       // +1/9!
        p = mulAdd (p, r2, -1.9841269841269841e-04);       // -1/7!
        p = mulAdd (p, r2,  8.3333333333333333e-03);       // +1/5!
        p = mulAdd (p, r2, -1.6666666666666666e-01);       // -1/3!
        return mulAdd (mul (r, r2), p, r);                 // r + r^3 * p
    }

    // cos(r) on |r| <= pi/4, Taylor to r^16.
    inline double polyCos (double r) noexcept
    {
        const double r2 = mul (r, r);
        double p =  4.7794773323873853e-14;                // +1/16!
        p = mulAdd (p, r2, -1.1470745597729725e-11);       // -1/14!
        p = mulAdd (p, r2,  2.0876756987868099e-09);       // +1/12!
        p = mulAdd (p, r2, -2.7557319223985893e-07);       // -1/10!
        p = mulAdd (p, r2,  2.4801587301587302e-05);       // +1/8!
        p = mulAdd (p, r2, -1.3888888888888889e-03);       // -1/6!
        p = mulAdd (p, r2,  4.1666666666666666e-02);       // +1/4!
        p = mulAdd (p, r2, -5.0000000000000000e-01);       // -1/2!
        return mulAdd (r2, p, 1.0);                        // 1 + r^2 * p
    }

    // Cody-Waite: x = k*(pi/2) + r with |r| <= pi/4. Returns the quadrant in `quadrant`.
    inline double reduce (double x, int& quadrant) noexcept
    {
        const double q = mul (x, k2OverPi);
        const double k = (double) (long long) (q + (q >= 0.0 ? 0.5 : -0.5));
        double r = mulAdd (-k, kPio2Hi,  x);
        r        = mulAdd (-k, kPio2Lo,  r);
        r        = mulAdd (-k, kPio2Lo2, r);
        quadrant = (int) (((long long) k) & 3);
        return r;
    }
} // namespace detail

// cos/sin for |x| < 2^24. Outside that, and for a non-finite argument, NaN — never a plausible number.
inline double cos (double x) noexcept
{
    if (! (std::fabs (x) < detail::kMaxArg)) return std::numeric_limits<double>::quiet_NaN();
    int q = 0;
    const double r = detail::reduce (x, q);
    switch (q)
    {
        case 0:  return  detail::polyCos (r);
        case 1:  return -detail::polySin (r);
        case 2:  return -detail::polyCos (r);
        default: return  detail::polySin (r);
    }
}

inline double sin (double x) noexcept
{
    if (! (std::fabs (x) < detail::kMaxArg)) return std::numeric_limits<double>::quiet_NaN();
    int q = 0;
    const double r = detail::reduce (x, q);
    switch (q)
    {
        case 0:  return  detail::polySin (r);
        case 1:  return  detail::polyCos (r);
        case 2:  return -detail::polySin (r);
        default: return -detail::polyCos (r);
    }
}

} // namespace felitronics::core::det

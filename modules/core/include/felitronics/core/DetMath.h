// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <bit>
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


//==============================================================================
// THE LOGARITHMIC AND EXPONENTIAL HALF.
//
// These are where the divergence is worst, and it is not the dB VALUES that make it dangerous — it is
// the dB THRESHOLDS. Measured over -120..+12 dB, 41 % of pow10 results differ between Apple's libm and
// both Linux ones; every silenceThresholdDb, quietThresholdDb, enterDb, exitDb and minDropDb in this
// library is built from one of them. A reported number that moves by an ulp moves by an ulp. A
// threshold that moves by an ulp moves a DECISION: a block lands inside the gate on one row and outside
// it on the next, and the two rows then disagree by a whole block of loudness, not by a bit.

namespace detail
{
    inline constexpr double kLog2E     = 1.44269504088896338700e+00;
    inline constexpr double kLog10E    = 4.34294481903251816668e-01;
    inline constexpr double kLog10_2Hi = 3.01029920578002929688e-01;   // 20 bits: e * Hi is EXACT for |e| < 2^32
    inline constexpr double kLog10_2Lo = 7.50859782655262350275e-08;
    inline constexpr double kLog2_10Hi = 3.32192808389663696289e+00;   // 26 bits, for the Dekker product
    inline constexpr double kLog2_10Lo = 1.09907253849796945458e-08;
    inline constexpr double kLn2       = 6.93147180559945286227e-01;
    inline constexpr double kSqrt2     = 1.41421356237309514547e+00;
    inline constexpr double kSplit     = 134217729.0;                  // 2^27 + 1, Veltkamp

    // atanh(f) * 2 on |f| <= (sqrt2-1)/(sqrt2+1) = 0.1716, series to f^21 (tail below 1e-17).
    inline double log1pSeries (double f) noexcept
    {
        const double f2 = mul (f, f);
        double p = 1.0 / 21.0;
        p = mulAdd (p, f2, 1.0 / 19.0);
        p = mulAdd (p, f2, 1.0 / 17.0);
        p = mulAdd (p, f2, 1.0 / 15.0);
        p = mulAdd (p, f2, 1.0 / 13.0);
        p = mulAdd (p, f2, 1.0 / 11.0);
        p = mulAdd (p, f2, 1.0 /  9.0);
        p = mulAdd (p, f2, 1.0 /  7.0);
        p = mulAdd (p, f2, 1.0 /  5.0);
        p = mulAdd (p, f2, 1.0 /  3.0);
        p = mulAdd (p, f2, 1.0);
        return mul (2.0, mul (f, p));
    }

    // x = m * 2^e with m in [sqrt2/2, sqrt2), so |(m-1)/(m+1)| <= 0.1716 and the series above converges.
    inline void decompose (double x, double& m, int& e) noexcept
    {
        const std::uint64_t b = std::bit_cast<std::uint64_t> (x);
        e = (int) ((b >> 52) & 0x7ff) - 1022;
        m = std::bit_cast<double> ((b & 0x800fffffffffffffull) | 0x3fe0000000000000ull);   // m in [0.5, 1)
        if (m < kSqrt2 * 0.5) { m = mul (m, 2.0); --e; }
    }

    // Veltkamp split: x = hi + lo with hi carrying 26 significant bits, so hi * (a 26-bit constant) is
    // EXACT. The subtractions cannot be contracted (no multiply-add pattern), so no pin is needed here.
    inline void split (double x, double& hi, double& lo) noexcept
    {
        const double t = mul (kSplit, x);
        hi = t - (t - x);
        lo = x - hi;
    }

    // 2^f on |f| <= 0.5, Taylor in f*ln2 to the 13th term.
    inline double exp2Frac (double f) noexcept
    {
        const double y = mul (f, kLn2);
        double p = 1.0 / 6227020800.0;                      // 1/13!
        p = mulAdd (p, y, 1.0 / 479001600.0);               // 1/12!
        p = mulAdd (p, y, 1.0 / 39916800.0);                // 1/11!
        p = mulAdd (p, y, 1.0 / 3628800.0);                 // 1/10!
        p = mulAdd (p, y, 1.0 / 362880.0);                  // 1/9!
        p = mulAdd (p, y, 1.0 / 40320.0);                   // 1/8!
        p = mulAdd (p, y, 1.0 / 5040.0);                    // 1/7!
        p = mulAdd (p, y, 1.0 / 720.0);                     // 1/6!
        p = mulAdd (p, y, 1.0 / 120.0);                     // 1/5!
        p = mulAdd (p, y, 1.0 / 24.0);                      // 1/4!
        p = mulAdd (p, y, 1.0 / 6.0);                       // 1/3!
        p = mulAdd (p, y, 0.5);                             // 1/2!
        p = mulAdd (p, y, 1.0);
        return mulAdd (p, y, 1.0);
    }

    // 2^k as an exact double, k an integer in [-1074, 1023]. Subnormal results are built by halving.
    inline double pow2i (int k) noexcept
    {
        if (k >= 1024)  return std::numeric_limits<double>::infinity();
        if (k >= -1022) return std::bit_cast<double> ((std::uint64_t) (k + 1023) << 52);
        if (k < -1074)  return 0.0;
        return std::bit_cast<double> ((std::uint64_t) 1 << (k + 1074));   // subnormal 2^k, exact
    }
} // namespace detail

// log2 / log10 for x > 0. Zero, negative and non-finite are refused, never answered with a plausible
// number: -inf for log(0) is the mathematically right answer, and it is given; a NaN argument stays NaN.
inline double log2 (double x) noexcept
{
    if (std::isnan (x)) return x;
    if (x < 0.0)  return std::numeric_limits<double>::quiet_NaN();
    if (x == 0.0) return -std::numeric_limits<double>::infinity();
    if (std::isinf (x)) return x;
    if (x < 2.2250738585072014e-308) return log2 (mul (x, 18014398509481984.0)) - 54.0;   // subnormal
    double m; int e;
    detail::decompose (x, m, e);
    const double f = (m - 1.0) / (m + 1.0);
    return mulAdd (detail::log1pSeries (f), detail::kLog2E, (double) e);
}

inline double log10 (double x) noexcept
{
    if (std::isnan (x)) return x;
    if (x < 0.0)  return std::numeric_limits<double>::quiet_NaN();
    if (x == 0.0) return -std::numeric_limits<double>::infinity();
    if (std::isinf (x)) return x;
    if (x < 2.2250738585072014e-308) return log10 (mul (x, 18014398509481984.0)) - 16.25561652641961;
    double m; int e;
    detail::decompose (x, m, e);
    const double f = (m - 1.0) / (m + 1.0);
    const double lm = mul (detail::log1pSeries (f), detail::kLog10E);
    // e * Hi is exact (20-bit Hi), so the only rounding is in the two sums.
    return mulAdd ((double) e, detail::kLog10_2Hi, mulAdd ((double) e, detail::kLog10_2Lo, lm));
}

inline double exp2 (double x) noexcept
{
    if (std::isnan (x)) return x;
    if (x >= 1024.0) return std::numeric_limits<double>::infinity();
    if (x <= -1075.0) return 0.0;
    const double kd = (double) (long long) (x + (x >= 0.0 ? 0.5 : -0.5));
    const double f  = x - kd;
    return mul (detail::exp2Frac (f), detail::pow2i ((int) kd));
}

// 10^x. The exponent is formed as an EXACT Dekker product before it reaches exp2: a plain
// exp2(x * log2(10)) loses the low bits of the product, and at 60 dB that is ~14 ulp in the result —
// which is the whole point, since this function's output is a threshold.
inline double pow10 (double x) noexcept
{
    if (std::isnan (x)) return x;
    if (x >= 308.26) return std::numeric_limits<double>::infinity();
    if (x <= -324.0) return 0.0;
    double xh, xl;
    detail::split (x, xh, xl);
    const double hi = mul (xh, detail::kLog2_10Hi);                        // exact: 26 bits * 26 bits
    const double lo = mulAdd (xh, detail::kLog2_10Lo,
                       mulAdd (xl, detail::kLog2_10Hi, mul (xl, detail::kLog2_10Lo)));
    const double t  = hi + lo;
    const double e  = (hi - t) + lo;                                       // the bits the sum dropped
    // 2^(t+e) = 2^t * 2^e, and |e| is tiny, so two terms of 2^e are more than enough.
    const double corr = mulAdd (mul (e, detail::kLn2), mulAdd (e, mul (0.5, detail::kLn2), 1.0), 1.0);
    return mul (exp2 (t), corr);
}

// tan from the same reduction. Poles are poles: tan(pi/2_d) is a huge finite number, as it should be.
inline double tan (double x) noexcept
{
    if (! (std::fabs (x) < detail::kMaxArg)) return std::numeric_limits<double>::quiet_NaN();
    int q = 0;
    const double r = detail::reduce (x, q);
    const double s = detail::polySin (r), c = detail::polyCos (r);
    return (q & 1) ? -(c / s) : (s / c);
}


// x^y for x > 0, formed as exp2(y * log2 x) with the product taken exactly (the same Dekker split as
// pow10 uses, and for the same reason). Only two call sites need a general power — the K-weighting
// shelf constants — and both have constant operands, so this is a design-time function, never a hot one.
inline double pow (double x, double y) noexcept
{
    if (std::isnan (x) || std::isnan (y)) return std::numeric_limits<double>::quiet_NaN();
    if (y == 0.0) return 1.0;
    if (x <= 0.0) return std::numeric_limits<double>::quiet_NaN();   // negative bases are not needed here
    const double l = log2 (x);
    double yh, yl, lh, ll;
    detail::split (y, yh, yl);
    detail::split (l, lh, ll);
    const double hi = mul (yh, lh);
    const double lo = mulAdd (yh, ll, mulAdd (yl, lh, mul (yl, ll)));
    const double t  = hi + lo;
    const double e  = (hi - t) + lo;
    const double corr = mulAdd (mul (e, detail::kLn2), mulAdd (e, mul (0.5, detail::kLn2), 1.0), 1.0);
    return mul (exp2 (t), corr);
}

} // namespace felitronics::core::det


namespace felitronics::core
{

//==============================================================================
// THE TWO COEFFICIENT-MATH POLICIES, and why they are a TYPE and not a flag.
//
// A filter's coefficients decide its bits; its bits decide a product's sound. eq::Svf is TabbyEQ's and
// OrbitCab's core, so its numbers must not move — while the offline analyzers that borrow the same
// filter need coefficients that are the same on every row. Both are true at once, so the choice travels
// with the TYPE: eq::Svf is and stays BasicSvf<core::SystemMath>, and an analyzer spells
// BasicSvf<core::DetMath> where it wants reproducibility.
//
// There is deliberately NO default on the offline side and no runtime switch: a compile-time flag would
// reroute every consumer in a build at once (and fcore_measure links both regimes into ONE binary — its
// `lufs` mode is the system meter and its `report` mode is the deterministic one), while a
// setParamsDet() would be a second copy of the coefficient formula, which is the drift these policies
// exist to prevent. static_assert in the suites pins both directions, so changing an alias is a
// deliberate act that must also edit a test.
struct SystemMath
{
    static double tan  (double x) noexcept { return std::tan (x); }
    static double sqrt (double x) noexcept { return std::sqrt (x); }
    static double pow10 (double x) noexcept { return std::pow (10.0, x); }
    static double log10 (double x) noexcept { return std::log10 (x); }
    static double pow  (double x, double y) noexcept { return std::pow (x, y); }
};

struct DetMath
{
    static double tan  (double x) noexcept { return det::tan (x); }
    static double sqrt (double x) noexcept { return std::sqrt (x); }   // IEEE-exact: no det version needed
    static double pow10 (double x) noexcept { return det::pow10 (x); }
    static double log10 (double x) noexcept { return det::log10 (x); }
    static double pow  (double x, double y) noexcept { return det::pow (x, y); }
};

} // namespace felitronics::core

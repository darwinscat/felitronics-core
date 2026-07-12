// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// Theory-first falsification suite for core::flushPoison (+ its flushDenormal superset contract).
// Every expectation below was derived from the DOCUMENTED contract in FlushToZero.h plus IEEE-754
// theory BEFORE the test was run, and none may be weakened to make a run pass:
//   (a) poison — any NaN (quiet/signalling, either sign, any payload) and ±Inf — maps to EXACT 0;
//   (b) finite |x| < 1e-15 (the documented threshold; covers ±0, ALL subnormals, FLT_MIN/DBL_MIN
//       neighbourhoods — 1e-15 sits far above both binary32 and binary64 subnormal ranges) → EXACT 0;
//   (c) finite |x| >= 1e-15 ("at or above ... passes through bit-identically") → IDENTITY, sign
//       preserved bit-for-bit because identity is bitwise;
//   (d) totality: for EVERY input bit pattern the result is finite;
//   (e) idempotence f(f(x)) == f(x) bitwise — follows from the mapping: 0 re-flushes to 0, an
//       identity survivor is still at/above the threshold;
//   (f) strict superset: on every FINITE input flushPoison agrees with flushDenormal bit-for-bit
//       ("on healthy streams it behaves exactly like flushDenormal — a drop-in strict superset").
// The sign of the flushed zero is NOT promised by the doc ("exact 0"); current behaviour (+0) is
// pinned separately as a documented-behaviour pin so an accidental change is still caught.
//
// Sweep design: exhaustive BY CLASS — every exponent field value × both signs × mantissa edge
// patterns (zero / min / mid / all-ones / randoms) covers every FP class boundary exactly (±0,
// denorm_min → max subnormal, FLT_MIN/DBL_MIN, every normal binade, ±Inf, sNaN/qNaN payload min/max),
// plus millions of dense random bit patterns. Set FELITRONICS_EXHAUSTIVE=1 to additionally sweep all
// 2^32 float bit patterns (seconds in Release; skipped by default to keep the sanitizer gate quick).

#include <felitronics_test.h>
#include <felitronics/core/FlushToZero.h>

#include <bit>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <string>

using namespace felitronics;

//==============================================================================
// Doc-derived oracles (written from the contract text, not from the implementation).

static float  oracleF (float x)  { if (std::isnan (x) || std::isinf (x)) return 0.0f; return std::fabs (x) < 1e-15f ? 0.0f : x; }
static double oracleD (double x) { if (std::isnan (x) || std::isinf (x)) return 0.0;  return std::fabs (x) < 1e-15  ? 0.0  : x; }

struct SweepCounters
{
    long long cases = 0, totality = 0, mapping = 0, idempotence = 0, superset = 0, plusZeroPin = 0;
};

static void checkOneF (std::uint32_t bits, SweepCounters& c)
{
    const float x = std::bit_cast<float> (bits);
    float y = x; core::flushPoison (y);
    ++c.cases;

    if (! std::isfinite (y)) ++c.totality;                                    // (d)

    const float e = oracleF (x);
    const bool mustZap = std::isnan (x) || std::isinf (x) || std::fabs (x) < 1e-15f;
    if (mustZap)
    {
        if (! (y == 0.0f)) ++c.mapping;                                       // (a)/(b): exact zero
        if (std::bit_cast<std::uint32_t> (y) != 0u) ++c.plusZeroPin;          // pin: +0, not -0
    }
    else if (std::bit_cast<std::uint32_t> (y) != std::bit_cast<std::uint32_t> (e))
        ++c.mapping;                                                          // (c): bit-identity

    float z = y; core::flushPoison (z);                                       // (e)
    if (std::bit_cast<std::uint32_t> (z) != std::bit_cast<std::uint32_t> (y)) ++c.idempotence;

    if (std::isfinite (x))                                                    // (f)
    {
        float d = x; core::flushDenormal (d);
        if (std::bit_cast<std::uint32_t> (d) != std::bit_cast<std::uint32_t> (y)) ++c.superset;
    }
}

static void checkOneD (std::uint64_t bits, SweepCounters& c)
{
    const double x = std::bit_cast<double> (bits);
    double y = x; core::flushPoison (y);
    ++c.cases;

    if (! std::isfinite (y)) ++c.totality;

    const double e = oracleD (x);
    const bool mustZap = std::isnan (x) || std::isinf (x) || std::fabs (x) < 1e-15;
    if (mustZap)
    {
        if (! (y == 0.0)) ++c.mapping;
        if (std::bit_cast<std::uint64_t> (y) != 0ull) ++c.plusZeroPin;
    }
    else if (std::bit_cast<std::uint64_t> (y) != std::bit_cast<std::uint64_t> (e))
        ++c.mapping;

    double z = y; core::flushPoison (z);
    if (std::bit_cast<std::uint64_t> (z) != std::bit_cast<std::uint64_t> (y)) ++c.idempotence;

    if (std::isfinite (x))
    {
        double d = x; core::flushDenormal (d);
        if (std::bit_cast<std::uint64_t> (d) != std::bit_cast<std::uint64_t> (y)) ++c.superset;
    }
}

static void assertCounters (const SweepCounters& c, const std::string& name)
{
    test::ok (c.totality    == 0, name + ": totality — output always finite (violations: "        + std::to_string (c.totality)    + " of " + std::to_string (c.cases) + ")");
    test::ok (c.mapping     == 0, name + ": doc mapping poison/sub-threshold->0, else identity (" + std::to_string (c.mapping)     + " of " + std::to_string (c.cases) + ")");
    test::ok (c.idempotence == 0, name + ": idempotence f(f(x))==f(x) bitwise ("                  + std::to_string (c.idempotence) + " of " + std::to_string (c.cases) + ")");
    test::ok (c.superset    == 0, name + ": flushPoison==flushDenormal bitwise on finite x ("     + std::to_string (c.superset)    + " of " + std::to_string (c.cases) + ")");
    test::ok (c.plusZeroPin == 0, name + ": PIN flushed zero is +0 (not contract - doc says 'exact 0' with no sign promise; pinned so a silent change is caught)");
}

static std::uint32_t xorshift32 (std::uint32_t& s) { s ^= s << 13; s ^= s >> 17; s ^= s << 5; return s; }
static std::uint64_t xorshift64 (std::uint64_t& s) { s ^= s << 13; s ^= s >> 7;  s ^= s << 17; return s; }

//==============================================================================
int main()
{
    std::printf ("FlushPoisonTheoryTests (theory-first falsification: contract derived before running)\n");

    //==========================================================================
    test::group ("float: named class table (poison / zero / denormal / FLT_MIN / threshold / normals)");
    {
        const float T = 1e-15f;   // the documented threshold literal

        // Poison classes: qNaN/sNaN, both payload signs, payload min/max; +/-Inf.
        const std::uint32_t poison[] = {
            0x7FC00000u,            // qNaN canonical
            0xFFC00000u,            // -qNaN
            0x7FC00001u,            // qNaN payload min
            0x7FFFFFFFu,            // qNaN payload max
            0xFFFFFFFFu,            // -qNaN payload max
            0x7F800001u,            // sNaN payload min
            0x7FBFFFFFu,            // sNaN payload max
            0xFF800001u,            // -sNaN payload min
            0xFFBFFFFFu,            // -sNaN payload max
            0x7F800000u,            // +Inf
            0xFF800000u,            // -Inf
        };
        for (auto b : poison)
        {
            float y = std::bit_cast<float> (b); core::flushPoison (y);
            test::ok (y == 0.0f && std::isfinite (y), "poison bits 0x" + std::to_string (b) + " -> exact 0");
        }

        // Sub-threshold classes -> exact 0: +/-0, full subnormal range boundaries, FLT_MIN neighbourhood.
        const float subThresh[] = {
            0.0f, -0.0f,
            std::numeric_limits<float>::denorm_min(),  -std::numeric_limits<float>::denorm_min(),
            std::bit_cast<float> (0x007FFFFFu),        std::bit_cast<float> (0x807FFFFFu),        // max subnormal
            std::numeric_limits<float>::min(),         -std::numeric_limits<float>::min(),        // FLT_MIN (still << 1e-15)
            std::nextafterf (std::numeric_limits<float>::min(), 1.0f),                            // FLT_MIN + 1 ulp
            1e-30f, -1e-30f, 1e-16f, -1e-16f,
            std::nextafterf (T, 0.0f), -std::nextafterf (T, 0.0f),                                // threshold - 1 ulp
        };
        for (float v : subThresh)
        {
            float y = v; core::flushPoison (y);
            test::ok (y == 0.0f, "sub-threshold " + std::to_string (v) + " -> exact 0");
        }

        // At/above threshold -> bit-identity (doc: "at or above 1e-15 passes through bit-identically").
        const float identity[] = {
            T, -T,                                                                                 // exactly the threshold (inclusive)
            std::nextafterf (T, std::numeric_limits<float>::infinity()),                           // threshold + 1 ulp
            -std::nextafterf (T, std::numeric_limits<float>::infinity()),
            1e-14f, -1e-14f, 1e-7f, 0.5f, -0.5f, 1.0f, -1.0f, 123456.789f,
            1e30f, -1e30f, std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(),
        };
        for (float v : identity)
        {
            float y = v; core::flushPoison (y);
            test::ok (std::bit_cast<std::uint32_t> (y) == std::bit_cast<std::uint32_t> (v),
                      "identity (>= threshold) " + std::to_string (v) + " bit-preserved incl. sign");
        }
    }

    //==========================================================================
    test::group ("double: named class table (poison / zero / denormal / DBL_MIN / threshold / normals)");
    {
        const double T = 1e-15;

        const std::uint64_t poison[] = {
            0x7FF8000000000000ull, 0xFFF8000000000000ull,                       // +/- qNaN canonical
            0x7FF8000000000001ull, 0x7FFFFFFFFFFFFFFFull, 0xFFFFFFFFFFFFFFFFull,
            0x7FF0000000000001ull, 0x7FF7FFFFFFFFFFFFull,                       // sNaN payload min/max
            0xFFF0000000000001ull, 0xFFF7FFFFFFFFFFFFull,
            0x7FF0000000000000ull, 0xFFF0000000000000ull,                       // +/- Inf
        };
        for (auto b : poison)
        {
            double y = std::bit_cast<double> (b); core::flushPoison (y);
            test::ok (y == 0.0 && std::isfinite (y), "double poison -> exact 0");
        }

        const double subThresh[] = {
            0.0, -0.0,
            std::numeric_limits<double>::denorm_min(), -std::numeric_limits<double>::denorm_min(),
            std::bit_cast<double> (0x000FFFFFFFFFFFFFull), std::bit_cast<double> (0x800FFFFFFFFFFFFFull), // max subnormal
            std::numeric_limits<double>::min(), -std::numeric_limits<double>::min(),                     // DBL_MIN
            std::nextafter (std::numeric_limits<double>::min(), 1.0),
            1e-300, -1e-300, 1e-16, -1e-16,
            std::nextafter (T, 0.0), -std::nextafter (T, 0.0),                                            // threshold - 1 ulp
        };
        for (double v : subThresh)
        {
            double y = v; core::flushPoison (y);
            test::ok (y == 0.0, "double sub-threshold -> exact 0");
        }

        const double identity[] = {
            T, -T,
            std::nextafter (T, std::numeric_limits<double>::infinity()),
            -std::nextafter (T, std::numeric_limits<double>::infinity()),
            1e-14, -1e-14, 0.5, -1.0, 6.02e23,
            std::numeric_limits<double>::max(), -std::numeric_limits<double>::max(),
        };
        for (double v : identity)
        {
            double y = v; core::flushPoison (y);
            test::ok (std::bit_cast<std::uint64_t> (y) == std::bit_cast<std::uint64_t> (v),
                      "double identity (>= threshold) bit-preserved incl. sign");
        }
    }

    //==========================================================================
    test::group ("float: exhaustive-by-class sweep — every exponent x sign x mantissa edge patterns");
    {
        SweepCounters c;
        std::uint32_t rng = 0xC0FFEE01u;
        for (std::uint32_t sign = 0; sign <= 1; ++sign)
            for (std::uint32_t exp = 0; exp <= 255; ++exp)
            {
                const std::uint32_t mants[] = { 0x000000u, 0x000001u, 0x400000u, 0x7FFFFFu,
                                                xorshift32 (rng) & 0x7FFFFFu, xorshift32 (rng) & 0x7FFFFFu,
                                                xorshift32 (rng) & 0x7FFFFFu, xorshift32 (rng) & 0x7FFFFFu };
                for (auto m : mants)
                    checkOneF ((sign << 31) | (exp << 23) | m, c);
            }
        assertCounters (c, "float class sweep");
    }

    test::group ("double: exhaustive-by-class sweep — every exponent x sign x mantissa edge patterns");
    {
        SweepCounters c;
        std::uint64_t rng = 0xDEADBEEFCAFE0001ull;
        for (std::uint64_t sign = 0; sign <= 1; ++sign)
            for (std::uint64_t exp = 0; exp <= 2047; ++exp)
            {
                const std::uint64_t mants[] = { 0x0000000000000ull, 0x0000000000001ull,
                                                0x8000000000000ull, 0xFFFFFFFFFFFFFull,
                                                xorshift64 (rng) & 0xFFFFFFFFFFFFFull,
                                                xorshift64 (rng) & 0xFFFFFFFFFFFFFull,
                                                xorshift64 (rng) & 0xFFFFFFFFFFFFFull,
                                                xorshift64 (rng) & 0xFFFFFFFFFFFFFull };
                for (auto m : mants)
                    checkOneD ((sign << 63) | (exp << 52) | m, c);
            }
        assertCounters (c, "double class sweep");
    }

    //==========================================================================
    test::group ("dense random bit patterns (4M float + 4M double)");
    {
        SweepCounters cf;
        std::uint32_t r32 = 0x12345678u;
        for (int i = 0; i < 4'000'000; ++i) checkOneF (xorshift32 (r32), cf);
        assertCounters (cf, "float random");

        SweepCounters cd;
        std::uint64_t r64 = 0x9E3779B97F4A7C15ull;
        for (int i = 0; i < 4'000'000; ++i) checkOneD (xorshift64 (r64), cd);
        assertCounters (cd, "double random");
    }

    //==========================================================================
    // Optional: ALL 2^32 float bit patterns (Release: seconds). Off by default so sanitizer runs stay quick.
    if (const char* e = std::getenv ("FELITRONICS_EXHAUSTIVE"); e != nullptr && e[0] == '1')
    {
        test::group ("float: FULL exhaustive sweep over all 2^32 bit patterns (FELITRONICS_EXHAUSTIVE=1)");
        SweepCounters c;
        std::uint32_t b = 0;
        do { checkOneF (b, c); } while (++b != 0u);
        assertCounters (c, "float FULL 2^32");
    }

    return test::report();
}

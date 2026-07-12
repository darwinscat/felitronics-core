// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// Compact adversarial pins for flushPoison: named IEEE-754 edge cases from the external case list.

#include <felitronics_test.h>
#include <felitronics/core/FlushToZero.h>

#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>

using felitronics::test::group;
using felitronics::test::ok;

namespace
{
template <typename T>
auto bitsOf (T x)
{
    if constexpr (sizeof (T) == sizeof (std::uint32_t)) return std::bit_cast<std::uint32_t> (x);
    else return std::bit_cast<std::uint64_t> (x);
}
} // namespace

int main()
{
    std::printf ("felitronics::core flushPoison FALSIFY tests\n");

    group ("DeepSeek 6.1 materialized: qNaN/sNaN payloads of both signs flush to exact +0");
    {
        const std::uint32_t fbits[] { 0x7FC12345u, 0xFFC12345u, 0x7FA00001u, 0xFFA00001u };
        bool all = true;
        for (auto b : fbits)
        {
            float x = std::bit_cast<float> (b);
            felitronics::core::flushPoison (x);
            all = all && (bitsOf (x) == 0u);
        }
        const std::uint64_t dbits[] { 0x7FF8000000012345ull, 0xFFF8000000012345ull,
                                      0x7FF0000000000001ull, 0xFFF0000000000001ull };
        for (auto b : dbits)
        {
            double x = std::bit_cast<double> (b);
            felitronics::core::flushPoison (x);
            all = all && (bitsOf (x) == 0ull);
        }
        ok (all, "NaN payload/sign variants map to exact +0");
    }

    group ("DeepSeek 6.2 materialized: threshold is strict < 1e-15; at-threshold survives");
    {
        float fBelow = std::nextafterf (1.0e-15f, 0.0f);
        float fAt = 1.0e-15f;
        float fAbove = std::nextafterf (1.0e-15f, std::numeric_limits<float>::infinity());
        felitronics::core::flushPoison (fBelow);
        felitronics::core::flushPoison (fAt);
        felitronics::core::flushPoison (fAbove);
        ok (bitsOf (fBelow) == 0u, "float just below threshold flushes");
        ok (bitsOf (fAt) == bitsOf (1.0e-15f), "float exactly at threshold is bit-preserved");
        ok (bitsOf (fAbove) == bitsOf (std::nextafterf (1.0e-15f, std::numeric_limits<float>::infinity())),
            "float just above threshold is bit-preserved");

        double dBelow = std::nextafter (1.0e-15, 0.0);
        double dAt = 1.0e-15;
        double dAbove = std::nextafter (1.0e-15, std::numeric_limits<double>::infinity());
        felitronics::core::flushPoison (dBelow);
        felitronics::core::flushPoison (dAt);
        felitronics::core::flushPoison (dAbove);
        ok (bitsOf (dBelow) == 0ull, "double just below threshold flushes");
        ok (bitsOf (dAt) == bitsOf (1.0e-15), "double exactly at threshold is bit-preserved");
        ok (bitsOf (dAbove) == bitsOf (std::nextafter (1.0e-15, std::numeric_limits<double>::infinity())),
            "double just above threshold is bit-preserved");
    }

    group ("DeepSeek 6.3 materialized: infinities and negative zero pin");
    {
        float pf = std::numeric_limits<float>::infinity();
        float nf = -std::numeric_limits<float>::infinity();
        float nz = -0.0f;
        felitronics::core::flushPoison (pf);
        felitronics::core::flushPoison (nf);
        felitronics::core::flushPoison (nz);
        ok (bitsOf (pf) == 0u && bitsOf (nf) == 0u, "float infinities flush to exact +0");
        ok (bitsOf (nz) == 0u, "float -0 flushes to +0 (documented-behaviour pin)");

        double pd = std::numeric_limits<double>::infinity();
        double nd = -std::numeric_limits<double>::infinity();
        double ndz = -0.0;
        felitronics::core::flushPoison (pd);
        felitronics::core::flushPoison (nd);
        felitronics::core::flushPoison (ndz);
        ok (bitsOf (pd) == 0ull && bitsOf (nd) == 0ull, "double infinities flush to exact +0");
        ok (bitsOf (ndz) == 0ull, "double -0 flushes to +0 (documented-behaviour pin)");
    }

    return felitronics::test::report();
}

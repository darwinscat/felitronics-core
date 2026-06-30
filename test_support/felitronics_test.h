// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

// Tiny JUCE-free, zero-dependency test harness shared by every felitronics-core module's self-tests
// (the `teq` discipline: measured audio == analytic curve). Each module's test .cpp owns its main():
//   int main() { ...checks...; return felitronics::test::report(); }

#include <cmath>
#include <cstdio>
#include <string>

namespace felitronics::test
{
    struct Stats { int checks = 0, failures = 0; };
    inline Stats& stats() { static Stats s; return s; }

    inline void ok (bool cond, const std::string& msg)
    {
        ++stats().checks;
        if (! cond) { ++stats().failures; std::printf ("    FAIL: %s\n", msg.c_str()); }
    }

    inline void approx (double got, double want, double tol, const std::string& msg)
    {
        ++stats().checks;
        if (std::fabs (got - want) > tol)
        {
            ++stats().failures;
            std::printf ("    FAIL: %s (got %.6g, want %.6g, tol %.3g)\n", msg.c_str(), got, want, tol);
        }
    }

    inline void group (const std::string& name) { std::printf ("  - %s\n", name.c_str()); }

    inline int report()
    {
        const auto& s = stats();
        std::printf ("\n%d checks, %d failures\n%s\n", s.checks, s.failures,
                     s.failures == 0 ? "ALL TESTS PASSED" : "TESTS FAILED");
        return s.failures == 0 ? 0 : 1;
    }
}

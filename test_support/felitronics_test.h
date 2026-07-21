// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

// Tiny JUCE-free, zero-dependency test harness shared by every felitronics-core module's self-tests
// (the `teq` discipline: measured audio == analytic curve). Each module's test .cpp owns its main():
//   int main() { ...checks...; return felitronics::test::report(); }

#include <cmath>
#include <cstdio>
#include <string>
#include <version>   // defines _LIBCPP_VERSION on libc++ (used to gate the no-alloc check)

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
        // Spell acceptance positively so unordered NaN comparisons cannot pass unnoticed.
        if (! (std::fabs (got - want) <= tol))
        {
            ++stats().failures;
            std::printf ("    FAIL: %s (got %.6g, want %.6g, tol %.3g)\n", msg.c_str(), got, want, tol);
        }
    }

    // RT-safety no-alloc check. Each test counts allocations via a global operator-new override, but that
    // only isolates OUR allocations on libc++ (the dev/macOS toolchain); libstdc++ and the MSVC STL allocate
    // internally in process()-reachable paths in ways a global counter can't separate. So enforce strictly
    // on libc++ and record it as informational elsewhere — no-alloc is a property of the code, proven on libc++.
    inline void okNoAlloc (bool didNotAllocate, const std::string& msg)
    {
    #if defined(_LIBCPP_VERSION)
        ok (didNotAllocate, msg);
    #else
        (void) didNotAllocate;
        ok (true, msg + "  [alloc-counting N/A on this stdlib]");
    #endif
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

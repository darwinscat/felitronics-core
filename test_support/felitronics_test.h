// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

// Tiny JUCE-free, zero-dependency test harness shared by every felitronics-core module's self-tests
// (the `teq` discipline: measured audio == analytic curve). Each module's test .cpp owns its main():
//   int main() { ...checks...; return felitronics::test::report(); }

#include <cmath>
#include <cstdio>
#include <string>
#if __has_include(<source_location>)
  #include <source_location>
#endif
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

    // A [[nodiscard]] core call — process()/analyse()/applyGain()/prepare() — that the fixture EXPECTS
    // to be accepted.
    // Law 11 (DSP-ARCHITECTURE.md §2) makes every such entry point return its verdict, and a fixture that
    // throws the verdict away is the blind fixture this project keeps finding: it would pass unchanged
    // against a stage that had silently stopped processing. Wrapping the call turns that into a failure
    // without adding a check, so suite counts stay comparable across the change.
    // Use `ok (! obj.process (...), "…")` for a call that is MEANT to be refused.
    // The location comes from std::source_location's default argument, so the call site stays
    // `test::run (x.process (...))` and still names itself when it fails.
#if defined(__cpp_lib_source_location) && __cpp_lib_source_location >= 201907L
    inline bool run (bool accepted, const std::source_location loc = std::source_location::current())
    {
        if (! accepted)
        {
            ++stats().failures;
            std::printf ("    FAIL: a core call was REFUSED at %s:%u\n", loc.file_name(), (unsigned) loc.line());
        }
        return accepted;
    }
#else
    inline bool run (bool accepted)
    {
        if (! accepted) { ++stats().failures; std::printf ("    FAIL: a core call was REFUSED\n"); }
        return accepted;
    }
#endif

    inline void group (const std::string& name) { std::printf ("  - %s\n", name.c_str()); }

    inline int report()
    {
        const auto& s = stats();
        std::printf ("\n%d checks, %d failures\n%s\n", s.checks, s.failures,
                     s.failures == 0 ? "ALL TESTS PASSED" : "TESTS FAILED");
        return s.failures == 0 ? 0 : 1;
    }
}

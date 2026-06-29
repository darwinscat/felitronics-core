// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko <oleh@darwinscat.com> & Alisa <alisa@darwinscat.com>. Part of TabbyEQ — see LICENSE.

#pragma once

// Tiny JUCE-free test harness so the teq module + its tests have ZERO framework deps.

#include <cmath>
#include <cstdio>
#include <string>

namespace teqtest
{
    struct Stats { int checks = 0, failures = 0; };
    inline Stats& stats() { static Stats s; return s; }

    inline void expectTrue (bool ok, const std::string& msg)
    {
        ++stats().checks;
        if (! ok) { ++stats().failures; std::printf ("    FAIL: %s\n", msg.c_str()); }
    }

    inline void expectNear (double got, double want, double tol, const std::string& msg)
    {
        ++stats().checks;
        if (std::fabs (got - want) > tol)
        {
            ++stats().failures;
            std::printf ("    FAIL: %s (got %.6g, want %.6g, tol %.3g)\n", msg.c_str(), got, want, tol);
        }
    }

    inline void group (const std::string& name) { std::printf ("  - %s\n", name.c_str()); }
}

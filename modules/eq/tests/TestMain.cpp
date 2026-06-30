// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// JUCE-free unit-test runner for the teq:: DSP core. Returns non-zero on any failure → CI gate.

#include "TestUtil.h"
#include <cstdio>

void runMatchedBiquadTests();
void runEqEngineTests();

int main()
{
    std::printf ("teq core tests\n");

    runMatchedBiquadTests();
    runEqEngineTests();

    const auto& s = teqtest::stats();
    std::printf ("\n%d checks, %d failures\n", s.checks, s.failures);
    std::printf ("%s\n", s.failures == 0 ? "ALL TESTS PASSED" : "TESTS FAILED");
    return s.failures == 0 ? 0 : 1;
}

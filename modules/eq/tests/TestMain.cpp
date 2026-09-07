// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// JUCE-free unit-test runner for the teq:: DSP core. Returns non-zero on any failure → CI gate.

#include "TestUtil.h"
#include <felitronics_test.h>   // the SHARED harness: felitronics::test::run() reports refused calls there
#include <cstdio>

void runMatchedBiquadTests();
void runEqEngineTests();

int main()
{
    std::printf ("teq core tests\n");

    runMatchedBiquadTests();
    runEqEngineTests();

    // Two counters, one verdict. The migrated teq suite counts in `teqtest`, but a REFUSED process()
    // call is recorded by the SHARED harness — and a refusal that only prints is exactly the silence
    // law 11 exists to remove.
    const auto& s = teqtest::stats();
    const int shared = felitronics::test::stats().failures;
    std::printf ("\n%d checks, %d failures (+%d from the shared harness)\n", s.checks, s.failures, shared);
    std::printf ("%s\n", (s.failures == 0 && shared == 0) ? "ALL TESTS PASSED" : "TESTS FAILED");
    return (s.failures == 0 && shared == 0) ? 0 : 1;
}

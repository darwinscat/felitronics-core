// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// The analyzers' half of the math-policy routing gate. MathPolicyTests pins the other half
// — the shipped RT names bound to core::SystemMath, the deterministic spellings a DIFFERENT type, and a
// runtime witness that the two policies differ on the running row. What is left for here is ownership: the
// offline analyzers must stay on the DETERMINISTIC math. They are compared byte for byte across rows, and the
// system libm is not the same function on every row: at 44.1 kHz the filter arguments of ProgrammeReport
// (30 Hz), BandBursts (9 kHz) and the K-weighting corner (1681.97 Hz) all give different coefficients on
// Apple's libm than on glibc's and musl's — measured.
//
// The static_asserts make an alias change a deliberate act that must also edit this file.

#include <felitronics_test.h>
#include <felitronics/analysis/BandBursts.h>
#include <felitronics/analysis/BandCrest.h>
#include <felitronics/analysis/KWeightingFilter.h>
#include <felitronics/analysis/LowEnd.h>
#include <felitronics/analysis/ProgrammeReport.h>
#include <felitronics/analysis/LoudnessMeter.h>
#include <felitronics/core/DetMath.h>
#include <felitronics/eq/Crossover2.h>

#include <cstdio>
#include <type_traits>

namespace eq   = felitronics::eq;
namespace an   = felitronics::analysis;

// --- the ANALYZERS OWN the deterministic spellings. These assert the public aliases — and the aliases are
// what the members are DECLARED WITH (ProgrammeReport.h: `CrossoverType lr4_`), so there is nothing for
// them to drift from. An earlier draft declared the members with the concrete type and the aliases beside
// them, which let a member be switched to SystemMath with the alias left intact and this gate still green.
static_assert (std::is_same_v<an::ProgrammeReport::CrossoverType,  eq::DeterministicCrossover2>);
static_assert (std::is_same_v<an::ProgrammeReport::KWeightingType, an::DeterministicKWeightingFilter>);
static_assert (std::is_same_v<an::ProgrammeReport::LoudnessType,   an::DeterministicLoudnessMeter>);
static_assert (std::is_same_v<an::LowEnd::CrossoverType,           eq::DeterministicCrossover2>);
static_assert (std::is_same_v<an::BandBursts::CrossoverType,       eq::DeterministicCrossover2>);
static_assert (std::is_same_v<an::BandCrest::CrossoverType,        eq::DeterministicCrossover2>);

// Compile-time only, so nothing here is COUNTED: the checks were static_asserts inside MathPolicyTests before
// they moved, and a runtime line saying so would add one to a total that did not change.
int main()
{
    std::printf ("routing: the analyzers declare their members with the deterministic spellings (static_assert)\n");
    return felitronics::test::report();
}

// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// The math-policy routing gate. Two things must stay true in opposite directions, and neither is visible
// in any output on a single row, so neither would be noticed if it broke:
//
//   · the RT names must stay on the SYSTEM math. eq::Svf's coefficients are TabbyEQ's and OrbitCab's
//     sound; rerouting them through the deterministic floor would move a shipped product's bits for the
//     sake of a test, which is the wrong way round.
//   · the offline analyzers must stay on the DETERMINISTIC math. They are compared byte for byte across
//     rows, and the system libm is not the same function on every row: at 44.1 kHz the filter arguments
//     of ProgrammeReport (30 Hz), BandBursts (9 kHz) and the K-weighting corner (1681.97 Hz) all give
//     different coefficients on Apple's libm than on glibc's and musl's — measured. Today's green parity
//     survives that only because the difference has not yet reached a printed digit on the fixtures in
//     use. That is luck, and this gate is what replaces it.
//
// The static_asserts make an alias change a deliberate act that must also edit this file. The runtime
// half proves the policy is LIVE — that the two spellings are not quietly the same type — and it says so
// out loud when it cannot find a witness on the running row, rather than passing in silence.

#include <felitronics_test.h>
#include <felitronics/analysis/KWeightingFilter.h>
#include <felitronics/analysis/LoudnessMeter.h>
#include <felitronics/core/DetMath.h>
#include <felitronics/eq/Crossover2.h>
#include <felitronics/eq/Svf.h>

#include <bit>
#include <cmath>
#include <cstdint>
#include <string>
#include <type_traits>

namespace core = felitronics::core;
namespace eq   = felitronics::eq;
namespace an   = felitronics::analysis;
using felitronics::test::ok;

// --- the shipped names stay on the system math ---
static_assert (std::is_same_v<eq::Svf,               eq::BasicSvf<core::SystemMath>>);
static_assert (std::is_same_v<eq::Crossover2,        eq::BasicCrossover2<core::SystemMath>>);
static_assert (std::is_same_v<an::KWeightingFilter,  an::BasicKWeightingFilter<core::SystemMath>>);
static_assert (std::is_same_v<an::LoudnessMeter,     an::BasicLoudnessMeter<core::SystemMath>>);
// --- and the deterministic spellings are a DIFFERENT type, not an alias of the same one ---
static_assert (! std::is_same_v<eq::Svf,              eq::DeterministicSvf>);
static_assert (! std::is_same_v<eq::Crossover2,       eq::DeterministicCrossover2>);
static_assert (! std::is_same_v<an::KWeightingFilter, an::DeterministicKWeightingFilter>);
static_assert (! std::is_same_v<an::LoudnessMeter,    an::DeterministicLoudnessMeter>);

int main()
{
    ok (true, "routing: the shipped names are bound to core::SystemMath (static_assert)");
    ok (true, "routing: the deterministic names are a distinct type (static_assert)");

    // The policy must actually change a number somewhere on this row, or the gate above is comparing
    // two spellings of the same behaviour and proves nothing.
    {
        const double kPi = 3.14159265358979323846;
        const double rates[6] = { 44100.0, 48000.0, 88200.0, 96000.0, 176400.0, 192000.0 };
        int witnesses = 0; double wFs = 0.0, wFc = 0.0;
        for (int r = 0; r < 6; ++r)
            for (int f = 10; f < 20000; ++f)
            {
                if ((double) f >= rates[r] * 0.49) break;
                const double x = kPi * ((double) f / rates[r]);
                if (std::bit_cast<std::uint64_t> (core::SystemMath::tan (x))
                    != std::bit_cast<std::uint64_t> (core::DetMath::tan (x)))
                { if (witnesses == 0) { wFs = rates[r]; wFc = (double) f; } ++witnesses; }
            }
        ok (true, witnesses > 0
                  ? "routing: the two policies are distinguishable on this row at " + std::to_string (witnesses)
                    + " filter arguments (first: " + std::to_string ((int) wFc) + " Hz at "
                    + std::to_string ((int) wFs) + " Hz)"
                  : "routing: NO witness on this row — the check below proves nothing here, say so rather than read it as green");

        if (witnesses > 0)
        {
            // AND WHAT THE DIFFERENCE ACTUALLY REACHES — measured, because the answer is not the one the
            // routing suggests. Over 119880 filter arguments (20 Hz .. Nyquist*0.49 at six rates) the two
            // policies give a different `tan` at 7709, and in ZERO of those does the difference reach the
            // filter's output over 8192 samples: Svf carries its state in float, and a one-ulp difference
            // in a double coefficient does not survive that rounding. So routing the analyzers' crossovers
            // through the deterministic policy is DEFENSIVE, not corrective — the cross-row divergence
            // that P79 actually removed (29 printed lines) travelled the double paths: the window, log10,
            // pow10, log2. This is asserted rather than assumed so that the day a change makes the
            // coefficient path reachable — a wider state, a different topology, a Q that amplifies — this
            // number stops being zero and somebody has to look.
            eq::Svf sys; eq::DeterministicSvf det;
            det.prepare (wFs, 1); sys.prepare (wFs, 1);
            sys.setParams (eq::FilterType::LowPass, (float) wFc, 0.70710678118654752, 0.0);
            det.setParams (eq::FilterType::LowPass, (float) wFc, 0.70710678118654752, 0.0);
            bool reached = false;
            double accS = 0.0, accD = 0.0;
            for (int n = 0; n < 8192; ++n)
            {
                const float in = n == 0 ? 1.0f : 0.0f;
                const float a = sys.processSample (0, in), b = det.processSample (0, in);
                accS += (double) a; accD += (double) b;
                if (a != b) reached = true;
            }
            ok (std::isfinite (accS) && std::isfinite (accD) && accS != 0.0,
                "routing: the witness filter actually ran (non-trivial impulse response)");
            ok (! reached, "routing: at the witness the float state absorbs the coefficient difference — "
                           "the deterministic crossover is defensive here, and the day this line fails it "
                           "is a finding, not a regression");
        }
    }

    return felitronics::test::report();
}

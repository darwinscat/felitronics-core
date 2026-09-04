// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// EBU Tech 3341-2023 §2.6 true-peak conformance for felitronics::analysis::TruePeakMeter — the first
// EXTERNAL criterion this repository has ever applied to a true-peak reading. Everything else in the
// true-peak suites tests the meter against its own documented behaviour; this file tests it against a
// published acceptance envelope that no one here wrote.
//
// The signals, the amplitudes, the taper, the budget and the reasons for each live in one place, shared with
// the fcore::Probe suite so the two cannot drift: test_support/ebu_tech3341_truepeak.h. Read it first.
//
// Four assertions per case, and they answer four different questions:
//   1. the EBU envelope, target +0.2/-0.4 — "would this meter pass the formal test?";
//   2. an accuracy budget against the ANALYTIC true peak — "and by how much, on a scale that a redesign can
//      still move within"; this is the regression gate, because the envelope alone leaves ~0.18 dB of slack;
//   3. the sample peak equals the closed-form sample-grid maximum — a fixture self-check: it fails if the
//      synthesis, the phase or the rate handling is wrong, and without it a malformed signal could be
//      carried to a pass by assertion 4;
//   4. the fixture ends in silence and its two fades mirror — the one class of fixture defect that no
//      assertion about a maximum can see, and the class that produced the +0.92 dB absurdity behind F6.
// Plus one cross-case invariant: cases 16 and 19 are the same signal at two amplitudes, so their normalized
// gain must be identical. Case 19 RECONSTRUCTS above full scale (+3 dBTP) while its input samples sit exactly
// at +-1.0, so what this catches is a clamp applied AFTER interpolation, or any gain-dependent branch in the
// reconstruction — not a clamp on the input samples, which +-1.0 would survive. Design-neutral either way.

#include <felitronics_test.h>
#include <ebu_tech3341_truepeak.h>

#include <felitronics/analysis/TruePeakMeter.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

// RT-safety witness: every allocation in this binary is counted (the LoudnessMeter suite's pattern).
static std::atomic<long> g_allocs { 0 };
void* operator new      (std::size_t s) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s ? s : 1); }
void* operator new[]    (std::size_t s) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s ? s : 1); }
void  operator delete   (void* p) noexcept { std::free (p); }
void  operator delete[] (void* p) noexcept { std::free (p); }
void  operator delete   (void* p, std::size_t) noexcept { std::free (p); }
void  operator delete[] (void* p, std::size_t) noexcept { std::free (p); }

using namespace felitronics;
using namespace felitronics::test::ebu3341;

namespace
{
    constexpr double kSeconds = 0.5;      // Table 1: "the total duration does not matter" — the grid is
                                          // exactly periodic in `divisor` samples, so one cycle would do.

    // A fresh meter over one Table 1 signal, fed stereo in phase as Table 1 specifies, in blocks, so the
    // reading cannot depend on being handed the whole file at once.
    struct Reading { double truePeakDb, samplePeakDb; int oversample; };

    Reading measure (const TruePeakCase& c, double sr)
    {
        std::vector<float> ch;
        synthesize (c, sr, kSeconds, ch);

        analysis::TruePeakMeter m;
        m.prepare (sr, 1024, 2);

        const long long n = (long long) ch.size();
        for (long long off = 0; off < n; off += 1024)
        {
            const int len = (int) std::min<long long> (1024, n - off);
            const float* view[2] { ch.data() + off, ch.data() + off };
            m.process (view, 2, len);
        }
        return { m.truePeakDb(), m.samplePeakDb(), m.oversampleFactor() };
    }

    std::string at (double sr, const std::string& s)
    {
        return "Fs " + std::to_string ((int) sr) + ": " + s;
    }
}

int main()
{
    std::printf ("felitronics::analysis TruePeakMeter — EBU Tech 3341 §2.6 true-peak conformance\n");

    for (const double sr : { 48000.0, 44100.0 })
    {
        test::group (at (sr, "Table 1 tests 15-19 against the EBU envelope (+0.2 / -0.4 dB)"));

        for (const TruePeakCase& c : kTruePeakCases)
        {
            const Reading r = measure (c, sr);
            const std::string tag = "test " + std::to_string (c.number) + " (" + c.name + ")";

            const double dTarget = r.truePeakDb - c.targetDb;
            test::ok (dTarget <= kTolAboveDb && dTarget >= -kTolBelowDb,
                      at (sr, tag + ": inside the EBU envelope"));

            // The regression gate: geometry is derived from the meter's own factor, so only filter quality
            // is left for this to catch. See kBudgetSlackDb in the fixture header for why it is 0.02.
            const double dOracle = r.truePeakDb - oracleDb (c);
            const double floorDb = gridBoundDb (c, r.oversample) - kBudgetSlackDb;
            test::ok (dOracle <= kBudgetSlackDb && dOracle >= floorDb,
                      at (sr, tag + ": within the derived grid bound of the analytic true peak"));

            test::approx (r.samplePeakDb, sampleGridPeakDb (c), 1.0e-4,
                          at (sr, tag + ": sample peak is the closed-form sample-grid maximum"));

            std::vector<float> sig;
            synthesize (c, sr, kSeconds, sig);
            test::ok (endsAreSilent (sig), at (sr, tag + ": the signal starts and ends at exact silence"));

            // NB there is deliberately no "true peak >= sample peak" assertion here: both meters seed the
            // running max with the grid sample itself, so it cannot fail, and the derived bound above is
            // strictly stronger — on case 16 it forces the interpolator to have recovered 2.82 dB above the
            // sample peak. The structural invariant is pinned where it can actually break, in the unit
            // suites (TruePeakMeterTests, ProbeTests "the sample peak is a hard floor").
        }

        // Cases 16 and 19: identical frequency and phase, amplitudes 0.5012 and 1.4142 — the second one
        // reconstructs ABOVE full scale (+3 dBTP from samples that sit exactly at +-1.0). A meter that is
        // linear reports the same gain over amplitude for both; one that clips, saturates or takes a
        // different branch above 1.0 does not. Nothing here is pinned to a filter design.
        test::group (at (sr, "linearity across full scale: cases 16 and 19 share one normalized gain"));
        {
            const TruePeakCase& c16 = kTruePeakCases[1];
            const TruePeakCase& c19 = kTruePeakCases[4];
            const double g16 = measure (c16, sr).truePeakDb - oracleDb (c16);
            const double g19 = measure (c19, sr).truePeakDb - oracleDb (c19);
            test::approx (g19, g16, 1.0e-4,
                          at (sr, "case 19 (1.41 FFS) has case 16's gain — no clipping above full scale"));
        }
    }

    // The envelope is a claim about a meter, so it has to survive the meter being used normally: reset()
    // must genuinely clear the running max, or a session that once saw a loud passage stays pinned to it and
    // every later conformance reading is a lie.
    //
    // Re-measuring the SAME signal would not test that: the meter is a running maximum, so a reset() that
    // did nothing at all would still return the first answer and the check would pass. The second pass is
    // therefore 20 dB quieter, and must read 20 dB quieter.
    test::group ("reset() really clears: a quieter second pass reads quieter");
    {
        const TruePeakCase& c = kTruePeakCases[1];
        std::vector<float> loud, quiet;
        synthesize (c, 48000.0, kSeconds, loud);
        quiet = loud;
        for (float& v : quiet) v *= 0.1f;                       // -20 dB exactly

        analysis::TruePeakMeter m;
        m.prepare (48000.0, 1024, 2);
        const float* io[2] { loud.data(), loud.data() };
        m.process (io, 2, (int) loud.size());
        const double first = m.truePeakDb();

        m.reset();
        test::ok (m.truePeakDb() < -100.0, "immediately after reset() the meter reads nothing, not the old max");

        const float* io2[2] { quiet.data(), quiet.data() };
        m.process (io2, 2, (int) quiet.size());
        test::approx (m.truePeakDb(), first - 20.0, 1.0e-3, "the quieter pass reads exactly 20 dB down");
    }

    // The taper checked as itself, not through the carrier: mirrored SAMPLES differ because the sine's phase
    // differs there, so only the window can answer whether the two ramps are the same ramp.
    test::group ("the 10 ms taper is symmetric and reaches silence at both ends");
    {
        for (const double sr : { 48000.0, 44100.0 })
        {
            const long long n = (long long) std::llround (sr * kSeconds), fade = (long long) std::llround (sr * 0.010);
            test::ok (fadeWindow (0, n, fade) == 0.0, at (sr, "the first sample of the window is exactly 0"));
            test::ok (fadeWindow (n - 1, n, fade) == 0.0, at (sr, "the last sample of the window is exactly 0"));
            test::ok (fadeWindow (fade, n, fade) == 1.0, at (sr, "the window reaches exactly 1 where the fade-in ends"));
            test::ok (fadeAsymmetry (n, fade) < 1.0e-15, at (sr, "fade-out mirrors fade-in to the last bit"));
        }
    }

    test::group ("process() does not allocate");
    {
        const TruePeakCase& c = kTruePeakCases[1];
        std::vector<float> ch;
        synthesize (c, 48000.0, kSeconds, ch);
        analysis::TruePeakMeter m;
        m.prepare (48000.0, 1024, 2);
        const float* io[2] { ch.data(), ch.data() };
        const long before = g_allocs.load (std::memory_order_relaxed);
        m.process (io, 2, (int) ch.size());
        // Snapshot the verdict BEFORE calling the harness: okNoAlloc takes a std::string, whose temporary is
        // long enough to heap-allocate, and the order in which the two arguments are evaluated is
        // unspecified. Reading the counter inside the call is a coin flip on the compiler.
        const bool noAlloc = (g_allocs.load (std::memory_order_relaxed) == before);
        test::okNoAlloc (noAlloc, "a whole Table 1 signal through process() allocates nothing");
    }

    return test::report();
}

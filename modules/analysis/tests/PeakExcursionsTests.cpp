// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.
//
// K10 — excursions over a delivery ceiling, on the reconstruction the certificate is issued on.
//
// THE ORACLES HERE ARE ARITHMETIC, NOT A SECOND IMPLEMENTATION. A sine of amplitude A at f over T seconds
// crosses a symmetric threshold exactly 2fT times, its crest's curvature gives its own frequency in closed
// form, and its excess per unit time does not depend on f at all. Each of those is computed in this file
// from the fixture's construction and compared with what the instrument says.

#include <felitronics_test.h>

#include <felitronics/analysis/PeakExcursions.h>
#include <felitronics/analysis/ReferenceTruePeakMeter.h>

#include <cmath>
#include <limits>
#include <vector>

namespace
{
using felitronics::analysis::PeakExcursions;
using felitronics::analysis::ReferenceTruePeakMeter;
namespace test = felitronics::test;
using test::ok;
using test::approx;

constexpr double kPi = 3.14159265358979323846;

std::vector<float> sine (double fs, double seconds, double hz, double dbfs)
{
    const auto n = (std::size_t) (fs * seconds);
    const double a = std::pow (10.0, dbfs / 20.0);
    std::vector<float> v (n);
    for (std::size_t i = 0; i < n; ++i) v[i] = (float) (a * std::sin (2.0 * kPi * hz * (double) i / fs));
    return v;
}

// One call, whole programme. The slicing variants are built where they are needed.
bool feed (PeakExcursions& pe, const std::vector<float>& x)
{
    const float* ch[1] { x.data() };
    if (! pe.process (ch, 1, (int) x.size())) return false;
    return pe.finish();
}
} // namespace

int main()
{
    std::printf ("felitronics::analysis PeakExcursions — K10\n");
    constexpr double kFs = 48000.0;

    // THE INSTRUMENT IS THE CERTIFICATE'S, asserted and not assumed. `PolyphaseOversampler`'s own default
    // is 64 taps per phase; taking it would build a meter that disagrees with the delivery certificate by
    // the instrument gap the true-peak tests pin, and nothing in a passing suite would have said so.
    static_assert (PeakExcursions::kFactor == ReferenceTruePeakMeter::kFactor, "factor must be the reference's");
    static_assert (PeakExcursions::kTapsPerPhase == ReferenceTruePeakMeter::kTapsPerPhase, "taps must be the reference's");

    test::group ("the reconstruction, not the sample grid — the excursion a grid finder cannot see");
    {
        // Two samples of 0.85: the grid says -1.41 dBFS and the reconstruction says about +0.42 dBTP. At a
        // ceiling of -1 dBTP a finder working on the samples reports NOTHING for an excursion of 1.4 dB.
        const int n = 4096;
        std::vector<float> x ((std::size_t) n, 0.0f);
        x[2400] = 0.85f; x[2401] = 0.85f;

        PeakExcursions::Params p; p.thresholdDbtp = -1.0;
        PeakExcursions pe; pe.setParams (p);
        ok (test::run (pe.prepare (kFs, 1)) && test::run (feed (pe, x)), "prepare+feed the pair");

        const double gridPeakDb = 20.0 * std::log10 (0.85);
        approx (gridPeakDb, -1.41, 0.01, "PRECONDITION: the sample peak is " + std::to_string (gridPeakDb) + " dBFS, UNDER the ceiling");
        ok (pe.runCount() == 1, "and the instrument finds exactly one run (got " + std::to_string (pe.runCount()) + ")");
        approx (20.0 * std::log10 (pe.reconstructedPeak()), 0.4224, 0.001,
                "the reconstruction reads " + std::to_string (20.0 * std::log10 (pe.reconstructedPeak())) + " dBTP");

        // AND IT IS THE SAME NUMBER THE CERTIFYING METER REPORTS, bit for bit — the whole point of pinning
        // the topology rather than taking the oversampler's default.
        ReferenceTruePeakMeter m;
        ok (test::run (m.prepare (kFs, n, 1)), "a reference meter over the same audio");
        const float* mc[1] { x.data() };
        ok (test::run (m.process (mc, 1, n)), "…fed");
        m.drain();
        ok (felitronics::core::exactlyEqual (pe.reconstructedPeak(), m.truePeakLinear()),
            "bit-identical to ReferenceTruePeakMeter::truePeakLinear()");
        approx (pe.samplePeakLinear(), (double) 0.85f, 1e-12,
                "and the sample peak is published separately, as the floor it is (the FIXTURE is a float,\n"
                "                 so the oracle is (double) 0.85f and not the literal 0.85 — 2.4e-8 apart)");

        const PeakExcursions::Run r = pe.run (0);
        // The two samples are at 2400 and 2401, so the reconstruction's crest sits between them: the run's
        // start divided by the factor must land inside that interval, which is the coordinate convention
        // working. The residual eighth of an input sample is the FIR's half-sample delay and is real.
        const double at = (double) r.startOs / (double) PeakExcursions::kFactor;
        ok (at > 2400.0 && at < 2401.0, "the run starts at input sample " + std::to_string (at) + ", between the two");
        ok (r.aboveOs == r.lengthOs, "with no merged gap inside it");
    }

    test::group ("a sine: the run count, the crest frequency and the dose are all arithmetic");
    {
        for (const double f : { 40.0, 60.0, 100.0 })
        {
            const auto x = sine (kFs, 1.0, f, -0.01);
            PeakExcursions::Params p; p.thresholdDbtp = -1.0;
            PeakExcursions pe; pe.setParams (p);
            if (! test::run (pe.prepare (kFs, 1)) || ! test::run (feed (pe, x))) continue;
            const std::string at = std::to_string ((int) f) + " Hz: ";

            // A SYMMETRIC THRESHOLD IS CROSSED TWICE PER CYCLE. One second of f Hz has exactly 2f crests
            // over it, and this count is what caught a real defect: the chunk bound was being compared
            // against a sample count updated AFTER the call, so the first chunk stopped early and the
            // answer was 118 where it had to be 120. The aggregate looked healthy — 117 of 118 runs read
            // the right crest frequency — and one run whose crest read twice the tone's was the only tell.
            ok (pe.runCount() == (std::int64_t) (2.0 * f),
                at + "exactly " + std::to_string ((std::int64_t) (2.0 * f)) + " runs, got " + std::to_string (pe.runCount()));

            // THE CREST'S OWN FREQUENCY, from curvature: e = A*w^2*D^2/8 for A*cos(wt), so w = sqrt(8e/(A D^2)).
            // Computed here from the run's published length and peak, and compared with the tone that made it.
            const PeakExcursions::Run r = pe.run (1);          // not run 0: at the programme's edge a crest
            const double d  = (double) r.lengthOs / ((double) PeakExcursions::kFactor * kFs);
            const double e  = r.peak - pe.thresholdLinear();
            const double hz = std::sqrt (8.0 * e / (r.peak * d * d)) / (2.0 * kPi);
            approx (r.crestHz (kFs, pe.thresholdLinear()), hz, 1e-9, at + "the accessor is that arithmetic");

            // AGAINST THE ESTIMATOR'S EXACT THEORY, not against the tone with a loose tolerance. A cosine
            // crest is deeper than the parabola that approximates it, so the estimate reads LOW by exactly
            // sqrt(2e/A) / acos(T/A) — 0.991 at one decibel over. Comparing with the tone at +-5 % would
            // pass a broken estimator; comparing with the bias's closed form pins it to four digits.
            const double A = r.peak, T = pe.thresholdLinear();
            const double bias = std::sqrt (2.0 * (A - T) / A) / std::acos (T / A);
            approx (hz, f * bias, 0.004 * f,
                    at + "reads " + std::to_string (hz) + " Hz where the parabolic estimator must read "
                       + std::to_string (f * bias) + " (bias " + std::to_string (bias) + ")");
            ok (bias < 1.0 && bias > 0.98, at + "…and that bias is under a percent at this depth");

            // …AND THE WHOLE PROGRAMME AGREES, not only one run: every run must land in ONE third-octave,
            // and it is the one holding the ESTIMATE — 40 Hz sits exactly on a bin edge and its estimate
            // lands in the bin below, which is the estimator being right and the tone being on a boundary.
            int want = -1;
            for (int b = 0; b < PeakExcursions::kCrestBins; ++b)
                if (PeakExcursions::crestBinLowHz (b) <= hz
                    && (b + 1 >= PeakExcursions::kCrestBins || PeakExcursions::crestBinLowHz (b + 1) > hz)) want = b;
            ok (want >= 0 && pe.crestBinCount (want) == pe.runCount(),
                at + "every run lands in one third-octave, the one at "
                   + std::to_string (PeakExcursions::crestBinLowHz (want)) + " Hz");
            approx (pe.crestBinDose (want) / pe.totalDose(), 1.0, 1e-12, at + "carrying all of the dose");
        }

        // DOSE PER SECOND DOES NOT DEPEND ON FREQUENCY, and that is the property that makes it comparable
        // between programmes: it is a time-average over the waveform's amplitude distribution, which a sine
        // has independently of how fast it traverses it. If dose were "peak times duration once per run" —
        // the other candidate — this would be false, because the count of runs scales with f.
        double d40 = 0.0, d100 = 0.0;
        for (const double f : { 40.0, 100.0 })
        {
            const auto x = sine (kFs, 1.0, f, -0.01);
            PeakExcursions::Params p; p.thresholdDbtp = -1.0;
            PeakExcursions pe; pe.setParams (p);
            if (! test::run (pe.prepare (kFs, 1)) || ! test::run (feed (pe, x))) continue;
            (f < 50.0 ? d40 : d100) = pe.totalDose();
        }
        ok (d40 > 0.0 && d100 > 0.0, "PRECONDITION: both measured a dose");
        approx (d100 / d40, 1.0, 2e-3, "40 Hz and 100 Hz carry the same dose per second ("
                                       + std::to_string (d40) + " against " + std::to_string (d100) + ")");
    }

    test::group ("merging moves the count and the percentile, and moves NOTHING else");
    {
        const auto x = sine (kFs, 1.0, 60.0, -0.01);
        PeakExcursions::Params zero, wide;
        zero.thresholdDbtp = wide.thresholdDbtp = -1.0;
        zero.mergeMs = 0.0;
        wide.mergeMs = 6.0;                       // over 1/(2*60) = 8.3 ms? no — under it, but well over the gaps
        PeakExcursions a, b;
        a.setParams (zero); b.setParams (wide);
        if (test::run (a.prepare (kFs, 1)) && test::run (feed (a, x))
            && test::run (b.prepare (kFs, 1)) && test::run (feed (b, x)))
        {
            ok (b.runCount() < a.runCount(), "a wider window fuses runs (" + std::to_string (a.runCount())
                                             + " against " + std::to_string (b.runCount()) + ")");
            ok (b.p90Ms() > a.p90Ms(), "and lengthens the percentile");
            // THE TWO FIXED POINTS. Both are accumulated over raw above-samples, so a consumer sweeping the
            // window keeps them still while everything else moves — which is what makes the sweep readable.
            ok (a.aboveOs() == b.aboveOs(), "the above-sample count is merge-free");
            ok (felitronics::core::exactlyEqual (a.occupancy(), b.occupancy()), "so occupancy is bit-identical");
            ok (felitronics::core::exactlyEqual (a.totalDose(), b.totalDose()), "and so is the total dose");
        }
    }

    test::group ("nothing over the ceiling: the trivially-satisfied rule, published so it cannot be read as one");
    {
        const auto x = sine (kFs, 0.5, 100.0, -12.0);          // nowhere near -1 dBTP
        PeakExcursions::Params p; p.thresholdDbtp = -1.0;
        PeakExcursions pe; pe.setParams (p);
        if (test::run (pe.prepare (kFs, 1)) && test::run (feed (pe, x)))
        {
            ok (pe.runCount() == 0, "no runs");
            ok (pe.aboveOs() == 0 && pe.totalDose() == 0.0, "nothing above, no dose");
            ok (pe.occupancy() == 0.0 && pe.runsPerMinute() == 0.0, "occupancy and rate read their canonical zero");
            // THE TRAP, STATED: a rule of the form "p90 <= 2 ms" is TRUE here, of a programme that never
            // went over at all. The field that separates "nothing happened" from "things happened briefly"
            // is the count, and it is published beside the percentile for exactly that reason.
            ok (pe.p90Ms() == 0.0, "the percentile is 0.0 — which satisfies 'p90 <= 2 ms' while meaning the opposite");
            ok (pe.runCount() == 0, "…and runCount() is what tells the two apart");
            ok (pe.valid(), "the measurement is valid: nothing over a ceiling is an answer, not a failure");
        }
    }

    test::group ("law 8a — the clock is the sample, never the call");
    {
        const auto x = sine (kFs, 0.4, 60.0, -0.01);
        PeakExcursions::Params p; p.thresholdDbtp = -1.0;
        PeakExcursions one, cut;
        one.setParams (p); cut.setParams (p);
        if (test::run (one.prepare (kFs, 1)) && test::run (feed (one, x)) && test::run (cut.prepare (kFs, 1)))
        {
            const int cuts[] = { 1, 4093, 17, 1024, 333, 2 };
            std::size_t at = 0; int k = 0; bool okAll = true;
            while (at < x.size())
            {
                const int m = (int) std::min<std::size_t> ((std::size_t) cuts[k++ % 6], x.size() - at);
                const float* ch[1] { x.data() + at };
                if (! cut.process (ch, 1, m)) { okAll = false; break; }
                at += (std::size_t) m;
            }
            ok (okAll && test::run (cut.finish()), "the same programme in irregular calls");
            ok (cut.runCount() == one.runCount() && cut.aboveOs() == one.aboveOs(), "the same counts");
            ok (felitronics::core::exactlyEqual (cut.totalDose(), one.totalDose()), "the same dose, bit for bit");
            ok (felitronics::core::exactlyEqual (cut.reconstructedPeak(), one.reconstructedPeak()), "the same peak");
            bool same = true;
            for (std::int64_t i = 0; i < one.storedRunCount(); ++i)
                if (cut.run (i).startOs != one.run (i).startOs || cut.run (i).lengthOs != one.run (i).lengthOs) same = false;
            ok (same, "and every run at the same coordinate");
        }
    }

    test::group ("law 11 — the list fills, every statistic does not");
    {
        const auto x = sine (kFs, 1.0, 60.0, -0.01);
        PeakExcursions::Params big, small;
        big.thresholdDbtp = small.thresholdDbtp = -1.0;
        small.maxRuns = 7;                                   // far under the 120 this programme makes
        PeakExcursions a, b;
        a.setParams (big); b.setParams (small);
        if (test::run (a.prepare (kFs, 1)) && test::run (feed (a, x))
            && test::run (b.prepare (kFs, 1)) && test::run (feed (b, x)))
        {
            ok (b.runCount() == a.runCount(), "the count keeps counting past the capacity ("
                                              + std::to_string (b.runCount()) + ")");
            ok (b.storedRunCount() == 7 && ! b.runsComplete(), "the list holds its prefix and says it is short");
            ok (a.runsComplete(), "while the roomy one is complete");
            // EVERY AGGREGATE IS AN ACCUMULATOR AND NEVER READS THE LIST, which is what keeps these exact.
            ok (felitronics::core::exactlyEqual (a.p90Ms(), b.p90Ms()), "the percentile is unchanged");
            ok (felitronics::core::exactlyEqual (a.occupancy(), b.occupancy()), "occupancy is unchanged");
            ok (felitronics::core::exactlyEqual (a.totalDose(), b.totalDose()), "the dose is unchanged");
            bool sameCls = true, sameCrest = true;
            for (int k = 0; k < PeakExcursions::kClasses; ++k)
                if (a.classCount (k) != b.classCount (k)) sameCls = false;
            for (int k = 0; k < PeakExcursions::kCrestBins; ++k)
                if (a.crestBinCount (k) != b.crestBinCount (k)) sameCrest = false;
            ok (sameCls, "every duration class is unchanged");
            ok (sameCrest, "and so is the crest-frequency histogram");
        }
    }

    test::group ("refusals — exactly the arguments prepare() will not take");
    {
        PeakExcursions::Params p;
        PeakExcursions::Storage st;
        ok (! PeakExcursions::storageFor (0.0, 1, p).ok, "a rate of zero");
        ok (! PeakExcursions::storageFor (std::numeric_limits<double>::quiet_NaN(), 1, p).ok, "a NaN rate");
        ok (! PeakExcursions::storageFor (1.0e9, 1, p).ok, "a rate past the ceiling");
        ok (! PeakExcursions::storageFor (kFs, 0, p).ok, "zero channels");
        ok (! PeakExcursions::storageFor (kFs, felitronics::core::kMaxChannels + 1, p).ok, "one channel too many");
        { auto q = p; q.thresholdDbtp = std::numeric_limits<double>::quiet_NaN();
          ok (! PeakExcursions::storageFor (kFs, 1, q).ok, "a NaN ceiling"); }
        { auto q = p; q.thresholdDbtp = -1.0e9;
          ok (! PeakExcursions::storageFor (kFs, 1, q).ok, "a ceiling past what det::pow10 keeps finite"); }
        { auto q = p; q.mergeMs = -1.0;
          ok (! PeakExcursions::storageFor (kFs, 1, q).ok, "a negative merge window"); }
        { auto q = p; q.mergeMs = std::numeric_limits<double>::infinity();
          ok (! PeakExcursions::storageFor (kFs, 1, q).ok, "an infinite one — refused BEFORE llround sees it"); }
        { auto q = p; q.classEdgesMs[2] = q.classEdgesMs[1];
          ok (! PeakExcursions::storageFor (kFs, 1, q).ok, "class edges that do not strictly increase"); }
        { auto q = p; q.classEdgesMs[0] = 0.0;
          ok (! PeakExcursions::storageFor (kFs, 1, q).ok, "a first edge of zero length"); }
        { auto q = p; q.maxRuns = -1;
          ok (! PeakExcursions::storageFor (kFs, 1, q).ok, "a negative capacity"); }
        ok (PeakExcursions::storageFor (kFs, 2, p).ok && PeakExcursions::storageFor (kFs, 2, p).bytes() > 0,
            "…while the defaults are accepted and priced");

        // Law 11b: a refused prepare leaves an object that cannot be read as a measurement.
        PeakExcursions bad;
        auto q = p; q.mergeMs = -1.0; bad.setParams (q);
        ok (! bad.prepare (kFs, 1), "prepare refuses what storageFor refuses");
        ok (! bad.isPrepared() && bad.reason() == PeakExcursions::Reason::NotFinished,
            "and the object is disarmed, not left holding a previous answer");

        // Law 11a: a narrower call than the first is refused rather than guessed at — the union at an
        // oversampled index needs every channel at that index.
        PeakExcursions two; two.setParams (p);
        if (test::run (two.prepare (kFs, 2)))
        {
            std::vector<float> l (2048, 0.1f), r2 (2048, 0.1f);
            const float* both[2] { l.data(), r2.data() };
            ok (two.process (both, 2, 2048), "two channels are accepted");
            const float* one[1] { l.data() };
            ok (! two.process (one, 1, 2048), "and a narrower call after them is REFUSED");
        }
    }

    test::group ("non-finite input is counted and located, and the measurement says so");
    {
        auto x = sine (kFs, 0.2, 60.0, -0.01);
        x[5000] = std::numeric_limits<float>::quiet_NaN();
        x[5001] = std::numeric_limits<float>::infinity();
        PeakExcursions::Params p; p.thresholdDbtp = -1.0;
        PeakExcursions pe; pe.setParams (p);
        if (test::run (pe.prepare (kFs, 1)) && test::run (feed (pe, x)))
        {
            ok (pe.nonFiniteSamples() == 2, "both were counted (got " + std::to_string (pe.nonFiniteSamples()) + ")");
            ok (pe.firstNonFiniteAt() == 5000, "and the first is located at " + std::to_string (pe.firstNonFiniteAt()));
            ok (pe.reason() == PeakExcursions::Reason::NonFiniteInput && ! pe.valid(),
                "the report refuses to call itself valid");
            ok (pe.runCount() > 0, "…while still publishing what it measured, which is the useful half");
        }
    }

    return test::report();
}

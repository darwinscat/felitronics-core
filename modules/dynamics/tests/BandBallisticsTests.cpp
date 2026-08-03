// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// JUCE-free self-tests for felitronics::dynamics::BandBallistics. Property tests against the ANALYTIC
// ring time of a 2nd-order bandpass — plus a reference NULL against a directly measured eq::Svf probe
// envelope, which is what falsifies the law rather than mirroring it.

#include <felitronics_test.h>
#include <felitronics/dynamics/BandBallistics.h>
#include <felitronics/eq/Svf.h>
#include <felitronics/dynamics/GainReductionFollower.h>
#include <felitronics/core/Math.h>

#include <cmath>
#include <cstdio>
#include <limits>
#include <vector>

using namespace felitronics;

// Measure the real 10-90% envelope rise of an SVF bandpass driven at its centre, tracked with a
// follower fast enough not to be what we are measuring. This is the ground truth the law must match.
static double measuredRiseMs (double fs, double fc, double Q)
{
    eq::Svf bp; bp.prepare (fs, 1);
    bp.setParams (eq::FilterType::BandPass, fc, Q, 0.0);

    const int n = (int) (fs * 3.0);
    std::vector<double> env ((size_t) n, 0.0);
    const double c = std::exp (-1.0 / (0.00005 * fs));      // 0.05 ms — well under any band here
    double e = 0.0;
    for (int i = 0; i < n; ++i)
    {
        const double x = std::sin (2.0 * core::kPi * fc * (double) i / fs);
        const double y = std::fabs ((double) bp.processSample (0, (float) x));
        e = y + c * (e - y);
        env[(size_t) i] = e;
    }
    double peak = 0.0;
    for (int i = n / 2; i < n; ++i) peak = std::fmax (peak, env[(size_t) i]);
    if (peak <= 0.0) return -1.0;

    int i10 = -1, i90 = -1;
    for (int i = 0; i < n; ++i)
    {
        if (i10 < 0 && env[(size_t) i] >= 0.1 * peak) i10 = i;
        if (i90 < 0 && env[(size_t) i] >= 0.9 * peak) { i90 = i; break; }
    }
    if (i10 < 0 || i90 < 0) return -1.0;
    return 1.0e3 * (double) (i90 - i10) / fs;
}

int main()
{
    std::printf ("felitronics::dynamics::BandBallistics tests\n");
    const double fs = 48000.0;
    using BB = dynamics::BandBallistics;

    // The whole point of the law: ring time is LINEAR in Q/fc (inverse bandwidth), not in the period.
    test::group ("ringMs == analytic ln9*Q/(pi*fc) well below Nyquist");
    {
        const double ln9 = 2.1972245773362196;
        for (double fc : { 100.0, 250.0, 1000.0 })
            for (double Q : { 3.5, 10.0, 40.0 })
            {
                const double analytic = 1.0e3 * ln9 * Q / (core::kPi * fc);
                test::approx (BB::ringMs (fs, fc, Q), analytic, analytic * 0.02,
                              "ring matches analytic form at fc/fs << 1");
            }
    }

    // A reference NULL against the actual filter — this is what a wrong law fails.
    test::group ("ringMs NULLs against a measured eq::Svf probe (incl. near Nyquist)");
    {
        struct Case { double fc, Q, tol; };
        for (const Case c : { Case { 1000.0, 10.0, 0.06 }, Case { 1000.0, 40.0, 0.06 },
                              Case {  100.0, 40.0, 0.06 }, Case { 7000.0, 40.0, 0.12 },
                              Case { 12000.0, 40.0, 0.15 } })
        {
            const double meas = measuredRiseMs (fs, c.fc, c.Q);
            const double law  = BB::ringMs (fs, c.fc, c.Q);
            test::ok (meas > 0.0, "probe rise measurable");
            test::approx (law, meas, meas * c.tol, "law matches the measured probe rise");
        }
    }

    // The analog form is not merely less accurate near Nyquist — it is wrong by a factor that grows.
    test::group ("digital-bandwidth form beats the analog one near Nyquist");
    {
        const double ln9 = 2.1972245773362196;
        const double fc = 16000.0, Q = 40.0;
        const double analog = 1.0e3 * ln9 * Q / (core::kPi * fc);
        const double meas   = measuredRiseMs (fs, fc, Q);
        const double law    = BB::ringMs (fs, fc, Q);
        test::ok (meas > 1.8 * analog, "analog form under-predicts by >1.8x at fc=16k/Q=40");
        test::ok (std::fabs (law - meas) < std::fabs (analog - meas), "digital form is closer");
    }

    // UNITS. ringMs is a 10-90% rise; GainReductionFollower reads ms as tau. compute() converts, so
    // the follower's OWN measured 10-90 must come out at attackFraction x the band's ring time. A
    // build that skips the ln9 conversion lands 2.197x slow and fails here — and nowhere else.
    test::group ("compute() composed with the real follower honours attackFraction");
    {
        // Q must be high enough that the RING term wins over the period floor (crossover is ~Q 23),
        // otherwise this would be measuring the floor and comparing it against the law.
        for (double fc : { 200.0, 1000.0, 4000.0 })
            for (double Q : { 30.0, 40.0 })
            {
                const double atkMs  = BB::compute (fs, fc, Q).attackMs;
                const double period = 2.5 * 1.0e3 / fc;
                if (atkMs <= 1.0 + 1.0e-9) continue;             // parked on the auto floor, not the law
                test::ok (atkMs > period, "the ring term, not the period floor, is what set this attack");

                dynamics::GainReductionFollower gr; gr.prepare (fs); gr.setTimes (atkMs, 200.0);
                const double target = -12.0;
                int i10 = -1, i90 = -1;
                for (int i = 0; i < (int) (fs * 2.0); ++i)
                {
                    const double v = gr.process ((float) target);
                    if (i10 < 0 && v <= 0.1 * target) i10 = i;
                    if (i90 < 0 && v <= 0.9 * target) { i90 = i; break; }
                }
                test::ok (i10 >= 0 && i90 > i10, "follower rise measurable");
                const double measured10to90 = 1.0e3 * (double) (i90 - i10) / fs;
                const double expected = 0.35 * BB::ringMs (fs, fc, Q);   // attackFraction x band ring
                test::approx (measured10to90, expected, expected * 0.12,
                              "follower reaches 10-90 in attackFraction x the band's own ring time");
            }
    }

    test::group ("monotonicity: attack falls with fc, rises with Q");
    {
        double prev = std::numeric_limits<double>::infinity();
        for (double fc : { 60.0, 120.0, 500.0, 2000.0, 8000.0 })
        {
            const double a = BB::compute (fs, fc, 8.0).attackMs;
            test::ok (a <= prev + 1.0e-9, "attack non-increasing in fc");
            prev = a;
        }
        prev = 0.0;
        for (double Q : { 0.5, 1.0, 4.0, 16.0, 40.0 })
        {
            const double a = BB::compute (fs, 500.0, Q).attackMs;
            test::ok (a >= prev - 1.0e-9, "attack non-decreasing in Q");
            prev = a;
        }
    }

    // Below Q = pi/ln9 ~ 1.43 the ring is shorter than a cycle; a sub-period attack would modulate
    // gain inside the waveform. The floor must bind there.
    test::group ("period floor binds at low Q");
    {
        const double fc = 200.0, Q = 0.5;
        const double a = BB::compute (fs, fc, Q).attackMs;
        test::ok (a >= 2.0 * 1.0e3 / fc - 1.0e-9, "attack at least ~2 periods at low Q");
        test::ok (a > BB::ringMs (fs, fc, Q) * 0.35, "floor, not the ring term, is what set it");
    }

    test::group ("rails hold, and the release hierarchy is respected");
    {
        const auto slow = BB::compute (fs, 100.0, 40.0);     // physics wants ~280 ms attack
        test::ok (slow.attackMs <= 300.0 + 1.0e-9, "attack ceiling respected");
        test::ok (slow.releaseMs <= 500.0 + 1.0e-9, "release ceiling respected — never seconds");
        const auto fast = BB::compute (fs, 12000.0, 0.7);
        test::ok (fast.attackMs >= 1.0 - 1.0e-9, "attack floor respected");
        test::ok (fast.releaseMs >= 12.0 - 1.0e-9, "release floor respected");
        for (double fc : { 80.0, 700.0, 6000.0 })
            for (double Q : { 0.7, 5.0, 30.0 })
            {
                const auto t = BB::compute (fs, fc, Q);
                test::ok (t.releaseMs >= t.attackMs, "release never shorter than attack");
            }
    }

    test::group ("deviation knob: x1 at centre, x1/4 and x4 at the ends");
    {
        test::approx (BB::knobMultiplier (0.5), 1.0, 1.0e-12, "centre is exactly auto");
        test::approx (BB::knobMultiplier (0.0), 0.25, 1.0e-12, "0 -> four times faster");
        test::approx (BB::knobMultiplier (1.0), 4.0, 1.0e-12, "1 -> four times slower");
        // Interior points too: three endpoints alone are satisfied by a straight line, which would
        // give completely the wrong feel everywhere between them.
        test::approx (BB::knobMultiplier (0.25), 0.5, 1.0e-12, "quarter travel is exactly half speed");
        test::approx (BB::knobMultiplier (0.75), 2.0, 1.0e-12, "three-quarter travel is exactly double");
        // Applied before the clamp, so a mid-band (off the rails) shows the multiplier exactly.
        const double autoA = BB::compute (fs, 500.0, 8.0, 0.5, 0.5).attackMs;
        const double slowA = BB::compute (fs, 500.0, 8.0, 1.0, 0.5).attackMs;
        test::approx (slowA / autoA, 4.0, 0.02, "knob scales attack off the rails");
    }

    // A release floor set too high pins every treble band to the same value and turns the release
    // knob into a placebo for exactly the de-esser case this module exists to serve.
    test::group ("the release knob is alive where a de-esser lives");
    {
        for (double fc : { 5000.0, 7000.0, 12000.0 })
            for (double Q : { 3.5, 10.0 })
            {
                const double fastRel = BB::compute (fs, fc, Q, 0.5, 0.0).releaseMs;
                const double slowRel = BB::compute (fs, fc, Q, 0.5, 1.0).releaseMs;
                test::ok (slowRel > fastRel * 1.5, "release knob still moves the release up here");
            }
    }

    test::group ("sample-rate invariance of the law");
    {
        for (double rate : { 44100.0, 96000.0 })
            for (double fc : { 500.0, 4000.0 })
            {
                const double meas = measuredRiseMs (rate, fc, 20.0);
                const double law  = BB::ringMs (rate, fc, 20.0);
                test::approx (law, meas, meas * 0.08, "law NULLs against the probe at this rate too");
            }
    }

    test::group ("a hostile Config cannot produce UB or escape sanity");
    {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        BB::Config bad;
        bad.attackMinMs  = 500.0; bad.attackMaxMs  = 2.0;      // inverted rails: std::clamp UB
        bad.releaseMinMs = 900.0; bad.releaseMaxMs = 10.0;
        bad.periodFloor  = nan;   bad.attackFraction = -1.0;
        bad.releaseFactor = nan;  bad.knobRange = nan;
        const auto t = BB::compute (fs, 1000.0, 4.0, 0.5, 0.5, bad);
        test::ok (std::isfinite (t.attackMs) && std::isfinite (t.releaseMs), "finite despite the rails");
        test::ok (t.releaseMs >= t.attackMs, "ordering invariant survives");
    }

    test::group ("release is never shorter than attack, even with opposing knobs");
    {
        for (double fc : { 80.0, 100.0, 1000.0 })
            for (double Q : { 10.0, 40.0 })
            {
                const auto t = BB::compute (fs, fc, Q, 0.0, 0.0);   // attack fast-ish, release fastest
                test::ok (t.releaseMs >= t.attackMs - 1.0e-9, "ordering holds at opposing knob extremes");
            }
    }

    test::group ("hostile inputs stay finite and inside the rails");
    {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        const double inf = std::numeric_limits<double>::infinity();
        for (double fc : { nan, inf, -inf, -100.0, 0.0, 1.0e9 })
            for (double Q : { nan, inf, -5.0, 0.0, 1.0e6 })
            {
                const auto t = BB::compute (fs, fc, Q, nan, inf);
                test::ok (std::isfinite (t.attackMs) && std::isfinite (t.releaseMs), "finite output");
                test::ok (t.attackMs >= 1.0 && t.attackMs <= 300.0, "attack inside rails");
                test::ok (t.releaseMs >= 12.0 && t.releaseMs <= 500.0, "release inside rails");
            }
        const auto z = BB::compute (0.0, 1000.0, 2.0);
        test::ok (std::isfinite (z.attackMs), "zero sample rate falls back rather than dividing by it");
    }

    return test::report();
}

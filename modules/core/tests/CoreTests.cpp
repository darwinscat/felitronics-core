// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// JUCE-free self-tests for felitronics::core (Smoother, Math, denormal flush).

#include <felitronics_test.h>
#include <felitronics/core/Config.h>
#include <felitronics/core/Math.h>
#include <felitronics/core/Smoother.h>
#include <felitronics/core/FlushToZero.h>
#include <felitronics/core/DelayLine.h>

#include <cmath>
#include <cstdio>
#include <limits>

using namespace felitronics;

int main()
{
    std::printf ("felitronics::core tests\n");

    // --- Config / sizes ---
    test::group ("Config");
    test::ok (core::kMaxChannels == 16, "kMaxChannels SSOT == 16");
    test::ok (sizeof (core::Sample) == sizeof (float), "Sample is float (Law 3)");

    // --- Smoother: snap is immediate; long run settles to target ---
    test::group ("Smoother snap + settle");
    {
        core::Smoother s; s.prepare (48000.0, 10.0); s.snap (1.0);
        test::approx (s.value(), 1.0, 1e-12, "snap sets value immediately");
        test::ok (s.settled(), "snapped == settled");
        s.setTarget (0.0);
        for (int i = 0; i < 48000; ++i) s.next();          // 1 s >> 10 ms time constant
        test::ok (s.settled(), "settled after 1 s");
        test::approx (s.value(), 0.0, 1e-6, "converged to target");
    }

    // --- advance(n) closed form == n stepwise next() ---
    test::group ("Smoother advance(n) == n*next()");
    {
        core::Smoother a, b; a.prepare (48000.0, 30.0); b.prepare (48000.0, 30.0);
        a.snap (0.0); b.snap (0.0); a.setTarget (1.0); b.setTarget (1.0);
        for (int i = 0; i < 100; ++i) a.next();
        b.advance (100);
        test::approx (b.value(), a.value(), 1e-9, "closed-form advance matches stepwise");
    }

    // --- LinearSmoother: bit-for-bit drop-in for juce::SmoothedValue<float, Linear> ---
    // (assertions ported from JUCE's own CommonSmoothedValueTests, so a swap can't change the feel)
    test::group ("LinearSmoother == juce::SmoothedValue (Linear)");
    {
        // initial state
        core::LinearSmoother sv;
        test::approx (sv.getCurrentValue(), sv.getTargetValue(), 0.0, "init: current == target");
        sv.getNextValue();
        test::ok (! sv.isSmoothing(), "init: not smoothing");

        // resetting + ramp arming
        core::LinearSmoother r (15.0f);
        r.reset (3);
        test::approx (r.getCurrentValue(), 15.0, 1e-6, "reset keeps current at 15");
        r.setTargetValue (16.0f);
        test::ok (r.getTargetValue() == 16.0f && r.getCurrentValue() == 15.0f && r.isSmoothing(),
                  "target armed, current not moved yet, smoothing");
        test::ok (r.getNextValue() > 15.0f, "first step moves toward target");
        r.reset (5);
        test::ok (r.getCurrentValue() == 16.0f && ! r.isSmoothing(), "reset snaps current to target + stops");
        r.setCurrentAndTargetValue (0.2f);
        test::approx (r.getNextValue(), 0.2, 1e-6, "snapped value held");
        test::ok (! r.isSmoothing(), "snap: not smoothing");

        // linear ramp SHAPE: 1 -> 2 over 12 steps, snap exactly on the 12th
        core::LinearSmoother b (1.0f);
        b.reset (12); b.setTargetValue (2.0f);
        float s[15]; for (int i = 0; i < 15; ++i) s[i] = b.getNextValue();
        test::approx (s[0],  1.0 + 1.0 / 12.0,  1e-6, "step 1 == 1 + 1/12 (fixed increment)");
        test::approx (s[10], 1.0 + 11.0 / 12.0, 1e-5, "step 11 == 1 + 11/12");
        test::ok     (s[10] < 2.0f, "still below target before the final step");
        test::approx (s[11], 2.0, 2e-7, "reaches exactly target on step 12");
        test::approx (s[14], 2.0, 0.0, "holds target after");

        // skip(n) == n * getNextValue()
        core::LinearSmoother g; g.reset (12); g.setCurrentAndTargetValue (1.0f); g.setTargetValue (2.0f);
        float ref[15]; for (int i = 0; i < 15; ++i) ref[i] = g.getNextValue();
        g.setCurrentAndTargetValue (1.0f); g.setTargetValue (2.0f);
        test::approx (g.skip (1), ref[0], 1e-5, "skip 1 == getNextValue #0");
        test::approx (g.skip (1), ref[1], 1e-5, "skip 1 == getNextValue #1");
        test::approx (g.skip (2), ref[3], 1e-5, "skip 2 == getNextValue #3");
        g.skip (3);
        test::approx (g.getCurrentValue(), ref[6], 1e-5, "after skip 3, current == #6");
        test::approx (g.skip (300), g.getTargetValue(), 0.0, "skip past the end snaps to target");

        // reset(sampleRate, seconds) uses floor(), matching JUCE
        core::LinearSmoother f; f.reset (1000.0, 0.0125);            // floor(12.5) == 12 steps
        f.setCurrentAndTargetValue (0.0f); f.setTargetValue (12.0f); // step == 1.0/sample
        test::approx (f.getNextValue(), 1.0, 1e-6, "floor(0.0125*1000)=12 steps -> +1.0/step");

        // juce::approximatelyEqual tolerance: a ~1-ULP target change is a NO-OP (no ramp restart / zipper),
        // exactly like juce::SmoothedValue::setTargetValue — automation that jitters the last bit stays quiet.
        core::LinearSmoother t (1.0f); t.reset (100);
        t.setTargetValue (std::nextafter (1.0f, 2.0f));             // one ULP up == exactly FLT_EPSILON away
        test::ok (! t.isSmoothing(), "1-ULP target change is ignored (juce tolerance) — no spurious ramp");
        t.setTargetValue (1.001f);                                 // a genuine change (>> epsilon*|v|)
        test::ok (t.isSmoothing(), "a real target change still arms the ramp");
    }

    // --- Math dB <-> gain ---
    test::group ("Math dB<->gain");
    test::approx (core::dbToGain (0.0),    1.0,    1e-12, "0 dB == unity");
    test::approx (core::dbToGain (6.0206), 2.0,    1e-3,  "+6.0206 dB ~ 2x");
    test::approx (core::gainToDb (2.0),    6.0206, 1e-3,  "2x ~ +6.0206 dB");

    // --- Software denormal flush (Law 8, portable) ---
    test::group ("flushDenormal");
    {
        float f = 1.0e-20f; core::flushDenormal (f); test::ok (f == 0.0f, "tiny float -> exact 0");
        float g = 0.5f;     core::flushDenormal (g); test::ok (g == 0.5f, "normal float kept");
        double d = 1.0e-20; core::flushDenormal (d); test::ok (d == 0.0,  "tiny double -> exact 0");
        double e = 0.5;     core::flushDenormal (e); test::ok (e == 0.5,  "normal double kept");
    }

    // --- ScopedFlushToZero: constructs/restores without crashing on any tier (smoke) ---
    test::group ("ScopedFlushToZero smoke");
    {
        { core::ScopedFlushToZero ftz; volatile float x = 1.0e-30f; x = x * 0.5f; (void) x; }
        test::ok (true, "scoped FTZ constructs + destructs (no-op on tiers without an FP-control reg)");
    }

    // --- DelayLine: passthrough at 0, exact integer delay, capacity edge ---
    test::group ("DelayLine");
    {
        core::DelayLine dl; dl.prepare (64);
        dl.setDelay (0);
        test::approx (dl.process (1.0f), 1.0, 1e-9, "delay 0 == passthrough");

        dl.reset (); dl.setDelay (5);
        float out[20];
        for (int i = 0; i < 20; ++i) out[i] = dl.process (i == 0 ? 1.0f : 0.0f);
        test::approx (out[5], 1.0, 1e-9, "impulse appears at exactly delay 5");
        test::ok (out[0] == 0.0f && out[4] == 0.0f, "nothing before the delay");

        dl.reset (); dl.setDelay (64);            // capacity edge
        float last = 0.0f;
        for (int i = 0; i < 65; ++i) last = dl.process (i == 0 ? 1.0f : 0.0f);
        test::approx (last, 1.0, 1e-9, "impulse appears at delay == capacity (64)");
    }

    // --- DelayLine: a DEFAULT (unprepared) line is a valid delay-0 passthrough, not an empty-buffer OOB ---
    test::group ("DelayLine: default (unprepared) is a safe passthrough");
    {
        core::DelayLine raw;                          // never prepared — buf_ defaults to one slot
        test::approx (raw.process (0.5f), 0.5, 1e-9, "unprepared process() = delay-0 passthrough (no empty-buffer OOB)");
        test::ok (raw.capacity() == 0, "unprepared capacity is 0");
    }

    // --- FALSIFICATION: a non-finite smoothing time must not poison the smoother state ---
    test::group ("Smoother non-finite time rejected");
    {
        core::Smoother s; s.prepare (48000.0, 30.0);
        s.setTimeMs (std::numeric_limits<double>::quiet_NaN());   // NaN tau → coeff must fall back, not go NaN
        s.snap (0.0); s.setTarget (1.0);
        const double v = s.next();
        test::ok (std::isfinite (v), "next() stays finite after setTimeMs(NaN)");
    }

    return test::report();
}

// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa. Part of felitronics-core — see LICENSE.

// JUCE-free self-tests for felitronics::core (Smoother, Math, denormal flush).

#include <felitronics_test.h>
#include <felitronics/core/Config.h>
#include <felitronics/core/Math.h>
#include <felitronics/core/Smoother.h>
#include <felitronics/core/FlushToZero.h>
#include <felitronics/core/DelayLine.h>

#include <cstdio>

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

    return test::report();
}

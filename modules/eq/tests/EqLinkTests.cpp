// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// felitronics::eq::EqLink (moved from OrbitAmp with its gate) — the console-grammar link: shelves with free
// corners, three bells, cut filters with real slopes. Plus the invariant that
// matters most: the curve the UI draws is the filter chain that actually runs.
//
// JUCE-free, like the code under test: links only felitronics::eq.

#include <felitronics/eq/EqLink.h>

#include <cmath>
#include <cstdio>
#include <vector>

using felitronics::eq::EqLink;

namespace
{
    int failures = 0;

    void check (const char* what, double got, double want, double tol)
    {
        const bool ok = std::fabs (got - want) <= tol;
        if (! ok)
            ++failures;

        std::printf ("%-52s %8.3f  (want %8.3f +-%.2f)  %s\n", what, got, want, tol, ok ? "ok" : "FAIL");
    }

    constexpr double sampleRate = 48000.0;
    constexpr double pi = 3.14159265358979323846;

    /** RMS of a sine after the link, in dB relative to the dry sine — the sound's own answer,
        to hold against magnitudeDb. */
    double processedDb (EqLink& eq, double freq)
    {
        eq.reset();

        constexpr int n = 48000;
        std::vector<float> data ((size_t) n);
        for (int i = 0; i < n; ++i)
            data[(size_t) i] = (float) std::sin (2.0 * pi * freq * i / sampleRate);

        float* chans[] = { data.data() };
        eq.process (chans, 1, n);

        // Skip the first quarter: filter transients and the level ramp settle there.
        double sum = 0.0;
        int counted = 0;
        for (int i = n / 4; i < n; ++i, ++counted)
            sum += (double) data[(size_t) i] * data[(size_t) i];

        const double rms = std::sqrt (sum / counted);
        return 20.0 * std::log10 (std::max (1.0e-9, rms / std::sqrt (0.5)));
    }
}

int main()
{
    std::printf ("felitronics::eq EqLink tests\n\n");

    EqLink eq;
    eq.prepare (sampleRate, 1);

    // A link at rest is bit-transparent — every band self-disables.
    check ("flat: 1 kHz", eq.magnitudeDb (1000.0), 0.0, 1.0e-9);
    check ("flat: 50 Hz", eq.magnitudeDb (50.0),   0.0, 1.0e-9);

    EqLink::Settings s;

    // Shelves with FREE corners — the plateau follows the corner you chose.
    s = {}; s.loDb = 12.0; s.loHz = 100.0; eq.setSettings (s);
    check ("lo +12 @100: 20 Hz (plateau)",   eq.magnitudeDb (20.0),    12.0, 0.6);
    check ("lo +12 @100: 10 kHz (clear)",    eq.magnitudeDb (10000.0),  0.0, 0.2);

    s = {}; s.loDb = 12.0; s.loHz = 400.0; eq.setSettings (s);
    check ("lo +12 @400: 400 Hz (half up)",  eq.magnitudeDb (400.0),    6.0, 1.5);

    s = {}; s.hiDb = -9.0; s.hiHz = 4000.0; eq.setSettings (s);
    check ("hi -9 @4k: 15 kHz (plateau)",    eq.magnitudeDb (15000.0), -9.0, 0.6);
    check ("hi -9 @4k: 200 Hz (clear)",      eq.magnitudeDb (200.0),    0.0, 0.3);

    // Bells at their centres; the narrow B3 only exists while switched on.
    s = {}; s.b1Db = -12.0; s.b1Hz = 400.0; eq.setSettings (s);
    check ("b1 -12 @400: centre",            eq.magnitudeDb (400.0),  -12.0, 0.3);

    s = {}; s.b2Db = 9.0; s.b2Hz = 2500.0; s.b2Q = 2.0; eq.setSettings (s);
    check ("b2 +9 @2.5k Q2: centre",         eq.magnitudeDb (2500.0),   9.0, 0.3);
    check ("b2 +9 @2.5k Q2: 250 Hz (clear)", eq.magnitudeDb (250.0),    0.0, 0.3);

    s = {}; s.b3Db = -15.0; s.b3Hz = 3000.0; s.b3Q = 8.0; eq.setSettings (s);
    check ("b3 OFF stays silent",            eq.magnitudeDb (3000.0),   0.0, 1.0e-9);
    s.b3On = true; eq.setSettings (s);
    check ("b3 ON -15 @3k Q8: centre",       eq.magnitudeDb (3000.0), -15.0, 0.5);
    check ("b3 ON: one octave out (narrow)", eq.magnitudeDb (1500.0),   0.0, 1.0);

    // The slope ladder, measured two octaves under an HPF at 200 Hz.
    const double slopeWant[] = { -12.3, -24.1, -36.3, -48.2, -96.3 };
    const int    slopes[]    = { 6, 12, 18, 24, 48 };

    for (int i = 0; i < 5; ++i)
    {
        s = {}; s.hpfOn = true; s.hpfHz = 200.0; s.hpfSlope = slopes[i];
        eq.setSettings (s);

        char what[64];
        std::snprintf (what, sizeof (what), "hpf %d dB/oct @200: 50 Hz", slopes[i]);
        check (what, eq.magnitudeDb (50.0), slopeWant[i], 3.0);
    }

    s = {}; s.lpfOn = true; s.lpfHz = 2000.0; s.lpfSlope = 48; eq.setSettings (s);
    check ("lpf 48 dB/oct @2k: 8 kHz",       eq.magnitudeDb (8000.0), -96.3, 4.0);
    check ("lpf 48 dB/oct @2k: 500 Hz",      eq.magnitudeDb (500.0),     0.0, 0.3);

    // The invariant: what the curve promises is what the sound does.
    s = {};
    s.loDb = 6.0; s.loHz = 150.0;
    s.b1Db = -9.0; s.b1Hz = 500.0; s.b1Q = 1.4;
    s.hiDb = 4.0; s.hiHz = 6000.0;
    s.hpfOn = true; s.hpfHz = 60.0; s.hpfSlope = 24;
    eq.setSettings (s);

    for (const double f : { 80.0, 500.0, 3000.0, 12000.0 })
    {
        char what[64];
        std::snprintf (what, sizeof (what), "curve == sound @ %.0f Hz", f);
        check (what, processedDb (eq, f), eq.magnitudeDb (f), 0.4);
    }

    std::printf ("\n%s\n", failures != 0 ? "FAILURES" : "all checks passed");
    return failures;
}

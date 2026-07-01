// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// JUCE-free golden tests for eq::Svf HPF/LPF — the exact filters OrbitCab's IRSlot de-JUCEs onto (replacing
// juce::dsp::StateVariableTPTFilter). The decisive property: the MEASURED magnitude matches the analytic
// 2nd-order Butterworth realised by the bilinear transform — a Zavalishin TPT SVF's effective analog
// frequency is tan(pi*f/fs), so the exact digital curve uses the prewarped ratio r = tan(pi*f/fs)/tan(pi*fc/fs)
// (NOT the raw f/fc — that only matches well below Nyquist). This survives after JUCE is removed; the direct
// null against juce lives in the orbitcab test-suite (where JUCE is still linked).

#include <felitronics_test.h>
#include <felitronics/eq/Svf.h>
#include <felitronics/core/Math.h>

#include <cmath>
#include <cstdio>

using namespace felitronics;

// |H(f)| of a prepared single-channel Svf: drive a unit sine at f (fs/f is an exact integer for every f
// below, so an integer number of periods is measured), discard the transient, RMS the steady state.
static double measureMag (eq::Svf& svf, double fs, double f)
{
    svf.reset();
    const int spp  = (int) std::lround (fs / f);         // samples/period — exact for the chosen f
    const int warm = 40 * spp, meas = 64 * spp;
    double acc = 0.0;
    for (int i = 0; i < warm + meas; ++i)
    {
        const double x = std::sin (2.0 * core::kPi * f * (double) i / fs);
        const float  y = svf.processSample (0, (float) x);
        if (i >= warm) acc += (double) y * (double) y;
    }
    return std::sqrt (acc / (double) meas) / std::sqrt (0.5);   // RMS(y) / RMS(unit sine) = |H|
}

static double toDb (double lin) { return 20.0 * std::log10 (lin > 1e-12 ? lin : 1e-12); }

// The BLT (prewarped) 2nd-order Butterworth the SVF realises exactly. r = tan-ratio; Q = 1/sqrt2 → +1 in the
// denominator's r^2 term drops out (Butterworth), leaving 1/sqrt(1+r^4) (LP) and r^2/sqrt(1+r^4) (HP).
static double warpRatio (double f, double fc, double fs)
{
    return std::tan (core::kPi * f / fs) / std::tan (core::kPi * fc / fs);
}
static double lpMag (double f, double fc, double fs) { const double r = warpRatio (f, fc, fs); return 1.0     / std::sqrt (1.0 + r*r*r*r); }
static double hpMag (double f, double fc, double fs) { const double r = warpRatio (f, fc, fs); return (r*r)   / std::sqrt (1.0 + r*r*r*r); }

int main()
{
    std::printf ("felitronics::eq Svf HPF/LPF golden tests\n");
    const double Q = 0.70710678118654752;                 // Butterworth (== juce setResonance(1/sqrt2))
    const double freqs[] = { 125.0, 250.0, 500.0, 1000.0, 2000.0, 4000.0, 8000.0 };

    for (double fs : { 44100.0, 48000.0, 96000.0 })
    {
        const double fc = 1000.0;
        eq::Svf lp; lp.prepare (fs, 1); lp.setParams (eq::FilterType::LowPass,  fc, Q, 0.0);
        eq::Svf hp; hp.prepare (fs, 1); hp.setParams (eq::FilterType::HighPass, fc, Q, 0.0);

        test::group (std::string ("Svf LP/HP == BLT 2nd-order Butterworth @ ") + std::to_string ((int) fs) + " Hz");
        for (double f : freqs)
        {
            if (fs == 44100.0 && f > 8000.0) continue;
            test::approx (toDb (measureMag (lp, fs, f)), toDb (lpMag (f, fc, fs)), 0.05, "LP |H| matches BLT Butterworth");
            test::approx (toDb (measureMag (hp, fs, f)), toDb (hpMag (f, fc, fs)), 0.05, "HP |H| matches BLT Butterworth");
        }
        // Exactly −3.01 dB at the cutoff (the prewarp makes the SVF cutoff sample-accurate).
        test::approx (toDb (measureMag (lp, fs, fc)), -3.0103, 0.03, "LP is −3 dB at the cutoff");
        test::approx (toDb (measureMag (hp, fs, fc)), -3.0103, 0.03, "HP is −3 dB at the cutoff");
    }

    // reset() clears the integrators → the next impulse starts from silence (no leaked tail).
    test::group ("Svf: reset() clears state; NaN/denormal-safe");
    {
        eq::Svf f; f.prepare (48000.0, 1); f.setParams (eq::FilterType::HighPass, 1000.0, Q, 0.0);
        for (int i = 0; i < 500; ++i) f.processSample (0, 1.0f);   // charge the state
        f.reset();
        const float y0 = f.processSample (0, 0.0f);
        test::ok (y0 == 0.0f, "after reset, a zero input yields exactly zero (state cleared)");

        f.reset();
        float acc = 0.0f; for (int i = 0; i < 64; ++i) acc += std::fabs (f.processSample (0, (i == 0) ? 1.0f : 0.0f));
        test::ok (std::isfinite (acc) && acc > 0.0f, "impulse response is finite + non-trivial");
        f.flushDenormals();
        test::ok (std::isfinite (f.processSample (0, 0.0f)), "flushDenormals leaves the filter finite");
    }

    return test::report();
}

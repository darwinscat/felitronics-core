// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// JUCE-free self-tests for felitronics::dynamics (EnvelopeFollower + GainComputer) — measured
// behaviour vs the analytic envelope / compressor-expander law.

#include <felitronics_test.h>
#include <felitronics/dynamics/EnvelopeFollower.h>
#include <felitronics/dynamics/GainComputer.h>
#include <felitronics/dynamics/ChannelLinker.h>
#include <felitronics/dynamics/GainReductionFollower.h>
#include <felitronics/core/Math.h>

#include <cmath>
#include <cstdio>

using namespace felitronics;

int main()
{
    std::printf ("felitronics::dynamics tests\n");
    const double fs = 48000.0;
    const double e1 = std::exp (1.0);

    // ---------------- EnvelopeFollower ----------------
    // Peak attack reaches 1 - 1/e of a unit step after exactly one attack time-constant.
    test::group ("EnvelopeFollower: peak attack 1-1/e @ attackMs");
    {
        dynamics::EnvelopeFollower ef; ef.setDetector (dynamics::Detector::Peak);
        ef.setTimes (10.0, 100.0); ef.prepare (fs);
        const int n = (int) std::lround (0.010 * fs);
        float v = 0.0f; for (int i = 0; i < n; ++i) v = ef.process (1.0f);
        test::approx (v, 1.0 - 1.0 / e1, 0.02, "peak env ~ 1-1/e after one attack TC");
    }

    // Peak release decays to ~1/e of the held value after one release time-constant.
    test::group ("EnvelopeFollower: peak release 1/e @ releaseMs");
    {
        dynamics::EnvelopeFollower ef; ef.setDetector (dynamics::Detector::Peak);
        ef.setTimes (1.0, 50.0); ef.prepare (fs);
        for (int i = 0; i < (int) (0.05 * fs); ++i) ef.process (1.0f);   // charge to ~1
        const float start = ef.envelope();
        const int n = (int) std::lround (0.050 * fs);
        float v = start; for (int i = 0; i < n; ++i) v = ef.process (0.0f);
        test::approx (v, start / e1, 0.03, "peak env decays to ~1/e after one release TC");
    }

    // RMS detector on a steady sine of amplitude A converges to A/sqrt(2).
    test::group ("EnvelopeFollower: RMS of a sine == A/sqrt2");
    {
        dynamics::EnvelopeFollower ef; ef.setDetector (dynamics::Detector::Rms);
        ef.setTimes (5.0, 5.0); ef.prepare (fs);
        const double f = 1000.0, A = 0.5;
        float v = 0.0f;
        for (int i = 0; i < (int) (0.5 * fs); ++i)
            v = ef.process ((float) (A * std::sin (2.0 * core::kPi * f * i / fs)));
        test::approx (v, A / std::sqrt (2.0), 0.02, "RMS env ~ A/sqrt2");
    }

    // Law 8: after a long silence the state flushes to EXACT zero (no sustained subnormals).
    test::group ("EnvelopeFollower: flushDenormals after long silence");
    {
        dynamics::EnvelopeFollower ef; ef.setDetector (dynamics::Detector::Peak);
        ef.setTimes (1.0, 1.0); ef.prepare (fs);
        for (int i = 0; i < (int) fs; ++i) ef.process (1.0f);
        for (int b = 0; b < 4000; ++b) { for (int i = 0; i < 64; ++i) ef.process (0.0f); ef.flushDenormals(); }
        test::ok (ef.envelope() == 0.0f, "env flushed to exact 0 after long silence");
    }

    // ---------------- GainComputer ----------------
    // Down-compress 2:1, thr -20: 10 dB over → -10·(1-1/2) = -5 dB; below thr → 0; at thr → 0 (hard knee).
    test::group ("GainComputer: down-compress 2:1");
    {
        dynamics::GainComputer gc; gc.setMode (dynamics::Mode::DownCompress);
        gc.setThresholdDb (-20.0); gc.setRatio (2.0); gc.setKneeDb (0.0); gc.setRangeDb (24.0);
        test::approx (gc.deltaDb (-10.0), -5.0, 1e-9, "10 dB over @2:1 -> -5 dB");
        test::approx (gc.deltaDb (-30.0),  0.0, 1e-9, "below threshold -> 0");
        test::approx (gc.deltaDb (-20.0),  0.0, 1e-9, "at threshold -> 0 (hard knee)");
    }

    // Range clamps the magnitude of the delta.
    test::group ("GainComputer: range clamp");
    {
        dynamics::GainComputer gc; gc.setMode (dynamics::Mode::DownCompress);
        gc.setThresholdDb (-40.0); gc.setRatio (8.0); gc.setKneeDb (0.0); gc.setRangeDb (6.0);
        test::approx (gc.deltaDb (0.0), -6.0, 1e-9, "huge overshoot clamps to -range");
    }

    // Up-compress lifts below threshold, transparent above.
    test::group ("GainComputer: up-compress");
    {
        dynamics::GainComputer gc; gc.setMode (dynamics::Mode::UpCompress);
        gc.setThresholdDb (-20.0); gc.setRatio (2.0); gc.setKneeDb (0.0); gc.setRangeDb (24.0);
        test::approx (gc.deltaDb (-30.0), +5.0, 1e-9, "10 dB under @2:1 -> +5 dB lift");
        test::approx (gc.deltaDb (-10.0),  0.0, 1e-9, "above threshold -> 0");
    }

    // Down-expand ducks below threshold by under·(ratio-1).
    test::group ("GainComputer: down-expand");
    {
        dynamics::GainComputer gc; gc.setMode (dynamics::Mode::DownExpand);
        gc.setThresholdDb (-40.0); gc.setRatio (2.0); gc.setKneeDb (0.0); gc.setRangeDb (60.0);
        test::approx (gc.deltaDb (-50.0), -10.0, 1e-9, "10 dB under, ratio 2 -> -10 dB");
        test::approx (gc.deltaDb (-30.0),   0.0, 1e-9, "above threshold -> 0");
    }

    // Soft knee: continuous through the corner (value matches at both edges + the analytic centre).
    test::group ("GainComputer: soft knee (C0/C1 corner)");
    {
        dynamics::GainComputer gc; gc.setMode (dynamics::Mode::DownCompress);
        gc.setThresholdDb (0.0); gc.setRatio (2.0); gc.setKneeDb (10.0); gc.setRangeDb (60.0);
        test::approx (gc.deltaDb (-5.0),  0.0,    1e-9, "lower knee edge -> 0");
        test::approx (gc.deltaDb (+5.0), -2.5,    1e-9, "upper knee edge -> hard value (-5*0.5)");
        test::approx (gc.deltaDb ( 0.0), -0.625,  1e-9, "knee centre quadratic (over=1.25, delta=-0.625)");
    }

    // Symmetric soft knee is INTENTIONAL (centred on threshold, smears ±knee/2) — pin the textbook
    // shape so a future "fix" can't silently make it one-sided. (The inactive-side tail can look like a
    // bug, but it is standard soft-knee behaviour.)
    test::group ("GainComputer: symmetric soft knee (intentional)");
    {
        dynamics::GainComputer gc; gc.setMode (dynamics::Mode::UpCompress);
        gc.setThresholdDb (0.0); gc.setRatio (2.0); gc.setKneeDb (10.0); gc.setRangeDb (60.0);
        test::approx (gc.deltaDb (+2.0), +0.225, 1e-9, "knee tail extends ~knee/2 ABOVE threshold");
        test::approx (gc.deltaDb (+5.0),  0.0,   1e-9, "zero beyond +knee/2");
        test::approx (gc.deltaDb (-2.0), +1.225, 1e-9, "fuller lift below threshold");
    }

    // ---------------- ChannelLinker ----------------
    test::group ("ChannelLinker (one linked level → same gain for all channels)");
    {
        float a[1] { 0.5f }, b[1] { -0.8f }; const float* ch2[2] { a, b };
        test::approx (dynamics::linkAmplitude (dynamics::LinkMode::Max, ch2, 2, 0), 0.8, 1e-6, "Max = loudest |ch|");
        test::approx (dynamics::linkAmplitude (dynamics::LinkMode::MeanPower, ch2, 2, 0),
                      std::sqrt ((0.25 + 0.64) / 2.0), 1e-6, "MeanPower = sqrt(mean ch^2)");
        float m[1] { -0.3f }; const float* ch1[1] { m };
        test::approx (dynamics::linkAmplitude (dynamics::LinkMode::Max, ch1, 1, 0), 0.3, 1e-6, "mono = |ch0|");
    }

    // ---------------- GainReductionFollower ----------------
    test::group ("GainReductionFollower attack/release ballistics on the GR");
    {
        dynamics::GainReductionFollower gr; gr.prepare (fs); gr.setTimes (10.0, 100.0);
        for (int i = 0; i < (int) std::lround (0.010 * fs); ++i) gr.process (-6.0f);   // one attack TC
        test::approx (gr.valueDb(), -6.0 * (1.0 - 1.0 / e1), 0.05, "attack reaches 1-1/e toward target");

        dynamics::GainReductionFollower g2; g2.prepare (fs); g2.setTimes (1.0, 50.0);
        for (int i = 0; i < (int) (0.05 * fs); ++i) g2.process (-6.0f);                // settle to ~-6
        const float start = g2.valueDb();
        for (int i = 0; i < (int) std::lround (0.050 * fs); ++i) g2.process (0.0f);     // one release TC
        test::approx (g2.valueDb(), start / e1, 0.1, "release decays to ~1/e toward 0");
    }

    return test::report();
}

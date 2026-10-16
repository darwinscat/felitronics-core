// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.
//
// dehum::HumRemover as a FILTER (the LTI configuration: one width, no tracking — dehumtest::lti): the null it digs at every harmonic, the notch shape against the exact
// transfer of the LMS loop (and the −3 dB point at ±BW/2), what it does to the partials of real notes, a
// sweep away from the harmonics, and the ring an impulse leaves. At 48 and 192 kHz, 50 and 60 Hz mains.
#include <felitronics_test.h>
#include "dehum_fixtures.h"

#include <cstdio>
#include <string>

using namespace dehumtest;
namespace dh = felitronics::dehum;
using felitronics::test::ok;
using felitronics::test::group;

//------------------------------------------------------------------------------------------------------
// NULL DEPTH. Hum of 24 harmonics, amplitude 1/k, random phases, fundamental at −30 dBFS, at exactly 50 or
// 60 Hz (the oscillator runs at the nominal frequency rounded to fs/2³², ~3e-6 Hz off). The level of each
// harmonic in the output over [settle, settle + 1 s] against its level in the input. Returns the worst.
static double nullDepth (double fs, dh::Mains m, double bw, double settleS, double limitDb, bool assert_)
{
    const double f0 = m == dh::Mains::Hz50 ? 50.0 : 60.0;
    Hum hum; hum.f0 = f0; hum.fs = fs; hum.seed (m == dh::Mains::Hz50 ? 11u : 23u);
    const std::int64_t S = (std::int64_t) (settleS * fs), W = (std::int64_t) fs;   // 1 s = whole periods
    std::vector<float> x ((std::size_t) (S + W + (std::int64_t) (0.05 * fs)));
    for (std::size_t n = 0; n < x.size(); ++n) x[n] = (float) hum.at ((std::int64_t) n);
    taper (x, fs);   // the tail end of the taper lies past the window
    const std::vector<float> in = x;

    dh::HumRemover<> h;
    h.setParams (lti (m, bw));   // params first: prepare applies them at once
    felitronics::test::run (h.prepare (fs, 512, 1));
    runMono (h, x);

    double worst = -400.0;
    std::string per;
    for (int k = 1; k <= 24; ++k)
    {
        const double d = dbOf (fitAmp (x, S, W, f0 * k, fs) / fitAmp (in, S, W, f0 * k, fs));
        worst = std::max (worst, d);
        if (k == 1 || k == 2 || k == 3 || k == 11 || k == 24) per += " k" + std::to_string (k) + " " + std::to_string ((int) std::lround (d));
    }
    std::printf ("      null %3.0f kHz %s BW %.1f after %.0f s: worst %.1f dB (%s )\n",
                 fs / 1000.0, f0 > 55.0 ? "60" : "50", bw, settleS, worst, per.c_str());
    if (assert_)
        ok (worst <= limitDb, "null " + std::to_string ((int) fs) + " Hz, mains " + std::to_string ((int) f0) + ", BW "
                              + std::to_string (bw) + ", after " + std::to_string (settleS) + " s: every harmonic <= "
                              + std::to_string (limitDb) + " dB (worst " + std::to_string (worst) + ")");
    return worst;
}

static void testNullDepth()
{
    group ("null depth per harmonic — 24 harmonics, 1/k, random phases, -30 dBFS");
    // The quiet/learn width (3 Hz, tau 0.11 s) is what a hum-only input runs at once the quiet detector
    // lands; in this version the width is set explicitly.
    for (double fs : { 48000.0, 192000.0 })
        for (dh::Mains m : { dh::Mains::Hz50, dh::Mains::Hz60 })
            nullDepth (fs, m, 3.0, 3.0, -50.0, true);
    // The playing width (0.5 Hz, tau 0.64 s): -50 dB needs 5.8 tau = 3.7 s, so after 3 s it reads about
    // -40 dB by construction; it is asserted after 6 s and the 3 s figure is printed.
    for (double fs : { 48000.0, 192000.0 })
    {
        nullDepth (fs, dh::Mains::Hz50, 0.5, 3.0, 0.0, false);
        nullDepth (fs, dh::Mains::Hz50, 0.5, 6.0, -50.0, true);
    }
}

//------------------------------------------------------------------------------------------------------
// NOTCH SHAPE. A tone at k·f_osc + d, settled for 10 s (15.7 tau at 0.5 Hz), fitted over 4 s, in vs out.
// Reference 1: the exact multi-notch transfer (all 24 notches) — 0.02 dB. Reference 2: the single-notch
// textbook shape, which ignores the neighbours and so only holds mid-comb (k = 11) — 0.05 dB. A second
// oracle (RMS ratio) must agree with the fit to 0.01 dB.
static void notchPoint (double fs, int k, double d, double bw, double& measured, double& exact)
{
    dh::HumRemover<> h;
    h.setParams (lti (dh::Mains::Hz50, bw));   // params first: prepare applies them at once
    felitronics::test::run (h.prepare (fs, 512, 1));
    const double fOsc = h.oscillatorHz(), f = k * fOsc + d;
    const std::int64_t S = (std::int64_t) (10.0 * fs), W = (std::int64_t) (4.0 * fs);
    std::vector<float> x ((std::size_t) (S + W + (std::int64_t) (0.02 * fs)));
    for (std::size_t n = 0; n < x.size(); ++n) x[n] = (float) (0.25 * std::sin (2.0 * kPi * std::fmod (f * (double) n / fs, 1.0) + 0.3));
    taper (x, fs);
    const std::vector<float> in = x;
    runMono (h, x);
    measured = dbOf (fitAmp (x, S, W, f, fs) / fitAmp (in, S, W, f, fs));
    exact    = exactGainDb (f, fOsc, h.activeHarmonics(), bw, fs);
    const double rmsDb = dbOf (rms (x, S, W) / rms (in, S, W));
    ok (std::fabs (measured - exact) <= 0.02,
        "notch " + std::to_string ((int) fs) + " Hz k" + std::to_string (k) + " d " + std::to_string (d) + " Hz: "
        + std::to_string (measured) + " dB vs exact transfer " + std::to_string (exact));
    ok (std::fabs (measured - rmsDb) <= 0.01, "  and the RMS oracle agrees (" + std::to_string (rmsDb) + ")");
    if (k == 11)
        ok (std::fabs (measured - singleNotchDb (d, bw)) <= 0.05,
            "  mid-comb, the single-notch shape holds (" + std::to_string (singleNotchDb (d, bw)) + ")");
}

static void testNotchShape()
{
    group ("notch shape vs BW: -3 dB at +-BW/2, against the exact transfer (BW 0.5 Hz)");
    const double bw = 0.5;
    for (double fs : { 48000.0, 192000.0 })
    {
        const std::vector<int> ks = fs > 96000.0 ? std::vector<int> { 11, 24 } : std::vector<int> { 1, 11, 24 };
        const std::vector<double> ds = fs > 96000.0 ? std::vector<double> { bw / 2, -bw / 2, bw }
                                                     : std::vector<double> { bw / 4, bw / 2, -bw / 2, bw, 3 * bw };
        for (int k : ks)
        {
            std::string line;
            for (double d : ds)
            {
                double m = 0.0, e = 0.0;
                notchPoint (fs, k, d, bw, m, e);
                char buf[64]; std::snprintf (buf, sizeof buf, " %+.3f:%.3f", d, m);
                line += buf;
            }
            std::printf ("      %3.0f kHz k%-2d offset:dB%s\n", fs / 1000.0, k, line.c_str());
        }
    }
}

//------------------------------------------------------------------------------------------------------
// PASSBAND: the partials of real notes (8 partials, B = 1e-4, 1/p, random phases), one note at a time,
// fitted jointly. Each partial must match the exact transfer (0.02 dB) — the canceller does nothing to a
// note that a 0.5 Hz notch comb would not — and a partial 1.2 Hz or more from every harmonic loses at most
// 0.2 dB. (Closer than that a 0.5 Hz notch takes more by its own shape: 0.67 dB at 0.61 Hz, 1.24 dB at
// 0.44 Hz — E4's 2nd and 4th partials against 660 and 1320 Hz on 60 Hz mains; printed below.)
static void testPassband()
{
    group ("passband: the partials of E2 / G3 / B3 / E4 (BW 0.5 Hz)");
    const double fs = 48000.0, bw = 0.5;
    const double notes[] = { 82.4069, 195.9977, 246.9417, 329.6276 };
    const char*  names[] = { "E2", "G3", "B3", "E4" };
    for (dh::Mains m : { dh::Mains::Hz50, dh::Mains::Hz60 })
    {
        int checked = 0;
        for (int ni = 0; ni < 4; ++ni)
        {
            dh::HumRemover<> h;
            h.setParams (lti (m, bw));   // params first: prepare applies them at once
            felitronics::test::run (h.prepare (fs, 512, 1));
            const double fOsc = h.oscillatorHz();
            const int    K    = h.activeHarmonics();
            std::vector<double> fr; double ph[8]; Lcg r { 77u + (std::uint32_t) ni };
            for (int p = 1; p <= 8; ++p) { fr.push_back (p * notes[ni] * std::sqrt (1.0 + 1e-4 * p * p)); ph[p - 1] = 2.0 * kPi * r.uniform(); }
            const std::int64_t S = (std::int64_t) (8.0 * fs), W = (std::int64_t) (4.0 * fs);
            std::vector<float> x ((std::size_t) (S + W + (std::int64_t) (0.02 * fs)));
            for (std::size_t n = 0; n < x.size(); ++n)
            {
                double v = 0.0;
                for (int p = 0; p < 8; ++p) v += 0.06 / (p + 1) * std::sin (2.0 * kPi * std::fmod (fr[(std::size_t) p] * (double) n / fs, 1.0) + ph[p]);
                x[n] = (float) v;
            }
            taper (x, fs);
            const std::vector<float> in = x;
            runMono (h, x);
            const auto ao = fitAmps (x, S, W, fr, fs), ai = fitAmps (in, S, W, fr, fs);
            for (int p = 0; p < 8; ++p)
            {
                const double f = fr[(std::size_t) p];
                double dist = 1e9;
                for (int k = 1; k <= K; ++k) dist = std::min (dist, std::fabs (f - k * fOsc));
                const double got = dbOf (ao[(std::size_t) p] / ai[(std::size_t) p]);
                const double ref = exactGainDb (f, fOsc, K, bw, fs);
                const std::string tag = std::string (names[ni]) + " p" + std::to_string (p + 1) + " " + std::to_string (f)
                                      + " Hz (" + (m == dh::Mains::Hz50 ? "50" : "60") + " Hz mains, " + std::to_string (dist) + " Hz from a harmonic)";
                ok (std::fabs (got - ref) <= 0.02, tag + ": " + std::to_string (got) + " dB = exact transfer " + std::to_string (ref));
                if (dist >= 1.2) ok (got >= -0.2, tag + ": within 0.2 dB");
                if (dist < 1.2) std::printf ("      near a harmonic: %s -> %.3f dB (single-notch theory %.3f)\n", tag.c_str(), got, singleNotchDb (dist, bw));
                ++checked;
            }
        }
        ok (checked == 32, "precondition: all 32 partials measured");
    }
}

//------------------------------------------------------------------------------------------------------
// SWEEP AWAY FROM THE HARMONICS: tones 3 Hz or more from every active harmonic, 20 Hz .. 19 kHz.
static void testSweep()
{
    group ("sweep away from harmonics: |delta| < 0.05 dB (BW 0.5 Hz)");
    const double fs = 48000.0, bw = 0.5;
    const double fr[] = { 20.0, 31.3, 77.7, 123.4, 271.9, 333.3, 517.1, 1013.0, 1777.7, 2600.0, 5000.0, 12000.0, 19000.0 };
    for (dh::Mains m : { dh::Mains::Hz50, dh::Mains::Hz60 })
    {
        double worst = 0.0;
        for (double f : fr)
        {
            dh::HumRemover<> h;
            h.setParams (lti (m, bw));   // params first: prepare applies them at once
            felitronics::test::run (h.prepare (fs, 512, 1));
            double dist = 1e9;
            for (int k = 1; k <= h.activeHarmonics(); ++k) dist = std::min (dist, std::fabs (f - k * h.oscillatorHz()));
            ok (dist >= 3.0, "precondition: " + std::to_string (f) + " Hz is 3 Hz or more from every harmonic");
            const std::int64_t S = (std::int64_t) (4.0 * fs), W = (std::int64_t) (2.0 * fs);
            std::vector<float> x ((std::size_t) (S + W + 960));
            for (std::size_t n = 0; n < x.size(); ++n) x[n] = (float) (0.25 * std::sin (2.0 * kPi * std::fmod (f * (double) n / fs, 1.0) + 0.7));
            taper (x, fs);
            const std::vector<float> in = x;
            runMono (h, x);
            const double d = dbOf (fitAmp (x, S, W, f, fs) / fitAmp (in, S, W, f, fs));
            worst = std::max (worst, std::fabs (d));
            ok (std::fabs (d) < 0.05, std::to_string (f) + " Hz: " + std::to_string (d) + " dB");
        }
        std::printf ("      %s Hz mains: worst |delta| %.4f dB over %zu tones\n", m == dh::Mains::Hz50 ? "50" : "60", worst, sizeof fr / sizeof fr[0]);
    }
}

//------------------------------------------------------------------------------------------------------
// IMPULSE RING. A canceller converged on hum is copied; one copy then sees a +1.0 impulse on top of the
// hum. The difference of the two outputs is the impulse response of the running stage: nothing before the
// impulse (bit-exact), the impulse itself, and a ring — one tone per harmonic, each at mu = 2π·BW/fs
// (-83.7 dB at 48 kHz, BW 0.5) — decaying with tau. Per harmonic it must be <= -80 dB re the impulse. The
// tones start in phase, so the time-domain peak of their sum is about K·mu; that figure is printed.
static void testImpulseRing()
{
    group ("impulse: nothing before the input, ring per harmonic <= -80 dB (48 kHz, BW 0.5 Hz)");
    const double fs = 48000.0, bw = 0.5;
    Hum hum; hum.fs = fs; hum.seed (5u);
    dh::HumRemover<> a;
    a.setParams (lti (dh::Mains::Hz50, bw));   // params first: prepare applies them at once
    felitronics::test::run (a.prepare (fs, 512, 1));
    const std::int64_t S = (std::int64_t) (6.0 * fs);
    std::vector<float> pre ((std::size_t) S);
    for (std::int64_t n = 0; n < S; ++n) pre[(std::size_t) n] = (float) hum.at (n);
    runMono (a, pre);
    dh::HumRemover<> b = a;   // no heap: a copy is an exact twin

    const int L = 48000, n0 = 333;
    std::vector<float> xr ((std::size_t) L), xi;
    for (int i = 0; i < L; ++i) xr[(std::size_t) i] = (float) hum.at (S + i);
    xi = xr;
    xi[(std::size_t) n0] += 1.0f;
    const float before = xi[(std::size_t) n0] - xr[(std::size_t) n0];
    // run both up to and including the impulse, then read the weight perturbation it caused
    { float* io[1] { xr.data() }; felitronics::test::run (a.process (io, 1, n0 + 1)); }
    { float* io[1] { xi.data() }; felitronics::test::run (b.process (io, 1, n0 + 1)); }
    double worstW = 0.0;
    for (int k = 1; k <= a.activeHarmonics(); ++k)
    {
        const auto wa = a.weight (0, k), wb = b.weight (0, k);
        worstW = std::max (worstW, std::hypot (wb.re - wa.re, wb.im - wa.im) / (double) before);
    }
    { float* io[1] { xr.data() + n0 + 1 }; felitronics::test::run (a.process (io, 1, L - n0 - 1)); }
    { float* io[1] { xi.data() + n0 + 1 }; felitronics::test::run (b.process (io, 1, L - n0 - 1)); }

    bool causal = true;
    for (int i = 0; i < n0; ++i) causal = causal && core::sameBits (xr[(std::size_t) i], xi[(std::size_t) i]);
    ok (causal, "nothing before the input: the outputs agree bit for bit up to the impulse");
    ok (std::fabs ((double) (xi[(std::size_t) n0] - xr[(std::size_t) n0]) - (double) before) < 1e-6, "the impulse itself passes at unity");

    std::vector<float> d ((std::size_t) L, 0.0f);
    double peak = 0.0;
    for (int i = n0 + 1; i < L; ++i) { d[(std::size_t) i] = xi[(std::size_t) i] - xr[(std::size_t) i]; peak = std::max (peak, (double) std::fabs (d[(std::size_t) i])); }
    // the emitted ring per harmonic, measured on the audio over its first 0.2 s (whole periods)
    std::vector<double> fr;
    for (int k = 1; k <= a.activeHarmonics(); ++k) fr.push_back (k * a.oscillatorHz());
    const auto amps = fitAmps (d, n0 + 1, 9600, fr, fs);
    double worstA = 0.0;
    for (double v : amps) worstA = std::max (worstA, v);
    const double mu = 2.0 * kPi * bw / fs;
    std::printf ("      mu %.1f dB; weight kick per harmonic %.1f dB; emitted ring per harmonic (first 0.2 s) %.1f dB;"
                 " time-domain peak of the sum %.1f dB\n", dbOf (mu), dbOf (worstW), dbOf (worstA), dbOf (peak));
    ok (dbOf (worstW) <= -80.0, "the impulse moves each harmonic's weight by <= -80 dB (" + std::to_string (dbOf (worstW)) + ")");
    ok (dbOf (worstA) <= -80.0, "the emitted ring per harmonic is <= -80 dB (" + std::to_string (dbOf (worstA)) + ")");
    ok (std::fabs (dbOf (worstW) - dbOf (mu)) < 0.1, "and the kick is mu, as the loop says");
}

int main()
{
    std::printf ("felitronics::dehum — HumRemover as a filter\n");
    testNullDepth();
    testNotchShape();
    testPassband();
    testSweep();
    testImpulseRing();
    return felitronics::test::report();
}

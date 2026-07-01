// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// Proves the matched (Vicanek) biquad family does its job: resonant filters and shelves
// pushed near Nyquist stay true to their analog prototype where the RBJ cookbook cramps.
// This is the correctness gate behind decision #3 (matched coeffs from day one).

#include "TestUtil.h"

#include <teq/MatchedBiquad.h>

#include <algorithm>
#include <cmath>
#include <cstdio>

using namespace teq;

namespace
{
    // Analog resonant-lowpass prototype magnitude at real frequency f (Hz). |H(i w0)| = Q.
    double analogLowpassMag (double f, double f0, double Q) noexcept
    {
        const double w0 = 2.0 * kPi * f0, w = 2.0 * kPi * f;
        const double re = w0 * w0 - w * w, im = w * w0 / Q;
        return (w0 * w0) / std::sqrt (re * re + im * im);
    }

    // Analog 2-pole Butterworth high-shelf magnitude (Vicanek 2poleShelvingFits eq. 1),
    // gain G at high frequencies, unity at DC. sg = G^(1/2) = (G^(1/4))^2.
    double analogHighShelfMag (double f, double f0, double G) noexcept
    {
        const double sg = std::sqrt (G), r2 = (f / f0) * (f / f0);
        const double numSq = (1.0 - sg * r2) * (1.0 - sg * r2) + 2.0 * sg * r2;
        const double denSq = (1.0 - r2 / sg) * (1.0 - r2 / sg) + 2.0 * r2 / sg;
        return std::sqrt (numSq / denSq);
    }

    double db (double x) noexcept { return 20.0 * std::log10 (std::max (1e-12, x)); }

    // Analog Butterworth band-stop magnitude (dB) for a prototype order m:
    //   |H(iΩ)|² = 1 / (1 + (BW·Ω / (Ω0² − Ω²))^{2m}),  Ω0 = 2π f0,  BW = Ω0/Q.
    // The independent reference the matched cascade is NULL-tested against.
    double analogBandStopDb (double f, double f0, double Q, int m) noexcept
    {
        const double Om0 = 2.0 * kPi * f0, BW = Om0 / Q, Om = 2.0 * kPi * f;
        const double x = BW * Om / (Om0 * Om0 - Om * Om);
        return 10.0 * std::log10 (1.0 / (1.0 + std::pow (std::fabs (x), 2.0 * (double) m)));
    }

    // |H(e^{jw})| in dB of a cascade of biquad sections (their product response).
    double cascadeDb (const BiquadCoeffs* s, int n, double w) noexcept
    {
        double m = 1.0;
        for (int i = 0; i < n; ++i) m *= s[i].magnitude (w);
        return 20.0 * std::log10 (std::max (1e-300, m));
    }
}

void runMatchedBiquadTests()
{
    using namespace teqtest;
    const double fs = 44100.0;
    auto w = [fs] (double f) { return 2.0 * kPi * f / fs; };

    group ("peaking: 0 dB is identity (flat) everywhere");
    {
        const auto c = matched::peakingDb (3000.0, fs, 2.0, 0.0);
        for (double f : { 50.0, 500.0, 3000.0, 9000.0, 18000.0 })
            expectNear (c.magnitudeDb (w (f)), 0.0, 1e-6, "flat at " + std::to_string ((int) f));
    }

    group ("peaking: nails the centre gain mid-band AND in the air band");
    {
        const auto mid = matched::peakingDb (1000.0, fs, 4.0, 6.0);
        expectNear (mid.magnitudeDb (w (1000.0)), 6.0, 0.02, "+6 dB @ 1k");
        const auto hi = matched::peakingDb (16000.0, fs, 4.0, 6.0);
        expectNear (hi.magnitudeDb (w (16000.0)), 6.0, 0.05, "+6 dB @ 16k air bell");
        expectTrue (hi.isStable(), "air bell stable");
    }

    group ("peaking: a cut is the exact mirror of the boost (symmetric, Pro-Q style)");
    {
        const double f0 = 100.0, Q = 2.0, gDb = 20.0;
        const auto boost = matched::peakingDb (f0, fs, Q,  gDb);
        const auto cut   = matched::peakingDb (f0, fs, Q, -gDb);
        expectNear (cut.magnitudeDb (w (f0)), -gDb, 0.02, "cut hits -20 dB at centre");
        expectTrue (cut.isStable(), "cut stable");
        double maxAsym = 0.0;
        for (double frac = 0.001; frac <= 0.49; frac += 0.002)
        {
            const double ww = 2.0 * kPi * frac;
            maxAsym = std::max (maxAsym, std::fabs (cut.magnitudeDb (ww) + boost.magnitudeDb (ww)));
        }
        std::printf ("      max |cutDb + boostDb| over band: %.6f dB\n", maxAsym);
        expectTrue (maxAsym < 1e-4, "cut(w) == -boost(w) everywhere (true mirror image)");
    }

    group ("notch / all-pass / tilt");
    {
        const auto n = matched::notch (1000.0, fs, 4.0);
        expectTrue (n.magnitudeDb (w (1000.0)) < -40.0,     "notch: deep null at f0");
        expectNear (n.magnitudeDb (0.0),        0.0, 0.05,  "notch: unity at DC");
        expectNear (n.magnitudeDb (w (100.0)),  0.0, 1.0,   "notch: ~unity a decade below");
        expectNear (n.magnitudeDb (w (10000.0)), 0.0, 1.0,  "notch: ~unity a decade above");
        expectTrue (n.isStable(),                           "notch stable");

        const auto ap = matched::allpass (2000.0, fs, 2.0);
        for (double f : { 50.0, 500.0, 2000.0, 9000.0, 18000.0 })
            expectNear (ap.magnitudeDb (w (f)), 0.0, 1e-4, "all-pass flat at " + std::to_string ((int) f));
        expectTrue (ap.isStable(),                         "all-pass stable");

        const double G = 12.0;                             // tilt = lowShelf(-G) * highShelf(+G)
        const auto lo = matched::lowShelfDb  (1000.0, fs, -G);
        const auto hi = matched::highShelfDb (1000.0, fs,  G);
        auto tiltDb = [&] (double f) { return lo.magnitudeDb (w (f)) + hi.magnitudeDb (w (f)); };
        expectNear (tiltDb (40.0),    -G,  1.0, "tilt: lows ~ -12 dB");
        expectNear (tiltDb (16000.0),  G,  1.0, "tilt: highs ~ +12 dB");
        expectNear (tiltDb (1000.0),   0.0, 1.5, "tilt: ~0 at the pivot");
    }

    group ("single notch: sweep f0 / Q — deep null, unity a decade out, log-symmetric skirts");
    {
        for (double f0 : { 60.0, 250.0, 1000.0, 4000.0, 12000.0 })
            for (double Q : { 0.5, 1.0, 4.0, 16.0 })
            {
                const auto c = matched::notch (f0, fs, Q);
                const std::string at = " (f0=" + std::to_string ((int) f0) + " Q=" + std::to_string (Q) + ")";
                expectTrue (c.isStable(),                        "swept-Q notch stable" + at);
                expectTrue (c.magnitudeDb (w (f0)) < -40.0,      "deep null at f0"       + at);
                // unity a decade out on BOTH sides (a decade above must stay below Nyquist).
                if (f0 * 10.0 < 0.45 * fs)
                    expectNear (c.magnitudeDb (w (f0 * 10.0)), 0.0, 1.0, "unity a decade above" + at);
                if (f0 / 10.0 > 20.0)
                    expectNear (c.magnitudeDb (w (f0 / 10.0)), 0.0, 1.0, "unity a decade below" + at);
                // skirts symmetric in LOG frequency: |H(f0·r)| ≈ |H(f0/r)| for a geometric detune r.
                // Holds away from Nyquist; a matched notch trades log-symmetry for analog-magnitude
                // fidelity as f→fs/2, so only assert it where both skirt points stay below 0.20·fs.
                for (double r : { 1.3, 2.0, 3.0 })
                    if (f0 * r < 0.20 * fs && f0 / r > 20.0)
                    {
                        const double up = c.magnitudeDb (w (f0 * r)), dn = c.magnitudeDb (w (f0 / r));
                        expectNear (up, dn, 0.6, "log-symmetric skirt r=" + std::to_string (r) + at);
                    }
            }
    }

    group ("variable-order band-stop: order 1 == matched::notch BIT-FOR-BIT (legacy sessions don't drift)");
    {
        BiquadCoeffs out[8];
        for (double f0 : { 50.0, 200.0, 1000.0, 5000.0, 18000.0 })
            for (double Q : { 0.3, 0.707, 2.0, 8.0, 30.0 })
            {
                const int n = matched::notchCascade (f0, fs, Q, 1, out);
                const auto ref = matched::notch (f0, fs, Q);
                const std::string at = " (f0=" + std::to_string ((int) f0) + " Q=" + std::to_string (Q) + ")";
                expectTrue (n == 1, "order-1 is a single section" + at);
                expectTrue (out[0].b0 == ref.b0 && out[0].b1 == ref.b1 && out[0].b2 == ref.b2
                         && out[0].a1 == ref.a1 && out[0].a2 == ref.a2, "coeffs bit-identical" + at);
            }
    }

    group ("variable-order band-stop: deep null + unity DC/Nyquist + skirt steepens monotonically with order");
    {
        for (double f0 : { 120.0, 1000.0, 6000.0 })
        {
            const double Q = 3.0;
            double prevSteep = -1e9;
            BiquadCoeffs out[8];
            for (int m : { 1, 2, 3, 4, 6, 8 })
            {
                const int n = matched::notchCascade (f0, fs, Q, m, out);
                const std::string at = " (f0=" + std::to_string ((int) f0) + " m=" + std::to_string (m) + ")";
                expectTrue (n == m,                                    "section count == order"    + at);
                bool stable = true; for (int i = 0; i < n; ++i) stable = stable && out[i].isStable();
                expectTrue (stable,                                    "all sections stable"       + at);
                expectTrue (cascadeDb (out, n, w (f0)) < -60.0,        "order-fold deep null at f0" + at);
                expectNear (cascadeDb (out, n, 0.0), 0.0, 0.02,        "unity at DC"               + at);
                expectNear (cascadeDb (out, n, kPi), 0.0, 0.5,         "≈unity at Nyquist"          + at);
                // steepness = dB recovered between a near-null detune and an octave out; grows with order.
                const double steep = cascadeDb (out, n, w (f0 * 2.0)) - cascadeDb (out, n, w (f0 * std::pow (2.0, 1.0 / 6.0)));
                if (m > 1) expectTrue (steep > prevSteep + 0.2, "skirt steeper than lower order" + at);
                prevSteep = steep;
            }
        }
    }

    group ("variable-order band-stop: reference-NULL vs analog Butterworth band-stop");
    {
        BiquadCoeffs out[8];
        // Mid-band notch: the matched cascade tracks the analog prototype tightly away from Nyquist.
        for (int m : { 2, 3, 4, 5, 7, 8 })
        {
            const double f0 = 1000.0, Q = 2.0;
            const int n = matched::notchCascade (f0, fs, Q, m, out);
            double maxErr = 0.0;
            for (double f = 100.0; f < 0.30 * fs; f *= 1.02)
                if (std::fabs (f - f0) > 0.06 * f0)                    // skip the ±ε neighbourhood of the −∞ null
                    maxErr = std::max (maxErr, std::fabs (cascadeDb (out, n, w (f)) - analogBandStopDb (f, f0, Q, m)));
            std::printf ("      m=%d mid-band max|cascade−analog| = %.4f dB\n", m, maxErr);
            expectTrue (maxErr < 0.10, "matched cascade == analog band-stop (m=" + std::to_string (m) + ")");
        }
        // High notch pushing skirts toward Nyquist: still close, looser tolerance (matched ≠ exact at fs/2).
        {
            const double f0 = 9000.0, Q = 3.0; const int m = 4;
            const int n = matched::notchCascade (f0, fs, Q, m, out);
            double maxErr = 0.0;
            for (double f = 2000.0; f < 0.45 * fs; f *= 1.02)
                if (std::fabs (f - f0) > 0.06 * f0)
                    maxErr = std::max (maxErr, std::fabs (cascadeDb (out, n, w (f)) - analogBandStopDb (f, f0, Q, m)));
            std::printf ("      high f0=9k m=4 broadband max|cascade−analog| = %.4f dB\n", maxErr);
            expectTrue (maxErr < 1.0, "high-f0 cascade tracks analog within 1 dB across band");
        }
    }

    group ("variable-order band-stop: stable + deep null across fs / f0 / Q / order grid");
    {
        int unstable = 0; double worstNull = -1e9;
        BiquadCoeffs out[8];
        for (double fsg : { 44100.0, 48000.0, 96000.0 })
            for (double Q : { 0.2, 0.5, 0.707, 2.0, 10.0, 30.0 })
                for (int m : { 2, 3, 4, 6, 8 })
                    for (double f0 = 20.0; f0 <= 0.45 * fsg; f0 *= 1.4)
                    {
                        const int n = matched::notchCascade (f0, fsg, Q, m, out);
                        for (int i = 0; i < n; ++i) if (! out[i].isStable()) ++unstable;
                        worstNull = std::max (worstNull, cascadeDb (out, n, 2.0 * kPi * f0 / fsg));
                    }
        std::printf ("      grid: %d unstable sections; worst null across grid = %.1f dB\n", unstable, worstNull);
        expectTrue (unstable == 0,       "every section stable across the whole grid");
        expectTrue (worstNull < -60.0,   "null stays ≥60 dB deep everywhere on the grid");
    }

    group ("variable-order band-stop: Nyquist edge — f0 near fs/2 stays deep + stable");
    {
        BiquadCoeffs out[8];
        for (double frac : { 0.40, 0.45, 0.48 })
            for (int m : { 2, 4, 8 })
            {
                const double f0 = frac * fs;
                const int n = matched::notchCascade (f0, fs, 3.0, m, out);
                const std::string at = " (f0=" + std::to_string ((int) f0) + " m=" + std::to_string (m) + ")";
                bool stable = true; for (int i = 0; i < n; ++i) stable = stable && out[i].isStable();
                expectTrue (stable,                              "sections stable near Nyquist" + at);
                expectTrue (cascadeDb (out, n, w (f0)) < -60.0,  "null survives near Nyquist"    + at);
                expectNear (cascadeDb (out, n, 0.0), 0.0, 0.05,  "unity at DC near Nyquist"      + at);
            }
    }

    group ("variable-order band-stop: exactly `order` sections for every order 1..8 (odd-order pairing)");
    {
        BiquadCoeffs out[8];
        for (int m = 1; m <= 8; ++m)
            for (double f0 : { 100.0, 1000.0, 9000.0 })
                for (double Q : { 0.5, 3.0, 12.0 })
                {
                    const int n = matched::notchCascade (f0, fs, Q, m, out);
                    const std::string at = " (order=" + std::to_string (m) + " f0=" + std::to_string ((int) f0) + ")";
                    expectTrue (n == m, "returns exactly `order` sections" + at);
                    bool stable = true; for (int i = 0; i < n; ++i) stable = stable && out[i].isStable();
                    expectTrue (stable, "all sections stable" + at);
                    expectTrue (cascadeDb (out, n, w (f0)) < -60.0, "deep null" + at);
                }
    }

    group ("variable-order band-stop: a cut-only filter never boosts above unity (maximally flat)");
    {
        BiquadCoeffs out[8];
        double worstOver = -1e9;
        for (double Q : { 0.05, 0.2, 0.707, 4.0, 40.0 })
            for (int m : { 2, 4, 6, 8 })
                for (double f0 : { 20.0, 200.0, 2000.0, 0.40 * fs, 0.48 * fs })
                {
                    const int n = matched::notchCascade (f0, fs, Q, m, out);
                    for (double f = 5.0; f < 0.5 * fs; f *= 1.03)
                        worstOver = std::max (worstOver, cascadeDb (out, n, w (f)));
                }
        std::printf ("      worst magnitude above 0 dB across the grid = %+.4f dB\n", worstOver);
        expectTrue (worstOver < 0.02, "band-stop never boosts (<= unity everywhere => maximally flat)");
    }

    group ("lowpass: matched hits analog |H(w0)| = Q where RBJ cramps");
    {
        const double f0 = 0.40 * fs, Q = 4.0;
        const double m = matched::lowpass (f0, fs, Q).magnitude (w (f0));
        const double r = rbj::lowpass     (f0, fs, Q).magnitude (w (f0));
        std::printf ("      |H(w0)|: Q=%.1f  matched=%.4f  rbj=%.4f\n", Q, m, r);
        expectNear (m, Q, 1e-3, "matched == Q at cutoff (by construction)");
        expectTrue (std::fabs (m - Q) < std::fabs (r - Q), "matched closer to analog than RBJ");
    }

    group ("lowpass: matched tracks the analog curve better than RBJ across the band");
    {
        const double f0 = 0.32 * fs, Q = 3.0;
        const auto m = matched::lowpass (f0, fs, Q);
        const auto r = rbj::lowpass     (f0, fs, Q);
        double mMax = 0.0, rMax = 0.0;
        for (double frac = 0.02; frac <= 0.47; frac += 0.005)
        {
            const double ww = 2.0 * kPi * frac, ref = analogLowpassMag (frac * fs, f0, Q);
            mMax = std::max (mMax, std::fabs (m.magnitude (ww) - ref));
            rMax = std::max (rMax, std::fabs (r.magnitude (ww) - ref));
        }
        std::printf ("      max |digital-analog|: matched=%.4f  rbj=%.4f\n", mMax, rMax);
        expectTrue (mMax < rMax, "matched max error < RBJ max error");
    }

    group ("high shelf: DC = 0 dB, centre = +G/2, Nyquist = +G");
    {
        const double f0 = 2000.0, gDb = 12.0;
        const auto c = matched::highShelfDb (f0, fs, gDb);
        expectNear (c.magnitudeDb (0.0),     0.0,  0.05, "unity at DC");
        expectNear (c.magnitudeDb (w (f0)),  6.0,  0.6,  "+6 dB at the corner");
        expectNear (c.magnitudeDb (kPi),    12.0,  0.3,  "+12 dB plateau at Nyquist");
        expectTrue (c.isStable(), "high shelf stable");
    }

    group ("low shelf: DC = +G, centre = +G/2, top = 0 dB");
    {
        const double f0 = 200.0, gDb = 12.0;
        const auto c = matched::lowShelfDb (f0, fs, gDb);
        expectNear (c.magnitudeDb (0.0),    12.0,  0.3,  "+12 dB at DC");
        expectNear (c.magnitudeDb (w (f0)),  6.0,  0.6,  "+6 dB at the corner");
        expectNear (c.magnitudeDb (kPi),     0.0,  0.3,  "unity at the top");
        expectTrue (c.isStable(), "low shelf stable");
    }

    group ("high shelf: matched tracks analog better than RBJ near Nyquist");
    {
        const double f0 = 10000.0, gDb = 15.0, G = std::pow (10.0, gDb / 20.0);
        const auto m = matched::highShelfDb (f0, fs, gDb);
        const auto r = rbj::highShelf      (f0, fs, G);
        double mMax = 0.0, rMax = 0.0;
        for (double frac = 0.02; frac <= 0.47; frac += 0.005)
        {
            const double ww = 2.0 * kPi * frac, refDb = db (analogHighShelfMag (frac * fs, f0, G));
            mMax = std::max (mMax, std::fabs (m.magnitudeDb (ww) - refDb));
            rMax = std::max (rMax, std::fabs (r.magnitudeDb (ww) - refDb));
        }
        std::printf ("      max |digital-analog| dB: matched=%.3f  rbj=%.3f\n", mMax, rMax);
        expectTrue (mMax < rMax, "matched max error (dB) < RBJ max error (dB)");
    }

    group ("resonant shelf (Q): plateau exact at any Q; corner == analog Q^2(A-1)^2+A; high Q overshoots");
    {
        const double f0 = 1000.0;
        for (double gDb : { 6.0, 12.0, -6.0, -12.0 })
        {
            const double A = std::pow (10.0, gDb / 20.0);
            for (double Q : { 0.5, 0.707, 2.0, 6.0 })
            {
                const auto lo = matched::lowShelfQDb  (f0, fs, gDb, Q);
                const auto hi = matched::highShelfQDb (f0, fs, gDb, Q);
                const std::string at = " (g=" + std::to_string ((int) gDb) + " Q=" + std::to_string (Q) + ")";

                // plateau is exact at DC / Nyquist regardless of Q
                expectNear (lo.magnitudeDb (0.0),  gDb, 0.05, "low shelf DC == gain"      + at);
                expectNear (lo.magnitudeDb (kPi),  0.0, 0.05, "low shelf Nyquist == 0 dB" + at);
                expectNear (hi.magnitudeDb (kPi),  gDb, 0.10, "high shelf Nyquist == gain"+ at);
                expectNear (hi.magnitudeDb (0.0),  0.0, 0.05, "high shelf DC == 0 dB"     + at);
                expectTrue (lo.isStable() && hi.isStable(), "resonant shelf stable" + at);

                // corner magnitude follows the analog prototype |H(w0)|^2 = Q^2 (A-1)^2 + A
                const double Ab = A >= 1.0 ? A : 1.0 / A;                       // boost-equivalent gain
                const double HcDb = (gDb >= 0.0 ? 1.0 : -1.0) * db (std::sqrt (Q * Q * (Ab - 1.0) * (Ab - 1.0) + Ab));
                expectNear (lo.magnitudeDb (w (f0)), HcDb, 0.20, "low shelf corner == analog"  + at);
                expectNear (hi.magnitudeDb (w (f0)), HcDb, 0.20, "high shelf corner == analog" + at);
            }
        }

        // high Q overshoots past the plateau; low Q stays under it (boost case)
        const auto loHiQ = matched::lowShelfQDb (f0, fs, 6.0, 6.0);
        const auto loLoQ = matched::lowShelfQDb (f0, fs, 6.0, 0.5);
        expectTrue (loHiQ.magnitudeDb (w (f0)) > 6.5,  "high-Q low shelf overshoots above the plateau");
        expectTrue (loLoQ.magnitudeDb (w (f0)) < 6.0,  "low-Q low shelf stays below the plateau");

        // high-Q cut dips below the plateau (mirror of the boost overshoot)
        const auto cutHiQ = matched::lowShelfQDb (f0, fs, -6.0, 6.0);
        expectTrue (cutHiQ.magnitudeDb (w (f0)) < -6.5, "high-Q low shelf cut dips below the plateau");

        // cut is the exact mirror of the boost at the same Q (away from the Nyquist edge case)
        const auto up = matched::lowShelfQDb (f0, fs,  12.0, 3.0);
        const auto dn = matched::lowShelfQDb (f0, fs, -12.0, 3.0);
        double maxAsym = 0.0;
        for (double frac = 0.002; frac <= 0.49; frac += 0.004)
            maxAsym = std::max (maxAsym, std::fabs (up.magnitudeDb (2.0 * kPi * frac) + dn.magnitudeDb (2.0 * kPi * frac)));
        std::printf ("      resonant shelf max |cut + boost|: %.6f dB\n", maxAsym);
        expectTrue (maxAsym < 1e-4, "resonant shelf cut == -boost (true mirror)");
    }

    group ("resonant shelf: stable across freq / Q / gain grid (incl. near Nyquist)");
    {
        for (double f0 : { 40.0, 200.0, 1000.0, 8000.0, 16000.0, 20000.0 })
            for (double Q : { 0.4, 0.707, 2.0, 6.0, 12.0 })
                for (double gDb : { -15.0, -6.0, 6.0, 15.0 })
                {
                    const std::string at = " (f0=" + std::to_string ((int) f0) + " Q=" + std::to_string (Q) + " g=" + std::to_string ((int) gDb) + ")";
                    expectTrue (matched::lowShelfQDb  (f0, fs, gDb, Q).isStable(), "low shelf Q unstable"  + at);
                    expectTrue (matched::highShelfQDb (f0, fs, gDb, Q).isStable(), "high shelf Q unstable" + at);
                }
    }

    group ("resonant shelf: plateaus stay EXACT near Nyquist at low Q + high gain (numerator feasibility clamp)");
    {
        // Regression for the numerator-feasibility bug: near Nyquist a low-Q, high-gain shelf drove the
        // single-biquad 3-point fit infeasible (W²+B2<0); the plateau then ballooned (a +30 dB high
        // shelf read +100 dB at DC → overload/NaN). The DC/Nyquist endpoints must be exact at any Q.
        for (double f0 : { 0.30 * fs, 0.40 * fs, 0.45 * fs, 0.49 * fs })
            for (double gDb : { 6.0, 15.0, 30.0 })
                for (double Q : { 0.05, 0.2, 0.5, 0.707 })
                {
                    const auto hs = matched::highShelfQDb (f0, fs, gDb, Q);
                    const auto ls = matched::lowShelfQDb  (f0, fs, gDb, Q);
                    const std::string at = " (f0=" + std::to_string ((int) f0) + " g=" + std::to_string ((int) gDb) + " Q=" + std::to_string (Q) + ")";
                    expectTrue (hs.isStable() && ls.isStable(),    "shelves stable near Nyquist"     + at);
                    expectNear (hs.magnitudeDb (0.0), 0.0, 0.1,    "high shelf DC plateau == 0 dB"    + at);
                    expectNear (hs.magnitudeDb (kPi), gDb, 0.1,    "high shelf Nyquist plateau == g"  + at);
                    expectNear (ls.magnitudeDb (kPi), 0.0, 0.1,    "low shelf Nyquist plateau == 0 dB"+ at);
                    expectNear (ls.magnitudeDb (0.0), gDb, 0.1,    "low shelf DC plateau == gain"     + at);
                }
    }

    group ("stability across a freq / Q / gain grid");
    {
        for (double f0 : { 40.0, 200.0, 1000.0, 8000.0, 18000.0, 21000.0 })
            for (double Q : { 0.3, 0.707, 2.0, 8.0 })
            {
                const std::string at = " (f0=" + std::to_string ((int) f0) + ", Q=" + std::to_string (Q) + ")";
                expectTrue (matched::peakingDb (f0, fs, Q, 6.0).isStable(), "peaking unstable" + at);
                expectTrue (matched::lowpass   (f0, fs, Q).isStable(),      "lowpass unstable" + at);
                expectTrue (matched::highpass  (f0, fs, Q).isStable(),      "highpass unstable" + at);
                expectTrue (matched::bandpass  (f0, fs, Q).isStable(),      "bandpass unstable" + at);
            }
        for (double f0 : { 60.0, 500.0, 3000.0, 12000.0, 18000.0 })
            for (double gDb : { -12.0, -3.0, 3.0, 12.0 })
            {
                const std::string at = " (f0=" + std::to_string ((int) f0) + ", g=" + std::to_string (gDb) + ")";
                expectTrue (matched::highShelfDb (f0, fs, gDb).isStable(), "high shelf unstable" + at);
                expectTrue (matched::lowShelfDb  (f0, fs, gDb).isStable(), "low shelf unstable" + at);
            }
    }
}

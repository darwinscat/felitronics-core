// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// The per-lane dynamic gain-delta seam (EqBand::setLaneDeltaDb). What matters here is not that a bell
// moves — Svf is tested elsewhere — but WHERE the delta sits: inside the lane, after the matched static
// sections and before the M/S fold. Everything below fails if it is anywhere else.

#include <felitronics_test.h>
#include <felitronics/eq/EqBand.h>
#include <felitronics/core/Math.h>

#include <cmath>
#include <cstdio>
#include <limits>
#include <vector>

using namespace felitronics;
using namespace felitronics::eq;

static BandParams bellPoint (double freq, double Q, double gainDb)
{
    BandParams p;
    p.on = true; p.type = FilterType::Bell;
    LaneParams& st = p.lane (Lane::Stereo);
    st.on = true; st.freq = freq; st.Q = Q; st.gainDb = gainDb;
    return p;
}

// Run a band over noise and return the two channels.
static void run (EqBand& b, std::vector<float>& L, std::vector<float>& R)
{
    float* ch[2] { L.data(), R.data() };
    b.processBlock (ch, 2, (int) L.size());
}

static void fillNoise (std::vector<float>& L, std::vector<float>& R, unsigned seed)
{
    unsigned s = seed;
    for (std::size_t i = 0; i < L.size(); ++i)
    {
        s = s * 1664525u + 1013904223u;
        L[i] = (float) ((double) (s >> 8) / 8388608.0 - 1.0) * 0.25f;
        s = s * 1664525u + 1013904223u;
        R[i] = (float) ((double) (s >> 8) / 8388608.0 - 1.0) * 0.25f;
    }
}

int main()
{
    std::printf ("felitronics::eq EqBand gain-delta tests\n");
    const double fs = 48000.0;
    const int    n  = 4096;

    // A band that never opted in must be untouched by a delta — the whole "additive, costs nothing"
    // claim rests on this.
    test::group ("dyn.on == false ignores the delta, bit for bit");
    {
        std::vector<float> aL (n), aR (n), bL (n), bR (n);
        fillNoise (aL, aR, 12345); bL = aL; bR = aR;

        EqBand ref; ref.prepare (fs, 2); ref.setParams (bellPoint (1000.0, 2.0, 4.0));
        run (ref, aL, aR);

        EqBand dut; dut.prepare (fs, 2); dut.setParams (bellPoint (1000.0, 2.0, 4.0));
        dut.setLaneDeltaDb (Lane::Stereo, -12.0);        // shouted at, and deaf by design
        run (dut, bL, bR);

        bool identical = true;
        for (int i = 0; i < n; ++i) identical = identical && (aL[i] == bL[i]) && (aR[i] == bR[i]);
        test::ok (identical, "output is bit-identical with dynamics off");
    }

    test::group ("a zero delta is transparent when dynamics IS on");
    {
        std::vector<float> aL (n), aR (n), bL (n), bR (n);
        fillNoise (aL, aR, 999); bL = aL; bR = aR;

        BandParams p = bellPoint (800.0, 1.5, -3.0);
        EqBand ref; ref.prepare (fs, 2); ref.setParams (p);
        run (ref, aL, aR);

        p.dyn.on = true;
        EqBand dut; dut.prepare (fs, 2); dut.setParams (p);
        dut.setLaneDeltaDb (Lane::Stereo, 0.0);
        run (dut, bL, bR);

        double worst = 0.0;
        for (int i = 0; i < n; ++i) worst = std::fmax (worst, (double) std::fabs (aL[i] - bL[i]));
        test::ok (worst < 1.0e-6, "unity-gain delta bell passes the signal through");
    }

    // The delta must land INSIDE the lane. If it were applied after the fold — or on the wrong axis —
    // a Mid-only delta would leak into Side.
    test::group ("a Mid-lane delta moves Mid and leaves the Side axis bit-exact");
    {
        BandParams p;
        p.on = true; p.type = FilterType::Bell; p.dyn.on = true;
        p.lane (Lane::Stereo).on = false;
        LaneParams& m = p.lane (Lane::Mid);
        m.on = true; m.freq = 1000.0; m.Q = 2.0; m.gainDb = 0.0;

        std::vector<float> aL (n), aR (n), bL (n), bR (n);
        fillNoise (aL, aR, 4242); bL = aL; bR = aR;

        EqBand ref; ref.prepare (fs, 2); ref.setParams (p);
        run (ref, aL, aR);

        EqBand dut; dut.prepare (fs, 2); dut.setParams (p);
        dut.setLaneDeltaDb (Lane::Mid, -9.0);
        run (dut, bL, bR);

        double midMoved = 0.0, sideDrift = 0.0;
        for (int i = 0; i < n; ++i)
        {
            const double refS = 0.5 * ((double) aL[i] - (double) aR[i]);
            const double dutS = 0.5 * ((double) bL[i] - (double) bR[i]);
            sideDrift = std::fmax (sideDrift, std::fabs (refS - dutS));
            midMoved  = std::fmax (midMoved,
                                   std::fabs (0.5 * ((double) aL[i] + (double) aR[i])
                                            - 0.5 * ((double) bL[i] + (double) bR[i])));
        }
        // Not bit-exactness: the fold writes L += dM and R += dM, and (L+dM)-(R+dM) does not cancel
        // exactly in float. What matters is the ORDER of magnitude — a delta leaking onto the wrong
        // axis would show up at the same scale as midMoved (1e-1), not at rounding noise (1e-8).
        test::ok (midMoved > 1.0e-4, "the Mid axis actually moved");
        test::ok (sideDrift < 1.0e-6, "the Side axis only carries rounding noise, not the delta");
        test::ok (midMoved > sideDrift * 1.0e3, "the Mid/Side separation is orders of magnitude, not luck");
    }

    // The static curve must remain the MATCHED design. Running the whole gain through one SVF instead
    // (what a simpler dynamic band does) would show up as a different response near Nyquist.
    test::group ("the static response is unchanged by opting into dynamics");
    {
        BandParams p = bellPoint (12000.0, 1.0, 8.0);    // high shelf-ish bell, where cramping shows
        EqBand a; a.prepare (fs, 2); a.setParams (p);
        p.dyn.on = true;
        EqBand b; b.prepare (fs, 2); b.setParams (p);
        b.setLaneDeltaDb (Lane::Stereo, 0.0);

        for (double f : { 6000.0, 12000.0, 18000.0, 23000.0 })
        {
            const double w = 2.0 * core::kPi * f / fs;
            const double da = 20.0 * std::log10 (std::abs (a.response (w, Axis::Mid)));
            const double db = 20.0 * std::log10 (std::abs (b.response (w, Axis::Mid)));
            test::approx (db, da, 1.0e-9, "matched static response is untouched by the delta seam");
        }
    }

    test::group ("hostile deltas are clamped, not propagated");
    {
        EqBand b; b.prepare (fs, 2);
        BandParams p = bellPoint (1000.0, 2.0, 0.0); p.dyn.on = true;
        b.setParams (p);
        b.setLaneDeltaDb (Lane::Stereo, std::numeric_limits<double>::quiet_NaN());
        test::approx (b.laneDeltaDb (Lane::Stereo), 0.0, 1.0e-12, "NaN delta becomes 0");
        b.setLaneDeltaDb (Lane::Stereo, 1.0e9);
        test::ok (std::isfinite (b.laneDeltaDb (Lane::Stereo)) && b.laneDeltaDb (Lane::Stereo) <= 60.0,
                  "absurd delta is clamped");

        std::vector<float> L (n), R (n);
        fillNoise (L, R, 7);
        run (b, L, R);
        bool finite = true;
        for (int i = 0; i < n; ++i) finite = finite && std::isfinite (L[i]) && std::isfinite (R[i]);
        test::ok (finite, "and the audio stays finite");
    }

    return test::report();
}

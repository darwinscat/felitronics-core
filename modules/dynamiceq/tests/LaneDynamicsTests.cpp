// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// felitronics::dynamiceq::LaneDynamics — the composition layer. The load-bearing test is LEVEL
// TRANSLATION INVARIANCE: the same programme at wildly different absolute levels must be processed
// the same way, because that is the entire promise of one shared setting across lanes 20+ dB apart.

#include <felitronics_test.h>
#include <felitronics/dynamiceq/LaneDynamics.h>
#include <felitronics/core/Math.h>

#include <cmath>
#include <cstdio>
#include <vector>

using namespace felitronics;

static eq::BandParams dynBell (double freq, double Q, double rangeDb, bool on = true)
{
    eq::BandParams p;
    p.on = true; p.type = eq::FilterType::Bell;
    eq::LaneParams& st = p.lane (eq::Lane::Stereo);
    st.on = true; st.freq = freq; st.Q = Q; st.gainDb = 0.0;
    p.dyn.on = on; p.dyn.rangeDb = rangeDb;
    return p;
}

static double t_elapsed (int samples, double fs) noexcept;

// Drive a tone at the band centre and report the DEEPEST delta reached during periodic bursts that
// sit `burstDb` above the programme. Bursts matter: a stationary tone becomes the programme and the
// band correctly idles, so measuring that would make every level look alike for the wrong reason.
static double peakDelta (double fs, double levelDb, double freq, double Q, double rangeDb,
                         double burstDb = 12.0, double sec = 8.0)
{
    eq::EqBand band; band.prepare (fs, 2);
    const eq::BandParams p = dynBell (freq, Q, rangeDb);
    band.setParams (p);

    dynamiceq::LaneDynamics dyn; dyn.prepare (fs, 2); dyn.setParams (p);

    const int block = 64, n = (int) (fs * sec);
    std::vector<float> L ((size_t) block), R ((size_t) block), sl ((size_t) block), sr ((size_t) block);
    long phase = 0;
    double deepest = 0.0;
    for (int done = 0; done < n; done += block)
    {
        for (int i = 0; i < block; ++i, ++phase)
        {
            const double t     = (double) phase / fs;
            const bool   burst = std::fmod (t, 1.0) > 0.75;          // 250 ms of burst per second
            const double amp   = core::dbToGain (levelDb + (burst ? burstDb : 0.0));
            const float  v     = (float) (amp * std::sin (2.0 * core::kPi * freq * t));
            L[(size_t) i] = R[(size_t) i] = sl[(size_t) i] = sr[(size_t) i] = v;
        }
        float* aud[2] { L.data(), R.data() };
        const float* sc[2] { sl.data(), sr.data() };
        dyn.processBand (aud, sc, 2, block, band);
        if (t_elapsed (done, fs) > 3.0)                              // let the estimator settle first
            deepest = std::fmin (deepest, dyn.deltaDb (eq::Lane::Stereo));
    }
    return deepest;
}

static double t_elapsed (int samples, double fs) noexcept { return (double) samples / fs; }

// Drive a steady tone and report the settled delta.
static double settledDelta (double fs, double levelDb, double freq, double Q, double rangeDb,
                            double sec = 6.0)
{
    eq::EqBand band; band.prepare (fs, 2);
    const eq::BandParams p = dynBell (freq, Q, rangeDb);
    band.setParams (p);

    dynamiceq::LaneDynamics dyn; dyn.prepare (fs, 2); dyn.setParams (p);

    const int block = 64, n = (int) (fs * sec);
    std::vector<float> L ((size_t) block), R ((size_t) block), sl ((size_t) block), sr ((size_t) block);
    const double amp = core::dbToGain (levelDb);
    long phase = 0;
    for (int done = 0; done < n; done += block)
    {
        for (int i = 0; i < block; ++i, ++phase)
        {
            const float v = (float) (amp * std::sin (2.0 * core::kPi * freq * (double) phase / fs));
            L[(size_t) i] = R[(size_t) i] = sl[(size_t) i] = sr[(size_t) i] = v;
        }
        float* aud[2] { L.data(), R.data() };
        const float* sc[2] { sl.data(), sr.data() };
        dyn.processBand (aud, sc, 2, block, band);
    }
    return dyn.deltaDb (eq::Lane::Stereo);
}

int main()
{
    std::printf ("felitronics::dynamiceq::LaneDynamics tests\n");
    const double fs = 48000.0;

    // THE promise: identical treatment at -50 / -30 / -10 dBFS. A fixed threshold would give three
    // completely different answers here.
    test::group ("level-translation invariance");
    {
        const double a = peakDelta (fs, -50.0, 3000.0, 4.0, -12.0);
        const double b = peakDelta (fs, -30.0, 3000.0, 4.0, -12.0);
        const double c = peakDelta (fs, -10.0, 3000.0, 4.0, -12.0);
        // First: the band must actually be WORKING, or "identical" would be a statement about zero.
        test::ok (a < -1.0, "the burst genuinely drives gain reduction");
        test::approx (b, a, 1.0, "-30 dBFS is treated like -50");
        test::approx (c, a, 1.0, "-10 dBFS is treated like -50");
    }

    // A steady tone eventually IS the programme, so an auto-threshold band settles near idle. That is
    // the documented semantic, not a bug: permanent problems are what the manual threshold is for.
    test::group ("a stationary tone settles toward idle under auto threshold");
    {
        const double d = settledDelta (fs, -20.0, 3000.0, 4.0, -12.0);
        test::ok (std::fabs (d) < 2.0, "settled delta is small once the tone became the norm");
        test::ok (std::isfinite (d), "and finite");
    }

    test::group ("dynamics off is bit-identical to a plain band");
    {
        const int n = 2048;
        std::vector<float> aL ((size_t) n), aR ((size_t) n), bL ((size_t) n), bR ((size_t) n);
        for (int i = 0; i < n; ++i)
        {
            const float v = (float) (0.2 * std::sin (2.0 * core::kPi * 1000.0 * (double) i / fs));
            aL[(size_t) i] = aR[(size_t) i] = bL[(size_t) i] = bR[(size_t) i] = v;
        }
        eq::BandParams p = dynBell (1000.0, 2.0, -12.0, /*on*/ false);
        p.lane (eq::Lane::Stereo).gainDb = 5.0;

        eq::EqBand ref; ref.prepare (fs, 2); ref.setParams (p);
        float* r[2] { aL.data(), aR.data() };
        ref.processBlock (r, 2, n);

        eq::EqBand dut; dut.prepare (fs, 2); dut.setParams (p);
        dynamiceq::LaneDynamics dyn; dyn.prepare (fs, 2); dyn.setParams (p);
        float* d[2] { bL.data(), bR.data() };
        const float* sc[2] { bL.data(), bR.data() };
        dyn.processBand (d, sc, 2, n, dut);

        bool identical = true;
        for (int i = 0; i < n; ++i) identical = identical && (aL[(size_t) i] == bL[(size_t) i]);
        test::ok (identical, "an opted-out point is untouched by the composition layer");
    }

    // Was vacuous once: with a POSITIVE range and an inverted DownCompress, far-below-threshold
    // silence already yields a zero target, so the test passed even with the gate removed. It now
    // drives the gate directly — a very quiet but NON-silent tail, where the estimator is held and
    // the relative level would otherwise read as a large excess.
    test::group ("a held estimator cannot manufacture gain reduction on a decaying tail");
    {
        eq::BandParams p = dynBell (3000.0, 4.0, +12.0);      // positive range: lift when loud
        eq::EqBand band; band.prepare (fs, 2); band.setParams (p);
        dynamiceq::LaneDynamics dyn; dyn.prepare (fs, 2); dyn.setParams (p);

        const int block = 64;
        std::vector<float> L ((size_t) block), R ((size_t) block);
        long phase = 0;
        auto feed = [&] (double amp, double sec)
        {
            for (int done = 0; done < (int) (fs * sec); done += block)
            {
                for (int i = 0; i < block; ++i, ++phase)
                {
                    const float v = (float) (amp * std::sin (2.0 * core::kPi * 3000.0 * (double) phase / fs));
                    L[(size_t) i] = R[(size_t) i] = v;
                }
                float* aud[2] { L.data(), R.data() };
                const float* sc[2] { L.data(), R.data() };
                dyn.processBand (aud, sc, 2, block, band);
            }
        };
        feed (core::dbToGain (-20.0), 5.0);       // programme: the estimator learns -20
        feed (core::dbToGain (-115.0), 3.0);      // a tail BELOW the activity floor, but NOT digital zero
        test::ok (std::fabs (dyn.deltaDb (eq::Lane::Stereo)) < 1.0,
                  "the gate holds the delta at 0 on a sub-floor tail");

        // The converse, so the gate is not merely "always off": genuine programme still engages.
        feed (core::dbToGain (-20.0), 4.0);
        feed (core::dbToGain (-4.0), 0.4);        // a real 16 dB event over the learned programme
        test::ok (dyn.deltaDb (eq::Lane::Stereo) > 1.0, "a real event still lifts (range is positive)");
    }

    // THE CADENCE BUG. GainReductionFollower's coefficients are PER SAMPLE. Advancing it once per
    // control chunk stretched every time constant by kControl — a 1 ms attack behaved like 16 ms —
    // and made the amount depend on the host block size, since a 17-sample block splits 16+1 while a
    // 64-sample one splits 16+16+16+16. Both properties are pinned here.
    test::group ("attack timing is real, and independent of the host block size");
    {
        // Time to reach 3 dB of reduction after a 20 dB step over the learned programme. An absolute
        // threshold, deliberately: it needs no reference run and cannot degenerate when the settled
        // value is small.
        auto attackMs = [&] (int block)
        {
            eq::BandParams p = dynBell (3000.0, 4.0, -12.0);
            eq::EqBand band; band.prepare (fs, 2); band.setParams (p);
            dynamiceq::LaneDynamics dyn; dyn.prepare (fs, 2); dyn.setParams (p);

            std::vector<float> L ((size_t) block), R ((size_t) block);
            long phase = 0;
            auto feed = [&] (double amp, double sec, bool timed) -> double
            {
                const int n = (int) (fs * sec);
                for (int done = 0; done < n; done += block)
                {
                    for (int i = 0; i < block; ++i, ++phase)
                    {
                        const float v = (float) (amp * std::sin (2.0 * core::kPi * 3000.0 * (double) phase / fs));
                        L[(size_t) i] = R[(size_t) i] = v;
                    }
                    float* aud[2] { L.data(), R.data() };
                    const float* sc[2] { L.data(), R.data() };
                    dyn.processBand (aud, sc, 2, block, band);
                    if (timed && dyn.deltaDb (eq::Lane::Stereo) <= -3.0)
                        return 1.0e3 * (double) (done + block) / fs;   // end of the chunk that got there
                }
                return -1.0;
            };
            feed (core::dbToGain (-30.0), 6.0, false);          // learn the programme
            return feed (core::dbToGain (-10.0), 0.5, true);    // step, and time the rise
        };

        const double a64 = attackMs (64), a17 = attackMs (17), a15 = attackMs (15);
        test::ok (a64 > 0.0 && a17 > 0.0 && a15 > 0.0, "attack measurable at every block size");
        // Auto attack at 3 kHz / Q 4 is a few ms. The pre-fix code advanced the follower once per
        // 16-sample chunk, stretching every constant by ~16x — this bound is what fails then.
        test::ok (a64 < 25.0, "attack is the order the ballistics asked for, not 16x slower");
        test::ok (std::fabs (a17 - a64) < 3.0, "a 17-sample block times like a 64-sample one");
        test::ok (std::fabs (a15 - a64) < 3.0, "a 15-sample block times like a 64-sample one");
    }

    // Every other test drives the ST lane only — a sign error in laneSignal(Side) would have passed
    // the whole suite. These exercise the derivations the product actually rests on.
    test::group ("the M/S and L/R lanes detect their own domain");
    {
        auto lanePeak = [&] (eq::Lane lane, bool sideOnlyContent)
        {
            eq::BandParams p;
            p.on = true; p.type = eq::FilterType::Bell;
            p.lane (eq::Lane::Stereo).on = false;
            eq::LaneParams& lp = p.lane (lane);
            lp.on = true; lp.freq = 3000.0; lp.Q = 4.0; lp.gainDb = 0.0;
            p.dyn.on = true; p.dyn.rangeDb = -12.0;

            eq::EqBand band; band.prepare (fs, 2); band.setParams (p);
            dynamiceq::LaneDynamics dyn; dyn.prepare (fs, 2); dyn.setParams (p);

            const int block = 64;
            std::vector<float> L ((size_t) block), R ((size_t) block);
            long phase = 0;
            double deepest = 0.0;
            const int total = (int) (fs * 8.0);
            for (int done = 0; done < total; done += block)
            {
                for (int i = 0; i < block; ++i, ++phase)
                {
                    const double t = (double) phase / fs;
                    const bool burst = std::fmod (t, 1.0) > 0.75;
                    const double amp = core::dbToGain (-30.0 + (burst ? 14.0 : 0.0));
                    const float v = (float) (amp * std::sin (2.0 * core::kPi * 3000.0 * t));
                    // Side-only content: L = +v, R = -v (mid cancels). Otherwise mono (side cancels).
                    L[(size_t) i] = v;
                    R[(size_t) i] = sideOnlyContent ? -v : v;
                }
                float* aud[2] { L.data(), R.data() };
                const float* sc[2] { L.data(), R.data() };
                dyn.processBand (aud, sc, 2, block, band);
                if ((double) done / fs > 3.0) deepest = std::fmin (deepest, dyn.deltaDb (lane));
            }
            return deepest;
        };

        test::ok (lanePeak (eq::Lane::Mid,  false) < -1.0, "a Mid lane engages on mono content");
        test::ok (lanePeak (eq::Lane::Mid,  true)  > -0.5, "and stays idle when the content is side-only");
        test::ok (lanePeak (eq::Lane::Side, true)  < -1.0, "a Side lane engages on side-only content");
        test::ok (lanePeak (eq::Lane::Side, false) > -0.5, "and stays idle when the content is mono");
        test::ok (lanePeak (eq::Lane::Left, false) < -1.0, "an L lane engages");
    }

    return test::report();
}

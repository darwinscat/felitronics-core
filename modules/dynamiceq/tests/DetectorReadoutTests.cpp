// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// JUCE-free self-tests for DynamicEqBand::detectorLevelDb(), the level the band's detector sees.
//   1. Fed to an INDEPENDENT GainComputer + GainReductionFollower, the readout reproduces the band's own
//      dynamicDeltaDb() bit for bit, sample by sample, across audio and pauses: it is the curve's argument.
//   2. A sine reads the analytic band-pass level through the envelope's steady state; silence lets it fall
//      to the floor.
//   3. A pause decays it along the envelope's own geometric release; reset() returns it to the floor a
//      cold envelope reads, exactly -180 dB.
//   4. Hostile input (NaN, inf, FLT_MAX) never makes it non-finite.
//   5. process(), a pause and the getter allocate nothing.

#include <felitronics_test.h>
#include <alloc_counter.h>   // installs the allocation counter: EVERY form of `new`, over-aligned included
#include <felitronics/core/Math.h>
#include <felitronics/dynamiceq/DynamicEqBand.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <string>
#include <vector>

using namespace felitronics;
using DB   = dynamiceq::DynamicEqBand;
using DP   = dynamiceq::DynamicEqBandParams;
using Mode = dynamiceq::DynamicEqMode;
using Det  = dynamics::Detector;
using Link = dynamics::LinkMode;

namespace
{
constexpr double kSr = 48000.0;

// |H| of the sidechain band-pass, from the textbook rather than from eq::Svf: the analog prototype
// k s / (s^2 + k s + 1), k = 1/Q, at the bilinear transform's prewarped frequency ratio.
double bandPassGain (double f, double fc, double Q)
{
    const double w = std::tan (core::kPi * f / kSr) / std::tan (core::kPi * fc / kSr), k = 1.0 / Q;
    return k * w / std::sqrt ((1.0 - w * w) * (1.0 - w * w) + (k * w) * (k * w));
}

// The steady envelope of a sine of amplitude a through the band's symmetric 2 ms level window: an RMS
// detector settles on a / sqrt 2, the rectifying peak follower on the rectified mean 2a / pi.
double steadyLevelDb (Det d, double a) { return 20.0 * std::log10 (d == Det::Rms ? a / std::sqrt (2.0) : 2.0 * a / core::kPi); }

// dB the readout loses per sample of zero input: the envelope's state is multiplied by the release
// coefficient of its 2 ms window, and an RMS level is the square root of that state.
double releaseDbPerSample (Det d) { return (d == Det::Rms ? 10.0 : 20.0) * std::log10 (std::exp (-1.0 / (0.002 * kSr))); }

const char* name (Det d)  { return d == Det::Rms ? "rms" : "peak"; }
const char* name (Link l) { return l == Link::Max ? "max" : "meanpower"; }
const char* name (Mode m) { return m == Mode::CutWhenLoud ? "cut-when-loud" : m == Mode::BoostWhenLoud ? "boost-when-loud" : "boost-when-quiet"; }

// What a prepared band that has processed nothing reads.
double coldReadout()
{
    DB b; test::run (b.prepare (kSr, 2));
    return b.detectorLevelDb();
}

void runTone (DB& b, double f, double a, int n, int& phase)
{
    std::vector<float> x (64);
    for (int o = 0; o < n; o += 64)
    {
        const int len = std::min (64, n - o);
        for (int j = 0; j < len; ++j, ++phase) x[(size_t) j] = (float) (a * std::sin (2.0 * core::kPi * f * phase / kSr));
        float* io[1] { x.data() };
        test::run (b.process (io, 1, len));
    }
}
} // namespace

int main()
{
    std::printf ("felitronics::dynamiceq detector readout tests\n");
    const double cold = coldReadout();

    // WHAT IS COMPARED: with coeffUpdatePeriod = 1 the band publishes its smoothed delta every sample, so
    // dynamicDeltaDb() after a one-sample call IS that sample's follower output. The reference takes ONLY
    // the readout, runs it through its own GainComputer (same curve) and its own GainReductionFollower
    // (same times, flushed at the same block ends), and must land on the same double — no tolerance.
    // The pauses start where the tone stops, so a readout that froze during one would hand the
    // reference a level ~20 dB too high while the band's own envelope keeps falling.
    test::group ("readout -> independent curve + ballistics == the band's own delta, bit for bit");
    for (Mode mode : { Mode::CutWhenLoud, Mode::BoostWhenLoud, Mode::BoostWhenQuiet })
    for (Det det : { Det::Rms, Det::Peak })
    for (Link link : { Link::Max, Link::MeanPower })
    {
        DP p; p.freq = 5000.0; p.Q = 2.0; p.mode = mode; p.thresholdDb = -30.0; p.ratio = 4.0; p.kneeDb = 6.0; p.rangeDb = 18.0;
        p.attackMs = 3.0; p.releaseMs = 40.0; p.detector = det; p.link = link; p.coeffUpdatePeriod = 1;
        DB b; test::run (b.prepare (kSr, 2)); b.setParams (p);

        dynamics::GainComputer gc;
        gc.setMode (mode == Mode::BoostWhenQuiet ? dynamics::Mode::UpCompress : dynamics::Mode::DownCompress);
        gc.setThresholdDb (p.thresholdDb); gc.setRatio (p.ratio); gc.setKneeDb (p.kneeDb); gc.setRangeDb (p.rangeDb);
        const float sign = mode == Mode::BoostWhenLoud ? -1.0f : 1.0f;
        dynamics::GainReductionFollower gr; gr.prepare (kSr); gr.setTimes (p.attackMs, p.releaseMs);

        int mismatches = 0, pauses = 0; double drove = 0.0;
        for (int i = 0; i < 9600; ++i)
        {
            const bool on    = i < 2400 || (i >= 4200 && i < 6400);
            const bool pause = (i >= 2400 && i < 3000) || (i >= 6400 && i < 6800);
            float l = on ? (float) (0.5 * std::sin (2.0 * core::kPi * 5000.0 * i / kSr)) : 0.0f;
            float r = on ? (float) (0.2 * std::sin (2.0 * core::kPi * 6300.0 * i / kSr)) : 0.0f;
            float* io[2] { &l, &r };
            test::run (b.process (io, pause ? 0 : 2, 1));
            pauses += pause ? 1 : 0;

            const float smooth = gr.process (sign * (float) gc.deltaDb (b.detectorLevelDb()));
            gr.flushDenormals();
            if (! core::exactlyEqual (b.dynamicDeltaDb(), (double) smooth)) ++mismatches;
            const double mag = std::fabs (b.dynamicDeltaDb());
            if (! (mag <= drove)) drove = mag;
        }
        test::ok (mismatches == 0 && pauses == 1000 && drove > 6.0 && drove <= 18.0,
                  std::string (name (mode)) + " / " + name (det) + " / " + name (link) + ": 9600 samples, 1000 of them paused, "
                  + std::to_string (mismatches) + " mismatches, delta reached " + std::to_string (drove) + " dB");
    }

    test::group ("a sine in the band reads its band-pass level; silence lets it fall to the floor");
    for (Det det : { Det::Rms, Det::Peak })
    for (double f : { 5000.0, 6500.0 })
    {
        DP p; p.freq = 5000.0; p.Q = 2.0; p.thresholdDb = 12.0; p.detector = det;    // the curve never acts: the readout is the detector's alone
        DB b; test::run (b.prepare (kSr, 1)); b.setParams (p);
        const double a = 0.5, want = steadyLevelDb (det, a * bandPassGain (f, 5000.0, 2.0));
        int phase = 0; runTone (b, f, a, 2400, phase);                              // 50 ms: the 2 ms window has long settled
        double worst = 0.0;
        for (int blk = 0; blk < 40; ++blk)                                          // then read it at 40 block ends
        {
            runTone (b, f, a, 64, phase);
            const double e = std::fabs (b.detectorLevelDb() - want);
            if (! (e <= worst)) worst = e;
        }
        // 0.1 dB is about twice the window's worst ripple: the 2 ms one-pole passes ~0.9 % of the 2f term
        // at 5 kHz, ~0.04 dB on the RMS level and, with the rectifier's higher harmonics, ~0.05 dB on the
        // peak follower.
        test::ok (worst < 0.1, std::string (name (det)) + " " + std::to_string ((int) f) + " Hz: reads "
                  + std::to_string (want) + " dB within " + std::to_string (worst) + " dB (bound 0.1)");

        const double r0 = b.detectorLevelDb();
        std::vector<float> z (48, 0.0f); float* io[1] { z.data() };
        for (int k = 0; k < 10; ++k) { std::fill (z.begin(), z.end(), 0.0f); test::run (b.process (io, 1, 48)); }
        const double fell = r0 - b.detectorLevelDb(), envFall = -480.0 * releaseDbPerSample (det);
        // The band-pass rings for a few samples after the tone stops (its decay is Q / (pi f), ~6 samples),
        // so the fall trails the bare envelope's by at most a few samples' worth of release.
        test::ok (fell > envFall - 1.0 && fell < envFall + 0.1, std::string (name (det)) + " " + std::to_string ((int) f)
                  + " Hz: 10 ms of zeros drops it " + std::to_string (fell) + " dB (the envelope alone: " + std::to_string (envFall) + ")");
        for (int k = 0; k < 200; ++k) { std::fill (z.begin(), z.end(), 0.0f); test::run (b.process (io, 1, 48)); }
        test::ok (core::exactlyEqual (b.detectorLevelDb(), cold), std::string (name (det)) + " " + std::to_string ((int) f)
                  + " Hz: 200 ms of zeros reach the cold floor");
    }

    test::group ("a pause decays it along the envelope's release; reset() returns it to the cold floor");
    {
        std::printf ("    cold floor: %.17g dB\n", cold);
        test::ok (core::exactlyEqual (cold, -180.0), "a cold band reads exactly -180 dB (the 1e-9 floor, narrowed to float)");
        DB once; test::run (once.prepare (kSr, 1));
        float zero = 0.0f; float* zio[1] { &zero };
        test::run (once.process (zio, 1, 1));
        test::ok (core::exactlyEqual (once.detectorLevelDb(), cold), "the floor is what one sample of digital silence hands the curve");

        for (Det det : { Det::Rms, Det::Peak })
        {
            DP p; p.freq = 5000.0; p.Q = 2.0; p.thresholdDb = -30.0; p.ratio = 4.0; p.detector = det;
            DB b; test::run (b.prepare (kSr, 1)); b.setParams (p);
            int phase = 0; runTone (b, 5000.0, 0.5, 4800, phase);
            const double r0 = b.detectorLevelDb();
            std::vector<float> x (64, 0.0f); float* io[1] { x.data() };
            test::run (b.process (io, 0, 480));
            test::approx (b.detectorLevelDb(), r0 + 480.0 * releaseDbPerSample (det), 0.01,
                          std::string (name (det)) + ": 480 paused samples decay it by the envelope's release, from " + std::to_string (r0) + " dB");
            test::run (b.process (io, 0, 48000));
            test::ok (core::exactlyEqual (b.detectorLevelDb(), cold), std::string (name (det)) + ": a 1 s pause reaches the cold floor");

            runTone (b, 5000.0, 0.5, 4800, phase);
            const double charged = b.detectorLevelDb();
            b.reset();
            test::ok (charged > -20.0 && core::exactlyEqual (b.detectorLevelDb(), cold),
                      std::string (name (det)) + ": reset() takes it from " + std::to_string (charged) + " dB to the cold floor");
            runTone (b, 5000.0, 0.5, 4800, phase);
            test::run (b.prepare (kSr, 1));
            test::ok (core::exactlyEqual (b.detectorLevelDb(), cold), std::string (name (det)) + ": so does a re-prepare");
        }
    }

    // The sidechain input is gated (non-finite -> 0, |x| clamped to 1e6) BEFORE the band-pass, so the
    // envelope stays bounded and the dB conversion never sees inf; a NaN could not reach it either. The
    // grid spans the Q clamp (1e-3) to a resonance far above any musical one.
    test::group ("hostile input never makes the readout non-finite");
    {
        const float inf = std::numeric_limits<float>::infinity(), qnan = std::numeric_limits<float>::quiet_NaN();
        const float hostile[] { qnan, inf, -inf, std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(), 1.0e30f, -1.0e30f, 1.0e-40f, 0.0f };
        bool finite = true; double loudest = -1.0e300;
        for (double Q : { 1.0e-3, 0.7, 40.0, 1.0e4 })
        for (Det det : { Det::Rms, Det::Peak })
        for (Link link : { Link::Max, Link::MeanPower })
        {
            DP p; p.freq = 5000.0; p.Q = Q; p.detector = det; p.link = link; p.thresholdDb = -30.0; p.ratio = 8.0; p.rangeDb = 24.0;
            DB b; test::run (b.prepare (kSr, 2)); b.setParams (p);
            std::vector<float> l (32), r (32);
            for (int blk = 0; blk < 100; ++blk)
            {
                for (int j = 0; j < 32; ++j)
                {
                    l[(size_t) j] = hostile[(size_t) ((blk * 32 + j) * 7 % 9)];
                    r[(size_t) j] = (j & 1) ? std::numeric_limits<float>::max() : -std::numeric_limits<float>::max();
                }
                float* io[2] { l.data(), r.data() };
                test::run (b.process (io, blk % 10 == 9 ? 0 : 2, 32));
                const double v = b.detectorLevelDb();
                if (! std::isfinite (v)) finite = false;
                else if (v > loudest) loudest = v;
            }
        }
        std::printf ("    loudest readout over the grid: %.3f dB\n", loudest);
        test::ok (finite && loudest < 200.0, "finite after every block, 32 rows (loudest " + std::to_string (loudest) + " dB)");
    }

    test::group ("no allocation across process(), a pause and the getter");
    {
        DP p; p.freq = 5000.0; p.thresholdDb = -30.0; p.ratio = 4.0;
        DB b; test::run (b.prepare (kSr, 2)); b.setParams (p);
        std::vector<float> l (512), r (512);
        for (int i = 0; i < 512; ++i) { l[(size_t) i] = (float) (0.5 * std::sin (2.0 * core::kPi * 5000.0 * i / kSr)); r[(size_t) i] = 0.5f * l[(size_t) i]; }
        float* io[2] { l.data(), r.data() };
        test::run (b.process (io, 2, 512));
        const long long before = alloc::count.load();
        test::run (b.process (io, 2, 512));
        test::run (b.process (io, 0, 512));
        volatile double sink = b.detectorLevelDb();
        const long long allocated = alloc::count.load() - before;   // read BEFORE the message string exists
        (void) sink;
        test::okNoAlloc (allocated == 0, "process(), a pause and detectorLevelDb() did not allocate");
    }

    return test::report();
}

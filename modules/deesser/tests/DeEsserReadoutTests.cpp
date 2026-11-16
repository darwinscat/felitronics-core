// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// JUCE-free self-tests for DeEsser::detectorLevelDb(), the level the de-esser's ACTIVE detector sees — the
// local split-band one, or the dynamic-EQ delegate's.
//   1. Fed to an INDEPENDENT GainComputer + GainReductionFollower, the readout reproduces gainReductionDb()
//      bit for bit, sample by sample, across audio and pauses, on every route.
//   2. A sine reads the analytic band-pass level through the envelope's steady state; silence lets it fall
//      to the floor.
//   3. A pause decays it along the envelope's own release; reset(), a re-prepare and a route switch return
//      it to the route's cold floor; a setParams() that does not switch leaves it alone.
//   4. Hostile input (NaN, inf, FLT_MAX) never makes it non-finite.
//   5. process(), a pause and the getter allocate nothing.

#include <felitronics_test.h>
#include <alloc_counter.h>   // installs the allocation counter: EVERY form of `new`, over-aligned included
#include <felitronics/core/Math.h>
#include <felitronics/deesser/DeEsser.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <string>
#include <vector>

using namespace felitronics;
using DE   = deesser::DeEsser;
using DP   = deesser::DeEsserParams;
using Mode = deesser::DeEsserMode;
using Det  = dynamics::Detector;
using Link = dynamics::LinkMode;

namespace
{
constexpr double kSr = 48000.0;

struct Route { const char* name; Mode mode; bool listen; bool delegated; };
constexpr Route kRoutes[] {
    { "split-band",          Mode::SplitBand, false, false },
    { "split-band + listen", Mode::SplitBand, true,  false },
    { "dynamic-eq",          Mode::DynamicEq, false, true  },   // the DynamicEqBand delegate's detector
    { "dynamic-eq + listen", Mode::DynamicEq, true,  false },   // listen takes the local detector
};

// |H| of the sidechain band-pass, from the textbook rather than from eq::Svf: the analog prototype
// k s / (s^2 + k s + 1), k = 1/Q, at the bilinear transform's prewarped frequency ratio.
double bandPassGain (double f, double fc, double Q)
{
    const double w = std::tan (core::kPi * f / kSr) / std::tan (core::kPi * fc / kSr), k = 1.0 / Q;
    return k * w / std::sqrt ((1.0 - w * w) * (1.0 - w * w) + (k * w) * (k * w));
}

// The steady envelope of a sine of amplitude a through the symmetric 2 ms level window both detectors
// use: RMS settles on a / sqrt 2, the rectifying peak follower on the rectified mean 2a / pi.
double steadyLevelDb (Det d, double a) { return 20.0 * std::log10 (d == Det::Rms ? a / std::sqrt (2.0) : 2.0 * a / core::kPi); }

// dB the readout loses per sample of zero input: the state is multiplied by the window's release
// coefficient, and an RMS level is the square root of that state.
double releaseDbPerSample (Det d) { return (d == Det::Rms ? 10.0 : 20.0) * std::log10 (std::exp (-1.0 / (0.002 * kSr))); }

const char* name (Det d)  { return d == Det::Rms ? "rms" : "peak"; }
const char* name (Link l) { return l == Link::Max ? "max" : "meanpower"; }

DP routeParams (const Route& r) { DP p; p.mode = r.mode; p.listen = r.listen; return p; }

// What a prepared de-esser on route r that has processed nothing reads.
double coldReadout (const Route& r)
{
    DE d; test::run (d.prepare (kSr, 512, 2)); d.setParams (routeParams (r));
    return d.detectorLevelDb();
}

void runTone (DE& d, double f, double a, int n, int& phase)
{
    std::vector<float> x (64);
    for (int o = 0; o < n; o += 64)
    {
        const int len = std::min (64, n - o);
        for (int j = 0; j < len; ++j, ++phase) x[(size_t) j] = (float) (a * std::sin (2.0 * core::kPi * f * phase / kSr));
        float* io[1] { x.data() };
        test::run (d.process (io, 1, len));
    }
}
} // namespace

int main()
{
    std::printf ("felitronics::deesser detector readout tests\n");
    const double coldLocal = coldReadout (kRoutes[0]), coldDelegate = coldReadout (kRoutes[2]);

    // WHAT IS COMPARED: the reference takes ONLY the readout, runs it through its own GainComputer (the
    // de-esser's curve) and its own GainReductionFollower (same times, flushed at the same block ends —
    // a block is one sample here), and must land on the same number as gainReductionDb(), no tolerance.
    // The local route publishes its gain every sample. The delegate publishes at control rate, every
    // 16th sample from a reset, so its expected meter is the reference sampled on that grid and held.
    // Pauses start where the tone stops, so a readout frozen during one hands the reference a level
    // ~20 dB too high while the detector's own envelope keeps falling.
    test::group ("readout -> independent curve + ballistics == gainReductionDb(), bit for bit");
    for (const Route& route : kRoutes)
    for (Det det : { Det::Rms, Det::Peak })
    for (Link link : { Link::Max, Link::MeanPower })
    {
        DP p = routeParams (route); p.fc = 7000.0; p.scQ = 2.0; p.thresholdDb = -30.0; p.ratio = 4.0; p.kneeDb = 6.0; p.rangeDb = 12.0;
        p.attackMs = 2.0; p.releaseMs = 40.0; p.detector = det; p.link = link;
        DE d; test::run (d.prepare (kSr, 512, 2)); d.setParams (p);

        dynamics::GainComputer gc; gc.setMode (dynamics::Mode::DownCompress);
        gc.setThresholdDb (p.thresholdDb); gc.setRatio (p.ratio); gc.setKneeDb (p.kneeDb); gc.setRangeDb (p.rangeDb);
        dynamics::GainReductionFollower gr; gr.prepare (kSr); gr.setTimes (p.attackMs, p.releaseMs);

        int mismatches = 0, pauses = 0; double drove = 0.0; float held = 0.0f;
        for (int i = 0; i < 9600; ++i)
        {
            const bool on    = i < 2400 || (i >= 4200 && i < 6400);
            const bool pause = (i >= 2400 && i < 3000) || (i >= 6400 && i < 6800);
            float l = on ? (float) (0.5 * std::sin (2.0 * core::kPi * 7000.0 * i / kSr)) : 0.0f;
            float r = on ? (float) (0.2 * std::sin (2.0 * core::kPi * 8500.0 * i / kSr)) : 0.0f;
            float* io[2] { &l, &r };
            test::run (d.process (io, pause ? 0 : 2, 1));
            pauses += pause ? 1 : 0;

            const float smooth = gr.process ((float) gc.deltaDb (d.detectorLevelDb()));
            gr.flushDenormals();
            if (! route.delegated || i % 16 == 0) held = smooth;
            if (! core::exactlyEqual (d.gainReductionDb(), (double) held)) ++mismatches;
            const double mag = std::fabs (d.gainReductionDb());
            if (! (mag <= drove)) drove = mag;
        }
        test::ok (mismatches == 0 && pauses == 1000 && drove > 6.0 && drove <= 12.0,
                  std::string (route.name) + " / " + name (det) + " / " + name (link) + ": 9600 samples, 1000 of them paused, "
                  + std::to_string (mismatches) + " mismatches, reduction reached " + std::to_string (drove) + " dB");
    }

    test::group ("a sine in the band reads its band-pass level; silence lets it fall to the floor");
    for (const Route& route : { kRoutes[0], kRoutes[2] })
    for (Det det : { Det::Rms, Det::Peak })
    for (double f : { 7000.0, 9000.0 })
    {
        DP p = routeParams (route); p.fc = 7000.0; p.scQ = 2.0; p.thresholdDb = 12.0; p.detector = det;   // the curve never acts
        DE d; test::run (d.prepare (kSr, 512, 1)); d.setParams (p);
        const double a = 0.5, want = steadyLevelDb (det, a * bandPassGain (f, 7000.0, 2.0));
        const std::string row = std::string (route.name) + " / " + name (det) + " / " + std::to_string ((int) f) + " Hz";
        int phase = 0; runTone (d, f, a, 2400, phase);                              // 50 ms: the 2 ms window has long settled
        double worst = 0.0;
        for (int blk = 0; blk < 40; ++blk)
        {
            runTone (d, f, a, 64, phase);
            const double e = std::fabs (d.detectorLevelDb() - want);
            if (! (e <= worst)) worst = e;
        }
        // 0.1 dB is about twice the 2 ms window's worst ripple at 2f (~0.05 dB, on the peak follower).
        test::ok (worst < 0.1, row + ": reads " + std::to_string (want) + " dB within " + std::to_string (worst) + " dB (bound 0.1)");

        const double r0 = d.detectorLevelDb();
        std::vector<float> z (48, 0.0f); float* io[1] { z.data() };
        for (int k = 0; k < 10; ++k) { std::fill (z.begin(), z.end(), 0.0f); test::run (d.process (io, 1, 48)); }
        const double fell = r0 - d.detectorLevelDb(), envFall = -480.0 * releaseDbPerSample (det);
        // The band-pass rings for a few samples after the tone stops, so the fall trails the bare
        // envelope's by at most a few samples' worth of release.
        test::ok (fell > envFall - 1.0 && fell < envFall + 0.1,
                  row + ": 10 ms of zeros drops it " + std::to_string (fell) + " dB (the envelope alone: " + std::to_string (envFall) + ")");
        for (int k = 0; k < 200; ++k) { std::fill (z.begin(), z.end(), 0.0f); test::run (d.process (io, 1, 48)); }
        test::ok (core::exactlyEqual (d.detectorLevelDb(), route.delegated ? coldDelegate : coldLocal), row + ": 200 ms of zeros reach the cold floor");
    }

    test::group ("a pause decays it along the release; reset, re-prepare and a route switch return it to the floor");
    {
        std::printf ("    cold floor: %.17g dB local, %.17g dB delegate\n", coldLocal, coldDelegate);
        // The local stage does not narrow its dB to float before the curve; the delegate does.
        test::approx (coldLocal, -180.00000025, 5.0e-9, "the split-band detector's cold floor is -180.00000025 dB");
        test::ok (core::exactlyEqual (coldDelegate, -180.0), "the dynamic-eq delegate's cold floor is exactly -180 dB");

        for (const Route& route : kRoutes)
        for (Det det : { Det::Rms, Det::Peak })
        {
            const double cold = route.delegated ? coldDelegate : coldLocal;
            const std::string row = std::string (route.name) + " / " + name (det);
            DP p = routeParams (route); p.fc = 7000.0; p.thresholdDb = -30.0; p.ratio = 4.0; p.detector = det;
            DE d; test::run (d.prepare (kSr, 512, 1)); d.setParams (p);
            int phase = 0; runTone (d, 7000.0, 0.5, 4800, phase);
            const double r0 = d.detectorLevelDb();
            std::vector<float> x (64, 0.0f); float* io[1] { x.data() };
            test::run (d.process (io, 0, 480));
            test::approx (d.detectorLevelDb(), r0 + 480.0 * releaseDbPerSample (det), 0.01,
                          row + ": 480 paused samples decay it by the envelope's release, from " + std::to_string (r0) + " dB");
            test::run (d.process (io, 0, 48000));
            test::ok (core::exactlyEqual (d.detectorLevelDb(), cold), row + ": a 1 s pause reaches the cold floor");

            runTone (d, 7000.0, 0.5, 4800, phase);
            const double charged = d.detectorLevelDb();
            DP q = p; q.thresholdDb = -20.0; d.setParams (q);                           // no route change: nothing is handed to the curve
            test::ok (core::exactlyEqual (d.detectorLevelDb(), charged), row + ": a setParams() that keeps the route leaves the readout alone");
            d.reset();
            test::ok (charged > -20.0 && core::exactlyEqual (d.detectorLevelDb(), cold),
                      row + ": reset() takes it from " + std::to_string (charged) + " dB to the cold floor");
            runTone (d, 7000.0, 0.5, 4800, phase);
            test::run (d.prepare (kSr, 512, 1));
            test::ok (core::exactlyEqual (d.detectorLevelDb(), cold), row + ": so does a re-prepare");

            // A route switch resets both detectors, so the readout starts from the NEW route's floor.
            runTone (d, 7000.0, 0.5, 4800, phase);
            DP s = q; s.mode = route.mode == Mode::SplitBand ? Mode::DynamicEq : Mode::SplitBand;
            d.setParams (s);
            const bool nowDelegated = s.mode == Mode::DynamicEq && ! s.listen;
            test::ok (core::exactlyEqual (d.detectorLevelDb(), nowDelegated ? coldDelegate : coldLocal),
                      row + ": switching the mode reads the new route's cold floor");
        }
    }

    // The sidechain input is gated (non-finite -> 0, |x| clamped to 1e6) BEFORE the band-pass on both
    // routes, so the envelope stays bounded and the dB conversion never sees inf.
    test::group ("hostile input never makes the readout non-finite");
    {
        const float inf = std::numeric_limits<float>::infinity(), qnan = std::numeric_limits<float>::quiet_NaN();
        const float hostile[] { qnan, inf, -inf, std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(), 1.0e30f, -1.0e30f, 1.0e-40f, 0.0f };
        bool finite = true; double loudest = -1.0e300; int rows = 0;
        for (const Route& route : kRoutes)
        for (double Q : { 1.0e-3, 2.0, 1.0e4 })
        for (Det det : { Det::Rms, Det::Peak })
        for (Link link : { Link::Max, Link::MeanPower })
        {
            DP p = routeParams (route); p.fc = 7000.0; p.scQ = Q; p.detector = det; p.link = link; p.thresholdDb = -30.0; p.ratio = 8.0; p.rangeDb = 24.0;
            DE d; test::run (d.prepare (kSr, 512, 2)); d.setParams (p); ++rows;
            std::vector<float> l (32), r (32);
            for (int blk = 0; blk < 100; ++blk)
            {
                for (int j = 0; j < 32; ++j)
                {
                    l[(size_t) j] = hostile[(size_t) ((blk * 32 + j) * 7 % 9)];
                    r[(size_t) j] = (j & 1) ? std::numeric_limits<float>::max() : -std::numeric_limits<float>::max();
                }
                float* io[2] { l.data(), r.data() };
                test::run (d.process (io, blk % 10 == 9 ? 0 : 2, 32));
                const double v = d.detectorLevelDb();
                if (! std::isfinite (v)) finite = false;
                else if (v > loudest) loudest = v;
            }
        }
        std::printf ("    loudest readout over the grid: %.3f dB\n", loudest);
        test::ok (finite && loudest < 200.0, "finite after every block, " + std::to_string (rows) + " rows (loudest " + std::to_string (loudest) + " dB)");
    }

    test::group ("no allocation across process(), a pause and the getter");
    for (const Route& route : kRoutes)
    {
        DP p = routeParams (route); p.fc = 7000.0; p.thresholdDb = -30.0; p.ratio = 4.0;
        DE d; test::run (d.prepare (kSr, 512, 2)); d.setParams (p);
        std::vector<float> l (512), r (512);
        for (int i = 0; i < 512; ++i) { l[(size_t) i] = (float) (0.5 * std::sin (2.0 * core::kPi * 7000.0 * i / kSr)); r[(size_t) i] = 0.5f * l[(size_t) i]; }
        float* io[2] { l.data(), r.data() };
        test::run (d.process (io, 2, 512));
        const long long before = alloc::count.load();
        test::run (d.process (io, 2, 512));
        test::run (d.process (io, 0, 512));
        volatile double sink = d.detectorLevelDb();
        const long long allocated = alloc::count.load() - before;   // read BEFORE the message string exists
        (void) sink;
        test::okNoAlloc (allocated == 0, std::string (route.name) + ": process(), a pause and detectorLevelDb() did not allocate");
    }

    return test::report();
}

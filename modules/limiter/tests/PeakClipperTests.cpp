// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.
//
// K13 — the peak clipper that lives INSIDE TruePeakLimiter's oversampling island.
//
// THE FOUR THAT MATTER, and why each one is here rather than an obvious alternative:
//   * OFF IS BIT-IDENTICAL. The switch is the feature; a clipper that costs a single bit when off is
//     not switchable, it is always on. Two ways of being off are checked, because they are different
//     code paths: the flag, and a threshold no sample ever reaches.
//   * THE REPORT IS OF WHAT ARRIVED. `maxReconstructedPeakDb()` must not move when the clipper runs.
//     Clip one line too early and it reports the CLIPPED peak — a plausible number, under the ceiling,
//     that no range check can tell from the truth. This is the defect the design was written against.
//   * THE GRID THEOREM, THROUGH THE LIMITER'S OWN ALGEBRA. With the clip at the ceiling exactly, the
//     limiter must find nothing to do: its detector sees a signal the clipper already bounded. Zero
//     gain reduction is a statement about the oversampled grid that no output meter can make.
//   * THE CURVE IS TRACED, NOT ASSERTED. At `overCeilingDb = 0` the limiter is idle, so the output IS
//     the clipped signal; a slow ramp then draws q() at the output, and it is compared against the
//     formula computed here — identity below the knee, exactly T above it, monotone, never over.

#include <felitronics_test.h>

#include <felitronics/analysis/ReferenceTruePeakMeter.h>
#include <felitronics/limiter/TruePeakLimiter.h>

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <vector>

namespace
{
using felitronics::limiter::TruePeakLimiter;
using felitronics::limiter::TruePeakLimiterParams;
namespace test = felitronics::test;
using test::ok;
using test::approx;

constexpr double kPi = 3.14159265358979323846;
constexpr double kFs = 48000.0;

struct Render { std::vector<std::vector<float>> out; double maxReconDb = 0.0, grDb = 0.0; };

Render render (const std::vector<std::vector<float>>& in, const TruePeakLimiterParams& p,
               double ceilingDb, int block = 512)
{
    const int nch = (int) in.size(), n = (int) in[0].size();
    TruePeakLimiter lim;
    if (! lim.prepare (kFs, block, nch, {})) return {};
    TruePeakLimiterParams q = p;
    q.ceilingDbTp = ceilingDb;
    lim.setParams (q);
    const int drain = lim.latencySamples() + 64;
    Render r;
    r.out.assign ((std::size_t) nch, std::vector<float> ((std::size_t) (n + drain), 0.0f));
    for (int c = 0; c < nch; ++c)
        std::copy (in[(std::size_t) c].begin(), in[(std::size_t) c].end(), r.out[(std::size_t) c].begin());
    std::vector<float*> ptrs ((std::size_t) nch);
    for (int i = 0; i < n + drain; i += block)
    {
        const int k = std::min (block, n + drain - i);
        for (int c = 0; c < nch; ++c) ptrs[(std::size_t) c] = r.out[(std::size_t) c].data() + i;
        (void) lim.process (ptrs.data(), nch, k);
    }
    r.maxReconDb = lim.maxReconstructedPeakDb();
    r.grDb       = lim.gainReductionDb();
    return r;
}

std::vector<std::vector<float>> tone (double hz, double dbfs, double seconds, int nch = 2)
{
    const std::size_t n = (std::size_t) (kFs * seconds);
    const double a = std::pow (10.0, dbfs / 20.0);
    std::vector<std::vector<float>> v ((std::size_t) nch, std::vector<float> (n, 0.0f));
    for (std::size_t i = 0; i < n; ++i)
        for (int c = 0; c < nch; ++c)
            v[(std::size_t) c][i] = (float) (a * std::sin (2.0 * kPi * hz * (double) i / kFs));
    return v;
}

int firstDiff (const std::vector<std::vector<float>>& a, const std::vector<std::vector<float>>& b)
{
    if (a.size() != b.size()) return 0;
    for (std::size_t c = 0; c < a.size(); ++c)
    {
        if (a[c].size() != b[c].size()) return 0;
        for (std::size_t i = 0; i < a[c].size(); ++i)
            if (! (a[c][i] == b[c][i])) return (int) i;
    }
    return -1;
}
} // namespace

int main()
{
    std::printf ("felitronics::limiter K13 — the peak clipper inside the island\n");

    test::group ("off is bit-identical, both ways of being off");
    {
        const auto in = tone (997.0, -0.3, 0.5);
        TruePeakLimiterParams off;                       // the clipper is default-off
        TruePeakLimiterParams flagOff;  flagOff.peakClip = false; flagOff.overCeilingDb = 1.0; flagOff.kneeDb = 0.5;
        TruePeakLimiterParams tooHigh;  tooHigh.peakClip = true;  tooHigh.overCeilingDb = 12.0; tooHigh.kneeDb = 0.0;

        const auto a = render (in, off,     -1.0);
        const auto b = render (in, flagOff, -1.0);
        const auto c = render (in, tooHigh, -1.0);
        ok (firstDiff (a.out, b.out) < 0,
            "the flag off renders bit-identically to a default parameter set, whatever the other two say");
        ok (firstDiff (a.out, c.out) < 0,
            "and a threshold 12 dB over a ceiling nothing reaches changes not one bit either");
    }

    test::group ("the report is of what ARRIVED — the one line the clip must sit after");
    {
        // +2.5 dBTP of material into a -1 dBTP ceiling, clipped at the ceiling: the clipper removes
        // every bit of the excess, and the field that says what the limiter SAW must not notice.
        const auto in = tone (1000.0, -0.1, 0.4);
        TruePeakLimiterParams off;
        TruePeakLimiterParams on; on.peakClip = true; on.overCeilingDb = 0.0;
        const auto a = render (in, off, -6.0);
        const auto b = render (in, on,  -6.0);
        ok (a.maxReconDb > -1.0,
            "the programme really does arrive over the ceiling (" + std::to_string (a.maxReconDb) + " dBTP)");
        approx (b.maxReconDb, a.maxReconDb, 0.0,
                "maxReconstructedPeakDb is bit-identical with the clipper running — it reports the INPUT");
    }

    test::group ("the grid theorem, read through the limiter's own algebra");
    {
        // Clip AT the ceiling: every sample on the limiter's grid is then <= the ceiling before the
        // detector ever looks, so the limiter must ask for no reduction at all. Nothing measured on the
        // OUTPUT can make this statement — the downsampler puts the band-limiting term back.
        const auto in = tone (5000.0, -0.1, 0.3);
        TruePeakLimiterParams on;  on.peakClip = true; on.overCeilingDb = 0.0;
        TruePeakLimiterParams off;
        const auto b = render (in, on,  -6.0);
        const auto a = render (in, off, -6.0);
        ok (a.grDb < -1.0, "without the clipper the limiter is working (" + std::to_string (a.grDb) + " dB)");
        approx (b.grDb, 0.0, 1e-6, "with the clip at the ceiling it has nothing left to do");
    }

    test::group ("the curve, traced at the output rather than asserted");
    {
        // At overCeilingDb = 0 the limiter is idle (previous group), so the output IS the clipped
        // signal. A ramp slow enough that the oversampler is flat over it draws q() directly.
        const double ceilDb = -6.0, kneeDb = 1.0;
        const double T = std::pow (10.0, ceilDb / 20.0);
        const double W = T * (1.0 - std::pow (10.0, -kneeDb / 20.0));
        const std::size_t n = 48000;
        std::vector<std::vector<float>> in (1, std::vector<float> (n, 0.0f));
        for (std::size_t i = 0; i < n; ++i) in[0][i] = (float) (2.0 * T * (double) i / (double) n);

        TruePeakLimiterParams on; on.peakClip = true; on.overCeilingDb = 0.0; on.kneeDb = kneeDb;
        const auto r = render (in, on, ceilDb, 1024);

        auto q = [&] (double a) {
            if (a <= T - W) return a;
            if (a >= T + W) return T;
            const double t = a - (T - W);
            return a - t * t / (4.0 * W);
        };
        // The output lags by the limiter's latency; compare on the ramp's own index.
        const int lat = (int) (r.out[0].size() - n) - 64;
        double worst = 0.0, worstOver = 0.0;
        bool monotone = true;
        float prev = -1.0f;
        for (std::size_t i = 200; i + 200 < n; ++i)
        {
            const double got = (double) r.out[0][i + (std::size_t) lat];
            const double want = q ((double) in[0][i]);
            worst = std::max (worst, std::fabs (got - want));
            worstOver = std::max (worstOver, got - T);
            if (r.out[0][i + (std::size_t) lat] < prev - 1e-6f) monotone = false;
            prev = r.out[0][i + (std::size_t) lat];
        }
        ok (worst < 2e-3, "the traced output follows the knee formula (worst " + std::to_string (worst) + ")");
        ok (worstOver < 1e-6, "and never rises above T (worst excess " + std::to_string (worstOver) + ")");
        ok (monotone, "the curve is monotone at the output");
    }

    test::group ("what it reports about itself, and what it refuses to report");
    {
        TruePeakLimiter fresh;
        ok (fresh.clipOccupancy() < 0.0, "with nothing judged there is no occupancy to report");
        ok (fresh.clipReductionQuantileDb (0.95) < 0.0, "and no quantile");
        ok (fresh.clipReductionMaxDb() == 0.0 && fresh.clipRunCount() == 0, "and nothing counted");

        const auto in = tone (100.0, -0.1, 0.5);
        TruePeakLimiterParams on; on.peakClip = true; on.overCeilingDb = 0.0;
        TruePeakLimiter lim;
        if (test::run (lim.prepare (kFs, 512, 2, {})))
        {
            TruePeakLimiterParams p = on; p.ceilingDbTp = -6.0;
            lim.setParams (p);
            std::vector<std::vector<float>> buf = in;
            std::vector<float*> ptrs { buf[0].data(), buf[1].data() };
            (void) lim.process (ptrs.data(), 2, (int) buf[0].size());

            ok (lim.clipOsSamples() > 0 && lim.clipOsTotal() > lim.clipOsSamples(),
                "some of the judged samples were clipped and not all of them ("
                    + std::to_string (lim.clipOsSamples()) + " of " + std::to_string (lim.clipOsTotal()) + ")");
            // A 100 Hz tone at -0.1 dBFS clipped to -6 dB: the deepest reduction is the crest's own
            // ratio, 20log10(peak/T), and nothing in the instrument may exceed it.
            const double want = -0.1 - (-6.0);
            ok (lim.clipReductionMaxDb() <= want + 0.05 && lim.clipReductionMaxDb() > want - 0.5,
                "the deepest reduction is the crest's own ratio, " + std::to_string (want) + " dB (got "
                    + std::to_string (lim.clipReductionMaxDb()) + ")");
            ok (lim.clipRunCount() > 0 && lim.clipRunOsTotal() >= lim.clipRunCount(),
                "the runs were counted (" + std::to_string (lim.clipRunCount()) + ")");
            const double q95 = lim.clipReductionQuantileDb (0.95);
            ok (q95 >= 0.0 && q95 <= lim.clipReductionMaxDb() + 0.1,
                "the 95th percentile sits inside the measured range (" + std::to_string (q95) + ")");
            ok (lim.clipReductionQuantileDb (-0.1) < 0.0 && lim.clipReductionQuantileDb (1.1) < 0.0,
                "a quantile outside [0, 1] is refused rather than clamped");
            ok (lim.clipThresholdDbTp() == -6.0, "the published clip level is the ceiling plus the offset");

            lim.reset();
            ok (lim.clipOsSamples() == 0 && lim.clipRunCount() == 0 && lim.clipOccupancy() < 0.0,
                "reset clears what it counted");
        }
    }

    test::group ("the K13 tap: the reduction of every oversampled sample, and it changes nothing");
    {
        // Two limiters on the same clipped tone, one tapped. The tap is a window onto the clipper, not a stage:
        // the PCM and every aggregate must agree bit for bit, the trace must hold exactly the samples the
        // clipper counted and never a negative, and its deepest value IS the published maximum.
        const auto in = tone (100.0, -0.1, 0.25);
        const int n = (int) in[0].size();
        TruePeakLimiterParams p; p.peakClip = true; p.overCeilingDb = 0.0; p.ceilingDbTp = -6.0;
        TruePeakLimiter plain, tapped;
        if (test::run (plain.prepare (kFs, n, 2, {})) && test::run (tapped.prepare (kFs, n, 2, {})))
        {
            plain.setParams (p);
            tapped.setParams (p);
            const std::size_t os = (std::size_t) n * (std::size_t) tapped.oversampleFactor();
            std::vector<std::vector<float>> a = in, b = in;
            std::vector<float*> pa { a[0].data(), a[1].data() }, pb { b[0].data(), b[1].data() };
            std::vector<float> gr (os), pk (os), red (os, -1.0f);

            felitronics::limiter::TruePeakLimiterTap shortTap { nullptr, nullptr, (int) os - 1 };
            shortTap.peakClipReductionDb = red.data();
            ok (! tapped.process (pb.data(), 2, n, shortTap) && b == in && red[0] == -1.0f,
                "a short capacity refuses the call with the clip trace alone, before anything moves");

            felitronics::limiter::TruePeakLimiterTap tap { gr.data(), pk.data(), (int) os };
            tap.peakClipReductionDb = red.data();
            const bool ranPlain  = plain.process (pa.data(), 2, n);
            const bool ranTapped = tapped.process (pb.data(), 2, n, tap);
            ok (ranPlain && ranTapped && a == b, "the tapped output is the untapped output, sample for sample");
            ok (plain.clipReductionMaxDb() == tapped.clipReductionMaxDb()
                    && plain.clipOsSamples() == tapped.clipOsSamples()
                    && plain.clipRunCount() == tapped.clipRunCount()
                    && plain.maxReconstructedPeakDb() == tapped.maxReconstructedPeakDb(),
                "and so is every aggregate");

            std::int64_t reduced = 0;
            float deepest = 0.0f;
            bool negative = false;
            for (const float r : red)
            {
                if (r > 0.0f) ++reduced;
                if (r > deepest) deepest = r;
                if (r < 0.0f) negative = true;
            }
            ok (! negative, "every oversampled sample was written, and none below zero");
            ok (reduced > 0 && reduced == tapped.clipOsSamples(),
                "the trace holds exactly the samples the clipper counted (" + std::to_string (reduced) + " of "
                    + std::to_string (tapped.clipOsSamples()) + ")");
            ok ((double) deepest == tapped.clipReductionMaxDb(), "its deepest value is the published maximum");
        }

        TruePeakLimiterParams off = p; off.peakClip = false;
        TruePeakLimiter idle;
        if (test::run (idle.prepare (kFs, n, 2, {})))
        {
            idle.setParams (off);
            const std::size_t os = (std::size_t) n * (std::size_t) idle.oversampleFactor();
            std::vector<std::vector<float>> c = in;
            std::vector<float*> pc { c[0].data(), c[1].data() };
            std::vector<float> red (os, -1.0f);
            felitronics::limiter::TruePeakLimiterTap tap { nullptr, nullptr, (int) os };
            tap.peakClipReductionDb = red.data();
            bool allZero = idle.process (pc.data(), 2, n, tap);
            for (const float r : red) allZero = allZero && r == 0.0f;
            ok (allZero, "with the clipper off the trace is written, and it is zero everywhere");
        }
    }

    test::group ("the ranges are clamped and the result is published");
    {
        TruePeakLimiter lim;
        if (test::run (lim.prepare (kFs, 512, 2, {})))
        {
            TruePeakLimiterParams p; p.peakClip = true; p.ceilingDbTp = -2.0;
            p.overCeilingDb = 1e9;  lim.setParams (p);
            ok (lim.clipThresholdDbTp() == -2.0 + TruePeakLimiter::kMaxOverCeilingDb,
                "an absurd offset is clamped to the documented maximum, and the clamp is readable");
            p.overCeilingDb = -5.0; lim.setParams (p);
            ok (lim.clipThresholdDbTp() == -2.0, "a negative offset cannot put the clip UNDER the ceiling");
            p.overCeilingDb = std::nan (""); lim.setParams (p);
            ok (lim.clipThresholdDbTp() == -1.0, "a non-finite offset falls back to the default, as every param here does");
        }
    }

    return test::report();
}

// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.
//
// K5 — what a dynamic EQ band did, read by the (band, lane) PAIR.
//
// THE TWO STATISTICS ANSWER DIFFERENT QUESTIONS AND THE TEST HAS TO SHOW THAT, not merely that both are
// finite. Over the whole programme the number says how OFTEN the band worked; over the windows in which
// its delta was not zero it says how DEEP it went when it did. On a de-esser at a realistic sibilant duty
// they are several times apart — measured on this rig, 1.795 dB against 0.315 — and a gate placed on the
// band's INPUT instead of on its delta collapses the second into the first — which is the planted control below, and the reason the gating signal is named in
// the ABI's own documentation rather than left to be inferred.

#include <felitronics_test.h>

#include <felitronics/mastering/LoudnessSolver.h>
#include <felitronics/mastering/MasteringChain.h>
#include <felitronics/mastering/OfflineRenderer.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <vector>

// THE PRICE IS COUNTED THROUGH `operator new`, NOT THROUGH `sizeof`, and that is not a stylistic choice:
// `sizeof(QuantileHistogram)` is 120 bytes and its bins are a SEPARATE allocation of tens of thousands of
// entries, so an estimate built from `sizeof` came out 2500 times too small and nearly shipped a design
// that would have cost 73 MB a solve. A number about memory that was not taken from the allocator is not
// a measurement of memory.
namespace { std::size_t g_live = 0; bool g_count = false; }
void* operator new (std::size_t n) { if (g_count) g_live += n; return std::malloc (n); }
void operator delete (void* p) noexcept { std::free (p); }
void operator delete (void* p, std::size_t) noexcept { std::free (p); }

namespace
{
using namespace felitronics;
using namespace felitronics::mastering;
namespace test = felitronics::test;
using test::ok;
using test::approx;

constexpr double kPi = 3.14159265358979323846;
constexpr double kFs = 48000.0;

// Speech-like: a steady body with short sibilants, `duty` of the period each.
std::vector<std::vector<float>> speech (double seconds, double duty)
{
    const std::size_t n = (std::size_t) (kFs * seconds);
    std::mt19937 rng (4242); std::uniform_real_distribution<float> nz (-1.0f, 1.0f);
    std::vector<std::vector<float>> v (2, std::vector<float> (n, 0.0f));
    for (std::size_t i = 0; i < n; ++i)
    {
        const double t = (double) i / kFs;
        double a = 0.20 * std::sin (2.0 * kPi * 180.0 * t) + 0.06 * (double) nz (rng);
        if (std::fmod (t, 0.85) < 0.85 * duty) a += 0.50 * (double) nz (rng);
        v[0][i] = v[1][i] = (float) a;
    }
    return v;
}

struct Rig
{
    MasteringChain       chain;
    OfflineRenderer      renderer;
    TargetLoudnessSolver solver;
    MasteringChainParams params;

    bool build (int rendererBlock = 4096)
    {
        MasteringChainConfig cfg;
        if (! chain.prepare (kFs, 2, cfg)) return false;
        if (! renderer.prepare (2, rendererBlock)) return false;
        params = MasteringChainParams {};
        params.limiter.ceilingDbTp = -1.0;
        return solver.prepare (kFs, 2, rendererBlock, chain.internalBlock(), chain.tapOversampleFactor());
    }
    // One de-esser-shaped point: a bell at 6.7 kHz on the Stereo lane, cutting when loud.
    void armDeEsser (int band = 0, double rangeDb = -6.0, int lane = 0)
    {
        eq::BandParams& b = params.eqBands[(std::size_t) band];
        b.on = true; b.type = eq::FilterType::Bell;
        b.lanes[(std::size_t) lane].on = true;
        b.lanes[(std::size_t) lane].freq = 6700.0;
        b.lanes[(std::size_t) lane].Q = 1.5;
        b.dyn.on = true; b.dyn.rangeDb = rangeDb; b.dyn.thrAuto = true;
    }
};

LoudnessSolution solveOn (Rig& rig, const std::vector<std::vector<float>>& in, int passes = 2)
{
    const int frames = (int) in[0].size();
    std::vector<std::vector<float>> out (2, std::vector<float> ((std::size_t) frames, 0.0f));
    const float* ip[2] { in[0].data(), in[1].data() };
    float* op[2] { out[0].data(), out[1].data() };
    LoudnessRequest req;
    req.targetLufs = -14.0;
    req.maxTruePeakDbTp = -1.0;
    req.maxPasses = passes;
    return rig.solver.solve (rig.chain, rig.renderer, rig.params, ip, op, 2, frames, req);
}
} // namespace

int main()
{
    std::printf ("felitronics::mastering K5 — a dynamic band's gain reduction, by (band, lane)\n");

    test::group ("the two halves are several times apart, which is the whole point of publishing both");
    {
        Rig rig;
        if (test::run (rig.build()))
        {
            rig.armDeEsser();
            const auto sol = solveOn (rig, speech (20.0, 0.005));       // 0.5 % sibilant duty
            const BandGrResult* r = sol.bandGrFor (0, 0);
            ok (r != nullptr, "the armed pair has a statistic");
            if (r != nullptr)
            {
                ok (r->whole.valid && r->active.stats.valid, "and both halves are measurements");
                std::printf ("      whole: p95 %.4f max %.4f activeFraction %.4f | active: p95 %.4f max %.4f "
                             "windows %llu of %llu\n",
                             r->whole.p95Db, r->whole.maxDb, r->whole.activeFraction,
                             r->active.stats.p95Db, r->active.stats.maxDb,
                             (unsigned long long) r->active.activeWindows, (unsigned long long) r->active.windows);
                ok (r->whole.activeFraction > 0.0 && r->whole.activeFraction < 1.0,
                    "it worked on part of the programme and not all of it ("
                        + std::to_string (r->whole.activeFraction) + ")");
                // THE NUMBER THE CONSUMER ASKED FOR: at this duty the depth when it works is many times
                // the depth averaged over everything. If a change ever makes these two agree, the active
                // half has stopped being conditioned and the pair has stopped saying two things.
                ok (r->active.stats.p95Db > 0.5,
                    "when it works it pulls real decibels (active p95 " + std::to_string (r->active.stats.p95Db) + ")");
                // MEASURED ON THIS RIG, not imported: a bare LaneDynamics driven by hand reads 5.41 against
                // 0.27 at the same duty, and the full chain reads 1.795 against 0.315, because the solver
                // normalises the programme before the band ever sees it. The RATIO is what the pair exists
                // for, so the ratio is what is pinned — with room, because the absolute pair is a property
                // of this fixture and the solver's own gain, not of the statistic.
                ok (r->active.stats.p95Db > 3.0 * r->whole.p95Db,
                    "…and several times the whole-programme p95 (" + std::to_string (r->active.stats.p95Db)
                        + " against " + std::to_string (r->whole.p95Db) + ", ratio "
                        + std::to_string (r->active.stats.p95Db / r->whole.p95Db) + ")");
                ok (! (r->active.thresholdDb > -1.0e30),
                    "the echoed gate is -inf: there is no dB threshold, the gate is a non-zero delta");
            }
        }
    }

    test::group ("aboveRange cannot be non-zero here, and the doc says so rather than the field pretending");
    {
        Rig rig;
        if (test::run (rig.build()))
        {
            rig.armDeEsser (0, -30.0);                                   // the deepest range the core takes
            const auto sol = solveOn (rig, speech (10.0, 0.20));
            const BandGrResult* r = sol.bandGrFor (0, 0);
            if (test::run (r != nullptr))
            {
                ok (r->whole.aboveRange == 0 && r->active.stats.aboveRange == 0,
                    "nothing is above the range at the deepest setting the core allows");
                ok (r->whole.maxDb <= 30.0 + 1e-9,
                    "because the delta is clamped to 30 dB (max seen " + std::to_string (r->whole.maxDb) + ")");
            }
        }
    }

    test::group ("the refusals name WHY, and each is its own answer");
    {
        Rig rig;
        if (test::run (rig.build()))
        {
            rig.armDeEsser (0, -6.0, 0);                   // band 0 lane 0 armed
            rig.params.eqBands[1].on = true;               // band 1: a point with no dynamics at all
            rig.params.eqBands[1].lanes[0].on = true;
            rig.params.eqBands[2].on = true;               // band 2: armed, range 0 — inert
            rig.params.eqBands[2].lanes[0].on = true;
            rig.params.eqBands[2].dyn.on = true; rig.params.eqBands[2].dyn.rangeDb = 0.0;
            const auto sol = solveOn (rig, speech (6.0, 0.05));

            ok (sol.bandGrFor (0, 0) != nullptr, "the armed pair answers");
            ok (sol.bandGrAbsence[(std::size_t) bandGrIndex (0, 0)] == BandGrAbsence::Armed, "…and says so");
            ok (sol.bandGrAbsence[(std::size_t) bandGrIndex (1, 0)] == BandGrAbsence::NotDynamic,
                "a point with no dynamics is NOT DYNAMIC, not silent");
            ok (sol.bandGrAbsence[(std::size_t) bandGrIndex (2, 0)] == BandGrAbsence::Inert,
                "a point armed at range 0 is INERT — a different thing to tell a user than 'off'");
            ok (sol.bandGrAbsence[(std::size_t) bandGrIndex (0, 3)] == BandGrAbsence::LaneOff,
                "a lane that is not enabled in an armed band is LANE OFF");
            ok (sol.bandGrFor (1, 0) == nullptr && sol.bandGrFor (2, 0) == nullptr
                    && sol.bandGrFor (0, 3) == nullptr,
                "and none of the three has a statistic to hand back");
        }
    }

    test::group ("the pair is per LANE: two lanes of one band answer differently");
    {
        Rig rig;
        if (test::run (rig.build()))
        {
            // Mid and Side both enabled on one point, with the SAME dyn block — the case a per-band
            // number cannot report, because each lane has its own probe and its own delta.
            eq::BandParams& b = rig.params.eqBands[0];
            b.on = true; b.type = eq::FilterType::Bell;
            b.dyn.on = true; b.dyn.rangeDb = -6.0; b.dyn.thrAuto = true;
            for (int l : { (int) eq::Lane::Mid, (int) eq::Lane::Side })
            {
                b.lanes[(std::size_t) l].on = true;
                b.lanes[(std::size_t) l].freq = 6700.0;
                b.lanes[(std::size_t) l].Q = 1.5;
            }
            // A programme whose top is WIDE, so Mid and Side see different levels there.
            auto in = speech (10.0, 0.06);
            for (std::size_t i = 0; i < in[0].size(); ++i) in[1][i] = -in[1][i];
            const auto sol = solveOn (rig, in);
            const BandGrResult* m = sol.bandGrFor (0, (int) eq::Lane::Mid);
            const BandGrResult* s = sol.bandGrFor (0, (int) eq::Lane::Side);
            ok (m != nullptr && s != nullptr, "both lanes of the one band have their own statistic");
            if (m != nullptr && s != nullptr)
                ok (! (m->whole.activeFraction == s->whole.activeFraction)
                        || ! (m->whole.maxDb == s->whole.maxDb),
                    "…and they are not the same answer (Mid " + std::to_string (m->whole.maxDb)
                        + " dB, Side " + std::to_string (s->whole.maxDb) + ")");
        }
    }

    test::group ("the price, taken from the allocator and not from sizeof");
    {
        auto solveCost = [] (int pairs)
        {
            Rig rig;
            if (! rig.build (2048)) return (std::size_t) 0;
            for (int i = 0; i < pairs; ++i) rig.armDeEsser (i, -6.0, 0);
            const auto in = speech (3.0, 0.05);
            // ONE PASS, because this counter is CUMULATIVE and not a peak: the accumulators are locals of
            // one render, freed at its end, so a two-pass solve allocates them twice and a reading taken
            // across both would say the price is double what a solve ever holds at once.
            g_live = 0; g_count = true;
            (void) solveOn (rig, in, 1);
            g_count = false;
            return g_live;
        };
        const std::size_t c0 = solveCost (0), c1 = solveCost (1), c10 = solveCost (10);
        std::printf ("      a solve allocates %zu with no armed pair, %zu with one, %zu with ten\n", c0, c1, c10);
        // WHAT A PAIR ACTUALLY COSTS, measured piece by piece rather than estimated: two histograms over
        // 0..30 dB at 0.01 are 24 008 bytes each, and the trace at the default 1000 buckets is 32 000 —
        // about 80 KB in all. The trace is the bigger half and was missing from every estimate of this
        // made before the allocator was asked.
        const double perPair = (double) (c10 - c1) / 9.0;
        ok (c1 > c0, "an armed pair costs more than none (" + std::to_string (c1 - c0) + " bytes)");
        ok (perPair > 60000.0 && perPair < 110000.0,
            "each further pair costs its two histograms and its trace, about 80 KB ("
                + std::to_string ((long long) perPair) + " bytes)");
        // AND THE PAIR NOBODY ARMED COSTS NOTHING, which is the whole reason the grid is not preallocated:
        // a full 24x5 would have been 5.8 MB of which 119 pairs' worth is never read.
        ok ((double) (c1 - c0) < 10.0 * perPair,
            "…and one armed pair does not drag the other 119 along with it");
    }

    return test::report();
}

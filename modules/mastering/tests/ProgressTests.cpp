// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#include <felitronics/mastering/DeliveredMastering.h>
#include <felitronics_test.h>

#include <cmath>
#include <cstdio>
#include <cstring>
#include <random>
#include <string>
#include <vector>

using namespace felitronics;
using namespace felitronics::mastering;
using felitronics::test::group;
using felitronics::test::ok;

namespace
{
constexpr double kPi = 3.14159265358979323846;
constexpr int    kNch = 2;
constexpr int    kBlock = 4096;

std::vector<float> programme (double fs, int frames, double peak, bool tone = false)
{
    std::vector<float> v ((std::size_t) (frames * kNch), 0.0f);
    std::mt19937 rng (4242u);
    std::uniform_real_distribution<float> u (-1.0f, 1.0f);
    float lp[kNch] = {};
    for (int i = 0; i < frames; ++i)
    {
        const double t = (double) i / fs;
        for (int c = 0; c < kNch; ++c)
        {
            double x = std::sin (2.0 * kPi * (440.0 + 3.0 * c) * t);
            if (! tone)
            {
                lp[c] = 0.9f * lp[c] + 0.1f * u (rng);
                const double macro = 0.35 + 0.65 * (0.5 + 0.5 * std::sin (2.0 * kPi * 0.2 * t));
                const double beat  = std::exp (-8.0 * std::fmod (t, 0.5));
                x = macro * (0.45 * (std::sin (2.0 * kPi * (110.0 + 3.0 * c) * t) + 0.55 * x) + 0.9 * beat * lp[c]);
                if (i % (int) (0.12 * fs) == 0) x += 1.6 * macro;
            }
            v[(std::size_t) (c * frames + i)] = (float) (peak * x);
        }
    }
    return v;
}

struct Planes
{
    const float* in[core::kMaxChannels] {};
    float* out[core::kMaxChannels] {};
    Planes (const std::vector<float>& i, long long inFrames, std::vector<float>& o, long long outFrames)
    {
        for (int c = 0; c < kNch; ++c) { in[c] = i.data() + c * inFrames; out[c] = o.data() + c * outFrames; }
    }
};

bool sameBits (const std::vector<float>& a, const std::vector<float>& b)
{
    return a.size() == b.size() && std::memcmp (a.data(), b.data(), a.size() * sizeof (float)) == 0;
}
bool sameBits (double a, double b) { return std::memcmp (&a, &b, sizeof a) == 0; }
bool sameRecord (const SolvePassRecord& a, const SolvePassRecord& b)
{
    return sameBits (a.gainDb, b.gainDb) && sameBits (a.ceilingDb, b.ceilingDb) && sameBits (a.integratedLufs, b.integratedLufs)
        && sameBits (a.truePeakDbTp, b.truePeakDbTp) && sameBits (a.plrDb, b.plrDb)
        && sameBits (a.limiterMaxGrDb, b.limiterMaxGrDb) && sameBits (a.loudnessRangeLu, b.loudnessRangeLu)
        && a.violated == b.violated;
}
bool sameSolution (const LoudnessSolution& a, const LoudnessSolution& b)
{
    bool same = a.status == b.status && a.binding == b.binding && a.alsoViolated == b.alsoViolated
             && a.passes == b.passes && a.logCount == b.logCount
             && sameBits (a.preLimiterGainDb, b.preLimiterGainDb) && sameBits (a.ceilingDbTp, b.ceilingDbTp)
             && sameBits (a.measured.integratedLufs, b.measured.integratedLufs)
             && sameBits (a.measured.truePeakDbTp, b.measured.truePeakDbTp)
             && sameBits (a.measured.loudnessRangeLu, b.measured.loudnessRangeLu)
             && sameBits (a.measured.limiter.meanDb, b.measured.limiter.meanDb)
             && sameBits (a.measured.limiter.p95Db, b.measured.limiter.p95Db)
             && sameBits (a.measured.compressor.maxDb, b.measured.compressor.maxDb)
             && sameBits (a.measured.limiterMaxReconstructedPeakDb, b.measured.limiterMaxReconstructedPeakDb)
             && a.limiterTrace.buckets == b.limiterTrace.buckets && a.compressorTrace.buckets == b.compressorTrace.buckets;
    for (int i = 0; same && i < a.logCount; ++i) same = sameRecord (a.log[i], b.log[i]);
    for (int i = 0; same && i < a.limiterTrace.buckets; ++i)
        same = sameBits (a.limiterTrace.bucket[i].maxDb, b.limiterTrace.bucket[i].maxDb)
            && sameBits (a.compressorTrace.bucket[i].meanDb, b.compressorTrace.bucket[i].meanDb);
    return same;
}

struct Seen { ProgressEvent e; bool hasRecord = false; SolvePassRecord rec; };
struct Recorder
{
    std::vector<Seen> seen;
    long long stopAt = -1;
    ProgressCallback callback()
    {
        return { [] (void* ctx, const ProgressEvent& e)
                 {
                     auto& r = *static_cast<Recorder*> (ctx);
                     Seen s; s.e = e; s.hasRecord = e.record != nullptr;
                     if (s.hasRecord) s.rec = *e.record;
                     r.seen.push_back (s);
                     return (long long) r.seen.size() - 1 != r.stopAt;
                 }, this };
    }
};

struct Stage { ProgressStage stage; int pass, maxPasses; std::vector<std::size_t> events; };

template <class FramesOf, class UnitsOf>
std::vector<Stage> checkShape (const std::string& label, const std::vector<Seen>& seen, FramesOf framesOf, UnitsOf unitsOf)
{
    std::vector<Stage> st;
    for (std::size_t i = 0; i < seen.size(); ++i)
    {
        if (seen[i].e.fraction == 0.0) st.push_back ({ seen[i].e.stage, seen[i].e.pass, seen[i].e.maxPasses, {} });
        if (! st.empty()) st.back().events.push_back (i);
    }
    int shape = 0, gaps = 0;
    double worst = 0.0;
    for (const Stage& s : st)
    {
        const bool render = s.stage == ProgressStage::SearchPass || s.stage == ProgressStage::FinalRender;
        const Seen& last = seen[s.events.back()];
        if (last.e.fraction != 1.0 || s.events.size() < 2) ++shape;
        if (render != last.hasRecord) ++shape;
        const double bound = (double) framesOf (s) / 100.0, units = (double) unitsOf (s);
        for (std::size_t k = 0; k < s.events.size(); ++k)
        {
            const Seen& x = seen[s.events[k]];
            if (x.e.stage != s.stage || x.e.pass != s.pass || x.e.maxPasses != s.maxPasses) ++shape;
            if (k + 1 < s.events.size() && x.hasRecord) ++shape;
            if (k == 0) continue;
            const double prev = seen[s.events[k - 1]].e.fraction;
            if (! (x.e.fraction > prev)) ++shape;
            const double gap = (x.e.fraction - prev) * units;
            worst = std::max (worst, gap / (double) framesOf (s));
            if (gap > bound * (1.0 + 1.0e-9)) ++gaps;
        }
    }
    ok (! st.empty() && shape == 0 && st.front().events.front() == 0,
        label + ": every stage opens at exactly 0, rises, closes at exactly 1 — a record on a render's close and nowhere else ("
        + std::to_string (st.size()) + " stages, " + std::to_string (seen.size()) + " events)");
    ok (gaps == 0, label + ": no two events of a stage more than 1 % of the programme apart (widest "
                   + std::to_string (100.0 * worst) + " %)");
    return st;
}

void checkRenders (const std::string& label, const std::vector<Seen>& seen, const std::vector<Stage>& st,
                   const LoudnessSolution& sol, int maxPasses)
{
    int renders = 0, numbering = 0, records = 0;
    for (const Stage& s : st)
    {
        if (s.stage != ProgressStage::SearchPass && s.stage != ProgressStage::FinalRender) continue;
        ++renders;
        const int bound = s.stage == ProgressStage::SearchPass ? maxPasses : maxPasses + 1;
        if (s.pass != renders || s.maxPasses != bound) ++numbering;
        if (s.stage == ProgressStage::FinalRender && &s != &st.back()) ++numbering;
        if (renders > sol.logCount || ! sameRecord (seen[s.events.back()].rec, sol.log[renders - 1])) ++records;
    }
    ok (renders == sol.passes && numbering == 0,
        label + ": renders numbered 1.." + std::to_string (sol.passes) + ", bound maxPasses (+1 for the re-render)");
    ok (records == 0, label + ": each render closes with its log entry, bit for bit");
}

struct Rig
{
    MasteringChain chain; OfflineRenderer renderer; TargetLoudnessSolver solver; MasteringChainParams params;
    bool build (double fs, int rendererBlock = kBlock)
    {
        MasteringChainConfig cfg;
        if (! chain.prepare (fs, kNch, cfg) || ! renderer.prepare (kNch, rendererBlock)) return false;
        params.compressor.thresholdDb = -20.0; params.compressor.ratio = 2.0;
        params.limiter.ceilingDbTp = -1.0;
        return solver.prepare (fs, kNch, rendererBlock, chain.internalBlock(), chain.tapOversampleFactor());
    }
};

struct Case { const char* name; bool tone; double peak; LoudnessRequest req; };
std::vector<Case> cases()
{
    Case music { "music -14 LUFS", false, 0.25, {} };
    music.req.targetLufs = -14.0; music.req.maxTruePeakDbTp = -1.0; music.req.maxPasses = 4;
    Case rerender { "re-render", true, 0.3, {} };
    rerender.req.targetLufs = -3.0; rerender.req.maxTruePeakDbTp = -1.0; rerender.req.maxPasses = 3;
    rerender.req.limiterGr.limitDb = 2.0;
    return { music, rerender };
}
}

static void testTheClock()
{
    group ("ProgressClock: the bound for any cutting of the work, and a stop is for good");
    std::mt19937 rng (7u);
    int badGaps = 0, badPieces = 0, badEnds = 0;
    for (long long units : { 1LL, 7LL, 199LL, 200LL, 201LL, 1000LL, 48000LL, 1234567LL })
    {
        Recorder r;
        ProgressClock clock (r.callback());
        const long long frames = units;
        if (! clock.begin (ProgressStage::LoudnessRange, 0, 0, units, frames)) ++badEnds;
        for (long long done = 0; done < units; )
        {
            const int cap = clock.piece (units - done);
            if (cap < 1 || cap > std::max (1LL, frames / 200)) ++badPieces;
            const int n = std::uniform_int_distribution<int> (1, cap) (rng);
            done += n;
            if (! clock.advance (n)) ++badEnds;
        }
        if (! clock.finish()) ++badEnds;
        for (std::size_t i = 1; i < r.seen.size(); ++i)
            if ((r.seen[i].e.fraction - r.seen[i - 1].e.fraction) * (double) units > (double) frames / 100.0 + 1.0e-9
                && (double) frames >= 100.0) ++badGaps;
        if (r.seen.front().e.fraction != 0.0 || r.seen.back().e.fraction != 1.0) ++badEnds;
    }
    ok (badPieces == 0, "a piece is at least 1 and at most half a percent of the frames");
    ok (badGaps == 0, "random pieces keep every gap within 1 % of the frames, at eight lengths");
    ok (badEnds == 0, "each stage opens at 0 and closes at 1");

    Recorder r; r.stopAt = 2;
    ProgressClock clock (r.callback());
    bool b = clock.begin (ProgressStage::Convert, 0, 0, 1000, 1000);
    int trueAfterStop = 0;
    for (int i = 0; i < 100; ++i) if (clock.advance (10) && clock.stopped()) ++trueAfterStop;
    const bool f = clock.finish();
    const bool b2 = clock.begin (ProgressStage::Convert, 0, 0, 1000, 1000);
    ok (b && ! f && ! b2 && clock.stopped() && r.seen.size() == 3 && trueAfterStop == 0,
        "after a false every call answers false and nothing more is sent — a new stage included");

    ProgressClock silent (ProgressCallback {});
    ok (silent.begin (ProgressStage::Convert, 0, 0, 10, 10) && silent.piece (4096) == 4096 && silent.advance (5)
        && silent.finish() && ! silent.stopped(), "with no callback nothing is cut and nothing stops");
}

static void testSolveReportsAndKeepsTheBits()
{
    group ("solve: the same bits with a callback, and the stages it reports");
    const double fs = 48000.0;
    for (const Case& k : cases())
    {
        const int n = (int) (6.0 * fs);
        const auto in = programme (fs, n, k.peak, k.tone);
        Rig rig;
        if (! test::run (rig.build (fs))) return;
        std::vector<float> plainOut (in.size()), heardOut (in.size());
        Planes pp (in, n, plainOut, n), ph (in, n, heardOut, n);
        const LoudnessSolution plain = rig.solver.solve (rig.chain, rig.renderer, rig.params, pp.in, pp.out, kNch, n, k.req);
        Recorder r;
        const LoudnessSolution heard = rig.solver.solve (rig.chain, rig.renderer, rig.params, ph.in, ph.out, kNch, n, k.req,
                                                         r.callback());
        const std::string label = k.name;
        ok (plain.status != MasteringSolveStatus::Cancelled && plain.passes >= 2,
            label + ": PRECONDITION: a verdict, and the search moved (" + std::to_string (plain.passes) + " renders)");
        ok (sameBits (plainOut, heardOut) && sameSolution (plain, heard), label + ": audio and solution bit for bit");
        const long long D = rig.chain.latencySamples();
        const auto st = checkShape (label, r.seen, [&] (const Stage&) { return (long long) n; },
                                    [&] (const Stage&) { return 2LL * n + D; });
        checkRenders (label, r.seen, st, heard, k.req.maxPasses);
        if (k.tone)
            ok (! st.empty() && st.back().stage == ProgressStage::FinalRender,
                label + ": PRECONDITION: the winner was re-rendered, so the FinalRender stage was reported");
    }
}

static void testAStopAtEveryStageBoundaryAndMiddle()
{
    group ("solve: a stop at the opening, the middle and the close of every stage — then the reference bits");
    const double fs = 48000.0;
    for (const Case& k : cases())
    {
        const int n = (int) (4.0 * fs);
        const auto in = programme (fs, n, k.peak, k.tone);
        Rig rig;
        if (! test::run (rig.build (fs))) return;
        std::vector<float> refOut (in.size()), out (in.size());
        Planes pr (in, n, refOut, n), po (in, n, out, n);
        Recorder ref;
        const LoudnessSolution reference = rig.solver.solve (rig.chain, rig.renderer, rig.params, pr.in, pr.out, kNch, n,
                                                             k.req, ref.callback());
        std::vector<Stage> st;
        for (std::size_t i = 0; i < ref.seen.size(); ++i)
        {
            if (ref.seen[i].e.fraction == 0.0) st.push_back ({ ref.seen[i].e.stage, ref.seen[i].e.pass, 0, {} });
            st.back().events.push_back (i);
        }
        int stops = 0, notCancelled = 0, sentAfter = 0, wrongPasses = 0, notBack = 0;
        for (std::size_t s = 0; s < st.size(); ++s)
            for (std::size_t at : { st[s].events.front(), st[s].events[st[s].events.size() / 2], st[s].events.back() })
            {
                Recorder r; r.stopAt = (long long) at;
                const LoudnessSolution stopped = rig.solver.solve (rig.chain, rig.renderer, rig.params, po.in, po.out, kNch,
                                                                   n, k.req, r.callback());
                ++stops;
                if (stopped.status != MasteringSolveStatus::Cancelled) ++notCancelled;
                if (r.seen.size() != at + 1) ++sentAfter;
                const int expected = st[s].pass - (at == st[s].events.front() ? 1 : 0);
                if (stopped.passes != expected) ++wrongPasses;
                std::fill (out.begin(), out.end(), 0.5f);
                const LoudnessSolution again = rig.solver.solve (rig.chain, rig.renderer, rig.params, po.in, po.out, kNch,
                                                                 n, k.req);
                if (! sameBits (out, refOut) || ! sameSolution (again, reference)) ++notBack;
            }
        const std::string label = k.name;
        ok (stops >= 9, label + ": PRECONDITION: " + std::to_string (stops) + " stops over " + std::to_string (st.size()) + " stages");
        ok (notCancelled == 0, label + ": every stop answers Cancelled");
        ok (sentAfter == 0, label + ": and nothing is sent after it");
        ok (wrongPasses == 0, label + ": `passes` counts the renders begun — none for a stop on a render's opening event");
        ok (notBack == 0, label + ": the same objects then solve to the reference bits, every time");
    }
}

static void testShortProgrammesCutTheRendererBlock()
{
    group ("a programme whose 1 % is shorter than the renderer's block: cut, and the same bits");
    const double fs = 48000.0;
    for (int n : { 3000, 30000, 150000 })
    {
        const auto in = programme (fs, n, 0.25);
        Rig rig;
        if (! test::run (rig.build (fs))) return;
        LoudnessRequest req; req.targetLufs = -14.0; req.maxTruePeakDbTp = -1.0;
        std::vector<float> a (in.size()), b (in.size());
        Planes pa (in, n, a, n), pb (in, n, b, n);
        const LoudnessSolution plain = rig.solver.solve (rig.chain, rig.renderer, rig.params, pa.in, pa.out, kNch, n, req);
        Recorder r;
        const LoudnessSolution heard = rig.solver.solve (rig.chain, rig.renderer, rig.params, pb.in, pb.out, kNch, n, req,
                                                         r.callback());
        const std::string label = std::to_string (n) + " frames";
        ok (plain.passes >= 1 && sameBits (a, b) && sameSolution (plain, heard), label + ": bit for bit");
        const long long D = rig.chain.latencySamples();
        checkShape (label, r.seen, [&] (const Stage&) { return (long long) n; }, [&] (const Stage&) { return 2LL * n + D; });
    }
}

static void testTheRange()
{
    group ("measureInputLoudnessRange: one stage, the same number, a stop anywhere");
    const double fs = 48000.0;
    const int n = (int) (12.0 * fs);
    const auto in = programme (fs, n, 0.25);
    Rig rig;
    if (! test::run (rig.build (fs))) return;
    const float* pl[kNch] = { in.data(), in.data() + n };
    double plain = -1.0, heard = -2.0;
    const bool a = rig.solver.measureInputLoudnessRange (pl, kNch, n, plain);
    Recorder r;
    const bool b = rig.solver.measureInputLoudnessRange (pl, kNch, n, heard, r.callback());
    ok (a && b && sameBits (plain, heard), "the same range, bit for bit (" + std::to_string (plain) + " LU)");
    const auto st = checkShape ("range", r.seen, [&] (const Stage&) { return (long long) n; },
                                [&] (const Stage&) { return (long long) n; });
    ok (st.size() == 1 && st[0].stage == ProgressStage::LoudnessRange && st[0].pass == 0 && st[0].maxPasses == 0,
        "one LoudnessRange stage, pass 0 of 0");

    int wrong = 0;
    for (std::size_t at : { std::size_t (0), r.seen.size() / 2, r.seen.size() - 1 })
    {
        Recorder s; s.stopAt = (long long) at;
        double v = 123.0;
        if (rig.solver.measureInputLoudnessRange (pl, kNch, n, v, s.callback()) || v != 123.0 || s.seen.size() != at + 1) ++wrong;
        double again = 0.0;
        if (! rig.solver.measureInputLoudnessRange (pl, kNch, n, again) || ! sameBits (again, plain)) ++wrong;
    }
    ok (wrong == 0, "a stop at the opening, the middle or the close answers false, leaves `out`, sends nothing after, "
                    "and the next measurement is the same bits");
}

static void testAStopEndsTheWorkThere()
{
    group ("a stop in the middle of a stage: the work after it is not done");
    {
        const double a = 44100.0, b = 48000.0;
        const long long n = (long long) (5.0 * a);
        const auto in = programme (a, (int) n, 0.25);
        const long long d = DeliveryConverter::deliveredFrames (a, b, n);
        DeliveryConverter conv;
        if (! test::run (conv.prepare (a, b, kNch, kBlock))) return;
        std::vector<float> out ((std::size_t) (d * kNch));
        Planes pl (in, n, out, d);
        Recorder all;
        ProgressClock full (all.callback());
        ok (conv.convert (pl.in, kNch, n, pl.out, d, full) && all.seen.size() > 10, "PRECONDITION: a converted programme");
        const long long mid = (long long) all.seen.size() / 2;
        std::fill (out.begin(), out.end(), 7.0f);
        Recorder r; r.stopAt = mid;
        ProgressClock clock (r.callback());
        const bool converted = conv.convert (pl.in, kNch, n, pl.out, d, clock);
        long long untouched = 0;
        for (long long i = d * 3 / 4; i < d; ++i) untouched += (out[(std::size_t) i] == 7.0f) ? 1 : 0;
        ok (! converted && clock.stopped() && (long long) r.seen.size() == mid + 1 && untouched == d - d * 3 / 4,
            "convert: stopped at fraction " + std::to_string (all.seen[(std::size_t) mid].e.fraction)
            + ", the last quarter of `out` is not written (" + std::to_string (untouched) + " of " + std::to_string (d - d * 3 / 4) + ")");
    }
    {
        const double fs = 48000.0;
        const int n = (int) (4.0 * fs);
        const auto in = programme (fs, n, 0.25);
        Rig rig;
        if (! test::run (rig.build (fs))) return;
        LoudnessRequest req; req.targetLufs = -14.0; req.maxTruePeakDbTp = -1.0;
        std::vector<float> out (in.size());
        Planes pl (in, n, out, n);
        Recorder all;
        (void) rig.solver.solve (rig.chain, rig.renderer, rig.params, pl.in, pl.out, kNch, n, req, all.callback());
        long long at = -1;
        for (std::size_t i = 0; i < all.seen.size() && at < 0; ++i)
            if (all.seen[i].e.pass == 1 && all.seen[i].e.fraction >= 0.2) at = (long long) i;
        std::fill (out.begin(), out.end(), 7.0f);
        Recorder r; r.stopAt = at;
        const LoudnessSolution stopped = rig.solver.solve (rig.chain, rig.renderer, rig.params, pl.in, pl.out, kNch, n, req, r.callback());
        long long untouched = 0;
        for (int c = 0; c < kNch; ++c)
            for (int i = n * 3 / 4; i < n; ++i) untouched += (pl.out[c][i] == 7.0f) ? 1 : 0;
        ok (at > 0 && stopped.status == MasteringSolveStatus::Cancelled && untouched == (long long) kNch * (n - n * 3 / 4),
            "solve: stopped in the first render at fraction " + std::to_string (at > 0 ? all.seen[(std::size_t) at].e.fraction : -1.0)
            + ", the last quarter of `out` is not written");
    }
}

static void testDelivered()
{
    group ("DeliveredMastering 44.1 -> 48 kHz: the conversion is a stage, the bits are kept, a stop there moves nothing");
    const double a = 44100.0, b = 48000.0;
    const long long n = (long long) (5.0 * a);
    const auto in = programme (a, (int) n, 0.25);
    const long long d = DeliveredMastering::deliveredFrames (a, b, n);
    MasteringChain chain; OfflineRenderer renderer; TargetLoudnessSolver solver; DeliveredMastering dm;
    MasteringChainParams params;
    const bool built = chain.prepare (b, kNch, MasteringChainConfig {}) && renderer.prepare (kNch, kBlock)
                    && solver.prepare (b, kNch, kBlock, chain.internalBlock(), chain.tapOversampleFactor())
                    && dm.prepare (a, b, kNch, kBlock);
    if (! test::run (built)) return;
    LoudnessRequest req; req.targetLufs = -14.0; req.maxTruePeakDbTp = -1.0;
    std::vector<float> plainOut ((std::size_t) (d * kNch)), heardOut (plainOut.size());
    Planes pp (in, n, plainOut, d), ph (in, n, heardOut, d);
    const LoudnessSolution plain = dm.solve (solver, chain, renderer, params, pp.in, kNch, n, pp.out, d, req);
    Recorder r;
    const LoudnessSolution heard = dm.solve (solver, chain, renderer, params, ph.in, kNch, n, ph.out, d, req, r.callback());
    ok (plain.passes >= 1 && sameBits (plainOut, heardOut) && sameSolution (plain, heard), "solve: bit for bit");
    const long long D = chain.latencySamples();
    auto framesOf = [&] (const Stage& s) { return s.stage == ProgressStage::Convert ? n : d; };
    auto unitsOf  = [&] (const Stage& s) { return s.stage == ProgressStage::Convert ? n : 2 * d + D; };
    const auto st = checkShape ("delivered solve", r.seen, framesOf, unitsOf);
    ok (! st.empty() && st[0].stage == ProgressStage::Convert, "solve: the conversion is the first stage");
    checkRenders ("delivered solve", r.seen, st, heard, req.maxPasses);

    int wrong = 0;
    for (std::size_t at : { st[0].events.front(), st[0].events[st[0].events.size() / 2], st[0].events.back() })
    {
        Recorder s; s.stopAt = (long long) at;
        const LoudnessSolution stopped = dm.solve (solver, chain, renderer, params, ph.in, kNch, n, ph.out, d, req, s.callback());
        if (stopped.status != MasteringSolveStatus::Cancelled || stopped.passes != 0 || s.seen.size() != at + 1) ++wrong;
        const LoudnessSolution again = dm.solve (solver, chain, renderer, params, ph.in, kNch, n, ph.out, d, req);
        if (! sameBits (heardOut, plainOut) || ! sameSolution (again, plain)) ++wrong;
    }
    ok (wrong == 0, "solve: a stop while converting is Cancelled with no render begun, and the next solve is the same bits");

    double lraPlain = -1.0, lraHeard = -2.0;
    Recorder q;
    const bool mp = dm.measureInputLoudnessRange (solver, pp.in, kNch, n, lraPlain);
    const bool mh = dm.measureInputLoudnessRange (solver, pp.in, kNch, n, lraHeard, q.callback());
    ok (mp && mh && sameBits (lraPlain, lraHeard), "range: bit for bit");
    const auto qs = checkShape ("delivered range", q.seen, [&] (const Stage& s) { return s.stage == ProgressStage::Convert ? n : d; },
                                [&] (const Stage& s) { return s.stage == ProgressStage::Convert ? n : d; });
    ok (qs.size() == 2 && qs[0].stage == ProgressStage::Convert && qs[1].stage == ProgressStage::LoudnessRange,
        "range: Convert, then LoudnessRange");
}

int main()
{
    std::printf ("felitronics::mastering — progress and cancellation (Progress.h)\n");
    testTheClock();
    testSolveReportsAndKeepsTheBits();
    testAStopAtEveryStageBoundaryAndMiddle();
    testShortProgrammesCutTheRendererBlock();
    testTheRange();
    testAStopEndsTheWorkThere();
    testDelivered();
    return felitronics::test::report();
}

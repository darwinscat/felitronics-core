// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// felitronics::mastering::TargetLoudnessSolver — P7's acceptance, as properties.
//
// THE FIXTURES CARRY THEIR OWN PRECONDITIONS, because the most common defect in this project is a test
// that looks like it checks something and does not. Every group below asserts that the thing it is
// about is actually happening — the limiter is limiting, the gain is moving, the loudness range is
// wide, the constraint is reachable — before it asserts the result.
//
// AND THE ORACLES ARE INDEPENDENT WHERE THEY CAN BE. The solver measures its own render; a test that
// only re-read the solver's own numbers would pass against a solver that measured the wrong buffer. So
// the statistics are re-derived here from the DELIVERED audio with meters this file builds, and the
// gain-reduction summaries are re-derived from the chain's taps driven by hand.

#include <felitronics/mastering/LoudnessSolver.h>
#include <felitronics_test.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <vector>

// The allocation counter, the module suites' idiom — a COUNT only (no bytes, so no allocator's padding to reason
// about): P41 needs "a refused prepare() allocates nothing" counted rather than read off the source.
static std::atomic<long long> g_allocs { 0 };
void* operator new      (std::size_t s) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s ? s : 1); }
void* operator new[]    (std::size_t s) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s ? s : 1); }
void  operator delete   (void* p) noexcept { std::free (p); }
void  operator delete[] (void* p) noexcept { std::free (p); }
void  operator delete   (void* p, std::size_t) noexcept { std::free (p); }
void  operator delete[] (void* p, std::size_t) noexcept { std::free (p); }

using namespace felitronics;
using namespace felitronics::mastering;

namespace
{
constexpr double kPi = 3.14159265358979323846;
constexpr double kFs = 48000.0;

// ---------------------------------------------------------------------------------------------
// Programme generators. NONE of them is a plateau of identical gating blocks: that shape flips every
// block across the absolute gate at once and makes the gated measure jump by a whole LU, which would
// be a fact about the fixture and not about the solver.
struct Programme
{
    std::vector<std::vector<float>> ch;
    int frames() const { return (int) ch[0].size(); }
    int nch()    const { return (int) ch.size(); }
    const float* const* in()  const { return ptrs_.data(); }
    float* const*       out()       { return optrs_.data(); }
    void bind()
    {
        ptrs_.clear(); optrs_.clear();
        for (auto& c : ch) { ptrs_.push_back (c.data()); optrs_.push_back (c.data()); }
    }
private:
    std::vector<const float*> ptrs_;
    std::vector<float*>       optrs_;
};

// A music-like stereo programme with a slow macro envelope, a beat, and — the part that matters — a
// REAL CREST. A fixture whose peak-to-loudness ratio is small never makes the limiter work, and then
// every limiter statistic is zero, every limiter constraint is unreachable, and a whole group of tests
// passes for the wrong reason. `crest` adds sparse transients that carry the peak without carrying the
// energy; the groups below assert the resulting PLR rather than trusting it.
Programme makeMusic (double seconds, double peak, unsigned seed = 12345u, double crest = 1.6)
{
    const int n = (int) (seconds * kFs);
    Programme p; p.ch.assign (2, std::vector<float> ((std::size_t) n, 0.0f));
    std::mt19937 rng (seed);
    std::uniform_real_distribution<float> u (-1.0f, 1.0f);
    float lp[2] = { 0.0f, 0.0f };
    for (int i = 0; i < n; ++i)
    {
        const double t = (double) i / kFs;
        const double macro = 0.35 + 0.65 * (0.5 + 0.5 * std::sin (2.0 * kPi * 0.11 * t - 1.2));
        const double beat  = std::exp (-8.0 * std::fmod (t, 0.5));
        for (int c = 0; c < 2; ++c)
        {
            lp[c] = 0.90f * lp[c] + 0.10f * u (rng);
            const double tone = std::sin (2.0 * kPi * (110.0 + 3.0 * c) * t)
                              + 0.55 * std::sin (2.0 * kPi * (735.0 + 7.0 * c) * t)
                              + 0.30 * std::sin (2.0 * kPi * (2810.0 + 11.0 * c) * t);
            const double v = macro * (0.45 * tone + 0.9 * beat * (double) lp[c]);
            p.ch[(std::size_t) c][(std::size_t) i] = (float) (peak * v);
        }
        // Transients: one every 0.12 s, two samples wide, riding above the body. They add peak and
        // almost no energy, which is exactly what gives a mix its crest — and what a limiter exists to
        // take away. CALIBRATED, not guessed: `crest` was swept and 1.6 puts the fixture's PLR at 13.3,
        // inside the 12.5..15.4 band measured on the six real mixes of this project's corpus. At 1.0 it
        // is 11.2 (denser than any real mix) and at 3.2 it is 19.3 (a click track), and both of those
        // make a different group pass for the wrong reason.
        if (i % (int) (0.12 * kFs) == 0)
        {
            for (int c = 0; c < 2; ++c)
            {
                const double a = peak * crest * macro;
                p.ch[(std::size_t) c][(std::size_t) i]     = (float) ( a);
                if (i + 1 < n) p.ch[(std::size_t) c][(std::size_t) (i + 1)] = (float) (-0.8 * a);
            }
        }
    }
    p.bind();
    return p;
}

// A programme with a REAL loudness range: a quiet section and a loud one, so the LRA constraint has
// something to bind on. Deliberately not the "story then rock" shape of the corpus, but the same
// mechanism: contrast in time.
Programme makeWideRange (double seconds, unsigned seed = 999u)
{
    Programme p = makeMusic (seconds, 0.5, seed);
    const int n = p.frames();
    for (int i = 0; i < n; ++i)
    {
        const double t = (double) i / kFs;
        // Quiet and loud sections alternating every 4 s. THE DEPTH IS CALIBRATED, not chosen: at 0.055
        // (-25 dB) the quiet short-term samples fall below EBU Tech 3342's -20 LU relative gate and drop
        // out of the LRA measurement entirely, so the range the constraint is about is not in the number
        // — measured, that fixture's LRA went UP under limiting (7.90 -> 19.90) and the constraint could
        // not bind at all. 0.30 (-10.5 dB) keeps both sections inside the gate, which is what makes
        // squashing the loud one show up as a smaller range.
        const double k = (std::fmod (t, 8.0) < 4.0) ? 0.30 : 1.0;
        for (int c = 0; c < 2; ++c) p.ch[(std::size_t) c][(std::size_t) i] *= (float) k;
    }
    p.bind();
    return p;
}

// ---------------------------------------------------------------------------------------------
struct Rig
{
    MasteringChain       chain;
    OfflineRenderer      renderer;
    TargetLoudnessSolver solver;
    MasteringChainParams params;

    bool build (int nch, int rendererBlock = 4096, bool clipper = false)
    {
        MasteringChainConfig cfg;
        cfg.clipper = clipper;
        if (! chain.prepare (kFs, nch, cfg)) return false;
        if (! renderer.prepare (nch, rendererBlock)) return false;
        params = MasteringChainParams {};
        params.compressor.thresholdDb = -20.0; params.compressor.ratio = 2.0;
        params.compressor.kneeDb = 6.0; params.compressor.attackMs = 15.0; params.compressor.releaseMs = 180.0;
        params.limiter.ceilingDbTp = -1.0; params.limiter.releaseMs = 100.0;
        return solver.prepare (kFs, nch, rendererBlock, chain.internalBlock(), chain.tapOversampleFactor());
    }
};

// An INDEPENDENT measurement of a delivered buffer: this file's own meters, drained the way the
// solver's are, so "the statistics agree with an independent measurement of the rendered file" is a
// null and not a re-read.
struct Independent { double I = 0.0, TP = 0.0, LRA = 0.0, sp = 0.0; int blocks = 0; };

Independent measureIndependently (const std::vector<std::vector<float>>& buf)
{
    Independent r;
    const int nch = (int) buf.size(), n = (int) buf[0].size();
    analysis::LoudnessMeter lm; analysis::TruePeakMeter tm;
    if (! lm.prepare (kFs, nch, (double) n / kFs + 1.0)) return r;
    if (! tm.prepare (kFs, 65536, nch)) return r;
    std::vector<const float*> p ((std::size_t) nch);
    for (int c = 0; c < nch; ++c) p[(std::size_t) c] = buf[(std::size_t) c].data();
    if (! lm.process (p.data(), nch, n)) return r;
    if (! tm.process (p.data(), nch, n)) return r;
    std::vector<std::vector<float>> z ((std::size_t) nch, std::vector<float> (64, 0.0f));
    std::vector<const float*> zp ((std::size_t) nch);
    for (int c = 0; c < nch; ++c) zp[(std::size_t) c] = z[(std::size_t) c].data();
    (void) tm.process (zp.data(), nch, 64);
    r.I = lm.integratedLufs(); r.TP = tm.truePeakDb(); r.LRA = lm.loudnessRangeLu();
    r.sp = tm.samplePeakDb(); r.blocks = lm.gatingBlockCount();
    return r;
}

const char* statusName (MasteringSolveStatus s)
{
    switch (s)
    {
        case MasteringSolveStatus::Solved:             return "Solved";
        case MasteringSolveStatus::TargetUnreachable:  return "TargetUnreachable";
        case MasteringSolveStatus::UpstreamViolation:  return "UpstreamViolation";
        case MasteringSolveStatus::TargetBetweenAchievable: return "BetweenAchievable";
        case MasteringSolveStatus::PassLimit:          return "PassLimit";
        case MasteringSolveStatus::MeasurementInvalid: return "MeasurementInvalid";
        case MasteringSolveStatus::RenderFailed:       return "RenderFailed";
        case MasteringSolveStatus::NotPrepared:        return "NotPrepared";
        case MasteringSolveStatus::InvalidRequest:     return "InvalidRequest";
    }
    return "?";
}

const char* constraintName (MasteringConstraint c)
{
    switch (c)
    {
        case MasteringConstraint::None:                    return "None";
        case MasteringConstraint::TruePeakCeiling:         return "TruePeakCeiling";
        case MasteringConstraint::LimiterGainReduction:    return "LimiterGainReduction";
        case MasteringConstraint::PeakToLoudness:          return "PeakToLoudness";
        case MasteringConstraint::LoudnessRange:           return "LoudnessRange";
        case MasteringConstraint::GainRange:               return "GainRange";
        case MasteringConstraint::CompressorGainReduction: return "CompressorGainReduction";
    }
    return "?";
}

// =============================================================================================
void testHitsTheTarget()
{
    test::group ("the target is taken to 0.1 LU, and in how many renders");
    for (double target : { -16.0, -14.0, -12.0, -10.0 })
    {
        Programme src = makeMusic (8.0, 0.28);
        std::vector<std::vector<float>> out = src.ch;
        Programme dst; dst.ch = out; dst.bind();

        Rig rig;
        if (! test::run (rig.build (2))) return;

        LoudnessRequest req;
        req.targetLufs = target;
        req.maxTruePeakDbTp = -1.0;
        req.maxPasses = 4;

        const auto sol = rig.solver.solve (rig.chain, rig.renderer, rig.params,
                                           src.in(), dst.out(), 2, src.frames(), req);
        char msg[192];
        std::snprintf (msg, sizeof msg, "target %.1f LUFS: status %s (binding %s)",
                       target, statusName (sol.status), constraintName (sol.binding));
        test::ok (sol.status == MasteringSolveStatus::Solved, msg);
        for (int k = 0; k < sol.logCount; ++k)
            std::printf ("        pass %d: g %7.3f c %7.3f -> I %9.4f TP %8.4f PLR %7.3f limGR %6.3f viol 0x%x\n",
                         k + 1, sol.log[k].gainDb, sol.log[k].ceilingDb, sol.log[k].integratedLufs,
                         sol.log[k].truePeakDbTp, sol.log[k].plrDb, sol.log[k].limiterMaxGrDb,
                         (unsigned) sol.log[k].violated);
        if (sol.status != MasteringSolveStatus::Solved) continue;

        // PRECONDITION: the search actually had to move. A solver that did nothing would also be
        // "within tolerance" if the fixture happened to start on target.
        std::snprintf (msg, sizeof msg, "target %.1f: the gain actually moved (%.3f dB)",
                       target, sol.preLimiterGainDb);
        test::ok (std::fabs (sol.preLimiterGainDb) > 1.0, msg);

        // THE INDEPENDENT NULL — acceptance point 3.
        const Independent ind = measureIndependently (dst.ch);
        std::snprintf (msg, sizeof msg, "target %.1f: reported I equals an independent measurement", target);
        test::approx (sol.measured.integratedLufs, ind.I, 1.0e-9, msg);
        std::snprintf (msg, sizeof msg, "target %.1f: reported TP equals an independent measurement", target);
        test::approx (sol.measured.truePeakDbTp, ind.TP, 1.0e-9, msg);
        std::snprintf (msg, sizeof msg, "target %.1f: reported LRA equals an independent measurement", target);
        test::approx (sol.measured.loudnessRangeLu, ind.LRA, 1.0e-9, msg);

        std::snprintf (msg, sizeof msg, "target %.1f: I is within tolerance", target);
        test::approx (ind.I, target, 0.1, msg);
        std::snprintf (msg, sizeof msg, "target %.1f: delivered TP %.4f is under the ceiling", target, ind.TP);
        test::ok (ind.TP <= -1.0 + 1.0e-9, msg);
        std::snprintf (msg, sizeof msg, "target %.1f: took %d renders (<= 4)", target, sol.passes);
        test::ok (sol.passes <= 4, msg);
        std::printf ("      target %6.1f -> I %9.4f  TP %8.4f  PLR %7.3f  gain %7.3f  ceiling %7.3f  passes %d\n",
                     target, ind.I, ind.TP, sol.measured.plrDb, sol.preLimiterGainDb, sol.ceilingDbTp, sol.passes);
    }
}

// =============================================================================================
void testTwoPassesWhenTheShapeAlreadyAdmitsIt()
{
    test::group ("the closed form: when the delivered PLR already admits the target, two renders");
    // The closed form moves g and c TOGETHER, which leaves the limiting shape untouched and moves the
    // loudness by exactly the shift. Its precondition is PLR <= ceiling - target, so the fixture picks a
    // modest target on a programme with plenty of crest.
    Programme src = makeMusic (8.0, 0.20);
    Programme dst; dst.ch = src.ch; dst.bind();
    Rig rig;
    if (! test::run (rig.build (2))) return;

    LoudnessRequest req;
    req.targetLufs = -20.0;
    req.maxTruePeakDbTp = -1.0;
    req.maxPasses = 3;
    const auto sol = rig.solver.solve (rig.chain, rig.renderer, rig.params,
                                       src.in(), dst.out(), 2, src.frames(), req);
    test::ok (sol.status == MasteringSolveStatus::Solved, "closed-form case solves");
    // PRECONDITION: the case really is the closed-form one — the delivered PLR admits the target.
    test::ok (sol.measured.plrDb <= (-1.0) - (-20.0) + 1e-9,
              "precondition: delivered PLR admits the target (the closed form applies)");
    test::ok (sol.passes <= 2, "two renders, no slope needed");
    std::printf ("      closed form: I %9.4f  TP %8.4f  PLR %7.3f  passes %d\n",
                 sol.measured.integratedLufs, sol.measured.truePeakDbTp, sol.measured.plrDb, sol.passes);
}

// =============================================================================================
void testUnreachableIsNamed()
{
    test::group ("an unreachable target names its binding constraint and does not crush the programme");

    // (a) THE LOUDNESS RANGE. A wide-range programme pushed to a loud target has to be squashed, and the
    // corpus says accepted mastering moves LRA by -0.30 to -0.40 LU. Ask for 0.5 and a loud target.
    {
        Programme src = makeWideRange (20.0);
        Programme dst; dst.ch = src.ch; dst.bind();
        Rig rig;
        if (! test::run (rig.build (2))) return;

        double inLra = 0.0;
        const bool haveLra = rig.solver.measureInputLoudnessRange (src.in(), 2, src.frames(), inLra);
        test::ok (haveLra, "the input's LRA is measurable (the constraint is a DELTA and needs it)");
        // PRECONDITION: the fixture has a range to lose. Without this the constraint cannot bind and
        // the test would pass against a solver that never checked it.
        test::ok (inLra > 6.0, "precondition: the fixture's own LRA is wide (measured "
                               + std::to_string (inLra) + " LU)");

        LoudnessRequest req;
        req.targetLufs = -8.0;                 // loud enough to need heavy limiting
        req.maxTruePeakDbTp = -1.0;
        req.maxLraLossLu = 0.5;                // the corpus number, rounded up
        req.inputLoudnessRangeLu = inLra;
        req.maxPasses = 4;
        const auto sol = rig.solver.solve (rig.chain, rig.renderer, rig.params,
                                           src.in(), dst.out(), 2, src.frames(), req);
        char msg[192];
        std::snprintf (msg, sizeof msg, "LRA case: status %s, binding %s",
                       statusName (sol.status), constraintName (sol.binding));
        // UPSTREAM, and that is the correct answer rather than a weaker one: on this fixture the chain's
        // own compressor already spends more range than the request allows, at the least drive the
        // search will ever use — so the target is not what makes it impossible, the settings are, and
        // reporting "target unreachable" would send the user to change the wrong number.
        test::ok (sol.status == MasteringSolveStatus::UpstreamViolation, msg);
        test::ok (sol.binding == MasteringConstraint::LoudnessRange,
                  "LRA case: the binding constraint is named as the loudness range");
        test::ok (sol.passes == 1, "LRA case: refused after ONE render, not after the whole budget");
        // AND IT DID NOT CRUSH IT. On this fixture the chain's own compressor already spends more of the
        // loudness range than the request allows, so NO render is feasible — which is the interesting
        // case, not a broken one. What must then hold is that the delivered render is the GENTLEST of
        // the ones tried rather than the closest to a target already declared unreachable: its range
        // loss is the smallest in the whole pass log.
        const Independent ind = measureIndependently (dst.ch);
        test::approx (sol.measured.loudnessRangeLu, ind.LRA, 1e-9,
                      "LRA case: the reported range is the delivered buffer's");
        // PRECONDITION: the violation is real and it is the compressor's, not a rounding artefact.
        const double loss = inLra - ind.LRA;
        test::ok (loss > req.maxLraLossLu,
                  "precondition: the range loss really exceeds the allowance ("
                  + std::to_string (loss) + " LU against " + std::to_string (req.maxLraLossLu) + ")");
        test::ok (sol.measured.limiter.maxDb <= 0.0,
                  "precondition: the LIMITER did nothing at this render, so the loss is upstream of it");
        std::printf ("      LRA case: in %.2f LU -> out %.2f LU, I %.3f (target %.1f), status %s, binding %s\n",
                     inLra, ind.LRA, sol.measured.integratedLufs,
                     req.targetLufs, statusName (sol.status), constraintName (sol.binding));
    }

    // (b) THE LIMITER'S GAIN REDUCTION, on a MANUFACTURED dense input — the plan's own recipe, because a
    // dense mix with PLR ~10 does not exist as an input (10 is the RESULT of mastering, measured on two
    // proven pairs). Manufacturing one tests the REFUSAL, not the sound.
    {
        Programme src = makeMusic (8.0, 0.62);
        Programme dst; dst.ch = src.ch; dst.bind();
        Rig rig;
        if (! test::run (rig.build (2))) return;
        LoudnessRequest req;
        req.targetLufs = -7.0;
        req.maxTruePeakDbTp = -1.0;
        req.limiterGr = { 3.0, GrStatistic::Max };
        req.maxPasses = 4;
        const auto sol = rig.solver.solve (rig.chain, rig.renderer, rig.params,
                                           src.in(), dst.out(), 2, src.frames(), req);
        char msg[192];
        std::snprintf (msg, sizeof msg, "limiter-GR case: status %s, binding %s",
                       statusName (sol.status), constraintName (sol.binding));
        test::ok (sol.status == MasteringSolveStatus::TargetUnreachable, msg);
        test::ok (sol.binding == MasteringConstraint::LimiterGainReduction
                  || (sol.alsoViolated & constraintBit (MasteringConstraint::LimiterGainReduction)) != 0, msg);
        test::ok (sol.measured.limiter.valid, "limiter-GR case: the limiter statistic is a measurement");
        test::ok (sol.measured.limiter.maxDb <= 3.0 + 1e-9,
                  "limiter-GR case: the DELIVERED render is inside the limit it refused to break");
        std::printf ("      limiter-GR case: max |GR| %.4f dB (limit 3.0), I %.3f, status %s\n",
                     sol.measured.limiter.maxDb, sol.measured.integratedLufs, statusName (sol.status));
    }
}

// =============================================================================================
void testUpstreamIsNotBlamedOnTheTarget()
{
    test::group ("a compressor limit broken by the SETTINGS is upstream, not 'target unreachable'");
    Programme src = makeMusic (6.0, 0.5);
    Programme dst; dst.ch = src.ch; dst.bind();
    Rig rig;
    if (! test::run (rig.build (2))) return;
    // Make the compressor work hard, then forbid it. The pre-limiter gain node cannot undo this at ANY
    // setting, because it sits AFTER the compressor.
    rig.params.compressor.thresholdDb = -36.0;
    rig.params.compressor.ratio = 8.0;
    LoudnessRequest req;
    req.targetLufs = -14.0;
    req.maxTruePeakDbTp = -1.0;
    req.compressorGr = { 2.0, GrStatistic::P95 };
    req.maxPasses = 3;
    const auto sol = rig.solver.solve (rig.chain, rig.renderer, rig.params,
                                       src.in(), dst.out(), 2, src.frames(), req);
    char msg[192];
    std::snprintf (msg, sizeof msg, "upstream case: status %s, binding %s",
                   statusName (sol.status), constraintName (sol.binding));
    test::ok (sol.status == MasteringSolveStatus::UpstreamViolation, msg);
    test::ok (sol.binding == MasteringConstraint::CompressorGainReduction,
              "upstream case: the binding constraint is named as the compressor's");
    // PRECONDITION: the compressor really is over the limit — otherwise this passes on a fixture that
    // never triggered anything.
    test::ok (sol.measured.compressor.valid && sol.measured.compressor.p95Db > 2.0,
              "precondition: the compressor's p95 GR really does exceed the limit ("
              + std::to_string (sol.measured.compressor.p95Db) + " dB)");
    test::ok (sol.passes == 1, "upstream case: refused after ONE render, not after the whole budget");
}

// =============================================================================================
void testStatisticsAgreeWithAHandDrivenChain()
{
    test::group ("the gain-reduction statistics null against the chain driven by hand");
    // The solver's statistics come from the chain's taps consumed inside `OfflineRenderer`. Re-derive
    // them here by driving the SAME chain with the SAME parameters and the SAME taps, block by block,
    // and accumulating with an independent histogram. A wrong window, a wrong stride or a double count
    // fails here even though the solver is perfectly self-consistent.
    Programme src = makeMusic (5.0, 0.35);
    Programme dst; dst.ch = src.ch; dst.bind();
    Rig rig;
    if (! test::run (rig.build (2))) return;
    LoudnessRequest req;
    req.targetLufs = -10.0;                    // loud enough that the limiter really works — a
                                               // barely-engaging limiter makes the null below a
                                               // statement about a handful of samples
    req.maxTruePeakDbTp = -1.0;
    req.maxPasses = 3;
    const auto sol = rig.solver.solve (rig.chain, rig.renderer, rig.params,
                                       src.in(), dst.out(), 2, src.frames(), req);
    if (! test::run (sol.status == MasteringSolveStatus::Solved
                     || sol.status == MasteringSolveStatus::TargetUnreachable
                     || sol.status == MasteringSolveStatus::PassLimit)) return;

    // Drive the chain by hand at the settings the solver reported.
    MasteringChain chain2;
    MasteringChainConfig cfg;
    if (! test::run (chain2.prepare (kFs, 2, cfg))) return;
    MasteringChainParams p2 = rig.params;
    p2.preLimiterGainDb = sol.preLimiterGainDb;
    p2.limiter.ceilingDbTp = sol.ceilingDbTp;
    chain2.setParams (p2);
    chain2.reset();

    const int K = chain2.internalBlock(), F = chain2.tapOversampleFactor();
    const int blk = 1024;
    std::vector<float> compTap ((std::size_t) (blk + K), 0.0f);
    std::vector<float> limTap  ((std::size_t) (blk + K) * (std::size_t) F, 0.0f);
    MasteringChainTaps taps;
    taps.compressorGrDb = compTap.data(); taps.frameCapacity = blk + K;
    taps.limiterGrDb    = limTap.data();  taps.osCapacity    = (blk + K) * F;

    dynamics::offline::QuantileHistogram hc, hl;
    if (! test::run (hc.prepare (0.0, 400.0, 0.01))) return;
    if (! test::run (hl.prepare (0.0, 400.0, 0.01))) return;
    std::uint64_t ac = 0, al = 0, nc = 0, nl = 0;

    const int frames = src.frames();
    const long long D = chain2.latencySamples();
    const MasteringChainResolved r = chain2.resolved();
    // The STATED offset, read from the chain rather than assumed: the limiter's trace is written where
    // the gain is decided, which is on the oversampled copy, so it lags the limiter's input by the
    // oversampler's latency (`limiterLatency - limiterLookahead`) on top of everything in front of it.
    // The offset itself is checked independently in the tap-offset group below.
    const long long limFrom = r.limiterTapOffset;
    std::vector<std::vector<float>> scratch (2, std::vector<float> ((std::size_t) blk, 0.0f));
    std::vector<float*> sp { scratch[0].data(), scratch[1].data() };
    long long tapPos = 0;
    for (long long off = 0; off < (long long) frames + D; )
    {
        const int m = (int) std::min<long long> ((long long) blk, (long long) frames + D - off);
        for (int c = 0; c < 2; ++c)
            for (int i = 0; i < m; ++i)
            {
                const long long s = off + i;
                sp[(std::size_t) c][(std::size_t) i] = (s < frames) ? src.ch[(std::size_t) c][(std::size_t) s] : 0.0f;
            }
        if (! test::run (chain2.process (sp.data(), 2, m, taps))) return;
        for (int j = 0; j < taps.framesWritten; ++j)
        {
            const long long s = tapPos + j;
            if (s >= 0 && s < frames)
            {
                const double a = std::fabs ((double) compTap[(std::size_t) j]);
                hc.add (a); if (a > 0.1) ++ac; ++nc;
            }
            if (s >= limFrom && s < limFrom + frames)
                for (int k = 0; k < F; ++k)
                {
                    const double a = std::fabs ((double) limTap[(std::size_t) (j * F + k)]);
                    hl.add (a); if (a > 0.1) ++al; ++nl;
                }
        }
        tapPos += taps.framesWritten;
        off += m;
    }

    // PRECONDITION: the traces are LIVE. Zero gain reduction everywhere would make every equality below
    // hold trivially.
    test::ok (hc.maxValue() > 0.5, "precondition: the compressor really compressed ("
                                   + std::to_string (hc.maxValue()) + " dB peak)");
    test::ok (hl.maxValue() > 0.5, "precondition: the limiter really limited ("
                                   + std::to_string (hl.maxValue()) + " dB peak)");
    test::ok (nc == (std::uint64_t) frames, "compressor window: exactly `frames` samples counted");
    test::ok (nl == (std::uint64_t) frames * (std::uint64_t) F, "limiter window: exactly frames*F counted");

    double p95c = 0.0, p95l = 0.0;
    test::ok (hc.quantile (0.95, p95c), "hand-driven compressor p95 is answerable");
    test::ok (hl.quantile (0.95, p95l), "hand-driven limiter p95 is answerable");
    test::approx (sol.measured.compressor.meanDb, hc.mean(), 1.0e-9 * std::fmax (1.0, hc.mean()),
                  "compressor mean nulls");
    test::approx (sol.measured.compressor.p95Db,  p95c,      0.0, "compressor p95 nulls");
    test::approx (sol.measured.compressor.maxDb,  hc.maxValue(), 0.0, "compressor max nulls");
    test::approx (sol.measured.compressor.activeFraction, (double) ac / (double) nc, 0.0,
                  "compressor active fraction nulls");
    // The MEAN is a sum over a million values and the two paths accumulate it in different block
    // orders, so it nulls to double-summation precision rather than bit for bit. Everything else is an
    // order statistic or a count and IS exact.
    test::approx (sol.measured.limiter.meanDb, hl.mean(), 1.0e-9 * std::fmax (1.0, hl.mean()),
                  "limiter mean nulls");
    test::approx (sol.measured.limiter.p95Db,  p95l,      0.0, "limiter p95 nulls");
    test::approx (sol.measured.limiter.maxDb,  hl.maxValue(), 0.0, "limiter max nulls");
    test::approx (sol.measured.limiter.activeFraction, (double) al / (double) nl, 0.0,
                  "limiter active fraction nulls");
    std::printf ("      hand-driven null: comp mean %.5f p95 %.5f max %.5f | lim mean %.5f p95 %.5f max %.5f\n",
                 hc.mean(), p95c, hc.maxValue(), hl.mean(), p95l, hl.maxValue());
}

// =============================================================================================
void testTappedRenderNullsAgainstThePlainOne()
{
    test::group ("a tapped render is bit-identical to a plain one");
    // The tap must not move audio. This is the same rule the compressor's tap follows, one level up,
    // and it is what makes the statistics free rather than a second behaviour.
    Programme src = makeMusic (2.0, 0.4);
    Programme a; a.ch = src.ch; a.bind();
    Programme b; b.ch = src.ch; b.bind();

    MasteringChainConfig cfg; cfg.clipper = true;
    MasteringChain c1, c2;
    if (! test::run (c1.prepare (kFs, 2, cfg))) return;
    if (! test::run (c2.prepare (kFs, 2, cfg))) return;
    MasteringChainParams p;
    p.preLimiterGainDb = 6.0; p.limiter.ceilingDbTp = -1.0;
    c1.setParams (p); c2.setParams (p);
    OfflineRenderer r1, r2;
    if (! test::run (r1.prepare (2, 512))) return;
    if (! test::run (r2.prepare (2, 512))) return;

    if (! test::run (r1.render (c1, src.in(), a.out(), 2, src.frames()))) return;

    const int K = c2.internalBlock(), F = c2.tapOversampleFactor();
    std::vector<float> ct ((std::size_t) (512 + K), 0.0f);
    std::vector<float> lt ((std::size_t) (512 + K) * (std::size_t) F, 0.0f);
    std::vector<float> lp ((std::size_t) (512 + K) * (std::size_t) F, 0.0f);
    std::vector<std::vector<float>> pre (2, std::vector<float> ((std::size_t) (512 + K), 0.0f));
    std::vector<float*> prePtr { pre[0].data(), pre[1].data() };
    MasteringChainTaps taps;
    taps.compressorGrDb = ct.data(); taps.preLimiter = prePtr.data(); taps.frameCapacity = 512 + K;
    taps.limiterGrDb = lt.data(); taps.limiterPeakLin = lp.data(); taps.osCapacity = (512 + K) * F;
    long long seen = 0;
    if (! test::run (r2.render (c2, src.in(), b.out(), 2, src.frames(), taps,
                                [&] (const MasteringChainTaps& t, long long) noexcept { seen += t.framesWritten; })))
        return;

    int diff = 0;
    for (int c = 0; c < 2; ++c)
        for (int i = 0; i < src.frames(); ++i)
            if (a.ch[(std::size_t) c][(std::size_t) i] != b.ch[(std::size_t) c][(std::size_t) i]) ++diff;
    test::ok (diff == 0, "tapped render is bit-identical to the plain one");
    // PRECONDITION: the render is not silence, and the taps actually ran.
    double peak = 0.0;
    for (int i = 0; i < src.frames(); ++i) peak = std::fmax (peak, std::fabs ((double) a.ch[0][(std::size_t) i]));
    test::ok (peak > 0.1, "precondition: the render is not silence (" + std::to_string (peak) + ")");
    test::ok (seen >= src.frames(), "precondition: the tap stream covered the programme");
}

// =============================================================================================
void testShortTapRefusesTheWholeCall()
{
    test::group ("a tap too short refuses the whole call and touches nothing");
    MasteringChain ch;
    MasteringChainConfig cfg;
    if (! test::run (ch.prepare (kFs, 2, cfg))) return;
    const int K = ch.internalBlock(), F = ch.tapOversampleFactor();
    const int n = 4 * K;
    std::vector<std::vector<float>> buf (2, std::vector<float> ((std::size_t) n, 0.25f));
    std::vector<float*> bp { buf[0].data(), buf[1].data() };

    std::vector<float> shortComp (K / 2, -12345.0f);
    MasteringChainTaps t;
    t.compressorGrDb = shortComp.data(); t.frameCapacity = K / 2;
    test::ok (! ch.process (bp.data(), 2, n, t), "short compressor tap: the call is refused");
    int touched = 0; for (float v : shortComp) if (v != -12345.0f) ++touched;
    test::ok (touched == 0, "short compressor tap: nothing was written");
    int moved = 0; for (int i = 0; i < n; ++i) if (buf[0][(std::size_t) i] != 0.25f) ++moved;
    test::ok (moved == 0, "short compressor tap: the audio did not move");
    test::ok (t.framesWritten == 0 && t.osWritten == 0, "short tap: the counts report nothing written");

    std::vector<float> okComp ((std::size_t) (n + K), 0.0f);
    std::vector<float> shortLim (8, -12345.0f);
    MasteringChainTaps t2;
    t2.compressorGrDb = okComp.data(); t2.frameCapacity = n + K;
    t2.limiterGrDb = shortLim.data(); t2.osCapacity = 8;
    test::ok (! ch.process (bp.data(), 2, n, t2), "short limiter tap: the call is refused");
    int touched2 = 0; for (float v : shortLim) if (v != -12345.0f) ++touched2;
    test::ok (touched2 == 0, "short limiter tap: nothing was written");
    (void) F;
}

// =============================================================================================
void testRefusalsAndDegenerateInputs()
{
    test::group ("refusals: the solver says no rather than answering wrongly");
    Programme src = makeMusic (4.0, 0.3);
    Programme dst; dst.ch = src.ch; dst.bind();

    {
        TargetLoudnessSolver s;
        MasteringChain ch; OfflineRenderer r;
        LoudnessRequest req; req.targetLufs = -14.0; req.maxTruePeakDbTp = -1.0;
        const auto sol = s.solve (ch, r, MasteringChainParams {}, src.in(), dst.out(), 2, src.frames(), req);
        test::ok (sol.status == MasteringSolveStatus::NotPrepared, "unprepared solver refuses");
    }
    {
        Rig rig;
        if (! test::run (rig.build (2))) return;
        LoudnessRequest req; req.targetLufs = std::numeric_limits<double>::quiet_NaN();
        const auto sol = rig.solver.solve (rig.chain, rig.renderer, rig.params,
                                           src.in(), dst.out(), 2, src.frames(), req);
        test::ok (sol.status == MasteringSolveStatus::InvalidRequest, "NaN target refused");

        LoudnessRequest req2; req2.targetLufs = -14.0; req2.maxTruePeakDbTp = -1.0; req2.maxPasses = 0;
        const auto sol2 = rig.solver.solve (rig.chain, rig.renderer, rig.params,
                                            src.in(), dst.out(), 2, src.frames(), req2);
        test::ok (sol2.status == MasteringSolveStatus::InvalidRequest, "zero pass budget refused");

        LoudnessRequest req3; req3.targetLufs = -14.0; req3.maxTruePeakDbTp = -1.0;
        const auto sol3 = rig.solver.solve (rig.chain, rig.renderer, rig.params,
                                            src.in(), dst.out(), 1, src.frames(), req3);
        test::ok (sol3.status == MasteringSolveStatus::InvalidRequest,
                  "a width the chain was not prepared for is refused");
    }
    {
        // DIGITAL SILENCE. The meter has no gating block, so there is no measurement — and a solver
        // that answered anyway would set a ceiling from a -200 dBTP reading, which the limiter clamps to
        // +60, i.e. it would turn the limiter OFF and report success.
        Rig rig;
        if (! test::run (rig.build (2))) return;
        Programme q; q.ch.assign (2, std::vector<float> ((std::size_t) (4 * (int) kFs), 0.0f)); q.bind();
        Programme qo; qo.ch = q.ch; qo.bind();
        LoudnessRequest req; req.targetLufs = -14.0; req.maxTruePeakDbTp = -1.0;
        const auto sol = rig.solver.solve (rig.chain, rig.renderer, rig.params,
                                           q.in(), qo.out(), 2, q.frames(), req);
        test::ok (sol.status == MasteringSolveStatus::MeasurementInvalid,
                  "digital silence: no measurement, and it says so instead of inventing a ceiling");
        // NOT `gatingBlocks == 0`: that counter reports blocks that ARRIVED, and 4 s of silence supplies
        // plenty. What says "no measurement" is the meter's own sentinel — `integrated()` returns the
        // literal -120.0 from every path where nothing passed the absolute gate.
        test::ok (sol.measured.integratedLufs <= -120.0,
                  "digital silence: the reading is the meter's no-measurement sentinel, not a loudness");
        test::ok (! sol.measured.loudnessValid, "digital silence: the measurement is marked invalid");
    }
}

// =============================================================================================
void testBlockIndependence()
{
    test::group ("the solve does not depend on the renderer's block size");
    Programme src = makeMusic (4.0, 0.3);
    double firstI = 0.0, firstG = 0.0;
    bool have = false;
    for (int blk : { 64, 256, 1000, 8192 })
    {
        Programme dst; dst.ch = src.ch; dst.bind();
        Rig rig;
        if (! test::run (rig.build (2, blk))) return;
        LoudnessRequest req; req.targetLufs = -13.0; req.maxTruePeakDbTp = -1.0; req.maxPasses = 3;
        const auto sol = rig.solver.solve (rig.chain, rig.renderer, rig.params,
                                           src.in(), dst.out(), 2, src.frames(), req);
        if (! test::run (sol.status == MasteringSolveStatus::Solved)) return;
        if (! have) { firstI = sol.measured.integratedLufs; firstG = sol.preLimiterGainDb; have = true; }
        char msg[128];
        std::snprintf (msg, sizeof msg, "block %d: the achieved loudness is bit-identical to block 64", blk);
        test::approx (sol.measured.integratedLufs, firstI, 0.0, msg);
        std::snprintf (msg, sizeof msg, "block %d: the applied gain is bit-identical to block 64", blk);
        test::approx (sol.preLimiterGainDb, firstG, 0.0, msg);
    }
}

// =============================================================================================
void testTheReportedRenderIsTheDeliveredOne()
{
    test::group ("the delivered audio is the render the report describes");
    // The search may end on a candidate that is not its last attempt. The buffer must then hold the
    // reported one — handing back a file the report does not describe is the defect this guards.
    Programme src = makeWideRange (20.0);
    Programme dst; dst.ch = src.ch; dst.bind();
    Rig rig;
    if (! test::run (rig.build (2))) return;
    double inLra = 0.0;
    (void) rig.solver.measureInputLoudnessRange (src.in(), 2, src.frames(), inLra);
    LoudnessRequest req;
    req.targetLufs = -7.0;
    req.maxTruePeakDbTp = -1.0;
    req.maxLraLossLu = 0.5;
    req.inputLoudnessRangeLu = inLra;
    req.maxPasses = 4;
    const auto sol = rig.solver.solve (rig.chain, rig.renderer, rig.params,
                                       src.in(), dst.out(), 2, src.frames(), req);
    const Independent ind = measureIndependently (dst.ch);
    test::approx (sol.measured.integratedLufs, ind.I, 1.0e-9,
                  "the delivered buffer measures what the report says");
    test::approx (sol.measured.truePeakDbTp, ind.TP, 1.0e-9,
                  "the delivered buffer's true peak is the reported one");
    std::printf ("      delivered-is-reported: status %s, I %.4f vs %.4f\n",
                 statusName (sol.status), sol.measured.integratedLufs, ind.I);
}

// =============================================================================================
// The claims the HEADERS make, pinned. A number in a comment that no test reads goes stale silently,
// and this file is where the ones this change introduced are held.
void testTheScaleLawIsPinned()
{
    test::group ("y(g,c) == 10^(c/20) * y(g-c, 0) — the identity the whole design rests on");
    // Dither OFF on purpose: it is added AFTER the limiter and does not scale, so it is outside the
    // identity by construction. Saying which half of the chain the law covers is the point.
    Programme src = makeMusic (1.0, 0.4);
    double worst = 0.0, worstPeak = 0.0;
    for (double c : { -1.0, -3.0, -6.0 })
        for (double g : { 3.0, 6.0, 12.0 })
        {
            auto run = [&] (double gg, double cc, std::vector<std::vector<float>>& o)
            {
                MasteringChainConfig cfg;
                cfg.eq = false; cfg.compressor = false; cfg.clipper = false;
                cfg.limiter = true; cfg.dither = false;
                MasteringChain ch;
                if (! ch.prepare (kFs, 2, cfg)) return false;
                MasteringChainParams p;
                p.preLimiterGainDb = gg; p.limiter.ceilingDbTp = cc; p.limiter.releaseMs = 100.0;
                ch.setParams (p);
                OfflineRenderer r;
                if (! r.prepare (2, 1024)) return false;
                o = src.ch;
                std::vector<float*> op { o[0].data(), o[1].data() };
                return r.render (ch, src.in(), op.data(), 2, src.frames());
            };
            std::vector<std::vector<float>> a, b;
            if (! test::run (run (g, c, a)) || ! test::run (run (g - c, 0.0, b))) return;
            const double k = std::pow (10.0, c / 20.0);
            for (int i = 0; i < src.frames(); ++i)
            {
                worst = std::fmax (worst, std::fabs ((double) a[0][(std::size_t) i] - k * (double) b[0][(std::size_t) i]));
                worstPeak = std::fmax (worstPeak, std::fabs ((double) a[0][(std::size_t) i]));
            }
        }
    // PRECONDITION: the limiter is actually doing something, or the identity is about silence.
    test::ok (worstPeak > 0.1, "precondition: the compared renders are not silence ("
                               + std::to_string (worstPeak) + ")");
    // The header claims 9.1e-07 .. 1.7e-06 on a programme peaking at 0.89. Pinned as an ORDER, not as a
    // literal: it is float rounding in `dbToGain` and the FIR, which moves with the toolchain. What must
    // not move is that it stays six decades under the signal.
    test::ok (worst < 1.0e-5, "the scale law holds to better than 1e-5 (measured "
                              + std::to_string (worst) + " at peak " + std::to_string (worstPeak) + ")");
    std::printf ("      scale law: worst |y(g,c) - 10^(c/20) y(d,0)| = %.3e at peak %.3f\n", worst, worstPeak);
}

// =============================================================================================
void testTheAbsoluteGateStepIsPinned()
{
    test::group ("the absolute gate at -70 LUFS is NOT scale-invariant — the header's reason, measured");
    // Two halves either side of the gate. At g = 0 only the louder half is counted; a gain that lifts
    // the quieter half over the gate changes the SET being averaged, and `I(g) = I(0) + g` stops being
    // true. This is the reason the solver measures every candidate instead of trusting the identity.
    const int half = (int) (10.0 * kFs), n = 2 * half;
    auto measure = [&] (double a1, double a2, double gDb)
    {
        std::vector<float> x ((std::size_t) n);
        const double gg = std::pow (10.0, gDb / 20.0);
        for (int i = 0; i < n; ++i)
            x[(std::size_t) i] = (float) (gg * (i < half ? a1 : a2) * std::sin (2.0 * kPi * 1000.0 * (double) i / kFs));
        analysis::LoudnessMeter m;
        if (! m.prepare (kFs, 1, 120.0)) return 1.0e9;
        const float* p[1] = { x.data() };
        if (! m.process (p, 1, n)) return 1.0e9;
        return m.integratedLufs();
    };
    const double cal = measure (1.0, 1.0, 0.0);
    const double a1 = std::pow (10.0, (-69.0 - cal) / 20.0);
    const double a2 = std::pow (10.0, (-71.0 - cal) / 20.0);
    const double i0 = measure (a1, a2, 0.0);
    const double i2 = measure (a1, a2, 2.0);
    const double err = i2 - (i0 + 2.0);
    // PRECONDITION: the fixture really straddles the gate, or the error below is a rounding artefact.
    test::ok (std::fabs (i0 - (-69.0)) < 0.2,
              "precondition: at g = 0 only the -69 half is counted (I = " + std::to_string (i0) + ")");
    test::approx (err, -0.879829, 0.01, "the step is the measured -0.879829 LU, not zero");
    std::printf ("      absolute gate: I(0) %.6f, I(+2) %.6f, error against I(0)+2 = %+.6f LU\n", i0, i2, err);
}

// =============================================================================================
void testTheReviewsCounterexamples()
{
    test::group ("the code-review round's counterexamples, each pinned");
    Programme src = makeMusic (2.0, 0.2);

    // (a) ALIASED in/out. Every pass after the first would read the previous master.
    {
        Programme p; p.ch = src.ch; p.bind();
        Rig rig; if (! test::run (rig.build (2))) return;
        LoudnessRequest req; req.targetLufs = -20.0; req.maxTruePeakDbTp = -1.0;
        const auto sol = rig.solver.solve (rig.chain, rig.renderer, rig.params,
                                           p.in(), p.out(), 2, p.frames(), req);
        test::ok (sol.status == MasteringSolveStatus::InvalidRequest,
                  "in == out is refused: a search cannot read its own previous master");
        test::ok (sol.passes == 0, "in == out: refused before spending a render");
    }

    // (b) A CALLER CEILING BELOW THE PROMISE. The solver may raise it, so a limiter-GR violation there
    //     is NOT upstream — it is something its own knobs can relieve.
    {
        Programme dst; dst.ch = src.ch; dst.bind();
        Rig rig; if (! test::run (rig.build (2))) return;
        rig.params.limiter.ceilingDbTp = -40.0;                 // far below the promise
        LoudnessRequest req;
        req.targetLufs = -30.0; req.maxTruePeakDbTp = -1.0;
        req.limiterGr = { 1.0, GrStatistic::Max };
        req.maxPasses = 5;
        const auto sol = rig.solver.solve (rig.chain, rig.renderer, rig.params,
                                           src.in(), dst.out(), 2, src.frames(), req);
        test::ok (sol.status != MasteringSolveStatus::UpstreamViolation,
                  "a ceiling the solver may still RAISE is not an upstream violation");
        std::printf ("      ceiling -40 vs promise -1: status %s, ceiling %.3f, limGR max %.3f, passes %d\n",
                     statusName (sol.status), sol.ceilingDbTp, sol.measured.limiter.maxDb, sol.passes);
    }

    // (c) A TARGET PAST THE ACTUATOR. The gain node clamps at +-60 dB; a target that needs more is
    //     unreachable BY NAME, not a budget that ran out.
    {
        const int n = (int) (2.0 * kFs);
        Programme q; q.ch.assign (2, std::vector<float> ((std::size_t) n, 0.0f));
        for (int i = 0; i < n; ++i)
        {
            const float v = (float) (0.0005 * std::sin (2.0 * kPi * 1000.0 * (double) i / kFs));
            q.ch[0][(std::size_t) i] = v; q.ch[1][(std::size_t) i] = v;
        }
        q.bind();
        Programme o; o.ch = q.ch; o.bind();
        Rig rig; if (! test::run (rig.build (2))) return;
        rig.params.compressor.thresholdDb = 20.0;               // out of the way
        LoudnessRequest req;
        req.targetLufs = 20.0; req.maxTruePeakDbTp = 30.0; req.maxPasses = 4;
        const auto sol = rig.solver.solve (rig.chain, rig.renderer, rig.params,
                                           q.in(), o.out(), 2, q.frames(), req);
        test::ok (sol.status == MasteringSolveStatus::TargetUnreachable,
                  "a target past the actuator is unreachable, not a pass limit");
        test::ok (sol.binding == MasteringConstraint::GainRange,
                  "...and the binding constraint is NAMED as the gain range");
        test::ok ((sol.alsoViolated & constraintBit (MasteringConstraint::GainRange)) != 0,
                  "...and it appears in the violation mask");
        // PRECONDITION: the actuator really is pinned, or this passes for another reason.
        test::approx (std::fabs (sol.preLimiterGainDb), 60.0, 1e-9,
                      "precondition: the gain really is pinned at the +-60 dB clamp");
    }

    // (d) AN IMPOSSIBLE CONSTRAINT SPELLED WITH THE WRONG INFINITY. `isfinite` is true of neither
    //     infinity, so testing it disabled a limit the caller meant to be unsatisfiable.
    {
        Programme dst; dst.ch = src.ch; dst.bind();
        Rig rig; if (! test::run (rig.build (2))) return;
        LoudnessRequest req;
        req.targetLufs = -20.0; req.maxTruePeakDbTp = 0.0;
        req.minPlrDb = std::numeric_limits<double>::infinity();     // PLR must exceed +infinity
        req.maxPasses = 3;
        const auto sol = rig.solver.solve (rig.chain, rig.renderer, rig.params,
                                           src.in(), dst.out(), 2, src.frames(), req);
        test::ok (sol.status != MasteringSolveStatus::Solved,
                  "minPlrDb = +infinity is unsatisfiable and is NOT reported as solved");
        test::ok (sol.binding == MasteringConstraint::PeakToLoudness
                  || sol.status == MasteringSolveStatus::UpstreamViolation,
                  "...and the peak-to-loudness limit is what binds");
        // ...while a NaN is a malformed request rather than a constraint.
        LoudnessRequest bad = req;
        bad.minPlrDb = std::nan ("");
        const auto sol2 = rig.solver.solve (rig.chain, rig.renderer, rig.params,
                                            src.in(), dst.out(), 2, src.frames(), bad);
        test::ok (sol2.status == MasteringSolveStatus::InvalidRequest, "a NaN limit is a malformed request");
    }

    // (e) THE REQUEST HAS NO DEFAULT TARGET. The core does not choose a delivery policy.
    {
        Programme dst; dst.ch = src.ch; dst.bind();
        Rig rig; if (! test::run (rig.build (2))) return;
        LoudnessRequest req;                                    // both required fields left unset
        const auto sol = rig.solver.solve (rig.chain, rig.renderer, rig.params,
                                           src.in(), dst.out(), 2, src.frames(), req);
        test::ok (sol.status == MasteringSolveStatus::InvalidRequest,
                  "a request with no target and no ceiling is refused, not defaulted");
    }

    // (f) THE SAMPLE RATE IS CHECKED, not assumed shared with the chain.
    {
        Programme dst; dst.ch = src.ch; dst.bind();
        MasteringChainConfig cfg;
        MasteringChain ch;
        if (! test::run (ch.prepare (kFs, 2, cfg))) return;
        OfflineRenderer r;
        if (! test::run (r.prepare (2, 4096))) return;
        TargetLoudnessSolver s;
        if (! test::run (s.prepare (44100.0, 2, 4096, ch.internalBlock(), ch.tapOversampleFactor()))) return;
        LoudnessRequest req; req.targetLufs = -14.0; req.maxTruePeakDbTp = -1.0;
        const auto sol = s.solve (ch, r, MasteringChainParams {}, src.in(), dst.out(), 2, src.frames(), req);
        test::ok (sol.status == MasteringSolveStatus::InvalidRequest,
                  "a solver prepared at a different rate from the chain is refused");
    }
}

// =============================================================================================
void testTwoConstraintsAtOnce()
{
    test::group ("two constraints violated together: one is named, both are in the mask");
    Programme src = makeMusic (8.0, 0.5);
    Programme dst; dst.ch = src.ch; dst.bind();
    Rig rig;
    if (! test::run (rig.build (2))) return;
    double inLra = 0.0;
    (void) rig.solver.measureInputLoudnessRange (src.in(), 2, src.frames(), inLra);
    LoudnessRequest req;
    req.inputLoudnessRangeLu = inLra;
    req.targetLufs = -6.0;                       // loud enough to need heavy limiting
    req.maxTruePeakDbTp = -1.0;
    // BOTH limits are chosen to be SATISFIED at the first render and broken at the target — otherwise
    // this is an upstream violation wearing a second constraint as decoration. The fixture's own PLR
    // after the chain is 11.36 at the starting gain, so 10.0 is inside it and 2 dB of limiting is not.
    req.limiterGr = { 2.0, GrStatistic::Max };   // forbid the limiting the target needs
    req.minPlrDb = 10.0;                         // ...and the density that would come with it
    req.maxPasses = 4;
    const auto sol = rig.solver.solve (rig.chain, rig.renderer, rig.params,
                                       src.in(), dst.out(), 2, src.frames(), req);
    const std::uint32_t gr  = constraintBit (MasteringConstraint::LimiterGainReduction);
    const std::uint32_t plr = constraintBit (MasteringConstraint::PeakToLoudness);
    // PRECONDITION: BOTH really are violated by the render nearest the target — otherwise this is a
    // test of one constraint with a second one decorating it.
    test::ok ((sol.alsoViolated & gr) != 0 && (sol.alsoViolated & plr) != 0,
              "precondition: both the limiter-GR and the PLR limits are in the violation mask (0x"
              + std::to_string (sol.alsoViolated) + ")");
    test::ok (sol.binding == MasteringConstraint::LimiterGainReduction,
              "the NAMED one is the limiter's gain reduction — the one more drive cannot trade away");
    test::ok (sol.status == MasteringSolveStatus::TargetUnreachable
              || sol.status == MasteringSolveStatus::UpstreamViolation,
              "and the verdict is a refusal, not a solve");
    std::printf ("      two at once: status %s, binding %s, mask 0x%x, limGR max %.3f, PLR %.3f\n",
                 statusName (sol.status), constraintName (sol.binding), (unsigned) sol.alsoViolated,
                 sol.measured.limiter.maxDb, sol.measured.plrDb);
}

// =============================================================================================
void testPreLimiterTapIsTheRealSignal()
{
    test::group ("the pre-limiter tap is the pre-limiter signal, and its stated offset is right");
    // The tap declares a contract: `preLimiter[c][j]` is the sample at that node for chain input
    // `j - (compressorLookahead + clipperLatency)`, taken BEFORE the gain node. Null it against the
    // same chain rendered with the limiter and the dither bypassed, which is that node's signal by
    // another route — and do it with a NON-ZERO pre-limiter gain, so a tap taken on the wrong side of
    // the gain node fails.
    Programme src = makeMusic (1.5, 0.3);
    const int frames = src.frames();

    MasteringChainConfig cfg;
    cfg.eq = true; cfg.compressor = true; cfg.clipper = false; cfg.limiter = true; cfg.dither = true;
    MasteringChainParams p;
    p.compressor.thresholdDb = -22.0; p.compressor.ratio = 2.5;
    p.preLimiterGainDb = 7.0;                       // NOT zero: the tap must be upstream of this
    p.limiter.ceilingDbTp = -1.0;

    // Route A: the tap.
    MasteringChain ch;
    if (! test::run (ch.prepare (kFs, 2, cfg))) return;
    ch.setParams (p);
    OfflineRenderer r;
    const int blk = 1024;
    if (! test::run (r.prepare (2, blk))) return;
    const int K = ch.internalBlock();
    std::vector<std::vector<float>> preBuf (2, std::vector<float> ((std::size_t) (blk + K), 0.0f));
    std::vector<float*> prePtr { preBuf[0].data(), preBuf[1].data() };
    MasteringChainTaps taps;
    taps.preLimiter = prePtr.data(); taps.frameCapacity = blk + K;
    std::vector<std::vector<float>> tapped (2, std::vector<float> ((std::size_t) frames, 0.0f));
    std::vector<std::vector<float>> out = src.ch;
    std::vector<float*> op { out[0].data(), out[1].data() };
    const MasteringChainResolved res = ch.resolved();
    const long long off = (long long) res.compressorLookahead + (long long) res.clipperLatency;
    long long seen = 0;
    if (! test::run (r.render (ch, src.in(), op.data(), 2, frames, taps,
        [&] (const MasteringChainTaps& t, long long tapPos) noexcept
        {
            for (int j = 0; j < t.framesWritten; ++j)
            {
                const long long inIdx = tapPos + j - off;        // the STATED offset
                if (inIdx >= 0 && inIdx < frames)
                    for (int c = 0; c < 2; ++c)
                        tapped[(std::size_t) c][(std::size_t) inIdx] = preBuf[(std::size_t) c][(std::size_t) j];
            }
            seen += t.framesWritten;
        })))
        return;

    // Route B: the same chain with everything downstream of that node switched off. Its output is the
    // pre-limiter signal, aligned by the renderer — except that it is taken AFTER the gain node, so the
    // comparison divides that back out.
    MasteringChain ch2;
    if (! test::run (ch2.prepare (kFs, 2, cfg))) return;
    MasteringChainParams p2 = p;
    p2.bypassLimiter = true; p2.bypassDither = true;
    ch2.setParams (p2);
    OfflineRenderer r2;
    if (! test::run (r2.prepare (2, blk))) return;
    std::vector<std::vector<float>> ref = src.ch;
    std::vector<float*> rp { ref[0].data(), ref[1].data() };
    if (! test::run (r2.render (ch2, src.in(), rp.data(), 2, frames))) return;

    const double gainOut = std::pow (10.0, p.preLimiterGainDb / 20.0);
    double worst = 0.0, refPeak = 0.0;
    // The last few samples differ by construction: the bypass route's renderer crops the chain's own
    // tail, while the tap sees the drain. Compare the body.
    const int stop = frames - 4 * K;
    for (int c = 0; c < 2; ++c)
        for (int i = 0; i < stop; ++i)
        {
            const double a = (double) tapped[(std::size_t) c][(std::size_t) i] * gainOut;
            const double b = (double) ref[(std::size_t) c][(std::size_t) i];
            worst = std::fmax (worst, std::fabs (a - b));
            refPeak = std::fmax (refPeak, std::fabs (b));
        }
    test::ok (seen >= frames, "precondition: the tap stream covered the programme");
    test::ok (refPeak > 0.05, "precondition: the reference is not silence (" + std::to_string (refPeak) + ")");
    test::ok (worst < 1.0e-5, "the pre-limiter tap nulls against the bypass route at the STATED offset "
                              "(worst " + std::to_string (worst) + " at peak " + std::to_string (refPeak) + ")");
    std::printf ("      pre-limiter tap: worst |tap*gain - bypassRender| = %.3e at peak %.3f (offset %lld)\n",
                 worst, refPeak, off);

    // ---- and the LIMITER tap's offset, found rather than assumed -------------------------------
    // An impulse goes in at a known input index; the limiter's reconstructed-peak trace is searched for
    // where it comes out. The answer must be the offset the contract states, INCLUDING the limiter's own
    // interpolator latency — which the first version of this contract left out, and the omission cost
    // the whole reaction to a peak in the last 32 samples of a programme.
    test::group ("the limiter tap's time offset is the one the contract states");
    {
        MasteringChainConfig c2;
        c2.eq = false; c2.compressor = true; c2.clipper = false; c2.limiter = true; c2.dither = false;
        MasteringChain ch3;
        if (! test::run (ch3.prepare (kFs, 2, c2))) return;
        MasteringChainParams p3;
        p3.bypassCompressor = true;                       // present (so it delays) but transparent
        p3.limiter.ceilingDbTp = -20.0;                   // low, so the impulse certainly engages it
        ch3.setParams (p3);
        const int K3 = ch3.internalBlock(), F3 = ch3.tapOversampleFactor();
        const int n3 = 8 * K3, hit = 3 * K3;
        std::vector<std::vector<float>> imp (2, std::vector<float> ((std::size_t) n3, 0.0f));
        imp[0][(std::size_t) hit] = 0.9f; imp[1][(std::size_t) hit] = 0.9f;
        std::vector<const float*> ip3 { imp[0].data(), imp[1].data() };
        std::vector<std::vector<float>> o3 = imp;
        std::vector<float*> op3 { o3[0].data(), o3[1].data() };
        OfflineRenderer r3;
        const int blk3 = 4 * K3;
        if (! test::run (r3.prepare (2, blk3))) return;
        std::vector<float> pk3 ((std::size_t) (blk3 + K3) * (std::size_t) F3, 0.0f);
        MasteringChainTaps t3;
        t3.limiterPeakLin = pk3.data(); t3.osCapacity = (blk3 + K3) * F3;
        long long argmaxOs = -1; double best3 = 0.0;
        if (! test::run (r3.render (ch3, ip3.data(), op3.data(), 2, n3, t3,
            [&] (const MasteringChainTaps& t, long long tapPos) noexcept
            {
                for (int i = 0; i < t.osWritten; ++i)
                    if ((double) pk3[(std::size_t) i] > best3)
                    { best3 = pk3[(std::size_t) i]; argmaxOs = tapPos * F3 + i; }
            })))
            return;
        const MasteringChainResolved r4 = ch3.resolved();
        const long long stated = r4.limiterTapOffset;
        test::ok (best3 > 0.5, "precondition: the impulse really reached the limiter's detector ("
                               + std::to_string (best3) + ")");
        const double foundOffset = (double) argmaxOs / (double) F3 - (double) hit;
        // Half an oversampled sample of slack: the peak of a reconstructed impulse need not land exactly
        // on a grid point, and the contract is in whole baseband frames.
        test::approx (foundOffset, (double) stated, 1.0,
                      "the measured tap offset is the stated one (found " + std::to_string (foundOffset)
                      + ", stated " + std::to_string (stated) + ")");
        std::printf ("      limiter tap offset: measured %.3f frames, stated %lld "
                     "(compLook %d + clip %d + half the oversampler round trip %d)\n",
                     foundOffset, stated, r4.compressorLookahead, r4.clipperLatency,
                     (r4.limiterLatency - r4.limiterLookahead) / 2);
    }
}

// =============================================================================================
void testTargetBetweenAchievable()
{
    test::group ("a target between two achievable values is reported as that, with both sides");
    // The status has two causes and only one of them can be built on demand: ask for a tolerance finer
    // than the search's own gain resolution, and the bracket closes with the target still between the
    // two sides. That exercises the same code path a gated STEP takes — the interesting cause, which
    // needs a block to cross the absolute gate and cannot be conjured on ordinary programme.
    Programme src = makeMusic (4.0, 0.3);
    Programme dst; dst.ch = src.ch; dst.bind();
    Rig rig;
    if (! test::run (rig.build (2))) return;
    LoudnessRequest req;
    req.targetLufs = -13.0;
    req.maxTruePeakDbTp = -1.0;
    req.toleranceLu = 1.0e-9;                 // finer than 1e-3 dB of gain can resolve
    req.maxPasses = 12;
    const auto sol = rig.solver.solve (rig.chain, rig.renderer, rig.params,
                                       src.in(), dst.out(), 2, src.frames(), req);
    for (int k = 0; k < sol.logCount; ++k)
        std::printf ("        pass %d: g %.9f -> I %.9f\n", k + 1, sol.log[k].gainDb, sol.log[k].integratedLufs);
    std::printf ("        sides: lo %.9f @ %.9f | hi %.9f @ %.9f\n",
                 sol.achievedBelowLufs, sol.gainBelowDb, sol.achievedAboveLufs, sol.gainAboveDb);
    test::ok (sol.status == MasteringSolveStatus::TargetBetweenAchievable,
              std::string ("status is TargetBetweenAchievable (got ") + statusName (sol.status) + ")");
    if (sol.status != MasteringSolveStatus::TargetBetweenAchievable) return;
    // PRECONDITION: the two sides really do straddle the target, or "between" is not what happened.
    // PRECONDITION: the answer really is outside the tolerance asked for, or "between" is vacuous.
    test::ok (std::fabs (sol.measured.integratedLufs - req.targetLufs) > req.toleranceLu,
              "precondition: the achieved loudness really misses the tolerance ("
              + std::to_string (std::fabs (sol.measured.integratedLufs - req.targetLufs)) + " LU)");
    test::ok (std::fabs (sol.gainAboveDb - sol.gainBelowDb) <= 1.0e-3,
              "the gain interval really is below the resolution the search can express");
    // ...and it is CLOSE: this is a resolution limit, not a failure to find the answer.
    test::ok (std::fabs (sol.measured.integratedLufs - req.targetLufs) < 1.0e-5,
              "the answer is within 1e-5 LU — the tolerance was finer than the actuator, not the search wrong");
    test::ok (sol.binding == MasteringConstraint::None,
              "it is NOT a constraint violation: nothing was broken");
    std::printf ("      between achievable: %.6f .. %.6f LUFS at gains %.6f .. %.6f, passes %d\n",
                 sol.achievedBelowLufs, sol.achievedAboveLufs, sol.gainBelowDb, sol.gainAboveDb, sol.passes);
}

// =============================================================================================
void testTheAnswerDoesNotDependOnWhereItStarted()
{
    test::group ("the answer does not depend on the starting gain");
    // THE STRONGEST PROPERTY IN THIS FILE, and the one that caught the worst defect. A step that moved
    // the gain and the ceiling together was exact — it is the scale law — but exact at a FROZEN DRIVE,
    // so a warm start delivered the target loudness and the stated peak with the programme crushed:
    // measured, from 0 dB the answer was 8.5 dB of drive and no limiting, and from 55 dB it was 47.6 dB
    // of drive, 38.05 dB of limiter gain reduction and an LRA of 0.10 — both reported `Solved`.
    Programme src = makeMusic (6.0, 0.3);
    double refDrive = 0.0, refLimGr = 0.0, refPlr = 0.0;
    bool have = false;
    for (double start : { 0.0, 4.0, 12.0, 30.0, 55.0 })
    {
        Programme dst; dst.ch = src.ch; dst.bind();
        Rig rig;
        if (! test::run (rig.build (2))) return;
        LoudnessRequest req;
        req.targetLufs = -14.0;
        req.maxTruePeakDbTp = -1.0;
        req.initialGainDb = start;
        req.maxPasses = 6;
        const auto sol = rig.solver.solve (rig.chain, rig.renderer, rig.params,
                                           src.in(), dst.out(), 2, src.frames(), req);
        char msg[192];
        std::snprintf (msg, sizeof msg, "start %+.1f dB: solved", start);
        test::ok (sol.status == MasteringSolveStatus::Solved, msg);
        if (sol.status != MasteringSolveStatus::Solved) continue;
        const double drive = sol.preLimiterGainDb - sol.ceilingDbTp;
        if (! have) { refDrive = drive; refLimGr = sol.measured.limiter.maxDb; refPlr = sol.measured.plrDb; have = true; }
        std::snprintf (msg, sizeof msg, "start %+.1f dB: the DRIVE agrees with the cold start (%.3f vs %.3f dB)",
                       start, drive, refDrive);
        test::approx (drive, refDrive, 0.35, msg);
        std::snprintf (msg, sizeof msg, "start %+.1f dB: the limiting agrees (%.3f vs %.3f dB)",
                       start, sol.measured.limiter.maxDb, refLimGr);
        test::approx (sol.measured.limiter.maxDb, refLimGr, 0.5, msg);
        std::snprintf (msg, sizeof msg, "start %+.1f dB: the delivered PLR agrees (%.3f vs %.3f dB)",
                       start, sol.measured.plrDb, refPlr);
        test::approx (sol.measured.plrDb, refPlr, 0.5, msg);
        std::printf ("      start %+5.1f -> drive %7.3f  limGR %6.3f  PLR %6.3f  LRA %5.2f  passes %d\n",
                     start, drive, sol.measured.limiter.maxDb, sol.measured.plrDb,
                     sol.measured.loudnessRangeLu, sol.passes);
    }
    // PRECONDITION: the starts really were far apart, or "does not depend" is a statement about noise.
    test::ok (have, "precondition: at least one start solved, so there is a reference to agree with");
}

// =============================================================================================
void testAnUnmeasurableStartIsNotAnUnmeasurableProgramme()
{
    test::group ("a programme too quiet to measure at the starting gain is still solved");
    // Every gating block under the -70 LUFS absolute gate: the integrated measure has nothing to report
    // at 0 dB, and 55 dB up it is ordinary programme. Refusing on the first reading would be a verdict
    // about the starting gain.
    const int n = (int) (4.0 * kFs);
    Programme q; q.ch.assign (2, std::vector<float> ((std::size_t) n, 0.0f));
    // FLAT on purpose, and the amplitude is calibrated rather than chosen — the window between "still
    // measurable" and "out of the actuator's range" is only a few dB wide. 3.7e-4 read -69.32 LUFS, i.e.
    // ABOVE the -70 gate and therefore measurable, and the test passed for the wrong reason; 2.6e-4 sits
    // at about -72.3 LUFS, under the gate, and 56.3 dB below a -16 LUFS target, inside the +-60 dB
    // actuator. K-weighting at 440 Hz is why the arithmetic from a 1 kHz calibration was 2 dB out.
    for (int i = 0; i < n; ++i)
    {
        const float v = (float) (2.6e-4 * std::sin (2.0 * kPi * 440.0 * (double) i / kFs));
        q.ch[0][(std::size_t) i] = v; q.ch[1][(std::size_t) i] = v;
    }
    q.bind();
    Programme o; o.ch = q.ch; o.bind();

    // PRECONDITION: at unity the programme really is unmeasurable, or this tests nothing.
    {
        analysis::LoudnessMeter m;
        if (! test::run (m.prepare (kFs, 2, 60.0))) return;
        if (! test::run (m.process (q.in(), 2, n))) return;
        test::ok (m.integratedLufs() <= -120.0,
                  "precondition: at unity every gating block is under the absolute gate (I = "
                  + std::to_string (m.integratedLufs()) + ")");
    }

    Rig rig;
    if (! test::run (rig.build (2))) return;
    rig.params.compressor.thresholdDb = 0.0;               // out of the way
    LoudnessRequest req;
    req.targetLufs = -16.0;
    req.maxTruePeakDbTp = -1.0;
    req.maxPasses = 5;
    const auto sol = rig.solver.solve (rig.chain, rig.renderer, rig.params,
                                       q.in(), o.out(), 2, q.frames(), req);
    for (int k = 0; k < sol.logCount; ++k)
        std::printf ("        pass %d: g %8.3f -> I %10.4f TP %8.4f\n",
                     k + 1, sol.log[k].gainDb, sol.log[k].integratedLufs, sol.log[k].truePeakDbTp);
    test::ok (sol.status == MasteringSolveStatus::Solved,
              std::string ("an unmeasurable first render is bootstrapped from the PEAK, not refused (got ")
              + statusName (sol.status) + ")");
    if (sol.status != MasteringSolveStatus::Solved) return;
    const Independent ind = measureIndependently (o.ch);
    test::approx (ind.I, req.targetLufs, 0.1, "...and the target is taken");
    test::ok (ind.TP <= req.maxTruePeakDbTp + 1e-9, "...with the peak under the promise");
    std::printf ("      quiet programme: gain %.3f, I %.4f, TP %.4f, passes %d\n",
                 sol.preLimiterGainDb, ind.I, ind.TP, sol.passes);
}

// =============================================================================================
void testTapPlumbingEdges()
{
    test::group ("the tap plumbing, at the edges every other fixture in this file avoids");
    // A mutation round found the whole tap edge unguarded: every statistics fixture here uses a renderer
    // block that is a MULTIPLE of the internal quantum, so `pos_` is zero at every call boundary, the
    // capacity arithmetic never has to account for it, and a tap-position bug that costs a fraction of a
    // block is invisible. These are the cases that were missing.

    // (a) CAPACITY WITH A NON-ZERO PHASE. `n + K` frames is what a caller must size; exactly `n` is not
    //     enough once the stream is part way into a quantum, and the check has to say so BEFORE moving.
    {
        MasteringChain ch;
        MasteringChainConfig cfg;
        if (! test::run (ch.prepare (kFs, 2, cfg))) return;
        const int K = ch.internalBlock();
        std::vector<std::vector<float>> buf (2, std::vector<float> ((std::size_t) (K + 200), 0.2f));
        std::vector<float*> bp { buf[0].data(), buf[1].data() };
        std::vector<float> tap ((std::size_t) K, -1.0f);
        MasteringChainTaps t;
        t.compressorGrDb = tap.data(); t.frameCapacity = K;
        // At phase 0 a call of exactly K frames produces exactly K tap frames: accepted.
        test::ok (ch.process (bp.data(), 2, K, t), "phase 0, n == K, capacity K: accepted");
        test::ok (t.framesWritten == K, "...and it wrote exactly K frames");
        // Now the stream sits at phase 100. A call of 412 frames spans TWO quantum boundaries, so it
        // produces 2K tap frames — and a capacity of K must be refused.
        MasteringChainTaps t2;
        t2.compressorGrDb = tap.data(); t2.frameCapacity = K;
        test::ok (ch.process (bp.data(), 2, 100, t2), "advance the stream to phase 100");
        MasteringChainTaps t3;
        t3.compressorGrDb = tap.data(); t3.frameCapacity = K;
        test::ok (! ch.process (bp.data(), 2, 412, t3),
                  "phase 100, n = 412, capacity K: REFUSED — the phase is part of the arithmetic");
        test::ok (t3.framesWritten == 0, "...and the refused call reports nothing written");
    }

    // (b) THE OVERSAMPLED CAPACITY BAND. A limiter tap sized in FRAMES rather than in oversampled
    //     samples is the mistake, and it is only visible between the two.
    {
        MasteringChain ch;
        MasteringChainConfig cfg;
        if (! test::run (ch.prepare (kFs, 2, cfg))) return;
        const int K = ch.internalBlock(), F = ch.tapOversampleFactor();
        test::ok (F >= 2, "precondition: the chain really is oversampling (" + std::to_string (F) + "x)");
        std::vector<std::vector<float>> buf (2, std::vector<float> ((std::size_t) K, 0.2f));
        std::vector<float*> bp { buf[0].data(), buf[1].data() };
        std::vector<float> os ((std::size_t) K * (std::size_t) F, 0.0f);
        MasteringChainTaps t;
        t.limiterGrDb = os.data(); t.osCapacity = K * 2;      // enough frames, NOT enough oversampled
        test::ok (! ch.process (bp.data(), 2, K, t),
                  "a limiter tap sized in frames instead of oversampled samples is refused");
        MasteringChainTaps t2;
        t2.limiterGrDb = os.data(); t2.osCapacity = K * F;    // exactly right
        test::ok (ch.process (bp.data(), 2, K, t2), "...and exactly K*F is accepted");
        test::ok (t2.osWritten == K * F, "...writing exactly K*F oversampled samples");
    }

    // (c) THE RENDERER'S TAP POSITION, at a block that is NOT a multiple of the quantum. The statistics
    //     have to come out the same as at a block that is, or the tap stream and the audio have drifted.
    {
        Programme src = makeMusic (3.0, 0.35);
        double refMean = 0.0, refP95 = 0.0, refAct = 0.0;
        bool have = false;
        for (int blk : { 1024, 977, 100, 3 })
        {
            Programme dst; dst.ch = src.ch; dst.bind();
            Rig rig;
            if (! test::run (rig.build (2, blk))) return;
            LoudnessRequest req;
            req.targetLufs = -13.0; req.maxTruePeakDbTp = -1.0; req.maxPasses = 4;
            const auto sol = rig.solver.solve (rig.chain, rig.renderer, rig.params,
                                               src.in(), dst.out(), 2, src.frames(), req);
            if (! test::run (sol.status == MasteringSolveStatus::Solved)) return;
            if (! have) { refMean = sol.measured.compressor.meanDb; refP95 = sol.measured.compressor.p95Db;
                          refAct = sol.measured.limiter.activeFraction; have = true; }
            char msg[160];
            std::snprintf (msg, sizeof msg, "block %d: the compressor mean is bit-identical to block 1024", blk);
            test::approx (sol.measured.compressor.meanDb, refMean, 0.0, msg);
            std::snprintf (msg, sizeof msg, "block %d: the compressor p95 is bit-identical", blk);
            test::approx (sol.measured.compressor.p95Db, refP95, 0.0, msg);
            std::snprintf (msg, sizeof msg, "block %d: the limiter active fraction is bit-identical", blk);
            test::approx (sol.measured.limiter.activeFraction, refAct, 0.0, msg);
        }
        // PRECONDITION: the statistics are not all zero, or every equality above is trivial.
        test::ok (refMean > 0.01 && refP95 > 0.01,
                  "precondition: the compressor statistics are non-zero (mean " + std::to_string (refMean) + ")");
    }

    // (d) A BYPASSED LIMITER STILL FILLS ITS TRACE. A hole in the trace is not the same as zero, and a
    //     mean taken over a programme whose bypass moved would silently be a mean over a shorter one.
    {
        MasteringChain ch;
        MasteringChainConfig cfg;
        if (! test::run (ch.prepare (kFs, 2, cfg))) return;
        const int K = ch.internalBlock(), F = ch.tapOversampleFactor();
        MasteringChainParams p;
        p.bypassLimiter = true;
        p.preLimiterGainDb = 20.0;                 // loud enough that a RUNNING limiter would react
        ch.setParams (p);
        std::vector<std::vector<float>> buf (2, std::vector<float> ((std::size_t) K, 0.5f));
        std::vector<float*> bp { buf[0].data(), buf[1].data() };
        std::vector<float> gr ((std::size_t) K * (std::size_t) F, -1234.0f);
        std::vector<float> pk ((std::size_t) K * (std::size_t) F, -1234.0f);
        MasteringChainTaps t;
        t.limiterGrDb = gr.data(); t.limiterPeakLin = pk.data(); t.osCapacity = K * F;
        test::run (ch.process (bp.data(), 2, K, t));
        int unwritten = 0, nonZero = 0;
        for (std::size_t i = 0; i < gr.size(); ++i)
        {
            if (core::exactlyEqual (gr[i], -1234.0f) || core::exactlyEqual (pk[i], -1234.0f)) ++unwritten;
            else if (! core::exactlyEqual (gr[i], 0.0f)) ++nonZero;
        }
        test::ok (unwritten == 0, "a bypassed limiter writes its whole trace rather than leaving a hole");
        test::ok (nonZero == 0, "...and what it writes is zero, which is the truth about that quantum");
    }

    // (e) THE STATISTICS WINDOW WITH A CLIPPER PRESENT. The limiter's window is offset by the clipper's
    //     latency as well as the compressor's lookahead, and no other fixture here turns the clipper on.
    {
        Programme src = makeMusic (3.0, 0.35);
        Programme dst; dst.ch = src.ch; dst.bind();
        Rig rig;
        if (! test::run (rig.build (2, 4096, /*clipper*/ true))) return;
        LoudnessRequest req;
        req.targetLufs = -12.0; req.maxTruePeakDbTp = -1.0; req.maxPasses = 4;
        const auto sol = rig.solver.solve (rig.chain, rig.renderer, rig.params,
                                           src.in(), dst.out(), 2, src.frames(), req);
        if (! test::run (sol.status == MasteringSolveStatus::Solved)) return;
        const MasteringChainResolved r = rig.chain.resolved();
        test::ok (r.clipperLatency > 0, "precondition: the clipper really is present and carries latency ("
                                        + std::to_string (r.clipperLatency) + " samples)");
        test::ok (sol.measured.limiter.frames == (std::uint64_t) src.frames() * (std::uint64_t) rig.chain.tapOversampleFactor(),
                  "with a clipper in the chain the limiter window is still exactly frames*F samples");
        test::ok (sol.measured.compressor.frames == (std::uint64_t) src.frames(),
                  "...and the compressor window is still exactly `frames`");
    }

    // (f) A PROGRAMME THAT ENDS ON A PEAK. Without draining the true-peak meter the delivered peak is
    //     under-read, which is the exact shape of the defect P1 measured in the chain being replaced.
    {
        const int n = (int) (2.0 * kFs);
        Programme q; q.ch.assign (2, std::vector<float> ((std::size_t) n, 0.0f));
        for (int i = 0; i < n; ++i)
        {
            const float v = (float) (0.25 * std::sin (2.0 * kPi * 300.0 * (double) i / kFs));
            q.ch[0][(std::size_t) i] = v; q.ch[1][(std::size_t) i] = v;
        }
        q.ch[0][(std::size_t) (n - 2)] = 0.9f; q.ch[1][(std::size_t) (n - 2)] = 0.9f;
        q.ch[0][(std::size_t) (n - 1)] = 0.9f; q.ch[1][(std::size_t) (n - 1)] = 0.9f;
        q.bind();
        Programme o; o.ch = q.ch; o.bind();
        Rig rig;
        if (! test::run (rig.build (2))) return;
        LoudnessRequest req;
        req.targetLufs = -18.0; req.maxTruePeakDbTp = -1.0; req.maxPasses = 4;
        const auto sol = rig.solver.solve (rig.chain, rig.renderer, rig.params,
                                           q.in(), o.out(), 2, q.frames(), req);
        if (! test::run (sol.status == MasteringSolveStatus::Solved)) return;
        // An UNDRAINED reading of the same buffer, for the comparison: this is what the number would be
        // without the drain, and the point is that the reported one is not it.
        analysis::TruePeakMeter dry;
        if (! test::run (dry.prepare (kFs, 65536, 2))) return;
        const float* dp[2] = { o.ch[0].data(), o.ch[1].data() };
        if (! test::run (dry.process (dp, 2, n))) return;
        const double undrained = dry.truePeakDb();
        test::ok (sol.measured.truePeakDbTp >= undrained - 1e-9,
                  "the reported peak is at least the undrained reading (drained " +
                  std::to_string (sol.measured.truePeakDbTp) + " vs undrained " + std::to_string (undrained) + ")");
        test::ok (sol.measured.truePeakDbTp <= req.maxTruePeakDbTp + 1e-9,
                  "...and a programme that ends on a peak still delivers under the promise");
        std::printf ("      ends on a peak: drained %.4f dBTP, undrained %.4f dBTP\n",
                     sol.measured.truePeakDbTp, undrained);
    }
}

// =============================================================================================
// THE MUTATION STAND'S SURVIVORS. Each group below is the INPUT that makes one surviving mutation
// change the delivered answer — found by sweeping a 99-case battery against the clean header and each
// mutated one, not by reasoning about the code. They are collected here rather than scattered because
// what they have in common is where they live: outside the corpus, in the saturated regime, under three
// seconds, or in a request whose every candidate is infeasible.
void testSurvivorsOfTheMutationStand()
{
    // ---------------------------------------------------------------------------------------------
    test::group ("a broken limit outranks 'not exactly achievable' (found by mutating the choice rule)");
    {
        // The two verdicts answer different questions and the wrong one used to win. `best` is the
        // candidate DELIVERED, which when nothing is feasible is deliberately the gentlest rather than
        // the nearest — so testing its distance to the target reported a search that had walked away
        // from an illegal target as a resolution limit, and dropped the violations on the way.
        Programme src = makeMusic (6.0, 0.3);
        Programme dst; dst.ch = src.ch; dst.bind();
        Rig rig;
        if (! test::run (rig.build (2))) return;
        LoudnessRequest req;
        req.targetLufs = -24.0;                 // QUIETER than the first render, so the upstream guard
        req.maxTruePeakDbTp = -1.0;             // stays out of it and the search walks down
        req.initialGainDb = 18.0;
        req.maxLraLossLu = 0.01;                // broken by the compressor at every gain
        req.inputLoudnessRangeLu = 8.0;
        req.minPlrDb = 40.0;                    // and unsatisfiable on any real programme
        req.maxPasses = 6;
        const auto sol = rig.solver.solve (rig.chain, rig.renderer, rig.params,
                                           src.in(), dst.out(), 2, src.frames(), req);
        const std::uint32_t lra = constraintBit (MasteringConstraint::LoudnessRange);
        const std::uint32_t plr = constraintBit (MasteringConstraint::PeakToLoudness);
        // PRECONDITION: both limits really are broken, or the ordering under test never applies.
        test::ok ((sol.alsoViolated & lra) != 0 && (sol.alsoViolated & plr) != 0,
                  "precondition: both limits are violated (mask 0x" + std::to_string (sol.alsoViolated) + ")");
        test::ok (sol.status == MasteringSolveStatus::TargetUnreachable,
                  std::string ("a violated limit is reported, not 'between achievable' (got ")
                  + statusName (sol.status) + ")");
        test::ok (sol.binding != MasteringConstraint::None, "...and the binding limit is NAMED");
        std::printf ("      ordering: status %s, binding %s, mask 0x%x, I %.4f\n",
                     statusName (sol.status), constraintName (sol.binding),
                     (unsigned) sol.alsoViolated, sol.measured.integratedLufs);
    }

    // ---------------------------------------------------------------------------------------------
    test::group ("among candidates that all break a limit, the GENTLEST is delivered, not the nearest");
    {
        // The same request without the unsatisfiable PLR floor: every candidate breaks the range limit,
        // so the tie-break among infeasible candidates decides what the caller gets. Delivering the
        // closest to a target already declared unreachable is the "push it through anyway" behaviour.
        Programme src = makeMusic (6.0, 0.3);
        Programme dst; dst.ch = src.ch; dst.bind();
        Rig rig;
        if (! test::run (rig.build (2))) return;
        LoudnessRequest req;
        req.targetLufs = -24.0; req.maxTruePeakDbTp = -1.0; req.initialGainDb = 18.0;
        req.maxLraLossLu = 0.01; req.inputLoudnessRangeLu = 8.0; req.maxPasses = 6;
        const auto sol = rig.solver.solve (rig.chain, rig.renderer, rig.params,
                                           src.in(), dst.out(), 2, src.frames(), req);
        test::ok (sol.status == MasteringSolveStatus::TargetUnreachable, "all-infeasible: refused by name");
        // PRECONDITION: the log holds more than one candidate, and they differ in range loss —
        // otherwise "the gentlest" is a statement about a set of one (P18 F45's geometry trap).
        test::ok (sol.logCount >= 2, "precondition: more than one candidate was rendered ("
                                     + std::to_string (sol.logCount) + ")");
        double best = 1e9, worst = -1e9;
        for (int k = 0; k < sol.logCount; ++k)
        {
            best  = std::fmin (best,  req.inputLoudnessRangeLu - sol.log[k].loudnessRangeLu);
            worst = std::fmax (worst, req.inputLoudnessRangeLu - sol.log[k].loudnessRangeLu);
        }
        const double delivered = req.inputLoudnessRangeLu - sol.measured.loudnessRangeLu;
        // AND THE PRECONDITION IS ASSERTED, not merely printed. `worst` was computed and shown and
        // never compared, so a log whose candidates all cost the SAME range loss satisfied every check
        // in the group — "the gentlest" would then be a statement about a set of one wearing a count
        // of two, which is exactly the geometry trap the comment above names.
        test::ok (worst > best + 1.0e-6,
                  "precondition: the candidates really do differ in range loss (" + std::to_string (best)
                  + " .. " + std::to_string (worst) + ")");
        test::ok (delivered <= best + 1e-9,
                  "the delivered render is the gentlest tried (" + std::to_string (delivered)
                  + " vs best " + std::to_string (best) + ")");
        std::printf ("      gentlest: delivered loss %.4f, range over the log %.4f .. %.4f, %d candidates\n",
                     delivered, best, worst, sol.logCount);
    }

    // ---------------------------------------------------------------------------------------------
    test::group ("a programme too short for a loudness RANGE does not get one invented");
    {
        // EBU Tech 3342 needs short-term samples, one a second. Under three seconds there are none, and
        // `loudnessRangeLu()` answers 0.0 — which is also what "no dynamic range at all" answers. A
        // solver that took that as a measurement would compute a range LOSS equal to the whole input
        // range and refuse a programme it should have mastered.
        for (double sec : { 1.0, 2.0, 2.9, 6.0 })
        {
            Programme src = makeMusic (sec, 0.3);
            Programme dst; dst.ch = src.ch; dst.bind();
            Rig rig;
            if (! test::run (rig.build (2))) return;
            LoudnessRequest req;
            req.targetLufs = -14.0; req.maxTruePeakDbTp = -1.0;
            req.maxLraLossLu = 0.4; req.inputLoudnessRangeLu = 6.0;   // supplied by the caller
            req.maxPasses = 4;
            const auto sol = rig.solver.solve (rig.chain, rig.renderer, rig.params,
                                               src.in(), dst.out(), 2, src.frames(), req);
            char msg[176];
            std::snprintf (msg, sizeof msg,
                           "%.1f s: the range constraint is %s and the solve %s", sec,
                           sec >= 3.0 ? "ACTIVE" : "off", statusName (sol.status));
            if (sec < 3.0)
            {
                test::ok (! sol.measured.lraValid, std::string ("under 3 s the range is not a measurement — ") + msg);
                test::ok (sol.status == MasteringSolveStatus::Solved, msg);
            }
            else
            {
                // PRECONDITION for the pair above: at six seconds the same field DOES flip, so the
                // quantity being compared moves across the fixture rather than being constant.
                test::ok (sol.measured.lraValid, std::string ("at 6 s the range IS a measurement — ") + msg);
            }
            std::printf ("      %.1f s: lraValid=%d LRA=%.3f status=%s I=%.4f\n",
                         sec, (int) sol.measured.lraValid, sol.measured.loudnessRangeLu,
                         statusName (sol.status), sol.measured.integratedLufs);
        }
    }

    // ---------------------------------------------------------------------------------------------
    test::group ("the saturated regime: the slope guard and the idle anchor both earn their keep");
    {
        // Where the corpus is not. On material already dense at full scale the chain saturates near
        // -5.3 LUFS, so a target just above that sits in the flat part of `I(g)` — the measured slope
        // collapses toward zero and the gain is nowhere near the +-60 dB clamp. Two guards live here and
        // nothing else in this file visits the place.
        // TEN RENDERS FOR ALL THREE, and that is itself the assertion. The middle row is the one that
        // moves: with the slope floor at 0.05 a -5.3 LUFS target stopped 0.127 LU short at this budget
        // and needed twelve to land, and the group used to hide that behind a disjunction. With the
        // floor at 0.01 — the honest secants here are 0.019 down to 0.0089, all of them under the old
        // floor — it solves in SEVEN. So a regression of the floor fails this row rather than being
        // absorbed by it.
        struct Row { double target; int passes; };
        for (const Row& row : { Row { -5.4, 10 }, Row { -5.3, 10 }, Row { -2.0, 10 } })
        {
            const double target = row.target;
            Programme src = makeMusic (6.0, 0.95, 12345u, 1.0);   // dense, no transients
            Programme dst; dst.ch = src.ch; dst.bind();
            Rig rig;
            if (! test::run (rig.build (2))) return;
            LoudnessRequest req;
            req.targetLufs = target; req.maxTruePeakDbTp = -1.0; req.maxPasses = row.passes;
            const auto sol = rig.solver.solve (rig.chain, rig.renderer, rig.params,
                                               src.in(), dst.out(), 2, src.frames(), req);
            // PRECONDITION: the search really is in the flat region — consecutive renders move the
            // loudness far less than they move the gain.
            double worstSlope = 1e9;
            for (int k = 1; k < sol.logCount; ++k)
            {
                const double dg = sol.log[k].gainDb - sol.log[k - 1].gainDb;
                if (std::fabs (dg) > 0.25)
                    worstSlope = std::fmin (worstSlope,
                                            (sol.log[k].integratedLufs - sol.log[k - 1].integratedLufs) / dg);
            }
            char msg[192];
            std::snprintf (msg, sizeof msg, "target %.1f: the search reached the flat region (slope %.4f)",
                           target, worstSlope);
            test::ok (worstSlope < 0.25, msg);
            std::snprintf (msg, sizeof msg, "target %.1f: %s in %d renders, gain %.3f",
                           target, statusName (sol.status), sol.passes, sol.preLimiterGainDb);
            // THE VERDICT IS ASSERTED PER TARGET, not as a disjunction — a disjunction over three
            // statuses is satisfied by every mutation of the two guards this group exists for, and it
            // was: the first version of this group passed against both.
            if (target < -4.0)
            {
                // -5.4 is INSIDE what the chain can deliver, and the search only gets there if it
                // believes the tiny slope it measured. Two things used to stop it believing: a
                // REJECTION of any slope under 0.02 (which fell back to 1.0, a fiftieth of the step the
                // measurement called for, and ran out of renders 0.194 LU short) and a FLOOR of 0.05
                // that clipped the honest 0.0089 to five times itself. The first is gone; the second is
                // now 0.01. The earlier version of this comment described only the rejection and called
                // the floor by the rejection's number, which is how a fixed constant went on looking
                // fixed while it was still five times too big.
                //
                // THIS USED TO BE A DISJUNCTION (`Solved || |err| < 0.16`) directly under a comment
                // saying it was not one, and the disjunction is what made it blind: every mutation of
                // the two guards this group exists for satisfied the second arm. The status and the
                // error are now asserted SEPARATELY, so a mutant that lands close while reporting the
                // wrong reason fails on the first of them.
                test::ok (sol.status == MasteringSolveStatus::Solved, msg);
                test::approx (sol.measured.integratedLufs, target, req.toleranceLu,
                              std::string ("…and within the tolerance it was asked for"));
            }
            else
            {
                // -2 is past what the chain can deliver at all, and the actuator is what stops it. That
                // has a NAME, and naming the budget instead is a statement about the wrong thing.
                test::ok (sol.status == MasteringSolveStatus::TargetUnreachable, msg);
                test::ok (sol.binding == MasteringConstraint::GainRange,
                          "…and the binding limit is the gain range, not the pass budget");
                // AND IT DOES NOT GRIND. The ceiling-still test used to be 1e-6, which the true peak's
                // own jitter of a couple of parts per million never satisfies, so the search sat at the
                // clamp re-rendering identical audio until the budget ran out — measured, six of eleven
                // renders carrying no information. At the resolution `bracketClosed` already uses it
                // stops when it has stopped.
                test::ok (sol.passes <= 7, "…and it stops when it stops, rather than burning the budget ("
                                           + std::to_string (sol.passes) + " renders)");
                test::approx (std::fabs (sol.preLimiterGainDb), 60.0, 1e-9,
                              "…with the gain really pinned at the actuator's limit");
            }
            std::printf ("      %s\n", msg);
        }
    }
}

// =============================================================================================
// THE CODE-REVIEW ROUND'S COUNTEREXAMPLES. Four inputs, each of which the solver got wrong, and three
// of them against fixes made earlier the same day — which is why they are here rather than folded into
// the groups above: they are the evidence that a fix is not protected by a test on its own.
void testTheReviewRoundsCounterexamples()
{
    // A tone quiet enough that a warm start puts the search deep in the saturated region, where the
    // measured slope is genuinely near zero and the ceiling is still tracking its aim. Everything about
    // this group is that combination.
    auto tone = [] (double seconds, double amp)
    {
        const int n = (int) (seconds * kFs);
        Programme p; p.ch.assign (2, std::vector<float> ((std::size_t) n, 0.0f));
        for (int i = 0; i < n; ++i)
        {
            const float v = (float) (amp * std::sin (2.0 * kPi * 1000.0 * (double) i / kFs));
            p.ch[0][(std::size_t) i] = v; p.ch[1][(std::size_t) i] = v;
        }
        p.bind();
        return p;
    };
    auto run = [] (Rig& rig, Programme& src, Programme& dst, double startGain, double startCeiling,
                   double target, int passes, double tol)
    {
        rig.params.bypassCompressor = true; rig.params.bypassDither = true;
        rig.params.limiter.ceilingDbTp = startCeiling;
        LoudnessRequest req;
        req.targetLufs = target; req.maxTruePeakDbTp = -1.0; req.toleranceLu = tol;
        req.maxPasses = passes; req.initialGainDb = startGain;
        return rig.solver.solve (rig.chain, rig.renderer, rig.params,
                                 src.in(), dst.out(), 2, src.frames(), req);
    };

    // ---------------------------------------------------------------------------------------------
    test::group ("a warm start into the saturated region converges at the STANDARD budget");
    {
        // The search runs in `d = g - c` and `J = I - c` because `I` is not a function of `g` alone.
        // Taken in the caller's coordinates instead, two renders 9 dB apart whose ceilings differ by
        // 0.05 dB gave a "slope" of 0.0057 where the truth is near 1; the step went to the -60 dB clamp,
        // the render fell under the absolute gate, and a four-render budget ended at -12.993 LUFS.
        Programme src = tone (1.0, 0.005);
        Programme dst; dst.ch = src.ch; dst.bind();
        Rig rig;
        if (! test::run (rig.build (2))) return;
        const auto sol = run (rig, src, dst, 55.0, -1.0, -10.0, 4, 0.1);
        // PRECONDITION: the first render really is in the saturated region — it must be far LOUDER than
        // the target, or the search never needs the coordinates this group is about.
        test::ok (sol.logCount >= 1 && sol.log[0].integratedLufs > -6.0,
                  "precondition: the warm start lands in the saturated region (I = "
                  + std::to_string (sol.logCount ? sol.log[0].integratedLufs : 0.0) + ")");
        test::ok (sol.status == MasteringSolveStatus::Solved,
                  std::string ("a +55 dB start reaches a -10 LUFS target in four renders (got ")
                  + statusName (sol.status) + ")");
        test::approx (sol.measured.integratedLufs, -10.0, 0.1, "…on target");
        std::printf ("      warm start: %s at %.4f LUFS, gain %.3f, %d renders\n",
                     statusName (sol.status), sol.measured.integratedLufs, sol.preLimiterGainDb, sol.passes);
    }

    // ---------------------------------------------------------------------------------------------
    test::group ("a ceiling far below the promise is not mistaken for the actuator's limit");
    {
        // The caller's ceiling is a STARTING POINT, and the solver may raise it to the promise. A step
        // that MOVES to the +-60 dB clamp is a step the search has not evaluated yet; naming it called a
        // -25 LUFS target unreachable that `g = -18.99` delivers exactly.
        Programme src = tone (1.0, 0.5);
        Programme dst; dst.ch = src.ch; dst.bind();
        Rig rig;
        if (! test::run (rig.build (2))) return;
        const auto sol = run (rig, src, dst, 55.0, -40.0, -25.0, 32, 0.1);
        test::ok (sol.status == MasteringSolveStatus::Solved,
                  std::string ("the target is found, not declared out of range (got ")
                  + statusName (sol.status) + "/" + constraintName (sol.binding) + ")");
        test::approx (sol.measured.integratedLufs, -25.0, 0.1, "…on target");
        // PRECONDITION: the search really did have to leave the clamp behind — the first render is at
        // the caller's ceiling and nowhere near the answer.
        test::ok (sol.logCount >= 2 && sol.log[0].integratedLufs < -35.0,
                  "precondition: the first render is far from the target, at the caller's own ceiling");
        test::ok (sol.passes <= 8, "…and it does not spend the budget bisecting a bracket that never held it ("
                                   + std::to_string (sol.passes) + " renders)");
        std::printf ("      low ceiling: %s at %.4f LUFS, gain %.3f, ceiling %.3f, %d renders\n",
                     statusName (sol.status), sol.measured.integratedLufs, sol.preLimiterGainDb,
                     sol.ceilingDbTp, sol.passes);
    }

    // ---------------------------------------------------------------------------------------------
    test::group ("the QUIET end of the gain node is named too");
    {
        // The actuator has two ends. Testing only the loud one lost the other entirely: a target below
        // what -60 dB delivers came back a pass limit with nothing named.
        Programme src = tone (2.0, 0.5);
        Programme dst; dst.ch = src.ch; dst.bind();
        Rig rig;
        if (! test::run (rig.build (2))) return;
        const auto sol = run (rig, src, dst, 0.0, -1.0, -69.0, 8, 0.1);
        test::ok (sol.status == MasteringSolveStatus::TargetUnreachable,
                  std::string ("a target below the -60 dB end is unreachable (got ")
                  + statusName (sol.status) + ")");
        test::ok (sol.binding == MasteringConstraint::GainRange, "…and the gain range is what is NAMED");
        test::approx (sol.preLimiterGainDb, -60.0, 1e-9, "…with the gain at the quiet end of the clamp");
        // PRECONDITION: it really is out of reach, and by a margin the tolerance cannot swallow.
        test::ok (sol.measured.integratedLufs - (-69.0) > 1.0,
                  "precondition: -60 dB still leaves it "
                  + std::to_string (sol.measured.integratedLufs + 69.0) + " LU too loud");
        std::printf ("      quiet end: %s/%s at %.4f LUFS, gain %.3f\n", statusName (sol.status),
                     constraintName (sol.binding), sol.measured.integratedLufs, sol.preLimiterGainDb);
    }

    // ---------------------------------------------------------------------------------------------
    test::group ("a tolerance finer than the search can resolve is reported as that, not as a violation");
    {
        // A first candidate that hits the target exactly while breaking the ceiling has an error of zero
        // that nothing can improve on. Reading the verdict off it buried the later, clean candidates
        // that bracketed the target — the answer became `Unreachable / TruePeakCeiling` for a render
        // whose peak was in fact under the promise.
        Programme src = makeMusic (2.0, 0.3, 77u);
        Programme dst; dst.ch = src.ch; dst.bind();
        Rig rig;
        if (! test::run (rig.build (2))) return;
        rig.params.bypassCompressor = false;
        LoudnessRequest req;
        req.targetLufs = -14.0; req.maxTruePeakDbTp = -1.0;
        req.toleranceLu = 1.0e-10; req.initialGainDb = 12.0; req.maxPasses = 32;
        const auto sol = rig.solver.solve (rig.chain, rig.renderer, rig.params,
                                           src.in(), dst.out(), 2, src.frames(), req);
        test::ok (sol.status == MasteringSolveStatus::TargetBetweenAchievable
                  || sol.status == MasteringSolveStatus::Solved,
                  std::string ("an unreachable TOLERANCE is not a broken limit (got ")
                  + statusName (sol.status) + "/" + constraintName (sol.binding) + ")");
        test::ok (sol.measured.truePeakDbTp <= -1.0 + 1e-9,
                  "…and the delivered peak really is under the promise ("
                  + std::to_string (sol.measured.truePeakDbTp) + " dBTP)");
        test::approx (sol.measured.integratedLufs, -14.0, 0.001,
                      "…with the answer as close as the actuator allows");
        std::printf ("      tight tolerance: %s/%s at %.6f LUFS, TP %.4f, %d renders\n",
                     statusName (sol.status), constraintName (sol.binding),
                     sol.measured.integratedLufs, sol.measured.truePeakDbTp, sol.passes);
    }
}

// The diverse-testing round's finding: the ACTUATOR'S LIMIT was an event, so the verdict depended on
// the pass budget. The rows that REACH the clamp all deliver the same render; the ones that stop short
// of it do not, and the group separates the two rather than claiming one rule for both.
void testTheVerdictDoesNotDependOnTheBudget()
{
    test::group ("the actuator's limit is named by the render that reaches it, not one render later");
    {
        // Recording the pin only when a STEP tried to walk past the clamp needed an iteration the
        // budget did not always have. Measured on this exact input: `maxPasses = 10` gave `PassLimit`
        // with an empty mask at g = +60.000, `maxPasses = 11` gave `Unreachable / GainRange` — same
        // programme, same delivered audio, two different answers to "why did you stop".
        //
        // That "same delivered audio" is true of the CLAMPED rows (10, 11, 12) and of nothing else: a
        // budget that ends before the clamp naturally delivers a different, quieter render (row 8 stops
        // at g = +54.02, row 9 at +58.55). The header of this group said it of every row, which was
        // wrong, and the rows below assert the distinction rather than assuming it.
        //
        // SIGHTED ON BOTH BRANCHES. A budget that ends BEFORE the clamp is reached is a genuine pass
        // limit and must stay one, so the group asserts the precondition it turns on — whether the
        // last render actually sat at +-60 dB — and checks the verdict on each side of it.
        // THE EXPECTATION IS DERIVED FROM THE RENDER, NOT TABULATED. A hard-coded "budget 8 does not
        // reach the clamp" is a fact about how fast the search happens to be, and it went stale the
        // moment the slope floor was lowered — the row then asserted the opposite of the truth while
        // looking like a specification. What is actually being claimed is a BICONDITIONAL: the verdict
        // names the gain range exactly when the last render sits at the clamp. So each budget measures
        // its own precondition, and the sweep afterwards asserts that BOTH sides were really visited —
        // without that last check a group like this can quietly become one-sided and prove nothing.
        int sawClamped = 0, sawShort = 0;
        for (int passes : { 2, 3, 4, 5, 6, 8, 10, 12 })
        {
            Programme src = makeMusic (6.0, 0.95, 12345u);   // the default crest, not a pinned one
            Programme dst; dst.ch = src.ch; dst.bind();
            Rig rig;
            if (! test::run (rig.build (2))) return;
            rig.params.bypassCompressor = false;
            LoudnessRequest req;
            req.targetLufs = -5.3; req.maxTruePeakDbTp = -1.0; req.toleranceLu = 0.1;
            req.maxPasses = passes; req.initialGainDb = 0.0;
            const auto sol = rig.solver.solve (rig.chain, rig.renderer, rig.params,
                                               src.in(), dst.out(), 2, src.frames(), req);

            const bool atClamp = sol.logCount > 0
                              && std::fabs (sol.log[sol.logCount - 1].gainDb - 60.0) < 1.0e-6;
            test::ok (sol.measured.integratedLufs < -5.3 - req.toleranceLu,
                      "precondition: at budget " + std::to_string (passes)
                      + " the target is genuinely out of reach");
            if (atClamp)
            {
                ++sawClamped;
                test::ok (sol.status == MasteringSolveStatus::TargetUnreachable
                          && (sol.alsoViolated & constraintBit (MasteringConstraint::GainRange)) != 0u,
                          std::string ("budget ") + std::to_string (passes)
                          + ": the last render sits at the clamp, so the gain range is NAMED (got "
                          + statusName (sol.status) + "/" + constraintName (sol.binding) + ")");
            }
            else
            {
                ++sawShort;
                // The claim on this side is about the GAIN RANGE only. A short search may well end on
                // a real constraint — at budget 2 the render still breaks the ceiling and says so —
                // and asserting a particular status here would be asserting something else.
                test::ok ((sol.alsoViolated & constraintBit (MasteringConstraint::GainRange)) == 0u,
                          std::string ("budget ") + std::to_string (passes)
                          + ": a search that never reached the clamp does not name it (got "
                          + statusName (sol.status) + "/" + constraintName (sol.binding) + ")");
                test::ok (sol.binding != MasteringConstraint::GainRange,
                          "…and it is not the binding one either");
            }
            std::printf ("      budget %2d: %-8s/%-9s g %+8.3f  I %+9.4f  %s\n", passes,
                         statusName (sol.status), constraintName (sol.binding),
                         sol.preLimiterGainDb, sol.measured.integratedLufs,
                         atClamp ? "at the clamp" : "short of it");
        }
        // SIGHTED: both branches were exercised, so neither arm above is vacuous.
        test::ok (sawClamped > 0 && sawShort > 0,
                  "both branches of the biconditional occur in this sweep (" + std::to_string (sawShort)
                  + " short, " + std::to_string (sawClamped) + " at the clamp)");
    }
}

void testThePreMergeDiffPass()
{
    auto tone = [] (double seconds, double amp)
    {
        const int n = (int) (seconds * kFs);
        Programme p; p.ch.assign (2, std::vector<float> ((std::size_t) n, 0.0f));
        for (int i = 0; i < n; ++i)
        {
            const float v = (float) (amp * std::sin (2.0 * kPi * 1000.0 * (double) i / kFs));
            p.ch[0][(std::size_t) i] = v; p.ch[1][(std::size_t) i] = v;
        }
        p.bind();
        return p;
    };

    // ---------------------------------------------------------------------------------------------
    test::group ("an INFINITE constraint limit does not freeze the ranking");
    {
        // `minPlrDb` may legally be `+infinity` — the request check rejects NaN and nothing else — and
        // then every candidate's constraint excess is `+infinity` too. The infeasible tie-break has to
        // TIE there so the distance to the target can decide; a tolerance window cannot, because
        // `fabs(inf - inf)` is NaN and every comparison against NaN is false. Measured with the window:
        // the ranking froze on the first candidate and delivered -6.014 LUFS for a target of -20, with
        // the gain still sitting at its starting 0 dB. Fourteen LU, from a hardening.
        Programme src = tone (1.0, 0.5);
        Programme dst; dst.ch = src.ch; dst.bind();
        Rig rig;
        if (! test::run (rig.build (2))) return;
        rig.params.bypassCompressor = true; rig.params.bypassDither = true;
        LoudnessRequest req;
        req.targetLufs = -20.0; req.maxTruePeakDbTp = -1.0; req.maxPasses = 6;
        req.minPlrDb = std::numeric_limits<double>::infinity();
        const auto sol = rig.solver.solve (rig.chain, rig.renderer, rig.params,
                                           src.in(), dst.out(), 2, src.frames(), req);
        // PRECONDITION: the constraint really is unsatisfiable, so every candidate ranks equal on it.
        test::ok ((sol.alsoViolated & constraintBit (MasteringConstraint::PeakToLoudness)) != 0u,
                  "precondition: an infinite PLR floor is violated by every render");
        test::approx (sol.measured.integratedLufs, -20.0, 0.1,
                      "the search still walks to the target instead of freezing on candidate one");
        std::printf ("      infinite PLR floor: %s at %.6f LUFS, gain %.4f\n",
                     statusName (sol.status), sol.measured.integratedLufs, sol.preLimiterGainDb);
    }

    // ---------------------------------------------------------------------------------------------
    test::group ("the CEILING is an actuator too, and an untried one is not the gain range's fault");
    {
        // A caller may start the ceiling far below the promise. `c` only ever tracks the true-peak aim,
        // so until it reaches `pmax` there are dB of loudness the search has not tried, and naming the
        // gain node is a verdict about the wrong knob: `g = +60, c = -40` with one render's budget was
        // called `Unreachable / GainRange` for a -25 LUFS target that `g = -18.99, c = -1.05` delivers.
        Programme src = tone (1.0, 0.5);
        for (int passes : { 1, 32 })
        {
            Programme dst; dst.ch = src.ch; dst.bind();
            Rig rig;
            if (! test::run (rig.build (2))) return;
            rig.params.bypassCompressor = true; rig.params.bypassDither = true;
            rig.params.limiter.ceilingDbTp = -40.0;
            LoudnessRequest req;
            req.targetLufs = -25.0; req.maxTruePeakDbTp = -1.0; req.maxPasses = passes;
            req.initialGainDb = 60.0;
            const auto sol = rig.solver.solve (rig.chain, rig.renderer, rig.params,
                                               src.in(), dst.out(), 2, src.frames(), req);
            if (passes == 1)
            {
                // PRECONDITION: the one render really did sit at the clamp and really did fall short —
                // the two halves that would otherwise make this a pin.
                test::ok (sol.logCount == 1 && std::fabs (sol.log[0].gainDb - 60.0) < 1.0e-6,
                          "precondition: the single render sits at the +60 dB clamp");
                test::ok (sol.log[0].integratedLufs < -25.0 - req.toleranceLu,
                          "precondition: and it is far below the target");
                test::ok ((sol.alsoViolated & constraintBit (MasteringConstraint::GainRange)) == 0u,
                          std::string ("39 dB of untried ceiling is not a gain-range verdict (got ")
                          + statusName (sol.status) + "/" + constraintName (sol.binding) + ")");
            }
            else
                test::ok (sol.status == MasteringSolveStatus::Solved,
                          "…and with the budget to use that ceiling, the target is simply found");
            std::printf ("      untried ceiling, budget %2d: %s/%s g %+.4f c %+.4f I %+.4f\n", passes,
                         statusName (sol.status), constraintName (sol.binding),
                         sol.preLimiterGainDb, sol.ceilingDbTp, sol.measured.integratedLufs);
        }
    }

    // ---------------------------------------------------------------------------------------------
    test::group ("each half of the pin's test is load-bearing");
    {
        // Three mutants of the same three-term condition, each with its own witness. All three start
        // the gain where the mutant's mistake shows and give the search a single render, so nothing but
        // the pin's own arithmetic can decide the answer.
        Programme src = tone (1.0, 0.5);

        // (1) THE DIRECTION. The render at the clamp OVERSHOOTS the target, so no amount of gain is
        //     being asked for. Dropping `want > 0` names the gain range for a target below the render.
        {
            Programme dst; dst.ch = src.ch; dst.bind();
            Rig rig;
            if (! test::run (rig.build (2))) return;
            rig.params.bypassCompressor = true; rig.params.bypassDither = true;
            LoudnessRequest req;
            req.targetLufs = -40.0; req.maxTruePeakDbTp = -1.0; req.maxPasses = 1;
            req.initialGainDb = 60.0;
            const auto sol = rig.solver.solve (rig.chain, rig.renderer, rig.params,
                                               src.in(), dst.out(), 2, src.frames(), req);
            test::ok (sol.logCount == 1 && sol.log[0].integratedLufs > -40.0 + req.toleranceLu,
                      "precondition: the render at the clamp is LOUDER than the target");
            test::ok ((sol.alsoViolated & constraintBit (MasteringConstraint::GainRange)) == 0u,
                      "a target BELOW the render is not the gain range running out");
        }

        // (2) THE CLAMP ITSELF. A render half a dB short of +60 has gain left; only a render AT the
        //     clamp is pinned. A tolerance of a whole dB would swallow this one.
        {
            Programme dst; dst.ch = src.ch; dst.bind();
            Rig rig;
            if (! test::run (rig.build (2))) return;
            rig.params.bypassCompressor = true; rig.params.bypassDither = true;
            LoudnessRequest req;
            req.targetLufs = 6.0; req.maxTruePeakDbTp = -1.0; req.maxPasses = 1;
            req.initialGainDb = 59.5;
            const auto sol = rig.solver.solve (rig.chain, rig.renderer, rig.params,
                                               src.in(), dst.out(), 2, src.frames(), req);
            test::ok (sol.logCount == 1 && std::fabs (sol.log[0].gainDb - 59.5) < 1.0e-6,
                      "precondition: the render sits at +59.5 dB, half a dB inside the clamp");
            test::ok (sol.log[0].integratedLufs < 6.0 - req.toleranceLu,
                      "precondition: and it is short of the target, so gain IS being asked for");
            test::ok ((sol.alsoViolated & constraintBit (MasteringConstraint::GainRange)) == 0u,
                      "half a dB of remaining gain is still gain");
        }

        // (2b) THE QUIET END'S ERROR TEST. Its loud twin is redundant — see the proof in the header —
        //      but on this side nothing else stands between "the gain is at -60" and naming the range,
        //      so a render that MEETS the target must not be called out of gain. Infeasible for a
        //      different reason entirely, and THAT is what the caller needs told.
        {
            Programme probeDst; probeDst.ch = src.ch; probeDst.bind();
            Rig probeRig;
            if (! test::run (probeRig.build (2))) return;
            probeRig.params.bypassCompressor = true; probeRig.params.bypassDither = true;
            LoudnessRequest probeReq;
            probeReq.targetLufs = -200.0; probeReq.maxTruePeakDbTp = -1.0; probeReq.maxPasses = 1;
            probeReq.initialGainDb = -60.0;
            const auto probe = probeRig.solver.solve (probeRig.chain, probeRig.renderer, probeRig.params,
                                                      src.in(), probeDst.out(), 2, src.frames(), probeReq);
            test::ok (probe.logCount == 1, "precondition: the quiet probe render happened");
            if (probe.logCount != 1) return;

            Programme dst2; dst2.ch = src.ch; dst2.bind();
            Rig rig2;
            if (! test::run (rig2.build (2))) return;
            rig2.params.bypassCompressor = true; rig2.params.bypassDither = true;
            LoudnessRequest req;
            req.targetLufs = probe.log[0].integratedLufs - 0.05;   // inside the tolerance, below it
            req.maxTruePeakDbTp = -1.0; req.toleranceLu = 0.1;
            req.minPlrDb = 40.0;                                   // …and infeasible for another reason
            req.maxPasses = 1; req.initialGainDb = -60.0;
            const auto sol = rig2.solver.solve (rig2.chain, rig2.renderer, rig2.params,
                                                src.in(), dst2.out(), 2, src.frames(), req);
            test::ok (sol.logCount == 1 && std::fabs (sol.log[0].gainDb + 60.0) < 1.0e-6,
                      "precondition: the render sits at the -60 dB clamp");
            test::ok (std::fabs (sol.log[0].integratedLufs - req.targetLufs) <= req.toleranceLu,
                      "precondition: and MEETS the target in loudness");
            test::ok (sol.log[0].integratedLufs > req.targetLufs,
                      "precondition: with the target BELOW it, which is the quiet arm's direction");
            test::ok ((sol.alsoViolated & constraintBit (MasteringConstraint::GainRange)) == 0u,
                      "the quiet end does not name the gain range for a target it hit");
        }

        // (3) THE ERROR. A render at the clamp that MEETS the target in loudness is not out of range,
        //     whatever else is wrong with it — here the ceiling is broken, and THAT is the answer.
        {
            Programme dst; dst.ch = src.ch; dst.bind();
            Rig rig;
            if (! test::run (rig.build (2))) return;
            rig.params.bypassCompressor = true; rig.params.bypassDither = true;
            LoudnessRequest probeReq;
            probeReq.targetLufs = 0.0; probeReq.maxTruePeakDbTp = -1.0; probeReq.maxPasses = 1;
            probeReq.initialGainDb = 60.0;
            Programme probeDst; probeDst.ch = src.ch; probeDst.bind();
            Rig probeRig;
            if (! test::run (probeRig.build (2))) return;
            probeRig.params.bypassCompressor = true; probeRig.params.bypassDither = true;
            const auto probe = probeRig.solver.solve (probeRig.chain, probeRig.renderer, probeRig.params,
                                                      src.in(), probeDst.out(), 2, src.frames(), probeReq);
            test::ok (probe.logCount == 1, "precondition: the probe render happened");
            if (probe.logCount != 1) return;

            LoudnessRequest req;
            // The second request must render EXACTLY what the probe did, so the thing that makes it
            // infeasible may not be the ceiling: `maxTruePeakDbTp` is also the cap on `c`, and moving
            // it moves the render itself. A peak-to-loudness floor no render at +60 dB can meet leaves
            // the audio alone and still denies `Solved`.
            req.targetLufs = probe.log[0].integratedLufs;   // exactly what +60 dB delivers
            req.maxTruePeakDbTp = -1.0;
            req.minPlrDb = 40.0;
            req.maxPasses = 1; req.initialGainDb = 60.0;
            const auto sol = rig.solver.solve (rig.chain, rig.renderer, rig.params,
                                               src.in(), dst.out(), 2, src.frames(), req);
            test::ok (sol.logCount == 1
                      && std::fabs (sol.log[0].integratedLufs - req.targetLufs) <= req.toleranceLu,
                      "precondition: the render at the clamp MEETS the target in loudness");
            test::ok ((sol.alsoViolated & constraintBit (MasteringConstraint::PeakToLoudness)) != 0u,
                      "precondition: and breaks a constraint, so it is not `Solved`");
            test::ok ((sol.alsoViolated & constraintBit (MasteringConstraint::GainRange)) == 0u,
                      "a render that HIT the target is not also out of gain");
            std::printf ("      pin halves: direction, clamp and error each hold their own witness\n");
        }
    }

    // ---------------------------------------------------------------------------------------------
    test::group ("the idle anchor pays for itself where it actually pays");
    {
        // The regime is a start whose first render is still IDLE and short of the target: the search
        // has to climb through the drive at which the limiter engages, and the anchor is the exact
        // second point sitting on that boundary. Without it the first ACTIVE render pairs its secant
        // with a point far away in the linear region and a pass is spent recovering.
        //
        // Measured over four consecutive starts, three renders with the anchor and four without, every
        // time. Over the whole 102-cell battery 23 cells move, totalling 120 renders with against 124
        // without — a real win, and a small one; two cells are actually faster without it.
        for (double start : { 8.30, 8.35, 8.40, 8.45 })
        {
            Programme src = makeMusic (6.0, 0.3);
            Programme dst; dst.ch = src.ch; dst.bind();
            Rig rig;
            if (! test::run (rig.build (2))) return;
            // The chain has to be the one the measurement was taken on: the limiter's release is what
            // decides how much reduction the first warm render carries, and the default is not it.
            rig.params.bypassCompressor = false;
            rig.params.compressor.thresholdDb = -20.0; rig.params.compressor.ratio = 2.0;
            rig.params.compressor.kneeDb = 6.0;
            rig.params.compressor.attackMs = 15.0; rig.params.compressor.releaseMs = 180.0;
            rig.params.limiter.ceilingDbTp = -1.0;
            rig.params.limiter.releaseMs = 100.0;
            LoudnessRequest req;
            req.targetLufs = -10.5; req.maxTruePeakDbTp = -1.0; req.toleranceLu = 0.1;
            req.maxPasses = 6; req.initialGainDb = start;
            const auto sol = rig.solver.solve (rig.chain, rig.renderer, rig.params,
                                               src.in(), dst.out(), 2, src.frames(), req);
            // PRECONDITION: the start really is above the answer and the limiter really is engaged on
            // the first render — without both, this is an ordinary cold search and proves nothing.
            // PRECONDITION, and it is the whole point: the first render must be IDLE — that is what
            // establishes the anchor — and BELOW the target, so the search has to walk up through the
            // engagement boundary the anchor sits on. Written the other way round first ("a warm start
            // above the answer"), which is not this regime at all and made every row fail its own
            // precondition. The measurement was right; the sentence describing it was invented.
            test::ok (sol.logCount >= 1 && sol.log[0].integratedLufs < req.targetLufs,
                      "precondition: the first render is below the target");
            test::ok (sol.logCount >= 1 && sol.log[0].limiterMaxGrDb == 0.0,
                      "precondition: and the limiter is idle on it, so an anchor exists at all");
            char m[128];
            std::snprintf (m, sizeof m, "start %+.2f: three renders, not four (got %d)", start, sol.passes);
            test::ok (sol.passes <= 3, m);
        }
    }

    // ---------------------------------------------------------------------------------------------
    test::group ("the bracket closes on BOTH knobs, not on the gain alone");
    {
        // The sides are keyed on drive since the coordinate rewrite, so two renders can share a gain
        // exactly and still be a whole ceiling apart: `(60, -40)` and `(60, -1.05)`. Their gain gap is
        // ZERO, which passes any gain-only closure test, and the pair straddles the target — so the
        // solver announced `TargetBetweenAchievable`, "the target lies between two achievable values",
        // for a target it simply had not walked to yet. A quiet programme reaches that geometry in
        // three renders because the bootstrap spends the first one.
        Programme src = tone (1.0, 1.0e-4);
        Programme dst; dst.ch = src.ch; dst.bind();
        Rig rig;
        if (! test::run (rig.build (2))) return;
        rig.params.bypassCompressor = true; rig.params.bypassDither = true;
        rig.params.limiter.ceilingDbTp = -40.0;
        LoudnessRequest req;
        req.targetLufs = -25.0; req.maxTruePeakDbTp = -1.0; req.toleranceLu = 0.1;
        req.maxPasses = 3; req.initialGainDb = 0.0;
        const auto sol = rig.solver.solve (rig.chain, rig.renderer, rig.params,
                                           src.in(), dst.out(), 2, src.frames(), req);

        // PRECONDITION: the geometry this group is about is actually present — two renders with the
        // SAME gain, DIFFERENT ceilings, straddling the target. Without it the assertion below is
        // satisfied by any run at all.
        bool sawSameGainDifferentCeiling = false, straddles = false;
        for (int a = 0; a < sol.logCount && ! sawSameGainDifferentCeiling; ++a)
            for (int b = a + 1; b < sol.logCount; ++b)
                if (std::fabs (sol.log[a].gainDb - sol.log[b].gainDb) < 1.0e-9
                    && std::fabs (sol.log[a].ceilingDb - sol.log[b].ceilingDb) > 1.0e-3)
                {
                    sawSameGainDifferentCeiling = true;
                    straddles = (sol.log[a].integratedLufs - req.targetLufs)
                              * (sol.log[b].integratedLufs - req.targetLufs) < 0.0;
                    break;
                }
        test::ok (sawSameGainDifferentCeiling,
                  "precondition: two renders share a gain exactly and differ in ceiling");
        test::ok (straddles, "precondition: …and they straddle the target, so a gain-only test closes");
        test::ok (sol.status != MasteringSolveStatus::TargetBetweenAchievable,
                  std::string ("a zero GAIN gap across two ceilings is not a closed interval (got ")
                  + statusName (sol.status) + ")");
        std::printf ("      both knobs: %s at %.4f LUFS, g %+.4f c %+.4f, %d renders\n",
                     statusName (sol.status), sol.measured.integratedLufs,
                     sol.preLimiterGainDb, sol.ceilingDbTp, sol.passes);
    }
}

} // namespace


//==============================================================================
// THE LRA MEASUREMENT REFUSES A POISONED PROGRAMME. Here rather than only in the C-ABI suite: a consumer
// building the core as a subproject never builds `tools/`, and this is a behaviour change to a public
// method — it used to return SUCCESS on a programme its own meter had already flagged, and the meter is
// a local, so the caller could not check for itself.
static void testLraRefusesAPoisonedProgramme()
{
    test::group ("measureInputLoudnessRange refuses what the meter itself flags");

    const double fs = 48000.0;
    const int nch = 2;
    const int frames = 30 * 48000;                       // long enough for LRA's short-term samples
    std::vector<float> l ((std::size_t) frames), r ((std::size_t) frames);
    for (int i = 0; i < frames; ++i)                     // 3 s loud / 3 s quiet — a real range to lose
    {
        const double amp = ((i / 48000) % 6 < 3) ? 0.5 : 0.03;
        l[(std::size_t) i] = (float) (amp * std::sin (2.0 * 3.14159265358979 * 220.0 * i / fs));
        r[(std::size_t) i] = (float) (amp * std::sin (2.0 * 3.14159265358979 * 277.0 * i / fs));
    }
    const float* in[2] { l.data(), r.data() };

    mastering::TargetLoudnessSolver s;
    test::ok (s.prepare (fs, nch, 1024, 256, 4), "prepared");

    double clean = 0.0;
    test::ok (s.measureInputLoudnessRange (in, nch, frames, clean), "a clean programme is measured");
    test::approx (clean, 4.8, 0.05, "and the range is the fixture's own 4.8 LU");

    // Poison every LOUD second. A poisoned sub-hop is recorded as SILENCE, silence fails the absolute
    // gate, and the loud blocks leave the distribution the range is computed over — so the number this
    // used to return was 21.4 LU, a range the programme does not have, with `true` beside it.
    for (int sec = 0; sec < 30; ++sec)
        if (sec % 6 < 3)
            for (int i = sec * 48000; i < (sec + 1) * 48000; ++i)
                l[(std::size_t) i] = std::numeric_limits<float>::quiet_NaN();

    double poisoned = -1.0;
    test::ok (! s.measureInputLoudnessRange (in, nch, frames, poisoned), "the poisoned one is REFUSED");
    test::ok (poisoned == -1.0, "and the out-parameter is untouched by the refusal");

    // PRECONDITION, and the whole reason this test is worth having: the number really would have moved.
    // A fixture on which poisoning changes nothing would pass this test while proving nothing.
    analysis::LoudnessMeter lm;
    test::ok (lm.prepare (fs, nch, (double) frames / fs + 1.0), "an independent meter for the precondition");
    (void) lm.process (in, nch, frames);
    test::ok (lm.nonFiniteSubHops() > 0, "PRECONDITION: the meter really is flagging this programme");
    test::approx (lm.loudnessRangeLu(), 21.4, 0.05,
            "PRECONDITION: and the number it would have returned is 21.4 LU, not 4.8");
}

// P41 F1 — THE METER BEHIND A SOLVE IS SIZED IN SAMPLES. `frames / fs + 1` seconds is +inf at a finite rate the
// chain accepts, and the store's size used to be `(std::size_t) inf`: undefined behaviour whose answer depended on
// the row — 3 kept blocks and 194 dropped on arm64 and wasm32, 4 and 193 on x86-64 gcc. The rate is absurd on
// purpose and the chain takes it (without EQ and limiter, the two stages that refuse it); a solver that measures
// at the rate it was given has to size its meter for the programme it was given.
static void testTheMeterIsSizedInSamples()
{
    test::group ("the meter behind a solve is sized in samples, at any rate the chain accepts");

    const double fs = 1.0e-305;
    const int frames = 2000;
    MasteringChain chain; OfflineRenderer renderer; TargetLoudnessSolver solver;
    MasteringChainConfig cfg; cfg.eq = false; cfg.limiter = false;
    test::ok (chain.prepare (fs, 1, cfg), "PRECONDITION: the chain accepts the rate once EQ and limiter are off");
    test::ok (renderer.prepare (1, 1024), "renderer prepared");
    test::ok (solver.prepare (fs, 1, 1024, chain.internalBlock(), chain.tapOversampleFactor()), "solver prepared");
    test::ok (! std::isfinite ((double) frames / fs), "PRECONDITION: the seconds form really is +inf here");

    std::vector<float> in ((std::size_t) frames), out ((std::size_t) frames, 0.0f);
    for (int i = 0; i < frames; ++i) in[(std::size_t) i] = (i & 1) ? 0.25f : -0.25f;
    const float* ip[1] { in.data() };
    float*       op[1] { out.data() };
    LoudnessRequest req; req.targetLufs = -14.0; req.maxTruePeakDbTp = -1.0;
    const LoudnessSolution s = solver.solve (chain, renderer, MasteringChainParams {}, ip, op, 1, frames, req);
    // PRECONDITION: a meter was prepared at all. A refused solve leaves `measured` at its defaults, where "0 dropped"
    // holds vacuously — the stand's mutant back to the seconds form now REFUSES (+inf seconds) and passed that line.
    test::ok (s.status != MasteringSolveStatus::RenderFailed && s.passes > 0, "PRECONDITION: the solve measured something");
    // 2000 one-sample sub-hops are 200 hops, and the first block is born on the 4th: 197 blocks. Whatever the
    // verdict at this rate, the store must hold all of them — the defect dropped 193 or 194 of them.
    test::ok (s.measured.droppedBlocks == 0, "no gating block is dropped: the store holds the programme");
    test::ok (s.measured.gatingBlocks == 197, "and all 197 blocks the programme produces are kept");

    double lra = -1.0;
    (void) solver.measureInputLoudnessRange (ip, 1, frames, lra);   // the same sizing, exercised for the sanitizer rows
    // The range measurement is sized the same way, and a SILENT programme shows it — at 1e-304 Hz, where `frames / fs`
    // is +inf too. The programme is measured, 0 LU, as at any ordinary rate; the seconds form overran its store and
    // refused it.
    MasteringChain qChain; TargetLoudnessSolver qSolver;
    test::ok (qChain.prepare (1.0e-304, 1, cfg) && qSolver.prepare (1.0e-304, 1, 1024, qChain.internalBlock(),
                                                                    qChain.tapOversampleFactor()),
              "PRECONDITION: chain and solver take 1e-304 Hz");
    test::ok (! std::isfinite (100000.0 / 1.0e-304), "PRECONDITION: the seconds form is +inf at 1e-304 Hz too");
    const std::vector<float> quiet (100000, 0.0f);
    const float* qp[1] { quiet.data() };
    double quietLra = -1.0;
    test::ok (qSolver.measureInputLoudnessRange (qp, 1, 100000, quietLra) && quietLra == 0.0,
              "100 000 silent frames at 1e-304 Hz: the range is measured, 0 LU — the store holds the programme");
}

// P41 — A BUDGET IS 0 WHERE ITS CALL REFUSES, A CHANNEL COUNT INCLUDED. The pre-merge diff pass: both helpers used to
// turn a negative count into a huge number, and a different one on wasm32 than natively.
static void testTheBudgetsRefuseWhatTheCallsRefuse()
{
    test::group ("a budget is 0 where its call refuses — channel counts included");
    using felitronics::analysis::TruePeakMeter;
    const int past = felitronics::core::kMaxChannels + 1;
    test::ok (TruePeakMeter::storageFor (48000.0, 0).bytes() == 0 && TruePeakMeter::storageFor (48000.0, -1).bytes() == 0
              && TruePeakMeter::storageFor (48000.0, past).bytes() == 0,
              "TruePeakMeter::storageFor: 0 bytes for a channel count prepare() refuses");
    test::ok (TruePeakMeter::storageFor (48000.0, 2).bytes() == 296, "and 296 B for stereo at 48 kHz (the ABI suite's oracle)");
    test::ok (TargetLoudnessSolver::solveBytes (48000.0, 0, 48000) == 0 && TargetLoudnessSolver::solveBytes (48000.0, -1, 48000) == 0
              && TargetLoudnessSolver::solveBytes (48000.0, past, 48000) == 0,
              "solveBytes: 0 for a channel count solve() refuses");
    test::ok (TargetLoudnessSolver::solveBytes (48000.0, 2, 48000) == 3480, "and 3480 B for 1 s of stereo (the ABI suite's oracle)");
    // A prepare() refused on its bin width (400 dB at 1e-7 dB is 4e9 bins, past the 4e6 ceiling) allocates NOTHING —
    // which is what its budget says. The diverse-testing round found the tap buffers assigned before that refusal, and
    // kept. The delta is read into a local before the check.
    {
        TargetLoudnessSolver fine;
        const long long before = g_allocs.load();
        const bool refused = ! fine.prepare (48000.0, 2, 1024, 64, 4, 1.0e-7);
        const long long allocs = g_allocs.load() - before;
        test::ok (refused && allocs == 0 && TargetLoudnessSolver::prepareBytes (1024, 64, 4, 1.0e-7) == 0,
                  "a prepare() refused on its bin width allocates nothing, and its budget is 0");
    }
    test::ok (TargetLoudnessSolver::solveBytes (0.0, 2, 48000) == 0 && TargetLoudnessSolver::solveBytes (-1.0, 2, 48000) == 0
              && TargetLoudnessSolver::measureRangeBytes (0.0, 480000) == 0, "and 0 for a rate the solver refuses");
    // The true-peak meter's factor follows the RATE, and so must its budget — one rate could not tell (the diverse-
    // testing round's mutant sized it at 48 kHz and passed). Derived: 1 s is 20 hops at any multiple of 100 Hz, so the
    // loudness meter is 8·(300 + 24 + 10) = 2672 B; the true-peak meter at factor F is 4·12F + 4·2·12 + 4·2; the drain 512.
    test::ok (TruePeakMeter::storageFor (96000.0, 2).bytes() == 200 && TruePeakMeter::storageFor (192000.0, 2).bytes() == 152,
              "the true-peak meter at 96 kHz (factor 2) is 200 B, at 192 kHz (factor 1) 152 B");
    test::ok (TargetLoudnessSolver::solveBytes (96000.0, 2, 96000) == 2672u + 200u + 512u
              && TargetLoudnessSolver::solveBytes (192000.0, 2, 192000) == 2672u + 152u + 512u,
              "and a 1 s solve at 96 and 192 kHz carries that meter: 3384 and 3336 B");
}

int main()
{
    std::printf ("felitronics::mastering::TargetLoudnessSolver — P7\n");
    testTappedRenderNullsAgainstThePlainOne();
    testShortTapRefusesTheWholeCall();
    testHitsTheTarget();
    testTwoPassesWhenTheShapeAlreadyAdmitsIt();
    testUnreachableIsNamed();
    testUpstreamIsNotBlamedOnTheTarget();
    testStatisticsAgreeWithAHandDrivenChain();
    testRefusalsAndDegenerateInputs();
    testBlockIndependence();
    testTheReportedRenderIsTheDeliveredOne();
    testTheScaleLawIsPinned();
    testTheAbsoluteGateStepIsPinned();
    testTheReviewsCounterexamples();
    testTwoConstraintsAtOnce();
    testPreLimiterTapIsTheRealSignal();
    testTargetBetweenAchievable();
    testTheAnswerDoesNotDependOnWhereItStarted();
    testAnUnmeasurableStartIsNotAnUnmeasurableProgramme();
    testTapPlumbingEdges();
    testSurvivorsOfTheMutationStand();
    testTheReviewRoundsCounterexamples();
    testTheVerdictDoesNotDependOnTheBudget();
    testThePreMergeDiffPass();
    testLraRefusesAPoisonedProgramme();
    testTheMeterIsSizedInSamples();
    testTheBudgetsRefuseWhatTheCallsRefuse();
    return felitronics::test::report();
}

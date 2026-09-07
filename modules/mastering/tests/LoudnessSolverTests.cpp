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
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

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
        case MasteringSolveStatus::GateStep:           return "GateStep";
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
        Programme src = makeMusic (24.0, 0.28);
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
    Programme src = makeMusic (20.0, 0.20);
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
        Programme src = makeWideRange (30.0);
        Programme dst; dst.ch = src.ch; dst.bind();
        Rig rig;
        if (! test::run (rig.build (2))) return;

        const bool haveLra = rig.solver.measureInputLra (src.in(), 2, src.frames());
        test::ok (haveLra, "the input's LRA is measurable (the constraint is a DELTA and needs it)");
        // PRECONDITION: the fixture has a range to lose. Without this the constraint cannot bind and
        // the test would pass against a solver that never checked it.
        test::ok (rig.solver.inputLoudnessRangeLu() > 6.0,
                  "precondition: the fixture's own LRA is wide (measured "
                  + std::to_string (rig.solver.inputLoudnessRangeLu()) + " LU)");

        LoudnessRequest req;
        req.targetLufs = -8.0;                 // loud enough to need heavy limiting
        req.maxTruePeakDbTp = -1.0;
        req.maxLraLossLu = 0.5;                // the corpus number, rounded up
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
        const double loss = rig.solver.inputLoudnessRangeLu() - ind.LRA;
        test::ok (loss > req.maxLraLossLu,
                  "precondition: the range loss really exceeds the allowance ("
                  + std::to_string (loss) + " LU against " + std::to_string (req.maxLraLossLu) + ")");
        test::ok (sol.measured.limiter.maxDb <= 0.0,
                  "precondition: the LIMITER did nothing at this render, so the loss is upstream of it");
        std::printf ("      LRA case: in %.2f LU -> out %.2f LU, I %.3f (target %.1f), status %s, binding %s\n",
                     rig.solver.inputLoudnessRangeLu(), ind.LRA, sol.measured.integratedLufs,
                     req.targetLufs, statusName (sol.status), constraintName (sol.binding));
    }

    // (b) THE LIMITER'S GAIN REDUCTION, on a MANUFACTURED dense input — the plan's own recipe, because a
    // dense mix with PLR ~10 does not exist as an input (10 is the RESULT of mastering, measured on two
    // proven pairs). Manufacturing one tests the REFUSAL, not the sound.
    {
        Programme src = makeMusic (16.0, 0.62);
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
    Programme src = makeMusic (12.0, 0.5);
    Programme dst; dst.ch = src.ch; dst.bind();
    Rig rig;
    if (! test::run (rig.build (2))) return;
    // Make the compressor work hard, then forbid it. The pre-limiter gain node cannot undo this at ANY
    // setting, because it sits AFTER the compressor.
    rig.params.compressor.thresholdDb = -36.0;
    rig.params.compressor.ratio = 8.0;
    LoudnessRequest req;
    req.targetLufs = -14.0;
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
    Programme src = makeMusic (10.0, 0.35);
    Programme dst; dst.ch = src.ch; dst.bind();
    Rig rig;
    if (! test::run (rig.build (2))) return;
    LoudnessRequest req;
    req.targetLufs = -12.0;
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
    const long long limFrom = (long long) r.compressorLookahead + (long long) r.clipperLatency;
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
    test::ok (hl.maxValue() > 0.05, "precondition: the limiter really limited ("
                                    + std::to_string (hl.maxValue()) + " dB peak)");
    test::ok (nc == (std::uint64_t) frames, "compressor window: exactly `frames` samples counted");
    test::ok (nl == (std::uint64_t) frames * (std::uint64_t) F, "limiter window: exactly frames*F counted");

    double p95c = 0.0, p95l = 0.0;
    test::ok (hc.quantile (0.95, p95c), "hand-driven compressor p95 is answerable");
    test::ok (hl.quantile (0.95, p95l), "hand-driven limiter p95 is answerable");
    test::approx (sol.measured.compressor.meanDb, hc.mean(), 1.0e-12, "compressor mean nulls");
    test::approx (sol.measured.compressor.p95Db,  p95c,      1.0e-12, "compressor p95 nulls");
    test::approx (sol.measured.compressor.maxDb,  hc.maxValue(), 1.0e-12, "compressor max nulls");
    test::approx (sol.measured.compressor.activeFraction, (double) ac / (double) nc, 1.0e-12,
                  "compressor active fraction nulls");
    test::approx (sol.measured.limiter.meanDb, hl.mean(), 1.0e-12, "limiter mean nulls");
    test::approx (sol.measured.limiter.p95Db,  p95l,      1.0e-12, "limiter p95 nulls");
    test::approx (sol.measured.limiter.maxDb,  hl.maxValue(), 1.0e-12, "limiter max nulls");
    test::approx (sol.measured.limiter.activeFraction, (double) al / (double) nl, 1.0e-12,
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
    Programme src = makeMusic (4.0, 0.4);
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
        LoudnessRequest req;
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

        LoudnessRequest req2; req2.maxPasses = 0;
        const auto sol2 = rig.solver.solve (rig.chain, rig.renderer, rig.params,
                                            src.in(), dst.out(), 2, src.frames(), req2);
        test::ok (sol2.status == MasteringSolveStatus::InvalidRequest, "zero pass budget refused");

        LoudnessRequest req3;
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
        LoudnessRequest req; req.targetLufs = -14.0;
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
    Programme src = makeMusic (8.0, 0.3);
    double firstI = 0.0, firstG = 0.0;
    bool have = false;
    for (int blk : { 64, 256, 1000, 8192 })
    {
        Programme dst; dst.ch = src.ch; dst.bind();
        Rig rig;
        if (! test::run (rig.build (2, blk))) return;
        LoudnessRequest req; req.targetLufs = -13.0; req.maxPasses = 3;
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
    Programme src = makeWideRange (24.0);
    Programme dst; dst.ch = src.ch; dst.bind();
    Rig rig;
    if (! test::run (rig.build (2))) return;
    (void) rig.solver.measureInputLra (src.in(), 2, src.frames());
    LoudnessRequest req;
    req.targetLufs = -7.0;
    req.maxLraLossLu = 0.5;
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

} // namespace

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
    return felitronics::test::report();
}

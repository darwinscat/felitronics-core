// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/dynamics/offline/DetectorAnalysis.h>

#include <cmath>
#include <limits>
#include <vector>

namespace felitronics::dynamics::offline
{

// Why the inverse problem can fail to have an answer, NAMED. A plausible number with a boolean beside
// it is the shape that gets ignored; a reason is the shape that gets handled.
enum class SolveStatus
{
    Solved,                 // a threshold was found whose measured amount is within tolerance
    ToleranceBelowResolution,  // the search closed its bracket and the statistic is still outside
                            // `toleranceDb` — which, the objective being CONTINUOUS in the threshold
                            // (see the note below), can only mean the tolerance asked for is finer than
                            // the measurement can resolve. It is a converged search reporting the limit
                            // of its own instrument, not a failure to converge, and not a property of
                            // the material. If it appears at a tolerance comfortably above
                            // `binDb + multiplier * resolutionDb`, the bracket is wrong — which is
                            // exactly how it was used to find a bracket bug.
    FloorDrivesTheAnswer,   // the answer processes frames sitting on `core::gainToDb`'s -240 dB clamp,
                            // i.e. it is steering by a FLOOR rather than by a measurement. Legitimate
                            // arithmetic, degenerate as advice, and true in every mode — it is asked of
                            // the curve, not inferred from which mode is running.
    NotPrepared,
    EmptyInput,             // no frames, or none survived the statistics gate
    InvalidTarget,          // amount non-finite or negative, quantile outside (0, 1]
    InvalidParams,          // non-finite curve parameters
    NoProcessingPossible,   // ratio == 1 (slope 0) or rangeDb == 0: the amount is identically zero
    TargetIsZero,           // degenerate but legitimate; the canonical no-processing threshold is returned
    TargetNotReachable,     // the most aggressive threshold in the bracket still falls short — see achievedDb
    QuantileOutOfRange,     // the measured quantile fell outside the histogram: no number to trust
    EvaluationLimit         // ran out of passes before the bracket closed
};

// WHICH FUNCTIONAL OF THE GAIN-REDUCTION TRACE IS BEING TARGETED. The search is valid for any of
// these and for any other that is monotone under POINTWISE DOMINATION of the trace — which is the
// actual content of the monotonicity argument, and is why the choice is a parameter rather than a
// hard-coded p95. It matters: two pieces of material with the same p95 can sound opposite (held at
// 3 dB all the time by a slow release, versus struck 5 % of the time), so the quantity a product
// should aim at is a product question, and a solver that only knew one would have to be rewritten to
// answer it. All three are in dB, so one tolerance covers them.
enum class TraceStatistic
{
    Quantile,   // the q-quantile of |gain change| — the classic "p95 GR"
    Mean,       // its arithmetic mean over the counted frames
    Max         // its maximum
};

struct ThresholdTarget
{
    TraceStatistic statistic = TraceStatistic::Quantile;
    double quantile      = 0.95;    // OF THE MAGNITUDE |gain change| — see DetectorAnalysis.h
    double amountDb      = 3.0;     // how much of it is wanted, as a positive number of dB
    double toleranceDb   = 0.05;    // accept when |achieved - amount| is within this
    double resolutionDb  = 1.0e-3;  // a CEILING on the bracket width the search settles for; a steep
                                    // curve tightens it further, because one dB of threshold is worth
                                    // `multiplier` dB of gain change
    int    maxEvaluations = 48;     // one evaluation is one full pass over the key
    double statisticsGateDb = -std::numeric_limits<double>::infinity();
    bool   useSeed       = true;    // bracket around the level-percentile estimate first
    double seedWindowDb  = 6.0;     // the half-width of that first bracket, doubled on each miss
};

struct ThresholdSolution
{
    double thresholdDb   = 0.0;     // the answer: a threshold VERIFIED to meet the target
    double achievedDb    = 0.0;     // the measured quantile of |gain change| at that threshold
    double bracketLoDb   = 0.0;     // the interval the search closed to; the answer's remaining slack
    double bracketHiDb   = 0.0;
    double seedDb        = 0.0;     // the level-percentile estimate, whether or not it was used
    double attainableDb  = 0.0;     // the most this material can be made to give at all
    double achievedBelowDb = 0.0;   // the statistic just the other side of the boundary: the size of
                                    // the step, and the only way to see a discontinuity for what it is
    int    evaluations   = 0;       // full passes over the key, the honest cost
    SolveStatus   status = SolveStatus::NotPrepared;
    EnvelopeStats envelope {};
};

//==============================================================================
// felitronics::dynamics::offline::ThresholdSolver — finds the threshold that delivers a wanted
// quantile of gain reduction on a given key, by running the real chain and measuring.
//
// THE SEARCH IS ON A PREDICATE, NOT ON A VALUE, and that is what decides the awkward cases for free.
// `f(threshold)` — the measured quantile of |gain change| — is monotone in the threshold (decreasing
// for DownCompress, increasing for UpCompress and DownExpand: the curve is monotone in
// `level - threshold`, and the follower preserves pointwise order, so the whole trace does and
// therefore so does every quantile of it). Bisecting the BOOLEAN `f(threshold) >= amount` converges to
// the boundary of the set that meets the target, so when the solution is not unique — a plateau, which
// a saturated range or a coarse material makes ordinary — the answer returned is the LEAST AGGRESSIVE
// threshold that still meets it. No tie-break rule to invent, no midpoint of a plateau to justify.
//
// THE SEED HAS TO BRACKET TO BE WORTH ANYTHING, and that is a fact about bisection rather than about
// this problem. Midpoint bisection is MINIMAX-OPTIMAL on the width of the interval, so no starting
// guess can improve its worst case; simply probing a seed instead of the midpoint splits the interval
// UNEVENLY and is, on average, worse — measured: one extra pass on all three test signals. What a seed
// can do is make the interval smaller to begin with. So the level-percentile estimate is used to open a
// narrow bracket around itself (`seedWindowDb`, doubling outward until it straddles the answer), and
// only then is the interval bisected. When the estimate is good that starts from 12 dB instead of 320
// and saves the log2 of the ratio; when it is bad the doubling walks out to the global bracket and the
// cost is the handful of probes that walk. Both outcomes are measured in the tests rather than assumed.
//
// THE OBJECTIVE IS ALSO CONTINUOUS IN THE THRESHOLD, and that is worth stating because it says which
// surprises are possible. The static curve is Lipschitz in `level - threshold` with constant
// `multiplier`; the follower is a convex combination and so is 1-Lipschitz in the pointwise supremum
// of its input; an order statistic is 1-Lipschitz in that same supremum. Composing,
// `|f(a) - f(b)| <= multiplier * |a - b|`. So a converged bracket of width `resolutionDb` pins the
// statistic to within `multiplier * resolutionDb` plus the histogram's half-bin either side — about
// 0.011 dB at the defaults — and there is no step for the search to fall into. Atoms in the material
// make the HISTOGRAM a step function of the value; they do not make `f` a step function of the
// threshold, and conflating the two is what once gave a real bracket bug an innocent-sounding name.
//
// MONOTONICITY IS EXACT IN ARITHMETIC AND APPROXIMATE IN FLOAT, and the contract says the second — but
// the size matters and the fixture has to be named. On ordinary music at 10/100 ms, sweeping the
// threshold in 1e-3 dB steps, the measured statistic never moves the wrong way at all. It takes an
// adversarial fixture to see it: 5 ms ballistics on a two-level signal (the F9 layout) at 1e-5 dB
// steps gives 2555 reversals in 4000 with a worst case of 2.7e-05 dB, and a 21-sample fixture with a
// 12.5 dB range gives 97 in 20000 at 9.5e-07 dB. The cause is rounding in the follower's
// `target + c*(current - target)` where `c` is close to 1. The default `resolutionDb` of 1e-3 sits an
// order of magnitude above the worst of that, and the search stops on a TOLERANCE and a RESOLUTION,
// never on an equality.
//
// EVERY PROBE STARTS FROM A RESET. A solve is an iteration, and between iterations the detector and
// the gain-reduction follower must not remember the previous threshold — `setParams()` is not a reset,
// and a solver that forgot this would return an answer that depends on the order the probes happened
// to be tried. Here the path is reset inside `evaluate()`, which is the only place a probe can happen.
//
// THE COST IS PASSES OVER THE KEY, and the number has to be the one a caller will actually see rather
// than the cost of the arithmetic inside it. Measured end to end — `solve()` timed, divided by the
// passes it reported — on a five-minute MONO 48 kHz tone with a 1.3 Hz amplitude envelope, Apple
// silicon, `-O2 -ffp-contract=on`: **319 ms per pass without the envelope cache, 189 ms with it**, so
// a 16-pass solve is 5.1 s or 3.0 s. The chain itself is only 73 ms of that; the rest is accumulating
// the statistic, one sample at a time. An earlier draft of this comment said 64 ms, which was the
// chain measured on its own — a real number for the wrong thing, and the reason the fixture is spelled
// out here. `begin()`/`step()` exist so a single-threaded host (the wasm tier has no other kind) can
// run ONE pass per turn and stay responsive, and `sol_.evaluations` never exceeds `maxEvaluations`.
// THE KEY IS BORROWED, not copied — it must outlive the solve.
class ThresholdSolver
{
public:
    // `maxCachedFrames` buys SPEED WITH MEMORY, and the trade is the caller's because only the caller
    // knows which is scarce. Zero (the default) keeps the solver O(1) in the length of the material and
    // re-runs the detector on every probe; a value covering the material records the detector envelope
    // once — four bytes a frame, 57.6 MB for five stereo minutes at 48 kHz — and every probe after the
    // first reads it. Measured on a five-minute track: 64 ms per pass re-running the detector against
    // 44 ms reading the cache. The two produce BIT-IDENTICAL traces by construction: both end in
    // `GainReductionPath::step()`, and there is no second implementation to drift.
    [[nodiscard]] bool prepare (double sampleRate, int maxCachedFrames = 0, double binDb = 0.01)
    {
        prepared_ = false;
        if (! (sampleRate > 0.0) || ! std::isfinite (sampleRate)) return false;
        if (maxCachedFrames < 0) return false;
        fs_ = sampleRate;
        if (! analyzer_.prepare (sampleRate, binDb)) return false;
        path_.prepare (sampleRate);
        // |gain change| is bounded by rangeDb, which GainComputer itself caps at 400 dB.
        if (! grHist_.prepare (0.0, 400.0, binDb)) return false;
        env_.assign ((std::size_t) maxCachedFrames, 0.0f);
        prepared_ = true;
        return true;
    }

    bool isPrepared() const noexcept { return prepared_; }

    ThresholdSolution solve (const float* const* key, int numKeyChannels, int numSamples,
                             const GainReductionParams& base, const ThresholdTarget& target)
    {
        if (! begin (key, numKeyChannels, numSamples, base, target)) return sol_;
        while (step()) {}
        return sol_;
    }

    // Sets the search up and does the fixed work: one envelope pass, both bracket endpoints measured.
    // Returns false when the answer is already decided (every status but `Solved` can be decided here).
    bool begin (const float* const* key, int numKeyChannels, int numSamples,
                const GainReductionParams& base, const ThresholdTarget& target)
    {
        sol_ = ThresholdSolution {};
        seeding_ = false; seedProbed_ = false; moreMoved_ = false; lessMoved_ = false;
        key_ = key; nKeyCh_ = numKeyChannels; n_ = numSamples; base_ = base; t_ = target;
        running_ = false;

        if (! prepared_)                                   { sol_.status = SolveStatus::NotPrepared;  return false; }
        if (key == nullptr || numKeyChannels < 1 || numSamples <= 0) { sol_.status = SolveStatus::EmptyInput; return false; }
        if (! std::isfinite (t_.amountDb) || t_.amountDb < 0.0
            || ! (t_.quantile > 0.0) || ! (t_.quantile <= 1.0)
            || ! std::isfinite (t_.toleranceDb) || t_.toleranceDb < 0.0
            || ! std::isfinite (t_.resolutionDb) || ! (t_.resolutionDb > 0.0)
            || t_.maxEvaluations < 1) { sol_.status = SolveStatus::InvalidTarget; return false; }
        if (! std::isfinite (base.ratio) || ! std::isfinite (base.kneeDb) || ! std::isfinite (base.rangeDb)
            || ! std::isfinite (base.attackMs) || ! std::isfinite (base.releaseMs))
                                                           { sol_.status = SolveStatus::InvalidParams; return false; }

        // --- one envelope pass: the bracket, the seed and the reported statistics all come from it ---
        analyzer_.setParams (base);
        analyzer_.setStatisticsGateDb (t_.statisticsGateDb);
        analyzer_.reset();
        cached_ = ((int) env_.size() >= numSamples);
        analyzer_.analyze (key, numKeyChannels, numSamples,
                           cached_ ? env_.data() : nullptr, (int) env_.size());
        sol_.envelope = analyzer_.stats();
        if (! sol_.envelope.valid)                         { sol_.status = SolveStatus::EmptyInput; return false; }

        // THE EFFECTIVE CURVE, read from a GainComputer configured exactly as the path's will be —
        // never re-derived. `setRangeDb` caps at 400 dB and `setRatio` floors at 1, so a solver doing
        // its own arithmetic on the raw fields diverges from the curve it is searching: dividing a raw
        // `rangeDb` of 1e308 by the slope overflowed the bracket to infinity, and the search then
        // returned a confident threshold of +Inf from entirely finite parameters.
        GainComputer curve;
        curve.setMode (base.mode); curve.setRatio (base.ratio);
        curve.setKneeDb (base.kneeDb); curve.setRangeDb (base.rangeDb);
        mult_ = curve.multiplier();
        const double rangeDb = curve.rangeDb();
        if (! (mult_ > 0.0) || ! (rangeDb > 0.0))          { sol_.status = SolveStatus::NoProcessingPossible; return false; }

        const double h  = curve.kneeDb() * 0.5;
        const double sat = (h > rangeDb / mult_ ? h : rangeDb / mult_) + 1.0;   // +1 dB of slack
        lowerIsMore_ = (base.mode == Mode::DownCompress);
        // THE BRACKET IS BUILT FROM THE LEVELS THAT OCCURRED, not from the ones that were counted. A
        // statistics gate removes frames from the DISTRIBUTION and never from the stream, so a gated
        // extreme is not the extreme the ballistics saw — measured, a bracket sized from the gated
        // pair failed to contain the answer at all, and the search then converged to the wrong end.
        // `more` is the endpoint that produces the most processing, `less` the one that produces none.
        more_ = lowerIsMore_ ? (sol_.envelope.minSeenDb - sat) : (sol_.envelope.maxSeenDb + sat);
        less_ = lowerIsMore_ ? (sol_.envelope.maxSeenDb + h + 1.0) : (sol_.envelope.minSeenDb - h - 1.0);
        if (! std::isfinite (more_) || ! std::isfinite (less_)) { sol_.status = SolveStatus::InvalidParams; return false; }
        sol_.seedDb = seedThreshold (h);

        // THE RESOLUTION HAS TO MATCH THE CURVE'S STEEPNESS. A dB of threshold is worth `mult_` dB of
        // gain change, so a fixed floor on the bracket width is a floor of `mult_ * width` on the
        // answer's error: at ratio 1000 in DownExpand the multiplier is 999, and a 1e-3 dB bracket
        // still left a whole dB of error — which the search then reported as a discontinuity that was
        // not there. Half a tolerance, converted through the multiplier, is the width that makes the
        // stated tolerance mean something; the caller's value is a ceiling on it, never a floor.
        resolution_ = t_.resolutionDb;
        const double needed = t_.toleranceDb / (2.0 * mult_);
        if (std::isfinite (needed) && needed < resolution_) resolution_ = needed;
        if (! (resolution_ > 1.0e-9)) resolution_ = 1.0e-9;

        if (! (t_.amountDb > 0.0))
        {
            // MEASURED, not asserted. "The gentlest threshold does nothing" is true by construction of
            // the endpoint — but a constant compared with itself is the shape of test that passes on
            // both sides of a behaviour change, and the same is true of a reported number that was
            // never taken. One pass is the whole cost.
            sol_.thresholdDb = less_;
            if (! evaluate (less_, sol_.achievedDb)) { sol_.status = SolveStatus::QuantileOutOfRange; return false; }
            sol_.bracketLoDb = sol_.bracketHiDb = less_;
            sol_.status = SolveStatus::TargetIsZero;
            return false;
        }

        // THE SATURATING ENDPOINT IS MEASURED, not assumed — that is the one that cannot be reasoned
        // about. The most a piece of material can give is not `rangeDb`: the follower starts at zero
        // and the passage may be too short or too sparse to charge it, so a target well inside the
        // range can still be unreachable, and that has to be its own answer rather than a failure to
        // converge. The other endpoint is not measured and does not need to be: at a threshold past
        // every level by half a knee, no sample reaches the curve at all, so the statistic there is
        // exactly zero by construction. Both are recorded below either way.
        double achievedMore = 0.0;
        if (! evaluate (more_, achievedMore))              { sol_.status = SolveStatus::QuantileOutOfRange; return false; }
        sol_.attainableDb = achievedMore;
        // THE PREDICATE IS `f >= amount - tolerance`, not `f >= amount`, and the two have to be the
        // same statement or the invariant is a fiction: accepting an endpoint that falls short by half
        // a tolerance while calling it "verified to meet the target" is the gap through which a Solved
        // could be returned for a target no threshold satisfies.
        if (achievedMore < t_.amountDb - t_.toleranceDb)
        {
            sol_.thresholdDb = more_; sol_.achievedDb = achievedMore;
            sol_.bracketLoDb = std::fmin (more_, less_); sol_.bracketHiDb = std::fmax (more_, less_);
            sol_.status = SolveStatus::TargetNotReachable;
            return false;
        }
        // Both endpoints are now KNOWN, and only one of them had to be measured: `less_` is the
        // threshold at which no sample reaches the knee at all, so its statistic is exactly zero by
        // construction rather than by measurement. Keeping both means `finish()` costs no further pass.
        achievedMore_ = achievedMore;
        achievedLess_ = 0.0;
        moreMoved_ = lessMoved_ = false;
        // NB when the saturated end already sits ON the target, nothing special happens: the search
        // below still walks it back to the gentlest threshold that holds, which is what searching a
        // predicate rather than a value buys.

        // The global bracket is now established: `more_` measures at or above the target, `less_` is
        // the no-processing end, which is zero by construction. Everything below only ever narrows it.
        seeding_ = t_.useSeed && std::isfinite (sol_.seedDb)
                   && sol_.seedDb > std::fmin (more_, less_) && sol_.seedDb < std::fmax (more_, less_);
        seedProbed_ = false;
        // Floored, so a degenerate window cannot make the doubling walk in place: at `denorm_min()`
        // the search spent its whole budget re-probing a single point.
        expand_ = (std::isfinite (t_.seedWindowDb) && t_.seedWindowDb > 0.0) ? t_.seedWindowDb : 6.0;
        if (expand_ < 4.0 * resolution_) expand_ = 4.0 * resolution_;

        running_ = true;
        return true;
    }

    // One bisection step = one full pass over the key. Returns false when the search is over.
    bool step()
    {
        if (! running_) return false;
        if (std::fabs (less_ - more_) <= resolution_)     { finish (true);  return false; }
        if (sol_.evaluations >= t_.maxEvaluations)        { finish (false); return false; }

        const double probe = nextProbe();
        double achieved = 0.0;
        if (! evaluate (probe, achieved)) { running_ = false; sol_.status = SolveStatus::QuantileOutOfRange; return false; }
        if (achieved >= t_.amountDb - t_.toleranceDb) { more_ = probe; achievedMore_ = achieved; moreMoved_ = true; }
        else                                          { less_ = probe; achievedLess_ = achieved; lessMoved_ = true; }
        return true;
    }

    const ThresholdSolution& solution() const noexcept { return sol_; }

    // The naive inversion, exposed so the acceptance can MEASURE how far it misses rather than argue
    // about it: the same detector, the same curve, the same quantile — with the ballistics removed.
    // That isolates the temporal error, which is the only thing bisection buys over it.
    double staticInversionThresholdDb (const GainReductionParams& base, double levelDb, double amountDb) const noexcept
    {
        const double h = (base.kneeDb > 0.0 ? base.kneeDb : 0.0) * 0.5;
        return invertCurve (base, levelDb, amountDb, h);
    }

private:
    // `closed` distinguishes the two ways the loop can end, and the distinction is the whole point:
    // a search that closed its bracket to `resolutionDb` HAS converged, and if the statistic is still
    // outside tolerance that is because the statistic STEPS — a histogram is a step function, and so is
    // a quantile of a trace with atoms in it. Reporting that as "did not converge" names the wrong
    // cause and sends the caller looking for a bug in the search.
    void finish (bool closed)
    {
        running_ = false;
        // NO FURTHER PASSES. Both endpoints' measurements were kept as they were made, so the answer
        // costs nothing to report — which is what makes `step()` honestly one pass and keeps the total
        // inside `maxEvaluations` instead of two over it. The answer is the endpoint KNOWN to satisfy
        // the predicate; `less_` was only ever shown to fall short, and reporting its value too is what
        // makes the size of a step visible rather than merely suspected.
        const double achieved = achievedMore_, below = achievedLess_;
        sol_.thresholdDb     = more_;
        sol_.achievedDb      = achieved;
        sol_.achievedBelowDb = below;
        sol_.bracketLoDb     = std::fmin (more_, less_);
        sol_.bracketHiDb     = std::fmax (more_, less_);

        // IS THE ANSWER STEERING BY THE FLOOR? Asked of the curve at the returned threshold, not
        // inferred from the mode: a frame at `gainToDb`'s -240 dB clamp is not a measurement, so an
        // answer that only works because such frames are being processed is arithmetic rather than
        // advice. A mode-shaped test got this wrong in exactly the ways a mode-shaped test does —
        // UpCompress happily returned Solved at -236 dB, lifting digital silence by the requested 3 dB.
        if (sol_.envelope.atFloorFrames > 0)
        {
            GainComputer c;
            c.setMode (base_.mode); c.setRatio (base_.ratio); c.setKneeDb (base_.kneeDb);
            c.setRangeDb (base_.rangeDb); c.setThresholdDb (more_);
            if (std::fabs (c.deltaDb (EnvelopeAnalyzer::kFloorDb)) > 0.0)
            { sol_.status = SolveStatus::FloorDrivesTheAnswer; return; }
        }
        if (std::fabs (achieved - t_.amountDb) <= t_.toleranceDb) { sol_.status = SolveStatus::Solved; return; }
        sol_.status = closed ? SolveStatus::ToleranceBelowResolution : SolveStatus::EvaluationLimit;
    }

    // Where to look next. While seeding, walk outward from the estimate until BOTH sides of the answer
    // have been found; after that, the midpoint, which is optimal on what is then a narrow interval.
    double nextProbe() noexcept
    {
        const double dirMore = lowerIsMore_ ? -1.0 : 1.0;    // the way that increases processing
        if (seeding_)
        {
            if (! seedProbed_) { seedProbed_ = true; return sol_.seedDb; }
            // One side is still sitting on the starting endpoint: step outward on that side, doubling.
            // Tracked with flags rather than by comparing doubles for equality — the equality would be
            // exact here, but it is the shape `-Wfloat-equal` exists to catch, and a flag says what is
            // meant.
            if (! moreMoved_ || ! lessMoved_)
            {
                const double d = (! moreMoved_) ? dirMore : -dirMore;
                const double probe = sol_.seedDb + d * expand_;
                expand_ *= 2.0;
                const double loB = std::fmin (more_, less_), hiB = std::fmax (more_, less_);
                if (! (probe > loB) || ! (probe < hiB)) { seeding_ = false; return 0.5 * (more_ + less_); }
                return probe;
            }
            seeding_ = false;                                // both sides found: bisect the narrow one
        }
        return 0.5 * (more_ + less_);
    }

    // One probe: the WHOLE chain over the WHOLE key, from a reset, into a fresh histogram.
    [[nodiscard]] bool evaluate (double thresholdDb, double& achievedDb)
    {
        GainReductionParams p = base_;
        p.thresholdDb = thresholdDb;
        path_.setParams (p);
        path_.reset();                                    // setParams() is NOT a reset — see the note above
        grHist_.reset();
        ++sol_.evaluations;

        // The gate reads the LEVEL, which does not move with the threshold, so the counted set is
        // identical at every probe and the objective stays monotone. Off is the common case and costs
        // nothing: no level below the -240 dB floor exists, so a gate there or lower excludes nothing.
        // `>=`, NOT `>`: the analyzer excludes a frame whose level is at or below the gate, so a gate
        // sitting exactly ON the -240 dB floor still excludes the floor frames. Reading that as "off"
        // made the solver measure a different set of frames from the analyzer that sized its bracket —
        // measured, a solver reporting 3.005 dB where the gated truth was 6.011.
        const bool gating = (t_.statisticsGateDb >= EnvelopeAnalyzer::kFloorDb);
        for (int i = 0; i < n_; ++i)
        {
            const float gr = cached_ ? path_.processLevel (env_[(std::size_t) i])
                                     : path_.process (key_, nKeyCh_, i);
            if (! gating || core::gainToDb (cached_ ? env_[(std::size_t) i] : path_.detectorLevel()) > t_.statisticsGateDb)
                grHist_.add (std::fabs ((double) gr));
        }
        path_.flushDenormals();
        if (grHist_.count() == 0) { achievedDb = 0.0; return true; }
        switch (t_.statistic)
        {
            case TraceStatistic::Mean:     achievedDb = grHist_.mean();     return true;
            case TraceStatistic::Max:      achievedDb = grHist_.maxValue(); return true;
            case TraceStatistic::Quantile: break;
        }
        return grHist_.quantile (t_.quantile, achievedDb);
    }

    // The seed: invert the STATIC curve at the level quantile. Cheap (no signal), and only ever a
    // starting point — see the counterexample in DetectorAnalysis.h for why it cannot be the answer.
    double seedThreshold (double h) const noexcept
    {
        double lq = 0.0;
        if (! analyzer_.histogram().quantile (t_.quantile, lq)) lq = sol_.envelope.p95Db;
        return invertCurve (base_, lq, t_.amountDb, h);
    }

    // Solve |deltaDb(levelDb ; threshold)| == amountDb for the threshold. Bisection on the curve alone
    // — it costs nothing and handles the soft knee, the range clamp and every mode without a special
    // case each.
    static double invertCurve (const GainReductionParams& base, double levelDb, double amountDb, double h) noexcept
    {
        GainComputer gc;
        gc.setMode (base.mode); gc.setRatio (base.ratio); gc.setKneeDb (base.kneeDb); gc.setRangeDb (base.rangeDb);
        const bool lowerIsMore = (base.mode == Mode::DownCompress);
        const double mult = (base.mode == Mode::DownExpand)
                                ? (base.ratio - 1.0) : (1.0 - 1.0 / (base.ratio > 1.0 ? base.ratio : 1.0));
        if (! (mult > 0.0)) return levelDb;
        const double sat = (h > base.rangeDb / mult ? h : base.rangeDb / mult) + 1.0;
        double more = lowerIsMore ? levelDb - sat : levelDb + sat;
        double less = lowerIsMore ? levelDb + h + 1.0 : levelDb - h - 1.0;
        for (int i = 0; i < 60; ++i)
        {
            const double mid = 0.5 * (more + less);
            gc.setThresholdDb (mid);
            if (std::fabs (gc.deltaDb (levelDb)) >= amountDb) more = mid; else less = mid;
        }
        return more;
    }

    EnvelopeAnalyzer   analyzer_;
    GainReductionPath  path_;
    QuantileHistogram  grHist_;
    std::vector<float> env_;        // the optional detector-envelope cache — see prepare()
    ThresholdSolution  sol_ {};
    ThresholdTarget    t_ {};
    GainReductionParams base_ {};

    const float* const* key_ = nullptr;
    int    nKeyCh_ = 0, n_ = 0;
    double fs_ = 48000.0, more_ = 0.0, less_ = 0.0, mult_ = 0.0;
    double expand_ = 6.0, resolution_ = 1.0e-3, achievedMore_ = 0.0, achievedLess_ = 0.0;
    bool   prepared_ = false, running_ = false, lowerIsMore_ = true, cached_ = false;
    bool   seeding_ = false, seedProbed_ = false, moreMoved_ = false, lessMoved_ = false;
};

} // namespace felitronics::dynamics::offline

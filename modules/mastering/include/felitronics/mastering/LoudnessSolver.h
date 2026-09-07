// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/analysis/LoudnessMeter.h>
#include <felitronics/analysis/TruePeakMeter.h>
#include <felitronics/dynamics/offline/Quantile.h>
#include <felitronics/mastering/MasteringChain.h>
#include <felitronics/mastering/OfflineRenderer.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

namespace felitronics::mastering
{

//==============================================================================
// felitronics::mastering::TargetLoudnessSolver — hit a target integrated loudness with a stated true-peak
// ceiling, in a bounded number of renders, and REFUSE rather than crush the programme when the target
// cannot be had without breaking a named limit.
//
// ====================================================================================
// THE ONE EQUATION — why this is a scalar search and not two loops chasing each other
// ====================================================================================
// The solver owns exactly two numbers: `preLimiterGainDb` (call it g) and `limiter.ceilingDbTp` (c).
// They look like two knobs and are not, and this is the load-bearing fact of the whole class. Inside
// `limiter::TruePeakLimiter` the reduction is `rawRedDb = min(0, c - smaxDb)`, and `smaxDb` is taken on
// the signal AFTER the gain node, so `smaxDb(10^(g/20) p) = smaxDb(p) + g`. Therefore
//
//     rawRedDb = min (0, -(g - c) - smaxDb(p))
//
// depends on g and c ONLY through the difference. Write `d = g - c` (the DRIVE, in dB above the
// ceiling); the release recursion `grDb = min(rawRed, grDb*relCoef)` reads only `rawRed`, so the whole
// gain trace is a function of d alone, and the output is
//
//     y(g, c) = 10^(c/20) * y(d, 0)
//
// i.e. c is a pure output SCALE and d is the entire SHAPE. Measured on this tree over a 3x3 grid of
// (g, c): max |y(g,c) - 10^(c/20) y(d,0)| = 9.1e-07 .. 1.7e-06 on a programme peaking at 0.89, which is
// float rounding in `dbToGain` and the FIR, not a structural gap. Consequences, and they are the design:
//
//     I(g, c)  = c + J(d)          achieved integrated loudness
//     TP(g, c) = c + T(d)          achieved true peak
//     PLR(d)   = T(d) - J(d)       peak-to-loudness — INDEPENDENT OF c
//
// So "hit the target AND stay under the ceiling" is not a pair of loops. It is
//
//     find the smallest d with PLR(d) <= ceiling - target;   then  c = target - J(d),  g = d + c.
//
// `PLR` is monotone non-increasing in d (more drive removes peaks and adds loudness), measured on three
// real mixes: 16.53 -> 6.95, 13.17 -> 7.29 and 15.87 -> 7.34 over an 18 dB sweep, with no reversal.
// SMALLEST d is not an aesthetic preference either: d IS how hard the limiter works, so the least d
// that satisfies the ceiling is the most transient the programme can keep at that loudness.
//
// THE DERATE IS NOT A CONSTANT HERE, and that is the point. `TruePeakLimiter` bounds its own F*fs grid
// exactly and overshoots the reconstructed peak between grid points — up to +1.22 dB in the worst
// non-degenerate case at 4x, and measured IN SITU on the corpus at only +0.0005 .. +0.1028 dB, because
// spectral tilt buys it cheaply. A fixed 1.2 dB derate would cost loudness on every track. Here the
// delivered peak is MEASURED and enters `PLR`, so the derate is whatever the material's is, on that
// material, and the ceiling that gets programmed into the limiter is an OUTPUT of the solve.
//
// ====================================================================================
// WHAT A PASS IS, and the honest cost — because the interesting claim is a cost claim
// ====================================================================================
// A PASS is one `OfflineRenderer::render` over the whole programme. `passes` in the result counts them
// and nothing else. Measured on this tree (60 s stereo, one render each), a render is NOT uniformly
// expensive: the limiter is 1417.7 ms of 1522.8, i.e. **93 %**; the dither 4 %; the EQ and the
// compressor together 3.3 %; the chain with the limiter and the dither bypassed is **3.7 %**. So a
// probe that re-runs only the limiter is not cheap, and a scheme that saves chain renders while
// running the limiter saves nothing. This class therefore minimises RENDERS, states its probe cost
// openly, and reports both numbers.
//
// ====================================================================================
// WHY THE PREDICATE IS NOT `I == target`
// ====================================================================================
// `analysis::LoudnessMeter`'s integrated measure is gated, and the ABSOLUTE gate at -70 LUFS is not
// scale-invariant (the relative one is). A block crossing it changes the SET being averaged, so `I` is
// piecewise-continuous in gain rather than smooth, and the step can be large: measured, a fixture of
// two halves at -69 and -71 LUFS reads -0.879829 LU away from `I(0) + g` the moment the quieter half
// enters. On real programme the step is bounded by one block joining a set of N: 10*log10((N+0.1)/(N+1)),
// which is -0.002 LU at 3.5 minutes and -0.096 LU at 4 seconds — and a sweep of three real mixes over
// 0..18 dB moved no blocks at all (the counts were 3168, 1437 and 3983 at every gain). The design
// answer is therefore not to assume smoothness: every candidate is MEASURED, the best FEASIBLE one is
// kept, and a target that lands inside a step is reported as `GateStep` with the two sides — never as
// `Solved`, and never as a constraint violation, because it is neither.
//
// ====================================================================================
// WHAT IT WILL NOT DO
// ====================================================================================
// The gain node it drives sits AFTER the compressor (`MasteringChain`: gate -> inputGain -> EQ ->
// mono-bass -> compressor -> clipper -> preLimiterGain -> limiter -> dither). So the compressor's gain
// reduction, and most of what the compressor does to the loudness range, DO NOT MOVE with the search.
// A compressor limit that the caller's settings already break is broken at every gain, and calling it
// "the reason the target is unreachable" would be a lie: it is reported as `UpstreamViolation`, with the
// measurement, before the search spends a pass on it. Moving that would mean solving for the compressor
// threshold, which is `dynamics::offline::ThresholdSolver`'s job and a different question.
//
// OFFLINE. `prepare()` allocates; `solve()` does not, and it is called from a worker rather than an
// audio callback. It drives `OfflineRenderer`, which drives `MasteringChain`, so every measurement it
// makes is one the chain's fixed internal quantum has already made independent of block size —
// measured: the achieved integrated loudness of the same programme is bit-identical at renderer block
// sizes 1, 64, 256, 977, 4096 and 65536.
//==============================================================================

// Which limit stopped the search. NAMED, because a bare `unreachable` flag is a sneeze: the chain this
// replaces raised one while standing 0.1 LU from its target, and a flag without the name of the binding
// limit would be exactly as useless.
enum class MasteringConstraint
{
    None = 0,
    TruePeakCeiling,          // the target cannot be reached without delivering above `maxTruePeakDbTp`
    LimiterGainReduction,     // ... without the limiter exceeding `limiterGrLimitDb`
    PeakToLoudness,           // ... without the delivered PLR falling below `minPlrDb`
    LoudnessRange,            // ... without LRA falling more than `maxLraLossLu` below the input's
    GainRange,                // ... without a gain outside what the chain accepts (+-60 dB)
    CompressorGainReduction   // UPSTREAM: broken by the caller's settings at every gain — see above
};

inline constexpr std::uint32_t constraintBit (MasteringConstraint c) noexcept
{
    return (c == MasteringConstraint::None) ? 0u : (1u << ((int) c - 1));
}

enum class MasteringSolveStatus
{
    Solved,                 // measured within `toleranceLu` of the target, ceiling held, nothing bound
    TargetUnreachable,      // a NAMED constraint binds; the result carries the best FEASIBLE render
    UpstreamViolation,      // a limit the search cannot move is already broken — see the note above
    GateStep,               // the target falls inside a discontinuity of the gated measure; both
                            // sides are reported and neither is within tolerance. Honest unreachability,
                            // and NOT a constraint violation
    PassLimit,              // ran out of `maxPasses` while still converging. Not a verdict about the
                            // material: it is a verdict about the budget, and it says so
    MeasurementInvalid,     // the meter could not answer (no gating block, dropped blocks, non-finite)
    RenderFailed,           // the chain or the renderer refused a call
    NotPrepared,
    InvalidRequest
};

// The pass budget's ceiling, named here because `LoudnessSolution` carries an array of that length and
// a struct cannot reach into the class that uses it.
struct TargetLoudnessSolverLimits { static constexpr int kMaxPasses = 32; };

// Which statistic a gain-reduction limit binds. Part of the limit's TYPE, never a hidden convention:
// "max 6 dB of limiting" and "p95 under 6 dB" are different products, and a field named for one while
// enforcing the other is the shape of defect this repository keeps closing.
enum class GrStatistic { Mean, P95, Max };

struct GainReductionLimit
{
    double      limitDb  = std::numeric_limits<double>::infinity();   // infinity = no limit
    GrStatistic statistic = GrStatistic::Max;
};

// A gain-reduction trace, summarised. `|GR|` throughout — the traces are SIGNED (reduction is negative),
// and a mean over signed values is a different number that nobody wants.
struct GainReductionStats
{
    double meanDb        = 0.0;
    double p95Db         = 0.0;
    double maxDb         = 0.0;
    double activeFraction = 0.0;    // fraction of samples with |GR| > activityThresholdDb
    std::uint64_t frames = 0;
    std::uint64_t nonFinite = 0;    // a poisoned trace: the numbers above are best effort, not measurement
    std::uint64_t aboveRange = 0;   // values past the histogram's top — quantiles are then not answerable
    bool valid = false;             // false ⇒ every number above is a placeholder, not a measurement
};

// Everything the request asked to be told, from ONE render. The three achieved numbers are measured on
// exactly the delivered frames; the two gain-reduction summaries are measured on exactly the tap frames
// that carry programme (each stage's own window — see MasteringChainTaps).
struct MasterMeasurement
{
    double integratedLufs   = 0.0;
    double truePeakDbTp     = 0.0;      // DRAINED: see the note in measure()
    double samplePeakDb     = 0.0;
    double loudnessRangeLu  = 0.0;
    double plrDb            = 0.0;      // truePeakDbTp - integratedLufs
    GainReductionStats compressor {};
    GainReductionStats limiter {};
    double limiterMaxReconstructedPeakDb = 0.0;   // the peak the limiter's own oversampler saw
    int    latencySamples   = 0;
    int    gatingBlocks     = 0;
    int    droppedBlocks    = 0;        // non-zero ⇒ the loudness numbers describe a PREFIX
    int    nonFiniteSubHops = 0;
    bool   loudnessValid    = false;
    bool   lraValid         = false;    // LRA needs short-term samples; a short programme has none, and
                                        // `loudnessRangeLu() == 0` cannot tell "no range" from "no data"
};

struct LoudnessRequest
{
    double targetLufs      = -14.0;
    double toleranceLu     =   0.1;
    double maxTruePeakDbTp =  -1.0;     // DELIVERED, measured. Not the limiter's setting — that is derived.
    // How far BELOW the promise to aim the delivered peak. Not decoration and not taste: the ceiling
    // loop drives the delivered peak toward its aim, and an aim of exactly `maxTruePeakDbTp` converges
    // to the boundary FROM ABOVE and never crosses it — measured, a solve stalled at -1.0000 dBTP with
    // the ceiling constraint still flagged on three consecutive renders. The aim has only to exceed the
    // change in the limiter's own between-grid overshoot across one correction step, which is a few
    // thousandths of a dB; 0.05 is two decades of slack and is still an order below the 0.5 dB a
    // mastering engineer would notice.
    double truePeakAimDb   =   0.05;

    // Constraints. A target that needs one of these broken is REFUSED with the name, not forced through.
    GainReductionLimit limiterGr    {};                                              // off by default
    GainReductionLimit compressorGr {};                                              // UPSTREAM — see above
    double minPlrDb     = -std::numeric_limits<double>::infinity();                  // delivered TP - I
    // A DELTA, not an absolute floor, and the difference is the whole calibration: the two proven
    // input->accepted-master pairs in this project's corpus move LRA by -0.40 and -0.30 LU, while the
    // ffmpeg chain being replaced collapses it 17.4 -> 2.9. An absolute floor cannot tell those apart —
    // the catalogue's own inputs run from 3.1 to 14.2 LU — and a delta can.
    double maxLraLossLu = std::numeric_limits<double>::infinity();

    // What counts as "the limiter was working" / "the compressor was working". A THRESHOLD, not a
    // comparison with zero: a release from 6 dB decays for tens of thousands of samples before it
    // reaches exactly 0, so `|GR| > 0` reports a duty cycle of the release time constant rather than of
    // the programme. 0.1 dB is an order below the loudness tolerance and is reported back in the result.
    double activityThresholdDb = 0.1;

    // Renders. Two is enough whenever the correction is EXACT (see the step rule in solve()) and three
    // is what a shape change costs, because a shape change needs a slope and a slope needs two real
    // measurements. Four leaves one spare for a programme that also has to bring its ceiling down.
    int    maxPasses = 4;
    double initialGainDb = std::numeric_limits<double>::quiet_NaN();   // NaN = use the params' own
};

// One render the search made. The whole trace is returned, not just the winner: a caller that has to
// explain "why is this track at -12.3 and not -12.0" cannot do it from a single row, and neither can a
// test. It is also what makes the pass COUNT auditable rather than a number to be trusted.
struct SolvePassRecord
{
    double gainDb = 0.0, ceilingDb = 0.0;
    double integratedLufs = 0.0, truePeakDbTp = 0.0, plrDb = 0.0;
    double limiterMaxGrDb = 0.0, loudnessRangeLu = 0.0;
    std::uint32_t violated = 0;
};

struct LoudnessSolution
{
    MasteringSolveStatus status = MasteringSolveStatus::NotPrepared;
    MasteringConstraint  binding = MasteringConstraint::None;
    std::uint32_t        alsoViolated = 0;      // bitmask over constraintBit(); the binding one included

    double preLimiterGainDb = 0.0;              // what was applied to produce the delivered render
    double ceilingDbTp      = 0.0;              // ... and the ceiling the limiter was actually given
    MasterMeasurement measured {};              // of the DELIVERED render, always — never of a probe

    int    passes = 0;                          // renders spent
    double activityThresholdDb = 0.1;           // echoed, because a fraction without its threshold is not a number

    // GateStep only: the two sides of the discontinuity the target fell into.
    double achievedBelowLufs = 0.0, achievedAboveLufs = 0.0;
    double gainBelowDb = 0.0, gainAboveDb = 0.0;

    // Every render, in order. `passes` is its length plus any re-render of the winner.
    SolvePassRecord log[TargetLoudnessSolverLimits::kMaxPasses] {};
    int logCount = 0;
};

class TargetLoudnessSolver
{
public:
    static constexpr double kMaxGainDb = 60.0;      // MasteringChain::kMaxGainDb — the search's actuator range
    static constexpr int    kMaxPasses = TargetLoudnessSolverLimits::kMaxPasses;

    // `maxFrames` and `maxChannels` size the tap buffers; `binDb` is the resolution every gain-reduction
    // quantile is reported to. The tap buffers are the whole allocation and they are per RENDERER BLOCK,
    // not per programme — the traces are consumed as they arrive, so a five-minute track costs the same
    // as a five-second one.
    [[nodiscard]] bool prepare (double sampleRate, int maxChannels, int rendererBlock,
                                int internalBlock, int oversampleFactor, double binDb = 0.01)
    {
        prepared_ = false;
        if (! (sampleRate > 0.0) || ! std::isfinite (sampleRate)) return false;
        if (maxChannels < 1 || maxChannels > core::kMaxChannels) return false;
        if (rendererBlock < 1 || internalBlock < 1 || oversampleFactor < 1) return false;
        if (! (binDb > 0.0) || ! std::isfinite (binDb)) return false;
        // A `process(n)` call writes K * floor((pos + n) / K) tap frames, which is at most n + K - 1.
        const long long cap = (long long) rendererBlock + (long long) internalBlock;
        if (cap > (long long) std::numeric_limits<int>::max() / (oversampleFactor + 1)) return false;

        fs_ = sampleRate;
        nch_ = maxChannels;
        frameCap_ = (int) cap;
        osCap_    = (int) (cap * (long long) oversampleFactor);
        compTap_.assign ((std::size_t) frameCap_, 0.0f);
        limTap_.assign  ((std::size_t) osCap_, 0.0f);
        limPeak_.assign ((std::size_t) osCap_, 0.0f);
        // The gain-reduction range: `GainComputer` caps its own range at 400 dB, and the limiter's is
        // bounded by its ceiling clamp. 400 covers both, and anything past it is COUNTED rather than
        // folded into the top bin, so a quantile that lands there answers `false` instead of lying.
        if (! compHist_.prepare (0.0, 400.0, binDb)) return false;
        if (! limHist_.prepare  (0.0, 400.0, binDb)) return false;
        prepared_ = true;
        return true;
    }

    bool isPrepared() const noexcept { return prepared_; }

    // Render `frames` of `in` into `out` at a gain and ceiling chosen to meet `req`. `params` is the
    // caller's whole parameter set; the solver overrides exactly `preLimiterGainDb` and
    // `limiter.ceilingDbTp` and leaves every other field alone.
    //
    // `params` is taken by value rather than read back from the chain ON PURPOSE: `MasteringChain::params()`
    // reports the last APPLIED set, which lags a `setParams()` by up to one internal quantum, and it
    // stores the caller's UNCLAMPED request while the chain applies a clamped one. A search that read
    // its own actuator through either of those would be measuring a number it did not apply.
    LoudnessSolution solve (MasteringChain& chain, OfflineRenderer& renderer,
                            MasteringChainParams params,
                            const float* const* in, float* const* out,
                            int numChannels, int frames, const LoudnessRequest& req)
    {
        LoudnessSolution sol;
        sol.activityThresholdDb = req.activityThresholdDb;
        if (! prepared_) { sol.status = MasteringSolveStatus::NotPrepared; return sol; }
        if (! chain.isPrepared() || numChannels != chain.numChannels() || numChannels > nch_
            || frames <= 0 || in == nullptr || out == nullptr)
            { sol.status = MasteringSolveStatus::InvalidRequest; return sol; }
        if (! std::isfinite (req.targetLufs) || ! std::isfinite (req.maxTruePeakDbTp)
            || ! std::isfinite (req.toleranceLu) || req.toleranceLu < 0.0
            || ! std::isfinite (req.activityThresholdDb) || req.activityThresholdDb < 0.0
            || req.maxPasses < 1 || req.maxPasses > kMaxPasses)
            { sol.status = MasteringSolveStatus::InvalidRequest; return sol; }
        // The tap buffers were sized for a geometry; a chain that does not match them would be measured
        // through a refused call, which is a silent zero rather than a statistic.
        if (chain.internalBlock() + renderer.blockSize() > frameCap_
            || (long long) (chain.internalBlock() + renderer.blockSize()) * chain.tapOversampleFactor() > (long long) osCap_)
            { sol.status = MasteringSolveStatus::InvalidRequest; return sol; }

        const double target = req.targetLufs;
        const double pmax   = req.maxTruePeakDbTp;
        // The limiter's ceiling is a MEANS; `maxTruePeakDbTp` is the promise. Starting above the promise
        // would ship a render whose only guard is the solver's own arithmetic, which is the failure mode
        // P1 measured in the chain this replaces. It starts at the tighter of the two and only ever falls.
        const double c0     = std::min (pmax, std::isfinite (params.limiter.ceilingDbTp)
                                                  ? params.limiter.ceilingDbTp : pmax);
        double g = std::isfinite (req.initialGainDb) ? req.initialGainDb : params.preLimiterGainDb;
        if (! std::isfinite (g)) g = 0.0;
        double c = c0;

        Best best;                          // the best FEASIBLE render seen, and the best of any kind
        double loG = 0.0, hiG = 0.0, loI = 0.0, hiI = 0.0;   // the two sides that bracket the target
        bool haveLo = false, haveHi = false;
        double prevG = 0.0, prevI = 0.0;    // the PREVIOUS render, for a local secant
        bool   havePrev = false;
        // THE IDLE ANCHOR, and it is a measurement rather than a model. While the limiter does not
        // engage, the chain from the gain node on is a plain multiply, so `I(g) = I1 + (g - g1)` holds
        // EXACTLY up to the gain at which the limiter starts working — and that gain is `c` minus the
        // reconstructed peak the limiter itself reports. So a single idle render hands the search a
        // second exact point sitting ON the boundary of the active region, which is the anchor a local
        // secant wants: the first active render then pairs with it instead of with a point far away in
        // the linear region. Measured, that pairing is worth a whole pass on a loud target.
        double anchorG = 0.0, anchorI = 0.0;
        bool   haveAnchor = false;
        const double aim = pmax - (std::isfinite (req.truePeakAimDb) && req.truePeakAimDb > 0.0
                                       ? req.truePeakAimDb : 0.0);

        for (int pass = 0; pass < req.maxPasses; ++pass)
        {
            g = std::clamp (g, -kMaxGainDb, kMaxGainDb);
            c = std::clamp (c, -kMaxGainDb, kMaxGainDb);
            params.preLimiterGainDb      = g;
            params.limiter.ceilingDbTp   = c;
            chain.setParams (params);

            MasterMeasurement m;
            if (! renderPass (chain, renderer, params, in, out, numChannels, frames, req, m))
                { sol.status = MasteringSolveStatus::RenderFailed; sol.passes = pass + 1; return sol; }
            ++sol.passes;

            if (! m.loudnessValid)
            {
                sol.status = MasteringSolveStatus::MeasurementInvalid;
                sol.measured = m; sol.preLimiterGainDb = g; sol.ceilingDbTp = c;
                return sol;
            }

            // UPSTREAM first, and before the search spends anything on it: the compressor's gain
            // reduction does not move with `preLimiterGainDb`, so a broken limit there is broken at
            // every gain, and reporting it as "the target is unreachable" would name the wrong thing.
            if (violates (m.compressor, req.compressorGr))
            {
                sol.status  = MasteringSolveStatus::UpstreamViolation;
                sol.binding = MasteringConstraint::CompressorGainReduction;
                sol.alsoViolated |= constraintBit (MasteringConstraint::CompressorGainReduction);
                sol.measured = m; sol.preLimiterGainDb = g; sol.ceilingDbTp = c;
                return sol;
            }

            const std::uint32_t viol = violatedMask (m, req);

            // UPSTREAM, PART TWO — and this one is not about the compressor's own limit. The search
            // makes a programme LOUDER by driving the limiter harder, and every one of these three gets
            // WORSE with drive: more limiting means more gain reduction, less peak-to-loudness and less
            // loudness range. So a limit already broken at the least drive the search will use is broken
            // at every drive it can reach, and naming it "the target is unreachable" points the user at
            // the target when the fault is in the chain's settings. Measured on the corpus: with a 0.4 LU
            // range allowance, this chain's own compressor spends 0.80 to 3.70 LU before the solver
            // applies a single dB, and every one of those tracks came back blaming the loudness target.
            // The guard is the FIRST render only, and only when the target needs MORE drive than it — a
            // quieter target unwinds the drive and can cure all three.
            if (pass == 0 && target > m.integratedLufs)
            {
                const std::uint32_t worsensWithDrive =
                    constraintBit (MasteringConstraint::LimiterGainReduction)
                  | constraintBit (MasteringConstraint::PeakToLoudness)
                  | constraintBit (MasteringConstraint::LoudnessRange);
                const std::uint32_t up = viol & worsensWithDrive;
                if (up != 0)
                {
                    sol.status  = MasteringSolveStatus::UpstreamViolation;
                    sol.binding = bindingOf (up);
                    sol.alsoViolated = viol;
                    sol.measured = m; sol.preLimiterGainDb = g; sol.ceilingDbTp = c;
                    if (sol.logCount < kMaxPasses)
                    {
                        SolvePassRecord& rec0 = sol.log[sol.logCount++];
                        rec0.gainDb = g; rec0.ceilingDb = c;
                        rec0.integratedLufs = m.integratedLufs; rec0.truePeakDbTp = m.truePeakDbTp;
                        rec0.plrDb = m.plrDb; rec0.limiterMaxGrDb = m.limiter.maxDb;
                        rec0.loudnessRangeLu = m.loudnessRangeLu; rec0.violated = viol;
                    }
                    return sol;
                }
            }

            if (sol.logCount < kMaxPasses)
            {
                SolvePassRecord& rec = sol.log[sol.logCount++];
                rec.gainDb = g; rec.ceilingDb = c;
                rec.integratedLufs = m.integratedLufs; rec.truePeakDbTp = m.truePeakDbTp;
                rec.plrDb = m.plrDb; rec.limiterMaxGrDb = m.limiter.maxDb;
                rec.loudnessRangeLu = m.loudnessRangeLu; rec.violated = viol;
            }
            const bool onTarget = std::fabs (m.integratedLufs - target) <= req.toleranceLu;
            const bool feasible = (viol == 0);

            best.offer (g, c, m, feasible, std::fabs (m.integratedLufs - target),
                        worstExcess (m, req), viol);
            if (m.integratedLufs <= target) { loG = g; loI = m.integratedLufs; haveLo = true; }
            else                            { hiG = g; hiI = m.integratedLufs; haveHi = true; }

            if (onTarget && feasible)
            {
                sol.status = MasteringSolveStatus::Solved;
                sol.preLimiterGainDb = g; sol.ceilingDbTp = c; sol.measured = m;
                return sol;
            }

            if (pass + 1 >= req.maxPasses) break;

            // --- choose the next (g, c) ------------------------------------------------------------
            const double shift = target - m.integratedLufs;
            double nextG = g, nextC = c;

            // THE EXACT STEP, and it is available more often than it looks. While the limiter is IDLE the
            // chain from the gain node on is a plain multiply, so `I(g)` has slope exactly 1 and the whole
            // correction is one addition — no slope estimate, no secant, no model. "Idle" is not assumed:
            // it is read off this render's own gain-reduction trace, and the range over which it stays
            // idle is read off the limiter's reconstructed peak, which the tap reports. Past that point
            // the limiter engages and the slope drops (measured on three real mixes: 0.99 at -14 LUFS,
            // 0.68 / 0.60 / 0.49 at -10), so the exact step is claimed only inside its own domain.
            const bool limiterIdle = m.limiter.valid && m.limiter.maxDb <= 0.0;
            const double headroomToEngage = c - m.limiterMaxReconstructedPeakDb;   // dB of gain left before
                                                                                   // the limiter starts working
            if (limiterIdle && std::isfinite (headroomToEngage))
            {
                anchorG = g + headroomToEngage;
                anchorI = m.integratedLufs + headroomToEngage;   // exact: the idle chain is a multiply
                haveAnchor = true;
            }
            if (limiterIdle && shift <= headroomToEngage)
            {
                nextG = g + shift;                                // exact
            }
            else if (shift < 0.0 && m.truePeakDbTp <= aim)
            {
                // GETTING QUIETER IS EXACT TOO, by a different route. `c` is a pure output scale, so
                // moving g and c DOWN together moves I and TP by exactly that amount while leaving the
                // limiting shape — and therefore PLR, LRA and every gain-reduction statistic — untouched.
                // It works downward and not upward for one reason and it is not symmetry: the ceiling is
                // the caller's delivery bound, and raising it above `maxTruePeakDbTp` would leave the
                // render with no guard at all. It is never raised, only lowered.
                nextG = g + shift;
                nextC = std::min (pmax, c + shift);
            }
            else
            {
                // The target is above what this drive delivers and the limiter is working, so the SHAPE
                // has to change and a slope is needed.
                //
                // THE LOCAL SECANT, and it has to be local. Taking it from the two points that BRACKET
                // the target uses whichever pair is widest, and the far one usually sits in the idle
                // region where the slope is exactly 1 — which over-estimates it and under-steps by the
                // same factor every time. Measured on this suite before the change: a target of -10 LUFS
                // walked -10.74, -10.33, -10.15 and ran out of budget, each step correcting only 55 % of
                // its own error, because the slope was read as 0.94 where it was 0.55. The two most
                // RECENT renders are the pair that describes where the search actually is.
                double s = 1.0;
                double pg = 0.0, pi = 0.0; bool pOk = false;
                if (havePrev) { pg = prevG; pi = prevI; pOk = true; }
                // Prefer the anchor when the previous render is still in the idle region: pairing an
                // active render with an idle one reads the slope as ~1 and under-steps by that factor.
                if (haveAnchor && (! pOk || pg < anchorG)) { pg = anchorG; pi = anchorI; pOk = true; }
                if (pOk && std::fabs (g - pg) > 1.0e-9)
                {
                    const double sl = (m.integratedLufs - pi) / (g - pg);
                    if (std::isfinite (sl) && sl > 0.02 && sl <= 1.2) s = sl;
                }
                else if (haveLo && haveHi && std::fabs (hiG - loG) > 1.0e-9)
                    s = (hiI - loI) / (hiG - loG);
                s = std::clamp (s, 0.05, 1.0);
                // A PLAIN SECANT STEP, and the reason there is no over-relaxation factor on it is worth
                // recording: `I(g)` is CONCAVE once the limiter engages (the slope falls from 1 to about
                // 0.2 as the drive rises — measured on three real mixes: 0.99 at -14 LUFS against
                // 0.68 / 0.60 / 0.49 at -10), so every extrapolated step lands SHORT and a deliberate
                // over-step to force a bracket looks obviously right. Measured, it is not: a 1.25x reach
                // moved a -12 LUFS solve from three renders to four, because the overshoot then had to be
                // walked back. The curvature costs a pass either way; paying it as an honest under-step
                // keeps the arithmetic explicable and leaves no tuned constant to go stale.
                nextG = g + shift / s;
                // ... and the ceiling TRACKS the aim, 1:1 and in both directions, because that is exactly
                // what `c` does to the delivered peak. It used to fall only, which ratchets: one
                // exploratory over-drive pushed it 0.5 dB below what the delivery needed and it never
                // came back, so the limiter kept working harder than the promise asked for ever after.
                // The cap at `pmax` is what keeps the promise; the floor is the search's own business.
                nextC = std::min (pmax, c + (aim - m.truePeakDbTp));
            }
            if (! std::isfinite (nextG) || ! std::isfinite (nextC)) break;
            nextC = std::min (nextC, pmax);      // the caller's ceiling is a bound, never a starting point
            // A step that cannot move the actuator cannot produce a new measurement: stop rather than
            // spend the remaining budget re-rendering the same point.
            if (std::fabs (nextG - g) < 1.0e-6 && std::fabs (nextC - c) < 1.0e-6) break;
            prevG = g; prevI = m.integratedLufs; havePrev = true;
            g = nextG; c = nextC;
        }

        // --- no candidate met the target -----------------------------------------------------------
        // THE GATE STEP IS ITS OWN ANSWER. Two renders straddling the target, both outside tolerance,
        // and an actuator gap too small to hold anything between them, is not a failure to converge and
        // not a constraint: it is the gated measure being discontinuous there. Saying so is the only
        // honest verdict, and it carries both sides so a caller can choose.
        if (haveLo && haveHi && std::fabs (hiG - loG) <= 1.0e-3
            && std::fabs (loI - target) > req.toleranceLu && std::fabs (hiI - target) > req.toleranceLu)
        {
            sol.status = MasteringSolveStatus::GateStep;
            sol.achievedBelowLufs = loI; sol.achievedAboveLufs = hiI;
            sol.gainBelowDb = loG;       sol.gainAboveDb = hiG;
        }
        else if (best.nearestViolated != 0)
        {
            sol.status  = MasteringSolveStatus::TargetUnreachable;
            sol.binding = bindingOf (best.nearestViolated);
            sol.alsoViolated = best.nearestViolated;
        }
        else
        {
            sol.status = MasteringSolveStatus::PassLimit;
        }

        if (best.have)
        {
            sol.preLimiterGainDb = best.g; sol.ceilingDbTp = best.c; sol.measured = best.m;
            // The DELIVERED render must be the one described. The last render written into `out` is the
            // last one attempted, which is not necessarily the best — re-render the reported one rather
            // than hand back a file the report does not describe.
            if (! best.isLast)
            {
                params.preLimiterGainDb    = best.g;
                params.limiter.ceilingDbTp = best.c;
                chain.setParams (params);
                MasterMeasurement again;
                if (! renderPass (chain, renderer, params, in, out, numChannels, frames, req, again))
                    { sol.status = MasteringSolveStatus::RenderFailed; return sol; }
                ++sol.passes;
                sol.measured = again;
            }
        }
        return sol;
    }

private:
    // Two candidates, and keeping them apart is what makes the verdict mean something.
    //   * `best`     — what gets DELIVERED. A feasible render always beats an infeasible one, whatever
    //                  their errors: handing back the render that broke the limit because it was 0.02 LU
    //                  closer is exactly the "crush it anyway" behaviour this class exists to refuse.
    //   * `nearest`  — the render closest to the target REGARDLESS of feasibility. It is the answer the
    //                  search would have delivered if nothing constrained it, so its violations are the
    //                  ones that stopped it. ORing violations across every render tried instead would
    //                  name constraints broken by an exploratory step the search had already left.
    struct Best
    {
        double g = 0.0, c = 0.0, err = std::numeric_limits<double>::infinity();
        double excess = std::numeric_limits<double>::infinity();
        MasterMeasurement m {};
        bool have = false, feasible = false, isLast = false;

        double nearestErr = std::numeric_limits<double>::infinity();
        std::uint32_t nearestViolated = 0;
        bool nearestHave = false;

        void offer (double gg, double cc, const MasterMeasurement& mm, bool feas, double e,
                    double exc, std::uint32_t viol) noexcept
        {
            if (! nearestHave || e < nearestErr) { nearestErr = e; nearestViolated = viol; nearestHave = true; }
            // WHEN NOTHING IS FEASIBLE, THE ANSWER IS THE GENTLEST RENDER, NOT THE CLOSEST ONE. Picking
            // the closest to a target that has already been declared unreachable delivers the most
            // crushed render there is AND a refusal to go with it — which is the "push it through
            // anyway" behaviour with a warning label. The tie-break among infeasible candidates is
            // therefore the WORST constraint excess, in the constraint's own units, and only then the
            // distance to the target.
            const bool better = ! have
                              || (feas && ! feasible)
                              || (feas == feasible && (feas ? (e < err)
                                                            : (exc < excess || (exc == excess && e < err))));
            isLast = better;
            if (! better) return;
            g = gg; c = cc; m = mm; err = e; excess = exc; have = true; feasible = feas;
        }
    };

    // How far past its limit the worst violated constraint is, in that constraint's own units (dB or
    // LU). Zero when nothing is violated. It is a RANKING, not a physical quantity — the units are only
    // comparable because every one of them is a logarithmic ratio and the caller set them all.
    double worstExcess (const MasterMeasurement& m, const LoudnessRequest& req) const noexcept
    {
        double e = 0.0;
        if (m.truePeakDbTp > req.maxTruePeakDbTp) e = std::fmax (e, m.truePeakDbTp - req.maxTruePeakDbTp);
        if (std::isfinite (req.limiterGr.limitDb) && m.limiter.valid)
        {
            const double v = (req.limiterGr.statistic == GrStatistic::Mean) ? m.limiter.meanDb
                           : (req.limiterGr.statistic == GrStatistic::P95)  ? m.limiter.p95Db : m.limiter.maxDb;
            if (v > req.limiterGr.limitDb) e = std::fmax (e, v - req.limiterGr.limitDb);
        }
        if (std::isfinite (req.minPlrDb) && m.plrDb < req.minPlrDb) e = std::fmax (e, req.minPlrDb - m.plrDb);
        if (std::isfinite (req.maxLraLossLu) && inputLraValid_ && m.lraValid)
        {
            const double loss = inputLraLu_ - m.loudnessRangeLu;
            if (loss > req.maxLraLossLu) e = std::fmax (e, loss - req.maxLraLossLu);
        }
        return e;
    }

    static MasteringConstraint bindingOf (std::uint32_t mask) noexcept
    {
        // ORDER IS AN ANSWER, not a tie-break, and it is ordered by WHAT A SEARCH CAN TRADE AWAY.
        // Reaching a louder target means more drive, and more drive REDUCES the true-peak excess (that
        // is what the limiter is for) while it INCREASES the limiter's gain reduction and REDUCES the
        // peak-to-loudness ratio and the loudness range. So when several are violated at once, the ones
        // the search cannot trade away are the reason it stopped, and a true-peak excess — always
        // fixable by bringing the ceiling down — is named last. Getting this backwards tells a user to
        // raise a ceiling that was never the problem.
        const MasteringConstraint order[] = { MasteringConstraint::LimiterGainReduction,
                                              MasteringConstraint::PeakToLoudness,
                                              MasteringConstraint::LoudnessRange,
                                              MasteringConstraint::GainRange,
                                              MasteringConstraint::TruePeakCeiling };
        for (MasteringConstraint k : order) if ((mask & constraintBit (k)) != 0) return k;
        return MasteringConstraint::None;
    }

    static bool violates (const GainReductionStats& s, const GainReductionLimit& lim) noexcept
    {
        if (! std::isfinite (lim.limitDb)) return false;
        if (! s.valid) return false;                 // an unanswerable statistic is not a violation
        const double v = (lim.statistic == GrStatistic::Mean) ? s.meanDb
                       : (lim.statistic == GrStatistic::P95)  ? s.p95Db : s.maxDb;
        return v > lim.limitDb;
    }

    std::uint32_t violatedMask (const MasterMeasurement& m, const LoudnessRequest& req) const noexcept
    {
        std::uint32_t v = 0;
        if (m.truePeakDbTp > req.maxTruePeakDbTp) v |= constraintBit (MasteringConstraint::TruePeakCeiling);
        if (violates (m.limiter, req.limiterGr))  v |= constraintBit (MasteringConstraint::LimiterGainReduction);
        if (std::isfinite (req.minPlrDb) && m.plrDb < req.minPlrDb)
            v |= constraintBit (MasteringConstraint::PeakToLoudness);
        // LRA is a DELTA against the input's, and it is only asked when both ends are measurements.
        if (std::isfinite (req.maxLraLossLu) && inputLraValid_ && m.lraValid
            && (inputLraLu_ - m.loudnessRangeLu) > req.maxLraLossLu)
            v |= constraintBit (MasteringConstraint::LoudnessRange);
        return v;
    }

    // One render, and every measurement of it. Everything here is measured on exactly the delivered
    // frames — the render is `frames` long by contract, and the drain is inside it.
    bool renderPass (MasteringChain& chain, OfflineRenderer& renderer, const MasteringChainParams& params,
                     const float* const* in, float* const* out, int nch, int frames,
                     const LoudnessRequest& req, MasterMeasurement& m)
    {
        (void) params;
        compHist_.reset();
        limHist_.reset();
        compActive_ = 0; limActive_ = 0; compFrames_ = 0; limFrames_ = 0;
        maxReconLin_ = 0.0f;
        if (! renderTapped (chain, renderer, in, out, nch, frames, req)) return false;

        m.latencySamples = chain.latencySamples();
        m.compressor = summarise (compHist_, compActive_, compFrames_);
        m.limiter    = summarise (limHist_,  limActive_,  limFrames_);
        m.limiterMaxReconstructedPeakDb = core::gainToDb ((double) maxReconLin_);
        return measure (out, nch, frames, m);
    }

    // The render, through `OfflineRenderer`'s OWN loop with a tap sink — not a second copy of the
    // alignment arithmetic. The sink below is the only thing this class adds to a plain render.
    bool renderTapped (MasteringChain& chain, OfflineRenderer& renderer,
                       const float* const* in, float* const* out, int nch, int frames,
                       const LoudnessRequest& req)
    {
        const int F = chain.tapOversampleFactor();
        // Each stage's own window into the tap stream — stated by MasteringChainTaps, read from the
        // chain, never guessed. Everything outside it is the chain's priming or its drain, both of
        // which are fed zeros and therefore produce no gain reduction and no peak; counting them would
        // dilute `mean` and the active fraction by exactly the ratio a short programme cannot afford.
        const MasteringChainResolved r = chain.resolved();
        const long long compTo  = frames;
        const long long limFrom = (long long) r.compressorLookahead + (long long) r.clipperLatency;
        const long long limTo   = limFrom + frames;
        const double activityDb = req.activityThresholdDb;

        MasteringChainTaps taps;
        taps.compressorGrDb = compTap_.data();
        taps.frameCapacity  = frameCap_;
        taps.limiterGrDb    = limTap_.data();
        taps.limiterPeakLin = limPeak_.data();
        taps.osCapacity     = osCap_;

        auto sink = [&] (const MasteringChainTaps& t, long long tapPos) noexcept
        {
            for (int j = 0; j < t.framesWritten; ++j)
            {
                const long long s = tapPos + j;
                if (s >= 0 && s < compTo)
                {
                    const double a = std::fabs ((double) compTap_[(std::size_t) j]);
                    compHist_.add (a);
                    if (a > activityDb) ++compActive_;
                    ++compFrames_;
                }
                if (s >= limFrom && s < limTo)
                    for (int k = 0; k < F; ++k)
                    {
                        const std::size_t idx = (std::size_t) j * (std::size_t) F + (std::size_t) k;
                        const double a = std::fabs ((double) limTap_[idx]);
                        limHist_.add (a);
                        if (a > activityDb) ++limActive_;
                        ++limFrames_;
                        const float pk = limPeak_[idx];
                        if (pk > maxReconLin_) maxReconLin_ = pk;
                    }
            }
        };
        return renderer.render (chain, in, out, nch, frames, taps, sink);
    }

    static GainReductionStats summarise (const dynamics::offline::QuantileHistogram& h,
                                         std::uint64_t active, std::uint64_t frames) noexcept
    {
        GainReductionStats s;
        s.frames = frames;
        s.nonFinite = h.nonFiniteCount();
        s.aboveRange = h.aboveRange();
        if (frames == 0 || h.count() == 0) return s;
        double p95 = 0.0;
        if (! h.quantile (0.95, p95)) return s;      // out of range: NOT reported as a plausible number
        s.meanDb = h.mean();
        s.p95Db  = p95;
        s.maxDb  = h.maxValue();
        s.activeFraction = (double) active / (double) frames;
        s.valid = (s.nonFinite == 0);
        return s;
    }

    // THE TWO RULES THE MEASURING RIG OWES, both measured rather than argued:
    //  * the loudness meter is sized for THIS programme and its `droppedBlocks()` is checked. A meter
    //    prepared for less silently answers about a prefix — measured, a 5 s meter given 20 s read
    //    -17.14 LUFS against the truth of -18.17, with the only outward sign a counter nobody reads.
    //  * the true-peak meter is DRAINED, and the drain is not given to the loudness meter. A programme
    //    ending on a peak under-reads without it: `[... 0, 1, 1]` reads +0.000000 dBTP undrained and
    //    +1.750350 once 8 zeros have gone through, and shipping the first number is precisely the
    //    defect P1 measured in the chain this replaces (17 of 36 rows above their own ceiling, reported
    //    as success).
    bool measure (float* const* out, int nch, int frames, MasterMeasurement& m)
    {
        analysis::LoudnessMeter  lm;
        analysis::TruePeakMeter  tm;
        const double seconds = (double) frames / fs_ + 1.0;
        if (! lm.prepare (fs_, nch, seconds)) return false;
        if (! tm.prepare (fs_, frames > 0 ? frames : 1, nch)) return false;
        const float* p[core::kMaxChannels] {};
        for (int c = 0; c < nch; ++c) p[c] = out[c];
        if (! lm.process (p, nch, frames)) return false;
        if (! tm.process (p, nch, frames)) return false;

        // 64 zeros: the meter's FIR holds 12 base-rate samples, and 8 were measured to be enough to
        // deliver the whole answer. Eight times that costs nothing and leaves no argument.
        drain_.assign ((std::size_t) nch * 64u, 0.0f);
        const float* z[core::kMaxChannels] {};
        for (int c = 0; c < nch; ++c) z[c] = drain_.data() + (std::size_t) c * 64u;
        if (! tm.process (z, nch, 64)) return false;

        m.gatingBlocks     = lm.gatingBlockCount();
        m.droppedBlocks    = lm.droppedBlocks();
        m.nonFiniteSubHops = lm.nonFiniteSubHops();
        m.integratedLufs   = lm.integratedLufs();
        m.loudnessRangeLu  = lm.loudnessRangeLu();
        m.truePeakDbTp     = tm.truePeakDb();
        m.samplePeakDb     = tm.samplePeakDb();
        m.plrDb            = m.truePeakDbTp - m.integratedLufs;
        // -120.0 EXACTLY is `LoudnessMeter`'s sentinel for "no gating block passed the absolute gate",
        // not a loudness: `integrated()` returns that literal from three different early exits. Treating
        // it as a measurement is how a solver ends up steering on digital silence — and the consequence
        // is not a wrong number but a DISABLED LIMITER, because a ceiling derived from a -200 dBTP peak
        // reading clamps to +60 dBTP and the true peak stops being bounded at all.
        m.loudnessValid    = (m.gatingBlocks > 0 && m.droppedBlocks == 0 && m.nonFiniteSubHops == 0
                              && m.integratedLufs > -120.0);
        // LRA needs short-term samples, one a second: a programme too short for them reports 0.0 LU,
        // which is also what "no dynamic range at all" reports. Saying which one it is is the only way
        // an LRA constraint can mean anything.
        m.lraValid         = m.loudnessValid && ((double) frames / fs_ >= 3.0);
        return true;
    }

public:
    // The input's loudness range, for the LRA constraint — which is a DELTA and therefore needs both
    // ends. Measured once, on the caller's input, before any render. Optional: an unset input LRA
    // simply switches that constraint off rather than making one up.
    [[nodiscard]] bool measureInputLra (const float* const* in, int nch, int frames)
    {
        inputLraValid_ = false;
        if (! prepared_ || in == nullptr || nch < 1 || nch > nch_ || frames <= 0) return false;
        analysis::LoudnessMeter lm;
        if (! lm.prepare (fs_, nch, (double) frames / fs_ + 1.0)) return false;
        if (! lm.process (in, nch, frames)) return false;
        if (lm.gatingBlockCount() <= 0 || lm.droppedBlocks() != 0) return false;
        if ((double) frames / fs_ < 3.0) return false;
        inputLraLu_ = lm.loudnessRangeLu();
        inputLraValid_ = true;
        return true;
    }

    double inputLoudnessRangeLu() const noexcept { return inputLraLu_; }
    bool   hasInputLoudnessRange() const noexcept { return inputLraValid_; }

private:
    double fs_ = 48000.0;
    int    nch_ = 0, frameCap_ = 0, osCap_ = 0;
    bool   prepared_ = false;

    std::vector<float> compTap_, limTap_, limPeak_, drain_;
    dynamics::offline::QuantileHistogram compHist_, limHist_;
    std::uint64_t compActive_ = 0, limActive_ = 0, compFrames_ = 0, limFrames_ = 0;
    float  maxReconLin_ = 0.0f;

    double inputLraLu_ = 0.0;
    bool   inputLraValid_ = false;
};

} // namespace felitronics::mastering

// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/analysis/LoudnessMeter.h>
#include <felitronics/analysis/ReferenceTruePeakMeter.h>
#include <felitronics/core/Config.h>
#include <felitronics/core/Math.h>
#include <felitronics/dynamics/offline/Quantile.h>
#include <felitronics/mastering/MasteringChain.h>
#include <felitronics/mastering/OfflineRenderer.h>
#include <felitronics/mastering/Planes.h>
#include <felitronics/mastering/Progress.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <utility>
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
// ceiling); the release recursion `grDb = min(rawRed, grDb*relCoef)` reads only `rawRed`, as do both envelopes of
// the limiter's `dualRelease`, so the whole gain trace is a function of d alone, and the output is
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
// WHAT THAT IS USED FOR HERE, precisely — because the identity is a REASON, not the algorithm, and an
// earlier draft of this comment described an algorithm this class does not implement. It is used for
// exactly two things:
//   * the loudness loop and the ceiling loop CANNOT chase each other, because they are not two loops.
//     Lowering the ceiling and raising the gain by the same amount leaves the limiting shape untouched,
//     so the ceiling correction is a trim rather than a competitor to the loudness correction;
//   * getting QUIETER is therefore an EXACT single step (see the step rule in solve()), with no slope
//     and no model, whenever the ceiling has room to come down with it.
// The search itself is a secant on MEASURED integrated loudness, not a closed-form inversion of `J`:
// the identity holds before the dither and while the gated block set does not move, and neither of
// those is guaranteed (the dither is after the limiter and does not scale; the absolute gate is not
// scale-invariant — see below). Every candidate is rendered and measured; the identity only ever
// chooses where to look next.
//
// `PLR(d)` being monotone non-increasing in d is measured rather than assumed — three real mixes,
// 16.53 -> 6.95, 13.17 -> 7.29 and 15.87 -> 7.34 over an 18 dB sweep, no reversal — and it is why the
// least drive that satisfies the ceiling is also the most transient the programme can keep at that
// loudness.
//
// THE DERATE IS NOT A CONSTANT HERE, and that is the point. `TruePeakLimiter` bounds its own F*fs grid
// exactly and overshoots the reconstructed peak between grid points — up to +1.22 dB in the worst
// non-degenerate case at 4x, and measured IN SITU on the corpus (with the pre-P62 meter) at only +0.0005 .. +0.1028 dB, because
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
// WHICH TRUE-PEAK METER, AND WHY IT IS THE EXPENSIVE ONE
// ====================================================================================
// `maxTruePeakDbTp` is a promise about the DELIVERED file, and a delivered file is certified by
// `analysis::ReferenceTruePeakMeter` — what `fcore_measure` and the browser's `fc_probe` report. So every
// render is read by that class and nothing else: the ceiling is aimed, feasibility is judged and
// `measured.truePeakDbTp` is reported on the same arithmetic the certificate runs, and the reported number
// IS the certificate — bit for bit, on the delivered samples. This class used to read with
// `analysis::TruePeakMeter`, the spec's short filter, and judged a render feasible by a number the
// certificate did not agree with: that meter stops interpolating at 2x / 1x on the high delivery rates, and
// on bright transients under-reads even at 4x. How far the two disagree is pinned by
// felitronics_truepeak_instrument_gap_tests. A wider aim margin was the alternative and was rejected: the
// disagreement depends on the material, so any margin is a guess about the worst programme and a tax on
// every other one.
// The price is the reference's filter on every pass, which is a fraction of the render the pass already
// is — and cheaper than any scheme that aims with one meter and verifies with the other, since that scheme
// needs an extra RENDER whenever the two disagree.
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
// OFFLINE, AND IT ALLOCATES — said plainly because the first draft of this line said it did not, and
// it does: every pass builds and prepares an `analysis::LoudnessMeter` sized for the programme, whose
// gating-block store is a function of the programme's length and therefore cannot live in `prepare()`
// where the length is unknown. That is fine here — this class is called from a worker, never from an
// audio callback — but it is not the RT-safe promise the chain underneath makes, and a caller must not
// read one as the other. What it does NOT do is allocate per BLOCK: the tap buffers and the render
// scratch are sized once. It drives `OfflineRenderer`, which drives `MasteringChain`, so every measurement it
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
    TargetBetweenAchievable,// the two sides of the smallest gain interval the search can still express
                            // BRACKET the target and both miss the tolerance. `achievedBelowLufs` and
                            // `achievedAboveLufs` carry them. It is honest unreachability and NOT a
                            // constraint violation: nothing was broken, the target simply is not an
                            // achievable value. Named for what is DETECTED rather than for a cause,
                            // because it has two — the gated measure stepping (a block crossing the
                            // absolute gate moves the SET being averaged, measured at 0.879829 LU on a
                            // fixture that straddles it), and a tolerance finer than the actuator's own
                            // resolution. The first is the interesting one and the second is the one a
                            // test can construct on demand
    PassLimit,              // ran out of `maxPasses` while still converging. Not a verdict about the
                            // material: it is a verdict about the budget, and it says so
    MeasurementInvalid,     // the meter could not answer (no gating block, dropped blocks, non-finite)
    RenderFailed,           // the chain or the renderer refused a call
    NotPrepared,
    InvalidRequest,
    Cancelled
};

// The pass budget's ceiling, named here because `LoudnessSolution` carries an array of that length and
// a struct cannot reach into the class that uses it.
struct TargetLoudnessSolverLimits { static constexpr int kMaxPasses = 32; };

// Which statistic a gain-reduction limit binds. Part of the limit's TYPE, never a hidden convention:
// "max 6 dB of limiting" and "p95 under 6 dB" are different products, and a field named for one while
// enforcing the other is the shape of defect this repository keeps closing.
enum class GrStatistic { Mean, P95, Max };

// `limitDb` OFF is `+infinity` and nothing else. `isfinite` looked like the right disabling test and is
// not: it is true of BOTH infinities and of NaN, so a caller expressing an unsatisfiable limit as
// `-infinity` silently switched the constraint OFF and got `Solved`. Measured: `minPlrDb = +infinity`
// — a peak-to-loudness ratio that must exceed infinity — came back Solved with nothing bound.
struct GainReductionLimit
{
    double      limitDb  = std::numeric_limits<double>::infinity();   // +infinity = no limit
    GrStatistic statistic = GrStatistic::Max;

    bool off()      const noexcept { return core::exactlyEqual (limitDb, std::numeric_limits<double>::infinity()); }
    bool malformed() const noexcept { return std::isnan (limitDb); }
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

// THE SAME TAPS, RESOLVED IN TIME (P59b). `GainReductionStats` says how much a stage worked over the whole programme;
// this says WHERE: the programme's frames cut into `buckets` uniform stretches, and per stretch the largest and the
// mean |GR|. A maximum per bucket and not a sample every N frames — a point sample steps over the peak it is meant to
// show. The same tap values the statistics see, in the same window (each stage's own tap offset), in the same units:
// frames for the compressor, frames × `tapOversampleFactor` sub-samples for the limiter.
//
// BUCKET k covers the programme frames [floor(k·F/B), floor((k+1)·F/B)) for F frames and B = min(grTraceBuckets, F)
// buckets — integer arithmetic, and never an empty bucket: a programme shorter than grTraceBuckets frames gets one
// bucket per frame, not trailing zeros that would read as "the stage did not work". Frame p of the programme is tap
// sample p + the stage's tap offset, so a bucket is in the time of the INPUT the gain was decided for, exactly as
// the statistics are; on the delivered-rate path (`DeliveredMastering`) that input is the converted programme, so
// the frames are delivered-rate frames.
struct GainReductionTraceBucket
{
    double        maxDb     = 0.0;  // largest finite |GR| in the bucket; 0 when it saw no finite sample
    double        meanDb    = 0.0;  // sequential sum of its finite |GR| / their count; 0 when it saw none
    std::uint64_t samples   = 0;    // tap samples that landed in it, finite or not
    std::uint64_t nonFinite = 0;    // ... of which were NaN or infinite: excluded from max and mean, and COUNTED,
                                    // because a max that skips a NaN silently reads 0 — "the stage was idle"
};

struct GainReductionTrace
{
    static constexpr int kDefaultBuckets = 1000;     // LoudnessRequest::grTraceBuckets' default
    static constexpr int kMaxBuckets     = 65536;    // the largest request admitted

    int           buckets   = 0;    // min(grTraceBuckets, programme frames) of the last render ATTEMPTED; 0 when none was
    // TRUE ONLY WHEN THE TRACE IS A MEASUREMENT: a render ran to its end for this solution, the stage's window saw
    // at least one sample, and none of them was non-finite. False on a refusal before any render (buckets 0), on a
    // render the renderer abandoned (RenderFailed from the render itself), and on a poisoned tap — in which case the
    // numbers are best effort and `nonFinite` says where. Like `GainReductionStats::valid`, not a nicety.
    bool          valid     = false;
    std::uint64_t samples   = 0, nonFinite = 0;
    std::vector<GainReductionTraceBucket> bucket;   // `buckets` entries

    // min(requested, frames); 0 when either is not positive.
    static int bucketsFor (int requested, int frames) noexcept
    {
        return frames > 0 && requested > 0 ? std::min (requested, frames) : 0;
    }
    // The bytes of that many buckets.
    static std::uint64_t bytesFor (int requested, int frames) noexcept
    {
        return (std::uint64_t) bucketsFor (requested, frames) * (std::uint64_t) sizeof (GainReductionTraceBucket);
    }
};

// THE TRACE'S THREE STEPS, as one small object the solver's tap sink drives — public so the one path no audio can
// reach through the solver (a non-finite tap: the chain sanitises its input) can be tested directly.
//   * construction — before a render: bucketsFor(requested, frames) buckets, `requested` clamped to [1, kMaxBuckets],
//     every bucket zero, `valid` false, the storage reused when it already holds that many;
//   * `add(frame, a)` — per tap sample, in stream order: `a` is |GR| in dB for programme frame `frame`; it counts the
//     sample, and either its non-finite count or its running max and SUM (the mean is divided out once, at the end);
//   * `finish()` — after a render that ran to its end: the means, the totals, and `valid`. A render that did not run
//     to its end never calls it, so its trace stays `valid == false` with whatever it had counted.
// The bucket of a frame is found by a cursor over the boundaries floor(k·F/B), advanced as the frames arrive — one
// division per bucket, not per sample. Frames outside [0, F) are ignored.
class GainReductionTraceBuilder
{
public:
    GainReductionTraceBuilder (GainReductionTrace& trace, int programmeFrames,
                               int requestedBuckets = GainReductionTrace::kDefaultBuckets)
        : t (trace), frames (programmeFrames > 0 ? (std::uint64_t) programmeFrames : 0u)
    {
        t.buckets   = GainReductionTrace::bucketsFor (std::clamp (requestedBuckets, 1, GainReductionTrace::kMaxBuckets),
                                                      programmeFrames);
        t.bucket.assign ((std::size_t) t.buckets, GainReductionTraceBucket {});
        t.valid     = false;
        t.samples   = 0;
        t.nonFinite = 0;
        end = t.buckets > 0 ? bucketEnd (0) : 0u;
    }

    void add (std::uint64_t frame, double a) noexcept
    {
        if (t.buckets <= 0 || frame >= frames) return;
        if (frame < start) { k = 0; start = 0; end = bucketEnd (0); }       // the stream never goes back; defensive
        while (frame >= end && k + 1 < (std::uint64_t) t.buckets) { ++k; start = end; end = bucketEnd (k); }
        GainReductionTraceBucket& b = t.bucket[(std::size_t) k];
        ++b.samples;
        if (! std::isfinite (a)) { ++b.nonFinite; return; }
        if (a > b.maxDb) b.maxDb = a;
        b.meanDb += a;                                                      // a running SUM until finish()
    }

    void finish() noexcept
    {
        std::uint64_t samples = 0, nonFinite = 0;
        for (int i = 0; i < t.buckets; ++i)
        {
            GainReductionTraceBucket& b = t.bucket[(std::size_t) i];
            const std::uint64_t finite = b.samples - b.nonFinite;
            b.meanDb = finite > 0 ? b.meanDb / (double) finite : 0.0;
            samples   += b.samples;
            nonFinite += b.nonFinite;
        }
        t.samples   = samples;
        t.nonFinite = nonFinite;
        t.valid     = samples > 0 && nonFinite == 0;
    }

private:
    // bucket k's end is bucket k+1's start: floor((k+1)·F/B), so a cursor that carries it over never recomputes a start
    std::uint64_t bucketEnd (std::uint64_t i) const noexcept { return (i + 1) * frames / (std::uint64_t) t.buckets; }

    GainReductionTrace& t;
    std::uint64_t frames;
    std::uint64_t k = 0, start = 0, end = 0;
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

// THE REQUEST HAS NO DEFAULT TARGET, and that is an architecture decision rather than an oversight.
// "-14 LUFS, -1 dBTP" is a delivery policy — a platform's, a label's, a taste — and the core is
// product-neutral by rule (docs/CORE-OVERVIEW.md: targets, curves and preset tables stay in the
// product). A core that shipped those numbers as defaults would be choosing the policy for every
// caller who forgot to, and forgetting is silent. So both are NaN and `solve()` refuses until the
// caller states them. `toleranceLu` DOES have a default, because it is a property of the measurement
// rather than of the product.
//
// "-1 dBTP" WITHOUT THE NAME OF AN INSTRUMENT IS HALF A PROMISE — one quantity, two definitions, the class of
// defect this core spent a sprint removing. The ceiling here is held as `analysis::ReferenceTruePeakMeter`
// reads it: the instrument `fcore_measure` certifies a file with, the one this class aims with, and so the
// one whose reading `measured.truePeakDbTp` IS. Like every BS.1770-class meter, that reference reads under
// the band-limited peak of the signal — by up to 0.33 dB on a full-band click, pinned in
// felitronics_truepeak_instrument_gap_tests — so another vendor's meter may read a delivered file higher than
// the promise. That is the method's property and a decision, not a gap to close here: moving the certifying
// instrument would move every certificate already issued.
struct LoudnessRequest
{
    double targetLufs      = std::numeric_limits<double>::quiet_NaN();   // REQUIRED
    double toleranceLu     =   0.1;
    double maxTruePeakDbTp = std::numeric_limits<double>::quiet_NaN();   // REQUIRED. DELIVERED, measured —
                                                                         // not the limiter's setting, which
                                                                         // this class derives.
    // How far BELOW the promise to aim the delivered peak. Not decoration and not taste: the ceiling
    // loop drives the delivered peak toward its aim, and an aim of exactly `maxTruePeakDbTp` converges
    // to the boundary FROM ABOVE and never crosses it — measured, a solve stalled at -1.0000 dBTP with
    // the ceiling constraint still flagged on three consecutive renders. The aim has only to exceed the
    // change in the limiter's own between-grid overshoot across one correction step, which is a few
    // thousandths of a dB; 0.05 is two decades of slack and is still an order below the 0.5 dB a
    // mastering engineer would notice.
    // THAT ARGUMENT HOLDS ONLY BECAUSE THE AIM AND THE PROMISE ARE READ BY ONE INSTRUMENT (P62). Before, the
    // solver read every render with `analysis::TruePeakMeter` while the delivered file was certified by the
    // reference, and the two disagree by more than this margin on bright material and at the high delivery
    // rates where the cheap meter stops interpolating — so a render the solver called feasible was delivered
    // above its promise. No margin fixes that, because the disagreement is the material's: see measure().
    double truePeakAimDb   =   0.05;

    // Constraints. A target that needs one of these broken is REFUSED with the name, not forced through.
    GainReductionLimit limiterGr    {};                                              // off by default
    GainReductionLimit compressorGr {};                                              // UPSTREAM — see above
    // NB with `MasteringChainParams::compressorMix < 1` this limits the COMPRESSED path's gain reduction,
    // not what reaches the output: at mix 0 a 26 dB reduction nobody hears is still refused against a
    // 3 dB limit. That is what the tap measures (MasteringChain.h), and a caller combining a GR budget
    // with parallel compression owns the translation.
    double minPlrDb     = -std::numeric_limits<double>::infinity();                  // delivered TP - I
    // A DELTA, not an absolute floor, and the difference is the whole calibration: the two proven
    // input->accepted-master pairs in this project's corpus move LRA by -0.40 and -0.30 LU, while the
    // ffmpeg chain being replaced collapses it 17.4 -> 2.9. An absolute floor cannot tell those apart —
    // the catalogue's own inputs run from 3.1 to 14.2 LU — and a delta can.
    double maxLraLossLu = std::numeric_limits<double>::infinity();
    // THE OTHER END OF THAT DELTA, and it lives in the REQUEST rather than in the solver because it is a
    // fact about THIS programme. Held as solver state it outlived the programme it was measured on:
    // solve A, then solve B without re-measuring, and B was judged against A's range — measured, 5.90
    // against 5.80 was enough to turn a healthy render into an `UpstreamViolation`. NaN means "not
    // supplied", which switches the range constraint off rather than inventing a number.
    // `TargetLoudnessSolver::measureInputLoudnessRange()` computes it; the caller passes it back in.
    double inputLoudnessRangeLu = std::numeric_limits<double>::quiet_NaN();

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

    // Buckets per gain-reduction trace of the solution: min(this, programme frames). 1..GainReductionTrace::kMaxBuckets,
    // else InvalidRequest.
    int    grTraceBuckets = GainReductionTrace::kDefaultBuckets;
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

    // RENDERS SPENT, all of them. `maxPasses` bounds the SEARCH; delivering the chosen candidate can
    // cost one more when the search did not end on it, so `passes` can be `maxPasses + 1` — and saying
    // so here is cheaper than a caller discovering it from a progress bar.
    int    passes = 0;
    double activityThresholdDb = 0.1;           // echoed, because a fraction without its threshold is not a number

    // TargetBetweenAchievable only: the two achievable values the target fell between.
    double achievedBelowLufs = 0.0, achievedAboveLufs = 0.0;
    double gainBelowDb = 0.0, gainAboveDb = 0.0;

    // Every render, in order, the delivery re-render of the winner included: `passes == logCount` while the log has room.
    // A re-render is the LAST record, repeating an earlier candidate's gain and ceiling.
    SolvePassRecord log[TargetLoudnessSolverLimits::kMaxPasses] {};
    int logCount = 0;

    // WHERE each stage worked, over the audio handed back in `out` (P59b) — see GainReductionTrace. Written by every
    // render and reset at the start of each, so it always describes the LAST render, and the last render is the one
    // in `out`: the search re-renders the winner when it did not end on it. It lives here, with the solution, and
    // not in the solver, because a solution outlives the next solve.
    GainReductionTrace compressorTrace {};
    GainReductionTrace limiterTrace {};
};

class TargetLoudnessSolver
{
public:
    static constexpr double kMaxGainDb = 60.0;      // MasteringChain::kMaxGainDb — the search's actuator range
    static constexpr int    kMaxPasses = TargetLoudnessSolverLimits::kMaxPasses;
    // The lowest rate the search measures at — the core's floor (P51). Below twice the K-weighting shelf the meters
    // this class builds are aliased and, in most of that range, unstable (a 0 dBFS sine read +3043 LUFS at 3300 Hz),
    // and below 8000 Hz what arrives is a
    // rate in kilohertz or a broken header: on origin/main a search handed 88.2 reported Solved at -14 LUFS.
    static constexpr double kMinSampleRate = core::kMinSampleRate;
    // And the highest — the chain's (P104). The search measures what a MasteringChain renders, at the chain's rate
    // (admits() checks that), so a rate the chain refuses has nothing to measure; before P104 prepare() took any
    // finite rate over the floor, 1e300 included, while solveBytes() answered 0 there.
    static constexpr double kMaxSampleRate = MasteringChain::kMaxSampleRate;

    // WHERE THIS CLASS STOPS REPORTING A dB AND STARTS REPORTING A SENTINEL. peakDb() below is the only
    // user; the constant is public so a test can pin WHERE it is, not merely that silence reads -200.
    // Digital silence sits below any plausible gate, so a test using silence alone cannot tell this value
    // from one ten times larger — an adversarial round raised it to 1e-9f and the whole suite stayed green
    // while the report and the certificate disagreed by 6 dB at a peak between the two. The float spelling
    // widened to double is deliberate and predates P80: it keeps the boundary exactly where it was.
    static constexpr double kPeakDbGate    = (double) 1.0e-10f;
    static constexpr double kPeakDbSilence = -200.0;

    // `maxFrames` and `maxChannels` size the tap buffers; `binDb` is the resolution every gain-reduction
    // quantile is reported to. The tap buffers are the whole allocation and they are per RENDERER BLOCK,
    // not per programme — the traces are consumed as they arrive, so a five-minute track costs the same
    // as a five-second one.
    TargetLoudnessSolver() noexcept { for (double& w : weights_) w = 1.0; }

    [[nodiscard]] bool prepare (double sampleRate, int maxChannels, int rendererBlock,
                                int internalBlock, int oversampleFactor, double binDb = 0.01)
    {
        prepared_ = false;
        if (! rateAdmitted (sampleRate)) return false;
        if (maxChannels < 1 || maxChannels > core::kMaxChannels) return false;
        if (rendererBlock < 1 || internalBlock < 1 || oversampleFactor < 1) return false;
        if (! (binDb > 0.0) || ! std::isfinite (binDb)) return false;
        int frameCap = 0, osCap = 0;
        if (! tapLayoutFor (rendererBlock, internalBlock, oversampleFactor, frameCap, osCap)) return false;
        // The histograms' refusal BEFORE the first write (law 11(b)): a prepare() refused on its bin width allocates
        // nothing, which is what prepareBytes() says of it. (The diverse-testing round: the tap buffers used to be
        // assigned first, and kept.)
        std::size_t bins = 0;
        if (! dynamics::offline::QuantileHistogram::binsFor (0.0, kGrRangeDb, binDb, bins)) return false;

        fs_ = sampleRate;
        nch_ = maxChannels;
        frameCap_ = frameCap;
        osCap_    = osCap;
        compTap_.assign ((std::size_t) frameCap_, 0.0f);
        limTap_.assign  ((std::size_t) osCap_, 0.0f);
        limPeak_.assign ((std::size_t) osCap_, 0.0f);
        // The gain-reduction range: `GainComputer` caps its own range at 400 dB, and the limiter's is
        // bounded by its ceiling clamp. 400 covers both, and anything past it is COUNTED rather than
        // folded into the top bin, so a quantile that lands there answers `false` instead of lying.
        if (! compHist_.prepare (0.0, kGrRangeDb, binDb)) return false;
        if (! limHist_.prepare  (0.0, kGrRangeDb, binDb)) return false;
        prepared_ = true;
        return true;
    }

    bool isPrepared() const noexcept { return prepared_; }
    // The rate the meters are built at — 0 until prepared. A caller that hands this solver a programme at another
    // rate (a converted one, say) checks it here; see `DeliveredMastering`.
    double sampleRate() const noexcept { return prepared_ ? fs_ : 0.0; }

    // EVERY VERDICT `solve()` REACHES BEFORE ITS FIRST PASS that does not look at the audio pointers: `NotPrepared`,
    // or `InvalidRequest` for a chain, a width, a length, a request, a rate or a tap geometry it will not search.
    // One definition, read by `solve()` itself — so a caller that has to spend memory BEFORE the search (a
    // converted programme, `DeliveredMastering`) can ask first and spend nothing on a call that would be refused.
    // True, and `why` untouched, when the call would reach a pass.
    [[nodiscard]] bool admits (const MasteringChain& chain, const OfflineRenderer& renderer, int numChannels, int frames,
                               const LoudnessRequest& req, MasteringSolveStatus& why) const noexcept
    {
        if (! prepared_) { why = MasteringSolveStatus::NotPrepared; return false; }
        why = MasteringSolveStatus::InvalidRequest;
        if (! chain.isPrepared() || numChannels != chain.numChannels() || numChannels > nch_ || frames <= 0) return false;
        if (! std::isfinite (req.targetLufs) || ! std::isfinite (req.maxTruePeakDbTp)
            || ! std::isfinite (req.toleranceLu) || req.toleranceLu < 0.0
            || ! std::isfinite (req.activityThresholdDb) || req.activityThresholdDb < 0.0
            || req.maxPasses < 1 || req.maxPasses > kMaxPasses
            || req.limiterGr.malformed() || req.compressorGr.malformed()
            || std::isnan (req.minPlrDb) || std::isnan (req.maxLraLossLu)
            || ! std::isfinite (req.truePeakAimDb) || req.truePeakAimDb < 0.0
            || ! traceBucketsAdmitted (req.grTraceBuckets)) return false;
        // The tap buffers were sized for a geometry; a chain that does not match them would be measured
        // through a refused call, which is a silent zero rather than a statistic.
        // THE RATE IS CHECKED, not assumed shared. The solver builds its own meters from `fs_`, and a
        // caller that passed 44100 here and 48000 to the chain would get a 48 kHz render measured on a
        // 44.1 kHz grid — every number plausible, every number wrong.
        if (! (std::fabs (chain.sampleRate() - fs_) < 1.0e-9)) return false;
        if (chain.internalBlock() + renderer.blockSize() > frameCap_
            || (long long) (chain.internalBlock() + renderer.blockSize()) * chain.tapOversampleFactor() > (long long) osCap_)
            return false;
        why = MasteringSolveStatus::Solved;     // not read by a caller on `true`; restored so `why` is untouched
        return true;
    }

    //==========================================================================================================
    // THE BUDGETS — what a call will ask the heap for, computed by the very functions the call sizes itself with
    // (tapLayoutFor, QuantileHistogram::binsFor, meterSamples, LoudnessMeter::storageFor, ReferenceTruePeakMeter::storageFor),
    // so a budget cannot drift from its allocation — exact for a FRESH object: one already prepared keeps whatever
    // storage still fits and asks nothing for it. REQUESTED bytes: allocator headers, alignment and fragmentation
    // are the caller's margin, and none of this is a promise that a heap can serve it. Static on purpose — a caller
    // budgets before it prepares anything, and the rate is an argument, not state.
    static constexpr double kGrRangeDb   = 400.0;    // the gain-reduction histograms' span — see prepare()

    // prepare(): the tap buffers and the two histograms. 0 where prepare() refuses the same arguments.
    static std::uint64_t prepareBytes (int rendererBlock, int internalBlock, int oversampleFactor, double binDb = 0.01) noexcept
    {
        if (rendererBlock < 1 || internalBlock < 1 || oversampleFactor < 1) return 0;
        int frameCap = 0, osCap = 0;
        if (! tapLayoutFor (rendererBlock, internalBlock, oversampleFactor, frameCap, osCap)) return 0;
        const std::uint64_t hist = dynamics::offline::QuantileHistogram::storageBytes (0.0, kGrRangeDb, binDb);
        if (hist == 0) return 0;
        return (std::uint64_t) sizeof (float) * ((std::uint64_t) frameCap + 2u * (std::uint64_t) osCap) + 2u * hist;
    }

    // solve(): its PEAK. Every pass builds a loudness meter and the reference true-peak meter and frees them at the
    // pass's end, so the peak is ONE pass. (The drain used to be a buffer of zeros the first solve allocated and later
    // ones reused; the reference meter drains from its own fixed array, so there is nothing left over.) 0 for a length
    // or a channel count solve() refuses before any pass.
    // Plus the solution's two traces of `grTraceBuckets` buckets, which the first render allocates; 0 for a bucket
    // count solve() refuses.
    static std::uint64_t solveBytes (double sampleRate, int numChannels, int frames, int grTraceBuckets) noexcept
    {
        if (frames <= 0 || numChannels < 1 || numChannels > core::kMaxChannels) return 0u;
        if (! traceBucketsAdmitted (grTraceBuckets)) return 0u;
        const std::uint64_t meter = meterBytes (sampleRate, frames);
        if (meter == 0) return 0u;       // the meter refuses its capacity: measure() stops before anything is allocated
        return meter + analysis::ReferenceTruePeakMeter::storageFor (sampleRate, frames, numChannels).bytes()
             + 2u * GainReductionTrace::bytesFor (grTraceBuckets, frames);
    }

    // The request's bucket count, as admits() judges it.
    static bool traceBucketsAdmitted (int grTraceBuckets) noexcept
    {
        return grTraceBuckets >= 1 && grTraceBuckets <= GainReductionTrace::kMaxBuckets;
    }

    // measureInputLoudnessRange(): one loudness meter — and NOTHING for a programme too short to have a range, which it
    // refuses before building one. (The code-review round: this budget used to promise a meter for a 1 s call that
    // allocates none.)
    static std::uint64_t measureRangeBytes (double sampleRate, int frames) noexcept
    {
        // The rate first: rangeMeasurable() divides by it. (The answer was 0 either way — meterBytes() refuses the
        // same rates — but a division by a zero rate is not a question this budget should have to ask.)
        return rateAdmitted (sampleRate) && rangeMeasurable (frames, sampleRate) ? meterBytes (sampleRate, frames) : 0u;
    }

    // Render `frames` of `in` into `out` at a gain and ceiling chosen to meet `req`. `params` is the
    // caller's whole parameter set; the solver overrides exactly `preLimiterGainDb` and
    // `limiter.ceilingDbTp` and leaves every other field alone.
    //
    // `params` is taken by value rather than read back from the chain ON PURPOSE: `MasteringChain::params()`
    // reports the last APPLIED set, which lags a `setParams()` by up to one internal quantum, and it
    // stores the caller's UNCLAMPED request while the chain applies a clamped one. A search that read
    // its own actuator through either of those would be measuring a number it did not apply.
    LoudnessSolution solve (MasteringChain& chain, OfflineRenderer& renderer,
                            const MasteringChainParams& params,
                            const float* const* in, float* const* out,
                            int numChannels, int frames, const LoudnessRequest& req)
    {
        return solve (chain, renderer, params, in, out, numChannels, frames, req, ProgressCallback {});
    }

    LoudnessSolution solve (MasteringChain& chain, OfflineRenderer& renderer,
                            MasteringChainParams params,
                            const float* const* in, float* const* out,
                            int numChannels, int frames, const LoudnessRequest& req,
                            const ProgressCallback& progress)
    {
        LoudnessSolution sol;
        sol.activityThresholdDb = req.activityThresholdDb;
        if (! admits (chain, renderer, numChannels, frames, req, sol.status)) return sol;
        // `in == out` IS REFUSED HERE, even though `OfflineRenderer` supports it. One render in place is
        // well defined; a SEARCH is not, because every pass after the first would read the previous
        // pass's master as its input. Measured: a 1 kHz tone solved to a reported -22.996 LUFS, and the
        // gain it returned applied to the untouched source gives -29.000 — the answer misses its own
        // programme by 6.0 LU, and every number in the report describes a programme the caller does not
        // have. Refused rather than copied: the copy is the caller's memory to spend, and only the
        // caller knows whether it can.
        // AND NOT ONLY CHANNEL AGAINST ITSELF: no output plane may touch ANY input plane, nor another output plane
        // (`planesUsable`, Planes.h). This used to test `in[c] == out[c]`, which let `out[0] = in[1]` through — the
        // render writes channel 0's master where channel 1 is read next pass, and the call answered an ordinary
        // verdict, at a plausible gain, over a master that is not the programme's — and `out[0] = out[1]`, where the
        // second channel's render overwrites the first and the search meters one channel twice (the suite's
        // witnesses, testCrossChannelAliasingIsRefused). A NULL PLANE is refused by the same predicate; it used to
        // reach the renderer and dereference it.
        if (! planesUsable (in, out, numChannels, frames, frames))
            { sol.status = MasteringSolveStatus::InvalidRequest; return sol; }

        const double target = req.targetLufs;
        const double pmax   = req.maxTruePeakDbTp;
        // THE CEILING IS THE SOLVER'S, THE PROMISE IS THE CALLER'S. `params.limiter.ceilingDbTp` is a
        // STARTING POINT; `maxTruePeakDbTp` is the bound, and the only thing this class guarantees about
        // the delivered file. The ceiling then TRACKS the aim in both directions, capped at the promise —
        // it can rise above what the caller set (a caller ceiling of -40 against a promise of -1 is
        // relieved by raising it, not by refusing) and it can fall below the promise by whatever the
        // material's between-grid overshoot turns out to be. What it can never do is exceed the promise:
        // shipping above a stated ceiling with the interface reporting success is the exact defect P1
        // measured 17 times in 36 in the chain this replaces.
        const double c0     = std::min (pmax, std::isfinite (params.limiter.ceilingDbTp)
                                                  ? params.limiter.ceilingDbTp : pmax);
        double g = std::isfinite (req.initialGainDb) ? req.initialGainDb : params.preLimiterGainDb;
        if (! std::isfinite (g)) g = 0.0;
        double c = c0;

        Best best;                          // the best FEASIBLE render seen, and the best of any kind
        double loD = 0.0, hiD = 0.0, loJ = 0.0, hiJ = 0.0;   // the two sides, in (drive, shape)
        double loG = 0.0, hiG = 0.0, loC = 0.0, hiC = 0.0, loI = 0.0, hiI = 0.0;   // ...and in the caller's units, for the report
        bool haveLo = false, haveHi = false;
        double prevD = 0.0, prevJ = 0.0;    // the PREVIOUS render, in (drive, shape) coordinates
        bool   havePrev = false;
        // Which end of the +-60 dB gain node the search wanted to pass, if any. Judged at the end.
        int    pinnedDir = 0;
        double lastG = 0.0, lastI = 0.0, lastC = 0.0; bool haveLastRender = false;
        // Set when the search stops because it CANNOT MOVE — the step it wants is below the resolution
        // the actuator can express — as opposed to running out of budget while still making progress.
        // The two are different answers and used to be the same one.
        bool   stoppedOnResolution = false;
        bool   bootstrapped = false;        // one attempt to bring an unmeasurable programme into range
        // THE IDLE ANCHOR, and it is a measurement rather than a model. While the limiter does not
        // engage, the chain from the gain node on is a plain multiply, so `J(d) = J1 + (d - d1)` holds
        // EXACTLY up to the DRIVE at which the limiter starts working — and that drive is `c` minus the
        // reconstructed peak the limiter itself reports. (In gain it was `I(g) = I1 + (g - g1)`, which
        // is the same statement only while `c` stands still; the rewrite moved it to drive and this
        // comment was left behind in the old coordinates.) So a single idle render hands the search a
        // second exact point sitting ON the boundary of the active region, which is the anchor a local
        // secant wants: the first active render then pairs with it instead of with a point far away in
        // the linear region.
        //
        // MEASURED, AND THE PLACE IT EARNS ITS KEEP MOVED once the other defects were fixed. It is not
        // the loud target (no change) and no longer the saturated one (a -2 LUFS target is now one
        // render FASTER without it). It is the WARM START just above the answer: four consecutive
        // starts from +8.30 to +8.45 dB toward a -10.5 LUFS target take three renders with the anchor
        // and four without, every time. Over the 102-cell battery, 23 cells move and the totals are
        // 120 renders with against 124 without -- a net win, and a small one. Two earlier numbers
        // written here were true of code that has since changed; this one is dated to the fixes above.
        double anchorD = 0.0, anchorJ = 0.0;
        bool   haveAnchor = false;
        DriveBound bound;
        const double aim = pmax - (std::isfinite (req.truePeakAimDb) && req.truePeakAimDb > 0.0
                                       ? req.truePeakAimDb : 0.0);

        ProgressClock clock (progress);
        const long long passUnits = 2LL * (long long) frames + (long long) chain.latencySamples();

        for (int pass = 0; pass < req.maxPasses; ++pass)
        {
            if (! clock.begin (ProgressStage::SearchPass, pass + 1, req.maxPasses, passUnits, frames))
                return cancelled (sol);
            g = std::clamp (g, -kMaxGainDb, kMaxGainDb);
            c = std::clamp (c, -kMaxGainDb, kMaxGainDb);
            params.preLimiterGainDb      = g;
            params.limiter.ceilingDbTp   = c;
            chain.setParams (params);

            MasterMeasurement m;
            if (! renderPass (chain, renderer, params, in, out, numChannels, frames, req, m, sol, clock))
            {
                sol.status = clock.stopped() ? MasteringSolveStatus::Cancelled : MasteringSolveStatus::RenderFailed;
                sol.passes = pass + 1;
                return sol;
            }
            ++sol.passes;

            // AN UNMEASURABLE FIRST RENDER IS NOT ALWAYS AN UNMEASURABLE PROGRAMME. A file quiet enough
            // that every gating block sits under the absolute gate has no integrated loudness at the
            // gain it was rendered at — and may have a perfectly ordinary one 55 dB up, which is inside
            // the actuator. Measured: a flat tone at -71.69 LUFS returned `MeasurementInvalid` from a
            // start of 0 dB and `Solved` in one render from a start of +55.7. Refusing on that is a
            // verdict about the starting gain, not about the material.
            //
            // The way out uses the measurement that DID work: the peak. Bring the sample peak to a
            // sensible distance under the promise and let the ordinary search take it from there. One
            // attempt only — if the peak cannot be read either, there really is nothing to measure.
            if (! m.loudnessValid && ! bootstrapped && m.samplePeakDb > -180.0 && pass + 1 < req.maxPasses)
            {
                bootstrapped = true;
                if (! keepRecord (sol, clock, g, c, m, 0u)) return cancelled (sol);
                // 12 dB under the promise: far enough below it that the limiter does not take over the
                // next measurement, high enough that an ordinary programme's blocks clear the -70 gate.
                g = std::clamp (g + ((pmax - 12.0) - m.samplePeakDb), -kMaxGainDb, kMaxGainDb);
                continue;
            }
            if (! m.loudnessValid)
            {
                sol.status = MasteringSolveStatus::MeasurementInvalid;
                sol.measured = m; sol.preLimiterGainDb = g; sol.ceilingDbTp = c;
                if (! keepRecord (sol, clock, g, c, m, 0u)) return cancelled (sol);
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
                if (! keepRecord (sol, clock, g, c, m, constraintBit (MasteringConstraint::CompressorGainReduction)))
                    return cancelled (sol);
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
            // ...and only when the solver's OWN knobs cannot relieve it. `c` is capped at the promise
            // but may still be raised toward it, and raising the ceiling is exactly what removes limiter
            // gain reduction. Measured: a caller ceiling of -40 dBTP against a promise of -1 made the
            // limiter pull 20 dB on the first render, and the guard called a violation "upstream" that
            // `g = -6.996, c = -1` reaches with no gain reduction at all.
            if (pass == 0 && target > m.integratedLufs && c >= pmax - 1.0e-9)
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
                    if (! keepRecord (sol, clock, g, c, m, viol)) return cancelled (sol);
                    return sol;
                }
            }

            if (! keepRecord (sol, clock, g, c, m, viol)) return cancelled (sol);
            const bool onTarget = std::fabs (m.integratedLufs - target) <= req.toleranceLu;
            const bool feasible = (viol == 0);
            bound.add (g, c, m, aim, pmax, req, viol);

            best.offer (g, c, m, feasible, std::fabs (m.integratedLufs - target),
                        worstExcess (m, req), viol);
            lastG = g; lastI = m.integratedLufs; lastC = c; haveLastRender = true;
            // THE BRACKET ONLY EVER TIGHTENS. Overwriting each side with the most RECENT render on it
            // is not a bracket: once the search converges from one side, the other side's record stays
            // where it was many dB ago, and the interval never closes — which makes the "the target
            // falls between two achievable values" verdict unreachable and turns it into a `PassLimit`.
            // `J(d)` is increasing, so the useful sides are the LARGEST DRIVE that undershoots and the
            // SMALLEST that overshoots -- drive, not gain, since the rewrite: the gain and the ceiling
            // are carried alongside only so the report can speak the caller's units.
            if (m.integratedLufs <= target)
            { if (! haveLo || (g - c) > loD) { loD = g - c; loJ = m.integratedLufs - c; loG = g; loC = c; loI = m.integratedLufs; haveLo = true; } }
            else
            { if (! haveHi || (g - c) < hiD) { hiD = g - c; hiJ = m.integratedLufs - c; hiG = g; hiC = c; hiI = m.integratedLufs; haveHi = true; } }

            if (onTarget && feasible)
            {
                sol.status = MasteringSolveStatus::Solved;
                sol.preLimiterGainDb = g; sol.ceilingDbTp = c; sol.measured = m;
                return sol;
            }

            if (pass + 1 >= req.maxPasses) break;

            // --- choose the next (g, c) ------------------------------------------------------------
            // THE SEARCH RUNS IN THE CHAIN'S OWN COORDINATES, `d = g - c` and `J = I - c`, AND IT USED
            // TO RUN IN `g` AND `I`. That was the same mistake the class's header exists to warn about,
            // made one level down: `I` is not a function of `g`, it is `c + J(g - c)`, and `c` MOVES
            // every pass while the ceiling tracks its aim. A secant taken across two renders whose
            // ceilings differ therefore measures nothing — it divides a change that is mostly `c` by a
            // change that is all `g`.
            //
            // Measured, and both of these are one defect: a 1 kHz tone at 0.005, started at +55 dB with
            // a target of -10 LUFS, gave two renders 9 dB apart whose ceilings differed by 0.05 dB; the
            // "slope" came out **0.0057** where the true one is near 1, the step went to the -60 dB
            // clamp, the render fell under the absolute gate, and a four-render budget ended at
            // **-12.993 LUFS instead of -10.000**. And a start of +55 dB with the caller's ceiling at
            // -40 spent twenty-two renders bisecting a bracket that never contained the answer, then
            // reported `Unreachable / GainRange` for a target that `g = -18.99, c = -1` delivers
            // exactly.
            //
            // In `(d, J)` both go away, because `J` is a function of `d` alone WHEREVER THE SCALE LAW
            // HOLDS. The exception is named rather than glossed: BS.1770's RELATIVE gate is
            // scale-invariant, but its ABSOLUTE gate at -70 LUFS is not, so two renders with the same
            // drive and different ceilings straddling that gate have different `J`. Measured on the
            // gate fixture, `(g, c) = (0, -3)` and `(2, -1)` — both `d = 3` — differ by 0.874 LU. The
            // consequence is bounded and local: it is a step whose secant is wrong near -70 LUFS, and
            // the physical bracket below is what keeps that from becoming a wild step. Nothing in this
            // file may claim the identity is unconditional; the header says the same at the top.
            //
            // The ceiling is chosen FIRST — it is the trim, and it is what makes `J` mean anything —
            // and the step is then taken on the shape.
            const double dNow = g - c;
            const double jNow = m.integratedLufs - c;
            // The ceiling tracks the aim, in both directions, capped at the promise.
            double nextC = std::min (pmax, c + (aim - m.truePeakDbTp));
            if (! std::isfinite (nextC)) nextC = c;
            double nextD = dNow;

            // "IDLE" HERE IS STRICTER THAN THE REQUEST'S `activityThresholdDb`, and deliberately so.
            // The two words look like they should agree -- `activeFraction` counts samples with
            // |GR| > `activityThresholdDb` -- but they gate different things. That one describes a
            // render for the CALLER; this one licenses the exact 1:1 step below, which is an identity
            // (`J(d) = J(d0) + (d - d0)`) and holds only while the chain is a plain multiply. A limiter
            // doing 0.09 dB is inactive by the report's standard and NOT a plain multiply.
            //
            // WHAT THIS DOES NOT CLAIM is that the loose form would be harmless. Two measurements of it
            // disagree, and the disagreement is between FIXTURES rather than between readings: on the
            // shipped test rig, a warm start above the target moves the delivered loudness by up to
            // 0.09 LU with the pass count unchanged; on a standalone probe over the same programme and
            // the same targets it is bit-identical. Both were run; neither is wrong. What settles the
            // choice is not inertness but the arithmetic: with a worst reduction of `h`, the exact step
            // is wrong by the loudness that reduction cost (call it e, in [0, h]) and the conservative
            // end of the bracket below is wrong by `h - e`. They are mirror images of one bound, and
            // which is smaller is a property of the material. The identity is the thing this code can
            // actually prove, so it is the thing it is allowed to assume.
            const bool limiterIdle = m.limiter.valid && m.limiter.maxDb <= 0.0;
            // dB of DRIVE left before the limiter starts working. `c - maxReconstructedPeak` is the same
            // number in either coordinate system: the limiter engages at `d = -R(p)`, and
            // `R(p) = maxReconstructedPeak - g`.
            const double headroomToEngage = c - m.limiterMaxReconstructedPeakDb;
            if (limiterIdle && std::isfinite (headroomToEngage))
            {
                anchorD = dNow + headroomToEngage;
                anchorJ = jNow + headroomToEngage;       // exact: while idle the chain is a multiply
                haveAnchor = true;
            }
            // A render the limiter will work on, after an idle one, is aimed `overshootAt` its drive under the aim.
            if (limiterIdle && std::isfinite (headroomToEngage) && (target - nextC) - jNow > headroomToEngage)
                nextC = std::min (nextC, aim - bound.overshootAt (dNow + (target - nextC) - jNow));
            const double jReq = target - nextC;          // the shape has to deliver this much
            const double dj   = jReq - jNow;

            if (limiterIdle && dj <= headroomToEngage)
            {
                nextD = dNow + dj;                       // exact, and it stays exact under a moving ceiling
            }
            else
            {
                // The local secant, in the coordinates where it is a slope. Prefer the anchor when the
                // previous render is still in the idle region: pairing an active render with an idle one
                // reads the slope as ~1 and under-steps by that factor.
                double s = 1.0;
                double pd = 0.0, pj = 0.0; bool pOk = false;
                if (havePrev) { pd = prevD; pj = prevJ; pOk = true; }
                if (haveAnchor && (! pOk || pd < anchorD)) { pd = anchorD; pj = anchorJ; pOk = true; }
                if (pOk && std::fabs (dNow - pd) > 1.0e-9)
                {
                    // NO SMALL-SLOPE REJECTION — the bound that stays is `sl <= 1.2`, which throws away
                    // a slope steeper than the scale law allows (`J` cannot outrun `d`), plus the
                    // non-positive and non-finite cases. A slope under 0.02 is not an
                    // implausible measurement, it is the TRUTH in the saturated region, where more drive
                    // buys almost no loudness; a floor of 0.02 there fell back to 1 and crept at a
                    // fiftieth of the step the measurement called for — measured, a -5.4 LUFS target ran
                    // out of ten renders 0.194 LU short and reported `PassLimit` where it now SOLVES.
                    const double sl = (jNow - pj) / (dNow - pd);
                    if (std::isfinite (sl) && sl > 0.0 && sl <= 1.2) s = sl;
                }
                else if (haveLo && haveHi && std::fabs (hiD - loD) > 1.0e-9)
                    s = (hiJ - loJ) / (hiD - loD);
                // THE FLOOR IS 0.01, NOT 0.05, AND THE DIFFERENCE IS MEASURED. In the saturated region
                // the honest secants run 0.019, 0.013, 0.011, 0.0089 — every one of them under 0.05, so
                // the old floor replaced the measurement with a step twenty times too small on exactly
                // the material where the search is already struggling. What makes a small slope safe to
                // believe is not the floor but the physical bracket below, which bounds the STEP; the
                // floor only has to stop a slope of zero from producing infinity.
                //
                // Measured over the four hardest cells of the corpus battery, floor 0.05 -> 0.01: every
                // one improves on BOTH counts, renders and accuracy — 8->7, 9->7, 10->7 renders, error
                // 0.075->0.054, 0.097->0.072, and a -5.3 LUFS target that was a `PassLimit` 0.122 LU
                // short now SOLVES 0.090 short. Nothing else in the battery's 102 cells moves.
                s = std::clamp (s, 0.01, 1.0);
                nextD = dNow + dj / s;
                // AND THE STEP IS BOUNDED BY PHYSICS — which is what the old `sl > 0.02` guard was
                // reaching for and getting wrong, because it bounded the SLOPE (a measurement) where
                // the thing that needed bounding was the STEP (a decision).
                //
                // Precisely, since the earlier wording claimed more than the code does: a floor on the
                // slope REMAINS — it is the `clamp` ten lines above, and it is the number that line
                // executes, not a second one written out here (this paragraph said `0.05` after that
                // line had already become `0.01`, in the same commit that lowered it) — and so
                // does an upper rejection (`sl <= 1.2` throws away a finite, positive 1.3). What the
                // rewrite removed is the 0.02 REJECTION threshold — a slope below it used to be
                // discarded in favour of the default 1.0, which is the fifty-fold overstep. The floor
                // that stayed is a bound on the step size; the bracket below is what makes it safe.
                //
                // Two facts bound it, and both are already measured every pass. Below the drive at which
                // the limiter engages the chain is a plain multiply, so `J` moves 1:1 there; above it the
                // limiter can take away anything at all. The engagement point is known on EVERY render,
                // active or not — `headroomToEngage = c - maxReconstructedPeak` is the distance to it —
                // so for a step DOWNWARD the answer is trapped between "the limiter gave back nothing on
                // the way down" and "it gave back all of it":
                //
                //     dNow + headroom + dj   <=   answer   <=   dNow + dj        (for dj < 0)
                //
                // and for a step UPWARD only the lower end survives, because the slope above engagement
                // may be anything down to zero. Without this, a slope measured in the saturated region —
                // genuinely near zero, and now correctly so — turned a -8.96 dB requirement into a 179 dB
                // step, landed on the -60 dB clamp, rendered digital silence and spent a pass climbing
                // back: a four-render budget ended at -12.993 LUFS against a target of -10.000. With it,
                // the same request solves in THREE -- the figure said four until the diff pass counted
                // it. Getting the direction of the bound wrong is its own
                // trap: clamping to the conservative end alone left a +55 dB warm start creeping 8 dB a
                // render and out of budget at -8.33.
                if (dj > 0.0)
                    nextD = std::fmax (nextD, dNow + dj);
                else if (dj < 0.0)
                {
                    const double a = dNow + dj, b = dNow + headroomToEngage + dj;
                    nextD = std::clamp (nextD, std::fmin (a, b), std::fmax (a, b));
                }
                // A BRACKET IS A GUARANTEE AND A SECANT IS NOT. The two sides are kept in `d` as well,
                // and which side each is on is re-decided against the CURRENT requirement — an entry
                // recorded under an older ceiling can otherwise sit on the wrong side of it.
                if (haveLo && haveHi && (loJ - jReq) * (hiJ - jReq) < 0.0)
                {
                    const double bLo = std::fmin (loD, hiD), bHi = std::fmax (loD, hiD);
                    if (! (nextD > bLo && nextD < bHi)) nextD = 0.5 * (bLo + bHi);
                }
            }
            // DriveBound::wants the step: `probe` replaces it, or `closed` ends the search.
            if (bound.wants (nextD))
            {
                bound.acted = true;
                if (best.have && best.feasible && bound.closed (best.m.integratedLufs, target, req.toleranceLu)) break;
                nextD = bound.probe (req, pass + 2 >= req.maxPasses);
                nextC = std::min (pmax, aim - bound.overshootAt (bound.cap.d));
            }

            double nextG = nextD + nextC;
            // `nextC` is finite and already `<= pmax` by construction above; only `nextG` can arrive
            // non-finite here, through `nextD`. The second half of this test and the `min (nextC, pmax)`
            // that used to follow it were both dead — kept only the live one.
            if (! std::isfinite (nextG)) break;
            const double clG = std::clamp (nextG, -kMaxGainDb, kMaxGainDb);
            const double clC = std::clamp (nextC, -kMaxGainDb, kMaxGainDb);
            // THE ACTUATOR'S LIMIT IS RECORDED IN BOTH DIRECTIONS AND JUDGED AT THE END, not the moment
            // an exploratory step touches it. Recording it immediately named `GainRange` for a target the
            // search simply had not looked at yet; testing only `shift > 0` lost the other bound
            // entirely, so a target 3 dB below what -60 dB delivers came back as a pass limit.
            const bool clamped = (std::fabs (clG - nextG) > 1.0e-9);
            // ALREADY AT THE CLAMP AND STILL ASKING TO GO PAST IT. That is the actuator's limit, and it
            // is judged on the GAIN alone: if the gain cannot rise and more loudness is wanted, the
            // ceiling cannot rescue it — lowering it makes the render quieter and raising it is capped
            // at the promise. Testing the pair (gain AND ceiling both still) instead let a ceiling that
            // was merely still tracking its aim hide the pin, and the verdict came back `PassLimit`.
            // The distinction from an exploratory touch is the `clG == g` half: a step that MOVES to the
            // clamp is a step the search has not evaluated yet, and naming it would be a verdict about a
            // target nobody has looked at — measured, that called a -25 LUFS target unreachable that
            // `g = -18.99` delivers exactly.
            // THE CEILING NEVER STANDS STILL AT 1e-6. The true peak jitters by a couple of parts in a
            // million from render to render, so `clC` keeps moving by ~2e-6 forever and this test never
            // fires: measured, six consecutive renders at g = +60, c = -1.11132 +- 2e-6 with identical
            // audio, five of eleven renders carrying no information at all. The resolution that matters
            // is the one `bracketClosed` already uses.
            if (std::fabs (clG - g) < 1.0e-6 && std::fabs (clC - c) < 1.0e-3)
            {
                if (! clamped) stoppedOnResolution = true;
                break;
            }
            nextG = clG; nextC = clC;
            prevD = g - c; prevJ = m.integratedLufs - c; havePrev = true;
            g = nextG; c = nextC;
        }

        // --- no candidate met the target -----------------------------------------------------------
        // A TARGET BETWEEN TWO ACHIEVABLE VALUES IS ITS OWN ANSWER. Two renders straddling the target,
        // both outside tolerance, and a gain gap too small to hold anything between them, is not a
        // failure to converge and not a constraint. Saying so is the only honest verdict, and it carries
        // both sides so a caller can choose which one to take.
        // ORDER: A BROKEN LIMIT OUTRANKS "NOT EXACTLY ACHIEVABLE", and it did not.
        //
        // The two verdicts answer different questions and the wrong one used to win. "The target lies
        // between two achievable values" is a statement about the SEARCH — it needs the candidate
        // NEAREST the target — while `best` is the candidate that gets DELIVERED, which when nothing is
        // feasible is deliberately the gentlest rather than the nearest. Testing `best.err` therefore
        // asked "how far is the render I am handing back from the target", got a large answer because
        // the search had walked away from a target it could not legally reach, and reported that
        // distance as a resolution limit — hiding the constraint violations entirely (`alsoViolated` is
        // only filled on the unreachable path).
        //
        // Measured, and this is a mutation that found a defect in the CLEAN code rather than a hole in
        // the suite: a request with `maxLraLossLu = 0.01` and `minPlrDb = 40` — both broken by every
        // candidate — came back `TargetBetweenAchievable`, binding `None`, mask `0x0`, delivering
        // -19.29 LUFS against a target of -24. The truth is that -24 IS reachable and costs the
        // loudness range and the peak-to-loudness ratio the caller forbade.
        //
        // So: constraints first, and the between-achievable test reads the NEAREST candidate's error.
        // THE ACTUATOR'S LIMIT IS A STATE, NOT AN EVENT. Recording it only when a step tried to walk
        // past the clamp made the verdict depend on the pass BUDGET rather than on the material: the
        // search renders at +60 dB, the budget ends, and the one extra iteration that would have
        // noticed the pin never runs. Measured on the suite's own generator, target -5.3 LUFS: with
        // `maxPasses = 10` the answer was `PassLimit` with an empty mask at g = +60.000, with 11 it was
        // `Unreachable / GainRange` -- same input, same delivered audio, two different verdicts.
        //
        // So it is decided here from what was actually RENDERED: the last render sat at the clamp, the
        // target is still outside tolerance, and closing the gap needs gain the range does not have.
        // Reading the last render rather than the next STEP also keeps the extrapolation out of it --
        // a step's direction can disagree with the measurement it was extrapolated from, and did.
        // JUDGED ON THE LAST RENDER, both halves. Reading `best.err` for "did we succeed" while
        // reading `lastI` for "which way is the target" mixes two different renders: `best` is the
        // candidate that gets DELIVERED, and when nothing is feasible that is deliberately the gentlest
        // rather than the nearest. A last render inside tolerance but breaking the ceiling would then
        // pick up `GainRange` on top of the true-peak violation that is the real answer.
        if (haveLastRender && best.have && std::fabs (target - lastI) > req.toleranceLu)
        {
            const double want = target - lastI;                 // > 0 asks for more gain
            // THE CEILING IS AN ACTUATOR TOO, and on the loud side it has not necessarily been tried.
            // `c` only ever tracks the true-peak aim, so a caller who starts it far below the promise
            // leaves real loudness on the table: raising `c` by one dB raises `I` by about one, and
            // until `c` reaches `pmax` the gain node is not the only thing that could close the gap.
            // Measured: a start of `g = +60, c = -40` with one render's budget was called
            // `Unreachable / GainRange` for a -25 LUFS target that `g = -18.99, c = -1.05` delivers
            // exactly — the search had 39 dB of untried ceiling and a verdict saying it had none.
            // On the QUIET side there is no such escape: `nextC` never lowers the ceiling to chase a
            // quiet target, so the gain node really is the only actuator and the pin stands.
            //
            // `want > 0.0` IS REDUNDANT ON THIS SIDE and kept only as the written intent: `nextC` is
            // `min(pmax, ...)`, so `c <= pmax` always and `ceilingLeft >= 0`, hence `want > ceilingLeft`
            // already implies `want > 0`. The mutation stand proved it — dropping the term changes no
            // answer on any input — and it is recorded here rather than removed because the quiet arm
            // below has no such implication and needs its own direction test to stay legible.
            const double ceilingLeft = pmax - lastC;
            if (lastG >=  kMaxGainDb - 1.0e-6 && want > 0.0 && want > ceilingLeft) pinnedDir = +1;
            if (lastG <= -kMaxGainDb + 1.0e-6 && want < 0.0) pinnedDir = -1;
        }

        // CLOSED ON BOTH KNOBS. The sides are keyed on drive now, but this test still compared gains
        // only, and two renders at the SAME gain with different ceilings have a gap of exactly zero:
        // `(60, -40)` and `(60, -1.05)` straddle a -25 LUFS target, pass the gap test, and report
        // `TargetBetweenAchievable` for a target that is simply reachable. The interval has to be small
        // in the coordinate the search actually moves in, which is both of them.
        const bool bracketClosed = haveLo && haveHi && std::fabs (hiG - loG) <= 1.0e-3
                                && std::fabs (hiC - loC) <= 1.0e-3
                                && std::fabs (loI - target) > req.toleranceLu
                                && std::fabs (hiI - target) > req.toleranceLu;
        const bool betweenAchievable = best.nearestViolated == 0
                                    && (bracketClosed
                                        || (stoppedOnResolution && best.nearestHave
                                            && best.nearestErr > req.toleranceLu));
        // Once the drive bound has acted: `TargetUnreachable`, `binding` from what `cap` broke.
        if (bound.acted)
        {
            std::uint32_t viol = best.nearestViolated | bound.capViol;
            if (pinnedDir != 0 && best.have && best.err > req.toleranceLu)
                viol |= constraintBit (MasteringConstraint::GainRange);
            sol.status  = MasteringSolveStatus::TargetUnreachable;
            sol.binding = bindingOf (bound.capViol);
            sol.alsoViolated = viol;
        }
        else if (betweenAchievable)
        {
            sol.status = MasteringSolveStatus::TargetBetweenAchievable;
            if (bracketClosed)
            {
                sol.achievedBelowLufs = loI; sol.achievedAboveLufs = hiI;
                sol.gainBelowDb       = loG; sol.gainAboveDb       = hiG;
            }
            else
            {
                // THE SEARCH COULD NOT MOVE FROM HERE, so "here" is the whole answer and both fields
                // carry it. Reporting the two sides of a WIDE bracket instead would be a different
                // claim — that the target lies between two gains a dB apart — and it was: the two
                // branches ran on different rows of the same suite (a bracket on the wasm tier, one
                // side on arm64 macOS), and only the second reported an interval the caller could act
                // on. The condition that stopped the search is the one that gets reported.
                sol.achievedBelowLufs = sol.achievedAboveLufs = best.m.integratedLufs;
                sol.gainBelowDb       = sol.gainAboveDb       = best.g;
            }
        }
        else if (best.nearestViolated != 0 || (pinnedDir != 0 && best.have && best.err > req.toleranceLu))
        {
            // THE ACTUATOR'S LIMIT IS JUDGED HERE, at the end, and in both directions. Recording it the
            // moment an exploratory step touched +-60 dB named `GainRange` for a target the search had
            // simply not looked at yet — measured, a start of +55 dB with the caller's ceiling at -40
            // reported `Unreachable / GainRange` for a target that `g = -18.99` delivers exactly. And
            // testing only the loud direction lost the quiet one: a target 3 dB below what -60 dB
            // delivers came back a pass limit with nothing named.
            std::uint32_t viol = best.nearestViolated;
            if (pinnedDir != 0 && best.have && best.err > req.toleranceLu)
                viol |= constraintBit (MasteringConstraint::GainRange);
            sol.status  = MasteringSolveStatus::TargetUnreachable;
            sol.binding = bindingOf (viol);
            sol.alsoViolated = viol;
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
            // `isLast` is set by the last `offer` that WON. A candidate re-offered on equal terms does
            // not win, so a search that ended on its own best point used to re-render it — a whole pass
            // for a buffer that already held the right audio.
            const bool alreadyDelivered = best.isLast
                                       || (std::fabs (g - best.g) < 1.0e-12 && std::fabs (c - best.c) < 1.0e-12);
            if (! alreadyDelivered)
            {
                if (! clock.begin (ProgressStage::FinalRender, sol.passes + 1, req.maxPasses + 1, passUnits, frames))
                    return cancelled (sol);
                params.preLimiterGainDb    = best.g;
                params.limiter.ceilingDbTp = best.c;
                chain.setParams (params);
                MasterMeasurement again;
                if (! renderPass (chain, renderer, params, in, out, numChannels, frames, req, again, sol, clock))
                {
                    if (clock.stopped()) { ++sol.passes; return cancelled (sol); }
                    sol.status = MasteringSolveStatus::RenderFailed;
                    return sol;
                }
                ++sol.passes;
                sol.measured = again;
                if (! keepRecord (sol, clock, best.g, best.c, again, violatedMask (again, req))) return cancelled (sol);
            }
        }
        return sol;
    }

private:
    static bool keepRecord (LoudnessSolution& sol, ProgressClock& clock, double g, double c,
                            const MasterMeasurement& m, std::uint32_t violated) noexcept
    {
        SolvePassRecord rec;
        rec.gainDb = g; rec.ceilingDb = c;
        rec.integratedLufs = m.integratedLufs; rec.truePeakDbTp = m.truePeakDbTp;
        rec.plrDb = m.plrDb; rec.limiterMaxGrDb = m.limiter.maxDb;
        rec.loudnessRangeLu = m.loudnessRangeLu; rec.violated = violated;
        if (sol.logCount < kMaxPasses) sol.log[sol.logCount++] = rec;
        return clock.finish (&rec);
    }

    // `sol` as Cancelled, moved into the return value.
    static LoudnessSolution cancelled (LoudnessSolution& sol) noexcept
    {
        sol.status = MasteringSolveStatus::Cancelled;
        return std::move (sol);
    }

    // Two candidates, and keeping them apart is what makes the verdict mean something.
    //   * `best`     — what gets DELIVERED. A feasible render always beats an infeasible one, whatever
    //                  their errors: handing back the render that broke the limit because it was 0.02 LU
    //                  closer is exactly the "crush it anyway" behaviour this class exists to refuse.
    //   * `nearest`  — the render closest to the target REGARDLESS of feasibility. It is the answer the
    //                  search would have delivered if nothing constrained it, so its violations are the
    //                  ones that stopped it. ORing violations across every render tried instead would
    //                  name constraints broken by an exploratory step the search had already left.
    // `exc == excess` COMPARES TWO DOUBLES FOR EQUALITY ON PURPOSE, and the obvious hardening is wrong.
    // Replacing it with a tolerance window (`fabs(exc - excess) <= 1e-12`) looks strictly safer and is
    // not: `minPlrDb` may legally be `+infinity` — the API rejects only NaN — and then every candidate's
    // excess is `+infinity`, `fabs(inf - inf)` is NaN, every comparison against it is false, and the
    // ranking freezes on whatever was offered first. Measured: target -20 LUFS came back at -6.014 with
    // the gain still at its starting 0 dB, fourteen LU out, where the equality delivers -20.000 exactly.
    // Equality is what makes the infinite case tie and fall through to the distance test.

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
                                                            : (exc < excess || (core::exactlyEqual (exc, excess) && e < err))));
            isLast = better;
            if (! better) return;
            g = gg; c = cc; m = mm; err = e; excess = exc; have = true; feasible = feas;
        }
    };

    // The limits that grow with drive, in the order of DriveBound's excesses.
    static constexpr std::uint32_t kDriveBound = constraintBit (MasteringConstraint::LimiterGainReduction)
                                               | constraintBit (MasteringConstraint::PeakToLoudness)
                                               | constraintBit (MasteringConstraint::LoudnessRange);
    static constexpr MasteringConstraint kDriveBoundConstraint[3] {
        MasteringConstraint::LimiterGainReduction, MasteringConstraint::PeakToLoudness, MasteringConstraint::LoudnessRange };

    // The overshoot `truePeakDbTp - c` assumed for a working render when none is measured at or under its drive.
    static constexpr double kFirstLimitingOvershootDb = 0.15;

    // THE DRIVE BOUND over the renders made, in drive `d = g - c`: `ok` is the largest drive that broke no `kDriveBound` limit,
    // `cap` the smallest that broke one, each with every limit's excess (> 0 where violatedMask() flags it, NaN when off or
    // unmeasured); `cap` also with what it broke and `capAimedI`, its loudness at the ceiling its true peak asks for. A bracket
    // is `ok` under `cap`. The bound only chooses renders; each is measured and judged like any other.
    struct DriveBound
    {
        struct End { double d = 0.0; double e[3] {}; bool have = false; };

        End ok, cap;
        std::uint32_t capViol = 0;
        double capAimedI = 0.0, engageD = 0.0;
        bool   okIdle = false, haveEngage = false;
        bool   acted = false;                        // it chose a render or ended the search
        int    lastSide = 0, streak = 0;             // the end last moved (-1 `ok`, +1 `cap`), and how often in a row
        double activeD[kMaxPasses] {}, activeOvershoot[kMaxPasses] {};
        int    nActive = 0;

        void add (double g, double c, const MasterMeasurement& m, double aim, double pmax, const LoudnessRequest& req,
                  std::uint32_t viol) noexcept
        {
            const double d = g - c;
            const bool idle = m.limiter.valid && m.limiter.maxDb <= 0.0;
            const double engage = g - m.limiterMaxReconstructedPeakDb;
            if (std::isfinite (engage)) { engageD = engage; haveEngage = true; }
            if (m.limiter.valid && m.limiter.maxDb > 0.0 && nActive < kMaxPasses)
            {
                activeD[nActive] = d;
                activeOvershoot[nActive] = std::fmax (0.0, m.truePeakDbTp - c);
                ++nActive;
            }
            const double nan = std::numeric_limits<double>::quiet_NaN();
            double e[3] = { nan, nan, nan };
            if (! req.limiterGr.off() && ! req.limiterGr.malformed() && m.limiter.valid)
                e[0] = ((req.limiterGr.statistic == GrStatistic::Mean) ? m.limiter.meanDb
                      : (req.limiterGr.statistic == GrStatistic::P95)  ? m.limiter.p95Db
                                                                        : m.limiter.maxDb) - req.limiterGr.limitDb;
            if (! core::exactlyEqual (req.minPlrDb, -std::numeric_limits<double>::infinity()))
                e[1] = req.minPlrDb - m.plrDb;
            if (! core::exactlyEqual (req.maxLraLossLu, std::numeric_limits<double>::infinity())
                && std::isfinite (req.inputLoudnessRangeLu) && m.lraValid)
                e[2] = (req.inputLoudnessRangeLu - m.loudnessRangeLu) - req.maxLraLossLu;

            const std::uint32_t broke = viol & kDriveBound;
            int moved = 0;
            if (broke == 0u)
            {
                if (! ok.have || d > ok.d) { set (ok, d, e); okIdle = idle; moved = -1; }
            }
            else if (! cap.have || d < cap.d)
            {
                set (cap, d, e); moved = +1;
                capViol   = broke;
                capAimedI = (m.integratedLufs - c) + std::min (pmax, c + (aim - m.truePeakDbTp));
            }
            if (moved == 0) return;
            streak   = (moved == lastSide) ? streak + 1 : 1;
            lastSide = moved;
        }

        // A bracket exists and `nextD` is not under `cap`.
        bool wants (double nextD) const noexcept
        {
            return ok.have && cap.have && ok.d < cap.d && nextD >= cap.d;
        }

        // `capAimedI` is under the target's tolerance and within `tol` of `feasibleI`.
        bool closed (double feasibleI, double target, double tol) const noexcept
        {
            return capAimedI < target - tol && capAimedI - feasibleI <= tol;
        }

        // Held `margin = min (toleranceLu, width / 2)` inside the bracket: the smallest regula-falsi root over what `cap` broke,
        // from `engageD` when `ok` is idle under it, the still end's excess halved `streak - 1` times; the midpoint when an
        // excess is not measured. `margin` lower for the last render allowed.
        double probe (const LoudnessRequest& req, bool lastRender) const noexcept
        {
            const double w = cap.d - ok.d, margin = std::min (req.toleranceLu, 0.5 * w);
            const double dl = (okIdle && haveEngage && engageD > ok.d && engageD < cap.d) ? engageD : ok.d;
            double next = cap.d;
            for (int k = 0; k < 3; ++k)
            {
                if ((capViol & constraintBit (kDriveBoundConstraint[k])) == 0u) continue;
                double el = ok.e[k], eh = cap.e[k];
                if (! (el <= 0.0 && eh > 0.0)) { next = ok.d + 0.5 * w; break; }
                if (streak >= 2 && lastSide < 0) eh = std::ldexp (eh, 1 - streak);
                if (streak >= 2 && lastSide > 0) el = std::ldexp (el, 1 - streak);
                next = std::fmin (next, dl - el * (cap.d - dl) / (eh - el));
            }
            if (lastRender) next -= margin;
            return std::clamp (next, ok.d + margin, cap.d - margin);
        }

        // The largest overshoot measured on a working render at or under drive `d`, else kFirstLimitingOvershootDb.
        double overshootAt (double d) const noexcept
        {
            double o = -1.0;
            for (int i = 0; i < nActive; ++i)
                if (activeD[i] <= d) o = std::fmax (o, activeOvershoot[i]);
            return o >= 0.0 ? o : kFirstLimitingOvershootDb;
        }

    private:
        static void set (End& s, double d, const double e[3]) noexcept
        {
            s.d = d;
            std::copy (e, e + 3, s.e);
            s.have = true;
        }
    };

    // How far past its limit the worst violated constraint is, in that constraint's own units (dB or
    // LU). Zero when nothing is violated. It is a RANKING, not a physical quantity — the units are only
    // comparable because every one of them is a logarithmic ratio and the caller set them all.
    double worstExcess (const MasterMeasurement& m, const LoudnessRequest& req) const noexcept
    {
        double e = 0.0;
        if (m.truePeakDbTp > req.maxTruePeakDbTp) e = std::fmax (e, m.truePeakDbTp - req.maxTruePeakDbTp);
        if (! req.limiterGr.off() && ! req.limiterGr.malformed() && m.limiter.valid)
        {
            const double v = (req.limiterGr.statistic == GrStatistic::Mean) ? m.limiter.meanDb
                           : (req.limiterGr.statistic == GrStatistic::P95)  ? m.limiter.p95Db : m.limiter.maxDb;
            if (v > req.limiterGr.limitDb) e = std::fmax (e, v - req.limiterGr.limitDb);
        }
        if (! core::exactlyEqual (req.minPlrDb, -std::numeric_limits<double>::infinity()) && m.plrDb < req.minPlrDb)
            e = std::fmax (e, req.minPlrDb - m.plrDb);
        if (! core::exactlyEqual (req.maxLraLossLu, std::numeric_limits<double>::infinity())
            && std::isfinite (req.inputLoudnessRangeLu) && m.lraValid)
        {
            const double loss = req.inputLoudnessRangeLu - m.loudnessRangeLu;
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
        if (lim.off()) return false;                 // +infinity, and ONLY +infinity, means "no limit"
        if (lim.malformed()) return false;           // NaN is refused up front; never silently permissive
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
        // OFF is `-infinity` for a FLOOR and `+infinity` for a CEILING — the sign is part of the
        // meaning, and testing `isfinite` throws it away in the direction that always says "satisfied".
        if (! core::exactlyEqual (req.minPlrDb, -std::numeric_limits<double>::infinity()) && m.plrDb < req.minPlrDb)
            v |= constraintBit (MasteringConstraint::PeakToLoudness);
        // LRA is a DELTA against the input's, and it is only asked when both ends are measurements.
        if (! core::exactlyEqual (req.maxLraLossLu, std::numeric_limits<double>::infinity())
            && std::isfinite (req.inputLoudnessRangeLu) && m.lraValid
            && (req.inputLoudnessRangeLu - m.loudnessRangeLu) > req.maxLraLossLu)
            v |= constraintBit (MasteringConstraint::LoudnessRange);
        return v;
    }

    // One render, and every measurement of it. Everything here is measured on exactly the delivered
    // frames — the render is `frames` long by contract, and the drain is inside it.
    bool renderPass (MasteringChain& chain, OfflineRenderer& renderer, const MasteringChainParams& params,
                     const float* const* in, float* const* out, int nch, int frames,
                     const LoudnessRequest& req, MasterMeasurement& m, LoudnessSolution& sol, ProgressClock& clock)
    {
        (void) params;
        compHist_.reset();
        limHist_.reset();
        compActive_ = 0; limActive_ = 0; compFrames_ = 0; limFrames_ = 0;
        maxReconLin_ = 0.0f;
        // The traces are reset with the histograms, at the same place and for the same reason: whatever they held
        // described a render that is about to be overwritten in `out`.
        GainReductionTraceBuilder compTrace (sol.compressorTrace, frames, req.grTraceBuckets);
        GainReductionTraceBuilder limTrace  (sol.limiterTrace, frames, req.grTraceBuckets);
        if (! renderTapped (chain, renderer, in, out, nch, frames, req, compTrace, limTrace, clock)) return false;
        compTrace.finish();
        limTrace.finish();

        m.latencySamples = chain.latencySamples();
        m.compressor = summarise (compHist_, compActive_, compFrames_);
        m.limiter    = summarise (limHist_,  limActive_,  limFrames_);
        // `gainToDbDet` for the same reason as peakDb() above, and this one is the stronger case of the two:
        // the field crosses the C ABI into the browser (fc_master_abi.h, fc_master.cpp) AND it is read back
        // as a DECISION — `headroomToEngage = ceiling - m.limiterMaxReconstructedPeakDb` a few hundred lines
        // up — so on the system spelling the solver could take a different branch on Apple than on the row
        // that rendered the same file. The LINEAR peak it converts is the limiter's own, and stays RT.
        m.limiterMaxReconstructedPeakDb = core::gainToDbDet ((double) maxReconLin_);
        return measure (out, nch, frames, m, clock);
    }

    // The render, through `OfflineRenderer`'s OWN loop with a tap sink — not a second copy of the
    // alignment arithmetic. The sink below is the only thing this class adds to a plain render.
    bool renderTapped (MasteringChain& chain, OfflineRenderer& renderer,
                       const float* const* in, float* const* out, int nch, int frames,
                       const LoudnessRequest& req, GainReductionTraceBuilder& compTrace, GainReductionTraceBuilder& limTrace,
                       ProgressClock& clock)
    {
        const int F = chain.tapOversampleFactor();
        // Each stage's own window into the tap stream — stated by MasteringChainTaps, read from the
        // chain, never guessed. Everything outside it is the chain's priming or its drain, and neither is
        // programme — which, and not that they are empty, is why they are left out. The drain in particular is
        // NOT free of gain reduction: a release still running when the programme ends runs on into it, and an
        // expanding or upward mode acts on the drain's silence itself. Counting outside the window would let the
        // chain's latency and the programme's last moments, not the programme, move the statistics — and a short
        // programme's the most.
        // EACH TAP'S WINDOW COMES FROM THE CHAIN, not from arithmetic here. The limiter's offset in
        // particular is not the sum of the stages in front of it — the trace is written where the gain
        // is decided, on the oversampled copy, so it lags by the UP leg of the limiter's own oversampler
        // as well. Deriving it here got it wrong twice (once omitting that term entirely, once counting
        // the whole round trip instead of half), and the cost of omitting it is the whole reaction to a
        // peak in the programme's last samples: measured, 10.9 dB of gain reduction reported as zero.
        const MasteringChainResolved r = chain.resolved();
        const long long compFrom = r.compressorTapOffset;
        const long long compTo   = compFrom + frames;
        const long long limFrom  = r.limiterTapOffset;
        const long long limTo    = limFrom + frames;
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
                if (s >= compFrom && s < compTo)
                {
                    const double a = std::fabs ((double) compTap_[(std::size_t) j]);
                    compHist_.add (a);
                    if (a > activityDb) ++compActive_;
                    ++compFrames_;
                    compTrace.add ((std::uint64_t) (s - compFrom), a);
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
                        limTrace.add ((std::uint64_t) (s - limFrom), a);
                    }
            }
        };
        return renderer.render (chain, in, out, nch, frames, taps, sink, &clock);
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
    //    +1.833993 drained (the reference, pinned in felitronics_reference_truepeak_tests; the cheap meter
    //    this class used before P62 read +1.750350 after 8 zeros), and shipping the first number is precisely the
    //    defect P1 measured in the chain this replaces — rows shipping ABOVE their own ceiling while
    //    the interface reports success. The COUNT deliberately does not live here: its owner is the
    //    baseline harness in another repository, it moves whenever that corpus does, and nothing in
    //    this tree can re-derive it. A number without a local owner rots and cannot be made not to.
    // THE RATE THIS CLASS MEASURES AT — the one test, read by prepare() and by every budget, so the two cannot disagree
    // about a rate. NaN fails both comparisons, -inf the first and +inf the second.
    static bool rateAdmitted (double sampleRate) noexcept
    {
        return sampleRate >= kMinSampleRate && sampleRate <= kMaxSampleRate;
    }

    // THE METER'S CAPACITY FOR A PROGRAMME, in samples: the programme plus one second of margin. The ONE place it is
    // decided — both meters this class builds size themselves with it. In samples, and not as `frames / fs + 1`
    // seconds, which was +inf at a finite rate the chain used to accept (1e-305 Hz with 2000 frames): the store's
    // size was then `(std::size_t) inf`, undefined behaviour, and the same ABI call kept 3 blocks on arm64 and wasm32
    // and 4 on x86-64 gcc. Since P51 no such rate reaches this class (`frames / fs` is at most INT_MAX / 8000 s), and since
    // P103 none reaches the meter either; LoudnessConformanceTests pins what is left of the property at the floor.
    static double meterSamples (int frames, double sampleRate) noexcept
    {
        return (double) frames + std::ceil (sampleRate);
    }

    // The loudness meter a programme of `frames` is measured with — the store both meters of this class are built with.
    // 0 for a rate prepare() refuses: the meter itself would take it (and read a rate <= 0 as 48 kHz), but no
    // measurement reaches a meter at it.
    static std::uint64_t meterBytes (double sampleRate, int frames) noexcept
    {
        if (! rateAdmitted (sampleRate)) return 0u;
        analysis::LoudnessMeter::Storage st;
        return analysis::LoudnessMeter::storageFor (sampleRate, meterSamples (frames, sampleRate), st) ? st.bytes() : 0u;
    }

    // EBU Tech 3342 needs short-term samples, one a second: under 3 s of programme there is no range to measure. ONE
    // rule, read by measureInputLoudnessRange() (which refuses before it builds a meter), by its budget, and by
    // `lraValid`.
    static bool rangeMeasurable (int frames, double sampleRate) noexcept
    {
        return frames > 0 && (double) frames / sampleRate >= 3.0;
    }

    // THE TAP BUFFERS' GEOMETRY, as the one function prepare() sizes them with. A `process(n)` call writes
    // K * floor((pos + n) / K) tap frames, which is at most n + K - 1.
    [[nodiscard]] static bool tapLayoutFor (int rendererBlock, int internalBlock, int oversampleFactor,
                                            int& frameCap, int& osCap) noexcept
    {
        const long long cap = (long long) rendererBlock + (long long) internalBlock;
        if (cap > (long long) std::numeric_limits<int>::max() / ((long long) oversampleFactor + 1)) return false;
        frameCap = (int) cap;
        osCap    = (int) (cap * (long long) oversampleFactor);
        return true;
    }

    // How the two peaks are written in dB — the form this class has always reported (it was `TruePeakMeter`'s):
    // `gainToDb` above 1e-10f — the float threshold, widened, so the boundary is the old one to the bit — and -200 for
    // anything quieter. Kept on purpose when the instrument changed (P62), so
    // the switch moves the READING and nothing about how silence is spelled to a caller that tests for it.
    // `gainToDbDet`, and it MUST move together with ReferenceTruePeakMeter::truePeakDb(). The note at the
    // top of this class says the reported number IS the certificate, bit for bit; the certificate is that
    // meter's dB getter, and this is the solver's. Converting one and not the other makes the two spellings
    // of one value disagree at the last ulp — which is not a theory: it is what
    // felitronics_delivered_ceiling_tests caught within one run of the first half of this change, on
    // 44100->88200, 48000->96000, 48000->192000 and 88200->192000.
    // AND THE EQUIVALENCE HAS A FLOOR, which "bit for bit" alone does not say. This function keeps its own
    // 1e-10f gate and its -200.0 sentinel; the meter's getter clamps at core::kGainToDbFloor (1e-12) and
    // reads -240 there. Below the float gate the two therefore differ BY DESIGN, and deliberately — the
    // -200 sentinel is this class's published answer for silence and LoudnessSolverTests pins it. The
    // identity is over peaks above that gate, which is every peak a delivered file has.
    static double peakDb (double lin) noexcept { return lin > kPeakDbGate ? core::gainToDbDet (lin) : kPeakDbSilence; }

    bool measure (float* const* out, int nch, int frames, MasterMeasurement& m, ProgressClock& clock)
    {
        analysis::LoudnessMeter           lm;
        analysis::ReferenceTruePeakMeter  tm;
        if (! lm.prepareForSamples (fs_, nch, meterSamples (frames, fs_))) return false;
        for (int c = 0; c < nch; ++c) lm.setChannelWeight (c, weights_[c]);
        if (! tm.prepare (fs_, frames > 0 ? frames : 1, nch)) return false;
        for (int off = 0; off < frames; )
        {
            const int n = clock.piece (frames - off);
            const float* p[core::kMaxChannels] {};
            for (int c = 0; c < nch; ++c) p[c] = out[c] + off;
            if (! lm.process (p, nch, n)) return false;
            if (! tm.process (p, nch, n)) return false;
            off += n;
            if (! clock.advance (n)) return false;
        }
        tm.drain();                         // its own kTapsPerPhase zeros: the whole FIR, and not given to `lm`

        m.gatingBlocks     = lm.gatingBlockCount();
        m.droppedBlocks    = lm.droppedBlocks();
        // THE uint64 FITS AN int HERE, and what bounds it is this function, not the meter's type. `lm` is a
        // local prepared above — prepareForSamples() ends in reset(), which zeroes the counter — and it sees
        // exactly `frames` samples, in one call or in the clock's pieces (the drain feeds `tm`, never `lm`). The
        // counter moves only in finishSubHop(), by at most one per sub-hop (its two increments are exclusive on
        // `poisoned`), and a sub-hop is at least one sample. So it cannot exceed `frames`, an int. Feeding
        // this meter more than `frames`, or widening `frames`, is what would make this cast wrong.
        m.nonFiniteSubHops = (int) lm.nonFiniteSubHops();
        // THIS MEASUREMENT MIXES THE TWO MATH POLICIES, on purpose and worth saying out loud. The two peak
        // fields below go through `peakDb`, which P80 put on `core::det`, because they are the certificate
        // and the certificate is compared bit for bit. The two loudness fields here come from
        // `analysis::LoudnessMeter`, which is `BasicLoudnessMeter<core::SystemMath>` — the SYSTEM policy —
        // so they are NOT the same bits on every row, and `plrDb` below subtracts one from the other.
        // That is consistent with what this class is measured against: tools/wasm/master-parity.mjs states
        // a TOLERANCE (1e-5 in sample value, 1e-3 dB in the reported numbers), not byte identity, because
        // the two roads render at two roundings by construction. ProgrammeReport, which IS byte-diffed,
        // uses `DeterministicLoudnessMeter` instead. Moving this one would change the solver's SEARCH, not
        // just its report, so it is a decision rather than a tidy-up — recorded here, not done in passing.
        m.integratedLufs   = lm.integratedLufs();
        m.loudnessRangeLu  = lm.loudnessRangeLu();
        m.truePeakDbTp     = peakDb (tm.truePeakLinear());
        m.samplePeakDb     = peakDb (tm.samplePeakLinear());
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
        m.lraValid         = m.loudnessValid && rangeMeasurable (frames, fs_);
        return true;
    }

public:
    // The input's loudness range, for the LRA constraint — which is a DELTA and therefore needs both
    // ends. STATELESS on purpose: it returns the number and the caller puts it in the request, so it
    // cannot outlive the programme it describes. Returns false, and leaves `out` alone, when the
    // programme is too short for the measure to mean anything (EBU Tech 3342 needs short-term samples,
    // one a second) or when the meter could not answer — because `loudnessRangeLu()` returning 0.0 is
    // also what "no dynamic range at all" returns, and a constraint cannot tell those apart.
    [[nodiscard]] bool measureInputLoudnessRange (const float* const* in, int nch, int frames,
                                                  double& out) const
    {
        return measureInputLoudnessRange (in, nch, frames, out, ProgressCallback {});
    }

    [[nodiscard]] bool measureInputLoudnessRange (const float* const* in, int nch, int frames,
                                                  double& out, const ProgressCallback& progress) const
    {
        if (! prepared_ || in == nullptr || nch < 1 || nch > nch_ || frames <= 0) return false;
        if (! rangeMeasurable (frames, fs_)) return false;
        analysis::LoudnessMeter lm;
        if (! lm.prepareForSamples (fs_, nch, meterSamples (frames, fs_))) return false;   // in samples: see meterSamples()
        for (int c = 0; c < nch; ++c) lm.setChannelWeight (c, weights_[c]);
        ProgressClock clock (progress);
        if (! clock.begin (ProgressStage::LoudnessRange, 0, 0, frames, frames)) return false;
        for (int off = 0; off < frames; )
        {
            const int n = clock.piece (frames - off);
            const float* p[core::kMaxChannels] {};
            for (int c = 0; c < nch; ++c) p[c] = in[c] + off;
            if (! lm.process (p, nch, n)) return false;
            off += n;
            if (! clock.advance (n)) return false;
        }
        if (! clock.finish()) return false;
        // `nonFiniteSubHops` belongs in this test and was missing from it, which made the function
        // report SUCCESS on a programme its own meter had already flagged as compromised — and the
        // meter is a local, so the caller could not check for itself. A poisoned 10 ms is recorded as
        // silence, silence fails the absolute gate, and the blocks it was in leave the distribution
        // the range is computed over. Measured on a 30 s programme alternating 3 s loud / 3 s quiet at
        // 48 kHz, with every LOUD second poisoned: 4.8000 LU clean against **21.4000 LU** poisoned,
        // both returned `true`. That number then travels into `LoudnessRequest::inputLoudnessRangeLu`
        // as the far end of the `maxLraLossLu` DELTA, so the constraint is judged against a range the
        // programme does not have. The same condition already guards `MasterMeasurement::loudnessValid`
        // below; the two now agree on what "measured" means.
        if (lm.gatingBlockCount() <= 0 || lm.droppedBlocks() != 0 || lm.nonFiniteSubHops() != 0) return false;
        out = lm.loudnessRangeLu();
        return true;
    }

    // BS.1770 CHANNEL WEIGHTS, forwarded to every meter this class builds. Default 1.0, which is correct
    // for mono and stereo and WRONG for surround: the standard weights Ls/Rs at 1.41 and excludes LFE
    // entirely. Without this the solver would happily steer a 16-channel layout by an LFE-only
    // programme that BS.1770 says has no measurable loudness at all. The host-layout-to-role mapping is
    // product glue and stays outside, exactly as `analysis::LoudnessMeter` says.
    void setChannelWeight (int c, double weight) noexcept
    {
        if (c >= 0 && c < core::kMaxChannels && std::isfinite (weight) && weight >= 0.0) weights_[c] = weight;
    }
    double channelWeight (int c) const noexcept
    {
        return (c >= 0 && c < core::kMaxChannels) ? weights_[c] : 0.0;
    }

private:
    double fs_ = 48000.0;
    int    nch_ = 0, frameCap_ = 0, osCap_ = 0;
    bool   prepared_ = false;

    std::vector<float> compTap_, limTap_, limPeak_;
    dynamics::offline::QuantileHistogram compHist_, limHist_;
    std::uint64_t compActive_ = 0, limActive_ = 0, compFrames_ = 0, limFrames_ = 0;
    float  maxReconLin_ = 0.0f;

    double weights_[core::kMaxChannels] { };
};

} // namespace felitronics::mastering

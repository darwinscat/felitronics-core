// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/analysis/KWeightingFilter.h>
#include <felitronics/analysis/LoudnessMeter.h>
#include <felitronics/analysis/ReferenceTruePeakMeter.h>
#include <felitronics/analysis/StereoColumns.h>
#include <felitronics/core/Config.h>
#include <felitronics/core/Math.h>
#include <felitronics/core/StateGrid.h>
#include <felitronics/eq/Crossover2.h>

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

namespace felitronics::analysis
{

//==============================================================================
// felitronics::analysis::ProgrammeReport — ONE offline pass over a delivered programme, and the continuous
// evidence it leaves: where the silence is, how much the channels differ, what sits under 30 Hz, how far
// the loudness travels. It reports COORDINATES AND NUMBERS, never a verdict: nothing here says "too loud",
// "this is a problem" or "you have rumble". Where the programme cannot support a number, the field says so
// — `valid = false` with a `reason` — instead of publishing a zero that reads as "fine".
//
// WHAT IT MEASURES
//   · DC OFFSET per channel — the mean sample value over the PROGRAMME SPAN (below).
//   · SILENCE at the head and at the tail, twice over: as a count of samples at or below a named THRESHOLD
//     (a parameter, default -96 dBFS), and as a count of EXACT digital silence, which needs no threshold at
//     all. Both are counts; neither is a verdict about whether that silence belongs there.
//   · THE TAIL — the energy of the last `tailWindowMs` of PROGRAMME relative to the programme's own mean
//     square, plus the last sample of the file and the last sample before the trailing digital silence. A
//     fade that was cut off leaves a tail ratio near 1 and a large last signal sample; a fade that
//     completed leaves both small. The instrument states the numbers; what they mean is the reader's.
//   · INFRA-LOW — the fraction of energy below `infraLowHz` (default 30 Hz), where "below" means WHAT THE
//     LR4 PASSES (see THE INFRA-LOW FRACTION IS A FILTER, NOT A BAND below). Not an ideal brick wall, and
//     the difference is large enough to mislead a reader who is not told.
//   · STEREO — the L/R level difference, whether the two channels are bit-identical, the integral
//     correlation over the whole programme (`analysis::StereoSums`, the same arithmetic
//     `fcore_measure correlation` runs), and the side/mid energy ratio.
//   · DYNAMICS — integrated loudness, reference true peak, PLR, LRA (EBU Tech 3342) and the percentiles of
//     the short-term loudness distribution, plus the crest factor.
//
// THE PROGRAMME SPAN, AND WHY EVERY MEAN IS TAKEN OVER IT. The span runs from the first to the last FRAME
// whose loudest finite present channel EXCEEDS `silenceThresholdDb`. Every mean square, every RMS, every
// crest factor and the DC offset are taken over the span, and the tail window ENDS at the last such frame
// rather than at the last sample of the file. The consequence is the property this report is built on:
// PADDING A MASTER DOES NOT MOVE ITS MEASUREMENTS. Pad the same music with five seconds of silence at both
// ends and the DC, the sample peak, the RMS, the crest factor, the infra-low fraction, the programme mean
// square, the tail ratio, the stereo relations and the reference TRUE PEAK all come back bit-identical.
// Measured over the whole file instead, every one of them moves for a reason that has nothing to do with
// the music: the mean square of a 60 s programme padded to 70 s falls by 0.67 dB, and the tail window lands
// entirely inside the padding and reads 0 — a truncated fade reporting as a perfect one.
//
// WHY THE THRESHOLD AND NOT EXACT ZERO, which reads cleaner and was the first design. An exact-zero span is
// a CLIFF exactly where masters live: a 24-bit dithered file — the dominant delivery format — has no exact
// zero anywhere, so its five "silent" seconds are ±2e-7 of dither, an exact-zero span swallows them, the
// RMS is diluted by 5.4 dB, and a programme cut mid-fade reads a tail ratio of 9.1e-12 instead of 1.0. So
// the cliff would sit under most real files and the invariance would apply to almost none of them. Against
// the threshold, |0| > thr is false for every threshold this class accepts, so digital padding stays
// bit-exactly invariant AND sub-threshold padding becomes invariant too.
//
// WHAT PADDING DOES MOVE, named rather than left to be discovered: `lastSample` (which is by definition the
// last sample of the FILE — `lastSignalSample` is its padding-invariant companion), and the LOUDNESS FAMILY
// (integrated, PLR, LRA, the short-term percentiles). BS.1770 anchors its 400 ms gating blocks and 3 s
// windows at the start of the stream, so a pad changes which windows straddle the music's edges even when
// it is a whole number of sub-hops. Anchoring the sub-meters at the first active frame would make these
// invariant too and is deliberately NOT done: it would put this report's loudness at odds with
// `fcore_measure lufs`, ffmpeg's ebur128 and every other EBU tool on the same file, which is a worse trade
// than a padding-sensitive LUFS. The suite pins the whole split as an allow-list.
//
// THE INFRA-LOW FRACTION IS A FILTER, NOT A BAND. The split is `eq::Crossover2` — a 4th-order
// Linkwitz-Riley, i.e. two cascaded Butterworth (Q = 1/√2) sections, so
//
//     |LP_LR4(f)|² = 1 / (1 + (f/fc)⁴)²        — NOT 1/(1 + (f/fc)⁸), which is 4th-order BUTTERWORTH
//
// and the two are not close. LR4 is -6 dB at its crossover (that is what makes low + high reconstruct
// flat), so with fc = 30 Hz a pure 30 Hz tone reads 25 % "below 30 Hz", not 50 %. A pure 40 Hz tone reads
// 5.78 %, a 60 Hz tone 0.346 %, and flat noise to Nyquist reads 0.104 % at 48 kHz (the filter's equivalent
// rectangular bandwidth is 30·3π/(8√2) = 24.99 Hz). So a few per cent "below 30 Hz" on a bass-heavy master
// is very often the FUNDAMENTAL at 40 Hz seen through the skirt of the filter, and not infrasound at all.
// The number is honest and useful; reading it as an ideal band is not. Reported per channel and over the
// programme. Two further properties of the definition, so neither is discovered as a surprise: the ratio
// is Σ lp4(x)² / Σ x² over the span, which equals the |H|²-weighted integral only in the limit — the
// filter's warm-up from `reset()` is inside the window and its ring-down past the last signal frame is
// outside it — and a 100 ms fixture at 40 Hz is four cycles, dominated by that warm-up rather than by the
// steady-state 5.78 %.
//
// "CANNOT SAY" IS A RESULT. A mono file does not get a correlation of +1.0 — it gets
// `valid = false, reason = MonoProgramme`. A programme shorter than the tail window does not get a ratio of
// exactly 1.0. A programme with fewer than two gated short-term observations does not get an LRA of 0.0,
// which is also the honest answer for a constant tone. An invalid field's `value` is a documented
// canonical +0.0 — never a NaN, which this repository does not use as an absence marker (a computed NaN's
// sign bit is not even portable: 0x7ff8… on arm64, 0xfff8… on x86-64).
//
// TWO PRIMITIVES ARE WRAPPED BECAUSE THEIR OWN DEFAULTS ARE VERDICTS, and that is worth naming:
//   · `StereoSums::correlation()` answers +1.0 when √(ΣLL·ΣRR) ≤ 1e-12 — "silence is neutral"
//     (StereoColumns.h:56). For a page drawing a meter that is right; for a report it says "perfectly
//     mono-compatible" about a programme that carries no stereo information at all. Here that case is
//     `valid = false, reason = SilentProgramme`, and the NUMBER, when valid, is that primitive's own, bit
//     for bit, so this report and `fcore_measure correlation` cannot disagree.
//   · `LoudnessMeter::integratedLufs()` answers -120.0 when no block passed the gates. That is a SENTINEL,
//     not a measurement, and it is invalid here.
//   · And LRA is NOT read from `LoudnessMeter::loudnessRangeLu()`, for the same class of reason: that
//     function answers 0.0 both for a constant tone (whose range really is 0 LU) and for a programme where
//     the gates left fewer than two short-term observations, and from outside the meter the two cannot be
//     told apart. EBU Tech 3342 is computed here from this class's own short-term series, where the gated
//     count is known exactly. The two agree to about 0.3 LU wherever both are defined — not bitwise, and
//     the reason is worth knowing: the meter sums its short-term window newest-first
//     (LoudnessMeter.h:325) and this class oldest-first, so the two orders leave different residues. The
//     suite pins the agreement at that tolerance.
//
// FORM. setParams / prepare / process / finish / reset, like `analysis::ClipDetector`. `process()` is
// READ-ONLY and allocates nothing; `finish()` drains the true-peak filter, closes the report and freezes
// it, is idempotent, and `process()` refuses until `reset()`. Every parameter is read at `prepare()` and
// not before — so the report can never depend on WHEN the caller called `setParams()`.
//
// LAW 8a — BIT-IDENTICAL UNDER ANY RE-SPLIT. With the same input bits, the same parameters, the same
// binary and the same per-sample timeline of channel presence, every field of the report — and every
// intermediate the trace hook below exposes — is bit-identical however the caller cuts the stream into
// `process()` calls. This is stronger than law 11(a), which only promises agreement with the caller
// chunking at the SAME boundaries (DSP-ARCHITECTURE.md:273), and it is what the design is arranged around:
//   · ONE integer clock. `totalSamples_` counts samples since `reset()`; the outer loop is over SAMPLES and
//     the inner over CHANNELS, because a cross-channel accumulator filled channel-outer adds in a different
//     order for a 2-sample call than for two 1-sample calls, and a different order is different bits.
//   · `maxBlock` sizes the scratch and NOTHING else — no window, no hop, no branch. A call longer than
//     `maxBlock` is consumed IN FULL; a refused call consumes nothing; `n == 0` changes nothing at all.
//   · Sub-hop boundaries are found by `totalSamples_ == nextSubHopEnd_`, never by a modulo.
//   · The denormal flush of the infra-low crossover rides `core::StateGrid` (64 samples of AUDIO time), not
//     the end of `process()` — `eq::Crossover2::processSample` does not clock itself (Crossover2.h:52), so
//     the owner must. The K-weighting flush rides the 10 ms sub-hop, the cadence `LoudnessMeter` documents
//     for the same filter (its binding constant is the 90 ms shelf, not the 2.85 s RLB).
//   · The sliding 3 s window is RECOMPUTED from its ring at each observation, never carried as
//     `sum += new - old`, which is deterministic but drifts over a million hops.
//   · Every accumulator is one fixed recurrence in `double`; the compensated (Neumaier) sums hold two
//     constant scalars and nothing that depends on where a call ended.
//   · `LoudnessMeter` and `ReferenceTruePeakMeter` both carry a falling-edge action keyed on the CALL's
//     width (their `ranNc_`). They are therefore driven here at a CONSTANT width through a scratch buffer,
//     with canonical zeros for a channel the caller stopped delivering, so that edge never fires. Feeding
//     silence is also the honest state for a whole-programme report: law 11c's "the channel heard
//     silence", and exactly what `ReferenceTruePeakMeter::drain()` feeds a stopped channel anyway.
// NOT PROMISED: bit-identity between platforms, toolchains or libms. `log10`, `pow` and `sqrt` are not
// bit-portable, and the percentile bin and the loudness conversions are computed through them.
//
// NON-FINITE INPUT. A non-finite sample is a HOLE: a canonical 0.0f goes into every filter, the sample
// enters no statistic, and it is counted per channel. What that invalidates is drawn along one line —
// A HOLE IS EXCLUDED WHERE EXCLUSION IS EXACT, AND FATAL WHERE STATE CARRIES IT FORWARD:
//   · DC, RMS, crest factor, the stereo sums and the tail window are sums over the finite samples; dropping
//     one term is exact, so they stay valid and the counts are published beside them.
//   · The infra-low fraction, the loudness family and the short-term percentiles run through recursive
//     filters that carry the substitution forward for their whole ring-down. Any non-finite input makes
//     them `valid = false, reason = NonFiniteInput`.
// AND FINITE INPUT IS NOT ENOUGH. `eq::Svf` updates its integrators as `ic1 = (float)(2·v1 - ic1)`, and
// that intermediate can leave the float range even where the state itself is not growing — a stream of
// finite FLT_MAX drives the 30 Hz crossover non-finite at sample 849 exactly, even with the grid flush
// running (a 1.7e38 DC input does NOT overflow, so "the state doubles every sample" would be the wrong
// mechanism to quote)
// — measured. The K-weighting biquads and `StereoSums`'s products overflow the same way. So every
// intermediate is tested for finiteness where it is produced, the event is counted, and the field it feeds
// becomes `valid = false, reason = NonFiniteIntermediate`. Nothing non-finite reaches an accumulator or
// the trace. This is also why `LoudnessMeter::nonFiniteSubHops()` is published even though this class
// sanitises what it feeds the meter: that counter still fires on a FINITE input that overflows the
// K-weighting, which an input-side count cannot see.
//
// CAPACITY IS DATA, NOT CONTROL. `LoudnessMeter` keeps every gating block to the end (its relative gate is
// a two-pass) and drops what will not fit, counting it in `droppedBlocks()`; the short-term series kept
// here for LRA has its own capacity. Both are sized from `maxDurationSec` and both are PUBLISHED by
// `storageFor()` before a byte is allocated (law 11d). Overflow never stops the measurement and never
// refuses a call — it invalidates exactly the fields it damages, with `reason = LoudnessCapacityExceeded`,
// and every counter keeps counting.

//==============================================================================
// The parameters. Each is read at `prepare()` and nowhere else: a parameter that took effect mid-stream
// would make the report a function of the caller's timeline, which is precisely what law 8a forbids.
struct ProgrammeReportParams
{
    // The level a sample must EXCEED to count as signal, in dBFS, compared against |x| per frame.
    // -96 dBFS = 1.5849e-5 linear. WHY NOT -90: one LSB of 16-bit PCM is 1/32768 = 3.0518e-5 = -90.309
    // dBFS, which is BELOW a -90 dBFS threshold — so a dithered 16-bit fade-out would be reported as
    // silence, and a 16-bit file whose left channel sits at exactly ±1 LSB would read as 100 % silent.
    // -96 dBFS sits 5.7 dB under that LSB, so every 16-bit dither tail is signal. At 20 bit (1 LSB =
    // -114.39 dBFS) and 24 bit (-138.47 dBFS) a dither tail is BELOW this default and does read as
    // silence: for a depth-independent answer use the digital-silence counts, which have no threshold.
    double silenceThresholdDb = -96.0;

    // The tail window, in milliseconds. It ENDS at the last signal frame, not at the last sample.
    double tailWindowMs = 100.0;

    // The infra-low crossover, in Hz. Must satisfy 1 <= infraLowHz <= 0.45·sampleRate — `eq::Svf` silently
    // CLAMPS a frequency outside its own range (Svf.h:61-62), and a clamped report would be a report about
    // a filter the caller did not ask for, so it is refused here instead.
    double infraLowHz = 30.0;

    // The programme length the loudness stores are sized for, in seconds. A longer programme keeps being
    // measured; the fields the overflow damages go invalid.
    double maxDurationSec = 3600.0;
};

// Why a field could not be measured. Never "what is wrong with the programme".
enum class ProgrammeReason : std::uint8_t
{
    None = 0,                      // the field is valid
    MonoProgramme,                 // fewer than two prepared channels: no stereo relation exists
    NoStereoFrames,                // prepared for stereo, but no frame carried both channels finite
    SilentProgramme,               // the quantity needs signal and the programme has none
    NoFiniteSamples,               // every sample of the region was a hole or absent
    ShorterThanTailWindow,         // the programme span is shorter than tailWindowMs
    ShorterThanLoudnessWindow,     // shorter than the window the measure is defined over
    TooFewObservations,            // fewer than two short-term observations survived the gates
    LoudnessCapacityExceeded,      // a store sized by maxDurationSec overflowed; the answer would be partial
    NonFiniteInput,                // a non-finite INPUT sample reached a filter that carries it forward
    NonFiniteIntermediate,         // a finite input overflowed a filter or an accumulator
    ZeroDenominator                // the ratio's denominator is zero
};

// A measured scalar. `value` is a canonical +0.0 whenever `valid` is false — never a NaN.
struct ProgrammeValue
{
    double          value  = 0.0;
    bool            valid  = false;
    ProgrammeReason reason = ProgrammeReason::None;
};

// A measured count. Counts are facts, but a count can still be untrustworthy — a silence run with a hole
// inside it is not a measured silence — so it carries the same verdict.
struct ProgrammeCount
{
    std::int64_t    value  = 0;
    bool            valid  = true;
    ProgrammeReason reason = ProgrammeReason::None;
};

//==============================================================================
// THE TRACE HOOK — a test surface, and it is here because the final report is NOT a sufficient gate for
// law 8a. This project has its own measured example: docs/LAW8-KWEIGHTING.md:78, where the final LUFS and
// dBTP matched BIT FOR BIT between two slicings while 5 of 97 intermediate block energies had moved. A
// re-splitting defect contracts inside a maximum, a percentile, a gate, a histogram and a rounding; a
// shifted window on a stationary tone gives the same spectrum; a lost frame and a duplicated identical one
// leave the mean alone. So the gate compares the whole SERIES, element by element, at the moment each one
// CLOSES — not what fell out of `process()`.
//
// The buffer is the caller's and is never resized: `setTraceBuffer()` takes a pointer and a capacity,
// `process()` writes up to that capacity and counts what did not fit. With no buffer set, nothing is
// recorded and nothing is paid. It allocates nothing, so a suite can hold the no-allocation check and the
// trace at the same time.
enum class ProgrammeTraceKind : std::uint8_t
{
    SubHopChannel = 1,   // per (closed sub-hop, channel): a = this hop's Σx², b = this hop's Σlp4(x)², c = its ΣK²
    SubHop        = 2,   // per closed sub-hop:            a = channel-weighted mean square, b = finite samples, c = 0
    ShortTerm     = 3,   // per short-term observation:    a = mean square (RAW, pre-log), b = histogram bin, c = 1 if abs-gated
    SilenceEdge   = 4,   // a predicate transition:        a = 0 thresholded / 1 digital, b = 1 entering signal, c = frame max |x|
    TailWindow    = 5    // once, from finish():           a = tail Σx², b = tail sample count, c = programme mean square
};

struct ProgrammeTraceEvent
{
    ProgrammeTraceKind kind    = ProgrammeTraceKind::SubHop;
    std::int32_t       channel = -1;          // -1 = the whole frame
    std::int64_t       index   = 0;           // sub-hop / observation / edge ordinal
    std::int64_t       begin   = 0;           // [begin, end) in samples since reset()
    std::int64_t       end     = 0;
    double             a       = 0.0;
    double             b       = 0.0;
    double             c       = 0.0;
};

//==============================================================================
class ProgrammeReport
{
public:
    static constexpr double kMinSampleRate      = 1000.0;
    static constexpr double kMaxSampleRate      = 768000.0;
    static constexpr int    kShortTermSubHops   = 300;    // 3 s of 10 ms sub-hops — the short-term window
    static constexpr int    kObservationHops    = 100;    // one short-term observation a second (libebur128's cadence)
    static constexpr int    kHistogramBins      = 1000;   // -70 … +30 LUFS in 0.1 LU bins (LoudnessMeter.h:366-385)
    static constexpr double kHistogramFloorLufs = -70.0;
    static constexpr double kHistogramBinLu     = 0.1;
    // The BS.1770 absolute gate, -70 LUFS, as ENERGY: 10^((-70 + 0.691)/10). Written as a literal rather
    // than computed with `std::pow` because it DECIDES INCLUSION, and a libm that disagrees by an ulp would
    // move a borderline observation in on one platform and out on another. Exact in binary64.
    static constexpr double kAbsoluteGateEnergy = 0x1.f791ec6e1d5b7p-24;   // == 1.1724653045822981e-07
    static constexpr double kRelativeGate       = 0.01;   // -20 LU, as an energy factor (EBU Tech 3342)
    static constexpr double kSilenceFloorLufs   = -120.0; // LoudnessMeter's own "no signal" sentinel

    //==========================================================================
    // THE REPORT. Produced by `finish()` and frozen there; before it, every field is invalid.
    struct PerChannel
    {
        ProgrammeValue dcOffset;            // mean sample value over the programme span
        ProgrammeValue samplePeak;          // max |x|, linear
        ProgrammeValue rms;                 // √(mean square over the span), linear
        ProgrammeValue crestFactorDb;       // 20·log10(samplePeak / rms)
        ProgrammeValue infraLowFraction;    // Σ lp4(x)² / Σ x² over the span — see the header
        ProgrammeValue lastSample;          // the sample at T-1, exactly as delivered
        ProgrammeValue lastSignalSample;    // the last sample this channel carried on a signal frame
        std::int64_t   finiteSamples    = 0;
        std::int64_t   nonFiniteSamples = 0;
        std::int64_t   absentSamples    = 0;  // frames where this channel was outside the call's width
    };

    struct Report
    {
        int            channels     = 0;
        std::int64_t   totalSamples = 0;
        double         sampleRate   = 0.0;

        std::array<PerChannel, core::kMaxChannels> channel {};

        // --- the whole programme ---
        ProgrammeValue samplePeak, rms, crestFactorDb, infraLowFraction, programmeMeanSquare;
        std::int64_t   programmeSpanSamples = 0;   // [firstSignal, lastSignal] inclusive, in frames
        std::int64_t   programmeSpanFinite  = 0;   // finite present samples inside it, all channels

        // --- silence ---
        ProgrammeCount leadingSilenceSamples, trailingSilenceSamples;      // at silenceThresholdDb
        std::int64_t   leadingDigitalSilenceSamples  = 0;                  // exact zeros: no threshold
        std::int64_t   trailingDigitalSilenceSamples = 0;
        std::int64_t   silenceEdges = 0;                                   // predicate transitions seen

        // --- the tail ---
        ProgrammeValue tailEnergyRatio;                 // mean square of the tail window / of the programme
        std::int64_t   tailWindowSamples = 0;           // the window's length in frames
        std::int64_t   tailWindowFinite  = 0;           // finite present samples inside it

        // --- stereo (channels 0 and 1) ---
        std::int64_t   stereoFrames             = 0;    // frames where both were present and finite
        std::int64_t   stereoBitIdenticalFrames = 0;    // …and their bit patterns matched
        std::int64_t   stereoValueEqualFrames   = 0;    // …and their VALUES matched (+0.0 == -0.0 here)
        ProgrammeValue stereoBalanceDb;                 // 10·log10(ΣLL / ΣRR)
        ProgrammeValue stereoCorrelation;               // ΣLR / √(ΣLL·ΣRR) — analysis::StereoSums
        ProgrammeValue stereoSideToMidRatio;            // Σ side² / Σ mid², the ½(L±R) convention

        // --- dynamics ---
        ProgrammeValue integratedLufs;                  // BS.1770 gated, analysis::LoudnessMeter
        ProgrammeValue truePeakDbtp;                    // analysis::ReferenceTruePeakMeter, 4× / 32 taps
        ProgrammeValue plrDb;                           // truePeakDbtp - integratedLufs
        ProgrammeValue lraLu;                           // EBU Tech 3342, P95-P10 of the gated short-term set
        ProgrammeValue shortTermP10, shortTermP50, shortTermP95;   // absolute-gated only
        ProgrammeValue shortTermSpreadLu;               // P95 - P10 of the same set
        std::int64_t   shortTermObservations      = 0;
        std::int64_t   shortTermGatedObservations = 0;  // …that passed the -70 LUFS absolute gate
        std::int64_t   uncoveredSubHopSamples     = 0;  // T - (closed sub-hops)·subHopSamples

        // --- damage and capacity ---
        std::int64_t   nonFiniteInputSamples        = 0;  // all channels
        std::int64_t   nonFiniteIntermediates       = 0;  // finite input that overflowed a filter or a sum
        std::int64_t   loudnessDroppedBlocks        = 0;  // LoudnessMeter::droppedBlocks()
        std::int64_t   loudnessNonFiniteSubHops     = 0;  // LoudnessMeter::nonFiniteSubHops()
        std::int64_t   shortTermDroppedObservations = 0;
        // NOT traceOverflow: how many events the caller's buffer could not hold is a property of that
        // buffer and not of the programme, and holding it here made the REPORT a function of the trace's
        // capacity — which the law-8a comparison then had to zero out by hand. It is `traceOverflow()`.

        // THE ONE ENUMERATION OF THE FIELDS, so that a consumer cannot fall behind the struct. The law-8a
        // gate compares through it and `fcore_measure report` prints through it, so a field added here is
        // automatically compared and automatically printed; a hand-written list in either place would
        // silently stop covering the new field, which is how a gate quietly stops being a gate.
        // `f (name, channel, value)` — channel is -1 for a whole-programme field.
        template <typename F>
        void visitValues (F&& f) const
        {
            for (int c = 0; c < channels; ++c)
            {
                const PerChannel& p = channel[(std::size_t) c];
                f ("dcOffset",         c, p.dcOffset);
                f ("samplePeak",       c, p.samplePeak);
                f ("rms",              c, p.rms);
                f ("crestFactorDb",    c, p.crestFactorDb);
                f ("infraLowFraction", c, p.infraLowFraction);
                f ("lastSample",       c, p.lastSample);
                f ("lastSignalSample", c, p.lastSignalSample);
            }
            f ("samplePeak",           -1, samplePeak);
            f ("rms",                  -1, rms);
            f ("crestFactorDb",        -1, crestFactorDb);
            f ("infraLowFraction",     -1, infraLowFraction);
            f ("programmeMeanSquare",  -1, programmeMeanSquare);
            f ("tailEnergyRatio",      -1, tailEnergyRatio);
            f ("stereoBalanceDb",      -1, stereoBalanceDb);
            f ("stereoCorrelation",    -1, stereoCorrelation);
            f ("stereoSideToMidRatio", -1, stereoSideToMidRatio);
            f ("integratedLufs",       -1, integratedLufs);
            f ("truePeakDbtp",         -1, truePeakDbtp);
            f ("plrDb",                -1, plrDb);
            f ("lraLu",                -1, lraLu);
            f ("shortTermP10",         -1, shortTermP10);
            f ("shortTermP50",         -1, shortTermP50);
            f ("shortTermP95",         -1, shortTermP95);
            f ("shortTermSpreadLu",    -1, shortTermSpreadLu);
        }

        // `f (name, channel, count)` — the integer half of the same enumeration.
        template <typename F>
        void visitCounts (F&& f) const
        {
            f ("channels",     -1, (std::int64_t) channels);
            f ("totalSamples", -1, totalSamples);
            for (int c = 0; c < channels; ++c)
            {
                const PerChannel& p = channel[(std::size_t) c];
                f ("finiteSamples",    c, p.finiteSamples);
                f ("nonFiniteSamples", c, p.nonFiniteSamples);
                f ("absentSamples",    c, p.absentSamples);
            }
            f ("programmeSpanSamples",          -1, programmeSpanSamples);
            f ("programmeSpanFinite",           -1, programmeSpanFinite);
            f ("leadingSilenceSamples",         -1, leadingSilenceSamples.value);
            f ("leadingSilenceValid",           -1, (std::int64_t) (leadingSilenceSamples.valid ? 1 : 0));
            f ("leadingSilenceReason",          -1, (std::int64_t) leadingSilenceSamples.reason);
            f ("trailingSilenceSamples",        -1, trailingSilenceSamples.value);
            f ("trailingSilenceValid",          -1, (std::int64_t) (trailingSilenceSamples.valid ? 1 : 0));
            f ("trailingSilenceReason",         -1, (std::int64_t) trailingSilenceSamples.reason);
            f ("leadingDigitalSilenceSamples",  -1, leadingDigitalSilenceSamples);
            f ("trailingDigitalSilenceSamples", -1, trailingDigitalSilenceSamples);
            f ("silenceEdges",                  -1, silenceEdges);
            f ("tailWindowSamples",             -1, tailWindowSamples);
            f ("tailWindowFinite",              -1, tailWindowFinite);
            f ("stereoFrames",                  -1, stereoFrames);
            f ("stereoBitIdenticalFrames",      -1, stereoBitIdenticalFrames);
            f ("stereoValueEqualFrames",        -1, stereoValueEqualFrames);
            f ("shortTermObservations",         -1, shortTermObservations);
            f ("shortTermGatedObservations",    -1, shortTermGatedObservations);
            f ("uncoveredSubHopSamples",        -1, uncoveredSubHopSamples);
            f ("nonFiniteInputSamples",         -1, nonFiniteInputSamples);
            f ("nonFiniteIntermediates",        -1, nonFiniteIntermediates);
            f ("loudnessDroppedBlocks",         -1, loudnessDroppedBlocks);
            f ("loudnessNonFiniteSubHops",      -1, loudnessNonFiniteSubHops);
            f ("shortTermDroppedObservations",  -1, shortTermDroppedObservations);
        }
    };

    //==========================================================================
    // WHAT prepare() ASKS THE HEAP FOR, from the function prepare() sizes and validates itself with, so the
    // budget and the allocation cannot drift (law 11d). Asked of a FRESH object; a prepared one keeps
    // storage that still fits. `ok == false` on exactly the arguments prepare() refuses, and then every
    // number below is zero. The 3 s sub-hop ring and the percentile histogram are NOT here: both are fixed
    // and inline, so they are not heap at all.
    struct Storage
    {
        bool          ok               = false;
        int           channels         = 0;
        std::int64_t  tailSamples      = 0;    // frames in the tail window
        std::int64_t  subHopSamples    = 0;    // lround(0.01·fs), at least 1 — the 10 ms sub-hop
        std::size_t   tailRingEnergies = 0;    // doubles
        std::size_t   tailRingCounts   = 0;    // int32
        std::size_t   tailPendingCounts = 0;   // int32 — the deferred quiet run's per-frame counts
        std::size_t   tailPendingEnergies = 0; // doubles — …and their energies, which are not always zero
        std::size_t   scratchFloats    = 0;    // maxBlock · channels — the constant-width feed
        std::size_t   shortTermEntries = 0;    // doubles: one short-term observation a second
        LoudnessMeter::Storage          loudness {};
        ReferenceTruePeakMeter::Storage truePeak {};

        // 64 bits because the product is the point: on wasm32 a size_t byte count wraps long before the
        // element counts do.
        std::uint64_t bytes() const noexcept
        {
            return (std::uint64_t) sizeof (double)       * (std::uint64_t) tailRingEnergies
                 + (std::uint64_t) sizeof (std::int32_t) * (std::uint64_t) tailRingCounts
                 + (std::uint64_t) sizeof (std::int32_t) * (std::uint64_t) tailPendingCounts
                 + (std::uint64_t) sizeof (double)       * (std::uint64_t) tailPendingEnergies
                 + (std::uint64_t) sizeof (float)        * (std::uint64_t) scratchFloats
                 + (std::uint64_t) sizeof (double)       * (std::uint64_t) shortTermEntries
                 + loudness.bytes() + truePeak.bytes();
        }
    };

    static Storage storageFor (double sampleRate, int maxBlock, int maxChannels,
                               const ProgrammeReportParams& p) noexcept
    {
        Storage st;
        if (! (sampleRate >= kMinSampleRate && sampleRate <= kMaxSampleRate)) return st;   // NaN fails too
        if (maxBlock < 1 || maxChannels < 1 || maxChannels > core::kMaxChannels) return st;
        if (! (p.tailWindowMs > 0.0 && p.tailWindowMs <= 60000.0)) return st;
        if (! (p.infraLowHz >= 1.0 && p.infraLowHz <= 0.45 * sampleRate)) return st;       // Svf would CLAMP
        if (! (p.silenceThresholdDb <= 0.0 && p.silenceThresholdDb >= -400.0)) return st;
        if (! (p.maxDurationSec > 0.0 && p.maxDurationSec <= 360000.0)) return st;

        const std::int64_t tail = (std::int64_t) std::lround (sampleRate * p.tailWindowMs / 1000.0);
        if (tail < 1) return st;
        const std::int64_t sub = std::max<std::int64_t> (1, (std::int64_t) std::lround (0.01 * sampleRate));

        const double maxSamples = p.maxDurationSec * sampleRate;
        if (! LoudnessMeter::storageFor (sampleRate, maxSamples, st.loudness)) return st;
        st.truePeak = ReferenceTruePeakMeter::storageFor (sampleRate, maxBlock, maxChannels);
        if (! st.truePeak.ok) return st;

        // One observation per kObservationHops sub-hops, the first at kShortTermSubHops. The +8 is margin,
        // not need — the shape LoudnessMeter::storageFor uses for its own short-term store.
        const double obs = std::ceil (maxSamples) / (double) (sub * (std::int64_t) kObservationHops);
        if (! (obs <= 1.0e8)) return st;                       // a count no sane store would hold

        st.ok               = true;
        st.channels         = maxChannels;
        st.tailSamples      = tail;
        st.subHopSamples    = sub;
        st.tailRingEnergies = (std::size_t) tail;
        st.tailRingCounts   = (std::size_t) tail;
        st.tailPendingCounts = (std::size_t) tail;
        st.tailPendingEnergies = (std::size_t) tail;
        st.scratchFloats    = (std::size_t) maxBlock * (std::size_t) maxChannels;
        st.shortTermEntries = (std::size_t) obs + 8u;
        return st;
    }

    void setParams (const ProgrammeReportParams& p) noexcept { params_ = p; }   // read at the next prepare()
    const ProgrammeReportParams& params() const noexcept { return params_; }

    // Law 11d: disarm, validate, size with the PUBLIC storageFor(), allocate, reset.
    [[nodiscard]] bool prepare (double sampleRate, int maxBlock, int maxChannels)
    {
        prepared_ = false;
        const Storage st = storageFor (sampleRate, maxBlock, maxChannels, params_);
        if (! st.ok) return false;

        sampleRate_    = sampleRate;
        maxBlock_      = maxBlock;
        channels_      = maxChannels;
        tailSamples_   = st.tailSamples;
        subHopSamples_ = st.subHopSamples;
        silenceThr_    = std::pow (10.0, params_.silenceThresholdDb / 20.0);

        tailSq_.assign (st.tailRingEnergies, 0.0);
        tailN_.assign (st.tailRingCounts, 0);
        tailPendN_.assign (st.tailPendingCounts, 0);
        tailPendSq_.assign (st.tailPendingEnergies, 0.0);
        scratch_.assign (st.scratchFloats, 0.0f);
        shortTerm_.assign (st.shortTermEntries, 0.0);

        if (! lm_.prepareForSamples (sampleRate, channels_, params_.maxDurationSec * sampleRate)) return false;
        if (! tp_.prepare (sampleRate, maxBlock, channels_)) return false;
        kw_.prepare (sampleRate, channels_);
        lr4_.prepare (sampleRate, channels_);
        lr4_.setFrequency ((float) params_.infraLowHz);

        prepared_ = true;
        reset();
        return true;
    }

    bool isPrepared() const noexcept { return prepared_; }
    bool isFinished() const noexcept { return finished_; }
    int  channels()   const noexcept { return channels_; }
    std::int64_t samplesProcessed()  const noexcept { return totalSamples_; }
    std::int64_t subHopSamples()     const noexcept { return subHopSamples_; }
    std::int64_t tailWindowSamples() const noexcept { return tailSamples_; }
    static constexpr int latencySamples() noexcept { return 0; }               // a read-only sink

    // Re-anchors EVERYTHING: the clock, the rings, the filters, the grid phase, the sub-meters, every
    // counter, every damage flag and the frozen report.
    void reset() noexcept
    {
        for (auto& c : chans_) c = ChannelState {};
        std::fill (tailSq_.begin(), tailSq_.end(), 0.0);
        std::fill (tailN_.begin(), tailN_.end(), 0);
        std::fill (tailPendN_.begin(), tailPendN_.end(), 0);
        std::fill (tailPendSq_.begin(), tailPendSq_.end(), 0.0);
        std::fill (subRing_.begin(), subRing_.end(), 0.0);

        totalSamples_ = 0;
        nextSubHopEnd_ = subHopSamples_;
        subHopIndex_ = 0; nextObservation_ = kShortTermSubHops; observationIndex_ = 0;
        tailWritten_ = 0; tailPendingZeros_ = 0;
        firstActive_ = -1; lastActive_ = -1; firstNonZero_ = -1; lastNonZero_ = -1;
        nfBeforeFirstActive_ = 0; nfAfterLastActive_ = 0;
        wasActive_ = false; wasNonZero_ = false; edges_ = 0;
        stereoFrames_ = 0; stereoBitEqual_ = 0; stereoValueEqual_ = 0;
        sums_ = StereoSums {};
        shortTermCount_ = 0; shortTermDropped_ = 0;
        nonFiniteInput_ = 0; nonFiniteIntermediate_ = 0; kwOverflow_ = 0;
        traceCount_ = 0; traceOverflow_ = 0;
        grid_.reset();
        kw_.reset();
        lr4_.reset();
        lm_.reset();
        tp_.reset();
        report_ = Report {};
        finished_ = false;
    }

    // The test surface — see THE TRACE HOOK above. `capacity` events at most, never resized, never
    // allocated. Passing nullptr turns recording off. It records nothing that has already happened, so a
    // fixture sets it before the first `process()`.
    void setTraceBuffer (ProgrammeTraceEvent* buffer, std::size_t capacity) noexcept
    {
        trace_    = buffer;
        traceCap_ = buffer != nullptr ? capacity : 0u;
    }
    std::int64_t traceCount()    const noexcept { return traceCount_; }
    std::int64_t traceOverflow() const noexcept { return traceOverflow_; }

    // READ-ONLY, allocation-free. Law 11: malformed → unprepared → finished → nch > maxChannels → a null
    // plane → n == 0 → run. A legal call longer than maxBlock is consumed IN FULL, in maxBlock pieces; a
    // refused one consumes nothing.
    [[nodiscard]] bool process (const float* const* in, int numChannels, int n) noexcept
    {
        if (numChannels < 0 || n < 0) return false;
        if (! prepared_ || finished_) return false;
        if (numChannels > channels_) return false;
        if (numChannels > 0)
        {
            if (in == nullptr) return false;
            for (int c = 0; c < numChannels; ++c) if (in[c] == nullptr) return false;
        }
        if (n == 0) return true;                                  // no samples: no time, no edge

        for (int off = 0; off < n; )
        {
            const int m = std::min (maxBlock_, n - off);
            runChunk (in, numChannels, off, m);
            off += m;
        }
        return true;
    }

    // End of stream. Drains the reference true-peak filter, closes the tail window, converts to dB, ratios
    // and percentiles exactly once, and freezes the report. Idempotent.
    void finish() noexcept
    {
        if (! prepared_ || finished_) return;
        tp_.drain();                                              // the last ~16 samples are still in the FIR
        buildReport();
        finished_ = true;
    }

    // The report. Meaningful after finish(); before it, every field is invalid and every count is zero.
    const Report& report() const noexcept { return report_; }

private:
    //==========================================================================
    // A Neumaier compensated sum: two constant scalars and one fixed recurrence, which is what law 8a's
    // "one recurrence per quantity" asks for — `ClipDetector.h:393` keeps its DC the same way. Legal
    // precisely because nothing about it depends on where the caller cut the stream.
    struct Sum
    {
        double s = 0.0, c = 0.0;
        inline void add (double x) noexcept
        {
            const double t = s + x;
            c += (std::fabs (s) >= std::fabs (x)) ? ((s - t) + x) : ((x - t) + s);
            s = t;
        }
        double value() const noexcept { return s + c; }
    };

    struct ChannelState
    {
        // THE THREE SPAN ACCUMULATORS, and all three are COMMITTED rather than summed over everything.
        // With an exact-zero span the distinction did not exist: a frame outside the span contributed
        // exactly 0 to Σx and Σx², so only the COUNT had to be corrected. A THRESHOLDED span changes that —
        // a sub-threshold frame is quiet, not silent, and its energy is real — so a numerator summed over
        // every sample beside a denominator counted over the span only is a mean square of one region
        // divided by the length of another. Measured when the span moved to the threshold: the RMS of a
        // 1e-4 tone disagreed with the whole-file reference in the 7th digit, which is the leak.
        Sum          dc;                      // Σ x over the span's finite present samples
        Sum          sumSq;                   // Σ x² over them
        Sum          lowSqSpan;               // Σ lp4(x)² over them
        double       dcPending    = 0.0;      // …the quiet run since the last active frame, not yet committed
        double       sumSqPending = 0.0;
        double       lowSqPending = 0.0;
        double       hopSq       = 0.0;       // Σ x² inside the sub-hop in progress — a trace quantity
        double       hopLowSq    = 0.0;       // Σ lp4(x)² inside it
        double       kwSq        = 0.0;       // the K-weighted square inside it
        double       peak        = 0.0;
        float        last        = 0.0f;
        bool         lastFinite  = false;
        float        lastSignal  = 0.0f;
        bool         lastSignalFinite = false;
        bool         hasLastSignal    = false;
        bool         lowFinite   = true;      // false once this channel's LR4 output overflowed
        std::int64_t finite = 0, nonFinite = 0, absent = 0;
        std::int64_t finiteLead  = 0;         // finite samples in frames before the first signal frame
        std::int64_t finiteTrail = 0;         // …since the last signal frame
    };

    //==========================================================================
    void runChunk (const float* const* in, int nch, int off, int m) noexcept
    {
        // The constant-width feed. Both sub-meters carry a falling edge keyed on the call's width; feeding
        // them a WIDTH THAT NEVER CHANGES is what keeps that edge from firing at a caller-chosen moment.
        // The samples are SANITISED here: `ReferenceTruePeakMeter` has no poison flush, so one +Inf would
        // sit in its polyphase ring and make the true peak +Inf for the rest of the programme with no way
        // back. What the substitution cannot hide is published beside the reading.
        for (int c = 0; c < channels_; ++c)
        {
            float* dst = scratch_.data() + (std::size_t) c * (std::size_t) maxBlock_;
            if (c < nch)
            {
                const float* src = in[c] + off;
                for (int i = 0; i < m; ++i) { const float v = src[i]; dst[i] = std::isfinite (v) ? v : 0.0f; }
            }
            else
            {
                for (int i = 0; i < m; ++i) dst[i] = 0.0f;
            }
        }

        for (int i = 0; i < m; ++i)
        {
            const std::int64_t t = totalSamples_;                 // the sample's absolute index
            double frameSumSq = 0.0, frameAbsMax = 0.0;
            int    frameFinite = 0;
            bool   frameNonFinite = false, frameNonZero = false;

            // This frame's per-channel contributions, held until the frame's SIGNAL predicate is known —
            // the span's trailing run is a property of the frame, not of a channel.
            double frameX[core::kMaxChannels] {}, frameSq[core::kMaxChannels] {};
            double frameLowSq[core::kMaxChannels] {};
            bool   framePresentFinite[core::kMaxChannels] {};

            for (int c = 0; c < channels_; ++c)                   // sample OUTSIDE, channel INSIDE
            {
                ChannelState& ch = chans_[(std::size_t) c];
                const bool  present = c < nch;
                const float raw     = present ? in[c][off + i] : 0.0f;
                const bool  finite  = present && std::isfinite (raw);
                const float x       = scratch_[(std::size_t) c * (std::size_t) maxBlock_ + (std::size_t) i];

                // The infra-low split and the K-weighting run on EVERY prepared channel at every sample —
                // an absent channel and a hole both feed them the canonical zero — so their state advances
                // on the sample clock and never on the call's.
                float low = 0.0f, high = 0.0f;
                lr4_.processSample (c, x, low, high);
                const double kwY = kw_.process (c, (double) x);

                if (present) { ch.last = raw; ch.lastFinite = finite; }
                if (! present) { ++ch.absent; continue; }
                if (! finite)
                {
                    ++ch.nonFinite; ++nonFiniteInput_;
                    // A hole is not a ZERO, so it ends a digital-silence run. It is not ACTIVE either —
                    // frameAbsMax is taken over the finite present channels only — so it cannot move the
                    // span, which is what makes the exclusion exact.
                    frameNonFinite = true; frameNonZero = true;
                    continue;
                }

                ++ch.finite; ++frameFinite;
                framePresentFinite[c] = true;
                const double d  = (double) x;
                const double sq = d * d;
                frameX[c]   = d;              // held until the frame's span predicate is known
                frameSq[c]  = sq;
                ch.hopSq   += sq;             // the sub-hop's own energy is not a span quantity
                frameSumSq += sq;

                // A finite input is not a finite intermediate: `Svf`'s `(float)(2·v - ic)` update can leave
                // the float range, and a stream of FLT_MAX drives this crossover non-finite at sample 849
                // exactly (measured). Caught where it is
                // produced, counted, and never allowed into an accumulator or the trace.
                const double lo = (double) low;
                if (std::isfinite (lo)) { const double l2 = lo * lo; ch.hopLowSq += l2; frameLowSq[c] = l2; }
                else                    { ++nonFiniteIntermediate_; ch.lowFinite = false; }
                if (std::isfinite (kwY)) ch.kwSq += kwY * kwY;
                else                     { ++nonFiniteIntermediate_; ++kwOverflow_; }

                const double a = std::fabs (d);
                if (a > ch.peak)     ch.peak = a;
                if (a > frameAbsMax) frameAbsMax = a;
                if (! core::exactlyEqual (x, 0.0f)) frameNonZero = true;
            }

            // --- the stereo pair: only frames that carried BOTH channels, finite ---
            if (channels_ >= 2 && nch >= 2)
            {
                const float l = in[0][off + i], r = in[1][off + i];
                if (std::isfinite (l) && std::isfinite (r))
                {
                    sums_.add (l, r);
                    ++stereoFrames_;
                    if (core::exactlyEqual (l, r)) ++stereoValueEqual_;
                    if (std::bit_cast<std::uint32_t> (l) == std::bit_cast<std::uint32_t> (r)) ++stereoBitEqual_;
                }
            }

            // --- THE SPAN, AND THE TWO SILENCE PREDICATES ---
            // THE SPAN IS THE THRESHOLDED REGION, and which predicate defines it is the decision that
            // decides whether this report's padding invariance applies to real deliverables. An exact-zero
            // span reads cleaner and is a CLIFF exactly where masters live: a 24-bit dithered file has no
            // exact zero anywhere, so its five seconds of "silent" lead-in are ±2e-7 of dither, an
            // exact-zero span swallows them, the RMS is diluted by 5.4 dB and a programme CUT mid-fade
            // reads tailEnergyRatio 9.1e-12 — a truncated fade reporting as a perfect one, which is the
            // one failure the span exists to prevent. Measured on this header before the change, on
            // material that is the dominant delivery format. Against the threshold instead: |0| > thr is
            // false for every threshold this class accepts, so DIGITAL padding stays bit-exactly invariant,
            // and sub-threshold padding becomes invariant too. The cost is named rather than hidden — at
            // the -96 dBFS default a 16-bit dither pad (±1 LSB = -90.3 dBFS) is ABOVE the threshold and
            // stays inside the span; a caller who wants it excluded raises the threshold, which is what a
            // named parameter with a documented default is for.
            const bool active = frameAbsMax > silenceThr_;
            if (active)
            {
                if (firstActive_ < 0) firstActive_ = t;
                lastActive_ = t;
                nfAfterLastActive_ = 0;
                for (int c = 0; c < channels_; ++c)
                {
                    ChannelState& ch = chans_[(std::size_t) c];
                    ch.finiteTrail = 0;                           // this frame is INSIDE the span
                    // COMMIT the run that just ended plus this frame: silence inside the programme is part
                    // of it, so its ring-down belongs to the span. Committing is not the same as
                    // subtracting a trailing run from a total, and the difference is not numerical taste —
                    // the crossover's ring-down past the LAST active frame can carry MORE low-band energy
                    // than the span does (a lone impulse is the extreme: one frame in the span, the whole
                    // ring-down outside it), so `total - trail` reaches exactly zero after 224 trailing
                    // zeros and goes NEGATIVE after 3601. A negative energy fraction is not a rounding
                    // error, it is a wrong measurement, and this shape cannot produce one.
                    ch.dc.add (ch.dcPending + frameX[c]);
                    ch.sumSq.add (ch.sumSqPending + frameSq[c]);
                    ch.lowSqSpan.add (ch.lowSqPending + frameLowSq[c]);
                    ch.dcPending = 0.0; ch.sumSqPending = 0.0; ch.lowSqPending = 0.0;
                    if (c < nch)
                    {
                        ch.lastSignal       = in[c][off + i];
                        ch.lastSignalFinite = std::isfinite (ch.lastSignal);
                        ch.hasLastSignal    = true;
                    }
                }
                tailFlushPendingZeros();
                tailWrite (frameSumSq, frameFinite);
            }
            else
            {
                if (frameNonFinite)
                {
                    if (firstActive_ < 0) ++nfBeforeFirstActive_;
                    ++nfAfterLastActive_;
                }
                for (int c = 0; c < channels_; ++c)
                {
                    ChannelState& ch = chans_[(std::size_t) c];
                    // ONLY PAST THE FIRST ACTIVE FRAME. Before it there is no span to commit into, and a
                    // pending run carried across it would put the LEADING quiet region inside the span —
                    // which with an exact-zero span was harmless (those frames contributed exactly 0) and
                    // with a thresholded one is a real leak: measured as a 1e-6 relative disagreement with
                    // the whole-file reference on a 1e-4 tone, whose sub-threshold zero crossings are the
                    // leading run. The oracle found this; no invariance check could have.
                    if (firstActive_ >= 0)
                    {
                        ch.dcPending    += frameX[c];             // committed only if an active frame follows
                        ch.sumSqPending += frameSq[c];
                        ch.lowSqPending += frameLowSq[c];
                    }
                    if (framePresentFinite[c])
                    {
                        ++ch.finiteTrail;
                        if (firstActive_ < 0) ++ch.finiteLead;
                    }
                }
                // A DEFERRED FRAME KEEPS ITS OWN SAMPLE COUNT. Its energy is not necessarily zero — a
                // sub-threshold frame is below the threshold, not silent — and its COUNT is not uniform
                // either, since a channel may be absent or a hole. Writing these back as count 0 would
                // shrink the tail window's denominator and INFLATE the tail mean square, a wrong number
                // rather than a missing one. Only the last `tailSamples_` frames of a run can survive into
                // the window, so a ring of that size holds every count the window can ever need.
                tailPendSq_[(std::size_t) (tailPendingZeros_ % tailSamples_)] = frameSumSq;
                tailPendN_[(std::size_t) (tailPendingZeros_ % tailSamples_)] = (std::int32_t) frameFinite;
                ++tailPendingZeros_;
            }
            if (active != wasActive_) { emitEdge (0, active, t, frameAbsMax); wasActive_ = active; }

            // THE DIGITAL PREDICATE is threshold-free evidence and now drives only its own two counts. A
            // non-finite sample is not a zero, so it ENDS a digital-silence run: a NaN in the last sample
            // of a file reports 0 trailing digital silence, with the hole counted beside it. That is
            // literal rather than helpful, and it is preferred to a count that pretends a hole is silence.
            if (frameNonZero)
            {
                if (firstNonZero_ < 0) firstNonZero_ = t;
                lastNonZero_ = t;
            }
            if (frameNonZero != wasNonZero_) { emitEdge (1, frameNonZero, t, frameAbsMax); wasNonZero_ = frameNonZero; }

            ++totalSamples_;

            // Maintenance rides AUDIO time, never the end of a call. `eq::Crossover2` hands out
            // `flushDenormals()` and keeps no cadence of its own (Crossover2.h:52) — the owner clocks it.
            if (grid_.advance (1)) lr4_.flushDenormals();

            if (totalSamples_ == nextSubHopEnd_) { closeSubHop(); nextSubHopEnd_ += subHopSamples_; }
        }

        // The two sub-meters, at the constant prepared width, over the SAME samples in the same order.
        // Both are bit-invariant under chunking at a fixed width — `ReferenceTruePeakMeter` because its
        // polyphase ring makes every output a function of the samples seen so far and a maximum is
        // associative, `LoudnessMeter` because its sub-hop boundary is its own sample clock.
        const float* planes[core::kMaxChannels] {};
        for (int c = 0; c < channels_; ++c)
            planes[(std::size_t) c] = scratch_.data() + (std::size_t) c * (std::size_t) maxBlock_;
        (void) lm_.process (planes, channels_, m);
        (void) tp_.process (planes, channels_, m);
    }

    //==========================================================================
    // The tail ring holds the last `tailSamples_` frames ENDING AT the last signal frame: trailing digital
    // silence is held back rather than written, so padding a master with zeros cannot push its own tail out
    // of the window. Interior silence IS written — it is part of the programme.
    void tailWrite (double sq, int count) noexcept
    {
        const std::size_t i = (std::size_t) (tailWritten_ % tailSamples_);
        tailSq_[i] = sq;
        tailN_[i]  = count;
        ++tailWritten_;
    }

    void tailFlushPendingZeros() noexcept
    {
        if (tailPendingZeros_ == 0) return;
        const std::int64_t keep = std::min (tailPendingZeros_, tailSamples_);
        tailWritten_ += tailPendingZeros_ - keep;                 // what a full ring would overwrite anyway
        for (std::int64_t j = tailPendingZeros_ - keep; j < tailPendingZeros_; ++j)
            tailWrite (tailPendSq_[(std::size_t) (j % tailSamples_)],
                       tailPendN_[(std::size_t) (j % tailSamples_)]);
        tailPendingZeros_ = 0;
    }

    //==========================================================================
    // A 10 ms sub-hop closes. The unit is BS.1770's, and it is the first energy boundary that does not
    // depend on how the caller chunks; the K-weighting flush lives here for the same reason
    // `LoudnessMeter` puts it here (docs/LAW8-KWEIGHTING.md).
    void closeSubHop() noexcept
    {
        const std::int64_t begin = totalSamples_ - subHopSamples_;
        double ms = 0.0;
        std::int64_t finiteInHop = 0;
        for (int c = 0; c < channels_; ++c)
        {
            ChannelState& ch = chans_[(std::size_t) c];
            emit (ProgrammeTraceKind::SubHopChannel, c, subHopIndex_, begin, totalSamples_,
                  ch.hopSq, ch.hopLowSq, ch.kwSq);
            ms += ch.kwSq / (double) subHopSamples_;               // unit channel weights, as LoudnessMeter's
            finiteInHop += ch.finite;
            ch.hopSq = 0.0; ch.hopLowSq = 0.0; ch.kwSq = 0.0;
        }
        if (! std::isfinite (ms)) { ms = 0.0; ++nonFiniteIntermediate_; ++kwOverflow_; }

        subRing_[(std::size_t) (subHopIndex_ % (std::int64_t) kShortTermSubHops)] = ms;
        ++subHopIndex_;
        kw_.flushDenormals();

        emit (ProgrammeTraceKind::SubHop, -1, subHopIndex_ - 1, begin, totalSamples_,
              ms, (double) finiteInHop, 0.0);

        if (subHopIndex_ == nextObservation_)
        {
            observeShortTerm();
            nextObservation_ += kObservationHops;
        }
    }

    // A short-term (3 s) observation. The window sum is RECOMPUTED from the ring — never carried as
    // `sum += new - old`, which is deterministic but drifts over a long programme (LoudnessMeter.h:320-327
    // recomputes for the same reason).
    void observeShortTerm() noexcept
    {
        double s = 0.0;
        for (std::int64_t j = 0; j < (std::int64_t) kShortTermSubHops; ++j)   // oldest first: one fixed order
        {
            const std::int64_t idx = subHopIndex_ - (std::int64_t) kShortTermSubHops + j;
            s += subRing_[(std::size_t) (idx % (std::int64_t) kShortTermSubHops)];
        }
        const double meanSquare = s / (double) kShortTermSubHops;
        const bool   gated = std::isfinite (meanSquare) && meanSquare >= kAbsoluteGateEnergy;

        if (shortTermCount_ < (std::int64_t) shortTerm_.size())
            shortTerm_[(std::size_t) shortTermCount_++] = std::isfinite (meanSquare) ? meanSquare : 0.0;
        else
            ++shortTermDropped_;                                  // capacity is data, never control

        const std::int64_t k = observationIndex_++;
        emit (ProgrammeTraceKind::ShortTerm, -1, k,
              totalSamples_ - (std::int64_t) kShortTermSubHops * subHopSamples_, totalSamples_,
              meanSquare, (double) binOf (meanSquare), gated ? 1.0 : 0.0);
    }

    void emitEdge (int predicate, bool entering, std::int64_t t, double frameMax) noexcept
    {
        ++edges_;
        emit (ProgrammeTraceKind::SilenceEdge, -1, edges_ - 1, t, t + 1,
              (double) predicate, entering ? 1.0 : 0.0, frameMax);
    }

    void emit (ProgrammeTraceKind kind, int channel, std::int64_t index, std::int64_t begin,
               std::int64_t end, double a, double b, double c) noexcept
    {
        if (trace_ == nullptr) return;
        if ((std::size_t) traceCount_ >= traceCap_) { ++traceOverflow_; return; }
        ProgrammeTraceEvent& e = trace_[(std::size_t) traceCount_++];
        e.kind = kind; e.channel = channel; e.index = index;
        e.begin = begin; e.end = end; e.a = a; e.b = b; e.c = c;
    }

    //==========================================================================
    static double lufsOf (double meanSquare) noexcept
    {
        return meanSquare > 1e-12 ? -0.691 + 10.0 * std::log10 (meanSquare) : kSilenceFloorLufs;
    }

    static int binOf (double meanSquare) noexcept
    {
        const double l = lufsOf (meanSquare);
        const int b = (int) ((l - kHistogramFloorLufs) / kHistogramBinLu);
        return b < 0 ? 0 : (b >= kHistogramBins ? kHistogramBins - 1 : b);
    }

    static ProgrammeValue good (double v) noexcept { return ProgrammeValue { v, true, ProgrammeReason::None }; }
    static ProgrammeValue bad (ProgrammeReason r) noexcept { return ProgrammeValue { 0.0, false, r }; }

    //==========================================================================
    void buildReport() noexcept
    {
        Report& R = report_;
        R = Report {};
        R.channels     = channels_;
        R.totalSamples = totalSamples_;
        R.sampleRate   = sampleRate_;

        R.nonFiniteInputSamples  = nonFiniteInput_;
        R.nonFiniteIntermediates = nonFiniteIntermediate_;
        R.loudnessDroppedBlocks  = (std::int64_t) lm_.droppedBlocks();
        R.loudnessNonFiniteSubHops = (std::int64_t) std::min<std::uint64_t> (
            lm_.nonFiniteSubHops(), (std::uint64_t) std::numeric_limits<std::int64_t>::max());
        R.shortTermDroppedObservations = shortTermDropped_;
        R.silenceEdges           = edges_;
        R.uncoveredSubHopSamples = totalSamples_ - subHopIndex_ * subHopSamples_;

        const bool anySignal = lastActive_ >= 0;                    // the span IS the thresholded region
        R.programmeSpanSamples = anySignal ? lastActive_ + 1 - firstActive_ : 0;

        // --- silence ---
        R.leadingDigitalSilenceSamples  = firstNonZero_ < 0 ? totalSamples_ : firstNonZero_;
        R.trailingDigitalSilenceSamples = lastNonZero_  < 0 ? totalSamples_ : totalSamples_ - 1 - lastNonZero_;
        R.leadingSilenceSamples.value   = firstActive_ < 0 ? totalSamples_ : firstActive_;
        R.trailingSilenceSamples.value  = lastActive_  < 0 ? totalSamples_ : totalSamples_ - 1 - lastActive_;
        if (nfBeforeFirstActive_ > 0)
        {
            R.leadingSilenceSamples.valid  = false;
            R.leadingSilenceSamples.reason = ProgrammeReason::NonFiniteInput;
        }
        if (nfAfterLastActive_ > 0)
        {
            R.trailingSilenceSamples.valid  = false;
            R.trailingSilenceSamples.reason = ProgrammeReason::NonFiniteInput;
        }

        // --- per channel, over the programme span ---
        const bool inputClean = nonFiniteInput_ == 0;
        Sum totalSq, totalLow;
        std::int64_t totalCount = 0;
        double totalPeak = 0.0;
        bool   allLowFinite = inputClean;

        for (int c = 0; c < channels_; ++c)
        {
            const ChannelState& ch = chans_[(std::size_t) c];
            PerChannel& P = R.channel[(std::size_t) c];
            P.finiteSamples    = ch.finite;
            P.nonFiniteSamples = ch.nonFinite;
            P.absentSamples    = ch.absent;
            if (! ch.lowFinite) allLowFinite = false;

            P.lastSample = ch.lastFinite
                         ? good ((double) ch.last)
                         : bad (totalSamples_ == 0 ? ProgrammeReason::NoFiniteSamples
                                                   : ProgrammeReason::NonFiniteInput);
            P.lastSignalSample = (ch.hasLastSignal && ch.lastSignalFinite)
                               ? good ((double) ch.lastSignal)
                               : bad (ch.hasLastSignal ? ProgrammeReason::NonFiniteInput
                                                       : ProgrammeReason::SilentProgramme);

            const std::int64_t spanCount = ch.finite - ch.finiteLead - ch.finiteTrail;
            // A MAXIMUM IS NOT EXCLUSION-EXACT: the sample dropped for being non-finite may have been the
            // largest one, so with nothing finite to compare there is no peak — and 0.0 published as valid
            // would be the comfortable zero this class forbids.
            P.samplePeak = ch.finite > 0 ? good (ch.peak)
                                         : bad (ch.nonFinite > 0 ? ProgrammeReason::NoFiniteSamples
                                                                 : ProgrammeReason::SilentProgramme);
            if (ch.peak > totalPeak) totalPeak = ch.peak;

            if (! anySignal || spanCount <= 0)
            {
                const ProgrammeReason r = anySignal ? ProgrammeReason::NoFiniteSamples
                                                    : ProgrammeReason::SilentProgramme;
                P.dcOffset = P.rms = P.crestFactorDb = P.infraLowFraction = bad (r);
                continue;
            }

            const double sq     = ch.sumSq.value();
            const double lowS   = ch.lowSqSpan.value();
            const double meanSq = sq / (double) spanCount;
            const double rms    = std::sqrt (meanSq);
            P.dcOffset      = good (ch.dc.value() / (double) spanCount);
            P.rms           = good (rms);
            P.crestFactorDb = rms > 0.0 ? good (20.0 * std::log10 (ch.peak / rms))
                                        : bad (ProgrammeReason::SilentProgramme);
            P.infraLowFraction = ! inputClean      ? bad (ProgrammeReason::NonFiniteInput)
                               : ! ch.lowFinite    ? bad (ProgrammeReason::NonFiniteIntermediate)
                               : ! (sq > 0.0)      ? bad (ProgrammeReason::SilentProgramme)
                               :                     good (lowS / sq);

            totalSq.add (sq);
            if (ch.lowFinite) totalLow.add (lowS);
            totalCount += spanCount;
        }

        R.programmeSpanFinite = totalCount;
        std::int64_t finiteAll = 0;
        for (int c = 0; c < channels_; ++c) finiteAll += chans_[(std::size_t) c].finite;
        R.samplePeak = finiteAll > 0 ? good (totalPeak)
                                     : bad (nonFiniteInput_ > 0 ? ProgrammeReason::NoFiniteSamples
                                                                : ProgrammeReason::SilentProgramme);

        const double progSq = totalSq.value();
        if (anySignal && totalCount > 0)
        {
            const double meanSq = progSq / (double) totalCount;
            const double rms    = std::sqrt (meanSq);
            R.programmeMeanSquare = good (meanSq);
            R.rms                 = good (rms);
            R.crestFactorDb = rms > 0.0 ? good (20.0 * std::log10 (totalPeak / rms))
                                        : bad (ProgrammeReason::SilentProgramme);
            R.infraLowFraction = ! inputClean   ? bad (ProgrammeReason::NonFiniteInput)
                               : ! allLowFinite ? bad (ProgrammeReason::NonFiniteIntermediate)
                               : ! (progSq > 0.0) ? bad (ProgrammeReason::SilentProgramme)
                               :                    good (totalLow.value() / progSq);
        }
        else
        {
            const ProgrammeReason r = anySignal ? ProgrammeReason::NoFiniteSamples
                                                : ProgrammeReason::SilentProgramme;
            R.programmeMeanSquare = R.rms = R.crestFactorDb = R.infraLowFraction = bad (r);
        }

        buildTail (R, progSq, totalCount);
        buildStereo (R);
        buildDynamics (R);
    }

    void buildTail (Report& R, double progSq, std::int64_t progCount) noexcept
    {
        R.tailWindowSamples = tailSamples_;
        if (lastActive_ < 0) { R.tailEnergyRatio = bad (ProgrammeReason::SilentProgramme); return; }

        const std::int64_t have = std::min (tailWritten_, tailSamples_);
        Sum tail; std::int64_t count = 0;
        for (std::int64_t k = tailWritten_ - have; k < tailWritten_; ++k)   // chronological, one fixed order
        {
            const std::size_t i = (std::size_t) (k % tailSamples_);
            tail.add (tailSq_[i]);
            count += (std::int64_t) tailN_[i];
        }
        R.tailWindowFinite = count;
        emit (ProgrammeTraceKind::TailWindow, -1, 0, tailWritten_ - have, tailWritten_,
              tail.value(), (double) count, progCount > 0 ? progSq / (double) progCount : 0.0);

        // THE SPAN, NOT THE FILE POSITION. `tailWritten_` is `lastActive_ + 1`, which counts the LEADING
        // digital silence too, so a file with 5 s of zeros and 10 ms of music would have passed this test
        // and then measured a window that is 99 % leading silence — frames that are outside the span every
        // other mean in this report is taken over. The window has to fit inside the programme.
        if (R.programmeSpanSamples < tailSamples_) { R.tailEnergyRatio = bad (ProgrammeReason::ShorterThanTailWindow); return; }
        if (count <= 0)                         { R.tailEnergyRatio = bad (ProgrammeReason::NoFiniteSamples); return; }
        if (progCount <= 0 || ! (progSq > 0.0)) { R.tailEnergyRatio = bad (ProgrammeReason::ZeroDenominator); return; }
        R.tailEnergyRatio = good ((tail.value() / (double) count) / (progSq / (double) progCount));
    }

    void buildStereo (Report& R) noexcept
    {
        R.stereoFrames             = stereoFrames_;
        R.stereoBitIdenticalFrames = stereoBitEqual_;
        R.stereoValueEqualFrames   = stereoValueEqual_;

        if (channels_ < 2)
        {
            R.stereoBalanceDb = R.stereoCorrelation = R.stereoSideToMidRatio
                = bad (ProgrammeReason::MonoProgramme);
            return;
        }
        if (stereoFrames_ == 0)
        {
            R.stereoBalanceDb = R.stereoCorrelation = R.stereoSideToMidRatio
                = bad (ProgrammeReason::NoStereoFrames);
            return;
        }
        if (! (std::isfinite (sums_.ll) && std::isfinite (sums_.rr) && std::isfinite (sums_.lr)
               && std::isfinite (sums_.mid) && std::isfinite (sums_.side)))
        {
            R.stereoBalanceDb = R.stereoCorrelation = R.stereoSideToMidRatio
                = bad (ProgrammeReason::NonFiniteIntermediate);
            return;
        }

        // THE PRIMITIVE'S OWN NUMBER, but not its own silence convention: `StereoSums::correlation()`
        // answers +1.0 below its 1e-12 denominator, and "perfectly in phase" is a verdict about a
        // programme that carries no stereo information at all.
        const double denom = std::sqrt (sums_.ll * sums_.rr);
        R.stereoCorrelation = denom > 1e-12 ? good (sums_.correlation())
                                            : bad (ProgrammeReason::SilentProgramme);
        R.stereoBalanceDb   = (sums_.ll > 0.0 && sums_.rr > 0.0)
                            ? good (10.0 * std::log10 (sums_.ll / sums_.rr))
                            : bad (ProgrammeReason::SilentProgramme);
        R.stereoSideToMidRatio = sums_.mid > 0.0 ? good (sums_.side / sums_.mid)
                                                 : bad (ProgrammeReason::SilentProgramme);
    }

    void buildDynamics (Report& R) noexcept
    {
        R.shortTermObservations = shortTermCount_ + shortTermDropped_;

        // What damages the recursive family, in the order the reasons should be read.
        const ProgrammeReason damage = nonFiniteInput_ > 0          ? ProgrammeReason::NonFiniteInput
                                     : lm_.nonFiniteSubHops() != 0  ? ProgrammeReason::NonFiniteIntermediate
                                     : kwOverflow_ > 0              ? ProgrammeReason::NonFiniteIntermediate
                                     :                                ProgrammeReason::None;
        const bool damaged = damage != ProgrammeReason::None;

        // --- true peak, from the REFERENCE meter: 4× oversampling, 32 taps per phase, drained ---
        const double tpLin = tp_.truePeakLinear();
        const ProgrammeValue tpDb = damaged                    ? bad (damage)
                                  : ! std::isfinite (tpLin)    ? bad (ProgrammeReason::NonFiniteIntermediate)
                                  // gainToDb FLOORS at kGainToDbFloor, so a peak at or under 1e-12 would
                                  // publish -240.0 dBTP as a valid reading rather than a floored one.
                                  : ! (tpLin > core::kGainToDbFloor) ? bad (ProgrammeReason::SilentProgramme)
                                  :                              good (core::gainToDb (tpLin));
        R.truePeakDbtp = tpDb;

        // --- integrated loudness: -120.0 is the meter's "nothing passed the gates" SENTINEL ---
        const double lufs = lm_.integratedLufs();
        const ProgrammeValue integrated =
              damaged                                      ? bad (damage)
            : lm_.droppedBlocks() != 0                     ? bad (ProgrammeReason::LoudnessCapacityExceeded)
            : ! std::isfinite (lufs)                       ? bad (ProgrammeReason::NonFiniteIntermediate)
            // -120.0 is the meter's sentinel for "nothing was averaged", and it has two causes that are
            // not the same answer: no gating block EXISTED (the programme is shorter than one 400 ms block)
            // or blocks existed and the gates rejected them all. `gatingBlockCount()` is public, so the two
            // are separable without touching that header.
            : lm_.gatingBlockCount() == 0                   ? bad (ProgrammeReason::ShorterThanLoudnessWindow)
            : core::exactlyEqual (lufs, kSilenceFloorLufs)  ? bad (ProgrammeReason::SilentProgramme)
            :                                                good (lufs);
        R.integratedLufs = integrated;

        R.plrDb = (tpDb.valid && integrated.valid) ? good (tpDb.value - integrated.value)
                : ! tpDb.valid                     ? bad (tpDb.reason)
                :                                    bad (integrated.reason);

        buildShortTerm (R, damaged, damage);
    }

    // The short-term distribution and EBU Tech 3342's LRA, both from the SAME stored series, so their
    // validity is exact rather than inferred from the programme's duration.
    void buildShortTerm (Report& R, bool damaged, ProgrammeReason damage) noexcept
    {
        const auto all = [&] (ProgrammeValue v)
        {
            R.lraLu = R.shortTermP10 = R.shortTermP50 = R.shortTermP95 = R.shortTermSpreadLu = v;
        };
        if (damaged)              { all (bad (damage)); return; }
        if (shortTermDropped_ > 0) { all (bad (ProgrammeReason::LoudnessCapacityExceeded)); return; }

        // Pass one: the absolute gate, and the mean the relative gate is taken from.
        int hist[kHistogramBins] = { 0 };
        double absSum = 0.0;
        int    absCount = 0;
        for (std::int64_t j = 0; j < shortTermCount_; ++j)
        {
            const double e = shortTerm_[(std::size_t) j];
            if (std::isfinite (e) && e >= kAbsoluteGateEnergy)
            {
                absSum += e; ++absCount;
                ++hist[(std::size_t) binOf (e)];
            }
        }
        R.shortTermGatedObservations = absCount;

        if (absCount < 2)
        {
            all (bad (shortTermCount_ < 2 ? ProgrammeReason::ShorterThanLoudnessWindow
                                          : ProgrammeReason::TooFewObservations));
            return;
        }
        R.shortTermP10 = good (percentile (hist, absCount, 0.10));
        R.shortTermP50 = good (percentile (hist, absCount, 0.50));
        R.shortTermP95 = good (percentile (hist, absCount, 0.95));
        R.shortTermSpreadLu = good (R.shortTermP95.value - R.shortTermP10.value);

        // Pass two: EBU Tech 3342 — the -20 LU relative gate, then P95 - P10 over what is left.
        const double relT = kRelativeGate * (absSum / (double) absCount);
        for (int b = 0; b < kHistogramBins; ++b) hist[b] = 0;
        int lraCount = 0;
        for (std::int64_t j = 0; j < shortTermCount_; ++j)
        {
            const double e = shortTerm_[(std::size_t) j];
            if (std::isfinite (e) && e >= kAbsoluteGateEnergy && e >= relT)
            {
                ++hist[(std::size_t) binOf (e)]; ++lraCount;
            }
        }
        R.lraLu = lraCount >= 2
                ? good (percentile (hist, lraCount, 0.95) - percentile (hist, lraCount, 0.10))
                : bad (ProgrammeReason::TooFewObservations);
    }

    // libebur128's rank rule over a fixed-edge histogram, returning the bin's lower bound — the same shape
    // `LoudnessMeter::lra()` uses, so the two answers are comparable. The fixed edges are what make the
    // percentile allocation-free; they are NOT what makes it re-split invariant — the SERIES being
    // identical is. Quantisation: each percentile is within one 0.1 LU bin, so a difference of two of them
    // is within 0.1 LU, and a distribution with few observations is coarse for the rank rule's own reason
    // (adjacent order statistics can be several LU apart) — hence the published observation counts.
    static double percentile (const int* hist, int total, double p) noexcept
    {
        const int rank = (int) ((double) (total - 1) * p + 0.5);
        int cum = 0, b = 0;
        for (; b < kHistogramBins; ++b) { cum += hist[b]; if (cum > rank) break; }
        if (b >= kHistogramBins) b = kHistogramBins - 1;
        return kHistogramFloorLufs + (double) b * kHistogramBinLu;
    }

    //==========================================================================
    ProgrammeReportParams params_ {};
    Report                report_ {};

    double       sampleRate_    = 48000.0;
    double       silenceThr_    = 1.5848931924611136e-05;      // 10^(-96/20)
    int          maxBlock_      = 0;
    int          channels_      = 0;
    std::int64_t tailSamples_   = 0;
    std::int64_t subHopSamples_ = 0;

    std::array<ChannelState, core::kMaxChannels> chans_ {};

    std::vector<double>       tailSq_;
    std::vector<std::int32_t> tailN_;
    std::vector<std::int32_t> tailPendN_;
    std::vector<double>       tailPendSq_;
    std::vector<float>        scratch_;
    std::vector<double>       shortTerm_;
    std::array<double, (std::size_t) kShortTermSubHops> subRing_ {};

    std::int64_t totalSamples_ = 0, nextSubHopEnd_ = 0;
    // No "filled" counter: the first observation is at sub-hop kShortTermSubHops, by which point the ring
    // holds exactly that many written sub-hops, so there is nothing for such a counter to guard. One was
    // here; a crew round mutated its comparison and the suite stayed green, because no audio input can
    // reach state nothing reads.
    std::int64_t subHopIndex_ = 0, nextObservation_ = 0, observationIndex_ = 0;
    std::int64_t tailWritten_ = 0, tailPendingZeros_ = 0;
    std::int64_t firstActive_ = -1, lastActive_ = -1, firstNonZero_ = -1, lastNonZero_ = -1;
    std::int64_t nfBeforeFirstActive_ = 0, nfAfterLastActive_ = 0, edges_ = 0;
    bool         wasActive_ = false, wasNonZero_ = false;
    std::int64_t stereoFrames_ = 0, stereoBitEqual_ = 0, stereoValueEqual_ = 0;
    std::int64_t shortTermCount_ = 0, shortTermDropped_ = 0;
    // Two overflow counters, not one: the crossover's overflow damages the infra-low fraction and the
    // K-weighting's damages the loudness family, and charging either to the other publishes a refusal whose
    // stated cause is a filter the field does not go through.
    std::int64_t nonFiniteInput_ = 0, nonFiniteIntermediate_ = 0, kwOverflow_ = 0;

    StereoSums             sums_ {};
    KWeightingFilter       kw_ {};
    eq::Crossover2         lr4_ {};
    core::StateGrid        grid_ {};
    LoudnessMeter          lm_ {};
    ReferenceTruePeakMeter tp_ {};

    ProgrammeTraceEvent* trace_      = nullptr;
    std::size_t          traceCap_   = 0;
    std::int64_t         traceCount_ = 0, traceOverflow_ = 0;

    bool prepared_ = false, finished_ = false;
};

} // namespace felitronics::analysis

// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/analysis/SpectrumFrames.h>
#include <felitronics/core/DetMath.h>
#include <felitronics/core/Config.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace felitronics::analysis
{

//==============================================================================
// felitronics::analysis::HumDetector — MAINS HUM: a narrow, STATIONARY line at 50 or 60 Hz with a comb of exact
// multiples, looked for only where the programme is QUIET, because music masks it.
//
// WHAT IS MEASURED, AND WHAT IS NOT CLAIMED. The report carries coordinates and continuous evidence: where the
// line sits (in Hz, interpolated inside the bin), how far it stands above the LOCAL background (in dB), which
// harmonics were seen, and which stretches of the programme were used (in samples). `mains` says a
// MAINS-COMPATIBLE STATIONARY LINE was measured at that nominal — it is not a causal claim about electricity.
// A synthesised 50.000 Hz pedal with exact harmonics can be sample-for-sample what an interference pickup
// leaves, and no audio-only instrument separates the two; ClipDetector states the same limit for a hard-driven
// clamp against a square wave ("one PCM, two histories", ClipDetector.h:62). What the instrument CAN do is
// refuse to guess, and say why.
//
// THE FOUR THINGS THAT SEPARATE HUM FROM A BASS NOTE — and the hostile notes are close: G1 = 49.0 Hz sits 1.0 Hz
// from 50, A#1 = 58.3 Hz and B1 = 61.7 Hz sit 1.7 Hz from 60.
//   1. HUM IS THERE IN THE QUIET. Frames are selected by their own loudness (below `quietThresholdDb`), so the
//      evidence comes from where the programme is not.
//   2. THE POSITION, MEASURED INSIDE THE BIN. A peak is accepted only if its INTERPOLATED position is within
//      `toleranceHz` (0.5 Hz) of 50 or 60. 0.5 Hz is what mains asks for and what the notes cannot have: a stiff
//      grid holds +-0.05 Hz and excurses to +-0.2 Hz, while the nearest note is 1.0 Hz away. It is also what
//      makes the resolution below non-negotiable — a bin centre would read "50.000" for anything in its bin.
//   3. IT STANDS STILL. The same line must be accepted in at least TWO quiet stretches, and the spread of every
//      position estimate — per stretch AND per frame — must stay inside `toleranceHz`. A note moves with the
//      music; one that merely passes through 50 Hz in one stretch gives one observation, and one observation is
//      not an answer (see SingleQuietStretch / the stationarity gate).
//   4. THE COMB IS REPORTED, NOT REQUIRED. Exact multiples corroborate, but a bass guitar's partials are also
//      near-exact multiples, so the comb is weak as a discriminator, and a well-filtered pickup can be
//      fundamental-only. Gating on a harmonic count would therefore buy little and cost a FALSE CLEAN, which is
//      the one outcome this instrument may not produce. `harmonicsObserved` is published so a consumer that
//      wants corroboration can ask for it. (Two crew seats argued for a minimum harmonic count; overruled for
//      that reason, and the objection is recorded in the task report.)
//
// WHY EVERY PEAK IN THE WINDOW IS ENUMERATED, NOT JUST THE STRONGEST. This was the design's one real defect and
// two crew seats found it independently. Reading only the argmax of the search window means a LOUDER note 1 Hz
// away steals the window, fails the tolerance, and the real hum beneath it is never examined — a false clean,
// with the resolution paid for and unused. So `locate()` enumerates every strict interior local maximum in the
// window, accepts those inside the tolerance, and keeps the strongest of them; the strongest local maximum
// overall is published separately as `windowPeak`, accepted or not, so contradictory evidence is never censored.
// This costs nothing in false positives: for ONE Hann-windowed tone the sampled skirt is MONOTONE away from the
// main lobe — sampling W(m+f) on the bin grid leaves the |sin(pi*f)| factor constant, so the magnitude falls as
// 1/|m+f|^3 with no interior maximum at all. (The "sidelobe maxima at 2.36 bins" of the CONTINUOUS kernel, which
// one crew seat warned would fake a 50 Hz peak from a 49 Hz note, do not survive sampling: the bin at +2.36 is
// below the bin at +1.36. Verified analytically and by the suite's monotone-skirt witness.) Two tones CAN
// interfere into an extra maximum; the prominence gate and the stationarity gate are what answer that.
//
// RESOLUTION IS A DURATION. To SEPARATE 49.0 from 50.0 Hz — not merely to notice something there — two Hann main
// lobes must sit >= 2 bins apart: at exactly 2 bins the sum of the two kernels is flat across the pair
// (W(0) + W(2) = 0.5 + 0 = 0.5 = 2*W(1) = 2*0.25), and closer than that the dip is gone. 1.0 Hz / 2 bins gives
// bin <= 0.5 Hz, i.e. N >= 96000, i.e. N = 2^17 at 48 kHz (bin 0.366 Hz, window 2.73 s, and 2.73 bins of
// separation for a 6 dB dip). At 2^16 the pair is 1.37 bins apart and merges. Hence `fftOrder = 0` = AUTO: the
// smallest order whose bin is <= 0.5 Hz at this rate. An explicitly coarser order is honoured and answers
// `valid = false, InsufficientResolution` — never a guessed 50. Zero-padding would sharpen the POSITION
// estimate and would not separate the pair; it is not used, and the two must not be confused.
//
// THE QUIET GATE IS THE PROGRAMME'S LOUDNESS, NOT THE HUM'S. A frame's Hann-weighted mean square comes out of
// its own bins by Parseval — for this producer's calibration, P[0] + P[N/2] + 2*sum_{k=1}^{N/2-1} P[k] is
// exactly sum(x^2 w^2) / sum(w^2) (SpectrumFramesTests.cpp:169 pins it) — with two subtractions:
//   * bins 0 and 1 are dropped. Under a PERIODIC Hann a constant d has energy at exactly k = 0, 1, N-1
//     (w = 0.5 - 0.25 e^{+j} - 0.25 e^{-j}), with P0 = 2d^2/3 and 2*P1 = d^2/3 summing to d^2, so in exact
//     arithmetic this removes a DC offset completely — a converter's DC must not read as loudness. A double
//     window through a double FFT leaves a rounding residue, so the removal is exact in theory, not bitwise.
//     Genuine sub-bin programme content goes with it, and a DRIFTING offset or a ramp does not.
//   * the candidate bands are dropped: +-`searchHz` around every h*50 and h*60 for h = 1..`maxHarmonic`, merged.
//     Without this the defect censors its own detection — a -45 dBFS hum makes every frame "not quiet" and the
//     answer becomes NoQuietStretch on the very files that have the problem (found by a crew seat). 16 windows
//     of 6 Hz out of the whole band change a programme measure by nothing; they change this one by everything.
// A frame carrying ANY hole is never quiet: a holed frame's power is deliberately zeroed by the producer
// (SpectrumFrames.h:231) and reading that as silence would manufacture quiet stretches out of dropouts.
// LIMIT, NAMED: the Hann window weights a frame's own edges to zero, so a transient at a frame edge is invisible
// to THAT frame. At the default hop = N/2 every sample sits at the peak weight of some frame, so the transient
// makes a neighbouring frame loud and breaks the run — except inside the first hop of the stream and past the
// last full frame, where fewer frames cover a sample.
//
// A QUIET STRETCH is a maximal run of consecutive SELECTED FRAME INDICES, published as
// [first frame start, last frame end) in samples. Runs separated by any unselected frame are separate stretches
// and are never glued. The interval is a FRAME-SUPPORT ENVELOPE, not a proof that every sample in it was quiet,
// and at hop < N/2 two distinct runs can report overlapping intervals — read `frames` beside the coordinates.
//
// STATIONARITY NEEDS TWO STRETCHES, AND ONE IS ITS OWN ANSWER. "Stands still between stretches" cannot be tested
// inside one, so exactly one usable stretch is `valid = false, SingleQuietStretch` — a named incompleteness, not
// a quietly lowered requirement, and not "no hum". The provisional candidate is still published (it is evidence,
// and the danger here runs one way only), together with `maxIntraStretchSpreadHz` from the individual FRAME
// estimates inside that stretch, so a consumer can see how still the line was without the instrument certifying
// it. Considered and rejected: cutting one long run into artificial "stretches" — that would answer the
// requirement by redefining it.
//
// A RATIO IS NOT EVIDENCE ON ITS OWN, so a line must also be LOUD ENOUGH TO BE A LINE. `minLevelDbfs`
// (-100 dBFS on the tone's mean square) is the second half of the acceptance test, and it is what keeps the
// instrument out of its own arithmetic: a programme that is one pure tone, or one constant sample value,
// has a spectrum of 1e-23 to 1e-35 away from that tone, and a local maximum of THAT residue stands 18 to
// 37 dB above the residue beside it. Measured: a pure 440 Hz tone reported 60 Hz mains at -230 dBFS, a
// pure 700 Hz tone reported 50 Hz at -228 dBFS, and a constant +1.0 reported 50 Hz at -343 dBFS — all with
// two eligible stretches, because a pure tone's frames are identical and so is their residue. -100 dBFS is
// under the noise floor of a 16-bit delivery and 130 decades above that residue.
//
// FOUR LIMITS THE ADVERSARIAL ROUND MEASURED AND THAT NO THRESHOLD CAN REMOVE. They are named here rather
// than papered over, and in every one the report publishes the evidence a consumer needs:
//   * A NOTE SHARP OF G1 ENTERS THE BAND. G1 is 1.0 Hz from 50 at A = 440, but 20 cents sharp it is 49.57,
//     0.43 Hz away — inside the tolerance. A band tuned that sharp, holding a dead-steady drone through two
//     quiet passages, is reported. What saves the realistic version is the stationarity gate: a free vibrato
//     of +-1 Hz at 0.5 Hz was refused (one stretch observation, frame spread 0.07 Hz).
//   * A STRONG NEIGHBOUR BIASES THE POSITION by up to ~0.07 bins, so a line within ~0.03 Hz of the
//     tolerance edge can fall on either side of it: 49.000 + 50.510 Hz reads as 50.484 and is accepted.
//     Every threshold has an edge; this one's width is 8 % of a bin.
//   * A REAL 50 OR 60 Hz LINE OF MUSICAL ORIGIN IS STILL A 50 OR 60 Hz LINE. A 20 Hz tone's third harmonic
//     is exactly 60 Hz; a 100 Hz tone tremoloed at 50 Hz literally contains 50 Hz. Both are reported, and
//     both are correct under what `mains` claims. `bandPowerSum()` is where a consumer sees the 20 Hz
//     fundamental that explains it.
//   * A SPEED-SHIFTED TRANSFER MOVES THE LINE OUT OF THE BAND. A -1.2 % tape transfer puts 50 Hz at 49.4,
//     0.6 Hz out, and no tolerance can reach it without reaching G1. The line is published anyway, as the
//     candidate's `windowPeak` (position and prominence, accepted or not) and counted in
//     `stretchOffTolerance` — which is exactly what those two fields are for.
// And one that is the seam's, not this instrument's: hum living only in the UNCOVERED TAIL (less than one
// window at the end of the programme) is not analysed at all. `tailUncoveredSamples` names how much.
//
// ONE PERIODOGRAM IS NOT A LINE, AND THAT IS WHY A STRETCH MUST HOLD TWO FRAMES. A single frame's
// periodogram has an exponential tail: the chance that some bin of a search window stands 10 dB over the
// paired median of its neighbours is small but not small enough, and if the two quiet stretches carry the
// SAME noise — a looped passage, a repeated room tone, a silent section duplicated by an editor — the
// excursion repeats exactly and passes the stationarity test with it. MEASURED, 1200 noise-only files per
// cell at 4 kHz (`minFramesPerObservation` = 1, i.e. the rule off), reported as mains:
//                                        10 dB   12 dB   14 dB   16 dB
//   one frame per stretch, independent        0       0       0       0
//   one frame per stretch, IDENTICAL noise   75       5       0       0
//   one frame per stretch, 16-bit dither     68       9       0       0
//   three frames per stretch, either way      0       0       0       0
// Two frames per stretch takes every one of those cells to 0 at every gate, and costs nothing: the weakest
// hum found is unchanged (2.2e-5 amplitude at 10 dB) for stretches that have two frames. Raising the gate
// instead would work — 14 dB clears the table too — but it would cost 4 dB of sensitivity on every file and
// leave the cause in place. So the gate stays at 10 dB and a stretch must show the line TWICE.
// The price is named rather than paid silently: fewer than two eligible stretches is
// `valid = false, StretchesTooShort`, never "no hum". At 48 kHz two frames is N + hop = 4.1 s of quiet.
// A LOOPED passage whose period puts a line inside the tolerance is NOT a false positive and is not treated
// as one: a buffer repeated every 997 samples at 4 kHz really does carry a stationary line at 60.18 Hz with
// an exact harmonic comb, and no spectral instrument can call that anything else. The header's first
// paragraph is the answer to it: `mains` is a mains-COMPATIBLE line, not a causal claim.
//
// AVERAGING FRAMES DOES NOT MAKE A LINE EASIER TO SEE, which is worth knowing before tuning anything here:
// it raises the line and its background together, so the peak-to-floor ratio barely moves with frame count.
// What averaging buys is the stability of the estimate and the accuracy of its position. SENSITIVITY comes
// from the WINDOW: a longer one narrows the bin, which lowers the background in that bin while the line
// keeps its power. Measured consequence — a stretch's observation and its frames' observations appear and
// disappear at the same amplitude, at 3 frames per stretch and at 16 alike.
//
// THE LOCAL BACKGROUND, AND WHY IT IS PAIRED. A quiet passage is usually falling LF rumble, and a plain median
// over +-`floorSpanHz` excluding +-`floorExcludeHz` then lands on the largest bin of the QUIETER side — about
// -0.3 dB on a 1/f^2 slope, and ~3.7 dB on a steep one, straight into the prominence. So the floor is the lower
// median of PAIRED geometric means sqrt(P[kp-d])*sqrt(P[kp+d]) over symmetric offsets d: for a locally
// log-linear background every pair estimates the level AT the peak, and the slope cancels instead of biasing.
// Two square roots per pair, not one over the product, so a tiny power cannot underflow to zero. `prominenceDb`
// is 10*log10(tonePower/floorPower), taken once at the end; the GATE is `tone >= floor * 10^(min/10)` with the
// factor built in prepare(), so no libm call decides an acceptance (ClipDetector.h:53's rule).
//
// `tonePower` IS THE TONE'S MEAN SQUARE, not the peak bin. This calibration puts A^2/6 in the peak bin of a
// bin-centred sine of amplitude A; the three-bin folded sum 2*(P[kp-1] + P[kp] + P[kp+1]) is A^2/2 exactly,
// which is the number a consumer means by "the level of the hum", so that is what is published.
// BUT THE PROMINENCE IS MEASURED IN ONE BIN'S BANDWIDTH — `peakBinPower / floorPower`, not tonePower over it.
// The two must share a bandwidth or the ratio carries a constant that has nothing to do with a tone: six bins
// against one is 7.8 dB of head start, and the suite caught exactly that — pure noise cleared a 10 dB gate at
// 119.9 Hz in both quiet stretches of a 49.0 Hz BASS fixture and was reported as 60 Hz mains, base h = 2.
//
// THE COMB, IN TWO PASSES. Pass A fixes the BASE from h = 1 (tolerance 0.5 Hz) or else h = 2 (tolerance 1.0 Hz —
// h*toleranceHz, because a real grid's h-th harmonic is h times as far off nominal): the first that yields an
// accepted, prominent line gives f0est = position / h. Pass B then searches every h in 1..`maxHarmonic` around
// h*f0est within `harmonicToleranceBins` bins, which is where the accuracy is: f0est carries ~0.02 bin of
// estimator error, so h*that is under 0.2 bin at h = 8, inside a 2-bin window. 100/120 Hz with NO fundamental —
// full-wave rectification — therefore reports base 2, `fundamentalObserved = false`, and a fundamental marked
// DERIVED; the unobserved line is never asserted as seen.
// LIMIT, NAMED: the base is looked for at h = 1 and h = 2 only, which is the brief's scope (the fundamental,
// and full-wave rectification's 100/120 Hz). A three-pulse rectifier's first ripple is 150/180 Hz and a
// three-phase six-pulse bridge's is 300/360 Hz, and 300 Hz is both 6*50 and 5*60 — which mains family produced
// it cannot be decided from audio without inventing the answer. So no `mains` is claimed there, but it is NOT
// reported as an absence either: when no base is found the harmonics are still scanned on the NOMINAL grid
// (tolerance h*toleranceHz), and >= 2 of them agreeing on one fundamental is `valid = false, CombWithoutBase`
// with `combFundamentalHz` and `lowestHarmonicObserved` naming what was seen. ONE prominent line at a multiple
// of 50 or 60 is not a comb and does not trigger it.
// LIMIT, NAMED: the comb is measured on the pooled average of every selected frame, so its components need not
// have CO-OCCURRED — two stretches of a bare 50 Hz line and one of a bare 150 Hz line pool into a comb. Per
// stretch, only the four fixed hypotheses (50 and 60 Hz at h = 1 and h = 2) are re-measured, and their per-
// stretch acceptance IS published, so co-occurrence is visible for the base and pooled for the rest.
//
// LAW 8a — bit-identical under ARBITRARY re-slicing into process() calls, and identical for any `maxBlock` at
// prepare(). What that costs in discipline: one integer clock (sample outside, channel inside, as
// ClipDetector.h:200); every completed frame consumed the instant tick() returns true, because the producer's
// power() is only readable until the next frame closes (SpectrumFrames.h:164) and reading once per process()
// would keep only the last one and make the answer a function of the block size; nothing closed, flushed or
// reset at a process() boundary; accumulators in double with one fixed recurrence each, bins ascending, frames
// in absolute order; ratios from accumulated numerator and denominator, never a mean of dB; one division by the
// integer frame count at the end. There is no IIR state and no per-sample recurrence in this instrument, so
// core::StateGrid has nothing to schedule here — the mutant that stands in for "maintenance moved to the call
// boundary" is the active stretch being closed at the end of process(), and the suite is red on it.
//
// CAPACITY IS DATA, NOT CONTROL (law 11d). `maxStretches` bounds the reportable stretch LIST and nothing else:
// the 257th stretch gets its absolute index, contributes every frame to every accumulator, is measured through
// the same reusable active buffer and updates every candidate summary — only its coordinate record is not
// stored, and `stretchesComplete` says so. A capacity of 256 and one of 512 must classify identically. The same
// holds for the frame trace.
//
// FORM. setParams / prepare / process / finish / reset, with Storage + storageFor() published before anything is
// allocated. process() is read-only and allocates nothing; the report is final after finish(), which is
// idempotent and closes the open stretch at its last selected frame's end — it invents no tail frame, because
// the producer emits none (SpectrumFrames.h:180). Every accessor is total: an index outside the report answers a
// zeroed value rather than reading past a vector.
//
// RT: OFFLINE (message thread). A 2^17 frame per 1.365 s of programme is a second or two for a five-minute file.
struct HumDetectorParams
{
    int    fftOrder = 0;                  // 0 = AUTO: the smallest order whose bin is <= kMaxBinHz at this rate
    int    hop = 0;                       // 0 = N/2
    double quietThresholdDb = -50.0;      // a frame is quiet below this (dB, mean square, candidate bands excluded)
    double toleranceHz = 0.5;             // distance from nominal AND the permitted spread of the estimates
    double searchHz = 3.0;                // the neighbourhood local maxima are enumerated in
    double floorSpanHz = 10.0;            // the local background's reach on each side of the peak
    double floorExcludeHz = 1.5;          // ... minus this much around the peak (the Hann main lobe plus margin)
    double minProminenceDb = 10.0;        // a line counts as OBSERVED this far above its local background
    double minLevelDbfs = -100.0;         // ... and only if the line itself is at least this loud (see the header)
    double harmonicToleranceBins = 2.0;   // pass B's window around h*f0est
    int    minFramesPerObservation = 2;   // frames a quiet stretch must hold before its line counts as seen
                                          // (see "ONE PERIODOGRAM IS NOT A LINE" in the header: measured)
    int    maxHarmonic = 8;               // h = 1..maxHarmonic
    int    maxStretches = 256;            // capacity of the reportable stretch list, per channel
    int    traceCapacity = 0;             // frame-trace records (0 = off); a test hook, bounded and published
};

enum class HumMains : std::uint8_t { None = 0, Hz50 = 1, Hz60 = 2 };

enum class HumReason : std::uint8_t
{
    Ok                     = 0,
    NotPrepared            = 1,  // prepare() has not been accepted
    NotFinished            = 2,  // read before finish(); the report is not final
    InsufficientResolution = 3,  // bin > kMaxBinHz: 49.0 and 50.0 Hz cannot be separated at this window
    ShorterThanWindow      = 4,  // the programme is shorter than one window, so there is no frame at all
    AllFramesHoled         = 5,  // every frame carried a non-finite sample or an absent channel
    NoQuietStretch         = 6,  // no frame was quiet — NOT "no hum"
    SingleQuietStretch     = 7,  // exactly one: "stands still between stretches" is untestable by construction
    StretchesTooShort      = 8,  // < 2 stretches hold minFramesPerObservation frames: see the header
    CandidateNotStationary = 9,  // a mains-compatible line WAS measured and did not stand still — not "clean"
    CombWithoutBase        = 10,  // >= 2 consistent harmonics of 50 or 60 Hz, none of them h = 1 or h = 2
};

// One measured spectral line. Evidence, never a verdict.
struct HumPeak
{
    bool   found     = false;   // a strict interior local maximum existed in the search window
    bool   prominent = false;   // ... and it stands minProminenceDb above its local background
    bool   accepted  = false;   // ... and its interpolated position is inside the tolerance
    int    bin       = -1;      // the maximum's bin, -1 when none
    double hz        = 0.0;     // interpolated position (three-bin parabola in dB), 0 when not found
    double tonePower = 0.0;     // linear, 2*(P[kp-1]+P[kp]+P[kp+1]) — the tone's MEAN SQUARE (A^2/2 for a sine)
    double peakBinPower = 0.0;  // linear, P[kp] alone — the one that shares a bandwidth with the background
    double floorPower = 0.0;    // linear, the paired lower median of the local background, 0 when unavailable
    double prominenceDb = 0.0;  // 10*log10(peakBinPower/floorPower), 0 when either is unavailable
};

struct HumHarmonic
{
    int     index  = 0;         // h
    bool    inBand = false;     // h*f0 plus its window fits below Nyquist and inside the accumulated band
    HumPeak peak {};
};

// A quiet stretch, in samples. The interval is the FRAME-SUPPORT ENVELOPE of the run (see the header).
struct HumStretch
{
    std::int64_t index       = 0;   // absolute index, so a truncated list still names what it holds
    std::int64_t startSample = 0;
    std::int64_t endSample   = 0;
    std::int64_t frames      = 0;   // selected frames in the run
};

// One mains hypothesis, measured whether or not it wins.
struct HumCandidate
{
    double  nominalHz = 0.0;
    bool    baseFound = false;
    int     baseHarmonic = 0;            // 1 or 2
    double  fundamentalHz = 0.0;         // base.hz / baseHarmonic
    bool    fundamentalObserved = false; // true only when the base IS h = 1
    HumPeak base {};                     // the base line on the pooled average
    HumPeak windowPeak {};               // the STRONGEST local maximum in the h = 1 window, accepted or not
    int     harmonicsObserved = 0;
    int     lowestHarmonicObserved = 0;  // 0 when none — names a comb whose base is out of this scope
    bool    combWithoutBase = false;     // >= 2 observed harmonics agreeing on a fundamental this scope missed
    double  combFundamentalHz = 0.0;     // ... implied by the lowest of them; NOT an observed line
    std::int64_t stretchObservations = 0;// stretches whose OWN average showed the base, accepted
    std::int64_t stretchOffTolerance = 0;// stretches with a prominent base-window peak OUTSIDE the tolerance
    std::int64_t frameObservations = 0;
    double  stretchSpreadHz = 0.0;       // max - min over stretch estimates
    double  frameSpreadHz = 0.0;         // max - min over every frame estimate, all stretches
    double  maxIntraStretchSpreadHz = 0.0;
    bool    stationary = false;
    bool    passed = false;
};

struct HumReport
{
    bool      valid  = false;
    HumReason reason = HumReason::NotPrepared;
    HumMains  mains  = HumMains::None;   // a mains-COMPATIBLE stationary line, not a causal claim
    HumPeak   line {};                   // the winner's base line; zeroed when ! valid
    double    fundamentalHz = 0.0;       // zeroed when ! valid
    bool      fundamentalObserved = false;
    bool      fundamentalDerived = false;// the fundamental was computed from h = 2, not seen
    int       baseHarmonic = 0;
    int       harmonicsObserved = 0;
    // --- evidence, ALWAYS published, whatever `valid` says ---
    std::int64_t frames = 0, finiteFrames = 0, holedFrames = 0, quietFrames = 0;
    std::int64_t silentFrames = 0;      // frames of EXACT digital silence: excluded, because they show nothing
    std::int64_t quietStretches = 0, eligibleStretches = 0, storedStretches = 0;
    bool         stretchesComplete = true;
    std::int64_t totalSamples = 0;
    std::int64_t tailUncoveredSamples = 0;
    double       binHz = 0.0;
    std::int64_t windowSamples = 0;
};

// The frame TRACE: one record per (frame, channel), filled at FRAME CLOSE. A bounded, published test hook —
// a final-report comparison cannot see a re-slicing defect that cancels in an average or a maximum
// (docs/LAW8-KWEIGHTING.md:78: identical LUFS, 5 of 97 block energies moved), so the suite compares this.
struct HumFrameTrace
{
    std::int64_t frameStart = 0, frameEnd = 0;
    std::int64_t stretchIndex = -1;      // -1 when the frame was not selected
    double       meanSquare = 0.0;       // raw linear, the quiet measure itself
    double       bandSum = 0.0;          // raw linear sum of this frame's accumulated band
    double       hypHz[4] {};            // per hypothesis (50/h1, 50/h2, 60/h1, 60/h2): position, raw
    double       hypTone[4] {};          // raw linear three-bin folded tone power
    double       hypFloor[4] {};         // raw linear local background
    std::int32_t hypBin[4] { -1, -1, -1, -1 };
    std::uint32_t hypAccepted = 0;       // bit h: that hypothesis was accepted in this frame
    std::int32_t channel = 0;
    std::uint8_t finite = 0, quiet = 0;
};

class HumDetector
{
    // Declared FIRST because Storage::bytes() measures them: a budget must be able to name its own rows.
    struct BinRange { int lo = 0, hi = -1; };

    // One fixed hypothesis' running summary: 50 or 60 Hz at h = 1 or h = 2.
    struct HypSummary
    {
        std::int64_t frameObs = 0, stretchObs = 0, stretchOffTol = 0;
        double frameMin = 0.0, frameMax = 0.0;
        double stretchMin = 0.0, stretchMax = 0.0;
        double maxIntra = 0.0;
        std::int64_t frameObsInStretch = 0;                   // reset when a stretch opens
        double frameMinInStretch = 0.0, frameMaxInStretch = 0.0;
    };

public:
    static constexpr double kNominal[2]        = { 50.0, 60.0 };
    static constexpr int    kCandidates        = 2;
    static constexpr int    kBaseHarmonics     = 2;    // the base is looked for at h = 1 and h = 2
    static constexpr int    kHypotheses        = kCandidates * kBaseHarmonics;
    // G1 = 49.0 Hz is the closest hostile note to 50 Hz, so 1.0 Hz is the separation that must be resolved.
    // The requirement is stated in BINS, because 2.0 bins is where two equal Hann main lobes stop having a dip
    // at all and 2.048 bins — what a 0.5 Hz bound admits at 4, 8, 16, 32 and 64 kHz — is measurably not enough:
    // a 49.0 Hz note over a 50.0 Hz hum 14 dB below it is found at 2.73 bins and MISSED at 2.048, and at 2.048
    // an equal-amplitude pair merges into one accepted peak at 49.66 Hz. 2.7 bins is the separation the
    // 48 kHz / 2^17 configuration actually has, and it is the one measured to work (also at 5:1 and 10:1).
    static constexpr double kMinNoteSeparationHz = 1.0;
    static constexpr double kMinSeparationBins = 2.7;
    static constexpr double kMaxBinHz          = kMinNoteSeparationHz / kMinSeparationBins;   // 0.370 Hz
    // ... and nothing is gained below this: a window finer than 0.01 Hz (100 s) measures the grid's own
    // wander rather than a line, and its background sort would dominate the runtime (order 22 at 1 kHz took
    // 114 s for two frames before this bound existed).
    static constexpr double kMinBinHz          = 0.01;
    static constexpr int    kMinHarmonic       = 1;
    static constexpr int    kMaxHarmonicLimit  = 64;
    static constexpr int    kMaxStretchesLimit = 1 << 20;
    static constexpr int    kMaxTraceLimit     = 1 << 20;
    static constexpr int    kMinFloorPairs     = 4;    // fewer than this and the background is unavailable
    static constexpr double kMaxSpanHz         = 100000.0;   // the largest any Hz parameter may be (see geometryFor)
    // Bounded like every other measurer in this module (ClipDetector.h:132): it also keeps binHz far from zero,
    // so every `(int) ceil(hz / binHz)` below is a small number rather than an undefined narrowing of +inf.
    static constexpr double kMinSampleRate     = 1000.0;
    static constexpr double kMaxSampleRate     = 768000.0;

    // The geometry every size and every window is derived from — ONE function, shared by storageFor() and
    // prepare(), so a budget and an allocation cannot drift apart. It is also pinned by hand-computed numbers in
    // the suite: a budget that sizes itself through the same function it is compared with proves nothing.
    struct Geometry
    {
        bool ok = false;
        int  order = 0, bins = 0;
        std::int64_t n = 0, hop = 0;
        double binHz = 0.0;
        int  bandLo = 0, bandHi = 0, bandBins = 0;           // the pooled average's bins
        int  stretchLo = 0, stretchHi = 0, stretchBins = 0;  // the reusable active-stretch buffer's bins
        int  floorSpanBins = 0, floorExcludeBins = 0, floorPairs = 0;
        int  searchBins = 0, spanHarmonics = 0;
        int  excludeCount = 0;                               // merged candidate bands kept out of the quiet gate
    };

    // Every bin index below goes through one of these: the clamp happens in DOUBLE and the conversion
    // afterwards, so a parameter of 1e20 cannot hand `(int)` a value outside its range. Converting first and
    // clamping after — `std::min(last, (int) std::ceil(hz / bin))` — is undefined for exactly those inputs,
    // and UBSan reported all three sites before this existed.
    static int binCeil (double hz, double bin, int lo, int hi) noexcept
    {
        const double v = std::ceil (hz / bin);
        if (! std::isfinite (v)) return hi;
        return (int) std::clamp (v, (double) lo, (double) hi);
    }
    static int binFloor (double hz, double bin, int lo, int hi) noexcept
    {
        const double v = std::floor (hz / bin);
        if (! std::isfinite (v)) return hi;
        return (int) std::clamp (v, (double) lo, (double) hi);
    }

    static Geometry geometryFor (double sampleRate, const HumDetectorParams& p) noexcept
    {
        Geometry g;
        if (! (sampleRate >= kMinSampleRate && sampleRate <= kMaxSampleRate)) return g;   // NaN fails too
        // Each Hz parameter is bounded as well as finite: `bin` can be as small as 1000/2^22, so an
        // unbounded numerator would make a bin index exceed `int` long before it reached infinity.
        if (! (p.toleranceHz > 0.0) || ! (p.toleranceHz <= kMaxSpanHz)) return g;
        if (! (p.searchHz >= p.toleranceHz) || ! (p.searchHz <= kMaxSpanHz)) return g;
        if (! (p.floorExcludeHz > 0.0) || ! (p.floorSpanHz > p.floorExcludeHz)) return g;
        if (! (p.floorSpanHz <= kMaxSpanHz)) return g;
        if (! (p.harmonicToleranceBins > 0.0) || ! (p.harmonicToleranceBins <= 1.0e6)) return g;
        if (! std::isfinite (p.quietThresholdDb) || ! std::isfinite (p.minProminenceDb)) return g;
        if (! std::isfinite (p.minLevelDbfs)) return g;
        if (p.maxHarmonic < kMinHarmonic || p.maxHarmonic > kMaxHarmonicLimit) return g;
        if (p.minFramesPerObservation < 1 || p.minFramesPerObservation > (1 << 20)) return g;
        if (p.maxStretches < 0 || p.maxStretches > kMaxStretchesLimit) return g;
        if (p.traceCapacity < 0 || p.traceCapacity > kMaxTraceLimit) return g;

        int order = p.fftOrder;
        if (order == 0)                                       // AUTO: the shortest window that resolves the pair
        {
            order = 4;
            while (order < 22 && sampleRate / (double) ((std::int64_t) 1 << order) > kMaxBinHz) ++order;
        }
        if (order < 4 || order > 22) return g;
        if (! (sampleRate / (double) ((std::int64_t) 1 << order) >= kMinBinHz)) return g;
        const std::int64_t n = (std::int64_t) 1 << order;
        const std::int64_t hop = p.hop == 0 ? n / 2 : (std::int64_t) p.hop;
        if (hop <= 0 || hop > n) return g;

        g.order = order;
        g.n = n;
        g.hop = hop;
        g.bins = (int) (n / 2 + 1);
        g.binHz = sampleRate / (double) n;
        const double bin = g.binHz;

        // The highest place any search or any background can reach: the largest fundamental the tolerance
        // admits times the harmonic count, plus pass B's window, plus the background's span, plus a bin.
        // The count is at least kBaseHarmonics, NOT `maxHarmonic`: the four fixed hypotheses are 50 and 60 Hz
        // at h = 1 and h = 2 whatever the report asks for, and `maxHarmonic = 1` used to cut the band at
        // 72 Hz — leaving the 100/120 Hz hypotheses unmeasurable and a 120 Hz comb reported as
        // `valid = true, mains = None`, which is the one answer this instrument may not give.
        const int spanHarmonics = std::max (p.maxHarmonic, kBaseHarmonics);
        const double topHz = (kNominal[1] + p.toleranceHz) * (double) spanHarmonics
                           + p.harmonicToleranceBins * bin + p.floorSpanHz + bin;
        const int last = g.bins - 2;                          // kp+1 must always be readable
        if (last < 3) return g;
        g.bandLo = 2;                                         // bins 0 and 1 are DC and its Hann skirt
        g.bandHi = binCeil (topHz, bin, g.bandLo, last);
        if (g.bandHi <= g.bandLo + 2) return g;
        g.bandBins = g.bandHi - g.bandLo + 1;
        g.spanHarmonics = spanHarmonics;

        // The active-stretch buffer covers every FIXED hypothesis: 50 and 60 Hz at h = 1 and h = 2, each with
        // the search neighbourhood and then the background's reach around a peak displaced to the window edge.
        const double margin = p.searchHz + p.floorSpanHz + bin;
        g.stretchLo = binFloor (kNominal[0] - margin, bin, g.bandLo, g.bandHi);
        g.stretchHi = binCeil (kNominal[1] * 2.0 + margin, bin, g.bandLo, g.bandHi);
        if (g.stretchHi <= g.stretchLo + 2) return g;
        g.stretchBins = g.stretchHi - g.stretchLo + 1;

        g.searchBins       = std::max (1, binCeil (p.searchHz, bin, 1, last));
        g.floorSpanBins    = std::max (2, binCeil (p.floorSpanHz, bin, 2, last));
        g.floorExcludeBins = std::max (1, binCeil (p.floorExcludeHz, bin, 1, last));
        g.floorPairs = g.floorSpanBins - g.floorExcludeBins;
        // FIX: geometry that cannot supply kMinFloorPairs pairs is REFUSED. Accepting it made every
        // prominence test fail for want of a background and the report come back `valid = true,
        // mains = None` on a textbook hum — a false clean produced by a configuration, not by the audio.
        if (g.floorPairs < kMinFloorPairs) return g;
        g.excludeCount = kCandidates * spanHarmonics;         // before merging; the merged count is <= this
        g.ok = true;
        return g;
    }

private:
    struct ChannelState
    {
        std::int64_t frames = 0, finiteFrames = 0, holedFrames = 0, quietFrames = 0, silentFrames = 0;
        std::int64_t stretchCount = 0;                 // absolute: keeps counting past the list's capacity
        std::int64_t eligibleStretches = 0;            // ... of which hold minFramesPerObservation frames
        std::int64_t stretchFrames = 0;
        std::int64_t stretchStart = 0, stretchEnd = 0;
        bool stretchOpen = false;
        HypSummary hyp[kHypotheses] {};
        HumCandidate cand[kCandidates] {};
        HumReason reason = HumReason::NotFinished;
        HumMains mains = HumMains::None;
        int winner = -1;
    };

public:
    // WHAT prepare() ASKS THE HEAP FOR (law 11d), from the function prepare() sizes itself with.
    struct Storage
    {
        bool ok = false;
        SpectrumFrames::Storage frames {};
        std::size_t bandDoubles     = 0;   // double, bandBins * channels
        std::size_t stretchDoubles  = 0;   // double, stretchBins * channels — ONE active buffer per channel
        std::size_t stretchEntries  = 0;   // HumStretch, maxStretches * channels
        std::size_t harmonicEntries = 0;   // HumHarmonic, maxHarmonic * kCandidates * channels
        std::size_t channelEntries  = 0;   // ChannelState, channels
        std::size_t traceEntries    = 0;   // HumFrameTrace, traceCapacity
        std::size_t floorDoubles    = 0;   // double, floorPairs (scratch, reused)
        std::size_t excludeEntries  = 0;   // BinRange, the merged candidate bands
        std::uint64_t bytes() const noexcept
        {
            return frames.bytes()
                 + (std::uint64_t) sizeof (double) * (bandDoubles + stretchDoubles + floorDoubles)
                 + (std::uint64_t) sizeof (HumStretch) * stretchEntries
                 + (std::uint64_t) sizeof (HumHarmonic) * harmonicEntries
                 + (std::uint64_t) sizeof (ChannelState) * channelEntries
                 + (std::uint64_t) sizeof (HumFrameTrace) * traceEntries
                 + (std::uint64_t) sizeof (BinRange) * excludeEntries;
        }
    };

    static Storage storageFor (double sampleRate, int maxChannels, const HumDetectorParams& p) noexcept
    {
        Storage s;
        const Geometry g = geometryFor (sampleRate, p);
        if (! g.ok) return s;
        if (maxChannels < 1 || maxChannels > core::kMaxChannels) return s;
        SpectrumFramesParams fp; fp.fftOrder = g.order; fp.hop = (int) g.hop;
        s.frames = SpectrumFrames::storageFor (sampleRate, maxChannels, fp);
        if (! s.frames.ok) return s;
        const std::size_t ch = (std::size_t) maxChannels;
        s.bandDoubles     = (std::size_t) g.bandBins * ch;
        s.stretchDoubles  = (std::size_t) g.stretchBins * ch;
        s.stretchEntries  = (std::size_t) p.maxStretches * ch;
        s.harmonicEntries = (std::size_t) p.maxHarmonic * (std::size_t) kCandidates * ch;
        s.channelEntries  = ch;
        s.traceEntries    = (std::size_t) p.traceCapacity;
        s.floorDoubles    = (std::size_t) g.floorPairs;
        s.excludeEntries  = (std::size_t) g.excludeCount;
        s.ok = true;
        return s;
    }

    // The PENDING parameters: they take effect at the next prepare() and nowhere else. Every stride, bound
    // and threshold that process(), finish() and the accessors use comes from the set prepare() COMMITTED, so
    // writing this under a live object cannot move an index or a threshold (it used to: maxHarmonic 1 -> 64
    // after prepare() then read entry 127 of a two-entry vector, and maxStretches 1 -> 64 wrote past its end).
    void setParams (const HumDetectorParams& p) noexcept { params_ = p; }
    const HumDetectorParams& params() const noexcept { return params_; }          // pending
    const HumDetectorParams& committedParams() const noexcept { return cfg_; }    // what prepare() accepted

    [[nodiscard]] bool prepare (double sampleRate, int /*maxBlock: nothing is sized by it*/, int maxChannels) noexcept
    {
        prepared_ = false;                                    // law 11b: disarm, validate, write
        const Storage st = storageFor (sampleRate, maxChannels, params_);
        if (! st.ok) return false;
        geom_ = geometryFor (sampleRate, params_);
        if (! geom_.ok) return false;
        cfg_ = params_;                                       // COMMITTED here; nothing else reads params_
        SpectrumFramesParams fp; fp.fftOrder = geom_.order; fp.hop = (int) geom_.hop;
        frames_.setParams (fp);
        if (! frames_.prepare (sampleRate, 0, maxChannels)) return false;

        sampleRate_ = sampleRate;
        channels_   = maxChannels;
        quietLin_   = core::det::pow10 (cfg_.quietThresholdDb / 10.0);
        promLin_    = core::det::pow10 (cfg_.minProminenceDb / 10.0);
        levelLin_   = core::det::pow10 (cfg_.minLevelDbfs / 10.0);
        harmTolHz_  = cfg_.harmonicToleranceBins * geom_.binHz;

        band_.assign (st.bandDoubles, 0.0);
        stretchBuf_.assign (st.stretchDoubles, 0.0);
        stretches_.assign (st.stretchEntries, HumStretch {});
        harmonics_.assign (st.harmonicEntries, HumHarmonic {});
        chans_.assign (st.channelEntries, ChannelState {});
        trace_.assign (st.traceEntries, HumFrameTrace {});
        floorScratch_.assign (st.floorDoubles, 0.0);
        exclude_.assign (st.excludeEntries, BinRange {});
        buildExclusions();
        prepared_ = true;
        reset();
        return true;
    }

    void reset() noexcept
    {
        frames_.reset();
        for (auto& v : band_) v = 0.0;
        for (auto& v : stretchBuf_) v = 0.0;
        for (auto& s : stretches_) s = HumStretch {};
        for (auto& h : harmonics_) h = HumHarmonic {};
        for (auto& c : chans_) c = ChannelState {};
        for (auto& t : trace_) t = HumFrameTrace {};
        traceCount_ = 0;
        finished_ = false;
        if (prepared_)
            for (int cand = 0; cand < kCandidates; ++cand)
                for (int c = 0; c < channels_; ++c)
                    candidateAt (c, cand).nominalHz = kNominal[cand];
    }

    static constexpr int latencySamples() noexcept { return 0; }          // a read-only sink

    // READ-ONLY. Law 11: malformed -> unprepared -> finished -> nch > maxChannels -> n == 0 -> run.
    [[nodiscard]] bool process (const float* const* in, int numChannels, int n) noexcept
    {
        if (numChannels < 0 || n < 0) return false;
        if (! prepared_ || finished_) return false;
        if (numChannels > channels_) return false;
        if (n == 0) return true;                              // an empty call changes NOTHING, stretches included
        for (int i = 0; i < n; ++i)
        {
            for (int c = 0; c < channels_; ++c)
            {
                const bool fed = c < numChannels;
                frames_.push (c, fed ? in[c][i] : 0.0f, fed);
            }
            if (frames_.tick())                               // consume IMMEDIATELY: power() lives until the next
                for (int c = 0; c < channels_; ++c)           // frame closes (SpectrumFrames.h:164)
                    onFrame (c);
        }
        return true;
    }

    // End of stream: close the open stretch at its last selected frame's end, then decide. Idempotent.
    void finish() noexcept
    {
        if (! prepared_ || finished_) return;
        frames_.finish();
        for (int c = 0; c < channels_; ++c)
        {
            if (chans_[(std::size_t) c].stretchOpen) closeStretch (c);
            decide (c);
        }
        finished_ = true;
    }

    // --- the report ---
    bool isPrepared() const noexcept { return prepared_; }
    bool isFinished() const noexcept { return finished_; }
    const Geometry& geometry() const noexcept { return geom_; }
    double binHz() const noexcept { return geom_.binHz; }
    double sampleRate() const noexcept { return sampleRate_; }
    std::int64_t windowSamples() const noexcept { return geom_.n; }
    std::int64_t hopSamples() const noexcept { return geom_.hop; }
    int channels() const noexcept { return channels_; }
    std::int64_t samplesProcessed() const noexcept { return frames_.totalSamples(); }

    HumReport report (int c) const noexcept
    {
        HumReport r;
        if (! prepared_ || ! has (c)) { r.reason = HumReason::NotPrepared; return r; }
        const ChannelState& s = chans_[(std::size_t) c];
        r.binHz = geom_.binHz;
        r.windowSamples = geom_.n;
        r.frames = s.frames;
        r.finiteFrames = s.finiteFrames;
        r.holedFrames = s.holedFrames;
        r.quietFrames = s.quietFrames;
        r.silentFrames = s.silentFrames;
        r.quietStretches = s.stretchCount;
        r.eligibleStretches = s.eligibleStretches;
        r.storedStretches = storedStretchCount (c);
        r.stretchesComplete = s.stretchCount <= (std::int64_t) cfg_.maxStretches;
        r.totalSamples = frames_.totalSamples();
        r.tailUncoveredSamples = frames_.tailUncoveredSamples();
        if (! finished_) { r.reason = HumReason::NotFinished; return r; }
        r.reason = s.reason;
        r.valid = s.reason == HumReason::Ok;
        if (! r.valid) return r;                              // conclusion fields stay at their canonical zero
        r.mains = s.mains;
        if (s.winner >= 0)
        {
            const HumCandidate& w = s.cand[(std::size_t) s.winner];
            r.line = w.base;
            r.fundamentalHz = w.fundamentalHz;
            r.fundamentalObserved = w.fundamentalObserved;
            r.fundamentalDerived = ! w.fundamentalObserved && w.baseFound;
            r.baseHarmonic = w.baseHarmonic;
            r.harmonicsObserved = w.harmonicsObserved;
        }
        return r;
    }

    // The channel the summary speaks for: the greatest aggregate prominence among channels that passed; ties by
    // the lowest index. When none passed, channel 0 — whose report then carries its own reason.
    int strongestChannel() const noexcept
    {
        int best = 0;
        double bestTone = 0.0, bestFloor = 0.0;
        bool any = false;
        for (int c = 0; c < channels_; ++c)
        {
            const ChannelState& s = chans_[(std::size_t) c];
            if (s.reason != HumReason::Ok || s.winner < 0) continue;
            const HumCandidate& w = s.cand[(std::size_t) s.winner];
            if (! (w.base.peakBinPower > 0.0) || ! (w.base.floorPower > 0.0)) continue;
            // compare peak/floor ratios by cross-multiplication, so no division decides the choice
            if (! any || w.base.peakBinPower * bestFloor > bestTone * w.base.floorPower)
            { any = true; best = c; bestTone = w.base.peakBinPower; bestFloor = w.base.floorPower; }
        }
        return best;
    }
    HumReport report() const noexcept { return report (strongestChannel()); }

    HumCandidate candidate (int c, int cand) const noexcept
    {
        if (! has (c) || cand < 0 || cand >= kCandidates) return HumCandidate {};
        return chans_[(std::size_t) c].cand[(std::size_t) cand];
    }
    HumHarmonic harmonic (int c, int cand, int h) const noexcept
    {
        if (! has (c) || cand < 0 || cand >= kCandidates) return HumHarmonic {};
        if (h < 1 || h > cfg_.maxHarmonic) return HumHarmonic {};
        return harmonics_[harmonicIndex (c, cand, h)];
    }
    // The winner's harmonics; an empty entry before finish() or when nothing passed.
    HumHarmonic harmonic (int c, int h) const noexcept
    {
        if (! has (c)) return HumHarmonic {};
        const int w = chans_[(std::size_t) c].winner;
        return w < 0 ? HumHarmonic {} : harmonic (c, w, h);
    }

    std::int64_t stretchCount (int c) const noexcept { return has (c) ? chans_[(std::size_t) c].stretchCount : 0; }
    std::int64_t storedStretchCount (int c) const noexcept
    {
        return has (c) ? std::min<std::int64_t> (chans_[(std::size_t) c].stretchCount, (std::int64_t) cfg_.maxStretches) : 0;
    }
    bool stretchesComplete (int c) const noexcept
    {
        return has (c) && chans_[(std::size_t) c].stretchCount <= (std::int64_t) cfg_.maxStretches;
    }
    HumStretch stretch (int c, std::int64_t i) const noexcept
    {
        if (! has (c) || i < 0 || i >= storedStretchCount (c)) return HumStretch {};
        return stretches_[(std::size_t) c * (std::size_t) cfg_.maxStretches + (std::size_t) i];
    }

    // The pooled power SUM over the accumulated band — every selected frame added once, nothing divided, so a
    // test can compare each bin bit-for-bit and a consumer can plot the evidence. Divide by selectedFrames(c)
    // for the average; a position or a prominence does not need it, an additive dB constant cancelling in both.
    // bandBinLo() is the bin the first element stands for.
    const double* bandPowerSum (int c) const noexcept
    {
        if (! has (c) || band_.empty()) return nullptr;
        return band_.data() + (std::size_t) c * (std::size_t) geom_.bandBins;
    }
    int bandBinLo() const noexcept { return geom_.bandLo; }
    int bandBins()  const noexcept { return geom_.bandBins; }
    std::int64_t selectedFrames (int c) const noexcept { return has (c) ? chans_[(std::size_t) c].quietFrames : 0; }

    // --- the frame trace (a bounded test hook; capacity is a parameter, overflow is data) ---
    std::int64_t traceCount() const noexcept { return traceCount_; }
    std::int64_t storedTraceCount() const noexcept { return std::min<std::int64_t> (traceCount_, (std::int64_t) trace_.size()); }
    bool traceComplete() const noexcept { return traceCount_ <= (std::int64_t) trace_.size(); }
    HumFrameTrace traceAt (std::int64_t i) const noexcept
    {
        return i >= 0 && i < storedTraceCount() ? trace_[(std::size_t) i] : HumFrameTrace {};
    }

private:
    bool has (int c) const noexcept { return prepared_ && c >= 0 && c < channels_; }
    HumCandidate& candidateAt (int c, int cand) noexcept { return chans_[(std::size_t) c].cand[(std::size_t) cand]; }
    std::size_t harmonicIndex (int c, int cand, int h) const noexcept
    {
        return ((std::size_t) c * (std::size_t) kCandidates + (std::size_t) cand) * (std::size_t) cfg_.maxHarmonic
             + (std::size_t) (h - 1);
    }
    static int hypIndex (int cand, int h) noexcept { return cand * kBaseHarmonics + (h - 1); }
    static double hypCentreHz (int cand, int h) noexcept { return kNominal[cand] * (double) h; }

    // The candidate bands the quiet gate ignores: +-searchHz around every h*50 and h*60, MERGED and ascending,
    // so the Parseval fold walks bins in one pass and skips them without ever revisiting a bin.
    void buildExclusions() noexcept
    {
        excludeUsed_ = 0;
        const int cap = (int) exclude_.size();
        for (int h = 1; h <= geom_.spanHarmonics; ++h)     // >= kBaseHarmonics, so 100/120 Hz is always covered
            for (int cand = 0; cand < kCandidates; ++cand)
            {
                const double f = hypCentreHz (cand, h);
                if ((f + cfg_.searchHz) / geom_.binHz < 2.0) continue;
                const int lo = binFloor (f - cfg_.searchHz, geom_.binHz, 2, geom_.bins - 1);
                const int hi = binCeil  (f + cfg_.searchHz, geom_.binHz, 2, geom_.bins - 1);
                if (hi < lo) continue;
                insertRange (lo, hi, cap);
            }
    }
    void insertRange (int lo, int hi, int cap) noexcept
    {
        // insertion into an ascending, disjoint list, merging on touch — O(k^2) over <= 2*maxHarmonic entries
        int at = 0;
        while (at < excludeUsed_ && exclude_[(std::size_t) at].hi + 1 < lo) ++at;
        int end = at;
        while (end < excludeUsed_ && exclude_[(std::size_t) end].lo <= hi + 1)
        {
            lo = std::min (lo, exclude_[(std::size_t) end].lo);
            hi = std::max (hi, exclude_[(std::size_t) end].hi);
            ++end;
        }
        const int removed = end - at;
        if (removed > 0)
        {
            for (int i = at; i + removed < excludeUsed_; ++i) exclude_[(std::size_t) i] = exclude_[(std::size_t) (i + removed)];
            excludeUsed_ -= removed;
        }
        if (excludeUsed_ >= cap) return;                       // cannot happen: cap == kCandidates*maxHarmonic
        for (int i = excludeUsed_; i > at; --i) exclude_[(std::size_t) i] = exclude_[(std::size_t) (i - 1)];
        exclude_[(std::size_t) at] = BinRange { lo, hi };
        ++excludeUsed_;
    }

    //--------------------------------------------------------------------------
    // ONE FRAME, ONE CHANNEL. Everything here reads the producer's power() for the frame that just closed.
    void onFrame (int c) noexcept
    {
        ChannelState& s = chans_[(std::size_t) c];
        ++s.frames;
        const std::int64_t start = frames_.frameStart(), end = frames_.frameEnd();
        const bool finite = frames_.frameFinite (c);
        HumFrameTrace tr;
        tr.frameStart = start; tr.frameEnd = end; tr.channel = c;
        tr.finite = (std::uint8_t) (finite ? 1 : 0);

        if (! finite)
        {
            ++s.holedFrames;
            if (s.stretchOpen) closeStretch (c);
            pushTrace (tr);
            return;
        }
        ++s.finiteFrames;
        const double* p = frames_.power (c);
        const double ms = quietMeasure (p);
        tr.meanSquare = ms;
        // A frame of EXACT digital silence carries no information about hum, so it is not evidence of its
        // absence: it is excluded like a hole rather than selected. Without this, a file whose quiet passages
        // are digital black reported `mains = None` over a -40 dBFS hum in the programme — a clean answer
        // built out of frames that could not have seen anything. Dither and denormals still count: they are
        // a background a line would stand above.
        // This IS an exact-zero predicate, and an exact-zero predicate is normally a cliff — dithered
        // silence has no exact zero, so it slips past. Here the cliff is covered from the other side: a
        // frame whose background is 1e-30 rather than 0 is selected, but no line inside it can be accepted
        // unless the line itself clears minLevelDbfs, so the pair leaves no gap. Widening this into a
        // threshold would instead start discarding the very quiet frames where a faint hum is visible.
        const bool silent = ! (ms > 0.0);
        if (silent) ++s.silentFrames;
        const bool quiet = ms < quietLin_ && ! silent;
        tr.quiet = (std::uint8_t) (quiet ? 1 : 0);
        if (! quiet)
        {
            if (s.stretchOpen) closeStretch (c);
            pushTrace (tr);
            return;
        }

        // A run of consecutive SELECTED frame indices. No gap check is needed and none is written: every
        // closed frame reaches this function, and every path but "quiet and finite" has already closed the
        // stretch, so two selected frames in one run are always consecutive. (A gap check here was dead code
        // — the mutation stand proved it by staying green when it was deleted.)
        if (! s.stretchOpen) openStretch (c, start);
        s.stretchEnd = end;
        ++s.stretchFrames;
        ++s.quietFrames;
        tr.stretchIndex = s.stretchCount;                      // the absolute index, capacity or no capacity

        // the pooled average and the active stretch, bins ASCENDING, one addition per bin per frame
        double* g = band_.data() + (std::size_t) c * (std::size_t) geom_.bandBins;
        double frameBandSum = 0.0;
        for (int k = geom_.bandLo; k <= geom_.bandHi; ++k)
        {
            const double v = p[k];
            g[k - geom_.bandLo] += v;
            frameBandSum += v;
        }
        tr.bandSum = frameBandSum;
        double* sb = stretchBuf_.data() + (std::size_t) c * (std::size_t) geom_.stretchBins;
        for (int k = geom_.stretchLo; k <= geom_.stretchHi; ++k) sb[k - geom_.stretchLo] += p[k];

        // per-FRAME estimates of the four fixed hypotheses: the intra-stretch stillness evidence, and the
        // spread gate a pair of averaged stretches cannot see (two symmetric sweeps average to the same 50.0)
        for (int cand = 0; cand < kCandidates; ++cand)
            for (int h = 1; h <= kBaseHarmonics; ++h)
            {
                const int hi = hypIndex (cand, h);
                const HumPeak pk = locate (p, 0, geom_.bins - 1, hypCentreHz (cand, h),
                                           (double) h * cfg_.toleranceHz, cfg_.searchHz, 1.0);
                tr.hypHz[hi] = pk.hz; tr.hypTone[hi] = pk.tonePower; tr.hypFloor[hi] = pk.floorPower;
                tr.hypBin[hi] = (std::int32_t) pk.bin;
                if (! pk.accepted) continue;
                tr.hypAccepted |= (std::uint32_t) 1u << hi;
                HypSummary& hs = s.hyp[(std::size_t) hi];
                const double f0 = pk.hz / (double) h;
                if (hs.frameObs == 0) { hs.frameMin = f0; hs.frameMax = f0; }
                else { hs.frameMin = std::min (hs.frameMin, f0); hs.frameMax = std::max (hs.frameMax, f0); }
                ++hs.frameObs;
                if (hs.frameObsInStretch == 0) { hs.frameMinInStretch = f0; hs.frameMaxInStretch = f0; }
                else { hs.frameMinInStretch = std::min (hs.frameMinInStretch, f0); hs.frameMaxInStretch = std::max (hs.frameMaxInStretch, f0); }
                ++hs.frameObsInStretch;
            }
        pushTrace (tr);
    }

    // The frame's Hann-weighted mean square, from its own bins, with DC and the candidate bands removed.
    // The fold is P[0] + P[N/2] + 2*sum_{1..N/2-1} P[k] == sum(x^2 w^2)/sum(w^2); bins 0 and 1 and the merged
    // exclusion ranges are the two subtractions the header explains. Ascending, one fixed recurrence.
    double quietMeasure (const double* p) const noexcept
    {
        const int last = geom_.bins - 1;                       // N/2, counted ONCE (it is its own mirror)
        double sum = 0.0;                                      // the interior bins, each worth 2 (its mirror)
        int k = 2;
        for (int e = 0; e < excludeUsed_ && k < last; ++e)
        {
            const BinRange& r = exclude_[(std::size_t) e];
            if (r.hi < k) continue;                            // entirely behind us (clamped to bin 2)
            const int stop = std::min (r.lo, last);
            for (; k < stop; ++k) sum += p[k];
            if (r.hi + 1 > k) k = r.hi + 1;
        }
        for (; k < last; ++k) sum += p[k];
        double total = 2.0 * sum;
        bool lastExcluded = false;
        for (int i = 0; i < excludeUsed_; ++i)
            if (last >= exclude_[(std::size_t) i].lo && last <= exclude_[(std::size_t) i].hi) lastExcluded = true;
        if (! lastExcluded) total += p[last];
        return total;
    }

    //--------------------------------------------------------------------------
    // ONE SPECTRAL LINE. `p` is indexed absolutely from bin `pLo` to `pHi`; the caller passes either a frame's
    // power (pLo = 0) or an accumulated band. Every strict interior local maximum in the search window is
    // enumerated; the strongest one INSIDE the tolerance is the answer, and the strongest one overall is the
    // evidence. Sums, not averages: an additive dB constant leaves the parabola's vertex and every ratio alone.
    // `radiusHz` is the ENUMERATION reach and `tolHz` the ACCEPTANCE distance. They are separate arguments
    // because pass B's window can legitimately be wider than searchHz (harmonicToleranceBins is a parameter);
    // one shared radius silently narrowed it, and a harmonic inside its configured tolerance was never looked at.
    HumPeak locate (const double* p, int pLo, int pHi, double centreHz, double tolHz, double radiusHz,
                    double frames, HumPeak* windowPeakOut = nullptr) const noexcept
    {
        HumPeak best {}, window {};
        const int lo = std::max (pLo + 1, binFloor (centreHz - radiusHz, geom_.binHz, 0, pHi));
        const int hi = std::min (pHi - 1, binCeil  (centreHz + radiusHz, geom_.binHz, 0, pHi));
        for (int k = lo; k <= hi; ++k)
        {
            const double a = p[k - 1 - pLo], b = p[k - pLo], cc = p[k + 1 - pLo];
            if (! (b > a) || ! (b >= cc)) continue;            // strict left, non-strict right: a two-bin
            if (! (b > 0.0)) continue;                         // plateau (a tone exactly between bins) is a peak
            HumPeak pk;
            pk.found = true;
            pk.bin = k;
            const double delta = parabolicDelta (a, b, cc);
            pk.hz = ((double) k + delta) * geom_.binHz;
            pk.tonePower = 2.0 * (a + b + cc);                 // the tone's mean square (A^2/2 for a sine)
            pk.peakBinPower = b;                               // one bin, to be compared with one bin
            // The background is NOT computed here. It is a sort over floorPairs elements, and computing one
            // per local maximum made a wide window quadratic: order 22 at 1 kHz spent 114 s on two frames.
            // Only two maxima can ever be reported — the strongest inside the tolerance and the strongest
            // overall — so the sort happens twice, below, after they are known.
            if (! window.found || pk.peakBinPower > window.peakBinPower) window = pk;
            if (std::fabs (pk.hz - centreHz) <= tolHz
                && (! best.found || pk.peakBinPower > best.peakBinPower)) best = pk;
        }
        for (HumPeak* q : { &best, &window })
        {
            if (! q->found) continue;
            q->floorPower = localFloor (p, pLo, pHi, q->bin);
            // Two tests, and the second is not optional. A RATIO alone certifies numerical residue: a pure
            // 440 Hz tone, or a constant +1.0, leaves a spectrum of ~1e-23 and ~1e-35 elsewhere, and a local
            // maximum of THAT stands 18-37 dB over the residue around it — three fixtures reported confident
            // 50/60 Hz mains from lines at -230 to -348 dBFS. A line must also be loud enough to be a line.
            q->prominent = q->peakBinPower > 0.0 && q->floorPower > 0.0
                        && q->peakBinPower >= q->floorPower * promLin_
                        && q->tonePower >= levelLin_ * (frames > 0.0 ? frames : 1.0);
        }
        best.accepted = best.prominent;                        // it was already inside the tolerance
        if (windowPeakOut != nullptr)
        {
            *windowPeakOut = window;
            if (window.found && window.peakBinPower > 0.0 && window.floorPower > 0.0)
                windowPeakOut->prominenceDb = 10.0 * core::det::log10 (window.peakBinPower / window.floorPower);
        }
        if (best.found && best.peakBinPower > 0.0 && best.floorPower > 0.0)
            best.prominenceDb = 10.0 * core::det::log10 (best.peakBinPower / best.floorPower);
        return best;
    }

    // Three-bin parabola on 10*log10(P). An additive constant cancels, so the accumulated SUM may be used
    // directly. Worst-case bias for a Hann main lobe is ~0.016 bin; the result is clamped to the bin it names.
    static double parabolicDelta (double a, double b, double cc) noexcept
    {
        if (! (a > 0.0) || ! (b > 0.0) || ! (cc > 0.0)) return 0.0;
        const double la = core::det::log10 (a), lb = core::det::log10 (b), lc = core::det::log10 (cc);
        const double den = la - 2.0 * lb + lc;
        if (! (den < 0.0)) return 0.0;                         // not a maximum in the log domain
        const double d = 0.5 * (la - lc) / den;
        if (! std::isfinite (d)) return 0.0;
        return std::clamp (d, -0.5, 0.5);
    }

    // A published linear power is a per-frame AVERAGE; the gates work on the accumulated sums, where the frame
    // count cancels in every ratio, so the single division happens here, once, when the number leaves the object.
    static void normalise (HumPeak& pk, std::int64_t count) noexcept
    {
        if (count <= 0) return;
        const double inv = 1.0 / (double) count;
        pk.tonePower *= inv;
        pk.peakBinPower *= inv;
        pk.floorPower *= inv;
    }

    // The local background at bin kp: the lower median of PAIRED geometric means over symmetric offsets, so a
    // spectral slope cancels instead of biasing the floor onto the quieter side.
    double localFloor (const double* p, int pLo, int pHi, int kp) const noexcept
    {
        int used = 0;
        for (int d = geom_.floorExcludeBins + 1; d <= geom_.floorSpanBins; ++d)
        {
            const int l = kp - d, r = kp + d;
            if (l < pLo || r > pHi) break;                     // shorten BOTH flanks symmetrically
            const double a = p[l - pLo], b = p[r - pLo];
            if (! (a > 0.0) || ! (b > 0.0)) continue;          // a zero or a non-finite bin is no evidence
            if (! std::isfinite (a) || ! std::isfinite (b)) continue;
            if (used >= (int) floorScratch_.size()) break;
            floorScratch_[(std::size_t) used] = std::sqrt (a) * std::sqrt (b);   // not sqrt(a*b): no underflow
            ++used;
        }
        if (used < kMinFloorPairs) return 0.0;                 // unavailable, and said so by being zero
        std::sort (floorScratch_.begin(), floorScratch_.begin() + used);
        return floorScratch_[(std::size_t) ((used - 1) / 2)];  // the LOWER median: no average of two elements
    }

    //--------------------------------------------------------------------------
    void openStretch (int c, std::int64_t start) noexcept
    {
        ChannelState& s = chans_[(std::size_t) c];
        s.stretchOpen = true;
        s.stretchStart = start;
        s.stretchEnd = start;
        s.stretchFrames = 0;
        for (auto& hs : s.hyp) { hs.frameObsInStretch = 0; hs.frameMinInStretch = 0.0; hs.frameMaxInStretch = 0.0; }
        double* sb = stretchBuf_.data() + (std::size_t) c * (std::size_t) geom_.stretchBins;
        for (int i = 0; i < geom_.stretchBins; ++i) sb[i] = 0.0;
    }

    // Fold the closed stretch into the fixed hypotheses' summaries and, if the list has room, record its
    // coordinates. The MEASUREMENT happens whatever the capacity — only the record may be dropped.
    void closeStretch (int c) noexcept
    {
        ChannelState& s = chans_[(std::size_t) c];
        if (! s.stretchOpen) return;
        s.stretchOpen = false;
        const std::int64_t index = s.stretchCount;
        ++s.stretchCount;
        const bool eligible = s.stretchFrames >= (std::int64_t) cfg_.minFramesPerObservation;
        if (eligible) ++s.eligibleStretches;
        const double* sb = stretchBuf_.data() + (std::size_t) c * (std::size_t) geom_.stretchBins;
        for (int cand = 0; cand < kCandidates; ++cand)
            for (int h = 1; h <= kBaseHarmonics; ++h)
            {
                const int hi = hypIndex (cand, h);
                HypSummary& hs = s.hyp[(std::size_t) hi];
                HumPeak win {};
                const HumPeak pk = locate (sb, geom_.stretchLo, geom_.stretchHi, hypCentreHz (cand, h),
                                           (double) h * cfg_.toleranceHz, cfg_.searchHz,
                                           (double) s.stretchFrames, &win);
                if (pk.accepted && eligible)
                {
                    const double f0 = pk.hz / (double) h;
                    if (hs.stretchObs == 0) { hs.stretchMin = f0; hs.stretchMax = f0; }
                    else { hs.stretchMin = std::min (hs.stretchMin, f0); hs.stretchMax = std::max (hs.stretchMax, f0); }
                    ++hs.stretchObs;
                }
                else if (win.found && win.prominent)
                    ++hs.stretchOffTol;                        // contradictory evidence, counted, not censored
                if (hs.frameObsInStretch >= 2)
                    hs.maxIntra = std::max (hs.maxIntra, hs.frameMaxInStretch - hs.frameMinInStretch);
            }
        if (index < (std::int64_t) cfg_.maxStretches)
            stretches_[(std::size_t) c * (std::size_t) cfg_.maxStretches + (std::size_t) index]
                = HumStretch { index, s.stretchStart, s.stretchEnd, s.stretchFrames };
    }

    void pushTrace (const HumFrameTrace& tr) noexcept
    {
        if (traceCount_ >= 0 && traceCount_ < (std::int64_t) trace_.size()) trace_[(std::size_t) traceCount_] = tr;
        ++traceCount_;                                         // keeps counting past the capacity
    }

    //--------------------------------------------------------------------------
    // THE DECISION, once, at finish(). Every candidate is measured whether or not it wins.
    void decide (int c) noexcept
    {
        ChannelState& s = chans_[(std::size_t) c];
        const double* g = band_.data() + (std::size_t) c * (std::size_t) geom_.bandBins;
        const bool haveAverage = s.quietFrames > 0;

        for (int cand = 0; cand < kCandidates; ++cand)
        {
            HumCandidate& k = s.cand[(std::size_t) cand];
            k.nominalHz = kNominal[cand];
            if (! haveAverage) continue;

            // pass A: the comb's base, from h = 1 then h = 2
            for (int h = 1; h <= kBaseHarmonics && ! k.baseFound; ++h)
            {
                HumPeak win {};
                const HumPeak pk = locate (g, geom_.bandLo, geom_.bandHi, hypCentreHz (cand, h),
                                           (double) h * cfg_.toleranceHz, cfg_.searchHz, (double) s.quietFrames, &win);
                if (h == 1) k.windowPeak = win;
                if (! pk.accepted) continue;
                k.baseFound = true;
                k.baseHarmonic = h;
                k.base = pk;
                k.fundamentalHz = pk.hz / (double) h;
                k.fundamentalObserved = h == 1;
            }
            normalise (k.base, s.quietFrames);
            normalise (k.windowPeak, s.quietFrames);
            // pass B: the comb. Around h*f0est when a base was found, else on the NOMINAL grid, so a ripple
            // whose lowest component is out of the base's scope still appears in the evidence.
            const double f0 = k.baseFound ? k.fundamentalHz : kNominal[cand];
            for (int h = 1; h <= cfg_.maxHarmonic; ++h)
            {
                HumHarmonic& hh = harmonics_[harmonicIndex (c, cand, h)];
                hh.index = h;
                const double centre = f0 * (double) h;
                const double tol = k.baseFound ? harmTolHz_ : (double) h * cfg_.toleranceHz;
                hh.inBand = centre > 0.0
                         && binCeil (centre + tol, geom_.binHz, 0, geom_.bandHi)
                              + geom_.floorExcludeBins + kMinFloorPairs <= geom_.bandHi;
                if (! hh.inBand) continue;
                hh.peak = locate (g, geom_.bandLo, geom_.bandHi, centre, tol, std::max (cfg_.searchHz, tol),
                                  (double) s.quietFrames);
                normalise (hh.peak, s.quietFrames);
                if (! hh.peak.accepted) continue;
                ++k.harmonicsObserved;
                if (k.lowestHarmonicObserved == 0) k.lowestHarmonicObserved = h;
            }

            // A comb whose base is out of this scope (h >= 3: a three-pulse rectifier's 150/180 Hz, a
            // six-pulse bridge's 300/360 Hz, a filtered hum) is not an absence either. Two observed
            // harmonics that agree on a fundamental are a comb; one prominent line at a multiple is not.
            if (! k.baseFound)
            {
                int obs = 0;
                double lo = 0.0, hi2 = 0.0;
                for (int h = 1; h <= cfg_.maxHarmonic; ++h)
                {
                    const HumHarmonic& hh = harmonics_[harmonicIndex (c, cand, h)];
                    if (! hh.peak.accepted) continue;
                    const double f = hh.peak.hz / (double) h;
                    if (obs == 0) { lo = f; hi2 = f; k.combFundamentalHz = f; }
                    else { lo = std::min (lo, f); hi2 = std::max (hi2, f); }
                    ++obs;
                }
                k.combWithoutBase = obs >= 2 && hi2 - lo <= cfg_.toleranceHz;
                if (! k.combWithoutBase) k.combFundamentalHz = 0.0;
            }

            // stationarity: >= 2 stretches showed the base, and every estimate — per stretch AND per frame —
            // sits inside one tolerance envelope. The frame spread is what catches two sweeps whose averages
            // agree; the stretch count is what "stands still BETWEEN stretches" means.
            if (! k.baseFound) continue;
            const HypSummary& hs = s.hyp[(std::size_t) hypIndex (cand, k.baseHarmonic)];
            k.stretchObservations = hs.stretchObs;
            k.stretchOffTolerance = hs.stretchOffTol;
            k.frameObservations = hs.frameObs;
            k.stretchSpreadHz = hs.stretchObs > 0 ? hs.stretchMax - hs.stretchMin : 0.0;
            k.frameSpreadHz = hs.frameObs > 0 ? hs.frameMax - hs.frameMin : 0.0;
            k.maxIntraStretchSpreadHz = hs.maxIntra;
            k.stationary = hs.stretchObs >= 2
                        && k.stretchSpreadHz <= cfg_.toleranceHz
                        && k.frameSpreadHz <= cfg_.toleranceHz;
            k.passed = k.stationary;
        }

        // the validity ladder — first match wins
        s.winner = -1;
        s.mains = HumMains::None;
        if (geom_.binHz > kMaxBinHz)            { s.reason = HumReason::InsufficientResolution; return; }
        if (s.frames == 0)                      { s.reason = HumReason::ShorterThanWindow; return; }
        if (s.finiteFrames == 0)                { s.reason = HumReason::AllFramesHoled; return; }
        if (s.stretchCount == 0)                { s.reason = HumReason::NoQuietStretch; return; }
        if (s.stretchCount == 1)                { s.reason = HumReason::SingleQuietStretch; return; }
        if (s.eligibleStretches < 2)            { s.reason = HumReason::StretchesTooShort; return; }
        // A candidate that was FOUND on the pooled average — an accepted, prominent line within the
        // tolerance — and then failed stationarity is not an absence. Three fixtures reached
        // `valid = true, mains = None` that way and all three read as clean: hum drifting 0.56 Hz between
        // stretches, hum present in one of two quiet stretches, and a line swept inside the tolerance band.
        // The evidence stays published in `candidate()`; the verdict does not.
        bool anyFound = false, anyPassed = false;
        for (int cand = 0; cand < kCandidates; ++cand)
        {
            anyFound = anyFound || s.cand[(std::size_t) cand].baseFound;
            anyPassed = anyPassed || s.cand[(std::size_t) cand].passed;
        }
        if (anyFound && ! anyPassed)             { s.reason = HumReason::CandidateNotStationary; return; }
        if (! anyFound
            && (s.cand[0].combWithoutBase || s.cand[1].combWithoutBase))
                                                 { s.reason = HumReason::CombWithoutBase; return; }
        s.reason = HumReason::Ok;
        for (int cand = 0; cand < kCandidates; ++cand)
        {
            const HumCandidate& k = s.cand[(std::size_t) cand];
            if (! k.passed) continue;
            if (s.winner < 0) { s.winner = cand; continue; }
            const HumCandidate& w = s.cand[(std::size_t) s.winner];
            // the greater aggregate prominence wins, by cross-multiplication; a tie keeps the lower nominal
            if (k.base.peakBinPower * w.base.floorPower > w.base.peakBinPower * k.base.floorPower) s.winner = cand;
        }
        if (s.winner >= 0) s.mains = s.winner == 0 ? HumMains::Hz50 : HumMains::Hz60;
    }

    HumDetectorParams params_ {};                          // pending: read only by prepare()/storageFor()
    HumDetectorParams cfg_ {};                             // committed: read by everything at runtime
    Geometry geom_ {};
    SpectrumFrames frames_;
    bool prepared_ = false, finished_ = false;
    double sampleRate_ = 48000.0;
    double quietLin_ = 0.0, promLin_ = 1.0, levelLin_ = 0.0, harmTolHz_ = 0.0;
    int channels_ = 0, excludeUsed_ = 0;
    std::int64_t traceCount_ = 0;
    std::vector<double> band_, stretchBuf_;
    mutable std::vector<double> floorScratch_;
    std::vector<HumStretch> stretches_;
    std::vector<HumHarmonic> harmonics_;
    std::vector<ChannelState> chans_;
    std::vector<HumFrameTrace> trace_;
    std::vector<BinRange> exclude_;
};

} // namespace felitronics::analysis

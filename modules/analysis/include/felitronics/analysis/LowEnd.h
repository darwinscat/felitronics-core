// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/analysis/SpectrumFrames.h>
#include <felitronics/core/DetMath.h>
#include <felitronics/core/Config.h>
#include <felitronics/core/StateGrid.h>
#include <felitronics/eq/Crossover2.h>
#include <felitronics/stereo/MidSide.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace felitronics::analysis
{

//==============================================================================
// felitronics::analysis::LowEnd — the two questions a LACQUER asks of a master's bottom end:
//
//   1. HOW WIDE IS THE BASS. A cutter head writes the mono sum LATERALLY and the difference VERTICALLY.
//      Vertical modulation is what lifts the stylus out of the groove, so out-of-phase low end is not a
//      matter of taste — the cut either survives it or it does not. This half splits the programme at
//      `crossoverHz` with an LR4 (eq::Crossover2), forms Mid/Side, and reports the SIDE share of the low
//      band's energy: over absolute 10 ms blocks, as a duration-weighted distribution, and integrated.
//   2. WHICH NOTE OWNS THE BOTTOM. 30..300 Hz folded to semitone bands, each band's energy in the
//      cuttable (Mid) and the vertical (Side) axis, the loudest band, and how far it stands above the
//      rest. A cutting engineer's first question about a dominant bass note is whether it is LATERAL,
//      and that falls out of the same table.
//
// WHAT IT REFUSES TO SAY. Nothing here is a verdict. The tonic of the key IS normally the loudest note
// in the bass — that is music, not a defect. A kick's fundamental is not a note at all and sits in the
// same bins; it is reported as the band it lands in, and nothing calls it wrong. Every threshold in
// this file DEFINES the instrument (a crossover frequency, a note range, a reference pitch) and is a
// named parameter with a documented default; no number in the report reads as good or bad.
//
// THE PUBLISHED WIDTH NUMBER IS A FRACTION, NOT A RATIO. The request was "Side-to-Mid energy ratio",
// and that ratio cannot be a number: L = -R is a real master (and the worst one for a lacquer), it
// makes Mid EXACTLY zero, and S/M is then +inf. So the report publishes both raw energies and the
// bounded SIDE ENERGY FRACTION f = S/(M+S) in [0, 1], which is monotone in r = S/M (f = r/(1+r)) and
// therefore orders programmes identically with no pole. r = f/(1-f) is one line away for a consumer
// that wants it. The scale: f = 0 is a mono bottom (perfectly lateral), f = 0.5 is a hard-panned OR an
// uncorrelated bottom (equal lateral and vertical energy — the fraction cannot tell those two apart;
// analysis::CorrelationMeter can), f = 1 is pure anti-phase, i.e. pure vertical.
//   * a MONO programme is not 0/0. Its Side is exactly 0 and its Mid is positive, so f is exactly 0 and
//     `widthValid()` is true — zero width is the RIGHT answer, not "undefined".
//   * A ZERO-ENERGY LOW BAND is the 0/0, and there a zero would lie, so it does not get one:
//     `widthReason() == NoEnergy`. Digital silence is the obvious case but NOT the only one: eq::Svf
//     keeps its state in FLOAT, so a programme below that state's own floor filters to exact zero even
//     though its raw energy is positive — measured, every sample at 2^-149 gives rawMidEnergy 2^-281
//     and lowMidEnergy exactly 0. NoEnergy therefore means "no energy reached the low band", which is
//     what the field can honestly claim; rawMidEnergy() is published beside it so the two are
//     distinguishable.
// NOTE this is an ENERGY fraction. `analysis::StereoSums::width` (StereoColumns.h:59) is an AMPLITUDE
// fraction, sqrt(S)/(sqrt(M)+sqrt(S)); it is a bit-exact port of a JavaScript spec and stays as it is.
// The amplitude form is the more sensitive of the two at small side levels (at S/M = -20 dB it reads
// 0.091 where the energy fraction reads 0.0099); both raw energies are published, so either is derivable.
//
// THE LOW BAND IS LR4-SHAPED, NOT A BRICK WALL, and the numbers say how much that matters. The LR4
// low-pass power response is |H|^2 = 1/(1+r^4)^2 with r = tan(pi*f/fs)/tan(pi*fc/fs) (the SVF's
// prewarped ratio, Svf.h:61 — NOT f/fc, which is the analogue prototype and is wrong near Nyquist).
// At fc = 120 Hz / 48 kHz that is -6.0 dB at 120 Hz, -24.6 dB at 240 Hz, -49.6 dB at 480 Hz and
// -73.7 dB at 1 kHz. So a programme whose wide content reaches down to 240 Hz WILL show side energy in
// this "low band" — against an equal-energy 60 Hz mono bass that reads f ~ 0.0039, and that is the
// measurement being correct, not leaking. The high band's own Mid/Side pair is published beside the low
// one precisely so a consumer can see where the width lives.
//   * DO NOT ADD THE TWO BANDS AND EXPECT THE INPUT. LR4 sums to an ALLPASS in amplitude, so
//     |H_lp|^2 + |H_hp|^2 = (1+r^8)/(1+r^4)^2, which is 1/2 at the crossover, not 1 (Crossover2.h:18).
//     The raw (unfiltered) Mid/Side energies are published as the full-band reference instead.
//   * AND THE INTEGRAL INCLUDES THE FILTER'S OWN STARTUP. The crossover begins at sample 0 with no
//     history, and while it charges it passes content the settled filter rejects. That is honest — the
//     programme really does start there — but it DOMINATES a near-zero side fraction: measured on a 4 s
//     file of mono 82 Hz bass plus anti-phase 900 Hz, the whole-file lowSideFraction() is 1.75e-6 while
//     the SETTLED value is 5.32e-8 (the analytic prediction is 5.3212e-8, matched to four figures from
//     20 ms on) — the first 10 ms block alone holds 96.65 % of the file's entire low Side energy, and
//     the first 100 ms 97.04 %. So the
//     integral is the answer for "what is on this record", and the 10 ms SERIES is the answer for "how
//     wide is the bass where it is playing": skip the first few blocks and the two agree. Nothing is
//     dropped here on the instrument's own initiative — the coordinates are published and the consumer
//     decides.
//
// THE SEMITONE FOLD. Bands are the MIDI notes whose centre f(n) = tuningHz * 2^((n-69)/12) lies inside
// [lowNoteHz, highNoteHz]; band n spans f(n) * 2^(-1/24) .. f(n) * 2^(+1/24). A frame's bins are
// integrated into a band with FRACTIONAL edge overlap in the POWER domain — the repository's convention
// (MultiResSpectrumPane.h:405) — under that convention's own assumption, that a bin's power is uniform
// across its cell [(k-1/2)*binHz, (k+1/2)*binHz].
//   * NOT by prefix sums. The pane needs them because a display queries arbitrary overlapping bands at
//     arbitrary frequencies; here the bands are FIXED and DISJOINT, so prepare() precomputes one flat
//     (bin, weight, moment-frequency) table and a frame is a single pass over it. And the prefix form
//     cannot carry a correct FIRST MOMENT: a fractional edge cell's power sits at the midpoint of the
//     OVERLAP, not at the bin's centre, and weighting P_k by the bin centre can place a centroid
//     outside its own band (half a bin is 0.183 Hz at 2^17/48 kHz, ~10 cents at 30 Hz).
//   * ONE-SIDED, FOLDED. SpectrumFrames does not fold (bin k is bin k), so a real tone's power splits
//     between k and N-k. The weights carry a factor 2 for 0 < k < N/2 and 1 for DC and Nyquist, which
//     makes a band's energy a genuine mean-square contribution: a full-scale sine of amplitude A inside
//     one band reads A^2/2.
//     THE FULL SUM IS THE WINDOW-WEIGHTED MEAN SQUARE, AND THE TABLE IS NOT THE FULL SUM. By Parseval
//     the folded sum over EVERY bin is sum(w*x)^2 / sum(w^2) — the mean square of x weighted by w^2,
//     published as frameEnergy(). Disjoint bands add up to that only if they COVER the spectrum, and the
//     published table covers [lowNoteHz, highNoteHz] alone: at the default range and 48 kHz that is
//     1.1345 % of Nyquist, so totalBandEnergy() is the energy of the NOTE RANGE and not of the frame.
//     bandRangeShare() is the ratio, and on a flat spectrum it equals exactly that covered fraction.
//   * RESOLUTION IS THE OBSERVATION LENGTH, and the default is the smallest order that has any. A
//     semitone at 30 Hz is 1.73 Hz wide; a Hann main lobe is 4 bins. At 48 kHz only fftOrder >= 17
//     (0.366 Hz bins, 4.87 bins per band) fits the lobe inside the band, and a tone at a band centre
//     then keeps 99.96 % of its power in its own band (measured against a direct DFT). Orders below
//     that do not resolve the bottom of the range and `underResolvedBands()` counts them rather than
//     hiding it. Zero-padding would not help: it interpolates a peak, it does not separate two tones.
//   * A TRANSIENT'S BAND ENERGY IS MODULATED UP TO 3.01 dB BY WHERE IT FALLS ON THE FRAME GRID. At the
//     default 50 % hop a click at a frame boundary is weighted w^2 = 1 by one frame and 0 by the next,
//     while a click a quarter-window later is weighted 0.5^2 by each of two frames: 1 against 0.5, or
//     exactly 3.0103 dB, over the same frame count. Hann is COLA at 50 % overlap for the WINDOW, not for
//     its square, so no hop choice removes this. It is deterministic and absolute — the grid does not
//     move with the caller's slicing, so law 8a is untouched — but it means the band energy of the click
//     half of a kick carries up to 3 dB of grid phase, and a consumer comparing two transient-heavy
//     programmes should know it. Steady tones are unaffected (they are present in every frame).
//   * WHAT STILL SPREADS IS THE BAND EDGE. A tone near a semitone boundary splits roughly 50/50 between
//     two bands, depressing the peak and raising its neighbour. The whole band table is published, so a
//     split is visible; the peak's centroid recovers a tuning offset well inside the band (true 10/20/30
//     cents read 10.0/19.9/28.6 at B0) and compresses toward the centre as it approaches +-50 cents.
//     It is a spectral centroid, not a tuner.
//   * THE SEMITONE GRID HAS A TILT, AND IT IS NAMED. Band width grows with frequency, so under a FLAT
//     spectrum the top band of the default range holds 3.17x (5.02 dB) the energy of the median band
//     with no note present at all. That is why the background is a median of DENSITIES (energy per Hz)
//     scaled to the peak band's width — under white noise that reads 1.0, as it should — and why the
//     argmax is published twice: `peakBand()` maximises ENERGY (what a tone does, since a tone's energy
//     is independent of its band's width) and `peakDensityBand()` maximises DENSITY (what noise does).
//     When the two disagree, the bottom end is noise-like rather than tonal, and that is evidence.
//   * AND A RANGE THAT HOLDS ONLY ROUND-OFF HOLDS NOTHING. Pure DC, pure Nyquist and a signal at the
//     float denormal floor all leave the 30..300 Hz table with a POSITIVE total that is pure transform
//     round-off: measured, DC gives a bandRangeShare of 1.3e-34 and Nyquist 1.9e-38, and `> 0` alone
//     then named F#3 as the dominant note of a signal that has no note at all. A double FFT resolves a
//     bin to about 1e-15 of the largest bin in amplitude, i.e. 1e-30 in power, so a share at or below
//     kNoteFloorShare = 1e-24 is the transform's own floor — six decades above it, and twenty-seven
//     below the 1.1e-7 a single impulse produces. Below it the answer is NoEnergy and a reason, not a
//     note. This is a threshold that DEFINES the instrument's resolution, published as a constant with
//     its derivation, not a verdict on the programme.
//   * THE DOMINANCE RATIO IS NOT A STORED FIELD, for the same reason S/M is not: the median energy of a
//     synthetic tone in digital silence is exactly 0 and the ratio is then +inf. The report publishes
//     `peakBandEnergy()`, `backgroundDensity()` and `peakBandWidthHz()` — divide them if you want it —
//     plus the bounded `peakShare()` = peak / total, and the runner-up band, which is what says whether
//     the "dominant" note has a rival (a split tone, or a two-note bass).
//
// THE PAIR IS CHANNEL 0 AND CHANNEL 1, AND NOTHING ELSE. L = channel 0, R = channel 1, or R = L for a
// mono object. Channels 2.. take no part in either half: this is an L/R lacquer instrument, and letting
// a surround or LFE channel decide the dominant note (or poison a frame) would measure a record nobody
// is cutting. `analysedChannels()` says so out loud.
//
// FORM (law 11d): setParams / prepare / process / finish / reset, every entry point [[nodiscard]] bool,
// and a public `Storage` + `storageFor()` published from the function prepare() sizes itself with.
// prepare() allocates; process() and finish() do not. OFFLINE (message thread): the default parameters
// ask the heap for ~7 MB.
//
// LAW 8a — bit-identical under ARBITRARY re-slicing, which is stronger than law 11(a)'s same-boundaries
// promise (DSP-ARCHITECTURE.md:273). One integer clock; the outer loop is over SAMPLES and the inner over
// channels; `maxBlock` sizes nothing at all; the 10 ms grid and the frame schedule are absolute; the
// denormal flush rides core::StateGrid rather than the end of process() (eq::Crossover2 does NOT run
// that cadence itself — Crossover2.h:58 hands `flushDenormals()` to its owner); and THE TAIL IS NOT
// TRANSFORMED — finish() invents no spectral frame and publishes `tailUncoveredSamples()` instead, while
// the time-domain half covers every sample up to T and closes its one partial block at [floor(T/B)*B, T).
struct LowEndParams
{
    double crossoverHz = 120.0;     // the LR4 split. Below it is the lacquer's phase question.
    // TWENTY SINCE K9, WHERE IT WAS THIRTY. E0 is 20.60 Hz and a lacquer's vertical hazard lives under B0,
    // so a table that starts at 30 could not see it. The price is four bands at the bottom narrower than a
    // Hann main lobe at 48 kHz — counted by underResolvedBands() and bounded by resolvedAboveHz(), not
    // hidden. Raising fftOrder to 18 would resolve them and was REFUSED: it doubles the window to 5.46 s,
    // which halves the frame count a duty is measured over and smears the intermittency duty exists to find.
    //
    // NOTE FOR A CALLER THAT ADDRESSES BANDS BY POSITION: it must not. Band 0 was MIDI 23 and is now MIDI
    // 16 — every index moved by seven semitones. `LowEndBand::midi` and `centreHz` are the addresses.
    double lowNoteHz   = 20.0;      // the semitone range: every note whose CENTRE lies in [low, high]
    double highNoteHz  = 300.0;     // defaults give MIDI 16..62 (E0..D4), 47 bands
    double tuningHz    = 440.0;     // A4. f(n) = tuningHz * 2^((n-69)/12), MIDI numbering (60 = C4)
    int    fftOrder    = 17;        // N = 1 << fftOrder. 17 is the smallest that resolves a semitone at
                                    // 30 Hz at 48 kHz — see "RESOLUTION" above.
    int    hop         = 0;         // 0 means N/2
    // K9. How far under a FRAME'S loudest band a band may sit and still count as present in that frame.
    // A parameter and not a constant because it is the consumer's hypothesis: 20 dB is its starting point,
    // and it will calibrate. Turned into a linear ratio ONCE, at prepare(), so no logarithm stands in a
    // per-frame decision.
    double dutyThresholdDb = 20.0;
    // K9. How many leading 10 ms blocks the HISTOGRAM ignores. Default 0, so nothing already shipped moves.
    //
    // ONLY THE HISTOGRAM, deliberately. The series, the integrals and the extrema stay complete: the
    // instrument's standing rule is that nothing is dropped on its own initiative — `worstFractionBlock()`
    // is a COORDINATE into a programme, and a coordinate that silently skipped a prefix would be a lie.
    // What the histogram is, and only it, is a population a caller reads a time fraction from, and a
    // caller that wants that fraction free of the crossover's own charge-up says so HERE, in a number it
    // chose. `settlingBlocks()` derives that number from the filter instead of leaving it to be guessed:
    // "the first 10 blocks" was a guess, and it is three to five times too many at 100–150 Hz and two
    // times too few at 20 Hz.
    int    skipBlocks = 0;
    int    maxBlocks   = 1 << 16;   // capacity of the stored 10 ms series: 10.9 min. Everything that is
                                    // not the series itself (integrals, histogram, extrema, counters)
                                    // keeps going past it — law 11, exhaustion is data.
};

// Why a field is not a number. Never NaN: an invalid field has a documented canonical value, and the
// reason is what carries the meaning.
enum class LowEndReason : std::uint8_t
{
    Ok                = 0,
    NotFinished       = 1,  // finish() has not been called yet
    NoFiniteSamples   = 2,  // not one sample of the L/R pair was usable
    NoEnergy          = 3,  // every relevant energy is exactly zero — digital silence, the 0/0
    ShorterThanWindow = 4,  // no spectral frame ever closed: the programme is shorter than the window
    NoUsableFrames    = 5,  // frames closed, but every one of them held a hole
    Overflowed        = 6,  // an accumulated band quantity left the finite range (see reduceBands)
};

// One 10 ms block of the LOW band. Raw energies; the fraction is derived so nothing is lost to it.
struct LowEndBlock
{
    std::int64_t index         = 0;
    std::int64_t start         = 0;     // first sample, in samples since reset()
    std::int64_t samples       = 0;     // the block's own length — the LAST block is short
    std::int64_t finiteSamples = 0;     // samples that reached the accumulators
    std::int64_t holes         = 0;     // samples that did not (non-finite, or a channel was absent)
    double       midEnergy     = 0.0;   // sum of m*m over the finite samples of the LOW band
    double       sideEnergy    = 0.0;   // sum of s*s, ditto
    bool         valid         = false; // no holes at all, and at least one sample

    // S/(M+S) in [0, 1]. Exactly 0 for a mono or a silent block; see the header note on the 0/0.
    double sideFraction() const noexcept
    {
        const double t = midEnergy + sideEnergy;
        return t > 0.0 ? sideEnergy / t : 0.0;
    }
    double energy() const noexcept { return midEnergy + sideEnergy; }
};

// One semitone band's accumulated spectrum, in mean-square units (see the FOLDED note above).
struct LowEndBand
{
    int          midi       = 0;
    double       centreHz   = 0.0;      // f(n), the nominal note
    double       widthHz    = 0.0;      // f(n)*2^(1/24) - f(n)*2^(-1/24)
    double       binsPerBand = 0.0;     // widthHz / binHz — under 4 and the Hann lobe does not fit
    double       midEnergy  = 0.0;      // the CUTTABLE axis, mean over the used frames
    double       sideEnergy = 0.0;      // the VERTICAL axis
    double       energy     = 0.0;      // midEnergy + sideEnergy
    double       density    = 0.0;      // energy / widthHz — the tilt-free quantity
    double       centroidHz = 0.0;      // energy-weighted mean frequency INSIDE the band
    double       centsOffset = 0.0;     // 1200*log2(centroidHz/centreHz), 0 when the band is empty
    double sideFraction() const noexcept
    {
        const double t = midEnergy + sideEnergy;
        return t > 0.0 ? sideEnergy / t : 0.0;
    }
};

// The law-8a observer. Fired at the moment a block CLOSES and at the moment a frame is CONSUMED — never
// at the exit of process(), which is the boundary the whole contract exists to be independent of.
// DIAGNOSTIC: the report is the product; this exists so a test can compare every intermediate, because
// comparing final reports is not enough (docs/LAW8-KWEIGHTING.md:78 — identical LUFS and dBTP over 5
// changed block energies out of 97). The pointers are scratch owned by the object: COPY what you need
// during the call, and never compare the addresses.
struct LowEndTrace
{
    enum class Kind : std::uint8_t { Block = 0, Frame = 1 };
    Kind         kind          = Kind::Block;
    std::int64_t index         = 0;         // block index, or frame index
    std::int64_t start         = 0;         // first sample covered
    std::int64_t end           = 0;         // one past the last
    bool         valid         = false;     // block: no holes. frame: both axes finite, so it was used
    std::int64_t finiteSamples = 0;         // Block only
    std::int64_t holes         = 0;         // Block only
    double       midEnergy     = 0.0;       // Block: the LOW band's raw mid energy over the block
    double       sideEnergy    = 0.0;       // Block: the LOW band's raw side energy
    double       frameEnergy   = 0.0;       // Frame: this frame's own window-weighted mean square
    const double* bandMid      = nullptr;   // Frame: bandCount raw per-frame band powers, mid axis
    const double* bandSide     = nullptr;   // Frame: ditto, side axis
    int          bandCount     = 0;         // Frame only
};

class LowEnd
{
private:
    // Declared first because the public Storage below sizes itself in terms of them.
    struct BlockRecord { double mid = 0.0, side = 0.0; std::int32_t samples = 0, finite = 0; };
    struct BinWeight   { double weight = 0.0, momentHz = 0.0; std::int32_t bin = 0; };

    static constexpr int    kAxes     = 2;                    // Mid and Side — the two things measured
    static constexpr int    kSums     = 6;                    // the six integral energies
    static constexpr double kSemiUp   = 1.0293022366434921;   // 2^( 1/24)
    static constexpr double kSemiDown = 0.9715319411536058;   // 2^(-1/24)

public:
    using TraceFn = void (*) (void* user, const LowEndTrace& t);

    static constexpr double kBlockMs      = 10.0;    // NOT a parameter. The grid is shared: it is
                                                     // lround(0.01*fs) samples, exactly the sub-hop
                                                     // analysis::LoudnessMeter builds (LoudnessMeter::storageFor),
                                                     // so the two instruments name the same intervals.
    static constexpr int    kHistogramBins = 100;    // the side fraction over [0, 1], fixed edges
    static constexpr double kMinSampleRate = core::kMinSampleRate;   // P51: the core's floor, one number
    static constexpr double kMaxSampleRate = 768000.0;
    static constexpr int    kMaxBlocksLimit = 1 << 24;
    static constexpr int    kLobeBins      = 4;      // a Hann main lobe, in bins
    static constexpr double kMinCrossoverHz = 1.0;   // eq::Svf clamps a cutoff below 1 Hz (Svf.h:61), so
                                                     // accepting less would make crossoverHz() report a
                                                     // filter that is not the one running
    static constexpr double kNoteFloorShare = 1.0e-24;   // below this share of frameEnergy() the note
                                                         // range holds only the double transform's own
                                                         // round-off — see the header note
    static constexpr double kMinNoteHz      = 1.0;   // below this the band-edge FREQUENCIES underflow the
                                                     // first moment to zero while the energies stay
                                                     // positive, and a centroid then leaves its own band

    //==============================================================================
    // WHAT prepare() ASKS THE HEAP FOR (law 11d), from the function prepare() sizes AND validates itself
    // with, so a caller budgeting memory reads the numbers the object is actually built from. The nested
    // SpectrumFrames ask is part of it — at the default order it is 5 505 040 bytes on its own and
    // leaving it out would make this a promise with no budget.
    struct Storage
    {
        bool ok = false;
        SpectrumFrames::Storage frames {};
        std::size_t blockRecords  = 0;      // BlockRecord, maxBlocks of them
        std::size_t bands         = 0;      // LowEndBand, the published table
        std::size_t binWeights    = 0;      // BinWeight, the precomputed fold table
        std::size_t bandBinCounts = 0;      // int, how many weights each band owns
        std::size_t accDoubles    = 0;      // 3 * bands: the mid / side / moment accumulators
        std::size_t traceDoubles  = 0;      // 2 * bands: one frame's raw band powers, for the observer
        std::size_t sortDoubles   = 0;      // bands: the median's scratch
        std::size_t dutyCounts    = 0;      // bands: how many frames each band was present in (int64)
        std::size_t dutyRatios    = 0;      // bands: the sum of its share of the frame max, over those frames
        std::int64_t blockSamples = 0;      // lround(0.01*fs), published because every coordinate uses it
        int bandCount             = 0;
        std::uint64_t bytes() const noexcept
        {
            return frames.bytes()
                 + (std::uint64_t) blockRecords * sizeof (BlockRecord)
                 + (std::uint64_t) bands * sizeof (LowEndBand)
                 + (std::uint64_t) binWeights * sizeof (BinWeight)
                 + (std::uint64_t) bandBinCounts * sizeof (int)
                 + (std::uint64_t) sizeof (double) * ((std::uint64_t) accDoubles + (std::uint64_t) traceDoubles
                                                    + (std::uint64_t) sortDoubles + (std::uint64_t) dutyRatios)
                 + (std::uint64_t) dutyCounts * sizeof (std::int64_t);
        }
    };

    // ok == false on exactly the arguments prepare() refuses, and then every count is zero.
    // A NaN parameter fails every one of these tests, which is the point: eq::Svf's setParams clamps a
    // frequency but walks a NaN straight into tan() and poisons its coefficients PERMANENTLY — no flush
    // and no reset() reaches that (Svf.h:44). It is refused here instead.
    static Storage storageFor (double sampleRate, int maxChannels, const LowEndParams& p) noexcept
    {
        Storage s;
        if (! (sampleRate >= kMinSampleRate && sampleRate <= kMaxSampleRate)) return s;
        if (maxChannels < 1 || maxChannels > core::kMaxChannels) return s;
        if (p.maxBlocks < 0 || p.maxBlocks > kMaxBlocksLimit) return s;
        // K9's two parameters, REFUSED HERE RATHER THAN IGNORED LATER — and they were neither, which is
        // the defect a consumer found by passing them. `skipBlocks = -1` compared as `blockIndex >= -1`,
        // true of every block, so the parameter was accepted and then did nothing; a caller reading its own
        // -1 back out of params() would have believed the histogram was cut. `dutyThresholdDb = -5` makes
        // the linear gate 10^(0.5) = 3.16, which no band can reach against its own frame's maximum, so
        // every duty read 0 and the measurement was silently empty. A parameter that is taken and then has
        // no effect is worse than one that is refused: the refusal is visible.
        if (p.skipBlocks < 0 || p.skipBlocks > kMaxBlocksLimit) return s;
        if (! (p.dutyThresholdDb >= 0.0) || ! (p.dutyThresholdDb <= 200.0)) return s;
        if (! (p.crossoverHz >= kMinCrossoverHz) || ! (p.crossoverHz <= 0.49 * sampleRate)) return s;   // Svf would clamp either end silently
        if (! (p.tuningHz >= kMinNoteHz) || ! (p.tuningHz < sampleRate)) return s;
        if (! (p.lowNoteHz >= kMinNoteHz) || ! (p.highNoteHz > p.lowNoteHz)) return s;
        // the top band must fit under Nyquist, or its energy would be a clipped band pretending to be whole
        if (! (p.highNoteHz * kSemiUp < 0.5 * sampleRate)) return s;

        const SpectrumFrames::Storage fs = SpectrumFrames::storageFor (sampleRate, kAxes, framesParams (p));
        if (! fs.ok) return s;

        int loMidi = 0, hiMidi = 0;
        if (! midiRange (p, loMidi, hiMidi)) return s;
        const int bands = hiMidi - loMidi + 1;
        if (bands < 2) return s;                                   // a "dominant note among one" is not a measurement

        const std::int64_t blockSamples = blockSamplesFor (sampleRate);
        if (blockSamples < 1) return s;

        const std::int64_t n = (std::int64_t) 1 << framesParams (p).fftOrder;
        const double binHz = sampleRate / (double) n;
        std::size_t weights = 0;
        for (int b = 0; b < bands; ++b)
        {
            int first = 0, count = 0;
            bandBins (noteHz (p, loMidi + b), binHz, (int) (n / 2 + 1), first, count);
            weights += (std::size_t) count;
        }

        s.ok = true;
        s.frames = fs;
        s.blockRecords  = (std::size_t) p.maxBlocks;
        s.bands         = (std::size_t) bands;
        s.binWeights    = weights;
        s.bandBinCounts = (std::size_t) bands;
        s.accDoubles    = 3u * (std::size_t) bands;
        s.traceDoubles  = 2u * (std::size_t) bands;
        s.sortDoubles   = (std::size_t) bands;
        s.dutyCounts    = (std::size_t) bands;
        s.dutyRatios    = (std::size_t) bands;
        s.blockSamples = blockSamples;
        s.bandCount    = bands;
        return s;
    }

    void setParams (const LowEndParams& p) noexcept { params_ = p; }   // structural: at the next prepare()

    // The observer is not a measurement parameter and may be set at any time; it is read, never stored
    // into the report. nullptr (the default) is no observer.
    void setTrace (TraceFn fn, void* user) noexcept { trace_ = fn; traceUser_ = user; }

    [[nodiscard]] bool prepare (double sampleRate, int /*maxBlock: nothing is sized by it*/, int maxChannels) noexcept
    {
        // Law 11b: disarm, validate, write — AND DISARM MEANS THE REPORT TOO, the way SourceForensics
        // already spells it. A refused prepare() after a finished measurement used to leave isFinished()
        // answering true and the previous report still readable, so an instance reconfigured with bad
        // arguments kept certifying the programme before it. Found by the release round; P71 had already
        // named the same shape a defect.
        prepared_ = false;
        finished_ = false;
        channels_ = 0;
        const Storage st = storageFor (sampleRate, maxChannels, params_);
        if (! st.ok) return false;

        frames_.setParams (framesParams (params_));
        if (! frames_.prepare (sampleRate, 0, kAxes)) return false;

        sampleRate_   = sampleRate;
        channels_     = maxChannels;
        blockSamples_ = st.blockSamples;
        bandCount_    = st.bandCount;
        binHz_        = frames_.binHz();
        midiLo_       = 0;
        {
            int lo = 0, hi = 0;
            if (! midiRange (params_, lo, hi)) return false;
            midiLo_ = lo;
        }

        blocks_.assign (st.blockRecords, BlockRecord {});
        bands_.assign (st.bands, LowEndBand {});
        weights_.assign (st.binWeights, BinWeight {});
        traceMid_.assign (st.bands, 0.0);
        traceSide_.assign (st.bands, 0.0);
        sort_.assign (st.sortDoubles, 0.0);
        bandCountBins_.assign (st.bandBinCounts, 0);
        accMid_.assign (st.bands, 0.0);
        accSide_.assign (st.bands, 0.0);
        accMoment_.assign (st.bands, 0.0);
        dutyCount_.assign (st.dutyCounts, 0);
        dutyRatio_.assign (st.dutyRatios, 0.0);
        // The threshold becomes a linear ratio HERE and never again: a decision taken per frame per band
        // must not carry a logarithm, and det::pow10 is the one this repo agrees on across libms.
        dutyShare_  = core::det::pow10 (-params_.dutyThresholdDb / 10.0);
        dutyFrames_ = 0;

        if (! buildBands()) return false;      // law 11b: prepared_ is still false here

        xover_.prepare (sampleRate, kAxes);                           // two axes: Mid and Side
        xover_.setFrequency ((float) params_.crossoverHz);

        prepared_ = true;
        reset();
        return true;
    }

    void reset() noexcept
    {
        totalSamples_ = 0;
        nextBlockEnd_ = blockSamples_;
        blockIndex_   = 0;
        blockCount_   = 0;
        finished_     = false;
        openMid_ = openSide_ = 0.0;
        openFinite_ = openHoles_ = 0;
        for (int i = 0; i < kSums; ++i) { sum_[i] = 0.0; comp_[i] = 0.0; }
        finiteSamples_ = holeSamples_ = nonFiniteSamples_ = filterNonFinite_ = absentSamples_ = 0;
        firstHole_ = lastHole_ = -1;
        for (int i = 0; i < kHistogramBins; ++i) hist_[i] = 0;
        skippedBlocks_ = 0;
        std::fill (dutyCount_.begin(), dutyCount_.end(), (std::int64_t) 0);
        std::fill (dutyRatio_.begin(), dutyRatio_.end(), 0.0);
        dutyFrames_ = 0;
        histSamples_ = 0;
        worstFrac_ = -1.0; worstFracBlock_ = -1; worstFracEnergy_ = 0.0;
        peakEnergy_ = -1.0; peakEnergyBlock_ = -1; peakEnergyFrac_ = 0.0;
        peakSide_ = -1.0; peakSideBlock_ = -1;
        peakSideAmp_ = 0.0; peakSideAmpAt_ = -1;
        usedFrames_ = 0; holedFrames_ = 0;
        accFrameEnergy_ = 0.0; frameEnergy_ = 0.0; bandRangeShare_ = 0.0; frameTotal_ = 0.0;
        for (std::size_t i = 0; i < accMid_.size(); ++i) { accMid_[i] = 0.0; accSide_[i] = 0.0; accMoment_[i] = 0.0; }
        for (auto& b : bands_) { b.midEnergy = 0.0; b.sideEnergy = 0.0; b.energy = 0.0; b.density = 0.0; b.centroidHz = 0.0; b.centsOffset = 0.0; }
        peakBand_ = -1; peakDensityBand_ = -1; secondBand_ = -1;
        backgroundDensity_ = 0.0; peakShare_ = 0.0; totalBandEnergy_ = 0.0;
        widthReason_ = LowEndReason::NotFinished;
        noteReason_  = LowEndReason::NotFinished;
        grid_.reset();
        xover_.reset();
        frames_.reset();
    }

    static constexpr int latencySamples() noexcept { return 0; }       // a read-only sink

    //==============================================================================
    // READ-ONLY. Law 11 order: malformed -> unprepared -> finished -> nch > maxChannels -> n == 0 -> run.
    // A legal call longer than any maxBlock is consumed WHOLE (nothing here is sized by a block length);
    // a refused call consumes nothing.
    [[nodiscard]] bool process (const float* const* in, int numChannels, int n) noexcept
    {
        if (numChannels < 0 || n < 0) return false;
        if (numChannels > 0 && in == nullptr) return false;            // malformed, before anything else
        if (! prepared_ || finished_) return false;
        if (numChannels > channels_) return false;
        if (n == 0) return true;

        for (int i = 0; i < n; ++i)
        {
            // --- the L/R pair, and whether this sample is a hole ---
            const bool fedL = numChannels > 0;
            const bool fedR = channels_ > 1 ? numChannels > 1 : fedL;    // mono: R IS L, so S is exactly 0
            const float l = fedL ? in[0][i] : 0.0f;
            const float r = channels_ > 1 ? (fedR ? in[1][i] : 0.0f) : l;
            const bool fed = fedL && fedR;
            bool hole = ! fed || ! std::isfinite (l) || ! std::isfinite (r);
            if (! fed) ++absentSamples_;                                // a channel the call did not carry
            else if (hole) ++nonFiniteSamples_;                         // a sample that was there and was not a number

            // --- Mid/Side FIRST, then the crossover. Encoding before filtering keeps the SIDE signal's
            // relative precision: filtering L and R separately would quantise two large correlated
            // histories into float SVF state and then subtract them, manufacturing side residue in
            // near-mono material. The two orders are mathematically identical (one linear filter, both
            // channels) and the choice is purely numerical. MidSide::encode can itself overflow a finite
            // l + r, so its outputs are checked before anything is squared. ---
            float m = 0.0f, s = 0.0f;
            if (! hole)
            {
                stereo::MidSide::encode (l, r, m, s);
                if (! std::isfinite (m) || ! std::isfinite (s)) { hole = true; m = 0.0f; s = 0.0f; }
            }

            float lowM = 0.0f, highM = 0.0f, lowS = 0.0f, highS = 0.0f;
            xover_.processSample (0, m, lowM, highM);
            xover_.processSample (1, s, lowS, highS);
            // A huge finite input could overflow INSIDE the filter. DEFENSIVE: no finite float input
            // found so far reaches it — the SVF's arithmetic is double and only its float state is at
            // risk, and at the largest input that survives MidSide::encode (l = -r = 1.7e38, so
            // s = 1.7e38) every output stays finite. Kept because the response is right if it ever does
            // fire, and checked HERE, on this sample, so the grid's flush below cannot erase a poison
            // this code has not counted. See filterNonFiniteSamples() for the case it still cannot see.
            if (! (std::isfinite (lowM) && std::isfinite (highM) && std::isfinite (lowS) && std::isfinite (highS)))
            {
                ++filterNonFinite_;
                hole = true;
                xover_.healPoison();
                lowM = highM = lowS = highS = 0.0f;
            }

            if (hole)
            {
                ++openHoles_; ++holeSamples_;
                if (firstHole_ < 0) firstHole_ = totalSamples_;
                lastHole_ = totalSamples_;
            }
            else
            {
                // promote BEFORE squaring: a float 3e38 squares to a float infinity
                const double dLowM = (double) lowM, dLowS = (double) lowS;
                const double dHighM = (double) highM, dHighS = (double) highS;
                const double dM = (double) m, dS = (double) s;
                addTo (0, dLowM * dLowM); addTo (1, dLowS * dLowS);
                addTo (2, dHighM * dHighM); addTo (3, dHighS * dHighS);
                addTo (4, dM * dM);        addTo (5, dS * dS);
                openMid_  += dLowM * dLowM;
                openSide_ += dLowS * dLowS;
                ++openFinite_; ++finiteSamples_;
                const double a = std::fabs (dLowS);                   // peak VERTICAL excursion of the low band
                if (a > peakSideAmp_) { peakSideAmp_ = a; peakSideAmpAt_ = totalSamples_; }
            }

            // --- the spectral axes see the RAW, UNFILTERED Mid and Side. The note range reaches 300 Hz,
            // 1.3 octaves above the crossover, so filtering here would erase most of it. ---
            frames_.push (0, m, ! hole);
            frames_.push (1, s, ! hole);

            // --- maintenance on the AUDIO clock, after the sample is complete. Never at the end of a
            // process() call: that boundary is the caller's, not the audio's (StateGrid.h:30). ---
            if (grid_.advance (1)) xover_.flushDenormals();

            ++totalSamples_;

            // A frame closes on the sample that completes it; the block closes on the sample that ends
            // it. When both land on the same sample the FRAME is emitted first — fixed, so a trace is
            // orderable.
            if (frames_.tick()) consumeFrame();
            if (totalSamples_ == nextBlockEnd_) { closeBlock (blockSamples_); nextBlockEnd_ += blockSamples_; }
        }
        return true;
    }

    // End of stream. Closes the one partial block at [floor(T/B)*B, T) and normalises it by its ACTUAL
    // length; emits NO spectral frame (the tail contract) and publishes tailUncoveredSamples() instead;
    // reduces the band table once and freezes. Idempotent.
    [[nodiscard]] bool finish() noexcept
    {
        if (! prepared_) return false;
        if (finished_) return true;                                    // idempotent, and not an error
        const std::int64_t openStart = blockIndex_ * blockSamples_;
        if (totalSamples_ > openStart) closeBlock (totalSamples_ - openStart);   // never a zero-length block
        frames_.finish();
        reduceBands();
        widthReason_ = finiteSamples_ == 0 ? LowEndReason::NoFiniteSamples
                     : (sum (0) + sum (1)) > 0.0 ? LowEndReason::Ok
                     : LowEndReason::NoEnergy;
        finished_ = true;
        return true;
    }

    //==============================================================================
    // --- the report: part 1, the width of the bass ---
    bool isFinished() const noexcept { return finished_; }
    bool isPrepared() const noexcept { return prepared_; }
    int  channels()  const noexcept { return channels_; }
    int  analysedChannels() const noexcept { return channels_ > 1 ? 2 : 1; }   // the L/R pair, and no more
    double sampleRate() const noexcept { return sampleRate_; }
    double crossoverHz() const noexcept { return (double) xover_.frequency(); }
    std::int64_t samplesProcessed() const noexcept { return totalSamples_; }
    std::int64_t blockSamples() const noexcept { return blockSamples_; }

    LowEndReason widthReason() const noexcept { return widthReason_; }
    bool widthValid() const noexcept { return widthReason_ == LowEndReason::Ok; }

    double lowMidEnergy()  const noexcept { return sum (0); }          // raw integrals, in samples^2
    double lowSideEnergy() const noexcept { return sum (1); }
    double highMidEnergy() const noexcept { return sum (2); }
    double highSideEnergy() const noexcept { return sum (3); }
    double rawMidEnergy()  const noexcept { return sum (4); }          // UNFILTERED — the full-band reference,
    double rawSideEnergy() const noexcept { return sum (5); }          // because low+high is an allpass, not the input
    double lowBandEnergy() const noexcept { return sum (0) + sum (1); }

    // S/(M+S) of the low band over the whole programme: the ENERGY-WEIGHTED answer. 0 for mono; 0 for
    // silence too, which is why widthReason() and not this number tells the two apart.
    double lowSideFraction() const noexcept { const double t = sum (0) + sum (1); return t > 0.0 ? sum (1) / t : 0.0; }

    // THE LOW BAND'S SHARE OF THE WHOLE PROGRAMME — (lowMid + lowSide) / (rawMid + rawSide). All four are
    // published already and the division is trivial; it is named here so that a caller neither recomposes
    // it from four scalars nor, far more likely, forgets what the numerator actually contains.
    //
    // IT IS NOT "THE ENERGY BELOW crossoverHz", AND THE DIFFERENCE IS LARGE ENOUGH TO INVERT A DECISION.
    // The low band is an LR4, not a brick wall, so this is the true share WEIGHTED BY THAT FILTER'S POWER
    // RESPONSE — and a fourth-order low-pass is -6 dB in power at its own corner, so content AT the
    // crossover counts at a QUARTER. Measured against a fixture whose share is exact by construction (a
    // sine of amplitude a carries a^2/2, so 3.000 % is arithmetic rather than a measurement), at 48 kHz
    // with fc = 30 Hz, and agreeing with the analytic response to four decimals at every row:
    //
    //     f / fc     0.17   0.33   0.40   0.60   0.67   0.83   1.00   1.33   2.00
    //     counted   0.998  0.976  0.950  0.784  0.697  0.455  0.250  0.058  0.0035
    //
    // The ratios are in units of fc and hold at any crossover: the response depends on f/fc alone this far
    // below Nyquist. So a real 3 % of infra energy reads as 2.85 % if it sits at 12 Hz and as 0.75 % if it
    // sits at 30 — the same physical fact, four times apart. A threshold carried over from a different
    // definition of "infra-low" will be wrong; one calibrated on real material against THIS number will
    // not, because the number is stable to four decimals against its own analytic form.
    //
    // WHAT IT IS FOR, now that the band table reaches 20 Hz: the table answers 20..300 Hz per semitone,
    // with duty and levelWhenOnDb, so the octave from E0 to B0 needs no filter at all. This ratio is the
    // right instrument for what the table CANNOT see — strictly below lowNoteHz — and the two together say
    // more than either alone: a share from the filter, and named notes with their occupancy from the bands.
    //
    // 0.0 for the 0/0 of a silent or unmeasured programme, like every other fraction here.
    double infraLowShare() const noexcept
    {
        const double lo = sum (0) + sum (1), raw = sum (4) + sum (5);
        return raw > 0.0 ? lo / raw : 0.0;
    }
    double highSideFraction() const noexcept { const double t = sum (2) + sum (3); return t > 0.0 ? sum (3) / t : 0.0; }
    double rawSideFraction() const noexcept { const double t = sum (4) + sum (5); return t > 0.0 ? sum (5) / t : 0.0; }

    std::int64_t finiteSamples() const noexcept { return finiteSamples_; }
    std::int64_t holeSamples()   const noexcept { return holeSamples_; }
    std::int64_t nonFiniteSamples() const noexcept { return nonFiniteSamples_; }
    std::int64_t absentSamples() const noexcept { return absentSamples_; }   // the L/R pair was not fully fed
    // A LOWER BOUND, not an exact count, and the mechanism is worth naming. eq::Svf computes in double
    // but stores its integrator state in FLOAT, and that state is updated as `2*v - ic`: a v just under
    // FLT_MAX therefore overflows the STATE while the output this code checks stays finite. The next
    // sample would expose it — an infinite state makes the next output non-finite, which is counted here
    // and healed — unless a StateGrid boundary falls in between, because flushDenormals() heals poison
    // as well as denormals (Svf.h:161). Measured reachable: at fs 1000, fc 400, mono +-1.7e38 at samples
    // 61..63 with the boundary at 64 — a rate below the core's floor since P51, and the same case at fs 8000,
    // fc 3200: the filter sees only tan(pi fc / fs), the same double for both pairs, and the grid period is 64
    // samples at any rate (P51 checked it: the two crossovers' outputs over that burst are bit-identical, and this
    // counter reads 0 at 8 kHz). NO MEASURED VALUE IS WRONG when that happens — the state is healed before any
    // sample consumes it, and every published energy stays finite — but this counter reads 0.
    // Detecting it properly needs the filter's state, which this instrument does not own.
    std::int64_t filterNonFiniteSamples() const noexcept { return filterNonFinite_; }
    // A hole feeds the documented canonical zero to the filters, so the LR4 state carries it for its
    // own settling time: blocks after lastHoleSample() are finite but not untouched. The coordinates are
    // published rather than a settling constant invented, so a consumer discounts exactly what it likes.
    std::int64_t firstHoleSample() const noexcept { return firstHole_; }
    std::int64_t lastHoleSample()  const noexcept { return lastHole_; }

    // --- the 10 ms series. Capacity exhaustion is DATA: the prefix is kept, the count keeps counting. ---
    std::int64_t blockCount() const noexcept { return blockCount_; }
    std::int64_t storedBlockCount() const noexcept { return std::min<std::int64_t> (blockCount_, (std::int64_t) blocks_.size()); }
    bool blocksComplete() const noexcept { return blockCount_ <= (std::int64_t) blocks_.size(); }
    LowEndBlock block (std::int64_t i) const noexcept
    {
        LowEndBlock b;
        if (i < 0 || i >= storedBlockCount()) return b;
        const BlockRecord& r = blocks_[(std::size_t) i];
        b.index = i;
        b.start = i * blockSamples_;
        b.samples = (std::int64_t) r.samples;
        b.finiteSamples = (std::int64_t) r.finite;
        b.holes = (std::int64_t) r.samples - (std::int64_t) r.finite;
        b.midEnergy = r.mid;
        b.sideEnergy = r.side;
        b.valid = r.finite > 0 && r.finite == r.samples;
        return b;
    }

    // The DURATION-weighted distribution of the block side fraction: fixed edges, bin j is
    // [j/100, (j+1)/100) and 1.0 lands in the last bin. Counted in SAMPLES, not blocks, so the short
    // final block weighs what it actually lasts; blocks with no energy are excluded (a silent block is
    // not a mono block) and histogramSamples() says how much was covered. It accumulates from block
    // zero and never stops, so it does not change with maxBlocks.
    std::int64_t histogram (int bin) const noexcept { return bin >= 0 && bin < kHistogramBins ? hist_[bin] : 0; }
    // Blocks the histogram did NOT count, because `skipBlocks` asked. Published rather than implied: a
    // population whose size a caller cannot see is a population it cannot divide by.
    std::int64_t skippedBlocks() const noexcept { return skippedBlocks_; }

    // HOW MANY 10 ms BLOCKS AN LR4 AT `crossoverHz` NEEDS TO FALL `dB` BELOW ITS OWN PEAK — the number
    // `skipBlocks` wants, derived from the filter rather than written down. Pure, static, allocation-free:
    // a caller computes it before prepare() and passes the answer back in.
    //
    // THE MODEL, stated so it can be argued with. LR4 is two cascaded Butterworth sections, so its poles
    // are a DOUBLE pair at real part -wc/sqrt(2): the transient envelope is t·exp(-t/tau) with
    // tau = sqrt(2)/(2·pi·fc). Normalised to its own peak (at t = tau) that is u·exp(1-u) with u = t/tau,
    // and the answer is the u where it reaches 10^(-dB/20). Solved here by bisection, which is
    // deterministic, needs no libm beyond an exponential and cannot diverge.
    //
    // u = 10.233 at -60 dB and 17.688 at -120, so at 48 kHz: 12 blocks at 20 Hz and 2 at 120 Hz for -60 dB.
    // The "ten blocks" this replaces was five times too many at 120 Hz and too few at 20.
    //
    // It is an ENVELOPE bound, not a promise about a particular programme: a filter driven by music is
    // never at its own impulse peak, so this is the pessimistic end. 0 for arguments the crossover itself
    // would refuse, and for a non-positive dB — a caller asking to skip nothing gets nothing skipped.
    static int settlingBlocks (double sampleRate, double crossoverHz, double dB) noexcept
    {
        if (! (sampleRate >= kMinSampleRate && sampleRate <= kMaxSampleRate)) return 0;
        if (! (crossoverHz >= kMinCrossoverHz) || ! (crossoverHz <= 0.49 * sampleRate)) return 0;
        if (! (dB > 0.0) || ! (dB <= 400.0)) return 0;
        const double target = core::det::pow10 (-dB / 20.0);
        double lo = 1.0, hi = 400.0;
        for (int i = 0; i < 200; ++i)                       // fixed count: the same work on every input
        {
            const double u = 0.5 * (lo + hi);
            const double v = u * core::det::exp2 ((1.0 - u) * 1.4426950408889634);   // exp(1-u)
            if (v > target) lo = u; else hi = u;
        }
        const double tau = 1.4142135623730951 / (2.0 * core::kPi * crossoverHz);
        const double blockSec = (double) std::max<std::int64_t> (1, (std::int64_t) std::lround (0.01 * sampleRate)) / sampleRate;
        const double blocks = std::ceil (0.5 * (lo + hi) * tau / blockSec);
        return (blocks >= 0.0 && blocks < 1.0e9) ? (int) blocks : 0;
    }
    std::int64_t histogramSamples() const noexcept { return histSamples_; }

    // The coordinates. Three different questions, three different blocks — a block can win any one of
    // them alone. Ties go to the EARLIEST block. -1 when no block qualified.
    std::int64_t worstFractionBlock() const noexcept { return worstFracBlock_; }
    double worstFraction() const noexcept { return worstFrac_ > 0.0 ? worstFrac_ : 0.0; }
    double worstFractionEnergy() const noexcept { return worstFracEnergy_; }     // weigh an accidental 1.0
    std::int64_t peakEnergyBlock() const noexcept { return peakEnergyBlock_; }
    double peakBlockEnergy() const noexcept { return peakEnergy_ > 0.0 ? peakEnergy_ : 0.0; }
    double peakEnergyBlockFraction() const noexcept { return peakEnergyFrac_; }
    // The greatest VERTICAL modulation, which is the quantity a cutting engineer asks for first and
    // which neither extremum above identifies.
    std::int64_t peakSideEnergyBlock() const noexcept { return peakSideBlock_; }
    double peakBlockSideEnergy() const noexcept { return peakSide_ > 0.0 ? peakSide_ : 0.0; }
    double peakLowSideAmplitude() const noexcept { return peakSideAmp_; }
    std::int64_t peakLowSideAmplitudeAt() const noexcept { return peakSideAmpAt_; }

    //==============================================================================
    // --- the report: part 2, the dominant low note ---
    LowEndReason noteReason() const noexcept { return noteReason_; }
    bool noteValid() const noexcept { return noteReason_ == LowEndReason::Ok; }
    int bandCount() const noexcept { return bandCount_; }
    LowEndBand band (int b) const noexcept { return b >= 0 && b < bandCount_ ? bands_[(std::size_t) b] : LowEndBand {}; }
    int underResolvedBands() const noexcept                    // bands narrower than a Hann main lobe
    {
        int k = 0;
        for (int b = 0; b < bandCount_; ++b) if (bands_[(std::size_t) b].binsPerBand < (double) kLobeBins) ++k;
        return k;
    }

    // WHERE THE UNDER-RESOLVED BANDS STOP. A count without a boundary is not actionable: a caller reading
    // `underResolvedBands() == 4` cannot tell WHICH four without walking the table. Band width grows
    // monotonically with centre frequency (widthHz = c·(kSemiUp - kSemiDown)), so the under-resolved set is
    // always a PREFIX and one index names it.
    //
    // THIS IS NOT A K9 PROBLEM, and that is the reason it is published rather than left implicit. At 96 kHz
    // the bin is 0.7324 Hz, and the OLD 30 Hz table already had NINE under-resolved bands, up to 49.0 Hz —
    // true of every 96 kHz file this class has ever measured, with nothing saying so.
    int firstResolvedBand() const noexcept { const int k = underResolvedBands(); return k < bandCount_ ? k : -1; }

    // The centre frequency at and above which a band is at least kLobeBins wide — the closed form of the
    // same criterion, so it is defined even where the table does not reach:
    //     resolvedAboveHz = kLobeBins·binHz / (kSemiUp - kSemiDown)
    // 25.356 Hz at 48 kHz / order 17, 23.298 at 44.1 kHz, 50.712 at 96 kHz, 12.678 at 48 kHz / order 18.
    // It is a NOMINAL threshold on the band CENTRE under the four-bin criterion, not a promise about a tone
    // near a band edge: a tone 50 cents off splits roughly 50/50 between two bands at any resolution.
    // 0.0 before prepare(), like every other geometry reading here.
    double resolvedAboveHz() const noexcept
    {
        return binHz_ > 0.0 ? (double) kLobeBins * binHz_ / (kSemiUp - kSemiDown) : 0.0;
    }
    static constexpr int lobeBins() noexcept { return kLobeBins; }   // so a caller never re-types the 4

    //==============================================================================
    // --- the report: part 2b, OCCUPANCY (K9). How OFTEN a band is there, not how loud it is on average ---
    //
    // WHY THE INTEGRAL IS NOT ENOUGH, which is the whole reason this exists. A sub that plays on an eighth
    // of the programme is 10·log10(1/8) = 9.03 dB down in the integral against one that plays throughout,
    // and a threshold set for the second misses the first. Duty separates "how loud when present" from
    // "how often present", and the band table already answers the first.
    //
    // DUTY IS A COUNT OF FRAMES, NOT A FRACTION OF TIME, and a caller must calibrate in those units. The
    // frames overlap — at the default hop each sample is in two — so one event can mark two frames, and
    // the quantum is 1/dutyFrames(). At 48 kHz and order 17 the hop is 1.365 s, so a three-minute
    // programme offers 130 frames: an event on "an eighth of the time" does NOT read 0.125.
    std::int64_t dutyFrames() const noexcept { return dutyFrames_; }          // the denominator, published
    std::int64_t dutyCount (int b) const noexcept                             // the numerator, published
    {
        return b >= 0 && b < bandCount_ ? dutyCount_[(std::size_t) b] : 0;
    }
    // In [0, 1]. Canonically 0 when no frame was counted — and then noteReason() says which kind of
    // nothing it was (ShorterThanWindow, NoUsableFrames, NoEnergy), so a caller is never left to guess.
    double duty (int b) const noexcept
    {
        return dutyFrames_ > 0 ? (double) dutyCount (b) / (double) dutyFrames_ : 0.0;
    }
    // HOW LOUD THE BAND IS WHEN IT IS ON, against the loudest band of the same frame, in dB — the mean
    // over the frames that counted it. In [-dutyThresholdDb, 0] by construction. This is the number duty
    // does NOT carry: two bands with the same duty can sit at 0 dB ("it IS the bass when it plays") and at
    // -19 dB ("it just scraped in"), and a consumer choosing a filter needs to tell those apart. 0.0 when
    // the band was never counted.
    //
    // The ratio is accumulated LINEARLY and the logarithm is taken here, once, on read: a dB is a view.
    double levelWhenOnDb (int b) const noexcept
    {
        const std::int64_t k = dutyCount (b);
        return k > 0 ? 10.0 * core::det::log10 (dutyRatio_[(std::size_t) b] / (double) k) : 0.0;
    }
    // …and the same number measured from the threshold rather than from the frame's peak, which is what a
    // caller asking "how much room did it have" wants. Zero means it sat exactly on the line.
    double marginWhenOnDb (int b) const noexcept
    {
        return dutyCount (b) > 0 ? levelWhenOnDb (b) + params_.dutyThresholdDb : 0.0;
    }
    double dutyThresholdDb() const noexcept { return params_.dutyThresholdDb; }   // the installed value

    // The lowest band present in at least `dutyMin` of the counted frames, with everything a caller needs
    // to act on it. `band == -1` is the canonical none — this file's own convention (an absent index is
    // -1, never a sentinel frequency) — and it covers both "no band qualified" and "nothing was measured";
    // dutyFrames() and noteReason() tell those apart without a second enum.
    struct LowestOccupied
    {
        int          band  = -1;
        int          midi  = 0;
        double       centreHz = 0.0;
        std::int64_t count = 0;
        double       duty  = 0.0;
        double       levelWhenOnDb = 0.0;
        double       marginWhenOnDb = 0.0;
    };
    // An ACCESSOR, taking the threshold, because the consumer sweeps it while calibrating and a field
    // frozen at finish() would force a re-run of the audio for every hypothesis. 47 comparisons.
    //
    // `count > 0` IS LOAD-BEARING: with dutyMin <= 0 the test `duty >= dutyMin` is true of a band present
    // in no frame at all, and the call would return band 0 "occupied 0 % of the time" — a number that
    // reads as a finding. A non-finite dutyMin fails every comparison and returns none, by the same rule.
    //
    // ±1 BAND BY CONSTRUCTION, and a caller must expect it: a tone 40 cents below a band centre leaves
    // about a third of its power in the band beneath, only ~5 dB down, which a 20 dB threshold admits. The
    // returned band's `centsOffset` in the table is the tell — near +50 cents means the energy sits at the
    // band's top edge and belongs to the note above.
    // THE SIDE FRACTION BELOW ANY CANDIDATE CROSSOVER, from the band table alone — no second pass over the
    // audio, no second filter, no extra state. The spectral axes are fed the RAW Mid and Side, BEFORE the
    // crossover (see process()), so the table is an unfiltered semitone spectrum of both axes and any LR4
    // low-pass can be applied to it afterwards as a weight:
    //
    //     sideFractionBelow(fc) = SUM_b |H(c_b;fc)|^2 · side_b  /  SUM_b |H(c_b;fc)|^2 · (mid_b + side_b)
    //
    // with the SAME response the real filter has — |H|^2 = 1/(1+r^4)^2, r = tan(pi·f/fs)/tan(pi·fc/fs).
    // The prewarped ratio is not optional: f/fc is the analogue prototype and is wrong near Nyquist, which
    // the header above says at length. `det::tan` for the same reason `det::pow10` is used elsewhere — a
    // published number must not depend on which libm the row was built against.
    //
    // WHAT IT IS AND IS NOT, MEASURED BOTH WAYS. It is a SWEEP instrument: it answers "where should the
    // crossover sit" for a whole curve at one call per point, which is what choosing a mono-bass frequency
    // actually asks. It is NOT the number the installed crossover produces.
    //
    // ON MATERIAL WHOSE LOW END LIVES INSIDE THE TABLE it is very nearly the same number: against a real
    // LR4 on mono 41 Hz plus anti-phase 98 and 220 Hz, the two agree to 4.7e-5 absolute at every crossover
    // from 60 to 300 Hz — far closer than the arithmetic needs to be for a question posed at 0.25.
    //
    // OUTSIDE THE TABLE IT IS BLIND, AND THE FAILURE IS NOT SMALL. The bands cover [lowNoteHz, highNoteHz]
    // and nothing else, so anti-phase energy under 20 Hz or above 300 Hz is absent from both sums while
    // the real low-pass passes it. Add anti-phase 15 Hz and 700 Hz to that fixture and at fc = 60 Hz this
    // function reads 0.0036 where the installed filter reads 0.488 — the opposite answer to a question
    // asked at 0.25, and the 15 Hz term is precisely the vertical hazard a lacquer cares about. Two more
    // differences are small beside that one and are named for completeness: the population is the USED
    // FRAMES, window-weighted, not every sample; and no filter start-up is carried.
    //
    // SO: sweep to CHOOSE a frequency, then install it and read lowSideFraction() for the answer. A caller
    // that wants one number and not a curve should not be here.
    //
    // 0.0 FOR THE 0/0 of a silent or unmeasured programme — the convention LowEndBand::sideFraction() uses.
    // But -1.0 FOR A FREQUENCY THIS CLASS WOULD REFUSE, and the difference matters: 0.0 is a legitimate
    // reading (a perfectly mono low end) and a caller that passed a bad frequency would have read it as
    // one. A fraction is in [0, 1], so a negative is unmistakably not an answer.
    double sideFractionBelow (double fc) const noexcept
    {
        if (! (fc >= kMinCrossoverHz) || ! (fc <= 0.49 * sampleRate_)) return -1.0;
        const double tc = core::det::tan (core::kPi * fc / sampleRate_);
        if (! (tc > 0.0) || ! std::isfinite (tc)) return -1.0;
        double num = 0.0, den = 0.0;
        for (int b = 0; b < bandCount_; ++b)
        {
            const LowEndBand& row = bands_[(std::size_t) b];
            const double r  = core::det::tan (core::kPi * row.centreHz / sampleRate_) / tc;
            const double r4 = r * r * r * r;
            const double d  = 1.0 + r4;
            if (! std::isfinite (d) || ! (d > 0.0)) continue;          // far above fc: weight is 0, skip it
            const double w  = 1.0 / (d * d);
            if (! std::isfinite (w)) continue;
            num += w * row.sideEnergy;
            den += w * (row.midEnergy + row.sideEnergy);
        }
        return (std::isfinite (num) && std::isfinite (den) && den > 0.0) ? num / den : 0.0;
    }

    LowestOccupied lowestOccupiedBand (double dutyMin) const noexcept
    {
        LowestOccupied r;
        for (int b = 0; b < bandCount_; ++b)
        {
            if (! (dutyCount_[(std::size_t) b] > 0) || ! (duty (b) >= dutyMin)) continue;
            r.band = b;
            r.midi = bands_[(std::size_t) b].midi;
            r.centreHz = bands_[(std::size_t) b].centreHz;
            r.count = dutyCount_[(std::size_t) b];
            r.duty  = duty (b);
            r.levelWhenOnDb  = levelWhenOnDb (b);
            r.marginWhenOnDb = marginWhenOnDb (b);
            break;
        }
        return r;
    }

    int peakBand() const noexcept { return peakBand_; }                  // argmax of ENERGY — what a tone does
    int peakDensityBand() const noexcept { return peakDensityBand_; }    // argmax of DENSITY — what noise does
    int secondBand() const noexcept { return secondBand_; }              // the runner-up by energy
    double backgroundDensity() const noexcept { return backgroundDensity_; }   // median density of the non-peak bands
    double peakBandEnergy() const noexcept { return peakBand_ >= 0 ? bands_[(std::size_t) peakBand_].energy : 0.0; }
    double peakBandWidthHz() const noexcept { return peakBand_ >= 0 ? bands_[(std::size_t) peakBand_].widthHz : 0.0; }
    double secondBandEnergy() const noexcept { return secondBand_ >= 0 ? bands_[(std::size_t) secondBand_].energy : 0.0; }
    double totalBandEnergy() const noexcept { return totalBandEnergy_; }
    // The frame's own window-weighted mean square, sum(w*x)^2/sum(w^2) by Parseval, folded over every
    // bin and averaged over the used frames. The denominator totalBandEnergy() is a share OF.
    double frameEnergy() const noexcept { return frameEnergy_; }
    // totalBandEnergy() / frameEnergy(): how much of the frame lies inside [lowNoteHz, highNoteHz].
    // On a FLAT spectrum this equals the fraction of Nyquist the bands cover (1.1345 % at the defaults
    // and 48 kHz), which is what makes it a usable sanity number rather than an abstraction.
    double bandRangeShare() const noexcept { return bandRangeShare_; }
    double peakShare() const noexcept { return peakShare_; }             // peak / total, in [0, 1]: pole-free
    // The dominance ratio is deliberately NOT a field: it is
    //   peakBandEnergy() / (backgroundDensity() * peakBandWidthHz()),
    // and the denominator is EXACTLY zero for a synthetic tone in digital silence. The three numbers are
    // published so a consumer divides them with its own eyes open, exactly as with S/M above.
    int    peakMidi()     const noexcept { return peakBand_ >= 0 ? bands_[(std::size_t) peakBand_].midi : 0; }
    double peakNoteHz()   const noexcept { return peakBand_ >= 0 ? bands_[(std::size_t) peakBand_].centreHz : 0.0; }
    double peakCentroidHz() const noexcept { return peakBand_ >= 0 ? bands_[(std::size_t) peakBand_].centroidHz : 0.0; }
    double peakCentsOffset() const noexcept { return peakBand_ >= 0 ? bands_[(std::size_t) peakBand_].centsOffset : 0.0; }
    // Is the dominant note CUTTABLE? S/(M+S) of the peak band alone.
    double peakBandSideFraction() const noexcept { return peakBand_ >= 0 ? bands_[(std::size_t) peakBand_].sideFraction() : 0.0; }

    std::int64_t usedFrames()  const noexcept { return usedFrames_; }
    std::int64_t holedFrames() const noexcept { return holedFrames_; }
    std::int64_t tailUncoveredSamples() const noexcept { return frames_.tailUncoveredSamples(); }
    std::int64_t windowSamples() const noexcept { return frames_.windowSamples(); }
    std::int64_t hopSamples()    const noexcept { return frames_.hopSamples(); }
    double binHz() const noexcept { return binHz_; }

    // The 12 pitch classes and the octave, so a caller prints "E2" without this header owning a string.
    // Scientific pitch: MIDI 60 is C4.
    static const char* pitchClassName (int midi) noexcept
    {
        static const char* const names[12] = { "C", "C#", "D", "D#", "E", "F", "F#", "G", "G#", "A", "A#", "B" };
        const int pc = ((midi % 12) + 12) % 12;
        return names[pc];
    }
    static int noteOctave (int midi) noexcept { return (int) std::floor ((double) midi / 12.0) - 1; }

using CrossoverType = eq::DeterministicCrossover2;   // asserted by the math-policy suite

private:
    static SpectrumFramesParams framesParams (const LowEndParams& p) noexcept
    {
        SpectrumFramesParams fp;
        fp.fftOrder = p.fftOrder;
        fp.hop = p.hop;
        return fp;
    }

    static std::int64_t blockSamplesFor (double sampleRate) noexcept
    {
        const double b = kBlockMs / 1000.0 * sampleRate;               // 0.01*fs, the shared 10 ms grid
        if (! (b >= 1.0) || ! (b < 1.0e15)) return 0;
        return (std::int64_t) std::lround (b);
    }

    static double noteHz (const LowEndParams& p, int midi) noexcept
    {
        return p.tuningHz * core::det::exp2 ((double) (midi - 69) / 12.0);
    }

    // The band set: every MIDI note whose CENTRE lies inside [lowNoteHz, highNoteHz]. Computed from the
    // parameters alone, once, so it is not a function of the data or of the slicing.
    static bool midiRange (const LowEndParams& p, int& loMidi, int& hiMidi) noexcept
    {
        const double a = 69.0 + 12.0 * core::det::log2 (p.lowNoteHz  / p.tuningHz);
        const double b = 69.0 + 12.0 * core::det::log2 (p.highNoteHz / p.tuningHz);
        if (! std::isfinite (a) || ! std::isfinite (b)) return false;
        if (! (a > -2000.0 && b < 2000.0)) return false;
        loMidi = (int) std::ceil (a);
        hiMidi = (int) std::floor (b);
        // The bound came through a logarithm, so a `high` one ulp below a note's centre rounds to that
        // centre and would admit a note ABOVE the range (measured: highNoteHz = nextafter(440, 0) still
        // included 440 Hz). Snap both ends against the SAME noteHz() that defines a band, so the rule
        // "every note whose centre lies in [low, high]" is exact by construction rather than to within
        // a rounding of its own bound.
        if (loMidi <= hiMidi && noteHz (p, loMidi) < p.lowNoteHz)  ++loMidi;
        if (loMidi <= hiMidi && noteHz (p, hiMidi) > p.highNoteHz) --hiMidi;
        return hiMidi >= loMidi;
    }

    // Which bins a band touches. Bin k's cell is [(k-1/2)*binHz, (k+1/2)*binHz], clipped to [0, fs/2],
    // so DC and Nyquist are half cells — that IS their one-sided weight.
    static void bandBins (double centreHz, double binHz, int bins, int& first, int& count) noexcept
    {
        first = 0; count = 0;
        if (! (binHz > 0.0) || bins < 1) return;
        const double lo = centreHz * kSemiDown, hi = centreHz * kSemiUp;
        int a = (int) std::floor (lo / binHz + 0.5);
        int b = (int) std::floor (hi / binHz + 0.5);
        if (a < 0) a = 0;
        if (b > bins - 1) b = bins - 1;
        if (b < a) return;
        first = a; count = b - a + 1;
    }

    [[nodiscard]] bool buildBands() noexcept
    {
        const int bins = frames_.bins();
        std::size_t at = 0;
        for (int b = 0; b < bandCount_; ++b)
        {
            const int midi = midiLo_ + b;
            const double c = noteHz (params_, midi);
            const double lo = c * kSemiDown, hi = c * kSemiUp;
            int first = 0, count = 0;
            bandBins (c, binHz_, bins, first, count);
            bandCountBins_[(std::size_t) b] = count;
            if (at + (std::size_t) count > weights_.size()) return false;   // storageFor() and this must agree
            for (int j = 0; j < count; ++j)
            {
                const int k = first + j;
                double cellLo = ((double) k - 0.5) * binHz_;
                double cellHi = ((double) k + 0.5) * binHz_;
                if (cellLo < 0.0) cellLo = 0.0;                        // bin 0 is a half cell
                const double nyq = 0.5 * sampleRate_;
                if (cellHi > nyq) cellHi = nyq;                        // and so is Nyquist
                const double ovLo = std::max (cellLo, lo), ovHi = std::min (cellHi, hi);
                const double cell = cellHi - cellLo;
                const double ov = ovHi > ovLo ? ovHi - ovLo : 0.0;
                // FOLD: a real tone's power splits between bin k and bin N-k, and only k is exposed, so
                // an interior bin counts twice. DC and Nyquist have no mirror and count once. With this
                // a band holding a full-scale sine of amplitude A reads A^2/2 — a true mean square.
                const double fold = (k == 0 || k == bins - 1) ? 1.0 : 2.0;
                BinWeight& w = weights_[at + (std::size_t) j];
                w.bin = k;
                w.weight = cell > 0.0 ? fold * ov / cell : 0.0;
                // the FIRST MOMENT of a partial cell sits at the midpoint of the OVERLAP, not at the
                // bin's centre — the correction that keeps a centroid inside its own band
                w.momentHz = ov > 0.0 ? 0.5 * (ovLo + ovHi) : 0.0;
            }
            at += (std::size_t) count;

            LowEndBand& row = bands_[(std::size_t) b];
            row.midi = midi;
            row.centreHz = c;
            row.widthHz = hi - lo;
            row.binsPerBand = binHz_ > 0.0 ? (hi - lo) / binHz_ : 0.0;
        }
        return at == weights_.size();
    }

    // Neumaier: two constant scalars per quantity, which law 7 permits explicitly and which
    // analysis::ClipDetector already uses for its DC sum (ClipDetector.h:393).
    void addTo (int i, double x) noexcept
    {
        const double s = sum_[i] + x;
        comp_[i] += std::fabs (sum_[i]) >= std::fabs (x) ? (sum_[i] - s) + x : (x - s) + sum_[i];
        sum_[i] = s;
    }
    double sum (int i) const noexcept { return sum_[i] + comp_[i]; }

    void closeBlock (std::int64_t samples) noexcept
    {
        const std::int64_t start = blockIndex_ * blockSamples_;
        if (blockCount_ < (std::int64_t) blocks_.size())
        {
            BlockRecord& r = blocks_[(std::size_t) blockCount_];
            r.mid = openMid_; r.side = openSide_;
            r.samples = (std::int32_t) samples;
            r.finite = (std::int32_t) openFinite_;
        }
        ++blockCount_;

        const double energy = openMid_ + openSide_;
        const double frac = energy > 0.0 ? openSide_ / energy : 0.0;
        if (energy > 0.0 && openFinite_ > 0)
        {
            int bin = (int) (frac * (double) kHistogramBins);
            if (bin < 0) bin = 0;
            if (bin >= kHistogramBins) bin = kHistogramBins - 1;
            if (blockIndex_ >= (std::int64_t) params_.skipBlocks)        // the leading blocks a caller asked to drop
            {
                hist_[bin] += openFinite_;                              // duration-weighted: samples, not blocks
                histSamples_ += openFinite_;
            }
            else ++skippedBlocks_;
            if (frac > worstFrac_) { worstFrac_ = frac; worstFracBlock_ = blockIndex_; worstFracEnergy_ = energy; }
            if (energy > peakEnergy_) { peakEnergy_ = energy; peakEnergyBlock_ = blockIndex_; peakEnergyFrac_ = frac; }
            if (openSide_ > peakSide_) { peakSide_ = openSide_; peakSideBlock_ = blockIndex_; }
        }

        if (trace_ != nullptr)
        {
            LowEndTrace t;
            t.kind = LowEndTrace::Kind::Block;
            t.index = blockIndex_;
            t.start = start;
            t.end = start + samples;
            t.finiteSamples = openFinite_;
            t.holes = openHoles_;
            t.valid = openFinite_ > 0 && openHoles_ == 0;
            t.midEnergy = openMid_;
            t.sideEnergy = openSide_;
            trace_ (traceUser_, t);
        }

        ++blockIndex_;
        openMid_ = openSide_ = 0.0;
        openFinite_ = openHoles_ = 0;
    }

    // One closed frame. A frame that holds a hole in EITHER axis is discarded whole and counted, so
    // every accumulated frame carries the same content.
    void consumeFrame() noexcept
    {
        const bool usable = frames_.frameFinite (0) && frames_.frameFinite (1);
        const double* pm = frames_.power (0);
        const double* ps = frames_.power (1);
        if (! usable || pm == nullptr || ps == nullptr)
        {
            ++holedFrames_;
            if (trace_ != nullptr) fireFrameTrace (false);
            return;
        }
        // The frame's OWN window-weighted mean square, folded over every bin — the honest denominator
        // for "how much of this frame is in the note range", and the number that says a range holding
        // only round-off holds nothing.
        double frameTotal = 0.0;
        {
            const int bins = frames_.bins();
            for (int k = 0; k < bins; ++k)
            {
                const double fold = (k == 0 || k == bins - 1) ? 1.0 : 2.0;
                frameTotal += (pm[k] + ps[k]) * fold;
            }
        }
        std::size_t at = 0;
        for (int b = 0; b < bandCount_; ++b)                            // bands ascending, bins ascending
        {
            const int count = bandCountBins_[(std::size_t) b];
            double em = 0.0, es = 0.0, mom = 0.0;
            for (int j = 0; j < count; ++j)
            {
                const BinWeight& w = weights_[at + (std::size_t) j];
                const double p = (pm[w.bin] + ps[w.bin]) * w.weight;
                em += pm[w.bin] * w.weight;
                es += ps[w.bin] * w.weight;
                mom += p * w.momentHz;
            }
            at += (std::size_t) count;
            if (! (std::isfinite (em) && std::isfinite (es) && std::isfinite (mom)))
            {
                ++holedFrames_;                                         // an overflowed transform is a hole too
                if (trace_ != nullptr) fireFrameTrace (false);
                return;
            }
            traceMid_[(std::size_t) b] = em;
            traceSide_[(std::size_t) b] = es;
            accMid_[(std::size_t) b] += em;
            accSide_[(std::size_t) b] += es;
            accMoment_[(std::size_t) b] += mom;
        }
        if (! std::isfinite (frameTotal)) { ++holedFrames_; if (trace_ != nullptr) fireFrameTrace (false); return; }
        accFrameEnergy_ += frameTotal;
        frameTotal_ = frameTotal;
        ++usedFrames_;
        countDuty (frameTotal);
        if (trace_ != nullptr) fireFrameTrace (true);
    }


    // DUTY: IN HOW MANY FRAMES WAS THIS BAND PRESENT. Called once per ACCEPTED frame, from the one place a
    // frame is accepted, so its clock is the frame schedule and law 8a holds by construction. The trace
    // arrays it reads are written unconditionally — the observer only decides whether they are REPORTED —
    // so duty costs a caller nothing and does not depend on one being installed.
    //
    // THE GATE IS THE WHOLE DESIGN. "Within dutyThresholdDb of the frame's loudest band" is `E_b >= max·q`,
    // and on a frame whose every band is zero that reads `0 >= 0` — TRUE for all of them. A silent
    // programme would report every band occupied in every frame, which is not a small error: it is the
    // exact opposite of the answer. A frame whose bands hold only the transform's own round-off is the
    // same failure wearing a number. So a frame enters the denominator only if its band total clears the
    // note floor, the same share of frameEnergy() the whole-programme report uses — one rule, two places.
    void countDuty (double frameTotal) noexcept
    {
        double bandTotal = 0.0, frameMax = 0.0;
        for (int b = 0; b < bandCount_; ++b)
        {
            const double e = traceMid_[(std::size_t) b] + traceSide_[(std::size_t) b];
            bandTotal += e;
            if (e > frameMax) frameMax = e;
        }
        if (! (frameMax > 0.0) || ! (bandTotal > kNoteFloorShare * frameTotal)) return;
        ++dutyFrames_;
        const double floorE = frameMax * dutyShare_;
        for (int b = 0; b < bandCount_; ++b)
        {
            const double e = traceMid_[(std::size_t) b] + traceSide_[(std::size_t) b];
            if (e >= floorE)
            {
                ++dutyCount_[(std::size_t) b];
                dutyRatio_[(std::size_t) b] += e / frameMax;   // linear; the dB is a VIEW, taken on read
            }
        }
    }

    void fireFrameTrace (bool used) noexcept
    {
        LowEndTrace t;
        t.kind = LowEndTrace::Kind::Frame;
        t.index = frames_.frameIndex();
        t.start = frames_.frameStart();
        t.end = frames_.frameEnd();
        t.valid = used;
        t.bandCount = bandCount_;
        t.frameEnergy = used ? frameTotal_ : 0.0;
        t.bandMid = used ? traceMid_.data() : nullptr;
        t.bandSide = used ? traceSide_.data() : nullptr;
        trace_ (traceUser_, t);
    }

    // Once, in finish(): the linear mean over the used frames, then the reductions. Power is averaged
    // LINEARLY and nothing is turned into dB in here at all.
    void reduceBands() noexcept
    {
        if (frames_.frameCount() <= 0) { noteReason_ = LowEndReason::ShorterThanWindow; return; }
        if (usedFrames_ <= 0)          { noteReason_ = LowEndReason::NoUsableFrames;    return; }
        const double inv = 1.0 / (double) usedFrames_;
        totalBandEnergy_ = 0.0;
        bool bad = false;
        for (int b = 0; b < bandCount_; ++b)
        {
            LowEndBand& row = bands_[(std::size_t) b];
            row.midEnergy = accMid_[(std::size_t) b] * inv;
            row.sideEnergy = accSide_[(std::size_t) b] * inv;
            row.energy = row.midEnergy + row.sideEnergy;
            row.density = row.widthHz > 0.0 ? row.energy / row.widthHz : 0.0;
            const double mom = accMoment_[(std::size_t) b] * inv;
            row.centroidHz = row.energy > 0.0 ? mom / row.energy : 0.0;
            row.centsOffset = row.centroidHz > 0.0 && row.centreHz > 0.0
                            ? 1200.0 * core::det::log2 (row.centroidHz / row.centreHz) : 0.0;
            totalBandEnergy_ += row.energy;
            if (! (std::isfinite (row.energy) && std::isfinite (row.density)
                   && std::isfinite (row.centroidHz) && std::isfinite (row.centsOffset)))
                bad = true;
        }
        // Each frame's contribution was checked finite before it was accumulated, and a float-fed band
        // cannot sum past ~1e83 over any realistic frame count — but "cannot" there is an argument about
        // magnitudes, and the contract is that a non-finite value never reaches a published field. So it
        // is CHECKED: peakShare_ would be inf/inf = NaN if a total overflowed, and that is precisely the
        // zero-that-reads-as-an-answer this instrument refuses to print.
        if (bad || ! std::isfinite (totalBandEnergy_))
        {
            for (auto& row : bands_) { row.centroidHz = 0.0; row.centsOffset = 0.0; }
            noteReason_ = LowEndReason::Overflowed;
            return;
        }
        frameEnergy_ = accFrameEnergy_ * inv;
        bandRangeShare_ = frameEnergy_ > 0.0 ? totalBandEnergy_ / frameEnergy_ : 0.0;
        // `> 0` alone is not "there is a note here": pure DC, pure Nyquist and a denormal-floor signal
        // all leave a positive total made of transform round-off, and naming a note from it is exactly
        // the number-that-reads-as-a-finding this instrument refuses to print.
        if (! (totalBandEnergy_ > 0.0) || ! (bandRangeShare_ > kNoteFloorShare))
        { noteReason_ = LowEndReason::NoEnergy; return; }

        // argmax twice, strictly, so a tie goes to the LOWEST band index
        peakBand_ = 0; peakDensityBand_ = 0;
        for (int b = 1; b < bandCount_; ++b)
        {
            if (bands_[(std::size_t) b].energy > bands_[(std::size_t) peakBand_].energy) peakBand_ = b;
            if (bands_[(std::size_t) b].density > bands_[(std::size_t) peakDensityBand_].density) peakDensityBand_ = b;
        }
        secondBand_ = -1;
        for (int b = 0; b < bandCount_; ++b)
        {
            if (b == peakBand_) continue;
            if (secondBand_ < 0 || bands_[(std::size_t) b].energy > bands_[(std::size_t) secondBand_].energy) secondBand_ = b;
        }
        peakShare_ = bands_[(std::size_t) peakBand_].energy / totalBandEnergy_;

        // the background: the MEDIAN DENSITY of the non-peak bands. Densities, not energies, because a
        // semitone band's width grows with frequency and the energies therefore tilt +5.02 dB across the
        // default range under a perfectly flat spectrum. Even count -> the mean of the two middles.
        int n = 0;
        for (int b = 0; b < bandCount_; ++b) if (b != peakBand_) sort_[(std::size_t) n++] = bands_[(std::size_t) b].density;
        if (n > 0)
        {
            std::sort (sort_.begin(), sort_.begin() + n);
            backgroundDensity_ = (n % 2) == 1 ? sort_[(std::size_t) (n / 2)]
                                              : 0.5 * (sort_[(std::size_t) (n / 2 - 1)] + sort_[(std::size_t) (n / 2)]);
        }
        noteReason_ = LowEndReason::Ok;
    }

    LowEndParams params_ {};
    bool prepared_ = false, finished_ = false;
    double sampleRate_ = 48000.0, binHz_ = 0.0;
    int channels_ = 0, bandCount_ = 0, midiLo_ = 0;
    std::int64_t blockSamples_ = 0;

    SpectrumFrames frames_;
    CrossoverType xover_;   // declared THROUGH the public alias, so the assertion cannot drift from the member   // deterministic coefficients (see core::DetMath)
    core::StateGrid grid_;

    // the clock, and the open block
    std::int64_t totalSamples_ = 0, nextBlockEnd_ = 0, blockIndex_ = 0, blockCount_ = 0;
    double openMid_ = 0.0, openSide_ = 0.0;
    std::int64_t openFinite_ = 0, openHoles_ = 0;

    // 0 lowMid, 1 lowSide, 2 highMid, 3 highSide, 4 rawMid, 5 rawSide — value + Neumaier compensation
    double sum_[kSums] {};
    double comp_[kSums] {};
    std::int64_t finiteSamples_ = 0, holeSamples_ = 0, nonFiniteSamples_ = 0, filterNonFinite_ = 0, absentSamples_ = 0;
    std::int64_t firstHole_ = -1, lastHole_ = -1;

    std::int64_t hist_[kHistogramBins] {};
    std::int64_t skippedBlocks_ = 0;
    std::int64_t histSamples_ = 0;
    double worstFrac_ = -1.0, worstFracEnergy_ = 0.0;
    std::int64_t worstFracBlock_ = -1;
    double peakEnergy_ = -1.0, peakEnergyFrac_ = 0.0;
    std::int64_t peakEnergyBlock_ = -1;
    double peakSide_ = -1.0;
    std::int64_t peakSideBlock_ = -1;
    double peakSideAmp_ = 0.0;
    std::int64_t peakSideAmpAt_ = -1;

    std::int64_t usedFrames_ = 0, holedFrames_ = 0;
    double accFrameEnergy_ = 0.0, frameEnergy_ = 0.0, bandRangeShare_ = 0.0, frameTotal_ = 0.0;
    int peakBand_ = -1, peakDensityBand_ = -1, secondBand_ = -1;
    double backgroundDensity_ = 0.0, peakShare_ = 0.0, totalBandEnergy_ = 0.0;
    LowEndReason widthReason_ = LowEndReason::NotFinished;
    LowEndReason noteReason_  = LowEndReason::NotFinished;

    TraceFn trace_ = nullptr;
    void* traceUser_ = nullptr;

    std::vector<BlockRecord> blocks_;
    std::vector<LowEndBand>  bands_;
    std::vector<BinWeight>   weights_;
    std::vector<int>         bandCountBins_;
    std::vector<double>      accMid_, accSide_, accMoment_, traceMid_, traceSide_, sort_;
    std::vector<std::int64_t> dutyCount_;   // K9: frames this band was present in
    std::vector<double>      dutyRatio_;    // …and the sum of its share of the frame max over those
    double                   dutyShare_ = 0.0;   // the threshold as a linear ratio, made once
    std::int64_t             dutyFrames_ = 0;    // frames that cleared the note floor and were counted
};

} // namespace felitronics::analysis

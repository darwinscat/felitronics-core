// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/analysis/SpectrumFrames.h>
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
//   * DIGITAL SILENCE is the 0/0, and it is the one place a zero would lie, so it does not get one:
//     `widthReason() == NoEnergy`.
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
//     20 ms on) — the first 10 ms alone holds 97.0 % of the file's entire low Side energy. So the
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
//     one band reads A^2/2, and disjoint bands add up to the windowed frame's mean square (Parseval).
//   * RESOLUTION IS THE OBSERVATION LENGTH, and the default is the smallest order that has any. A
//     semitone at 30 Hz is 1.73 Hz wide; a Hann main lobe is 4 bins. At 48 kHz only fftOrder >= 17
//     (0.366 Hz bins, 4.87 bins per band) fits the lobe inside the band, and a tone at a band centre
//     then keeps 99.96 % of its power in its own band (measured against a direct DFT). Orders below
//     that do not resolve the bottom of the range and `underResolvedBands()` counts them rather than
//     hiding it. Zero-padding would not help: it interpolates a peak, it does not separate two tones.
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
    double lowNoteHz   = 30.0;      // the semitone range: every note whose CENTRE lies in [low, high]
    double highNoteHz  = 300.0;     // defaults give MIDI 23..62 (B0..D4), 40 bands
    double tuningHz    = 440.0;     // A4. f(n) = tuningHz * 2^((n-69)/12), MIDI numbering (60 = C4)
    int    fftOrder    = 17;        // N = 1 << fftOrder. 17 is the smallest that resolves a semitone at
                                    // 30 Hz at 48 kHz — see "RESOLUTION" above.
    int    hop         = 0;         // 0 means N/2
    int    maxBlocks   = 1 << 16;   // capacity of the stored 10 ms series: 10.9 min. Everything that is
                                    // not the series itself (integrals, histogram, extrema, counters)
                                    // keeps going past it — law 11, exhaustion is data.
};

// Why a field is not a number. Never NaN: an invalid field has a documented canonical value, and the
// reason is what carries the meaning.
enum class LowEndReason : std::uint8_t
{
    Ok = 0,
    NotFinished,            // finish() has not been called yet
    NoFiniteSamples,        // not one sample of the L/R pair was usable
    NoEnergy,              // every relevant energy is exactly zero — digital silence, the 0/0
    ShorterThanWindow,      // no spectral frame ever closed: the programme is shorter than the window
    NoUsableFrames          // frames closed, but every one of them held a hole
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
                                                     // analysis::LoudnessMeter builds (LoudnessMeter.h:110),
                                                     // so the two instruments name the same intervals.
    static constexpr int    kHistogramBins = 100;    // the side fraction over [0, 1], fixed edges
    static constexpr double kMinSampleRate = 1000.0;
    static constexpr double kMaxSampleRate = 768000.0;
    static constexpr int    kMaxBlocksLimit = 1 << 24;
    static constexpr int    kLobeBins      = 4;      // a Hann main lobe, in bins

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
                                                    + (std::uint64_t) sortDoubles);
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
        if (! (p.crossoverHz > 0.0) || ! (p.crossoverHz <= 0.49 * sampleRate)) return s;   // Svf would clamp it silently
        if (! (p.tuningHz > 0.0) || ! (p.tuningHz < sampleRate)) return s;
        if (! (p.lowNoteHz > 0.0) || ! (p.highNoteHz > p.lowNoteHz)) return s;
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
        prepared_ = false;                                             // law 11b: disarm, validate, write
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

        buildBands();

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
        histSamples_ = 0;
        worstFrac_ = -1.0; worstFracBlock_ = -1; worstFracEnergy_ = 0.0;
        peakEnergy_ = -1.0; peakEnergyBlock_ = -1; peakEnergyFrac_ = 0.0;
        peakSide_ = -1.0; peakSideBlock_ = -1;
        peakSideAmp_ = 0.0; peakSideAmpAt_ = -1;
        usedFrames_ = 0; holedFrames_ = 0;
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
            // A huge finite input can overflow INSIDE the filter. Checked here, on this sample, so the
            // grid's flush below (which also heals poison — Svf.h:161) can never erase an uncounted one.
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
    double highSideFraction() const noexcept { const double t = sum (2) + sum (3); return t > 0.0 ? sum (3) / t : 0.0; }
    double rawSideFraction() const noexcept { const double t = sum (4) + sum (5); return t > 0.0 ? sum (5) / t : 0.0; }

    std::int64_t finiteSamples() const noexcept { return finiteSamples_; }
    std::int64_t holeSamples()   const noexcept { return holeSamples_; }
    std::int64_t nonFiniteSamples() const noexcept { return nonFiniteSamples_; }
    std::int64_t absentSamples() const noexcept { return absentSamples_; }   // the L/R pair was not fully fed
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

    int peakBand() const noexcept { return peakBand_; }                  // argmax of ENERGY — what a tone does
    int peakDensityBand() const noexcept { return peakDensityBand_; }    // argmax of DENSITY — what noise does
    int secondBand() const noexcept { return secondBand_; }              // the runner-up by energy
    double backgroundDensity() const noexcept { return backgroundDensity_; }   // median density of the non-peak bands
    double peakBandEnergy() const noexcept { return peakBand_ >= 0 ? bands_[(std::size_t) peakBand_].energy : 0.0; }
    double peakBandWidthHz() const noexcept { return peakBand_ >= 0 ? bands_[(std::size_t) peakBand_].widthHz : 0.0; }
    double secondBandEnergy() const noexcept { return secondBand_ >= 0 ? bands_[(std::size_t) secondBand_].energy : 0.0; }
    double totalBandEnergy() const noexcept { return totalBandEnergy_; }
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
        return p.tuningHz * std::exp2 ((double) (midi - 69) / 12.0);
    }

    // The band set: every MIDI note whose CENTRE lies inside [lowNoteHz, highNoteHz]. Computed from the
    // parameters alone, once, so it is not a function of the data or of the slicing.
    static bool midiRange (const LowEndParams& p, int& loMidi, int& hiMidi) noexcept
    {
        const double a = 69.0 + 12.0 * std::log2 (p.lowNoteHz  / p.tuningHz);
        const double b = 69.0 + 12.0 * std::log2 (p.highNoteHz / p.tuningHz);
        if (! std::isfinite (a) || ! std::isfinite (b)) return false;
        if (! (a > -2000.0 && b < 2000.0)) return false;
        loMidi = (int) std::ceil (a);
        hiMidi = (int) std::floor (b);
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

    void buildBands() noexcept
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
            hist_[bin] += openFinite_;                                  // duration-weighted: samples, not blocks
            histSamples_ += openFinite_;
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
        ++usedFrames_;
        if (trace_ != nullptr) fireFrameTrace (true);
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
                            ? 1200.0 * std::log2 (row.centroidHz / row.centreHz) : 0.0;
            totalBandEnergy_ += row.energy;
        }
        if (! (totalBandEnergy_ > 0.0)) { noteReason_ = LowEndReason::NoEnergy; return; }

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
    eq::Crossover2 xover_;
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
};

} // namespace felitronics::analysis

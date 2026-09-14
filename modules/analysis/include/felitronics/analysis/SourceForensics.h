// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/analysis/SpectrumFrames.h>
#include <felitronics/core/Config.h>
#include <felitronics/core/Math.h>

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

namespace felitronics::analysis
{

//==============================================================================
// felitronics::analysis::SourceForensics — WHAT THE FILE ACTUALLY WAS, as far as the samples can prove it.
// Two families of evidence, nothing else, and NO attribution:
//
//   1. THE SPECTRAL WALL. A long-term (Welch) mean power spectrum, then the sharp upper edge in it: where
//      it is (in Hz AND as a fraction of this rate's Nyquist), how deep the drop is, how wide the
//      transition is, and separately how far up the spectrum is EMPTY (the upsampling tell).
//   2. THE SAMPLE GRID. Per sample, no FFT: the coarsest dyadic grid 2^-k every finite sample lies on, and
//      how many distinct sample values the stream carries.
//
// WHAT THIS INSTRUMENT REFUSES TO DO. It never says "this is an MP3", "this is upsampled", "this is really
// 16-bit". It publishes coordinates and continuous evidence; the thresholds it does carry are the ones that
// DEFINE the measurement (what counts as a sharp edge, what counts as empty, how wide a plateau is
// sampled), each a named parameter with a documented default. Where the evidence does not reach, the answer
// is `valid = false` plus a reason — never a number that reads as "nothing found".
//
// THE TWO WAYS A NAIVE VERSION OF THIS LIES CONFIDENTLY, and what is done about them:
//
//   · A BAND LIMIT IS A SUFFIX PROPERTY, NOT A LOCAL ONE. The obvious construction — find the steepest
//     local descent, then measure the level just above it — reports a deep NOTCH as a wall: a band-stop, a
//     comb filter or a room null gives a 100 dB "drop" and a one-cell "transition" while full-power
//     spectrum resumes 400 Hz higher. Worse, the notch WINS the search and the real wall above it is then
//     never even evaluated. So the floor of the drop here is `maxAbove` — the loudest cell ANYWHERE above
//     the candidate (a suffix maximum) — and the search maximises that conservative drop. A notch scores
//     ~0 dB because the spectrum comes back; a genuinely band-limited file scores its whole depth. The
//     drop against the LOCAL floor is published beside it as `localDropDb`, and their difference — how far
//     the spectrum comes BACK above the edge — as `recoveryDb`. This also settles the gentle-rolloff case without any special
//     rule: on a monotone 10 dB/kHz slope the loudest cell above a boundary is the cell right at it, so the
//     drop is ~10 dB and the edge is simply not sharp.
//     A suffix MAXIMUM is a blunt instrument, though: one narrow sustained line above the edge — a 19 kHz
//     switching whine, a leaked pilot tone — would be the loudest thing above it and would collapse the
//     drop of a real wall to nothing, and then the argmax moves to an unrelated boundary and publishes ITS
//     frequency. Measured on a -9 dB/octave programme walled at 16 kHz: one 19 kHz line at -45 dB relative
//     to the midband is already louder than the 15 kHz cell, and the reported cutoff jumps from 16 kHz to
//     2.5 kHz. Two guards, and both are needed: the search runs on a 3-cell MEDIAN-FILTERED copy of the
//     cells, and the floor is not the suffix maximum but the (t+1)-th LARGEST cell of the suffix, where
//     t forgives at most `exemptCells` cells and never more than a tenth of the suffix — so the exemption
//     self-disables near Nyquist, where the suffix is short and no forgiveness is wanted. The line is
//     never lost: `sufMaxPower` and `strictDropDb` publish the unforgiven numbers, `exemptedCells` says
//     how many cells were skipped, `recoveryDb` is measured from the STRICT maximum, and the "empty above"
//     test runs on per-BIN maxima and sees every line.
//   · AN EDGE THAT IS NOT THERE MUST NOT BE GIVEN A FREQUENCY. A monotone-falling full-band spectrum has
//     its largest plateau-to-suffix ratio at the LOWEST candidate, because a plateau span fixed in Hz
//     spans 1.6 octaves at 1.3 kHz and 0.19 at 16 kHz; a flat one has it at the highest, where the suffix
//     is shortest. Either way an argmax always returns something. So `minDropDb` is not only a
//     post-hoc flag: a report whose best descent does not reach it is `valid = false`, reason
//     `ShallowerThanMinDrop`, with every number still published as evidence. Measured margin: on white
//     noise the largest drop an argmax can find over 482 candidates is 1.7 dB at one frame and 0.6 dB at
//     sixteen, against a 24 dB default — estimator noise cannot manufacture an edge.
//     Hence three deliberately different resolutions, each published: the edge geometry at cell
//     resolution on median-filtered cells, the emptiness at bin resolution on raw means, and
//     `peakCellPower` — the reference the emptiness threshold hangs from — as the loudest RAW cell. That
//     reference is published, and it needs to be: a louder one raises the threshold and makes emptiness
//     EASIER to claim, and a cell being a 17-bin mean does not stop a dominant tone from setting it.
//   · A POSITION NEAR NYQUIST CANNOT BE ATTRIBUTED. A wall at 0.9-0.95 of Nyquist is where a converter's
//     or an export's anti-alias filter sits, where a 320 kbit/s codec sits, and where a genuinely
//     band-limited master sits. That band is ambiguous BY NATURE, so `nearNyquist` marks it and the
//     instrument stops there. (At 44.1 kHz the default 0.80 is 17.6 kHz and at 48 kHz 19.2 kHz, while
//     128 kbit/s at 16 kHz is 0.73 / 0.67 and stays below it.) `nearNyquist == false` is NOT a claim that
//     something compressed the file: a deliberate 16 kHz low-pass is always possible, and this instrument
//     cannot tell the two apart either. It reports the number.
//
// WHAT THE SAMPLE GRID CAN AND CANNOT PROVE. `gridExponent` k is a fact about the float stream AS IT
// ARRIVED: every finite sample is an exact multiple of 2^-k. From it, `minExactPcmBits = k + 1` is the
// SHORTEST normalised PCM word (values i/2^(b-1), |value| <= 1) that holds every observed sample exactly.
// Read it in that direction only:
//   · It does NOT bound the source's word length from above. A 24-bit file carrying a 16-bit master padded
//     with zeros is indistinguishable from a 16-bit file — that is the whole point of the measurement, and
//     also its limit: the coarse grid is always compatible with storage on a finer one.
//   · The ABSENCE of zero low bits proves nothing at all. 16-bit dithered up into 24, or normalised by a
//     non-dyadic gain, fills the low bits and reports "<= 24 bits" (or no PCM grid at all), which is the
//     honest answer, not "24 bits, really".
//   · A grid finer than 2^-23 is not a PCM word this measurement will name: `pcmCompatible = false`,
//     reason `FinerThan24BitGrid`. A decoded lossy file, anything through a float gain, and 32-bit PCM all
//     land there — 32-bit PCM at EVERY level, measured: float32's top binade [0.5, 1) is spaced 2^-24, one
//     bit finer than a 24-bit grid, so near full scale a 32-bit stream reads k = 24 and lower down k = 31.
//     What DOES read as 24-bit is a 32-bit CONTAINER carrying 24-bit content (every code a multiple of
//     256): k = 23, and a container carrying 16-bit content reads k = 15. Container and content are
//     different questions and only the content is in the samples. (io::Wav decodes 32-bit PCM to double
//     and keeps what a float32 path here would have discarded — see io/Wav.h.)
//   · A sample OUTSIDE [-1, +1) fits no normalised PCM word at any depth — a b-bit word carries
//     i/2^(b-1) for i in [-2^(b-1), 2^(b-1) - 1], so -1 is representable and +1 is not — so
//     `outsidePcmRange` withholds `minExactPcmBits` while `gridExponent` stays valid. Exactly +1.0 is
//     therefore not a PCM word's sample, however routinely a float master's peak sits there.
//
// FORM. setParams / prepare / process / finish / reset, as analysis::ClipDetector. process() is READ-ONLY
// and allocates nothing; `maxBlock` sizes NOTHING (a report must be identical for prepare(sr, 64, nch) and
// prepare(sr, 8192, nch), and a legal call longer than maxBlock is consumed whole). One integer clock, the
// outer loop over samples and the inner over channels. Non-finite samples are holes: counted, never
// measured, never guessed — and they poison exactly the spectral frames that contain them (SpectrumFrames
// does that part). A channel a call does not carry (index >= numChannels, law 11a) is a hole for those
// samples too, counted separately. finish() closes the measurement and freezes the report; process() then
// refuses until reset().
//
// LAW 8a. With the same input bits, parameters, binary and FFT backend, every published field is
// bit-identical under ARBITRARY re-slicing of the stream into process() calls — stronger than law 11(a),
// which only promises agreement with the caller's own boundaries. Bit-exactness across platforms or FFT
// backends is NOT promised. What that costs here: one add per usable frame into a per-channel Neumaier
// per-bin accumulator in frame order; no per-call subtotals and no averaging of per-call averages; power
// averaged LINEARLY with dB taken once at the end; the aggregate spectrum rebuilt from the finished
// per-channel MEANS in ascending channel order (pooling raw sums instead would weight channels by how long
// each was present); quantiles by integer index into a freshly copied scratch, never by permuting the
// accumulator; and the whole cell/edge analysis deferred to finish(), so reading the report cannot change
// it. The capacity of the distinct-value set is published before allocation and its exhaustion is DATA
// (a count plus `distinctComplete`), never a refusal in the middle of a call.
//
// WHAT `transitionHz` CAN MEAN AT BEST. The measured width is bounded below by the analysis window's own
// skirt, so it is a property of the geometry as much as of the file. For a brickwalled passband under a
// periodic Hann window the integrated leakage past d bins falls as 1/(15*pi^2*d^5), which inverts to
//     transitionHz >= binHz * (10^(D/10) / (15*pi^2))^(1/5)        D = the drop being measured, in dB
// (verified against the exact kernel integral beyond d = 10). At fftOrder 14 / 48 kHz that is ~90 Hz for a
// 96 dB drop and ~170 Hz for 110 dB — far under the 1500 Hz default. It is NOT negligible at every legal
// setting: at fftOrder 12 a perfect brickwall over a -140 dB floor measures ~2.7 kHz and would be refused
// `sharp` by its own window. A caller who lowers fftOrder must raise maxTransitionHz with it.
//
// TWO NUMBERS THIS FILE DOES NOT INVENT. The quantile rule is NEAREST-RANK — the order statistic of
// 1-based rank ceil(q*n), clamped to [1, n] — which is the rule dynamics::offline::QuantileHistogram
// already names for this repository (Quantile.h:22), NOT analysis::LoudnessMeter's LRA rule (which takes
// the bin's LOWER bound at rank (N-1)*p + 0.5, faithful to libebur128 and correct there because LRA is a
// DIFFERENCE of two percentiles, so that half-BIN offset cancels), and NOT a third rule of its own. It is
// computed by integer arithmetic, which simply removes the question of whether ceil(q*n) in floating point
// can land on the wrong side of an integer (for these quantiles it does not, checked to n = 20000 — the
// integer form needs no such check). And every dB here is 10*log10 of a power RATIO, inline:
// core::gainToDb is an AMPLITUDE function, 20*log10 with its argument floored at 1e-12, so it would double
// every number and clamp any ratio below 1e-12 to -240 dB.
//
// WHY NO CELL CAN BE NaN, so that a sort of them is not undefined behaviour. Only hole-free frames of
// finite floats are transformed, and for those |X_k| <= N*max|x| <= 2^22 * 3.4e38, so a bin power is under
// 2e90 and the accumulated sum overflows only past ~1e200 frames. The frame producer's own hole gate is
// what keeps that true if a backend ever changes underneath.
//
// RT: an OFFLINE instrument (message thread). prepare() allocates; process/finish/reset do not.
//==============================================================================

enum class ForensicsReason : std::uint8_t
{
    Ok = 0,
    NotFinished,            // read before finish() — the edge analysis runs once, there
    NoChannel,             // asked about a channel outside the prepared width
    ShorterThanWindow,     // the programme is shorter than one FFT window: there is no spectral frame at all
    AllFramesHoled,        // frames existed, but every one of them contained a hole
    NoSpectralEnergy,      // the usable frames carry no power anywhere (digital silence is FINITE, not a hole)
    InsufficientSpan,      // the spectrum is too narrow to hold one plateau span plus a candidate
    NoDownwardEvidence,    // no boundary has more power below it than anywhere above: a flat spectrum has no edge
    ShallowerThanMinDrop,  // the best descent in the search band does not reach minDropDb — there is no edge here
    NoNonZeroSample,       // every finite sample was exactly zero, so no sample witnesses a grid
    FinerThan24BitGrid,    // the grid is finer than 2^-23: no <= 24-bit PCM word holds the stream
    OutsidePcmRange,       // a sample outside [-1, +1) fits no normalised PCM word at any depth
};

// The canonical moments inside process() at which the report's state changes in a way a re-slicing could
// move. They exist so the invariance can be certified from OUTSIDE (see setObserver).
enum class ForensicsEvent : std::uint8_t
{
    FrameClosed = 0,        // a spectral frame closed and was accumulated (or refused); channel == -1
    FirstOffGrid,           // the first sample of `channel` that lies off the <= 24-bit PCM grid
    DistinctSaturated,      // `channel`'s distinct-value set met a new value it had no room for
};

// One measured upper edge of a power spectrum. Every field is published whatever the verdict: the flags are
// the instrument's own definitions, the numbers are the evidence.
struct SpectralWall
{
    bool valid = false;                     // an edge was measured at all; when false, read `reason`
    ForensicsReason reason = ForensicsReason::NotFinished;
    bool sharp = false;                     // the EDGE is sharp by this instrument's definition — never what made it
    bool nearNyquist = false;               // at/above nearNyquistFraction of Nyquist: the position cannot be attributed

    double cutoffHz = 0.0;                  // the band limit: the last cell edge still within transitionStartDb of the plateau
    double cutoffFractionOfNyquist = 0.0;
    double steepestHz = 0.0;                // the boundary that maximised the conservative drop
    double transitionEndHz = 0.0;           // the first cell edge at which the level has reached the local floor
    double transitionHz = 0.0;              // transitionEndHz - cutoffHz; 0 means "narrower than one cell"
    bool   transitionClipped = false;       // a scan ran into its span bound: transitionHz is a LOWER bound
    bool   truncatedAtNyquist = false;      // ...and it ran out of SPECTRUM, not out of span

    // Linear mean power, read from the median-filtered cells the decision itself used — publishing
    // numbers other than the ones that decided would be worse than publishing these.
    double plateauPower = 0.0;              // the median cell of the plateau span below the STEEPEST boundary
    double floorLocalPower = 0.0;           // the median cell of the floor span just above it
    double maxAbovePower = 0.0;             // the loudest cell ANYWHERE above it — the conservative floor
    double sufMaxPower = 0.0;               // the suffix maximum of the RAW cells: nothing above the edge is
                                            // forgiven, neither the rank nor the 3-cell median filter
    int    exemptedCells = 0;               // how many of the loudest cells above the edge the floor skipped
    double dropDb = 0.0;                    // 10log10(plateau / maxAbove); exactly +inf when nothing is above
    double strictDropDb = 0.0;              // ...and the same against sufMaxPower, forgiving nothing
    double localDropDb = 0.0;               // 10log10(plateau / floorLocal): the LOCAL drop, which a notch also has
    // 10log10(sufMax / floorLocal): how far the loudest thing above the edge stands above the local floor
    // just above it, measured from the STRICT maximum so a forgiven line still shows. Read it against its
    // own baseline, not against zero: a random stopband floor's peak-to-median spread over a 2 kHz span is
    // 10-15 dB here by itself (measured), while a sustained line or a genuine recovery stands 30-70 dB up.
    double recoveryDb = 0.0;
    double plateauSpreadDb = 0.0;           // 90th-to-10th-percentile spread of the span below cutoffHz: how flat the plateau was
    double steepnessDbPerOctave = 0.0;      // dropDb per octave of the transition; 0.0 when undefined

    // The second edge, measured the same way over the candidates ABOVE this one — far enough above that its
    // own plateau span is clear of this transition, or it would measure this plateau against that floor and
    // invent an edge. Nested edges are real: a 16 kHz codec wall inside a 20.5 kHz export anti-alias filter.
    // WHICH of the two becomes the primary is decided by the conservative drop and nothing else: with a
    // deep shelf between them the inner edge wins (its plateau is the programme itself) and the outer is
    // published here; with a shallow shelf the OUTER edge wins instead, and then the inner one is not
    // reported at all, because a second search below the primary is structurally useless — everything
    // above such a candidate includes the primary's own plateau, so its conservative drop is ~0 dB.
    bool   secondValid = false;
    bool   secondSharp = false;
    double secondCutoffHz = 0.0;
    double secondDropDb = 0.0;
    double secondTransitionHz = 0.0;

    // "Empty above": tested on the per-BIN suffix maximum, not on cell means — one bin 95 dB down inside a
    // 17-bin cell averages to 107 dB down and would pass a cell-mean test that this one fails.
    // A steep low-pass, a dark master and an upsampled file all produce this — it says only that nothing
    // above emptyAboveHz comes within emptyDb of the loudest cell. Read it beside `sharp` and `dropDb`.
    bool   emptyAboveValid = false;
    double emptyAboveHz = 0.0;              // at BIN resolution (finer than the edge coordinates)
    double emptyAboveFractionOfNyquist = 0.0;
    double emptyThresholdPower = 0.0;       // peakCellPower * 10^(-emptyDb/10), the threshold as applied
    // The reference the emptiness threshold hangs from: the loudest RAW cell. Read it, because the test is
    // relative to it in the direction that matters — a LOUDER reference raises the threshold and makes
    // emptiness EASIER to claim, and a cell is a 17-bin mean, which a dominant tone can still set.
    double peakCellPower = 0.0;

    // the geometry this edge was measured with, published so the resolution is never guessed
    double binHz = 0.0, cellHz = 0.0, nyquistHz = 0.0;
    double searchFromHz = 0.0, searchToHz = 0.0;   // the EFFECTIVE searched band (a plateau span must fit below)
    std::int64_t framesUsed = 0, framesHoled = 0;
};

// What the samples themselves prove about the grid they sit on.
struct SampleGrid
{
    bool valid = false;                     // a grid was witnessed (at least one finite, non-zero sample)
    ForensicsReason reason = ForensicsReason::NotFinished;   // Ok only when minExactPcmBits is readable too
    int  gridExponent = 0;                  // k in 0..149: every finite sample is an exact multiple of 2^-k
    // k <= 23 AND every sample inside [-1, +1): some normalised PCM word of at most 24 bits holds them all.
    // The range is part of the claim, not decoration: a b-bit word carries i/2^(b-1) for i in
    // [-2^(b-1), 2^(b-1) - 1], so -1 IS representable and +1 is not, at any depth.
    bool pcmCompatible = false;
    bool outsidePcmRange = false;           // a sample at or above +1, or below -1
    int  minExactPcmBits = 0;               // k + 1, the SHORTEST exact word; 0 when not readable
    double absPeak = 0.0, sampleMin = 0.0, sampleMax = 0.0;   // the range the claim rests on

    std::int64_t nonZeroSamples = 0, zeroSamples = 0;        // finite; zeros lie on every grid and witness none
    std::int64_t nonFiniteSamples = 0;                       // fed, but NaN or infinite
    std::int64_t absentSamples = 0;                          // the call did not carry this channel (law 11a)
    std::int64_t offGridSamples = 0;                         // finite, non-zero, own k > 23
    std::int64_t firstOffGridSample = -1;                    // its coordinate in samples since reset()
    // gridExponent is a MAX, so one dithered sample, one fade sample or one denormal out of a reverb tail
    // sets it for the whole programme. That is arithmetically right and forensically thin, so the
    // distribution is published beside it (gridExponentHistogram()) and this is where the max was first
    // reached. A 16-bit master with a 24-bit fade is two spikes in that histogram, not one number.
    std::int64_t firstMaxGridSample = -1;

    // EXACT while complete; otherwise a lower bound (the truth is >= this + 1). Exact zero is a value and
    // is counted; -0.0 is the same value as +0.0 and is not counted twice.
    std::int64_t distinctValues = 0;
    bool distinctComplete = false;

    // The same fact read the other way round: how many LOW bits of a `containerBits`-wide PCM word are
    // always zero. A 16-bit master padded into a 24-bit container answers 8 — and answers it about the
    // CONTENT, so it says the container is wider than the content needs, never that the container is 16.
    // 0 when the word length is not readable, or when the container is no wider than the content.
    int alwaysZeroLowBits (int containerBits) const noexcept
    {
        return minExactPcmBits > 0 && containerBits > minExactPcmBits ? containerBits - minExactPcmBits : 0;
    }
};

struct SourceForensicsParams
{
    // --- the spectrum ---
    int fftOrder = 14;                  // N = 1 << fftOrder. 14 at 48 kHz = 0.34 s, 2.93 Hz bins — a codec's
                                        // transition is 200-600 Hz wide, so this resolves it many times over.
    int hop = 0;                        // 0 means N/2. Must satisfy 0 < hop <= N.
    double cellWidthHz = 50.0;          // the spectrum is averaged into cells of about this width, and every
                                        // edge coordinate is a cell edge. The EFFECTIVE width is a whole
                                        // number of bins and is published (cellHz).
    // --- the edge search ---
    double searchFromHz = 1000.0;       // below this an edge is a high-pass, not a band limit. The effective
                                        // start is also pushed up by plateauSpanHz and is published.
    double plateauSpanHz = 2000.0;      // the reference level below a candidate is the MEDIAN cell of this
                                        // span: a median survives a transition occupying up to half of it,
                                        // where a 90th percentile would compare passband peaks with
                                        // stopband minima and inflate the drop.
    double floorSpanHz = 2000.0;        // the LOCAL floor above a candidate, same construction
    double transitionStartDb = 3.0;     // the transition begins where the level leaves the plateau by this much
    double transitionEndDb = 3.0;       // ...and ends where it has come within this much of the local floor
    double minDropDb = 24.0;            // `sharp` needs at least this conservative drop
    double maxTransitionHz = 1500.0;    // ...within at most this width, and not a clipped (lower-bound) one
    int exemptCells = 2;                // the floor above an edge may skip this many of the loudest cells (and
                                        // never more than a tenth of the suffix): a narrow sustained line above
                                        // a real wall must not hide it. 0 restores the strict suffix maximum.
    double nearNyquistFraction = 0.80;  // at/above this fraction of Nyquist the position cannot be attributed
    double emptyDb = 100.0;             // "empty" is this far below the loudest cell
    double emptyMinHz = 200.0;          // ...and must hold over at least this much spectrum (never fewer than
                                        // two bins), so the Nyquist bin alone cannot declare the band empty
    // --- the sample grid ---
    int maxDistinctValues = 1 << 16;    // counted EXACTLY up to here (65536 = every 16-bit code); past it the
                                        // count is a lower bound. The table is sized so this limit sits at or
                                        // below 3/4 load, and that size is published before it is allocated.
};

class SourceForensics
{
public:
    static constexpr double kMinSampleRate = 1000.0;
    static constexpr double kMaxSampleRate = 768000.0;
    static constexpr int    kMinFftOrder   = 8;
    static constexpr int    kMaxFftOrder   = 22;
    static constexpr int    kPcmGridK      = 23;        // a <= 24-bit normalised PCM word means k <= 23
    static constexpr int    kMaxGridK      = 149;       // the smallest positive float is 2^-149
    static constexpr int    kMaxDistinctLimit = 1 << 22;
    static constexpr int    kMaxExemptCells = 16;
    static constexpr double kMaxDb = 400.0;             // a dB parameter past this is a mistake, not a setting
    static constexpr int    kSuffixExemptFraction = 10; // t is also capped at 1/10 of the suffix
    // The plateau reference is the median; the same sort also yields the spread that says how flat it was.
    static constexpr int kMedianNum = 1, kMedianDen = 2;
    static constexpr int kSpreadLoNum = 1, kSpreadLoDen = 10;
    static constexpr int kSpreadHiNum = 9, kSpreadHiDen = 10;

    // An observer of the canonical moments INSIDE process(). It exists because the invariance of this report
    // cannot be certified from the outside: a block-shaped caller only ever looks BETWEEN calls, and the
    // defect that moves a frame onto a call boundary is invisible there — docs/LAW8-KWEIGHTING.md:78 has the
    // measured example, where the integrated result matched bit-for-bit while 5 of 97 intermediate block
    // energies had moved. The type is `noexcept` because process() is: a throwing callback would terminate.
    // The observer must not re-enter this object.
    using Observer = void (*) (void* user, const SourceForensics& self, ForensicsEvent ev, int channel) noexcept;

    // The geometry every sizing and every coordinate is derived from — ONE function, called by both
    // storageFor() and prepare(), so the published budget and the built object cannot disagree.
    struct Geometry
    {
        bool ok = false;
        std::int64_t n = 0, hop = 0;
        int bins = 0;
        double binHz = 0.0, cellHz = 0.0;
        int binsPerCell = 0, cellCount = 0;
        int plateauCells = 0, floorCells = 0;
        int firstCandidate = 0;             // the lowest cell boundary a plateau span fits below
        int emptyMinBins = 0;
        int exemptCells = 0;
        double minDropRatio = 0.0;          // 10^(minDropDb/10): the candidate filter, not just a flag
        int distinctLimit = 0;
        std::size_t tableSlots = 0;         // a power of two, >= 4/3 * distinctLimit
        double startRatio = 0.0, endRatio = 0.0, emptyRatio = 0.0;   // the dB parameters as power ratios
    };

    static Geometry geometryFor (double sampleRate, int maxChannels, const SourceForensicsParams& p) noexcept
    {
        Geometry g;
        if (! (sampleRate >= kMinSampleRate && sampleRate <= kMaxSampleRate)) return g;   // NaN fails too
        if (maxChannels < 1 || maxChannels > core::kMaxChannels) return g;
        if (p.fftOrder < kMinFftOrder || p.fftOrder > kMaxFftOrder) return g;
        if (p.hop < 0) return g;
        // RANGES, not merely finiteness. A finite 1e300 cell width divided by binHz and narrowed to `int`
        // is undefined behaviour, which is the class this project has already paid for twice (a gate on
        // isfinite where a gate on the RANGE was needed; an early return replaced by a clamp). Every
        // frequency parameter is bounded by the sample rate, every dB parameter by kMaxDb.
        if (! (p.cellWidthHz > 0.0 && p.cellWidthHz <= sampleRate)) return g;
        if (! (p.searchFromHz >= 0.0 && p.searchFromHz <= sampleRate)) return g;
        if (! (p.plateauSpanHz > 0.0 && p.plateauSpanHz <= sampleRate)) return g;
        if (! (p.floorSpanHz > 0.0 && p.floorSpanHz <= sampleRate)) return g;
        if (! (p.transitionStartDb >= 0.0 && p.transitionStartDb <= kMaxDb)) return g;
        if (! (p.transitionEndDb >= 0.0 && p.transitionEndDb <= kMaxDb)) return g;
        if (! (p.minDropDb >= 0.0 && p.minDropDb <= kMaxDb)) return g;
        if (! (p.maxTransitionHz > 0.0 && p.maxTransitionHz <= sampleRate)) return g;
        if (! (p.nearNyquistFraction > 0.0 && p.nearNyquistFraction <= 1.0)) return g;
        if (! (p.emptyDb > 0.0 && p.emptyDb <= kMaxDb)) return g;
        if (! (p.emptyMinHz > 0.0 && p.emptyMinHz <= sampleRate)) return g;
        if (p.exemptCells < 0 || p.exemptCells > kMaxExemptCells) return g;
        if (p.maxDistinctValues < 1 || p.maxDistinctValues > kMaxDistinctLimit) return g;

        g.n = (std::int64_t) 1 << p.fftOrder;
        g.hop = p.hop == 0 ? g.n / 2 : (std::int64_t) p.hop;
        if (g.hop < 1 || g.hop > g.n) return g;
        g.bins = (int) (g.n / 2 + 1);
        g.binHz = sampleRate / (double) g.n;
        g.binsPerCell = (int) (p.cellWidthHz / g.binHz);
        if (g.binsPerCell < 1) g.binsPerCell = 1;                       // a cell is never narrower than a bin
        if (g.binsPerCell > g.bins) g.binsPerCell = g.bins;
        g.cellHz = (double) g.binsPerCell * g.binHz;
        // Cells tile the bins, and the TOP cell ABSORBS the remainder instead of standing as a partial cell
        // of its own. Dropping the remainder would hide exactly the top of the spectrum this instrument is
        // looking at; keeping it as its own cell is worse still, and measurably so. At order 12 / 48 kHz the
        // remainder is ONE bin — Nyquist — whose level is structurally below its neighbours', so that cell
        // became the `maxAbove` of every candidate below it and full-band noise reported a 25 dB "wall" at
        // exactly Nyquist with a 0 Hz transition. A cell of 17 bins and a cell of 1 are not comparable
        // measurements, and a suffix maximum compares them.
        g.cellCount = g.bins / g.binsPerCell;
        if (g.cellCount < 1) g.cellCount = 1;
        g.plateauCells = (int) (p.plateauSpanHz / g.cellHz);
        if (g.plateauCells < 1) g.plateauCells = 1;
        g.floorCells = (int) (p.floorSpanHz / g.cellHz);
        if (g.floorCells < 1) g.floorCells = 1;
        const int fromCell = (int) std::ceil (p.searchFromHz / g.cellHz);
        g.firstCandidate = std::max (g.plateauCells, fromCell);         // a FULL plateau span must fit below
        // CEIL, not floor: this is a minimum-width GUARD, and truncating it weakens the very promise it
        // makes (at fftOrder 9 a floor gives 187.5 Hz where 200 was asked for). A SPAN may round down —
        // plateauSpanCells() and floorSpanCells() publish what those became — but a guard may not.
        g.emptyMinBins = (int) std::ceil (p.emptyMinHz / g.binHz);
        // At least TWO bins, whatever emptyMinHz rounds to. With one, the claim can be made about the
        // Nyquist bin alone and reads "empty above 24000 Hz" on a 48 kHz file — vacuously true and
        // published as if it meant something (found at fftOrder 8, where binHz is 187.5 and a 200 Hz
        // minimum floors to a single bin). Two bins also makes the reported coordinate STRICTLY below
        // Nyquist, which is a property the suite asserts directly.
        if (g.emptyMinBins < 2) g.emptyMinBins = 2;
        g.exemptCells = p.exemptCells;
        g.minDropRatio = std::pow (10.0, p.minDropDb / 10.0);
        g.distinctLimit = p.maxDistinctValues;
        g.tableSlots = core::offline::nextPow2 ((std::size_t) g.distinctLimit * 4u / 3u + 1u);
        g.startRatio = std::pow (10.0, -p.transitionStartDb / 10.0);
        g.endRatio   = std::pow (10.0,  p.transitionEndDb   / 10.0);
        g.emptyRatio = std::pow (10.0, -p.emptyDb / 10.0);
        g.ok = true;
        return g;
    }

    // WHAT prepare() ASKS THE HEAP FOR (law 11d), from the function prepare() sizes itself with — the frame
    // producer's own budget included, since this object owns one.
    struct Storage
    {
        bool ok = false;
        SpectrumFrames::Storage frames {};
        std::size_t powerSum   = 0;      // double, bins * channels
        std::size_t powerComp  = 0;      // double, bins * channels — the Neumaier compensation
        std::size_t meanBins   = 0;      // double, bins — one channel's mean, or the aggregate, at a time
        std::size_t cells      = 0;      // double, cellCount — the raw cell means
        std::size_t smoothCells = 0;     // double, cellCount — their 3-cell median, what the edge search reads
        std::size_t suffixMax  = 0;      // double, cellCount — the forgiving floor the search uses
        std::size_t strictMax  = 0;      // double, cellCount — and the raw, unforgiving one, published beside it
        std::size_t sortScratch = 0;     // double, max(plateauCells, floorCells)
        std::size_t tableSlots = 0;      // uint32, slots * channels
        std::uint64_t bytes() const noexcept
        {
            return frames.bytes()
                 + (std::uint64_t) sizeof (double) * (powerSum + powerComp + meanBins + cells + smoothCells
                                                      + suffixMax + strictMax + sortScratch)
                 + (std::uint64_t) sizeof (std::uint32_t) * tableSlots;
        }
    };

    static Storage storageFor (double sampleRate, int maxChannels, const SourceForensicsParams& p) noexcept
    {
        Storage s;
        const Geometry g = geometryFor (sampleRate, maxChannels, p);
        if (! g.ok) return s;
        SpectrumFramesParams fp;
        fp.fftOrder = p.fftOrder;
        fp.hop = p.hop;
        s.frames = SpectrumFrames::storageFor (sampleRate, maxChannels, fp);
        if (! s.frames.ok) return s;
        const std::size_t ch = (std::size_t) maxChannels;
        s.powerSum    = (std::size_t) g.bins * ch;
        s.powerComp   = (std::size_t) g.bins * ch;
        s.meanBins    = (std::size_t) g.bins;
        s.cells       = (std::size_t) g.cellCount;
        s.smoothCells = (std::size_t) g.cellCount;
        s.suffixMax   = (std::size_t) g.cellCount;
        s.strictMax   = (std::size_t) g.cellCount;
        s.sortScratch = (std::size_t) std::max (g.plateauCells, g.floorCells);
        s.tableSlots  = g.tableSlots * ch;
        s.ok = true;
        return s;
    }

    void setParams (const SourceForensicsParams& p) noexcept { params_ = p; }   // takes effect at the next prepare()
    const SourceForensicsParams& params() const noexcept { return params_; }
    // What the CURRENT measurement is being made with — the snapshot prepare() took. setParams() between
    // prepare() and finish() must not move a published number, so the edge analysis reads this and never
    // `params_` (found in the code-review round: minDropDb was read at finish(), so raising it afterwards
    // turned a measured 79 dB edge from valid to invalid without a new prepare()).
    const SourceForensicsParams& activeParams() const noexcept { return active_; }

    void setObserver (Observer fn, void* user) noexcept { obs_ = fn; obsUser_ = user; }

    [[nodiscard]] bool prepare (double sampleRate, int /*maxBlock: nothing is sized by it*/, int maxChannels) noexcept
    {
        // Law 11b: disarm, validate, write — and DISARM MEANS THE REPORT TOO. A failed prepare() after a
        // finished measurement used to leave isPrepared() false while wall(0).valid stayed true, because
        // the channel count and the finished flag survived (code-review round).
        prepared_ = false;
        finished_ = false;
        channels_ = 0;
        const Geometry g = geometryFor (sampleRate, maxChannels, params_);
        const Storage st = storageFor (sampleRate, maxChannels, params_);
        if (! g.ok || ! st.ok) return false;
        SpectrumFramesParams fp;
        fp.fftOrder = params_.fftOrder;
        fp.hop = params_.hop;
        frames_.setParams (fp);
        if (! frames_.prepare (sampleRate, 0, maxChannels)) return false;
        geo_ = g;
        active_ = params_;                                               // the snapshot the report is made with
        sampleRate_ = sampleRate;
        channels_ = maxChannels;
        bins_ = g.bins;
        tableMask_ = g.tableSlots - 1u;
        sum_.assign (st.powerSum, 0.0);
        comp_.assign (st.powerComp, 0.0);
        mean_.assign (st.meanBins, 0.0);
        cells_.assign (st.cells, 0.0);
        smooth_.assign (st.smoothCells, 0.0);
        suffix_.assign (st.suffixMax, 0.0);
        strict_.assign (st.strictMax, 0.0);
        sortScratch_.assign (st.sortScratch, 0.0);
        table_.assign (st.tableSlots, kEmptySlot);
        prepared_ = true;
        reset();
        return true;
    }

    void reset() noexcept
    {
        t_ = 0;
        finished_ = false;
        aggFramesUsed_ = 0;
        aggFramesHoled_ = 0;
        for (auto& v : sum_)  v = 0.0;
        for (auto& v : comp_) v = 0.0;
        for (auto& v : table_) v = kEmptySlot;
        for (int c = 0; c < core::kMaxChannels; ++c)
        {
            chans_[c] = Channel {};
            walls_[c] = SpectralWall {};
        }
        walls_[core::kMaxChannels] = SpectralWall {};
        frames_.reset();
    }

    static constexpr int latencySamples() noexcept { return 0; }          // a read-only sink

    // READ-ONLY. Law 11's order: malformed → unprepared → finished → nch > maxChannels → n == 0 → run.
    // A legal call of ANY length is consumed whole (maxBlock sizes nothing); a refused one consumes nothing.
    [[nodiscard]] bool process (const float* const* in, int numChannels, int n) noexcept
    {
        if (numChannels < 0 || n < 0) return false;
        if (! prepared_ || finished_) return false;
        if (numChannels > channels_) return false;
        if (n == 0) return true;
        for (int i = 0; i < n; ++i)
        {
            // The sample-domain family and the frame producer advance TOGETHER, sample by sample. Doing all
            // of a call's grid work first and its frames afterwards would show a frame callback at sample t
            // the counters of samples that come after t, and re-slicing at t would then change the trace.
            for (int c = 0; c < channels_; ++c)
            {
                const bool fed = c < numChannels;
                const float x = fed ? in[c][i] : 0.0f;
                observeSample (c, x, fed);
                frames_.push (c, x, fed);
            }
            const bool closed = frames_.tick();
            ++t_;                                                        // one clock, advanced exactly once
            if (closed) closeFrame();
        }
        return true;
    }

    // End of stream. Runs the edge analysis ONCE (it is a function of the accumulated spectrum, not of the
    // stream's boundaries), publishes the uncovered tail, freezes. Idempotent; allocates nothing.
    void finish() noexcept
    {
        if (! prepared_ || finished_) return;
        frames_.finish();
        for (int c = 0; c < channels_; ++c)
        {
            fillChannelMean (c);
            walls_[c] = measureWall (chans_[c].framesUsed, chans_[c].framesHoled);
        }
        (void) fillAggregateMean();                                  // aggFramesUsed_ == 0 -> the no-frames report
        walls_[channels_] = measureWall (aggFramesUsed_, aggFramesHoled_);
        finished_ = true;
    }

    // --- the report (the edge analysis after finish(); the counters live while streaming) ---
    bool isPrepared() const noexcept { return prepared_; }
    bool isFinished() const noexcept { return finished_; }
    int  channels() const noexcept { return channels_; }
    std::int64_t samplesProcessed() const noexcept { return t_; }
    std::int64_t tailUncoveredSamples() const noexcept { return prepared_ ? frames_.tailUncoveredSamples() : 0; }

    // Total accessors: a channel outside the prepared width answers an empty report with reason NoChannel
    // rather than reading past the array (the P58 crew round's rule for every published accessor).
    SpectralWall wall (int c) const noexcept
    {
        if (c < 0 || c >= channels_) { SpectralWall w; w.reason = ForensicsReason::NoChannel; return w; }
        return walls_[c];
    }
    // The file's edge: the per-channel MEANS averaged in ascending channel order (not the pooled sums,
    // which would weight a channel by how many frames it was present for).
    SpectralWall wall() const noexcept
    {
        if (! prepared_) { SpectralWall w; w.reason = ForensicsReason::NoChannel; return w; }
        return walls_[channels_];
    }
    // Computed from the counters on every read, so it is LIVE while streaming and final afterwards (the
    // counters cannot move once finish() has refused further calls). That is also what lets a frame
    // observer see the sample-domain family's state at the exact sample a frame closed.
    SampleGrid sampleGrid (int c) const noexcept
    {
        if (c < 0 || c >= channels_) { SampleGrid g; g.reason = ForensicsReason::NoChannel; return g; }
        return measureGrid (c);
    }

    // --- the raw evidence, published for a trace and for a cross-toolchain diff ---
    // The accumulated per-bin power SUM (the Neumaier main term). Linear, pre-division, pre-dB.
    const double* powerSum (int c) const noexcept
    {
        if (c < 0 || c >= channels_ || sum_.empty()) return nullptr;
        return sum_.data() + (std::size_t) c * (std::size_t) bins_;
    }
    const double* powerSumCompensation (int c) const noexcept
    {
        if (c < 0 || c >= channels_ || comp_.empty()) return nullptr;
        return comp_.data() + (std::size_t) c * (std::size_t) bins_;
    }
    // The Welch mean of bin `bin` — the sum and its compensation over the channel's usable frame count.
    double meanPower (int c, int bin) const noexcept
    {
        if (c < 0 || c >= channels_ || bin < 0 || bin >= bins_ || chans_[c].framesUsed <= 0) return 0.0;
        const std::size_t i = (std::size_t) c * (std::size_t) bins_ + (std::size_t) bin;
        return (sum_[i] + comp_[i]) / (double) chans_[c].framesUsed;
    }
    // kPcmGridK + 3 counters: k = 0..23 exactly, then 24, then one bucket for every finer grid (25..149).
    // A max is one sample; this is the distribution behind it.
    const std::int64_t* gridExponentHistogram (int c) const noexcept
    {
        return c >= 0 && c < channels_ ? chans_[c].kHist : nullptr;
    }
    static constexpr int gridExponentBuckets() noexcept { return kPcmGridK + 3; }
    std::int64_t framesUsed (int c) const noexcept { return c >= 0 && c < channels_ ? chans_[c].framesUsed : 0; }
    std::int64_t framesHoled (int c) const noexcept { return c >= 0 && c < channels_ ? chans_[c].framesHoled : 0; }
    bool lastFrameUsed (int c) const noexcept { return c >= 0 && c < channels_ && chans_[c].lastUsed; }
    // The frame producer, for a trace: the frame just closed (index, start, finite flags, raw power bins).
    const SpectrumFrames& frames() const noexcept { return frames_; }

    // --- the geometry, published so nothing about the resolution is guessed ---
    int bins() const noexcept { return bins_; }
    int cellCount() const noexcept { return geo_.cellCount; }
    int binsPerCell() const noexcept { return geo_.binsPerCell; }
    double binHz() const noexcept { return geo_.binHz; }
    double cellHz() const noexcept { return geo_.cellHz; }
    double sampleRate() const noexcept { return sampleRate_; }
    std::int64_t windowSamples() const noexcept { return geo_.n; }
    std::int64_t hopSamples() const noexcept { return geo_.hop; }
    // The lowest boundary the search CONSIDERS (a full plateau span must fit below it, which can push it
    // above the requested searchFromHz). A reported cutoffHz may sit up to one plateau span BELOW this: the
    // band limit is found by scanning down from the winning candidate.
    double searchFromHz() const noexcept { return (double) geo_.firstCandidate * geo_.cellHz; }
    double searchToHz() const noexcept { return (double) std::max (0, geo_.cellCount - 1) * geo_.cellHz; }
    int plateauSpanCells() const noexcept { return geo_.plateauCells; }   // the EFFECTIVE span, in cells
    int floorSpanCells() const noexcept { return geo_.floorCells; }
    int exemptCells() const noexcept { return geo_.exemptCells; }
    int distinctLimit() const noexcept { return geo_.distinctLimit; }
    std::size_t distinctTableSlots() const noexcept { return geo_.tableSlots; }

private:
    static constexpr std::uint32_t kEmptySlot = 0xFFFFFFFFu;   // a quiet-NaN pattern: never a finite sample

    struct Channel
    {
        std::int64_t framesUsed = 0, framesHoled = 0;
        bool lastUsed = false;
        int  maxK = -1;                                        // no grid witnessed yet
        std::int64_t firstMaxK = -1;
        std::int64_t kHist[kPcmGridK + 3] {};                  // k = 0..23, then 24, then "25 or finer"
        double absPeak = 0.0, minVal = 0.0, maxVal = 0.0;
        bool rangeSeen = false;
        std::int64_t nonZero = 0, zero = 0, nonFinite = 0, absent = 0, offGrid = 0;
        std::int64_t firstOffGrid = -1;
        std::int64_t distinct = 0;
        bool distinctComplete = true;
    };

    // One measured edge, before it is written into the report.
    struct Edge
    {
        bool valid = false;
        ForensicsReason reason = ForensicsReason::NoDownwardEvidence;
        int jc = 0, js = 0, je = 0;
        double plateau = 0.0, floorLocal = 0.0, maxAbove = 0.0, sufMax = 0.0, spreadLo = 0.0, spreadHi = 0.0;
        int exempted = 0;
        bool clipped = false, truncated = false;
    };

    // The coarsest dyadic grid a finite, non-zero float lies on: the smallest k >= 0 with x*2^k an integer.
    // Built from ClipDetector's private helper (ClipDetector.h:327) and differs in one documented way: that
    // one answers a sentinel just past its own cap for every subnormal, which is all a capped comparison
    // needs, while this one PUBLISHES k as a number, so a sentinel would be a lie. The exact range is
    // 0..149 — the smallest positive subnormal is 2^-149. (A shared primitive is the right home for this;
    // ClipDetector is read-only in this branch, so the duplication is deliberate and noted.)
    static int gridExponentOf (float x) noexcept
    {
        const std::uint32_t bits = std::bit_cast<std::uint32_t> (x);
        const int ex = (int) ((bits >> 23) & 0xFFu);
        const std::uint32_t frac = bits & 0x7FFFFFu;
        if (ex == 0) return kMaxGridK - std::countr_zero (frac);         // subnormal: frac * 2^-149, frac != 0
        return std::max (0, 150 - ex - std::countr_zero (frac | 0x800000u));
    }

    // -0.0 and +0.0 are the same VALUE and must not count twice.
    static std::uint32_t valueKey (float x) noexcept
    {
        const std::uint32_t bits = std::bit_cast<std::uint32_t> (x);
        return bits == 0x80000000u ? 0u : bits;
    }

    // The median of three, by min/max alone: exact, no libm, no sort, no tie to resolve.
    static double med3 (double a, double b, double c) noexcept
    {
        return std::max (std::min (a, b), std::min (std::max (a, b), c));
    }

    static std::uint32_t mix (std::uint32_t x) noexcept                  // a 32-bit finaliser: no libm, exact
    {
        x ^= x >> 16; x *= 0x7feb352du;
        x ^= x >> 15; x *= 0x846ca68bu;
        x ^= x >> 16;
        return x;
    }

    void notify (ForensicsEvent ev, int c) const noexcept { if (obs_ != nullptr) obs_ (obsUser_, *this, ev, c); }

    // One sample of one channel, at coordinate t_. Non-finite values are gated BEFORE any arithmetic.
    void observeSample (int c, float x, bool fed) noexcept
    {
        Channel& ch = chans_[c];
        if (! fed) { ++ch.absent; return; }
        if (! std::isfinite (x)) { ++ch.nonFinite; return; }
        const double a = std::fabs ((double) x);
        if (a > ch.absPeak) ch.absPeak = a;
        const double v = (double) x;
        if (! ch.rangeSeen) { ch.minVal = v; ch.maxVal = v; ch.rangeSeen = true; }
        else { if (v < ch.minVal) ch.minVal = v; if (v > ch.maxVal) ch.maxVal = v; }
        insertValue (c, valueKey (x));
        if (core::exactlyEqual (x, 0.0f)) { ++ch.zero; return; }         // on every grid, witness to none
        ++ch.nonZero;
        const int k = gridExponentOf (x);
        ++ch.kHist[k <= kPcmGridK + 1 ? k : kPcmGridK + 2];
        if (k > ch.maxK) { ch.maxK = k; ch.firstMaxK = t_; }
        if (k > kPcmGridK)
        {
            ++ch.offGrid;
            if (ch.firstOffGrid < 0) { ch.firstOffGrid = t_; notify (ForensicsEvent::FirstOffGrid, c); }
        }
    }

    // Open addressing, linear probing, fixed capacity. At the limit LOOKUPS CONTINUE, so a stream that
    // repeats exactly `distinctLimit` values forever stays exactly countable; the set is marked incomplete
    // only by a genuinely NEW value it has no room for. The probe is bounded by the table, and the limit
    // sits at or below 3/4 load, so an empty slot always exists.
    void insertValue (int c, std::uint32_t key) noexcept
    {
        Channel& ch = chans_[c];
        // Once the set is known incomplete, both published fields are frozen — the count at the limit and
        // the flag — so every further probe is a proven no-op. Skipping them keeps a float stream O(1)
        // instead of paying an unsuccessful linear probe per sample for the rest of the programme.
        if (! ch.distinctComplete) return;
        std::uint32_t* tab = table_.data() + (std::size_t) c * geo_.tableSlots;
        std::size_t i = (std::size_t) mix (key) & tableMask_;
        for (std::size_t probe = 0; probe < geo_.tableSlots; ++probe)
        {
            const std::uint32_t slot = tab[i];
            if (slot == kEmptySlot)
            {
                if (ch.distinct >= (std::int64_t) geo_.distinctLimit)
                {
                    if (ch.distinctComplete)
                    {
                        ch.distinctComplete = false;
                        notify (ForensicsEvent::DistinctSaturated, c);
                    }
                    return;
                }
                tab[i] = key;
                ++ch.distinct;
                return;
            }
            if (slot == key) return;
            i = (i + 1u) & tableMask_;
        }
    }

    // A frame closed on this very sample: accumulate it, or refuse it whole. One add per usable frame, in
    // frame order, ascending bins, into a per-channel accumulator — the order that makes the sum's bits a
    // function of the samples alone.
    void closeFrame() noexcept
    {
        for (int c = 0; c < channels_; ++c)
        {
            Channel& ch = chans_[c];
            ch.lastUsed = frames_.frameFinite (c);
            if (! ch.lastUsed) { ++ch.framesHoled; continue; }
            ++ch.framesUsed;
            const double* p = frames_.power (c);
            double* s = sum_.data() + (std::size_t) c * (std::size_t) bins_;
            double* k = comp_.data() + (std::size_t) c * (std::size_t) bins_;
            for (int b = 0; b < bins_; ++b)
            {
                const double v = p[(std::size_t) b];
                const double sb = s[(std::size_t) b];
                const double tt = sb + v;                                 // Neumaier
                k[(std::size_t) b] += std::fabs (sb) >= std::fabs (v) ? (sb - tt) + v : (v - tt) + sb;
                s[(std::size_t) b] = tt;
            }
        }
        notify (ForensicsEvent::FrameClosed, -1);
    }

    // --- finish(): the sample-grid report ---
    SampleGrid measureGrid (int c) const noexcept
    {
        const Channel& ch = chans_[c];
        SampleGrid g;
        g.nonZeroSamples = ch.nonZero;
        g.zeroSamples = ch.zero;
        g.nonFiniteSamples = ch.nonFinite;
        g.absentSamples = ch.absent;
        g.offGridSamples = ch.offGrid;
        g.firstOffGridSample = ch.firstOffGrid;
        g.firstMaxGridSample = ch.firstMaxK;
        g.distinctValues = ch.distinct;
        g.distinctComplete = ch.distinctComplete;
        g.absPeak = ch.absPeak;
        g.sampleMin = ch.minVal;
        g.sampleMax = ch.maxVal;
        g.outsidePcmRange = ch.rangeSeen && ! (ch.minVal >= -1.0 && ch.maxVal < 1.0);
        if (ch.maxK < 0) { g.reason = ForensicsReason::NoNonZeroSample; return g; }
        g.valid = true;
        g.gridExponent = ch.maxK;
        g.pcmCompatible = ch.maxK <= kPcmGridK && ! g.outsidePcmRange;
        if (ch.maxK > kPcmGridK) { g.reason = ForensicsReason::FinerThan24BitGrid; return g; }
        if (g.outsidePcmRange)  { g.reason = ForensicsReason::OutsidePcmRange; return g; }
        g.minExactPcmBits = ch.maxK + 1;
        g.reason = ForensicsReason::Ok;
        return g;
    }

    // --- finish(): the spectrum ---
    void fillChannelMean (int c) noexcept
    {
        const std::int64_t f = chans_[c].framesUsed;
        const std::size_t base = (std::size_t) c * (std::size_t) bins_;
        if (f <= 0) { for (int b = 0; b < bins_; ++b) mean_[(std::size_t) b] = 0.0; return; }
        // DIVISION, not a multiplication by a reciprocal: meanPower() below publishes (sum + comp) / frames,
        // and a report derived from x * (1/f) would be a report about numbers the accessor never showed.
        for (int b = 0; b < bins_; ++b)
            mean_[(std::size_t) b] = (sum_[base + (std::size_t) b] + comp_[base + (std::size_t) b]) / (double) f;
    }

    // The aggregate: the mean of the per-channel MEANS, in ascending channel order. Returns how many
    // channels took part. Pooling the raw sums instead would weight a channel by its usable duration — a
    // channel present for one frame and one present for a hundred would not contribute equally.
    int fillAggregateMean() noexcept
    {
        for (int b = 0; b < bins_; ++b) mean_[(std::size_t) b] = 0.0;
        int usable = 0;
        aggFramesUsed_ = 0;
        aggFramesHoled_ = 0;
        for (int c = 0; c < channels_; ++c)
        {
            const std::int64_t f = chans_[c].framesUsed;
            aggFramesHoled_ += chans_[c].framesHoled;
            if (f <= 0) continue;
            ++usable;
            aggFramesUsed_ += f;
            const std::size_t base = (std::size_t) c * (std::size_t) bins_;
            for (int b = 0; b < bins_; ++b)
                mean_[(std::size_t) b] += (sum_[base + (std::size_t) b] + comp_[base + (std::size_t) b]) / (double) f;
        }
        if (usable > 0)
            for (int b = 0; b < bins_; ++b) mean_[(std::size_t) b] /= (double) usable;
        return usable;
    }

    void fillGeometry (SpectralWall& w) const noexcept
    {
        w.binHz = geo_.binHz;
        w.cellHz = geo_.cellHz;
        w.nyquistHz = sampleRate_ * 0.5;
        w.searchFromHz = searchFromHz();
        w.searchToHz = searchToHz();
    }

    // Everything below reads mean_ (one channel's mean, or the aggregate) and writes one SpectralWall.
    SpectralWall measureWall (std::int64_t framesUsedIn, std::int64_t framesHoledIn) noexcept
    {
        SpectralWall w;
        fillGeometry (w);
        w.framesUsed = framesUsedIn;
        w.framesHoled = framesHoledIn;
        if (framesUsedIn <= 0)
        {
            w.reason = frames_.frameCount() == 0 ? ForensicsReason::ShorterThanWindow
                                                 : ForensicsReason::AllFramesHoled;
            return w;
        }
        // cells: the mean over each cell's own bins, so the top partial cell is a mean too
        double peakCell = 0.0;
        for (int j = 0; j < geo_.cellCount; ++j)
        {
            const int from = j * geo_.binsPerCell;
            const int to = j == geo_.cellCount - 1 ? bins_ : from + geo_.binsPerCell;   // the top cell takes the rest
            double acc = 0.0;
            for (int b = from; b < to; ++b) acc += mean_[(std::size_t) b];
            const double v = acc / (double) (to - from);
            cells_[(std::size_t) j] = v;
            if (v > peakCell) peakCell = v;
        }
        w.peakCellPower = peakCell;
        // A 3-cell median of the cells, and the suffix maximum OF THAT: an isolated narrow line above the
        // edge is removed, a real recovery (many cells wide) is not. The first and last cells keep their
        // raw value — a median of three needs three.
        for (int j = 0; j < geo_.cellCount; ++j)
            smooth_[(std::size_t) j] = (j == 0 || j == geo_.cellCount - 1)
                                     ? cells_[(std::size_t) j]
                                     : med3 (cells_[(std::size_t) (j - 1)], cells_[(std::size_t) j],
                                             cells_[(std::size_t) (j + 1)]);
        // The floor above each boundary, in one backward pass: `strict_` is the suffix maximum of the RAW
        // cells — forgiving nothing, not even the median filter, since a 3-cell median is itself
        // forgiveness and an "unforgiven" number computed from filtered cells would be neither (found in
        // the code-review round: it read 27 dB of recovery where the raw suffix says 63) — and `suffix_`
        // is the (t+1)-th largest of the FILTERED cells, t forgiving at most exemptCells of them and at
        // most a tenth of the suffix. The top-(exemptCells+1) values are held in a fixed array of at most
        // 17 doubles, kept descending by one insertion step — no allocation, no sort, no tie to resolve.
        {
            const int cap = geo_.exemptCells + 1;
            double top[kMaxExemptCells + 1] {};
            int held = 0;
            double rawRun = 0.0;
            for (int j = geo_.cellCount - 1; j >= 0; --j)
            {
                if (cells_[(std::size_t) j] > rawRun) rawRun = cells_[(std::size_t) j];
                strict_[(std::size_t) j] = rawRun;
                const double v = smooth_[(std::size_t) j];
                if (held < cap) { top[held] = v; ++held; }
                else if (v > top[cap - 1]) top[cap - 1] = v;
                for (int i = held - 1; i > 0 && top[i] > top[i - 1]; --i) std::swap (top[i], top[i - 1]);
                const int above = geo_.cellCount - j;
                int t = (above - 1) / kSuffixExemptFraction;
                if (t > geo_.exemptCells) t = geo_.exemptCells;
                if (t > held - 1) t = held - 1;
                suffix_[(std::size_t) j] = top[t];
            }
        }
        // "empty above", on the per-BIN suffix maximum
        measureEmptyAbove (w, peakCell);

        if (geo_.firstCandidate > geo_.cellCount - 1) { w.reason = ForensicsReason::InsufficientSpan; return w; }
        if (! (peakCell > 0.0)) { w.reason = ForensicsReason::NoSpectralEnergy; return w; }

        Edge e;
        if (! findEdge (geo_.firstCandidate, geo_.cellCount, e)) { w.reason = e.reason; return w; }
        writeEdge (w, e);
        // An argmax always returns something, so the instrument's own minimum drop is the gate on whether
        // there is an edge AT ALL — otherwise a monotone-falling programme is handed the lowest boundary in
        // the band and a flat one the highest, each with a frequency that means nothing. The numbers stay
        // published; only `valid` is withheld, with the reason that says why.
        if (! (w.dropDb >= active_.minDropDb)) { w.reason = ForensicsReason::ShallowerThanMinDrop; return w; }
        // the second edge: the same measurement, above this one's transition and clear of it by a plateau span
        Edge e2;
        const int secondFrom = std::max (e.je + geo_.plateauCells, geo_.firstCandidate);
        if (secondFrom < geo_.cellCount && findEdge (secondFrom, geo_.cellCount, e2)
            && e2.plateau >= e2.maxAbove * geo_.minDropRatio)
        {
            w.secondValid = true;
            w.secondCutoffHz = (double) e2.js * geo_.cellHz;
            w.secondDropDb = powerRatioDb (e2.plateau, e2.maxAbove);
            w.secondTransitionHz = (double) (e2.je - e2.js) * geo_.cellHz;
            w.secondSharp = ! e2.clipped && w.secondDropDb >= active_.minDropDb
                          && w.secondTransitionHz <= active_.maxTransitionHz;
        }
        w.valid = true;
        w.reason = ForensicsReason::Ok;
        return w;
    }

    void measureEmptyAbove (SpectralWall& w, double peakCell) const noexcept
    {
        if (! (peakCell > 0.0)) return;
        const double thr = peakCell * geo_.emptyRatio;
        w.emptyThresholdPower = thr;
        double run = 0.0;
        int edge = -1;
        for (int b = bins_ - 1; b >= 0; --b)
        {
            if (mean_[(std::size_t) b] > run) run = mean_[(std::size_t) b];
            if (run <= thr && bins_ - b >= geo_.emptyMinBins) edge = b;   // the LOWEST bin that still holds
        }
        if (edge < 0) return;
        w.emptyAboveValid = true;
        w.emptyAboveHz = (double) edge * geo_.binHz;
        w.emptyAboveFractionOfNyquist = w.nyquistHz > 0.0 ? w.emptyAboveHz / w.nyquistHz : 0.0;
    }

    // The edge search over the candidate boundaries [jFrom, jTo). The metric is the CONSERVATIVE drop —
    // the plateau below against the loudest cell ANYWHERE above — compared by cross-multiplication, so
    // there is no division and an empty suffix is simply the largest possible drop. A candidate must have
    // energy below it and must actually descend, or a flat spectrum would be handed a cutoff.
    bool findEdge (int jFrom, int jTo, Edge& e) noexcept
    {
        int best = -1;
        double bestP = 0.0, bestM = 0.0;
        for (int j = jFrom; j < jTo; ++j)
        {
            const double p = pick (loadSpan (j - geo_.plateauCells, j), kMedianNum, kMedianDen);
            const double m = suffix_[(std::size_t) j];
            if (! (p > 0.0) || ! (p > m)) continue;
            if (best < 0 || p * bestM > bestP * m) { best = j; bestP = p; bestM = m; }   // ties keep the lower
        }
        if (best < 0) { e.reason = ForensicsReason::NoDownwardEvidence; return false; }

        e.jc = best;
        e.maxAbove = bestM;
        e.sufMax = strict_[(std::size_t) best];
        const int above = geo_.cellCount - best;
        e.exempted = std::min (geo_.exemptCells, (above - 1) / kSuffixExemptFraction);
        // the plateau reference and its own spread, from ONE sort of the span below the winner
        e.plateau = pick (loadSpan (best - geo_.plateauCells, best), kMedianNum, kMedianDen);
        const int floorTo = std::min (best + geo_.floorCells, geo_.cellCount);
        e.floorLocal = pick (loadSpan (best, floorTo), kMedianNum, kMedianDen);

        // DOWN from the winner to the last cell still within transitionStartDb of the plateau
        const int lowBound = best - geo_.plateauCells;
        int js = -1;
        for (int j = best; j > lowBound; --j)
            if (smooth_[(std::size_t) (j - 1)] >= e.plateau * geo_.startRatio) { js = j; break; }
        if (js < 0) { e.clipped = true; js = lowBound; }
        // UP from there to the first cell that has reached the local floor
        int je = -1;
        for (int j = js; j < floorTo; ++j)
            if (smooth_[(std::size_t) j] <= e.floorLocal * geo_.endRatio) { je = j; break; }
        if (je < 0)
        {
            e.clipped = true;
            je = floorTo;
            e.truncated = floorTo >= geo_.cellCount;
        }
        e.js = js;
        e.je = je;
        // The plateau's own flatness, measured over the span below the BAND LIMIT rather than below the
        // steepest boundary: the winner usually sits just above the edge, so its own span straddles the
        // transition and its spread would say "the span was mixed" instead of "the plateau was ragged".
        const std::size_t ns = loadSpan (js - geo_.plateauCells, js);
        e.spreadLo = pick (ns, kSpreadLoNum, kSpreadLoDen);
        e.spreadHi = pick (ns, kSpreadHiNum, kSpreadHiDen);
        e.valid = true;
        e.reason = ForensicsReason::Ok;
        return true;
    }

    void writeEdge (SpectralWall& w, const Edge& e) const noexcept
    {
        w.cutoffHz = (double) e.js * geo_.cellHz;
        w.cutoffFractionOfNyquist = w.nyquistHz > 0.0 ? w.cutoffHz / w.nyquistHz : 0.0;
        w.steepestHz = (double) e.jc * geo_.cellHz;
        w.transitionEndHz = (double) e.je * geo_.cellHz;
        w.transitionHz = w.transitionEndHz - w.cutoffHz;
        w.transitionClipped = e.clipped;
        w.truncatedAtNyquist = e.truncated;
        w.plateauPower = e.plateau;
        w.floorLocalPower = e.floorLocal;
        w.maxAbovePower = e.maxAbove;
        w.sufMaxPower = e.sufMax;
        w.exemptedCells = e.exempted;
        w.dropDb = powerRatioDb (e.plateau, e.maxAbove);
        w.strictDropDb = powerRatioDb (e.plateau, e.sufMax);
        w.localDropDb = powerRatioDb (e.plateau, e.floorLocal);
        // recovery: how far the spectrum comes back above the edge. Both floors zero means nothing is up
        // there at all and nothing recovers — the documented 0.0, rather than an inf-minus-inf NaN.
        w.recoveryDb = (! (e.sufMax > 0.0) && ! (e.floorLocal > 0.0)) ? 0.0
                     : powerRatioDb (e.sufMax, e.floorLocal);
        w.plateauSpreadDb = powerRatioDb (e.spreadHi, e.spreadLo);
        // dB per octave of the transition. A transition narrower than one cell spans no octaves at all, so
        // a real drop over it is infinitely steep — the documented value there is +inf, not a 0.0 that
        // reads as "flat", and 0/0 is never produced.
        w.steepnessDbPerOctave = (w.transitionEndHz > w.cutoffHz && w.cutoffHz > 0.0 && std::isfinite (w.dropDb))
                               ? w.dropDb / std::log2 (w.transitionEndHz / w.cutoffHz)
                               : (w.dropDb > 0.0 ? std::numeric_limits<double>::infinity() : 0.0);
        w.sharp = ! w.transitionClipped && w.dropDb >= active_.minDropDb
                && w.transitionHz <= active_.maxTransitionHz;
        w.nearNyquist = w.cutoffFractionOfNyquist >= active_.nearNyquistFraction;
    }

    // 10log10(a/b) — a power ratio, taken ONCE at the end of a linear average. b == 0 with a > 0 is
    // legitimately +inf (that IS the value); 0/0 is never asked, the callers exclude it.
    static double powerRatioDb (double a, double b) noexcept { return 10.0 * std::log10 (a / b); }

    // Sort a COPY of cells [from, to) into the scratch and answer how many landed there; every quantile of
    // that span is then an integer index into it (pick()), so one span costs one sort and the accumulator
    // is never permuted. 0 means the span is unusable, and every caller treats that as "no energy": the
    // range is empty, or — a path the sizing makes unreachable, never silently truncated — longer than the
    // scratch prepare() published.
    std::size_t loadSpan (int from, int to) noexcept
    {
        const int lo = std::max (0, from);
        const int hi = std::min (to, geo_.cellCount);
        const int count = hi - lo;
        if (count <= 0 || (std::size_t) count > sortScratch_.size()) return 0;
        const std::size_t n = (std::size_t) count;
        for (std::size_t i = 0; i < n; ++i) sortScratch_[i] = smooth_[(std::size_t) lo + i];
        std::sort (sortScratch_.begin(), sortScratch_.begin() + (std::ptrdiff_t) n);
        return n;
    }

    // NEAREST-RANK: the order statistic of 1-based rank ceil(q*n) in the loaded span, by integer
    // arithmetic. The repository's one named quantile rule (dynamics/offline/Quantile.h:22), not
    // LoudnessMeter's LRA rule and not a third one; a value that actually occurred, never a point invented
    // between two that did.
    double pick (std::size_t n, int num, int den) const noexcept
    {
        if (n == 0) return 0.0;
        const std::size_t rank = (n * (std::size_t) num + (std::size_t) den - 1u) / (std::size_t) den;
        const std::size_t idx = rank < 1u ? 0u : std::min (rank, n) - 1u;
        return sortScratch_[idx];
    }

    SourceForensicsParams params_ {}, active_ {};
    Geometry geo_ {};
    SpectrumFrames frames_ {};
    bool prepared_ = false, finished_ = false;
    double sampleRate_ = 48000.0;
    int channels_ = 0, bins_ = 0;
    std::int64_t t_ = 0;
    std::int64_t aggFramesUsed_ = 0, aggFramesHoled_ = 0;
    std::size_t tableMask_ = 0;
    Observer obs_ = nullptr;
    void* obsUser_ = nullptr;
    Channel chans_[core::kMaxChannels] {};
    SpectralWall walls_[core::kMaxChannels + 1] {};            // one per channel, plus the aggregate
    std::vector<double> sum_, comp_, mean_, cells_, smooth_, suffix_, strict_, sortScratch_;
    std::vector<std::uint32_t> table_;
};

} // namespace felitronics::analysis

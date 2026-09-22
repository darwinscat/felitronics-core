// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/core/Config.h>
#include <felitronics/core/DetMath.h>
#include <felitronics/core/Math.h>
#include <felitronics/core/StateGrid.h>
#include <felitronics/eq/Crossover2.h>
#include <felitronics/oversampling/PolyphaseOversampler.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace felitronics::analysis
{

//==================================================================================================
// felitronics::analysis::BandCrest — the peak-to-RMS of a programme, per 400 ms block and per band.
//
// WHAT IT IS FOR. An automatic mastering mode needs a MICRO-DYNAMICS damage counter: run this on the input,
// run it again on the delivered output, and the damage is the crest a band LOST, block by block. Macro damage
// is LRA (it exists); timbre is a band energy shift (it does not yet); this is the third counter.
//
// WHAT A BAND CREST CAN AND CANNOT MEAN. It is the peak-to-RMS of the signal seen through a NAMED filter. It is
// never a statement about "the content between 2 and 6 kHz": the band is a dome, and an LR4 skirt is only
// 24 dB/octave down, so a band's peak can be a neighbour's transient leaking in. In the DIFFERENCE of two runs
// through the SAME filter that leakage largely cancels; in one run's absolute number it does not, and this
// header says so rather than letting a reader assume otherwise.
//
//==================================================================================================
// THE FOUR DECISIONS, and each one is a place where the obvious choice is wrong
//==================================================================================================
//
// 1. RECONSTRUCTED PEAK, NOT SAMPLE PEAK — and the reason is not the sample rate. Even with both runs on the
//    SAME grid, a sample peak is biased in ONE direction: the input carries inter-sample peaks and the output
//    has been through a true-peak limiter whose whole job is to remove them. A sample-peak crest would
//    under-read the input's crest and therefore UNDER-REPORT the loss — the one number this class exists to
//    produce. The reconstruction is the certifying instrument's own: 4x, 32 taps per phase, the same
//    `oversampling::PolyphaseOversampler` topology `ReferenceTruePeakMeter` is built from, so the full-band
//    peak this class reports can be nulled against the certificate rather than merely believed.
//
// 2. ONE RECONSTRUCTION, THEN THE SPLIT AT 4x — not a split at the base rate followed by five
//    reconstructions. It is cheaper (one interpolation instead of five) and it makes the band magnitudes
//    nearly rate-independent as a side effect: an LR4 corner is prewarped by `tan(pi*f/fs)`, so at the base
//    rate the skirt one octave out moves with the sample rate, and a comparison across a delivery resampler
//    would read that movement as damage.
//
// 3. THE BANDS ARE INDEPENDENT FILTERS OF THE INPUT, not a reconstructing tree:
//
//        low     = LP4(f0)
//        lowMid  = HP4(f0) -> LP4(f1)
//        highMid = HP4(f1) -> LP4(f2)
//        high    = HP4(f2)
//        full    = the reconstructed stream itself
//
//    A serial tree (`eq::MultibandSplitter`) is cheaper and sums back to the input, but it pays for that with
//    ALLPASS compensation on every band — and an allpass changes the waveform, which is to say it changes the
//    PEAK, which is to say it changes the very quantity being measured. Nothing here ever needs the bands to
//    sum, so nothing here pays for summing. The bands do NOT reconstruct the input and are not meant to: at an
//    LR crossover both branches sit at 1/2 amplitude and their powers sum to 1/2, not to 1.
//
// 4. THE ACTIVITY MASK IS BUILT FROM THE INPUT ONLY. This is the decision that can quietly invalidate
//    everything else. A mask computed on each side separately, or on the output, can drop exactly the blocks
//    the chain damaged most — a loud block that the chain crushed may fall below a gate the input version
//    passed, and the damage then leaves the population along with it. So the consumer builds ONE mask, on the
//    source, and applies it to both sides; this class publishes the mask and its counts so that the population
//    behind a number is visible rather than assumed.
//
//    The gate itself is a CONJUNCTION, because each half alone fails in its own direction:
//      * an absolute programme floor (default -70 dBFS on the full band) — without it a quiet passage is
//        judged, and a relative loudness gate would throw away exactly the quiet passages where a compressor
//        does its work;
//      * a band SHARE floor (default -40 dB of the full band's energy) — without it a band with no content
//        contributes the crest of its own noise floor, and an input that is digital black in a band against an
//        output carrying dither reads as several dB of "damage" that is nothing of the kind.
//    Share is level-invariant, so it does not drop a quiet passage's minor bands the way an absolute band
//    floor would.
//
//==================================================================================================
// TIME, EXACTLY
//==================================================================================================
// The hop grid is `LoudnessMeter`'s, and it is built the way that meter builds it — out of `kSubHopsPerHop`
// sub-hops of `lround(0.01*fs)` samples, not by rounding the whole hop. The two agree wherever `0.01*fs` is an
// integer and part where it is not: at 22050 Hz a rounded hop is 2205 samples and the meter's is 2210. A block
// is `blockHops` hops (400 ms at the defaults), one block every hop once the first `blockHops` are in, so
// block j covers input samples [j*H, (j+blockHops)*H) and block j of one run pairs with block j of another.
//
// WITH ONE NAMED DIFFERENCE: this analyzer closes the hop the programme ENDS INSIDE and the meter does not, so
// on a programme that does not end on a hop boundary there is one more block here. That is deliberate — the
// last hop is the end of the file, and a master's loudest material often sits there — and it is stated rather
// than left for a consumer to discover from a count that is one off.
//
// THE INTERPOLATOR'S DELAY IS HALF AN OVERSAMPLED SAMPLE, and it is named rather than rounded away in silence.
// The polyphase FIR is symmetric over `factor*tapsPerPhase` oversampled taps, so its group delay is
// `(factor*tapsPerPhase - 1)/2` of them — MEASURED, not assumed: an impulse at base index 64 comes back
// symmetric about oversampled index 319.5, bit for bit (the suite pins it). With the shipped 4x/32 that is
// 63.5 oversampled samples. The integer part is skipped exactly; the remaining half sample is 5.2 microseconds
// at 48 kHz, 13 parts per million of a 400 ms block, and it is the same half sample on both sides of any
// comparison.
//
//==================================================================================================
// WHAT IS STORED, AND WHY IT IS LINEAR
//==================================================================================================
// One cell per hop per band: the reconstructed PEAK (linear) and the sum of squares. Blocks are derived when
// read — a block's peak is the max of its hops' peaks and its mean square the sum over its hops — and dB is a
// VIEW taken with one `core::det::log10` at the end. Storing dB would put a logarithm in the store and lose
// the one property the comparison needs: formed from linear cells, a crest LOSS is
// `20*log10((peakA*rmsB)/(rmsA*peakB))`, one logarithm and invariant to any gain applied to either side.
//
// Ten minutes at the defaults is 6000 hops x 5 bands x 16 B = 469 KB. It is not reduced online, and cannot be:
// the loss is a per-cell difference between two runs, and no single run can form it.
//==================================================================================================

struct BandCrestParams
{
    // The three corners, as a pair of Linkwitz-Riley corners each. Named in the contract, adjustable so that a
    // consumer is not forced to re-implement the analyzer to ask a different question.
    double bandEdgeHz[3] = { 120.0, 2000.0, 6000.0 };
    double hopMs         = 100.0;   // the grid's step; the block is `blockHops` of these
    int    blockHops     = 4;       // 4 x 100 ms = the 400 ms momentary window

    // THE GATE, both halves — see decision 4. In dB; turned into linear ratios once, in prepare(), so no
    // logarithm takes part in a decision.
    double programmeFloorDb = -70.0;   // the full band's mean square, absolute
    double bandShareFloorDb = -40.0;   // the band's share of the full band's energy
};

// Why a quantity has no number. `None` is the only value that means the quantity is present.
enum class BandCrestInvalid : std::uint8_t
{
    None           = 0,
    NoBlock        = 1,   // the programme never reached `blockHops` hops: there is no 400 ms window
    NonFiniteInput = 2    // a non-finite sample, or an overflowed filter, occurred in the programme
};

class BandCrest
{
public:
    // THE FILTERS ARE THE DETERMINISTIC ONES, and the alias is here so a static_assert can say so. `eq::Crossover2`
    // is `BasicCrossover2<core::SystemMath>`, whose prewarp is `std::tan` — and `std::tan` does not agree
    // between libms: a review measured Apple's against glibc/musl/UCRT at 120 Hz on 44.1 kHz and found them
    // apart, and one ulp in `g` moves the SVF's coefficients. On a 6.3-minute mix that showed up as 12 block
    // rows of 3752 differing between the two spellings. Every other offline analyzer here already took this
    // road (BandBursts, LowEnd, ProgrammeReport) and MathPolicyTests asserts it of each; this one was written
    // without it, which is the same class of defect the wasm tier caught in K11's dB threshold.
    using CrossoverType = eq::DeterministicCrossover2;   // asserted by the math-policy suite
    using SvfType       = eq::DeterministicSvf;          // what the cascades below are actually built from
    // `eq::Crossover2`'s own Q, so the cascades here ARE its cascades: 1/sqrt(2) per section, two sections,
    // which is a 4th-order Linkwitz-Riley.
    static constexpr double kCrossoverQ = 0.7071067811865476;

    // low, lowMid, highMid, high, full — `full` last so that a loop over the four BANDS is `b < kFull`.
    static constexpr int kBands = 5;
    static constexpr int kLow = 0, kLowMid = 1, kHighMid = 2, kHigh = 3, kFull = 4;

    // THE RECONSTRUCTION IS THE CERTIFICATE'S. Both constants are `ReferenceTruePeakMeter`'s, spelled here so
    // that a test can assert they are the same pair rather than hope so: the full-band peak is then the same
    // quantity the certifying meter reports, and can be nulled against it.
    static constexpr int kFactor       = 4;
    static constexpr int kTapsPerPhase = 32;

    // The interpolator's group delay in OVERSAMPLED samples, doubled so it is an integer: the FIR is symmetric
    // over `kFactor*kTapsPerPhase` taps, so the delay is `(kFactor*kTapsPerPhase - 1)/2` = 63.5 at the shipped
    // topology. Kept as a doubled integer because the half is real and rounding it away silently is how an
    // instrument acquires a bias nobody can find later.
    static constexpr int kDelayOsX2 = kFactor * kTapsPerPhase - 1;

    static constexpr int kChunk = 1024;                 // the interpolator's scratch, in base samples

    // The corner must clear the oversampled Nyquist with the same margin `eq::Svf` is honest over.
    static constexpr double kMaxCornerFraction = 0.49;

    //----------------------------------------------------------------------------------------------
    // WHAT prepare() ALLOCATES, from the function prepare() sizes and validates itself with (law 11d).
    struct Storage
    {
        bool        ok          = false;
        int         hopSamples  = 0;
        std::size_t hopCapacity = 0;                    // hops the cell store can hold
        std::size_t cellEntries = 0;                    // hopCapacity * kBands * 2 doubles
        std::size_t countEntries = 0;                   // one sample count per hop — the exact denominator
        std::size_t scratchOs   = 0;                    // the interpolator's output, oversampled floats
        std::size_t channels    = 0;

        std::uint64_t bytes() const noexcept
        {
            return (std::uint64_t) cellEntries  * sizeof (double)
                 + (std::uint64_t) countEntries * (sizeof (std::uint64_t) + 2u * sizeof (double))
                 + (std::uint64_t) scratchOs   * sizeof (float)
                 // THE INTERPOLATORS THEMSELVES, not only their buffers: `os_.resize(maxChannels)` asks the
                 // heap for `maxChannels * sizeof(PolyphaseOversampler)` and that was missing from the demand,
                 // which under-reported by 168 B a channel at EVERY length. A budget is a promise about what
                 // the call asks for, so the container it asks for counts.
                 + (std::uint64_t) channels    * ((std::uint64_t) oneOversamplerBytes()
                                                  + sizeof (oversampling::PolyphaseOversampler));
        }

        static std::uint64_t oneOversamplerBytes() noexcept
        {
            oversampling::PolyphaseOversampler::Storage one {};
            if (! oversampling::PolyphaseOversampler::storageFor (kFactor, 1, kTapsPerPhase, one)) return 0;
            return one.bytes();
        }
    };

    // `maxSamples` is the programme this preparation can hold cells for. Blocks past it are COUNTED and not
    // kept (`droppedHops()`), the rule `LoudnessMeter` already follows — an instrument that silently measured
    // the first part of a programme and reported it as the whole is the failure this counts instead.
    // The rate ceiling every sibling analyzer carries (BandBursts, ClipDetector, HumDetector, LowEnd,
    // ProgrammeReport, SourceForensics). This one was written without it, and the consequence was not merely a
    // silly number being accepted: `llround` SATURATES on arm64 and WRAPS on x86-64 glibc, so `fs = 1e308` was
    // refused on one row and accepted on another with a 10-sample hop — a refusal set that differs between
    // platforms, which is half a parity contract gone. Measured on both rows before this was written.
    static constexpr double kMaxSampleRate = 768000.0;
    static constexpr double kMaxHopMs      = 10000.0;    // ten seconds of hop is already absurd; the point is
                                                         // that the bound exists BEFORE `llround` sees it

    static Storage storageFor (double sampleRate, int maxChannels, long long maxSamples,
                               const BandCrestParams& p) noexcept
    {
        Storage st;
        if (! (sampleRate >= core::kMinSampleRate) || ! (sampleRate <= kMaxSampleRate)) return st;   // NaN fails both
        if (maxChannels < 1 || maxChannels > core::kMaxChannels) return st;
        if (maxSamples < 0) return st;
        // BOUNDED BEFORE THE CONVERSION, not after: `llround(1e299/10)` is undefined, and what it does in
        // practice differs by row. A comparison, not a conversion, is what makes this portable.
        // AND THE FLOOR IS THE QUANTUM, not 1 ms. A hop is built from 10 ms sub-hops, so `hopMs` is rounded to
        // a multiple of 10 and anything under 10 becomes 10 — which means a caller passing 1 would read its own
        // 1 back out of `params()` while the measurement ran at 10. The parameter refuses rather than echoes a
        // number that describes nothing.
        if (! (p.hopMs >= 10.0) || ! (p.hopMs <= kMaxHopMs)) return st;
        if (p.blockHops < 1 || p.blockHops > 64) return st;
        if (! cornersAdmitted (sampleRate, p)) return st;

        // THE HOP IS THE LOUDNESS METER'S OWN ARITHMETIC, not a formula that agrees with it at the rates
        // anyone happened to test. It builds a hop out of `kSubHopsPerHop` sub-hops of `lround(0.01*fs)`
        // samples; rounding the whole hop instead gives a different integer wherever 0.01*fs is not one —
        // at 22050 Hz that is 2210 samples against 2205, and the header's claim that block j is the meter's
        // block j would have been false there. A review found it; the test that claimed the identity ran
        // only at 44.1, 48 and 96 kHz, where the two agree.
        const double subHop = 0.01 * sampleRate;
        const long long sub = std::max (1LL, (long long) std::llround (subHop));
        const long long hopsPerStep = std::max (1LL, (long long) std::llround (p.hopMs / 10.0));
        const long long hop = sub * hopsPerStep;
        if (hop < 1 || hop > 0x7FFFFFFFLL) return st;
        // The hops a programme of `maxSamples` can close, plus the partial one `finish()` closes.
        const long long hops = maxSamples / hop + 1;
        if (hops < 1 || hops > (long long) (1u << 28)) return st;

        st.ok          = true;
        st.hopSamples  = (int) hop;
        st.hopCapacity = (std::size_t) hops;
        st.cellEntries  = (std::size_t) hops * (std::size_t) kBands * 2u;
        st.countEntries = (std::size_t) hops;   // and, at the same count, the base-rate sample peaks
        // THE SCRATCH IS PER CHANNEL, and that is forced by the hop clock rather than chosen: a hop boundary
        // falls inside a chunk, so the accumulation has to walk the oversampled samples IN ORDER with every
        // channel available at each index. Interpolating one channel at a time into one buffer would have
        // attributed every sample after a boundary to the hop before it — up to `kChunk` of them.
        st.scratchOs   = (std::size_t) kChunk * (std::size_t) kFactor * (std::size_t) maxChannels;
        st.channels    = (std::size_t) maxChannels;
        return st;
    }

    // Every corner in order, each clearing the OVERSAMPLED Nyquist — which is the rate the filters run at, and
    // the reason a 6 kHz corner is admitted at base rates a base-rate split would refuse.
    static bool cornersAdmitted (double sampleRate, const BandCrestParams& p) noexcept
    {
        const double osRate = sampleRate * (double) kFactor;
        double last = 0.0;
        for (const double f : p.bandEdgeHz)
        {
            if (! std::isfinite (f) || ! (f > last)) return false;
            if (! (f <= kMaxCornerFraction * osRate)) return false;
            last = f;
        }
        return true;
    }

    void setParams (const BandCrestParams& p) noexcept { pending_ = p; }   // taken at the next prepare()
    // WHAT THE LAST prepare() TOOK, not what the next one would: a reader asking an object what it measured
    // with must not be told what it is about to measure with. `pendingParams()` answers the other question.
    const BandCrestParams& params() const noexcept { return p_; }
    const BandCrestParams& pendingParams() const noexcept { return pending_; }

    [[nodiscard]] bool prepare (double sampleRate, int maxChannels, long long maxSamples) noexcept
    {
        // LAW 11(b): DISARM FIRST, validate, then write. A refused prepare() leaves the object unusable rather
        // than leaving the previous preparation standing and answering process() calls.
        // LAW 11(b), AND ITS OTHER HALF. Disarming `prepared_` stops `process()` and stopped nothing else: a
        // refused prepare() left the PREVIOUS run's blocks, crests and validity readable, and `bandCrestLoss`
        // would happily compare them. The header claimed the object was unusable; only one of its doors was
        // shut. A done-flag is not a validity flag — every reader is keyed on this now.
        prepared_ = false;
        hops_ = 0; dropped_ = 0; samples_ = 0; baseHops_ = 0; finished_ = false;
        nonFinite_ = 0; firstNonFiniteAt_ = -1; narrowed_ = 0; widest_ = 0; ranNc_ = 0;
        const Storage st = storageFor (sampleRate, maxChannels, maxSamples, pending_);
        if (! st.ok) return false;

        fs_         = sampleRate;
        p_          = pending_;
        hopSamples_ = st.hopSamples;
        channels_   = maxChannels;
        hopCap_     = st.hopCapacity;

        cells_.assign (st.cellEntries, 0.0);
        counts_.assign (st.countEntries, 0u);
        basePeak_.assign (st.countEntries, 0.0);
        baseMs_.assign   (st.countEntries, 0.0);
        scratch_.assign (st.scratchOs, 0.0f);
        os_.resize ((std::size_t) maxChannels);
        for (auto& o : os_) if (! o.prepare (kFactor, 1, kTapsPerPhase)) return false;

        // SIX CASCADES, NOT FIVE CROSSOVERS. A `Crossover2` computes BOTH of its outputs and has no one-sided
        // form, so four of the five here threw one away: 8 of 20 filter evaluations a sample were work nobody
        // read. The cascades below are the SAME arithmetic — `Crossover2` is exactly `lp2(lp1(x))` and
        // `hp2(hp1(x))` at Q = 1/sqrt(2), with no correction term — so this is bit-identical and 40 % cheaper,
        // which matters because an automatic search runs this once per candidate. The suite nulls each band
        // against a cascade it builds itself, so the identity is asserted rather than asserted to be obvious.
        const double osRate = fs_ * (double) kFactor;
        const auto arm = [&] (SvfType& a, SvfType& b2, eq::FilterType t, double hz)
        {
            // THE ARGUMENTS ARE `Crossover2`'s, TO THE TYPE. It stores its corner as a float and then passes
            // `kQ` as a DOUBLE; rounding Q to float here made the coefficients differ and the rework stopped
            // being the same filter — caught by diffing the output against the previous build rather than by
            // reasoning about it, which is the only way this class of claim is worth anything.
            a.prepare (osRate, maxChannels);  a.setParams  (t, (double) (float) hz, kCrossoverQ, 0.0);
            b2.prepare (osRate, maxChannels); b2.setParams (t, (double) (float) hz, kCrossoverQ, 0.0);
        };
        arm (lo0a_, lo0b_, eq::FilterType::LowPass,  p_.bandEdgeHz[0]);   // band `low`
        arm (hi0a_, hi0b_, eq::FilterType::HighPass, p_.bandEdgeHz[0]);   // ... and the head of `lowMid`
        arm (m1a_,  m1b_,  eq::FilterType::LowPass,  p_.bandEdgeHz[1]);   // band `lowMid`
        arm (hi1a_, hi1b_, eq::FilterType::HighPass, p_.bandEdgeHz[1]);   // the head of `highMid`
        arm (m2a_,  m2b_,  eq::FilterType::LowPass,  p_.bandEdgeHz[2]);   // band `highMid`
        arm (hi2a_, hi2b_, eq::FilterType::HighPass, p_.bandEdgeHz[2]);   // band `high`

        // THE GATES BECOME RATIOS HERE, so that no logarithm takes part in a decision later.
        floorMs_    = core::det::pow10 (p_.programmeFloorDb / 10.0);
        shareRatio_ = core::det::pow10 (p_.bandShareFloorDb / 10.0);

        prepared_ = true;
        reset();
        return true;
    }

    bool isPrepared() const noexcept { return prepared_ && ! cells_.empty(); }

    void reset() noexcept
    {
        std::fill (cells_.begin(), cells_.end(), 0.0);
        std::fill (counts_.begin(), counts_.end(), 0u);
        std::fill (basePeak_.begin(), basePeak_.end(), 0.0);
        std::fill (baseMs_.begin(), baseMs_.end(), 0.0);
        for (SvfType* f : filters()) f->reset();
        for (auto& o : os_) o.reset();
        grid_.reset();
        hops_ = 0; dropped_ = 0; samples_ = 0; nonFinite_ = 0; firstNonFiniteAt_ = -1;
        osSkipped_ = 0; osInHop_ = 0; osMeasured_ = 0; ranNc_ = 0; finished_ = false; curCount_ = 0;
        narrowed_ = 0; widest_ = 0;
        baseHops_ = 0; baseInHop_ = 0; baseCur_ = 0.0; baseSum_ = 0.0; baseCount_ = 0;
        for (int b = 0; b < kBands; ++b) { curPeak_[b] = 0.0; curSum_[b] = 0.0; }
    }

    //----------------------------------------------------------------------------------------------
    [[nodiscard]] bool process (const float* const* in, int numChannels, int n) noexcept
    {
        if (! prepared_ || finished_) return false;
        if (numChannels < 0 || numChannels > channels_ || n < 0) return false;
        if (n == 0) return true;
        if (in == nullptr) return false;
        for (int c = 0; c < numChannels; ++c) if (in[c] == nullptr) return false;
        if (numChannels > widest_) widest_ = numChannels;
        if (numChannels < widest_) narrowed_ += n;
        ranNc_ = numChannels;

        for (int off = 0; off < n; )
        {
            const int m = std::min (kChunk, n - off);
            // Interpolate every channel first, then walk the oversampled samples IN ORDER. The channel loop is
            // inside the sample loop, which is this repository's idiom and is also what the hop clock needs.
            for (int c = 0; c < numChannels; ++c)
            {
                const float* src[1] { in[c] + off };
                float*       dst[1] { planeFor (c) };
                os_[(std::size_t) c].upsample (src, 1, m, dst);
            }
            // THE SAMPLE-PEAK FLOOR, on its own clock. The interpolator reads UNDER the truth on a full-band
            // click — `ReferenceTruePeakMeter` says so and carries the same floor — and an under-read peak on
            // the INPUT side under-reports the crest it lost, which is the one error this class is built to
            // avoid. The base samples cannot be recovered from the oversampled stream (the FIR's half-sample
            // delay means no output phase reproduces them), so they are measured here, against a hop clock of
            // their own. Both clocks partition the SAME input time, so hop h is hop h in either; `finish()`
            // checks that they agree rather than assuming it.
            for (int i = 0; i < m; ++i)
            {
                float g = 0.0f;
                for (int c = 0; c < numChannels; ++c)
                {
                    const float a = std::fabs (in[c][off + i]);
                    if (std::isfinite (a) && a > g) g = a;
                }
                if ((double) g > baseCur_) baseCur_ = (double) g;
                for (int c = 0; c < numChannels; ++c)
                {
                    const double x = (double) in[c][off + i];
                    if (std::isfinite (x)) { baseSum_ += x * x; ++baseCount_; }
                }
                if (++baseInHop_ >= hopSamples_)
                {
                    baseInHop_ = 0;
                    if ((std::size_t) baseHops_ < hopCap_)
                    {
                        basePeak_[(std::size_t) baseHops_] = baseCur_;
                        baseMs_[(std::size_t) baseHops_]   = baseCount_ > 0 ? baseSum_ / (double) baseCount_ : 0.0;
                    }
                    ++baseHops_;
                    baseCur_ = 0.0; baseSum_ = 0.0; baseCount_ = 0;
                }
            }
            walkOs (m * kFactor, numChannels);
            samples_ += m;
            off      += m;
        }
        return true;
    }

    // Drain the interpolator and close the partial hop. After this the store describes the whole programme and
    // `process()` refuses — a measurement that could be appended to is a measurement of nothing in particular.
    void finish() noexcept
    {
        if (! prepared_ || finished_) return;
        // The interpolator still holds `kTapsPerPhase` base samples of programme. Push silence through it, as
        // `ReferenceTruePeakMeter::drain` does, so the last transient's reconstruction is measured rather than
        // left in the ring.
        const float zeros[kTapsPerPhase] {};
        for (int c = 0; c < ranNc_; ++c)
        {
            const float* src[1] { zeros };
            float*       dst[1] { planeFor (c) };
            os_[(std::size_t) c].upsample (src, 1, kTapsPerPhase, dst);
        }
        // THE DRAIN CARRIES THE PROGRAMME'S TAIL, NOT NEW TIME. The interpolator lags by `kDelayOs`, so at the
        // end of `process()` exactly that many oversampled samples of the programme are still inside it; the
        // drain pushes them out. Taking the WHOLE drain would measure `kTapsPerPhase` base samples of silence
        // as programme and close a hop that never happened — the two clocks caught it, 31 hops against 30.
        // The window is the programme's own: `samples_ * kFactor` measurable samples, no more.
        const long long room = (long long) samples_ * (long long) kFactor - osMeasured_;
        walkOs (kTapsPerPhase * kFactor, ranNc_, room > 0 ? room : 0);
        // THE PARTIAL LAST HOP IS AN ENTRY LIKE ANY OTHER, closed over its own length — the rule the loudness
        // meter follows. `osInHop_ > 0` is the whole test: a hop that closed exactly on the drain's last
        // sample has already been stored and must not be stored again as an empty one.
        if (osInHop_ > 0) closeHop();
        // The base clock's partial hop, closed like the oversampled one.
        if (baseInHop_ > 0)
        {
            if ((std::size_t) baseHops_ < hopCap_)
            {
                basePeak_[(std::size_t) baseHops_] = baseCur_;
                baseMs_[(std::size_t) baseHops_]   = baseCount_ > 0 ? baseSum_ / (double) baseCount_ : 0.0;
            }
            ++baseHops_;
            baseInHop_ = 0; baseCur_ = 0.0; baseSum_ = 0.0; baseCount_ = 0;
        }
        finished_ = true;
    }

    //----------------------------------------------------------------------------------------------
    // READING. Everything below is derived from the cells; nothing is stored in dB.
    double       sampleRate() const noexcept { return fs_; }
    int          hopSamples() const noexcept { return hopSamples_; }
    int          blockHops()  const noexcept { return p_.blockHops; }
    int          channels()   const noexcept { return ranNc_; }
    long long    hopCount()   const noexcept { return hops_; }
    // The base-rate clock's own hop count. It must equal `hopCount()`; a test asserts it rather than a comment
    // claiming it, because the two are counted by different code over the same time.
    long long    basePeakHops() const noexcept { return baseHops_; }
    long long    droppedHops() const noexcept { return dropped_; }
    long long    samplesProcessed() const noexcept { return samples_; }
    long long    nonFiniteSamples() const noexcept { return nonFinite_; }

    // A CHANNEL THAT VANISHES MID-PROGRAMME CHANGES THE MEASUREMENT, and silently unless it is counted. The
    // peak is a maximum over the channels PRESENT and the mean square is divided by the samples actually
    // accumulated, so a stretch fed at one channel where the rest was fed at two is measured over a different
    // width — the numbers stay finite and plausible and nothing in them says the width moved. This is the
    // count: base samples processed at fewer channels than the widest this run has seen. `BandBursts` counts
    // its absent channels for the same reason and this is that rule kept.
    long long    narrowedSamples() const noexcept { return narrowed_; }
    int          widestChannels() const noexcept { return widest_; }
    long long    firstNonFiniteAt() const noexcept { return firstNonFiniteAt_; }

    long long blockCount() const noexcept
    {
        if (! prepared_) return 0;                     // a refused preparation answers nothing, not the last run
        const long long b = hops_ - (long long) p_.blockHops + 1;
        return b > 0 ? b : 0;
    }

    BandCrestInvalid invalidReason() const noexcept
    {
        if (nonFinite_ > 0) return BandCrestInvalid::NonFiniteInput;
        if (blockCount() == 0) return BandCrestInvalid::NoBlock;
        return BandCrestInvalid::None;
    }
    bool valid() const noexcept { return invalidReason() == BandCrestInvalid::None; }

    // Block j covers input samples [j*H, (j+blockHops)*H). Out of range answers 0 rather than reading memory.
    double blockPeakLin (long long block, int band) const noexcept
    {
        if (! inRange (block, band)) return 0.0;
        double pk = 0.0;
        for (int k = 0; k < p_.blockHops; ++k) pk = std::max (pk, cellPeak (block + k, band));
        return pk;
    }

    // THE FULL BAND'S MEAN SQUARE IS THE BASE-RATE ONE, and that is a correction rather than a choice. Its
    // PEAK is the reconstruction (max of the sample peak and the 4x one, the certifying meter's own quantity),
    // but its mean square used to come from the oversampled stream — and the 32-tap interpolator droops above
    // about 0.375*fs. The two together made a HYBRID: a pure 0.45*fs tone read a crest of 9.03 dB where the
    // truth is 3.0103, and 23.8 dB at 0.49*fs. Measured, not derived. The per-band crests were never affected
    // — their peak and their mean square come from the same stream — so only this one number was wrong, and it
    // fed the absolute gate and `programmeMeanSquareDb` with it.
    //
    // The SHARE gate still divides by the oversampled full-band power (`osMeanSq`), because a band's power is
    // an oversampled quantity and a ratio across two domains is not a share of anything.
    double blockMeanSq (long long block, int band) const noexcept
    {
        if (! inRange (block, band)) return 0.0;
        if (band == kFull)
        {
            double s = 0.0; long long have = 0;
            for (int k = 0; k < p_.blockHops; ++k)
                if (block + k < baseHops_) { s += baseMs_[(std::size_t) (block + k)]; ++have; }
            return have > 0 ? s / (double) have : 0.0;
        }
        // THE DENOMINATOR IS COUNTED, NOT COMPUTED. The last hop is partial — the drain closes it over its own
        // length — so a denominator derived from the hop length would understate exactly one block's mean
        // square, and it would be the block at the end of the programme, where a master's loudest material
        // often sits. The count is what was actually accumulated, channels included.
        double s = 0.0;
        std::uint64_t n = 0;
        for (int k = 0; k < p_.blockHops; ++k) { s += cellSum (block + k, band); n += counts_[(std::size_t) (block + k)]; }
        return n > 0u ? s / (double) n : 0.0;
    }

    // ONE LOGARITHM, from the linear cell. 20log10(peak) - 10log10(meanSq), formed as a single ratio so the
    // two halves cannot be rounded apart. Silence answers the floor rather than -inf.
    double blockCrestDb (long long block, int band) const noexcept
    {
        const double pk = blockPeakLin (block, band), ms = blockMeanSq (block, band);
        if (! (pk > 0.0) || ! (ms > 0.0)) return 0.0;
        return 10.0 * core::det::log10 ((pk * pk) / ms);
    }

    // THE MASK, and it is a property of THIS run — the consumer builds it on the source and applies it to both
    // sides (see decision 4). Both halves of the conjunction, so a caller can see which one refused.
    bool blockActive (long long block, int band) const noexcept
    {
        if (! inRange (block, band)) return false;
        if (! (blockMeanSq (block, kFull) >= floorMs_)) return false;      // the ABSOLUTE gate, base-rate
        if (band == kFull) return true;
        const double os = osMeanSq (block);                                 // the SHARE gate, one domain
        return os > 0.0 && blockMeanSq (block, band) >= shareRatio_ * os;
    }

    long long activeBlocks (int band) const noexcept
    {
        long long n = 0;
        for (long long j = 0, e = blockCount(); j < e; ++j) if (blockActive (j, band)) ++n;
        return n;
    }

    // THE PROGRAMME'S OWN LEVEL, in the units `programmeFloorDb` is compared against — so a caller can set a
    // RELATIVE floor without converting between two quantities that are not the same one.
    //
    // WHY IT EXISTS. A caller wanting "40-odd dB below the programme" naturally reaches for integrated
    // loudness, and `I - 42` is the obvious spelling. But `I` is K-weighted and this gate is not: measured on
    // three 20 s fixtures, `I` minus this number runs from -0.15 dB on bass-heavy material to +5.87 dB on
    // bright material — the offset is a function of the SPECTRUM, so a floor derived that way is a floor that
    // moves with the mix, and a calibration made on one kind of material would silently be wrong on another.
    // This is the same quantity, by the same code path, so the difference is zero by construction.
    //
    // THE GATE HERE IS A FIXED -70 dBFS, NOT `programmeFloorDb`, and that is the point rather than a detail:
    // a caller derives its floor FROM this number, so averaging over the population that floor selects would
    // be a threshold chasing its own tail. Two stages, an absolute gate and then a relative one computed from
    // what it admitted, is exactly BS.1770's own construction for the integrated measure — this is not a new
    // idea, it is that idea in unweighted units.
    //
    // 0 blocks above the gate answers `kSilenceDb`, which is a sentinel and not a level.
    static constexpr double kProgrammeGateDb = -70.0;
    static constexpr double kSilenceDb       = -200.0;

    double programmeMeanSquareDb() const noexcept
    {
        const double gate = core::det::pow10 (kProgrammeGateDb / 10.0);
        double sum = 0.0; long long n = 0;
        for (long long j = 0, e = blockCount(); j < e; ++j)
        {
            const double ms = blockMeanSq (j, kFull);
            if (std::isfinite (ms) && ms >= gate) { sum += ms; ++n; }
        }
        if (n == 0) return kSilenceDb;
        return 10.0 * core::det::log10 (sum / (double) n);
    }

    // The whole programme's reconstructed peak, which is the quantity `ReferenceTruePeakMeter` certifies — the
    // one number here with an instrument OUTSIDE this class to be nulled against.
    double fullBandPeakLin() const noexcept
    {
        double pk = 0.0;
        for (long long h = 0; h < hops_; ++h) pk = std::max (pk, cellPeak (h, kFull));
        return pk;
    }

private:
    bool inRange (long long block, int band) const noexcept
    {
        return band >= 0 && band < kBands && block >= 0 && block < blockCount();
    }
    double cellPeak (long long hop, int band) const noexcept
    {
        const double os = cells_[(std::size_t) ((hop * kBands + band) * 2 + 0)];
        // The floor applies to the FULL band only: a band signal exists at the oversampled rate and nowhere
        // else, so there is no sample of it to be a floor.
        if (band != kFull || hop >= baseHops_) return os;
        const double sp = basePeak_[(std::size_t) hop];
        return sp > os ? sp : os;
    }
    double cellSum (long long hop, int band) const noexcept
    {
        return cells_[(std::size_t) ((hop * kBands + band) * 2 + 1)];
    }

    // The full band's OVERSAMPLED mean square — the share gate's denominator, in the bands' own domain.
    double osMeanSq (long long block) const noexcept
    {
        double s = 0.0; std::uint64_t n = 0;
        for (int k = 0; k < p_.blockHops; ++k)
        { s += cellSum (block + k, kFull); n += counts_[(std::size_t) (block + k)]; }
        return n > 0u ? s / (double) n : 0.0;
    }

    // One channel's plane inside the shared oversampled scratch.
    float* planeFor (int c) noexcept
    {
        return scratch_.data() + (std::size_t) c * (std::size_t) (kChunk * kFactor);
    }

    // The oversampled samples of one chunk, IN ORDER, every channel at each index. The bands are filtered for
    // EVERY sample, the interpolator's ramp-up included — the filters' state has to be right when the
    // measurement starts — and accumulated only once the integer part of the delay has gone by.
    // `limit` caps how many MEASURABLE oversampled samples this call may take; -1 is no cap.
    void walkOs (int osSamples, int nch, long long limit = -1) noexcept
    {
        const long long hopOs = (long long) hopSamples_ * (long long) kFactor;
        long long taken = 0;
        for (int i = 0; i < osSamples; ++i)
        {
            bool measuring = (osSkipped_ >= kDelayOs);
            if (! measuring) ++osSkipped_;
            else if (limit >= 0 && taken >= limit) measuring = false;   // past the programme's own time

            for (int c = 0; c < nch; ++c)
            {
                const float x = planeFor (c)[i];
                // NON-FINITE IS GATED, NOT FLUSHED: an infinity in a recursive filter is not healable, and a
                // hop it poisoned is not a measurement. It is replaced by silence, counted, and located.
                float v = x;
                if (! std::isfinite (v))
                {
                    v = 0.0f;
                    if (nonFinite_ != ~0LL) ++nonFinite_;
                    if (firstNonFiniteAt_ < 0) firstNonFiniteAt_ = samples_;
                }
                const float lo0 = lo0b_.processSample (c, lo0a_.processSample (c, v));     // LP4(f0) = `low`
                const float hi0 = hi0b_.processSample (c, hi0a_.processSample (c, v));     // HP4(f0)
                const float m1  = m1b_ .processSample (c, m1a_ .processSample (c, hi0));   // -> LP4(f1) = `lowMid`
                const float hi1 = hi1b_.processSample (c, hi1a_.processSample (c, v));     // HP4(f1)
                const float m2  = m2b_ .processSample (c, m2a_ .processSample (c, hi1));   // -> LP4(f2) = `highMid`
                const float hi2 = hi2b_.processSample (c, hi2a_.processSample (c, v));     // HP4(f2) = `high`
                if (! measuring) continue;

                const float band[kBands] { lo0, m1, m2, hi2, v };
                for (int b = 0; b < kBands; ++b)
                {
                    const float a = std::fabs (band[b]);
                    // A FILTER CAN OVERFLOW FROM A FINITE INPUT (`eq::Svf` narrows its state to float), so the
                    // OUTPUT is gated too — the rule `BandBursts` applies, and for the same reason.
                    if (! std::isfinite (a)) { if (nonFinite_ != ~0LL) ++nonFinite_; continue; }
                    if ((double) a > curPeak_[b]) curPeak_[b] = (double) a;
                    curSum_[b] += (double) a * (double) a;
                }
                ++curCount_;
            }

            if (measuring)
            {
                ++taken; ++osMeasured_;
                if (++osInHop_ >= hopOs) { osInHop_ = 0; closeHop(); }
            }
            // Denormal maintenance on the OVERSAMPLED clock, which is this bank's own clock — ONE SAMPLE AT A
            // TIME, which is `core::StateGrid::advance`'s precondition and what every other analyzer here
            // does. Advancing it by a whole chunk made the flush land at a different point for different call
            // sizes, and two cells then differed between a 1-sample call and a 4096-sample one.
            if (grid_.advance (1)) flushAll();
        }
    }

    void closeHop() noexcept
    {
        if ((std::size_t) hops_ < hopCap_)
        {
            for (int b = 0; b < kBands; ++b)
            {
                cells_[(std::size_t) ((hops_ * kBands + b) * 2 + 0)] = curPeak_[b];
                cells_[(std::size_t) ((hops_ * kBands + b) * 2 + 1)] = curSum_[b];
            }
            counts_[(std::size_t) hops_] = curCount_;
            ++hops_;
        }
        else ++dropped_;                                   // past the capacity: counted, not kept
        for (int b = 0; b < kBands; ++b) { curPeak_[b] = 0.0; curSum_[b] = 0.0; }
        curCount_ = 0;
    }

    // The integer part of the interpolator's delay, in OVERSAMPLED samples, rounded DOWN. The half that is
    // left over is the residual the header names — 63 of 63.5 — and it is the same half on both sides of any
    // comparison, so it cannot become a difference.
    static constexpr int kDelayOs = kDelayOsX2 / 2;

    BandCrestParams pending_ {}, p_ {};
    bool      prepared_ = false, finished_ = false;
    double    fs_ = 0.0;
    int       hopSamples_ = 0, channels_ = 0, ranNc_ = 0;
    std::size_t hopCap_ = 0;
    long long hops_ = 0, dropped_ = 0, samples_ = 0, nonFinite_ = 0, firstNonFiniteAt_ = -1;
    int           osSkipped_ = 0;
    long long     osInHop_ = 0, osMeasured_ = 0;
    std::uint64_t curCount_ = 0;      // the samples accumulated into the hop being built, all channels
    long long     baseHops_ = 0, baseInHop_ = 0, narrowed_ = 0;
    int           widest_ = 0;
    double        baseCur_ = 0.0, baseSum_ = 0.0;
    long long     baseCount_ = 0;
    double    curPeak_[kBands] {}, curSum_[kBands] {};
    double    floorMs_ = 0.0, shareRatio_ = 0.0;

    std::vector<double>        cells_, basePeak_, baseMs_;
    std::vector<std::uint64_t> counts_;
    std::vector<float>  scratch_;
    std::vector<oversampling::PolyphaseOversampler> os_;
    SvfType lo0a_, lo0b_, hi0a_, hi0b_, m1a_, m1b_, hi1a_, hi1b_, m2a_, m2b_, hi2a_, hi2b_;

    std::array<SvfType*, 12> filters() noexcept
    {
        return { &lo0a_, &lo0b_, &hi0a_, &hi0b_, &m1a_, &m1b_, &hi1a_, &hi1b_, &m2a_, &m2b_, &hi2a_, &hi2b_ };
    }
    void flushAll() noexcept { for (SvfType* f : filters()) f->flushDenormals(); }
    core::StateGrid grid_;
};

//==================================================================================================
// THE PAIRED COMPARISON — what K1 exists to produce. A per-run crest table is descriptive; the DAMAGE is the
// crest a band lost between the source and the master, block by block.
//
// THE POPULATION IS THE SOURCE'S MASK, and only the source's. A mask taken from the output, or from both, can
// drop exactly the blocks the chain damaged most: a loud block the chain crushed can fall below a gate its
// input passed, and the damage then leaves the population with it. So the caller passes the source as `in` and
// the master as `out`, and the mask comes from `in` alone.
//
// THE LOSS IS ONE LOGARITHM OF A RATIO OF LINEAR CELLS:
//
//     loss(j) = 10*log10( (peakIn^2 * meanSqOut) / (peakOut^2 * meanSqIn) )
//
// which is `crestIn - crestOut` written so that no pair of dB values is subtracted. It is INVARIANT to any gain
// applied to either side — a master is louder than its source by construction, and a loss that moved with the
// gain would measure the gain.
//
// FOUR NUMBERS, NOT ONE, and the reason is that the obvious one is blind to the instrument doing the damage.
// A limiter works on the loudest few percent of blocks. If 4 % of the blocks lose 6 dB of crest, the p95 of the
// positive part reads ZERO — the statistic is 0 exactly where the damage is concentrated. So:
//   * `cvar95` — the mean of the worst 5 % — is the number that sees concentrated damage, and is meant to be
//     the one a cost function reads;
//   * `p95` stays beside it, because it is what a reader expects and its disagreement with `cvar95` is itself
//     the signal that the damage is concentrated;
//   * `over1Db` / `over3Db` / `over6Db` count the blocks past each threshold, which is the same fact as a count
//     rather than a quantile;
//   * `p5` is the NEGATIVE side, and it is not symmetry for its own sake: a limiter that pulls a kick down for
//     tens of milliseconds lowers a band's RMS while its peak elsewhere in the block is untouched, so the crest
//     RISES. Clipped to its positive part that damage is invisible; `p5` is where it shows.
//
// A crest loss with no PEAK loss is a floor RISE, not a flattening — a clipper putting a steady harmonic into a
// band that held only a click reads the same way here. The band's own peak and mean-square shifts are published
// so a consumer can tell the two apart rather than attribute both to micro-dynamics.
struct BandCrestLoss
{
    long long blocks   = 0;     // blocks the two runs have in common
    long long inActive = 0;     // ... of which the SOURCE's mask admits, which is the population
    long long usable   = 0;     // ... of which both sides carry a finite peak and energy
    long long outSilent = 0;    // ... and of the admitted ones, those whose MASTER band is digital silence.
                                // They are not "no data": they are an infinite loss, and folding them into a
                                // quantile would put an infinity in it while dropping them silently makes a
                                // muted band read as a clean one. Counted here, decided by the consumer.
    long long lagBlocks = 0;    // the block lag at which the two full-band peak series agree best. NOT a
                                // guarantee of alignment: it sees whole blocks, so a delay of a few
                                // milliseconds reads 0 — and a few milliseconds is enough to fabricate damage
                                // (a review measured 5 ms as up to 1.1 dB of CVaR95 and a 17 dB maximum, the
                                // order of a real limiter's). Alignment is the caller's precondition; this is
                                // the coarse half of it, published so a GROSS mistake cannot pass unseen.
    double    p50Db = 0.0, p95Db = 0.0, cvar95Db = 0.0, meanDb = 0.0, maxDb = 0.0;   // of max(0, loss)
    double    p5Db  = 0.0;      // the signed loss at the 5th percentile: negative means crest was GAINED
    long long over1Db = 0, over3Db = 0, over6Db = 0;
    double    peakShiftDb = 0.0, levelShiftDb = 0.0;   // median shifts, so a floor rise can be told apart
    bool      valid = false;    // false => every number above is a placeholder, not a measurement
};

// `scratch` is the caller's, and it is an argument rather than a local so that the allocation is visible: it
// needs `min(in.blockCount(), out.blockCount())` doubles and is resized here if it is short.
inline BandCrestLoss bandCrestLoss (const BandCrest& in, const BandCrest& out, int band,
                                    std::vector<double>& scratch)
{
    BandCrestLoss r;
    if (band < 0 || band >= BandCrest::kBands) return r;
    // THE TWO RUNS MUST BE COMPARABLE, and nothing upstream checks it. A loss is a per-block difference, so
    // two runs on different grids are two different questions: a review fed the same audio to both slots with
    // the master on a 50 ms hop and got CVaR95 7.8 dB and `valid` true — a plausible number for a comparison
    // that never happened. A mismatch is a refusal here, not a number with a caveat.
    // `core::exactlyEqual` because the comparison IS exact here and says so: two runs are the same question
    // only if they were asked with the same numbers, and a tolerance would let a grid that differs in its last
    // bit pass as the same grid. This is the one place a float comparison is not a mistake, and the helper is
    // how this tree spells that.
    if (! core::exactlyEqual (in.sampleRate(), out.sampleRate())
        || in.hopSamples() != out.hopSamples() || in.blockHops() != out.blockHops()
        || in.channels() != out.channels()) return r;
    for (int k = 0; k < 3; ++k)
        if (! core::exactlyEqual (in.params().bandEdgeHz[k], out.params().bandEdgeHz[k])) return r;
    const long long n = std::min (in.blockCount(), out.blockCount());
    r.blocks = n;
    if (n <= 0) return r;

    // THE COARSE ALIGNMENT CHECK — where the two full-band peak series agree best, in blocks. See the note on
    // `lagBlocks`: this catches a gross mistake and cannot catch a small one.
    {
        // NO EVIDENCE IS NOT A LAG. `best` began at -1 and every correlation of a silent side is 0, so the
        // FIRST lag tried won and a muted master read "misaligned by 8 blocks". The correlation at lag 0 is
        // the reference now: another lag has to BEAT it, so a flat or empty series answers 0.
        const long long span = std::min<long long> (8, n - 1);
        double best = 0.0; long long bestLag = 0;
        for (long long j = 0; j < n; ++j)
            best += in.blockPeakLin (j, BandCrest::kFull) * out.blockPeakLin (j, BandCrest::kFull);
        for (long long lag = -span; lag <= span; ++lag)
        {
            double acc = 0.0;
            for (long long j = 0; j < n; ++j)
            {
                const long long k = j + lag;
                if (k < 0 || k >= n) continue;
                acc += in.blockPeakLin (j, BandCrest::kFull) * out.blockPeakLin (k, BandCrest::kFull);
            }
            if (acc > best) { best = acc; bestLag = lag; }
        }
        r.lagBlocks = bestLag;
    }
    if (scratch.size() < (std::size_t) n) scratch.resize ((std::size_t) n);

    std::size_t k = 0;
    double sum = 0.0;
    std::vector<double> shiftPeak, shiftLevel;
    shiftPeak.reserve ((std::size_t) n); shiftLevel.reserve ((std::size_t) n);
    for (long long j = 0; j < n; ++j)
    {
        if (! in.blockActive (j, band)) continue;       // THE SOURCE'S MASK, and only the source's
        ++r.inActive;
        const double pA = in.blockPeakLin (j, band),  mA = in.blockMeanSq (j, band);
        const double pB = out.blockPeakLin (j, band), mB = out.blockMeanSq (j, band);
        // A MASTER BAND THAT IS DIGITAL SILENCE is counted rather than skipped in silence: the loss there is
        // infinite, and an infinity cannot enter a quantile — but a band the chain muted must not read as a
        // band it left alone, which is what dropping the block quietly would do.
        if (! (pA > 0.0) || ! (mA > 0.0)) continue;
        if (! (pB > 0.0) || ! (mB > 0.0)) { ++r.outSilent; continue; }
        const double ratio = (pA * pA * mB) / (pB * pB * mA);
        if (! (ratio > 0.0) || ! std::isfinite (ratio)) continue;
        const double loss = 10.0 * core::det::log10 (ratio);
        if (! std::isfinite (loss)) continue;
        scratch[k++] = loss;
        const double pos = loss > 0.0 ? loss : 0.0;
        sum += pos;
        if (pos > r.maxDb) r.maxDb = pos;
        if (loss > 1.0) ++r.over1Db;
        if (loss > 3.0) ++r.over3Db;
        if (loss > 6.0) ++r.over6Db;
        shiftPeak.push_back (20.0 * core::det::log10 (pB / pA));
        shiftLevel.push_back (10.0 * core::det::log10 (mB / mA));
    }
    r.usable = (long long) k;
    if (k == 0) return r;

    std::sort (scratch.begin(), scratch.begin() + (std::ptrdiff_t) k);
    // NEAREST RANK, on the SIGNED losses, then the positive part is taken — the two commute for a quantile
    // because max(0, .) is monotone, and doing it in this order keeps `p5Db` signed.
    auto at = [&] (double q)
    {
        auto i = (std::size_t) std::llround (q * (double) (k - 1));
        if (i >= k) i = k - 1;
        return scratch[i];
    };
    r.p50Db = std::max (0.0, at (0.50));
    r.p95Db = std::max (0.0, at (0.95));
    r.p5Db  = at (0.05);
    r.meanDb = sum / (double) k;
    // CVaR95 — the mean of the worst 5 %, and at least one entry, so a short programme answers its own maximum
    // rather than nothing.
    const std::size_t tail = std::max<std::size_t> (1, (std::size_t) ((double) k * 0.05 + 0.5));
    double t = 0.0;
    for (std::size_t i = k - tail; i < k; ++i) t += std::max (0.0, scratch[i]);
    r.cvar95Db = t / (double) tail;

    std::sort (shiftPeak.begin(), shiftPeak.end());
    std::sort (shiftLevel.begin(), shiftLevel.end());
    r.peakShiftDb  = shiftPeak[shiftPeak.size() / 2];
    r.levelShiftDb = shiftLevel[shiftLevel.size() / 2];
    r.valid = in.valid() && out.valid();
    return r;
}

} // namespace felitronics::analysis

// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/core/Math.h>
#include <felitronics/dynamics/GainComputer.h>
#include <felitronics/dynamics/GainReductionFollower.h>
#include <felitronics/dynamics/LinkedDetector.h>

#include <cmath>

namespace felitronics::dynamics
{

// The parameters of the whole gain-reduction path, detector included. Derives from DetectorParams for
// the reason that one derives from anything here: so that handing this object to a LinkedDetector is a
// base-class binding rather than a field mapping somebody has to keep in step. `CompressorParams`
// derives from THIS in turn and adds only what lives outside the path (makeup, lookahead).
struct GainReductionParams : DetectorParams
{
    Mode mode = Mode::DownCompress;

    double thresholdDb = -18.0;
    double ratio       = 2.0;
    double kneeDb      = 6.0;
    double rangeDb     = 60.0;

    double attackMs    = 10.0;                    // ballistics on the GAIN REDUCTION (not the level)
    double releaseMs   = 100.0;
};

//==============================================================================
// felitronics::dynamics::GainReductionPath — key frame → SIGNED gain delta in dB, as ONE object:
//
//     LinkedDetector (gate → link → envelope) → gainToDb → GainComputer (curve) → GR ballistics
//
// It exists for the same reason LinkedDetector does, one stage further along. P2 made "the offline
// analysis ran the same DETECTOR" a structural fact; the three lines downstream of it — the dB
// conversion, the stateless curve, and the cast to float that feeds the follower — were still written
// out inside `Compressor::process`, so anything that wanted to reproduce the compressor's gain
// reduction offline had to copy them, in the right order, with the right types, and keep copying every
// future fix. That is exactly the drift LinkedDetector was created to remove, stopped half way.
//
// WHAT IS AND IS NOT IN HERE. In: the detector, the curve, the ballistics. Out: makeup (a decision
// about output level, not about how much to compress), lookahead and the delay line (a decision about
// alignment), and applying the gain to anything at all. That split is what lets an offline pass drive
// this object over a key and get, sample for sample, the same gain reduction the running compressor
// computes — while the compressor keeps owning what it does with it.
//
// THE dB CONVERSION IS PART OF THE PATH, and deliberately so. `core::gainToDb` clamps its argument at
// 1e-12, i.e. digital silence enters the curve as -240 dB. That is a FLOOR, not a measurement, and it
// has a consequence worth stating where it can be read: a threshold BELOW -240 dB makes silence an
// ACTIVE sample. Measured, with `thresholdDb = -300`, ratio 4 and a 60 dB range: digital silence earns
// 45 dB of gain reduction. That is arithmetic, not a defect — but it is why an offline solver may not
// reason about the curve in the linear domain, where no such floor exists.
//
// EXACTNESS. Same object, same parameters, same sample rate, same flush cadence, same binary ⇒ the
// same bits. Across toolchains the desktop tier declares `-ffp-contract=on` and libm is not
// bit-portable, so a pass compiled elsewhere agrees to rounding, not to the last bit — the same
// honest reading LinkedDetector states for the envelope.
//
// RT-safe: no allocation, lock, IO or throw. State is the detector's one float plus the follower's.
class GainReductionPath
{
public:
    void prepare (double sampleRate) noexcept
    {
        fs_ = (std::isfinite (sampleRate) && sampleRate > 0.0) ? sampleRate : 48000.0;
        det_.prepare (fs_);
        grf_.prepare (fs_);
        apply();
        reset();
    }

    void setParams (const GainReductionParams& p) noexcept { p_ = p; apply(); }

    const GainReductionParams& params()     const noexcept { return p_; }
    double                     sampleRate() const noexcept { return fs_; }

    void reset() noexcept { det_.reset(); grf_.reset(); }

    // One frame of the KEY at `sampleIndex` → the smoothed, SIGNED gain delta in dB. Negative is a
    // reduction (DownCompress, DownExpand), positive a boost (UpCompress).
    // PRECONDITION: `key` is non-null and `numKeyChannels >= 1` — the same contract LinkedDetector
    // states, for the same reason: this is the primitive, and resolving "no key" into a source is the
    // composite's job.
    inline float process (const float* const* key, int numKeyChannels, int sampleIndex) noexcept
    {
        return step (det_.process (key, numKeyChannels, sampleIndex));
    }

    // Same, for a caller that already holds one mono key sample.
    inline float processSample (float keySample) noexcept { return step (det_.processSample (keySample)); }

    // The chain FROM THE DETECTOR'S OUTPUT ONWARD, for a caller that already has the linked level —
    // an offline search that recorded the envelope once and now varies the curve over it. It is the
    // same `step()` the two entry points above end in, which is the point: a cached envelope and a
    // live key are the same arithmetic, not two implementations that agree today.
    // NB it does NOT advance the detector; a caller mixing this with `process()` would be driving two
    // different histories.
    inline float processLevel (float linkedLevel) noexcept { return step (linkedLevel); }

    float valueDb()       const noexcept { return grf_.valueDb(); }   // last gain delta (signed dB)
    float detectorLevel() const noexcept { return det_.level(); }     // linked level, linear amplitude

    // LAW 11c — advance the WHOLE path through `n` samples of DIGITAL SILENCE. Bit for bit what
    // `processSample(0.0f)` repeated `n` times does, and it is that same call that runs in phase one, so
    // there is no second implementation of the arithmetic to drift.
    //
    // The linked level of an all-zero frame is exactly +0.0f at EVERY width, width zero included
    // (`detail::linkAmplitudeImpl` returns 0 for `numChannels <= 0`, `|get(0)|` for 1, a max or a
    // `sqrt(0/nc)` above that), so `processSample(0.0f)` is what a silent block does to this object no
    // matter how wide the caller's silence was. That is what makes the pause well posed at width zero.
    //
    // TWO PHASES, and the split is where the money is. While the detector level is still ABOVE
    // `core::kGainToDbFloor` the curve's input moves, so the honest call runs — a log10 and a knee per
    // sample. At or below the floor `gainToDb` clamps, so the curve's argument, and therefore its output,
    // is the SAME NUMBER for the rest of the pause: compute it once and hand it to the follower, which
    // still picks its own attack/release branch per sample. The level cannot climb back out on its own
    // (the silent recurrence is `env *= c`, non-increasing for `c` in [0,1]), so the phase boundary is
    // crossed once and never re-crossed. Measured at 48 kHz ON THE DETECTOR THIS PATH ACTUALLY RUNS, the
    // expensive phase is ~1.75x shorter than the settling horizon behind it: 13 046 samples against
    // 23 393 at a 5 ms Rms window, 255 745 against 448 304 at 100 ms, 2 452 365 against 4 265 082 at 1 s.
    // In PEAK mode phase 1 is ONE step, because `LinkedDetector` makes Peak instant (times 0, 0).
    // An earlier draft said "3.5-19x" here. Those digits are real but they belong to a Peak follower WITH
    // a release, which this path never has — its Peak is instant and its Rms is symmetric on the POWER,
    // which halves the crossing. The 19 was two horizons divided by each other and meant nothing at all.
    void advanceSilence (int n) noexcept
    {
        int i = 0;
        // THE PREDICATE IS SPELLED IN DOUBLE because `gainToDb`'s clamp is: it tests `gain > 1.0e-12`
        // on the float level PROMOTED to double, and this is the same test on the same value.
        // BE HONEST ABOUT WHAT THAT BUYS, because the first version of this comment overstated it: it
        // claimed the two spellings straddle a band of levels, and they do not. `(float) 1e-12` is
        // 9.99999996e-13, BELOW the clamp, so no float at all separates `(double) x > 1e-12` from
        // `x > (float) 1e-12` — the float spelling merely runs a few more honest iterations before it
        // stops, and is equally correct. The cast is kept because it makes the predicate the clamp's own
        // test rather than one that happens to agree; it is not load-bearing, and a mutation of it is an
        // equivalent mutant, which is how the stand classifies it.
        for (; i < n && (double) det_.level() > core::kGainToDbFloor; ++i)
        {
            const float e0 = det_.stateWord(), g0 = grf_.valueDb();
            (void) processSample (0.0f);
            if (core::sameBits (det_.stateWord(), e0) && core::sameBits (grf_.valueDb(), g0)) return;
            if (! std::isfinite (det_.stateWord()) && ! std::isfinite (e0)) break;   // let phase 2 finish it
        }
        if (i >= n) return;

        // The level is pinned under the floor: ONE evaluation of the same expression `step()` uses, then
        // the follower alone. Spelling it through `gc_`/`gainToDb` rather than through a hard-coded
        // -240 dB keeps the mode, threshold, ratio, knee and range in the answer — a threshold below the
        // floor makes digital silence an ACTIVE sample, which the header says and a constant would lose.
        const float target = (float) gc_.deltaDb (core::gainToDb (det_.level()));
        grf_.advanceConstant (target, n - i);
        det_.advanceSilence (n - i);
    }

    // The curve itself, for the callers that have to READ it — a static auto-makeup asks the curve
    // what a 0 dBFS signal would get, and an offline solver asks it where the knee starts.
    const GainComputer& curve() const noexcept { return gc_; }

    // Law 8: call once per block. Both members hold decaying recursive state.
    void flushDenormals() noexcept { det_.flushDenormals(); grf_.flushDenormals(); }

private:
    // The three lines this class exists to own. THE ORDER AND THE TYPES ARE THE CONTRACT: the level
    // enters the curve as a double (gainToDb returns one, and the curve's knee arithmetic is done in
    // double), and the resulting delta is narrowed to float exactly once, before the follower — which
    // carries float state. Reproducing this in a second place is what this class prevents.
    inline float step (float level) noexcept
    {
        return grf_.process ((float) gc_.deltaDb (core::gainToDb (level)));
    }

    void apply() noexcept
    {
        gc_.setMode (p_.mode);
        gc_.setThresholdDb (p_.thresholdDb);
        gc_.setRatio (p_.ratio);
        gc_.setKneeDb (p_.kneeDb);
        gc_.setRangeDb (p_.rangeDb);
        grf_.setTimes (p_.attackMs, p_.releaseMs);   // non-finite times → instant (the follower's guard)
        det_.setParams (p_);                          // SLICING, on purpose — see GainReductionParams
    }

    LinkedDetector        det_;
    GainComputer          gc_;
    GainReductionFollower grf_;
    GainReductionParams   p_ {};
    double                fs_ = 48000.0;
};

} // namespace felitronics::dynamics

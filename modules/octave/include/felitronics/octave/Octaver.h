// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/core/DetMath.h>
#include <felitronics/core/FlushToZero.h>
#include <felitronics/core/Math.h>
#include <felitronics/core/StateGrid.h>
#include <felitronics/dynamics/EnvelopeFollower.h>
#include <felitronics/eq/Svf.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>

namespace felitronics::octave
{

//==============================================================================
// felitronics::octave::Octaver — the bass amp's sub-octave: the analog octaver of the OC-2 family, not a pitch
// shifter. It does not analyse a window and resynthesise; it follows the string's zero crossings one by one and
// divides them, so nothing waits for anything and the dry is never delayed.
//
// THE TRACKER. The playing channels are summed (scaled by 1/nch, so a mono source split to two channels tracks
// exactly as it does in mono) and band-limited to the bass fundamental range: a 2-pole high pass at 20 Hz (low B
// is 30.9 Hz) and a 2-pole low pass at 400 Hz, then a one-pole at 50 Hz that tilts the band toward its lowest
// partial. A comparator with hysteresis squares that up — and its thresholds are 0.88 of each POLARITY's own
// decaying peak, not a fixed level: a sawtooth's filtered wave swings twice as far one way as the other, and one
// symmetric threshold sits where the short side never reaches, so the comparator sticks (measured: a 30.9 Hz saw
// made two edges in three seconds). Thresholds that close to the peaks are what keeps a harmonic-rich note from
// double-triggering: a partial makes a second bump inside a period, and only the period's tallest bump reaches 88 %
// of a peak that has decayed for half a period (40 ms release). A rising edge sooner than 1/480 s after the last
// one is ignored — no bass fundamental is that high — which is what keeps a slap's click out of the divider.
//
// THE DIVIDER. Each accepted rising edge flips a flip-flop: −1 octave. Every second flip flips another: −2. The
// squares change sign mid-sample where the crossing actually lay (the sample carries the box average of the
// step), so the sub's period does not wobble on the sample grid.
//
// THE VOICE. Each square is scaled by the band's peak over the last input cycle — let down only at an edge, so a
// steady note gives a steady sub with no follower ripple on it, and up at once when the cycle in progress is
// already louder, so an accent is not a cycle late — glided over 10 ms, and by a lock gain, then rounded by its own 2-pole
// low pass (TONE) and only THEN by its level: a level move rides a smooth signal, and a level of zero is an exact
// zero with nothing ringing behind it. At level 1 the sub's peak matches the band's.
//
// THE GATE. The sub plays only while the band is above the gate (a 1 ms / 30 ms peak follower, 6 dB of hysteresis,
// 50 ms hold — longer than the lowest note's period, so a follower's ripple on a decaying low B cannot chatter it)
// AND an edge has come within the last 1/22 s. Either one failing fades the sub out over 4 ms; an edge with the gate
// open fades it back in. Between notes the sub is an exact zero.
//
// MONOPHONIC. A chord has no one fundamental: the comparator follows whatever its sum crosses — usually the lowest
// note when it is the loudest in the band, otherwise the beat pattern between the notes — and the sub comes out as
// that, often a growl. A single note at a time is the instrument's use; a double stop with a low root mostly
// tracks the root.
//
// CONTRACT (DSP-ARCHITECTURE.md §2).
//   · prepare(sampleRate, maxBlock, maxChannels) sets every coefficient; process()/reset() never allocate, lock or
//     throw. maxChannels is 1 or 2 and the rate 8 kHz … 768 kHz — anything else is REFUSED, and a refused prepare
//     DISARMS the instance: process() refuses until a prepare succeeds (law 11b).
//   · process() is [[nodiscard]] bool (law 11): malformed, unprepared, or wider than prepared is refused before
//     anything moves; any n >= 0 is processed in full (no block-sized scratch, so maxBlock is not a limit).
//   · The sub is ONE signal added to every channel; the dry is scaled by its own level and is otherwise untouched,
//     sample for sample: latencySamples() == 0.
//   · A NARROWER call is legal: the tracker follows the channels still playing (law 11a — a stopped channel is no
//     longer in the sum; the sum's memory is shared and stays). nch == 0 with n > 0 means EVERY channel stopped:
//     the tracker and both voices hold nothing but the stopped audio, so they are dropped — the sub cannot ring
//     into the return — while the three level glides, the stage's only shared ballistics, travel on exactly as
//     through silence (law 11c).
//   · reset() is a stream RESTART: tracker and voices dropped, the glides re-snapped to their targets on the next
//     call, the audio-time grid re-anchored.
//   · AT REST NOTHING IS COMPUTED. When both sub levels stand at zero and the dry at one, the stage drops its tracker
//     at the next audio-time grid point (law 8a, every 64 samples) and from then on does not touch the buffer at
//     all: bit-exact dry. A fresh stage set that way is at rest from its first call. Turning a level up starts the
//     tracker cold; it locks on the next rising edge.
//   · Law 8: every recursive state (the filters, the followers, the per-cycle level) is flushed to zero below 1e-15
//     on the audio-time grid (core::StateGrid), so the output does not depend on how the caller slices the stream.
//   · POISON: the tracker never sees a non-finite sample (it reads such a sample as zero) or one beyond ±4 (+12 dBFS),
//     so its state stays finite by construction and a wild sample cannot make a wild sub; the dry path is not this stage's to sanitise — a NaN in comes out where
//     it went in. Knob values are bounded: a non-finite one is ignored, a finite one clamped.
//   · THREADS: prepare, reset, the setters and process belong to ONE thread.
//   · Math: coefficients through core::det (exp2, pow10, tan) — the same bits on every row; the gate's follower is
//     dynamics::EnvelopeFollower, whose coefficient is the system exp (manifest: retain-rt).
//==============================================================================

/** The knobs, in the units the stage takes. `Octaver::setParams` applies them through the same guarded setters a
    caller may use one by one. `paramsForAmount()` (OctaveLaw.h) is the one-slider face of these. */
struct OctaverParams
{
    float sub1   = 1.0f;     ///< −1 octave level, 0 … 2 (1: the sub's peak matches the band's)
    float sub2   = 0.0f;     ///< −2 octave level, 0 … 2
    float dry    = 1.0f;     ///< dry level, 0 … 1
    float toneHz = 120.0f;   ///< the sub's low pass, 30 … 2000 Hz — lower is rounder
    float gateDb = -50.0f;   ///< tracking gate, dBFS of the 20 … 400 Hz band, −90 … −20 — lower tracks quieter notes
};

namespace detail
{
    inline float clampf (float lo, float hi, float v) noexcept { return v < lo ? lo : (hi < v ? hi : v); }

    inline constexpr double kLog2E = 1.4426950408889634074;

    /** e^x, the same bits on every row. */
    inline double detExp (double x) noexcept { return core::det::exp2 (core::det::mul (x, kLog2E)); }

    /** A level glide that LANDS: a counted linear ramp from where the level stood to its target. The value is a
        function of the count, never an accumulation, so it is the same bits however the ramp is sliced, it can be
        advanced over a pause in one step, and at the end it is the target exactly — a level of zero is a zero. */
    struct Ramp
    {
        float from = 0.0f, to = 0.0f;
        int   pos = 1, len = 1;
        float inv = 1.0f;

        void setLength (int n) noexcept { len = std::max (1, n); inv = 1.0f / (float) len; pos = len; from = to; }
        void snap() noexcept { from = to; pos = len; }
        void setTarget (float v) noexcept
        {
            if (core::sameBits (v, to)) return;
            from = value();
            to   = v;
            pos  = 0;
        }
        float value() const noexcept { return pos >= len ? to : from + (to - from) * ((float) pos * inv); }
        float next() noexcept { if (pos < len) ++pos; return value(); }
        void  advance (int n) noexcept { pos = n >= len - pos ? len : pos + n; }
        bool  landedAt (float v) const noexcept { return pos >= len && core::sameBits (to, v); }
    };
} // namespace detail

class Octaver
{
public:
    Octaver() noexcept { setParams (OctaverParams {}); }

    static constexpr int    kMaxChannels   = 2;
    static constexpr double kMinSampleRate = 8000.0, kMaxSampleRate = 768000.0;

    static constexpr float kMaxSub = 2.0f, kMaxDry = 1.0f;
    static constexpr float kMinToneHz = 30.0f, kMaxToneHz = 2000.0f;
    static constexpr float kMinGateDb = -90.0f, kMaxGateDb = -20.0f;

    // The tracker, fixed by measurement (README.md): the band, the tilt, the hysteresis and its release.
    static constexpr double kBandLowHz  = 20.0;
    static constexpr double kBandHighHz = 400.0;
    static constexpr double kTiltHz     = 50.0;
    static constexpr float  kHysteresis = 0.88f;
    static constexpr double kThresholdReleaseSeconds = 0.040;
    static constexpr double kMinTrackHz = 22.0;    // an edge later than 1/22 s after the last: the note is lost
    static constexpr double kMaxTrackHz = 480.0;   // an edge sooner than 1/480 s after the last: not a fundamental
    static constexpr float  kInputLimit = 4.0f;    // what the tracker reads: +12 dBFS; beyond is clamped

    // The voice and the gate.
    static constexpr double kLevelGlideSeconds  = 0.010;  // the per-cycle level's glide between cycles
    static constexpr double kLockRampSeconds    = 0.004;  // the sub's fade in/out on lock and gate
    static constexpr double kGateAttackMs       = 1.0, kGateReleaseMs = 30.0;
    static constexpr double kGateHoldSeconds    = 0.050;
    static constexpr float  kGateCloseRatio     = 0.5f;   // the gate closes 6 dB under where it opens
    static constexpr float  kFloorUnderGate     = 0.01f;  // the comparator's floor: 40 dB under the gate
    static constexpr double kLevelRampSeconds   = 0.020;  // sub1 / sub2 / dry glides

    /** Sets every coefficient for `maxChannels` (1 or 2) at `sampleRate` (8 kHz … 768 kHz). `maxBlock` (>= 1) is
        recorded, not a limit. A refused prepare DISARMS the instance. */
    [[nodiscard]] bool prepare (double sampleRate, int maxBlock, int maxChannels) noexcept
    {
        prepared_ = false;   // law 11b: whatever happens next, the old preparation no longer stands

        if (! (sampleRate >= kMinSampleRate && sampleRate <= kMaxSampleRate) || maxBlock < 1
            || maxChannels < 1 || maxChannels > kMaxChannels)
            return false;

        fs_       = sampleRate;
        maxBlock_ = maxBlock;
        maxCh_    = maxChannels;

        hpBand_.prepare (fs_, 1);
        lpBand_.prepare (fs_, 1);
        hpBand_.setParams (eq::FilterType::HighPass, kBandLowHz,  kButterworthQ, 0.0);
        lpBand_.setParams (eq::FilterType::LowPass,  kBandHighHz, kButterworthQ, 0.0);
        tone1_.prepare (fs_, 1);
        tone2_.prepare (fs_, 1);
        applyTone();

        gateEnv_.prepare (fs_);
        gateEnv_.setTimes (kGateAttackMs, kGateReleaseMs);
        gateEnv_.setDetector (dynamics::Detector::Peak);
        applyGate();

        tiltCoeff_ = (float) (1.0 - detail::detExp (-2.0 * core::kPi * kTiltHz / fs_));
        peakRel_   = (float) detail::detExp (-1.0 / (kThresholdReleaseSeconds * fs_));
        ampCoeff_  = (float) (1.0 - detail::detExp (-1.0 / (kLevelGlideSeconds * fs_)));
        minPeriod_ = std::max (1, (int) std::floor (fs_ / kMaxTrackHz));
        maxPeriod_ = (int) std::ceil (fs_ / kMinTrackHz);
        lockLen_   = std::max (1, (int) std::lround (kLockRampSeconds * fs_));
        lockStep_  = 1.0f / (float) lockLen_;
        gateHoldLen_ = std::max (1, (int) std::lround (kGateHoldSeconds * fs_));

        const int rampLen = std::max (1, (int) std::lround (kLevelRampSeconds * fs_));
        sub1_.setLength (rampLen);
        sub2_.setLength (rampLen);
        dry_.setLength (rampLen);

        cold_ = false;   // force the drop: a re-prepare holds nothing it heard before
        dropTracker();
        grid_.reset();
        primed_   = false;
        prepared_ = true;
        return true;
    }

    bool isPrepared()  const noexcept { return prepared_; }
    int  maxChannels() const noexcept { return prepared_ ? maxCh_ : 0; }
    int  maxBlock()    const noexcept { return prepared_ ? maxBlock_ : 0; }
    static constexpr int latencySamples() noexcept { return 0; }

    /** A stream restart: the tracker and the voices dropped, the glides re-snapped on the next call, the grid
        re-anchored. Cheap and idempotent. */
    void reset() noexcept
    {
        dropTracker();
        grid_.reset();
        primed_ = false;
    }

    void setParams (const OctaverParams& p) noexcept
    {
        setSub1 (p.sub1);
        setSub2 (p.sub2);
        setDry (p.dry);
        setToneHz (p.toneHz);
        setGateDb (p.gateDb);
    }

    /** −1 octave level, 0 … 2; glides over kLevelRampSeconds. */
    void setSub1 (float level) noexcept { if (std::isfinite (level)) sub1_.setTarget (detail::clampf (0.0f, kMaxSub, level)); }
    /** −2 octave level, 0 … 2; glides over kLevelRampSeconds. */
    void setSub2 (float level) noexcept { if (std::isfinite (level)) sub2_.setTarget (detail::clampf (0.0f, kMaxSub, level)); }
    /** Dry level, 0 … 1; glides over kLevelRampSeconds. */
    void setDry (float level) noexcept { if (std::isfinite (level)) dry_.setTarget (detail::clampf (0.0f, kMaxDry, level)); }

    /** The sub's low pass, 30 … 2000 Hz. Both voices' filters take it at once (a TPT filter does not click). */
    void setToneHz (float hz) noexcept
    {
        if (! std::isfinite (hz))
            return;
        const float want = detail::clampf (kMinToneHz, kMaxToneHz, hz);
        if (core::sameBits (want, toneHz_))
            return;
        toneHz_ = want;
        applyTone();
    }

    /** The tracking gate, dBFS of the 20 … 400 Hz band, −90 … −20. */
    void setGateDb (float db) noexcept
    {
        if (! std::isfinite (db))
            return;
        const float want = detail::clampf (kMinGateDb, kMaxGateDb, db);
        if (core::sameBits (want, gateDb_))
            return;
        gateDb_ = want;
        applyGate();
    }

    /** In place. Channels at index >= nch are not touched (they stopped — law 11a). */
    [[nodiscard]] bool process (float* const* io, int numChannels, int numSamples) noexcept
    {
        if (numChannels < 0 || numSamples < 0 || ! prepared_ || numChannels > maxCh_)
            return false;   // law 11: refused before anything moves
        if (numSamples == 0)
            return true;    // the only true no-op: no time, no edge
        if (numChannels > 0 && io == nullptr)
            return false;

        if (! primed_)
        {
            // A fresh or restarted stage starts ON its levels rather than fading in from wherever they stood.
            sub1_.snap(); sub2_.snap(); dry_.snap();
            primed_ = true;
        }

        if (numChannels == 0)
        {
            // Law 11a/11c — every channel stopped: what the tracker and the voices hold is the stopped audio, so it
            // goes; the level glides are shared ballistics and travel on as through silence.
            dropTracker();
            sub1_.advance (numSamples); sub2_.advance (numSamples); dry_.advance (numSamples);
            grid_.skip (numSamples);
            return true;
        }

        const float invN = 1.0f / (float) numChannels;   // exact for 1 and 2

        for (int done = 0; done < numSamples; )
        {
            const int seg = grid_.segment (numSamples - done);

            if (! (cold_ && atRestLevels()))
            {
                cold_ = false;
                runSegment (io, numChannels, done, seg, invN);
            }
            // else: AT REST — the dry stands at one and no sub plays; the buffer is not touched.

            done += seg;
            if (grid_.advance (seg))
                maintain();
        }
        return true;
    }

    /** The glided levels where they stand now — what the next sample would use before its own step. */
    float sub1Level() const noexcept { return sub1_.value(); }
    float sub2Level() const noexcept { return sub2_.value(); }
    float dryLevel()  const noexcept { return dry_.value(); }

    /** True while the sub plays: the gate is open and an edge came within the last 1/22 s. */
    bool isLocked() const noexcept { return gateOpen_ && sinceEdge_ <= maxPeriod_; }
    /** True while the stage is at rest and computes nothing (see the contract). */
    bool isAtRest() const noexcept { return cold_ && atRestLevels(); }
    /** Accepted rising edges since the tracker last started cold — a tap for tests and a tuner light. */
    std::int64_t edgeCount() const noexcept { return edges_; }

    /** True when any of the stage's own recursive words sits in the subnormal band — the law-8 witness for the
        tests. (The four eq::Svf filters keep their words private; eq's own suite and LAW8-AUDIT.md answer for
        them, and this stage flushes them on the same grid.) */
    bool anySubnormalState() const noexcept
    {
        const auto subn = [] (float v) { return std::fpclassify (v) == FP_SUBNORMAL; };
        return subn (tilt_) || subn (prevT_) || subn (peakPos_) || subn (peakNeg_) || subn (cycleMax_)
            || subn (amp_) || subn (ampTarget_) || subn (gateEnv_.stateWord());
    }

private:
    static constexpr double kButterworthQ = 0.70710678118654752;

    bool atRestLevels() const noexcept
    {
        return sub1_.landedAt (0.0f) && sub2_.landedAt (0.0f) && dry_.landedAt (1.0f);
    }

    void applyTone() noexcept
    {
        tone1_.setParams (eq::FilterType::LowPass, (double) toneHz_, kButterworthQ, 0.0);
        tone2_.setParams (eq::FilterType::LowPass, (double) toneHz_, kButterworthQ, 0.0);
    }

    void applyGate() noexcept
    {
        gateOn_  = (float) core::det::pow10 ((double) gateDb_ / 20.0);
        gateOff_ = gateOn_ * kGateCloseRatio;
        floor_   = gateOn_ * kFloorUnderGate;
    }

    /** The tracker and both voices back to a stage that has heard nothing: the stage is COLD. */
    void dropTracker() noexcept
    {
        if (cold_)
            return;   // idempotent: a pause calls this every block, and cold is cold

        hpBand_.reset(); lpBand_.reset(); tone1_.reset(); tone2_.reset();
        gateEnv_.reset();
        tilt_ = prevT_ = peakPos_ = peakNeg_ = cycleMax_ = amp_ = ampTarget_ = 0.0f;
        high_ = false;
        q1_ = q2_ = -1.0f;
        sinceEdge_ = maxPeriod_ + 1;
        lockPos_ = 0;
        gateOpen_ = false;
        gateHold_ = 0;
        edges_ = 0;
        cold_  = true;
    }

    /** Once per grid period of audio time: law 8's flush, and — if the levels stand at rest — going cold. */
    void maintain() noexcept
    {
        if (cold_)
            return;   // nothing ran since the stage went cold: there is nothing to flush

        hpBand_.flushDenormals(); lpBand_.flushDenormals(); tone1_.flushDenormals(); tone2_.flushDenormals();
        gateEnv_.flushDenormals();
        core::flushPoison (tilt_);    core::flushPoison (prevT_);
        core::flushPoison (peakPos_); core::flushPoison (peakNeg_);
        core::flushPoison (cycleMax_);
        core::flushPoison (amp_);     core::flushPoison (ampTarget_);

        if (atRestLevels())
            dropTracker();   // both sub levels are an exact zero: nothing the tracker holds can be heard again
    }

    void runSegment (float* const* io, int nch, int offset, int len, float invN) noexcept
    {
        for (int k = offset; k < offset + len; ++k)
        {
            float sum = 0.0f;
            for (int ch = 0; ch < nch; ++ch)
            {
                const float v = io[ch][k];
                sum += std::isfinite (v) ? detail::clampf (-kInputLimit, kInputLimit, v) : 0.0f;
            }

            const float sub = tick (sum * invN);
            const float d   = dry_.next();

            for (int ch = 0; ch < nch; ++ch)
                io[ch][k] = d * io[ch][k] + sub;
        }
    }

    /** One sample of the tracker and the voices; returns the sub. */
    float tick (float x) noexcept
    {
        // The bass fundamental range.
        const float b = lpBand_.processSample (0, hpBand_.processSample (0, x));

        // The gate: open above gateOn_, held open while above gateOff_ and for gateHoldLen_ after.
        const float g = gateEnv_.process (b);
        if (g >= gateOn_)
        {
            gateOpen_ = true;
            gateHold_ = gateHoldLen_;
        }
        else if (gateOpen_)
        {
            if (g >= gateOff_)        gateHold_ = gateHoldLen_;
            else if (gateHold_ > 0)   --gateHold_;
            else                      gateOpen_ = false;
        }

        // The comparator's input: the band tilted toward its lowest partial.
        tilt_ += tiltCoeff_ * (b - tilt_);
        const float t = tilt_;

        // Each polarity's own peak: instant up, kThresholdReleaseSeconds down.
        peakPos_ = t  > peakPos_ ? t  : peakPos_ * peakRel_;
        peakNeg_ = -t > peakNeg_ ? -t : peakNeg_ * peakRel_;
        const float hi = std::max (kHysteresis * peakPos_, floor_);
        const float lo = std::max (kHysteresis * peakNeg_, floor_);

        const float ab = std::fabs (b);
        cycleMax_ = ab > cycleMax_ ? ab : cycleMax_;

        if (sinceEdge_ <= maxPeriod_)
        {
            ++sinceEdge_;
            if (sinceEdge_ > maxPeriod_)
                cycleMax_ = ab;   // the note is lost: the next note's level starts from here
        }

        bool  edge = false;
        float frac = 1.0f;
        if (! high_)
        {
            if (t > hi)
            {
                high_ = true;
                if (sinceEdge_ >= minPeriod_)
                {
                    edge = true;
                    const float span = t - prevT_;   // > 0: t crossed hi from at or below it
                    frac = span > 0.0f ? detail::clampf (0.0f, 1.0f, (t - hi) / span) : 1.0f;
                    sinceEdge_ = 0;
                }
            }
        }
        else if (t < -lo)
        {
            high_ = false;
        }
        prevT_ = t;

        // The divider; at the edge sample the square carries the box average of its step.
        float s1 = q1_, s2 = q2_;
        if (edge)
        {
            ++edges_;
            q1_ = -q1_;
            const float box = 2.0f * frac - 1.0f;
            s1 = q1_ * box;
            if (q1_ > 0.0f)
            {
                q2_ = -q2_;
                s2  = q2_ * box;
            }
            else
            {
                s2 = q2_;
            }

            ampTarget_ = cycleMax_;   // the band's peak over the cycle that just ended
            cycleMax_  = 0.0f;
        }
        if (sinceEdge_ > maxPeriod_)
            ampTarget_ = 0.0f;
        else if (cycleMax_ > ampTarget_)
            ampTarget_ = cycleMax_;   // the cycle in progress is louder already: up now, not at its end
        amp_ += ampCoeff_ * (ampTarget_ - amp_);

        const bool locked = gateOpen_ && sinceEdge_ <= maxPeriod_;
        if (locked) { if (lockPos_ < lockLen_) ++lockPos_; }
        else        { if (lockPos_ > 0)        --lockPos_; }
        const float lock = lockPos_ == lockLen_ ? 1.0f : (float) lockPos_ * lockStep_;

        const float a  = amp_ * lock;
        const float y1 = tone1_.processSample (0, s1 * a);
        const float y2 = tone2_.processSample (0, s2 * a);
        return sub1_.next() * y1 + sub2_.next() * y2;
    }

    // The tracker's filters, the voices' tone and the gate's follower (core primitives, deterministic coefficients).
    eq::DeterministicSvf hpBand_, lpBand_, tone1_, tone2_;
    dynamics::EnvelopeFollower gateEnv_;

    // Knobs.
    detail::Ramp sub1_, sub2_, dry_;
    float toneHz_ = 120.0f, gateDb_ = -50.0f;
    float gateOn_ = 0.0f, gateOff_ = 0.0f, floor_ = 0.0f;

    // Coefficients and counts at the prepared rate.
    double fs_ = 48000.0;
    float  tiltCoeff_ = 0.0f, peakRel_ = 0.0f, ampCoeff_ = 0.0f, lockStep_ = 1.0f;
    int    minPeriod_ = 1, maxPeriod_ = 1, lockLen_ = 1, gateHoldLen_ = 1;

    // The tracker.
    float tilt_ = 0.0f, prevT_ = 0.0f, peakPos_ = 0.0f, peakNeg_ = 0.0f;
    float cycleMax_ = 0.0f, amp_ = 0.0f, ampTarget_ = 0.0f;
    float q1_ = -1.0f, q2_ = -1.0f;
    bool  high_ = false, gateOpen_ = false;
    int   sinceEdge_ = 2, lockPos_ = 0, gateHold_ = 0;
    std::int64_t edges_ = 0;

    core::StateGrid grid_;
    bool cold_ = false, primed_ = false, prepared_ = false;
    int  maxBlock_ = 0, maxCh_ = 0;
};

} // namespace felitronics::octave

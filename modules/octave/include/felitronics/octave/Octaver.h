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
#include <atomic>
#include <cmath>
#include <cstdint>
#include <initializer_list>

namespace felitronics::octave
{

//==============================================================================
// felitronics::octave::Octaver — the bass amp's sub-octave: the analog octaver of the OC-2 family, not a pitch
// shifter. It does not analyse a window and resynthesise; it follows the string's zero crossings one by one and
// divides them, so nothing waits for anything and the dry is never touched.
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
// THE NOTE'S LIFE. An edge COUNTS only with the gate open. A lost note is ACQUIRED again only on crossings that are
// fresh twice over: the band's peak over the last half to whole longest period stands at the gate's opening level
// (not the follower at the instant of the crossing — its ripple kept a steady note at 0.84 … 1.28 times the level
// from ever acquiring, and its tail after a loud note let the hum under it in), and the raw input itself has just
// swung that way — past a quarter of its own peak since the loss, DC taken off at 5 Hz — since the last counted
// crossing of the other direction (filter ringing after a source has stopped is not the source swinging). What
// plays under the gate, a DI's hum included, can neither prove a note nor set its period; the comparator's floor
// rides 20 dB under the band's own level, so hum under a sounding note never squares up at all.
//
// The note's PERIOD is the longest of the last three intervals between counted edges, and it sets the deadlines
// only once CONFIRMED — two intervals in a row within a quarter of each other; a note's first interval can be a
// partial attack cycle (813 samples on a low B whose period is 1555), and a deadline drawn from it lost the note
// before its next true edge. A confirmed note's next edge comes within about 1 % of its period; one a quarter period
// LATE means the string is dying, and from that sample the sub fades out over 4 ms — AFTER the filters, so it ends
// in an exact zero and the high pass is not handed a step to ring on — and fades back in on the next edge. A muted
// note therefore does not hold its last cycle as a DC plateau. Two periods without an edge (at most 1/22 s; 1/22 s
// while unconfirmed) and the note is LOST: the thresholds drop to nothing (a quieter next note is not locked out by
// a louder one's), and the next note must prove itself with TWO counted rises AND TWO counted falls before the sub
// plays again. One cycle of anything — the smoothest click is a single cycle with ramps — has at most two of one and
// one of the other, so no click makes a sub (measured on 7188cfc, which asked for two rises: 0.19 from a smooth one). Dying moves the fade only, never the thresholds: re-anchoring them on a
// late edge was tried, and on a weak fundamental it fed itself (one double trigger halved the period, the next true
// edge looked late, the lowered thresholds took the 2nd harmonic's bump). Nor is it a release of the level BEFORE
// the filters: that hands the sub's high pass a falling step, which comes out as one more pulse.
//
// THE REOPEN. While the sub is faded out its voices hold nothing (reset once, then idle) and its level is SET to the
// cycle in progress, not glided. When the fade reopens, each voice is primed where a long half-cycle of the other
// polarity would have left it, so its first output is a full step at the level of the note it reopens on — on an edge
// the edge itself; for a note proven on its second fall, a step there (its first half-cycle a little short) rather
// than half a cycle of near-silence waiting for the next edge. Measured on 7188cfc: a loud low E then a G2 20 dB
// quieter reopened at 6.29 times the G2's own level, from the level and the filter state the low E had left.
//
// THE DIVIDER. Each counted rising edge flips a flip-flop: −1 octave. Every second flip flips another: −2. The
// squares change sign mid-sample where the crossing actually lay (the sample carries the box average of the
// step), so the sub's period does not wobble on the sample grid.
//
// THE VOICE. Each square is scaled by the band's peak over the last input cycle — let down only at an edge, so a
// steady note gives a steady sub with no follower ripple on it, and up at once when the cycle in progress is already
// louder — glided over 10 ms; then rounded by its own 2-pole low pass (TONE), cleared of what no amp or speaker
// should be asked for by a 2-pole high pass (SUB HP, 25 Hz: the −1 of the B string's lowest frets is
// 15.4 … 20.6 Hz), and only THEN set by its level and the lock gain: a level of zero is an exact zero with nothing
// ringing behind it. At level 1 the sub's square swings to the band's peak before the filters. TONE and SUB HP glide
// to a new corner in log frequency over 20 ms, re-designed EVERY SAMPLE while they travel: swapping a filter's
// coefficients under live state in one step is a click (1075 times the held filter's curvature for 120 → 2000 Hz in
// mid-edge, measured on 3971592), and swapping them every 16 samples was a staircase of small ones (7.25 % per step).
//
// THE GATE. The sub plays only while the band is above the gate (a 1 ms / 30 ms peak follower, 6 dB of hysteresis,
// 50 ms hold — longer than the lowest note's period, so a follower's ripple on a decaying low B cannot chatter it)
// AND the note is locked and not dying. Between notes the sub is an exact zero. A gate down at −90 dBFS opens on a
// DI's own noise floor and the tracker follows the hiss — the range stops at −90 for that reason; −50 is the default.
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
//   · THE DRY IS NEVER TOUCHED: each output sample is its input sample plus the sub, the same sub on every channel
//     (latencySamples() == 0). There is no dry level — the amp has its own mix.
//   · A NARROWER call is legal: the tracker follows the channels still playing (law 11a — a stopped channel is no
//     longer in the sum). The tracker and the voices are ONE shared, linked instrument, and nch == 0 with n > 0 is a
//     pause that runs them exactly as n samples of silence would (law 11c) — the gate holds, glides travel, a note
//     dies — and stops computing the moment a silent sample would change no state word at all (the silent fixed
//     point), so a long pause costs its first second and a half at most.
//   · reset() is a stream RESTART: tracker and voices cleared, the levels and the corners re-snapped to their
//     targets on the next call, the audio-time grid re-anchored.
//   · IDLE NOTHING IS COMPUTED. When both sub levels stand at exactly zero the stage drops its tracker at the next
//     audio-time grid point (law 8a, every 64 samples) and from then on does not touch the buffer at all: the dry
//     bit for bit, isIdle() true. A fresh stage set that way is idle from its first call.
//   · A WAKE FROM IDLE IS A RESTART: the grid is re-anchored, the corners land on their targets, the levels glide up
//     from zero, and the tracker starts from silence through a 10 ms input door — it locks on the note's second
//     edge.
//   · Law 8: every recursive state is flushed to zero below 1e-15 on the audio-time grid (core::StateGrid), so the
//     output does not depend on how the caller slices the stream.
//   · POISON: the tracker never sees a non-finite sample (it reads such a sample as zero) or one beyond ±4 (+12 dBFS),
//     so its state stays finite by construction and a wild sample cannot make a wild sub; the dry path is not this
//     stage's to sanitise — a NaN in comes out where it went in. Knob values are bounded: a non-finite one is
//     ignored, a finite one clamped.
//   · THREADS: prepare, reset, the setters, process and the taps (sub1Level, sub2Level, the corners, edgeCount,
//     isIdle) belong to ONE thread. isLocked() — a tuner light — is a relaxed atomic any thread may read.
//   · Math: coefficients through core::det (exp2, log2, pow10, tan) — the same bits on every row; the gate's
//     follower is dynamics::EnvelopeFollower, whose coefficient is the system exp (manifest: retain-rt).
//==============================================================================

/** The knobs, in the units the stage takes. `Octaver::setParams` applies them through the same guarded setters a
    caller may use one by one. `paramsForAmount()` (OctaveLaw.h) is the one-slider face of these: it OWNS `sub1`
    and `sub2` (the levels now) and writes them from the menu's `sub1Top` / `sub2Top`, which it never changes — so the
    stage's current params can be handed back to it as the menu. */
struct OctaverParams
{
    float sub1Top = 1.25f;    ///< MENU: the −1 level at the slider's top (read by the law; the stage ignores it)
    float sub2Top = 0.0f;     ///< MENU: the −2 level at the slider's top (read by the law; the stage ignores it)
    float sub1    = 1.25f;    ///< −1 octave level now, 0 … 2 (1: the square swings to the band's peak) — the law's
    float sub2    = 0.0f;     ///< −2 octave level now, 0 … 2                                          — the law's
    float toneHz  = 120.0f;   ///< the sub's low pass, 30 … 2000 Hz — lower is rounder
    float subHpHz = 25.0f;    ///< the sub's high pass, 10 … 80 Hz — what no amp or speaker should be asked for
    float gateDb  = -50.0f;   ///< tracking gate, dBFS of the 20 … 400 Hz band, −90 … −20 — lower tracks quieter notes
};

namespace detail
{
    inline float clampf (float lo, float hi, float v) noexcept { return v < lo ? lo : (hi < v ? hi : v); }

    inline constexpr double kLog2E = 1.4426950408889634074;

    /** e^x, the same bits on every row. */
    inline double detExp (double x) noexcept { return core::det::exp2 (core::det::mul (x, kLog2E)); }

    /** A glide that LANDS: a counted linear ramp from where the value stood to its target. The value is a function
        of the count, never an accumulation, so it is the same bits however the ramp is sliced, it can be advanced
        over a stretch in one step, and at the end it is the target exactly — a level of zero is a zero. */
    struct Ramp
    {
        float from = 0.0f, to = 0.0f;
        int   pos = 1, len = 1;
        float inv = 1.0f;

        void setLength (int n) noexcept { len = std::max (1, n); inv = 1.0f / (float) len; pos = len; from = to; }
        void snap() noexcept { from = to; pos = len; }
        void snapTo (float v) noexcept { from = to = v; pos = len; }
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
        bool  landed() const noexcept { return pos >= len; }
        bool  landedAt (float v) const noexcept { return pos >= len && core::sameBits (to, v); }
    };

    /** +0.0f exactly (not −0, not a hair): the value every word of the silent fixed point holds. */
    inline bool plusZero (float v) noexcept { return core::sameBits (v, 0.0f); }
} // namespace detail

class Octaver
{
public:
    static constexpr int    kMaxChannels   = 2;
    static constexpr double kMinSampleRate = 8000.0, kMaxSampleRate = 768000.0;

    static constexpr float kMaxSub = 2.0f;
    static constexpr float kMinToneHz = 30.0f, kMaxToneHz = 2000.0f;
    static constexpr float kMinSubHpHz = 10.0f, kMaxSubHpHz = 80.0f;
    static constexpr float kMinGateDb = -90.0f, kMaxGateDb = -20.0f;

    // The tracker, fixed by measurement (README.md): the band, the tilt, the hysteresis and its release.
    static constexpr double kBandLowHz  = 20.0;
    static constexpr double kBandHighHz = 400.0;
    static constexpr double kTiltHz     = 50.0;
    static constexpr float  kHysteresis = 0.88f;
    static constexpr double kThresholdReleaseSeconds = 0.040;
    static constexpr double kMinTrackHz = 22.0;    // 1/22 s between edges: lost, whatever the period was
    static constexpr double kMaxTrackHz = 480.0;   // an edge sooner than 1/480 s after the last: not a fundamental
    static constexpr float  kInputLimit = 4.0f;    // what the tracker reads: +12 dBFS; beyond is clamped
    static constexpr float  kFloorUnderBand = 0.1f; // the comparator's floor rides 20 dB under the band's own level
    static constexpr double kSourceHpHz     = 5.0;  // the raw input's own swing is read with its DC taken off here
    static constexpr float  kSourceSwing    = 0.25f; // ...and must pass a quarter of its peak since the loss

    // The note's life, in its own periods.
    static constexpr int kDyingQuarters = 5;       // an edge 5/4 of a period late: the note is dying, the sub fades
    static constexpr int kLostPeriods   = 2;       // two periods without one: the note is lost (capped at 1/22 s)

    // The voice and the gate.
    static constexpr double kLevelGlideSeconds   = 0.010;  // the per-cycle level's glide between cycles
    static constexpr double kLostReleaseSeconds  = 0.005;  // ...and its release once the note is lost (behind the fade)
    static constexpr double kLockRampSeconds     = 0.004;  // the sub's fade in/out on lock and gate
    static constexpr double kGateAttackMs        = 1.0, kGateReleaseMs = 30.0;
    static constexpr double kGateHoldSeconds     = 0.050;
    static constexpr float  kGateCloseRatio      = 0.5f;   // the gate closes 6 dB under where it opens
    static constexpr float  kFloorUnderGate      = 0.01f;  // the comparator's floor: 40 dB under the gate
    static constexpr double kLevelRampSeconds    = 0.020;  // sub1 / sub2 glides
    static constexpr double kDoorSeconds         = 0.010;  // the tracker's input door on a wake
    static constexpr double kCornerGlideSeconds  = 0.020;  // TONE / SUB HP glides, in log frequency, re-designed per sample

    Octaver() noexcept { setParams (OctaverParams {}); }

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

        for (auto* f : { &hpBand_, &lpBand_, &tone1_, &tone2_, &hp1_, &hp2_ })
            f->prepare (fs_, 1);
        hpBand_.setParams (eq::FilterType::HighPass, kBandLowHz,  kButterworthQ, 0.0);
        lpBand_.setParams (eq::FilterType::LowPass,  kBandHighHz, kButterworthQ, 0.0);

        gateEnv_.prepare (fs_);
        gateEnv_.setTimes (kGateAttackMs, kGateReleaseMs);
        gateEnv_.setDetector (dynamics::Detector::Peak);
        applyGate();

        tiltCoeff_  = (float) (1.0 - detail::detExp (-2.0 * core::kPi * kTiltHz / fs_));
        peakRel_    = (float) detail::detExp (-1.0 / (kThresholdReleaseSeconds * fs_));
        ampCoeff_   = (float) (1.0 - detail::detExp (-1.0 / (kLevelGlideSeconds * fs_)));
        lostCoeff_  = (float) (1.0 - detail::detExp (-1.0 / (kLostReleaseSeconds * fs_)));
        sourceCoeff_ = (float) (1.0 - detail::detExp (-2.0 * core::kPi * kSourceHpHz / fs_));
        minPeriod_  = std::max (1, (int) std::floor (fs_ / kMaxTrackHz));
        maxPeriod_  = (int) std::ceil (fs_ / kMinTrackHz);
        freshLen_   = std::max (1, (maxPeriod_ + 1) / 2);
        lockLen_    = std::max (1, (int) std::lround (kLockRampSeconds * fs_));
        lockStep_   = 1.0f / (float) lockLen_;
        gateHoldLen_ = std::max (1, (int) std::lround (kGateHoldSeconds * fs_));

        sub1_.setLength ((int) std::lround (kLevelRampSeconds * fs_));
        sub2_.setLength ((int) std::lround (kLevelRampSeconds * fs_));
        door_.setLength ((int) std::lround (kDoorSeconds * fs_));
        toneLog_.setLength ((int) std::lround (kCornerGlideSeconds * fs_));
        hpLog_.setLength ((int) std::lround (kCornerGlideSeconds * fs_));

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

    /** A stream restart: the tracker and the voices cleared, the levels and corners re-snapped on the next call,
        the grid re-anchored. Cheap and idempotent. */
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
        setToneHz (p.toneHz);
        setSubHpHz (p.subHpHz);
        setGateDb (p.gateDb);
    }

    /** −1 octave level, 0 … 2; glides over kLevelRampSeconds. Non-finite: ignored. */
    void setSub1 (float level) noexcept { if (std::isfinite (level)) sub1_.setTarget (detail::clampf (0.0f, kMaxSub, level) + 0.0f); }
    /** −2 octave level, 0 … 2; glides over kLevelRampSeconds. Non-finite: ignored. */
    void setSub2 (float level) noexcept { if (std::isfinite (level)) sub2_.setTarget (detail::clampf (0.0f, kMaxSub, level) + 0.0f); }

    /** The sub's low pass, 30 … 2000 Hz; glides there in log frequency over kCornerGlideSeconds. */
    void setToneHz (float hz) noexcept
    {
        if (std::isfinite (hz))
            toneLog_.setTarget ((float) core::det::log2 ((double) detail::clampf (kMinToneHz, kMaxToneHz, hz)));
    }

    /** The sub's high pass, 10 … 80 Hz; glides there in log frequency over kCornerGlideSeconds. */
    void setSubHpHz (float hz) noexcept
    {
        if (std::isfinite (hz))
            hpLog_.setTarget ((float) core::det::log2 ((double) detail::clampf (kMinSubHpHz, kMaxSubHpHz, hz)));
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

    /** In place: each sample gets the sub added. Channels at index >= nch are not touched (law 11a). */
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
            // A fresh or restarted stage starts ON its levels and corners rather than gliding in from wherever
            // they stood.
            sub1_.snap(); sub2_.snap(); toneLog_.snap(); hpLog_.snap();
            applyCorners();
            primed_ = true;
        }

        const float invN = numChannels > 0 ? 1.0f / (float) numChannels : 0.0f;   // exact for 1 and 2

        for (int done = 0; done < numSamples; )
        {
            if (cold_)
            {
                if (atRestLevels())
                    break;   // IDLE: the buffer is not touched and nothing is computed
                wake();
            }

            const int seg = grid_.segment (numSamples - done);   // segments end on the 64-sample grid of audio time

            if (numChannels == 0 && trackerAtRest())
            {
                // The silent fixed point: a silent sample would change no tracker word, so only the counted
                // things move — exactly as far as `seg` silent samples would move them. Once they have landed too
                // (and the levels are not on their way to idle, which the next grid point must see), nothing moves
                // but the grid's phase, for the rest of the pause.
                if (sub1_.landed() && sub2_.landed() && door_.landed() && toneLog_.landed() && hpLog_.landed()
                    && ! atRestLevels())
                {
                    grid_.skip (numSamples - done);
                    break;
                }
                sub1_.advance (seg); sub2_.advance (seg); door_.advance (seg);
                if (! (toneLog_.landed() && hpLog_.landed()))
                {
                    toneLog_.advance (seg); hpLog_.advance (seg);   // the filters hold nothing: only where they land
                    applyCorners();                                  // matters, not the way there
                }
            }
            else
            {
                runSegment (io, numChannels, done, seg, invN);
            }

            done += seg;
            if (grid_.advance (seg))
                maintain();
        }

        locked_.store (isLockedNow(), std::memory_order_relaxed);
        return true;
    }

    /** The glided levels where they stand now — taps for tests (owner thread). */
    float sub1Level() const noexcept { return sub1_.value(); }
    float sub2Level() const noexcept { return sub2_.value(); }
    /** The corners where their glides stand now, Hz (owner thread). */
    float toneHzNow()  const noexcept { return (float) core::det::exp2 ((double) toneLog_.value()); }
    float subHpHzNow() const noexcept { return (float) core::det::exp2 ((double) hpLog_.value()); }

    /** True while the sub plays (locked, not dying) — as of the end of the last process() call. A relaxed atomic:
        any thread. */
    bool isLocked() const noexcept { return locked_.load (std::memory_order_relaxed); }
    /** True while the stage is idle and computes nothing (see the contract). Owner thread. */
    bool isIdle() const noexcept { return cold_ && atRestLevels(); }
    /** Accepted rising edges with the gate open since the tracker last started cold. Owner thread. */
    std::int64_t edgeCount() const noexcept { return edges_; }

    /** True when any of the stage's own recursive words sits in the subnormal band — the law-8 witness for the
        tests. (The six eq::Svf filters keep their words private; eq's own suite and LAW8-AUDIT.md answer for
        them, and this stage flushes them on the same grid.) */
    bool anySubnormalState() const noexcept
    {
        const auto subn = [] (float v) { return std::fpclassify (v) == FP_SUBNORMAL; };
        return subn (tilt_) || subn (prevT_) || subn (peakPos_) || subn (peakNeg_) || subn (cycleMax_)
            || subn (amp_) || subn (ampTarget_) || subn (gateEnv_.stateWord()) || subn (source_)
            || subn (freshNow_) || subn (freshBefore_);
    }

private:
    static constexpr double kButterworthQ = 0.70710678118654752;

    bool atRestLevels() const noexcept { return sub1_.landedAt (0.0f) && sub2_.landedAt (0.0f); }

    bool isLockedNow() const noexcept
    {
        return gateOpen_ && edgesSinceLoss_ >= 2 && fallsSinceLoss_ >= 2 && sinceEdge_ <= lockLimit_ && sinceEdge_ < dyingAt_;
    }

    void applyGate() noexcept
    {
        gateOn_  = (float) core::det::pow10 ((double) gateDb_ / 20.0);
        gateOff_ = gateOn_ * kGateCloseRatio;
        floor_   = gateOn_ * kFloorUnderGate;
    }

    /** The voices' corners from where their glides stand. */
    void applyCorners() noexcept
    {
        const double tone = core::det::exp2 ((double) toneLog_.value());
        const double hp   = core::det::exp2 ((double) hpLog_.value());
        tone1_.setParams (eq::FilterType::LowPass,  tone, kButterworthQ, 0.0);
        tone2_.setParams (eq::FilterType::LowPass,  tone, kButterworthQ, 0.0);
        hp1_.setParams   (eq::FilterType::HighPass, hp,   kButterworthQ, 0.0);
        hp2_.setParams   (eq::FilterType::HighPass, hp,   kButterworthQ, 0.0);
    }

    /** One sample of the corner glides: re-designed every sample while they travel, so the filter moves on a curve
        and never on stairs. */
    void stepCorners() noexcept
    {
        if (toneLog_.landed() && hpLog_.landed())
            return;
        (void) toneLog_.next();
        (void) hpLog_.next();
        applyCorners();
    }

    /** Out of idle: a restart. The grid is re-anchored, the corners land (a cold filter has nothing to retune
        under), and the tracker hears the input through its door. */
    void wake() noexcept
    {
        cold_ = false;
        grid_.reset();
        toneLog_.snap(); hpLog_.snap();
        applyCorners();
        door_.snapTo (0.0f);
        door_.setTarget (1.0f);
    }

    /** The tracker and both voices back to a stage that has heard nothing: the stage is COLD. */
    void dropTracker() noexcept
    {
        if (cold_)
            return;   // idempotent

        for (auto* f : { &hpBand_, &lpBand_, &tone1_, &tone2_, &hp1_, &hp2_ })
            f->reset();
        gateEnv_.reset();
        tilt_ = prevT_ = peakPos_ = peakNeg_ = cycleMax_ = amp_ = ampTarget_ = 0.0f;
        high_ = false;
        q1_ = q2_ = -1.0f;
        period_ = lockLimit_ = maxPeriod_;
        intervals_[0] = intervals_[1] = intervals_[2] = 0;
        dyingAt_ = maxPeriod_ + 2;   // unreachable: the counter stops at maxPeriod_ + 1
        fallsSinceLoss_ = 0;
        sawLow_ = false;
        source_ = sourcePeak_ = 0.0f;
        sourceSide_ = 0;
        sourceRose_ = sourceFell_ = false;
        freshNow_ = freshBefore_ = 0.0f;
        freshPos_ = 0;
        sinceEdge_ = maxPeriod_ + 1;
        edgesSinceLoss_ = 0;
        lockPos_ = 0;
        gateOpen_ = false;
        gateHold_ = 0;
        edges_ = 0;
        voicesQuiet_ = true;
        door_.snapTo (0.0f);
        cold_ = true;
    }

    /** The silent fixed point: a silent sample would leave every tracker and voice word exactly as it is. The
        decaying words must sit on +0; the cycle's running peak only grows, so a silent sample leaves it alone, and
        it need only be a value the grid's flush would not move either. */
    bool trackerAtRest() const noexcept
    {
        using detail::plusZero;
        return plusZero (tilt_) && plusZero (prevT_) && plusZero (peakPos_) && plusZero (peakNeg_)
            && (plusZero (cycleMax_) || cycleMax_ >= 1.0e-15f) && plusZero (amp_)
            && plusZero (ampTarget_) && plusZero (gateEnv_.stateWord())
            && plusZero (source_) && plusZero (freshNow_) && plusZero (freshBefore_)
            && (plusZero (sourcePeak_) || sourcePeak_ >= 1.0e-15f)
            && ! gateOpen_ && gateHold_ == 0 && lockPos_ == 0 && voicesQuiet_ && edgesSinceLoss_ == 0
            && sinceEdge_ > maxPeriod_
            && hpBand_.isAtRest() && lpBand_.isAtRest() && tone1_.isAtRest() && tone2_.isAtRest()
            && hp1_.isAtRest() && hp2_.isAtRest();
    }

    /** Once per grid period of audio time: law 8's flush, and — if both levels stand at zero — going idle. */
    void maintain() noexcept
    {
        if (cold_)
            return;

        for (auto* f : { &hpBand_, &lpBand_, &tone1_, &tone2_, &hp1_, &hp2_ })
            f->flushDenormals();
        gateEnv_.flushDenormals();
        for (float* v : { &tilt_, &prevT_, &peakPos_, &peakNeg_, &cycleMax_, &amp_, &ampTarget_, &source_, &sourcePeak_, &freshNow_, &freshBefore_ })
            core::flushPoison (*v);

        if (atRestLevels())
            dropTracker();   // both sub levels are an exact zero: nothing the tracker holds can be heard again
    }

    /** `len` samples from `offset`. nch == 0 is a pause: silence in, nothing out. */
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

            for (int ch = 0; ch < nch; ++ch)
                io[ch][k] += sub;
        }
    }

    /** One sample of the tracker and the voices; returns the sub. */
    float tick (float x) noexcept
    {
        stepCorners();

        // The bass fundamental range, through the wake's door.
        const float xin = x * door_.next();
        const float b   = lpBand_.processSample (0, hpBand_.processSample (0, xin));

        // The SOURCE's own swing, before any filter (DC taken off at 5 Hz, one pole: it does not ring). A note is only
        // acquired on crossings the input itself has just made: past +/- the gate's closing level, in the same
        // direction, since the last counted crossing of the other one. Filter ringing after a source has stopped is
        // not the source swinging (measured on 0faf723: one smooth 41.2 Hz cycle, its second rise supplied by the
        // band's ringing 2.6 ms after the input ended, 0.338 of sub for 50 ms) — and one cycle of anything swings at
        // most twice one way and once the other.
        // The swing must pass a quarter of the source's own peak since the note was lost (and at least the gate's
        // closing level): the high pass's own tail after a lone cycle is a few per cent of it, a note's swings pass a
        // quarter of its peak every cycle.
        source_ += sourceCoeff_ * (xin - source_);
        const float xs = xin - source_, axs = std::fabs (xs);
        sourcePeak_ = axs > sourcePeak_ ? axs : sourcePeak_;
        const float swing = std::max (gateOff_, kSourceSwing * sourcePeak_);
        if (xs > swing && sourceSide_ <= 0)       { sourceSide_ =  1; sourceRose_ = true; }
        else if (xs < -swing && sourceSide_ >= 0) { sourceSide_ = -1; sourceFell_ = true; }

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
        const float t = tilt_, nt = -t;

        // Each polarity's own peak: instant up, kThresholdReleaseSeconds down.
        peakPos_ = t  > peakPos_ ? t  : peakPos_ * peakRel_;
        peakNeg_ = nt > peakNeg_ ? nt : peakNeg_ * peakRel_;

        const float ab = std::fabs (b);
        cycleMax_ = ab > cycleMax_ ? ab : cycleMax_;

        // The band's FRESH level: its peak over the last half to whole longest period (two buckets). A note is
        // acquired only while this stands at the gate's opening level — not the follower at the instant of a crossing,
        // whose ripple kept a steady note at 0.84 ... 1.28 times the level from ever acquiring, and whose tail after a
        // loud note let the hum under it be acquired 153 ms later (both measured on 0faf723).
        freshNow_ = ab > freshNow_ ? ab : freshNow_;
        if ((freshNow_ > 0.0f || freshBefore_ > 0.0f) && ++freshPos_ >= freshLen_)
        {
            freshBefore_ = freshNow_;
            freshNow_    = 0.0f;
            freshPos_    = 0;
        }
        const bool fresh = (freshNow_ > freshBefore_ ? freshNow_ : freshBefore_) >= gateOn_;

        if (sinceEdge_ <= maxPeriod_)
        {
            ++sinceEdge_;

            if (sinceEdge_ == dyingAt_)
                cycleMax_ = ab;   // a quarter period late: the note is dying, and its next level is measured from here

            if (sinceEdge_ == lockLimit_ + 1)
            {
                // Lost: the next note starts from nothing and must prove itself with two edges.
                edgesSinceLoss_ = 0;
                peakPos_ = peakNeg_ = 0.0f;
                high_ = false;
                cycleMax_ = ab;
                intervals_[0] = intervals_[1] = intervals_[2] = 0;
                fallsSinceLoss_ = 0;
                sawLow_ = false;
                sourceSide_ = 0;
                sourceRose_ = sourceFell_ = false;
                sourcePeak_ = 0.0f;
            }
        }
        const bool dyingBefore = sinceEdge_ >= dyingAt_;

        // The thresholds, and their floor: a fixed one under the gate, and one 20 dB under the band's own level — so
        // hum or hiss under a note never squares up as edges of its own, at any level of the note.
        const float fl = std::max (floor_, kFloorUnderBand * g);
        const float hi = std::max (kHysteresis * peakPos_, fl);
        const float lo = std::max (kHysteresis * peakNeg_, fl);

        bool  edge = false;
        float frac = 1.0f;
        if (! high_)
        {
            // An edge COUNTS only with the gate open — and a note is only ACQUIRED (the two-edge proof) above the
            // gate's opening level: what plays under the gate, hum included, can neither prove a note nor set its
            // period. A crossing that may not count leaves the comparator LOW, so the rise that opens the gate is the
            // note's first edge, not the one a cycle later.
            if (t > hi && gateOpen_ && (edgesSinceLoss_ >= 2 || (fresh && sourceRose_)))
            {
                high_ = true;
                sourceRose_ = false;
                if (sinceEdge_ >= minPeriod_)
                {
                    edge = true;
                    const float span = t - prevT_;   // > 0 when t crossed hi from at or below it on this sample
                    frac = span > 0.0f ? detail::clampf (0.0f, 1.0f, (t - hi) / span) : 1.0f;
                }
            }
        }
        else if (t < -lo)
        {
            high_ = false;
        }

        // ...and its FALLS, counted the same way (gate open, above the gate while acquiring), once per swing. A note
        // opens the fade only after two rises AND two falls since it was lost: one cycle of anything — the smoothest
        // click is a single cycle — has at most two of one and one of the other, so it never opens it.
        if (t < -lo && ! sawLow_ && gateOpen_ && (fallsSinceLoss_ >= 2 || (fresh && sourceFell_)))
        {
            sawLow_ = true;
            sourceFell_ = false;
            if (fallsSinceLoss_ < 2)
                ++fallsSinceLoss_;
        }
        prevT_ = t;

        // The divider; at the edge sample the square carries the box average of its step.
        float s1 = q1_, s2 = q2_;
        if (edge)
        {
            if (edgesSinceLoss_ == 0)
            {
                edgesSinceLoss_ = 1;   // the first edge after a loss: no period yet
            }
            else
            {
                edgesSinceLoss_ = 2;
                intervals_[2] = intervals_[1];
                intervals_[1] = intervals_[0];
                intervals_[0] = sinceEdge_;
            }

            // The LOSS deadline comes from a CONFIRMED period only — two intervals in a row that agree within a
            // quarter. A note's first interval can be a partial attack cycle (measured: 813 samples on a low B whose
            // period is 1555), and a loss drawn from it threw the note away before its next true edge. The DYING
            // deadline (a fade, which the next edge reopens without a new proof) runs from the longest interval seen
            // so far even before that: a two-lobe click — two edges, then nothing — fades after 5/4 of its own
            // interval instead of playing for 1/22 s.
            const int a0 = intervals_[0], a1 = intervals_[1];
            const bool confirmed = a0 > 0 && a1 > 0 && 4 * std::abs (a0 - a1) <= std::max (a0, a1);
            period_    = std::max (a0, std::max (a1, intervals_[2]));
            lockLimit_ = confirmed ? std::min (maxPeriod_, kLostPeriods * period_) : maxPeriod_;
            dyingAt_   = confirmed ? period_ * kDyingQuarters / 4 + 1 : maxPeriod_ + 2;
            sinceEdge_ = 0;
            sawLow_    = false;
            ++edges_;

            q1_ = -q1_;
            const float box = 2.0f * frac - 1.0f;
            s1 = q1_ * box;
            if (q1_ > 0.0f)
            {
                q2_ = -q2_;
                s2  = q2_ * box;
            }

            ampTarget_ = cycleMax_;   // the band's peak over the cycle that just ended
            cycleMax_  = 0.0f;
        }

        // The level: up at once within a cycle, down at an edge, released once the note is lost. While the sub is
        // faded out it is not glided but SET — the voice that reopens speaks at the note it reopens on, never at the
        // one before — and so it is at an edge that ends a dying stretch.
        const bool lost  = edgesSinceLoss_ == 0;
        const bool muted = lockPos_ == 0;
        if (lost)
            ampTarget_ = 0.0f;
        else if (cycleMax_ > ampTarget_)
            ampTarget_ = cycleMax_;
        if (muted || (edge && dyingBefore))
            amp_ = ampTarget_;
        else
            amp_ += (lost ? lostCoeff_ : ampCoeff_) * (ampTarget_ - amp_);

        // The voices: tone, sub HP, then the level and the lock fade (locked and not dying) — so a zero of either is
        // an exact zero. Faded out, they hold nothing: they reopen from rest on the edge that reopens them, and no
        // filter carries the last note's level into the next.
        const bool open = isLockedNow();
        float y1 = 0.0f, y2 = 0.0f;
        if (muted && ! open)
        {
            if (! voicesQuiet_)
            {
                tone1_.reset(); tone2_.reset(); hp1_.reset(); hp2_.reset();
                voicesQuiet_ = true;
            }
        }
        else
        {
            if (voicesQuiet_)
            {
                // Reopening: each voice is put where a long half-cycle of the OTHER polarity would have left it, so
                // its first output is a full step of the note it reopens on — on an edge that is the edge itself;
                // off one (a note proven on its second fall, mid half-cycle) the sub enters with a step there, its
                // first half-cycle a little short, rather than waiting half a cycle for the next edge to be heard.
                tone1_.primeDc (0, -q1_ * amp_); hp1_.primeDc (0, -q1_ * amp_);
                tone2_.primeDc (0, -q2_ * amp_); hp2_.primeDc (0, -q2_ * amp_);
            }
            voicesQuiet_ = false;
            y1 = hp1_.processSample (0, tone1_.processSample (0, s1 * amp_));
            y2 = hp2_.processSample (0, tone2_.processSample (0, s2 * amp_));
        }

        if (open) { if (lockPos_ < lockLen_) ++lockPos_; }
        else      { if (lockPos_ > 0)        --lockPos_; }
        const float lock = lockPos_ == lockLen_ ? 1.0f : (float) lockPos_ * lockStep_;

        return (sub1_.next() * y1 + sub2_.next() * y2) * lock;
    }

    // The tracker's filters, the voices' tone and high pass, the gate's follower.
    eq::DeterministicSvf hpBand_, lpBand_, tone1_, tone2_, hp1_, hp2_;
    dynamics::EnvelopeFollower gateEnv_;

    // Knobs.
    detail::Ramp sub1_, sub2_, door_, toneLog_, hpLog_;
    float gateDb_ = -50.0f;
    float gateOn_ = 0.0f, gateOff_ = 0.0f, floor_ = 0.0f;

    // Coefficients and counts at the prepared rate.
    double fs_ = 48000.0;
    float  tiltCoeff_ = 0.0f, peakRel_ = 0.0f, ampCoeff_ = 0.0f, lostCoeff_ = 0.0f, lockStep_ = 1.0f;
    int    minPeriod_ = 1, maxPeriod_ = 1, lockLen_ = 1, gateHoldLen_ = 1;

    // The tracker.
    float tilt_ = 0.0f, prevT_ = 0.0f, peakPos_ = 0.0f, peakNeg_ = 0.0f;
    float cycleMax_ = 0.0f, amp_ = 0.0f, ampTarget_ = 0.0f;
    float q1_ = -1.0f, q2_ = -1.0f;
    bool  high_ = false, gateOpen_ = false, voicesQuiet_ = true, sawLow_ = false;
    float source_ = 0.0f, sourcePeak_ = 0.0f, freshNow_ = 0.0f, freshBefore_ = 0.0f, sourceCoeff_ = 0.0f;
    int   sourceSide_ = 0, freshPos_ = 0, freshLen_ = 1;
    bool  sourceRose_ = false, sourceFell_ = false;
    int   intervals_[3] { 0, 0, 0 };   // the last three edge intervals, newest first (0: none yet)
    int   fallsSinceLoss_ = 0;
    int   sinceEdge_ = 2, period_ = 1, lockLimit_ = 1, dyingAt_ = 3, edgesSinceLoss_ = 0, lockPos_ = 0, gateHold_ = 0;
    std::int64_t edges_ = 0;

    std::atomic<bool> locked_ { false };
    core::StateGrid grid_;
    bool cold_ = false, primed_ = false, prepared_ = false;
    int  maxBlock_ = 0, maxCh_ = 0;
};

} // namespace felitronics::octave

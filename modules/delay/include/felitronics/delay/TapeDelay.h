// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/core/DetMath.h>
#include <felitronics/core/FlushToZero.h>
#include <felitronics/delay/DelayLaw.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

namespace felitronics::delay
{

//==============================================================================
// felitronics::delay::TapeDelay — the echo before the space: repeats of what the preamp made, which a
// reverb then rooms. Moved here from OrbitAmp (its `DelayStage`) so the amp plugins and the browser demo
// run ONE implementation. The signal path is the app's, expression for expression; the transcendentals are
// core::det rather than the system libm, so every row computes the same bits from the same arithmetic — the
// distance to the app's own class is measured by a parity harness and is ulp-level (see README.md).
//
// A TAPE HEAD, NOT A SAMPLER. The delay TIME is a target the line glides toward, and while it glides the
// read head moves at the wrong speed — the repeats bend in pitch, the way a tape echo bends when its motor
// is dragged. Turning the knob IS the effect; there is no crossfade hiding the move.
//
// THE LOOP IS DARK ON PURPOSE. Every pass through the line — the first echo included — goes through a
// one-pole low pass and a soft saturator, the record head's own voice: repeats darken and compress as they
// recirculate instead of piling up as clean copies. DARK is the low-pass corner; the saturator is fixed by
// taste, not a knob.
//
// OFFSET spreads the wet: one channel's repeats arrive later than the other's by a signed handful of
// milliseconds — the whole comb shifts as one, a constant width, not a progressive drift. The offset line
// glides like the time does, so turning it bends too. Dry is untouched: unity always, the wet ADDED.
//
// UP TO kMaxHeads PLAYBACK HEADS on the one tape, each at a POSITION that is a fraction of TIME in (0, 1], with
// a LEVEL and a PAN. TIME stays the one motor: a head stands at position x the motor's glided time, so turning
// TIME moves every head with it, and a head moved on its own glides under TIME's law and bends that repeat. The
// LOOP does not read a head: REPEATS feeds back from 1x TIME, so the heads' figure repeats once per TIME and
// decays with REPEATS — heads at 1/3, 2/3 and 1 play a triplet, and the third echo can be the loudest. PAN is
// a BALANCE on the head's own lines, at constant power with the centre at unity on both sides (+3 dB on its
// side at hard pan); on dual mono — one guitar on both lines, as the amps feed it — heads panned apart are a
// ping-pong. A mono call aims every head at its level, pan aside. The default layout, one head at 1x, level
// 1, centre, is the one-head stage to the bit (TapeDelayParityTests, against a frozen copy).
//
// CONTRACT (DSP-ARCHITECTURE.md §2).
//   · prepare(sampleRate, maxBlock, maxChannels) allocates every line; process()/reset() never allocate,
//     lock or throw. maxChannels is 1 or 2 and the rate 1 kHz … 768 kHz — anything else is REFUSED, and a
//     refused prepare DISARMS the instance: process() refuses until a prepare succeeds (law 11b).
//   · process() is [[nodiscard]] bool (law 11): an unprepared instance, a width above the prepared one, or a
//     negative extent is refused before anything moves; any n >= 0 is processed in full (the stage has no
//     block-sized scratch, so maxBlock is not a limit here).
//   · A NARROWER call is legal; a channel at an index >= nch with n > 0 has STOPPED (law 11a): its line,
//     offset line and loop filter are dropped once, on that edge, and its offset head snaps to its target
//     when it returns. nch == 0 with n > 0 is a gap: the shared motor glides on through it as through
//     silence (law 11c) and every channel's memory is dropped.
//   · reset() is a stream RESTART: the lines cleared, the write clocks re-anchored at zero and the heads
//     snapped to their targets on the next call — what follows a reset does not depend on what came before.
//   · FEED (process's last argument) is the insert's stand-by, and it is a DOOR, not a switch: what enters the
//     record head is ramped linearly to nothing, or back, over kFeedRampSeconds — a flip on a loud note would
//     otherwise be recorded as a step that every repeat replays as a tick. Only the input is ramped, never the
//     recirculation, and a held `feed` records exactly what the switch did. The door is shared by the channels;
//     a pause moves it as silence would (law 11c), and after prepare() or reset() it stands where the next
//     call's `feed` puts it, with no travel, so a fresh stage does not fade in its first note.
//   · Law 8: the loop filter's state and the offset glide are flushed to zero below 1e-15 PER SAMPLE
//     (core::flushDenormal), so no path reaches the subnormal band — the line is written only through the
//     loop filter, and a glide whose target is zero would otherwise park a few subnormals short of it. Per
//     sample, so the flush cannot depend on how the caller slices the stream. The motor and the offset heads
//     glide in double and LAND: once an update no longer moves a head, it is put on its target (LAW8-AUDIT.md,
//     item 6) — the app's float heads parked up to 1/(2·coeff) ulps short, 2.3 ms at 6 s. The playback heads'
//     positions glide the same way (their floor is the smallest normal float, so no glide heads for zero),
//     every contribution a head scales on its way into the offset line is flushed like the loop filter, and a
//     gain's ramp value is floored at the smallest normal float (a ramp restarted again and again from where the
//     last one stood shrinks geometrically, and nothing else would stop it).
//   · THE HEADS' CONTROLS ARE SAMPLE-CLOCKED like the motor: a position glides once per sample, LEVEL and PAN
//     reach the head as two side gains, each on a counted kHeadRampSeconds ramp that lands on its target, and a
//     pause steps all of them as silence would. A head that is silent stands on its target, so turning one on
//     plays it where it was sent instead of sweeping it in. The gains are aimed at the call's width: a mono call
//     ramps side 0 to the head's level, a stereo call back to its balance — at the edge, not as a switch. After
//     prepare() or reset() the positions snap with the motor, and the gains — with the positions again, since no
//     head has been heard until then — on the first call that carries a channel (a pause has no width to aim
//     them for).
//   · POISON: a non-finite INPUT sample is not recorded (the dry keeps it — it is not this stage's to
//     change); a non-finite loop state at the end of a call restarts the line from silence. Knob values are
//     bounded before they reach any arithmetic: a non-finite one is ignored, a finite one is clamped. A head's
//     three fields are bounded each on its own, and a head index off the array is ignored.
//   · THREADS: prepare, reset, the setters and process belong to ONE thread (the audio thread, or whichever
//     owns the stage). tailSeconds, shownTimeMs, shownOffsetMs, shownHeadMs and inputPeak are relaxed atomics
//     that any other thread may read.
//   · Math: core::det (exp2, log2, pow) for the coefficients and the tail, and a det-built tanh per sample in
//     the loop — ~5 ns a call, measured against the system tanh's ~2: one call per sample per channel. The pan
//     law is rational (one correctly rounded sqrt, in the setter), so it needs no det call.
//==============================================================================

/** How many playback heads one tape carries. Every per-head state is an array this long, so a tape with more
    heads is this number and nothing else. At namespace scope: the knob struct below is declared before the class. */
inline constexpr int kMaxHeads = 4;

/** One playback head. A default-constructed head is OFF (level 0). */
struct TapeDelayHead
{
    float position = 1.0f;   ///< a fraction of TIME, (0, 1]: 1 is the full delay, 1/3 a third of it
    float level    = 0.0f;   ///< linear, 0 … 1; 0 is off
    float pan      = 0.0f;   ///< the balance on the head's own lines, -1 left … +1 right, constant power
};

using TapeDelayHeads = std::array<TapeDelayHead, kMaxHeads>;

/** The knobs, in the units the stage takes — the WHOLE state: `TapeDelay::setParams` applies every field, the
    heads included, in this order, through the same guarded setters a caller may use one by one. A struct built
    without heads carries the default layout, so a caller that sets the heads with `setHead` and then hands over a
    fresh struct every block hands the default layout back every block: use one door or the other. */
struct TapeDelayParams
{
    float timeMs   = 350.0f;   ///< the delay's destination, ms (1 … 6100); free, or from `timeMsFor()`
    float repeats  = 0.35f;    ///< feedback 0..1 — the dark filter and the saturator keep even 1 bounded
    float darkHz   = 4000.0f;  ///< the loop low-pass corner, Hz (`darkHzForPercent()` maps the DARK knob)
    float offsetMs = 0.0f;     ///< signed: positive holds the RIGHT repeats back, negative the left (±31)
    float mix      = 0.25f;    ///< 0 = dry only, 1 = the heads' repeats added at their own levels; dry never moves
    /// the playback heads; the default is ONE head at 1x TIME, level 1, centre — the one-head stage to the bit
    TapeDelayHeads heads { { TapeDelayHead { 1.0f, 1.0f, 0.0f } } };
};

namespace detail
{
    /** Equal within one relative float epsilon (or the smallest normal, absolutely) — the guard every setter
        uses, so a knob that re-sends the value it already holds costs nothing. Non-finite: exact only. */
    inline bool nearlyEqual (float a, float b) noexcept
    {
        if (! (std::isfinite (a) && std::isfinite (b)))   // exact only: the same infinity, never a NaN
            return std::isinf (a) && std::isinf (b) && std::signbit (a) == std::signbit (b);

        const float diff = std::abs (a - b);
        return diff <= std::numeric_limits<float>::min()
            || diff <= std::numeric_limits<float>::epsilon() * std::max (std::abs (a), std::abs (b));
    }

    inline float clampf (float lo, float hi, float v) noexcept { return v < lo ? lo : (hi < v ? hi : v); }

    /** EXACTLY equal, for the heads' guards (both values finite and bounded by then). nearlyEqual's absolute
        branch would take a request for exactly 0 for its 1e-38 neighbour and keep the neighbour, and a level
        that never reaches 0 never turns its head off. Spelt without == (-Wfloat-equal); ±0 are one value. */
    inline bool sameValue (float a, float b) noexcept { return ! (a < b) && ! (b < a); }
    inline bool sameValue (double a, double b) noexcept { return ! (a < b) && ! (b < a); }

    inline constexpr double kLog2E = 1.4426950408889634074;

    /** e^x, the same bits on every row: core::det::exp2 of a pinned product. */
    inline double detExp (double x) noexcept { return core::det::exp2 (core::det::mul (x, kLog2E)); }

    /** tanh, the same bits on every row, for the loop's saturator. Below 1/16 a Taylor polynomial (its next
        term is < 1e-14 relative there), above it 1 − 2/(e^{2|x|} + 1), which no longer cancels; every
        multiply-add pinned against contraction (core::det::mulAdd). Accurate far past float in both. */
    inline float detTanh (float xf) noexcept
    {
        const double x = (double) xf;
        const double a = std::fabs (x);

        if (a < 0.0625)
        {
            const double x2 = core::det::mul (x, x);
            double p = core::det::mulAdd (x2, 62.0 / 2835.0, -17.0 / 315.0);
            p = core::det::mulAdd (x2, p, 2.0 / 15.0);
            p = core::det::mulAdd (x2, p, -1.0 / 3.0);
            return (float) core::det::mulAdd (core::det::mul (x, x2), p, x);
        }

        if (a > 20.0)
            return x > 0.0 ? 1.0f : -1.0f;

        const double e = core::det::exp2 (core::det::mul (2.0 * a, kLog2E));
        const double t = 1.0 - 2.0 / (e + 1.0);
        return (float) (x < 0.0 ? -t : t);
    }
} // namespace detail

class TapeDelay
{
public:
    static constexpr int    kMaxChannels      = 2;
    static constexpr double kMaxDelaySeconds  = 6.1;     // 1/1 at 40 BPM, with margin
    static constexpr double kMaxOffsetSeconds = 0.031;   // the offset knob's reach
    static constexpr float  kLineMaxMs = (float) (kMaxDelaySeconds * 1000.0);   // the STAGE's clamp on TIME, 6100 ms
                                                                                 // (the KNOB's range is DelayLaw's kTimeMaxMs)
    static constexpr float  kMinDarkHz = 1.0f, kMaxDarkHz = 100000.0f;          // a corner outside is clamped
    static constexpr double kMinSampleRate = 1000.0, kMaxSampleRate = 768000.0;

    /** How long `feed` takes to open or close the record head. A linear ramp of T is the switch's step smeared
        over T: its spectrum is the step's times sinc(f·T), so above 1/T it falls at 12 dB/oct instead of 6. At
        5 ms that corner is 200 Hz and the edge is 36 dB under the step's by 4 kHz; on a full-scale tone the
        repeats' largest step falls from 13-16x the tone's own to 0.98-1.05x (TapeDelayTests). Short on purpose:
        the door still answers the hand, and a note struck on the switch loses only its first 5 ms to the fade.
        Counted in samples at the prepared rate (240 at 48 kHz, 5 at 1 kHz). */
    static constexpr double kFeedRampSeconds = 0.005;

    static constexpr int kMaxHeads = delay::kMaxHeads;

    /** A head's position floor: the smallest normal float, so a position never heads for zero (a glide toward
        zero would park subnormals short of it). Any position under one sample of delay plays at one sample —
        the motor's own floor — which at 48 kHz is 20.8 µs. */
    static constexpr float kMinHeadPosition = std::numeric_limits<float>::min();

    /** How long a head's LEVEL or PAN takes to arrive: each of its two side gains travels a counted linear ramp of
        this length and lands on its target exactly (960 samples at 48 kHz). Linear in GAIN: a jump from hard left
        to hard right dips the head's power by up to 3 dB mid-ramp; a hand on the knob moves far less a ramp. */
    static constexpr double kHeadRampSeconds = 0.02;

    TapeDelay() noexcept
    {
        setHeads (TapeDelayParams {}.heads);
        for (int k = 0; k < kMaxHeads; ++k)
            headPos_[k] = (double) headWant_[k].position;
        publishHeadTargets();
        refreshTail();
    }

    /** Allocates every line for `maxChannels` (1 or 2) at `sampleRate` (1 kHz … 768 kHz). `maxBlock` (>= 1) is
        recorded, not a limit: the stage has no block-sized scratch. A refused prepare DISARMS the instance. */
    [[nodiscard]] bool prepare (double newSampleRate, int maxBlock, int maxChannels)
    {
        prepared_ = false;   // law 11b: whatever happens next, the old preparation no longer stands

        if (! (newSampleRate >= kMinSampleRate && newSampleRate <= kMaxSampleRate) || maxBlock < 1
            || maxChannels < 1 || maxChannels > kMaxChannels)
            return false;

        sampleRate = newSampleRate;
        maxBlock_  = maxBlock;
        maxCh_     = maxChannels;

        // Sized for the ranges, not the moment: the slowest division at the slowest BPM (a whole note at
        // 40) is six seconds, and the offset rides on top.
        const auto lineMax = (std::size_t) std::ceil (kMaxDelaySeconds * sampleRate) + 8;
        const auto offMax  = (std::size_t) std::ceil (kMaxOffsetSeconds * sampleRate) + 8;

        for (std::size_t ch = 0; ch < (std::size_t) kMaxChannels; ++ch)
        {
            const bool used = (int) ch < maxCh_;
            line[ch].assign (used ? lineMax : 0, 0.0f);
            offLine[ch].assign (used ? offMax : 0, 0.0f);
            lpfState[ch]   = 0.0f;
            live_[ch]      = false;
            offPrimed_[ch] = true;
        }

        lineN_ = (int) lineMax;
        offN_  = (int) offMax;

        writePos = 0;
        offPos   = 0;
        primed   = false;
        cleared  = true;

        // The glides, seeded: a fresh line starts ON its targets rather than sweeping in from wherever the
        // last session left the heads.
        currentDelay = targetDelaySamples();
        currentOff[0] = currentOff[1] = 0.0;

        updateCoeffs();
        refreshTail();

        // The door's travel; `primed` is false, so the first call settles it on its own `feed`.
        feedLen_  = std::max (1, (int) (kFeedRampSeconds * sampleRate + 0.5));
        feedStep_ = 1.0f / (float) feedLen_;

        // The heads, seeded like the motor: on their targets, the gains aimed on the first call that carries
        // a channel.
        headRampLen_  = std::max (1, (int) (kHeadRampSeconds * sampleRate + 0.5));
        headRampStep_ = 1.0 / (double) headRampLen_;
        for (int k = 0; k < kMaxHeads; ++k)
        {
            headPos_[k] = (double) headWant_[k].position;
            for (int s = 0; s < 2; ++s)
            {
                headGain_[k][s] = headGainFrom_[k][s] = headGainTo_[k][s] = 0.0f;
                headGainLeft_[k][s] = 0;
            }
        }
        gainsPrimed_ = false;
        headsMono_   = false;
        publishHeadTargets();

        prepared_ = true;
        return true;
    }

    bool isPrepared()  const noexcept { return prepared_; }
    int  maxChannels() const noexcept { return prepared_ ? maxCh_ : 0; }
    int  maxBlock()    const noexcept { return prepared_ ? maxBlock_ : 0; }
    static constexpr int latencySamples() noexcept { return 0; }

    /** A stream restart: the lines cleared, the write clocks back at zero, the heads re-snapped and the door
        settled on the next call. Cheap-idempotent — a stage out of the rig may call it every block; the clear
        runs once. */
    void reset() noexcept
    {
        if (cleared)
            return;   // an off block calls this every block; the memset must not be per-block

        for (std::size_t ch = 0; ch < (std::size_t) kMaxChannels; ++ch)
        {
            std::fill (line[ch].begin(), line[ch].end(), 0.0f);
            std::fill (offLine[ch].begin(), offLine[ch].end(), 0.0f);
            lpfState[ch] = 0.0f;
            live_[ch]    = false;
        }

        writePos = 0;   // law 11a: a restart re-anchors the clocks, so the read fractions after it do not
        offPos   = 0;   // depend on how long the stream ran before it
        primed  = false;
        cleared = true;
        gainsPrimed_ = false;   // the heads snap with the motor, their gains on the next call that carries a channel
    }

    /** The whole state: the five knobs in the order the app sets them, then the heads. */
    void setParams (const TapeDelayParams& p) noexcept
    {
        setTimeMs (p.timeMs);
        setRepeats (p.repeats);
        setDarkHz (p.darkHz);
        setOffsetMs (p.offsetMs);
        setMix (p.mix);
        setHeads (p.heads);
    }

    /** Head `k` (0 … kMaxHeads-1; any other index is ignored). Each field is bounded on its own — a non-finite one
        leaves that field as it was, a finite one is clamped: POSITION to [kMinHeadPosition, 1], LEVEL to [0, 1]
        (under 1e-15 it is 0: off), PAN to [-1, 1]. Re-sending what the head already holds changes nothing; a new
        position glides there under TIME's law, a new level or pan ramps there over kHeadRampSeconds. */
    void setHead (int k, const TapeDelayHead& h) noexcept
    {
        if (k < 0 || k >= kMaxHeads)
            return;

        TapeDelayHead want = headWant_[k];
        if (std::isfinite (h.position))
            want.position = detail::clampf (kMinHeadPosition, 1.0f, h.position);
        if (std::isfinite (h.level))
        {
            want.level = detail::clampf (0.0f, 1.0f, h.level);
            core::flushDenormal (want.level);   // a level of 1e-38 is a head that never turns off
        }
        if (std::isfinite (h.pan))
            want.pan = detail::clampf (-1.0f, 1.0f, h.pan);

        const bool moved = ! detail::sameValue (want.position, headWant_[k].position);
        const bool aimed = ! detail::sameValue (want.level, headWant_[k].level)
                        || ! detail::sameValue (want.pan, headWant_[k].pan);
        headWant_[k] = want;

        if (aimed)
        {
            aimHead (k);

            // The tail follows what the heads add up to — and only that: a silent head's pan, or the default
            // re-sent, moves no estimate.
            const double peak = headPeakGain();
            if (! detail::sameValue (peak, headPeak_))
            {
                headPeak_ = peak;
                refreshTail();
            }
        }
        if ((moved || aimed) && ! primed)
            publishHeadTargets();   // a running stage publishes where its heads stand, after each call
    }

    /** Every head at once — `setHead` for each, in index order. */
    void setHeads (const TapeDelayHeads& heads) noexcept
    {
        for (int k = 0; k < kMaxHeads; ++k)
            setHead (k, heads[(std::size_t) k]);
    }

    /** The head's destination, in milliseconds — free, or `timeMsFor()`. The line glides there; big moves
        bend. Clamped to 1 … 6100 ms (the line's reach); a non-finite value is ignored. */
    void setTimeMs (float ms) noexcept
    {
        if (! std::isfinite (ms))
            return;

        const float want = detail::clampf (1.0f, kLineMaxMs, ms);

        if (! detail::nearlyEqual (want, timeMs))
        {
            timeMs = want;
            refreshTail();
        }

        shownTimeMs_.store (timeMs, std::memory_order_relaxed);   // process() overwrites with the glide
        if (! primed)
            publishHeadTargets();   // the heads stand at a fraction of TIME, so the picture's heads move with it
    }

    /** Feedback, 0..1. The dark filter and the saturator live inside the loop, so even 1 is a long
        compressed bloom rather than a runaway. */
    void setRepeats (float amount) noexcept
    {
        if (! std::isfinite (amount))
            return;

        const float want = detail::clampf (0.0f, 1.0f, amount);

        if (detail::nearlyEqual (want, repeats))
            return;

        repeats = want;
        refreshTail();
    }

    /** The loop low-pass corner — LOWER is darker, and every pass darkens again. Clamped to 1 Hz … 100 kHz. */
    void setDarkHz (float hz) noexcept
    {
        if (! std::isfinite (hz))
            return;

        const float want = detail::clampf (kMinDarkHz, kMaxDarkHz, hz);

        if (! detail::nearlyEqual (want, darkHz))
        {
            darkHz = want;
            updateCoeffs();
        }
    }

    /** The stereo shift of the wet, signed: positive holds the RIGHT repeats back, negative the left. Dry
        never moves. */
    void setOffsetMs (float ms) noexcept
    {
        if (! std::isfinite (ms))
            return;

        const float before = offsetMs;
        offsetMs = detail::clampf (-1000.0f * (float) kMaxOffsetSeconds, 1000.0f * (float) kMaxOffsetSeconds, ms);
        shownOffsetMs_.store (offsetMs, std::memory_order_relaxed);

        if (! detail::nearlyEqual (before, offsetMs))
            refreshTail();
    }

    /** 0 = fully dry, 1 = the heads' repeats added at their own levels (the default head's at unity). Dry never
        moves. */
    void setMix (float newMix) noexcept
    {
        if (std::isfinite (newMix))
            mix = detail::clampf (0.0f, 1.0f, newMix);
    }

    /** `feed` false is the INSERT'S BYPASS: the record head stops taking new signal, but the recirculation
        keeps going round, so the repeats already in the line ring out into the dry instead of being chopped.
        The line is only cleared by `reset()` — the difference between standing a pedal by and unplugging it.
        A flip is not instant: the input fades out (or in) over kFeedRampSeconds, from wherever a reversed flip
        left it. In place; channels at index >= 2 are never part of a call (a wider call is refused). */
    [[nodiscard]] bool process (float* const* channels, int numChannels, int numSamples, bool feed = true) noexcept
    {
        if (! prepared_ || numChannels < 0 || numSamples < 0 || numChannels > maxCh_)
            return false;   // law 11: refused before anything moves
        if (numSamples == 0)
            return true;    // the only true no-op: no time, no edge
        if (numChannels > 0 && channels == nullptr)
            return false;

        const int nch = numChannels;

        // Law 11a — the channels this call does not carry have stopped: their memory goes, once.
        for (int ch = nch; ch < maxCh_; ++ch)
            stopChannel (ch);

        // Waking up, the heads SNAP to their targets: a first block should not spend six seconds sweeping
        // in from wherever the knob stood when the block was last on.
        const float tgtDelay = targetDelaySamples();
        if (! primed)
        {
            currentDelay = tgtDelay;
            currentOff[0] = offsetMs < 0.0f ? -offsetMs * 0.001f * (float) sampleRate : 0.0f;
            currentOff[1] = offsetMs > 0.0f ?  offsetMs * 0.001f * (float) sampleRate : 0.0f;
            offPrimed_[0] = offPrimed_[1] = true;
            feedPos_ = feed ? feedLen_ : 0;   // ...and the door stands where this call wants it: no fade-in
            for (int k = 0; k < kMaxHeads; ++k)
                headPos_[k] = (double) headWant_[k].position;   // the playback heads snap with the motor
            primed  = true;
        }
        cleared = false;

        // The heads' gains are aimed at the width of a call that carries a channel; a pause has none.
        if (nch > 0)
            aimWidth (nch == 1);

        // AT REST: every head on its target and every ramp landed. The setters run between calls and nothing
        // inside one moves a target, so heads at rest stay at rest to the call's end, and stepHeads() would be
        // an exact no-op on every sample of it (a snap onto the target a head stands on, a glide of zero that
        // lands where it is, ramps that are done) — skipped. And when the one head sounding is a 1x head at a
        // gain of exactly 1 on every side the call reads — the default layout — the heads play the loop's own
        // tap, which is what the general sum would make of it.
        const bool atRest = headsAtRest();
        const bool plain  = atRest && plainLoopHead (nch);

        const float tgtOff[2] = {
            offsetMs < 0.0f ? -offsetMs * 0.001f * (float) sampleRate : 0.0f,
            offsetMs > 0.0f ?  offsetMs * 0.001f * (float) sampleRate : 0.0f,
        };

        // A channel coming back after it stopped starts its offset head ON its target, as a fresh line does.
        for (int ch = 0; ch < nch; ++ch)
        {
            if (! offPrimed_[ch])
            {
                currentOff[ch] = tgtOff[ch];
                offPrimed_[ch] = true;
            }
            live_[ch] = true;
        }

        if (nch == 0)
        {
            // Law 11c — a pause is silence for the shared motor: it glides on as it would through n samples
            // of digital silence. Every channel stopped above, so there is no line to read or write.
            for (int i = 0; i < numSamples; ++i)
            {
                glide (currentDelay, tgtDelay, glideCoeff);
                if (! atRest)
                    stepHeads();   // ...and so do the heads' positions and gains
            }

            feedPos_ = feedAfter (numSamples, feed);   // the door travels on, as through silence
            writePos = (int) (((long long) writePos + numSamples) % lineN_);
            offPos   = (int) (((long long) offPos   + numSamples) % offN_);
            envIn_.store (0.0f, std::memory_order_relaxed);
            publishHeads();
            return true;
        }

        const int lineN  = lineN_;
        const int offN   = offN_;
        const int feedTo = feed ? feedLen_ : 0;

        // The comb's pulse listens at the door: the block peak of what the line is being FED, taken before
        // the loop overwrites the channel in place.
        {
            float pk = 0.0f;
            for (int i = 0; i < numSamples; ++i)
                pk = std::max (pk, std::abs (channels[0][i]));
            envIn_.store (pk, std::memory_order_relaxed);
        }

        for (int i = 0; i < numSamples; ++i)
        {
            // The motor: one exponential glide shared by both channels, so the comb never splits while it
            // travels.
            glide (currentDelay, tgtDelay, glideCoeff);
            const float dly = detail::clampf (1.0f, (float) (lineN - 4), (float) currentDelay);

            // The playback heads: one step of their controls, then where each sounding head reads — once per
            // sample, shared by the channels like the motor. A head at 1x reads where the loop does, so it is
            // the loop's own tap; any other stands at position x the motor's glided time, floored at one sample
            // like the motor, and reads with the loop tap's own arithmetic below.
            if (! atRest)
                stepHeads();
            // Per sample, as locals: measured faster with four heads than one set per call, and no slower on the
            // default path.
            int   h0[kMaxHeads] {}, h1[kMaxHeads] {};
            float hf[kMaxHeads] {};
            bool  hLoop[kMaxHeads] {};
            for (int k = 0; k < kMaxHeads && ! plain; ++k)
            {
                hLoop[k] = headPos_[k] >= 1.0;   // the glide never overshoots its target, and no target is past 1
                if (hLoop[k] || ! (headGain_[k][0] > 0.0f || headGain_[k][1] > 0.0f))
                    continue;

                const float hd = detail::clampf (1.0f, (float) (lineN - 4), (float) (headPos_[k] * currentDelay));
                float hpos = (float) writePos - hd;
                if (hpos < 0.0f)
                    hpos += (float) lineN;
                const int hwhole = (int) hpos;   // A HAIR SHORT OF THE END, as for the loop tap below
                hf[k] = hpos - (float) hwhole;
                h0[k] = hwhole == lineN ? 0 : hwhole;
                h1[k] = h0[k] + 1 == lineN ? 0 : h0[k] + 1;
            }

            // The door: one sample of travel toward `feed`, shared by both channels like the motor.
            if (feedPos_ != feedTo)
                feedPos_ += feedPos_ < feedTo ? 1 : -1;
            const bool  feedOpen = feedPos_ == feedLen_;
            const bool  feedShut = feedPos_ == 0;
            const float feedGain = (float) feedPos_ * feedStep_;   // read only strictly between the ends

            for (int ch = 0; ch < nch; ++ch)
            {
                const float x = channels[ch][i];
                float* const ln = line[(std::size_t) ch].data();
                float* const ol = offLine[(std::size_t) ch].data();

                // The play head, fractional: this read at a moving distance IS the repitch.
                float pos = (float) writePos - dly;
                if (pos < 0.0f)
                    pos += (float) lineN;

                // A HAIR SHORT OF THE END IS THE END. A read a hair behind a write head that has just wrapped
                // lands a hair short of the line's length, and a float that size cannot hold the hair: it
                // rounds to the length itself, one past the last sample. That position is the line's first
                // sample, and it is read as one.
                const int   whole = (int) pos;
                const float f  = pos - (float) whole;
                const int   i0 = whole == lineN ? 0 : whole;
                const int   i1 = i0 + 1 == lineN ? 0 : i0 + 1;
                const float tap = ln[i0] * (1.0f - f)
                                + ln[i1] * f;

                // What the heads play on this channel, read before the record head writes: each sounding head's
                // tap on this side, in index order. A contribution at a gain of exactly 1 is the tap itself — the
                // default head is the one-head stage's `tap`, bit for bit — and any other is scaled, then flushed
                // (law 8: a landed ramp still scales, and a tiny product would reach the offset line's band). The
                // sum starts from the first contribution; with none, the heads play +0 and the loop runs on. A
                // mono call reads side 0, which aimWidth has aimed at each head's level.
                float heard = plain ? tap : 0.0f;
                bool  sounding = false;
                for (int k = 0; k < kMaxHeads && ! plain; ++k)
                {
                    const float g = headGain_[k][ch];
                    if (! (g > 0.0f))
                        continue;

                    const float ht = hLoop[k] ? tap
                                              : ln[h0[k]] * (1.0f - hf[k])
                                              + ln[h1[k]] * hf[k];
                    float c = ht;
                    if (! (g >= 1.0f && g <= 1.0f))
                    {
                        c = ht * g;
                        core::flushDenormal (c);
                    }
                    heard = sounding ? heard + c : c;
                    sounding = true;
                }

                // The record head: input plus the recirculation, darkened and pressed — the first echo
                // already wears one pass, the n-th wears n. A sample that is not a number is not recorded;
                // the dry keeps it, it is not ours. The input comes through the door — the plain sample when
                // it is open and a plain zero when it is shut, so a held `feed` records the same bits it
                // always did, scaled only while the door travels. The recirculation does not pass the door.
                const float in = feedOpen ? x : (feedShut ? 0.0f : x * feedGain);
                float v = (std::isfinite (x) ? in : 0.0f) + tap * repeats;
                lpfState[ch] += lpfCoeff * (v - lpfState[ch]);
                core::flushDenormal (lpfState[ch]);   // law 8: the only way into the line is through here
                v = detail::detTanh (lpfState[ch] * satDrive) * satNorm;
                ln[writePos] = v;

                // The offset line rides the WET only — one channel's comb held back whole. It runs on its OWN
                // write position: it is a short line, and the long line's clock would walk past its end.
                glide (currentOff[ch], tgtOff[ch], offGlideCoeff);
                core::flushDenormal (currentOff[ch]);   // law 8: a glide to zero parks subnormals short of it
                ol[offPos] = heard;

                float opos = (float) offPos - detail::clampf (0.0f, (float) (offN - 4), (float) currentOff[ch]);
                if (opos < 0.0f)
                    opos += (float) offN;

                // ...the same edge: a knob's zero can arrive as a hair of a sample (0.000000447 ms from a host
                // on arm64), so a wrap of this line would read past its end — see A HAIR SHORT OF THE END.
                const int   owhole = (int) opos;
                const float of = opos - (float) owhole;
                const int   o0 = owhole == offN ? 0 : owhole;
                const int   o1 = o0 + 1 == offN ? 0 : o0 + 1;
                const float wet = ol[o0] * (1.0f - of)
                                + ol[o1] * of;

                channels[ch][i] = x + wet * mix;
            }

            writePos = writePos + 1 == lineN ? 0 : writePos + 1;
            offPos   = offPos   + 1 == offN  ? 0 : offPos   + 1;
        }

        // POISON. The line is recursive, and a NaN or an infinity recorded into it would go round for ever.
        // The door above keeps one from being recorded; this is the proof nothing else got in. Every write
        // into the line goes through the loop filter's state, so a finite state at the end of the call means
        // a clean line, and anything else restarts it from silence.
        if (! std::isfinite (lpfState[0]) || ! std::isfinite (lpfState[1]))
            reset();

        publishHeads();
        return true;
    }

    /** HOW LONG THE LINE GOES ON SOUNDING after the last note went in, in seconds — what a plugin tells a
        host that renders offline. A recirculating line loses `repeats` of itself every time round; the
        echoes are discrete, and the n-th pass leaves at n·t carrying A·g^(n-1) — A the heads' summed first-echo
        level (refreshTail), 1 for the default layout — so the first under a thousandth is
        `1 + ln(0.001/A)/ln(g)` of them, at least one. At zero it is one echo; at one it never decays and the answer is a
        million seconds — the caller caps it. Deliberately the LONG estimate: the darkening and the
        saturator bleed the loop further. Correct from construction on, before any knob is touched. */
    float tailSeconds()   const noexcept { return tailSec_.load (std::memory_order_relaxed); }

    /** The picture's taps — the glided time and offset while the stage runs (the set targets while it does
        not), and the block peak of what the line was fed. */
    float shownTimeMs()   const noexcept { return shownTimeMs_.load (std::memory_order_relaxed); }
    float shownOffsetMs() const noexcept { return shownOffsetMs_.load (std::memory_order_relaxed); }
    float inputPeak()     const noexcept { return envIn_.load (std::memory_order_relaxed); }

    /** Where playback head `k` stands, in ms of delay before OFFSET — its glided position times the glided time,
        floored at one sample, as it reads; its target while the stage does not run. For the picture: each tap is
        its own atomic, not one snapshot of all four. An index off the array reads 0. */
    float shownHeadMs (int k) const noexcept
    {
        return k >= 0 && k < kMaxHeads ? shownHeadMs_[k].load (std::memory_order_relaxed) : 0.0f;
    }

    /** Where the heads stand, in samples — the motor's delay and channel `ch`'s offset (0 or 1). The tests read
        them to see a glide LAND on its target, bit for bit. */
    double headSamples() const noexcept { return currentDelay; }
    double offsetHeadSamples (int ch) const noexcept { return currentOff[ch < 0 ? 0 : (ch > 1 ? 1 : ch)]; }

    /** Playback head `k`'s position (a fraction of TIME) and its side gain `side` (0 or 1) as the controls stand
        now — the tests' witnesses for the glide and the ramps landing. Between a restart and the first call that
        carries a channel the gains are not aimed yet: after prepare() they read 0, after reset() what they were;
        a pause in that window puts the positions on their targets. An index off the array reads 0. */
    double headPosition (int k) const noexcept { return k >= 0 && k < kMaxHeads ? headPos_[k] : 0.0; }
    float headGain (int k, int side) const noexcept
    {
        return k >= 0 && k < kMaxHeads && side >= 0 && side < 2 ? headGain_[k][side] : 0.0f;
    }

    /** True when any state sits in the subnormal band — the law-8 witness for the tests: the lines and the offset
        lines, the loop filters, the motor and the offset glides, and per playback head its position, its gains as
        heard, both ends of each ramp and its balance targets. */
    bool anySubnormalState() const noexcept
    {
        const auto sub = [] (auto v) { return std::fpclassify (v) == FP_SUBNORMAL; };
        for (std::size_t ch = 0; ch < (std::size_t) kMaxChannels; ++ch)
        {
            if (sub (lpfState[ch]) || sub (currentOff[ch]))
                return true;
            for (float v : line[ch])
                if (sub (v))
                    return true;
            for (float v : offLine[ch])
                if (sub (v))
                    return true;
        }
        for (int k = 0; k < kMaxHeads; ++k)
            for (int s = 0; s < 2; ++s)
                if (sub (headPos_[k]) || sub (headGain_[k][s]) || sub (headGainFrom_[k][s]) || sub (headGainTo_[k][s])
                    || sub (headGainWant_[k][s]))
                    return true;
        return sub (currentDelay);
    }

private:
    /** One step of an exponential glide that LANDS. The heads glide in DOUBLE: in float the update stops
        moving the value once (target − current)·coeff falls under half an ulp of `current`, which is
        1/(2·coeff) ulps short — measured on the app's float motor, 112 samples (2.3 ms) short of a 6 s target
        at 48 kHz, parked there for good. In double that distance is ~1e-6 of a sample, so when the update no
        longer moves the value, putting it on the target is a step nobody can hear. */
    static void glide (double& current, double target, double coeff) noexcept
    {
        const double next = current + (target - current) * coeff;
        current = std::fabs (next - current) > 0.0 ? next : target;
    }

    /** Where the door stands `n` samples on, travelling toward `feed`'s end and stopping there. */
    int feedAfter (int n, bool feed) const noexcept
    {
        const int to   = feed ? feedLen_ : 0;
        const int step = std::min (n, std::abs (to - feedPos_));
        return feedPos_ < to ? feedPos_ + step : feedPos_ - step;
    }

    float targetDelaySamples() const noexcept
    {
        return (float) ((double) timeMs * 0.001 * sampleRate);
    }

    void updateCoeffs() noexcept
    {
        // The motor's inertia: a 150 ms time constant (37 % of the distance still to go after it; the landing
        // comes much later, by glide()'s rule) — enough travel to hear the bend.
        glideCoeff    = 1.0 - detail::detExp (-1.0 / (0.15 * sampleRate));
        offGlideCoeff = 1.0 - detail::detExp (-1.0 / (0.03 * sampleRate));
        lpfCoeff      = 1.0f - (float) detail::detExp (-kTwoPi * (double) darkHz / sampleRate);
    }

    void stopChannel (int ch) noexcept
    {
        if (ch < 0 || ch >= kMaxChannels || ! live_[ch])   // the bound is for the optimiser's array-bounds check
            return;

        std::fill (line[(std::size_t) ch].begin(), line[(std::size_t) ch].end(), 0.0f);
        std::fill (offLine[(std::size_t) ch].begin(), offLine[(std::size_t) ch].end(), 0.0f);
        lpfState[ch]   = 0.0f;
        live_[ch]      = false;
        offPrimed_[ch] = false;
    }

    void publishHeads() noexcept
    {
        // Where the heads actually STAND, for the picture: the glided values, so the comb slides with the
        // motor instead of jumping to the target.
        shownTimeMs_.store ((float) (currentDelay / sampleRate * 1000.0), std::memory_order_relaxed);
        shownOffsetMs_.store ((float) ((currentOff[1] - currentOff[0]) / sampleRate * 1000.0),
                              std::memory_order_relaxed);

        // ...and each playback head where it reads: the same product and floor as the read itself.
        const float hi = (float) (lineN_ - 4);
        for (int k = 0; k < kMaxHeads; ++k)
        {
            const float hd = detail::clampf (1.0f, hi, (float) (headPos_[k] * currentDelay));
            shownHeadMs_[k].store ((float) ((double) hd / sampleRate * 1000.0), std::memory_order_relaxed);
        }
    }

    /** Before the stage runs, the picture shows where the heads are SENT: position x TIME. */
    void publishHeadTargets() noexcept
    {
        for (int k = 0; k < kMaxHeads; ++k)
            shownHeadMs_[k].store ((float) ((double) headWant_[k].position * (double) timeMs), std::memory_order_relaxed);
    }

    /** One sample of every playback head's controls — the audio loop and a pause both run THIS, so a pause moves
        them exactly as silence would. A head silent coming into the sample stands on its target; a sounding one
        glides there under TIME's law and lands. Then each side gain takes one step of its ramp. Until the first
        call that carries a channel aims the gains (after prepare or reset), there is nothing to hear: the heads
        stand on their targets and the ramps wait. */
    void stepHeads() noexcept
    {
        for (int k = 0; k < kMaxHeads; ++k)
        {
            const double want = (double) headWant_[k].position;
            if (! gainsPrimed_ || ! (headGain_[k][0] > 0.0f || headGain_[k][1] > 0.0f))
                headPos_[k] = want;
            else
                glide (headPos_[k], want, glideCoeff);

            if (! gainsPrimed_)
                continue;
            for (int s = 0; s < 2; ++s)
            {
                if (headGainLeft_[k][s] == 0)
                    continue;   // landed: the gain IS its target, no arithmetic
                const int left = --headGainLeft_[k][s];
                const double from = (double) headGainFrom_[k][s], to = (double) headGainTo_[k][s];
                float g = left == 0 ? headGainTo_[k][s]
                                    : (float) core::det::mulAdd (to - from,
                                                                 (double) (headRampLen_ - left) * headRampStep_,
                                                                 from);

                // Law 8. One ramp never steps under target/N, and every target is zero or normal (a level of
                // 1e-15 one float short of hard pan is ~4e-23). But a ramp restarted from where the last one
                // stood — a level sent up and back down before every sample — shrinks the gain geometrically, and
                // only this floor stops it at the edge of the subnormal band. Not the house 1e-15: the pan law
                // makes gains that small on purpose, and a ramp must land on the law's bits. (A gain that small
                // plays nothing either way — every scaled contribution is flushed at 1e-15.)
                if (std::fabs (g) < std::numeric_limits<float>::min())
                    g = 0.0f;
                headGain_[k][s] = g;
            }
        }
    }

    /** Every head standing on its target and every ramp landed — see AT REST in process(). */
    bool headsAtRest() const noexcept
    {
        for (int k = 0; k < kMaxHeads; ++k)
            if (! detail::sameValue (headPos_[k], (double) headWant_[k].position)
                || headGainLeft_[k][0] != 0 || headGainLeft_[k][1] != 0)
                return false;
        return true;
    }

    /** Exactly one head sounds on the sides a call of `nch` channels reads (side 0 alone in a mono call), it reads
        where the loop does (position 1), and it plays them at a gain of exactly 1; every other head is silent on
        them. Then the heads' sum is the loop's tap, unscaled. */
    bool plainLoopHead (int nch) const noexcept
    {
        int unity = 0;
        for (int k = 0; k < kMaxHeads; ++k)
        {
            bool silent = true, one = headPos_[k] >= 1.0;
            for (int s = 0; s < nch && s < 2; ++s)
            {
                const float g = headGain_[k][s];
                silent = silent && ! (g > 0.0f);
                one    = one && g >= 1.0f && g <= 1.0f;
            }
            if (one)
                ++unity;
            else if (! silent)
                return false;
        }
        return unity == 1;
    }

    /** A, the tail's first-echo level: the sum over the heads of level x the louder side gain — at least the
        level, so it covers a mono call too. Exactly 1 for the default layout. */
    double headPeakGain() const noexcept
    {
        double a = 0.0;
        for (int k = 0; k < kMaxHeads; ++k)
            a += (double) std::max (headGainWant_[k][0], headGainWant_[k][1]);
        return a;
    }

    /** Side `s`'s destination for head `k` at the width the gains are aimed for: in a mono call side 0 plays
        the head's level, pan aside; otherwise each side its balance. */
    float gainTarget (int k, int s) const noexcept
    {
        return s == 0 && headsMono_ ? headWant_[k].level : headGainWant_[k][s];
    }

    /** Send side `s` of head `k` on a fresh ramp from where it stands — unless it is already headed there. */
    void retarget (int k, int s) noexcept
    {
        const float to = gainTarget (k, s);
        if (detail::sameValue (to, headGainTo_[k][s]))
            return;   // a re-send must not restart the ramp: that would hang the sound on the caller's blocks
        headGainFrom_[k][s] = headGain_[k][s];
        headGainTo_[k][s]   = to;
        headGainLeft_[k][s] = headRampLen_;
    }

    /** A head's level or pan moved: its two balance gains, from the rational constant-power law
        gL = (1 - p)/sqrt(1 + p²), gR = (1 + p)/sqrt(1 + p²) — gL² + gR² = 2 for every p, exactly (1, 1) at the
        centre and (√2, 0) at the left end, mirror-exact, every operation correctly rounded (p·p is exact in
        double), so the same bits on every row with no det call. It is a rotation by atan(p) rather than p·45°,
        at most 4.1° from the sine law. Then the ramps, if the gains are aimed. */
    void aimHead (int k) noexcept
    {
        const double p = (double) headWant_[k].pan, lv = (double) headWant_[k].level;
        const double r = std::sqrt (1.0 + p * p);
        headGainWant_[k][0] = (float) (lv * ((1.0 - p) / r));
        headGainWant_[k][1] = (float) (lv * ((1.0 + p) / r));

        if (gainsPrimed_)
        {
            retarget (k, 0);
            retarget (k, 1);
        }
    }

    /** The call carries a channel: the gains are aimed for its width. The first such call after prepare() or
        reset() puts every gain ON its target (the stage is waking, like the motor); after that a change of width
        sends side 0 to its new target on a ramp — the mono level or the balance — so the edge has no step. */
    void aimWidth (bool mono) noexcept
    {
        if (! gainsPrimed_)
        {
            headsMono_ = mono;
            for (int k = 0; k < kMaxHeads; ++k)
            {
                // No head has been heard since the restart, so each is SILENT up to here and stands on its target
                // — a position sent after a pause primed the motor would otherwise be glided to from where the
                // pause left it, at full level from the first sample.
                headPos_[k] = (double) headWant_[k].position;
                for (int s = 0; s < 2; ++s)
                {
                    headGain_[k][s] = headGainFrom_[k][s] = headGainTo_[k][s] = gainTarget (k, s);
                    headGainLeft_[k][s] = 0;
                }
            }
            gainsPrimed_ = true;
            return;
        }

        if (mono == headsMono_)
            return;
        headsMono_ = mono;
        for (int k = 0; k < kMaxHeads; ++k)
            retarget (k, 0);
    }

    // The saturator, fixed by taste: gentle at echo level, a press at full recirculation.
    static constexpr float  satDrive = 1.2f;
    static constexpr float  satNorm  = 1.0f / satDrive;
    static constexpr double kTwoPi   = 6.283185307179586476925286766559;

    /** Worked out where the numbers it depends on are set, rather than where it is asked: the asker may be
        any thread and these are the audio thread's. Every setter checks first, so an untouched knob costs no
        logarithm. log2 for ln: the ratio is the same in any base.
        THE HEADS. No head plays past 1x TIME, so the last echo of each pass is no later than the one head's; but
        heads stacked on one place add, and hard pan is +3 dB on its side. So the first echo starts at A, the sum
        over heads of level × the louder side gain, and the count is the passes until A·g^(n-1) falls under a
        thousandth: 1 + log2(0.001/A)/log2(g), at least one. The default layout's A is exactly 1, and 0.001/1 is
        0.001: the one-head stage's estimate, to the bit. */
    void refreshTail() noexcept
    {
        const double t   = (double) timeMs * 0.001;
        const double g   = (double) repeats;
        const double off = std::abs ((double) offsetMs) * 0.001;   // one side arrives this much later

        const double a = headPeak_;

        tailSec_.store ((float) (off + (g <= 0.0 ? t
                                      : g >= 1.0 ? 1.0e6            // never decays; the caller caps it
                                      : a > 0.001 ? t * (1.0 + core::det::log2 (0.001 / a) / core::det::log2 (g))
                                                  : t)),            // the first echo is already under a thousandth
                        std::memory_order_relaxed);
    }

    std::atomic<float> tailSec_       { 0.0f };   // set by the constructor, from the formula
    std::atomic<float> shownTimeMs_   { 350.0f };
    std::atomic<float> shownOffsetMs_ { 0.0f };
    std::atomic<float> envIn_         { 0.0f };

    float timeMs   = 350.0f;
    float repeats  = 0.35f;
    float darkHz   = 4000.0f;
    float offsetMs = 0.0f;
    float mix      = 0.25f;

    double sampleRate = 48000.0;
    double glideCoeff = 0.0, offGlideCoeff = 0.0;   // the heads' glides, in double
    float  lpfCoeff = 1.0f;

    std::array<std::vector<float>, kMaxChannels> line, offLine;
    int   lineN_ = 0, offN_ = 0;
    int   writePos = 0, offPos = 0;
    double currentDelay = 0.0;                        // the heads glide in double — see glide()
    double currentOff[kMaxChannels] { 0.0, 0.0 };
    float lpfState[kMaxChannels]   { 0.0f, 0.0f };
    // The door, 0 = shut … feedLen_ = open. COUNTED, not glided: it lands on its ends exactly, where a glide
    // toward zero would park subnormals short of it (law 8) and a door a hair open would still record.
    int   feedLen_ = 1, feedPos_ = 0;
    float feedStep_ = 1.0f;                            // 1 / feedLen_
    bool  live_[kMaxChannels]      { false, false };   // carried by a call since the last clear (law 11a)
    bool  offPrimed_[kMaxChannels] { true, true };     // false: the offset head snaps on the channel's return
    bool  primed = false, cleared = true;
    bool  prepared_ = false;
    int   maxBlock_ = 0, maxCh_ = 0;

    // The playback heads — RAW arrays, as the rest of the state is: the container type of an operand on the
    // default path was measured to move FMA contraction under gcc (arm64, -ffp-contract=on).
    TapeDelayHead headWant_[kMaxHeads] {};           // the targets, canonical: bounded and flushed by setHead
    float  headGainWant_[kMaxHeads][2] {};           // each target's balance gains, level x the pan law
    double headPos_[kMaxHeads] {};                   // the glided positions, in double like the motor
    float  headGain_[kMaxHeads][2] {};               // the side gains as heard; side 0 is the mono gain in a mono call
    float  headGainFrom_[kMaxHeads][2] {}, headGainTo_[kMaxHeads][2] {};   // each ramp's ends
    int    headGainLeft_[kMaxHeads][2] {};           // samples left on each ramp; 0 = landed on headGainTo_
    int    headRampLen_ = 1;
    double headRampStep_ = 1.0;                      // 1 / headRampLen_
    double headPeak_ = 0.0;                          // headPeakGain() as the tail last read it
    bool   gainsPrimed_ = false;                     // false: the next call that carries a channel puts them on target
    bool   headsMono_ = false;                       // the width the gains are aimed at: the last call that had one
    std::atomic<float> shownHeadMs_[kMaxHeads] {};
};

} // namespace felitronics::delay

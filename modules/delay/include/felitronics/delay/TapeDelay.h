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
//   · Law 8: the loop filter's state and the offset glide are flushed to zero below 1e-15 PER SAMPLE
//     (core::flushDenormal), so no path reaches the subnormal band — the line is written only through the
//     loop filter, and a glide whose target is zero would otherwise park a few subnormals short of it. Per
//     sample, so the flush cannot depend on how the caller slices the stream. Both heads glide in double and
//     LAND: once an update no longer moves a head, it is put on its target (LAW8-AUDIT.md, item 6) — the
//     app's float heads parked up to 1/(2·coeff) ulps short, 2.3 ms at 6 s.
//   · POISON: a non-finite INPUT sample is not recorded (the dry keeps it — it is not this stage's to
//     change); a non-finite loop state at the end of a call restarts the line from silence. Knob values are
//     bounded before they reach any arithmetic: a non-finite one is ignored, a finite one is clamped.
//   · THREADS: prepare, reset, the setters and process belong to ONE thread (the audio thread, or whichever
//     owns the stage). tailSeconds, shownTimeMs, shownOffsetMs and inputPeak are relaxed atomics that any
//     other thread may read.
//   · Math: core::det (exp2, log2, pow) for the coefficients and the tail, and a det-built tanh per sample in
//     the loop — ~5 ns a call, measured against the system tanh's ~2: one call per sample per channel.
//==============================================================================

/** The knobs, in the units the stage takes. `TapeDelay::setParams` applies them in this order, through
    the same guarded setters a caller may use one by one. */
struct TapeDelayParams
{
    float timeMs   = 350.0f;   ///< the head's destination, ms (1 … 6100); free, or from `timeMsFor()`
    float repeats  = 0.35f;    ///< feedback 0..1 — the dark filter and the saturator keep even 1 bounded
    float darkHz   = 4000.0f;  ///< the loop low-pass corner, Hz (`darkHzForPercent()` maps the DARK knob)
    float offsetMs = 0.0f;     ///< signed: positive holds the RIGHT repeats back, negative the left (±31)
    float mix      = 0.25f;    ///< 0 = dry only, 1 = the repeats added at unity; dry never moves
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
    static constexpr float  kMaxTimeMs = (float) (kMaxDelaySeconds * 1000.0);   // a longer TIME is this
    static constexpr float  kMinDarkHz = 1.0f, kMaxDarkHz = 100000.0f;          // a corner outside is clamped
    static constexpr double kMinSampleRate = 1000.0, kMaxSampleRate = 768000.0;

    TapeDelay() noexcept { refreshTail(); }

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
        prepared_ = true;
        return true;
    }

    bool isPrepared()  const noexcept { return prepared_; }
    int  maxChannels() const noexcept { return prepared_ ? maxCh_ : 0; }
    int  maxBlock()    const noexcept { return prepared_ ? maxBlock_ : 0; }
    static constexpr int latencySamples() noexcept { return 0; }

    /** A stream restart: the lines cleared, the write clocks back at zero, the heads re-snapped on the next
        call. Cheap-idempotent — a stage out of the rig may call it every block; the clear runs once. */
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
    }

    /** All five knobs, in the order the app sets them. */
    void setParams (const TapeDelayParams& p) noexcept
    {
        setTimeMs (p.timeMs);
        setRepeats (p.repeats);
        setDarkHz (p.darkHz);
        setOffsetMs (p.offsetMs);
        setMix (p.mix);
    }

    /** The head's destination, in milliseconds — free, or `timeMsFor()`. The line glides there; big moves
        bend. Clamped to 1 … 6100 ms (the line's reach); a non-finite value is ignored. */
    void setTimeMs (float ms) noexcept
    {
        if (! std::isfinite (ms))
            return;

        const float want = detail::clampf (1.0f, kMaxTimeMs, ms);

        if (! detail::nearlyEqual (want, timeMs))
        {
            timeMs = want;
            refreshTail();
        }

        shownTimeMs_.store (timeMs, std::memory_order_relaxed);   // process() overwrites with the glide
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

    /** 0 = fully dry, 1 = the repeats added at unity. Dry never moves. */
    void setMix (float newMix) noexcept
    {
        if (std::isfinite (newMix))
            mix = detail::clampf (0.0f, 1.0f, newMix);
    }

    /** `feed` false is the INSERT'S BYPASS: the record head stops taking new signal, but the recirculation
        keeps going round, so the repeats already in the line ring out into the dry instead of being chopped.
        The line is only cleared by `reset()` — the difference between standing a pedal by and unplugging it.
        In place; channels at index >= 2 are never part of a call (a wider call is refused). */
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
            primed  = true;
        }
        cleared = false;

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
                glide (currentDelay, tgtDelay, glideCoeff);

            writePos = (int) (((long long) writePos + numSamples) % lineN_);
            offPos   = (int) (((long long) offPos   + numSamples) % offN_);
            envIn_.store (0.0f, std::memory_order_relaxed);
            publishHeads();
            return true;
        }

        const int lineN = lineN_;
        const int offN  = offN_;

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

                // The record head: input plus the recirculation, darkened and pressed — the first echo
                // already wears one pass, the n-th wears n. A sample that is not a number is not recorded;
                // the dry keeps it, it is not ours.
                float v = (feed && std::isfinite (x) ? x : 0.0f) + tap * repeats;
                lpfState[ch] += lpfCoeff * (v - lpfState[ch]);
                core::flushDenormal (lpfState[ch]);   // law 8: the only way into the line is through here
                v = detail::detTanh (lpfState[ch] * satDrive) * satNorm;
                ln[writePos] = v;

                // The offset line rides the WET only — one channel's comb held back whole. It runs on its OWN
                // write position: it is a short line, and the long line's clock would walk past its end.
                glide (currentOff[ch], tgtOff[ch], offGlideCoeff);
                core::flushDenormal (currentOff[ch]);   // law 8: a glide to zero parks subnormals short of it
                ol[offPos] = tap;

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
        echoes are discrete, and the n-th leaves at n·t carrying g^(n-1), so the first under a thousandth is
        `1 + ln(0.001)/ln(g)` of them. At zero it is one echo; at one it never decays and the answer is a
        million seconds — the caller caps it. Deliberately the LONG estimate: the darkening and the
        saturator bleed the loop further. Correct from construction on, before any knob is touched. */
    float tailSeconds()   const noexcept { return tailSec_.load (std::memory_order_relaxed); }

    /** The picture's taps — the glided time and offset while the stage runs (the set targets while it does
        not), and the block peak of what the line was fed. */
    float shownTimeMs()   const noexcept { return shownTimeMs_.load (std::memory_order_relaxed); }
    float shownOffsetMs() const noexcept { return shownOffsetMs_.load (std::memory_order_relaxed); }
    float inputPeak()     const noexcept { return envIn_.load (std::memory_order_relaxed); }

    /** True when any recursive state sits in the subnormal band — the law-8 witness for the tests. */
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
        }
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

    float targetDelaySamples() const noexcept
    {
        return (float) ((double) timeMs * 0.001 * sampleRate);
    }

    void updateCoeffs() noexcept
    {
        // The motor's inertia: ~150 ms to the new time — enough travel to hear the bend.
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
    }

    // The saturator, fixed by taste: gentle at echo level, a press at full recirculation.
    static constexpr float  satDrive = 1.2f;
    static constexpr float  satNorm  = 1.0f / satDrive;
    static constexpr double kTwoPi   = 6.283185307179586476925286766559;

    /** Worked out where the numbers it depends on are set, rather than where it is asked: the asker may be
        any thread and these are the audio thread's. Every setter checks first, so an untouched knob costs no
        logarithm. log2 for ln: the ratio is the same in any base. */
    void refreshTail() noexcept
    {
        const double t   = (double) timeMs * 0.001;
        const double g   = (double) repeats;
        const double off = std::abs ((double) offsetMs) * 0.001;   // one side arrives this much later

        tailSec_.store ((float) (off + (g <= 0.0 ? t
                                      : g >= 1.0 ? 1.0e6            // never decays; the caller caps it
                                                 : t * (1.0 + core::det::log2 (0.001) / core::det::log2 (g)))),
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
    bool  live_[kMaxChannels]      { false, false };   // carried by a call since the last clear (law 11a)
    bool  offPrimed_[kMaxChannels] { true, true };     // false: the offset head snaps on the channel's return
    bool  primed = false, cleared = true;
    bool  prepared_ = false;
    int   maxBlock_ = 0, maxCh_ = 0;
};

} // namespace felitronics::delay

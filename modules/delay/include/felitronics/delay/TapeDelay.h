// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

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
// run ONE implementation; the arithmetic is the app's, expression for expression, and a parity harness
// against the app's own class pins that (bit-identical on the rows it was run on — see README.md).
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
//     lock or throw. maxChannels is 1 or 2 — a wider prepare is REFUSED, not clamped (law 11b).
//   · process() is [[nodiscard]] bool (law 11): an unprepared instance, a width above the prepared one, or a
//     negative extent is refused before anything moves; any n >= 0 is processed in full (the stage has no
//     block-sized scratch, so maxBlock is not a limit here).
//   · A NARROWER call is legal; a channel at an index >= nch with n > 0 has STOPPED (law 11a): its line,
//     offset line and loop filter are dropped once, on that edge, and its offset head snaps to its target
//     when it returns. nch == 0 with n > 0 is a gap: the shared motor glides on through it as through
//     silence (law 11c) and every channel's memory is dropped.
//   · Law 8: the loop filter's state and the offset glide are flushed to zero below 1e-30 PER SAMPLE, so no
//     path reaches the subnormal band — the line is written only through the loop filter, and a glide whose
//     target is zero would otherwise park a few subnormals short of it. Per sample, so the flush cannot
//     depend on how the caller slices the stream. Neither flush can move a sample by more than ~1e-29.
//   · POISON: a non-finite INPUT sample is not recorded (the dry keeps it — it is not this stage's to
//     change); a non-finite loop state at the end of a call restarts the line from silence.
//   · The picture taps (tailSeconds, shownTimeMs, shownOffsetMs, inputPeak) are relaxed atomics: the audio
//     thread writes them, any other thread may read them.
//   · libm: std::exp (coefficients), std::tanh (per sample), std::log (the tail estimate) and std::pow (the
//     DARK mapping) are the SYSTEM functions, deliberately — they are the app's, and parity with the shipped
//     sound is the point. Listed in tools/lint/det-math-manifest.txt as retain-rt.
//==============================================================================

/** The knobs, in the units the stage takes. `TapeDelay::setParams` applies them in this order, through
    the same guarded setters a caller may use one by one. */
struct TapeDelayParams
{
    float timeMs   = 350.0f;   ///< the head's destination, ms (>= 1); free, or `divisionMs()` for a synced one
    float repeats  = 0.35f;    ///< feedback 0..1 — the dark filter and the saturator keep even 1 bounded
    float darkHz   = 4000.0f;  ///< the loop low-pass corner, Hz (`darkHzForPercent()` maps the DARK knob)
    float offsetMs = 0.0f;     ///< signed: positive holds the RIGHT repeats back, negative the left (±31)
    float mix      = 0.25f;    ///< 0 = dry only, 1 = the repeats added at unity; dry never moves
};

//------------------------------------------------------------------------------
// The parameter law of the shipped amps — the sync ladder and the DARK knob — so a host that drives this
// stage (a plugin, the web demo) lands on the same numbers the app does.
//------------------------------------------------------------------------------

/** The sync ladder, longest first, dotted and triplet beside each plain value. */
inline constexpr int kNumDivisions = 13;
inline constexpr const char* kDivisionNames[kNumDivisions] = { "1/1", "1/2.", "1/2", "1/2T", "1/4.", "1/4",
                                                               "1/4T", "1/8.", "1/8", "1/8T", "1/16.", "1/16",
                                                               "1/16T" };
/** The same ladder in quarter notes: ms = beats × 60000 / BPM. */
inline constexpr float kDivisionBeats[kNumDivisions] = { 4.0f, 3.0f, 2.0f, 4.0f / 3.0f, 1.5f, 1.0f,
                                                         2.0f / 3.0f, 0.75f, 0.5f, 1.0f / 3.0f, 0.375f,
                                                         0.25f, 1.0f / 6.0f };
inline constexpr int   kDivisionDefault = 5;                       // 1/4
inline constexpr float kBpmMin = 40.0f, kBpmMax = 240.0f, kBpmDefault = 120.0f;
inline constexpr float kTimeMinMs = 20.0f, kTimeMaxMs = 2000.0f;   // the free TIME knob's range
inline constexpr float kOffsetMaxMs = 30.0f;                       // the OFFSET knob's reach
inline constexpr float kDarkLoHz = 800.0f, kDarkHiHz = 16000.0f;   // DARK's floor and ceiling
inline constexpr float kDarkDefaultPct = 45.0f;                    // ~4 kHz — the tape-ish middle

/** A synced division's time at `bpm`, in ms. The division index is clamped to the ladder; `bpm` must be
    positive (the caller chooses whose tempo — the host's or its own field). */
inline float divisionMs (int division, double bpm) noexcept
{
    const int d = division < 0 ? 0 : (division >= kNumDivisions ? kNumDivisions - 1 : division);
    return (float) ((double) kDivisionBeats[d] * 60000.0 / bpm);
}

/** DARK is a percent of darkness, so the knob's direction tells the truth: clockwise darkens. The corner
    falls log-evenly from the bright ceiling to the dark floor, every degree darkening by the same ear-step. */
inline float darkHzForPercent (float percent) noexcept
{
    return kDarkHiHz * std::pow (kDarkLoHz / kDarkHiHz, percent * 0.01f);
}

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

    /** Law 8, per sample: exact zero below 1e-30. Not the house 1e-15 (core::flushDenormal): the loop
        integrates what it flushes, so the threshold sets how far the deep tail departs from the unflushed
        reference — 1.5e-15 at 1e-15, measured. 1e-30 is still eight decades clear of the subnormal floor,
        and no product in the loop (coefficients >= ~1e-4 at any sane rate) can reach it. */
    inline void flushTiny (float& x) noexcept { if (std::fabs (x) < 1.0e-30f) x = 0.0f; }
} // namespace detail

class TapeDelay
{
public:
    static constexpr int    kMaxChannels      = 2;
    static constexpr double kMaxDelaySeconds  = 6.1;     // 1/1 at 40 BPM, with margin
    static constexpr double kMaxOffsetSeconds = 0.031;   // the offset knob's reach
    static constexpr float  kMinDarkHz = 1.0f, kMaxDarkHz = 100000.0f;   // a corner outside is clamped

    /** Allocates every line for `maxChannels` (1 or 2) at `sampleRate`. `maxBlock` (>= 1) is recorded, not
        a limit: the stage has no block-sized scratch. A refused prepare leaves the instance as it was. */
    [[nodiscard]] bool prepare (double newSampleRate, int maxBlock, int maxChannels)
    {
        if (! (newSampleRate > 0.0) || ! std::isfinite (newSampleRate) || maxBlock < 1
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
            lpfState[ch] = 0.0f;
            live_[ch]    = false;
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
        currentOff[0] = currentOff[1] = 0.0f;

        updateCoeffs();
        prepared_ = true;
        return true;
    }

    bool isPrepared()  const noexcept { return prepared_; }
    int  maxChannels() const noexcept { return prepared_ ? maxCh_ : 0; }
    int  maxBlock()    const noexcept { return prepared_ ? maxBlock_ : 0; }
    static constexpr int latencySamples() noexcept { return 0; }

    /** A stream restart: the lines cleared and the heads re-snapped on the next call. Cheap-idempotent — a
        stage out of the rig may call it every block; the clear runs once. */
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

    /** The head's destination, in milliseconds — free or from `divisionMs()`. The line glides there; big
        moves bend. Below 1 ms is 1 ms; a non-finite value is ignored. */
    void setTimeMs (float ms) noexcept
    {
        if (! std::isfinite (ms))
            return;

        const float want = std::max (1.0f, ms);

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

    /** The loop low-pass corner — LOWER is darker, and every pass darkens again. */
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
                currentDelay += (tgtDelay - currentDelay) * glideCoeff;

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
            currentDelay += (tgtDelay - currentDelay) * glideCoeff;
            const float dly = detail::clampf (1.0f, (float) (lineN - 4), currentDelay);

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
                detail::flushTiny (lpfState[ch]);   // law 8: the only way into the line is through here
                v = std::tanh (lpfState[ch] * satDrive) * satNorm;
                ln[writePos] = v;

                // The offset line rides the WET only — one channel's comb held back whole. It runs on its OWN
                // write position: it is a short line, and the long line's clock would walk past its end.
                currentOff[ch] += (tgtOff[ch] - currentOff[ch]) * offGlideCoeff;
                detail::flushTiny (currentOff[ch]);   // law 8: a glide to zero parks subnormals short of it
                ol[offPos] = tap;

                float opos = (float) offPos - detail::clampf (0.0f, (float) (offN - 4), currentOff[ch]);
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
        saturator bleed the loop further. */
    float tailSeconds()   const noexcept { return tailSec_.load (std::memory_order_relaxed); }

    /** The picture's taps — the glided time and offset while the stage runs (the set targets while it does
        not), and the block peak of what the line was fed. */
    float shownTimeMs()   const noexcept { return shownTimeMs_.load (std::memory_order_relaxed); }
    float shownOffsetMs() const noexcept { return shownOffsetMs_.load (std::memory_order_relaxed); }
    float inputPeak()     const noexcept { return envIn_.load (std::memory_order_relaxed); }

    /** True when any recursive state sits in the subnormal band — the law-8 witness for the tests. */
    bool anySubnormalState() const noexcept
    {
        const auto sub = [] (float v) { return std::fpclassify (v) == FP_SUBNORMAL; };
        for (int ch = 0; ch < kMaxChannels; ++ch)
        {
            if (sub (lpfState[ch]) || sub (currentOff[ch]))
                return true;
            for (float v : line[(std::size_t) ch])
                if (sub (v))
                    return true;
        }
        return sub (currentDelay);
    }

private:
    float targetDelaySamples() const noexcept
    {
        return (float) ((double) timeMs * 0.001 * sampleRate);
    }

    void updateCoeffs() noexcept
    {
        // The motor's inertia: ~150 ms to the new time — enough travel to hear the bend.
        glideCoeff    = 1.0f - std::exp ((float) (-1.0 / (0.15 * sampleRate)));
        offGlideCoeff = 1.0f - std::exp ((float) (-1.0 / (0.03 * sampleRate)));
        lpfCoeff      = 1.0f - std::exp ((float) (-kTwoPi * (double) darkHz / sampleRate));
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
        logarithm. */
    void refreshTail() noexcept
    {
        const double t   = (double) timeMs * 0.001;
        const double g   = (double) repeats;
        const double off = std::abs ((double) offsetMs) * 0.001;   // one side arrives this much later

        tailSec_.store ((float) (off + (g <= 0.0 ? t
                                      : g >= 1.0 ? 1.0e6            // never decays; the caller caps it
                                                 : t * (1.0 + std::log (0.001) / std::log (g)))),
                        std::memory_order_relaxed);
    }

    std::atomic<float> tailSec_       { 0.35f };
    std::atomic<float> shownTimeMs_   { 350.0f };
    std::atomic<float> shownOffsetMs_ { 0.0f };
    std::atomic<float> envIn_         { 0.0f };

    float timeMs   = 350.0f;
    float repeats  = 0.35f;
    float darkHz   = 4000.0f;
    float offsetMs = 0.0f;
    float mix      = 0.25f;

    double sampleRate = 48000.0;
    float  glideCoeff = 0.0f, offGlideCoeff = 0.0f, lpfCoeff = 1.0f;

    std::array<std::vector<float>, kMaxChannels> line, offLine;
    int   lineN_ = 0, offN_ = 0;
    int   writePos = 0, offPos = 0;
    float currentDelay = 0.0f;
    float currentOff[kMaxChannels] { 0.0f, 0.0f };
    float lpfState[kMaxChannels]   { 0.0f, 0.0f };
    bool  live_[kMaxChannels]      { false, false };   // carried by a call since the last clear (law 11a)
    bool  offPrimed_[kMaxChannels] { true, true };     // false: the offset head snaps on the channel's return
    bool  primed = false, cleared = true;
    bool  prepared_ = false;
    int   maxBlock_ = 0, maxCh_ = 0;
};

} // namespace felitronics::delay

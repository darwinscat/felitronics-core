// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/core/DetMath.h>
#include <felitronics/core/StateGrid.h>
#include <felitronics/eq/Svf.h>
#include <felitronics/reverb/Freeverb.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <vector>

namespace felitronics::reverb
{

//==============================================================================
// felitronics::reverb::AmpReverb — the reverb tail of the amp plugins: the last of the amp's own layers
// before the signal leaves for the speaker. Moved here from OrbitAmp (its `ReverbStage`, which sat on the
// framework's Freeverb and therefore had to live in the product); it now runs on `reverb::Freeverb`, so the
// plugins and the browser demo share one implementation. The arithmetic is the app's, expression for
// expression, and a parity harness against the app's own class pins that (see README.md).
//
// Mix is the hero, the CHARACTER is the title, and everything else is a late refinement: DECAY scales the
// chosen character's tail rather than replacing it, PREDELAY holds the tail back so the attack stays dry,
// and the HPF cleans the WET only — a low tail is mud, and every layer after this one can only pass it on.
//
// THE DRY PATH IS THE AMP'S LAW: an amp's reverb control ADDS the tank return to the dry signal — it does
// not crossfade away from it. Dry stays at unity at every setting; the whole wet chain (room, modulation,
// predelay, HPF) runs on a scratch copy and is added on top.
//
// Two characters are more than presets: AMBIENCE is the room you don't hear as an effect; MODULATED wears a
// slow chorus on its tail (the Lexicon dress).
//
// CONTRACT (DSP-ARCHITECTURE.md §2).
//   · prepare(sampleRate, maxBlock, maxChannels) allocates; process()/reset() never allocate, lock or throw.
//     maxChannels is 1 or 2 and the rate 1 kHz … 768 kHz — anything else is REFUSED, and a refused prepare
//     DISARMS the instance: process() refuses until a prepare succeeds (law 11b). A prepared room starts AT
//     its set parameters: no glide in from anywhere.
//   · process() is [[nodiscard]] bool (law 11). Any n >= 0 is processed IN FULL, in chunks of
//     min(maxBlock, 8192) — the scratch is capped, so a whole-file maxBlock costs nothing — the same as the
//     caller chunking at those boundaries (11a); `addedWet` then holds the LAST chunk.
//   · A NARROWER call is legal: one channel runs the room's mono path. A channel at an index >= nch with
//     n > 0 has STOPPED (11a): its combs, allpasses, predelay, chorus line and HPF state are dropped once, on
//     that edge. nch == 0 with n > 0 is a gap (11c): the room's parameter glides and the chorus LFO run on
//     through it as through silence, and every channel's memory is dropped.
//   · reset() is a stream RESTART: every line cleared, every clock re-anchored (the line positions, the chorus
//     LFO's phase, the law-8 grid) and every parameter glide landed on its target — what follows a reset does
//     not depend on what came before it.
//   · Law 8: the room flushes its own state to zero below 1e-15 per sample (Freeverb.h); the predelay and
//     chorus lines are pure delays; the HPF's integrators are flushed on `core::StateGrid` (audio time, 64
//     samples) and healed at the end of every call.
//   · POISON: a non-finite input sample never enters the wet chain (the dry keeps it); a block whose tail is
//     not finite is dropped whole and the room restarts from silence, the dry passing untouched. Knob values
//     are bounded before they reach any arithmetic: a non-finite one is ignored, a finite one clamped.
//   · THREADS: prepare, reset, the setters and process belong to ONE thread (the audio thread, or whichever
//     owns the stage). tailSeconds() is a relaxed atomic any thread may read; addedWet()/addedLength() are the
//     process thread's — read them there, between calls.
//   · Math: core::det (log2 for DECAY and the tail, sin for the chorus LFO — per sample, ~5 ns a call), so
//     every row computes the same coefficients and the same LFO.
//==============================================================================

enum class Character { ambience, room, hall, plate, spring, modulated };

inline constexpr int kNumCharacters = 6;
inline constexpr const char* kCharacterNames[kNumCharacters] = { "Ambience", "Room", "Hall", "Plate", "Spring",
                                                                 "Modulated" };

/** The knobs, in the units the stage takes. `AmpReverb::setParams` applies them in this order, through the
    same guarded setters a caller may use one by one. */
struct AmpReverbParams
{
    Character character  = Character::room;
    float     mix        = 0.2f;     ///< 0 = dry only, 1 = the tail added at unity; dry never moves
    float     decay      = 1.0f;     ///< scales the character's decay, 0.5..2
    float     predelayMs = 0.0f;     ///< 0..100 ms
    float     hpfHz      = 120.0f;   ///< the tail's own high-pass corner (the amps offer 40..500 Hz)
};

class AmpReverb
{
public:
    static constexpr int kMaxChannels = 2;

    static constexpr double kMinSampleRate = 1000.0, kMaxSampleRate = 768000.0;
    static constexpr int    kMaxChunk      = 8192;   // the scratch's cap: longer calls are chunked here

    AmpReverb() noexcept { refreshTail(); }

    [[nodiscard]] bool prepare (double newSampleRate, int maxBlock, int maxChannels)
    {
        prepared_ = false;   // law 11b: whatever happens next, the old preparation no longer stands

        if (! (newSampleRate >= kMinSampleRate && newSampleRate <= kMaxSampleRate)
            || maxBlock < 1 || maxChannels < 1 || maxChannels > kMaxChannels)
            return false;

        sampleRate = newSampleRate;
        maxBlock_  = std::min (maxBlock, kMaxChunk);
        maxCh_     = maxChannels;

        // The room's parameters FIRST, then its prepare: prepare() lands every glide on the targets set
        // here, so a fresh room starts at its character instead of gliding in from Freeverb's published
        // initial parameters (whose dry gain would leak the input into the tail for the first 10 ms).
        hpf.prepare (sampleRate, kMaxChannels);
        applyHpf();
        apply();
        (void) room_.prepare (sampleRate);   // cannot refuse: the rate was checked above
        room_.reset();

        for (auto& b : wet)
            b.assign ((std::size_t) maxBlock_, 0.0f);
        for (auto& b : wetOut)
            b.assign ((std::size_t) maxBlock_, 0.0f);

        const auto preMax = (std::size_t) std::ceil (0.1 * sampleRate) + 8;   // 100 ms of predelay
        const auto modMax = (std::size_t) std::ceil (0.02 * sampleRate) + 8;  // 12 ms base + depth

        for (std::size_t ch = 0; ch < (std::size_t) kMaxChannels; ++ch)
        {
            preLine[ch].assign (preMax, 0.0f);
            modLine[ch].assign (modMax, 0.0f);
            live_[ch] = false;
        }

        prePos = modPos = 0;
        lfoPhase = 0.0f;
        addedCount = 0;

        hpf.reset();
        grid_.reset();

        cleared   = true;
        prepared_ = true;
        return true;
    }

    bool isPrepared()  const noexcept { return prepared_; }
    int  maxChannels() const noexcept { return prepared_ ? maxCh_ : 0; }
    int  maxBlock()    const noexcept { return prepared_ ? maxBlock_ : 0; }
    static constexpr int latencySamples() noexcept { return 0; }

    /** A stream restart. Idempotent on purpose: a chain calls this every block a room is out of the rig, and
        the clear underneath is eight combs, four allpasses and two delay lines per channel — about a hundred
        and fifty kilobytes — which must not be a per-block memset. */
    void reset() noexcept
    {
        room_.landGlides();   // a restart starts AT the set parameters, like a fresh prepare — even on a cleared room

        if (cleared)
            return;

        cleared = true;
        room_.reset();
        hpf.reset();
        grid_.reset();
        for (std::size_t ch = 0; ch < (std::size_t) kMaxChannels; ++ch)
        {
            std::fill (preLine[ch].begin(), preLine[ch].end(), 0.0f);
            std::fill (modLine[ch].begin(), modLine[ch].end(), 0.0f);
            live_[ch] = false;
        }
        prePos = modPos = 0;   // law 11a: the clocks are re-anchored with the lines
        lfoPhase = 0.0f;
    }

    /** All five knobs, in the order the app sets them. */
    void setParams (const AmpReverbParams& p) noexcept
    {
        setCharacter (p.character);
        setMix (p.mix);
        setDecay (p.decay);
        setPredelayMs (p.predelayMs);
        setHpfHz (p.hpfHz);
    }

    void setCharacter (Character c) noexcept
    {
        if ((int) c < 0 || (int) c >= kNumCharacters || c == character)
            return;

        character = c;
        apply();
    }

    /** 0 = fully dry, 1 = the tail added at unity. Dry never moves. */
    void setMix (float newMix) noexcept
    {
        if (std::isfinite (newMix))
            mix = clampf (0.0f, 1.0f, newMix);
    }

    /** Scales the character's decay, 0.5..2 — the character stays the voice, this is its breath. */
    void setDecay (float scale) noexcept
    {
        if (! std::isfinite (scale))
            return;

        const float s = clampf (0.5f, 2.0f, scale);
        if (! detail::nearlyEqual (s, decayScale))
        {
            decayScale = s;
            apply();
        }
    }

    /** How long the attack stays dry before the tail arrives, 0..100 ms. */
    void setPredelayMs (float ms) noexcept
    {
        if (! std::isfinite (ms))
            return;

        const float want = clampf (0.0f, 100.0f, ms);

        if (detail::nearlyEqual (want, predelayMs))
            return;

        predelayMs = want;
        refreshTail();
    }

    /** The tail's own high-pass — the WET only, always in: a low tail is mud in any speaker, and at a 40 Hz
        floor the filter is as good as air. The corner is bounded by the filter itself to [1 Hz, 0.49·rate]. */
    void setHpfHz (float hz) noexcept
    {
        if (! std::isfinite (hz))
            return;

        if (! detail::nearlyEqual (hz, hpfHz))
        {
            hpfHz = hz;
            applyHpf();
        }
    }

    /** `feed` false is the INSERT'S BYPASS: nothing new enters the wet chain, so the room that is already
        ringing decays into the dry rather than being cut. The dry is untouched either way; this stage only
        ever ADDS. */
    [[nodiscard]] bool process (float* const* channels, int numChannels, int numSamples, bool feed = true) noexcept
    {
        if (! prepared_ || numChannels < 0 || numSamples < 0 || numChannels > maxCh_)
            return false;   // law 11: refused before anything moves
        if (numSamples == 0)
            return true;    // the only true no-op
        if (numChannels > 0 && channels == nullptr)
            return false;

        // Law 11a — the channels this call does not carry have stopped: their memory goes, once.
        for (int ch = numChannels; ch < maxCh_; ++ch)
            stopChannel (ch);

        if (numChannels == 0)
        {
            pause (numSamples);
            return true;
        }

        for (int done = 0; done < numSamples;)
        {
            const int len = std::min (maxBlock_, numSamples - done);
            float* chunk[kMaxChannels] = { channels[0] + done, numChannels > 1 ? channels[1] + done : nullptr };
            processChunk (chunk, numChannels, len, feed);
            done += len;
        }

        return true;
    }

    /** What the stage ADDED in the last chunk it ran, per channel (0 or 1) — for a picture that shows the
        door and the tail. */
    const float* addedWet (int ch) const noexcept { return wetOut[(std::size_t) std::clamp (ch, 0, 1)].data(); }

    /** How many samples of `addedWet` the last chunk wrote (0 after a pause). A reader must not go past it. */

    /** The room itself, read-only — its glide state is what a pause has to move as silence would. */
    const Freeverb& room() const noexcept { return room_; }

    /** True when the room, the predelay or the chorus line holds a subnormal — the law-8 witness for the tests.
        (The HPF's two integrators per channel are the eq module's; its own flush is tested there.) */
    bool anySubnormalState() const noexcept
    {
        const auto sub = [] (float v) { return std::fpclassify (v) == FP_SUBNORMAL; };
        for (std::size_t ch = 0; ch < (std::size_t) kMaxChannels; ++ch)
        {
            for (float v : preLine[ch]) if (sub (v)) return true;
            for (float v : modLine[ch]) if (sub (v)) return true;
        }
        return room_.anySubnormalState();
    }
    int addedLength() const noexcept { return addedCount; }

    /** HOW LONG THE ROOM GOES ON SOUNDING after the last note went in, in seconds — what a plugin tells a host
        that renders offline. The longest comb is the right channel's, 1617 + 23 samples of 44.1 kHz whatever
        the rate (the tunings scale with it); each pass loses `roomSize × 0.28 + 0.7` of itself, so reaching a
        thousandth takes that comb times `1 + ln(0.001)/ln(g)` — a HALL at DECAY ×2 is a shade under ten
        seconds. The damping takes the top off faster, so this is the long side, the side to be wrong on. */
    float tailSeconds() const noexcept { return tailSec.load (std::memory_order_relaxed); }

private:
    static float clampf (float lo, float hi, float v) noexcept { return v < lo ? lo : (hi < v ? hi : v); }

    void stopChannel (int ch) noexcept
    {
        if (ch < 0 || ch >= kMaxChannels || ! live_[ch])   // the bound is for the optimiser's array-bounds check
            return;

        room_.clearChannel (ch);
        hpf.resetChannel (ch);
        std::fill (preLine[(std::size_t) ch].begin(), preLine[(std::size_t) ch].end(), 0.0f);
        std::fill (modLine[(std::size_t) ch].begin(), modLine[(std::size_t) ch].end(), 0.0f);
        live_[ch] = false;
    }

    /** Law 11c: n samples of digital silence at no width. Every channel has stopped (its memory is gone), so
        what is left to move is shared: the room's glides, the chorus LFO and the line clocks. */
    void pause (int n) noexcept
    {
        cleared = false;
        addedCount = 0;
        room_.advanceGlides (n);

        if (character == Character::modulated)
        {
            const float inc = (float) (0.45 / sampleRate);
            for (int i = 0; i < n; ++i)
            {
                lfoPhase += inc;
                if (lfoPhase >= 1.0f)
                    lfoPhase -= 1.0f;
            }
            modPos = (int) (((long long) modPos + n) % (long long) modLine[0].size());
        }

        prePos = (int) (((long long) prePos + n) % (long long) preLine[0].size());
        grid_.skip (n);
    }

    void processChunk (float* const* channels, int nch, int n, bool feed) noexcept
    {
        cleared = false;   // there is state in here again; the next reset() has work to do
        addedCount = n;

        for (int ch = 0; ch < nch; ++ch)
            live_[ch] = true;   // per chunk: a poison restart below clears the marks mid-call

        // The wet copy: the whole wet chain runs on it, and the dry never enters. A sample that is not a
        // number does not enter it either — the dry keeps whatever it was handed.
        for (int ch = 0; ch < nch; ++ch)
        {
            float* w = wet[(std::size_t) ch].data();
            const float* d = channels[ch];

            for (int i = 0; i < n; ++i)
                w[i] = feed && std::isfinite (d[i]) ? d[i] : 0.0f;
        }

        // The room itself, 100% wet (the wet level 1/3 undoes Freeverb's ×3).
        if (nch >= 2)
            (void) room_.processStereo (wet[0].data(), wet[1].data(), n);   // prepared, planes ours: cannot refuse
        else
            (void) room_.processMono (wet[0].data(), n);

        // The slow chorus on the tail — a modulated delay, the Lexicon dress.
        if (character == Character::modulated)
            modulateTail (nch, n);

        // Predelay, then the tail's own high-pass, then the mix — into `wetOut`, not yet added.
        const int preSamples = (int) ((double) predelayMs * 0.001 * sampleRate);

        for (int i = 0; i < n; ++i)
        {
            const int preN = (int) preLine[0].size();
            const int rd   = (prePos - preSamples + preN) % preN;

            for (int ch = 0; ch < nch; ++ch)
            {
                float* w = wet[(std::size_t) ch].data();
                preLine[(std::size_t) ch][(std::size_t) prePos] = w[i];
                float s = preLine[(std::size_t) ch][(std::size_t) rd];

                s = hpf.processSample (ch, s);

                wetOut[(std::size_t) ch][(std::size_t) i] = s * mix;
            }

            prePos = (prePos + 1) % (int) preLine[0].size();

            if (grid_.advance (1))
                hpf.flushDenormals();   // law 8, on audio time
        }

        hpf.healPoison();   // law 8a(a): the poison half runs at the end of every call

        // POISON. Every element of the room is recursive, and a single NaN or infinity in any of them would
        // reproduce itself for ever. The door above keeps a bad sample from walking in; this catches one the
        // room makes of its own. The whole chunk's tail is checked BEFORE any of it reaches the dry: a tail
        // that is not a number is dropped, the room restarts from silence and the dry passes untouched.
        // `x - x` is 0 for every finite x and NaN for anything else, so one running sum answers.
        float poison = 0.0f;
        for (int ch = 0; ch < nch; ++ch)
            for (int i = 0; i < n; ++i)
                poison += wetOut[(std::size_t) ch][(std::size_t) i] - wetOut[(std::size_t) ch][(std::size_t) i];

        if (! std::isfinite (poison))
        {
            reset();

            for (int ch = 0; ch < nch; ++ch)
                std::fill_n (wetOut[(std::size_t) ch].begin(), n, 0.0f);

            return;
        }

        for (int ch = 0; ch < nch; ++ch)
            for (int i = 0; i < n; ++i)
                channels[ch][i] += wetOut[(std::size_t) ch][(std::size_t) i];
    }

    /** The room the character asks for, as DECAY leaves it — the one place this is worked out. */
    float roomSizeNow() const noexcept
    {
        float size = 0.35f;

        switch (character)
        {
            // The room you don't hear as an effect: tiny, dark, gone before the next note.
            case Character::ambience:  size = 0.15f; break;
            case Character::room:      size = 0.35f; break;
            case Character::hall:      size = 0.85f; break;
            case Character::plate:     size = 0.60f; break;
            case Character::spring:    size = 0.25f; break;
            // The dressed tail: a big smooth room the chorus will ride.
            case Character::modulated: size = 0.70f; break;
        }

        // DECAY breathes through the room size: half-to-double maps to a quarter of the scale either way,
        // clamped clear of runaway.
        return clampf (0.05f, 0.98f, size + (float) core::det::log2 ((double) decayScale) * 0.25f);
    }

    void refreshTail() noexcept
    {
        const double g    = (double) roomSizeNow() * 0.28 + 0.7;   // Freeverb's room scaling
        const double comb = (1617.0 + 23.0) / 44100.0;               // the right channel's longest comb, s
        // The energy of the k-th pass sits in the interval that ENDS at (k+1) combs, so the last one still
        // above a thousandth runs one comb longer than the envelope's answer.
        tailSec.store ((float) ((double) predelayMs * 0.001
                                  + comb * (1.0 + core::det::log2 (0.001) / core::det::log2 (g))),   // any base
                       std::memory_order_relaxed);
    }

    void apply() noexcept
    {
        FreeverbParams p;
        p.roomSize = roomSizeNow();

        switch (character)
        {
            case Character::ambience:  p.damping = 0.85f; p.width = 0.90f; break;
            case Character::room:      p.damping = 0.50f; p.width = 0.80f; break;
            case Character::hall:      p.damping = 0.30f; p.width = 1.00f; break;
            case Character::plate:     p.damping = 0.15f; p.width = 1.00f; break;
            case Character::spring:    p.damping = 0.70f; p.width = 0.45f; break;
            case Character::modulated: p.damping = 0.25f; p.width = 1.00f; break;
        }

        refreshTail();

        // The wet path is OURS: the room runs fully wet (1/3 undoes its ×3), the dry never enters it, and
        // the ADD at the mix happens in process().
        p.dryLevel   = 0.0f;
        p.wetLevel   = 1.0f / 3.0f;
        p.freezeMode = 0.0f;

        room_.setParams (p);
    }

    void applyHpf() noexcept
    {
        hpf.setParams (eq::FilterType::HighPass, (double) hpfHz, 0.7071, 0.0);
    }

    /** The tail through a slowly breathing delay — depth and rate fixed by taste, not knobs. */
    void modulateTail (int nch, int n) noexcept
    {
        const int   modN   = (int) modLine[0].size();
        const float baseS  = (float) (0.012 * sampleRate);
        const float depthS = (float) (0.003 * sampleRate);
        const float inc    = (float) (0.45 / sampleRate);   // Hz

        for (int i = 0; i < n; ++i)
        {
            lfoPhase += inc;
            if (lfoPhase >= 1.0f)
                lfoPhase -= 1.0f;

            for (int ch = 0; ch < nch; ++ch)
            {
                modLine[(std::size_t) ch][(std::size_t) modPos] = wet[(std::size_t) ch][(std::size_t) i];

                // Quadrature between the channels, so the dress swirls instead of pumping.
                const float ph  = lfoPhase + (ch == 1 ? 0.25f : 0.0f);
                const float lfo = (float) core::det::sin ((double) (ph * kTwoPiF));
                const float dly = baseS + depthS * lfo;

                double pos = (double) modPos - (double) dly;
                while (pos < 0.0) pos += (double) modN;

                const int i0 = (int) pos;
                const int i1 = (i0 + 1) % modN;
                const float f = (float) (pos - (double) i0);
                wet[(std::size_t) ch][(std::size_t) i] = modLine[(std::size_t) ch][(std::size_t) i0] * (1.0f - f)
                                                       + modLine[(std::size_t) ch][(std::size_t) i1] * f;
            }

            modPos = (modPos + 1) % modN;
        }
    }

    static constexpr float kTwoPiF = 6.283185307179586476925286766559f;

    Freeverb room_;
    eq::Svf  hpf;
    core::StateGrid grid_;

    Character character  = Character::room;
    float     mix        = 0.2f;
    float     decayScale = 1.0f;
    std::atomic<float> tailSec { 0.0f };   // see tailSeconds(); set by the constructor, from the formula
    float     predelayMs = 0.0f;
    float     hpfHz      = 120.0f;

    double sampleRate = 48000.0;

    std::array<std::vector<float>, kMaxChannels> wet, wetOut;
    std::array<std::vector<float>, kMaxChannels> preLine, modLine;
    bool  live_[kMaxChannels] { false, false };   // carried by a call since the last clear (law 11a)
    bool  cleared   = true;    // see reset()
    bool  prepared_ = false;
    int   addedCount = 0;      // see addedLength()
    int   prePos = 0, modPos = 0;
    int   maxBlock_ = 0, maxCh_ = 0;
    float lfoPhase = 0.0f;
};

} // namespace felitronics::reverb

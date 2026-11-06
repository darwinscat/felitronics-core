// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/core/BypassFade.h>
#include <felitronics/dynamics/CharacterComp.h>
#include <felitronics/dynamics/Compressor.h>
#include <felitronics/dynamics/GainComputer.h>
#include <felitronics/eq/BassCuts.h>   // CutBiquad — the detector's key

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <vector>

namespace felitronics::basscomp
{

using core::BypassFade;
using dynamics::CharacterComp;
using eq::CutBiquad;

/** THE BASS AMP'S COMP: felitronics-core's compressor, EACH CONTROL ONE JOB. Moved here from OrbitAmp
    (src/core/BassComp.h), its arithmetic unchanged; what it used to write into the plug-in's picture is now a
    read-out of the block (picture()), which the plug-in draws from and a page need not.

      - THE THRESHOLD, dB (-48..0), its own: the line on the picture is where the detector works, in either character.
      - THE RATIO, 1.5:1..10:1, its own.
      - THE KNOB IS BLEND — parallel compression: 0 is no compressor at all (once glided there nothing is computed and
        the signal passes bit for bit), 1 all of the signal compressed, between the compressed signal over the clean one,
        sample-aligned (no lookahead, no latency).
      - MAKE-UP is automatic and LEVEL-NEUTRAL (the owner, 31.10: BLEND changes evenness, not loudness): it gives back
        what the squeeze has TAKEN FROM THIS PLAYING — its average over `makeupAverageMs`, followed only while a note
        sounds and held through the pauses (the owner, 02.11: a threshold still over every note gave make-up from the
        reference note's curve and drove the amp into distortion). A threshold no note reaches gives nothing back; a
        steady note keeps its level at any BLEND. It glides as the squeeze does — up at the attack's pace, down at the
        release's. (`makeupDbFor`, the curve's take from the reference note, stays for the picture and the gates.)
      - Every move GLIDES: the blend over 20 ms for its whole travel, the threshold over 20 ms for its 48 dB, the ratio
        over 20 ms for 1.5:1..10:1 — a drag, a recall or automation is never a step.
      - CHARACTER chooses the detector: SMOOTH hears the level's RMS over 10 ms, PUNCHY its peaks. (The old law raised
        PUNCHY's threshold 4 dB so the two played about as loud; with a threshold of its own the line says where it works,
        and PUNCHY at the same line squeezes harder — the player sets the line; the make-up reads the reference note as
        each detector does, so the two play it as loud.) OPTO, FET and SQUASH are our own (the owner, 01.11 — not copies;
        dynamics::CharacterComp): OPTO fed back on the RMS, its release in two stages with a memory of the squeeze, clean; FET fast
        on the peaks, a grain of odd harmonics growing with the squeeze; SQUASH fed back on the peaks, its attack letting
        the pluck through, its release long, the top let down as it squeezes. Each its own times (attackTable,
        releaseTable); a switch between core's compressor and ours crosses over in characterFadeMs, both fed.

    STANDBY (the strip's arrow) is the house's crossfade (BypassFade, 15 ms), handed in per block. While the link stands
    by with BLEND up the compressor is still FED, so switching it in fades to a compressor already settled. Only with
    BLEND at zero does it stop, its memory cleared then. A BLEND recalled up from zero starts a FRESH compressor, fed and
    not heard for three of its attacks and `warmMs`, then faded in. A PREPARE is not a recall: a chain that arrives with
    BLEND up plays it from its first sample (snapAmount).

    The detector hears the bass through a 30 Hz high-pass, 12 dB/oct (the key): the lowest fundamentals do not pump the
    rest. ATTACK slow 30 / MED 15 / fast 5 ms; RELEASE AUTO — 80 ms after a short peak growing toward 400 ms as the squeeze
    is held — or FAST 80 / SLOW 400 ms. ITS TIME IS ITS OWN: the controls and the release reach the compressor on a fixed
    cadence counted in samples from the prepare — the same playing gives the same samples whatever the host's blocks.

    Mono, as the bass's chain is. Real-time safe: prepare() allocates, process() does not.

    THE PICTURE (picture()): while `pictured` is set and COMP is heard, a block leaves what a picture of it draws — the
    signal it was handed, the squeezed signal before the make-up, the same blended as the sound is, and the squeeze per
    sample — valid until the next process(); the first block of a new run of them says the trace breaks before it.
    Standing by, at BLEND 0, or not pictured, there is none: the picture stands still. Nothing of it reaches the sound. */
class BassComp
{
public:
    enum class Attack    { slow, med, fast };
    enum class Release   { automatic, fast, slow };
    enum class Character { smooth, punchy, opto, fet, squash };
    static constexpr int characterCount = 5;

    static constexpr double kneeDb            = 8.0;
    static constexpr double attackSlowMs      = 30.0,  attackMedMs = 15.0, attackFastMs = 5.0;
    static constexpr double releaseFastMs     = 80.0,  releaseSlowMs = 400.0;
    static constexpr double amountGlideMs     = 20.0;   // the knob's whole travel
    static constexpr double holdMs            = 500.0;  // AUTO's memory of the squeeze
    static constexpr double makeupAverageMs   = 400.0;  // the make-up: the squeeze taken, averaged over this...
    static constexpr double presenceDb        = -50.0;  // ...while a note sounds, its peak over this (dBFS)
    static constexpr double presenceReleaseMs = 150.0;  // how long a note's peak is held for it
    static constexpr double holdFullDb        = 6.0;    // ...held this deep, the release is at its slowest
    static constexpr double cadenceMs         = 1.0;    // how often the knob and the release reach the compressor
    static constexpr double warmMs            = 30.0;   // a fresh start is fed, not heard, three attacks and this long...
    static constexpr double warmFadeMs        = 15.0;   // ...then fades in over BypassFade's length
    static constexpr double characterFadeMs   = 15.0;   // a switch between core's compressor and ours

    /** EACH CHARACTER'S TIMES, ms: ATTACK slow / med / fast, RELEASE fast / slow (AUTO between them as the hold says) —
        OPTO's release is its memory's, its first stage optoFirstMs. */
    static constexpr double attackTable[characterCount][3]  = { { attackSlowMs, attackMedMs, attackFastMs },
                                                                { attackSlowMs, attackMedMs, attackFastMs },
                                                                { 20.0, 10.0, 5.0 }, { 3.0, 1.0, 0.3 }, { 12.0, 6.0, 3.0 } };
    static constexpr double releaseTable[characterCount][2] = { { releaseFastMs, releaseSlowMs }, { releaseFastMs, releaseSlowMs },
                                                                { 400.0, 2000.0 }, { 50.0, 500.0 }, { 250.0, 900.0 } };
    static constexpr double optoFirstMs = 60.0, optoChargeMs = 800.0;
    static constexpr double fetGrain = 0.16, squashGrain = 0.08, squashDarkHz = 3000.0;

    /** Whether the character is ours (dynamics::CharacterComp) rather than core's compressor; fed back; on the RMS. */
    static bool ours     (Character c) noexcept { return c == Character::opto || c == Character::fet || c == Character::squash; }
    static bool fedBack  (Character c) noexcept { return c == Character::opto || c == Character::squash; }
    static bool hearsRms (Character c) noexcept { return c == Character::smooth || c == Character::opto; }

    /** Our character's voice at an ATTACK and a release time (ms, as AUTO has it now). */
    static CharacterComp::Voice voiceOf (Character c, Attack a, double releaseMs) noexcept
    {
        CharacterComp::Voice v;
        v.attackMs  = attackMsOf (c, a);
        v.releaseMs = releaseMs;
        v.rms       = hearsRms (c);
        v.feedback  = fedBack (c);
        if (c == Character::opto)
        {
            v.releaseMs = optoFirstMs;
            v.memoryMs  = releaseMs;
            v.chargeMs  = optoChargeMs;
        }
        else if (c == Character::fet)
            v.grain = fetGrain;
        else if (c == Character::squash)
        {
            v.grain  = squashGrain;
            v.darkHz = squashDarkHz;
        }
        return v;
    }


    /** The key's high-pass: this one, or another for a comparison — the first version's was 40 Hz, two
        sections (24 dB/oct). */
    struct KeyCut
    {
        double hz       = 30.0;
        int    sections = 1;    // 12 dB/oct each
    };


    /** EACH CONTROL ONE JOB (the owner, 31.10): the THRESHOLD is the threshold alone, in dB; the RATIO the ratio alone;
        the knob is BLEND — parallel compression, 0 none at all, 1 all of the signal compressed. Make-up is automatic,
        from the threshold and the ratio. */
    static constexpr double thresholdMinDb = -48.0, thresholdMaxDb = 0.0, defaultThresholdDb = -24.0;
    static constexpr double ratioMin = 1.5, ratioMax = 10.0, defaultRatio = 4.0;

    /** THE MAKE-UP'S REFERENCE: a steady note at this RMS, dBFS, comes out at its own level at any BLEND. */
    static constexpr double referenceDb = -18.0;

    /** HOW EACH DETECTOR READS THE REFERENCE NOTE as its squeeze settles, dB over the note's RMS — measured on core's own
        compressor (steady notes A1..A4 at -18 dBFS RMS, -24 dB, 4:1: the middle of their spread), by
        [CHARACTER smooth, punchy][ATTACK slow, med, fast][RELEASE fast 80 ms, slow 400 ms]. The ballistics sit on the
        squeeze, so a steady note's settles where each cycle's attack meets its release: SMOOTH's RMS reads the note
        itself, within a fifth of a decibel; PUNCHY's peaks a steady tone's crest (3 dB) less what its squeeze lets go
        between them — much with a fast release, little with a fast attack. AUTO's release lies between, by the same
        fraction of its hold that apply() gives its time. The key's 30 Hz high-pass reads B0 and E1 low: those notes come
        up 1.5..3 dB at the full BLEND (the gate says how much). OPTO, FET and SQUASH measured the same way, on our
        compressor (dynamics::CharacterComp), FET's and SQUASH's colour in what they read; OPTO on its own fed-back curve.
        A LIMIT, KNOWN: the table is calibrated at the defaults (-24 dB, 4:1). How far a detector's squeeze lets go between
        peaks scales with how deep it squeezes, so away from them the reference note drifts — at the deep corner (-48 dB,
        10:1) PUNCHY with a fast release comes up about 1 dB at the full BLEND, half that at half; under the knee (-12 dB,
        2:1) there is no make-up and nothing to drift; at the gentlest ratio, 1.5:1, it stays within 0.4 dB at every
        threshold. The gate prints the corners and holds them within 1.2 dB. */
    static constexpr double readsDb[characterCount][3][2] = { { { -0.07, 0.11 }, { 0.01, 0.15 }, { 0.12, 0.17 } },
                                                              { {  0.29, 1.74 }, { 1.05, 2.11 }, { 1.88, 2.47 } },
                                                              { {  0.03, 0.08 }, { 0.09, 0.13 }, { 0.13, 0.14 } },
                                                              { {  1.97, 2.62 }, { 2.43, 2.74 }, { 2.65, 2.80 } },
                                                              { {  1.94, 2.45 }, { 2.27, 2.61 }, { 2.50, 2.71 } } };

    /** ...for a release `slowness` of the way from FAST (0) to SLOW (1). */
    static double readingDb (Character character, Attack attack, double slowness) noexcept
    {
        const auto& r = readsDb[std::clamp ((int) character, 0, characterCount - 1)][std::clamp ((int) attack, 0, 2)];
        return r[0] + (r[1] - r[0]) * std::clamp (slowness, 0.0, 1.0);
    }

    /** The make-up, dB, at a threshold and a ratio, in a character, an attack and a release's slowness: what core's curve
        (the compressor's own, its knee and range) takes from the reference note as that detector reads it — so the
        note at the reference keeps its level. 0 where the reference lies under the knee. */
    static double makeupDbFor (double thresholdDb, double ratio, Character character, Attack attack = Attack::med,
                               double slowness = 1.0) noexcept
    {
        if (fedBack (character))   // the fed-back curve, its knee wider (CharacterComp::staticDeltaDb)
            return -CharacterComp::staticDeltaDb (true, referenceDb + readingDb (character, attack, slowness), thresholdDb,
                                                  ratio, kneeDb);
        felitronics::dynamics::GainComputer curve;
        curve.setMode (felitronics::dynamics::Mode::DownCompress);
        curve.setThresholdDb (thresholdDb);
        curve.setRatio (ratio);
        curve.setKneeDb (kneeDb);
        curve.setRangeDb (felitronics::dynamics::CompressorParams {}.rangeDb);   // the compressor's, as apply() leaves it
        return -curve.deltaDb (referenceDb + readingDb (character, attack, slowness));
    }

    static double attackMsOf (Attack a) noexcept
    {
        return a == Attack::slow ? attackSlowMs : a == Attack::fast ? attackFastMs : attackMedMs;
    }

    static double attackMsOf (Character c, Attack a) noexcept
    {
        return attackTable[std::clamp ((int) c, 0, characterCount - 1)][std::clamp ((int) a, 0, 2)];
    }

    /** The release, ms, a `slowness` of the way from the character's FAST (0) to its SLOW (1). */
    static double releaseMsOf (Character c, double slowness) noexcept
    {
        const auto& r = releaseTable[std::clamp ((int) c, 0, characterCount - 1)];
        return r[0] + (r[1] - r[0]) * std::clamp (slowness, 0.0, 1.0);
    }

    /** The static curve the character's detector works on, dB (0 or below) — the picture's square draws it. */
    static double staticDeltaDb (Character c, double levelDb, double thresholdDb, double ratio) noexcept
    {
        return CharacterComp::staticDeltaDb (fedBack (c), levelDb, thresholdDb, ratio, kneeDb);
    }

    [[nodiscard]] bool prepare (double sampleRate, int maxBlock) { return prepare (sampleRate, maxBlock, KeyCut {}); }

    [[nodiscard]] bool prepare (double sampleRate, int maxBlock, KeyCut cut)
    {
        fs        = sampleRate > 0.0 ? sampleRate : 48000.0;
        block     = std::max (1, maxBlock);
        cadence   = std::max (1, (int) std::lround (fs * cadenceMs * 0.001));
        glideStep = 1.0 / (amountGlideMs * 0.001 * fs);
        ratioStep = (ratioMax - ratioMin) * (double) cadence * glideStep;                      // 1.5:1 to 10:1 in 20 ms
        thrStep   = (thresholdMaxDb - thresholdMinDb) * (double) cadence * glideStep;          // the whole travel in 20 ms
        warmSamples = (3.0 * attackMedMs + warmMs) * 0.001 * fs;   // the attack's own, from the first apply
        warmFade    = warmFadeMs * 0.001 * fs;
        holdK     = 1.0 - std::exp (-1.0 / (holdMs * 0.001 * fs));
        takenK    = 1.0 - std::exp (-1.0 / (makeupAverageMs * 0.001 * fs));
        envDecay  = std::exp (-1.0 / (presenceReleaseMs * 0.001 * fs));
        presence  = std::pow (10.0, presenceDb / 20.0);
        key.assign ((size_t) block, 0.0f);
        dry.assign ((size_t) block, 0.0f);
        work.assign ((size_t) block, 0.0f);
        gr.assign ((size_t) block, 0.0f);
        alt.assign ((size_t) block, 0.0f);
        grAlt.assign ((size_t) block, 0.0f);
        shown.assign ((size_t) block, 0.0f);
        path.assign ((size_t) block, 0.0f);

        keySections = std::clamp (cut.sections, 1, (int) keyCut.size());
        for (auto& b : keyCut)
            b.design (true, cut.hz, fs);

        ready = comp.prepare (fs, block, 1, 0.0);
        ours_.prepare (fs);
        famStep = 1.0 / (characterFadeMs * 0.001 * fs);
        reset();
        return ready;
    }

    /** Silence: no squeeze held, the cadence started again. */
    void reset() noexcept
    {
        comp.reset();
        ours_.reset();
        oursRunning = false;
        famValid    = false;
        for (auto& b : keyCut)
            b.reset();
        heldDb      = 0.0;
        clock       = 0;
        sinceStart  = 0;
        makeupNowDb = 0.0;
        takenDb     = 0.0;
        signalEnv   = 0.0;
        snapMakeup  = false;
        running     = false;
        paramsValid = false;
        tapWrote    = false;
        grDb.store (0.0f, std::memory_order_relaxed);
    }

    /** The knob lands at `amount` without a glide — a prepare's first block. */
    void snapAmount (double amount) noexcept
    {
        amt        = std::clamp (amount, 0.0, 1.0);
        snapMakeup = amt > 0.0;
        if (amt > 0.0)
            sinceStart = 1LL << 40;   // past the warm: what a prepare finds up plays at once
    }

    /** One block of the mono chain, in place. `link` is the strip's STANDBY crossfade for this block; the threshold
        (dB) and the ratio are the hand's, glided to (apply). */
    void process (float* x, int n, double amountTarget, const BypassFade::Span& link,
                  Attack attack, Release release, Character character,
                  double thresholdDb = defaultThresholdDb, double ratio = defaultRatio) noexcept
    {
        shownNow = false;   // this call's picture, if it leaves one (below)
        const double target = std::clamp (amountTarget, 0.0, 1.0);
        thrWant   = std::isfinite (thresholdDb) ? std::clamp (thresholdDb, thresholdMinDb, thresholdMaxDb) : thrWant;
        ratioWant = std::isfinite (ratio) ? std::clamp (ratio, ratioMin, ratioMax) : ratioWant;

        // The knob at zero and glided there: nothing at all — and the memory cleared, now that nothing of it
        // is in the sound. NOTHING WRITTEN for the picture either: a COMP that does nothing stands still (the owner, 31.10).
        if (! (target > 0.0) && ! (amt > 0.0))
        {
            if (running)
                reset();
            grDb.store (0.0f, std::memory_order_relaxed);
            tapWrote = false;
            return;
        }
        if (! ready || x == nullptr || n <= 0 || n > block)
            return;

        running = true;
        const bool heard = link.moving || link.to > 0.0f;   // standing by and settled: fed, not heard

        std::copy (x, x + n, dry.data());
        std::copy (x, x + n, work.data());
        std::copy (x, x + n, alt.data());

        const int   ramp     = std::clamp (link.ramp, 1, n);
        const float linkStep = (link.to - link.from) / (float) ramp;

        for (int i = 0; i < n;)
        {
            // On the cadence, counted from the prepare: the knob as it has glided and the release as AUTO reads
            // the hold, handed to the compressor.
            const int phase = (int) (clock % (long long) cadence);
            if (phase == 0)
                apply (attack, release, character);

            const int len = std::min (n - i, cadence - phase);

            for (int k = i; k < i + len; ++k)
            {
                float v = work[(size_t) k];
                for (int s = 0; s < keySections; ++s)
                    v = keyCut[(size_t) s].process (v);
                key[(size_t) k] = v;
            }

            // OUR CHARACTERS beside core's compressor, both fed while a switch crosses between them; idle, ours takes
            // over the squeeze core's holds, without a step.
            const bool oursWanted = famTarget > 0.0 || fam > 0.0;
            if (oursWanted && ! oursRunning)
                ours_.seed (-comp.gainReductionDb());
            oursRunning = oursWanted;

            float* io[] = { work.data() + i };
            const float* keyIn[] = { key.data() + i };
            (void) comp.process (io, 1, len, keyIn, 1, felitronics::dynamics::GainReductionTap { gr.data() + i, len });

            if (oursRunning)
                ours_.process (alt.data() + i, key.data() + i, len, grAlt.data() + i);

            for (int k = i; k < i + len; ++k)
            {
                if (oursRunning)
                {
                    fam = famTarget > fam ? std::min (famTarget, fam + famStep) : std::max (famTarget, fam - famStep);
                    work[(size_t) k] += (float) fam * (alt[(size_t) k] - work[(size_t) k]);
                    gr[(size_t) k]   += (float) fam * (grAlt[(size_t) k] - gr[(size_t) k]);
                }
                heldDb += holdK * (std::max (0.0, -(double) gr[(size_t) k]) - heldDb);
                amt = target > amt ? std::min (target, amt + glideStep) : std::max (target, amt - glideStep);

                // THE MAKE-UP IS WHAT THE SQUEEZE HAS TAKEN FROM THIS PLAYING: its average over makeupAverageMs, followed
                // only while a note sounds (its peak over presenceDb) and held through the pauses — so a threshold over
                // every note gives nothing back and drives nothing harder into the amp (the owner, 02.11).
                signalEnv = std::max ((double) std::abs (dry[(size_t) k]), signalEnv * envDecay);
                if (signalEnv > presence)
                    takenDb += takenK * (std::max (0.0, -(double) gr[(size_t) k]) - takenDb);
                makeupTargetDb = takenDb;
                makeupNowDb += (makeupTargetDb - makeupNowDb) * (makeupTargetDb > makeupNowDb ? makeupUpK : makeupDownK);
                const float squeezed = work[(size_t) k];   // after the squeeze, before the make-up: the picture's
                work[(size_t) k] *= (float) std::exp (makeupNowDb * 0.11512925464970228);   // ln 10 / 20

                if (! heard)
                    continue;   // standing by: nothing of it in the sound, and the picture is not written (below)

                // BLEND — parallel compression: the compressed signal over the clean one, sample-aligned (no lookahead),
                // in the share the knob has glided to; continuous from its zero by construction. A JUMP to zero is the
                // link's own fade's (the processor fades the link toward "on AND up"), out in 15 ms, before this glide
                // lands and the stage stops.
                const double fresh = std::clamp (((double) (sinceStart + (k - i)) - warmSamples) / warmFade, 0.0, 1.0);
                const double blend = amt * fresh;
                const float  g     = link.moving ? (k < ramp ? link.from + linkStep * (float) k : link.to) : link.to;
                const float  d     = dry[(size_t) k];
                const float  wet   = work[(size_t) k];
                x[k] = d + (wet - d) * (float) (blend * (double) g);
                // THE PICTURE'S OUT: the same blend of the squeezed signal BEFORE the make-up — so the gap under IN is what
                // COMP took (the owner, 31.10: no make-up number anywhere). The sound is x, unchanged.
                shown[(size_t) k] = d + (squeezed - d) * (float) (blend * (double) g);
                path[(size_t) k]  = squeezed;   // AFTER COMP itself: under the input by exactly the squeeze
            }

            i          += len;
            clock      += len;
            sinceStart += len;
        }

        // The indicator shows the squeeze that is heard: none while the link stands by, fed or not.
        grDb.store (heard ? gr[(size_t) n - 1] : 0.0f, std::memory_order_relaxed);

        // THE PICTURE only while COMP is heard (the owner, 31.10): standing by — its strip arrow off, its fade done — it
        // writes nothing, and the picture stands still, dimmed with its block; switched in again, it is alive, on a trace
        // of its own. IN is the signal COMP is handed — the full band, as AFTER and AFTER + BLEND (picture()).
        const bool write = pictured && heard;
        if (write)
        {
            picturedBlock = { dry.data(), path.data(), shown.data(), gr.data(), n, hearsRms (character), ! tapWrote };
            shownNow = true;
        }
        tapWrote = write;
    }

    /** Whether the block's picture is wanted — on screen (the owner sets it per block). */
    bool pictured = false;

    /** WHAT A PICTURE OF THE LAST BLOCK DRAWS: IN, the signal COMP was handed; AFTER, the squeezed signal before the
        make-up; AFTER + BLEND, the same blended as the sound is; the squeeze per sample, dB (0 or below) — `n` samples
        each, `rms` whether the detector heard the RMS, `breaksTrace` whether the trace breaks before it (the block before
        left no picture). */
    struct Picture
    {
        const float* in = nullptr;
        const float* after = nullptr;
        const float* afterBlend = nullptr;
        const float* grDb = nullptr;
        int  n = 0;
        bool rms = false;
        bool breaksTrace = false;
    };

    /** The last process()'s picture — null where it left none (not pictured, standing by, at BLEND 0, refused). Valid
        until the next process(). */
    const Picture* picture() const noexcept { return shownNow ? &picturedBlock : nullptr; }

    /** The ratio the compressor works at now, gliding to or holding the hand's (a gate reads it). */
    double ratioNow() const noexcept { return ratioAt; }

    /** The make-up applied now, dB, as it has glided (a gate reads it). */
    double makeupNow() const noexcept { return makeupNowDb; }

    /** The last block's squeeze, dB per sample (0 or below) — core's own, before the make-up (a gate reads it). */
    const float* squeezeDb() const noexcept { return gr.data(); }

    /** The squeeze now, in dB (0 or below) — for the indicator beside the knob. */
    std::atomic<float> grDb { 0.0f };

    /** Where the knob has glided to — for a gate. */
    double amountNow() const noexcept { return amt; }

private:
    void apply (Attack attack, Release release, Character character) noexcept
    {
        // The threshold and the ratio, each its own control, GLIDING to where they are asked — a drag or a recall is a
        // glide, never a step; a fresh compressor starts on them.
        if (! paramsValid)
        {
            thrAt   = thrWant;
            ratioAt = ratioWant;
        }
        else
        {
            thrAt   = thrWant > thrAt ? std::min (thrWant, thrAt + thrStep) : std::max (thrWant, thrAt - thrStep);
            ratioAt = ratioWant > ratioAt ? std::min (ratioWant, ratioAt + ratioStep) : std::max (ratioWant, ratioAt - ratioStep);
        }

        // The release's slowness, FAST 0 .. SLOW 1 — AUTO's by how deep the squeeze has been held.
        const double slowness  = release == Release::fast ? 0.0 : release == Release::slow ? 1.0 : std::clamp (heldDb / holdFullDb, 0.0, 1.0);
        const double releaseMs = releaseMsOf (character, slowness);

        // Which compressor is heard — core's or ours — a fresh start on it at once, a switch crossing over.
        famTarget = ours (character) ? 1.0 : 0.0;
        if (! famValid)
        {
            fam      = famTarget;
            famValid = true;
        }
        if (ours (character))
            ours_.set (voiceOf (character, attack, releaseMs), thrAt, ratioAt, kneeDb);

        felitronics::dynamics::CompressorParams p;
        p.thresholdDb = thrAt;   // the line is where the detector works, in either character (no offset any more)

        // The paces the make-up glides at; its number is what the squeeze takes from the playing (process).
        makeupUpK      = 1.0 - std::exp (-1.0 / (attackMsOf (character, attack) * 0.001 * fs));
        makeupDownK    = 1.0 - std::exp (-1.0 / (releaseMs * 0.001 * fs));
        warmSamples    = (3.0 * attackMsOf (character, attack) + warmMs) * 0.001 * fs;
        if (snapMakeup)
        {
            makeupNowDb = takenDb;
            snapMakeup  = false;
        }
        p.ratio       = ratioAt;
        p.kneeDb      = kneeDb;
        p.attackMs    = attackMsOf (character, attack);
        p.releaseMs   = releaseMs;
        p.detector    = hearsRms (character) ? felitronics::dynamics::Detector::Rms : felitronics::dynamics::Detector::Peak;
        p.rmsWindowMs = 10.0;
        p.autoMakeup  = false;   // the make-up is ours (above)
        p.makeupDb    = 0.0;
        p.lookaheadMs = 0.0;

        if (! paramsValid || ! same (p, set))
        {
            comp.setParams (p);
            set = p;
            paramsValid = true;
        }
    }

    static bool same (const felitronics::dynamics::CompressorParams& a,
                      const felitronics::dynamics::CompressorParams& b) noexcept
    {
        const auto near = [] (double u, double v, double within) { return std::abs (u - v) < within; };
        return near (a.thresholdDb, b.thresholdDb, 1.0e-9) && near (a.ratio, b.ratio, 1.0e-9)
            && near (a.attackMs, b.attackMs, 1.0e-9) && near (a.releaseMs, b.releaseMs, 1.0) && a.detector == b.detector;
    }

    felitronics::dynamics::Compressor       comp;
    felitronics::dynamics::CompressorParams set;
    CharacterComp ours_;
    double fam = 0.0, famTarget = 0.0, famStep = 0.001;   // 0 core's compressor heard, 1 ours
    bool   famValid = false, oursRunning = false;
    std::array<CutBiquad, 2> keyCut;
    int    keySections = 1;
    std::vector<float> key, dry, work, alt, gr, grAlt, shown, path;
    double fs = 48000.0, heldDb = 0.0, amt = 0.0, glideStep = 0.001, holdK = 0.001;
    double thrWant = defaultThresholdDb, thrAt = defaultThresholdDb, thrStep = 0.05;
    double ratioWant = defaultRatio, ratioAt = defaultRatio, ratioStep = 0.4;
    long long clock = 0, sinceStart = 0;
    double warmSamples = 960.0, warmFade = 720.0;
    double makeupNowDb = 0.0, makeupTargetDb = 0.0, makeupUpK = 0.001, makeupDownK = 0.001;
    double takenDb = 0.0, takenK = 0.0001, signalEnv = 0.0, envDecay = 0.999, presence = 0.003;   // the make-up's measure
    bool   snapMakeup = false;
    int    block = 0, cadence = 48;
    bool   ready = false, running = false, paramsValid = false;
    bool   tapWrote = false;   // a picture left last block: not, and the next one breaks its trace
    bool   shownNow = false;   // this call left one (picture())
    Picture picturedBlock;
};

} // namespace felitronics::basscomp

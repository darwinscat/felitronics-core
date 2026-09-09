// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.
#pragma once
// felitronics::rigplayer — ONE DEVICE OF A PACK, PLAYING. Hand it the pack's structures
// (namz::rig — the manifest as the canonical reader returns it) and a way to get a file's bytes, tell it
// where the knobs stand, feed it audio. It does, in signal order:
//
//     in → pack in → dry copy → pre tone (bands, curve) → chain trim → ┬→ trim A → model A → align ─┐
//                                                                      └→ trim B → model B → align ─┴→ mix
//        → post tone (bands, curve) → dry/wet blend → pack out → out
//
//   • WHICH files play is namz::rig's policy (a turned control is law) joined to ModelBlend.h's pair
//     along the gain dial — RigSelection.h, pure.
//   • WHEN a model may be replaced, and how loudly each is heard, is felitronics::nam::BlendLaw —
//     run here once per block on the audio thread, the only writer of the weight. The message thread
//     posts a request and answers load asks; two one-way mailboxes, no locks. A load itself — bytes
//     fetched, a network built and warmed, some twenty milliseconds for a WaveNet — is a JOB the
//     player hands out (takeLoadJob) and the host runs where it likes (run, any thread) before
//     bringing it back (deliver): no thread of its own, no policy, no hiccup on the drawing thread.
//   • A slot that has stood silent under an unchanged request for kColdAfterSeconds goes COLD and its
//     model stops being run on audio — a dial parked on a capture costs one network, not two. It is
//     still HANDED the call, at width zero, until the stage has drained the silence it owes (law 11c);
//     that is one bounded drain per sleep and nothing after it, and without it the slot replayed what
//     it was holding when it fell asleep — 0.500000 out of digital silence. The law owns the flag
//     (BlendState::cold) and wakes the slot, warm-up first, the moment the request changes; the model
//     stays where it is throughout. slotCold() and coldBlocks() read it out.
//   • A tone knob is applied as the pack describes it — a curve becomes a minimum-phase FIR, bands
//     become biquads — on the side of the models its circuit puts it. ToneKnobs.h, pure.
//   • A blend knob mixes the DI the models were fed, through the dry path's own response. BlendKnob.h.
//   • THE PACK'S OWN TWO LEVELS BRACKET EVERYTHING IT DOES. `chain[].input_db` is the level of the
//     guitar into this device, and it goes FIRST — ahead of the dry copy, because a blend knob mixes
//     one guitar with itself and both ends must be fed the same signal. `chain[].output_db` is the
//     level the device leaves at, and it goes LAST, after the mix, so balancing one pack against
//     another cannot move the blend the pack was told to hold. Neither is behind setInputTrims():
//     they are the author's statement about the pack, not a listener's option.
//   • The slot trims are the files' own `input_db` — a linked setting plays its neighbour softer, and
//     this is the one thing setInputTrims() switches off — and the chain trim is the extension past
//     the ends of the captured range, which the dry path never saw either.
//   • The models' offsets are measured from the bytes (AlignmentTable.h) by the host, on any thread,
//     and handed in; each slot's delay travels with its model and lands when the slot is silent.
//
// NOTHING OF THE LIBRARY IS KNOWN HERE — no take, no database, no file system, no JUCE. The bytes come
// through ModelSource by `files[].id`, when the law asks and never before: a swept dial fetches each
// capture once. This is the unit that leaves for the plugin as one directory.
//
// Threads, exactly as NamStage: prepare() / load() / the knob setters / service() / takeLoadJob() /
// deliver() on the message thread; run() on any thread but the audio one; process() alone on the audio
// thread — it allocates nothing, locks nothing, touches no file. The read-outs are atomics.

#include <felitronics/rigplayer/AlignmentTable.h>
#include <felitronics/rigplayer/BlendKnob.h>
#include <felitronics/rigplayer/RigSelection.h>
#include <felitronics/rigplayer/ToneKnobs.h>

#include <felitronics/convolution/CabConvolver.h>
#include <felitronics/core/DryAligner.h>
#include <felitronics/core/StreamResampler.h>
#include <felitronics/core/StateGrid.h>
#include <felitronics/eq/MatchedBiquad.h>
#include <felitronics/lineareq/MagnitudeCurve.h>
#include <felitronics/nam/BlendLaw.h>
#include <felitronics/nam/NamStage.h>
#include <namz_rig.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace felitronics::rigplayer {

class RigPlayer {
public:
    static constexpr int kMaxChannels = 2;
    // Biquads per side. Every band of every `sections` knob on that side is one; a pack carrying more
    // than this on one side would be a first, and the rest are left out (bandsDropped() says so).
    static constexpr int kMaxBands = 24;

    // HOW LONG A COEFFICIENT MOVE TAKES, and it is READ OFF THE BAND rather than chosen or inherited.
    //
    // It used to be the host's block: the coefficients travelled linearly across whatever `process()`
    // was handed, so the same turn of the same knob took 0.33 ms at a block of 16 and 85 ms at 4096 —
    // 180 SECONDS on the whole-file call an offline driver makes. A spread of 500000x is not a
    // behaviour, so there was nothing to preserve; what there was instead was a question with no
    // number in the file, and this is that number.
    //
    // WHAT SETS IT. A biquad cannot be moved faster than it forgets: its poles sit at radius
    // r = sqrt(|a2|), the impulse response decays as r^n, and a move completed well inside that decay
    // asks the filter to be somewhere its own state has not caught up to. So ONE ring time — the
    // filter's own -60 dB memory — is the length, and the two rails below are what keep that honest.
    //
    // THERE IS NO KNEE, AND THE FIRST VERSION OF THIS NOTE CLAIMED ONE. The transient a fast move adds
    // over the ideal is +14.88 dB at L = 1, +3.13 at 64, +0.93 at 256, +0.50 at 512, +0.13 at 2048 —
    // which LOOKS like it flattens, and does not: as a linear excess it is L x excess = 27.8, 28.9,
    // 30.3, 30.9, i.e. a plain 1/L all the way down. So no length is "where the artefact stops"; every
    // doubling halves it, for ever, and the number is a CHOSEN TOLERANCE dressed as physics if it is
    // stated any other way. The adiabatic error behaves the same, ~3 dB per doubling with no knee
    // either.
    //
    // WHAT THE MULTIPLIER IS, then: a responsiveness policy, and ONE ring is defended UNIVERSALLY rather
    // than from the only pack that exists — appealing to that pack is over-fitting, and an earlier
    // version of this note did exactly that on a reading taken at a BLIND PROBE FREQUENCY. It said the
    // pack's bass shelf "overshoots 0.00 dB at every length"; that was measured at 90 Hz and 1400 Hz,
    // both outside the shelf's transition band. Inside it the same move reads +3.04 dB at L = 1, +1.41
    // at 64 and +0.20 at one ring (436 samples) on a 334 Hz probe, and the OTHER direction (+6 -> -12)
    // reads +17.19 dB at 42 Hz on an instantaneous jump. The pack's own knob was never free.
    //
    // The universal statement, measured across every shape tried (Q 1.4 bell, the pack shelf, a 400 Hz
    // Q 0.3 bell, a 300 Hz Q 0.5, a 1 kHz Q 0.7 at -24 dB, a 100 Hz Q 8, a Q 2, an 80 Hz Q 0.3 shelf):
    // AT EXACTLY ONE RING the peak excess over the ideal move on a FULL-SWING SLAM is +0.24 to +1.06 dB.
    // That is a property of the ramp-to-ring ratio, not of a pack, because the excess falls as 1/L.
    // Five rings buys 5x less (+0.05 to +0.25) for 5x the lag. One ring also lands the pack's bass knob
    // at 9.1 ms, which is what a 512-sample host was already producing — a tie-breaker, not the reason.
    //
    // ⚠️ AND IT DOES NOT COVER THE PACK'S SWITCH, whatever an earlier version of this note said. That
    // switch's two positions carry ZERO and ONE section, so flipping it changes the BAND COUNT, and
    // `takeBands` starts an appearing band AT its coefficients with reset state — measured through the
    // player, the output matches a static "smooth" render from ELEVEN samples after the flip, where a
    // 64-sample ramp would need at least 64. An appearance is a hard step today. Ramping it from
    // identity (a stable endpoint, and the stability triangle is convex so the linear path stays
    // stable) and draining a disappearing band back to identity is a separate policy — P26 follow-up.
    //
    // THE FLOOR IS A RAIL, NOT A GUARANTEE, and the difference is the part that is easy to overclaim.
    // Every band in the reference pack (namz conformance `rig-measured`: a 167 Hz Q 0.58 low shelf
    // travelling [-12, +6] dB with its corner over [120, 240] Hz, a 1077 Hz Q 0.65 tilt, a 10 kHz
    // Q 0.8 switch position) is a shelf or tilt with Q BELOW ONE, and on those the adiabatic error is
    // already -47 dB at an INSTANTANEOUS jump. But 64 samples does NOT bound the artefact in general:
    // a schema-legal 5 kHz bell at Q 1.4 swept +12 -> -12 dB is handed exactly 64 by the rule above and
    // overshoots by 0.67 dB there (1.49 dB at the worst edit phase), falling as 1/L to 0.11 dB at 512.
    // 64 is where the lag budget starts, chosen so an ordinary knob is not made sluggish; a move large
    // and sharp enough will exceed it, and that is stated rather than hidden.
    //
    // THE CEILING IS A LAG BUDGET, not a DSP one, which is why it is a named constant and not a
    // formula: a band that rings 176 ms (a 100 Hz bell at Q 8 — the pack schema permits it, nothing
    // in `namz_rig_load.h` bounds `hz` or `q` beyond "positive") would ask for 880 ms, and a knob that
    // trails the hand by a second is broken in a different way. Past the ceiling the move is honestly
    // too fast for the filter and the overshoot above is what it costs.
    // THE RAILS ARE TIME, and the first version had them in SAMPLES — which is the same defect as the
    // one this whole change removes, one scale down. A ring is a duration: a 167 Hz shelf rings 9.08 ms
    // at every rate. Rails in samples do not: 2048 is 46.4 ms at 44.1 kHz, 42.7 at 48, **21.3 at 96 and
    // 10.7 at 192**, so the same pack gets a quarter of the budget on a 192 kHz session, and a 100 Hz
    // Q 2 bell swept +-12 dB costs +1.5 dB of overshoot at 42.7 ms against **+4.7 dB at 10.7**.
    static constexpr double kBandRampMinSeconds =   64.0 / 48000.0;   // 1.33 ms — the responsiveness rail
    static constexpr double kBandRampMaxSeconds = 2048.0 / 48000.0;   // 42.7 ms — the lag budget, the only taste in here
    static constexpr double kBandRampRings = 1.0;   // ONE ring; see the note above for why not five

    // `sampleRate` has no default on purpose: a caller that forgets it would silently get the 48 kHz
    // rails at every rate, which is the defect the note above describes.
    static int bandRampLength(const felitronics::rigplayer::SectionBiquad& from,
                              const felitronics::rigplayer::SectionBiquad& to,
                              double sampleRate) {
        const double fs   = sampleRate > 0.0 ? sampleRate : 48000.0;
        const int kBandRampMin = std::max(1,  (int) (kBandRampMinSeconds * fs));
        const int kBandRampMax = std::max(kBandRampMin, (int) (kBandRampMaxSeconds * fs));
        // THE DOMINANT ROOT, and `sqrt(|a2|)` is NOT it. a2 is the pole PRODUCT, so its square root is
        // the GEOMETRIC MEAN of the two magnitudes — right only when they are equal, which is exactly
        // the complex-conjugate case and no other. A real pair is what a LOW-Q design gives, and packs
        // are full of low-Q shelves. `matched::highShelfQDb(1 kHz, 48 kHz, +12 dB, Q 0.05)` has real
        // poles at 0.99346 and 0.07343 — BOTH POSITIVE, and it is the FAST one that the mean drags the
        // answer down to: the mean says 5.28 samples where the slow pole rings for 1052.78, a 199x
        // underestimate. Its -12 dB sibling (poles 0.99670 and 0.39809) misses by 140x: 14.95 against
        // 2090.39. Two shapes, two different factors — quoting one shape's pair with the other's factor
        // is how the first draft of this note was wrong.
        auto ring = [kBandRampMax] (double a1, double a2) {
            const double disc = a1 * a1 - 4.0 * a2;
            const double rho  = disc >= 0.0 ? (std::fabs(a1) + std::sqrt(disc)) * 0.5   // real pair: the larger root
                                            : std::sqrt(std::fabs(a2));                 // complex pair: both are sqrt(a2)
            if (! (rho > 0.0)) return 0.0;                     // no pole at all: nothing to outrun
            if (! (rho < 1.0)) return (double) kBandRampMax;    // marginal or worse: take the ceiling
            return std::log(1.0e-3) / std::log(rho);            // samples to -60 dB of the SLOWEST mode
        };
        const double n = kBandRampRings * std::max(ring(from.a1, from.a2), ring(to.a1, to.a2));
        if (! (n > (double) kBandRampMin)) return kBandRampMin;
        if (n >= (double) kBandRampMax)    return kBandRampMax;
        return (int) n;
    }

    // The curve form's FIR: 1024 taps designed through 8192 points, the numbers the bench auditions with.
    static constexpr int kFirTaps = 1024, kFirDesign = 8192;
    static constexpr int kMaxDelay = felitronics::nam::kBlendMaxDelay;
    // How long a slot may stand at exactly zero, its request unchanged, before its model stops being
    // run (BlendLaw.h, "and one economy"). Long against a hand — a pause this long IS the dial at rest
    // — and short against a session, so the saving is there whenever nobody is turning. A host may set
    // its own (setColdAfterSeconds); zero or less = never.
    static constexpr double kColdAfterSeconds = 2.0;

    // 🔴 THE CEILING ON A HOST RATE, and every buffer derived from one is bounded by it. Far above any
    // audio rate — 3 MHz is sixteen times the highest a DAW offers — and its job is not to judge taste
    // but to keep `(int) f(sampleRate)` inside int and the rings it sizes inside memory. A rate above
    // it falls back to 48 kHz the same way a negative or non-finite one does.
    //
    // The VALUE is the house convention, not a new number: dynamics::Compressor::kMaxSampleRate and
    // limiter::TruePeakLimiter::kMaxSampleRate are both 3.0e6 with the same stated reason, and
    // eq::EqBand::prepare refuses past the same figure as a literal. That it now has a FOURTH spelling
    // is the very disease this branch exists to treat, and consolidating the four is its own task —
    // rigplayer links neither dynamics nor limiter, so it cannot simply ask.
    static constexpr double kMaxSampleRate = 3.0e6;

    // The ONE place a host rate is judged, so that nothing downstream re-decides it. A rate outside
    // the accepted range becomes the factory rate rather than a refusal, which is what this class did
    // for a non-positive rate before and what stereo::MonoBass and the dynamics detectors do too.
    //
    // Spelled as a RANGE and positively — the idiom TruePeakLimiter documents: `> 0.0` is false for a
    // NaN and `<= kMaxSampleRate` is false for an infinity, so between them they exclude everything
    // std::isfinite would have, and an isfinite here would be a clause no test could distinguish.
    static double usableSampleRate (double hostSR) noexcept
    {
        return (hostSR > 0.0 && hostSR <= kMaxSampleRate) ? hostSR : 48000.0;
    }

    // 🔴 THE DRY ALIGNER'S CAPACITY — ASKED, NOT DERIVED, and the whole of it is now one question
    // put to the module that owns the answer. Pure and public for the reason rateMatch() is: it is the
    // only way to pin the number at rates this repository does not itself run.
    //
    // 🔴 WHAT USED TO BE HERE AND WHY IT IS GONE. This was
    // `max(256, ceil(pairDelayHostSamples(fs, kModelSampleRate)) + 2)`, and both halves of that
    // existed for ONE reason, named in its own comment: NamStage::install() accepted a model within
    // half a hertz of the CURRENT run rate and prepare() then adopted it, so repeated loads walked the
    // run rate away from the factory value without bound. A lower run rate is a longer round trip, so a
    // capacity derived from 48 kHz could come up a slot short of a delay the stage really reports —
    // and DryAligner clamps silently. The floor and the spare slot were a MITIGATION with a measured
    // price, not a proof: the first silent clamp arrived after 82256 half-hertz steps at a 48 kHz host,
    // 68511 at 96 kHz, 41021 at 192 kHz, 560 at 384 kHz and 72 at a 3 MHz one.
    //
    // P38 closed the walk at its source — the gate's reference is the constant now, so an accepted
    // model's rate lies in [kModelSampleRate - kModelRateTolerance, kModelSampleRate +
    // kModelRateTolerance] for the life of the process — and a mitigation outlives its defect only as
    // a number nobody can re-derive. So it goes, together with the fix rather than after it.
    //
    // What replaces it is not a smaller guess but the contract's own consequence:
    // NamStage::maxLatencySamples(fs) is an upper bound on what any ACCEPTED model can report at this
    // host, +1 because DryAligner's usable range is capacity-1 (it clamps to [0, capacity-1], and does
    // it in silence — which is exactly how a downstream repository shipped a bypass path that
    // under-delayed every host rate above 48 kHz). One reason for the "+1", none for anything else, and
    // no copy of 48000 or of the half hertz on this side of the seam. The bound's floating-point
    // environment precondition applies here too: capacity prepared to nearest at
    // nextafter(96748.9921875, 0) carries 96, but a model prepared upward asks for 97 and clamps.
    // See NamStage::maxLatencySamples and review plan 2; this function supplies no extra margin.
    static int dryAlignerCapacity (double hostSR) noexcept
    {
        return felitronics::nam::NamStage::maxLatencySamples (usableSampleRate (hostSR)) + 1;
    }

    RigPlayer() = default;
    RigPlayer(const RigPlayer&) = delete;
    RigPlayer& operator=(const RigPlayer&) = delete;

    // ---------------------------------------------------------------- message thread: setup ----

    // Sizes every buffer and (re)builds anything designed for a rate. Never while process() runs.
    // LAW 11(b): BINDING, and it used to CLAMP — `prepare(48000, 256, 4)` was accepted as a 2-channel
    // player, after which every `process(io, 4, n)` was refused forever and the caller was told only at
    // the second call. That is the defect the law's own text describes, in this file.
    [[nodiscard]] bool prepare(double sampleRate, int maxBlock, int numChannels) {
        prepared_ = false;                       // law 11(b): disarm FIRST, then validate, then write
        if (numChannels < 1 || numChannels > kMaxChannels) return false;
        if (maxBlock < 1) return false;
        // 🔴 IN RANGE, not merely finite — and the difference is the whole point, because an earlier
        // version of this line said `isfinite` and that is NOT the property the code below depends on.
        // The dry-aligner capacity is `(int) ceil(<a function of fs_>)`, and an out-of-range float→int
        // conversion is undefined: `isfinite` lets 1e300 through, and 1e300 converts just as badly as
        // an infinity does. Measured through this very function with UBSan, before this line was
        // widened: prepare(1e300, 64, 2) returned TRUE while firing twice — "6.66667e+296 is outside
        // the range of representable values of type 'int'", then "signed integer overflow: 2147483647
        // + 2". On the base commit, where the capacity was a literal, neither fired: the exposure came
        // in with the computation, and half of it survived the fix that was written for it.
        //
        // The second face of the same line is not undefined at all and is worse for being legal: the
        // capacity used to be a CONSTANT and is now a function of the argument, so the ring it sizes
        // grew without an upper bound. Measured on the same probe — heap requested inside one
        // prepare() call: 0 bytes at 1e11 on the base commit, 533 333 608 (508.6 MiB) here.
        //
        // kMaxSampleRate closes both, and it is the house number rather than a new one. Spelled as a
        // RANGE and positively, which is the same idiom TruePeakLimiter documents: `sampleRate > 0.0`
        // is false for a NaN and `sampleRate <= kMaxSampleRate` is false for an infinity, so the two
        // comparisons already exclude everything `std::isfinite` would have — an explicit isfinite
        // here would be a clause no test could ever distinguish.
        fs_       = usableSampleRate (sampleRate);
        maxBlock_ = maxBlock;
        channels_ = numChannels;   // validated above — law 11(b) forbids the clamp that was here
        coldAfter_.store(coldSamples(coldSeconds_), std::memory_order_release);
        for (auto& n : nam_) n.prepare(fs_, maxBlock_);
        // Not normalised: a tone curve's broadband level is part of what the pack says, and a dry path's
        // level rides in its gain. Reference-unity RMS is for cabinets, which these are not.
        prepared_ = false;                       // a refused sub-prepare leaves the player unprepared
        for (auto& f : fir_) if (! f.prepare(fs_, maxBlock_, channels_, 0.1, false)) return false;
        if (! dry_.prepare(fs_, maxBlock_, channels_, 0.1, false)) return false;
        // 🔴 THE DRY PATH HAS TO BE DELAYED BY WHAT THE WET PATH COSTS, or the blend is a COMB FILTER.
        // The mix in process() sums `a` (through the models, hence through their rate-match) with `d`
        // (the DI, through the dry FIR only). Those two are misaligned by exactly latencySamples(), and
        // summing them notches at fs/(2·D) with every odd multiple above it.
        //
        // This was ALREADY wrong before P34 and nobody had put a number on it: with the cubic's 3.84
        // samples the first null sat at 5742 Hz. P34's 61.4 samples move it to 359 Hz — the body of a
        // guitar, not a phasey top — and the host cannot fix it, because latencySamples() reports the
        // whole player's PDC outward while this notch is INTERNAL to the blend. Measured on a 50/50
        // blend, worst dip across 100 Hz … 10 kHz: -50 dB unaligned against 0.00 dB aligned.
        //
        // Capacity: it must EXCEED the delay it will ever hold — DryAligner clamps to [0, capacity-1]
        // and does it SILENTLY, which is exactly how a downstream repository shipped a bypass path
        // that under-delayed every host rate above 48 kHz.
        //
        // 🔴 THE PREVIOUS VERSION OF THESE LINES WAS THE SAME MISTAKE, WRITTEN BY THE FIX FOR IT. It
        // said "the largest this can report is a 192 kHz host against a 48 kHz model — D·(1 + 4) = 160
        // — so 256 covers it with margin", and that sentence is a RESTATEMENT of the round-trip
        // formula sitting in front of a silent clamp. It was true when written and would have stopped
        // being true the moment the kernel length became a function of the ratio. Ask instead — and
        // the asking now lives in dryAlignerCapacity(), where it can be pinned at rates this file
        // never runs, because inline arithmetic here is exactly what could not be.
        dryLatency_.prepare(channels_, maxBlock_, dryAlignerCapacity (fs_));
        for (int c = 0; c < kMaxChannels; ++c) {
            slotB_[c].assign((std::size_t) maxBlock_, 0.0f);
            dryBuf_[c].assign((std::size_t) maxBlock_, 0.0f);
            spare_[c].assign((std::size_t) maxBlock_, 0.0f);
            for (auto& t : lagTail_) t[c].fill(0.0f);
        }
        for (auto& side : bq_) for (auto& band : side) for (auto& b : band) b.reset();
        for (auto& g : bandGrid_) g.reset();             // a stream restart re-anchors the audio-time grid
        // ...AND RETIRES THE BANDS THEMSELVES, which cancelling the ramp alone does not. `bandTo_` and
        // `bandCur_` survive a prepare() otherwise, `rebuildBands` republishes the same COUNT, and
        // `takeBands` re-initialises only the bands beyond that count — so the first block after a
        // restart GLIDES from the previous stream's coefficients to the new ones. At a new sample rate
        // that is worse than stale: a 167 Hz shelf's 48 kHz coefficients ARE a 334 Hz shelf at 96 kHz,
        // and the band slides down from there over a ring. Measured on a dial published before
        // prepare(): 5951 samples differ, worst 0.54. Zeroing the count makes every band re-appear AT
        // its coefficients on the first block, which is what a restart means. prepare() is
        // contractually never concurrent with process().
        for (int s = 0; s < 2; ++s) {
            bandRt_[s].count = 0;
            for (int k = 0; k < kMaxBands; ++k) bandPos_[s][k] = bandLen_[s][k] = 0;
        }
        curIn_    = inGain_.load(std::memory_order_relaxed);
        curOut_   = outGain_.load(std::memory_order_relaxed);
        curChain_ = chainGain_.load(std::memory_order_relaxed);
        for (int i = 0; i < 2; ++i) curSlot_[i] = slotGain_[i].load(std::memory_order_relaxed);
        curDry_ = dryGain_.load(std::memory_order_relaxed);
        curWet_ = wetGain_.load(std::memory_order_relaxed);
        prepared_ = true;
        if (loaded_) {                                   // the FIRs and the bands were designed for the old rate
            for (int s = 0; s < 2; ++s) { rebuildCurves(s); rebuildBands(s); }
            rebuildDry();
        }
        return true;
    }
    bool   prepared()   const { return prepared_; }
    double sampleRate() const { return fs_; }
    int    channels()   const { return channels_; }

    // The first stage of the chain this player can run. False = nothing to play (no NAM stage, or one
    // with no files) and the player is left empty.
    bool load(const namz::rig::Rig& rig, ModelSource source) {
        for (const auto& st : rig.chain)
            if (st.kind == namz::rig::StageKind::Nam) return load(st, std::move(source));
        unload();
        return false;
    }

    bool load(const namz::rig::Stage& stage, ModelSource source) {
        unload();
        if (stage.kind != namz::rig::StageKind::Nam || stage.device.files.empty()) return false;
        stage_  = stage;
        publishLevels();          // the new pack's levels stand from here, not from the end of load()
        source_ = std::move(source);
        const auto& d = stage_.device;
        // One model per DISTINCT file: several settings pointing at one file — the pack's link — share
        // its bytes, its id and its delay; only their trims differ, and those belong to the slot.
        fileModel_.assign(d.files.size(), -1);
        for (std::size_t i = 0; i < d.files.size(); ++i) {
            const auto& id = d.files[i].id;
            int m = -1;
            for (std::size_t k = 0; k < models_.size(); ++k) if (models_[k].id == id) { m = (int) k; break; }
            if (m < 0) { models_.push_back({ id, {} }); m = (int) models_.size() - 1; }
            fileModel_[i] = m;
        }
        // Where every knob starts: the pack says, the player never invents (namz::rig::defaultSettings,
        // and `default`-else-`reference` for the linear knobs).
        axes_ = namz::rig::defaultSettings(d);
        if (const auto* dial = crossfadeDial(d)) {
            dial_ = dial->name;
            dialSweep_ = dial->sweep;
            const int deg = degreesOf(axes_[dial_]);
            dialDeg_ = deg >= 0 ? (double) deg : 0.0;
        }
        // …and the defaults are per CONTROL, not a promise that the combination was captured. If no
        // file stands there, pin the first control at its default and let resolve() find the closest
        // captured combination — the pack opens on something sounding, never on a silent panel.
        if (select(d, axes_, dial_, dialDeg_, shape_, topExtend_).fileA < 0)
            for (const auto& c : d.controls)
                if (namz::rig::resolve(d, axes_, c.name, axes_[c.name]) != nullptr) break;
        for (const auto& t : tones())      toneAt_[t.name]  = toneStart(t);
        for (const auto& b : stage_.blend) blendAt_[b.name] = blendStart(b);
        // The models' offsets, from the pack when it carries them — measured once, when it was built.
        // A pack without them leaves the table empty; the host may measure and setAlignment().
        align_ = AlignmentTable::fromDevice(d);
        loaded_ = true;
        if (prepared_) {
            for (int s = 0; s < 2; ++s) { rebuildCurves(s); rebuildBands(s); }
            rebuildDry();
        }
        apply();
        return true;
    }

    // Forget the device: both models cleared, the law told to start over, every filter a wire.
    void unload() {
        loaded_ = false;
        ++gen_;                                          // a job out for this pack comes back to be dropped
        setRequest(0, 0, 0.0f);
        // A landing (or a failure) published and not yet consumed dies with its pack — and it must die
        // BEFORE the forget is posted: consumed after it, it would land a stale model in a law that was
        // just wiped and a stage that was just emptied, and an empty stage carried is the raw DI.
        landFlag_.store(false, std::memory_order_release);
        landFail_.store(0, std::memory_order_release);
        forget_.store(true, std::memory_order_release);
        loadAsk_.store(0, std::memory_order_release);    // …and an ask for it is void with it
        for (auto& n : nam_) n.clearModel();
        models_.clear(); fileModel_.clear();
        axes_.clear(); dial_.clear(); dialSweep_ = 0; dialDeg_ = 0.0;
        toneAt_.clear(); blendAt_.clear(); toneOverride_.clear();
        sel_ = {}; plan_ = {};
        stage_ = {};
        source_ = nullptr;
        align_ = {};
        for (int s = 0; s < 2; ++s) {
            curveSum_[s].clear();
            firBypass_[s].store(true, std::memory_order_release);
            publishBands(s, BandSet {});
        }
        bandsDropped_ = 0;
        dryActive_.store(false, std::memory_order_release);
        dryFirBypass_.store(true, std::memory_order_release);
        dryGain_.store(0.0f, std::memory_order_release);
        wetGain_.store(1.0f, std::memory_order_release);
        // …and with no stage there is nothing to be fed or delivered: an empty player is a wire,
        // whatever hand the host is holding. The hand itself is kept — it is the bench's, not the
        // pack's — and lands again on the next load().
        inGain_.store(1.0f, std::memory_order_release);
        outGain_.store(1.0f, std::memory_order_release);
        chainGain_.store(1.0f, std::memory_order_release);
        for (auto& g : slotGain_) g.store(1.0f, std::memory_order_release);
        for (auto& p : pendDelay_) p.store(0, std::memory_order_release);
        for (auto& c : slotCold_) c.store(false, std::memory_order_release);
    }
    bool loaded() const { return loaded_; }
    const namz::rig::Stage& stage() const { return stage_; }

    /// THE HOST'S OWN HAND ON THE PACK'S TWO LEVELS — the WHOLE number it wants this device played
    /// at, never a difference from what the pack states. Absent (the default) leaves the pack's own
    /// number standing, which is what a plugin wants and why a plugin needs none of this.
    ///
    /// A host with a fader used to subtract: it sent `hand - what the pack says`, and applied that
    /// outside the player. That cannot be made correct. The subtraction needs to know which pack is
    /// loaded RIGHT NOW, and a host reads that from its own document — which changes when somebody
    /// edits it, when a rebuild is in flight, when a device is switched mid-build. Every one of those
    /// left the level wrong, silently, in the same class as the double application this pair of keys
    /// was added to end. The player is the one applying the level and the only place that cannot
    /// disagree with itself about which pack it holds, so the subtraction lives here or nowhere.
    ///
    /// Message thread, like every other setter. Survives a load: the hand belongs to the bench and
    /// not to the pack, so the next pack is played at the same number without the host restating it.
    void setHostInputDb (std::optional<double> db) { hostIn_  = db; publishLevels(); }
    void setHostOutputDb(std::optional<double> db) { hostOut_ = db; publishLevels(); }
    std::optional<double> hostInputDb()  const { return hostIn_; }
    std::optional<double> hostOutputDb() const { return hostOut_; }
    // The pack's own two levels, for a host that prints or draws them. These are what the player IS
    // applying: a host with a fader of its own applies only its deviation from them, or the level
    // lands twice — which is exactly the bug that made these keys necessary.
    double stageInputDb()  const { return stage_.inputDb; }
    double stageOutputDb() const { return stage_.outputDb; }

    // ---------------------------------------------------------------- message thread: the panel ----

    // A dial, in degrees of its rotation. The crossfade dial moves continuously; any other captured dial
    // selects its nearest captured position; a tone or blend dial moves its own filter. False = no such
    // dial on this device, or a value it cannot take.
    bool setDial(const std::string& control, double degrees) {
        if (! loaded_) return false;
        if (! dial_.empty() && control == dial_) {
            dialDeg_ = std::clamp(degrees, 0.0, (double) dialSweep_);
            apply();
            return true;
        }
        for (const auto& c : stage_.device.controls)
            if (c.name == control) {
                const auto v = nearestValue(c, degrees);
                return ! v.empty() && setSwitch(control, v);
            }
        for (const auto& t : tones())
            if (t.name == control && t.sweep > 0) {
                toneAt_[control] = degreeValue(degrees, t.sweep);
                if (prepared_) rebuildKnob(t);
                return true;
            }
        for (const auto& b : stage_.blend)
            if (b.name == control && b.sweep > 0) {
                blendAt_[control] = degreeValue(degrees, b.sweep);
                applyBlendGains();
                return true;
            }
        return false;
    }

    // A control set to one of its values. For a captured axis this is namz::rig::resolve: the turned
    // control is pinned, everything else stays where the hand left it, and the closest captured
    // combination is chosen — the crossfade dial keeps its angle and follows along the new knots. False
    // = no such control, or a value nothing was captured at (then nothing changes).
    bool setSwitch(const std::string& control, const std::string& value) {
        if (! loaded_) return false;
        for (const auto& c : stage_.device.controls)
            if (c.name == control) {
                if (! dial_.empty() && control == dial_) {
                    const int deg = degreesOf(value);
                    return deg >= 0 && setDial(control, (double) deg);
                }
                if (namz::rig::resolve(stage_.device, axes_, control, value) == nullptr) return false;
                apply();
                return true;
            }
        for (const auto& t : tones())
            if (t.name == control) {
                // A SWITCH ONLY HAS THE POSITIONS IT DECLARES. Accepting any word and storing it read
                // back as a setting the knob does not have, while the sound was the reference — the
                // caller was told `true` and heard nothing. A dial keeps taking any degree: every angle
                // of its travel is playable, swept or not.
                if (t.sweep <= 0 && ! t.positions.empty()) {
                    bool found = false;
                    for (const auto& p : t.positions) if (p.value == value) { found = true; break; }
                    if (! found) return false;
                }
                toneAt_[control] = value;
                if (prepared_) rebuildKnob(t);
                return true;
            }
        for (const auto& b : stage_.blend)
            if (b.name == control) {
                blendAt_[control] = value;
                applyBlendGains();
                return true;
            }
        return false;
    }

    // What the panel reads back. `settings()` is the captured combination — for the crossfade dial it
    // names the knot at or below the hand, the angle itself is `dialDegrees()`.
    const Settings&    settings()    const { return axes_; }
    const std::string& dialName()    const { return dial_; }
    int                dialSweep()   const { return dialSweep_; }
    double             dialDegrees() const { return dialDeg_; }
    // A tone or blend knob's position as the pack spells it ("150", or a switch's token). Empty = no
    // such knob.
    std::string knobValue(const std::string& control) const {
        if (const auto t = toneAt_.find(control); t != toneAt_.end()) return t->second;
        if (const auto b = blendAt_.find(control); b != blendAt_.end()) return b->second;
        if (const auto a = axes_.find(control); a != axes_.end()) return a->second;
        return {};
    }
    const Selection& selection() const { return sel_; }
    const SlotPlan&  plan()      const { return plan_; }

    // TONE HANDED IN BESIDE THE MANIFEST — the same structures the pack carries, from another source: a
    // bench fitting bands by ear hands them here and hears them at once, then packs and hears them from
    // the file. ONE path applies tone; only where the knob's description came from differs. A knob named
    // here replaces the pack's block of that name (its `positions` or `sections`, whichever it had); a
    // name the pack has no block for is a new knob, starting at its default. The pack's own blocks are
    // untouched and return with clearToneOverride().
    void setToneOverride(std::vector<namz::rig::Tone> tone) {
        toneOverride_ = std::move(tone);
        for (const auto& t : tones()) if (! toneAt_.count(t.name)) toneAt_[t.name] = toneStart(t);
        if (prepared_ && loaded_) for (int s = 0; s < 2; ++s) { rebuildCurves(s); rebuildBands(s); }
    }
    void clearToneOverride() {
        if (toneOverride_.empty()) return;
        toneOverride_.clear();
        if (prepared_ && loaded_) for (int s = 0; s < 2; ++s) { rebuildCurves(s); rebuildBands(s); }
    }
    bool toneOverridden() const { return ! toneOverride_.empty(); }
    // The tone knobs as they PLAY: the pack's, with any handed-in block of the same name in its place,
    // and the handed-in knobs the pack does not have after them.
    std::vector<namz::rig::Tone> tones() const {
        std::vector<namz::rig::Tone> out;
        for (const auto& t : stage_.tone) {
            const namz::rig::Tone* use = &t;
            for (const auto& o : toneOverride_) if (o.name == t.name) { use = &o; break; }
            out.push_back(*use);
        }
        for (const auto& o : toneOverride_)
            if (std::none_of(stage_.tone.begin(), stage_.tone.end(), [&o](const namz::rig::Tone& t) { return t.name == o.name; }))
                out.push_back(o);
        return out;
    }

    // THE SHAPE OF THE HANDOVER — where the 50/50 lands between two captures and how wide the fade is.
    // The default is the original law (midpoint, full span); the bench turns these, a plugin need not.
    void setBlendShape(BlendShape s) { shape_ = s; if (loaded_) apply(); }
    BlendShape blendShape() const { return shape_; }
    // How much harder the top capture is driven per degree past it (ModelBlend.h explains the guess).
    void setTopExtendDbPerDeg(double v) { topExtend_ = v; if (loaded_) apply(); }
    double topExtendDbPerDeg() const { return topExtend_; }

    // Apply each model's loudness tag. Off, the models play raw, as the hardware returned them.
    void setNormalize(bool on) { normalize_.store(on, std::memory_order_release); }
    bool normalize() const { return normalize_.load(std::memory_order_acquire); }

    // Apply the pack's per-file input trims (`input_db` — a linked setting plays its neighbour
    // softer). OFF feeds every capture at unity instead: for a library shot at one honest level
    // the stated attenuations are somebody else's story, and into a nonlinear model a few dB
    // less in is a lot less out. Read on the audio side, so the toggle lands on the next block.
    void setInputTrims(bool on) { inputTrims_.store(on, std::memory_order_release); }
    bool inputTrims() const { return inputTrims_.load(std::memory_order_acquire); }

    // The rest after which a silent slot's model is no longer run (kColdAfterSeconds). A slot that
    // sounds at all is never cold; between two captures both models run, on one they do not. Zero or
    // less = never. Takes effect from the next block; a slot already asleep stays asleep until woken.
    void   setColdAfterSeconds(double seconds) {
        coldSeconds_ = seconds;
        coldAfter_.store(coldSamples(seconds), std::memory_order_release);
    }
    double coldAfterSeconds() const { return coldSeconds_; }

    // The models' measured offsets (AlignmentTable.h), for a pack that does not carry its own. A delay
    // travels with a model and is applied at the one instant its slot is silent — weight exactly zero
    // — never as a splice on a live signal. So a table handed in BEFORE playing lands with the first
    // loads; one that arrives mid-mix waits until the dial visits a knot. A host that has the bytes in
    // hand should measure before it plays.
    void setAlignment(AlignmentTable table) {
        align_ = std::move(table);
        if (loaded_) stageDelays();
    }
    const AlignmentTable& alignment() const { return align_; }
    bool alignmentFromPack() const { return align_.fromPack && ! align_.empty(); }

    // Regular message-thread housekeeping: frees models the audio thread has stepped past, retries a
    // filter the convolver rejected mid-fade. Call from a timer, a few times a second at least — and
    // ask for the load job after it (takeLoadJob), because the law waits for that.
    void service() {
        for (auto& n : nam_) n.collectGarbage();
        for (auto& f : fir_) f.flushPending();
        dry_.flushPending();
    }

    // ---------------------------------------------------------------- message thread: loading ----

    // A MODEL THE LAW ASKED FOR, as work for the host to run where it likes — a worker thread, or right
    // here. The job carries everything the work needs (the source, the bytes if already fetched, the
    // numbers the stages were prepared with) and nothing of the player, so run() is a pure function of
    // it. At most ONE job is out at a time, and a job taken MUST come back through deliver() — with a
    // null model if it failed — or the law waits for it forever.
    struct LoadJob {
        int           slot = 0;
        std::uint64_t model = 0;                                   // BlendModelId: models_ index + 1
        std::string   fileId;                                      // the pack's `files[].id`
        std::shared_ptr<const std::vector<std::byte>> bytes;       // null = not fetched yet
        ModelSource   source;
        double        sampleRate = 48000.0;
        int           maxBlock = 512;
        std::uint64_t generation = 0;                              // the load() it belongs to
    };
    struct Loaded {
        LoadJob job;
        std::shared_ptr<const std::vector<std::byte>> bytes;       // the job's own, or what run() fetched
        bool    fetched = false;                                   // …and whether it had to
        felitronics::nam::NamStage::PreparedModel model;           // null = the bytes were not a model
    };

    std::optional<LoadJob> takeLoadJob() {
        if (jobOut_ || ! loaded_) return std::nullopt;
        const auto ask = loadAsk_.load(std::memory_order_acquire);
        if (ask == 0) return std::nullopt;
        const auto id = ask >> 1;
        if (id == 0 || id > models_.size()) {                      // an ask that outlived its pack
            loadAsk_.store(0, std::memory_order_release);
            landFail_.store((int) (ask & 1u) + 1, std::memory_order_release);
            return std::nullopt;
        }
        const auto& m = models_[(std::size_t) id - 1];
        LoadJob job;
        job.slot = (int) (ask & 1u); job.model = id; job.fileId = m.id; job.bytes = m.bytes;
        job.source = source_; job.sampleRate = fs_; job.maxBlock = maxBlock_; job.generation = gen_;
        jobOut_ = true;
        return job;
    }

    // Any thread but the audio one. The bytes through the source when the job has none, then the heavy
    // half of a load (NamStage::prepareModel). Touches nothing of any player.
    static Loaded run(LoadJob job) {
        Loaded out;
        out.bytes = job.bytes;
        if (out.bytes == nullptr && job.source) {
            out.bytes = std::make_shared<const std::vector<std::byte>>(job.source(job.fileId));
            out.fetched = true;
        }
        if (out.bytes != nullptr && ! out.bytes->empty())
            out.model = felitronics::nam::NamStage::prepareModel(out.bytes->data(), out.bytes->size(),
                                                                 job.sampleRate, job.maxBlock);
        out.job = std::move(job);
        return out;
    }

    // Message thread. The light half: the bytes kept (fetched once per file, however many times the
    // dial crosses it), the model installed in its slot, the law told it landed — or failed. A job from
    // a pack that is gone is dropped whole: the law was told to start over when the pack left.
    void deliver(Loaded loaded) {
        jobOut_ = false;
        auto& job = loaded.job;
        if (job.generation != gen_ || ! loaded_) return;
        if (loaded.fetched && job.model >= 1 && job.model <= models_.size()) {
            auto& m = models_[(std::size_t) job.model - 1];
            if (m.bytes == nullptr) { m.bytes = loaded.bytes; ++loads_; }
        }
        loadSlot(job.slot, std::move(loaded.model), job.model);
    }

    // service() with every job due run right here, on this thread — for a host without a worker, and
    // for a test. The same path on one thread, not a second one.
    void serviceHere() {
        service();
        while (auto job = takeLoadJob()) deliver(run(std::move(*job)));
    }

    // Host-rate latency: the models' own rate-matching (the filters are minimum-phase, the alignment
    // delays are relative and the reference model carries none).
    int latencySamples() const {
        return std::max(nam_[0].latencySamples(), nam_[1].latencySamples());
    }

    // ---- read-outs of the audio thread, for a strip or a dump ----
    float liveMix() const { return liveMix_.load(std::memory_order_acquire); }        // applied weight of slot 1
    // Applied blend gains — what the audio thread multiplied by on the last block, not what was asked for.
    // Same read-out contract as liveMix() above: a strip can show the gain that is actually in the sound
    // while the ramp is still travelling, and it is the only place a law-8 regression on these ramps is
    // visible at all (the stuck value is a subnormal that adds nothing to a normal sample).
    float liveDry() const { return liveDry_.load(std::memory_order_acquire); }
    float liveWet() const { return liveWet_.load(std::memory_order_acquire); }
    // The file id held in a slot right now ("" = nothing yet).
    std::string heldFileId(int slot) const {
        const auto id = held_[(std::size_t) (slot & 1)].load(std::memory_order_relaxed);
        return id == 0 || id > models_.size() ? std::string() : models_[(std::size_t) id - 1].id;
    }
    int appliedSlotDelay(int slot) const { return slotDelay_[(std::size_t) (slot & 1)].load(std::memory_order_relaxed); }
    // THE LOUDNESS TAG OF WHAT IS SOUNDING — of the slot the sound actually leaves by, which is not
    // always slot 0. Slots are handed out by the PARITY of a capture's place on the dial (slotPlan,
    // RigSelection.h), so on an odd rung the whole sound comes out of slot 1 while slot 0 holds the
    // silent neighbour: a face that read slot 0 named the wrong capture, or warned "plays raw" about a
    // model that carries a tag.
    //
    // DURING A CROSSFADE THERE IS NO SINGLE HONEST NUMBER, and this does not invent one. `db` and
    // `tagged` are the slot carrying MOST of the sound, and `blended` says the other slot holds a
    // DIFFERENT capture and is audible beside it — so a face can say "one of two" instead of stating a
    // level the sound has not got. At the ends of the dial both slots hold the same capture, and there
    // `blended` is false however the weight sits.
    //
    // The weight is the APPLIED one (liveMix()), not the one asked for: while a model is still loading
    // or warming the slot that IS the sound is the old one, and its tag is what a face should show.
    // `tagged` is false both for a model without a tag and for a slot with no model at all —
    // heldFileId(slot) tells those two apart. This is the model mix alone: a blend knob at its dry end
    // takes every model out of the sound and says nothing here. Message thread.
    struct SoundingLoudness {
        double db      = 0.0;      // the tag of the slot carrying most of the sound, in dB
        bool   tagged  = false;    // …and whether that slot carries one at all
        bool   blended = false;    // …and whether a second, different capture is audible beside it
        int    slot    = 0;        // which slot that was
    };
    SoundingLoudness soundingLoudness() const {
        // unload() empties both stages at once while the audio thread clears its own numbers a block
        // later: with no pack there is nothing sounding to carry a tag.
        if (! loaded_) return {};
        const float mix = liveMix_.load(std::memory_order_acquire);
        const auto  h0  = held_[0].load(std::memory_order_relaxed);
        const auto  h1  = held_[1].load(std::memory_order_relaxed);
        const int   s   = mix > 0.5f ? 1 : 0;
        return { nam_[s].modelLoudness(), nam_[s].modelHasLoudness(),
                 mix > 0.0f && mix < 1.0f && h0 != 0 && h1 != 0 && h0 != h1, s };
    }
    // Instrumentation, for a bench's dump: blocks in which a slot was still warming, weight moves of
    // more than 2 % inside one block and the biggest of them — a whole swing inside one block is a step
    // in all but name. Counted since the last clearCounters().
    int   warmBlocks()  const { return warmBlocks_.load(std::memory_order_relaxed); }
    int   mixJumps()    const { return mixJumps_.load(std::memory_order_relaxed); }
    float biggestJump() const { return biggestJump_.load(std::memory_order_relaxed); }
    // Whether a slot is COLD right now — held, silent, its model no longer run on AUDIO
    // (kColdAfterSeconds; it is still clocked at width zero until its drain is spent) — and the blocks
    // it has slept since the last clearCounters(), for the dump beside warmBlocks().
    bool  slotCold(int slot)   const { return slotCold_[(std::size_t) (slot & 1)].load(std::memory_order_relaxed); }
    int   coldBlocks(int slot) const { return coldBlocks_[(std::size_t) (slot & 1)].load(std::memory_order_relaxed); }
    void  clearCounters() {
        warmBlocks_.store(0, std::memory_order_relaxed);
        mixJumps_.store(0, std::memory_order_relaxed);
        biggestJump_.store(0.0f, std::memory_order_relaxed);
        for (auto& c : coldBlocks_) c.store(0, std::memory_order_relaxed);
    }
    long long modelLoads() const { return loads_; }                                   // how many the source was asked for
    int bandsDropped() const { return bandsDropped_; }                                // sections past kMaxBands

    // The tone as it stands, for drawing: the summed curve of the curve-form knobs on a side (dB on
    // `commonGrid()`, empty = nothing on that side) and the bands of the section-form knobs.
    const std::vector<double>& commonGrid() const { return commonGrid_; }
    const std::vector<double>& curveDb(int side) const { return curveSum_[side & 1]; }
    bool curveActive(int side) const { return ! firBypass_[side & 1].load(std::memory_order_acquire); }   // …and whether it is a FIR or a wire
    const std::vector<felitronics::rigplayer::SectionBiquad>& bands(int side) const { return bandsShown_[side & 1]; }
    BlendGains blendGains() const { return { (double) dryGain_.load(std::memory_order_relaxed),
                                             (double) wetGain_.load(std::memory_order_relaxed) }; }

    // ---------------------------------------------------------------------- the audio thread ----

    // In place. `numChannels` up to what prepare() was given (fewer is fine: the missing planes are
    // silent). Any `numSamples`; longer than maxBlock is walked in pieces.
    //
    // THE MODELS RUN ON THE PLANES THEY ARE GIVEN, not on the width the player was prepared for.
    // A host prepared for a stereo bus that plays one plane — a mono chain — used to push the silent
    // second plane through both networks as well: four WaveNet passes a block for one channel of
    // sound, and a quarter of a core gone to nothing. The convolvers still take the prepared width
    // (the convolvers take the prepared width and REFUSE any other — law 11c), and the spare planes stay
    // zeroed for them; only the expensive part — the two models — shrinks to what is playing.
    [[nodiscard]] bool process(float* const* io, int numChannels, int numSamples) {
        if (numChannels < 0 || numSamples < 0) return false;         // malformed — law 11
        if (! prepared_ || io == nullptr) return false;
        if (numChannels > channels_) return false;                   // width is a LIMIT — law 11(b)
        if (numSamples == 0) return true;               // law 11(d): no samples, no time, no edge
        // `numChannels == 0` is NOT short-circuited here: it is a GAP, and the pipeline below is what
        // makes the gap real — the per-slot lagTail_ fill on the falling edge, and the three convolvers,
        // which run on the zeroed spare planes and go on decaying. Returning early left both slots
        // frozen: measured 0.2317 out of DIGITAL SILENCE (-12.7 dBFS) on the return, on both slots.
        const int nch = numChannels;
        bool ok = true;                                              // the stages' verdicts, ANDed

        // A PLANE THAT STOPS PLAYING AND PLAYS AGAIN. The models run on the planes they are given, and so
        // do the per-slot alignment delay lines beside them — `for (int c = 0; c < nch; …)` below. A plane
        // the host stops handing over therefore keeps its `lagTail_` frozen rather than draining it, and
        // hands it back when the host widens again. The tail is exactly the slot's delay in samples, so its
        // loudness is whatever the plane was carrying — a ceiling, not a sample: it sits after the slot's
        // gain and before the blend, so it cannot exceed amplitude x slot weight, and a sweep of all 218
        // leaving phases at sample resolution reaches 99.99% of that. Out of DIGITAL SILENCE: up to 0.25,
        // which is -12.0 dBFS, in the first three samples of the return, none of it on the plane that stayed.
        //
        // 🔴 THAT ARITHMETIC IS THIS LINE'S OWN AND IT DOES NOT COVER THE MODELS' HALF — a whole paragraph
        // here used to spend it on both, and it was wrong by a factor of 3.8. Both are stated now because
        // the two ceilings are DIFFERENT SHAPES, not one number applied twice:
        //
        //     the delay line   only ONE slot carries the delay (delayOf() is maxLag − lag), and the tail is
        //                      a literal replay, so the operator's gain is exactly 1: A · w = 0.5 · 0.5.
        //     the models       BOTH slots freeze at once, so the weights SUM rather than pick — and the
        //                      frozen rate-matcher is not a replay: a phase row of the polyphase table is
        //                      normalised by its SUM, not by its modulus, so its peak gain exceeds one.
        //
        // In one formula, with `g_s` the frozen return operator of slot s (the capture's impulse response,
        // convolved with the rate-matcher pair where one is installed):
        //
        //          ceiling  =  A · Σ_s w_s · ‖g_s‖₁
        //
        // and it is ATTAINABLE, because the whole path is linear in the pre-gap input for a Linear capture:
        // the return is `y[n] = Σ_k h_n[k]·x[−k]`, so `A·max_n ‖h_n‖₁` is reached exactly by the input
        // `x[−k] = A·sign(h_n[k])`. Measured that way, `A` = 0.5, both slots at 0.5, a unit capture:
        //
        //     host      44100     48000    88200    96000   176400   192000
        //     ‖h‖₁    1.898766  1.000000 1.685674 1.612676 1.505512 1.492565   (48 kHz: no rate-matcher)
        //     ceiling 0.949383  0.500000 0.842837 0.806338 0.752756 0.746283
        //     reached  100.00%   100.00%  100.00%  100.00%  100.00%  100.00%   by the sign pattern
        //     a sine     55.4%    99.91%    62.4%    65.0%    69.0%    70.0%   swept over EVERY leaving phase
        //
        // The last row is why the number above is not the ceiling of this half: the phase sweep that closed
        // the delay line reaches 55 % here, and would have published 0.5257. It agreed there because THAT
        // operator is a pure delay — a property of the operator, not of the method.
        //
        // 🟢 THE MODELS' HALF IS CLOSED, in nam::NamStage rather than here: a lane the host stops handing
        // over is fed the digital silence it is actually receiving, for exactly as long as its state can
        // still be heard (its receptive field, plus the rate-matcher's tap windows at each end). Out of
        // digital silence the return is now EXACTLY zero at 8 · 22.05 · 44.1 · 48 · 88.2 · 96 · 176.4 ·
        // 192 kHz, with a memoryless capture and with a real 6332-sample WaveNet, at gap widths 0 and 1
        // and on both slots. The SLEEPING slot — this file's other way of not clocking a stage — is closed
        // beside the models below. What is not: an LSTM's recurrent cell has no finite drain, so its
        // number is a bound on NAM's own heuristic (0.419 against 0.023 for a lane clocked throughout).
        // The comment beside those lines says the history is "advanced every block including at zero" —
        // true for the planes that are playing, which is exactly the gap. Dropping it here matches what
        // this file already does in three other places where a slot's delay line stops being valid: a
        // landing, a retime, and falling asleep all fill it. The convolvers and the tone are not touched:
        // they are handed the PREPARED width with the spare planes zeroed, so they never stop running.
        if (nch < ranNch_)
            for (int i = 0; i < 2; ++i)
                for (int c = nch; c < ranNch_; ++c) lagTail_[i][(std::size_t) c].fill(0.0f);
        ranNch_ = nch;
        const bool norm = normalize_.load(std::memory_order_acquire);
        float* a[kMaxChannels]; float* b[kMaxChannels]; float* d[kMaxChannels];
        for (int c = 0; c < kMaxChannels; ++c) { b[c] = slotB_[c].data(); d[c] = dryBuf_[c].data(); }

        int done = 0;
        while (done < numSamples) {
            const int count = std::min(maxBlock_, numSamples - done);
            for (int c = 0; c < channels_; ++c) {
                if (c < nch && io[c] != nullptr) a[c] = io[c] + done;
                else { a[c] = spare_[c].data(); std::fill(a[c], a[c] + count, 0.0f); }
            }
            for (int s = 0; s < 2; ++s) takeBands(s);

            // THE PACK'S INPUT LEVEL, and it goes first: before the dry block, before the pre-model
            // tone, before everything. It is the level of the GUITAR into this device, and a blend
            // knob mixes one guitar with itself — feed the models less and the dry side the same as
            // ever and the mix becomes two instruments at two volumes. `extendDb` is a different
            // animal and stays in chainGain_ below: a trick played inside the pack past the top
            // capture, which the dry path leaving at the input jack never saw.
            rampInto(a, nch, count, inGain_.load(std::memory_order_acquire), curIn_);

            // Keep the DRY block before anything touches it: the same DI the models are fed, which is
            // the whole economy of a blend. Taken ahead of the pre-model tone — the hardware's dry path
            // leaves at the input jack.
            const bool mixDry = dryActive_.load(std::memory_order_acquire);

            // 🔴 THE DRY PATH IS DELAYED BY WHAT THE WET PATH COSTS, and the delay is taken UNCONDITIONALLY.
            // The blend below sums `a` (through the models, hence through their rate-match) with `d`; those
            // two are misaligned by exactly latencySamples() unless something holds the dry back, and
            // summing them notches at fs/(2·D). See prepare() for the numbers and the history.
            // Advancing only while the blend is ON would be the classic cold-ring bug — DryAligner.h says
            // it in as many words: a ring fed only while a stage runs is COLD the moment it is first read
            // and emits its latency in zeros. Turning the dry knob up would then start with D samples of
            // silence in the dry leg. The copy costs one pass over the block when the blend is off, which
            // is nothing beside two neural models.
            for (int c = 0; c < channels_; ++c) std::copy(a[c], a[c] + count, d[c]);
            dryLatency_.advance((const float* const*) d, channels_, count, latencySamples());
            for (int c = 0; c < channels_; ++c) std::copy_n(dryLatency_.delayed(c), count, d[c]);

            // The knobs that sit BEFORE the distortion in the hardware, shaping what gets distorted.
            runBands(0, a, count);
            if (! firBypass_[0].load(std::memory_order_acquire)) ok = fir_[0].process(a, channels_, count) && ok;

            // The chain's trim: past the top capture it rises with the angle, and a fast hand moves
            // tens of degrees between two events — so it is smoothed like every other gain here.
            rampInto(a, nch, count, chainGain_.load(std::memory_order_acquire), curChain_);

            // THE LAW, once per block, and the only writer of the weight. Everything it needs arrives
            // through atomics; everything it asks for leaves the same way.
            if (forget_.exchange(false, std::memory_order_acquire)) { blend_ = {}; coldRt_[0] = coldRt_[1] = false; }
            if (const int f = landFail_.exchange(0, std::memory_order_acquire); f != 0)
                felitronics::nam::blendLoadFailed(blend_, f - 1);
            if (landFlag_.exchange(false, std::memory_order_acquire)) {
                const int ls = landSlot_.load(std::memory_order_relaxed);
                const int nd = landDelay_.load(std::memory_order_relaxed);
                slotDelay_[(std::size_t) (ls & 1)].store(nd, std::memory_order_relaxed);
                pendDelay_[(std::size_t) (ls & 1)].store(nd, std::memory_order_release);
                for (auto& t : lagTail_[(std::size_t) (ls & 1)]) t.fill(0.0f);
                felitronics::nam::blendLanded(blend_, ls, landModel_.load(std::memory_order_relaxed),
                                              landWarm_.load(std::memory_order_relaxed));
            }
            // ACQUIRE FIRST: the writer publishes the two ids and then releases on the target.
            felitronics::nam::BlendRequest req;
            req.targetB = (double) reqTarget_.load(std::memory_order_acquire);
            req.want[0] = reqWant_[0].load(std::memory_order_relaxed);
            req.want[1] = reqWant_[1].load(std::memory_order_relaxed);
            felitronics::nam::BlendPolicy policy;
            policy.coldAfterSamples = coldAfter_.load(std::memory_order_relaxed);
            const auto law = felitronics::nam::blendStep(blend_, req, count, policy);
            if (law.load.wanted && loadAsk_.load(std::memory_order_relaxed) == 0)
                loadAsk_.store((law.load.model << 1) | (std::uint64_t) (law.load.slot & 1), std::memory_order_release);
            for (int i = 0; i < 2; ++i) {
                held_[(std::size_t) i].store(blend_.held[i], std::memory_order_relaxed);
                slotCold_[(std::size_t) i].store(blend_.cold[i], std::memory_order_relaxed);
                // FALLING ASLEEP CLEARS THE DELAY LINE, exactly as a landing does (above): the tail would
                // otherwise hold the slot's last samples from before the rest, and a model that declares
                // no field (need == 0) is heard on the very block it wakes in — with those samples first.
                if (blend_.cold[i] && ! coldRt_[i])
                    for (auto& t : lagTail_[(std::size_t) i]) t.fill(0.0f);
                coldRt_[i] = blend_.cold[i];
                if (blend_.cold[i]) coldBlocks_[(std::size_t) i].fetch_add(1, std::memory_order_relaxed);
            }
            // …and a staged retime lands the same way a model does: only where the slot is silent.
            for (int i = 0; i < 2; ++i) {
                const double w = i == 0 ? 1.0 - law.endB : law.endB;
                const int want = pendDelay_[(std::size_t) i].load(std::memory_order_acquire);
                if (w <= 0.0 && want != slotDelay_[(std::size_t) i].load(std::memory_order_relaxed)) {   // exactly zero: the law clamps to [0, 1]
                    slotDelay_[(std::size_t) i].store(want, std::memory_order_relaxed);
                    for (auto& t : lagTail_[(std::size_t) i]) t.fill(0.0f);
                }
            }
            liveMix_.store((float) law.endB, std::memory_order_release);
            if (blend_.fed[0] < blend_.need[0] || blend_.fed[1] < blend_.need[1])
                warmBlocks_.fetch_add(1, std::memory_order_relaxed);
            if (const float jump = (float) std::abs(law.endB - law.beginB); jump > 0.02f) {
                mixJumps_.fetch_add(1, std::memory_order_relaxed);
                float prev = biggestJump_.load(std::memory_order_relaxed);
                while (jump > prev && ! biggestJump_.compare_exchange_weak(prev, jump, std::memory_order_relaxed)) {}
            }

            // THE SECOND CAPTURE'S COPY is taken here, after the pre-model tone and the chain trim, so
            // both models are fed the identical signal — and only now do the two part company, each
            // through its own trim: a linked setting plays its neighbour's model with less going in.
            // A COLD SLOT IS NOT RUN ON AUDIO — not copied into, not trimmed, not modelled, not delayed,
            // not even mixed: the law holds its weight at exactly zero for as long as it sleeps, so the
            // other slot IS the sound, and its buffer — whatever it holds, however old — is left out
            // rather than multiplied by zero (a NaN times zero is a NaN, for ever). What it IS handed is
            // a width-zero call (below), so the stage can drain the silence the slot is receiving; that
            // is law 11(d)'s clock-only call and it stops of its own accord. Its trim ramp resumes from
            // where it stopped and settles within its block, like any other gain change; its delay line
            // was cleared on the way to sleep. (The rail is checked as well as the flag: the flag says
            // what the law decided, the weights say what this block sounds like.)
            const bool run[2] { ! blend_.cold[0], ! blend_.cold[1] };
            const bool aAlone = ! run[1] && law.beginB <= 0.0 && law.endB <= 0.0;
            const bool bAlone = ! run[0] && law.beginB >= 1.0 && law.endB >= 1.0;
            const bool trims = inputTrims_.load(std::memory_order_acquire);
            if (run[1]) for (int c = 0; c < nch; ++c) std::copy(a[c], a[c] + count, b[c]);
            if (run[0]) rampInto(a, nch, count, trims ? slotGain_[0].load(std::memory_order_acquire) : 1.0f, curSlot_[0]);
            if (run[1]) rampInto(b, nch, count, trims ? slotGain_[1].load(std::memory_order_acquire) : 1.0f, curSlot_[1]);
            // 🔴 A SLEEPING SLOT IS STILL HANDED THE CALL, at width ZERO. Skipping it entirely is the
            // SECOND way this file stops clocking a NamStage — the first is a plane the host takes away
            // — and it has the same consequence: the slot's models and rate-matchers freeze holding
            // whatever was playing when the law put it to sleep, and a wake hands it back. Measured on a
            // three-knot device whose 240 capture has a 2001-sample memory: the dial parked at 150 until
            // slot 0 slept, eight blocks of DIGITAL SILENCE, then a turn to 240 — per-block peaks
            // 0.114 · 0.22 · 0.366 · 0.486 · 0.5 · 0.5 · 0.5 · 0.5, i.e. **0.500000 out of digital
            // silence** for exactly one receptive field, under the law's own 0.25-per-call ramp. A
            // TWO-knot device cannot show it: both slots end up holding the same capture, so the turn
            // SWAPS the backend rather than waking it, and a swap has nothing to leak. That is why the
            // first fixture read a clean zero, and it is a fixture fault, not a property of the code.
            // Width zero is law 11(d)'s clock-only call: the stage drains the lane and stops, so this
            // costs one bounded drain per sleep and nothing thereafter — a sleeping slot is still free.
            // Measured on a real 6332-sample capture: 12.9 ms of CPU at 48 kHz and 9.3 at 44.1 for the
            // whole sleep, and 0.0000 ms per block once the debt is spent.
            // ⚠️ ONE CONTRACT CONSEQUENCE, stated because nobody would look for it: a cold slot's verdict
            // is now ANDed too, and `NamStage::process` refuses for an UNPREPARED backend. A slot whose
            // live re-prepare was refused (a low-memory `configureRates`) used to be silent-and-accepted
            // while it slept and now fails every block instead. That is law 11's answer — a call that
            // cannot be honoured says so — and it is not constructible in a test without exhausting
            // memory, so it is written here rather than gated.
            ok = nam_[0].process(a, run[0] ? nch : 0, count, norm) && ok;
            ok = nam_[1].process(b, run[1] ? nch : 0, count, norm) && ok;
            // Align BEFORE the weights: during a ramp the two gains must sum to one at the SAME instant.
            // Each slot carries its own history per channel, advanced every block including at zero.
            for (int c = 0; c < nch; ++c) {
                if (run[0]) felitronics::nam::blendDelay(a[c], lagTail_[0][(std::size_t) c].data(), kMaxDelay,
                                                         slotDelay_[0].load(std::memory_order_acquire), count);
                if (run[1]) felitronics::nam::blendDelay(b[c], lagTail_[1][(std::size_t) c].data(), kMaxDelay,
                                                         slotDelay_[1].load(std::memory_order_acquire), count);
                if (aAlone)      { /* slot 1 asleep at zero: `a` is the whole sound */ }
                else if (bAlone) std::copy(b[c], b[c] + count, a[c]);
                else             felitronics::nam::blendMix(a[c], b[c], count, law.beginB, law.endB);
                // …and the law's own gain, which is below one only while nothing has been fed yet.
                if (law.beginGain < 1.0 || law.endGain < 1.0) {
                    const float dg = (float) (law.endGain - law.beginGain) / (float) count;
                    float g = (float) law.beginGain;
                    for (int i = 0; i < count; ++i) { g += dg; a[c][i] *= g; }
                }
            }

            // The knobs AFTER the distortion, shaping what came out.
            runBands(1, a, count);
            if (! firBypass_[1].load(std::memory_order_acquire)) ok = fir_[1].process(a, channels_, count) && ok;

            // The blend, as the hardware sums it: the dry path through its own response, then the two
            // gains the pack states for this position.
            if (mixDry) {
                if (! dryFirBypass_.load(std::memory_order_acquire)) ok = dry_.process(d, channels_, count) && ok;
                const float wantDry = dryGain_.load(std::memory_order_acquire);
                const float wantWet = wetGain_.load(std::memory_order_acquire);
                const float decay   = std::exp(-(float) count / (float) (0.010 * fs_));
                const float endDry  = rampEnd(wantDry, curDry_, decay);   // law 8 — see rampEnd
                const float endWet  = rampEnd(wantWet, curWet_, decay);
                const float stepDry = (endDry - curDry_) / (float) count;
                const float stepWet = (endWet - curWet_) / (float) count;
                for (int c = 0; c < channels_; ++c) {
                    float gd = curDry_, gw = curWet_;
                    for (int i = 0; i < count; ++i) {
                        gd += stepDry; gw += stepWet;
                        a[c][i] = a[c][i] * gw + d[c][i] * gd;
                    }
                }
                curDry_ = endDry; curWet_ = endWet;
                liveDry_.store(endDry, std::memory_order_release);
                liveWet_.store(endWet, std::memory_order_release);
            }

            // THE PACK'S OUTPUT LEVEL, and it goes last: after the mix. One number for the whole
            // stage, and applied any earlier it would ride the wet side alone and move a blend the
            // pack states for this position.
            rampInto(a, nch, count, outGain_.load(std::memory_order_acquire), curOut_);
            done += count;
        }
        return ok;
    }

private:
    // ------------------------------------------------------------------ the decision, applied ----

    static std::string degreeValue(double degrees, int sweep) {
        return std::to_string((int) std::lround(std::clamp(degrees, 0.0, (double) sweep)));
    }

    static float linOf(double db) { return std::abs(db) < 0.0005 ? 1.0f : (float) std::pow(10.0, db / 20.0); }

    /// WHAT THE PLAYER IS APPLYING, decided in ONE place: the host's hand if it has one, else the
    /// pack's own number. Called from apply(), from the host setters, and from load() before anything
    /// slow — a pack swap used to publish unity here (unload's reset) and the real gains only at the
    /// end of load(), so a callback landing in between heard neither pack's level.
    void publishLevels() {
        inGain_ .store(linOf(hostIn_  ? *hostIn_  : stage_.inputDb),  std::memory_order_release);
        outGain_.store(linOf(hostOut_ ? *hostOut_ : stage_.outputDb), std::memory_order_release);
    }

    std::uint64_t modelIdOf(int file) const {
        return file < 0 || (std::size_t) file >= fileModel_.size() ? 0 : (std::uint64_t) fileModel_[(std::size_t) file] + 1;
    }

    const std::string& fileIdOfModel(std::uint64_t id) const {
        static const std::string none;
        return id == 0 || id > models_.size() ? none : models_[(std::size_t) id - 1].id;
    }

    // In HOST samples. A pack's numbers are in the model's own rate, which the stage knows once the
    // model is loaded; a measured table carries the rate it was measured at.
    int delayOfModel(std::uint64_t id, const felitronics::nam::NamStage& st) const {
        // Same reason as warmFor(): a pack's lag is written in the model's OWN rate, and for an
        // untagged capture that rate is not the -1 the stage reports — it is the rate the stage will
        // actually run it at. AlignmentTable falls back to scale 1.0 when handed a non-positive rate,
        // which for an untagged model at 96 kHz would apply a 48 kHz lag unscaled.
        const double runSR = felitronics::nam::NamStage::rateMatch(fs_, st.modelSampleRate()).modelRunSR;
        return align_.delayOf(fileIdOfModel(id), fs_, runSR);
    }

    // The whole decision for the panel as it stands, handed to the audio thread as a request. What the
    // law does with it — in what order, and whether a model may be replaced yet — is its business.
    void apply() {
        const auto& d = stage_.device;
        publishLevels();   // ahead of the early return: they are the stage's, not the selection's
        sel_  = select(d, axes_, dial_, dialDeg_, shape_, topExtend_);
        plan_ = slotPlan(sel_, d);
        if (plan_.file[0] < 0) { setRequest(0, 0, 0.0f); return; }
        setRequest(modelIdOf(plan_.file[0]), modelIdOf(plan_.file[1]), (float) plan_.targetB);
        for (int i = 0; i < 2; ++i) slotGain_[i].store(linOf(plan_.inputDb[i]), std::memory_order_release);
        chainGain_.store(linOf(sel_.extendDb), std::memory_order_release);
        stageDelays();
        // The dial's own entry names the knot at (or below) the hand, so a later resolve() on another
        // control scores "keep what the user had" against the capture actually sounding.
        if (! dial_.empty() && sel_.fileA >= 0)
            if (const auto v = d.files[(std::size_t) sel_.fileA].settings.find(dial_);
                v != d.files[(std::size_t) sel_.fileA].settings.end())
                axes_[dial_] = v->second;
    }

    // For the slots as loaded now; a model still to come brings its own delay when it lands.
    void stageDelays() {
        for (int i = 0; i < 2; ++i)
            pendDelay_[(std::size_t) i].store(delayOfModel(modelIdOf(plan_.file[i]), nam_[i]), std::memory_order_release);
    }

    void setRequest(std::uint64_t want0, std::uint64_t want1, float targetB) {
        reqWant_[0].store(want0, std::memory_order_relaxed);
        reqWant_[1].store(want1, std::memory_order_relaxed);
        reqTarget_.store(std::clamp(targetB, 0.0f, 1.0f), std::memory_order_release);
    }

    // Put this model in that slot and tell the law it landed, with the delay that travels with the
    // model — applied when it lands, the one instant the slot is guaranteed silent. Only the ask this
    // answers is cleared, and a failure carries its own slot: sharing either cell stranded the other
    // slot for the rest of a session.
    bool loadSlot(int slot, felitronics::nam::NamStage::PreparedModel prepared, std::uint64_t model) {
        auto& stage = nam_[(std::size_t) (slot & 1)];
        const bool ok = prepared != nullptr && stage.install(std::move(prepared));
        const int delay = ok ? delayOfModel(model, stage) : 0;    // after the load: the rate is the model's
        if (const auto ask = loadAsk_.load(std::memory_order_acquire); ask != 0 && (int) (ask & 1u) == slot)
            loadAsk_.store(0, std::memory_order_release);
        if (! ok) { landFail_.store(slot + 1, std::memory_order_release); return false; }
        landSlot_.store(slot, std::memory_order_relaxed);
        landModel_.store(model, std::memory_order_relaxed);
        landWarm_.store(warmFor(stage), std::memory_order_relaxed);
        landDelay_.store(std::clamp(delay, 0, kMaxDelay), std::memory_order_relaxed);
        landFlag_.store(true, std::memory_order_release);
        return true;
    }

    // How many samples this model owes before it may be heard, in THIS rate: its receptive field,
    // scaled — a 96 kHz host feeds twice as many to fill the same network — plus the RATE-MATCH DELAY,
    // plus one block, so a slot is never marked audible for a block it is still short in.
    //
    // 🔴 THE LATENCY TERM IS NOT DECORATION, and it used to be missing. When the stage is resampling,
    // the first latencySamples() host samples out of it are the network's response to the resampler's
    // own leading zeros, not to the signal — so a full receptive field of REAL material needs that many
    // more. The old expression relied on `+ maxBlock_` to absorb it silently, which worked only while
    // that delay was 3.84 samples: P34's 64-tap kernel makes it 61.4 at 44.1 kHz, and at the 32- and
    // 64-sample blocks live rigs run the slack stopped covering it. The effect is small against a
    // ~6300-sample receptive field (a crossfade starting up to 1.4 ms early), but the sentence above
    // was a GUARANTEE, and a guarantee that quietly stopped holding is worse than a smaller number.
    long long warmFor(const felitronics::nam::NamStage& st) const {
        const int pre = st.prewarmSamples();
        if (pre <= 0) return 0;
        // 🔴 ASK for the rate the model will be RUN at, do not read the rate it REPORTS. An untagged
        // capture reports -1 and NamStage runs it at kModelSampleRate anyway; the previous line here
        // restated the normalisation with a different answer — `scale = 1.0` — and therefore
        // under-warmed such a model by a whole receptive field at a 96 kHz host, which is exactly the
        // guarantee the comment above says must not quietly stop holding.
        const double mr = felitronics::nam::NamStage::rateMatch(fs_, st.modelSampleRate()).modelRunSR;
        const double scale = (mr > 0.0 && fs_ > 0.0) ? fs_ / mr : 1.0;
        return (long long) std::ceil((double) pre * scale) + (long long) st.latencySamples() + maxBlock_;
    }

    // The rest before a slot goes cold, in this rate's samples, for the law. Zero = never.
    long long coldSamples(double seconds) const {
        return seconds > 0.0 ? (long long) std::llround(seconds * fs_) : 0;
    }

    // ------------------------------------------------------------------------------- the tone ----

    // Which side of the models a knob sits on: 0 before, 1 after, -1 for a placement this player does
    // not know — such a knob is skipped, as the format says.
    static int sideOf(const namz::rig::Tone& t) {
        return t.placement == "pre" ? 0 : t.placement == "post" ? 1 : -1;
    }

    void rebuildKnob(const namz::rig::Tone& t) {
        const int side = sideOf(t);
        if (side < 0) return;
        // Bands, whether they travel with a dial or stand still at a switch's position; else the curve.
        if (! t.sections.empty() || bandsPerPosition(t)) rebuildBands(side); else rebuildCurves(side);
    }

    // Every curve-form knob on one side, summed (cascaded linear filters multiply, which in dB is a
    // sum) on one common grid — each knob may ship its own — and designed as ONE minimum-phase FIR.
    // Flat, or nothing there, is a bypass rather than a convolution with an impulse.
    void rebuildCurves(int side) {
        std::vector<double> sum(commonGrid_.size(), 0.0);
        bool any = false;
        for (const auto& t : tones()) {
            // A knob whose positions carry BANDS ships no curve and no grid: it belongs to the other
            // path, and asking this one for its curve would ask an empty grid for a value.
            if (sideOf(t) != side || t.positions.empty() || bandsPerPosition(t)) continue;
            const auto grid  = gridOf(t);
            const auto curve = curveAt(t, grid, toneAt_[t.name]);
            if (curve.empty()) continue;
            any = true;
            for (std::size_t k = 0; k < sum.size(); ++k)
                sum[k] += felitronics::lineareq::curveDbAt(curve, grid, commonGrid_[k]);
        }
        curveSum_[side] = any ? sum : std::vector<double> {};
        auto taps = any ? felitronics::lineareq::magnitudeCurveToFir(sum, commonGrid_, fs_, kFirTaps, kFirDesign)
                        : std::vector<float> {};
        firTaps_[side] = std::move(taps);
        if (! firTaps_[side].empty() && prepared_) {
            const float* t[1] { firTaps_[side].data() };
            fir_[side].loadIR(t, 1, (int) firTaps_[side].size(), fs_);
            firBypass_[side].store(false, std::memory_order_release);
        } else
            firBypass_[side].store(true, std::memory_order_release);
    }

    // Every band of every section-form knob on one side, at its knob's position — in a fixed order,
    // so band k keeps its state across knob moves.
    void rebuildBands(int side) {
        BandSet set;
        bandsShown_[side].clear();
        int dropped = 0;
        for (const auto& t : tones()) {
            if (sideOf(t) != side) continue;
            // A DIAL's bands travel with its rotation; a SWITCH's stand whole at the position it is on.
            // One list either way, in a fixed order, so band k keeps its state across knob moves.
            const auto qs = ! t.sections.empty() ? sectionsAt(t, toneAt_[t.name], fs_)
                          : bandsPerPosition(t)  ? sectionsAtValue(t, toneAt_[t.name], fs_)
                                                 : std::vector<felitronics::rigplayer::SectionBiquad> {};
            for (const auto& q : qs) {
                if (set.count < kMaxBands) { set.c[set.count++] = q; bandsShown_[side].push_back(q); }
                else ++dropped;
            }
        }
        bandsDropped_ = dropped;
        publishBands(side, set);
    }

    // The first blend knob of the stage — a pedal has one. Its dry path's shape is designed once; the
    // gains move with the knob and never rebuild a filter.
    const namz::rig::Blend* blendKnob() const { return stage_.blend.empty() ? nullptr : &stage_.blend.front(); }

    void rebuildDry() {
        const auto* b = blendKnob();
        dryTaps_.clear();
        if (b != nullptr) {
            const auto grid = gridOf(*b);
            dryTaps_ = felitronics::lineareq::magnitudeCurveToFir(dryCurve(*b, grid), grid, fs_, kFirTaps, kFirDesign);
        }
        if (! dryTaps_.empty() && prepared_) {
            const float* t[1] { dryTaps_.data() };
            dry_.loadIR(t, 1, (int) dryTaps_.size(), fs_);
            dryFirBypass_.store(false, std::memory_order_release);
        } else
            dryFirBypass_.store(true, std::memory_order_release);
        applyBlendGains();
        dryActive_.store(b != nullptr, std::memory_order_release);
    }

    void applyBlendGains() {
        const auto* b = blendKnob();
        const BlendGains g = b != nullptr ? blendGainsAt(*b, blendAt_[b->name]) : BlendGains {};
        dryGain_.store((float) std::clamp(g.dry, -4.0, 4.0), std::memory_order_release);
        wetGain_.store((float) std::clamp(g.wet,  0.0, 4.0), std::memory_order_release);
    }

    // ------------------------------------------------------------ the bands, across the threads ----

    struct BandSet {
        int count = 0;
        felitronics::rigplayer::SectionBiquad c[kMaxBands];
    };
    // A small pool instead of a heap: the message thread fills a free set and publishes its index; the
    // audio thread takes the index, copies the set and frees it. At most one pending and one being
    // filled at a time, so four never run out — and nothing is allocated or freed on the audio thread.
    struct BandSlot { BandSet set; std::atomic<bool> busy { false }; };
    static constexpr int kBandPool = 4;

    void publishBands(int side, const BandSet& set) {
        auto& pool = bandPool_[side];
        for (int i = 0; i < kBandPool; ++i) {
            if (pool[i].busy.load(std::memory_order_acquire)) continue;
            pool[i].set = set;
            pool[i].busy.store(true, std::memory_order_release);
            const int old = bandNext_[side].exchange(i, std::memory_order_acq_rel);
            if (old >= 0) pool[old].busy.store(false, std::memory_order_release);   // never taken: superseded
            return;
        }
    }

    // Audio thread, at the top of every block.
    void takeBands(int side) {
        const int i = bandNext_[side].exchange(-1, std::memory_order_acq_rel);
        if (i < 0) return;
        auto& rt = bandRt_[side];
        const auto& in = bandPool_[side][i].set;
        // A band that appears (a load) starts AT its coefficients — there is nothing to ramp from; one
        // that is already running ramps to the new ones across the next block.
        for (int k = rt.count; k < in.count; ++k) {
            bandCur_[side][k] = bandFrom_[side][k] = bandTo_[side][k] = in.c[k];
            bandPos_[side][k] = bandLen_[side][k] = 0;         // arrived: nothing to travel
            for (auto& b : bq_[side][k]) b.reset();
        }
        rt = in;
        bandPool_[side][i].busy.store(false, std::memory_order_release);
    }

    // The same five numbers, bit for bit: a band that did not move ramps nothing.
    static bool same(const felitronics::rigplayer::SectionBiquad& x, const felitronics::rigplayer::SectionBiquad& y) {
        return std::memcmp(&x.b0, &y.b0, sizeof(double)) == 0 && std::memcmp(&x.b1, &y.b1, sizeof(double)) == 0
            && std::memcmp(&x.b2, &y.b2, sizeof(double)) == 0 && std::memcmp(&x.a1, &y.a1, sizeof(double)) == 0
            && std::memcmp(&x.a2, &y.a2, sizeof(double)) == 0;
    }

    // Every band of a side, on every channel. Coefficients that changed travel linearly over a length
    // THE BAND ITSELF DECLARES (see bandRampLength) — a knob dragged sixty times a second never steps a
    // filter, and it now does not step it differently on a different host either. The ramp carries
    // across calls: `bandPos_` counts SAMPLES, so where the caller cuts the stream changes nothing.
    //
    // THE SEGMENT LOOP IS ON THE OUTSIDE, and the first version of this had it on the inside — all the
    // audio processed first, then a loop that counted grid boundaries and flushed at each. That flushes
    // the same FINAL state eight times for a 512-sample call and never once inside a whole-file call,
    // which is the very defect the grid exists to remove, dressed as its fix. Found by the diff-pass
    // consilium; the suite could not see it, because a mutation that removes the flush entirely fails
    // while one that merely misplaces it does not.
    void runBands(int side, float* const* planes, int count) {
        auto& rt = bandRt_[side];

        // A NEW TARGET STARTS A NEW RAMP, from wherever the last one had got to — a hand that moves
        // again mid-travel is the ordinary case, not the exception. This is an ARRIVAL, so it is read
        // once per call, at the call boundary, exactly where the caller put it.
        for (int k = 0; k < rt.count; ++k)
            if (! same(rt.c[k], bandTo_[side][k])) {
                bandFrom_[side][k] = bandCur_[side][k];
                bandTo_[side][k]   = rt.c[k];
                bandLen_[side][k]  = bandRampLength(bandFrom_[side][k], rt.c[k], fs_);
                bandPos_[side][k]  = 0;
            }

        auto& grid = bandGrid_[side];
        for (int off = 0; off < count; ) {
            const int seg = grid.segment(count - off);
            for (int k = 0; k < rt.count; ++k) runBandSegment(side, k, planes, off, seg);
            off += seg;
            // LAW 8 ON THE AUDIO-TIME GRID (core/StateGrid.h, law 8a), which means HERE — after the
            // segment that ends on a boundary and before the next one, not after the call.
            if (grid.advance(seg))
                for (int k = 0; k < rt.count; ++k)
                    for (int c = 0; c < channels_; ++c) bq_[side][k][(std::size_t) c].flushDenormals();
        }

        // The poison half keeps the call's clock — see eq::Biquad::healPoison(). Invisible while finite.
        for (int k = 0; k < rt.count; ++k)
            for (int c = 0; c < channels_; ++c) bq_[side][k][(std::size_t) c].healPoison();
    }

    // One band, one grid segment: `n` samples starting at `off`, every channel, walking the ramp.
    void runBandSegment(int side, int k, float* const* planes, int off, int n) {
        auto& bands = bq_[side][k];
        if (bandPos_[side][k] >= bandLen_[side][k]) {           // arrived: no interpolation at all
            bandCur_[side][k] = bandTo_[side][k];
            for (int c = 0; c < channels_; ++c) {
                auto& bq = bands[(std::size_t) c];
                bq.c = bandTo_[side][k];
                float* x = planes[c] + off;
                for (int i = 0; i < n; ++i) x[i] = bq.processSample(x[i]);
            }
            return;
        }
        const auto& from = bandFrom_[side][k];
        const auto& to   = bandTo_[side][k];
        const int    len = bandLen_[side][k];
        const double inv = 1.0 / (double) len;
        const int    p0  = bandPos_[side][k];
        for (int c = 0; c < channels_; ++c) {
            auto& bq = bands[(std::size_t) c];
            float* x = planes[c] + off;
            int p = p0;
            for (int i = 0; i < n; ++i) {
                if (p < len) {
                    ++p;
                    const double t = (double) p * inv;
                    bq.c.b0 = from.b0 + (to.b0 - from.b0) * t;
                    bq.c.b1 = from.b1 + (to.b1 - from.b1) * t;
                    bq.c.b2 = from.b2 + (to.b2 - from.b2) * t;
                    bq.c.a1 = from.a1 + (to.a1 - from.a1) * t;
                    bq.c.a2 = from.a2 + (to.a2 - from.a2) * t;
                } else bq.c = to;
                x[i] = bq.processSample(x[i]);
            }
            bandCur_[side][k] = bq.c;                           // every channel walks the same path
        }
        bandPos_[side][k] = std::min(p0 + n, len);
    }

    // Law 8. Every gain here approaches its target asymptotically — `end = want + (current-want)*decay`
    // — so it never arrives: the residual parks in a fixed-point band and stays. For want == 0 that band
    // is subnormal (at 128 samples / 48 kHz decay = 0.766, so every k*u with k <= 0.5/(1-decay) ~= 2.1
    // maps to itself; measured: frozen at 2 ulp, stall beginning ~1 s in), and the audio then pays for it
    // per sample: `gd += stepDry; a[c][i]*gw + d[c][i]*gd` is ~3 assisting operations per sample per
    // channel, and the mix loop is gated on `dryActive_` — a dry IR being LOADED — not on either gain.
    //
    // An EXACT zero is reachable and the repo's own fixture ships one: `BlendKnob::linOf` returns 0.0 for
    // any level at or below -120 dB, which is how a pack spells "this path is off" (RigPlayerTests.cpp's
    // rig: wetDb -120 at the dry end, dryDb -120 at the wet end). What does NOT stall is the initialised
    // or reset 0 — 0 -> 0 stays exactly 0; the stall needs an audible position first and a -120 dB one
    // after it. The second condition below covers the non-zero targets the same way core::Smoother does.
    static float rampEnd(float want, float current, float decay) {
        const float d  = (current - want) * decay;
        const float nx = want + d;
        return (std::fabs(d) < 1e-30f || felitronics::core::exactlyEqual(nx, current)) ? want : nx;
    }

    // One linear ramp per block toward an exponential ~10 ms endpoint: smooth moves without an exp()
    // per sample. `current` is audio-thread-only.
    static bool isOne(float g) { const float one = 1.0f; return std::memcmp(&g, &one, sizeof(float)) == 0; }   // bit for bit, as same()
    void rampInto(float* const* planes, int nch, int count, float want, float& current) const {
        const float decay = std::exp(-(float) count / (float) (0.010 * fs_));
        const float end   = rampEnd(want, current, decay);
        const float step  = (end - current) / (float) count;
        if (! isOne(current) || ! isOne(want))
            for (int c = 0; c < nch; ++c) {
                float g = current;
                float* x = planes[c];
                for (int i = 0; i < count; ++i) { g += step; x[i] *= g; }
            }
        current = end;
    }

    // ------------------------------------------------------------------------------------ state ----

    // message thread
    struct Model { std::string id; std::shared_ptr<const std::vector<std::byte>> bytes; };   // null = not fetched yet
    std::vector<Model> models_;                        // BlendModelId = index + 1; 0 = nothing
    std::uint64_t      gen_ = 0;                       // bumped by every unload(): a job carries the one it was taken under
    bool               jobOut_ = false;                // a job is with the host; nothing else is handed out
    std::vector<int>   fileModel_;                     // Device::files index → models_ index
    namz::rig::Stage   stage_;
    ModelSource        source_;
    bool               loaded_ = false;
    Settings           axes_;
    std::string        dial_;
    int                dialSweep_ = 0;
    double             dialDeg_ = 0.0;
    std::map<std::string, std::string> toneAt_, blendAt_;
    std::vector<namz::rig::Tone> toneOverride_;        // tone handed in beside the manifest; empty = none
    Selection          sel_;
    SlotPlan           plan_;
    BlendShape         shape_ {};
    double             topExtend_ = kTopExtendDbPerDeg;
    AlignmentTable     align_;
    long long          loads_ = 0;
    int                bandsDropped_ = 0;
    double             coldSeconds_ = kColdAfterSeconds;
    // A 1/12-octave grid from 20 Hz to 20 kHz, on which every curve-form knob is summed whatever grid
    // it shipped on.
    std::vector<double> commonGrid_ = felitronics::lineareq::logFreqGrid(20.0, 20000.0, 121);
    std::vector<double> curveSum_[2];
    std::vector<float>  firTaps_[2], dryTaps_;
    std::vector<felitronics::rigplayer::SectionBiquad> bandsShown_[2];

    // both, by contract
    double fs_ = 48000.0;
    int    maxBlock_ = 0, channels_ = 1;
    bool   prepared_ = false;

    // the mailboxes (message → audio, audio → message)
    felitronics::nam::NamStage nam_[2];
    std::atomic<std::uint64_t> reqWant_[2] { 0, 0 };
    std::atomic<float>         reqTarget_ { 0.0f };
    std::atomic<std::uint64_t> loadAsk_ { 0 };          // (model << 1) | slot; 0 = nothing asked
    std::atomic<bool>          landFlag_ { false };
    std::atomic<int>           landSlot_ { 0 };
    std::atomic<std::uint64_t> landModel_ { 0 };
    std::atomic<long long>     landWarm_ { 0 };
    std::atomic<int>           landDelay_ { 0 };
    std::atomic<int>           landFail_ { 0 };          // slot + 1; 0 = none
    std::atomic<bool>          forget_ { false };
    std::atomic<int>           pendDelay_[2] { 0, 0 };
    std::atomic<int>           slotDelay_[2] { 0, 0 };
    std::atomic<std::uint64_t> held_[2] { 0, 0 };
    std::atomic<float>         liveMix_ { 0.0f };
    std::atomic<float>         liveDry_ { 0.0f }, liveWet_ { 1.0f };
    std::atomic<int>           warmBlocks_ { 0 }, mixJumps_ { 0 };
    std::atomic<float>         biggestJump_ { 0.0f };
    std::atomic<long long>     coldAfter_ { 0 };          // coldSeconds_ at fs_, for the law; 0 = never
    std::atomic<bool>          slotCold_[2] { false, false };
    std::atomic<int>           coldBlocks_[2] { 0, 0 };
    // The pack's own levels, and the reason they are two gains of their own rather than one more
    // term in chainGain_: inGain_ has to reach the dry side of a blend, and chainGain_ is applied
    // after the dry block has already been taken.
    std::atomic<float>         inGain_ { 1.0f }, outGain_ { 1.0f };
    std::optional<double>      hostIn_, hostOut_;      // message thread only; empty = the pack's own
    std::atomic<float>         chainGain_ { 1.0f };
    std::atomic<float>         slotGain_[2] { 1.0f, 1.0f };
    // A model's `metadata.loudness` tag is a CONTRACT, not a listener's option: off, every capture
    // plays at whatever level the hardware happened to give, and two packs cannot be compared at all.
    std::atomic<bool>          normalize_ { true };
    std::atomic<bool>          inputTrims_ { true };
    std::atomic<bool>          firBypass_[2] { true, true };
    std::atomic<bool>          dryActive_ { false }, dryFirBypass_ { true };
    std::atomic<float>         dryGain_ { 0.0f }, wetGain_ { 1.0f };
    BandSlot                   bandPool_[2][kBandPool];
    std::atomic<int>           bandNext_[2] { -1, -1 };

    // audio thread only
    felitronics::nam::BlendState blend_ {};
    bool coldRt_[2] {};                                // the law's cold flags as of the last block, to see a slot fall asleep
    std::array<float, (std::size_t) kMaxDelay> lagTail_[2][kMaxChannels] {};
    int ranNch_ = 0;                        // planes that advanced state on the previous call
    float curIn_ = 1.0f, curOut_ = 1.0f;
    float curChain_ = 1.0f, curSlot_[2] { 1.0f, 1.0f }, curDry_ = 0.0f, curWet_ = 1.0f;
    BandSet                    bandRt_[2];
    felitronics::rigplayer::SectionBiquad   bandCur_[2][kMaxBands];    // what the filters are running NOW
    felitronics::rigplayer::SectionBiquad   bandFrom_[2][kMaxBands];   // where the ramp in flight started
    felitronics::rigplayer::SectionBiquad   bandTo_[2][kMaxBands];     // ...and where it is going
    int                                     bandPos_[2][kMaxBands] {}; // samples of it consumed
    int                                     bandLen_[2][kMaxBands] {}; // its declared length (bandRampLength)
    felitronics::core::StateGrid            bandGrid_[2];              // law 8 on audio time, per side
    felitronics::eq::Biquad    bq_[2][kMaxBands][kMaxChannels];
    felitronics::convolution::CabConvolver fir_[2], dry_;
    felitronics::core::DryAligner dryLatency_;   // holds the dry path back by the models' rate-match
    std::vector<float> slotB_[kMaxChannels], dryBuf_[kMaxChannels], spare_[kMaxChannels];
};

} // namespace felitronics::rigplayer
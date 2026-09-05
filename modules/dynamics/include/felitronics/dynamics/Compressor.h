// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/core/Config.h>
#include <felitronics/core/Math.h>
#include <felitronics/core/DelayLine.h>
#include <felitronics/core/FlushToZero.h>
#include <felitronics/dynamics/GainComputer.h>
#include <felitronics/dynamics/GainReductionFollower.h>
#include <felitronics/dynamics/ChannelLinker.h>
#include <felitronics/dynamics/EnvelopeFollower.h>   // Detector enum (Peak / Rms)
#include <felitronics/dynamics/LinkedDetector.h>     // DetectorParams + the detector path itself

#include <cmath>
#include <vector>

namespace felitronics::dynamics
{

// Derives from DetectorParams (`detector`, `link`, `rmsWindowMs`) rather than repeating those three
// fields, so that handing a CompressorParams to a LinkedDetector is a base-class binding and not a
// mapping somebody has to keep in step. Field access is unchanged — `p.detector`, `p.link`,
// `p.rmsWindowMs` all still name the same things. It stays an aggregate (a public, non-virtual base is
// allowed in one since C++17), so `CompressorParams p;` then assigning members works exactly as before,
// and it stays trivially copyable. THIS IS A SOURCE-COMPATIBILITY BREAK, and the header is the place to
// say so plainly rather than to reassure: POSITIONAL brace initialisation shifts by three fields, a
// a designated initialiser naming a base member (`{.detector = ...}`) no longer compiles, and
// `std::is_standard_layout_v` becomes false (`is_aggregate_v` and `is_trivially_copyable_v` stay true,
// and positional brace initialisation still works both flat and with the base nested — compile-checked,
// both spellings). Nothing in this repository does any of those; a consumer outside it might.
struct CompressorParams : DetectorParams
{
    Mode mode = Mode::DownCompress;

    double thresholdDb = -18.0;
    double ratio       = 2.0;
    double kneeDb      = 6.0;
    double rangeDb     = 60.0;

    double attackMs    = 10.0;                    // ballistics on the GAIN REDUCTION (not the level)
    double releaseMs   = 100.0;

    double makeupDb    = 0.0;
    bool   autoMakeup  = false;                   // adds a STATIC auto-makeup (see autoMakeupDb)
    double lookaheadMs = 0.0;
};

// A STATIC auto-makeup (dB): half-compensate the gain reduction a 0 dBFS signal would receive. Static
// (a function of the curve only) — NOT signal-following (which would be a second compressor and pump).
// It reads the CURVE, so it assumes the key and the programme share a scale; an external key at a
// different level makes it over- or under-correct, and the product owns that decision.
inline double autoMakeupDb (const GainComputer& gc) noexcept { return -0.5 * gc.deltaDb (0.0); }

//==============================================================================
// felitronics::dynamics::Compressor — a broadband, channel-LINKED compressor that composes the dynamics
// primitives in the clean topology:
//     LinkedDetector (gate → link → envelope) → stateless GainComputer (curve) → GR ballistics (timing)
//     → makeup → apply to the LOOKAHEAD-delayed signal.
// Smoothing the gain reduction (not the detector level) keeps the attack/release independent of the
// knee/ratio. Product-neutral: no sidechain EQ, params system, GUI, metering policy, dry/wet — those
// stay in the product (it's the "core gives primitives + a broadband composite; products compose" line).
// RT-safe: prepare() allocates; process() does no alloc/lock/throw. Reported latency = lookahead.
//
// ====================================================================================
// THE EXTERNAL KEY (sidechain), and the contract the caller owes it
// ====================================================================================
// `process()` has a five-argument form that feeds the DETECTOR from a separate buffer while the gain is
// still applied to `io`. That is how a mastering bus stops its kick and bass from ducking the whole mix:
// the product high-passes a copy of the bus and passes it as the key. THE FILTER IS NOT IN HERE and
// will not be — ADR §4, "the core gives primitives; products compose them" — so the key arrives already
// filtered, built from `eq::MatchedBiquad` or anything else the product likes.
//
//   * TIME BASE. `key[c][i]` must be the SAME INSTANT as `io[c][i]` at the input. The compressor
//     supplies the lookahead itself, by delaying `io` and not the key, so key sample `i` sets the gain
//     applied to programme sample `i - latencySamples()`. Do NOT pre-advance or pre-delay the key by
//     this module's latency: a host compensates a plug-in's reported latency downstream of it, so both
//     buses arrive on the same timeline. (Hosts have not always got the key bus right — Pro Tools
//     before 2021.6 and Logic before 10.8 did not compensate it — but that is a wrapper's problem to
//     detect, not something a DSP primitive can guess.)
//   * KEY-FILTER GROUP DELAY eats lookahead. A minimum-phase HPF is not a constant delay, and it is
//     not free either: a 2nd-order Butterworth high-pass at fc has a group delay of sqrt(2)/(2*pi*fc)
//     AT the corner — 0.75 ms, or 36 samples at 48 kHz, for a 300 Hz corner — falling as 1/f^2 above
//     it, which is why it is normally left uncompensated for the transients a key exists to pass. How
//     much is a property of the filter the product chose, not something this module knows. A
//     LINEAR-PHASE key filter is a different
//     matter: an N-tap symmetric FIR delays the key by exactly (N-1)/2 samples, which subtracts that
//     much from the effective lookahead and can make it negative. Use a minimum-phase key filter, or
//     delay `io` to match and declare the extra latency.
//   * CHANNEL COUNT is the key's own. `numKeyChannels` has nothing to do with `numChannels`: one mono
//     kick bus keying a stereo master is `numKeyChannels == 1`, and it is read for every sample of the
//     block. Only `key[0 .. numKeyChannels-1]` are touched, and the count must be the true length of
//     that array — an array length is not observable from in here. NB `LinkMode::MeanPower` divides by
//     the key's channel count, so a mono key and the same key duplicated into two channels agree, while
//     a mono key and a stereo key whose second channel is silent differ by 3 dB. That is arithmetic, not
//     a bug, and it is the reason the count travels with the pointer.
//   * `key == nullptr`, OR A COUNT OF ZERO OR LESS, is the SELF-KEYED path, identical to the
//     three-argument form. Passing `io` itself as the key is also exactly identical, bit for bit: at
//     each sample every detector read happens before any write. NB FOR A PLUG-IN WRAPPER: a host that
//     reports zero sidechain channels because the bus is not connected therefore gets SELF-KEYING, not
//     silence — which is the safe default, but it is not "external sidechain engaged and receiving
//     nothing". A product that wants the second must say so itself.
//   * THERE IS NO FOUR-ARGUMENT FORM, deliberately. A bare pointer with no count is the shape that
//     reads one element past a mono key array on a stereo bus, silently; here it does not compile.
//   * ALIASING IS EXACT OR NOT AT ALL. `key[c] == io[c]` is the self-keyed identity above. A key that
//     OVERLAPS the programme at a shift (`io[0] - 1`, say) reads samples this call has already
//     overwritten, and two `io` planes pointing at the same buffer get the gain and the delay applied
//     twice. Neither is diagnosable from in here; both are the caller's to avoid.
//   * A POISONED KEY FILTER IS UPSTREAM OF THE GATE. The gate bounds what reaches the detector, not
//     what reaches the product's own filter: one Inf into a recursive HPF makes every sample it emits
//     afterwards a NaN, which this module then reads — correctly — as silence, for good. A product that
//     filters the key owns guarding that filter's input and state.
//
// ====================================================================================
// THE DETECTOR GATE — what it protects, and what it deliberately does not touch
// ====================================================================================
// Every key sample passes `detectorGate()` (non-finite → 0, magnitude clamped to 1e6) BEFORE the link.
// This is not decoration: the detector is recursive, and `state += c*(x - state)` never decays a NaN
// out, so ONE bad key sample used to be permanent. Measured on the unguarded code, from a single +Inf
// (or a NaN, or a finite 1e20 — which squares to +Inf in float): DownCompress became a permanent
// bypass, DownExpand a permanent -rangeDb of silence, and UpCompress a permanent +60 dB. The gate is
// per SAMPLE and not on the linked result because `MeanPower` squares inside the link, so a gate on the
// result would see the overflowed frame as silence — a loud key read as nothing.
// BIT-TRANSPARENT for any finite key sample within ±1e6 (+120 dBFS): it takes neither branch, so every
// real audio sample passes through untouched with 120 dB to spare. ABOVE that bound the two detectors
// part company, and it is worth being exact about it. `Peak` is still output-transparent once the curve
// has saturated, because the level is instantaneous and both the clamped and unclamped spikes land on
// -rangeDb. `Rms` is NOT: the window stores x^2/tau and carries it forward, so what the clamp bounds is
// the MEMORY of the spike, not just its instant. Measured with one key sample of 1e8, a 5 ms window and
// the defaults: gated and ungated gain reduction differ in 46010 of 48000 samples, by up to 16.65 dB at
// 4:1 and 17.45 dB at 20:1, where `Peak` differs in none. `MeanPower` behaves like `Rms`.
// THE AUDIO PATH IS NOT GATED. A compressor is not a sanitiser: `io` carries no recursive state (the
// delay line is a ring, not a feedback path), so a non-finite programme sample passes through and out,
// which is the honest answer and leaves bounding the signal to the limiter that owns that promise.
//
// ====================================================================================
// WHAT CHANGES OF SHAPE DO, all of them deliberate
// ====================================================================================
//   * MORE CHANNELS THAN PREPARED are REFUSED — the buffer is returned untouched. This module's one
//     structural promise is that every channel gets the SAME gain; processing a prefix breaks exactly
//     that, and silently: measured with a stereo buffer on a mono-prepared instance, the second channel
//     came out neither gained nor delayed, i.e. a level mismatch plus a lookahead-sized time skew — a
//     comb-filtered fold-down that can ship unnoticed. Delivering no compression cannot. Prepare for the
//     most channels you will ever pass, and check what prepare() returned.
//   * A CHANGE OF CHANNEL COUNT resets the delay lines of the channels that were idle, and nothing
//     else. Only `delays[]` is per-channel; the detector and the gain reduction are shared and keep
//     advancing, so they are not stale and clearing them would be a gain jump the change never asked
//     for. Without this, stereo → mono → stereo re-emitted the returning channel's pre-mono audio:
//     measured at 0.75 out of digital silence, 417 ms after the fact.
//   * A CHANGE OF THE RESULTING SAMPLE DELAY — the rounded, clamped one, not every touch of
//     `lookaheadMs` — clears the delay lines and moves `latencySamples()`, which makes it a
//     resynchronisation event rather than automation. Moving the read pointer of a live ring re-emits
//     audio already delivered (measured: 48 samples of 0.9 into silence) or skips over some; a cleared
//     line emits the gap as zeros instead, which is the bounded failure. Set it once if you can.
class Compressor
{
public:
    // The bounds prepare() enforces, public so a wrapper that computes sizes from the SAME arguments
    // can reject an out-of-range value before its own arithmetic rather than after ours.
    static constexpr double kMaxSampleRate  = 3.0e6;    // far above any audio rate; keeps every derived
                                                        // size inside int with the lookahead capacity
    static constexpr double kMaxLookaheadMs = 250.0;    // five times any musical setting; bounds what one
                                                        // prepare() can be asked to allocate per channel

    // Returns false and leaves the compressor UNPREPARED on a configuration it cannot honour; process()
    // then does nothing at all rather than half-processing. A false return is the only way to learn
    // this, so check it. What it prevents is not `std::ceil(NaN)` — that just returns NaN — but the
    // out-of-range floating-to-INTEGER conversions downstream of it, which are undefined.
    [[nodiscard]] bool prepare (double sampleRate, int maxBlock, int maxChannels, double maxLookaheadMs = 50.0)
    {
        prepared_ = false;                                     // any early return below leaves it unprepared
        // Spelled positively so NaN fails: `sampleRate <= 0.0` is FALSE for a NaN, which is how one
        // reaches the undefined conversions above.
        if (! (sampleRate > 0.0) || ! std::isfinite (sampleRate) || sampleRate > kMaxSampleRate) return false;
        if (maxBlock < 1) return false;                        // no scratch is sized by it, but a caller
                                                               // that passes 0 has a bug worth reporting
        // Refused rather than clamped: a silently reduced channel count would leave the surplus channels
        // refused at every process() call, i.e. a compressor that does nothing, discovered much later.
        if (maxChannels < 1 || maxChannels > core::kMaxChannels) return false;
        if (! std::isfinite (maxLookaheadMs) || maxLookaheadMs < 0.0 || maxLookaheadMs > kMaxLookaheadMs) return false;

        fs    = sampleRate;
        maxCh = maxChannels;
        maxLookSamples = (int) std::ceil (maxLookaheadMs * 0.001 * fs);
        delays.assign ((std::size_t) maxCh, core::DelayLine {});
        for (auto& d : delays) d.prepare (maxLookSamples);
        det.prepare (fs);
        grFollower.prepare (fs);
        lookSamples = -1;                                      // force apply() to size the fresh lines
        apply (params);
        reset();
        prepared_ = true;
        return true;
    }

    void reset() noexcept
    {
        for (auto& d : delays) d.reset();
        det.reset();
        grFollower.reset();
        lastNc_ = 0;
    }

    void   setParams (const CompressorParams& p) noexcept { params = p; apply (p); }
    int    latencySamples()  const noexcept { return prepared_ ? lookSamples : 0; }
    double gainReductionDb() const noexcept { return grFollower.valueDb(); }   // for metering
    float  detectorLevel()   const noexcept { return det.level(); }            // linked level, linear
    bool   isPrepared()      const noexcept { return prepared_; }

    // Audio thread, in place, self-keyed (the detector reads the programme). RT-safe.
    void process (float* const* channels, int numChannels, int numSamples) noexcept
    {
        process (channels, numChannels, numSamples, nullptr, 0);
    }

    // Audio thread, in place, with an EXTERNAL KEY — read the header note above for the time-base and
    // channel-count contract. `key == nullptr` or `numKeyChannels <= 0` means self-keyed. RT-safe.
    void process (float* const* channels, int numChannels, int numSamples,
                  const float* const* key, int numKeyChannels) noexcept
    {
        if (! prepared_ || numChannels <= 0 || numSamples <= 0) return;
        if (numChannels > maxCh) return;                       // refuse — see "MORE CHANNELS THAN PREPARED"
        const int nc = numChannels;

        // A channel that sat out blocks holds `lookSamples` of audio from before it left. Zero only
        // those lines; the shared detector and gain reduction are still tracking a real signal.
        if (lastNc_ > 0 && nc > lastNc_)
            for (int c = lastNc_; c < nc; ++c) delays[(std::size_t) c].reset();
        lastNc_ = nc;

        const float* const* detCh = key;
        int                 detNc = numKeyChannels;
        if (key == nullptr || numKeyChannels <= 0) { detCh = channels; detNc = nc; }

        for (int i = 0; i < numSamples; ++i)
        {
            // --- detector: gate → link → envelope, in the object P3 can build from the same params ---
            const float level = det.process (detCh, detNc, i);

            // --- curve → GR ballistics → gain ---
            const double levelDb    = core::gainToDb (level);
            const float  targetDb   = (float) gc.deltaDb (levelDb);
            const float  grDb       = grFollower.process (targetDb);
            // THE SUM is what reaches dbToGain, and it is bounded here and nowhere else. An UpCompress
            // curve sitting at +rangeDb plus a large makeup is +800 dB, i.e. 1e40 — +Inf once cast to
            // float — and `0.0f * Inf` is NaN, so a silent passage came out poisoned from entirely
            // finite parameters (measured: 64 of 64 output samples). Two comparisons, and
            // bit-transparent for any total inside ±400 dB, a gain of 1e±20 no real setting approaches.
            const double totalDb    = (double) grDb + makeupAppliedDb;
            const float  gain       = (float) core::dbToGain (totalDb < -kMaxGainDb ? -kMaxGainDb : (totalDb > kMaxGainDb ? kMaxGainDb : totalDb));

            // --- apply to the lookahead-delayed signal (read input before overwriting in place) ---
            for (int c = 0; c < nc; ++c)
            {
                const float x = channels[c][i];
                const float y = delays[(std::size_t) c].process (x);
                channels[c][i] = y * gain;
            }
        }

        det.flushDenormals();
        grFollower.flushDenormals();
    }

private:
    // With kMaxSampleRate and kMaxLookaheadMs the delay ring is at most 750001 slots, and DelayLine
    // casts its own size to int internally — so the pair above is what keeps that conversion in range.
    static constexpr double kMaxGainDb      = 400.0;    // gain reduction PLUS makeup, the only thing that
                                                        // reaches dbToGain; 400 dB is a gain of 1e20

    void apply (const CompressorParams& p) noexcept
    {
        gc.setMode (p.mode);
        gc.setThresholdDb (p.thresholdDb);
        gc.setRatio (p.ratio);
        gc.setKneeDb (p.kneeDb);
        gc.setRangeDb (p.rangeDb);
        grFollower.setTimes (p.attackMs, p.releaseMs);   // non-finite times → instant (the follower's coeff guard)

        // SLICING, on purpose: the compressor's own parameter object IS the detector's. Nothing maps
        // fields across, so nothing can fall out of step when a detector field is added. Non-finite
        // windows are handled inside the follower's coefficient guard, same as the times above.
        det.setParams (p);

        // A non-finite lookahead would feed lround() undefined behaviour; clamp in the DOUBLE domain
        // first, because std::lround(1e300) is out of range for a long before any int clamp can help.
        const double lookMs = std::isfinite (p.lookaheadMs) ? (p.lookaheadMs < 0.0 ? 0.0 : (p.lookaheadMs > kMaxLookaheadMs ? kMaxLookaheadMs : p.lookaheadMs)) : 0.0;
        int newLook = (int) std::lround (lookMs * 0.001 * fs);
        if (newLook < 0) newLook = 0;
        if (newLook > maxLookSamples) newLook = maxLookSamples;

        // Changing the delay of a LIVE ring moves only the read pointer: raising it re-emits audio
        // already delivered, lowering it skips over some. Clear instead — a hole of zeros is the
        // bounded failure, and the reported latency has moved anyway, which is a resync either way.
        const bool lookChanged = (newLook != lookSamples);
        lookSamples = newLook;
        for (auto& d : delays) d.setDelay (lookSamples);
        if (lookChanged) for (auto& d : delays) d.reset();

        // Not bounded HERE, deliberately: what has to stay inside float is the SUM of this and the
        // gain reduction, and process() bounds that in one place. Two clamps would mask each other —
        // and did, until a mutation of either survived the suite because the other was still standing.
        makeupAppliedDb = (std::isfinite (p.makeupDb) ? p.makeupDb : 0.0) + (p.autoMakeup ? autoMakeupDb (gc) : 0.0);
    }

    double fs = 48000.0;
    CompressorParams params;

    GainComputer          gc;
    LinkedDetector        det;
    GainReductionFollower grFollower;
    std::vector<core::DelayLine> delays;

    int    lookSamples = 0, maxLookSamples = 0, maxCh = 0, lastNc_ = 0;
    double makeupAppliedDb = 0.0;
    bool   prepared_ = false;
};

} // namespace felitronics::dynamics

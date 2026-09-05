// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/core/Config.h>
#include <felitronics/core/DryAligner.h>
#include <felitronics/core/Math.h>
#include <felitronics/dither/Dither.h>
#include <felitronics/dynamics/Compressor.h>
#include <felitronics/eq/EqEngine.h>
#include <felitronics/eq/MatchedBiquad.h>
#include <felitronics/limiter/TruePeakLimiter.h>
#include <felitronics/saturation/Saturator.h>
#include <felitronics/stereo/MonoBass.h>

#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

namespace felitronics::mastering
{

//==============================================================================
// TOPOLOGY, fixed for the life of a prepared stream — which stages EXIST, and the three numbers that
// decide how much latency they cost. None of these can be a per-block parameter: each one moves
// latencySamples(), and a moving latency is a host resynchronisation event rather than automation.
// This is the split `TruePeakLimiterConfig` introduced in P12, applied one level up.
//
// A stage that is absent here is not allocated, not called, and not in the latency sum. A stage that
// is PRESENT can still be bypassed per block (MasteringChainParams), and that never moves latency.
struct MasteringChainConfig
{
    // THE INTERNAL QUANTUM — the load-bearing decision of this whole module, and the only reason
    // criterion "same input, same output, whatever the caller's block size" is a THEOREM here rather
    // than a hope. The chain buffers, and every stage is always called with exactly this many samples,
    // no matter how the caller cuts the stream. See the block-invariance note below for the measured
    // reason that is necessary. Costs exactly `internalBlock` samples of latency, and — measured —
    // no CPU: the whole chain runs at 2.26-2.30 %RT flat for every K from 16 to 4096.
    int internalBlock = 256;

    // Which stages exist. gain nodes are free and always present.
    bool eq         = true;
    bool monoBass   = false;      // OFF by default: it is a decision about the low end, not a default
    bool compressor = true;
    bool clipper    = false;      // OFF by default: it is where loudness is won, and it costs latency
    bool limiter    = true;
    bool dither     = true;

    // Latency-bearing topology of the stages above.
    double compressorLookaheadMs = 1.0;
    double limiterLookaheadMs    = 1.0;
    int    oversampleFactor      = 4;    // shared by clipper and limiter
    // 64 rather than the stages' own default of 32, and the reason is measured. The 0.90*Nyquist
    // prototype is a ROUND TRIP, so its loss doubles in dB, and clipper + limiter in series double it
    // again. At 44.1 kHz with both stages at 32 taps the chain costs -1.549 dB at 17.6 kHz, -6.033 at
    // 18.5 and -16.131 at 19.4. At 64 taps: +0.000 / -0.610 / -10.182. A signature like that is not a
    // rounding error on a mastering chain; +64 samples of total latency buys it back.
    int    tapsPerPhase          = 64;

    // The compressor's key filter. > 0 Hz builds a minimum-phase high-passed copy of the compressor's
    // own input and feeds it as the external key (the input P2 added). 0 means self-keyed, which is
    // bit-identical to the three-argument form.
    double sidechainHpfHz = 0.0;
};

// PER-BLOCK parameters. Nothing here moves latencySamples(). Each stage's own parameter type appears
// verbatim — no field is re-declared, so nothing can fall out of step when a stage grows one. That is
// the same rule `CompressorParams : GainReductionParams : DetectorParams` follows; inheritance is not
// available at this level (seven bases, and `releaseMs` exists in both the compressor's and the
// limiter's), so the mechanism is containment by value.
struct MasteringChainParams
{
    double inputGainDb      = 0.0;    // ahead of everything: sets the compressor's operating point
    double preLimiterGainDb = 0.0;    // after the clipper: moves loudness WITHOUT re-compressing

    eq::BandParams                 eqBands[eq::EqEngine::kMaxBands] {};
    stereo::MonoBassParams         monoBass {};
    dynamics::CompressorParams     compressor {};   // lookaheadMs IGNORED — it is topology, see prepare()
    saturation::Saturator::Params  clipper {};
    limiter::TruePeakLimiterParams limiter {};
    dither::DitherParams           dither {};

    // Runtime bypass of a PRESENT stage. Latency-neutral by construction.
    bool bypassEq = false, bypassMonoBass = false, bypassCompressor = false;
    bool bypassClipper = false, bypassLimiter = false, bypassDither = false;
};

// What the chain ACTUALLY applied, after every stage's own clamps and refusals. `configure` in the
// C-ABI has to hand these back — the form's nerd block shows them — and a caller that asked for
// something outside a stage's range can see what it got instead of guessing.
struct MasteringChainResolved
{
    int    latencySamples        = 0;
    int    internalBlock         = 0;
    int    compressorLookahead   = 0;
    int    clipperLatency        = 0;
    int    limiterLatency        = 0;
    int    limiterLookahead      = 0;
    int    oversampleFactor      = 0;
    double limiterCeilingDbTp    = 0.0;
    double limiterReleaseMs      = 0.0;
    stereo::MonoBassParams monoBass {};
};

//==============================================================================
// felitronics::mastering::MasteringChain — the mastering signal chain as ONE streaming, RT-safe,
// block-independent object:
//
//     gate -> inputGain -> EQ -> [M/S mono-bass] -> compressor (opt. keyed) -> [soft clipper]
//          -> preLimiterGain -> true-peak limiter -> dither
//
// Nothing here is new DSP. Every stage is a module that already ships and is already tested; this is
// the composition, and the composition is where the defects live. `OfflineRenderer` is a thin wrapper
// over this — "run to the end, flush, cut the latency" — and not the other way round, so a live
// preview gets the same object a file render uses.
//
// ====================================================================================
// BLOCK INVARIANCE, and why it needs an internal quantum rather than a promise
// ====================================================================================
// The requirement is that the same input and the same parameters give the same output whatever the
// caller's block size. Five of the six stages break that on their own, and NOT in a way any tolerance
// can be written around. Everything below is measured on this tree, not argued:
//
//   * `eq::EqBand` and `stereo::MonoBass` flush their filter state once per process() CALL. State
//     below 1e-15 becomes exact zero at an instant the CALLER chose. Renders at block 4096 and block 1
//     differ in 4721 and 8908 samples respectively.
//   * That is NOT confined to a decaying tail. A band sitting at 0 dB keeps state right at the
//     threshold, so a 50000-sample tone diverges in 2083 of its samples, starting at sample 24.
//   * `dither::Dither` then AMPLIFIES it by eight orders of magnitude. Its auto-blank compares the
//     input to zero EXACTLY, and the EQ flush is what decides whether the tail is exactly zero. One
//     partition exports digital black, the other keeps emitting dither noise: 3.576e-07, three LSB of
//     a 24-bit master, for the whole tail. (With autoBlank off the quantiser absorbs the difference
//     entirely — which is what proves the auto-blank is the amplifier.)
//   * `dynamics::Compressor` is not exempt either, though its own header used to say so. With the RMS
//     window set so the follower coefficient is exactly 0.5 and a key of [1.3e-15, 1.5e-12], the power
//     state after one sample is 8.45e-31 — below the 1e-30 flush floor — while the amplitude after two
//     is 1.06e-12, ABOVE the 1e-12 floor `core::gainToDb` clamps at. One 2-sample call gives 0.478396237,
//     two 1-sample calls give 0.478396297. Half an LSB of 24-bit, out of a floor that was supposed to
//     hide it.
//   * And a whole-file call — the natural shape for an offline render — opens a LAW 8 hole:
//     `core::FlushToZero` reasons that "a block is too short to re-traverse the gap", which is false
//     for a big block. Measured on a 38000-sample tail: 37678 samples subnormal at one call for the
//     whole file, against 0 at block 64. That is the 10-100x stall Law 8 exists to prevent.
//
// So the chain does not hand the caller's block boundaries to anything. It accumulates into a fixed
// `internalBlock` and calls every stage with exactly that, always. Every per-call behaviour — the two
// filter flushes, the compressor's, the limiter's, the EQ's once-per-block coefficient recompute, and
// whatever a future stage does — is then clocked by the STREAM, not by the caller. Measured: 0
// differing samples between callers using 1, 1021 and 4096 samples, at K = 64, 128 and 512, over the
// full chain with 24-bit dither and auto-blank on. It also closes the Law 8 hole for free.
//
// WHAT IT COSTS: exactly `internalBlock` samples of latency, declared in latencySamples() like any
// other. Nothing in CPU — measured 2.26-2.30 %RT for the whole chain at every K from 16 to 4096, i.e.
// flat, so K is chosen for latency and never for speed.
//
// WHAT IT DOES NOT COVER: a parameter change still lands where the caller puts it. Changes are
// deferred to the next quantum boundary, so a change made at stream position p always takes effect at
// the same place regardless of blocking — but a caller that pushes params at different stream
// positions is asking for different renders, which is not this module's to hide. Set them once before
// the render and the question does not arise.
//
// ====================================================================================
// BYPASS, and why "just don't call it" is the wrong default
// ====================================================================================
// A stage is PRESENT (config) or absent; a present stage is ACTIVE or bypassed (params). Latency is
// the sum over PRESENT stages and never moves, which is what makes a bypass toggle safe in a host.
//
// Bypass may not be spelled with neutral-looking numbers. Measured: a limiter with its ceiling at
// +60 dBTP still costs the 0.90*Nyquist round trip (-0.775 dB at 0.40*fs) and its 79 samples; a
// saturator at driveDb = 0 is not linear either, because `WaveShaper` floors the drive at 1e-4.
//
// Which mechanism each stage gets is decided by measurement, not by uniformity:
//   * COMPRESSOR — its own parameters give an exactly transparent, WARM bypass: `ratio = 1` makes the
//     curve's slope exactly 0, so the gain is exactly 1.0f and the signal comes out of the lookahead
//     ring untouched, sign of zero included. Measured bit-exact against the input delayed by its own
//     latency. The detector keeps tracking, so un-bypassing does not jump. Nothing to align.
//   * CLIPPER and LIMITER — skipped, with a `core::DryAligner` holding their PDC. The clipper does have
//     a `mix = 0` bypass that is bit-exact for ordinary audio, but it normalises -0.0f to +0.0f
//     (`dry + 0.0f * wet`), which a bit-exact null test would report, and it pays the whole oversampled
//     wet path for nothing. The limiter has no bypass at all.
//   * EQ, MONO-BASS, DITHER — zero latency, so bypass is simply not calling them. They are reset on the
//     transition, which is exactly what `EqBand` and `MonoBass` do to themselves when they go idle.
// A skipped stage is reset when the bypass flag changes, in either direction. Without that a stage
// that sat out re-emits audio from before the gap — the defect class measured at +19.76 dB over the
// ceiling in the limiter (#119), 0.75 out of digital silence in the compressor (#120), 0.93 in the
// saturator and +7.39 dBFS in the EQ (both still open, see .private/p6-findings.md).
//
// ====================================================================================
// THE CHANNEL COUNT IS EXACT
// ====================================================================================
// prepare() takes the channel count, and process() REFUSES any other, touching neither the buffer nor
// any state. It is not clamped to a prefix and the count is not allowed to vary. Two reasons, both
// concrete: the stages disagree — the compressor refuses a wider call, the limiter, saturator and EQ
// process a PREFIX and leave the rest untouched, mono-bass ignores anything that is not exactly two —
// so reconciling them from out here is precisely the field-mapping that falls out of step; and the
// C-ABI this feeds (`fc_master_create(sampleRate, channels)`) fixes the count at creation anyway.
// This turns the "stereo -> mono -> stereo" sequence into a provable property: a refused call is
// indistinguishable from one never made, so [A, refused, B] is bit-identical to [A, B].
//
// ====================================================================================
// THE GATE
// ====================================================================================
// One `isfinite`/clamp on every input sample, ahead of everything, in the shape `saturation::Saturator`
// and `limiter::TruePeakLimiter` already use. It is not decoration: `eq::Biquad` clears poison only at
// the end of a call, so one Inf reaching the EQ costs the rest of the quantum, and with every stage
// bypassed an Inf would otherwise reach the output untouched. With the gate, the chain's behaviour on
// a bad sample is ONE rule that does not depend on which stages are on: feeding a NaN is bit-identical
// to feeding the sanitised value. Bit-transparent for any finite sample within +-1e6.
//
// RT-safe: prepare() allocates, process()/flush() do not allocate, lock or throw, and accept any block
// length. `sizeof(eq::EqEngine)` is ~324 KB, so the engine is held behind a pointer and this object
// stays small enough to put on a stack.
class MasteringChain
{
public:
    static constexpr int    kMaxInternalBlock = 8192;
    static constexpr int    kMinInternalBlock = 8;
    static constexpr double kMaxSampleRate    = 3.0e6;
    static constexpr double kMaxGainDb        = 60.0;    // both gain nodes; beyond this is not a trim

    // Returns false and leaves the chain UNPREPARED on anything it cannot honour. A false return is
    // the only way to learn that, so check it. Spelled positively so a NaN fails: two of the stages
    // below accept a NaN sample rate through `sampleRate <= 0.0` and go on to emit NaN, so the chain
    // validates once, here, rather than trusting them.
    [[nodiscard]] bool prepare (double sampleRate, int numChannels, const MasteringChainConfig& config = {})
    {
        prepared_ = false;                                     // any early return leaves it unprepared
        if (! (sampleRate > 0.0) || ! std::isfinite (sampleRate) || sampleRate > kMaxSampleRate) return false;
        if (numChannels < 1 || numChannels > core::kMaxChannels) return false;
        if (config.internalBlock < kMinInternalBlock || config.internalBlock > kMaxInternalBlock) return false;
        if (config.oversampleFactor < 2 || config.tapsPerPhase < 4) return false;
        if (! std::isfinite (config.compressorLookaheadMs) || config.compressorLookaheadMs < 0.0) return false;
        if (! std::isfinite (config.limiterLookaheadMs) || config.limiterLookaheadMs < 0.0) return false;
        if (! std::isfinite (config.sidechainHpfHz) || config.sidechainHpfHz < 0.0
            || config.sidechainHpfHz >= 0.5 * sampleRate) return false;
        // REFUSED, not silently ignored. `stereo::MonoBass` leaves a non-stereo buffer untouched, so a
        // mono chain that reported this stage as enabled would be reporting a stage that does nothing —
        // the class of silent no-op this plan keeps closing.
        if (config.monoBass && numChannels != 2) return false;

        cfg_ = config;
        fs_  = sampleRate;
        nch_ = numChannels;
        K_   = config.internalBlock;

        // One buffer, not two: process() SWAPS the caller's samples with the quantum buffer, so the
        // slot a sample is read from is the slot the next input goes into. Latency is exactly K.
        fifo_.assign ((std::size_t) K_ * (std::size_t) nch_, 0.0f);
        pos_ = 0;

        if (cfg_.eq)
        {
            eq_ = std::make_unique<eq::EqEngine>();
            eq_->prepare (fs_, K_, nch_);
        }
        else eq_.reset();

        if (cfg_.monoBass) monoBass_.prepare (fs_, K_, nch_);

        int compLat = 0;
        if (cfg_.compressor)
        {
            // maxLookaheadMs is what sizes the ring, so it must be at least what the config asks for —
            // otherwise the compressor CLAMPS the lookahead and reports a latency smaller than the one
            // this chain would have computed. That mismatch is exactly why latency is read back below
            // rather than computed here.
            if (! comp_.prepare (fs_, K_, nch_, std::max (cfg_.compressorLookaheadMs, 1.0))) return false;
            dynamics::CompressorParams cp;
            cp.lookaheadMs = cfg_.compressorLookaheadMs;
            comp_.setParams (cp);
            compLat = comp_.latencySamples();
        }

        int clipLat = 0;
        if (cfg_.clipper)
        {
            if (! sat_.prepare (fs_, K_, nch_, cfg_.oversampleFactor, cfg_.tapsPerPhase)) return false;
            clipLat = sat_.latencySamples();
            alignClip_.prepare (nch_, K_, clipLat + 2);        // capacity must EXCEED the delay it will hold
        }

        int limLat = 0;
        if (cfg_.limiter)
        {
            limiter::TruePeakLimiterConfig lc;
            lc.lookaheadMs      = cfg_.limiterLookaheadMs;
            lc.oversampleFactor = cfg_.oversampleFactor;
            lc.tapsPerPhase     = cfg_.tapsPerPhase;
            if (! lim_.prepare (fs_, K_, nch_, lc)) return false;
            limLat = lim_.latencySamples();
            alignLim_.prepare (nch_, K_, limLat + 2);
        }

        if (cfg_.dither) dith_.prepare (fs_, K_, nch_);

        if (cfg_.sidechainHpfHz > 0.0)
        {
            keyBuf_.assign ((std::size_t) K_ * (std::size_t) nch_, 0.0f);
            const eq::BiquadCoeffs hc = eq::matched::highpass (cfg_.sidechainHpfHz, fs_, 0.70710678118654752);
            for (int c = 0; c < nch_; ++c) hpf_[c].setCoeffs (hc);
        }
        else keyBuf_.clear();

        // READ BACK, never computed. A chain that derives the number itself can agree with its own
        // aligners while disagreeing with the stage: `Compressor::prepare(.., maxLookaheadMs = 50)`
        // caps the lookahead at 2400 samples, so a chain asking for 60 ms and computing lround(2880)
        // would hold BOTH its latency and its aligner at 2880 and pass a bypass null test while the
        // active chain sat 480 samples out of alignment.
        latency_ = K_ + compLat + clipLat + limLat;

        pendingParams_ = params_;
        paramsDirty_   = true;
        reset();
        prepared_ = true;
        return true;
    }

    void reset() noexcept
    {
        std::fill (fifo_.begin(), fifo_.end(), 0.0f);
        pos_ = 0;
        if (eq_) eq_->reset();
        if (cfg_.monoBass)   monoBass_.reset();
        if (cfg_.compressor) comp_.reset();
        if (cfg_.clipper)  { sat_.reset();  alignClip_.reset(); }
        if (cfg_.limiter)  { lim_.reset();  alignLim_.reset(); }
        if (cfg_.dither)     dith_.reset();
        for (int c = 0; c < core::kMaxChannels; ++c) hpf_[c].reset();
        std::fill (keyBuf_.begin(), keyBuf_.end(), 0.0f);
        // A stream restart has to restore the PARAMETER smoothers too, and summing the stages' own
        // resets does not do it: `EqBand::reset()` clears filter state and deliberately leaves the
        // freq/Q/gain smoothers where they are. Measured before this line existed — set a band to 0 dB,
        // render, set it to +9 dB (a 30 ms ramp), render 10 ms of it, reset, render again: the result
        // differed from a freshly prepared chain by up to 0.51, full scale, because the ramp simply
        // resumed. That would have made `OfflineRenderer`'s "two renders of the same input are
        // bit-identical" false for any programme whose parameters moved during the first pass.
        paramsDirty_ = true;
        forceSnap_   = true;
        bypassKnown_ = false;
    }

    // Takes effect at the next internal quantum boundary, which is what keeps a parameter change from
    // depending on where the caller happened to cut the stream. `compressor.lookaheadMs` is ignored —
    // it is topology (it moves latency), and it is overwritten from the config.
    void setParams (const MasteringChainParams& p) noexcept
    {
        pendingParams_ = p;
        paramsDirty_   = true;
    }

    const MasteringChainParams& params() const noexcept { return params_; }

    int  latencySamples() const noexcept { return prepared_ ? latency_ : 0; }
    int  numChannels()    const noexcept { return prepared_ ? nch_ : 0; }
    int  internalBlock()  const noexcept { return prepared_ ? K_ : 0; }
    bool isPrepared()     const noexcept { return prepared_; }

    MasteringChainResolved resolved() const noexcept
    {
        MasteringChainResolved r;
        if (! prepared_) return r;
        r.latencySamples      = latency_;
        r.internalBlock       = K_;
        r.compressorLookahead = cfg_.compressor ? comp_.latencySamples() : 0;
        r.clipperLatency      = cfg_.clipper ? sat_.latencySamples() : 0;
        r.limiterLatency      = cfg_.limiter ? lim_.latencySamples() : 0;
        r.limiterLookahead    = cfg_.limiter ? lim_.lookaheadSamples() : 0;
        r.oversampleFactor    = cfg_.limiter ? lim_.oversampleFactor()
                                             : (cfg_.clipper ? cfg_.oversampleFactor : 0);
        r.limiterCeilingDbTp  = cfg_.limiter ? lim_.effectiveCeilingDbTp() : 0.0;
        r.limiterReleaseMs    = cfg_.limiter ? lim_.effectiveReleaseMs() : 0.0;
        r.monoBass            = cfg_.monoBass ? monoBass_.params() : stereo::MonoBassParams { false, 0.0f, 0.0f };
        return r;
    }

    // Audio thread, in place, planar. RT-safe. `numSamples` may be anything at all, including a whole
    // file. Returns false — TOUCHING NOTHING — if the chain is unprepared or the channel count is not
    // the prepared one; a refused call is indistinguishable from one never made.
    bool process (float* const* io, int numChannels, int numSamples) noexcept
    {
        if (! prepared_ || io == nullptr) return false;
        if (numChannels != nch_) return false;
        if (numSamples < 0) return false;
        if (numSamples == 0) return true;

        for (int off = 0; off < numSamples; )
        {
            const int take = std::min (K_ - pos_, numSamples - off);
            for (int c = 0; c < nch_; ++c)
            {
                float* slot = fifo_.data() + (std::size_t) c * (std::size_t) K_ + (std::size_t) pos_;
                std::swap_ranges (slot, slot + take, io[c] + off);
            }
            pos_ += take;
            off  += take;
            if (pos_ == K_) { runQuantum(); pos_ = 0; }
        }
        return true;
    }

    // Drain the chain: writes exactly min(latencySamples(), capacity) frames and returns that count.
    // It is not a second code path — it IS process() over that many zeros, which is what makes "the
    // tail is not lost" a definition rather than a promise, and what a null test can be written
    // against. `render(x)` equals `process(x followed by latencySamples() zeros)` with the first
    // `latencySamples()` output samples dropped.
    int flush (float* const* out, int numChannels, int capacity) noexcept
    {
        if (! prepared_ || out == nullptr || numChannels != nch_ || capacity <= 0) return 0;
        const int n = std::min (latency_, capacity);
        if (n <= 0) return 0;
        for (int c = 0; c < nch_; ++c) std::fill (out[c], out[c] + n, 0.0f);
        (void) process (out, numChannels, n);
        return n;
    }

private:
    // One quantum: exactly K_ samples, every stage, always. This is the only place a stage is called.
    void runQuantum() noexcept
    {
        float* ch[core::kMaxChannels] {};
        for (int c = 0; c < nch_; ++c) ch[c] = fifo_.data() + (std::size_t) c * (std::size_t) K_;

        if (paramsDirty_) { applyParams(); paramsDirty_ = false; }

        // --- the gate, ahead of everything (see the header note) -------------------------------
        for (int c = 0; c < nch_; ++c)
            for (int i = 0; i < K_; ++i)
            {
                const float v = ch[c][i];
                ch[c][i] = std::clamp (std::isfinite (v) ? v : 0.0f, -1.0e6f, 1.0e6f);
            }

        applyGain (ch, inputGain_);

        // --- EQ (zero latency: bypass is simply not calling it) --------------------------------
        if (eq_ != nullptr)
        {
            if (! params_.bypassEq) eq_->process (ch, nch_, K_);
            else if (bypassChanged_.eq) eq_->reset();
        }

        // --- M/S mono-bass (zero latency) ------------------------------------------------------
        if (cfg_.monoBass)
        {
            if (! params_.bypassMonoBass) monoBass_.process (ch, nch_, K_);
            else if (bypassChanged_.monoBass) monoBass_.reset();
        }

        // --- compressor: WARM bypass through its own curve, so nothing has to be aligned --------
        if (cfg_.compressor)
        {
            if (! keyBuf_.empty())
            {
                // The key is the compressor's OWN input, same instant, minimum-phase high-passed. It is
                // deliberately not pre-shifted: the compressor supplies the lookahead by delaying the
                // programme, so key sample i sets the gain applied to programme sample i - latency.
                // The filter's group delay (sqrt(2)/(2*pi*fc) at the corner) eats into the lookahead and
                // is NOT declared latency — it is frequency-dependent and cannot be one number.
                const float* key[core::kMaxChannels] {};
                for (int c = 0; c < nch_; ++c)
                {
                    float* k = keyBuf_.data() + (std::size_t) c * (std::size_t) K_;
                    for (int i = 0; i < K_; ++i) k[i] = hpf_[c].processSample (ch[c][i]);
                    hpf_[c].flushDenormals();
                    key[c] = k;
                }
                comp_.process (ch, nch_, K_, key, nch_);
            }
            else comp_.process (ch, nch_, K_);
        }

        // --- soft clipper: skipped when bypassed, its PDC held by the aligner -------------------
        if (cfg_.clipper)
        {
            // The aligner is advanced whether or not the stage runs — a ring fed only while bypassed is
            // cold at the moment it is first read and emits its latency in zeros (DryAligner.h says so
            // in as many words). ONE reset on any change of the flag, before the branch: skipping a
            // stage freezes its history, and a frozen oversampler replays pre-gap audio on re-entry.
            alignClip_.advance ((const float* const*) ch, nch_, K_, sat_.latencySamples());
            if (bypassChanged_.clipper) sat_.reset();
            if (! params_.bypassClipper) sat_.process (ch, nch_, K_);
            else for (int c = 0; c < nch_; ++c) std::copy_n (alignClip_.delayed (c), K_, ch[c]);
        }

        applyGain (ch, preLimGain_);

        // --- true-peak limiter: the one stage with no bypass of its own -------------------------
        if (cfg_.limiter)
        {
            alignLim_.advance ((const float* const*) ch, nch_, K_, lim_.latencySamples());
            if (bypassChanged_.limiter) lim_.reset();
            if (! params_.bypassLimiter) lim_.process (ch, nch_, K_);
            else for (int c = 0; c < nch_; ++c) std::copy_n (alignLim_.delayed (c), K_, ch[c]);
        }

        // --- dither, last, and only when it is not bypassed --------------------------------------
        if (cfg_.dither && ! params_.bypassDither) dith_.process (ch, nch_, K_);

        bypassChanged_ = {};
    }

    void applyGain (float* const* ch, float g) noexcept
    {
        if (core::exactlyEqual (g, 1.0f)) return;              // 0 dB is a bit-exact no-op, not a multiply
        for (int c = 0; c < nch_; ++c)
            for (int i = 0; i < K_; ++i) ch[c][i] *= g;
    }

    struct BypassFlags { bool eq = false, monoBass = false, compressor = false, clipper = false, limiter = false, dither = false; };

    void applyParams() noexcept
    {
        const MasteringChainParams& p = pendingParams_;

        if (! bypassKnown_)
        {
            bypassChanged_ = { true, true, true, true, true, true };   // first quantum: settle everything
            bypassKnown_   = true;
        }
        else
        {
            bypassChanged_.eq         = bypassChanged_.eq         || (p.bypassEq         != params_.bypassEq);
            bypassChanged_.monoBass   = bypassChanged_.monoBass   || (p.bypassMonoBass   != params_.bypassMonoBass);
            bypassChanged_.compressor = bypassChanged_.compressor || (p.bypassCompressor != params_.bypassCompressor);
            bypassChanged_.clipper    = bypassChanged_.clipper    || (p.bypassClipper    != params_.bypassClipper);
            bypassChanged_.limiter    = bypassChanged_.limiter    || (p.bypassLimiter    != params_.bypassLimiter);
            bypassChanged_.dither     = bypassChanged_.dither     || (p.bypassDither     != params_.bypassDither);
        }
        params_ = p;

        inputGain_  = gainOf (p.inputGainDb);
        preLimGain_ = gainOf (p.preLimiterGainDb);

        if (eq_)
        {
            // Make the write below behave like the FIRST write after prepare(), which snaps. `EqBand`
            // snaps a lane whose previous state was OFF and ramps one that stays on, so writing every
            // lane off first turns the real write into a snap. There is no other way in: `initialized`
            // is only cleared by `EqBand::prepare`, and reaching it means `EqEngine::prepare`, which
            // allocates. Two stores per band, on a reset, and nothing at all on the streaming path.
            if (forceSnap_)
            {
                for (int i = 0; i < eq::EqEngine::kMaxBands; ++i)
                {
                    eq::BandParams off = p.eqBands[i];
                    off.on = false;
                    for (eq::LaneParams& l : off.lanes) l.on = false;
                    eq_->setBand (i, off);
                }
                forceSnap_ = false;
            }
            for (int i = 0; i < eq::EqEngine::kMaxBands; ++i) eq_->setBand (i, p.eqBands[i]);
        }
        if (cfg_.monoBass) monoBass_.setParams (p.monoBass);

        if (cfg_.compressor)
        {
            dynamics::CompressorParams cp = p.compressor;
            cp.lookaheadMs = cfg_.compressorLookaheadMs;       // topology, never a per-block field
            if (p.bypassCompressor)
            {
                // The warm bypass, and it is exact rather than approximately transparent: `ratio = 1`
                // gives the curve a slope of exactly 0, so the delta is exactly 0 dB and the gain is
                // exactly 1.0f. The signal leaves the lookahead ring untouched, sign of zero included,
                // while the detector keeps tracking — so coming back out of bypass does not jump.
                cp.mode       = dynamics::Mode::DownCompress;
                cp.ratio      = 1.0;
                cp.makeupDb   = 0.0;
                cp.autoMakeup = false;
            }
            comp_.setParams (cp);
        }

        if (cfg_.clipper) sat_.setParams (p.clipper);
        if (cfg_.limiter) lim_.setParams (p.limiter);
        if (cfg_.dither)  dith_.setParams (p.dither);
    }

    static float gainOf (double db) noexcept
    {
        const double d = std::isfinite (db) ? std::clamp (db, -kMaxGainDb, kMaxGainDb) : 0.0;
        return (float) core::dbToGain (d);
    }

    double fs_ = 48000.0;
    int    nch_ = 0, K_ = 0, pos_ = 0, latency_ = 0;
    bool   prepared_ = false, paramsDirty_ = true, bypassKnown_ = false, forceSnap_ = true;

    MasteringChainConfig cfg_ {};
    MasteringChainParams params_ {}, pendingParams_ {};
    BypassFlags          bypassChanged_ {};

    float inputGain_ = 1.0f, preLimGain_ = 1.0f;

    std::vector<float> fifo_, keyBuf_;

    std::unique_ptr<eq::EqEngine> eq_;                 // ~324 KB — behind a pointer so this object is stack-sized
    stereo::MonoBass              monoBass_;
    dynamics::Compressor          comp_;
    saturation::Saturator         sat_;
    limiter::TruePeakLimiter      lim_;
    dither::Dither                dith_;
    core::DryAligner              alignClip_, alignLim_;
    eq::Biquad                    hpf_[core::kMaxChannels] {};
};

} // namespace felitronics::mastering

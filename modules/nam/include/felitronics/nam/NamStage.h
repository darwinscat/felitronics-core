// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <cstddef>
#include <memory>

//==============================================================================
// felitronics::nam::NamStage — a neural amp (NAM) stage that runs IN FRONT of cab convolution
// (signal path: input → AMP → cab). Same seam discipline as CabConvolver: the public API speaks
// raw float* + sizes, and the Neural Amp Modeler inference engine (Eigen + nlohmann/json) is hidden
// in the .cpp behind a pImpl — so NAM never leaks into the rest of the Core, and the dependency
// stays swappable.
//
// Threading mirrors the IR path exactly:
//   • prepare() / loadModelFromMemory() / clearModel() / collectGarbage() — message thread.
//   • process() — the ONLY thing the audio thread calls. 🔴 It never allocates, locks,
//     does IO or throws. A freshly-loaded model is atomic-swapped into the live pointer;
//     the replaced model is parked and freed on the message thread (collectGarbage) only
//     once the audio thread has provably stepped past it — so no use-after-free and the
//     audio thread never deletes.
//
// The no-allocation process() guarantee applies to architectures whose pinned NAM implementation
// preallocates its work in Reset (Linear and WaveNet). LSTM at this pin returns an owning dynamic
// Eigen hidden-state vector per sample; ConvNet constructs dynamic Eigen temporaries per block.
// Those inherited OrbitCab behaviors are tracked upstream and are not rejected at load, preserving
// byte-compatibility with OrbitCab.
//
// The swap machinery itself (atomic live pointer, block-counter retire, message-thread GC) is
// the shared felitronics::neural::NeuralStage — extracted from OrbitCab's original AmpStage;
// what lives here is the NAM-specific backend (dual-instance true stereo, loudness makeup, .namz
// unpack, rate-match). This class is the stable public seam, so product adapters can simply alias it.
//
// Channels: a MONO stream (1 ch) runs ONE model instance on the single channel; a STEREO stream
// (2 ch) runs TWO independent instances of the SAME capture (true stereo — L/R independent). Both
// instances are always built on load, and prepare() configures both per-channel resamplers, so a
// host switching the layout mono<->stereo (a re-prepare) is safe.
//==============================================================================
namespace felitronics::nam
{

class NamStage
{
public:
    NamStage();
    ~NamStage();

    // Allocate the mono scratch for this stream and (re)configure a live model for the
    // new sample-rate / block size. Message/host thread (prepareToPlay) — never the audio
    // thread (it can allocate + prewarm the network).
    void prepare (double sampleRate, int maxBlock);
    void reset();

    // 🔴 RT-safe, in place. No model loaded → clean passthrough (no-op, and an ACCEPTED call).
    // `normalize` applies the model's loudness makeup (output normalisation) when the model
    // carries a loudness tag — brings raw model output to a consistent reference level.
    // Law 11: `numSamples` is any length (chunked internally); `numChannels` is 1 or 2 and anything
    // else is REFUSED whole — false means nothing was touched. See DSP-ARCHITECTURE.md §2 law 11.
    [[nodiscard]] bool process (float* const* io, int numChannels, int numSamples, bool normalize) noexcept;

    //--- model lifecycle (message thread) ----------------------------------------
    // Build a NAM model from raw .nam bytes off the audio thread and atomic-swap it in.
    // Returns false (and leaves the current model untouched) on a bad / unsupported /
    // non-mono model. The replaced model is reclaimed later via collectGarbage().
    // A load (or clear) that was ACCEPTED is guaranteed to apply: immediately in the normal
    // case, or — if the bounded swap-retire queue is momentarily full (audio frozen across
    // many swaps) — on the next collectGarbage()/prepare() drain. The info getters below
    // update only when it actually lands, so they always describe the model that is live.
    // trimDb = per-model output offset (dB) folded into the loudness-normalisation makeup.
    // (Bytes only — file I/O stays in the adapter layer, never in pure-DSP core.)
    bool   loadModelFromMemory (const void* data, std::size_t size, float trimDb = 0.0f);

    //--- the same load, in two halves ----------------------------------------------
    // The heavy half — the bytes parsed, both instances built, the rate-match and the prewarm
    // prepared — done ANYWHERE BUT the audio thread and touching no stage at all, so a host can run
    // it on a worker while its drawing thread stays free: a WaveNet costs some twenty milliseconds
    // here, which a hand on a knob feels as a hiccup. `sampleRate` and `maxBlock` are what the stage
    // was prepared with; a model prepared for other numbers is prepared again by install(), on that
    // thread, at the old cost. Null = a bad or unsupported model, the cases loadModelFromMemory()
    // refuses — except the rate contract, which only a stage can judge (install() refuses those).
    // loadModelFromMemory() IS prepareModel() followed by install(): one path, two entry points.
    struct Prepared;
    struct PreparedDeleter { void operator() (Prepared*) const noexcept; };
    using PreparedModel = std::unique_ptr<Prepared, PreparedDeleter>;
    static PreparedModel prepareModel (const void* data, std::size_t size, double sampleRate, int maxBlock,
                                       float trimDb = 0.0f);
    // The light half: message thread, a pointer swap. False — and the stage untouched — for a null
    // handle or a model whose native rate breaks the rate contract (see loadModelFromMemory). Accepted
    // = guaranteed to apply, exactly as a load: immediately, or on the next drain if the retire queue
    // is momentarily full.
    bool   install (PreparedModel model);
    void   clearModel();
    bool   collectGarbage();          // free models retired by a swap, once audio moved past them,
                                      // and land any load/clear deferred by a full retire queue.
                                      // Returns true when a DEFERRED intent landed on this call —
                                      // model state (and possibly latencySamples()) changed, so the
                                      // caller re-reports host PDC like after a normal load.

    bool   hasModel() const;
    double modelSampleRate() const;   // the model's expected sample rate (<= 0 if none / unknown)
    double modelLoudness()   const;   // the model's tagged loudness in dB (0 if none / no model)
    bool   modelHasLoudness() const;  // whether the loaded model carries a loudness tag
    int    latencySamples()  const;   // host-rate latency from rate-matching (0 if none / not resampling)

    //--- THE RATE-MATCH, ANSWERED WITHOUT A MODEL ---------------------------------------------
    // 🔴 WHY THIS IS PUBLIC AND STATIC. Three facts decide what a rate-match costs — the run rate a
    // model will actually get, whether a resampler is installed at all, and the host-rate delay if one
    // is — and until now all three lived inside an INSTANCE method that needs a loaded, prepared
    // model. A consumer that must size a fixed delay line BEFORE any model exists therefore rewrote
    // the arithmetic from a comment, and rewrote it wrong: a downstream bypass path capped itself at
    // 64 samples on a formula two kernel generations stale, and silently under-delayed every host rate
    // above 48 kHz. Restating this is the defect; asking is the fix.
    //
    // The ORDER of the three is part of the contract and not an implementation detail: the model rate
    // is normalised FIRST, and the gate then sees the normalised value. An untagged model (rate <= 0)
    // therefore runs at the default and is NOT resampled at a default-rate host — reverse the two and
    // that case changes.
    // The rate a model runs at when it reports none, and in practice the rate almost every model
    // runs at: install() refuses a tagged model whose rate differs from the stage's current run rate
    // by more than half a hertz, and that rate is only ever re-derived from a model that already
    // passed the same check — so a fresh stage accepts 48 kHz captures and refuses 44.1 and 96 kHz
    // ones outright (measured: every public route was tried).
    //
    // ⚠️ IT IS NOT A PROVABLE CEILING, and an earlier version of this comment claimed it was. The
    // check is a TOLERANCE, not equality: a model at 48000.5 is accepted, prepare() then adopts that
    // rate, and the next half-hertz step is accepted against the new one. Measured: 65 such steps
    // walked the run rate from 48000 to 47967.5. A LOWER run rate means a LONGER round trip, so
    // sizing a buffer from this constant alone can come up short — floor it at what already shipped.
    //
    // It is public because a consumer sizing a delay line before any model exists has to start
    // somewhere, and until now it could not even name this number: a downstream repository invented a
    // "lowest pack rate" of 8 kHz to stand in for it, and sized itself wrong.
    static constexpr double kModelSampleRate = 48000.0;

    struct RateMatch
    {
        double modelRunSR     = kModelSampleRate;   // the rate the model would be run at, normalised
        bool   resampling     = false;              // whether a rate-matcher would be installed
        int    latencySamples = 0;                  // host-rate latency it costs; 0 when none is
    };

    // Pure function of the two rates: no state, no model, no allocation. `modelSR` is the rate a model
    // REPORTS (<= 0 meaning "unknown"), not one already normalised.
    //
    // It answers for the rates it is GIVEN. It does not know whether a stage would accept a model at
    // that rate — install() has its own contract — so rateMatch(h, 44100) describes what 44.1 kHz
    // would cost, not a configuration a fresh stage can reach.
    //
    // Precondition, documented rather than enforced because this extraction promises that no number
    // moves: hostSR positive and finite. Outside that the answer is whatever the arithmetic gives
    // (h = 0 reports 32 against a real 64; h = inf reports -1), exactly as it did before.
    static RateMatch rateMatch (double hostSR, double modelSR) noexcept;
    // How many samples this model must be FED before its output means anything — its receptive field.
    // A network with empty buffers describes the silence it was born into for exactly this long, so
    // anything that fades a freshly loaded model in has to run it silently for this many samples
    // first. NAM answers for convnet/lstm and a container forwards to its submodel, but WaveNet does
    // not override it — so for a WaveNet capture this is computed from the config instead of trusted
    // as zero. 0 = no model, or an architecture nothing here can read.
    int    prewarmSamples()  const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl;
};

} // namespace felitronics::nam

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
    // thread, at the old cost. Null = a bad or unsupported model — EVERY case loadModelFromMemory()
    // refuses, the rate contract included: acceptsModelRate() reads no stage state, so this half can
    // and does judge it. What that saves is the PREWARM and the rate configuration, not the parse: the
    // tag is read off the built network (`GetExpectedSampleRate`), so both instances exist by the time
    // the contract can be asked. An earlier wording here said "the network is never built", which is
    // simply false and was caught by a review round rather than by a test.
    // loadModelFromMemory() IS prepareModel() followed by install(): one path, two entry points.
    struct Prepared;
    struct PreparedDeleter { void operator() (Prepared*) const noexcept; };
    using PreparedModel = std::unique_ptr<Prepared, PreparedDeleter>;
    static PreparedModel prepareModel (const void* data, std::size_t size, double sampleRate, int maxBlock,
                                       float trimDb = 0.0f);
    // The light half: message thread, a pointer swap. False — and the stage untouched — for a null
    // handle, and for the one case left that a handle cannot answer on its own: a re-preparation for
    // THIS stage's host rate and block that the backend cannot honour (see prepare()'s refusal). The
    // rate CONTRACT is settled in prepareModel(), so "wrong model rate" is no longer among the reasons;
    // "a non-null handle is always installable" would be a stronger claim than the code makes and an
    // earlier draft of this line made it. Accepted = guaranteed to apply,
    // exactly as a load: immediately, or on the next drain if the retire queue is momentarily full.
    // A handle prepared for a different host rate or block size is re-prepared here, and the test for
    // that is EXACT: a tolerance there let a backend keep a rate-match computed for a host it is not
    // installed into, and report 0 samples of latency where the policy charges 64 (and 64 where the
    // policy charges 0 — both measured). Re-preparing a stateful LSTM also prewarms its existing cell
    // again, even if both hosts use the direct path. A mismatched split load can therefore start from
    // a different state than a fused load; Reset is not an idempotent state restoration in pinned NAM.
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
    // The rate a model runs at when it reports none, and — within kModelRateTolerance — the rate every
    // model this stage will accept runs at. It is a PROVABLE CEILING now, and the two sentences that
    // used to stand here saying otherwise are the subject of this fix rather than a caveat to it.
    //
    // 🔴 IT IS A FIXED WINDOW, AND IT USED TO BE A MOVING ONE — that was the whole defect. The gate
    // compared a model's tag against the rate of whatever was live at the last prepare(), and prepare()
    // ADOPTED the accepted tag, so a fuzz on equality was a random walk: each accepted load moved the
    // reference the next load is judged against. Measured on the base commit, half-hertz steps down,
    // host 48000, loop capped at 5000:
    //
    //     load only .................... 1 step      ← nothing else moves the reference
    //     load + process() ............. 1 step      ← audio is NOT the clock
    //     load + prepare() ............. 66 steps    ← stopped by kMaxRetiredModels, not by rates
    //     load + process() + prepare() . 5000 steps, no refusal → run rate 45500.0
    //
    // So the walk was clocked by prepare(), and audio only drained the retire queue — which corrects
    // BOTH numbers this comment used to carry ("65 steps", and "run audio and the walk does not stop").
    // The cost was not an exotic one either: a stage walked to 47900 REFUSES an ordinary 48000 capture,
    // measured. The window moved; it never widened.
    //
    // The reference is now the constant itself, so no sequence of loads can move it, and a consumer may
    // size a buffer from this number — which is what it is public for. Ask maxLatencySamples() rather
    // than deriving one: a downstream repository invented a "lowest pack rate" of 8 kHz to stand in for
    // this constant and sized itself wrong, and deriving from the NOMINAL rate is a second way to get
    // the same answer wrong, because the accepted window's LOW edge is what costs the most.
    static constexpr double kModelSampleRate = 48000.0;

    // The half-hertz of tag noise the gate forgives. It is the SAME number rateMatch() uses to decide
    // whether a resampler is worth installing, and deliberately so: "this stage accepts a model exactly
    // when a factory-rate host would not resample it" is one rule with one owner, and acceptsModelRate()
    // below is that sentence spelled as code rather than a second copy of the 0.5.
    //
    // 🔴 AND rateMatch() SPENDS THIS CONSTANT RATHER THAN ITS OWN LITERAL, which it did not at first.
    // The mutation stand is what said so: moving this number to 0.6 changed nothing anywhere, because
    // the only thing reading it was maxLatencySamples() while the gate still carried a private 0.5.
    // A constant that names a rule nobody consults is a restatement waiting to drift, and this one
    // would have drifted silently in the direction that widens the accepted window.
    static constexpr double kModelRateTolerance = 0.5;

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
    // moves — and stated as the range it actually KEEPS, since "positive and finite" was measured to be
    // wider than the arithmetic supports: hostSR in (0, 3.22e12]. Outside that the answer is whatever
    // the arithmetic gives, exactly as it did before the extraction, and there are four regimes rather
    // than the two the .cpp used to name — they are enumerated with their thresholds at rateMatch()'s
    // definition. The two that bite: past ~3.22e12 the narrowing of lround's long to int invents a
    // plausible positive answer with no flag raised, and past ~1.38e22 (an infinity included) the
    // answer is whatever that platform's lround saturates to, which is NOT the same on all of them —
    // measured on four rows, see the .cpp. h = 0 reports 32 against a real 64.
    static RateMatch rateMatch (double hostSR, double modelSR) noexcept;

    // 🔴 THE RATE CONTRACT, AS A PURE PREDICATE ON THE MODEL'S OWN TAG. True = this stage will take a
    // model reporting `modelSR`; false = prepareModel() returns null for it and no stage anywhere will
    // accept it. It depends on NOTHING but the argument — that is the fix, and it is why the check now
    // lives in prepareModel() rather than install(): a rule that reads no stage state is not a stage's
    // to judge, and judging it in the heavy half saves a WaveNet's twenty milliseconds of PREWARM on a
    // model that was never going to load. Not the parse and not the network: the tag is read off the
    // built instances, so those exist by the time this can be asked — a claim that said otherwise stood
    // here until a review round measured it.
    //
    // A model that reports NO rate (<= 0) is accepted and runs at kModelSampleRate — rateMatch() owns
    // that normalisation, and this asks it rather than repeating it. The window is therefore
    // [kModelSampleRate - kModelRateTolerance, kModelSampleRate + kModelRateTolerance], closed at both
    // ends, for the life of the process.
    //
    // ⚠️ THE ARGUMENT IS A double AND THE EDGES NEED IT. Measured: 48000.5001 is refused, but
    // `(double)(float) 48000.5001` is exactly 48000.5 and is ACCEPTED — a tag that passes through a
    // float anywhere upstream collapses onto the edge and widens the effective window by about two
    // thousandths of a hertz. NAM reports the tag as a double and nothing in this repository narrows
    // it; a consumer that does is choosing a slightly different window.
    static bool acceptsModelRate (double modelSR) noexcept;

    // An upper bound on latencySamples() over every model this stage would ACCEPT at this host rate —
    // the number a consumer sizing a fixed delay line actually needs, and the reason it must not derive
    // one itself. A LOWER model rate is a LONGER round trip, so the bound is taken at the accepted
    // window's low edge, not at the nominal rate: at a 3 MHz host that difference is 0.0208 samples,
    // which is nothing until it lands on the wrong side of a rounding boundary and a ring comes up one
    // slot short.
    //
    // A BOUND rather than an attained maximum, and it has TWO gaps rather than the one an earlier
    // wording admitted: at hostSR exactly kModelSampleRate no accepted model resamples at all, so the
    // true maximum there is 0 and this still answers 64; and at any host too slow for a backend to be
    // prepared at all (see prepare()'s refusal) nothing runs, so nothing attains it. It is a policy
    // answer about rates, not a promise about a particular prepared instance. See the definition for
    // why the bound is taken on the geometry instead of on rateMatch().
    //
    // Same precondition as rateMatch(): hostSR in (0, 3.22e12]. And one more that is easy to miss —
    // THE ANSWER IS EVALUATED IN THE CALLER'S FLOATING-POINT ENVIRONMENT. A consumer that sizes a ring
    // here and asks a MODEL for its latency on another thread gets one number only if both threads
    // round the same way: at hostSR = nextafter(96748.9921875, 0) the geometry is 96.499999999999986
    // to nearest and exactly 96.5 upward, so the two answers are 96 and 97 and a ring sized by the
    // first is one short of a model prepared under the second. Nothing in this repository changes the
    // rounding mode; a consumer that does owes itself a spare slot.
    //
    // A caller sizing a ring wants this PLUS ONE, because a delay line's usable range is capacity-1
    // (core::DryAligner clamps to [0, capacity-1], silently). That "+1" has exactly one reason and
    // belongs to the ring, so it is not folded in here.
    static int maxLatencySamples (double hostSR) noexcept;
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

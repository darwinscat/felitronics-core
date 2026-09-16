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
//   • process() and reset() — what the audio thread calls. 🔴 Neither allocates, locks,
//     does IO or throws — for the architectures named in the next paragraph, which is the
//     same carve-out for both. A freshly-loaded model is atomic-swapped into the live pointer;
//     the replaced model is parked and freed on the message thread (collectGarbage) only
//     once the audio thread has provably stepped past it — so no use-after-free and the
//     audio thread never deletes. process() is O(the block); reset() is NOT — see its own
//     note below: it costs a whole drain length of inference per lane that has played.
//
// The no-allocation guarantee — for process() AND for reset() — applies to architectures whose pinned
// NAM implementation preallocates its work in Reset (Linear and WaveNet). LSTM at this pin returns an
// owning dynamic Eigen hidden-state vector per sample; ConvNet constructs dynamic Eigen temporaries per
// block. Those inherited OrbitCab behaviors are tracked upstream and are not rejected at load,
// preserving byte-compatibility with OrbitCab. reset() INHERITS that carve-out — it runs the same
// process path — rather than adding one, and at this pin it did not spend it: a stereo LSTM restart of
// 48 000 samples allocated nothing through either gate (the house operator-new counter or Eigen's own
// EIGEN_RUNTIME_NO_MALLOC), on the minimal fixture and on a real lstm.nam. So the carve-out is what
// upstream's shape PERMITS, measured as unspent here; a restart is not a new exposure to it.
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
    //
    // 🔴 WHAT IT PROMISES ABOUT STATE, AND IT IS THE SAME SENTENCE reset() MAKES: a prepared stage holds
    // no audio the caller fed. Nothing played before this call can be heard after it — at any rate, on
    // any shape, and whether or not the rate or the block size actually changed. "Prepared" and "just
    // constructed" name one state, so the two verbs of this class do not have to be read against each
    // other. The exceptions are reset()'s three, word for word, because this IS reset(): a recurrent
    // cell, memory the ledger cannot see past its allowance (see reset()), and NAM's own
    // partitioned-FFT clock. (A capture whose conditioner is a model of its own was the second one until
    // P87 put conditioners in the ledger.)
    //
    // It is NOT "silence comes out" and not, in general, bit-identity with a stage prepared a moment
    // ago — the same two qualifications reset() carries below, for the same reasons.
    //
    // ⚠️ WHAT IT COSTS. Where it charges anything, the charge is reset()'s: one lane's whole drain
    // length of inference per lane that may still be holding audio — fed since it was last emptied —
    // 3.77 ms for a real Standard WaveNet at a 64-sample block, per lane. A lane whose falling-edge
    // drain ran to the end HAS been emptied and is not charged (P90: it used to be, a second time) —
    // unless the capture is RECURRENT, which no finite drain empties, so it is charged every time (see
    // reset() below). This call is already the expensive one (it allocates and prewarms the
    // network), it is the message thread's, and it has no callback to miss. A FIRST prepare after a
    // load costs nothing at all — nothing has been fed, so there is nothing owed — and neither does the
    // re-prepare a model swap performs, for the same reason. The lane a permanently mono host never
    // hands over is never charged either. For an architecture whose own Reset already PREWARMS (every
    // WaveNet), this drain is a second pass over the field and roughly doubles the call: a stereo real
    // Standard measured 6.5 ms before and 13.1 ms after at 48 kHz.
    //
    // ⚠️ AND THAT SECOND PASS IS SKIPPABLE ONLY PER ARCHITECTURE — measured, and registered as P98 rather
    // than taken. With the drain removed outright, every capture in NAM's own example set keeps
    // independence at exactly 0 (and its 3.58e-07 residue against a fresh stage vanishes), because a
    // WaveNet's `SetMaxBufferSize` zeroes its rings; this file's own fixtures then leak on 72 of 96 grid
    // cells, worst 0.567861, because a `Buffer`-based capture (Linear, ConvNet) keeps its window through
    // `Reset`. The real example set contains no such capture anywhere in any tree, so a check run only on
    // real captures would have called a blanket skip safe. The predicate belongs to the receptive-field
    // registry, which already walks the tree.
    //
    // ⚠️ AND THE COST NOW SCALES WITH THE HOST RATE, which it did not before and which no caller would
    // guess. The drain is the field converted into HOST samples, so a rate-matched capture is charged
    // `ceil((field + ring + taps) · hostSR/modelRunSR)` of them: measured on a delay(514) capture,
    // stereo, 0.3 ms at 48 kHz, 7.2 ms at 192 kHz, 12.5 ms at 3 MHz, 174 ms at 1e8 Hz and 1.6 s at
    // 1e9 Hz, and at the precondition's own ceiling the length clamps at INT_MAX − 1 per lane. Inside
    // the rates a host offers this is nothing; a caller that prepares this stage at an arbitrary rate
    // is buying inference proportional to it. `rigplayer::RigPlayer` is not exposed to it — its
    // usableSampleRate() substitutes 48 kHz outside (0, 3e6] — but a direct consumer of this class is.
    //
    // 🔴 THIS USED TO BE THE OTHER HALF OF THE DEFECT P47 CLOSED IN reset(), and it is where the number
    // was first measured: a tone, then this call, then digital silence at full width answered
    // 0.224604502320 — −12.97 dBFS — on a dense 2001-tap capture at 48 kHz, because `::nam::DSP::Reset`
    // calls SetMaxBufferSize and then a prewarm that is ZERO samples for a `Linear`, so `Buffer`'s
    // per-channel window survived. There is deliberately no predicate on what changed: a re-prepare at
    // the SAME rate and block is the common case (a host's buffer-size slider moves more often than its
    // rate one), and it is exactly the case that leaked loudest.
    void prepare (double sampleRate, int maxBlock);

    // 🔴 THE STREAM RESTART, AND IT MEANS IT. Audio thread (or any thread with audio stopped): every
    // lane that carried audio is fed the digital silence it still owes, here and in full, until its
    // state is provably the state of a lane that was silent all along.
    //
    // WHAT IS PROMISED IS INDEPENDENCE: nothing the CALLER fed before the restart can be heard after it.
    // That one is EXACT and is what the suite pins — two stages fed different audio before the restart
    // answer the next programme with the same bits, at every rate, on every shape.
    //
    // It is NOT "silence comes out": a capture answers digital zero with whatever its own biases make of
    // it, fresh and restarted alike — measured on the shipped examples at 48 kHz, 0.001195220510 for a
    // real Standard and 9.266554832458 for the A2-max feature set. Exact zero is a property of the
    // bias-free fixtures, which is what lets them witness the promise.
    //
    // And it is not, in general, "bit-identical to a stage prepared a moment ago": NAM's answer depends
    // on how the stream is CUT INTO CALLS (the same property that moves a decaying cell's first sample
    // when `maxModelFrames` changes — see configureRates), and a restart's chunking is its own. Measured
    // against a stage prepared a moment ago: exactly 0 for every fixture in the suite and for a real
    // slimmable WaveNet, and 1.037e-06 for a real Standard at blocks 64…512, where the restart's last
    // chunk is short. That residue is not audio this failed to flush — independence is exactly 0 for the
    // same capture — and feeding the debt in whole prepared blocks removes it, which is a decision about
    // the ledger rather than a patch.
    //
    // ⚠️ IT IS NOT O(THE BLOCK), and it is the only call here that is not. The price is one lane's whole
    // DRAIN LENGTH of inference per dirty lane — the model's field, plus the partitioned-FFT ring any
    // `Linear` capture is charged, plus each rate-matcher leg's own tap window, so it is bigger than
    // prewarmSamples() and that getter is not an estimate of it (a 2001-tap Linear reports 2000 and is
    // charged 4048; an untagged LSTM reports 1 and is charged 24 000). On an M-series core, per lane, a
    // real Standard WaveNet costs 3.77 ms at a 64-sample block — 282 % of that callback — 3.46 ms at
    // 256, 3.43 at 512; a real LSTM 1.3 ms; a dense 2001-tap Linear 0.13 ms — and it grows faster than
    // linearly as the block SHRINKS, because NAM's per-call overhead is paid `debt / maxBlock` times: a
    // real Standard costs 3.61 ms per lane at block 256 and 19.47 ms at block 1. It is allocation-,
    // lock- and throw-free exactly where process() is — the same code path, so the same carve-out, and
    // that carve-out is wider than the architecture names above suggest: `wavenet_a2_max.nam`, a WaveNet
    // in NAM's own example set, allocates 4 times per sample in process() (1024 allocations for one
    // stereo 256-block, measured) and therefore inside a restart too. So it is safe to call from the
    // audio thread and it is NOT free there:
    // at a small block a host that calls it mid-stream buys one late callback. The natural place is
    // where prepareToPlay is — and the debt is re-armed only by audio actually being FED, so a second
    // restart with nothing in between costs nothing, and a mono host pays for one lane, not two. (For a
    // RECURRENT capture it is deliberately NOT idempotent — see below: each restart spends the heuristic
    // again, because such a lane never becomes provably clean.)
    //
    // ⚠️ A RECURRENT CAPTURE IS THE NAMED EXCEPTION (DSP-ARCHITECTURE.md, law 11a): no finite length of
    // silence empties an LSTM cell, so this spends NAM's own half-second heuristic and leaves what that
    // leaves — measured 300 samples differing from a fresh instance, worst 1.49e-07, on a real capture,
    // and 0.419413 on the repository's deliberately slow fixture. A recurrent lane is therefore never
    // marked clean: every restart spends the heuristic again, which is what NAM's own Reset does too.
    //
    // ⚠️ AND IT FLUSHES WHAT THE LEDGER CAN SEE — which now includes a capture's CONDITIONER. A
    // `config.condition_dsp` is a whole model of its own and NAM builds it with `get_dsp` like any
    // other; its memory used to be outside the ledger, because NAM answers 0 for a `Linear` conditioner
    // and `detail::receptiveFieldFromConfig` walked `submodels`, not `condition_dsp`. Such a capture was
    // under-flushed here by exactly as much as law 11a's falling edge under-drained it — measured on a
    // WaveNet with a 2001-sample conditioner, 0.905147969723 after a full drain and 0.905148267746 after
    // a restart. One defect in one ledger, and it was corrected in the LEDGER, so both readers moved
    // together: the registry now adds the conditioner's memory to the network's, in series, and the same
    // fixture answers digital silence with the model's own silence state through either verb.
    // ⚠️ AND WHERE THE LEDGER CANNOT PLACE SOMETHING, IT CHARGES ONE ALLOWANCE — NEVER ZERO, NEVER THE
    // FACE VALUE (P92). The measured case was NAM's slimmable wrapper, whose real config sits under
    // `config.model` where nothing read it: 0 samples flushed for a model reaching 2046, on NAM's own
    // shipped `slimmable_wavenet.nam` rewrapped. An unread config, a value that cannot be read, or a
    // reading the ledger sets aside now adds `detail::kUnreadShapeCeiling` (48 000 samples) ONCE to what
    // it did read, plus the 2048-sample ring — measured through this call at 39.7 ms per lane (256 block)
    // on the most expensive real capture rewrapped, 44.5 ms at 64, and 4.5 ms on the wrapped slimmable.
    // A dead number it CANNOT place cannot inflate it; a dead number it DOES place still can, as before
    // (a wrapped form's own decoy stack) — and a LIVE memory it cannot place, longer than the allowance,
    // drains short by the difference. Both doors are named in `ReceptiveField.h`, and which one stays open
    // is a registered policy question. The rule fires on none of the 1229 distinct captures on the
    // author's machine.
    //
    // ⚠️ AND IT CANNOT REWIND A THIRD PARTY'S CLOCK. NAM's partitioned `Linear` engine counts every
    // sample the instance has ever seen and decides from it where the next programme falls against its
    // partition boundaries; rewinding that means re-configuring the engine, which allocates. The residue
    // is 1.1e-07 against a stage prepared a moment ago and EXACTLY ZERO against one clocked to the same
    // point — the engine's own arithmetic, not state this stage kept.
    void reset() noexcept;

    // 🔴 RT-safe, in place. No model loaded → clean passthrough (no-op, and an ACCEPTED call).
    // `normalize` applies the model's loudness makeup (output normalisation) when the model
    // carries a loudness tag — brings raw model output to a consistent reference level.
    // Law 11: `numSamples` is any length (chunked internally); `numChannels` is 1 or 2, and 0 is the
    // law-11(d) CLOCK-ONLY call — accepted, `io` may be null, and the lanes the caller stopped handing
    // over are fed the digital silence they are receiving (see drainedSamples()). A width above 2 is
    // REFUSED whole — false means nothing was touched. See DSP-ARCHITECTURE.md §2 law 11.
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
    // How many samples this model must be FED before its output means anything — its MEMORY. A network
    // with empty buffers describes the silence it was born into for exactly this long, so anything that
    // fades a freshly loaded model in has to run it silently for this many samples first, and anything
    // that drains a lane the caller stopped feeding owes it the same.
    //
    // Where the number comes from, because three sources disagree and the answer is the largest: NAM's
    // own `GetPrewarmSamples()` answers for ConvNet, LSTM and a plain WaveNet, and a container forwards
    // to its ACTIVE submodel; a SlimmableWavenet answers zero (`wavenet/slimmable.h`), and `Linear`
    // inherits the base class's zero — so the config is read here as well, one dilated convolution at a
    // time, or the plain `receptive_field` where the architecture simply states it. 0 = no model, a
    // capture with no memory (a one-tap gain), or an architecture nothing here can read.
    //
    // ⚠️ It is a MEMORY, not a tap count: an N-tap impulse response reaches back N−1 samples, so a
    // one-tap capture answers 0. This used to answer 0 for EVERY Linear capture, whatever its length.
    int    prewarmSamples()  const;

    // 🔴 THE DRAIN'S ODOMETER — how many samples of digital silence this stage has fed to lanes the
    // caller STOPPED handing over, since it was CONSTRUCTED. It is not cleared by a load, a clear or a
    // prepare: it belongs to the stage, not to the model, so a caller comparing two moments subtracts. A lane the host takes away is not
    // skipped: its network window and its two rate-matchers would otherwise freeze and be replayed on
    // the return (measured 0.518588 out of digital silence at 44.1 kHz, and 0.499533 at 48 kHz with a
    // 2001-tap capture, where no rate-matcher exists at all). It is fed silence instead, for as long as
    // its state can still be heard, and then it STOPS — a permanently mono host pays for one network,
    // not two.
    //
    // This exists to be TESTED, not to be acted on: past the debt the lane's output is zero whether it
    // is still being clocked or not, so "and then it stops" has no witness in the audio and a drain that
    // ran for ever would look identical. Counted per chunk on the audio thread with a relaxed store.
    long long drainedSamples() const;

    // 🔴 THE RESTART'S OWN ODOMETER — how many samples of digital silence `reset()` has spent, since
    // this stage was CONSTRUCTED. Kept apart from drainedSamples() on purpose: the two mechanisms feed
    // silence for different reasons (a lane the caller stopped handing over, against a stream restart),
    // and one number for both would be an oracle that cannot say which of them moved.
    //
    // Like the drain's, this exists to be TESTED: "the full length for a lane that was playing, the
    // remainder for one mid-drain, and NOTHING for a lane that never carried audio" has no witness in
    // the audio — a restart that ran a second network for a mono host sounds exactly like one that did
    // not, and costs 132 ms of a real WaveNet per call.
    long long clearedSamples() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl;
};

} // namespace felitronics::nam

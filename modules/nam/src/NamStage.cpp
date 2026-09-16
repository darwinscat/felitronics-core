// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#include <felitronics/nam/NamStage.h>

#include <NAM/dsp.h>        // NAM_SAMPLE (float, via NAM_SAMPLE_FLOAT) + class ::nam::DSP
#include "ReceptiveField.h"

#include <NAM/get_dsp.h>    // ::nam::get_dsp(path|json)

#include <felitronics/core/StreamResampler.h>
#include <felitronics/neural/NeuralStage.h>   // the shared swap-safe holder extracted from AmpStage

#define NAMZ_IMPLEMENTATION
#include <namz.h>           // load path accepts BOTH raw .nam JSON and packed .namz

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <vector>

namespace felitronics::nam
{

namespace
{
    // Output-normalisation base reference. Each .nam tags its measured output loudness (dB); we
    // bring it to this target so swapping tubes doesn't jump the level. A PER-MODEL trim (passed
    // into the load) is added on top. Final level is still governed by output gain + auto-level.
    constexpr float  kNormTargetDb   = -18.0f;
    inline float dbToGain (float db) { return std::pow (10.0f, db * (1.0f / 20.0f)); }

    // Factory captures are conventionally trained at 48 kHz; the core StreamResampler converts
    // the host stream to/from this so an untagged model still runs at the historical native rate.
    // The one place this number lives is NamStage::kModelSampleRate (public, so consumers can name
    // it). This alias exists only so the lines below stay readable.
    constexpr double kModelSampleRate = NamStage::kModelSampleRate;

    // NeuralStage bounds its retire queue (the pre-extraction code kept an unbounded vector) and
    // REFUSES a swap/clear against a full one. 64 pending models is unreachable in normal use —
    // the queue only backs up while the audio thread is frozen across that many swaps — and the
    // pending-intent slot in NamStage::Impl makes even that case lossless.
    constexpr int kMaxRetiredModels = 64;

//==============================================================================
// NamBackend — the NAM half of the split. The generic model-swap machinery (atomic live pointer,
// block-counter retire, message-thread GC) lives in felitronics::neural::NeuralStage; what stays
// local is everything NAM-specific, shaped as a felitronics::neural::Inference backend. One backend
// instance = one loaded capture: TWO independent ::nam::DSP instances of the SAME capture (true
// stereo — L/R independent; a mono stream uses just instance 0), the −18 dB loudness makeup
// (per-model trim folded in), and the per-channel StreamResampler rate-match state + scratch.
// NeuralStage swaps WHOLE instances of this class, so the resampler/scratch state travels with the
// model it belongs to.
class NamBackend
{
public:
    // Both instances always exist (built by the loader) and prepare() configures both per-channel
    // resamplers, so a host switching the layout mono<->stereo (a re-prepare) is safe.
    // The normalize flag is owned by NamStage::Impl: the audio thread stores the per-call `normalize`
    // argument there right before NeuralStage::process() reaches this backend (same thread →
    // sequenced), keeping NamStage's per-block flag without widening the shared Inference seam. It is
    // BOUND when the backend meets its stage (bindNormalize, before it goes live) — a backend can be
    // built and prepared with no stage in sight; unbound, it plays raw.
    NamBackend (std::unique_ptr<::nam::DSP> m0, std::unique_ptr<::nam::DSP> m1,
                float trimDb, int prewarmFromConfig = 0, int drainTailFromConfig = 0,
                bool recurrentFromConfig = false)
    {
        prewarm    = prewarmFromConfig;       // …raised below by NAM's own answer, where it has one
        drainTail_ = drainTailFromConfig;     // the partitioned-FFT ring, which is not the field
        recurrent_ = recurrentFromConfig;
        inst[0] = std::move (m0);
        inst[1] = std::move (m1);

        expectedSR = inst[0]->GetExpectedSampleRate();
        prewarm = std::max (prewarm, inst[0]->GetPrewarmSamples());
        if (inst[0]->HasLoudness())
        {
            loudnessDb  = inst[0]->GetLoudness();
            hasLoudness = true;
            makeup = dbToGain ((float) (kNormTargetDb - loudnessDb) + trimDb);
        }
        else
        {
            loudnessDb = 0.0; hasLoudness = false;
            makeup = dbToGain (trimDb);
        }
    }

    // Retire-ledger hook (see NamStage::Impl::retiredLedger). Attached ONLY when this instance is
    // handed to a swap that goes live — a superseded, never-live pending instance stays detached,
    // so its death doesn't skew the count. The dtor runs on the message thread only (NeuralStage
    // frees retired backends in collectGarbage / its own teardown, never on the audio thread).
    void attachRetireLedger (int* ledger) noexcept { retireLedger = ledger; }
    // Message thread, on a backend that is NOT live (the loader, before the swap).
    void bindNormalize (const std::atomic<bool>& flag) noexcept { normalize = &flag; }
    // The drain's own odometer, published for a test to READ rather than infer. "It drains, and then it
    // STOPS" has no other witness: past the debt the lane's output is zero either way, so a mutation that
    // never stops draining is invisible in the audio and shows up only as CPU. Counted per CHUNK, not per
    // sample — a counter inside the sample loop is how a gate killed vectorisation once already.
    void bindDrainCounter (std::atomic<long long>& counter) noexcept { drained_ = &counter; }
    // …and the RESTART's own, kept apart from the drain's on purpose: the two mechanisms spend silence
    // for different reasons (a falling edge against a stream restart), and one counter for both would be
    // an oracle that cannot say which of them moved. See NamStage::clearedSamples().
    void bindClearCounter (std::atomic<long long>& counter) noexcept { cleared_ = &counter; }
    double preparedSampleRate() const noexcept { return hostSR; }
    int    preparedMaxBlock()   const noexcept { return maxBlock; }
    ~NamBackend() noexcept { if (retireLedger != nullptr) --(*retireLedger); }

    //--- felitronics::neural::Inference seam --------------------------------------
    // Message thread, with audio stopped (live instance, via NeuralStage::prepare) or on a
    // not-yet-live instance (the loader). Decide the run rate, (re)configure the per-channel
    // resamplers/scratch, and Reset both instances to it (alloc + prewarm — fine off the audio
    // thread). Both lanes are always configured — see the ctor note on mono<->stereo re-prepare.
    //
    // The Inference concept REQUIRES prepare() to be noexcept, yet the body genuinely allocates
    // (resampler/scratch resizes + the NAM Reset prewarm) — sizes depend on the host's future
    // sampleRate/maxBlock, so they can't be preallocated in the ctor. A low-memory failure is
    // therefore caught HERE instead of std::terminate()-ing the host: the backend marks itself
    // unprepared and latencySamples() drops to 0. The loader checks prepared() and fails the load
    // honestly; a LIVE backend whose re-prepare fails leaves the caller's audio ALONE.
    //
    // 🔴 AND IT REFUSES THAT CALL RATHER THAN CLAIMING IT — the two halves of that sentence are not the
    // same and an earlier wording here said only the second ("degrades to a clean passthrough"), which
    // reads as "returns true". It does not: process() returns FALSE for an unprepared backend, which is
    // law 11's answer for a call that cannot be honoured, and a composite that ANDs its stages'
    // verdicts (rigplayer does) therefore reports the refused block outward. The buffer is untouched
    // either way; what the verdict adds is that the caller can tell.
    void prepare (double sampleRate, int maxBlockIn, int /*maxChannels*/) noexcept
    {
        prepared_ = false;
        try
        {
            hostSR   = sampleRate;
            maxBlock = std::max (1, maxBlockIn);
            // 🔴 THE RAW RATE, NOT A NORMALISED ONE. This line used to spell
            // `expectedSR > 0.0 ? expectedSR : kModelSampleRate` — a second copy of the normalisation
            // that rateMatch() already owns — and the mutation stand proved it was not decoration:
            // with the default pre-applied here, breaking the default INSIDE rateMatch changed nothing
            // for an untagged model, so the suite could not see it. Hand over what the model actually
            // reported and let the one owner decide.
            if (! configureRates (expectedSR))
                return;                    // prepared_ stays false — see configureRates()
            for (auto& m : inst)
                if (m) m->Reset (modelRunSR, maxModelFrames);
            // 🔴 AND GROW WHAT NAM GROWS LAZILY, HERE, WHERE ALLOCATING IS LEGAL. `Buffer::_update_buffers_`
            // resizes its per-channel window on demand inside process(), and `Reset` pre-grows it only via
            // prewarm() — which runs `GetPrewarmSamples()` samples, and that is ZERO for a Linear capture.
            // So the first audio call on a lane allocated, and until this task that was invisible because
            // instance 1 was never touched on a mono host. The drain touches it, which turned a dormant
            // hole into a live one: measured 4 allocations in the first width-1 call against 2 with the
            // drain off, i.e. two per instance. One block of silence through each such instance grows
            // everything the audio path can ask for — it is exactly what NAM's own prewarm does for the
            // architectures that declare one, and the state it leaves is the silence state either way.
            for (auto& m : inst)
                if (m != nullptr && m->GetPrewarmSamples() <= 0)
                {
                    std::fill (ch[0].modelIn.begin(), ch[0].modelIn.end(), 0.0f);
                    NAM_SAMPLE* in [1] { ch[0].modelIn .data() };
                    NAM_SAMPLE* out[1] { ch[0].modelOut.data() };
                    m->process (in, out, maxModelFrames);
                }
            prepared_ = true;
        }
        catch (...) {}   // bad_alloc (or a throwing NAM Reset): stay unprepared, never crash
        // 🔴 AND THE PREPARATION ENDS IN THE RESTART — ALWAYS, not only when one was asked for. This is
        // the whole of P85 and it is three lines because P24 and P47 already wrote the mechanism: the
        // debt `configureRates` just re-charged above is exactly what this stream left in the network,
        // and reset() is the one thing that spends it. A prepare() is where a host restarts a stream —
        // every path into this function is a host changing its rate or its block size, and a driver
        // stops the stream to do either — so "prepared" and "just constructed" now mean the same state,
        // and the two verbs of this class say one thing. What it replaces is 0.224604502320 out of
        // DIGITAL SILENCE (−12.97 dBFS) on a dense 2001-tap capture: `::nam::DSP::Reset` calls
        // SetMaxBufferSize and a prewarm that is ZERO samples for a `Linear`, so `Buffer`'s per-channel
        // window survived the call — and the grow loop above pushes `maxModelFrames` zeros, which at a
        // 64-sample block is 80 against a 2001-sample window.
        //
        // NO PREDICATE ON "WHAT CHANGED", deliberately: a re-prepare at the SAME rate and block is the
        // common case (a buffer-size slider moves more often than a rate one) and is exactly where the
        // old stream's tail was loudest — measured 0.468718945980 there, against 0.469410002232 across
        // a rate change. A predicate would have fixed the rarer case and left the frequent one.
        //
        // IT IS FREE WHERE IT IS CALLED FROM ANYWHERE ELSE, which is what makes this safe to put here
        // rather than at the four call sites. A backend reached by `prepareModel()`, by `install()`'s
        // re-prepare, or as `NamStage::prepare()`'s parked `pendingBackend` has never been fed, so
        // `everFed_` is false on both lanes, `configureRates` charged it nothing, and this spends
        // nothing — a model change goes through prepare() once, or twice when `install()` finds the host's
        // numbers moved between the halves, and pays for neither. The one backend that can be dirty is
        // the LIVE one, and it is reached exactly once per host prepare.
        //
        // ⚠️ WHERE IT IS NOT FREE IT IS NOT HIDDEN EITHER: for an architecture whose own Reset already
        // prewarms (every WaveNet), this drain is a SECOND pass over the field and roughly DOUBLES the
        // call — a stereo real Standard measured 6.8 ms before and 12.6 ms after at 48 kHz. That is the
        // message thread, with no callback to miss, and it is the price of one rule instead of two.
        // Skipping it where NAM's own prewarm provably covers the ledger is a real optimisation and is
        // registered as one rather than taken here: it is a predicate, and it would put back exactly
        // the kind of "this case does not need it" reasoning this fix exists to remove.
        //
        // OUTSIDE the try on purpose: reset() is noexcept, so a throw inside it terminates rather than
        // landing in that catch, and putting it there would suggest otherwise. Skipped when this prepare
        // REFUSED, because chunking by a `maxBlock` a refused prepare already stored is a heap-buffer-
        // overflow (see reset()); the next prepare that succeeds restarts instead, owed or not — which
        // is why the `restartOwed_` bit P47 needed here is gone rather than kept as a second answer.
        //
        // ⚠️ AND THAT IS A NEW EXPOSURE ON THIS PATH, named rather than hidden: the drain runs NAM's own
        // `process()`, and for the architectures whose process() allocates (`wavenet_a2_max.nam`, LSTM
        // and ConvNet temporaries — the same carve-out reset() already inherits) a `bad_alloc` here
        // TERMINATES, where the identical failure a few lines above, inside `Reset`/the prewarm/the grow
        // loop, is caught and leaves the backend honestly unprepared. It is exactly the exposure those
        // captures already carry on the audio thread every block, so it is not a new class; what is new
        // is that a host's prepare() can now meet it. Catching it here would be worse than the risk: a
        // half-spent drain caught and swallowed leaves a lane the ledger calls clean and the network
        // does not, which is the one thing this whole verb exists to make impossible.
        if (prepared_)
            reset();
    }

    bool prepared() const noexcept { return prepared_; }

    // 🔴 RT-safe, in place — the audio thread reaches this via NeuralStage::process() on the
    // live instance. Never allocates, locks, does IO or throws. An unprepared backend (failed
    // low-memory prepare — see above) is a clean passthrough.
    // LAW 11: the length is a CAPACITY, so a call longer than maxBlock is CHUNKED, not clamped. It used
    // to be `n = std::min (numSamples, maxBlock)`, and the tail past maxBlock then came out of the amp
    // BIT-IDENTICAL TO ITS INPUT — measured on a Linear model prepared for 64 and called with 512: 448
    // of 512 samples never met the model at all. Simply dropping the clamp would have been far worse than
    // the defect: `processChannel` copies n samples into `modelIn`, which is `maxModelFrames` long, and
    // NAM's own buffers are sized from the prepared block too — an unclamped n is a heap overflow here
    // and a resize (an allocation, on the audio thread) inside NAM.
    //
    // The loop lives HERE and not in NeuralStage::process on purpose: the live model is resolved ONCE
    // per host call, above us. Chunking a level up would re-resolve it per chunk, so a model swap landing
    // mid-buffer could put half a block through the old capture and half through the new one.
    // This is exactly what rigplayer::RigPlayer already does around its own NamStage instances.
    [[nodiscard]] bool process (float* const* io, int numChannels, int numSamples) noexcept
    {
        if (numChannels < 0 || numSamples < 0) return false;
        if (! prepared_) return false;
        if (numChannels > 2) return false;                 // a NAM capture is mono or true-stereo; nothing else
        if (numSamples == 0) return true;                  // law 11(d): no samples, no time, no edge

        const float g = (normalize != nullptr && normalize->load (std::memory_order_relaxed)) ? makeup : 1.0f;
        for (int off = 0; off < numSamples; )
        {
            const int n = std::min (numSamples - off, maxBlock);
            for (int c = 0; c < 2; ++c)                                        // mono → 1 instance; stereo → 2 independent ones
            {
                if (c < numChannels)
                {
                    processChannel (ch[c], inst[c].get(), io[c] + off, n, g);
                    drain_[c] = drainSamples_;             // a lane that is playing owes a full drain when it stops
                    everFed_[c] = true;                    // …and is a lane a later prepare must charge
                }
                else if (drain_[c] > 0)
                {
                    // 🔴 THE PAUSE IS SILENCE, AND FOR THIS STAGE THAT IS PER LANE — law 11c/11a. A lane the
                    // host stops handing over used to be SKIPPED, so its network window and its two
                    // rate-matchers froze and were replayed on the return: measured through
                    // rigplayer::RigPlayer, 0.518588 out of DIGITAL SILENCE at 44.1 kHz (-5.70 dBFS, 125
                    // host samples) with a memoryless capture, and 0.499533 with a 2001-tap one at 48 kHz,
                    // where no rate-matcher is installed at all — the two halves are independent and the
                    // second is the one a 48 kHz fixture cannot see. Feeding it the digital silence it is
                    // actually receiving takes both to EXACTLY zero (measured on the whole audio rate grid).
                    // The scratch is this stage's own: `io` need not have a plane here at all, and at
                    // numChannels == 0 it may legally be null.
                    // BOUNDED, not forever: past drainSamples_ the lane's state IS the silence state, so
                    // further zeros change nothing a caller can hear and the expensive part — the network —
                    // stops. That is what keeps a permanently mono host at one model instead of two.
                    const int d = std::min (n, drain_[c]);
                    std::fill (hush_.data(), hush_.data() + d, 0.0f);
                    processChannel (ch[c], inst[c].get(), hush_.data(), d, g);
                    drain_[c] -= d;
                    if (drained_ != nullptr) drained_->fetch_add ((long long) d, std::memory_order_relaxed);
                    // 🔴 …AND A DRAIN THAT RAN ALL THE WAY MAKES THE LANE CLEAN IN THE LEDGER TOO, which
                    // until P90 it did not. `everFed_` was set by every fed chunk above and cleared ONLY
                    // by reset(), so a lane that had just spent its whole debt here was still marked
                    // "may be holding audio" and `configureRates` charged it a FULL drain again at the
                    // next prepare() — the stage paying twice for one departure. P85 registered it as a
                    // cost and measured it on a fixture (2562 spent, 5124 charged); on the real captures
                    // the wasted lane is 4093 samples for `wavenet_a1_standard`, 6347 for `A2` and 2047
                    // for `slimmable_wavenet`, and the re-prepare that spends it measures 14.7 ms
                    // against the 7.4 ms the same call costs with nothing owed.
                    //
                    // IT IS THE SAME CLAIM reset() ALREADY RESTS ON, not a new one: `drain_[c] == 0`
                    // means exactly "this lane has been fed the silence it owed", and reset() reads that
                    // very counter (`owed = drain_[c]`) to decide it has nothing to spend. What reset()
                    // does that this path does not is clear the two rate-matcher legs — and that is NOT
                    // a hole here, for two reasons. The first is that this is a PAUSE and not a restart
                    // (law 11c): the stream continued, with zeros, so the legs' sub-sample phase is
                    // exactly where it would have been, and for a lane that comes BACK that continuity
                    // is the correct state — re-anchoring it here would be the 1.039e-06 divergence
                    // reset() documents, introduced deliberately. The second is that it does not matter
                    // to the only reader: `configureRates` has already re-derived both legs, coefficients
                    // AND state, before it consults `everFed_` at all, and reset() clears them
                    // unconditionally whether or not it spends a sample.
                    //
                    // A RECURRENT CAPTURE KEEPS ITS FLAG, and the exclusion is copied from reset()'s own
                    // line rather than reasoned about again: no finite length of silence empties an LSTM
                    // cell, so `drain_ == 0` is not cleanliness for one and it never becomes provably
                    // clean. Measured, the lstm fixture drains 24000 here and is charged 24000 again —
                    // and must be.
                    if (drain_[c] == 0 && ! recurrent_) everFed_[c] = false;
                }
            }
            off += n;                                      // `off += maxBlock` could step past INT_MAX
        }
        return true;
    }

    // 🔴 THE STREAM RESTART, AND IT IS THE ONE VERB HERE THAT KEEPS ITS WORD. It used to be empty, and
    // the comment that stood here recorded why with a number: `::nam::DSP::Reset` calls
    // SetMaxBufferSize and then prewarm(), and `Linear` overrides neither the prewarm (the base class
    // answers 0) nor the input buffer — `Buffer`'s per-channel window survives. Measured on a dense
    // 2001-tap capture at 48 kHz: a tone, then `NamStage::prepare()`, then digital silence at full
    // width returns 0.224604502320; through `reset()`, the same. P24 closed the half that was a falling
    // edge — a lane the host STOPS handing over is fed the silence it is receiving — and left this half,
    // the lane that is PRESENT, whose stale window speaks into the caller's own samples.
    //
    // WHAT IT PROMISES — the whole statement, with its numbers and its three exceptions, is on the public
    // declaration in NamStage.h. In short: every lane that carried audio is fed the digital silence it
    // still owes, here and in full, until its state is provably the state of a lane that was silent all
    // along. The promise is INDEPENDENCE — nothing the CALLER fed before the restart can be heard after
    // it — which is exact and is what the suite pins. It is NOT "silence comes out" (a capture answers
    // digital zero with its own biases: 0.001195220510 for a real Standard, 9.266554832458 for the
    // A2-max feature set), and it is not, in general, bit-identity with a stage prepared a moment ago,
    // because NAM's answer depends on how the stream is cut into calls: 1.037e-06 on a real Standard
    // where this restart's last chunk is short, and exactly 0 on a real slimmable at the same blocks.
    //
    // WHAT IT COSTS, because this is the one house verb here that is not O(block): one lane's whole
    // DRAIN LENGTH of inference per dirty lane — the field, plus the FFT ring a Linear is charged, plus
    // each leg's tap window. Real Standard WaveNet on an M-series core, per lane: 3.77 ms at a
    // 64-sample block (282 % of the callback), 3.62 at 128, 3.46 at 256, 3.43 at 512, 3.33 at 1024;
    // the A2 container 3.0-4.0 ms; a real LSTM 1.3 ms; a dense 2001-tap Linear 0.13 ms. It allocates,
    // locks, blocks and throws exactly where `process()` does — i.e. nowhere, for the architectures
    // whose NAM implementation preallocates — which is NOT every Linear and WaveNet: `wavenet_a2_max.nam`
    // from NAM's own examples allocates 4 times per sample in process() and therefore in here too
    // (measured), beside the LSTM/ConvNet exception the header already names. There is no cheaper
    // exact mechanism: asking NAM to Reset with its prewarm off zeroes the Conv1D rings in 0.014 ms and
    // still misses the prepared state by 4089 samples (worst 0.324) because that state is a PREWARMED
    // one, and on a Linear with the FFT engine it allocates 46 times.
    //
    // AND IT IS IDEMPOTENT, which is what keeps that price a one-off: the debt is re-armed only by audio
    // actually being fed, so a second restart with nothing in between spends nothing, and a mono host
    // pays for one lane rather than two.
    //
    // ⚠️ WHAT IT DOES NOT REACH, both named with numbers in the header: a recurrent cell (nothing finite
    // empties one) and NAM's own partitioned-FFT clock (1.1e-07 against a stage prepared a moment ago,
    // exactly zero against one clocked to the same point). A capture whose CONDITIONER is a model of its
    // own used to belong on that list — its memory was outside the ledger, the same hole that
    // under-drained it at a falling edge, 0.905147969723 either way — and no longer does: the ledger
    // counts it, so both readers of it do.
    void reset() noexcept
    {
        // NOTHING, on a backend whose preparation was refused — not even the ledgers. `prepare()` writes
        // `hostSR`/`maxBlock` before `configureRates` can refuse (law 11b's disarm-first leaves the
        // object unusable, not consistent), so `maxBlock` can already be the NEW value while `hush_` and
        // `modelIn` are still the old ones: chunking by it would `std::fill` past the end of `hush_` and
        // copy past the end of `modelIn` — measured under ASan as a heap-buffer-overflow, a WRITE 0
        // bytes past the 1024-byte scratch. A first prepare that refuses leaves `hush_` empty outright.
        // The ledgers must stay too: they are the only record the next SUCCESSFUL prepare re-charges the
        // lanes from.
        //
        // 🔴 AND THE REQUEST NEED NOT BE REMEMBERED, which it was until P85 — a `restartOwed_` bit set
        // here and read at the end of prepare(). Every successful prepare now restarts whether one was
        // asked for or not, so the bit had one answer for two questions and is gone. The sequence it was
        // written for is still closed, and by the stronger rule: "play, a refused prepare, reset(), a
        // prepare that succeeds, play again" used to hand back ~242 samples of a delay(514) capture (the
        // rest of its window having been pushed out by prepare()'s own scratch warm-up), and the
        // succeeding prepare flushes it now with no reset() in the sequence at all.
        if (! prepared_)
            return;

        for (int c = 0; c < 2; ++c)
        {
            // A lane with no instance cannot be flushed, so its ledger is left exactly as it is rather
            // than marked clean — the same rule as the unprepared backend above. UNREACHABLE as the
            // class stands (prepareModel refuses a model missing either instance, so a prepared backend
            // has both), and kept honest rather than asserted, exactly like the swap refusal in
            // tryApplyPending below. The mutation stand says so too: removing it changes nothing any
            // test can see, which is the signature of a branch no caller can enter.
            if (inst[c] == nullptr)
                continue;

            // 🔴 THE DEBT IS THE DIRT — for a capture whose memory is FINITE. `drain_[c]` is exactly how
            // much more silence this lane needs before its state is the silence state: the full length
            // while it is playing (process() re-arms it on every chunk it feeds), the remainder while it
            // is mid-drain, and zero for a lane that has never carried audio — which is why a mono host
            // does not run a second network here, the hole `everFed_` closed for prepare().
            //
            // ⚠️ AND FOR A RECURRENT CELL IT IS NOT THE DIRT, which is the one place this reading breaks
            // and three crew seats went for it: an LSTM lane that has already SPENT its drain reads
            // `drain_ == 0` and is still not empty — the fixture leaves 0.419413 there. So a recurrent
            // lane that ever played is charged the WHOLE heuristic again at every restart, which is also
            // what NAM's own Reset does (prewarm, unconditionally, half a second of it), and it never
            // becomes provably clean: `everFed_` stays set for one, because nothing finite empties it.
            int owed = recurrent_ ? (everFed_[c] ? drainSamples_ : 0) : drain_[c];
            while (owed > 0)
            {
                const int d = std::min (maxBlock, owed);
                // REFILLED EVERY CHUNK, AND EVERY LANE: `processChannel` writes the model's OUTPUT back
                // into what it was handed, so a buffer filled once would feed the network its own answer
                // from the chunk before. The gain is 1.0f and not the live makeup — every sample of this
                // goes into the scratch and nothing reads it.
                std::fill (hush_.data(), hush_.data() + d, 0.0f);
                processChannel (ch[c], inst[c].get(), hush_.data(), d, 1.0f);
                owed -= d;
                if (cleared_ != nullptr) cleared_->fetch_add ((long long) d, std::memory_order_relaxed);
            }
            // …AND THE RATE-MATCHERS ARE RE-PRIMED, not re-designed. A restart re-anchors the audio-time
            // clocks — the same reason `eq::EqBand::reset()` re-anchors its StateGrid — and for these two
            // legs that clock is the sub-sample phase of the model grid: leave it and the stage answers
            // the next programme through a different fractional alignment than a fresh one, measured
            // 1.039e-06 at 44.1 kHz. `reset (rates, capacity)` is the wrong tool for it here (it
            // re-derives 513 x 64 coefficients, 1.77 ms, and is not noexcept); `clearAudioState()` is
            // the state alone. Both legs, whether or not one is installed: an identity leg is inert, and
            // a clear that depends on which path ran is a clear with two answers.
            ch[c].down.clearAudioState();
            ch[c].up  .clearAudioState();
            drain_[c] = 0;
            if (! recurrent_) everFed_[c] = false;
        }
    }

    // Host-rate latency the rate-matcher introduces (0 when not resampling).
    //
    // 🔴 THIS LINE HAS BEEN WRONG TWICE, BOTH TIMES BY RESTATING SOMEBODY ELSE'S ARITHMETIC. It was
    // `ceil(3*hostSR/modelRunSR) + 3` once — a guess at "~3 samples of lookahead per stage", 2.16
    // samples out — and P32 replaced it with the cubic kernel's real geometry, which P34 then made
    // stale again by swapping the kernel. So it does not compute anything now: `configureRates()`
    // asked `NamStage::rateMatch()` once, and this reports what it was told. The derivation lives
    // where the geometry lives, in core::StreamResampler; the gate and the rounding live in
    // rateMatch(); and there is exactly one copy of each.
    //
    // The number IS acted on outside this repository — OrbitCab and orbit-amp delay their dry/bypass
    // path by it — so it is an audio-alignment figure there, not only PDC.
    int latencySamples() const noexcept
    {
        return prepared_ ? rm_.latencySamples : 0;   // an unprepared backend passes through
    }

    //--- model info (read by the loader for NamStage's UI-mirror atomics) ----------
    double reportedSampleRate()  const noexcept { return expectedSR; }
    int    reportedPrewarm()     const noexcept { return prewarm; }
    double reportedLoudnessDb()  const noexcept { return loudnessDb; }
    bool   reportedHasLoudness() const noexcept { return hasLoudness; }

private:
    // Decide the run rate + (re)configure per-channel resamplers/scratch (prepare() only — never
    // while this instance is live and audio runs). modelRunSR = the loaded model's native rate.
    // FALSE = this configuration cannot be honoured; the caller must leave the backend unprepared.
    [[nodiscard]] bool configureRates (double modelSR)
    {
        // ONE call decides all three: the run rate, whether a resampler is installed, and what it
        // costs. Recomputing any of them here is how they drifted apart before.
        rm_         = NamStage::rateMatch (hostSR, modelSR);
        modelRunSR  = rm_.modelRunSR;
        resampling  = rm_.resampling;
        // 🔴 THE SCRATCH SERVES TWO PATHS, AND ONE RATIO CANNOT SIZE BOTH. It used to be
        // `ceil(maxBlock * modelRunSR / max(8000, hostSR)) + 16`, and that number was wrong at BOTH
        // ends — each way silent, each way measured:
        //
        //   • TOO SMALL, past the end of the heap. Without a resampler the model is clocked by the
        //     HOST: processChannel() copies a whole chunk — up to maxBlock host samples — into this
        //     buffer, and the ratio never enters. A tag half a hertz BELOW the host is inside the
        //     acceptance window (the gate and the resampler gate are the same half hertz), so the
        //     ratio is just under 1 and the buffer comes out just under maxBlock. Measured with ASan
        //     on the base commit, host 48000, tag 47999.5, maxBlock 2000000: heap-buffer-overflow,
        //     a WRITE of 8000000 bytes into a 7999984-byte region. The threshold is maxBlock >=
        //     34*hostSR — 34 seconds of audio in one call, which law 11(a) explicitly invites an
        //     offline caller to pass and which this very function's own comment (see process()) says
        //     it is guarding against.
        //   • TOO SMALL AGAIN, and inaudibly at first, below 8 kHz: `max(8000, hostSR)` is an
        //     ASSUMED host rate standing in for the real one, so a slower host's block converts to
        //     more model frames than fit and produceAvailable() drops the surplus without a word.
        //     Measured on a unity model, 100 Hz tone: -1.22 dB at a 6 kHz host, -3.00 at 4 kHz,
        //     -6.05 at 2 kHz, -9.09 at 1 kHz, and exactly 0.00 at 8 kHz and above — the floor's own
        //     edge. rigplayer::RigPlayer::usableSampleRate accepts every one of those rates.
        //
        // Ask the path that runs for its own number, with the REAL host rate. The bound
        // that keeps the arithmetic honest is not a floor on the rate but a refusal: a frame count
        // that will not fit in an int cannot be honoured, and this backend already has an observable
        // way to say so — stay unprepared, which makes the loader fail the load and a live instance
        // degrade to a clean passthrough. That is law 11(b) at the only place in this class that can
        // still refuse; a silent substitute rate is what it replaces.
        // ASK THE PATH THAT WILL ACTUALLY RUN — `resampling` is decided above and does not change until
        // the next prepare(), so exactly one of these two numbers is the requirement and the other is
        // irrelevant. Taking the max of both instead looks safer and is not: at a non-positive host
        // rate the ratio is not a number, and a max() would quietly hand the CONVERTED path the direct
        // path's answer — smaller than what the old floored expression gave, which is a regression
        // wearing a fix's clothes. Here that case has no honest number, so it is refused instead.
        //
        // 🔴 THE CONVERTED BRANCH NEEDS ITS HOST RATE STATED, NOT INFERRED FROM THE RESULT. A guard that
        // only looks at the frame count lets a NEGATIVE host through whenever the slack outweighs the
        // (negative) converted term: at host -2000000, block 512, `ceil(512 * 48000 / -2000000) + 16`
        // is 4 — a positive, perfectly representable four-frame scratch for a rate that cannot convert
        // anything. -48000 happens to give -496 and is refused, which is exactly how a wrong rule looks
        // right on the cell you tried. Found by a pre-merge diff round; the claim it falsified was mine.
        if (rm_.resampling && ! (hostSR > 0.0))
            return false;
        // Keep +16 on the direct path: maxModelFrames is also NAM's Reset block, and whole-block
        // prewarming makes its value observable with a stateful capture. This preserves the old block
        // when model and host rates are equal. It does NOT preserve every integer-tagged capture at
        // every host: host 47999.75, tag 48000, block 512 used 529 and now uses 528. A decaying-cell LSTM
        // then starts at 0.1600718498 instead of 0.1597609967. The general direct-path change is
        // maxBlock - ceil(maxBlock * (modelRunSR / hostSR)), not always one frame and not restricted to
        // tags above the nominal rate. See review-p38/REPORT.md, plan 1.
        const double frames = rm_.resampling
                                  ? std::ceil ((double) maxBlock * (modelRunSR / hostSR)) + 16.0
                                  : (double) maxBlock + 16.0;
        // The ceiling is not INT_MAX but (INT_MAX-16)/2, because this number is doubled and offset
        // again on the next line (`maxModelFrames * 2 + 16`) — a bound that only keeps its own
        // conversion legal is a bound that overflows one line later. maxBlock rides the same
        // arithmetic through `down.reset`, so it is held to the same ceiling.
        constexpr double kMaxFrames = (double) ((std::numeric_limits<int>::max() - 16) / 2);
        if (! (frames >= 1.0 && frames <= kMaxFrames && (double) maxBlock <= kMaxFrames))
            return false;                 // refused — the caller leaves this backend unprepared
        maxModelFrames = (int) frames;
        for (auto& c : ch)
        {
            c.down.reset (hostSR,    modelRunSR, maxBlock * 2 + 16);
            c.up  .reset (modelRunSR, hostSR,    maxModelFrames * 2 + 16);
            c.modelIn .assign ((size_t) maxModelFrames, 0.0f);
            c.modelOut.assign ((size_t) maxModelFrames, 0.0f);
        }
        // HOW LONG AN ABSENT LANE IS FED SILENCE — see process(). It is the point past which the lane's
        // state is the state of a lane that was silent all along, and it is a SUM over the three things
        // that hold samples, each counted in ITS OWN rate:
        //
        //     direct path      D = prewarm                                        (host samples)
        //     rate-matched     D = kTaps + ceil((prewarm + kTaps) · hostSR/modelRunSR)
        //
        // The leading kTaps flushes the DOWN leg in HOST samples; the ceil term drives `prewarm` model
        // frames of zeros through the network and another kTaps through the UP leg, both MODEL frames,
        // hence the conversion.
        //
        // 🔴 WHAT A LEG ACTUALLY NEEDS IS kTaps − 1 = 63 OF ITS OWN INPUTS, and this line said 32 until a
        // crew round did the arithmetic. produceAvailable emits the head at `i` only while
        // `i < len − kHalf`, and the window is `buf[i−kBehind … i+kHalf]`, so the LAST head that still
        // reads input sample Q is `i = Q + kBehind` and it runs only once `len ≥ Q + kBehind + kHalf + 1
        // = Q + 64`. Sixty-three further inputs, not thirty-two — so kTaps is ONE sample of margin here,
        // not a doubling, and the honest reading of this term is "the tap window, rounded up to the
        // power of two it is built from". Measured over 96 cells (three shapes x eight rates x blocks
        // 7/64/256/1024) the slack is 1 sample at 8 k / 22.05 k / 44.1 k and up to 8 at 192 k, and the
        // need is block-INDEPENDENT.
        //
        // A third kTaps of slack stood at the end of this expression and is gone: the mutation stand
        // could not tell it from nothing, and a border with no witness is a number nobody can reproduce.
        // What the stand CAN tell is the leading term — drop it and a memoryless capture at a host rate
        // well below the model's drains ceil(kTaps·h/m) samples, which is under 63 the moment h/m falls
        // below about 0.98, which is why 8 kHz and 22.05 kHz are in the suite's rate grid.
        // `prewarm` is the model's receptive field (detail::receptiveFieldFromConfig, raised by NAM's own
        // GetPrewarmSamples) — the same number a caller warms a fresh model for, asked once, here.
        // ⚠️ For an LSTM this is a BOUND ON THE HEURISTIC, not on the memory: a recurrent cell never
        // reaches the silence state exactly, and NAM's own answer there is half a second of samples.
        {
            constexpr double taps = (double) felitronics::core::StreamResampler::kTaps;
            // 🔴 A RECURRENT CELL HAS NO FLUSH LENGTH, so what is spent for one is NAM's own half-second
            // heuristic — taken at the rate the model is RUN at rather than at the tag it reports, which
            // it may not have. An untagged LSTM answers `GetPrewarmSamples() == 1` and a drain sized
            // from that is ONE SAMPLE: measured 0.499222 out of digital silence, where the same model
            // WITH a tag reads 0.419413 and a lane clocked through the whole gap reads 0.022842. So this
            // closes the hole in the heuristic — untagged now behaves as tagged — and NOT the gap between
            // either of them and the truth: half a second is 1.1 time constants of a 22 000-sample cell,
            // and no finite drain closes a recurrent state. The test says exactly that rather than
            // promising a zero the arithmetic cannot keep.
            const double field = recurrent_ ? std::fmax ((double) prewarm, 0.5 * modelRunSR)
                                            : (double) prewarm;
            const double d = rm_.resampling
                                 ? taps + std::ceil ((field + (double) drainTail_ + taps) * (hostSR / modelRunSR))
                                 : field + (double) drainTail_;
            constexpr double kMaxDrain = (double) (std::numeric_limits<int>::max() - 1);
            drainSamples_ = (d >= 1.0) ? (int) std::min (d, kMaxDrain) : 0;
        }
        // 🔴 A PREPARE IS A DEPARTURE OF EVERY LANE THAT WAS PLAYING, so those owe a full drain afterwards
        // — NOT zero, which is what "the Reset cleared it" would imply and what an earlier draft of this
        // line said. `DSP::Reset` does not clear a Linear capture's window at all (`Buffer::_input_buffers`
        // survive `SetMaxBufferSize`): measured on the delay fixture, a tone then `prepare()` then silence
        // replays 0.500000, and on a dense 2001-tap kernel 0.224604502320.
        // 🔴 AND SINCE P85 THIS LINE IS WHAT THE TAIL OF prepare() SPENDS. The debt charged here is the
        // whole of what the old stream left, at the NEW rates, and prepare() ends by feeding it (see the
        // reset() at the end of prepare()) — so the numbers above are what this used to hand back and no
        // longer does. The charge still has to happen HERE rather than inside that restart, because the
        // length is a property of the configuration this function has just decided; and it still has to
        // be a FULL drain rather than a remainder, for a lane that was away with a half-spent debt.
        // …and a lane that was NEVER fed owes nothing, which is not decoration: without `everFed_` a fresh
        // load arms lane 1 on a MONO host and runs a second network for a whole drain — 132 ms of a real
        // WaveNet, per load, for a window that is already the silence state NAM zero-filled it with.
        // 🔴 AND "NEVER FED" NOW HAS A SECOND WAY IN: a lane a `reset()` emptied is provably in that same
        // state, so the restart clears `everFed_` for it and this line charges it nothing — which is what
        // keeps a re-prepare after a restart free. A RECURRENT capture is the exception at both ends: the
        // restart leaves its flag set, so this line goes on charging it. See reset().
        for (int c = 0; c < 2; ++c) drain_[c] = everFed_[c] ? drainSamples_ : 0;
        hush_.assign ((size_t) maxBlock, 0.0f);            // ALLOC here (message thread) — never in process
        return true;
    }

    // Per-channel resampler state + model-rate scratch (used only when resampling).
    struct Ch
    {
        felitronics::core::StreamResampler down, up;       // host->model, model->host
        std::vector<float> modelIn, modelOut;
    };

    // RT-safe: run one channel through its model instance, with optional resampling, applying makeup.
    void processChannel (Ch& c, ::nam::DSP* m, float* io, int n, float g) noexcept
    {
        if (! resampling)
        {
            // copy → model scratch (NAM may not allow in==out) → process → back, with makeup.
            std::copy (io, io + n, c.modelIn.data());
            NAM_SAMPLE* in [1] = { c.modelIn.data() };
            NAM_SAMPLE* out[1] = { c.modelOut.data() };
            m->process (in, out, n);
            for (int i = 0; i < n; ++i) io[i] = c.modelOut[i] * g;
            return;
        }

        // host → model (variable count), run NAM, model → host (exactly n).
        c.down.feed (io, n);
        const int mFrames = c.down.produceAvailable (c.modelIn.data(), maxModelFrames);
        if (mFrames > 0)
        {
            NAM_SAMPLE* in [1] = { c.modelIn.data() };
            NAM_SAMPLE* out[1] = { c.modelOut.data() };
            m->process (in, out, mFrames);
            c.up.feed (c.modelOut.data(), mFrames);
        }
        c.up.produceExact (io, n);
        for (int i = 0; i < n; ++i) io[i] *= g;
    }

    std::unique_ptr<::nam::DSP> inst[2];
    const std::atomic<bool>*  normalize = nullptr;   // NamStage::Impl's per-call flag (bindNormalize)
    std::atomic<long long>*   drained_  = nullptr;   // …and its silence odometer (bindDrainCounter)
    std::atomic<long long>*   cleared_  = nullptr;   // …and the restart's own (bindClearCounter)
    int*  retireLedger = nullptr;                    // attached only once live (see attachRetireLedger)
    bool  prepared_    = false;                      // false until prepare() fully succeeded

    float  makeup      = 1.0f;
    double expectedSR  = 0.0;
    double loudnessDb  = 0.0;
    bool   hasLoudness = false;
    int    prewarm     = 0;

    double hostSR   = 48000.0;
    int    maxBlock = 512;
    bool   resampling = false;              // hostSR != model native rate
    double modelRunSR = kModelSampleRate;   // the rate the NAM instances are Reset to / run at
    NamStage::RateMatch rm_ { kModelSampleRate, false, 0 };   // decided once per prepare(), reported after
    int    maxModelFrames = 1024;
    int    drainSamples_  = 0;              // how long an absent lane is fed silence (configureRates)
    int    drainTail_     = 0;              // …plus what the engine holds past the field (Linear's FFT ring)
    bool   recurrent_     = false;          // …and whether the field is a bound at all (LSTM)
    int    drain_[2] { 0, 0 };              // …and how much of that each lane still owes
    bool   everFed_[2] { false, false };    // …and whether it may still be holding some (see reset())
    std::vector<float> hush_;               // the silence an absent lane is fed, and where its output goes
    Ch     ch[2];
};

static_assert (felitronics::neural::Inference<NamBackend>,
               "NamBackend must satisfy the felitronics::neural process-only inference seam");

} // namespace

//==============================================================================
struct NamStage::Impl
{
    // UI mirrors (message thread reads these) — published ONLY when a swap/clear actually lands
    // (swap first, mirrors second), so they always describe the model that is audibly live.
    std::atomic<int>    prewarmSamples { 0 };
    std::atomic<double> expectedSR  { 0.0 };
    std::atomic<double> loudnessDb  { 0.0 };
    std::atomic<bool>   hasLoudness { false };

    // Per-call `normalize` handoff to the live backend (see NamStage::process / NamBackend ctor).
    std::atomic<bool> normalize { true };

    // How many samples of silence the live backend has fed to absent lanes — see drainedSamples().
    std::atomic<long long> drained { 0 };

    // …and how many a reset() has spent restarting them — see clearedSamples(). Two mechanisms, two
    // odometers: past the debt the audio is the same either way, so a single number could not say
    // which of them moved.
    std::atomic<long long> cleared { 0 };

    // EXACT mirror of NeuralStage's internal retire-queue count (which it doesn't expose): +1 when
    // a successful swap/clear retires the live model; -1 from the dtor of every once-live backend
    // when the GC frees it (instances are attached to the ledger only when they go live, so a
    // superseded pending instance doesn't skew it; the live one freed at stage teardown decrements
    // a dying int — harmless). Needed because swapPrepared() takes ownership and destroys the
    // handed-in instance when it refuses: tryApplyPending() must KNOW a swap can't refuse before
    // handing over a load it is not allowed to lose. Message thread only.
    int retiredLedger = 0;

    // The swap-safe holder: atomic live-pointer swap, block-counter retire, message-thread GC.
    // Declared AFTER every member the backends touch — the normalize flag they point at and the
    // ledger their dtors decrement — because members destruct in reverse order: the stage, and
    // every backend it still owns, must die first.
    //
    // Refusal contract (verified against NeuralStage.h v0.8.0 — settled, don't re-litigate):
    // a swap/clear against a FULL retire queue returns false BEFORE the live-pointer exchange,
    // so the LIVE model is neither replaced nor deleted (no use-after-free, audio keeps playing
    // it). The only thing freed on refusal is the handed-in candidate — by its own unique_ptr,
    // on the message thread. That loss-on-refusal is exactly why the pending machinery below
    // never attempts a swap it cannot prove will succeed.
    felitronics::neural::NeuralStage<NamBackend, kMaxRetiredModels> stage;

    // Deferred intent, last-wins (message thread only): set when the retire queue was full at
    // swap/clear time. nullptr backend = a pending CLEAR (mirroring stage.clear() ==
    // swapPrepared (nullptr)). Retried by tryApplyPending() on every message-thread touch —
    // load / clear / prepare / the periodic collectGarbage() tick — until it lands.
    bool pendingActive = false;
    std::unique_ptr<NamBackend> pendingBackend;

    double hostSR   = 48000.0;
    int    maxBlock = 512;
    // There is deliberately no `modelRunSR` here any more. It existed to be the reference the
    // mid-stream rate contract compared against, and being a MUTABLE reference for a tolerance is
    // exactly what let the run rate walk: kModelSampleRate is the reference now, and it cannot move.

    void publishMirrors (double sr, double ldb, bool hl, int prewarm = 0)
    {
        prewarmSamples.store (prewarm, std::memory_order_relaxed);
        expectedSR.store (sr, std::memory_order_relaxed);
        loudnessDb.store (ldb, std::memory_order_relaxed);
        hasLoudness.store (hl, std::memory_order_relaxed);
    }

    // Message thread. Drain the retire queue, then try to land the deferred intent (if any).
    // A pending CLEAR is lossless to attempt (a refused clear() has no side effects — see the
    // refusal contract above), so it is simply retried. A pending SWAP is attempted ONLY when
    // NeuralStage's own precondition guarantees success: live == nullptr can never refuse;
    // otherwise the ledger must show room in the retire queue.
    //
    // Returns true when an intent LANDED in this call — the model state (and with it possibly
    // latencySamples()) just changed. Load/clear/prepare callers already re-report host PDC via
    // their normal paths and may ignore the result; collectGarbage() forwards it, so the
    // processor's timer re-reports PDC when a DEFERRED swap/clear lands between those calls.
    bool tryApplyPending()
    {
        stage.collectGarbage();            // frees retired models → their dtors drop the ledger
        if (! pendingActive)
            return false;

        if (pendingBackend == nullptr)     // pending CLEAR
        {
            const bool hadModel = stage.hasModel();
            if (! stage.clear())
                return false;              // queue still full — keep the intent for the next tick
            if (hadModel)
                ++retiredLedger;           // the cleared model just entered the retire queue
            pendingActive = false;
            publishMirrors (0.0, 0.0, false);
            return true;
        }

        // swapPrepared() refuses only when a model is live AND the retire queue is full; the
        // ledger mirrors that queue exactly, so this test equals the core's own precondition.
        if (stage.hasModel() && retiredLedger >= kMaxRetiredModels)
            return false;                  // would lose the load — keep it parked instead

        const bool   hadModel = stage.hasModel();
        const double sr  = pendingBackend->reportedSampleRate();
        const double ldb = pendingBackend->reportedLoudnessDb();
        const bool   hl  = pendingBackend->reportedHasLoudness();
        const int    pw  = pendingBackend->reportedPrewarm();
        pendingBackend->attachRetireLedger (&retiredLedger);
        const bool ok = stage.swapPrepared (std::move (pendingBackend));
        pendingActive = false;
        if (! ok)                          // unreachable given the precondition above; kept honest —
            return false;                  // mirrors stay untouched if it ever fired
        if (hadModel)
            ++retiredLedger;               // the replaced model just entered the retire queue
        publishMirrors (sr, ldb, hl, pw);
        return true;
    }
};

//==============================================================================
NamStage::NamStage()  : impl (std::make_unique<Impl>()) {}
NamStage::~NamStage() = default;   // NeuralStage frees the live + retired models

void NamStage::prepare (double sampleRate, int maxBlock)
{
    impl->hostSR   = sampleRate;
    impl->maxBlock = std::max (1, maxBlock);

    // Audio is stopped here. NeuralStage re-prepares the live backend in place (it captures the live
    // pointer ONCE, so a concurrent message-thread swap can't split the rate-config from the Reset).
    // What this used to do as well was ADOPT the live model's rate as the reference for the next
    // load's rate check — which made every accepted load move the goalposts. Nothing is adopted now.
    if (impl->pendingBackend != nullptr)                                // a parked load must follow the
        impl->pendingBackend->prepare (sampleRate, impl->maxBlock, 2);  // new rates before it goes live
    impl->stage.prepare ({ sampleRate, impl->maxBlock, 2 });
    impl->tryApplyPending();   // if the retire queue drained since the deferral, land the intent now
}

// THE STREAM RESTART — the contract and the price are on NamBackend::reset(); this is the wiring.
// It reaches only the LIVE backend, which is the only one that can be dirty: a parked load has been
// prepared and never fed, and NeuralStage resolves the live pointer once (acquire), so a swap landing
// beside this restarts one instance or the other and never half of each. It does NOT step the block
// counter — that counter is what keeps the garbage collector off an instance the audio thread may be
// inside, and a restart is not a block. It does not touch the pending intent either: a parked load is
// something the caller asked for, not something the previous stream left behind.
void NamStage::reset() noexcept { impl->stage.reset(); }

//==============================================================================
bool NamStage::process (float* const* io, int numChannels, int numSamples, bool normalize) noexcept
{
    // The shared Inference seam is process(io, nc, n) — the per-block `normalize` flag travels via
    // an atomic the live backend reads inside the SAME call (same thread → sequenced). Block
    // counting + the live-model resolve live in NeuralStage; no model → clean passthrough.
    if (io == nullptr && numChannels > 0 && numSamples > 0) return false;
    impl->normalize.store (normalize, std::memory_order_relaxed);
    return impl->stage.process (io, numChannels, numSamples);
}



//==============================================================================
// A model built and prepared with no stage in sight: everything a load costs, minus the swap.
struct NamStage::Prepared
{
    std::unique_ptr<NamBackend> backend;
};

void NamStage::PreparedDeleter::operator() (Prepared* p) const noexcept { delete p; }

NamStage::PreparedModel NamStage::prepareModel (const void* data, std::size_t size,
                                                double sampleRate, int maxBlock, float trimDb)
{
    // Accept BOTH raw .nam JSON and packed .namz (weights as float32). A packed blob is unpacked
    // to the equivalent JSON first, then the existing parse→get_dsp path runs unchanged — .namz
    // is bit-exact to the float32 the engine computes, so the model is identical. Any thread but
    // the audio thread: nothing here touches a stage, and the alloc/parse cost is the point.
    constexpr std::size_t kMaxUnpackedNamBytes = 64u * 1024u * 1024u;   // zip-bomb guard on unpack

    std::unique_ptr<::nam::DSP> m0, m1;
    int prewarmFromConfig   = 0;
    int drainTailFromConfig = 0;
    bool recurrentFromConfig = false;
    try
    {
        std::vector<std::uint8_t> unpacked;            // owns reconstructed JSON iff input was packed
        const char* begin = static_cast<const char*> (data);
        const char* end   = begin + size;
        if (namz::isNamz (data, size))
        {
            unpacked = namz::unpack (data, size, kMaxUnpackedNamBytes);
            if (unpacked.empty())
                return nullptr;
            begin = reinterpret_cast<const char*> (unpacked.data());
            end   = begin + unpacked.size();
        }
        auto j = nlohmann::json::parse (begin, end);
        m0 = ::nam::get_dsp (j);            // two independent instances of the same capture
        m1 = ::nam::get_dsp (j);
        prewarmFromConfig   = detail::receptiveFieldFromConfig (j);
        drainTailFromConfig = detail::partitionedTailSamples (j);
        recurrentFromConfig = detail::isRecurrent (j);
    }
    catch (...) { return nullptr; }

    if (m0 == nullptr || m1 == nullptr)
        return nullptr;
    // prototype: mono amp captures only
    if (m0->NumInputChannels() != 1 || m0->NumOutputChannels() != 1)
        return nullptr;
    const double modelSR = m0->GetExpectedSampleRate();
    // THE RATE CONTRACT, and it is settled HERE because it reads nothing but this number. It used to
    // live in install(), judged against the rate of whatever was live at the last prepare() — a
    // reference that MOVED as loads were accepted, which is the ratchet this task closed. Judging it
    // here spares a doomed model its PREWARM, which is the expensive half of everything below; it does
    // not spare the parse or the two get_dsp() calls above, because the tag is read off the built
    // network and there is nothing to judge until then.
    if (! acceptsModelRate (modelSR))
        return nullptr;

    // Build + prepare the new backend while it is NOT live (alloc + prewarm is fine here — and a
    // low-memory failure fails the LOAD, never the host: prepare() self-catches, see NamBackend).
    std::unique_ptr<Prepared> out;
    try
    {
        out = std::make_unique<Prepared>();
        out->backend = std::make_unique<NamBackend> (std::move (m0), std::move (m1), trimDb, prewarmFromConfig,
                                                     drainTailFromConfig, recurrentFromConfig);
    }
    catch (...) { return nullptr; }
    out->backend->prepare (sampleRate, std::max (1, maxBlock), 2);
    if (! out->backend->prepared())
        return nullptr;
    return PreparedModel (out.release());
}

bool NamStage::install (PreparedModel model)
{
    if (model == nullptr || model->backend == nullptr)
        return false;

    // The rate contract was settled in prepareModel() — acceptsModelRate() reads no stage state, so a
    // handle that exists is one this stage takes. Nothing about a rate is decided here any more.
    auto backend = std::move (model->backend);
    backend->bindNormalize (impl->normalize);
    backend->bindDrainCounter (impl->drained);
    backend->bindClearCounter (impl->cleared);
    // Prepared for other numbers than this stage runs at: the same work again, here — the rare case
    // of a rate change between the two halves, at the cost a one-call load always paid.
    //
    // 🔴 THE TEST IS EXACT, AND IT USED TO BE A HALF-HERTZ TOLERANCE. Two tolerances of the same size
    // do not compose: this one decided whether the backend's rate-match is still VALID, while
    // rateMatch() decided what that rate-match IS, and a host rate landing between them left a backend
    // answering for a rate it is not installed into. Measured on the base commit — backend prepared for
    // host 48000.4 installed into a stage at 48000.6: accepted, reports 0 samples where the policy
    // charges 64; the mirror (prepared 48000.6, stage 48000.4) reports 64 where the policy charges 0.
    // Equality is the only spelling under which "what the stage reports IS rateMatch(hostSR, tag)"
    // holds with no argument, and it costs a re-prepare only when a host hands two different doubles
    // for one rate — the cost the one-call load pays every time anyway.
    // Re-preparation also changes model state: NAM's LSTM Reset adds another prewarm to the existing
    // cell. Prepared at 48000.1 and installed at 48000.2 (block 512), the decaying-cell fixture starts
    // at 0.0548291542 instead of 0.1600718498, although both configurations report zero latency.
    // A matching prepareModel+install does one prewarm; a mismatch does two. See plan 1 in the review.
    if (backend->preparedSampleRate() != impl->hostSR || backend->preparedMaxBlock() != impl->maxBlock)
        backend->prepare (impl->hostSR, impl->maxBlock, 2);
    // 🔴 AND THE VERDICT IS CHECKED WHETHER OR NOT WE RE-PREPARED. It used to be checked only inside
    // that branch, which was safe by an invariant nobody stated: `prepareModel()` returns null for a
    // backend that failed to prepare, so an unprepared one could not reach here. `preparedSampleRate()`
    // reports the rate the backend was ASKED for, not one it honoured, so a handle prepared for exactly
    // this stage's numbers takes the skip path — and if that preparation had failed, the skip would
    // install a passthrough that reports a model and no latency. Nothing reaches it today; the check
    // costs a load of an already-loaded bool and stops depending on a guarantee two functions away.
    if (! backend->prepared())
        return false;

    // Accepted — from here the load is GUARANTEED to apply. Park it as the (single, last-wins)
    // pending intent and try to land it now: the normal path lands immediately; only a full
    // retire queue (audio frozen across kMaxRetiredModels swaps) defers it to a later drain
    // (collectGarbage / prepare / the next load or clear). Mirrors update only when it lands.
    impl->pendingBackend = std::move (backend);
    impl->pendingActive  = true;
    impl->tryApplyPending();
    return true;
}

bool NamStage::loadModelFromMemory (const void* data, std::size_t size, float trimDb)
{
    // Both halves, here, on the caller's thread — the message thread, by this class's contract.
    return install (prepareModel (data, size, impl->hostSR, impl->maxBlock, trimDb));
}

void NamStage::clearModel()
{
    // Last-wins intent: a clear supersedes any parked (never-applied) load. Mirrors are zeroed
    // only when the clear actually LANDS (swap first, mirrors second) — a deferred clear keeps
    // reporting the model that is still audibly live, instead of lying "empty".
    impl->pendingBackend.reset();
    impl->pendingActive = true;
    impl->tryApplyPending();
}

bool NamStage::collectGarbage()
{
    // Drains the retire queue, then lands any deferred swap/clear. True = a deferred intent
    // landed HERE (between load/clear calls) — the caller must re-report host PDC exactly like
    // after a normal load, because the landing can change latencySamples().
    return impl->tryApplyPending();
}

bool   NamStage::hasModel()         const { return impl->stage.hasModel(); }
double NamStage::modelSampleRate()  const { return impl->expectedSR.load  (std::memory_order_relaxed); }
double NamStage::modelLoudness()    const { return impl->loudnessDb.load  (std::memory_order_relaxed); }
bool   NamStage::modelHasLoudness() const { return impl->hasLoudness.load (std::memory_order_relaxed); }
// THE ONE ANSWER. Three facts, each written down exactly once, and the order between them is part of
// the contract rather than an accident of how the old code happened to be laid out.
//
//  1. NORMALISE. A model that reports no rate (<= 0) runs at the factory rate. This happens FIRST, so
//     the gate below compares against the rate the model will ACTUALLY run at. Reverse the two and an
//     untagged model at a 48 kHz host would be judged against 0 and come out "resampling".
//  2. GATE. A resampler is installed only past half a hertz of difference. Below that the rates are
//     the same clock as far as anything audible is concerned, and installing a 64-tap kernel to
//     convert 48000 to 48000.4 would cost 64 samples of delay to no purpose.
//  3. DERIVE, and only if one is installed. The geometry is core's — `pairDelayHostSamples` — and
//     the rounding is ours: the true delay is fractional and a host wants an integer, so round to
//     nearest. The residual is at most half a sample (worst on the shipped grid: 0.40 at 44.1 kHz,
//     whose first comb notch against an undelayed dry path sits at 55 kHz, out of band).
//
// 🔴 NOT GUARDED, DELIBERATELY — and what an absurd hostSR does has now been measured rather than
// reasoned about, because two earlier wordings of this paragraph were wrong about it. There are FOUR
// regimes against a 48 kHz model, not two (m = the model rate, the delay is 32 + 32·h/m):
//
//   h negative      → a perfectly finite, perfectly useless number (-48000 gives exactly 0).
//   h NaN           → fails the gate, returns 0 (the comparison itself raises FE_INVALID).
//   h > ~3.22e12    → lround is fine, but NARROWING ITS long TO int silently invents a plausible
//                     positive answer: h = 1e18 reports 1842981579 samples of latency, no flag raised.
//   h > ~1.38e22    → the value passes out of long's range too: lround saturates and FE_INVALID is
//                     raised. An INFINITE host lands in this same regime, which is why the earlier
//                     claim that "only an infinite host" gets this far was false: a finite 1e23
//                     reaches it identically, on every row measured.
//
// 🔴 AND THE ANSWER IN THE LAST TWO REGIMES IS PLATFORM-SPECIFIC, so no number is quoted for it here.
// Measured on four rows rather than reasoned about, because an earlier draft of this paragraph quoted
// "-1" and that is one toolchain's answer out of three:
//
//   arm64 macOS / x86-64 macOS (Apple libm) : lround saturates to LONG_MAX  → (int) = -1
//   x86-64 Debian (gcc 14 + glibc)          : lround saturates to LONG_MIN  → (int) =  0
//   x86-64 Windows (MSVC 19.44 + UCRT)      : long is 32 BITS, so lround saturates far earlier —
//                                             even the 1e18 case, where all three POSIX rows agree on
//                                             1842981579, reads 0 there.
//
// The two Mac rows and the Debian row share an ISA in one pairing and a toolchain in the other, which
// is what identifies libm rather than the ISA as the thing that differs. On the SHIPPED grid — 16 host
// rates x 3 model rates, up to the 3 MHz ceiling rigplayer now enforces — all four rows are
// byte-identical, so this divergence lives strictly outside the documented domain.
//
// None of that is new: the base commit's latencySamples() computed `lround(d + d*hostSR/modelRunSR)`
// with the same types and the same gate, so every one of the four regimes predates this extraction.
// Adding a guard would be a BEHAVIOUR change wearing a refactor's clothes, while this commit promises
// that no number moves. The one input that could divide by zero cannot reach the division: step 1
// turns a non-positive model rate into the factory rate. A real contract for absurd host rates is a
// separate question — and the consumer that actually sizes a buffer from this, rigplayer, no longer
// depends on the answer: it bounds its host rate before it asks.
NamStage::RateMatch NamStage::rateMatch (double hostSR, double modelSR) noexcept
{
    RateMatch r {};
    r.modelRunSR = (modelSR > 0.0 ? modelSR : kModelSampleRate);
    r.resampling = std::abs (hostSR - r.modelRunSR) > kModelRateTolerance;
    r.latencySamples = r.resampling
        ? (int) std::lround (felitronics::core::StreamResampler::pairDelayHostSamples (hostSR, r.modelRunSR))
        : 0;
    return r;
}

// The rate contract, spelled ONCE and through the owner of the tolerance rather than beside it: this
// stage takes a model exactly when a factory-rate host would run it without a resampler. The `<= 0`
// normalisation is rateMatch()'s too, so an untagged capture is accepted here for the same reason it
// runs at kModelSampleRate there — one sentence, one place, no second copy to fall out of step.
bool NamStage::acceptsModelRate (double modelSR) noexcept
{
    return ! rateMatch (kModelSampleRate, modelSR).resampling;
}

// The worst an accepted model can cost this host — the number a consumer sizes a fixed delay line
// from. `pairDelayHostSamples` is monotonically DECREASING in the model rate (the return leg is
// converted at hostSR/modelSR) and the accepted window's low edge is kModelSampleRate -
// kModelRateTolerance, so one evaluation there dominates every accepted model and no search is needed.
// Asking at the NOMINAL rate instead — which is what a consumer does when it derives its own number
// from kModelSampleRate — reads up to 0.0208 samples short at a 3 MHz host, and a shortfall of any
// size is a delay line that clamps in silence.
//
// 🔴 IT ASKS THE GEOMETRY, NOT rateMatch(), AND THAT IS DELIBERATE. rateMatch's answer is zero
// whenever the two rates are within half a hertz, and that gate depends on the MODEL rate as well as
// the host's — so the window's low edge is not always the worst cell of rateMatch: at a 47999.5 Hz host
// the low edge costs nothing while the window's HIGH edge resamples and costs 64. Taking the geometry
// at the low edge is >= every accepted model's reported latency at every host, gate or no gate, which
// is the property a ring needs. As a BOUND it holds everywhere; as an attained maximum it has two
// exceptions, and both are stated because "exact" was claimed here once with only the first in mind:
// a host at exactly kModelSampleRate, where no accepted model resamples at all and the true maximum is
// 0 (64 slots of a ring, paid so the rule has no exception in the code); and any host so slow that no
// backend can be PREPARED for it at all, where the stage stays at zero because it never runs. This is
// a policy answer about rates, not a promise about a particular prepared instance.
int NamStage::maxLatencySamples (double hostSR) noexcept
{
    return (int) std::lround (felitronics::core::StreamResampler::pairDelayHostSamples (
                                  hostSR, kModelSampleRate - kModelRateTolerance));
}

int    NamStage::latencySamples()   const { return impl->stage.latencySamples(); }
int    NamStage::prewarmSamples()   const { return impl->prewarmSamples.load (std::memory_order_relaxed); }
long long NamStage::drainedSamples() const { return impl->drained.load (std::memory_order_relaxed); }
long long NamStage::clearedSamples() const { return impl->cleared.load (std::memory_order_relaxed); }

} // namespace felitronics::nam

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
                float trimDb, int prewarmFromConfig = 0)
    {
        prewarm = prewarmFromConfig;          // …raised below by NAM's own answer, where it has one
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
            prepared_ = true;
        }
        catch (...) {}   // bad_alloc (or a throwing NAM Reset): stay unprepared, never crash
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
        if (numChannels == 0 || numSamples == 0) return true;

        const float g = (normalize != nullptr && normalize->load (std::memory_order_relaxed)) ? makeup : 1.0f;
        for (int off = 0; off < numSamples; )
        {
            const int n = std::min (numSamples - off, maxBlock);
            processChannel (ch[0], inst[0].get(), io[0] + off, n, g);          // mono track → 1 instance
            if (numChannels > 1)                                               // stereo → 2 independent instances
                processChannel (ch[1], inst[1].get(), io[1] + off, n, g);
            off += n;                                      // `off += maxBlock` could step past INT_MAX
        }
        return true;
    }

    void reset() noexcept {}   // transient state cleared by prepare()'s Reset on the next play

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

void NamStage::reset() {}   // transient state cleared by prepare()'s Reset on the next play

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
    int prewarmFromConfig = 0;
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
        prewarmFromConfig = detail::receptiveFieldFromConfig (j);
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
        out->backend = std::make_unique<NamBackend> (std::move (m0), std::move (m1), trimDb, prewarmFromConfig);
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

} // namespace felitronics::nam

// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/core/Fft.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cstring>
#include <vector>

namespace felitronics::convolution
{

//==============================================================================
// felitronics::convolution::MatrixConvolverNupc — the NON-UNIFORM (Gardner) sibling of MatrixConvolver:
// the SAME public 2×2-operator API on one canonical raw-L/R history, but backed by the block-INDEPENDENT
// non-uniform partition schedule (a time-domain head P0 + growing overlap-save stages) proven in
// NonUniformConvolver. Drop-in for `MatrixConvolver<Fft>`: identical prepare/setOperator/stageOperator/
// publishStaged/setIr/process/reset/isBusy/latencySamples()==0 signatures, so lineareq's
// `using Conv = MatrixConvolver<AudioFft>` becomes `= MatrixConvolverNupc<AudioFft>` with no other change.
//
// COMPLETE. ALL topologies (mono / LRDiag / MSDiag / Full) + the CLICK-FREE 2-slot smoothstep WARM crossfade
// swap — a full drop-in for MatrixConvolver.
//   • MSDiag forms the ½(X_L ± X_R) spectral views PER STAGE from the per-stage L/R FDLs (macViewStage), MACs
//     bank_M / bank_S, and decodes yL = yM+yS, yR = yM−yS; Full does the 4-bank cross sums per stage.
//   • setOperator stages an operator into the inactive slot; process() crossfades the OLD and NEW operators'
//     FULL outputs (head + every stage's tail) by a smoothstep weight, BOTH computed from the SAME shared warm
//     FDL (chunkStage recomputes both slots' tails as stages wrap). Every swap uses the same short smoothstep
//     crossfade (crossfadeSamples): a cold-started FDL already produces the EXACT causal convolution, so the fade
//     only masks the silence→convolution onset (no long "cold prime"). isBusy() is true while fading (host coalesces).
//
// STRUCTURE (channel-indexed per-stage, mirrors MatrixConvolver's decoupled helpers replicated per stage):
//   • the IR is tiled by a head [0,P0) + stages s of block B_s, count C_s, offset_s==B_s (zero-latency
//     invariant, recurrence B_{s+1}=(C_s+1)·B_s); the default schedule caps doubling at B_max=2048.
//   • per stage: an audio FFT (size 2·B_s), a message-thread build FFT, and a PER-CHANNEL input frame + FDL
//     (the raw-L and raw-R past input spectra at that stage's size). LRDiag reads L-FDL for yL, R-FDL for yR;
//     MSDiag/Full (later) read BOTH per stage. IR-independent, so a topology change is just a bank swap.
//   • per operator slot: per-bank {head taps h0, per-stage IR partition spectra} + per-stage per-output-channel
//     cached tails. Two slots for the click-free swap (Phase 4).
//   • ONE shared scratch set sized to the largest stage — stages fire SEQUENTIALLY in the per-sample advance
//     loop (no concurrency), so a max-sized acc/ifft/input scratch is reused across them.
//
// COST NOTE. Two parts, and only one of them follows maxIrSamples.
//   • FIXED BY THE SCHEDULE, regardless of the live IR: every stage takes its input sample, the head dot runs
//     over P0 taps, and each stage FFTs its input and inverse-FFTs each output on its own boundary. This keeps
//     every stage's FDL warm so a swap to a LONGER IR finds real history (the MatrixConvolver invariant).
//     maxIrSamples sets how many STAGES there are, and the count stops growing at the B_max stage: with a
//     128-sample head, five stages for any cap above 2048 samples, whatever the cap. A larger cap only adds
//     partitions to the last stage.
//   • FOLLOWS THE LOADED IR: the spectral MACs, the bulk of the work for a long IR, loop over each bank's
//     `activeParts` — the partitions the current IR reaches — never over the stage's full count (macBankStage /
//     macViewStage). A short IR in a long cap pays for its own length plus the fixed part above.
// A future "optimize short IRs" pass must NOT simply gate stages off — that would cold-start the FDL and break
// warm swaps.
//
// MEMORY NOTE. The same split. The shared FDL is sized to maxIrSamples at prepare(), because it is ONE history
// that both slots read: a new operator answers the last maxIrSamples of input at once, which is what keeps a
// swap click-free and bit-exact. The bank spectra follow the LOADED IR: stageOperator() grows each bank, per
// stage, to the partitions the IR reaches, and allocates only the banks the topology reads (mono 1, LRDiag and
// MSDiag 2, Full 4). They GROW ONLY: a shorter IR reuses the storage a longer one left, so each slot holds its
// own high-water mark until the next prepare(), which releases it. A caller that must not allocate after its
// prepare() reserves up front instead (reserve()). With the default schedule a sample of cap costs ~8 B of FDL
// per channel (~16 B stereo); a sample of loaded IR costs ~8 B per bank, in each slot that holds it.
//
// RT-safe: prepare()/setOperator() allocate (message thread); process() is alloc/lock/throw-free, zero latency.
template <core::fft::RealFftBackend Fft = core::fft::DefaultRealFft>
class MatrixConvolverNupc
{
public:
    enum class Topology { MSDiag, LRDiag, Full };

    static constexpr int kMaxBanks  = 4;
    static constexpr int kMaxStages = 16;
    static constexpr int kMaxPartition = 1 << 24;   // hard cap so the FFT size 2*B never overflows signed int (far above any real IR)
    static constexpr int kMaxIrSamples = 1 << 24;   // hard cap so a stage offset (stored as int) can't overflow (16M smp ≫ any IR)
    static constexpr int kDefaultMaxBlock = 2048;   // B_max cap: flat mean at any B_max, smaller worst-buffer spike than 4096

    static int numBanksFor (Topology t) noexcept { return t == Topology::Full ? 4 : 2; }

    // Signature-compatible with MatrixConvolver::prepare. `partitionSize` is the time-domain HEAD size P0
    // (lineareq passes 128); B_max is an internal default. crossfadeSamples is the warm anti-click fade length —
    // use >= ~32 for a click-free ramp (below that the smoothstep degenerates toward a 1-sample hard switch; the
    // floor is 1, matching MatrixConvolver). lineareq passes ~20 ms, well above that.
    bool prepare (int partitionSize, int maxIrSamples, int crossfadeSamples, int numChannels = 2)
    {
        prepared_ = false;
        if (numChannels < 1 || numChannels > 2) return false;
        if (! core::fft::isPow2 (partitionSize) || partitionSize > kMaxPartition) return false;
        if (maxIrSamples > kMaxIrSamples) return false;
        if (maxIrSamples < 0) maxIrSamples = 0;
        channels_ = numChannels;
        mono_     = (channels_ == 1);
        P0_ = partitionSize;
        headMask_ = P0_ - 1;

        const int bMax = std::max (P0_, kDefaultMaxBlock);
        numStages_ = buildSchedule (P0_, bMax, maxIrSamples);
        if (numStages_ < 0) return false;                       // schedule can't cover maxIrSamples in kMaxStages

        warmXfade_ = crossfadeSamples < 1 ? 1 : crossfadeSamples;

        int maxSpecF = 1, maxN = 1, maxB = 1;
        for (int st = 0; st < numStages_; ++st) { maxSpecF = std::max (maxSpecF, spec_[st].specF); maxN = std::max (maxN, spec_[st].N); maxB = std::max (maxB, spec_[st].B); }
        inputSpec_.assign ((std::size_t) maxSpecF, 0.0f);
        viewSpec_.assign  ((std::size_t) maxSpecF, 0.0f);
        acc_.assign       ((std::size_t) maxSpecF, 0.0f);
        ifftOut_.assign   ((std::size_t) maxN,     0.0f);
        buildBuf_.assign  ((std::size_t) maxN,     0.0f);
        tmpTailA_.assign  ((std::size_t) maxB,     0.0f);
        tmpTailB_.assign  ((std::size_t) maxB,     0.0f);

        for (int ch = 0; ch < channels_; ++ch) { headHist_[ch].assign ((std::size_t) P0_, 0.0f); headPos_[ch] = 0; }

        for (int st = 0; st < numStages_; ++st)
        {
            History& h = hist_[(std::size_t) st];
            const int N = spec_[st].N, specF = spec_[st].specF, C = spec_[st].C;
            if (! h.fft.prepare (N)) return false;
            if (! h.buildFft.prepare (N)) return false;         // separate message-thread FFT (Phase-4 build/audio race-free)
            for (int ch = 0; ch < channels_; ++ch)
            {
                h.frame[ch].assign ((std::size_t) N, 0.0f);
                // Released first, so a re-prepare to a SMALLER cap (a lower host rate) gives the memory back;
                // assign() alone keeps the old capacity.
                core::fft::AlignedVector<float>().swap (h.fdl[ch]);
                h.fdl[ch].assign   ((std::size_t) C * (std::size_t) specF, 0.0f);
            }
            h.phase = 0; h.fdlPos = 0;
        }
        // What the new geometry no longer reads is released too: a stage the new schedule dropped, and a channel
        // a mono re-prepare dropped (a stereo history at the old cap would otherwise stay held beside it).
        for (int st = 0; st < kMaxStages; ++st)
            for (int ch = (st < numStages_ ? channels_ : 0); ch < 2; ++ch)
            {
                core::fft::AlignedVector<float>().swap (hist_[(std::size_t) st].fdl[ch]);
                core::fft::AlignedVector<float>().swap (hist_[(std::size_t) st].frame[ch]);
            }
        for (int ch = channels_; ch < 2; ++ch) std::vector<float>().swap (headHist_[ch]);

        // The head taps of every bank a slot can hold (a few hundred bytes each): a slot nothing was staged
        // into still runs its head dot. The bank SPECTRA start empty and are grown by stageOperator() to the IR
        // it is handed — see the MEMORY NOTE above.
        const int nb = mono_ ? 1 : kMaxBanks;
        for (int k = 0; k < 2; ++k)
        {
            Slot& sl = slot_[k];
            sl.topo = Topology::LRDiag; sl.numBanks = 0;
            for (int b = 0; b < kMaxBanks; ++b)
            {
                if (b < nb) sl.banks[b].h0.assign ((std::size_t) P0_, 0.0f);
                else        std::vector<float>().swap (sl.banks[b].h0);
                for (int st = 0; st < kMaxStages; ++st)
                {
                    core::fft::AlignedVector<float>().swap (sl.banks[b].irSpec[(std::size_t) st]);   // released
                    sl.banks[b].activeParts[st] = 0;
                }
            }
            for (int st = 0; st < numStages_; ++st)
                for (int ch = 0; ch < 2; ++ch) sl.tail[(std::size_t) st][(std::size_t) ch].assign ((std::size_t) spec_[st].B, 0.0f);
        }

        cur_ = 0; xfadePos_ = 0; xfadeLen_ = warmXfade_;
        state_.store (0, std::memory_order_relaxed);
        prepared_ = true;
        return true;
    }

    // 🔴 THE AUDIO THE CALLER FED, AND NOTHING ELSE — every buffer below is written ONLY inside
    // process(), so this touches no field the message thread owns and cancels no swap. That is the
    // whole difference from reset(), and it is what lets a composite restart its history on the audio
    // thread WITHOUT writing the two fields setOperator() reads (`cur_` and `state_`), and with a fade
    // left running rather than ended.
    //
    // Added for a composite's stream restart, against a measured sequence: publish a tone curve and
    // restart before the 50 ms crossfade ends. reset() below used to drop the incoming operator on the
    // floor there — `cur_` stayed on the OLD one and `CabConvolver::pendingRetry_` is already false after
    // a successful publish, so nothing ever re-staged it, and the knob move was lost until the next knob
    // move. P88 closed that at the source: reset() now ADOPTS the publication. This verb remains the half
    // that decides nothing and leaves a fade running — which is also why it does NOT give law 11a
    // independence mid-fade (CabConvolver::clearAudioState(), which forwards here, has the number); a
    // restart that owes independence calls reset().
    void clearAudioState() noexcept
    {
        for (int ch = 0; ch < channels_; ++ch) { std::fill (headHist_[ch].begin(), headHist_[ch].end(), 0.0f); headPos_[ch] = 0; }
        for (int st = 0; st < numStages_; ++st)
        {
            History& h = hist_[(std::size_t) st];
            for (int ch = 0; ch < channels_; ++ch) { std::fill (h.frame[ch].begin(), h.frame[ch].end(), 0.0f); std::fill (h.fdl[ch].begin(), h.fdl[ch].end(), 0.0f); }
            h.phase = 0; h.fdlPos = 0;
        }
        // BOTH slots' tails, including the incoming operator's: a fade that survives this call keeps
        // computing that operator's tail from `h.fdl`, which is now zeros, so leaving its tail behind
        // would feed the new operator the old stream through the back door.
        for (int k = 0; k < 2; ++k)
            for (int st = 0; st < numStages_; ++st)
                for (int ch = 0; ch < 2; ++ch) std::fill (slot_[k].tail[(std::size_t) st][(std::size_t) ch].begin(), slot_[k].tail[(std::size_t) st][(std::size_t) ch].end(), 0.0f);
    }

    // Flush the running history and KEEP THE OPERATOR THE CALLER LAST PUBLISHED — including one still
    // crossfading in, which this verb used to drop (see clearAudioState() above for the sequence that made
    // the loss visible, and ConvolutionEngine::reset() for the full reasoning: an accepted publication
    // is never lost, law 11a independence, and the price at the seam). A swap in flight is ENDED IN FAVOUR
    // OF THE NEW operator; one merely STAGED (stageOperator() without publishStaged(), state_ == 0) is left
    // exactly where it is, for the later publish to find. On the body this replaces, the suite's
    // independence group reads the OLD operator in every cell, a worst sample 1.927e+00 from the settled restart.
    // Audio thread (or externally synced); must not race setOperator().
    void reset() noexcept
    {
        clearAudioState();
        const int s = state_.load (std::memory_order_acquire);   // pairs with publishStaged()'s release store
        if (s == 0) return;                                      // nothing in flight, and nothing is stored:
                                                                 // an unconditional store(0) would wipe a
                                                                 // publication landing between load and store
        cur_ = 1 - cur_;                                         // the published slot becomes the live one
        xfadePos_ = 0;
        state_.store (0, std::memory_order_release);             // release: setOperator() reads `cur_` after
                                                                 // acquiring this store
    }

    static constexpr int latencySamples() noexcept { return 0; }
    bool isBusy() const noexcept { return state_.load (std::memory_order_acquire) != 0; }
    int  numChannels() const noexcept { return channels_; }

    // Message thread — single producer, must not race reset(). Stage a whole operator into the inactive slot,
    // then publish; the next process() picks it up and crossfades it in over `crossfadeSamples` (see the
    // state 1 -> 2 handoff in process() below — the warm swap is implemented, not pending).
    bool setOperator (Topology topo, const float* const* banks, int numBanks, int len)
    {
        if (! stageOperator (topo, banks, numBanks, len)) return false;
        publishStaged();
        return true;
    }

    // Build the operator into the inactive slot WITHOUT publishing (the surround ST-only path stages every
    // instance then publishes back-to-back). All topologies (mono / LRDiag / MSDiag / Full) are supported.
    //
    // ALLOCATES when the IR reaches further than this slot's banks have held before (the MEMORY NOTE above),
    // and only then. That is safe on THIS slot alone: with state_ == 0 — checked first — the audio thread reads
    // nothing of the inactive slot. Every read of `slot_[1 - cur_]` in process() is under state 1 (priming the
    // tails) or state 2 (the fade), both of which refuse here, and the fade's end stores `cur_` before the
    // release store of 0 that this acquire pairs with. Exhaustion throws out of here (law 11d) with the slot
    // still holding what it held: every bank grows before any is rebuilt or relabelled.
    bool stageOperator (Topology topo, const float* const* banks, int numBanks, int len)
    {
        if (! prepared_ || banks == nullptr || state_.load (std::memory_order_acquire) != 0) return false;
        const Topology t = mono_ ? Topology::LRDiag : topo;
        const int nb = mono_ ? 1 : numBanksFor (t);
        if (numBanks < nb) return false;
        // A null bank pointer is only rejected when it would be dereferenced (len > 0); (nullptr, len<=0) is the
        // documented "clear to silence" case buildBankSlice handles without a deref — matching MatrixConvolver.
        for (int b = 0; b < nb; ++b) if (banks[b] == nullptr && len > 0) return false;
        const int stg = 1 - cur_;
        Slot& sl = slot_[stg];
        for (int b = 0; b < nb; ++b) growBank (sl.banks[b], len);
        sl.topo = t; sl.numBanks = nb;
        for (int b = 0; b < nb; ++b) buildBankSlice (sl.banks[b], banks[b], len);
        return true;
    }

    void publishStaged() noexcept { state_.store (1, std::memory_order_release); }

    // Grow BOTH slots' banks for `topo` (prepared mono: the one bank) to an IR of `len` taps NOW, so that no later
    // swap up to that length and topology allocates. For a caller whose own contract is that its prepare()
    // allocates everything — lineareq's EQs, whose IR is always the cap's length and which promise a setBands()
    // that allocates nothing. The same thread contract as prepare(): no process() may run concurrently, because
    // this grows the LIVE slot too. Grow-only, like stageOperator(); exhaustion throws (law 11d). False when
    // unprepared.
    bool reserve (Topology topo, int len)
    {
        if (! prepared_) return false;
        const int nb = mono_ ? 1 : numBanksFor (topo);
        for (int k = 0; k < 2; ++k)
            for (int b = 0; b < nb; ++b) growBank (slot_[k].banks[b], len);
        return true;
    }

    // Mono convenience: broadcast a single IR onto both output channels (or the one bank when prepared mono).
    bool setIr (const float* ir, int len)
    {
        const float* both[2] { ir, ir };
        return setOperator (Topology::LRDiag, both, mono_ ? 1 : 2, len);
    }

    // Audio thread. Planar in/out may alias (in-place). RT-safe, zero latency. Needs EXACTLY channels_
    // planes — law 11(b)+(c); any other width is refused and RETURNS false.
    [[nodiscard]] bool process (const float* const* in, float* const* out, int numChannelsToProcess, int n) noexcept
    {
        if (numChannelsToProcess < 0 || n < 0) return false;
        if (! prepared_) return false;
        // The WIDTH is validated before the empty-length shortcut: `process(io, channels_ + 1, 0)` is a
        // malformed call that happens to have nothing to do, and answering "accepted" about it would tell
        // the caller its geometry is fine right up until the first non-empty block. Law 11's check order.
        // LAW 11(b)+(c): the matrix width is EXACT. A NARROW call is refused because a 2x2 matrix needs
        // BOTH input planes to compute either output — refused OBSERVABLY now, where it used to return
        // void, so a caller that pre-zeroed `out` got DIGITAL SILENCE and no way to find out (P18 F35,
        // measured through CabConvolver as a 9.12e-02 "leak" that was really forty unprocessed blocks
        // and a tail). A WIDER call is refused for the ordinary reason: the planes past `channels_` have
        // no operator to run and would have gone out DRY while the first two were convolved.
        if (numChannelsToProcess != channels_) return false;
        if (n == 0) return true;                                       // no samples: no time, no edge

        int s = state_.load (std::memory_order_acquire);
        if (s == 1)   // begin the crossfade. A cold-started FDL already yields the EXACT causal convolution, so ONE
        {             // short smoothstep fade (warmXfade_) masks only the silence→convolution onset — no long cold prime.
            xfadePos_ = 0;
            xfadeLen_ = warmXfade_;
            for (int st = 0; st < numStages_; ++st)         // prime the new slot's tails from the per-stage FDLs
            {
                int base = hist_[(std::size_t) st].fdlPos - 1; if (base < 0) base += spec_[st].C;
                computeStageTails (slot_[1 - cur_], st, base);
            }
            state_.store (2, std::memory_order_relaxed);
            s = 2;
        }

        if (s != 2) { processRange (in, out, 0, n); return true; }   // Idle → single live operator
        processFade (in, out, n);
        return true;
    }

private:
    //==========================================================================
    struct StageSpec { int B = 0, C = 0, offset = 0, specF = 0, N = 0; };

    // Per-stage input history — channel-indexed raw L/R past-input spectra at this stage's block size.
    struct History
    {
        Fft fft;        // audio-thread FFT (size 2B)
        Fft buildFft;   // message-thread build FFT (separate — Phase-4 build/audio race-free)
        core::fft::AlignedVector<float> frame[2];   // 2B per channel: [prev | cur]
        core::fft::AlignedVector<float> fdl[2];     // C×specF per channel: ring of past raw-input spectra
        int phase = 0, fdlPos = 0;
    };

    // An IR bank: head taps + per-stage tail partition spectra (built from a time-domain IR).
    struct Bank
    {
        std::vector<float> h0;                                                   // P0 head taps (time domain)
        std::array<core::fft::AlignedVector<float>, kMaxStages> irSpec;          // per stage: C_s×specF_s
        int activeParts[kMaxStages] {};
    };

    // An operator slot: topology + banks + the cached per-stage per-output-channel tails.
    struct Slot
    {
        Topology topo = Topology::LRDiag;
        int numBanks = 0;
        Bank banks[kMaxBanks];
        std::array<std::array<std::vector<float>, 2>, kMaxStages> tail;   // tail[stage][outCh] — B_s each
    };

    // Capped-doubling schedule: head P0 + doubling stages (C=1) to B_max, then one uniform B_max stage covering
    // the rest. Fills spec_[], returns the stage count, or -1 if it can't cover maxIr within kMaxStages.
    int buildSchedule (int P0, int maxBlock, int maxIr) noexcept
    {
        int n = 0;
        long long cum = P0;
        int B = P0;
        while (cum < (long long) maxIr && n < kMaxStages)
        {
            int C;
            if (B < maxBlock) C = 1;
            else { long long rem = (long long) maxIr - cum; C = (int) ((rem + B - 1) / B); if (C < 1) C = 1; }
            spec_[n] = { B, C, (int) cum, Fft::spectrumFloats (2 * B), 2 * B };   // offset == cum == B (invariant)
            cum += (long long) C * B;
            ++n;
            if (B >= maxBlock) break;                                            // tail stage absorbed the remainder
            B <<= 1;
        }
        return (cum < (long long) maxIr) ? -1 : n;
    }

    // How many of stage `st`'s partitions a `len`-tap IR reaches: those whose first tap lies inside the IR,
    // never more than the stage holds. So an IR longer than the schedule convolves only what the schedule covers:
    // maxIrSamples rounded up to whole partitions of the last stage, exactly as before.
    int partsFor (int st, int len) const noexcept
    {
        const long long rem = (long long) len - spec_[st].offset;
        if (rem <= 0) return 0;
        return (int) std::min<long long> ((rem + spec_[st].B - 1) / spec_[st].B, spec_[st].C);
    }

    // Grow one bank's spectra so a `len`-tap IR fits, stage by stage. Never shrinks. Message thread, the inactive
    // slot only (stageOperator). reserve() first: it allocates exactly what is asked, where a growing resize()
    // alone may take up to twice that.
    void growBank (Bank& bank, int len)
    {
        for (int st = 0; st < numStages_; ++st)
        {
            const std::size_t need = (std::size_t) partsFor (st, len) * (std::size_t) spec_[st].specF;
            core::fft::AlignedVector<float>& v = bank.irSpec[(std::size_t) st];
            if (v.size() >= need) continue;
            v.reserve (need);
            v.resize (need, 0.0f);
        }
    }

    // Build one bank's head + per-stage partition spectra from a time-domain IR slice. Message thread; the
    // bank was grown for `len` first (growBank), so every partition written here has its storage.
    void buildBankSlice (Bank& bank, const float* ir, int len) noexcept
    {
        if (len < 0) len = 0;
        for (int i = 0; i < P0_; ++i) bank.h0[(std::size_t) i] = (i < len) ? ir[i] : 0.0f;
        for (int st = 0; st < numStages_; ++st)
        {
            const int B = spec_[st].B, offset = spec_[st].offset, specF = spec_[st].specF, N = spec_[st].N;
            const int parts = partsFor (st, len);
            for (int j = 0; j < parts; ++j)
            {
                const int base = offset + j * B;
                std::fill (buildBuf_.begin(), buildBuf_.begin() + N, 0.0f);
                for (int i = 0; i < B; ++i) { const int src = base + i; buildBuf_[(std::size_t) i] = (src < len) ? ir[src] : 0.0f; }
                hist_[(std::size_t) st].buildFft.forward (buildBuf_.data(), &bank.irSpec[(std::size_t) st][(std::size_t) j * (std::size_t) specF]);
            }
            bank.activeParts[st] = parts;
        }
    }

    float headDot (const Bank& bank, int ch) const noexcept
    {
        float h = 0.0f;
        const int pos = headPos_[ch], mask = headMask_;
        for (int i = 0; i < P0_; ++i) h += bank.h0[(std::size_t) i] * headHist_[ch][(std::size_t) ((pos - i) & mask)];
        return h;
    }

    void clearAcc (int st) noexcept { std::fill (acc_.begin(), acc_.begin() + spec_[st].specF, 0.0f); }

    // acc_ += Σ_j srcFdl[(base-j) wrap] .* bank.irSpec[st][j]  — one channel's stage history against one bank.
    void macBankStage (const core::fft::AlignedVector<float>& srcFdl, const Bank& bank, int st, int base) noexcept
    {
        const int C = spec_[st].C, specF = spec_[st].specF;
        for (int j = 0; j < bank.activeParts[st]; ++j)
        {
            int idx = base - j; if (idx < 0) idx += C;
            hist_[(std::size_t) st].fft.spectralMultiplyAdd (&srcFdl[(std::size_t) idx * (std::size_t) specF],
                                                             &bank.irSpec[(std::size_t) st][(std::size_t) j * (std::size_t) specF], acc_.data());
        }
    }

    // acc_ += Σ_j (½(X_L ± X_R))[base-j] .* bank.irSpec[st][j] — the on-the-fly per-stage M/S spectral view
    // (sign +→M, −→S), formed elementwise over the stage's spectrum floats by FFT layout-linearity.
    void macViewStage (const Bank& bank, int st, int base, float sign) noexcept
    {
        const int C = spec_[st].C, specF = spec_[st].specF;
        const core::fft::AlignedVector<float>& fL = hist_[(std::size_t) st].fdl[0];
        const core::fft::AlignedVector<float>& fR = hist_[(std::size_t) st].fdl[1];
        for (int j = 0; j < bank.activeParts[st]; ++j)
        {
            int idx = base - j; if (idx < 0) idx += C;
            const float* xl = &fL[(std::size_t) idx * (std::size_t) specF];
            const float* xr = &fR[(std::size_t) idx * (std::size_t) specF];
            for (int i = 0; i < specF; ++i) viewSpec_[(std::size_t) i] = 0.5f * (xl[i] + sign * xr[i]);
            hist_[(std::size_t) st].fft.spectralMultiplyAdd (viewSpec_.data(),
                                                             &bank.irSpec[(std::size_t) st][(std::size_t) j * (std::size_t) specF], acc_.data());
        }
    }

    void invAccToStage (std::vector<float>& tail, int st) noexcept
    {
        hist_[(std::size_t) st].fft.inverse (acc_.data(), ifftOut_.data());
        const int B = spec_[st].B;
        for (int i = 0; i < B; ++i) tail[(std::size_t) i] = ifftOut_[(std::size_t) (B + i)];   // overlap-save: last B
    }

    // Recompute one stage's cached tails for a slot from the FDL history ending at partition `base`.
    void computeStageTails (Slot& sl, int st, int base) noexcept
    {
        if (mono_)
        {
            clearAcc (st); macBankStage (hist_[(std::size_t) st].fdl[0], sl.banks[0], st, base); invAccToStage (sl.tail[(std::size_t) st][0], st);
            return;
        }
        switch (sl.topo)
        {
            case Topology::LRDiag:
                clearAcc (st); macBankStage (hist_[(std::size_t) st].fdl[0], sl.banks[0], st, base); invAccToStage (sl.tail[(std::size_t) st][0], st);
                clearAcc (st); macBankStage (hist_[(std::size_t) st].fdl[1], sl.banks[1], st, base); invAccToStage (sl.tail[(std::size_t) st][1], st);
                break;

            // NB (numerical): decoding M/S → L/R in float means the Side channel of a near-mono input
            // (x_L ≈ x_R ⇒ S ≈ 0) has a large RELATIVE error (~eps/‖S‖·‖M‖) — inherent to ANY L/R-output M/S
            // process (a time-domain M/S convolver shows the same), not the spectral view. The L/R OUTPUTS stay
            // accurate to float noise (~2e-7 rel) and the residual floor is the Mid level (~−136 dB), inaudible.
            case Topology::MSDiag:
                clearAcc (st); macViewStage (sl.banks[0], st, base, +1.0f); invAccToStage (tmpTailA_, st);   // yM tail
                clearAcc (st); macViewStage (sl.banks[1], st, base, -1.0f); invAccToStage (tmpTailB_, st);   // yS tail
                for (int i = 0; i < spec_[st].B; ++i)
                {
                    const float m = tmpTailA_[(std::size_t) i], sd = tmpTailB_[(std::size_t) i];   // m = yM, sd = yS
                    sl.tail[(std::size_t) st][0][(std::size_t) i] = m + sd;   // decode yL = yM + yS
                    sl.tail[(std::size_t) st][1][(std::size_t) i] = m - sd;   //        yR = yM − yS
                }
                break;

            case Topology::Full:
                clearAcc (st); macBankStage (hist_[(std::size_t) st].fdl[0], sl.banks[0], st, base); macBankStage (hist_[(std::size_t) st].fdl[1], sl.banks[1], st, base); invAccToStage (sl.tail[(std::size_t) st][0], st);   // yL = LL∗xL + LR∗xR
                clearAcc (st); macBankStage (hist_[(std::size_t) st].fdl[0], sl.banks[2], st, base); macBankStage (hist_[(std::size_t) st].fdl[1], sl.banks[3], st, base); invAccToStage (sl.tail[(std::size_t) st][1], st);   // yR = RL∗xL + RR∗xR
                break;
        }
    }

    // One stage's chunk boundary: FFT each channel's frame into the FDL, recompute the live operator's tails,
    // shift cur→prev, advance the ring.
    void chunkStage (int st, bool fading) noexcept
    {
        History& h = hist_[(std::size_t) st];
        const int B = spec_[st].B, specF = spec_[st].specF;
        for (int ch = 0; ch < channels_; ++ch)
        {
            h.fft.forward (h.frame[ch].data(), inputSpec_.data());
            std::memcpy (&h.fdl[ch][(std::size_t) h.fdlPos * (std::size_t) specF], inputSpec_.data(), (std::size_t) specF * sizeof (float));
        }
        computeStageTails (slot_[cur_], st, h.fdlPos);
        if (fading) computeStageTails (slot_[1 - cur_], st, h.fdlPos);   // keep the incoming operator's tail coherent through the fade
        for (int ch = 0; ch < channels_; ++ch)
            for (int i = 0; i < B; ++i) h.frame[ch][(std::size_t) i] = h.frame[ch][(std::size_t) (B + i)];
        if (++h.fdlPos >= spec_[st].C) h.fdlPos = 0;
    }

    // The FULL output (head + every stage's tail at its phase) for ONE slot, per the topology. Called once
    // (idle) or twice (during the crossfade — old + new slot from the SAME shared FDL) per sample.
    void computeSlotOutputs (const Slot& sl, float& oL, float& oR) const noexcept
    {
        if (mono_)
        {
            float y = headDot (sl.banks[0], 0);
            for (int st = 0; st < numStages_; ++st) y += sl.tail[(std::size_t) st][0][(std::size_t) hist_[(std::size_t) st].phase];
            oL = y; oR = 0.0f;
            return;
        }
        float yl = 0.0f, yr = 0.0f;
        switch (sl.topo)
        {
            case Topology::LRDiag:
                yl = headDot (sl.banks[0], 0); yr = headDot (sl.banks[1], 1);
                break;
            case Topology::MSDiag:
            {
                const std::vector<float>& hM = sl.banks[0].h0;
                const std::vector<float>& hS = sl.banks[1].h0;
                const int pos0 = headPos_[0], pos1 = headPos_[1], mask = headMask_;
                float headM = 0.0f, headS = 0.0f;                          // form m/s on the fly by linearity
                for (int i = 0; i < P0_; ++i)
                {
                    const float l = headHist_[0][(std::size_t) ((pos0 - i) & mask)];
                    const float r = headHist_[1][(std::size_t) ((pos1 - i) & mask)];
                    headM += hM[(std::size_t) i] * 0.5f * (l + r);
                    headS += hS[(std::size_t) i] * 0.5f * (l - r);
                }
                yl = headM + headS; yr = headM - headS;                    // tails already decoded to yL/yR
                break;
            }
            case Topology::Full:
            {
                const std::vector<float>& LL = sl.banks[0].h0, & LR = sl.banks[1].h0, & RL = sl.banks[2].h0, & RR = sl.banks[3].h0;
                const int pos0 = headPos_[0], pos1 = headPos_[1], mask = headMask_;
                float hLL = 0.0f, hLR = 0.0f, hRL = 0.0f, hRR = 0.0f;
                for (int i = 0; i < P0_; ++i)
                {
                    const float l = headHist_[0][(std::size_t) ((pos0 - i) & mask)];
                    const float r = headHist_[1][(std::size_t) ((pos1 - i) & mask)];
                    hLL += LL[(std::size_t) i] * l; hLR += LR[(std::size_t) i] * r;
                    hRL += RL[(std::size_t) i] * l; hRR += RR[(std::size_t) i] * r;
                }
                yl = hLL + hLR; yr = hRL + hRR;
                break;
            }
        }
        for (int st = 0; st < numStages_; ++st)
        {
            const int ph = hist_[(std::size_t) st].phase;
            yl += sl.tail[(std::size_t) st][0][(std::size_t) ph];
            yr += sl.tail[(std::size_t) st][1][(std::size_t) ph];
        }
        oL = yl; oR = yr;
    }

    // Feed one input sample into the head rings + every stage's current frame half (per channel).
    void writeInputs (const float* const* in, int s) noexcept
    {
        for (int ch = 0; ch < channels_; ++ch)
        {
            const float x = in[ch][s];
            headPos_[ch] = (headPos_[ch] + 1) & headMask_;
            headHist_[ch][(std::size_t) headPos_[ch]] = x;
            for (int st = 0; st < numStages_; ++st) hist_[(std::size_t) st].frame[ch][(std::size_t) (spec_[st].B + hist_[(std::size_t) st].phase)] = x;
        }
    }

    // Advance every stage's phase; a chunk boundary recomputes tails (both slots while fading) + warms the FDL.
    void advanceStages (bool fading) noexcept
    {
        for (int st = 0; st < numStages_; ++st)
            if (++hist_[(std::size_t) st].phase == spec_[st].B) { hist_[(std::size_t) st].phase = 0; chunkStage (st, fading); }
    }

    // Idle path: samples [a,b) through the single live operator.
    void processRange (const float* const* in, float* const* out, int a, int b) noexcept
    {
        const Slot& sl = slot_[cur_];
        for (int s = a; s < b; ++s)
        {
            writeInputs (in, s);
            float oL = 0.0f, oR = 0.0f; computeSlotOutputs (sl, oL, oR);
            out[0][s] = oL; if (! mono_) out[1][s] = oR;
            advanceStages (false);
        }
    }

    // Crossfade path: blend the old (cur_) and new (1-cur_) operators by a smoothstep weight — both computed
    // from the SAME shared FDL — then finish the block on the single new operator once the fade ends.
    void processFade (const float* const* in, float* const* out, int n) noexcept
    {
        for (int s = 0; s < n; ++s)
        {
            const float t    = (float) xfadePos_ / (float) (xfadeLen_ > 1 ? xfadeLen_ - 1 : 1);   // 0→1, exactly 1 on the last sample
            const float wNew = t * t * (3.0f - 2.0f * t);   // smoothstep: zero slope at both ends (click-free)
            const float wOld = 1.0f - wNew;
            writeInputs (in, s);
            float oLo = 0.0f, oRo = 0.0f, oLn = 0.0f, oRn = 0.0f;
            computeSlotOutputs (slot_[cur_],     oLo, oRo);
            computeSlotOutputs (slot_[1 - cur_], oLn, oRn);
            out[0][s] = oLo * wOld + oLn * wNew;
            if (! mono_) out[1][s] = oRo * wOld + oRn * wNew;
            advanceStages (true);
            if (++xfadePos_ >= xfadeLen_)
            {
                cur_ = 1 - cur_;                                // new operator live; old now free to re-stage
                state_.store (0, std::memory_order_release);
                processRange (in, out, s + 1, n);              // finish the block on the single new operator
                return;
            }
        }
    }

    std::array<History, kMaxStages> hist_ {};
    StageSpec spec_[kMaxStages] {};
    Slot slot_[2];
    std::vector<float> headHist_[2];                               // P0 ring per channel (time domain)
    core::fft::AlignedVector<float> inputSpec_, viewSpec_, acc_, ifftOut_, buildBuf_;   // shared scratch (seam), sized to the largest stage
    std::vector<float> tmpTailA_, tmpTailB_;                                            // max-B M/S decode scratch (time domain)
    static_assert (std::atomic<int>::is_always_lock_free, "state_ must be lock-free — it is read/written on the audio thread");
    std::atomic<int> state_ { 0 };                                 // 0 Idle · 1 Pending (staged) · 2 Crossfading
    bool prepared_ = false;
    bool mono_ = false;
    int P0_ = 0, headMask_ = 0, headPos_[2] { 0, 0 };
    int numStages_ = 0, channels_ = 2;
    int cur_ = 0, xfadePos_ = 0, xfadeLen_ = 1, warmXfade_ = 1;
};

} // namespace felitronics::convolution

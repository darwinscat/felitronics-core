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
// COST NOTE. EVERY stage runs every sample and FFTs on its boundary REGARDLESS of the live IR length — this
// keeps every stage's FDL warm so a swap to a LONGER IR finds real history (the MatrixConvolver invariant). So
// CPU is set by maxIrSamples, not the current IR: a short-IR operator pays full price (unlike NonUniformConvolver,
// which skips inactive stages). A future "optimize short IRs" pass must NOT simply gate stages off — that would
// cold-start the FDL and break warm swaps.
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
                h.fdl[ch].assign   ((std::size_t) C * (std::size_t) specF, 0.0f);
            }
            h.phase = 0; h.fdlPos = 0;
        }

        const int nb = mono_ ? 1 : kMaxBanks;
        for (int k = 0; k < 2; ++k)
        {
            Slot& sl = slot_[k];
            sl.topo = Topology::LRDiag; sl.numBanks = 0;
            for (int b = 0; b < nb; ++b)
            {
                sl.banks[b].h0.assign ((std::size_t) P0_, 0.0f);
                for (int st = 0; st < numStages_; ++st)
                {
                    sl.banks[b].irSpec[(std::size_t) st].assign ((std::size_t) spec_[st].C * (std::size_t) spec_[st].specF, 0.0f);
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

    // Flush the running history (keeps the live operator). Audio thread; must not race setOperator().
    void reset() noexcept
    {
        for (int ch = 0; ch < channels_; ++ch) { std::fill (headHist_[ch].begin(), headHist_[ch].end(), 0.0f); headPos_[ch] = 0; }
        for (int st = 0; st < numStages_; ++st)
        {
            History& h = hist_[(std::size_t) st];
            for (int ch = 0; ch < channels_; ++ch) { std::fill (h.frame[ch].begin(), h.frame[ch].end(), 0.0f); std::fill (h.fdl[ch].begin(), h.fdl[ch].end(), 0.0f); }
            h.phase = 0; h.fdlPos = 0;
        }
        for (int k = 0; k < 2; ++k)
            for (int st = 0; st < numStages_; ++st)
                for (int ch = 0; ch < 2; ++ch) std::fill (slot_[k].tail[(std::size_t) st][(std::size_t) ch].begin(), slot_[k].tail[(std::size_t) st][(std::size_t) ch].end(), 0.0f);
        xfadePos_ = 0;                                 // cancel any pending/active fade
        state_.store (0, std::memory_order_relaxed);   // keep cur_ (the live operator)
    }

    static constexpr int latencySamples() noexcept { return 0; }
    bool isBusy() const noexcept { return state_.load (std::memory_order_acquire) != 0; }
    int  numChannels() const noexcept { return channels_; }

    // Message thread — single producer, must not race reset(). Stage a whole operator into the inactive slot,
    // then publish; process() goes live INSTANTLY (Phase 2 — the click-free crossfade lands in Phase 4).
    bool setOperator (Topology topo, const float* const* banks, int numBanks, int len)
    {
        if (! stageOperator (topo, banks, numBanks, len)) return false;
        publishStaged();
        return true;
    }

    // Build the operator into the inactive slot WITHOUT publishing (the surround ST-only path stages every
    // instance then publishes back-to-back). All topologies (mono / LRDiag / MSDiag / Full) are supported.
    bool stageOperator (Topology topo, const float* const* banks, int numBanks, int len)
    {
        if (! prepared_ || banks == nullptr || state_.load (std::memory_order_acquire) != 0) return false;
        const Topology t = mono_ ? Topology::LRDiag : topo;
        const int nb = mono_ ? 1 : numBanksFor (t);
        if (numBanks < nb) return false;
        const int stg = 1 - cur_;
        Slot& sl = slot_[stg];
        sl.topo = t; sl.numBanks = nb;
        // A null bank pointer is only rejected when it would be dereferenced (len > 0); (nullptr, len<=0) is the
        // documented "clear to silence" case buildBankSlice handles without a deref — matching MatrixConvolver.
        for (int b = 0; b < nb; ++b) { if (banks[b] == nullptr && len > 0) return false; buildBankSlice (sl.banks[b], banks[b], len); }
        return true;
    }

    void publishStaged() noexcept { state_.store (1, std::memory_order_release); }

    // Mono convenience: broadcast a single IR onto both output channels (or the one bank when prepared mono).
    bool setIr (const float* ir, int len)
    {
        const float* both[2] { ir, ir };
        return setOperator (Topology::LRDiag, both, mono_ ? 1 : 2, len);
    }

    // Audio thread. Planar in/out may alias (in-place). RT-safe, zero latency. Needs >= channels_ planes.
    void process (const float* const* in, float* const* out, int numChannelsToProcess, int n) noexcept
    {
        if (! prepared_ || n <= 0 || numChannelsToProcess < channels_) return;

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

        if (s != 2) { processRange (in, out, 0, n); return; }   // Idle → single live operator
        processFade (in, out, n);
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

    // Build one bank's head + per-stage partition spectra from a time-domain IR slice. Message thread.
    void buildBankSlice (Bank& bank, const float* ir, int len) noexcept
    {
        if (len < 0) len = 0;
        for (int i = 0; i < P0_; ++i) bank.h0[(std::size_t) i] = (i < len) ? ir[i] : 0.0f;
        for (int st = 0; st < numStages_; ++st)
        {
            const int B = spec_[st].B, C = spec_[st].C, offset = spec_[st].offset, specF = spec_[st].specF, N = spec_[st].N;
            int parts = 0;
            for (int j = 0; j < C; ++j)
            {
                const int base = offset + j * B;
                if (base >= len) break;                                          // this + later partitions all-zero
                std::fill (buildBuf_.begin(), buildBuf_.begin() + N, 0.0f);
                for (int i = 0; i < B; ++i) { const int src = base + i; buildBuf_[(std::size_t) i] = (src < len) ? ir[src] : 0.0f; }
                hist_[(std::size_t) st].buildFft.forward (buildBuf_.data(), &bank.irSpec[(std::size_t) st][(std::size_t) j * (std::size_t) specF]);
                parts = j + 1;
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

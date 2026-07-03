// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/core/Fft.h>
#include <felitronics/core/Config.h>

#include <algorithm>
#include <atomic>
#include <cstring>
#include <vector>

namespace felitronics::convolution
{

//==============================================================================
// felitronics::convolution::MatrixConvolver — a CLICK-FREE stereo 2×2 OPERATOR convolver on ONE canonical
// raw-L/R history. It is the matrix sibling of ConvolutionEngine: same proven zero-latency partitioning
// (head + partitioned tail) and the same one-shared-history / two-slot crossfade, but with cross-input
// routing so it can render the full 2×2 EQ transfer matrix that split-lane FIR EQ needs.
//
// ONE CANONICAL FDL BASIS — RAW L/R. The per-channel input frequency-domain delay line always stores the
// FFTs of past raw L and R input blocks (IR-independent, exactly like ConvolutionEngine's FDL). Every
// topology reads THAT shared history, which is what lets a topology CHANGE mid-playback (2-bank → 4-bank)
// crossfade with no cold-start tail — both operators are computed from the same warm FDL.
//
// THREE ROUTING TOPOLOGIES over that shared L/R history (an "operator" = {topology, IR bank set}):
//   • MSDiag (2 banks M,S): forms the M/S SPECTRAL VIEWS on the fly from the L/R FDL by linearity of the
//       FFT — X_M(k) = ½(X_L+X_R), X_S(k) = ½(X_L−X_R) — convolves bank_M against the M view and bank_S
//       against the S view, decodes in time yL = yM+yS, yR = yM−yS. Same math + 2-convolution cost as a
//       classic time-domain M/S encode→convolve→decode, but the FDL basis stays raw L/R.
//   • LRDiag (2 banks L,R): direct — yL = bank_L ∗ xL, yR = bank_R ∗ xR.
//   • Full   (4 banks LL,LR,RL,RR): yL = bank_LL∗xL + bank_LR∗xR, yR = bank_RL∗xL + bank_RR∗xR (≈2× cost).
//
// ATOMIC OPERATOR SWAP with a warm-history crossfade. setOperator() stages a whole operator (topology +
// its banks) into the inactive slot; process() crossfades the OUTPUTS of the old and new operators, both
// computed from the SAME shared FDL. It never mixes an old bank with a new one, and a topology change is
// just another operator swap. The very first activation (empty history) uses a long fade to mask the cold
// prime, exactly like ConvolutionEngine; every later swap uses the short anti-click fade.
//
// MONO (numChannels==1) degenerates to ONE FDL + ONE bank (topology ignored): yL = bank0 ∗ xL.
//
// CONTRACTS mirror ConvolutionEngine: prepare() allocates (message thread, RT-UNSAFE); setOperator()/setIr()
// are RT-UNSAFE single-producer message-thread (must not race reset()); process() is RT-safe (no
// alloc/lock/throw), zero latency; isBusy() is true while a swap is crossfading (the host coalesces);
// reset() flushes history (re-arms a cold prime) but keeps the live operator. The FFT is the compile-time
// backend seam (no hot-path vtable). Stereo only (2 output channels); the class is intentionally not
// channel-parameterised — the 2×2 matrix IS its reason to exist.
template <core::fft::RealFftBackend Fft = core::fft::DefaultRealFft>
class MatrixConvolver
{
public:
    enum class Topology { MSDiag, LRDiag, Full };

    // Bank layout per topology (indices into the operator's bank set):
    //   MSDiag → { M, S }      LRDiag → { L, R }      Full → { LL, LR, RL, RR }      mono → { single }
    static constexpr int kMaxBanks = 4;

    static int numBanksFor (Topology t) noexcept { return t == Topology::Full ? 4 : 2; }

    // partitionSize P (pow2; FFT size = 2P). maxIrSamples sizes the partition arrays. crossfadeSamples is
    // the SHORT (warm) anti-click fade; the cold first-prime length is derived internally. numChannels is 1
    // (mono) or 2 (stereo). Message thread; allocates.
    bool prepare (int partitionSize, int maxIrSamples, int crossfadeSamples, int numChannels = 2)
    {
        prepared_ = false;                            // any early return below leaves the engine unprepared
        if (numChannels < 1 || numChannels > 2) return false;
        if (! core::fft::isPow2 (partitionSize)) return false;
        channels_ = numChannels;
        mono_     = (channels_ == 1);
        P_ = partitionSize;
        N_ = 2 * P_;
        if (! fft_.prepare (N_)) return false;
        if (! buildFft_.prepare (N_)) return false;   // separate message-thread FFT for the IR build (no race with the audio fft_)
        specF_    = Fft::spectrumFloats (N_);
        maxParts_ = (maxIrSamples > P_) ? ((maxIrSamples - P_ + P_ - 1) / P_) : 0;

        warmXfade_ = crossfadeSamples < 1 ? 1 : crossfadeSamples;
        coldXfade_ = maxParts_ > 0 ? std::max (warmXfade_, maxParts_ * P_) : warmXfade_;

        inputSpec_.assign ((std::size_t) specF_, 0.0f);
        viewSpec_.assign  ((std::size_t) specF_, 0.0f);
        acc_.assign       ((std::size_t) specF_, 0.0f);
        ifftOut_.assign   ((std::size_t) N_,     0.0f);
        tmpTailA_.assign  ((std::size_t) P_,     0.0f);
        tmpTailB_.assign  ((std::size_t) P_,     0.0f);

        frameL_.assign ((std::size_t) N_, 0.0f);
        fdlL_.assign   ((std::size_t) maxParts_ * (std::size_t) specF_, 0.0f);
        if (! mono_)
        {
            frameR_.assign ((std::size_t) N_, 0.0f);
            fdlR_.assign   ((std::size_t) maxParts_ * (std::size_t) specF_, 0.0f);
        }
        for (int k = 0; k < 2; ++k) slot_[k].prepare (P_, maxParts_, specF_, mono_ ? 1 : kMaxBanks);

        cur_ = 0; xfadePos_ = 0; xfadeLen_ = warmXfade_;
        phase_ = 0; fdlPos_ = 0; warmSamples_ = 0;
        state_.store (0, std::memory_order_relaxed);
        prepared_ = true;
        return true;
    }

    // Flush the running history (next swap re-primes cold) but KEEP the live operator. Cancels a pending
    // swap. Audio thread (or externally synced); must not run concurrently with setOperator().
    void reset() noexcept
    {
        std::fill (frameL_.begin(), frameL_.end(), 0.0f);
        std::fill (fdlL_.begin(),   fdlL_.end(),   0.0f);
        if (! mono_) { std::fill (frameR_.begin(), frameR_.end(), 0.0f); std::fill (fdlR_.begin(), fdlR_.end(), 0.0f); }
        for (int k = 0; k < 2; ++k) slot_[k].resetTails();
        xfadePos_ = 0; phase_ = 0; fdlPos_ = 0; warmSamples_ = 0;   // keep cur_ (the live operator)
        state_.store (0, std::memory_order_relaxed);
    }

    static constexpr int latencySamples() noexcept { return 0; }
    bool isBusy() const noexcept { return state_.load (std::memory_order_acquire) != 0; }
    int  numChannels() const noexcept { return channels_; }

    // Message thread — SINGLE producer; must not run concurrently with reset(). Stage a whole operator
    // (topology + its banks) into the INACTIVE slot; process() crossfades old→new over the shared warm
    // history. `banks` holds numBanksFor(topo) IR pointers in the layout above; every bank shares `len`.
    // Returns false if a swap is already pending/crossfading (the caller coalesces with the latest snapshot).
    bool setOperator (Topology topo, const float* const* banks, int numBanks, int len)
    {
        if (! stageOperator (topo, banks, numBanks, len)) return false;
        publishStaged();
        return true;
    }

    // Two-phase variant: build the operator into the inactive slot WITHOUT publishing. A host running one
    // convolver per channel (the surround ST-only path) stages every instance first (the expensive bank
    // FFTs), then publishes them back-to-back — so all channels enter their crossfade within nanoseconds
    // of each other instead of skewing by whole audio blocks while later channels are still building.
    bool stageOperator (Topology topo, const float* const* banks, int numBanks, int len)
    {
        if (! prepared_ || banks == nullptr || state_.load (std::memory_order_acquire) != 0) return false;
        const int stg = 1 - cur_;                                          // inactive slot
        const int nb  = mono_ ? 1 : std::min (numBanks, numBanksFor (topo));
        if (nb < (mono_ ? 1 : numBanksFor (topo))) return false;           // not enough banks supplied
        slot_[stg].topo     = mono_ ? Topology::LRDiag : topo;             // mono: topology is inert
        slot_[stg].numBanks = mono_ ? 1 : numBanksFor (topo);
        for (int b = 0; b < slot_[stg].numBanks; ++b)
            slot_[stg].banks[b].build (banks[b], len, P_, maxParts_, specF_, buildFft_);
        return true;
    }
    // KNOWN residual: across N instances the publishes are N independent atomic stores — if the audio
    // thread interleaves mid-loop, one channel's fade can start up to ONE block before another's. With
    // identical IRs on every channel and a ≥20 ms smoothstep fade this is inaudible; a shared fade clock
    // across instances was judged over-engineering (revisit only if a real bed swap ever images).
    void publishStaged() noexcept { state_.store (1, std::memory_order_release); }   // → Pending

    // Mono convenience: a single IR broadcast onto the one bank.
    bool setIr (const float* ir, int len)
    {
        const float* one[1] { ir };
        return setOperator (Topology::LRDiag, one, 1, len);
    }

    // Audio thread. Planar `in`/`out` may alias (in-place). RT-safe, zero latency. The matrix width is
    // FIXED at prepare() (1 or 2) — the 2×2 routing produces two outputs from two inputs — so process()
    // needs at least `channels_` planes and always operates on exactly that many; an under-width call is a
    // safe no-op (the documented contract: re-prepare on a bus change, as the adapter does).
    void process (const float* const* in, float* const* out, int numChannelsToProcess, int n) noexcept
    {
        if (! prepared_ || n <= 0 || numChannelsToProcess < channels_) return;

        int s = state_.load (std::memory_order_acquire);
        if (s == 1)   // begin crossfade — length by warmth (cold first-prime vs short warm swap)
        {
            xfadePos_ = 0;
            const int newParts = slot_[1 - cur_].maxNumParts();
            xfadeLen_ = (warmSamples_ >= coldXfade_) ? warmXfade_
                                                     : std::clamp ((newParts + 1) * P_, warmXfade_, coldXfade_);
            if (maxParts_ > 0)                                             // prime the new slot's tail from the warm FDL
            {
                int base = fdlPos_ - 1; if (base < 0) base += maxParts_;
                computeTails (slot_[1 - cur_], base);
            }
            else slot_[1 - cur_].resetTails();
            state_.store (2, std::memory_order_relaxed);
            s = 2;
        }

        if (s != 2) { processRange (in, out, 0, n); return; }             // Idle → single active operator
        processFade (in, out, n);
    }

private:
    //==========================================================================
    // A spectral IR bank: head taps (partition 0, time domain) + tail partition spectra + partition count.
    // Holds IR data only — reset() does NOT clear it (history lives in the FDL / frames, not here).
    struct Bank
    {
        std::vector<float> h0;        // P: direct (head) taps
        std::vector<float> irSpec;    // maxParts × specF: tail partition spectra
        int numParts = 0;

        void prepare (int P, int maxParts, int specF)
        {
            h0.assign     ((std::size_t) P, 0.0f);
            irSpec.assign ((std::size_t) maxParts * (std::size_t) specF, 0.0f);
            numParts = 0;
        }

        // Message thread: build this bank's head + tail spectra from a time-domain IR (allocates a scratch).
        void build (const float* ir, int len, int P, int maxParts, int specF, Fft& fft)
        {
            if (P <= 0) return;                        // unprepared — also guards the /P below
            if (len < 0) len = 0;
            for (int i = 0; i < P; ++i) h0[(std::size_t) i] = (i < len) ? ir[i] : 0.0f;

            const int tailLen = (len > P) ? (len - P) : 0;
            int parts = (tailLen > 0) ? ((tailLen + P - 1) / P) : 0;
            if (parts > maxParts) parts = maxParts;

            std::vector<float> part ((std::size_t) (2 * P), 0.0f);
            for (int j = 0; j < parts; ++j)
            {
                std::fill (part.begin(), part.end(), 0.0f);
                for (int i = 0; i < P; ++i)
                {
                    const int src = P + j * P + i;
                    part[(std::size_t) i] = (src < len) ? ir[src] : 0.0f;
                }
                fft.forward (part.data(), &irSpec[(std::size_t) j * (std::size_t) specF]);
            }
            numParts = parts;
        }
    };

    // An operator slot: topology + its banks + the cached tail outputs (per output channel yL/yR).
    struct Slot
    {
        Topology topo = Topology::MSDiag;
        int numBanks = 0;
        Bank banks[kMaxBanks];
        std::vector<float> tailL, tailR;   // P each — this operator's cached tail for the NEXT chunk

        void prepare (int P, int maxParts, int specF, int nBanks)
        {
            for (int b = 0; b < nBanks; ++b) banks[b].prepare (P, maxParts, specF);   // mono: bank 0 only — a
            // mono operator never touches banks[1..3] (setOperator caps numBanks), so their spectra stay unallocated
            tailL.assign ((std::size_t) P, 0.0f);
            tailR.assign ((std::size_t) P, 0.0f);
        }
        void resetTails() noexcept
        {
            std::fill (tailL.begin(), tailL.end(), 0.0f);
            std::fill (tailR.begin(), tailR.end(), 0.0f);
        }
        int maxNumParts() const noexcept
        {
            int m = 0; for (int b = 0; b < numBanks; ++b) m = std::max (m, banks[b].numParts); return m;
        }
    };

    static float headDot (const std::vector<float>& h0, const float* fr, int P) noexcept
    {
        float h = 0.0f;
        for (int i = 0; i < P; ++i) h += h0[(std::size_t) i] * fr[-i];
        return h;
    }

    void clearAcc() noexcept { std::fill (acc_.begin(), acc_.end(), 0.0f); }

    // IFFT the accumulated spectrum, overlap-save the last P into `dst`.
    void invAccTo (std::vector<float>& dst) noexcept
    {
        fft_.inverse (acc_.data(), ifftOut_.data());
        for (int i = 0; i < P_; ++i) dst[(std::size_t) i] = ifftOut_[(std::size_t) (P_ + i)];
    }

    // acc_ += Σ_j srcFdl[(base-j) wrap] .* bank.irSpec[j]  (direct — L or R history against one bank)
    void macBank (const std::vector<float>& srcFdl, const Bank& bank, int base) noexcept
    {
        for (int j = 0; j < bank.numParts; ++j)
        {
            int idx = base - j; if (idx < 0) idx += maxParts_;
            fft_.spectralMultiplyAdd (&srcFdl[(std::size_t) idx * (std::size_t) specF_],
                                      &bank.irSpec[(std::size_t) j * (std::size_t) specF_], acc_.data());
        }
    }

    // acc_ += Σ_j (½(X_L ± X_R))[base-j] .* bank.irSpec[j]  — the on-the-fly M/S spectral view (sign +→M, −→S).
    void macView (const Bank& bank, int base, float sign) noexcept
    {
        for (int j = 0; j < bank.numParts; ++j)
        {
            int idx = base - j; if (idx < 0) idx += maxParts_;
            const float* xl = &fdlL_[(std::size_t) idx * (std::size_t) specF_];
            const float* xr = &fdlR_[(std::size_t) idx * (std::size_t) specF_];
            for (int i = 0; i < specF_; ++i) viewSpec_[(std::size_t) i] = 0.5f * (xl[i] + sign * xr[i]);
            fft_.spectralMultiplyAdd (viewSpec_.data(),
                                      &bank.irSpec[(std::size_t) j * (std::size_t) specF_], acc_.data());
        }
    }

    // Fill a slot's cached tail outputs (yL/yR) from the FDL history ending at partition `base`, per topology.
    void computeTails (Slot& sl, int base) noexcept
    {
        if (mono_)
        {
            clearAcc(); macBank (fdlL_, sl.banks[0], base); invAccTo (sl.tailL);
            return;
        }
        switch (sl.topo)
        {
            case Topology::MSDiag:
                clearAcc(); macView (sl.banks[0], base, +1.0f); invAccTo (tmpTailA_);   // yM tail
                clearAcc(); macView (sl.banks[1], base, -1.0f); invAccTo (tmpTailB_);   // yS tail
                for (int i = 0; i < P_; ++i)
                {
                    const float m = tmpTailA_[(std::size_t) i], s = tmpTailB_[(std::size_t) i];
                    sl.tailL[(std::size_t) i] = m + s;   // decode yL = yM + yS
                    sl.tailR[(std::size_t) i] = m - s;   //        yR = yM − yS
                }
                break;

            case Topology::LRDiag:
                clearAcc(); macBank (fdlL_, sl.banks[0], base); invAccTo (sl.tailL);
                clearAcc(); macBank (fdlR_, sl.banks[1], base); invAccTo (sl.tailR);
                break;

            case Topology::Full:
                clearAcc(); macBank (fdlL_, sl.banks[0], base); macBank (fdlR_, sl.banks[1], base); invAccTo (sl.tailL);
                clearAcc(); macBank (fdlL_, sl.banks[2], base); macBank (fdlR_, sl.banks[3], base); invAccTo (sl.tailR);
                break;
        }
    }

    // One sample's output from a slot: head (partition 0, time domain) + the cached tail at `phase`.
    void computeOutputs (const Slot& sl, int phase, const float* frL, const float* frR,
                         float& outL, float& outR) const noexcept
    {
        if (mono_)
        {
            outL = headDot (sl.banks[0].h0, frL, P_) + sl.tailL[(std::size_t) phase];
            outR = 0.0f;
            return;
        }
        switch (sl.topo)
        {
            case Topology::MSDiag:
            {
                const std::vector<float>& hM = sl.banks[0].h0;
                const std::vector<float>& hS = sl.banks[1].h0;
                float headM = 0.0f, headS = 0.0f;                         // form m/s on the fly by linearity
                for (int i = 0; i < P_; ++i)
                {
                    const float l = frL[-i], rr = frR[-i];
                    headM += hM[(std::size_t) i] * 0.5f * (l + rr);
                    headS += hS[(std::size_t) i] * 0.5f * (l - rr);
                }
                outL = headM + headS + sl.tailL[(std::size_t) phase];     // tailL already = yM+yS
                outR = headM - headS + sl.tailR[(std::size_t) phase];
                break;
            }
            case Topology::LRDiag:
                outL = headDot (sl.banks[0].h0, frL, P_) + sl.tailL[(std::size_t) phase];
                outR = headDot (sl.banks[1].h0, frR, P_) + sl.tailR[(std::size_t) phase];
                break;

            case Topology::Full:
            {
                const std::vector<float>& LL = sl.banks[0].h0;
                const std::vector<float>& LR = sl.banks[1].h0;
                const std::vector<float>& RL = sl.banks[2].h0;
                const std::vector<float>& RR = sl.banks[3].h0;
                float hLL = 0.0f, hLR = 0.0f, hRL = 0.0f, hRR = 0.0f;
                for (int i = 0; i < P_; ++i)
                {
                    const float l = frL[-i], rr = frR[-i];
                    hLL += LL[(std::size_t) i] * l; hLR += LR[(std::size_t) i] * rr;
                    hRL += RL[(std::size_t) i] * l; hRR += RR[(std::size_t) i] * rr;
                }
                outL = hLL + hLR + sl.tailL[(std::size_t) phase];
                outR = hRL + hRR + sl.tailR[(std::size_t) phase];
                break;
            }
        }
    }

    // A full P-sample chunk just completed. FFT each input frame once into the SHARED FDL, recompute the
    // active operator's tail (and the incoming operator's too while fading), then shift current→previous.
    void chunkAll (bool fading) noexcept
    {
        fft_.forward (frameL_.data(), inputSpec_.data());
        if (maxParts_ > 0)
            std::memcpy (&fdlL_[(std::size_t) fdlPos_ * (std::size_t) specF_], inputSpec_.data(),
                         (std::size_t) specF_ * sizeof (float));
        if (! mono_)
        {
            fft_.forward (frameR_.data(), inputSpec_.data());
            if (maxParts_ > 0)
                std::memcpy (&fdlR_[(std::size_t) fdlPos_ * (std::size_t) specF_], inputSpec_.data(),
                             (std::size_t) specF_ * sizeof (float));
        }

        if (maxParts_ > 0)
        {
            computeTails (slot_[cur_], fdlPos_);
            if (fading) computeTails (slot_[1 - cur_], fdlPos_);
        }
        else
        {
            slot_[cur_].resetTails();
            if (fading) slot_[1 - cur_].resetTails();
        }

        for (int i = 0; i < P_; ++i) frameL_[(std::size_t) i] = frameL_[(std::size_t) (P_ + i)];
        if (! mono_) for (int i = 0; i < P_; ++i) frameR_[(std::size_t) i] = frameR_[(std::size_t) (P_ + i)];
        if (maxParts_ > 0) { if (++fdlPos_ >= maxParts_) fdlPos_ = 0; }
    }

    // Idle path: samples [a,b) through the single active operator cur_.
    void processRange (const float* const* in, float* const* out, int a, int b) noexcept
    {
        const Slot& sl = slot_[cur_];
        for (int s = a; s < b; ++s)
        {
            frameL_[(std::size_t) (P_ + phase_)] = in[0][s];
            if (! mono_) frameR_[(std::size_t) (P_ + phase_)] = in[1][s];
            const float* frL = &frameL_[(std::size_t) (P_ + phase_)];
            const float* frR = mono_ ? nullptr : &frameR_[(std::size_t) (P_ + phase_)];
            float oL = 0.0f, oR = 0.0f;   // zero-init: the mono path leaves oR unwritten (GCC -Wmaybe-uninitialized at -O2)
            computeOutputs (sl, phase_, frL, frR, oL, oR);
            out[0][s] = oL;
            if (! mono_) out[1][s] = oR;
            if (++phase_ == P_) { phase_ = 0; chunkAll (false); }
            if (warmSamples_ < coldXfade_) ++warmSamples_;
        }
    }

    // Crossfade path: blend the old operator (cur_) and the new (1-cur_) by a smoothstep weight, both
    // computed from the SAME shared FDL, then finish the block on the single new operator once the fade ends.
    void processFade (const float* const* in, float* const* out, int n) noexcept
    {
        const Slot& sOld = slot_[cur_];
        const Slot& sNew = slot_[1 - cur_];
        for (int s = 0; s < n; ++s)
        {
            const float t    = (float) xfadePos_ / (float) (xfadeLen_ > 1 ? xfadeLen_ - 1 : 1);   // 0→1, exactly 1 on the last sample
            const float wNew = t * t * (3.0f - 2.0f * t);   // smoothstep: zero slope at both ends (click-free)
            const float wOld = 1.0f - wNew;
            frameL_[(std::size_t) (P_ + phase_)] = in[0][s];
            if (! mono_) frameR_[(std::size_t) (P_ + phase_)] = in[1][s];
            const float* frL = &frameL_[(std::size_t) (P_ + phase_)];
            const float* frR = mono_ ? nullptr : &frameR_[(std::size_t) (P_ + phase_)];
            float oLo = 0.0f, oRo = 0.0f, oLn = 0.0f, oRn = 0.0f;   // zero-init (see processRange)
            computeOutputs (sOld, phase_, frL, frR, oLo, oRo);
            computeOutputs (sNew, phase_, frL, frR, oLn, oRn);
            out[0][s] = oLo * wOld + oLn * wNew;
            if (! mono_) out[1][s] = oRo * wOld + oRn * wNew;
            if (++phase_ == P_) { phase_ = 0; chunkAll (true); }
            if (warmSamples_ < coldXfade_) ++warmSamples_;
            if (++xfadePos_ >= xfadeLen_)
            {
                cur_ = 1 - cur_;                                // new operator live; old now free to re-stage
                state_.store (0, std::memory_order_release);    // → Idle
                processRange (in, out, s + 1, n);               // finish the block on the single new operator
                return;
            }
        }
    }

    Fft fft_;                                            // audio-thread FFT (chunkAll / computeTails)
    Fft buildFft_;                                        // message-thread FFT (Bank::build) — separate, no race
    Slot slot_[2];
    std::vector<float> frameL_, frameR_;                  // 2P each: [prev | current] input accumulators
    std::vector<float> fdlL_, fdlR_;                      // maxParts × specF: rings of past raw L/R input spectra
    std::vector<float> inputSpec_, viewSpec_, acc_, ifftOut_, tmpTailA_, tmpTailB_;   // shared FFT scratch
    std::atomic<int> state_ { 0 };                        // 0 Idle · 1 Pending (staged) · 2 Crossfading
    bool prepared_ = false;
    bool mono_ = false;
    int P_ = 0, N_ = 0, specF_ = 0, maxParts_ = 0, channels_ = 2;
    int cur_ = 0;                                         // active slot (0/1)
    int phase_ = 0, fdlPos_ = 0;                          // shared per-chunk timing (both channels lockstep)
    int xfadePos_ = 0, xfadeLen_ = 1, warmXfade_ = 1, coldXfade_ = 1;
    long long warmSamples_ = 0;                           // samples processed (saturates at coldXfade_) → warm test
};

} // namespace felitronics::convolution

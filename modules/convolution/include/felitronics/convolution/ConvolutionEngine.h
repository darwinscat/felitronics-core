// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa. Part of felitronics-core — see LICENSE.

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
// felitronics::convolution::ConvolutionEngine — the production, swap-safe, CLICK-FREE convolver,
// MULTI-CHANNEL with a single LOCKSTEP crossfade. Zero latency (head + partitioned-tail, like
// PartitionedConvolver, which it reimplements internally with dual coefficients).
//
// DESIGN B — ONE SHARED warm history, TWO IR slots. Per channel there is a SINGLE running input
// history (frame_ + the frequency-domain delay line `fdl`) — and that history is IR-INDEPENDENT (the
// FDL stores FFTs of past INPUT blocks). On top of it sit TWO IR coefficient slots. setIr() stages the
// new IR into the INACTIVE slot and process() crossfades old→new over a SHORT window against the SAME
// warm history. So a fresh IR is response-correct from the first block (no per-swap re-prime): an
// interactive IR change — an EQ band drag — lands within that short fade, not a full-FIR-length ramp.
//
// COST: ~1× convolution normally; 2× only DURING the brief fade (the forward FFT into the FDL is shared
// once; only the cheap per-partition spectral MAC and the direct head dot run twice). The two slots ARE
// the swap double-buffer, so this is race-free WITHOUT keeping a second history warm: the message thread
// writes the inactive slot only while Idle and publishes with a release store; the audio thread reads
// both slots only after acquiring that store (during the fade), and commits the slot flip at fade end.
//
// COLD START: on the very first activation the shared history is empty, so that ONE swap uses a long
// fade (≈ the tail length) to mask the cold prime; once the history has filled, every swap uses the
// short (anti-click) fade. Detected by an audio-thread sample counter (no message-thread race).
//
// LOCKSTEP: one crossfade position advanced once per sample-frame across all channels, so a stereo IR
// swap (Mid IR ch0, Side IR ch1) can never move the image or decorrelate L/R. Race-free via a 3-state
// atomic (0 Idle / 1 Pending / 2 Crossfading). The audio thread NEVER allocates, locks, or blocks; the
// adapter coalesces setIr() while isBusy(). `MaxChannels` is the compile-time channel bound (2 = stereo
// default; a mono consumer uses <Fft, 1>). The FFT is the compile-time backend seam (no hot-path vtable).
template <core::fft::RealFftBackend Fft = core::fft::DefaultRealFft, int MaxChannels = 2>
class ConvolutionEngine
{
public:
    // partitionSize P (pow2; FFT size = 2P). maxIrSamples sizes the partition arrays. crossfadeSamples is
    // the SHORT (warm) anti-click fade; the cold first-prime length is derived internally. Message thread.
    bool prepare (int partitionSize, int maxIrSamples, int crossfadeSamples, int numChannels = 1)
    {
        if (numChannels < 1 || numChannels > MaxChannels) return false;
        if (! core::fft::isPow2 (partitionSize)) return false;
        channels_ = numChannels;
        P_ = partitionSize;
        N_ = 2 * P_;
        if (! fft_.prepare (N_)) return false;
        if (! buildFft_.prepare (N_)) return false;   // separate FFT for the message-thread IR build — no race with the audio fft_
        specF_    = Fft::spectrumFloats (N_);
        maxParts_ = (maxIrSamples > P_) ? ((maxIrSamples - P_ + P_ - 1) / P_) : 0;

        warmXfade_ = crossfadeSamples < 1 ? 1 : crossfadeSamples;
        coldXfade_ = maxParts_ > 0 ? maxParts_ * P_ : warmXfade_;   // cold prime ≈ the tail length

        inputSpec_.assign ((std::size_t) specF_, 0.0f);
        acc_.assign       ((std::size_t) specF_, 0.0f);
        ifftOut_.assign   ((std::size_t) N_,     0.0f);
        for (int c = 0; c < MaxChannels; ++c) chan_[c].prepare (P_, N_, specF_, maxParts_);

        cur_ = 0; xfadePos_ = 0; xfadeLen_ = warmXfade_;
        phase_ = 0; fdlPos_ = 0; warmSamples_ = 0;
        state_.store (0, std::memory_order_relaxed);
        return true;
    }

    // Clears the running history (so the next swap re-primes cold) but KEEPS the current IR slot — reset
    // means "flush the tail", not "revert the EQ". Cancels any pending swap. Audio thread (or externally
    // synced); must not run concurrently with setIr().
    void reset() noexcept
    {
        for (int c = 0; c < channels_; ++c) chan_[c].reset();
        xfadePos_ = 0; phase_ = 0; fdlPos_ = 0; warmSamples_ = 0;   // keep cur_ (the live IR slot)
        state_.store (0, std::memory_order_relaxed);
    }

    static constexpr int latencySamples() noexcept { return 0; }
    bool isBusy() const noexcept { return state_.load (std::memory_order_acquire) != 0; }
    int  numChannels() const noexcept { return channels_; }

    // Message thread — SINGLE producer; must not run concurrently with reset(). Per-channel IR (a
    // true-stereo IR); if fewer channels are supplied than configured, the last is broadcast to the
    // remainder. Stages into the INACTIVE slot (keeps the shared warm history). Returns false if a swap
    // is already pending/crossfading (the caller coalesces with the latest snapshot).
    bool setIr (const float* const* irPerCh, int nch, int len)
    {
        if (nch < 1 || state_.load (std::memory_order_acquire) != 0) return false;   // busy
        const int stg = 1 - cur_;                                                    // inactive slot
        for (int c = 0; c < channels_; ++c)
            chan_[c].buildIr (stg, irPerCh[c < nch ? c : nch - 1], len, P_, maxParts_, specF_, buildFft_);
        state_.store (1, std::memory_order_release);                                 // → Pending (publishes)
        return true;
    }

    // Message thread. Mono IR broadcast to ALL channels.
    bool setIr (const float* ir, int len)
    {
        const float* one[1] { ir };
        return setIr (one, 1, len);
    }

    // Audio thread. Planar; `in`/`out` may alias (in-place). RT-safe, zero latency.
    void process (const float* const* in, float* const* out, int numChannelsToProcess, int n) noexcept
    {
        const int nc = numChannelsToProcess < channels_ ? numChannelsToProcess : channels_;
        if (nc <= 0 || n <= 0) return;

        int s = state_.load (std::memory_order_acquire);
        if (s == 1)   // begin crossfade — length by warmth (cold first-prime vs short warm swap)
        {
            xfadePos_ = 0;
            xfadeLen_ = (warmSamples_ >= coldXfade_) ? warmXfade_ : coldXfade_;
            for (int c = 0; c < nc; ++c) primeTail (chan_[c], 1 - cur_);   // make the new slot's tail valid NOW from the warm FDL (else its zeroed tail leaks into the blend for ≤P samples)
            state_.store (2, std::memory_order_relaxed);
            s = 2;
        }

        if (s != 2) { processRange (in, out, nc, 0, n); return; }   // Idle → single active slot
        processFade (in, out, nc, n);
    }

    // Mono convenience (keeps a single-channel consumer on the 3-arg process()).
    void process (const float* in, float* out, int n) noexcept
    {
        const float* ins[1]  { in };
        float*       outs[1] { out };
        process (ins, outs, 1, n);
    }

private:
    //==========================================================================
    // Per-channel: ONE running input history shared by TWO IR coefficient slots.
    struct Chan
    {
        std::vector<float> frame;            // 2P: [previous chunk | current chunk] input accumulator
        std::vector<float> fdl;              // maxParts × specF: ring of past INPUT spectra (IR-independent)
        std::vector<float> pendingTail[2];   // P each: cached tail output per slot
        std::vector<float> h0[2];            // P each: direct (head) taps per slot
        std::vector<float> irSpec[2];        // maxParts × specF each: tail partition spectra per slot
        int numParts[2] { 0, 0 };

        void prepare (int P, int N, int specF, int maxParts)
        {
            frame.assign ((std::size_t) N, 0.0f);
            fdl.assign   ((std::size_t) maxParts * (std::size_t) specF, 0.0f);
            for (int k = 0; k < 2; ++k)
            {
                pendingTail[k].assign ((std::size_t) P, 0.0f);
                h0[k].assign          ((std::size_t) P, 0.0f);
                irSpec[k].assign      ((std::size_t) maxParts * (std::size_t) specF, 0.0f);
                numParts[k] = 0;
            }
        }

        void reset() noexcept
        {
            std::fill (frame.begin(), frame.end(), 0.0f);
            std::fill (fdl.begin(),   fdl.end(),   0.0f);
            for (int k = 0; k < 2; ++k) std::fill (pendingTail[k].begin(), pendingTail[k].end(), 0.0f);
        }

        // Message thread: build slot k's head + tail spectra (allocates a scratch). Zeroes that slot's
        // cached tail — it is recomputed from the warm shared FDL on the next chunk after the swap.
        void buildIr (int k, const float* ir, int len, int P, int maxParts, int specF, Fft& fft)
        {
            if (len < 0) len = 0;
            for (int i = 0; i < P; ++i) h0[k][(std::size_t) i] = (i < len) ? ir[i] : 0.0f;

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
                fft.forward (part.data(), &irSpec[k][(std::size_t) j * (std::size_t) specF]);
            }
            numParts[k] = parts;
            std::fill (pendingTail[k].begin(), pendingTail[k].end(), 0.0f);
        }
    };

    static float headDot (const std::vector<float>& h0, const float* fr, int P) noexcept
    {
        float h = 0.0f;
        for (int i = 0; i < P; ++i) h += h0[(std::size_t) i] * fr[-i];
        return h;
    }

    // One channel's tail MAC for slot k: Σ_j FDL[fdlPos-j] · irSpec[k][j] → IFFT → pendingTail[k].
    void macSlot (Chan& ch, int k) noexcept
    {
        if (ch.numParts[k] <= 0) { std::fill (ch.pendingTail[k].begin(), ch.pendingTail[k].end(), 0.0f); return; }
        std::fill (acc_.begin(), acc_.end(), 0.0f);
        for (int j = 0; j < ch.numParts[k]; ++j)
        {
            int idx = fdlPos_ - j;
            if (idx < 0) idx += maxParts_;
            fft_.spectralMultiplyAdd (&ch.fdl[(std::size_t) idx * (std::size_t) specF_],
                                      &ch.irSpec[k][(std::size_t) j * (std::size_t) specF_], acc_.data());
        }
        fft_.inverse (acc_.data(), ifftOut_.data());
        for (int i = 0; i < P_; ++i) ch.pendingTail[k][(std::size_t) i] = ifftOut_[(std::size_t) (P_ + i)];  // overlap-save: last P
    }

    // Prime slot k's pendingTail from the CURRENT warm FDL — base = the LAST-written chunk (fdlPos_-1),
    // matching the active slot's cached tail. Called at fade start so a freshly-staged slot's tail is
    // valid immediately, instead of zero for ≤P samples (which would leak into the blend). Audio thread.
    void primeTail (Chan& ch, int k) noexcept
    {
        if (ch.numParts[k] <= 0 || maxParts_ <= 0) { std::fill (ch.pendingTail[k].begin(), ch.pendingTail[k].end(), 0.0f); return; }
        int base = fdlPos_ - 1; if (base < 0) base += maxParts_;
        std::fill (acc_.begin(), acc_.end(), 0.0f);
        for (int j = 0; j < ch.numParts[k]; ++j)
        {
            int idx = base - j;
            if (idx < 0) idx += maxParts_;
            fft_.spectralMultiplyAdd (&ch.fdl[(std::size_t) idx * (std::size_t) specF_],
                                      &ch.irSpec[k][(std::size_t) j * (std::size_t) specF_], acc_.data());
        }
        fft_.inverse (acc_.data(), ifftOut_.data());
        for (int i = 0; i < P_; ++i) ch.pendingTail[k][(std::size_t) i] = ifftOut_[(std::size_t) (P_ + i)];
    }

    // A full P-sample chunk just completed for every channel. FFT each channel's frame into the SHARED
    // FDL once, MAC the active slot (and the other slot while fading), then shift current→previous.
    void chunkAll (int nc, bool fading) noexcept
    {
        const int other = 1 - cur_;
        for (int c = 0; c < nc; ++c)
        {
            Chan& ch = chan_[c];
            fft_.forward (ch.frame.data(), inputSpec_.data());
            if (maxParts_ > 0)
            {
                std::memcpy (&ch.fdl[(std::size_t) fdlPos_ * (std::size_t) specF_], inputSpec_.data(),
                             (std::size_t) specF_ * sizeof (float));
                macSlot (ch, cur_);
                if (fading) macSlot (ch, other);
            }
            else
            {
                std::fill (ch.pendingTail[cur_].begin(), ch.pendingTail[cur_].end(), 0.0f);
                if (fading) std::fill (ch.pendingTail[other].begin(), ch.pendingTail[other].end(), 0.0f);
            }
            for (int i = 0; i < P_; ++i) ch.frame[(std::size_t) i] = ch.frame[(std::size_t) (P_ + i)];  // current → previous
        }
        if (maxParts_ > 0) { if (++fdlPos_ >= maxParts_) fdlPos_ = 0; }
    }

    // Idle path: process samples [a, b) through the single active slot cur_.
    void processRange (const float* const* in, float* const* out, int nc, int a, int b) noexcept
    {
        for (int s = a; s < b; ++s)
        {
            for (int c = 0; c < nc; ++c)
            {
                Chan& ch = chan_[c];
                ch.frame[(std::size_t) (P_ + phase_)] = in[c][s];
                const float* fr = &ch.frame[(std::size_t) (P_ + phase_)];
                out[c][s] = headDot (ch.h0[cur_], fr, P_) + ch.pendingTail[cur_][(std::size_t) phase_];
            }
            if (++phase_ == P_) { phase_ = 0; chunkAll (nc, false); }
            if (warmSamples_ < coldXfade_) ++warmSamples_;
        }
    }

    // Crossfade path: blend slot cur_ (out) and 1-cur_ (in) by a smoothstep weight, then finish the
    // block on the single new slot once the fade completes.
    void processFade (const float* const* in, float* const* out, int nc, int n) noexcept
    {
        const int other = 1 - cur_;
        for (int s = 0; s < n; ++s)
        {
            const float t    = (float) xfadePos_ / (float) (xfadeLen_ > 1 ? xfadeLen_ - 1 : 1);   // 0→1, reaches exactly 1 on the last fade sample (no step at the hand-off)
            const float wNew = t * t * (3.0f - 2.0f * t);   // smoothstep: zero slope at both ends (click-free)
            const float wOld = 1.0f - wNew;
            for (int c = 0; c < nc; ++c)
            {
                Chan& ch = chan_[c];
                ch.frame[(std::size_t) (P_ + phase_)] = in[c][s];
                const float* fr  = &ch.frame[(std::size_t) (P_ + phase_)];
                const float  oOld = headDot (ch.h0[cur_],  fr, P_) + ch.pendingTail[cur_] [(std::size_t) phase_];
                const float  oNew = headDot (ch.h0[other], fr, P_) + ch.pendingTail[other][(std::size_t) phase_];
                out[c][s] = oOld * wOld + oNew * wNew;
            }
            if (++phase_ == P_) { phase_ = 0; chunkAll (nc, true); }
            if (warmSamples_ < coldXfade_) ++warmSamples_;
            if (++xfadePos_ >= xfadeLen_)
            {
                cur_ = other;                                   // new slot live; old now free to re-stage
                state_.store (0, std::memory_order_release);    // → Idle
                processRange (in, out, nc, s + 1, n);           // finish the block on the new single slot
                return;
            }
        }
    }

    Fft fft_;                                            // audio-thread FFT (chunkAll / macSlot / primeTail)
    Fft buildFft_;                                        // message-thread FFT (buildIr) — separate, so an IR build never races the audio FFT
    Chan chan_[MaxChannels];
    std::vector<float> inputSpec_, acc_, ifftOut_;        // shared FFT scratch (one channel at a time)
    std::atomic<int> state_ { 0 };                        // 0 Idle · 1 Pending (staged) · 2 Crossfading
    int P_ = 0, N_ = 0, specF_ = 0, maxParts_ = 0, channels_ = 1;
    int cur_ = 0;                                         // active slot (0/1)
    int phase_ = 0, fdlPos_ = 0;                          // shared per-chunk timing (all channels lockstep)
    int xfadePos_ = 0, xfadeLen_ = 1, warmXfade_ = 1, coldXfade_ = 1;
    long long warmSamples_ = 0;                           // samples processed (saturates at coldXfade_) → warm test
};

} // namespace felitronics::convolution

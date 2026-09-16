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
// COLD START: on the very first activation the shared history is empty — but a cold-started overlap-save
// FDL already produces the EXACT causal convolution (empty history = "no input before t=0"), so there is
// nothing to mask. Every swap, first or later, uses the same short (anti-click) crossfade.
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
    // the anti-click crossfade length, used for every swap (first activation and warm swaps alike). Message thread.
    bool prepare (int partitionSize, int maxIrSamples, int crossfadeSamples, int numChannels = 1)
    {
        prepared_ = false;                            // any early return below leaves the engine unprepared (setIr/process reject)
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

        inputSpec_.assign ((std::size_t) specF_, 0.0f);
        acc_.assign       ((std::size_t) specF_, 0.0f);
        ifftOut_.assign   ((std::size_t) N_,     0.0f);
        for (int c = 0; c < MaxChannels; ++c) chan_[c].prepare (P_, N_, specF_, maxParts_);

        cur_ = 0; xfadePos_ = 0; xfadeLen_ = warmXfade_;
        phase_ = 0; fdlPos_ = 0;
        state_.store (0, std::memory_order_relaxed);
        ranNc_ = 0;                                   // a re-prepare may NARROW channels_; the ledger of
                                                      // "what ran last call" belongs to the old topology
        prepared_ = true;                             // fully built — only now may setIr()/process() run
        return true;
    }

    // 🔴 THE AUDIO THE CALLER FED, AND NOTHING ELSE — the running history, and no decision about which
    // IR is live. Nothing below is touched by setIr(), the one message-thread verb this can run beside
    // (buildIr() no longer writes a cached tail — see there), so it races nothing and cancels no swap: a
    // fade in flight goes on from where it was, now over an empty history.
    // BOTH slots' cached tails go, not just the live one: while a fade runs, the other slot's tail is
    // blended into every output sample, so leaving it would replay pre-restart audio through the incoming
    // operator for up to P samples.
    // WHAT THIS VERB DOES NOT GIVE is independence from WHERE in a fade it was called — the fade it leaves
    // running is state the next programme hears. A restart that owes law 11a calls reset(), which ends the
    // fade. Measured by this class's own suite, a clear 128 and a clear 256 samples into a 512-sample fade
    // answer the same programme differently in 383 of 1600 samples (worst 1.6556e-01), where a reset()
    // with the swap in flight answers exactly as a reset() after it settled (the independence group).
    void clearAudioState() noexcept
    {
        for (int c = 0; c < channels_; ++c) chan_[c].reset();
        ranNc_ = 0;                               // nothing has run, so nothing can be stopping
        phase_ = 0; fdlPos_ = 0;
    }

    // Clears the running history (so the next swap fades in from silence) but KEEPS THE IR THE CALLER LAST
    // PUBLISHED — reset() means "flush the tail", not "revert the EQ", and THAT PROMISE COVERS A PUBLICATION
    // THE AUDIO THREAD HAS NOT PICKED UP YET: a restart never loses an accepted publication (law 11e;
    // prepare() is a re-initialisation and discards everything, by contract). So a swap in flight is ENDED
    // IN FAVOUR OF THE NEW OPERATOR, not abandoned — `cur_` flips only at fade end, so keeping it would put
    // the PREVIOUS IR back, and nothing would ever re-stage the published one (setIr() already said true,
    // and a consumer's retry flag is clear after a successful publish). On the body this replaces, the
    // suite's independence group reads the OLD operator in every cell {Pending, Crossfading} x
    // {mono, stereo} x both slot parities, a worst sample 7.494e-01 away from the settled restart.
    //
    // WHAT MAKES ADOPTION RIGHT IS LAW 11a INDEPENDENCE, NOT CLICK-FREEDOM: a restart must answer the next
    // programme with bits that do not depend on what came before, and a half-finished fade is precisely
    // such a dependency — two engines holding the same published operator, one mid-fade and one settled,
    // would answer differently for up to `xfadeLen_` samples.
    // THE PRICE IS AT THE SEAM AND IT IS NOT ZERO. Flushing the history is itself a cut — the previous
    // stream's tail stops mid-decay — so the first sample after reset() steps whichever slot is live: on
    // DC 0.5 into a settled 700-tap operator, 1.8203e-01 with nothing in flight. Adopting puts the new head
    // tap where the old one was at that sample, so it moves the step by AT MOST |Δh0|·|x|, in EITHER
    // direction: 7.3203e-01 for a pair whose head tap flips +0.70 -> -0.40, 1.3933e-01 for a +1 dB
    // broadband move and 2.2018e-01 for a -1 dB one, and the flush exactly for a shape change that leaves
    // the head tap alone. That difference IS the change the caller asked for; refusing to make it does not
    // avoid a step, it plays the wrong operator and never plays the right one. The house precedent is
    // `eq::EqBand::reset()`, which SNAPS a pending design onto the target rather than dropping it or playing
    // out its ramp.
    // Audio thread (or externally synced); must not run concurrently with setIr().
    void reset() noexcept
    {
        clearAudioState();
        // The acquire pairs with setIr()'s release store, so a slot adopted at Pending is fully built. (The
        // 1 -> 2 store in process() is relaxed: a reset() on the audio thread makes this load on that same
        // thread, and one that is "externally synced" gets its ordering from that synchronisation.)
        const int s = state_.load (std::memory_order_acquire);
        if (s == 0) return;                 // nothing in flight — and NOTHING IS STORED: an unconditional
                                            // store(0) would wipe a publication landing between the load
                                            // and the store, the same loss through a race instead of a
                                            // sequence. At Idle the flush above is all the work there is.
        cur_ = 1 - cur_;                    // the published slot becomes the live one
        xfadePos_ = 0;                      // hygiene: process() re-arms the clock at the next fade start
        // RELEASE, not relaxed: setIr() reads `cur_` after acquiring this store, and `cur_` is written just
        // above. Relaxed would let the message thread see Idle with a stale `cur_` and build the next IR
        // into the slot that is now LIVE, overwriting coefficients the audio thread is reading.
        // With the store made only here and no cached tail written off the audio thread, reset() no longer
        // races a single-producer loader at all: ThreadSanitizer, a loader publishing in a loop beside
        // process() / reset() / clearAudioState(), reports no race on this class — nor on its two swap-safe
        // siblings once they carry the same rule — against 54–62 across the three before P88 (Apple clang
        // with -fno-builtin — without it a std::fill write is not instrumented and the run is blind to the
        // tail write — and gcc 14 alike). That is a measurement, not yet the contract above, which stays until
        // a sanitizer row carries it. The store's CONDITION is gated by a real-thread test (law 11e).
        state_.store (0, std::memory_order_release);
    }

    static constexpr int latencySamples() noexcept { return 0; }
    bool isBusy() const noexcept { return state_.load (std::memory_order_acquire) != 0; }
    int  numChannels() const noexcept { return channels_; }

    // Channels that advanced state on the previous accepted call (see the note in process()).
    int  ranNc_ = 0;

    // Message thread — SINGLE producer; must not run concurrently with reset(). Per-channel IR (a
    // true-stereo IR); if fewer channels are supplied than configured, the last is broadcast to the
    // remainder. Stages into the INACTIVE slot (keeps the shared warm history). Returns false if a swap
    // is already pending/crossfading (the caller coalesces with the latest snapshot).
    bool setIr (const float* const* irPerCh, int nch, int len)
    {
        if (nch < 1 || ! prepared_ || state_.load (std::memory_order_acquire) != 0) return false;   // unprepared / busy
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
    // Law 11 (DSP-ARCHITECTURE.md §2).
    [[nodiscard]] bool process (const float* const* in, float* const* out, int numChannelsToProcess, int n) noexcept
    {
        if (numChannelsToProcess < 0 || n < 0) return false;
        if (! prepared_) return false;                              // never prepared — channel buffers are empty
        if (numChannelsToProcess > channels_) return false;         // width is a LIMIT — law 11(b)
        if (n == 0) return true;                                    // no samples: no time, no edge
        const int nc = numChannelsToProcess;

        // A channel that stops being asked for keeps a full frame, FDL and pending tail — frozen, not
        // decayed — and replays them when it is asked for again. Measured through CabConvolver, stereo ->
        // mono -> stereo: 9.12e-02 (-20.8 dBFS) on the returning channel out of DIGITAL SILENCE. Per
        // channel only: `phase_`, `fdlPos_` and `xfadePos_` below are the engine's shared block position
        // and are supposed to keep running for the channels that stayed.
        for (int c = nc; c < ranNc_; ++c) chan_[c].reset();
        ranNc_ = nc;
        if (nc == 0) return true;                                   // law 11(d): the edge above IS the work

        int s = state_.load (std::memory_order_acquire);
        if (s == 1)   // begin the crossfade. A cold-started FDL already yields the EXACT causal convolution, so ONE
        {             // short smoothstep fade (warmXfade_) masks only the silence→convolution onset — no long cold prime.
            xfadePos_ = 0;
            xfadeLen_ = warmXfade_;
            for (int c = 0; c < nc; ++c) primeTail (chan_[c], 1 - cur_);   // make the new slot's tail valid NOW from the FDL (else its zeroed tail leaks into the blend for ≤P samples)
            state_.store (2, std::memory_order_relaxed);
            s = 2;
        }

        if (s != 2) { processRange (in, out, nc, 0, n); return true; }   // Idle → single active slot
        processFade (in, out, nc, n);
        return true;
    }

    // Mono convenience (keeps a single-channel consumer on the 3-arg process()).
    [[nodiscard]] bool process (const float* in, float* out, int n) noexcept
    {
        const float* ins[1]  { in };
        float*       outs[1] { out };
        return process (ins, outs, 1, n);
    }

private:
    //==========================================================================
    // Per-channel: ONE running input history shared by TWO IR coefficient slots.
    struct Chan
    {
        core::fft::AlignedVector<float> frame;       // 2P: [prev | current] input accumulator (seam: forward input)
        core::fft::AlignedVector<float> fdl;         // maxParts × specF: ring of past INPUT spectra (seam)
        std::vector<float> pendingTail[2];   // P each: cached tail output per slot (time domain — plain)
        std::vector<float> h0[2];            // P each: direct (head) taps per slot (time domain — plain)
        core::fft::AlignedVector<float> irSpec[2];   // maxParts × specF each: tail partition spectra (seam)
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

        // Message thread: build slot k's head + tail spectra (allocates a scratch). It does NOT touch that
        // slot's cached tail, which is audio-thread state: primeTail() overwrites it whole at fade start
        // for every channel being processed, and a channel that is NOT being processed has had it zeroed
        // by the drop-out path in process() (or by prepare()) and is cold, for which zero is the right
        // value. Zeroing it here instead put a MESSAGE-thread write into a history buffer that the audio
        // thread also writes — the same 0.0f, so nothing could be lost, but a data race all the same, and
        // one that made clearAudioState()'s promise untrue: with this write put back, ThreadSanitizer reports
        // it against the history clear under a loader stress (2 reports a run, Apple clang -fno-builtin and
        // gcc 14 alike). Measured byte-identical against the body before the change over width transitions
        // x {no reset, reset at Idle} x block sizes.
        void buildIr (int k, const float* ir, int len, int P, int maxParts, int specF, Fft& fft)
        {
            if (P <= 0) return;                 // unprepared engine — nothing to build (also guards the /P below)
            if (len < 0) len = 0;
            for (int i = 0; i < P; ++i) h0[k][(std::size_t) i] = (i < len) ? ir[i] : 0.0f;

            const int tailLen = (len > P) ? (len - P) : 0;
            int parts = (tailLen > 0) ? ((tailLen + P - 1) / P) : 0;
            if (parts > maxParts) parts = maxParts;

            core::fft::AlignedVector<float> part ((std::size_t) (2 * P), 0.0f);   // forward input — aligned
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
    Chan chan_[(std::size_t) MaxChannels];               // size_t cast: GCC -Wsign-conversion flags a dependent int bound
    core::fft::AlignedVector<float> inputSpec_, acc_, ifftOut_;   // shared FFT scratch — SIMD-aligned (seam)
    std::atomic<int> state_ { 0 };                        // 0 Idle · 1 Pending (staged) · 2 Crossfading
    bool prepared_ = false;                               // true ONLY after a fully-successful prepare(); setIr()/process() reject until then (a partial/failed prepare leaves it false even if P_>0)
    int P_ = 0, N_ = 0, specF_ = 0, maxParts_ = 0, channels_ = 1;
    int cur_ = 0;                                         // active slot (0/1)
    int phase_ = 0, fdlPos_ = 0;                          // shared per-chunk timing (all channels lockstep)
    int xfadePos_ = 0, xfadeLen_ = 1, warmXfade_ = 1;
};

} // namespace felitronics::convolution

// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <algorithm>   // std::copy, std::min
#include <atomic>
#include <cstddef>
#include <cstdint>

namespace felitronics::analysis
{

//==============================================================================
// felitronics::analysis::RollingSpectrumTap — a lock-free SPSC analyser tap whose snapshot CADENCE
// (hop) is decoupled from its analysis WINDOW SIZE (N). The audio thread push()es samples into a
// rolling ring of the maximum window (kMaxSize = 1<<MaxOrder); on demand it publishIfDue()s the most
// recent N = (1<<order) samples into a single-slot immutable mailbox the GUI drains with tryPull().
//
// Why not the plain SpectrumTap? SpectrumTapT publishes exactly one window every kSize samples — hop
// is welded to the window, so a bigger window means a slower frame rate (an 8192-pt window at 48 kHz
// would refresh at ~5.8 fps at EVERY resolution). Here the consumer chooses order and hop per block:
// ONE ring serves any order ≤ MaxOrder at a steady UI-rate cadence (overlapping windows for large N,
// gapped windows for small N — both correct for a spectrum display), with a single write per sample
// and no per-order ring duplication.
//
// Threading contract — tear-free by the same ownership handoff SpectrumTap uses, extended with the
// window metadata:
//   • The ring (`ring`, `wpos`, `lastPublish`, `lastOrder`) is touched ONLY by the audio thread; the
//     GUI never reads it. push() and publishIfDue() are both audio-thread — the snapshot copy reads
//     the ring on the same thread that writes it, so it can never tear against a concurrent push().
//   • The mailbox (`data`, `frameSize`, `frameOrder`) is written by the audio thread ONLY while
//     `ready == false`, and read by the GUI ONLY while `ready == true`. The acquire/release `ready`
//     flag hands off exclusive ownership of the mailbox — there is never concurrent access to it.
//   • The request (`order`, `hopSamples`) is passed BY VALUE into publishIfDue from a single atomic
//     the consumer reads ONCE per block, so the audio thread never splits a resolution change across
//     the tap-select and the window-size decision (no torn request).
//
// RT-safe: fixed-size storage, bounded copies (≤ kMaxSize floats), no alloc / lock / IO / throw.
// JUCE-free (`std::copy`). Cache-line separation (alignas 64) keeps the producer's per-sample ring
// write off the cache line the GUI polls for `ready`/`data`.
//
// A resolution switch = the consumer changing the `order` it passes. publishIfDue force-publishes a
// fresh frame the moment `order` differs from the last published one (it does not wait for the next
// hop), and the GUI discards any frame whose reported order ≠ the one it now wants — so a live switch
// is click-free: never a blank frame, never a wrong-size (garbage) frame. The ring is kept warm by
// the continuous feed, so the first frame at the new order is already a full, valid window.
template <int MaxOrder = 13>
struct RollingSpectrumTapT
{
    static constexpr int kMaxOrder = MaxOrder;
    static constexpr int kMaxSize  = 1 << MaxOrder;   // power of two → mask wrap
    static_assert (MaxOrder >= 1 && MaxOrder <= 30, "RollingSpectrumTap: order out of range");

    //==============================================================================
    // Audio thread.

    // Append one sample to the rolling history. One store + a monotonic counter bump; no branch.
    void push (float s) noexcept
    {
        ring[(std::size_t) (wpos & (std::uint64_t) (kMaxSize - 1))] = s;
        ++wpos;
    }

    // Publish the most-recent (1<<order) samples into the mailbox IF a hop has elapsed OR `order`
    // changed since the last published frame (a resolution switch forces an immediate fresh frame).
    // No-op until the ring holds ≥ that many samples (avoids a short/garbage first frame after reset),
    // and a no-op if the GUI hasn't drained the previous frame (the slot stays owned by the reader —
    // the forced state persists, so the next block retries). `hopSamples` is the consumer's cadence
    // target in samples (e.g. round(fs / uiRateHz)); it is clamped to ≥ 1.
    void publishIfDue (int order, int hopSamples) noexcept
    {
        if (order < 0)          order = 0;
        if (order > kMaxOrder)  order = kMaxOrder;

        const std::uint64_t N = (std::uint64_t) 1 << order;
        if (wpos < N) return;                                        // not enough history yet

        const bool          forced = (order != lastOrder);           // resolution switch → publish now
        const std::uint64_t hop    = (std::uint64_t) (hopSamples > 1 ? hopSamples : 1);
        if (! forced && (wpos - lastPublish) < hop) return;          // hop not elapsed

        if (ready.load (std::memory_order_acquire)) return;          // reader still owns the slot

        // Copy the most-recent N samples in chronological order, splitting on the ring wrap.
        const std::uint64_t start = wpos - N;                        // oldest sample of the window
        const std::size_t   mask  = (std::size_t) (kMaxSize - 1);
        const std::size_t   s0    = (std::size_t) (start & (std::uint64_t) mask);
        const std::size_t   nN    = (std::size_t) N;
        const std::size_t   first = std::min (nN, (std::size_t) kMaxSize - s0);
        std::copy (ring + s0, ring + s0 + first, data);
        if (nN > first) std::copy (ring, ring + (nN - first), data + first);

        frameSize   = (int) N;
        frameOrder  = order;
        lastPublish = wpos;
        lastOrder   = order;
        ready.store (true, std::memory_order_release);
    }

    // Restart the PRODUCER only. Deliberately does NOT touch the mailbox (`ready`/`data`/frameSize/
    // frameOrder): reset() runs on the audio (prepareToPlay) thread while tryPull() may be mid-copy on
    // the message thread, so clearing `ready` here would revoke the reader's ownership and let the next
    // publishIfDue overwrite `data` underneath it — a torn frame / data race. Leaving the handshake to
    // the reader keeps reset() race-free: any in-flight frame is drained once by the reader, then the
    // restarted ring (lastOrder = -1 re-arms the forced first frame) publishes fresh once warm again.
    void reset() noexcept
    {
        wpos = 0;
        lastPublish = 0;
        lastOrder = -1;                                              // force the first frame after restart
    }

    //==============================================================================
    // GUI thread.

    // If a frame is ready, copy its `frameSize` samples into dst (dst must hold ≥ kMaxSize floats),
    // report the order it was captured at via outOrder, re-arm the slot, and return true. Else false.
    // The caller compares outOrder against the resolution it currently wants and discards a mismatch.
    bool tryPull (float* dst, int& outOrder) noexcept
    {
        if (! ready.load (std::memory_order_acquire)) return false;
        const int n = frameSize;                                     // published before the release store
        outOrder = frameOrder;
        std::copy (data, data + n, dst);
        ready.store (false, std::memory_order_release);
        return true;
    }

    //==============================================================================
    // Producer-only rolling history (audio thread) — its own cache line, away from the GUI mailbox.
    alignas (64) float         ring[(std::size_t) kMaxSize] {};
    std::uint64_t              wpos        = 0;   // total samples written (monotonic); index = wpos & (kMaxSize-1)
    std::uint64_t              lastPublish = 0;   // wpos at the last successful publish
    int                        lastOrder   = -1;  // order of the last publish (−1 → force the first frame)

    // Cross-thread mailbox — its own cache line. Written by audio only while !ready, read by GUI only
    // while ready; the flag is the ownership handoff.
    alignas (64) std::atomic<bool> ready { false };
    alignas (64) float             data[(std::size_t) kMaxSize] {};
    int                            frameSize  = 0;   // N of the published frame (1<<frameOrder)
    int                            frameOrder = 0;   // order the frame was captured at
};

// Drop-in default: a 16384-sample (order-14) rolling history — enough for any analyser window
// 1024 / 2048 / 4096 / 8192 / 16384 without a fork.
using RollingSpectrumTap = RollingSpectrumTapT<14>;

} // namespace felitronics::analysis

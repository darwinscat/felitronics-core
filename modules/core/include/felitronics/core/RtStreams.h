// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

//==============================================================================
// felitronics::core — the RT buffer-swap discipline, made a type. PROMOTED from OrbitCapture
// (audio/RtStreams.h, de-monolith step 8) when a second capture product needed the same plumbing,
// then hardened behind an adversarial review (codex + deepseek): zero-length publishes no longer
// arm an empty looping buffer, the recorder's length is release-published so a harvest never reads
// an unsynchronized sample, the audition buffer is FIXED-SIZE (no reliance on assign-not-
// reallocating), and int lengths are guarded at reserve time.
//
// The audio thread streams these buffers while a mode/playing atomic is up. The load-bearing rule:
//
//     flag DOWN  →  mutate the buffer WITHIN its fixed/reserved storage  →  flag UP
//
// SINGLE-WRITER contract (all three types): exactly ONE message-thread writer and ONE RT reader
// (RecStream: one RT writer, one message-thread harvester). These are not general SPSC queues —
// they are the family's capture-app swap discipline with its exact semantics.
//
// READER RULES (the product's audio callback):
//   • Gate every read on the flag with memory_order_acquire — a relaxed gate can observe the flag
//     up while missing the writer's len/pos resets (a real head-drop on ARM, found by the crew).
//   • Re-load the flag each BLOCK (not each sample); never cache buf.data() across blocks.
//
// KNOWN, DELIBERATE GAP (crew consensus, fix deferred): a publish()/refill() that lands while an
// already-admitted callback block is still streaming is a formal data race on the sample memory —
// the flag went down after the reader's block-start check. In practice the window is one audio
// block (~10 ms) hit only by user-action-rate publishes, and every float is written/read whole;
// but it is UB by the letter, so treat it as a bound to engineer away, not a licence: the planned
// fix is a reader-acknowledged epoch (the reader bumps a counter when it observes the flag down;
// the writer waits for the bump before mutating). That needs reader cooperation and lands as its
// own change; until then, writers SHOULD leave a block-length pause between flag-down and mutate
// where the product allows it.
//==============================================================================

#include <atomic>
#include <cassert>
#include <cstddef>
#include <climits>
#include <cstring>
#include <vector>

namespace felitronics::core {

// A pre-rendered audition clip: the message thread publishes a whole render, the RT callback
// streams buf[0..len) (loop-aware). The buffer is FIXED-SIZE from reserve() on — publish() copies
// into it and never touches the vector's shape, so the RT pointer is valid by construction (no
// reliance on assign()-within-capacity implementation behaviour). The reader must use `len`,
// never buf.size().
struct AuditionStream {
    std::vector<float> buf;
    std::atomic<bool> playing { false }, loop { false };
    std::atomic<int>  pos { 0 }, len { 0 };

    void reserve(std::size_t cap) {
        assert(cap <= (std::size_t)INT_MAX);                // len/pos are int (crew: overflow guard)
        buf.assign(cap, 0.0f);                              // fixed storage; size() stays == cap
    }
    void stop() { playing.store(false, std::memory_order_release); }

    // stop → copy within the fixed storage (truncating) → start. A ZERO-length render never arms:
    // the reader's loop path resets pos to 0 and reads buf[0] — on an armed empty clip that is an
    // out-of-bounds read (crew finding, pinned in the tests).
    void publish(const float* data, std::size_t n) {
        playing.store(false, std::memory_order_release);
        const std::size_t nn = n < buf.size() ? n : buf.size();
        if (nn > 0) std::memcpy(buf.data(), data, nn * sizeof(float));
        len.store((int)nn); pos.store(0);
        if (nn > 0) playing.store(true, std::memory_order_release);
    }
};

// The live-convolver input stream. mode: 0 off · 1 DI clip through the IR · 2 live input
// through the IR. `buf`/`pos` are only read in mode 1; the mode atomic also gates the
// convolver + live-input path, so it stays the single flag for the whole conv route.
// NB the reader takes the clip LENGTH from buf.size(), so this type keeps assign() semantics —
// the no-reallocation-within-capacity reliance is implementation behaviour (universal in
// libc++/libstdc++/MSVC), documented here rather than assumed silently.
struct ConvStream {
    std::vector<float> buf;
    std::atomic<int> mode { 0 };
    std::atomic<int> pos { 0 };

    void reserve(std::size_t cap) {
        assert(cap <= (std::size_t)INT_MAX);                // pos is int (crew: overflow guard)
        buf.reserve(cap);
    }
    void setMode(int m) { mode.store(m, std::memory_order_release); }
    void stop() { setMode(0); }

    // Refill the DI clip. The caller must have taken the mode DOWN first (stop()) — the RT
    // side may otherwise stream the buffer mid-assign. Truncates to the reserved capacity.
    void refill(const float* data, std::size_t n) {
        assert(mode.load(std::memory_order_acquire) == 0);  // never mutate under a live reader
        const std::size_t nn = n < buf.capacity() ? n : buf.capacity();
        buf.assign(data, data + (std::ptrdiff_t)nn);
        pos.store(0);
    }
};

// The live-input recorder: the RT callback appends dry input samples while `recording` is up;
// the message thread harvests buf[0..len) only after taking it down. `reserve` RESIZES (the RT
// side writes by index), so it must run while the device is stopped (audioDeviceAboutToStart).
//
// Harvest pattern (message thread):  stop();  n = len.load();  read buf[0..n).
// push() release-stores len, so the harvester's load (seq_cst ⊇ acquire) synchronizes with the
// last push and every harvested sample is visible — a relaxed len would let len=i+1 become
// visible before buf[i] on ARM (crew finding). The RT gate on `recording` must be ACQUIRE: it
// pairs with start()'s release so the len=0 rewind is seen before the first push (a relaxed gate
// can drop the head of a take by seeing its own stale full length).
struct RecStream {
    std::vector<float> buf;
    std::atomic<bool> recording { false };
    std::atomic<int>  len { 0 };

    void reserve(std::size_t cap) {
        assert(cap <= (std::size_t)INT_MAX);                // len is int (crew: overflow guard)
        buf.assign(cap, 0.0f);
    }
    void start() { len.store(0); recording.store(true, std::memory_order_release); }
    void stop()  { recording.store(false, std::memory_order_release); }
    bool full() const { return (std::size_t)len.load() >= buf.size(); }

    inline void push(float s) {                          // RT side: append while space remains
        const int i = len.load(std::memory_order_relaxed);
        if ((std::size_t)i < buf.size()) {
            buf[(std::size_t)i] = s;
            len.store(i + 1, std::memory_order_release); // publish the sample WITH the new length
        }
    }
};

} // namespace felitronics::core

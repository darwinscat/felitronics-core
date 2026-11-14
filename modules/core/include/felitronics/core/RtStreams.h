// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

//==============================================================================
// felitronics::core — the RT buffer-swap discipline, made a type. PROMOTED from OrbitCapture
// (audio/RtStreams.h) when a second capture product needed the same plumbing. Three streams share one
// guard: AuditionStream (a pre-rendered clip the audio thread plays), ConvStream (the clip or the live
// input fed to a convolver) and RecStream (the audio thread records, the controller harvests).
//
// THE GATE. Each stream is guarded by one atomic word, `RtGate`: the MODE in the low 31 bits (0 = down)
// and a BUSY bit the audio thread holds while it is inside a block it was admitted to. The word changes
// ONLY by read-modify-write, so it has one modification order and every change continues the release
// sequence of the change before it.
//
//   audio thread   a Block opens:   admitted only by a CAS that sees the mode up and sets BUSY (acquire)
//                  a Block closes:  stores what it advanced (position, length), then clears BUSY, and the
//                                   mode too when the clip ended, in one release RMW
//   controller     try…():          takes the mode down and learns, in that same RMW, whether a block is
//                                   open. Open: refuse and mutate nothing. Closed: mutate, then raise the
//                                   mode (release) if the call starts something.
//
// So no sample the controller writes is in use by an admitted block, no store a block makes at its close
// lands after a re-arm, and a block admitted after the re-arm sees the whole new content. There is no
// counter and no wait. A stopped device holds no block, so every try lands at once. A try refused because
// a block was open lands on a retry after that block closes: nothing can open a block while the mode is
// down.
//
// THREADS. ONE controller thread calls every member that is not a Block's, tryReserve() included. ONE
// audio thread opens Blocks. isPlaying(), isRecording(), mode(), recorded() and `loop` may be read from
// any thread. A device's lifecycle callbacks (JUCE's audioDeviceAboutToStart / audioDeviceStopped) are not
// the controller: JUCE runs them on the thread that adds or removes the callback, or wherever the backend
// restarts the device, and documents no thread, so a call made there can overlap the controller's. Size the
// streams from the controller.
//
// THE AUDIO THREAD'S RULES.
//   • Open each stream's Block once per callback, read and write the stream only through it, and let it
//     close (leave its scope) before the callback returns. Never keep its pointer past it.
//   • Keep a Block's scope to the loop that needs it. BUSY lasts as long as the scope, and every try made
//     meanwhile is refused and must be retried.
//
// THE CONTROLLER'S RULES.
//   • Every try…() is [[nodiscard]]. A refusal mutated nothing and left the stream stopped; keep what
//     you meant to publish and try again from a timer (the house's ConvolutionEngine::setIr + isBusy()
//     pattern).
//
// FAILURE IS VISIBLE, NEVER A RACE. A Block that never closes (a thread killed inside a callback) keeps
// BUSY set: every later try that needs the stream to itself is refused (stop() and trySetMode(0) still take
// the mode down) and every later Block is empty. Silence, not torn audio.
//==============================================================================

#include <atomic>
#include <climits>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <vector>

namespace felitronics::core {

static_assert (std::atomic<std::uint32_t>::is_always_lock_free, "RtGate needs a lock-free 32-bit atomic");
static_assert (std::atomic<int>::is_always_lock_free && std::atomic<bool>::is_always_lock_free,
               "the streams' readings need lock-free atomics");

// The guard itself, public so a product can gate its own swaps of plain memory the same way. `Atom` is the
// test seam: any type with std::atomic<std::uint32_t>'s load / compare_exchange_strong / fetch_and (with
// explicit orders) and a constructor from std::uint32_t. Products use `RtGate`.
//
// The controller's sequence is   if (gate.quiesce()) { mutate; (void) gate.arm (mode); }   and a TRUE
// quiesce() before the mutation is the caller's precondition: it is what orders the last block's accesses
// before the controller's. arm() lands only on a word that is exactly 0 (down, no block open), so it refuses
// an arm over a raised mode or an open block; it cannot tell whether the quiesce() before it was true.
template <class Atom = std::atomic<std::uint32_t>>
class BasicRtGate {
public:
    static constexpr std::uint32_t kBusy     = 0x80000000u;
    static constexpr std::uint32_t kModeMask = 0x7fffffffu;

    // AUDIO THREAD. Admits a block: returns the mode it was admitted at, or 0 (not admitted — the mode is
    // down, or a block of this gate is already open). Lock-free, not wait-free: the CAS retries only when the
    // controller changed the word in between, which a product does at user rate.
    [[nodiscard]] std::uint32_t enter() noexcept {
        std::uint32_t v = s_.load (std::memory_order_acquire);
        while ((v & kModeMask) != 0u && (v & kBusy) == 0u)
            if (s_.compare_exchange_strong (v, v | kBusy, std::memory_order_acq_rel, std::memory_order_acquire))
                return v & kModeMask;
        return 0u;
    }
    // AUDIO THREAD. Closes an admitted block: clears BUSY, and the mode too when `down` (the clip ended).
    // Everything the block stored before this call is published by it.
    void leave (bool down) noexcept {
        (void) s_.fetch_and (down ? 0u : kModeMask, std::memory_order_release);
    }

    // CONTROLLER. Takes the mode down. True when no block is open: from here until arm(), the guarded
    // memory belongs to the controller alone. False when a block is open: mutate nothing, retry later.
    [[nodiscard]] bool quiesce() noexcept {
        return (s_.fetch_and (kBusy, std::memory_order_acq_rel) & kBusy) == 0u;
    }
    // CONTROLLER. Raises the mode (low 31 bits; 0 raises nothing). Lands only on a word that is exactly 0.
    [[nodiscard]] bool arm (std::uint32_t mode) noexcept {
        std::uint32_t expected = 0u;
        const std::uint32_t m = mode & kModeMask;
        return m != 0u && s_.compare_exchange_strong (expected, m, std::memory_order_release,
                                                      std::memory_order_relaxed);
    }

    // ANY THREAD: readings.
    [[nodiscard]] std::uint32_t mode() const noexcept { return s_.load (std::memory_order_acquire) & kModeMask; }
    [[nodiscard]] bool busy() const noexcept { return (s_.load (std::memory_order_acquire) & kBusy) != 0u; }

private:
    Atom s_ { 0u };
};
using RtGate = BasicRtGate<>;

namespace rtstreams_detail {
    // Capacities are int on the audio thread's side; a larger request is clamped, never refused, so a false
    // from a try always means "a block was open".
    inline std::size_t clampCapacity (std::size_t cap) noexcept {
        return cap < (std::size_t) INT_MAX ? cap : (std::size_t) INT_MAX;
    }
    inline int clampPos (int pos, int len) noexcept { return pos < 0 ? 0 : (pos > len ? len : pos); }
}

// A pre-rendered audition clip. The controller publishes a whole render into FIXED storage (sized by
// tryReserve, never reshaped by a publish); the audio thread streams data()[pos, len()) through a Block,
// looping when `loop` is set, and calls finish() when a one-shot clip has ended.
template <class Atom = std::atomic<std::uint32_t>>
class BasicAuditionStream {
public:
    std::atomic<bool> loop { false };                     // any thread; the audio thread reads it per block

    // CONTROLLER. Sizes the storage (allocates), stops, empties. Refused while a block is open.
    [[nodiscard]] bool tryReserve (std::size_t capacity) {
        const std::size_t cap = rtstreams_detail::clampCapacity (capacity);
        if (! gate_.quiesce()) return false;
        len_ = 0; pos_ = 0;                                // first: an assign that throws leaves an EMPTY clip
        buf_.assign (cap, 0.0f);
        return true;
    }
    // CONTROLLER. Stops: the stream goes down at once; a block already open finishes its loop.
    void stop() noexcept { (void) gate_.quiesce(); }
    // CONTROLLER. Copies a render into the storage (truncating to it), rewinds and plays it. A ZERO-length
    // render (or a null pointer) leaves the stream stopped and empty: an armed empty clip would let a loop
    // path read data()[0] out of bounds. Refused, with nothing changed but the stop, while a block is open.
    [[nodiscard]] bool tryPublish (const float* data, std::size_t n) noexcept {
        if (! gate_.quiesce()) return false;
        const std::size_t nn = data == nullptr ? 0 : (n < buf_.size() ? n : buf_.size());
        if (nn > 0) std::memcpy (buf_.data(), data, nn * sizeof (float));
        len_ = (int) nn; pos_ = 0;
        if (nn > 0) (void) gate_.arm (1u);
        return true;
    }
    [[nodiscard]] bool isPlaying() const noexcept { return gate_.mode() != 0u; }  // any thread
    [[nodiscard]] std::size_t capacity() const noexcept { return buf_.size(); }   // controller

    // AUDIO THREAD. Empty (false) unless the clip plays. While it is open, data()[0, len()) is the clip and
    // `pos` the read position to advance; on close the stream keeps pos, clamped to [0, len()].
    class Block {
    public:
        int pos = 0;
        explicit operator bool() const noexcept { return s_ != nullptr; }
        [[nodiscard]] const float* data() const noexcept { return data_; }
        [[nodiscard]] int len() const noexcept { return len_; }
        void finish() noexcept { ended_ = true; }         // the clip ended: the stream goes down on close
        Block (const Block&) = delete;
        Block& operator= (const Block&) = delete;
        ~Block() {
            if (s_ == nullptr) return;
            s_->pos_ = rtstreams_detail::clampPos (pos, len_);
            s_->gate_.leave (ended_);
        }
    private:
        friend class BasicAuditionStream;
        explicit Block (BasicAuditionStream* s) noexcept : s_ (s) {
            if (s_ == nullptr) return;
            pos = s_->pos_; len_ = s_->len_; data_ = s_->buf_.data();
        }
        BasicAuditionStream* s_;
        const float* data_ = nullptr;
        int len_ = 0;
        bool ended_ = false;
    };
    [[nodiscard]] Block beginBlock() noexcept { return Block (gate_.enter() != 0u ? this : nullptr); }

private:
    std::vector<float> buf_;
    int len_ = 0, pos_ = 0;                                // guarded by the gate, so plain
    BasicRtGate<Atom> gate_;
};

// The live-convolver input. mode: 0 off · kClip (1) the clip through the convolver · kLive (2) the live
// input through it — any positive mode is accepted, the meaning is the product's. data()/len()/pos are
// the clip's, used in kClip. The storage is FIXED (sized by tryReserve); len() is the clip, never the
// storage.
template <class Atom = std::atomic<std::uint32_t>>
class BasicConvStream {
public:
    static constexpr int kOff = 0, kClip = 1, kLive = 2;

    // CONTROLLER. Sizes the clip storage (allocates), stops, empties. Refused while a block is open.
    [[nodiscard]] bool tryReserve (std::size_t capacity) {
        const std::size_t cap = rtstreams_detail::clampCapacity (capacity);
        if (! gate_.quiesce()) return false;
        len_ = 0; pos_ = 0;                                // first: an assign that throws leaves an EMPTY clip
        buf_.assign (cap, 0.0f);
        return true;
    }
    // CONTROLLER. Stops, copies the clip (truncating to the storage), rewinds, and leaves the route OFF.
    // Raise it with trySetMode, which cannot be refused right after a true tryRefill: no block can open
    // while the mode is down. Refused while a block is open.
    [[nodiscard]] bool tryRefill (const float* data, std::size_t n) noexcept {
        if (! gate_.quiesce()) return false;
        const std::size_t nn = data == nullptr ? 0 : (n < buf_.size() ? n : buf_.size());
        if (nn > 0) std::memcpy (buf_.data(), data, nn * sizeof (float));
        len_ = (int) nn; pos_ = 0;
        return true;
    }
    // CONTROLLER. A positive mode switches the route through a quiesce, so a block that ends the clip cannot
    // take the new mode down behind it; the clip position is kept. A mode <= 0 stops, and is never refused.
    [[nodiscard]] bool trySetMode (int m) noexcept {
        if (m <= 0) { stop(); return true; }
        if (! gate_.quiesce()) return false;
        (void) gate_.arm ((std::uint32_t) m);
        return true;
    }
    void stop() noexcept { (void) gate_.quiesce(); }                               // controller
    [[nodiscard]] int mode() const noexcept { return (int) gate_.mode(); }       // any thread
    [[nodiscard]] std::size_t capacity() const noexcept { return buf_.size(); }   // controller

    // AUDIO THREAD. Empty unless the route is up; mode() is the mode it was admitted at. Close it before the
    // convolver runs: BUSY should span the copy loop only.
    class Block {
    public:
        int pos = 0;
        explicit operator bool() const noexcept { return s_ != nullptr; }
        [[nodiscard]] int mode() const noexcept { return mode_; }
        [[nodiscard]] const float* data() const noexcept { return data_; }
        [[nodiscard]] int len() const noexcept { return len_; }
        void finish() noexcept { ended_ = true; }         // the route ends: the stream goes down on close
        Block (const Block&) = delete;
        Block& operator= (const Block&) = delete;
        ~Block() {
            if (s_ == nullptr) return;
            s_->pos_ = rtstreams_detail::clampPos (pos, len_);
            s_->gate_.leave (ended_);
        }
    private:
        friend class BasicConvStream;
        Block (BasicConvStream* s, int mode) noexcept : s_ (s), mode_ (mode) {
            if (s_ == nullptr) return;
            pos = s_->pos_; len_ = s_->len_; data_ = s_->buf_.data();
        }
        BasicConvStream* s_;
        int mode_ = 0;
        const float* data_ = nullptr;
        int len_ = 0;
        bool ended_ = false;
    };
    [[nodiscard]] Block beginBlock() noexcept {
        const int m = (int) gate_.enter();
        return Block (m != 0 ? this : nullptr, m);
    }

private:
    std::vector<float> buf_;
    int len_ = 0, pos_ = 0;                                // guarded by the gate, so plain
    BasicRtGate<Atom> gate_;
};

// The live-input recorder: the audio thread appends through a Block while recording is up; the controller
// harvests the take. tryHarvest stops and copies buf[0, n) only once no block is open, so n is FINAL; a
// tryStart rewinds only then, so a block of the old take can neither write into the new one nor undo its
// rewind.
template <class Atom = std::atomic<std::uint32_t>>
class BasicRecStream {
public:
    // CONTROLLER. Sizes the storage (allocates), stops, discards the take. Refused while a block is open.
    [[nodiscard]] bool tryReserve (std::size_t capacity) {
        const std::size_t cap = rtstreams_detail::clampCapacity (capacity);
        if (! gate_.quiesce()) return false;
        n_ = 0; shown_.store (0, std::memory_order_relaxed); // first: an assign that throws leaves an EMPTY take
        buf_.assign (cap, 0.0f);
        return true;
    }
    // CONTROLLER. Rewinds and records. Refused while a block is open.
    [[nodiscard]] bool tryStart() noexcept {
        if (! gate_.quiesce()) return false;
        n_ = 0; shown_.store (0, std::memory_order_relaxed);
        (void) gate_.arm (1u);
        return true;
    }
    void stop() noexcept { (void) gate_.quiesce(); }                               // controller
    // CONTROLLER (allocates `take`). Stops; once no block is open, copies the whole take into `take` and
    // returns true. Refused, with `take` untouched, while a block is open: that block may still append.
    [[nodiscard]] bool tryHarvest (std::vector<float>& take) {
        if (! gate_.quiesce()) return false;
        take.assign (buf_.begin(), buf_.begin() + n_);
        return true;
    }
    [[nodiscard]] bool isRecording() const noexcept { return gate_.mode() != 0u; } // any thread
    // ANY THREAD: the length recorded as of the last closed block, a reading for display only.
    [[nodiscard]] std::size_t recorded() const noexcept {
        return (std::size_t) shown_.load (std::memory_order_relaxed);
    }
    [[nodiscard]] std::size_t capacity() const noexcept { return buf_.size(); }   // controller
    [[nodiscard]] bool full() const noexcept { return recorded() >= buf_.size(); } // controller

    // AUDIO THREAD. Empty unless recording. push() appends while space remains.
    class Block {
    public:
        explicit operator bool() const noexcept { return s_ != nullptr; }
        void push (float v) noexcept { if (i_ < cap_) d_[i_++] = v; }
        [[nodiscard]] bool full() const noexcept { return i_ >= cap_; }
        Block (const Block&) = delete;
        Block& operator= (const Block&) = delete;
        ~Block() {
            if (s_ == nullptr) return;
            s_->n_ = i_;
            s_->shown_.store (i_, std::memory_order_relaxed);
            s_->gate_.leave (false);
        }
    private:
        friend class BasicRecStream;
        explicit Block (BasicRecStream* s) noexcept : s_ (s) {
            if (s_ == nullptr) return;
            i_ = s_->n_; cap_ = (int) s_->buf_.size(); d_ = s_->buf_.data();
        }
        BasicRecStream* s_;
        float* d_ = nullptr;
        int i_ = 0, cap_ = 0;
    };
    [[nodiscard]] Block beginBlock() noexcept { return Block (gate_.enter() != 0u ? this : nullptr); }

private:
    std::vector<float> buf_;
    int n_ = 0;                                            // guarded by the gate, so plain
    std::atomic<int> shown_ { 0 };                         // the display's copy of n_
    BasicRtGate<Atom> gate_;
};

using AuditionStream = BasicAuditionStream<>;
using ConvStream     = BasicConvStream<>;
using RecStream      = BasicRecStream<>;

} // namespace felitronics::core

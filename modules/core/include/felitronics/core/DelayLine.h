// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace felitronics::core
{

//==============================================================================
// felitronics::core::DelayLine — a fixed-size single-channel integer delay line (RT-safe ring). The
// generic building block for LOOKAHEAD (delay the signal by the lookahead while a detector reads ahead),
// PDC, and short delays. prepare() allocates; process()/reset are alloc-free. One per channel.
class DelayLine
{
public:
    // WHAT prepare() ASKS THE HEAP FOR (law 11d) — the one function it sizes itself with, so a caller
    // budgeting memory reads the number the ring is actually built from and the two cannot drift. Asked
    // of a FRESH line; one already holding at least this much keeps it and asks for nothing.
    struct Storage
    {
        std::size_t samples = 1;      // ring slots: the delay PLUS ONE, so a delay of `cap` is addressable
        std::uint64_t bytes() const noexcept { return (std::uint64_t) sizeof (float) * (std::uint64_t) samples; }
    };

    // Negative reads as zero, exactly as prepare() reads it: the budget of a call is storageFor() with the
    // SAME arguments.
    static Storage storageFor (int maxDelaySamples) noexcept
    {
        Storage s;
        s.samples = (std::size_t) (maxDelaySamples < 0 ? 0 : maxDelaySamples) + 1u;
        return s;
    }

    // What a DEFAULT-CONSTRUCTED line asks for — the 1-slot passthrough below. A bank that default-builds
    // n lines and then prepares them pays this n times in blocks that are FREED again, so a budget that
    // models such a bank has to carry n transients that are never live at the end. The sizing constructor
    // below is how a bank avoids owing them at all.
    static constexpr std::uint64_t constructBytes() noexcept { return (std::uint64_t) sizeof (float); }

    DelayLine() = default;

    // BUILT AT ITS FINAL SIZE, in ONE allocation — for the per-channel banks in Compressor, Saturator and
    // TruePeakLimiter, which used to `assign(n, DelayLine{})`: a temporary, n copies of it, and then n
    // reallocations at the real size. Same end state, three fewer allocation classes for a budget to model.
    explicit DelayLine (int maxDelaySamples)
        : buf_ (storageFor (maxDelaySamples).samples, 0.0f),
          cap_ (maxDelaySamples < 0 ? 0 : maxDelaySamples) {}

    void prepare (int maxDelaySamples)
    {
        cap_ = maxDelaySamples < 0 ? 0 : maxDelaySamples;
        buf_.assign (storageFor (maxDelaySamples).samples, 0.0f);   // +1 so a delay of `cap_` is addressable
        // THE DELAY IS RE-CLAMPED INTO THE NEW CAPACITY, and that repairs an invariant this call used to
        // break: `setDelay` holds `delay_` inside [0, cap_], but a prepare() to a SMALLER capacity left the
        // old delay standing, after which `process()` computes `pos_ - delay_` and adds the ring length
        // ONCE — still negative, and the read is out of bounds. `prepare(8); setDelay(8); prepare(2);` is
        // the whole reproduction, and it predates the delay bank below.
        if (delay_ > cap_) delay_ = cap_;
        reset();
    }

    void reset() noexcept { std::fill (buf_.begin(), buf_.end(), 0.0f); pos_ = 0; }

    void setDelay (int d) noexcept { delay_ = d < 0 ? 0 : (d > cap_ ? cap_ : d); }
    int  delay()    const noexcept { return delay_; }
    int  capacity() const noexcept { return cap_; }

    // Push one sample, return the sample `delay()` samples ago (delay 0 = passthrough). RT-safe.
    inline float process (float x) noexcept
    {
        const int n = (int) buf_.size();
        buf_[(std::size_t) pos_] = x;
        int rd = pos_ - delay_;
        if (rd < 0) rd += n;
        const float y = buf_[(std::size_t) rd];
        if (++pos_ >= n) pos_ = 0;
        return y;
    }

private:
    std::vector<float> buf_ { 0.0f };   // default-VALID: a 1-slot delay-0 passthrough until prepare() (no empty-buffer OOB, no hot-path branch)
    int cap_ = 0, delay_ = 0, pos_ = 0;
};

//==============================================================================
// A PER-CHANNEL BANK OF DELAY LINES, sized in place — the shape `dynamics::Compressor`,
// `saturation::Saturator` and `limiter::TruePeakLimiter` all need, in one place so their budgets read the
// same arithmetic their preparations run.
//
// It replaces `bank.assign (n, DelayLine {})`, which cost FOUR allocation classes where one will do: a
// temporary line, n copies of it, and then n reallocations when each copy was prepared at the real size —
// three of them transient, and every one of them something a budget had to model without any of it being
// live at the end. Here a fresh bank builds each line AT its final size and a bank that already has the
// right shape asks the heap for nothing at all, which is what makes re-preparing a chain at the same
// geometry cost zero.
//
// IT KEEPS MORE THAN THE `assign` IT REPLACES DID, AND EXACTLY WHERE IS WORTH BEING PRECISE ABOUT — the
// obvious statement is wrong. `assign(n, DelayLine {})` did NOT hand a narrower bank's memory back in
// general: it copy-assigns into the elements that survive, and a `vector` copy-assigned from a shorter one
// keeps its capacity, so the old form retained just as this one does. The two diverge when the bank GROWS
// PAST ITS CAPACITY: `assign` reallocated the outer vector and built fresh lines from the temporary,
// destroying the old lines and their rings with them, while `reserve` + `emplace_back` MOVES the existing
// lines into the new buffer and their rings travel with them. Measured on that edge:
// `Compressor::prepare(48000, 64, 1, 50 ms)` then `(48000, 64, 2, 1 ms)` holds 9 880 B where it held 472,
// and `TruePeakLimiter::prepare(48000, 65536, 1)` then `(48000, 256, 2)` holds 1 133 236 B where it held
// 88 756. Law 11d's budgets are stated for exactly this object ("one already prepared keeps storage that
// still fits"), and it is the same property that makes a re-preparation free — but a consumer that must
// give a large geometry back destroys the stage rather than re-preparing it wider.

inline void prepareDelayBank (std::vector<DelayLine>& bank, std::size_t lines, int maxDelaySamples)
{
    bank.reserve (lines);                                            // no-op once the bank is wide enough
    while (bank.size() > lines) bank.pop_back();
    while (bank.size() < lines) bank.emplace_back (maxDelaySamples);  // one allocation, at the final size
    for (DelayLine& d : bank)
    {
        d.prepare (maxDelaySamples);                                 // and the lines that were already there
        // THE TAP GOES BACK TO ZERO, because the `assign(n, DelayLine {})` this replaces built every line
        // FRESH and a fresh line has no delay. Keeping a re-used line's tap would make "same end state"
        // false in the one field a preparation does not otherwise touch — and every caller here sets the
        // tap immediately afterwards, so this is the end state they already get, now stated rather than
        // relied upon.
        d.setDelay (0);
    }
}

// What prepareDelayBank() asks the heap for on a FRESH bank: the vector itself plus each line's ring.
inline std::uint64_t delayBankBytes (std::size_t lines, int maxDelaySamples) noexcept
{
    return (std::uint64_t) lines * ((std::uint64_t) sizeof (DelayLine)
                                    + DelayLine::storageFor (maxDelaySamples).bytes());
}

} // namespace felitronics::core

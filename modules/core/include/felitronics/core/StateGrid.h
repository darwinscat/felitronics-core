// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

namespace felitronics::core
{

//==============================================================================
// THE AUDIO-TIME MAINTENANCE GRID — where periodic state maintenance happens, instead of "once per
// process() call".
//
// WHY THIS EXISTS. Law 8's denormal flush was written as "once per block, after the loop, exactly like
// every other kernel in the core". That cadence is chosen by the CALLER, not by the audio, and it fails
// in both directions at once:
//   * the output becomes a function of how the host sliced the stream — a numerical event lands wherever
//     the caller cut, so two renders of the same programme differ (measured on `eq::EqEngine`: 38522 of
//     40000 tail samples differ between a whole-file call and one-sample calls);
//   * on a call that spans a whole file the flush never fires INSIDE it, so the subnormal stall law 8
//     exists to prevent happens in full (37678 of 38000 tail samples subnormal, 10-100x on any CPU
//     without hardware FTZ), and a single non-finite input sample is never healed — measured, one +Inf
//     into an EQ poisons 479900 of the next 480000 samples on a whole-file call against 1 sample when
//     the same stream is fed one sample at a time.
// Two modules in this repository had already refused the per-call pattern in writing, each for its own
// kernel: `analysis::LoudnessMeter` flushes on its 10 ms sub-hop (docs/LAW8-KWEIGHTING.md: "the end of
// process() is wherever the CALLER chose to cut the stream") and `saturation::Saturator` flushes its DC
// blocker per sample. This class generalises what they concluded, so the next kernel is not written to
// the old pattern.
//
// WHAT IT GUARANTEES. `phase_` counts AUDIO samples across calls, so the boundaries fall on the same
// absolute sample indices however the caller slices the stream. A kernel that does all of its periodic
// maintenance at those boundaries is bit-identical under arbitrary re-slicing FROM THE SAME reset().
//
// WHAT IT DOES NOT. It says nothing about events that ARRIVE per call — a parameter write, a bypass
// toggle, a channel-count change. Those are the caller's own timeline and land where the caller puts
// them; the grid makes the maintenance deterministic, not the automation.
//
// THE PERIOD. `kPeriod = 64` samples is chosen, not derived from law 8 — law 8 permits anything in
// [1, ~128] once a boundary is guaranteed to arrive at all, because a state below the 1e-15 flush
// threshold is 23 decades above the subnormal floor and cannot outlive one period. What 64 buys:
//   * subnormal exposure per silence event stays at or below what today's smallest realistic host block
//     already gives (measured worst run inside one period, real pole from 1e-15: 8 samples at K=32,
//     16 at K=64, 180 at K=256, 8116 at K=8192);
//   * it is the smallest block a live rig runs, and a power of two, so every host block that is a
//     multiple of it (64, 128, 256, 512, ...) is a whole number of periods — the boundaries land on the
//     call boundaries and no segment is ever split;
//   * where a kernel also puts its CONTROL rate on this grid it sets the redesign budget: 750 ticks/s
//     at 48 kHz.
// A stall CANNOT be dismissed by pole radius, which is the argument this file used to invite: a DF2T
// biquad with a1 = -65/128, a2 = 9/128 (poles at radius 0.265) walks from 1e-15 to the exact nonzero
// fixed point (z1, z2) = (2^-149, -0) in 53 silent updates and stays there for ever. It is harmless
// HERE only because the next boundary always comes.
//
// OVERSAMPLED KERNELS COUNT THEIR OWN CLOCK. 64 base-rate frames is 512 recursive updates at 8x, so a
// kernel whose state advances at an internal rate drives its grid from THAT rate, not from the frames
// its caller passed.
class StateGrid
{
public:
    static constexpr int kPeriod = 64;

    // A stream restart re-anchors the grid: the first sample after it is grid sample 0.
    void reset() noexcept { phase_ = 0; }

    // How many samples may be processed before the next boundary. Precondition: remaining > 0; the
    // result is then in [1, kPeriod] and never 0, so a segment loop always makes progress.
    int segment (int remaining) const noexcept
    {
        const int room = kPeriod - phase_;
        return remaining < room ? remaining : room;
    }

    // Advance by a segment. Precondition: 0 < n <= segment(...). Returns true exactly when the segment
    // ENDS on a boundary — i.e. when maintenance is due after processing it.
    bool advance (int n) noexcept
    {
        phase_ += n;
        if (phase_ < kPeriod) return false;
        phase_ = 0;
        return true;
    }

    // Advance by `n` samples of audio that produced no maintenance — a stage that was bypassed or ran
    // nothing still consumed audio TIME, and freezing the phase there would make every later boundary
    // depend on how long the bypass lasted. Any size; a non-positive `n` moves nothing, which is not
    // decoration — a negative one would leave the phase negative and `segment()` would then hand out
    // MORE than a period, past the caller's own buffer.
    void skip (int n) noexcept { if (n > 0) phase_ = (phase_ + n % kPeriod) % kPeriod; }

    int phase() const noexcept { return phase_; }

private:
    int phase_ = 0;
};

} // namespace felitronics::core

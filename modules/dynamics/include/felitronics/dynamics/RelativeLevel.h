// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/core/FlushToZero.h>

#include <algorithm>
#include <cmath>

namespace felitronics::dynamics
{

//==============================================================================
// felitronics::dynamics::RelativeLevel — a slow estimate of a signal's PROGRAM level, so a dynamics
// setting can mean the same thing at any gain staging. Feed it a detector envelope; ask it how far
// the current level sits ABOVE that program level. That difference — not an absolute dBFS number —
// is what a "one knob" dynamic EQ drives its gain computer with.
//
// WIRING (deliberate): do NOT rewrite a GainComputer's threshold from this. SET that threshold to 0
// once, and feed it `relativeDb (levelDb)`. The transfer is identical, but nothing then updates a
// threshold at control rate, so a whole class of unit/coefficient bug cannot exist. The steady-state
// gain reduction also becomes a one-line query on the computer itself.
//
//   🔴 GainComputer's DEFAULT threshold is -18 dB, not 0. A consumer that forgets the one-time
//      setThresholdDb(0) gets a permanent gain reduction at relative level 0 (-9 dB for the default
//      2:1) with nothing in the signal to explain it.
//   🔴 When activity() is false, BYPASS the gain computer and target a zero delta. Feeding it a
//      relative level of 0 is NOT neutral: a soft knee is centred on the threshold, so 0 still asks
//      for knee/8 * slope of reduction.
//
// WHAT IT CAN AND CANNOT DO. It observes departures from a signal's OWN recent norm — sibilance, a
// resonance ringing on one note, a boom on one chord. It cannot mark a permanently excessive level
// as wrong: given enough time that level IS the norm. This is a property of same-signal estimation,
// not a defect to engineer around; a fixed threshold is the tool for permanent problems.
//
// RECTIFICATION BIAS — why the averaging is SYMMETRIC. The obvious design is a fast-rise/slow-fall
// follower ("programme meters do that"). Measured, it is a trap: with different up/down constants,
// an observation that swings symmetrically about a steady level takes BIGGER steps up than down, so
// the estimate settles ABOVE the truth. At a +-12 dB swing — ordinary music — that bias measured
// +3.9 dB, i.e. the threshold silently drifts up and the band under-processes exactly the dynamic
// material it exists for. And it buys almost nothing: through a 400 ms inter-phrase gap the
// asymmetric estimator sagged 4.37 dB where the symmetric one sagged 4.80 dB. So: ONE symmetric
// constant, and gap protection done honestly —
//
// GAP SLOWDOWN. An observation more than `gapBelowDb` under the estimate is not a quieter programme,
// it is a pause or a decay tail. There the fall constant is stretched by `gapFactor`, so short gaps
// barely move the estimate while a genuinely quieter section still arrives (slowly) instead of
// latching forever. Inside the normal swing of a programme the estimator stays exactly symmetric, so
// no bias is reintroduced.
//
// The rest exists because it breaks otherwise:
//   • SLEW CAP, BOTH WAYS — a one-pole crosses many dB almost instantly on a 100x step. Which
//     direction is dangerous depends on the consumer's mode (a sinking estimate slams a
//     down-compressor on the next chorus hit; a hanging one rails an up-compressor), so cap both.
//   • FLOOR: HOLD, NOT DECAY — a silent passage must not drag the estimate into the noise and
//     range-cap the first note back. `floorDb` is kept >= `clampLoDb` (sanitised): a programme
//     BELOW the clamp could never be represented anyway, and allowing that grey zone would break
//     translation invariance and leave a hole in the activity gate.
//   • ESTIMATE CLAMP — so it cannot learn a noise floor as "programme".
//   • SEED — the first valid observation sets the estimate outright; crawling from zero would
//     out-mass the signal for several time constants. reset() also ARMS THE FAST WINDOW, because a
//     detector envelope is still climbing its own attack ramp at the first control tick, so the very
//     first observation reads low; the window lets that self-correct in a fraction of a second
//     instead of a slew-capped second and a half.
//   • FAST WINDOW ON RETUNE — a large fc/Q change invalidates what was learned, but reseeding on
//     every parameter write would make the estimator instantaneous while a user DRAGS a node. So a
//     retune opens a bounded window of faster adaptation instead. Both the time constant and the
//     slew cap scale by `fastFactor`, which keeps the crossover between them invariant (the cap
//     binds when the error exceeds slew*tau) — scaling only one would slew-starve the window.
//   • activity() — false while held, WITH HYSTERESIS: a tail hovering at the floor would otherwise
//     flip the flag every control tick and the consumer's gate would chatter. The consumer must gate
//     its gain computer's INPUT on this, not merely stop updating: an up-compressing band otherwise
//     sees a huge "under" during silence and rails its boost onto the noise floor.
//
// 🔴 RT-safe: prepare()/reset() are for stream (re)starts; accumulate()/update() are allocation-,
// lock-, IO- and exception-free. accumulate() is per sample and cheap (one multiply-add); update()
// carries the transcendentals and is meant for control rate.
class RelativeLevel
{
public:
    struct Params
    {
        double timeMs        = 2000.0;   // SYMMETRIC averaging constant (see RECTIFICATION BIAS)
        double gapBelowDb    = 15.0;     // an observation this far under the estimate reads as a gap
        double gapFactor     = 8.0;      // ... and the estimate falls that much slower through it
        double slewDbPerSec  = 12.0;     // hard bidirectional cap on estimate motion
        double clampLoDb     = -80.0;    // the estimate itself stays inside these
        double clampHiDb     =   0.0;
        double floorDb       = -80.0;    // below this: hold + report inactive. Kept >= clampLoDb.
        double activityHystDb = 6.0;     // re-activate only this far above the floor (no chatter)
        double fastWindowMs  = 500.0;    // bounded fast-adaptation window after a retune / reset
        double fastFactor    = 8.0;      // how much faster inside that window (time AND slew)
    };

    void prepare (double sampleRate) noexcept
    {
        fs_ = (std::isfinite (sampleRate) && sampleRate > 0.0) ? sampleRate : 48000.0;
        reset();
    }

    void reset() noexcept
    {
        acc_ = 0.0; accN_ = 0; estDb_ = p_.clampLoDb; seeded_ = false; active_ = false;
        armFastWindow();                 // the first observation is a ramp artefact — let it correct
    }

    void setParams (const Params& p) noexcept { p_ = sanitize (p); }
    const Params& params() const noexcept { return p_; }

    // One detector-envelope sample (LINEAR AMPLITUDE, e.g. EnvelopeFollower::process). We accumulate
    // mean SQUARE: decimating the raw probe instead would beat against the carrier. Non-finite
    // samples are skipped for the MEAN but still counted as elapsed time by the caller's numSamples,
    // which is correct — time passed even if that sample was garbage.
    inline void accumulate (float envAmplitude) noexcept
    {
        const double a = (double) envAmplitude;                  // widen BEFORE squaring
        if (std::isfinite (a)) { acc_ += a * a; ++accN_; }
    }

    // Advance the estimate over the `numSamples` of wall clock just accumulated. Control rate.
    inline void update (int numSamples) noexcept
    {
        if (numSamples <= 0 || accN_ <= 0) { acc_ = 0.0; accN_ = 0; return; }

        const double meanSq = acc_ / (double) accN_;
        acc_ = 0.0; accN_ = 0;

        // 10*log10 of a mean SQUARE == 20*log10 of the rms amplitude — the same units the consumer
        // feeds GainComputer, so relativeDb() subtracts like for like. Floored, never log(0).
        const double obsDb = 10.0 * std::log10 (meanSq > 1.0e-20 ? meanSq : 1.0e-20);

        // Hysteretic activity: leave at floorDb, return only floorDb + activityHystDb.
        const double onDb = p_.floorDb + p_.activityHystDb;
        if (active_ ? (obsDb < p_.floorDb) : (obsDb < onDb)) { active_ = false; return; }   // HOLD
        active_ = true;

        const double dt = (double) numSamples / fs_;
        const bool   fast = fastLeft_ > 0;

        if (! seeded_) { estDb_ = clampEstimate (obsDb); seeded_ = true; }                  // SEED
        else
        {
            // Symmetric by default; only a genuine GAP stretches the constant (see header).
            const bool   gap    = obsDb < estDb_ - p_.gapBelowDb;
            const double tauSec = std::max (1.0e-4, p_.timeMs * 1.0e-3
                                                    * (gap ? p_.gapFactor : 1.0)
                                                    / (fast ? p_.fastFactor : 1.0));

            // Interval-correct coefficient. A PER-SAMPLE coefficient used on a call that covers K
            // samples would silently make the follower K times slower and alias the probe.
            const double coeff  = std::exp (-dt / tauSec);
            double       target = obsDb + coeff * (estDb_ - obsDb);

            const double maxStep = std::max (0.0, p_.slewDbPerSec * (fast ? p_.fastFactor : 1.0)) * dt;
            target = std::clamp (target, estDb_ - maxStep, estDb_ + maxStep);   // BOTH directions

            estDb_ = clampEstimate (target);
        }

        // Spent only while actually tracking, so a retune during a pause still gets its window when
        // the signal returns. The window bounds ADAPTATION time, not wall clock.
        if (fast) fastLeft_ = std::max (0, fastLeft_ - numSamples);
    }

    // How far `levelDb` sits above the programme level — what the gain computer consumes. Meaningless
    // before the first observation, so it reads 0 (transparent) rather than a wild excess.
    inline double relativeDb (double levelDb) const noexcept
    {
        return seeded_ ? levelDb - estDb_ : 0.0;
    }

    double programDb() const noexcept { return estDb_; }
    bool   activity()  const noexcept { return active_; }
    bool   seeded()    const noexcept { return seeded_; }

    // A MATERIAL fc/Q change — not a gain change, and not every parameter write. Opens the bounded
    // fast-adaptation window; deliberately does NOT reseed, so a drag cannot flatten the estimator.
    void retuned() noexcept { armFastWindow(); }

    // Law 8. Belt and braces: accumulate() already rejects non-finite input, and a sum of squares of
    // normal floats cannot itself reach the double subnormal range — this guards the path where a
    // consumer reaches in, not the arithmetic here.
    void flushDenormals() noexcept { core::flushPoison (acc_); }

private:
    static double finite (double v, double fallback) noexcept { return std::isfinite (v) ? v : fallback; }

    void armFastWindow() noexcept { fastLeft_ = (int) std::lround (p_.fastWindowMs * 1.0e-3 * fs_); }

    static Params sanitize (const Params& in) noexcept
    {
        Params p;
        p.timeMs         = std::max (1.0, finite (in.timeMs, 2000.0));
        p.gapBelowDb     = std::max (0.0, finite (in.gapBelowDb, 15.0));
        p.gapFactor      = std::max (1.0, finite (in.gapFactor, 8.0));
        p.slewDbPerSec   = std::max (0.0, finite (in.slewDbPerSec, 12.0));
        p.clampLoDb      =                finite (in.clampLoDb, -80.0);
        p.clampHiDb      = std::max (p.clampLoDb, finite (in.clampHiDb, 0.0));
        // A programme under clampLoDb cannot be represented, so tracking there would pin the estimate
        // and silently break translation invariance — and leave the activity gate blind exactly where
        // an up-compressor would rail onto the noise floor. Hold instead.
        p.floorDb        = std::max (p.clampLoDb, finite (in.floorDb, -80.0));
        p.activityHystDb = std::max (0.0, finite (in.activityHystDb, 6.0));
        p.fastWindowMs   = std::max (0.0, finite (in.fastWindowMs, 500.0));
        p.fastFactor     = std::max (1.0, finite (in.fastFactor, 8.0));
        return p;
    }

    double clampEstimate (double db) const noexcept { return std::clamp (db, p_.clampLoDb, p_.clampHiDb); }

    Params p_;
    double fs_    = 48000.0;
    double acc_   = 0.0;      // running sum of squares since the last update()
    int    accN_  = 0;
    double estDb_ = -80.0;
    bool   seeded_ = false, active_ = false;
    int    fastLeft_ = 0;     // samples of ADAPTATION remaining in the fast window
};

} // namespace felitronics::dynamics

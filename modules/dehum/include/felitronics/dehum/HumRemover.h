// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/core/Config.h>
#include <felitronics/core/DetMath.h>
#include <felitronics/core/FlushToZero.h>
#include <felitronics/core/Math.h>
#include <felitronics/core/Smoother.h>
#include <felitronics/core/StateGrid.h>

#include <cmath>
#include <cstdint>
#include <limits>

namespace felitronics::dehum
{

//==============================================================================
// felitronics::dehum::HumRemover — removes mains hum (50/60 Hz and its harmonics) from an INPUT signal.
//
// WHAT IT IS. A phase-coherent synchronous canceller. One mains phase drives a reference phasor per
// harmonic, r_k = e^{j·k·θ}; every channel keeps one complex weight per harmonic, w_k, and subtracts its
// model of the hum, Re(Σ w_k·r_k), from the input. The weights follow the LMS rule w_k += μ·e·conj(r_k)
// on the residual e. Zero latency, nothing before the input, no heap.
//
// IT IS A NOTCH COMB, EXACTLY (Widrow 1975). With |r| = 1 and a fixed μ the loop is linear and time
// invariant: per harmonic H(z) = (z² − 2z·cos ω_k + 1) / (z² − (2 − μ)·z·cos ω_k + (1 − μ)), poles at
// radius √(1 − μ) ≈ 1 − μ/2, so the −3 dB width is BW = μ·fs/(2π) Hz and the time constant τ = 1/(π·BW).
// μ is set from a width in hertz as μ = 2π·BW/fs, so the width is a CONSTANT ABSOLUTE number of hertz at
// every harmonic and every rate — not a Q. That is the point: a note's partial 0.6 Hz from 660 Hz loses the
// same fraction of a dB as one 0.6 Hz from 60, and widening the upper notches (a constant Q) would eat every
// E's 660 Hz partial. The playing width (`bandwidthHz`, 0.5 Hz) costs ≈ −0.97 dB on a partial 0.5 Hz from a
// harmonic and −0.03 dB at 3 Hz. The nulls sit EXACTLY on k·f_osc whatever the other harmonics do.
//
// TWO WIDTHS. While the input is QUIET the canceller learns at `learnBandwidthHz` (3 Hz, τ 0.11 s); the
// moment it is not, at `bandwidthHz` (0.5 Hz, τ 0.64 s), so a held note is touched only as a 0.5 Hz notch
// touches it. Switching μ moves no weight — the correction, and so the output, is continuous through the
// switch; only the rate of learning changes. QUIET is three conditions together:
//   (a) the last 0.5 s window was hum: its broadband power is within 12 dB of the strongest hum family's
//       Goertzel power (bins {50, 100} and {60, 120} Hz, summed over channels as POWERS);
//   (b) the 20 ms-smoothed broadband power is within 12 dB of the FLOOR — the minimum of that smoothed power
//       over the last 4 s (sixteen 0.25 s hops);
//   (c) no per-sample trigger in the last 50 ms: a sample 12 dB over the floor's sine peak (x² > 32·floor),
//       or — while the canceller runs — a RESIDUAL sample 12 dB over the residual's own 20 ms-smoothed level
//       (e² > 32·max(smoothed e², 1e-4·floor)). The residual is what a pluck shows in first (the hum is
//       already gone from it), so quiet ends on the first few samples of an attack and an attack is never
//       learned at the fast width; a stationary residual (an unconverged model, a ringing ghost) never
//       trips it. The 50 ms latch covers the smoothed power's lag. (a) keeps a held note from becoming
//       "the floor" after 4 s; (b) and (c) react within a period. ACTIVE (a pluck, reported by active()) is the brief's: the 64-sample
//       power above floor + 12 dB for at least 20 ms.
//
// KNOWING THE MAINS. Off / Hz50 / Hz60 / Auto.
//   * AUTO decides the family from the bins on rectangular 0.5 s windows hopped every 0.25 s on the grid: a
//     family must win by more than 6 dB in 3 consecutive QUIET windows (≈ 1 s of hum). Ambiguous → abstain
//     (dry). Decided only in quiet, so it is sticky while playing; re-decided by the same rule.
//   * COARSE frequency: the winning family's stronger bin, phase advance between windows one hop apart
//     (±2 Hz unambiguous), taken only between two quiet windows and only for a family that wins by 6 dB. Auto starts the canceller at it; a fixed
//     mode jumps to it once, the first time it exists ("acquisition").
//   * FINE (an FLL, only in quiet, and only while the running family is the one the last window's bins
//     heard win — a fixed mode on the other family's hum would otherwise chase beats): the weights' own rotation, Δf = Σ k·Im(w_k·conj w_k') / (2πT·Σ k²|w_k|²)
//     over each grid period T, weighted by k·|w_k|² so a weak fundamental does not matter. Gain
//     K = (π/2)·learnBandwidthHz per second: the weights follow the hum's phase through the LMS loop's own
//     lag τ, so the frequency loop is second order, and K·τ = 1/2 is a damping of 0.707. A gain of 0.1 per
//     period (K ≈ 75/s, ζ ≈ 0.18) was measured against it: under −60 dBFS hiss the estimate jitters 2.7×
//     more (0.55 vs 0.20 mHz rms), while on a 0.0025 Hz/s mains ramp it lags less (0.10 vs 0.57 mHz).
//     Slewed to ≤ 0.05 Hz/s, bounded to nominal ± `lockRangeHz` (1 Hz; 0 = run at the nominal frequency).
//     A frequency change moves only the uint32 increment — the phase is continuous, nothing steps. Held
//     whenever the input is not quiet, so a note beside a harmonic cannot pull it.
//   * A change of the family or of the active harmonic set ramps the correction out over 20 ms, then — at
//     the next grid boundary — re-arms with zero weights. Zero weights are a zero correction: no step.
//
// THE GHOST A NOTE LEAVES. While a note sounds Δ Hz from harmonic k, the notch's state carries the note's
// beat: amplitude A·(BW/2)/√(Δ² + (BW/2)²) in w_k. When the note stops, that state rings out as a tone at the
// harmonic, decaying with τ — an LTI fact of any notch of this width, not a defect of the learner. With the
// two widths it decays at τ 0.64 s until the input is quiet again (the first clear 0.5 s window, 0.5–0.75 s
// after the note), then at τ 0.11 s. Measured (48 kHz, a 330 Hz tone 20 Hz from k7 = 350): it starts at
// −38.7 dB re the tone (the LTI value is −38.1), the input is quiet again 0.51 s after the stop, and the ghost
// is below −60 dB re the tone 0.70 s after the stop — against 1.69 s at the playing width alone. The switch
// back to the learn width is a μ change only: the correction's second difference stays at −83 dBFS through
// it, the same as without the switch (HumRemoverTrackingTests).
//
// NUMERICS (law 3 carve-out, stated): weights and references in DOUBLE, samples in float. Float weights
// stall where an increment drops below half an ulp of the weight (the design study read −63 dB re the hum at
// 48 kHz, −51 at 192 kHz). The mains phase is a uint32 advanced by round(f·2³²/fs); every 64 samples
// (core::StateGrid) each reference is RE-ANCHORED from that phase through core::det sin/cos — k·phase wraps
// exactly in uint32 arithmetic — and between anchors it advances by one complex rotation in double. All
// clocks, decisions and tracker steps live on the grid or on the sample, so the output is BIT-IDENTICAL
// under any re-slicing of the stream (law 8a). The coarse estimate uses std::atan2 (libm), so rows agree
// within tolerance, not bit for bit.
// LAW 10: the weight update is spelled so nothing in it can contract under ANY mode: the product μ·e·r is
// ROUNDED TO FLOAT and stored, then added to the double weight (detail::adapt). The rest of the sample path
// is written one operation per statement, so the repository's -ffp-contract=on fuses nothing there either.
//
// SILENCE, POISON, STOPPED CHANNELS:
//   * An exact-zero input sample ADAPTS NOTHING and feeds no detector — digital zero is the absence of a
//     measurement. After kBlackRun (64) consecutive exact zeros the correction is MUTED too, so digital black
//     in is bit-exact black out; weights, floor and decisions hold. The first 63 zeros still carry the
//     correction (the anti-hum of the stream that just stopped), then it cuts.
//   * A non-finite input sample passes through untouched and adapts / feeds nothing. Law 8/8a: on the grid
//     each weight is flushed (|w| < 1e-15 → 0), healed (non-finite → 0) and bounded (beyond kWeightCeiling,
//     +60 dBFS, it is not a hum model → 0); the non-finite half also runs at the end of every call.
//   * Law 11a: a channel at index >= nch STOPPED for that call. Its weights are a MODEL of the hum on that
//     input, not audio memory, so they are KEPT; its zero-run counter advances (it is receiving silence).
//   * Law 11c: a clock-only call (nch == 0) is digital silence for every shared clock — the grid, the mains
//     and bin phases, the hop and window counters, the ramp and a pending re-arm advance exactly as through
//     n zeros at a live width; weights and detector statistics hold, as they do through zeros.
//
// WHAT IT CANNOT DO: remove hum that is made after this stage's input (e.g. inside a physical amplifier
// downstream), harmonics above the active set, or hum whose phase moves faster than τ lets the weights follow.
//
// RT-safe: no allocation, lock, IO or throw anywhere after construction — there is no heap at all. Footprint
// sizeof(HumRemover): weights and their previous-period copy 4·kMaxChannels·MaxHarmonics doubles (32 KiB at
// 64), plus ≈ 3 KiB of detector state.
//==============================================================================

enum class Mains { Off, Hz50, Hz60, Auto };

struct HumRemoverParams
{
    Mains  mains            = Mains::Auto;
    int    harmonics        = 24;            // [1, MaxHarmonics]; harmonics above maxHz or 0.9·Nyquist are skipped
    double bandwidthHz      = 0.5;           // notch −3 dB width while playing; [kMinBandwidthHz, kMaxBandwidthHz]
    double learnBandwidthHz = 3.0;           // notch width while the input is quiet; same range
    double maxHz            = 2000.0;        // highest harmonic frequency processed
    double lockRangeHz      = 1.0;           // tracker bound around 50/60; [0, 2]; 0 = the nominal frequency
};

namespace detail
{
    // The LMS step w += muErr·conj(r), i.e. re += muErr·rRe, im −= muErr·rIm. Each product is ROUNDED TO
    // FLOAT and stored before the double add, so no contraction mode can fuse it (law 10): an FMA would
    // have to skip the float rounding. The increment keeps 24 bits of its own precision; the weight keeps
    // 53, which is what defeats the float-weight stall.
    inline void adapt (double& wRe, double& wIm, double muErr, double rRe, double rIm) noexcept
    {
        const float dRe = (float) (muErr * rRe);
        const float dIm = (float) (muErr * rIm);
        wRe += (double) dRe;
        wIm -= (double) dIm;
    }

    // e^{j·2π·u/2³²}. The uint32 is the exact phase (wrapped by unsigned arithmetic); it is read as a signed
    // angle in (−π, π] so det::sin/cos see the smallest argument.
    inline void phasor (std::uint32_t u, double& re, double& im) noexcept
    {
        constexpr double kRadPerUnit = 2.0 * core::kPi / 4294967296.0;
        const double a = (double) (std::int32_t) u * kRadPerUnit;
        re = core::det::cos (a);
        im = core::det::sin (a);
    }
} // namespace detail

template <int MaxHarmonics = 64>
class HumRemover
{
    static_assert (MaxHarmonics >= 1, "HumRemover needs at least one harmonic");

public:
    static constexpr int    kMaxHarmonics   = MaxHarmonics;
    static constexpr int    kBlackRun       = 64;       // exact zeros in a row that mute the correction
    static constexpr double kRampSeconds    = 0.020;    // correction ramp-out on a family / harmonic-set change
    static constexpr double kMinBandwidthHz = 0.01;
    static constexpr double kMaxBandwidthHz = 20.0;
    static constexpr double kMaxLockRangeHz = 2.0;
    static constexpr double kWeightCeiling  = 1.0e3;    // +60 dBFS: beyond it a weight is not a hum model
    static constexpr double kHopSeconds     = 0.25;     // detector hop; a window is two hops
    static constexpr int    kFloorHops      = 16;       // the floor is the minimum over this many hops (4 s)
    static constexpr int    kWinStreak      = 3;        // consecutive quiet winning windows that decide Auto
    static constexpr double kQuietRatio     = 16.0;     // 12 dB, in power
    static constexpr double kWinRatio       = 3.981071705534973;   // 6 dB, in power
    static constexpr double kHitRatio       = 32.0;     // x² vs floor: 12 dB over the floor's sine peak
    static constexpr double kSmoothSeconds  = 0.020;    // the broadband power's smoothing
    static constexpr double kActiveSeconds  = 0.020;    // above floor + 12 dB this long = active
    static constexpr double kHoldSeconds    = 0.050;    // a trigger keeps the playing width at least this long
    static constexpr double kSlewHzPerSec   = 0.05;     // the tracker's slew bound

    struct Weight { double re, im; };

    // Law 11b: disarm, validate, write. Refuses fs < core::kMinSampleRate (or non-finite), maxBlock < 1 and
    // channels outside [1, core::kMaxChannels]. A prepare is a restart (law 11a): weights, clocks, detector.
    [[nodiscard]] bool prepare (double sampleRate, int maxBlock, int maxChannels) noexcept
    {
        prepared_ = false;
        if (! (sampleRate >= core::kMinSampleRate) || ! std::isfinite (sampleRate)) return false;
        if (maxBlock < 1) return false;
        if (maxChannels < 1 || maxChannels > core::kMaxChannels) return false;
        fs_       = sampleRate;
        maxBlock_ = maxBlock;
        maxCh_    = maxChannels;
        ramp_.reset (fs_, kRampSeconds);
        periodSec_     = (double) kPeriod / fs_;
        smoothA_       = 1.0 - core::det::exp2 (-periodSec_ / (kSmoothSeconds * 0.6931471805599453));
        const long long hp = std::llround (kHopSeconds * fs_ / (double) kPeriod);
        hopPeriods_    = hp < 1 ? 1 : (int) hp;
        hopSec_        = (double) (hopPeriods_ * kPeriod) / fs_;
        activeSamples_ = (int) std::llround (kActiveSeconds * fs_);
        holdSamples_   = (int) std::llround (kHoldSeconds * fs_);
        for (int b = 0; b < 4; ++b) binInc_[b] = incFor (kBinHz[b]);
        prepared_ = true;
        reset();
        return true;
    }

    bool isPrepared()  const noexcept { return prepared_; }
    int  maxChannels() const noexcept { return prepared_ ? maxCh_ : 0; }
    int  maxBlock()    const noexcept { return prepared_ ? maxBlock_ : 0; }
    static constexpr int latencySamples() noexcept { return 0; }

    // A stream restart: zero weights, re-anchored clocks, a fresh detector (Auto undecided), the current
    // params applied at once (no ramp), zero-run counters cleared.
    void reset() noexcept
    {
        grid_.reset();
        anchorPhase_ = 0;
        for (int b = 0; b < 4; ++b) binAnchor_[b] = 0;
        for (int c = 0; c < core::kMaxChannels; ++c) zeros_[c] = 0;
        resetDetector();
        applyWanted();
    }

    // The widths apply at once (μ moves no weight); a family or harmonic-set change ramps out and re-arms;
    // a narrower lock range pulls the estimate inside it at once (phase-continuous).
    void setParams (const HumRemoverParams& p) noexcept
    {
        const bool wasOff = params_.mains == Mains::Off;
        params_.mains            = sanitizeMains (p.mains);
        params_.harmonics        = p.harmonics < 1 ? 1 : (p.harmonics > MaxHarmonics ? MaxHarmonics : p.harmonics);
        params_.bandwidthHz      = sanitizeBw (p.bandwidthHz, params_.bandwidthHz);
        params_.learnBandwidthHz = sanitizeBw (p.learnBandwidthHz, params_.learnBandwidthHz);
        params_.maxHz            = (std::isfinite (p.maxHz) && p.maxHz >= 0.0) ? p.maxHz : params_.maxHz;
        params_.lockRangeHz      = std::isfinite (p.lockRangeHz)
                                       ? (p.lockRangeHz < 0.0 ? 0.0 : (p.lockRangeHz > kMaxLockRangeHz ? kMaxLockRangeHz : p.lockRangeHz))
                                       : params_.lockRangeHz;
        if (! prepared_) return;
        updateRates();
        if (wasOff && params_.mains != Mains::Off) resetDetector();
        if (running_) setFrequency (clampToLock (fHat_));
        retarget();
    }

    const HumRemoverParams& params() const noexcept { return params_; }

    // Law 11: malformed → unprepared → nch > maxChannels → n == 0 → nch == 0 (clock only) → run.
    [[nodiscard]] bool process (core::Sample* const* io, int nch, int n) noexcept
    {
        if (n < 0 || nch < 0) return false;
        if (! prepared_)      return false;
        if (nch > maxCh_)     return false;
        if (n == 0)           return true;

        // Law 11a: every plane at index >= nch stopped for these n samples and is receiving digital silence.
        for (int c = nch; c < maxCh_; ++c)
            zeros_[c] = (n >= kBlackRun - zeros_[c]) ? kBlackRun : zeros_[c] + n;

        for (int done = 0; done < n;)
        {
            const int seg = grid_.segment (n - done);
            if (nch > 0 && running_)      runSegment (io, nch, done, seg);
            else if (nch > 0 && detOn())  detectSegment (io, nch, done, seg);
            else                          clockSegment (seg);
            if (grid_.advance (seg)) atGridBoundary();
            done += seg;
        }
        healPoison();   // law 8a(a): the poison half keeps its own clock — the end of every call
        return true;
    }

    // ---- inspection (process thread) -------------------------------------------------------------
    bool   running()         const noexcept { return running_; }
    int    activeHarmonics() const noexcept { return nActive_; }
    // The tracker's estimate of the mains frequency, or 0 while dry.
    double mainsHz()         const noexcept { return running_ ? fHat_ : 0.0; }
    // The frequency the oscillator actually runs at (the estimate rounded to fs/2³²), or 0 while dry.
    double oscillatorHz()    const noexcept { return running_ ? (double) inc_ * fs_ / 4294967296.0 : 0.0; }
    // What Auto has decided (Hz50 / Hz60), or Off while it abstains — whatever the mode.
    Mains  detected()        const noexcept { return autoFamily_ == 50 ? Mains::Hz50 : (autoFamily_ == 60 ? Mains::Hz60 : Mains::Off); }
    bool   quiet()           const noexcept { return quietPeriod_; }
    bool   active()          const noexcept { return activeRun_ >= activeSamples_; }
    Weight weight (int channel, int harmonic) const noexcept
    {
        if (channel < 0 || channel >= core::kMaxChannels || harmonic < 1 || harmonic > MaxHarmonics) return { 0.0, 0.0 };
        return { wRe_[channel][harmonic - 1], wIm_[channel][harmonic - 1] };
    }

private:
    static constexpr int    kPeriod  = core::StateGrid::kPeriod;
    static constexpr double kBinHz[4] = { 50.0, 100.0, 60.0, 120.0 };
    static constexpr double kNone    = std::numeric_limits<double>::max();

    static Mains sanitizeMains (Mains m) noexcept
    {
        switch (m)
        {
            case Mains::Off:  return Mains::Off;
            case Mains::Hz50: return Mains::Hz50;
            case Mains::Hz60: return Mains::Hz60;
            case Mains::Auto: return Mains::Auto;
        }
        return Mains::Off;
    }

    static double sanitizeBw (double v, double fallback) noexcept
    {
        if (! std::isfinite (v)) return fallback;
        return v < kMinBandwidthHz ? kMinBandwidthHz : (v > kMaxBandwidthHz ? kMaxBandwidthHz : v);
    }

    std::uint32_t incFor (double hz) const noexcept { return (std::uint32_t) std::llround (hz * 4294967296.0 / fs_); }
    // The detector listens whenever a mode wants it — and while the canceller still runs (a ramp-out to Off).
    bool detOn() const noexcept { return params_.mains != Mains::Off || running_; }

    int wantedFamily() const noexcept
    {
        switch (params_.mains)
        {
            case Mains::Off:  return 0;
            case Mains::Hz50: return 50;
            case Mains::Hz60: return 60;
            case Mains::Auto: return autoFamily_;
        }
        return 0;
    }

    int activeFor (int family) const noexcept
    {
        if (family <= 0) return 0;
        const double f0  = (double) family;
        const double top = params_.maxHz < 0.45 * fs_ ? params_.maxHz : 0.45 * fs_;
        const double fit = std::floor (top / f0);
        const int    k   = fit < (double) params_.harmonics ? (int) fit : params_.harmonics;
        return k < 0 ? 0 : k;
    }

    double clampToLock (double f) const noexcept
    {
        const double lo = (double) family_ - params_.lockRangeHz, hi = (double) family_ + params_.lockRangeHz;
        return f < lo ? lo : (f > hi ? hi : f);
    }

    void updateRates() noexcept
    {
        muSlow_    = 2.0 * core::kPi * params_.bandwidthHz / fs_;
        muFast_    = 2.0 * core::kPi * params_.learnBandwidthHz / fs_;
        fllGainT_  = 0.5 * core::kPi * params_.learnBandwidthHz * periodSec_;   // K·T, K·τ_learn = 1/2
        slewT_     = kSlewHzPerSec * periodSec_;
    }

    // Wanted vs applied: anything structural pending ramps the correction out (or keeps it out).
    void retarget() noexcept
    {
        const int fam = wantedFamily();
        pending_ = ! (fam == family_ && activeFor (fam) == nActive_);
        ramp_.setTargetValue (pending_ || ! running_ ? 0.0f : 1.0f);
    }

    // The oscillator at f: only the increment and the per-sample rotations change, never the phase.
    void setFrequency (double f) noexcept
    {
        fHat_ = f;
        const std::uint32_t ni = f > 0.0 ? incFor (f) : 0u;
        if (ni == inc_ && rotsValid_) return;
        inc_ = ni;
        detail::phasor (inc_, rotRe_[0], rotIm_[0]);
        for (int k = 1; k < nActive_; ++k)
        {
            const double a = rotRe_[k - 1] * rotRe_[0];
            const double b = rotIm_[k - 1] * rotIm_[0];
            const double c = rotRe_[k - 1] * rotIm_[0];
            const double d = rotIm_[k - 1] * rotRe_[0];
            rotRe_[k] = a - b;
            rotIm_[k] = c + d;
        }
        rotsValid_ = true;
    }

    // Re-arm on the current params and decision: zero weights, the family's frequency (the coarse estimate
    // when there is a fresh one), the gain at its resting value.
    void applyWanted() noexcept
    {
        updateRates();
        family_   = wantedFamily();
        nActive_  = activeFor (family_);
        running_  = nActive_ > 0;
        acquired_ = false;
        double f  = (double) family_;
        const int s = family_ == 60 ? 1 : 0;
        if (running_ && params_.lockRangeHz > 0.0 && coarseValid_[s]) { f = clampToLock (coarse_[s]); acquired_ = true; }
        rotsValid_ = false;
        setFrequency (running_ ? f : 0.0);
        for (int c = 0; c < core::kMaxChannels; ++c)
            for (int k = 0; k < MaxHarmonics; ++k) { wRe_[c][k] = 0.0; wIm_[c][k] = 0.0; pRe_[c][k] = 0.0; pIm_[c][k] = 0.0; }
        refsValid_ = false;
        pending_   = false;
        psE_ = 0.0; psEValid_ = false; periodSumE_ = 0.0;
        ramp_.setCurrentAndTargetValue (running_ ? 1.0f : 0.0f);
    }

    void resetDetector() noexcept
    {
        for (int c = 0; c < core::kMaxChannels; ++c)
            for (int b = 0; b < 4; ++b)
            { hopRe_[c][b] = 0.0; hopIm_[c][b] = 0.0; prevHopRe_[c][b] = 0.0; prevHopIm_[c][b] = 0.0; winRe_[c][b] = 0.0; winIm_[c][b] = 0.0; }
        hopPow_ = 0.0; prevHopPow_ = 0.0; periodSum_ = 0.0;
        hopCount_ = 0; hopsSeen_ = 0;
        prevWinValid_ = false; prevWinQuiet_ = false;
        ps_ = 0.0; psValid_ = false; hopMin_ = kNone;
        for (int i = 0; i < kFloorHops; ++i) floorRing_[i] = 0.0;
        floorCount_ = 0; floorPos_ = 0; floor_ = 0.0; floorValid_ = false;
        hit_ = false; hitHold_ = 0; activeRun_ = 0; quietWin_ = false; quietPeriod_ = false; thr2_ = kNone; thrE2_ = kNone;
        streak_ = 0; streakFamily_ = 0; autoFamily_ = 0; winner_ = 0;
        coarse_[0] = coarse_[1] = 0.0; coarseValid_[0] = coarseValid_[1] = false;
    }

    // The references at grid phase p: the anchor of the current period, rotated p times. Lazy, so a stretch
    // nobody listens to (black, a pause) costs no rotation — and the arithmetic is the same sequence either way.
    void syncRefs (int p) noexcept
    {
        if (! refsValid_)
        {
            for (int k = 0; k < nActive_; ++k)
                detail::phasor ((std::uint32_t) (k + 1) * anchorPhase_, rRe_[k], rIm_[k]);
            refStep_   = 0;
            refsValid_ = true;
        }
        for (; refStep_ < p; ++refStep_)
            for (int k = 0; k < nActive_; ++k)
                rotate (rRe_[k], rIm_[k], rotRe_[k], rotIm_[k]);
    }

    // The detector's four bin phasors, the same way.
    void syncBins (int p) noexcept
    {
        if (! binsValid_)
        {
            for (int b = 0; b < 4; ++b)
            {
                detail::phasor (binAnchor_[b], bRe_[b], bIm_[b]);
                detail::phasor (binInc_[b], bRotRe_[b], bRotIm_[b]);
            }
            binStep_   = 0;
            binsValid_ = true;
        }
        for (; binStep_ < p; ++binStep_)
            for (int b = 0; b < 4; ++b)
                rotate (bRe_[b], bIm_[b], bRotRe_[b], bRotIm_[b]);
    }

    static void rotate (double& re, double& im, double c, double s) noexcept
    {
        const double a = re * c;
        const double b = im * s;
        const double d = re * s;
        const double e = im * c;
        re = a - b;
        im = d + e;
    }

    double correction (int c) const noexcept
    {
        double acc = 0.0;
        for (int k = 0; k < nActive_; ++k)
        {
            const double a = wRe_[c][k] * rRe_[k];
            const double b = wIm_[c][k] * rIm_[k];
            const double t = a - b;
            acc += t;
        }
        return acc;
    }

    // The detector's view of one sample of every live channel: broadband power, the per-sample trigger, the
    // four bins. Zeros and non-finite samples are no measurement.
    void observe (core::Sample* const* io, int nch, int idx, int p) noexcept
    {
        for (int c = 0; c < nch; ++c)
        {
            const core::Sample x = io[c][idx];
            if (core::exactlyEqual (x, 0.0f) || ! std::isfinite (x)) continue;
            const double xd = (double) x;
            const double x2 = xd * xd;
            periodSum_ += x2;
            if (x2 > thr2_) hit_ = true;
            syncBins (p);
            for (int b = 0; b < 4; ++b)
            {
                const double a = xd * bRe_[b];
                const double s = xd * bIm_[b];
                hopRe_[c][b] += a;
                hopIm_[c][b] -= s;
            }
        }
    }

    void runSegment (core::Sample* const* io, int nch, int off, int seg) noexcept
    {
        const int  p0  = grid_.phase();
        const bool det = detOn();
        for (int i = 0; i < seg; ++i)
        {
            const double g = (double) ramp_.getNextValue();
            if (det) observe (io, nch, off + i, p0 + i);
            for (int c = 0; c < nch; ++c)
            {
                core::Sample& s = io[c][off + i];
                const core::Sample x = s;
                const bool zero = core::exactlyEqual (x, 0.0f);
                if (zero)
                {
                    if (zeros_[c] < kBlackRun) ++zeros_[c];
                    if (zeros_[c] >= kBlackRun) continue;             // digital black: muted, weights hold
                }
                else
                    zeros_[c] = 0;
                if (! std::isfinite (x)) continue;                   // passes untouched, adapts nothing

                syncRefs (p0 + i);
                const double corr = correction (c);
                const double gc   = g * corr;
                s = (core::Sample) ((double) x - gc);
                if (zero) continue;                                  // a zero is no measurement

                const double err = (double) x - corr;
                if (det)
                {
                    const double e2 = err * err;
                    periodSumE_ += e2;
                    if (e2 > thrE2_) hit_ = true;
                }
                const double mu    = (quietPeriod_ && ! hit_) ? muFast_ : muSlow_;
                const double muErr = mu * err;
                for (int k = 0; k < nActive_; ++k)
                    detail::adapt (wRe_[c][k], wIm_[c][k], muErr, rRe_[k], rIm_[k]);
            }
        }
    }

    // Dry, but listening (Auto before a decision): the detector sees the input, the output is untouched.
    void detectSegment (core::Sample* const* io, int nch, int off, int seg) noexcept
    {
        const int p0 = grid_.phase();
        for (int i = 0; i < seg; ++i)
        {
            if (ramp_.isSmoothing()) (void) ramp_.getNextValue();
            observe (io, nch, off + i, p0 + i);
        }
    }

    // Nothing is listening: only the ramp moves, sample by sample exactly as in runSegment.
    void clockSegment (int seg) noexcept
    {
        if (ramp_.isSmoothing())
            for (int i = 0; i < seg; ++i) (void) ramp_.getNextValue();
    }

    void atGridBoundary() noexcept
    {
        const std::uint32_t P = (std::uint32_t) kPeriod;
        anchorPhase_ += P * inc_;                                   // exact, mod 2³²
        for (int b = 0; b < 4; ++b) binAnchor_[b] += P * binInc_[b];
        refsValid_ = false;
        binsValid_ = false;
        if (detOn()) endPeriod();
        for (int c = 0; c < maxCh_; ++c)
            for (int k = 0; k < nActive_; ++k)
            {
                double& re = wRe_[c][k];
                double& im = wIm_[c][k];
                if (! (std::fabs (re) <= kWeightCeiling && std::fabs (im) <= kWeightCeiling)) { re = 0.0; im = 0.0; }
                core::flushDenormal (re);
                core::flushDenormal (im);
                pRe_[c][k] = re;
                pIm_[c][k] = im;
            }
        if (pending_ && ! ramp_.isSmoothing() && core::exactlyEqual (ramp_.getCurrentValue(), 0.0f))
            applyWanted();
    }

    void endPeriod() noexcept
    {
        const bool   wasQuiet = quietPeriod_ && ! hit_;       // the whole period ran at the learn width
        const double sum      = periodSum_;
        const double sumE     = periodSumE_;
        periodSum_  = 0.0;
        periodSumE_ = 0.0;
        hopPow_   += sum;
        if (sum > 0.0)                                         // an all-zero period (black, a pause) holds
        {
            const double p = sum / (double) kPeriod;
            if (! psValid_) { ps_ = p; psValid_ = true; }
            else
            {
                const double d = p - ps_;
                const double s = smoothA_ * d;
                ps_ += s;
            }
            const bool above = floorValid_ && p > kQuietRatio * floor_;
            activeRun_ = above ? (activeRun_ < activeSamples_ ? activeRun_ + kPeriod : activeRun_) : 0;
            if (ps_ < hopMin_) hopMin_ = ps_;
            if (running_ && wasQuiet && winner_ == family_ && params_.lockRangeHz > 0.0) trackFrequency();
        }
        if (sumE > 0.0)                                        // the residual's smoothed level (running only)
        {
            const double pe = sumE / (double) kPeriod;
            if (! psEValid_) { psE_ = pe; psEValid_ = true; }
            else
            {
                const double d = pe - psE_;
                const double s = smoothA_ * d;
                psE_ += s;
            }
        }
        if (++hopCount_ >= hopPeriods_) endHop();
        if (hit_) hitHold_ = holdSamples_;
        else if (hitHold_ > 0) hitHold_ -= kPeriod;
        quietPeriod_ = quietWin_ && floorValid_ && ! (ps_ > kQuietRatio * floor_) && hitHold_ <= 0;
        hit_  = false;
        thr2_ = floorValid_ ? kHitRatio * floor_ : kNone;
        if (psEValid_ && floorValid_)
        {
            const double lo = 1e-4 * floor_;
            thrE2_ = kHitRatio * (psE_ > lo ? psE_ : lo);
        }
        else thrE2_ = thr2_;
    }

    // The FLL: the weights' rotation over this period, weighted by k·|w|², one gain-and-slew step.
    void trackFrequency() noexcept
    {
        double num = 0.0, den = 0.0;
        for (int c = 0; c < maxCh_; ++c)
            for (int k = 0; k < nActive_; ++k)
            {
                const double kk = (double) (k + 1);
                const double a  = wIm_[c][k] * pRe_[c][k];
                const double b  = wRe_[c][k] * pIm_[c][k];
                const double x  = a - b;                       // Im(w · conj w')
                const double t  = kk * x;
                num += t;
                const double m1 = wRe_[c][k] * wRe_[c][k];
                const double m2 = wIm_[c][k] * wIm_[c][k];
                const double m  = m1 + m2;
                const double kq = kk * kk;
                const double v  = kq * m;
                den += v;
            }
        if (! (den > 1e-16)) return;
        const double delta = num / (2.0 * core::kPi * periodSec_ * den);
        double step = fllGainT_ * delta;
        if (step >  slewT_) step =  slewT_;
        if (step < -slewT_) step = -slewT_;
        setFrequency (clampToLock (fHat_ + step));
    }

    void endHop() noexcept
    {
        hopCount_ = 0;
        if (hopsSeen_ >= 1)
        {
            const double nW = 2.0 * (double) (hopPeriods_ * kPeriod);
            double pw[4] = { 0.0, 0.0, 0.0, 0.0 };
            for (int c = 0; c < maxCh_; ++c)
                for (int b = 0; b < 4; ++b)
                {
                    const double re = prevHopRe_[c][b] + hopRe_[c][b];
                    const double im = prevHopIm_[c][b] + hopIm_[c][b];
                    const double a = re * re;
                    const double d = im * im;
                    const double m = a + d;
                    pw[b] += m;
                    nowRe_[c][b] = re;
                    nowIm_[c][b] = im;
                }
            const double norm  = 2.0 / (nW * nW);
            const double fam50 = (pw[0] + pw[1]) * norm;
            const double fam60 = (pw[2] + pw[3]) * norm;
            const double pm    = (prevHopPow_ + hopPow_) / nW;
            const double hum   = fam50 > fam60 ? fam50 : fam60;
            const bool   quiet = pm > 0.0 && pm < kQuietRatio * hum;

            const bool wins[2] = { fam50 > kWinRatio * fam60, fam60 > kWinRatio * fam50 };
            for (int s = 0; s < 2; ++s)                        // coarse frequency, per family, only where it is the hum
            {
                coarseValid_[s] = false;
                if (! (quiet && wins[s] && prevWinValid_ && prevWinQuiet_)) continue;
                const int b = pw[2 * s] >= pw[2 * s + 1] ? 2 * s : 2 * s + 1;
                double cr = 0.0, ci = 0.0;
                for (int c = 0; c < maxCh_; ++c)
                {
                    const double a = nowRe_[c][b] * winRe_[c][b];
                    const double d = nowIm_[c][b] * winIm_[c][b];
                    const double e = nowIm_[c][b] * winRe_[c][b];
                    const double f = nowRe_[c][b] * winIm_[c][b];
                    const double r = a + d;
                    const double q = e - f;
                    cr += r;
                    ci += q;
                }
                const double c2a = cr * cr;
                const double c2b = ci * ci;
                if (! (c2a + c2b > 0.0)) continue;
                const double off = std::atan2 (ci, cr) / (2.0 * core::kPi * hopSec_);
                coarse_[s]      = (kBinHz[b] + off) / (double) ((b & 1) + 1);
                coarseValid_[s] = true;
            }

            const int winner = wins[0] ? 50 : (wins[1] ? 60 : 0);
            winner_ = winner;
            if (quiet && winner != 0) { streak_ = winner == streakFamily_ ? streak_ + 1 : 1; streakFamily_ = winner; }
            else                      { streak_ = 0; streakFamily_ = 0; }
            if (streak_ >= kWinStreak && winner != autoFamily_)
            {
                autoFamily_ = winner;
                if (params_.mains == Mains::Auto) retarget();
            }
            if (running_ && ! acquired_ && params_.lockRangeHz > 0.0 && coarseValid_[family_ == 60 ? 1 : 0])
            {
                setFrequency (clampToLock (coarse_[family_ == 60 ? 1 : 0]));
                acquired_ = true;
            }
            quietWin_ = quiet;
            for (int c = 0; c < maxCh_; ++c)
                for (int b = 0; b < 4; ++b) { winRe_[c][b] = nowRe_[c][b]; winIm_[c][b] = nowIm_[c][b]; }
            prevWinValid_ = true;
            prevWinQuiet_ = quiet;
        }
        if (hopMin_ < kNone)                                   // the floor: minimum over the last kFloorHops
        {
            floorRing_[floorPos_] = hopMin_;
            floorPos_ = (floorPos_ + 1) % kFloorHops;
            if (floorCount_ < kFloorHops) ++floorCount_;
            double m = floorRing_[0];
            for (int i = 1; i < floorCount_; ++i) if (floorRing_[i] < m) m = floorRing_[i];
            floor_ = m;
            floorValid_ = true;
        }
        hopMin_ = kNone;
        for (int c = 0; c < maxCh_; ++c)
            for (int b = 0; b < 4; ++b)
            { prevHopRe_[c][b] = hopRe_[c][b]; prevHopIm_[c][b] = hopIm_[c][b]; hopRe_[c][b] = 0.0; hopIm_[c][b] = 0.0; }
        prevHopPow_ = hopPow_;
        hopPow_     = 0.0;
        if (hopsSeen_ < 2) ++hopsSeen_;
    }

    void healPoison() noexcept
    {
        for (int c = 0; c < maxCh_; ++c)
            for (int k = 0; k < nActive_; ++k)
                if (! (std::isfinite (wRe_[c][k]) && std::isfinite (wIm_[c][k]))) { wRe_[c][k] = 0.0; wIm_[c][k] = 0.0; }
    }

    // ---- canceller ----
    // Weights per channel, per harmonic (index k-1), and their value at the last grid boundary (the FLL's).
    double wRe_[core::kMaxChannels][MaxHarmonics] {};
    double wIm_[core::kMaxChannels][MaxHarmonics] {};
    double pRe_[core::kMaxChannels][MaxHarmonics] {};
    double pIm_[core::kMaxChannels][MaxHarmonics] {};
    double rRe_[MaxHarmonics] {}, rIm_[MaxHarmonics] {};
    double rotRe_[MaxHarmonics] {}, rotIm_[MaxHarmonics] {};
    int    zeros_[core::kMaxChannels] {};

    // ---- detector ----
    double hopRe_[core::kMaxChannels][4] {}, hopIm_[core::kMaxChannels][4] {};
    double prevHopRe_[core::kMaxChannels][4] {}, prevHopIm_[core::kMaxChannels][4] {};
    double winRe_[core::kMaxChannels][4] {}, winIm_[core::kMaxChannels][4] {};   // the previous window
    double nowRe_[core::kMaxChannels][4] {}, nowIm_[core::kMaxChannels][4] {};   // scratch: this window
    double bRe_[4] {}, bIm_[4] {}, bRotRe_[4] {}, bRotIm_[4] {};
    std::uint32_t binInc_[4] {}, binAnchor_[4] {};
    double floorRing_[kFloorHops] {};
    double hopPow_ = 0.0, prevHopPow_ = 0.0, periodSum_ = 0.0;
    double ps_ = 0.0, hopMin_ = kNone, floor_ = 0.0, thr2_ = kNone;
    double psE_ = 0.0, periodSumE_ = 0.0, thrE2_ = kNone;
    double coarse_[2] {};
    int    hopCount_ = 0, hopsSeen_ = 0, hopPeriods_ = 1, binStep_ = 0;
    int    floorCount_ = 0, floorPos_ = 0, activeRun_ = 0, activeSamples_ = 960, hitHold_ = 0, holdSamples_ = 2400;
    int    streak_ = 0, streakFamily_ = 0, autoFamily_ = 0, winner_ = 0;   // winner_: the last window's family
    bool   coarseValid_[2] {};
    bool   prevWinValid_ = false, prevWinQuiet_ = false, psValid_ = false, floorValid_ = false, psEValid_ = false;
    bool   hit_ = false, quietWin_ = false, quietPeriod_ = false, binsValid_ = false;

    HumRemoverParams     params_ {};
    core::StateGrid      grid_ {};
    core::LinearSmoother ramp_ { 0.0f };
    std::uint32_t        anchorPhase_ = 0;   // the mains phase at the start of the current grid period
    std::uint32_t        inc_         = 0;   // per-sample phase increment, round(f·2³²/fs)
    double               fs_          = 48000.0;
    double               periodSec_   = 64.0 / 48000.0;
    double               hopSec_      = 0.25;
    double               smoothA_     = 0.0;
    double               muSlow_      = 0.0, muFast_ = 0.0, fllGainT_ = 0.0, slewT_ = 0.0;
    double               fHat_        = 0.0;
    int                  family_      = 0;   // 0 (dry), 50 or 60 — what the canceller runs
    int                  nActive_     = 0;
    int                  refStep_     = 0;
    int                  maxCh_       = 0;
    int                  maxBlock_    = 0;
    bool                 refsValid_   = false;
    bool                 rotsValid_   = false;
    bool                 running_     = false;
    bool                 pending_     = false;
    bool                 acquired_    = false;
    bool                 prepared_    = false;
};

} // namespace felitronics::dehum

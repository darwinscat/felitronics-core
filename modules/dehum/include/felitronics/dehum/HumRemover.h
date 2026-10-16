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
//   (b) the canceller's RESIDUAL (20 ms-smoothed) is within 12 dB of the RESIDUAL FLOOR — the minimum of
//       that power over the last 4 s (sixteen 0.25 s hop minima), which may RISE by at most 0.5 dB per hop
//       (2 dB/s) and restarts on every re-arm, zeroing or restore. After cancellation true quiet is the
//       noise floor; a note is not, however stationary — a −40 dBFS note over a −90 dB residual floor is
//       not quiet for ≈ 19 s. A LOUDER HUM is told apart from a note: the running family's bins up 3 dB in a
//       quiet, trigger-free window AND the coarse estimate within 0.05 Hz of ours → the residual floor
//       restarts and quiet returns (measured: 0.51 s after a +10 dB step). The input floor no longer gates
//       quiet; it keeps active() and the trigger used before the canceller runs. Smoothed powers are seeded
//       from a settled level, never from one 64-sample period (a waveform trough there set every trigger off);
//   (c) no per-sample trigger in the last 50 ms: a sample 12 dB over the floor's sine peak (x² > 32·floor),
//       or — while the canceller runs — a RESIDUAL sample 12 dB over the residual's own 20 ms-smoothed level
//       (e² > 32·max(smoothed e², 1e-4·floor)). The residual is what a pluck shows in first (the hum is
//       already gone from it), so quiet ends on the first few samples of an attack and an attack is never
//       learned at the fast width; a stationary residual (an unconverged model, a ringing ghost) never
//       trips it. The 50 ms latch covers the smoothed power's lag, and a trigger also voids the last window's
//       verdict at once (a window that held a trigger is never quiet). (a) needs hum in the k1/k2 bins — a hum
//       without them never qualifies; (b) keeps a note from becoming quiet; (c) reacts within a period. ACTIVE (a pluck, reported by active()) is the brief's: the 64-sample
//       power above floor + 12 dB for at least 20 ms.
//
// KNOWING THE MAINS. Off / Hz50 / Hz60 / Auto.
//   * AUTO decides the family from the bins on rectangular 0.5 s windows hopped every 0.25 s on the grid: a
//     family must win by more than 6 dB in 3 consecutive QUIET windows (≈ 1 s of hum). Ambiguous → abstain
//     (dry). Decided only in quiet, so it is sticky while playing. To move OFF a decided family, that family
//     must also have LOST (its bins 10 dB under their level when last confirmed) and the windows must be
//     free of triggers: a 60/120 Hz tone over a 50 Hz hum that is still there flips nothing; a real swap
//     re-decides (measured 1.26 s after it).
//   * COARSE frequency: the winning family's stronger bin, phase advance between windows one hop apart
//     (±2 Hz unambiguous), taken only between two quiet windows and only for a family that wins by 6 dB. Auto starts the canceller at it; a fixed
//     mode jumps to it once, the first time it exists ("acquisition").
//   * FINE (an FLL, only in quiet, and only while the running family is the one the last window's bins
//     heard win, while its coarse estimate agrees with ours within 0.05 Hz — a note at a bin drags the coarse
//     estimate off and so holds the FLL — and not in the first hop after fresh weights; NOT gated by the
//     residual floor, since a residual grown by the FLL's own lag must not hold the FLL that cures it; a
//     fixed mode on the other family's hum would otherwise chase beats): the weights' own rotation, Δf = Σ k·Im(w_k·conj w_k') / (2πT·Σ k²|w_k|²)
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
// harmonic, decaying with τ — an LTI fact of any notch of this width, not a defect of the learner. The ghost
// is RESIDUAL, so quiet (b) does not return while it rings, and it decays at the playing width's τ 0.64 s:
// NOTE SAFETY WINS over a faster clean-up (a ghost and a soft note's own tail look the same to the stage).
// Measured (48 kHz, a 330 Hz tone 20 Hz from k7 = 350): it starts at −38.7 dB re the tone (the LTI value is
// −38.1) and is below −60 dB re the tone 1.56 s after the stop (1.62 s with one width) — it was 0.70 s when
// quiet admitted the ghost, and with it every soft note (HumRemoverTrackingTests).
//
// WHILE PLAYING THE FREQUENCY IS HELD. The FLL runs only in quiet, so a long stretch of playing with no quiet
// gap runs at the frequency it had when the playing began. If the mains drifts at r Hz/s for T s, harmonic k
// is off by Δ_k = k·r·T, and the canceller (still adapting at the playing width) leaves it at
// 10·log10(Δ²/(Δ² + (BW/2)²)) re uncancelled. Worked: r = 1 mHz/s for 60 s (0.06 Hz, a large excursion) leaves
// k = 1 at −12.6 dB, k = 4 at −3.2 dB, k = 10 at −0.7 dB — the upper harmonics are back. The first quiet gap
// re-locks it: once one 0.5 s window is clear (0.5–0.75 s after the playing stops) the FLL closes at up to
// 0.05 Hz/s, so 0.06 Hz takes about 1.2 s of quiet; the coarse estimate is taken once per arm, so a drift
// accumulated while playing is closed by the FLL alone.
//
// LEARN (beginLearn / learnStatus / endLearn): a short, deliberate quiet — the player mutes the strings and
// asks the stage to learn. The model is snapshotted, ramped out, zeroed, and relearnt at the learn width
// with the family Auto's windows hear (in any mode); the correction stays live. Timeline from beginLearn:
// the family is decided in 0.5–0.75 s (3 windows; at once if the running family is confirmed), the weights
// converge, and the FLOOR WINDOW opens at max(1.0 s, arm + 0.75 s) — 0.75 s is 7 τ at the 3 Hz learn width,
// −61 dB of the initial error; a cold start arms at ≈ 1.0 s, so its window opens at ≈ 1.75 s. The report's
// floorAfterDb is the peak |out| over the window — what a gate threshold should sit above. 2 s is the
// useful minimum. Anything a per-sample
// trigger makes it DISTURBED and freezes adaptation from that sample — and so does a floor window that is not
// a floor: a hop whose lowest residual is 6 dB over the lowest hop of the window (a swell, a note creeping in
// under every trigger). A LEARN pressed after the rig changed (more noise, louder hum) is not Disturbed: the
// rule compares the window with itself, not with the stage's past. A note ringing through the whole learn
// that falls (decays) is not refused; floorAfterDb, the window's peak, reports it — high, so a gate set from
// it errs safe. A learn that decided no family, or learned a hum below kNothingDb (−90 dB), is NOTHING HEARD. Only a
// committed, Learned learn is kept, and the report says whether it was (kept: false for a learn of the other
// family in a fixed mode, or in Off); anything else restores the snapshot bit for bit through a ramp. A floor
// window shorter than kLearnFloorMinSec (one hop) is TooShort. A LEARN begun while a cancelled one still ramps
// out keeps the original snapshot. Before the input floor exists (a cold learn's first 0.25 s) the trigger
// works from the settled input level, so an early pluck is Disturbed too. endLearn with no learn running
// says NotLearning; prepare()/reset() are restarts and cancel a learn. THREADING: setParams, beginLearn and
// endLearn belong on the process thread, or strictly between process calls — nothing here is atomic.
//
// TELEMETRY (process thread): detected(), mainsHz(), oscillatorHz(), locked() (the FLL's last steps under
// 5 mHz for a hop, or no tracking asked for), humDb() (10·log10 Σ|w_k|², max over channels), harmonicDb(k),
// quiet(), active().
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
//     measurement. Over a run of exact zeros the correction FADES OUT (a raised cosine, one step per zero,
//     64 steps), so from zero #64 on it is MUTED and digital black in is bit-exact black out; weights, floor
//     and decisions hold. On the way back it fades in the same way over the first 64 real samples (from
//     wherever a shorter run left it). Neither edge steps: the correction's second difference stays under
//     −60 dBFS on a −30 dBFS hum at a region end and a region start, abrupt or faded (measured ≈ −70).
//     Price: the first 64 samples after black are cancelled only partly.
//   * A non-finite input sample passes through untouched and adapts / feeds nothing. Law 8/8a: on the grid
//     each weight is flushed (|w| < 1e-15 → 0), healed (non-finite → 0) and bounded (beyond kWeightCeiling,
//     +60 dBFS, it is not a hum model → 0); the non-finite half also runs at the end of every call.
//   * Law 11a: a channel at index >= nch STOPPED for that call. Its weights are a MODEL of the hum on that
//     input, not audio memory, so they are KEPT; its zero-run counter and black fade advance (it is receiving
//     silence), so it returns faded in, as from any black.
//   * Law 11c: a clock-only call (nch == 0) is digital silence for every shared clock — the grid, the mains
//     and bin phases, the hop and window counters, the ramp and a pending re-arm advance exactly as through
//     n zeros at a live width; weights and detector statistics hold, as they do through zeros.
//
// WHAT IT CANNOT DO: remove hum that is made after this stage's input (e.g. inside a physical amplifier
// downstream), harmonics above the active set, or hum whose phase moves faster than τ lets the weights follow.
//
// RT-safe: no allocation, lock, IO or throw anywhere after construction — there is no heap at all. Footprint
// sizeof(HumRemover): weights, their previous-period copy and the LEARN snapshot, 6·kMaxChannels·MaxHarmonics
// doubles (48 KiB at 64, 18 KiB at 24), plus ≈ 7 KiB of detector state and references.
//==============================================================================

enum class Mains { Off, Hz50, Hz60, Auto };

struct HumRemoverParams
{
    Mains  mains            = Mains::Auto;
    int    harmonics        = 24;            // [1, MaxHarmonics]; a harmonic that could pass maxHz or 0.9·Nyquist at the top of the lock range is skipped
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
    static constexpr double kResidualRise   = 1.1220184543019633;   // the residual floor rises <= 0.5 dB per hop (2 dB/s)

    struct Weight { double re, im; };

    static constexpr double kSilentDb          = -200.0;   // what an empty level reads as
    static constexpr double kNothingDb         = -90.0;    // a learned hum below this is "nothing heard"
    static constexpr double kLearnSettleSec    = 1.0;      // the floor window opens no earlier than this...
    static constexpr double kLearnConvergeSec  = 0.75;     // ...nor earlier than this after the canceller armed (7 τ at 3 Hz)
    static constexpr double kLearnFloorMinSec  = 0.25;     // a floor window shorter than one hop is TooShort

    enum class LearnStatus  { Idle, Running, Disturbed, NothingHeard };
    enum class LearnVerdict { Learned, Disturbed, NothingHeard, TooShort, NotLearning };
    struct LearnReport
    {
        Mains        mains        = Mains::Off;   // the family the learn heard (Off: none)
        double       hz           = 0.0;          // the tracker's estimate at the end
        double       humDb        = kSilentDb;    // the learned hum: 10·log10 Σ|w_k|², dB re a full-scale sine
        double       harmonicDb[MaxHarmonics] {}; // each harmonic's amplitude, dBFS (index k-1)
        double       floorAfterDb = kSilentDb;    // peak |out| over the settled window, dBFS
        double       settledAtMs  = -1.0;         // when the settled window opened, ms after beginLearn (-1: never)
        LearnVerdict verdict      = LearnVerdict::NotLearning;
        bool         kept         = false;        // the stage now runs the learned model (commit, Learned, the mode keeps the family)
    };

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
        settlePeriods_ = (int) std::llround (3.0 * kSmoothSeconds * fs_ / (double) kPeriod) + 1;
        for (int i = 0; i <= kBlackRun; ++i)                  // the black fade: a raised cosine over 64 samples
            blackFade_[i] = 0.5 - 0.5 * core::det::cos (core::kPi * (double) i / (double) kBlackRun);
        blackFade_[0] = 0.0; blackFade_[kBlackRun] = 1.0;
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
        for (int c = 0; c < core::kMaxChannels; ++c) { zeros_[c] = 0; blackGain_[c] = kBlackRun; }
        learning_ = false; swapPhase_ = 0; learnFamily_ = 0; bkLive_ = false; freqPending_ = false;
        binsValid_ = false; binStep_ = 0;                      // every lazily-anchored phasor re-anchors (law 11a)
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
        if (running_) { fHat_ = clampToLock (fHat_); freqPending_ = true; }   // takes effect on the next boundary
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
        {
            zeros_[c]     = (n >= kBlackRun - zeros_[c]) ? kBlackRun : zeros_[c] + n;
            blackGain_[c] = (n >= blackGain_[c]) ? 0 : blackGain_[c] - n;    // its correction faded with the zeros
        }

        for (int done = 0; done < n;)
        {
            const int seg = grid_.segment (n - done);
            if (nch > 0 && running_)      runSegment (io, nch, done, seg);
            else if (nch > 0 && detOn())  detectSegment (io, nch, done, seg);
            else                          clockSegment (seg);
            if (learning_) learnElapsed_ += seg;
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

    // Telemetry (process thread). Levels are the model's, max over channels.
    bool   locked() const noexcept { return running_ && (params_.lockRangeHz <= 0.0 || lockedPeriods_ >= hopPeriods_); }
    double humDb() const noexcept
    {
        double m = 0.0;
        for (int c = 0; c < maxCh_; ++c)
        {
            double sum = 0.0;
            for (int k = 0; k < nActive_; ++k) { const double a = wRe_[c][k] * wRe_[c][k]; const double b = wIm_[c][k] * wIm_[c][k]; sum += a; sum += b; }
            if (sum > m) m = sum;
        }
        return m > 0.0 ? 10.0 * std::log10 (m) : kSilentDb;
    }
    double harmonicDb (int k) const noexcept
    {
        if (k < 1 || k > nActive_) return kSilentDb;
        double m = 0.0;
        for (int c = 0; c < maxCh_; ++c) { const double a = std::hypot (wRe_[c][k - 1], wIm_[c][k - 1]); if (a > m) m = a; }
        return m > 0.0 ? 20.0 * std::log10 (m) : kSilentDb;
    }

    // ---- LEARN ---------------------------------------------------------------------------------------
    // beginLearn: snapshot the model, ramp the correction out (20 ms, nothing adapts meanwhile), zero the
    // weights, and learn at the learn width with the family decided by Auto's windows whatever the mode
    // says; the correction stays live. A per-sample trigger or a period 12 dB over the floor makes it
    // DISTURBED and freezes adaptation from that sample on. endLearn(commit) reports; it keeps the learned
    // model only for commit && Learned, otherwise it ramps out, restores the snapshot bit for bit (weights,
    // frequency, family, Auto's decision) and ramps it back in, adapting nothing until done (swapping()).
    // In a fixed mode a committed learn of the other family re-arms on the mode's family.
    void beginLearn() noexcept
    {
        if (! prepared_ || learning_) return;
        if (swapPhase_ == 1 && swapAction_ == 1)
        {
            // a cancelled learn is still ramping out to its restore: the snapshot it would restore IS the model
            // from before any learn, so keep it, and let the pending swap zero for this learn instead
            learning_ = true; learnStatus_ = LearnStatus::Running; learnElapsed_ = 0; learnFamily_ = 0;
            learnArmAt_ = -1; peakOut_ = 0.0; streak_ = 0; streakFamily_ = 0; lwMin_ = kNone;
            swapAction_ = 0;
            return;
        }
        for (int c = 0; c < core::kMaxChannels; ++c)
            for (int k = 0; k < MaxHarmonics; ++k) { bkRe_[c][k] = wRe_[c][k]; bkIm_[c][k] = wIm_[c][k]; }
        bkFamily_ = family_; bkRunning_ = running_; bkAcquired_ = acquired_; bkAuto_ = autoFamily_;
        bkLocked_ = lockedPeriods_; bkHz_ = fHat_; bkActive_ = nActive_;
        bkInc_ = inc_; bkAnchor_ = anchorPhase_; bkLive_ = true;   // the frame the snapshot's weights belong to
        learning_ = true; learnStatus_ = LearnStatus::Running; learnElapsed_ = 0; learnFamily_ = 0;
        learnArmAt_ = -1; peakOut_ = 0.0; streak_ = 0; streakFamily_ = 0; lwMin_ = kNone;
        if (running_ && ramp_.getCurrentValue() > 0.0f) { swapPhase_ = 1; swapAction_ = 0; ramp_.setTargetValue (0.0f); }
        else zeroForLearn();
    }
    bool learning() const noexcept { return learning_; }
    bool swapping() const noexcept { return swapPhase_ != 0; }
    LearnStatus learnStatus() const noexcept
    {
        if (! learning_) return LearnStatus::Idle;
        if (learnStatus_ == LearnStatus::Disturbed) return LearnStatus::Disturbed;
        if ((double) learnElapsed_ >= kLearnSettleSec * fs_ && (learnFamily_ == 0 || humDb() < kNothingDb)) return LearnStatus::NothingHeard;
        return LearnStatus::Running;
    }
    LearnReport endLearn (bool commit) noexcept
    {
        LearnReport r;
        if (! learning_) return r;                             // NotLearning: there was nothing to end
        const LearnStatus st = learnStatus();
        r.mains = learnFamily_ == 50 ? Mains::Hz50 : (learnFamily_ == 60 ? Mains::Hz60 : Mains::Off);
        r.hz    = running_ ? fHat_ : 0.0;
        r.humDb = humDb();
        for (int k = 0; k < MaxHarmonics; ++k) r.harmonicDb[k] = harmonicDb (k + 1);
        const long long ss = settleStart();
        const bool seen = learnArmAt_ >= 0 && learnElapsed_ >= ss + (long long) std::llround (kLearnFloorMinSec * fs_);
        r.settledAtMs  = seen ? 1000.0 * (double) ss / fs_ : -1.0;
        r.floorAfterDb = seen && peakOut_ > 0.0 ? 20.0 * std::log10 (peakOut_) : kSilentDb;
        r.verdict = st == LearnStatus::Disturbed ? LearnVerdict::Disturbed
                  : (st == LearnStatus::NothingHeard || learnFamily_ == 0 || r.humDb < kNothingDb) ? LearnVerdict::NothingHeard
                  : ! seen ? LearnVerdict::TooShort : LearnVerdict::Learned;
        learning_ = false;
        const int modeFamily = params_.mains == Mains::Hz50 ? 50 : (params_.mains == Mains::Hz60 ? 60 : (params_.mains == Mains::Auto ? learnFamily_ : 0));
        r.kept = commit && r.verdict == LearnVerdict::Learned && modeFamily == learnFamily_;
        if (commit && r.verdict == LearnVerdict::Learned)
        {
            if (params_.mains == Mains::Auto) autoFamily_ = learnFamily_;
            bkLive_ = false;
            if (swapPhase_ == 0) retarget();
        }
        else if (running_ && ramp_.getCurrentValue() > 0.0f) { swapPhase_ = 1; swapAction_ = 1; pending_ = false; ramp_.setTargetValue (0.0f); }
        else restoreBackup();
        return r;
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
    bool detOn() const noexcept { return params_.mains != Mains::Off || running_ || learning_; }

    int wantedFamily() const noexcept
    {
        if (learning_ && learnFamily_ != 0) return learnFamily_;
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
        const double f0  = (double) family + params_.lockRangeHz;   // the cut-off holds wherever the tracker goes
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
        if (swapPhase_ != 0) return;                           // a LEARN swap in flight reconciles when it ends
        const int fam = wantedFamily();
        pending_ = ! (fam == family_ && activeFor (fam) == nActive_);
        ramp_.setTargetValue (pending_ || ! running_ ? 0.0f : 1.0f);
    }

    // The oscillator at f: only the increment and the per-sample rotations change, never the phase. Called
    // ONLY on a grid boundary (after the anchor advance) or on a restart — a period that ran p samples at one
    // increment must be advanced at that increment; a change asked for between calls waits (freqPending_).
    void setFrequency (double f) noexcept
    {
        fHat_ = f;
        setIncrement (f > 0.0 ? incFor (f) : 0u);
    }

    void setIncrement (std::uint32_t ni) noexcept
    {
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
        rotsValid_ = false; freqPending_ = false;
        setFrequency (running_ ? f : 0.0);
        for (int c = 0; c < core::kMaxChannels; ++c)
            for (int k = 0; k < MaxHarmonics; ++k) { wRe_[c][k] = 0.0; wIm_[c][k] = 0.0; pRe_[c][k] = 0.0; pIm_[c][k] = 0.0; }
        refsValid_ = false;
        pending_   = false;
        periodSumE_ = 0.0; seedResidual();
        lockedPeriods_ = 0;
        resetResidualFloor();
        if (learning_) learnArmAt_ = running_ ? learnElapsed_ : -1;
        ramp_.setCurrentAndTargetValue (running_ ? 1.0f : 0.0f);
    }

    long long settleStart() const noexcept
    {
        const long long a = (long long) std::llround (kLearnSettleSec * fs_);
        const long long b = learnArmAt_ + (long long) std::llround (kLearnConvergeSec * fs_);
        return a > b ? a : b;
    }

    void zeroForLearn() noexcept
    {
        for (int c = 0; c < core::kMaxChannels; ++c)
            for (int k = 0; k < MaxHarmonics; ++k) { wRe_[c][k] = 0.0; wIm_[c][k] = 0.0; pRe_[c][k] = 0.0; pIm_[c][k] = 0.0; }
        seedResidual(); hitHold_ = 0; lockedPeriods_ = 0;
        resetResidualFloor();
        learnArmAt_ = running_ ? learnElapsed_ : -1;
        swapPhase_ = 0;
        if (running_) ramp_.setCurrentAndTargetValue (1.0f);
    }

    void restoreBackup() noexcept
    {
        for (int c = 0; c < core::kMaxChannels; ++c)
            for (int k = 0; k < MaxHarmonics; ++k) { wRe_[c][k] = bkRe_[c][k]; wIm_[c][k] = bkIm_[c][k]; pRe_[c][k] = bkRe_[c][k]; pIm_[c][k] = bkIm_[c][k]; }
        family_ = bkFamily_; nActive_ = bkActive_; running_ = bkRunning_; acquired_ = bkAcquired_; autoFamily_ = bkAuto_;
        lockedPeriods_ = bkLocked_;
        // the snapshot's own phase frame: its anchor advanced at its own increment through the learn, so the
        // weights come back paired with the references they were learned against (exact, uint32)
        anchorPhase_ = bkAnchor_; fHat_ = bkHz_; freqPending_ = false; bkLive_ = false;
        streak_ = 0; streakFamily_ = 0;                        // Auto re-decides only on windows after the restore
        rotsValid_ = false;
        setIncrement (bkInc_);
        refsValid_ = false; seedResidual(); pending_ = false;
        resetResidualFloor();
        if (running_) { ramp_.setCurrentAndTargetValue (0.0f); ramp_.setTargetValue (1.0f); swapPhase_ = 2; }
        else          { ramp_.setCurrentAndTargetValue (0.0f); swapPhase_ = 0; retarget(); }
    }

    // After cancellation true quiet is the noise floor; a note is not. Within 12 dB of the residual floor.
    bool residualQuiet() const noexcept { return ! running_ || ! rfValid_ || ! psEValid_ || ! (psE_ > kQuietRatio * rf_); }

    // With zero (or freshly restored) weights the residual IS the input, so its level starts at the input's
    // settled level, not at one period's power — a trough there would set every trigger off on the next peak.
    void seedResidual() noexcept
    {
        fllHold_ = hopPeriods_ * kPeriod;   // fresh weights are small and moving: the FLL waits one hop for them
        psE_ = ps_; psEValid_ = psValid_; psEAge_ = settlePeriods_;
        if (psEValid_ && floorValid_) { const double lo = 1e-4 * floor_; thrE2_ = kHitRatio * (psE_ > lo ? psE_ : lo); }
        else thrE2_ = thr2_;
    }

    void resetResidualFloor() noexcept { rfCount_ = 0; rfPos_ = 0; rfValid_ = false; rf_ = 0.0; hopMinE_ = kNone; }

    void resetDetector() noexcept
    {
        famRef_[0] = famRef_[1] = 0.0; hitInHop_ = false; hitPrevHop_ = false; lwMin_ = kNone; coarseOk_ = false; psAge_ = 0; psEAge_ = 0;
        resetResidualFloor();
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
            if (! running_ && x2 > thr2_) hit_ = true;       // while the canceller runs, the residual's trigger rules
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
        const int       p0  = grid_.phase();
        const bool      det = detOn();
        const long long ss  = learning_ && learnArmAt_ >= 0 ? settleStart() : -1;
        for (int i = 0; i < seg; ++i)
        {
            const double g = (double) ramp_.getNextValue();
            const bool   floorWindow = ss >= 0 && learnElapsed_ + i >= ss;
            if (det) observe (io, nch, off + i, p0 + i);
            for (int c = 0; c < nch; ++c)
            {
                core::Sample& s = io[c][off + i];
                const core::Sample x = s;
                const bool zero = core::exactlyEqual (x, 0.0f);
                if (zero)
                {
                    if (zeros_[c] < kBlackRun) ++zeros_[c];
                    if (blackGain_[c] > 0) --blackGain_[c];           // the correction fades out over the run...
                    if (zeros_[c] >= kBlackRun) continue;             // ...and from zero #64: muted, bit-exact
                }
                else
                {
                    zeros_[c] = 0;
                    if (! std::isfinite (x)) continue;               // passes untouched, adapts nothing
                    if (blackGain_[c] < kBlackRun) ++blackGain_[c];   // back from (partial) black: fades in
                }

                syncRefs (p0 + i);
                const double corr = correction (c);
                const double gb   = g * blackFade_[blackGain_[c]];
                const double gc   = gb * corr;
                s = (core::Sample) ((double) x - gc);
                if (floorWindow) { const double a = std::fabs ((double) s); if (a > peakOut_) peakOut_ = a; }
                if (zero) continue;                                  // a zero is no measurement

                const double err = (double) x - corr;
                if (det)
                {
                    const double e2 = err * err;
                    periodSumE_ += e2;
                    if (e2 > thrE2_) hit_ = true;
                }
                if (learning_ && hit_) learnStatus_ = LearnStatus::Disturbed;   // frozen from this very sample
                chanLive_[c] = true;
                const double mu    = swapPhase_ != 0 ? 0.0
                                   : learning_ ? (learnStatus_ == LearnStatus::Disturbed ? 0.0 : muFast_)
                                   : ((quietPeriod_ && ! hit_) ? muFast_ : muSlow_);
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
        if (bkLive_) bkAnchor_ += P * bkInc_;                       // the LEARN snapshot's frame, at its own rate
        for (int b = 0; b < 4; ++b) binAnchor_[b] += P * binInc_[b];
        refsValid_ = false;
        binsValid_ = false;
        if (detOn()) endPeriod();
        for (int c = 0; c < core::kMaxChannels; ++c) chanLive_[c] = false;
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
        if (swapPhase_ == 1 && ! ramp_.isSmoothing() && core::exactlyEqual (ramp_.getCurrentValue(), 0.0f))
        {
            if (swapAction_ == 0) zeroForLearn(); else restoreBackup();
        }
        else if (swapPhase_ == 2 && ! ramp_.isSmoothing()) { swapPhase_ = 0; retarget(); }
        else if (swapPhase_ == 0 && pending_ && ! ramp_.isSmoothing() && core::exactlyEqual (ramp_.getCurrentValue(), 0.0f))
            applyWanted();
        if (freqPending_) { freqPending_ = false; if (running_) setFrequency (fHat_); }
    }

    void endPeriod() noexcept
    {
        // the FLL's own gate: the window was hum, no trigger, and the bins' coarse estimate agrees with ours —
        // not the residual criterion (a residual grown by the FLL's own lag must not hold the FLL that cures it)
        const bool   wasQuiet = quietWin_ && floorValid_ && hitHold_ <= 0 && ! hit_ && coarseOk_ && fllHold_ <= 0;
        if (fllHold_ > 0) fllHold_ -= kPeriod;
        const double sum      = periodSum_;
        const double sumE     = periodSumE_;
        periodSum_  = 0.0;
        periodSumE_ = 0.0;
        hopPow_   += sum;
        if (sum > 0.0)                                         // an all-zero period (black, a pause) holds
        {
            const double p = sum / (double) kPeriod;
            if (! psValid_) { ps_ = p; psValid_ = true; psAge_ = 0; }
            else
            {
                const double d = p - ps_;
                const double s = smoothA_ * d;
                ps_ += s;
            }
            if (psAge_ < settlePeriods_) ++psAge_;             // one period's power is not a level: settle first
            const bool above = floorValid_ && p > kQuietRatio * floor_;
            activeRun_ = above ? (activeRun_ < activeSamples_ ? activeRun_ + kPeriod : activeRun_) : 0;
            if (psAge_ >= settlePeriods_ && ps_ < hopMin_) hopMin_ = ps_;
            if (running_ && wasQuiet && winner_ == family_ && swapPhase_ == 0 && params_.lockRangeHz > 0.0) trackFrequency();
        }
        if (sumE > 0.0)                                        // the residual's smoothed level (running only)
        {
            const double pe = sumE / (double) kPeriod;
            if (! psEValid_) { psE_ = pe; psEValid_ = true; psEAge_ = 0; }
            else
            {
                const double d = pe - psE_;
                const double s = smoothA_ * d;
                psE_ += s;
            }
            if (psEAge_ < settlePeriods_) ++psEAge_;
            else if (psE_ < hopMinE_) hopMinE_ = psE_;
        }
        if (hit_) { hitInHop_ = true; quietWin_ = false; }   // a trigger voids the last window's verdict at once
        if (hit_ && learning_) learnStatus_ = LearnStatus::Disturbed;   // also before the canceller runs (cold learn)
        if (++hopCount_ >= hopPeriods_) endHop();
        if (hit_) hitHold_ = holdSamples_;
        else if (hitHold_ > 0) hitHold_ -= kPeriod;
        quietPeriod_ = quietWin_ && floorValid_ && hitHold_ <= 0 && residualQuiet();
        hit_  = false;
        thr2_ = floorValid_ ? kHitRatio * floor_ : (psValid_ && psAge_ >= settlePeriods_ ? kHitRatio * ps_ : kNone);
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
            for (int k = 0; k < nActive_ && chanLive_[c]; ++k)   // only channels that adapted this period
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
        lockedPeriods_ = std::fabs (delta) < 0.005 ? (lockedPeriods_ < hopPeriods_ ? lockedPeriods_ + 1 : lockedPeriods_) : 0;
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
            // A window counts for a decision only if it was quiet, free of triggers, and — to move off a decided
            // family — that family has LOST (its bins 10 dB under their level when it was last confirmed): a tone
            // at 60 Hz over a 50 Hz hum that is still there wins nothing.
            const bool   trigFree = ! hitInHop_ && ! hitPrevHop_;
            const double famOf[2] = { fam50, fam60 };
            const int    decided  = learning_ ? learnFamily_ : autoFamily_;
            const int    ds       = decided == 60 ? 1 : 0;
            const bool   lost     = decided != 0 && famRef_[ds] > 0.0 && famOf[ds] < 0.1 * famRef_[ds];
            if (quiet && winner != 0 && (decided == 0 || winner == decided || (lost && trigFree)))
                 { streak_ = winner == streakFamily_ ? streak_ + 1 : 1; streakFamily_ = winner; }
            else { streak_ = 0; streakFamily_ = 0; }
            // A LOUDER HUM is not a note: the running family's bins up 3 dB, at the frequency we track — the
            // residual floor restarts, so quiet returns while the canceller re-converges at the learn width.
            if (running_ && quiet && trigFree && winner == family_)
            {
                const int fs2 = family_ == 60 ? 1 : 0;
                if (famRef_[fs2] > 0.0 && famOf[fs2] > 2.0 * famRef_[fs2] && coarseValid_[fs2] && std::fabs (coarse_[fs2] - fHat_) < 0.05)
                    resetResidualFloor();
            }
            if (quiet && trigFree && winner != 0) famRef_[winner == 60 ? 1 : 0] = famOf[winner == 60 ? 1 : 0];
            if (learning_)                                      // a learn decides with Auto's windows, whatever the mode
            {
                if (streak_ >= kWinStreak && winner != learnFamily_) { learnFamily_ = winner; retarget(); }
            }
            else if (streak_ >= kWinStreak && winner != autoFamily_)
            {
                autoFamily_ = winner;
                if (params_.mains == Mains::Auto) retarget();
            }
            {
                const int cs = family_ == 60 ? 1 : 0;
                coarseOk_ = running_ && coarseValid_[cs] && std::fabs (coarse_[cs] - fHat_) < 0.05;
            }
            if (running_ && ! acquired_ && params_.lockRangeHz > 0.0 && coarseValid_[family_ == 60 ? 1 : 0])
            {
                setFrequency (clampToLock (coarse_[family_ == 60 ? 1 : 0]));
                acquired_ = true;
                coarseOk_ = true;
            }
            quietWin_ = quiet && trigFree;                     // a window that held a trigger is not quiet
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
        // LEARN's floor window must be a floor: a hop whose LOWEST residual is 6 dB over the lowest hop seen in
        // the window (a swell, a note creeping in under every trigger) disturbs the learn. Hop minima, so a
        // spike cannot trip it and a converging (falling) residual never does.
        if (learning_ && learnArmAt_ >= 0 && learnElapsed_ >= settleStart() && hopMinE_ < kNone)
        {
            const double lo = floorValid_ ? 1e-6 * floor_ : 0.0;
            if (lwMin_ < kNone && hopMinE_ > 4.0 * (lwMin_ > lo ? lwMin_ : lo)) learnStatus_ = LearnStatus::Disturbed;
            if (hopMinE_ < lwMin_) lwMin_ = hopMinE_;
        }
        if (hopMinE_ < kNone)                                  // the residual floor: same, its rise capped per hop
        {
            rfRing_[rfPos_] = hopMinE_;
            rfPos_ = (rfPos_ + 1) % kFloorHops;
            if (rfCount_ < kFloorHops) ++rfCount_;
            double m = rfRing_[0];
            for (int i = 1; i < rfCount_; ++i) if (rfRing_[i] < m) m = rfRing_[i];
            const double capped = rf_ * kResidualRise;
            rf_ = (rfValid_ && m > capped) ? capped : m;
            rfValid_ = true;
        }
        hopMinE_ = kNone;
        hitPrevHop_ = hitInHop_;
        hitInHop_   = false;
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
    int    blackGain_[core::kMaxChannels] {};   // 0..kBlackRun: the black fade's position per channel
    bool   chanLive_[core::kMaxChannels] {};    // adapted during this grid period (the FLL sums only these)
    double blackFade_[kBlackRun + 1] {};

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
    int    lockedPeriods_ = 0;
    std::uint32_t bkInc_ = 0, bkAnchor_ = 0;
    double rfRing_[kFloorHops] {}, rf_ = 0.0, hopMinE_ = kNone, famRef_[2] {}, lwMin_ = kNone;
    int    rfCount_ = 0, rfPos_ = 0, psAge_ = 0, psEAge_ = 0, settlePeriods_ = 46;
    bool   coarseOk_ = false;
    int    fllHold_ = 0;
    bool   rfValid_ = false, hitInHop_ = false, hitPrevHop_ = false;
    bool   bkLive_ = false, freqPending_ = false;

    // ---- LEARN ----
    double      bkRe_[core::kMaxChannels][MaxHarmonics] {};   // the snapshot endLearn(false) restores
    double      bkIm_[core::kMaxChannels][MaxHarmonics] {};
    double      bkHz_ = 0.0, peakOut_ = 0.0;
    int         bkFamily_ = 0, bkActive_ = 0, bkAuto_ = 0, bkLocked_ = 0, learnFamily_ = 0, swapPhase_ = 0, swapAction_ = 0;
    bool        bkRunning_ = false, bkAcquired_ = false, learning_ = false;
    long long   learnElapsed_ = 0, learnArmAt_ = -1;
    LearnStatus learnStatus_ = LearnStatus::Idle;
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

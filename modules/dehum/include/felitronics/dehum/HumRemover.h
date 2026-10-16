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
// μ is set from `bandwidthHz` as μ = 2π·BW/fs, so the width is a CONSTANT ABSOLUTE number of hertz at every
// harmonic and every rate — not a Q. That is the point: a note's partial 0.6 Hz from 660 Hz loses the same
// fraction of a dB as one 0.6 Hz from 60, and widening the upper notches (a constant Q) would eat every E's
// 660 Hz partial. The default 0.5 Hz costs ≈ −0.97 dB on a partial 0.5 Hz from a harmonic and −0.03 dB at
// 3 Hz. The nulls sit EXACTLY on k·f_osc whatever the other harmonics do (each notch's zeros are its own).
// The price of a narrow notch is that the frequency must be known to millihertz; until the tracker lands,
// 50 and 60 run at the nominal frequency, rounded to the oscillator's resolution (fs/2³² Hz).
//
// MODES. Off / Hz50 / Hz60 / Auto. AUTO IS DRY IN THIS VERSION — the mains detector that decides it is not
// here yet, and a stage that guessed would inject the wrong comb; Auto passes the input through untouched,
// exactly like Off. `learnBandwidthHz` is stored for the quiet-time learner and is not used yet.
// A change of the mode or of the active harmonic set ramps the correction out over 20 ms (a
// core::LinearSmoother), then — at the next law-8a grid boundary — re-arms: the weights are zeroed and the
// new set applies. Zero weights make a zero correction, so the gain returns to 1 without a step. A change of
// `bandwidthHz` alone applies at once: μ moves no weight, so it cannot step the output.
//
// NUMERICS (law 3 carve-out, stated): weights and references in DOUBLE, samples in float. Float weights
// stall where an increment drops below half an ulp of the weight (the design study read −63 dB re the hum at
// 48 kHz, −51 at 192 kHz). The mains phase is a uint32 advanced by round(f·2³²/fs); every 64 samples (core::StateGrid)
// each reference is RE-ANCHORED from that phase through core::det sin/cos — k·phase wraps exactly in uint32
// arithmetic — and between anchors it advances by one complex rotation in double. The anchors are at fixed
// audio-time positions, so the output is BIT-IDENTICAL under any re-slicing of the stream (law 8a).
// LAW 10: the weight update is spelled so nothing in it can contract under ANY mode: the product μ·e·r is
// ROUNDED TO FLOAT and stored, then added to the double weight (detail::adapt). An FMA would have to skip
// that rounding, which no contraction mode may do. The correction sum and the rotation are written one
// operation per statement, so the repository's -ffp-contract=on fuses nothing there either.
//
// SILENCE, POISON, STOPPED CHANNELS:
//   * An exact-zero input sample ADAPTS NOTHING — digital zero is the absence of a measurement, not a
//     measurement of no hum. After kBlackRun (64) consecutive exact zeros the correction is MUTED too, so
//     digital black in is bit-exact black out, and the weights hold. The first 63 zeros still carry the
//     correction (the anti-hum of the stream that just stopped), then it cuts.
//   * A non-finite input sample passes through untouched and adapts nothing, so a NaN/Inf never reaches
//     the model. Law 8/8a: on the grid each weight is flushed (|w| < 1e-15 → 0), healed (non-finite → 0) and
//     bounded (beyond kWeightCeiling, i.e. +60 dBFS, it is not a hum model → 0); the non-finite half also
//     runs at the end of every call. The weights decay only by adapting to real, non-black quiet.
//   * Law 11a: a channel at index >= nch STOPPED for that call. Its weights are a MODEL of the hum on that
//     input, not audio memory, so they are KEPT; the channel is receiving digital silence, so its zero-run
//     counter advances (after 64 stopped samples it returns muted until a non-zero sample arrives).
//   * Law 11c: a clock-only call (nch == 0) is digital silence for the shared clocks — the grid, the mains
//     phase, the 20 ms ramp and a pending re-arm advance exactly as they would through n zeros at a live
//     width; the weights hold, as they do through zeros. The two are bit-identical.
//
// WHAT IT CANNOT DO: remove hum that is made after this stage's input (e.g. inside a physical amplifier
// downstream), harmonics above the active set, or hum whose phase moves faster than τ lets the weights follow.
//
// RT-safe: no allocation, lock, IO or throw anywhere after construction — there is no heap at all; the
// footprint is sizeof(HumRemover) (≈ 2·16·MaxHarmonics doubles of weights + 4·MaxHarmonics of references).
//==============================================================================

enum class Mains { Off, Hz50, Hz60, Auto };

struct HumRemoverParams
{
    Mains  mains            = Mains::Auto;   // Auto is dry in this version (see the header comment)
    int    harmonics        = 24;            // [1, MaxHarmonics]; harmonics above maxHz or 0.9·Nyquist are skipped
    double bandwidthHz      = 0.5;           // notch −3 dB width while playing; [kMinBandwidthHz, kMaxBandwidthHz]
    double learnBandwidthHz = 3.0;           // the quiet/learn width — stored, not used yet
    double maxHz            = 2000.0;        // highest harmonic frequency processed
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
    static constexpr double kRampSeconds    = 0.020;    // correction ramp-out on a mode / harmonic-set change
    static constexpr double kMinBandwidthHz = 0.01;
    static constexpr double kMaxBandwidthHz = 20.0;
    static constexpr double kWeightCeiling  = 1.0e3;    // +60 dBFS: beyond it a weight is not a hum model

    struct Weight { double re, im; };

    // Law 11b: disarm, validate, write. Refuses fs < core::kMinSampleRate (or non-finite), maxBlock < 1 and
    // channels outside [1, core::kMaxChannels]. A prepare is a restart (law 11a): weights, clocks, ramp.
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
        prepared_ = true;
        reset();
        return true;
    }

    bool isPrepared()  const noexcept { return prepared_; }
    int  maxChannels() const noexcept { return prepared_ ? maxCh_ : 0; }
    int  maxBlock()    const noexcept { return prepared_ ? maxBlock_ : 0; }
    static constexpr int latencySamples() noexcept { return 0; }

    // A stream restart: zero weights, re-anchored clocks (phase 0 at the grid origin), the current params
    // applied at once (no ramp), zero-run counters cleared.
    void reset() noexcept
    {
        grid_.reset();
        anchorPhase_ = 0;
        for (int c = 0; c < core::kMaxChannels; ++c) zeros_[c] = 0;
        applyWanted();
    }

    // Takes effect at once for the widths; a mode or harmonic-set change ramps out and re-arms (see above).
    void setParams (const HumRemoverParams& p) noexcept
    {
        params_.mains            = sanitizeMains (p.mains);
        params_.harmonics        = p.harmonics < 1 ? 1 : (p.harmonics > MaxHarmonics ? MaxHarmonics : p.harmonics);
        params_.bandwidthHz      = sanitizeBw (p.bandwidthHz, params_.bandwidthHz);
        params_.learnBandwidthHz = sanitizeBw (p.learnBandwidthHz, params_.learnBandwidthHz);
        params_.maxHz            = (std::isfinite (p.maxHz) && p.maxHz >= 0.0) ? p.maxHz : params_.maxHz;
        if (! prepared_) return;
        mu_ = 2.0 * core::kPi * params_.bandwidthHz / fs_;
        const int wantedActive = activeFor (params_.mains);
        pending_ = ! (params_.mains == mode_ && wantedActive == nActive_);
        ramp_.setTargetValue (pending_ || ! running_ ? 0.0f : 1.0f);
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
            if (running_ && nch > 0) runSegment (io, nch, done, seg);
            else                     clockSegment (seg);
            if (grid_.advance (seg)) atGridBoundary();
            done += seg;
        }
        healPoison();   // law 8a(a): the poison half keeps its own clock — the end of every call
        return true;
    }

    // ---- inspection (process thread) -------------------------------------------------------------
    bool   running()         const noexcept { return running_; }
    int    activeHarmonics() const noexcept { return nActive_; }
    // The frequency the canceller runs at — the nominal one rounded to fs/2³² — or 0 while dry.
    double mainsHz()         const noexcept { return running_ ? (double) inc_ * fs_ / 4294967296.0 : 0.0; }
    Weight weight (int channel, int harmonic) const noexcept
    {
        if (channel < 0 || channel >= core::kMaxChannels || harmonic < 1 || harmonic > MaxHarmonics) return { 0.0, 0.0 };
        return { wRe_[channel][harmonic - 1], wIm_[channel][harmonic - 1] };
    }

private:
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

    // Auto is dry until the detector exists, so only 50 and 60 have a frequency.
    static double nominalHz (Mains m) noexcept
    {
        return m == Mains::Hz50 ? 50.0 : (m == Mains::Hz60 ? 60.0 : 0.0);
    }

    int activeFor (Mains m) const noexcept
    {
        const double f0 = nominalHz (m);
        if (! (f0 > 0.0)) return 0;
        const double top = params_.maxHz < 0.45 * fs_ ? params_.maxHz : 0.45 * fs_;
        const double fit = std::floor (top / f0);
        const int    k   = fit < (double) params_.harmonics ? (int) fit : params_.harmonics;
        return k < 0 ? 0 : k;
    }

    // Re-arm on the CURRENT params: zero weights, the new frequency and set, the gain at its resting value.
    void applyWanted() noexcept
    {
        mu_       = 2.0 * core::kPi * params_.bandwidthHz / fs_;
        mode_     = params_.mains;
        nActive_  = activeFor (mode_);
        running_  = nActive_ > 0;
        inc_      = running_ ? (std::uint32_t) std::llround (nominalHz (mode_) * 4294967296.0 / fs_) : 0u;
        for (int k = 0; k < nActive_; ++k)
            detail::phasor ((std::uint32_t) (k + 1) * inc_, rotRe_[k], rotIm_[k]);
        for (int c = 0; c < core::kMaxChannels; ++c)
            for (int k = 0; k < MaxHarmonics; ++k) { wRe_[c][k] = 0.0; wIm_[c][k] = 0.0; }
        refsValid_ = false;
        pending_   = false;
        ramp_.setCurrentAndTargetValue (running_ ? 1.0f : 0.0f);
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
            {
                const double re = rRe_[k], im = rIm_[k], c = rotRe_[k], s = rotIm_[k];
                const double a = re * c;
                const double b = im * s;
                const double d = re * s;
                const double e = im * c;
                rRe_[k] = a - b;
                rIm_[k] = d + e;
            }
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

    void runSegment (core::Sample* const* io, int nch, int off, int seg) noexcept
    {
        const int p0 = grid_.phase();
        for (int i = 0; i < seg; ++i)
        {
            const double g = (double) ramp_.getNextValue();
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

                const double err   = (double) x - corr;
                const double muErr = mu_ * err;
                for (int k = 0; k < nActive_; ++k)
                    detail::adapt (wRe_[c][k], wIm_[c][k], muErr, rRe_[k], rIm_[k]);
            }
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
        anchorPhase_ += (std::uint32_t) core::StateGrid::kPeriod * inc_;   // exact, mod 2³²
        refsValid_ = false;
        for (int c = 0; c < maxCh_; ++c)
            for (int k = 0; k < nActive_; ++k)
            {
                double& re = wRe_[c][k];
                double& im = wIm_[c][k];
                if (! (std::fabs (re) <= kWeightCeiling && std::fabs (im) <= kWeightCeiling)) { re = 0.0; im = 0.0; }
                core::flushDenormal (re);
                core::flushDenormal (im);
            }
        if (pending_ && ! ramp_.isSmoothing() && core::exactlyEqual (ramp_.getCurrentValue(), 0.0f))
            applyWanted();
    }

    void healPoison() noexcept
    {
        for (int c = 0; c < maxCh_; ++c)
            for (int k = 0; k < nActive_; ++k)
                if (! (std::isfinite (wRe_[c][k]) && std::isfinite (wIm_[c][k]))) { wRe_[c][k] = 0.0; wIm_[c][k] = 0.0; }
    }

    // Weights per channel, per harmonic (index k-1). References and per-sample rotations are shared.
    double wRe_[core::kMaxChannels][MaxHarmonics] {};
    double wIm_[core::kMaxChannels][MaxHarmonics] {};
    double rRe_[MaxHarmonics] {}, rIm_[MaxHarmonics] {};
    double rotRe_[MaxHarmonics] {}, rotIm_[MaxHarmonics] {};
    int    zeros_[core::kMaxChannels] {};

    HumRemoverParams    params_ {};
    core::StateGrid     grid_ {};
    core::LinearSmoother ramp_ { 0.0f };
    std::uint32_t       anchorPhase_ = 0;   // the mains phase at the start of the current grid period
    std::uint32_t       inc_         = 0;   // per-sample phase increment, round(f·2³²/fs)
    double              fs_          = 48000.0;
    double              mu_          = 0.0;
    Mains               mode_        = Mains::Off;
    int                 nActive_     = 0;
    int                 refStep_     = 0;
    int                 maxCh_       = 0;
    int                 maxBlock_    = 0;
    bool                refsValid_   = false;
    bool                running_     = false;
    bool                pending_     = false;
    bool                prepared_    = false;
};

} // namespace felitronics::dehum

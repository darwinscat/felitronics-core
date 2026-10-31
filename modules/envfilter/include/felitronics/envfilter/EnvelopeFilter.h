// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/core/DetMath.h>
#include <felitronics/core/FlushToZero.h>
#include <felitronics/core/Math.h>
#include <felitronics/core/StateGrid.h>
#include <felitronics/dynamics/EnvelopeFollower.h>
#include <felitronics/eq/Svf.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <limits>

namespace felitronics::envfilter
{

//==============================================================================
// felitronics::envfilter::EnvelopeFilter — the funk envelope filter (the Mu-Tron III family), voiced for bass:
// a resonant filter whose cutoff is driven by how hard the string is hit. The note opens it, its decay closes
// it, and the follower's release is the "wow" between the two.
//
// THE PATH, per sample.
//   · DETECTOR — the loudest |x| of the call's channels (a non-finite sample counts as silence) into a peak
//     dynamics::EnvelopeFollower (attack, release). Linked: one envelope moves both channels' filter.
//   · THE SWEEP — the envelope in dB, pushed by SENSITIVITY, places the cutoff on a log-frequency line between
//     the range's two ends: the sweep is fully open when envelope + sensitivity reaches 0 dBFS and closed
//     kSpanDb below that. UP rests at the low end and a hit opens it; DOWN rests at the high end and a hit
//     pushes it down. The position goes through ONE kSmoothSeconds one-pole on log2(Hz) — the single place
//     every move of the cutoff passes through: the playing, the slider, a range edit, a direction flip — and
//     the filter is redesigned EVERY sample from it, so the cutoff never steps.
//   · THE CLEAN LOWS — a 2-pole high pass (eq::DeterministicSvf) at keepLowsHz gives H, and the lows are
//     x − H: complementary BY CONSTRUCTION, the two add back to x. `keepLows` (0..1) takes that share out of
//     the filter's input and adds it back clean, so the note keeps its weight under a band-pass or a heavy
//     resonance. Not eq::Crossover2, on purpose: an LR4 pair adds up to an ALLPASS that is −1 at the crossover,
//     so a dry/wet blend at mix 0.5 would notch the crossover frequency out of the bass. With x − H the wet's
//     neutral is x itself and the mix crossfades two signals that agree wherever the filter does nothing. The
//     price is a soft split: the clean share rolls off at 6 dB/oct above keepLowsHz, and it lifts the band
//     below it a little. The high pass is CRITICALLY DAMPED (Q 0.5) because that is the flattest complement a
//     2-pole split has: |x − H|² = (1 + 4w²)/(1 + w²)², at most +1.25 dB (at 0.71 × keepLowsHz), where a
//     Butterworth split lifts +2.1 dB — measured in EnvelopeFilterTests.
//   · THE FILTER — a Cytomic TPT state-variable filter, eq::Svf's update line for line, but in a kernel of its
//     own (detail::SweptSvf) for two reasons: ONE pair of integrators gives both taps, so a mode change
//     crossfades the low-pass and band-pass outputs of one state (kModeSeconds) instead of switching filters;
//     and the per-sample redesign needs only g = tan(π·fc/fs) and k = 1/Q, where eq::Svf::setParams designs a
//     gain term too. Both taps are scaled by 1/√max(1, Q), so the resonant peak grows as √Q in either mode
//     and the two modes sit at one level.
//   · THE MIX — out = x + mix·(wet − x), wet = keepLows·lows + filter. Dry is the input itself.
//
// AMOUNT 0 IS NOT THERE. When `mix` stands at 0 — set so, or its glide has landed there — the stage SLEEPS:
// process() returns with the buffer untouched (bit-exact dry, NaN payloads included) and computes nothing.
// Waking up starts it fresh: the follower empty, the cutoff at rest, the mix gliding in from 0. The one-slider
// law that maps a face's slider onto these knobs is EnvelopeFilterLaw.h.
//
// CONTRACT (DSP-ARCHITECTURE.md §2).
//   · prepare(sampleRate, maxBlock, maxChannels): maxChannels 1 or 2, the rate 8 kHz … 768 kHz, maxBlock >= 1
//     (recorded, not a limit: the stage has no block-sized scratch). Anything else is REFUSED, and a refused
//     prepare DISARMS the instance (law 11b). prepare() performs a reset().
//   · process() is [[nodiscard]] bool (law 11): unprepared, a width above the prepared one, a negative extent,
//     or a null buffer with nch > 0 is refused before anything moves; any n >= 0 is processed in full.
//   · A NARROWER call is legal; a channel at an index >= nch with n > 0 has STOPPED (law 11a): its filter and
//     its split are cleared once, on that edge. nch == 0 with n > 0 is a gap: the SHARED ballistics — the
//     follower, the cutoff smoother, the knob glides and the mode door — run through it exactly as through n
//     samples of digital silence (law 11c), and stop computing once they stand still.
//   · reset() is a stream RESTART: every state cleared, the audio-time grid re-anchored, and the next call
//     starts the glides ON their targets — a fresh stage does not fade in.
//   · KNOBS GLIDE. mix, sensitivity, resonance and keepLows glide over kGlideSeconds (one-pole, in double,
//     and they LAND — LAW8-AUDIT.md item 6); the cutoff glides through the smoother above; the mode crossfades
//     over kModeSeconds, counted, so it lands exactly. The attack/release times and keepLowsHz retune at once:
//     they are menu settings, not something a hand rides.
//   · Law 8: the follower, both filters' integrators and the split's are flushed at the house 1e-15 on the
//     core::StateGrid (64 samples of AUDIO time, law 8a), so the output does not depend on how the caller
//     slices the stream. The glides land rather than creep, so none of them reaches the subnormal band.
//   · POISON: a non-finite INPUT sample is not fed to the detector, the split or the filter (they hear silence
//     there); the output sample keeps it — it is the caller's, not this stage's to change. A non-finite filter
//     state at the end of a call (a finite input can still overflow it) restarts that channel from silence.
//     Knob values are bounded before any arithmetic: a non-finite one is ignored, a finite one is clamped.
//   · THREADS: prepare, reset, the setters and process belong to ONE thread. shownCutoffHz is a relaxed atomic
//     that any other thread may read.
//   · Math: core::det (log2, exp2, tan) for the per-sample sweep, so every row maps an envelope to the same
//     cutoff; the follower's two coefficients are dynamics::EnvelopeFollower's own (system exp, per setTimes).
//==============================================================================

enum class Mode      { lowPass, bandPass };   ///< which tap of the swept filter is heard
enum class Direction { up, down };            ///< up: rests low, a hit opens it; down: rests high, a hit closes it

/** Every knob, in the units the stage takes. `EnvelopeFilter::setParams` applies them through the same guarded
    setters a caller may use one by one. */
struct EnvelopeFilterParams
{
    float     sensitivityDb = 6.0f;      ///< how hard the envelope pushes the sweep, dB (−24 … +36)
    float     attackMs      = 3.0f;      ///< the follower's attack (0.1 … 100 ms)
    float     releaseMs     = 120.0f;    ///< the follower's release — the "wow" (5 … 2000 ms)
    float     sweepLoHz     = 150.0f;    ///< one end of the sweep (20 … 12000 Hz; the ends may come in any order)
    float     sweepHiHz     = 3000.0f;   ///< the other end
    Direction direction     = Direction::up;
    Mode      mode          = Mode::lowPass;
    float     resonance     = 4.0f;      ///< the swept filter's Q (0.5 … 12)
    float     keepLows      = 1.0f;      ///< 0 = the whole signal is filtered, 1 = everything under keepLowsHz stays clean
    float     keepLowsHz    = 100.0f;    ///< where the clean lows end (30 … 300 Hz)
    float     mix           = 1.0f;      ///< 0 = dry (and nothing runs), 1 = all filter
};

namespace detail
{
    inline float clampf (float lo, float hi, float v) noexcept { return v < lo ? lo : (hi < v ? hi : v); }

    /** Equal within one relative float epsilon — the guard that makes a re-sent knob cost nothing. */
    inline bool nearlyEqual (float a, float b) noexcept
    {
        const float diff = std::fabs (a - b);
        return diff <= std::numeric_limits<float>::min()
            || diff <= std::numeric_limits<float>::epsilon() * std::max (std::fabs (a), std::fabs (b));
    }

    inline constexpr double kLog2E = 1.4426950408889634074;

    /** e^x, the same bits on every row. */
    inline double detExp (double x) noexcept { return core::det::exp2 (core::det::mul (x, kLog2E)); }

    /** The swept filter: eq::Svf's Cytomic TPT update, line for line, with BOTH taps out of one pair of
        integrators — low pass (v2) and band pass (v1, whose centre gain is Q). Coefficients from g and k only,
        so a per-sample redesign is three multiplies and a divide. State float, arithmetic double, as eq::Svf. */
    struct SweptSvf
    {
        static constexpr int kChannels = 2;

        double a1 = 1.0, a2 = 0.0, a3 = 0.0;
        float  ic1[kChannels] {}, ic2[kChannels] {};

        void design (double g, double k) noexcept
        {
            a1 = 1.0 / (1.0 + g * (g + k));
            a2 = g * a1;
            a3 = g * a2;
        }

        void tick (int ch, double u, double& lp, double& bp) noexcept
        {
            const double s1 = (double) ic1[ch], s2 = (double) ic2[ch];
            const double v3 = u - s2;
            const double v1 = a1 * s1 + a2 * v3;
            const double v2 = s2 + a2 * s1 + a3 * v3;
            ic1[ch] = (float) (2.0 * v1 - s1);
            ic2[ch] = (float) (2.0 * v2 - s2);
            lp = v2;
            bp = v1;
        }

        void clear (int ch) noexcept { ic1[ch] = 0.0f; ic2[ch] = 0.0f; }

        /** Law 8, on the owner's grid: the two words of a channel are one state — if either is not finite both
            go (eq::Svf's rule, for eq::Svf's reason), else each is flushed under 1e-15. */
        void flush (int nch) noexcept
        {
            for (int c = 0; c < nch; ++c)
            {
                if (! std::isfinite (ic1[c]) || ! std::isfinite (ic2[c])) { clear (c); continue; }
                core::flushDenormal (ic1[c]);
                core::flushDenormal (ic2[c]);
            }
        }

        /** The poison half alone, at the end of every call: cannot change a bit while the state is finite. */
        void heal (int nch) noexcept
        {
            for (int c = 0; c < nch; ++c)
                if (! std::isfinite (ic1[c]) || ! std::isfinite (ic2[c]))
                    clear (c);
        }
    };
} // namespace detail

class EnvelopeFilter
{
public:
    static constexpr int    kMaxChannels   = 2;
    static constexpr double kMinSampleRate = 8000.0, kMaxSampleRate = 768000.0;

    static constexpr float kMinSensDb  = -24.0f,  kMaxSensDb  = 36.0f;
    static constexpr float kMinAttack  = 0.1f,    kMaxAttack  = 100.0f;    // ms
    static constexpr float kMinRelease = 5.0f,    kMaxRelease = 2000.0f;   // ms
    static constexpr float kMinHz      = 20.0f,   kMaxHz      = 12000.0f;  // the sweep's ends
    static constexpr float kMinQ       = 0.5f,    kMaxQ       = 12.0f;
    static constexpr float kMinKeepHz  = 30.0f,   kMaxKeepHz  = 300.0f;

    /** How far below full opening the sweep is closed: the envelope's last 30 dB under (0 dBFS − sensitivity)
        walk the cutoff across the whole range, so a note decaying by 30 dB closes it all the way. */
    static constexpr double kSpanDb = 30.0;

    /** The cutoff's one smoother. Short: it must not blunt the attack (the follower's own is 3 ms), and long
        enough to take most of the rectifier's ripple off a low note — at a 41 Hz E the follower droops ~0.9 dB
        between the 82 Hz peaks of |x|, and this pole is down 9 dB there. */
    static constexpr double kSmoothSeconds = 0.004;
    static constexpr double kGlideSeconds  = 0.020;   ///< the knobs' glide (one-pole time constant)
    static constexpr double kModeSeconds   = 0.010;   ///< the low-pass ↔ band-pass crossfade (linear, counted)

    /** Never above this share of the rate, whatever the range says: the prewarp's tan() runs away at Nyquist. */
    static constexpr double kMaxCutoffShare = 0.45;

    EnvelopeFilter() noexcept
    {
        follower_.setDetector (dynamics::Detector::Peak);
        applyTimes();
        lfLo_ = core::det::log2 ((double) loHz_);
        lfHi_ = core::det::log2 ((double) hiHz_);
        lf_   = lfLo_;
    }

    [[nodiscard]] bool prepare (double newSampleRate, int maxBlock, int maxChannels) noexcept
    {
        prepared_ = false;   // law 11b: whatever happens next, the old preparation no longer stands

        if (! (newSampleRate >= kMinSampleRate && newSampleRate <= kMaxSampleRate) || maxBlock < 1
            || maxChannels < 1 || maxChannels > kMaxChannels)
            return false;

        fs_       = newSampleRate;
        maxBlock_ = maxBlock;
        maxCh_    = maxChannels;

        follower_.prepare (fs_);
        applyTimes();
        split_.prepare (fs_, maxCh_);
        split_.setParams (eq::FilterType::HighPass, (double) keepHz_, kSplitQ, 0.0);

        glideC_  = 1.0 - detail::detExp (-1.0 / (kGlideSeconds * fs_));
        smoothC_ = 1.0 - detail::detExp (-1.0 / (kSmoothSeconds * fs_));
        modeLen_ = std::max (1, (int) (kModeSeconds * fs_ + 0.5));
        fMax_    = kMaxCutoffShare * fs_;

        cleared_ = false;
        reset();
        prepared_ = true;
        return true;
    }

    bool isPrepared()  const noexcept { return prepared_; }
    int  maxChannels() const noexcept { return prepared_ ? maxCh_ : 0; }
    int  maxBlock()    const noexcept { return prepared_ ? maxBlock_ : 0; }
    static constexpr int latencySamples() noexcept { return 0; }

    /** A stream restart: every state cleared, the grid re-anchored, and the glides put on their targets by the
        next call. Cheap-idempotent: a stage out of the rig may call it every block. */
    void reset() noexcept
    {
        primed_ = false;
        if (cleared_)
            return;

        follower_.reset();
        split_.reset();
        for (int ch = 0; ch < kMaxChannels; ++ch)
        {
            svf_.clear (ch);
            live_[ch] = false;
        }
        grid_.reset();
        cleared_ = true;
    }

    void setParams (const EnvelopeFilterParams& p) noexcept
    {
        setSensitivityDb (p.sensitivityDb);
        setTimes (p.attackMs, p.releaseMs);
        setRange (p.sweepLoHz, p.sweepHiHz);
        setDirection (p.direction);
        setMode (p.mode);
        setResonance (p.resonance);
        setKeepLows (p.keepLows);
        setKeepLowsHz (p.keepLowsHz);
        setMix (p.mix);
    }

    void setSensitivityDb (float db) noexcept
    {
        if (std::isfinite (db))
            sensTgt_ = (double) detail::clampf (kMinSensDb, kMaxSensDb, db);
    }

    /** The follower's attack and release. Retunes at once (two exponentials, only when a value changes). */
    void setTimes (float attackMs, float releaseMs) noexcept
    {
        const float a = std::isfinite (attackMs)  ? detail::clampf (kMinAttack, kMaxAttack, attackMs)    : atkMs_;
        const float r = std::isfinite (releaseMs) ? detail::clampf (kMinRelease, kMaxRelease, releaseMs) : relMs_;
        if (detail::nearlyEqual (a, atkMs_) && detail::nearlyEqual (r, relMs_))
            return;
        atkMs_ = a;
        relMs_ = r;
        applyTimes();
    }

    /** The two ends of the sweep, in Hz, in either order. A move re-aims the cutoff, which then travels there
        through the smoother. */
    void setRange (float aHz, float bHz) noexcept
    {
        const float a = std::isfinite (aHz) ? detail::clampf (kMinHz, kMaxHz, aHz) : loHz_;
        const float b = std::isfinite (bHz) ? detail::clampf (kMinHz, kMaxHz, bHz) : hiHz_;
        const float lo = std::min (a, b), hi = std::max (a, b);
        if (detail::nearlyEqual (lo, loHz_) && detail::nearlyEqual (hi, hiHz_))
            return;
        loHz_ = lo;
        hiHz_ = hi;
        lfLo_ = core::det::log2 ((double) loHz_);
        lfHi_ = core::det::log2 ((double) hiHz_);
    }

    void setDirection (Direction d) noexcept { up_ = d != Direction::down; }
    void setMode (Mode m) noexcept { lowPass_ = m != Mode::bandPass; }

    void setResonance (float q) noexcept
    {
        if (std::isfinite (q))
            qTgt_ = (double) detail::clampf (kMinQ, kMaxQ, q);
    }

    void setKeepLows (float amount) noexcept
    {
        if (std::isfinite (amount))
            keepTgt_ = (double) detail::clampf (0.0f, 1.0f, amount);
    }

    /** Where the clean lows end. Retunes the split at once — a menu setting. */
    void setKeepLowsHz (float hz) noexcept
    {
        if (! std::isfinite (hz))
            return;
        const float want = detail::clampf (kMinKeepHz, kMaxKeepHz, hz);
        if (detail::nearlyEqual (want, keepHz_))
            return;
        keepHz_ = want;
        split_.setParams (eq::FilterType::HighPass, (double) keepHz_, kSplitQ, 0.0);
    }

    /** 0 = dry, and the stage sleeps once the glide lands there; 1 = all filter. */
    void setMix (float m) noexcept
    {
        if (std::isfinite (m))
            mixTgt_ = (double) detail::clampf (0.0f, 1.0f, m);
    }

    /** In place. Channels at index >= nch with n > 0 have stopped (law 11a); nch == 0 is a gap (law 11c). */
    [[nodiscard]] bool process (float* const* channels, int numChannels, int numSamples) noexcept
    {
        if (! prepared_ || numChannels < 0 || numSamples < 0 || numChannels > maxCh_)
            return false;   // law 11: refused before anything moves
        if (numSamples == 0)
            return true;    // the only true no-op: no time, no edge
        if (numChannels > 0 && channels == nullptr)
            return false;

        const int nch = numChannels;

        for (int ch = nch; ch < maxCh_; ++ch)   // law 11a — the channels this call does not carry have stopped
            stopChannel (ch);

        if (! primed_)
            prime();

        if (asleep_)
        {
            if (! (mixTgt_ > 0.0))
                return true;   // AMOUNT 0: the buffer untouched, nothing computed
            wake();
        }

        cleared_ = false;
        for (int ch = 0; ch < nch; ++ch)
            live_[ch] = true;

        int done = 0;
        while (done < numSamples)
        {
            const int seg = grid_.segment (numSamples - done);

            for (int i = done; i < done + seg; ++i)
            {
                // The detector hears the loudest finite sample of the call; a non-finite one is silence to it.
                float key = 0.0f;
                for (int ch = 0; ch < nch; ++ch)
                {
                    const float a = std::fabs (channels[ch][i]);
                    if (std::isfinite (a) && a > key)
                        key = a;
                }

                const bool moving = controlStep (key);

                if (! (mix_ > 0.0) && ! (mixTgt_ > 0.0))
                {
                    fallAsleep();   // the mix has landed on 0: from this sample on the call is dry
                    publish();
                    return true;
                }

                if (nch == 0 && ! moving && core::sameBits (follower_.stateWord(), 0.0f))
                {
                    // A gap whose shared state stands still: silence would not move it either, and every
                    // channel has stopped, so a flush would find nothing — the rest is audio time passing.
                    grid_.skip (numSamples - done);
                    publish();
                    return true;
                }

                for (int ch = 0; ch < nch; ++ch)
                    channels[ch][i] = channelStep (ch, channels[ch][i]);
            }

            if (grid_.advance (seg))
                flushState();   // law 8a: clocked by audio time, not by the call
            done += seg;
        }

        svf_.heal (maxCh_);    // the poison half keeps the call's clock (law 8a)
        split_.healPoison();
        publish();
        return true;
    }

    /** Where the cutoff stands, in Hz, at the end of the last call (at rest while the stage sleeps) — for the
        picture, from any thread. */
    float shownCutoffHz() const noexcept { return shownHz_.load (std::memory_order_relaxed); }

    /** True while amount 0 holds the stage out of the path (nothing is computed). */
    bool isAsleep() const noexcept { return primed_ && asleep_; }

    /** Where the knob glides stand — the tests read them to see a glide LAND on its target. */
    double currentMix() const noexcept { return mix_; }
    double currentResonance() const noexcept { return q_; }

    /** True when any recursive state sits in the subnormal band — the law-8 witness for the tests. */
    bool anySubnormalState() const noexcept
    {
        const auto sub = [] (auto v) { return std::fpclassify (v) == FP_SUBNORMAL; };
        for (int ch = 0; ch < kMaxChannels; ++ch)
            if (sub (svf_.ic1[ch]) || sub (svf_.ic2[ch]))
                return true;
        return sub (follower_.stateWord()) || sub (mix_) || sub (sens_) || sub (keep_) || sub (q_) || sub (lf_);
    }

private:
    static constexpr double kSplitQ    = 0.5;   // critically damped: the flattest complement (see the header)
    static constexpr double kPi        = 3.14159265358979323846;
    static constexpr double kDbPerLog2 = 6.0205999132796239042;    // 20·log10(2)
    static constexpr float  kEnvFloor  = 1.0e-9f;                  // −180 dBFS: below this the sweep is at rest
    static constexpr double kLandTol   = 1.0e-9;                   // a glide this close (× max(1,|target|)) lands

    /** One step of an exponential glide that LANDS: on its target once the update no longer moves it, or once
        it is within kLandTol of it (a glide toward zero would otherwise creep for ever). Returns whether the
        value changed. */
    static bool glide (double& cur, double target, double c) noexcept
    {
        const double before = cur;
        const double next = cur + (target - cur) * c;
        const bool   far  = std::fabs (target - next) > kLandTol * std::max (1.0, std::fabs (target));
        cur = (far && std::fabs (next - cur) > 0.0) ? next : target;
        return std::fabs (cur - before) > 0.0;
    }

    void applyTimes() noexcept { follower_.setTimes ((double) atkMs_, (double) relMs_); }

    /** The sweep's target, in log2(Hz), for an envelope. */
    double cutoffTarget (float env) const noexcept
    {
        double p = 0.0;
        if (env > kEnvFloor)
        {
            p = 1.0 + (kDbPerLog2 * core::det::log2 ((double) env) + sens_) / kSpanDb;
            p = p < 0.0 ? 0.0 : (p > 1.0 ? 1.0 : p);
        }
        const double span = lfHi_ - lfLo_;
        return up_ ? lfLo_ + p * span : lfHi_ - p * span;
    }

    void design() noexcept
    {
        double f = core::det::exp2 (lf_);
        if (f > fMax_)
            f = fMax_;
        cutoff_ = f;
        svf_.design (core::det::tan (kPi * f / fs_), k_);
    }

    void refreshQ() noexcept
    {
        k_    = 1.0 / q_;
        comp_ = q_ > 1.0 ? 1.0 / std::sqrt (q_) : 1.0;
    }

    /** One sample of the SHARED ballistics: the knob glides, the mode door, the follower and the cutoff. Every
        sample of a call runs it, at any width, and a gap runs it with the detector at digital zero (law 11c).
        Returns whether anything other than the follower moved. */
    bool controlStep (float key) noexcept
    {
        const bool mixMoved  = glide (mix_, mixTgt_, glideC_);
        const bool sensMoved = glide (sens_, sensTgt_, glideC_);
        const bool keepMoved = glide (keep_, keepTgt_, glideC_);
        const bool qMoved    = glide (q_, qTgt_, glideC_);
        if (qMoved)
            refreshQ();

        bool doorMoved = false;
        const int modeTo = lowPass_ ? modeLen_ : 0;
        if (modePos_ != modeTo)
        {
            modePos_ += modePos_ < modeTo ? 1 : -1;
            wLp_ = (double) modePos_ / (double) modeLen_;
            doorMoved = true;
        }

        const float env = follower_.process (key);
        const bool lfMoved = glide (lf_, cutoffTarget (env), smoothC_);
        if (lfMoved || qMoved)
            design();

        return mixMoved || sensMoved || keepMoved || qMoved || doorMoved || lfMoved;
    }

    /** One channel's sample: the split, the swept filter, the clean lows back, the mix. A non-finite input is
        not this stage's to change, and not for its state to hear: the filters run on a plain zero there (one
        path, so the state is exactly what a zero would have left) and the sample goes out as it came. */
    float channelStep (int ch, float xin) noexcept
    {
        const bool   finite = std::isfinite (xin);
        const float  xs   = finite ? xin : 0.0f;
        const float  h    = split_.processSample (ch, xs);
        const double x    = (double) xs;
        const double lows = x - (double) h;              // the complement: lows + h == x
        const double u    = x - keep_ * lows;            // what the filter hears: h + (1 - keepLows)*lows
        double lp = 0.0, bp = 0.0;
        svf_.tick (ch, u, lp, bp);
        const double y   = comp_ * (wLp_ * lp + (1.0 - wLp_) * bp);
        const double wet = keep_ * lows + y;
        return finite ? (float) (x + mix_ * (wet - x)) : xin;
    }

    void stopChannel (int ch) noexcept
    {
        if (ch < 0 || ch >= kMaxChannels || ! live_[ch])
            return;
        svf_.clear (ch);
        split_.resetChannel (ch);
        live_[ch] = false;
    }

    /** The glides on their targets, the follower empty, the cutoff at rest — a fresh stream. */
    void prime() noexcept
    {
        mix_  = mixTgt_;
        sens_ = sensTgt_;
        keep_ = keepTgt_;
        q_    = qTgt_;
        refreshQ();
        modePos_ = lowPass_ ? modeLen_ : 0;
        wLp_     = (double) modePos_ / (double) modeLen_;
        follower_.reset();
        lf_ = cutoffTarget (0.0f);
        design();
        grid_.reset();
        asleep_ = ! (mixTgt_ > 0.0);
        primed_ = true;
    }

    /** Out of sleep: a fresh stream, but the mix comes in from 0 — the slider's first move glides. */
    void wake() noexcept
    {
        for (int ch = 0; ch < kMaxChannels; ++ch)
            svf_.clear (ch);
        split_.reset();
        prime();
        mix_    = 0.0;
        asleep_ = false;
    }

    void fallAsleep() noexcept
    {
        follower_.reset();
        split_.reset();
        for (int ch = 0; ch < kMaxChannels; ++ch)
        {
            svf_.clear (ch);
            live_[ch] = false;
        }
        lf_ = cutoffTarget (0.0f);
        design();
        asleep_ = true;
    }

    void flushState() noexcept
    {
        follower_.flushDenormals();
        svf_.flush (maxCh_);
        split_.flushDenormals();
    }

    void publish() noexcept { shownHz_.store ((float) cutoff_, std::memory_order_relaxed); }

    std::atomic<float> shownHz_ { 150.0f };

    // The knobs' targets, as set (bounded).
    double mixTgt_ = 1.0, sensTgt_ = 6.0, keepTgt_ = 1.0, qTgt_ = 4.0;
    float  atkMs_ = 3.0f, relMs_ = 120.0f;
    float  loHz_ = 150.0f, hiHz_ = 3000.0f, keepHz_ = 100.0f;
    double lfLo_ = 0.0, lfHi_ = 0.0;   // log2 of the ends, set by the constructor and setRange
    bool   up_ = true, lowPass_ = true;

    // Where the glides stand (double, and they land).
    double mix_ = 1.0, sens_ = 6.0, keep_ = 1.0, q_ = 4.0;
    double k_ = 0.25, comp_ = 0.5;
    double lf_ = 0.0, cutoff_ = 150.0;
    int    modePos_ = 1, modeLen_ = 1;
    double wLp_ = 1.0;

    double fs_ = 48000.0, fMax_ = 21600.0;
    double glideC_ = 1.0, smoothC_ = 1.0;

    dynamics::EnvelopeFollower follower_;
    eq::DeterministicSvf       split_;
    detail::SweptSvf           svf_;
    core::StateGrid            grid_;

    bool live_[kMaxChannels] { false, false };   // carried by a call since the last clear (law 11a)
    bool primed_ = false, asleep_ = false, cleared_ = true;
    bool prepared_ = false;
    int  maxBlock_ = 0, maxCh_ = 0;
};

} // namespace felitronics::envfilter

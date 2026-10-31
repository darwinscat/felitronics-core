// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/core/DetMath.h>
#include <felitronics/core/FlushToZero.h>
#include <felitronics/core/Math.h>
#include <felitronics/core/StateGrid.h>
#include <felitronics/dynamics/EnvelopeFollower.h>
#include <felitronics/eq/Crossover2.h>
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
//   · DETECTOR — the loudest finite |sample| of the KEY into a peak dynamics::EnvelopeFollower (attack,
//     release). The key is the programme itself, or an external one (process's last two arguments): in an amp
//     the filter may sit after a compressor and an octaver, and a key taken before them lets it follow the
//     player's hand rather than a squashed envelope and a square sub. Linked: one envelope moves every channel.
//   · THE SWEEP — the envelope in dB, pushed by SENSITIVITY, places the cutoff on a log-frequency line between
//     the range's ends: fully open when envelope + sensitivity reaches 0 dBFS, closed kSpanDb below. UP rests at
//     the low end and a hit opens it; DOWN rests at the high end and a hit pushes it down. The PLAYING reaches
//     the cutoff through one kSmoothSeconds one-pole on log2(Hz); the hand's edits — the direction (a 0..1
//     blend of the two rests), the range's ends, the sensitivity — glide over kGlideSeconds before it, so a
//     direction flip travels as a sweep and is not heard as a hit. The filter is redesigned EVERY sample.
//   · THE FILTER — a Cytomic TPT state-variable filter, eq::Svf's update line for line, in a kernel of its own
//     (detail::SweptSvf): the per-sample redesign needs only g = tan(π·fc/fs) and k = 1/Q, where
//     eq::Svf::setParams designs a gain term too. Both taps are scaled by 1/√max(1, Q): at the cutoff the two
//     are then equal at every Q (Q/√Q), and the peak grows as √Q. The modes SOUND alike in level from about
//     Q 2 up, where the peak dominates; at low Q the low pass also carries its whole passband and is louder. The
//     PRICE of the scaling is the low pass's passband: −6.0 dB at Q 4, −10.8 dB at Q 12 (pinned by the tests).
//   · THE CLEAN LOWS — `keepLows` (0..1) takes the band under keepLowsHz out of the filter and adds it back
//     clean, so the note keeps its weight. HOW depends on the mode, because a protected band only survives if
//     the filtered part cannot arrive opposite in phase to it:
//       LOW PASS — the lows are x − HP(x) (a critically damped high pass, eq::DeterministicSvf), an exact
//         complement of the input, and the dry is x itself. The low pass's passband is in phase with x, so the
//         two can only add: measured over cutoff × Q × corner × mix, the band under the corner never drops more
//         than 0.07 dB. (The complement lifts it a little: at most +1.25 dB, at 0.71 × the corner.)
//       BAND PASS — that complement does NOT survive a band pass: its ±90° skirts rotate the filtered half
//         against the clean one, and a band pass parked just above the corner took a 98 Hz note down 13.6 dB.
//         So this mode splits with eq::DeterministicCrossover2, an LR4 pair whose two halves are IN PHASE at
//         every frequency: the band pass's output never has a negative real part, so it can only add to the low
//         half, and the band under the corner keeps at least the LR4 low share (−0.5 dB at half the corner,
//         −1.9 at 0.7, −6 at the corner). Its dry is the same LR4 allpass (low + high), so a mix never notches
//         the crossover the way x against an allpass would.
//     The two modes run side by side and the mode door crossfades their outputs (and their dries), so a mode
//     change is a crossfade of two live filters, never a restart. The corner glides (kGlideSeconds, on
//     log2 Hz) with both splits redesigned while it moves.
//   · THE MIX — out = dry + mix·(wet − dry), then through the ENGAGE DOOR: x + door·(that − x).
//
// AMOUNT 0 IS NOT THERE. When `mix` is set to 0 the engage door closes — a COUNTED kDoorSeconds ramp that lands
// exactly — and from the sample it lands the stage is IDLE: process() returns with the buffer untouched (bit-exact
// dry, NaN payloads included) and computes nothing. A wake is a RESTART: every state from silence, the grid
// re-anchored, the glides on their targets, the door opening from dry over kDoorSeconds and the mix gliding in
// from 0 (the filters starting from silence on a sounding note are faded in, not heard). The one-slider law
// that maps a face's slider onto these knobs is EnvelopeFilterLaw.h.
//
// CONTRACT (DSP-ARCHITECTURE.md §2).
//   · prepare(sampleRate, maxBlock, maxChannels): maxChannels 1 or 2, the rate 8 kHz … 768 kHz, maxBlock >= 1
//     (recorded, not a limit: the stage has no block-sized scratch). Anything else is REFUSED, and a refused
//     prepare DISARMS the instance (law 11b). prepare() performs a reset().
//   · process() is [[nodiscard]] bool (law 11): unprepared, a width above the prepared one, a negative extent,
//     or a null buffer with nch > 0 is refused before anything moves; any n >= 0 is processed in full.
//   · THE KEY (the five-argument form, the dynamics house shape): `key[0 .. numKeyChannels-1]`, each `n` long,
//     the SAME instant as `io`. `key == nullptr` or a count <= 0 is the self-keyed path, identical to the
//     three-argument form; `key` may be `io` itself (identical too). Any count; the loudest finite sample of the
//     key's channels is the detector's input.
//   · A NARROWER call is legal; a channel at an index >= nch with n > 0 has STOPPED (law 11a): its filters and
//     splits are cleared once, on that edge. nch == 0 with n > 0 is a gap: the SHARED ballistics — the
//     follower, the cutoff smoother, the glides and both doors — run through it exactly as through n samples of
//     digital silence at a live width (law 11c). AN EXTERNAL KEY IS STILL CONSUMED: the programme stopped, not
//     the key. A keyless gap stops computing once the shared state stands still.
//   · reset() is a stream RESTART: every state cleared, the audio-time grid re-anchored, and the next call
//     starts the glides ON their targets and the engage door open — a fresh stage does not fade in.
//   · GLIDES. mix, sensitivity, resonance, keepLows, the direction, the range's ends and keepLowsHz glide over
//     kGlideSeconds (one-pole, in double, and they LAND — LAW8-AUDIT.md item 6); the mode door and the engage
//     door are counted and land exactly. Only the follower's attack/release retune at once.
//   · Law 8: the follower, the swept filters' integrators and both splits are flushed at the house 1e-15 on the
//     core::StateGrid (64 samples of AUDIO time, law 8a), so the output does not depend on how the caller
//     slices the stream. The glides land rather than creep, so none of them reaches the subnormal band.
//   · POISON: a non-finite INPUT sample is not fed to the detector, the splits or the filters (they hear silence
//     there); the output sample keeps it — it is the caller's, not this stage's to change. A non-finite filter
//     state at the end of a call (a finite input can still overflow it) restarts that lane from silence, so the
//     next call is clean even when the overflow came after the call's last grid boundary. Knob values are
//     bounded before any arithmetic: a non-finite one is ignored, a finite one is clamped.
//   · THREADS: prepare, reset, the setters and process belong to ONE thread. shownCutoffHz is a relaxed atomic
//     that any other thread may read; every other getter is the audio thread's (the tests read them).
//   · Math: the transcendentals of the per-sample sweep (log2, exp2, tan) are core::det, the same bits on every
//     row — deliberately per sample, as delay::TapeDelay's tanh is (DetMath.h's cost note). The arithmetic
//     between them follows law 10 (a contracted a·b + c may move a last bit), and the follower's coefficients
//     are dynamics::EnvelopeFollower's own (system exp), so the cutoff is not claimed bit-equal across rows.
//==============================================================================

enum class Mode      { lowPass, bandPass };   ///< which filter is heard
enum class Direction { up, down };            ///< up: rests low, a hit opens it; down: rests high, a hit closes it

/** Every knob, in the units the stage takes. `EnvelopeFilter::setParams` applies them through the same guarded
    setters a caller may use one by one. `sensitivityDb` and `mix` are the ONE-SLIDER LAW's (EnvelopeFilterLaw.h):
    a face that drives the stage through the law must not expose them. */
struct EnvelopeFilterParams
{
    float     sensitivityDb = 6.0f;      ///< LAW-OWNED. How hard the envelope pushes the sweep, dB (−24 … +36)
    float     attackMs      = 3.0f;      ///< the follower's attack (0.1 … 100 ms)
    float     releaseMs     = 120.0f;    ///< the follower's release — the "wow" (5 … 2000 ms)
    float     sweepLoHz     = 150.0f;    ///< one end of the sweep (20 … 12000 Hz; the ends may come in any order)
    float     sweepHiHz     = 3000.0f;   ///< the other end
    Direction direction     = Direction::up;
    Mode      mode          = Mode::lowPass;
    float     resonance     = 4.0f;      ///< the swept filter's Q (0.5 … 12)
    float     keepLows      = 1.0f;      ///< 0 = the whole signal is filtered, 1 = everything under keepLowsHz stays clean
    float     keepLowsHz    = 150.0f;    ///< the clean lows' corner (30 … 300 Hz) — a fifth above G2
    float     mix           = 1.0f;      ///< LAW-OWNED. 0 = dry (the stage goes idle), 1 = all filter
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

    /** The swept filter: eq::Svf's Cytomic TPT update, line for line, giving both taps — low pass (v2) and band
        pass (v1, whose centre gain is Q) — from one pair of integrators per LANE. One coefficient set serves every
        lane: the stage runs the low-pass path in lanes 0..1 and the band-pass path in lanes 2..3, both at the same
        cutoff. Coefficients from g and k only, so a per-sample redesign is three multiplies and a divide. State
        float, arithmetic double, as eq::Svf. */
    struct SweptSvf
    {
        static constexpr int kLanes = 4;

        double a1 = 1.0, a2 = 0.0, a3 = 0.0;
        float  ic1[kLanes] {}, ic2[kLanes] {};

        void design (double g, double k) noexcept
        {
            a1 = 1.0 / (1.0 + g * (g + k));
            a2 = g * a1;
            a3 = g * a2;
        }

        void tick (int lane, double u, double& lp, double& bp) noexcept
        {
            const double s1 = (double) ic1[lane], s2 = (double) ic2[lane];
            const double v3 = u - s2;
            const double v1 = a1 * s1 + a2 * v3;
            const double v2 = s2 + a2 * s1 + a3 * v3;
            ic1[lane] = (float) (2.0 * v1 - s1);
            ic2[lane] = (float) (2.0 * v2 - s2);
            lp = v2;
            bp = v1;
        }

        void clear (int lane) noexcept { ic1[lane] = 0.0f; ic2[lane] = 0.0f; }

        bool finite (int lane) const noexcept { return std::isfinite (ic1[lane]) && std::isfinite (ic2[lane]); }

        /** Law 8, on the owner's grid: the two words of a lane are one state — if either is not finite both go
            (eq::Svf's rule, for eq::Svf's reason), else each is flushed under 1e-15. */
        void flush() noexcept
        {
            for (int c = 0; c < kLanes; ++c)
            {
                if (! finite (c)) { clear (c); continue; }
                core::flushDenormal (ic1[c]);
                core::flushDenormal (ic2[c]);
            }
        }

        /** The poison half alone, at the end of every call: cannot change a bit while the state is finite. */
        void heal() noexcept
        {
            for (int c = 0; c < kLanes; ++c)
                if (! finite (c))
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

    /** The PLAYING's smoother on the cutoff. Short: it must not blunt the attack (the follower's own is 3 ms),
        and long enough to take some of the rectifier's ripple off a low note — at a 41 Hz E the follower droops
        ~0.9 dB between the 82 Hz peaks of |x|, and this pole is down 7.2 dB there. */
    static constexpr double kSmoothSeconds = 0.004;
    static constexpr double kGlideSeconds  = 0.020;   ///< the hand's glides (one-pole time constant)
    static constexpr double kModeSeconds   = 0.010;   ///< the low-pass ↔ band-pass crossfade (linear, counted)
    static constexpr double kDoorSeconds   = 0.010;   ///< the engage door, in and out (linear, counted)

    /** Never above this share of the rate, whatever the range says: the prewarp's tan() runs away at Nyquist. */
    static constexpr double kMaxCutoffShare = 0.45;

    EnvelopeFilter() noexcept
    {
        follower_.setDetector (dynamics::Detector::Peak);
        applyTimes();
        lfLoTgt_ = lfLo_ = core::det::log2 ((double) loHz_);
        lfHiTgt_ = lfHi_ = core::det::log2 ((double) hiHz_);
        kfTgt_   = kf_   = core::det::log2 ((double) keepHz_);
        lf_ = lfLo_;
        refreshFloor();
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
        splitLp_.prepare (fs_, maxCh_);
        splitBp_.prepare (fs_, maxCh_);

        glideC_  = 1.0 - detail::detExp (-1.0 / (kGlideSeconds * fs_));
        smoothC_ = 1.0 - detail::detExp (-1.0 / (kSmoothSeconds * fs_));
        modeLen_ = std::max (1, (int) (kModeSeconds * fs_ + 0.5));
        doorLen_ = std::max (1, (int) (kDoorSeconds * fs_ + 0.5));
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
        clearAll();
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

    /** The two ends of the sweep, in Hz, in either order. Each end glides to its new place. */
    void setRange (float aHz, float bHz) noexcept
    {
        const float a = std::isfinite (aHz) ? detail::clampf (kMinHz, kMaxHz, aHz) : loHz_;
        const float b = std::isfinite (bHz) ? detail::clampf (kMinHz, kMaxHz, bHz) : hiHz_;
        const float lo = std::min (a, b), hi = std::max (a, b);
        if (detail::nearlyEqual (lo, loHz_) && detail::nearlyEqual (hi, hiHz_))
            return;
        loHz_ = lo;
        hiHz_ = hi;
        lfLoTgt_ = core::det::log2 ((double) loHz_);
        lfHiTgt_ = core::det::log2 ((double) hiHz_);
    }

    /** Glides: the rest moves from one end to the other as a sweep, not a jump. */
    void setDirection (Direction d) noexcept { dirTgt_ = d == Direction::down ? 1.0 : 0.0; }
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

    /** The clean lows' corner. Glides (on log2 Hz), with both splits redesigned while it moves. */
    void setKeepLowsHz (float hz) noexcept
    {
        if (! std::isfinite (hz))
            return;
        const float want = detail::clampf (kMinKeepHz, kMaxKeepHz, hz);
        if (detail::nearlyEqual (want, keepHz_))
            return;
        keepHz_ = want;
        kfTgt_  = core::det::log2 ((double) keepHz_);
    }

    /** 0 = dry: the engage door closes and the stage goes idle; 1 = all filter. */
    void setMix (float m) noexcept
    {
        if (std::isfinite (m))
            mixTgt_ = (double) detail::clampf (0.0f, 1.0f, m);
    }

    /** In place, self-keyed. Channels at index >= nch with n > 0 have stopped (law 11a); nch == 0 is a gap. */
    [[nodiscard]] bool process (float* const* channels, int numChannels, int numSamples) noexcept
    {
        return process (channels, numChannels, numSamples, nullptr, 0);
    }

    /** In place, keyed: the detector reads `key[0 .. numKeyChannels-1]` (each numSamples long, the same instant
        as `channels`); `key == nullptr` or `numKeyChannels <= 0` is the self-keyed form. */
    [[nodiscard]] bool process (float* const* channels, int numChannels, int numSamples,
                                const float* const* key, int numKeyChannels) noexcept
    {
        if (! prepared_ || numChannels < 0 || numSamples < 0 || numChannels > maxCh_)
            return false;   // law 11: refused before anything moves
        if (numSamples == 0)
            return true;    // the only true no-op: no time, no edge
        if (numChannels > 0 && channels == nullptr)
            return false;

        const int  nch   = numChannels;
        const bool keyed = key != nullptr && numKeyChannels > 0;

        for (int ch = nch; ch < maxCh_; ++ch)   // law 11a — the channels this call does not carry have stopped
            stopChannel (ch);

        if (! primed_)
            prime();

        if (idle_)
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
                // The detector hears the loudest finite sample of its key; a non-finite one is silence to it.
                float k = 0.0f;
                if (keyed)
                {
                    for (int c = 0; c < numKeyChannels; ++c)
                    {
                        const float a = std::fabs (key[c][i]);
                        if (std::isfinite (a) && a > k)
                            k = a;
                    }
                }
                else
                {
                    for (int ch = 0; ch < nch; ++ch)
                    {
                        const float a = std::fabs (channels[ch][i]);
                        if (std::isfinite (a) && a > k)
                            k = a;
                    }
                }

                const bool moving = controlStep (k);

                if (doorPos_ == 0 && ! (mixTgt_ > 0.0))
                {
                    goIdle();   // the engage door has landed shut: from this sample on the call is dry
                    publish();
                    return true;
                }

                if (nch == 0 && ! keyed && ! moving && core::sameBits (follower_.stateWord(), 0.0f))
                {
                    // A keyless gap whose shared state stands still: silence would not move it either, and
                    // every channel has stopped, so a flush would find nothing — the rest is audio time passing.
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

        svf_.heal();    // the poison half keeps the call's clock (law 8a)
        splitLp_.healPoison();
        splitBp_.healPoison();
        publish();
        return true;
    }

    /** Where the cutoff stands, in Hz, at the end of the last call (at rest while idle) — for the picture, from
        any thread. */
    float shownCutoffHz() const noexcept { return shownHz_.load (std::memory_order_relaxed); }

    /** True while amount 0 holds the stage out of the path (nothing is computed). */
    bool isIdle() const noexcept { return primed_ && idle_; }

    /** Where the glides stand — the tests read them to see a glide LAND on its target. Audio thread only. */
    double currentMix() const noexcept { return mix_; }
    double currentResonance() const noexcept { return q_; }

    /** True when every swept-filter lane holds a finite state — the poison witness for the tests. */
    bool filterStateIsFinite() const noexcept
    {
        for (int c = 0; c < detail::SweptSvf::kLanes; ++c)
            if (! svf_.finite (c))
                return false;
        return true;
    }

    /** True when any recursive state sits in the subnormal band — the law-8 witness for the tests. */
    bool anySubnormalState() const noexcept
    {
        const auto sub = [] (auto v) { return std::fpclassify (v) == FP_SUBNORMAL; };
        for (int c = 0; c < detail::SweptSvf::kLanes; ++c)
            if (sub (svf_.ic1[c]) || sub (svf_.ic2[c]))
                return true;
        return sub (follower_.stateWord()) || sub (mix_) || sub (sens_) || sub (keep_) || sub (q_) || sub (lf_)
            || sub (dir_) || sub (lfLo_) || sub (lfHi_) || sub (kf_);
    }

private:
    static constexpr double kSplitQ    = 0.5;                      // critically damped: the flattest complement
    static constexpr double kPi        = 3.14159265358979323846;
    static constexpr double kDbPerLog2 = 6.0205999132796239042;    // 20·log10(2)
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

    /** One step of a counted door toward its end; returns whether it moved. */
    static bool step (int& pos, int to) noexcept
    {
        if (pos == to)
            return false;
        pos += pos < to ? 1 : -1;
        return true;
    }

    void applyTimes() noexcept { follower_.setTimes ((double) atkMs_, (double) relMs_); }

    /** The envelope under which the sweep is at rest whatever log2 would say — so a decaying tail, or a long
        gap, costs no logarithm. p <= 0 exactly when envDb + sensitivity <= −kSpanDb. */
    void refreshFloor() noexcept
    {
        envFloor_ = (float) core::det::pow10 (-(kSpanDb + sens_) / 20.0);
    }

    /** The sweep's target, in log2(Hz), for an envelope. */
    double cutoffTarget (float env) const noexcept
    {
        double p = 0.0;
        if (env > envFloor_)
        {
            p = 1.0 + (kDbPerLog2 * core::det::log2 ((double) env) + sens_) / kSpanDb;
            p = p < 0.0 ? 0.0 : (p > 1.0 ? 1.0 : p);
        }
        const double span = lfHi_ - lfLo_;
        const double upAt = lfLo_ + p * span, downAt = lfHi_ - p * span;
        return upAt + dir_ * (downAt - upAt);
    }

    void design() noexcept
    {
        double f = core::det::exp2 (lf_);
        if (f > fMax_)
            f = fMax_;
        cutoff_ = f;
        svf_.design (core::det::tan (kPi * f / fs_), k_);
    }

    void designSplits() noexcept
    {
        const double hz = core::det::exp2 (kf_);
        splitLp_.setParams (eq::FilterType::HighPass, hz, kSplitQ, 0.0);
        splitBp_.setFrequency ((float) hz);
    }

    void refreshQ() noexcept
    {
        k_    = 1.0 / q_;
        comp_ = q_ > 1.0 ? 1.0 / std::sqrt (q_) : 1.0;
    }

    /** One sample of the SHARED ballistics: the glides, both doors, the follower and the cutoff. Every sample of
        a call runs it, at any width, and a gap runs it with the detector at its key — digital zero when keyless
        (law 11c). Returns whether anything other than the follower moved. */
    bool controlStep (float k) noexcept
    {
        const bool mixMoved  = glide (mix_, mixTgt_, glideC_);
        const bool sensMoved = glide (sens_, sensTgt_, glideC_);
        const bool keepMoved = glide (keep_, keepTgt_, glideC_);
        const bool qMoved    = glide (q_, qTgt_, glideC_);
        const bool dirMoved  = glide (dir_, dirTgt_, glideC_);
        const bool loMoved   = glide (lfLo_, lfLoTgt_, glideC_);
        const bool hiMoved   = glide (lfHi_, lfHiTgt_, glideC_);
        const bool kfMoved   = glide (kf_, kfTgt_, glideC_);
        if (qMoved)
            refreshQ();
        if (sensMoved)
            refreshFloor();
        if (kfMoved)
            designSplits();

        const bool modeMoved = step (modePos_, lowPass_ ? modeLen_ : 0);
        if (modeMoved)
            wLp_ = (double) modePos_ / (double) modeLen_;
        const bool doorMoved = step (doorPos_, mixTgt_ > 0.0 ? doorLen_ : 0);
        if (doorMoved)
            door_ = (double) doorPos_ / (double) doorLen_;

        const float env = follower_.process (k);
        const bool lfMoved = glide (lf_, cutoffTarget (env), smoothC_);
        if (lfMoved || qMoved)
            design();

        return mixMoved || sensMoved || keepMoved || qMoved || dirMoved || loMoved || hiMoved || kfMoved
            || modeMoved || doorMoved || lfMoved;
    }

    /** One channel's sample: both modes' splits and filters, the clean lows back, the mode door, the mix, the
        engage door. A non-finite input is not this stage's to change, and not for its state to hear: everything
        runs on a plain zero there (one path, so the state is exactly what a zero would have left) and the
        sample goes out as it came. */
    float channelStep (int ch, float xin) noexcept
    {
        const bool   finite = std::isfinite (xin);
        const float  xs = finite ? xin : 0.0f;
        const double x  = (double) xs;

        // LOW PASS: the complement split, the dry is x.
        const float  hc   = splitLp_.processSample (ch, xs);
        const double lowC = x - (double) hc;                    // lowC + hc == x
        double lp = 0.0, unusedBp = 0.0;
        svf_.tick (ch, x - keep_ * lowC, lp, unusedBp);
        const double wetLp = keep_ * lowC + comp_ * lp;

        // BAND PASS: the LR4 split, in phase at every frequency; the dry is its allpass.
        float l4 = 0.0f, h4 = 0.0f;
        splitBp_.processSample (ch, xs, l4, h4);
        double unusedLp = 0.0, bp = 0.0;
        svf_.tick (2 + ch, (double) h4 + (1.0 - keep_) * (double) l4, unusedLp, bp);
        const double wetBp = keep_ * (double) l4 + comp_ * bp;
        const double dryBp = (double) l4 + (double) h4;

        const double wet = wetBp + wLp_ * (wetLp - wetBp);
        const double dry = dryBp + wLp_ * (x - dryBp);
        const double mixed = dry + mix_ * (wet - dry);
        return finite ? (float) (x + door_ * (mixed - x)) : xin;
    }

    void stopChannel (int ch) noexcept
    {
        if (ch < 0 || ch >= kMaxChannels || ! live_[ch])
            return;
        svf_.clear (ch);
        svf_.clear (2 + ch);
        splitLp_.resetChannel (ch);
        splitBp_.resetChannel (ch);
        live_[ch] = false;
    }

    void clearAll() noexcept
    {
        follower_.reset();
        splitLp_.reset();
        splitBp_.reset();
        for (int c = 0; c < detail::SweptSvf::kLanes; ++c)
            svf_.clear (c);
        for (int ch = 0; ch < kMaxChannels; ++ch)
            live_[ch] = false;
    }

    /** The glides on their targets, the follower empty, the cutoff at rest, the engage door open — a fresh
        stream. */
    void prime() noexcept
    {
        mix_  = mixTgt_;
        sens_ = sensTgt_;
        keep_ = keepTgt_;
        q_    = qTgt_;
        dir_  = dirTgt_;
        lfLo_ = lfLoTgt_;
        lfHi_ = lfHiTgt_;
        kf_   = kfTgt_;
        refreshQ();
        refreshFloor();
        designSplits();
        modePos_ = lowPass_ ? modeLen_ : 0;
        wLp_     = (double) modePos_ / (double) modeLen_;
        idle_    = ! (mixTgt_ > 0.0);
        doorPos_ = idle_ ? 0 : doorLen_;
        door_    = (double) doorPos_ / (double) doorLen_;
        follower_.reset();
        lf_ = cutoffTarget (0.0f);
        design();
        grid_.reset();
        primed_ = true;
    }

    /** Out of idle: a RESTART — every state from silence, the grid re-anchored, the glides on their targets —
        with the engage door opening from dry over kDoorSeconds and the mix gliding in from 0, so the filters'
        own start from silence is faded in rather than heard. */
    void wake() noexcept
    {
        clearAll();
        prime();
        doorPos_ = 0;
        door_    = 0.0;
        mix_     = 0.0;
        idle_    = false;
    }

    void goIdle() noexcept
    {
        clearAll();
        lf_ = cutoffTarget (0.0f);
        design();
        idle_ = true;
    }

    void flushState() noexcept
    {
        follower_.flushDenormals();
        svf_.flush();
        splitLp_.flushDenormals();
        splitBp_.flushDenormals();
    }

    void publish() noexcept { shownHz_.store ((float) cutoff_, std::memory_order_relaxed); }

    std::atomic<float> shownHz_ { 150.0f };

    // The knobs' targets, as set (bounded).
    double mixTgt_ = 1.0, sensTgt_ = 6.0, keepTgt_ = 1.0, qTgt_ = 4.0, dirTgt_ = 0.0;
    double lfLoTgt_ = 0.0, lfHiTgt_ = 0.0, kfTgt_ = 0.0;      // log2 Hz, set by the constructor and the setters
    float  atkMs_ = 3.0f, relMs_ = 120.0f;
    float  loHz_ = 150.0f, hiHz_ = 3000.0f, keepHz_ = 150.0f;
    bool   lowPass_ = true;

    // Where the glides stand (double, and they land).
    double mix_ = 1.0, sens_ = 6.0, keep_ = 1.0, q_ = 4.0, dir_ = 0.0;
    double lfLo_ = 0.0, lfHi_ = 0.0, kf_ = 0.0;
    double k_ = 0.25, comp_ = 0.5;
    double lf_ = 0.0, cutoff_ = 150.0;
    float  envFloor_ = 0.0f;
    int    modePos_ = 1, modeLen_ = 1, doorPos_ = 1, doorLen_ = 1;
    double wLp_ = 1.0, door_ = 1.0;

    double fs_ = 48000.0, fMax_ = 21600.0;
    double glideC_ = 1.0, smoothC_ = 1.0;

    dynamics::EnvelopeFollower      follower_;
    eq::DeterministicSvf            splitLp_;   // the low-pass mode's complement split
    eq::DeterministicCrossover2     splitBp_;   // the band-pass mode's LR4 split
    detail::SweptSvf                svf_;
    core::StateGrid                 grid_;

    bool live_[kMaxChannels] { false, false };   // carried by a call since the last clear (law 11a)
    bool primed_ = false, idle_ = false, cleared_ = true;
    bool prepared_ = false;
    int  maxBlock_ = 0, maxCh_ = 0;
};

} // namespace felitronics::envfilter

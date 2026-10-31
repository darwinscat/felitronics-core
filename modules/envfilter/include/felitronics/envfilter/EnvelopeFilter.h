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
//     The two modes run side by side, so a mode change is a crossfade of two live filters, never a restart. The
//     corner glides (kGlideSeconds, on log2 Hz) with both splits redesigned while it moves.
//   · THE MIX — each mode's out = its dry + mix·(its wet − its dry).
//
// THE DOORS — HOW THE STAGE GETS FROM ONE PHASE TO THE OTHER WITHOUT A HOLE. Bypass and the low pass are
// referenced to x; the band pass is referenced to the LR4 allpass, which is −1 at the corner. A crossfade between
// the two families has a hole there — at mix 0.1 a tone at the corner fell 7.5 dB through a 10 ms door. So no door
// ever crosses families. The x family has one door: the LOW-PASS DOOR (x ↔ the low pass, kDoorSeconds). The band
// pass is reached by a LADDER of three counted stages, each phase-continuous: (1) a kRungSeconds crossfade from
// the x-family signal to the same signal through the LR4 pair's own allpass parked at the CEILING (kLadderOctaves
// above the corner, where it is ≈ +1 over the protected bass); (2) that allpass's corner sweeps down to the
// clean-lows corner over kCornerSweepSeconds (an exponential approach on log2 Hz) — a pure allpass, so the level
// holds while the phase turns, which is heard as a small pitch bend (≤ ~48 cents, see kCornerSweepSeconds); (3) a
// kRungSeconds crossfade from it to the band pass, which now shares its phase. Down is the same in reverse. The
// low-pass door stays open until the band pass is fully in, and the ladder comes down only once it is open again,
// so a mode change never passes through dry. Measured over tones at and above the corner, cutoff at it and 10×
// above, Q 4 and 12, mix 0.1 and 1, both modes: the weakest whole cycle through a wake, a disengage and either
// mode change sits within 0.92 dB of the lower settled level. What remains is the two taps' own difference: at
// the cutoff they are 90° apart, so the crossing (a rung) dips a tone there by ~3 dB (pinned).
//
// AMOUNT 0 IS NOT THERE. When `mix` is set to 0 the doors close — the low-pass door, or the ladder down
// (2·kRungSeconds + kCornerSweepSeconds = 250 ms), all counted, landing exactly — and from the sample the last one lands
// the stage is IDLE: process() returns with the buffer untouched (bit-exact dry, NaN payloads included) and
// computes nothing. A wake is a RESTART: every state from silence, the grid re-anchored, the glides on their
// targets, the doors opening from dry and the mix gliding in from 0 (the filters starting from silence on a
// sounding note are faded in, not heard). The one-slider law that maps a face's slider onto these knobs is
// EnvelopeFilterLaw.h.
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
//     starts the glides ON their targets and the doors where the mode puts them — a fresh stage does not fade in.
//   · GLIDES. mix, sensitivity, resonance, keepLows, the direction, the range's ends and keepLowsHz glide over
//     kGlideSeconds (one-pole, in double, and they LAND — LAW8-AUDIT.md item 6); the doors are counted and land
//     exactly. Only the follower's attack/release retune at once.
//   · Law 8: the follower, the swept filters' integrators, both splits and the ladder's allpass are flushed at the house 1e-15 on the
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
    /** One step of a counted door toward `to`; returns whether it moved. */
    inline bool stepToward (int& pos, int to) noexcept
    {
        if (pos == to)
            return false;
        pos += pos < to ? 1 : -1;
        return true;
    }

    /** THE DOOR MACHINE, one sample (see THE DOORS in the stage's header). `lp` is the low-pass door (0 … door),
        `lad` the band pass's ladder (0 … 2·rung + sweep: rung 1, the sweep, rung 2). The ladder climbs whenever the
        band pass is wanted. It comes down once the low-pass door stands where it is going — or at once from
        STRICTLY INSIDE rung 2, where the low-pass door is frozen and waiting for it would wait for ever (a mode
        change to low pass, or the slider to 0, in the band pass's last rung froze both doors in a partial blend
        that never idled); rung 2 crosses within one phase family whatever the low-pass door holds. At the TOP it
        still waits: the low-pass door can move there, and opening it first is what keeps band pass -> low pass
        from passing through dry. The low-pass door stays open until the band pass is fully in, and moves only
        while the band pass is out of the picture (below rung 2) or wholly in it. Pure, so the tests sweep EVERY
        state against every target: none is stuck. */
    inline void stepDoors (int& lp, int& lad, int door, int rung, int sweep, bool engaged, bool lowPass,
                           bool& lpMoved, bool& ladMoved) noexcept
    {
        const int top   = 2 * rung + sweep;
        const int lpTo  = engaged && (lowPass || (lp == door && lad < top)) ? door : 0;
        const int ladTo = engaged && ! lowPass ? top : 0;
        ladMoved = false;
        lpMoved  = false;
        if (ladTo > lad || lp == lpTo || (lad > rung + sweep && lad < top))
            ladMoved = stepToward (lad, ladTo);
        if (lad <= rung + sweep || lad == top)
            lpMoved = stepToward (lp, lpTo);
    }
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
    static constexpr double kDoorSeconds   = 0.010;   ///< the low-pass door's travel (linear, counted, lands exactly)
    /** The band pass's LADDER (see THE DOORS above). Its allpass turns a bass tone's phase by up to ~320° (−96° at
        82 Hz, −180° at the 150 Hz corner, −319° at 600 Hz), and turning it is a PITCH BEND: Δf = turn rate / 2π.
        The first ladder (0.45 fs ceiling, a 50 ms smoothstep sweep, 10 ms rungs) bent 82/150/300 Hz tones by
        309/339/292 cents. What sets the bend is the time the turn is spread over and how evenly: the ceiling sits
        kLadderOctaves above the corner (so the sweep spends its time where bass phases turn, not in the treble),
        the corner APPROACHES its target exponentially (shape kSweepShape: the rate falls with the distance left,
        so the tones nearest the corner, which turn most, see the slowest travel), and the rungs are long enough
        that the ceiling allpass's own small phase (≈ 2√2·f/ceiling) does not bend on the first rung. Measured
        (EnvelopeFilterTests, 44.1/48/96 kHz): at most ~48 cents on a wake or a disengage, for 250 ms of ladder;
        a click on the first rung loses 0.3–0.7 dB. Longer halves it again, at the cost of a slower band pass — the
        owner's call by ear (the render tool's door demo). */
    static constexpr double kCornerSweepSeconds = 0.200;
    static constexpr double kRungSeconds        = 0.025;  ///< the ladder's first and last crossfades
    static constexpr double kLadderOctaves      = 4.0;    ///< the ladder's ceiling above the corner (capped at 0.45 fs)
    static constexpr double kSweepShape         = 4.0;    ///< the exponential approach's rate (e-folds over the sweep)

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
        doorLen_ = std::max (1, (int) (kDoorSeconds * fs_ + 0.5));
        rungLen_ = std::max (1, (int) (kRungSeconds * fs_ + 0.5));
        sweepLen_ = std::max (1, (int) (kCornerSweepSeconds * fs_ + 0.5));
        kfMax_   = core::det::log2 (kMaxCutoffShare * fs_);
        apSvf_.prepare (fs_, maxCh_);
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

    /** 0 = dry: the doors close and the stage goes idle; 1 = all filter. */
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

                if (lpPos_ == 0 && ladPos_ == 0 && ! (mixTgt_ > 0.0))
                {
                    goIdle();   // the last door has landed shut: from this sample on the call is dry
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
        apSvf_.healPoison();
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

    /** The doors, in samples — the tests read them to see every edit settle. Audio thread only. The low-pass door
        runs 0 … doorLength(), the band pass's ladder 0 … ladderLength(); steady: bypass idle, low pass
        (doorLength, 0), band pass (0, ladderLength). */
    int lowPassDoor()  const noexcept { return lpPos_; }
    int ladder()       const noexcept { return ladPos_; }
    int doorLength()   const noexcept { return doorLen_; }
    int ladderLength() const noexcept { return 2 * rungLen_ + sweepLen_; }

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
    static constexpr double kApQ       = 0.70710678118654752440;   // the LR4 pair's sum is this allpass
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

    /** The ladder's allpass — the LR4 pair's own (2nd order, Q 1/√2) — at the ceiling while the ladder is below
        its sweep, at the clean-lows corner above it, eased (smoothstep, on log2 Hz) in between. */
    void designAllpass() noexcept
    {
        const int    into = std::min (std::max (ladPos_ - rungLen_, 0), sweepLen_);
        const double c = (double) into / (double) sweepLen_;
        const double e = ease (c);
        const double span = std::min (kLadderOctaves, std::max (kfMax_ - kf_, 0.0));   // the ceiling, above the corner
        apSvf_.setParams (eq::FilterType::AllPass, core::det::exp2 (kf_ + (1.0 - e) * span), kApQ, 0.0);
    }

    /** The sweep's shape, 0 -> 1: an exponential approach, normalised to land on 1 — fast while the corner is far
        from its target, slow as it arrives. Per sample only while the ladder sweeps. */
    static double ease (double c) noexcept
    {
        return (1.0 - detail::detExp (-kSweepShape * c)) / (1.0 - detail::detExp (-kSweepShape));
    }

    void refreshDoors() noexcept
    {
        lp_   = (double) lpPos_ / (double) doorLen_;
        g_    = (double) std::min (ladPos_, rungLen_) / (double) rungLen_;
        beta_ = (double) std::max (ladPos_ - rungLen_ - sweepLen_, 0) / (double) rungLen_;
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

        const bool engaged = mixTgt_ > 0.0;
        bool lpMoved = false, ladMoved = false;
        detail::stepDoors (lpPos_, ladPos_, doorLen_, rungLen_, sweepLen_, engaged, lowPass_, lpMoved, ladMoved);
        if (ladMoved || lpMoved)
            refreshDoors();
        if (ladMoved || kfMoved)
            designAllpass();

        const float env = follower_.process (k);
        const bool lfMoved = glide (lf_, cutoffTarget (env), smoothC_);
        if (lfMoved || qMoved)
            design();

        return mixMoved || sensMoved || keepMoved || qMoved || dirMoved || loMoved || hiMoved || kfMoved
            || ladMoved || lpMoved || lfMoved;
    }

    /** One channel's sample: both modes' splits and filters, the clean lows back, the mix, the doors. A non-finite input is not this stage's to change, and not for its state to hear: everything
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

        // THE DOORS: the low-pass door from x, the ladder's allpass over that, the band pass across.
        const double mixedLp = x + mix_ * (wetLp - x);
        const double mixedBp = dryBp + mix_ * (wetBp - dryBp);
        const double rx = x + lp_ * (mixedLp - x);
        const double ap = (double) apSvf_.processSample (ch, (float) rx);
        const double xr = rx + g_ * (ap - rx);
        return finite ? (float) (xr + beta_ * (mixedBp - xr)) : xin;
    }

    void stopChannel (int ch) noexcept
    {
        if (ch < 0 || ch >= kMaxChannels || ! live_[ch])
            return;
        svf_.clear (ch);
        svf_.clear (2 + ch);
        splitLp_.resetChannel (ch);
        splitBp_.resetChannel (ch);
        apSvf_.resetChannel (ch);
        live_[ch] = false;
    }

    void clearAll() noexcept
    {
        follower_.reset();
        splitLp_.reset();
        splitBp_.reset();
        apSvf_.reset();
        for (int c = 0; c < detail::SweptSvf::kLanes; ++c)
            svf_.clear (c);
        for (int ch = 0; ch < kMaxChannels; ++ch)
            live_[ch] = false;
    }

    /** The glides on their targets, the follower empty, the cutoff at rest, the doors where the mode puts them —
        a fresh stream. */
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
        idle_   = ! (mixTgt_ > 0.0);
        lpPos_  = ! idle_ && lowPass_ ? doorLen_ : 0;
        ladPos_ = ! idle_ && ! lowPass_ ? 2 * rungLen_ + sweepLen_ : 0;
        refreshDoors();
        designSplits();
        designAllpass();
        follower_.reset();
        lf_ = cutoffTarget (0.0f);
        design();
        grid_.reset();
        primed_ = true;
    }

    /** Out of idle: a RESTART — every state from silence, the grid re-anchored, the glides on their targets —
        with the doors opening from dry and the mix gliding in from 0, so the filters'
        own start from silence is faded in rather than heard. */
    void wake() noexcept
    {
        clearAll();
        prime();
        lpPos_  = 0;
        ladPos_ = 0;
        refreshDoors();
        designAllpass();
        mix_    = 0.0;
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
        apSvf_.flushDenormals();
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
    int    lpPos_ = 1, ladPos_ = 0, doorLen_ = 1, rungLen_ = 1, sweepLen_ = 1;   // the doors, counted
    double lp_ = 1.0, g_ = 0.0, beta_ = 0.0;                        // ...as weights
    double kfMax_ = 14.4;   // log2 of 0.45 fs, the ladder allpass's highest ceiling (set by prepare)

    double fs_ = 48000.0, fMax_ = 21600.0;
    double glideC_ = 1.0, smoothC_ = 1.0;

    dynamics::EnvelopeFollower      follower_;
    eq::DeterministicSvf            splitLp_;   // the low-pass mode's complement split
    eq::DeterministicCrossover2     splitBp_;   // the band-pass mode's LR4 split
    eq::DeterministicSvf            apSvf_;     // the ladder's allpass (the LR4 pair's own)
    detail::SweptSvf                svf_;
    core::StateGrid                 grid_;

    bool live_[kMaxChannels] { false, false };   // carried by a call since the last clear (law 11a)
    bool primed_ = false, idle_ = false, cleared_ = true;
    bool prepared_ = false;
    int  maxBlock_ = 0, maxCh_ = 0;
};

} // namespace felitronics::envfilter

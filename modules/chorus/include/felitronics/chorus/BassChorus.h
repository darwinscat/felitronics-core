// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/core/DetMath.h>
#include <felitronics/core/FlushToZero.h>
#include <felitronics/core/StateGrid.h>
#include <felitronics/eq/Svf.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

namespace felitronics::chorus
{

//==============================================================================
// felitronics::chorus::BassChorus — a bass chorus: the dry plus a short, slowly swinging copy of itself,
// with the bass rule that the LOW END stays dry, centred and in phase. The one-slider law is ChorusLaw.h.
//
// THE DRY IS NEVER SPLIT. A chorus is the dry plus a VOICE, and the bass rule is a rule about the voice: it
// is made from the input's Linkwitz–Riley high band (two Butterworth high passes at the crossover, LR4's
// high half). The dry goes out as it came in, sample for sample, at every setting — there is no crossover in
// its path, so it has no crossover phase either, and no switch from "split" to "unsplit" when the effect comes
// in. What the voice takes is LR4's high half and what it leaves is LR4's low half: the two sum flat (LP4 + HP4
// is an allpass), which BassChorusTests measures on the stage's own voice.
//
// ATTENUATED, NOT EXCLUDED. The voice's band falls at LR4's 24 dB/oct under the crossover: -6 dB at it,
// -24.6 dB an octave below, -48 two octaves below. So the low end is kept out of the voice by that much, not
// absolutely, and the crossover's floor is 100 Hz: there the low B and E fundamentals (30.9, 41.2 Hz) enter the
// voice at least 31 dB down, and the Side energy an octave under the crossover stays under -30 dB of the dry's
// at every menu extreme (measured in BassChorusTests, the swing's sidebands included).
//
// THE VOICE. The high band is written into a line and read back by a head that swings around `delayMs` by
// `depthMs` on a SINE: the swing is the chorus's pitch wobble, peak frequency ratio 1 ± 2π·rate·depth. Sine
// rather than triangle because a triangle's corners are pitch STEPS, heard on a sustained bass note as a hiccup
// at every turn. The read is a 4-point cubic Hermite (Catmull–Rom), in double — -0.04 dB at 5 kHz and -0.5 dB
// at 10 kHz at 48 kHz in the worst fraction, where a linear read loses -0.47 and -2.0 dB and its loss swings
// with the LFO. On a stereo call the two channels carry two voices whose LFOs are in QUADRATURE (90°, as
// reverb::AmpReverb's tail): opposite phases would spread wider, but their mono sum is a standing comb that
// breathes like a flanger, while quadrature voices still swing when a PA or a DI sums the pair. A mono call is
// one voice — an adapter that wants two puts the mono signal in two DISTINCT buffers. The voice is ADDED: `mix`
// 1 is the voice at the dry's level, so the voiced band rises ~+3 dB (up to +6 on coherent peaks) at full mix.
//
// IDLE. At mix 0 — the slider's amount 0 — the stage is IDLE once its mix ramp has landed: process() returns
// with the buffer untouched and computes nothing. A WAKE from idle is a RESTART: the stage comes back exactly as
// reset() leaves it (lines and filters cleared, LFO, write clock and law-8 grid at zero, every knob on its
// target) and its voice opens through the 10 ms door.
//
// CONTRACT (DSP-ARCHITECTURE.md §2).
//   · prepare(sampleRate, maxBlock, maxChannels) allocates the lines; process()/reset()/the setters never
//     allocate, lock or throw. maxChannels is 1 or 2 and the rate 8 kHz … 768 kHz — anything else is
//     REFUSED, and a refused prepare DISARMS the instance until a prepare succeeds (law 11b). maxBlock (>= 1)
//     is recorded, not a limit: the stage has no block-sized scratch, so any n >= 0 is processed in full.
//   · process() is [[nodiscard]] bool (law 11): a negative extent, an unprepared instance or a width above
//     the prepared one is refused before anything moves; n == 0 is the only no-op.
//   · A NARROWER call is legal; a channel at an index >= nch with n > 0 has STOPPED (law 11a): its line, its
//     filter state and its DOOR are dropped once, on that edge, so it comes back as a new stream — its voice
//     fades in over the door again, while the lane that kept playing and the shared LFO go on untouched. nch == 0
//     with n > 0 is a gap: every lane stops, and the shared clocks — the LFO, the glides, the mix ramp, the law-8
//     grid — move exactly as through n samples of silence (law 11c).
//   · reset() is a stream RESTART (see IDLE). Cheap-idempotent: the clear runs once.
//   · GLIDES: no knob ever steps the head, which is what a click in a modulated line is. The head (depth and
//     centre) glides CRITICALLY DAMPED — two landing one-poles in cascade, 75 and 80 ms each — so its SPEED
//     starts and ends at zero too: a one-pole starts at full speed, a pitch jump (6 % for a 3 ms depth gliding at
//     50 ms, measured 1.1x to 3.7x past the smooth bound in BassChorusTests). Rate (50 ms) and crossover (30 ms)
//     are one-poles: neither moves the head's speed by a step; the crossover's coefficients follow per sample.
//     The swing is held 1 ms clear of the write head by a SMOOTH minimum (a 1 ms quadratic knee), so a centre
//     gliding into or out of that limit bends the head instead of kinking its speed (a hard min was a ~118-cent
//     step at an LFO crest). The centre's floor, 4 ms, is the law's full depth plus that 1 ms: the slider is
//     never capped by it.
//   · MIX is a COUNTED linear ramp of kMixRampSeconds that lands exactly — on 0, which is what idles the stage.
//   · THE DOOR: each lane's voice INPUT opens linearly over kDoorSeconds after prepare(), reset(), a wake or the
//     lane's own return, so the filters and the line start from a fade instead of from a step in mid-note.
//   · Law 8: the filters' state is flushed on core::StateGrid (audio time, 64 samples); what is written into
//     the line is flushed per sample (core::flushDenormal) and the glides land, so nothing parks in the
//     subnormal band. The line is not recursive: there is no feedback to sustain one.
//   · POISON: a non-finite INPUT sample is not voiced (the dry keeps it — it is not this stage's to change);
//     a filter that overflows writes silence into the line, and its state is healed at the end of the call. The
//     read runs in double and is bounded to the float range, so a line holding a near-FLT_MAX burst replays it
//     as at most ±inf on a finite dry sample, never NaN. Knob values are bounded before they reach arithmetic:
//     a non-finite one is ignored, a finite one clamped.
//   · THREADS: prepare, reset, the setters, params() and process belong to ONE thread (whichever owns the
//     stage). What a face prints comes from the law's pure functions (ChorusLaw.h), not from the stage.
//   · Math: core::det (sin and cos for the LFO — per sample, shared by the voices; exp2 for the glide
//     coefficients) and the eq module's DetMath Svf for the crossover.
//==============================================================================

/** The knobs, in the units the stage takes. `BassChorus::setParams` applies them through the same guarded
    setters a caller may use one by one. ChorusLaw.h says which belong to the slider and which to the menu. */
struct BassChorusParams
{
    float depthMs     = 1.5f;     ///< peak swing of the head around delayMs, ms (0 … 10; eased to <= delayMs - 1)
    float mix         = 0.5f;     ///< the voice's level: 0 = the dry alone (the stage idles), 1 = the dry's level
    float rateHz      = 0.7f;     ///< the LFO, Hz (0.05 … 10)
    float crossoverHz = 200.0f;   ///< below it the voice falls at 24 dB/oct, Hz (100 … 1000)
    float delayMs     = 8.0f;     ///< the centre of the head, ms (4 … 25)
};

namespace detail
{
    inline float clampf (float lo, float hi, float v) noexcept { return v < lo ? lo : (hi < v ? hi : v); }

    /** Equality without -Wfloat-equal: the values compared here are finite by construction. */
    inline bool same (double a, double b) noexcept { return ! (a < b) && ! (b < a); }

    inline constexpr double kLog2E = 1.4426950408889634074;
    inline constexpr double kTwoPi = 6.283185307179586476925286766559;

    /** e^x, the same bits on every row. */
    inline double detExp (double x) noexcept { return core::det::exp2 (core::det::mul (x, kLog2E)); }

    /** min(a, b) with a quadratic knee of width `k`: exactly the min once a and b are `k` apart, never above it,
        and continuous in its first derivative — so a limit that is reached while both sides glide bends the
        result instead of kinking it. */
    inline double smoothMin (double a, double b, double k) noexcept
    {
        const double h = std::max (k - std::fabs (a - b), 0.0) / k;
        return std::min (a, b) - h * h * k * 0.25;
    }

    /** A one-pole glide in double that LANDS: once the residual is under `land` (in the glide's own units, far
        below anything audible) or a step no longer moves the value, the value is put on its target — the two
        conditions core::Smoother documents. `keep` is the residual's factor per step, e^(-1/(tau·fs)). */
    struct Glide
    {
        double cur = 0.0, tgt = 0.0, keep = 0.0, land = 1.0e-15;

        void snap() noexcept { cur = tgt; }
        bool moving() const noexcept { return ! same (cur, tgt); }

        void step() noexcept
        {
            if (! moving())
                return;
            const double d  = keep * (cur - tgt);
            const double nx = tgt + d;
            cur = (std::fabs (d) < land || same (nx, cur)) ? tgt : nx;
        }
    };

    /** Two landing one-poles in cascade: a critically damped glide, whose speed starts and ends at zero. `cur` is
        the output, `mid` the first stage, `tgt` where both are going. */
    struct Glide2
    {
        double cur = 0.0, mid = 0.0, tgt = 0.0, keep = 0.0, land = 1.0e-15;

        void snap() noexcept { cur = mid = tgt; }
        bool moving() const noexcept { return ! same (cur, tgt) || ! same (mid, tgt); }

        void step() noexcept
        {
            if (! moving())
                return;
            Glide first  { mid, tgt, keep, land };
            first.step();
            mid = first.cur;
            Glide second { cur, mid, keep, land };
            second.step();
            cur = second.cur;
        }
    };

    /** A COUNTED linear ramp: `len` samples from wherever it stands to its target, landing on it exactly. A new
        target restarts it from where it stands; the same target leaves it alone. */
    struct Ramp
    {
        double cur = 0.0, to = 0.0, inc = 0.0;
        int    left = 0, len = 1;

        void snap() noexcept { cur = to; left = 0; }
        bool moving() const noexcept { return left > 0; }

        void setTarget (double t) noexcept
        {
            if (same (t, to))
                return;
            to   = t;
            left = len;
            inc  = (to - cur) / (double) len;
        }

        void step() noexcept
        {
            if (left <= 0)
                return;
            --left;
            cur = left == 0 ? to : cur + inc;
        }
    };
} // namespace detail

//==============================================================================
class BassChorus
{
public:
    static constexpr int    kMaxChannels    = 2;
    static constexpr double kMinSampleRate  = 8000.0, kMaxSampleRate = 768000.0;
    static constexpr float  kHeadroomMs     = 1.0f;    ///< the head is held this far clear of the write head
    static constexpr float  kCapKneeMs      = 1.0f;    ///< the width of the smooth minimum that holds it there
    static constexpr float  kMinDelayMs     = 4.0f,  kMaxDelayMs     = 25.0f;   ///< floor = the law's full depth + headroom
    static constexpr float  kMaxDepthMs     = 10.0f;
    static constexpr float  kMinRateHz      = 0.05f, kMaxRateHz      = 10.0f;
    static constexpr float  kMinCrossoverHz = 100.0f, kMaxCrossoverHz = 1000.0f;

    /** How long each lane's voice input takes to open after prepare(), reset(), a wake or the lane's return. */
    static constexpr double kDoorSeconds = 0.010;

    /** The mix's counted ramp. */
    static constexpr double kMixRampSeconds = 0.020;

    /** The glides' time constants, seconds — for depth and centre, of EACH of the two stages. */
    static constexpr double kDepthGlide = 0.075, kDelayGlide = 0.080, kRateGlide = 0.050, kCrossoverGlide = 0.030;

    /** HOW LONG THE VOICE GOES ON after the input stops, at any setting, down to -100 dB of the input's peak: the
        farthest the head reaches (35 ms) plus the high pass's ring at the lowest crossover. Two identical
        Butterworth sections are a double pole pair whose ring falls as (1 + t/τ)·e^(-t/τ), τ = √2/(2π·fc) —
        2.25 ms at 100 Hz — and that is under 2e-6 by 16τ (36 ms): 71 ms in all. Measured at the lowest crossover
        and through a crossover glide in BassChorusTests. */
    static constexpr float tailSeconds() noexcept
    {
        constexpr double tau = 1.4142135623730951 / (detail::kTwoPi * (double) kMinCrossoverHz);
        return (float) ((double) (kMaxDelayMs + kMaxDepthMs) * 0.001 + 16.0 * tau);
    }

    static constexpr int latencySamples() noexcept { return 0; }

    BassChorus() noexcept = default;

    /** Allocates the lines for `maxChannels` (1 or 2) at `sampleRate` (8 kHz … 768 kHz). A refused prepare
        DISARMS the instance. A successful one is a restart (law 11a): no audio fed before it survives. */
    [[nodiscard]] bool prepare (double sampleRate, int maxBlock, int maxChannels)
    {
        prepared_ = false;   // law 11b: whatever happens next, the old preparation no longer stands

        if (! (sampleRate >= kMinSampleRate && sampleRate <= kMaxSampleRate) || maxBlock < 1
            || maxChannels < 1 || maxChannels > kMaxChannels)
            return false;

        fs_       = sampleRate;
        maxBlock_ = maxBlock;
        maxCh_    = maxChannels;

        // A power of two, so a read wraps with a mask: the farthest head plus the interpolator's reach.
        const auto need = (std::size_t) std::ceil ((double) (kMaxDelayMs + kMaxDepthMs) * fs_ / 1000.0) + 8u;
        std::size_t len = 64;
        while (len < need)
            len <<= 1;

        for (std::size_t ch = 0; ch < (std::size_t) kMaxChannels; ++ch)
            line_[ch].assign ((int) ch < maxCh_ ? len : 0u, 0.0f);

        mask_    = (int) len - 1;
        maxHead_ = (double) (len - 4u);

        hp_.prepare (fs_, 2 * maxCh_);   // lane ch is the first section, lane ch + maxCh_ the second

        // The landings, each far under what its knob can make heard: the heads 1e-9 of a sample, the rate
        // 1e-12 Hz, the crossover 1e-3 Hz (a landing there ends ~200 ms of per-sample redesign sooner).
        depth_ = { 0.0, 0.0, 0.0, keepFor (kDepthGlide), 1.0e-9 };
        delay_ = { 0.0, 0.0, 0.0, keepFor (kDelayGlide), 1.0e-9 };
        rate_  = { 0.0, 0.0, keepFor (kRateGlide),      1.0e-12 };
        xover_ = { 0.0, 0.0, keepFor (kCrossoverGlide), 1.0e-3 };
        mix_.len  = std::max (1, (int) (kMixRampSeconds * fs_ + 0.5));
        headroom_ = (double) kHeadroomMs * fs_ / 1000.0;
        knee_     = (double) kCapKneeMs * fs_ / 1000.0;

        doorLen_  = std::max (1, (int) (kDoorSeconds * fs_ + 0.5));
        doorStep_ = 1.0f / (float) doorLen_;

        retarget();
        dirty_ = true;   // the lines were just built; restart() clears them only when they were written
        restart();

        prepared_ = true;
        return true;
    }

    bool isPrepared()  const noexcept { return prepared_; }
    int  maxChannels() const noexcept { return prepared_ ? maxCh_ : 0; }
    int  maxBlock()    const noexcept { return prepared_ ? maxBlock_ : 0; }

    /** A stream restart: the lines and filters cleared, the write clock, the LFO and the grid back at zero,
        every knob on its target, every door shut (each opens over the next kDoorSeconds of its lane). */
    void reset() noexcept { restart(); }

    /** All five knobs. For the face's slider: `setParams (paramsForAmount (amount, params()))` (ChorusLaw.h). */
    void setParams (const BassChorusParams& p) noexcept
    {
        setDepthMs (p.depthMs);
        setMix (p.mix);
        setRateHz (p.rateHz);
        setCrossoverHz (p.crossoverHz);
        setDelayMs (p.delayMs);
    }

    void setDepthMs     (float ms) noexcept { if (std::isfinite (ms)) { depthMs_ = detail::clampf (0.0f, kMaxDepthMs, ms); retarget(); } }
    void setMix         (float m)  noexcept { if (std::isfinite (m))  { mix01_   = detail::clampf (0.0f, 1.0f, m); retarget(); } }
    void setRateHz      (float hz) noexcept { if (std::isfinite (hz)) { rateHz_  = detail::clampf (kMinRateHz, kMaxRateHz, hz); retarget(); } }
    void setCrossoverHz (float hz) noexcept { if (std::isfinite (hz)) { xoverHz_ = detail::clampf (kMinCrossoverHz, kMaxCrossoverHz, hz); retarget(); } }
    void setDelayMs     (float ms) noexcept { if (std::isfinite (ms)) { delayMs_ = detail::clampf (kMinDelayMs, kMaxDelayMs, ms); retarget(); } }

    /** The knobs as set (clamped) — the owning thread's, like the setters. */
    BassChorusParams params() const noexcept { return { depthMs_, mix01_, rateHz_, xoverHz_, delayMs_ }; }

    /** True while the stage computes nothing: mix is 0 and its ramp has landed there. */
    bool isIdle() const noexcept { return idle_; }

    /** In place. Channels at index >= 2 are never part of a call (a wider call is refused). */
    [[nodiscard]] bool process (float* const* io, int numChannels, int numSamples) noexcept
    {
        if (numChannels < 0 || numSamples < 0 || ! prepared_ || numChannels > maxCh_)
            return false;   // law 11: refused before anything moves
        if (numSamples == 0)
            return true;    // the only true no-op: no time, no edge
        if (numChannels > 0 && io == nullptr)
            return false;

        if (idle_)
        {
            if (! (mix_.to > 0.0))
            {
                grid_.skip (numSamples);   // amount 0: nothing computed, but audio time still passed
                return true;
            }
            restart();                     // a wake is a restart
        }

        // Law 11a — the lanes this call does not carry have stopped: their memory and their door go, once.
        for (int ch = numChannels; ch < maxCh_; ++ch)
            stopChannel (ch);
        for (int ch = 0; ch < numChannels; ++ch)
            live_[(std::size_t) ch] = true;
        if (numChannels > 0)
            dirty_ = true;

        int done = 0;
        while (done < numSamples)
        {
            const int seg = grid_.segment (numSamples - done);
            run (io, numChannels, done, seg);
            done += seg;
            if (grid_.advance (seg))
                hp_.flushDenormals();      // law 8, on audio time
        }

        // The poison half of the filters' flush runs on every call (law 8a (a)); a finite state is untouched.
        hp_.healPoison();

        if (! mix_.moving() && ! (mix_.to > 0.0))
            idle_ = true;   // landed on zero: from the next call on, nothing is computed

        return true;
    }

    /** True when the line or a glide holds a subnormal — the law-8 witness for the tests. (The Svf's integrators
        are the eq module's, under its own flush contract.) */
    bool anySubnormalState() const noexcept
    {
        const auto sub = [] (auto v) { return std::fpclassify (v) == FP_SUBNORMAL; };
        for (const auto& ln : line_)
            for (float v : ln)
                if (sub (v))
                    return true;
        return sub (mix_.cur) || sub (depth_.cur) || sub (depth_.mid) || sub (delay_.cur) || sub (delay_.mid)
            || sub (rate_.cur) || sub (xover_.cur);
    }

private:
    double keepFor (double tauSeconds) const noexcept { return detail::detExp (-1.0 / (tauSeconds * fs_)); }

    /** The ramps' and glides' targets from the knobs, in the units the loop runs in. */
    void retarget() noexcept
    {
        mix_.setTarget ((double) mix01_);
        depth_.tgt = (double) depthMs_ * fs_ / 1000.0;
        delay_.tgt = (double) delayMs_ * fs_ / 1000.0;
        rate_.tgt  = (double) rateHz_;
        xover_.tgt = (double) xoverHz_;
    }

    /** Both sections share one design, so one Svf carries them as two lanes per channel: one tan per update. */
    void applyCrossover() noexcept
    {
        constexpr double kQ = 0.7071067811865476;   // Butterworth: two of them in cascade are LR4's high half
        hp_.setParams (eq::FilterType::HighPass, xover_.cur, kQ, 0.0);
    }

    /** prepare(), reset() and a wake: what follows does not depend on anything that came before. */
    void restart() noexcept
    {
        if (dirty_)
            for (auto& ln : line_)
                std::fill (ln.begin(), ln.end(), 0.0f);
        dirty_ = false;
        hp_.reset();
        live_.fill (false);
        doorPos_.fill (0);
        writePos_ = 0;
        phase_    = 0.0;
        mix_.snap(); depth_.snap(); rate_.snap(); delay_.snap(); xover_.snap();
        applyCrossover();
        grid_.reset();
        idle_ = ! (mix_.to > 0.0);
    }

    void stopChannel (int ch) noexcept
    {
        if (ch < 0 || ch >= kMaxChannels || ! live_[(std::size_t) ch])
            return;
        std::fill (line_[(std::size_t) ch].begin(), line_[(std::size_t) ch].end(), 0.0f);
        hp_.resetChannel (ch);
        hp_.resetChannel (ch + maxCh_);
        doorPos_[(std::size_t) ch] = 0;   // it comes back as a new stream, through its door
        live_[(std::size_t) ch] = false;
    }

    /** Samples [off, off + n) of the call, n within one grid segment. nch == 0 runs only the shared clocks. */
    void run (float* const* io, int nch, int off, int n) noexcept
    {
        constexpr double kFloatMax = (double) std::numeric_limits<float>::max();
        const int    mask   = mask_;
        const int    lane2  = maxCh_;
        const double invFs  = 1.0 / fs_;
        const bool   xoverMoved = xover_.moving();

        for (int i = off; i < off + n; ++i)
        {
            // The shared clocks: the glides, the mix ramp, the LFO — one per stage, so the voices never drift
            // apart. The crossover's coefficients follow its glide sample by sample: a coarser step is a small
            // jump in the voice, and on a tone near the corner it is heard as a tick.
            mix_.step(); depth_.step(); rate_.step(); delay_.step();
            if (xover_.moving())
            {
                xover_.step();
                if (nch > 0)
                    applyCrossover();
            }
            phase_ += rate_.cur * invFs;
            if (phase_ >= 1.0)
                phase_ -= 1.0;

            if (nch == 0)
                continue;   // a gap: nothing to voice, the clocks only

            const float  mix   = (float) mix_.cur;
            const double swing = detail::smoothMin (depth_.cur, delay_.cur - headroom_, knee_);
            const double arg   = detail::kTwoPi * phase_;
            const double lfo[kMaxChannels] = { core::det::sin (arg), nch > 1 ? core::det::cos (arg) : 0.0 };

            for (int ch = 0; ch < nch; ++ch)
            {
                float* const x  = io[ch];
                float* const ln = line_[(std::size_t) ch].data();
                const float  in = x[i];
                int&         door = doorPos_[(std::size_t) ch];

                // The voice's input: LR4's high half of the sample, through the lane's door. A sample that is not
                // a number is not voiced — the dry keeps it.
                float v = std::isfinite (in) ? in : 0.0f;
                if (door < doorLen_)
                    v *= (float) ++door * doorStep_;
                float h = hp_.processSample (ch + lane2, hp_.processSample (ch, v));
                if (! std::isfinite (h))
                    h = 0.0f;               // an overflowing filter writes silence; its state heals at the call's end
                core::flushDenormal (h);    // law 8: the only way into the line
                ln[writePos_] = h;

                // The head: the centre plus the swing, the second voice a quarter period on.
                const double head  = std::clamp (delay_.cur + swing * lfo[ch], 2.0, maxHead_);
                const int    whole = (int) head;
                const double t  = head - (double) whole;
                const int    r  = writePos_ - whole;
                const double p0 = ln[(r + 1) & mask];   // one sample newer than the head
                const double p1 = ln[r & mask];
                const double p2 = ln[(r - 1) & mask];
                const double p3 = ln[(r - 2) & mask];

                // 4-point cubic Hermite (Catmull–Rom) from p1 toward p2, in double: a line near FLT_MAX would
                // overflow these sums in float and come out NaN.
                const double c1 = 0.5 * (p2 - p0);
                const double c2 = p0 - 2.5 * p1 + 2.0 * p2 - 0.5 * p3;
                const double c3 = 0.5 * (p3 - p0) + 1.5 * (p1 - p2);
                const double voice = std::clamp (((c3 * t + c2) * t + c1) * t + p1, -kFloatMax, kFloatMax);

                if (mix > 0.0f)
                    x[i] = in + mix * (float) voice;   // at mix 0 the dry is not even added to: its bits stay its own
            }

            writePos_ = (writePos_ + 1) & mask;
        }

        if (nch == 0)
        {
            writePos_ = (int) (((long long) writePos_ + n) & (long long) mask);
            if (xoverMoved)
                applyCrossover();   // a gap moves the glide sample by sample; the coefficients follow once
        }
    }

    double fs_ = 48000.0;
    int    maxBlock_ = 0, maxCh_ = 0;
    bool   prepared_ = false;

    // The knobs, as set.
    float depthMs_ = 1.5f, mix01_ = 0.5f, rateHz_ = 0.7f, xoverHz_ = 200.0f, delayMs_ = 8.0f;

    // The ramp and the glides, in the loop's units: mix 0..1, depth and delay in samples, rate and crossover in Hz.
    detail::Ramp   mix_;                  // counted: lands on an exact 0, which is what idles the stage
    detail::Glide2 depth_, delay_;        // critically damped: the head's speed never steps
    detail::Glide  rate_, xover_;
    double headroom_ = 48.0;              // kHeadroomMs in samples
    double knee_     = 48.0;              // kCapKneeMs in samples
    double maxHead_  = 60.0;              // the safety clamp on the head, inside the line with the interpolator's reach

    std::array<std::vector<float>, kMaxChannels> line_;
    int    mask_ = 0, writePos_ = 0;
    double phase_ = 0.0;                  // the LFO, 0..1
    eq::BasicSvf<core::DetMath> hp_;      // LR4's high half: lanes 0..maxCh_-1 the first section, the rest the second
    core::StateGrid grid_;

    int   doorLen_ = 1;                   // each lane's door, counted: lands on open exactly
    float doorStep_ = 1.0f;
    std::array<int, kMaxChannels>  doorPos_ { 0, 0 };
    std::array<bool, kMaxChannels> live_ { false, false };   // carried by a call since the last clear (law 11a)
    bool idle_  = true;
    bool dirty_ = false;   // the lines were written since they were last cleared
};

} // namespace felitronics::chorus

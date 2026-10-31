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
#include <vector>

namespace felitronics::chorus
{

//==============================================================================
// felitronics::chorus::BassChorus — a bass chorus: the dry plus a short, slowly swinging copy of itself,
// with the bass rule that the LOW END stays dry, centred and in phase.
//
// THE DRY IS NEVER SPLIT. A chorus is the dry plus a VOICE, and the bass rule is a rule about the voice: it
// is made from the input's Linkwitz–Riley high band (two Butterworth high passes at the crossover, LR4's
// high half), so nothing below the crossover is ever delayed, modulated or spread. The dry goes out as it
// came in, sample for sample, at every setting — there is no crossover in its path, so it has no crossover
// phase either, and no switch from "split" to "unsplit" when the effect comes in. What the voice takes is
// LR4's high half and what it leaves alone is the LR4 low half's band: the two are complementary (LP4 + HP4 is
// an allpass, flat in magnitude), which BassChorusTests measures on the stage's own voice.
//
// THE VOICE. The high band is written into a line and read back by a head that swings around `delayMs` by
// `depthMs` on a SINE: the swing is the chorus's pitch wobble, peak ratio 1 + 2π·rate·depth. Sine rather than
// triangle because a triangle's corners are pitch STEPS, heard on a sustained bass note as a hiccup at every
// turn; the sine's deviation is a smooth cosine. The read is a 4-point cubic Hermite (Catmull–Rom) — -0.04 dB
// at 5 kHz and -0.5 dB at 10 kHz at 48 kHz in the worst fraction, where a linear read loses -0.6 and -2 dB
// and its loss swings with the fraction, i.e. with the LFO, into a 2·rate amplitude flutter of the treble.
// On a stereo call the two channels carry two voices whose LFOs are in opposite phase, so one side's copy is
// sharp while the other's is flat and the pair spreads left and right; the lows stay mono because neither
// voice has any. A mono call is one voice. The voice is ADDED: `mix` 1 is the voice at the dry's level.
//
// THE ONE-SLIDER LAW (`law::paramsFor`): `amount` 0..1 drives depth (0 … 3 ms) and mix (0 … 1) together,
// linearly; rate, crossover and the centre delay are the menu's and the slider does not touch them. Amount 0
// is mix 0, and at mix 0 the stage is IDLE: process() returns with the buffer untouched and computes nothing.
// Coming back from idle the voice restarts from silence — its line and filters cleared, its LFO at phase 0 —
// while mix glides up from 0, so a slider pulled off zero fades the chorus in.
//
// CONTRACT (DSP-ARCHITECTURE.md §2).
//   · prepare(sampleRate, maxBlock, maxChannels) allocates the lines; process()/reset()/the setters never
//     allocate, lock or throw. maxChannels is 1 or 2 and the rate 8 kHz … 768 kHz — anything else is
//     REFUSED, and a refused prepare DISARMS the instance until a prepare succeeds (law 11b). maxBlock (>= 1)
//     is recorded, not a limit: the stage has no block-sized scratch, so any n >= 0 is processed in full.
//   · process() is [[nodiscard]] bool (law 11): a negative extent, an unprepared instance or a width above
//     the prepared one is refused before anything moves; n == 0 is the only no-op.
//   · A NARROWER call is legal; a channel at an index >= nch with n > 0 has STOPPED (law 11a): its line and
//     its filter state are dropped once, on that edge. nch == 0 with n > 0 is a gap: every channel stops,
//     and the shared clocks — the LFO, the five glides, the door, the
//     law-8 grid — move exactly as through n samples of silence (law 11c).
//   · reset() is a stream RESTART: the lines and filters cleared, the write clock and the LFO back at zero,
//     every glide landed on its target. What follows a reset does not depend on what came before.
//     Cheap-idempotent: the clear runs once.
//   · GLIDES: every knob moves through a glide in double that LANDS — a knob never steps the head, which is
//     what a click in a modulated line is. The head (depth, centre) and the mix glide CRITICALLY DAMPED, two
//     one-poles in cascade (each 40, 80 and 10 ms), so their SPEED starts and ends at zero too: a one-pole
//     starts at full speed, and a head that suddenly moves is a pitch jump — 6 % for a 3 ms depth gliding at
//     50 ms, heard as a blip the instant the slider is yanked (measured: a third difference 1.1x to 3.7x past
//     the smooth bound in BassChorusTests). Rate (50 ms) and crossover (30 ms) are one-poles: neither moves
//     the head's speed by a step. The crossover's coefficients follow its glide per sample. The swing is held inside the centre by 1 ms per sample (`depth <= delay - 1 ms`), so no pair of
//     glides can bring the head nearer the write head than that.
//   · THE DOOR: after prepare(), reset() or a wake from idle the voice's INPUT opens linearly over
//     kDoorSeconds, so the filters and the line start from a fade instead of from a step in mid-note.
//   · Law 8: the filters' state is flushed on core::StateGrid (audio time, 64 samples); what is written into
//     the line is flushed per sample (core::flushDenormal) and the glides land, so nothing parks in the
//     subnormal band. The line is not recursive: there is no feedback to sustain one.
//   · POISON: a non-finite INPUT sample is not voiced (the dry keeps it — it is not this stage's to change);
//     a filter that overflows writes silence into the line, and its state is healed at the end of the call.
//     Knob values are bounded before they reach arithmetic: a non-finite one is ignored, a finite one clamped.
//   · THREADS: prepare, reset, the setters and process belong to ONE thread (whichever owns the stage).
//   · Math: core::det (sin for the LFO — one call per sample, shared by both voices; exp2 for the glide
//     coefficients; log2 for the published deviation) and the eq module's DetMath Svf for the crossover.
//==============================================================================

/** The knobs, in the units the stage takes. `BassChorus::setParams` applies them through the same guarded
    setters a caller may use one by one. The defaults are the law's menu at amount 0.5. */
struct BassChorusParams
{
    float depthMs     = 1.5f;     ///< peak swing of the head around delayMs, ms (0 … 10; held <= delayMs - 1)
    float mix         = 0.5f;     ///< the voice's level: 0 = the dry alone (the stage idles), 1 = the dry's level
    float rateHz      = 0.7f;     ///< the LFO, Hz (0.05 … 10) — the menu's
    float crossoverHz = 200.0f;   ///< below it nothing is voiced, Hz (40 … 1000) — the menu's
    float delayMs     = 8.0f;     ///< the centre of the head, ms (2 … 25)
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

    /** A one-pole glide in double that LANDS: once the residual is under `land` (in the glide's own units, far
        below anything audible) or a step no longer moves the value, the value is put on its target — the two
        conditions core::Smoother documents. `keep` is the residual's factor per step, e^(-1/(tau·fs)). */
    struct Glide
    {
        double cur = 0.0, tgt = 0.0, keep = 0.0, land = 1.0e-15;

        void snap() noexcept { cur = tgt; }
        bool moving() const noexcept { return ! same (cur, tgt); }

        void stepBy (double k) noexcept
        {
            if (! moving())
                return;
            const double d  = k * (cur - tgt);
            const double nx = tgt + d;
            cur = (std::fabs (d) < land || same (nx, cur)) ? tgt : nx;
        }

        void step() noexcept { stepBy (keep); }
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
} // namespace detail

//==============================================================================
/** The one-slider law: what the face's single CHORUS slider sets. */
namespace law
{
    inline constexpr float kDepthMsAtFull = 3.0f;   ///< the swing at amount 1: 22.7 cents peak at the 0.7 Hz default
    inline constexpr float kMixAtFull     = 1.0f;   ///< the voice at the dry's level at amount 1 — the classic 50/50

    /** 0..1; a non-finite amount is no effect at all. */
    inline float clampAmount (float amount) noexcept
    {
        return std::isfinite (amount) ? detail::clampf (0.0f, 1.0f, amount) : 0.0f;
    }

    inline float depthMsFor (float amount) noexcept { return kDepthMsAtFull * clampAmount (amount); }
    inline float mixFor     (float amount) noexcept { return kMixAtFull * clampAmount (amount); }

    /** The menu's knobs (rate, crossover, centre delay) as given, with depth and mix from `amount`. */
    inline BassChorusParams paramsFor (float amount, BassChorusParams menu = {}) noexcept
    {
        menu.depthMs = depthMsFor (amount);
        menu.mix     = mixFor (amount);
        return menu;
    }
} // namespace law

//==============================================================================
class BassChorus
{
public:
    static constexpr int    kMaxChannels    = 2;
    static constexpr double kMinSampleRate  = 8000.0, kMaxSampleRate = 768000.0;
    static constexpr float  kMinDelayMs     = 2.0f,  kMaxDelayMs     = 25.0f;
    static constexpr float  kMaxDepthMs     = 10.0f;
    static constexpr float  kHeadroomMs     = 1.0f;   ///< the head never comes nearer the write head than this
    static constexpr float  kMinRateHz      = 0.05f, kMaxRateHz      = 10.0f;
    static constexpr float  kMinCrossoverHz = 40.0f, kMaxCrossoverHz = 1000.0f;

    /** How long the voice's input takes to open after prepare(), reset() or a wake from idle. */
    static constexpr double kDoorSeconds = 0.010;

    /** The glides' time constants, seconds — for mix, depth and centre delay, of EACH of the two stages. */
    static constexpr double kMixGlide = 0.010, kDepthGlide = 0.040, kRateGlide = 0.050,
                            kDelayGlide = 0.080, kCrossoverGlide = 0.030;

    /** The longest the voice sounds after the input stops: the farthest the head can reach. */
    static constexpr float tailSeconds() noexcept { return (kMaxDelayMs + kMaxDepthMs) * 0.001f; }
    static constexpr int   latencySamples() noexcept { return 0; }

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

        // The landings, each far under what its knob can make heard: mix 1e-15 (it must reach an exact 0 for
        // the stage to idle), the heads 1e-9 of a sample, the rate 1e-12 Hz, the crossover 1e-6 Hz.
        mix_   = { 0.0, 0.0, 0.0, keepFor (kMixGlide),   1.0e-15 };
        depth_ = { 0.0, 0.0, 0.0, keepFor (kDepthGlide), 1.0e-9 };
        delay_ = { 0.0, 0.0, 0.0, keepFor (kDelayGlide), 1.0e-9 };
        rate_  = { 0.0, 0.0, keepFor (kRateGlide),      1.0e-12 };
        xover_ = { 0.0, 0.0, keepFor (kCrossoverGlide), 1.0e-6 };
        headroom_ = (double) kHeadroomMs * fs_ / 1000.0;

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
        every glide on its target, the door shut (it opens over the next kDoorSeconds of a running stage). */
    void reset() noexcept { restart(); }

    void setParams (const BassChorusParams& p) noexcept
    {
        setDepthMs (p.depthMs);
        setMix (p.mix);
        setRateHz (p.rateHz);
        setCrossoverHz (p.crossoverHz);
        setDelayMs (p.delayMs);
    }

    /** The face's slider: depth and mix from `law::paramsFor`, the menu's knobs left where they are. A
        non-finite amount is ignored. */
    void setAmount (float amount) noexcept
    {
        if (! std::isfinite (amount))
            return;
        setDepthMs (law::depthMsFor (amount));
        setMix (law::mixFor (amount));
    }

    void setDepthMs     (float ms) noexcept { if (std::isfinite (ms)) { depthMs_ = detail::clampf (0.0f, kMaxDepthMs, ms); retarget(); } }
    void setMix         (float m)  noexcept { if (std::isfinite (m))  { mix01_   = detail::clampf (0.0f, 1.0f, m); retarget(); } }
    void setRateHz      (float hz) noexcept { if (std::isfinite (hz)) { rateHz_  = detail::clampf (kMinRateHz, kMaxRateHz, hz); retarget(); } }
    void setCrossoverHz (float hz) noexcept { if (std::isfinite (hz)) { xoverHz_ = detail::clampf (kMinCrossoverHz, kMaxCrossoverHz, hz); retarget(); } }
    void setDelayMs     (float ms) noexcept { if (std::isfinite (ms)) { delayMs_ = detail::clampf (kMinDelayMs, kMaxDelayMs, ms); retarget(); } }

    BassChorusParams params() const noexcept { return { depthMs_, mix01_, rateHz_, xoverHz_, delayMs_ }; }

    /** The set depth and rate as the voice's peak pitch deviation, in cents (the swing held inside the centre
        delay as the stage holds it): 1200·log2(1 + 2π·rate·depth). What a face may print. */
    float pitchDeviationCents() const noexcept
    {
        const double depthS = (double) std::min (depthMs_, delayMs_ - kHeadroomMs) * 0.001;
        return (float) (1200.0 * core::det::log2 (1.0 + detail::kTwoPi * (double) rateHz_ * depthS));
    }

    /** True while the stage computes nothing: mix is 0 and has landed there. */
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
            if (! (mix_.tgt > 0.0))
            {
                grid_.skip (numSamples);   // amount 0: nothing computed, but audio time still passed
                return true;
            }
            wake();
        }

        // Law 11a — the channels this call does not carry have stopped: their memory goes, once.
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
                maintain();
        }

        // The poison half of the filters' flush runs on every call (law 8a (a)); a finite state is untouched.
        hp_.healPoison();

        if (! mix_.moving() && ! (mix_.tgt > 0.0))
            idle_ = true;   // landed on zero: from the next call on, nothing is computed

        return true;
    }

    /** True when the line or a glide holds a subnormal — the law-8 witness for the tests. */
    bool anySubnormalState() const noexcept
    {
        const auto sub = [] (auto v) { return std::fpclassify (v) == FP_SUBNORMAL; };
        for (const auto& ln : line_)
            for (float v : ln)
                if (sub (v))
                    return true;
        return sub (mix_.cur) || sub (mix_.mid) || sub (depth_.cur) || sub (depth_.mid) || sub (delay_.cur)
            || sub (delay_.mid) || sub (rate_.cur) || sub (xover_.cur);
    }

private:
    double keepFor (double tauSeconds) const noexcept { return detail::detExp (-1.0 / (tauSeconds * fs_)); }

    /** The glides' targets from the knobs, in the units the loop runs in. */
    void retarget() noexcept
    {
        mix_.tgt   = (double) mix01_;
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

    void clearLines() noexcept
    {
        if (dirty_)
            for (auto& ln : line_)
                std::fill (ln.begin(), ln.end(), 0.0f);
        dirty_ = false;
        hp_.reset();
        live_.fill (false);
        writePos_ = 0;   // the clocks re-anchor, so what follows does not depend on how long the stream ran
        phase_    = 0.0;
        doorPos_  = 0;
    }

    /** prepare() and reset(): everything on its target, the mix included. */
    void restart() noexcept
    {
        clearLines();
        mix_.snap(); depth_.snap(); rate_.snap(); delay_.snap(); xover_.snap();
        applyCrossover();
        grid_.reset();
        idle_ = ! (mix_.tgt > 0.0);
    }

    /** Off idle: the voice starts from silence, its knobs on their targets, while the mix glides up from 0. */
    void wake() noexcept
    {
        clearLines();
        depth_.snap(); rate_.snap(); delay_.snap(); xover_.snap();
        applyCrossover();
        mix_.cur = mix_.mid = 0.0;
        idle_ = false;
    }

    void stopChannel (int ch) noexcept
    {
        if (ch < 0 || ch >= kMaxChannels || ! live_[(std::size_t) ch])
            return;
        std::fill (line_[(std::size_t) ch].begin(), line_[(std::size_t) ch].end(), 0.0f);
        hp_.resetChannel (ch);
        hp_.resetChannel (ch + maxCh_);
        live_[(std::size_t) ch] = false;
    }

    /** On the grid: the filters' law-8 flush. */
    void maintain() noexcept { hp_.flushDenormals(); }

    /** Samples [off, off + n) of the call, n within one grid segment. nch == 0 runs only the shared clocks. */
    void run (float* const* io, int nch, int off, int n) noexcept
    {
        const int    mask   = mask_;
        const int    lane2  = maxCh_;
        const double invFs  = 1.0 / fs_;
        const bool   xoverMoved = xover_.moving();

        for (int i = off; i < off + n; ++i)
        {
            // The shared clocks: the glides, the LFO, the door — one per stage, so the voices never drift apart.
            // The crossover's coefficients follow its glide sample by sample: a coarser step is a small jump in
            // the voice, and on a tone near the corner it is heard as a tick.
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
            if (doorPos_ < doorLen_)
                ++doorPos_;

            if (nch == 0)
                continue;   // a gap: nothing to voice, the clocks only

            const bool   doorOpen = doorPos_ >= doorLen_;
            const float  doorGain = (float) doorPos_ * doorStep_;
            const float  mix      = (float) mix_.cur;
            const double swing    = std::min (depth_.cur, delay_.cur - headroom_);
            const double lfo      = core::det::sin (detail::kTwoPi * phase_);

            for (int ch = 0; ch < nch; ++ch)
            {
                float* const x  = io[ch];
                float* const ln = line_[(std::size_t) ch].data();
                const float  in = x[i];

                // The voice's input: LR4's high half of the sample, through the door. A sample that is not a
                // number is not voiced — the dry keeps it.
                const float v = std::isfinite (in) ? (doorOpen ? in : in * doorGain) : 0.0f;
                float h = hp_.processSample (ch + lane2, hp_.processSample (ch, v));
                if (! std::isfinite (h))
                    h = 0.0f;               // an overflowing filter writes silence; its state heals at the call's end
                core::flushDenormal (h);    // law 8: the only way into the line
                ln[writePos_] = h;

                // The head: the centre plus the swing, the second voice on the opposite phase.
                const double head = std::clamp (delay_.cur + (ch == 0 ? swing : -swing) * lfo, 2.0, maxHead_);
                const int    whole = (int) head;
                const float  t  = (float) (head - (double) whole);
                const int    r  = writePos_ - whole;
                const float  p0 = ln[(r + 1) & mask];   // one sample newer than the head
                const float  p1 = ln[r & mask];
                const float  p2 = ln[(r - 1) & mask];
                const float  p3 = ln[(r - 2) & mask];

                // 4-point cubic Hermite (Catmull–Rom) from p1 toward p2.
                const float c1 = 0.5f * (p2 - p0);
                const float c2 = p0 - 2.5f * p1 + 2.0f * p2 - 0.5f * p3;
                const float c3 = 0.5f * (p3 - p0) + 1.5f * (p1 - p2);
                const float voice = ((c3 * t + c2) * t + c1) * t + p1;

                if (mix > 0.0f)
                    x[i] = in + mix * voice;   // at mix 0 the dry is not even added to: its bits stay its own
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

    // The glides, in the loop's units: mix 0..1, depth and delay in samples, rate and crossover in Hz.
    detail::Glide2 mix_, depth_, delay_;   // critically damped: the head's speed and the mix's slope never step
    detail::Glide  rate_, xover_;
    double headroom_    = 48.0;  // kHeadroomMs in samples
    double maxHead_     = 60.0;  // the safety clamp on the head, inside the line with the interpolator's reach

    std::array<std::vector<float>, kMaxChannels> line_;
    int    mask_ = 0, writePos_ = 0;
    double phase_ = 0.0;         // the LFO, 0..1
    eq::BasicSvf<core::DetMath> hp_;   // LR4's high half: lanes 0..maxCh_-1 the first section, the rest the second
    core::StateGrid grid_;

    int   doorLen_ = 1, doorPos_ = 0;   // the door, counted: lands on open exactly
    float doorStep_ = 1.0f;

    std::array<bool, kMaxChannels> live_ { false, false };   // carried by a call since the last clear (law 11a)
    bool idle_  = true;
    bool dirty_ = false;   // the lines were written since they were last cleared
};

} // namespace felitronics::chorus

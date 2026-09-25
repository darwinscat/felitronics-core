// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/stereo/MidSide.h>
#include <felitronics/eq/Crossover2.h>
#include <felitronics/eq/MatchedBiquad.h>
#include <felitronics/core/Smoother.h>
#include <felitronics/core/StateGrid.h>

#include <algorithm>
#include <cstdint>
#include <cmath>

namespace felitronics::stereo
{

//==============================================================================
// felitronics::stereo::MonoBass — "elliptical EQ" / bass mono-maker: collapse the low end to mono below a
// crossover (vinyl/LP lacquer cutting — out-of-phase lows cause vertical stylus excursion — and general
// mastering tightening). Zero added latency, minimum-phase.
//
// TOPOLOGY (verified by math + the self-tests). Operate ONLY on the Side of a Mid/Side split, so the Mid
// (the mono content you keep) passes UNFILTERED — no transient smear on the kept bass; only the
// low-frequency Side (stereo) energy is removed. A full L/R crossover would needlessly filter the mono
// bass too. The side rolls off at 24 dB/oct below fc (−6 dB at fc), so deep below fc the output tends to
// L == R by construction — the tool can only IMPROVE mono compatibility, never harm it.
//
//   M = ½(L+R),  S = ½(L-R)
//   S_wet = lowWidth·LP4(S) + HP4(S)      // LP4/HP4 = eq::Crossover2, a 4th-order Linkwitz-Riley split
//   S'    = xf·S + (1-xf)·S_wet           // xf = the full-wide bypass crossfade (see PARAMS below)
//   L' = M + S',  R' = M − S'
//
// WHY LINKWITZ-RILEY, fixed at LR4 / 24 dB/oct — consilium-settled; don't re-litigate without new math:
//  * The LR4 identity LP4+HP4 = allpass — (1+s⁴)/(s²+√2s+1)², |·| = (1+ω⁴)/(1+ω⁴) = 1, the branches
//    exactly IN-PHASE — is what makes the `lowWidth` taper bump-free: the blended side magnitude is
//    (lowWidth + r⁴)/(1 + r⁴) with r = tan(πf/fs)/tan(πfc/fs) (the BLT warp the TPT SVF realises exactly).
//  * matched::highpass/lowpass (the Nyquist-matched EQ filters) were REJECTED: they are magnitude fits,
//    not a complementary pair — their sum is not that allpass, so the blend would ripple; and Nyquist
//    accuracy buys nothing at a ~120 Hz crossover.
//  * An allpass-subtraction variant was REJECTED: AP4−LP4 has magnitude 1.5 at fc (a boost, not a reject).
//  * A crossfade blend w·S + (1−w)·HP4(S) was REJECTED: S and HP4(S) are anti-phase at fc (HP4(jω₀) = −½),
//    so that blend NULLS the side at fc when w = ⅓.
//  * LR2 / 12 dB/oct was REJECTED: LP2+HP2 = (1+s²)/(s+1)² has a NULL at fc — its flat sum needs a
//    polarity flip (LP−HP), which breaks S' = w·LP + HP. LR8 / 48 dB/oct would be admissible (the in-phase
//    allpass sum holds) but is YAGNI until a product asks for it.
//
// PARAMS. `lowWidth` ∈ [0,1] — 0 = fully mono below fc (default), 1 = full stereo — is LINEAR-smoothed
// (~20 ms) so live automation doesn't click. Full-wide is special: the wet path at lowWidth=1 is the LR4
// allpass — magnitude-flat but phase-rotated (−1 at fc) — so a hard switch to raw bypass there would
// click. Instead lowWidth=1 also ramps the dry/wet crossfade `xf`→1 (same 20 ms), and only once xf has
// SETTLED at 1.0 does process() take the bit-exact early-return bypass. Entering ANY bypass (settled
// full-wide / disabled / non-stereo) resets the crossover state, so re-entry never replays stale tails —
// it restarts the filters from zero, a bounded click-free settle (tested). Setters reject non-finite values
// and clamp (a stray host NaN can't poison the SVF state — std::clamp would pass it through). reset() snaps
// the smoothers to their targets (snap-on-load: no ramp on session recall). `frequency()` reads back the
// clamped value ([20, 0.45·fs], re-clamped if prepare() lowers fs).
//
// LIVE MOVES GLIDE, and the first write of a stream SNAPS. Measured on a 103.7 Hz tone with L and R 90° apart
// (a large Side) through the mastering chain, max|Δ²y| where the change reaches the output against -84 dBFS
// for the steady tone: the crossover frequency 60 -> 250 Hz stepped at -45.3 dBFS (its SVF coefficients were
// set at once), `enabled` false -> true at -48.7 and true -> false at -21.5 (a hard switch with a crossover
// reset). So:
//  * THE CROSSOVER FREQUENCY rides a one-pole `core::Smoother` (kFreqSmoothMs = 30 ms, the EQ's own), advanced
//    by 64 samples at every `core::StateGrid` boundary with the four SVFs redesigned there when it moved —
//    EqBand's idiom, and a TPT SVF is the topology made for a moving coefficient.
//  * `enabled` IS THE xf CROSSFADE: disabling fades the Side to dry exactly as a full-wide lowWidth does, and
//    the bit-exact bypass follows once it has settled; enabling fades back in from a crossover restarted at
//    zero. It used to be a hard toggle "for parity with StereoWidth"; the mastering chain's live preview is
//    the product that disagrees.
//  * THE AIR SHELF'S CORNER rides a Smoother like the crossover's, and its `enabled` glides the plateau to and
//    from 0 dB (where the shelf is skipped) instead of switching it.
// Every one of those advances on audio time — per sample, or at grid boundaries counted in samples, bypassed
// stretches included — so a stream cut any other way renders the same bits (law 8a). And every write that
// lands before the first sample after prepare()/reset() SNAPS, lowWidth and the air plateau included, which
// they did not: a write between prepare() and the first sample used to glide those two from the prepared
// values, so "configure, then prepare" and "prepare, then configure" rendered differently (measured through
// the mastering chain, 71 529 and 59 030 samples of a one-second render). They are one behaviour now — the
// configure-then-prepare one, whose bits are unchanged.
//
// STRICTLY STEREO: process() touches the buffer ONLY when numChannels == 2 — a mono or surround bus
// passes through whole (treating a stereo pair inside a wider layout is a routing decision the host
// owns, not this class). RT-safe: no alloc/lock/throw in process(); state = one eq::Crossover2 (the Side
// is one channel) + two LinearSmoothers.
// THE STAGE OWNS THE TYPE OF ITS OWN PARAMETERS. Added so a composite that drives this stage passes
// one object through instead of copying three fields across, which is the shape that falls out of step
// the first time a fourth setting appears here (the reason `CompressorParams` derives from
// `GainReductionParams` rather than repeating it). Purely additive: `setParams` calls the three
// setters below in the order a caller would, so every clamp, non-finite rejection and smoothing
// decision is unchanged and the result is bit-identical to setting them by hand — pinned in the suite.
struct MonoBassParams
{
    bool  enabled     = true;
    float frequencyHz = 120.0f;
    float lowWidth    = 0.0f;      // 0 = mono below fc, 1 = full stereo (settles into bypass)
};

//==============================================================================
// K14 — "stereo air": a high shelf on SIDE ONLY, riding in the SAME M/S island this class already
// opens. A separate stage would have opened a second one, and the encode/decode round trip is not the
// identity in float (see the note in process()), so a second one would perturb the programme twice for
// no reason. It is here, and not in an EQ band, because an EQ band works on L and R.
//
// WHAT IT DOES TO THE SOUND, stated where a caller will read it:
//  * THE MONO FOLD DOES NOT CHANGE. With m = (l+r)/2 and l = m+s, (l'+r')/2 = m for ANY Side
//    processing, phase included. What grows is the GAP: stereo gains top end, mono does not. On
//    uncorrelated highs a +3 dB plateau takes the stereo-to-mono drop from -3.01 dB to -4.76.
//  * ON ANTI-PHASE HIGHS EVERY WIDTH NUMBER IS BLIND. Width reads 1.000 before and 1.000 after while
//    the Side energy grows by the full plateau, because the fold was already empty. Read the band
//    ENERGIES, not the fraction — which is why the chain publishes all three.
//  * IT BLEEDS A HARD-PANNED TOP INTO THE OTHER CHANNEL. With content only in L, the output R is
//    about -15 dB of L with INVERTED polarity above the corner, because r' = m - H·s. No L/R shelf
//    does this; it is a property of acting on Side, not a defect.
//  * `gainDb` IS THE PLATEAU, not the gain at the corner: the shelf reaches half of it AT `frequencyHz`
//    and the rest above. At a 12 kHz corner and 44.1 kHz even the analogue prototype only reaches
//    +2.75 dB by Nyquist, so the band energy moves by x1.74 rather than x2 for a +3 dB request.
struct StereoAirParams
{
    bool  enabled     = false;     // OFF, so a parameter set written before this existed renders as it did
    float frequencyHz = 6000.0f;   // the shelf corner; clamped to [3000, min(12000, 0.45*fs)]
    float gainDb      = 0.0f;      // the PLATEAU, clamped to [0, kMaxAirDb]; 0 skips the filter entirely
};

class MonoBass
{
public:
    static constexpr double kSmoothingMs = 20.0;    // click-free lowWidth automation + the bypass crossfade
    static constexpr double kFreqSmoothMs = 30.0;   // the crossover and the air corner, one-pole on the grid
    static constexpr float  kMinFreq     = 20.0f;

    // Law 11: a default-constructed MonoBass is NOT a valid configuration, whatever its member defaults
    // suggest — `xo_` has no coefficients until this runs, so the fold does not fold. Measured on a pure
    // SIDE signal at lowWidth 0: a default object passes the side band at -6.02 dB where a prepared one
    // kills it to -54.22 — 48.199 dB of difference at 30 Hz, 6.02 at 120.
    // LAW 11(b): it TAKES a width, so the width is binding — "a module that ignores an argument is
    // lying about its contract" was true of this one. Disarm first, validate, then write.
    [[nodiscard]] bool prepare (double sampleRate, int /*maxBlock*/ = 0, int maxChannels = 2) noexcept
    {
        prepared_ = false;
        if (maxChannels < 1 || maxChannels > 2) return false;
        fs_ = (std::isfinite (sampleRate) && sampleRate > 0.0) ? sampleRate : 48000.0;   // inf/NaN/≤0 -> 48 kHz
        fresh_ = true;                              // a preparation is a stream restart: the corners design at once,
        xoHz_  = -1.0f;                             // and at THIS rate, whatever the last one designed
        xo_.prepare (fs_, 1);
        widthSm_.reset (fs_, kSmoothingMs * 0.001);
        xfSm_.reset (fs_, kSmoothingMs * 0.001);
        freqSm_.prepare (fs_, kFreqSmoothMs);
        airFreqSm_.prepare (fs_, kFreqSmoothMs);
        applyFrequency();
        airSm_.reset (fs_, kSmoothingMs * 0.001);
        clampAirFrequency();                        // the ceiling is 0.45*fs, so it moves with the rate
        wxM_.prepare (fs_, 1); wxSb_.prepare (fs_, 1); wxSa_.prepare (fs_, 1);
        retuneWidthBand();
        reset();
        prepared_ = true;
        return true;
    }

    void reset() noexcept                           // snap smoothers to targets (settled, no glide) + clear filter state
    {
        widthSm_.setCurrentAndTargetValue (lowWidth_);
        xfSm_.setCurrentAndTargetValue (xfTarget());
        // The corners land where they were written and are designed there — a glide in flight does not resume.
        freqSm_.snap ((double) freq_);
        designXo();
        airFreqSm_.snap ((double) airFreq_);
        airFreqCur_ = airFreq_;
        xo_.reset();
        grid_.reset();                              // a stream restart re-anchors the maintenance grid
        bypassed_ = false;
        // K14. The shelf's ramp snaps like the others, its state clears like the crossover's, and the
        // width interval starts over: a total carried across a restart is not a measurement of either
        // side of it.
        airSm_.setCurrentAndTargetValue (airTarget());
        airShelf_.reset();
        airBypassed_ = core::exactlyEqual (airTarget(), 0.0f);
        airDesignedDb_ = airDesignedHz_ = -1.0f;
        retuneWidthBand();
        fresh_ = true;                              // the next writes, up to the first sample, snap
    }

    // Glides through the xf crossfade (see LIVE MOVES GLIDE) — or lands at once before the stream's first sample.
    void setEnabled (bool e) noexcept
    {
        enabled_ = e;
        if (fresh_) xfSm_.setCurrentAndTargetValue (xfTarget());
        else        xfSm_.setTargetValue (xfTarget());
    }

    void setFrequency (float hz) noexcept
    {
        if (! std::isfinite (hz)) return;           // reject NaN/inf — keep the last good value
        freq_ = hz;
        applyFrequency();
    }

    void setLowWidth (float w) noexcept             // 0 = mono below fc, 1 = full stereo (settles into bypass)
    {
        if (! std::isfinite (w)) return;
        lowWidth_ = std::clamp (w, 0.0f, 1.0f);
        if (fresh_)
        {
            widthSm_.setCurrentAndTargetValue (lowWidth_);
            xfSm_.setCurrentAndTargetValue (xfTarget());
            return;
        }
        widthSm_.setTargetValue (lowWidth_);
        xfSm_.setTargetValue (xfTarget());
    }

    // The three setters above as one object — see MonoBassParams. Order matters only in that it is the
    // order a caller would use; none of the three interacts with another.
    void setParams (const MonoBassParams& p) noexcept
    {
        setEnabled   (p.enabled);
        setFrequency (p.frequencyHz);
        setLowWidth  (p.lowWidth);
    }

    // The RESOLVED values, after this class's own clamps and rejections — `frequencyHz` comes back
    // clamped to [20, 0.45*fs] and a non-finite write comes back as the last good value. A composite
    // that has to report what it actually applied reads this, rather than echoing what it was handed.
    MonoBassParams params() const noexcept { return { enabled_, freq_, lowWidth_ }; }

    // K14 — the Side shelf. Same house rule as everything else here: non-finite is REFUSED and the last
    // good value stands, the ranges are clamped, and `air()` reads back what was actually applied.
    void setAir (const StereoAirParams& p) noexcept
    {
        // THE CORNER IS COMPARED AGAINST WHAT IT WAS BEFORE THIS WRITE. It used to be captured after the
        // assignment, so the width-measurement band was re-tuned only when the CLAMP moved the corner, never
        // when the caller did: a live 6000 -> 8000 Hz move kept summing across both bands (measured:
        // airJudgedSamples() read 2000 where the new band had judged 1000). The note at retuneWidthBand() was
        // always the contract; this is the code agreeing with it.
        const float wasFreq = airFreq_;
        airEnabled_ = p.enabled;
        if (std::isfinite (p.frequencyHz)) airFreq_ = p.frequencyHz;
        if (std::isfinite (p.gainDb))      airDb_   = std::clamp (p.gainDb, 0.0f, kMaxAirDb);
        clampAirFrequency();
        if (! core::exactlyEqual (wasFreq, airFreq_)) retuneWidthBand();
        // THE SMOOTHER CARRIES THE dB, not the coefficients: a +3 dB step on Side is a 0.41*S jump, and
        // a shelf redesigned between two settled values is still a discontinuity in the output. `enabled`
        // rides the same smoother now — its target is the plateau, or 0 dB (where the shelf is skipped) — and
        // the corner rides its own, on the grid.
        if (fresh_)
        {
            airSm_.setCurrentAndTargetValue (airTarget());
            airFreqSm_.snap ((double) airFreq_);
            airFreqCur_ = airFreq_;
            return;
        }
        airSm_.setTargetValue (airTarget());
        airFreqSm_.setTarget ((double) airFreq_);
    }
    StereoAirParams air() const noexcept { return { airEnabled_, airFreq_, airDb_ }; }

    // ==============================================================================================
    // K14 — WHAT THE SHELF DID TO THE TOP, measured where it acted and on the band it acted on.
    // ==============================================================================================
    // THREE ENERGIES, NOT A FRACTION, and the reason is a case a fraction cannot report: on anti-phase
    // highs the width reads 1.000 before and 1.000 after while the Side energy grows by the whole
    // plateau, because the mono fold was already empty. From these three a caller gets the width on any
    // convention it likes, the fold loss 10log10(M/(M+S)), and the stereo-level rise — all of which the
    // fraction alone cannot give back.
    //
    // THE BAND IS AN LR4 HIGH-PASS AT THE SHELF'S OWN CORNER, so content AT the corner counts a quarter
    // (-6 dB) and the shelf also acts a little below it: this is a WEIGHTING, not a brick wall, and a
    // caller comparing it with a textbook "energy above 6 kHz" will be wrong by the skirt.
    //
    // ONLY WHILE THE AIR IS ENGAGED. Three LR4 pairs per sample is real work, and with the tool off
    // there is nothing to report anyway — `airJudgedSamples() == 0` says exactly that, and is not a
    // measurement of zero.
    double       airMidEnergy()        const noexcept { return wSumM_; }
    double       airSideEnergyBefore() const noexcept { return wSumSb_; }
    double       airSideEnergyAfter()  const noexcept { return wSumSa_; }
    std::int64_t airJudgedSamples()    const noexcept { return wSamples_; }
    // The PAGE's own convention, named so nobody puts two different "widths" side by side: the amplitude
    // fraction sqrt(S)/(sqrt(M)+sqrt(S)), which is what stereo-meter.js computes. -1.0, never 0.0, when
    // there is no energy to judge — 0.0 is a legitimate reading (an exactly mono top).
    static double airWidth (double mid, double side) noexcept
    {
        if (! (mid >= 0.0) || ! (side >= 0.0) || (mid + side) <= 0.0) return -1.0;
        const double m = std::sqrt (mid), s = std::sqrt (side);
        return (m + s) > 0.0 ? s / (m + s) : -1.0;
    }
    double airWidthBefore() const noexcept { return airWidth (wSumM_, wSumSb_); }
    double airWidthAfter()  const noexcept { return airWidth (wSumM_, wSumSa_); }
    static constexpr float kMaxAirDb = 6.0f;
    static constexpr float kMinAirFreq = 3000.0f, kMaxAirFreq = 12000.0f;

    bool  isEnabled() const noexcept { return enabled_; }
    float frequency() const noexcept { return freq_; }
    float lowWidth()  const noexcept { return lowWidth_; }
    static constexpr int latencySamples() noexcept { return 0; }

    // Stereo, in place: io[0]=L, io[1]=R. Bypass (buffer untouched) when disabled, numChannels != 2, or
    // settled full-wide. Entering bypass clears the crossover so re-entry starts from silence, not stale tails.
    // Law 11 (DSP-ARCHITECTURE.md §2): this is a fixed STEREO stage, so 2 is both its maximum and the
    // only width it can act on; narrower is a documented passthrough, wider is refused.
    [[nodiscard]] bool process (float* const* io, int numChannels, int n) noexcept
    {
        if (numChannels < 0 || n < 0) return false;
        if (! prepared_) return false;                  // see prepare(): the crossover is not built yet
        if (numChannels > 2) return false;              // width is a LIMIT — law 11(b)
        if (n == 0) return true;                        // no samples: no time, no edge, nothing at all —
                                                        // the bypass edge below used to fire even here
        fresh_ = false;                                 // the stream has started: a write glides from here on
        // TWO TOOLS SHARE THIS ISLAND NOW, so the gate asks about BOTH. The island is skipped only when
        // neither has anything to do; with the bass settled full-wide and the air ramping, the round trip
        // still has to run, and the old spelling would have skipped it and dropped the shelf silently.
        // A DISABLED bass is one whose crossfade has settled at dry — `enabled` rides xf now.
        const bool mbIdle = ! xfSm_.isSmoothing() && core::exactlyEqual (xfSm_.getCurrentValue(), 1.0f);
        if (numChannels != 2 || (mbIdle && airIdle()))
        {
            if (! bypassed_)    { bypassed_ = true;    xo_.reset(); }
            if (! airBypassed_) { airBypassed_ = true; airShelf_.reset(); }
            skipTime (n);                   // bypassed audio is still audio TIME — the grid and the corners move on
            return true;
        }
        // …and each tool owns its own latch. Resetting the crossover because the AIR toggled would click
        // the bass; resetting the shelf because the bass settled would click the top.
        bool mbDone = mbIdle;
        if (mbDone) { if (! bypassed_) { bypassed_ = true; xo_.reset(); } }
        else bypassed_ = false;
        float* L = io[0];
        float* R = io[1];
        for (int i = 0; i < n; ++i)
        {
            // SETTLING INTO FULL-WIDE IS A SAMPLE EVENT, not a call event. `xfSm_` arrives inside this
            // loop, and testing for it only at the top of the next call left every sample in between
            // going through the M/S round trip — which is NOT the identity in float (0.5(L+R) + 0.5(L-R)
            // rounds twice), so the output depended on where the caller cut: measured 1 LSB of 24 bit
            // (5.96e-08) on 22953 of 336000 samples across the re-slicing sweep. Same test, same
            // one-shot reset, moved to the clock it belongs on.
            // …and it retires the BASS, not the island: the air may still be working, and leaving the
            // loop here would drop its shelf for the rest of the block — silently, with no refusal and
            // nothing in the report. Only when BOTH are done is there nothing left to do.
            if (! mbDone && ! xfSm_.isSmoothing() && core::exactlyEqual (xfSm_.getCurrentValue(), 1.0f))
            {
                mbDone = true;
                if (! bypassed_) { bypassed_ = true; xo_.reset(); }
                if (airIdle()) { skipTime (n - i); return true; }
            }
            float m, s; MidSide::encode (L[i], R[i], m, s);
            float sOut = s;
            if (! mbDone)
            {
                const float w  = widthSm_.getNextValue();
                const float xf = xfSm_.getNextValue();
                float lp, hp; xo_.processSample (0, s, lp, hp);
                const float wet = w * lp + hp;               // side magnitude (w + r⁴)/(1+r⁴) — bump-free (LR4 in-phase)
                sOut = xf * s + (1.0f - xf) * wet;           // ≠ wet only while fading into/out of full-wide
            }
            // K14. NOT a multiply by one and NOT a unity biquad when the gain is zero: `highShelfDb`
            // substitutes g = 1.00001 for a zero request (it is designed as a +8.7e-5 dB shelf), and even
            // forced identity coefficients are not transparent — the biquad computes 1.0*x + 0.0 in double
            // and turns -0.0f into +0.0f, measured, which is the very defect that disqualified the
            // clipper's mix = 0 bypass. So zero is a BRANCH, taken on the smoothed and clamped value.
            // `enabled` rides the plateau's smoother (its target is 0 dB when off), so the shelf runs for as long
            // as that smoother is away from 0 — a switch-off fades instead of stepping.
            {
                const bool  moving = airSm_.isSmoothing();
                const float db     = airSm_.getNextValue();
                if (moving || ! core::exactlyEqual (db, 0.0f))
                {
                    designAir (db, moving);
                    airBypassed_ = false;
                    sOut = airShelf_.processSample (sOut);
                }
                else if (! airBypassed_) { airBypassed_ = true; airShelf_.reset(); }
            }
            // The measurement, on the band the shelf acts on, taking `s` and `sOut` while both are in hand.
            if (airEnabled_)
            {
                float lp, hp;
                wxM_ .processSample (0, m,    lp, hp); wSumM_  += (double) hp * (double) hp;
                wxSb_.processSample (0, s,    lp, hp); wSumSb_ += (double) hp * (double) hp;
                wxSa_.processSample (0, sOut, lp, hp); wSumSa_ += (double) hp * (double) hp;
                ++wSamples_;
            }
            MidSide::decode (m, sOut, L[i], R[i]);
            // LAW 8 on the AUDIO-TIME grid, not at the end of the call: the flush zeroes state, so
            // putting it where the caller happened to cut made the output a function of the host's
            // block size (measured on this stage: 37180 of 40000 tail samples differ between a
            // whole-file call and one-sample calls, and a whole-file call never flushed at all).
            // One increment and a compare per sample; the flush itself still runs once per period.
            if (grid_.advance (1)) { xo_.flushDenormals(); airShelf_.flushDenormals(); tickCorners(); }
        }
        // The poison half stays per call — see eq::Biquad::healPoison(). Invisible on a finite stream.
        xo_.healPoison();
        airShelf_.healPoison();
        return true;
    }

private:
    // The clamp, then either the design at once (before the stream's first sample) or a new target for the
    // corner's smoother, which the grid walks to it.
    void applyFrequency() noexcept
    {
        const float hi = std::max (kMinFreq, (float) (0.45 * fs_));   // keep lo <= hi even at absurdly low fs (clamp UB otherwise)
        freq_ = std::clamp (freq_, kMinFreq, hi);
        if (fresh_) { freqSm_.snap ((double) freq_); designXo(); }
        else        freqSm_.setTarget ((double) freq_);
    }

    // The crossover at the smoother's value — redesigned only when that value moved, so a settled corner costs
    // a comparison per grid boundary and nothing else.
    void designXo() noexcept
    {
        const float f = (float) freqSm_.value();
        if (core::exactlyEqual (f, xoHz_)) return;
        xoHz_ = f;
        xo_.setFrequency (f);
    }

    // ONE GRID BOUNDARY for the two corners: each smoother advances by exactly one period (EqBand's rule — a
    // partial period must never advance control, or the glide would be back on the caller's clock), and the
    // crossover is redesigned if its corner moved. The shelf picks its corner up at its next design.
    void tickCorners() noexcept
    {
        freqSm_.advance (core::StateGrid::kPeriod);
        airFreqSm_.advance (core::StateGrid::kPeriod);
        designXo();
        airFreqCur_ = (float) airFreqSm_.value();
    }

    // Audio time that ran no audio — a bypassed island, a mono stretch, a clock-only call — still carries the
    // corners: every boundary it crosses is a tick, exactly as the audio path would have taken it. Bounded by the
    // smoothers' own settling (they arrive exactly and then cost nothing), not by the length of the call.
    void skipTime (int n) noexcept
    {
        long long boundaries = ((long long) grid_.phase() + (long long) n) / core::StateGrid::kPeriod;
        grid_.skip (n);
        const bool moving = ! (core::exactlyEqual (freqSm_.value(), freqSm_.targetValue())
                               && core::exactlyEqual (airFreqSm_.value(), airFreqSm_.targetValue()));
        if (! moving) return;
        for (; boundaries > 0; --boundaries)
        {
            freqSm_.advance (core::StateGrid::kPeriod);
            airFreqSm_.advance (core::StateGrid::kPeriod);
            if (core::exactlyEqual (freqSm_.value(), freqSm_.targetValue())
                && core::exactlyEqual (airFreqSm_.value(), airFreqSm_.targetValue())) break;
        }
        designXo();
        airFreqCur_ = (float) airFreqSm_.value();
    }

    float xfTarget()  const noexcept { return (! enabled_ || lowWidth_ >= 1.0f) ? 1.0f : 0.0f; }
    float airTarget() const noexcept { return airEnabled_ ? airDb_ : 0.0f; }

    void clampAirFrequency() noexcept
    {
        // The shelf has no Nyquist clamp of its own, unlike the crossover, so the ceiling is applied
        // here and `air()` reports it. At the core's 8 kHz floor this collapses to [3000, 3600], and
        // lo <= hi still holds — the same shape applyFrequency() uses for the crossover.
        const float hi = std::max (kMinAirFreq, std::min (kMaxAirFreq, (float) (0.45 * fs_)));
        airFreq_ = std::clamp (airFreq_, kMinAirFreq, hi);
    }

    // IDLE IS NOT "gainDb == 0": a ramp still moving is still an operation. Same test StereoWidth makes
    // before it skips its own round trip, and the same reason.
    bool airIdle() const noexcept
    {
        return ! airSm_.isSmoothing() && core::exactlyEqual (airSm_.getCurrentValue(), 0.0f);
    }

    // WHILE THE RAMP MOVES the design is quantised to 0.01 dB, so a 6 dB move costs 600 redesigns
    // rather than one per sample; at REST it is designed at the exact value, so `air()` and the sound
    // agree. The quantisation is a function of the smoothed value alone, which is a function of audio
    // time alone — so it does not depend on where the caller cut the block (law 8a).
    void designAir (float db, bool moving) noexcept
    {
        const float use = moving ? std::round (db * 100.0f) * 0.01f : db;
        if (core::exactlyEqual (use, airDesignedDb_) && core::exactlyEqual (airFreqCur_, airDesignedHz_)) return;
        airDesignedDb_ = use;
        airDesignedHz_ = airFreqCur_;
        airShelf_.setCoeffs (eq::matched::highShelfDb ((double) airFreqCur_, fs_, (double) use));
    }

    double fs_ = 48000.0;
    bool   prepared_ = false;   // law 11: the crossover has no coefficients before prepare()
    bool   fresh_ = true;       // no sample since prepare()/reset(): a write SNAPS instead of gliding
    float  freq_ = 120.0f, lowWidth_ = 0.0f;
    float  xoHz_ = -1.0f;       // the corner the crossover is designed at (-1: none yet)
    bool   enabled_ = true, bypassed_ = false;
    core::Smoother freqSm_, airFreqSm_;             // the two corners' glides, on the grid
    float  airFreqCur_ = 6000.0f;                   // the air corner the shelf is designed at, as the grid walks it
    eq::Crossover2 xo_;                             // the Side-channel LR4 split (the primitive extracted from here, reused back)
    core::StateGrid grid_;                          // law 8 on audio time, never on the caller's block
    core::LinearSmoother widthSm_ { 0.0f }, xfSm_ { 0.0f };

    // K14. `airBypassed_` is the shelf's OWN latch: the crossover must not be reset because the air
    // toggled, and the shelf must not be reset because the bass settled. Two tools, two latches.
    // The measurement's own band, three LR4 pairs at the shelf's corner: Mid, Side before, Side after.
    // Reset together with their sums whenever the corner moves — the interval belongs to one band, and a
    // total accumulated across two of them is not a measurement of either.
    void retuneWidthBand() noexcept
    {
        wxM_.setFrequency (airFreq_); wxSb_.setFrequency (airFreq_); wxSa_.setFrequency (airFreq_);
        wxM_.reset(); wxSb_.reset(); wxSa_.reset();
        wSumM_ = wSumSb_ = wSumSa_ = 0.0; wSamples_ = 0;
    }

    bool   airEnabled_ = false, airBypassed_ = true;
    float  airFreq_ = 6000.0f, airDb_ = 0.0f;
    float  airDesignedDb_ = -1.0f, airDesignedHz_ = -1.0f;   // -1: nothing designed yet
    eq::Biquad airShelf_ {};
    core::LinearSmoother airSm_ { 0.0f };
    eq::Crossover2 wxM_, wxSb_, wxSa_;
    double wSumM_ = 0.0, wSumSb_ = 0.0, wSumSa_ = 0.0;
    std::int64_t wSamples_ = 0;
};

} // namespace felitronics::stereo

// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

namespace felitronics::reverb
{

//==============================================================================
// felitronics::reverb::Freeverb — Jezar's Freeverb (Dreampoint, 2000; released to the public domain):
// per channel, eight lowpass-feedback combs in parallel into four allpasses in series, the right channel's
// delays spread 23 samples past the left's. Written for this repository from the published algorithm.
//
// WHY THIS PARTICULAR FREEVERB. It reproduces, sample for sample, the behaviour of the Freeverb the shipped
// amps were voiced on (the one the JUCE framework ships as juce::Reverb), so the amps keep their sound when
// they move onto this core. That behaviour was established by MEASURING that implementation's output against
// this one's — impulses, noise, sweeps, tails and parameter moves at 44.1/48/96 kHz — not by copying it:
//   · the published tunings (combs 1116 1188 1277 1356 1422 1491 1557 1617, allpasses 556 441 341 225, spread
//     23) are scaled to the rate in INTEGER arithmetic, `size = (int) sampleRate * tuning / 44100`;
//   · the input gain is 0.015 (0 while frozen); damping is `damping × 0.4`, the comb feedback
//     `roomSize × 0.28 + 0.7` (0 and 1 while frozen); wet is `wetLevel × 3` split as
//     `wet1 = 0.5·wet·(1 + width)`, `wet2 = 0.5·wet·(1 − width)`; dry is `dryLevel × 2`; frozen means
//     `freezeMode >= 0.5`;
//   · damping, feedback, dry, wet1 and wet2 each glide LINEARLY to a new value over floor(0.01 · rate)
//     samples, one step per sample, restarting from where they stand when a new value arrives; a value equal
//     to the one already targeted (within a float epsilon) is not a new value. prepare() lands every glide on
//     its target. A fresh instance holds the published initial parameters (room 0.5, damping 0.5, wet 0.33,
//     dry 0.4, width 1) as if prepared at 44.1 kHz — so the first setParams() after the first prepare() glides
//     FROM them, as the measured implementation does;
//   · no state is quantised: the measured implementation carries every comb and allpass value at full float
//     precision on arm64 (a tail fed 1e-7 matched bit for bit), so neither does this one.
//
// LAW 8. Left at that, the combs' damping filters and every recirculating comb and allpass value decay
// geometrically into the subnormal band on silence — and stay there on any tier without hardware FTZ (wasm
// has none). So each comb's filter state and every value written into a comb or an allpass is flushed to
// exact zero below 1e-30 PER SAMPLE, inside the kernel. Per sample because the buffers are the state: a value
// written a hair above zero comes round again a comb length later, and only a flush at the write keeps every
// stored value either zero or at least 1e-30, eight decades clear of the subnormal floor — no product in the
// kernel (factors 0.5, 1 - damp >= 0.6, feedback >= 0.7) can then land in the band.
// WHY 1e-30 AND NOT THE HOUSE 1e-15 (core::flushDenormal). The flushed values are not independent: a decaying
// tail is flushed at every one of the eight combs on every sample, all with the same sign, and the loop
// integrates them over its time constant (a comb length over 1 - feedback, ~63000 samples for a long HALL).
// Measured against the reference at 1e-15: 1.2e-10 of difference on a HALL impulse — still -199 dBFS, but
// five decades above the threshold that caused it. At 1e-30 the same accumulation stays near 1e-25.
// It is the only way this Freeverb departs from the measured one.
//
// Real-time rule: prepare() allocates; process*/reset/setParams never allocate, lock or throw. No libm call
// (std::floor of the glide length in prepare() is exact).
//==============================================================================

namespace detail
{
    /** Law 8, per sample: exact zero below 1e-30 — see the header for why not 1e-15. */
    inline void flushTiny (float& x) noexcept { if (std::fabs (x) < 1.0e-30f) x = 0.0f; }

    /** Equal within one relative float epsilon (or the smallest normal, absolutely); non-finite: exact only.
        The reverb's guard for "the same value again" — in a glide target and in every AmpReverb setter. */
    inline bool nearlyEqual (float a, float b) noexcept
    {
        if (! (std::isfinite (a) && std::isfinite (b)))   // exact only: the same infinity, never a NaN
            return std::isinf (a) && std::isinf (b) && std::signbit (a) == std::signbit (b);

        const float diff = std::abs (a - b);
        return diff <= std::numeric_limits<float>::min()
            || diff <= std::numeric_limits<float>::epsilon() * std::max (std::abs (a), std::abs (b));
    }
} // namespace detail

struct FreeverbParams
{
    float roomSize   = 0.5f;    ///< 0..1
    float damping    = 0.5f;    ///< 0..1
    float wetLevel   = 0.33f;   ///< 0..1
    float dryLevel   = 0.4f;    ///< 0..1
    float width      = 1.0f;    ///< 0..1
    float freezeMode = 0.0f;    ///< >= 0.5 freezes the tail (no input, no loss)
};

class Freeverb
{
public:
    static constexpr int kNumCombs     = 8;
    static constexpr int kNumAllpasses = 4;
    static constexpr int kStereoSpread = 23;
    static constexpr int kCombTunings[kNumCombs]         = { 1116, 1188, 1277, 1356, 1422, 1491, 1557, 1617 };
    static constexpr int kAllpassTunings[kNumAllpasses]  = { 556, 441, 341, 225 };

    static constexpr float  kFixedGain     = 0.015f;
    static constexpr float  kScaleWet      = 3.0f;
    static constexpr float  kScaleDry      = 2.0f;
    static constexpr float  kScaleDamp     = 0.4f;
    static constexpr float  kScaleRoom     = 0.28f;
    static constexpr float  kOffsetRoom    = 0.7f;
    static constexpr double kGlideSeconds  = 0.01;
    static constexpr double kTuningRate    = 44100.0;

    Freeverb() noexcept
    {
        setParams (FreeverbParams {});
        landGlides (kTuningRate);
    }

    /** Sizes every comb and allpass for `sampleRate` (allocates) and clears them; every parameter glide
        lands on its target. A non-positive or non-finite rate is refused and nothing changes. */
    [[nodiscard]] bool prepare (double sampleRate)
    {
        if (! (sampleRate > 0.0) || ! std::isfinite (sampleRate) || sampleRate > 2.0e9)
            return false;

        const long long rate = (long long) sampleRate;

        for (std::size_t ch = 0; ch < 2; ++ch)
        {
            const int spread = ch == 0 ? 0 : kStereoSpread;

            for (std::size_t j = 0; j < (std::size_t) kNumCombs; ++j)
                combs[ch][j].resize (sizeFor (rate, kCombTunings[j] + spread));

            for (std::size_t j = 0; j < (std::size_t) kNumAllpasses; ++j)
                allpasses[ch][j].resize (sizeFor (rate, kAllpassTunings[j] + spread));
        }

        landGlides (sampleRate);
        prepared_ = true;
        return true;
    }

    bool isPrepared() const noexcept { return prepared_; }

    /** Clears every comb and allpass of both channels. The parameter glides are not touched. */
    void reset() noexcept
    {
        clearChannel (0);
        clearChannel (1);
    }

    /** Clears one channel's combs and allpasses, leaving the other bit-exact. */
    void clearChannel (int ch) noexcept
    {
        if (ch < 0 || ch > 1)
            return;

        for (auto& c : combs[(std::size_t) ch])     c.clear();
        for (auto& a : allpasses[(std::size_t) ch]) a.clear();
    }

    /** New targets; each value glides there over floor(0.01 · rate) samples. */
    void setParams (const FreeverbParams& p) noexcept
    {
        const bool  frozen = p.freezeMode >= 0.5f;
        const float wet    = p.wetLevel * kScaleWet;

        dry_.setTarget (p.dryLevel * kScaleDry);
        wet1_.setTarget (0.5f * wet * (1.0f + p.width));
        wet2_.setTarget (0.5f * wet * (1.0f - p.width));
        gain_ = frozen ? 0.0f : kFixedGain;

        if (frozen)
        {
            damp_.setTarget (0.0f);
            feedback_.setTarget (1.0f);
        }
        else
        {
            damp_.setTarget (p.damping * kScaleDamp);
            feedback_.setTarget (p.roomSize * kScaleRoom + kOffsetRoom);
        }
    }

    /** Both channels' combs hear `(left + right) · gain`; each side's tail is mixed back with width. */
    void processStereo (float* left, float* right, int numSamples) noexcept
    {
        if (! prepared_ || left == nullptr || right == nullptr)
            return;

        for (int i = 0; i < numSamples; ++i)
        {
            const float input = (left[i] + right[i]) * gain_;
            const float damp  = damp_.next();
            const float fb    = feedback_.next();

            float outL = 0.0f, outR = 0.0f;

            for (std::size_t j = 0; j < (std::size_t) kNumCombs; ++j)
            {
                outL += combs[0][j].process (input, damp, fb);
                outR += combs[1][j].process (input, damp, fb);
            }

            for (std::size_t j = 0; j < (std::size_t) kNumAllpasses; ++j)
            {
                outL = allpasses[0][j].process (outL);
                outR = allpasses[1][j].process (outR);
            }

            const float dry  = dry_.next();
            const float wet1 = wet1_.next();
            const float wet2 = wet2_.next();

            left[i]  = outL * wet1 + outR * wet2 + left[i]  * dry;
            right[i] = outR * wet1 + outL * wet2 + right[i] * dry;
        }
    }

    /** The left channel's combs and allpasses only; wet2 has no other side to take from. */
    void processMono (float* samples, int numSamples) noexcept
    {
        if (! prepared_ || samples == nullptr)
            return;

        for (int i = 0; i < numSamples; ++i)
        {
            const float input = samples[i] * gain_;
            const float damp  = damp_.next();
            const float fb    = feedback_.next();

            float out = 0.0f;

            for (std::size_t j = 0; j < (std::size_t) kNumCombs; ++j)
                out += combs[0][j].process (input, damp, fb);

            for (std::size_t j = 0; j < (std::size_t) kNumAllpasses; ++j)
                out = allpasses[0][j].process (out);

            const float dry  = dry_.next();
            const float wet1 = wet1_.next();   // wet2 does not move on a mono sample — as measured

            samples[i] = out * wet1 + samples[i] * dry;
        }
    }

    /** `n` samples of parameter glide with no audio — a pause in the stream, during which the glides keep
        moving exactly as they would through silence. */
    void advanceGlides (int numSamples) noexcept
    {
        for (int i = 0; i < numSamples; ++i)
        {
            (void) damp_.next(); (void) feedback_.next();
            (void) dry_.next();  (void) wet1_.next(); (void) wet2_.next();
        }
    }

    /** Longest comb of the right channel at the prepared rate, in samples (0 before prepare). */
    int longestCombSamples() const noexcept { return (int) combs[1][(std::size_t) kNumCombs - 1].buf.size(); }

private:

    static std::size_t sizeFor (long long rate, int tuning) noexcept
    {
        const long long s = rate * (long long) tuning / 44100;
        return (std::size_t) std::max (1LL, s);
    }

    /** One linear glide: `steps` per move, restarting from where it stands. */
    struct Glide
    {
        float current = 0.0f, target = 0.0f, step = 0.0f;
        int   steps = 0, countdown = 0;

        void land (int newSteps) noexcept { steps = newSteps; current = target; countdown = 0; }

        void setTarget (float v) noexcept
        {
            if (detail::nearlyEqual (v, target))
                return;

            if (steps <= 0)
            {
                current = target = v;
                countdown = 0;
                return;
            }

            target    = v;
            countdown = steps;
            step      = (target - current) / (float) countdown;
        }

        float next() noexcept
        {
            if (countdown <= 0)
                return target;

            --countdown;
            if (countdown > 0)
                current += step;
            else
                current = target;

            return current;
        }

    };

    struct Comb
    {
        std::vector<float> buf;
        std::size_t idx  = 0;
        float       last = 0.0f;

        void resize (std::size_t n)
        {
            if (n != buf.size())
            {
                buf.assign (n, 0.0f);
                idx = 0;
            }
            clear();
        }

        void clear() noexcept
        {
            std::fill (buf.begin(), buf.end(), 0.0f);
            last = 0.0f;
        }

        float process (float input, float damp, float feedback) noexcept
        {
            const float out = buf[idx];
            last = out * (1.0f - damp) + last * damp;
            detail::flushTiny (last);
            float w = input + last * feedback;
            detail::flushTiny (w);
            buf[idx] = w;
            if (++idx == buf.size())
                idx = 0;
            return out;
        }
    };

    struct Allpass
    {
        std::vector<float> buf;
        std::size_t idx = 0;

        void resize (std::size_t n)
        {
            if (n != buf.size())
            {
                buf.assign (n, 0.0f);
                idx = 0;
            }
            clear();
        }

        void clear() noexcept { std::fill (buf.begin(), buf.end(), 0.0f); }

        float process (float input) noexcept
        {
            const float held = buf[idx];
            float w = input + held * 0.5f;
            detail::flushTiny (w);
            buf[idx] = w;
            if (++idx == buf.size())
                idx = 0;
            return held - input;
        }
    };

    void landGlides (double sampleRate) noexcept
    {
        const int steps = (int) std::floor (kGlideSeconds * sampleRate);
        damp_.land (steps); feedback_.land (steps);
        dry_.land (steps);  wet1_.land (steps); wet2_.land (steps);
    }

    std::array<std::array<Comb, kNumCombs>, 2>        combs;
    std::array<std::array<Allpass, kNumAllpasses>, 2> allpasses;
    Glide damp_, feedback_, dry_, wet1_, wet2_;
    float gain_     = kFixedGain;
    bool  prepared_ = false;
};

} // namespace felitronics::reverb

// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/core/Config.h>
#include <felitronics/core/Math.h>
#include <felitronics/core/DelayLine.h>
#include <felitronics/core/FlushToZero.h>
#include <felitronics/dynamics/GainComputer.h>
#include <felitronics/dynamics/GainReductionFollower.h>
#include <felitronics/dynamics/ChannelLinker.h>
#include <felitronics/dynamics/EnvelopeFollower.h>   // Detector enum (Peak / Rms)

#include <cmath>
#include <vector>

namespace felitronics::dynamics
{

struct CompressorParams
{
    Detector detector    = Detector::Peak;       // Peak = instant |x|; Rms = short-window power
    LinkMode link        = LinkMode::Max;         // how channels collapse to one gain
    Mode     mode        = Mode::DownCompress;

    double thresholdDb = -18.0;
    double ratio       = 2.0;
    double kneeDb      = 6.0;
    double rangeDb     = 60.0;

    double attackMs    = 10.0;                    // ballistics on the GAIN REDUCTION (not the level)
    double releaseMs   = 100.0;
    double rmsWindowMs = 5.0;                     // only used when detector == Rms

    double makeupDb    = 0.0;
    bool   autoMakeup  = false;                   // adds a STATIC auto-makeup (see autoMakeupDb)
    double lookaheadMs = 0.0;
};

// A STATIC auto-makeup (dB): half-compensate the gain reduction a 0 dBFS signal would receive. Static
// (a function of the curve only) — NOT signal-following (which would be a second compressor and pump).
inline double autoMakeupDb (const GainComputer& gc) noexcept { return -0.5 * gc.deltaDb (0.0); }

//==============================================================================
// felitronics::dynamics::Compressor — a broadband, channel-LINKED compressor that composes the dynamics
// primitives in the clean topology:
//     instant linked detector → stateless GainComputer (curve) → GR ballistics (timing) → makeup
//     → apply to the LOOKAHEAD-delayed signal.
// Smoothing the gain reduction (not the detector level) keeps the attack/release independent of the
// knee/ratio. Product-neutral: no sidechain EQ, params system, GUI, metering policy, dry/wet — those
// stay in the product (it's the "core gives primitives + a broadband composite; products compose" line).
// RT-safe: prepare() allocates; process() does no alloc/lock/throw. Reported latency = lookahead.
class Compressor
{
public:
    void prepare (double sampleRate, int /*maxBlock*/, int maxChannels, double maxLookaheadMs = 50.0)
    {
        fs = sampleRate;
        const int mc = maxChannels < 1 ? 1 : (maxChannels > core::kMaxChannels ? core::kMaxChannels : maxChannels);
        maxLookSamples = (int) std::ceil (maxLookaheadMs * 0.001 * fs);
        delays.assign ((std::size_t) mc, core::DelayLine {});
        for (auto& d : delays) d.prepare (maxLookSamples);
        grFollower.prepare (fs);
        apply (params);
        reset();
    }

    void reset() noexcept
    {
        for (auto& d : delays) d.reset();
        grFollower.reset();
        rmsState = 0.0f;
    }

    void   setParams (const CompressorParams& p) noexcept { params = p; apply (p); }
    int    latencySamples()  const noexcept { return lookSamples; }
    double gainReductionDb() const noexcept { return grFollower.valueDb(); }   // for metering

    // Audio thread, in place. RT-safe.
    void process (float* const* channels, int numChannels, int numSamples) noexcept
    {
        const int nc = numChannels < (int) delays.size() ? numChannels : (int) delays.size();
        if (nc <= 0) return;

        for (int i = 0; i < numSamples; ++i)
        {
            // --- detector: one linked level (+ optional RMS time-averaging) ---
            const float linked = linkAmplitude (params.link, channels, nc, i);
            float level;
            if (params.detector == Detector::Rms)
            {
                const float pw = linked * linked;
                rmsState = pw + rmsCoeff * (rmsState - pw);
                level = std::sqrt (rmsState);
            }
            else
            {
                level = linked;
            }

            // --- curve → GR ballistics → gain ---
            const double levelDb    = core::gainToDb (level);
            const float  targetDb   = (float) gc.deltaDb (levelDb);
            const float  grDb       = grFollower.process (targetDb);
            const float  gain       = (float) core::dbToGain ((double) grDb + makeupAppliedDb);

            // --- apply to the lookahead-delayed signal (read input before overwriting in place) ---
            for (int c = 0; c < nc; ++c)
            {
                const float x = channels[c][i];
                const float y = delays[(std::size_t) c].process (x);
                channels[c][i] = y * gain;
            }
        }

        grFollower.flushDenormals();
        core::flushDenormal (rmsState);
    }

private:
    void apply (const CompressorParams& p) noexcept
    {
        gc.setMode (p.mode);
        gc.setThresholdDb (p.thresholdDb);
        gc.setRatio (p.ratio);
        gc.setKneeDb (p.kneeDb);
        gc.setRangeDb (p.rangeDb);
        grFollower.setTimes (p.attackMs, p.releaseMs);   // non-finite times → instant (the follower's coeff guard)

        // Non-finite params fall back (house rule): a NaN window/lookahead/makeup would leave a NaN coeff,
        // feed lround() UB, or poison the per-sample gain.
        const double t = (std::isfinite (p.rmsWindowMs) ? p.rmsWindowMs : 0.0) * 0.001;
        rmsCoeff = (t <= 0.0 || fs <= 0.0) ? 0.0f : (float) std::exp (-1.0 / (t * fs));

        lookSamples = (int) std::lround ((std::isfinite (p.lookaheadMs) ? p.lookaheadMs : 0.0) * 0.001 * fs);
        if (lookSamples < 0) lookSamples = 0;
        if (lookSamples > maxLookSamples) lookSamples = maxLookSamples;
        for (auto& d : delays) d.setDelay (lookSamples);

        makeupAppliedDb = (std::isfinite (p.makeupDb) ? p.makeupDb : 0.0) + (p.autoMakeup ? autoMakeupDb (gc) : 0.0);
    }

    double fs = 48000.0;
    CompressorParams params;

    GainComputer          gc;
    GainReductionFollower  grFollower;
    std::vector<core::DelayLine> delays;

    float  rmsState = 0.0f, rmsCoeff = 0.0f;
    int    lookSamples = 0, maxLookSamples = 0;
    double makeupAppliedDb = 0.0;
};

} // namespace felitronics::dynamics

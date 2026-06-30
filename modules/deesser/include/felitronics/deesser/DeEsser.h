// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/eq/Crossover2.h>
#include <felitronics/eq/Svf.h>
#include <felitronics/dynamics/EnvelopeFollower.h>
#include <felitronics/dynamics/ChannelLinker.h>
#include <felitronics/dynamics/GainComputer.h>
#include <felitronics/dynamics/GainReductionFollower.h>
#include <felitronics/dynamiceq/DynamicEqBand.h>
#include <felitronics/core/Math.h>
#include <felitronics/core/Config.h>

#include <algorithm>
#include <cmath>

namespace felitronics::deesser
{

//==============================================================================
// felitronics::deesser::DeEsser — tame sibilance. Two topologies:
//   DynamicEq  (default) — a CutWhenLoud Bell at fc (delegates to dynamiceq::DynamicEqBand). SURGICAL: only
//              the sibilant band moves, the rest of the top stays — mastering-safe.
//   SplitBand            — an LR4 eq::Crossover2 at fc; the high band is DUCKED (the low band untouched);
//              the allpass-flat sum makes it transparent when idle. The classic clean de-esser,
//              but ducking the WHOLE band can dull air on a full mix.
//
// Detection is a BANDPASS sidechain at fc (NOT a high-pass) so loud non-sibilant HF — cymbals, air — doesn't
// false-trigger; `rangeDb` caps the cut to avoid lisping. One linked detector → one gain on both channels.
// `listen` solos the sidechain so you can hear what's being detected. RT-safe; zero latency.
enum class DeEsserMode { DynamicEq, SplitBand };

struct DeEsserParams
{
    DeEsserMode mode = DeEsserMode::DynamicEq;
    double fc = 7000.0, scQ = 2.0;                 // sidechain band + crossover/EQ centre
    double thresholdDb = -30.0, ratio = 3.0, kneeDb = 6.0, rangeDb = 8.0;
    double attackMs = 2.0, releaseMs = 90.0;
    dynamics::Detector detector = dynamics::Detector::Rms;
    dynamics::LinkMode link     = dynamics::LinkMode::Max;
    bool listen = false;                           // solo the sidechain (hear the sibilance the detector sees)
};

class DeEsser
{
public:
    void prepare (double sampleRate, int maxBlock, int maxChannels) noexcept
    {
        fs_ = sampleRate > 0.0 ? sampleRate : 48000.0;
        channels_ = std::clamp (maxChannels, 1, core::kMaxChannels);
        side_.prepare (fs_, channels_);
        xover_.prepare (fs_, channels_);
        env_.prepare (fs_); gr_.prepare (fs_);
        (void) maxBlock;
        deq_.prepare (fs_, channels_);
        apply (params_);
        reset();
    }

    void reset() noexcept { side_.reset(); xover_.reset(); env_.reset(); gr_.reset(); deq_.reset(); grDb_ = 0.0f; }
    void setParams (const DeEsserParams& p) noexcept
    {
        const bool relatch = (p.mode != params_.mode) || (p.listen != params_.listen);   // topology / signal-path change
        params_ = p; apply (p);
        if (relatch) reset();   // clear the now-inactive path's stale state → no thump when toggling mode/listen live
    }
    static constexpr int latencySamples() noexcept { return 0; }
    double gainReductionDb() const noexcept { return grDb_; }

    void process (float* const* io, int numChannels, int n) noexcept
    {
        const int nc = std::min (numChannels, channels_);
        if (nc <= 0) return;

        // DynamicEq mode = the surgical dynamic-EQ band (no split, no listen detour).
        if (params_.mode == DeEsserMode::DynamicEq && ! params_.listen)
        {
            deq_.process (io, numChannels, n);
            grDb_ = (float) deq_.dynamicDeltaDb();
            return;
        }

        for (int i = 0; i < n; ++i)
        {
            // detector: per-channel sidechain BandPass → one linked level
            float sc[core::kMaxChannels]; float linked = 0.0f, sq = 0.0f;
            for (int c = 0; c < nc; ++c)
            {
                sc[c] = side_.processSample (c, io[c][i]);
                const float a = std::fabs (sc[c]);
                if (link_ == dynamics::LinkMode::Max) { if (a > linked) linked = a; } else sq += sc[c] * sc[c];
            }
            if (link_ == dynamics::LinkMode::MeanPower) linked = std::sqrt (sq / (float) nc);

            const float lvl = env_.process (linked);
            const float gr  = gr_.process ((float) gc_.deltaDb ((double) core::gainToDb ((double) std::max (lvl, 1.0e-9f))));   // ≤ 0 (cut)
            grDb_ = gr;
            const float gain = (float) core::dbToGain ((double) gr);

            for (int c = 0; c < nc; ++c)
            {
                if (params_.listen) { io[c][i] = sc[c]; continue; }
                float lo, hi; xover_.processSample (c, io[c][i], lo, hi);
                io[c][i] = lo + gain * hi;                          // duck only the high band
            }
        }
        side_.flushDenormals(); xover_.flushDenormals(); env_.flushDenormals(); gr_.flushDenormals();
    }

private:
    void apply (const DeEsserParams& p) noexcept
    {
        side_.setParams (eq::FilterType::BandPass, p.fc, p.scQ, 0.0);
        xover_.setFrequency ((float) p.fc);
        env_.setDetector (p.detector); env_.setTimes (2.0, 2.0);
        gc_.setMode (dynamics::Mode::DownCompress);
        gc_.setThresholdDb (p.thresholdDb); gc_.setRatio (p.ratio); gc_.setKneeDb (p.kneeDb); gc_.setRangeDb (p.rangeDb);
        gr_.setTimes (p.attackMs, p.releaseMs);
        link_ = p.link;

        dynamiceq::DynamicEqBandParams dp;
        dp.type = eq::FilterType::Bell; dp.freq = p.fc; dp.Q = p.scQ; dp.mode = dynamiceq::DynamicEqMode::CutWhenLoud;
        dp.thresholdDb = p.thresholdDb; dp.ratio = p.ratio; dp.kneeDb = p.kneeDb; dp.rangeDb = p.rangeDb;
        dp.attackMs = p.attackMs; dp.releaseMs = p.releaseMs; dp.detector = p.detector; dp.link = p.link;
        deq_.setParams (dp);
    }

    double fs_ = 48000.0;
    int channels_ = 2;
    DeEsserParams params_;

    eq::Svf side_;                               // sidechain BandPass
    eq::Crossover2 xover_;                       // LR4 split (SplitBand mode)
    dynamics::EnvelopeFollower env_;
    dynamics::GainComputer gc_;
    dynamics::GainReductionFollower gr_;
    dynamiceq::DynamicEqBand deq_;               // DynamicEq mode delegate

    dynamics::LinkMode link_ = dynamics::LinkMode::Max;
    float grDb_ = 0.0f;
};

} // namespace felitronics::deesser

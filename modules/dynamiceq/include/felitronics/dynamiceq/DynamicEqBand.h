// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/eq/Svf.h>
#include <felitronics/dynamics/EnvelopeFollower.h>
#include <felitronics/dynamics/ChannelLinker.h>
#include <felitronics/dynamics/GainComputer.h>
#include <felitronics/dynamics/GainReductionFollower.h>
#include <felitronics/core/Math.h>
#include <felitronics/core/Config.h>

#include <algorithm>
#include <cmath>

namespace felitronics::dynamiceq
{

//==============================================================================
// felitronics::dynamiceq::DynamicEqBand — an EQ band whose gain moves with a level detector: tame a
// resonance only when it rings, de-ess, dynamic shelf. The cornerstone "surgical-but-musical" mastering tool.
//
// Detector: a SIDECHAIN eq::Svf BandPass at the band's freq/Q isolates the band's energy (independent of the
// band's own moving gain) → EnvelopeFollower level → dynamics::GainComputer curve → dynamics::GainReduction-
// Follower ballistics (smooth the GAIN, not the level — the Compressor's lesson). One linked detector → one
// gain on all channels (image-preserving). Audio: an eq::Svf Bell/shelf whose gain = staticGainDb + the
// dynamic delta. The Cytomic gain enters the damping (k = 1/(Q·A)) so a gain change needs a full coeff
// recompute → do it at CONTROL RATE (every `coeffUpdatePeriod` samples, K=16 ≈ 0.33 ms) off the smoothed
// gain, keeping the SAME bell shape as the static EQ (a Regalia–Mitra constant-Q form was rejected: cheaper
// but an asymmetric bell that wouldn't match eq::Svf::Bell).
//
//   Mode             GainComputer       sign   result
//   CutWhenLoud      DownCompress        +1    loud → cut   (de-ess / tame resonance)
//   BoostWhenLoud    DownCompress        −1    loud → boost
//   BoostWhenQuiet   UpCompress          +1    quiet → boost
//
// RT-safe: prepare() allocates; process() is alloc/lock/throw-free, in place. Zero latency.
enum class DynamicEqMode { CutWhenLoud, BoostWhenLoud, BoostWhenQuiet };

struct DynamicEqBandParams
{
    eq::FilterType type = eq::FilterType::Bell;
    double freq = 1000.0, Q = 2.0, staticGainDb = 0.0;
    double thresholdDb = -24.0, ratio = 2.0, kneeDb = 6.0, rangeDb = 12.0;
    double attackMs = 10.0, releaseMs = 120.0;
    DynamicEqMode    mode     = DynamicEqMode::CutWhenLoud;
    dynamics::Detector detector = dynamics::Detector::Rms;
    dynamics::LinkMode link     = dynamics::LinkMode::Max;
    int coeffUpdatePeriod = 16;                 // control-rate gain recompute (samples)
};

class DynamicEqBand
{
public:
    void prepare (double sampleRate, int maxChannels) noexcept
    {
        fs_ = sampleRate > 0.0 ? sampleRate : 48000.0;
        channels_ = std::clamp (maxChannels, 1, core::kMaxChannels);
        side_.prepare (fs_, channels_);
        audio_.prepare (fs_, channels_);
        env_.prepare (fs_);
        gr_.prepare (fs_);
        apply (params_);
        reset();
    }

    void reset() noexcept
    {
        side_.reset(); audio_.reset(); env_.reset(); gr_.reset();
        ksamp_ = 0; curGainDb_ = params_.staticGainDb; updateAudio();
    }

    void setParams (const DynamicEqBandParams& p) noexcept { params_ = p; apply (p); }
    static constexpr int latencySamples() noexcept { return 0; }
    double dynamicDeltaDb() const noexcept { return curGainDb_ - params_.staticGainDb; }   // for metering

    void process (float* const* io, int numChannels, int n) noexcept
    {
        const int nc = std::min (numChannels, channels_);
        if (nc <= 0) return;
        for (int i = 0; i < n; ++i)
        {
            // detector: per-channel sidechain BandPass → one linked level
            float linked = 0.0f;
            if (link_ == dynamics::LinkMode::Max)
            {
                for (int c = 0; c < nc; ++c) { const float a = std::fabs (side_.processSample (c, io[c][i])); if (a > linked) linked = a; }
            }
            else
            {
                float sq = 0.0f; for (int c = 0; c < nc; ++c) { const float bp = side_.processSample (c, io[c][i]); sq += bp * bp; }
                linked = std::sqrt (sq / (float) nc);
            }
            const float lvl   = env_.process (linked);
            const float lvlDb = (float) core::gainToDb ((double) std::max (lvl, 1.0e-9f));
            const float delta = sign_ * (float) gc_.deltaDb ((double) lvlDb);
            const float smooth = gr_.process (delta);

            if (ksamp_ == 0) { curGainDb_ = (double) params_.staticGainDb + (double) smooth; updateAudio(); }
            if (++ksamp_ >= K_) ksamp_ = 0;

            for (int c = 0; c < nc; ++c) io[c][i] = audio_.processSample (c, io[c][i]);
        }
        env_.flushDenormals(); gr_.flushDenormals(); side_.flushDenormals(); audio_.flushDenormals();
    }

private:
    void updateAudio() noexcept { audio_.setParams (params_.type, params_.freq, params_.Q, curGainDb_); }

    void apply (const DynamicEqBandParams& p) noexcept
    {
        side_.setParams (eq::FilterType::BandPass, p.freq, p.Q, 0.0);
        env_.setDetector (p.detector);
        env_.setTimes (2.0, 2.0);                                 // a short, symmetric level window; the GR follower does the musical timing
        gc_.setThresholdDb (p.thresholdDb); gc_.setRatio (p.ratio); gc_.setKneeDb (p.kneeDb); gc_.setRangeDb (p.rangeDb);
        gr_.setTimes (p.attackMs, p.releaseMs);
        link_ = p.link;
        K_ = std::max (1, p.coeffUpdatePeriod);
        switch (p.mode)
        {
            case DynamicEqMode::CutWhenLoud:    gc_.setMode (dynamics::Mode::DownCompress); sign_ =  1.0f; break;
            case DynamicEqMode::BoostWhenLoud:  gc_.setMode (dynamics::Mode::DownCompress); sign_ = -1.0f; break;
            case DynamicEqMode::BoostWhenQuiet: gc_.setMode (dynamics::Mode::UpCompress);   sign_ =  1.0f; break;
        }
        updateAudio();
    }

    double fs_ = 48000.0;
    int channels_ = 2;
    DynamicEqBandParams params_;

    eq::Svf side_, audio_;                       // sidechain BandPass · audio Bell/shelf
    dynamics::EnvelopeFollower    env_;
    dynamics::GainComputer        gc_;
    dynamics::GainReductionFollower gr_;

    dynamics::LinkMode link_ = dynamics::LinkMode::Max;
    float  sign_ = 1.0f;
    double curGainDb_ = 0.0;
    int    ksamp_ = 0, K_ = 16;
};

} // namespace felitronics::dynamiceq

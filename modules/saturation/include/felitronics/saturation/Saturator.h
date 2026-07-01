// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/saturation/WaveShaper.h>
#include <felitronics/oversampling/PolyphaseOversampler.h>
#include <felitronics/core/Config.h>
#include <felitronics/core/Math.h>

#include <algorithm>
#include <cmath>
#include <vector>

namespace felitronics::saturation
{

//==============================================================================
// felitronics::saturation::Saturator — the production soft-saturation stage (the mastering "glue/color"
// before the make-loud gain). Wraps a stateless WaveShaper curve in OVERSAMPLING (default 4×) so the new
// harmonics fold above the base Nyquist instead of aliasing back, runs a DC blocker inside the oversampled
// region (the asymmetric curve's even harmonics shift the mean), then returns to base rate and applies
// DRIVE-COMPENSATION + a linear dry/wet + output trim.
//
// Gain-staging (a reference tool reverted its saturator twice over this): the curve
// is peak-normalised (|x|≤1 → |y|≤1), `autoComp` undoes the small-signal-gain bump so loudness doesn't jump
// with drive, and the dry/wet blend is LINEAR (convex combination of bounded signals → peak-safe). Place it
// BEFORE the loudness gain + final limiter.
//
// RT-safe: prepare() does all allocation; process() is alloc/lock/throw-free, in place. Params apply per
// block (mastering params are near-static; a parameter smoother can wrap setParams() for automation). n
// must be ≤ the maxBlock passed to prepare(). With oversampling the stage reports a round-trip latency.
class Saturator
{
public:
    struct Params
    {
        WaveShaper::Shape shape = WaveShaper::Shape::Tanh;
        float driveDb   = 3.0f;    // 0..36 (mastering 1..6); k = dbToGain(driveDb) − 1
        float bias      = 0.0f;    // −0.5..0.5 — Asym only, sets the even-harmonic content
        float mix       = 1.0f;    // 0..1 dry/wet
        float outputDb  = 0.0f;    // −24..24 post trim
        float autoComp  = 0.5f;    // 0..1 — drive-compensation amount (1 = unity small-signal gain)
        float dcBlockHz = 10.0f;   // DC blocker corner (in the oversampled domain)
    };

    bool prepare (double sampleRate, int maxBlock, int maxChannels, int oversampleFactor = 4, int tapsPerPhase = 32)
    {
        prepared_ = false;                                             // any early return below leaves it unprepared
        if (sampleRate <= 0.0 || maxBlock < 1) return false;
        fs_       = sampleRate;
        maxBlock_ = maxBlock;
        channels_ = std::clamp (maxChannels, 1, core::kMaxChannels);
        os_       = (oversampleFactor >= 2) ? oversampleFactor : 1;
        if (os_ > 1 && ! ovs_.prepare (os_, channels_, tapsPerPhase)) return false;

        osBuf_.assign  ((std::size_t) channels_ * (std::size_t) (maxBlock_ * os_), 0.0f);
        wetBuf_.assign ((std::size_t) channels_ * (std::size_t) maxBlock_,         0.0f);
        osPtrs_.assign ((std::size_t) channels_, nullptr);
        wetPtrs_.assign((std::size_t) channels_, nullptr);
        dcX1_.assign   ((std::size_t) channels_, 0.0f);
        dcY1_.assign   ((std::size_t) channels_, 0.0f);
        applyParams();
        prepared_ = true;                                              // fully built — process() may now run
        return true;
    }

    void reset() noexcept
    {
        if (os_ > 1) ovs_.reset();
        std::fill (dcX1_.begin(), dcX1_.end(), 0.0f);
        std::fill (dcY1_.begin(), dcY1_.end(), 0.0f);
    }

    int  latencySamples() const noexcept { return os_ > 1 ? ovs_.latencySamples() : 0; }

    void setParams (const Params& p) noexcept { params_ = p; applyParams(); }

    // In place, planar. RT-safe. n must be ≤ maxBlock.
    void process (float* const* io, int numChannels, int n) noexcept
    {
        const int nc = std::min (numChannels, channels_);
        if (! prepared_ || nc <= 0 || n <= 0 || n > maxBlock_) return;   // unprepared / failed-prepare / oversized block → no OOB
        const int osN = n * os_;
        for (int c = 0; c < nc; ++c)
        {
            osPtrs_[(std::size_t) c]  = &osBuf_[(std::size_t) c * (std::size_t) (maxBlock_ * os_)];
            wetPtrs_[(std::size_t) c] = &wetBuf_[(std::size_t) c * (std::size_t) maxBlock_];
        }

        // 1) to the oversampled domain (io is only READ here, so it still holds the dry signal below)
        if (os_ > 1) ovs_.upsample (io, nc, n, osPtrs_.data());
        else for (int c = 0; c < nc; ++c) std::copy (io[c], io[c] + n, osPtrs_[(std::size_t) c]);

        // 2) waveshape per channel in the oversampled domain; DC-block only the asymmetric curve (the
        //    even-harmonic offset). Symmetric curves are zero-mean → no blocker → no needless phase shift.
        for (int c = 0; c < nc; ++c)
        {
            float* b = osPtrs_[(std::size_t) c];
            if (dcEnabled_)
            {
                float x1 = dcX1_[(std::size_t) c], y1 = dcY1_[(std::size_t) c];
                for (int i = 0; i < osN; ++i)
                {
                    const float w  = shaper_.processSample (b[i]);
                    const float dc = w - x1 + dcR_ * y1;     // one-pole DC blocker  H(z)=(1-z⁻¹)/(1-R·z⁻¹)
                    x1 = w; y1 = dc;
                    b[i] = dc;
                }
                dcX1_[(std::size_t) c] = x1; dcY1_[(std::size_t) c] = y1;
            }
            else
            {
                for (int i = 0; i < osN; ++i) b[i] = shaper_.processSample (b[i]);
            }
        }

        // 3) back to base rate (wet)
        if (os_ > 1) ovs_.downsample (osPtrs_.data(), nc, n, wetPtrs_.data());
        else for (int c = 0; c < nc; ++c) std::copy (osPtrs_[(std::size_t) c], osPtrs_[(std::size_t) c] + n, wetPtrs_[(std::size_t) c]);

        // 4) drive-compensate the wet, then a linear (peak-safe) dry/wet, then output trim
        for (int c = 0; c < nc; ++c)
            for (int i = 0; i < n; ++i)
            {
                const float dry = io[c][i];
                const float wet = comp_ * wetPtrs_[(std::size_t) c][i];
                io[c][i] = outGain_ * ((1.0f - mix_) * dry + mix_ * wet);
            }
    }

private:
    void applyParams() noexcept
    {
        shaper_.setShape (params_.shape);
        shaper_.setBias  (params_.bias);
        shaper_.setDrive (core::dbToGain (params_.driveDb) - 1.0f);     // driveDb 0 → k≈0 (linear)
        comp_    = std::pow (std::max (1.0e-6f, shaper_.slopeAtZero()), -std::clamp (params_.autoComp, 0.0f, 1.0f));
        mix_     = std::clamp (params_.mix, 0.0f, 1.0f);
        outGain_ = core::dbToGain (params_.outputDb);
        const double fsOs = fs_ * (double) os_;
        const double fc   = std::clamp ((double) params_.dcBlockHz, 0.0, 0.49 * fsOs);
        dcR_ = (fc <= 0.0) ? 0.0f : (float) std::exp (-2.0 * core::kPi * fc / fsOs);
        dcEnabled_ = (params_.dcBlockHz > 0.0f) && (params_.shape == WaveShaper::Shape::Asym);
    }

    Params      params_ {};
    WaveShaper  shaper_ {};
    oversampling::PolyphaseOversampler ovs_;

    double fs_ = 0.0;
    int    maxBlock_ = 0, channels_ = 0, os_ = 1;
    float  comp_ = 1.0f, dcR_ = 0.0f, mix_ = 1.0f, outGain_ = 1.0f;
    bool   dcEnabled_ = false;
    bool   prepared_  = false;                             // true only after a fully-successful prepare()

    std::vector<float>  osBuf_, wetBuf_;
    std::vector<float*> osPtrs_, wetPtrs_;
    std::vector<float>  dcX1_, dcY1_;
};

} // namespace felitronics::saturation

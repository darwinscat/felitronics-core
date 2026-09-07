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
    [[nodiscard]] bool prepare (double sampleRate, int maxChannels) noexcept
    {
        prepared_ = false;                                        // any early return below leaves it unprepared
        if (maxChannels < 1 || maxChannels > core::kMaxChannels) return false;   // law 11(b): BINDING
        fs_ = sampleRate > 0.0 ? sampleRate : 48000.0;
        channels_ = maxChannels;
        side_.prepare (fs_, channels_);
        audio_.prepare (fs_, channels_);
        env_.prepare (fs_);
        gr_.prepare (fs_);
        apply (params_);
        reset();
        prepared_ = true;
        return true;
    }

    bool isPrepared() const noexcept { return prepared_; }

    void reset() noexcept
    {
        side_.reset(); audio_.reset(); env_.reset(); gr_.reset();
        ksamp_ = 0; curGainDb_ = params_.staticGainDb; updateAudio();
        ranNc_ = 0;   // nothing has run, so nothing can be stopping (see dropStoppedChannels)
    }

    void setParams (const DynamicEqBandParams& p) noexcept { params_ = p; apply (p); }
    static constexpr int latencySamples() noexcept { return 0; }
    double dynamicDeltaDb() const noexcept { return curGainDb_ - params_.staticGainDb; }   // for metering

    // Law 11 (DSP-ARCHITECTURE.md §2).
    [[nodiscard]] bool process (float* const* io, int numChannels, int n) noexcept
    {
        if (numChannels < 0 || n < 0) return false;
        if (! prepared_) return false;
        if (numChannels > channels_) return false;                // width is a LIMIT — law 11(b)
        if (n == 0) return true;                                  // no samples: no time, no edge
        const int nc = numChannels;
        dropStoppedChannels (nc);                                 // law 11(d): the edge is clocked by n
        if (nc == 0) { advanceSilence (n); return true; }         // law 11c: a pause IS silence
        for (int i = 0; i < n; ++i)
        {
            // detector: per-channel sidechain BandPass → one linked level
            float linked = 0.0f;
            if (link_ == dynamics::LinkMode::Max)
            {
                // Gated on the filter's INPUT — see the note in DeEsser: the sidechain filter is
                // recursive, so gating its output would only change the sign of the defect.
                for (int c = 0; c < nc; ++c) { const float a = std::fabs (side_.processSample (c, dynamics::detectorGate (io[c][i]))); if (a > linked) linked = a; }
            }
            else
            {
                float sq = 0.0f; for (int c = 0; c < nc; ++c) { const float bp = side_.processSample (c, dynamics::detectorGate (io[c][i])); sq += bp * bp; }
                linked = std::sqrt (sq / (float) nc);
            }
            sharedStep (linked);

            for (int c = 0; c < nc; ++c) io[c][i] = audio_.processSample (c, io[c][i]);
        }
        env_.flushDenormals(); gr_.flushDenormals(); side_.flushDenormals(); audio_.flushDenormals();
        return true;
    }

private:
    // LAW 11c — `n` samples of DIGITAL SILENCE through the SHARED half of the loop above. The per-channel
    // sidechain columns were just dropped by `dropStoppedChannels(0)`, so a silent band's linked level is
    // exactly +0.0f: the loop's own MeanPower branch cannot be reused for it — `std::sqrt (sq / (float) nc)`
    // is 0/0 at width zero — so the level is PINNED rather than computed, which is the same number by a
    // route that has no division in it. Freezing instead held -21.774 dB of dynamic delta through a
    // second of gap where silence releases to -2.985 (and to -0.000 through ten seconds).
    //
    // THE CONTROL COUNTER IS PART OF THE CLOCK. `ksamp_` selects which samples redesign the audio filter,
    // and a pause that did not advance it would put the band's coefficient grid out of phase with the
    // stream for ever after. It runs honestly in the loop below; once the CONTINUOUS state has stopped
    // moving, the only thing left to do is the counter, and that closes in the obvious form — but only
    // AFTER at least one honest step, because `apply()` can lower `K_` without normalising `ksamp_`, and
    // `(ksamp_ + rem) % K_` disagrees with iterating whenever `ksamp_ >= K_` (from ksamp_ = 10, K_ = 4 the
    // loop reaches 0 in one step; the formula says 3). One honest step maps any counter into [0, K_), and
    // the loop below has always run one by the time the shortcut is reached. The arithmetic is widened
    // before the add for the reason `core::StateGrid::skip` widens its own: `ksamp_ + rem` is an int.
    void advanceSilence (int n) noexcept
    {
        int i = 0;
        for (; i < n; ++i)
        {
            const float e0 = env_.stateWord(), g0 = gr_.valueDb();
            const double c0 = curGainDb_;
            sharedStep (0.0f);   // a silent band links to exactly +0.0f at every width
            if (core::sameBits (e0, env_.stateWord()) && core::sameBits (g0, gr_.valueDb())
                && core::exactlyEqual (c0, curGainDb_))
            { ++i; break; }
            if (! std::isfinite (env_.stateWord()) && ! std::isfinite (e0)
                && ! std::isfinite (gr_.valueDb()) && ! std::isfinite (g0)) { ++i; break; }
        }
        const long long rem = (long long) n - (long long) i;
        if (rem > 0)
        {
            // A coefficient write lands iff some remaining index enters the body at counter 0. The first
            // such index is 1 when ksamp_ == 0 and K_ - ksamp_ + 1 otherwise; it is idempotent here
            // because `smooth` no longer moves, so applying it once is applying it every time.
            if (ksamp_ == 0 || rem > (long long) K_ - (long long) ksamp_) { curGainDb_ = (double) params_.staticGainDb + (double) gr_.valueDb(); updateAudio(); }
            ksamp_ = (int) (((long long) ksamp_ + rem) % (long long) K_);
        }
        env_.flushDenormals(); gr_.flushDenormals(); side_.flushDenormals(); audio_.flushDenormals();
    }

    // ONE sample of the SHARED body — everything downstream of the link and upstream of the per-channel
    // audio filter. It is written once and called from both loops, so the pause and the audio cannot end
    // up with two spellings of the same arithmetic that merely agree today; the narrowing of `lvlDb` to
    // float BEFORE the curve is one of the differences between this stage, `Compressor` and `DeEsser`
    // that a shared "silence primitive" would have quietly flattened.
    void sharedStep (float linked) noexcept
    {
        const float lvl   = env_.process (linked);
        const float lvlDb = (float) core::gainToDb ((double) std::max (lvl, 1.0e-9f));
        const float delta = sign_ * (float) gc_.deltaDb ((double) lvlDb);
        const float smooth = gr_.process (delta);
        if (ksamp_ == 0) { curGainDb_ = (double) params_.staticGainDb + (double) smooth; updateAudio(); }
        if (++ksamp_ >= K_) ksamp_ = 0;
    }

    // Both filters are per channel and advance only for c < nc; the envelope and the GR follower they feed
    // are shared and keep running. A channel that leaves and RETURNS therefore re-enters with a sidechain
    // column and an audio column frozen from before the gap — measured, on DIGITAL SILENCE, an 0.388 audio
    // tail eleven samples in and 8.85 dB of unearned reduction. Per channel only: the shared detector path
    // legitimately follows whatever channels are actually present.
    void dropStoppedChannels (int nc) noexcept
    {
        for (int c = nc; c < ranNc_; ++c) { side_.resetChannel (c); audio_.resetChannel (c); }
        ranNc_ = nc;
    }

    void updateAudio() noexcept { audio_.setParams (params_.type, params_.freq, params_.Q, curGainDb_); }

    static double finite (double v, double fallback) noexcept { return std::isfinite (v) ? v : fallback; }

    void apply (const DynamicEqBandParams& p) noexcept
    {
        // Sanitize the fields that feed the SVFs / gain math directly — a NaN freq reaches tan() and poisons
        // both filters. threshold/ratio/knee/range are guarded by the GainComputer; the followers' coeff()
        // guards attack/release. Sanitized into params_ so updateAudio()/dynamicDeltaDb() see safe values.
        params_.freq         = finite (p.freq, 1000.0);
        params_.Q            = finite (p.Q, 2.0);
        params_.staticGainDb = finite (p.staticGainDb, 0.0);
        side_.setParams (eq::FilterType::BandPass, params_.freq, params_.Q, 0.0);
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
    bool prepared_ = false;                                       // true only after a prepare() that succeeded
    DynamicEqBandParams params_;

    eq::Svf side_, audio_;                       // sidechain BandPass · audio Bell/shelf
    int     ranNc_ = 0;                          // channels that advanced state on the previous call
    dynamics::EnvelopeFollower    env_;
    dynamics::GainComputer        gc_;
    dynamics::GainReductionFollower gr_;

    dynamics::LinkMode link_ = dynamics::LinkMode::Max;
    float  sign_ = 1.0f;
    double curGainDb_ = 0.0;
    int    ksamp_ = 0, K_ = 16;
};

} // namespace felitronics::dynamiceq

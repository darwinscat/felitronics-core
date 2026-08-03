// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/eq/EqBand.h>
#include <felitronics/eq/Svf.h>
#include <felitronics/dynamics/EnvelopeFollower.h>
#include <felitronics/dynamics/GainComputer.h>
#include <felitronics/dynamics/GainReductionFollower.h>
#include <felitronics/dynamics/RelativeLevel.h>
#include <felitronics/dynamics/BandBallistics.h>
#include <felitronics/core/Math.h>

#include <algorithm>
#include <cmath>

namespace felitronics::dynamiceq
{

//==============================================================================
// felitronics::dynamiceq::LaneDynamics — drives the dynamics of ONE eq::EqBand point across its
// placement lanes. This is the composition layer: eq owns the filter and a per-lane gain-delta seam
// (eq::EqBand::setLaneDeltaDb), dynamics owns the detector/computer/ballistics primitives, and
// NEITHER knows about the other. Everything product-shaped lives here.
//
// AUTO-FIRST, by design. The user sets ONE number — dyn.rangeDb, whose SIGN is the direction. Ratio
// and knee are fixed internally; threshold is relative to the lane's own programme level; attack and
// release come from the band's own fc/Q. Nothing else is exposed, and none of it is guesswork the
// user has to sanity-check.
//
// SHARED SETTINGS, PER-LANE STATE. All lanes of a point share dyn, but each gets its own probe,
// follower, programme estimate and delta — because each lane has its own freq/Q and its own signal.
// That is only sound because the threshold is auto-relative: an absolute dBFS threshold shared
// between a Mid and a Side lane, 20+ dB apart, would be meaningless.
//
// THE SIDECHAIN IS THE SECTION INPUT, not the band's own input. In a 24-point chain each point's
// input is the previous points' OUTPUT, so probing that would let earlier points' moving deltas
// modulate later detectors at overlapping frequencies — chained pumping the slow estimator cannot
// absorb, because the feed-through is fast. Every probe taps the same pre-EQ signal instead.
//
// CONTROL RATE. Deltas are pushed every `kControl` samples; processBand() splits the buffer for you
// so the audio and the deltas stay in step. Feeding a whole 512-sample block at once would quantise
// a 1 ms attack to 10 ms.
//
// 🔴 RT-safe: prepare() allocates nothing beyond its members; process paths are alloc/lock/throw-free.
class LaneDynamics
{
public:
    static constexpr int kControl = 16;      // samples between delta updates (~0.33 ms at 48 k)

    void prepare (double sampleRate, int maxChannels) noexcept
    {
        fs_ = (std::isfinite (sampleRate) && sampleRate > 0.0) ? sampleRate : 48000.0;
        ch_ = std::clamp (maxChannels, 1, core::kMaxChannels);
        for (int i = 0; i < eq::kNumLanes; ++i)
        {
            LaneState& s = st_[i];
            s.probe.prepare (fs_, i == 0 ? ch_ : 1);   // only the ST lane is multi-channel
            s.env.prepare (fs_);
            s.env.setDetector (dynamics::Detector::Rms);
            s.env.setTimes (2.0, 2.0);                 // a short, symmetric LEVEL window: the musical
                                                       // timing lives in the GR follower, so the knee
                                                       // cannot warp the attack
            s.rel.prepare (fs_);
            s.gr.prepare (fs_);
            s.gc.setMode (dynamics::Mode::DownCompress);   // direction is the sign of range, applied outside
            s.gc.setRatio (kRatio);
            s.gc.setThresholdDb (0.0);                     // RelativeLevel wiring: never rewritten
        }
        reset();
    }

    void reset() noexcept
    {
        for (auto& s : st_)
        {
            s.probe.reset(); s.env.reset(); s.rel.reset(); s.gr.reset();
            s.deltaDb = 0.0;
            // Invalidate the retune sentinels too. Svf::prepare() resets STATE but keeps coefficients,
            // which were baked at the old sample rate (g = tan(pi*f/fs)); if the host re-prepares at a
            // new rate and the adapter re-sends identical params, "unchanged" would skip the redesign
            // and leave the probe listening an octave away from the band it drives.
            s.freq = -1.0; s.Q = -1.0;
        }
        engaged_ = false;
    }

    // Point parameters. A MATERIAL fc/Q change re-arms that lane's estimator adaptation window — a
    // gain change deliberately does not, and neither does any other write, or dragging a node would
    // flatten the estimator into an instantaneous follower.
    void setParams (const eq::BandParams& p) noexcept
    {
        for (int i = 0; i < eq::kNumLanes; ++i)
        {
            const eq::LaneParams& lp = p.lanes[i];
            LaneState& s = st_[i];

            const bool retuned = ! nearlyEqual (lp.freq, s.freq) || ! nearlyEqual (lp.Q, s.Q);
            if (retuned)
            {
                const double f = probeFreq (lp.freq), q = probeQ (lp.Q);
                s.probe.setParams (eq::FilterType::BandPass, f, q, 0.0);
                s.rel.retuned();
                const auto t = dynamics::BandBallistics::compute (fs_, f, q, p.dyn.atk, p.dyn.rel);
                s.gr.setTimes (t.attackMs, t.releaseMs);
                s.freq = lp.freq; s.Q = lp.Q;   // raw, so retune detection tracks what the caller sent
            }
            else if (! nearlyEqual (p.dyn.atk, dynAtk_) || ! nearlyEqual (p.dyn.rel, dynRel_))
            {
                const auto t = dynamics::BandBallistics::compute (fs_, probeFreq (lp.freq), probeQ (lp.Q),
                                                                  p.dyn.atk, p.dyn.rel);
                s.gr.setTimes (t.attackMs, t.releaseMs);
            }

            const double range = std::fabs (p.dyn.rangeDb);
            s.gc.setRangeDb (range);
            // Knee scaled off range so a small range keeps a proportionate corner instead of a 6 dB
            // smear, and offset by knee/2 so "idle on stationary programme" is EXACT: the steady-state
            // reduction is -slope*kneeOver(-offset), which is zero only once offset >= knee/2.
            const double knee = std::min (6.0, 1.5 * range);
            s.gc.setKneeDb (knee);
            s.offsetDb = knee * 0.5 + kHeadroomDb;
        }
        dynAtk_ = p.dyn.atk; dynRel_ = p.dyn.rel;
        dyn_ = p.dyn;
        sign_ = (p.dyn.rangeDb < 0.0) ? 1.0f : -1.0f;   // range<0 = cut when loud (DownCompress as-is)
    }

    // Audio + sidechain in, band driven out. `audio` is processed IN PLACE by `band`; `sidechain` is
    // the EQ section's common input. Split into control-rate chunks so a 1 ms attack means 1 ms.
    void processBand (float* const* audio, const float* const* sidechain,
                      int numChannels, int numSamples, eq::EqBand& band) noexcept
    {
        if (numSamples <= 0 || numChannels <= 0) return;
        const int nc = std::clamp (numChannels, 1, ch_);

        // Disengaging must not leave the band frozen mid-duck: zero the seams and drop the
        // detector/programme/GR state once, on the edge. Without this, toggling dynamics off during a
        // loud passage and back on during a quiet one replays the old gain reduction onto material
        // that asked for nothing.
        if (! dyn_.on || sidechain == nullptr)
        {
            if (engaged_) { disengage (band); engaged_ = false; }
            band.processBlock (audio, nc, numSamples);
            return;
        }
        engaged_ = true;

        int done = 0;
        while (done < numSamples)
        {
            const int n = std::min (kControl, numSamples - done);
            float* aud[core::kMaxChannels];
            const float* sc[core::kMaxChannels];
            for (int c = 0; c < nc; ++c) { aud[c] = audio[c] + done; sc[c] = sidechain[c] + done; }

            // AUDIO FIRST, then detect. Detecting this chunk before processing it would let a
            // transient at sample 15 alter output sample 0 — up to 15 samples of undeclared
            // look-ahead in a plugin that reports zero latency. Running the band on the delta derived
            // from the PREVIOUS chunk keeps the path causal; the cost is one control period (0.33 ms)
            // of delay on the gain, which the ballistics dwarf.
            band.processBlock (aud, nc, n);
            advance (sc, nc, n, band);
            done += n;
        }
    }

    // Live gain reduction of one lane (dB, signed) — for metering. Lock-free by being a plain read on
    // the audio thread's own value; a host/UI reader should publish it through its own atomic.
    double deltaDb (eq::Lane l) const noexcept { return st_[(std::size_t) l].deltaDb; }

private:
    struct LaneState
    {
        eq::Svf                         probe;      // sidechain band-pass at this lane's freq/Q
        dynamics::EnvelopeFollower      env;
        dynamics::RelativeLevel         rel;
        dynamics::GainComputer          gc;
        dynamics::GainReductionFollower gr;
        double freq = -1.0, Q = -1.0;               // last applied, for retune detection
        double offsetDb = 0.0;                      // knee/2 + headroom: makes "idle" exact
        double deltaDb = 0.0;
    };

    static constexpr double kRatio      = 4.0;   // fixed: the knob is range, not ratio
    static constexpr double kHeadroomDb = 3.0;   // how far above "normal" a peak must sit to engage

    static bool nearlyEqual (double a, double b) noexcept
    {
        return std::isfinite (a) && std::isfinite (b) && std::fabs (a - b) <= 1.0e-9 * std::fmax (1.0, std::fabs (a));
    }

    // The SAME rails eq::EqBand applies to itself, so the detector, the ballistics and the filter can
    // never describe three different bands — and so a NaN freq cannot reach tan() and poison the
    // follower permanently (env = in + c*(env-in) never recovers from NaN).
    double probeFreq (double f) const noexcept
    {
        return std::clamp (std::isfinite (f) ? f : 1000.0, 10.0, 0.49 * fs_);
    }
    static double probeQ (double q) noexcept
    {
        return std::clamp (std::isfinite (q) ? q : 1.0, 0.05, 40.0);
    }

    // Leave the band exactly as an opted-out one: no residual delta, no state that could resume.
    void disengage (eq::EqBand& band) noexcept
    {
        for (int i = 0; i < eq::kNumLanes; ++i)
        {
            st_[i].deltaDb = 0.0;
            st_[i].gr.reset(); st_[i].env.reset(); st_[i].probe.reset(); st_[i].rel.reset();
            band.setLaneDeltaDb ((eq::Lane) i, 0.0);
        }
    }

    // One control-rate chunk: run every running lane's detector over the SECTION INPUT and push its
    // delta into the band.
    void advance (const float* const* sc, int nc, int n, eq::EqBand& band) noexcept
    {
        const eq::BandParams& p = band.params();
        for (int i = 0; i < eq::kNumLanes; ++i)
        {
            const eq::Lane l = (eq::Lane) i;
            if (! laneRuns (p, l, nc)) { st_[i].deltaDb = 0.0; band.setLaneDeltaDb (l, 0.0); continue; }

            LaneState& s = st_[i];
            float linked = 0.0f;
            for (int k = 0; k < n; ++k)
            {
                // Channel LINKING is mandatory on the ST lane: probing each channel independently
                // would compress L and R by different amounts and the image would wander on every
                // sibilant. The single-signal lanes derive their own domain first.
                float e = 0.0f;
                if (l == eq::Lane::Stereo)
                {
                    for (int c = 0; c < nc; ++c)
                        e = std::fmax (e, std::fabs (s.probe.processSample (c, sc[c][k])));
                }
                else
                {
                    const float x = laneSignal (l, sc, nc, k);
                    e = std::fabs (s.probe.processSample (0, x));
                }
                linked = s.env.process (e);
                s.rel.accumulate (linked);
            }
            s.rel.update (n);

            const double levelDb = core::gainToDb (std::fmax ((double) linked, 1.0e-9));
            // Gate the COMPUTER'S INPUT, not merely the estimator: with a soft knee a relative level
            // of 0 still asks for knee/8 * slope, so an up-lifting band would ride the noise floor up
            // through every pause.
            const double target = s.rel.activity()
                                ? (double) sign_ * s.gc.deltaDb (s.rel.relativeDb (levelDb) - s.offsetDb)
                                : 0.0;

            // The GR follower's coefficients are PER SAMPLE, so it must be advanced once per sample.
            // Calling it once per control chunk would stretch every time constant by kControl — a
            // nominal 1 ms attack would behave like 16 ms, and the amount would depend on the host's
            // block size (a 17-sample block gives a different factor than a 64-sample one). The
            // target is what updates at control rate; the smoothing does not.
            float smoothed = (float) s.deltaDb;
            for (int k = 0; k < n; ++k) smoothed = s.gr.process ((float) target);
            s.deltaDb = (double) smoothed;
            band.setLaneDeltaDb (l, s.deltaDb);

            s.probe.flushDenormals(); s.env.flushDenormals(); s.rel.flushDenormals(); s.gr.flushDenormals();
        }
    }

    static bool laneRuns (const eq::BandParams& p, eq::Lane l, int nc) noexcept
    {
        if (! p.on || p.bypass) return false;
        // EqBand's swept branch runs its search SVF and never applies the delta, so a detector here
        // would drive a seam the audio ignores — and deltaDb() would report gain reduction that is
        // not happening. Same predicate EqBand::sweptActive uses, via the shared free function, so
        // the two gates cannot drift apart.
        if (p.swept && eq::onlyStereoEnabled (p)) return false;
        const eq::LaneParams& lp = p.lane (l);
        if (! lp.on || lp.bypass) return false;
        return l == eq::Lane::Stereo || nc == 2;      // L/R/M/S are stereo-only, as in EqBand
    }

    static float laneSignal (eq::Lane l, const float* const* sc, int nc, int k) noexcept
    {
        const float L = sc[0][k];
        const float R = nc > 1 ? sc[1][k] : L;
        switch (l)
        {
            case eq::Lane::Left:  return L;
            case eq::Lane::Right: return R;
            case eq::Lane::Mid:   return 0.5f * (L + R);
            case eq::Lane::Side:  return 0.5f * (L - R);
            default:              return L;
        }
    }

    double        fs_ = 48000.0;
    int           ch_ = 2;
    eq::DynParams dyn_;
    double        dynAtk_ = 0.5, dynRel_ = 0.5;
    float         sign_ = 1.0f;
    bool          engaged_ = false;      // was the dynamics path live last call? (edge detect)
    LaneState     st_[eq::kNumLanes];
};

} // namespace felitronics::dynamiceq

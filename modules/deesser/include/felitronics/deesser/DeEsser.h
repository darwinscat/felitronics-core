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
    [[nodiscard]] bool prepare (double sampleRate, int maxBlock, int maxChannels) noexcept
    {
        prepared_ = false;                                       // any early return below leaves it unprepared
        if (maxChannels < 1 || maxChannels > core::kMaxChannels) return false;   // law 11(b): BINDING
        fs_ = sampleRate > 0.0 ? sampleRate : 48000.0;
        channels_ = maxChannels;
        side_.prepare (fs_, channels_);
        xover_.prepare (fs_, channels_);
        env_.prepare (fs_); gr_.prepare (fs_);
        (void) maxBlock;
        if (! deq_.prepare (fs_, channels_)) return false;
        apply (params_);
        reset();
        prepared_ = true;
        return true;
    }

    bool isPrepared() const noexcept { return prepared_; }

    void reset() noexcept
    {
        side_.reset(); xover_.reset(); env_.reset(); gr_.reset(); deq_.reset(); grDb_ = 0.0f;
        ranNc_ = 0;   // nothing has run, so nothing can be stopping (see dropStoppedChannels)
    }
    void setParams (const DeEsserParams& p) noexcept
    {
        const bool relatch = (p.mode != params_.mode) || (p.listen != params_.listen);   // topology / signal-path change
        params_ = p; apply (p);
        if (relatch) reset();   // clear the now-inactive path's stale state → no thump when toggling mode/listen live
    }
    static constexpr int latencySamples() noexcept { return 0; }
    double gainReductionDb() const noexcept { return grDb_; }

    // Law 11 (DSP-ARCHITECTURE.md §2) — see the law for what each refusal means.
    [[nodiscard]] bool process (float* const* io, int numChannels, int n) noexcept
    {
        if (numChannels < 0 || n < 0) return false;
        if (! prepared_) return false;
        if (numChannels > channels_) return false;                // width is a LIMIT — law 11(b)
        if (n == 0) return true;                                  // no samples: no time, no edge
        const int nc = numChannels;
        dropStoppedChannels (nc);                                 // law 11(d): the edge is clocked by n
        // ...and the gap has to reach the INNER band too: its own per-channel SVFs freeze the same way.
        // Measured on a 7 kHz tone through the default DynamicEq mode: 0.182 out of DIGITAL SILENCE
        // (-14.8 dBFS) when this returned early instead of passing the zero width down.
        //
        // LAW 11c — AND IT HAS TO REACH THE HALF THAT IS ACTUALLY RUNNING. This used to forward the gap
        // to `deq_` whatever the mode was, which was harmless while both halves merely froze and stops
        // being harmless the moment they advance: in SplitBand (or with `listen` on) the live detector is
        // the LOCAL one below — `side_`, `env_`, `gr_`, `xover_` — and `deq_` is not in the signal path
        // at all, so forwarding only there would have advanced the idle child and left the working
        // detector frozen. Measured on the freeze this replaces, SplitBand at 7 kHz: -8.000 dB of gain
        // reduction held through a second of gap where silence releases to -1.104, and to -0.000 through
        // ten. The meter is written on both routes, because a stale `gainReductionDb()` after a gap is
        // the same defect one storey up.
        if (nc == 0)
        {
            if (params_.mode == DeEsserMode::DynamicEq && ! params_.listen)
            {
                const bool ok = deq_.process (io, 0, n);
                grDb_ = (float) deq_.dynamicDeltaDb();
                return ok;
            }
            advanceSilence (n);
            return true;
        }

        // DynamicEq mode = the surgical dynamic-EQ band (no split, no listen detour).
        if (params_.mode == DeEsserMode::DynamicEq && ! params_.listen)
        {
            const bool ok = deq_.process (io, numChannels, n);
            grDb_ = (float) deq_.dynamicDeltaDb();
            return ok;
        }

        for (int i = 0; i < n; ++i)
        {
            // detector: per-channel sidechain BandPass → one linked level
            float sc[core::kMaxChannels]; float linked = 0.0f, sq = 0.0f;
            for (int c = 0; c < nc; ++c)
            {
                // GATED ON THE FILTER'S INPUT, not on its output. The filter is recursive, so a
                // poisoned sample reaching it is permanent — and a gate placed AFTER it would read
                // the resulting NaN as silence, turning "noise for ever" into "silence for ever"
                // rather than fixing anything. Measured: gate-after diverges from a sanitised
                // reference by 1222..5533 samples and 0.2..171 dB depending on the caller's block
                // size; gate-before diverges by zero, at every block size.
                sc[c] = side_.processSample (c, dynamics::detectorGate (io[c][i]));
                const float a = std::fabs (sc[c]);
                if (link_ == dynamics::LinkMode::Max) { if (a > linked) linked = a; } else sq += sc[c] * sc[c];
            }
            if (link_ == dynamics::LinkMode::MeanPower) linked = std::sqrt (sq / (float) nc);

            const float gain = sharedStep (linked);

            for (int c = 0; c < nc; ++c)
            {
                if (params_.listen) { io[c][i] = sc[c]; continue; }
                float lo, hi; xover_.processSample (c, io[c][i], lo, hi);
                io[c][i] = lo + gain * hi;                          // duck only the high band
            }
        }
        side_.flushDenormals(); xover_.flushDenormals(); env_.flushDenormals(); gr_.flushDenormals();
        return true;
    }

private:
    // LAW 11c — `n` samples of DIGITAL SILENCE through the SplitBand/listen detector. The per-channel
    // sidechain columns were dropped by `dropStoppedChannels(0)` just above, so the linked level of a
    // silent band is exactly +0.0f; it is PINNED rather than run through the loop's own link, whose
    // MeanPower branch would divide by a channel count of zero. Every remaining line is the line the
    // audio loop runs, in the same order and with the same types — note this stage does NOT narrow the
    // dB to float before the curve, where `DynamicEqBand` does, and that difference is preserved here.
    // `xover_` is per channel and has no columns left to advance, so only its flush is owed.
    void advanceSilence (int n) noexcept
    {
        for (int i = 0; i < n; ++i)
        {
            const float e0 = env_.stateWord(), g0 = gr_.valueDb();
            (void) sharedStep (0.0f);       // a silent band links to exactly +0.0f at every width
            if (core::sameBits (e0, env_.stateWord()) && core::sameBits (g0, gr_.valueDb())) break;
            if (! std::isfinite (env_.stateWord()) && ! std::isfinite (e0)
                && ! std::isfinite (gr_.valueDb()) && ! std::isfinite (g0)) break;
        }
        side_.flushDenormals(); xover_.flushDenormals(); env_.flushDenormals(); gr_.flushDenormals();
    }

    // ONE sample of the SHARED body — the detector level in, the duck gain out, the meter written on the
    // way. Written once and called from both loops. NB this stage does NOT narrow the dB to float before
    // the curve where `DynamicEqBand` does, and it floors the level at 1e-9 where `GainReductionPath`
    // floors at 1e-12: three neighbouring stages, three spellings, which is exactly why the silent path
    // reuses this one instead of growing a fourth.
    float sharedStep (float linked) noexcept
    {
        const float lvl = env_.process (linked);
        const float gr  = gr_.process ((float) gc_.deltaDb ((double) core::gainToDb ((double) std::max (lvl, 1.0e-9f))));   // ≤ 0 (cut)
        grDb_ = gr;
        return (float) core::dbToGain ((double) gr);
    }

    static double finite (double v, double fallback) noexcept { return std::isfinite (v) ? v : fallback; }

    void apply (const DeEsserParams& p) noexcept
    {
        // Sanitize what feeds the SVF/crossover directly (a NaN fc reaches tan()); threshold/ratio/knee/range
        // are guarded by the GainComputer, attack/release by the followers' coeff(). Written back to params_
        // so a later relatch/reset sees safe values too.
        params_.fc  = finite (p.fc, 7000.0);
        params_.scQ = finite (p.scQ, 2.0);
        side_.setParams (eq::FilterType::BandPass, params_.fc, params_.scQ, 0.0);
        xover_.setFrequency ((float) params_.fc);
        env_.setDetector (p.detector); env_.setTimes (2.0, 2.0);
        gc_.setMode (dynamics::Mode::DownCompress);
        gc_.setThresholdDb (p.thresholdDb); gc_.setRatio (p.ratio); gc_.setKneeDb (p.kneeDb); gc_.setRangeDb (p.rangeDb);
        gr_.setTimes (p.attackMs, p.releaseMs);
        link_ = p.link;

        dynamiceq::DynamicEqBandParams dp;
        dp.type = eq::FilterType::Bell; dp.freq = params_.fc; dp.Q = params_.scQ; dp.mode = dynamiceq::DynamicEqMode::CutWhenLoud;
        dp.thresholdDb = p.thresholdDb; dp.ratio = p.ratio; dp.kneeDb = p.kneeDb; dp.rangeDb = p.rangeDb;
        dp.attackMs = p.attackMs; dp.releaseMs = p.releaseMs; dp.detector = p.detector; dp.link = p.link;
        deq_.setParams (dp);
    }

    double fs_ = 48000.0;
    int channels_ = 2;
    bool prepared_ = false;                                       // true only after a prepare() that succeeded
    DeEsserParams params_;

    // The sidechain band-pass and the crossover are both PER CHANNEL and both advance only for c < nc,
    // while the envelope and the GR follower they feed are shared. In SplitBand mode a channel that
    // leaves and returns therefore replays two frozen recursions into the new stream — measured 1.706e-01
    // (-15.4 dBFS) on the returning channel at sample 0, out of DIGITAL SILENCE, against 0.000 for the
    // channel that stayed. DynamicEq mode was already safe because it delegates to DynamicEqBand, which
    // carries its own ledger; this suite's own tests never processed more than one channel, which is why
    // a green suite said nothing about it. The shared detector half is deliberately untouched: it is
    // supposed to follow whichever channels are actually present.
    void dropStoppedChannels (int nc) noexcept
    {
        for (int c = nc; c < ranNc_; ++c) { side_.resetChannel (c); xover_.resetChannel (c); }
        ranNc_ = nc;
    }

    int     ranNc_ = 0;                          // channels that advanced state on the previous call
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

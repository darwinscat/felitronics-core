// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/core/Math.h>
#include <felitronics/dynamics/ChannelLinker.h>
#include <felitronics/dynamics/GainReductionPath.h>
#include <felitronics/dynamics/offline/Quantile.h>

#include <cmath>
#include <cstdint>
#include <limits>

namespace felitronics::dynamics::offline
{

//==============================================================================
// DETECTOR-DOMAIN ANALYSIS — the statistics of the envelope THE COMPRESSOR WILL ACTUALLY SEE, and the
// threshold that delivers a wanted amount of gain reduction on a given piece of material.
//
// WHY THIS EXISTS. The predecessor pipeline took the 95th percentile of a K-weighted three-second
// short-term LUFS and put it in the threshold of an unweighted 5 ms RMS detector. Different weighting,
// different window, different crest factor: a category error, and the pile of correction offsets on top
// of it was the symptom. The fix is to measure in the detector's own domain — same `Detector`, same
// `rmsWindowMs`, same `LinkMode`, same already-filtered key — which since P2 is a structural fact
// rather than a convention: `GainReductionParams` derives from `DetectorParams`, so the analysis and
// the compressor are configured by the same object, not by two that must be kept in step.
//
// AND WHY THE OBVIOUS INVERSION IS WRONG. It is tempting to invert the static curve at a level
// percentile: pick the level the signal exceeds 5 % of the time, ask the curve what reduction that
// earns, done. That cannot work, and the reason is topological rather than numerical: THE BALLISTICS
// SIT AFTER THE CURVE and integrate over time, while a percentile of the level cannot see time at all.
// Measured counterexample, five seconds of Peak-detector levels, 20 % at -4 dB and 80 % at -40 dB,
// threshold -20, ratio 4, hard knee, attack = release = 5 ms: laid out as ONE CONTIGUOUS SECOND the p95
// gain reduction is 11.9999 dB; laid out as EVERY FIFTH SAMPLE it is 2.4200 dB. The two signals have
// bit-identical level percentiles. No function of those percentiles can tell them apart.
//
// So the threshold is found by BISECTION ON THE MEASURED QUANTILE, running the whole
// `LinkedDetector → GainComputer → GainReductionFollower` chain — as one `GainReductionPath`, the same
// object the compressor runs — over the key at each probe. The level percentile survives only as a
// SEED. Measured on an 18-track release corpus at three targets (1, 3 and 6 dB — 54 rows), the
// static-curve inversion in its strongest honest form misses by up to 1.39 dB, mean 0.49, and 41 of
// the 54 rows miss by more than 0.3 dB. Swept instead across detector, attack, release and ratio on
// one signal (another 54 rows), the worst miss is 2.19 dB and 34 rows exceed 0.3. Two different
// sweeps, so neither number is the other's; both say the same thing. That is what the apparatus buys.
//
// WHAT IS MEASURED, stated rather than implied:
//   * THE INPUT IS THE KEY, already filtered by the product, exactly as the compressor receives it.
//     No sidechain EQ in here (ADR §4), and the same filtered key must go to both.
//   * EVERY FRAME ADVANCES THE STATE. A statistics gate excludes samples from the DISTRIBUTION; it
//     never skips the detector or the follower, because the state a skipped sample would have left is
//     the state the next counted one starts from.
//   * THE GATE IS ON THE DETECTOR LEVEL, not on the resulting gain reduction, and that is what keeps
//     the solver's objective monotone: the level does not depend on the threshold, so the set of
//     counted samples is the same at every probe. Gating on the RESULT would make the objective
//     select its own sample set and the bisection would be searching a moving target. Default off.
//   * SILENCE IS INCLUDED BY DEFAULT, and it moves the answer: with half the material silent, the p95
//     of the whole is the p90 of the part that sounds. Excluding it is a PRODUCT decision with a
//     number attached (the gate), not a default the core may take on the product's behalf.
//   * -240 dB IS A FLOOR, NOT A MEASUREMENT. `core::gainToDb` clamps at 1e-12, so digital silence
//     enters as -240 dB; `atFloorFrames` says how much of the answer is that floor. And the floor has
//     teeth: a threshold BELOW -240 dB makes silence an ACTIVE sample (measured: 45 dB of reduction on
//     digital silence at threshold -300, ratio 4, range 60), which is why the solver reasons in dB
//     throughout and never in the linear domain, where no such floor exists.
//   * THE STATISTIC IS THE MAGNITUDE |grDb|. On a signed gain reduction the 95th percentile picks the
//     LEAST processed end, which is the opposite of what anyone asking for it means. For
//     `DownCompress` and `DownExpand` the magnitude is the reduction; for `UpCompress` it is a boost,
//     and the honest general name is the amount of processing.
//   * IDENTITY IS WITHIN A BUILD. Same object, same parameters, same rate, same binary ⇒ the same
//     bits; across toolchains the desktop tier declares `-ffp-contract=on` and libm is not
//     bit-portable, so it is equality to rounding.
//
// OFFLINE / MESSAGE THREAD ONLY: prepare() allocates. Nothing here is RT-safe and nothing here needs
// to be — the compressor is what runs on the audio thread.
//==============================================================================

// The statistics of the detector envelope. Levels are dB relative to full scale.
struct EnvelopeStats
{
    double p50Db = 0.0, p90Db = 0.0, p95Db = 0.0, p99Db = 0.0;
    double minDb = 0.0, maxDb = 0.0;              // of the COUNTED frames

    // The extremes of every frame the detector saw, gate or no gate. They are NOT the same as the pair
    // above and the difference is load-bearing: an excluded frame still drives the follower, so
    // anything sizing a search over thresholds has to reason about the levels that actually occurred,
    // not about the subset somebody chose to count. A bracket built from the gated extremes can fail to
    // contain the answer — measured, and by a wide margin.
    double minSeenDb = 0.0, maxSeenDb = 0.0;

    // Crest factor of the GATED, LINKED KEY — peak over RMS of what enters the envelope follower, not
    // of what leaves it. The crest of a smoothed envelope is a property of the smoothing at least as
    // much as of the signal, and would report a different number for the same music at a different
    // `rmsWindowMs`; this one does not.
    double crestDb = 0.0, peakDb = 0.0, rmsDb = 0.0;

    std::uint64_t frames        = 0;   // frames that entered the distribution (after the gate)
    std::uint64_t framesSeen    = 0;   // frames processed, gate or no gate
    std::uint64_t atFloorFrames = 0;   // frames whose level was AT the -240 dB floor: censored
    bool          valid         = false;
};

//==============================================================================
// felitronics::dynamics::offline::EnvelopeAnalyzer — streams a key through the compressor's own
// detector and reports the statistics of its envelope. O(1) in the length of the material: it holds
// histograms, not the signal. `analyze()` may be called repeatedly to feed the key in chunks; the
// state carries across, exactly as it would in the compressor.
class EnvelopeAnalyzer
{
public:
    static constexpr double kFloorDb = -240.0;    // core::gainToDb's clamp, spelled out

    // `binDb` is the resolution every reported percentile carries. Allocates, and REFUSES a resolution
    // it cannot honour over the whole level range rather than covering part of it.
    [[nodiscard]] bool prepare (double sampleRate, double binDb = 0.01)
    {
        prepared_ = false;
        if (! (sampleRate > 0.0) || ! std::isfinite (sampleRate)) return false;
        fs_ = sampleRate;
        det_.prepare (fs_);
        // Levels run from the -240 dB floor to +120 dBFS, which is where the detector gate clamps the
        // key (ChannelLinker's 1e6). Nothing can be reported outside that by construction.
        if (! hist_.prepare (kFloorDb, 120.0, binDb)) return false;
        reset();
        prepared_ = true;
        return true;
    }

    void setParams (const DetectorParams& p) noexcept { det_.setParams (p); }
    const DetectorParams& params() const noexcept { return det_.params(); }

    // Excludes frames whose DETECTOR LEVEL is at or below this from the distribution — they still
    // advance the detector. -infinity (the default) counts everything.
    void setStatisticsGateDb (double gateDb) noexcept { gateDb_ = gateDb; }
    double statisticsGateDb() const noexcept { return gateDb_; }

    void reset() noexcept
    {
        det_.reset(); hist_.reset();
        peak_ = 0.0; sumSq_ = 0.0; linked_ = 0; seen_ = 0; atFloor_ = 0;
        minSeen_ = 0.0; maxSeen_ = 0.0;
    }

    // One chunk of the key. PRECONDITION: `key` non-null, `numKeyChannels >= 1`, `numSamples >= 0`.
    //
    // `envOut` optionally receives the LINKED DETECTOR LEVEL of every frame — the linear amplitude, the
    // exact float the compressor's detector produces. It is written only when `envCapacity` covers the
    // whole chunk; a short buffer writes NOTHING rather than a prefix, for the reason the compressor
    // refuses a short gain-reduction tap: a partial trace looks like data. Recording it costs nothing
    // extra here, and it is what lets a search reuse one detector pass across many probes.
    // INDEXED WITHIN THE CHUNK: `envOut[0]` is this call's first frame, not the stream's. A caller
    // feeding the key in pieces advances the pointer itself, the same way it advances the key.
    void analyze (const float* const* key, int numKeyChannels, int numSamples,
                  float* envOut = nullptr, int envCapacity = 0) noexcept
    {
        if (! prepared_ || key == nullptr || numKeyChannels < 1 || numSamples <= 0) return;
        if (envOut != nullptr && envCapacity < numSamples) envOut = nullptr;
        for (int i = 0; i < numSamples; ++i)
        {
            // The crest is taken on the linked key BEFORE the follower — same gate, same link, so it
            // is the same quantity the detector integrates, just unsmoothed.
            const double a = (double) linkAmplitudeGated (det_.params().link, key, numKeyChannels, i);
            if (a > peak_) peak_ = a;
            sumSq_ += a * a;
            ++linked_;

            const float  level   = det_.process (key, numKeyChannels, i);
            if (envOut != nullptr) envOut[i] = level;
            const double levelDb = core::gainToDb (level);
            if (seen_ == 0) { minSeen_ = maxSeen_ = levelDb; }
            else { if (levelDb < minSeen_) minSeen_ = levelDb; if (levelDb > maxSeen_) maxSeen_ = levelDb; }
            ++seen_;
            if (levelDb <= kFloorDb) ++atFloor_;
            if (levelDb > gateDb_) hist_.add (levelDb);
        }
        det_.flushDenormals();
    }

    EnvelopeStats stats() const noexcept
    {
        EnvelopeStats s;
        s.framesSeen = seen_; s.atFloorFrames = atFloor_; s.frames = hist_.count();
        s.minSeenDb = minSeen_; s.maxSeenDb = maxSeen_;
        if (hist_.count() == 0) return s;
        // `valid` MEANS the numbers below are numbers. A quantile that could not be reported is not a
        // zero — reporting one as if it were is how a histogram too coarse for its range once produced
        // a confident p95 of 0 dB for a signal at -20.
        bool okAll = hist_.quantile (0.50, s.p50Db);
        okAll = hist_.quantile (0.90, s.p90Db) && okAll;
        okAll = hist_.quantile (0.95, s.p95Db) && okAll;
        okAll = hist_.quantile (0.99, s.p99Db) && okAll;
        if (! okAll) return s;
        s.minDb = hist_.minValue(); s.maxDb = hist_.maxValue();
        s.peakDb = core::gainToDb (peak_);
        s.rmsDb  = core::gainToDb (linked_ > 0 ? std::sqrt (sumSq_ / (double) linked_) : 0.0);
        s.crestDb = s.peakDb - s.rmsDb;
        s.valid = true;
        return s;
    }

    // The fraction of counted frames whose level is above `levelDb`, to within one histogram bin.
    double fractionAboveDb (double levelDb) const noexcept { return hist_.fractionAbove (levelDb); }

    const QuantileHistogram& histogram() const noexcept { return hist_; }
    bool isPrepared() const noexcept { return prepared_; }

private:
    LinkedDetector    det_;
    QuantileHistogram hist_;
    double fs_ = 48000.0, gateDb_ = -std::numeric_limits<double>::infinity();
    double peak_ = 0.0, sumSq_ = 0.0, minSeen_ = 0.0, maxSeen_ = 0.0;
    std::uint64_t linked_ = 0, seen_ = 0, atFloor_ = 0;
    bool prepared_ = false;
};

} // namespace felitronics::dynamics::offline

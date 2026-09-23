// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.
#pragma once

#include <felitronics/core/Config.h>
#include <felitronics/core/DetMath.h>
#include <felitronics/core/Math.h>
#include <felitronics/oversampling/PolyphaseOversampler.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

namespace felitronics::analysis
{

//======================================================================================================
// HOW FAR, HOW OFTEN AND FOR HOW LONG A RENDER GOES OVER A CEILING — the measurement behind the choice
// between a clipper and a limiter. The caller supplies a render made WITHOUT the limiter and the ceiling
// it intends to deliver at; this reports the excursions that would have to be removed.
//
// WHY THE RECONSTRUCTION AND NOT THE SAMPLES, which is the first thing to get right. A dBTP ceiling is a
// statement about the RECONSTRUCTED signal: it exists between the samples. Two samples of 0.85 have a
// sample peak of -1.41 dBFS and a reference true peak of +0.42 dBTP, so a finder working on the grid
// reports NOTHING for an excursion of 1.4 dB. The under-read of a single crest of local frequency f is
// 20*log10(cos(pi*f/fs)): 0.000 dB at 60 Hz, -0.019 at 1 kHz, -2.011 at 10 kHz, -5.105 at 15 kHz (48 kHz).
// Bass crests survive the grid intact and transients do not — so a grid instrument would report bass
// excursions and hide the fast ones, which is precisely backwards for this decision. Every run here is
// found on the same 4x / 32-taps-per-phase reconstruction the delivery certificate is issued on.
//
// THE CONSTANTS ARE `ReferenceTruePeakMeter`'s, SPELLED. `PolyphaseOversampler`'s own default is 64 taps
// per phase; taking it would build an instrument that disagrees with the certificate by the instrument
// gap the true-peak tests pin. The pair is asserted equal in the suite rather than trusted.
//
// AND THE CERTIFICATE HAS A SAMPLE-PEAK FLOOR THIS CANNOT REPRODUCE. `truePeakLinear()` is
// `max(reconstructed, samplePeak)`, so a certificate can exceed the ceiling because of a raw SAMPLE while
// the reconstruction never does. A run has no meaning in that case — there is no interval over the
// ceiling, only one grid point at it — so none is invented: `samplePeakLinear()` and `reconstructedPeak()`
// are published separately and a caller comparing this instrument with the certificate can see which of
// the two carried the reading.
//
// WHAT A RUN IS. Per oversampled sample the excess is `e(k) = max over channels of (|y_c(k)| - T, 0)`,
// and a sample is ABOVE when `e(k) > 0` (strictly). Runs are the UNION ACROSS CHANNELS: a stereo-linked
// limiter acts once on the programme when either channel exceeds, and per-channel runs would count a
// centred kick twice and make the count depend on panning. `aboveOs(c)` is published per channel beside
// it so a one-sided problem stays visible.
//
// MERGING, AND WHAT IT DOES NOT TOUCH. Two excursions separated by a gap of at most `mergeMs` are one
// run. That parameter is meant to be swept, so `runCount`, the duration histogram and the percentile all
// move with it — but `occupancy` and the total dose DO NOT: both are accumulated over raw above-samples
// and are merge-free by construction. A consumer sweeping the window therefore keeps two fixed points.
//   * THE WINDOW HAS A PHYSICAL CEILING. Two crests of a tone at f are 1/(2f) apart, so a window of
//     4.2 ms fuses every half-cycle of 120 Hz bass into a single run and the percentile jumps by an order
//     of magnitude. Keep `mergeMs` well under `500/f_low` milliseconds; the class does not refuse it,
//     because the number is a property of the material and not of the instrument.
//
// DOSE IS PER-SAMPLE LINEAR EXCESS, `sum of e(k)` over the run, and that choice is load-bearing. It is
// additive across runs and programmes, it is merge-invariant in total (a gap contributes zero), and it is
// what a limiter's gain computer actually removes. "Peak times duration, once per run" over-reads a
// parabolic crest by 3/2 and a triangular one by 2, and counts merged gaps as excess. dB excess is not
// additive at all: a 2 dB run is not two 1 dB runs' worth of removed material.
//
// THE PERCENTILE NEVER READS THE RUN LIST. It comes from a fixed-time duration histogram, so it stays
// exact-to-the-bin after the list is full — law 11's "exhaustion is data" applied to a statistic rather
// than to a coordinate. Same for every class count, every dose and both ceiling densities.
//
// WHAT IT MEANS WHEN THERE ARE NO RUNS. Everything reads its canonical zero and `runCount()` is 0. A rule
// of the form "p90 <= 2 ms" is satisfied TRIVIALLY by a programme with no excursions at all, so a caller
// must test `runCount() > 0` beside it: "nothing went over" and "things went over but briefly" are
// different findings and are published differently.
class PeakExcursions
{
public:
    // THE RECONSTRUCTION IS THE CERTIFICATE'S — see the header note. Not the oversampler's defaults.
    static constexpr int kFactor       = 4;
    static constexpr int kTapsPerPhase = 32;

    // The interpolator's group delay in OVERSAMPLED samples, doubled so it stays an integer: the FIR is
    // symmetric over `kFactor*kTapsPerPhase` taps, so the delay is (128-1)/2 = 63.5 oversampled samples,
    // which is 15.875 INPUT samples and not 63.5 of them. Kept doubled because the half is real: rounding
    // it away is how an instrument acquires a bias nobody can find afterwards.
    static constexpr int kDelayOsX2 = kFactor * kTapsPerPhase - 1;
    static constexpr int kDelayOs   = kDelayOsX2 / 2;

    static constexpr int kChunk = 1024;                  // the interpolator's scratch, in base samples

    static constexpr double kMaxSampleRate = 768000.0;   // the ceiling every sibling analyzer carries
    static constexpr double kMaxMergeMs    = 1000.0;
    static constexpr double kMaxClassMs    = 1000.0;
    static constexpr int    kClassEdges    = 4;          // …so five classes: four edges and "longer"
    static constexpr int    kClasses       = kClassEdges + 1;

    // THE DURATION HISTOGRAM IS IN TIME, NOT IN OVERSAMPLED SAMPLES, and the difference is not cosmetic.
    // Bins of oversampled samples would cover 2.7 ms at 768 kHz — under the longest class edge — so the
    // percentile would saturate on a legal rate. 1/32 ms bins over 2048 bins is 64 ms at every rate.
    static constexpr int    kDurationBins  = 2048;
    static constexpr double kDurationBinMs = 1.0 / 32.0;

    // Local-maximum levels, for the ceiling densities. 0.05 dB bins so the 0.2 dB window is exactly four
    // of them, over [-100, +12] dBFS: 2240 bins.
    static constexpr int    kCeilingBins   = 2240;
    static constexpr double kCeilingBinDb  = 0.05;
    static constexpr double kCeilingTopDb  = 12.0;
    static constexpr double kCeilingSpanDb = 112.0;

    // Crest frequency, in third-octaves from 20 Hz: 20 * 2^(k/3) for k in [0, 30) reaches 20 kHz.
    static constexpr int    kCrestBins     = 30;
    static constexpr double kCrestBaseHz   = 20.0;

    struct Params
    {
        // The ceiling this render is judged against, in dBTP on the reconstruction above. Finite and
        // within +-200 dB so that `det::pow10` of a twentieth of it is finite.
        double thresholdDbtp = -1.0;
        // Gaps of at most this are bridged. 0 means no merging at all — every above-span is its own run.
        double mergeMs       = 1.0;
        // Duration class edges, strictly increasing, in milliseconds. The last class is "longer than the
        // last edge". They sit at the same scale as `mergeMs` on purpose, and the consequence is worth
        // saying out loud: the shortest class is the one merging empties.
        double classEdgesMs[kClassEdges] = { 0.5, 1.0, 2.0, 5.0 };
        // Capacity of the run list. Every AGGREGATE survives its exhaustion; only per-run coordinates
        // past the prefix are lost, and `runsComplete()` says so.
        int    maxRuns       = 1 << 16;
    };

    enum class Reason : std::uint8_t
    {
        Ok             = 0,
        NotFinished    = 1,   // finish() has not been called
        Empty          = 2,   // no sample was measured
        NonFiniteInput = 3,   // a non-finite sample reached the instrument; see nonFiniteSamples()
    };

    // One excursion. Coordinates are OVERSAMPLED samples on the programme's own grid — the interpolator's
    // 63.5-sample lag is already taken off, so `startOs / kFactor` is the input-sample position and the
    // remaining eighth of an input sample is a real residual, named rather than rounded away.
    struct Run
    {
        std::int64_t startOs  = 0;     // first above-sample, oversampled, since reset()
        std::int64_t lengthOs = 0;     // first above to last above INCLUSIVE — merged gaps count here
        std::int64_t aboveOs  = 0;     // …of which actually above the ceiling; equals lengthOs when mergeMs is 0
        double       peak     = 0.0;   // the largest |y| in the run, linear
        double       dose     = 0.0;   // sum over above-samples of (|y| - T), linear
        bool closedByFinish   = false; // the programme ended while it was still above
        // The crest's own frequency, from its SHAPE and nothing else. A crest of peak A over a ceiling T
        // by e = A - T for a duration D is, to second order, A*cos(wt): e = A*w^2*D^2/8, so
        // w = sqrt(8e / (A*D^2)). No filter, no phase, no window, and gain-invariant because e/A is a
        // ratio. Measured on the reconstruction at 48 kHz: a 60 Hz sine reads 59.3 Hz, a 40 Hz sine 39.6,
        // and a two-sample click 8.40 kHz at amplitude 0.85 — 11.40 kHz at 1.0. A CLICK HAS NO ONE ANSWER:
        // its excess over the ceiling grows with its amplitude while its width does not, and e/A is what
        // the formula reads, so the figure is meaningless unless the amplitude is given with it. Both
        // numbers are pinned in the tests; the 10.8 kHz that stood here was carried over from elsewhere
        // and belonged to neither.
        //
        // IT READS THE CREST'S SHAPE, NOT A SPECTRUM, and the difference shows on a crest that has been
        // flattened by earlier saturation: it reads low, because a flat top IS slow. For the decision this
        // instrument serves that is the right answer — a clipper would mangle it — but it is not a claim
        // about the band the energy sits in.
        //
        // AND IT READS LOW BY A KNOWN AMOUNT, because a cosine's crest is deeper than the parabola that
        // approximates it. Exactly: the true half-width solves A*cos(wD/2) = T, so wD/2 = acos(T/A), while
        // the parabola gives sqrt(2e/A). The ratio is the estimator's own bias and depends only on how far
        // the crest goes over:
        //
        //     excess    0.1 dB   0.3    0.5    1.0    2.0    3.0    6.0
        //     reads     0.999    0.997  0.995  0.991  0.982  0.975  0.955   of the true frequency
        //
        // So a 40 Hz tone one decibel over reads 39.63 Hz, measured 39.6. It is monotone, it is under a
        // percent at the depths this instrument is used at, and it is NOT corrected here: the correction
        // would need the true crest shape, which is the thing being estimated. A consumer calibrating a
        // boundary on its own corpus measures through the same bias and is unaffected; one carrying a
        // boundary in from elsewhere should know the estimate sits a little low.
        double crestHz (double sampleRate, double thresholdLinear) const noexcept
        {
            const double d = (double) lengthOs / ((double) kFactor * sampleRate);
            const double e = peak - thresholdLinear;
            if (! (d > 0.0) || ! (e > 0.0) || ! (peak > 0.0)) return 0.0;
            const double w = std::sqrt (8.0 * e / (peak * d * d));
            return std::isfinite (w) ? w / (2.0 * core::kPi) : 0.0;
        }
    };

    //--------------------------------------------------------------------------------------------------
    // WHAT prepare() ASKS THE HEAP FOR (law 11d): the budget is the allocation, and the container the
    // interpolators live in counts as much as their buffers do.
    struct Storage
    {
        bool        ok        = false;
        std::size_t runs      = 0;     // Run records
        std::size_t scratchOs = 0;     // interpolator output, oversampled floats
        std::size_t channels  = 0;

        std::uint64_t bytes() const noexcept
        {
            return (std::uint64_t) runs      * sizeof (Run)
                 + (std::uint64_t) scratchOs * sizeof (float)
                 + (std::uint64_t) channels  * ((std::uint64_t) oneOversamplerBytes()
                                              + sizeof (oversampling::PolyphaseOversampler));
        }

        static std::uint64_t oneOversamplerBytes() noexcept
        {
            oversampling::PolyphaseOversampler::Storage one {};
            if (! oversampling::PolyphaseOversampler::storageFor (kFactor, 1, kTapsPerPhase, one)) return 0;
            return one.bytes();
        }
    };

    // ok == false on exactly the arguments prepare() refuses, and then every count is zero. A NaN
    // parameter fails every one of these comparisons, which is the point.
    static Storage storageFor (double sampleRate, int maxChannels, const Params& p) noexcept
    {
        Storage st;
        if (! (sampleRate >= core::kMinSampleRate) || ! (sampleRate <= kMaxSampleRate)) return st;
        if (maxChannels < 1 || maxChannels > core::kMaxChannels) return st;
        if (! (p.thresholdDbtp >= -200.0) || ! (p.thresholdDbtp <= 200.0)) return st;
        // BOUNDED BEFORE ANY CONVERSION. `llround` saturates on arm64 and wraps on x86-64, so a refusal
        // set decided after the conversion differs between rows — half a parity contract gone.
        if (! (p.mergeMs >= 0.0) || ! (p.mergeMs <= kMaxMergeMs)) return st;
        if (p.maxRuns < 0) return st;
        double prev = 0.0;
        for (int k = 0; k < kClassEdges; ++k)
        {
            const double e = p.classEdgesMs[k];
            if (! (e > prev) || ! (e <= kMaxClassMs)) return st;   // strictly increasing, finite, bounded
            prev = e;
        }
        st.ok        = true;
        st.runs      = (std::size_t) p.maxRuns;
        st.scratchOs = (std::size_t) maxChannels * (std::size_t) (kChunk * kFactor);
        st.channels  = (std::size_t) maxChannels;
        return st;
    }

    void setParams (const Params& p) noexcept { params_ = p; }     // structural: read at the next prepare()
    const Params& params() const noexcept { return params_; }

    // Law 11b: DISARMED FIRST, so a refused prepare leaves an object that cannot be read as a measurement
    // rather than one carrying the previous run's numbers.
    [[nodiscard]] bool prepare (double sampleRate, int maxChannels) noexcept
    {
        prepared_ = false;
        finished_ = false;
        const Storage st = storageFor (sampleRate, maxChannels, params_);
        if (! st.ok) return false;

        sampleRate_  = sampleRate;
        maxChannels_ = maxChannels;
        thresholdLinear_ = core::det::pow10 (params_.thresholdDbtp / 20.0);
        if (! std::isfinite (thresholdLinear_)) return false;

        // Every decision below is on an INTEGER oversampled clock, converted once, here.
        const double osRate = sampleRate * (double) kFactor;
        mergeOs_ = (std::int64_t) std::llround (params_.mergeMs * osRate / 1000.0);
        if (mergeOs_ < 0) mergeOs_ = 0;
        for (int k = 0; k < kClassEdges; ++k)
            classEdgeOs_[k] = (std::int64_t) std::llround (params_.classEdgesMs[k] * osRate / 1000.0);

        runs_.assign (st.runs, Run {});
        scratch_.assign (st.scratchOs, 0.0f);
        os_.resize ((std::size_t) maxChannels);
        for (auto& o : os_) if (! o.prepare (kFactor, 1, kTapsPerPhase)) return false;

        prepared_ = true;
        reset();
        return true;
    }

    void reset() noexcept
    {
        finished_   = false;
        ranNc_      = 0;
        samples_    = 0;
        osMeasured_ = 0;
        osSkipped_  = 0;
        runCount_   = 0;
        aboveOsTotal_ = 0;
        totalDose_  = 0.0;
        maxExcess_  = 0.0;
        reconPeak_  = 0.0;
        samplePeak_ = 0.0;
        nonFinite_  = 0;
        firstNonFiniteAt_ = -1;
        open_       = false;
        openStart_ = openLastAbove_ = 0;
        openAbove_ = 0;
        openPeak_  = 0.0;
        openDose_  = 0.0;
        prev1_ = prev2_ = 0.0;
        haveTwo_ = false;
        for (int c = 0; c < core::kMaxChannels; ++c) aboveOsCh_[c] = 0;
        for (int k = 0; k < kClasses; ++k) { classCount_[k] = 0; classDose_[k] = 0.0; }
        for (int b = 0; b < kDurationBins; ++b) durHist_[b] = 0;
        for (int b = 0; b < kCeilingBins; ++b) ceilHist_[b] = 0;
        for (int b = 0; b < kCrestBins; ++b) { crestCount_[b] = 0; crestDose_[b] = 0.0; }
        ceilingMaxima_ = 0;
        for (auto& o : os_) o.resetChannel (0);
    }

    //--------------------------------------------------------------------------------------------------
    // READ-ONLY, allocation-free. Law 11a: a call narrower than the first is REFUSED rather than guessed
    // at — the union at an oversampled index needs every channel at that index, and a stopped channel's
    // drained tail arrives 63 oversampled samples out of phase with the others, so a "reset the missing
    // lane" rule would publish a union over two different times.
    [[nodiscard]] bool process (const float* const* in, int numChannels, int numSamples) noexcept
    {
        if (! prepared_ || finished_ || in == nullptr) return false;
        if (numChannels < 1 || numChannels > maxChannels_) return false;
        if (numSamples < 0) return false;
        if (ranNc_ == 0) ranNc_ = numChannels;
        else if (numChannels != ranNc_) return false;
        for (int c = 0; c < numChannels; ++c) if (in[c] == nullptr) return false;
        if (numSamples == 0) return true;

        for (int off = 0; off < numSamples; off += kChunk)
        {
            const int m = std::min (kChunk, numSamples - off);
            for (int c = 0; c < numChannels; ++c)
            {
                const float* src[1] { in[c] + off };
                float*       dst[1] { planeFor (c) };
                os_[(std::size_t) c].upsample (src, 1, m, dst);
            }
            // The sample peak on its own clock: the interpolator's half-sample delay means no output
            // phase reproduces the input samples, so the floor the certificate carries cannot be read
            // back out of the oversampled stream and is measured here.
            for (int i = 0; i < m; ++i)
                for (int c = 0; c < numChannels; ++c)
                {
                    const float x = in[c][off + i];
                    if (! std::isfinite (x))
                    {
                        if (nonFinite_ != ~std::uint64_t {}) ++nonFinite_;
                        if (firstNonFiniteAt_ < 0) firstNonFiniteAt_ = samples_ + (std::int64_t) i;
                        continue;
                    }
                    const double a = std::fabs ((double) x);
                    if (a > samplePeak_) samplePeak_ = a;
                }
            // THE COUNT MOVES FIRST, and that is not a tidy-up. `consumeOs` bounds itself against
            // `samples_ * kFactor`; with the increment after the call, the FIRST chunk was bounded against
            // ZERO and stopped after `kTapsPerPhase * kFactor` oversampled samples — it swallowed input
            // samples 32 to 1024 of every programme. The aggregate looked healthy (117 of 118 runs read the
            // right crest frequency), which is how a defect of this shape hides; what showed it was one run
            // whose crest frequency was twice the tone's, because a truncated crest is a faster crest.
            samples_ += m;
            consumeOs (m * kFactor, numChannels);
        }
        return true;
    }

    // Idempotent. Pushes the interpolator's own tail out and takes exactly the programme's share of it —
    // taking the WHOLE drain would measure `kTapsPerPhase` base samples of silence as programme.
    [[nodiscard]] bool finish() noexcept
    {
        if (! prepared_) return false;
        if (finished_) return true;
        if (ranNc_ > 0)
        {
            float zeros[kTapsPerPhase] {};
            for (int c = 0; c < ranNc_; ++c)
            {
                const float* src[1] { zeros };
                float*       dst[1] { planeFor (c) };
                os_[(std::size_t) c].upsample (src, 1, kTapsPerPhase, dst);
            }
            const std::int64_t room = samples_ * (std::int64_t) kFactor - osMeasured_;
            if (room > 0) consumeOs ((int) std::min<std::int64_t> (room, kTapsPerPhase * kFactor), ranNc_);
        }
        if (open_) closeRun (true);
        finished_ = true;
        return true;
    }

    //--------------------------------------------------------------------------------------------------
    // --- the report ---
    bool  isPrepared() const noexcept { return prepared_; }
    bool  isFinished() const noexcept { return finished_; }
    Reason reason() const noexcept
    {
        if (! finished_) return Reason::NotFinished;
        if (osMeasured_ <= 0) return Reason::Empty;
        if (nonFinite_ > 0) return Reason::NonFiniteInput;
        return Reason::Ok;
    }
    bool valid() const noexcept { return reason() == Reason::Ok; }

    double sampleRate() const noexcept { return sampleRate_; }
    int    channels()   const noexcept { return ranNc_; }
    double thresholdLinear() const noexcept { return thresholdLinear_; }
    std::int64_t samplesProcessed() const noexcept { return samples_; }
    // Oversampled samples that belonged to the programme. `samples_ * kFactor` exactly, once finished.
    std::int64_t measuredOs() const noexcept { return osMeasured_; }

    // THE TWO PEAKS, SEPARATELY, because the certificate is the larger of them and an instrument that
    // published only one would disagree with it without saying where.
    double reconstructedPeak() const noexcept { return reconPeak_; }
    double samplePeakLinear() const noexcept { return samplePeak_; }
    double truePeakLinear()   const noexcept { return std::max (reconPeak_, samplePeak_); }

    std::int64_t runCount()       const noexcept { return runCount_; }
    std::int64_t storedRunCount() const noexcept { return std::min<std::int64_t> (runCount_, (std::int64_t) runs_.size()); }
    bool         runsComplete()   const noexcept { return runCount_ <= (std::int64_t) runs_.size(); }
    Run run (std::int64_t i) const noexcept
    {
        return i >= 0 && i < storedRunCount() ? runs_[(std::size_t) i] : Run {};
    }

    // MERGE-FREE, both of them: accumulated over raw above-samples, so a consumer sweeping `mergeMs`
    // keeps two fixed points while the count and the percentile move.
    std::int64_t aboveOs() const noexcept { return aboveOsTotal_; }
    std::int64_t aboveOs (int c) const noexcept
    {
        return c >= 0 && c < core::kMaxChannels ? aboveOsCh_[c] : 0;
    }
    double occupancy() const noexcept
    {
        return osMeasured_ > 0 ? (double) aboveOsTotal_ / (double) osMeasured_ : 0.0;
    }
    double runsPerMinute() const noexcept
    {
        const double sec = osMeasured_ > 0 ? (double) osMeasured_ / ((double) kFactor * sampleRate_) : 0.0;
        return sec > 0.0 ? (double) runCount_ * 60.0 / sec : 0.0;
    }

    double totalDose() const noexcept { return totalDose_; }
    double maxExcess() const noexcept { return maxExcess_; }
    std::int64_t classCount (int k) const noexcept { return k >= 0 && k < kClasses ? classCount_[k] : 0; }
    double       classDose  (int k) const noexcept { return k >= 0 && k < kClasses ? classDose_[k]  : 0.0; }

    // NEAREST-RANK, FROM THE HISTOGRAM AND NEVER FROM THE LIST, so it stays exact-to-the-bin after the
    // list is full. r = ceil(0.9 N) in integers. For N <= 10 that is N itself, so the answer IS the
    // longest run — said here because "the 90th percentile of three runs" is a sentence that reads as
    // more than it is. 0.0 at N = 0, with runCount() beside it to tell that from "short runs".
    double p90Ms() const noexcept
    {
        const std::int64_t n = runCount_;
        if (n <= 0) return 0.0;
        const std::int64_t rank = (9 * n + 9) / 10;
        std::int64_t seen = 0;
        for (int b = 0; b < kDurationBins; ++b)
        {
            seen += durHist_[b];
            if (seen >= rank) return (double) (b + 1) * kDurationBinMs;
        }
        return (double) kDurationBins * kDurationBinMs;
    }
    bool p90Saturated() const noexcept
    {
        std::int64_t seen = 0;
        for (int b = 0; b < kDurationBins - 1; ++b) seen += durHist_[b];
        return runCount_ > 0 && seen < runCount_;
    }
    std::int64_t durationBin (int b) const noexcept { return b >= 0 && b < kDurationBins ? durHist_[b] : 0; }

    // THE CREST'S FREQUENCY, aggregated. Third-octaves from 20 Hz, by count and by dose, so a consumer
    // sweeps "is this bass" on its own corpus instead of receiving a boundary baked in at 120 Hz. The
    // histogram is an accumulator and survives the run list's exhaustion.
    std::int64_t crestBinCount (int b) const noexcept { return b >= 0 && b < kCrestBins ? crestCount_[b] : 0; }
    double       crestBinDose  (int b) const noexcept { return b >= 0 && b < kCrestBins ? crestDose_[b]  : 0.0; }
    static double crestBinLowHz (int b) noexcept
    {
        return b >= 0 && b < kCrestBins ? kCrestBaseHz * core::det::exp2 ((double) b / 3.0) : 0.0;
    }
    // Dose whose crest frequency is under `hz`, over the total — the number `lowShare` was meant to be,
    // built out of shape rather than out of a filter. Exact to the third-octave boundary below `hz`.
    double doseShareBelow (double hz) const noexcept
    {
        if (! (hz > 0.0) || ! (totalDose_ > 0.0)) return 0.0;
        double d = 0.0;
        for (int b = 0; b < kCrestBins; ++b)
            if (crestBinLowHz (b + 1 < kCrestBins ? b + 1 : b) <= hz || crestBinLowHz (b) < hz) d += crestDose_[b];
        return d / totalDose_;
    }

    // HOW TIGHTLY THE PEAKS SIT UNDER THE CEILING — the tell for a source that was true-peak LIMITED
    // rather than clipped: its maxima cluster at the ceiling with no flat tops for a clipping test to
    // find. The denominator is every local maximum of the reconstruction, which is dominated by the
    // musical waveform, so a dynamic programme reads low for reasons that are not about peak control;
    // `ceilingDensityAbove` takes the same histogram with a floor under it.
    std::int64_t ceilingMaxima() const noexcept { return ceilingMaxima_; }
    // THE BINS RUN DOWNWARD IN LEVEL — `b = floor((kCeilingTopDb - dB) / kCeilingBinDb)` — so `top` is the
    // LOUDEST populated bin and every quieter maximum sits at a LARGER index. The first version of this
    // walked `b = 0 .. top`, which is the range ABOVE the loudest maximum and is empty by definition, so
    // numerator and denominator were both `ceilHist_[top]` and the answer was EXACTLY 1 on every
    // programme. A consumer measured 1 on 102 files and said so; nothing here had a test.
    //
    // -1.0, NOT 0.0, for a request this cannot answer — a bad argument, or a call before a measurement.
    // A density lives in [0, 1] and 0.0 is a legitimate reading ("no maximum sits near the loudest"), so
    // the two must not share a value. Same rule, and the same reason, as LowEnd::sideFractionBelow.
    double ceilingDensity (double withinDb = 0.2) const noexcept
    {
        return ceilingDensityAbove ((double) kCeilingSpanDb, withinDb);
    }
    double ceilingDensityAbove (double minusDb, double withinDb = 0.2) const noexcept
    {
        if (! (withinDb > 0.0) || ! (minusDb > 0.0)) return -1.0;
        if (ceilingMaxima_ <= 0) return -1.0;
        const int top = ceilingTopBin();
        if (top < 0) return -1.0;
        // BOTH WIDTHS ARE CLAMPED IN DOUBLE BEFORE THEY BECOME AN int. The harness asks for 1e9 dB to mean
        // "every maximum", and 1e9 / 0.05 is 2e10 — not representable in an int, so the conversion alone
        // would be undefined. The histogram is kCeilingBins wide; nothing past it exists to count.
        const double binsPerDb = 1.0 / kCeilingBinDb;
        const int inner = (int) std::ceil  (std::fmin (withinDb * binsPerDb, (double) kCeilingBins));
        const int span  = (int) std::floor (std::fmin (minusDb  * binsPerDb, (double) kCeilingBins));
        std::int64_t num = 0, den = 0;
        for (int b = top; b < kCeilingBins && b - top <= span; ++b)
        {
            den += ceilHist_[b];
            if (b - top < inner) num += ceilHist_[b];
        }
        return den > 0 ? (double) num / (double) den : -1.0;
    }
    std::int64_t ceilingBin (int b) const noexcept { return b >= 0 && b < kCeilingBins ? ceilHist_[b] : 0; }
    static double ceilingBinTopDb (int b) noexcept { return kCeilingTopDb - (double) b * kCeilingBinDb; }

    std::uint64_t nonFiniteSamples() const noexcept { return nonFinite_; }
    std::int64_t  firstNonFiniteAt() const noexcept { return firstNonFiniteAt_; }

private:
    float* planeFor (int c) noexcept
    {
        return scratch_.data() + (std::size_t) c * (std::size_t) (kChunk * kFactor);
    }

    // One pass over `n` oversampled samples of the scratch planes. The clock is the oversampled index,
    // which is a function of the samples fed and not of how they were sliced — law 8a by construction.
    void consumeOs (int n, int numChannels) noexcept
    {
        for (int k = 0; k < n; ++k)
        {
            double mag = 0.0;
            for (int c = 0; c < numChannels; ++c)
            {
                const float v = planeFor (c)[k];
                const double a = std::isfinite (v) ? std::fabs ((double) v) : 0.0;
                if (a > mag) mag = a;
                if (a > thresholdLinear_)
                {
                    ++aboveOsCh_[c];
                }
            }
            // THE FIRST kDelayOs OUTPUTS ARE THE FILTER'S RAMP-UP, not programme time. Skipping them is
            // what puts `startOs` on the programme's own grid.
            if (osSkipped_ < kDelayOs) { ++osSkipped_; continue; }
            if (osMeasured_ >= samples_ * (std::int64_t) kFactor + (std::int64_t) kTapsPerPhase * kFactor) return;

            if (mag > reconPeak_) reconPeak_ = mag;
            noteLocalMax (mag);

            const double e = mag - thresholdLinear_;
            const bool above = e > 0.0;
            if (above)
            {
                ++aboveOsTotal_;
                totalDose_ += e;
                if (e > maxExcess_) maxExcess_ = e;
                if (! open_)
                {
                    open_ = true;
                    openStart_ = osMeasured_;
                    openAbove_ = 0;
                    openPeak_  = 0.0;
                    openDose_  = 0.0;
                }
                openLastAbove_ = osMeasured_;
                ++openAbove_;
                openDose_ += e;
                if (mag > openPeak_) openPeak_ = mag;
            }
            else if (open_ && osMeasured_ - openLastAbove_ > mergeOs_)
            {
                closeRun (false);
            }
            ++osMeasured_;
        }
    }

    // A local maximum of the reconstruction, from a three-sample window. It lags the stream by one sample,
    // which costs nothing here because the histogram is an aggregate with no coordinate.
    //
    // THE RISE IS STRICT AND THE FALL IS NOT, so a PLATEAU counts ONCE — at its first sample. With `>=` on
    // both sides every sample of a flat top is a maximum, and a flat top is exactly what this instrument
    // is pointed at: a clipped programme's L-sample plateau then contributed L counts at the ceiling and
    // the density read the DURATION of the flatness as its crowding. It is not a corner case either — the
    // reconstruction's group delay is 63.5 oversampled samples, a half sample, so the crest of an isolated
    // impulse falls BETWEEN two samples that are then bit-identical by symmetry, and even a lone spike
    // counted twice.
    void noteLocalMax (double mag) noexcept
    {
        if (haveTwo_ && prev1_ > prev2_ && prev1_ >= mag && prev1_ > 0.0)
        {
            const double db = 20.0 * core::det::log10 (prev1_);
            int b = (int) std::floor ((kCeilingTopDb - db) / kCeilingBinDb);
            if (b < 0) b = 0;
            if (b >= kCeilingBins) b = kCeilingBins - 1;
            ++ceilHist_[b];
            ++ceilingMaxima_;
        }
        prev2_ = prev1_;
        prev1_ = mag;
        haveTwo_ = true;
    }

    int ceilingTopBin() const noexcept
    {
        for (int b = 0; b < kCeilingBins; ++b) if (ceilHist_[b] > 0) return b;
        return -1;
    }

    void closeRun (bool byFinish) noexcept
    {
        const std::int64_t len = openLastAbove_ - openStart_ + 1;
        if (runCount_ < (std::int64_t) runs_.size())
        {
            Run& r = runs_[(std::size_t) runCount_];
            r.startOs  = openStart_;
            r.lengthOs = len;
            r.aboveOs  = openAbove_;
            r.peak     = openPeak_;
            r.dose     = openDose_;
            r.closedByFinish = byFinish;
        }
        ++runCount_;

        // EVERY AGGREGATE IS AN ACCUMULATOR AND NEVER READS THE LIST, which is what keeps the percentile,
        // the classes and the crest histogram exact after the list is full.
        const double ms = (double) len * 1000.0 / ((double) kFactor * sampleRate_);
        int db = (int) std::floor (ms / kDurationBinMs);
        if (db < 0) db = 0;
        if (db >= kDurationBins) db = kDurationBins - 1;
        ++durHist_[db];

        int cls = kClasses - 1;
        for (int k = 0; k < kClassEdges; ++k) if (len <= classEdgeOs_[k]) { cls = k; break; }
        ++classCount_[cls];
        classDose_[cls] += openDose_;

        Run tmp; tmp.lengthOs = len; tmp.peak = openPeak_;
        const double hz = tmp.crestHz (sampleRate_, thresholdLinear_);
        if (hz > 0.0)
        {
            int cb = (int) std::floor (3.0 * core::det::log2 (hz / kCrestBaseHz));
            if (cb < 0) cb = 0;
            if (cb >= kCrestBins) cb = kCrestBins - 1;
            ++crestCount_[cb];
            crestDose_[cb] += openDose_;
        }
        open_ = false;
    }

    Params params_ {};
    bool   prepared_ = false, finished_ = false;
    double sampleRate_ = 48000.0, thresholdLinear_ = 0.0;
    int    maxChannels_ = 0, ranNc_ = 0;
    std::int64_t mergeOs_ = 0, classEdgeOs_[kClassEdges] {};

    std::vector<Run>   runs_;
    std::vector<float> scratch_;
    std::vector<oversampling::PolyphaseOversampler> os_;

    std::int64_t samples_ = 0, osMeasured_ = 0, osSkipped_ = 0, runCount_ = 0, aboveOsTotal_ = 0;
    std::int64_t aboveOsCh_[core::kMaxChannels] {};
    double totalDose_ = 0.0, maxExcess_ = 0.0, reconPeak_ = 0.0, samplePeak_ = 0.0;
    std::uint64_t nonFinite_ = 0;
    std::int64_t  firstNonFiniteAt_ = -1;

    bool open_ = false;
    std::int64_t openStart_ = 0, openLastAbove_ = 0, openAbove_ = 0;
    double openPeak_ = 0.0, openDose_ = 0.0;

    double prev1_ = 0.0, prev2_ = 0.0;
    bool   haveTwo_ = false;

    std::int64_t classCount_[kClasses] {};
    double       classDose_[kClasses] {};
    std::int64_t durHist_[kDurationBins] {};
    std::int64_t ceilHist_[kCeilingBins] {};
    std::int64_t ceilingMaxima_ = 0;
    std::int64_t crestCount_[kCrestBins] {};
    double       crestDose_[kCrestBins] {};
};

} // namespace felitronics::analysis

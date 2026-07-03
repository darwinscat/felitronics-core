// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/core/Config.h>
#include <felitronics/core/Math.h>
#include <felitronics/core/DelayLine.h>
#include <felitronics/core/FlushToZero.h>
#include <felitronics/oversampling/PolyphaseOversampler.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

namespace felitronics::limiter
{

namespace detail
{
    // O(1)-amortized sliding-window maximum (monotonic-decreasing deque). Front = max over the last W.
    // prepare() fixes the CAPACITY; setWindow() picks the effective window ≤ capacity, RT-safe (no alloc),
    // so the limiter can retune the lookahead per setParams() without a re-prepare. Expiry runs BEFORE the
    // insert — insert-first can overwrite the head (the current max) once a strictly-decreasing run fills
    // the deque, under-reading the max exactly when a decaying peak is still inside the lookahead. Indices
    // are int64 (a `long` is 32-bit on Windows → wraps after ~3 h of oversampled pushes).
    class SlidingMax
    {
    public:
        void prepare (int maxWindow)
        {
            cap = maxWindow < 1 ? 1 : maxWindow;
            v.assign ((std::size_t) cap, 0.0f); ix.assign ((std::size_t) cap, 0);
            W = cap;
            reset();
        }

        void reset() noexcept { head = tail = count = 0; n = 0; }

        void setWindow (int w) noexcept { W = w < 1 ? 1 : (w > cap ? cap : w); }   // stale entries expire on the next pushes

        float push (float x) noexcept
        {
            while (count > 0 && ix[(std::size_t) head] <= n - (std::int64_t) W) { head = (head + 1) % cap; --count; }
            while (count > 0) { const int b = (tail - 1 + cap) % cap; if (v[(std::size_t) b] <= x) { tail = b; --count; } else break; }
            v[(std::size_t) tail] = x; ix[(std::size_t) tail] = n; tail = (tail + 1) % cap; ++count;
            ++n;
            return v[(std::size_t) head];
        }

    private:
        int cap = 1, W = 1; std::vector<float> v; std::vector<std::int64_t> ix; int head = 0, tail = 0, count = 0; std::int64_t n = 0;
    };
}

struct TruePeakLimiterParams
{
    double ceilingDbTp     = -1.0;     // output true-peak ceiling (dBTP)
    double releaseMs       = 100.0;
    double lookaheadMs     = 1.0;
    int    oversampleFactor = 4;       // fixed at prepare(); changing it needs a re-prepare
};

//==============================================================================
// felitronics::limiter::TruePeakLimiter — a brickwall true-peak limiter that GUARANTEES the output
// true peak ≤ ceiling (within the down-sampler's tiny ripple). It oversamples, limits at the
// oversampled rate (so inter-sample peaks are real samples), then downsamples (Option B — the only
// path that truly bounds true-peak; applying a fast-moving gain at baseband would itself create ISP).
//
// Gain law: a sliding-window MAXIMUM of the channel-linked oversampled peak over the lookahead window
// → required gain = ceiling / slidingMax (so EVERY sample in the window, including the one being output
// `lookahead` samples behind, is ≤ ceiling) → instant attack + rate-limited (release) recovery. This is
// the structure that makes the bound provable, not heuristic.
//
// RT-safe: prepare() allocates; process() does no alloc/lock/throw. Reported latency = oversampler
// round-trip + lookahead (baseband samples). One linked gain for all channels (no image shift).
class TruePeakLimiter
{
public:
    void prepare (double sampleRate, int maxBlock, int maxChannels, int tapsPerPhase = 32)
    {
        prepared_ = false;                                     // any early return below leaves it unprepared
        if (sampleRate <= 0.0 || maxBlock < 1) return;         // invalid stream config
        fs    = sampleRate;
        maxCh = maxChannels < 1 ? 1 : (maxChannels > core::kMaxChannels ? core::kMaxChannels : maxChannels);
        maxBlock_ = maxBlock;
        tpp   = tapsPerPhase;
        F     = params.oversampleFactor < 2 ? 2 : params.oversampleFactor;
        if (! os.prepare (F, maxCh, tpp)) return;              // oversampler rejected (e.g. tapsPerPhase<4) → stay unprepared

        osBuf.assign ((std::size_t) maxCh, std::vector<float> ((std::size_t) maxBlock * F, 0.0f));
        osPtrs.assign ((std::size_t) maxCh, nullptr);

        const int maxLookOS = (int) std::ceil (kMaxLookaheadMs * 0.001 * fs) * F;
        osDelays.assign ((std::size_t) maxCh, core::DelayLine {});
        for (auto& d : osDelays) d.prepare (maxLookOS);
        slide.prepare (maxLookOS + 1);

        apply (params);
        reset();
        prepared_ = true;                                      // fully built — process() may now run
    }

    void reset() noexcept
    {
        for (auto& b : osBuf) std::fill (b.begin(), b.end(), 0.0f);
        for (auto& d : osDelays) d.reset();
        slide.reset();
        grDb = 0.0f;
    }

    void setParams (const TruePeakLimiterParams& p) noexcept { params = p; apply (p); }   // factor change needs re-prepare
    int    latencySamples()  const noexcept { return os.latencySamples() + lookBaseband; }
    double gainReductionDb() const noexcept { return grDb; }

    // Audio thread, in place (baseband). RT-safe.
    void process (float* const* channels, int numChannels, int numSamples) noexcept
    {
        const int nc = numChannels < maxCh ? numChannels : maxCh;
        if (! prepared_ || nc <= 0 || numSamples <= 0 || numSamples > maxBlock_) return;   // unprepared / oversized → no OOB in osBuf
        const int osN = numSamples * F;

        for (int c = 0; c < nc; ++c) osPtrs[(std::size_t) c] = osBuf[(std::size_t) c].data();
        os.upsample (channels, nc, numSamples, osPtrs.data());

        for (int i = 0; i < osN; ++i)
        {
            float linkedPeak = 0.0f;
            for (int c = 0; c < nc; ++c) { const float a = std::fabs (osBuf[(std::size_t) c][(std::size_t) i]); if (a > linkedPeak) linkedPeak = a; }

            const float  smax    = slide.push (linkedPeak);
            const double smaxDb  = core::gainToDb (smax);
            double rawRedDb = ceilingDb - smaxDb;
            if (rawRedDb > 0.0) rawRedDb = 0.0;

            grDb = std::min ((float) rawRedDb, grDb * relCoef);   // instant attack, exponential release toward 0 dB
            const float gain = (float) core::dbToGain ((double) grDb);

            for (int c = 0; c < nc; ++c)
            {
                const float x = osBuf[(std::size_t) c][(std::size_t) i];
                osBuf[(std::size_t) c][(std::size_t) i] = osDelays[(std::size_t) c].process (x) * gain;
            }
        }

        os.downsample ((const float* const*) osPtrs.data(), nc, numSamples, channels);
        core::flushDenormal (grDb);
    }

private:
    static constexpr double kMaxLookaheadMs = 20.0;

    void apply (const TruePeakLimiterParams& p) noexcept
    {
        // Non-finite params fall back to the defaults (house rule) — a NaN ceiling/lookahead would poison
        // grDb / feed lround() UB.
        ceilingDb    = std::isfinite (p.ceilingDbTp) ? p.ceilingDbTp : -1.0;
        const double lookMs = std::isfinite (p.lookaheadMs) ? p.lookaheadMs : 0.0;
        lookBaseband = (int) std::lround (lookMs * 0.001 * fs);
        if (lookBaseband < 0) lookBaseband = 0;
        const int maxLookBb = (int) std::ceil (kMaxLookaheadMs * 0.001 * fs);
        if (lookBaseband > maxLookBb) lookBaseband = maxLookBb;
        const int lookOS = lookBaseband * F;
        for (auto& d : osDelays) d.setDelay (lookOS);
        // The max window must equal the ACTUAL lookahead (+1 for the sample being output). The capacity-wide
        // window held the gain down for the full 20 ms after every peak — a release that ignored the knob.
        slide.setWindow (lookOS + 1);

        const double relMs = std::isfinite (p.releaseMs) ? p.releaseMs : 100.0;
        const double t = relMs * 0.001 * fs * F;         // release time constant at the OS rate
        relCoef = (t <= 0.0) ? 0.0f : (float) std::exp (-1.0 / t);
    }

    double fs = 48000.0;
    int maxCh = 0, tpp = 32, F = 4, maxBlock_ = 0;
    bool prepared_ = false;                     // true only after a fully-successful prepare()
    TruePeakLimiterParams params;

    oversampling::PolyphaseOversampler os;
    std::vector<std::vector<float>>    osBuf;
    std::vector<float*>                osPtrs;
    std::vector<core::DelayLine>       osDelays;
    detail::SlidingMax                 slide;

    double ceilingDb = -1.0;
    int    lookBaseband = 0;
    float  relCoef = 0.0f, grDb = 0.0f;
};

} // namespace felitronics::limiter

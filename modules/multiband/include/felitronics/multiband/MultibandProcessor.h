// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/eq/MultibandSplitter.h>
#include <felitronics/core/DelayLine.h>
#include <felitronics/core/Config.h>

#include <algorithm>
#include <array>
#include <vector>

namespace felitronics::multiband
{

//==============================================================================
// felitronics::multiband::MultibandProcessor<BandProcessor, MaxBands> — split a signal into N LR4 bands
// (eq::MultibandSplitter), run ONE independent BandProcessor per band, and recombine. The generic engine
// behind a multiband compressor / saturator / width.
//
// BandProcessor concept: process(float* const* io, int nch, int n) [in place] · int latencySamples() · reset().
// prepare() is given a callable that prepares each band instance (their prepare signatures differ).
//
// Per-band bypass + solo. A global parallel dry/wet whose DRY is the band-sum ALLPASS RECONSTRUCTION (delayed
// to match the wet's latency) — NOT the raw input, which would comb against the allpass-phase wet. Bands
// are aligned to the max band latency via per-band
// DelayLines so a band with lookahead doesn't desync the sum. RT-safe: prepare() allocates, process() is
// alloc/lock/throw-free. n ≤ maxBlock.
template <class BandProcessor, int MaxBands = 4>
class MultibandProcessor
{
    static_assert (MaxBands >= 2 && MaxBands <= 8, "MultibandProcessor supports 2..8 bands");

public:
    // `maxAlignSamples` bounds the per-band latency-alignment delay (≥ the largest band lookahead in samples).
    template <class PrepareBand>
    bool prepare (double sampleRate, int maxBlock, int maxChannels, int maxAlignSamples, PrepareBand&& prepareBand)
    {
        if (sampleRate <= 0.0 || maxBlock < 1) return false;
        fs_ = sampleRate; maxBlock_ = maxBlock;
        channels_ = std::clamp (maxChannels, 1, core::kMaxChannels);
        if (maxAlignSamples < 0) maxAlignSamples = 0;

        splitter_.prepare (sampleRate, channels_);
        bandBuf_.assign ((std::size_t) MaxBands * (std::size_t) channels_ * (std::size_t) maxBlock_, 0.0f);
        dryBuf_.assign  ((std::size_t) channels_ * (std::size_t) maxBlock_, 0.0f);
        for (int b = 0; b < MaxBands; ++b) { bandPtrs_[(std::size_t) b].assign ((std::size_t) channels_, nullptr); prepareBand (proc_[(std::size_t) b]); }

        align_.assign ((std::size_t) MaxBands * (std::size_t) channels_, core::DelayLine {});
        for (auto& d : align_) d.prepare (maxAlignSamples);
        dryDelay_.assign ((std::size_t) channels_, core::DelayLine {});
        for (auto& d : dryDelay_) d.prepare (maxAlignSamples);

        refreshLatency(); reset();
        return true;
    }

    void reset() noexcept
    {
        splitter_.reset();
        for (auto& p : proc_) p.reset();
        for (auto& d : align_) d.reset();
        for (auto& d : dryDelay_) d.reset();
    }

    bool setNumBands (int n) noexcept { if (! splitter_.setNumBands (n)) return false; refreshLatency(); reset(); return true; }
    void setCrossovers (const float* hz, int count) noexcept { splitter_.setCrossovers (hz, count); }
    int  numBands() const noexcept { return splitter_.numBands(); }

    BandProcessor&       band (int b)       noexcept { return proc_[(std::size_t) std::clamp (b, 0, MaxBands - 1)]; }
    const BandProcessor& band (int b) const noexcept { return proc_[(std::size_t) std::clamp (b, 0, MaxBands - 1)]; }

    // Set a band's params + re-derive latency (a lookahead change shifts alignment). OOB index → ignored.
    template <class Params>
    void setBandParams (int b, const Params& p) noexcept { if ((unsigned) b < (unsigned) MaxBands) { proc_[(std::size_t) b].setParams (p); refreshLatency(); } }

    void setBandBypass (int b, bool x) noexcept { if ((unsigned) b < (unsigned) MaxBands) { bypass_[(std::size_t) b] = x; refreshLatency(); } }
    void setBandSolo   (int b, bool x) noexcept { if ((unsigned) b < (unsigned) MaxBands) solo_[(std::size_t) b] = x; }
    void setMix        (float m) noexcept       { mix_ = std::clamp (m, 0.0f, 1.0f); }
    void refreshLatency() noexcept { recomputeLatency(); }

    int  latencySamples() const noexcept { return latency_; }

    void process (float* const* io, int numChannels, int n) noexcept
    {
        const int nc = std::min (numChannels, channels_);
        const int nb = splitter_.numBands();
        if (nc <= 0 || n <= 0 || n > maxBlock_) return;

        // 1) split into per-band planar buffers; capture the allpass-reconstructed dry for the parallel mix
        float tmp[(std::size_t) MaxBands] {};
        for (int c = 0; c < nc; ++c)
        {
            float* dry = dryData (c);
            for (int i = 0; i < n; ++i)
            {
                splitter_.splitSample (c, io[c][i], tmp);
                dry[i] = splitter_.sumSample (tmp);
                for (int b = 0; b < nb; ++b) bandData (b, c)[i] = tmp[b];
            }
        }

        // 2) process each band in place (bypassed bands keep their split signal)
        for (int b = 0; b < nb; ++b)
        {
            if (bypass_[(std::size_t) b]) continue;
            for (int c = 0; c < nc; ++c) bandPtrs_[(std::size_t) b][(std::size_t) c] = bandData (b, c);
            proc_[(std::size_t) b].process (bandPtrs_[(std::size_t) b].data(), nc, n);
        }

        // 3) align every band to the max band latency
        for (int b = 0; b < nb; ++b)
            for (int c = 0; c < nc; ++c)
            {
                core::DelayLine& d = align_[(std::size_t) (b * channels_ + c)];
                float* x = bandData (b, c);
                for (int i = 0; i < n; ++i) x[i] = d.process (x[i]);
            }

        // 4) sum (honouring solo) + global dry/wet against the latency-matched allpass dry
        bool anySolo = false; for (int b = 0; b < nb; ++b) anySolo |= solo_[(std::size_t) b];
        for (int c = 0; c < nc; ++c)
        {
            float* dry = dryData (c);
            for (int i = 0; i < n; ++i)
            {
                for (int b = 0; b < nb; ++b) tmp[b] = (! anySolo || solo_[(std::size_t) b]) ? bandData (b, c)[i] : 0.0f;
                const float wet = splitter_.sumSample (tmp);
                const float d   = dryDelay_[(std::size_t) c].process (dry[i]);   // always advance the delay
                // Solo is a monitor: hear ONLY the soloed bands (full), bypassing the parallel dry — else the
                // full-band dry would leak the non-soloed bands back in at mix < 1.
                io[c][i] = anySolo ? wet : (d + mix_ * (wet - d));
            }
        }
    }

private:
    void recomputeLatency() noexcept
    {
        const int nb = splitter_.numBands();
        latency_ = 0;
        for (int b = 0; b < nb; ++b) latency_ = std::max (latency_, bypass_[(std::size_t) b] ? 0 : proc_[(std::size_t) b].latencySamples());
        for (int b = 0; b < nb; ++b)
        {
            const int bl = bypass_[(std::size_t) b] ? 0 : proc_[(std::size_t) b].latencySamples();
            const int extra = std::max (0, latency_ - bl);
            for (int c = 0; c < channels_; ++c) align_[(std::size_t) (b * channels_ + c)].setDelay (extra);
        }
        for (auto& d : dryDelay_) d.setDelay (latency_);
    }

    float* bandData (int b, int c) noexcept { return bandBuf_.data() + ((std::size_t) b * (std::size_t) channels_ + (std::size_t) c) * (std::size_t) maxBlock_; }
    float* dryData  (int c)        noexcept { return dryBuf_.data()  + (std::size_t) c * (std::size_t) maxBlock_; }

    double fs_ = 48000.0;
    int maxBlock_ = 0, channels_ = 0, latency_ = 0;
    float mix_ = 1.0f;

    eq::MultibandSplitter<MaxBands> splitter_;
    std::array<BandProcessor, (std::size_t) MaxBands> proc_ {};
    std::array<bool, (std::size_t) MaxBands> bypass_ {}, solo_ {};
    std::array<std::vector<float*>, (std::size_t) MaxBands> bandPtrs_;
    std::vector<float> bandBuf_, dryBuf_;
    std::vector<core::DelayLine> align_, dryDelay_;
};

} // namespace felitronics::multiband

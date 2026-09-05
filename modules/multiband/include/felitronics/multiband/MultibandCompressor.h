// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/multiband/MultibandProcessor.h>
#include <felitronics/dynamics/Compressor.h>

#include <felitronics/core/Config.h>

#include <algorithm>
#include <cmath>

namespace felitronics::multiband
{

//==============================================================================
// felitronics::multiband::MultibandCompressor — a per-band broadband dynamics::Compressor on the LR4
// splitter (the cornerstone mastering tool). Per-band DETECTION (each band compresses its own level); a
// linked/broadband sidechain is a future feature (the Compressor now HAS an external detector input,
// but routing one key through the LR4 split is a multiband decision this wrapper has not taken). NO
// −6 dB crossover threshold compensation — thresholds are per-band-relative (the detector sees the actual
// band signal), the predictable convention. Per-band params / bypass / solo + a global parallel dry/wet.
// RT-safe (prepare() allocates, process() is alloc-free). Reported latency = the max per-band lookahead.
template <int MaxBands = 4>
class MultibandCompressor
{
public:
    [[nodiscard]] bool prepare (double sampleRate, int maxBlock, int maxChannels, double maxLookaheadMs = 50.0)
    {
        if (! (sampleRate > 0.0) || ! std::isfinite (sampleRate)) return false;
        // Finite is not enough: `maxLookaheadMs = 1e300` is finite, and the `(int) std::ceil(...)`
        // two lines down is then an out-of-range floating-to-integer conversion — undefined, and
        // reached BEFORE any band's prepare() could refuse it. Validate against the band's own bound.
        if (! std::isfinite (maxLookaheadMs) || maxLookaheadMs < 0.0
            || maxLookaheadMs > dynamics::Compressor::kMaxLookaheadMs) return false;
        if (sampleRate > dynamics::Compressor::kMaxSampleRate) return false;
        // Clamp ONCE and hand the same number to both layers. MultibandProcessor clamps its own copy to
        // core::kMaxChannels, so passing the raw count to the band lambda used to prepare each
        // Compressor for a wider layout than the buffers around it — invisible while Compressor::prepare
        // returned void, a hard prepare() failure now that it can refuse.
        const int mc = std::clamp (maxChannels, 1, core::kMaxChannels);
        const int maxAlign = (int) std::ceil (maxLookaheadMs * 0.001 * sampleRate) + 1;
        return mb_.prepare (sampleRate, maxBlock, mc, maxAlign,
                            [&] (dynamics::Compressor& c) { return c.prepare (sampleRate, maxBlock, mc, maxLookaheadMs); });
    }

    void reset() noexcept { mb_.reset(); }

    bool setNumBands (int n) noexcept { return mb_.setNumBands (n); }
    void setCrossovers (const float* hz, int count) noexcept { mb_.setCrossovers (hz, count); }
    int  numBands() const noexcept { return mb_.numBands(); }

    void setBandParams (int b, const dynamics::CompressorParams& p) noexcept { mb_.setBandParams (b, p); }
    void setBandBypass (int b, bool x) noexcept { mb_.setBandBypass (b, x); }
    void setBandSolo   (int b, bool x) noexcept { mb_.setBandSolo (b, x); }
    void setMix        (float m) noexcept       { mb_.setMix (m); }

    int    latencySamples() const noexcept { return mb_.latencySamples(); }
    double bandGainReductionDb (int b) noexcept { return mb_.band (b).gainReductionDb(); }

    void process (float* const* io, int numChannels, int n) noexcept { mb_.process (io, numChannels, n); }

private:
    MultibandProcessor<dynamics::Compressor, MaxBands> mb_;
};

} // namespace felitronics::multiband

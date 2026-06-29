// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/multiband/MultibandProcessor.h>
#include <felitronics/dynamics/Compressor.h>

#include <cmath>

namespace felitronics::multiband
{

//==============================================================================
// felitronics::multiband::MultibandCompressor — a per-band broadband dynamics::Compressor on the LR4
// splitter (the cornerstone mastering tool). Per-band DETECTION (each band compresses its own level); a
// linked/broadband sidechain is a future feature (the Compressor has no external detector input yet). NO
// −6 dB crossover threshold compensation — thresholds are per-band-relative (the detector sees the actual
// band signal), the predictable convention. Per-band params / bypass / solo + a global parallel dry/wet.
// RT-safe (prepare() allocates, process() is alloc-free). Reported latency = the max per-band lookahead.
template <int MaxBands = 4>
class MultibandCompressor
{
public:
    bool prepare (double sampleRate, int maxBlock, int maxChannels, double maxLookaheadMs = 50.0)
    {
        const int maxAlign = (int) std::ceil (maxLookaheadMs * 0.001 * sampleRate) + 1;
        return mb_.prepare (sampleRate, maxBlock, maxChannels, maxAlign,
                            [&] (dynamics::Compressor& c) { c.prepare (sampleRate, maxBlock, maxChannels, maxLookaheadMs); });
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

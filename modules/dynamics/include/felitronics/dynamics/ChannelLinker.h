// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <cmath>

namespace felitronics::dynamics
{

// How a multichannel detector collapses to ONE level so every channel gets the SAME gain (a stereo
// compressor must not gain channels independently — that shifts the image).
enum class LinkMode
{
    Max,        // max|ch| — image-preserving, reacts to the loudest channel (limiter-ish, firmer)
    MeanPower   // sqrt(mean(ch^2)) — gentler, channel-count-invariant, pumps less on one-sided transients
};

// One linked detector level (linear amplitude) from the channel frame at `sampleIndex`. Stateless.
inline float linkAmplitude (LinkMode mode, const float* const* channels, int numChannels, int sampleIndex) noexcept
{
    if (numChannels <= 1) return std::fabs (channels[0][sampleIndex]);

    if (mode == LinkMode::Max)
    {
        float mx = 0.0f;
        for (int c = 0; c < numChannels; ++c)
        {
            const float a = std::fabs (channels[c][sampleIndex]);
            if (a > mx) mx = a;
        }
        return mx;
    }

    float s = 0.0f;
    for (int c = 0; c < numChannels; ++c) { const float v = channels[c][sampleIndex]; s += v * v; }
    return std::sqrt (s / (float) numChannels);
}

} // namespace felitronics::dynamics

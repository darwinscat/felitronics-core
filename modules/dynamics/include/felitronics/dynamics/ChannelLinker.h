// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <algorithm>
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

// The DETECTOR GATE. A key sample is a MEASUREMENT, and one bad measurement must not be able to reach
// recursive state — `state += c*(x - state)` never decays a NaN out, so a single poisoned sample is
// permanent. Two things get in:
//   * non-finite -> 0. An Inf reaches the follower as `in`, and even an INSTANT coefficient poisons it
//     (`0.0f * (env - Inf)` is `0*Inf` = NaN, not 0);
//   * |x| > 1e6 (+120 dBFS) -> clamped, because `MeanPower` and the RMS detector both SQUARE the sample:
//     a finite 1e20 squares to +Inf in float and lands in exactly the case above. With the clamp each
//     lane contributes at most 1e12 to the MeanPower sum, so the count would have to pass 3e26 lanes
//     before that sum could overflow — the key channel count is deliberately unbounded, so the bound
//     that matters is per-lane, not per-frame.
// Bit-transparent by construction: any finite sample within +-1e6 — which is every real audio sample by
// 120 dB — takes neither branch and comes out unchanged. Same shape and same bound as the gate in
// limiter::TruePeakLimiter and saturation::Saturator, deliberately.
inline float detectorGate (float x) noexcept
{
    return std::clamp (std::isfinite (x) ? x : 0.0f, -1.0e6f, 1.0e6f);
}

namespace detail
{
    // The link arithmetic, written ONCE. `get(c)` yields channel c's sample for this frame — raw for
    // linkAmplitude, gated for linkAmplitudeGated — so the two spellings cannot drift apart into two
    // different detectors, which is the whole point of P2 (the offline analysis pass must be able to
    // reproduce the running compressor's envelope EXACTLY, not approximately).
    template <class Fetch>
    inline float linkAmplitudeImpl (LinkMode mode, int numChannels, Fetch&& get) noexcept
    {
        // A count of ZERO must not read plane 0. `<= 1` used to fold 0 and negatives into "mono", so
        // `(key, 0)` — the natural spelling of "no key" for a caller holding the detector directly —
        // dereferenced a pointer it was told nothing was behind. Unreachable through Compressor, which
        // resolves "no key" into a source before it gets here; reachable by anything else, and the
        // offline analysis pass that comes next drives this object on its own.
        if (numChannels <= 0) return 0.0f;
        if (numChannels == 1) return std::fabs (get (0));

        if (mode == LinkMode::Max)
        {
            float mx = 0.0f;
            for (int c = 0; c < numChannels; ++c)
            {
                const float a = std::fabs (get (c));
                if (a > mx) mx = a;
            }
            return mx;
        }

        float s = 0.0f;
        for (int c = 0; c < numChannels; ++c) { const float v = get (c); s += v * v; }
        return std::sqrt (s / (float) numChannels);
    }
}

// One linked detector level (linear amplitude) from the channel frame at `sampleIndex`. Stateless.
inline float linkAmplitude (LinkMode mode, const float* const* channels, int numChannels, int sampleIndex) noexcept
{
    return detail::linkAmplitudeImpl (mode, numChannels, [&] (int c) noexcept { return channels[c][sampleIndex]; });
}

// Same, with every sample passed through detectorGate() FIRST. Gating after the link would be too late:
// `MeanPower` squares inside the link, so a large-but-finite frame overflows to +Inf there and a gate on
// the RESULT would then read that as invalid and substitute silence — a loud frame mistaken for nothing.
// This is the entry point a detector uses; `linkAmplitude` stays the raw primitive.
inline float linkAmplitudeGated (LinkMode mode, const float* const* channels, int numChannels, int sampleIndex) noexcept
{
    return detail::linkAmplitudeImpl (mode, numChannels, [&] (int c) noexcept { return detectorGate (channels[c][sampleIndex]); });
}

} // namespace felitronics::dynamics

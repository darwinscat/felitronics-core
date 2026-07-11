// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/core/Fft.h>                  // DefaultRealFft — the analysis FFT (scalar reference)
#include <felitronics/measurement/LevelProbe.h>    // THE contract: kShapeHz — the reference's power spectrum

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

//==============================================================================
// felitronics::measurement — reference-unity IR normalization: the ONE gain that makes an IR pass
// the family reference stimulus (LevelProbe.h) at unity RMS. The other half of the family "unity"
// contract — the probe MEASURES a stage against the reference, this NORMALIZES an IR to it, and
// both take their numbers from LevelProbe.h, so "unity" means the same thing family-wide.
//
// WHY (from OrbitCab's cab loader): vendor cab IRs are peak-normalized bandpasses, so their
// passband sits well above unity (+13…+18 dB measured across a factory set) — convolving them RAW
// makes the wet path that much hotter than the dry. Normalized, a cab contributes TONE, not gain:
// browsing compares timbre at matched loudness and an output trim stays inside its range. Apply
// the returned gain to the FINAL taps (post-trim, post-resample) — ONE common gain across channels
// (stereo imaging untouched). NOT for reverb/decay IRs: RMS-normalizing a mostly-decayed tail
// blows up the wet gain — leave those at their authored level.
//
// HOW: 1 / (reference RMS gain) of the IR: G² = Σ w(f)·P(f) / Σ w(f) over the positive-frequency
// bins (DC excluded), where P(f) is the channel-mean power response |H(f)|² and
// w(f) = 1 / (1 + (f/levelprobe::kShapeHz)²) is the one-pole-shaped reference's power spectrum.
// Frequency domain (one real FFT per channel) — equals convolving the reference noise through the
// IR and reading the RMS ratio, without needing a signal.
//
// OFFLINE, MESSAGE-THREAD-ONLY (allocates; one FFT per channel per call — cheap enough that a
// trim-handle drag re-runs it per mouse move). FLOAT taps/gain ON PURPOSE (documented exception to
// the module's double convention): this is the exact arithmetic products shipped — the returned
// gain defines what users HEAR at unity, so its outputs are golden and must not drift.
//==============================================================================
namespace felitronics::measurement
{

// Clamp: a pathological IR can't blast or vanish. ±30 dB clears every REAL library measured
// (186 commercial IRs: reference gains +12…+24.1 dB — 96 kHz packs run the hottest) with
// headroom, so the guard only ever bites genuine garbage.
constexpr float  kIrNormMinDb  = -30.0f;
constexpr float  kIrNormMaxDb  = 30.0f;
constexpr double kIrRefFloorDb = -60.0;   // near-silent IR ⇒ keep g = 1 (don't amplify garbage)
// Analysis window: the first second. An IR's tail past that carries so little band energy that it
// moves the reference gain by < 0.1 dB (measured on a factory set down to 10% trims) — the cap
// keeps interactive re-normalization (trim drags) light.
constexpr double kIrAnalysisSeconds = 1.0;

// The reference-unity gain for `numChannels` planar channels of `numSamples` taps at `sampleRate`
// (the rate the IR will be CONVOLVED at — measure the final, resampled taps). Returns 1 on a
// near-silent IR (below kIrRefFloorDb) or degenerate input — non-positive sizes/rate, a null
// channel array, or ANY null channel pointer (a null channel is a caller bug, not a silent
// channel: don't guess a gain from what can't be measured — leave the IR alone, the same stance
// as the near-silent floor). Otherwise the clamped 1/G above.
inline float referenceUnityGain (const float* const* channels, int numChannels, int numSamples,
                                 double sampleRate)
{
    if (channels == nullptr || numChannels <= 0 || numSamples <= 0 || sampleRate <= 0.0)
        return 1.0f;
    for (int c = 0; c < numChannels; ++c)                        // a null channel poisons the whole call
        if (channels[c] == nullptr) return 1.0f;

    const int cap = std::min (numSamples, (int) std::lround (kIrAnalysisSeconds * sampleRate));
    int N = 256; while (N < cap * 2 && N < (1 << 21)) N <<= 1;   // dense DTFT sampling of the IR
    core::fft::DefaultRealFft fft;
    if (! fft.prepare (N)) return 1.0f;

    std::vector<float> padded ((std::size_t) N, 0.0f);
    std::vector<float> spec ((std::size_t) core::fft::DefaultRealFft::spectrumFloats (N));
    std::vector<double> power ((std::size_t) (N / 2 + 1), 0.0);
    for (int c = 0; c < numChannels; ++c)
    {
        std::fill (padded.begin(), padded.end(), 0.0f);
        std::copy (channels[c], channels[c] + cap, padded.begin());
        fft.forward (padded.data(), spec.data());
        power[0]                     += (double) spec[0] * spec[0];                 // DC (unused below)
        power[(std::size_t) (N / 2)] += (double) spec[1] * spec[1];                 // Nyquist
        for (int k = 1; k < N / 2; ++k)
            power[(std::size_t) k] += (double) spec[(std::size_t) (2 * k)] * spec[(std::size_t) (2 * k)]
                                    + (double) spec[(std::size_t) (2 * k + 1)] * spec[(std::size_t) (2 * k + 1)];
    }

    double num = 0.0, den = 0.0;
    const double chInv = 1.0 / std::max (1, numChannels);
    for (int k = 1; k <= N / 2; ++k)                             // DC excluded: not audio
    {
        const double f = (double) k * sampleRate / N;
        const double w = 1.0 / (1.0 + (f / levelprobe::kShapeHz) * (f / levelprobe::kShapeHz));
        num += w * power[(std::size_t) k] * chInv;
        den += w;
    }
    const double gSq = den > 0.0 ? num / den : 0.0;
    if (gSq < std::pow (10.0, kIrRefFloorDb / 10.0))             // near-silent IR: leave it alone
        return 1.0f;
    const float gDb = (float) (-10.0 * std::log10 (gSq));
    return std::pow (10.0f, std::clamp (gDb, kIrNormMinDb, kIrNormMaxDb) / 20.0f);
}

// Mono convenience — same contract, one channel.
inline float referenceUnityGain (const float* taps, int numSamples, double sampleRate)
{
    const float* channels[1] { taps };
    return referenceUnityGain (channels, 1, numSamples, sampleRate);
}

} // namespace felitronics::measurement

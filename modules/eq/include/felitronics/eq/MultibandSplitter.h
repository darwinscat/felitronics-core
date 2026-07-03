// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/eq/Crossover2.h>
#include <felitronics/eq/Svf.h>

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace felitronics::eq
{

//==============================================================================
// felitronics::eq::MultibandSplitter<MaxBands> — split a signal into N frequency bands with FLAT
// reconstruction (Σ bands = the input through an allpass — |·|=1, not bit-transparent; phase is rotated).
// The foundation for multiband compressor / saturator / width.
//
// Naive cascading of 2-way LR4 crossovers is NOT flat: the low band of crossover-1 skips crossover-2's
// phase, so the sum combs at f1. Fix: pass each already-split lower band
// through the LR4 ALLPASS of every LATER crossover, using the identity LP4+HP4 = allpass. Then the bands
// telescope: Σ = A1·A2·…·A(N-1)·x, |·|=1.
//
//   B0 = L1·A2·A3·…·A(N-1)          (lowest band, compensated by all later crossovers)
//   Bi = H1·…·Hi·L(i+1)·A(i+2)·…    (0 < i < N-1)
//   B(N-1) = H1·H2·…·H(N-1)         (highest band, no compensation)
//
// Optimisation we verified: that LR4 allpass is only 2ND-ORDER —
// LP4+HP4 = (s⁴+1)/(s²+√2s+1)² = (s²−√2s+1)/(s²+√2s+1) — so each compensation is ONE `Svf` in AllPass mode
// (NOT a 4-SVF Crossover, NOT a wrong 2nd-order cascade). Zero latency, per-channel state, RT-safe.
template <int MaxBands = 4>
class MultibandSplitter
{
    static_assert (MaxBands >= 2 && MaxBands <= 8, "MultibandSplitter supports 2..8 bands");

public:
    static constexpr int kMaxCrossovers = MaxBands - 1;

    void prepare (double sampleRate, int numChannels) noexcept
    {
        fs_ = sampleRate;
        for (int k = 0; k < kMaxCrossovers; ++k)
        {
            xover_[k].prepare (sampleRate, numChannels);
            for (int j = 0; j < kMaxCrossovers; ++j) comp_[k][j].prepare (sampleRate, numChannels);
        }
        defaultFreqs();
        applyFreqs();
        reset();
    }

    void reset() noexcept
    {
        for (int k = 0; k < kMaxCrossovers; ++k)
        {
            xover_[k].reset();
            for (int j = 0; j < kMaxCrossovers; ++j) comp_[k][j].reset();
        }
    }

    bool setNumBands (int bands) noexcept
    {
        if (bands < 2 || bands > MaxBands) return false;
        numBands_ = bands; numXovers_ = bands - 1;
        defaultFreqs(); applyFreqs(); reset();   // topology change → clear stale filter state
        return true;
    }

    int numBands() const noexcept { return numBands_; }

    // `hz` ascending, length = numBands-1. Sanitised (finite), clamped, kept strictly increasing — and the
    // lower bound is itself capped at the ceiling so std::clamp can never see lo > hi (UB) near Nyquist.
    void setCrossovers (const float* hz, int count) noexcept
    {
        const int n = std::min (count, numXovers_);
        const float ceil = (float) (0.45 * fs_);
        float prev = 10.0f;
        for (int k = 0; k < n; ++k)
        {
            float f  = std::isfinite (hz[k]) ? hz[k] : prev * 2.0f;        // NaN/Inf guard
            const float lo = std::min (prev * 1.0001f, ceil);             // never lo > hi
            f = std::clamp (f, lo, ceil);
            freqs_[k] = f; prev = f;
        }
        for (int k = n; k < numXovers_; ++k)                              // repair any untouched suffix
        {
            const float lo = std::min (prev * 1.0001f, ceil);
            freqs_[k] = std::clamp (freqs_[k], lo, ceil); prev = freqs_[k];
        }
        applyFreqs();
    }

    // One sample, one channel → numBands outputs (already compensated, so sum is flat).
    void splitSample (int ch, float x, float* bandOut) noexcept
    {
        float residual = x;
        for (int k = 0; k < numXovers_; ++k)
        {
            float low, high;
            xover_[k].processSample (ch, residual, low, high);
            for (int j = k + 1; j < numXovers_; ++j) low = comp_[k][j].processSample (ch, low);  // allpass of later crossovers
            bandOut[k] = low;
            residual = high;
        }
        bandOut[numXovers_] = residual;
    }

    // Reconstruct — just addition (the bands are pre-compensated).
    float sumSample (const float* bandIn) const noexcept
    {
        float s = 0.0f;
        for (int b = 0; b < numBands_; ++b) s += bandIn[b];
        return s;
    }

    void flushDenormals() noexcept
    {
        for (int k = 0; k < kMaxCrossovers; ++k)
        {
            xover_[k].flushDenormals();
            for (int j = 0; j < kMaxCrossovers; ++j) comp_[k][j].flushDenormals();
        }
    }

    static constexpr int latencySamples() noexcept { return 0; }

private:
    void defaultFreqs() noexcept   // geometric spread ~60 Hz..12 kHz so any N is usable out of the box
    {
        for (int k = 0; k < numXovers_; ++k)
            freqs_[k] = (float) (60.0 * std::pow (200.0, (double) (k + 1) / (double) (numXovers_ + 1)));
    }

    void applyFreqs() noexcept
    {
        for (int k = 0; k < numXovers_; ++k)
        {
            xover_[k].setFrequency (freqs_[k]);
            for (int i = 0; i < k; ++i) comp_[i][k].setParams (FilterType::AllPass, freqs_[k], kQ, 0.0);  // band i < crossover k
        }
    }

    static constexpr double kQ = 0.7071067811865476;
    double fs_ = 48000.0;
    int numBands_ = MaxBands, numXovers_ = MaxBands - 1;
    // size_t casts on the template-dependent bounds: GCC's -Wsign-conversion flags a dependent int bound.
    float freqs_[(std::size_t) kMaxCrossovers] = {};
    Crossover2 xover_[(std::size_t) kMaxCrossovers];
    Svf comp_[(std::size_t) kMaxCrossovers][(std::size_t) kMaxCrossovers];   // comp_[band][crossover]: allpass for band < crossover
};

} // namespace felitronics::eq

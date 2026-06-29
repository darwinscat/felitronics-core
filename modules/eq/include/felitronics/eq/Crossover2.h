// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/eq/Svf.h>

namespace felitronics::eq
{

//==============================================================================
// felitronics::eq::Crossover2 — a 4th-order Linkwitz-Riley 2-way crossover (low | high) built on the
// Cytomic SVF: low = LP4(x), high = HP4(x), each a cascade of two Butterworth (Q=1/√2) sections. The LR4
// identity low + high = allpass (flat magnitude, in-phase) is what lets a MultibandSplitter reconstruct
// flat. Minimum-phase, ZERO latency, per-channel state. This is the reusable splitter primitive extracted
// from stereo::MonoBass.
//
// Note (used by MultibandSplitter): that allpass is only 2nd-order — LP4+HP4 = (s⁴+1)/(s²+√2s+1)² =
// (s²−√2s+1)/(s²+√2s+1) — so band-compensation needs a single `Svf` in AllPass mode at the same fc/Q, not
// a second Crossover2.
class Crossover2
{
public:
    void prepare (double sampleRate, int numChannels) noexcept
    {
        lp1_.prepare (sampleRate, numChannels); lp2_.prepare (sampleRate, numChannels);
        hp1_.prepare (sampleRate, numChannels); hp2_.prepare (sampleRate, numChannels);
        setFrequency (freq_);
    }

    void reset() noexcept { lp1_.reset(); lp2_.reset(); hp1_.reset(); hp2_.reset(); }

    void setFrequency (float hz) noexcept
    {
        freq_ = hz;
        lp1_.setParams (FilterType::LowPass,  freq_, kQ, 0.0); lp2_.setParams (FilterType::LowPass,  freq_, kQ, 0.0);
        hp1_.setParams (FilterType::HighPass, freq_, kQ, 0.0); hp2_.setParams (FilterType::HighPass, freq_, kQ, 0.0);
    }

    float frequency() const noexcept { return freq_; }
    static constexpr int latencySamples() noexcept { return 0; }

    void processSample (int ch, float x, float& low, float& high) noexcept
    {
        low  = lp2_.processSample (ch, lp1_.processSample (ch, x));
        high = hp2_.processSample (ch, hp1_.processSample (ch, x));
    }

    void flushDenormals() noexcept { lp1_.flushDenormals(); lp2_.flushDenormals(); hp1_.flushDenormals(); hp2_.flushDenormals(); }

private:
    static constexpr double kQ = 0.7071067811865476;   // 1/√2 Butterworth → cascade = 4th-order Linkwitz-Riley
    float freq_ = 1000.0f;
    Svf lp1_, lp2_, hp1_, hp2_;
};

} // namespace felitronics::eq

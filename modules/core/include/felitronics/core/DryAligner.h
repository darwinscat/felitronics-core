// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <algorithm>
#include <cstddef>
#include <vector>

namespace felitronics::core
{

//==============================================================================
// felitronics::core::DryAligner — emits the input delayed by a per-block-VARIABLE number of samples,
// through a fixed-capacity per-channel circular ring fed EVERY block (so it stays warm and the delay
// tap can change block-to-block — a stage toggling between modes, or a model load moving a rate-match
// latency — with no reallocation and no glitch). Its job: hold a bypassed / powered-off stage at the
// ACTIVE stage's PDC, so toggling that stage's power never changes the plugin's reported latency (no
// host re-sync GAP) and an off↔active crossfade blends time-aligned signals (no comb / level JUMP).
//
// Same integer-ring algorithm as core::DelayLine, deliberately a SEPARATE type — the shape differs:
// DelayLine is one channel, one sample per call, with a delay you set and hold; DryAligner is
// N channels, one block per call, tap passed per block, and it STAGES the delayed copy into internal
// scratch so the caller can snapshot the raw input BEFORE an in-place stage overwrites the buffer and
// still read the aligned dry afterwards. Folding either into the other would cost the per-sample API
// or the stage-then-read semantics.
//
// Bit-exact: the delayed sample is a pure copy of the input from exactly D samples earlier — no
// arithmetic, no interpolation. D == 0 → an exact identity copy.
//
// Contract: advance() with numChannels ≤ the prepared channel count and numSamples ≤ the prepared
// maxBlock; call it BEFORE any in-place stage overwrites `io`; delayed(ch) stays valid until the
// next advance()/prepare().
//
// 🔴 RT: advance() never allocates, locks or throws — both buffers are sized once in prepare().
//==============================================================================
class DryAligner
{
public:
    // `capacity` must exceed any latency the tap will ever request (the delay is clamped to
    // [0, capacity-1]). Size it to the largest stage latency plus margin.
    void prepare (int numChannels, int maxBlock, int capacity)
    {
        channels_ = std::max (1, numChannels);
        cap_      = std::max (2, capacity);
        maxBlock_ = std::max (1, maxBlock);
        ring_.assign    ((std::size_t) channels_ * (std::size_t) cap_,      0.0f);
        scratch_.assign ((std::size_t) channels_ * (std::size_t) maxBlock_, 0.0f);
        pos_ = 0;
    }

    void reset() noexcept
    {
        std::fill (ring_.begin(),    ring_.end(),    0.0f);
        std::fill (scratch_.begin(), scratch_.end(), 0.0f);
        pos_ = 0;
    }

    int capacity()    const noexcept { return cap_; }
    int numChannels() const noexcept { return channels_; }

    // Stage the internal scratch = `io` delayed by `delaySamples` (0 → an exact copy) and advance the
    // ring by numSamples. Read the result via delayed(ch). `io` may be the caller's live buffer — call
    // this BEFORE any in-place stage overwrites it.
    void advance (const float* const* io, int numChannels, int numSamples, int delaySamples) noexcept
    {
        const int C = cap_;
        const int D = std::clamp (delaySamples, 0, C - 1);
        for (int ch = 0; ch < numChannels; ++ch)
        {
            float* r = ring_.data()    + (std::size_t) ch * (std::size_t) C;
            float* s = scratch_.data() + (std::size_t) ch * (std::size_t) maxBlock_;
            const float* x = io[ch];
            int wp = pos_;
            for (int i = 0; i < numSamples; ++i)
            {
                r[wp] = x[i];                       // write freshest input
                int rp = wp - D; if (rp < 0) rp += C;
                s[i] = r[rp];                       // read the sample D samples ago (D == 0 → == input)
                if (++wp >= C) wp = 0;
            }
        }
        pos_ = (pos_ + numSamples) % C;             // every channel advanced identically → one commit
    }

    const float* delayed (int channel) const noexcept
    {
        return scratch_.data() + (std::size_t) channel * (std::size_t) maxBlock_;
    }

private:
    std::vector<float> ring_    { 0.0f, 0.0f };   // circular delay history (capacity ≥ any stage latency);
    std::vector<float> scratch_ { 0.0f };         //   default-VALID mono delay-0 state until prepare()
    int channels_ = 1, cap_ = 2, maxBlock_ = 1;
    int pos_ = 0;                                 // shared write cursor into ring_
};

} // namespace felitronics::core

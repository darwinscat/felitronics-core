// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/multiband/MultibandProcessor.h>
#include <felitronics/stereo/StereoWidth.h>

namespace felitronics::multiband
{

//==============================================================================
// felitronics::multiband::MultibandWidth — per-band stereo width: a stereo::StereoWidth on every LR4 band of
// the splitter (low-narrow / high-wide is the mastering norm). PURE COMPOSITION — zero new DSP (DSP council).
//
// Because StereoWidth touches ONLY the Side (the Mid passes untouched at unity output gain), the GLOBAL
// mono-fold invariant survives ANY per-band widening: ½(L+R) folds to the allpass-reconstructed input Mid,
// Σ band-Mids = allpass(M), no matter the individual side gains — so widening can never collapse the mono
// sum. All widths = 1 → the splitter's allpass reconstruction (magnitude-flat, not raw bit-exact — the same
// contract as Crossover2 / MonoBass). Differing per-band widths reconstruct cleanly: the in-phase LR4 split
// makes a per-band side-gain a smooth frequency-dependent width, never a comb. Mono (< 2 ch) is a TRUE
// passthrough (no allpass phase). Per-band output gain is deliberately NOT exposed — it would break the
// mono-fold guarantee. Zero latency. RT-safe: prepare() allocates, process() is alloc/lock/throw-free.
template <int MaxBands = 4>
class MultibandWidth
{
public:
    bool prepare (double sampleRate, int maxBlock, int maxChannels)
    {
        return mb_.prepare (sampleRate, maxBlock, maxChannels, 0,            // StereoWidth is zero-latency → no align delay
                            [&] (stereo::StereoWidth& w) { w.prepare (sampleRate, maxBlock, maxChannels); });
    }

    void reset() noexcept { mb_.reset(); }

    bool setNumBands (int n) noexcept { return mb_.setNumBands (n); }
    void setCrossovers (const float* hz, int count) noexcept { mb_.setCrossovers (hz, count); }
    int  numBands() const noexcept { return mb_.numBands(); }

    void  setBandWidth (int b, float width) noexcept { if ((unsigned) b < (unsigned) MaxBands) mb_.band (b).setWidth (width); }   // 0=mono, 1=neutral, 2=wide; OOB ignored
    float bandWidth (int b) const noexcept { return ((unsigned) b < (unsigned) MaxBands) ? mb_.band (b).width() : 1.0f; }
    void  setBandBypass (int b, bool x) noexcept { mb_.setBandBypass (b, x); }
    void  setBandSolo   (int b, bool x) noexcept { mb_.setBandSolo (b, x); }
    void  setMix        (float m) noexcept       { mb_.setMix (m); }

    int latencySamples() const noexcept { return mb_.latencySamples(); }     // 0 (LR4 split + StereoWidth)

    void process (float* const* io, int numChannels, int n) noexcept
    {
        if (numChannels < 2) return;                                         // mono → true passthrough (no allpass phase)
        mb_.process (io, numChannels, n);
    }

private:
    MultibandProcessor<stereo::StereoWidth, MaxBands> mb_;
};

} // namespace felitronics::multiband

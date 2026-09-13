// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/core/Config.h>
#include <felitronics/core/Math.h>
#include <felitronics/oversampling/PolyphaseOversampler.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <vector>

namespace felitronics::analysis
{

//==============================================================================
// felitronics::analysis::ReferenceTruePeakMeter — the core's REFERENCE true peak: the instrument a delivered
// file is CERTIFIED with, and therefore the one anything that promises a delivered ceiling has to aim with.
//
// The filter: `oversampling::PolyphaseOversampler` at 4x, 32 taps per phase (a 128-tap Kaiser-sinc prototype,
// cutoff 0.90 x base Nyquist, beta 9), at EVERY sample rate; the running maximum of |x| over the 4x stream,
// floored at the sample peak. This is the arithmetic `fcore::Probe` (tools/fcore_probe.h) has always run —
// `fcore_measure` and the browser's `fc_probe` report it, it is what the site's baseline harness and the EBU
// Tech 3341 runs were judged by — moved here so that a second caller cannot be a second copy. The probe now
// measures THROUGH this class, so "the solver's reading equals `fcore_measure` on the delivered file" is a
// statement about one class and one set of doubles, not about two implementations that happen to agree.
//
// WHY THE CORE HOLDS TWO TRUE-PEAK METERS, AND WHICH ONE TO REACH FOR. `analysis::TruePeakMeter` is the spec's
// short filter (12 taps per phase, cutoff at the base Nyquist, factor 4/2/1 by rate): cheap, built for a live
// chain. It is a different filter and reads a different number, and on high delivery rates it stops
// interpolating altogether (2x at 88.2-176.4 kHz, a plain sample peak at >= 176.4). How far the two readings
// sit apart — on which material, at which rate — is not stated here: it is measured and pinned by
// `felitronics_truepeak_instrument_gap_tests`, the only owner of that number. The rule is the design, not the
// number: a PROMISE about a delivered file is aimed and certified with this class; a display is free to use
// the cheap one.
//
// NOR IS THIS THE TRUE PEAK OF THE SIGNAL — it is the reference's reading of it. Its 0.90 x Nyquist cutoff
// with 32 taps per phase droops near the top of the band, so on content the filter attenuates (a full-band
// click, a tone above ~0.4 fs) a band-limited oracle reads higher than this does. That is a property of the
// certifying instrument, which is chosen by contract, and it is pinned in the same gap suite (its groups "the worst
// cell against the truth" and "the reference reads under the truth on a full-band click") rather than argued away here.
//
// THE DRAIN IS THE CALLER'S TO ASK FOR, and asking is not optional for a whole-file reading. The FIR is causal
// with a group delay of (128-1)/2 oversampled samples, so the reconstruction of the last ~16 input samples has
// not left the filter when the input ends. The sample-peak floor does not stand in for it: a buffer ending on ten
// samples at 0.95 reads 0.95 undrained and 1.062468 drained — an over the floor cannot see (the 4x maximum alone
// is 0.015 undrained). `drain()` pushes `kTapsPerPhase` samples of digital silence through every prepared channel
// — exactly what `fcore::Probe::finish()` always did — and after it the history IS that silence, so a stream that
// continues is simply a stream with 32 zeros in it. One more drain adds nothing: the ring is already zero.
//
// A CHANNEL THAT STOPS IS DRAINED, NOT DROPPED. Law 11 (DSP-ARCHITECTURE.md §2) makes a narrower call legal and says
// a stopped channel must not replay its history when it returns (11a) — for a delay line, drop it. For a running
// MAXIMUM, dropping is wrong: the samples still inside the filter were submitted, and their reconstruction is part
// of the reading. Dropped, a stereo stream whose right channel ends on a peak and is then followed by one mono
// block lost an over (+0.53 dBTP read as -0.45). So the falling edge drains the stopped channel — the same 32 zeros
// `drain()` feeds — which is what 11c asks of a pause (the channel heard silence) and leaves its history silent for
// 11a. READ-ONLY: `io` is sampled, never written. RT-safe: prepare() allocates, process()/drain() do not.
class ReferenceTruePeakMeter
{
public:
    static constexpr int kFactor       = 4;     // \ THE reference. Changing either moves every certificate the
    static constexpr int kTapsPerPhase = 32;    // / core and the site have issued — see the header above
    // The internal walk: process() takes any `n` and upsamples it in pieces of this many frames into one shared
    // scratch, so the scratch is a constant 16 KB rather than a function of the caller's block. Chunking cannot
    // change a bit: the oversampler's ring makes each output a pure function of the samples seen so far.
    static constexpr int kChunk = 1024;

    // WHAT prepare() ASKS THE HEAP FOR (law 11d): one `PolyphaseOversampler` per channel at the reference
    // topology, and the scratch. Asked of a FRESH meter; a prepared one keeps storage that still fits. All zeros
    // (`ok == false`) exactly where prepare() refuses the same arguments — which are its arguments, maxBlock included.
    struct Storage
    {
        bool          ok = false;
        std::uint64_t oversamplerBytes = 0;     // all channels together
        std::size_t   scratchFloats    = 0;
        std::uint64_t bytes() const noexcept
        {
            return oversamplerBytes + (std::uint64_t) sizeof (float) * (std::uint64_t) scratchFloats;
        }
    };
    static Storage storageFor (double sampleRate, int maxBlock, int maxChannels) noexcept
    {
        Storage st;
        if (! validRate (sampleRate) || maxBlock < 1 || maxChannels < 1 || maxChannels > core::kMaxChannels) return st;
        oversampling::PolyphaseOversampler::Storage one;
        if (! oversampling::PolyphaseOversampler::storageFor (kFactor, 1, kTapsPerPhase, one)) return st;
        st.oversamplerBytes = one.bytes() * (std::uint64_t) maxChannels;
        st.scratchFloats    = (std::size_t) kChunk * (std::size_t) kFactor;
        st.ok = true;
        return st;
    }

    // Neither the rate nor the block shapes the filter or the storage — the reference is 4x at every rate and walks
    // any `n` in kChunk pieces — but both are arguments, and law 11b makes every argument binding: a rate no audio
    // stream can have, or a block of no samples, is refused rather than ignored.
    [[nodiscard]] bool prepare (double sampleRate, int maxBlock, int maxChannels)
    {
        prepared_ = false;
        const Storage st = storageFor (sampleRate, maxBlock, maxChannels);
        if (! st.ok) return false;
        for (int c = 0; c < maxChannels; ++c)
            if (! os_[(std::size_t) c].prepare (kFactor, 1, kTapsPerPhase)) return false;
        scratch_.assign (st.scratchFloats, 0.0f);
        channels_ = maxChannels;
        reset();
        prepared_ = true;
        return true;
    }

    bool isPrepared() const noexcept { return prepared_; }

    void reset() noexcept
    {
        for (int c = 0; c < channels_; ++c) os_[(std::size_t) c].reset();
        maxOs_ = samplePeak_ = blockMax_ = 0.0;
        ranNc_ = 0;
    }

    static constexpr int latencySamples() noexcept { return 0; }

    // The reading: the running maximum since reset(), never below the sample peak — the reconstruction passes
    // through the samples by construction, and the floor makes a filter-side mistake (a mis-sized drain, a
    // wrong prototype) impossible to hide under it. LINEAR is the bit-exactness surface; dB is derived.
    double truePeakLinear()   const noexcept { return std::max (maxOs_, samplePeak_); }
    double samplePeakLinear() const noexcept { return samplePeak_; }
    // The maximum of the last accepted process() call alone, oversampled and grid, for a caller that asks what
    // one block contained. It includes whatever of the previous block the FIR was still carrying.
    double truePeakLinearBlock() const noexcept { return blockMax_; }
    double truePeakDb()      const noexcept { return core::gainToDb (truePeakLinear()); }
    double truePeakDbBlock() const noexcept { return core::gainToDb (blockMax_); }

    [[nodiscard]] bool process (const float* const* io, int numChannels, int n) noexcept
    {
        if (numChannels < 0 || n < 0) return false;
        if (! prepared_) return false;
        if (numChannels > channels_) return false;                   // law 11(b): the width is a LIMIT
        if (n == 0) return true;                                     // no samples: no time, no edge
        const int nc = numChannels;
        if (nc > 0)                                                  // a plane that is not there is refused
        {                                                            // BEFORE the edge moves anything
            if (io == nullptr) return false;
            for (int c = 0; c < nc; ++c) if (io[c] == nullptr) return false;
        }
        double blockMax = 0.0;
        for (int c = nc; c < ranNc_; ++c) blockMax = std::max (blockMax, drainChannel (c));   // these channels stopped
        ranNc_ = nc;
        if (nc == 0) { maxOs_ = std::max (maxOs_, blockMax); blockMax_ = blockMax; return true; }

        for (int c = 0; c < nc; ++c)
        {
            const float* x = io[c];
            float grid = 0.0f;
            for (int i = 0; i < n; ++i) grid = std::max (grid, std::fabs (x[i]));
            samplePeak_ = std::max (samplePeak_, (double) grid);
            blockMax    = std::max (blockMax, (double) grid);

            for (int off = 0; off < n; )                                   // advanced by what was taken: no
            {                                                               // `off + kChunk` past INT_MAX
                const int m = std::min (kChunk, n - off);
                blockMax = std::max (blockMax, upsampleMax (c, x + off, m));
                off += m;
            }
        }
        maxOs_    = std::max (maxOs_, blockMax);
        blockMax_ = blockMax;
        return true;
    }

    // See the header: `kTapsPerPhase` zeros through every PREPARED channel, folded into the running maximum.
    // A channel that never ran, or stopped and was drained then, has a zero ring and contributes exactly zero.
    void drain() noexcept
    {
        if (! prepared_) return;
        for (int c = 0; c < channels_; ++c)
            maxOs_ = std::max (maxOs_, drainChannel (c));
    }

private:
    static bool validRate (double fs) noexcept { return fs > 0.0 && std::isfinite (fs); }

    double drainChannel (int c) noexcept
    {
        const float zeros[kTapsPerPhase] {};
        return upsampleMax (c, zeros, kTapsPerPhase);
    }

    // One piece of one channel through its oversampler; the maximum |x| of the 4x output. Taken as a float
    // maximum and widened once, which is the same number as widening every sample first (a maximum of floats is
    // one of them) and a loop a compiler can vectorise. NaN cannot win a std::max whose left side is not NaN.
    double upsampleMax (int c, const float* x, int m) noexcept
    {
        const float* in[1] { x };
        float*       out[1] { scratch_.data() };
        os_[(std::size_t) c].upsample (in, 1, m, out);
        float mx = 0.0f;
        const int upN = m * kFactor;
        for (int k = 0; k < upN; ++k) mx = std::max (mx, std::fabs (scratch_[(std::size_t) k]));
        return (double) mx;
    }

    std::array<oversampling::PolyphaseOversampler, core::kMaxChannels> os_ {};
    std::vector<float> scratch_;
    double maxOs_ = 0.0, samplePeak_ = 0.0, blockMax_ = 0.0;
    int    channels_ = 0, ranNc_ = 0;
    bool   prepared_ = false;
};

} // namespace felitronics::analysis

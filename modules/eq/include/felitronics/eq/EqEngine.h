// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/eq/EqBand.h>
#include <felitronics/eq/EqTypes.h>
#include <felitronics/analysis/SpectrumTap.h>

#include <algorithm>
#include <array>
#include <vector>
#include <atomic>
#include <cmath>
#include <cmath>
#include <complex>

namespace felitronics::eq
{

//==============================================================================
// teq::EqEngine — a bank of N EQ bands in series, plus pre/post spectrum taps and a
// magnitude-response readout for the GUI curve. Framework-agnostic: it processes raw float
// buffers and owns no parameter system. A host adapter (e.g. a JUCE AudioProcessor) maps its
// parameters into BandParams via setBand().
//
// Threading: prepare()/setBand()/process() run on the SAME thread (typically the audio thread —
// the host adapter reads its atomic parameters there and feeds setBand() at the top of process()),
// or are synchronised externally; the engine takes no internal lock. process() does no
// alloc/lock/IO. setSpectrumActive() is safe from any thread (atomic). The host should enable
// FTZ/DAZ on the audio thread (e.g. juce::ScopedNoDenormals); the engine also flushes per block.
class EqEngine
{
public:
    static constexpr int kMaxBands    = 24;
    static constexpr int kMaxChannels = EqBand::kMaxChannels;

    // REFUSES a configuration it cannot honour. `sampleRate` reaches every coefficient in every band, so
    // a zero or non-finite one does not degrade the EQ, it poisons it — measured, prepare(NaN, 64, 2)
    // then a 0.25 input gives 63 non-finite samples out of 64. This was the second half of P6 F12:
    // Saturator's identical hole was closed then and this one was not, so only the mastering chain, which
    // validates at its own входе, was protected — a direct consumer of `eq` got NaN in silence.
    // [[nodiscard]] and a refusal, matching TruePeakLimiter (P12) and Compressor (P2), rather than the
    // silent substitution some void prepare()s in this repo do: a caller that believes it prepared and
    // processes at the wrong rate is the failure those two exist to prevent. Spelled positively so NaN fails.
    [[nodiscard]] bool prepare (double sampleRate, int maxBlock, int numChannels) noexcept
    {
        prepared_ = false;                              // any early return below leaves the engine unprepared
        // Same domain test as the bands below, and for the same reason: a band clamps its frequency to
        // [10 Hz, 0.49*fs], so a rate under 20.41 Hz hands `std::clamp` a lo above its hi. The engine
        // repeats it rather than relying on the bands because it is the module's front door and refuses
        // before allocating the scratch buffer. Spelled positively so NaN fails.
        if (! (std::isfinite (sampleRate) && 0.49 * sampleRate >= 10.0 && sampleRate <= 3.0e6)) return false;
        fs = sampleRate;
        if (numChannels < 1 || numChannels > kMaxChannels) return false;   // law 11(b): BINDING — a width
        ch = numChannels;                                                  // it clamped was a width it lied about
        maxBlock_ = maxBlock > 0 ? maxBlock : 0;
        // Sidechain scratch: the SECTION INPUT, preserved before any band touches the signal. A
        // dynamics layer must detect on this and not on a band's own input — in a series chain that
        // input is the previous bands' OUTPUT, so their moving deltas would modulate later detectors
        // at overlapping frequencies and the chain would pump. Allocated here, never in process().
        // maxBlock 0 (or a consumer that never asks) costs nothing.
        scratch_.assign ((std::size_t) (maxBlock_ * ch), 0.0f);
        for (int c = 0; c < kMaxChannels; ++c)
            scPtr_[c] = (c < ch && maxBlock_ > 0) ? scratch_.data() + (std::size_t) c * (std::size_t) maxBlock_
                                                  : nullptr;
        scValid_ = 0;
        bool ok = true;
        for (auto& b : bands) ok = b.prepare (fs, ch) && ok;   // every band, then the verdict — never short-circuit
        if (! ok) return false;
        reset();
        prepared_ = true;
        return true;
    }

    // Capture this block's SECTION INPUT into the scratch buffer, then hand it back. Call at the top
    // of a block, BEFORE processing any band; the pointers stay valid until the next capture. Returns
    // nullptr if prepare() was given no maxBlock or the block is larger than promised — a caller that
    // gets nullptr must skip its dynamics rather than detect on the wrong signal.
    const float* const* captureSectionInput (const float* const* channels, int numChannels, int numSamples) noexcept
    {
        // LAW 11: this is the second two-phase API the law names, and it owes the same two things as
        // NoiseGate::analyse. (1) The WIDTH is a limit, not a clamp: a capture wider than the prepared
        // engine used to be silently narrowed, so a consumer asking for 3 planes got 2 and no word of
        // it. (2) A REFUSED capture must be inert — clobbering `scValid_` on refusal destroyed a
        // capture the previous call had legitimately produced (measured: a valid length of 32 became 0
        // when a 65-sample capture was refused against a 64-sample capacity).
        if (! prepared_) return nullptr;
        if (numChannels < 0 || numSamples < 0) return nullptr;
        if (numChannels > ch || numSamples > maxBlock_ || scratch_.empty()) return nullptr;
        // `n == 0` is the ONE true no-op — the same clause this branch's sibling (NoiseGate::analyse)
        // was corrected for in this very change. Zeroing the bookkeeping here destroyed a capture the
        // previous call had legitimately produced: a valid 32 became 0.
        if (numSamples == 0) return nullptr;
        const int nc = numChannels;
        if (nc == 0) { scValid_ = scNc_ = 0; return nullptr; }   // a real, EMPTY capture: no columns
        for (int c = 0; c < nc; ++c)
        {
            scPtr_[c] = scratch_.data() + (std::size_t) c * (std::size_t) maxBlock_;   // a capture wider than
            std::copy (channels[c], channels[c] + numSamples, scPtr_[c]);              // the last one restores
        }
        // The capture has a WIDTH as well as a length, and only the length used to be recorded. A narrow
        // capture after a wide one left the columns past nc pointing at the PREVIOUS block's audio, so a
        // consumer reading one column too far detected on a signal that is not this block's — silently, and
        // on data that looks perfectly plausible. Columns outside the capture are handed back as nullptr
        // instead: the same refusal this method already makes for a block that is too long.
        for (int c = nc; c < kMaxChannels; ++c) scPtr_[c] = nullptr;
        scValid_ = numSamples;
        scNc_    = nc;
        return scPtr_;
    }

    int  sectionInputSamples()  const noexcept { return scValid_; }
    int  sectionInputChannels() const noexcept { return scNc_; }   // width of the last capture, 0 if none

    // Direct band access, so a composition layer can interleave "compute this band's deltas" with
    // "run this band" — the order dynamics requires, and one this engine deliberately does not
    // hard-code, since it knows nothing about what drives it.
    EqBand&       bandAt (int i)       noexcept { return bands[(size_t) std::clamp (i, 0, kMaxBands - 1)]; }
    const EqBand& bandAt (int i) const noexcept { return bands[(size_t) std::clamp (i, 0, kMaxBands - 1)]; }
    static constexpr int bandCount() noexcept { return kMaxBands; }   // not "numBands": several statics here already take that as a PARAMETER

    void reset() noexcept
    {
        scValid_ = scNc_ = 0;   // a stream restart invalidates any captured section input
        for (auto& b : bands) b.reset();
        inTap.reset();
        outTap.reset();
    }

    // A STOP, not a restart: clear what the previous audio left behind and leave the parameter epoch
    // alone — see `EqBand::clearAudioState()`. This is what a consumer that skips the engine for a while
    // wants at the bypass edge, and it is what `reset()` used to do before reset() became a real stream
    // restart. Reaching for reset() there would now also snap every ramp in flight and make the next
    // parameter write a hard step.
    void clearAudioState() noexcept
    {
        scValid_ = scNc_ = 0;   // the captured section input belongs to audio that will not be continued
        for (auto& b : bands) b.clearAudioState();
        inTap.reset();
        outTap.reset();
    }

    void setBand (int i, const BandParams& p) noexcept
    {
        if (i >= 0 && i < kMaxBands) bands[(size_t) i].setParams (p);
    }

    BandParams band (int i) const noexcept
    {
        return (i >= 0 && i < kMaxBands) ? bands[(size_t) i].params() : BandParams {};
    }

    void setSpectrumActive (bool active) noexcept { spectrumOn.store (active, std::memory_order_relaxed); }

    // Audio thread. In-place over `numChannels` planar buffers of `numSamples`.
    // Law 11 (DSP-ARCHITECTURE.md §2). The engine used to drop a zero-width call outright while the very
    // same band, driven directly, advanced its audio-time grid — so one EqBand behaved two ways depending
    // on which side of this method it was called from. Measured on one 500 -> 4000 Hz glide and 10240
    // samples of (nch = 0) calls: the two had diverged by 11.08 dB. The bands are now always driven,
    // because a call carrying samples spent that much audio time whoever was listening.
    [[nodiscard]] bool process (float* const* channels, int numChannels, int numSamples) noexcept
    {
        if (numChannels < 0 || numSamples < 0) return false;   // malformed
        if (! prepared_) return false;                         // a refused configuration has nothing to run
        if (numChannels > ch) return false;                    // width is a LIMIT — law 11(b)
        if (numSamples == 0) return true;                      // no samples: no time, no edge, nothing
        const int nc = numChannels;

        // The taps read channel 0, so they run only when there IS one; the bands run either way.
        const bool spec = spectrumOn.load (std::memory_order_relaxed) && nc > 0;
        if (spec) for (int n = 0; n < numSamples; ++n) inTap.push (channels[0][n]);
        bool ok = true;
        for (auto& b : bands) ok = b.processBlock (channels, nc, numSamples) && ok;
        if (spec) for (int n = 0; n < numSamples; ++n) outTap.push (channels[0][n]);
        return ok;
    }

    // Total magnitude (dB) of all active bands on one stereo axis at a real frequency — best-effort
    // LIVE readout (reads the bands' current coefficients; may briefly race a concurrent param change
    // on the audio thread). Defaults to the Mid axis (== the L=R response when no band is split). For a
    // guaranteed race-free GUI curve use the static magnitudeDbFor() below.
    double magnitudeDb (double freqHz, Axis a = Axis::Mid) const noexcept
    {
        const double w = 2.0 * kPi * freqHz / fs;
        std::complex<double> h { 1.0, 0.0 };
        for (auto& b : bands) h *= b.response (w, a);
        return 20.0 * std::log10 (std::max (1e-9, std::abs (h)));
    }

    void magnitudeResponse (const double* freqs, double* outDb, int n) const noexcept
    {
        for (int i = 0; i < n; ++i) outDb[i] = magnitudeDb (freqs[i]);
    }

    // Race-free GUI curve: total magnitude (dB) from a caller-owned BandParams array (touches no
    // engine state). Pass the same params you fed setBand(), plus sampleRate(). The plain overload is
    // the Mid axis (== the L=R response when no band is split); pass an Axis to pick another.
    static double magnitudeDbFor (const BandParams* bandsIn, int numBands, double freqHz, double fs) noexcept
    {
        return magnitudeDbFor (bandsIn, numBands, freqHz, fs, Axis::Mid);
    }

    static double magnitudeDbFor (const BandParams* bandsIn, int numBands, double freqHz, double fs, Axis a) noexcept
    {
        const double w = 2.0 * kPi * freqHz / fs;
        return 20.0 * std::log10 (std::max (1e-9, std::abs (compositeResponse (bandsIn, numBands, fs, w, a))));
    }

    // Fill out[0..n-1] with the LINEAR-magnitude composite on a uniform grid f[k] = k·fs/(2·(n-1)) —
    // i.e. k=0 at DC … k=n-1 at Nyquist (so for an FFT of size N pass n = N/2 + 1). Chosen stereo axis,
    // race-free. A host builds a linear-phase FIR from this (zero-phase magnitude → IFFT → window).
    static void magnitudeGridFor (const BandParams* bandsIn, int numBands, double fs,
                                  float* out, int n, Axis a) noexcept
    {
        const double nyq = 0.5 * fs;
        for (int k = 0; k < n; ++k)
        {
            const double f = (n > 1) ? nyq * (double) k / (double) (n - 1) : 0.0;
            const double w = 2.0 * kPi * f / fs;
            out[k] = (float) std::abs (compositeResponse (bandsIn, numBands, fs, w, a));
        }
    }

    // Fill the FOUR entries of the ZERO-PHASE 2×2 matrix (matrixResponseZeroPhase) on the same uniform
    // grid as magnitudeGridFor. The entries are REAL but SIGNED (the off-diagonals can go negative) — a
    // linear-phase FIR host inverse-transforms each to a symmetric IR (a signed real spectrum stays
    // symmetric, keeping the sign). Only the Full topology needs all four; the diagonal topologies use
    // magnitudeGridFor on the axis composites. Race-free (caller-owned snapshot).
    static void matrixGridZeroPhase (const BandParams* bandsIn, int numBands, double fs,
                                     float* outLL, float* outLR, float* outRL, float* outRR, int n) noexcept
    {
        const double nyq = 0.5 * fs;
        for (int k = 0; k < n; ++k)
        {
            const double f = (n > 1) ? nyq * (double) k / (double) (n - 1) : 0.0;
            const double w = 2.0 * kPi * f / fs;
            const ResponseMatrix H = matrixResponseZeroPhase (bandsIn, numBands, fs, w);
            outLL[k] = (float) H.hLL.real();  outLR[k] = (float) H.hLR.real();
            outRL[k] = (float) H.hRL.real();  outRR[k] = (float) H.hRR.real();
        }
    }

    double sampleRate() const noexcept { return fs; }

    SpectrumTap& inputTap()  noexcept { return inTap; }
    SpectrumTap& outputTap() noexcept { return outTap; }

private:
    double fs = 44100.0;
    int    ch = 2;
    std::atomic<bool> spectrumOn { false };

    std::array<EqBand, kMaxBands> bands;
    std::vector<float> scratch_;                 // section-input capture (see captureSectionInput)
    float*             scPtr_[kMaxChannels] {};
    bool               prepared_ = false;          // true only after a fully-successful prepare()
    int                maxBlock_ = 0, scValid_ = 0, scNc_ = 0;
    SpectrumTap inTap, outTap;
};

} // namespace felitronics::eq

// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/eq/EqBand.h>
#include <felitronics/eq/EqTypes.h>
#include <felitronics/analysis/SpectrumTap.h>

#include <algorithm>
#include <array>
#include <atomic>
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

    void prepare (double sampleRate, int /*maxBlock*/, int numChannels) noexcept
    {
        fs = sampleRate;
        ch = numChannels < 1 ? 1 : (numChannels > kMaxChannels ? kMaxChannels : numChannels);
        for (auto& b : bands) b.prepare (fs, ch);
        reset();
    }

    void reset() noexcept
    {
        for (auto& b : bands) b.reset();
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
    void process (float* const* channels, int numChannels, int numSamples) noexcept
    {
        const int nc = numChannels < ch ? numChannels : ch;
        if (nc <= 0 || numSamples <= 0) return;

        const bool spec = spectrumOn.load (std::memory_order_relaxed);
        if (spec) for (int n = 0; n < numSamples; ++n) inTap.push (channels[0][n]);
        for (auto& b : bands) b.processBlock (channels, nc, numSamples);
        if (spec) for (int n = 0; n < numSamples; ++n) outTap.push (channels[0][n]);
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
    SpectrumTap inTap, outTap;
};

} // namespace felitronics::eq

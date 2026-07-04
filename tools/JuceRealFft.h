// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// BENCHMARK-ONLY JUCE/vDSP reference FFT backend. Deliberately NOT part of the core build — the fftbench tool
// includes it only when FELITRONICS_BENCH_JUCE=1 (which FetchContent's JUCE). Lets the bench run
// MatrixConvolver<JuceRealFft> — our identical partitioned engine with JUCE's FFT (vDSP on Apple) — beside the
// pffft backend. NOTE: the spectral MAC here is scalar, so this is an FFT-backend proxy, not a stand-in for
// juce::dsp::Convolution's own optimized convolver.

#pragma once

#include <felitronics/core/Fft.h>

#include <juce_dsp/juce_dsp.h>

#include <algorithm>
#include <cstring>
#include <memory>
#include <vector>

namespace felitronics::bench
{

class JuceRealFft
{
public:
    static constexpr int spectrumFloats (int n) noexcept { return 2 * n; }   // juce real-only: interleaved complex, 2N floats

    bool prepare (int n) noexcept
    {
        fft_.reset();
        scratch_.clear();
        n_ = 0;
        if (! core::fft::isPow2 (n) || n < 4) return false;
        int order = 0;
        while ((1 << order) < n) ++order;
        try
        {
            fft_ = std::make_unique<juce::dsp::FFT> (order);
            scratch_.assign (static_cast<std::size_t> (2 * n), 0.0f);
        }
        catch (...) { fft_.reset(); scratch_.clear(); return false; }
        n_ = n;
        return true;
    }

    int size() const noexcept { return n_; }

    void forward (const float* real, float* spec) noexcept
    {
        if (fft_ == nullptr || n_ <= 0) return;
        std::fill (scratch_.begin(), scratch_.end(), 0.0f);
        std::memcpy (scratch_.data(), real, static_cast<std::size_t> (n_) * sizeof (float));
        fft_->performRealOnlyForwardTransform (scratch_.data(), false);
        std::memcpy (spec, scratch_.data(), static_cast<std::size_t> (2 * n_) * sizeof (float));
    }

    void inverse (const float* spec, float* real) noexcept
    {
        if (fft_ == nullptr || n_ <= 0) return;
        std::memcpy (scratch_.data(), spec, static_cast<std::size_t> (2 * n_) * sizeof (float));
        fft_->performRealOnlyInverseTransform (scratch_.data());
        std::memcpy (real, scratch_.data(), static_cast<std::size_t> (n_) * sizeof (float));
    }

    void spectralMultiplyAdd (const float* a, const float* b, float* acc) const noexcept
    {
        if (n_ <= 0) return;
        for (int i = 0; i < n_; ++i)
        {
            const float ar = a[2 * i], ai = a[2 * i + 1], br = b[2 * i], bi = b[2 * i + 1];
            acc[2 * i]     += ar * br - ai * bi;
            acc[2 * i + 1] += ar * bi + ai * br;
        }
    }

private:
    std::unique_ptr<juce::dsp::FFT> fft_;
    std::vector<float> scratch_;
    int n_ = 0;
};

static_assert (core::fft::RealFftBackend<JuceRealFft>, "JuceRealFft must satisfy RealFftBackend");

} // namespace felitronics::bench

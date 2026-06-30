// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/core/Fft.h>

#include <algorithm>
#include <cmath>
#include <vector>

namespace felitronics::lineareq
{

//==============================================================================
// felitronics::lineareq::MixedPhaseFir — render a magnitude response |H(ω)| into an FIR whose phase is a
// CONTINUOUS blend between zero (linear) and the minimum phase: φ(ω) = k · φ_min(ω), k ∈ [0,1].
//
//   k = 0  → zero-phase (symmetric) FIR — full linear phase: flat phase, max pre-ringing, N/2 group delay.
//   k = 1  → minimum-phase FIR — causal, NO pre-ring, ~0 bulk delay, but the full (analog-like) phase shift.
//   0<k<1  → "Natural phase": most of linear's flat phase, but the impulse is shifted forward so pre-ringing
//            and bulk delay drop sharply. The mastering middle ground.
//
// METHOD — cepstral, exact (DSP council: codex + deepseek, both verified). The cepstrum is the log-spectrum
// domain, so a LINEAR blend of cepstra is a linear blend of (log|H| + jφ): the magnitude term log|H| is the
// same in both endpoints, so it is preserved EXACTLY, while the phase term blends 0 → φ_min linearly. Steps:
//   1. logM = ln(max(|H|, floor))               (zero-phase spectrum)         → inverse FFT → real cepstrum c
//   2. blend-fold c: causal side ×(1+k), anti-causal side ×(1−k)              (= (1−k)·c_even + k·c_minphase)
//   3. forward FFT → imaginary part is φ = k·φ_min ; rebuild H = |H|·e^{jφ}   → inverse FFT → impulse h[0..D)
// The design size D must be ≫ the kept IR length (cepstral time-aliasing on steep filters/notches — use ≥ 8×).
// No centred window (it would gut the front-loaded impulse); the consumer truncates with a TAIL taper + a
// causal bulk shift, and reports that shift as latency.
//
// RT-UNSAFE (message thread): prepare()/build() allocate / run FFTs. Pairs with ConvolutionEngine for the swap.
template <class Fft = core::fft::DefaultRealFft>
class MixedPhaseFir
{
public:
    // designSize D: a power of two, ≥ ~8× the FIR length you intend to keep. Message thread.
    bool prepare (int designSize)
    {
        if (! core::fft::isPow2 (designSize)) return false;
        D_ = designSize;
        if (! fft_.prepare (D_)) return false;
        specF_ = Fft::spectrumFloats (D_);
        ceps_.assign ((std::size_t) D_, 0.0f);
        spec_.assign ((std::size_t) specF_, 0.0f);
        h_.assign    ((std::size_t) D_, 0.0f);
        return true;
    }

    int size() const noexcept { return D_; }

    // mag: |H| at D/2+1 points (DC … Nyquist), linear scale (≥ 0). k ∈ [0,1] phase blend. Returns the
    // internal D-sample impulse whose DFT magnitude == mag and phase == k·φ_min (circular; the consumer
    // shifts it causal). RT-UNSAFE.
    const float* build (const float* mag, float k) noexcept
    {
        const int H = D_ / 2;
        k = std::clamp (k, 0.0f, 1.0f);
        constexpr float floorMag = 1.0e-5f;

        // 1. log-magnitude as a zero-phase spectrum → real cepstrum (inverse is 1/D-normalised).
        spec_[0] = std::log (std::max (mag[0],             floorMag));   // DC   (real)
        spec_[1] = std::log (std::max (mag[(std::size_t) H], floorMag)); // Nyq  (real)
        for (int b = 1; b < H; ++b)
        {
            spec_[(std::size_t) (2 * b)]     = std::log (std::max (mag[(std::size_t) b], floorMag));
            spec_[(std::size_t) (2 * b + 1)] = 0.0f;
        }
        fft_.inverse (spec_.data(), ceps_.data());           // c[0..D), real + even

        // 2. blend-fold: c_blend = (1−k)·c_even + k·c_minphase. c is even (c[D−n]==c[n]); the minimum-phase
        //    fold doubles the causal side and zeroes the anti-causal, so the blend is:
        //      n∈[1,D/2)  → (1+k)·c[n]      (causal)
        //      n∈(D/2,D)  → (1−k)·c[n]      (anti-causal)   ; c[0], c[D/2] unchanged.
        for (int n = 1; n < H; ++n)
        {
            const float c = ceps_[(std::size_t) n];          // == ceps_[D-n] (even)
            ceps_[(std::size_t) n]        = c * (1.0f + k);
            ceps_[(std::size_t) (D_ - n)] = c * (1.0f - k);
        }

        // 3. forward FFT of the blended cepstrum → log H = log|H| + j·(k·φ_min); take the phase, rebuild
        //    H = |H|·e^{jφ}, inverse → the mixed-phase impulse. (DC/Nyquist are forced real.)
        fft_.forward (ceps_.data(), spec_.data());           // spec_ = [DC, Nyq, re,im, …]
        spec_[0] = mag[0];                                   // DC   (phase 0)
        spec_[1] = mag[(std::size_t) H];                     // Nyq  (phase 0)
        for (int b = 1; b < H; ++b)
        {
            const float phi = spec_[(std::size_t) (2 * b + 1)];   // imag of FFT(c_blend) = k·φ_min
            const float m   = mag[(std::size_t) b];
            spec_[(std::size_t) (2 * b)]     = m * std::cos (phi);
            spec_[(std::size_t) (2 * b + 1)] = m * std::sin (phi);
        }
        fft_.inverse (spec_.data(), h_.data());              // h[0..D): mixed-phase impulse (circular)
        return h_.data();
    }

private:
    Fft fft_;
    int D_ = 0, specF_ = 0;
    std::vector<float> ceps_, spec_, h_;
};

} // namespace felitronics::lineareq

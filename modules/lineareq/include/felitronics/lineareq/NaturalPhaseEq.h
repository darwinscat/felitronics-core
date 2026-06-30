// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/eq/EqEngine.h>
#include <felitronics/convolution/ConvolutionEngine.h>
#include <felitronics/lineareq/MixedPhaseFir.h>
#include <felitronics/core/Fft.h>
#include <felitronics/core/Config.h>

#include <algorithm>
#include <cmath>
#include <vector>

namespace felitronics::lineareq
{

//==============================================================================
// felitronics::lineareq::NaturalPhaseEq — the "Natural phase" rendering of an eq::EqEngine bank: the
// mastering middle ground between Zero-Latency (minimum-phase IIR) and Linear (zero-phase FIR).
//
// It builds a MIXED-PHASE FIR (MixedPhaseFir, cepstral blend φ = k·φ_min) of the bank's composite
// magnitude, per Mid/Side axis, and convolves with the same click-free ConvolutionEngine as LinearPhaseEq.
// At k≈0.5 the phase stays close to flat (most of linear's benefit) while pre-ringing and bulk delay drop
// sharply — so the reported PDC latency is only (1−k)·L/2, FAR below linear's L/2, and transients are not
// smeared by symmetric pre-ring.
//
//  • Design FFT D = 8·L masks cepstral time-aliasing on steep filters (council: codex + deepseek).
//  • The mixed-phase impulse is shifted causal by `bulkDelay = (1−k)·L/2` (that shift IS the latency) and
//    truncated to L taps with a TAIL taper (NO centred window — that would gut the front-loaded impulse).
//  • A flat EQ renders a unit impulse at `bulkDelay` ⇒ exact unity-gain pass-through at the reported latency.
//
// Same threading contract as LinearPhaseEq: setBands()/buildFir()/prepare() are RT-UNSAFE (message thread,
// host-serialised); process() is RT-safe. To change quality or the phase blend k, re-prepare().
class NaturalPhaseEq
{
public:
    static constexpr int kNumQuality = 4;                       // kept-IR length (L) per quality — lighter than Linear
    static int firSizeForQuality (int q) noexcept
    {
        static constexpr int sizes[kNumQuality] = { 2048, 4096, 8192, 16384 };
        return sizes[std::clamp (q, 0, kNumQuality - 1)];
    }

    // quality picks L (the convolution IR length). k ∈ [0,1] blends phase (0 linear … 1 minimum). RT-UNSAFE.
    bool prepare (double sampleRate, int maxBlock, int numChannels, int quality, float k = 0.5f) noexcept
    {
        fs_       = sampleRate > 0.0 ? sampleRate : 48000.0;
        maxBlock_ = std::max (1, maxBlock);
        channels_ = std::clamp (numChannels, 1, core::kMaxChannels);
        k_        = std::clamp (k, 0.0f, 1.0f);
        L_        = firSizeForQuality (quality);
        D_        = 8 * L_;                                      // cepstral design size (≥ 8× the kept IR)
        bulkDelay_ = (int) std::lround ((1.0f - k_) * (float) L_ * 0.5f);   // causal shift = reported latency

        if (! mp_.prepare (D_)) return false;
        magBuf_.assign ((std::size_t) (D_ / 2 + 1), 0.0f);
        firMid_.assign ((std::size_t) L_, 0.0f);
        firSide_.assign((std::size_t) L_, 0.0f);
        computeTaper();

        int part = 64; while (part < maxBlock_) part <<= 1;     // partition ≥ maxBlock, pow2
        const int warmXfade = std::max (2 * part, (int) std::lround (0.02 * fs_));   // short anti-click fade (design B)
        if (! conv_.prepare (part, L_, warmXfade, 2)) return false;   // Mid IR (ch0) + Side IR (ch1)
        return true;
    }

    void reset() noexcept { conv_.reset(); }

    int  firSize() const noexcept { return L_; }
    int  latencySamples() const noexcept { return bulkDelay_; }   // the causal bulk shift (convolver adds 0)
    bool isBusy() const noexcept { return conv_.isBusy(); }
    float blend() const noexcept { return k_; }

    // RT-UNSAFE (message thread): rebuild the Mid+Side mixed-phase FIRs from a snapshot, hand to the
    // convolver for a click-free swap. False if a swap is still crossfading (host coalesces with latest).
    bool setBands (const eq::BandParams* bands, int numBands) noexcept
    {
        buildFir (bands, numBands, false, firMid_.data());     // Mid axis
        buildFir (bands, numBands, true,  firSide_.data());    // Side axis
        const float* irs[2] { firMid_.data(), firSide_.data() };
        return conv_.setIr (irs, 2, L_);
    }

    // RT-UNSAFE: build ONE axis's mixed-phase FIR into out[0..L). Exposed for tests / a host that drives
    // the convolver itself. magnitude → MixedPhaseFir (phase k·φ_min) → causal shift by bulkDelay + taper.
    void buildFir (const eq::BandParams* bands, int numBands, bool side, float* out) noexcept
    {
        eq::EqEngine::magnitudeGridFor (bands, numBands, fs_, magBuf_.data(), D_ / 2 + 1, side);
        const float* h = mp_.build (magBuf_.data(), k_);       // D-point mixed-phase impulse (circular)

        // Extract L causal taps: the peak sits at h[0], pre-ring wraps to h[D-1..]; shift right by bulkDelay
        // so out[bulkDelay] = h[0]. Tail (+ light head) taper suppresses truncation ripple.
        for (int i = 0; i < L_; ++i)
        {
            int idx = i - bulkDelay_;
            if (idx < 0) idx += D_;
            out[i] = h[(std::size_t) idx] * taper_[(std::size_t) i];
        }
    }

    // RT-safe (audio thread), in place. Stereo (≥2 ch) → M/S; mono → the Mid IR. `n` ≤ maxBlock.
    void process (float* const* io, int numChannels, int n) noexcept
    {
        if (n <= 0) return;
        if (numChannels >= 2)
        {
            float* L = io[0];
            float* R = io[1];
            for (int i = 0; i < n; ++i) { const float m = 0.5f * (L[i] + R[i]), s = 0.5f * (L[i] - R[i]); L[i] = m; R[i] = s; }
            conv_.process (io, io, 2, n);
            for (int i = 0; i < n; ++i) { const float l = L[i] + R[i], r = L[i] - R[i]; L[i] = l; R[i] = r; }
        }
        else if (numChannels == 1)
        {
            conv_.process (io, io, 1, n);
        }
    }

private:
    // Tukey-ish taper: a short rising head (suppress any clipped pre-ring) + a longer falling tail
    // (suppress post-ring truncation). NOT a centred window — the bulk of the impulse stays untouched.
    void computeTaper() noexcept
    {
        taper_.assign ((std::size_t) L_, 1.0f);
        const int head = std::min (bulkDelay_ / 2, L_ / 16);                 // ≤ half the pre-ring region
        const int tail = std::max (1, L_ / 8);                              // last ~12.5%
        for (int i = 0; i < head; ++i)
            taper_[(std::size_t) i] = 0.5f * (1.0f - std::cos (core::kPi * (float) i / (float) head));
        for (int i = 0; i < tail; ++i)
        {
            const int idx = L_ - 1 - i;
            taper_[(std::size_t) idx] = 0.5f * (1.0f - std::cos (core::kPi * (float) i / (float) tail));
        }
    }

    double fs_ = 48000.0;
    int maxBlock_ = 512, channels_ = 2;
    int L_ = 4096, D_ = 32768, bulkDelay_ = 1024;
    float k_ = 0.5f;

    MixedPhaseFir<core::fft::DefaultRealFft> mp_;                            // cepstral mixed-phase FIR designer
    convolution::ConvolutionEngine<core::fft::DefaultRealFft, 2> conv_;      // Mid (ch0) + Side (ch1), click-free swap

    std::vector<float> magBuf_, firMid_, firSide_, taper_;
};

} // namespace felitronics::lineareq

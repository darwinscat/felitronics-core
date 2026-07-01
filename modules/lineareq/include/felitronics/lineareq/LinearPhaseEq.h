// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/eq/EqEngine.h>
#include <felitronics/convolution/ConvolutionEngine.h>
#include <felitronics/core/Fft.h>
#include <felitronics/core/Config.h>

#include <algorithm>
#include <cmath>
#include <vector>

namespace felitronics::lineareq
{

//==============================================================================
// felitronics::lineareq::LinearPhaseEq — a linear-phase rendering of an eq::EqEngine bank. From the bank's
// composite MAGNITUDE response (eq::EqEngine::magnitudeGridFor, per Mid/Side axis) it builds a symmetric
// zero-phase FIR and convolves: the EQ's amplitude shape, ZERO phase shift, constant group delay = N/2.
//
// Built to FIX the plugin prototype it's extracted from (which stuttered
// rebuilding a huge FIR on every node drag):
//   • CORE owns NO thread. setBands()/buildFir() are RT-UNSAFE (FFT/IFFT/window) — the HOST calls them off
//     the audio thread and DEBOUNCES (rebuild on drag-settle, not per move). process() is RT-safe. setBands()
//     returns false while a swap is mid-crossfade (isBusy()) → the host coalesces with the latest snapshot.
//   • The FIR is N+1 symmetric taps (type-I linear phase): zero-phase magnitude grid → core::fft inverse
//     (1/N-normalised) → fftshift → 4-term Blackman-Harris (peaks EXACTLY 1.0 at the centre tap N/2). So a
//     flat EQ renders an exact unit impulse at N/2 — no run-time gain hack needed.
//   • The convolution is convolution::ConvolutionEngine (zero-latency, partitioned), 2 banks: the Mid IR on
//     channel 0, the Side IR on channel 1, with ONE lockstep click-free crossfade so an IR swap never moves
//     the stereo image. Stereo in → M/S encode (½ convention) → convolve → decode; mono → the Mid IR.
//
// Latency = N/2 (the FIR group delay; the convolver adds 0) — 2048…65536 samples, so this is an OFFLINE /
// MASTERING tool, not for live monitoring. CONTRACTS: setBands()/buildFir()/prepare() are message-thread
// (the host serialises them — a single producer; concurrent setBands() is NOT safe); reset() is audio-thread
// (or externally synced) and cancels any pending swap. To change quality (N) re-prepare(). RT-safe: process()
// never allocates/locks/throws.
class LinearPhaseEq
{
public:
    static constexpr int kNumQuality = 5;                       // Low / Medium / High / Very High / Maximum
    static int firSizeForQuality (int q) noexcept
    {
        // Power-of-two FIR lengths (our FFT is radix-2). Latency = N/2 ≈ FabFilter's linear-phase ladder
        // — 2048 / 4096 / 8192 / 16384 / 65536 vs FabFilter's 3072 / 5120 / 9216 / 17408 / 66560 @ 44.1k:
        // comparable, converging to equal at Maximum.
        static constexpr int sizes[kNumQuality] = { 4096, 8192, 16384, 32768, 131072 };
        return sizes[std::clamp (q, 0, kNumQuality - 1)];
    }

    // RT-UNSAFE (message thread): allocates the FIR/FFT/convolver. quality ∈ [0,3] picks the FIR length.
    bool prepare (double sampleRate, int maxBlock, int numChannels, int quality) noexcept
    {
        prepared_ = false;                                     // any early return below leaves it unprepared
        fs_ = sampleRate > 0.0 ? sampleRate : 48000.0;
        maxBlock_ = std::max (1, maxBlock);
        channels_ = std::clamp (numChannels, 1, core::kMaxChannels);
        N_ = firSizeForQuality (quality);

        if (! buildFft_.prepare (N_)) return false;
        magBuf_.assign ((std::size_t) (N_ / 2 + 1), 0.0f);
        spec_.assign   ((std::size_t) N_, 0.0f);
        time_.assign   ((std::size_t) N_, 0.0f);
        firMid_.assign ((std::size_t) (N_ + 1), 0.0f);
        firSide_.assign((std::size_t) (N_ + 1), 0.0f);
        computeWindow();

        int part = 64; while (part < maxBlock_) part <<= 1;     // partition ≥ maxBlock, pow2
        // SHORT crossfade for interactive swaps. ConvolutionEngine (design B) keeps ONE warm input history
        // shared by two IR slots, so a swapped IR is response-correct immediately — only a short anti-click
        // fade is needed (≥ the partition P, to cover the one-chunk pendingTail recompute). The engine
        // derives its OWN long fade for the COLD first activation, so an EQ band drag lands in ~tens of ms,
        // not a full FIR (the old N-long fade is what made interactive dragging lag ~seconds). isBusy() is
        // then true only for that short window → the host's coalescing rebuild stays responsive.
        const int warmXfade = std::max (2 * part, (int) std::lround (0.02 * fs_));   // ≥ 2P and ≥ ~20 ms
        if (! conv_.prepare (part, N_ + 1, warmXfade, 2)) return false;   // Mid IR (ch0) + Side IR (ch1)
        prepared_ = true;                                      // fully built — setBands()/buildFir() may now run
        return true;
    }

    void reset() noexcept { conv_.reset(); }

    int  firSize() const noexcept { return N_; }
    int  latencySamples() const noexcept { return N_ / 2; }     // symmetric FIR group delay (convolver adds 0)
    bool isBusy() const noexcept { return conv_.isBusy(); }     // a swap is mid-crossfade → coalesce setBands()

    // RT-UNSAFE (message thread, off the audio thread): rebuild the Mid+Side FIRs from a caller-owned
    // BandParams snapshot and hand them to the convolver for a click-free swap. Returns false if a swap is
    // still crossfading (isBusy()) — the host should retry with the LATEST snapshot once it's free.
    bool setBands (const eq::BandParams* bands, int numBands) noexcept
    {
        if (! prepared_) return false;                         // unprepared — the FIR/scratch buffers are empty
        buildFir (bands, numBands, false, firMid_.data());     // Mid axis
        buildFir (bands, numBands, true,  firSide_.data());    // Side axis
        const float* irs[2] { firMid_.data(), firSide_.data() };
        return conv_.setIr (irs, 2, N_ + 1);
    }

    // RT-UNSAFE: build ONE axis's symmetric zero-phase FIR into out[0..N] (N+1 taps). Exposed for tests /
    // a host that wants to drive the convolver itself.
    void buildFir (const eq::BandParams* bands, int numBands, bool side, float* out) noexcept
    {
        if (! prepared_) return;                               // unprepared — magBuf_/spec_/time_/window_ are empty
        const int N = N_;
        eq::EqEngine::magnitudeGridFor (bands, numBands, fs_, magBuf_.data(), N / 2 + 1, side);

        spec_[0] = magBuf_[0];                                  // DC (real)
        spec_[1] = magBuf_[(std::size_t) (N / 2)];              // Nyquist (real)
        for (int k = 1; k < N / 2; ++k)                         // zero phase: imag = 0
        {
            spec_[(std::size_t) (2 * k)]     = magBuf_[(std::size_t) k];
            spec_[(std::size_t) (2 * k + 1)] = 0.0f;
        }
        buildFft_.inverse (spec_.data(), time_.data());        // time_[0..N-1], 1/N-normalised, even-symmetric about 0

        for (int i = 0; i <= N; ++i)                            // fftshift (peak → N/2) + window → N+1 symmetric taps
            out[i] = time_[(std::size_t) ((i + N / 2) % N)] * window_[(std::size_t) i];
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
            conv_.process (io, io, 2, n);                       // ch0 = Mid·firMid, ch1 = Side·firSide (in place)
            for (int i = 0; i < n; ++i) { const float l = L[i] + R[i], r = L[i] - R[i]; L[i] = l; R[i] = r; }
        }
        else if (numChannels == 1)
        {
            conv_.process (io, io, 1, n);                       // mono → the Mid IR (bank 0)
        }
    }

private:
    void computeWindow() noexcept                              // 4-term Blackman-Harris over N+1 taps (peak 1.0 at N/2)
    {
        window_.assign ((std::size_t) (N_ + 1), 0.0f);
        constexpr double a0 = 0.35875, a1 = 0.48829, a2 = 0.14128, a3 = 0.01168;
        const double M1 = (double) N_;                          // (taps - 1) = N
        for (int i = 0; i <= N_; ++i)
        {
            const double x = (double) i / M1;
            window_[(std::size_t) i] = (float) (a0 - a1 * std::cos (2.0 * core::kPi * x)
                                                   + a2 * std::cos (4.0 * core::kPi * x)
                                                   - a3 * std::cos (6.0 * core::kPi * x));
        }
    }

    double fs_ = 48000.0;
    int maxBlock_ = 512, channels_ = 2, N_ = 16384;
    bool prepared_ = false;                                    // true only after a fully-successful prepare()

    core::fft::DefaultRealFft buildFft_;                        // size-N IFFT for the FIR design
    convolution::ConvolutionEngine<core::fft::DefaultRealFft, 2> conv_;   // Mid (ch0) + Side (ch1), click-free swap

    std::vector<float> magBuf_, spec_, time_, firMid_, firSide_, window_;
};

} // namespace felitronics::lineareq

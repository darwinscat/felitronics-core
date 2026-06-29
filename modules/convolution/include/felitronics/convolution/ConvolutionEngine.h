// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/convolution/PartitionedConvolver.h>

#include <atomic>

namespace felitronics::convolution
{

//==============================================================================
// felitronics::convolution::ConvolutionEngine — the production, swap-safe, CLICK-FREE convolver,
// MULTI-CHANNEL with a single LOCKSTEP crossfade.
//
// Holds 2 banks × N channels of PartitionedConvolver. setIr() (message thread) builds the new IR into the
// INACTIVE bank for every channel; process() (audio thread) CROSSFADES old→new over `crossfadeSamples`,
// feeding BOTH banks the same input so the new bank's frequency-domain delay line primes during the fade.
//
// LOCKSTEP is the point: there is ONE crossfade position shared by all channels, advanced once per sample
// (not once per channel), so a stereo IR swap can never move the image or decorrelate L/R — every channel
// is at the exact same fade coefficient on every sample. (N independent engines would race: setIr on L then
// R could start the two fades on different blocks.) Race-free via a 3-state atomic (Idle / Pending /
// Crossfading): the message thread only ever touches the inactive bank while Idle, and publishes the build
// with a release store the audio thread acquires. Zero latency.
//
// The audio thread NEVER allocates, locks, or blocks. The adapter coalesces setIr() while isBusy(). For a
// fully-primed swap the crossfade should be ≥ the IR length; cab IRs prime within JUCE's 50 ms fade.
// To match juce::dsp::Convolution exactly the adapter pre-applies JUCE's IR gains (Normalise::yes →
// 0.125/sqrt(max-channel ΣIR²); Normalise::no → irSr/hostSr) and resamples to the host rate first.
//
// `MaxChannels` is the compile-time channel bound (2 = stereo, the plugin default); a mono-only consumer
// uses <Fft, 1>. The FFT is the compile-time seam backend (a JUCE adapter plugs juce::dsp::FFT for desktop
// speed; the scalar reference is the default + the test/wasm path).
template <core::fft::RealFftBackend Fft = core::fft::DefaultRealFft, int MaxChannels = 2>
class ConvolutionEngine
{
public:
    // numChannels defaults to 1 so a mono consumer keeps the 3-arg prepare().
    bool prepare (int partitionSize, int maxIrSamples, int crossfadeSamples, int numChannels = 1)
    {
        if (numChannels < 1 || numChannels > MaxChannels) return false;
        channels = numChannels;
        for (int b = 0; b < 2; ++b)
            for (int c = 0; c < channels; ++c)
                if (! eng[b][c].prepare (partitionSize, maxIrSamples)) return false;
        xfadeLen = crossfadeSamples < 1 ? 1 : crossfadeSamples;
        active = 0; xfadePos = 0;
        state.store (0, std::memory_order_relaxed);
        return true;
    }

    void reset() noexcept
    {
        for (int b = 0; b < 2; ++b)
            for (int c = 0; c < channels; ++c)
                eng[b][c].reset();
        active = 0; xfadePos = 0;
        state.store (0, std::memory_order_relaxed);
    }

    static constexpr int latencySamples() noexcept { return 0; }
    bool isBusy() const noexcept { return state.load (std::memory_order_acquire) != 0; }
    int  numChannels() const noexcept { return channels; }

    // Message thread. Mono IR broadcast to ALL channels (matches juce Stereo::yes: a mono IR on L & R).
    // Returns false if a swap is already pending/crossfading (the caller coalesces with the latest IR).
    bool setIr (const float* ir, int len)
    {
        if (state.load (std::memory_order_acquire) != 0) return false;   // busy
        const int inactive = 1 - active;
        for (int c = 0; c < channels; ++c)
            eng[inactive][c].setIr (ir, len);
        state.store (1, std::memory_order_release);                      // → Pending (publishes the build)
        return true;
    }

    // Message thread. Per-channel IR (a true-stereo IR). irPerCh[c] feeds channel c; if fewer channels are
    // supplied than configured, the last one is broadcast to the remainder. False if busy.
    bool setIr (const float* const* irPerCh, int nch, int len)
    {
        if (nch < 1 || state.load (std::memory_order_acquire) != 0) return false;
        const int inactive = 1 - active;
        for (int c = 0; c < channels; ++c)
            eng[inactive][c].setIr (irPerCh[c < nch ? c : nch - 1], len);
        state.store (1, std::memory_order_release);
        return true;
    }

    // Audio thread. Planar; `in`/`out` may alias (in-place). RT-safe, zero latency.
    void process (const float* const* in, float* const* out, int numChannelsToProcess, int n) noexcept
    {
        const int nc = numChannelsToProcess < channels ? numChannelsToProcess : channels;

        int s = state.load (std::memory_order_acquire);
        if (s == 1) { xfadePos = 0; state.store (2, std::memory_order_relaxed); s = 2; }   // begin crossfade

        if (s != 2)
        {
            for (int c = 0; c < nc; ++c)
                eng[active][c].process (in[c], out[c], n);
            return;
        }

        // Crossfading: ONE shared position, advanced once per sample → all channels at the same fade t.
        const int other = 1 - active;
        int pos = xfadePos;
        for (int i = 0; i < n; ++i)
        {
            const float t = (float) pos / (float) xfadeLen;
            for (int c = 0; c < nc; ++c)
            {
                float a, b;
                eng[active][c].process (&in[c][i], &a, 1);
                eng[other][c].process  (&in[c][i], &b, 1);
                out[c][i] = a * (1.0f - t) + b * t;                      // linear blend (smooth, click-free)
            }

            if (++pos >= xfadeLen)
            {
                active = other;
                state.store (0, std::memory_order_release);              // → Idle (new bank live; old now free)
                const int rem = n - i - 1;
                if (rem > 0)
                    for (int c = 0; c < nc; ++c)
                        eng[active][c].process (&in[c][i + 1], &out[c][i + 1], rem);
                xfadePos = 0;
                return;
            }
        }
        xfadePos = pos;
    }

    // Mono convenience (keeps a single-channel consumer on the 3-arg process()).
    void process (const float* in, float* out, int n) noexcept
    {
        const float* ins[1]  { in };
        float*       outs[1] { out };
        process (ins, outs, 1, n);
    }

private:
    PartitionedConvolver<Fft> eng[2][MaxChannels];   // [bank][channel]; bank = active / inactive
    std::atomic<int> state { 0 };                    // 0 Idle · 1 Pending (built, awaiting audio) · 2 Crossfading
    int active = 0, xfadePos = 0, xfadeLen = 1, channels = 1;
};

} // namespace felitronics::convolution

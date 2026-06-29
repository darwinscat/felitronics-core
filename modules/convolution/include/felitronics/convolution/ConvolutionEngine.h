// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/convolution/PartitionedConvolver.h>

#include <atomic>

namespace felitronics::convolution
{

//==============================================================================
// felitronics::convolution::ConvolutionEngine — the production, swap-safe, CLICK-FREE convolver.
//
// Holds TWO PartitionedConvolvers. setIr() (message thread) builds the new IR into the INACTIVE engine;
// process() (audio thread) CROSSFADES old→new over `crossfadeSamples`, feeding BOTH engines the same
// input so the new engine's frequency-domain delay line primes during the fade. Race-free via a 3-state
// atomic (Idle / Pending / Crossfading): the message thread only ever touches the inactive engine while
// Idle, and publishes the build with a release store the audio thread acquires. Zero latency.
//
// The audio thread NEVER allocates, locks, or blocks. The adapter coalesces setIr() while isBusy()
// (mirrors OrbitCab's existing reload-coalescing). For a fully-primed swap the crossfade should be ≥ the
// IR length — short cab IRs prime within a ~30–60 ms fade; a very long IR gets a small, masked transient.
//
// The FFT is the compile-time seam backend (a JUCE adapter plugs juce::dsp::FFT for desktop speed; the
// scalar reference is the default + the test/wasm path).
template <core::fft::RealFftBackend Fft = core::fft::DefaultRealFft>
class ConvolutionEngine
{
public:
    bool prepare (int partitionSize, int maxIrSamples, int crossfadeSamples)
    {
        if (! eng[0].prepare (partitionSize, maxIrSamples)) return false;
        if (! eng[1].prepare (partitionSize, maxIrSamples)) return false;
        xfadeLen = crossfadeSamples < 1 ? 1 : crossfadeSamples;
        active = 0; xfadePos = 0;
        state.store (0, std::memory_order_relaxed);
        return true;
    }

    void reset() noexcept
    {
        eng[0].reset(); eng[1].reset();
        active = 0; xfadePos = 0;
        state.store (0, std::memory_order_relaxed);
    }

    static constexpr int latencySamples() noexcept { return 0; }
    bool isBusy() const noexcept { return state.load (std::memory_order_acquire) != 0; }

    // Message thread. Builds the new IR into the inactive engine + arms a crossfade. Returns false if a
    // swap is already pending/crossfading (the caller coalesces — retries with the latest IR).
    bool setIr (const float* ir, int len)
    {
        if (state.load (std::memory_order_acquire) != 0) return false;   // busy
        eng[1 - active].setIr (ir, len);                                 // build the inactive engine (it is idle)
        state.store (1, std::memory_order_release);                      // → Pending (publishes the build)
        return true;
    }

    // Audio thread. In place not required. RT-safe.
    void process (const float* in, float* out, int n) noexcept
    {
        int s = state.load (std::memory_order_acquire);
        if (s == 1) { xfadePos = 0; state.store (2, std::memory_order_relaxed); s = 2; }   // begin crossfade
        if (s != 2) { eng[active].process (in, out, n); return; }

        const int other = 1 - active;
        for (int i = 0; i < n; ++i)
        {
            float a, b;
            eng[active].process (&in[i], &a, 1);
            eng[other].process  (&in[i], &b, 1);
            const float t = (float) xfadePos / (float) xfadeLen;
            out[i] = a * (1.0f - t) + b * t;                             // equal-? linear blend (smooth, click-free)

            if (++xfadePos >= xfadeLen)
            {
                active = other;
                state.store (0, std::memory_order_release);              // → Idle (new engine live; old now free)
                if (i + 1 < n) eng[active].process (in + i + 1, out + i + 1, n - i - 1);
                return;
            }
        }
    }

private:
    PartitionedConvolver<Fft> eng[2];
    std::atomic<int> state { 0 };       // 0 Idle · 1 Pending (built, awaiting audio) · 2 Crossfading
    int active = 0, xfadePos = 0, xfadeLen = 1;
};

} // namespace felitronics::convolution

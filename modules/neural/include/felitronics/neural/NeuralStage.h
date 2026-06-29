// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/core/Config.h>
#include <felitronics/neural/Inference.h>

#include <array>
#include <atomic>
#include <cstdint>
#include <memory>

namespace felitronics::neural
{

//==============================================================================
// felitronics::neural::NeuralStage<Backend> — the swap-safe holder that generalizes orbitcab's
// cab::AmpStage. The adapter builds + prepares a new model INSTANCE off the audio thread and hands it in
// via swapPrepared(); process() (audio) uses the live instance and NEVER allocates, locks, or DELETES;
// the replaced instance is retired and freed later on the message thread (collectGarbage()), only once
// the audio thread has provably stepped into a LATER block — so no use-after-free and no audio-thread
// delete. No model loaded → in-place passthrough. Vtable-free: `Backend` is the compile-time type, the
// runtime-swapped thing is just an instance of it.
//
// Threading: process()/reset() = audio thread. prepare()/swapPrepared()/clear()/collectGarbage() =
// message (control) thread, serialized with each other.
template <Inference Backend, int MaxRetired = 8>
class NeuralStage
{
public:
    struct Spec { double sampleRate = 48000.0; int maxBlock = 512; int maxChannels = 2; };

    void prepare (Spec s) noexcept
    {
        spec_ = s;
        Backend* p = live_.load (std::memory_order_acquire);
        if (p) p->prepare (s.sampleRate, s.maxBlock, s.maxChannels);
        latency_.store (p ? p->latencySamples() : 0, std::memory_order_release);
    }

    // Message thread. `next` is already built + prepared. Returns false (and keeps the live model) if
    // the retire queue is full — the caller drains it via collectGarbage() and retries.
    bool swapPrepared (std::unique_ptr<Backend> next) noexcept
    {
        collectGarbage();
        if (live_.load (std::memory_order_acquire) != nullptr && retiredCount_ >= MaxRetired)
            return false;
        const int latency = next ? next->latencySamples() : 0;
        Backend* old = live_.exchange (next.release(), std::memory_order_acq_rel);
        latency_.store (latency, std::memory_order_release);
        retire (old);
        return true;
    }

    bool clear() noexcept { return swapPrepared (nullptr); }
    bool hasModel() const noexcept { return live_.load (std::memory_order_acquire) != nullptr; }

    // Audio thread, in place. RT-safe. No model → passthrough (leaves `io` unchanged).
    void process (float* const* io, int numChannels, int numSamples) noexcept
    {
        audioBlock_.fetch_add (1, std::memory_order_acq_rel);
        Backend* p = live_.load (std::memory_order_acquire);
        if (p) p->process (io, numChannels, numSamples);
    }

    void reset() noexcept { if (Backend* p = live_.load (std::memory_order_acquire)) p->reset(); }
    int  latencySamples() const noexcept { return latency_.load (std::memory_order_acquire); }

    // Message thread. Frees retired instances whose retire-block the audio thread has stepped PAST
    // (now > retiredAtBlock — the `>` is load-bearing: `>=` would free a model the audio may still be
    // inside, a use-after-free). Never called on the audio thread.
    void collectGarbage() noexcept
    {
        const std::uint64_t now = audioBlock_.load (std::memory_order_acquire);
        int w = 0;
        for (int r = 0; r < retiredCount_; ++r)
        {
            if (now > retired_[(std::size_t) r].atBlock)
                retired_[(std::size_t) r].object.reset();           // delete — message thread only
            else
            {
                if (w != r) retired_[(std::size_t) w] = std::move (retired_[(std::size_t) r]);
                ++w;
            }
        }
        retiredCount_ = w;
    }

private:
    struct Retired { std::unique_ptr<Backend> object; std::uint64_t atBlock = 0; };

    void retire (Backend* old) noexcept
    {
        if (! old) return;
        retired_[(std::size_t) retiredCount_++] =
            Retired { std::unique_ptr<Backend> (old), audioBlock_.load (std::memory_order_acquire) };
    }

    Spec spec_ {};
    std::atomic<Backend*>      live_ { nullptr };
    std::atomic<std::uint64_t> audioBlock_ { 0 };
    std::atomic<int>           latency_ { 0 };
    std::array<Retired, (std::size_t) MaxRetired> retired_ {};
    int retiredCount_ = 0;                                          // message thread only
};

} // namespace felitronics::neural

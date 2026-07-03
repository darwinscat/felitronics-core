// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <algorithm>   // std::copy
#include <atomic>
#include <cstddef>

namespace felitronics::analysis
{

// The FFT analysis size is a shared contract between the audio-side tap and the GUI-side FFT.
constexpr int kSpectrumFftOrder = 11;                 // 2048-point default
constexpr int kSpectrumFftSize  = 1 << kSpectrumFftOrder;

//==============================================================================
// felitronics::analysis::SpectrumTap — one mono capture window with an atomic ready handshake (SPSC,
// lock-free). The audio thread push()es samples; once a window fills it snapshots into `data` (unless
// the reader hasn't consumed the previous frame yet) and flags `ready`; the GUI tryPull()s the latest
// frame and re-arms. RT-safe: fixed-size copies, no alloc / lock. JUCE-free (`std::copy`).
//
// Consolidates two diverged copies into one canonical primitive: cab::SpectrumTap (a JUCE-coupled
// copy using juce::FloatVectorOperations) and teq::SpectrumTap. Optimization over both: the
// per-sample-written producer cursor (`fifo`/`idx`) is cache-line separated (alignas) from the
// GUI-polled handshake (`ready`/`data`), so the audio thread's per-sample `idx` write never
// false-shares the cache line the GUI polls — a desktop hot-path win.
//
// Templated on the FFT order (default 2048) so a consumer can pick a different analysis resolution
// without forking; the `SpectrumTap` alias below is the drop-in default.
template <int FftOrder = kSpectrumFftOrder>
struct SpectrumTapT
{
    static constexpr int kOrder = FftOrder;
    static constexpr int kSize  = 1 << FftOrder;

    // Producer-only hot state (the audio thread writes `idx` every sample). NB: explicit size_t casts on
    // the template-dependent array bounds — GCC's -Wsign-conversion flags a dependent int bound.
    alignas (64) float fifo[(std::size_t) kSize] {};
    int idx = 0;

    // Cross-thread handshake (`ready`) + the snapshot the GUI reads (`data`), each on its own cache
    // line so they don't false-share the producer's per-sample `idx` write.
    alignas (64) std::atomic<bool> ready { false };
    alignas (64) float data[(std::size_t) kSize] {};

    void reset() noexcept { idx = 0; ready.store (false, std::memory_order_release); }

    // GUI thread: if a frame is ready, copy it into dst[kSize] and re-arm; else return false.
    // Preferred read path — never hands out the internal buffer.
    bool tryPull (float* dst) noexcept
    {
        if (! ready.load (std::memory_order_acquire)) return false;
        std::copy (data, data + kSize, dst);
        ready.store (false, std::memory_order_release);
        return true;
    }

    // Audio thread: push one sample. On a full window, snapshot into `data` + flag ready (unless the
    // reader hasn't consumed the previous frame), then wrap.
    void push (float s) noexcept
    {
        if (idx == kSize)
        {
            if (! ready.load (std::memory_order_acquire))
            {
                std::copy (fifo, fifo + kSize, data);
                ready.store (true, std::memory_order_release);
            }
            idx = 0;
        }
        fifo[idx++] = s;
    }
};

// Drop-in default (2048-point) — what both plugins use.
using SpectrumTap = SpectrumTapT<kSpectrumFftOrder>;

} // namespace felitronics::analysis

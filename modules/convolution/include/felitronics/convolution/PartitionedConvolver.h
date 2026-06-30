// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/core/Fft.h>

#include <algorithm>
#include <cstring>
#include <vector>

namespace felitronics::convolution
{

//==============================================================================
// felitronics::convolution::PartitionedConvolver — ZERO-LATENCY partitioned FFT convolution.
//
// IR = head (IR[0:P]) + tail (IR[P:M]). The head is convolved DIRECTLY in the time domain → immediate,
// sample-accurate output. The tail is a uniform-partitioned overlap-save FFT convolution whose inherent
// P-sample latency is EXACTLY absorbed by the head offset (the tail's earliest contribution is at delay
// P), so total added latency is 0. The partition size P is fixed at prepare(); an internal 2P frame
// absorbs ARBITRARY host block sizes (the partition step is decoupled from the host block).
//
// RT-safe: prepare()/setIr() allocate (message thread); process() does NO alloc/lock/throw. The FFT is a
// compile-time backend (the seam) — no vtable in the hot path.
//
// SPIKE SCOPE: uniform partitioning (not Gardner non-uniform) — adequate for short IRs (guitar cab /
// short reverb) on desktop, which is the funded target. Long-IR efficiency (growing tail partitions) is
// a later production optimization behind this same API. A live IR swap here is INSTANTANEOUS; a
// click-free production swap needs a crossfade (run old+new over ~10-20 ms).
template <core::fft::RealFftBackend Fft = core::fft::DefaultRealFft>
class PartitionedConvolver
{
public:
    // partitionSize P (pow2; FFT size = 2P). maxIrSamples sizes the partition arrays. Message thread.
    bool prepare (int partitionSize, int maxIrSamples)
    {
        if (! core::fft::isPow2 (partitionSize)) return false;
        P_ = partitionSize;
        N_ = 2 * P_;
        if (! fft_.prepare (N_)) return false;
        specF_    = Fft::spectrumFloats (N_);
        maxParts_ = (maxIrSamples > P_) ? ((maxIrSamples - P_ + P_ - 1) / P_) : 0;

        frame_.assign       ((std::size_t) N_, 0.0f);
        h0_.assign          ((std::size_t) P_, 0.0f);
        pendingTail_.assign ((std::size_t) P_, 0.0f);
        inputSpec_.assign   ((std::size_t) specF_, 0.0f);
        acc_.assign         ((std::size_t) specF_, 0.0f);
        ifftOut_.assign     ((std::size_t) N_, 0.0f);
        fdl_.assign         ((std::size_t) maxParts_ * (std::size_t) specF_, 0.0f);
        irSpec_.assign      ((std::size_t) maxParts_ * (std::size_t) specF_, 0.0f);

        numParts_ = 0;
        reset();
        return true;
    }

    void reset() noexcept
    {
        std::fill (frame_.begin(),       frame_.end(),       0.0f);
        std::fill (pendingTail_.begin(), pendingTail_.end(), 0.0f);
        std::fill (fdl_.begin(),         fdl_.end(),         0.0f);
        phase_  = 0;
        fdlPos_ = 0;
    }

    int  partitionSize() const noexcept { return P_; }
    static constexpr int latencySamples() noexcept { return 0; }

    // Set the IR (copies, builds the head + the partition spectra). Message thread (allocates a scratch).
    // SPIKE: resets the running state (instantaneous swap). Production: crossfade instead.
    void setIr (const float* ir, int len)
    {
        if (len < 0) len = 0;
        for (int i = 0; i < P_; ++i) h0_[(std::size_t) i] = (i < len) ? ir[i] : 0.0f;

        const int tailLen = (len > P_) ? (len - P_) : 0;
        int parts = (tailLen > 0) ? ((tailLen + P_ - 1) / P_) : 0;
        if (parts > maxParts_) parts = maxParts_;

        std::vector<float> part ((std::size_t) N_, 0.0f);
        for (int j = 0; j < parts; ++j)
        {
            std::fill (part.begin(), part.end(), 0.0f);
            for (int i = 0; i < P_; ++i)
            {
                const int src = P_ + j * P_ + i;
                part[(std::size_t) i] = (src < len) ? ir[src] : 0.0f;
            }
            fft_.forward (part.data(), &irSpec_[(std::size_t) j * (std::size_t) specF_]);
        }
        numParts_ = parts;
        reset();
    }

    // Audio thread. out[n] = (in * IR)[n] (causal linear convolution). RT-safe.
    void process (const float* in, float* out, int n) noexcept
    {
        for (int s = 0; s < n; ++s)
        {
            frame_[(std::size_t) (P_ + phase_)] = in[s];

            // direct head: sum_{i=0}^{P-1} h0[i] * x[t-i]  (reads current+previous chunk in `frame_`)
            const float* fr = &frame_[(std::size_t) (P_ + phase_)];
            float head = 0.0f;
            for (int i = 0; i < P_; ++i) head += h0_[(std::size_t) i] * fr[-i];

            out[s] = head + pendingTail_[(std::size_t) phase_];

            if (++phase_ == P_) { phase_ = 0; processChunk(); }
        }
    }

private:
    // A full P-sample chunk just arrived (frame_ = [prev chunk | current chunk]). FFT it, push to the
    // frequency-domain delay line, MAC against the IR tail partitions, IFFT → next chunk's pendingTail.
    void processChunk() noexcept
    {
        fft_.forward (frame_.data(), inputSpec_.data());

        if (numParts_ > 0)
        {
            std::memcpy (&fdl_[(std::size_t) fdlPos_ * (std::size_t) specF_], inputSpec_.data(),
                         (std::size_t) specF_ * sizeof (float));

            std::fill (acc_.begin(), acc_.end(), 0.0f);
            for (int j = 0; j < numParts_; ++j)
            {
                int idx = fdlPos_ - j;
                if (idx < 0) idx += numParts_;
                fft_.spectralMultiplyAdd (&fdl_[(std::size_t) idx * (std::size_t) specF_],
                                          &irSpec_[(std::size_t) j * (std::size_t) specF_], acc_.data());
            }

            fft_.inverse (acc_.data(), ifftOut_.data());
            for (int i = 0; i < P_; ++i) pendingTail_[(std::size_t) i] = ifftOut_[(std::size_t) (P_ + i)];  // overlap-save: last P

            if (++fdlPos_ >= numParts_) fdlPos_ = 0;
        }
        else
        {
            std::fill (pendingTail_.begin(), pendingTail_.end(), 0.0f);
        }

        for (int i = 0; i < P_; ++i) frame_[(std::size_t) i] = frame_[(std::size_t) (P_ + i)];  // current → previous
    }

    Fft fft_;
    int P_ = 0, N_ = 0, specF_ = 0, maxParts_ = 0, numParts_ = 0, phase_ = 0, fdlPos_ = 0;
    std::vector<float> frame_, h0_, pendingTail_, inputSpec_, acc_, ifftOut_, fdl_, irSpec_;
};

} // namespace felitronics::convolution

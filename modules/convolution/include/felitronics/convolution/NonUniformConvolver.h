// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/core/Fft.h>

#include <algorithm>
#include <array>
#include <cstring>
#include <vector>

namespace felitronics::convolution
{

//==============================================================================
// felitronics::convolution::NonUniformConvolver — MONO, TRUE sample-zero-latency, NON-UNIFORM partitioned FFT
// convolution (Gardner 1995). Phase-1 primitive: the proof that a non-uniform partition schedule undercuts the
// fixed-P=128 PartitionedConvolver / MatrixConvolver at every host block while staying block-INDEPENDENT and
// keeping true sample-zero-latency (unlike juce::dsp::Convolution, which is block-granular-zero-latency and only
// cheap when its P = the host block). Matrix routing (MSDiag/LRDiag/Full) + the click-free warm swap come in
// later phases on top of this; LRDiag is proven here as two independent instances.
//
// THE SCHEDULE. The IR is tiled by a HEAD + a list of overlap-save STAGES of GROWING block size:
//   head  [0, P0)                       — direct time-domain FIR, immediate, sample-accurate.
//   stage s: block B_s, count C_s        — a uniform-partitioned overlap-save conv of the IR slice
//            covers [offset_s, offset_s + C_s·B_s), FFT size 2·B_s, its own frequency-domain delay line (depth C_s).
// ZERO-LATENCY INVARIANT: offset_s == B_s for every stage — i.e. each stage begins exactly where the head + all
// earlier stages stopped, and that cumulative coverage equals the stage's own block size. This makes a stage's
// inherent B_s-sample overlap-save latency exactly absorbed by the earlier coverage (Gardner's condition), so
// TOTAL added latency is 0. The invariant is the recurrence  B_{s+1} = (C_s + 1)·B_s,  B_0 = P0. Pure octave
// doubling is the C_s≡1 special case (B doubles each stage); the efficient schedule uses a few DISTINCT block
// sizes each REPEATED (e.g. head 128 + 128×3 + 512×7 + 4096×31 for a 131072-tap IR) and caps the largest FFT.
//
// WHY NON-UNIFORM. A fixed small P pays O(L/P) spectral MACs (1023 partitions at P=128). JUCE's fixed large
// P=maxBlock pays few MACs but forces its partition = the host block (block-granular latency, cost explodes at
// small blocks). Non-uniform keeps small partitions ONLY in the near field (for zero latency) and large ones in
// the tail (for cheap MACs), independent of the host block. See docs/PERF-CONVOLVER-JUCE-GAP.md.
//
// COST NOTE (RT-safety, Phase-1 scope). Each stage fires its FFT synchronously when its block completes; stages
// of block B_s fire every B_s samples, so their firings COINCIDE every lcm(B_s) = B_max samples — a per-buffer
// CPU spike (the mean %RT hides it). Correctness/latency are unaffected (the work always finishes within the
// hop it is due). Spreading the large FFTs over their deadline (time-distributed scheduling) is a later-phase
// RT-hardening; capping B_max bounds the worst single FFT meanwhile.
//
// RT-safe: prepare()/setIr() allocate (message thread); process() does NO alloc/lock/throw, zero latency. The
// FFT is the compile-time backend seam — no vtable in the hot path. SPIKE SCOPE (like PartitionedConvolver): a
// live IR swap here is INSTANTANEOUS (setIr resets running state); the click-free crossfade is a later phase.
template <core::fft::RealFftBackend Fft = core::fft::DefaultRealFft>
class NonUniformConvolver
{
public:
    static constexpr int kMaxStages = 16;   // P0>=64, B_max<=2^21 → <=15 doubling stages + 1 tail; generous
    static constexpr int kMaxPartition = 1 << 24;   // hard cap so the FFT size 2*B never overflows signed int (far above any real IR)
    static constexpr int kMaxIrSamples = 1 << 24;   // hard cap so a cumulative stage offset (stored as int) can't overflow (16M smp ≫ any IR)

    // One tunable schedule step: `count` partitions of block size `blockSize` (pow2). The head size P0 and the
    // steps must satisfy the zero-latency recurrence B_{s+1} = (count_s + 1)·B_s with B_0 = P0 (validated in
    // prepare). A caller tunes the near-field/tail split for a platform; buildCappedSchedule() derives a
    // sensible default (doubling to B_max, then one uniform B_max stage covering the rest).
    struct ScheduleStep { int blockSize; int count; };

    // Default schedule: head P0 + doubling P0,2·P0,… (each ×1) until the block reaches maxBlock, then ONE uniform
    // maxBlock stage covering the remainder up to maxIrSamples. Returns the step count (<=kMaxStages), 0 if the
    // whole IR fits in the head. P0 and maxBlock must be pow2 with P0<=maxBlock.
    static int buildCappedSchedule (int P0, int maxBlock, int maxIrSamples,
                                    std::array<ScheduleStep, kMaxStages>& out) noexcept
    {
        int n = 0;
        long long cum = P0;                                   // head covers [0, P0)
        int B = P0;
        while (cum < (long long) maxIrSamples && n < kMaxStages)
        {
            if (B < maxBlock)
            {
                out[(std::size_t) n++] = { B, 1 };            // one doubling partition: covers [cum, cum+B)
                cum += B;
                B <<= 1;
            }
            else                                              // B == maxBlock: absorb the whole remainder here
            {
                const long long remaining = (long long) maxIrSamples - cum;
                int C = (int) ((remaining + B - 1) / B);      // ceil
                if (C < 1) C = 1;
                out[(std::size_t) n++] = { B, C };
                cum += (long long) C * B;
                break;
            }
        }
        return n;
    }

    // Convenience: prepare with the default capped-doubling schedule. headSize and maxBlock are pow2 with
    // maxBlock >= headSize (the doubling reaches maxBlock exactly, then caps). Leaves the engine UNPREPARED
    // on an invalid argument (house rule: a failed re-prepare must not keep a stale plan).
    bool prepare (int headSize, int maxBlock, int maxIrSamples)
    {
        prepared_ = false;
        if (! core::fft::isPow2 (headSize) || ! core::fft::isPow2 (maxBlock) || maxBlock < headSize) return false;
        if (headSize > kMaxPartition || maxBlock > kMaxPartition) return false;
        if (maxIrSamples > kMaxIrSamples) return false;
        std::array<ScheduleStep, kMaxStages> steps {};
        const int n = buildCappedSchedule (headSize, maxBlock, maxIrSamples, steps);
        return prepareWithSchedule (headSize, steps.data(), n, maxIrSamples);
    }

    // Full control: an explicit schedule for bench tuning. headSize P0 (pow2) is the time-domain head; `steps`
    // holds the growing overlap-save stages (each blockSize pow2). Validates the zero-latency recurrence
    // B_{s+1}=(count_s+1)·B_s with B_0=P0 and rejects any violation. maxIrSamples sizes the partition arrays.
    // Message thread; allocates. Returns false (leaving the engine UNPREPARED) on any invalid argument.
    bool prepareWithSchedule (int headSize, const ScheduleStep* steps, int numSteps, int maxIrSamples)
    {
        prepared_ = false;                                    // any early return below leaves it unprepared
        if (! core::fft::isPow2 (headSize) || headSize > kMaxPartition) return false;
        if (numSteps < 0 || numSteps > kMaxStages) return false;
        if (numSteps > 0 && steps == nullptr) return false;
        if (maxIrSamples > kMaxIrSamples) return false;
        if (maxIrSamples < 0) maxIrSamples = 0;

        // Validate the WHOLE schedule — recurrence AND coverage — BEFORE touching any state, so a bad schedule
        // can't half-build. cum = cumulative coverage = the required offset of the next stage.
        long long cum = headSize;
        for (int s = 0; s < numSteps; ++s)
        {
            const int B = steps[s].blockSize;
            const int C = steps[s].count;
            if (! core::fft::isPow2 (B) || C < 1 || B > kMaxPartition) return false;
            if ((long long) B != cum) return false;           // zero-latency invariant: offset_s == B_s
            cum += (long long) C * B;
        }
        if (cum < (long long) maxIrSamples) return false;     // the schedule MUST cover maxIrSamples — else setIr
                                                              // would silently drop every tap past the last stage

        P0_ = headSize;
        headMask_ = P0_ - 1;
        h0_.assign ((std::size_t) P0_, 0.0f);
        headHist_.assign ((std::size_t) P0_, 0.0f);
        headPos_ = 0;

        numStages_ = numSteps;
        long long off = P0_;
        for (int s = 0; s < numStages_; ++s)
        {
            if (! stages_[(std::size_t) s].prepare (steps[s].blockSize, steps[s].count, (int) off)) return false;
            off += (long long) steps[s].count * steps[s].blockSize;
        }
        numActiveStages_ = 0;                                 // set by setIr

        prepared_ = true;
        return true;
    }

    void reset() noexcept
    {
        std::fill (headHist_.begin(), headHist_.end(), 0.0f);
        headPos_ = 0;
        for (int s = 0; s < numStages_; ++s) stages_[(std::size_t) s].reset();
    }

    static constexpr int latencySamples() noexcept { return 0; }
    int  headSize()  const noexcept { return P0_; }
    int  numStages() const noexcept { return numStages_; }

    // Set the IR (copies, builds the head taps + every stage's partition spectra). Message thread (the stages'
    // build FFTs allocate a scratch). SPIKE: resets the running state (instantaneous swap). Production: crossfade.
    // NOT concurrency-safe with process(): setIr() reuses each stage's audio-thread fft_ and calls reset(), so it
    // must NOT run while process() is on another thread — this primitive is a "swap while stopped" API. For a
    // click-free LIVE swap on the audio thread use MatrixConvolverNupc (separate build FFT + a 2-slot crossfade).
    void setIr (const float* ir, int len)
    {
        if (! prepared_) return;
        if (len < 0) len = 0;
        for (int i = 0; i < P0_; ++i) h0_[(std::size_t) i] = (i < len) ? ir[i] : 0.0f;

        numActiveStages_ = 0;
        for (int s = 0; s < numStages_; ++s)
        {
            const bool anyActive = stages_[(std::size_t) s].setSlice (ir, len);
            if (anyActive) numActiveStages_ = s + 1;          // stages are offset-ordered: last active bounds the loop
        }
        reset();
    }

    // Audio thread. out[n] = (in * IR)[n] (causal linear convolution). RT-safe, zero latency. May alias in==out.
    void process (const float* in, float* out, int n) noexcept
    {
        if (! prepared_) return;
        for (int s = 0; s < n; ++s)
        {
            const float x = in[s];

            // Direct time-domain head: y = Σ_{i=0}^{P0-1} h0[i]·x[t-i] (immediate; reads a P0-deep ring).
            headPos_ = (headPos_ + 1) & headMask_;
            headHist_[(std::size_t) headPos_] = x;
            float y = 0.0f;
            for (int i = 0; i < P0_; ++i)
                y += h0_[(std::size_t) i] * headHist_[(std::size_t) ((headPos_ - i) & headMask_)];

            // Each active stage: feed the sample and add its pending overlap-save tail at the stage's phase.
            for (int st = 0; st < numActiveStages_; ++st)
                y += stages_[(std::size_t) st].tapAndTail (x);

            out[s] = y;

            // Advance stage phases AFTER the output is formed (a chunk boundary recomputes the NEXT tail).
            for (int st = 0; st < numActiveStages_; ++st)
                stages_[(std::size_t) st].advance();
        }
    }

private:
    //==========================================================================
    // One overlap-save STAGE: a uniform-partitioned convolution of the IR slice [offset, offset + C·B), block
    // size B, FFT size 2B, own frequency-domain delay line (depth C). It is the PartitionedConvolver tail
    // machinery WITHOUT a head — the head + earlier stages already cover [0, offset). Since offset == B (the
    // zero-latency invariant), partition j covers IR[offset + jB : offset + (j+1)B) and is MAC'd against the
    // input spectrum from j chunks ago; the earliest (j=0) contributes at delay B — exactly one overlap-save hop.
    struct Stage
    {
        Fft fft_;                                             // per-stage FFT (size 2B; each stage a distinct size)
        core::fft::AlignedVector<float> frame_;               // 2B: [prev B | cur B] input accumulator (seam)
        core::fft::AlignedVector<float> fdl_;                 // C × specF: ring of past input spectra (seam)
        core::fft::AlignedVector<float> irSpec_;              // C × specF: the slice's partition spectra (seam)
        core::fft::AlignedVector<float> inputSpec_, acc_, ifftOut_, buildBuf_;   // scratch (seam-crossing)
        std::vector<float> pendingTail_;                      // B: tail for the next B output samples (time domain)
        int B_ = 0, C_ = 0, specF_ = 0, offset_ = 0;
        int activeParts_ = 0, phase_ = 0, fdlPos_ = 0;

        bool prepare (int B, int C, int offset) noexcept
        {
            if (! core::fft::isPow2 (B) || C < 1) return false;
            B_ = B; C_ = C; offset_ = offset;
            const int N = 2 * B_;
            if (! fft_.prepare (N)) return false;
            specF_ = Fft::spectrumFloats (N);

            frame_.assign     ((std::size_t) N, 0.0f);
            inputSpec_.assign ((std::size_t) specF_, 0.0f);
            acc_.assign       ((std::size_t) specF_, 0.0f);
            ifftOut_.assign   ((std::size_t) N, 0.0f);
            buildBuf_.assign  ((std::size_t) N, 0.0f);
            fdl_.assign       ((std::size_t) C_ * (std::size_t) specF_, 0.0f);
            irSpec_.assign    ((std::size_t) C_ * (std::size_t) specF_, 0.0f);
            pendingTail_.assign ((std::size_t) B_, 0.0f);
            activeParts_ = 0; phase_ = 0; fdlPos_ = 0;
            return true;
        }

        // Build this stage's partition spectra from the IR slice [offset, offset + C·B). A partition is ACTIVE
        // only if its slice start is within `len`; returns true if any partition is active (the stage does work).
        bool setSlice (const float* ir, int len) noexcept
        {
            int parts = 0;
            for (int j = 0; j < C_; ++j)
            {
                const int base = offset_ + j * B_;            // absolute IR index of this partition's first tap
                if (base >= len) break;                       // this and all later partitions are all-zero
                std::fill (buildBuf_.begin(), buildBuf_.end(), 0.0f);
                for (int i = 0; i < B_; ++i)
                {
                    const int src = base + i;
                    buildBuf_[(std::size_t) i] = (src < len) ? ir[src] : 0.0f;   // partition in [0,B); [B,2B) stays 0
                }
                fft_.forward (buildBuf_.data(), &irSpec_[(std::size_t) j * (std::size_t) specF_]);
                parts = j + 1;
            }
            activeParts_ = parts;
            return parts > 0;
        }

        void reset() noexcept
        {
            std::fill (frame_.begin(),       frame_.end(),       0.0f);
            std::fill (fdl_.begin(),         fdl_.end(),         0.0f);
            std::fill (pendingTail_.begin(), pendingTail_.end(), 0.0f);
            phase_ = 0; fdlPos_ = 0;
        }

        // Feed one input sample into the current frame half and return the tail sample due at this phase.
        float tapAndTail (float x) noexcept
        {
            frame_[(std::size_t) (B_ + phase_)] = x;
            return pendingTail_[(std::size_t) phase_];
        }

        // Advance the phase; on a chunk boundary FFT the frame, MAC the partitions, IFFT into the next tail.
        void advance() noexcept
        {
            if (++phase_ != B_) return;
            phase_ = 0;
            processChunk();
        }

        void processChunk() noexcept
        {
            fft_.forward (frame_.data(), inputSpec_.data());  // 2B frame → spectrum
            if (activeParts_ > 0)                             // (an all-inactive stage is never in the processed
            {                                                 //  range, but guard the MAC/IFFT like PartitionedConvolver)
                std::memcpy (&fdl_[(std::size_t) fdlPos_ * (std::size_t) specF_], inputSpec_.data(),
                             (std::size_t) specF_ * sizeof (float));

                std::fill (acc_.begin(), acc_.end(), 0.0f);
                for (int j = 0; j < activeParts_; ++j)
                {
                    int idx = fdlPos_ - j; if (idx < 0) idx += C_;
                    fft_.spectralMultiplyAdd (&fdl_[(std::size_t) idx * (std::size_t) specF_],
                                              &irSpec_[(std::size_t) j * (std::size_t) specF_], acc_.data());
                }
                fft_.inverse (acc_.data(), ifftOut_.data());
                for (int i = 0; i < B_; ++i) pendingTail_[(std::size_t) i] = ifftOut_[(std::size_t) (B_ + i)];   // overlap-save: last B

                if (++fdlPos_ >= C_) fdlPos_ = 0;
            }
            for (int i = 0; i < B_; ++i) frame_[(std::size_t) i] = frame_[(std::size_t) (B_ + i)];   // cur → prev
        }
    };

    std::array<Stage, kMaxStages> stages_ {};
    std::vector<float> h0_, headHist_;                        // P0 each: head taps + input ring (time domain)
    int P0_ = 0, headMask_ = 0, headPos_ = 0;
    int numStages_ = 0, numActiveStages_ = 0;
    bool prepared_ = false;
};

} // namespace felitronics::convolution

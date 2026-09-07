// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/core/Config.h>
#include <felitronics/mastering/MasteringChain.h>

#include <algorithm>
#include <vector>

namespace felitronics::mastering
{

//==============================================================================
// felitronics::mastering::OfflineRenderer — the file-render wrapper over the streaming chain: run to
// the end, drain the tail, cut the latency back off. A THIN layer, deliberately: the chain is the
// object, and a live preview drives the same one. An offline monolith with the loop inside would have
// to be taken apart the first time a plug-in wanted it.
//
// THE CONTRACT IS A FORMULA, not a description. Let `y` be what the chain emits for the input `x`
// followed by `D = chain.latencySamples()` zeros. Then
//
//     out[n] = y[n + D]   for 0 <= n < frames
//
// so the output has exactly as many frames as the input, sample n of the output is sample n of the
// input processed, and the last `D` frames — the ones that only leave the chain after the input has
// ended — are in it. That last clause is the whole reason this class exists: the ffmpeg chain this
// replaces never emitted them, measured, so a click 4 ms before the end of a file disappeared and one
// 5 ms before it survived. Restating that as arithmetic is what makes it testable rather than hoped
// for; a suite can compute the right-hand side independently and null against it.
//
// It works with frames < D (the whole output then comes out of the drain), with frames == 0, and with
// `in == out` — every block is copied into internal scratch before anything is written back, and the
// write cursor trails the read cursor by D in any case.
//
// EVERY RENDER STARTS FROM A FULL RESET, including the dither's RNG, so two renders of the same input
// with the same parameters are bit-identical. P7 needs that: it renders, measures, adjusts one gain
// and renders again, and a second pass that started from the first pass's tail state would make its
// measurement mean something else. NB it also means a render is not a continuation — the chain is a
// streaming object and this is a whole-programme operation on it.
//
// RT: this is an OFFLINE class and does not pretend otherwise — it is called from a worker, not an
// audio callback. It does not allocate inside `render()` (the scratch is sized in `prepare()`), which
// matters because the chain underneath is RT-safe and a test that counts allocations over a whole
// render should see none.
class OfflineRenderer
{
public:
    // `blockSize` is the renderer's own read/write granularity and has NO effect on the result — the
    // chain re-blocks everything to its internal quantum anyway, which is exactly what makes this
    // parameter free to choose. Pinned in the suite: rendering the same programme at 1, 977 and 65536
    // is bit-identical.
    [[nodiscard]] bool prepare (int maxChannels, int blockSize)
    {
        // LAW 11(b), AND ITS MISSING HALF: a refused prepare() writes NOTHING **and leaves the object
        // unusable**. The two are only compatible in this order — DISARM first, validate, then write —
        // because validating first means returning before the disarm, which leaves the previous
        // preparation standing and answering process() calls. This is Saturator's idiom, made general.
        block_ = 0;                              // 0 == unprepared; render() refuses on it Validating one argument,
        // storing it, and then refusing on the next left the new WIDTH standing beside the old buffer —
        // and render() sizes its scratch pointers from the width. ASan: heap-buffer-overflow, a WRITE
        // four bytes past a 1024-byte region, after prepare(1, 256) then prepare(2, 0).
        if (maxChannels < 1 || maxChannels > core::kMaxChannels) return false;
        if (blockSize < 1) return false;
        maxCh_ = maxChannels;
        block_ = blockSize;
        scratch_.assign ((std::size_t) maxCh_ * (std::size_t) block_, 0.0f);
        return true;
    }

    int maxChannels() const noexcept { return maxCh_; }
    int blockSize()   const noexcept { return block_; }

    // Render `frames` of `in` into `out`. Returns false — having written nothing — if the chain is not
    // prepared, if the channel count is not the chain's, or if this renderer was not prepared for that
    // many channels. `in` and `out` may be the same buffers.
    bool render (MasteringChain& chain, const float* const* in, float* const* out,
                 int numChannels, int frames)
    {
        if (block_ < 1) return false;                   // a refused prepare() leaves it unusable
        if (! chain.isPrepared() || numChannels != chain.numChannels()) return false;
        if (numChannels < 1 || numChannels > maxCh_ || scratch_.empty()) return false;
        if (frames < 0) return false;
        if (frames > 0 && (in == nullptr || out == nullptr)) return false;

        chain.reset();
        const long long D     = chain.latencySamples();
        const long long total = (long long) frames + D;

        float* sp[core::kMaxChannels] {};
        for (int c = 0; c < numChannels; ++c)
            sp[c] = scratch_.data() + (std::size_t) c * (std::size_t) block_;

        for (long long off = 0; off < total; )
        {
            const int m = (int) std::min<long long> ((long long) block_, total - off);

            // Read the whole slice into scratch BEFORE writing anything back, so `in == out` is safe
            // even at D == 0. Past the end of the input the chain is fed zeros — which is what "flush"
            // means here, and why there is no separate tail code path to get wrong.
            for (int c = 0; c < numChannels; ++c)
                for (int i = 0; i < m; ++i)
                {
                    const long long s = off + i;
                    sp[c][i] = (s < (long long) frames) ? in[c][s] : 0.0f;
                }

            if (! chain.process (sp, numChannels, m)) return false;

            // Streaming sample (off + i) carries input sample (off + i - D). Everything before 0 is the
            // chain's own priming and is dropped; everything from `frames` on is past the end.
            for (int c = 0; c < numChannels; ++c)
                for (int i = 0; i < m; ++i)
                {
                    const long long o = off + i - D;
                    if (o >= 0 && o < (long long) frames) out[c][o] = sp[c][i];
                }

            off += m;
        }
        return true;
    }

private:
    std::vector<float> scratch_;
    int maxCh_ = 0, block_ = 0;
};

} // namespace felitronics::mastering

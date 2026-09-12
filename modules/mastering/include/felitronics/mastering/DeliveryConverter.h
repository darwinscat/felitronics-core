// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/core/Config.h>
#include <felitronics/core/DeliveryResampler.h>

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <vector>

namespace felitronics::mastering
{

//==============================================================================
// felitronics::mastering::DeliveryConverter — SRC FIRST, as one whole-programme operation. It turns a
// programme at the SOURCE rate into the same programme at the DELIVERY rate, so that everything after it
// — `OfflineRenderer::render`, `TargetLoudnessSolver::solve`, the range measurement — runs unchanged on a
// chain prepared at the delivery rate, and the dither stays the last thing to touch the samples.
//
// WHY SRC IS FIRST, and not up for debate: resampling AFTER the true-peak limiter raises inter-sample
// peaks the limiter never saw, so -1 dBTP stops being a promise exactly at the conversion. Everything
// that measures or limits has to see the delivered rate.
//
// WHY THIS IS A CORE CLASS AND NOT A FEW LINES IN THE C ABI FACADE. Two paths render a programme — the C
// ABI and the direct C++ call `fcore_master`'s selftest nulls it against, bit for bit — and the delicate
// part of a delivered render is the arithmetic around the converter: the exact length, the latency trim,
// the drain. Written twice by hand, the two copies can agree with each other while both are wrong. So
// the arithmetic lives here once and both paths call it.
//
// THE CONTRACT IS A FORMULA. For a programme of `inFrames` frames at `inRate`:
//
//     outFrames = ceil (inFrames * deliveryRate / inRate)            in EXACT integer arithmetic
//     out[n]    = y[T0 + n]      T0 = round (D),  D = the converter's exact rational latency
//
// where y is everything the converter emits for the input followed by as much silence as it takes. So
// output sample n carries input time (n - d) / deliveryRate with d = D - T0 in [-0.5, 0.5]: a constant
// sub-sample offset of the whole file, stated rather than hidden, because an integer trim of a rational
// delay cannot remove it.
//
// 🔴 THE LENGTH IS NOT `std::ceil` OF A DOUBLE. 147 frames at 44.1 -> 48 kHz is exactly 160 frames, and
// `std::ceil (147.0 * 48000.0 / 44100.0)` is 161, because the quotient is 160.00000000000003 in double.
// It comes from the reduced ratio L:M instead: (inFrames * L + M - 1) / M.
//
// 🔴 THE DRAIN IS FED UNTIL IT IS ENOUGH, not computed once. `DeliveryResampler::flush()` pushes the
// cascade's delay rounded up, which for a delivered programme can still be one output short of
// T0 + outFrames (176,401 frames at 176.4 -> 44.1 kHz: 44,220 emitted, 44,221 needed). So `convert()`
// never flushes; it feeds silence in blocks until the converter has produced what the formula needs.
//
// EQUAL RATES ARE A COPY: outFrames == inFrames, D = 0, and the samples are the caller's bits.
//==============================================================================
class DeliveryConverter
{
public:
    // The exact delivered length. -1 for rates the resampler would refuse, or a result past `long long`.
    [[nodiscard]] static long long deliveredFrames (double inRate, double deliveryRate, long long inFrames) noexcept
    {
        if (inFrames < 0) return -1;
        core::DeliveryResampler::Params p;
        p.inRate = inRate; p.outRate = deliveryRate;
        const auto pl = core::DeliveryResampler::plan (p);
        if (! pl.ok) return -1;
        if (pl.identity) return inFrames;
        const long long a = (long long) std::llround (inRate), b = (long long) std::llround (deliveryRate);
        const long long g = std::gcd (a, b);
        const long long L = b / g, M = a / g;
        if (inFrames > (LLONG_MAX - M) / L) return -1;
        return (inFrames * L + M - 1) / M;
    }

    // WHAT prepare() ASKS THE HEAP FOR (law 11d): the resampler plus one block of output staging per channel.
    struct Storage
    {
        bool ok = false;
        std::uint64_t resampler = 0, staging = 0;
        std::uint64_t bytes() const noexcept { return ok ? resampler + staging : 0u; }
    };
    [[nodiscard]] static Storage storageFor (double inRate, double deliveryRate, int numChannels, int block) noexcept
    {
        Storage st;
        core::DeliveryResampler::Params p;
        p.inRate = inRate; p.outRate = deliveryRate;
        const auto rs = core::DeliveryResampler::storageFor (p, numChannels, block);
        if (! rs.ok) return st;
        const long long perCall = outputBoundFor (rs.plan, block);
        if (perCall <= 0 || perCall > (1 << 26)) return st;
        st.resampler = rs.bytes();
        st.staging   = (std::uint64_t) sizeof (float) * (std::uint64_t) numChannels * (std::uint64_t) perCall
                     + (std::uint64_t) sizeof (float) * (std::uint64_t) numChannels * (std::uint64_t) block;   // silence
        st.ok = true;
        return st;
    }
    [[nodiscard]] static std::uint64_t prepareBytes (double inRate, double deliveryRate, int numChannels, int block) noexcept
    {
        return storageFor (inRate, deliveryRate, numChannels, block).bytes();
    }

    [[nodiscard]] bool prepare (double inRate, double deliveryRate, int numChannels, int block)
    {
        prepared_ = false;
        const Storage st = storageFor (inRate, deliveryRate, numChannels, block);
        if (! st.ok) return false;
        core::DeliveryResampler::Params p;
        p.inRate = inRate; p.outRate = deliveryRate;
        if (! src_.prepare (p, numChannels, block)) return false;
        inRate_ = inRate; deliveryRate_ = deliveryRate; nch_ = numChannels; block_ = block;
        perCall_ = (int) outputBoundFor (src_.currentPlan(), block);
        staging_.assign ((std::size_t) numChannels * (std::size_t) perCall_, 0.0f);
        silence_.assign ((std::size_t) numChannels * (std::size_t) block, 0.0f);
        prepared_ = true;
        return true;
    }

    bool isPrepared() const noexcept { return prepared_; }
    double latencyOutputSamples() const noexcept { return src_.latencyOutputSamples(); }
    long long trimSamples() const noexcept { return std::llround (src_.latencyOutputSamples()); }
    const core::DeliveryResampler::Plan& plan() const noexcept { return src_.currentPlan(); }

    // Convert a whole programme. `out` must hold exactly deliveredFrames(inFrames) frames per channel, and
    // `outFrames` must be that number — a caller that computed its own is refused, not trusted. Every call
    // starts from a reset, so two conversions of the same programme are bit-identical. No allocation.
    [[nodiscard]] bool convert (const float* const* in, int numChannels, long long inFrames,
                                float* const* out, long long outFrames) noexcept
    {
        if (! prepared_ || numChannels != nch_ || inFrames < 0) return false;
        if (outFrames != deliveredFrames (inRate_, deliveryRate_, inFrames)) return false;
        if (outFrames > 0 && (in == nullptr || out == nullptr)) return false;
        for (int c = 0; c < numChannels && outFrames > 0; ++c)
            if (in[c] == nullptr || out[c] == nullptr) return false;
        if (outFrames == 0) return true;

        src_.reset();
        const long long T0 = src_.currentPlan().identity ? 0 : trimSamples();
        long long emitted = 0;                       // converter outputs seen so far
        auto take = [&] (int got) noexcept
        {
            for (int k = 0; k < got; ++k, ++emitted)
            {
                const long long n = emitted - T0;
                if (n < 0 || n >= outFrames) continue;
                for (int c = 0; c < numChannels; ++c)
                    out[c][n] = staging_[(std::size_t) c * (std::size_t) perCall_ + (std::size_t) k];
            }
        };
        float* sp[core::kMaxChannels] {};
        for (int c = 0; c < numChannels; ++c) sp[c] = staging_.data() + (std::size_t) c * (std::size_t) perCall_;

        for (long long off = 0; off < inFrames; off += block_)
        {
            const int m = (int) std::min<long long> (block_, inFrames - off);
            const float* ip[core::kMaxChannels] {};
            for (int c = 0; c < numChannels; ++c) ip[c] = in[c] + off;
            int got = 0;
            if (! src_.process (ip, numChannels, m, sp, perCall_, got)) return false;
            take (got);
        }
        const float* zp[core::kMaxChannels] {};
        for (int c = 0; c < numChannels; ++c) zp[c] = silence_.data() + (std::size_t) c * (std::size_t) block_;
        // The drain: silence until T0 + outFrames outputs exist. Bounded — every block of silence yields
        // at least block*L/M outputs, and the shortfall after the programme is about one latency.
        for (long long guard = 0; emitted < T0 + outFrames; ++guard)
        {
            if (guard > 64 + (T0 + outFrames) / std::max (1, block_)) return false;
            int got = 0;
            if (! src_.process (zp, numChannels, block_, sp, perCall_, got)) return false;
            take (got);
        }
        return true;
    }

private:
    static long long outputBoundFor (const core::DeliveryResampler::Plan& pl, long long n) noexcept
    {
        if (pl.identity) return n;
        long long b = n;
        for (int s = 0; s < pl.count; ++s) b = (b * pl.stage[s].L) / pl.stage[s].M + 1;
        return b;
    }

    core::DeliveryResampler src_;
    std::vector<float> staging_, silence_;
    double inRate_ = 0.0, deliveryRate_ = 0.0;
    int nch_ = 0, block_ = 0, perCall_ = 0;
    bool prepared_ = false;
};

} // namespace felitronics::mastering

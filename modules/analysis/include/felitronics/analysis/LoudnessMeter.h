// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/core/Config.h>
#include <felitronics/analysis/KWeightingFilter.h>

#include <algorithm>
#include <cmath>
#include <vector>

namespace felitronics::analysis
{

//==============================================================================
// felitronics::analysis::LoudnessMeter — ITU-R BS.1770-4 / EBU R128 loudness: MOMENTARY (400 ms),
// SHORT-TERM (3 s) and INTEGRATED (gated). LUFS = -0.691 + 10·log10(Σ weightₖ·meanSquareₖ) of the
// K-weighted signal.
//
// Energy is accumulated in 10 ms SUB-HOPS: momentary is the mean of the last 40, short-term of the last 300,
// so either window lands within 10 ms of any event. EBU Tech 3341's file-based cases 10 and 13 slide a 3 s /
// 400 ms tone in 150 ms / 20 ms steps and expect the maximum reading to be the tone ±0.1 LU at EVERY offset —
// a window that only moves in 100 ms steps cannot do that (a 400 ms burst 40 ms off its grid reads 0.45 LU
// low). The integrated measure keeps its 100 ms HOP: every tenth sub-hop closes a 400 ms gating block
// (75 % overlap), and integrated = the gated mean over those blocks — absolute gate at -70 LUFS, then a
// -10 LU relative gate (a two-pass over all absolute-gated blocks — the threshold moves as more program
// arrives, so it must NOT be a one-pass running sum).
//
// RT-safe: prepare() allocates the sub-hop ring + the integrated-block buffer; process() only indexes them
// (no alloc/lock/throw). The gated measure keeps every block's energy to the end (the relative gate is a
// two-pass), so the block store is sized for maxDurationSec of program at the prepared rate; blocks past it
// are counted in droppedBlocks() and not kept — size for the longest program and that reads 0. Channel
// weights default to 1.0 (correct for mono/stereo); set them per the BS.1770 roles (Ls/Rs = 1.41, LFE
// excluded) for surround — the host-layout→role mapping is product glue.
class LoudnessMeter
{
public:
    void prepare (double sampleRate, int numChannels, double maxDurationSec = 3600.0)
    {
        prepared_ = false;
        fs = sampleRate > 0.0 ? sampleRate : 48000.0;                  // fs<=0 → subSamples 0 → /0 in finishSubHop
        ch = numChannels < 1 ? 1 : (numChannels > kMaxChannels ? kMaxChannels : numChannels);
        kw.prepare (fs, ch);
        subSamples = std::max (1, (int) std::lround (0.01 * fs));      // 10 ms; a gating hop is ten of them
        for (int c = 0; c < kMaxChannels; ++c) w[c] = 1.0;
        subRing.assign (kSubRing, 0.0);
        // Sized by HOPS at this rate, not by seconds: a hop is 10 × lround (0.01·fs) samples, which is 100 ms
        // only where fs is a multiple of 100 — elsewhere a per-second count drifts from the blocks that
        // actually arrive. +4 blocks / +8 samples of slack cover the first block's onset and the 1 s cadence.
        const double hops = std::ceil (std::max (0.0, maxDurationSec) * fs) / (double) (subSamples * kSubHopsPerHop);
        blockE.assign ((std::size_t) hops + 4, 0.0);
        stE.assign ((std::size_t) (hops / 10.0) + 8, 0.0);              // 1 short-term sample/s for LRA
        reset();
        prepared_ = true;
    }

    void reset() noexcept
    {
        kw.reset();
        for (int c = 0; c < kMaxChannels; ++c) subSumSq[c] = 0.0;
        subCount = 0; subWrite = 0; subFilled = 0; subInHop = 0; blockCount = 0; droppedBlocks_ = 0;
        stCount = 0; stSince = 0;
        std::fill (subRing.begin(), subRing.end(), 0.0);
    }

    void setChannelWeight (int c, double weight) noexcept { if (c >= 0 && c < kMaxChannels) w[c] = weight; }

    void process (const float* const* channels, int numChannels, int n) noexcept
    {
        if (! prepared_) return;                                       // unprepared — subRing/blockE/stE empty
        const int nc = numChannels < ch ? numChannels : ch;
        for (int i = 0; i < n; ++i)
        {
            for (int c = 0; c < nc; ++c) { const double y = kw.process (c, (double) channels[c][i]); subSumSq[c] += y * y; }
            if (++subCount >= subSamples) finishSubHop (nc);
        }
    }

    double momentaryLufs()  const noexcept { return lufsOf (meanLastSubHops (kMomentarySubHops)); }  // 400 ms
    double shortTermLufs()  const noexcept { return lufsOf (meanLastSubHops (kSubRing)); }           // 3 s
    double integratedLufs() const noexcept { return integrated(); }                                  // gated
    double loudnessRangeLu() const noexcept { return lra(); }                                        // EBU Tech 3342 (P95−P10)

    // Gating blocks that arrived past maxDurationSec and were not kept. Non-zero means integratedLufs() and
    // loudnessRangeLu() describe the first maxDurationSec of the program only — a caller that must not lose a
    // block sizes prepare() for its longest program and checks this reads 0.
    int droppedBlocks() const noexcept { return droppedBlocks_; }

    // The 400 ms gating blocks' K-weighted mean-square energies, in arrival order, exactly as finishHop()
    // recorded them — BEFORE either gate. This is what a cross-toolchain bit-exactness check compares:
    // integratedLufs() is DISCONTINUOUS in these energies (a block landing within ~1e-12 of a gate flips its
    // inclusion and moves the reading by ~0.01 dB), so it cannot carry a bit-identity claim, while the vector
    // itself is continuous in the input samples and can. Raw energy, not dB, deliberately — a dB accessor
    // would route the comparison back through log10, whose 1-ulp disagreement between libms is precisely what
    // such a check must not inherit. Valid for gatingBlockCount() entries until the next
    // process()/reset()/prepare(); the count is 0 (and the pointer may be null) before prepare().
    int           gatingBlockCount()    const noexcept { return blockCount; }
    const double* gatingBlockEnergies() const noexcept { return blockE.data(); }

private:
    static constexpr int kMaxChannels      = core::kMaxChannels;
    static constexpr int kSubHopsPerHop    = 10;                                    // 10 × 10 ms = the 100 ms gating hop
    static constexpr int kMomentarySubHops = 40;                                    // 400 ms
    static constexpr int kSubRing          = 300;                                   // 3 s — the short-term window, and the ring

    void finishSubHop (int nc) noexcept
    {
        double subMS = 0.0;
        for (int c = 0; c < nc; ++c) subMS += w[c] * (subSumSq[c] / (double) subSamples);
        subRing[(std::size_t) subWrite] = subMS;
        subWrite = (subWrite + 1) % kSubRing;
        if (subFilled < kSubRing) ++subFilled;
        if (++subInHop >= kSubHopsPerHop) { subInHop = 0; finishHop(); }
        // Clear EVERY channel's accumulator, not just c < nc: a channel that vanishes mid-hop (host drops
        // the channel count) must not park its partial energy and leak it into a later hop when it returns.
        for (int c = 0; c < kMaxChannels; ++c) subSumSq[c] = 0.0;
        subCount = 0;
    }

    // Every 100 ms: a 400 ms gating block for the integrated measure and, once 3 s are in, a short-term
    // sample every 1 s for LRA.
    void finishHop() noexcept
    {
        if (subFilled >= kMomentarySubHops)                                         // a 400 ms block every 100 ms
        {
            if (blockCount < (int) blockE.size()) blockE[(std::size_t) blockCount++] = meanLastSubHops (kMomentarySubHops);
            else ++droppedBlocks_;                                                  // past maxDurationSec: counted, not kept
        }
        if (subFilled >= kSubRing)                                                  // a 3 s short-term sample every 1 s (LRA)
        {
            if (stSince == 0 && stCount < (int) stE.size()) stE[(std::size_t) stCount++] = meanLastSubHops (kSubRing);
            if (++stSince >= 10) stSince = 0;                                        // first at 3 s, then every 10 hops (libebur128)
        }
    }

    double meanLastSubHops (int k) const noexcept
    {
        const int kk = k < subFilled ? k : subFilled;
        if (kk <= 0) return 0.0;
        double s = 0.0;
        for (int j = 0; j < kk; ++j) { const int idx = (subWrite - 1 - j + kSubRing) % kSubRing; s += subRing[(std::size_t) idx]; }
        return s / kk;
    }

    static double lufsOf (double meanSquare) noexcept { return meanSquare > 1e-12 ? -0.691 + 10.0 * std::log10 (meanSquare) : -120.0; }

    double integrated() const noexcept
    {
        if (blockCount <= 0) return -120.0;
        const double absT = std::pow (10.0, (-70.0 + 0.691) / 10.0);               // energy for -70 LUFS
        double sum = 0.0; int cnt = 0;
        for (int j = 0; j < blockCount; ++j) if (blockE[(std::size_t) j] > absT) { sum += blockE[(std::size_t) j]; ++cnt; }
        if (cnt == 0) return -120.0;
        const double relT = 0.1 * (sum / cnt);                                     // -10 LU relative to the abs-gated mean
        double s2 = 0.0; int c2 = 0;
        for (int j = 0; j < blockCount; ++j) { const double z = blockE[(std::size_t) j]; if (z > absT && z > relT) { s2 += z; ++c2; } }
        return c2 > 0 ? lufsOf (s2 / c2) : -120.0;
    }

    // LRA (EBU Tech 3342) = P95 − P10 of the gated 3 s short-term loudness distribution. Two-pass gate over
    // the 1 s-cadence stE[] energies: absolute −70 LUFS, then −20 LU below the energy-mean of the abs-gated
    // set; percentiles via a fixed 0.1 LU histogram (−70..+30 LUFS), libebur128-faithful (no sort, no alloc).
    double lra() const noexcept
    {
        if (stCount <= 0) return 0.0;
        const double absT = std::pow (10.0, (-70.0 + 0.691) / 10.0);                // energy for −70 LUFS
        double sum = 0.0; int cnt = 0;
        for (int j = 0; j < stCount; ++j) if (stE[(std::size_t) j] >= absT) { sum += stE[(std::size_t) j]; ++cnt; }
        if (cnt == 0) return 0.0;
        const double relT = 0.01 * (sum / cnt);                                     // −20 LU relative (energy ×0.01)

        constexpr int kBins = 1000;                                                 // −70..+30 LUFS, 0.1 LU bins
        int hist[kBins] = { 0 }; int total = 0;
        for (int j = 0; j < stCount; ++j)
        {
            const double e = stE[(std::size_t) j];
            if (e >= absT && e >= relT)                                             // ≥ gates (libebur128-faithful)
            {
                int b = (int) ((lufsOf (e) + 70.0) * 10.0);                         // 0.1 LU bins from −70 LUFS
                b = b < 0 ? 0 : (b >= kBins ? kBins - 1 : b);
                ++hist[b]; ++total;
            }
        }
        if (total < 2) return 0.0;
        auto pct = [&] (double p) {                                                 // libebur128 rank rule, bin lower bound
            const int rank = (int) ((double) (total - 1) * p + 0.5);
            int cum = 0, b = 0;
            for (; b < kBins; ++b) { cum += hist[b]; if (cum > rank) break; }
            return -70.0 + (double) (b < kBins ? b : kBins - 1) * 0.1;
        };
        return pct (0.95) - pct (0.10);                                             // the 0.05-LU bin offset cancels
    }

    double fs = 48000.0; int ch = 2, subSamples = 480;
    bool prepared_ = false;                     // true only after prepare() (subRing/blockE/stE allocated)
    KWeightingFilter kw;
    double w[kMaxChannels] {};
    double subSumSq[kMaxChannels] {};
    int subCount = 0;
    std::vector<double> subRing;                                                    // 300 × 10 ms sub-hop energies
    int subWrite = 0, subFilled = 0, subInHop = 0;
    std::vector<double> blockE;
    int blockCount = 0, droppedBlocks_ = 0;
    std::vector<double> stE;                                                        // 3 s short-term energies @1 s (LRA)
    int stCount = 0, stSince = 0;
};

} // namespace felitronics::analysis

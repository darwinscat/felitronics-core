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
// K-weighted signal. 100 ms hops; integrated = the gated mean over 400 ms blocks (75 % overlap):
// absolute gate at -70 LUFS, then a -10 LU relative gate (a two-pass over all absolute-gated blocks —
// the threshold moves as more program arrives, so it must NOT be a one-pass running sum).
//
// RT-safe: prepare() allocates the hop ring + the integrated-block buffer; process() only indexes them
// (no alloc/lock/throw). Channel weights default to 1.0 (correct for mono/stereo); set them per the
// BS.1770 roles (Ls/Rs = 1.41, LFE excluded) for surround — the host-layout→role mapping is product glue.
class LoudnessMeter
{
public:
    void prepare (double sampleRate, int numChannels, double maxDurationSec = 3600.0)
    {
        fs = sampleRate;
        ch = numChannels < 1 ? 1 : (numChannels > kMaxChannels ? kMaxChannels : numChannels);
        kw.prepare (fs, ch);
        hopSamples = (int) std::lround (0.1 * fs);                      // 100 ms
        for (int c = 0; c < kMaxChannels; ++c) w[c] = 1.0;
        hopRing.assign (kHopRing, 0.0);
        blockE.assign ((std::size_t) ((int) std::ceil (maxDurationSec * 10.0) + 4), 0.0);
        stE.assign ((std::size_t) ((int) std::ceil (maxDurationSec) + 8), 0.0);     // 1 short-term sample/s for LRA
        reset();
    }

    void reset() noexcept
    {
        kw.reset();
        for (int c = 0; c < kMaxChannels; ++c) hopSumSq[c] = 0.0;
        hopCount = 0; hopWrite = 0; hopFilled = 0; blockCount = 0;
        stCount = 0; stSince = 0;
        std::fill (hopRing.begin(), hopRing.end(), 0.0);
    }

    void setChannelWeight (int c, double weight) noexcept { if (c >= 0 && c < kMaxChannels) w[c] = weight; }

    void process (const float* const* channels, int numChannels, int n) noexcept
    {
        const int nc = numChannels < ch ? numChannels : ch;
        for (int i = 0; i < n; ++i)
        {
            for (int c = 0; c < nc; ++c) { const double y = kw.process (c, (double) channels[c][i]); hopSumSq[c] += y * y; }
            if (++hopCount >= hopSamples) finishHop (nc);
        }
    }

    double momentaryLufs()  const noexcept { return lufsOf (meanLastHops (4)); }    // 400 ms
    double shortTermLufs()  const noexcept { return lufsOf (meanLastHops (30)); }   // 3 s
    double integratedLufs() const noexcept { return integrated(); }                 // gated
    double loudnessRangeLu() const noexcept { return lra(); }                       // EBU Tech 3342 (P95−P10)

private:
    static constexpr int kMaxChannels = core::kMaxChannels;
    static constexpr int kHopRing = 30;                                             // 30 × 100 ms = 3 s

    void finishHop (int nc) noexcept
    {
        double hopMS = 0.0;
        for (int c = 0; c < nc; ++c) hopMS += w[c] * (hopSumSq[c] / (double) hopSamples);
        hopRing[(std::size_t) hopWrite] = hopMS;
        hopWrite = (hopWrite + 1) % kHopRing;
        if (hopFilled < kHopRing) ++hopFilled;
        if (hopFilled >= 4 && blockCount < (int) blockE.size())                     // a 400 ms block every 100 ms
            blockE[(std::size_t) blockCount++] = meanLastHops (4);
        if (hopFilled >= kHopRing)                                                  // a 3 s short-term sample every 1 s (LRA)
        {
            if (stSince == 0 && stCount < (int) stE.size()) stE[(std::size_t) stCount++] = meanLastHops (kHopRing);
            if (++stSince >= 10) stSince = 0;                                        // first at 3 s, then every 10 hops (libebur128)
        }
        for (int c = 0; c < nc; ++c) hopSumSq[c] = 0.0;
        hopCount = 0;
    }

    double meanLastHops (int k) const noexcept
    {
        const int kk = k < hopFilled ? k : hopFilled;
        if (kk <= 0) return 0.0;
        double s = 0.0;
        for (int j = 0; j < kk; ++j) { const int idx = (hopWrite - 1 - j + kHopRing) % kHopRing; s += hopRing[(std::size_t) idx]; }
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

    double fs = 48000.0; int ch = 2, hopSamples = 4800;
    KWeightingFilter kw;
    double w[kMaxChannels] {};
    double hopSumSq[kMaxChannels] {};
    int hopCount = 0;
    std::vector<double> hopRing;
    int hopWrite = 0, hopFilled = 0;
    std::vector<double> blockE;
    int blockCount = 0;
    std::vector<double> stE;                                                        // 3 s short-term energies @1 s (LRA)
    int stCount = 0, stSince = 0;
};

} // namespace felitronics::analysis

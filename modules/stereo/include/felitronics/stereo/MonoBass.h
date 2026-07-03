// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/stereo/MidSide.h>
#include <felitronics/eq/Svf.h>

#include <algorithm>
#include <cmath>

namespace felitronics::stereo
{

//==============================================================================
// felitronics::stereo::MonoBass — "elliptical EQ" / bass mono-maker: collapse the low end to mono below a
// cutoff (vinyl/LP lacquer cutting — out-of-phase lows cause vertical stylus excursion — and general
// mastering tightening). Zero added latency, minimum-phase.
//
// Design — verified by math + the self-tests. An allpass-subtraction variant was checked and REJECTED:
// AP4−LP4 has magnitude 1.5 at the cutoff, i.e. it boosts, not rejects.
// Operate ONLY on the Side of a Mid/Side split, so the Mid (the mono content you keep) passes UNFILTERED —
// no transient smear on the kept bass; only the low-frequency Side (stereo) energy is removed. A full L/R
// crossover would needlessly filter the mono bass too.
//
//   M = ½(L+R),  S = ½(L-R)
//   S' = lowWidth·LP4(S) + HP4(S)
//   L' = M + S',  R' = M − S'
//
// LP4/HP4 are 4th-order Linkwitz-Riley (two cascaded Butterworth Q=1/√2 sections each). The LR4 identity
// LP4 + HP4 = allpass (proven flat: |LP4+HP4| = (1+ω⁴)/(1+ω⁴) = 1) makes the `lowWidth` taper exact:
// lowWidth=0 → S'=HP4(S) (mono below fc); lowWidth=1 → S'=allpass(S) (side magnitude preserved → we bypass
// entirely there for bit-transparency). Above fc, HP4≈1 so the full stereo image is intact.
//
// RT-safe: no alloc/lock/throw in process(); four SVF instances hold the only state (the Side is one channel).
class MonoBass
{
public:
    void prepare (double sampleRate, int /*maxBlock*/ = 0, int /*maxChannels*/ = 2) noexcept
    {
        fs_ = sampleRate > 0.0 ? sampleRate : 48000.0;
        lp1_.prepare (fs_, 1); lp2_.prepare (fs_, 1);
        hp1_.prepare (fs_, 1); hp2_.prepare (fs_, 1);
        applyFilters();
        reset();
    }

    void reset() noexcept { lp1_.reset(); lp2_.reset(); hp1_.reset(); hp2_.reset(); }

    void setEnabled (bool e) noexcept { enabled_ = e; }
    // Setters reject non-finite values (std::clamp passes NaN through — it would poison the SVF coefficients).
    void setFrequency (float hz) noexcept { if (! std::isfinite (hz)) return; freq_ = std::clamp (hz, 20.0f, (float) (0.45 * fs_)); applyFilters(); }
    void setLowWidth (float w) noexcept { if (! std::isfinite (w)) return; lowWidth_ = std::clamp (w, 0.0f, 1.0f); }   // 0 = mono below fc, 1 = full stereo

    bool  isEnabled() const noexcept { return enabled_; }
    float frequency() const noexcept { return freq_; }
    float lowWidth()  const noexcept { return lowWidth_; }
    static constexpr int latencySamples() noexcept { return 0; }

    // Stereo, in place: io[0]=L, io[1]=R. Bypass (signal untouched) when disabled, < 2 channels, or fully wide.
    void process (float* const* io, int numChannels, int n) noexcept
    {
        if (numChannels < 2 || ! enabled_ || lowWidth_ >= 1.0f - 1.0e-6f) return;
        float* L = io[0];
        float* R = io[1];
        for (int i = 0; i < n; ++i)
        {
            float m, s; MidSide::encode (L[i], R[i], m, s);
            const float lp = lp2_.processSample (0, lp1_.processSample (0, s));   // LP4(S)
            const float hp = hp2_.processSample (0, hp1_.processSample (0, s));   // HP4(S)
            const float sOut = lowWidth_ * lp + hp;                               // lowWidth·LP4 + HP4
            MidSide::decode (m, sOut, L[i], R[i]);
        }
        lp1_.flushDenormals(); lp2_.flushDenormals(); hp1_.flushDenormals(); hp2_.flushDenormals();
    }

private:
    void applyFilters() noexcept   // coefficients only — no state reset, so a freq move doesn't click
    {
        constexpr double Q = 0.7071067811865476;   // 1/√2 Butterworth → cascade = 4th-order Linkwitz-Riley
        lp1_.setParams (eq::FilterType::LowPass,  freq_, Q, 0.0);
        lp2_.setParams (eq::FilterType::LowPass,  freq_, Q, 0.0);
        hp1_.setParams (eq::FilterType::HighPass, freq_, Q, 0.0);
        hp2_.setParams (eq::FilterType::HighPass, freq_, Q, 0.0);
    }

    double fs_ = 48000.0;
    float  freq_ = 150.0f, lowWidth_ = 0.0f;
    bool   enabled_ = true;
    eq::Svf lp1_, lp2_, hp1_, hp2_;   // Side-channel LR4 low-pass + high-pass cascades
};

} // namespace felitronics::stereo

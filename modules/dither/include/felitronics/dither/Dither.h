// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/core/Config.h>
#include <felitronics/core/Math.h>

#include <algorithm>
#include <cmath>
#include <cstdint>

namespace felitronics::dither
{

//==============================================================================
// felitronics::dither::Dither — the LAST stage before bit-depth reduction on export (float → 16/20/24-bit):
// TPDF dither (eliminates quantization-noise modulation/distortion) + optional error-feedback noise shaping
// (pushes the residual quantization noise out of the ear's sensitive band). The structure was
// cross-checked against the literature and the sign/curve confirmed by the self-tests.
//
//   lsb   = 1 / 2^(bits-1)                                   (float audio in [-1, 1) → 2^bits codes)
//   u[n]  = x[n] − Σ_k h_k · e[n−k]·lsb                       (error feedback; SUBTRACT → NTF(z) = 1 − H(z))
//   v[n]  = u[n] + d[n],  d = (½LSB-uniform + ½LSB-uniform)   (TPDF, peak ±1 LSB → noise-modulation-free)
//   y[n]  = quantize(v[n]) = round(v·2^(bits-1)) · lsb        (mid-tread, clamped to signed-PCM range)
//   e[n]  = (y[n] − u[n]) / lsb                               (TOTAL added noise = dither + quant residue)
//
//   Y(z) = X(z) + (1 − H(z))·E(z)   →   ALL the added noise E (the dither INCLUDED) is shaped by NTF = 1 − H.
//   Feeding back y−u (not y−v) is what lets the curve lower the dither floor too — shape only the quantizer
//   residue and the flat dither sets an unshakable in-band floor (the psychoacoustic curves then do nothing).
//
// Shaping (NTF = 1 − H, H = Σ h_k z^−k, h[0] multiplies e[n−1]):
//   None           — flat TPDF only.
//   Weighted       — H = 1.5 z^−1 − 0.5 z^−2 → NTF = (1−z^−1)(1−0.5 z^−1): a gentle, sample-rate-agnostic
//                    high-pass with a DC zero. THE MASTERING DEFAULT — benign, no overload risk.
//   Psychoacoustic — the 9-tap "F-weighted" curve (Wannamaker/Lipshitz, designed for 44.1 kHz): a deeper
//                    notch around the 3–4 kHz hearing peak. Opt-in; best at 44.1 kHz.
//
// TPDF source: one PCG32 stream per channel (SplitMix64-seeded, per-channel decorrelated → independent L/R
// dither), two draws per sample. Auto-blanking: after `autoBlankSamples` of EXACT digital silence a channel
// outputs true zero (preserves digital black on exported tails) and clears its shaper; instant resume on
// signal. Non-finite input → 0; output codes hard-clamped. bits ≥ 32 → bypass (32-bit float export).
//
// RT-safe: no alloc/lock/throw in process(); fixed per-channel state (no heap); zero latency.
enum class NoiseShaping { None, Weighted, Psychoacoustic };

struct DitherParams
{
    int bits = 24;                                   // 16 / 20 / 24 target; ≥ 32 → bypass (float export)
    NoiseShaping shaping = NoiseShaping::Weighted;
    std::uint64_t seed = 0x853c49e6748fea9bULL;      // base seed; per-channel decorrelated internally
    bool autoBlank = true;                           // mute dither on sustained exact digital silence
    int  autoBlankSamples = 4096;                    // consecutive zero samples before a channel blanks
};

class Dither
{
public:
    void prepare (double sampleRate, int /*maxBlock*/, int maxChannels) noexcept
    {
        fs_ = sampleRate > 0.0 ? sampleRate : 48000.0;
        channels_ = std::clamp (maxChannels, 1, core::kMaxChannels);
        apply (params_);
        reset();
    }

    void reset() noexcept
    {
        for (int c = 0; c < core::kMaxChannels; ++c)
        {
            seedChannel (c);
            ch_[c].clearShaper();
            ch_[c].blank = 0;
        }
    }

    void setParams (const DitherParams& p) noexcept
    {
        const bool seedChanged    = (p.seed != params_.seed);
        const bool shapingChanged = (p.shaping != params_.shaping) || (p.bits != params_.bits);
        params_ = p; apply (p);
        if (seedChanged)         reset();                 // new stream → reseed + clear
        else if (shapingChanged) for (int c = 0; c < core::kMaxChannels; ++c) ch_[c].clearShaper();
    }

    static constexpr int latencySamples() noexcept { return 0; }

    void process (float* const* io, int numChannels, int n) noexcept
    {
        if (bits_ >= 32 || numChannels <= 0) return;      // float export → bypass
        const int nc = std::min (numChannels, core::kMaxChannels);   // state exists for every channel → never silently skip
        for (int c = 0; c < nc; ++c)
        {
            Channel& st = ch_[c];
            float* x = io[c];
            for (int i = 0; i < n; ++i)
            {
                float in = x[i];
                if (! std::isfinite (in)) in = 0.0f;       // a stray NaN/inf must not poison the export

                if (autoBlank_ && core::exactlyEqual (in, 0.0f))   // intentional exact == (digital black only)
                {
                    if (st.blank < blankSamples_) ++st.blank;
                    if (st.blank >= blankSamples_) { st.clearShaper(); x[i] = 0.0f; continue; }   // digital black
                }
                else st.blank = 0;

                double fb = 0.0;                            // Σ h_k e[n−k]  (LSB units)
                for (int k = 0; k < order_; ++k) fb += coeffs_[k] * st.e[k];
                const double u = (double) in - fb * lsb_;   // SUBTRACT → NTF = 1 − H
                const double d = tpdf (st) * lsb_;          // ±1 LSB triangular
                const double v = u + d;

                // CLAMP IN THE DOUBLE DOMAIN, BEFORE THE CONVERSION. The clamp used to come after it,
                // which is too late: a floating-to-integer conversion whose value is outside the target
                // range is undefined, and `v` is only bounded by what the caller sent. Non-finite input
                // is already zeroed above, but a large FINITE one is not — nothing here clamps
                // magnitude — so a stage upstream running without a limiter (a bypassed one, plus a
                // makeup gain) reaches it: at 24 bits `scale_` is 8388608, so any |v| past ~1.1e12
                // overflows int64 before the old clamp could act. Bit-identical inside the code range,
                // where the bounds are exactly representable in double and the branch is not taken.
                const double scaled = std::floor (v * scale_ + 0.5);   // mid-tread round-half-up
                const std::int64_t code = scaled <= (double) minCode_ ? minCode_
                                        : scaled >= (double) maxCode_ ? maxCode_
                                                                      : (std::int64_t) scaled;
                const double y = (double) code * lsb_;

                double eLsb = (y - u) * scale_;             // TOTAL added noise d+q (LSB units) → shapes the dither too
                eLsb = std::clamp (eLsb, -8.0, 8.0);        // defensive: bound the 9th-order feedback
                for (int k = order_ - 1; k > 0; --k) st.e[k] = st.e[k - 1];
                if (order_ > 0) st.e[0] = eLsb;

                x[i] = (float) y;
            }
        }
    }

private:
    static constexpr int kMaxOrder = 9;

    struct Channel
    {
        std::uint64_t state = 0, inc = 0;                  // PCG32
        double e[kMaxOrder] = {};                          // error history, LSB units (e[0] = e[n−1])
        int blank = 0;
        void clearShaper() noexcept { for (double& v : e) v = 0.0; }
    };

    // --- PCG32 (O'Neill) — one stream per channel, two draws per sample for the TPDF pair ---
    static inline std::uint32_t pcg (Channel& s) noexcept
    {
        const std::uint64_t old = s.state;
        s.state = old * 6364136223846793005ULL + s.inc;
        const std::uint32_t xs  = (std::uint32_t) (((old >> 18) ^ old) >> 27);
        const std::uint32_t rot = (std::uint32_t) (old >> 59);
        return (xs >> rot) | (xs << ((0u - rot) & 31u));
    }
    static inline double uni (Channel& s) noexcept { return (double) (pcg (s) >> 8) * 0x1.0p-24 - 0.5; }   // [-0.5, 0.5)
    static inline double tpdf (Channel& s) noexcept { return uni (s) + uni (s); }                          // triangular [-1, 1)

    static inline std::uint64_t splitmix64 (std::uint64_t& x) noexcept
    {
        std::uint64_t z = (x += 0x9E3779B97F4A7C15ULL);
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
        return z ^ (z >> 31);
    }

    void seedChannel (int c) noexcept
    {
        constexpr std::uint64_t kGolden = 0x9E3779B97F4A7C15ULL;   // uint64_t, not ULL literal: keeps the ^/* in one type (GCC -Wsign-conversion)
        std::uint64_t s = params_.seed ^ (kGolden * (std::uint64_t) (c + 1));
        const std::uint64_t a = splitmix64 (s), b = splitmix64 (s);
        ch_[c].state = 0u; ch_[c].inc = (b << 1) | 1u;     // inc must be odd
        pcg (ch_[c]); ch_[c].state += a; pcg (ch_[c]);     // PCG seeding ritual
    }

    void apply (const DitherParams& p) noexcept
    {
        bits_ = p.bits;
        const int b = std::clamp (p.bits, 2, 31);
        scale_ = (double) ((std::int64_t) 1 << (b - 1));
        lsb_   = 1.0 / scale_;
        minCode_ = -((std::int64_t) 1 << (b - 1));
        maxCode_ =  ((std::int64_t) 1 << (b - 1)) - 1;
        autoBlank_ = p.autoBlank;
        blankSamples_ = std::max (1, p.autoBlankSamples);

        switch (p.shaping)
        {
            case NoiseShaping::None:    order_ = 0; break;
            case NoiseShaping::Weighted:
                order_ = 2; coeffs_[0] = 1.5; coeffs_[1] = -0.5; break;
            case NoiseShaping::Psychoacoustic:
            {
                order_ = 9;
                static constexpr double f9[9] = { 2.412, -3.370, 3.937, -4.174, 3.353, -2.205, 1.281, -0.569, 0.0847 };
                for (int k = 0; k < 9; ++k) coeffs_[k] = f9[k];
                break;
            }
        }
    }

    double fs_ = 48000.0;
    int channels_ = 2;
    DitherParams params_;

    int bits_ = 24, order_ = 2;
    double scale_ = 8388608.0, lsb_ = 1.0 / 8388608.0;
    std::int64_t minCode_ = -8388608, maxCode_ = 8388607;
    double coeffs_[kMaxOrder] = { 1.5, -0.5 };
    bool autoBlank_ = true;
    int  blankSamples_ = 4096;

    Channel ch_[core::kMaxChannels];
};

} // namespace felitronics::dither

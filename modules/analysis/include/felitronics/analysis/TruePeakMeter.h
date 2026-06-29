// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/core/Config.h>
#include <felitronics/core/Math.h>

#include <algorithm>
#include <cmath>
#include <vector>

namespace felitronics::analysis
{

namespace detail
{
    inline double tpBesselI0 (double x) noexcept
    {
        double sum = 1.0, term = 1.0;
        const double y = x * x * 0.25;
        for (int k = 1; k < 64; ++k) { term *= y / ((double) k * (double) k); sum += term; if (term < 1e-17 * sum) break; }
        return sum;
    }
}

//==============================================================================
// felitronics::analysis::TruePeakMeter — inter-sample (true) peak per ITU-R BS.1770-4 / EBU R128 (dBTP):
// oversample, then take the max of the reconstructed waveform — the inter-sample peaks a raw sample meter
// misses. A READ-ONLY sink (it never touches the audio), so it adds ZERO latency; the interpolation FIR's
// group delay only postpones *discovery* of a peak, which is irrelevant to a running max (DSP council:
// codex + deepseek).
//
// The interpolator is a Kaiser-windowed-sinc polyphase FIR (factor × tapsPerPhase), each phase normalised to
// unity DC so a constant reads exactly and the sample grid is reproduced — the meter therefore can never
// under-read the sample peak. BS.1770-4 specifies a particular table; the Recommendation also allows an
// equivalent filter, so we DESIGN one meeting the envelope (flat pass-band, ≥ ~60 dB image rejection) and
// VERIFY it against the canonical fs/4 inter-sample test rather than trusting a copied coefficient list.
//
// Oversampling factor by source rate (targets the spec's ~176.4/192 kHz analysis rate): < 88.2 kHz → 4×,
// < 176.4 kHz → 2×, else 1× (the grid already resolves the peak). It is a PURE MAX — no integration, no
// gating, no K-weighting (that is loudness/LUFS, a separate meter). truePeakDb()/samplePeakDb() are
// free-running maxima since reset(); an optional hold+decay ballistic (displayTruePeakDb()) is for a moving
// display only and never feeds the authoritative max. RT-safe: prepare() allocates; process() does not.
struct TruePeakMeterParams
{
    int    oversample    = 0;        // 0 = auto by sample rate; else force 1 / 2 / 4
    double holdMs        = 1000.0;   // display ballistic only (does NOT affect truePeakDb())
    double decayDbPerSec = 20.0;     // display ballistic fall rate
};

class TruePeakMeter
{
public:
    static constexpr int kTapsPerPhase = 12;     // 48-tap prototype at 4× — matches the spec filter length

    void prepare (double sampleRate, int /*maxBlock*/, int maxChannels) noexcept
    {
        fs_ = sampleRate > 0.0 ? sampleRate : 48000.0;
        channels_ = std::clamp (maxChannels, 1, core::kMaxChannels);
        chooseFactor();
        designFilter();
        hist_.assign ((std::size_t) channels_ * (std::size_t) kTapsPerPhase, 0.0f);
        pos_.assign  ((std::size_t) channels_, 0);
        applyBallistics();
        reset();
    }

    // NB: a factor change re-designs the FIR (allocates) → call from the prepare/stopped context, not the
    // audio thread. Ballistic/seed-only changes are RT-safe.
    void setParams (const TruePeakMeterParams& p) noexcept
    {
        const bool refactor = (p.oversample != params_.oversample);
        params_ = p;
        if (refactor) { chooseFactor(); designFilter(); reset(); }   // reset → no mixed-filter maxima
        applyBallistics();
    }

    void reset() noexcept
    {
        std::fill (hist_.begin(), hist_.end(), 0.0f);
        std::fill (pos_.begin(),  pos_.end(),  0);
        truePeakLin_ = samplePeakLin_ = blockTpLin_ = 0.0f;
        holdLin_ = 0.0f; holdCount_ = 0;
    }

    static constexpr int latencySamples() noexcept { return 0; }
    int    oversampleFactor() const noexcept { return L_; }
    double truePeakLinear()   const noexcept { return truePeakLin_; }
    double truePeakDb()       const noexcept { return toDb (truePeakLin_); }     // authoritative max since reset
    double samplePeakDb()     const noexcept { return toDb (samplePeakLin_); }
    double truePeakDbBlock()  const noexcept { return toDb (blockTpLin_); }      // max of the last process() block
    double displayTruePeakDb()const noexcept { return toDb (holdLin_); }         // hold/decay ballistic (display)

    // READ-ONLY: io is sampled, never modified.
    void process (const float* const* io, int numChannels, int n) noexcept
    {
        const int nc = std::min (numChannels, channels_);
        if (nc <= 0) return;
        float blockMax = 0.0f;
        for (int c = 0; c < nc; ++c)
        {
            int pos = pos_[(std::size_t) c];
            float* h = &hist_[(std::size_t) c * (std::size_t) kTapsPerPhase];
            for (int m = 0; m < n; ++m)
            {
                float x = io[c][m];
                if (! std::isfinite (x) || std::fabs (x) < 1.0e-30f) x = 0.0f;   // guard NaN/inf, kill denormals
                const float ax = std::fabs (x);
                if (ax > samplePeakLin_) samplePeakLin_ = ax;

                h[pos] = x; if (++pos >= kTapsPerPhase) pos = 0;     // pos-1 = newest x[m]
                float tp = ax;                                        // the grid sample itself (→ TP ≥ sample peak)
                for (int p = 0; p < L_ && L_ > 1; ++p)
                {
                    float acc = 0.0f;
                    for (int k = 0; k < kTapsPerPhase; ++k)
                    {
                        int hi = pos - 1 - k; if (hi < 0) hi += kTapsPerPhase;
                        acc += proto_[(std::size_t) (k * L_ + p)] * h[hi];
                    }
                    const float v = std::fabs ((float) L_ * acc);     // ×L restores the polyphase gain
                    if (v > tp) tp = v;
                }
                if (tp > blockMax) blockMax = tp;

                if (holdSamples_ > 0)                                 // display ballistic (per sample, base rate)
                {
                    if (tp >= holdLin_) { holdLin_ = tp; holdCount_ = holdSamples_; }
                    else if (holdCount_ > 0) --holdCount_;
                    else { holdLin_ *= decayMul_; if (holdLin_ < 1.0e-20f) holdLin_ = 0.0f; }   // no denormal tail
                }
            }
            pos_[(std::size_t) c] = pos;
        }
        blockTpLin_ = blockMax;
        if (blockMax > truePeakLin_) truePeakLin_ = blockMax;
    }

private:
    void chooseFactor() noexcept
    {
        if (params_.oversample == 1 || params_.oversample == 2 || params_.oversample == 4) { L_ = params_.oversample; return; }
        L_ = (fs_ < 88200.0) ? 4 : (fs_ < 176400.0 ? 2 : 1);          // reach the spec's ~176.4/192 kHz analysis rate
    }

    void designFilter() noexcept
    {
        const int N = L_ * kTapsPerPhase;
        proto_.assign ((std::size_t) std::max (1, N), 0.0f);
        if (L_ <= 1) return;                                          // no interpolation needed at ≥ 176.4 kHz
        const double fc   = 0.5 / (double) L_;                        // cutoff at the BASE Nyquist (full true-peak bandwidth)
        const double cen  = (double) (N - 1) * 0.5;
        const double beta = 8.0;                                      // Kaiser: ~ −70 dB images, < 0.1 dB pass-band ripple
        const double i0b  = detail::tpBesselI0 (beta);
        double sum = 0.0;
        for (int i = 0; i < N; ++i)
        {
            const double xx   = (double) i - cen;
            const double sinc = (std::fabs (xx) < 1e-9) ? (2.0 * fc) : std::sin (2.0 * core::kPi * fc * xx) / (core::kPi * xx);
            const double r    = (double) (2 * i - (N - 1)) / (double) (N - 1);
            const double win  = detail::tpBesselI0 (beta * std::sqrt (std::max (0.0, 1.0 - r * r))) / i0b;
            const double v    = sinc * win;
            proto_[(std::size_t) i] = (float) v; sum += v;
        }
        const float inv = (float) (1.0 / sum);                        // Σ → 1 → each phase ≈ 1/L → unity DC after ×L
        for (auto& v : proto_) v *= inv;
    }

    void applyBallistics() noexcept
    {
        holdSamples_ = (int) std::max (0.0, params_.holdMs * 0.001 * fs_);
        decayMul_    = params_.decayDbPerSec > 0.0 ? (float) std::pow (10.0, -params_.decayDbPerSec / (20.0 * fs_)) : 1.0f;
    }

    static double toDb (float lin) noexcept { return lin > 1.0e-10f ? core::gainToDb ((double) lin) : -200.0; }

    double fs_ = 48000.0;
    int channels_ = 2, L_ = 4;
    TruePeakMeterParams params_;

    std::vector<float> proto_;                  // N = L*tapsPerPhase taps, Σ = 1
    std::vector<float> hist_;                   // per-channel ring (tapsPerPhase)
    std::vector<int>   pos_;

    float truePeakLin_ = 0.0f, samplePeakLin_ = 0.0f, blockTpLin_ = 0.0f;
    float holdLin_ = 0.0f; int holdCount_ = 0, holdSamples_ = 0;
    float decayMul_ = 1.0f;
};

} // namespace felitronics::analysis

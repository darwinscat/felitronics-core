// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/core/Config.h>
#include <felitronics/core/Math.h>

#include <algorithm>
#include <cmath>
#include <vector>

namespace felitronics::oversampling
{

namespace detail
{
    inline double besselI0 (double x) noexcept
    {
        double sum = 1.0, term = 1.0;
        const double y = x * x * 0.25;
        for (int k = 1; k < 64; ++k) { term *= y / ((double) k * (double) k); sum += term; if (term < 1e-17 * sum) break; }
        return sum;
    }
}

//==============================================================================
// felitronics::oversampling::PolyphaseOversampler — integer up/down sampling by a polyphase
// windowed-sinc (Kaiser) FIR. Linear phase → a constant group delay (reported). For inter-sample-peak
// detection (the true-peak limiter) and, with up+down around a process, alias-free nonlinear processing.
// RT-safe: prepare() allocates the FIR + per-channel histories; up/downsample() do no alloc/lock/throw.
// Scalar reference (correctness-first); a SIMD backend can replace it later behind the same API.
class PolyphaseOversampler
{
public:
    // factor 2/4/8; tapsPerPhase = FIR taps per polyphase branch (filter length = factor*tapsPerPhase).
    bool prepare (int factor, int maxChannels, int tapsPerPhase = 32)
    {
        if (factor < 2 || tapsPerPhase < 4) return false;
        L = factor; tpp = tapsPerPhase; N = L * tpp;
        channels_ = maxChannels < 1 ? 1 : (maxChannels > core::kMaxChannels ? core::kMaxChannels : maxChannels);
        designFilter();
        upHist.assign   ((std::size_t) channels_ * tpp, 0.0f);
        downHist.assign ((std::size_t) channels_ * N,   0.0f);
        upPos.assign    ((std::size_t) channels_, 0);
        downPos.assign  ((std::size_t) channels_, 0);
        return true;
    }

    void reset() noexcept
    {
        std::fill (upHist.begin(),   upHist.end(),   0.0f);
        std::fill (downHist.begin(), downHist.end(), 0.0f);
        std::fill (upPos.begin(),    upPos.end(),    0);
        std::fill (downPos.begin(),  downPos.end(),  0);
    }

    int factor() const noexcept { return L; }
    // Group delay of ONE filter pass, in OVERSAMPLED samples ((N-1)/2 for a linear-phase FIR).
    double filterLatencyOversampled() const noexcept { return (double) (N - 1) * 0.5; }
    // Round-trip (up THEN down) latency in BASEBAND samples. EXACT integer: the two (N-1)/2 OS group
    // delays plus the decimation phase (L-1) combine to (N - L)/L = tpp - 1 baseband samples.
    // Unprepared (tpp == 0) reports 0, not -1 — hosts query latency before prepare().
    int latencySamples() const noexcept { return tpp > 0 ? tpp - 1 : 0; }

    // n baseband samples per channel → n*L oversampled samples. out[ch] holds n*L.
    void upsample (const float* const* in, int channels, int n, float* const* out) noexcept
    {
        const int nc = channels < channels_ ? channels : channels_;
        for (int c = 0; c < nc; ++c)
        {
            int   pos  = upPos[(std::size_t) c];
            float* h   = &upHist[(std::size_t) c * tpp];
            for (int m = 0; m < n; ++m)
            {
                h[pos] = in[c][m];
                if (++pos >= tpp) pos = 0;                       // pos-1 = newest x[m]
                for (int p = 0; p < L; ++p)
                {
                    float acc = 0.0f;
                    for (int k = 0; k < tpp; ++k)
                    {
                        int hi = pos - 1 - k; if (hi < 0) hi += tpp;
                        acc += proto[(std::size_t) (k * L + p)] * h[hi];
                    }
                    out[c][m * L + p] = (float) L * acc;          // *L restores the zero-stuff gain
                }
            }
            upPos[(std::size_t) c] = pos;
        }
    }

    // n*L oversampled samples per channel → n baseband samples (lowpass + decimate). out[ch] holds n.
    void downsample (const float* const* in, int channels, int n, float* const* out) noexcept
    {
        const int nc = channels < channels_ ? channels : channels_;
        for (int c = 0; c < nc; ++c)
        {
            int   pos = downPos[(std::size_t) c];
            float* h  = &downHist[(std::size_t) c * N];
            for (int m = 0; m < n; ++m)
            {
                for (int p = 0; p < L; ++p) { h[pos] = in[c][m * L + p]; if (++pos >= N) pos = 0; }
                float acc = 0.0f;
                for (int j = 0; j < N; ++j) { int hi = pos - 1 - j; if (hi < 0) hi += N; acc += proto[(std::size_t) j] * h[hi]; }
                out[c][m] = acc;                                  // proto sums to 1 → unity passband
            }
            downPos[(std::size_t) c] = pos;
        }
    }

private:
    void designFilter()
    {
        proto.assign ((std::size_t) N, 0.0f);
        const double fc   = 0.5 / (double) L * 0.90;              // cutoff (cycles/OS-sample), guard below baseband Nyquist
        const double cen  = (double) (N - 1) * 0.5;
        const double beta = 9.0;                                  // Kaiser ~ -90 dB stopband
        const double i0b  = detail::besselI0 (beta);
        double sum = 0.0;
        for (int i = 0; i < N; ++i)
        {
            const double x    = (double) i - cen;
            const double sinc = (std::fabs (x) < 1e-9) ? (2.0 * fc)
                                                       : std::sin (2.0 * core::kPi * fc * x) / (core::kPi * x);
            const double r    = (double) (2 * i - (N - 1)) / (double) (N - 1);   // ∈ [-1,1]
            const double win  = detail::besselI0 (beta * std::sqrt (std::max (0.0, 1.0 - r * r))) / i0b;
            const double v    = sinc * win;
            proto[(std::size_t) i] = (float) v;
            sum += v;
        }
        const float inv = (float) (1.0 / sum);                   // normalize Σ → 1 (unity DC)
        for (auto& v : proto) v *= inv;
    }

    int L = 0, tpp = 0, N = 0, channels_ = 0;
    std::vector<float> proto;                  // N taps, Σ = 1
    std::vector<float> upHist, downHist;       // per-channel ring histories
    std::vector<int>   upPos, downPos;
};

} // namespace felitronics::oversampling

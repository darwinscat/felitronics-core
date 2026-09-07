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
//
// WHY THE DEFAULT IS 64, AND WHAT SETS IT. The cutoff is FIXED at 0.90 x baseband Nyquist
// (designFilter() below), so the transition band has to fit between 0.45 fs and the fold at 0.50 fs,
// and tapsPerPhase is the only thing that decides whether it does. The design DECLARES its own target
// one line down — `beta = 9.0` is a ~90 dB Kaiser stopband — and everything above 0.50 fs folds
// straight back into the audio band, so that number is the filter's job rather than a nicety.
// Worst |H| over the whole fold region [0.5, factor/2] fs — the WORST of factors 2, 4 and 8, from a
// 400k-point scan (they agree to 0.5 dB, and the worst is at 2x above 60 taps):
//
//        tpp      32       40       48       56       60       64       72       96      128
//        dB    -26.9    -37.0    -51.1    -76.2    -90.7    -90.5    -91.9    -94.3    -97.9
//
// So 60 is where the filter first DELIVERS the design it declares, and above it the stopband is on the
// window's floor — further taps buy pass-band width, not rejection. THE DEFAULT USED TO BE 32, i.e. the
// prototype delivered 27 dB where it promised 90.
//
// End to end, and READ ITS SCOPE: on a Saturator (0.17 fs tone, tanh at +24 dB, 4x) the 3rd harmonic
// lands at 0.51 fs — inside the transition band — and folded back to 0.49 fs at -44.7 dBc, where 64
// taps put it at -94.0. That is the component the TAPS own. It is NOT a claim about total aliasing: a
// tanh has infinitely many harmonics and the ones above the OS Nyquist fold inside the oversampled
// domain, which no decimation filter reaches, so the total non-harmonic energy of a hard-driven
// waveshaper barely moves (-31.97 -> -31.27 dBc at +24 dB of drive — marginally WORSE, since a flatter
// pass band also delivers what had already folded). Total aliasing is the FACTOR's axis.
// Pass-band flatness is a COROLLARY of the same width, not a second requirement: two oversampled
// stages in series (the real assembly — a clipper in front of a limiter) then droop under 0.04 dB to
// 0.41 fs, against -3.25 dB at 32. The pass-band edge follows 0.45 fs - c/tpp with c = 2.385 (0.1 dB),
// 2.505 (0.05 dB), 2.723 (0.01 dB), so a WIDER flat band is bought hyperbolically: 0.42 fs costs 80
// taps, 0.43 costs 128, 0.44 costs 239. Do not buy it that way — 20 kHz at 44.1 kHz is 0.4535 fs, ABOVE
// the fixed cutoff, so no tapsPerPhase reaches it and more taps make it WORSE (-27.4 dB at 32 against
// -31.1 at 64). Widening the audible band is the CUTOFF's axis, not this one.
// Cost, both linear in tapsPerPhase: the FIR (measured 2.09x for the clipper+limiter pair, 2.13 %RT ->
// 4.44 %RT at 48 kHz stereo 4x) and the reported latency (tapsPerPhase - 1 baseband samples).
class PolyphaseOversampler
{
public:
    // The shipped topology decision, named so a consumer can state "the core's default" rather than
    // re-spell the number — and so a test can assert the DEFAULT moved, not just that some number did.
    static constexpr int kDefaultTapsPerPhase = 64;
    // BOTH factors of N = factor*tapsPerPhase are bounded, HERE rather than at each caller, because the
    // product is what overflows and either argument alone can do it. tapsPerPhase = INT_MAX overflowed the
    // multiply below and then threw a length_error out of assign() — a terminate under the wasm tier's
    // -fno-exceptions. TruePeakLimiter::prepare refused both at its own gate (its contract, kept), but
    // Saturator passes an unbounded oversampleFactor straight through, and `prepare(INT_MAX, 1)` on the
    // default taps is UBSan-confirmed signed overflow followed by that same length_error. The ceiling is
    // generous rather than tight: the largest factor anything in the tree asks for is 32 (PowerAmpStage's
    // alias-free reference), and 64 x 1024 taps is still a trivially small allocation.
    static constexpr int kMaxTapsPerPhase = 1024;
    static constexpr int kMaxFactor       = 64;

    // factor 2/4/8; tapsPerPhase = FIR taps per polyphase branch (filter length = factor*tapsPerPhase).
    // The default is the topology decision above; pass it explicitly to pin one against this default.
    bool prepare (int factor, int maxChannels, int tapsPerPhase = kDefaultTapsPerPhase)
    {
        if (factor < 2 || factor > kMaxFactor) return false;
        if (tapsPerPhase < 4 || tapsPerPhase > kMaxTapsPerPhase) return false;
        L = factor; tpp = tapsPerPhase; N = L * tpp;
        channels_ = maxChannels < 1 ? 1 : (maxChannels > core::kMaxChannels ? core::kMaxChannels : maxChannels);
        designFilter();
        upHist.assign   ((std::size_t) channels_ * (std::size_t) tpp, 0.0f);
        downHist.assign ((std::size_t) channels_ * (std::size_t) N,   0.0f);
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

    // Clear ONE channel's histories, leaving every other channel BIT-EXACT — for an owner whose
    // channel stopped being fed and will be fed again. The ring POSITIONS go with the samples: they
    // are state like any other, and a channel whose ring is zeroed while its cursor stays put sits at
    // a phase reset() can never produce, so "this column now equals a fresh one" would be false about
    // a column that looks clean. Unprepared (channels_ == 0) is a no-op, not a bad index.
    void resetChannel (int c) noexcept
    {
        if (c < 0 || c >= channels_) return;
        std::fill_n (upHist.begin()   + (std::ptrdiff_t) c * (std::ptrdiff_t) tpp, tpp, 0.0f);
        std::fill_n (downHist.begin() + (std::ptrdiff_t) c * (std::ptrdiff_t) N,   N,   0.0f);
        upPos[(std::size_t) c]   = 0;
        downPos[(std::size_t) c] = 0;
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
            float* h   = &upHist[(std::size_t) c * (std::size_t) tpp];
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
            float* h  = &downHist[(std::size_t) c * (std::size_t) N];
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

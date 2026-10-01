// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/codecgrid/Windows.h>   // window::detail::bessel0
#include <felitronics/core/DetMath.h>
#include <felitronics/core/Math.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

//==============================================================================
// felitronics::codecgrid::BackResampler — takes a stretch of PCM BACK to the rate a codec ran at.
//
// A file at 44.1 kHz may have been coded at 48 kHz and converted afterwards (every Opus file is coded at 48 kHz;
// one corpus of "unmastered" 44.1 kHz files turned out to carry an AAC grid at 48). The grid is then not on the
// file's samples at all. Converted back with a ZERO-PHASE rational resampler it is again: output sample j sits at
// input time j * down / up exactly, so sample 0 stays sample 0 and a frame grid that was on whole samples at the
// codec's rate returns to whole samples. A twentieth of a sample off costs two thirds of the evidence, so the
// alignment is the whole point and the filter's quality is not: what is read afterwards is floored 80 dB under
// each coefficient's own level.
//
// THE FILTER IS THE ONE THE MEASUREMENT WAS MADE WITH, deliberately and exactly: a Kaiser (beta 14) windowed sinc
// of 20 * max (up, down) + 1 taps at the common rate, cut off at the lower of the two Nyquist frequencies and
// normalised to unity gain at DC — scipy.signal.resample_poly's default design with that window. It is short
// (about 21 input samples) and its transition is wide; the thresholds of this module were set on stretches that
// went through it, so a better filter here would be a different instrument.
//
// THIS IS NOT core::DeliveryResampler AND MUST NOT BECOME IT. That one is a streaming converter with a 140 dB
// stopband and a latency that is a real number of samples; this one is a one-shot over a few seconds, addressed
// by output index, with no latency at all.
//==============================================================================
namespace felitronics::codecgrid
{

class BackResampler
{
public:
    static constexpr double kKaiserBeta = 14.0;
    static constexpr int kHalfLengthPerRate = 10;

    // The lowest terms of rate conversion `from` -> `to` for the two rates this module pairs.
    // ok == false for any other pair.
    struct Ratio { bool ok = false; int up = 1, down = 1; };
    static constexpr Ratio ratioFor (int from, int to) noexcept
    {
        if (from == 44100 && to == 48000) return { true, 160, 147 };
        if (from == 48000 && to == 44100) return { true, 147, 160 };
        return {};
    }

    static constexpr std::size_t tapsFor (int up, int down) noexcept
    {
        return (std::size_t) (2 * kHalfLengthPerRate * (up > down ? up : down) + 1);
    }

    [[nodiscard]] bool prepare (int up, int down)
    {
        up_ = 0;
        if (up < 1 || down < 1 || up > 4096 || down > 4096) return false;
        const int rate = up > down ? up : down;
        const int half = kHalfLengthPerRate * rate, taps = 2 * half + 1;
        std::vector<double> h ((std::size_t) taps);
        const double denom = window::detail::bessel0 (kKaiserBeta), fc = 1.0 / (double) rate;
        double sum = 0.0;
        for (int n = 0; n < taps; ++n)
        {
            const double m = (double) (n - half);
            const double r = m / (double) half;                                      // -1 .. 1
            const double w = window::detail::bessel0 (kKaiserBeta * std::sqrt (1.0 - r * r)) / denom;
            const double a = core::kPi * fc * m;
            const double sinc = n == half ? 1.0 : core::det::sin (a) / a;
            h[(std::size_t) n] = w * sinc;
            sum += h[(std::size_t) n];
        }
        taps_.assign ((std::size_t) taps, 0.0f);
        for (int n = 0; n < taps; ++n) taps_[(std::size_t) n] = (float) ((double) up * h[(std::size_t) n] / sum);
        up_ = up;
        down_ = down;
        half_ = half;
        return true;
    }

    bool prepared() const noexcept { return up_ > 0; }

    // How many output samples `n` input samples make: ceil (n * up / down).
    std::int64_t outputLength (std::int64_t n) const noexcept
    {
        return up_ > 0 && n > 0 ? (n * (std::int64_t) up_ + (std::int64_t) down_ - 1) / (std::int64_t) down_ : 0;
    }

    // out[j - first] for j in [first, first + count): y[j] = sum_i x[i] h[j * down - i * up + half], the input
    // taken as zero outside [0, n).
    void resample (const float* x, int n, std::int64_t first, int count, float* out) const noexcept
    {
        if (up_ <= 0) return;
        const std::int64_t up = up_, down = down_, half = half_;
        const float* h = taps_.data();
        for (int c = 0; c < count; ++c)
        {
            const std::int64_t t = (first + (std::int64_t) c) * down;                 // the output's time, in common-rate samples
            std::int64_t lo = floorDiv (t - half + up - 1, up);                        // ceil ((t - half) / up)
            std::int64_t hi = floorDiv (t + half, up);
            if (lo < 0) lo = 0;
            if (hi > (std::int64_t) n - 1) hi = (std::int64_t) n - 1;
            float acc = 0.0f;
            for (std::int64_t i = lo; i <= hi; ++i) acc += x[i] * h[t - i * up + half];
            out[c] = acc;
        }
    }

private:
    static constexpr std::int64_t floorDiv (std::int64_t a, std::int64_t b) noexcept
    {
        const std::int64_t q = a / b;
        return (a % b != 0 && (a < 0) != (b < 0)) ? q - 1 : q;
    }

    int up_ = 0, down_ = 1, half_ = 0;
    std::vector<float> taps_;
};

} // namespace felitronics::codecgrid

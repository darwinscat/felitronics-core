// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

// An INDEPENDENT true-peak oracle for suites that must prove a bound rather than observe a number.
//
// WHY THIS EXISTS AT ALL. `felitronics::limiter::TruePeakLimiter` enforces its ceiling on an internal
// oversampled grid built by `oversampling::PolyphaseOversampler`. Measuring the result with a fresh
// instance of THAT SAME CLASS — which is what the limiter suite did, and what the migration plan's
// "measure with an independent oversampler" instruction was originally satisfied by — is not an
// independent measurement: it shares the prototype design, the Kaiser β, the 0.90 × Nyquist cutoff and
// the polyphase geometry with the thing under test. Measured on the witnesses in
// `truepeak_witnesses.h`, a fresh 8× polyphase reads 0.27 dB BELOW the true value on near-Nyquist
// transient material and 0.042 dB below it on EBU Tech 3341 test 16 — biased toward PASS, which is the
// one direction a ceiling proof cannot tolerate.
//
// THE PRIMARY ORACLE IS THE DEFINITION, NOT A BETTER FILTER. `truePeakDbFft()` interpolates by
// spectral zero-padding in double precision: for a buffer that begins and ends in silence, periodic
// sinc interpolation IS band-limited reconstruction, so there is no filter design to get wrong and no
// pass-band to droop. Measured against the five analytically-known EBU Tech 3341 Table 1 signals it
// reads the analytic amplitude to 9e-6 dB — four orders of magnitude below the tightest assertion any
// caller makes here. `selfCheckAgainstEbu()` runs exactly that comparison, so a drifted oracle fails
// loudly instead of silently moving everyone's budget.
//
// `truePeakDbSinc()` is the cross-check and deliberately shares nothing with either the primary oracle
// or the code under test: a different window (Blackman, not rectangular-in-frequency), a different
// cutoff (the base Nyquist, not 0.90 of it), a different evaluation (direct convolution on a
// sub-sample grid, not a transform). Two constructions that agree to 1e-4 dB are evidence; one
// construction is a hope.
//
// `truePeakDbPolyphase()` is kept ONLY so a suite can print the status-quo reading side by side with
// the truth and LABEL it. Every true-peak number in this repository is required to name the path that
// produced it, because the core holds two true-peak designs that disagree by up to 0.13 dB on nothing
// more than the sample-grid offset of one band-limited signal (see `analysis/TruePeakMeter.h`).
//
// OFFLINE, ALLOCATING, TEST-ONLY. This is not an RT meter and must never be used as one.

#include <felitronics/core/OfflineFft.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>
#include <utility>
#include <vector>

namespace felitronics::test::tp
{

inline constexpr double kPi = 3.14159265358979323846;

// Band-limited peak by spectral zero-padding, in dBFS/dBTP. `pad` is the interpolation factor.
//
// THE ORACLE IS ALSO A GRID, AND THAT BITES ON EXACTLY THE SIGNALS THIS FILE EXISTS FOR. Zero-padding
// evaluates the reconstruction at pad points per sample — finely, but discretely — so for a tone at a
// rational fs*p/q the oracle's own grid visits a finite phase set and can miss the crest by the same
// geometry the code under test is being measured for. Measured, not theorised: at pad = 8 a 4fs/11
// tone read 0.003 dB of excess where the truth is 0.089, because 11*8/gcd(4,88) = 22 puts the oracle's
// phases exactly where the 2x limiter's are. The oracle's own under-read is 20*log10(cos(pi/M_pad))
// with M_pad = q*pad/gcd(p, q*pad), so the default of 32 keeps it at or below 0.006 dB for every tone
// this repository asserts on — five times under the tightest tolerance any caller uses. Raise it, do
// not lower it, and never measure a rational tone at a pad that shares its denominator's factors.
//
// (For broadband material there is no such conspiracy: a click or noise has no single phase to hide in,
// and 16 is plentiful. The default is set for the harder case.)
//
// The Nyquist bin is SPLIT in half across +fs/2 and -fs/2. Placing it at only one of the two makes the
// interpolated signal complex, and taking the real part of that silently scales the Nyquist component
// by cos() of a phase nobody chose — an error that stays invisible on any signal without energy at
// exactly Nyquist, which is to say on every signal except the one that would catch it.
//
// The buffer MUST begin and end in silence (every witness here is tapered): spectral zero-padding
// interpolates the PERIODIC extension, so a buffer that ends mid-waveform gets a discontinuity at the
// wrap and the oracle would honestly report the overshoot of a step the signal never contained.
inline double truePeakDbFft (const std::vector<float>& x, int pad = 32)
{
    // `pad` is rounded UP to a power of two and floored at 2, because both are load-bearing rather than
    // tidiness: `fftInplace` is radix-2 and a pad of 3 would hand it a length that is not a power of two,
    // and at pad == 1 the two Nyquist destinations below collide — the second write lands on the first
    // and silently halves any energy at exactly Nyquist. Neither shows up on a signal without Nyquist
    // content, which is to say on every signal except the one that would catch it.
    if (x.size() < 2) return -std::numeric_limits<double>::infinity();
    pad = (int) core::offline::nextPow2 ((std::size_t) (pad < 2 ? 2 : pad));
    const std::size_t n = core::offline::nextPow2 (x.size());
    if (n == 0) return -std::numeric_limits<double>::infinity();
    std::vector<std::complex<double>> a (n, std::complex<double> {});
    for (std::size_t i = 0; i < x.size(); ++i) a[i] = (double) x[i];
    core::offline::fftInplace (a, -1);

    std::vector<std::complex<double>> b ((std::size_t) pad * n, std::complex<double> {});
    for (std::size_t k = 0; k < n / 2; ++k)            b[k] = a[k];                            // positive freqs
    for (std::size_t k = n / 2 + 1; k < n; ++k)        b[(std::size_t) pad * n - (n - k)] = a[k];   // negative freqs
    b[n / 2]                            = a[n / 2] * 0.5;
    b[(std::size_t) pad * n - n / 2]    = a[n / 2] * 0.5;
    core::offline::fftInplace (b, +1);

    const double sc = 1.0 / (double) n;               // fftInplace is unnormalised
    double mx = 0.0;
    for (const auto& v : b) mx = std::max (mx, std::fabs (v.real()) * sc);
    return 20.0 * std::log10 (mx > 1e-15 ? mx : 1e-15);
}

// Cross-check oracle: direct Blackman-windowed-sinc reconstruction in double, cutoff at the BASE
// Nyquist (not 0.90 of it — that guard band is the code-under-test's choice, not the truth's), taps
// 2*half+1, evaluated at `sub` positions per sample around the `topK` LARGEST local maxima of |x|.
// Restricting the search that way is safe here and not a shortcut: a band-limited signal's continuous
// maximum lies within one sample of a sample-grid local maximum, and it lies near the largest ones —
// an inter-sample peak hides at most ~4 dB behind its neighbouring sample, far less than the spread
// between the top of the signal and its 16th-loudest crest. It is also why this is the CROSS-check and
// the transform-based oracle above, which examines every position, is the judge.
inline double truePeakDbSinc (const std::vector<float>& x, int half = 512, int sub = 128, int topK = 16)
{
    if (x.empty()) return -std::numeric_limits<double>::infinity();
    const double n1 = 2.0 * (double) half;
    auto kern = [&] (double t) -> double
    {
        if (std::fabs (t) > (double) half) return 0.0;
        const double s = (std::fabs (t) < 1e-12) ? 1.0 : std::sin (kPi * t) / (kPi * t);
        const double u = (t + (double) half) / n1;
        return s * (0.42 - 0.5 * std::cos (2.0 * kPi * u) + 0.08 * std::cos (4.0 * kPi * u));
    };
    double smax = 0.0;
    for (float v : x) smax = std::max (smax, (double) std::fabs (v));
    double best = smax;
    const double th = smax * 0.5;
    std::vector<std::pair<double, int>> cand;
    for (int i = 1; i + 1 < (int) x.size(); ++i)
    {
        const double a = std::fabs ((double) x[(std::size_t) i]);
        if (a < th) continue;
        if (a < std::fabs ((double) x[(std::size_t) (i - 1)]) || a < std::fabs ((double) x[(std::size_t) (i + 1)])) continue;
        cand.emplace_back (a, i);
    }
    std::sort (cand.begin(), cand.end(), [] (const auto& l, const auto& r) { return l.first > r.first; });
    if ((int) cand.size() > topK) cand.resize ((std::size_t) topK);
    for (const auto& c : cand)
    {
        const int i = c.second;
        for (int s = -sub; s <= sub; ++s)
        {
            const double t = (double) s / (double) sub;
            double acc = 0.0;
            for (int k = -half; k <= half; ++k)
            {
                const int j = i + k;
                if (j < 0 || j >= (int) x.size()) continue;
                acc += (double) x[(std::size_t) j] * kern (t - (double) k);
            }
            best = std::max (best, std::fabs (acc));
        }
    }
    return 20.0 * std::log10 (best > 1e-15 ? best : 1e-15);
}

// The SAMPLE peak, for the one thing it is good for: separating a downsampler's overshoot (visible in
// the samples) from grid geometry (visible only between them).
inline double samplePeakDb (const std::vector<float>& x)
{
    double m = 0.0;
    for (float v : x) m = std::max (m, (double) std::fabs (v));
    return 20.0 * std::log10 (m > 1e-15 ? m : 1e-15);
}

} // namespace felitronics::test::tp

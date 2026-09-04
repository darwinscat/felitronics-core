// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

// Synthesized witnesses for proving a TRUE-PEAK CEILING, and the closed forms that say what each one
// is allowed to deliver. Written for `limiter::TruePeakLimiter`; the compressor work needs the same
// transient-dense material, which is why these live in test_support rather than in one suite.
//
// WHY SYNTHESIZED. The material that broke the predecessor (ffmpeg `alimiter` shipping +2.7 dB above
// its own ceiling on a click train under ~+88 dB of makeup) is gitignored audio in another repository,
// and this repo must not carry audio at all. The precedent is `ebu_tech3341_truepeak.h`, which
// synthesizes the EBU Tech 3341 signals instead of shipping 87 MB of third-party goldens and
// reproduces the official readings to 0.0005 dB. The taper here is that header's `fadeWindow()` — one
// copy of the rule, not two, because the fixture bug this project has already paid for three times is
// two copies of a specification drifting apart.
//
// ============================================================================================
// WHAT THE MEASUREMENTS SAID, AND WHERE THE TASK'S OWN PREMISE WAS WRONG
// ============================================================================================
// The brief for this fixture asked for "content above 0.8 x Nyquist, where both interpolators are
// weakest". That is true of the METERS and false of the LIMITER, and building the witness on it would
// have produced a test that passes for the wrong reason. The limiter up- and downsamples through its
// own Kaiser prototype cut at 0.90 x Nyquist, so its response in BASE-rate terms is 0.00 dB at
// 0.36 fs, -0.40 at 0.40, -4.0 at 0.44, -6.0 at 0.45 (identical at every factor — the transition width
// is constant in base-rate units). Content up there is therefore ATTENUATED ON THE WAY IN, cannot
// reach the ceiling on the way out, and measures 6 dB BELOW the closed form rather than above it:
// a 0.45 fs tone hot enough to demand 13 dB of limiting was delivered at -6.0 dB re ceiling. The
// witness that actually finds the ceiling's weak point is a tone at a SMALL-DENOMINATOR fraction of
// the rate, well inside the pass band, because what matters is not how high the frequency is but HOW
// FEW DISTINCT PHASES it visits on the limiter's internal grid.
//
// THE GEOMETRY LAW (derived, then measured to 0.002 dB — see gridPhaseCount/gridBreachDb below).
// The limiter's bound is exact on its own F x fs grid. A steady tone at f = fs*p/q visits
// M = q*F/gcd(p, q*F) distinct phases on that grid. The limiter compares MAGNITUDES, which folds the
// phase circle onto a half-turn: for even M the folded phases sit 2*pi/M apart, for odd M they sit
// pi/M apart. The crest can hide half a spacing away, so the grid under-reads it by cos(pi/M') with
// M' = M for even M and 2M for odd M — and the delivered true peak exceeds the ceiling by exactly
// that. Measured: fs/3 -> +1.2499 / +0.3020 / +0.0757 dB at F = 2 / 4 / 8 against a closed form of
// +1.2494 / +0.3011 / +0.0746. The residual is at most 0.0011 dB.
//
// (Two wrong forms were rejected on measurement, and both are easy to re-derive by accident. The
// half-step form -20*log10(cos(pi*f/(F*fs))) is the p/gcd == 1 special case only; it predicts
// -2.38 dB at 0.45 fs and F=2 where the truth is -0.027 before the transition band even applies.
// Dropping the odd-M doubling predicts +1.94 dB for 2fs/5 at F=2 where the measurement, corrected for
// the 0.40 fs droop, gives +0.44.)
// ============================================================================================

#include "ebu_tech3341_truepeak.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <vector>

namespace felitronics::test::tpw
{

inline constexpr double kPi = 3.14159265358979323846;

// The denominator M' of the breach formula for a tone at fs*p/q on an F*fs grid. NOT literally the
// number of distinct folded phase positions — for odd M that count is M and this returns 2M, because
// what the formula needs is the HALF-spacing, and folding an odd-M phase set halves its spacing.
inline int gridPhaseCount (int p, int q, int factor)
{
    const int m = q * factor / std::gcd (p, q * factor);
    return (m % 2 == 0) ? m : 2 * m;
}

// How far above a ceiling enforced on the F*fs grid the DELIVERED true peak of that tone may legally
// sit. Derived from the factor and the tone, not pinned: a redesign that changes the factor moves this
// assertion with it, the way ebu3341::gridBoundDb() takes the factor from the meter under test.
inline double gridBreachDb (int p, int q, int factor)
{
    return -20.0 * std::log10 (std::cos (kPi / (double) gridPhaseCount (p, q, factor)));
}

// ============================================================================================
// THE SECOND MECHANISM, AND THE ONE THAT DOES NOT GO AWAY BY OVERSAMPLING HARDER
// ============================================================================================
// Grid geometry above is only half the delivered excess, and it is the half that shrinks with the
// factor. The other half comes from the gain sequence itself: the attack is instantaneous, so the
// limited product is NOT band-limited even when the input was, and the downsampler re-band-limits it
// on the way out with a step overshoot no grid argument covers.
//
// It behaves nothing like the grid term, and two of the three things one would GUESS about it are
// false. Delivered excess in dB at a -1 dBTP ceiling, measured with the FFT oracle:
//
//   depth (white noise, peak dBFS)      -1     +2     +5    +10    +20      release
//     F=4                             0.435  0.678  0.919  0.920  0.921      50 ms
//     F=8                             0.332  0.699  0.867  0.868  0.869      50 ms
//     F=4                             0.197  0.435  0.435  0.919  0.920     100 ms
//
//   spectral tilt (one-pole, at +10 dBFS)   flat   0.5    0.8    0.9   0.95   0.98
//     F=4, release 50 ms                    0.920  0.320  0.174  0.100  0.090  0.043
//     F=4, release  1 ms                    0.989  0.686  0.344  0.241  0.163  0.118
//
//   * NOT bought by oversampling. 8x costs twice the work, cuts the grid term by four, and moves this
//     by 0.05 dB.
//   * NOT bought by a slow release either, which is the guess that cost the most: the term SATURATES
//     near 0.92 (4x) / 0.87 (8x) once there is about 5 dB of gain reduction, and it does that at a
//     100 ms release as readily as at 1 ms. Release only helps while the limiter is barely working.
//   * It IS bought by spectral tilt, and steeply: the same depth on roughly music-like material
//     (one-pole 0.9) costs 0.10 dB where flat noise costs 0.92. So the worst case here is a bright,
//     dense, hard-limited master — not an ordinary one.
// A single hard EDGE, by contrast, is nearly free (0.004 dB at a 100 ms release, 0.04 at 1 ms): it is
// the continuous re-modulation of dense material that costs, not any one transition.
//
// Hence ONE envelope rather than the release-split pair an earlier draft carried — that split asserted
// a dependence the measurements above show does not exist, and it passed only because its dense rows
// sat at about 1 dB of gain reduction. It is an ENVELOPE over the witnesses in this header at
// lookahead >= 1 ms and release > 0, not a theorem; the degenerate parameters (either one at zero)
// leave it behind and are pinned separately in the suite.
inline constexpr double kModulationEnvelopeDb = 1.15;   // worst measured component: 1.06 (8x, click at 0.1 ms)

// The whole delivered allowance for a factor: the grid's closed form at the worst tone the round trip
// passes flat, plus the modulation envelope. (Content above ~0.40 fs has a LARGER grid term — a
// 0.35-0.45 fs transient sits 0.43 dB over on the 4x grid against the tone term's 0.30 — but the
// prototype's own droop removes it before delivery, so it cannot reach the output. Widen the pass band
// and this term has to be re-derived.)
inline double deliveredBudgetDb (int factor)
{
    return gridBreachDb (1, 3, factor) + kModulationEnvelopeDb;
}

// The taper, borrowed rather than re-spelled. 10 ms at both ends, raised cosine, indexed from the end
// on the way out so the two ramps mirror exactly (see the EBU header for why that distinction is not
// cosmetic). Both oracles in `truepeak_oracle.h` REQUIRE the ends to be silent.
inline double taper (long long i, long long n, double sampleRate)
{
    return ebu3341::fadeWindow (i, n, (long long) std::llround (sampleRate * 0.010));
}

//==============================================================================
// W1 — the click train. The predecessor's condition: dense transients, a huge crest factor, and a
// monstrous makeup gain. One click is a band-limited impulse with a flat spectrum from DC to
// `bandEdge` cycles/sample, windowed to +-`support` samples, and centred `offset` samples off the
// integer grid. At offset 0.5 the sample grid lands as far from the crest as it can: the sample peak
// falls ~3.1 dB below the true peak, so a limiter that bounds SAMPLES (which is precisely what
// `alimiter` does) ships that 3.1 dB straight to the file.
//
// The 0.45 default band edge is the widest that survives the limiter's own pass band; going wider
// only adds content the round trip removes, which is why the "above 0.8 x Nyquist" brief would have
// measured the filter instead of the ceiling.
inline std::vector<float> clickTrain (double sampleRate, double seconds, double periodMs,
                                      double peakDb, double offset = 0.5,
                                      double bandEdge = 0.45, int support = 256)
{
    const long long n = (long long) std::llround (sampleRate * seconds);
    std::vector<float> x ((std::size_t) (n < 0 ? 0 : n), 0.0f);
    if (n <= 0) return x;
    const double amp  = std::pow (10.0, peakDb / 20.0);
    const long long step = std::max (1LL, (long long) std::llround (sampleRate * periodMs * 0.001));
    for (long long c = step; c + step < n; c += step)
        for (long long i = c - support; i <= c + support; ++i)
        {
            if (i < 0 || i >= n) continue;
            const double t = (double) i - ((double) c + offset);
            const double s = (std::fabs (t) < 1e-12) ? 1.0
                           : std::sin (2.0 * kPi * bandEdge * t) / (2.0 * kPi * bandEdge * t);
            const double u = t / (double) support;                       // Blackman envelope on the pulse
            const double w = std::fabs (u) >= 1.0 ? 0.0
                           : (0.42 + 0.5 * std::cos (kPi * u) + 0.08 * std::cos (2.0 * kPi * u));
            x[(std::size_t) i] += (float) (amp * s * w);
        }
    for (long long i = 0; i < n; ++i) x[(std::size_t) i] = (float) (x[(std::size_t) i] * taper (i, n, sampleRate));
    return x;
}

//==============================================================================
// W2 — the grid-worst tone. `phaseTurns` is scanned by the caller: the crest must sit half a folded
// phase step from the limiter's grid, and the internal grid is offset from the input grid by half an
// oversampled sample (the prototype's group delay (N-1)/2 is a half-integer), so the worst phase is
// not derivable from p/q alone and is FOUND, per the fixture rule that a grid landing on the crest
// makes every implementation read the truth and the test blind.
inline std::vector<float> gridTone (double sampleRate, double seconds, int p, int q,
                                    double peakDb, double phaseTurns)
{
    const long long n = (long long) std::llround (sampleRate * seconds);
    std::vector<float> x ((std::size_t) (n < 0 ? 0 : n), 0.0f);
    if (n <= 0) return x;
    const double amp = std::pow (10.0, peakDb / 20.0);
    const double phi = 2.0 * kPi * phaseTurns;
    for (long long i = 0; i < n; ++i)
        x[(std::size_t) i] = (float) (taper (i, n, sampleRate) * amp
                                      * std::sin (2.0 * kPi * (double) p * (double) i / (double) q + phi));
    return x;
}

//==============================================================================
// W3 — dense broadband. The witness for the OTHER mechanism: the gain sequence has an instant attack,
// so the limited product is not band-limited, and the downsampler re-band-limits it on the way out
// with a 9.2 % step overshoot that no grid argument covers. A deterministic LCG so a failure is
// reproducible; one-pole lowpassed when `lpCoef` > 0 so the suite can carry a "realistic master" row
// next to the flat-spectrum worst case.
// `peakDb` is the amplitude the uniform generator reaches, i.e. very nearly the SAMPLE peak — it is not
// an RMS figure (uniform noise sits 4.77 dB below its peak in RMS, and calling this rmsDb, as an earlier
// draft did, would have described every row 4.77 dB quieter than it is and made a heavily-limited case
// look like a bypass).
inline std::vector<float> denseNoise (double sampleRate, double seconds, double peakDb,
                                      std::uint64_t seed = 0x9E3779B97F4A7C15ULL, double lpCoef = 0.0)
{
    const long long n = (long long) std::llround (sampleRate * seconds);
    std::vector<float> x ((std::size_t) (n < 0 ? 0 : n), 0.0f);
    if (n <= 0) return x;
    const double amp = std::pow (10.0, peakDb / 20.0);
    std::uint64_t s = seed ? seed : 1ULL;
    double z = 0.0;
    for (long long i = 0; i < n; ++i)
    {
        s = s * 6364136223846793005ULL + 1442695040888963407ULL;
        const double u = (double) ((s >> 33) & 0x7FFFFF) / (double) 0x400000 - 1.0;   // [-1,1)
        z = lpCoef > 0.0 ? (lpCoef * z + (1.0 - lpCoef) * u) : u;
        const double g = lpCoef > 0.0 ? 1.0 / std::sqrt ((1.0 - lpCoef) / (1.0 + lpCoef)) : 1.0;
        x[(std::size_t) i] = (float) (taper (i, n, sampleRate) * amp * z * g);
    }
    return x;
}

//==============================================================================
// W4 — the hard edge. A tone already at the ceiling, then a step to `stepDb` above it. Band-limited
// pulses cannot make the gain move discontinuously — their own tails pre-announce them inside the
// lookahead window — so they under-test the modulation mechanism. A step does.
inline std::vector<float> hardEdge (double sampleRate, double seconds, double toneHz,
                                    double baseDb, double stepDb, double stepAtSeconds)
{
    const long long n = (long long) std::llround (sampleRate * seconds);
    std::vector<float> x ((std::size_t) (n < 0 ? 0 : n), 0.0f);
    if (n <= 0) return x;
    const long long at = (long long) std::llround (sampleRate * stepAtSeconds);
    const double a0 = std::pow (10.0, baseDb / 20.0), a1 = std::pow (10.0, (baseDb + stepDb) / 20.0);
    for (long long i = 0; i < n; ++i)
        x[(std::size_t) i] = (float) (taper (i, n, sampleRate) * (i < at ? a0 : a1)
                                      * std::sin (2.0 * kPi * toneHz * (double) i / sampleRate));
    return x;
}

//==============================================================================
// W5 — the plateau. A long stretch held far over the ceiling at a low frequency, ending abruptly, with
// the edge placed at a sub-sample offset the caller sweeps. This is the adversary for the downsampler's
// SIGNED step response: the limited product sits flat near the ceiling and then stops, and the decimation
// FIR's overshoot (9.2 %, intrinsic to a sharp cutoff — a Kaiser window does not remove it) lands on top
// of the ceiling rather than on top of a decaying tail. It is the witness that separates "the release is
// too fast" from "the gain path itself overshoots": at a release of 0 it delivers +1.48 dB over, and at
// 0.1 ms the same construction gives +0.69.
inline std::vector<float> plateau (double sampleRate, double seconds, double toneHz, double peakDb,
                                   double fromSeconds, double toSeconds, double edgeOffset)
{
    const long long n = (long long) std::llround (sampleRate * seconds);
    std::vector<float> x ((std::size_t) (n < 0 ? 0 : n), 0.0f);
    if (n <= 0) return x;
    const double amp = std::pow (10.0, peakDb / 20.0);
    const double a = sampleRate * fromSeconds + edgeOffset, b = sampleRate * toSeconds + edgeOffset;
    for (long long i = 0; i < n; ++i)
    {
        const double t = (double) i;
        const double v = (t > a && t < b) ? amp * std::sin (2.0 * kPi * toneHz * t / sampleRate + 1.0) : 0.0;
        x[(std::size_t) i] = (float) (v * taper (i, n, sampleRate));
    }
    return x;
}

//==============================================================================
// FIXTURE SELF-CHECKS. Every assertion in the suites below looks at a MAXIMUM, so a malformed
// envelope — the fixture failure this project has already paid for — survives all of them. These do
// not look at maxima.
// Both oracles interpolate the PERIODIC extension, so what they need is a wrap that joins silence to
// silence smoothly. Two things are therefore checked, and the second is not "guard samples of exact
// zero" — a correct raised-cosine taper reaches zero only AT the endpoint (with zero derivative, which
// is what actually makes the wrap smooth), so demanding a run of exact zeros would reject the right
// fixture and invite someone to "fix" the taper into a worse one.
//   * the endpoints are exactly zero — this is the check that catches the off-by-one fade-out the EBU
//     header describes, where the last sample sits at 1.1e-5 of amplitude and the whole ramp is shifted;
//   * the guard interval is at least 60 dB down — this catches a grossly malformed envelope, the kind
//     that ends at a tenth of full scale and makes the oracle measure its own wrap.
inline bool endsAreSilent (const std::vector<float>& x, int guard = 8)
{
    if ((int) x.size() < 4 * guard) return false;
    if (x.front() != 0.0f || x.back() != 0.0f) return false;
    double pk = 0.0;
    for (float v : x) pk = std::max (pk, (double) std::fabs (v));
    const double lim = pk * 1.0e-3;
    for (int i = 0; i < guard; ++i)
        if (std::fabs ((double) x[(std::size_t) i]) > lim
            || std::fabs ((double) x[x.size() - 1 - (std::size_t) i]) > lim) return false;
    return true;
}

inline double crestFactorDb (const std::vector<float>& x)
{
    if (x.empty()) return 0.0;
    double pk = 0.0, s2 = 0.0;
    for (float v : x) { pk = std::max (pk, (double) std::fabs (v)); s2 += (double) v * (double) v; }
    const double rms = std::sqrt (s2 / (double) x.size());
    return 20.0 * std::log10 ((pk > 1e-30 ? pk : 1e-30) / (rms > 1e-30 ? rms : 1e-30));
}

} // namespace felitronics::test::tpw

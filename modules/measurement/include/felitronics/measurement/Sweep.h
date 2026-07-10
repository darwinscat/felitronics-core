// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

//==============================================================================
// felitronics::measurement — exponential sine sweep (ESS / Farina) + matched inverse filter.
// OFFLINE, MESSAGE-THREAD-ONLY (double, allocates). See Convolve.h for the precision rationale.
//
//   Sweep:   x(t) = sin( (w1*T/R) * (exp((t/T)*R) - 1) ),  t in [0,T],  R = ln(w2/w1)
//   Inverse: time-reversed sweep with a +6 dB/oct amplitude envelope so that x ⊛ inv ≈ δ (whitened),
//            normalized so the resulting δ has unit peak → deconvolution has unity gain.
//
// The k-th harmonic response of a memoryless nonlinearity appears time-ADVANCED by dt_k = L*ln(k)
// (L = T/R) relative to the linear IR — this lets the linear IR be windowed out clean (see Deconvolve.h).
//==============================================================================

#include <felitronics/measurement/Convolve.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace felitronics::measurement
{

struct SweepSpec
{
    double f1              = 20.0;      // start freq, Hz
    double f2              = 20000.0;   // stop freq, Hz
    double durationSeconds = 5.0;       // sweep duration, s
    double sampleRate      = 48000.0;   // Hz
    double amplitude       = 0.5;       // linear peak (~ -6 dBFS)
    double fadeSeconds     = 0.05;      // raised-cosine fade in/out, s
    double tailSeconds     = 1.0;       // trailing silence, s (room/cab decay capture)
};

struct Sweep
{
    std::vector<double> signal;         // played sweep (sweep proper + tail silence)
    std::vector<double> inverse;        // Farina inverse filter (length == sweepLen)
    SweepSpec           spec;           // AS-SANITIZED (clamped) — the truth the caller must trust
    std::size_t         sweepLen = 0;   // samples of the sweep proper (excludes tail)
    double              harmonicL = 0.0; // L = T/R — the harmonic-advance time constant
};

// Build the sweep + matched inverse. Clamps/heals every parameter (house rule: reject non-finite,
// clamp to a valid finite domain) — guards f1=0→R=∞, f1==f2→R=0, dur=0, f2>Nyquist. The returned
// `spec` reflects the clamps.
inline Sweep makeSweep (SweepSpec s)
{
    if (! (s.sampleRate >= 3000.0)) s.sampleRate = 48000.0;   // floor first so clamp lo<=hi holds
    const double nyq = 0.5 * s.sampleRate;
    if (! (s.f1 > 0.0)) s.f1 = 20.0;
    s.f1 = std::clamp (s.f1, 1.0, nyq - 2.0);
    if (! (s.f2 > s.f1)) s.f2 = s.f1 * 10.0;
    s.f2 = std::clamp (s.f2, s.f1 * 1.0001 + 1.0, nyq * 0.999);
    if (! (s.durationSeconds > 0.0)) s.durationSeconds = 1.0;
    s.durationSeconds = std::max (s.durationSeconds, 16.0 / s.sampleRate);
    if (! (s.amplitude >= 0.0)) s.amplitude = 0.5;
    s.fadeSeconds = std::clamp (s.fadeSeconds, 0.0, 0.49 * s.durationSeconds);
    s.tailSeconds = std::max (s.tailSeconds, 0.0);

    Sweep out;
    out.spec = s;
    const double w1 = 2.0 * core::kPi * s.f1;
    const double w2 = 2.0 * core::kPi * s.f2;
    const double T  = s.durationSeconds;
    const double R  = std::log (w2 / w1);
    const std::size_t N     = (std::size_t) std::llround (T * s.sampleRate);
    const std::size_t Ntail = (std::size_t) std::llround (s.tailSeconds * s.sampleRate);
    const std::size_t Nf    = (std::size_t) std::llround (s.fadeSeconds * s.sampleRate);
    out.sweepLen  = N;
    out.harmonicL = T / R;

    std::vector<double> x (N + Ntail, 0.0);
    for (std::size_t n = 0; n < N; ++n)
    {
        const double t   = (double) n / s.sampleRate;
        const double phi = (w1 * T / R) * (std::exp ((t / T) * R) - 1.0);
        double env = 1.0;
        if (Nf > 0)
        {
            if (n < Nf)             env = 0.5 * (1.0 - std::cos (core::kPi * (double) n / (double) Nf));
            else if (n >= N - Nf)   env = 0.5 * (1.0 - std::cos (core::kPi * (double) (N - 1 - n) / (double) Nf));
        }
        x[n] = s.amplitude * std::sin (phi) * env;
    }
    out.signal = std::move (x);

    // Inverse: reverse the sweep, apply +6 dB/oct boost (envelope exp((t/T)*R)).
    std::vector<double> inv (N);
    for (std::size_t n = 0; n < N; ++n)
    {
        const double tfrac = (double) n / (double) N;
        inv[N - 1 - n] = out.signal[n] * std::exp (tfrac * R);
    }

    // Normalize so (sweep_proper ⊛ inv) has unit peak → deconvolution is unity gain.
    std::vector<double> xp (out.signal.begin(), out.signal.begin() + (std::ptrdiff_t) N);
    const auto g = convolve (xp, inv);
    double peak = 0.0;
    for (double v : g) peak = std::max (peak, std::fabs (v));
    const double scale = (peak > 0.0) ? 1.0 / peak : 1.0;
    for (double& v : inv) v *= scale;
    out.inverse = std::move (inv);
    return out;
}

// Time-advance (in samples) of the k-th harmonic image relative to the linear IR. k <= 1 → 0.
inline double harmonicAdvanceSamples (const Sweep& sw, int k) noexcept
{
    return (k > 1) ? sw.harmonicL * std::log ((double) k) * sw.spec.sampleRate : 0.0;
}

} // namespace felitronics::measurement

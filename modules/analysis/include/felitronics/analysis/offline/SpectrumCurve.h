// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

//==============================================================================
// felitronics::analysis::offline — OFFLINE, MESSAGE-THREAD-ONLY display/analysis curves (allocates,
// double FFT; NOT RT-safe — distinct from the rest of `analysis`, which is RT metering). The RT path
// uses the float `felitronics::core::fft` seam; these take a whole IR and are called off the audio
// thread to draw a spectrum. Extracted from OrbitCapture's SpectrumView.
//
//   logMagnitudeCurve — a 1/N-octave RMS-power-smoothed magnitude curve on a log-f grid, in dB.
//                       The cabinet's EQ curve at a glance.
//   interferenceDb    — where a multi-mic blend cancels (phase eats the frequency) or reinforces,
//                       relative to an incoherent (phase-blind) power sum of the mics.
//
// Both use the shared offline double FFT (felitronics::core::offline). RMS-power averaging is
// energy-preserving: 10·log10(mean(|X|²)) == 20·log10(rms|X|).
//==============================================================================

#include <felitronics/core/OfflineFft.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <vector>

namespace felitronics::analysis::offline
{

// Controls for logMagnitudeCurve. Defaults reproduce OrbitCapture's 20 Hz–20 kHz, 256-point,
// 1/12-octave curve. fHi is clamped to Nyquist internally.
struct LogCurveSpec
{
    double      fLo             = 20.0;      // grid low edge, Hz (> 0)
    double      fHi             = 20000.0;   // grid high edge, Hz (clamped to Nyquist)
    int         points          = 256;       // grid points (>= 2; smaller is clamped up)
    double      octaveDivisions = 12.0;      // 1/N-octave smoothing band (12 → 1/12-oct); > 0
    std::size_t minNfft         = 8192;      // FFT floor so short IRs still get fine bins
    bool        normalize       = true;      // subtract the peak → 0 dB max (skipped for a silent curve)
};

// One mic's magnitude curve + whether it counts in the blend, for interferenceDb. `db` must be a
// logMagnitudeCurve computed with normalize=false (shared absolute reference).
struct MicCurveView
{
    std::span<const double> db;
    bool                    audible = true;
};

namespace detail
{
template <class T>
inline std::vector<double> logMagCurveImpl (std::span<const T> ir, double sr, LogCurveSpec s)
{
    if (ir.size() < 8 || ! (sr > 0.0) || ! std::isfinite (sr)) return {};
    s.points = std::max (2, s.points);
    if (! (s.fLo > 0.0) || ! (s.fHi > s.fLo)) return {};          // need an ascending, positive grid
    if (! (s.octaveDivisions > 0.0) || ! std::isfinite (s.octaveDivisions)) s.octaveDivisions = 12.0;   // heal <=0 / NaN / inf
    s.octaveDivisions = std::clamp (s.octaveDivisions, 1.0, 96.0); // finite band → halfBand never inf (an inf→int cast is UB)
    const double nyq = 0.5 * sr;
    s.fHi = std::min (s.fHi, nyq);                                // never ask for bins past Nyquist
    if (! (s.fHi > s.fLo)) return {};                            // Nyquist collapsed the grid (sr too low)

    std::vector<double> x (ir.size());                           // to double + neutralize non-finite input
    for (std::size_t i = 0; i < ir.size(); ++i) { const double v = (double) ir[i]; x[i] = std::isfinite (v) ? v : 0.0; }

    const std::size_t np2 = core::offline::nextPow2 (x.size());
    if (np2 == 0) return {};                                      // absurd IR length overflowed pow2
    const std::size_t nfft = core::offline::nextPow2 (std::max (s.minNfft, np2));  // pow2 → binHz matches magSpectrum's actual transform
    if (nfft == 0) return {};                                     // minNfft overflowed pow2
    const auto M = core::offline::magSpectrum (x, nfft);
    if (M.size() < 2) return {};
    const double binHz    = sr / (double) nfft;
    const double halfBand = std::pow (2.0, 1.0 / (2.0 * s.octaveDivisions));
    const int    lastBin  = (int) M.size() - 1;

    std::vector<double> curve ((std::size_t) s.points);
    for (int i = 0; i < s.points; ++i)
    {
        const double f = s.fLo * std::pow (s.fHi / s.fLo, (double) i / (double) (s.points - 1));
        int blo = (int) std::floor (f / halfBand / binHz);
        int bhi = (int) std::ceil  (f * halfBand / binHz);
        blo = std::min (std::max (1, blo), lastBin);              // skip DC bin (0); never past the top bin (no inverted band)
        bhi = std::min (lastBin, std::max (bhi, blo));
        double sum = 0.0; int c = 0;
        for (int b = blo; b <= bhi; ++b) { const double m = M[(std::size_t) b]; sum += m * m; ++c; }
        curve[(std::size_t) i] = (c > 0 && std::isfinite (sum)) ? 10.0 * std::log10 (sum / (double) c + 1e-20) : -200.0;
    }
    if (s.normalize)
    {
        double pk = std::numeric_limits<double>::lowest();
        for (double v : curve) pk = std::max (pk, v);
        if (std::isfinite (pk) && pk > -199.0)                    // a silent curve would flat-line at 0 dB — leave it at the floor
            for (double& v : curve) v -= pk;
    }
    return curve;
}
} // namespace detail

// 1/N-octave RMS-power-smoothed magnitude of an IR on a log-f grid, in dB. OFFLINE (allocates, double
// FFT). Returns spec.points values (empty on a degenerate request: <8 samples, sr<=0/NaN, fLo<=0,
// fHi<=fLo). Non-finite input samples are treated as 0. See LogCurveSpec.
inline std::vector<double> logMagnitudeCurve (std::span<const float> ir, double sampleRate, LogCurveSpec spec = {})
{ return detail::logMagCurveImpl<float> (ir, sampleRate, spec); }

inline std::vector<double> logMagnitudeCurve (std::span<const double> ir, double sampleRate, LogCurveSpec spec = {})
{ return detail::logMagCurveImpl<double> (ir, sampleRate, spec); }

// Interference in dB per grid point: coherentBlendDb[p] − 10·log10(Σ incoherent mic power). Negative
// = phase cancellation eats the frequency; positive = reinforcement. Returns 0 dB where no mic is
// audible (no signal → no interference). Output length = coherentBlendDb.size(); point p reads
// mic.db[p] where present.
//
// PRECONDITIONS: every input curve (coherent + each mic) MUST be a logMagnitudeCurve with
// normalize=FALSE — one shared absolute reference; a normalized curve makes the subtraction
// meaningless. Mic curves must carry the SAME per-mic gains/filters used to build the coherent blend.
// The result is PURE dB — any display gating (fading the tint in a filter's stopband, mapping to an
// alpha) is the caller's concern (multiplying a dB delta by a linear [0,1] gate is not valid dB math).
inline std::vector<double> interferenceDb (std::span<const double>       coherentBlendDb,
                                           std::span<const MicCurveView> mics)
{
    std::vector<double> out (coherentBlendDb.size(), 0.0);
    for (std::size_t p = 0; p < coherentBlendDb.size(); ++p)
    {
        double pw = 0.0; int nAud = 0;
        for (const auto& m : mics)
            if (m.audible && p < m.db.size()) { pw += std::pow (10.0, m.db[p] / 10.0); ++nAud; }
        out[p] = (nAud > 0 && pw > 1e-19) ? coherentBlendDb[p] - 10.0 * std::log10 (pw) : 0.0;
    }
    return out;
}

} // namespace felitronics::analysis::offline

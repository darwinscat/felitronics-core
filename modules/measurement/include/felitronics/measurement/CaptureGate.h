// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

//==============================================================================
// felitronics::measurement — input-validation gate BEFORE deconvolution. Garbage-in (clipped / no
// sweep / non-finite) makes plausible-but-wrong IRs that pollute the catalog. Because a capture tool
// OWNS the sweep + sample rate, the checks are cheap and reliable. OFFLINE, MESSAGE-THREAD-ONLY.
//
// User-facing reject STRINGS stay app-side (localization/tone don't belong in core) — the gate returns
// a GateReject enum. `sweepConfidence` is a PRESENCE ratio (deconv peak / early floor), NOT an SNR;
// the honest band-limited measurement SNR is `snrDb`. Thresholds ship as STARTING POINTS — the real
// values must be calibrated on hardware (mic self-noise, room).
//==============================================================================

#include <felitronics/measurement/Convolve.h>
#include <felitronics/measurement/PeakClip.h>
#include <felitronics/measurement/Sweep.h>

#include <cmath>
#include <cstddef>
#include <span>
#include <vector>

namespace felitronics::measurement
{

struct GateConfig
{
    double clipLevel          = 0.999;    // full-scale threshold for the flat-top run
    int    clipRunSamples     = 3;        // consecutive full-scale samples that count as real saturation
    double minSweepConfidence = 1000.0;   // CALIBRATION STARTING POINT (see file header)
    double bandLoHz           = 100.0;    // honest-SNR band — guitar-cab default; a room tool passes 20/20k
    double bandHiHz           = 6000.0;
};

enum class GateReject { None, EmptyRecording, NonFinite, Clipped, SweepNotDetected };

struct GateReport
{
    bool        ok             = true;
    GateReject  reason         = GateReject::None;
    bool        clipped        = false;
    double      peakDbfs       = -120.0;
    int         clipRun        = 0;       // longest flat-top run at full scale
    bool        sweepPresent   = false;
    double      sweepConfidence = 0.0;    // deconv peak / early floor — a PRESENCE ratio, NOT an SNR
    double      snrDb          = 0.0;     // honest band-limited measurement SNR (what REW reports)
    bool        snrValid       = false;   // false when the recording was too short to measure SNR (snrDb is meaningless)
    std::size_t sweepLag       = 0;
};

// Band-limited [loHz, hiHz] energy of a window (sum of |X[k]|² in-band). The guitar-cab working band by default.
inline double bandEnergy (std::span<const double> x, double sr, double loHz, double hiHz)
{
    if (x.empty() || ! (sr > 0.0)) return 0.0;
    std::size_t nf = 1; while (nf < x.size() && nf != 0) nf <<= 1; nf = std::max<std::size_t> (nf, 4096);
    const auto M = magSpectrum (x, nf);
    const double binHz = sr / (double) nf, nyq = 0.5 * sr;
    // Clamp the band to a finite [0, Nyquist] BEFORE the float→size_t cast (casting a negative / NaN / inf
    // double to an unsigned is UB; a raw negative loHz would wrap to a huge index).
    const double lof = std::clamp (std::isfinite (loHz) ? loHz : 0.0, 0.0, nyq);
    const double hif = std::clamp (std::isfinite (hiHz) ? hiHz : nyq, 0.0, nyq);
    const std::size_t lo = (std::size_t) std::ceil (lof / binHz), hi = (std::size_t) std::floor (hif / binHz);
    double e = 0.0;
    for (std::size_t i = lo; i <= hi && i < M.size(); ++i) e += M[i] * M[i];
    return e;
}

// Reject a recording that must NOT reach deconvolution. clipped = a flat-top RUN at full scale (real
// saturation), NOT mere peak proximity, so a hot-but-clean −0.1 dBFS take is not false-rejected.
inline GateReport gateRecording (std::span<const double> rec, const Sweep& sw, const GateConfig& cfg = {})
{
    GateReport g;
    if (rec.empty()) { g.ok = false; g.reason = GateReject::EmptyRecording; return g; }

    // House rule: reject non-finite BEFORE any math (else one NaN makes the deconv all-NaN → every
    // comparison false → a wrong "sweep not detected" verdict). Honest failure with the true reason.
    // The scan itself is the shared standalone pass (PeakClip.h).
    const auto pc = scanPeakClip (rec, { cfg.clipLevel, cfg.clipRunSamples });
    if (pc.nonFinite) { g.ok = false; g.reason = GateReject::NonFinite; return g; }
    g.peakDbfs = pc.peakDbfs;
    g.clipRun  = pc.clipRun;
    g.clipped  = pc.clipped;

    const auto d = convolve (rec, sw.inverse);
    if (d.empty()) { g.ok = false; g.reason = GateReject::SweepNotDetected; return g; }   // malformed sweep (empty inverse)
    double dpk = 0.0; std::size_t dl = 0;
    for (std::size_t i = 0; i < d.size(); ++i)
    {
        const double a = std::fabs (d[i]);
        if (a > dpk) { dpk = a; dl = i; }
    }
    // floor from an early window silent for a real sweep (only ~f2/f1-order harmonics reach index 0).
    const std::size_t fw = std::min<std::size_t> (2000, d.size() / 8 + 1);
    double fe = 0.0;
    for (std::size_t i = 0; i < fw; ++i) fe += d[i] * d[i];
    const double floor = std::sqrt (fe / (double) std::max<std::size_t> (1, fw)) + 1e-12;
    g.sweepConfidence = dpk / floor;
    g.sweepLag        = dl;
    g.sweepPresent    = (g.sweepConfidence >= cfg.minSweepConfidence);

    // HONEST measurement SNR on the RAW recording (no deconv gain): band-limited signal power during
    // the sweep vs the noise floor in the decayed tail. This is what REW reports.
    const std::size_t L = rec.size();
    if (L >= 4096)
    {
        auto slice = [&] (double a, double b)
        {
            std::vector<double> w;
            const std::size_t i0 = (std::size_t) (a * (double) L), i1 = (std::size_t) (b * (double) L);
            for (std::size_t i = i0; i < i1 && i < L; ++i) w.push_back (rec[i]);
            return w;
        };
        const auto sig = slice (0.15, 0.55);   // solidly inside the sweep
        const auto nse = slice (0.90, 1.00);   // solidly inside the decayed tail
        const double sp = bandEnergy (sig, sw.spec.sampleRate, cfg.bandLoHz, cfg.bandHiHz)
                        / (double) std::max<std::size_t> (1, sig.size());
        const double np = bandEnergy (nse, sw.spec.sampleRate, cfg.bandLoHz, cfg.bandHiHz)
                        / (double) std::max<std::size_t> (1, nse.size()) + 1e-30;
        g.snrDb    = 10.0 * std::log10 ((sp + 1e-30) / np);
        g.snrValid = true;                                 // measured (else snrDb stays 0 + snrValid=false)
    }

    g.ok = ! g.clipped && g.sweepPresent;
    if (g.clipped)             g.reason = GateReject::Clipped;
    else if (! g.sweepPresent) g.reason = GateReject::SweepNotDetected;
    return g;
}

} // namespace felitronics::measurement

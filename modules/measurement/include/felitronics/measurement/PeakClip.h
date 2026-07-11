// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

//==============================================================================
// felitronics::measurement — standalone peak / clipping scan of a recording. The clipping verdict is
// a flat-top RUN at full scale (real saturation), NOT mere peak proximity, so a hot-but-clean
// −0.1 dBFS take is not false-rejected. Extracted from CaptureGate's gateRecording so consumers that
// own no sweep (e.g. a NAM reamp-capture take) can gate quality without dragging in Convolve/Sweep.
// OFFLINE, MESSAGE-THREAD-ONLY. Zero dependencies beyond <cmath>.
//==============================================================================

#include <cmath>
#include <span>

namespace felitronics::measurement
{

struct PeakClipConfig
{
    double clipLevel      = 0.999;   // full-scale threshold for the flat-top run
    int    clipRunSamples = 3;       // consecutive full-scale samples that count as real saturation
};

struct PeakClipReport
{
    bool   nonFinite = false;        // a NaN/Inf was seen — peak/clip fields are then meaningless
    bool   clipped   = false;
    double peakDbfs  = -120.0;
    int    clipRun   = 0;            // longest flat-top run at full scale
};

// One pass over the recording: non-finite detection + peak + longest flat-top run. An empty
// recording reports the quiet defaults (peak −120 dBFS, not clipped) — emptiness is the caller's
// rejection to make (CaptureGate keeps its own EmptyRecording reason).
inline PeakClipReport scanPeakClip (std::span<const double> rec, const PeakClipConfig& cfg = {})
{
    PeakClipReport r;
    double pk = 0.0;
    int run = 0;
    for (double v : rec)
    {
        if (! std::isfinite (v)) { r.nonFinite = true; continue; }   // fabs(NaN) comparisons are all-false anyway; flag honestly
        const double a = std::fabs (v);
        pk = std::max (pk, a);
        if (a >= cfg.clipLevel) { if (++run > r.clipRun) r.clipRun = run; }
        else run = 0;
    }
    r.peakDbfs = (pk > 0.0) ? 20.0 * std::log10 (pk) : -120.0;
    r.clipped  = (r.clipRun >= cfg.clipRunSamples);
    return r;
}

} // namespace felitronics::measurement

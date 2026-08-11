// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <algorithm>
#include <cmath>

//==============================================================================
// Incubated inside TabbyEQ as `eqview`, moved home once OrbitAmp became the second consumer. The
// discipline that made the move a file move still holds: data in, values out — no JUCE, no plugin
// includes, no component state.
//==============================================================================
// PlotMap — the log-frequency / dB coordinate map of an EQ plot, and the SINGLE home of the
// freq↔px / dB↔px math. The owning view builds one from its current geometry whenever it needs to
// convert; painting, node placement and hit-testing all share the same instance, so they can never
// disagree about where a frequency lives on screen.
//
// Two vertical scales share the surface on purpose:
//  • the CURVE scale (dbToY/yToDb) spans ±dbRange over [0, plotBottom] — plotBottom stops short of
//    the full height when the owner reserves a bottom strip, and values beyond ±dbRange keep
//    mapping linearly PAST those bounds (a deep notch deliberately overshoots the nominal scale and
//    dives for the window bottom — do not clamp);
//  • the SPECTRUM scale (specDbToY) spans [specTop, specBottom] dBFS over the FULL height, clamped.

namespace felitronics::analysis
{

struct PlotMap
{
    float  width      = 0.0f;                    // full plot width, px
    float  height     = 0.0f;                    // full plot height, px (spectrum scale spans this)
    float  plotBottom = 0.0f;                    // bottom of the curve/node area, px (≤ height)
    double freqMin    = 20.0,  freqMax = 28000.0; // log X axis range, Hz
    double dbRange    = 12.0;                    // visible ±dB window of the curve Y scale
    double specTop    = 6.0,   specBottom = -90.0; // spectrum dBFS scale

    float freqToX (double f) const noexcept
    {
        const double lo = std::log10 (freqMin), hi = std::log10 (freqMax);
        return (float) ((std::log10 (std::clamp (f, freqMin, freqMax)) - lo) / (hi - lo) * (double) width);
    }

    double xToFreq (float x) const noexcept
    {
        const double lo = std::log10 (freqMin), hi = std::log10 (freqMax);
        const double t  = std::clamp ((double) x / std::max (1.0, (double) width), 0.0, 1.0);
        return std::pow (10.0, lo + t * (hi - lo));
    }

    float dbToY (double db) const noexcept
    {
        return (float) ((0.5 - db / (2.0 * dbRange)) * (double) plotBottom);
    }

    double yToDb (float y) const noexcept
    {
        return (0.5 - (double) y / std::max (1.0, (double) plotBottom)) * (2.0 * dbRange);
    }

    float specDbToY (double db) const noexcept
    {
        const double t = std::clamp ((specTop - db) / (specTop - specBottom), 0.0, 1.0);
        return (float) (t * (double) height);
    }
};

} // namespace felitronics::analysis

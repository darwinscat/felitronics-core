// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.
#pragma once
// felitronics::rigplayer — ONE PARAMETRIC BAND AS A BIQUAD, the formulas NAMZ-FORMAT.md spells out for
// `sections`. This is the half of core/ToneSections.h that a PLAYER needs: a shape, a frequency, a Q
// and a gain go in, five coefficients come out. The other half — finding those numbers in a measured
// ladder — is the bench's and stays there. Split so the pack player can leave with the code that plays
// a section and without the code that fits one.
//
// RBJ cookbook, the same formulas felitronics::eq speaks, so a curve designed on the capture side and a
// biquad run on the player side describe one filter. The coefficient struct IS felitronics::eq's, so
// the player runs these through its Biquad without copying five numbers into a twin.

#include <felitronics/eq/MatchedBiquad.h>

#include <algorithm>
#include <cmath>

namespace felitronics::rigplayer {

enum class SectionKind { LowShelf, Peak, HighShelf, Tilt };

inline const char* nameOf(SectionKind k) {
    return k == SectionKind::LowShelf ? "low shelf" : k == SectionKind::Peak ? "bell"
         : k == SectionKind::Tilt     ? "tilt"      : "high shelf";
}

using SectionBiquad = felitronics::eq::BiquadCoeffs;

inline SectionBiquad designSection(SectionKind kind, double hz, double gainDb, double q,
                                   double sampleRate, double pivot = 0.5) {
    const double fs = sampleRate > 0.0 ? sampleRate : 48000.0;
    hz = std::clamp(hz, 10.0, fs * 0.49);
    q  = std::clamp(q, 0.05, 20.0);
    const double w0 = 2.0 * 3.14159265358979323846 * hz / fs;
    const double cw = std::cos(w0), sw = std::sin(w0);
    const double A = std::pow(10.0, gainDb / 40.0);
    const double alpha = sw / (2.0 * q);
    double b0, b1, b2, a0, a1, a2;
    // Peak, low shelf, and then BOTH the high shelf and the tilt, which are the same coefficients:
    // the tilt only differs by the broadband trim applied at the bottom of this function.
    if (kind == SectionKind::Peak) {
        b0 = 1 + alpha * A; b1 = -2 * cw; b2 = 1 - alpha * A;
        a0 = 1 + alpha / A; a1 = -2 * cw; a2 = 1 - alpha / A;
    } else if (kind == SectionKind::LowShelf) {
        const double s = 2.0 * std::sqrt(A) * alpha;
        b0 =     A * ((A + 1) - (A - 1) * cw + s);
        b1 = 2 * A * ((A - 1) - (A + 1) * cw);
        b2 =     A * ((A + 1) - (A - 1) * cw - s);
        a0 =         ((A + 1) + (A - 1) * cw + s);
        a1 = -2 *    ((A - 1) + (A + 1) * cw);
        a2 =         ((A + 1) + (A - 1) * cw - s);
    } else {
        const double s = 2.0 * std::sqrt(A) * alpha;
        b0 =      A * ((A + 1) + (A - 1) * cw + s);
        b1 = -2 * A * ((A - 1) + (A + 1) * cw);
        b2 =      A * ((A + 1) + (A - 1) * cw - s);
        a0 =          ((A + 1) - (A - 1) * cw + s);
        a1 =  2 *     ((A - 1) - (A + 1) * cw);
        a2 =          ((A + 1) - (A - 1) * cw - s);
    }
    SectionBiquad q2;
    q2.b0 = b0 / a0; q2.b1 = b1 / a0; q2.b2 = b2 / a0; q2.a1 = a1 / a0; q2.a2 = a2 / a0;
    // A TILT IS THAT HIGH SHELF WITH THE FLOOR TAKEN OUT FROM UNDER IT. The numerator is one scalar,
    // so a broadband trim costs no second filter and no second state: the shelf still lifts the top by
    // g, and the trim drops everything by pivot*g, which leaves the two arms straddling the hinge.
    // This is why `tilt` is a KIND and not a field beside the pack's gain — the level that makes a
    // tone control a tone control rather than a volume control is inside the filter, where it travels.
    if (kind == SectionKind::Tilt) {
        const double trim = std::pow(10.0, -std::clamp(pivot, 0.0, 1.0) * gainDb / 20.0);
        q2.b0 *= trim; q2.b1 *= trim; q2.b2 *= trim;
    }
    return q2;
}

// The band's response at one frequency, dB. What the fit is scored against and what the player draws,
// so a fit measured one way and played another cannot happen: there is one magnitude.
inline double sectionMagnitudeDb(const SectionBiquad& q, double hz, double sampleRate) {
    const double fs = sampleRate > 0.0 ? sampleRate : 48000.0;
    const double w = 2.0 * 3.14159265358979323846 * hz / fs;
    return 20.0 * std::log10(std::max(1e-12, q.magnitude(w)));
}

} // namespace felitronics::rigplayer
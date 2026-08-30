// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.
#pragma once
// felitronics::rigplayer — WHICH TWO MODELS SOUND AT THIS KNOB POSITION, AND IN WHAT PROPORTION.
//
// A gain knob is captured at a handful of angles and the hand stops anywhere. The rule here is the
// simplest one that can work, and it is deliberately not a measurement:
//
//     the two neighbouring captures play together, mixed by ANGLE.
//
// No decibels are fitted, nothing has to be trained before the dial works, and no position can fall
// into a hole. It rests on a fact about how this library is shot rather than on a theory: the angles
// are dense enough that stepping between them without any crossfade already sounds acceptable, so
// anything smooth between two of them is at least as good. A fitted drive curve would buy accuracy
// that nobody listens for, at hours of inference per device.
//
// TWO EXCEPTIONS, both asked for by the bench:
//
// OUTSIDE THE CAPTURED RANGE, THE NEAREST CAPTURE IS FED DIFFERENTLY. The dial turns its whole
// travel whatever was shot on it: below the lowest capture that one plays SOFTER, above the highest
// it plays HARDER, both at the rate below. Downwards is the safe direction — less signal into a model
// keeps it inside what it was trained on — and upwards is admitted extrapolation, which exists
// because a dial that physically reaches 17h with its last capture at 15h has a dead top, the exact
// complaint this mechanism was built to cure. The rate is a starting guess for the ear, not a reading.

#include <algorithm>
#include <cstddef>
#include <span>
#include <vector>

namespace felitronics::rigplayer {

// How much harder the top capture is driven per degree past it. A guess, and printed in the player's
// dump so it can be argued with: measured 30-degree spans in this library run 0.5 to 6.5 dB, and the
// ones at the TOP of a dial are the small ones — a saturating stage stops answering. So the rate is
// taken from that end, not from the average, and a full two hours past the last capture adds ~3 dB.
inline constexpr double kTopExtendDbPerDeg = 0.05;

struct BlendPick {
    int    lo = -1, hi = -1;   // indices into the knot list; equal = one model alone; -1 = nothing to play
    double mixHi = 0.0;        // weight of `hi`; `lo` gets 1 - mixHi
    double extendDb = 0.0;     // extra drive into the model, above the top knot only
};

// WHERE THE HANDOVER SITS AND HOW WIDE IT IS — the two numbers the bench turns.
//
// The original law faded across the WHOLE gap, so the handover was always at the midpoint and always
// as wide as the gap allowed. Both are voicing decisions, not facts: on a dial whose captures are 30
// degrees apart, "the next model starts a third of the way up, and takes ten degrees to arrive" is a
// different instrument from "half way, and instantly" — and only an ear can say which is the pedal.
//
//   point  0.01 … 0.99  where the 50/50 lands, as a fraction of the gap between the two captures.
//   width  0 … 1        how much of the room on EACH SIDE of that point the fade uses. Each side is
//                       measured against its own distance to its own capture, so the two halves are
//                       free to differ: at a point of 10% and a width of 90% the fade starts just
//                       above the lower capture and finishes just short of the upper one, arriving
//                       early and taking its time. Measuring both halves against the NEARER capture —
//                       which is what this did first — made 90% mean a band of ±9% and left nine
//                       tenths of the span unused (Oleh, off the picture). Zero is a hard switch AT
//                       the point; no width can reach past a capture into the next gap.
//
// The defaults reproduce the original law exactly (midpoint, full span), so a caller that says nothing
// gets what it always got.
struct BlendShape {
    double point = 0.5;
    double width = 1.0;
};

// `knots` are the angles that HAVE a model, ascending. `deg` is where the knob stands.
inline BlendPick pickBlend(std::span<const double> knots, double deg,
                           double topExtendDbPerDeg = kTopExtendDbPerDeg,
                           BlendShape shape = {}) {
    BlendPick p;
    if (knots.empty()) return p;
    const int last = (int) knots.size() - 1;
    if (deg <= knots.front()) {
        p.lo = p.hi = 0;
        p.extendDb = (deg - knots.front()) * topExtendDbPerDeg;    // negative: fed softer
        return p;
    }
    if (deg >= knots.back()) {
        p.lo = p.hi = last;
        p.extendDb = (deg - knots.back()) * topExtendDbPerDeg;
        return p;
    }
    for (int i = 0; i < last; ++i)
        if (deg >= knots[(std::size_t) i] && deg < knots[(std::size_t) i + 1]) {
            const double span = knots[(std::size_t) i + 1] - knots[(std::size_t) i];
            p.lo = i; p.hi = i + 1;
            // A span of zero would be two captures at one angle — take the lower and do not divide.
            const double t = span > 0.0 ? (deg - knots[(std::size_t) i]) / span : 0.0;
            const double pt = std::clamp(shape.point, 0.01, 0.99);
            const double k  = std::clamp(shape.width, 0.0, 1.0);
            const double wl = k * pt, wr = k * (1.0 - pt);       // each side stretches into its own room
            // Zero width is a switch, not a division by nothing: below the point the lower capture, at
            // or above it the upper one. Otherwise two ramps meeting at exactly 0.5 on the point — the
            // kink IS the decision, and either half may be the long one.
            p.mixHi = (wl <= 0.0 && wr <= 0.0) ? (t < pt ? 0.0 : 1.0)
                    : t < pt ? (wl <= 0.0 ? 0.0 : std::clamp(0.5 * (t - (pt - wl)) / wl, 0.0, 0.5))
                             : (wr <= 0.0 ? 1.0 : std::clamp(0.5 + 0.5 * (t - pt) / wr, 0.5, 1.0));
            return p;
        }
    p.lo = p.hi = last;
    return p;
}

// The same pick with the crossfade taken away: whichever neighbour is nearer plays alone. This is
// what the dial did before, and it stays reachable so the two can be compared by ear rather than
// argued about — the crossfade has to earn its place.
inline BlendPick stepped(BlendPick p) {
    if (p.lo < 0) return p;
    if (p.mixHi >= 0.5) p.lo = p.hi; else p.hi = p.lo;
    p.mixHi = 0.0;
    return p;
}

} // namespace felitronics::rigplayer
// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.
#pragma once
// felitronics::rigplayer — A BLEND KNOB OF THE PACK, AT ONE POSITION. A dry/wet mix is neither a
// captured axis nor a filter: the dry signal IS the DI the model is being fed, so the pack describes the
// mix and the player performs it (namz_rig.h, NAMZ-FORMAT.md):
//
//     out = polarity · dry_gain(pos) · (DI * dry) + wet_gain(pos) · model(DI)
//
// `dry` is the dry path's own response, shipped as SHAPE (`dry_db`, mean removed) and LEVEL
// (`dry_level_db`) — both, because neither can be recovered from the other and nobody can re-measure a
// pedal they do not hold. The gains per position are shipped explicitly; `law` is provenance only.
// Pure, std-only.

#include <felitronics/rigplayer/ToneKnobs.h>     // the trust and floor rules — one set for every curve of the pack

#include <namz_rig.h>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace felitronics::rigplayer {

using namz::rig::Blend;

struct BlendGains {
    double dry = 0.0;       // linear, WITH the dry path's level and the box's polarity folded in
    double wet = 1.0;       // linear
};

inline double linOf(double db) { return db <= -120.0 ? 0.0 : std::pow(10.0, db / 20.0); }

inline bool blendNorm(const Blend& b, const std::string& value, double& norm) {
    if (b.sweep <= 0) return false;
    const int deg = namz::rig::detail::degrees(value);
    if (deg < 0) return false;
    norm = std::clamp((double) deg / (double) b.sweep, 0.0, 1.0);
    return true;
}

inline std::string blendStart(const Blend& b) { return b.defaultValue.empty() ? b.reference : b.defaultValue; }

inline std::vector<double> gridOf(const Blend& b) {
    return b.grid.points > 0 ? felitronics::lineareq::logFreqGrid(b.grid.fLo, b.grid.fHi, b.grid.points)
                             : std::vector<double> {};
}

// The pair of gains at a dial's rotation: the two shipped positions either side, interpolated in LINEAR
// amplitude — a pot is a divider, and a straight line in decibels from -120 to 0 would put the middle of
// the travel sixty decibels down. Outside the shipped range the nearest position holds. A knob that
// ships no positions is full wet: exactly what its weights encode.
inline BlendGains blendGainsAtNorm(const Blend& b, double norm) {
    BlendGains g;
    if (b.positions.empty()) return g;
    const auto dryOf = [&](const namz::rig::BlendPosition& p) { return linOf(p.dryDb) * linOf(b.dryLevelDb); };
    const auto wetOf = [](const namz::rig::BlendPosition& p) { return linOf(p.wetDb); };
    const auto& ps = b.positions;                                      // ascending by norm, by the producer
    if (norm <= ps.front().norm || ps.size() == 1) { g.dry = dryOf(ps.front()); g.wet = wetOf(ps.front()); }
    else if (norm >= ps.back().norm)                { g.dry = dryOf(ps.back());  g.wet = wetOf(ps.back()); }
    else {
        std::size_t hi = 1;
        while (hi + 1 < ps.size() && ps[hi].norm < norm) ++hi;
        const auto& a = ps[hi - 1];
        const auto& c = ps[hi];
        const double span = c.norm - a.norm;
        const double t = span > 1e-12 ? (norm - a.norm) / span : 0.0;
        g.dry = dryOf(a) + t * (dryOf(c) - dryOf(a));
        g.wet = wetOf(a) + t * (wetOf(c) - wetOf(a));
    }
    if (b.polarity < 0) g.dry = -g.dry;
    return g;
}

// A switch: that position's own pair; a value the pack never measured plays full wet.
inline BlendGains blendGainsAtValue(const Blend& b, const std::string& value) {
    BlendGains g;
    for (const auto& p : b.positions)
        if (p.value == value) {
            g.dry = linOf(p.dryDb) * linOf(b.dryLevelDb) * (b.polarity < 0 ? -1.0 : 1.0);
            g.wet = linOf(p.wetDb);
            return g;
        }
    return g;
}

inline BlendGains blendGainsAt(const Blend& b, const std::string& value) {
    double norm = 0.0;
    return blendNorm(b, value, norm) ? blendGainsAtNorm(b, norm) : blendGainsAtValue(b, value);
}

// The dry path's SHAPE on the grid, with the pack's trusted band and this player's floor applied — the
// level rides in the gains above. Empty = a wire (no filter to convolve), which a flat shape also is.
inline std::vector<double> dryCurve(const Blend& b, const std::vector<double>& grid) {
    if (b.dryDb.empty() || b.dryDb.size() != grid.size()) return {};
    return heldToFloor(heldToTrust(b.dryDb, b.trusted), grid);
}

} // namespace felitronics::rigplayer
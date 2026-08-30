// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.
#pragma once
// felitronics::rigplayer — WHICH FILES OF THE PACK SOUND AT THIS PANEL, and in what proportion.
//
// Two policies meet here and neither is written here:
//
//   the COMBINATION — which switch positions, which captured axes — is namz::rig's: a turned control
//   is law, everything else stays where the hand left it, and the closest captured combination wins
//   (namz::rig::resolve). The player does not get a second opinion about that.
//
//   the PAIR along the gain dial is ModelBlend.h's: the two captures either side of the angle play
//   together, mixed by angle, and past the ends the nearest one is fed differently (pickBlend).
//
// What this file adds is only the join: the files of the pack that lie along the dial at the panel's
// OTHER settings are the knots, sorted by their angle, and the pair is picked among them. A file that
// several settings point at — the pack's spelling of a link — is a knot at each of those settings, with
// its own `input_db`; it plays its neighbour's weights fed softer, and the dial fades into it instead
// of falling into a hole.
//
// std-only and pure, so every rule about who sounds is provable without a sound card.

#include <felitronics/rigplayer/ModelBlend.h>

#include <namz_rig.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

namespace felitronics::rigplayer {

using namz::rig::Control;
using namz::rig::Device;
using namz::rig::FileEntry;
using namz::rig::Settings;

// "150" → 150; anything that is not a bare integer of degrees → -1. namz::rig's own reading of a dial
// position, so the player and the format agree about what a position IS.
inline int degreesOf(const std::string& v) { return namz::rig::detail::degrees(v); }

// The dial the crossfade runs along: the first Gain-role control that turns and has at least one
// position in degrees. nullptr = this device has no such dial, and every panel is one file alone.
inline const Control* crossfadeDial(const Device& d) {
    for (const auto& c : d.controls) {
        if (c.role != namz::rig::Role::Gain || c.sweep <= 0) continue;
        for (const auto& v : c.values)
            if (degreesOf(v) >= 0) return &c;
    }
    return nullptr;
}

// A dial that is NOT the crossfade axis selects like a switch: the captured position nearest the
// angle. Empty when the control has no position in degrees at all.
inline std::string nearestValue(const Control& c, double deg) {
    std::string best;
    double bestDist = 1.0e300;
    for (const auto& v : c.values) {
        const int d = degreesOf(v);
        if (d < 0) continue;
        const double dist = std::abs((double) d - deg);
        if (dist < bestDist) { bestDist = dist; best = v; }
    }
    return best;
}

struct Knot {
    int    file = -1;      // index into Device::files
    double deg  = 0.0;     // where on the dial it was captured
};

// The files along `dial` whose every OTHER setting equals `axes`, ascending by angle. Two files at one
// angle cannot both be knots — the later one in the index is dropped, so the list stays a function.
inline std::vector<Knot> knotsAlong(const Device& d, const Settings& axes, const std::string& dial) {
    std::vector<Knot> out;
    for (std::size_t i = 0; i < d.files.size(); ++i) {
        const auto& f = d.files[i];
        const auto dv = f.settings.find(dial);
        if (dv == f.settings.end()) continue;
        const int deg = degreesOf(dv->second);
        if (deg < 0) continue;
        bool same = true;
        for (const auto& c : d.controls) {
            if (c.name == dial) continue;
            const auto fv = f.settings.find(c.name);
            const auto av = axes.find(c.name);
            const std::string fs = fv != f.settings.end() ? fv->second : std::string();
            const std::string as = av != axes.end() ? av->second : std::string();
            if (fs != as) { same = false; break; }
        }
        if (! same) continue;
        if (std::any_of(out.begin(), out.end(), [deg](const Knot& k) { return (int) std::lround(k.deg) == deg; }))
            continue;
        out.push_back({ (int) i, (double) deg });
    }
    std::sort(out.begin(), out.end(), [](const Knot& a, const Knot& b) { return a.deg < b.deg; });
    return out;
}

struct Selection {
    std::vector<Knot> knots;          // along the dial at this panel, ascending
    BlendPick pick;                   // lo/hi INTO `knots`
    int    fileA = -1, fileB = -1;    // Device::files indices; equal = one model alone; -1 = nothing
    double mixB = 0.0;                // the weight of fileB; fileA gets 1 - mixB
    double extendDb = 0.0;            // the chain's trim past the ends of the captured range
    double deg = 0.0;                 // where the dial stands
    bool   pair() const { return fileA >= 0 && fileB >= 0 && fileA != fileB; }
};

// The whole decision for one panel: `axes` are the switch and non-crossfade-dial settings (the
// crossfade dial's own entry is ignored), `deg` is the crossfade dial's angle. With no dial, or no knot
// at this panel, the file whose settings equal `axes` plays alone.
inline Selection select(const Device& d, const Settings& axes, const std::string& dial, double deg,
                        BlendShape shape = {}, double topExtendDbPerDeg = kTopExtendDbPerDeg) {
    Selection s;
    s.deg = deg;
    if (! dial.empty()) s.knots = knotsAlong(d, axes, dial);
    if (s.knots.empty()) {
        if (const auto* f = d.find(axes)) s.fileA = s.fileB = (int) (f - d.files.data());
        return s;
    }
    std::vector<double> degs;
    degs.reserve(s.knots.size());
    for (const auto& k : s.knots) degs.push_back(k.deg);
    s.pick = pickBlend(degs, deg, topExtendDbPerDeg, shape);
    if (s.pick.lo < 0) return s;
    s.fileA    = s.knots[(std::size_t) s.pick.lo].file;
    s.fileB    = s.knots[(std::size_t) s.pick.hi].file;
    s.mixB     = s.pick.mixHi;
    s.extendDb = s.pick.extendDb;
    return s;
}

// WHICH SLOT EACH CAPTURE LIVES IN — by the PARITY of its position on the dial, never by which one is
// nearer. Nearness swaps the roles at the halfway mark, which is the one place both models are at
// full-ish weight; the freshly-loaded one then speaks from an empty memory while it is half the sound,
// and that is the crackle. By parity, a capture keeps its slot for as long as it is audible at all:
// crossing a capture replaces only the FAR one, whose weight is exactly zero at that moment.
//
// AT THE ENDS OF THE DIAL BOTH SLOTS ASK FOR THE SAME CAPTURE, rather than one of them asking for
// nothing. The law can only repair a slot it has been given a model for: told "nothing", it marks that
// slot wrong for ever, and in some orders of events the weight settles on it and stays.
struct SlotPlan {
    int    file[2] { -1, -1 };        // Device::files index per slot; -1 = nothing to play
    double targetB = 0.0;             // the weight slot 1 should end up with
    double deg[2] { 0.0, 0.0 };       // the angle each slot's capture was shot at
    double inputDb[2] { 0.0, 0.0 };   // …and how much softer it is fed (the file's own `input_db`)
};

inline SlotPlan slotPlan(const Selection& s, const Device& d) {
    SlotPlan p;
    if (s.fileA < 0) return p;
    const bool loEven = s.pick.lo < 0 || (s.pick.lo % 2) == 0;
    p.file[0] = loEven ? s.fileA : s.fileB;
    p.file[1] = loEven ? s.fileB : s.fileA;
    p.targetB = loEven ? s.mixB : 1.0 - s.mixB;
    if (s.pick.lo >= 0) {
        const double degLo = s.knots[(std::size_t) s.pick.lo].deg;
        const double degHi = s.knots[(std::size_t) s.pick.hi].deg;
        p.deg[0] = loEven ? degLo : degHi;
        p.deg[1] = loEven ? degHi : degLo;
    }
    for (int i = 0; i < 2; ++i)
        if (p.file[i] >= 0 && (std::size_t) p.file[i] < d.files.size())
            p.inputDb[i] = d.files[(std::size_t) p.file[i]].inputDb;
    return p;
}

} // namespace felitronics::rigplayer
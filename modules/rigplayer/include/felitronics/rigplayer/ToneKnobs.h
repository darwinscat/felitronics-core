// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.
#pragma once
// felitronics::rigplayer — A TONE KNOB OF THE PACK, AT ONE POSITION. The pack ships a linear knob
// in ONE of two forms (namz_rig.h, NAMZ-FORMAT.md): `positions`, the measured ladder — a dB table per
// swept position on a log grid — or `sections`, the same knob as parametric bands whose gain travels
// with the dial. Both are anchored at `reference`, the position every model was captured at, where the
// knob is flat by construction. This file turns either form into what the player runs — a curve on the
// grid, to become a FIR, or a set of biquads — and nothing else. Pure, std-only, testable.
//
// The rules applied are the format's, stated where each is applied: which positions to interpolate
// between and against what, where a curve is held rather than followed, how t runs from the reference to
// each stop and what follows it.

#include <felitronics/rigplayer/SectionBiquad.h>

#include <felitronics/lineareq/MagnitudeCurve.h>
#include <namz_rig.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

namespace felitronics::rigplayer {

using namz::rig::Tone;

// THE BOTTOM OF THE SWEEP IS NOT A MEASUREMENT, so a curve is held at 50 Hz on the way down rather than
// followed into it. Proof from the ladders on disk: every position is shot at three drive levels, and a
// tone stack cannot change with drive — so the disagreement between those shots IS the noise. Worst case
// over fourteen positions of one device: 0.2 dB above 1 kHz, 1.9 at 50-60 Hz, 7.0 at 30-40 and 15.1 dB
// at 20-30 — while the knobs themselves move about 10 dB down there. Left alone a +13 dB spike at 22 Hz
// went into the audio ahead of the cabinet: cone travel and headroom spent on a number nobody measured.
// The top needs no such guard, so the high edge is the grid's own end. The same floor the capture app's
// bench holds, carried here so the pack plays the way it was auditioned.
inline constexpr double kCurveFloorHz = 50.0;

// A dial's position → its rotation 0..1. False for a switch (no sweep) or a value that is not degrees.
inline bool toneNorm(const Tone& t, const std::string& value, double& norm) {
    if (t.sweep <= 0) return false;
    const int deg = namz::rig::detail::degrees(value);
    if (deg < 0) return false;
    norm = std::clamp((double) deg / (double) t.sweep, 0.0, 1.0);
    return true;
}

// Where a player starts the knob: `default` when stated, else `reference`. Never invented.
inline std::string toneStart(const Tone& t) { return t.defaultValue.empty() ? t.reference : t.defaultValue; }

inline std::vector<double> gridOf(const Tone& t) {
    return t.grid.points > 0 ? felitronics::lineareq::logFreqGrid(t.grid.fLo, t.grid.fHi, t.grid.points)
                             : std::vector<double> {};
}

// ---------------------------------------------------------------------------------- the curve form

// The pack's `trusted` band, applied on the grid — the INDICES are the authority, the Hz are the
// same band for humans. Outside the band the curve is held at the nearest edge; a band that was tested
// and failed everywhere (two or more levels, empty band) applies nothing; a band never tested trusts
// the whole grid.
inline std::vector<double> heldToTrust(std::vector<double> db, const namz::rig::ToneTrust& tr) {
    if (db.empty()) return db;
    const bool tested = tr.levels >= 2;
    const bool stated = tr.loHz > 0.0 || tr.hiHz > 0.0 || tr.loIndex != 0 || tr.hiIndex != 0;
    if (! tested || ! stated) return db;
    const int n = (int) db.size();
    const bool failedEverywhere = tr.loIndex > tr.hiIndex || (tr.hiHz > 0.0 && tr.hiHz <= tr.loHz);
    if (failedEverywhere) return std::vector<double>(db.size(), 0.0);
    const int lo = std::clamp(tr.loIndex, 0, n - 1), hi = std::clamp(tr.hiIndex, 0, n - 1);
    for (int i = 0; i < lo; ++i)      db[(std::size_t) i] = db[(std::size_t) lo];
    for (int i = hi + 1; i < n; ++i)  db[(std::size_t) i] = db[(std::size_t) hi];
    return db;
}

// …and the floor above, which is this player's, not the pack's.
inline std::vector<double> heldToFloor(const std::vector<double>& db, const std::vector<double>& grid) {
    return grid.empty() ? db : felitronics::lineareq::heldOutsideBand(db, grid, kCurveFloorHz, grid.back());
}

// The curve a DIAL has at `norm` (0..1 of its rotation): the two measured positions that bracket it,
// interpolated against their own `norm`s, never against their index — the producer states where each
// was, and the producer wins. Outside the measured range the nearest curve is held. Empty = flat.
inline std::vector<double> curveAtNorm(const Tone& t, const std::vector<double>& grid, double norm) {
    if (t.positions.empty() || grid.empty()) return {};
    std::vector<std::vector<double>> curves;
    std::vector<double> norms;
    for (const auto& p : t.positions) {
        if (p.db.size() != grid.size()) continue;              // the loader refuses these; belt and braces
        curves.push_back(heldToTrust(p.db, t.trusted));
        norms.push_back(p.norm);
    }
    if (curves.empty()) return {};
    return heldToFloor(felitronics::lineareq::curveAtPosition(curves, norms, norm), grid);
}

// The curve a SWITCH has at `value`: that position's own, and nothing between. Empty = flat, which is
// also the answer for a value the pack never swept — the reference tone is what every knob promises.
inline std::vector<double> curveAtValue(const Tone& t, const std::vector<double>& grid,
                                        const std::string& value) {
    for (const auto& p : t.positions)
        if (p.value == value && p.db.size() == grid.size())
            return heldToFloor(heldToTrust(p.db, t.trusted), grid);
    return {};
}

// Whichever the knob is. Degrees for a dial, the value itself for a switch.
inline std::vector<double> curveAt(const Tone& t, const std::vector<double>& grid, const std::string& value) {
    double norm = 0.0;
    return toneNorm(t, value, norm) ? curveAtNorm(t, grid, norm) : curveAtValue(t, grid, value);
}

// ------------------------------------------------------------------------------- the sections form

// THE TRAVEL LAW. t is the knob's rotation normalised so that `reference` is 0, the minus stop (the
// dial's 0) is -1 and the plus stop (`sweep`) is +1 — two straight segments meeting at zero, not one
// line through the two stops. A function of POSITION and of nothing else: never recovered from a gain.
inline double travelOf(const Tone& t, double norm) {
    if (t.sweep <= 0) return 0.0;
    const int ref = namz::rig::detail::degrees(t.reference);
    if (ref < 0 || ref > t.sweep) return 0.0;                   // the loader drops such a knob; be safe
    const double normRef = (double) ref / (double) t.sweep;
    norm = std::clamp(norm, 0.0, 1.0);
    if (std::abs(norm - normRef) < 1e-12) return 0.0;
    if (norm < normRef) return normRef > 0.0 ? (norm - normRef) / normRef : 0.0;
    return normRef < 1.0 ? (norm - normRef) / (1.0 - normRef) : 0.0;
}

inline felitronics::rigplayer::SectionKind kindOf(namz::rig::SectionKind k) {
    switch (k) {
        case namz::rig::SectionKind::LowShelf:  return felitronics::rigplayer::SectionKind::LowShelf;
        case namz::rig::SectionKind::Bell:      return felitronics::rigplayer::SectionKind::Peak;
        case namz::rig::SectionKind::HighShelf: return felitronics::rigplayer::SectionKind::HighShelf;
        case namz::rig::SectionKind::Tilt:      return felitronics::rigplayer::SectionKind::Tilt;
    }
    return felitronics::rigplayer::SectionKind::Peak;
}

// One band at travel `t`: the gain runs linearly in dB from zero at the reference to `range_db` at the
// stop on t's side; the frequency heads for `hz_at` on that side geometrically and the Q for `q_at`
// linearly, by |t|. A side that is not stated (0) does not move.
inline felitronics::rigplayer::SectionBiquad sectionAt(const namz::rig::Section& s, double t, double sampleRate) {
    t = std::clamp(t, -1.0, 1.0);
    const double a = std::abs(t);
    const double gainDb = t < 0.0 ? -t * s.dbAtMin : t * s.dbAtMax;
    const double hzEnd  = t < 0.0 ? s.hzAtMin : s.hzAtMax;
    const double qEnd   = t < 0.0 ? s.qAtMin  : s.qAtMax;
    double hz = s.hz, q = s.q;
    if (hzEnd > 0.0 && s.hz > 0.0) hz = s.hz * std::pow(hzEnd / s.hz, a);
    if (qEnd > 0.0)                q  = s.q + (qEnd - s.q) * a;
    return felitronics::rigplayer::designSection(kindOf(s.kind), hz, gainDb, q, sampleRate, s.pivot);
}

// Every band of a `sections` knob at `value` (degrees). A knob whose bands cannot be built is empty —
// the reference position plays, which is what dropping a tone knob has always meant.
inline std::vector<felitronics::rigplayer::SectionBiquad> sectionsAt(const Tone& t, const std::string& value,
                                                        double sampleRate) {
    std::vector<felitronics::rigplayer::SectionBiquad> out;
    double norm = 0.0;
    if (t.sections.empty() || ! toneNorm(t, value, norm)) return out;
    const double tr = travelOf(t, norm);
    out.reserve(t.sections.size());
    for (const auto& s : t.sections) out.push_back(sectionAt(s, tr, sampleRate));
    return out;
}

} // namespace felitronics::rigplayer
// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.
#pragma once
// felitronics::lineareq — a MAGNITUDE TABLE turned into a filter you can actually run.
//
// The table is dB against a log-spaced frequency grid: what a measurement produces, and what a pack
// ships. That is deliberately not a filter — no bands, no coefficients, nothing tied to a sample rate
// — so building one is the reader's job, and this is the reference way of doing it.
//
// Minimum phase, via MixedPhaseFir at k = 1: the magnitude is reproduced exactly, there is no pre-ring
// and no bulk delay, so switching a control (or A/B-ing against the hardware through a reamp loop)
// does not shift the audio in time. A linear-phase design would need half its length as latency —
// about 10 ms for a filter long enough to shape 80 Hz.
//
// Why not a fitted 1-pole, which is what this used to be: measured on real hardware, a Big Muff tone
// control is a bass-cut/treble-boost blend, roughly 30 dB of tilt across the band, and the best 1-pole
// fit missed it by 15 to 65 dB. The curve is the only honest description, so the curve is what gets built.
//
// Message thread only: allocates, runs FFTs. The result is handed to a convolver.
//
// NO PRODUCT POLICY LIVES HERE. Trust thresholds, how wide a band must reproduce before it counts,
// what a pack is allowed to claim — those are opinions about a format, and a library holding someone's
// opinion is worse than a copy of the code. The trusted band arrives as two numbers because the CALLER
// knows what a trusted band is and this file does not.

#include <felitronics/lineareq/MixedPhaseFir.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace felitronics::lineareq {

// The exact frequencies a curve's dB values belong to, in order. No defaults on purpose: how many
// points and from where to where is the format's choice, and a library that guessed it would be
// stating someone else's opinion as a fact.
inline std::vector<double> logFreqGrid(double fLo, double fHi, int points) {
    const int pts = std::max(2, points);
    std::vector<double> f((std::size_t) pts);
    for (int i = 0; i < pts; ++i)
        f[(std::size_t) i] = fLo * std::pow(fHi / fLo, (double) i / (double) (pts - 1));
    return f;
}

// dB at `hz`, read off a curve sampled on `freqHz` (ascending, log-spaced). Interpolated in log
// frequency, because that is the axis the curve was measured on; held flat past either end, where the
// measurement says nothing and an extrapolation would be invention.
inline double curveDbAt(const std::vector<double>& db, const std::vector<double>& freqHz, double hz) {
    if (db.empty() || freqHz.empty()) return 0.0;
    if (hz <= freqHz.front()) return db.front();
    if (hz >= freqHz.back())  return db.back();
    const auto up = std::lower_bound(freqHz.begin(), freqHz.end(), hz);
    const std::size_t hi = (std::size_t) (up - freqHz.begin());
    const std::size_t lo = hi > 0 ? hi - 1 : 0;
    if (hi >= db.size() || lo >= db.size()) return db.back();
    const double fLo = freqHz[lo], fHi = freqHz[hi];
    if (! (fHi > fLo)) return db[lo];
    const double t = std::log(hz / fLo) / std::log(fHi / fLo);
    return db[lo] + t * (db[hi] - db[lo]);
}

// Hold the curve flat outside the band it was shown to be a filter in. Not a roll-off and not a cut:
// beyond the trusted edge nobody knows what the control does, and the nearest thing anybody does know
// is the value at that edge. A measurement that did not reproduce must not be applied — but it also
// must not be replaced with an invention.
inline std::vector<double> heldOutsideBand(const std::vector<double>& db, const std::vector<double>& freqHz,
                                           double loHz, double hiHz) {
    if (db.size() != freqHz.size()) return db;
    // An empty band is a STATEMENT, not a missing one: the curve reproduced at no frequency, so there
    // is nothing to apply anywhere. Returning it unmodified — which is what a lenient guard did — turns
    // the strongest possible warning into the widest possible permission. The caller distinguishes
    // "never tested" from "tested and failed" before it gets here; both ends of that decision are
    // outside this file, which knows only what two numbers mean.
    if (hiHz > 0.0 && loHz > 0.0 && hiHz <= loHz) return std::vector<double>(db.size(), 0.0);
    if (! (hiHz > loHz)) return db;
    std::vector<double> out = db;
    for (std::size_t i = 0; i < out.size(); ++i)                 // below the band: hold its low edge
        if (freqHz[i] >= loHz) { for (std::size_t k = 0; k < i; ++k) out[k] = db[i]; break; }
    for (std::size_t i = out.size(); i-- > 0;)                   // above it: hold its high edge
        if (freqHz[i] <= hiHz) { for (std::size_t k = i + 1; k < out.size(); ++k) out[k] = db[i]; break; }
    return out;
}

// The curve as a minimum-phase FIR at `sampleRate`. `taps` is what gets kept; `designSize` must be a
// power of two and comfortably longer (the cepstral design aliases in time otherwise — the module asks
// for 8x). An empty or flat curve returns {} — the caller reads that as "nothing to apply".
inline std::vector<float> magnitudeCurveToFir(const std::vector<double>& dbIn,
                                              const std::vector<double>& freqHz,
                                              double sampleRate, int taps = 1024,
                                              int designSize = 8192,
                                              double trustLoHz = 0.0, double trustHiHz = 0.0) {
    if (dbIn.size() < 2 || dbIn.size() != freqHz.size() || ! (sampleRate > 0.0)) return {};
    const auto db = heldOutsideBand(dbIn, freqHz, trustLoHz, trustHiHz);
    if (taps < 16 || designSize < 8 * taps) return {};

    // A curve that does nothing is not worth a convolution: a reference position is flat by
    // construction, and that is the common case while a control sits where it was captured.
    double worst = 0.0;
    for (const double v : db) worst = std::max(worst, std::abs(v));
    if (worst < 0.05) return {};

    MixedPhaseFir<> design;
    if (! design.prepare(designSize)) return {};

    const int H = designSize / 2;
    std::vector<float> mag((std::size_t) H + 1, 1.0f);
    for (int b = 0; b <= H; ++b) {
        const double hz = (double) b * sampleRate / (double) designSize;
        mag[(std::size_t) b] = (float) std::pow(10.0, curveDbAt(db, freqHz, hz) / 20.0);
    }

    const float* h = design.build(mag.data(), 1.0f);          // k = 1 → minimum phase
    if (h == nullptr) return {};

    // Keep the front, taper the tail. A minimum-phase impulse is front-loaded, so the cut costs almost
    // nothing; the taper is what stops the truncation itself from ringing.
    std::vector<float> fir((std::size_t) taps, 0.0f);
    const int fade = std::max(1, taps / 4);
    for (int i = 0; i < taps; ++i) {
        float w = 1.0f;
        if (i >= taps - fade) {
            const float t = (float) (i - (taps - fade)) / (float) fade;
            w = 0.5f * (1.0f + std::cos(3.14159265358979f * t));
        }
        fir[(std::size_t) i] = h[(std::size_t) i] * w;
    }
    return fir;
}

// The curve a control has at an arbitrary setting, interpolating between the two measured positions
// that bracket it — the whole point of measuring instead of capturing is that the knob does not have
// to land on one.
//
// Against `norms`, NOT against the array index. They agree only when the measured positions are evenly
// spaced, and there is no reason they should be: measuring 0/30/150/300 on a 300-degree dial (denser
// where the pot acts fastest) gives norms 0, 0.1, 0.5, 1.0, and at a knob halfway the two readings are
// four to six decibels apart on a Big Muff-class tilt. The producer states `norm`; the producer wins.
//
// `norms` must be ascending and the same length as `positions`. Outside the measured range the nearest
// curve is held — never extrapolated, which on a 30 dB tilt would invent decibels wholesale.
inline std::vector<double> curveAtPosition(const std::vector<std::vector<double>>& positions,
                                           const std::vector<double>& norms, double norm) {
    if (positions.empty()) return {};
    if (norms.size() != positions.size() || positions.size() == 1) return positions.front();
    if (norm <= norms.front()) return positions.front();
    if (norm >= norms.back())  return positions.back();
    std::size_t hi = 1;
    while (hi + 1 < norms.size() && norms[hi] < norm) ++hi;
    const std::size_t lo = hi - 1;
    const double span = norms[hi] - norms[lo];
    const double t = span > 1e-12 ? (norm - norms[lo]) / span : 0.0;
    const auto& a = positions[lo];
    const auto& b = positions[hi];
    std::vector<double> out(std::min(a.size(), b.size()));
    for (std::size_t i = 0; i < out.size(); ++i) out[i] = a[i] + t * (b[i] - a[i]);   // in dB, per point
    return out;
}

} // namespace felitronics::lineareq

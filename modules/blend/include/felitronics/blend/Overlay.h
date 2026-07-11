// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

//==============================================================================
// felitronics::blend — the mix-view overlay FACADE (reuse audit N6): everything a product needs
// to DRAW a multi-mic blend in one call — per-mic magnitude curves (post filters + gain), the
// blend curve (post Master chain, so its rolloff shows), and the interference column
// (where phase eats or reinforces vs the incoherent power sum). Extracted from OrbitCapture's
// Mixer-tab overlay so OrbitCab's Blend Lens (and any future mix view) computes the SAME picture
// instead of re-deriving the pipeline from the blend/analysis primitives.
//
// Display semantics baked in (each was a hard-won product decision — do not "fix" them here):
//   • The interference is computed against the PRE-Master-filter blend: a Master HPF/LPF rolloff
//     is not phase cancellation, and using the post-filter curve painted whole stopbands red.
//   • The interference is then FADED by the Master filter's |H| (a display weight, deliberately
//     not valid dB math): where the Master removed the signal there is nothing to tint.
//   • Muted / off-solo mics keep their curve slot (the UI may still list them) but are flagged
//     not audible and excluded from the incoherent reference — flags come from the same
//     solo/mute rules the audible blend uses (channelAudible).
//   • Curves are absolute dB (normalize is forced OFF): a shared display reference must be the
//     CALLER's, anchored once per take/session — per-render re-normalizing hides gain moves.
//   • Non-finite strip/Master numbers are HEALED to their defaults at the boundary (the family's
//     offline convention) — a NaN gain must not paint a NaN picture.
//   • The Master |H| fade curve is measured on an impulse of irs[0].size() samples — for very
//     short IRs its low-frequency read is coarse. Same as the product original; display-only.
//
// Message-thread only (allocates, double precision). All maths comes from felitronics::blend
// (processedMic / blendIrs / applyFilters) + felitronics::analysis::offline (logMagnitudeCurve /
// interferenceDb); this header adds NO new numerics — the NULL test pins it to the hand-composed
// pipeline at machine epsilon.
//==============================================================================

#include <felitronics/analysis/offline/SpectrumCurve.h>
#include <felitronics/blend/IrBlend.h>

#include <algorithm>
#include <cmath>
#include <span>
#include <vector>

namespace felitronics::blend {

struct OverlaySpec {
    analysis::offline::LogCurveSpec curve;   // grid/smoothing; `normalize` is ignored (forced OFF)
};

struct Overlay {
    std::vector<std::vector<double>> mic;    // [m] dB curve of mic m POST its filters + gain
    std::vector<bool> audible;               // [m] counts in the blend (solo/mute rules)
    std::vector<double> blend;               // dB curve of the audible blend POST Master chain
    std::vector<double> interference;        // dB vs incoherent sum, pre-Master basis, |H|-faded
};

namespace detail {
    inline double finiteOr(double v, double fallback) { return std::isfinite(v) ? v : fallback; }
    inline Filter healFilter(Filter f, const Filter& def) {
        f.hz = finiteOr(f.hz, def.hz);
        return f;
    }
    // Heal non-finite params to their defaults (offline convention) — a NaN gain paints no picture.
    inline StripParams heal(StripParams s) {
        const StripParams d;
        s.gainDb = finiteOr(s.gainDb, d.gainDb); s.phaseDeg = finiteOr(s.phaseDeg, d.phaseDeg);
        s.shiftMs = finiteOr(s.shiftMs, d.shiftMs);
        s.hpf = healFilter(s.hpf, d.hpf); s.lpf = healFilter(s.lpf, d.lpf);
        return s;
    }
    inline MasterParams heal(MasterParams m) {
        const MasterParams d;
        m.gainDb = finiteOr(m.gainDb, d.gainDb);
        m.hpf = healFilter(m.hpf, d.hpf); m.lpf = healFilter(m.lpf, d.lpf);
        return m;
    }
}

// One call per mixer change. Empty/degenerate inputs (no irs, ANY empty ir, strip-count mismatch,
// sr <= 0) → an empty Overlay; callers must treat an empty `blend` as "nothing to draw".
inline Overlay makeOverlay(std::span<const std::vector<float>> irs,
                           std::span<const StripParams> stripsIn, const MasterParams& masterIn,
                           double sr, const OverlaySpec& spec = {}) {
    namespace off = analysis::offline;
    Overlay out;
    if (irs.empty() || stripsIn.size() != irs.size() || sr <= 0.0 || !std::isfinite(sr)) return out;
    for (const auto& ir : irs) if (ir.empty()) return out;       // a partial overlay is worse than none
    std::vector<StripParams> strips(stripsIn.begin(), stripsIn.end());
    for (auto& s : strips) s = detail::heal(s);                  // offline convention: heal, don't paint NaN
    const MasterParams master = detail::heal(masterIn);
    auto cs = spec.curve;
    cs.normalize = false;                                        // absolute dB — see the header note

    out.mic.resize(irs.size());
    out.audible.resize(irs.size());
    for (std::size_t m = 0; m < irs.size(); ++m) {               // per-mic curve WITH its filters + gain
        std::vector<float> gi = processedMic(irs[m], strips[m], sr);
        const double g = std::pow(10.0, strips[m].gainDb / 20.0);
        for (auto& v : gi) v *= (float)g;                        // EXACT original expression (single-rounding
                                                                 // float×float — a double-rounded variant broke
                                                                 // bit-parity with the app; crew finding)
        out.mic[m]     = off::logMagnitudeCurve(std::span<const float>(gi), sr, cs);
        out.audible[m] = channelAudible(strips, m);
    }

    const std::vector<float> post = blendIrs(irs, strips, master, sr, true);    // WITH Master → rolloff shows
    out.blend = off::logMagnitudeCurve(std::span<const float>(post), sr, cs);
    if (out.blend.empty()) return out;

    const std::vector<float> pre = blendIrs(irs, strips, master, sr, false);    // interference basis (see note)
    const auto rawC = off::logMagnitudeCurve(std::span<const float>(pre), sr, cs);
    std::vector<double> coherent(out.blend.size());
    for (std::size_t p = 0; p < coherent.size(); ++p)
        coherent[p] = p < rawC.size() ? rawC[p] : out.blend[p];

    std::vector<off::MicCurveView> mv; mv.reserve(out.mic.size());
    for (std::size_t m = 0; m < out.mic.size(); ++m)
        mv.push_back({ std::span<const double>(out.mic[m]), out.audible[m] });
    out.interference = off::interferenceDb(coherent, mv);

    if (hpActiveSlope(master.hpf) || lpActiveSlope(master.lpf)) {               // |H| fade in the Master stopband
        std::vector<float> imp(irs[0].size(), 0.0f);
        if (!imp.empty()) imp[0] = 1.0f;
        applyFilters(imp, master.hpf, master.lpf, sr);
        const auto mc = off::logMagnitudeCurve(std::span<const float>(imp), sr, cs);   // dB, ~0 in the passband
        for (std::size_t p = 0; p < out.interference.size() && p < mc.size(); ++p)
            out.interference[p] *= std::pow(10.0, std::min(0.0, mc[p]) / 20.0);
    }
    return out;
}

} // namespace felitronics::blend

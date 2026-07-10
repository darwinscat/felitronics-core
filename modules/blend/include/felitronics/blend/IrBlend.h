// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

//==============================================================================
// felitronics::blend — the mic-blend engine (JUCE-free, headless-testable). A PURE FUNCTION of params +
// IRs: per-strip phase → HPF/LPF → shift, then Σ(gain·strip), then the Master bus. Solo/mute (cross-strip
// semantics) lives here. Promoted VERBATIM from OrbitCapture's app-local core/BlendEngine.h (de-monolith
// step 3 → promoted step 4). The bit-load-bearing no-op guards (|shift|>0.005, |mg−1|>1e-9) are ported
// VERBATIM from the old CaptureComponent so the blend is byte-identical.
//==============================================================================

#include <felitronics/blend/BlendParams.h>
#include <felitronics/blend/BlendKernels.h>

#include <vector>
#include <span>
#include <cstddef>
#include <cmath>

namespace felitronics::blend
{

inline bool anySolo (std::span<const StripParams> strips) noexcept {
    for (const auto& s : strips) if (s.solo) return true;
    return false;
}

// Solo overrides mute: if any strip is soloed, only soloed strips are audible; else audible = !mute.
inline bool channelAudible (std::span<const StripParams> strips, std::size_t m) noexcept {
    if (m >= strips.size()) return false;
    return anySolo (strips) ? strips[m].solo : ! strips[m].mute;
}

// A strip's enabled HPF then LPF (each at its active slope), applied in place.
inline void applyFilters (std::vector<float>& v, const Filter& hpf, const Filter& lpf, double sr) {
    if (const int hs = hpActiveSlope (hpf)) applyBlendSlope (v, sr, hpf.hz, true,  hs);
    if (const int ls = lpActiveSlope (lpf)) applyBlendSlope (v, sr, lpf.hz, false, ls);
}

// One channel's IR after ITS strip: phase rotation → HPF/LPF → time-shift. NO gain (summed at blend time).
inline std::vector<float> processedMic (const std::vector<float>& ir, const StripParams& s, double sr) {
    std::vector<float> v = ir;
    rotatePhase (v, s.phaseDeg);
    applyFilters (v, s.hpf, s.lpf, sr);
    if (std::abs (s.shiftMs) > 0.005) v = shiftFrac (v, s.shiftMs / 1000.0 * sr);
    return v;
}

// Σ over audible strips of (gain · processedMic), then the Master bus (gain + filters) if applyMaster.
// Returns {} on an empty / size-mismatched set (matches the old computeBlend guard).
inline std::vector<float> blendIrs (std::span<const std::vector<float>> irs,
                                    std::span<const StripParams> strips,
                                    const MasterParams& master, double sr, bool applyMaster = true) {
    if (irs.empty() || strips.size() != irs.size()) return {};
    const std::size_t len = irs[0].size();
    std::vector<float> b (len, 0.0f);
    for (std::size_t m = 0; m < irs.size(); ++m) {
        if (! channelAudible (strips, m)) continue;
        const double g = std::pow (10.0, strips[m].gainDb / 20.0);
        const auto v = processedMic (irs[m], strips[m], sr);
        for (std::size_t i = 0; i < len && i < v.size(); ++i) b[i] += (float) (g * v[i]);
    }
    if (applyMaster) {
        const double mg = std::pow (10.0, master.gainDb / 20.0);
        if (std::abs (mg - 1.0) > 1e-9) for (auto& v : b) v *= (float) mg;
        applyFilters (b, master.hpf, master.lpf, sr);
    }
    return b;
}

} // namespace felitronics::blend

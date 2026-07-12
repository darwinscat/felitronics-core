// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

//==============================================================================
// felitronics::measurement — post-processing of a derived IR: onset detection (two detectors:
// analysis-onset detectOnset + HEAD-trim detectLeadingSilence), onset-trim, peak gain.
// OFFLINE, MESSAGE-THREAD-ONLY (double, allocates; detectLeadingSilence is float by contract —
// golden trim points — and allocation-free). Non-destructive by intent — the caller decides
// whether/where to apply a gain (the raw recording keeps its true level).
//==============================================================================

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <span>
#include <vector>

namespace felitronics::measurement
{

struct OnsetResult
{
    std::size_t onset     = 0;   // trim point (attack preserved by preRoll)
    std::size_t peakIndex = 0;
    double      peak      = 0.0;
};

// Onset detection robust to leading floor / latency-gap / ring: find the PEAK in [searchStart, end),
// then walk BACK to where the signal drops below a peak-relative threshold (the attack), then back off
// preRoll. Searching back from the peak (not forward from the start) makes it immune to any pre-onset
// floor above the threshold (deconv skirt, harmonic-tail residue, mic self-noise in the latency gap).
inline OnsetResult detectOnset (std::span<const double> ir,
                                double thresholdDb = -40.0,
                                std::size_t preRoll = 8,
                                std::size_t searchStart = 0) noexcept
{
    OnsetResult r;
    if (searchStart >= ir.size()) { r.onset = ir.empty() ? 0 : ir.size() - 1; return r; }
    for (std::size_t i = searchStart; i < ir.size(); ++i)
    {
        const double a = std::fabs (ir[i]);
        if (a > r.peak) { r.peak = a; r.peakIndex = i; }
    }
    if (r.peak <= 0.0) { r.onset = searchStart; return r; }
    const double thr = r.peak * std::pow (10.0, thresholdDb / 20.0);
    std::size_t on = r.peakIndex;
    while (on > searchStart && std::fabs (ir[on - 1]) >= thr) --on;
    r.onset = (on > searchStart + preRoll) ? on - preRoll : searchStart;
    return r;
}

// The complementary detector to detectOnset — a conservative LEADING-SILENCE measure for a destructive
// HEAD trim (from OrbitCab, where the IR slot's trim and the waveform HEAD indicator share it so the
// two can't drift apart). Forward scan for the FIRST sample (on any channel) with |x| > 0.001·peak —
// peak taken across ALL channels — minus a ~0.2 ms pre-roll so the transient's leading edge isn't
// clipped; returns 0 ("don't trim") unless more than ~0.5 ms would be gained.
//
// Two detectors, two jobs:
//   * detectOnset (above)  — ANALYSIS onset: walks BACK from the peak, so a pre-onset floor above the
//                            threshold can't fool it — but that same walk can step over real low-level
//                            early energy. Right for measurement/alignment, wrong for destructive trims.
//   * detectLeadingSilence — HEAD trim: never eats early energy (everything it removes is ≤ 0.1 % of
//                            peak on EVERY channel, and the pre-roll + minimum-lead gate keep it shy);
//                            a noise floor above 0.1 % of peak stops it early — it under-trims, never
//                            over-trims. Idempotent: applying the returned trim and re-detecting gives 0.
//
// FINGERPRINT-SACRED (float by contract — golden trim points): consumers pin shipped HEAD-trim indices
// to these exact values, so the arithmetic below must stay bit-identical — float compares with a strict
// '>', TRUNCATING (int) casts for the pre-roll and the gate (44.1 kHz → 8 samples, not round()'s 9).
// channels = per-channel read pointers (AudioBuffer-shaped, e.g. getArrayOfReadPointers()); reads the
// first numSamples of each channel.
inline int detectLeadingSilence (const float* const* channels, int numChannels, int numSamples,
                                 double sampleRate) noexcept
{
    const int total = numSamples;
    const int nch   = numChannels;
    if (total <= 0 || nch <= 0 || channels == nullptr)
        return 0;

    float peak = 0.0f;                                  // == max getMagnitude across channels: the
    for (int ch = 0; ch < nch; ++ch)                    // max/|x| fold is exact, so bit-identical
        for (int i = 0; i < total; ++i)
            peak = std::max (peak, std::abs (channels[ch][i]));
    if (peak <= 0.0f)
        return 0;

    const float thresh = 0.001f * peak;
    int onset = total;                                  // min over channels of the first supra-threshold
    for (int ch = 0; ch < nch && onset > 0; ++ch)       // index (each channel scans up to the current min)
    {
        const float* d = channels[ch];
        for (int i = 0; i < onset; ++i)
            if (std::abs (d[i]) > thresh) { onset = i; break; }
    }

    const int preRoll = (int) (0.0002 * sampleRate);    // ~0.2 ms (truncated)
    const int lead    = std::max (0, onset - preRoll);
    return lead > (int) (0.0005 * sampleRate) ? lead : 0;
}

// Copy the IR from the onset forward. len == 0 → to end; else exactly len (zero-padded past the source).
inline std::vector<double> trimToOnset (std::span<const double> ir, const OnsetResult& on, std::size_t len = 0)
{
    const std::size_t start = on.onset;
    const std::size_t n = (len > 0) ? len : (ir.size() > start ? ir.size() - start : 0);
    std::vector<double> out (n, 0.0);
    for (std::size_t k = 0; k < n; ++k)
    {
        const std::size_t idx = start + k;
        if (idx < ir.size()) out[k] = ir[idx];
    }
    return out;
}

// Gain that brings the signal's peak to targetPeak (1.0 if silent). The caller applies it where it wants.
inline double peakGain (std::span<const double> x, double targetPeak = 0.98) noexcept
{
    double pk = 0.0;
    for (double v : x) if (std::isfinite (v)) pk = std::max (pk, std::fabs (v));   // ignore poisoned samples
    return (pk > 0.0) ? targetPeak / pk : 1.0;
}

inline void applyGain (std::span<double> x, double g) noexcept { for (double& v : x) v *= g; }

} // namespace felitronics::measurement

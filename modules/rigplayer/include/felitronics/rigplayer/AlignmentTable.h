// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.
#pragma once
// felitronics::rigplayer — HOW FAR APART THE MODELS OF ONE DEVICE LAND. Two captures of one pedal
// come out of training shifted by a few samples; crossfaded as they are, the pair combs. The reading is
// a relation between two files, re-measured whenever a pack is built — which is whenever a model could
// have changed — and THE PACK CARRIES IT: `lag_samples` on every file entry (NAMZ-FORMAT.md), measured
// by the capture side with every model in hand. A player takes it from there (fromDevice) and delays
// the earlier models so any pair lines up.
//
// For a pack that does not carry it, the same measurement is here: half a second of broadband through
// each model, a normalized cross-correlation against the first — a pure function over bytes, heavy
// (NAM inference per file), for the host to run on whatever thread it likes and hand to the player.
// Nothing here touches audio hardware or the player's state.

#include <felitronics/rigplayer/ModelAlignment.h>

#include <felitronics/nam/BlendLaw.h>       // kBlendMaxDelay — the most a slot can be delayed
#include <felitronics/nam/NamStage.h>
#include <namz_rig.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <map>
#include <string>
#include <vector>

namespace felitronics::rigplayer {

// The bytes of one file of the pack, by its `files[].id`. The pack lives wherever the host keeps it —
// a zip, a folder, a database row — and the player never finds out. Empty = not available.
using ModelSource = std::function<std::vector<std::byte>(const std::string& fileId)>;

/// HALF A SECOND OF BROADBAND, generated rather than fetched. It never sounds and nobody keeps it: it
/// is pushed through two models to see how far apart their outputs land, and for that any signal with
/// energy everywhere will do. Deterministic, so two runs of the same check agree. A quiet probe on
/// purpose: a model is a nonlinearity, and two of them only shift alike while neither is being driven
/// into its own compression.
inline std::vector<float> alignmentProbe(double sampleRate = 48000.0) {
    const double fs = sampleRate > 0.0 ? sampleRate : 48000.0;
    const int len = (int) std::lround(fs * 0.5);
    std::vector<float> out((std::size_t) len);
    const double f1 = 40.0, f2 = 10000.0, T = (double) len / fs;
    const double k = std::log(f2 / f1) / T;
    for (int n = 0; n < len; ++n) {
        const double t = (double) n / fs;
        const double ph = 2.0 * 3.14159265358979323846 * f1 * (std::exp(k * t) - 1.0) / k;
        out[(std::size_t) n] = (float) (0.05 * std::sin(ph));
    }
    return out;
}

// A file's output on the probe, at `sampleRate`, in 512-sample blocks. Empty when the bytes do not load.
inline std::vector<float> probeThrough(const std::vector<std::byte>& model, const std::vector<float>& probe,
                                       double sampleRate) {
    std::vector<float> y;
    if (model.empty()) return y;
    felitronics::nam::NamStage m;
    m.prepare(sampleRate, 512);
    if (! m.loadModelFromMemory(model.data(), model.size())) return y;
    y = probe;
    for (std::size_t p = 0; p < y.size(); p += 512) {
        const int n = (int) std::min<std::size_t>(512, y.size() - p);
        float* io[1] { y.data() + p };
        m.process(io, 1, n, false);
    }
    return y;
}

struct AlignmentTable {
    // Each file's LAG against the first measurable file, in samples. Positive = this file comes out
    // later than it.
    std::map<std::string, int> lagByFile;
    // The rate those samples are in: the rate the probe ran at, or 0 for a pack's own numbers, which
    // are in each MODEL's rate (the pack does not say it; the player learns it from the model it has
    // loaded). A player converts to its own rate with the three-argument delayOf.
    double sampleRate = 0.0;
    bool fromPack = false;

    bool empty() const { return lagByFile.empty(); }

    // How much to DELAY this file so it lines up with the latest one, in the table's own samples:
    // every model is held back to the slowest, so no delay is ever negative and the pair sums in
    // phase. Unknown file = no delay.
    int delayOf(const std::string& fileId) const {
        const auto me = lagByFile.find(fileId);
        if (me == lagByFile.end()) return 0;
        int latest = me->second;
        for (const auto& [id, lag] : lagByFile) latest = std::max(latest, lag);
        return std::max(0, latest - me->second);
    }

    // …and in HOST samples, given the rate the number is in (the table's, else the model's own).
    int delayOf(const std::string& fileId, double hostRate, double modelRate) const {
        const double at = sampleRate > 0.0 ? sampleRate : modelRate;
        const double scale = (at > 0.0 && hostRate > 0.0) ? hostRate / at : 1.0;
        return std::clamp((int) std::lround((double) delayOf(fileId) * scale), 0, felitronics::nam::kBlendMaxDelay);
    }

    // THE PACK'S OWN READING. Every distinct file of the stage must carry `lag_samples`, or the stage
    // is not measured and the table is empty: a pair with one side unknown cannot be lined up, and a
    // half-measured table applied as if whole delays some models against nothing.
    static AlignmentTable fromDevice(const namz::rig::Device& d) {
        AlignmentTable out;
        out.fromPack = true;
        for (const auto& f : d.files) {
            if (f.id.empty()) continue;
            if (! f.lagSamples) return AlignmentTable {};
            out.lagByFile[f.id] = *f.lagSamples;
        }
        return out;
    }
};

// The table for a device: every DISTINCT file once (several settings may point at one file), against the
// first one whose bytes load. A file with no reading — bytes that will not load, an output with no
// energy — is simply absent from the table and plays undelayed.
inline AlignmentTable measureAlignment(const namz::rig::Device& d, const ModelSource& source,
                                       double sampleRate = 48000.0) {
    AlignmentTable out;
    out.sampleRate = sampleRate > 0.0 ? sampleRate : 48000.0;
    if (! source) return out;
    const auto probe = alignmentProbe(sampleRate);
    std::vector<std::string> ids;
    for (const auto& f : d.files)
        if (! f.id.empty() && std::find(ids.begin(), ids.end(), f.id) == ids.end()) ids.push_back(f.id);
    std::vector<float> ref;
    std::string refId;
    for (const auto& id : ids) {
        const auto y = probeThrough(source(id), probe, sampleRate);
        if (y.empty()) continue;
        if (ref.empty()) { ref = y; refId = id; out.lagByFile[id] = 0; continue; }
        const auto b = felitronics::rigplayer::bestLag(ref, y, felitronics::rigplayer::kSearch);
        if (b.corr < 0.0) continue;                          // no reading is not a complaint
        out.lagByFile[id] = b.lag;
    }
    return out;
}

} // namespace felitronics::rigplayer
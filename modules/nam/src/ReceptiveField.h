// SPDX-License-Identifier: MIT
#pragma once
// HOW MUCH SIGNAL A MODEL NEEDS BEFORE ITS OUTPUT MEANS ANYTHING. A network started with empty
// buffers describes the silence it was born into for as long as its memory reaches back, so anything
// that fades a freshly loaded model in has to wait this long first.
//
// NAM answers this itself for convnet and lstm, and a container forwards to its active submodel —
// but the SlimmableContainer of WaveNets that nearly every capture is answers zero
// (`wavenet/slimmable.h`), and `Linear` inherits the base class's zero. Hence this: the receptive
// field read straight off the config — one dilated convolution at a time where there are layers to
// walk, and the plain `receptive_field` number where the architecture simply states it.

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cstddef>
#include <limits>
#include <string>

namespace felitronics::nam::detail
{
inline int receptiveFieldFromConfig (const nlohmann::json& model);

inline int receptiveFieldOfLayers (const nlohmann::json& cfg)
{
    if (! cfg.contains ("layers") || ! cfg["layers"].is_array()) return 0;
    long long total = 0;
    for (const auto& grp : cfg["layers"])
    {
        if (! grp.contains ("dilations") || ! grp["dilations"].is_array()) continue;
        const auto& ds = grp["dilations"];
        // 🔴 TWO SPELLINGS, AND NAM ACCEPTS BOTH — `kernel_sizes` per layer, or the LEGACY single
        // `kernel_size` applied to every layer (v0.5.4 `wavenet/model.cpp:914-945`, which refuses a
        // config carrying both). Reading only the array made this answer ZERO for a legacy capture,
        // and that is not a display defect: a SlimmableWavenet answers zero from NAM as well
        // (`wavenet/slimmable.h:66`), so nothing at all knew the field and an absent lane drained for
        // no samples. Measured by a crew round on a loaded slimmable model: 0.462117 out of digital
        // silence with the scalar spelling, zero with the array.
        const bool perLayer = grp.contains ("kernel_sizes") && grp["kernel_sizes"].is_array();
        const bool single   = ! perLayer && grp.contains ("kernel_size")
                           && grp["kernel_size"].is_number_integer();
        if (! perLayer && ! single) continue;
        const long long ks1 = single ? grp["kernel_size"].get<long long>() : 0;
        for (std::size_t i = 0; i < ds.size(); ++i)
        {
            if (perLayer && i >= grp["kernel_sizes"].size()) break;
            const long long k = perLayer ? grp["kernel_sizes"][i].get<long long>() : ks1;
            if (! ds[i].is_number_integer()) continue;
            total += (long long) ds[i].get<long long>() * (k - 1);
        }
    }
    // 🔴 THE CAP IS ONLY THERE TO KEEP THE int CONVERSION LEGAL, and it used to be 1<<20 — about
    // 22 seconds of field at 48 kHz, which looked like "nobody has a longer one". That was safe while
    // this number was only REPORTED; it is now also SPENT, as the length an absent lane is fed silence
    // (NamStage's drainSamples_), and a cap on a number that is spent is a silent UNDER-DRAIN: a stack
    // reaching back further than the cap replays the difference. Nine scalars build one — a single
    // layer, kernel 2, `dilations:[2000000]` — so this is constructible, not hypothetical. NAM's own
    // GetPrewarmSamples() answers 2000001 for exactly that model and the max() below already took it,
    // which is why the defect is invisible on a WaveNet and REAL on the paths where NAM answers zero.
    return total > 0 ? (int) std::min (total + 1, (long long) std::numeric_limits<int>::max()) : 0;
}

// 🔴 AND THE ARCHITECTURE THAT SIMPLY DECLARES IT. A Linear capture is an impulse response, and its
// config carries `receptive_field` as a plain number — no layers to walk. NAM does not answer for it
// either (`Linear : Buffer : DSP` inherits `GetPrewarmSamples() { return 0; }`), so before this the
// whole path reported ZERO for a capture with a field of two thousand samples: measured, every
// `receptive_field` of 1 / 2 / 65 / 257 / 2001 gave `NamStage::prewarmSamples() == 0`. That number is
// what a caller uses to size a warm-up or to drain a lane that stopped being fed, and a zero there is
// not "no memory", it is "nobody asked".
//
// 🔴 TAPS ARE NOT MEMORY, AND THIS FUNCTION OWES MEMORY. `receptive_field` counts TAPS:
// `y[n] = Σ_{k<RF} h[k]·x[n−k]`, so the last real input at time T is out of the window from output
// T+RF onward — RF−1 further samples. A one-tap capture is a GAIN and has no memory at all, and
// reporting 1 for it is not a harmless rounding: `rigplayer::RigPlayer::warmFor()` returns early on
// `pre <= 0`, so a 1 there turns a memoryless capture's warm-up from nothing into
// `1 + latency + maxBlock` and a slot that wakes goes silent for a block. That regression was caught
// by the repository's own suite (RigPlayerTests, "a player that sleeps sounds bit-identically to one
// that never does") on a crew round, not by reasoning. The layers path above returns TAPS as well
// (`total + 1`), and NAM's ConvNet does the same (`1 + Σ dilations`); that one sample of margin is
// harmless there and is deliberately left alone rather than swept into this change.
inline int declaredReceptiveField (const nlohmann::json& cfg)
{
    // is_number_integer() rather than get<>() in a try: this whole file is called from inside
    // prepareModel's catch-all, so a throw here would turn a legal-but-odd config into a REFUSED load.
    if (! cfg.contains ("receptive_field") || ! cfg["receptive_field"].is_number_integer()) return 0;
    const long long rf = cfg["receptive_field"].get<long long>();
    return rf > 1 ? (int) std::min (rf - 1, (long long) std::numeric_limits<int>::max()) : 0;
}

// 🔴 AND THE PART OF THE MEMORY THAT IS NOT THE FIELD AT ALL. A Linear capture is not always run tap
// by tap: NAM picks a partitioned-FFT convolution for anything past 256 taps unless the config says
// otherwise (`implementation`, default `auto` — NAM v0.5.4 `linear.cpp:14-17, 100-108`), and that
// engine holds input SPECTRA in a ring, so a lane goes on emitting for up to two of its blocks after
// its field has been flushed. Measured by a crew round on a dense 2001-tap kernel at 48 kHz: a drain
// of exactly the field left 1.909e-08 on 46 samples of the return, and the same kernel with
// `"implementation":"direct"` left nothing — i.e. a sparse or direct fixture is blind to it.
// The block is 256 / 512 / 1024 taps for a field of ≤ 2048 / ≤ 8192 / more, so two of the largest is
// the bound for every case at this pin. It is a RESTATEMENT of a third-party constant and is written
// as one: named, pinned to the tag, and gated by a test with a dense kernel, because there is no
// owner here to ask.
inline bool anySubmodel (const nlohmann::json& model, bool (*pred) (const nlohmann::json&));

inline bool isLinearArchitecture (const nlohmann::json& model)
{
    return model.contains ("architecture") && model["architecture"].is_string()
        && model["architecture"].get<std::string>() == "Linear";
}

inline int partitionedTailSamples (const nlohmann::json& model)
{
    // A CONTAINER IS ASKED THROUGH, exactly as the field is. `receptiveFieldFromConfig` recurses into
    // `submodels` and this must too, or a container of Linear captures reports the top-level
    // architecture (SlimmableContainer), is charged no ring, and drains short on whichever submodel is
    // speaking — measured on a loaded one, 1.48e-08 out of digital silence.
    // It is charged whatever the `implementation` key says, including "direct", which cannot use a ring
    // at all: an over-drain of 2048 samples is a few milliseconds of one lane once per departure, and
    // the alternative is this file second-guessing NAM's own `auto` rule (`linear.cpp:100-108`), which
    // is exactly the restatement rule 9u exists against.
    return (isLinearArchitecture (model) || anySubmodel (model, isLinearArchitecture)) ? 2 * 1024 : 0;
}

// Whether the architecture carries a RECURRENT cell, whose state no finite length of silence empties.
// It is asked because NAM's own answer for one is a heuristic on a tag that may be absent: an LSTM
// with no `sample_rate` reports `GetPrewarmSamples() == 1` (`lstm.cpp:125-131`: 0.5 × −1 ≤ 0 → 1), so
// a drain sized from it is ONE SAMPLE. Measured by a crew round on a τ ≈ 22 000-sample cell: tagged,
// the truncated drain left 0.419 against 0.023 for a lane clocked throughout; untagged, 0.499.
inline bool isLstmArchitecture (const nlohmann::json& model)
{
    return model.contains ("architecture") && model["architecture"].is_string()
        && model["architecture"].get<std::string>() == "LSTM";
}

inline bool isRecurrent (const nlohmann::json& model)
{
    return isLstmArchitecture (model) || anySubmodel (model, isLstmArchitecture);
}

// True when ANY submodel, at any depth, satisfies `pred` — the container switches by level and any of
// them can be the one speaking, so the answer has to be the worst case over all of them, which is the
// same rule receptiveFieldFromConfig applies to the field.
inline bool anySubmodel (const nlohmann::json& model, bool (*pred) (const nlohmann::json&))
{
    if (! model.contains ("config")) return false;
    const auto& cfg = model["config"];
    if (! cfg.contains ("submodels") || ! cfg["submodels"].is_array()) return false;
    for (const auto& sub : cfg["submodels"])
        if (sub.contains ("model") && (pred (sub["model"]) || anySubmodel (sub["model"], pred)))
            return true;
    return false;
}

inline int receptiveFieldFromConfig (const nlohmann::json& model)
{
    if (! model.contains ("config")) return 0;
    const auto& cfg = model["config"];
    // A container holds several models and switches between them by level; any of them can be the
    // one speaking, so the longest memory is the one that has to be waited out.
    if (cfg.contains ("submodels") && cfg["submodels"].is_array())
    {
        int worst = 0;
        for (const auto& sub : cfg["submodels"])
            if (sub.contains ("model"))
                worst = std::max (worst, receptiveFieldFromConfig (sub["model"]));
        return worst;
    }
    // The layers are asked FIRST and the declaration is the fallback, not the other way round: a
    // WaveNet config may carry both, and the dilated stack is the one that describes what it does.
    if (const int fromLayers = receptiveFieldOfLayers (cfg); fromLayers > 0) return fromLayers;
    return declaredReceptiveField (cfg);
}
} // namespace felitronics::nam::detail

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
//
// 🔴 WHAT THIS FILE PROMISES IS AN UPPER BOUND ON THE MEMORY OF THE WHOLE MODEL, NOT THE FIELD OF ITS
// STACK. A .nam config nests whole models in two places, and NAM builds both by handing the sub-node
// straight back to `get_dsp()`:
//   `config.submodels[i].model` — a container, one of which is SPEAKING          (container.cpp:163)
//   `config.condition_dsp`      — the CONDITIONER, which is always RUNNING       (wavenet/model.cpp:844)
// Those two are the complete list at this pin: they are the only `get_dsp()` calls in the library
// outside `get_dsp.cpp` itself, and the third one (`wavenet/slimmable.cpp:442`) rebuilds that SAME
// `condition_dsp` node. A container is answered for with the WORST of its submodels, because any of
// them can be the one playing. A conditioner is answered for with a SUM, because it is in SERIES:
// `_process_condition` runs the raw input through it and the layer arrays are then processed against
// its output (`wavenet/model.cpp:699-729, :749-761`), so the network's window reaches back through
// the conditioner's. NAM's own arithmetic composes them the same way — `mPrewarmSamples` starts at
// the conditioner's prewarm and ADDS the stack (`wavenet/model.cpp:616-620`).
//
// 🔴 THE SUM IS AN UPPER BOUND AND NOT THE EXACT REACH, AND THE SLACK IS NAMED. The condition enters
// each layer AFTER that layer's dilated convolution (`z = conv(input) + input_mixin(condition)`,
// `wavenet/model.cpp:196-204`; the mixin is a memoryless Conv1x1, `wavenet/detail.h:47-49`), so the
// deepest path out of the conditioner skips the FIRST convolution's own lookback L₀. The true reach
// is `Mₒ + H + max(0, M_c − L₀)` and this file answers `Mₒ + H + M_c + 1`: over by `L₀ + 1`, never
// short. That formula is not reasoning, it is measured — an impulse against digital silence, bit
// compared, four shapes including the branch where the conditioner is SHORTER than L₀:
//     dilations {1,64,256}, no conditioner                      reach  321, answered  322
//     …with a 1000-sample conditioner                           reach 1320, answered 1322
//     …and a head kernel of 9 on top of that                    reach 1328, answered 1330
//     dilations {500,64} with a 100-sample conditioner          reach  564, answered  665
// (the last one is `max(0, M_c − L₀) == 0`: a conditioner that reaches back less far than the first
// convolution adds nothing to the model's memory, and the over-estimate there is the whole M_c.)
// `conv_pre_film` would put L₀ in series too (`wavenet/model.cpp:172-177`), which only makes the sum
// tighter, never wrong.
//
// With SEVERAL layer arrays the bound is looser again, and by the same amount NAM's own arithmetic is:
// an array's head outputs are copied into the NEXT array's head accumulator (`model.cpp:434-447`), so
// the head rechannels are in SERIES, while the residual path into the next array's convolutions
// BYPASSES this array's head. The true reach across a junction is `M_A + max(H_A, M_B) + H_B`, and both
// this file and `mPrewarmSamples` sum the two terms instead of taking the larger. Measured on three
// two-array stacks — {d=1,hk=2}+{d=1,hk=2} reaches 3 and both answer 5; {d=1,hk=64}+{d=100,hk=2}
// reaches 102 and both answer 166; {d=1,hk=100}+{d=8,hk=2} reaches 101 and both answer 110. Over in
// every case, and never looser than the number NAM would have used on its own.
//
// ⚠️ ONE SHAPE THIS FILE STILL ANSWERS ZERO FOR — registered rather than guessed at, and named as ONE
// rather than as the end of a list, because nobody has proved the list closed. The HYBRID SLIMMABLE
// WRAPPER: there is no registered `"SlimmableWavenet"` architecture, so the route in is `"WaveNet"`,
// whose parser delegates when a TOP-LEVEL `config.layers[i].slimmable.method` marker is present
// (`wavenet/model.cpp:1205-1229`), and the parser it delegates to then unwraps the REAL config from
// `config.model` (`wavenet/slimmable.cpp:543`). A config carrying both the marker and `config.model`
// therefore loads, and everything this file reads is in the half it does not look at. The hole is the
// model's WHOLE field, not a corner of it: measured on one whose inner stack has a single dilation of
// 200, NAM answers 0, this file answers 0, and the impulse reaches sample 200. Reading it means
// restating NAM's own dispatch heuristic in order to know when `config.model` is the config at all,
// which is the restatement rule 9u exists against — so it is a number in the plan, not a branch here.
// No capture in NAM's own `example_models/` is that shape.

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>

namespace felitronics::nam::detail
{
inline int receptiveFieldFromConfig (const nlohmann::json& model);

// 🔴 A JSON NUMBER READ AS A SAMPLE COUNT — AS NAM READS IT, which is `get<int>()` on every one of
// these keys, and that fixes all three of the decisions here.
//
// NOT `is_number_integer()`: `get<int>()` static_casts any arithmetic node, so a capture spelling its
// field `2001.0` LOADS and a narrow guard would answer zero for it and drain for nothing (measured —
// the narrow guard also broke the `dilations` reader against base, which never guarded at all: base
// answered 3 for `"dilations":[2.0]` and the narrow version 0).
//
// A BOOLEAN IS A NUMBER TO `get<int>()` AND NOT TO `is_number()`, and that gap was a real under-drain.
// Measured on nlohmann 3.12: `get<int>(true)` is 1 while `get<long long>(true)` and `get<double>(true)`
// THROW. So NAM loads `"dilations":[true,true]` as `[1,1]` — verified, its own prewarm answers 3 — and
// this file read it as nothing: on a slimmable stack, where NAM's answer is zero too, that left a lane
// replaying two samples at 0.1875 out of digital silence.
//
// AND A VALUE THAT DOES NOT FIT AN `int` IS REFUSED, because the model NAM built does not contain it.
// NAM's dilation IS an int; `"dilations":[4294967396]` builds a network whose dilation is 100, and a
// file like that loads. Casting a float out of range is undefined behaviour besides ([conv.fpint]/1),
// with a hard-fail UBSan job downstream. Refusing reads as ABSENT — the same answer this file already
// gives a field spelled as a string — and NAM's own number then covers whatever it did build. The
// alternative, clamping to INT_MAX, is an upper bound that is true and useless: it spends 2 147 483 646
// samples of inference, a measured 23.7 s of synchronous audio-thread work, for a model that remembers
// a hundred samples.
inline bool readCount (const nlohmann::json& value, long long& out)
{
    constexpr long long kMax = (long long) std::numeric_limits<int>::max();
    constexpr long long kMin = (long long) std::numeric_limits<int>::min();
    if (value.is_boolean()) { out = value.get<bool>() ? 1 : 0; return true; }
    if (! value.is_number()) return false;
    if (value.is_number_integer() || value.is_number_unsigned())
    {
        if (value.is_number_unsigned() && value.get<std::uint64_t>() > (std::uint64_t) kMax) return false;
        const long long v = value.get<long long>();
        if (v > kMax || v < kMin) return false;
        out = v;
        return true;
    }
    const double d = value.get<double>();
    if (! (d >= (double) kMin && d <= (double) kMax)) return false;   // a NaN fails both, which is the point
    out = (long long) d;
    return true;
}

// Every count this file spends is bounded by the int it must survive as. Clamping BEFORE the products
// and the sums below is what keeps them from signed overflow, which is undefined behaviour and not a
// large number.
inline long long clampCount (long long v)
{
    return std::min (std::max (v, 0LL), (long long) std::numeric_limits<int>::max());
}

// 🔴 AND THE HEAD IS MEMORY TOO. A layer array ends in `_head_rechannel`, a causal Conv1D whose
// kernel is the layer's own `head.kernel_size`, and NAM charges `kernel − 1` for it on top of the
// dilations (`LayerArray::get_receptive_field`, v0.5.4 `wavenet/model.cpp:417-423`; parsed at
// `:882-899`, where the LEGACY `head_size` spelling means an implicit kernel of 1 and therefore
// nothing). The real captures carry it: `example_models/A2.nam` spells `head.kernel_size` 16 in both
// submodels, so this file answered 6332 where NAM answers 6347 and the ARCHITECTURE reaches back 6346.
// (The architecture, not the impulse: that capture's weights attenuate, and its measured impulse
// response dies at sample 5899 — which is why the gate for this row is NAM's own number and not a
// measurement, the one place in this file where those two differ.)
// That shortfall is normally hidden by `NamStage.cpp:78`, which raises this number with NAM's own —
// but a SlimmableWavenet answers zero (`wavenet/slimmable.h:66`) and there is nothing to raise it
// with. Measured on one: a single dilation-1 layer with `head.kernel_size` 16 loads, NAM answers 0,
// this file answered 2, and the impulse reaches sample 16.
inline int headRechannelReach (const nlohmann::json& grp)
{
    // A group with no nested `head` object is the legacy spelling — `head_size` + `head_bias`, kernel
    // 1, no memory. `is_number()`, not `is_number_integer()`, for the reason declaredReceptiveField
    // states below: NAM's `get<int>()` takes `16.0` and a narrow guard would silently answer zero.
    if (! grp.contains ("head") || ! grp["head"].is_object()) return 0;
    const auto& head = grp["head"];
    if (! head.contains ("kernel_size")) return 0;
    long long k = 0;
    if (! readCount (head["kernel_size"], k)) return 0;
    return k > 1 ? (int) clampCount (k - 1) : 0;
}

// …AND THE POST-STACK HEAD, which is a different key and a different shape: `config.head` is a chain
// of causal convolutions over the accumulated head outputs, `receptive_field()` is `1 + Σ(kᵢ − 1)`
// and NAM charges it minus one (`wavenet/model.cpp:58-67`, charged at `:620`, parsed at `:1154-1186`).
// `"head": null` means ABSENT, exactly as NAM reads it (`:1154`).
inline long long postStackHeadReach (const nlohmann::json& cfg)
{
    if (! cfg.contains ("head") || ! cfg["head"].is_object()) return 0;
    const auto& head = cfg["head"];
    if (! head.contains ("kernel_sizes") || ! head["kernel_sizes"].is_array()) return 0;
    long long total = 0;
    for (const auto& entry : head["kernel_sizes"])
    {
        long long k = 0;
        if (readCount (entry, k))
            total = clampCount (total + clampCount (k - 1));
    }
    return total;
}

inline int receptiveFieldOfLayers (const nlohmann::json& cfg)
{
    if (! cfg.contains ("layers") || ! cfg["layers"].is_array()) return 0;
    long long total = 0;
    for (const auto& grp : cfg["layers"])
    {
        // The head rechannel is charged for the ARRAY, not for a layer, so it is added before the two
        // `continue`s below, which only decide whether this group's CONVOLUTIONS can be read. That
        // placement is DEFENSIVE and is said so rather than claimed: no config that `get_dsp` accepts
        // distinguishes it, because NAM throws on every shape that would reach a `continue` — absent or
        // non-array `dilations` (`model.cpp` inserts a null and the `std::vector<int>` conversion
        // throws) and a missing or non-numeric `kernel_size` alike. It costs nothing and it means a
        // future spelling this reader cannot parse still pays for the head it can.
        total = clampCount (total + headRechannelReach (grp));
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
        long long ks1 = 0;
        const bool single = ! perLayer && grp.contains ("kernel_size")
                         && readCount (grp["kernel_size"], ks1);   // readCount: NAM takes 2.0 as well
        if (! perLayer && ! single) continue;
        for (std::size_t i = 0; i < ds.size(); ++i)
        {
            if (perLayer && i >= grp["kernel_sizes"].size()) break;
            long long k = ks1, d = 0;
            if (perLayer && ! readCount (grp["kernel_sizes"][i], k)) continue;
            if (! readCount (ds[i], d)) continue;
            total = clampCount (total + clampCount (d) * clampCount (k - 1));
        }
    }
    total = clampCount (total + postStackHeadReach (cfg));
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

// 🔴 AND THE ARCHITECTURE THAT KEEPS ITS STACK AT THE TOP LEVEL. A ConvNet is a chain of dilated
// convolutions with a kernel NAM hard-codes to 2 ("HACK 2 kernel", v0.5.4 `convnet.cpp:56-57`), and
// its field is `1 + Σ dilations` off a TOP-LEVEL `dilations` array (`convnet.cpp:200-202`) — the key
// a WaveNet keeps inside its `layers` groups. Nothing here read that shape, so every ConvNet answered
// ZERO: measured on `dilations:[1,2,4,8]`, this file answered 0 where NAM answers 16 and the impulse
// reaches sample 15. NAM's own answer hides it at the top level and behind a plain WaveNet's
// conditioner; it does not hide it behind a SlimmableWavenet, which answers zero for everything.
inline int convNetField (const nlohmann::json& cfg)
{
    // A WaveNet keeps its dilations one level down, inside `layers`; reading a top-level `dilations`
    // on a config whose `layers` this file can walk would be reading a key NAM never looks at there.
    // 🔴 A NON-EMPTY ARRAY, and not merely the KEY's presence: NAM's ConvNet parser never reads
    // `layers` (`convnet.cpp:349-355`), so a ConvNet carrying a dead `"layers": []` had its whole stack
    // suppressed here — the same shape of mistake as a chain that lets one reader silence another.
    // Measured on one that loads: as a slimmable capture's conditioner, where NAM answers zero and this
    // file is the only answer, the lane drained 102 for a model reaching 131 and handed back 29 samples.
    if (cfg.contains ("layers") && cfg["layers"].is_array() && ! cfg["layers"].empty()) return 0;
    if (! cfg.contains ("dilations") || ! cfg["dilations"].is_array()) return 0;
    long long total = 0;
    for (const auto& entry : cfg["dilations"])
    {
        long long d = 0;
        if (readCount (entry, d))
            total = clampCount (total + clampCount (d));
    }
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
    // 🔴 is_number(), NOT is_number_integer(), and the difference is measured. A guard is wanted at all
    // because this file runs inside prepareModel's catch-all and a throw there turns a legal-but-odd
    // config into a REFUSED load. But `2.0` is a json number that is NOT an integer, and NAM's own
    // parser takes it — `get<int>()` static_casts any arithmetic node — so a capture spelling its field
    // `2001.0` LOADS and would then drain for nothing. The narrow guard also broke the `dilations`
    // reader against base, which never guarded at all: base answered 3 for `"dilations":[2.0]` and the
    // narrow version answered 0. Caught by a pre-merge round, not by the suite.
    if (! cfg.contains ("receptive_field")) return 0;
    long long rf = 0;
    if (! readCount (cfg["receptive_field"], rf)) return 0;
    return rf > 1 ? (int) clampCount (rf - 1) : 0;
}

// THE CONDITIONER, AND THE SHAPE GATE THAT SAYS WHEN IT IS ONE. `condition_dsp` is read by exactly
// one parser in the library — `nam::wavenet::parse_config_json` (v0.5.4 `wavenet/model.cpp:841-852`)
// — and it is DEAD JSON on every other architecture: a Linear (`linear.cpp:324`), an LSTM
// (`lstm.cpp:188`), a ConvNet (`convnet.cpp:349`) or a container (`container.cpp:146-168`) carrying
// the key builds nothing from it, so charging one there would be inventing memory that no instance
// has. The gate is the same shape NAM dispatches on, a `layers` ARRAY, which is also what keeps a
// CONTAINER's own sibling `condition_dsp` uncharged — a decision the panel was unanimous on.
// `null` means absent, exactly as NAM reads it (`:841`, `!is_null()`).
inline const nlohmann::json* conditionerOf (const nlohmann::json& cfg)
{
    if (! cfg.contains ("layers") || ! cfg["layers"].is_array()) return nullptr;
    if (! cfg.contains ("condition_dsp")) return nullptr;
    const auto& cond = cfg["condition_dsp"];
    return cond.is_object() ? &cond : nullptr;
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
inline bool anyNestedModel (const nlohmann::json& model, bool (*pred) (const nlohmann::json&));

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
    // …AND SO IS A CONDITIONER, for the same reason and with the same engine: the conditioner NAM
    // builds is a real `nam::Linear` (`wavenet/model.cpp:844`), `_configure_implementation` runs in its
    // constructor and again on every `SetMaxBufferSize`, which the WaveNet forwards to it
    // (`wavenet/model.cpp:658`), and its ring feeds the layer arrays through the mixin. One charge is
    // enough for the whole tree: a Linear is a leaf, and the rings do not compose in series here.
    // It is charged whatever the `implementation` key says, including "direct", which cannot use a ring
    // at all: an over-drain of 2048 samples is a few milliseconds of one lane once per departure, and
    // the alternative is this file second-guessing NAM's own `auto` rule (`linear.cpp:100-108`), which
    // is exactly the restatement rule 9u exists against.
    return (isLinearArchitecture (model) || anyNestedModel (model, isLinearArchitecture)) ? 2 * 1024 : 0;
}

// Whether the architecture carries a RECURRENT cell, whose state no finite length of silence empties.
// It is asked because NAM's own answer for one is a heuristic on a tag that may be absent: an LSTM
// with no `sample_rate` reports `GetPrewarmSamples() == 1` (`lstm.cpp:125-131`: 0.5 × −1 ≤ 0 → 1), so
// a drain sized from it is ONE SAMPLE. Measured by a crew round on a τ ≈ 22 000-sample cell: tagged,
// the truncated drain left 0.419 against 0.023 for a lane clocked throughout; untagged, 0.499.
// A CONDITIONER COUNTS: an LSTM conditioner's cell enters every layer's `z` through the memoryless
// mixin, so the model as a whole carries recurrent state and no finite silence empties it either.
inline bool isLstmArchitecture (const nlohmann::json& model)
{
    return model.contains ("architecture") && model["architecture"].is_string()
        && model["architecture"].get<std::string>() == "LSTM";
}

inline bool isRecurrent (const nlohmann::json& model)
{
    return isLstmArchitecture (model) || anyNestedModel (model, isLstmArchitecture);
}

// True when ANY model NESTED in this one, at any depth, satisfies `pred` — the two nesting keys are
// the ones named at the top of this file, and BOTH are walked here because a predicate that reached
// only half the tree is the same defect in a different function. The container switches by level and
// any of them can be the one speaking; the conditioner is always running. Either way the answer has
// to be the worst case over the tree, which is the same rule receptiveFieldFromConfig applies to the
// field.
inline bool anyNestedModel (const nlohmann::json& model, bool (*pred) (const nlohmann::json&))
{
    if (! model.contains ("config")) return false;
    const auto& cfg = model["config"];
    if (const nlohmann::json* cond = conditionerOf (cfg))
        if (pred (*cond) || anyNestedModel (*cond, pred))
            return true;
    if (! cfg.contains ("submodels") || ! cfg["submodels"].is_array()) return false;
    for (const auto& sub : cfg["submodels"])
        if (sub.contains ("model") && (pred (sub["model"]) || anyNestedModel (sub["model"], pred)))
            return true;
    return false;
}

inline int receptiveFieldFromConfig (const nlohmann::json& model)
{
    if (! model.contains ("config")) return 0;
    const auto& cfg = model["config"];
    // A container holds several models and switches between them by level; any of them can be the
    // one speaking, so the longest memory is the one that has to be waited out.
    long long worst = 0;
    if (cfg.contains ("submodels") && cfg["submodels"].is_array())
        for (const auto& sub : cfg["submodels"])
            if (sub.contains ("model"))
                worst = std::max (worst, (long long) receptiveFieldFromConfig (sub["model"]));
    // 🔴 AND THE CONTAINER BRANCH DOES NOT RETURN HERE, because NAM dispatches on the `architecture`
    // STRING and never on shape (`get_dsp.cpp:261`): a `"WaveNet"` carrying a stray `submodels` array
    // LOADS, and every key this file reads is in the half a `return` would skip. Measured on one — a
    // slimmable WaveNet, so NAM's own answer is zero and there is nothing to raise this with — the
    // early return answered 1 for a model reaching back 4200. Taking the worst of the two instead is
    // also what keeps the three functions answering about the SAME tree: `anyNestedModel` walks the
    // conditioner BEFORE it looks at `submodels`, so a return here made `isRecurrent` charge a
    // conditioner that the field ignored. Measured on a container carrying a stray `layers` key and a
    // dead LSTM conditioner: the field answered 501 and `isRecurrent` answered TRUE, which costs the
    // reader a 24 000-sample floor and a restart that is never idempotent again, for a tree whose only
    // instance reaches back 500.
    //
    // The layers are asked FIRST and the two declarations are the fallback, not the other way round: a
    // WaveNet config may carry `receptive_field` as well, and the dilated stack is the one that
    // describes what it does — maxing there would drain a real capture for a stale number. The other
    // two are MAXED against each other rather than chained, because a chain lets the ConvNet reader
    // SUPPRESS the declared one: a `Linear` carrying a stray `dilations` array loads (its parser reads
    // neither key, `linear.cpp:306-316`) and a first-non-zero chain answered 2 for a 4999-sample
    // impulse response that answered 4999 without the stray key. A max cannot take anything away.
    long long own = receptiveFieldOfLayers (cfg);
    if (own == 0) own = std::max ((long long) convNetField (cfg), (long long) declaredReceptiveField (cfg));
    // …AND THE CONDITIONER IS ADDED TO IT, because it is in series (see the top of this file). The sum
    // is taken in long long: both addends are already individually clamped to INT_MAX, so an int
    // addition here is undefined behaviour for a config a unit test can write down.
    const nlohmann::json* cond = conditionerOf (cfg);
    const long long total = std::max (worst,
                                      own + (cond != nullptr ? (long long) receptiveFieldFromConfig (*cond) : 0LL));
    return (int) std::min (total, (long long) std::numeric_limits<int>::max());
}
} // namespace felitronics::nam::detail

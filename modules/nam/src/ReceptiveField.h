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
// walk, the plain `receptive_field` number where the architecture simply states it, and a named
// CEILING where the config is a shape this file cannot place at all.
//
// 🔴 WHAT THIS FILE PROMISES IS AN UPPER BOUND ON THE MEMORY OF THE WHOLE MODEL, NOT THE FIELD OF ITS
// STACK. A .nam config nests whole models in two places, and NAM builds both by handing the sub-node
// straight back to `get_dsp()`:
//   `config.submodels[i].model` — a container, one of which is SPEAKING          (container.cpp:163)
//   `config.condition_dsp`      — the CONDITIONER, which is always RUNNING       (wavenet/model.cpp:844)
// Those two are the complete list at this pin: they are the only `get_dsp()` calls in the library
// outside `get_dsp.cpp` itself, and the third one (`wavenet/slimmable.cpp:442`) rebuilds that SAME
// `condition_dsp` node. (A THIRD place a whole config can hang — the slimmable wrapper's `config.model`,
// handed to `parse_config_json` rather than to `get_dsp()` — is the reason for the ceiling rule below,
// and it is deliberately NOT added to this list by name: see `forEachUnplacedConfig`.)
// A container is answered for with the WORST of its submodels, because any of
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
// 🔴 AND WHAT IT ANSWERS FOR A CONFIG IT CANNOT PLACE: A CEILING, NEVER ZERO. This file promises a
// BOUND and not an exact value, so "I do not know" has a natural answer and it is not zero. The costs
// are asymmetric and so is the rule: an understated number is the tail of the previous sound coming
// out of digital silence, an overstated one is inference nobody hears. `kUnreadShapeCeiling` below is
// that answer, with its measured price beside it.
//
// WHAT MADE THE RULE NECESSARY — measured on a REAL SHIPPED CAPTURE, not a toy. There is no registered
// `"SlimmableWavenet"` architecture, so the route in is `"WaveNet"`, whose parser delegates when a
// TOP-LEVEL `config.layers[i].slimmable.method` marker is present (`wavenet/model.cpp:1205-1229`), and
// the parser it delegates to then takes the REAL config from `config.model` when that key is there and
// from the config itself when it is not (`wavenet/slimmable.cpp:544`). A config carrying both the
// marker and `config.model` therefore loads, and everything this file used to read was in the half NAM
// ignores. Take NAM's own `example_models/slimmable_wavenet.nam`, move its config under `config.model`,
// leave a top-level `layers` carrying nothing but the marker, and it is the same capture making the
// same sound: the impulse reaches sample 2046 either way, NAM's own `GetPrewarmSamples()` is `return 0`
// for the architecture (`wavenet/slimmable.h:66`) so `NamStage.cpp:78`'s max() raises nothing, and this
// file answered 2047 for the flat file and ZERO for the wrapped one. All three readers were blind
// together: with an LSTM conditioner inside the wrapper `isRecurrent` went true → FALSE while the
// impulse never died, and with a dense `Linear` conditioner the ring went 2048 → 0 and the field
// 4047 → 0 for a model reaching 4042.
//
// AND THE RULE IS NOT A PATCH FOR THAT ONE ADDRESS. The same class — a config whose memory this file
// cannot derive, answered with a number anyway — had further members, all measured on loaded models:
//   a `Linear` of 5000 dense taps carrying a READABLE stray `layers` array answered 2 for an impulse
//   reaching 4999, because the stack reading silenced the declared field;
//   a `ConvNet` carrying a NON-EMPTY dead `layers` array answered 0 where NAM's field is 256, because
//   `convNetField` suppressed itself for it;
//   a SLIMMABLE WaveNet spelling a dilation past the int (`4294967396`, which NAM builds as 100)
//   answered 0, because a refused value read as absent.
// One rule closes all of them, and it is stated where the numbers meet (`placedField` and the block
// above `kUnreadShapeCeiling`): WHAT THIS FILE PLACES, IT TRUSTS; FOR WHAT IT CANNOT PLACE — AN UNREAD
// CONFIG, A REFUSED VALUE, A READING IT SETS ASIDE — IT CHARGES ONE ALLOWANCE, NEVER THE FACE VALUE AND
// NEVER NOTHING. The review rounds that shaped it measured both failures of the simpler rules: face
// value turned a few dead bytes into a half-hour `reset()`, and an allowance charged per node turned a
// 1.6 MB file into INT_MAX.

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
// with a hard-fail UBSan job downstream. The alternative, clamping to INT_MAX, is an upper bound that is
// true and useless: it spends 2 147 483 646 samples of inference, a measured 23.7 s of synchronous
// audio-thread work, for a model that remembers a hundred samples.
// 🔴 BUT A REFUSAL IS NOT "ABSENT" (P92). This used to say it was, "and NAM's own number then covers
// whatever it did build" — which is false exactly where this file is the only answer: a SLIMMABLE
// WaveNet spelling `"dilations":[4294967396]` loads, NAM builds dilation 100, the impulse reaches 100,
// NAM answers `return 0`, and this file answered 0. A value that is THERE and cannot be read is "I do
// not know", and every caller below passes a `refused` flag that turns it into the allowance
// (`kUnreadShapeCeiling`) — never the value, never nothing.
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

// `readCount`, and a refusal of a value that is PRESENT is written down — see the note above it.
inline bool readCountOrNote (const nlohmann::json& value, long long& out, bool* refused)
{
    if (readCount (value, out)) return true;
    if (refused != nullptr) *refused = true;
    return false;
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
inline int headRechannelReach (const nlohmann::json& grp, bool* refused = nullptr)
{
    // A group with no nested `head` object is the legacy spelling — `head_size` + `head_bias`, kernel
    // 1, no memory. `is_number()`, not `is_number_integer()`, for the reason declaredReceptiveField
    // states below: NAM's `get<int>()` takes `16.0` and a narrow guard would silently answer zero.
    if (! grp.contains ("head") || ! grp["head"].is_object()) return 0;
    const auto& head = grp["head"];
    if (! head.contains ("kernel_size")) return 0;
    long long k = 0;
    if (! readCountOrNote (head["kernel_size"], k, refused)) return 0;
    return k > 1 ? (int) clampCount (k - 1) : 0;
}

// …AND THE POST-STACK HEAD, which is a different key and a different shape: `config.head` is a chain
// of causal convolutions over the accumulated head outputs, `receptive_field()` is `1 + Σ(kᵢ − 1)`
// and NAM charges it minus one (`wavenet/model.cpp:58-67`, charged at `:620`, parsed at `:1154-1186`).
// `"head": null` means ABSENT, exactly as NAM reads it (`:1154`).
inline long long postStackHeadReach (const nlohmann::json& cfg, bool* refused = nullptr)
{
    if (! cfg.contains ("head") || ! cfg["head"].is_object()) return 0;
    const auto& head = cfg["head"];
    if (! head.contains ("kernel_sizes") || ! head["kernel_sizes"].is_array()) return 0;
    long long total = 0;
    for (const auto& entry : head["kernel_sizes"])
    {
        long long k = 0;
        if (readCountOrNote (entry, k, refused))
            total = clampCount (total + clampCount (k - 1));
    }
    return total;
}

inline int receptiveFieldOfLayers (const nlohmann::json& cfg, bool* refused = nullptr)
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
        total = clampCount (total + headRechannelReach (grp, refused));
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
                         && readCountOrNote (grp["kernel_size"], ks1, refused);   // NAM takes 2.0 as well
        if (! perLayer && ! single) continue;
        for (std::size_t i = 0; i < ds.size(); ++i)
        {
            if (perLayer && i >= grp["kernel_sizes"].size()) break;
            long long k = ks1, d = 0;
            if (perLayer && ! readCountOrNote (grp["kernel_sizes"][i], k, refused)) continue;
            if (! readCountOrNote (ds[i], d, refused)) continue;
            total = clampCount (total + clampCount (d) * clampCount (k - 1));
        }
    }
    total = clampCount (total + postStackHeadReach (cfg, refused));
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
inline int convNetField (const nlohmann::json& cfg, bool* refused = nullptr)
{
    // 🔴 NO `layers` GATE HERE ANY MORE — the decision it made is taken where the three readings meet.
    // A gate stood here: "a WaveNet keeps its dilations inside `layers`, so a top-level `dilations`
    // beside real layers is a key NAM never looks at there". P87 took it back for an EMPTY `layers` (a
    // ConvNet carrying `"layers": []` drained 102 for a model reaching 131, as a slimmable capture's
    // conditioner); it still fired for a NON-EMPTY one, and NAM's ConvNet parser never reads `layers`
    // (`convnet.cpp:326-338`), so a ConvNet carrying `"layers":[{}]` LOADS and this answered 0 where
    // NAM's field is 256. A reader silencing another is the chain mistake inside one function; the
    // junction in `placedField` now says what happens when two readings disagree.
    if (! cfg.contains ("dilations") || ! cfg["dilations"].is_array()) return 0;
    long long total = 0;
    for (const auto& entry : cfg["dilations"])
    {
        long long d = 0;
        if (readCountOrNote (entry, d, refused))
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
inline int declaredReceptiveField (const nlohmann::json& cfg, bool* refused = nullptr)
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
    if (! readCountOrNote (cfg["receptive_field"], rf, refused)) return 0;
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

//======================================================================================================
// 🔴 WHAT THIS FILE ANSWERS WHERE IT DOES NOT KNOW: WHAT IT READ, PLUS ONE ALLOWANCE. NEVER ZERO.
//
// THE RULE, in one line: WHAT THIS FILE PLACES, IT TRUSTS; FOR ANYTHING IT CANNOT PLACE, IT CHARGES ONE
// ALLOWANCE — NOT THE THING'S FACE VALUE, AND NOT NOTHING. "Cannot place" is three events, all read off
// the config and none off NAM's dispatch:
//   an UNPLACED CONFIG — an object carrying a model's vocabulary under a key this file does not read
//     (`forEachUnplacedConfig`); the measured case is the slimmable wrapper's `config.model`;
//   a REFUSED VALUE — a count that is there and cannot be read (`readCountOrNote`);
//   a DISCARDED READING — a stack was read, and a declared or top-level field beside it was set aside
//     (the junction in `placedField`).
// The answer is then `placed + kUnreadShapeCeiling`, charged ONCE for the whole tree.
//
// WHY ONCE, AND WHY ADDED RATHER THAN MAXED — each is the answer to a measured failure of the other.
//   A floor at EVERY node that carried an unplaced child made the allowance ADDITIVE across siblings:
//   a 30-byte dead `{"layers":[],"m":{"layers":[]}}` cost a full 48 000 each, so 50 000 of them — a
//   1.6 MB file NAM loads in microseconds — drained INT_MAX, which is about half an hour per lane of a
//   real Standard inside `reset()`. Ignorance of N things is one ignorance, not N.
//   A floor as a MAX at the root lets a large known part SWALLOW the allowance: a read stack of
//   100 000 in series with an unreadable stage answers 100 000, i.e. the unknown stage is charged
//   ZERO — the one answer the rule forbids. Added, it is 148 000.
//
// WHY NOT FACE VALUE. An unplaced or set-aside number may be DEAD — NAM never builds it — and a dead
// number is not bounded by anything: `"receptive_field": 2147483647` beside a real Standard's stack is
// a file NAM loads unchanged, and trusting it made `reset()` spend about half an hour per lane. NAM
// allocates a LIVE stack of that length and refuses the load; a dead one costs NAM nothing. This file
// cannot tell the two apart without restating NAM's dispatch (rule 9u), so it trusts neither.
//
// 🔴 AND THE RULE HAS TWO DOORS, BOTH MEASURED, AND NEITHER CAN BE SHUT WITHOUT OPENING THE OTHER.
// (A differential fuzz of 2 000 loadable configs against NAM found no short answer outside these and
// the recurrent exception; these are what it found.)
//   D1 — a LIVE memory this file cannot place is charged the allowance, so one LONGER than the allowance
//     drains short by the difference: a wrapped model whose inner field is 60 000 drains 50 048 (the
//     row that measures it leaked 9 952 samples); a 60 001-tap `Linear` carrying a dead readable
//     `layers` array has its declared field SET ASIDE for the stack and drains 50 050. The longest field
//     of any real capture is 6347.
//   D2 — a DEAD number this file PLACES is trusted at face value, as it was before this rule: a lower
//     reading with no stack beside it (the two lower readings are maxed, as P87 ratified — a dead
//     `dilations:[2e9]` beside a `Linear`'s declared field), the wrapped form's own decoy stack (its
//     top-level `layers` are placed, and NAM builds from `config.model` instead — a decoy spelling
//     `dilations:[2000000]` drains two million samples), a `layers` array on an architecture that never
//     reads it. Base answered the same numbers for all of them.
// Treating EVERY uncertain config as unplaceable (dropping its own reading too) shuts D2 and widens D1
// to a live stack that merely has a dead sibling; trusting face value shuts D1 and reopens the
// half-hour `reset()` on 30 dead bytes. Which door stays open is a policy decision, and it is
// registered as one.
//
// THE NUMBER. 48 000 samples, a POLICY constant: 7.6x the longest field any real capture has
// (A2.nam's 6347), one second at 48 kHz and a quarter of one at 192 kHz — which in TIME is still 1.9x
// A2's 132 ms, because a capture built to cover the same milliseconds at a higher rate needs
// proportionally more samples. Against NAM's own half-second heuristic for an unbounded memory it is
// 2x at 48 kHz, 1x at 96 kHz and HALF at 192 kHz — said plainly, because a rounded comparison in a
// comment is this project's documented systemic leak.
//
// THE PRICE, MEASURED, because it is SPENT and not displayed. `configureRates` sizes `drainSamples_`
// from this number and `reset()` drives that many samples of silence through the network — and
// `reset()` is documented callable from the AUDIO thread, so this lands there, once per departure, per
// lane. Timed through `NamStage::reset()` itself on NAM's shipped captures, each made slimmable at full
// width and then WRAPPED, so the ledger charges 0 + 48 000 + the 2048 ring = 50 048 samples; M-series
// core, 48 kHz, one mono restart, median of seven:
//     capture                    as shipped             wrapped, block 256    block 64    block 1024
//     `wavenet_a1_standard.nam`  4093 samples  3.44 ms  39.7 ms               44.5 ms     37.5 ms
//     `A2.nam` (its submodel)    6347 samples  3.06 ms  24.1 ms
//     `slimmable_wavenet.nam`    2047 samples  0.19 ms   4.5 ms
// So a real capture arriving in a shape this file cannot place pays about 40 ms per lane plus what it
// placed — some 80 ms for a stereo restart, fifteen 256-sample callbacks at 48 kHz. That is a dropout,
// and it is the price of the rule; the alternative was 2046 samples of the previous sound.
//
// WHY NOT RATE-AWARE. Two of NAM's nine shipped captures — `my_model.nam` and
// `wavenet_a1_standard.nam` — carry NO `sample_rate` key at all (their top-level keys are
// `architecture, config, version, weights`), so a rate-aware allowance needs an absolute constant for
// the missing case anyway, which is two numbers where one does. And this file runs at
// `NamStage.cpp:831`, BEFORE `acceptsModelRate` and before `modelRunSR` exists.
constexpr int kUnreadShapeCeiling = 48000;

// 🔴 HOW MANY UNPLACED CONFIGS DEEP THE RECURRENCE WALK FOLLOWS ONE INTO ANOTHER — and ONLY those hops
// are counted. It is a crash guard and it is measured. The field never descends into an unplaced node
// (it charges the allowance instead), so the one NEW recursion this file has is `isRecurrent` /
// `partitionedTailSamples` looking through unplaced nodes for an architecture; its frame is about
// 500 B per hop at -O2 and 1.0-1.2 KB at -O0 (Apple clang), so an unguarded walk died at ~1000 hops on
// a 512 KiB thread — a plain `std::thread`'s stack on macOS — and 32 hops is about 16 KB (-O2) / 37 KB
// (-O0): 14x under that stack even unoptimised. NAM itself unwraps `config.model` exactly ONCE
// (`wavenet/slimmable.cpp:544` is not recursive), and the deepest nesting of any kind in the author's
// corpus of 1229 distinct captures is ONE, so 32 hops is 32x anything that exists.
// Past it the walk stops: the answer carries the allowance already (the first unplaced node set it),
// and recurrence is never assumed (see `isRecurrent`).
// 🔴 THE AXES BASE WALKED — `condition_dsp` and `submodels` — ARE NOT COUNTED, deliberately. A guard
// there TRUNCATED a model NAM loads: 65 plain WaveNets chained through `condition_dsp` into a
// 60 001-tap `Linear` drained 50 243 samples under a guard where base drained 62 243, and a restarted
// lane replayed 1.95484 out of digital silence. Those axes recurse exactly as base did, and so they
// carry base's exposure, which this file does not close: a 100 000-level chain overflows this walk on a
// small stack. In the load path NAM gets there first — `get_dsp` copies the config RECURSIVELY before
// this file runs (`get_dsp.cpp:146`), and on a 512 KiB thread that copy takes the host down at some
// 2 000 levels (134 KB of JSON) through `NamStage::prepareModel`. Registered, not closed here.
constexpr int kMaxUnplacedHops = 32;

// THE ПРИЗНАК: an object that DESCRIBES A MODEL, judged only by the vocabulary THIS file answers in.
// It deliberately knows nothing about which parser NAM would hand the node to — that heuristic is
// upstream's, a second copy of it would have to track upstream forever, and rule 9u exists because this
// project has paid for such copies three times. All this says is: *here is a description of memory, and
// I cannot place it in the topology I understand.*
// `model` is in the list because it is NAM's OWN name for a nested node in both places it nests one
// (`container.cpp:160`, `wavenet/slimmable.cpp:544`); `config` catches a whole model node. Nothing
// here is keyed on the slimmable marker or on `slice_channels_uniform` — rename the wrapper's key and
// the признак still fires on its shape.
inline bool looksLikeModelDescription (const nlohmann::json& v)
{
    if (! v.is_object()) return false;
    return (v.contains ("config")        && v["config"].is_object())
        || (v.contains ("model")         && v["model"].is_object())
        || (v.contains ("layers")        && v["layers"].is_array())
        || (v.contains ("dilations")     && v["dilations"].is_array())
        || (v.contains ("submodels")     && v["submodels"].is_array())
        || (v.contains ("condition_dsp") && v["condition_dsp"].is_object())
        ||  v.contains ("receptive_field");
}

// The keys this file ACCOUNTS FOR at each position it walks. Everything else at those positions is a
// place a model could hang that nothing here reads. Stated as what IS read, because the other list has
// no end.
inline bool isReadConfigKey (const std::string& k)
{
    return k == "layers" || k == "dilations" || k == "head" || k == "receptive_field"
        || k == "condition_dsp" || k == "submodels";
}

inline bool isReadLayerKey (const std::string& k)
{
    return k == "dilations" || k == "kernel_sizes" || k == "kernel_size" || k == "head";
}

inline bool isReadSubmodelKey (const std::string& k)
{
    return k == "model";
}

namespace unplaced
{
// One key at one position. An object under an unread key is visited when it carries the vocabulary —
// or when the key is `config`, because a model node's config is its config by the FILE FORMAT this file
// already reads everywhere (`model["config"]`), whatever it happens to contain: an LSTM's config has no
// vocabulary word at all, and a node hiding one under a vocabulary-less config was walked past. An
// ARRAY under an unread key is the `submodels` shape by another name, so its elements are asked the
// same question.
inline void visitKey (const std::string& key, const nlohmann::json& value, bool (*isRead) (const std::string&),
                      void (*visit) (const nlohmann::json&, void*), void* ctx)
{
    if (isRead (key)) return;
    if (looksLikeModelDescription (value) || (key == "config" && value.is_object())) { visit (value, ctx); return; }
    if (value.is_array())
        for (const auto& entry : value)
            if (looksLikeModelDescription (entry)) visit (entry, ctx);
}

inline void visitEntries (const nlohmann::json& cfg, const char* arrayKey, bool (*isRead) (const std::string&),
                          void (*visit) (const nlohmann::json&, void*), void* ctx)
{
    if (! cfg.contains (arrayKey) || ! cfg[arrayKey].is_array()) return;
    for (const auto& entry : cfg[arrayKey])
    {
        if (! entry.is_object()) continue;     // a non-object entry has no keys; see the row in the suite
        for (const auto& item : entry.items())
            visitKey (item.key(), item.value(), isRead, visit, ctx);
    }
}
} // namespace unplaced

// 🔴 THE ONE WALK OVER THE CONFIGS THIS FILE CANNOT PLACE — and it being ONE is the point, not tidiness.
// The field's allowance, the ring and the recurrence all reach the unplaced nodes through this function,
// so they cannot disagree about WHICH nodes exist. Three functions answering differently about one model
// is the class this whole line of work exists to close, and the wrapper was a measured instance of it:
// field blind, ring blind, recurrence blind, all at once.
//
// WHERE IT LOOKS. At a config's own keys, and at the keys of every entry of the two arrays this file
// walks (`layers`, `submodels`). A description of memory hanging off any of those is one this file is
// not reading. It does NOT descend into arbitrary sub-objects: `metadata` — which sits on the MODEL
// node, a sibling of `config`, and is never visited from a model this file places — carries a real
// capture's free-form user tree (`training.data.latency.calibration...`), and a hunt through unschema'd
// JSON for words that look like a stack would invent memory for data nothing builds. `config` has a
// schema and `metadata` has not; that is the line, and it is structural.
// Two review seats argued against scanning the ENTRIES, on the grounds that a model there would be a
// new NAM dispatch and therefore a number in the plan; they are scanned anyway, because that argument
// was equally true of `config.model` until it was measured, and a false fire costs the allowance, which
// is the side the rule has already chosen.
//
// A config that is not an object is not iterated. NAM cannot load one (its parsers index the config by
// key and throw), so the answer is never spent — and an ARRAY config is the one shape where iterating
// would find something: its elements would read as unplaced configs.
//
// THE COST TODAY IS ZERO AND IT IS MEASURED: over the 1229 distinct captures on the author's machine —
// WaveNets, LSTMs and SlimmableContainers, no Linear and no ConvNet among them — this walk visits
// nothing at all. The only object-valued config key that occurs is `condition_dsp`, which is read; every
// object in a layer entry is a layer FEATURE (`head`, `head1x1`, `layer1x1`, eight `*_film` objects,
// `activation`, `slimmable`) carrying no vocabulary word.
inline void forEachUnplacedConfig (const nlohmann::json& cfg,
                                   void (*visit) (const nlohmann::json&, void*), void* ctx)
{
    if (! cfg.is_object()) return;
    for (const auto& item : cfg.items())
        unplaced::visitKey (item.key(), item.value(), isReadConfigKey, visit, ctx);
    unplaced::visitEntries (cfg, "layers",    isReadLayerKey,    visit, ctx);
    unplaced::visitEntries (cfg, "submodels", isReadSubmodelKey, visit, ctx);
}
//======================================================================================================

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
inline bool anyNestedInConfig (const nlohmann::json& cfg, bool (*pred) (const nlohmann::json&), int hops);
inline bool carriesAllowance (const nlohmann::json& model);

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
    // 🔴 …AND IT IS CHARGED WHEREVER THE FIELD CARRIES THE ALLOWANCE, taken from the field's OWN walk
    // (`carriesAllowance`) rather than from a second one that tries to agree with it. The ring is ADDED
    // to the field (`NamStage.cpp:543`, `field + drainTail_`), so the allowance does not cover it, and a
    // node this file could not place may be a `Linear` whose architecture nobody can see. 2048 samples is
    // under 2 ms of a real Standard's inference (1.6 ms at a 256 block, 1.8 at 64); a first draft that
    // asked a separate scan here diverged from the field twice (a submodel's unplaced node; the far side
    // of a guard).
    return (isLinearArchitecture (model) || anyNestedModel (model, isLinearArchitecture)
            || carriesAllowance (model)) ? 2 * 1024 : 0;
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

// 🔴 IT IS NEVER ASSUMED — the one place the allowance rule stops — BUT IT IS LOOKED FOR THROUGH NODES
// THIS FILE COULD NOT PLACE. Recurrence is a claim about the KIND of state ("no finite silence empties
// it"), and the arithmetic says assuming it buys nothing anyway: `configureRates` spends
// `fmax(prewarm, 0.5 * modelRunSR)` for a recurrent capture (`NamStage.cpp:539`), and with the
// allowance in `prewarm` that is the allowance for every model rate up to 96 kHz. What the flag WOULD
// buy is permanent: a recurrent lane is re-charged the whole drain on every `reset()` and its
// `everFed_` is never cleared (`NamStage.cpp:369, :393`), so every later `prepare()` re-charges it too.
// But an `architecture:"LSTM"` FOUND under an unplaced key is evidence, not assumption, and the wrapper
// that made this rule carries its conditioner exactly there: measured, a slimmable wrapper with an LSTM
// conditioner never goes silent, and read as not recurrent it would be marked clean after one drain.
// The price is stated rather than hidden: an LSTM-shaped object under a key NAM never builds — dead
// JSON at this pin — makes the lane re-drain on every restart. And a recurrent model this file cannot
// see at all gets the allowance, which is 2x NAM's own heuristic at 48 kHz, and is then marked clean.
inline bool isRecurrent (const nlohmann::json& model)
{
    return isLstmArchitecture (model) || anyNestedModel (model, isLstmArchitecture);
}

// True when ANY model nested in this one satisfies `pred` — through the two nesting keys named at the
// top of this file, and through every config this file could not place, because a predicate that
// reached only part of the tree is the same defect in a different function.
inline bool anyNestedInModel (const nlohmann::json& model, bool (*pred) (const nlohmann::json&), int hops)
{
    return model.contains ("config") ? anyNestedInConfig (model["config"], pred, hops) : false;
}

inline bool anyNestedModel (const nlohmann::json& model, bool (*pred) (const nlohmann::json&))
{
    return anyNestedInModel (model, pred, 0);
}

namespace unplaced
{
struct PredScan { bool (*pred) (const nlohmann::json&); int hops; bool hit; };

inline void scanPred (const nlohmann::json& v, void* ctx)
{
    auto* s = static_cast<PredScan*> (ctx);
    if (s->hit || s->hops >= kMaxUnplacedHops) return;
    // 🔴 ONE READING OF AN UNPLACED NODE, AND IT COVERS BOTH THINGS THE NODE MIGHT BE. It may be a MODEL
    // NODE (`architecture` + `config`) or a RAW CONFIG — the wrapper's `config.model` is the latter
    // (`wavenet/slimmable.cpp:544` hands it straight to `parse_config_json`). Deciding which NAM would
    // build is restating NAM's dispatch, so it is walked as a raw config, full stop: if it is really a
    // model node, its `config` key is visited as an unplaced child one hop down (`visitKey`). `pred` is
    // applied directly first, because a raw config has no `architecture` for it to read.
    // 🔴 AND NOT TWO READINGS, which is what a first draft did, and it was a HANG. Reading the node as a
    // model (through `v["config"]`) AND as a raw config (which reaches `v["config"]` again as an
    // unplaced child) visits every level by two paths, and the work grew like Fibonacci: measured on
    // `{"layers":[1],"config":{"layers":[1],"config":…}}`, 1.07 ms at 16 levels, 7.93 ms at 20, 19.1 ms
    // at 22 — about 1.6x per level, which is of the order of 10^7 seconds at 64 levels, from a file of
    // a few kilobytes. One reading visits each node once.
    if (s->pred (v) || anyNestedInConfig (v, s->pred, s->hops + 1))
        s->hit = true;
}

struct AnyScan { bool hit; };
inline void scanAny (const nlohmann::json&, void* ctx) { static_cast<AnyScan*> (ctx)->hit = true; }
} // namespace unplaced

inline bool anyNestedInConfig (const nlohmann::json& cfg, bool (*pred) (const nlohmann::json&), int hops)
{
    if (const nlohmann::json* cond = conditionerOf (cfg))
        if (pred (*cond) || anyNestedInModel (*cond, pred, hops))
            return true;
    if (cfg.contains ("submodels") && cfg["submodels"].is_array())
        for (const auto& sub : cfg["submodels"])
            if (sub.contains ("model") && (pred (sub["model"]) || anyNestedInModel (sub["model"], pred, hops)))
                return true;
    unplaced::PredScan s { pred, hops, false };
    forEachUnplacedConfig (cfg, unplaced::scanPred, &s);
    return s.hit;
}

// The FIELD of what this file can place, over the two axes base walked — and whether anything in the
// tree could not be placed. It never descends into an unplaced node: that node is the allowance.
inline long long placedField (const nlohmann::json& cfg, bool& unknown)
{
    // A container holds several models and switches between them by level; any of them can be the
    // one speaking, so the longest memory is the one that has to be waited out.
    long long worst = 0;
    if (cfg.contains ("submodels") && cfg["submodels"].is_array())
        for (const auto& sub : cfg["submodels"])
            if (sub.contains ("model") && sub["model"].contains ("config"))
                worst = std::max (worst, placedField (sub["model"]["config"], unknown));
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
    // 🔴 THE STACK IS READ FIRST AND A DECLARATION BESIDE IT IS NOT TRUSTED — BUT NOT IGNORED EITHER.
    // A WaveNet config may carry `receptive_field` or a top-level `dilations` as well, and NAM's WaveNet
    // parser reads neither; a `Linear` or a `ConvNet` may carry `layers`, and theirs reads no `layers`.
    // The stack is the reading this file takes, and a lower reading beside it is DISCARDED — which was
    // the defect: a 5000-tap `Linear` carrying a readable stray `layers` array answered 2 for an impulse
    // reaching 4999, on a loaded model. Trusting the larger instead (a max) closed that and opened a
    // worse door: `"receptive_field": 2147483647` beside a real Standard's stack loads unchanged, and
    // `reset()` then ran for about half an hour per lane. So a discarded reading is neither: it is a thing
    // this file could not place, and it costs the allowance. The two LOWER readings are maxed against
    // each other as P87 ratified — a chain there let the ConvNet reader silence the declared one.
    bool refused = false;
    const long long stack = receptiveFieldOfLayers (cfg, &refused);
    const long long lower = std::max ((long long) convNetField (cfg, &refused),
                                      (long long) declaredReceptiveField (cfg, &refused));
    const long long own = stack > 0 ? stack : lower;
    if (stack > 0 && lower > 0) unknown = true;        // a reading set aside
    if (refused) unknown = true;                       // a value there that could not be read
    // …AND ANY CONFIG THIS FILE CANNOT PLACE, which it does not read at all: see the rule above.
    {
        unplaced::AnyScan u { false };
        forEachUnplacedConfig (cfg, unplaced::scanAny, &u);
        if (u.hit) unknown = true;
    }
    // …AND THE CONDITIONER IS ADDED TO IT, because it is in series (see the top of this file). The sums
    // are taken in long long and clamped: every addend is already an int's worth at most.
    long long condField = 0;
    if (const nlohmann::json* cond = conditionerOf (cfg))
        if (cond->contains ("config"))
            condField = placedField ((*cond)["config"], unknown);
    return std::max (worst, clampCount (own + condField));
}

inline int receptiveFieldFromConfig (const nlohmann::json& model)
{
    if (! model.contains ("config")) return 0;
    bool unknown = false;
    const long long placed = placedField (model["config"], unknown);
    return (int) clampCount (placed + (unknown ? (long long) kUnreadShapeCeiling : 0LL));
}

// Whether the field carries the allowance — asked OF THE FIELD'S OWN WALK, so the ring cannot disagree
// with it.
inline bool carriesAllowance (const nlohmann::json& model)
{
    if (! model.contains ("config")) return false;
    bool unknown = false;
    (void) placedField (model["config"], unknown);
    return unknown;
}
} // namespace felitronics::nam::detail

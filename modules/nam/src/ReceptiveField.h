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
// from the config itself when it is not (`wavenet/slimmable.cpp:543`). A config carrying both the
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
// cannot derive, answered with a number anyway — had two further members on `main`, both measured on
// loaded models with real weights:
//   a `Linear` of 5000 dense taps carrying a READABLE stray `layers` array answered 2 for an impulse
//   reaching 4999, because `own == 0` was a CHAIN and a non-zero stack reading suppressed the declared
//   field (the very shape the comment in `fieldOfConfig` used to condemn between the other two);
//   a `ConvNet` carrying a NON-EMPTY dead `layers` array answered 0 where NAM's field is 256, because
//   `convNetField`'s suppression had been fixed for the EMPTY array only.
// Both are closed here by one rule rather than two branches: THE THREE SOURCES ARE MAXED, NEVER
// CHAINED, and no reader may silence another. See `fieldOfConfig`.

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
    // 🔴 NO `layers` GATE, AND THE REASON IS THE MAX THAT CALLS THIS. A gate stood here — "a WaveNet
    // keeps its dilations inside `layers`, so a top-level `dilations` beside real layers is a key NAM
    // never looks at there" — and it was a READER SILENCING ANOTHER, the same mistake as a chain. P87
    // took it back for an EMPTY `layers` (a ConvNet carrying `"layers": []` drained 102 for a model
    // reaching 131, as a slimmable capture's conditioner); it still fired for a NON-EMPTY one, and NAM's
    // ConvNet parser reads `layers` at no length at all (`convnet.cpp:326-338`), so a ConvNet carrying
    // `"layers":[{}]` LOADS and this answered 0 where NAM's field is 256. The gate existed to keep a
    // WaveNet's stack from being counted twice; `receptiveFieldFromConfig` now takes the MAX of the
    // readers, which cannot count anything twice, so there is nothing left for a gate to protect.
    // What moves is a synthetic row — a WaveNet with a 2-sample stack and a stray top-level
    // `dilations:[1000]` now answers 1001 — and no real capture: across 2369 of them no config carries a
    // top-level `dilations`.
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

//======================================================================================================
// 🔴 THE ANSWER FOR A CONFIG THIS FILE CANNOT PLACE, AND WHAT IT COSTS.
//
// THE NUMBER. 48 000 samples. It is a POLICY constant and not a derivation, so it is named, measured
// and justified rather than explained: it is 7.6x the longest field any real capture has (A2.nam's
// 6347), one second at 48 kHz and a quarter of one at 192 kHz — which in TIME is still 1.9x A2's
// 132 ms, because a capture built to cover the same milliseconds at a higher rate needs proportionally
// more samples of field.
//
// THE PRICE, MEASURED, because it is SPENT and not displayed. `configureRates` sizes `drainSamples_`
// from this number and `reset()` drives that many samples of silence through the network — and
// `reset()` is documented callable from the AUDIO thread, so this lands there, once per departure, per
// lane. Timed through `NamStage::reset()` itself on NAM's shipped captures, each made slimmable at full
// width and then WRAPPED, so the ledger charges 48 000 + the 2048 ring = 50 048 samples; M-series core,
// 48 kHz, one mono restart, median of seven:
//     capture                    as shipped             wrapped, block 256    block 64    block 1024
//     `wavenet_a1_standard.nam`  4093 samples  3.44 ms  39.7 ms               44.5 ms     37.5 ms
//     `A2.nam` (its submodel)    6347 samples  3.06 ms  24.1 ms
//     `slimmable_wavenet.nam`    2047 samples  0.19 ms   4.5 ms
// So about 40 ms per lane is the worst a real capture would pay if it arrived in a shape this file
// cannot place — some 80 ms for a stereo restart, fifteen 256-sample callbacks at 48 kHz. That is a
// dropout, and it is the price of the rule; the alternative was 2046 samples of the previous sound.
// The comparison against NAM's own half-second heuristic for an unbounded memory holds only up to
// 96 kHz: 48 000 is 2x it at 48 kHz, 1x at 96 kHz and HALF of it at 192 kHz. Said plainly rather than
// rounded, because a stale comparison in a comment is this project's documented systemic leak.
//
// WHY NOT RATE-AWARE. Two of NAM's nine shipped captures — `my_model.nam` and
// `wavenet_a1_standard.nam` — carry NO `sample_rate` key at all (verified: their top-level keys are
// `architecture, config, version, weights`), so a rate-aware ceiling needs an absolute constant for the
// missing case anyway, which is two numbers where one does. And this file runs at
// `NamStage.cpp:831`, BEFORE `acceptsModelRate` and before `modelRunSR` exists.
// WHY NOT GEOMETRIC. Nothing in a config bounds a dilated architecture: nine scalars write
// `dilations:[2000000]` (see the cap note in `receptiveFieldOfLayers`). There is no geometric bound on
// a config you could not read.
constexpr int kUnreadShapeCeiling = 48000;

// 🔴 AND HOW DEEP THIS FILE WILL FOLLOW A CONFIG BEFORE IT CALLS THAT "CANNOT PLACE" TOO.
// This is a CRASH GUARD and it is measured. The walks below recurse once per nesting level, and a
// `.nam` is user data: nlohmann 3.12 parses and destroys ITERATIVELY (its parser drives an explicit
// `states` vector), so a 100 000-deep file arrives intact and only this file's recursion dies on it.
// Measured on the code as it stood, nesting through `condition_dsp` alone: 60 000 levels survive on
// the 8 MiB main stack and 80 000 SIGSEGVs; on a 512 KiB worker stack — which is what a plain
// `std::thread` gets, and `prepareModel` is "any thread but the audio thread" — 4 000 survive and
// 6 000 die, on about 205 KB of JSON. `prepareModel`'s `catch (...)` cannot catch a SIGSEGV and the
// 64 MiB unpack guard is three orders of magnitude above the trigger.
// That exposure predates this file's unplaced-config walk but is WIDENED by it, because an unplaced
// key name is chosen by whoever wrote the file rather than by NAM — and NAM itself unwraps `config.model`
// exactly ONCE (`wavenet/slimmable.cpp:543` is not recursive), so a chain of dead `model` keys is a
// file NAM loads happily and this file used to be asked to walk. 64 is 32x the deepest nesting any real
// capture has (2: a container's submodel that carries a conditioner) and about 7.7 KB of stack, three
// orders of magnitude below the smallest stack this ships on. Past it the answer is the ceiling, which
// is the same thing this file says about every other shape it cannot read.
constexpr int kMaxConfigNesting = 64;

// THE ПРИЗНАК: an object that DESCRIBES A MODEL, judged only by the vocabulary THIS file answers in.
// It deliberately knows nothing about which parser NAM would hand the node to — that heuristic is
// upstream's, a second copy of it would have to track upstream forever, and rule 9u exists because this
// project has paid for such copies three times. All this says is: *here is a description of memory, and
// I cannot place it in the topology I understand.*
//
// `model` is in the list because it is NAM's OWN name for a nested node in both places it nests one
// (`container.cpp:161`, `slimmable.cpp:543`); `config` catches a whole model node. Nothing here is
// keyed on the slimmable marker, on `slice_channels_uniform`, or on where NAM looks — remove the
// wrapper's key name from a config and the признак still fires on its shape.
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

// The keys this file ACCOUNTS FOR at a config, and at one of its layer-array entries. Everything else
// at those two positions is a place a model could hang that nothing here walks. Stated as what IS read
// rather than as what is not, because the second list has no end.
inline bool isReadConfigKey (const std::string& k)
{
    return k == "layers" || k == "dilations" || k == "head" || k == "receptive_field"
        || k == "condition_dsp" || k == "submodels";
}

inline bool isReadLayerKey (const std::string& k)
{
    return k == "dilations" || k == "kernel_sizes" || k == "kernel_size" || k == "head";
}

// 🔴 THE ONE WALK OVER THE CONFIGS THIS FILE CANNOT PLACE — and it being ONE is the point, not tidiness.
// `receptiveFieldFromConfig`, `isRecurrent` and `partitionedTailSamples` all reach the unplaced nodes
// through this function, so they cannot disagree about WHICH nodes exist. Three functions answering
// differently about one model is the class this whole line of work exists to close, and the hybrid was
// a measured instance of it: field blind, ring blind, recurrence blind, all at once.
//
// WHERE IT LOOKS. At a config's own keys, and at the keys of each entry in its `layers` array. Those
// are the two positions this file walks, and a description of memory hanging off either is one it is
// not reading. It does NOT descend into arbitrary sub-objects: `metadata` — which sits on the MODEL
// node, a sibling of `config`, and is never visited from here — carries a real capture's free-form
// user tree (`training.data.latency.calibration...`), and a scan that went hunting through unschema'd
// JSON for words that look like a stack would invent memory for data nothing builds. `config` has a
// schema and `metadata` has not; that is the line, and it is structural rather than a survey result.
//
// THE COST TODAY IS ZERO AND IT IS MEASURED, NOT ASSUMED: over 2369 real captures on this machine —
// 589 112 model nodes and 393 082 layer entries, NAM's own examples and namz's conformance vectors
// included — this walk visits nothing at all. The only object-valued config key that ever occurs is
// `condition_dsp`, which is read; every layer-entry object is a layer FEATURE (`head1x1`, `layer1x1`,
// the eight `*_film` objects, `activation`, `slimmable`) and not one of them carries the vocabulary.
// Scanning the layer entries was argued against twice in the round on the grounds that a model there
// would be a new NAM dispatch and therefore a number in the plan; it is scanned anyway, because that
// argument was equally true of `config.model` until it was measured, and because a false fire here is
// an over-charge, which is the side the rule has already chosen.
inline void forEachUnplacedConfig (const nlohmann::json& cfg,
                                   void (*visit) (const nlohmann::json&, void*), void* ctx)
{
    if (! cfg.is_object()) return;                       // `config: null / [] / 7` is not iterated
    for (const auto& item : cfg.items())
    {
        if (isReadConfigKey (item.key())) continue;
        if (looksLikeModelDescription (item.value())) { visit (item.value(), ctx); continue; }
        // …and an ARRAY of them, which is the `submodels` shape under a name this file does not know.
        if (item.value().is_array())
            for (const auto& entry : item.value())
                if (looksLikeModelDescription (entry)) visit (entry, ctx);
    }
    if (! cfg.contains ("layers") || ! cfg["layers"].is_array()) return;
    for (const auto& grp : cfg["layers"])
    {
        if (! grp.is_object()) continue;
        for (const auto& item : grp.items())
            if (! isReadLayerKey (item.key()) && looksLikeModelDescription (item.value()))
                visit (item.value(), ctx);
    }
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
inline bool anyNestedInConfig (const nlohmann::json& cfg, bool (*pred) (const nlohmann::json&), int depth);
inline bool anyUnplacedConfig (const nlohmann::json& model);

inline bool isLinearArchitecture (const nlohmann::json& model)
{
    return model.contains ("architecture") && model["architecture"].is_string()
        && model["architecture"].get<std::string>() == "Linear";
}

inline int partitionedTailSamples (const nlohmann::json& model)
{
    // 🔴 AND IT IS CHARGED FOR A SHAPE THIS FILE COULD NOT PLACE, unconditionally. The ring is ADDED
    // to the field (`NamStage.cpp:543`, `field + drainTail_`), so the ceiling does not cover it: a
    // wrapper hiding a `Linear` would take the ceiling for the field and still lose the ring. 2048
    // samples is 43 µs; the alternative is a lane that drains its whole ceiling and still emits.
    // It is also what keeps the THREE readers moving together on the признак, which is the property
    // the joint rows in the suite exist to hold.
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
    return (isLinearArchitecture (model) || anyNestedModel (model, isLinearArchitecture)
            || anyUnplacedConfig (model)) ? 2 * 1024 : 0;
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

// 🔴 AND IT IS DELIBERATELY *NOT* FORCED TRUE FOR A SHAPE THIS FILE COULD NOT PLACE, which is the one
// place the ceiling rule stops. Recurrence is a claim about the KIND of state — "no finite silence
// empties it" — and not about a magnitude, and the arithmetic says it buys nothing here anyway:
// `configureRates` spends `fmax(prewarm, 0.5 * modelRunSR)` for a recurrent capture
// (`NamStage.cpp:539`), and with `prewarm` at the 48 000 ceiling that is 48 000 for EVERY model rate
// up to 96 kHz — i.e. for every capture that exists. What the flag WOULD buy is permanent: a recurrent
// lane is re-charged the whole drain on every `reset()` and its `everFed_` is never cleared
// (`NamStage.cpp:369, :393`), so every later `prepare()` re-charges it too, forever. Zero drain length,
// unbounded repeated cost. So it stays a READ fact — and "read" still means read through the unplaced
// node: the walker below goes there, so an `architecture:"LSTM"` under a wrapper is found. What is left
// uncovered is a recurrent model whose config says so in a vocabulary this file has never seen; that
// one gets 48 000 samples, which is twice NAM's own heuristic at 48 kHz, and then is marked clean.
inline bool isRecurrent (const nlohmann::json& model)
{
    return isLstmArchitecture (model) || anyNestedModel (model, isLstmArchitecture);
}

// True when ANY model NESTED in this one, at any depth, satisfies `pred` — the nesting keys are the
// ones named at the top of this file, and ALL of them are walked here because a predicate that reached
// only part of the tree is the same defect in a different function. The container switches by level and
// any of them can be the one speaking; the conditioner is always running; a config this file could not
// place may hold either. Either way the answer has to be the worst case over the tree, which is the
// same rule receptiveFieldFromConfig applies to the field.
inline bool anyNestedInModel (const nlohmann::json& model, bool (*pred) (const nlohmann::json&), int depth);

inline bool anyNestedModel (const nlohmann::json& model, bool (*pred) (const nlohmann::json&))
{
    return anyNestedInModel (model, pred, 0);
}

namespace unplaced
{
// The visitor's context. A function pointer plus a `void*` rather than a template, so the walk stays
// ONE function with ONE definition of "unplaced" for every reader that uses it.
struct PredScan { bool (*pred) (const nlohmann::json&); int depth; bool hit; };

inline void scanPred (const nlohmann::json& v, void* ctx)
{
    auto* s = static_cast<PredScan*> (ctx);
    if (s->hit) return;
    // 🔴 ONE READING OF AN UNPLACED NODE, AND IT COVERS BOTH THINGS THE NODE MIGHT BE. An unplaced
    // object may be a MODEL NODE (`architecture` + `config`) or a RAW CONFIG — the wrapper's
    // `config.model` is the latter (`wavenet/slimmable.cpp:543` hands it straight to
    // `parse_config_json`, which reads `model_json["layers"]`). Deciding which NAM would build is
    // restating NAM's dispatch, so it is read as a raw config, full stop: if it is really a model node,
    // its own `config` key is not one this file reads at a config, so the walk reaches it as an
    // unplaced child one level down. `pred` is applied directly first, because a raw config has no
    // `architecture` for it to read and a model node does.
    // 🔴 AND NOT TWO READINGS, which is what this first said and it was a HANG. Reading the node both
    // as a model (through `v["config"]`) and as a raw config (which reaches `v["config"]` again as an
    // unplaced child) visits every level twice by two paths, and the work grows like Fibonacci: measured
    // on `{"layers":[1],"config":{"layers":[1],"config":…}}`, 1.07 ms at 16 levels, 7.93 ms at 20,
    // 19.1 ms at 22 — about 1.6x per level, which is some 10^6 seconds at the 64-level guard, from a
    // file of 1.5 KB. One reading visits each node once.
    if (s->pred (v) || anyNestedInConfig (v, s->pred, s->depth + 1))
        s->hit = true;
}
} // namespace unplaced

inline bool anyNestedInConfig (const nlohmann::json& cfg, bool (*pred) (const nlohmann::json&), int depth)
{
    if (depth > kMaxConfigNesting) return false;   // past the guard nothing is READ; the FIELD answers the ceiling
    if (const nlohmann::json* cond = conditionerOf (cfg))
        if (pred (*cond) || anyNestedInModel (*cond, pred, depth + 1))
            return true;
    if (cfg.contains ("submodels") && cfg["submodels"].is_array())
        for (const auto& sub : cfg["submodels"])
            if (sub.contains ("model") && (pred (sub["model"]) || anyNestedInModel (sub["model"], pred, depth + 1)))
                return true;
    unplaced::PredScan s { pred, depth, false };
    forEachUnplacedConfig (cfg, unplaced::scanPred, &s);
    return s.hit;
}

inline bool anyNestedInModel (const nlohmann::json& model, bool (*pred) (const nlohmann::json&), int depth)
{
    return model.contains ("config") ? anyNestedInConfig (model["config"], pred, depth) : false;
}

inline int fieldOfConfig (const nlohmann::json& cfg, int depth, bool& ceilinged);

inline int receptiveFieldFromConfig (const nlohmann::json& model)
{
    bool ceilinged = false;
    return model.contains ("config") ? fieldOfConfig (model["config"], 0, ceilinged) : 0;
}

// 🔴 WHETHER THE FIELD TOOK THE CEILING ANYWHERE IN THE TREE — asked OF THE FIELD'S OWN WALK, not of a
// second walk that tries to agree with it. The first version of this was a separate scan of the top
// config only, and it diverged on both axes a separate scan can: a container whose SUBMODEL carried an
// unplaced node read 48 000 for the field and 0 for the ring, and past the nesting guard the field said
// "I stopped reading" while the ring, having stopped reading too, said "I read, and there is no Linear"
// — measured on a `condition_dsp` chain: ring 2048 at 65 levels, 0 at 66. Deriving the answer from the
// function that decides the ceiling makes that divergence unwritable rather than untested.
inline bool anyUnplacedConfig (const nlohmann::json& model)
{
    bool ceilinged = false;
    if (model.contains ("config")) (void) fieldOfConfig (model["config"], 0, ceilinged);
    return ceilinged;
}

namespace unplaced
{
// Tallies what the configs this file could not place are worth, and remembers that there WAS one.
struct FieldScan { long long total; int depth; bool* ceilinged; bool local; };

inline void scanField (const nlohmann::json& v, void* ctx)
{
    auto* s = static_cast<FieldScan*> (ctx);
    s->local = true;                 // THIS node carries one: floor THIS node's answer
    *s->ceilinged = true;            // …and the tree took the ceiling somewhere: the ring reads this
    // ONE reading, for the reason and with the measured hang stated at `scanPred`: a model node's
    // `config` is reached as an unplaced child of the node read as a raw config.
    // ACROSS nodes, and against what the top level read, the composition is a SUM. Two seats of the
    // round argued for a max on the grounds that NAM builds exactly one of the two halves at this pin
    // (`wavenet/slimmable.cpp:543` picks one) — which is true, and is precisely the upstream fact rule
    // 9u forbids depending on. Not knowing whether an unplaced node is an ALTERNATIVE to what was read
    // or a STAGE IN SERIES with it, this takes series, because series is the worse of the two. It costs
    // nothing on the shape that made the rule: the wrapper's decoy top level reads 0.
    s->total = clampCount (s->total + (long long) fieldOfConfig (v, s->depth + 1, *s->ceilinged));
}
} // namespace unplaced

inline int fieldOfConfig (const nlohmann::json& cfg, int depth, bool& ceilinged)
{
    // 🔴 PAST THE NESTING GUARD THE ANSWER IS THE CEILING, not zero and not a crash — see
    // `kMaxConfigNesting`. "I stopped reading" is a case of "I do not know", and the rule for that is
    // already written. This also changes a RELEASED answer for a config nested deeper than 64 through
    // `condition_dsp` or `submodels`: it used to read its real (small) number and now reads 48 000.
    // The deepest nesting any real capture has is 2, and the same shape at 6 000 deep used to take the
    // process down on a 512 KiB thread, so the shift is from a crash-or-a-number to a bound.
    if (depth > kMaxConfigNesting) { ceilinged = true; return kUnreadShapeCeiling; }
    // A container holds several models and switches between them by level; any of them can be the
    // one speaking, so the longest memory is the one that has to be waited out.
    long long worst = 0;
    if (cfg.contains ("submodels") && cfg["submodels"].is_array())
        for (const auto& sub : cfg["submodels"])
            if (sub.contains ("model"))
                worst = std::max (worst, (long long) (sub["model"].contains ("config")
                                                          ? fieldOfConfig (sub["model"]["config"], depth + 1, ceilinged) : 0));
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
    // 🔴 THE THREE SOURCES ARE MAXED, NEVER CHAINED, AND NO READER MAY SILENCE ANOTHER.
    // This line used to read `own = layers; if (own == 0) own = max(convNet, declared);` — a chain,
    // and the comment that stood here condemned exactly that pattern between the other two while
    // leaving it in place ABOVE them. It cost two measured leaks on loaded models with real weights:
    //   a `Linear` of 5000 dense taps carrying a READABLE stray `layers` array (its parser reads
    //   neither `layers` nor `dilations`, `linear.cpp:306-316`, so the file LOADS) answered 2 where
    //   the impulse reaches 4999 — a drain of 2 + 2048 against a model that needs 4999, i.e. 2949
    //   samples of the previous sound handed back out of digital silence;
    //   a `Linear` carrying a stray `dilations` array answered 2 for the same reason before P87 split
    //   the lower two — the row is still in the suite, and this closes the half above it.
    // The chain's stated justification was that "maxing would drain a real capture for a stale number":
    // under the ratified asymmetry that is the CHEAP side, and it is also empty in fact — across 2369
    // real captures on this machine (589 112 model nodes) no config carries `receptive_field` or a
    // top-level `dilations` at all, so no shipped capture can pay for the max. What moves instead is
    // a synthetic row: a WaveNet carrying a stale `receptive_field` of 99999 beside a 9-sample stack
    // now answers 99998 rather than 9.
    const long long own = std::max ((long long) receptiveFieldOfLayers (cfg),
                                    std::max ((long long) convNetField (cfg),
                                              (long long) declaredReceptiveField (cfg)));
    // …AND THE CONDITIONER IS ADDED TO IT, because it is in series (see the top of this file). The sum
    // is taken in long long: both addends are already individually clamped to INT_MAX, so an int
    // addition here is undefined behaviour for a config a unit test can write down.
    const nlohmann::json* cond = conditionerOf (cfg);
    const long long condField = cond != nullptr && cond->contains ("config")
                              ? (long long) fieldOfConfig ((*cond)["config"], depth + 1, ceilinged) : 0LL;
    // …AND THE CONFIGS THIS FILE COULD NOT PLACE. See `forEachUnplacedConfig`. The ceiling is a FLOOR
    // on the whole answer and not an alternative to it: a node that reads as something is charged that
    // something AND at least the ceiling, so the branch can never be short, and there is no
    // "was it readable?" predicate to get wrong — `receptiveFieldOfLayers` returns 0 both for "there is
    // nothing here" and for "there are layers whose dilations I could not read", so that predicate does
    // not exist to be asked. The price is measured beside the constant: on the wrapped
    // `slimmable_wavenet.nam` that made the rule, the ceiling and its ring cost 4.5 ms of one lane per departure.
    unplaced::FieldScan u { 0, depth, &ceilinged, false };
    forEachUnplacedConfig (cfg, unplaced::scanField, &u);
    long long total = std::max (worst, clampCount (own + condField + u.total));
    if (u.local) total = std::max (total, (long long) kUnreadShapeCeiling);
    return (int) std::min (total, (long long) std::numeric_limits<int>::max());
}
} // namespace felitronics::nam::detail

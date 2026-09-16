// SPDX-License-Identifier: MIT
// How long a model has to be fed before its output means anything. The number matters because
// anything that fades a freshly loaded model in — a crossfade between two captures of one device —
// hears the difference: a network still describing the silence it was born into is not the device.
//
// These read the config's SHAPE only, never its weights, which is what lets a real architecture be
// checked here without shipping a real model.

#include "../src/ReceptiveField.h"

#include <felitronics_test.h>

#include <atomic>
#include <chrono>
#include <string>
#include <thread>

using felitronics::test::ok;
using felitronics::test::group;
using felitronics::nam::detail::receptiveFieldFromConfig;
using felitronics::nam::detail::receptiveFieldOfLayers;
using felitronics::nam::detail::partitionedTailSamples;
using felitronics::nam::detail::isRecurrent;
using felitronics::nam::detail::convNetField;
using felitronics::nam::detail::anyNestedModel;
using felitronics::nam::detail::kUnreadShapeCeiling;
using felitronics::nam::detail::kMaxUnplacedHops;

namespace {

// One dilated stack, written the way a .nam does.
nlohmann::json wavenet(std::vector<int> kernels, std::vector<int> dilations) {
    return { { "architecture", "WaveNet" },
             { "config", { { "layers", nlohmann::json::array({
                   { { "kernel_sizes", kernels }, { "dilations", dilations } } }) } } } };
}

// The same stack with a CONDITIONER — a whole model of its own, which is how NAM spells it: the node
// is handed straight back to `get_dsp` (v0.5.4 wavenet/model.cpp:844), so it has a model's shape.
nlohmann::json conditionedBy(const nlohmann::json& conditioner, std::vector<int> kernels,
                             std::vector<int> dilations) {
    auto model = wavenet(std::move(kernels), std::move(dilations));
    model["config"]["condition_dsp"] = conditioner;
    return model;
}

nlohmann::json linear(int receptiveField) {
    return { { "architecture", "Linear" }, { "config", { { "receptive_field", receptiveField } } } };
}

// NAM's own shipped `example_models/slimmable_wavenet.nam`, config verbatim (weights are not read here).
// Its field by hand, NOT by the function under test: kernel 3 charges 2 per unit of dilation, the
// dilations sum to 1023, so 2·1023 = 2046 samples of memory, and the file's usual +1 of margin: 2047.
nlohmann::json realSlimmableConfig() {
    return { { "layers", nlohmann::json::array({
                 { { "input_size", 1 }, { "condition_size", 1 }, { "head_size", 1 }, { "channels", 3 },
                   { "kernel_size", 3 },
                   { "dilations", nlohmann::json::array({ 1, 2, 4, 8, 16, 32, 64, 128, 256, 512 }) },
                   { "activation", "ReLU" }, { "gated", false }, { "head_bias", false },
                   { "slimmable", { { "method", "slice_channels_uniform" },
                                    { "kwargs", { { "allowed_channels", nlohmann::json::array({ 1, 2, 3 }) } } } } } } }) },
             { "head", nullptr }, { "head_scale", 0.02 } };
}

nlohmann::json flatSlimmable(const nlohmann::json& cfg) {
    return { { "version", "0.7.0" }, { "architecture", "WaveNet" }, { "config", cfg }, { "sample_rate", 48000 } };
}

// THE WRAPPER, the shape that LOADS and that this file read as zero: a decoy top-level `layers` carrying
// only the marker NAM's WaveNet parser delegates on, and the real config one level down under `key`.
// `key` is a parameter on purpose — "model" is where NAM looks today, and the rule must not care.
nlohmann::json wrapped(const nlohmann::json& inner, const std::string& key = "model") {
    nlohmann::json cfg = { { "layers", nlohmann::json::array({
                               { { "slimmable", { { "method", "slice_channels_uniform" } } } } }) } };
    cfg[key] = inner;
    return flatSlimmable(cfg);
}

// Counts how often a walk asks its predicate — the deterministic witness for "linear, not exponential".
std::atomic<long long> predCalls { 0 };
bool countingNever(const nlohmann::json&) { predCalls.fetch_add(1); return false; }

// Nesting chains, built as TEXT and parsed, because that is how a hostile file arrives: nlohmann parses
// and destroys iteratively, so only the reader's recursion is under test here.
// The axis P92 ADDED: dead `model` keys, which NAM unwraps exactly once and then ignores.
std::string wrapperChain(int depth) {
    std::string s = R"({"architecture":"WaveNet","config":)";
    for (int i = 0; i < depth; ++i) s += R"({"layers":[{"kernel_size":2,"dilations":[1]}],"model":)";
    s += R"({"layers":[{"kernel_size":2,"dilations":[1]}]})";
    for (int i = 0; i < depth; ++i) s += "}";
    return s + "}";
}
// An LSTM model node `hops` unplaced configs down: {"layers":[], "m": {"layers":[], "m": … LSTM}}.
std::string lstmUnderUnplaced(int hops) {
    std::string s = R"({"architecture":"WaveNet","config":{"layers":[])";
    for (int i = 1; i < hops; ++i) s += R"(,"m":{"layers":[])";
    s += R"(,"m":{"architecture":"LSTM","config":{"hidden_size":3}})";
    for (int i = 1; i < hops; ++i) s += "}";
    return s + "}}";
}
// A plain chain of WaveNets through `condition_dsp`, ending in a Linear of `taps` — a model NAM loads.
std::string conditionersInto(int depth, const std::string& leaf) {
    std::string s;
    for (int i = 0; i < depth; ++i) s += R"({"architecture":"WaveNet","config":{"layers":[1],"condition_dsp":)";
    s += leaf;
    for (int i = 0; i < depth; ++i) s += "}}";
    return s;
}
// The shape that made the first draft of the walk EXPONENTIAL: every level a node that is both a raw
// config and a model node, so a walk that read it both ways visited each level by two paths.
std::string doubleReadingChain(int depth) {
    std::string s = R"({"architecture":"WaveNet","config":{"layers":[1],"model":)";
    for (int i = 0; i < depth; ++i) s += R"({"layers":[1],"config":)";
    s += R"({"layers":[1]})";
    for (int i = 0; i < depth; ++i) s += "}";
    return s + "}}";
}

} // namespace

int main() {
    std::printf("felitronics::nam receptive-field tests\n");

    group("a dilated stack reaches back the sum of its taps");
    {
        // One layer, kernel 2, dilation 1 sees this sample and the one before it.
        ok(receptiveFieldFromConfig(wavenet({ 2 }, { 1 })) == 2, "kernel 2, dilation 1 -> 2 samples");
        // …and the LEGACY spelling, a single `kernel_size` for every layer. NAM takes either
        // (v0.5.4 wavenet/model.cpp:914-945) and refuses a config carrying both; reading only the array
        // answered ZERO for a legacy capture, which on a SlimmableWavenet — whose own answer is also
        // zero — meant no drain at all: measured 0.462117 out of digital silence.
        const auto legacy = [](int kernel, std::vector<int> dilations) {
            return nlohmann::json { { "architecture", "WaveNet" },
                                    { "config", { { "layers", nlohmann::json::array({
                                          { { "kernel_size", kernel }, { "dilations", dilations } } }) } } } };
        };
        ok(receptiveFieldFromConfig(legacy(2, { 1 })) == 2, "…the legacy single kernel_size reads the same");
        ok(receptiveFieldFromConfig(legacy(2, { 1, 3, 7 })) == 12, "…and applies to EVERY layer: 1+3+7 taps back");
        ok(receptiveFieldFromConfig(legacy(3, { 4 })) == 9, "…kernel 3, dilation 4 -> 4*2 + 1, as the array form");
        ok(receptiveFieldFromConfig(wavenet({ 3 }, { 4 })) == 9, "kernel 3, dilation 4 -> 4*2 + 1");
        ok(receptiveFieldFromConfig(wavenet({ 2, 2 }, { 1, 2 })) == 4, "…and layers add: 1 + 2 + 1");
    }

    group("the real thing: a SlimmableContainer of WaveNets");
    {
        // The architecture every capture in this project uses, config shape verbatim. Two submodels,
        // each the same stack; the answer is 6332 samples — 132 ms at 48 kHz, which is two and a half
        // times the fifty milliseconds a player might have guessed.
        const std::vector<int> ks { 6,6,6,6,6,6,6, 6,6,6,6,6,6,6, 15,15, 6,6,6,6,6,6,6 };
        const std::vector<int> ds { 1,3,7,17,41,101,239, 1,3,7,17,41,101,239, 1,13, 1,3,7,17,41,101,239 };
        const auto one = wavenet(ks, ds);
        ok(receptiveFieldFromConfig(one) == 6332, "one WaveNet stack reaches back 6332 samples");

        nlohmann::json container = { { "architecture", "SlimmableContainer" },
                                     { "config", { { "submodels", nlohmann::json::array({
                                           { { "max_value", 0.5 }, { "model", one } },
                                           { { "max_value", 1.0 }, { "model", one } } }) } } } };
        ok(receptiveFieldFromConfig(container) == 6332, "…and a container answers for its submodels");
    }

    group("a container answers with its LONGEST memory");
    {
        // Which submodel speaks depends on the signal's level, so any of them can be the one playing
        // when the fade starts. Waiting out the shortest would fade in a model that is still empty.
        nlohmann::json mixed = { { "architecture", "SlimmableContainer" },
                                 { "config", { { "submodels", nlohmann::json::array({
                                       { { "model", wavenet({ 2 }, { 1 }) } },
                                       { { "model", wavenet({ 3 }, { 100 }) } } }) } } } };
        ok(receptiveFieldFromConfig(mixed) == 201, "the deepest submodel sets the wait");
    }

    group("an architecture that DECLARES its field is read, not guessed at");
    {
        // This used to answer 0 with the note "the caller falls back to what NAM itself reports", and
        // NAM reports zero too: `Linear : Buffer : DSP` inherits `GetPrewarmSamples() { return 0; }`.
        // So the whole path answered 0 for an impulse response two thousand taps long — measured
        // through NamStage::prewarmSamples() at receptive_field 1 / 2 / 65 / 257 / 2001, all zero.
        // TAPS, NOT MEMORY: `y[n] = sum_{k<RF} h[k]x[n-k]`, so RF taps reach back RF-1 samples and a
        // ONE-tap capture is a gain with no memory at all. Reporting 1 for that is not harmless —
        // rigplayer::RigPlayer::warmFor() returns early on `pre <= 0`, so a 1 turns a memoryless
        // capture's warm-up into `1 + latency + maxBlock` and a woken slot goes silent for a block.
        ok(receptiveFieldFromConfig({ { "architecture", "Linear" }, { "config", { { "receptive_field", 1 } } } }) == 0,
           "a ONE-tap Linear capture is a gain: no memory, and the warm-up must stay at zero");
        ok(receptiveFieldFromConfig({ { "architecture", "Linear" }, { "config", { { "receptive_field", 3 } } } }) == 2,
           "a Linear capture's declared receptive_field is read — as TAPS, so the memory is one less");
        ok(receptiveFieldFromConfig({ { "architecture", "Linear" }, { "config", { { "receptive_field", 2001 } } } }) == 2000,
           "…at any length");
        // The number is SPENT now (NamStage feeds an absent lane silence for this long), so a cap on it
        // is a silent under-drain rather than a tidy display. Only the int conversion is guarded.
        ok(receptiveFieldFromConfig({ { "architecture", "Linear" }, { "config", { { "receptive_field", 2000000 } } } }) == 1999999,
           "…and it is NOT capped at 1<<20 = 1048576, which is where a spent number becomes a defect");
        ok(receptiveFieldOfLayers(wavenet({ 2 }, { 2000000 })["config"]) == 2000001,
           "…the same for a dilated stack that reaches back further than the old cap");
        // 🔴 A config may carry both, and the stack is the reading — BUT THE DECLARATION IT SETS ASIDE
        // COSTS THE ALLOWANCE (P92, moved on purpose; this row asserted 9). The stack describes what a
        // WaveNet does, and a declared number beside it may be dead; trusting it at face value (a max)
        // was measured to be a door: `receptive_field: 2147483647` beside a real Standard's stack loads
        // unchanged and made `reset()` run for about half an hour per lane. Ignoring it (the old chain) was
        // the other door, the row below. So a reading set aside is a thing this file could not place.
        nlohmann::json both = wavenet({ 3 }, { 4 });
        both["config"]["receptive_field"] = 99999;
        ok(receptiveFieldFromConfig(both) == 9 + 48000,
           "a declared number beside a stack is not trusted and not ignored: the stack's 9 plus the"
           " allowance — read " + std::to_string(receptiveFieldFromConfig(both)) + " (was 9 before P92)");
        both["config"]["receptive_field"] = 2147483647;
        ok(receptiveFieldFromConfig(both) == 9 + 48000,
           "…and a DEAD one at INT_MAX costs the same allowance, not 2^31 samples of audio-thread work");
        {
            nlohmann::json strayLayers = linear(5000);
            strayLayers["config"]["layers"] = nlohmann::json::array({
                { { "kernel_size", 2 }, { "dilations", nlohmann::json::array({ 1 }) } } });
            ok(receptiveFieldFromConfig(strayLayers) == 2 + 48000,
               "…which is what closes a Linear carrying a READABLE stray `layers` array: its parser reads"
               " neither key, the file loads, the chain answered 2 for a 4999-sample memory, and the"
               " declared field it sets aside now costs the allowance — read "
               + std::to_string(receptiveFieldFromConfig(strayLayers)));
        }
        ok(receptiveFieldFromConfig({ { "architecture", "Linear" }, { "config", { { "receptive_field", 0 } } } }) == 0,
           "a declared zero is zero");
        ok(partitionedTailSamples({ { "architecture", "Linear" } }) == 2048
               && partitionedTailSamples({ { "architecture", "WaveNet" } }) == 0
               && partitionedTailSamples(nlohmann::json::object()) == 0,
           "a Linear capture is charged its partitioned-FFT ring (2 x the largest block, NAM v0.5.4"
           " linear.cpp:14-17); a sample-by-sample architecture is not");
        ok(isRecurrent({ { "architecture", "LSTM" } })
               && ! isRecurrent({ { "architecture", "WaveNet" } })
               && ! isRecurrent(nlohmann::json::object()),
           "…and only LSTM is recurrent, where no finite length of silence empties the cell");
        // A CONTAINER IS ASKED THROUGH for both, the same way the field is: the top-level architecture
        // is SlimmableContainer, and answering for THAT charges no ring and no floor to a container of
        // Linear or LSTM submodels — whichever of them is speaking would then drain short.
        const auto container = [](const nlohmann::json& inner) {
            return nlohmann::json { { "architecture", "SlimmableContainer" },
                                    { "config", { { "submodels", nlohmann::json::array({
                                          { { "model", nlohmann::json { { "architecture", "WaveNet" } } } },
                                          { { "model", inner } } }) } } } };
        };
        ok(partitionedTailSamples(container({ { "architecture", "Linear" } })) == 2048,
           "a container holding a Linear submodel is charged the ring");
        ok(isRecurrent(container({ { "architecture", "LSTM" } })),
           "…and one holding an LSTM submodel is recurrent");
        ok(partitionedTailSamples(container({ { "architecture", "WaveNet" } })) == 0
               && ! isRecurrent(container({ { "architecture", "WaveNet" } })),
           "…and one holding neither is charged neither");
        ok(receptiveFieldFromConfig({ { "architecture", "Linear" }, { "config", { { "receptive_field", -5 } } } }) == 0,
           "…and a negative one is not carried into a length");
        ok(receptiveFieldFromConfig({ { "architecture", "Linear" }, { "config", { { "receptive_field", "long" } } } }) == 48000,
           "…and one that is not a number is refused rather than thrown on — this file runs inside"
           " prepareModel's catch-all, where a throw would REFUSE the load — and a value that is THERE"
           " and cannot be read is the allowance, not zero (P92; this row asserted 0)");
        // 🔴 BUT A FLOAT SPELLING IS A NUMBER. `2001.0` is not `is_number_integer()`, and NAM's own
        // parser takes it (`get<int>()` static_casts any arithmetic node), so a guard that demands an
        // integer refuses a capture that LOADS — and it drains for nothing. Measured: with the narrow
        // guard the dilated reader answered 0 where base answered 3.
        ok(receptiveFieldFromConfig({ { "architecture", "Linear" }, { "config", { { "receptive_field", 2001.0 } } } }) == 2000,
           "a field spelled 2001.0 is read, because that is a model NAM loads");
        ok(receptiveFieldOfLayers(wavenet({ 2 }, { 1 })["config"]) == 2, "…and the array form still reads");
        {
            nlohmann::json legacyFloat = { { "architecture", "WaveNet" },
                                           { "config", { { "layers", nlohmann::json::array({
                                                 { { "kernel_size", 2.0 },
                                                   { "dilations", nlohmann::json::array({ 1 }) } } }) } } } };
            ok(receptiveFieldFromConfig(legacyFloat) == 2,
               "…and the LEGACY single kernel_size spelled as a float, which is its own guard: narrowing"
               " it back to an integer passed every other row here");
        }
        {
            nlohmann::json floaty = { { "architecture", "WaveNet" },
                                      { "config", { { "layers", nlohmann::json::array({
                                            { { "kernel_sizes", nlohmann::json::array({ 2.0 }) },
                                              { "dilations",    nlohmann::json::array({ 2.0 }) } } }) } } } };
            ok(receptiveFieldFromConfig(floaty) == 3, "…and so does a dilated stack spelled in floats,"
                                                      " which is what the base did before this file guarded");
        }
    }

    group("an architecture nothing here can read says so, rather than guessing");
    {
        ok(receptiveFieldFromConfig({ { "architecture", "WaveNet" } }) == 0, "no config at all -> 0");
        ok(receptiveFieldFromConfig({ { "architecture", "WaveNet" }, { "config", { { "layers", 7 } } } }) == 0,
           "…and a layers field that is not an array is refused, not indexed");
        ok(receptiveFieldFromConfig({ { "architecture", "Linear" }, { "config", nlohmann::json::object() } }) == 0,
           "…and a config with neither layers nor a declared field is still zero");
    }


    group("a CONDITIONER is a whole model, and its memory is IN SERIES with the network's");
    {
        // 🔴 WHAT THIS CLOSES, and it is one hole with two readers. `config.condition_dsp` is a whole
        // model — NAM hands the node straight back to `get_dsp` (v0.5.4 wavenet/model.cpp:844) — and
        // nothing counted its memory: NAM answers ZERO for a `Linear` conditioner (the base class's
        // `GetPrewarmSamples`), and this file walked `submodels`, not `condition_dsp`. Law 11a's drain
        // and P47's stream restart both spend this number, so both were short by exactly the same
        // amount: measured through NamStage on a WaveNet with a 2001-sample Linear conditioner,
        // 0.905147969723 after a "full" drain against 0.905148267746 after a restart.
        //
        // SUM, NOT WORST-OF, and the difference is the number rather than the wording. The conditioner's
        // OUTPUT is the network's conditioning INPUT (`_process_condition` then the layer arrays against
        // its output, wavenet/model.cpp:699-729, :749-761), so the two memories compose in SERIES —
        // which is also how NAM's own arithmetic composes them (`mPrewarmSamples` starts at the
        // conditioner's prewarm and ADDS the stack, :616-620). Measured from outside this file, on a
        // loaded two-layer stack with a 2500-sample conditioner: the impulse's last non-zero sample is
        // 5000, where the same stack without a conditioner reaches 2501 and a worst-of would answer 2502.
        ok(receptiveFieldFromConfig(conditionedBy(linear(2002), { 2 }, { 1 })) == 2003,
           "the network's own 2 plus the conditioner's 2001 — a walk that DISCARDS the conditioner's"
           " field answers 2, and a worst-of answers 2001");
        ok(receptiveFieldFromConfig(conditionedBy(wavenet({ 2 }, { 100 }), { 2 }, { 1 })) == 103,
           "…and a WaveNet conditioner is read the same way: 2 + 101, each keeping its own sample of margin");
        // A conditioner may have a conditioner: the recursion is through receptiveFieldFromConfig, not
        // through a hand-rolled read of the conditioner's layers.
        ok(receptiveFieldFromConfig(conditionedBy(conditionedBy(linear(2002), { 2 }, { 1 }), { 2 }, { 1 })) == 2005,
           "…and a conditioner's own conditioner is counted too: 2 + (2 + 2001)");
        // …and it may be a CONTAINER, which is answered for with its worst submodel and then added.
        {
            nlohmann::json held = { { "architecture", "SlimmableContainer" },
                                    { "config", { { "submodels", nlohmann::json::array({
                                          { { "model", linear(66) } },
                                          { { "model", linear(2002) } } }) } } } };
            ok(receptiveFieldFromConfig(conditionedBy(held, { 2 }, { 1 })) == 2003,
               "…and a conditioner that is a CONTAINER contributes its longest submodel: 2 + 2001");
        }
        // …and a conditioned WaveNet can itself be a container's submodel.
        {
            nlohmann::json container = { { "architecture", "SlimmableContainer" },
                                         { "config", { { "submodels", nlohmann::json::array({
                                               { { "model", wavenet({ 2 }, { 1 }) } },
                                               { { "model", conditionedBy(linear(2002), { 2 }, { 1 }) } } }) } } } };
            ok(receptiveFieldFromConfig(container) == 2003,
               "…and a container answers through to a submodel's conditioner");
        }
        // 🔴 THE OTHER TWO FUNCTIONS OF THE REGISTRY READ THE SAME BRANCH, and that is not tidiness:
        // three functions giving divergent answers about ONE model is the class this house has been
        // bitten by before. A conditioner that is a Linear is charged the partitioned-FFT ring — its
        // instance owns the same engine and its ring feeds the layer arrays through the mixin — and a
        // conditioner that is an LSTM makes the whole capture recurrent, because the cell's state enters
        // every layer through a memoryless Conv1x1 and no finite silence empties it.
        ok(partitionedTailSamples(conditionedBy(linear(2002), { 2 }, { 1 })) == 2048,
           "a Linear CONDITIONER is charged the ring, exactly as a Linear submodel is");
        ok(isRecurrent(conditionedBy({ { "architecture", "LSTM" } }, { 2 }, { 1 })),
           "…and an LSTM CONDITIONER makes the capture recurrent");
        ok(partitionedTailSamples(conditionedBy(wavenet({ 2 }, { 1 }), { 2 }, { 1 })) == 0
               && ! isRecurrent(conditionedBy(wavenet({ 2 }, { 1 }), { 2 }, { 1 })),
           "…and a conditioner that is neither is charged neither");
        // THE JOINT ROW — one model read by all three functions at once, which is what a mutation that
        // teaches only ONE of them fails. A conditioner that is a container of a Linear and an LSTM
        // makes every one of the three answers move.
        {
            nlohmann::json mixed = { { "architecture", "SlimmableContainer" },
                                     { "config", { { "submodels", nlohmann::json::array({
                                           { { "model", linear(66) } },
                                           { { "model", { { "architecture", "LSTM" } } } } }) } } } };
            const auto model = conditionedBy(mixed, { 2, 2 }, { 1, 3 });
            ok(receptiveFieldFromConfig(model) == 70 && isRecurrent(model) && partitionedTailSamples(model) == 2048,
               "one model, three answers, all three through the conditioner: field "
               + std::to_string(receptiveFieldFromConfig(model)) + ", recurrent "
               + std::to_string((int) isRecurrent(model)) + ", ring "
               + std::to_string(partitionedTailSamples(model)));
        }
        // 🔴 AND A CONDITIONER INSIDE A CONTAINER'S SUBMODEL, for ALL THREE — which the rows above do not
        // reach: they put a container under a conditioner, and this is the other order. A walker that
        // recursed into `submodels` but not into THOSE submodels' conditioners passed every other row
        // here, so a container whose active WaveNet carries a dense Linear conditioner would lose its
        // 2048-sample ring (measured on one: the ledger reads 2102 for a model reaching 2386, and only
        // the ring covers the difference) and one whose WaveNet carries an LSTM would read as finite.
        {
            nlohmann::json held = { { "architecture", "SlimmableContainer" },
                                    { "config", { { "submodels", nlohmann::json::array({
                                          { { "model", wavenet({ 2 }, { 1 }) } },
                                          { { "model", conditionedBy(linear(2002), { 2 }, { 100 }) } } }) } } } };
            ok(receptiveFieldFromConfig(held) == 2102 && partitionedTailSamples(held) == 2048 && ! isRecurrent(held),
               "a container reached THROUGH to its submodel's conditioner, by all three: field "
               + std::to_string(receptiveFieldFromConfig(held)) + ", ring "
               + std::to_string(partitionedTailSamples(held)));
            nlohmann::json heldRecurrent = { { "architecture", "SlimmableContainer" },
                                             { "config", { { "submodels", nlohmann::json::array({
                                                   { { "model", wavenet({ 2 }, { 1 }) } },
                                                   { { "model", conditionedBy({ { "architecture", "LSTM" } },
                                                                              { 2 }, { 100 }) } } }) } } } };
            ok(isRecurrent(heldRecurrent),
               "…and an LSTM conditioner one level down still makes the whole tree recurrent");
        }
        // …AND THE WORST-OF IS A WORST-OF IN BOTH DIRECTIONS. A container carrying a stray `layers` key
        // and a dead conditioner SHORTER than its deepest submodel: an implementation that let the
        // own+conditioner sum WIN whenever it is non-zero, instead of taking the larger, answers 2 here.
        {
            nlohmann::json shallowSibling = { { "architecture", "SlimmableContainer" },
                                              { "config", { { "submodels", nlohmann::json::array({
                                                    { { "model", wavenet({ 2 }, { 500 }) } } }) },
                                                            { "layers", nlohmann::json::array() },
                                                            { "condition_dsp", linear(3) } } } };
            ok(receptiveFieldFromConfig(shallowSibling) == 501,
               "…and a sibling conditioner SHORTER than the deepest submodel does not replace it: "
               + std::to_string(receptiveFieldFromConfig(shallowSibling)));
        }
        // 🔴 AND THE SHAPES THAT MUST NOT BE CHARGED, because charging them would be inventing memory
        // no instance has. `condition_dsp` is read by exactly ONE parser in the library —
        // `nam::wavenet::parse_config_json` — so the key is DEAD on a Linear, an LSTM, a ConvNet and on
        // a CONTAINER, none of which build anything from it. The gate is the same shape NAM dispatches
        // on: a `layers` array.
        {
            nlohmann::json deadOnLinear = linear(3);
            deadOnLinear["config"]["condition_dsp"] = linear(2002);
            ok(receptiveFieldFromConfig(deadOnLinear) == 2,
               "a `condition_dsp` on a Linear config is DEAD JSON — NAM's Linear parser never reads it,"
               " so neither does this");
            nlohmann::json deadOnContainer = { { "architecture", "SlimmableContainer" },
                                               { "config", { { "submodels", nlohmann::json::array({
                                                     { { "model", wavenet({ 2 }, { 1 }) } } }) },
                                                             { "condition_dsp", linear(2002) } } } };
            ok(receptiveFieldFromConfig(deadOnContainer) == 2 && partitionedTailSamples(deadOnContainer) == 0
                   && ! isRecurrent(deadOnContainer),
               "…and a container's own sibling `condition_dsp` is dead the same way: its parser reads"
               " `submodels` and nothing else");
        }
        // 🔴 …AND THE CONFIGS THAT CARRY BOTH KEYS, which is where a SHAPE gate and NAM's own dispatch
        // part company. NAM dispatches on the `architecture` STRING (get_dsp.cpp:261) and never on
        // shape, so every one of these LOADS. What is asserted is that the three functions answer about
        // the SAME tree — a divergence between them is the class this task exists to close, and the
        // container branch returning early made exactly that: the field ignored a conditioner that
        // `isRecurrent` was charging.
        {
            nlohmann::json bothKeys = { { "architecture", "SlimmableContainer" },
                                        { "config", { { "submodels", nlohmann::json::array({
                                              { { "model", wavenet({ 2 }, { 500 }) } } }) },
                                                      { "layers", nlohmann::json::array() },
                                                      { "condition_dsp", linear(2002) } } } };
            ok(receptiveFieldFromConfig(bothKeys) == 2001 && partitionedTailSamples(bothKeys) == 2048
                   && ! isRecurrent(bothKeys),
               "a container carrying a stray `layers` key is answered about CONSISTENTLY by all three:"
               " field " + std::to_string(receptiveFieldFromConfig(bothKeys)) + ", ring "
               + std::to_string(partitionedTailSamples(bothKeys)) + " — the field used to return from"
               " `submodels` before it ever looked, while isRecurrent was already charging the same node");
            // …and the stack of a model carrying a stray `submodels` is not hidden by it either. A
            // SLIMMABLE WaveNet, because that is where NAM answers zero and this file is the only answer:
            // the early return read 1 for a model whose impulse reaches 4000.
            nlohmann::json strayContainer = wavenet({ 2 }, { 4000 });
            strayContainer["config"]["submodels"] = nlohmann::json::array({ { { "model", linear(100) } } });
            ok(receptiveFieldFromConfig(strayContainer) == 4001,
               "…and a model carrying a stray `submodels` array still answers for its own stack: "
               + std::to_string(receptiveFieldFromConfig(strayContainer)));
        }
        {
            nlohmann::json nulled = wavenet({ 2 }, { 1 });
            nulled["config"]["condition_dsp"] = nullptr;
            ok(receptiveFieldFromConfig(nulled) == 2 && partitionedTailSamples(nulled) == 0,
               "…and `\"condition_dsp\": null` means ABSENT, which is how NAM reads it too");
        }
        // THE SUM IS TAKEN IN long long. Both addends are already clamped to INT_MAX on their own, so
        // an int addition here is undefined behaviour for a config a test can write down — and this is
        // that config.
        ok(receptiveFieldFromConfig(conditionedBy(linear(2000000001), { 2 }, { 2000000000 }))
               == std::numeric_limits<int>::max(),
           "…and two saturating fields do not wrap when they are added");
    }

    group("the HEAD is memory too, and it is a branch this file did not read");
    {
        // 🔴 A layer array ends in `_head_rechannel`, a causal Conv1D whose kernel is the layer's own
        // `head.kernel_size`, and NAM charges `kernel − 1` for it on top of the dilations
        // (v0.5.4 wavenet/model.cpp:417-423, parsed at :882-899). The real captures carry it:
        // `example_models/A2.nam` spells it 16 in both submodels. Measured from outside this file: a
        // single dilation-1 layer with `head.kernel_size` 16 has an impulse response reaching sample 16,
        // where the same layer with the legacy `head_size` spelling reaches 1.
        //
        // WHY IT MATTERS EVEN THOUGH NAM USUALLY ANSWERS: `NamStage` raises this number with NAM's own,
        // and NAM's own includes the head — but a SlimmableWavenet answers ZERO for everything
        // (wavenet/slimmable.h:66) and there is nothing to raise it with. Measured on one that loads:
        // NAM 0, this file 2 before, 17 after, impulse reach 16. It is also what keeps the CONDITIONER
        // sum an upper bound at all: the conditioner's memory enters a layer AFTER that layer's own
        // convolution, so what the sum has to cover is `Mₒ + H + M_c − L₀`, and an unread H makes it
        // short by `H − L₀ − 1` however correctly the conditioner is walked.
        const auto headed = [](int kernel, std::vector<int> dilations) {
            return nlohmann::json { { "architecture", "WaveNet" },
                                    { "config", { { "layers", nlohmann::json::array({
                                          { { "kernel_size", 2 }, { "dilations", dilations },
                                            { "head", { { "out_channels", 1 }, { "kernel_size", kernel },
                                                        { "bias", false } } } } }) } } } };
        };
        ok(receptiveFieldFromConfig(headed(16, { 1 })) == 17,
           "a head kernel of 16 reaches back 15 further samples: 1 + 15, plus the file's own sample of margin");
        ok(receptiveFieldFromConfig(headed(2, { 1 })) == 3,
           "…a head kernel of 2 reaches back ONE further sample, which is the boundary of the `k > 1` guard");
        ok(receptiveFieldFromConfig(headed(1, { 1 })) == 2,
           "…a head kernel of 1 is a rechannel and has no memory, which is what the LEGACY `head_size`"
           " spelling means");
        ok(receptiveFieldFromConfig(wavenet({ 2 }, { 1 })) == 2,
           "…and the legacy spelling itself is unchanged: nothing is charged where there is no `head` object");
        ok(receptiveFieldFromConfig(headed(16, { 1, 100 })) == 117,
           "…and it is charged for the ARRAY, once, on top of every dilation in it");
        // …AND THE POST-STACK HEAD, a different key with a different shape: `config.head` is a chain of
        // causal convolutions and NAM charges `Σ(kᵢ − 1)` (wavenet/model.cpp:58-67, :620, parsed :1154-1186).
        {
            nlohmann::json posted = wavenet({ 2 }, { 1 });
            posted["config"]["head"] = { { "channels", 1 }, { "out_channels", 1 },
                                         { "kernel_sizes", nlohmann::json::array({ 4, 7 }) },
                                         { "activation", "Tanh" } };
            ok(receptiveFieldFromConfig(posted) == 11, "a post-stack head of kernels {4,7} adds 3 + 6");
            posted["config"]["head"] = nullptr;
            ok(receptiveFieldFromConfig(posted) == 2, "…and `\"head\": null` means ABSENT, as NAM reads it");
        }
        // The shape every capture in this project actually is, with its head: A2.nam's stack answers
        // 6332 without it and 6347 with — which is EXACTLY what NAM answers for the same file
        // (measured through `::nam::get_dsp`), and the model's impulse reaches 5899 of it.
        {
            const std::vector<int> ks { 6,6,6,6,6,6,6, 6,6,6,6,6,6,6, 15,15, 6,6,6,6,6,6,6 };
            const std::vector<int> ds { 1,3,7,17,41,101,239, 1,3,7,17,41,101,239, 1,13, 1,3,7,17,41,101,239 };
            nlohmann::json a2 = wavenet(ks, ds);
            a2["config"]["layers"][0]["head"] = { { "out_channels", 1 }, { "kernel_size", 16 }, { "bias", true } };
            ok(receptiveFieldFromConfig(a2) == 6347,
               "the real A2 stack with its 16-tap head: 6331 + 15 + 1, which is NAM's own answer for the"
               " shipped file to the sample");
        }
    }

    group("a ConvNet declares its stack at the TOP level, and nothing here was reading it");
    {
        // 🔴 A ConvNet is a chain of dilated convolutions with a kernel NAM hard-codes to 2
        // ("HACK 2 kernel", v0.5.4 convnet.cpp:56-57) and a field of `1 + Σ dilations`
        // (convnet.cpp:200-202) read off a TOP-LEVEL `dilations` array — the key a WaveNet keeps one
        // level down, inside its `layers` groups. This file read neither, so EVERY ConvNet answered 0.
        // Measured on `dilations:[1,2,4,8]`: this file answered 0, NAM answers 16, the impulse reaches
        // sample 15. Masked by NAM at the top level and behind a plain WaveNet's conditioner; not
        // masked behind a SlimmableWavenet, which answers zero for everything.
        const auto convnet = [](std::vector<int> dilations) {
            return nlohmann::json { { "architecture", "ConvNet" },
                                    { "config", { { "channels", 1 }, { "dilations", dilations } } } };
        };
        ok(receptiveFieldFromConfig(convnet({ 1, 2, 4, 8 })) == 16, "1 + 1 + 2 + 4 + 8, as NAM computes it");
        ok(receptiveFieldFromConfig(convnet({ })) == 0, "…and an empty stack is zero, not one");
        ok(receptiveFieldFromConfig(conditionedBy(convnet({ 1, 2, 4, 8 }), { 2 }, { 1 })) == 18,
           "…and a ConvNet CONDITIONER is added like any other: 2 + 16");
        // 🔴 A top-level `dilations` beside real `layers` is set aside, and COSTS THE ALLOWANCE (P92,
        // moved on purpose; this row asserted 2). NAM's WaveNet parser never looks there, and the stack
        // is the reading; but a reading set aside is one this file cannot place — see the declared row
        // above for why neither trusting it nor ignoring it is safe.
        {
            nlohmann::json both = wavenet({ 2 }, { 1 });
            both["config"]["dilations"] = nlohmann::json::array({ 1000 });
            ok(receptiveFieldFromConfig(both) == 2 + 48000,
               "…and a stray top-level `dilations` beside real `layers` is the stack plus the allowance: "
               + std::to_string(receptiveFieldFromConfig(both)) + " (was 2 before P92)");
            // …and the gate that used to suppress this reader for a NON-EMPTY `layers` is gone: NAM's
            // ConvNet parser never reads `layers`, so a ConvNet carrying `"layers":[{}]` LOADS, and this
            // answered 0 where NAM's field is 256. An entry with no dilations reads as no stack, so
            // nothing is set aside and nothing is charged.
            nlohmann::json deadNonEmpty = convnet({ 1, 2, 4, 8, 16, 32, 64, 128 });
            deadNonEmpty["config"]["layers"] = nlohmann::json::array({ nlohmann::json::object() });
            ok(receptiveFieldFromConfig(deadNonEmpty) == 256 && partitionedTailSamples(deadNonEmpty) == 0,
               "…which is what closes a ConvNet carrying a NON-EMPTY dead `layers`: 1 + 255, as NAM computes"
               " it, where the gate answered 0 — read " + std::to_string(receptiveFieldFromConfig(deadNonEmpty)));
        }
        ok(convNetField(nlohmann::json::object()) == 0, "…and a config with no stack at all is still zero");
        // 🔴 AND A DEAD `layers` KEY MUST NOT SUPPRESS IT. NAM's ConvNet parser never reads `layers`
        // (`convnet.cpp:349-355`), so a ConvNet carrying `"layers": []` had its whole stack silenced
        // here — the same shape of mistake as a chain letting one reader speak over another. Measured on
        // one that loads, as a slimmable capture's conditioner where NAM answers zero and this file is
        // the only answer: the lane drained 102 for a model reaching 131 and handed back 29 samples.
        {
            nlohmann::json deadLayers = convnet({ 1, 2, 4, 8, 16 });
            deadLayers["config"]["layers"] = nlohmann::json::array();
            ok(receptiveFieldFromConfig(deadLayers) == 32,
               "a ConvNet carrying a dead `layers: []` still answers for its own stack: "
               + std::to_string(receptiveFieldFromConfig(deadLayers)));
        }
        // 🔴 A NUMBER THAT DOES NOT FIT IS REFUSED, NOT CAST. `is_number()` is deliberately wide (NAM's
        // own `get<int>()` takes `2001.0`), but static_casting a float that is out of the integer's
        // range is UNDEFINED BEHAVIOUR, and this repository's CI has a hard-fail UBSan job.
        // 🔴 AND A REFUSAL IS THE ALLOWANCE, NOT "ABSENT" (P92; these rows asserted the value without the
        // refused entry). "Absent" relied on NAM's own number covering what NAM built, and a slimmable
        // WaveNet answers `return 0`: `dilations:[4294967396]` loads as dilation 100, reaches 100, and
        // this file answered 0 for it.
        ok(receptiveFieldFromConfig({ { "architecture", "Linear" }, { "config", { { "receptive_field", 1e300 } } } }) == 48000,
           "a field of 1e300 is refused, not cast — the cast is UB — and charged the allowance");
        {
            nlohmann::json hugeHead = { { "architecture", "WaveNet" },
                                        { "config", { { "layers", nlohmann::json::array({
                                              { { "kernel_size", 2 }, { "dilations", nlohmann::json::array({ 1 }) },
                                                { "head", { { "out_channels", 1 }, { "kernel_size", 1e300 },
                                                            { "bias", false } } } } }) } } } };
            ok(receptiveFieldFromConfig(hugeHead) == 2 + 48000,
               "…and so is a head kernel of 1e300: the stack's 2 plus the allowance");
            nlohmann::json hugeConv = { { "architecture", "ConvNet" },
                                        { "config", { { "channels", 1 },
                                                      { "dilations", nlohmann::json::array({ 1e300, 4 }) } } } };
            ok(receptiveFieldFromConfig(hugeConv) == 5 + 48000,
               "…and one dilation of 1e300 is skipped, not carried, and costs the allowance beside its neighbour's 5");
            // …and a product that would overflow is CLAMPED, not wrapped. 2e9 x 2e9 is 4e18, past int
            // and past a careless long long accumulation; a wrap here would answer a NEGATIVE field.
            nlohmann::json overflowing = { { "architecture", "WaveNet" },
                                           { "config", { { "layers", nlohmann::json::array({
                                                 { { "kernel_sizes", nlohmann::json::array({ 2000000001 }) },
                                                   { "dilations", nlohmann::json::array({ 2000000000 }) } } }) } } } };
            ok(receptiveFieldFromConfig(overflowing) == std::numeric_limits<int>::max(),
               "…and a dilation-by-kernel product past the int is clamped, never wrapped: "
               + std::to_string(receptiveFieldFromConfig(overflowing)));
            // 🔴 …AND A DILATION THAT DOES NOT FIT AN `int` IS REFUSED, not clamped, because the model
            // NAM BUILT does not contain it: NAM's dilation is an int, and `4294967396` builds a network
            // whose dilation is 100 — a file that loads. Clamping to INT_MAX would be an upper bound
            // that is true and useless, and it is SPENT: 2 147 483 646 samples of inference, a measured
            // 23.7 s of synchronous audio-thread work in reset(), for a model that remembers a hundred.
            // Refused is the allowance: NAM's own answer covers what it built only where NAM answers.
            nlohmann::json pastTheInt = { { "architecture", "WaveNet" },
                                          { "config", { { "layers", nlohmann::json::array({
                                                { { "kernel_size", 2 },
                                                  { "dilations", nlohmann::json::array({ 4294967396LL, 4 }) } } }) } } } };
            ok(receptiveFieldFromConfig(pastTheInt) == 5 + 48000,
               "a dilation past the int is refused, its neighbour still counts, and the refusal costs the allowance: "
               + std::to_string(receptiveFieldFromConfig(pastTheInt)));
            // 🔴 AND A BOOLEAN IS A NUMBER TO NAM. `get<int>(true)` is 1 — measured — so NAM loads
            // `"dilations":[true,true]` as `[1,1]` and answers 3 for it, while `is_number()` says no and
            // `get<long long>()` THROWS. Read as nothing, that was a real under-drain on a slimmable
            // stack, where NAM's own answer is zero too: two samples at 0.1875 out of digital silence.
            nlohmann::json booleans = { { "architecture", "WaveNet" },
                                        { "config", { { "layers", nlohmann::json::array({
                                              { { "kernel_size", 2 },
                                                { "dilations", nlohmann::json::array({ true, true }) } } }) } } } };
            ok(receptiveFieldFromConfig(booleans) == 3,
               "…and a dilation spelled `true` is the 1 that NAM builds from it: "
               + std::to_string(receptiveFieldFromConfig(booleans)));
            booleans["config"]["layers"][0]["dilations"] = nlohmann::json::array({ false, true });
            ok(receptiveFieldFromConfig(booleans) == 2,
               "…and one spelled `false` is the 0 NAM builds from it, not a 1: "
               + std::to_string(receptiveFieldFromConfig(booleans)));
            // 🔴 AND A NEGATIVE PAST THE INT IS REFUSED TOO, not read as a negative (P92). NAM's `get<int>()`
            // wraps `-4294967196` to 100 — a slimmable carrying it LOADS and reaches 100 — so reading it as a
            // negative number and clamping it to 0 would answer zero for that memory. Found by a mutant
            // that dropped the lower bound and survived the suite; `-1e300` is the same bound on the float
            // path, where dropping it is an undefined cast besides.
            nlohmann::json negativePast = { { "architecture", "WaveNet" },
                                            { "config", { { "layers", nlohmann::json::array({
                                                  { { "kernel_size", 2 },
                                                    { "dilations", nlohmann::json::array({ -4294967196LL }) } } }) } } } };
            ok(receptiveFieldFromConfig(negativePast) == 48000,
               "a dilation below the int is refused and costs the allowance: "
               + std::to_string(receptiveFieldFromConfig(negativePast)));
            ok(receptiveFieldFromConfig({ { "architecture", "Linear" }, { "config", { { "receptive_field", -1e300 } } } }) == 48000,
               "…and so is a field of -1e300, whose cast would be undefined");
            // …AND THE KEY GUARDS INSIDE THE HEAD READERS AND THE PER-LAYER KERNELS, each of which stands
            // between a const `operator[]` and a key or an index that is not there — undefined behaviour,
            // not a zero. Each was removed by a mutant that the suite did not notice.
            nlohmann::json headNoKernel = wavenet({ 2 }, { 1 });
            headNoKernel["config"]["layers"][0]["head"] = { { "out_channels", 1 } };
            ok(receptiveFieldFromConfig(headNoKernel) == 2,
               "a layer head with no `kernel_size` is the legacy kernel of 1, read without indexing it");
            nlohmann::json postNoKernels = wavenet({ 2 }, { 1 });
            postNoKernels["config"]["head"] = { { "channels", 1 } };
            ok(receptiveFieldFromConfig(postNoKernels) == 2,
               "…and a post-stack head with no `kernel_sizes` adds nothing, read without indexing it");
            ok(receptiveFieldFromConfig(wavenet({ 2 }, { 1, 2 })) == 2,
               "…and `kernel_sizes` shorter than `dilations` stops at the last kernel instead of reading past it");
            // …and an ACCUMULATION past it: four layer arrays that each clamp to INT_MAX still sum
            // inside the type, because the running total is clamped after every one.
            nlohmann::json accumulating = { { "architecture", "WaveNet" },
                                            { "config", { { "layers", nlohmann::json::array({
                                                  { { "kernel_sizes", nlohmann::json::array({ 2000000001LL }) },
                                                    { "dilations",    nlohmann::json::array({ 2000000000LL }) } },
                                                  { { "kernel_sizes", nlohmann::json::array({ 2000000001LL }) },
                                                    { "dilations",    nlohmann::json::array({ 2000000000LL }) } },
                                                  { { "kernel_sizes", nlohmann::json::array({ 2000000001LL }) },
                                                    { "dilations",    nlohmann::json::array({ 2000000000LL }) } },
                                                  { { "kernel_sizes", nlohmann::json::array({ 2000000001LL }) },
                                                    { "dilations",    nlohmann::json::array({ 2000000000LL }) } } }) } } } };
            ok(receptiveFieldFromConfig(accumulating) == std::numeric_limits<int>::max(),
               "…and four of them accumulate to the ceiling rather than through it: "
               + std::to_string(receptiveFieldFromConfig(accumulating)));
        }
        // 🔴 AND IT MUST NOT SUPPRESS THE DECLARED FIELD. A `Linear` carrying a stray `dilations` array
        // LOADS — its parser reads neither key (linear.cpp:306-316) — and a first-non-zero chain made
        // this reader answer for the stray key INSTEAD of the real one: measured through NamStage,
        // 2 for a 4999-sample impulse response that answers 4999 without the stray key, i.e. a lane
        // draining 2050 where it needs 4999. The three sources are the stack, then the WORST of the
        // other two; a chain can take a number away and a max cannot.
        {
            nlohmann::json strayDilations = linear(5000);
            strayDilations["config"]["dilations"] = nlohmann::json::array({ 1 });
            ok(receptiveFieldFromConfig(strayDilations) == 4999,
               "a Linear capture carrying a stray `dilations` array still answers its DECLARED field: "
               + std::to_string(receptiveFieldFromConfig(strayDilations)));
        }
    }

    group("P92: what this file cannot place costs ONE allowance — never the face value, never zero");
    {
        // 🔴 WHAT THE REGISTRY PROMISES, FIRST. An UPPER BOUND on the memory of the whole model. What it
        // PLACES it trusts; for anything it cannot place — an unread config, a refused value, a reading
        // it sets aside — it adds ONE allowance for the whole tree. Not zero: zero is the previous sound
        // coming out of silence. Not the face value: a dead number is unbounded and is spent on the audio
        // thread. The number is a POLICY literal; a legal change to it must update these rows on purpose.
        ok(kUnreadShapeCeiling == 48000 && kMaxUnplacedHops == 32,
           "the allowance is 48 000 samples, and the recurrence walk follows 32 unplaced hops");

        // THE SHAPE THAT MADE THE RULE — NAM's own shipped capture, rewrapped. It LOADS (the WaveNet
        // parser delegates on the top-level marker and the delegate reads `config.model`), it is the same
        // network, and its impulse reaches sample 2046 either way. Before P92: 2047 flat, 0 wrapped.
        const auto flat = flatSlimmable(realSlimmableConfig());
        const auto wrap = wrapped(realSlimmableConfig());
        ok(receptiveFieldFromConfig(flat) == 2047 && partitionedTailSamples(flat) == 0 && ! isRecurrent(flat),
           "the FLAT shipped slimmable is unchanged: field 2047 (2*1023 + 1, by hand), no ring, not recurrent");
        ok(receptiveFieldFromConfig(wrap) == 48000 && partitionedTailSamples(wrap) == 2048 && ! isRecurrent(wrap),
           "the WRAPPED one answers the allowance and the ring, where it answered 0 — read "
           + std::to_string(receptiveFieldFromConfig(wrap)) + ", ring " + std::to_string(partitionedTailSamples(wrap)));

        // 🔴 THE JOINT ROW — one model, all three readers, one assertion. A fix made in ONE of the three
        // functions fails it: the wrapper hides a conditioner that is a container of a Linear and an LSTM.
        // Before P92 all three read nothing: 0, 0, false.
        {
            nlohmann::json inner = realSlimmableConfig();
            inner["condition_dsp"] = { { "architecture", "SlimmableContainer" },
                                       { "config", { { "submodels", nlohmann::json::array({
                                             { { "model", linear(66) } },
                                             { { "model", { { "architecture", "LSTM" } } } } }) } } } };
            const auto model = wrapped(inner);
            ok(receptiveFieldFromConfig(model) == 48000 && partitionedTailSamples(model) == 2048 && isRecurrent(model),
               "one wrapped model, three answers, all three moved by the wrapper: field "
               + std::to_string(receptiveFieldFromConfig(model)) + ", ring "
               + std::to_string(partitionedTailSamples(model)) + ", recurrent "
               + std::to_string((int) isRecurrent(model)));
            // …and the recurrence is READ, not assumed: the same wrapper without the LSTM is not recurrent.
            nlohmann::json plainInner = realSlimmableConfig();
            plainInner["condition_dsp"] = linear(66);
            ok(! isRecurrent(wrapped(plainInner)),
               "a wrapper is NOT declared recurrent for being a wrapper — recurrence stays a read fact");
        }

        // 🔴 ONCE — measured failure of the per-node floor: a 30-byte dead sibling cost a full allowance
        // EACH, and 50 000 of them (a 1.6 MB file NAM loads) drained INT_MAX.
        for (const int n : { 1, 10, 1000 })
        {
            nlohmann::json many = wavenet({ 2 }, { 100 });
            for (int i = 0; i < n; ++i)
                many["config"]["dead" + std::to_string(i)] = { { "layers", nlohmann::json::array() },
                                                               { "m", { { "layers", nlohmann::json::array() } } } };
            ok(receptiveFieldFromConfig(many) == 101 + 48000,
               std::to_string(n) + " unplaced siblings cost the allowance ONCE: "
               + std::to_string(receptiveFieldFromConfig(many)));
        }
        // 🔴 ADDED, NOT MAXED — measured failure of a max at the root: a known 100 001 swallowed the
        // allowance of an unreadable stage in series with it, i.e. charged that stage zero.
        {
            nlohmann::json known = wavenet({ 2 }, { 100000 });
            known["config"]["stage"] = { { "layers", nlohmann::json::array() } };
            ok(receptiveFieldFromConfig(known) == 100001 + 48000,
               "a read stack of 100 001 plus an unreadable stage is 148 001, not 100 001: "
               + std::to_string(receptiveFieldFromConfig(known)));
        }
        // 🔴 NOT THE FACE VALUE — measured failure of trusting it: a dead key costs nothing to NAM and
        // was spent as 2^31 samples here. This is also the rule's stated PRICE: a LIVE wrapped model
        // longer than the allowance is under-drained by the difference. No real capture is (6347 max).
        {
            nlohmann::json deep = { { "layers", nlohmann::json::array({
                                        { { "kernel_size", 2 }, { "dilations", nlohmann::json::array({ 100000 }) } } }) } };
            ok(receptiveFieldFromConfig(wrapped(deep)) == 48000,
               "a wrapped stack of 100 001 is charged the allowance, not its face value — the price of a"
               " rule that a dead key cannot inflate: " + std::to_string(receptiveFieldFromConfig(wrapped(deep))));
            // …and the same price on the other path the rule does not trust: a declared field SET ASIDE
            // beside a (dead) stack. A 60 001-tap Linear with a readable stray `layers` array drains the
            // stack's 2 plus the allowance — short of its 60 000 by 9 950, measured (door D1).
            nlohmann::json setAside = linear(60001);
            setAside["config"]["layers"] = nlohmann::json::array({
                { { "kernel_size", 2 }, { "dilations", nlohmann::json::array({ 1 }) } } });
            ok(receptiveFieldFromConfig(setAside) == 2 + 48000,
               "a declared field longer than the allowance, set aside beside a stack, is charged the allowance"
               " (door D1, registered): " + std::to_string(receptiveFieldFromConfig(setAside)));
            // 🔴 AND THE OTHER DOOR (D2), pinned so that changing the policy moves this row on purpose: a
            // DEAD number the file PLACES is trusted at face value, as base did. The wrapped form's decoy
            // `layers` are placed — NAM builds from `config.model` — so a decoy spelling a two-million
            // dilation drains two million samples, plus the allowance for the unplaced config beside it.
            nlohmann::json decoy = wrapped(realSlimmableConfig());
            decoy["config"]["layers"][0]["kernel_size"] = 2;
            decoy["config"]["layers"][0]["dilations"] = nlohmann::json::array({ 2000000 });
            ok(receptiveFieldFromConfig(decoy) == 2000001 + 48000,
               "a dead decoy stack in the wrapped form is placed and trusted (door D2, registered — base"
               " answered 2 000 001): " + std::to_string(receptiveFieldFromConfig(decoy)));
            nlohmann::json deadHuge = wavenet({ 2 }, { 100 });
            deadHuge["config"]["notes"] = { { "receptive_field", 2147483647 } };
            ok(receptiveFieldFromConfig(deadHuge) == 101 + 48000,
               "…and a dead `receptive_field` of INT_MAX under an unread key costs the same allowance: "
               + std::to_string(receptiveFieldFromConfig(deadHuge)));
        }

        // A NODE THAT CARRIES THE VOCABULARY AND READS AS NOTHING is the plainest "I do not know". Where
        // zero would be honest (`receptive_field: 1` is a gain) the allowance is an over-charge, and the
        // rule has already chosen that side.
        ok(receptiveFieldFromConfig(wrapped({ { "layers", nlohmann::json::array() } })) == 48000
               && receptiveFieldFromConfig(wrapped({ { "receptive_field", 1 } })) == 48000
               && receptiveFieldFromConfig(wrapped({ { "dilations", nlohmann::json::array({ "x" }) } })) == 48000,
           "an unplaced node that reads as NOTHING costs the allowance, not zero");

        // A REFUSED VALUE on the architecture where nothing else answers: a slimmable dilation past the
        // int, which NAM builds as 100 and whose impulse reaches 100. It answered 0.
        {
            nlohmann::json refused = realSlimmableConfig();
            refused["layers"][0]["dilations"] = nlohmann::json::array({ 4294967396LL });
            ok(receptiveFieldFromConfig(flatSlimmable(refused)) == 48000
                   && partitionedTailSamples(flatSlimmable(refused)) == 2048,
               "a slimmable spelling a dilation past the int costs the allowance and the ring — it answered 0");
        }

        // 🔴 THE CLASS, NOT THE ADDRESS. Nothing below is spelled the way NAM spells its wrapper today.
        {
            ok(receptiveFieldFromConfig(wrapped(realSlimmableConfig(), "future_wrapper")) == 48000,
               "the same config under a key NAM does not use still costs the allowance");
            nlohmann::json asArray = flatSlimmable({ { "layers", nlohmann::json::array() } });
            asArray["config"]["variants"] = nlohmann::json::array({ realSlimmableConfig(), realSlimmableConfig() });
            ok(receptiveFieldFromConfig(asArray) == 48000 && partitionedTailSamples(asArray) == 2048,
               "…and so does an ARRAY of configs under such a key — the `submodels` shape by another name");
            nlohmann::json asNode = flatSlimmable({ { "layers", nlohmann::json::array() } });
            asNode["config"]["alt"] = { { "architecture", "LSTM" }, { "config", { { "hidden_size", 3 } } } };
            ok(receptiveFieldFromConfig(asNode) == 48000 && isRecurrent(asNode),
               "…and a whole MODEL NODE under such a key is walked for its architecture");
            nlohmann::json outer = flatSlimmable({ { "layers", nlohmann::json::array() } });
            outer["config"]["outer"] = { { "model", realSlimmableConfig() } };
            ok(receptiveFieldFromConfig(outer) == 48000,
               "…and a wrapper of a wrapper, whose only vocabulary word is `model`");
            nlohmann::json hidden = flatSlimmable({ { "layers", nlohmann::json::array() } });
            hidden["config"]["x"] = { { "architecture", "WaveNet" },
                                      { "config", { { "weird", { { "architecture", "LSTM" },
                                                                 { "config", { { "hidden_size", 3 } } } } } } } };
            ok(isRecurrent(hidden),
               "…and a model node's config is walked even when it carries no vocabulary word — an LSTM one"
               " level further down is still found");
            nlohmann::json inLayer = flat;
            inLayer["config"]["layers"][0]["sidechain"] = { { "layers", nlohmann::json::array({
                { { "kernel_size", 2 }, { "dilations", nlohmann::json::array({ 70000 }) } } }) } };
            ok(receptiveFieldFromConfig(inLayer) == 2047 + 48000,
               "…and a config hanging off a LAYER ENTRY is the allowance on top of the stack it rides on: "
               + std::to_string(receptiveFieldFromConfig(inLayer)));
            nlohmann::json inLayerArray = flat;
            inLayerArray["config"]["layers"][0]["sidechains"] = nlohmann::json::array({
                { { "layers", nlohmann::json::array() } } });
            ok(receptiveFieldFromConfig(inLayerArray) == 2047 + 48000,
               "…and an ARRAY of them there too");
            nlohmann::json inSubmodel = { { "architecture", "SlimmableContainer" },
                                          { "config", { { "submodels", nlohmann::json::array({
                                                { { "max_value", 1.0 }, { "model", linear(2) },
                                                  { "extra", { { "layers", nlohmann::json::array() } } } } }) } } } };
            ok(receptiveFieldFromConfig(inSubmodel) == 1 + 48000 && partitionedTailSamples(inSubmodel) == 2048,
               "…and one hanging off a SUBMODEL ENTRY beside its `model`");
        }

        // 🔴 AND WHERE IT MUST NOT FIRE — identity, because a false fire is ~40 ms of a real capture.
        {
            // Every object a real layer entry carries is a layer FEATURE, and none is a model.
            nlohmann::json featured = flat;
            auto& grp = featured["config"]["layers"][0];
            grp["head1x1"] = { { "active", false }, { "out_channels", 1 }, { "groups", 1 } };
            grp["layer1x1"] = { { "active", true }, { "groups", 1 } };
            for (const char* film : { "conv_pre_film", "conv_post_film", "input_mixin_pre_film",
                                      "input_mixin_post_film", "activation_pre_film", "activation_post_film",
                                      "layer1x1_post_film", "head1x1_post_film" })
                grp[film] = { { "active", false }, { "shift", true }, { "groups", 1 } };
            grp["activation"] = { { "type", "PReLU" } };
            ok(receptiveFieldFromConfig(featured) == 2047 && partitionedTailSamples(featured) == 0,
               "the layer features of a real A2-class entry (head1x1, layer1x1, eight FiLMs, an activation"
               " object) do not fire it");
            // `metadata` sits on the MODEL node, beside `config`, and is never visited — even when a
            // user's free-form tree happens to use the vocabulary.
            nlohmann::json tagged = flat;
            tagged["metadata"] = { { "layers", nlohmann::json::array({ "di", "amp" }) },
                                   { "training", { { "model", { { "layers", nlohmann::json::array() } } } } } };
            ok(receptiveFieldFromConfig(tagged) == 2047 && partitionedTailSamples(tagged) == 0,
               "a model's `metadata` carrying vocabulary words is not read");
            // The post-stack `head` is a key this file reads, so a stray `layers` inside it is not a model.
            nlohmann::json headed = wavenet({ 2 }, { 1 });
            headed["config"]["head"] = { { "kernel_sizes", nlohmann::json::array({ 4 }) },
                                         { "layers", nlohmann::json::array() } };
            ok(receptiveFieldFromConfig(headed) == 5, "a `layers` key inside the post-stack head does not fire it");
            // A ConvNet's `activation` is an OBJECT at CONFIG level, and it is not a model.
            nlohmann::json conv = { { "architecture", "ConvNet" },
                                    { "config", { { "channels", 1 }, { "dilations", nlohmann::json::array({ 1, 2, 4, 8 }) },
                                                  { "activation", { { "type", "ReLU" } } } } } };
            ok(receptiveFieldFromConfig(conv) == 16 && partitionedTailSamples(conv) == 0,
               "a ConvNet's object-valued `activation` does not fire it");
            // Configs that are not objects, and a model that is not one, are refused rather than thrown on.
            bool quiet = true;
            for (const auto& cfg : { nlohmann::json(nullptr), nlohmann::json::array(), nlohmann::json(7) })
            {
                const nlohmann::json m = { { "architecture", "WaveNet" }, { "config", cfg } };
                quiet = quiet && receptiveFieldFromConfig(m) == 0 && partitionedTailSamples(m) == 0 && ! isRecurrent(m);
            }
            const auto topArray = nlohmann::json::array({ 1, 2 });
            quiet = quiet && receptiveFieldFromConfig(topArray) == 0 && partitionedTailSamples(topArray) == 0;
            nlohmann::json scalarWrap = flatSlimmable({ { "layers", nlohmann::json::array() }, { "model", 7 } });
            quiet = quiet && receptiveFieldFromConfig(scalarWrap) == 0;
            ok(quiet, "`config: null / [] / 7`, a top-level array and a scalar `model` all answer 0 without throwing");
            // …and the two shapes where NOT iterating is a decision rather than a no-op: an ARRAY config
            // and an ARRAY layer entry. NAM cannot load either (its parsers index them by key and throw),
            // so the answer is never spent; iterating would read their elements as unplaced configs.
            const nlohmann::json arrayConfig = { { "architecture", "WaveNet" },
                                                 { "config", nlohmann::json::array({
                                                       { { "layers", nlohmann::json::array() } } }) } };
            const nlohmann::json arrayEntry = { { "architecture", "WaveNet" },
                                                { "config", { { "layers", nlohmann::json::array({
                                                      nlohmann::json::array({ { { "layers", nlohmann::json::array() } } }) }) } } } };
            ok(receptiveFieldFromConfig(arrayConfig) == 0 && receptiveFieldFromConfig(arrayEntry) == 0,
               "an array config and an array layer entry are not iterated: "
               + std::to_string(receptiveFieldFromConfig(arrayConfig)) + " / "
               + std::to_string(receptiveFieldFromConfig(arrayEntry)));
        }

        // 🔴 THE AXES BASE WALKED ARE NOT TRUNCATED. A first draft guarded ALL nesting at 64 levels and
        // answered the allowance past it — which is SHORT of a model NAM loads: 65 WaveNets chained
        // through `condition_dsp` into a 60 001-tap Linear drained 50 243 where base drained 62 243, and
        // a restarted lane replayed 1.95484 out of silence. These axes are read as base read them.
        ok(receptiveFieldFromConfig(nlohmann::json::parse(conditionersInto(65, R"({"architecture":"Linear","config":{"receptive_field":60001}})"))) == 60000,
           "65 conditioners into a 60 001-tap Linear are read to the end: 60 000, as base answered");
        ok(isRecurrent(nlohmann::json::parse(conditionersInto(70, R"({"architecture":"LSTM","config":{"hidden_size":3}})"))),
           "…and an LSTM at the end of 70 conditioners is still recurrent, as base answered");
        {
            std::string nested;
            for (int i = 0; i < 70; ++i) nested += R"({"architecture":"SlimmableContainer","config":{"submodels":[{"model":)";
            nested += R"({"architecture":"LSTM","config":{"hidden_size":3}})";
            for (int i = 0; i < 70; ++i) nested += "}]}}";
            ok(isRecurrent(nlohmann::json::parse(nested)),
               "…and so is one at the end of 70 nested containers — neither axis base walked is counted");
        }

        // 🔴 THE UNPLACED AXIS IS GUARDED, and only it — the one recursion P92 added. Within the guard an
        // LSTM is found; past it the walk stops and recurrence is not assumed; the field carries the
        // allowance either way. And a 100 000-hop file is walked on a worker thread without a crash: a
        // std::thread's default stack is 512 KiB on macOS, 1 MiB on Windows, 8 MiB on glibc, and the
        // unguarded walk died at ~1000 hops on the smallest of them.
        {
            const auto near = nlohmann::json::parse(lstmUnderUnplaced(kMaxUnplacedHops));
            const auto far  = nlohmann::json::parse(lstmUnderUnplaced(kMaxUnplacedHops + 1));
            ok(isRecurrent(near) && ! isRecurrent(far)
                   && receptiveFieldFromConfig(near) == 48000 && receptiveFieldFromConfig(far) == 48000
                   && partitionedTailSamples(far) == 2048,
               "an LSTM 32 unplaced hops down is found, one 33 hops down is not assumed, and both cost the"
               " allowance and the ring");
            int field = -1, ring = -1;
            bool rec = true;
            std::thread worker([&] {
                const nlohmann::json hostile = nlohmann::json::parse(wrapperChain(100000));
                field = receptiveFieldFromConfig(hostile);
                ring  = partitionedTailSamples(hostile);
                rec   = isRecurrent(hostile);
            });
            worker.join();
            ok(field == 2 + 48000 && ring == 2048 && ! rec,
               "100 000 dead `model` keys, on a worker thread: no crash, the top stack's 2 plus the allowance —"
               " read " + std::to_string(field));
        }

        // 🔴 AND LINEAR, NOT EXPONENTIAL. A first draft read an unplaced node twice — as a model node
        // through its `config`, and as a raw config whose `config` key is itself unplaced — so every
        // level was reached by two paths and the work grew about 1.6x per level. Asked two ways: by
        // COUNTING the predicate calls (deterministic), and by bounding the time within the guard, where
        // the doubled walk takes seconds and the single one microseconds.
        {
            const nlohmann::json chain30 = nlohmann::json::parse(doubleReadingChain(30));
            predCalls.store(0);
            (void) anyNestedModel(chain30, countingNever);
            ok(predCalls.load() <= 2 * 30 + 8,
               "the walk asks each level ONCE: " + std::to_string(predCalls.load())
               + " predicate calls over 30 levels (the doubled walk asks about 1.6^30)");
            const nlohmann::json chain32 = nlohmann::json::parse(doubleReadingChain(kMaxUnplacedHops));
            const auto t0 = std::chrono::steady_clock::now();
            const bool r32 = isRecurrent(chain32);
            const double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
            ok(! r32 && ms < 250.0 && receptiveFieldFromConfig(chain32) == 48000,
               "…and 32 hops take " + std::to_string(ms) + " ms (the doubled walk takes seconds here)");
        }
    }

    return felitronics::test::report();
}

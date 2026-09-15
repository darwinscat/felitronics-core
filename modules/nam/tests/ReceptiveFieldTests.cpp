// SPDX-License-Identifier: MIT
// How long a model has to be fed before its output means anything. The number matters because
// anything that fades a freshly loaded model in — a crossfade between two captures of one device —
// hears the difference: a network still describing the silence it was born into is not the device.
//
// These read the config's SHAPE only, never its weights, which is what lets a real architecture be
// checked here without shipping a real model.

#include "../src/ReceptiveField.h"

#include <felitronics_test.h>

using felitronics::test::ok;
using felitronics::test::group;
using felitronics::nam::detail::receptiveFieldFromConfig;
using felitronics::nam::detail::receptiveFieldOfLayers;
using felitronics::nam::detail::partitionedTailSamples;
using felitronics::nam::detail::isRecurrent;
using felitronics::nam::detail::convNetField;

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
        // A config may carry both; the stack is the one that describes what the network does.
        nlohmann::json both = wavenet({ 3 }, { 4 });
        both["config"]["receptive_field"] = 99999;
        ok(receptiveFieldFromConfig(both) == 9, "layers WIN over a declared number when both are present");
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
        ok(receptiveFieldFromConfig({ { "architecture", "Linear" }, { "config", { { "receptive_field", "long" } } } }) == 0,
           "…and one that is not a number is refused rather than thrown on: this file runs inside"
           " prepareModel's catch-all, where a throw would REFUSE the load");
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
        // A WaveNet's dilations live inside `layers`; a top-level `dilations` key on a config that HAS
        // layers is a key NAM never looks at there, and reading it would double-count the stack.
        {
            nlohmann::json both = wavenet({ 2 }, { 1 });
            both["config"]["dilations"] = nlohmann::json::array({ 1000 });
            ok(receptiveFieldFromConfig(both) == 2,
               "…and a stray top-level `dilations` beside real `layers` is not read: NAM's WaveNet parser"
               " never looks there");
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
        // range is UNDEFINED BEHAVIOUR, and this repository's CI has a hard-fail UBSan job. A refusal
        // reads as absent — the same answer this file already gives a field spelled as a string.
        ok(receptiveFieldFromConfig({ { "architecture", "Linear" }, { "config", { { "receptive_field", 1e300 } } } }) == 0,
           "a field of 1e300 is refused, not cast: the cast is UB and the answer would be a number"
           " nobody can trust");
        {
            nlohmann::json hugeHead = { { "architecture", "WaveNet" },
                                        { "config", { { "layers", nlohmann::json::array({
                                              { { "kernel_size", 2 }, { "dilations", nlohmann::json::array({ 1 }) },
                                                { "head", { { "out_channels", 1 }, { "kernel_size", 1e300 },
                                                            { "bias", false } } } } }) } } } };
            ok(receptiveFieldFromConfig(hugeHead) == 2, "…and so is a head kernel of 1e300");
            nlohmann::json hugeConv = { { "architecture", "ConvNet" },
                                        { "config", { { "channels", 1 },
                                                      { "dilations", nlohmann::json::array({ 1e300, 4 }) } } } };
            ok(receptiveFieldFromConfig(hugeConv) == 5, "…and one dilation of 1e300 is skipped, not carried");
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
            // Refused reads as absent, and NAM's own answer then covers what it did build.
            nlohmann::json pastTheInt = { { "architecture", "WaveNet" },
                                          { "config", { { "layers", nlohmann::json::array({
                                                { { "kernel_size", 2 },
                                                  { "dilations", nlohmann::json::array({ 4294967396LL, 4 }) } } }) } } } };
            ok(receptiveFieldFromConfig(pastTheInt) == 5,
               "a dilation past the int is refused and its neighbour still counts: "
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

    return felitronics::test::report();
}

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

namespace {

// One dilated stack, written the way a .nam does.
nlohmann::json wavenet(std::vector<int> kernels, std::vector<int> dilations) {
    return { { "architecture", "WaveNet" },
             { "config", { { "layers", nlohmann::json::array({
                   { { "kernel_sizes", kernels }, { "dilations", dilations } } }) } } } };
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

    return felitronics::test::report();
}

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

    group("an architecture nothing here can read says so, rather than guessing");
    {
        ok(receptiveFieldFromConfig({ { "architecture", "Linear" }, { "config", { { "receptive_field", 3 } } } }) == 0,
           "no layers array -> 0, and the caller falls back to what NAM itself reports");
        ok(receptiveFieldFromConfig({ { "architecture", "WaveNet" } }) == 0, "no config at all -> 0");
        ok(receptiveFieldFromConfig({ { "architecture", "WaveNet" }, { "config", { { "layers", 7 } } } }) == 0,
           "…and a layers field that is not an array is refused, not indexed");
    }

    return felitronics::test::report();
}

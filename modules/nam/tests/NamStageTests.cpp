// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// Deterministic, file-free NAM backend tests. The tiny Linear model is an analytic FIR fixture and,
// because its parser is registered from linear.cpp, also guards the WHOLE_ARCHIVE link contract.

#include <felitronics_test.h>
#include <felitronics/nam/NamStage.h>

#include <namz.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

static std::atomic<long> g_allocs { 0 };
void* operator new      (std::size_t size) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (size ? size : 1); }
void* operator new[]    (std::size_t size) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (size ? size : 1); }
void  operator delete   (void* p) noexcept { std::free (p); }
void  operator delete[] (void* p) noexcept { std::free (p); }
void  operator delete   (void* p, std::size_t) noexcept { std::free (p); }
void  operator delete[] (void* p, std::size_t) noexcept { std::free (p); }

namespace
{
std::string firModel (const char* sampleRate = "48000", const char* metadata = "")
{
    std::string json =
        R"({"version":"0.5.0","architecture":"Linear","config":{"receptive_field":3,"bias":true,"implementation":"direct"},"weights":[0.5,-0.25,0.125,0.1])";
    if (sampleRate != nullptr)
        json += std::string (R"(,"sample_rate":)") + sampleRate;
    if (metadata[0] != '\0')
        json += std::string (R"(,"metadata":)") + metadata;
    json += '}';
    return json;
}

std::string gainModel (const char* sampleRate = "48000", const char* metadata = "")
{
    std::string json =
        R"({"version":"0.5.0","architecture":"Linear","config":{"receptive_field":1,"bias":false,"implementation":"direct"},"weights":[1.0])";
    if (sampleRate != nullptr)
        json += std::string (R"(,"sample_rate":)") + sampleRate;
    if (metadata[0] != '\0')
        json += std::string (R"(,"metadata":)") + metadata;
    json += '}';
    return json;
}

std::string scalarModel (int id, bool withLoudness)
{
    std::string json =
        R"({"version":"0.5.0","architecture":"Linear","config":{"receptive_field":1,"bias":false,"implementation":"direct"},"weights":[)"
        + std::to_string (0.01 * (double) id) + R"(],"sample_rate":48000)";
    if (withLoudness)
        json += R"(,"metadata":{"loudness":)" + std::to_string (-30.0 - (double) id) + '}';
    json += '}';
    return json;
}

void putU32 (std::vector<std::uint8_t>& bytes, std::uint32_t value)
{
    bytes.push_back ((std::uint8_t) value);
    bytes.push_back ((std::uint8_t) (value >> 8));
    bytes.push_back ((std::uint8_t) (value >> 16));
    bytes.push_back ((std::uint8_t) (value >> 24));
}

std::vector<std::uint8_t> namzV1GainModel()
{
    // Fixed v1 wire fixture: unlike namz::pack (which emits v2), v1 begins its body immediately
    // after the 8-byte header and has no display-metadata block or meta-length field.
    const std::string skeleton =
        R"({"version":"0.5.0","architecture":"Linear","config":{"receptive_field":1,"bias":false,"implementation":"direct"},"weights":0,"sample_rate":48000})";
    std::vector<std::uint8_t> bytes { 'N', 'A', 'M', 'Z', 1, 0, 0, 0 };
    putU32 (bytes, (std::uint32_t) skeleton.size());
    bytes.insert (bytes.end(), skeleton.begin(), skeleton.end());
    putU32 (bytes, 1);       // one weights array
    putU32 (bytes, 1);       // containing one float
    const float gain = 0.625f;
    const auto* payload = reinterpret_cast<const std::uint8_t*> (&gain);
    bytes.insert (bytes.end(), payload, payload + sizeof (gain));
    return bytes;
}

bool load (felitronics::nam::NamStage& stage, const std::string& json, float trimDb = 0.0f)
{
    return stage.loadModelFromMemory (json.data(), json.size(), trimDb);
}

std::vector<float> runMono (felitronics::nam::NamStage& stage, const std::vector<float>& input,
                            bool normalize)
{
    std::vector<float> output = input;
    float* io[1] { output.data() };
    stage.process (io, 1, (int) output.size(), normalize);
    return output;
}

bool allFinite (const std::vector<float>& values)
{
    for (float value : values)
        if (! std::isfinite (value))
            return false;
    return true;
}
} // namespace

int main()
{
    using namespace felitronics;
    std::printf ("felitronics::nam NamStage tests\n");

    test::group ("no-model passthrough");
    {
        nam::NamStage stage;
        stage.prepare (48000.0, 512);
        std::vector<float> left (512), right (512);
        for (int i = 0; i < 512; ++i)
        {
            left[(std::size_t) i] = 0.3f * std::sin (0.01f * (float) i);
            right[(std::size_t) i] = -left[(std::size_t) i];
        }
        const auto leftBefore = left;
        const auto rightBefore = right;
        float* io[2] { left.data(), right.data() };
        stage.process (io, 2, 512, true);
        test::ok (left == leftBefore && right == rightBefore, "buffer is untouched without a model");
        test::ok (stage.latencySamples() == 0, "no model reports zero latency");
    }

    test::group ("analytic Linear model impulse response");
    {
        nam::NamStage stage;
        stage.prepare (48000.0, 64);
        const auto json = firModel();
        test::ok (load (stage, json), "minimal Linear .nam loads from memory");
        test::ok (stage.hasModel(), "loaded model is live synchronously");
        std::vector<float> impulse (16, 0.0f);
        impulse[0] = 1.0f;
        const auto output = runMono (stage, impulse, false);
        test::approx (output[0], 0.6, 1.0e-6, "tap 0 plus bias");
        test::approx (output[1], -0.15, 1.0e-6, "tap 1 plus bias");
        test::approx (output[2], 0.225, 1.0e-6, "tap 2 plus bias");
        test::approx (output[8], 0.1, 1.0e-6, "tail is the configured bias");
    }

    test::group ("packed .namz is sample-identical to raw JSON");
    {
        const auto json = firModel();
        const auto packed = namz::pack (json.data(), json.size());
        test::ok (! packed.empty(), "namz::pack returns a packed model");

        nam::NamStage rawStage, packedStage;
        rawStage.prepare (48000.0, 128);
        packedStage.prepare (48000.0, 128);
        test::ok (load (rawStage, json), "raw model loads");
        test::ok (packedStage.loadModelFromMemory (packed.data(), packed.size()), "packed model loads");

        std::vector<float> input (128);
        for (int i = 0; i < 128; ++i)
            input[(std::size_t) i] = 0.4f * std::sin (0.07f * (float) i);
        const auto raw = runMono (rawStage, input, false);
        const auto zipped = runMono (packedStage, input, false);
        test::ok (raw == zipped, ".nam and .namz produce byte-identical samples");
    }

    test::group ("namz v1 wire compatibility");
    {
        const auto packed = namzV1GainModel();
        nam::NamStage stage;
        stage.prepare (48000.0, 32);
        test::ok (stage.loadModelFromMemory (packed.data(), packed.size()),
                  "fixed v1 fixture without a metadata block loads");
        const std::vector<float> input (16, 0.4f);
        const auto output = runMono (stage, input, false);
        test::approx (output[7], 0.25, 1.0e-7, "v1 payload weight reaches the NAM model unchanged");
    }

    test::group ("loudness makeup and per-model trim");
    {
        const auto json = gainModel ("48000", R"({"loudness":-20.0})");
        nam::NamStage rawStage, normalizedStage;
        rawStage.prepare (48000.0, 64);
        normalizedStage.prepare (48000.0, 64);
        test::ok (load (rawStage, json, 6.0f) && load (normalizedStage, json, 6.0f),
                  "tagged model loads with trim");
        test::ok (normalizedStage.modelHasLoudness(), "loudness tag is mirrored");
        test::approx (normalizedStage.modelLoudness(), -20.0, 1.0e-9, "loudness value is mirrored");

        const std::vector<float> input (64, 0.25f);
        const auto raw = runMono (rawStage, input, false);
        const auto normalized = runMono (normalizedStage, input, true);
        const double expected = std::pow (10.0, 8.0 / 20.0);   // (-18 - -20) + trimDb(6)
        test::approx (normalized[20] / raw[20], expected, 1.0e-5,
                      "-18 dB target makeup and trimDb fold into one gain");
    }

    test::group ("model-rate contract");
    {
        nam::NamStage stage;
        stage.prepare (48000.0, 128);
        const auto live = gainModel ("48000");
        test::ok (load (stage, live), "48 kHz baseline model is live before the refusal");
        const std::vector<float> probe { -0.3f, 0.2f, 0.75f, -0.125f };
        const auto beforeRefusal = runMono (stage, probe, false);

        const auto wrongRate = gainModel ("44100");
        test::ok (! load (stage, wrongRate), "44.1 kHz tagged model is refused after a 48 kHz prepare");
        test::ok (stage.hasModel() && stage.modelSampleRate() == 48000.0,
                  "refused load leaves the prior 48 kHz model and mirrors live");
        const auto afterRefusal = runMono (stage, probe, false);
        test::ok (beforeRefusal == afterRefusal,
                  "the previous model remains audibly byte-identical after a refused load");

        stage.clearModel();
        float advance = 0.0f;
        float* advanceIo[1] { &advance };
        stage.process (advanceIo, 1, 1, false);
        stage.collectGarbage();
        const auto untagged = gainModel (nullptr);
        test::ok (load (stage, untagged), "untagged model is accepted under the historical 48 kHz fallback");
        test::ok (stage.modelSampleRate() <= 0.0, "untagged model reports an unknown sample rate");
    }

    test::group ("bounded retire queue preserves accepted last-wins intents and live mirrors");
    {
        nam::NamStage stage;
        stage.prepare (48000.0, 16);
        bool allAccepted = true;
        for (int id = 1; id <= 72; ++id)
        {
            allAccepted = load (stage, scalarModel (id, (id % 2) == 0)) && allAccepted;
            test::ok (stage.modelSampleRate() == 48000.0,
                      "sample-rate mirror always describes the live model while loads are frozen");
            if (id >= 65)
                test::ok (! stage.modelHasLoudness(),
                          "parked tagged loads do not overwrite the untagged 65th live mirror");
        }
        test::ok (allAccepted, "more than 70 frozen-audio loads are accepted by the lossless contract");
        test::ok (stage.hasModel() && ! stage.modelHasLoudness(),
                  "after 64 retirements the 65th model remains live and later loads are parked");

        float probe = 1.0f;
        float* io[1] { &probe };
        stage.process (io, 1, 1, false);
        test::approx (probe, 0.65, 1.0e-6, "the frozen live model is audibly the 65th accepted load");
        test::ok (stage.collectGarbage(), "one audio block plus garbage collection lands the pending intent");
        test::ok (stage.modelHasLoudness() && stage.modelSampleRate() == 48000.0,
                  "latest pending model mirrors publish only when that model lands");
        test::approx (stage.modelLoudness(), -102.0, 1.0e-9,
                      "the latest (72nd), not the first parked load, wins");
        probe = 1.0f;
        stage.process (io, 1, 1, false);
        test::approx (probe, 0.72, 1.0e-6, "the latest pending model is audible after the drain");
        stage.collectGarbage();

        // The next load first reclaims the now-audio-safe prior model. Sixty-four frozen replacements fill the
        // bounded queue again; clearModel() must then park, retain live mirrors, and land after drain.
        for (int id = 100; id <= 163; ++id)
            allAccepted = load (stage, scalarModel (id, false)) && allAccepted;
        test::ok (allAccepted && stage.hasModel() && ! stage.modelHasLoudness(),
                  "retire queue is full again with an untagged live model");
        stage.clearModel();
        test::ok (stage.hasModel() && stage.modelSampleRate() == 48000.0 && ! stage.modelHasLoudness(),
                  "a deferred clear leaves the live model and mirrors truthful");
        probe = 0.0f;
        stage.process (io, 1, 1, false);
        test::ok (stage.collectGarbage(), "deferred clear lands after one block drains the queue");
        test::ok (! stage.hasModel() && stage.modelSampleRate() == 0.0
                  && ! stage.modelHasLoudness() && stage.modelLoudness() == 0.0,
                  "landed clear atomically publishes empty mirrors");
    }

    test::group ("true-stereo instances are independent and equal independent mono runs");
    {
        nam::NamStage stereo, monoLeft, monoRight;
        stereo.prepare (48000.0, 32);
        monoLeft.prepare (48000.0, 32);
        monoRight.prepare (48000.0, 32);
        const auto json = firModel();
        test::ok (load (stereo, json) && load (monoLeft, json) && load (monoRight, json),
                  "FIR fixture loads into stereo and independent mono stages");

        std::vector<float> left (16, 0.0f), right (16, 0.0f);
        left[0] = 1.0f;
        const auto leftInput = left;
        const auto rightInput = right;
        float* stereoIo[2] { left.data(), right.data() };
        stereo.process (stereoIo, 2, 16, false);
        const auto independentLeft = runMono (monoLeft, leftInput, false);
        const auto independentRight = runMono (monoRight, rightInput, false);

        bool biasOnly = true;
        for (float value : right)
            biasOnly = biasOnly && value == 0.1f;
        test::ok (biasOnly, "silent R receives exactly its own bias-only response with no L-tap crosstalk");
        test::ok (left == independentLeft && right == independentRight,
                  "one stereo run equals two independent mono runs sample-for-sample");
    }

    test::group ("96 kHz host rate matching");
    {
        nam::NamStage stage;
        stage.prepare (96000.0, 256);
        const auto json = gainModel();
        test::ok (load (stage, json), "48 kHz model loads on a 96 kHz host");
        test::ok (stage.latencySamples() == 9,
                  "latency pins ceil(3 * 96000 / 48000) + 3 exactly");

        std::vector<float> block (256);
        bool finite = true;
        float peak = 0.0f;
        for (int pass = 0; pass < 20; ++pass)
        {
            for (int i = 0; i < 256; ++i)
                block[(std::size_t) i] = 0.25f * std::sin (0.03f * (float) (pass * 256 + i));
            float* io[1] { block.data() };
            stage.process (io, 1, 256, false);
            finite = finite && allFinite (block);
            for (float value : block) peak = std::max (peak, std::fabs (value));
        }
        test::ok (finite, "resampled processing stays finite");
        test::ok (peak > 0.05f && peak < 1.0f, "resampled unity model output remains sane");

        auto processLayout = [&stage] (int channels)
        {
            std::vector<float> left (256, 0.1f), right (256, -0.1f);
            float* io[2] { left.data(), right.data() };
            stage.process (io, channels, 256, false);
            return allFinite (left) && (channels == 1 || allFinite (right));
        };
        const bool monoFirst = processLayout (1);
        stage.prepare (96000.0, 256);
        const bool stereo = processLayout (2);
        const int stereoLatency = stage.latencySamples();
        stage.prepare (96000.0, 256);
        const bool monoAgain = processLayout (1);
        test::ok (monoFirst && stereo && monoAgain,
                  "mono-to-stereo-to-mono re-prepare stays finite in both resampler lanes");
        test::ok (stereoLatency == 9 && stage.latencySamples() == 9,
                  "layout re-prepare leaves the pinned host latency stable");
    }

    test::group ("process is RT no-alloc");
    {
        nam::NamStage stage;
        stage.prepare (48000.0, 512);
        const auto json = gainModel();
        test::ok (load (stage, json), "RT fixture model loads");
        std::vector<float> left (512, 0.2f), right (512, -0.15f);
        float* io[2] { left.data(), right.data() };
        stage.process (io, 2, 512, false);        // warm every process-reachable container first
        const long before = g_allocs.load (std::memory_order_relaxed);
        stage.process (io, 2, 512, false);
        stage.process (io, 2, 512, true);
        test::okNoAlloc (g_allocs.load (std::memory_order_relaxed) == before,
                         "NamStage::process performs no heap allocation");
    }

    return test::report();
}

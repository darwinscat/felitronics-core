// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// Deterministic, file-free NAM backend tests. The tiny Linear model is an analytic FIR fixture and,
// because its parser is registered from linear.cpp, also guards the WHOLE_ARCHIVE link contract.

#include <felitronics_test.h>
#include <felitronics/nam/NamStage.h>
#include <felitronics/core/StreamResampler.h>   // the delay geometry is ASKED of the class, never restated

#include <namz.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <string>
#include <thread>
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
constexpr double kPi = 3.14159265358979323846;   // house convention: M_PI is not portable (MSVC)

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

// 🔴 THE ONE STATEFUL CAPTURE IN THIS FILE, AND IT EXISTS BECAUSE EVERYTHING ELSE HERE IS NOT. Every
// other model below is a `Linear` — an FIR whose state after any prewarm on silence is exactly zeros —
// so no fixture built on them can see how LONG the prewarm ran. That length is decided by the block
// handed to NAM's `Reset`, which is `maxModelFrames`, so a whole class of sizing changes was invisible
// to this suite: a mutation round ran it twice against a survivor and twice reported "equivalent",
// which was WRONG and only became visible once this fixture existed.
//
// One LSTM cell, arranged so its state on silence neither settles nor explodes:
//   W = 0 (4x2), b = [i,f,g,o] = [0, 10, 0, 0], h0 = 0, c0 = 1, head weight 1, head bias 0.
// With W = 0 and silence, `ifgo` IS `b`, so the recurrence collapses to c <- sigmoid(10)*c and
// h = sigmoid(0)*tanh(c) = 0.5*tanh(c): a decay with a time constant of ~20000 samples, which is the
// same order as the half second NAM prewarms an LSTM for. The first output after the prewarm is
// therefore 0.5*tanh(sigmoid(10)^N) with N the prewarm length — a direct readout of the Reset block.
//
// PORTABILITY, because a pinned value has to survive three libms: the sigmoid's argument is CONSTANT,
// so the state is one libm result raised to an integer power by repeated multiplication. A one-ulp
// difference in that single `exp` drifts the product by ~1e-11 over 24000 steps, four orders below the
// ~1e-3 the mutations here move — hence the 1e-6 tolerances used against these pins.
std::string lstmModel (const char* sampleRate = "48000")
{
    return std::string (R"({"version":"0.5.0","architecture":"LSTM","config":{"num_layers":1,)"
                        R"("input_size":1,"hidden_size":1},)"
                        R"("weights":[0,0,0,0,0,0,0,0,0,10,0,0,0,1,1,0],"sample_rate":)")
           + sampleRate + "}";
}

// 🔴 AND THE SECOND ORACLE, OF A DIFFERENT CONSTRUCTION — the house rule, and here it is not decoration.
// The group below pins the model's outputs as LITERALS, which is the right shape for a value that must
// not move when the code does; this one COMPUTES what the same silence should leave behind, from the
// two facts the fixture is built on (sigmoid(10) per sample, 0.5*tanh at the head) and a prewarm length
// the caller states. A pin catches a number that moved; this catches a number that moved for a reason
// the pin's author did not have in mind — and a pre-merge round wrote it precisely because the pins
// alone let three mutants through.
float stateAfterSilence (int samples)
{
    const float forget = 1.0f / (1.0f + std::exp (-10.0f));
    float cell = 1.0f;
    for (int i = 0; i < samples; ++i) cell *= forget;
    return 0.5f * std::tanh (cell);
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

// A DENSE Linear capture, deterministic: an LCG, one step per tap. Dense is the point — NAM runs a
// Linear capture past 256 taps through a partitioned-FFT engine (`implementation` defaults to `auto`),
// whose ring holds input SPECTRA past the field, and a SPARSE kernel cannot see what that leaves
// behind: measured, this kernel drained by exactly its field left 1.909e-08 on 46 samples of the
// return, and a single-tap kernel at the same length left nothing.
std::string denseLinearModel (int taps, const char* implementation)
{
    std::string w = "[";
    std::uint32_t state = 7u;
    char buf[32];
    for (int i = 0; i < taps; ++i)
    {
        state = state * 1664525u + 1013904223u;
        const double v = ((double) (state >> 8) / 16777216.0 - 0.5) * 2.0 / std::sqrt ((double) taps);
        std::snprintf (buf, sizeof buf, "%.7g", v);
        w += buf;
        if (i + 1 < taps) w += ',';
    }
    w += ']';
    std::string json = R"({"version":"0.5.0","architecture":"Linear","config":{"receptive_field":)"
                     + std::to_string (taps) + R"(,"bias":false)";
    if (implementation != nullptr) json += std::string (R"(,"implementation":")") + implementation + R"(")";
    return json + R"(},"weights":)" + w + R"(,"sample_rate":48000})";
}

// A one-cell LSTM with a SLOW forget gate: sigmoid(10) = 0.99995, so the cell's time constant is about
// 22 000 samples. Weight order is the one lstmModel in the RT-alloc suite already relies on — W is
// 4x2 row-major over gates i,f,g,o and columns (x, h), then b[i,f,g,o], then h0, c0, then the head
// weight and bias. `tagged` decides whether the model carries a `sample_rate`, which is the whole
// point: NAM sizes an LSTM's prewarm as 0.5 x the TAG, and an untagged one therefore answers 1.
std::string slowLstmModel (bool tagged)
{
    std::string json = R"({"version":"0.5.0","architecture":"LSTM","config":{"num_layers":1,)"
                       R"("input_size":1,"hidden_size":1},"weights":[3,0, 0,0, 0.002,0, 0,0,)"
                       R"(  0,10,0,0,  0,0, 1,0])";
    if (tagged) json += R"(,"sample_rate":48000)";
    return json + "}";
}

// A Linear capture that hands back what came in `d` samples ago: `d` zero taps and then a one. It is
// the cheapest fixture with a LONG memory, and it is also the one NAM answers zero about — the field is
// a plain number in the config and nothing but detail::declaredReceptiveField reads it.
std::string delayModel (int d)
{
    std::string w = "[";
    for (int i = 0; i < d; ++i) w += "0.0,";
    w += "1.0]";
    return R"({"version":"0.5.0","architecture":"Linear","config":{"receptive_field":)" + std::to_string (d + 1)
         + R"(,"bias":false,"implementation":"direct"},"weights":)" + w + R"(,"sample_rate":48000})";
}

// A REAL WaveNet whose memory is as long as it says. One one-channel layer, kernel 2, one dilation —
// the weight order is rechannel; dilated conv (OLDEST tap, then newest) + bias; condition mixin;
// residual 1x1 + bias; head rechannel; head scale. Putting the 1 on the OLDEST tap is what makes the
// network actually reach back `dilation` samples: the shipped [1,0,0,…] fixture has BOTH convolution
// taps at zero, so it declares a field of thousands and forgets after one sample. Measured: with
// [1,0,…] the impulse response's last non-zero sample is 0; with [1,1,0,…] it is `dilation`.
std::string waveNetDelayModel (int dilation)
{
    return std::string (R"({"version":"0.5.0","architecture":"WaveNet","config":{"layers":[{"input_size":1,)")
         + R"("condition_size":1,"head_size":1,"head_bias":false,"channels":1,"kernel_size":2,"dilations":[)"
         + std::to_string (dilation)
         + R"(],"activation":"Tanh","gated":false}],"head_scale":1.0},"weights":[1,1,0,0,1,0,0,1,1],)"
         + R"("sample_rate":48000})";
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
    felitronics::test::run (stage.process (io, 1, (int) output.size(), normalize));
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
        felitronics::test::run (stage.process (io, 2, 512, true));
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
        felitronics::test::run (stage.process (advanceIo, 1, 1, false));
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
        felitronics::test::run (stage.process (io, 1, 1, false));
        test::approx (probe, 0.65, 1.0e-6, "the frozen live model is audibly the 65th accepted load");
        test::ok (stage.collectGarbage(), "one audio block plus garbage collection lands the pending intent");
        test::ok (stage.modelHasLoudness() && stage.modelSampleRate() == 48000.0,
                  "latest pending model mirrors publish only when that model lands");
        test::approx (stage.modelLoudness(), -102.0, 1.0e-9,
                      "the latest (72nd), not the first parked load, wins");
        probe = 1.0f;
        felitronics::test::run (stage.process (io, 1, 1, false));
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
        felitronics::test::run (stage.process (io, 1, 1, false));
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
        felitronics::test::run (stereo.process (stereoIo, 2, 16, false));
        const auto independentLeft = runMono (monoLeft, leftInput, false);
        const auto independentRight = runMono (monoRight, rightInput, false);

        bool biasOnly = true;
        for (float value : right)
            biasOnly = biasOnly && value == 0.1f;
        test::ok (biasOnly, "silent R receives exactly its own bias-only response with no L-tap crosstalk");
        test::ok (left == independentLeft && right == independentRight,
                  "one stereo run equals two independent mono runs sample-for-sample");
    }

    // LAW 11(a) — the length is a CAPACITY. This used to be `n = std::min (numSamples, maxBlock)`, so
    // everything past the prepared block came out of the amp BIT-IDENTICAL to its input: measured on
    // this same Linear FIR, prepared for 64 and called with 512, 448 of 512 samples never met the model.
    // Simply dropping the clamp would have been worse than the defect — processChannel copies n samples
    // into a `maxModelFrames` scratch and NAM's own buffers are sized from the prepared block, so an
    // unclamped n is a heap overflow here and a resize (an allocation) inside NAM, on the audio thread.
    test::group ("law 11a: a call longer than maxBlock is CHUNKED, not clamped");
    {
        const int MB = 64, N = 517;                       // 517 is deliberately not a multiple of 64
        nam::NamStage one, many;
        one.prepare (48000.0, MB);
        many.prepare (48000.0, MB);
        const auto json = firModel();
        test::ok (load (one, json) && load (many, json), "precondition: the FIR fixture loads into both");

        std::vector<float> a ((std::size_t) N, 0.0f), b;
        for (int i = 0; i < N; ++i) a[(std::size_t) i] = (i % 41 == 0) ? 0.9f : 0.05f * (float) std::sin (0.037 * i);
        const auto in = a;
        b = a;

        float* pa[1] { a.data() };
        felitronics::test::run (one.process (pa, 1, N, false));
        for (int off = 0; off < N; )
        {
            const int m = std::min (N - off, MB);
            float* pb[1] { b.data() + off };
            felitronics::test::run (many.process (pb, 1, m, false));
            off += m;
        }

        int touched = 0; double spread = 0.0;
        for (int i = 0; i < N; ++i)
        {
            if (a[(std::size_t) i] != in[(std::size_t) i]) ++touched;
            spread = std::max (spread, (double) std::fabs (a[(std::size_t) i] - in[(std::size_t) i]));
        }
        test::ok (spread > 1e-3, "precondition: the model MOVED the signal (max |out-in| = " + std::to_string (spread) + ")");
        test::ok (touched == N, "every one of the " + std::to_string (N) + " samples met the model (was: 64 of 517)");
        test::ok (a == b, "...and one long call is bit-identical to the caller's own 64-sample calls");
        test::ok (! one.process (pa, 3, N, false), "a 3-channel call is REFUSED — a NAM capture is mono or true-stereo");
        test::ok (! one.process (pa, 1, -1, false), "a negative length is REFUSED");
        test::ok (one.process (pa, 0, N, false) && one.process (pa, 1, 0, false), "degenerate calls are accepted no-ops");
    }

    test::group ("latency is a MEASUREMENT, not a formula — the reported number matches the real delay");
    {
        // The old `ceil(3*hostSR/modelRunSR) + 3` was a guess and it was 2.16 samples long at 44.1 kHz;
        // the tests pinned the FORMULA, so nothing caught it. This one measures the delay the stage
        // ACTUALLY has and compares. Geometry, asked of the class rather than restated here (restating
        // it is how the last one went stale): reset() leaves kTaps leading zeros with pos = kHalf, so
        // each stage delays by exactly D = StreamResampler::delayInputSamples() of its OWN input
        // samples -> the round trip is D·(1 + hostSR/modelRunSR) host samples. With the P34 sinc kernel
        // D = 32, i.e. 61.4 at 44.1 kHz and 96.0 at 96 kHz; with the cubic it was D = 2 and 3.8375.
        //
        // TWO instruments, because each is blind where the other sees:
        //  * PHASE gives the fraction exactly but is periodic — it cannot tell a delay from that delay
        //    plus a whole tone period. A crew round proved that is not theoretical: a mutation that
        //    inserted a real 512-sample FIFO in front of the resampling branch left this test printing
        //    3.8375, because the tone period was exactly 512 samples.
        //  * ONSET (a single impulse) gives the integer unambiguously but is a poor fraction oracle
        //    through a DECIMATING stage — an impulse is not DC and the stage does not preserve it
        //    (measured DC sums 0.999 / 0.743 / 2.000 at 44.1 / 88.2 / 96 kHz).
        // So the onset resolves which period the phase belongs to, and the phase supplies the fraction.
        //
        // 🔴 THE WINDOW IS ITSELF A GRID and the first version of this test tripped on that too: with a
        // window that is not a whole number of tone periods, the discarded sum(sin(2Wn-WD)) term does
        // not cancel and reads as ~0.1 samples of phantom delay (6.11 where the geometry says 6.00).
        // So: f = hostSR/512 makes the tone period EXACTLY 512 samples at any host rate, and the window
        // is 150528 = 512*294 = 147*1024 samples — a whole number of tone periods AND of the
        // resampler's own 147-sample modulation period, both also whole multiples of the 64 block.
        auto phaseDelay = [] (double hostSR, int block, int channels)
        {
            nam::NamStage stage;
            stage.prepare (hostSR, block);
            const auto json = gainModel();
            if (! load (stage, json)) return -1.0e9;
            stage.prepare (hostSR, block);
            const double f = hostSR / 512.0;                       // period = 512 samples, exactly
            const double W = 2.0 * 3.14159265358979323846 * f / hostSR;
            const int skip = 76800, span = 150528;                 // 512*294 and 147*1024
            double sc = 0.0, ss = 0.0;
            std::vector<float> l ((std::size_t) block), r ((std::size_t) block);
            const int probe = channels - 1;                        // measure the LAST lane, so a stereo
            for (int off = 0; off < skip + span; off += block)     // run cannot pass on lane 0 alone
            {
                for (int i = 0; i < block; ++i)
                {
                    l[(std::size_t) i] = (float) std::sin (W * (off + i));
                    r[(std::size_t) i] = (float) std::sin (W * (off + i));
                }
                float* io[2] { l.data(), r.data() };
                felitronics::test::run (stage.process (io, channels, block, false));
                const float* probed = (probe == 0) ? l.data() : r.data();
                if (off >= skip)
                    for (int i = 0; i < block; ++i)
                    {
                        sc += (double) probed[(std::size_t) i] * std::cos (W * (off + i));
                        ss += (double) probed[(std::size_t) i] * std::sin (W * (off + i));
                    }
            }
            return std::atan2 (-sc, ss) / W;                       // in [-256, 256) host samples
        };

        // ONSET: where a single impulse comes out. Unambiguous integer, no periodicity to alias.
        auto onsetDelay = [] (double hostSR, int block)
        {
            nam::NamStage stage;
            stage.prepare (hostSR, block);
            const auto json = gainModel();
            if (! load (stage, json)) return -1;
            stage.prepare (hostSR, block);
            const int at = 4000;
            int peak = -1; double pv = 0.0; int idx = 0;
            std::vector<float> b ((std::size_t) block);
            for (int off = 0; off < at + 4096; off += block)
            {
                for (int i = 0; i < block; ++i) b[(std::size_t) i] = ((off + i) == at) ? 1.0f : 0.0f;
                float* io[1] { b.data() };
                felitronics::test::run (stage.process (io, 1, block, false));
                for (int i = 0; i < block; ++i, ++idx)
                {
                    const double v = std::fabs ((double) b[(std::size_t) i]);
                    if (v > pv && idx > at - 50) { pv = v; peak = idx; }
                }
            }
            return peak - at;
        };

        // precondition: the instruments read ZERO where there is no resampler at all (NamStage.cpp
        // engages it only when |hostSR - modelRunSR| > 0.5), so neither can be reading its own bias.
        test::approx (phaseDelay (48000.0, 64, 1), 0.0, 0.05,
                      "precondition: at the model's own rate the measured phase delay is 0.00 samples");
        test::ok (onsetDelay (48000.0, 64) == 0,
                  "precondition: and the impulse comes straight back out, 0 samples late");

        // Every expectation below is ASKED of the one function that owns the composition, so the table
        // cannot drift away from the kernel the way the hand-written 3.8375 did — and, unlike the
        // previous version of these two lines, it cannot drift away from the SHAPE of the composition
        // either. Writing `kD + kD*host/48000` here was itself a restatement: it would keep passing
        // after a change that made the two legs different, which is exactly what the open kTaps-scaling
        // item does.
        //
        // 🔴 THIS IS NOT THE TEST COMPARING THE CODE WITH ITSELF. What is asserted below is that the
        // number the stage REPORTS matches a delay measured out of the SIGNAL — carrier phase, with the
        // whole-period ambiguity resolved by an impulse onset. The helper only supplies the expectation
        // for that independent measurement; if it lied, the signal would disagree with it.
        // 🔴 FRACTIONAL, and the distinction is not pedantry — it cost a red suite while this was
        // being written. The delay measured out of the SIGNAL is fractional (61.4000 at 44.1 kHz);
        // the number the stage REPORTS is that value rounded (61). They are two different quantities
        // and the tests below need both: the fractional one to compare against the measurement, the
        // rounded one to compare against the report. Asking the wrong one of the two functions is a
        // restatement in a new costume, so each is asked where it belongs.
        auto geoAt = [] (double host)
        {
            return felitronics::core::StreamResampler::pairDelayHostSamples (host, 48000.0);
        };

        struct Case { double host; };
        // 32 kHz is in the list ON PURPOSE: it is the row where lround and ceil disagree most visibly
        // (geometry 53.3333 -> 53 against 54). Without it the choice of rounding is untested, which a
        // crew mutation proved by swapping lround for ceil and surviving.
        // 🔴 THE ONLY NON-TAUTOLOGICAL ORACLE IN THIS FILE, and the list grew for that reason. Since
        // the composition became one function, every test that compares the stage's report against
        // "the geometry" is comparing that function with itself and is green by construction. What is
        // NOT tautological is this: a delay measured out of the SIGNAL — carrier phase for the
        // fraction, impulse onset to resolve which whole period it belongs to — against the number the
        // stage reports. A wrong body in the one function would be invisible everywhere else and
        // visible here. So the measured list carries the rates that matter rather than four of them.
        for (const Case cc : { Case { 44100.0 }, Case { 96000.0 }, Case { 88200.0 }, Case { 32000.0 },
                               Case { 22050.0 }, Case { 64000.0 }, Case { 176400.0 }, Case { 192000.0 } })
        {
            const struct { double host, expect; } c { cc.host, geoAt (cc.host) };
            const int on = onsetDelay (c.host, 64);
            const double ph = phaseDelay (c.host, 64, 1);
            const double per = 512.0;
            double d = ph;                                          // resolve the period with the onset
            while (d < (double) on - per * 0.5) d += per;
            while (d > (double) on + per * 0.5) d -= per;
            std::printf ("      host %7.0f: onset +%d, phase-resolved delay %.4f samples (geometry %.4f)\n",
                         c.host, on, d, c.expect);
            test::ok (std::abs ((double) on - std::floor (c.expect + 0.5)) <= 1.0,
                      "host " + std::to_string ((int) c.host) + ": the IMPULSE comes out where the geometry "
                      "says, so no whole periods are hiding in the phase reading");
            test::approx (d, c.expect, 0.05,
                          "host " + std::to_string ((int) c.host) + ": the delay IS D*(1 + host/model)");
            nam::NamStage st;
            st.prepare (c.host, 64);
            const auto json = gainModel();
            test::ok (load (st, json), "model loads for the reported-latency comparison");
            st.prepare (c.host, 64);
            test::ok (std::fabs ((double) st.latencySamples() - d) <= 0.5,
                      "host " + std::to_string ((int) c.host) + ": latencySamples() = "
                      + std::to_string (st.latencySamples()) + " is the NEAREST integer to the measured "
                      + std::to_string (d) + " (ceil would report "
                      + std::to_string ((int) std::ceil (c.expect)) + ")");
        }

        // The right lane must rate-match too: a mutation that left instance 1's resamplers at the
        // identity ratio passed everything, because nothing measured the stereo delay. Phase ALONE is
        // not enough there either — it is periodic, so an extra whole tone period on the right lane
        // would be invisible exactly as it was on the left. Both instruments, both lanes.
        test::approx (phaseDelay (44100.0, 64, 2), geoAt (44100.0), 0.05,
                      "the RIGHT channel of a stereo call has the same measured phase delay as the left");
        {
            nam::NamStage st;
            st.prepare (44100.0, 64);
            const auto json = gainModel();
            test::ok (load (st, json), "model loads for the stereo onset");
            st.prepare (44100.0, 64);
            const int at = 4000;
            int peakR = -1; double pv = 0.0; int idx = 0;
            std::vector<float> l (64, 0.0f), r (64, 0.0f);
            for (int off = 0; off < at + 4096; off += 64)
            {
                for (int i = 0; i < 64; ++i)
                { l[(std::size_t) i] = 0.0f; r[(std::size_t) i] = ((off + i) == at) ? 1.0f : 0.0f; }
                float* io[2] { l.data(), r.data() };
                felitronics::test::run (st.process (io, 2, 64, false));
                for (int i = 0; i < 64; ++i, ++idx)
                {
                    const double v = std::fabs ((double) r[(std::size_t) i]);
                    if (v > pv && idx > at - 50) { pv = v; peakR = idx; }
                }
            }
            test::ok ((double) (peakR - at) == std::floor (geoAt (44100.0) + 0.5),
                      "…and the RIGHT lane's impulse comes out where the geometry says (+"
                      + std::to_string (peakR - at) + "), so no whole periods hide there either");
        }

        // Rounding is asserted as an INVARIANT over a sweep, not as a list of rates, because the one
        // case a list always misses is the exact half — and 🔴 WHERE THAT HALF LIVES MOVED WITH THE
        // KERNEL. With D = 2 the halves sat at h = 24000k − 36000 (12000, 36000, 60000, 84000 …); with
        // D = 32 they sit at h = 750·(2k+1) (750, 2250, … 44250, 45750 …), and 60000 now gives exactly
        // 72.0 — a whole number, which is to say the old exact-half case stopped being one. A sweep that
        // simply kept its rate list would have gone quietly blind, so 44250 is in the list on purpose.
        // Reporting an integer for a fractional delay costs at most half a sample; the old guessed
        // formula cost up to 3.3.
        {
            double worst = 0.0; double worstAt = 0.0;
            for (const double host : { 8000.0, 11025.0, 12000.0, 16000.0, 22050.0, 24000.0, 32000.0,
                                       36000.0, 44100.0, 44250.0, 47999.0, 48001.0, 60000.0, 64000.0,
                                       84000.0, 88200.0, 96000.0, 176400.0, 192000.0 })
            {
                nam::NamStage st;
                st.prepare (host, 64);
                const auto json = gainModel();
                if (! load (st, json)) { test::ok (false, "model loads at every swept rate"); break; }
                st.prepare (host, 64);
                // 🔴 THE ONE PLACE A RESTATEMENT IS THE POINT, and a crew round proved it by injecting
                // a +1 error into the geometry that this line — when it asked core — waved through
                // while the OLD tests caught it. An oracle's whole job is to disagree with the thing it
                // measures, so it must not share its arithmetic. Everything else in this file asks;
                // this computes, independently, from the two facts the header states: every stage
                // delays kHalf of its own input samples, and the return leg is converted at h/m.
                const double kD  = (double) felitronics::core::StreamResampler::kHalf;
                const double geo = kD + kD * host / 48000.0;
                const double err = std::fabs ((double) st.latencySamples() - geo);
                if (err > worst) { worst = err; worstAt = host; }
            }
            std::printf ("      worst reported-vs-geometry error over 19 host rates: %.4f samples (at %.0f Hz)\n",
                         worst, worstAt);
            // ⚠️ WHAT THIS LINE IS AND IS NOT — and this paragraph was itself WRONG for one commit,
            // which is the reason it now names its own history. It used to say both sides came from
            // the same function and that the check was therefore "a statement about ROUNDING ONLY".
            // That was true of a version in which `geo` ASKED core for the fractional value; a crew
            // round proved what it cost by injecting a +1 error into the geometry, which that version
            // waved through while the older tests caught it. `geo` was restored to an independent
            // computation — the two facts the header states, nothing shared with the thing measured —
            // and the paragraph explaining why it need not be independent was left standing behind it.
            //
            // What it is NOW: an independent oracle, so it catches BOTH a wrong composition and a
            // wrong rounding step. Measured at the close of that round, this line among them: a +1
            // error in the geometry fails 43 checks across the suite at nam=ON pffft=ON. (The commit
            // that restored this oracle published "21", which was already wrong when written — three
            // independent runs put it at 31 on that tree. A count is a measurement and goes stale like
            // any other, so it carries its configuration and its date of measurement or it carries
            // nothing.)
            test::ok (worst <= 0.5 + 1e-9,
                      "over 19 host rates including every exact-half case, the REPORTED INTEGER is the "
                      "nearest one to the fractional geometry (worst " + std::to_string (worst)
                      + ") — a claim about the rounding step, not about the formula");
            // …and at an exact half the BOUND accepts either neighbour, so pin the RULE itself: lround
            // takes halves away from zero. 44250 Hz gives exactly 61.5 with D = 32 (32 + 32·44250/48000
            // = 32 + 29.5). The rate that used to serve here, 60 kHz, now gives a whole 72.0.
            {
                nam::NamStage half;
                half.prepare (44250.0, 64);
                const auto json = gainModel();
                test::ok (load (half, json), "model loads at the exact-half rate");
                half.prepare (44250.0, 64);
                test::ok (half.latencySamples() == 62,
                          "at 44250 Hz the geometry is exactly 61.5 and the reported number is 62 — halves "
                          "go away from zero, which the <=0.5 bound alone would not pin");
            }
        }

        // The gate itself: NamStage.cpp engages the resampler only past |hostSR - modelRunSR| > 0.5,
        // and nothing tested either side of that edge — a mutation widening it to 10 Hz survived.
        {
            nam::NamStage near, past;
            near.prepare (48000.4, 64); past.prepare (48001.0, 64);
            const auto json = gainModel();
            test::ok (load (near, json) && load (past, json), "models load either side of the resampling gate");
            near.prepare (48000.4, 64); past.prepare (48001.0, 64);
            test::ok (near.latencySamples() == 0,
                      "0.4 Hz off the model rate is INSIDE the gate: no resampler, no latency");
            // 🔴 AND 0.6 — the OUTSIDE edge, a hair past it. This cell exists because a crew round
            // widened the gate from `> 0.5` all the way to `> 0.9` and the mutant SURVIVED: the only
            // outside pin was 48001.0, a whole hertz out, so every threshold between 0.5 and 1.0 was
            // free. The paragraph below already claimed 48000.6 was in the fixture — it was not, and a
            // comment naming a value the code does not have is how a hole stays open.
            {
                nam::NamStage past2;
                past2.prepare (48000.6, 64);
                const auto j3 = gainModel();
                test::ok (load (past2, j3), "a model loads on a stage half a hertz past the gate");
                past2.prepare (48000.6, 64);
                test::ok (past2.latencySamples() == 64,
                          "0.6 Hz off the model rate is OUTSIDE the gate: a resampler is installed and "
                          "costs " + std::to_string (past2.latencySamples()) + " samples. Widening the "
                          "threshold anywhere into (0.5, 0.6] now fails here");
            }
            // 🔴 …AND A HAIR ABOVE THE BOUNDARY, not merely near it. 0.6 leaves the whole interval
            // (0.5, 0.6) free, and a later round widened the gate to 0.55 and walked through. The fix
            // is not another arbitrary distance: pin the predicate IMMEDIATELY above the edge, and
            // every `> 0.5 + eps` mutation with eps >= 1e-4 dies at once.
            {
                nam::NamStage hair;
                hair.prepare (48000.5001, 64);
                const auto j4 = gainModel();
                test::ok (load (hair, j4), "a model loads a ten-thousandth of a hertz past the gate");
                hair.prepare (48000.5001, 64);
                test::ok (hair.latencySamples() == 64,
                          "48000.5001 is OUTSIDE the gate and costs "
                          + std::to_string (hair.latencySamples()) + " samples, while 48000.5 exactly "
                          "costs none — the predicate is `> 0.5` and nothing wider");
            }
            // 🔴 EXACTLY 0.5 — the boundary itself, which nothing tested. A mutation flipping `> 0.5`
            // to `>= 0.5` survived the whole suite because 48000.4 and 48001.0 stood well clear of the
            // edge without standing on it, and the gate's own predicate is only visible AT it.
            {
                nam::NamStage edge;
                edge.prepare (48000.5, 64);
                const auto j2 = gainModel();
                test::ok (load (edge, j2), "model loads exactly half a hertz off the model rate");
                edge.prepare (48000.5, 64);
                test::ok (edge.latencySamples() == 0,
                          "EXACTLY 0.5 Hz off is INSIDE the gate — the predicate is strict `> 0.5`, so "
                          "the boundary belongs to the no-resampler side; `>= 0.5` would report "
                          + std::to_string (nam::NamStage::rateMatch (48000.5, 48000.0).latencySamples));
            }

            // …and the NORMALISATION, which nothing tested either: two mutations survived here, one
            // moving the default rate an untagged model runs at, one asking the gate about the raw
            // reported rate instead of the normalised one. Both are invisible until an UNTAGGED model
            // is asked what it costs, because a tagged one normalises to itself.
            {
                const auto untagged = gainModel (nullptr);          // no sample_rate field at all
                nam::NamStage a48, a441;
                a48.prepare (48000.0, 64);   test::ok (load (a48, untagged), "untagged model loads at 48 kHz");
                a48.prepare (48000.0, 64);
                a441.prepare (44100.0, 64);  test::ok (load (a441, untagged), "…and at 44.1 kHz");
                a441.prepare (44100.0, 64);
                test::ok (a48.latencySamples() == 0,
                          "an UNTAGGED model runs at the factory rate, so at a 48 kHz host it is not "
                          "resampled at all — which is only true if the rate is normalised BEFORE the "
                          "gate sees it. Ask the raw -1 instead and this reports "
                          + std::to_string (nam::NamStage::rateMatch (48000.0, 48000.0).latencySamples));
                test::ok (a441.latencySamples() == nam::NamStage::rateMatch (44100.0, nam::NamStage::kModelSampleRate).latencySamples
                          && a441.latencySamples() > 0,
                          "…and at 44.1 kHz it costs exactly what a model tagged at the factory rate "
                          "costs (" + std::to_string (a441.latencySamples()) + ") — which pins WHICH "
                          "rate the default is, not merely that there is one");
            }

            test::ok (past.latencySamples() == 64,
                      "1.0 Hz off it is OUTSIDE: the resampler engages and reports its full 64 samples — "
                      "D·(1 + 48001/48000) rounds to 64, and the DISCONTINUITY at the gate is the point: "
                      "0.4 Hz costs nothing and 1.0 Hz costs 64 samples of PDC");
        }
    }

    test::group ("THE TWO FUNCTIONS THEMSELVES — pinned directly, because nothing else varies them");
    {
        // 🔴 WHY THIS GROUP EXISTS. A diverse-testing round replaced `modelRunSR` inside
        // pairDelayHostSamples with the literal 48000 and the mutant survived every suite in the
        // repository — because every call site in the tree passes 48000, and every model that goes live
        // in practice runs at 48000. (Not "every model that can": install() compares with a HALF-HERTZ
        // TOLERANCE, so a fresh stage accepts 48000.5 and prepare() then adopts it — the same
        // over-claim the header carried until it was corrected there. Correcting one copy of a
        // sentence and not the other is this branch's own subject matter.) The same round
        // showed rateMatch's normalisation path is only ever reached through an instance, so an
        // untagged rate never reaches it directly either.
        //
        // The whole point of these two functions is to be called by consumers this repository does not
        // contain, at rates it does not itself use. So they are pinned as PURE FUNCTIONS here, with
        // values computed by hand from the two documented facts — kHalf per leg, the return leg
        // converted at h/m — and not by asking the code.
        using felitronics::core::StreamResampler;
        struct G { double h, m, want; };
        for (const G g : { G { 96000.0,  32000.0, 128.0 },      // 32 + 32·3
                           G { 44100.0,  44100.0,  64.0 },      // identity ratio: still two legs
                           G { 48000.0,  96000.0,  48.0 },      // 32 + 32·0.5
                           G { 22050.0,  44100.0,  48.0 },      // …the same ratio, different rates
                           G { 192000.0, 48000.0, 160.0 },      // 32 + 32·4
                           G { 44100.0,  48000.0,  61.4 } })
            test::approx (StreamResampler::pairDelayHostSamples (g.h, g.m), g.want, 1e-9,
                          "pairDelayHostSamples(" + std::to_string ((int) g.h) + ", "
                          + std::to_string ((int) g.m) + ") = " + std::to_string (g.want));

        test::ok (StreamResampler::pairDelayHostSamples (96000.0, 32000.0)
                  != StreamResampler::pairDelayHostSamples (32000.0, 96000.0),
                  "…and the two arguments are NOT interchangeable — a swapped call at a consumer is a "
                  "different number, not a reciprocal one, which is why the order is in the name");

        // 🔴 THE DEFAULTS OF THE STRUCT ITSELF, which nothing else in the tree varies either. They are
        // not decoration: RateMatch gained member initialisers in the same commit that pinned these
        // functions, and that SILENTLY changed what `RateMatch r {}` means for every consumer —
        // modelRunSR read 0 before it and reads the factory rate after, and the type stopped being
        // trivially default-constructible. A crew round mutated the default back to 0 and the mutant
        // survived the whole repository. The new value is the right one — a rate-match that has not
        // been computed yet describes a stage running at the factory rate, not one running at 0 Hz —
        // but "right and unannounced and untested" is how the restatements this branch removes got in.
        {
            const nam::NamStage::RateMatch d {};
            test::ok (d.modelRunSR == nam::NamStage::kModelSampleRate && ! d.resampling
                          && d.latencySamples == 0,
                      "a default-constructed RateMatch reads {" + std::to_string (d.modelRunSR)
                      + ", false, 0} — the factory rate, not zero");
        }

        // rateMatch: the three facts. One of these cells an instance really cannot reach — a rate no
        // stage would accept — and the OTHER one it reaches every day, which an earlier version of this
        // comment got backwards while the paragraph eight lines above had it right. prepare() hands
        // configureRates() the RAW reported rate (NamStage.cpp, and that was the whole point of the fix
        // there), so an untagged model arrives here as modelSR = -1 and takes the normalisation path.
        // Two spellings of one fact in one group, already disagreeing — inside the group written to
        // stop exactly that.
        struct R { double h, m; double runSR; bool res; int lat; const char* what; };
        for (const R r : { R { 44100.0,     -1.0, 48000.0, true,  61, "unknown rate normalises to the factory one" },
                           R { 44100.0,      0.0, 48000.0, true,  61, "…and so does zero" },
                           R { 48000.0,     -1.0, 48000.0, false,  0, "…which then is NOT resampled at a 48 kHz host" },
                           R { 48000.0,  44100.0, 44100.0, true,  67, "a 44.1 kHz model would cost 67, whether or not a stage would take it" },
                           R { 96000.0,  96000.0, 96000.0, false,  0, "equal rates: no resampler, no latency" },
                           R { 48000.5,  48000.0, 48000.0, false,  0, "exactly half a hertz is inside the gate" } })
        {
            const auto got = nam::NamStage::rateMatch (r.h, r.m);
            test::ok (got.modelRunSR == r.runSR && got.resampling == r.res && got.latencySamples == r.lat,
                      std::string (r.what) + " — got {" + std::to_string (got.modelRunSR) + ", "
                      + (got.resampling ? "true" : "false") + ", " + std::to_string (got.latencySamples) + "}");
        }
    }

    test::group ("🔴 THE RATE CONTRACT IS A FIXED WINDOW — the reference cannot move (P38)");
    {
        // WHAT THE CONTRACT PROMISES, as the property and not as the threshold: in every reachable
        // state, a model this stage holds was tagged within kModelRateTolerance of kModelSampleRate —
        // the SAME window on the first load and on the ten-thousandth. It used to be a window around
        // the last ACCEPTED tag, which prepare() then adopted, so the window WALKED. Measured on the
        // base commit, half-hertz steps down, host 48000: {load} 1 step, {load, process} 1 step,
        // {load, prepare} 66 (stopped by the retire queue, not by rates), {load, process, prepare}
        // 5000 with no refusal at all — run rate 45500.0. The magnitude that settles it is the
        // excursion sup|runRate - kModelSampleRate|: UNBOUNDED before, kModelRateTolerance after.
        using nam::NamStage;

        // 1. The predicate itself, at both edges. By hand: the window is closed at both ends because
        //    the gate it delegates to refuses only PAST half a hertz, and a model that reports no rate
        //    normalises INTO the window rather than being judged against zero.
        struct A { double tag; bool want; const char* why; };
        for (const A a : { A { 48000.0,    true,  "the factory rate" },
                           A { 48000.5,    true,  "the high edge is INSIDE — the gate refuses past 0.5, not at it" },
                           A { 47999.5,    true,  "…and so is the low edge" },
                           A { 48000.5001, false, "a ten-thousandth of a hertz past the high edge" },
                           A { 47999.4999, false, "…and past the low edge" },
                           A { 44100.0,    false, "44.1 kHz, outright" },
                           A { 96000.0,    false, "96 kHz, outright" },
                           A { 0.0,        true,  "a model reporting no rate runs at the factory one" },
                           A { -1.0,       true,  "…which is what an untagged capture actually reports" } })
            test::ok (NamStage::acceptsModelRate (a.tag) == a.want,
                      "acceptsModelRate(" + std::to_string (a.tag) + ") = "
                      + (NamStage::acceptsModelRate (a.tag) ? "true" : "false") + " — " + a.why);
        // 🔴 AND THE EDGE IS PINNED AT THE EDGE, not a ten-thousandth of a hertz past it. A pre-merge
        // round found that a tolerance of 0.50001 — an excess ten times smaller than the cells
        // above resolve — passes every one of them while admitting a tag of 48000.500005. `nextafter`
        // is the only spelling that says "closed HERE"; 48000.5001 says "closed somewhere below here".
        test::ok (NamStage::acceptsModelRate (48000.5) && NamStage::acceptsModelRate (47999.5)
                      && ! NamStage::acceptsModelRate (std::nextafter (48000.5, 1.0e9))
                      && ! NamStage::acceptsModelRate (std::nextafter (47999.5, 0.0)),
                  "the window is closed at exactly [47999.5, 48000.5]: one ulp past either edge is "
                  "refused; the first representable tags beyond both edges are excluded");
        test::ok (NamStage::acceptsModelRate (std::numeric_limits<double>::quiet_NaN()),
                  "a NaN tag is not a rate at all: it takes the untagged door, exactly as rateMatch's "
                  "normalisation says it must — stated because the guard is `> 0.0`, which is FALSE "
                  "for a NaN, and that is easy to read the other way round");

        // 2. The same four edges through the real load path, BOTH entry points — the predicate being
        //    right is worth nothing if the loader asks it a different question.
        for (const A a : { A { 48000.5,    true,  "" }, A { 47999.5,    true,  "" },
                           A { 48000.5001, false, "" }, A { 47999.4999, false, "" } })
        {
            char tag[32];
            std::snprintf (tag, sizeof tag, "%.10g", a.tag);
            const auto json = gainModel (tag);
            NamStage one;   one.prepare (48000.0, 64);
            const bool fused = one.loadModelFromMemory (json.data(), json.size());
            NamStage two;   two.prepare (48000.0, 64);
            auto handle = NamStage::prepareModel (json.data(), json.size(), 48000.0, 64);
            const bool split = handle != nullptr && two.install (std::move (handle));
            test::ok (fused == a.want && split == a.want,
                      std::string ("a model tagged ") + tag + " loads=" + (fused ? "1" : "0")
                          + " through loadModelFromMemory and " + (split ? "1" : "0")
                          + " through prepareModel+install, want " + (a.want ? "1" : "0")
                          + " — and the two halves agree, which is what makes the gate's new home safe");
        }

        // 3. THE LADDER, which is the defect itself: the four modes of the probe, as a test. A single
        //    cell cannot see this — the first step is ACCEPTED in every one of them, and it is the
        //    SECOND that separates a fixed window from a walking one.
        {
            struct M { bool audio, prep; const char* name; };
            std::vector<float> l (64, 0.1f), r (64, 0.1f);
            float* io[2] { l.data(), r.data() };
            for (const M m : { M { false, false, "load only" },
                               M { true,  false, "load + audio" },
                               M { false, true,  "load + prepare" },
                               M { true,  true,  "load + audio + prepare (what a DAW does)" } })
            {
                NamStage s;
                s.prepare (48000.0, 64);
                double asked = 48000.0;
                int steps = 0;
                for (int i = 0; i < 64; ++i)          // 64 is far past the 1 a fixed window allows
                {
                    char tag[32];
                    std::snprintf (tag, sizeof tag, "%.10g", asked - 0.5);
                    const auto json = gainModel (tag);
                    if (! s.loadModelFromMemory (json.data(), json.size()))
                        break;
                    ++steps;
                    asked -= 0.5;
                    if (m.audio) test::run (s.process (io, 2, 64, false));
                    if (m.prep)  s.prepare (48000.0, 64);
                }
                test::ok (steps == 1 && NamStage::acceptsModelRate (s.modelSampleRate()),
                          std::string (m.name) + ": " + std::to_string (steps)
                              + " accepted half-hertz steps (want 1 — the first is inside the window, "
                                "the second is not), run rate now " + std::to_string (s.modelSampleRate()));
            }
        }

        // 4. …AND THE COST OF THE WALK WAS THE OPPOSITE OF WHAT IT LOOKED LIKE. A walked stage did not
        //    accept MORE, it accepted ELSEWHERE: at 47900 it refused an ordinary 48000 capture. So the
        //    fix cannot be checked only by what it now refuses — check what it keeps.
        {
            NamStage s;
            s.prepare (48000.0, 64);
            std::vector<float> l (64, 0.1f), r (64, 0.1f);
            float* io[2] { l.data(), r.data() };
            int accepted = 0;
            for (int i = 0; i < 200; ++i)             // 200 steps: enough to reach 47900 on the base
            {
                char tag[32];
                std::snprintf (tag, sizeof tag, "%.10g", 48000.0 - 0.5 * (i + 1));
                const auto json = gainModel (tag);
                if (s.loadModelFromMemory (json.data(), json.size()))
                    ++accepted;                       // no break: every step is ATTEMPTED, which is what
                test::run (s.process (io, 2, 64, false));   // the claim below says, and a loop that
                s.prepare (48000.0, 64);              // stopped at the first refusal would not say it
            }
            test::ok (accepted == 1, "precondition: of 200 attempted half-hertz steps exactly "
                                         + std::to_string (accepted) + " was accepted — the first, which "
                                         "is inside the window; the fixture is live and the walk is not");
            const auto plain = gainModel ("48000");
            test::ok (s.loadModelFromMemory (plain.data(), plain.size()) && s.modelSampleRate() == 48000.0,
                      "…and after all 200 an ordinary 48000 capture still loads — on the base commit "
                      "the stage had walked to 47900 by then and refused it");
        }

        // 5. THE INVARIANT, over sequences rather than cells: whatever order the public verbs come in,
        //    a held model is inside the window and the reported latency IS the policy's answer for the
        //    stage's own host rate. The second half is what F21 broke, and it is checked after EVERY
        //    step rather than at the end, because a wrong latency heals on the next re-prepare.
        {
            const double hosts[] = { 48000.0, 44100.0, 96000.0, 48000.4, 48000.6, 192000.0 };
            // 🔴 THE ALPHABET IS THE FIXTURE HERE, and a diverse-testing round showed the first draft's
            // was blind: every tag in it lay INSIDE the window, so on the base commit — where the
            // reference walked — the walk could only ever reach the same three rates and the window
            // check could not fail. "47999" and "48001" are what make it bite: on the base, a held
            // 47999.5 plus a prepare() moves the reference and 47999 is then accepted, which is exactly
            // the state this invariant says is unreachable.
            const char*  tags[]  = { "48000", "47999.5", "48000.5", "47999", "48001", "44100", "96000",
                                     nullptr, "0" };
            std::uint32_t rng = 0x38u;
            auto next = [&rng] { rng = rng * 1664525u + 1013904223u; return rng >> 16; };
            std::vector<float> l (64, 0.05f), r (64, 0.05f);
            float* io[2] { l.data(), r.data() };
            NamStage s;
            double host = 48000.0;
            s.prepare (host, 64);
            int held = 0, window = 0, latency = 0;
            for (int step = 0; step < 400; ++step)
            {
                // An alphabet alone does not exercise a sequence. With the original seed,
                // main passed both invariants in all 128 held states despite the added tags.
                // Force load(low edge), prepare, load(outside) before the random suffix.
                if (step == 0 || step == 2)
                {
                    const auto json = gainModel (step == 0 ? "47999.5" : "47999");
                    (void) s.loadModelFromMemory (json.data(), json.size());
                }
                else if (step == 1) s.prepare (host, 64);
                else switch (next() % 5u)
                {
                    case 0: host = hosts[next() % 6u]; s.prepare (host, 64); break;
                    case 1: { const char* t = tags[next() % 9u];
                              const auto json = gainModel (t);
                              (void) s.loadModelFromMemory (json.data(), json.size()); } break;
                    case 2: test::run (s.process (io, 2, 64, false)); break;
                    case 3: s.clearModel(); break;
                    default: (void) s.collectGarbage(); break;
                }
                if (s.hasModel())
                {
                    ++held;
                    if (NamStage::acceptsModelRate (s.modelSampleRate())) ++window;
                    if (s.latencySamples() == NamStage::rateMatch (host, s.modelSampleRate()).latencySamples)
                        ++latency;
                }
            }
            test::ok (held > 100, "precondition: the sequence checked more than 100 held-model states — "
                                  + std::to_string (held) + " of 400 steps, so the two checks below ran");
            test::ok (window == held, "…in every one of those states the held model is inside the window ("
                                      + std::to_string (window) + "/" + std::to_string (held) + ")");
            test::ok (latency == held, "…and in every one the reported latency IS rateMatch(hostSR, tag) ("
                                       + std::to_string (latency) + "/" + std::to_string (held)
                                       + ") — the property F21 broke");
        }

        // 6. F21 SPELLED OUT AS CELLS, because the invariant above would pass on a stage that simply
        //    re-prepared everything always. These are the cells where the two half-hertz tolerances
        //    straddle: a backend prepared for one host installed into a stage running another.
        {
            struct C { double preparedFor, stageHost; int want; const char* why; };
            for (const C c : { C { 48000.4, 48000.6, 64, "0.2 Hz apart: the skip used to keep a rate-match computed for the WRONG host — reported 0" },
                               C { 48000.3, 48000.7, 64, "0.4 Hz apart: the same door, one step wider" },
                               C { 48000.6, 48000.4,  0, "the MIRROR — this one used to over-report, charging 64 where the policy charges none" },
                               C { 48000.0, 48000.6, 64, "0.6 Hz apart: outside the old tolerance, so this cell was RIGHT before and must stay right" },
                               C { 48000.0, 48000.5,  0, "…and inside it, where the answers happen to agree" },
                               // 🔴 THE CELL THAT SEPARATES "EXACT" FROM "A SMALLER TOLERANCE". A
                               // pre-merge round pointed out that every row above also passes with a
                               // 1e-6 tolerance in place of equality — the rows are all further apart
                               // than that. These two straddle the policy's own boundary by a
                               // quarter-millionth of a hertz on either side, so a 1e-6 test misses them.
                               // The adjacent-double rows below resolve the boundary further.
                               C { 48000.49999975, 48000.50000025, 64,
                                   "half a millionth of a hertz apart, straddling the gate: exact "
                                   "reconfigures; a 1e-6 tolerance does not" },
                               C { 48000.5, std::nextafter (48000.5, 1.0e9), 64,
                                   "adjacent representable hosts straddle the gate" },
                               C { std::nextafter (48000.5, 1.0e9), 48000.5, 0,
                                   "the adjacent-host mirror must remove the converter" } })
            {
                NamStage s;
                s.prepare (c.stageHost, 64);
                const auto json = gainModel ("48000");
                auto handle = NamStage::prepareModel (json.data(), json.size(), c.preparedFor, 64);
                const bool ok = handle != nullptr && s.install (std::move (handle));
                test::ok (ok && s.latencySamples() == c.want,
                          "prepared for host " + std::to_string (c.preparedFor) + ", installed into "
                              + std::to_string (c.stageHost) + ": reports "
                              + std::to_string (s.latencySamples()) + ", want " + std::to_string (c.want)
                              + " — " + c.why);
            }
        }

        // 7. maxLatencySamples: the consequence the contract owes a consumer sizing a ring. Values by
        //    hand from the two documented facts — kHalf per leg, the return leg converted at h/m — at
        //    the window's LOW edge 47999.5, rounded to nearest, and NOT by asking the code.
        struct L { double h; int want; const char* why; };
        for (const L x : { L {  48000.0,   64, "32 + 32*48000/47999.5 = 64.000333" },
                           L {  44100.0,   61, "32 + 29.400306 = 61.400306" },
                           L {  96000.0,   96, "32 + 64.000667 = 96.000667" },
                           L { 192000.0,  160, "32 + 128.001333" },
                           L { 384000.0,  288, "32 + 256.002667" },
                           L { 3.0e6,    2032, "32 + 2000.020834 — the house ceiling for a host rate" },
                           L { 2999249.0, 2032, "🔴 THE CELL THAT SEPARATES THE TWO DERIVATIONS: asking at "
                                                "the NOMINAL 48000 gives 2031 here, and a ring sized from that "
                                                "clamps an accepted model by one sample, in silence" } })
            test::ok (NamStage::maxLatencySamples (x.h) == x.want,
                      "maxLatencySamples(" + std::to_string (x.h) + ") = "
                          + std::to_string (NamStage::maxLatencySamples (x.h)) + ", want "
                          + std::to_string (x.want) + " — " + x.why);
        test::ok (NamStage::maxLatencySamples (2999249.0)
                      > nam::NamStage::rateMatch (2999249.0, nam::NamStage::kModelSampleRate).latencySamples,
                  "…and that cell is a real difference, not a restatement: the nominal rate answers "
                  + std::to_string (nam::NamStage::rateMatch (2999249.0, nam::NamStage::kModelSampleRate).latencySamples)
                  + " where the window's low edge answers "
                          + std::to_string (NamStage::maxLatencySamples (2999249.0)));
        // Exact half-integers at the low accepted tag expose even a tiny upward nudge
        // of the tag used to size the ring. A sweep of ordinary rates cannot see it.
        for (const double h : { 48749.4921875, 96748.9921875, 191248.0078125, 2999218.7578125 })
            test::ok (NamStage::maxLatencySamples (h)
                          >= NamStage::rateMatch (h, 47999.5).latencySamples,
                      "the bound covers a low-edge model at a half-integer rounding boundary");
        // …and the bound really is a BOUND: sweep the window against it at every shipped host rate.
        {
            int checked = 0, covered = 0, hosts = 0, tightHosts = 0;
            for (const double h : { 8000.0, 11025.0, 16000.0, 22050.0, 32000.0, 44100.0, 47999.0, 48000.0,
                                    48000.4, 48000.6, 48001.0, 88200.0, 96000.0, 176400.0, 192000.0,
                                    352800.0, 384000.0, 2999249.0, 3.0e6 })
            {
                const int bound = NamStage::maxLatencySamples (h);
                bool attainedHere = false;
                for (int k = 0; k <= 20; ++k)             // the whole window, in twentieths
                {
                    const double m = 47999.5 + 0.05 * k;
                    const int rep = NamStage::rateMatch (h, m).latencySamples;
                    ++checked;
                    if (rep <= bound) ++covered;
                    if (rep == bound) attainedHere = true;
                }
                ++hosts;
                if (attainedHere) ++tightHosts;
            }
            test::ok (covered == checked, "the bound covers every model rate in the window at every "
                                          "shipped host: " + std::to_string (covered) + "/"
                                          + std::to_string (checked));
            // 🔴 TIGHT PER HOST, NOT IN AGGREGATE. A pre-merge round pointed out that "attained
            // somewhere" is satisfied by inflating the bound at every host but one, which is precisely
            // the failure a tightness check exists to catch. The exception is the single host where the
            // whole window is inside the resampling gate — hostSR exactly kModelSampleRate — and it is
            // named rather than tolerated.
            test::ok (tightHosts == hosts - 1,
                      "…and the bound is attained AT EVERY HOST but one: " + std::to_string (tightHosts)
                          + " of " + std::to_string (hosts) + ", the exception being hostSR exactly "
                            "48000, where no accepted model resamples and the true maximum is 0");
        }
    }

    test::group ("🔴 THE MODEL SCRATCH IS SIZED FOR THE PATH THAT RUNS, not for a ratio (P38)");
    {
        using nam::NamStage;
        // The scratch used to be `ceil(maxBlock * modelRunSR / max(8000, hostSR)) + 16`, and that one
        // number was asked to serve two paths that need different ones. Both halves were measured
        // wrong, and both were silent.

        // 🔴 HALF ONE — PAST THE END OF THE HEAP. With no resampler the model is clocked by the HOST,
        // so a whole chunk of maxBlock host samples is copied into the scratch and the ratio never
        // enters. A tag half a hertz BELOW the host is inside the acceptance window (the acceptance
        // gate and the resampling gate are the same half hertz), so the ratio is just under one and
        // the scratch came out just under maxBlock. Short by `maxBlock*(h-m)/h - 16` samples, so the
        // first failing block is 34*hostSR = 34 SECONDS of audio in one call — which law 11(a)
        // explicitly invites an offline caller to pass. Measured on the base commit under ASan at
        // maxBlock 2000000: heap-buffer-overflow, a WRITE of 8000000 bytes into a 7999984-byte region.
        // This cell is at the first undersized block rather than comfortably past it, so it also pins
        // WHERE the threshold is; the CI sanitizer job builds nam, so the mutant dies there by name.
        {
            const int blk = 1632000;                      // 34 * 48000, by the derivation above
            NamStage s;
            s.prepare (48000.0, blk);
            const auto json = gainModel ("47999.5");      // the window's low edge: accepted, not resampled
            test::ok (s.loadModelFromMemory (json.data(), json.size()),
                      "precondition: the low-edge tag loads at a 34-second block");
            test::ok (s.latencySamples() == 0,
                      "precondition: and it runs with NO resampler — the direct path is the one under "
                      "test, not the converted one");
            std::vector<float> in ((std::size_t) blk);
            for (int i = 0; i < blk; ++i)                 // a signal whose every sample is distinct, so
                in[(std::size_t) i] = (float) ((i % 2039) - 1019) * (1.0f / 2048.0f);   // a stale slot shows
            std::vector<float> io = in;
            float* p[1] { io.data() };
            test::run (s.process (p, 1, blk, false));
            std::size_t wrong = 0, first = 0;
            for (std::size_t i = 0; i < in.size(); ++i)
                if (io[i] != in[i]) { if (wrong == 0) first = i; ++wrong; }
            test::ok (wrong == 0, "a unity capture returns a 34-second block sample for sample — "
                                  + std::to_string (wrong) + " samples differ"
                                  + (wrong ? ", first at " + std::to_string (first) : std::string()));
        }

        // 🔴 HALF THREE — AND WHAT REPLACED THE ASSUMED RATE IS A REFUSAL, NOT ANOTHER SUBSTITUTE.
        // A host rate that cannot size a conversion has no honest frame count, and there are two ways
        // to answer that: invent one — which is what `max(8000, hostSR)` did, and it cost the samples
        // measured below — or say so. This says so, through the mechanism the class already had for a
        // preparation it cannot honour: the backend stays unprepared, so the LOAD FAILS visibly at the
        // call, with the stage untouched and still passing audio through.
        //
        // BEHAVIOUR CHANGE, stated rather than left to be discovered: at these rates a load used to
        // SUCCEED and then produce garbage or silence, and at the far end it converted an out-of-range
        // double to an int, which is undefined. Nothing shipped reaches them — rigplayer maps every
        // host outside (0, 3e6] to the factory rate before the stage sees it — and with a 512-sample
        // block the arithmetic lower limit is about 0.023 Hz; that is not a promise that allocation
        // succeeds or that StreamResampler supports the ratio (its limit is 1e6:1).
        {
            struct R { double fs; bool want; const char* why; };
            for (const R r : { R { 48000.0, true,  "the factory host, for contrast" },
                               R {     0.0, false, "zero: there is no ratio, so there is no frame count" },
                               R {-48000.0, false, "negative: the same, and it used to be floored to 8 kHz" },
                               // 🔴 THE NEGATIVE HOST THAT THE ARITHMETIC ALONE LETS THROUGH, and the
                               // reason the guard is written as a rate test rather than a frame-count
                               // test. At -2000000 the converted term is -12 and the slack of 16 makes
                               // it a positive, perfectly representable 4 — so a rule that inspects only
                               // the result accepts it. -48000 gives -496 and is refused, which is how a
                               // wrong rule looks right on the cell you happened to try; a mutation
                               // round deleted the guard and this suite stayed green without this cell.
                               R {-2000000.0, false, "the slack outweighs the negative converted term: 4 frames for a rate that converts nothing" },
                               R {   -1.0e7, false, "…and it is not one freak value either" },
                               R {   1e-30, false, "a host so slow one block converts to more frames than the arithmetic holds" },
                               R {    0.01, false, "…and the boundary is a RATE, not a special value: 0.01 Hz is refused" },
                               R {    1.0,  true,  "…while 1 Hz is absurd but sizeable, and IS honoured" } })
            {
                nam::NamStage s;
                s.prepare (r.fs, 512);
                const auto json = gainModel ("48000");
                const bool ok = s.loadModelFromMemory (json.data(), json.size());
                test::ok (ok == r.want, "a load at host " + std::to_string (r.fs) + " returns "
                                            + (ok ? "true" : "false") + ", want "
                                            + (r.want ? "true" : "false") + " — " + r.why);
                std::vector<float> b (8, 0.25f);
                float* io[1] { b.data() };
                const bool ran = s.process (io, 1, 8, false);
                // A REFUSED load leaves the stage empty, and an empty stage is a passthrough — that is
                // what makes the refusal safe rather than merely honest. An ACCEPTED one is only asked
                // to run: at a 1 Hz host it runs THROUGH a rate-matcher at a ratio of 48000, and
                // demanding 0.25 back from that would be asserting the resampler's arithmetic in a
                // group about refusals. (The first draft did demand it, and this row is where it fell.)
                test::ok (ran && (ok || b[0] == 0.25f),
                          std::string ("…the stage answers process() ") + (ok ? "with its model" : "")
                              + (ok ? "" : "and a refused load leaves the audio alone"));
            }
            // 🔴 AND THE BLOCK RIDES THE SAME ARITHMETIC, so it is held to the same ceiling — which is
            // not decoration either: `maxBlock * 2 + 16` is handed to the down-converter, and past
            // (INT_MAX-16)/2 that expression overflowed, quietly, on the base commit. The refusal costs
            // nothing to test because it happens BEFORE a byte is allocated.
            {
                nam::NamStage big;
                big.prepare (48000.0, 1073741816);        // (INT_MAX - 16)/2 + 1
                const auto json = gainModel ("48000");
                test::ok (! big.loadModelFromMemory (json.data(), json.size()),
                          "a block one past (INT_MAX-16)/2 is refused rather than doubled into an "
                          "overflow one line later");
                nam::NamStage ok2;
                ok2.prepare (48000.0, 1024);
                test::ok (ok2.loadModelFromMemory (json.data(), json.size()),
                          "…precondition: an ordinary block at the same host still loads, so the row "
                          "above is the ceiling and not a broken fixture");
            }

            // 🔴 AND THE SAME REFUSAL REACHES install(), WHICH IS WHERE IT WAS UNGATED. A diverse-testing
            // round deleted install()'s verdict check and the whole suite stayed green: the check had
            // been guarded only by an invariant two functions away (prepareModel returns null for a
            // backend that failed), and nothing asserted the case where the RE-preparation is the one
            // that fails. A handle prepared for a host this stage can honour, installed into a stage
            // whose host it cannot, must be refused with the stage left empty — not installed as a
            // model that reports itself present and passes audio through.
            for (const double bad : { 0.0, -48000.0, 1e-30 })
            {
                nam::NamStage s;
                s.prepare (bad, 64);
                const auto json = gainModel ("48000");
                auto handle = nam::NamStage::prepareModel (json.data(), json.size(), 48000.0, 64);
                test::ok (handle != nullptr,
                          "precondition: the handle itself prepares fine at 48000 — the refusal below "
                          "belongs to the STAGE's host rate, not to the model");
                test::ok (! s.install (std::move (handle)) && ! s.hasModel(),
                          "a handle prepared at 48000 and installed into a stage at host "
                              + std::to_string (bad) + " is refused, and the stage stays empty");
                std::vector<float> b (8, 0.25f);
                float* io[1] { b.data() };
                test::ok (s.process (io, 1, 8, false) && b[0] == 0.25f,
                          "…and an empty stage passes its audio through untouched");
            }

            // 🔴 THE BLOCK'S OWN CEILING, ON THE BRANCH WHERE IT IS THE ONLY THING THAT BINDS. At a host
            // far above the model rate the converted frame count is tiny — 68 at a 1e12 Hz host — so
            // `frames <= kMaxFrames` passes and only the clause about maxBlock stops `maxBlock * 2 + 16`
            // from overflowing. The cell above (a 48 kHz host) cannot reach this clause at all, because
            // there `frames` IS maxBlock. Named by a diverse-testing round, which deleted the clause and
            // watched the suite stay green.
            {
                nam::NamStage far;
                far.prepare (1e12, 1073741816);
                const auto json = gainModel ("48000");
                test::ok (! far.loadModelFromMemory (json.data(), json.size()),
                          "a 1.07e9-sample block at a 1e12 Hz host is refused — the converted count is "
                          "68 and passes, so this is the maxBlock clause and nothing else");
                nam::NamStage near;
                near.prepare (1e12, 1024);
                test::ok (near.loadModelFromMemory (json.data(), json.size()),
                          "…precondition: the same absurd host with an ordinary block still loads, so "
                          "the row above is the block ceiling and not the host");
            }

            // The precondition without which the rows above would be worthless: an INFINITE host is a
            // different animal and is deliberately NOT refused — its ratio is zero, which is a number,
            // and what it then REPORTS belongs to rateMatch's documented platform-specific regimes.
            // Stated so nobody reads this group as "absurd rates are refused".
            nam::NamStage e;
            e.prepare (std::numeric_limits<double>::infinity(), 512);
            const auto json = gainModel ("48000");
            test::ok (e.loadModelFromMemory (json.data(), json.size()),
                      "an infinite host still loads — the refusal is about a frame count that cannot "
                      "be represented, not about a rate looking unreasonable");
        }

        // 🔴 HALF TWO — AN ASSUMED HOST RATE, AND THE SAMPLES IT LOSES. `max(8000, hostSR)` put a
        // substitute rate into the sizing, so below 8 kHz a block converted to more model frames than
        // fit and produceAvailable() dropped the surplus without a word. Measured on the base commit
        // with a unity capture and a 100 Hz tone: -1.22 dB at a 6 kHz host, -3.00 at 4 kHz, -6.05 at
        // 2 kHz, -9.09 at 1 kHz — and exactly 0.00 at 8 kHz, the floor's own edge, which is what
        // identified the floor rather than the resampler as the cause. Every one of those rates is
        // accepted by rigplayer::RigPlayer::usableSampleRate.
        {
            struct H { double fs; double floorEdge; };
            double worst = 0.0;
            double atFloor = 0.0;
            for (const double fs : { 1000.0, 2000.0, 4000.0, 6000.0, 8000.0 })
            {
                NamStage s;
                s.prepare (fs, 512);
                const auto json = gainModel ("48000");
                test::ok (s.loadModelFromMemory (json.data(), json.size()),
                          "precondition: a factory capture loads at a " + std::to_string ((int) fs)
                              + " Hz host");
                double si = 0.0, so = 0.0, corr = 0.0;
                int counted = 0;
                std::vector<float> l (512);
                std::vector<float> whole, back;          // the full input and the full output, so the
                whole.reserve (512 * 200); back.reserve (512 * 200);   // phase check can ALIGN them
                float* p[1] { l.data() };
                for (int off = 0; off < 512 * 200; off += 512)
                {
                    for (int i = 0; i < 512; ++i)
                        l[(std::size_t) i] = (float) std::sin (2.0 * kPi * 100.0 * (off + i) / fs);
                    whole.insert (whole.end(), l.begin(), l.end());
                    if (off >= 512 * 100)
                        for (int i = 0; i < 512; ++i) si += (double) l[(std::size_t) i] * l[(std::size_t) i];
                    test::run (s.process (p, 1, 512, false));
                    back.insert (back.end(), l.begin(), l.end());
                    if (off >= 512 * 100)
                    {
                        for (int i = 0; i < 512; ++i) so += (double) l[(std::size_t) i] * l[(std::size_t) i];
                        counted += 512;
                    }
                }
                const double lossDb = 10.0 * std::log10 ((so / counted) / (si / counted));
                // 🔴 ENERGY IS BLIND TO SIGN, and a pre-merge round said so: squaring the samples makes a
                // polarity inversion identical to the truth. A correlation closes it — ALIGNED, because
                // the stage is resampling here and reports 33 samples of it at a 1 kHz host, and an
                // unaligned correlation reads NEGATIVE on a 100 Hz tone whose period is ten samples.
                // (The first draft of this check did exactly that and failed on correct code at four of
                // the five rates, which is what a phase-blind oracle looks like from the other side.)
                const int lat = s.latencySamples();
                for (std::size_t i = (std::size_t) (512 * 100 + lat); i < back.size(); ++i)
                    corr += (double) back[i] * (double) whole[i - (std::size_t) lat];
                test::ok (corr > 0.0, "…and the output is in PHASE with the input at "
                                          + std::to_string ((int) fs) + " Hz once its own "
                                          + std::to_string (lat) + " samples of rate-match delay are "
                                          "taken out (correlation " + std::to_string (corr)
                                          + ") — the energy figure above cannot tell an inversion "
                                            "from a match");
                if (fs == 8000.0) atFloor = lossDb;
                else if (std::fabs (lossDb) > std::fabs (worst)) worst = lossDb;   // largest EXCURSION,
                                                                                   // whatever its sign
                test::ok (lossDb > -0.5,
                          "a 100 Hz tone through a unity capture at a " + std::to_string ((int) fs)
                              + " Hz host loses " + std::to_string (lossDb)
                              + " dB, want better than -0.5 (the base commit lost up to -9.09)");
            }
            // 🔴 THE PRECONDITION THIS GROUP NEEDS IS ABOUT THE SWEEP'S REACH, NOT ABOUT ITS SIGN. A first
            // draft asserted that the worst rate below 8 kHz still reads BELOW zero, and a
            // diverse-testing round measured why that is not a property of anything under test: what is
            // left at 1 kHz is the resampler kernel's own passband droop (-0.069 dB there, -0.045 at
            // 2 kHz, -0.003 at 4 kHz and +0.0005 at 6 kHz — it changes SIGN across the sweep), so a
            // flatter kernel would fail the precondition on perfectly correct code. What must be true
            // is that the sweep straddles the rate the removed floor stood at; the residual belongs to
            // the kernel and is asserted only as "small", by the -0.5 dB rows above.
            test::ok (atFloor > -0.01 && std::fabs (worst) < 0.5,
                      "precondition: the sweep straddles the old 8 kHz floor — 8 kHz itself reads "
                      + std::to_string (atFloor) + " dB and the extreme below it reads "
                      + std::to_string (worst) + ", both inside the kernel's own residual, so a fixture "
                      "that only ran at or above the floor could not have seen the -9.09 dB at all");
        }
    }

    test::group ("\U0001f534 THE RESET BLOCK IS PART OF THE CONTRACT — read out through a STATEFUL capture");
    {
        using nam::NamStage;
        // WHY THIS GROUP EXISTS, and it is not a nicety. `maxModelFrames` is not only a buffer size: it
        // is the block handed to NAM's `Reset`, and NAM prewarms in WHOLE blocks of it. So it decides
        // how many samples of silence a capture is warmed with, and for anything with recurrent state
        // that decides the first real output. Every other model in this file is memoryless, so this
        // was invisible: a mutation stand ran the "+16 slack" and the "max() instead of the branch"
        // mutants against the whole suite and called both EQUIVALENT — and both of those verdicts were
        // wrong, which only showed once lstmModel() existed. A pre-merge round is what asked for it.

        // Read the first output that carries model state: at a resampling host the leading
        // latencySamples() outputs are the converter's own zeros.
        auto firstOut = [] (double host, int blk, const char* tag)
        {
            NamStage s;
            s.prepare (host, blk);
            const auto json = lstmModel (tag);
            if (! s.loadModelFromMemory (json.data(), json.size()))
                return std::numeric_limits<double>::quiet_NaN();
            const int lat = s.latencySamples();
            std::vector<float> whole;
            for (int k = 0; k < 4; ++k)
            {
                std::vector<float> b ((std::size_t) blk, 0.0f);
                float* io[1] { b.data() };
                test::run (s.process (io, 1, blk, false));
                whole.insert (whole.end(), b.begin(), b.end());
            }
            return (double) whole[(std::size_t) lat];
        };

        // PRECONDITION — the instrument reads the RESET BLOCK, not just "something". Three block sizes
        // at one host must give three DIFFERENT answers; if they did not, every assertion below would
        // hold vacuously and the group would be decoration.
        const double a64 = firstOut (48000.0, 64, "48000");
        const double a512 = firstOut (48000.0, 512, "48000");
        const double a4096 = firstOut (48000.0, 4096, "48000");
        test::ok (std::isfinite (a64) && std::isfinite (a512) && std::isfinite (a4096)
                      && a64 != a512 && a512 != a4096 && a64 != a4096,
                  "precondition: the fixture FEELS the Reset block — blocks 64/512/4096 read "
                  + std::to_string (a64) + " / " + std::to_string (a512) + " / "
                  + std::to_string (a4096) + ", three different numbers");

        // THE PINS. Values by construction: N = ceil(24000 / R) * R prewarm samples with R the Reset
        // block, and the output is 0.5*tanh(sigmoid(10)^N). R = maxBlock + 16 on the direct path, so
        // 80 / 528 / 4112 here. Tolerance 1e-6 — see lstmModel()'s note on why that survives three libms.
        test::approx (a64,   0.1620296240, 1e-6, "48000 Hz, block 64: Reset block 80, prewarm 24000");
        test::approx (a512,  0.1600718498, 1e-6, "48000 Hz, block 512: Reset block 528, prewarm 24288");
        test::approx (a4096, 0.1574926972, 1e-6, "48000 Hz, block 4096: Reset block 4112, prewarm 24672");

        // 🔴 THE SURVIVOR THIS GROUP WAS WRITTEN FOR. On the direct path the model is clocked by the
        // HOST and the ratio never enters, so the Reset block must NOT depend on the model's tag. The
        // mutant that sizes with `max(maxBlock, ceil(maxBlock*m/h) + 16)` breaks exactly this: at a tag
        // above the nominal rate it reads one frame more (529 against 528 at block 512), which moves
        // the prewarm and moves this number. That mutant survived two full mutation rounds before this
        // assertion existed.
        for (const int blk : { 64, 512, 4096 })
        {
            const double nominal = firstOut (48000.0, blk, "48000");
            const double high    = firstOut (48000.0, blk, "48000.5");
            const double low     = firstOut (48000.0, blk, "47999.5");
            test::ok (nominal == high && nominal == low,
                      "at a 48 kHz host, block " + std::to_string (blk)
                          + ", the Reset block does NOT depend on the tag: 48000 / 48000.5 / 47999.5 "
                            "all read " + std::to_string (nominal)
                          + " (mutant: " + std::to_string (high) + " for the high tag)");
        }

        // …AND THE CONVERTED PATH KEEPS ITS OWN SLACK, which is the other verdict this fixture
        // overturned. At a 96 kHz host the converted count is ceil(512*48000/96000) + 16 = 272; drop
        // the slack and it is 256, a different prewarm and a different number. Measured on the mutant:
        // this cell moves by ~5e-4, five hundred times the tolerance.
        test::approx (firstOut (96000.0, 512, "48000"), 0.1602614224, 1e-6,
                      "96 kHz host, block 512: the converted Reset block is 272, not 256");
        test::approx (firstOut (44100.0, 512, "48000"), 0.1610591710, 1e-6,
                      "44.1 kHz host, block 512: converted Reset block 574");
    }

    test::group ("P38 stateful prewarm and live preparation refusal");
    {
        using nam::NamStage;
        const auto json = lstmModel();
        for (const int block : { 64, 256, 512, 1024, 4096 })
        {
            NamStage s; s.prepare (48000.0, block);
            test::ok (s.loadModelFromMemory (json.data(), json.size()), "stateful fixture loads");
            float x = 0.0f; float* io[] { &x };
            test::run (s.process (io, 1, 1, false));
            const int resetBlock = block + 16;
            const int warmed = ((24000 + resetBlock - 1) / resetBlock) * resetBlock;
            test::approx (x, stateAfterSilence (warmed + 1), 2e-6,
                          "one direct preparation advances exactly the scheduled whole blocks");
        }
        {
            NamStage s; s.prepare (44100.0, 512);
            test::ok (s.loadModelFromMemory (json.data(), json.size()), "converted stateful fixture loads");
            std::vector<float> x (512, 0.0f); float* io[] { x.data() };
            test::run (s.process (io, 1, 512, false));
            // Signal pin includes the 574-frame Reset block and the converter's transient.
            test::approx (x.back(), 0.1577584296, 2e-6, "converted preparation preserves its state trajectory");
        }
        {
            // Current behavior: NAM's LSTM Reset prewarms the existing cell again. An exact
            // host mismatch now reaches this even when both hosts use the direct path.
            NamStage s; s.prepare (48000.2, 512);
            auto handle = NamStage::prepareModel (json.data(), json.size(), 48000.1, 512);
            test::ok (s.install (std::move (handle)), "same-path host mismatch installs");
            float x = 0.0f; float* io[] { &x };
            test::run (s.process (io, 1, 1, false));
            test::approx (x, stateAfterSilence (2 * 24288 + 1), 2e-6,
                          "host mismatch adds a second prewarm even without a latency change");
        }
        {
            NamStage s; s.prepare (48000.0, 512);
            auto handle = NamStage::prepareModel (json.data(), json.size(), 48000.0, 64);
            test::ok (s.install (std::move (handle)), "block mismatch installs");
            float x = 0.0f; float* io[] { &x };
            test::run (s.process (io, 1, 1, false));
            test::approx (x, stateAfterSilence (24000 + 24288 + 1), 2e-6,
                          "block mismatch advances the state through the second prewarm");
        }
        {
            NamStage s; s.prepare (47999.75, 512);
            test::ok (s.loadModelFromMemory (json.data(), json.size()), "integer tag at fractional host loads");
            float x = 0.0f; float* io[] { &x };
            test::run (s.process (io, 1, 1, false));
            test::approx (x, stateAfterSilence (24288 + 1), 2e-6,
                          "the direct Reset block follows the host block even when the tag is higher");
        }
        {
            NamStage s; s.prepare (48000.0, 64);
            const auto unity = gainModel();
            test::ok (s.loadModelFromMemory (unity.data(), unity.size()), "live refusal starts with a model");
            s.prepare (0.0, 64);
            float x = 0.25f; float* io[] { &x };
            test::ok (! s.process (io, 1, 1, false) && x == 0.25f
                          && s.hasModel() && s.latencySamples() == 0,
                      "failed live reprepare refuses processing, retains the model, and preserves audio");
            s.prepare (48000.0, 64);
            test::ok (s.process (io, 1, 1, false) && x == 0.25f, "valid preparation recovers the live model");
        }
        {
            // Compare identical architectures and equal-width tags. No wall-clock threshold:
            // rejection before backend construction/Reset must avoid their allocations.
            const auto rejected = lstmModel ("96000");
            auto before = g_allocs.load();
            auto acceptedHandle = NamStage::prepareModel (json.data(), json.size(), 48000.0, 512);
            const auto acceptedAllocs = g_allocs.load() - before;
            before = g_allocs.load();
            auto rejectedHandle = NamStage::prepareModel (rejected.data(), rejected.size(), 48000.0, 512);
            const auto rejectedAllocs = g_allocs.load() - before;
            test::ok (acceptedHandle != nullptr && rejectedHandle == nullptr, "allocation comparison has both outcomes");
            test::ok (rejectedAllocs < acceptedAllocs, "rate rejection avoids backend/prewarm allocations");
        }
    }

    test::group ("🔴 RATE-CHANGE NULLS — the resampler must not survive its own reconfiguration");
    {
        // THE MAIN DEFENCE OF P34. NamStage engages the rate-matcher only past
        // |hostSR - modelRunSR| > 0.5, so at a 48 kHz host with a 48 kHz capture there is no resampler
        // in the path at all — processChannel takes its `if (! resampling)` branch and StreamResampler
        // is never called. A kernel swap therefore has to leave this output BIT-IDENTICAL, and any
        // divergence means the change reached somewhere it was not aimed.
        //
        // The checksum below was computed on the BASE commit (main = 52582a8, the Catmull-Rom kernel)
        // with this same fixture and compared against the same computation after the swap: identical.
        // It is a FNV-1a over the raw bit patterns of every output sample, so it cannot be satisfied by
        // anything short of bit equality — a one-ulp difference in one sample changes it completely.
        //
        // The signal deliberately covers what a rate-match would disturb if it were wrongly engaged:
        // a sweep through the top octave (where the two kernels differ by 5.5 dB), a DC step, silence,
        // and a full-scale impulse. It is 4096 samples through the FIR fixture, which is the model with
        // memory — a gain model would hide a one-sample misalignment.
        auto checksum = [] (double hostSR)
        {
            nam::NamStage stage;
            stage.prepare (hostSR, 64);
            const auto json = firModel();
            if (! load (stage, json)) return (std::uint64_t) 0;
            stage.prepare (hostSR, 64);

            std::vector<float> in (4096, 0.0f);
            for (int i = 0; i < 2048; ++i)                       // sweep 8 kHz -> 22 kHz
            {
                const double t = (double) i / 2048.0;
                const double f = 8000.0 + 14000.0 * t;
                in[(std::size_t) i] = 0.5f * (float) std::sin (2.0 * 3.14159265358979323846 * f * i / hostSR);
            }
            for (int i = 2048; i < 2560; ++i) in[(std::size_t) i] = 0.75f;      // DC step
            for (int i = 2560; i < 3072; ++i) in[(std::size_t) i] = 0.0f;       // silence
            in[3072] = 1.0f;                                                    // full-scale impulse

            const auto out = runMono (stage, in, false);
            std::uint64_t h = 1469598103934665603ull;                           // FNV-1a offset basis
            for (float v : out)
            {
                std::uint32_t bits = 0;
                std::memcpy (&bits, &v, sizeof (bits));
                for (int byte = 0; byte < 4; ++byte)
                {
                    h ^= (std::uint64_t) ((bits >> (8 * byte)) & 0xffu);
                    h *= 1099511628211ull;
                }
            }
            return h;
        };

        // 🔴 WHY THE LITERAL HASH IS NOT ASSERTED HERE. The cross-tree comparison IS the acceptance and
        // it was done: base main = 52582a8 (Catmull-Rom) and this branch both render 0x8f19a4552add60ff
        // at 48 kHz on this machine — bit for bit — while 44.1 kHz moves from 0x8a79981a41b078b2 to
        // 0x9ba23a94952a4fab, which is the divergence the fix is FOR. But that hash is a NAM render
        // through Eigen: it is not expected to survive a change of toolchain or FMA contraction, so
        // pinning the constant would buy a red CI row on another platform and prove nothing extra.
        // What IS pinned below is the same defence in a platform-independent form.
        const std::uint64_t got = checksum (48000.0);
        std::printf ("      48 kHz FNV-1a over the whole render: 0x%016llx  (base main: 0x8f19a4552add60ff)\n",
                     (unsigned long long) got);

        // Precondition, so the lines below cannot pass by measuring nothing: at 44.1 kHz, where the
        // resampler IS engaged, the same render must hash differently. Without this an all-zero render
        // would satisfy any null test forever.
        const std::uint64_t off = checksum (44100.0);
        std::printf ("      44.1 kHz, where the resampler IS in the path: 0x%016llx\n", (unsigned long long) off);
        test::ok (off != got,
                  "precondition: the instrument is not blind — at 44.1 kHz, where the rate-match runs, "
                  "the same render hashes differently");

        // THE PORTABLE GATE. The failure this whole item defends against is the kernel leaking into a
        // path it does not belong on, and the only route it could take is STATE: a stage that has been
        // configured for resampling and then re-prepared at the model's own rate. configureRates()
        // resets both StreamResamplers on every prepare, so a stage that has been through 44.1 kHz and
        // back to 48 must render exactly what a virgin one does — bit for bit, on any toolchain.
        {
            nam::NamStage viaResampling;
            viaResampling.prepare (44100.0, 64);
            const auto json = firModel();
            test::ok (load (viaResampling, json), "model loads on the stage that starts at 44.1 kHz");
            viaResampling.prepare (44100.0, 64);

            viaResampling.prepare (48000.0, 64);                                   // …then come back to 48 kHz

            // 🔴 NOTE ON WHAT IS DELIBERATELY *NOT* DONE HERE, because getting it wrong cost a round.
            // The obvious stronger version — run audio at 44.1 kHz before re-preparing — measures the
            // wrong thing: the two stages then differ from sample 0, and they differ IDENTICALLY on
            // main with the Catmull-Rom kernel (checked, cross-tree). That divergence belongs to the
            // NAM model's own re-prewarm across a prepare(), not to the resampler, and asserting it
            // here would pin someone else's behaviour onto this item. Recorded as a finding instead.
            // What this gate does assert is the part that IS about the resampler: a stage whose
            // StreamResamplers were configured for a 44.1 kHz ratio and then re-configured for 48 kHz
            // must render exactly what one that never saw another rate does.

            std::vector<float> in (1024, 0.0f);
            for (int i = 0; i < 512; ++i)
                in[(std::size_t) i] = 0.5f * (float) std::sin (2.0 * 3.14159265358979323846 * 19000.0 * i / 48000.0);
            in[600] = 1.0f;

            nam::NamStage virginStage;
            virginStage.prepare (48000.0, 64);
            test::ok (load (virginStage, json), "…and on a stage that has never seen another rate");
            virginStage.prepare (48000.0, 64);

            const auto a = runMono (viaResampling, in, false);
            const auto c2 = runMono (virginStage, in, false);
            bool identical = (a.size() == c2.size());
            for (std::size_t i = 0; identical && i < a.size(); ++i)
                identical = (std::memcmp (&a[i], &c2[i], sizeof (float)) == 0);
            test::ok (identical && a.size() > 900,
                      "a stage that has RUN the rate-matcher and been re-prepared at 48 kHz renders "
                      "bit-identically to one that never did — the kernel cannot reach the path where "
                      "|hostSR - modelRunSR| <= 0.5, by state or otherwise");
            test::ok (viaResampling.latencySamples() == 0,
                      "…and it reports no latency there either, which is the same claim in the number "
                      "the host acts on");
        }

        // 🔴 AND THE FORM THAT ACTUALLY CATCHES A STALE TABLE. A diverse-testing round mutated
        // StreamResampler::reset() to early-return when its table was already populated — a resampler
        // that keeps designing for the ratio it saw first — and that mutant passed 345 of 345 checks
        // across this repository, INCLUDING the gate above. The gate above cannot see it: at 48 kHz the
        // resampler is never invoked, so `viaResampling ≡ virgin` by construction whatever the table
        // holds. What catches it is a rate change that lands somewhere the resampler DOES run.
        //
        // The fixture is the GAIN model on purpose. The FIR model has two samples of history that
        // survive NAM's own Reset across a prepare() — verified cross-tree, it does the same on main
        // with the Catmull-Rom kernel — so a FIR-based version of this test would fail for a reason
        // that is not the resampler's and is not this PR's to fix. A memoryless model removes that
        // term and leaves only the question being asked.
        {
            const auto json = gainModel();
            auto renderAt = [&json] (double firstRate, bool runAudio, double finalRate)
            {
                nam::NamStage st;
                st.prepare (firstRate, 64);
                if (! load (st, json)) return std::vector<float>{};
                st.prepare (firstRate, 64);
                if (runAudio)
                {
                    std::vector<float> warm (512, 0.3f);
                    float* wio[1] { warm.data() };
                    felitronics::test::run (st.process (wio, 1, 512, false));
                }
                st.prepare (finalRate, 64);
                std::vector<float> in (1024, 0.0f);
                for (int i = 0; i < 1024; ++i)
                    in[(std::size_t) i] = 0.5f * (float) std::sin (2.0 * 3.14159265358979323846 * 6000.0 * i / finalRate);
                return runMono (st, in, false);
            };
            struct Case { double first; bool audio; double final_; const char* what; };
            for (const Case c : { Case { 44100.0, true,  48000.0, "ran at 44.1 kHz, then re-prepared at 48" },
                                  Case { 44100.0, true,  44100.0, "ran at 44.1 kHz, then re-prepared at 44.1" },
                                  Case { 44100.0, true,  96000.0, "ran at 44.1 kHz, then re-prepared at 96" },
                                  Case { 48000.0, false, 44100.0, "loaded at 48 kHz, then prepared at 44.1" } })
            {
                const auto viaOther = renderAt (c.first, c.audio, c.final_);
                const auto virginAt = renderAt (c.final_, false, c.final_);
                bool identical = (! viaOther.empty() && viaOther.size() == virginAt.size());
                for (std::size_t i = 0; identical && i < viaOther.size(); ++i)
                    identical = (std::memcmp (&viaOther[i], &virginAt[i], sizeof (float)) == 0);
                test::ok (identical,
                          std::string ("a stage that ") + c.what + " renders bit-identically to one "
                          "prepared there directly — reset() re-DESIGNS the kernel, it does not top it up");
            }
        }
    }

    test::group ("the rate-match COST, measured through the real plumbing (not a replica of it)");
    {
        // felitronics_core_streamresampler_lptv_tests pins the round-trip carrier and worst phase, but it
        // builds its own copy of NamStage::processChannel's call pattern. A crew round pointed out what
        // that misses: change the PRIMING here — one extra produceExact pad at startup, a different
        // capacity, a reordered feed — and the core suite stays green while the shipped carrier moves by
        // up to 8.7 dB, because the cascade's coherent gain depends on how the two stages are aligned.
        // So the same number is measured once more through the REAL stage, with a unity model in it.
        //
        // The class is linear, so cos and sin through two identical stages combine into the response to a
        // complex exponential — a per-sample complex gain, no bucketing. The window is a whole number of
        // the resampler's 147-sample modulation periods, which is what makes the mean the coherent term.
        auto carrierDb = [] (double f, double hostSR, int block)
        {
            nam::NamStage c, s2;
            c.prepare (hostSR, block); s2.prepare (hostSR, block);
            const auto json = gainModel();
            if (! load (c, json) || ! load (s2, json)) return 1.0e9;
            c.prepare (hostSR, block); s2.prepare (hostSR, block);
            const double W = 2.0 * 3.14159265358979323846 * f / hostSR;
            const int skip = 20000, span = 147 * 400;
            double re = 0.0, im = 0.0; int cnt = 0;
            std::vector<float> bc ((std::size_t) block), bs ((std::size_t) block);
            for (int off = 0; off < skip + span + block; off += block)
            {
                for (int i = 0; i < block; ++i)
                {
                    bc[(std::size_t) i] = (float) std::cos (W * (off + i));
                    bs[(std::size_t) i] = (float) std::sin (W * (off + i));
                }
                float* ic[1] { bc.data() }; float* is[1] { bs.data() };
                felitronics::test::run (c.process (ic, 1, block, false));
                felitronics::test::run (s2.process (is, 1, block, false));
                for (int i = 0; i < block; ++i)
                {
                    const int m = off + i;
                    if (m < skip || cnt >= span) continue;
                    const double cr = std::cos (-W * m), sr = std::sin (-W * m);
                    const double a = (double) bc[(std::size_t) i], b = (double) bs[(std::size_t) i];
                    re += a * cr - b * sr;                  // (a + i b) * e^{-i W m}
                    im += a * sr + b * cr;
                    ++cnt;
                }
            }
            return 20.0 * std::log10 (std::max (std::hypot (re, im) / (double) cnt, 1e-30));
        };

        // liveness: at the model's own rate there is no resampler, so the same instrument must read 0.00.
        test::approx (carrierDb (17640.0, 48000.0, 64), 0.0, 0.01,
                      "precondition: at 48 kHz the instrument reads 0.00 dB — there is no resampler to read");
        // 🔴 THE NUMBERS IN THIS TABLE ARE THE ACCEPTANCE OF P34, measured where it actually ships.
        // `was` is what this same instrument read through the same plumbing with the Catmull-Rom cubic;
        // `want` is the sinc. The core suite reproduces both to 1e-6 with its own replica of the call
        // pattern, and THAT agreement is a second claim worth having: it says the priming in
        // configureRates and the priming in the replica are the same priming.
        struct Row { double f, was, want; };
        for (const Row r : { Row {10000.0, -0.64, -0.000054}, Row {15000.0, -2.59, +0.000003},
                             Row {17640.0, -4.17, +0.000160}, Row {20000.0, -5.48, -0.013301} })
        {
            const double got = carrierDb (r.f, 44100.0, 64);
            std::printf ("      %5.0f Hz through the real NamStage at 44.1 kHz: %+9.6f dB (cubic read %.2f)\n",
                         r.f, got, r.was);
            test::approx (got, r.want, 0.002,
                          std::to_string ((int) r.f) + " Hz: the SHIPPED stage costs what the core suite says");
            test::ok (std::fabs (got) < std::fabs (r.was) * 0.01,
                      std::to_string ((int) r.f) + " Hz: …and that is at least 40 dB less carrier droop than "
                      "the cubic cost at the same point of the same chain");
        }
    }

    test::group ("96 kHz host rate matching");
    {
        nam::NamStage stage;
        stage.prepare (96000.0, 256);
        const auto json = gainModel();
        test::ok (load (stage, json), "48 kHz model loads on a 96 kHz host");
        // `delayInputSamples() * 3.0` used to stand here with a comment deriving D·(1 + 96000/48000)
        // = 3D. That arithmetic is only true while the two legs of the round trip are the same length,
        // and the open kTaps-scaling item makes them different — so the test would have kept passing
        // through a change it exists to notice. Ask instead.
        const int kGeo96 = nam::NamStage::rateMatch (96000.0, 48000.0).latencySamples;
        test::ok (stage.latencySamples() == kGeo96,
                  "a prepared stage reports what rateMatch() says for its rates (" + std::to_string (kGeo96)
                  + "), and nothing here retypes the composition");

        std::vector<float> block (256);
        bool finite = true;
        float peak = 0.0f;
        for (int pass = 0; pass < 20; ++pass)
        {
            for (int i = 0; i < 256; ++i)
                block[(std::size_t) i] = 0.25f * std::sin (0.03f * (float) (pass * 256 + i));
            float* io[1] { block.data() };
            felitronics::test::run (stage.process (io, 1, 256, false));
            finite = finite && allFinite (block);
            for (float value : block) peak = std::max (peak, std::fabs (value));
        }
        test::ok (finite, "resampled processing stays finite");
        test::ok (peak > 0.05f && peak < 1.0f, "resampled unity model output remains sane");

        auto processLayout = [&stage] (int channels)
        {
            std::vector<float> left (256, 0.1f), right (256, -0.1f);
            float* io[2] { left.data(), right.data() };
            felitronics::test::run (stage.process (io, channels, 256, false));
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
        test::ok (stereoLatency == kGeo96 && stage.latencySamples() == kGeo96,
                  "layout re-prepare leaves the pinned host latency stable");
    }

    test::group ("the load in two halves: prepared anywhere, installed on the message thread");
    {
        nam::NamStage stage;
        stage.prepare (48000.0, 512);
        const auto json = firModel();
        auto prepared = nam::NamStage::prepareModel (json.data(), json.size(), 48000.0, 512);
        test::ok (prepared != nullptr, "prepared with no stage in sight");
        test::ok (! stage.hasModel(), "…and the stage has nothing yet");
        nam::NamStage ref;                             // the reference: the same bytes, the one-call load
        ref.prepare (48000.0, 512);
        test::ok (load (ref, json), "reference loads");
        test::ok (stage.install (std::move (prepared)), "installed");
        test::ok (stage.hasModel() && prepared == nullptr, "…the stage holds it and the handle is spent");
        std::vector<float> input (64, 0.0f);
        input[0] = 1.0f;
        const auto a = runMono (stage, input, false), b = runMono (ref, input, false);
        bool same = a.size() == b.size();
        for (std::size_t i = 0; same && i < a.size(); ++i) same = std::abs (a[i] - b[i]) < 1e-7f;
        test::ok (same, "…and it is the model loadModelFromMemory gives, sample for sample");
        test::ok (stage.modelSampleRate() == 48000.0 && stage.prewarmSamples() == ref.prewarmSamples()
                      && stage.latencySamples() == ref.latencySamples(),
                  "the mirrors read the installed model");

        test::ok (! stage.install (nullptr), "a null handle is refused");
        test::ok (stage.hasModel(), "…and the stage keeps its model");

        auto other = nam::NamStage::prepareModel (json.data(), json.size(), 96000.0, 256);
        test::ok (other != nullptr && stage.install (std::move (other)), "prepared for 96k/256, installed into a 48k/512 stage");
        const auto c = runMono (stage, input, false);
        same = c.size() == b.size();
        for (std::size_t i = 0; same && i < c.size(); ++i) same = std::abs (c[i] - b[i]) < 1e-7f;
        test::ok (same && stage.latencySamples() == ref.latencySamples(),
                  "…prepared again for the stage's own numbers: the same output, the same latency");

        const std::string junk = "{not a model";
        test::ok (nam::NamStage::prepareModel (junk.data(), junk.size(), 48000.0, 512) == nullptr, "junk prepares to nothing");

        // 🔴 THIS PAIR USED TO SAY THE OPPOSITE, AND THE OLD WORDING WAS THE DEFECT'S OWN VOICE: "a 96 kHz
        // model prepares — no stage was there to judge it". It was true only because the rate contract
        // compared against a rate a stage CARRIED, and a rate a stage carries is a rate loads can move.
        // The contract reads nothing but the tag now, so the heavy half judges it and the network is
        // paid: the refusal arrives earlier, at the same place a non-mono or corrupt model is refused.
        // MEASURED, not assumed, and stated as narrowly as it was measured: every consumer in this tree
        // and every one found outside it reaches the load through `loadModelFromMemory` (orbitcab
        // `CabEngine.h:107,163`, orbit-amp `NamBench.cpp:45`) or through
        // `prepared != nullptr && install(...)` (`RigPlayer.h:1114`), and both fold a null handle and a
        // refused install into one outcome. A loader that reported "bad file" for one and "wrong rate"
        // for the other WOULD see the difference; none exists today.
        const auto json96 = firModel ("96000");
        auto wrong = nam::NamStage::prepareModel (json96.data(), json96.size(), 48000.0, 512);
        test::ok (wrong == nullptr, "a 96 kHz model no longer prepares at all — the rate contract needs "
                                    "no stage to judge it, so it is settled before the prewarm is paid");
        test::ok (! stage.install (std::move (wrong)) && stage.hasModel(),
                  "…and installing the null it produced is refused, keeping the model, as a load would");
        const auto j441 = firModel ("44100");
        test::ok (nam::NamStage::prepareModel (j441.data(), j441.size(), 48000.0, 512) == nullptr
                      && ! load (stage, j441) && stage.hasModel(),
                  "…and 44.1 kHz the same way through both entry points, with the stage untouched");

        nam::NamStage::PreparedModel fromWorker;      // the point of the split: another thread does the work
        std::thread ([&] { fromWorker = nam::NamStage::prepareModel (json.data(), json.size(), 48000.0, 512); }).join();
        test::ok (fromWorker != nullptr && stage.install (std::move (fromWorker)), "prepared on a worker thread, installed here");
        const auto d = runMono (stage, input, false);
        same = d.size() == b.size();
        for (std::size_t i = 0; same && i < d.size(); ++i) same = std::abs (d[i] - b[i]) < 1e-7f;
        test::ok (same, "…the same model again");
    }

    test::group ("law 11a: a lane that stops being fed brings nothing back with it");
    {
        // 🔴 EVERY per-lane thing this stage holds, and there are THREE of them: the network's own
        // window, and the two core::StreamResamplers that stand either side of it when the host rate is
        // not the model's. A lane the host stops handing over used to be skipped whole, so all three
        // froze and were REPLAYED on the return. Measured before the fix, worst |out| out of DIGITAL
        // SILENCE / tail in host samples, through rigplayer::RigPlayer:
        //
        //     memoryless capture   44.1k 0.518588/125 · 48k 0.000000/0 · 88.2k 0.332768/182
        //                          96k 0.500179/193 · 176.4k 0.354183/300 · 192k 0.453147/319
        //     2001-tap capture     44.1k 0.499446/1962 · 48k 0.499533/2003 · 96k 0.500179/4193
        //
        // 🔴 AND THAT IS WHY THIS SWEEPS BOTH AXES. The two halves are INDEPENDENT and each has a
        // fixture that cannot see it: at 48 kHz no rate-matcher is installed at all, so a suite pinned
        // there (RigPlayerTests was) measures only the model's memory; and a capture whose field is one
        // sample measures only the resamplers. The 48 kHz row with a long capture is the cell that was
        // missing, and it was NOT clean — 0.499533 with a 2003-sample tail.
        //
        // The ceiling is NOT this measurement. The whole return path is linear in the pre-gap input for
        // a Linear capture, so `y[n] = sum_k h_n[k]·x[-k]` and the ceiling is `A·max_n ||h_n||_1`,
        // attained by the sign pattern of the worst row. At 44.1 kHz through the player that is
        // 0.949383 (-0.45 dBFS), attained to 100.00 %, where a 220 Hz sine swept over ALL 201 leaving
        // phases reaches 55 % of it. See rigplayer's own group for the derivation.
        // `field` is the MEMORY the stage should report (taps − 1 for a Linear capture, and taps for a
        // dilated stack, which keeps NAM's own +1 of margin); `reach` is how far the impulse response
        // is measured to go, which is what proves the fixture is not blind.
        struct Shape { const char* name; std::string json; int field; int reach; };
        const Shape shapes[] {
            // A Linear capture DECLARES its field and NAM answers zero for it — the path that reported
            // prewarmSamples() == 0 for a 2001-tap impulse response until detail::declaredReceptiveField.
            { "Linear delay(512)", delayModel (512), 512, 512 },
            // …and a real WaveNet with the real captures' field. A dilated tap costs the same nine
            // scalars at any distance, so 6332 samples of memory is a nine-number fixture. The weight
            // sits on the OLDEST tap deliberately: with the shipped [1,0,0,…] shape both convolution
            // taps are zero and the network's effective memory is ONE sample, whatever it declares.
            { "WaveNet field 6332", waveNetDelayModel (6331), 6332, 6331 },
        };

        for (const auto& shape : shapes)
        {
            // PRECONDITION — the fixture's memory is MEASURED, not declared. This is the check the
            // shipped WaveNet fixture would fail: it declares 6332 and forgets after one sample.
            {
                nam::NamStage stage;
                stage.prepare (48000.0, 8192);
                test::ok (load (stage, shape.json), std::string (shape.name) + " loads");
                test::ok (stage.prewarmSamples() == shape.field,
                          std::string (shape.name) + " reports its field: " + std::to_string (shape.field));
                std::vector<float> x (8192, 0.0f);
                float* io[1] { x.data() };
                x[0] = 0.5f;
                felitronics::test::run (stage.process (io, 1, 8192, false));
                int last = -1;
                for (int i = 0; i < 8192; ++i) if (x[(std::size_t) i] != 0.0f) last = i;
                test::ok (last >= shape.reach,
                          std::string ("precondition: ") + shape.name + " really remembers that far — its"
                          " impulse response reaches sample " + std::to_string (last)
                          + ", at least " + std::to_string (shape.reach));
                std::fill (x.begin(), x.end(), 0.0f);
                felitronics::test::run (stage.process (io, 1, 8192, false));
                double zero = 0.0;
                for (float v : x) zero = std::fmax (zero, (double) std::fabs (v));
                test::ok (zero == 0.0,
                          std::string ("precondition: ") + shape.name + " answers digital silence with"
                          " digital silence, so a non-zero return can only be the frozen state");
            }

            // 8000 and 22050 are not decoration: at a host rate well BELOW the model's, the ratio makes
            // the converted term small, and what keeps the DOWN leg honest is its own tap window in host
            // samples. Drop that term and these two rows are the ones that notice.
            for (const double fs : { 8000.0, 22050.0, 44100.0, 48000.0, 88200.0, 96000.0, 176400.0, 192000.0 })
            {
                // PRECONDITION — the grid straddles the rate-match gate rather than assuming it does.
                // ASK the owner; a literal here would be the restatement rule 9u exists against.
                const bool resampled = nam::NamStage::rateMatch (fs, nam::NamStage::kModelSampleRate).resampling;
                test::ok (resampled == (fs != nam::NamStage::kModelSampleRate),
                          "precondition: at " + std::to_string ((int) fs) + " Hz a rate-matcher is "
                          + (resampled ? "INSTALLED" : "not installed"));

                for (const int gapWidth : { 0, 1 })
                {
                    constexpr int kBlk = 256;
                    nam::NamStage stage;
                    stage.prepare (fs, kBlk);
                    if (! load (stage, shape.json)) { test::ok (false, "gap fixture loads"); continue; }

                    std::vector<float> l ((std::size_t) kBlk), r ((std::size_t) kBlk);
                    float* io[2] { l.data(), r.data() };
                    const int fill = (int) std::ceil ((double) shape.reach * fs / 48000.0) + 4 * kBlk;
                    double phase = 0.0, charged = 0.0;
                    for (int n = 0; n < fill; n += kBlk)
                    {
                        for (int i = 0; i < kBlk; ++i)
                        {
                            const float v = (float) (0.5 * std::sin (phase));
                            phase += 2.0 * kPi * 220.0 / fs;
                            l[(std::size_t) i] = v; r[(std::size_t) i] = v;
                        }
                        felitronics::test::run (stage.process (io, 2, kBlk, false));
                        for (float v : r) charged = std::fmax (charged, (double) std::fabs (v));
                    }
                    test::ok (charged > 0.1, "precondition: the lane under test really was playing at "
                                             + std::to_string ((int) fs) + " Hz");

                    // THE GAP. Width 1 leaves lane 0 playing (zeros) and takes lane 1 away; width 0
                    // takes both. Long enough that a lane still holding anything has time to say so.
                    const int gap = fill + 8 * kBlk;
                    for (int n = 0; n < gap; n += kBlk)
                    {
                        std::fill (l.begin(), l.end(), 0.0f); std::fill (r.begin(), r.end(), 0.0f);
                        felitronics::test::run (stage.process (io, gapWidth, kBlk, false));
                    }

                    double worst = 0.0;
                    for (int n = 0; n < fill + 4 * kBlk; n += kBlk)
                    {
                        std::fill (l.begin(), l.end(), 0.0f); std::fill (r.begin(), r.end(), 0.0f);
                        felitronics::test::run (stage.process (io, 2, kBlk, false));
                        for (float v : l) worst = std::fmax (worst, (double) std::fabs (v));
                        for (float v : r) worst = std::fmax (worst, (double) std::fabs (v));
                    }
                    test::ok (worst == 0.0,
                              std::string ("silence in, EXACT zero out after a ")
                              + (gapWidth == 1 ? "narrow" : "zero-width") + " gap — " + shape.name
                              + " at " + std::to_string ((int) fs) + " Hz");
                }
            }
        }
    }

    test::group ("law 11a: what a drain of exactly the receptive field does NOT flush");
    {
        // TWO things reach past the field, and each has a fixture that hides it.
        //
        // 1. THE PARTITIONED-FFT ENGINE. `implementation` defaults to `auto`, which is FFT past 256
        //    taps — so this is the configuration a real IR-as-NAM capture ships, not an exotic one. Its
        //    ring holds input spectra for a couple of its own blocks (256 / 512 / 1024 taps for a field
        //    of <= 2048 / <= 8192 / more, NAM v0.5.4 linear.cpp:14-17), so a lane goes on emitting after
        //    its field is clean. Only a DENSE kernel shows it: with a single tap all three
        //    implementations read zero.
        for (const double fs : { 48000.0, 44100.0 })       // …and at a rate where a rate-matcher IS installed
        for (const char* impl : { "fft", (const char*) nullptr })
        {
            nam::NamStage stage;
            stage.prepare (fs, 256);
            const auto json = denseLinearModel (2001, impl);
            const std::string what = std::string (impl != nullptr ? "\"fft\"" : "no `implementation` key")
                                   + " at " + std::to_string ((int) fs) + " Hz";
            if (! load (stage, json)) { test::ok (false, "dense Linear loads with " + what); continue; }
            std::vector<float> l (256), r (256);
            float* io[2] { l.data(), r.data() };
            double phase = 0.0, charged = 0.0;
            for (int k = 0; k < 40; ++k)
            {
                for (int i = 0; i < 256; ++i)
                {
                    const float v = (float) (0.5 * std::sin (phase));
                    phase += 2.0 * kPi * 220.0 / fs;
                    l[(std::size_t) i] = v; r[(std::size_t) i] = v;
                }
                felitronics::test::run (stage.process (io, 2, 256, false));
                for (float v : r) charged = std::fmax (charged, (double) std::fabs (v));
            }
            test::ok (charged > 0.1, "precondition: the dense kernel really sounds with " + what);
            for (int k = 0; k < 100; ++k)
            {
                std::fill (l.begin(), l.end(), 0.0f); std::fill (r.begin(), r.end(), 0.0f);
                felitronics::test::run (stage.process (io, 1, 256, false));
            }
            double worst = 0.0;
            for (int k = 0; k < 40; ++k)
            {
                std::fill (l.begin(), l.end(), 0.0f); std::fill (r.begin(), r.end(), 0.0f);
                felitronics::test::run (stage.process (io, 2, 256, false));
                for (float v : r) worst = std::fmax (worst, (double) std::fabs (v));
            }
            test::ok (worst == 0.0, "a dense 2001-tap capture with " + what + " returns EXACT zero — the"
                                    " field alone leaves 1.909e-08 on 46 samples");
        }

        // 2. A RECURRENT CELL, where there is no flush length at all — the drain is NAM's own
        //    half-second heuristic and nothing more. What IS closed here is the heuristic's own hole: it
        //    is `0.5 x the model's TAG`, so an untagged LSTM answers ONE sample and would drain for one.
        //    The gate is therefore "the untagged model behaves as the tagged one", not "exact zero":
        //    exact zero is unreachable for a cell whose time constant is 22 000 samples, and saying
        //    otherwise would be a promise the arithmetic cannot keep. Measured before the floor:
        //    untagged 0.499222, tagged 0.419413, a lane clocked throughout 0.022842.
        double leak[2] { 0.0, 0.0 };
        int    reported[2] { 0, 0 };
        for (int tagged = 0; tagged < 2; ++tagged)
        {
            nam::NamStage stage;
            stage.prepare (48000.0, 256);
            const auto json = slowLstmModel (tagged != 0);
            if (! load (stage, json)) { test::ok (false, "slow LSTM loads"); continue; }
            reported[tagged] = stage.prewarmSamples();
            std::vector<float> l (256), r (256);
            float* io[2] { l.data(), r.data() };
            double phase = 0.0;
            for (int k = 0; k < 400; ++k)                    // ~2.1 s: the cell settles (tau ~ 22 000)
            {
                for (int i = 0; i < 256; ++i)
                {
                    const float v = (float) (0.5 * std::sin (phase));
                    phase += 2.0 * kPi * 220.0 / 48000.0;
                    l[(std::size_t) i] = v; r[(std::size_t) i] = v;
                }
                felitronics::test::run (stage.process (io, 2, 256, false));
            }
            for (int k = 0; k < 375; ++k)                    // 2 s of gap, lane 1 away
            {
                std::fill (l.begin(), l.end(), 0.0f); std::fill (r.begin(), r.end(), 0.0f);
                felitronics::test::run (stage.process (io, 1, 256, false));
            }
            std::fill (l.begin(), l.end(), 0.0f); std::fill (r.begin(), r.end(), 0.0f);
            felitronics::test::run (stage.process (io, 2, 256, false));
            for (float v : r) leak[tagged] = std::fmax (leak[tagged], (double) std::fabs (v));
        }
        test::ok (reported[0] == 1 && reported[1] == 24000,
                  "precondition: NAM answers 1 for the untagged LSTM and 24000 for the tagged one, so the"
                  " two would drain by a factor of 24000 if the drain trusted that number");
        test::ok (leak[0] == leak[1],
                  "…and they leak the SAME, because the drain floors at half a second of the RUN rate: "
                  + std::to_string (leak[0]) + " against " + std::to_string (leak[1]));
        test::ok (leak[1] < 0.45,
                  "…and that is a bound on the heuristic, not on the cell: " + std::to_string (leak[1])
                  + " where a lane clocked through the whole gap reads 0.022842 — no finite drain closes"
                    " a recurrent state, and this test says so rather than promising zero");
    }

    test::group ("law 11a: a loaded CONTAINER, and the partition the small fixture cannot reach");
    {
        // 🔴 A CONTAINER IS ASKED THROUGH — and the synthetic JSON tests in ReceptiveFieldTests cannot
        // establish that, because they never LOAD anything. A crew round measured both holes on a real
        // model: a SlimmableContainer wrapping the dense 2001-tap capture returned 1.48e-08 after a gap
        // (its Linear submodel's FFT ring, uncharged because the top-level architecture is not Linear),
        // and one wrapping the untagged LSTM went from 0.419115 to 0.499275 (the recurrent floor, same
        // reason). Both are runtime classification, so both are pinned here rather than there.
        // A submodel is a WHOLE model spec (NAM v0.5.4 `container.cpp:157-166`: "has architecture,
        // config, weights, etc."), and the container's own `weights` array is empty.
        const auto container = [] (const std::string& inner) {
            return R"({"version":"0.5.0","architecture":"SlimmableContainer","config":{"submodels":[)"
                   R"({"max_value":1.0,"model":)" + inner + R"(}]},"weights":[],"sample_rate":48000})";
        };
        const auto json = container (denseLinearModel (2001, nullptr));
        nam::NamStage stage;
        stage.prepare (48000.0, 256);
        if (! load (stage, json))
            test::ok (false, "a SlimmableContainer of a dense Linear capture loads");
        else
        {
            test::ok (stage.prewarmSamples() == 2000, "…and the container reports its submodel's field: "
                                                      + std::to_string (stage.prewarmSamples()));
            std::vector<float> l (256), r (256);
            float* io[2] { l.data(), r.data() };
            double phase = 0.0, charged = 0.0;
            for (int k = 0; k < 40; ++k)
            {
                for (int i = 0; i < 256; ++i)
                {
                    const float v = (float) (0.5 * std::sin (phase));
                    phase += 2.0 * kPi * 220.0 / 48000.0;
                    l[(std::size_t) i] = v; r[(std::size_t) i] = v;
                }
                felitronics::test::run (stage.process (io, 2, 256, false));
                for (float v : r) charged = std::fmax (charged, (double) std::fabs (v));
            }
            test::ok (charged > 0.1, "precondition: the container really sounds");
            for (int k = 0; k < 100; ++k)
            {
                std::fill (l.begin(), l.end(), 0.0f); std::fill (r.begin(), r.end(), 0.0f);
                felitronics::test::run (stage.process (io, 1, 256, false));
            }
            double worst = 0.0; bool finite = true;
            for (int k = 0; k < 40; ++k)
            {
                std::fill (l.begin(), l.end(), 0.0f); std::fill (r.begin(), r.end(), 0.0f);
                felitronics::test::run (stage.process (io, 2, 256, false));
                for (float v : r) { worst = std::fmax (worst, (double) std::fabs (v)); finite = finite && std::isfinite (v); }
            }
            test::ok (finite && worst == 0.0, "a CONTAINER's submodel is charged its ring too — the field"
                                              " alone left 1.48e-08");
        }

        // THE LARGE PARTITION. NAM's Linear FFT block is 256 / 512 / 1024 taps for a field of <= 2048 /
        // <= 8192 / more (v0.5.4 linear.cpp:14-31), so the 2001-tap fixture above only ever exercises
        // the SMALLEST one: a crew round shrank the charge from 2·1024 to 512 and it survived. A capture
        // past 8192 taps runs the largest block, where the tail beyond the field reaches 2·1024 − 2.
        {
            nam::NamStage stage;
            stage.prepare (48000.0, 256);
            const auto big = denseLinearModel (8193, nullptr);
            if (! load (stage, big)) { test::ok (false, "an 8193-tap dense capture loads"); }
            else
            {
                std::vector<float> l (256), r (256);
                float* io[2] { l.data(), r.data() };
                double phase = 0.0, charged = 0.0;
                for (int k = 0; k < 80; ++k)
                {
                    for (int i = 0; i < 256; ++i)
                    {
                        const float v = (float) (0.5 * std::sin (phase));
                        phase += 2.0 * kPi * 220.0 / 48000.0;
                        l[(std::size_t) i] = v; r[(std::size_t) i] = v;
                    }
                    felitronics::test::run (stage.process (io, 2, 256, false));
                    for (float v : r) charged = std::fmax (charged, (double) std::fabs (v));
                }
                test::ok (charged > 0.1, "precondition: the 8193-tap capture sounds");
                for (int k = 0; k < 200; ++k)
                {
                    std::fill (l.begin(), l.end(), 0.0f); std::fill (r.begin(), r.end(), 0.0f);
                    felitronics::test::run (stage.process (io, 1, 256, false));
                }
                double worst = 0.0; bool finite = true;
                for (int k = 0; k < 60; ++k)
                {
                    std::fill (l.begin(), l.end(), 0.0f); std::fill (r.begin(), r.end(), 0.0f);
                    felitronics::test::run (stage.process (io, 2, 256, false));
                    for (float v : r) { worst = std::fmax (worst, (double) std::fabs (v)); finite = finite && std::isfinite (v); }
                }
                test::ok (finite && worst == 0.0, "…and the LARGEST FFT partition is covered by the same"
                                                  " 2 x 1024, which is what makes 2048 the number and not 512");
            }
        }
    }

    test::group ("law 11a: a prepare is a departure of EVERY lane, so both owe a drain afterwards");
    {
        // 🔴 THE DEBT MUST NOT BE STRANDED. `configureRates` rebuilds both rate-matchers clean and Resets
        // both networks — and a Reset does NOT clear a Linear capture's window (`Buffer::_input_buffers`
        // survive `SetMaxBufferSize`, and Linear's prewarm is the base class's zero): measured on a dense
        // 2001-tap capture, a tone then `prepare()` then digital silence at FULL width returns
        // 0.224604502320. So a host that changes its buffer size while a lane is AWAY would, if the
        // counters were cleared to zero there, hand that lane's tone back on the widen — measured
        // 0.224604502320 with `direct` and 0.072609648108 with the FFT engine before the counters were
        // charged instead. Both lanes owe a FULL drain after a prepare; that is what this pins.
        for (const double fs : { 48000.0, 44100.0 })
        {
            nam::NamStage stage;
            stage.prepare (fs, 256);
            const auto json = denseLinearModel (2001, "direct");
            if (! load (stage, json)) { test::ok (false, "dense Linear loads"); continue; }
            std::vector<float> l (256), r (256);
            float* io[2] { l.data(), r.data() };
            double phase = 0.0, charged = 0.0;
            for (int k = 0; k < 40; ++k)
            {
                for (int i = 0; i < 256; ++i)
                {
                    const float v = (float) (0.5 * std::sin (phase));
                    phase += 2.0 * kPi * 220.0 / fs;
                    l[(std::size_t) i] = v; r[(std::size_t) i] = v;
                }
                felitronics::test::run (stage.process (io, 2, 256, false));
                for (float v : r) charged = std::fmax (charged, (double) std::fabs (v));
            }
            test::ok (charged > 0.1, "precondition: the capture sounds at " + std::to_string ((int) fs) + " Hz");
            std::fill (l.begin(), l.end(), 0.0f); std::fill (r.begin(), r.end(), 0.0f);
            felitronics::test::run (stage.process (io, 1, 100, false));      // …100 samples of the drain spent
            stage.prepare (fs, 256);                                         // …and NOW the host re-prepares
            for (int k = 0; k < 100; ++k)
            {
                std::fill (l.begin(), l.end(), 0.0f); std::fill (r.begin(), r.end(), 0.0f);
                felitronics::test::run (stage.process (io, 1, 256, false));
            }
            double worst = 0.0;
            for (int k = 0; k < 40; ++k)
            {
                std::fill (l.begin(), l.end(), 0.0f); std::fill (r.begin(), r.end(), 0.0f);
                felitronics::test::run (stage.process (io, 2, 256, false));
                for (float v : r) worst = std::fmax (worst, (double) std::fabs (v));
            }
            test::ok (worst == 0.0, "a prepare() in the middle of a drain does not strand it — "
                                    + std::to_string ((int) fs) + " Hz");
        }
    }

    test::group ("law 11a: the drain's LENGTH, its lifecycle, and the shapes that hide it");
    {
        // 🔴 THE ODOMETER, because the audio cannot say this. Past the debt an absent lane's output is
        // zero whether it is still being clocked or not, so "it drains, and then it STOPS" has no
        // witness in the sound: a mutation that never decrements the debt, or rounds it up to the whole
        // chunk, or doubles it, is INAUDIBLE. A crew round ran exactly those three and all three
        // survived a suite of 960 checks. drainedSamples() is what closes them.
        {
            nam::NamStage stage;
            stage.prepare (48000.0, 256);
            // 🔴 THE DEBT IS DELIBERATELY NOT A MULTIPLE OF THE BLOCK. delayModel(514) is 515 taps, so
            // 514 samples of memory, and the drain is 514 + 2048 = 2562 against a 256-sample block: a
            // mutation that rounds each drained chunk up to the whole call then spends 2816, and the
            // odometer says so. With delayModel(512) the debt was 2560 = ten blocks exactly and that
            // mutation was invisible — measured, it survived the suite.
            const auto json = delayModel (514);
            test::ok (load (stage, json), "the odometer fixture loads");
            test::ok (stage.drainedSamples() == 0, "nothing is owed before a lane has ever played");
            std::vector<float> l (256, 0.1f), r (256, 0.1f);
            float* io[2] { l.data(), r.data() };
            for (int k = 0; k < 40; ++k) felitronics::test::run (stage.process (io, 1, 256, false));
            test::ok (stage.drainedSamples() == 0,
                      "…nor on a MONO host, where lane 1 has never carried audio and owes nothing: a debt"
                      " armed for a lane that never played costs a real WaveNet 132 ms per load for a"
                      " window NAM already zero-filled");
            felitronics::test::run (stage.process (io, 2, 256, false));
            test::ok (stage.drainedSamples() == 0, "…nor while both lanes are playing");
            // THE NUMBER, with its derivation rather than a call to the code that computes it: the host
            // rate IS the model rate here, so no rate-matcher is installed and the drain is the field
            // plus the partitioned-FFT ring a Linear capture is charged — 514 + 2·1024 = 2562 per lane,
            // and the gap below takes BOTH lanes away, so 5124. (Legal changes to either term must
            // update this literal ON PURPOSE; that is what a literal is for in an oracle.)
            //
            // …and a ZERO-LENGTH call in the middle of it spends nothing and CLEARS nothing: law 11(d)
            // is "no samples, no time, no edge", and a mutation that drops the debt there replays
            // 0.470000 on the return. It is cut in HERE, with the drain unfinished, because after the
            // debt is spent there is nothing left for it to clear.
            for (int k = 0; k < 2; ++k)
            {
                std::fill (l.begin(), l.end(), 0.0f); std::fill (r.begin(), r.end(), 0.0f);
                felitronics::test::run (stage.process (io, 0, 256, false));
            }
            felitronics::test::run (stage.process (io, 0, 0, false));
            felitronics::test::run (stage.process (io, 1, 0, false));
            for (int k = 0; k < 58; ++k)
            {
                std::fill (l.begin(), l.end(), 0.0f); std::fill (r.begin(), r.end(), 0.0f);
                felitronics::test::run (stage.process (io, 0, 256, false));
            }
            test::ok (stage.drainedSamples() == 5124,
                      "the drain is EXACTLY the field plus the FFT ring, per lane: 2 x (514 + 2048) = 5124,"
                      " read " + std::to_string (stage.drainedSamples()));
            const long long settled = stage.drainedSamples();
            for (int k = 0; k < 60; ++k) felitronics::test::run (stage.process (io, 0, 256, false));
            test::ok (stage.drainedSamples() == settled, "…and it STOPS: a hundred more blocks of gap owe nothing");

            // A REFUSED CALL MOVES NOTHING AT ALL, the debt included — law 11's "the refused call is
            // indistinguishable from one never made". A width above 2 and a negative length are the two
            // ways in, and neither may spend a sample of the drain.
            test::ok (! stage.process (io, 3, 256, false) && ! stage.process (io, 2, -1, false)
                          && ! stage.process (io, -1, 256, false),
                      "a malformed call is refused");
            test::ok (stage.drainedSamples() == settled, "…and a refused call spends none of the debt");

            // A ZERO-LENGTH CALL DOES NOT SPEND OR CLEAR THE DEBT — law 11(d) is "no samples, no time,
            // no edge", and a mutation clearing the debt there replays 0.470000 on the return.
            felitronics::test::run (stage.process (io, 0, 0, false));
            felitronics::test::run (stage.process (io, 1, 0, false));
            test::ok (stage.drainedSamples() == settled, "a zero-length call after it is no time either");

            // WIDTHS OUT OF ORDER: 2 -> 1 -> 0 -> 1 -> 2, which is what a host does when it reconfigures
            // a bus twice in a row. Lane 1 leaves at the 1, lane 0 at the 0, and each owes from ITS OWN
            // departure — the debts are per lane and do not share a clock.
            for (int k = 0; k < 4; ++k)
            {
                std::fill (l.begin(), l.end(), 0.3f); std::fill (r.begin(), r.end(), 0.3f);
                felitronics::test::run (stage.process (io, 2, 256, false));
            }
            const long long beforeMixed = stage.drainedSamples();
            felitronics::test::run (stage.process (io, 1, 256, false));    // lane 1 leaves: 256 of its debt
            felitronics::test::run (stage.process (io, 0, 256, false));    // …lane 0 too, and lane 1 goes on
            test::ok (stage.drainedSamples() == beforeMixed + 256 + 2 * 256,
                      "widths out of order: each lane owes from its OWN departure — "
                      + std::to_string (stage.drainedSamples() - beforeMixed) + " samples over three lane-blocks");

            // A SECOND DEPARTURE RE-ARMS. Feeding a lane again is what owes the next drain; a mutation
            // that arms the debt only once replays 0.5 on the second return.
            for (int k = 0; k < 20; ++k)
            {
                std::fill (l.begin(), l.end(), 0.2f); std::fill (r.begin(), r.end(), 0.2f);
                felitronics::test::run (stage.process (io, 2, 256, false));
            }
            const long long beforeSecond = stage.drainedSamples();
            for (int k = 0; k < 60; ++k)
            {
                std::fill (l.begin(), l.end(), 0.0f); std::fill (r.begin(), r.end(), 0.0f);
                felitronics::test::run (stage.process (io, 0, 256, false));
            }
            test::ok (stage.drainedSamples() - beforeSecond == settled,
                      "…and a SECOND departure owes the same again, exactly: "
                      + std::to_string (stage.drainedSamples() - beforeSecond));
        }

        // A LOAD LANDING WHILE A LANE IS AWAY. The debt belongs to the backend, and a load REPLACES it,
        // so the arriving instance owes nothing — its window is the zeros NAM filled it with. What must
        // not happen is the arriving model speaking the DEPARTED one's audio on the widen, and what must
        // also not happen is the new backend inheriting a debt it cannot owe.
        {
            nam::NamStage stage;
            stage.prepare (48000.0, 256);
            test::ok (load (stage, delayModel (514)), "the first capture loads");
            std::vector<float> l (256, 0.4f), r (256, 0.4f);
            float* io[2] { l.data(), r.data() };
            for (int k = 0; k < 20; ++k) felitronics::test::run (stage.process (io, 2, 256, false));
            std::fill (l.begin(), l.end(), 0.0f); std::fill (r.begin(), r.end(), 0.0f);
            felitronics::test::run (stage.process (io, 1, 256, false));      // lane 1 leaves, mid-debt
            const long long owed = stage.drainedSamples();
            test::ok (owed > 0, "precondition: lane 1 really was draining when the load landed");
            test::ok (load (stage, delayModel (300)), "…and a DIFFERENT capture lands while it is away");
            for (int k = 0; k < 40; ++k)
            {
                std::fill (l.begin(), l.end(), 0.0f); std::fill (r.begin(), r.end(), 0.0f);
                felitronics::test::run (stage.process (io, 1, 256, false));
            }
            test::ok (stage.drainedSamples() == owed,
                      "the arriving backend owes NOTHING — its window is the zeros it was built with, and"
                      " the departed one's debt went with it");
            double worst = 0.0; bool finite = true;
            for (int k = 0; k < 40; ++k)
            {
                std::fill (l.begin(), l.end(), 0.0f); std::fill (r.begin(), r.end(), 0.0f);
                felitronics::test::run (stage.process (io, 2, 256, false));
                for (float v : r) { worst = std::fmax (worst, (double) std::fabs (v)); finite = finite && std::isfinite (v); }
            }
            test::ok (finite && worst == 0.0, "…and the widen brings back FINITE silence, not the capture"
                                              " that left");
        }

        // SHAPES THE MAIN GROUP CANNOT SEE, each one a mutation that survived it:
        //   · call lengths that are not the prepared block, including 1 and a remainder — a drain that
        //     skips chunks shorter than the block replays 0.5 at blocks 1/17/63/255;
        //   · a NON-UNITY makeup, because a drain skipped whenever the gain differs replays 1.0 with
        //     `normalize` on and a loudness tag two doublings away from the reference;
        //   · FINITENESS, because std::fmax IGNORES a NaN — writing NaNs into the whole first returning
        //     chunk passed 960 checks, since a peak taken with fmax stays 0.
        for (const int call : { 1, 17, 63, 255, 256, 257, 1000 })
            for (const bool normalise : { false, true })
            {
                nam::NamStage stage;
                stage.prepare (44100.0, 256);
                // -24 dB of tagged loudness against the -18 dB reference is a makeup of +6 dB, so the
                // drained lane and the live one are NOT running the same gain — the mutation this is for
                // skips the drain exactly when they differ.
                std::string taps = "[1.0";
                for (int i = 0; i < 64; ++i) taps += ",0.0";          // 65 taps: 64 samples of memory
                const auto json = R"({"version":"0.5.0","architecture":"Linear","config":)"
                                  R"({"receptive_field":65,"bias":false,"implementation":"direct"},)"
                                  R"("weights":)" + taps + R"(],"sample_rate":48000,)"
                                  R"("metadata":{"loudness":-24.0}})";
                if (! load (stage, json)) { test::ok (false, "the makeup fixture loads"); continue; }
                std::vector<float> l ((std::size_t) call), r ((std::size_t) call);
                float* io[2] { l.data(), r.data() };
                double phase = 0.0, charged = 0.0;
                for (int n = 0; n < 4096; n += call)
                {
                    for (int i = 0; i < call; ++i)
                    {
                        const float v = (float) (0.4 * std::sin (phase));
                        phase += 2.0 * kPi * 220.0 / 44100.0;
                        l[(std::size_t) i] = v; r[(std::size_t) i] = v;
                    }
                    felitronics::test::run (stage.process (io, 2, call, normalise));
                    for (float v : r) charged = std::fmax (charged, (double) std::fabs (v));
                }
                test::ok (charged > 0.1, "precondition: the makeup fixture sounds at call length "
                                         + std::to_string (call));
                for (int n = 0; n < 16384; n += call)
                {
                    std::fill (l.begin(), l.end(), 0.0f); std::fill (r.begin(), r.end(), 0.0f);
                    felitronics::test::run (stage.process (io, 1, call, normalise));
                }
                double worst = 0.0; bool finite = true;
                for (int n = 0; n < 8192; n += call)
                {
                    std::fill (l.begin(), l.end(), 0.0f); std::fill (r.begin(), r.end(), 0.0f);
                    felitronics::test::run (stage.process (io, 2, call, normalise));
                    for (float v : r) { worst = std::fmax (worst, (double) std::fabs (v)); finite = finite && std::isfinite (v); }
                    for (float v : l) { worst = std::fmax (worst, (double) std::fabs (v)); finite = finite && std::isfinite (v); }
                }
                test::ok (finite && worst == 0.0,
                          std::string ("silence in, FINITE exact zero out at call length ")
                          + std::to_string (call) + (normalise ? " with the makeup ON" : " with it off"));
            }
    }

    test::group ("process is RT no-alloc");
    {
        // 🔴 TWO RATES, AND THE SECOND ONE IS THE WHOLE POINT. This test prepared only at 48 kHz, where
        // `resampling` is false and processChannel takes its early branch — so the single check that
        // NamStage::process never allocates was STRUCTURALLY BLIND to the resampler, before this change
        // and after it. A crew round confirmed the hole by mutation: a resampler that allocates a copy
        // of its table inside every decimating callback survived the entire suite. 44.1 kHz is the
        // configuration a live rig actually runs, and it is the one where the table is read.
        //
        // The FIRST call after a prepare is included in the measured window on purpose: the table is
        // built in reset(), but any lazy initialisation anywhere would land exactly there.
        for (const double rate : { 48000.0, 44100.0 })
        {
            nam::NamStage stage;
            stage.prepare (rate, 512);
            const auto json = gainModel();
            test::ok (load (stage, json), "RT fixture model loads at " + std::to_string ((int) rate));
            std::vector<float> left (512, 0.2f), right (512, -0.15f);
            float* io[2] { left.data(), right.data() };
            felitronics::test::run (stage.process (io, 2, 512, false));    // warm every process-reachable container
            const long before = g_allocs.load (std::memory_order_relaxed);
            felitronics::test::run (stage.process (io, 2, 512, false));
            felitronics::test::run (stage.process (io, 2, 512, true));
            test::okNoAlloc (g_allocs.load (std::memory_order_relaxed) == before,
                             "NamStage::process performs no heap allocation at "
                             + std::to_string ((int) rate) + " Hz"
                             + (rate == 48000.0 ? " (no resampler in the path)" : " (resampler ACTIVE)"));
        }

        // 🔴 AND THE FIRST CALL, WHICH THIS GROUP USED TO WARM AWAY. Every row here processed a block
        // before starting the counter — "warm every process-reachable container" — so the one place NAM
        // grows a buffer was structurally invisible: `Buffer::_update_buffers_` resizes its per-channel
        // window on DEMAND inside process(), and `Reset` pre-grows it only through prewarm(), which runs
        // `GetPrewarmSamples()` samples — zero for a Linear capture. Measured before the fix: 4
        // allocations in the first width-1 call, two per instance, and instance 1's were NEW, because
        // the drain is the first thing that ever touched it on a mono host. Counted from the very first
        // call now, at both widths and both rates, and for both engines a Linear capture can pick.
        for (const double rate : { 48000.0, 44100.0 })
            for (const char* impl : { "direct", (const char*) nullptr })
                for (const int width : { 2, 1, 0 })
                {
                    nam::NamStage stage;
                    stage.prepare (rate, 512);
                    const auto json = impl != nullptr ? delayModel (512) : denseLinearModel (2001, nullptr);
                    test::ok (load (stage, json), "first-call fixture loads");
                    std::vector<float> left (512, 0.2f), right (512, -0.15f);
                    float* io[2] { left.data(), right.data() };
                    const long before = g_allocs.load (std::memory_order_relaxed);
                    felitronics::test::run (stage.process (io, width, 512, false));
                    felitronics::test::run (stage.process (io, width, 512, false));
                    test::okNoAlloc (g_allocs.load (std::memory_order_relaxed) == before,
                                     std::string ("the FIRST call after a prepare allocates nothing — width ")
                                     + std::to_string (width) + ", "
                                     + (impl != nullptr ? "direct" : "the FFT engine") + ", "
                                     + std::to_string ((int) rate) + " Hz");
                }

        // 🔴 AND THE THIRD BRANCH: THE DRAIN. A lane the host stops handing over is now fed digital
        // silence through the same processChannel, which is a code path the two rows above never enter —
        // exactly the blindness the 44.1 kHz row was added for, one branch further in. The capture has a
        // 512-sample memory on purpose, so the drain is long enough to be running throughout the
        // measured window rather than finishing inside the first block. Counted, not read.
        for (const double rate : { 48000.0, 44100.0 })
        {
            nam::NamStage stage;
            stage.prepare (rate, 512);
            const auto json = delayModel (512);
            test::ok (load (stage, json), "drain fixture model loads at " + std::to_string ((int) rate));
            std::vector<float> left (512, 0.2f), right (512, -0.15f);
            float* io[2] { left.data(), right.data() };
            const long before = g_allocs.load (std::memory_order_relaxed);   // …counted from the FIRST drain
            felitronics::test::run (stage.process (io, 2, 512, false));
            felitronics::test::run (stage.process (io, 0, 512, false));
            felitronics::test::run (stage.process (io, 1, 512, false));    // lane 1 drains beside a live lane 0
            felitronics::test::run (stage.process (nullptr, 0, 512, false));   // …and with no buffers at all
            felitronics::test::run (stage.process (io, 2, 512, false));
            test::okNoAlloc (g_allocs.load (std::memory_order_relaxed) == before,
                             "…nor when an absent lane is being DRAINED at "
                             + std::to_string ((int) rate) + " Hz"
                             + (rate == 48000.0 ? " (no resampler in the path)" : " (resampler ACTIVE)"));
        }
    }

    return test::report();
}

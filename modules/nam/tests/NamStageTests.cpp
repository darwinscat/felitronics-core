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
            // ⚠️ WHAT THIS LINE IS AND IS NOT, said plainly because it changed meaning. Both sides now
            // come from the same function — the stage reports rateMatch(), and `geo` asks core for the
            // fractional value rateMatch() rounds — so this can no longer catch a wrong composition.
            // It is a statement about ROUNDING ONLY: that the reporting step is a round-to-nearest and
            // not a floor, a ceil, or a truncation, across every rate including the exact halves. That
            // is worth keeping (a `ceil` here would pass a ≤1.0 bound and fail this one), and the group
            // above is where a wrong composition is caught, by measurement.
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
            // 🔴 EXACTLY 0.5 — the boundary itself, which nothing tested. A mutation flipping `> 0.5`
            // to `>= 0.5` survived the whole suite because 48000.4 and 48000.6 straddle the edge
            // without standing on it, and the gate's own predicate is only visible AT it.
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
        // repository — because every call site in the tree passes 48000, and every model that can go
        // live runs at 48000. The function's second argument was, in effect, untested. The same round
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

        // rateMatch: the three facts, including the two an instance cannot reach — an unknown rate and
        // a rate no stage would accept.
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

        const auto json96 = firModel ("96000");
        auto wrong = nam::NamStage::prepareModel (json96.data(), json96.size(), 48000.0, 512);
        test::ok (wrong != nullptr, "a 96 kHz model prepares — no stage was there to judge it");
        test::ok (! stage.install (std::move (wrong)) && stage.hasModel(),
                  "…and a stage running at 48 kHz refuses it, keeping its model, as a load would");

        nam::NamStage::PreparedModel fromWorker;      // the point of the split: another thread does the work
        std::thread ([&] { fromWorker = nam::NamStage::prepareModel (json.data(), json.size(), 48000.0, 512); }).join();
        test::ok (fromWorker != nullptr && stage.install (std::move (fromWorker)), "prepared on a worker thread, installed here");
        const auto d = runMono (stage, input, false);
        same = d.size() == b.size();
        for (std::size_t i = 0; same && i < d.size(); ++i) same = std::abs (d[i] - b[i]) < 1e-7f;
        test::ok (same, "…the same model again");
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
    }

    return test::report();
}

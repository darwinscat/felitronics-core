// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// Golden coverage for the product-level cab wrapper: reference-unity gain math, direct-convolution
// parity, mono broadcast / true-stereo publication, verbatim reverb loading AT THE HOST'S OWN RATE, and
// latest-wins retry. And (P67) the rate match with its tolerance, a one-tap IR off the host rate, and
// loads that stage nothing. And (P68) the rate factor an un-normalized load carries when it IS resampled:
// the same convolution gain at every host rate, the gain reported for it, and what the normalized path
// reports so that "the normalized path did not move" is a measurement and not a promise.

#include <felitronics_test.h>
#include <felitronics/convolution/CabConvolver.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdio>
#include <limits>
#include <string>
#include <vector>

using felitronics::convolution::CabConvolver;

namespace
{
struct Lcg
{
    unsigned long long state;
    float next()
    {
        state = state * 6364136223846793005ULL + 1442695040888963407ULL;
        return (float) ((state >> 40) & 0xFFFFFFULL) / 8388608.0f - 1.0f;
    }
};

std::vector<float> directConvolve (const std::vector<float>& input, const std::vector<float>& ir)
{
    std::vector<float> out (input.size(), 0.0f);
    for (std::size_t n = 0; n < input.size(); ++n)
    {
        double sum = 0.0;
        const std::size_t taps = std::min (ir.size(), n + 1);
        for (std::size_t k = 0; k < taps; ++k)
            sum += (double) input[n - k] * (double) ir[k];
        out[n] = (float) sum;
    }
    return out;
}

double maxDifference (const std::vector<float>& a, const std::vector<float>& b)
{
    if (a.size() != b.size())
        return std::numeric_limits<double>::infinity();

    double result = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i)
    {
        if (! std::isfinite (a[i]) || ! std::isfinite (b[i]))
            return std::numeric_limits<double>::infinity();
        result = std::max (result, std::fabs ((double) a[i] - (double) b[i]));
    }
    return result;
}

void pumpCrossfade (CabConvolver& convolver, int channels = 2)
{
    std::vector<float> left (512, 0.0f), right (512, 0.0f);
    float* io[2] { left.data(), right.data() };
    for (int i = 0; i < 12; ++i)
        felitronics::test::run (convolver.process (io, channels, 512));    // 6144 samples > the fixed 50 ms fade at 48 kHz
}

void renderVariableBlocks (CabConvolver& convolver, std::vector<float>& left, std::vector<float>& right)
{
    static constexpr int blocks[] { 1, 17, 64, 128, 333, 512 };
    std::size_t pos = 0;
    int blockIndex = 0;
    while (pos < left.size())
    {
        const int remaining = (int) (left.size() - pos);
        const int count = std::min (blocks[blockIndex % 6], remaining);
        float* io[2] { left.data() + pos, right.data() + pos };
        felitronics::test::run (convolver.process (io, 2, count));
        pos += (std::size_t) count;
        ++blockIndex;
    }
}

std::vector<float> decayingIr (int length)
{
    Lcg random { 12345 };
    std::vector<float> ir ((std::size_t) length);
    for (int i = 0; i < length; ++i)
    {
        const float envelope = std::exp (-3.0f * (float) i / (float) length);
        ir[(std::size_t) i] = envelope * (0.5f * random.next() + (i == 0 ? 1.0f : 0.0f));
    }
    return ir;
}

// Independent reference for the loader's gain: direct double-precision DTFT on a linear grid,
// deliberately unlike the shipping real-FFT implementation. Agreement pins the formula itself.
// `bins` samples [0, sampleRate/2]; its DEFAULT is a bin count, so the spacing in HERTZ widens with the
// rate, and the weight w(f) = 1/(1 + (f/2 kHz)^2) concentrates the integrand below a couple of kHz where
// that spacing lives. At 48 kHz the two grids (this one and the loader's zero-padded FFT) agree to 0.02 dB
// and the sibling check allows 0.15; at 192 kHz the same 4096 bins are four times coarser in hertz and the
// disagreement grows to 0.18 dB, which is not the loader drifting. A caller comparing across rates passes
// a count PROPORTIONAL to the rate, so the grid is the same one in hertz at every rate.
double referenceNormalizationGainDb (const std::vector<float>& ir, double sampleRate, int bins = 4096)
{
    constexpr double pi = 3.14159265358979323846;
    double numerator = 0.0;
    double denominator = 0.0;
    for (int k = 1; k <= bins; ++k)
    {
        const double frequency = 0.5 * sampleRate * (double) k / (double) bins;
        const auto step = std::exp (std::complex<double> (0.0, -2.0 * pi * frequency / sampleRate));
        std::complex<double> phase (1.0, 0.0);
        std::complex<double> response (0.0, 0.0);
        for (float tap : ir)
        {
            response += (double) tap * phase;
            phase *= step;
        }
        const double ratio = frequency / CabConvolver::kIrRefShapeHz;
        const double weight = 1.0 / (1.0 + ratio * ratio);
        numerator += weight * std::norm (response);
        denominator += weight;
    }
    return -10.0 * std::log10 (std::max (1.0e-12, numerator / denominator));
}

// P68 — THE ORACLE FOR THE RATE FACTOR, and it is computed OUTSIDE the thing it measures: the magnitude
// response at a PHYSICAL frequency, straight from the definition, sharing not a line with the resampler
// or with the loader. Two IRs at two different rates are comparable through this and through nothing
// else — a bin index means a different frequency at each rate, and a broadband energy sum would average
// away exactly the frequency-dependent residue this has to see (so the assertions below read it at five
// separate frequencies rather than once over a band).
double magnitudeAt (const std::vector<float>& ir, double sampleRate, double hz)
{
    constexpr double pi = 3.14159265358979323846;
    const auto step = std::exp (std::complex<double> (0.0, -2.0 * pi * hz / sampleRate));
    std::complex<double> phase (1.0, 0.0), response (0.0, 0.0);
    for (const float tap : ir) { response += (double) tap * phase; phase *= step; }
    return std::abs (response);
}

// A fixture with THREE properties it needs and one it must not have. It carries real content at the
// frequencies the oracle reads (a DC term plus three tones), it decays so its tail edge is 60 dB down
// before the taps end, and — the one that matters here — it opens with `lead` samples of silence.
// ITS LENGTH IS A DURATION AND ITS LEAD IS A TAP COUNT, and neither is a matter of taste: the decay
// has to finish inside the taps at EVERY source rate (a fixed 4096 taps is 85 ms at 48 kHz but 43 at
// 96, which cut this tail off at -30 dB and read 1.4e-2 dB of "error" that was nothing but the step at
// the end), while the lead answers to the window's ONE-SIDED reach, `IrResampleConfig::halfTaps` = 32
// INPUT SAMPLES — the same count at every rate, and exactly where the measured loss reaches zero (a
// lone impulse at 96 -> 48 kHz loses 2.34 dB at sample 0, 0.25 dB eight samples in, 0.05 at sixteen and
// nothing at all from thirty-one on). Twice that is used below, which is slack, not a second rule.
// WITHOUT THAT LEAD the measurement is not about the rate factor at all: a resample drops the kernel's
// pre-ringing that would fall before output sample 0, which costs an IR whose sample 0 is full scale
// 2.34 dB (a lone impulse at 96 -> 48 kHz), and on THIS fixture 0.082 dB at 10 kHz with no lead against
// 2.8e-4 dB with the lead — while the worst cell of the grid that uses it reads 5.65e-4 dB. That loss is
// a property of the resample, older than P68 and untouched by it; a fixture that starts at full scale
// would charge it to the rate factor and read a hundred times the real error.
std::vector<float> rateFixture (double sampleRate, double seconds, int lead)
{
    constexpr double pi = 3.14159265358979323846;
    const int body = (int) std::lround (seconds * sampleRate);
    std::vector<float> ir ((std::size_t) (lead + body), 0.0f);
    for (int i = 0; i < body; ++i)
    {
        const double t = (double) i / sampleRate;
        const double envelope = std::exp (-(double) i / (0.0125 * sampleRate));       // 12.5 ms, in SECONDS
        ir[(std::size_t) (lead + i)] = (float) (envelope * (0.6 + 0.3 * std::cos (2.0 * pi * 800.0 * t)
                                                                + 0.25 * std::cos (2.0 * pi * 4500.0 * t)
                                                                + 0.2 * std::cos (2.0 * pi * 9000.0 * t)));
    }
    return ir;
}
} // namespace

int main()
{
    using namespace felitronics;
    std::printf ("felitronics::convolution CabConvolver tests\n");

    test::group ("reference-unity normalization math and guards");
    {
        auto loadDelta = [] (float amplitude)
        {
            std::vector<float> ir (256, 0.0f);
            ir[0] = amplitude;
            CabConvolver convolver;
            felitronics::test::run (convolver.prepare (48000.0, 512, 2, 0.05));
            const float* banks[1] { ir.data() };
            convolver.loadIR (banks, 1, (int) ir.size(), 48000.0);
            return convolver.irNormalizationGainDb();
        };

        test::approx (loadDelta (1.0f), 0.0, 0.001, "unit delta is already reference-unity");
        test::approx (loadDelta (100.0f), CabConvolver::kIrNormMinDb, 0.001,
                      "hot delta is clamped at -30 dB");
        test::approx (loadDelta (0.01f), CabConvolver::kIrNormMaxDb, 0.001,
                      "quiet-but-valid delta is clamped at +30 dB");
        test::approx (loadDelta (0.0001f), 0.0, 0.001,
                      "near-silent delta stays at unity by the -60 dB power floor");
    }

    test::group ("golden convolution and mono-broadcast parity");
    {
        const auto ir = decayingIr (257);
        CabConvolver convolver;
        felitronics::test::run (convolver.prepare (48000.0, 512, 2, 0.05));
        const float* banks[1] { ir.data() };
        convolver.loadIR (banks, 1, (int) ir.size(), 48000.0);
        test::approx (convolver.irNormalizationGainDb(), referenceNormalizationGainDb (ir, 48000.0),
                      0.15, "applied gain matches an independent weighted-DTFT reference");
        pumpCrossfade (convolver);
        test::ok (! convolver.isBusy(), "initial 50 ms publication completes");
        convolver.reset();

        constexpr int numSamples = 4096;
        Lcg leftRandom { 999 }, rightRandom { 55555 };
        std::vector<float> left ((std::size_t) numSamples), right ((std::size_t) numSamples);
        for (float& v : left)  v = 0.4f * leftRandom.next();
        for (float& v : right) v = 0.3f * rightRandom.next();
        const auto leftIn = left;
        const auto rightIn = right;
        const auto audibleIr = convolver.stagedTaps()[0];
        const auto leftReference = directConvolve (leftIn, audibleIr);
        const auto rightReference = directConvolve (rightIn, audibleIr);
        renderVariableBlocks (convolver, left, right);

        test::ok (maxDifference (left, leftReference) < 2.0e-4,
                  "left output matches naive direct convolution");
        test::ok (maxDifference (right, rightReference) < 2.0e-4,
                  "mono IR broadcasts to right and matches direct convolution");
        test::ok (CabConvolver::latencySamples() == 0, "cab path remains sample-zero-latency");
    }

    test::group ("resampled loads publish their final staged taps");
    {
        auto verifyRate = [] (double irRate, int sourceLength, int expectedLength, const char* lengthMessage,
                              const char* parityMessage)
        {
            const auto ir = decayingIr (sourceLength);
            CabConvolver convolver;
            felitronics::test::run (convolver.prepare (48000.0, 512, 2, 0.05));
            const float* banks[1] { ir.data() };
            convolver.loadIR (banks, 1, (int) ir.size(), irRate);
            const auto& staged = convolver.stagedTaps();
            test::ok (staged.size() == 1 && (int) staged[0].size() == expectedLength, lengthMessage);

            pumpCrossfade (convolver);
            convolver.reset();
            constexpr int numSamples = 3072;
            Lcg leftRandom { 76123 }, rightRandom { 991 };
            std::vector<float> left ((std::size_t) numSamples), right ((std::size_t) numSamples);
            for (float& value : left)  value = 0.25f * leftRandom.next();
            for (float& value : right) value = 0.25f * rightRandom.next();
            const auto leftReference = directConvolve (left, staged[0]);
            const auto rightReference = directConvolve (right, staged[0]);
            renderVariableBlocks (convolver, left, right);
            test::ok (maxDifference (left, leftReference) < 2.0e-4
                      && maxDifference (right, rightReference) < 2.0e-4, parityMessage);
        };

        verifyRate (96000.0, 600, 300, "96 kHz IR stages at half length on a 48 kHz host",
                    "96-to-48 resample, normalize, and publish matches direct convolution of staged taps");
        verifyRate (44100.0, 600, (int) std::llround (600.0 * 48000.0 / 44100.0),
                    "44.1 kHz IR length scales by 48/44.1 on a 48 kHz host",
                    "44.1-to-48 resample, normalize, and publish matches direct convolution of staged taps");
    }

    test::group ("prepare-time maximum IR length is an audible truncation cap");
    {
        constexpr int capSamples = 256;
        std::vector<float> ir (512, 0.0f);
        ir[0] = 1.0f;
        ir[200] = -0.375f;     // still inside the scheduled IR span
        ir[400] = 0.75f;       // retained in stagedTaps(), deliberately outside the audible schedule

        CabConvolver convolver;
        felitronics::test::run (convolver.prepare (48000.0, 128, 2, (double) capSamples / 48000.0, false));
        const float* banks[1] { ir.data() };
        convolver.loadIR (banks, 1, (int) ir.size(), 48000.0);
        test::ok (convolver.stagedTaps()[0] == ir,
                  "loader retains the full staged IR while the fixed NUPC schedule enforces the cap");
        pumpCrossfade (convolver);
        convolver.reset();

        std::vector<float> left (640, 0.0f), right (640, 0.0f);
        left[0] = 1.0f;
        right[0] = 1.0f;
        const auto input = left;
        std::vector<float> capped (ir.begin(), ir.begin() + capSamples);
        const auto reference = directConvolve (input, capped);
        renderVariableBlocks (convolver, left, right);
        test::ok (maxDifference (left, reference) < 2.0e-5
                  && maxDifference (right, reference) < 2.0e-5,
                  "convolution contains the in-cap tap and omits the staged out-of-cap tail");
    }

    test::group ("true-stereo publication preserves one common normalization gain");
    {
        const auto leftIr = decayingIr (193);
        std::vector<float> rightIr = leftIr;
        for (float& v : rightIr) v *= 0.5f;
        CabConvolver convolver;
        felitronics::test::run (convolver.prepare (48000.0, 256, 2, 0.05));
        const float* banks[2] { leftIr.data(), rightIr.data() };
        convolver.loadIR (banks, 2, (int) leftIr.size(), 48000.0);
        const auto& staged = convolver.stagedTaps();
        test::ok (staged.size() == 2, "two input banks stage as LRDiag true stereo");
        test::approx (staged[1][17] / staged[0][17], 0.5, 1.0e-6,
                      "one common gain preserves the L/R tap ratio");
        pumpCrossfade (convolver);
        convolver.reset();
        std::vector<float> left (256, 0.0f), right (256, 0.0f);
        left[0] = 1.0f;
        right[0] = 1.0f;
        float* io[2] { left.data(), right.data() };
        felitronics::test::run (convolver.process (io, 2, 256));
        bool lrDiag = true;
        for (std::size_t i = 0; i < leftIr.size(); ++i)
            lrDiag = lrDiag && std::fabs (left[i] - staged[0][i]) < 2.0e-4f
                            && std::fabs (right[i] - staged[1][i]) < 2.0e-4f;
        test::ok (lrDiag, "rendered stereo impulse uses the distinct L/R diagonal banks");
    }

    test::group ("normalize=false is a verbatim reverb path AT THE HOST'S OWN RATE");
    {
        std::vector<float> ir { 0.75f, -0.25f, 0.125f, 0.0625f };
        CabConvolver convolver;
        felitronics::test::run (convolver.prepare (48000.0, 128, 2, 0.05, false));
        const float* banks[1] { ir.data() };
        convolver.loadIR (banks, 1, (int) ir.size(), 48000.0);
        test::approx (convolver.irNormalizationGain(), 1.0, 0.0, "normalization gain is exactly one");
        test::ok (convolver.stagedTaps()[0] == ir, "staged taps are byte-verbatim");
    }

    test::group ("load during crossfade retries and coalesces latest-wins");
    {
        std::vector<float> first  { 1.0f, 0.25f, 0.0f, 0.0f };
        std::vector<float> second { 0.0f, 1.0f, 0.0f, 0.0f };
        std::vector<float> latest { 0.25f, -0.5f, 0.75f, -1.0f };
        CabConvolver convolver;
        felitronics::test::run (convolver.prepare (48000.0, 128, 2, 0.05, false));
        const float* a[1] { first.data() };
        const float* b[1] { second.data() };
        const float* c[1] { latest.data() };
        convolver.loadIR (a, 1, (int) first.size(), 48000.0);

        float left[64] {}, right[64] {};
        float* io[2] { left, right };
        felitronics::test::run (convolver.process (io, 2, 64));       // begins the first publication's 50 ms crossfade
        convolver.loadIR (b, 1, (int) second.size(), 48000.0);
        test::ok (convolver.hasPending(), "load is retained when NUPC rejects it mid-crossfade");
        convolver.loadIR (c, 1, (int) latest.size(), 48000.0);
        test::ok (convolver.hasPending() && convolver.stagedTaps()[0] == latest,
                  "a newer rejected load replaces the pending taps");

        pumpCrossfade (convolver);
        test::ok (convolver.flushPending(), "pending latest load publishes once the first fade is idle");
        test::ok (! convolver.hasPending(), "retry queue is empty after publication");
        pumpCrossfade (convolver);
        convolver.reset();

        std::vector<float> impulse (64, 0.0f), other (64, 0.0f);
        impulse[0] = 1.0f;
        other[0] = 1.0f;
        float* render[2] { impulse.data(), other.data() };
        felitronics::test::run (convolver.process (render, 2, (int) impulse.size()));
        bool latestWon = true;
        for (std::size_t i = 0; i < latest.size(); ++i)
            latestWon = latestWon && impulse[i] == latest[i] && other[i] == latest[i];
        test::ok (latestWon, "retried operator is the latest IR and broadcasts to both channels");
    }

    // P67. The witness that the resampler did NOT run is the taps themselves: even at a ratio of exactly one it
    // band-limits at 0.95 of Nyquist and moves every tap of this noise-like IR, so a verbatim copy, to the bit,
    // cannot come out of it. The rates are LITERALS, so the documented one part per million is pinned from both
    // sides (0.9 ppm verbatim, 1.1 ppm resampled) rather than following whatever the constant says; and it is
    // RELATIVE — half a ppm at 192 kHz is 0.096 Hz, which a tolerance in hertz small enough for 48000.0000001
    // would resample.
    test::group ("P67 — rates within kRateMatchTolerance are ONE rate: the taps go in verbatim; past it they are resampled");
    {
        const auto ir = decayingIr (600);
        const auto stagedAt = [&ir] (double hostRate, double irRate, bool normalize, float* gain = nullptr)
        {
            CabConvolver convolver;
            felitronics::test::run (convolver.prepare (hostRate, 512, 2, 0.05, normalize));
            const float* banks[1] { ir.data() };
            convolver.loadIR (banks, 1, (int) ir.size(), irRate);
            if (gain != nullptr) *gain = convolver.irNormalizationGain();
            return convolver.stagedTaps().size() == 1 ? convolver.stagedTaps()[0] : std::vector<float> {};
        };
        test::approx (CabConvolver::kRateMatchTolerance, 1.0e-6, 0.0, "the tolerance is one part per million");
        test::ok (stagedAt (48000.0000001, 48000.0, false) == ir,
                  "a host at 48000.0000001 loads a 48 kHz IR verbatim, to the bit (the resampler did not run)");
        test::ok (stagedAt (48000.0, 48000.0000001, false) == ir, "and a 48000.0000001 IR on a 48 kHz host, the same");
        test::ok (stagedAt (48000.0, 48000.0432, false) == ir && stagedAt (48000.0432, 48000.0, false) == ir,
                  "0.9 ppm at 48 kHz (0.0432 Hz), both ways: verbatim");
        test::ok (stagedAt (192000.0, 192000.096, false) == ir && stagedAt (192000.096, 192000.0, false) == ir,
                  "0.5 ppm at 192 kHz (0.096 Hz), both ways: verbatim — relative, not in hertz");
        // "EVERY TAP MOVED" STOPPED BEING A WITNESS WHEN P68 LANDED, and this is the repair. These two
        // loads are normalize=false PAST the tolerance, so they now also carry the rate factor
        // (0.999998927 at 1.1 ppm) — and that scalar alone moves all 600 taps of this fixture, so the old
        // check would have passed with the resampler deleted. A gain can only SCALE; a resample changes
        // the SHAPE. So the witness is now "no single scalar maps the input onto these taps", which is
        // strictly stronger than what was here and cannot be satisfied by any future gain either.
        const auto notAScaledCopy = [&ir] (const std::vector<float>& taps)
        {
            if (taps.size() != ir.size()) return false;
            std::size_t loudest = 0;
            for (std::size_t i = 0; i < ir.size(); ++i)
                if (std::fabs (ir[i]) > std::fabs (ir[loudest])) loudest = i;
            if (! (std::fabs (ir[loudest]) > 0.0f)) return false;
            const double k = (double) taps[loudest] / (double) ir[loudest];      // the best a pure gain could do
            double worst = 0.0;
            for (std::size_t i = 0; i < ir.size(); ++i)
            {
                if (! std::isfinite (taps[i])) return false;
                worst = std::max (worst, std::fabs ((double) taps[i] - k * (double) ir[i]));
            }
            return worst > 1.0e-6 * (double) std::fabs (ir[loudest]);
        };
        test::ok (notAScaledCopy (stagedAt (48000.0, 48000.0528, false)) && notAScaledCopy (stagedAt (48000.0528, 48000.0, false)),
                  "1.1 ppm at 48 kHz (0.0528 Hz), both ways: resampled — the taps are not the input times any one number");
        float gain = 0.0f;
        const auto normalized = stagedAt (48000.0000001, 48000.0, true, &gain);
        bool scaled = normalized.size() == ir.size() && gain > 0.0f;
        for (std::size_t i = 0; scaled && i < ir.size(); ++i) scaled = normalized[i] == ir[i] * gain;
        test::ok (scaled, "normalize=true at 48000.0000001: exactly the input times the one normalization gain");
    }

    // P67. An IR rate that is not a positive finite number is UNKNOWN, and an unknown rate loads the taps as is:
    // the samples are fine, only the metadata is broken, and refusing would play silence. One rule for NaN, zero,
    // negative and both infinities — on main +inf was the one that differed (sent to the resampler, which gave
    // back nothing, and the load was dropped without a word).
    test::group ("P67 — an unknown IR rate (NaN, zero, negative, +inf, -inf) loads the taps as is, and they play");
    {
        const std::vector<float> ir { 0.5f, -0.25f, 0.125f, 0.0625f, -0.03125f, 0.015625f };   // inside the head: plays exactly
        const double inf = std::numeric_limits<double>::infinity();
        struct Unknown { double rate; const char* name; };
        for (const Unknown& unknown : { Unknown { std::numeric_limits<double>::quiet_NaN(), "NaN" }, Unknown { 0.0, "zero" },
                                        Unknown { -48000.0, "negative" }, Unknown { inf, "+inf" }, Unknown { -inf, "-inf" } })
            for (bool normalize : { false, true })
            {
                CabConvolver convolver;
                felitronics::test::run (convolver.prepare (44100.0, 128, 2, 0.05, normalize));
                const float* banks[1] { ir.data() };
                convolver.loadIR (banks, 1, (int) ir.size(), unknown.rate);
                const float gain = convolver.irNormalizationGain();
                pumpCrossfade (convolver);
                const bool flushed = convolver.flushPending();
                const auto& staged = convolver.stagedTaps();
                bool asIs = staged.size() == 1 && staged[0].size() == ir.size()
                         && (normalize ? std::isfinite (gain) && gain > 0.0f && gain != 1.0f : gain == 1.0f);
                for (std::size_t i = 0; asIs && i < ir.size(); ++i) asIs = staged[0][i] == ir[i] * gain;

                convolver.reset();
                std::vector<float> left (16, 0.0f), right (16, 0.0f);
                left[0] = right[0] = 1.0f;
                float* io[2] { left.data(), right.data() };
                felitronics::test::run (convolver.process (io, 2, 16));
                bool plays = asIs;
                for (std::size_t i = 0; plays && i < left.size(); ++i)
                    plays = left[i] == (i < ir.size() ? ir[i] * gain : 0.0f) && right[i] == left[i];
                test::ok (asIs && flushed && ! convolver.hasPending() && plays,
                          std::string ("an IR at rate ") + unknown.name + (normalize ? ", normalized" : "")
                              + ": the taps load as is, publish and play");
            }
    }

    // P67. A one-tap IR off the host rate used to resample to NO taps, and the load silently did nothing: the
    // previous cabinet kept playing. Now it stages its one tap, publishes it, and that is what plays.
    test::group ("P67 — a one-tap IR off the host rate stages one tap and is what plays");
    {
        struct Rates { double host, ir; };
        for (const Rates& rates : { Rates { 44100.0, 96000.0 }, Rates { 44100.0, 192000.0 }, Rates { 48000.0, 176400.0 },
                                    Rates { 8000.0, 384000.0 } })
        {
            CabConvolver convolver;
            felitronics::test::run (convolver.prepare (rates.host, 128, 2, 0.05));       // normalize=true: the cabinet path
            const std::vector<float> previous { 0.5f, 0.25f, -0.125f, 0.0625f };
            const float* first[1] { previous.data() };
            convolver.loadIR (first, 1, (int) previous.size(), rates.host);
            pumpCrossfade (convolver);

            const std::vector<float> tap { 0.8f };
            const float* banks[1] { tap.data() };
            convolver.loadIR (banks, 1, 1, rates.ir);
            pumpCrossfade (convolver);
            const bool flushed = convolver.flushPending();
            const auto& staged = convolver.stagedTaps();
            const bool one = staged.size() == 1 && staged[0].size() == 1 && std::isfinite (staged[0][0]) && staged[0][0] != 0.0f;
            char msg[160];
            std::snprintf (msg, sizeof msg, "one tap at %.0f Hz on a %.0f Hz host stages one tap and publishes it",
                           rates.ir, rates.host);
            test::ok (one && flushed && ! convolver.hasPending() && ! convolver.isBusy(), msg);

            convolver.reset();
            std::vector<float> left (16, 0.0f), right (16, 0.0f);
            left[0] = right[0] = 1.0f;
            float* io[2] { left.data(), right.data() };
            felitronics::test::run (convolver.process (io, 2, 16));
            bool plays = one;
            for (std::size_t i = 0; plays && i < left.size(); ++i)
                plays = left[i] == (i == 0 ? staged[0][0] : 0.0f) && right[i] == left[i];
            test::ok (plays, std::string (msg) + " — and it is what plays; the previous IR is gone");
        }
    }

    // P67. A load that stages nothing used to return AFTER overwriting the retained taps, so a load still
    // pending from mid-crossfade was retried with its OLD length over an emptied or narrowed store: a read past
    // the end of a vector, or a retry refused forever with isBusy() stuck true. Each such load must now leave
    // the pending one exactly as it was — its taps staged, its retry publishing, its IR the one that plays.
    // The malformed calls go last: before the fix they crash outright, and the earlier cases fail readably.
    test::group ("P67 — a load that stages nothing is ignored whole: the pending load survives it and plays");
    {
        const std::vector<float> first { 1.0f, 0.25f, 0.0f, 0.0f };
        const std::vector<float> pendingLeft { 0.0f, 0.25f, 0.0f, 0.0f }, pendingRight { 0.0f, 0.0f, 0.125f, 0.0f };
        enum class Nothing { zeroLength, refusedRate, negativeLength, nullArray, nullPlane, nullSecondPlane };
        // normalize=true gives the pending load a gain that is not 1, so a load that resets the gain on its way
        // out is caught as well as one that touches the taps.
        const auto survives = [&] (bool stereoPending, bool normalize, Nothing nothing, const std::string& what)
        {
            CabConvolver convolver;
            felitronics::test::run (convolver.prepare (44100.0, 128, 2, 0.05, normalize));
            const float* a[1] { first.data() };
            convolver.loadIR (a, 1, (int) first.size(), 44100.0);
            float l[64] {}, r[64] {};
            float* io[2] { l, r };
            felitronics::test::run (convolver.process (io, 2, 64));                    // the first fade is in flight
            const float* b[2] { pendingLeft.data(), pendingRight.data() };
            convolver.loadIR (b, stereoPending ? 2 : 1, 4, 44100.0);                   // rejected mid-fade -> pending
            const bool pendingBefore = convolver.hasPending();
            const float gain = convolver.irNormalizationGain();
            const float gainDb = convolver.irNormalizationGainDb();

            const std::vector<float> one { 0.8f };
            const float* c[1] { one.data() };
            const float* nullPlane[1] { nullptr };
            switch (nothing)
            {
                case Nothing::zeroLength:     convolver.loadIR (c, 1, 0, 44100.0); break;
                case Nothing::refusedRate:    convolver.loadIR (c, 1, 1, 1.0e-5); break;     // x 4.41e9: past INT_MAX
                case Nothing::negativeLength: convolver.loadIR (c, 1, -1, 44100.0); break;
                case Nothing::nullArray:      convolver.loadIR (nullptr, 1, 4, 44100.0); break;
                case Nothing::nullPlane:      convolver.loadIR (nullPlane, 1, 4, 44100.0); break;
                case Nothing::nullSecondPlane:
                {
                    const float* halfStereo[2] { one.data(), nullptr };
                    convolver.loadIR (halfStereo, 2, 1, 44100.0);
                    break;
                }
            }
            const auto scaledBy = [gain] (const std::vector<float>& v) { auto s = v; for (float& x : s) x *= gain; return s; };
            const auto& staged = convolver.stagedTaps();
            const bool retained = pendingBefore && convolver.hasPending()
                               && (normalize ? std::isfinite (gain) && gain > 0.0f && gain != 1.0f : gain == 1.0f)
                               && staged.size() == (stereoPending ? 2u : 1u) && staged[0] == scaledBy (pendingLeft)
                               && (! stereoPending || staged[1] == scaledBy (pendingRight))
                               && convolver.irNormalizationGain() == gain && convolver.irNormalizationGainDb() == gainDb;
            test::ok (retained, what + ": the pending load, its staged taps and its gain are untouched");

            pumpCrossfade (convolver);
            test::ok (convolver.flushPending() && ! convolver.hasPending(), what + ": its retry publishes");
            pumpCrossfade (convolver);
            convolver.reset();
            std::vector<float> left (8, 0.0f), right (8, 0.0f);
            left[0] = right[0] = 1.0f;
            float* render[2] { left.data(), right.data() };
            felitronics::test::run (convolver.process (render, 2, 8));
            bool plays = true;
            for (std::size_t i = 0; i < pendingLeft.size(); ++i)
                plays = plays && left[i] == pendingLeft[i] * gain
                              && right[i] == (stereoPending ? pendingRight[i] : pendingLeft[i]) * gain;
            test::ok (plays, what + ": and the pending IR is what plays");
        };
        survives (false, false, Nothing::zeroLength,  "a zero-length load over a pending mono load");
        survives (false, true,  Nothing::zeroLength,  "a zero-length load over a pending NORMALIZED load");
        survives (true,  false, Nothing::zeroLength,  "a zero-length MONO load over a pending STEREO load");
        survives (true,  true,  Nothing::refusedRate, "a load at 1e-5 Hz (its resample would be past INT_MAX) over a pending normalized stereo load");

        // And the load that used to stage nothing — one tap at 96 kHz — now stages, so it is the LATEST and wins.
        {
            CabConvolver convolver;
            felitronics::test::run (convolver.prepare (44100.0, 128, 2, 0.05, false));
            const float* a[1] { first.data() };
            convolver.loadIR (a, 1, (int) first.size(), 44100.0);
            float l[64] {}, r[64] {};
            float* io[2] { l, r };
            felitronics::test::run (convolver.process (io, 2, 64));
            const float* b[1] { pendingLeft.data() };
            convolver.loadIR (b, 1, 4, 44100.0);
            const std::vector<float> tap { 0.8f };
            const float* c[1] { tap.data() };
            convolver.loadIR (c, 1, 1, 96000.0);
            const bool latestStaged = convolver.hasPending() && convolver.stagedTaps().size() == 1
                                   && convolver.stagedTaps()[0].size() == 1;
            pumpCrossfade (convolver);
            const bool flushed = convolver.flushPending();
            pumpCrossfade (convolver);
            convolver.reset();
            std::vector<float> left (8, 0.0f), right (8, 0.0f);
            left[0] = right[0] = 1.0f;
            float* render[2] { left.data(), right.data() };
            felitronics::test::run (convolver.process (render, 2, 8));
            const bool plays = latestStaged && left[0] == convolver.stagedTaps()[0][0] && left[0] != 0.0f
                            && left[1] == 0.0f && right[0] == left[0];
            test::ok (latestStaged && flushed && plays,
                      "one tap at 96 kHz over a pending load stages as the latest, publishes, and is what plays");
        }

        survives (false, true,  Nothing::negativeLength, "a negative-length load over a pending load");
        survives (false, false, Nothing::nullArray,      "a load with no channel array over a pending load");
        survives (false, true,  Nothing::nullPlane,      "a load with a null plane over a pending load");
        survives (true,  false, Nothing::nullPlane,      "a MONO load with a null plane over a pending STEREO load");
        survives (false, true,  Nothing::nullSecondPlane, "a STEREO load whose second plane is null over a pending load");
    }

    // P69 — law 11(b) reaches the RATE too. prepare() used to answer a rate it had not been given with the
    // factory 48 kHz (`sampleRate > 0.0 ? sampleRate : 48000.0`) and then size a crossfade and an IR budget
    // from the answer, reporting success; and an infinite rate made two float-to-int conversions undefined
    // ((long long) ceil(maxIrSeconds * hostSr_), and lround(0.05 * hostSr_)). A rate outside (0, kMaxSampleRate]
    // is refused here now, and a refused prepare leaves the object unusable — which process() has to report.
    test::group ("P69 — prepare refuses a host rate it cannot honour, and a refusal disarms the convolver");
    {
        const double inf = std::numeric_limits<double>::infinity(), nan = std::numeric_limits<double>::quiet_NaN();
        struct Case { double rate; const char* what; };
        const Case bad[] { { 0.0, "zero" }, { -48000.0, "negative" }, { nan, "NaN" }, { inf, "+inf" }, { -inf, "-inf" },
                           { CabConvolver::kMaxSampleRate * 2.0, "twice the ceiling" },
                           { std::nextafter (CabConvolver::kMaxSampleRate, inf), "one ulp past the ceiling" },
                           { 1.0e300, "1e300 — finite, and what an isfinite gate would have let through" } };
        int refused = 0, disarmed = 0;
        for (const Case& c : bad)
        {
            CabConvolver cab;
            if (! cab.prepare (c.rate, 256, 2, 4.0, true)) ++refused;
            float l[8] {}, r[8] {};
            float* io[2] { l, r };
            if (! cab.process (io, 2, 8)) ++disarmed;                      // law 11: refused prepare => unusable
        }
        test::ok (refused == (int) (sizeof bad / sizeof bad[0]), "every rate outside (0, 3 MHz] is refused");
        test::ok (disarmed == (int) (sizeof bad / sizeof bad[0]), "and process() refuses afterwards — the object is not armed");

        test::ok (CabConvolver::kMaxSampleRate == 3.0e6,
                  "the ceiling is the house figure — Compressor, TruePeakLimiter, EqBand and RigPlayer spell the same one");
        CabConvolver ok1, ok2;
        // 0.001 s, not the 4 s default: what is under test is the RATE, and four seconds AT THE CEILING is a
        // twelve-million-tap schedule — a gigabyte of partitions to prove a comparison.
        test::ok (ok1.prepare (CabConvolver::kMaxSampleRate, 256, 2, 0.001, true), "the ceiling itself is honoured");
        test::ok (ok2.prepare (44100.0, 256, 2, 4.0, true), "and so is an ordinary rate");

        // A REFUSED re-prepare must not leave the PREVIOUS one usable either.
        CabConvolver rearmed;
        const bool first = rearmed.prepare (48000.0, 256, 2, 0.25, true);
        const bool second = rearmed.prepare (inf, 256, 2, 0.25, true);
        float l[8] {}, r[8] {};
        float* io[2] { l, r };
        test::ok (first && ! second && ! rearmed.process (io, 2, 8),
                  "a prepared convolver re-prepared with an infinite rate is left unusable, not left as it was");
    }

    // The same law, the same call, the other argument: `std::max (0.0, NaN)` returns its FIRST operand, so a
    // NaN duration prepared successfully with a zero IR budget — a convolver holding only the backend's
    // 128-sample head. +inf means "as much as the backend can hold", which is what every huge finite value
    // already meant, and the saturation that gives it that meaning happens in double, BEFORE the cast that
    // was undefined for it.
    //
    // THE SATURATION IS TESTED THROUGH `maxIrSamplesFor` AND NOT THROUGH `prepare`, ON PURPOSE: a budget of
    // kMaxIrSamples is a 1.34 GB partition schedule (measured with /usr/bin/time -l), and this suite also
    // runs on the wasm tier, whose heap default is 2 GB. The expected values here come from outside the
    // function — 4 s at 48 kHz is 192000 by multiplication, and the ceiling is the backend's own published
    // constant — so this is not the budget being asked to confirm itself.
    test::group ("P69 — prepare refuses a maxIrSeconds that is not a non-negative number, and saturates the rest");
    {
        const double inf = std::numeric_limits<double>::infinity(), nan = std::numeric_limits<double>::quiet_NaN();
        CabConvolver a, b, c;
        test::ok (! a.prepare (48000.0, 256, 2, nan, true), "a NaN duration is refused");
        test::ok (! b.prepare (48000.0, 256, 2, -1.0, true), "a negative duration is refused");
        test::ok (! c.prepare (48000.0, 256, 2, -inf, true), "-inf is refused");
        float l[8] {}, r[8] {};
        float* io[2] { l, r };
        test::ok (! a.process (io, 2, 8) && ! b.process (io, 2, 8) && ! c.process (io, 2, 8),
                  "and each refusal leaves its convolver unusable");

        constexpr int ceiling = felitronics::convolution::MatrixConvolverNupc<felitronics::convolution::CabConvFft>::kMaxIrSamples;
        test::ok (CabConvolver::maxIrSamplesFor (48000.0, 4.0) == 192000, "4 s at 48 kHz is 192000 samples of budget");
        test::ok (CabConvolver::maxIrSamplesFor (48000.0, 0.0) == 0, "zero seconds is zero samples");
        CabConvolver zeroSeconds;
        test::ok (zeroSeconds.prepare (48000.0, 256, 2, 0.0, false),
                  "and a zero-second budget PREPARES — the backend keeps its own head, which is what it did before");
        test::ok (CabConvolver::maxIrSamplesFor (48000.0, inf) == ceiling,
                  "+inf saturates to the backend's ceiling instead of converting an infinity to an integer");
        test::ok (CabConvolver::maxIrSamplesFor (48000.0, 1.0e300) == ceiling,
                  "and so does 1e300 — the value whose cast was undefined before the clamp moved ahead of it");
        test::ok (CabConvolver::maxIrSamplesFor (3.0e6, 1.0e12) == ceiling, "at the rate ceiling too");
        test::ok (CabConvolver::maxIrSamplesFor (48000.0, 1.0 / 48000.0) == 1, "and one sample is one sample");
        test::ok (CabConvolver::maxIrSamplesFor (48000.0, 1.25 / 48000.0) == 2,
                  "a part-sample duration rounds UP — every whole-sample case alone cannot tell ceil from floor");
    }

    // A refused prepare has to leave the WHOLE object unusable, not just its flag. The pending-retry
    // geometry outlived a refusal, and a retry that survives an unusable object can never be published:
    // flushPending() returns on !prepared_, so isBusy() answered true for the rest of the object's life
    // and the host's reload poll spun on it forever. Four calls reach it.
    test::group ("P69 — a refused prepare clears the pending retry it can no longer publish");
    {
        std::vector<float> first { 1.0f, 0.25f, 0.0f, 0.0f }, second { 0.0f, 1.0f, 0.0f, 0.0f };
        CabConvolver cab;
        felitronics::test::run (cab.prepare (48000.0, 128, 2, 0.05, false));
        const float* a[1] { first.data() };
        const float* b[1] { second.data() };
        cab.loadIR (a, 1, (int) first.size(), 48000.0);
        float left[64] {}, right[64] {};
        float* io[2] { left, right };
        felitronics::test::run (cab.process (io, 2, 64));                  // begins the crossfade
        cab.loadIR (b, 1, (int) second.size(), 48000.0);                   // rejected mid-fade, retained
        const bool pendingBefore = cab.hasPending();
        const bool refused = ! cab.prepare (std::numeric_limits<double>::infinity(), 128, 2, 0.05, false);
        test::ok (pendingBefore && refused && ! cab.hasPending(),
                  "the retry a refused prepare can never publish does not survive it");
    }

    // `maxIrSamplesFor` is PUBLIC, so it answers for every double, not only for the ones prepare() lets
    // through. A comment is not a precondition: `std::min (NaN, ceiling)` keeps the NaN, and converting
    // that to int is the same undefined conversion this whole group exists to have closed.
    test::group ("P69 — the budget helper is total: no argument converts a non-number to an int");
    {
        const double inf2 = std::numeric_limits<double>::infinity(), nan2 = std::numeric_limits<double>::quiet_NaN();
        test::ok (CabConvolver::maxIrSamplesFor (nan2, 1.0) == 0, "a NaN rate gives no budget");
        test::ok (CabConvolver::maxIrSamplesFor (48000.0, nan2) == 0, "nor does a NaN duration");
        test::ok (CabConvolver::maxIrSamplesFor (-inf2, 1.0) == 0 && CabConvolver::maxIrSamplesFor (48000.0, -inf2) == 0,
                  "nor -inf on either side");
        test::ok (CabConvolver::maxIrSamplesFor (-48000.0, 1.0) == 0 && CabConvolver::maxIrSamplesFor (48000.0, -1.0) == 0,
                  "nor a negative rate or duration");
        test::ok (CabConvolver::maxIrSamplesFor (0.0, 4.0) == 0, "and a zero rate is zero samples, not a conversion");
    }

    // P69 — the analysis window and the transform that carries it. `normalizationGain` takes the first
    // second of the IR (`cap`) but sizes its FFT at `N`, which stops doubling at 1<<21; above a megasample
    // of window the two part company and the copy ran off the end of the buffer. Measured on the base
    // commit under AddressSanitizer as a 12 MB heap-buffer-overflow WRITE at a 4 MHz host; this is the
    // smallest input that reaches it — one sample past 2^21, at a rate one sample past 2^21 — so the
    // overwrite there is a single float, which only a sanitizer or a hardened allocator will notice. The
    // IR budget is kept tiny on purpose: it is the NORMALIZATION path under test, not the partition build.
    test::group ("P69 — the reference-gain analysis window cannot outrun its own transform");
    {
        constexpr int len = (1 << 21) + 1;
        CabConvolver cab;
        const bool armed = cab.prepare ((double) len, 64, 1, 0.001, true);
        std::vector<float> ir ((std::size_t) len, 0.0f);
        // The impulse is NOT at index 0, and that is the point: an amplitude of 2 anywhere inside the window
        // gives |H(f)| = 2 at every frequency and so a reference gain of exactly 0.5, while a window
        // truncated to nothing — or to one sample — sees only silence, falls through the near-silent floor
        // and returns unity. A pure impulse at index 0 could not tell those apart.
        ir[1000] = 2.0f;
        ir[(std::size_t) len - 1] = 0.25f;                                 // one past the window's last sample
        const float* planes[1] { ir.data() };
        cab.loadIR (planes, 1, len, (double) len);                         // same rate: no resample, straight to normalization
        const float g = cab.irNormalizationGain();
        // The VALUE, not merely finiteness: an analysis that copied nothing would return unity and pass a
        // "positive and finite" check while hiding exactly the defect this group exists for.
        test::ok (armed && std::fabs ((double) g - 0.5) < 1.0e-3,
                  "a 2 097 153-sample IR at a 2 097 153 Hz host normalizes to 0.5 without writing past the transform");
    }

    // P68 — AN IR'S AUTHORED LEVEL IS A CONVOLUTION GAIN, NOT A TAP AMPLITUDE, and a convolution gain is a
    // sum over taps: proportional to how many of them fit into a second. `resampleIr` preserves the
    // amplitude (unity DC per output tap), so before this the same file convolved outSr/inSr times as loud
    // — a 48 kHz spring measured +6.02 dB on a 96 kHz host and -0.74 dB on a 44.1 one, for one Mix knob.
    // Only normalize=false could show it: the reference-unity path measures the FINAL taps and swallows
    // the factor whole.
    //
    // THE TOLERANCE IS THE KERNEL'S, not a taste. What is left after the factor is taken out is the
    // 64-tap Kaiser (beta = 8) windowed sinc's own passband ripple: its stopband is
    // 8/0.1102 + 8.7 = 81.3 dB, so delta = 8.6e-5 and the ripple is 20*log10(1 + delta) = 0.00075 dB
    // one-sided. Its column sums (the quantity the factor actually inverts) were measured across 64
    // phases of every standard rate pair and stray at most 1.6e-4 dB from the ratio; the worst reading
    // over this grid is 5.65e-4 dB. So 0.01 dB is 13x the kernel's ripple and 18x the worst reading,
    // while the weakest mutation below moves a cell by 0.736 dB — 74 tolerances away. (The Kaiser
    // figure is an engineering formula, not a proof, which is why the measured column sums stand beside
    // it rather than behind it.)
    test::group ("P68 — normalize=false: the authored convolution gain survives the host's clock");
    {
        // TWO BANKS ON HALF THE CELLS, and not for symmetry: with every cell mono, "apply the factor only
        // when nch == 1" is a change to shipped sound that passes this whole file. The second bank carries
        // its own signal (the fixture's mirror image), so a mutation that scales one plane and not the
        // other cannot hide behind a broadcast either.
        struct Cell { double irSr, hostSr; int channels; };
        const Cell cells[] {
            { 48000.0,  44100.0, 1 }, { 48000.0,  88200.0, 2 }, { 48000.0, 96000.0, 1 }, { 48000.0, 192000.0, 2 },
            { 96000.0,  48000.0, 2 }, { 96000.0,  44100.0, 1 },                 // ratio < 1: the OTHER side
            { 44100.0,  48000.0, 2 }, { 44100.0, 192000.0, 1 },                 // a source off the family grid
        };
        const double probes[] { 0.0, 100.0, 1000.0, 5000.0, 10000.0 };          // all well inside the passband
        constexpr double tolerance = 0.01;                                       // dB — derived above

        // WHAT THIS GROUP DOES NOT CERTIFY, said out loud so nobody reads it as more than it is. It reads
        // MAGNITUDE, so a change that moved the taps in time (a reversal, a shift) is invisible here — the
        // resampler's own shift-invariance and edge-tap groups own that. It reads up to 10 kHz, so the
        // transition band, the images and the anti-aliasing are invisible too — the resampler's "a
        // cabinet's top octave keeps its own level" group owns those. What is left is exactly what the
        // rate factor is: a flat level, and whether it is the IR's own at every host rate.

        double worst = 0.0;
        double worstIr = 0.0, worstHost = 0.0, worstHz = 0.0;
        bool gainsExact = true, allFinite = true;
        for (const Cell& cell : cells)
        {
            const auto source = rateFixture (cell.irSr, 0.1, 2 * felitronics::convolution::IrResampleConfig{}.halfTaps);
            std::vector<float> second (source.size());
            for (std::size_t i = 0; i < source.size(); ++i) second[i] = -0.5f * source[source.size() - 1 - i];
            CabConvolver convolver;
            felitronics::test::run (convolver.prepare (cell.hostSr, 512, cell.channels, 4.0, /*normalize*/ false));
            const float* banks[2] { source.data(), second.data() };
            convolver.loadIR (banks, cell.channels, (int) source.size(), cell.irSr);

            const float applied = convolver.irNormalizationGain();
            gainsExact = gainsExact && applied == (float) (cell.irSr / cell.hostSr)
                                    // ...and the dB reading is the SAME number, not a floored stand-in
                                    && std::fabs ((double) convolver.irNormalizationGainDb()
                                                  - 20.0 * std::log10 ((double) applied)) < 1.0e-4;
            for (int c = 0; c < cell.channels; ++c)
            {
                const auto& staged = convolver.stagedTaps()[(std::size_t) c];
                const auto& plane  = c == 0 ? source : second;
                for (const double hz : probes)
                {
                    const double moved = 20.0 * std::log10 (magnitudeAt (staged, cell.hostSr, hz)
                                                            / magnitudeAt (plane, cell.irSr, hz));
                    allFinite = allFinite && std::isfinite (moved);
                    if (std::fabs (moved) > worst)
                    { worst = std::fabs (moved); worstIr = cell.irSr; worstHost = cell.hostSr; worstHz = hz; }
                }
            }
        }
        std::printf ("    %d rate pairs x %d frequencies: worst |level change| = %.2e dB at %.0f -> %.0f Hz host, "
                     "%.0f Hz (tolerance %.2f; an uncompensated load moves the weakest cell by 0.736)\n",
                     (int) (sizeof (cells) / sizeof (cells[0])), (int) (sizeof (probes) / sizeof (probes[0])),
                     worst, worstIr, worstHost, worstHz, tolerance);
        test::ok (allFinite && worst <= tolerance,
                  "the response at every physical frequency is the IR's own, at every host rate");
        test::ok (gainsExact, "irNormalizationGain() reports exactly irSr/hostSr on every resampled reverb load");
    }

    // The other half of the same contract — and the half that pins the DOMAIN of the change. A load that
    // is not resampled must still hand the taps over untouched, to the bit, with a gain of exactly one:
    // the rate factor is 1 there by construction, and "exactly 1.0f" is what says it was never applied
    // as an approximation of 1.
    test::group ("P68 — a load that is not resampled is still byte-verbatim with a gain of exactly one");
    {
        const std::vector<float> ir { 0.75f, -0.25f, 0.125f, 0.0625f, -0.5f };
        // 48000.0432 is 0.9 ppm off — inside kRateMatchTolerance, so no resample runs — and it is there
        // because 48000.0000001 cannot do its job alone: (float)(48000/48000.0000001) IS 1.0f, so float
        // rounding would hand the assertion its answer and a missing `resample ?` gate would pass. At 0.9
        // ppm the factor would be 0.99999911f, and only the gate can make it one.
        const double asIs[] { 48000.0, 48000.0000001, 48000.0432,               // the same rate, three ways
                              std::numeric_limits<double>::quiet_NaN(), 0.0, -48000.0,
                              std::numeric_limits<double>::infinity(),
                              -std::numeric_limits<double>::infinity() };      // ...and every unknown one
        bool verbatim = true, unity = true;
        for (const double rate : asIs)
        {
            CabConvolver convolver;
            felitronics::test::run (convolver.prepare (48000.0, 128, 2, 0.05, /*normalize*/ false));
            const float* banks[1] { ir.data() };
            convolver.loadIR (banks, 1, (int) ir.size(), rate);
            verbatim = verbatim && convolver.stagedTaps()[0] == ir;
            unity    = unity && convolver.irNormalizationGain() == 1.0f
                             && convolver.irNormalizationGainDb() == 0.0f;
        }
        test::ok (verbatim, "the host's own rate, one within the tolerance, and five unusable ones: taps byte-verbatim");
        test::ok (unity, "and the applied gain is exactly 1.0f (0.0 dB) on every one of them");
    }


    // P68 — THE dB DIAGNOSTIC'S FLOOR, which this change made reachable. `irNormalizationGainDb()` floors
    // its argument so a zero cannot read as -inf; the floor was 1e-6, which sat BELOW every gain the old
    // code could apply (the normalization clamps at -30 dB = 0.0316) and ABOVE what the rate factor can
    // reach. A file claiming 0.024 Hz is not a real file, but it is a rate the loader accepts — P67 settled
    // that a finite positive rate is a KNOWN rate — and it resamples one tap into two million, applying
    // 5.0e-7. The reading said -120.0000 dB for a gain of -126.0206.
    // THIS IS THE MOST EXPENSIVE CASE IN THE FILE AND ITS COST IS MEASURED, NOT GUESSED: 1.6 s of the
    // binary's 2.1 s on macOS/arm64 clang, and the whole file runs 3.5 s in the CHECKED wasm row
    // (SAFE_HEAP + assertions, which instruments every load and store across 1.28e8 window-tap
    // iterations) — the row worth naming, because it is the one a local release build never exercises.
    // There is no cheaper way to stand on that floor: a gain under 1e-6 needs an output of inLen/g taps
    // and inLen is already 1, so the million is arithmetic, not a choice.
    test::group ("P68 — an extreme accepted rate reports its true dB, not the anti-infinity floor");
    {
        std::vector<float> one { 1.0f };
        CabConvolver convolver;
        felitronics::test::run (convolver.prepare (48000.0, 128, 1, 0.05, /*normalize*/ false));
        const float* banks[1] { one.data() };
        convolver.loadIR (banks, 1, 1, 0.024);
        const double applied = (double) convolver.irNormalizationGain();
        const double reported = (double) convolver.irNormalizationGainDb();
        std::printf ("    0.024 Hz -> 48 kHz: %d taps, gain %.6e, %.4f dB (the floor would say -120.0000)\n",
                     (int) convolver.stagedTaps()[0].size(), applied, reported);
        test::ok (convolver.stagedTaps()[0].size() == 2000000, "one tap at 0.024 Hz stages two million");
        test::approx (applied, 0.024 / 48000.0, 1.0e-13, "the applied gain is the rate ratio, 5e-7");
        test::approx (reported, 20.0 * std::log10 (0.024 / 48000.0), 1.0e-3,
                      "and its dB reading is -126.0206, the gain that was actually applied");
    }

    // P68 — A GAIN THAT IS NOT A NUMBER MUST NOT READ AS A NUMBER. `std::max(a, b)` is `(a < b) ? b : a`
    // and every comparison against a NaN is false, so `max(floor, NaN)` returns the FLOOR: the dB accessor
    // answered a plausible -120.0000 for a gain that is a NaN, on the NORMALIZED path, on both sides of
    // this change (the floor moved, so the plausible lie moved with it — -120 before, -400 after, which is
    // how a "nothing changed here" claim gets falsified by an input nobody fixtured). A NaN tap reaches
    // this from a float32 WAV, and the linear accessor has always told the truth about it.
    test::group ("P68 — a NaN gain reports a NaN, not the anti-infinity floor dressed as a level");
    {
        std::vector<float> ir (256, 0.25f);
        ir[0] = std::numeric_limits<float>::quiet_NaN();
        CabConvolver convolver;
        felitronics::test::run (convolver.prepare (48000.0, 128, 1, 0.05, /*normalize*/ true));
        const float* banks[1] { ir.data() };
        convolver.loadIR (banks, 1, (int) ir.size(), 48000.0);
        test::ok (std::isnan (convolver.irNormalizationGain()), "a NaN tap makes the measured gain a NaN");
        test::ok (std::isnan (convolver.irNormalizationGainDb()),
                  "and the dB reading is a NaN too — not -120, and not -400 either");
    }

    // P68 — AND THE NORMALIZED PATH DID NOT MOVE, which is a claim and not an assumption. It cannot be
    // shown by the reference-unity check above: that one runs at the host's own rate, where the rate
    // factor is 1 and any pre-scaling would be invisible. It cannot be shown by re-normalizing either —
    // reference-unity is HOMOGENEOUS, so scaling the taps by k before measuring them multiplies the
    // measured gain by 1/k and the product comes out the same to within the last bits. (A mutation that
    // applies the rate factor on the normalize=true path as well passed every other check in this file.)
    // What DOES see it: recompute the gain here, from the resampled taps, with the loader's own resampler
    // but nothing else of the loader's, and ask whether the number the loader reports is the gain of
    // THOSE taps — not of taps somebody scaled first. At 192 kHz a pre-scaled measurement reads 12.04 dB
    // away from the true one in whichever direction the scaling went, which is 120 times this group's
    // tolerance; where that lands relative to the +-30 dB clamp depends on the IR, and for this quiet
    // fixture it lands nowhere near it. The clamp is not what catches the mutation — the 12 dB is.
    test::group ("P68 — a resampled normalize=true load reports the gain of the taps it actually plays");
    {
        // A SMOOTH, QUIET fixture, and both adjectives are load-bearing. Smooth (a one-pole decay, so
        // |H(f)| has no peaks) because the two grids here are different — the loader's zero-padded FFT
        // against a linear DTFT — and a peaky spectrum makes them disagree by more than the thing being
        // measured. Quiet (x 0.05) because the loudness this is checking has to stay INSIDE the +-30 dB
        // clamp at every rate: at 192 kHz the density adds 12 dB to the reading, and a fixture at the
        // fixture's natural level clamps there, which would compare two guards instead of two gains.
        std::vector<float> source ((std::size_t) (64 + 1440), 0.0f);
        for (int i = 0; i < 1440; ++i)
            source[(std::size_t) (64 + i)] = (float) (0.05 * std::exp (-(double) i / 48.0));
        bool matched = true, unclamped = true;
        for (const double hostSr : { 44100.0, 88200.0, 96000.0, 192000.0 })
        {
            const auto resampled = felitronics::convolution::resampleIr (source, 48000.0, hostSr);
            const double expected = referenceNormalizationGainDb (resampled, hostSr,
                                        (int) std::lround (4096.0 * hostSr / 48000.0));   // constant in HERTZ

            CabConvolver convolver;
            felitronics::test::run (convolver.prepare (hostSr, 512, 1, 4.0, /*normalize*/ true));
            const float* banks[1] { source.data() };
            convolver.loadIR (banks, 1, (int) source.size(), 48000.0);
            const double reported = (double) convolver.irNormalizationGainDb();

            std::printf ("    %7.0f Hz host: reported %+8.4f dB, recomputed %+8.4f dB, difference %+.4f\n",
                         hostSr, reported, expected, reported - expected);
            // 0.1 dB, and the headroom here is 1.9x, not the 13-74x of the first group — say so rather
            // than let that paragraph's numbers be read as covering this one. The two grids sit a
            // SYSTEMATIC 0.043-0.052 dB apart once the spacing is constant in hertz (it was 0.02 -> 0.18
            // and rate-dependent before that), which is grid geometry, not drift: gcc/glibc/x86-64
            // reproduces all four figures to the digit. What the budget has to clear is the mutation —
            // the rate factor applied here as well — and that moves the reading by 0.736 dB at the
            // weakest rate, 14 budgets away.
            matched = matched && std::fabs (reported - expected) <= 0.1;
            unclamped = unclamped && std::fabs (reported) < 29.99;             // the guard must not be in play
        }
        test::ok (matched && unclamped,
                  "at 44.1 / 88.2 / 96 / 192 kHz the reported gain is the reference-unity gain of the resampled taps");
    }

    return test::report();
}

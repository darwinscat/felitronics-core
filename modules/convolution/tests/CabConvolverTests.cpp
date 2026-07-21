// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// Golden coverage for the product-level cab wrapper: reference-unity gain math, direct-convolution
// parity, mono broadcast / true-stereo publication, verbatim reverb loading, and latest-wins retry.

#include <felitronics_test.h>
#include <felitronics/convolution/CabConvolver.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdio>
#include <limits>
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
        convolver.process (io, channels, 512);    // 6144 samples > the fixed 50 ms fade at 48 kHz
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
        convolver.process (io, 2, count);
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
double referenceNormalizationGainDb (const std::vector<float>& ir, double sampleRate)
{
    constexpr int bins = 4096;
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
            convolver.prepare (48000.0, 512, 2, 0.05);
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
        convolver.prepare (48000.0, 512, 2, 0.05);
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
            convolver.prepare (48000.0, 512, 2, 0.05);
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
        convolver.prepare (48000.0, 128, 2, (double) capSamples / 48000.0, false);
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
        convolver.prepare (48000.0, 256, 2, 0.05);
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
        convolver.process (io, 2, 256);
        bool lrDiag = true;
        for (std::size_t i = 0; i < leftIr.size(); ++i)
            lrDiag = lrDiag && std::fabs (left[i] - staged[0][i]) < 2.0e-4f
                            && std::fabs (right[i] - staged[1][i]) < 2.0e-4f;
        test::ok (lrDiag, "rendered stereo impulse uses the distinct L/R diagonal banks");
    }

    test::group ("normalize=false is a verbatim reverb path");
    {
        std::vector<float> ir { 0.75f, -0.25f, 0.125f, 0.0625f };
        CabConvolver convolver;
        convolver.prepare (48000.0, 128, 2, 0.05, false);
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
        convolver.prepare (48000.0, 128, 2, 0.05, false);
        const float* a[1] { first.data() };
        const float* b[1] { second.data() };
        const float* c[1] { latest.data() };
        convolver.loadIR (a, 1, (int) first.size(), 48000.0);

        float left[64] {}, right[64] {};
        float* io[2] { left, right };
        convolver.process (io, 2, 64);       // begins the first publication's 50 ms crossfade
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
        convolver.process (render, 2, (int) impulse.size());
        bool latestWon = true;
        for (std::size_t i = 0; i < latest.size(); ++i)
            latestWon = latestWon && impulse[i] == latest[i] && other[i] == latest[i];
        test::ok (latestWon, "retried operator is the latest IR and broadcasts to both channels");
    }

    return test::report();
}

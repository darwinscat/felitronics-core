// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// A channel that stops being split and is split again. Behind one `c < nc` sits the entire crossover
// tree for that column — an Svf in every Linkwitz-Riley crossover and in every allpass compensator —
// plus its per-band alignment delay lines. None of it decays while the channel is away.
//
// Measured before this suite, stereo -> mono -> stereo with DIGITAL SILENCE on the return: 1.55e-01,
// which is -16.2 dBFS. The band processors are not the subject here: each carries its own ledger, or
// needs none, and this suite uses a pass-through band so only the splitter and the aligners are in play.

#include <felitronics_test.h>
#include <felitronics/multiband/MultibandProcessor.h>

#include <cmath>
#include <vector>

using namespace felitronics;
using felitronics::test::ok;
using felitronics::test::group;

// A band that does nothing, so anything heard on the return came from the splitter or the aligners.
struct PassBand
{
    bool prepare() noexcept { return true; }
    void reset() noexcept {}
    void process (float* const*, int, int) noexcept {}
    int  latencySamples() const noexcept { return 0; }
};

// A band whose latency is set per INSTANCE, because the aligners hold `maxLatency - ownLatency` samples:
// with every band reporting the same number every line is zero-length and holds nothing. A first version
// of this suite used a zero-latency band, a second used one constant latency for all four, and both were
// blind to the alignment delays for the same reason.
struct LateBand
{
    int lat = 0;
    bool prepare (int l) noexcept { lat = l; return true; }
    void reset() noexcept {}
    void process (float* const*, int, int) noexcept {}
    int  latencySamples() const noexcept { return lat; }
};

static double peakOf (const std::vector<float>& v)
{
    double m = 0.0; for (float x : v) m = std::fmax (m, (double) std::fabs (x)); return m;
}

static bool bitEqual (const std::vector<float>& a, const std::vector<float>& b)
{
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i) if (! (a[i] == b[i])) return false;
    return true;
}

static void fillNoise (std::vector<float>& v, unsigned seed)
{
    unsigned s = seed;
    for (auto& x : v) { s = s * 1664525u + 1013904223u; x = (float) ((int) (s >> 8) % 2000 - 1000) * 0.0005f; }
}

int main()
{
    std::printf ("felitronics::multiband — execution-gate state\n");

    const int N = 128;

    group ("channel gate: a returning channel emits nothing from digital silence");
    {
        multiband::MultibandProcessor<PassBand, 4> m;
        ok (m.prepare (48000.0, N, 2, 0, [] (PassBand& b) { return b.prepare(); }),
            "precondition: the processor prepared for two channels");

        std::vector<float> L ((std::size_t) N, 0.0f), R ((std::size_t) N, 0.0f);
        float* io[2] { L.data(), R.data() };

        double charged = 0.0;
        for (int k = 0; k < 40; ++k)
        {
            std::fill (L.begin(), L.end(), 0.0f); fillNoise (R, 7u + (unsigned) k);
            felitronics::test::run (m.process (io, 2, N));
            charged = std::fmax (charged, peakOf (R));
        }
        ok (charged > 0.1, "precondition: the right channel really was being split");

        for (int k = 0; k < 40; ++k) { std::fill (L.begin(), L.end(), 0.0f); felitronics::test::run (m.process (io, 1, N)); }

        double worst = 0.0;
        for (int k = 0; k < 8; ++k)
        {
            std::fill (L.begin(), L.end(), 0.0f); std::fill (R.begin(), R.end(), 0.0f);
            felitronics::test::run (m.process (io, 2, N));
            worst = std::fmax (worst, peakOf (R));
        }
        ok (worst == 0.0, "silence in, exact zero out on the returned channel (was 1.55e-01 = -16.2 dBFS)");
    }

    group ("isolation: the channel that never left keeps its history bit-exact");
    {
        multiband::MultibandProcessor<PassBand, 4> dut, ref;
        ok (dut.prepare (48000.0, N, 2, 0, [] (PassBand& b) { return b.prepare(); })
            && ref.prepare (48000.0, N, 2, 0, [] (PassBand& b) { return b.prepare(); }),
            "precondition: both processors prepared");

        std::vector<float> d0 ((std::size_t) N), d1 ((std::size_t) N), r0 ((std::size_t) N), r1 ((std::size_t) N);
        float* dio[2] { d0.data(), d1.data() };
        float* rio[2] { r0.data(), r1.data() };

        bool equal = true; double energy = 0.0;
        for (int k = 0; k < 120; ++k)
        {
            fillNoise (d0, 51u + (unsigned) k);
            std::fill (d1.begin(), d1.end(), 0.0f);
            r0 = d0; r1 = d1;
            felitronics::test::run (dut.process (dio, (k >= 40 && k < 80) ? 1 : 2, N));
            felitronics::test::run (ref.process (rio, 2, N));
            equal = equal && bitEqual (d0, r0);
            energy += peakOf (d0);
        }
        ok (energy > 1.0, "precondition: the surviving channel carried signal");
        ok (equal, "channel 0 is bit-equal to the run where nothing ever left");
    }

    group ("the per-band alignment delays are cleared too");
    {
        // A band that declares latency makes the processor align the others through DelayLines, one per
        // band and channel. Those hold audio, and they are per channel, so a returning channel replays
        // them exactly like the crossover columns. Only a band with non-zero latency puts them in play.
        multiband::MultibandProcessor<LateBand, 4> m;
        int next = 0;
        ok (m.prepare (48000.0, N, 2, 256, [&next] (LateBand& b) { return b.prepare (64 * next++); }),
            "precondition: the processor prepared with room to align 256 samples");
        ok (m.latencySamples() >= 192, "precondition: the bands report DIFFERENT latencies, so the aligners hold audio");

        std::vector<float> L ((std::size_t) N, 0.0f), R ((std::size_t) N, 0.0f);
        float* io[2] { L.data(), R.data() };

        double charged = 0.0;
        for (int k = 0; k < 40; ++k)
        {
            std::fill (L.begin(), L.end(), 0.0f); fillNoise (R, 23u + (unsigned) k);
            felitronics::test::run (m.process (io, 2, N));
            charged = std::fmax (charged, peakOf (R));
        }
        ok (charged > 0.1, "precondition: the right channel really was being split and aligned");

        for (int k = 0; k < 40; ++k) { std::fill (L.begin(), L.end(), 0.0f); felitronics::test::run (m.process (io, 1, N)); }

        double worst = 0.0;
        for (int k = 0; k < 8; ++k)
        {
            std::fill (L.begin(), L.end(), 0.0f); std::fill (R.begin(), R.end(), 0.0f);
            felitronics::test::run (m.process (io, 2, N));
            worst = std::fmax (worst, peakOf (R));
        }
        ok (worst == 0.0, "silence in, exact zero out — the alignment lines came back empty too");
    }

    return felitronics::test::report();
}

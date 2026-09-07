// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// A convolver channel that stops being asked for and is asked for again. It holds a whole frame, an FDL
// of past input spectra and a pending tail — none of which decays while the channel is absent. Measured
// before this suite, stereo -> mono -> stereo with DIGITAL SILENCE on the return: 7.44e-01, which is
// -2.6 dBFS, the loudest leak of this class found anywhere in the repository.
//
// The engine's own block position (`phase_`, `fdlPos_`, the crossfade) is deliberately untouched: it is
// shared by every channel and is supposed to keep running for the ones that stayed.

#include <felitronics_test.h>
#include <felitronics/convolution/ConvolutionEngine.h>

#include <cmath>
#include <vector>

using namespace felitronics;
using felitronics::test::ok;
using felitronics::test::group;

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

static std::vector<float> decayIr (int len)
{
    std::vector<float> ir ((std::size_t) len, 0.0f);
    ir[0] = 1.0f;
    for (int i = 1; i < len; ++i) ir[(std::size_t) i] = (float) (0.5 * std::exp (-i / 60.0));
    return ir;
}

int main()
{
    std::printf ("felitronics::convolution — execution-gate state\n");

    const int N = 128;
    const std::vector<float> ir = decayIr (512);

    group ("channel gate: a returning channel emits nothing from digital silence");
    {
        convolution::ConvolutionEngine<> e;
        ok (e.prepare (128, 1024, 64, 2), "precondition: the engine prepared for two channels");
        ok (e.setIr (ir.data(), (int) ir.size()), "precondition: an impulse response was accepted");

        std::vector<float> L ((std::size_t) N, 0.0f), R ((std::size_t) N, 0.0f);
        float* out[2] { L.data(), R.data() };
        const float* in[2] { L.data(), R.data() };

        double charged = 0.0;
        for (int k = 0; k < 40; ++k)
        {
            std::fill (L.begin(), L.end(), 0.0f); fillNoise (R, 5u + (unsigned) k);
            felitronics::test::run (e.process (in, out, 2, N));
            charged = std::fmax (charged, peakOf (R));
        }
        ok (charged > 0.1, "precondition: the right channel really was convolving");

        for (int k = 0; k < 40; ++k) { std::fill (L.begin(), L.end(), 0.0f); felitronics::test::run (e.process (in, out, 1, N)); }

        double worst = 0.0;
        for (int k = 0; k < 8; ++k)
        {
            std::fill (L.begin(), L.end(), 0.0f); std::fill (R.begin(), R.end(), 0.0f);
            felitronics::test::run (e.process (in, out, 2, N));
            worst = std::fmax (worst, peakOf (R));
        }
        ok (worst == 0.0, "silence in, exact zero out on the returned channel (was 7.44e-01 = -2.6 dBFS)");
    }

    group ("isolation: the channel that never left keeps its history bit-exact");
    {
        convolution::ConvolutionEngine<> dut, ref;
        ok (dut.prepare (128, 1024, 64, 2) && ref.prepare (128, 1024, 64, 2), "precondition: both engines prepared");
        ok (dut.setIr (ir.data(), (int) ir.size()) && ref.setIr (ir.data(), (int) ir.size()), "precondition: both loaded");

        std::vector<float> d0 ((std::size_t) N), d1 ((std::size_t) N), r0 ((std::size_t) N), r1 ((std::size_t) N);
        float* dout[2] { d0.data(), d1.data() };
        const float* din[2] { d0.data(), d1.data() };
        float* rout[2] { r0.data(), r1.data() };
        const float* rin[2] { r0.data(), r1.data() };

        bool equal = true; double energy = 0.0;
        for (int k = 0; k < 120; ++k)
        {
            fillNoise (d0, 41u + (unsigned) k);
            std::fill (d1.begin(), d1.end(), 0.0f);
            r0 = d0; r1 = d1;
            felitronics::test::run (dut.process (din, dout, (k >= 40 && k < 80) ? 1 : 2, N));
            felitronics::test::run (ref.process (rin, rout, 2, N));
            equal = equal && bitEqual (d0, r0);
            energy += peakOf (d0);
        }
        ok (energy > 1.0, "precondition: the surviving channel carried signal");
        ok (equal, "channel 0 is bit-equal to the run where nothing ever left");
    }

    return felitronics::test::report();
}

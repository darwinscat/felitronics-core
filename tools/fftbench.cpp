// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// Acceptance benchmark for the #1-core-debt fix. Measures the partitioned convolver's %RT (CPU time / audio
// time) at ~LinearPhaseEq-Maximum FIR length (131072 taps), across host block sizes, for the OLD partition
// policy (P >= host block, which made cost EXPLODE with the block) vs the NEW fixed P=128 (block-INDEPENDENT),
// on the scalar reference vs the pffft SIMD backend, in BOTH the common LRDiag (2-conv stereo) and the
// worst-case Full (4-conv 2x2 matrix) routings. Empirical proof of CLAUDE.md's acceptance: "comparable perf
// vs a SIMD reference at block 2048 and larger".
//
// Dev tool (top-level builds only). NOT real-time (mt19937 + printf); build Release, run manually. Numbers are
// steady-state (the timed window starts after the cold FDL prime + a cache warm-up). Run with
// -DFELITRONICS_WITH_PFFFT=ON for the SIMD column.

#include <felitronics/convolution/MatrixConvolver.h>
#include <felitronics/core/Fft.h>
#include <felitronics/core/FlushToZero.h>
#if defined(FELITRONICS_WITH_PFFFT)
 #include <felitronics/fftpffft/PffftRealFft.h>
#endif

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

using namespace felitronics;

// Steady-state %RT of a stereo MatrixConvolver<Backend> at partition P and host block `block`, for `nBanks`
// IRs of length irLen routed by topology `topoSel` (0 = LRDiag, 1 = Full).
template <class Backend>
static double benchRT (int topoSel, int nBanks, int irLen, int P, int block, double fs, double warmSec, double measSec)
{
    using MC = convolution::MatrixConvolver<Backend>;
    MC mc;
    const int warmXfade = std::max (2 * P, (int) std::lround (0.02 * fs));
    if (! mc.prepare (P, irLen, warmXfade, 2)) return -1.0;

    std::mt19937 rng (12345);
    std::uniform_real_distribution<float> u (-1.0f, 1.0f);
    std::vector<std::vector<float>> irs ((std::size_t) nBanks, std::vector<float> ((std::size_t) irLen));
    for (auto& ir : irs)
        for (int i = 0; i < irLen; ++i)
            ir[(std::size_t) i] = 0.02f * std::exp (-3.0f * (float) i / (float) irLen) * u (rng);   // decaying, healthy (no denormals)
    std::vector<const float*> banks ((std::size_t) nBanks);
    for (int b = 0; b < nBanks; ++b) banks[(std::size_t) b] = irs[(std::size_t) b].data();

    using Topo = typename MC::Topology;
    mc.setOperator (topoSel == 0 ? Topo::LRDiag : Topo::Full, banks.data(), nBanks, irLen);

    std::vector<float> L ((std::size_t) block), R ((std::size_t) block);
    for (int i = 0; i < block; ++i) { L[(std::size_t) i] = 0.3f * u (rng); R[(std::size_t) i] = 0.3f * u (rng); }

    core::ScopedFlushToZero ftz;   // a real audio host runs with FTZ/DAZ — measure under the same FP mode
    auto once = [&] { float* io[2] { L.data(), R.data() }; mc.process (io, io, 2, block); };

    for (int i = 0, w = (int) (warmSec * fs / block); i < w; ++i) once();   // pass the cold prime + warm caches

    const int measBlocks = std::max (1, (int) (measSec * fs / block));
    const auto t0 = std::chrono::steady_clock::now();
    for (int i = 0; i < measBlocks; ++i) once();
    const auto t1 = std::chrono::steady_clock::now();

    const double cpu   = std::chrono::duration<double> (t1 - t0).count();
    const double audio = (double) measBlocks * (double) block / fs;
    return 100.0 * cpu / audio;
}

int main()
{
    const double fs = 48000.0;
    const int irLen = 131072;   // ~ LinearPhaseEq Maximum (N=131072; the EQ passes N+1 taps — the +1 is perf-negligible)
    const int blocks[] = { 256, 512, 1024, 2048, 4096, 8192 };
    const struct { int sel, nBanks; const char* name; } topos[] = {
        { 0, 2, "LRDiag  (stereo, 2 convolutions — the common case)" },
        { 1, 4, "Full    (2x2 matrix, 4 convolutions — worst-case routing)" },
    };

    std::printf ("felitronics-core #1-debt acceptance bench — MatrixConvolver, %d-tap IR @ %.0f Hz\n", irLen, fs);
    std::printf ("%%RT = CPU / audio time; lower is better (100%% = the whole real-time budget). Steady-state, FTZ on.\n");
#if defined(FELITRONICS_WITH_PFFFT)
    std::printf ("pffft SIMD width = %d (4 = SSE/NEON active)\n", fftpffft::PffftRealFft::simdWidth());
#else
    std::printf ("(configure -DFELITRONICS_WITH_PFFFT=ON for the SIMD column)\n");
#endif

    for (const auto& tp : topos)
    {
        std::printf ("\n== %s ==\n", tp.name);
#if defined(FELITRONICS_WITH_PFFFT)
        std::printf ("  host block |  OLD scalar (P>=block)  |  NEW scalar (P=128)  |  NEW pffft (P=128)\n");
        std::printf ("  ----------- ------------------------- ---------------------- --------------------\n");
#else
        std::printf ("  host block |  OLD scalar (P>=block)  |  NEW scalar (P=128)\n");
        std::printf ("  ----------- ------------------------- ----------------------\n");
#endif
        for (int b : blocks)
        {
            int oldP = 64; while (oldP < b) oldP <<= 1;   // the pre-PR4 policy: P = pow2 >= host block
            const double oldS = benchRT<core::fft::ScalarRadix2Real> (tp.sel, tp.nBanks, irLen, oldP, b, fs, 3.0, 2.0);
            const double newS = benchRT<core::fft::ScalarRadix2Real> (tp.sel, tp.nBanks, irLen, 128,  b, fs, 3.0, 2.0);
#if defined(FELITRONICS_WITH_PFFFT)
            const double newP = benchRT<fftpffft::PffftRealFft>       (tp.sel, tp.nBanks, irLen, 128,  b, fs, 3.0, 2.0);
            std::printf ("  %9d  |  %7.2f%%  (P=%5d)  |  %7.2f%%           |  %7.2f%%\n", b, oldS, oldP, newS, newP);
#else
            std::printf ("  %9d  |  %7.2f%%  (P=%5d)  |  %7.2f%%\n", b, oldS, oldP, newS);
#endif
        }
    }
    std::printf ("\nOLD scalar grows with the host block (the shipped bottleneck); NEW is block-independent;\n");
    std::printf ("NEW pffft is lowest — and the fix holds even in Full (worst-case 4-conv) routing.\n");
    return 0;
}

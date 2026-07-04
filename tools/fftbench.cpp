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
#if defined(FELITRONICS_BENCH_JUCE)
 #include "JuceRealFft.h"   // bench-only JUCE/vDSP FFT reference (FetchContent'd, never in core)
 #include <juce_events/juce_events.h>   // MessageManager — juce::dsp::Convolution's loader needs it in a headless tool
 #include <thread>
#endif

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
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

    std::vector<float> inL ((std::size_t) block), inR ((std::size_t) block);
    for (int i = 0; i < block; ++i) { inL[(std::size_t) i] = 0.3f * u (rng); inR[(std::size_t) i] = 0.3f * u (rng); }
    std::vector<float> L ((std::size_t) block), R ((std::size_t) block);

    core::ScopedFlushToZero ftz;   // a real audio host runs with FTZ/DAZ — measure under the same FP mode
    // Refill a FRESH bounded input each block. process() is in-place, so without this the output feeds back and
    // degenerates (grows/decays), which skews an input-sensitive convolver like JUCE's. The memcpy is identical
    // overhead for every backend and negligible vs the convolution.
    auto once = [&]
    {
        std::memcpy (L.data(), inL.data(), (std::size_t) block * sizeof (float));
        std::memcpy (R.data(), inR.data(), (std::size_t) block * sizeof (float));
        float* io[2] { L.data(), R.data() }; mc.process (io, io, 2, block);
    };

    for (int i = 0, w = (int) (warmSec * fs / block); i < w; ++i) once();   // pass the cold prime + warm caches

    const int measBlocks = std::max (1, (int) (measSec * fs / block));
    const auto t0 = std::chrono::steady_clock::now();
    for (int i = 0; i < measBlocks; ++i) once();
    const auto t1 = std::chrono::steady_clock::now();

    const double cpu   = std::chrono::duration<double> (t1 - t0).count();
    const double audio = (double) measBlocks * (double) block / fs;
    return 100.0 * cpu / audio;
}

#if defined(FELITRONICS_BENCH_JUCE)
// The FAIR head-to-head: JUCE's OWN partitioned convolver (juce::dsp::Convolution) at the same IR/block,
// vs our MatrixConvolver<Pffft>. loadImpulseResponse is async (background thread) → poll getCurrentIRSize()
// until the IR is live before warming/measuring. Returns -1 if the IR never loaded.
static double benchJuceConvolution (int irLen, int block, double fs, double warmSec, double measSec, int& latencyOut)
{
    juce::dsp::Convolution conv;
    juce::dsp::ProcessSpec spec { fs, (juce::uint32) block, 2 };
    conv.prepare (spec);

    juce::AudioBuffer<float> ir (2, irLen);
    std::mt19937 rng (12345);
    std::uniform_real_distribution<float> u (-1.0f, 1.0f);
    for (int ch = 0; ch < 2; ++ch)
    { float* p = ir.getWritePointer (ch); for (int i = 0; i < irLen; ++i) p[i] = 0.02f * std::exp (-3.0f * (float) i / (float) irLen) * u (rng); }
    conv.loadImpulseResponse (std::move (ir), fs, juce::dsp::Convolution::Stereo::yes,
                              juce::dsp::Convolution::Trim::no, juce::dsp::Convolution::Normalise::no);

    juce::AudioBuffer<float> inBuf (2, block), outBuf (2, block);
    for (int ch = 0; ch < 2; ++ch) { float* p = inBuf.getWritePointer (ch); for (int i = 0; i < block; ++i) p[i] = 0.3f * u (rng); }
    outBuf.clear();
    // NON-replacing: read the constant inBuf, write outBuf — the input never degenerates into feedback.
    auto once = [&] { juce::dsp::AudioBlock<float> ib (inBuf), ob (outBuf);
                      juce::dsp::ProcessContextNonReplacing<float> ctx (ib, ob); conv.process (ctx); };

    core::ScopedFlushToZero ftz;
    // getCurrentIRSize() starts at 1 and only reaches irLen once the background loader has published — pump the
    // message loop + process until it does (or give up → -1).
    for (int g = 0; conv.getCurrentIRSize() < irLen && g < 3000; ++g)
    {
        once();                                                            // the background loader publishes; process() swaps it in
        std::this_thread::sleep_for (std::chrono::milliseconds (1));
    }
    if (conv.getCurrentIRSize() < irLen) return -1.0;                      // IR never fully loaded (headless)
    latencyOut = (int) conv.getLatency();

    for (int i = 0, w = (int) (warmSec * fs / block); i < w; ++i) once();
    const int measBlocks = std::max (1, (int) (measSec * fs / block));
    const auto t0 = std::chrono::steady_clock::now();
    for (int i = 0; i < measBlocks; ++i) once();
    const auto t1 = std::chrono::steady_clock::now();
    return 100.0 * std::chrono::duration<double> (t1 - t0).count() / ((double) measBlocks * (double) block / fs);
}

// DECISIVE sanity check: does juce::dsp::Convolution actually convolve the FULL irLen-tap IR (so its low %RT
// is legit), or is it silently truncating/approximating? Feed a unit impulse → the output IS the impulse
// response, which must equal the loaded IR (latency 0, Normalise::no). Compare sample-by-sample + by energy.
static void juceCorrectnessCheck (int irLen, double fs)
{
    const int block = 512;
    juce::dsp::Convolution conv;
    conv.prepare ({ fs, (juce::uint32) block, 1u });

    juce::AudioBuffer<float> ir (1, irLen);
    std::mt19937 rng (777);
    std::uniform_real_distribution<float> u (-1.0f, 1.0f);
    std::vector<float> irRef ((std::size_t) irLen);
    { float* p = ir.getWritePointer (0);
      for (int i = 0; i < irLen; ++i) { irRef[(std::size_t) i] = 0.02f * std::exp (-3.0f * (float) i / (float) irLen) * u (rng); p[i] = irRef[(std::size_t) i]; } }
    conv.loadImpulseResponse (std::move (ir), fs, juce::dsp::Convolution::Stereo::no,
                              juce::dsp::Convolution::Trim::no, juce::dsp::Convolution::Normalise::no);

    juce::AudioBuffer<float> buf (1, block);
    auto proc = [&] { juce::dsp::AudioBlock<float> ab (buf); juce::dsp::ProcessContextReplacing<float> ctx (ab); conv.process (ctx); };
    for (int g = 0; conv.getCurrentIRSize() < irLen && g < 20000; ++g) { buf.clear(); proc(); std::this_thread::sleep_for (std::chrono::milliseconds (1)); }
    if (conv.getCurrentIRSize() < irLen) { std::printf ("  [correctness] IR did not load — skip\n"); return; }
    for (int i = 0; i < 40; ++i) { buf.clear(); proc(); }   // drain JUCE's ~50 ms IR CROSSFADE (getCurrentIRSize flips before it ends — codex)
    conv.reset();                                            // clear the FDL so the impulse response is measured clean

    std::vector<float> out; out.reserve ((std::size_t) (irLen + block));
    bool firstBlock = true;
    while ((int) out.size() < irLen)
    {
        buf.clear();
        if (firstBlock) { buf.setSample (0, 0, 1.0f); firstBlock = false; }   // unit impulse at sample 0
        proc();
        for (int i = 0; i < block && (int) out.size() < irLen; ++i) out.push_back (buf.getSample (0, i));
    }
    double maxErr = 0.0, irPeak = 0.0, eIr = 0.0, eOut = 0.0;
    for (int i = 0; i < irLen; ++i)
    {
        maxErr = std::max (maxErr, (double) std::fabs (out[(std::size_t) i] - irRef[(std::size_t) i]));
        irPeak = std::max (irPeak, (double) std::fabs (irRef[(std::size_t) i]));
        eIr  += (double) irRef[(std::size_t) i] * irRef[(std::size_t) i];
        eOut += (double) out[(std::size_t) i] * out[(std::size_t) i];
    }
    std::printf ("  [correctness] JUCE impulse-response vs the %d-tap IR: maxErr=%.2e (IR peak %.2e, rel %.2e), energy out/IR=%.4f\n",
                 irLen, maxErr, irPeak, maxErr / (irPeak + 1e-30), eOut / (eIr + 1e-30));
    std::printf ("  [correctness] => %s\n", maxErr < 1e-3 * irPeak
                 ? "JUCE reproduces the FULL IR — its low %RT is a real full-length convolution"
                 : "MISMATCH — JUCE is NOT convolving the full IR (its low %RT would be bogus)");
}
#endif

int main()
{
#if defined(FELITRONICS_BENCH_JUCE)
    juce::MessageManager::getInstance();   // this thread becomes JUCE's message thread → the Convolution loader runs
#endif
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

#if defined(FELITRONICS_BENCH_JUCE)
    std::printf ("\n== JUCE correctness probe (does juce::dsp::Convolution really convolve all %d taps?) ==\n", irLen);
    juceCorrectnessCheck (irLen, fs);
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

#if defined(FELITRONICS_BENCH_JUCE)
    // THE FAIR HEAD-TO-HEAD: JUCE's own convolver vs our MatrixConvolver<Pffft>, same IR, LRDiag/stereo.
    std::printf ("\n== JUCE head-to-head (same 131072-tap stereo IR) ==\n");
    std::printf ("  NB juce::dsp::Convolution reports ZERO latency too — it is UNIFORM-partitioned (partition =\n");
    std::printf ("  host block), so a bigger block => far fewer partitions => cheaper. Ours is a fixed P=128\n");
    std::printf ("  head+tail (block-INDEPENDENT) — the flat ~2%% is its 1023 partitions + an O(P) scalar head.\n");
    std::printf ("  host block |  juce::dsp::Convolution (its latency)  |  our MatrixConvolver<pffft> (0 latency)\n");
    std::printf ("  ----------- --------------------------------------- ----------------------------------------\n");
    for (int b : blocks)
    {
        int lat = 0;
        const double j = benchJuceConvolution (irLen, b, fs, 3.0, 2.0, lat);
        const double p = benchRT<fftpffft::PffftRealFft> (0, 2, irLen, 128, b, fs, 3.0, 2.0);
        if (j < 0.0) std::printf ("  %9d  |  (IR load failed)                     |  %7.2f%%\n", b, p);
        else         std::printf ("  %9d  |  %6.2f%%  (lat %6d smp = %6.1f ms)  |  %7.2f%%  (0.0 ms)\n",
                                  b, j, lat, 1000.0 * lat / fs, p);
    }
    // Same-engine FFT-backend reference (JUCE's FFT/vDSP in OUR engine; scalar MAC + 2N layout → a handicap,
    // NOT juce::dsp::Convolution — shown only to isolate the FFT).
    std::printf ("\n== JUCE-FFT-in-our-engine (P=128, LRDiag; scalar MAC — proxy only) ==\n");
    for (int b : blocks)
        std::printf ("  %9d  |  %7.2f%%\n", b, benchRT<bench::JuceRealFft> (0, 2, irLen, 128, b, fs, 3.0, 2.0));
#endif
    return 0;
}

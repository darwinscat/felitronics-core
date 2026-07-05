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
#include <felitronics/convolution/MatrixConvolverNupc.h>
#include <felitronics/convolution/NonUniformConvolver.h>
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
#include <cstdlib>   // getenv (FCORE_SWEEP_ONLY)
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

// Steady-state %RT of the Phase-1 NonUniformConvolver<Backend> in the common LRDiag case = TWO mono instances
// (L, R). Head P0 + capped-doubling to Bmax, irLen taps, host block `block`. Returns the MEAN %RT; writes the
// WORST single-block %RT to maxBlockPctOut (the mean hides the periodic spike when every stage's FFT coincides
// — the per-buffer cost a real-time thread must survive). block-INDEPENDENT + TRUE sample-zero-latency.
template <class Backend>
static double benchNU (int irLen, int P0, int Bmax, int block, double fs, double warmSec, double measSec, double& maxBlockPctOut)
{
    using NUc = convolution::NonUniformConvolver<Backend>;
    NUc l, r;
    if (! l.prepare (P0, Bmax, irLen) || ! r.prepare (P0, Bmax, irLen)) { maxBlockPctOut = -1.0; return -1.0; }

    std::mt19937 rng (12345);
    std::uniform_real_distribution<float> u (-1.0f, 1.0f);
    std::vector<float> irL ((std::size_t) irLen), irR ((std::size_t) irLen);
    for (int i = 0; i < irLen; ++i)
    {
        irL[(std::size_t) i] = 0.02f * std::exp (-3.0f * (float) i / (float) irLen) * u (rng);
        irR[(std::size_t) i] = 0.02f * std::exp (-3.0f * (float) i / (float) irLen) * u (rng);
    }
    l.setIr (irL.data(), irLen);
    r.setIr (irR.data(), irLen);

    std::vector<float> inL ((std::size_t) block), inR ((std::size_t) block), L ((std::size_t) block), R ((std::size_t) block);
    for (int i = 0; i < block; ++i) { inL[(std::size_t) i] = 0.3f * u (rng); inR[(std::size_t) i] = 0.3f * u (rng); }

    core::ScopedFlushToZero ftz;
    auto once = [&]
    {
        std::memcpy (L.data(), inL.data(), (std::size_t) block * sizeof (float));
        std::memcpy (R.data(), inR.data(), (std::size_t) block * sizeof (float));
        l.process (L.data(), L.data(), block);
        r.process (R.data(), R.data(), block);
    };

    for (int i = 0, w = (int) (warmSec * fs / block); i < w; ++i) once();   // pass the cold prime + warm caches

    const int measBlocks = std::max (1, (int) (measSec * fs / block));
    const double audioPerBlock = (double) block / fs;

    // Pass 1 — clean whole-loop timing for the MEAN (no per-block clock overhead).
    const auto t0 = std::chrono::steady_clock::now();
    for (int i = 0; i < measBlocks; ++i) once();
    const auto t1 = std::chrono::steady_clock::now();
    const double mean = 100.0 * std::chrono::duration<double> (t1 - t0).count() / ((double) measBlocks * audioPerBlock);

    // Pass 2 — per-block timing for the WORST buffer (captures the coincident-FFT spike).
    double maxBlk = 0.0;
    for (int i = 0; i < measBlocks; ++i)
    {
        const auto b0 = std::chrono::steady_clock::now();
        once();
        const auto b1 = std::chrono::steady_clock::now();
        maxBlk = std::max (maxBlk, std::chrono::duration<double> (b1 - b0).count() / audioPerBlock);
    }
    maxBlockPctOut = 100.0 * maxBlk;
    return mean;
}

// The SHIPPING matrix convolver: a single stereo MatrixConvolverNupc<Backend> routed by `topo` (0 LRDiag,
// 1 MSDiag, 2 Full). This is what lineareq now convolves on. Mean %RT + worst single-block %RT. block-INDEPENDENT.
template <class Backend>
static double benchMatrixNupc (int topoSel, int nBanks, int irLen, int block, double fs, double warmSec, double measSec, double& maxBlockPctOut)
{
    using MC = convolution::MatrixConvolverNupc<Backend>;
    MC mc;
    const int warmXfade = std::max (256, (int) std::lround (0.02 * fs));
    if (! mc.prepare (128, irLen, warmXfade, 2)) { maxBlockPctOut = -1.0; return -1.0; }

    std::mt19937 rng (12345);
    std::uniform_real_distribution<float> u (-1.0f, 1.0f);
    std::vector<std::vector<float>> irs ((std::size_t) nBanks, std::vector<float> ((std::size_t) irLen));
    for (auto& ir : irs) for (int i = 0; i < irLen; ++i) ir[(std::size_t) i] = 0.02f * std::exp (-3.0f * (float) i / (float) irLen) * u (rng);
    std::vector<const float*> banks ((std::size_t) nBanks);
    for (int b = 0; b < nBanks; ++b) banks[(std::size_t) b] = irs[(std::size_t) b].data();
    using Topo = typename MC::Topology;
    mc.setOperator (topoSel == 0 ? Topo::LRDiag : topoSel == 1 ? Topo::MSDiag : Topo::Full, banks.data(), nBanks, irLen);

    std::vector<float> inL ((std::size_t) block), inR ((std::size_t) block), L ((std::size_t) block), R ((std::size_t) block);
    for (int i = 0; i < block; ++i) { inL[(std::size_t) i] = 0.3f * u (rng); inR[(std::size_t) i] = 0.3f * u (rng); }

    core::ScopedFlushToZero ftz;
    auto once = [&]
    {
        std::memcpy (L.data(), inL.data(), (std::size_t) block * sizeof (float));
        std::memcpy (R.data(), inR.data(), (std::size_t) block * sizeof (float));
        float* io[2] { L.data(), R.data() }; mc.process (io, io, 2, block);
    };

    for (int i = 0, w = (int) (warmSec * fs / block); i < w; ++i) once();   // pass the cold prime + warm caches

    const int measBlocks = std::max (1, (int) (measSec * fs / block));
    const double audioPerBlock = (double) block / fs;
    const auto t0 = std::chrono::steady_clock::now();
    for (int i = 0; i < measBlocks; ++i) once();
    const auto t1 = std::chrono::steady_clock::now();
    const double mean = 100.0 * std::chrono::duration<double> (t1 - t0).count() / ((double) measBlocks * audioPerBlock);

    double maxBlk = 0.0;
    for (int i = 0; i < measBlocks; ++i) { const auto b0 = std::chrono::steady_clock::now(); once(); const auto b1 = std::chrono::steady_clock::now(); maxBlk = std::max (maxBlk, std::chrono::duration<double> (b1 - b0).count() / audioPerBlock); }
    maxBlockPctOut = 100.0 * maxBlk;
    return mean;
}

#if defined(FELITRONICS_BENCH_JUCE)
// The FAIR head-to-head: JUCE's OWN convolver (juce::dsp::Convolution). prepMaxBlock is what the plugin declares
// at prepare() (ProcessSpec::maximumBlockSize — this is what JUCE sizes its partitioning + reports its latency
// from); hostBlock is what the host actually delivers to process() (hostBlock <= prepMaxBlock). They are usually
// equal (the head-to-head passes them equal), but DECOUPLING them exposes the real-world case where a host
// declares a large max yet delivers smaller/variable buffers. Async IR load → poll getCurrentIRSize(). -1 if it
// never loaded.
static double benchJuceConvolution (int irLen, int prepMaxBlock, int hostBlock, double fs, double warmSec, double measSec, int& latencyOut)
{
    juce::dsp::Convolution conv;
    juce::dsp::ProcessSpec spec { fs, (juce::uint32) prepMaxBlock, 2 };
    conv.prepare (spec);

    juce::AudioBuffer<float> ir (2, irLen);
    std::mt19937 rng (12345);
    std::uniform_real_distribution<float> u (-1.0f, 1.0f);
    for (int ch = 0; ch < 2; ++ch)
    { float* p = ir.getWritePointer (ch); for (int i = 0; i < irLen; ++i) p[i] = 0.02f * std::exp (-3.0f * (float) i / (float) irLen) * u (rng); }
    conv.loadImpulseResponse (std::move (ir), fs, juce::dsp::Convolution::Stereo::yes,
                              juce::dsp::Convolution::Trim::no, juce::dsp::Convolution::Normalise::no);

    juce::AudioBuffer<float> inBuf (2, hostBlock), outBuf (2, hostBlock);
    for (int ch = 0; ch < 2; ++ch) { float* p = inBuf.getWritePointer (ch); for (int i = 0; i < hostBlock; ++i) p[i] = 0.3f * u (rng); }
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

    for (int i = 0, w = (int) (warmSec * fs / hostBlock); i < w; ++i) once();
    const int measBlocks = std::max (1, (int) (measSec * fs / hostBlock));
    const auto t0 = std::chrono::steady_clock::now();
    for (int i = 0; i < measBlocks; ++i) once();
    const auto t1 = std::chrono::steady_clock::now();
    return 100.0 * std::chrono::duration<double> (t1 - t0).count() / ((double) measBlocks * (double) hostBlock / fs);
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
    const bool sweepOnly = std::getenv ("FCORE_SWEEP_ONLY") != nullptr;   // FCORE_SWEEP_ONLY=1 → run only the DAW-buffer sweep
    const int blocks[] = { 256, 512, 1024, 2048, 4096, 8192 };
    // NUPC + JUCE head-to-head also cover the SMALL blocks that low-latency / live rigs (guitar amp, monitoring)
    // actually run — that is exactly NUPC's win regime, so it must be measured, not extrapolated.
    const int nuBlocks[] = { 64, 128, 256, 512, 1024, 2048, 4096, 8192 };
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
  if (! sweepOnly) {
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

    //==========================================================================================================
    // Phase-1 NonUniformConvolver (Gardner): head 128 + capped-doubling FFT tail (default Bmax=4096), LRDiag =
    // 2 mono instances. TRUE sample-zero-latency, block-INDEPENDENT. The uniform MatrixConvolver above uses
    // 1023 partitions (P=128); NUPC uses ~5 doubling + 1 uniform stage → far cheaper spectral MAC. The
    // max/block column exposes the periodic coincident-FFT spike (a per-buffer RT cost the mean hides).
    {
        const int P0 = 128, Bmax = 4096;
        std::printf ("\n== NonUniformConvolver — LRDiag (head %d + capped-doubling to Bmax=%d) ==\n", P0, Bmax);
#if defined(FELITRONICS_WITH_PFFFT)
        std::printf ("  host block |  NUPC scalar (mean)  |  NUPC pffft (mean)  |  NUPC pffft (max/block)\n");
        std::printf ("  ----------- ---------------------- --------------------- ------------------------\n");
        for (int b : nuBlocks)
        {
            double sMax = 0.0, pMax = 0.0;
            const double s = benchNU<core::fft::ScalarRadix2Real> (irLen, P0, Bmax, b, fs, 3.0, 2.0, sMax);
            const double p = benchNU<fftpffft::PffftRealFft>       (irLen, P0, Bmax, b, fs, 3.0, 2.0, pMax);
            std::printf ("  %9d  |  %7.2f%%           |  %7.2f%%          |  %7.2f%%\n", b, s, p, pMax);
        }
        // Bmax sweep — find the platform optimum (Fable: octave-doubling is ~26% suboptimal; the near-field/tail
        // split is tunable). At the SMALL block (64 — our realtime target) the coincident-FFT spike matters most,
        // so watch max/block, not just mean: a big B_max lands a big FFT in a tiny callback.
        for (int probeBlk : { 64, 512 })
        {
            std::printf ("\n  Bmax sweep (NUPC pffft, host block %d):  ", probeBlk);
            for (int Bm : { 1024, 2048, 4096, 8192 })
            {
                double mx = 0.0; const double m = benchNU<fftpffft::PffftRealFft> (irLen, P0, Bm, probeBlk, fs, 3.0, 2.0, mx);
                std::printf ("Bmax=%d: %.2f%% (max %.2f%%)   ", Bm, m, mx);
            }
        }
        std::printf ("\n");
#else
        std::printf ("  host block |  NUPC scalar (mean)  |  NUPC scalar (max/block)\n");
        std::printf ("  ----------- ---------------------- -------------------------\n");
        for (int b : nuBlocks)
        {
            double sMax = 0.0;
            const double s = benchNU<core::fft::ScalarRadix2Real> (irLen, P0, Bmax, b, fs, 3.0, 2.0, sMax);
            std::printf ("  %9d  |  %7.2f%%           |  %7.2f%%\n", b, s, sMax);
        }
        std::printf ("  (configure -DFELITRONICS_WITH_PFFFT=ON for the SIMD NUPC column)\n");
#endif
        std::printf ("  → block-INDEPENDENT mean (~0.9%%) + zero latency at EVERY host block from ONE prepare.\n");
        std::printf ("    JUCE is zero-latency too, but its cost TRACKS the actual block (cheap only at big blocks);\n");
        std::printf ("    NUPC's edge is a FLAT cost and the small-block regime — see the head-to-head below.\n");
    }

    //==========================================================================================================
    // The SHIPPING matrix convolver — MatrixConvolverNupc, the class lineareq now convolves on. A SINGLE stereo
    // instance on the raw-L/R history, routed per topology (LRDiag / MSDiag / Full). block-INDEPENDENT, 0 latency.
    {
        std::printf ("\n== MatrixConvolverNupc — the shipping matrix convolver (head 128, Bmax=2048, %d-tap IR) ==\n", irLen);
#if defined(FELITRONICS_WITH_PFFFT)
        std::printf ("  host block |  LRDiag pffft mean(max) |  MSDiag pffft mean(max) |  Full pffft mean(max)\n");
        std::printf ("  ----------- ------------------------- ------------------------- -----------------------\n");
        for (int b : nuBlocks)
        {
            double m0 = 0.0, m1 = 0.0, m2 = 0.0;
            const double lr = benchMatrixNupc<fftpffft::PffftRealFft> (0, 2, irLen, b, fs, 3.0, 2.0, m0);
            const double ms = benchMatrixNupc<fftpffft::PffftRealFft> (1, 2, irLen, b, fs, 3.0, 2.0, m1);
            const double fu = benchMatrixNupc<fftpffft::PffftRealFft> (2, 4, irLen, b, fs, 3.0, 2.0, m2);
            std::printf ("  %9d  |  %6.2f%% (%5.2f%%)       |  %6.2f%% (%5.2f%%)       |  %6.2f%% (%5.2f%%)\n", b, lr, m0, ms, m1, fu, m2);
        }
#else
        std::printf ("  (configure -DFELITRONICS_WITH_PFFFT=ON for the SIMD matrix-NUPC columns)\n");
#endif
        std::printf ("  → ONE stereo instance, every topology FLAT + true sample-zero-latency — this is what the EQ ships on.\n");
    }
  }   // end if (! sweepOnly) — the non-JUCE %RT tables

#if defined(FELITRONICS_BENCH_JUCE)
  if (! sweepOnly) {
    // THE HEAD-TO-HEAD nose-to-nose with JUCE, on the SHIPPING matrix convolver + the old fixed-P=128 one.
    std::printf ("\n== HEAD-TO-HEAD: juce::dsp::Convolution vs our convolvers (131072-tap stereo LRDiag) ==\n");
    std::printf ("  All three report ZERO latency. JUCE's cost falls with the block (bigger partitions); ours are\n");
    std::printf ("  FLAT (block-independent). MatrixConvolverNupc is the shipping matrix convolver.\n");
    std::printf ("  host block |  juce::dsp::Convolution  |  MatrixConvolver<pffft> (old) |  MatrixConvolverNupc<pffft>\n");
    std::printf ("  ----------- -------------------------- ------------------------------- ---------------------------\n");
    for (int b : nuBlocks)
    {
        int lat = 0;
        const double j  = benchJuceConvolution (irLen, b, b, fs, 3.0, 2.0, lat);            // oracle-tuned (JUCE's best case)
        const double mc = benchRT<fftpffft::PffftRealFft> (0, 2, irLen, 128, b, fs, 3.0, 2.0);   // old fixed-P=128 matrix convolver
        double nMax = 0.0;
        const double nu = benchMatrixNupc<fftpffft::PffftRealFft> (0, 2, irLen, b, fs, 3.0, 2.0, nMax);
        if (j < 0.0) std::printf ("  %9d  |  (IR load failed)        |  %7.2f%%                     |  %6.2f%% (max %5.2f%%)\n", b, mc, nu, nMax);
        else std::printf ("  %9d  |  %6.2f%% (%6.1f ms lat)  |  %7.2f%%                     |  %6.2f%% (max %5.2f%%)\n",
                          b, j, 1000.0 * lat / fs, mc, nu, nMax);
    }
    std::printf ("  → MatrixConvolverNupc: FLAT ~0.9%%, 3-9x cheaper than JUCE at the 64-128 blocks live rigs run,\n");
    std::printf ("    ~2.5x cheaper than the old fixed-P=128 matrix convolver, true sample-zero-latency, block-independent.\n");
    // MEASURED (not assumed): the head-to-head oracle-tunes JUCE's maxBlock = the actual host block (its best
    // case). A plugin declares ONE maximumBlockSize at prepare(); the host may then deliver SMALLER/variable
    // buffers. Prepare JUCE for maxBlock=8192, feed it a smaller block: JUCE STAYS zero-latency, but its cost
    // tracks the ACTUAL block — fed 256 it pays ~2.8% (worse than its own 256-tuned 2.06%), while NUPC is flat
    // ~0.9% from one prepare. (Correcting an earlier assumption that JUCE would gain latency here — it does not.)
    std::printf ("\n== JUCE prepared maxBlock=8192, then fed a SMALLER host block (real-world variable buffers) ==\n");
    std::printf ("  host block |  juce::dsp::Convolution(maxBlock=8192)  |  NUPC<pffft> (one prepare, 0 latency)\n");
    std::printf ("  ----------- --------------------------------------- ----------------------------------------\n");
    for (int hb : { 256, 1024, 8192 })
    {
        int lat = 0;
        const double j = benchJuceConvolution (irLen, 8192, hb, fs, 3.0, 2.0, lat);
        double nuMax = 0.0;
        const double nu = benchNU<fftpffft::PffftRealFft> (irLen, 128, 4096, hb, fs, 3.0, 2.0, nuMax);
        if (j < 0.0) std::printf ("  %9d  |  (IR load failed)                     |  %6.2f%%  (0.0 ms)\n", hb, nu);
        else         std::printf ("  %9d  |  %6.2f%%  (lat %6d smp = %6.1f ms)  |  %6.2f%%  (max %5.2f%%, 0.0 ms)\n",
                                  hb, j, lat, 1000.0 * lat / fs, nu, nuMax);
    }
    std::printf ("  → JUCE's cost follows the ACTUAL block (small block = expensive, always); NUPC is flat. A\n");
    std::printf ("    host running small/variable buffers is cheaper + steadier on NUPC, both at zero latency.\n");
  }   // end if (! sweepOnly) — the JUCE head-to-head tables

    // REAL DAW buffer sizes — including the NON-power-of-two ones a host actually offers (96 / 160 / 192 / 992…).
    // NUPC is FLAT at every one (block-independent handles ANY n); JUCE's cost tracks the buffer. This is the
    // block-independence claim measured at the sizes users pick, not just powers of two.
    std::printf ("\n== Real DAW buffer sweep (EVERY size a host offers) — juce::dsp::Convolution vs MatrixConvolverNupc, LRDiag ==\n");
    std::printf ("  buffer,ms,juce_pct,nupc_pct,nupc_maxpct   (CSV: every 32-sample DAW buffer; * suffix = non-pow2)\n");
    std::vector<int> dawBufs { 16, 32, 64 };
    for (int b = 96; b <= 2048; b += 32) dawBufs.push_back (b);   // the real DAW ladder: every 32 samples
    // warm MUST outlast NUPC's cold-prime crossfade (coldXfade_ = max_s(C_s*B_s) ≈ 2.5 s for Bmax=2048 over a
    // 131072-tap IR) or the measure window still double-convolves the warm+cold slots and INFLATES the mean by
    // ~0.35% — matches the 3.0 s head-to-head only once the crossfade has fully settled.
    for (int b : dawBufs)
    {
        int lat = 0;
        const double j  = benchJuceConvolution (irLen, b, b, fs, 3.0, 1.0, lat);
        double nMax = 0.0;
        const double nu = benchMatrixNupc<fftpffft::PffftRealFft> (0, 2, irLen, b, fs, 3.0, 1.0, nMax);
        std::printf ("  %d%s,%.2f,%.2f,%.2f,%.2f\n", b, (b & (b - 1)) ? "*" : "", 1000.0 * b / fs, j, nu, nMax);
    }
    std::printf ("  → NUPC dead-flat at EVERY buffer; JUCE a sawtooth that exceeds real-time at <= 32 samples.\n");

    if (! sweepOnly)
    {
        // Same-engine FFT-backend reference (JUCE's FFT/vDSP in OUR engine; scalar MAC + 2N layout → a handicap,
        // NOT juce::dsp::Convolution — shown only to isolate the FFT).
        std::printf ("\n== JUCE-FFT-in-our-engine (P=128, LRDiag; scalar MAC — proxy only) ==\n");
        for (int b : blocks)
            std::printf ("  %9d  |  %7.2f%%\n", b, benchRT<bench::JuceRealFft> (0, 2, irLen, 128, b, fs, 3.0, 2.0));
    }
#endif
    return 0;
}

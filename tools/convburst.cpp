// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// convburst — the WORST process() call of MatrixConvolverNupc, not its mean. A partitioned convolver's mean can
// be small while one call in every B_max samples does a whole stage's work; a host does not care about the mean,
// it misses a deadline on the worst call. This measures, for a Full (true-stereo) IR of a given length under a
// 20 s cap at 48 kHz, the worst call out of a swap and the worst call from a publish until the swap's fade has
// ended (the start of a swap is where the incoming operator is primed), at host blocks of 64, 128 and 256, each
// as a share of the block's real-time budget. The pffft rows appear with -DFELITRONICS_WITH_PFFFT=ON.
//
//   fcore_convburst [seconds of audio per row, default 20]
//
// Timing, not a test: CI gates the same property structurally (the work per call, counted through an FFT
// backend that counts). Run it on the machines whose budgets matter; a laptop's frequency scaling and an
// unpinned hybrid core (P or E) both widen the spread, so read the worst over a long run, not one number.
#include <felitronics/convolution/MatrixConvolverNupc.h>
#if defined(FELITRONICS_WITH_PFFFT)
 #include <felitronics/fftpffft/PffftRealFft.h>
#endif
#include <felitronics/core/FlushToZero.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

namespace fc = felitronics::convolution;

namespace
{
struct Lcg { unsigned long long s; float next() { s = s * 6364136223846793005ULL + 1442695040888963407ULL; return (float) ((s >> 40) & 0xffff) / 32768.0f - 1.0f; } };

struct Worst { double meanUs = 0.0, p999Us = 0.0, maxUs = 0.0; long calls = 0; };

// The maximum, and the 99.9th percentile beside it: on a laptop an isolated maximum can be the machine (an
// interrupt, a frequency step) rather than the code, and the two together say which.
Worst worstOf (std::vector<double> t)
{
    Worst w; w.calls = (long) t.size();
    if (t.empty()) return w;
    double sum = 0.0;
    for (double v : t) sum += v;
    w.meanUs = sum / (double) t.size();
    std::sort (t.begin(), t.end());
    w.p999Us = t[(std::size_t) std::min ((double) t.size() - 1.0, std::floor (0.999 * (double) t.size()))];
    w.maxUs = t.back();
    return w;
}

template <class Fft>
void row (const char* backend, double sr, double capSec, double irSec, int block, double seconds)
{
    using C = fc::MatrixConvolverNupc<Fft>;
    const int cap = (int) std::ceil (capSec * sr), len = (int) std::lround (irSec * sr);
    const int xfade = (int) std::lround (0.05 * sr);   // CabConvolver's 50 ms
    C mc;
    if (! mc.prepare (128, cap, xfade, 2)) { std::printf ("prepare refused\n"); return; }
    std::vector<float> ir[4];
    Lcg r { 1234 };
    for (auto& h : ir)
    {
        h.resize ((std::size_t) len);
        float env = 0.5f;
        const float k = std::pow (1.0e-3f, 1.0f / (float) std::max (1, len));   // -60 dB over the IR
        for (float& v : h) { v = env * r.next(); env *= k; }
    }
    const float* banks[4] { ir[0].data(), ir[1].data(), ir[2].data(), ir[3].data() };

    std::vector<float> xl ((std::size_t) block), xr ((std::size_t) block), yl ((std::size_t) block), yr ((std::size_t) block);
    auto call = [&]() -> double
    {
        for (int i = 0; i < block; ++i) { xl[(std::size_t) i] = 0.25f * r.next(); xr[(std::size_t) i] = 0.25f * r.next(); }
        const float* in[2] { xl.data(), xr.data() };
        float* out[2] { yl.data(), yr.data() };
        const auto t0 = std::chrono::steady_clock::now();
        const bool ok = mc.process (in, out, 2, block);
        const auto t1 = std::chrono::steady_clock::now();
        if (! ok) std::abort();
        return std::chrono::duration<double, std::micro> (t1 - t0).count();
    };

    (void) mc.setOperator (C::Topology::Full, banks, 4, len);
    for (long s = 0; s < (long) cap + 4 * xfade + 8192; s += block) (void) call();   // fill the history, settle

    std::vector<double> idle;
    const long calls = (long) (seconds * sr / block);
    idle.reserve ((std::size_t) calls);
    for (long c = 0; c < calls; ++c) idle.push_back (call());

    std::vector<double> swap;
    long busyCalls = 0, busyMax = 0;
    const int swaps = std::max (8, (int) (seconds * sr / (xfade + 8192)));
    for (int k = 0; k < swaps; ++k)
    {
        (void) mc.setOperator (C::Topology::Full, banks, 4, len);   // the message thread's work: untimed
        long busy = 0;
        do { swap.push_back (call()); ++busy; } while (mc.isBusy());
        busyCalls += busy; busyMax = std::max (busyMax, busy);
        const int gap = 1 + (int) ((unsigned) (r.next() * 1000.0f + 1000.0f) % 2048u);   // the next publish at another phase
        for (int s = 0; s < gap; s += block) (void) call();
    }

    const Worst a = worstOf (idle), f = worstOf (swap);
    const double budget = 1e6 * block / sr, ms = 1e3 * block / sr;
    std::printf ("  %-6s IR %6.3f s  block %4d (budget %6.1f us) | idle mean %6.1f p99.9 %7.1f max %7.1f us (%5.1f%%) | swap mean %6.1f p99.9 %7.1f max %7.1f us (%5.1f%%) | busy %5.1f ms mean, %5.1f max\n",
                 backend, irSec, block, budget, a.meanUs, a.p999Us, a.maxUs, 100.0 * a.maxUs / budget,
                 f.meanUs, f.p999Us, f.maxUs, 100.0 * f.maxUs / budget,
                 ms * (double) busyCalls / (double) swaps, ms * (double) busyMax);
}
} // namespace

int main (int argc, char** argv)
{
    felitronics::core::ScopedFlushToZero ftz;
    const double seconds = argc > 1 ? std::atof (argv[1]) : 20.0;
    std::printf ("MatrixConvolverNupc worst process() call: Full, cap 20 s, 48 kHz, 50 ms fade — out of a swap, and from publish to the fade's end\n");
    // 2.773 s and 2.816 s are the threshold's edge for Full at 48 kHz: 133120 taps (P = 64, 252 MACs a hop, the old
    // path — the worst unspread swap in the tree, lineareq's maximum) and 135168 (P = 65, 256, spread).
    for (double irSec : { 20.0, 4.0, 133120.0 / 48000.0, 135168.0 / 48000.0, 2.0 })
        for (int block : { 64, 128, 256 })
        {
#if defined(FELITRONICS_WITH_PFFFT)
            row<felitronics::fftpffft::PffftRealFft> ("pffft", 48000.0, 20.0, irSec, block, seconds);
#endif
            if (irSec >= 20.0) row<felitronics::core::fft::DefaultRealFft> ("scalar", 48000.0, 20.0, irSec, block, seconds);
        }
    return 0;
}

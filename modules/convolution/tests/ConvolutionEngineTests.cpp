// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa. Part of felitronics-core — see LICENSE.

// JUCE-free self-tests for the production ConvolutionEngine: a CLICK-FREE IR swap (no derivative spike)
// that CONVERGES to the new IR's response, zero latency, and no-allocation-in-process().

#include <felitronics_test.h>
#include <felitronics/convolution/ConvolutionEngine.h>

#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

static std::atomic<long> g_allocs { 0 };
void* operator new      (std::size_t s) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s ? s : 1); }
void* operator new[]    (std::size_t s) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s ? s : 1); }
void  operator delete   (void* p) noexcept { std::free (p); }
void  operator delete[] (void* p) noexcept { std::free (p); }
void  operator delete   (void* p, std::size_t) noexcept { std::free (p); }
void  operator delete[] (void* p, std::size_t) noexcept { std::free (p); }

using namespace felitronics;

struct Lcg { unsigned long s; float next() { s = s * 6364136223846793005ULL + 1442695040888963407ULL; return (float) ((s >> 40) & 0xffff) / 32768.0f - 1.0f; } };

static double maxDeriv (const std::vector<float>& y, int from, int to)
{
    double m = 0.0;
    for (int i = from + 1; i < to && i < (int) y.size(); ++i) m = std::max (m, (double) std::fabs (y[i] - y[i - 1]));
    return m;
}

int main()
{
    std::printf ("felitronics::convolution ConvolutionEngine tests\n");
    const double sr = 48000.0;
    const int P = 64, irMax = 400, irLen = 200, xfade = 256, n = 3000;

    Lcg r { 31 };
    std::vector<float> irA (irLen), irB (irLen);
    for (auto& v : irA) v = 0.1f * r.next();
    for (auto& v : irB) v = 0.1f * r.next();
    std::vector<float> x (n);
    const double f = 500.0;
    for (int i = 0; i < n; ++i) x[i] = (float) (0.6 * std::sin (2.0 * core::kPi * f * i / sr));   // smooth → a click would show

    // --- click-free swap that converges to IR B ---
    test::group ("ConvolutionEngine click-free swap → converges to new IR");
    {
        convolution::ConvolutionEngine<> eng; test::ok (eng.prepare (P, irMax, xfade), "prepare");
        test::ok (convolution::ConvolutionEngine<>::latencySamples() == 0, "zero latency");

        std::vector<float> y (n, 0.0f);
        eng.setIr (irA.data(), irLen);                 // first load (fades in from silence)
        eng.process (x.data(), y.data(), 1400);        // settle on IR A
        const bool swapOk = eng.setIr (irB.data(), irLen);
        eng.process (x.data() + 1400, y.data() + 1400, n - 1400);   // crossfade A→B + settle
        test::ok (swapOk, "swap accepted while idle");

        // continuity: the derivative across the swap must not spike vs the steady region (no click).
        const double steady = maxDeriv (y, 900, 1300);
        const double swapRgn = maxDeriv (y, 1400, 1400 + xfade + 100);
        test::ok (swapRgn < 3.0 * steady + 1e-6, "no click — swap derivative ~ steady derivative");

        // convergence: well after the swap, output == a fresh IR-B convolver fed the same input.
        convolution::PartitionedConvolver<> ref; ref.prepare (P, irMax); ref.setIr (irB.data(), irLen);
        std::vector<float> yref (n, 0.0f); ref.process (x.data(), yref.data(), n);
        double maxErr = 0.0; for (int i = 2200; i < n; ++i) maxErr = std::max (maxErr, (double) std::fabs (y[i] - yref[i]));
        test::ok (maxErr < 3e-3, "converges to the new IR's steady response");
    }

    // --- no allocation during process() (incl. crossing a swap) ---
    test::group ("ConvolutionEngine no-alloc in process()");
    {
        convolution::ConvolutionEngine<> eng; eng.prepare (P, irMax, xfade);
        std::vector<float> in (512, 0.2f), out (512, 0.0f);
        eng.setIr (irA.data(), irLen);
        eng.process (in.data(), out.data(), 512);      // consume the initial fade-in
        eng.setIr (irB.data(), irLen);                 // arm a swap (build is message-thread, before the snapshot)
        const long before = g_allocs.load();
        eng.process (in.data(), out.data(), 512);      // crosses the crossfade
        eng.process (in.data(), out.data(), 512);
        const long after = g_allocs.load();
        test::ok (after == before, "process() performed zero heap allocations (even across a swap)");
    }

    return test::report();
}

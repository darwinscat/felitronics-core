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

    // --- STEREO lockstep: identical input + broadcast IR ⇒ L and R are bit-identical THROUGH a swap ---
    // (the whole point of one shared xfadePos; a per-channel race would decorrelate L/R during the fade)
    test::group ("ConvolutionEngine stereo lockstep (L==R through a swap)");
    {
        convolution::ConvolutionEngine<> eng;                      // MaxChannels defaults to 2
        test::ok (eng.prepare (P, irMax, xfade, 2), "prepare stereo (2 ch)");
        test::ok (eng.numChannels() == 2, "configured for 2 channels");

        std::vector<float> l (n, 0.0f), rr (n, 0.0f);
        eng.setIr (irA.data(), irLen);                             // broadcast mono IR to L & R
        {
            const float* in[2]  { x.data(), x.data() };            // same input on both channels
            float*       out[2] { l.data(), rr.data() };
            eng.process (in, out, 2, 1400);                        // settle on A
        }
        eng.setIr (irB.data(), irLen);                             // arm swap on both channels at once
        {
            const float* in[2]  { x.data() + 1400, x.data() + 1400 };
            float*       out[2] { l.data() + 1400, rr.data() + 1400 };
            eng.process (in, out, 2, n - 1400);                    // crossfade A→B
        }
        double maxLR = 0.0;
        for (int i = 0; i < n; ++i) maxLR = std::max (maxLR, (double) std::fabs (l[i] - rr[i]));
        test::ok (maxLR == 0.0, "L and R bit-identical every sample (incl. crossfade) → lockstep");
    }

    // --- stereo broadcast == the mono engine, bit-for-bit (incl. the fade-in from silence) ---
    test::group ("ConvolutionEngine stereo broadcast == mono");
    {
        convolution::ConvolutionEngine<> st; st.prepare (P, irMax, xfade, 2);
        convolution::ConvolutionEngine<> mo; mo.prepare (P, irMax, xfade, 1);
        std::vector<float> sl (n, 0.0f), srr (n, 0.0f), mout (n, 0.0f);
        const float* sin[2] { x.data(), x.data() }; float* sout[2] { sl.data(), srr.data() };
        st.setIr (irA.data(), irLen); st.process (sin, sout, 2, n);
        mo.setIr (irA.data(), irLen); mo.process (x.data(), mout.data(), n);
        double mErr = 0.0; for (int i = 0; i < n; ++i) mErr = std::max (mErr, (double) std::fabs (sl[i] - mout[i]));
        test::ok (mErr == 0.0, "broadcast stereo channel matches the mono engine bit-for-bit");
    }

    // --- per-channel (true-stereo) IR: L←irA, R←irB ⇒ distinct outputs that each converge ---
    test::group ("ConvolutionEngine per-channel IR");
    {
        convolution::ConvolutionEngine<> eng; eng.prepare (P, irMax, xfade, 2);
        std::vector<float> l (n, 0.0f), rr (n, 0.0f);
        const float* irs[2] { irA.data(), irB.data() };
        eng.setIr (irs, 2, irLen);                                 // L gets A, R gets B (one armed swap)
        const float* in[2] { x.data(), x.data() }; float* out[2] { l.data(), rr.data() };
        eng.process (in, out, 2, n);

        convolution::PartitionedConvolver<> rA, rB; rA.prepare (P, irMax); rB.prepare (P, irMax);
        rA.setIr (irA.data(), irLen); rB.setIr (irB.data(), irLen);
        std::vector<float> yA (n, 0.0f), yB (n, 0.0f); rA.process (x.data(), yA.data(), n); rB.process (x.data(), yB.data(), n);
        double eL = 0.0, eR = 0.0, lr = 0.0;
        for (int i = 2200; i < n; ++i)
        {
            eL = std::max (eL, (double) std::fabs (l[i]  - yA[i]));
            eR = std::max (eR, (double) std::fabs (rr[i] - yB[i]));
            lr = std::max (lr, (double) std::fabs (l[i]  - rr[i]));
        }
        test::ok (eL < 3e-3, "L converges to IR-A response");
        test::ok (eR < 3e-3, "R converges to IR-B response");
        test::ok (lr > 1e-3, "L and R differ (true-stereo IR actually applied per channel)");
    }

    return test::report();
}

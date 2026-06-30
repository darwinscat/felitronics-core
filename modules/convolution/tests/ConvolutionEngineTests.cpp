// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// JUCE-free self-tests for the production ConvolutionEngine: a CLICK-FREE IR swap (no derivative spike)
// that CONVERGES to the new IR's response, zero latency, and no-allocation-in-process().

#include <felitronics_test.h>
#include <felitronics/convolution/ConvolutionEngine.h>
#include <felitronics/convolution/PartitionedConvolver.h>   // reference convolver for the null tests

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

struct Lcg { unsigned long long s; float next() { s = s * 6364136223846793005ULL + 1442695040888963407ULL; return (float) ((s >> 40) & 0xffff) / 32768.0f - 1.0f; } };

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
        test::okNoAlloc (after == before, "process() performed zero heap allocations (even across a swap)");
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

    // --- DESIGN B: a WARM swap settles within the SHORT fade, not a full FIR length ---
    // (the whole point of the shared-FDL rework: dragging a band must not lag ~seconds)
    test::group ("ConvolutionEngine warm swap settles within the short fade (not the FIR length)");
    {
        const int Pw = 64, irMaxW = 2048, irLenW = 1500, warmXf = 128, nw = 2900;
        const int coldXf = ((irMaxW - Pw + Pw - 1) / Pw) * Pw;     // = maxParts·P (the cold-prime length)
        Lcg rw { 7 };
        std::vector<float> iA (irLenW), iB (irLenW), xw (nw);
        for (auto& v : iA) v = 0.05f * rw.next();
        for (auto& v : iB) v = 0.05f * rw.next();
        for (auto& v : xw) v = 0.5f  * rw.next();

        convolution::ConvolutionEngine<> eng; eng.prepare (Pw, irMaxW, warmXf);
        std::vector<float> yw (nw, 0.0f);
        const int sw = coldXf + 200;                               // swap well after the shared history is warm
        eng.setIr (iA.data(), irLenW); eng.process (xw.data(), yw.data(), sw);
        const bool warmOk = eng.setIr (iB.data(), irLenW);
        eng.process (xw.data() + sw, yw.data() + sw, nw - sw);
        test::ok (warmOk, "warm swap accepted");

        convolution::PartitionedConvolver<> refB; refB.prepare (Pw, irMaxW); refB.setIr (iB.data(), irLenW);
        std::vector<float> yB (nw, 0.0f); refB.process (xw.data(), yB.data(), nw);
        const int settled = sw + warmXf + Pw + 4;                  // just past the SHORT fade + one chunk
        double e = 0.0; for (int i = settled; i < nw; ++i) e = std::max (e, (double) std::fabs (yw[i] - yB[i]));
        test::ok (e < 2e-3, "matches the new IR right after the short fade — no N-length lag");
    }

    // --- DESIGN B: cold first prime is LONG; subsequent warm swaps are SHORT (isBusy timing) ---
    test::group ("ConvolutionEngine cold prime long, warm swap short");
    {
        const int Pc = 64, irMaxC = 1024, irLenC = 700, warmXf = 128;
        const int coldXf = ((irMaxC - Pc + Pc - 1) / Pc) * Pc;
        Lcg rc { 99 };
        std::vector<float> iA (irLenC), iB (irLenC), xc (4000);
        for (auto& v : iA) v = 0.05f * rc.next();
        for (auto& v : iB) v = 0.05f * rc.next();
        for (auto& v : xc) v = 0.5f  * rc.next();

        convolution::ConvolutionEngine<> eng; eng.prepare (Pc, irMaxC, warmXf);
        std::vector<float> y (4000, 0.0f); int pos = 0;
        auto run = [&] (int k) { eng.process (xc.data() + pos, y.data() + pos, k); pos += k; };
        eng.setIr (iA.data(), irLenC);
        run (warmXf + 8);   test::ok (eng.isBusy(),   "cold first prime still busy past the short-fade length");
        run (coldXf);       test::ok (! eng.isBusy(), "cold prime finished (≈ tail length)");
        test::ok (eng.setIr (iB.data(), irLenC), "warm swap accepted");
        run (warmXf + Pc + 8); test::ok (! eng.isBusy(), "warm swap finished within the short fade");
    }

    // --- DESIGN B: the swap PRIMES the new slot's tail (no ≤P dip from a zeroed tail) ---
    test::group ("ConvolutionEngine swap primes the new slot tail (no dip)");
    {
        const int Ps = 64, irMaxS = 400, irLenS = 200, warmXf = 128;
        const int coldXf = ((irMaxS - Ps + Ps - 1) / Ps) * Ps;
        std::vector<float> zA (irLenS, 0.0f);                       // IR A: silence
        std::vector<float> dB (irLenS, 0.0f); dB[Ps] = 1.0f;        // IR B: one tail tap → delay by P (steady tail = 1 for DC in)
        std::vector<float> ones (4000, 1.0f), y (4000, 0.0f);

        convolution::ConvolutionEngine<> eng; eng.prepare (Ps, irMaxS, warmXf);
        int pos = 0; auto run = [&] (int k) { eng.process (ones.data() + pos, y.data() + pos, k); pos += k; };
        eng.setIr (zA.data(), irLenS); run (coldXf + 4 * Ps);       // warm the FDL with the DC input (output stays 0)
        const int sw = pos;
        eng.setIr (dB.data(), irLenS);                             // warm swap → short fade; B's tail must be primed to 1
        run (warmXf + 2 * Ps);

        double e = 0.0;
        for (int m = 0; m < warmXf; ++m)
        {
            const float t = (float) m / (float) (warmXf - 1);
            const float expect = t * t * (3.0f - 2.0f * t);        // smoothstep weight on B (whose primed tail = 1)
            e = std::max (e, (double) std::fabs (y[sw + m] - expect));
        }
        test::ok (e < 2e-3, "blend follows smoothstep with B's tail live from sample 0 (a zeroed tail would dip ~0.5)");
    }

    // --- DESIGN B: reset() flushes history + re-arms a cold prime, but KEEPS the live IR ---
    test::group ("ConvolutionEngine reset keeps the live IR but re-arms cold");
    {
        const int Pr = 64, irMaxR = 1024, irLenR = 700, warmXf = 128;
        const int coldXf = ((irMaxR - Pr + Pr - 1) / Pr) * Pr;
        Lcg rr3 { 5 };
        std::vector<float> iA (irLenR), iB (irLenR), iC (irLenR), xr (4000);
        for (auto& v : iA) v = 0.05f * rr3.next();
        for (auto& v : iB) v = 0.05f * rr3.next();
        for (auto& v : iC) v = 0.05f * rr3.next();
        for (auto& v : xr) v = 0.5f  * rr3.next();
        iA[0] = 0.70f; iB[0] = -0.40f;                              // distinct head taps → tells which IR is live

        convolution::ConvolutionEngine<> eng; eng.prepare (Pr, irMaxR, warmXf);
        std::vector<float> y (4000, 0.0f); int pos = 0;
        auto run = [&] (int k) { eng.process (xr.data() + pos, y.data() + pos, k); pos += k; };
        eng.setIr (iA.data(), irLenR); run (coldXf + 200);          // cold-prime A → warm
        test::ok (eng.setIr (iB.data(), irLenR), "swap to B accepted (warm)");
        run (warmXf + Pr + 8);
        test::ok (! eng.isBusy(), "warm swap to B done — B is the live slot (cur_ flipped)");

        eng.reset();                                                // flush history; B must stay live
        float imp[1] { 1.0f }, oimp[1] { 0.0f };
        eng.process (imp, oimp, 1);                                 // first post-reset sample = head[live]·δ = iB[0]
        test::ok (std::fabs (oimp[0] - iB[0]) < 1e-6, "reset kept the live IR (B's head), did not revert to A");
        test::ok (std::fabs (oimp[0] - iA[0]) > 1e-3, "  …and it is NOT A");
        test::ok (eng.setIr (iC.data(), irLenR), "post-reset swap accepted");
        run (warmXf + 8);
        test::ok (eng.isBusy(), "post-reset swap is COLD again (busy past the short fade) → history re-armed");
    }

    return test::report();
}

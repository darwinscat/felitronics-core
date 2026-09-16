// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// JUCE-free self-tests for MatrixConvolverNupc — the non-uniform (Gardner) sibling of MatrixConvolver: all four
// topologies (mono / LRDiag / MSDiag / Full) + the CLICK-FREE 2-slot warm crossfade swap. Reference-NULL against
// the proven PartitionedConvolver, but — like MatrixConvolver's own tests — in a SETTLED window past the cold-
// prime crossfade (the first activation ramps the output in over the slowest stage's FDL fill, click-free). The
// crossfade uses a small maxIrSamples in the swap tests so the cold prime is short and a WARM swap is cheap to
// reach. Proves: per-topology routing == the time-domain reference; channel isolation; zero latency (mid-stream
// impulse, un-shifted); the cold prime ramps in click-free; a warm swap is click-free AND honours the warm
// history (differs from a cold-started instance); topology switches both directions; no-alloc; in-place.

#include <felitronics_test.h>
#include <alloc_counter.h>   // installs the allocation counter: EVERY form of `new`, over-aligned included
#include <felitronics/convolution/MatrixConvolverNupc.h>
#include <felitronics/convolution/PartitionedConvolver.h>
#include <felitronics/core/Math.h>   // core::sameBits — a bitwise compare, not ==

#include <cmath>
#include <cstdio>
#include <cstdlib>
#if defined(_WIN32)
 #include <malloc.h>
#endif
#include <string>
#include <utility>
#include <vector>

using namespace felitronics;
using MCN = convolution::MatrixConvolverNupc<>;

struct Lcg { unsigned long long s; float next() { s = s * 6364136223846793005ULL + 1442695040888963407ULL; return (float) ((s >> 40) & 0xffff) / 32768.0f - 1.0f; } };

static std::string sci (double v) { char b[32]; std::snprintf (b, sizeof b, "%.3e", v); return b; }

static std::vector<float> convRef (const std::vector<float>& x, const std::vector<float>& h, int maxIr)
{
    convolution::PartitionedConvolver<> pc; pc.prepare (128, maxIr); pc.setIr (h.data(), (int) h.size());
    std::vector<float> y (x.size(), 0.0f); felitronics::test::run (pc.process (x.data(), y.data(), (int) x.size()));
    return y;
}
static double maxDiff (const std::vector<float>& a, const std::vector<float>& b, int from, int to)
{ double m = 0.0; for (int i = from; i < to; ++i) m = std::max (m, (double) std::fabs (a[(std::size_t) i] - b[(std::size_t) i])); return m; }
static double peak (const std::vector<float>& a) { double m = 0.0; for (float v : a) m = std::max (m, (double) std::fabs (v)); return m; }
static double maxDeriv (const std::vector<float>& y, int from, int to)
{ double m = 0.0; for (int i = from + 1; i < to; ++i) m = std::max (m, (double) std::fabs (y[(std::size_t) i] - y[(std::size_t) (i - 1)])); return m; }

static void runStereo (MCN& mc, const std::vector<float>& xL, const std::vector<float>& xR, std::vector<float>& yL, std::vector<float>& yR, int block)
{
    yL.assign (xL.size(), 0.0f); yR.assign (xR.size(), 0.0f);
    for (int i = 0; i < (int) xL.size(); i += block) { const int m = std::min (block, (int) xL.size() - i); const float* in[2] { &xL[(std::size_t) i], &xR[(std::size_t) i] }; float* out[2] { &yL[(std::size_t) i], &yR[(std::size_t) i] }; felitronics::test::run (mc.process (in, out, 2, m)); }
}
static void runMono (MCN& mc, const std::vector<float>& x, std::vector<float>& y, int block)
{
    y.assign (x.size(), 0.0f);
    for (int i = 0; i < (int) x.size(); i += block) { const int m = std::min (block, (int) x.size() - i); const float* in[1] { &x[(std::size_t) i] }; float* out[1] { &y[(std::size_t) i] }; felitronics::test::run (mc.process (in, out, 1, m)); }
}

int main()
{
    std::printf ("felitronics::convolution MatrixConvolverNupc tests (all topologies + click-free crossfade)\n");
    // Modest maxIr → a short cold prime; `settled` is safely past the first-activation crossfade ramp.
    const int maxIr = 8192, settled = 7000;

    // --- per-topology routing == the time-domain reference (in the settled window, past the cold prime) ---
    test::group ("routing == reference (mono / LRDiag / MSDiag / Full), settled window");
    {
        const int L = 5000, N = L + 9000;
        Lcg r { 33 };
        std::vector<float> b0 ((std::size_t) L), b1 ((std::size_t) L), b2 ((std::size_t) L), b3 ((std::size_t) L);
        for (auto& v : b0) v = 0.05f * r.next(); for (auto& v : b1) v = 0.03f * r.next(); for (auto& v : b2) v = 0.03f * r.next(); for (auto& v : b3) v = 0.05f * r.next();
        std::vector<float> xL ((std::size_t) N), xR ((std::size_t) N); for (auto& v : xL) v = 0.3f * r.next(); for (auto& v : xR) v = 0.4f * r.next();

        { MCN mc; test::ok (mc.prepare (128, maxIr, 128, 1), "prepare mono"); mc.setIr (b0.data(), L);
          std::vector<float> y; runMono (mc, xL, y, 256);
          const std::vector<float> ref = convRef (xL, b0, maxIr);
          test::ok (maxDiff (y, ref, settled, N) / (peak (ref) + 1e-30) < 2e-4, "mono == b0 ∗ xL"); }

        { MCN mc; mc.prepare (128, maxIr, 128, 2); const float* bk[2] { b0.data(), b1.data() };
          test::ok (mc.setOperator (MCN::Topology::LRDiag, bk, 2, L), "setOperator LRDiag");
          std::vector<float> yL, yR; runStereo (mc, xL, xR, yL, yR, 256);
          const std::vector<float> rL = convRef (xL, b0, maxIr), rR = convRef (xR, b1, maxIr);
          test::ok (maxDiff (yL, rL, settled, N) / (peak (rL) + 1e-30) < 2e-4, "LRDiag yL == b0 ∗ xL");
          test::ok (maxDiff (yR, rR, settled, N) / (peak (rR) + 1e-30) < 2e-4, "LRDiag yR == b1 ∗ xR"); }

        { MCN mc; mc.prepare (128, maxIr, 128, 2); const float* bk[2] { b0.data(), b1.data() };
          test::ok (mc.setOperator (MCN::Topology::MSDiag, bk, 2, L), "setOperator MSDiag");
          std::vector<float> yL, yR; runStereo (mc, xL, xR, yL, yR, 256);
          std::vector<float> mm ((std::size_t) N), ss ((std::size_t) N); for (int i = 0; i < N; ++i) { mm[(std::size_t) i] = 0.5f * (xL[(std::size_t) i] + xR[(std::size_t) i]); ss[(std::size_t) i] = 0.5f * (xL[(std::size_t) i] - xR[(std::size_t) i]); }
          const std::vector<float> yM = convRef (mm, b0, maxIr), yS = convRef (ss, b1, maxIr);
          std::vector<float> rL ((std::size_t) N), rR ((std::size_t) N); for (int i = 0; i < N; ++i) { rL[(std::size_t) i] = yM[(std::size_t) i] + yS[(std::size_t) i]; rR[(std::size_t) i] = yM[(std::size_t) i] - yS[(std::size_t) i]; }
          test::ok (maxDiff (yL, rL, settled, N) / (peak (rL) + 1e-30) < 2e-4, "MSDiag yL == yM+yS (per-stage M/S views)");
          test::ok (maxDiff (yR, rR, settled, N) / (peak (rR) + 1e-30) < 2e-4, "MSDiag yR == yM−yS"); }

        { MCN mc; mc.prepare (128, maxIr, 128, 2); const float* bk[4] { b0.data(), b1.data(), b2.data(), b3.data() };
          test::ok (mc.setOperator (MCN::Topology::Full, bk, 4, L), "setOperator Full");
          std::vector<float> yL, yR; runStereo (mc, xL, xR, yL, yR, 256);
          const std::vector<float> cLL = convRef (xL, b0, maxIr), cLR = convRef (xR, b1, maxIr), cRL = convRef (xL, b2, maxIr), cRR = convRef (xR, b3, maxIr);
          std::vector<float> rL ((std::size_t) N), rR ((std::size_t) N); for (int i = 0; i < N; ++i) { rL[(std::size_t) i] = cLL[(std::size_t) i] + cLR[(std::size_t) i]; rR[(std::size_t) i] = cRL[(std::size_t) i] + cRR[(std::size_t) i]; }
          test::ok (maxDiff (yL, rL, settled, N) / (peak (rL) + 1e-30) < 2e-4, "Full yL == LL∗xL + LR∗xR");
          test::ok (maxDiff (yR, rR, settled, N) / (peak (rR) + 1e-30) < 2e-4, "Full yR == RL∗xL + RR∗xR"); }
    }

    // --- setIr() on a STEREO instance broadcasts the mono IR to BOTH channels (regression: it returned false) ---
    test::group ("setIr broadcasts to both channels on a stereo instance");
    {
        const int L = 4000, N = L + 9000;
        Lcg r { 91 };
        std::vector<float> h ((std::size_t) L); for (auto& v : h) v = 0.05f * r.next();
        std::vector<float> xL ((std::size_t) N), xR ((std::size_t) N); for (auto& v : xL) v = 0.30f * r.next(); for (auto& v : xR) v = 0.35f * r.next();
        MCN mc; test::ok (mc.prepare (128, maxIr, 128, 2), "prepare stereo");
        test::ok (mc.setIr (h.data(), L), "setIr accepted on a stereo instance");   // regression: used to return false (numBanks 1 < 2)
        std::vector<float> yL, yR; runStereo (mc, xL, xR, yL, yR, 256);
        const std::vector<float> rL = convRef (xL, h, maxIr), rR = convRef (xR, h, maxIr);
        test::ok (maxDiff (yL, rL, settled, N) / (peak (rL) + 1e-30) < 2e-4, "stereo setIr: yL == h ∗ xL");
        test::ok (maxDiff (yR, rR, settled, N) / (peak (rR) + 1e-30) < 2e-4, "stereo setIr: yR == h ∗ xR");
    }

    // --- F1 regression: a render SHORTER than the (removed) cold-prime length must be UN-attenuated. The old cold
    //     fade scaled to the IR's reach (~6000 samples here); for a render of 4000 it never completed, so the WHOLE
    //     output stayed on the attenuating ramp (~10 dB down, worse earlier). The cold-started FDL is already the
    //     exact causal convolution, so now only the short warm fade (128) applies and the output is exact past it. ---
    test::group ("short render is un-attenuated (no long cold fade)");
    {
        const int bigIr = 131072, L = 6000, Nr = 4000, xf = 128;   // IR reach (6000) > render (4000): old fade never finished
        Lcg r { 4242 };
        std::vector<float> h ((std::size_t) L); for (auto& v : h) v = 0.05f * r.next();
        std::vector<float> x ((std::size_t) Nr); for (auto& v : x) v = 0.3f * r.next();
        MCN mc; test::ok (mc.prepare (128, bigIr, xf, 1), "prepare mono, 131072-tap capacity");
        mc.setIr (h.data(), L);
        std::vector<float> y; runMono (mc, x, y, 256);
        const std::vector<float> ref = convRef (x, h, bigIr);
        const double rel = maxDiff (y, ref, xf + 256, Nr) / (peak (ref) + 1e-30);   // exact past the SHORT fade only
        test::ok (rel < 2e-4, "output == reference right after the short fade (the whole render would be attenuated on the old cold ramp)");
    }

    // --- channel isolation (LRDiag): R input all-zero → R output exactly silent (blend of two zero paths) ---
    test::group ("channel isolation: R input zero → R output silent");
    {
        const int L = 4000, N = L + 9000;
        Lcg r { 5 };
        std::vector<float> hL ((std::size_t) L), hR ((std::size_t) L); for (auto& v : hL) v = 0.05f * r.next(); for (auto& v : hR) v = 0.05f * r.next();
        std::vector<float> xL ((std::size_t) N), xR ((std::size_t) N, 0.0f); for (auto& v : xL) v = 0.3f * r.next();
        MCN mc; mc.prepare (128, maxIr, 128, 2); const float* bk[2] { hL.data(), hR.data() }; mc.setOperator (MCN::Topology::LRDiag, bk, 2, L);
        std::vector<float> yL, yR; runStereo (mc, xL, xR, yL, yR, 128);
        test::ok (peak (yR) < 1e-9, "yR exactly silent (no cross-channel bleed, even through the crossfade)");
        test::ok (peak (yL) > 1e-2, "yL responds to xL");
    }

    // --- zero latency: a mid-stream impulse (past the cold prime) reproduces the IR with NO shift ---
    test::group ("TRUE zero latency (mid-stream impulse == IR, un-shifted)");
    {
        const int L = 1500, t0 = 7000, N = t0 + L + 200;
        Lcg r { 909 };
        std::vector<float> h ((std::size_t) L); for (auto& v : h) v = 0.05f * r.next();
        std::vector<float> x ((std::size_t) N, 0.0f); for (int i = 0; i < t0; ++i) x[(std::size_t) i] = 0.2f * r.next();   // warm past the cold prime
        x[(std::size_t) t0] = 1.0f;                                                                                       // then an isolated impulse
        MCN mc; mc.prepare (128, maxIr, 128, 1); mc.setIr (h.data(), L);
        std::vector<float> y; runMono (mc, x, y, 100);
        // subtract the no-impulse response so only the impulse's contribution remains
        std::vector<float> x2 = x; x2[(std::size_t) t0] = 0.0f;
        MCN mc2; mc2.prepare (128, maxIr, 128, 1); mc2.setIr (h.data(), L); std::vector<float> y2; runMono (mc2, x2, y2, 100);
        double e = 0.0; for (int i = 0; i < L; ++i) e = std::max (e, (double) std::fabs ((y[(std::size_t) (t0 + i)] - y2[(std::size_t) (t0 + i)]) - h[(std::size_t) i]));
        test::ok (e < 1e-4, "impulse response == h, first sample at t0 → zero added latency");
    }

    // --- the cold prime ramps in click-free (first activation: no discontinuity at the start) ---
    test::group ("cold prime is click-free (smooth ramp-in)");
    {
        const int L = 3000, N = 12000;
        Lcg r { 4141 };
        std::vector<float> h ((std::size_t) L); for (auto& v : h) v = 0.05f * r.next();
        std::vector<float> x ((std::size_t) N); for (auto& v : x) v = 0.3f * r.next();
        MCN mc; mc.prepare (128, maxIr, 128, 1); mc.setIr (h.data(), L);
        std::vector<float> y; runMono (mc, x, y, 256);
        const double steady = maxDeriv (y, settled, settled + 500);
        const double onset  = maxDeriv (y, 1, 400);                          // right at the ramp-in
        test::ok (onset < 4.0 * steady + 1e-6, "no click at the cold-prime onset (smoothstep ramp)");
    }

    // --- click-free WARM swap that honours the warm history (small maxIr → cheap warm-up) ---
    test::group ("warm swap: click-free + warm history");
    {
        const int mIr = 4096, L = 2000, N = 20000, T = 9000;   // T is well past the first-activation fade → a plain warm swap
        Lcg r { 3131 };
        std::vector<float> h1 ((std::size_t) L), h2 ((std::size_t) L); for (auto& v : h1) v = 0.05f * r.next(); for (auto& v : h2) v = 0.05f * r.next();
        std::vector<float> x ((std::size_t) N); for (auto& v : x) v = 0.3f * r.next();
        MCN mc; mc.prepare (128, mIr, 128, 1); mc.setIr (h1.data(), L);
        std::vector<float> y ((std::size_t) N, 0.0f); int swapAt = -1;
        for (int i = 0; i < N; i += 256)
        {
            if (swapAt < 0 && i >= T) { if (mc.setIr (h2.data(), L)) swapAt = i; }
            const int m = std::min (256, N - i); const float* in[1] { &x[(std::size_t) i] }; float* out[1] { &y[(std::size_t) i] }; felitronics::test::run (mc.process (in, out, 1, m));
        }
        test::ok (swapAt > 0, "warm swap accepted mid-stream");
        const double steady = maxDeriv (y, swapAt - 800, swapAt - 100);
        const double across = maxDeriv (y, swapAt - 10, swapAt + 400);
        test::ok (across < 4.0 * steady + 1e-6, "no click across the warm swap");
        const std::vector<float> ref2 = convRef (x, h2, mIr);
        test::ok (maxDiff (y, ref2, swapAt + 600, N) / (peak (ref2) + 1e-30) < 2e-4, "post-swap output == h2 ∗ x (warm)");

        // WARM proof: a cold instance fed only post-swap samples differs from the warm one where the pre-swap
        // history still matters (within one IR length of the swap, past the short fade).
        MCN cold; cold.prepare (128, mIr, 128, 1); cold.setIr (h2.data(), L);
        const int tail = N - swapAt; std::vector<float> cx ((std::size_t) tail), cy;
        for (int i = 0; i < tail; ++i) cx[(std::size_t) i] = x[(std::size_t) (swapAt + i)];
        runMono (cold, cx, cy, 256);
        double warmDiff = 0.0; for (int i = 300; i < 900; ++i) warmDiff = std::max (warmDiff, (double) std::fabs (y[(std::size_t) (swapAt + i)] - cy[(std::size_t) i]));
        test::ok (warmDiff > 1e-3, "warm swap ≠ a cold-started instance (the incoming operator read the warm FDL)");
    }

    // --- topology switch, both directions (warm, click-free) ---
    test::group ("topology switch MSDiag↔Full (warm, click-free)");
    {
        const int mIr = 4096, L = 2000, N = 22000, T = 9000;
        Lcg r { 909 };
        std::vector<float> hM ((std::size_t) L), hS ((std::size_t) L), LL ((std::size_t) L), LR ((std::size_t) L), RL ((std::size_t) L), RR ((std::size_t) L);
        for (auto& v : hM) v = 0.04f * r.next(); for (auto& v : hS) v = 0.04f * r.next();
        for (auto& v : LL) v = 0.04f * r.next(); for (auto& v : LR) v = 0.03f * r.next(); for (auto& v : RL) v = 0.03f * r.next(); for (auto& v : RR) v = 0.04f * r.next();
        std::vector<float> xL ((std::size_t) N), xR ((std::size_t) N); for (auto& v : xL) v = 0.3f * r.next(); for (auto& v : xR) v = 0.4f * r.next();
        for (int dir = 0; dir < 2; ++dir)   // dir 0: MSDiag→Full ; dir 1: Full→MSDiag
        {
            MCN mc; mc.prepare (128, mIr, 128, 2);
            if (dir == 0) { const float* b[2] { hM.data(), hS.data() }; mc.setOperator (MCN::Topology::MSDiag, b, 2, L); }
            else          { const float* b[4] { LL.data(), LR.data(), RL.data(), RR.data() }; mc.setOperator (MCN::Topology::Full, b, 4, L); }
            std::vector<float> yL ((std::size_t) N, 0.0f), yR ((std::size_t) N, 0.0f); int swapAt = -1;
            for (int i = 0; i < N; i += 256)
            {
                if (swapAt < 0 && i >= T)
                {
                    bool ok;
                    if (dir == 0) { const float* b[4] { LL.data(), LR.data(), RL.data(), RR.data() }; ok = mc.setOperator (MCN::Topology::Full, b, 4, L); }
                    else          { const float* b[2] { hM.data(), hS.data() };                       ok = mc.setOperator (MCN::Topology::MSDiag, b, 2, L); }
                    if (ok) swapAt = i;
                }
                const int m = std::min (256, N - i); const float* in[2] { &xL[(std::size_t) i], &xR[(std::size_t) i] }; float* out[2] { &yL[(std::size_t) i], &yR[(std::size_t) i] }; felitronics::test::run (mc.process (in, out, 2, m));
            }
            test::ok (swapAt > 0, dir == 0 ? "MSDiag→Full accepted" : "Full→MSDiag accepted");
            const double steady = maxDeriv (yL, swapAt - 800, swapAt - 100);
            const double across = maxDeriv (yL, swapAt - 10, swapAt + 400);
            test::ok (across < 4.0 * steady + 1e-6, dir == 0 ? "no click MSDiag→Full" : "no click Full→MSDiag");
            // post-switch converges to the destination topology
            std::vector<float> rL ((std::size_t) N), rR ((std::size_t) N);
            if (dir == 0) { const std::vector<float> cLL = convRef (xL, LL, mIr), cLR = convRef (xR, LR, mIr), cRL = convRef (xL, RL, mIr), cRR = convRef (xR, RR, mIr); for (int i = 0; i < N; ++i) { rL[(std::size_t) i] = cLL[(std::size_t) i] + cLR[(std::size_t) i]; rR[(std::size_t) i] = cRL[(std::size_t) i] + cRR[(std::size_t) i]; } }
            else          { std::vector<float> mm ((std::size_t) N), ss ((std::size_t) N); for (int i = 0; i < N; ++i) { mm[(std::size_t) i] = 0.5f * (xL[(std::size_t) i] + xR[(std::size_t) i]); ss[(std::size_t) i] = 0.5f * (xL[(std::size_t) i] - xR[(std::size_t) i]); } const std::vector<float> yM = convRef (mm, hM, mIr), yS = convRef (ss, hS, mIr); for (int i = 0; i < N; ++i) { rL[(std::size_t) i] = yM[(std::size_t) i] + yS[(std::size_t) i]; rR[(std::size_t) i] = yM[(std::size_t) i] - yS[(std::size_t) i]; } }
            test::ok (maxDiff (yL, rL, swapAt + 600, N) / (peak (rL) + 1e-30) < 2e-4 && maxDiff (yR, rR, swapAt + 600, N) / (peak (rR) + 1e-30) < 2e-4, dir == 0 ? "converges to Full" : "converges to MSDiag");
        }
    }

    // --- reset() in the MIDDLE of a crossfade ENDS it cleanly (the one fade path with no other coverage) ---
    // P88: "ends" is the fade, not the operator — reset() adopts the published one (the group below proves
    // which operator is live and that nothing of the old stream survives). Here: idle afterwards, finite,
    // and a fresh operator is accepted.
    test::group ("reset() mid-crossfade ends the fade cleanly");
    {
        const int mIr = 4096, L = 2000, N = 4000, xf = 1024;   // a long WARM fade so 512 samples lands mid-crossfade
        Lcg r { 321 };
        std::vector<float> h ((std::size_t) L); for (auto& v : h) v = 0.05f * r.next();
        std::vector<float> x ((std::size_t) N); for (auto& v : x) v = 0.3f * r.next();
        MCN mc; mc.prepare (128, mIr, xf, 1); mc.setIr (h.data(), L);
        std::vector<float> y ((std::size_t) N, 0.0f);
        for (int i = 0; i < 512; i += 256) { const int m = std::min (256, 512 - i); const float* in[1] { &x[(std::size_t) i] }; float* out[1] { &y[(std::size_t) i] }; felitronics::test::run (mc.process (in, out, 1, m)); }   // 512 samples into the 1024-sample fade
        test::ok (mc.isBusy(), "busy mid-crossfade");
        mc.reset();
        test::ok (! mc.isBusy(), "reset() ends the crossfade (idle)");
        std::vector<float> z (1000, 0.1f), zo (1000, 0.0f); const float* zin[1] { z.data() }; float* zout[1] { zo.data() }; felitronics::test::run (mc.process (zin, zout, 1, 1000));
        bool finite = true; for (float v : zo) if (! std::isfinite (v)) finite = false;
        test::ok (finite, "output finite + stable after a mid-fade reset");
        test::ok (mc.setIr (h.data(), L), "a fresh operator is accepted after the reset");
    }

    // --- 🔴 P88: reset() KEEPS AN OPERATOR PUBLISHED BUT NOT YET ADOPTED (the first instance of the class) ---
    // P85 closed the composite's path here by adding clearAudioState() and routing rigplayer::RigPlayer
    // through it; the VERB itself still dropped the publication, and lineareq — the live product path —
    // reaches this class through `reset()`. On that body this group reads the old operator in every cell,
    // a worst sample 1.927e+00 from the settled restart. Full reasoning in ConvolutionEngine::reset().
    // The oracle is law 11a INDEPENDENCE — two convolvers fed DIFFERENT audio of DIFFERENT length, one
    // restarted with the publication in flight and one after it settled, answer the next programme with
    // identical bits — for every TOPOLOGY, with four distinct banks and different L/R inputs (with L == R
    // an MSDiag side is silent and a routing slip after the adoption is invisible — a review round showed
    // `slot_[cur_].topo = LRDiag` surviving the first version of this group), and at BOTH slot parities,
    // since `cur_ = 0` is right half the time. The LRDiag cells are certified from OUTSIDE the class too.
    test::group ("MatrixConvolverNupc reset keeps a published-but-unadopted operator (independence, every topology)");
    {
        const int Pk = 128, irMaxK = 4096, irLenK = 2000, xfK = 1024;
        const int settleK = xfK + 4 * Pk + irMaxK;
        const int progN   = 8000;
        Lcg rk { 90210 };
        // operators 0 = A, 1 = A2 (the parity pre-swap), 2 = B (published), 3 = C (published after); 4 banks each
        std::vector<float> ops[4][4];
        for (int o = 0; o < 4; ++o)
            for (int b = 0; b < 4; ++b)
            {
                ops[o][b].assign ((std::size_t) irLenK, 0.0f);
                for (auto& v : ops[o][b]) v = 0.05f * rk.next();
                ops[o][b][0] = 0.1f * (float) (1 + o) * (float) (1 + b) * (((o + b) & 1) != 0 ? -1.0f : 1.0f);
            }
        std::vector<float> preD[2] { std::vector<float> ((std::size_t) 12000), std::vector<float> ((std::size_t) 12000) };
        std::vector<float> preR[2] { std::vector<float> ((std::size_t) 12600), std::vector<float> ((std::size_t) 12600) };
        std::vector<float> prog[2] { std::vector<float> ((std::size_t) progN), std::vector<float> ((std::size_t) progN) };
        for (auto& ch : preD) for (auto& v : ch) v = 0.5f * rk.next();
        for (auto& ch : preR) for (auto& v : ch) v = 0.4f * rk.next();
        for (auto& ch : prog) for (auto& v : ch) v = 0.3f * rk.next();

        // mode 0 = restart at Pending · 1 = restart mid-crossfade · 2 = the reference (publication settled)
        const auto render = [&] (int nch, MCN::Topology topo, int parity, int mode,
                                 std::vector<float>& L, std::vector<float>& R, bool* busyOut, bool* acceptOut)
        {
            MCN mc;
            felitronics::test::run (mc.prepare (Pk, irMaxK, xfK, nch));
            const std::vector<float>* pre = (mode == 2) ? preR : preD;
            const int preN = (int) pre[0].size();
            std::vector<float> jl ((std::size_t) preN, 0.0f), jr ((std::size_t) preN, 0.0f);
            int at = 0;
            const auto feed = [&] (int k)                          // wraps around the pre-roll if asked for more
            {
                for (int done = 0; done < k; )
                {
                    const int m = std::min (k - done, preN - at);
                    const float* in[2]  { pre[0].data() + at, pre[1].data() + at };
                    float*       out[2] { jl.data(), jr.data() };
                    felitronics::test::run (mc.process (in, out, nch, m));
                    done += m; at = (at + m) % preN;
                }
            };
            const auto publish = [&] (int o)
            {
                const float* b[4] { ops[o][0].data(), ops[o][1].data(), ops[o][2].data(), ops[o][3].data() };
                return mc.setOperator (topo, b, nch == 1 ? 1 : MCN::numBanksFor (topo), irLenK);
            };
            if (parity == 1) { test::ok (publish (1), "parity pre-swap accepted"); feed (settleK); }
            test::ok (publish (0), "operator A accepted");
            feed (preN);
            test::ok (publish (2), "operator B published");
            if (mode == 1) feed (xfK / 2);
            if (mode == 2) feed (settleK);
            mc.reset();
            if (busyOut != nullptr) *busyOut = mc.isBusy();
            L.assign ((std::size_t) progN, 0.0f); R.assign ((std::size_t) progN, 0.0f);
            {
                const float* in[2]  { prog[0].data(), prog[1].data() };
                float*       out[2] { L.data(), R.data() };
                felitronics::test::run (mc.process (in, out, nch, progN));
            }
            if (acceptOut != nullptr) *acceptOut = publish (3);
        };

        using Topo = MCN::Topology;
        const std::pair<Topo, const char*> topos[] { { Topo::LRDiag, "LRDiag" }, { Topo::MSDiag, "MSDiag" }, { Topo::Full, "Full" } };
        for (int nch = 1; nch <= 2; ++nch)
            for (const auto& [topo, topoName] : topos)
            {
                if (nch == 1 && topo != Topo::LRDiag) continue;        // mono: the topology is inert
                for (int parity = 0; parity <= 1; ++parity)
                {
                    std::vector<float> ref[2];
                    render (nch, topo, parity, 2, ref[0], ref[1], nullptr, nullptr);
                    const std::string cell = "  [" + std::to_string (nch) + " ch, " + topoName + ", slot parity " + std::to_string (parity) + "]";
                    for (int mode = 0; mode <= 1; ++mode)
                    {
                        std::vector<float> dut[2]; bool busy = true, accepts = false;
                        render (nch, topo, parity, mode, dut[0], dut[1], &busy, &accepts);
                        bool same = true; double worst = 0.0;
                        for (int c = 0; c < nch; ++c)
                            for (int i = 0; i < progN; ++i)
                            {
                                same  = same && core::sameBits (dut[c][(std::size_t) i], ref[c][(std::size_t) i]);
                                worst = std::max (worst, (double) std::fabs (dut[c][(std::size_t) i] - ref[c][(std::size_t) i]));
                            }
                        test::ok (same, std::string ("restarted at ") + (mode == 0 ? "Pending" : "Crossfading")
                                        + ", bit-identical to a settled restart" + cell + " (worst " + sci (worst) + ")");
                        test::ok (! busy,  "  …and Idle right after reset()" + cell);
                        test::ok (accepts, "  …and a new publication is accepted at once" + cell);

                        if (topo == Topo::LRDiag)          // LRDiag: yL = bank0 ∗ xL, yR = bank1 ∗ xR — checkable from outside
                            for (int c = 0; c < nch; ++c)
                            {
                                const std::vector<float> yB = convRef (prog[c], ops[2][c], irMaxK), yA = convRef (prog[c], ops[0][c], irMaxK);
                                const double eB = maxDiff (dut[c], yB, 0, progN) / (peak (yB) + 1e-30), eA = maxDiff (dut[c], yA, 0, progN) / (peak (yA) + 1e-30);
                                test::ok (eB < 2e-4, "  …and channel " + std::to_string (c) + " answers as the PUBLISHED operator" + cell + " (" + sci (eB) + ")");
                                test::ok (eA > 1e-2,  "  …and not as the previous one" + cell + " (" + sci (eA) + ")");
                            }
                    }
                }
            }

        // A STAGED but NOT PUBLISHED operator (state_ == 0) is none of reset()'s business: `cur_` must not
        // move, the banks must survive, and the later publishStaged() must still find it and fade it in.
        for (int nch = 1; nch <= 2; ++nch)
        {
            MCN mc;
            felitronics::test::run (mc.prepare (Pk, irMaxK, xfK, nch));
            std::vector<float> jl (preD[0].size(), 0.0f), jr (preD[0].size(), 0.0f);
            {
                const float* b[4] { ops[0][0].data(), ops[0][1].data(), ops[0][2].data(), ops[0][3].data() };
                test::ok (mc.setOperator (Topo::LRDiag, b, nch == 1 ? 1 : 2, irLenK), "A accepted");
                const float* in[2]  { preD[0].data(), preD[1].data() };
                float*       out[2] { jl.data(), jr.data() };
                felitronics::test::run (mc.process (in, out, nch, (int) preD[0].size()));
            }
            {
                const float* b[4] { ops[2][0].data(), ops[2][1].data(), ops[2][2].data(), ops[2][3].data() };
                test::ok (mc.stageOperator (Topo::LRDiag, b, nch == 1 ? 1 : 2, irLenK), "B staged, not published");
            }
            mc.reset();
            test::ok (! mc.isBusy(), "a staged-not-published operator leaves reset() idle  [" + std::to_string (nch) + " ch]");
            mc.publishStaged();
            std::vector<float> y[2] { std::vector<float> ((std::size_t) progN, 0.0f), std::vector<float> ((std::size_t) progN, 0.0f) };
            {
                const float* in[2]  { prog[0].data(), prog[1].data() };
                float*       out[2] { y[0].data(), y[1].data() };
                felitronics::test::run (mc.process (in, out, nch, progN));
            }
            for (int c = 0; c < nch; ++c)
            {
                const std::vector<float> yB = convRef (prog[c], ops[2][c], irMaxK);
                const double tail = maxDiff (y[c], yB, xfK + 1200, progN) / (peak (yB) + 1e-30);
                test::ok (tail < 2e-4, "the publish that follows the restart still lands on the staged operator  ["
                                          + std::to_string (nch) + " ch, channel " + std::to_string (c) + "] (" + sci (tail) + ")");
            }
        }
    }

    // --- 🔴 P88: clearAudioState() — the history alone, the swap left running FROM WHERE IT WAS ---
    // Two convolvers with different pasts, cleared at the same fade position, answer the same bits; the fade
    // then ENDS exactly when it would have (a review round showed `xfadePos_ = 0` inside this verb surviving
    // the first version of the group), lands on the published operator, and — the one thing this verb does
    // not give — the answer still depends on WHERE in the fade it was called. That is reset()'s promise.
    test::group ("MatrixConvolverNupc clearAudioState leaves the swap running and the history gone");
    {
        const int Pc = 128, irMaxC = 4096, irLenC = 2000, xfC = 4096;
        const int progC = 12000;
        Lcg rc { 5150 };
        std::vector<float> oA[2], oB[2];
        for (int b = 0; b < 2; ++b)
        {
            oA[b].assign ((std::size_t) irLenC, 0.0f); for (auto& v : oA[b]) v = 0.05f * rc.next();
            oB[b].assign ((std::size_t) irLenC, 0.0f); for (auto& v : oB[b]) v = 0.05f * rc.next();
        }
        std::vector<float> preD[2] { std::vector<float> ((std::size_t) 12000), std::vector<float> ((std::size_t) 12000) };
        std::vector<float> preR[2] { std::vector<float> ((std::size_t) 12600), std::vector<float> ((std::size_t) 12600) };
        std::vector<float> prog[2] { std::vector<float> ((std::size_t) progC), std::vector<float> ((std::size_t) progC) };
        for (auto& ch : preD) for (auto& v : ch) v = 0.5f * rc.next();
        for (auto& ch : preR) for (auto& v : ch) v = 0.4f * rc.next();
        for (auto& ch : prog) for (auto& v : ch) v = 0.3f * rc.next();

        // renders `prog` after a clear `intoFade` samples into the fade; reports busy just before and just
        // after the sample on which the fade must end
        const auto render = [&] (int nch, const std::vector<float>* pre, int intoFade, std::vector<float>* out,
                                 bool* busyAfterClear, bool* busyOneShort, bool* busyOnTime)
        {
            MCN mc;
            felitronics::test::run (mc.prepare (Pc, irMaxC, xfC, nch));
            const int preN = (int) pre[0].size();
            std::vector<float> jl ((std::size_t) preN, 0.0f), jr ((std::size_t) preN, 0.0f);
            const auto feed = [&] (int from, int k)
            {
                const float* in[2]  { pre[0].data() + from, pre[1].data() + from };
                float*       o[2]   { jl.data(), jr.data() };
                felitronics::test::run (mc.process (in, o, nch, k));
            };
            { const float* b[2] { oA[0].data(), oA[1].data() }; test::ok (mc.setOperator (MCN::Topology::LRDiag, b, nch == 1 ? 1 : 2, irLenC), "A accepted"); }
            feed (0, preN - xfC);
            { const float* b[2] { oB[0].data(), oB[1].data() }; test::ok (mc.setOperator (MCN::Topology::LRDiag, b, nch == 1 ? 1 : 2, irLenC), "B published"); }
            feed (preN - xfC, intoFade);
            mc.clearAudioState();
            *busyAfterClear = mc.isBusy();
            out[0].assign ((std::size_t) progC, 0.0f); out[1].assign ((std::size_t) progC, 0.0f);
            const int remaining = xfC - intoFade;                   // the fade ends on this sample, not later
            const auto run = [&] (int from, int k)
            {
                const float* in[2]  { prog[0].data() + from, prog[1].data() + from };
                float*       o[2]   { out[0].data() + from, out[1].data() + from };
                felitronics::test::run (mc.process (in, o, nch, k));
            };
            run (0, remaining - 1);          *busyOneShort = mc.isBusy();
            run (remaining - 1, 1);          *busyOnTime   = mc.isBusy();
            run (remaining, progC - remaining);
        };
        for (int nch = 1; nch <= 2; ++nch)
        {
            const std::string cell = "  [" + std::to_string (nch) + " ch]";
            std::vector<float> a[2], b[2], late[2];
            bool a0 = false, a1 = false, a2 = true, b0 = false, b1 = false, b2 = true, l0 = false, l1 = false, l2 = true;
            render (nch, preD, xfC / 4, a, &a0, &a1, &a2);
            render (nch, preR, xfC / 4, b, &b0, &b1, &b2);
            render (nch, preD, xfC / 2, late, &l0, &l1, &l2);
            test::ok (a0 && b0 && l0, "clearAudioState() leaves the crossfade in flight (still busy)" + cell);
            test::ok (a1 && b1 && l1 && ! a2 && ! b2 && ! l2, "…and it ends on the sample it was always going to end on, not later" + cell);
            bool same = true; double worst = 0.0, apart = 0.0;
            for (int c = 0; c < nch; ++c)
                for (int i = 0; i < progC; ++i)
                {
                    same  = same && core::sameBits (a[c][(std::size_t) i], b[c][(std::size_t) i]);
                    worst = std::max (worst, (double) std::fabs (a[c][(std::size_t) i] - b[c][(std::size_t) i]));
                    apart = std::max (apart, (double) std::fabs (a[c][(std::size_t) i] - late[c][(std::size_t) i]));
                }
            test::ok (same, "two histories, one fade position, identical bits" + cell + " (worst " + sci (worst) + ")");
            for (int c = 0; c < nch; ++c)
            {
                const std::vector<float>& p = prog[c];
                const std::vector<float> yB = convRef (p, oB[c], irMaxC);
                const double tail = maxDiff (a[c], yB, xfC + 2000, progC) / (peak (yB) + 1e-30);
                test::ok (tail < 2e-4, "the surviving fade completes onto the published operator, channel " + std::to_string (c) + cell + " (" + sci (tail) + ")");
            }
            test::ok (apart > 1e-3, "…but the answer depends on WHERE in the fade it was called — reset() is the verb that does not" + cell + " (" + sci (apart) + ")");
        }
    }

    // --- in-place aliasing (in == out) matches out-of-place for every topology (both go through the same fade) ---
    test::group ("in-place aliasing (in == out): mono / LRDiag / MSDiag / Full");
    {
        const int L = 3000, N = L + 6000;
        Lcg r { 4040 };
        std::vector<float> b0 ((std::size_t) L), b1 ((std::size_t) L), b2 ((std::size_t) L), b3 ((std::size_t) L);
        for (auto& v : b0) v = 0.05f * r.next(); for (auto& v : b1) v = 0.03f * r.next(); for (auto& v : b2) v = 0.03f * r.next(); for (auto& v : b3) v = 0.05f * r.next();
        std::vector<float> xL ((std::size_t) N), xR ((std::size_t) N); for (auto& v : xL) v = 0.3f * r.next(); for (auto& v : xR) v = 0.4f * r.next();
        for (auto topo : { MCN::Topology::LRDiag, MCN::Topology::MSDiag, MCN::Topology::Full })
        {
            const int nb = MCN::numBanksFor (topo);
            const float* bk[4] { b0.data(), b1.data(), b2.data(), b3.data() };
            MCN a; a.prepare (128, maxIr, 128, 2); a.setOperator (topo, bk, nb, L);
            std::vector<float> yL, yR; runStereo (a, xL, xR, yL, yR, 128);
            MCN b; b.prepare (128, maxIr, 128, 2); b.setOperator (topo, bk, nb, L);
            std::vector<float> iL = xL, iR = xR;
            for (int i = 0; i < N; i += 128) { const int m = std::min (128, N - i); const float* in[2] { &iL[(std::size_t) i], &iR[(std::size_t) i] }; float* out[2] { &iL[(std::size_t) i], &iR[(std::size_t) i] }; felitronics::test::run (b.process (in, out, 2, m)); }
            test::ok (maxDiff (iL, yL, 0, N) < 1e-6 && maxDiff (iR, yR, 0, N) < 1e-6, std::string ("in-place bit-matches out-of-place"));
        }
    }

    // --- no allocation in process(), incl. across a crossfade (mono + stereo Full) ---
    test::group ("no allocation in process() (incl. across a crossfade)");
    {
        const int L = 4000;
        std::vector<float> h1 ((std::size_t) L, 0.001f), h2 ((std::size_t) L, -0.001f);
        MCN mc; mc.prepare (128, maxIr, 128, 1); mc.setIr (h1.data(), L);
        std::vector<float> x (2048, 0.2f), y (2048, 0.0f);
        const float* in[1] { x.data() }; float* out[1] { y.data() };
        felitronics::test::run (mc.process (in, out, 1, 2048));
        mc.setIr (h2.data(), L);                                     // stage a swap → next process crossfades
        const long long before = alloc::count.load();
        felitronics::test::run (mc.process (in, out, 1, 2048));                              // inside the crossfade (blends both slots)
        felitronics::test::run (mc.process (in, out, 1, 2048));
        test::okNoAlloc (alloc::count.load() == before, "mono process() zero heap allocations across a crossfade");

        std::vector<float> f0 ((std::size_t) L, 0.001f), f1 ((std::size_t) L, 0.0005f), f2 ((std::size_t) L, 0.0005f), f3 ((std::size_t) L, -0.001f);
        MCN mf; mf.prepare (128, maxIr, 128, 2);
        { const float* bk[4] { f0.data(), f1.data(), f2.data(), f3.data() }; mf.setOperator (MCN::Topology::Full, bk, 4, L); }
        std::vector<float> xl (2048, 0.2f), xr (2048, -0.1f); const float* sin[2] { xl.data(), xr.data() }; float* sout[2] { xl.data(), xr.data() };
        felitronics::test::run (mf.process (sin, sout, 2, 2048));
        const long long before2 = alloc::count.load();
        felitronics::test::run (mf.process (sin, sout, 2, 2048));
        test::okNoAlloc (alloc::count.load() == before2, "stereo Full process() zero heap allocations");
    }

    // --- API surface + prepare guards ---
    test::group ("API parity + prepare guards");
    {
        test::ok (MCN::latencySamples() == 0, "latencySamples()==0");
        test::ok (MCN::numBanksFor (MCN::Topology::Full) == 4 && MCN::numBanksFor (MCN::Topology::LRDiag) == 2, "numBanksFor parity");
        MCN mc;
        const float* b[2] { nullptr, nullptr };
        test::ok (! mc.setOperator (MCN::Topology::LRDiag, b, 2, 100), "setOperator before prepare rejected");
        test::ok (! mc.prepare (127, maxIr, 128, 2), "non-pow2 head rejected");
        test::ok (! mc.prepare (128, maxIr, 128, 3), "numChannels=3 rejected");
        test::ok (  mc.prepare (128, maxIr, 128, 2), "prepare after rejected calls");
        test::ok (! mc.isBusy(), "idle after prepare");
        std::vector<float> ir (200, 0.01f); const float* one[2] { ir.data(), ir.data() };
        test::ok (mc.setOperator (MCN::Topology::LRDiag, one, 2, 200), "setOperator accepted (idle)");
        std::vector<float> l (8, 0.1f), rr (8, 0.1f); const float* in[2] { l.data(), rr.data() }; float* out[2] { l.data(), rr.data() };
        felitronics::test::run (mc.process (in, out, 2, 8));                                  // begins the cold fade → busy
        test::ok (mc.isBusy(), "busy during the cold prime crossfade");
        const float* b2[2] { ir.data(), ir.data() };
        test::ok (! mc.setOperator (MCN::Topology::LRDiag, b2, 2, 200), "second operator rejected while fading (host coalesces)");
    }

    return test::report();
}

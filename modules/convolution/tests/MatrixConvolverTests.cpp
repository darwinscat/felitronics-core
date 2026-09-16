// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// JUCE-free self-tests for the MatrixConvolver: the 2×2 OPERATOR convolver on one canonical raw-L/R
// history. Reference-NULL against the proven PartitionedConvolver: each routing topology is nulled against
// an INDEPENDENT time-domain computation of the same math, so the re-plumb (raw-L/R FDL + on-the-fly M/S
// views + cross-input sums) is proven, not assumed. Plus: an atomic topology SWITCH mid-playback proven
// click-free AND warm (its incoming operator carries the input's past, unlike a cold-started instance),
// and no allocation in process() even across a swap.

#include <felitronics_test.h>
#include <alloc_counter.h>   // installs the allocation counter: EVERY form of `new`, over-aligned included
#include <felitronics/convolution/MatrixConvolver.h>
#include <felitronics/convolution/PartitionedConvolver.h>   // reference convolver for the null tests
#include <felitronics/core/Math.h>   // core::sameBits — a bitwise compare, not ==

#include <cmath>
#include <cstdio>
#include <cstdlib>
#if defined(_WIN32)
 #include <malloc.h>   // _aligned_malloc / _aligned_free (MSVC has no posix_memalign)
#endif
#include <string>
#include <utility>
#include <vector>

using namespace felitronics;
using MC = convolution::MatrixConvolver<>;

struct Lcg { unsigned long long s; float next() { s = s * 6364136223846793005ULL + 1442695040888963407ULL; return (float) ((s >> 40) & 0xffff) / 32768.0f - 1.0f; } };

static std::string sci (double v) { char b[32]; std::snprintf (b, sizeof b, "%.3e", v); return b; }

// One mono reference convolution (in ∗ ir) via the proven PartitionedConvolver, from zero history.
static std::vector<float> convRef (const std::vector<float>& ir, const std::vector<float>& in, int P, int irMax)
{
    convolution::PartitionedConvolver<> pc; pc.prepare (P, irMax); pc.setIr (ir.data(), (int) ir.size());
    std::vector<float> out (in.size(), 0.0f);
    felitronics::test::run (pc.process (in.data(), out.data(), (int) in.size()));
    return out;
}

static double maxDiff (const std::vector<float>& a, const std::vector<float>& b, int from, int to)
{
    double m = 0.0; for (int i = from; i < to; ++i) m = std::max (m, (double) std::fabs (a[(std::size_t) i] - b[(std::size_t) i])); return m;
}
static double maxDeriv (const std::vector<float>& y, int from, int to)
{
    double m = 0.0; for (int i = from + 1; i < to; ++i) m = std::max (m, (double) std::fabs (y[(std::size_t) i] - y[(std::size_t) (i - 1)])); return m;
}

int main()
{
    std::printf ("felitronics::convolution MatrixConvolver tests\n");
    const int P = 64, irMax = 400, len = 200, xfade = 128, n = 4000;
    const int settled = 900;                              // past the short crossfade + tail fill + margin

    Lcg r { 2718 };
    auto mkIr = [&] { std::vector<float> v ((std::size_t) len); for (auto& x : v) x = 0.15f * r.next(); return v; };
    std::vector<float> irM = mkIr(), irS = mkIr(), irL = mkIr(), irR = mkIr();
    std::vector<float> irLL = mkIr(), irLR = mkIr(), irRL = mkIr(), irRR = mkIr();

    std::vector<float> xL (n), xR (n);
    for (int i = 0; i < n; ++i) { xL[(std::size_t) i] = (float) (0.5 * std::sin (2.0 * core::kPi * 500.0 * i / 48000.0)); xR[(std::size_t) i] = 0.4f * r.next(); }

    auto runStereo = [&] (MC& mc, std::vector<float>& oL, std::vector<float>& oR)
    {
        oL.assign ((std::size_t) n, 0.0f); oR.assign ((std::size_t) n, 0.0f);
        for (int i = 0; i < n; ++i) { oL[(std::size_t) i] = xL[(std::size_t) i]; oR[(std::size_t) i] = xR[(std::size_t) i]; }
        for (int o = 0; o < n; o += 512) { float* io[2] { oL.data() + o, oR.data() + o }; felitronics::test::run (mc.process (io, io, 2, std::min (512, n - o))); }
    };

    // --- LRDiag == two independent mono convolutions (direct L/R routing) ---
    test::group ("MatrixConvolver LRDiag == two independent mono convolutions");
    {
        MC mc; test::ok (mc.prepare (P, irMax, xfade, 2), "prepare stereo");
        const float* banks[2] { irL.data(), irR.data() };
        test::ok (mc.setOperator (MC::Topology::LRDiag, banks, 2, len), "setOperator LRDiag");
        std::vector<float> oL, oR; runStereo (mc, oL, oR);
        const std::vector<float> refL = convRef (irL, xL, P, irMax), refR = convRef (irR, xR, P, irMax);
        test::ok (maxDiff (oL, refL, settled, n) < 1e-4, "yL == irL ∗ xL (independent of xR)");
        test::ok (maxDiff (oR, refR, settled, n) < 1e-4, "yR == irR ∗ xR (independent of xL)");
    }

    // --- MSDiag == time-domain encode→conv→decode (the spectral-view re-plumb proof) ---
    test::group ("MatrixConvolver MSDiag == encode→conv→decode reference");
    {
        MC mc; mc.prepare (P, irMax, xfade, 2);
        const float* banks[2] { irM.data(), irS.data() };
        test::ok (mc.setOperator (MC::Topology::MSDiag, banks, 2, len), "setOperator MSDiag");
        std::vector<float> oL, oR; runStereo (mc, oL, oR);
        std::vector<float> m (n), s (n);
        for (int i = 0; i < n; ++i) { m[(std::size_t) i] = 0.5f * (xL[(std::size_t) i] + xR[(std::size_t) i]); s[(std::size_t) i] = 0.5f * (xL[(std::size_t) i] - xR[(std::size_t) i]); }
        const std::vector<float> yM = convRef (irM, m, P, irMax), yS = convRef (irS, s, P, irMax);
        std::vector<float> refL (n), refR (n);
        for (int i = 0; i < n; ++i) { refL[(std::size_t) i] = yM[(std::size_t) i] + yS[(std::size_t) i]; refR[(std::size_t) i] = yM[(std::size_t) i] - yS[(std::size_t) i]; }
        const double eL = maxDiff (oL, refL, settled, n), eR = maxDiff (oR, refR, settled, n);
        std::printf ("      MSDiag view-vs-time-encode max err L=%.2e R=%.2e\n", eL, eR);
        test::ok (eL < 1e-5, "yL == (yM+yS) — M/S spectral view == time-domain encode (float-exact)");
        test::ok (eR < 1e-5, "yR == (yM−yS) — the on-the-fly view re-plumb is exact");
    }

    // --- Full == direct 4-conv cross sums (cross-input routing + polarity) ---
    test::group ("MatrixConvolver Full == 4-conv cross-sum reference");
    {
        MC mc; mc.prepare (P, irMax, xfade, 2);
        const float* banks[4] { irLL.data(), irLR.data(), irRL.data(), irRR.data() };
        test::ok (mc.setOperator (MC::Topology::Full, banks, 4, len), "setOperator Full (4 banks)");
        std::vector<float> oL, oR; runStereo (mc, oL, oR);
        const std::vector<float> cLL = convRef (irLL, xL, P, irMax), cLR = convRef (irLR, xR, P, irMax);
        const std::vector<float> cRL = convRef (irRL, xL, P, irMax), cRR = convRef (irRR, xR, P, irMax);
        std::vector<float> refL (n), refR (n);
        for (int i = 0; i < n; ++i) { refL[(std::size_t) i] = cLL[(std::size_t) i] + cLR[(std::size_t) i]; refR[(std::size_t) i] = cRL[(std::size_t) i] + cRR[(std::size_t) i]; }
        test::ok (maxDiff (oL, refL, settled, n) < 5e-4, "yL == LL∗xL + LR∗xR");
        test::ok (maxDiff (oR, refR, settled, n) < 5e-4, "yR == RL∗xL + RR∗xR");
        // off-diagonal actually routes: dropping xR must change yL (LR∗xR term is real)
        double lrEnergy = 0.0; for (int i = settled; i < n; ++i) lrEnergy += std::fabs (cLR[(std::size_t) i]);
        test::ok (lrEnergy > 1e-2, "the LR cross term is non-trivial (Full genuinely mixes channels)");
    }

    // --- topology SWITCH mid-playback (2-bank MSDiag → 4-bank Full) : click-free + WARM history ---
    test::group ("MatrixConvolver topology switch MSDiag→Full: click-free + warm history");
    {
        const int T = 2000;                              // switch point (well past the cold prime → warm)
        // A: warm switch. Run MSDiag to T, then switch to Full mid-stream, keep streaming.
        MC A; A.prepare (P, irMax, xfade, 2);
        { const float* b[2] { irM.data(), irS.data() }; A.setOperator (MC::Topology::MSDiag, b, 2, len); }
        std::vector<float> aL (n), aR (n);
        for (int i = 0; i < n; ++i) { aL[(std::size_t) i] = xL[(std::size_t) i]; aR[(std::size_t) i] = xR[(std::size_t) i]; }
        bool switched = false;
        for (int o = 0; o < n; o += 256)
        {
            if (! switched && o >= T)
            { const float* b[4] { irLL.data(), irLR.data(), irRL.data(), irRR.data() };
              if (A.setOperator (MC::Topology::Full, b, 4, len)) switched = true; }
            float* io[2] { aL.data() + o, aR.data() + o }; felitronics::test::run (A.process (io, io, 2, std::min (256, n - o)));
        }
        test::ok (switched, "topology switch accepted mid-stream");

        // click-free: the derivative across the switch stays near the steady-state derivative.
        const double steady = maxDeriv (aL, T - 400, T - 50);
        const double across = maxDeriv (aL, T, T + xfade + 300);
        test::ok (across < 4.0 * steady + 1e-6, "no click across the topology change");

        // correctness: after the fade, A == the Full operator fed the WHOLE stream (warm history honoured).
        const std::vector<float> cLL = convRef (irLL, xL, P, irMax), cLR = convRef (irLR, xR, P, irMax);
        const std::vector<float> cRL = convRef (irRL, xL, P, irMax), cRR = convRef (irRR, xR, P, irMax);
        std::vector<float> fullL (n), fullR (n);
        for (int i = 0; i < n; ++i) { fullL[(std::size_t) i] = cLL[(std::size_t) i] + cLR[(std::size_t) i]; fullR[(std::size_t) i] = cRL[(std::size_t) i] + cRR[(std::size_t) i]; }
        const int after = T + xfade + P + 8;
        test::ok (maxDiff (aL, fullL, after, n) < 5e-4, "A converges to the true Full response (switched correctly)");

        // WARM proof: a COLD instance fed ONLY the post-switch samples differs from A right after the switch —
        // A's incoming Full operator reads the shared warm FDL (the input's past), the cold one cannot.
        MC B; B.prepare (P, irMax, xfade, 2);
        { const float* b[4] { irLL.data(), irLR.data(), irRL.data(), irRR.data() }; B.setOperator (MC::Topology::Full, b, 4, len); }
        const int tail = n - T;
        std::vector<float> bL (tail), bR (tail);
        for (int i = 0; i < tail; ++i) { bL[(std::size_t) i] = xL[(std::size_t) (T + i)]; bR[(std::size_t) i] = xR[(std::size_t) (T + i)]; }
        for (int o = 0; o < tail; o += 256) { float* io[2] { bL.data() + o, bR.data() + o }; felitronics::test::run (B.process (io, io, 2, std::min (256, tail - o))); }
        // compare A[T + after..] to B[after..] over a window where the pre-switch tail still matters
        double warmDiff = 0.0;
        for (int i = xfade + P; i < xfade + P + 400; ++i)
            warmDiff = std::max (warmDiff, (double) std::fabs (aL[(std::size_t) (T + i)] - bL[(std::size_t) i]));
        test::ok (warmDiff > 1e-2, "A (warm) ≠ a cold instance fed only post-switch samples → shared warm FDL");
    }

    // --- MONO degenerate == a single PartitionedConvolver ---
    test::group ("MatrixConvolver mono == one PartitionedConvolver");
    {
        MC mc; test::ok (mc.prepare (P, irMax, xfade, 1), "prepare mono"); test::ok (mc.numChannels() == 1, "1 channel");
        test::ok (mc.setIr (irM.data(), len), "setIr (mono convenience)");
        std::vector<float> y (n); for (int i = 0; i < n; ++i) y[(std::size_t) i] = xL[(std::size_t) i];
        for (int o = 0; o < n; o += 512) { float* io[1] { y.data() + o }; felitronics::test::run (mc.process (io, io, 1, std::min (512, n - o))); }
        const std::vector<float> ref = convRef (irM, xL, P, irMax);
        test::ok (maxDiff (y, ref, settled, n) < 1e-4, "mono output == irM ∗ xL");
    }

    // --- no allocation in process(), even across a topology-changing swap ---
    test::group ("MatrixConvolver no-alloc in process() (incl. across a swap)");
    {
        MC mc; mc.prepare (P, irMax, xfade, 2);
        { const float* b[2] { irM.data(), irS.data() }; mc.setOperator (MC::Topology::MSDiag, b, 2, len); }
        std::vector<float> l (512, 0.2f), rr (512, -0.1f); float* io[2] { l.data(), rr.data() };
        felitronics::test::run (mc.process (io, io, 2, 512));                      // consume the initial fade-in
        { const float* b[4] { irLL.data(), irLR.data(), irRL.data(), irRR.data() }; mc.setOperator (MC::Topology::Full, b, 4, len); }
        const long long before = alloc::count.load();
        felitronics::test::run (mc.process (io, io, 2, 512));                      // crosses the crossfade (MSDiag→Full)
        felitronics::test::run (mc.process (io, io, 2, 512));
        test::okNoAlloc (alloc::count.load() == before, "process() performed zero heap allocations across a topology swap");
    }

    // --- 🔴 P88: reset() KEEPS AN OPERATOR PUBLISHED BUT NOT YET ADOPTED (the third instance of the class) ---
    // Same defect and same fix as ConvolutionEngine::reset(), which carries the full reasoning: `cur_` flips
    // only at fade end, so a reset() that wiped `state_` and kept `cur_` put the PREVIOUS operator back and
    // nothing re-staged the published one. On that body this group reads the old operator in every cell,
    // a worst sample 1.017e+00 from the settled restart.
    // The oracle is law 11a INDEPENDENCE — two convolvers fed DIFFERENT audio of DIFFERENT length, one
    // restarted with the publication in flight and one after it settled, answer the next programme with
    // identical bits — for every TOPOLOGY, with four distinct banks and different L/R inputs (with L == R
    // an MSDiag side is silent and a routing slip after the adoption is invisible — a review round showed
    // `slot_[cur_].topo = LRDiag` surviving the first version of this group), and at BOTH slot parities,
    // since `cur_ = 0` is right half the time. The LRDiag cells are certified from OUTSIDE the class too.
    test::group ("MatrixConvolver reset keeps a published-but-unadopted operator (independence, every topology)");
    {
        const int Pk = 64, irMaxK = 1024, irLenK = 700, xfK = 128;
        const int settleK = xfK + Pk + 8;
        const int progN   = xfK + 17 * Pk;
        Lcg rk { 4242 };
        // operators 0 = A, 1 = A2 (the parity pre-swap), 2 = B (published), 3 = C (published after); 4 banks each
        std::vector<float> ops[4][4];
        for (int o = 0; o < 4; ++o)
            for (int b = 0; b < 4; ++b)
            {
                ops[o][b].assign ((std::size_t) irLenK, 0.0f);
                for (auto& v : ops[o][b]) v = 0.05f * rk.next();
                ops[o][b][0] = 0.1f * (float) (1 + o) * (float) (1 + b) * (((o + b) & 1) != 0 ? -1.0f : 1.0f);
            }
        const std::size_t preDN = (std::size_t) (((irMaxK - Pk + Pk - 1) / Pk) * Pk + 700);
        std::vector<float> preD[2] { std::vector<float> (preDN), std::vector<float> (preDN) };
        const std::size_t preRN = (std::size_t) (((irMaxK - Pk + Pk - 1) / Pk) * Pk + 1100);
        std::vector<float> preR[2] { std::vector<float> (preRN), std::vector<float> (preRN) };
        std::vector<float> prog[2] { std::vector<float> ((std::size_t) progN), std::vector<float> ((std::size_t) progN) };
        for (auto& ch : preD) for (auto& v : ch) v = 0.5f * rk.next();
        for (auto& ch : preR) for (auto& v : ch) v = 0.4f * rk.next();
        for (auto& ch : prog) for (auto& v : ch) v = 0.3f * rk.next();

        // mode 0 = restart at Pending · 1 = restart mid-crossfade · 2 = the reference (publication settled)
        const auto render = [&] (int nch, MC::Topology topo, int parity, int mode,
                                 std::vector<float>& L, std::vector<float>& R, bool* busyOut, bool* acceptOut)
        {
            MC mc;
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
                return mc.setOperator (topo, b, nch == 1 ? 1 : MC::numBanksFor (topo), irLenK);
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

        using Topo = MC::Topology;
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
                                const std::vector<float> yB = convRef (ops[2][c], prog[c], Pk, irMaxK), yA = convRef (ops[0][c], prog[c], Pk, irMaxK);
                                const double eB = maxDiff (dut[c], yB, 0, progN), eA = maxDiff (dut[c], yA, 0, progN);
                                test::ok (eB < 1e-6, "  …and channel " + std::to_string (c) + " answers as the PUBLISHED operator" + cell + " (" + sci (eB) + ")");
                                test::ok (eA > 1e-2,  "  …and not as the previous one" + cell + " (" + sci (eA) + ")");
                            }
                    }
                }
            }

        // A STAGED but NOT PUBLISHED operator (state_ == 0) is none of reset()'s business: `cur_` must not
        // move, the banks must survive, and the later publishStaged() must still find it and fade it in.
        for (int nch = 1; nch <= 2; ++nch)
        {
            MC mc;
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
                const std::vector<float> yB = convRef (ops[2][c], prog[c], Pk, irMaxK);
                const double tail = maxDiff (y[c], yB, xfK + 2 * Pk, progN);
                test::ok (tail < 1e-6, "the publish that follows the restart still lands on the staged operator  ["
                                          + std::to_string (nch) + " ch, channel " + std::to_string (c) + "] (" + sci (tail) + ")");
            }
        }
    }

    // --- 🔴 P88: clearAudioState() — the history alone, the swap left running FROM WHERE IT WAS ---
    // Two convolvers with different pasts, cleared at the same fade position, answer the same bits; the fade
    // then ENDS exactly when it would have (a review round showed `xfadePos_ = 0` inside this verb surviving
    // the first version of the group), lands on the published operator, and — the one thing this verb does
    // not give — the answer still depends on WHERE in the fade it was called. That is reset()'s promise.
    test::group ("MatrixConvolver clearAudioState leaves the swap running and the history gone");
    {
        const int Pc = 64, irMaxC = 1024, irLenC = 700, xfC = 512;
        const int progC = xfC + 17 * Pc;
        Lcg rc { 8181 };
        std::vector<float> oA[2], oB[2];
        for (int b = 0; b < 2; ++b)
        {
            oA[b].assign ((std::size_t) irLenC, 0.0f); for (auto& v : oA[b]) v = 0.05f * rc.next();
            oB[b].assign ((std::size_t) irLenC, 0.0f); for (auto& v : oB[b]) v = 0.05f * rc.next();
        }
        const std::size_t preDN = (std::size_t) (((irMaxC - Pc + Pc - 1) / Pc) * Pc + 500);
        std::vector<float> preD[2] { std::vector<float> (preDN), std::vector<float> (preDN) };
        const std::size_t preRN = (std::size_t) (((irMaxC - Pc + Pc - 1) / Pc) * Pc + 900);
        std::vector<float> preR[2] { std::vector<float> (preRN), std::vector<float> (preRN) };
        std::vector<float> prog[2] { std::vector<float> ((std::size_t) progC), std::vector<float> ((std::size_t) progC) };
        for (auto& ch : preD) for (auto& v : ch) v = 0.5f * rc.next();
        for (auto& ch : preR) for (auto& v : ch) v = 0.4f * rc.next();
        for (auto& ch : prog) for (auto& v : ch) v = 0.3f * rc.next();

        // renders `prog` after a clear `intoFade` samples into the fade; reports busy just before and just
        // after the sample on which the fade must end
        const auto render = [&] (int nch, const std::vector<float>* pre, int intoFade, std::vector<float>* out,
                                 bool* busyAfterClear, bool* busyOneShort, bool* busyOnTime)
        {
            MC mc;
            felitronics::test::run (mc.prepare (Pc, irMaxC, xfC, nch));
            const int preN = (int) pre[0].size();
            std::vector<float> jl ((std::size_t) preN, 0.0f), jr ((std::size_t) preN, 0.0f);
            const auto feed = [&] (int from, int k)
            {
                const float* in[2]  { pre[0].data() + from, pre[1].data() + from };
                float*       o[2]   { jl.data(), jr.data() };
                felitronics::test::run (mc.process (in, o, nch, k));
            };
            { const float* b[2] { oA[0].data(), oA[1].data() }; test::ok (mc.setOperator (MC::Topology::LRDiag, b, nch == 1 ? 1 : 2, irLenC), "A accepted"); }
            feed (0, preN - xfC);
            { const float* b[2] { oB[0].data(), oB[1].data() }; test::ok (mc.setOperator (MC::Topology::LRDiag, b, nch == 1 ? 1 : 2, irLenC), "B published"); }
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
                const std::vector<float> yB = convRef (oB[c], p, Pc, irMaxC);
                const double tail = maxDiff (a[c], yB, xfC + 2 * Pc, progC);
                test::ok (tail < 1e-6, "the surviving fade completes onto the published operator, channel " + std::to_string (c) + cell + " (" + sci (tail) + ")");
            }
            test::ok (apart > 1e-3, "…but the answer depends on WHERE in the fade it was called — reset() is the verb that does not" + cell + " (" + sci (apart) + ")");
        }
    }

    // --- swap coalescing + unprepared guards ---
    test::group ("MatrixConvolver swap coalescing + unprepared guards");
    {
        MC fresh;
        const float* b[2] { irM.data(), irS.data() };
        test::ok (! fresh.setOperator (MC::Topology::MSDiag, b, 2, len), "setOperator before prepare() rejected");
        float l[16] {}, rr[16] {}; float* io[2] { l, rr };
        test::ok (! fresh.process (io, io, 2, 16), "process() before prepare() is REFUSED (law 11)");
        test::ok (fresh.prepare (P, irMax, xfade, 2), "prepare after the rejected call");
        test::ok (fresh.setOperator (MC::Topology::MSDiag, b, 2, len), "first operator accepted (idle)");
        std::vector<float> l2 (512, 0.1f), r2 (512, 0.1f); float* io2[2] { l2.data(), r2.data() };
        felitronics::test::run (fresh.process (io2, io2, 2, 8));                  // begin the (cold) fade → busy
        test::ok (fresh.isBusy(), "busy during the cold prime");
        const float* b2[4] { irLL.data(), irLR.data(), irRL.data(), irRR.data() };
        test::ok (! fresh.setOperator (MC::Topology::Full, b2, 4, len), "second operator rejected while busy (host coalesces)");
    }

    return test::report();
}

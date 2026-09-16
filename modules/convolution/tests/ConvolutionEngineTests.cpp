// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// JUCE-free self-tests for the production ConvolutionEngine: a CLICK-FREE IR swap (no derivative spike)
// that CONVERGES to the new IR's response, zero latency, and no-allocation-in-process().

#include <felitronics_test.h>
#include <alloc_counter.h>   // installs the allocation counter: EVERY form of `new`, over-aligned included
#include <felitronics/convolution/ConvolutionEngine.h>
#include <felitronics/convolution/PartitionedConvolver.h>   // reference convolver for the null tests
#include <felitronics/core/Math.h>                         // core::sameBits — a bitwise compare, not ==

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <utility>
#if defined(_WIN32)
 #include <malloc.h>   // _aligned_malloc / _aligned_free (MSVC has no posix_memalign)
#endif
#include <vector>

using namespace felitronics;

struct Lcg { unsigned long long s; float next() { s = s * 6364136223846793005ULL + 1442695040888963407ULL; return (float) ((s >> 40) & 0xffff) / 32768.0f - 1.0f; } };

static std::string sci (double v) { char b[32]; std::snprintf (b, sizeof b, "%.3e", v); return b; }

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
        felitronics::test::run (eng.process (x.data(), y.data(), 1400));        // settle on IR A
        const bool swapOk = eng.setIr (irB.data(), irLen);
        felitronics::test::run (eng.process (x.data() + 1400, y.data() + 1400, n - 1400));   // crossfade A→B + settle
        test::ok (swapOk, "swap accepted while idle");

        // continuity: the derivative across the swap must not spike vs the steady region (no click).
        const double steady = maxDeriv (y, 900, 1300);
        const double swapRgn = maxDeriv (y, 1400, 1400 + xfade + 100);
        test::ok (swapRgn < 3.0 * steady + 1e-6, "no click — swap derivative ~ steady derivative");

        // convergence: well after the swap, output == a fresh IR-B convolver fed the same input.
        convolution::PartitionedConvolver<> ref; ref.prepare (P, irMax); ref.setIr (irB.data(), irLen);
        std::vector<float> yref (n, 0.0f); felitronics::test::run (ref.process (x.data(), yref.data(), n));
        double maxErr = 0.0; for (int i = 2200; i < n; ++i) maxErr = std::max (maxErr, (double) std::fabs (y[i] - yref[i]));
        test::ok (maxErr < 3e-3, "converges to the new IR's steady response");
    }

    // --- no allocation during process() (incl. crossing a swap) ---
    test::group ("ConvolutionEngine no-alloc in process()");
    {
        convolution::ConvolutionEngine<> eng; eng.prepare (P, irMax, xfade);
        std::vector<float> in (512, 0.2f), out (512, 0.0f);
        eng.setIr (irA.data(), irLen);
        felitronics::test::run (eng.process (in.data(), out.data(), 512));      // consume the initial fade-in
        eng.setIr (irB.data(), irLen);                 // arm a swap (build is message-thread, before the snapshot)
        const long long before = alloc::count.load();
        felitronics::test::run (eng.process (in.data(), out.data(), 512));      // crosses the crossfade
        felitronics::test::run (eng.process (in.data(), out.data(), 512));
        const long long after = alloc::count.load();
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
            felitronics::test::run (eng.process (in, out, 2, 1400));                        // settle on A
        }
        eng.setIr (irB.data(), irLen);                             // arm swap on both channels at once
        {
            const float* in[2]  { x.data() + 1400, x.data() + 1400 };
            float*       out[2] { l.data() + 1400, rr.data() + 1400 };
            felitronics::test::run (eng.process (in, out, 2, n - 1400));                    // crossfade A→B
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
        st.setIr (irA.data(), irLen); felitronics::test::run (st.process (sin, sout, 2, n));
        mo.setIr (irA.data(), irLen); felitronics::test::run (mo.process (x.data(), mout.data(), n));
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
        felitronics::test::run (eng.process (in, out, 2, n));

        convolution::PartitionedConvolver<> rA, rB; rA.prepare (P, irMax); rB.prepare (P, irMax);
        rA.setIr (irA.data(), irLen); rB.setIr (irB.data(), irLen);
        std::vector<float> yA (n, 0.0f), yB (n, 0.0f); felitronics::test::run (rA.process (x.data(), yA.data(), n)); felitronics::test::run (rB.process (x.data(), yB.data(), n));
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
        eng.setIr (iA.data(), irLenW); felitronics::test::run (eng.process (xw.data(), yw.data(), sw));
        const bool warmOk = eng.setIr (iB.data(), irLenW);
        felitronics::test::run (eng.process (xw.data() + sw, yw.data() + sw, nw - sw));
        test::ok (warmOk, "warm swap accepted");

        convolution::PartitionedConvolver<> refB; refB.prepare (Pw, irMaxW); refB.setIr (iB.data(), irLenW);
        std::vector<float> yB (nw, 0.0f); felitronics::test::run (refB.process (xw.data(), yB.data(), nw));
        const int settled = sw + warmXf + Pw + 4;                  // just past the SHORT fade + one chunk
        double e = 0.0; for (int i = settled; i < nw; ++i) e = std::max (e, (double) std::fabs (yw[i] - yB[i]));
        test::ok (e < 2e-3, "matches the new IR right after the short fade — no N-length lag");
    }

    // --- every activation (the first prime AND a warm swap) uses the SAME short fade — no long cold prime ---
    // (a cold-started FDL already yields the exact causal convolution, so there is nothing to mask by attenuating)
    test::group ("ConvolutionEngine first activation uses the short fade (like a warm swap)");
    {
        const int Pc = 64, irMaxC = 1024, irLenC = 700, warmXf = 128;
        Lcg rc { 99 };
        std::vector<float> iA (irLenC), iB (irLenC), xc (4000);
        for (auto& v : iA) v = 0.05f * rc.next();
        for (auto& v : iB) v = 0.05f * rc.next();
        for (auto& v : xc) v = 0.5f  * rc.next();

        convolution::ConvolutionEngine<> eng; eng.prepare (Pc, irMaxC, warmXf);
        std::vector<float> y (4000, 0.0f); int pos = 0;
        auto run = [&] (int k) { felitronics::test::run (eng.process (xc.data() + pos, y.data() + pos, k)); pos += k; };
        eng.setIr (iA.data(), irLenC);
        run (warmXf / 2);      test::ok (eng.isBusy(),   "first activation busy mid-fade");
        run (warmXf + Pc + 8); test::ok (! eng.isBusy(), "first activation finished within the short fade (no long cold prime)");
        test::ok (eng.setIr (iB.data(), irLenC), "warm swap accepted");
        run (warmXf + Pc + 8); test::ok (! eng.isBusy(), "warm swap finished within the same short fade");
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
        int pos = 0; auto run = [&] (int k) { felitronics::test::run (eng.process (ones.data() + pos, y.data() + pos, k)); pos += k; };
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

    // --- reset() flushes the running history but KEEPS the live IR; the next swap uses the short fade ---
    test::group ("ConvolutionEngine reset keeps the live IR, next swap short");
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
        auto run = [&] (int k) { felitronics::test::run (eng.process (xr.data() + pos, y.data() + pos, k)); pos += k; };
        eng.setIr (iA.data(), irLenR); run (coldXf + 200);          // run A a while to warm the shared history
        test::ok (eng.setIr (iB.data(), irLenR), "swap to B accepted (warm)");
        run (warmXf + Pr + 8);
        test::ok (! eng.isBusy(), "warm swap to B done — B is the live slot (cur_ flipped)");

        eng.reset();                                                // flush history; B must stay live
        float imp[1] { 1.0f }, oimp[1] { 0.0f };
        felitronics::test::run (eng.process (imp, oimp, 1));                                 // first post-reset sample = head[live]·δ = iB[0]
        test::ok (std::fabs (oimp[0] - iB[0]) < 1e-6, "reset kept the live IR (B's head), did not revert to A");
        test::ok (std::fabs (oimp[0] - iA[0]) > 1e-3, "  …and it is NOT A");
        test::ok (eng.setIr (iC.data(), irLenR), "post-reset swap accepted");
        run (warmXf + Pr + 8);
        test::ok (! eng.isBusy(), "post-reset swap settles within the SHORT fade (history was flushed, not re-armed long)");
    }

    // --- 🔴 P88: reset() KEEPS AN OPERATOR PUBLISHED BUT NOT YET ADOPTED — the property, both parities ---
    // reset() promises "flush the tail", not "revert the EQ", and that promise covers a publication the
    // audio thread has not picked up yet (law 11e). The body before P88 kept `cur_` on the OLD slot and
    // wiped `state_`, so a published operator was dropped with nothing left to re-stage it — setIr() had
    // already said true. On that body this group reads the OLD operator in every cell, a worst sample
    // 7.494e-01 from the settled restart. The group above tests the promise only from Idle,
    // the one state where it cannot fail.
    // THE ORACLE IS LAW 11a INDEPENDENCE, not a head tap: two engines fed DIFFERENT audio for DIFFERENT
    // lengths — so `phase_` and `fdlPos_` differ too — one restarted with the publication still in flight
    // and one restarted after it settled, must answer the next programme with IDENTICAL BITS. A head-tap
    // probe alone passes a restart that forgets the FDL or the frame; this does not, and the programme is
    // longer than the whole FDL span so a single stale partition shows. The channels get DIFFERENT inputs,
    // so a restart that mixed them up would show too. Both slot parities run, because `cur_ = 0` instead
    // of `cur_ = 1 - cur_` is right half the time — the mutation stand catches that mutant only at parity 1.
    // Every cell is also certified from OUTSIDE the class against a PartitionedConvolver holding the
    // PUBLISHED IR, and checked NOT to be the previous one.
    test::group ("ConvolutionEngine reset keeps a published-but-unadopted IR (law 11a independence)");
    {
        const int Pk = 64, irMaxK = 1024, irLenK = 700, xfK = 128;
        const int coldK   = ((irMaxK - Pk + Pk - 1) / Pk) * Pk;   // maxParts·P — the whole FDL span
        const int settleK = xfK + Pk + 8;
        const int progN   = xfK + 17 * Pk;                        // past the FDL span: a stale partition shows
        Lcg rk { 1234 };
        std::vector<float> iA (irLenK), iA2 (irLenK), iB (irLenK), iC (irLenK);
        for (auto& v : iA)  v = 0.05f * rk.next();
        for (auto& v : iA2) v = 0.05f * rk.next();
        for (auto& v : iB)  v = 0.05f * rk.next();
        for (auto& v : iC)  v = 0.05f * rk.next();
        iA[0] = 0.70f; iA2[0] = 0.31f; iB[0] = -0.40f;            // distinct heads: says which operator is live
        // [0]/[1] = left/right; D = the engine under test, R = the reference, fed different audio
        std::vector<float> preD[2] { std::vector<float> (coldK + 700),  std::vector<float> (coldK + 700) };
        std::vector<float> preR[2] { std::vector<float> (coldK + 1100), std::vector<float> (coldK + 1100) };
        std::vector<float> prog[2] { std::vector<float> (progN), std::vector<float> (progN) };
        for (auto& ch : preD) for (auto& v : ch) v = 0.5f * rk.next();
        for (auto& ch : preR) for (auto& v : ch) v = 0.4f * rk.next();
        for (auto& ch : prog) for (auto& v : ch) v = 0.3f * rk.next();

        // mode 0 = restart with the publication PENDING · 1 = restart mid-CROSSFADE · 2 = the reference,
        // which lets the publication settle before restarting.
        const auto render = [&] (int nch, int parity, int mode, std::vector<float>& L, std::vector<float>& R,
                                 bool* busyOut, bool* acceptOut)
        {
            convolution::ConvolutionEngine<> e;
            felitronics::test::run (e.prepare (Pk, irMaxK, xfK, nch));
            const std::vector<float>* pre = (mode == 2) ? preR : preD;
            const int blk = 128;
            std::vector<float> jl (pre[0].size() + (std::size_t) blk, 0.0f), jr (jl.size(), 0.0f);
            int fed = 0;
            const auto feed = [&] (int k)                          // rendered into scratch — only `prog` is compared
            {
                for (int done = 0; done < k; )
                {
                    const int m = std::min (blk, k - done);
                    const int at = (fed + done) % (int) pre[0].size();
                    const int mm = std::min (m, (int) pre[0].size() - at);
                    const float* in[2]  { pre[0].data() + at, pre[1].data() + at };
                    float*       out[2] { jl.data(), jr.data() };
                    felitronics::test::run (e.process (in, out, nch, mm));
                    done += mm;
                }
                fed += k;
            };
            if (parity == 1) { test::ok (e.setIr (iA2.data(), irLenK), "parity pre-swap accepted"); feed (settleK); }
            test::ok (e.setIr (iA.data(), irLenK), "operator A accepted");
            feed ((int) pre[0].size());                            // warm the whole FDL and settle on A
            test::ok (e.setIr (iB.data(), irLenK), "operator B published");
            if (mode == 1) feed (xfK / 2);                         // … and picked up: mid-crossfade
            if (mode == 2) feed (settleK);                         // … and settled: the reference state
            e.reset();
            if (busyOut != nullptr) *busyOut = e.isBusy();          // Idle immediately: no swap is left owing
            L.assign (prog[0].size(), 0.0f); R.assign (prog[0].size(), 0.0f);
            for (int done = 0; done < progN; )
            {
                const int m = std::min (blk, progN - done);
                const float* in[2]  { prog[0].data() + done, prog[1].data() + done };
                float*       out[2] { L.data() + done, R.data() + done };
                felitronics::test::run (e.process (in, out, nch, m));
                done += m;
            }
            if (acceptOut != nullptr) *acceptOut = e.setIr (iC.data(), irLenK);   // a NEW publication, at once
        };

        const auto reference = [&] (const std::vector<float>& ir, int ch)
        {
            convolution::PartitionedConvolver<> pc; pc.prepare (Pk, irMaxK); pc.setIr (ir.data(), irLenK);
            std::vector<float> y (prog[ch].size(), 0.0f);
            felitronics::test::run (pc.process (prog[ch].data(), y.data(), progN));
            return y;
        };
        const std::vector<float> yB[2] { reference (iB, 0), reference (iB, 1) };
        const std::vector<float> yA[2] { reference (iA, 0), reference (iA, 1) };

        for (int nch = 1; nch <= 2; ++nch)
            for (int parity = 0; parity <= 1; ++parity)
            {
                std::vector<float> refL, refR;
                render (nch, parity, 2, refL, refR, nullptr, nullptr);
                const std::string cell = "  [" + std::to_string (nch) + " ch, slot parity " + std::to_string (parity) + "]";
                for (int mode = 0; mode <= 1; ++mode)
                {
                    std::vector<float> dut[2]; bool busy = true, accepts = false;
                    render (nch, parity, mode, dut[0], dut[1], &busy, &accepts);
                    const std::vector<float>* ref[2] { &refL, &refR };
                    bool same = true; double worst = 0.0, eB = 0.0, eA = 0.0;
                    for (int c = 0; c < nch; ++c)
                        for (int i = 0; i < progN; ++i)
                        {
                            const float d = dut[c][(std::size_t) i];
                            same  = same && core::sameBits (d, (*ref[c])[(std::size_t) i]);
                            worst = std::max (worst, (double) std::fabs (d - (*ref[c])[(std::size_t) i]));
                            eB    = std::max (eB, (double) std::fabs (d - yB[c][(std::size_t) i]));
                            eA    = std::max (eA, (double) std::fabs (d - yA[c][(std::size_t) i]));
                        }
                    const std::string what = (mode == 0 ? "Pending" : "Crossfading");
                    test::ok (same, "restarted at " + what + ", bit-identical to a settled restart" + cell + " (worst " + sci (worst) + ")");
                    test::ok (eB < 1e-6, "  …and it answers as the PUBLISHED IR" + cell + " (" + sci (eB) + ")");
                    test::ok (eA > 1e-2, "  …and not as the previous one" + cell + " (" + sci (eA) + ")");
                    test::ok (! busy, "  …and Idle right after reset(), so the consumer may publish again" + cell);
                    test::ok (accepts, "  …and a new publication is accepted at once" + cell);
                }
            }
    }

    // --- 🔴 P88: clearAudioState() is the history ALONE — it decides nothing about the operator ---
    // The verb that leaves a fade running (reset() ends it) — running FROM WHERE IT WAS: the fade ends on the
    // sample it was always going to end on (`xfadePos_ = 0` inside this verb would restart it, and a review
    // round showed that surviving the first version of the sibling groups). It must still erase every trace
    // of the old stream, INCLUDING the incoming slot's cached tail, which is blended into every output sample
    // while the fade runs: two engines with different pasts, cleared at the same fade position, answer the
    // same bits.
    // What it does NOT give is independence from WHERE in the fade it was called — that is reset()'s promise,
    // and the last check says so with a number, because a composite choosing between the two verbs needs it.
    test::group ("ConvolutionEngine clearAudioState leaves the swap alone and the history gone");
    {
        const int Pc = 64, irMaxC = 1024, irLenC = 700, xfC = 512;   // a long fade: the restart lands inside it
        const int coldC = ((irMaxC - Pc + Pc - 1) / Pc) * Pc;
        const int progC = xfC + 17 * Pc;
        Lcg rc { 99 };
        std::vector<float> iA (irLenC), iB (irLenC);
        for (auto& v : iA) v = 0.05f * rc.next();
        for (auto& v : iB) v = 0.05f * rc.next();
        iA[0] = 0.70f; iB[0] = -0.40f;
        std::vector<float> preD[2] { std::vector<float> (coldC + 500), std::vector<float> (coldC + 500) };
        std::vector<float> preR[2] { std::vector<float> (coldC + 900), std::vector<float> (coldC + 900) };
        std::vector<float> prog[2] { std::vector<float> (progC), std::vector<float> (progC) };
        for (auto& ch : preD) for (auto& v : ch) v = 0.5f * rc.next();
        for (auto& ch : preR) for (auto& v : ch) v = 0.4f * rc.next();
        for (auto& ch : prog) for (auto& v : ch) v = 0.3f * rc.next();

        const auto render = [&] (int nch, const std::vector<float>* pre, int intoFade, std::vector<float>* out,
                                 bool* busyOut, bool* busyOneShort, bool* busyOnTime)
        {
            convolution::ConvolutionEngine<> e;
            felitronics::test::run (e.prepare (Pc, irMaxC, xfC, nch));
            std::vector<float> jl (pre[0].size(), 0.0f), jr (pre[0].size(), 0.0f);
            const auto feed = [&] (int from, int k)
            {
                const float* in[2]  { pre[0].data() + from, pre[1].data() + from };
                float*       o[2]   { jl.data(), jr.data() };
                felitronics::test::run (e.process (in, o, nch, k));
            };
            test::ok (e.setIr (iA.data(), irLenC), "A accepted");
            feed (0, (int) pre[0].size() - xfC);
            test::ok (e.setIr (iB.data(), irLenC), "B published");
            feed ((int) pre[0].size() - xfC, intoFade);             // into the fade
            e.clearAudioState();
            if (busyOut != nullptr) *busyOut = e.isBusy();           // the fade is NOT cancelled
            out[0].assign ((std::size_t) progC, 0.0f); out[1].assign ((std::size_t) progC, 0.0f);
            const int remaining = xfC - intoFade;                    // the fade ends on this sample, not later
            const auto run = [&] (int from, int k)
            {
                const float* in[2]  { prog[0].data() + from, prog[1].data() + from };
                float*       o[2]   { out[0].data() + from, out[1].data() + from };
                felitronics::test::run (e.process (in, o, nch, k));
            };
            run (0, remaining - 1);   *busyOneShort = e.isBusy();
            run (remaining - 1, 1);   *busyOnTime   = e.isBusy();
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
            int apartCount = 0;
            for (int c = 0; c < nch; ++c)
                for (int i = 0; i < progC; ++i)
                {
                    same   = same && core::sameBits (a[c][(std::size_t) i], b[c][(std::size_t) i]);
                    worst  = std::max (worst, (double) std::fabs (a[c][(std::size_t) i] - b[c][(std::size_t) i]));
                    apart  = std::max (apart, (double) std::fabs (a[c][(std::size_t) i] - late[c][(std::size_t) i]));
                    if (! core::sameBits (a[c][(std::size_t) i], late[c][(std::size_t) i])) ++apartCount;
                }
            std::printf ("      [%d ch] clearAudioState at %d vs %d samples into a %d-sample fade: %d of %d samples differ, worst %.4e\n",
                         nch, xfC / 4, xfC / 2, xfC, apartCount, nch * progC, apart);
            test::ok (same, "two histories, one fade position, identical bits — no trace of either stream" + cell + " (worst " + sci (worst) + ")");
            // …and the fade it left alone still lands on the published operator, on both channels.
            double tailErr = 0.0;
            for (int c = 0; c < nch; ++c)
            {
                convolution::PartitionedConvolver<> refB; refB.prepare (Pc, irMaxC); refB.setIr (iB.data(), irLenC);
                std::vector<float> yB ((std::size_t) progC, 0.0f);
                felitronics::test::run (refB.process (prog[c].data(), yB.data(), progC));
                for (int i = xfC + 2 * Pc; i < progC; ++i) tailErr = std::max (tailErr, (double) std::fabs (a[c][(std::size_t) i] - yB[(std::size_t) i]));
            }
            test::ok (tailErr < 1e-6, "the surviving fade completes onto the published operator" + cell + " (" + sci (tailErr) + ")");
            test::ok (apart > 1e-3, "…but it is NOT independent of where in the fade it was called — reset() is" + cell + " (" + sci (apart) + ")");
        }
    }

    // --- 🔴 P88: the price of adoption is at the SEAM — bounded by the head-tap difference, either way ---
    // The flush itself cuts the stream (the previous tail stops mid-decay), so the first sample after reset()
    // steps whichever slot is live. Adoption replaces the old head tap with the new one at that sample, so
    // on DC it moves the step by AT MOST |Δh0|·|x| — in EITHER direction: a +1 dB move happens to shrink it
    // here and a -1 dB move grows it (a review round caught the text claiming an interactive change never
    // grows it). The first check is the one that tells adoption from the old drop: the first post-reset
    // sample IS the published operator's head.
    test::group ("ConvolutionEngine reset seam: the step moves by at most the head-tap difference");
    {
        const int Ps = 64, irMaxS = 1024, irLenS = 700, xfS = 128;
        const int coldS = ((irMaxS - Ps + Ps - 1) / Ps) * Ps;
        const float dcIn = 0.5f;
        Lcg rs { 5 };
        std::vector<float> iA (irLenS); for (auto& v : iA) v = 0.05f * rs.next();
        iA[0] = 0.70f;
        const auto scaled = [&] (float g) { auto t = iA; for (auto& v : t) v *= g; return t; };
        std::vector<float> sameHead = iA;                         // a shape change that leaves the head tap alone
        for (int i = 1; i < irLenS; ++i) sameHead[(std::size_t) i] += 0.002f * (float) std::sin (0.01 * i);
        std::vector<float> flipped = iA; flipped[0] = -0.40f;     // two ARBITRARY operators
        const auto step = [&] (const std::vector<float>* publish, float* firstOut)
        {
            convolution::ConvolutionEngine<> e;
            felitronics::test::run (e.prepare (Ps, irMaxS, xfS, 1));
            std::vector<float> dc ((std::size_t) coldS + 1600, dcIn), y (dc.size(), 0.0f);
            e.setIr (iA.data(), irLenS);
            felitronics::test::run (e.process (dc.data(), y.data(), (int) dc.size()));
            const double last = y[dc.size() - 1];
            if (publish != nullptr) test::ok (e.setIr (publish->data(), irLenS), "seam publish accepted");
            e.reset();
            float in = dcIn, out = 0.0f;
            felitronics::test::run (e.process (&in, &out, 1));
            *firstOut = out;
            return std::fabs ((double) out - last);
        };
        float f0 = 0.0f;
        const double flushOnly = step (nullptr, &f0);
        std::printf ("      seam step, flush alone: %.4e\n", flushOnly);
        test::ok (flushOnly > 1e-2, "the flush alone already steps at the seam (" + sci (flushOnly) + ")");
        const std::vector<float> up = scaled (1.122f), down = scaled (0.891f);
        const std::pair<const char*, const std::vector<float>*> pairs[] {
            { "same head tap", &sameHead }, { "+1 dB", &up }, { "-1 dB", &down }, { "head tap flipped", &flipped } };
        for (const auto& [label, ir] : pairs)
        {
            float first = 0.0f;
            const double st = step (ir, &first);
            const double bound = std::fabs ((double) (*ir)[0] - (double) iA[0]) * dcIn;
            std::printf ("      seam step, adopting (%s): %.4e\n", label, st);
            test::ok (core::sameBits (first, (*ir)[0] * dcIn), std::string ("first sample after reset() is the PUBLISHED head (") + label + ")");
            test::ok (std::fabs (st - flushOnly) <= bound + 1e-6, std::string ("the step moves by at most |Δh0|·x (") + label + ": "
                                                                    + sci (st) + " vs flush " + sci (flushOnly) + ", bound " + sci (bound) + ")");
        }
    }

    // --- P88: no allocation in reset() / clearAudioState() (both are audio-thread verbs) ---
    test::group ("ConvolutionEngine reset/clearAudioState allocate nothing");
    {
        convolution::ConvolutionEngine<> eng; eng.prepare (P, irMax, xfade, 2);
        std::vector<float> in (512, 0.2f), outL (512, 0.0f), outR (512, 0.0f);
        const float* ins[2]  { in.data(), in.data() };
        float*       outs[2] { outL.data(), outR.data() };
        eng.setIr (irA.data(), irLen);
        felitronics::test::run (eng.process (ins, outs, 2, 512));
        eng.setIr (irB.data(), irLen);                      // a publication in flight, so the adopt path runs
        const long long before = alloc::count.load();
        eng.reset();
        eng.clearAudioState();
        const long long after = alloc::count.load();
        test::okNoAlloc (after == before, "reset() and clearAudioState() performed zero heap allocations");
    }

    // --- unprepared-engine guards (regression: a consumer's ctor pushed an IR before prepare()) ---
    // setIr() before prepare() hit buildIr with P_=0 → `(tailLen + P - 1) / P` divided by zero. ARM64
    // SDIV returns 0 (silent no-op) so macOS was fine; x86-64 IDIV raises #DE → the host process crashed
    // inside the consumer's constructor. setIr()/process() must now REJECT until a successful prepare().
    test::group ("ConvolutionEngine rejects setIr()/process() before prepare()");
    {
        convolution::ConvolutionEngine<> eng;                        // default-constructed, NOT prepared (P_=0)
        std::vector<float> ir (irLen), in (256, 0.3f), out (256, 0.0f);
        for (auto& v : ir) v = 0.1f * r.next();

        test::ok (! eng.isBusy(), "fresh engine is idle");
        test::ok (! eng.setIr (ir.data(), irLen), "setIr() before prepare() returns false (no div-by-0, no crash)");
        test::ok (! eng.isBusy(), "  …and it did NOT arm a swap (state stays idle)");
        test::ok (! eng.process (in.data(), out.data(), 256), "process() before prepare() is REFUSED (law 11)");
        bool untouched = true; for (float v : out) untouched = untouched && (v == 0.0f);
        test::ok (untouched, "process() before prepare() is a no-op (output left untouched)");

        test::ok (eng.prepare (P, irMax, xfade), "prepare() after the rejected loads");
        test::ok (eng.setIr (ir.data(), irLen), "setIr() now accepted once prepared");
        felitronics::test::run (eng.process (in.data(), out.data(), 256));
        bool finite = true; for (float v : out) finite = finite && std::isfinite (v);
        test::ok (finite, "process() produces finite output after prepare");
    }

    // --- a FAILED prepare() must stay unprepared (partial init: P_>0 but FFT setup failed) ---
    // prepare(P=1) passes isPow2(1) and sets P_=1, but the 2-point FFT prepare fails (size < 4) and
    // prepare() returns false. A plain `P_<=0` guard would NOT catch this (P_==1); the prepared_ flag
    // must, so setIr() still rejects rather than writing into the empty (unallocated) channel buffers.
    test::group ("ConvolutionEngine: a failed prepare() stays unprepared");
    {
        convolution::ConvolutionEngine<> eng;
        std::vector<float> ir (irLen, 0.05f);
        test::ok (! eng.prepare (1, irMax, xfade), "prepare(P=1) fails (FFT size 2 < 4)");
        test::ok (! eng.setIr (ir.data(), irLen), "setIr() still rejected after a failed prepare (no OOB into empty h0)");
    }

    return test::report();
}

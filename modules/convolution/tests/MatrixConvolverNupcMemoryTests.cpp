// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// MatrixConvolverNupc's bank spectra follow the LOADED IR: each slot grows, per bank and per stage, to the
// partitions the IR it is handed reaches, and only for the banks its topology reads — while the shared history
// (the FDL) stays sized to the cap. This suite proves what that change could break, and prints what it buys.
//
//   1. AS BEFORE, TO THE BIT. Against a FROZEN copy of the class from before the change
//      (MatrixConvolverNupcFrozen.h), through sequences of operators short, long, short, long — every topology,
//      uneven and fixed block slicing, a swap tried mid-fade, reset() and clearAudioState() mid-fade, a staged
//      operator published late — every output sample of every channel is compared BITWISE, fades included.
//   2. GROWTH IS COMPLETE. A cap-length IR loaded after a short one convolves its whole length: against a
//      direct convolution of the whole input (the history is shared, so the new IR answers input from before
//      the swap), and its LAST tap answers an impulse.
//   3. MEMORY FOLLOWS THE IR. What stageOperator() requests equals an oracle rebuilt here from the schedule's
//      documented rule, per bank and per stage; a shorter or equal IR requests nothing; the topology decides
//      the bank count (mono 1, LRDiag and MSDiag 2, Full 4); prepare() requests exactly the frozen class's
//      prepare() minus the bank spectra.
//   4. process() allocates nothing across any of it.
//   5. A TABLE of bytes per instance, before and after, from FRESH instances (both slots loaded, which is
//      where a live instance settles after its second swap).

#include <felitronics_test.h>
#include <alloc_counter.h>   // installs the allocation counter: EVERY form of `new`, over-aligned included
#include <felitronics/convolution/MatrixConvolverNupc.h>
#include <felitronics/convolution/CabConvolver.h>   // CabConvFft (the backend products ship) and maxIrSamplesFor
#include <felitronics/core/Math.h>                  // core::sameBits — a bitwise compare, not ==
#include "MatrixConvolverNupcFrozen.h"

#include <algorithm>
#include <atomic>
#include <new>
#include <cmath>
#include <cstdio>
#include <string>
#include <type_traits>
#include <vector>

using namespace felitronics;

namespace
{

struct Lcg
{
    unsigned long long s;
    unsigned next() { s = s * 6364136223846793005ULL + 1442695040888963407ULL; return (unsigned) (s >> 33); }
    float bipolar() { return (float) (next() & 0xffffu) / 32768.0f - 1.0f; }
};

// A dense IR of `len` taps: noise under a multiplicative decay (no libm), distinct per seed.
std::vector<float> denseIr (int len, unsigned long long seed)
{
    std::vector<float> h ((std::size_t) std::max (len, 0));
    Lcg r { seed };
    float env = 0.9f;
    const float k = len > 0 ? 1.0f - 4.0f / (float) (len + 64) : 1.0f;
    for (float& v : h) { v = env * r.bipolar(); env *= k; }
    return h;
}

enum class Topo { Mono, LR, MS, Full };

template <class C> typename C::Topology topologyOf (Topo t)
{
    switch (t)
    {
        case Topo::MS:   return C::Topology::MSDiag;
        case Topo::Full: return C::Topology::Full;
        case Topo::Mono:
        case Topo::LR:   return C::Topology::LRDiag;
    }
    return C::Topology::LRDiag;
}

int banksOf (Topo t) { return t == Topo::Mono ? 1 : (t == Topo::Full ? 4 : 2); }

//==============================================================================
// 1 + 4. The lockstep run: the shipping class and the frozen one, fed the same calls.

enum class MidFade { None, TrySwap, Reset, Clear, StageLate, StageResetLate };

// newCap >= 0 re-prepares both instances to that cap (same head) before the swap; newCh > 0 re-prepares them to
// that width as well (the cap stays where it is unless newCap says otherwise).
struct Step { Topo topo; int len; MidFade mid; int newCap = -1; int newCh = 0; };

struct Outcome
{
    long long compared = 0, mismatches = 0, firstMismatch = -1;
    long long processAllocs = 0;
    int callsDiffer = 0, swaps = 0, midActions = 0;
    double peak = 0.0;
    bool finite = true;
};

template <class Fft>
Outcome lockstep (int nch0, int head, int cap, int xfade, const std::vector<Step>& steps, int samplesPerStep,
                  int fixedBlock, unsigned long long seed)
{
    using Now = convolution::MatrixConvolverNupc<Fft>;
    using Was = convolution::frozen::MatrixConvolverNupc<Fft>;
    Outcome r;
    Now now; Was was;
    int nch = nch0, capNow = cap;
    if (! now.prepare (head, cap, xfade, nch) || ! was.prepare (head, cap, xfade, nch)) { r.callsDiffer = -1; return r; }

    Lcg rng { seed };
    constexpr int kMaxChunk = 4096;
    std::vector<float> x[2], yNow[2], yWas[2];
    for (int c = 0; c < 2; ++c) { x[c].assign (kMaxChunk, 0.0f); yNow[c].assign (kMaxChunk, 0.0f); yWas[c].assign (kMaxChunk, 0.0f); }
    long long t = 0;

    for (std::size_t k = 0; k < steps.size(); ++k)
    {
        const Step& st = steps[k];
        std::vector<float> irs[4];
        const float* banks[4] {};
        for (int b = 0; b < 4; ++b) { irs[b] = denseIr (st.len, seed * 977ULL + k * 31ULL + (unsigned long long) b); banks[b] = st.len > 0 ? irs[b].data() : nullptr; }
        if (st.newCap >= 0 || st.newCh > 0)
        {
            if (st.newCap >= 0) capNow = st.newCap;
            if (st.newCh > 0)   nch = st.newCh;
            if (! now.prepare (head, capNow, xfade, nch) || ! was.prepare (head, capNow, xfade, nch)) ++r.callsDiffer;
        }
        const int nb = 4;   // always four bank pointers: a topology that reads fewer must ignore the rest
        const bool late = st.mid == MidFade::StageLate || st.mid == MidFade::StageResetLate;
        const bool a = late ? now.stageOperator (topologyOf<Now> (st.topo), banks, nb, st.len)
                            : now.setOperator   (topologyOf<Now> (st.topo), banks, nb, st.len);
        const bool b = late ? was.stageOperator (topologyOf<Was> (st.topo), banks, nb, st.len)
                            : was.setOperator   (topologyOf<Was> (st.topo), banks, nb, st.len);
        if (a != b) ++r.callsDiffer;
        if (a) ++r.swaps;

        int done = 0;
        bool midDone = st.mid == MidFade::None;
        while (done < samplesPerStep)
        {
            // The first block after a swap is never longer than a third of the fade, so the mid-fade action
            // lands inside the fade at every slicing; every other block is the fixed size or a random one.
            int m = fixedBlock > 0 ? fixedBlock : 1 + (int) (rng.next() % 777u);
            if (done == 0) m = std::min (m, std::max (1, xfade / 3));
            m = std::min ({ m, samplesPerStep - done, kMaxChunk });
            for (int c = 0; c < nch; ++c) for (int i = 0; i < m; ++i) x[c][(std::size_t) i] = 0.5f * rng.bipolar();
            const float* in[2] { x[0].data(), x[1].data() };
            float* on[2] { yNow[0].data(), yNow[1].data() };
            float* ow[2] { yWas[0].data(), yWas[1].data() };

            const long long c0 = alloc::count.load();
            const bool pa = now.process (in, on, nch, m);
            const bool pb = was.process (in, ow, nch, m);
            r.processAllocs += alloc::count.load() - c0;
            if (pa != pb || ! pa) ++r.callsDiffer;

            for (int c = 0; c < nch; ++c)
                for (int i = 0; i < m; ++i)
                {
                    const float u = yNow[c][(std::size_t) i], v = yWas[c][(std::size_t) i];
                    ++r.compared;
                    if (! core::sameBits (u, v)) { if (r.firstMismatch < 0) r.firstMismatch = t + i; ++r.mismatches; }
                    if (! std::isfinite (u)) r.finite = false;
                    r.peak = std::max (r.peak, (double) std::fabs (u));
                }
            done += m; t += m;

            if (! midDone)
            {
                midDone = true;
                switch (st.mid)
                {
                    case MidFade::None: break;
                    case MidFade::TrySwap:   // refused by both while the fade runs, before any work
                    {
                        if (! now.isBusy() || ! was.isBusy()) { ++r.callsDiffer; break; }
                        const long long g0 = alloc::count.load();
                        const bool ra = now.setOperator (topologyOf<Now> (st.topo), banks, nb, st.len);
                        const bool rb = was.setOperator (topologyOf<Was> (st.topo), banks, nb, st.len);
                        if (ra || rb || alloc::count.load() != g0) ++r.callsDiffer;
                        ++r.midActions;
                        break;
                    }
                    case MidFade::Reset:          now.reset();           was.reset();           ++r.midActions; break;
                    case MidFade::Clear:          now.clearAudioState(); was.clearAudioState(); ++r.midActions; break;
                    case MidFade::StageLate:      now.publishStaged();   was.publishStaged();   ++r.midActions; break;
                    case MidFade::StageResetLate: now.reset();           was.reset();
                                                  now.publishStaged();   was.publishStaged();   ++r.midActions; break;
                }
            }
        }
        if (now.isBusy() != was.isBusy()) ++r.callsDiffer;
    }
    return r;
}

// Short, long, short, long with growth in between — then past the cap, cleared, and back.
std::vector<Step> shortLong (Topo t, int cap)
{
    return { { t, 150, MidFade::None },          { t, cap - 1000, MidFade::TrySwap }, { t, 700, MidFade::Reset },
             { t, cap, MidFade::Clear },         { t, 90, MidFade::StageLate },       { t, cap + 9000, MidFade::None },
             { t, 0, MidFade::TrySwap },         { t, 4097, MidFade::StageResetLate },{ t, 2049, MidFade::None },
             { t, cap / 2, MidFade::Reset },     { t, -7, MidFade::None },            { t, 128, MidFade::TrySwap },
             { t, 129, MidFade::None },          { t, 2048, MidFade::Clear },         { t, cap - 2048, MidFade::None } };
}

// Every topology in one instance, so a Full operator lands on banks a LRDiag one never grew, and back.
std::vector<Step> mixed (int cap)
{
    return { { Topo::LR, 150, MidFade::None },        { Topo::LR, cap - 1000, MidFade::TrySwap }, { Topo::Full, 700, MidFade::None },
             { Topo::Full, cap, MidFade::Reset },     { Topo::MS, 300, MidFade::Clear },          { Topo::MS, cap + 9000, MidFade::None },
             { Topo::LR, 0, MidFade::StageLate },     { Topo::Full, 2049, MidFade::None },        { Topo::LR, cap - 1, MidFade::StageResetLate },
             { Topo::Full, 1, MidFade::TrySwap },     { Topo::MS, cap, MidFade::None },           { Topo::Full, cap / 3, MidFade::None } };
}

// Re-prepared in the middle — to half the cap, to twice it, back, and to no tail at all — between swaps.
std::vector<Step> withReprepares (int cap)
{
    return { { Topo::LR, 150, MidFade::None },             { Topo::Full, cap, MidFade::None },
             { Topo::MS, 700, MidFade::None, cap / 2 },    { Topo::Full, cap, MidFade::Reset },
             { Topo::Full, 2 * cap, MidFade::None, 2 * cap }, { Topo::LR, 300, MidFade::Clear, cap },
             { Topo::Full, cap - 1, MidFade::None, 0 },    { Topo::LR, 129, MidFade::TrySwap, cap },
             { Topo::Full, cap, MidFade::None } };
}

// Re-prepared across widths: stereo to mono and back, with and without a cap change — the path where a mono
// re-prepare releases the right channel's history.
std::vector<Step> withWidthChanges (int cap)
{
    return { { Topo::LR, 150, MidFade::None },                 { Topo::Full, cap, MidFade::None },
             { Topo::LR, 700, MidFade::None, -1, 1 },          { Topo::LR, cap, MidFade::Reset },
             { Topo::Full, cap, MidFade::None, -1, 2 },        { Topo::MS, 300, MidFade::Clear },
             { Topo::LR, cap - 1, MidFade::None, cap / 2, 1 }, { Topo::Full, cap, MidFade::TrySwap, cap, 2 } };
}

template <class Fft>
void checkLockstep (const char* backend, const std::string& what, int nch, const std::vector<Step>& steps, int fixedBlock,
                    unsigned long long seed, int head = 128, int cap = 12000)
{
    constexpr int xfade = 400;
    const int perStep = 2 * cap + 3000;   // the fade, then every stage wraps at least once
    const Outcome r = lockstep<Fft> (nch, head, cap, xfade, steps, perStep, fixedBlock, seed);
    long long want = 0;
    for (int width = nch; const Step& s : steps) { if (s.newCh > 0) width = s.newCh; want += (long long) perStep * width; }
    int wantMid = 0; for (const Step& s : steps) if (s.mid != MidFade::None) ++wantMid;
    const std::string tag = std::string (backend) + " " + what + (fixedBlock > 0 ? ", block " + std::to_string (fixedBlock) : ", uneven blocks");
    test::ok (r.callsDiffer == 0, tag + ": every call answered as the frozen class answered (" + std::to_string (r.callsDiffer) + " differ)");
    test::ok (r.compared == want, tag + ": every sample compared (" + std::to_string (r.compared) + " of " + std::to_string (want) + ")");
    test::ok (r.mismatches == 0, tag + ": bit-exact against the frozen class (" + std::to_string (r.mismatches) + " differ, first at "
                                 + std::to_string (r.firstMismatch) + ")");
    test::ok (r.swaps == (int) steps.size() && r.midActions == wantMid,
              tag + ": every swap accepted and every mid-fade action taken (" + std::to_string (r.swaps) + ", " + std::to_string (r.midActions) + ")");
    test::ok (r.finite && r.peak > 0.05, tag + ": the output is finite and not silence (peak " + std::to_string (r.peak) + ")");
    test::okNoAlloc (r.processAllocs == 0, tag + ": process() allocated nothing (" + std::to_string (r.processAllocs) + ")");
}

//==============================================================================
// 2. Growth is complete — against a direct convolution.

struct Tap { int at; float v; };

// A sparse cap-length IR: a tap on each side of every stage boundary, the very last tap, and a scatter.
std::vector<Tap> sparseIr (int cap, unsigned long long seed)
{
    std::vector<Tap> h;
    Lcg r { seed };
    for (int edge = 128; edge < cap; edge *= 2) { h.push_back ({ edge - 1, 0.3f * r.bipolar() }); h.push_back ({ edge, 0.3f * r.bipolar() }); }
    for (int edge = 4096; edge < cap; edge += 2048) { h.push_back ({ edge - 1, 0.2f * r.bipolar() }); h.push_back ({ edge, 0.2f * r.bipolar() }); }
    for (int i = 0; i < 40; ++i) h.push_back ({ (int) (r.next() % (unsigned) cap), 0.1f * r.bipolar() });
    h.push_back ({ 0, 0.25f });
    h.push_back ({ cap - 1, 0.5f });   // the last tap: an impulse must find it
    return h;
}

std::vector<float> dense (const std::vector<Tap>& taps, int len)
{
    std::vector<float> h ((std::size_t) len, 0.0f);
    for (const Tap& t : taps) h[(std::size_t) t.at] += t.v;
    return h;
}

// y[n] = Σ h[k] x[n-k] over the sparse taps, in double.
double directAt (const std::vector<Tap>& taps, const std::vector<float>& x, long long n)
{
    double y = 0.0;
    for (const Tap& t : taps) if (n - t.at >= 0) y += (double) t.v * (double) x[(std::size_t) (n - t.at)];
    return y;
}

template <class Fft>
void checkGrowth (const char* backend, Topo topo)
{
    using C = convolution::MatrixConvolverNupc<Fft>;
    constexpr int cap = 12000, xfade = 64;
    const int nch = topo == Topo::Mono ? 1 : 2;
    const int nb  = banksOf (topo);
    C mc;
    test::run (mc.prepare (128, cap, xfade, nch));

    // The input: noise, silence long enough to empty the history, one impulse on the LEFT input, silence past
    // the last tap. The left impulse reaches the left output through bank 0 in every topology (in MSDiag
    // through M and S together), so the last tap's answer there is 0.5 everywhere.
    const int noiseLen = 4 * cap, impulseAt = noiseLen + cap + 100, total = impulseAt + cap + 500;
    std::vector<float> x[2];
    Lcg rng { 7 };
    for (int c = 0; c < 2; ++c)
    {
        x[c].assign ((std::size_t) total, 0.0f);
        for (int i = 0; i < noiseLen; ++i) x[c][(std::size_t) i] = 0.5f * rng.bipolar();
    }
    x[0][(std::size_t) impulseAt] = 1.0f;
    // MSDiag's banks see the mid and side views of the input, ½(L ± R).
    std::vector<float> mid ((std::size_t) total), side ((std::size_t) total);
    for (int i = 0; i < total; ++i)
    {
        mid[(std::size_t) i]  = 0.5f * (x[0][(std::size_t) i] + x[1][(std::size_t) i]);
        side[(std::size_t) i] = 0.5f * (x[0][(std::size_t) i] - x[1][(std::size_t) i]);
    }

    // Two short operators (300 taps) land in slots 1 and 0; the cap-length one then lands in slot 1 again, so
    // it GROWS a slot that already held a short IR rather than filling an empty one.
    std::vector<float> shortIr[4];
    for (int b = 0; b < 4; ++b) shortIr[b] = denseIr (300, 100ULL + (unsigned long long) b);
    std::vector<Tap> longTaps[4];
    std::vector<float> longIr[4];
    for (int b = 0; b < 4; ++b) { longTaps[b] = sparseIr (cap, 200ULL + (unsigned long long) b); longIr[b] = dense (longTaps[b], cap); }
    const float* sb[4] { shortIr[0].data(), shortIr[1].data(), shortIr[2].data(), shortIr[3].data() };
    const float* lb[4] { longIr[0].data(), longIr[1].data(), longIr[2].data(), longIr[3].data() };

    std::vector<float> y[2];
    for (int c = 0; c < 2; ++c) y[c].assign ((std::size_t) total, 0.0f);
    const int swapShort2 = cap, swapLong = 2 * cap;
    int pos = 0;
    auto runTo = [&] (int end)
    {
        while (pos < end)
        {
            const int m = std::min (333, end - pos);
            const float* in[2] { &x[0][(std::size_t) pos], &x[1][(std::size_t) pos] };
            float* out[2] { &y[0][(std::size_t) pos], &y[1][(std::size_t) pos] };
            test::run (mc.process (in, out, nch, m));
            pos += m;
        }
    };
    const bool s1 = mc.setOperator (topologyOf<C> (topo), sb, nb, 300);   // into slot 1
    runTo (swapShort2);
    const bool s2 = mc.setOperator (topologyOf<C> (topo), sb, nb, 300);   // into slot 0
    runTo (swapLong);
    const bool s3 = mc.setOperator (topologyOf<C> (topo), lb, nb, cap);   // into slot 1, which held 300 taps
    runTo (total);
    test::ok (s1 && s2 && s3, std::string (backend) + " growth: three swaps accepted");

    // The topology's own arithmetic, in double over the sparse taps.
    auto want = [&] (int c, long long n) -> double
    {
        switch (topo)
        {
            case Topo::Full: return c == 0 ? directAt (longTaps[0], x[0], n) + directAt (longTaps[1], x[1], n)
                                           : directAt (longTaps[2], x[0], n) + directAt (longTaps[3], x[1], n);
            case Topo::MS:   return c == 0 ? directAt (longTaps[0], mid, n) + directAt (longTaps[1], side, n)
                                           : directAt (longTaps[0], mid, n) - directAt (longTaps[1], side, n);
            case Topo::Mono:
            case Topo::LR:   return directAt (longTaps[c], x[c], n);
        }
        return 0.0;
    };

    // Settled from the end of the long IR's fade: the exact convolution of the WHOLE input with it — input from
    // before the swap included, because the history is shared.
    double err = 0.0, ref = 0.0;
    for (long long n = swapLong + xfade; n < total; ++n)
        for (int c = 0; c < nch; ++c)
        {
            const double w = want (c, n);
            ref = std::max (ref, std::fabs (w));
            err = std::max (err, std::fabs (w - (double) y[c][(std::size_t) n]));
        }
    const char* name = topo == Topo::Mono ? "mono" : topo == Topo::LR ? "LRDiag" : topo == Topo::MS ? "MSDiag" : "Full";
    std::printf ("      %s %s growth: max error %.3e against a direct convolution (peak %.3e)\n", backend, name, err, ref);
    test::ok (ref > 0.1 && err < 5e-6 * ref, std::string (backend) + " " + name + ": a cap-length IR after a short one is the direct convolution of the whole input");

    // The last tap answers the impulse: a partition missing from the grown bank would leave this at zero.
    const long long at = impulseAt + cap - 1;
    const double lastWant = want (0, at), lastGot = (double) y[0][(std::size_t) at];
    test::ok (std::fabs (lastWant) > 0.3 && std::fabs (lastGot - lastWant) < 5e-6,
              std::string (backend) + " " + name + ": the IR's last tap answers an impulse (" + std::to_string (lastGot) + ", want "
              + std::to_string (lastWant) + ")");
}

//==============================================================================
// 3. Memory follows the IR — an oracle rebuilt from the schedule's documented rule.

struct StageO { int B, C, offset, specF; };

// The rule the class documents: a head of P0, then stages of ONE partition doubling from B = P0, each at
// offset B, up to B_max = max(P0, 2048); the B_max stage takes as many partitions as the rest of the cap needs.
template <class Fft>
std::vector<StageO> scheduleOracle (int P0, int cap)
{
    std::vector<StageO> s;
    long long covered = P0;
    int B = P0;
    const int bMax = std::max (P0, 2048);
    while (covered < (long long) cap)
    {
        const int C = B < bMax ? 1 : (int) (((long long) cap - covered + B - 1) / B);
        s.push_back ({ B, C, (int) covered, Fft::spectrumFloats (2 * B) });
        covered += (long long) C * B;
        if (B >= bMax) break;
        B *= 2;
    }
    return s;
}

int partsO (const StageO& s, int len)
{
    if (len <= s.offset) return 0;
    return (int) std::min<long long> (((long long) len - s.offset + s.B - 1) / s.B, s.C);
}

// Bytes ONE bank requests growing from holding `had` taps to holding `len`: a stage that grows is
// re-reserved WHOLE (exactly the new size); a stage that does not grow requests nothing.
long long growthO (const std::vector<StageO>& s, int had, int len)
{
    long long b = 0;
    for (const StageO& st : s)
        if (partsO (st, len) > partsO (st, had)) b += (long long) partsO (st, len) * st.specF * (long long) sizeof (float);
    return b;
}

long long fullBanksO (const std::vector<StageO>& s)   // one bank at the cap
{
    long long b = 0;
    for (const StageO& st : s) b += (long long) st.C * st.specF * (long long) sizeof (float);
    return b;
}

// Silence through a convolver until its fade is over. The buffers are made once, outside every counting window;
// drain() returns what its process() call allocated, which must be nothing.
struct Drainer
{
    std::vector<float> z, w;
    explicit Drainer (int n) : z ((std::size_t) n, 0.0f), w ((std::size_t) n, 0.0f) {}
    template <class C> long long drain (C& mc, int nch)
    {
        const float* in[2] { z.data(), z.data() };
        float* out[2] { w.data(), w.data() };
        const long long c0 = alloc::count.load();
        test::run (mc.process (in, out, nch, (int) z.size()));
        return alloc::count.load() - c0;
    }
};

template <class Fft>
void checkMemory (const char* backend)
{
    using C = convolution::MatrixConvolverNupc<Fft>;
    using W = convolution::frozen::MatrixConvolverNupc<Fft>;
    const int cap = 240000, xfade = 64, P0 = 128;
    const std::vector<StageO> sch = scheduleOracle<Fft> (P0, cap);
    const std::string be (backend);

    const std::vector<float> ir = denseIr (cap + 5000, 9);
    const float* b4[4] { ir.data(), ir.data(), ir.data(), ir.data() };

    // prepare(): the frozen class's request minus the bank spectra it no longer pre-sizes (2 slots × 4 banks).
    for (int nch = 1; nch <= 2; ++nch)
    {
        const long long w0 = alloc::bytes.load();
        { W was; test::run (was.prepare (P0, cap, xfade, nch)); }
        const long long wasBytes = alloc::bytes.load() - w0;
        const long long n0 = alloc::bytes.load();
        { C now; test::run (now.prepare (P0, cap, xfade, nch)); }
        const long long nowBytes = alloc::bytes.load() - n0;
        const long long spectra = 2LL * (nch == 1 ? 1 : 4) * fullBanksO (sch);
        test::ok (nowBytes == wasBytes - spectra, be + " " + std::to_string (nch) + "ch prepare() requests the frozen class's bytes minus the bank spectra ("
                                                  + std::to_string (nowBytes) + " == " + std::to_string (wasBytes) + " - " + std::to_string (spectra) + ")");
    }

    // The frozen class never allocates at a swap: the whole saving is the shipping class's.
    {
        W was; test::run (was.prepare (P0, cap, xfade, 2));
        const long long a0 = alloc::count.load();
        const bool ok = was.setOperator (W::Topology::Full, b4, 4, cap);
        const long long calls = alloc::count.load() - a0;   // read BEFORE test::ok: its message allocates, and
                                                            // gcc builds that argument first
        test::ok (ok && calls == 0, be + " the frozen class's setOperator() allocates nothing (it pre-sized at the cap)");
    }

    // One stereo instance through a sequence of swaps. Slots alternate: 1, 0, 1, 0, ...
    C mc;
    test::run (mc.prepare (P0, cap, xfade, 2));
    Drainer dr (xfade + 1);
    int held[2][4] {};   // taps each slot's bank was grown for
    long long processAllocs = 0;
    struct Swap { Topo topo; int len; const char* what; };
    const Swap seq[] = {
        { Topo::LR,   61440,      "slot 1, LRDiag 1.28 s, empty: two banks, each the oracle's size" },
        { Topo::LR,   61440,      "slot 0, the same, empty: the same again" },
        { Topo::LR,   10000,      "slot 1, a shorter IR: nothing" },
        { Topo::LR,   61440,      "slot 0, back to its own high-water mark: nothing" },
        // Less than double the stage it grows: a growing resize() alone would reserve twice the old size here.
        { Topo::LR,   100000,     "slot 1, 1.6 times longer: the grown stage re-reserved at exactly its new size" },
        { Topo::LR,   100000,     "slot 0, the same" },
        { Topo::Full, 61440,      "slot 1, Full at the same length: only the two banks LRDiag never read" },
        { Topo::MS,   30000,      "slot 0, MSDiag, shorter, on banks grown for LRDiag: nothing" },
        { Topo::Full, cap,        "slot 1, Full at the cap: the stages that grow, four banks, each re-reserved whole" },
        { Topo::Full, cap,        "slot 0, Full at the cap: its two LRDiag banks grow, its two empty ones fill" },
        { Topo::Full, cap + 5000, "slot 1, past the cap, which it holds: nothing (an IR is cut at the cap)" },
        { Topo::LR,   0,          "slot 0, a zero-length IR: nothing" },
        { Topo::Full, 2049,       "slot 1, a short Full IR after the cap: nothing" },
    };
    int slot = 1;
    for (const Swap& s : seq)
    {
        const int nb = banksOf (s.topo);
        long long want = 0;
        for (int b = 0; b < nb; ++b) { want += growthO (sch, held[slot][b], s.len); held[slot][b] = std::max (held[slot][b], std::min (s.len, cap)); }
        const long long c0 = alloc::count.load(), y0 = alloc::bytes.load();
        const bool ok = mc.setOperator (topologyOf<C> (s.topo), b4, nb, s.len);
        const long long got = alloc::bytes.load() - y0, calls = alloc::count.load() - c0;
        test::ok (ok && got == want && (want != 0 || calls == 0),
                  be + " " + s.what + " (" + std::to_string (got) + " B requested, oracle " + std::to_string (want) + ")");
        processAllocs += dr.drain (mc, 2);
        test::ok (! mc.isBusy(), be + " the fade is over before the next swap");
        slot = 1 - slot;
    }
    test::okNoAlloc (processAllocs == 0, be + " process() allocated nothing across the growth sequence (" + std::to_string (processAllocs) + ")");

    // Mono: one bank.
    {
        C m;
        test::run (m.prepare (P0, cap, xfade, 1));
        const long long y0 = alloc::bytes.load();
        const bool ok = m.setIr (ir.data(), 61440);
        const long long got = alloc::bytes.load() - y0;
        test::ok (ok && got == growthO (sch, 0, 61440), be + " mono: one bank, the oracle's size (" + std::to_string (got) + " B)");
    }

    // prepare() releases what the previous run held. A request counter cannot see a free, but it can see what a
    // release makes necessary: re-prepared at the same geometry, the instance requests its whole history again
    // (assign() alone would reuse the old storage and request nothing), and the next swap requests the oracle's
    // bank bytes again. Re-prepared to a smaller cap, it requests the smaller history exactly.
    {
        auto fdlO = [] (const std::vector<StageO>& s, int nch)
        {
            long long b = 0;
            for (const StageO& st : s) b += (long long) nch * st.C * st.specF * (long long) sizeof (float);
            return b;
        };
        const int small = 48000;
        const std::vector<StageO> schSmall = scheduleOracle<Fft> (P0, small);
        C r;
        test::run (r.prepare (P0, cap, xfade, 2));
        Drainer dr3 (xfade + 1);
        for (int i = 0; i < 2; ++i) { test::run (r.setOperator (C::Topology::Full, b4, 4, cap)); (void) dr3.drain (r, 2); }
        long long y0 = alloc::bytes.load();
        test::run (r.prepare (P0, cap, xfade, 2));
        const long long again = alloc::bytes.load() - y0;
        y0 = alloc::bytes.load();
        const bool ok1 = r.setOperator (C::Topology::Full, b4, 4, cap);
        const long long swapAgain = alloc::bytes.load() - y0;
        (void) dr3.drain (r, 2);
        y0 = alloc::bytes.load();
        test::run (r.prepare (P0, small, xfade, 2));
        const long long down = alloc::bytes.load() - y0;
        y0 = alloc::bytes.load();
        const bool ok2 = r.setOperator (C::Topology::Full, b4, 4, small);
        const long long swapDown = alloc::bytes.load() - y0;
        test::ok (again == fdlO (sch, 2), be + " a re-prepare at the same geometry requests the whole history again (" + std::to_string (again)
                                          + " B, oracle " + std::to_string (fdlO (sch, 2)) + ")");
        test::ok (ok1 && swapAgain == 4 * growthO (sch, 0, cap), be + " ...and the next swap requests the bank bytes again (" + std::to_string (swapAgain) + " B)");
        test::ok (down == fdlO (schSmall, 2), be + " a re-prepare to a smaller cap requests the smaller history exactly (" + std::to_string (down) + " B)");
        test::ok (ok2 && swapDown == 4 * growthO (schSmall, 0, small), be + " ...and its first swap the smaller banks (" + std::to_string (swapDown) + " B)");
    }

    // reserve(): both slots' banks for a topology, now — and then no swap up to that length allocates, whatever
    // its topology. This is what lineareq's EQs call in their prepare().
    {
        C m;
        test::run (m.prepare (P0, cap, xfade, 2));
        long long y0 = alloc::bytes.load();
        const bool ok = m.reserve (C::Topology::Full, cap);
        const long long got = alloc::bytes.load() - y0, want = 2LL * 4 * growthO (sch, 0, cap);
        test::ok (ok && got == want, be + " reserve(Full, cap): both slots' four banks, the oracle's size (" + std::to_string (got) + " B, oracle "
                                     + std::to_string (want) + ")");
        Drainer dr2 (xfade + 1);
        const Topo ts[] { Topo::LR, Topo::Full, Topo::MS, Topo::Full, Topo::LR, Topo::MS };
        const int lens[] { 61440, cap, 1000, cap + 5000, cap, 0 };
        long long later = 0, drains = 0;
        int accepted = 0;
        for (int i = 0; i < 6; ++i)
        {
            const long long c0 = alloc::count.load();
            accepted += m.setOperator (topologyOf<C> (ts[i]), b4, banksOf (ts[i]), lens[i]) ? 1 : 0;
            later += alloc::count.load() - c0;
            drains += dr2.drain (m, 2);
        }
        test::ok (accepted == 6 && later == 0 && drains == 0,
                  be + " after reserve(Full, cap), six swaps of every topology up to the cap allocate nothing (" + std::to_string (later) + ")");

        C mono;
        test::run (mono.prepare (P0, cap, xfade, 1));
        y0 = alloc::bytes.load();
        const bool okMono = mono.reserve (C::Topology::Full, cap);
        const long long gotMono = alloc::bytes.load() - y0, wantMono = 2LL * growthO (sch, 0, cap);
        test::ok (okMono && gotMono == wantMono, be + " reserve() on a mono instance: one bank in each slot (" + std::to_string (gotMono) + " B)");

        C unprepared;
        const long long c0 = alloc::count.load();
        const bool refused = ! unprepared.reserve (C::Topology::Full, cap);
        const long long calls = alloc::count.load() - c0;
        test::ok (refused && calls == 0, be + " reserve() before prepare() is refused and requests nothing");
    }

    // A refused swap touches no storage: a null bank with taps to read is refused before anything grows.
    {
        C m;
        test::run (m.prepare (P0, cap, xfade, 2));
        const float* holes[4] { ir.data(), nullptr, ir.data(), ir.data() };
        const long long c0 = alloc::count.load();
        const bool refused = ! m.setOperator (C::Topology::Full, holes, 4, 61440);
        const long long calls = alloc::count.load() - c0;
        test::ok (refused && calls == 0, be + " a swap refused for a null bank requests nothing (" + std::to_string (calls) + ")");
    }
}

//==============================================================================
// 3b. Exhaustion while growing. Law 11d makes it fatal rather than a refusal and promises nothing after it; this
// pins what the class does give: the throw leaves the call (it is not a terminate at buildBankSlice's noexcept),
// nothing of it reaches the audio, and the next swap builds a correct operator. Every allocation the growth makes
// is failed in turn, then the run is compared to the bit against a frozen instance that never failed.

#if defined(__cpp_exceptions) || defined(_CPPUNWIND)
std::atomic<int> failCountdown { -1 };
void failAfterCountdown() noexcept
{
    if (failCountdown.fetch_sub (1) <= 1) { alloc::failNext.store (true); return; }   // the NEXT request fails
    alloc::onNext.store (&failAfterCountdown);                                       // not yet: re-arm
}

// Arms a failure on the k-th (0-based) allocation from now.
void failAllocation (int k)
{
    if (k == 0) { alloc::failNext.store (true); return; }
    failCountdown.store (k);
    alloc::onNext.store (&failAfterCountdown);
}

template <class Fft>
void checkExhaustion (const char* backend)
{
    using C = convolution::MatrixConvolverNupc<Fft>;
    using W = convolution::frozen::MatrixConvolverNupc<Fft>;
    constexpr int cap = 12000, xfade = 64;
    const std::vector<float> shortIr = denseIr (500, 41), longIr = denseIr (cap, 42);
    const float* sb[4] { shortIr.data(), shortIr.data(), shortIr.data(), shortIr.data() };
    const float* lb[4] { longIr.data(), longIr.data(), longIr.data(), longIr.data() };

    // How many allocations growing from a short LRDiag to a cap-length Full asks for — the sweep's length.
    int requests = 0;
    {
        C mc; test::run (mc.prepare (128, cap, xfade, 2));
        Drainer dr (xfade + 1);
        test::run (mc.setOperator (C::Topology::LRDiag, sb, 2, 500)); (void) dr.drain (mc, 2);
        test::run (mc.setOperator (C::Topology::LRDiag, sb, 2, 500)); (void) dr.drain (mc, 2);   // both slots hold 500
        const long long c0 = alloc::count.load();
        test::run (mc.setOperator (C::Topology::Full, lb, 4, cap));
        requests = (int) (alloc::count.load() - c0);
    }

    // Two ways on from a failed growth, each swept over every request:
    //   RETRY — the instance plays on (a third of the run) before the swap is tried again, then the whole run is
    //           compared with a frozen instance that only ever made the successful swap;
    //   PUBLISH THE STAGED ONE — an operator staged (not published) before the failed swap is published after it:
    //           the failed growth must have left that staging intact — its banks AND its topology label.
    // The staged IR is as long as what the slot held (500 taps), so staging it grows nothing and the failed swap
    // after it makes exactly the requests counted above; its taps differ, so a slot it did not survive in shows.
    const std::vector<float> stagedIr = denseIr (500, 44);
    const float* tb[4] { stagedIr.data(), stagedIr.data(), stagedIr.data(), stagedIr.data() };
    int thrown = 0, recovered = 0;
    long long mismatches = 0, stagedMismatches = 0;
    std::vector<float> x[2], yNow[2], yWas[2];
    constexpr int run = 3 * cap;
    Lcg rng { 43 };
    for (int c = 0; c < 2; ++c) { x[c].resize (run); yNow[c].assign (run, 0.0f); yWas[c].assign (run, 0.0f); for (float& v : x[c]) v = 0.5f * rng.bipolar(); }
    auto play = [&] (C& now, W& was, int from, int to)
    {
        for (int i = from; i < to; i += 512)
        {
            const int m = std::min (512, to - i);
            const float* in[2] { &x[0][(std::size_t) i], &x[1][(std::size_t) i] };
            float* on[2] { &yNow[0][(std::size_t) i], &yNow[1][(std::size_t) i] };
            float* ow[2] { &yWas[0][(std::size_t) i], &yWas[1][(std::size_t) i] };
            test::run (now.process (in, on, 2, m)); test::run (was.process (in, ow, 2, m));
        }
    };
    auto differ = [&]
    {
        long long d = 0;
        for (int c = 0; c < 2; ++c) for (int i = 0; i < run; ++i) if (! core::sameBits (yNow[c][(std::size_t) i], yWas[c][(std::size_t) i])) ++d;
        return d;
    };
    for (int variant = 0; variant < 2; ++variant)
        for (int k = 0; k < requests; ++k)
        {
            C now; W was;
            test::run (now.prepare (128, cap, xfade, 2)); test::run (was.prepare (128, cap, xfade, 2));
            Drainer dr (xfade + 1);
            for (int i = 0; i < 2; ++i)
            {
                test::run (now.setOperator (C::Topology::LRDiag, sb, 2, 500)); (void) dr.drain (now, 2);
                test::run (was.setOperator (W::Topology::LRDiag, sb, 2, 500)); (void) dr.drain (was, 2);
            }
            if (variant == 1) { test::run (now.stageOperator (C::Topology::LRDiag, tb, 2, 500)); test::run (was.stageOperator (W::Topology::LRDiag, tb, 2, 500)); }
            failAllocation (k);
            try { (void) now.setOperator (C::Topology::Full, lb, 4, cap); }
            catch (const std::bad_alloc&) { ++thrown; }
            alloc::failNext.store (false); alloc::onNext.store (nullptr);
            if (variant == 0)
            {
                play (now, was, 0, run / 3);                                     // plays on after the failure
                const bool a = now.setOperator (C::Topology::Full, lb, 4, cap);
                const bool b = was.setOperator (W::Topology::Full, lb, 4, cap);
                if (a && b) ++recovered;
                play (now, was, run / 3, run);
                mismatches += differ();
            }
            else
            {
                now.publishStaged(); was.publishStaged();
                play (now, was, 0, run);
                stagedMismatches += differ();
            }
        }
    const std::string be (backend);
    std::printf ("      %s exhaustion: %d requests in the growth, each failed in turn twice: %d thrown, %d recovered, %lld + %lld samples differ\n",
                 backend, requests, thrown, recovered, mismatches, stagedMismatches);
    test::ok (requests >= 8, be + " exhaustion: the growth makes the requests the sweep walks (" + std::to_string (requests) + ")");
    test::ok (thrown == 2 * requests, be + " exhaustion: every failed request leaves the call as bad_alloc (" + std::to_string (thrown) + " of " + std::to_string (2 * requests) + ")");
    test::ok (recovered == requests && mismatches == 0,
              be + " exhaustion: the instance plays on, the next swap is accepted, and the run is bit-exact against a frozen instance that never failed ("
              + std::to_string (mismatches) + " samples differ)");
    test::ok (stagedMismatches == 0,
              be + " exhaustion: an operator staged before the failed swap publishes intact, bit-exact against the frozen class ("
              + std::to_string (stagedMismatches) + " samples differ)");
}
#endif

//==============================================================================
// 5. The table: bytes per instance, before and after, from FRESH instances. "After" is an instance whose both
// slots have been loaded with the IR (where a live one settles after its second swap); "before" is the frozen
// class, whose bytes are all requested by prepare() — it is checked above that its swaps request none.

template <class Fft>
long long frozenBytes (int cap, int nch)
{
    using W = convolution::frozen::MatrixConvolverNupc<Fft>;
    const long long y0 = alloc::bytes.load();
    W was;
    test::run (was.prepare (128, cap, 64, nch));
    return alloc::bytes.load() - y0;
}

template <class Fft>
long long nowBytes (int cap, int nch, Topo topo, const std::vector<float>& ir, int len)
{
    using C = convolution::MatrixConvolverNupc<Fft>;
    const float* b4[4] { ir.data(), ir.data(), ir.data(), ir.data() };
    Drainer dr (65);
    // Only prepare() and the two swaps are counted: the drains between them allocate nothing, which the
    // sequence above checks, and the Drainer's own buffers are made before the first window.
    long long got = 0;
    C mc;
    long long y0 = alloc::bytes.load();
    test::run (mc.prepare (128, cap, 64, nch));
    got += alloc::bytes.load() - y0;
    bool ok = true;
    for (int k = 0; k < 2; ++k)
    {
        y0 = alloc::bytes.load();
        ok = mc.setOperator (topologyOf<C> (topo), b4, banksOf (topo), len) && ok;
        got += alloc::bytes.load() - y0;
        ok = dr.drain (mc, nch) == 0 && ! mc.isBusy() && ok;
    }
    return ok ? got : -1;
}

template <class Fft>
void memoryTable (const char* backend)
{
    std::printf ("      bytes per instance, MB (1e6 B), fresh instances, both slots loaded — backend %s\n", backend);
    std::printf ("      %-34s %9s %9s %9s %9s\n", "case: before -> after", "44.1 kHz", "48 kHz", "96 kHz", "192 kHz");
    struct Case { const char* name; double capSec; Topo topo; double irSec; };
    const Case cases[] = {
        { "amp, cap 5 s, LRDiag, 1.28 s IR",  5.0, Topo::LR,   1.28 },
        { "amp, cap 5 s, LRDiag, 5 s IR",     5.0, Topo::LR,   5.0  },
        { "cab, cap 4 s, LRDiag, 1.28 s IR",  4.0, Topo::LR,   1.28 },
        { "reverb, cap 20 s, Full, 2 s IR",  20.0, Topo::Full, 2.0  },
        { "reverb, cap 20 s, Full, 4 s IR",  20.0, Topo::Full, 4.0  },
        { "reverb, cap 20 s, Full, 20 s IR", 20.0, Topo::Full, 20.0 },
    };
    const double rates[] = { 44100.0, 48000.0, 96000.0, 192000.0 };
    bool neverMore = true, sameAtCap = true, allMeasured = true;
    for (const Case& c : cases)
    {
        char line[256];
        int at = std::snprintf (line, sizeof line, "      %-34s", c.name);
        for (double sr : rates)
        {
            const int cap = convolution::CabConvolver::maxIrSamplesFor (sr, c.capSec);
            const int len = (int) std::lround (c.irSec * sr);
            const std::vector<float> ir = denseIr (len, 5);
            const long long was = frozenBytes<Fft> (cap, 2);
            const long long now = nowBytes<Fft> (cap, 2, c.topo, ir, len);
            if (now < 0) allMeasured = false;
            if (now > was) neverMore = false;
            if (c.topo == Topo::Full && len >= cap && now != was) sameAtCap = false;
            if (at > 0 && at < (int) sizeof line)
                at += std::snprintf (line + at, sizeof line - (std::size_t) at, " %4.1f>%4.1f", (double) was / 1e6, (double) now / 1e6);
        }
        std::printf ("%s\n", line);
    }
    test::ok (allMeasured, std::string (backend) + " table: every cell's swaps accepted");
    test::ok (neverMore, std::string (backend) + " table: no cell holds more than before");
    test::ok (sameAtCap, std::string (backend) + " table: a Full IR at the cap holds exactly what the frozen class held");
}

template <class Fft>
void suite (const char* backend)
{
    test::group (std::string ("as before, to the bit — against the frozen class (") + backend + ")");
    checkLockstep<Fft> (backend, "mono",   1, shortLong (Topo::Mono, 12000), 0, 11);
    checkLockstep<Fft> (backend, "LRDiag", 2, shortLong (Topo::LR,   12000), 0, 12);
    checkLockstep<Fft> (backend, "MSDiag", 2, shortLong (Topo::MS,   12000), 0, 13);
    checkLockstep<Fft> (backend, "Full",   2, shortLong (Topo::Full, 12000), 0, 14);
    checkLockstep<Fft> (backend, "mixed topologies", 2, mixed (12000), 0,    15);
    checkLockstep<Fft> (backend, "mixed topologies", 2, mixed (12000), 64,   16);
    checkLockstep<Fft> (backend, "mixed topologies", 2, mixed (12000), 4096, 17);
    checkLockstep<Fft> (backend, "mixed topologies", 2, mixed (12000), 1,    18);
    checkLockstep<Fft> (backend, "Full, cap 2049 (the first with a B_max stage)", 2, shortLong (Topo::Full, 2049), 0, 19, 128, 2049);
    checkLockstep<Fft> (backend, "LRDiag, cap 129",     2, shortLong (Topo::LR, 129), 0, 20, 128, 129);
    checkLockstep<Fft> (backend, "mono, cap 0 (a head and no stage)", 1, shortLong (Topo::Mono, 0), 0, 21, 128, 0);
    checkLockstep<Fft> (backend, "Full, head 64",       2, shortLong (Topo::Full, 7000), 0, 22, 64, 7000);
    checkLockstep<Fft> (backend, "mixed, head 512",     2, mixed (30000), 0, 23, 512, 30000);
    checkLockstep<Fft> (backend, "MSDiag, head 16 (pffft's smallest real FFT is 32)", 2, shortLong (Topo::MS, 3000), 0, 24, 16, 3000);
    checkLockstep<Fft> (backend, "through re-prepares", 2, withReprepares (12000), 0, 25);
    checkLockstep<Fft> (backend, "through re-prepares across widths", 2, withWidthChanges (12000), 0, 26);

    test::group (std::string ("growth is complete — against a direct convolution (") + backend + ")");
    checkGrowth<Fft> (backend, Topo::Mono);
    checkGrowth<Fft> (backend, Topo::LR);
    checkGrowth<Fft> (backend, Topo::MS);
    checkGrowth<Fft> (backend, Topo::Full);

    test::group (std::string ("memory follows the IR — against the schedule's rule (") + backend + ")");
    checkMemory<Fft> (backend);

#if defined(__cpp_exceptions) || defined(_CPPUNWIND)
    test::group (std::string ("exhaustion while growing — every request failed in turn (") + backend + ")");
    checkExhaustion<Fft> (backend);
#endif

    test::group (std::string ("memory per instance, before and after (") + backend + ")");
    memoryTable<Fft> (backend);
}

} // namespace

int main()
{
    std::printf ("felitronics::convolution MatrixConvolverNupc — bank spectra follow the loaded IR\n");
    suite<core::fft::DefaultRealFft> ("scalar");
    if constexpr (! std::is_same_v<convolution::CabConvFft, core::fft::DefaultRealFft>)
        suite<convolution::CabConvFft> ("pffft");
    return test::report();
}

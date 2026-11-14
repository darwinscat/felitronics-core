// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// MatrixConvolverNupc SPREADS its last (B_max) stage for an operator long enough to need it: the partitions
// j >= 1 of the next boundary's tail are accumulated through the hop on a schedule of sample positions, and the
// boundary adds j = 0 and inverts. A spread operator picked up over a warm history is prepared for one hop before
// its fade starts. This suite proves what that could break, and prints what it bounds.
//
//   1. A COUNTING BACKEND wraps the shipping FFT and counts forward / inverse / multiply-add calls by transform
//      size. Staging a Full operator transforms exactly the partitions the schedule gives it.
//   2. THE BURST EXISTS (the positive control): the frozen class — the engine before the change — does a whole
//      B_max stage, P partitions of four banks, inside one call of 64 samples, and the counter sees it.
//   3. THE SPREAD BOUND at the shipping threshold, a 20 s Full IR under a 20 s cap at 48 kHz, blocks 1 to 4096:
//      no call does more than its share of the hop's units plus two, idle, nor three times that through a swap;
//      every whole hop does exactly P units — nothing skipped, nothing twice.
//   4. NULLS against an EXTERNAL reference built from frozen instances that hold each operator from the start:
//      the expected output is their smoothstep blend, the fade starting one hop after the pickup (warm and
//      spread) or at the pickup (cold, or unspread). Every topology, pickups at phases 0, 1, B/2 and B-1, fades
//      of 1, B-1, B, B+1 and 3B samples, slot reuse, spread <-> unspread, a topology change. Negative controls:
//      the fade at the pickup (the old start) fails every warm spread case, a shifted start fails every cold one.
//   5. SLICING INDEPENDENCE: one scenario of swaps and restarts, bitwise the same at blocks 1, 64, 333, 4096 and
//      random 1..777, every topology.
//   6. RESTARTS against the same reference: clearAudioState() while priming starts the fade on the next sample;
//      reset() while priming or fading makes the published operator live at once; an operator published and
//      adopted by reset() before any process() plays correctly and spreads (no burst afterwards).
//   7. BELOW THE THRESHOLD the shipping class IS the frozen one: bitwise through swaps, clears, resets and an
//      adoption, the same transforms and multiply-adds in every call.
//   8. maxSwapSamples() bounds every swap above, and a warm spread swap takes exactly crossfade + B.
//   9. process() allocates nothing, in any of it.
//  10. SEEDED RANDOM EVENTS on six geometries (heads 64 to 2048, a B_max stage of 2 to 12 partitions, a partial last
//      one, a single-stage schedule): publications, clears, resets and adoptions at random samples, aimed at the
//      edges — a pickup B-1, B or B+1 samples after the history restarted, a clear on the publish sample, at the
//      fade's start and just before it, a publication on the sample after a fade ends. Each scenario nulls against
//      the frozen references, is busy exactly as documented, and is bitwise the same at random 1..777 and at 4096.

#include <felitronics_test.h>
#include <alloc_counter.h>   // installs the allocation counter: EVERY form of `new`, over-aligned included
#include <felitronics/convolution/MatrixConvolverNupc.h>
#include <felitronics/convolution/CabConvolver.h>   // CabConvFft, the backend products ship
#include <felitronics/core/Math.h>                  // core::sameBits — a bitwise compare, not ==
#include "MatrixConvolverNupcFrozen.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <string>
#include <type_traits>
#include <vector>

using namespace felitronics;

namespace
{

constexpr int kB        = convolution::MatrixConvolverNupc<>::kDefaultMaxBlock;   // the spread stage's hop
constexpr int kLgN      = 12;                                                       // log2 of its transform size
static_assert ((1 << kLgN) == 2 * kB, "the B_max stage transforms 2 * B_max");
constexpr int kHead     = 128;
constexpr int kSmallCap = 6 * kB + kB;      // a B_max stage of six partitions
constexpr int kBigCap   = 960000;           // 20 s at 48 kHz
constexpr double kTol   = 2e-6;             // of the expected output's peak
constexpr double kPeakFloor = 0.5;          // the expected output must be this loud, or a silent engine would pass

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

long long ceilDiv (long long a, long long b) { return (a + b - 1) / b; }

// The schedule's rule at head 128 (the class's comments): stages of 128 .. 1024 taps hold one partition each, so
// the B_max stage starts at tap 2048 and holds ceil((cap - 2048) / 2048) partitions; an IR reaches those whose
// first tap it covers. The spread rule: P >= 2 partitions and (P - 1) * banks >= the threshold.
int stageParts (int cap) { return cap <= kB ? 0 : (int) ceilDiv (cap - kB, kB); }
int partsOf (int len, int cap) { return len <= kB ? 0 : std::min ((int) ceilDiv (len - kB, kB), stageParts (cap)); }
bool spreads (int len, int nb, int cap, int threshold)
{
    const int P = partsOf (len, cap);
    return P >= 2 && (P - 1) * nb >= threshold;
}

//==============================================================================
// 1. The counting backend: the shipping one underneath, every call counted by transform size.

struct FftCounts
{
    long long fwd[32] {}, inv[32] {}, mac[32] {};
    void clear() noexcept { *this = FftCounts {}; }
    bool operator== (const FftCounts&) const = default;
};

int log2Of (int n) noexcept { int k = 0; while (k < 31 && (1 << k) < n) ++k; return k; }

template <class Inner>
struct CountingFft
{
    static inline FftCounts counts {};   // single-threaded tests: plain statics
    static int spectrumFloats (int n) noexcept { return Inner::spectrumFloats (n); }
    bool prepare (int n) noexcept { lg = log2Of (n); return inner.prepare (n); }
    void forward (const float* r, float* w) noexcept { ++counts.fwd[lg]; inner.forward (r, w); }
    void inverse (const float* s, float* w) noexcept { ++counts.inv[lg]; inner.inverse (s, w); }
    void spectralMultiplyAdd (const float* a, const float* b, float* acc) noexcept { ++counts.mac[lg]; inner.spectralMultiplyAdd (a, b, acc); }
    Inner inner {};
    int lg = 0;
};

//==============================================================================
// Operators, signals, scenarios.

enum class Topo { Mono, LR, MS, Full };
int banksOf (Topo t) { return t == Topo::Mono ? 1 : (t == Topo::Full ? 4 : 2); }
int channelsOf (Topo t) { return t == Topo::Mono ? 1 : 2; }
const char* nameOf (Topo t)
{
    switch (t)
    {
        case Topo::Mono: return "mono";
        case Topo::LR:   return "LRDiag";
        case Topo::MS:   return "MSDiag";
        case Topo::Full: return "Full";
    }
    return "?";
}

struct Op { Topo topo; int len; std::array<std::vector<float>, 4> ir; };

Op makeOp (Topo t, int len, unsigned long long seed)
{
    Op o { t, len, {} };
    for (int b = 0; b < banksOf (t); ++b) o.ir[(std::size_t) b] = denseIr (len, seed * 131ULL + (unsigned long long) b);
    return o;
}

template <class C> bool load (C& c, const Op& o)
{
    const float* banks[4] {};
    for (std::size_t b = 0; b < 4; ++b) banks[b] = o.ir[b].empty() ? nullptr : o.ir[b].data();
    typename C::Topology t = C::Topology::LRDiag;
    if (o.topo == Topo::MS)   t = C::Topology::MSDiag;
    if (o.topo == Topo::Full) t = C::Topology::Full;
    return c.setOperator (t, banks, banksOf (o.topo), o.len);
}

struct Signal { std::array<std::vector<float>, 2> x; long long n = 0; };

Signal noise (long long n, unsigned long long seed)
{
    Signal s;
    s.n = n;
    Lcg r { seed };
    for (auto& v : s.x) { v.resize ((std::size_t) n); for (float& f : v) f = 0.5f * r.bipolar(); }
    return s;
}

const float* at (const std::vector<float>& v, long long i) { return v.data() + (std::ptrdiff_t) i; }
float*       at (std::vector<float>& v, long long i)       { return v.data() + (std::ptrdiff_t) i; }

// What happens before the sample at `pos`: a publication (setOperator), a clearAudioState(), a reset(), or a
// publication adopted by a reset() before any process() sees it.
enum class Ev { Publish, Clear, Reset, PublishReset };
struct Event { long long pos; Ev kind; int op; };

struct Blocks
{
    int fixed = 0;               // > 0: every call this long (split at events); 0: random 1..maxRandom
    int maxRandom = 97;
    bool fineWhileBusy = false;  // one sample per call while a swap is in flight: busy time exact to the sample
    unsigned long long seed = 1;
};

struct Swap { long long pickup = 0, busy = -1; };   // busy: from the pickup until isBusy() is false
struct Trace
{
    std::array<std::vector<float>, 2> y;
    std::vector<Swap> swaps;
    long long allocs = 0;        // inside process() only
    bool accepted = true;        // every setOperator() and process() accepted
};

struct NoHook { void operator() (long long, int) const noexcept {} };

// The device under test over a scenario. The pickup of a publication is the first sample of the next call, which
// starts at the event's position: every call is split there.
template <class C, class Hook = NoHook>
Trace drive (C& c, int nch, const std::vector<Op>& ops, const std::vector<Event>& evs, const Signal& s, const Blocks& bl,
             Hook hook = {})
{
    Trace t;
    for (auto& v : t.y) v.assign ((std::size_t) s.n, 0.0f);
    Lcg rng { bl.seed };
    std::size_t e = 0;
    int open = -1;
    for (long long pos = 0; pos < s.n; )
    {
        for (; e < evs.size() && evs[e].pos == pos; ++e)
        {
            const Event& ev = evs[e];
            if (ev.kind == Ev::Clear) { c.clearAudioState(); continue; }
            if (ev.kind == Ev::Reset)
            {
                c.reset();
                if (open >= 0) { t.swaps[(std::size_t) open].busy = pos - t.swaps[(std::size_t) open].pickup; open = -1; }
                continue;
            }
            t.accepted = load (c, ops[(std::size_t) ev.op]) && t.accepted;
            t.swaps.push_back ({ pos, -1 });
            open = (int) t.swaps.size() - 1;
            if (ev.kind == Ev::PublishReset) { c.reset(); t.swaps.back().busy = 0; open = -1; }
        }
        long long stop = s.n;
        if (e < evs.size()) stop = std::min (stop, evs[e].pos);
        const int want = (bl.fineWhileBusy && c.isBusy()) ? 1
                       : (bl.fixed > 0 ? bl.fixed : 1 + (int) (rng.next() % (unsigned) bl.maxRandom));
        const int m = (int) std::min<long long> (want, stop - pos);
        const float* in[2] { at (s.x[0], pos), at (s.x[1], pos) };
        float* out[2] { at (t.y[0], pos), at (t.y[1], pos) };
        const long long a0 = alloc::count.load();
        const bool ok = c.process (in, out, nch, m);
        t.allocs += alloc::count.load() - a0;
        t.accepted = ok && t.accepted;
        hook (pos, m);
        pos += m;
        if (open >= 0 && ! c.isBusy()) { t.swaps[(std::size_t) open].busy = pos - t.swaps[(std::size_t) open].pickup; open = -1; }
    }
    return t;
}

// One frozen instance per operator, holding it from t = 0 over the same input; every clear, reset or adoption in
// the scenario restarts its history at the same sample (they are idle there, so it changes nothing else).
struct Refs { std::vector<std::array<std::vector<float>, 2>> y; bool ok = true; };

template <class Fft>
Refs references (int nch, int cap, int L, const std::vector<Op>& ops, const std::vector<Event>& evs, const Signal& s,
                 int head = kHead)
{
    using F = convolution::frozen::MatrixConvolverNupc<Fft>;
    Refs r;
    r.y.resize (ops.size());
    for (std::size_t k = 0; k < ops.size(); ++k)
    {
        F f;
        r.ok = f.prepare (head, cap, L, nch) && load (f, ops[k]) && r.ok;
        for (auto& v : r.y[k]) v.assign ((std::size_t) s.n, 0.0f);
        std::size_t e = 0;
        for (long long pos = 0; pos < s.n; )
        {
            for (; e < evs.size() && evs[e].pos == pos; ++e)
                if (evs[e].kind != Ev::Publish) f.clearAudioState();
            long long stop = std::min<long long> (s.n, pos + 512);
            if (e < evs.size()) stop = std::min (stop, evs[e].pos);
            const float* in[2] { at (s.x[0], pos), at (s.x[1], pos) };
            float* out[2] { at (r.y[k][0], pos), at (r.y[k][1], pos) };
            r.ok = f.process (in, out, nch, (int) (stop - pos)) && r.ok;
            pos = stop;
        }
    }
    return r;
}

// The engine's smoothstep weight `since` samples into a fade of L (exactly 1 on its last sample).
float smooth (long long since, int L)
{
    const float t = (float) since / (float) (L > 1 ? L - 1 : 1);
    return t * t * (3.0f - 2.0f * t);
}

// Where a fade starts. Documented: one hop after the pickup for a spread operator over a warm history, else at
// the pickup; a clear while priming starts it at once. The others are wrong on purpose (negative controls).
enum class Start { Documented, AtPickup, AlwaysShifted, ClearKeepsStart };

struct Null
{
    double err = 0.0, peak = 0.0;
    bool finite = true, valid = true;
    int warmSpread = 0, coldSpread = 0, unspread = 0, primingCleared = 0, resets = 0;
    std::vector<long long> busy;   // per publication: from the pickup to the fade's end, or the reset that ended it
    double rel() const { return peak > 0.0 ? err / peak : 1e300; }
};

// The expected output, rebuilt from the documented timeline and the references alone — nothing of the engine.
Null null (const Trace& tr, const Refs& refs, const std::vector<Op>& ops, const std::vector<Event>& evs, int nch, int cap,
           int L, int threshold, Start mode)
{
    Null r;
    int live = -1, inc = -1;
    long long ts = 0, pickup = 0, lastClear = 0;
    std::size_t e = 0;
    const long long n = (long long) tr.y[0].size();
    if (n == 0 || refs.y.size() != ops.size()) { r.valid = false; return r; }
    for (long long pos = 0; pos < n; ++pos)
    {
        for (; e < evs.size() && evs[e].pos == pos; ++e)
        {
            const Event& ev = evs[e];
            if (ev.kind == Ev::Clear)
            {
                lastClear = pos;
                if (inc >= 0 && pos < ts && mode != Start::ClearKeepsStart) { ts = pos; ++r.primingCleared; }
                continue;
            }
            if (ev.kind == Ev::Reset)
            {
                lastClear = pos;
                if (inc >= 0) { live = inc; inc = -1; r.busy.back() = pos - pickup; ++r.resets; }
                continue;
            }
            if (inc >= 0) r.valid = false;                 // the engine refuses a publication while busy
            if (live < 0)                                  // the first operator fades in from silence, cold, exactly
            {                                              // as in its own reference
                if (pos != 0 || ev.kind != Ev::Publish) r.valid = false;
                live = ev.op; r.busy.push_back (L);
                continue;
            }
            if (pos < L) r.valid = false;                  // every reference is still fading in from silence there
            if (ev.kind == Ev::PublishReset) { lastClear = pos; live = ev.op; r.busy.push_back (0); ++r.resets; continue; }
            const Op& o = ops[(std::size_t) ev.op];
            const bool sp = spreads (o.len, nch == 1 ? 1 : banksOf (o.topo), cap, threshold);
            const bool warm = pos - lastClear >= kB;       // a B_max boundary has passed since the history restarted
            if (! sp) ++r.unspread; else if (warm) ++r.warmSpread; else ++r.coldSpread;
            const bool shift = mode == Start::AtPickup ? false : (mode == Start::AlwaysShifted ? sp : (sp && warm));
            inc = ev.op; pickup = pos; ts = pos + (shift ? kB : 0); r.busy.push_back (-1);
        }
        const bool fading = inc >= 0 && pos >= ts;
        const float w = fading ? smooth (pos - ts, L) : 0.0f;
        const std::size_t i = (std::size_t) pos;
        for (std::size_t ch = 0; ch < (std::size_t) nch; ++ch)
        {
            double want = refs.y[(std::size_t) live][ch][i];
            if (fading) want = want * (double) (1.0f - w) + (double) refs.y[(std::size_t) inc][ch][i] * (double) w;
            const float got = tr.y[ch][i];
            if (! std::isfinite (got)) r.finite = false;
            r.err  = std::max (r.err, std::fabs ((double) got - want));
            r.peak = std::max (r.peak, std::fabs (want));
        }
        if (fading && pos - ts == L - 1) { live = inc; inc = -1; r.busy.back() = pos + 1 - pickup; }
    }
    return r;
}

std::string num (double v) { char b[32]; std::snprintf (b, sizeof b, "%.3g", v); return b; }

struct Tally
{
    double worst = 0.0, minPeak = 1e300, minNeg = 1e300, minShifted = 1e300;
    int cases = 0, warm = 0, cold = 0, unspread = 0, exactWarm = 0;
    long long allocs = 0;
};

// Groups 4 and 8 for one scenario: the null, its negative controls, and every swap's duration.
template <class Fft>
void nullCase (const std::string& tag, int nch, int L, const std::vector<Op>& ops, const std::vector<Event>& evs,
               long long total, Tally& tl, const Blocks& bl)
{
    using C = convolution::MatrixConvolverNupc<Fft, 1>;
    const Signal s = noise (total, 0x51ULL + (unsigned long long) total);
    C c;
    const bool prepared = c.prepare (kHead, kSmallCap, L, nch);
    const Trace tr = prepared ? drive (c, nch, ops, evs, s, bl) : Trace {};
    const Refs refs = references<Fft> (nch, kSmallCap, L, ops, evs, s);
    const Null n = null (tr, refs, ops, evs, nch, kSmallCap, L, 1, Start::Documented);
    test::ok (prepared && tr.accepted && refs.ok && n.valid && n.finite && n.peak >= kPeakFloor && n.rel() <= kTol,
              tag + ": nulls against the frozen references (error " + num (n.rel()) + " of peak " + num (n.peak) + ")");
    if (n.warmSpread > 0)
    {
        const Null neg = null (tr, refs, ops, evs, nch, kSmallCap, L, 1, Start::AtPickup);
        test::ok (neg.rel() > kTol, tag + ": negative control — the fade at the pickup does NOT null (" + num (neg.rel()) + ")");
        tl.minNeg = std::min (tl.minNeg, neg.rel());
    }
    if (n.coldSpread > 0)
    {
        const Null sh = null (tr, refs, ops, evs, nch, kSmallCap, L, 1, Start::AlwaysShifted);
        test::ok (sh.rel() > kTol, tag + ": negative control — a cold swap's fade one hop late does NOT null (" + num (sh.rel()) + ")");
        tl.minShifted = std::min (tl.minShifted, sh.rel());
    }
    bool busyOk = n.busy.size() == tr.swaps.size();
    for (std::size_t i = 0; busyOk && i < n.busy.size(); ++i)
    {
        busyOk = tr.swaps[i].busy >= 0 && tr.swaps[i].busy == n.busy[i] && tr.swaps[i].busy <= c.maxSwapSamples();
        if (n.busy[i] == L + kB) ++tl.exactWarm;
    }
    test::ok (busyOk && c.maxSwapSamples() == L + kB,
              tag + ": every swap busy exactly as documented, within maxSwapSamples() = " + std::to_string (c.maxSwapSamples()));
    tl.worst = std::max (tl.worst, n.rel()); tl.minPeak = std::min (tl.minPeak, n.peak);
    ++tl.cases; tl.warm += n.warmSpread; tl.cold += n.coldSpread; tl.unspread += n.unspread;
    tl.allocs += tr.allocs;
}

//==============================================================================
// 1. The counter is live and keyed by transform size.

template <class Fft>
void countingIsLive (const char* be)
{
    using CF = CountingFft<Fft>;
    convolution::MatrixConvolverNupc<CF, 1> c;
    test::run (c.prepare (kHead, kSmallCap, 64, 2));
    const Op o = makeOp (Topo::Full, kSmallCap, 3);
    CF::counts.clear();
    test::run (load (c, o));
    const FftCounts k = CF::counts;
    // Head 128: stages of 128, 256, 512, 1024 (transforms 256 .. 2048) hold one partition, the B_max stage six.
    const bool ok = k.fwd[8] == 4 && k.fwd[9] == 4 && k.fwd[10] == 4 && k.fwd[11] == 4 && k.fwd[kLgN] == 4 * partsOf (kSmallCap, kSmallCap)
                 && k.inv[kLgN] == 0 && k.mac[kLgN] == 0;
    test::ok (ok, std::string (be) + " staging a Full operator transforms 4 x the schedule's partitions at each size (B_max: "
                  + std::to_string (k.fwd[kLgN]) + ")");
}

//==============================================================================
// 2. The positive control: the frozen class bursts, and the counter sees it.

template <class Fft>
long long burstExists (const char* be, const Op& a)
{
    using CF = CountingFft<Fft>;
    convolution::frozen::MatrixConvolverNupc<CF> f;
    const int P = partsOf (kBigCap, kBigCap), nb = 4;
    test::run (f.prepare (kHead, kBigCap, 960, 2));
    test::run (load (f, a));
    std::array<std::vector<float>, 2> x, y;
    Lcg r { 9 };
    for (std::size_t ch = 0; ch < 2; ++ch) { x[ch].resize (64); y[ch].assign (64, 0.0f); for (float& v : x[ch]) v = 0.5f * r.bipolar(); }
    const float* in[2] { x[0].data(), x[1].data() };
    float* out[2] { y[0].data(), y[1].data() };
    long long pos = 0, maxCall = 0, hopSum = 0, allocs = 0;
    int hopsWrong = 0, hops = 0;
    for (; pos < 2LL * kB; pos += 64) test::run (f.process (in, out, 2, 64));
    for (; pos < 5LL * kB; )
    {
        CF::counts.clear();
        const long long a0 = alloc::count.load();
        const bool ok = f.process (in, out, 2, 64);
        allocs += alloc::count.load() - a0;
        test::run (ok);
        pos += 64;
        maxCall = std::max (maxCall, CF::counts.mac[kLgN]);
        hopSum += CF::counts.mac[kLgN];
        if (pos % kB == 0) { if (hopSum != (long long) nb * P) ++hopsWrong; ++hops; hopSum = 0; }
    }
    test::ok (maxCall >= (long long) nb * P && hopsWrong == 0 && hops == 3,
              std::string (be) + " frozen, 20 s Full: one call of 64 does " + std::to_string (maxCall) + " B_max MACs (>= P*4 = "
              + std::to_string (nb * P) + "), every hop exactly P*4");
    std::printf ("    %s frozen burst: %lld MACs at N = %d in one call of 64 samples (P = %d)\n", be, maxCall, 2 * kB, P);
    return allocs;
}

//==============================================================================
// 3. The spread bound at the shipping threshold.

long long callBound (int m, int U, int nb) { return (long long) nb * (ceilDiv ((long long) U * m, kB) + ceilDiv (m, kB) + 1); }

template <class Fft>
long long spreadBound (const char* be, const Op& a, const Op& b)
{
    using CF = CountingFft<Fft>;
    using C = convolution::MatrixConvolverNupc<CF>;
    constexpr int L = 960, nb = 4;
    const int P = partsOf (kBigCap, kBigCap), U = P - 1;
    test::ok (spreads (kBigCap, nb, kBigCap, C::kSpreadFromMacs), std::string (be) + " a 20 s Full IR is spread at the shipping threshold");
    std::array<std::vector<float>, 2> x, y;
    Lcg r { 10 };
    for (std::size_t ch = 0; ch < 2; ++ch) { x[ch].resize (4096); y[ch].assign (4096, 0.0f); for (float& v : x[ch]) v = 0.5f * r.bipolar(); }
    long long allocs = 0;
    for (const int n : { 1, 64, 128, 256, 4096 })
    {
        C c;
        bool fine = c.prepare (kHead, kBigCap, L, 2) && load (c, a);
        long long pos = 0, idleMax = 0, swapMax = 0;
        int idleOver = 0, swapOver = 0, hopsWrong = 0, hops = 0, busyWrong = 0, swaps = 0;
        auto call = [&] (int m) -> long long
        {
            const float* in[2] { x[0].data(), x[1].data() };
            float* out[2] { y[0].data(), y[1].data() };
            CF::counts.clear();
            const long long a0 = alloc::count.load();
            const bool ok = c.process (in, out, 2, m);
            allocs += alloc::count.load() - a0;
            fine = ok && fine;
            pos += m;
            return CF::counts.mac[kLgN];
        };
        auto idle = [&] (int m)
        {
            const long long macs = call (m);
            if (m == n) idleMax = std::max (idleMax, macs);   // the printed maximum is of full calls; every call is bounded
            if (macs > callBound (m, U, nb)) ++idleOver;
            return macs;
        };
        while (pos < 2LL * kB) (void) call ((int) std::min<long long> (n, 2LL * kB - pos));   // past the fade-in, warm
        long long hopSum = 0, mark = pos;
        const long long idleEnd = pos + 4LL * kB;   // `pos` moves inside call(): a while, not a for
        while (pos < idleEnd)
        {
            hopSum += idle ((int) std::min<long long> (n, idleEnd - pos));
            if (pos % kB == 0)
            {
                const long long whole = (pos - mark) / kB;
                if (hopSum != whole * nb * P) ++hopsWrong;
                hops += (int) whole; hopSum = 0; mark = pos;
            }
        }
        const Op* next = &b;
        for (const int phase : { 0, 1, 1000, kB - 1 })
        {
            const int d = (int) (((long long) phase - pos % kB + kB) % kB);
            if (d > 0) (void) idle (d);
            fine = ! c.isBusy() && load (c, *next) && fine;
            next = next == &b ? &a : &b;
            const long long pickup = pos;
            long long lastBusy = 0, freeAt = -1;
            while (c.isBusy() && pos - pickup < 4LL * kB)
            {
                const long long macs = call (n);
                swapMax = std::max (swapMax, macs);
                if (macs > 3 * callBound (n, U, nb)) ++swapOver;
                if (c.isBusy()) lastBusy = pos - pickup; else freeAt = pos - pickup;
            }
            // Group 8: the swap ends at pickup + crossfade + B, which the call grid brackets (exactly, at block 1).
            const long long want = L + kB;
            if (! (lastBusy < want && freeAt >= want && (n > 1 || freeAt == want) && want == c.maxSwapSamples())) ++busyWrong;
            ++swaps;
            const long long hopEnd = pos + kB;
            while (pos < hopEnd) (void) idle ((int) std::min<long long> (n, hopEnd - pos));
        }
        const std::string tag = std::string (be) + " block " + std::to_string (n);
        test::ok (fine, tag + ": every call and swap accepted");
        test::ok (idleOver == 0, tag + ": idle, no call over nb*(ceil(U*n/B) + 2) (max " + std::to_string (idleMax) + ", " + std::to_string (idleOver) + " over)");
        test::ok (swapOver == 0, tag + ": through a swap, no call over 3*nb*(ceil(U*n/B) + 2) (max " + std::to_string (swapMax) + ", " + std::to_string (swapOver) + " over)");
        test::ok (hopsWrong == 0 && hops == 4, tag + ": every whole hop exactly nb*P = " + std::to_string (nb * P) + " MACs (" + std::to_string (hops) + " hops)");
        test::ok (busyWrong == 0 && swaps == 4, tag + ": each warm spread swap busy for crossfade + B = maxSwapSamples()");
        std::printf ("    %s block %4d: idle max %5lld MACs/call (bound %5lld), swap max %5lld (bound %5lld)\n", be, n, idleMax,
                     callBound (n, U, nb), swapMax, 3 * callBound (n, U, nb));
    }
    return allocs;
}

//==============================================================================
// 4. Nulls on small geometries, MatrixConvolverNupc<Fft, 1> (spread from P >= 2).

template <class Fft>
long long nulls (const char* be)
{
    Tally tl;
    const std::string b = be;
    const Blocks bl { 0, 97, true, 7 };
    const int kL = 200;
    // Pickups at phases 0, 1, B/2 and B-1 of the B_max stage, every topology.
    for (const Topo tp : { Topo::Mono, Topo::LR, Topo::MS, Topo::Full })
        for (const int ph : { 0, 1, kB / 2, kB - 1 })
        {
            const long long pick = 2LL * kB + ph;
            const std::vector<Op> ops { makeOp (tp, kSmallCap, 11), makeOp (tp, 9000, 12) };
            nullCase<Fft> (b + " " + nameOf (tp) + ", pickup at phase " + std::to_string (ph), channelsOf (tp), kL, ops,
                           { { 0, Ev::Publish, 0 }, { pick, Ev::Publish, 1 } }, pick + 2LL * kB + kL + 300, tl, bl);
        }
    // Fade lengths around the hop.
    for (const Topo tp : { Topo::LR, Topo::Full })
        for (const int L : { 1, kB - 1, kB, kB + 1, 3 * kB })
        {
            const long long pick = (long long) kB * (1 + ceilDiv (L, kB)) + 777;
            const std::vector<Op> ops { makeOp (tp, kSmallCap, 21), makeOp (tp, kSmallCap - 700, 22) };
            nullCase<Fft> (b + " " + nameOf (tp) + ", crossfade " + std::to_string (L), 2, L, ops,
                           { { 0, Ev::Publish, 0 }, { pick, Ev::Publish, 1 } }, pick + kB + L + kB + 100, tl, bl);
        }
    // Slot reuse A -> B -> A -> B, a fade shorter than a hop.
    for (const Topo tp : { Topo::Mono, Topo::Full })
    {
        const int L = 300;
        const std::vector<Op> ops { makeOp (tp, kSmallCap, 31), makeOp (tp, 11000, 32) };
        const long long p1 = 3000, p2 = p1 + kB + L + 111, p3 = p2 + kB + L + 517;
        nullCase<Fft> (b + " " + nameOf (tp) + ", slot reuse A B A B", channelsOf (tp), L, ops,
                       { { 0, Ev::Publish, 0 }, { p1, Ev::Publish, 1 }, { p2, Ev::Publish, 0 }, { p3, Ev::Publish, 1 } },
                       p3 + kB + L + 900, tl, bl);
    }
    // Spread <-> unspread (an unspread operator: one B_max partition, len <= 4096 here).
    {
        const int L = 400;
        const std::vector<Op> ops { makeOp (Topo::Full, 3000, 41), makeOp (Topo::Full, kSmallCap, 42) };
        const long long p1 = 3000, p2 = p1 + kB + L + 200, p3 = p2 + L + 900;
        nullCase<Fft> (b + " Full, unspread -> spread -> unspread -> spread", 2, L, ops,
                       { { 0, Ev::Publish, 0 }, { p1, Ev::Publish, 1 }, { p2, Ev::Publish, 0 }, { p3, Ev::Publish, 1 } },
                       p3 + kB + L + 500, tl, bl);
        const std::vector<Op> mops { makeOp (Topo::Mono, kSmallCap, 43), makeOp (Topo::Mono, 4096, 44) };
        nullCase<Fft> (b + " mono, spread -> unspread -> spread", 1, L, mops,
                       { { 0, Ev::Publish, 0 }, { p1, Ev::Publish, 1 }, { p1 + L + 1000, Ev::Publish, 0 } },
                       p1 + L + 1000 + kB + L + 500, tl, bl);
    }
    // A topology change in the swap: LRDiag -> Full -> MSDiag -> LRDiag (unspread).
    {
        const int L = 250;
        const std::vector<Op> ops { makeOp (Topo::LR, kSmallCap, 51), makeOp (Topo::Full, 10000, 52), makeOp (Topo::MS, kSmallCap, 53),
                                    makeOp (Topo::LR, 2500, 54) };
        const long long p1 = 2500, p2 = p1 + kB + L + 333, p3 = p2 + kB + L + 1;
        nullCase<Fft> (b + " LRDiag -> Full -> MSDiag -> LRDiag", 2, L, ops,
                       { { 0, Ev::Publish, 0 }, { p1, Ev::Publish, 1 }, { p2, Ev::Publish, 2 }, { p3, Ev::Publish, 3 } },
                       p3 + L + 2 * kB, tl, bl);
    }
    // Cold: before the first B_max boundary, and again after a clearAudioState().
    for (const Topo tp : { Topo::Mono, Topo::MS, Topo::Full })
    {
        const int L = 150;
        const std::vector<Op> ops { makeOp (tp, kSmallCap, 61), makeOp (tp, kSmallCap, 62) };
        nullCase<Fft> (b + " " + nameOf (tp) + ", cold swaps", channelsOf (tp), L, ops,
                       { { 0, Ev::Publish, 0 }, { 1000, Ev::Publish, 1 }, { 5000, Ev::Clear, -1 }, { 6000, Ev::Publish, 0 } },
                       6000 + 2 * kB + L, tl, bl);
    }
    test::ok (tl.warm >= 30 && tl.cold >= 6 && tl.unspread >= 3, b + " the nulls cover warm, cold and unspread swaps ("
              + std::to_string (tl.warm) + " / " + std::to_string (tl.cold) + " / " + std::to_string (tl.unspread) + ")");
    std::printf ("    %s nulls: %d cases (%d warm spread, %d cold spread, %d unspread swaps), worst error %.3g of peak (tol %.0e), "
                 "peaks >= %.3g; negative controls >= %.3g (fade at pickup), >= %.3g (cold fade shifted); %d swaps busy exactly crossfade + B\n",
                 be, tl.cases, tl.warm, tl.cold, tl.unspread, tl.worst, kTol, tl.minPeak, tl.minNeg, tl.minShifted, tl.exactWarm);
    return tl.allocs;
}

//==============================================================================
// 5 + 6. One scenario of swaps and restarts: slicing independence, then the restarts against the references.

std::vector<Event> restartScenario()
{
    return {
        { 0,     Ev::Publish, 0 },        // S1, cold
        { 3000,  Ev::Publish, 1 },        // S2, warm spread: fade from 5048
        { 6000,  Ev::Publish, 2 },        // U, unspread: fade at once
        { 7100,  Ev::Publish, 0 },        // S1, warm spread: priming until 9148 ...
        { 7900,  Ev::Clear,  -1 },        // ... cut by a clear: the fade starts at 7900
        { 9000,  Ev::Publish, 1 },        // S2, cold (1100 samples after the clear)
        { 12500, Ev::Publish, 0 },        // S1, warm spread ...
        { 13300, Ev::Reset,  -1 },        // ... reset while priming: S1 live at once
        { 16000, Ev::Publish, 1 },        // S2, warm spread: fade from 18048 ...
        { 18200, Ev::Reset,  -1 },        // ... reset while fading: S2 live at once
        { 21000, Ev::PublishReset, 0 },   // S1 published and adopted before any process()
        { 23500, Ev::Publish, 2 },        // U, warm but unspread
        { 24600, Ev::Publish, 1 },        // S2, warm spread: fade from 26648 ...
        { 26900, Ev::Clear,  -1 },        // ... a clear mid-fade leaves the fade running
    };
}
constexpr long long kRestartTotal = 28500;
constexpr int kRestartL = 500;

template <class Fft>
long long slicingAndRestarts (const char* be)
{
    using C = convolution::MatrixConvolverNupc<Fft, 1>;
    const std::vector<Event> evs = restartScenario();
    long long allocs = 0;
    for (const Topo tp : { Topo::Mono, Topo::LR, Topo::MS, Topo::Full })
    {
        const int nch = channelsOf (tp);
        const std::vector<Op> ops { makeOp (tp, kSmallCap, 71), makeOp (tp, 9000, 72), makeOp (tp, 3000, 73) };
        const Signal s = noise (kRestartTotal, 0x77ULL);
        const std::string tag = std::string (be) + " " + nameOf (tp);
        std::vector<Trace> runs;
        for (const Blocks bl : { Blocks { 1, 1, false, 1 }, Blocks { 64, 1, false, 1 }, Blocks { 333, 1, false, 1 },
                                 Blocks { 4096, 1, false, 1 }, Blocks { 0, 777, false, 5 } })
        {
            C c;
            test::run (c.prepare (kHead, kSmallCap, kRestartL, nch));
            runs.push_back (drive (c, nch, ops, evs, s, bl));
            allocs += runs.back().allocs;
        }
        long long differ = 0;
        bool accepted = true;
        for (const Trace& t : runs)
        {
            accepted = accepted && t.accepted;
            for (std::size_t ch = 0; ch < (std::size_t) nch; ++ch)
                for (std::size_t i = 0; i < (std::size_t) kRestartTotal; ++i)
                    if (! core::sameBits (t.y[ch][i], runs[0].y[ch][i])) ++differ;
        }
        test::ok (accepted && differ == 0, tag + ": bitwise the same at blocks 1, 64, 333, 4096 and random 1..777 ("
                  + std::to_string (differ) + " samples differ)");

        // Group 6 on the block-1 run (busy times exact to the sample).
        const Refs refs = references<Fft> (nch, kSmallCap, kRestartL, ops, evs, s);
        const Null n = null (runs[0], refs, ops, evs, nch, kSmallCap, kRestartL, 1, Start::Documented);
        test::ok (refs.ok && n.valid && n.finite && n.peak >= kPeakFloor && n.rel() <= kTol && n.primingCleared == 1 && n.resets == 3,
                  tag + ": restarts null against the references (error " + num (n.rel()) + " of peak " + num (n.peak) + ")");
        const Null keep = null (runs[0], refs, ops, evs, nch, kSmallCap, kRestartL, 1, Start::ClearKeepsStart);
        test::ok (keep.rel() > kTol, tag + ": negative control — a clear that left the priming running does NOT null (" + num (keep.rel()) + ")");
        bool busyOk = n.busy.size() == runs[0].swaps.size();
        for (std::size_t i = 0; busyOk && i < n.busy.size(); ++i)
            busyOk = runs[0].swaps[i].busy == n.busy[i] && n.busy[i] <= kRestartL + kB;
        test::ok (busyOk, tag + ": every swap of the scenario busy exactly as documented, within maxSwapSamples()");
        if (tp == Topo::Full)
            std::printf ("    %s restarts: error %.3g of peak %.3g; negative control (clear keeps the priming) %.3g\n", be, n.rel(), n.peak, keep.rel());
    }

    // The adopted operator spreads: no burst after reset() adopts it (block 64, Full, P = 6: a whole stage would be 24 MACs).
    {
        using CF = CountingFft<Fft>;
        convolution::MatrixConvolverNupc<CF, 1> c;
        const int P = partsOf (kSmallCap, kSmallCap), U = P - 1, nb = 4;
        const std::vector<Op> ops { makeOp (Topo::Full, kSmallCap, 81), makeOp (Topo::Full, kSmallCap, 82) };
        const long long adopt = 3000, total = adopt + 5LL * kB;
        const Signal s = noise (total, 0x88ULL);
        std::vector<long long> perCall ((std::size_t) total, -1);   // MACs of the call starting at a position
        test::run (c.prepare (kHead, kSmallCap, kRestartL, 2));
        CF::counts.clear();
        const Trace t = drive (c, 2, ops, { { 0, Ev::Publish, 0 }, { adopt, Ev::PublishReset, 1 } }, s, Blocks { 64, 1, false, 1 },
                               [&] (long long pos, int) { perCall[(std::size_t) pos] = CF::counts.mac[kLgN]; CF::counts.clear(); });
        allocs += t.allocs;
        long long maxCall = 0, hopSum = 0;
        int hopsWrong = 0, hops = 0;
        for (long long p = adopt; p < total; ++p)
        {
            const long long m = perCall[(std::size_t) p];
            if (m < 0) continue;
            maxCall = std::max (maxCall, m);
            if (p >= adopt + kB) hopSum += m;
            if (p >= adopt + kB && (p + 64 - adopt) % kB == 0) { if (hopSum != (long long) nb * P) ++hopsWrong; ++hops; hopSum = 0; }
        }
        test::ok (t.accepted && maxCall <= callBound (64, U, nb) && hopsWrong == 0 && hops == 4,
                  std::string (be) + " an operator adopted by reset() spreads: max " + std::to_string (maxCall) + " B_max MACs per call of 64 (bound "
                  + std::to_string (callBound (64, U, nb)) + ", a whole stage " + std::to_string (nb * P) + "), every hop exactly nb*P");
    }
    return allocs;
}

//==============================================================================
// 7. Below the threshold: the frozen class, to the bit and to the call.

template <class Fft>
long long belowThreshold (const char* be)
{
    using CF = CountingFft<Fft>;
    using C = convolution::MatrixConvolverNupc<CF>;
    using F = convolution::frozen::MatrixConvolverNupc<CF>;
    constexpr int L = 960, len = 96000;   // 2 s at 48 kHz
    const int P = partsOf (len, kBigCap);
    test::ok (! spreads (len, 4, kBigCap, C::kSpreadFromMacs) && ! spreads (len, 2, kBigCap, C::kSpreadFromMacs)
              && spreads (kBigCap, 4, kBigCap, C::kSpreadFromMacs),
              std::string (be) + " 2 s operators are under the threshold (P = " + std::to_string (P) + ": " + std::to_string ((P - 1) * 4)
              + " < " + std::to_string (C::kSpreadFromMacs) + ") at a cap that can spread");
    const std::vector<Op> ops { makeOp (Topo::Full, len, 91), makeOp (Topo::Full, len, 92), makeOp (Topo::MS, len, 93) };
    // A clear mid-fade, a reset mid-fade, an adoption by reset() before any process(), a clear on a publish sample.
    const std::vector<Event> evs { { 0, Ev::Publish, 0 },    { 3000, Ev::Publish, 1 }, { 3500, Ev::Clear, -1 },
                                   { 5737, Ev::Publish, 2 }, { 6200, Ev::Reset, -1 },  { 8000, Ev::PublishReset, 1 },
                                   { 9056, Ev::Publish, 0 }, { 12000, Ev::Publish, 2 }, { 12000, Ev::Clear, -1 } };
    constexpr long long total = 16000;
    const Signal s = noise (total, 0x99ULL);
    C c; F f;
    bool fine = c.prepare (kHead, kBigCap, L, 2) && f.prepare (kHead, kBigCap, L, 2);
    std::array<std::vector<float>, 2> yc, yf;
    for (std::size_t ch = 0; ch < 2; ++ch) { yc[ch].assign ((std::size_t) total, 0.0f); yf[ch].assign ((std::size_t) total, 0.0f); }
    Lcg rng { 3 };
    std::size_t e = 0;
    long long allocs = 0, pickup = -1, busyMax = 0;
    int callsDiffer = 0, busyDiffer = 0, calls = 0;
    for (long long pos = 0; pos < total; )
    {
        for (; e < evs.size() && evs[e].pos == pos; ++e)
        {
            const Event& ev = evs[e];
            if (ev.kind == Ev::Clear) { c.clearAudioState(); f.clearAudioState(); continue; }
            if (ev.kind == Ev::Reset) { c.reset(); f.reset(); pickup = -1; continue; }
            fine = load (c, ops[(std::size_t) ev.op]) && load (f, ops[(std::size_t) ev.op]) && fine;
            pickup = pos;
            if (ev.kind == Ev::PublishReset) { c.reset(); f.reset(); pickup = -1; }
        }
        long long stop = total;
        if (e < evs.size()) stop = std::min (stop, evs[e].pos);
        const int m = (int) std::min<long long> (1 + (int) (rng.next() % 500u), stop - pos);
        const float* in[2] { at (s.x[0], pos), at (s.x[1], pos) };
        float* oc[2] { at (yc[0], pos), at (yc[1], pos) };
        float* of[2] { at (yf[0], pos), at (yf[1], pos) };
        CF::counts.clear();
        fine = f.process (in, of, 2, m) && fine;
        const FftCounts was = CF::counts;
        CF::counts.clear();
        const long long a0 = alloc::count.load();
        const bool ok = c.process (in, oc, 2, m);
        allocs += alloc::count.load() - a0;
        fine = ok && fine;
        if (! (CF::counts == was)) ++callsDiffer;
        pos += m; ++calls;
        if (c.isBusy() != f.isBusy()) ++busyDiffer;
        if (pickup >= 0 && ! c.isBusy()) { busyMax = std::max (busyMax, pos - pickup); pickup = -1; }
    }
    long long differ = 0;
    for (std::size_t ch = 0; ch < 2; ++ch)
        for (std::size_t i = 0; i < (std::size_t) total; ++i) if (! core::sameBits (yc[ch][i], yf[ch][i])) ++differ;
    test::ok (fine && differ == 0, std::string (be) + " under the threshold: bitwise the frozen class through six publications, two clears, a reset and an adoption (" + std::to_string (differ) + " samples differ)");
    test::ok (callsDiffer == 0, std::string (be) + " ...the same forward / inverse / multiply-add counts at every size in every call ("
              + std::to_string (callsDiffer) + " of " + std::to_string (calls) + " differ)");
    test::ok (busyDiffer == 0 && busyMax <= c.maxSwapSamples() && c.maxSwapSamples() == L + kB,
              std::string (be) + " ...busy in lockstep with it, every swap within maxSwapSamples() (" + std::to_string (busyMax) + " <= "
              + std::to_string (c.maxSwapSamples()) + ")");
    return allocs;
}

//==============================================================================

// THE THRESHOLD'S EDGE, at the shipping value: one partition below it runs the old path (the burst at the boundary,
// an immediate fade), one partition more spreads (no burst, the fade a hop later). Full: (P-1)*4 = 252 against 256;
// LRDiag: (P-1)*2 = 254 against 256 — the 264192 taps CabConvolver's comment states. A new operator published
// while the spread one is being prepared is refused, and the convolver says it is busy.
template <class Fft>
long long thresholdEdge (const char* be)
{
    using C = convolution::MatrixConvolverNupc<CountingFft<Fft>>;
    constexpr int xfade = 64;
    long long allocs = 0;
    struct Edge { Topo topo; int cap, below, at; };
    const Edge edges[] { { Topo::Full, kB + 65 * kB, kB + 64 * kB, kB + 65 * kB },
                         { Topo::LR,   kB + 129 * kB, kB + 128 * kB, kB + 129 * kB } };
    for (const Edge& e : edges)
        for (int len : { e.below, e.at })
        {
            const bool spread = len == e.at;
            const int nb = banksOf (e.topo), P = partsOf (len, e.cap);
            const std::string tag = std::string (be) + " " + nameOf (e.topo) + ", P " + std::to_string (P) + " ((P-1)*" + std::to_string (nb)
                                    + " = " + std::to_string ((P - 1) * nb) + ")";
            test::ok (spreads (len, nb, e.cap, C::kSpreadFromMacs) == spread, tag + ": the oracle puts it on the side the class must");
            C mc;
            test::run (mc.prepare (kHead, e.cap, xfade, 2));
            const Op a = makeOp (e.topo, len, 77), b = makeOp (e.topo, len, 78);
            test::run (load (mc, a));
            const Signal x = noise (9LL * kB, 79);   // warm 3 hops, steady 2, then room for a swap of a fade and a hop
            std::vector<float> yl ((std::size_t) kB), yr ((std::size_t) kB);
            long long pos = 0, maxMacs = 0;
            auto call = [&] (int n)
            {
                const float* in[2] { &x.x[0][(std::size_t) pos], &x.x[1][(std::size_t) pos] };
                float* out[2] { yl.data(), yr.data() };
                const long long m0 = CountingFft<Fft>::counts.mac[kLgN], c0 = alloc::count.load();
                test::run (mc.process (in, out, 2, n));
                allocs += alloc::count.load() - c0;
                maxMacs = std::max (maxMacs, CountingFft<Fft>::counts.mac[kLgN] - m0);
                pos += n;
            };
            while (pos < 3LL * kB) call (64);                      // warm: the first fade over, the history full
            maxMacs = 0;
            while (pos < 5LL * kB) call (64);                      // steady: two whole hops
            const long long bound = callBound (64, P - 1, nb);
            if (spread) test::ok (maxMacs <= bound, tag + ": spread, no call above the bound (" + std::to_string (maxMacs) + " <= " + std::to_string (bound) + ")");
            else        test::ok (maxMacs >= (long long) P * nb, tag + ": not spread, the boundary does the whole stage (" + std::to_string (maxMacs) + ")");
            test::run (load (mc, b));
            long long busy = 0;
            bool refused = true;
            do
            {
                call (1); ++busy;
                if (busy == 1 || busy == kB / 2)                   // just after the pickup, and mid-preparation
                {
                    const bool busyNow = mc.isBusy();
                    refused = refused && busyNow && ! load (mc, a);
                }
            } while (mc.isBusy() && pos < 9LL * kB);
            test::ok (! mc.isBusy(), tag + ": the swap ended inside the signal");
            const long long want = spread ? xfade + kB : xfade;
            test::ok (busy == want, tag + ": the swap is busy " + std::to_string (busy) + " samples (want " + std::to_string (want) + ")");
            test::ok (refused, tag + ": a publish while the swap is in flight is refused, and isBusy() says so");
        }
    return allocs;
}

//==============================================================================
// 10. Seeded random events, MatrixConvolverNupc<Fft, 1>, against the same references and the same timeline.

struct Gen
{
    Lcg r;
    unsigned u (unsigned m) { return r.next() % m; }
};

// How often the scenarios hit the edges they are aimed at (gated below, so a generator that drifted cannot pass).
struct Cover
{
    int coldEdge[3] {};   // a spread operator picked up B-1, B, B+1 samples after the history restarted
    int clearOnPublish = 0, clearAtStart = 0, clearBeforeStart = 0, publishAfterFade = 0;
    int resetPriming = 0, resetFading = 0, adoptions = 0, scenarios = 0;
};

long long randomGap (Gen& g)
{
    switch (g.u (8))
    {
        case 0:  return 0;
        case 1:  return 1;
        case 2:  return kB - 1;
        case 3:  return kB;
        case 4:  return kB + 1;
        case 5:  return (long long) g.u (64);
        default: return (long long) g.u (2 * kB);
    }
}

// Events placed by the documented timeline (the rule null() checks): a publication only once the last swap has
// ended — its fade starting one hop after the pickup for a spread operator over a warm history, else at once, and
// at once from a clear made while it was being prepared.
std::vector<Event> randomScenario (Gen& g, const std::vector<Op>& ops, int nch, int cap, int L, long long& total, Cover& cv)
{
    std::vector<Event> evs { { 0, Ev::Publish, 0 } };
    long long pos = L, lastClear = 0, ts = 0, end = 0, pickup = 0;
    bool busy = false;
    int live = 0;
    auto publish = [&] (long long p)
    {
        int op = (int) g.u ((unsigned) ops.size());
        if (op == live) op = (op + 1) % (int) ops.size();
        const Op& o = ops[(std::size_t) op];
        const bool sp = spreads (o.len, nch == 1 ? 1 : banksOf (o.topo), cap, 1);
        const long long since = p - lastClear;
        if (sp && since >= kB - 1 && since <= kB + 1) ++cv.coldEdge[since - (kB - 1)];
        ts = p + ((sp && since >= kB) ? kB : 0); end = ts + L; pickup = p; busy = true; live = op;
        evs.push_back ({ p, Ev::Publish, op });
    };
    const int nEv = 3 + (int) g.u (6);
    for (int i = 0; i < nEv; ++i)
    {
        if (busy)
        {
            const unsigned k = g.u (10);
            if (k < 5)   // a clear (k < 3) or a reset inside the swap
            {
                long long p = pos + 1 + (long long) g.u ((unsigned) std::max<long long> (end - pos - 1, 1));
                const unsigned w = g.u (4);
                if (k < 3 && w == 0 && pos == pickup)     { p = pickup; ++cv.clearOnPublish; }   // state 1: not picked up yet
                else if (w == 1 && ts > pos)              { p = ts; if (k < 3) ++cv.clearAtStart; }
                else if (k < 3 && w == 2 && ts - 1 > pos) { p = ts - 1; ++cv.clearBeforeStart; }
                evs.push_back ({ p, k < 3 ? Ev::Clear : Ev::Reset, -1 });
                if (k >= 3 && p < end) { if (p < ts) ++cv.resetPriming; else ++cv.resetFading; }
                if (k >= 3 || p >= end) busy = false;
                else if (p < ts) { ts = p; end = p + L; }
                lastClear = p; pos = p;
                continue;
            }
            pos = end; busy = false;                                   // wait it out
            if (g.u (3) == 0) { publish (pos); ++cv.publishAfterFade; continue; }   // on the sample after the fade
        }
        const unsigned k = g.u (12);
        const long long p = pos + randomGap (g);
        if (k < 2) { evs.push_back ({ p, k == 0 ? Ev::Clear : Ev::Reset, -1 }); lastClear = p; pos = p; }
        else if (k == 2)
        {
            const int op = (int) g.u ((unsigned) ops.size());
            evs.push_back ({ p, Ev::PublishReset, op }); lastClear = p; live = op; pos = p; ++cv.adoptions;
        }
        else if (k == 3)   // the cold edge: a clear, then a publication B-1, B or B+1 samples on
        {
            evs.push_back ({ p, Ev::Clear, -1 }); lastClear = p;
            pos = p + kB - 1 + (long long) g.u (3);
            publish (pos);
        }
        else { publish (p); pos = p; }
    }
    total = (busy ? std::max (pos, end) : pos) + kB + 300 + (long long) g.u (700);   // past the last swap's first boundary
    ++cv.scenarios;
    return evs;
}

template <class Fft>
long long randomEvents (const char* be, Cover& cv)
{
    using C = convolution::MatrixConvolverNupc<Fft, 1>;
    struct Geo { int head, cap, cases; const char* what; };
    // The time-domain head dominates the cost at heads 1024 and 2048 (a tap per sample per bank): fewer cases there.
    // The events, the timeline and every decision the engine takes are the backend's business only in arithmetic, and
    // the scalar FFT is several times slower than pffft: when it is not the backend products ship, it runs the first
    // third of the scenarios.
    const Geo geos[] { { 128, 7 * kB, 40, "head 128, cap 7B (C = 6)" },         { 128, 3 * kB, 40, "head 128, cap 3B (C = 2)" },
                       { 128, 4 * kB + 5, 40, "head 128, cap 4B+5 (C = 4, a partial last partition)" },
                       { 64, 13 * kB, 30, "head 64, cap 13B (C = 12)" },        { 1024, 3 * kB - 7, 5, "head 1024, cap 3B-7 (C = 2)" },
                       { 2048, 5 * kB, 3, "head 2048, cap 5B (a single stage, C = 4)" } };
    const bool shipping = std::is_same_v<Fft, convolution::CabConvFft>;
    static constexpr int Ls[] { 1, 2, 37, 300, kB - 1, kB, kB + 1, 2500 };
    long long allocs = 0;
    int gi = 0;
    for (const Geo& geo : geos)
    {
        ++gi;
        int bad = 0, printed = 0;
        double worst = 0.0, minPeak = 1e300;
        long long samples = 0, swaps = 0;
        const int cases = shipping ? geo.cases : (geo.cases + 2) / 3;
        for (int ci = 0; ci < cases; ++ci)
        {
            Gen g { { 0xF00DULL * (unsigned long long) gi + (unsigned long long) ci * 7919ULL + 17ULL } };
            const int nch = g.u (4) == 0 ? 1 : 2;
            const int L = Ls[g.u (8)];
            const int nOps = 2 + (int) g.u (3);
            std::vector<Op> ops;
            for (int k = 0; k < nOps; ++k)
            {
                Topo tp = Topo::Mono;
                if (nch == 2) { const unsigned t = g.u (3); tp = t == 0 ? Topo::LR : (t == 1 ? Topo::MS : Topo::Full); }
                static constexpr int fr[] { 100, 3000, 4096, 4097, 6145, 0, 0, 0 };
                int len = fr[g.u (8)];
                if (len == 0) len = kB + 1 + (int) g.u ((unsigned) (geo.cap + 600 - kB - 1));
                if (k == 0) len = geo.cap;   // the first operator is the cap's length: every partition, loud
                ops.push_back (makeOp (tp, len, 1000ULL * (unsigned long long) gi + 10ULL * (unsigned long long) ci + (unsigned long long) k));
            }
            long long total = 0;
            const std::vector<Event> evs = randomScenario (g, ops, nch, geo.cap, L, total, cv);
            const Signal s = noise (total, 0x51ULL + (unsigned long long) total + (unsigned long long) ci);
            C c1, c2, c3;
            const bool prepared = c1.prepare (geo.head, geo.cap, L, nch) && c2.prepare (geo.head, geo.cap, L, nch)
                               && c3.prepare (geo.head, geo.cap, L, nch);
            if (! prepared) { ++bad; continue; }
            const Trace t1 = drive (c1, nch, ops, evs, s, Blocks { 1, 1, false, 1 });   // busy exact to the sample
            const Trace t2 = drive (c2, nch, ops, evs, s, Blocks { 0, 777, false, 99 + (unsigned long long) ci });
            const Trace t3 = drive (c3, nch, ops, evs, s, Blocks { 4096, 1, false, 1 });
            allocs += t1.allocs + t2.allocs + t3.allocs;
            const Refs refs = references<Fft> (nch, geo.cap, L, ops, evs, s, geo.head);
            const Null n = null (t1, refs, ops, evs, nch, geo.cap, L, 1, Start::Documented);
            bool busyOk = n.busy.size() == t1.swaps.size();
            for (std::size_t i = 0; busyOk && i < n.busy.size(); ++i)
                busyOk = t1.swaps[i].busy >= 0 && t1.swaps[i].busy == n.busy[i] && t1.swaps[i].busy <= c1.maxSwapSamples();
            long long differ = 0;
            for (std::size_t ch = 0; ch < (std::size_t) nch; ++ch)
                for (std::size_t i = 0; i < (std::size_t) total; ++i)
                    differ += (core::sameBits (t2.y[ch][i], t1.y[ch][i]) ? 0 : 1) + (core::sameBits (t3.y[ch][i], t1.y[ch][i]) ? 0 : 1);
            const bool good = t1.accepted && t2.accepted && t3.accepted && refs.ok && n.valid && n.finite && n.peak >= kPeakFloor
                           && n.rel() <= kTol && busyOk && differ == 0;
            if (! good)
            {
                ++bad;
                if (printed++ < 3)
                    std::printf ("    BAD %s %s case %d: nch %d L %d, accepted %d/%d/%d, refs %d, valid %d, finite %d, peak %.3g, error %.3g, busy %d, %lld differ\n",
                                 be, geo.what, ci, nch, L, t1.accepted, t2.accepted, t3.accepted, refs.ok, n.valid, n.finite, n.peak, n.rel(),
                                 busyOk, differ);
            }
            worst = std::max (worst, n.rel()); minPeak = std::min (minPeak, n.peak);
            samples += total; swaps += (long long) t1.swaps.size();
        }
        test::ok (bad == 0, std::string (be) + " " + geo.what + ": " + std::to_string (cases) + " random scenarios null (worst "
                  + num (worst) + " of peak), busy exact, bitwise at 777 and 4096 (" + std::to_string (bad) + " bad)");
        std::printf ("    %s %-52s %3d scenarios, %4lld swaps, %8lld samples: worst error %.3g of peak (tol %.0e), peaks >= %.3g\n",
                     be, geo.what, cases, swaps, samples, worst, kTol, minPeak);
    }
    return allocs;
}

template <class Fft>
void suite (const char* backend)
{
    const Op a = makeOp (Topo::Full, kBigCap, 1001), b = makeOp (Topo::Full, kBigCap, 1002);
    const std::string be = backend;

    test::group ("the counting backend sees every transform, by size (" + be + ")");
    countingIsLive<Fft> (backend);

    test::group ("the burst exists — the frozen class, 20 s Full under a 20 s cap (" + be + ")");
    const long long a2 = burstExists<Fft> (backend, a);

    test::group ("the spread bound — shipping threshold, 20 s Full, blocks 1 to 4096 (" + be + ")");
    const long long a3 = spreadBound<Fft> (backend, a, b);

    test::group ("nulls against the frozen references — every topology, phase, fade length, slot reuse (" + be + ")");
    const long long a4 = nulls<Fft> (backend);

    test::group ("slicing independence and restarts (" + be + ")");
    const long long a5 = slicingAndRestarts<Fft> (backend);

    test::group ("below the threshold — the frozen class, to the bit and to the call (" + be + ")");
    const long long a7 = belowThreshold<Fft> (backend);

    test::group ("the threshold's edge — one partition either side, Full and LRDiag (" + be + ")");
    const long long a8 = thresholdEdge<Fft> (backend);

    test::group ("seeded random events — six geometries, aimed at the edges (" + be + ")");
    Cover cv;
    const long long a10 = randomEvents<Fft> (backend, cv);
    test::ok (cv.coldEdge[0] > 0 && cv.coldEdge[1] > 0 && cv.coldEdge[2] > 0 && cv.clearOnPublish > 0 && cv.clearAtStart > 0
              && cv.clearBeforeStart > 0 && cv.publishAfterFade > 0 && cv.resetPriming > 0 && cv.resetFading > 0 && cv.adoptions > 0,
              be + " the random scenarios hit every edge they aim at");
    std::printf ("    %s %d scenarios: pickups at B-1 / B / B+1 after a restart %d / %d / %d; clears on the publish sample %d, at the fade's "
                 "start %d, just before %d; publications right after a fade %d; resets priming %d, fading %d; adoptions %d\n",
                 backend, cv.scenarios, cv.coldEdge[0], cv.coldEdge[1], cv.coldEdge[2], cv.clearOnPublish, cv.clearAtStart, cv.clearBeforeStart,
                 cv.publishAfterFade, cv.resetPriming, cv.resetFading, cv.adoptions);

    test::group ("process() allocates nothing (" + be + ")");
    test::okNoAlloc (a8 == 0, be + " threshold-edge runs (" + std::to_string (a8) + ")");
    test::okNoAlloc (a10 == 0, be + " random-event runs (" + std::to_string (a10) + ")");
    test::okNoAlloc (a2 == 0, be + " frozen burst run (" + std::to_string (a2) + ")");
    test::okNoAlloc (a3 == 0, be + " spread bound runs (" + std::to_string (a3) + ")");
    test::okNoAlloc (a4 == 0, be + " null runs (" + std::to_string (a4) + ")");
    test::okNoAlloc (a5 == 0, be + " slicing and restart runs (" + std::to_string (a5) + ")");
    test::okNoAlloc (a7 == 0, be + " below-threshold runs (" + std::to_string (a7) + ")");
}

} // namespace

int main()
{
    std::printf ("felitronics::convolution MatrixConvolverNupc — the B_max stage spread across the hop\n");
    suite<core::fft::DefaultRealFft> ("scalar");
    if constexpr (! std::is_same_v<convolution::CabConvFft, core::fft::DefaultRealFft>)
        suite<convolution::CabConvFft> ("pffft");
    return test::report();
}

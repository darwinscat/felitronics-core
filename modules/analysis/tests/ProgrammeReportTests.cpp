// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// ProgrammeReport self-tests. Three kinds of evidence, deliberately not one:
//
//   1. CLOSED-FORM ORACLES, computed outside this object and outside this repository's code. A pure sine's
//      crest factor is 20·log10(√2) = 3.0103 dB; an LR4 at 30 Hz passes |1/(1+(f/30)⁴)²| of a tone's power,
//      so 40 Hz reads 5.78 % and 60 Hz 0.346 %; a 1 kHz sine at −23 dBFS reads −23.0 LUFS (EBU Tech 3341
//      case 1); a programme built so the last 100 ms carries exactly four times the mean square reads a
//      tail ratio of 4. These are the checks that survive a schedule which is deterministically WRONG.
//   2. AN INDEPENDENT WHOLE-FILE REFERENCE — the same quantities computed by a second program: one flat
//      pass over the planes, no streaming, no scratch, no ring, no absent-channel path, and percentiles
//      taken by SORTING rather than from a histogram. It nulls the streaming machinery, not itself.
//   3. RE-SPLIT INVARIANCE (law 8a) over the report AND over the whole frame trace, bit for bit.
//      Invariance alone certifies self-consistency and nothing else — a deterministically wrong schedule is
//      equally wrong at every slicing — which is why (1) and (2) sit beside it and not behind it.
//
// THE MUTATION STAND — 23 mutants, run on this suite (2026-09-15, Release, AppleClang 21, 48 kHz
// fixtures). **22 DIE. The one survivor is equivalent by construction** and is named below rather than
// left looking like a hole. Eleven of the twenty-two were SURVIVORS on an earlier run of this stand and
// are dead only because two crew rounds went looking; each one's gate carries a comment at the group that
// kills it, saying what it survived.
//
//   THE SCHEDULE                                                                         checks failed
//     M1  the sub-hop closes at the end of process(), not on its audio-time boundary            50
//     M2  Crossover2::flushDenormals() moves off core::StateGrid to the end of process()        14
//     M4  the observation stride becomes 101 sub-hops instead of 100                            15
//     M14 the 10 ms sub-hop is floor(0.01·fs) instead of lround                                  4
//   THE SPAN
//     M16 the span's end excludes the last active frame                                          6
//     M17 the LEADING quiet run is committed into the span                                      12
//     M18 the span goes back to the exact-zero predicate                                        17
//     M22 DC and Σx² are summed over every sample instead of over the span                      12
//   THE TAIL
//     M3  the tail window ends at T instead of at the last active frame                          8
//     M7  the deferred run's skip (`tailWritten_ += pending - keep`) is dropped                  2
//     M9  the deferred frames are written back with a sample count of 0                          7
//     M9b …or with an energy of 0 (a sub-threshold frame is quiet, not silent)                   1
//     M10 the minimum tail length tests the FILE position, not the programme span                1
//     M13 the low-band pending run is never committed on an active frame                         5
//     M15 the tail ring is read one entry short                                                 12
//   THE STATISTICS AND THE VERDICTS
//     M5a binOf's TOP clamp is off by one                                                        2
//     M5b binOf's BOTTOM clamp is off by one                                                     6
//     M6  the two silence predicates' trace ids are swapped                                      2
//     M8  the damage-reason ladder in buildDynamics() is reordered                                3
//     M19 samplePeak stays valid with nothing finite to take a maximum over                      1
//     M20 the true peak is gated at 0 instead of at gainToDb's 1e-12 floor                       1
//     M21 "no gating block existed" is reported as SilentProgramme                               1
//
//   SURVIVED, AND WHY THAT IS NOT A HOLE:
//     M12 the absolute gate's `>=` becoming `>`. The two differ only on an observation whose mean square
//         is EXACTLY 1.1724653045822981e-07; no float programme through two recursive filters produces
//         that value, and constructing it would need a seam inside the gate. `>=` is libebur128's
//         spelling and is kept for that reason, not because a test forced it.
//     (An earlier stand also listed `percentile`'s own top clamp. It is unreachable — the histogram holds
//      `total` entries and `rank <= total - 1`, so the running sum passes `rank` at a populated bin — so
//      it is not counted as a mutant here; the LIVE clamp is binOf's, which M5a and M5b cover.)
//
//   **WHERE THEY DIED IS THE POINT: not one of the 22 produced a single "report bit-identical" failure.**
//   Under M1 and M2 all 19 slicings and all 8 maxBlock values still read a BIT-IDENTICAL finished report.
//   Comparing the finished report is not a gate for law 8a — docs/LAW8-KWEIGHTING.md:78 said so from a
//   measurement, and this is that measurement again. What killed them was the trace, the closed forms and
//   the whole-file reference; M17 and M22, both real defects, were caught by the reference alone.
//
//   Each mutant was built with its object file DELETED first and the build output checked for a
// `Building CXX` line, then the pristine header was restored and re-run to confirm it is green: make's
// mtime granularity is one second, so a header edited inside the same second as the previous build does
// not rebuild, and the stand would then report the OLD binary's result as the mutant's.

#include <felitronics_test.h>
#include <alloc_counter.h>   // installs the allocation counter: EVERY form of `new`, over-aligned included
#include <felitronics/analysis/BandCrest.h>
#include <felitronics/analysis/ProgrammeReport.h>

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <random>
#include <string>
#include <vector>

using namespace felitronics;
using PR     = analysis::ProgrammeReport;
using Params = analysis::ProgrammeReportParams;
using Reason = analysis::ProgrammeReason;
using Event  = analysis::ProgrammeTraceEvent;
using Planes = std::vector<std::vector<float>>;

static constexpr double kTwoPi = 6.283185307179586;
static constexpr double kFs    = 48000.0;

// ============================================================================================== helpers
static std::uint64_t b64 (double d) noexcept { return std::bit_cast<std::uint64_t> (d); }

static std::vector<const float*> ptrs (const Planes& p, std::size_t at)
{
    std::vector<const float*> v (p.size());
    for (std::size_t c = 0; c < p.size(); ++c) v[c] = p[c].data() + at;
    return v;
}

struct Out
{
    bool               prepared = false;
    bool               accepted = true;
    PR::Report         R {};
    std::vector<Event> trace;
    std::int64_t       traceOverflow = 0;   // an accessor, not a report field: it is the BUFFER's property
    long               allocsInProcess = 0;
};

// One streamed run. `split` is cycled for the call lengths; `prepCh` defaults to the plane count; `feedCh`
// is the width handed to process() (default = prepCh), so a narrower feed exercises the absent-channel path.
static Out run (const Planes& p, double fs, int maxBlock, const std::vector<int>& split,
                Params pm = {}, int prepCh = -1, int feedCh = -1, std::size_t traceCap = 1u << 16)
{
    Out o;
    const int nplanes = (int) p.size();
    const int pc = prepCh < 0 ? nplanes : prepCh;
    const int fc = feedCh < 0 ? std::min (pc, nplanes) : feedCh;

    std::vector<Event> buf (traceCap);                    // allocated BEFORE process(): never inside it
    PR a;
    a.setParams (pm);
    if (! a.prepare (fs, maxBlock, pc)) return o;
    o.prepared = true;
    a.setTraceBuffer (buf.data(), buf.size());

    const std::size_t frames = p.empty() ? 0u : p[0].size();
    const long long before = alloc::count.load();
    std::size_t at = 0, k = 0;
    while (at < frames)
    {
        const int n = (int) std::min<std::size_t> ((std::size_t) split[k++ % split.size()], frames - at);
        const auto v = ptrs (p, at);
        if (! a.process (v.data(), fc, n)) { o.accepted = false; break; }
        at += (std::size_t) n;
    }
    a.finish();
    o.allocsInProcess = alloc::count.load() - before;
    o.R = a.report();
    o.traceOverflow = a.traceOverflow();
    o.trace.assign (buf.begin(), buf.begin() + (std::size_t) a.traceCount());
    return o;
}

// --- the two field enumerations, materialised so two reports can be walked in lockstep ---
struct FieldV { std::string name; int ch; analysis::ProgrammeValue v; };
struct FieldC { std::string name; int ch; std::int64_t v; };

static std::vector<FieldV> valuesOf (const PR::Report& R)
{
    std::vector<FieldV> out;
    R.visitValues ([&] (const char* n, int c, const analysis::ProgrammeValue& v)
                   { out.push_back (FieldV { n, c, v }); });
    return out;
}
static std::vector<FieldC> countsOf (const PR::Report& R)
{
    std::vector<FieldC> out;
    R.visitCounts ([&] (const char* n, int c, std::int64_t v) { out.push_back (FieldC { n, c, v }); });
    return out;
}

// RELATIVE agreement with the whole-file reference, and why it is not bit-equality. The two programs sum
// in DIFFERENT ORDERS on purpose: the subject commits each quiet run into a compensated accumulator when
// the next active frame arrives, the reference walks the span straight through. Demanding the same bits
// would force the oracle to copy the subject's schedule, and an oracle that copies the schedule cannot
// test the schedule — which is the whole reason it exists. So: agreement to a relative tolerance, and
// EXACT only where exactness is genuinely claimed (a maximum, and the correlation, which both compute with
// the same shipped `StereoSums`).
static void relClose (double got, double want, double rel, const std::string& msg)
{
    const double tol = rel * std::max (std::fabs (want), std::fabs (got));
    test::ok (std::fabs (got - want) <= tol,
              msg + " (got " + std::to_string (got) + ", want " + std::to_string (want) + ")");
}

// BIT comparison, field by field — not memcmp (padding) and not `!=` (which cannot tell −0.0 from +0.0).
static bool sameReport (const PR::Report& a, const PR::Report& b, std::string& why)
{
    const auto va = valuesOf (a), vb = valuesOf (b);
    if (va.size() != vb.size()) { why = "value field count"; return false; }
    for (std::size_t i = 0; i < va.size(); ++i)
    {
        if (va[i].name != vb[i].name || va[i].ch != vb[i].ch) { why = "field order at " + va[i].name; return false; }
        if (b64 (va[i].v.value) != b64 (vb[i].v.value) || va[i].v.valid != vb[i].v.valid
            || va[i].v.reason != vb[i].v.reason)
        {
            why = va[i].name + "[" + std::to_string (va[i].ch) + "]";
            return false;
        }
    }
    const auto ca = countsOf (a), cb = countsOf (b);
    if (ca.size() != cb.size()) { why = "count field count"; return false; }
    for (std::size_t i = 0; i < ca.size(); ++i)
        if (ca[i].name != cb[i].name || ca[i].ch != cb[i].ch || ca[i].v != cb[i].v)
        {
            why = ca[i].name + "[" + std::to_string (ca[i].ch) + "] "
                + std::to_string (ca[i].v) + " vs " + std::to_string (cb[i].v);
            return false;
        }
    if (b64 (a.sampleRate) != b64 (b.sampleRate)) { why = "sampleRate"; return false; }
    return true;
}

static bool sameTrace (const std::vector<Event>& a, const std::vector<Event>& b, std::string& why)
{
    if (a.size() != b.size())
    {
        why = "trace length " + std::to_string (a.size()) + " vs " + std::to_string (b.size());
        return false;
    }
    for (std::size_t i = 0; i < a.size(); ++i)
    {
        const Event& x = a[i]; const Event& y = b[i];
        if (x.kind != y.kind || x.channel != y.channel || x.index != y.index || x.begin != y.begin
            || x.end != y.end || b64 (x.a) != b64 (y.a) || b64 (x.b) != b64 (y.b) || b64 (x.c) != b64 (y.c))
        {
            why = "trace[" + std::to_string (i) + "] kind " + std::to_string ((int) x.kind)
                + " ch " + std::to_string (x.channel) + " idx " + std::to_string (x.index);
            return false;
        }
    }
    return true;
}

// ============================================================================================== signals
static Planes silence (int nch, std::size_t n) { return Planes ((std::size_t) nch, std::vector<float> (n, 0.0f)); }

static Planes tone (int nch, std::size_t n, double hz, double amp, double fs = kFs, double phase = 0.0)
{
    Planes p ((std::size_t) nch, std::vector<float> (n));
    for (std::size_t i = 0; i < n; ++i)
    {
        const float v = (float) (amp * std::sin (kTwoPi * hz * (double) i / fs + phase));
        for (auto& ch : p) ch[i] = v;
    }
    return p;
}

// Non-stationary material with a fixed seed: bursts, gaps, transients and a slow drift. A steady tone is a
// weak witness of invariance — it looks the same through every window.
static Planes programme (int nch, std::size_t n, unsigned seed)
{
    Planes p ((std::size_t) nch, std::vector<float> (n, 0.0f));
    std::mt19937 rng (seed);
    std::uniform_real_distribution<double> u (-1.0, 1.0);
    std::size_t i = 0;
    while (i < n)
    {
        const std::size_t len = 200u + (std::size_t) (rng() % 9000u);
        const bool quiet = (rng() % 5u) == 0u;
        const double amp = quiet ? 1.0e-6 : 0.02 + 0.6 * ((double) (rng() % 1000u) / 1000.0);
        const double hz  = 30.0 + (double) (rng() % 4000u);
        for (std::size_t j = 0; j < len && i < n; ++j, ++i)
            for (std::size_t c = 0; c < p.size(); ++c)
            {
                const double env = std::exp (-4.0 * (double) j / (double) len);
                const double s = amp * env * std::sin (kTwoPi * hz * (double) i / kFs + 0.7 * (double) c)
                               + 0.004 * u (rng);
                p[c][i] = (float) s;
            }
        if ((rng() % 3u) == 0u)                          // a hard transient right on a boundary
            for (std::size_t c = 0; c < p.size() && i > 0; ++c) p[c][i - 1] = (c == 0 ? 0.9f : -0.85f);
    }
    return p;
}

static Planes pad (const Planes& p, std::size_t lead, std::size_t trail)
{
    Planes q (p.size());
    for (std::size_t c = 0; c < p.size(); ++c)
    {
        q[c].assign (lead, 0.0f);
        q[c].insert (q[c].end(), p[c].begin(), p[c].end());
        q[c].insert (q[c].end(), trail, 0.0f);
    }
    return q;
}

// ==================================================================== the independent whole-file reference
// A SECOND PROGRAM for the same quantities: one flat pass, no streaming state, no scratch, no ring, and
// percentiles taken by SORTING. It shares the two filter primitives — those ARE the definition — but not
// the schedule, the span bookkeeping, the tail placement or the statistics.
struct Ref
{
    std::vector<double> dc, rms, peak, crest, infra;
    double   progMeanSquare = 0.0, tailRatio = 0.0, crestAll = 0.0, rmsAll = 0.0, peakAll = 0.0;
    double   corr = 0.0, balanceDb = 0.0, sideToMid = 0.0;
    double   p10 = 0.0, p50 = 0.0, p95 = 0.0, lra = 0.0;
    std::int64_t span = 0, leadDigital = 0, trailDigital = 0, leadThr = 0, trailThr = 0;
    std::int64_t observations = 0, gated = 0;
    std::vector<double> obsSeries;                  // the raw short-term mean squares, pre-log
    bool     tailValid = false, statsValid = false;
};

static Ref reference (const Planes& p, double fs, const Params& pm, int channels)
{
    Ref r;
    const std::size_t T = p.empty() ? 0u : p[0].size();
    const std::size_t nch = (std::size_t) channels;
    const double thr = std::pow (10.0, pm.silenceThresholdDb / 20.0);
    const auto present = [&] (std::size_t c) { return c < p.size(); };

    // --- the exact-zero bounds: threshold-free evidence, and they drive ONLY the digital counts ---
    std::int64_t firstZ = -1, lastZ = -1;
    for (std::size_t i = 0; i < T; ++i)
    {
        bool nonZero = false;
        for (std::size_t c = 0; c < nch; ++c)
            if (present (c) && ! (p[c][i] == 0.0f)) nonZero = true;      // a NaN is not a zero either
        if (nonZero) { if (firstZ < 0) firstZ = (std::int64_t) i; lastZ = (std::int64_t) i; }
    }
    r.leadDigital  = firstZ < 0 ? (std::int64_t) T : firstZ;
    r.trailDigital = lastZ  < 0 ? (std::int64_t) T : (std::int64_t) T - 1 - lastZ;

    // --- THE SPAN is the thresholded region: the first to the last frame whose loudest finite present
    //     channel exceeds the threshold. Same definition as the subject's, derived here from the parameter
    //     rather than read from its state.
    std::int64_t first = -1, last = -1;
    for (std::size_t i = 0; i < T; ++i)
    {
        double mx = 0.0;
        for (std::size_t c = 0; c < nch; ++c)
            if (present (c) && std::isfinite (p[c][i])) mx = std::max (mx, std::fabs ((double) p[c][i]));
        if (mx > thr) { if (first < 0) first = (std::int64_t) i; last = (std::int64_t) i; }
    }
    r.span     = last < 0 ? 0 : last + 1 - first;
    r.leadThr  = first < 0 ? (std::int64_t) T : first;
    r.trailThr = last  < 0 ? (std::int64_t) T : (std::int64_t) T - 1 - last;

    // --- per channel over the span, plus the infra-low split driven flat with its own grid ---
    r.dc.assign (nch, 0.0); r.rms.assign (nch, 0.0); r.peak.assign (nch, 0.0);
    r.crest.assign (nch, 0.0); r.infra.assign (nch, 0.0);
    double totSq = 0.0, totLow = 0.0;
    std::int64_t totCount = 0;

    eq::Crossover2 xo;
    xo.prepare (fs, channels);
    xo.setFrequency ((float) pm.infraLowHz);
    core::StateGrid grid;
    std::vector<double> chSq (nch, 0.0), chLow (nch, 0.0), chDc (nch, 0.0);
    std::vector<std::int64_t> chCount (nch, 0);
    for (std::size_t i = 0; i < T; ++i)
    {
        const bool inSpan = first >= 0 && (std::int64_t) i >= first && (std::int64_t) i <= last;
        for (std::size_t c = 0; c < nch; ++c)
        {
            const float x = present (c) && std::isfinite (p[c][i]) ? p[c][i] : 0.0f;
            float lo = 0.0f, hi = 0.0f;
            xo.processSample ((int) c, x, lo, hi);
            if (! present (c) || ! std::isfinite (p[c][i])) continue;
            r.peak[c] = std::max (r.peak[c], std::fabs ((double) x));
            if (! inSpan) continue;
            chDc[c] += (double) x;
            chSq[c] += (double) x * (double) x;
            chLow[c] += (double) lo * (double) lo;
            ++chCount[c];
        }
        if (grid.advance (1)) xo.flushDenormals();
    }
    for (std::size_t c = 0; c < nch; ++c)
    {
        if (chCount[c] <= 0) continue;
        r.dc[c]  = chDc[c] / (double) chCount[c];
        r.rms[c] = std::sqrt (chSq[c] / (double) chCount[c]);
        r.crest[c] = r.rms[c] > 0.0 ? 20.0 * std::log10 (r.peak[c] / r.rms[c]) : 0.0;
        r.infra[c] = chSq[c] > 0.0 ? chLow[c] / chSq[c] : 0.0;
        totSq += chSq[c]; totLow += chLow[c]; totCount += chCount[c];
        r.peakAll = std::max (r.peakAll, r.peak[c]);
    }
    if (totCount > 0)
    {
        r.progMeanSquare = totSq / (double) totCount;
        r.rmsAll = std::sqrt (r.progMeanSquare);
        r.crestAll = r.rmsAll > 0.0 ? 20.0 * std::log10 (r.peakAll / r.rmsAll) : 0.0;
    }

    // --- the tail window, placed from the last signal frame ---
    const std::int64_t N = (std::int64_t) std::lround (fs * pm.tailWindowMs / 1000.0);
    // The window must fit inside the PROGRAMME, not merely inside the file: `last + 1 >= N` would admit a
    // window made almost entirely of the LEADING silence. (The first version of this reference had exactly
    // the implementation's own bug here, which is how a shared mistake hides — the review round found it in
    // both at once.)
    if (last >= 0 && (last + 1 - first) >= N)
    {
        double s = 0.0; std::int64_t cnt = 0;
        for (std::int64_t i = last + 1 - N; i <= last; ++i)
            for (std::size_t c = 0; c < nch; ++c)
                if (present (c) && std::isfinite (p[c][(std::size_t) i]))
                { const double x = (double) p[c][(std::size_t) i]; s += x * x; ++cnt; }
        if (cnt > 0 && r.progMeanSquare > 0.0)
        {
            r.tailRatio = (s / (double) cnt) / r.progMeanSquare;
            r.tailValid = true;
        }
    }

    // --- stereo, straight from the two planes ---
    if (nch >= 2 && p.size() >= 2)
    {
        analysis::StereoSums s;
        for (std::size_t i = 0; i < T; ++i)
            if (std::isfinite (p[0][i]) && std::isfinite (p[1][i])) s.add (p[0][i], p[1][i]);
        r.corr = s.correlation();
        if (s.ll > 0.0 && s.rr > 0.0) r.balanceDb = 10.0 * std::log10 (s.ll / s.rr);
        if (s.mid > 0.0) r.sideToMid = s.side / s.mid;
    }

    // --- the short-term series, its own K-weighting, and percentiles BY SORTING ---
    // THE CADENCE IS WRITTEN OUT HERE ON PURPOSE. The first version of this reference read
    // `PR::kShortTermSubHops` and `PR::kObservationHops` from the header, and a crew testing round showed
    // what that costs: changing the observation stride left all 350 checks GREEN, because the "independent"
    // reference imported the mutated constant and moved with it instead of opposing it. An oracle that reads
    // the subject's own parameters is the subject. BOTH numbers are written out here: 300 sub-hops is the 3 s
    // window and 10 is the 100 ms step, and EBU Tech 3342 §3.1 is the source of both — "a sliding
    // analysis-window of length 3 seconds" and "a minimum block overlap of 2.9 s … i.e. >=10 Hz sampling".
    // The stride was 100 until v0.43, one second, which came from libebur128 and not from the standard; this
    // comment used to credit 3342 with it.
    const int kRefWindowHops = 300;
    const int kRefObsStride  = 10;
    const double kRefAbsGate = 1.1724653045822981e-07;      // 10^((-70 + 0.691)/10), the BS.1770 gate
    const std::int64_t sub = std::max<std::int64_t> (1, (std::int64_t) std::lround (0.01 * fs));
    analysis::KWeightingFilter kw;
    kw.prepare (fs, channels);
    std::vector<double> hops;
    std::vector<double> acc (nch, 0.0);
    std::int64_t inHop = 0;
    for (std::size_t i = 0; i < T; ++i)
    {
        for (std::size_t c = 0; c < nch; ++c)
        {
            const float x = present (c) && std::isfinite (p[c][i]) ? p[c][i] : 0.0f;
            const double y = kw.process ((int) c, (double) x);
            acc[c] += y * y;
        }
        if (++inHop == sub)
        {
            double ms = 0.0;
            for (std::size_t c = 0; c < nch; ++c) { ms += acc[c] / (double) sub; acc[c] = 0.0; }
            hops.push_back (ms);
            kw.flushDenormals();
            inHop = 0;
        }
    }
    std::vector<double> obs;
    for (std::size_t h = (std::size_t) kRefWindowHops; h <= hops.size(); h += (std::size_t) kRefObsStride)
    {
        double s = 0.0;
        for (std::size_t j = h - (std::size_t) kRefWindowHops; j < h; ++j) s += hops[j];
        obs.push_back (s / (double) kRefWindowHops);
    }
    r.observations = (std::int64_t) obs.size();
    r.obsSeries = obs;

    const auto lufs = [] (double ms) { return ms > 1e-12 ? -0.691 + 10.0 * std::log10 (ms) : -120.0; };
    std::vector<double> gatedL;
    double absSum = 0.0;
    for (double e : obs) if (e >= kRefAbsGate) { gatedL.push_back (lufs (e)); absSum += e; }
    r.gated = (std::int64_t) gatedL.size();
    const auto pct = [] (std::vector<double> v, double q)
    {
        std::sort (v.begin(), v.end());
        const int rank = (int) ((double) ((int) v.size() - 1) * q + 0.5);
        return v[(std::size_t) rank];
    };
    if (gatedL.size() >= 2)
    {
        r.statsValid = true;
        r.p10 = pct (gatedL, 0.10); r.p50 = pct (gatedL, 0.50); r.p95 = pct (gatedL, 0.95);
        const double relT = 0.01 * (absSum / (double) gatedL.size());
        std::vector<double> lraL;
        for (double e : obs) if (e >= kRefAbsGate && e >= relT) lraL.push_back (lufs (e));
        if (lraL.size() >= 2) r.lra = pct (lraL, 0.95) - pct (lraL, 0.10);
    }
    return r;
}

// ============================================================================================== groups
static void testStorageAndRefusals()
{
    test::group ("storageFor / prepare: every argument is binding, and the budget IS the allocation");
    {
        Params pm;
        const auto st = PR::storageFor (kFs, 512, 2, pm);
        test::ok (st.ok, "a sane configuration is accepted");
        test::ok (st.subHopSamples == 480, "the 10 ms sub-hop at 48 kHz is 480 samples");
        test::ok (st.tailSamples == 4800, "a 100 ms tail window at 48 kHz is 4800 frames");
        test::ok (st.scratchFloats == 512u * 2u, "the scratch is maxBlock x channels and nothing else");
        test::ok (st.bytes() > 0, "the published budget is not zero");

        test::ok (! PR::storageFor (0.0, 512, 2, pm).ok,            "a zero rate is refused");
        test::ok (! PR::storageFor (std::numeric_limits<double>::quiet_NaN(), 512, 2, pm).ok, "a NaN rate is refused");
        // P51 — the floor. This class holds the K-weighting shelf, whose pole leaves the unit circle below 3364 Hz,
        // and it took 1000 Hz until the core's floor replaced its own copy.
        test::ok (PR::kMinSampleRate == 8000.0 && PR::storageFor (8000.0, 512, 2, pm).ok,
                  "the floor is 8000 Hz, and 8000 itself is accepted");
        test::ok (! PR::storageFor (std::nextafter (8000.0, 0.0), 512, 2, pm).ok && ! PR::storageFor (3300.0, 512, 2, pm).ok
                  && ! PR::storageFor (1000.0, 512, 2, pm).ok && ! PR::storageFor (44.1, 512, 2, pm).ok,
                  "one ulp under it, 3300 (the shelf past Nyquist), 1000 (the old floor) and 44.1 are refused");
        {
            PR q;
            test::ok (q.prepare (kFs, 512, 2), "PRECONDITION: a prepared report");
            test::ok (! q.prepare (std::nextafter (8000.0, 0.0), 512, 2) && ! q.isPrepared(),
                      "and a refused rate disarms it (law 11b) — it does not stay on its 48 kHz build");
        }
        test::ok (! PR::storageFor (kFs, 0, 2, pm).ok,              "maxBlock 0 is refused");
        test::ok (! PR::storageFor (kFs, 512, 0, pm).ok,            "0 channels is refused");
        test::ok (! PR::storageFor (kFs, 512, core::kMaxChannels + 1, pm).ok, "too many channels is refused");

        Params bad = pm; bad.infraLowHz = 0.5;
        test::ok (! PR::storageFor (kFs, 512, 2, bad).ok, "an infra-low frequency Svf would CLAMP is refused, not clamped");
        bad = pm; bad.infraLowHz = 0.46 * kFs;
        test::ok (! PR::storageFor (kFs, 512, 2, bad).ok, "an infra-low frequency past 0.45 fs is refused");
        bad = pm; bad.tailWindowMs = 0.0;
        test::ok (! PR::storageFor (kFs, 512, 2, bad).ok, "a zero tail window is refused");
        bad = pm; bad.maxDurationSec = 0.0;
        test::ok (! PR::storageFor (kFs, 512, 2, bad).ok, "a zero capacity is refused");
        bad = pm; bad.silenceThresholdDb = 6.0;
        test::ok (! PR::storageFor (kFs, 512, 2, bad).ok, "a silence threshold above full scale is refused");

        // THE BUDGET IS THE FORMULA, re-derived here from the same published numbers: if bytes() ever stops
        // summing what storageFor() states, a caller budgeting memory reads a number the allocation does not
        // honour. The two sub-meters' own byte counts are pinned by their own suites; what is checked here is
        // that this one includes them rather than quietly dropping them.
        const std::uint64_t hand = (std::uint64_t) sizeof (double) * (std::uint64_t) st.tailRingEnergies
                                 + (std::uint64_t) sizeof (std::int32_t) * (std::uint64_t) st.tailRingCounts
                                 + (std::uint64_t) sizeof (std::int32_t) * (std::uint64_t) st.tailPendingCounts
                                 + (std::uint64_t) sizeof (double) * (std::uint64_t) st.tailPendingEnergies
                                 + (std::uint64_t) sizeof (float) * (std::uint64_t) st.scratchFloats
                                 + (std::uint64_t) sizeof (double) * (std::uint64_t) st.shortTermEntries
                                 + st.loudness.bytes() + st.truePeak.bytes();
        test::ok (st.bytes() == hand, "bytes() is the sum of the published parts, the sub-meters included");
        test::ok (st.loudness.bytes() > 0 && st.truePeak.bytes() > 0,
                  "…and neither sub-meter's budget is silently zero");

        // A FRESH object allocates; a prepared one asked for the same shape again allocates nothing more,
        // which is what "a prepared one keeps storage that still fits" has to mean.
        PR a;
        a.setParams (pm);
        const long long before = alloc::count.load();
        test::run (a.prepare (kFs, 512, 2));
        const long long first = alloc::count.load() - before;
        test::ok (first > 0, "a fresh prepare() allocates (" + std::to_string (first) + " requests)");
        // THE DELTA IS CAPTURED BEFORE test::ok IS CALLED, and that is not style. `test::ok` takes a
        // std::string, so the message allocates; the order in which a compiler evaluates the two arguments
        // is unspecified, and it DIFFERS — spelled inline, this check read 0 under Apple clang/libc++ and 1
        // under gcc 14.2/libstdc++ on deb, i.e. it failed on one row for a reason that has nothing to do
        // with prepare(). An isolated probe with no strings in it reads 0 on both, which is what is
        // actually being asserted here.
        const long long mid = alloc::count.load();
        test::run (a.prepare (kFs, 512, 2));
        const long long again = alloc::count.load() - mid;
        test::ok (again == 0, "the same prepare() again asks the heap for nothing (got "
                              + std::to_string (again) + ")");

        // AND A REFUSED prepare() LEAVES NOTHING HALF-CONFIGURED (law 11b). The object was prepared and
        // working; a parameter it must refuse has to put it back to unprepared rather than leave the old
        // configuration live under new parameters.
        Params broken = pm; broken.infraLowHz = 0.5;
        a.setParams (broken);
        test::ok (! a.prepare (kFs, 512, 2), "a parameter out of range is refused");
        test::ok (! a.isPrepared(), "…and the object is left UNPREPARED, not half-configured");
        const auto tonePlanes = tone (2, 256, 1000.0, 0.25);
        const auto tonePtrs = ptrs (tonePlanes, 0);
        test::ok (! a.process (tonePtrs.data(), 2, 256), "…so process() refuses until a prepare() succeeds");
        a.setParams (pm);
        test::run (a.prepare (kFs, 512, 2));
        test::ok (a.isPrepared(), "a good prepare() after a refused one works");
    }

    test::group ("law 11: the call is a request against a prepared capacity");
    {
        const auto p = tone (2, 4096, 1000.0, 0.25);
        const auto v = ptrs (p, 0);
        PR a;
        test::ok (! a.process (v.data(), 2, 16), "process() before prepare() is refused");
        test::run (a.prepare (kFs, 256, 2));
        test::ok (! a.process (v.data(), 2, -1), "a negative n is refused");
        test::ok (! a.process (v.data(), -1, 16), "a negative width is refused");
        test::ok (! a.process (v.data(), 3, 16), "a width past the prepared one is refused");
        test::ok (! a.process (nullptr, 2, 16), "a null plane array is refused");
        const float* holed[2] { p[0].data(), nullptr };
        test::ok (! a.process (holed, 2, 16), "a null plane is refused");
        test::ok (a.process (v.data(), 2, 0) && a.samplesProcessed() == 0, "n == 0 is the one true no-op");
        test::run (a.process (v.data(), 2, 4096));
        test::ok (a.samplesProcessed() == 4096, "a call 16x longer than maxBlock is consumed IN FULL");
        a.finish();
        test::ok (! a.process (v.data(), 2, 16), "process() after finish() is refused");
        const auto first = a.report();
        a.finish();
        std::string why;
        test::ok (sameReport (first, a.report(), why), "finish() is idempotent (" + why + ")");
        a.reset();
        test::ok (a.samplesProcessed() == 0 && ! a.isFinished(), "reset() re-anchors the clock and unfreezes");
        test::run (a.process (v.data(), 2, 4096));
        a.finish();
        test::ok (sameReport (first, a.report(), why), "reset() then the same stream gives the same report (" + why + ")");
    }
}

static void testClosedForms()
{
    test::group ("closed forms: a pure sine's crest factor is 20 log10 sqrt(2) = 3.0103 dB");
    {
        // 480 samples per period at 100 Hz / 48 kHz — a whole number of periods, so the mean is exactly 0
        // in theory and the peak is reached.
        const auto p = tone (2, 48000, 100.0, 0.5);
        const auto o = run (p, kFs, 1024, { 1024 });
        test::ok (o.prepared && o.accepted, "the sine ran");
        test::approx (o.R.crestFactorDb.value, 3.0103, 0.01, "crest factor of a sine");
        test::approx (o.R.channel[0].rms.value, 0.5 / std::sqrt (2.0), 1e-4, "RMS of a 0.5 sine is 0.5/sqrt2");
        test::approx (o.R.channel[0].samplePeak.value, 0.5, 1e-3, "the sample peak is the amplitude");
        test::approx (o.R.channel[0].dcOffset.value, 0.0, 1e-6, "a whole number of periods has no DC");
        test::approx (o.R.stereoCorrelation.value, 1.0, 1e-12, "L == R reads correlation +1");
        test::ok (o.R.stereoBitIdenticalFrames == 48000, "every frame is bit-identical L to R");
        test::approx (o.R.stereoSideToMidRatio.value, 0.0, 1e-12, "no side energy in a mono-in-stereo file");
    }

    test::group ("closed forms: the infra-low fraction is the LR4's OWN |H|^2, not a brick wall");
    {
        // |LP_LR4(f)|^2 = 1/(1+(f/fc)^4)^2 at fc = 30 Hz. This is the check that dies if the filter is ever
        // swapped for a 4th-order Butterworth, whose 1/(1+(f/fc)^8) reads 9.10 % at 40 Hz against LR4's
        // 5.78 % — and the same mistake reads 50 % at the crossover where LR4 reads 25 %.
        // The tolerance is PURELY RELATIVE. An absolute floor beside it would have made the high-frequency
        // cells vacuous: at 100 Hz the steady state is 6.46e-5, so a +5e-4 floor accepts every answer
        // between "nothing" and eight times too much, and the cell tests the tolerance instead of the
        // filter. Which is also why the table STOPS at 60 Hz — see the note under it.
        struct Cell { double hz, want; };
        for (const Cell& cell : { Cell { 20.0, 0.697311 }, Cell { 30.0, 0.250000 },
                                  Cell { 40.0, 0.057771 }, Cell { 60.0, 0.003460 } })
        {
            const auto p = tone (1, (std::size_t) (kFs * 8.0), cell.hz, 0.5);
            const auto o = run (p, kFs, 4096, { 4096 });
            const double got = o.R.channel[0].infraLowFraction.value;
            test::ok (o.R.channel[0].infraLowFraction.valid, "the fraction is valid at "
                                                             + std::to_string ((int) cell.hz) + " Hz");
            test::ok (std::fabs (got - cell.want) <= 0.05 * cell.want,
                      "a " + std::to_string ((int) cell.hz) + " Hz tone reads |LR4|^2 = "
                      + std::to_string (cell.want) + " to 5 % (got " + std::to_string (got) + ")");
        }
        // WHY THE TABLE STOPS AT 60 Hz, and it is a property of the measurement rather than a gap in the
        // test: the ratio is Σ lp4(x)² / Σ x² over a FINITE window, and the filter's warm-up from reset()
        // is inside it. The warm-up contributes a roughly fixed amount of low-band energy, so it is
        // negligible against a large steady state and dominant against a small one. Measured on an 8 s
        // fixture: 100 Hz reads 1.55e-4 where the steady state is 6.46e-5 — 2.4x, all of it transient. A
        // cell up there would be pinning the onset, not |H|². The header says the same thing in words.
        {
            const auto p = tone (1, (std::size_t) (kFs * 8.0), 100.0, 0.5);
            const auto o = run (p, kFs, 4096, { 4096 });
            const double got = o.R.channel[0].infraLowFraction.value;
            const double steady = 1.0 / ((1.0 + std::pow (100.0 / 30.0, 4.0)) * (1.0 + std::pow (100.0 / 30.0, 4.0)));
            test::ok (got > steady && got < 20.0 * steady,
                      "at 100 Hz the reading is ABOVE the steady state and within a decade of it — the "
                      "warm-up, named rather than tuned away (got " + std::to_string (got)
                      + ", steady " + std::to_string (steady) + ")");
        }
    }

    test::group ("closed forms: EBU Tech 3341 case 1 — a 1 kHz sine at -23 dBFS reads -23.0 LUFS");
    {
        const double amp = std::pow (10.0, -23.0 / 20.0) * std::sqrt (2.0);   // -23 dBFS RMS
        const auto p = tone (2, (std::size_t) (kFs * 10.0), 1000.0, amp);
        const auto o = run (p, kFs, 2048, { 2048 });
        test::ok (o.R.integratedLufs.valid, "integrated loudness is a measurement, not the -120 sentinel");
        test::approx (o.R.integratedLufs.value, -23.0 + 10.0 * std::log10 (2.0), 0.2,
                      "two identical channels at -23 dBFS each read -19.99 LUFS (the channel sum)");
        test::ok (o.R.truePeakDbtp.valid, "the reference true peak is valid");
        test::approx (o.R.truePeakDbtp.value, 20.0 * std::log10 (amp), 0.15, "the true peak is the amplitude");
        test::ok (o.R.plrDb.valid, "PLR is valid when both its terms are");
        test::approx (o.R.plrDb.value, o.R.truePeakDbtp.value - o.R.integratedLufs.value, 1e-12,
                      "PLR is exactly the difference of the two published numbers");
    }

    test::group ("closed forms: a tail window built to carry 4x the programme's mean square reads 4");
    {
        // 2 s at 0.1, then exactly 100 ms at 0.2 — mean square 0.01 over 1.9 s and 0.04 over 0.1 s, so the
        // programme mean square is (1.9*0.01 + 0.1*0.04)/2.0 and the ratio is 0.04 over that.
        const std::size_t body = (std::size_t) (kFs * 1.9), tailN = (std::size_t) (kFs * 0.1);
        Planes p (1, std::vector<float> (body + tailN, 0.1f));
        for (std::size_t i = body; i < body + tailN; ++i) p[0][i] = 0.2f;
        const auto o = run (p, kFs, 1000, { 1000 });
        const double meanSq = (1.9 * 0.01 + 0.1 * 0.04) / 2.0;
        test::ok (o.R.tailEnergyRatio.valid, "the tail ratio is valid on a programme longer than the window");
        test::approx (o.R.tailEnergyRatio.value, 0.04 / meanSq, 1e-6, "the tail ratio is 0.04 / the mean square");
        test::approx (o.R.channel[0].lastSample.value, 0.2, 1e-7, "the last sample is published as delivered");
    }

    test::group ("closed forms: the stereo relations");
    {
        const std::size_t n = 24000;
        Planes p (2, std::vector<float> (n));
        for (std::size_t i = 0; i < n; ++i)
        {
            const double s = 0.4 * std::sin (kTwoPi * 220.0 * (double) i / kFs);
            p[0][i] = (float) s;
            p[1][i] = (float) (0.5 * s);                            // R is L at -6.02 dB
        }
        const auto o = run (p, kFs, 777, { 777 });
        test::approx (o.R.stereoBalanceDb.value, 10.0 * std::log10 (4.0), 1e-6,
                      "R at half of L reads +6.0206 dB of balance");
        test::approx (o.R.stereoCorrelation.value, 1.0, 1e-9, "a scaled copy is still fully correlated");
        test::ok (o.R.stereoValueEqualFrames < o.R.stereoFrames, "the channels are not value-identical");

        Planes q (2, std::vector<float> (n));
        for (std::size_t i = 0; i < n; ++i)
        {
            const double s = 0.4 * std::sin (kTwoPi * 220.0 * (double) i / kFs);
            q[0][i] = (float) s; q[1][i] = (float) -s;
        }
        const auto o2 = run (q, kFs, 777, { 777 });
        test::approx (o2.R.stereoCorrelation.value, -1.0, 1e-9, "an inverted channel reads correlation -1");
        test::ok (o2.R.stereoSideToMidRatio.valid == false || o2.R.stereoSideToMidRatio.value > 1e6,
                  "all side and no mid: the ratio is either refused or enormous");
    }

    test::group ("closed forms: LRA over a programme that steps 10 LU");
    {
        // 20 s at -23 dBFS then 20 s at -33 dBFS: the gated short-term distribution has two clusters 10 LU
        // apart, so P95 - P10 is 10 LU to within the relative gate and the 0.1 LU bin.
        const double a1 = std::pow (10.0, -23.0 / 20.0) * std::sqrt (2.0);
        const double a2 = std::pow (10.0, -33.0 / 20.0) * std::sqrt (2.0);
        const std::size_t half = (std::size_t) (kFs * 20.0);
        Planes p (2, std::vector<float> (half * 2));
        for (std::size_t i = 0; i < half * 2; ++i)
        {
            const double amp = i < half ? a1 : a2;
            const float v = (float) (amp * std::sin (kTwoPi * 1000.0 * (double) i / kFs));
            p[0][i] = v; p[1][i] = v;
        }
        Params pm; pm.maxDurationSec = 120.0;
        const auto o = run (p, kFs, 4096, { 4096 }, pm);
        test::ok (o.R.lraLu.valid, "LRA is valid with this many gated observations");
        test::approx (o.R.lraLu.value, 10.0, 0.6, "a 10 LU step reads about 10 LU of range");
        test::ok (o.R.shortTermGatedObservations > 20, "the gated observation count is published ("
                                                       + std::to_string (o.R.shortTermGatedObservations) + ")");
        test::approx (o.R.shortTermSpreadLu.value, 10.0, 0.6, "the ungated-by-relative spread agrees here");

        // AND THE SHIPPED METER AGREES. This class computes EBU Tech 3342 from its own series so that the
        // gated count is exactly known; the number itself must still be the meter's.
        analysis::LoudnessMeter lm;
        test::run (lm.prepare (kFs, 2, 120.0));
        for (std::size_t at = 0; at < p[0].size(); at += 4096)
        {
            const int n = (int) std::min<std::size_t> (4096, p[0].size() - at);
            const float* ch[2] { p[0].data() + at, p[1].data() + at };
            test::run (lm.process (ch, 2, n));
        }
        test::approx (o.R.lraLu.value, lm.loudnessRangeLu(), 0.3,
                      "on THIS programme the LRA here and LoudnessMeter::loudnessRangeLu() are the same number "
                      "(long steady steps, where the sampling cadence cannot matter — see the group below)");

        // …AND THE CHECK ABOVE, ALONE, IS NARROWER THAN THE SENTENCE IT STANDS UNDER. Its programme is long
        // steady steps, where every window that is not straddling a step reads the same value and the sampling
        // CADENCE cannot matter. That is exactly why it kept passing through the window in which the two
        // classes had different cadences: K12 gave LoudnessMeter the 10 Hz EBU Tech 3342 §3.1 requires while
        // this class still took one observation a second, and on a square envelope whose states last exactly
        // the 3 s window they read 20.0 LU against 9.5 — 10.5 LU apart, both labelled Tech 3342, with this
        // group's ancestor measuring the gap and calling it an open defect.
        //
        // v0.43 CLOSED IT: `kObservationHops` is 10, both series sample at 10 Hz, and what was a gap detector
        // is now the guard for the equality. Held on the fixture that can SEE a cadence, at the same 0.3 LU
        // the steady one uses — that tolerance is summation order (this class oldest-first, the meter
        // newest-first), and nothing else should fit inside it.
        {
            const double fs = kFs;
            const int frames = (int) (fs * 12.0);
            const std::size_t state = (std::size_t) (fs * 3.0);
            std::vector<std::vector<float>> q (2, std::vector<float> ((std::size_t) frames, 0.0f));
            for (int i = 0; i < frames; ++i)
            {
                const double s = std::sin (2.0 * core::kPi * 1000.0 * (double) i / fs);
                const double g = (((std::size_t) i / state) % 2 == 0) ? 0.25 : 0.025;   // 20 dB apart
                q[0][(std::size_t) i] = q[1][(std::size_t) i] = (float) (g * s);
            }
            Params wpm; wpm.maxDurationSec = 120.0;
            const auto w = run (q, fs, 4096, { 4096 }, wpm);
            analysis::LoudnessMeter ref;
            test::run (ref.prepareForSamples (fs, 2, (double) frames));
            const float* rp[2] { q[0].data(), q[1].data() };
            test::run (ref.process (rp, 2, frames));
            test::ok (w.R.lraLu.valid, "PRECONDITION: the report publishes an LRA for this programme");
            // The fixture is live as a cadence probe only if a 1 Hz grid would answer differently — 20.0 LU,
            // which is also the envelope's full swing and so the ceiling. Asserting the pair is BELOW that
            // keeps the precondition honest without pinning a number no shipped build can produce.
            test::ok (w.R.lraLu.value < 12.0 && ref.loudnessRangeLu() < 12.0,
                      "PRECONDITION: neither reads the 20 LU swing a 1 Hz grid lands on (report "
                      + std::to_string (w.R.lraLu.value) + ", meter " + std::to_string (ref.loudnessRangeLu()) + ")");
            test::approx (w.R.lraLu.value, ref.loudnessRangeLu(), 0.3,
                          "on an envelope that CAN tell two cadences apart, the LRA here and "
                          "LoudnessMeter::loudnessRangeLu() are still the same number");
        }
    }
}

static void testNegatives()
{
    test::group ("the histogram's edges: the top clamp, and the absolute gate's own constant");
    {
        // The gate is a LITERAL in the header rather than a std::pow call, because it decides inclusion and
        // a libm that disagrees by an ulp would move a borderline observation in on one row and out on
        // another. Pin it against the formula it is derived from.
        const double formula = std::pow (10.0, (-70.0 + 0.691) / 10.0);
        test::ok (std::fabs (PR::kAbsoluteGateEnergy - formula) <= 1.0e-22,
                  "kAbsoluteGateEnergy is 10^((-70 + 0.691)/10) to within an ulp");

        // A programme loud enough to leave the histogram's +30 LUFS top: amplitude 60 over two channels
        // reads -0.691 + 10 log10(2 * 60^2/2) = +34.87 LUFS, so every observation lands in the last bin and
        // the percentile is that bin's lower bound, -70 + 999*0.1 = +29.9. Nothing else in this suite ever
        // reaches an edge of the histogram, and a crew round showed the clamp unpinned.
        Params pm; pm.maxDurationSec = 60.0;
        const auto o = run (tone (2, (std::size_t) (kFs * 20.0), 1000.0, 60.0), kFs, 4096, { 4096 }, pm);
        test::ok (o.R.shortTermP95.valid, "an enormous programme still produces percentiles");
        test::approx (o.R.shortTermP95.value, 29.9, 1e-9, "P95 clamps to the TOP bin's lower bound, +29.9");
        test::approx (o.R.shortTermP10.value, 29.9, 1e-9, "…and so does P10, the whole distribution being there");
        test::approx (o.R.lraLu.value, 0.0, 1e-9, "…so the range across one bin is 0 LU");
    }

    test::group ("the two silence predicates transition at DIFFERENT samples, and the trace says which is which");
    {
        // 0.2 s of exact zeros, 0.2 s at 1e-7 (non-zero, but far under the -96 dBFS = 1.585e-5 threshold),
        // 0.2 s of real signal. The DIGITAL predicate turns on at 9600; the THRESHOLDED one not until
        // 19200. So the two edges are distinguishable by position, which is what pins the predicate ids —
        // a crew round swapped them and every other check in this file stayed green.
        const std::size_t q = (std::size_t) (kFs * 0.2);
        Planes p (1, std::vector<float> (q * 3, 0.0f));
        for (std::size_t i = q; i < 2 * q; ++i) p[0][i] = 1.0e-7f;
        for (std::size_t i = 2 * q; i < 3 * q; ++i)
            p[0][i] = (float) (0.3 * std::cos (kTwoPi * 440.0 * (double) (i - 2 * q) / kFs));   // starts AT 0.3
        const auto o = run (p, kFs, 512, { 512 });
        bool digitalOn = false, thresholdOn = false;
        for (const Event& e : o.trace)
            if (e.kind == analysis::ProgrammeTraceKind::SilenceEdge && e.b > 0.5)
            {
                if (e.a > 0.5 && e.begin == (std::int64_t) q)         digitalOn = true;    // predicate id 1
                if (e.a < 0.5 && e.begin == (std::int64_t) (2 * q))   thresholdOn = true;  // predicate id 0
            }
        test::ok (digitalOn, "the DIGITAL predicate (id 1) turns on where the zeros end, at 9600");
        test::ok (thresholdOn, "the THRESHOLDED predicate (id 0) turns on where the level clears -96 dBFS, at 19200");
        test::ok (o.R.leadingDigitalSilenceSamples == (std::int64_t) q, "leading digital silence is the first fifth");
        test::ok (o.R.leadingSilenceSamples.value == (std::int64_t) (2 * q),
                  "leading thresholded silence is the first TWO fifths");
    }

    test::group ("NEGATIVE 1: digital-silence padding must not lose the tail, or move ANY measurement");
    {
        const std::size_t body = (std::size_t) (kFs * 1.9), tailN = (std::size_t) (kFs * 0.1);
        Planes bare (2, std::vector<float> (body + tailN, 0.1f));
        for (std::size_t i = body; i < body + tailN; ++i) { bare[0][i] = 0.2f; bare[1][i] = 0.2f; }
        const auto padded = pad (bare, (std::size_t) (kFs * 5.0), (std::size_t) (kFs * 5.0));

        const auto a = run (bare, kFs, 1024, { 1024 });
        const auto b = run (padded, kFs, 1024, { 1024 });

        test::ok (a.R.tailEnergyRatio.valid && b.R.tailEnergyRatio.valid, "both tails are measurable");
        test::ok (b64 (a.R.tailEnergyRatio.value) == b64 (b.R.tailEnergyRatio.value),
                  "the tail ratio is BIT-identical with 5 s of zeros at each end");
        test::ok (b64 (a.R.rms.value) == b64 (b.R.rms.value), "the RMS is bit-identical");
        test::ok (b64 (a.R.crestFactorDb.value) == b64 (b.R.crestFactorDb.value), "the crest factor is bit-identical");
        test::ok (b64 (a.R.channel[0].dcOffset.value) == b64 (b.R.channel[0].dcOffset.value), "the DC is bit-identical");
        test::ok (b64 (a.R.programmeMeanSquare.value) == b64 (b.R.programmeMeanSquare.value),
                  "the programme mean square is bit-identical");
        test::ok (a.R.programmeSpanSamples == b.R.programmeSpanSamples, "the span is the same length");
        test::ok (b.R.leadingDigitalSilenceSamples == (std::int64_t) (kFs * 5.0)
                  && b.R.trailingDigitalSilenceSamples == (std::int64_t) (kFs * 5.0),
                  "…and the padding is REPORTED, not absorbed");
        test::ok (a.R.leadingDigitalSilenceSamples == 0 && a.R.trailingDigitalSilenceSamples == 0,
                  "the unpadded file reports no digital silence");

        // AND NOW EVERY FIELD, against an explicit ALLOW-LIST of what may move — so the claim is measured
        // rather than asserted about the handful of fields somebody thought to check, and so a field added
        // later that turns out to be padding-sensitive fails here instead of quietly widening the claim.
        //
        // WHAT MAY MOVE, AND WHY IT IS NOT A DEFECT:
        //   lastSample        — the last sample OF THE FILE, which is the pad's zero. That is its
        //                       definition; `lastSignalSample` is the padding-invariant companion and is
        //                       NOT on this list.
        //   the loudness family — integrated, PLR, LRA and the short-term percentiles. BS.1770 anchors its
        //                       400 ms gating blocks and 3 s windows at the START OF THE STREAM, so a pad
        //                       shifts which windows straddle the music's edges even when the pad is a
        //                       whole number of sub-hops (5 s is exactly 500). Anchoring the sub-meters at
        //                       the first active frame instead would make these invariant too, and was not
        //                       done: it would put this report's loudness numbers at odds with
        //                       `fcore_measure lufs`, ffmpeg's ebur128 and every other EBU tool on the same
        //                       file, which is a worse trade than a padding-sensitive LUFS.
        // Everything else — DC, peak, RMS, crest, infra-low, the programme mean square, the tail ratio, the
        // stereo relations, and the reference TRUE PEAK (a maximum, to which zeros add nothing) — is
        // bit-identical.
        {
            static const char* kMayMove[] { "lastSample", "integratedLufs", "plrDb", "lraLu",
                                            "shortTermP10", "shortTermP50", "shortTermP95",
                                            "shortTermSpreadLu" };
            const auto va = valuesOf (a.R), vb = valuesOf (b.R);
            test::ok (va.size() == vb.size(), "both reports enumerate the same fields");
            std::string moved;
            int unexpected = 0, invariant = 0;
            for (std::size_t i = 0; i < va.size() && i < vb.size(); ++i)
            {
                bool allowed = false;
                for (const char* n : kMayMove) if (va[i].name == n) allowed = true;
                const bool same = b64 (va[i].v.value) == b64 (vb[i].v.value)
                               && va[i].v.valid == vb[i].v.valid && va[i].v.reason == vb[i].v.reason;
                if (allowed) continue;
                if (same) ++invariant;
                else { ++unexpected; moved += " " + va[i].name + "[" + std::to_string (va[i].ch) + "]"; }
            }
            test::ok (unexpected == 0, "every field outside the allow-list is BIT-identical under padding;"
                                       " moved:" + moved);
            test::ok (invariant >= 15, "…and that is most of the report, not a corner of it ("
                                       + std::to_string (invariant) + " fields)");
        }
    }

    test::group ("NEGATIVE 2: asymmetric content has a large mean on a short window and none on the programme");
    {
        // Half-wave-ish asymmetry whose polarity flips every 0.5 s: any 100 ms window has a big mean, the
        // programme has none. A DC measured over a window instead of the programme would read ~0.15.
        const std::size_t n = (std::size_t) (kFs * 4.0), flip = (std::size_t) (kFs * 0.5);
        Planes p (1, std::vector<float> (n));
        for (std::size_t i = 0; i < n; ++i)
        {
            const double s = std::sin (kTwoPi * 200.0 * (double) i / kFs);
            const double sat = s > 0.0 ? 0.6 * s : 0.15 * s;                 // one-sided saturation
            p[0][i] = (float) (((i / flip) % 2u == 0u) ? sat : -sat);        // …with alternating polarity
        }
        const auto o = run (p, kFs, 640, { 640 });

        double windowMean = 0.0;
        const std::size_t w = (std::size_t) (kFs * 0.1);
        for (std::size_t i = 0; i < w; ++i) windowMean += (double) p[0][i];
        windowMean /= (double) w;
        test::ok (std::fabs (windowMean) > 0.05, "precondition — a 100 ms window's mean is large ("
                                                 + std::to_string (windowMean) + ")");
        test::ok (std::fabs (o.R.channel[0].dcOffset.value) < 1e-3,
                  "the programme's DC is ~0 (" + std::to_string (o.R.channel[0].dcOffset.value) + ")");
    }

    test::group ("NEGATIVE 3: a mono programme gets no polarity, and no fictitious correlation");
    {
        const auto p = tone (1, 48000, 440.0, 0.3);
        const auto o = run (p, kFs, 512, { 512 });
        test::ok (! o.R.stereoCorrelation.valid && o.R.stereoCorrelation.reason == Reason::MonoProgramme,
                  "a one-channel programme refuses the correlation with MonoProgramme");
        test::ok (b64 (o.R.stereoCorrelation.value) == b64 (0.0), "…and its canonical value is +0.0, not +1.0");
        test::ok (! o.R.stereoBalanceDb.valid && ! o.R.stereoSideToMidRatio.valid,
                  "balance and side/mid are refused too");
        test::ok (o.R.stereoFrames == 0, "no stereo frame was counted");

        // Prepared for stereo and then fed one channel for the whole programme: the relation still does not
        // exist, and the reason says which of the two situations it was.
        const auto o2 = run (p, kFs, 512, { 512 }, Params {}, 2, 1);
        test::ok (! o2.R.stereoCorrelation.valid && o2.R.stereoCorrelation.reason == Reason::NoStereoFrames,
                  "prepared stereo but never fed stereo reads NoStereoFrames");
        test::ok (o2.R.channel[1].absentSamples == 48000, "…and channel 1's absence is counted");
    }
}

static void testInvariance()
{
    test::group ("law 8a: the report and the WHOLE TRACE are bit-identical under any re-splitting");
    {
        // 4 s of non-stationary stereo with digital-silence padding, so the fixture crosses every boundary
        // the schedule has: sub-hop ends, StateGrid ticks, both silence predicates, and one short-term
        // observation at 3 s plus one at 4 s.
        auto p = pad (programme (2, (std::size_t) (kFs * 3.6), 20260914u),
                      (std::size_t) (kFs * 0.2), (std::size_t) (kFs * 0.2));
        const int maxBlock = 1024;
        const int H = 480;                                   // the 10 ms sub-hop at 48 kHz
        const int G = core::StateGrid::kPeriod;              // 64
        const int W = 300 * H;                               // the 3 s short-term window

        const auto base = run (p, kFs, maxBlock, { (int) p[0].size() });
        test::ok (base.prepared && base.accepted, "the whole-file call ran");
        test::ok (base.trace.size() > 1000, "the trace is substantial (" + std::to_string (base.trace.size()) + " events)");
        test::ok (base.traceOverflow == 0, "the trace buffer held every event");

        // THE CADENCE, ASSERTED. 3 s of sub-hops, then one observation per 100 ms HOP — ten sub-hops, the
        // >=10 Hz EBU Tech 3342 §3.1 requires: observation k closes at sample (300 + 10k) * subHopSamples.
        // Nothing else in the suite pins the stride, and a crew round showed a stride one off surviving every
        // other check in this file. (It read 100 — one second, libebur128's — until v0.43.)
        {
            const std::int64_t H64 = (std::int64_t) H;
            std::int64_t k = 0;
            for (const Event& e : base.trace)
                if (e.kind == analysis::ProgrammeTraceKind::ShortTerm)
                {
                    test::ok (e.index == k, "observation " + std::to_string (k) + " is numbered in order");
                    test::ok (e.end == (300 + 10 * k) * H64,
                              "observation " + std::to_string (k) + " closes at (300 + 10k) sub-hops = "
                              + std::to_string ((300 + 10 * k) * H64) + " (got " + std::to_string (e.end) + ")");
                    test::ok (e.end - e.begin == 300 * H64, "…and spans exactly 3 s of sub-hops");
                    ++k;
                }
            test::ok (k >= 2, "the fixture produced at least two observations (" + std::to_string (k) + ")");
        }

        const std::vector<std::vector<int>> splits {
            { 1 }, { 2 }, { 3 }, { 7 }, { 8191 },
            { G - 1 }, { G }, { G + 1 },
            { H - 1 }, { H }, { H + 1 },
            { W - 1 }, { W }, { W + 1 },
            { maxBlock }, { maxBlock + 1 }, { 5 * maxBlock + 3 },
            { 1, 480, 63, 1, 4096, 17 },                     // ragged, cycled
            { H, 1, H - 1, 2, G, 3 }
        };
        int differing = 0;
        for (const auto& s : splits)
        {
            const auto o = run (p, kFs, maxBlock, s);
            std::string why;
            const bool sameR = sameReport (base.R, o.R, why);
            const bool sameT = sameTrace (base.trace, o.trace, why);
            if (! sameR || ! sameT) ++differing;
            test::ok (sameR, "report bit-identical at split " + std::to_string (s[0]) + " (" + why + ")");
            test::ok (sameT, "trace bit-identical at split " + std::to_string (s[0]) + " (" + why + ")");
        }
        test::ok (differing == 0, "every slicing agrees");

        // …and boundaries placed exactly on the transitions, from a fixed seed.
        std::mt19937 rng (7u);
        for (int trial = 0; trial < 6; ++trial)
        {
            std::vector<int> s;
            for (int k = 0; k < 40; ++k)
            {
                static const int anchors[] { 1, 2, 63, 64, 65, 479, 480, 481, 1023, 1024, 1025, 9600 };
                s.push_back (anchors[rng() % (sizeof anchors / sizeof anchors[0])]);
            }
            const auto o = run (p, kFs, maxBlock, s);
            std::string why;
            test::ok (sameReport (base.R, o.R, why) && sameTrace (base.trace, o.trace, why),
                      "anchored random slicing " + std::to_string (trial) + " agrees (" + why + ")");
        }
    }

    test::group ("law 8a: maxBlock is scratch and nothing else — a different one reads the same report");
    {
        auto p = programme (2, (std::size_t) (kFs * 3.5), 4242u);
        const auto base = run (p, kFs, 64, { 997 });
        for (int mb : { 1, 2, 63, 65, 256, 1024, 4096, (int) p[0].size() })
        {
            const auto o = run (p, kFs, mb, { 997 });
            std::string why;
            test::ok (sameReport (base.R, o.R, why) && sameTrace (base.trace, o.trace, why),
                      "maxBlock " + std::to_string (mb) + " reads the same (" + why + ")");
        }
    }

    test::group ("the denormal-flush cadence is OBSERVABLE, so the gate for it is not theatre");
    {
        // A 2^-50 impulse then 4800 zeros. The LR4's float state decays through the 1e-15 flush threshold
        // during the silence, so WHERE the flush happens decides whether the low-band energy of those
        // sub-hops is exactly 0 or ~1e-30 — and this fixture is what makes the mutant that moves the flush
        // to the end of process() fail. It is checked here as a PRECONDITION on the fixture (the tiny tail
        // energy exists at all) so a future change that made it unobservable is caught as a failure and not
        // as a silent weakening.
        // THE POSITIVE ASSERTION ABOUT THE CADENCE, and it is an oracle rather than a precondition: the grid
        // flush zeroes the crossover's state the moment it falls under 1e-15, so the low-band energy of the
        // sub-hops after the impulse is EXACTLY zero — the decay is CUT. Move the flush to the end of
        // process() and a whole-file call never reaches a boundary inside itself, so the state rings on and
        // every one of those sub-hops carries ~1e-30 instead. That is the difference this line reads.
        for (int amplitudeCase = 0; amplitudeCase < 2; ++amplitudeCase)
        {
            // 0.5 s = 50 sub-hops of room: at 30 Hz the LR4 rings down at about 133 nepers a second, so a
            // 1e-9 impulse needs ~104 ms to reach the 1e-15 flush threshold on its own. A 100 ms fixture is
            // one sub-hop short of showing the cut and would have read as "the gate does not fire".
            Planes q (1, std::vector<float> (24000, 0.0f));
            q[0][0] = amplitudeCase == 0 ? std::ldexp (1.0f, -50) : 1.0e-9f;
            const auto w  = run (q, kFs, 24000, { 24000 });     // maxBlock = the whole call: a flush moved
            const auto s1 = run (q, kFs, 24000, { 1 });         // to the end of process() never fires inside
            std::string w2;
            test::ok (sameReport (w.R, s1.R, w2) && sameTrace (w.trace, s1.trace, w2),
                      "the impulse fixture is split-invariant at amplitude case "
                      + std::to_string (amplitudeCase) + " (" + w2 + ")");

            int lastNonZeroHop = -1, firstZeroAfter = -1;
            double lowSeen = 0.0;
            for (const Event& e : w.trace)
                if (e.kind == analysis::ProgrammeTraceKind::SubHopChannel)
                {
                    if (e.b > 0.0) { lastNonZeroHop = (int) e.index; lowSeen = std::max (lowSeen, e.b); }
                    else if (lastNonZeroHop >= 0 && firstZeroAfter < 0) firstZeroAfter = (int) e.index;
                }
            test::ok (lastNonZeroHop >= 0, "the impulse put energy in the low band (case "
                                           + std::to_string (amplitudeCase) + ")");
            test::ok (firstZeroAfter > lastNonZeroHop,
                      "…and a later sub-hop reads EXACTLY zero: the grid flush cut the decay (case "
                      + std::to_string (amplitudeCase) + ", cut after hop "
                      + std::to_string (lastNonZeroHop) + ")");
            // WHERE the cut falls. A 2^-50 (8.88e-16) impulse puts the crossover's state UNDER the 1e-15
            // threshold at once, so the first grid boundary — sample 63, inside sub-hop 0 — zeroes it: the
            // last sub-hop carrying any low-band energy is 0, and that is arithmetic, not a measurement.
            // A 1e-9 impulse has 6 decades = 13.8 nepers to fall. The single-section estimate is
            // omega/2Q = 2*pi*30/1.414 = 133.29 nepers a second, i.e. 103.7 ms, and it is TOO SHORT: LR4 is
            // two cascaded sections with the SAME pole pair, so the envelope is t*exp(-sigma*t) rather than
            // exp(-sigma*t), and the second section lags the first. Measured: the last non-zero low-band
            // sample is 7103, i.e. 148 ms, in sub-hop 14. The bound
            // here is loose on purpose — the gate is not this number but the line above, which the mutant
            // fails outright by never cutting the decay at all (it reads hop 49, the end of the fixture).
            const int wantMax = amplitudeCase == 0 ? 0 : 20;
            test::ok (lastNonZeroHop <= wantMax,
                      "…at the hop the flush threshold puts it at, not wherever the call happened to end "
                      "(case " + std::to_string (amplitudeCase) + ": hop "
                      + std::to_string (lastNonZeroHop) + " <= " + std::to_string (wantMax) + ")");
            test::ok (lowSeen > 0.0 && lowSeen < 1e-12,
                      "…in the tiny range the 1e-15 flush threshold decides (" + std::to_string (lowSeen) + ")");
        }
    }
}

static void testAgainstReference()
{
    test::group ("the independent whole-file reference: a second program for the same numbers");
    {
        struct Case { const char* name; Planes p; int ch; };
        std::vector<Case> cases;
        cases.push_back (Case { "non-stationary stereo", programme (2, (std::size_t) (kFs * 6.0), 11u), 2 });
        cases.push_back (Case { "padded mono", pad (programme (1, (std::size_t) (kFs * 5.0), 12u),
                                                    12345u, 54321u), 1 });
        cases.push_back (Case { "quiet stereo", tone (2, (std::size_t) (kFs * 5.0), 63.0, 1e-4), 2 });
        cases.push_back (Case { "loud transient stereo", programme (2, (std::size_t) (kFs * 4.5), 13u), 2 });
        {
            // AN INTERIOR SILENCE GAP LONGER THAN THE TAIL WINDOW, then a burst SHORTER than it, then
            // trailing silence. The tail window therefore reaches back across the gap, and it is the only
            // fixture whose window is part signal and part interior silence. Two defects live exactly here:
            // the deferred silent frames' own sample counts (write them back as 0 and the window's
            // denominator shrinks and the tail mean square inflates), and the skip in
            // `tailWritten_ += pending - keep` (drop it and the window is read from the wrong ring
            // positions). Both survived every other check in this file.
            const std::size_t head = (std::size_t) (kFs * 0.3), gap = (std::size_t) (kFs * 0.3);
            const std::size_t burst = (std::size_t) (kFs * 0.05), tailPad = (std::size_t) (kFs * 0.4);
            Planes g (2, std::vector<float> (head + gap + burst + tailPad, 0.0f));
            for (std::size_t i = 0; i < head; ++i)
            {
                const float v = (float) (0.4 * std::sin (kTwoPi * 180.0 * (double) i / kFs));
                g[0][i] = v; g[1][i] = 0.7f * v;
            }
            for (std::size_t i = 0; i < burst; ++i)
            {
                const float v = (float) (0.25 * std::sin (kTwoPi * 300.0 * (double) i / kFs));
                g[0][head + gap + i] = v; g[1][head + gap + i] = -v;
            }
            cases.push_back (Case { "interior gap wider than the tail window", g, 2 });
        }

        Params pm;
        for (const auto& c : cases)
        {
            const auto o = run (c.p, kFs, 1024, { 1000, 1, 4096, 63 }, pm);
            const Ref r = reference (c.p, kFs, pm, c.ch);
            const std::string at = std::string (" [") + c.name + "]";

            test::ok (o.R.programmeSpanSamples == r.span, "the span agrees" + at);
            test::ok (o.R.leadingDigitalSilenceSamples == r.leadDigital, "leading digital silence agrees" + at);
            test::ok (o.R.trailingDigitalSilenceSamples == r.trailDigital, "trailing digital silence agrees" + at);
            test::ok (o.R.leadingSilenceSamples.value == r.leadThr, "leading thresholded silence agrees" + at);
            test::ok (o.R.trailingSilenceSamples.value == r.trailThr, "trailing thresholded silence agrees" + at);

            for (int ci = 0; ci < c.ch; ++ci)
            {
                const auto& P = o.R.channel[(std::size_t) ci];
                const std::string atc = at + " ch" + std::to_string (ci);
                // The DC of a symmetric programme is a near-perfect cancellation, so the two orders leave
                // different residues and a relative tolerance on it would be meaningless. What is
                // meaningful is that both are zero AGAINST THE SIGNAL: a billionth of the channel's RMS.
                test::ok (std::fabs (P.dcOffset.value - r.dc[(std::size_t) ci])
                          <= 1.0e-9 * std::max (r.rms[(std::size_t) ci], 1.0e-30),
                          "DC agrees to 1e-9 of the signal level" + atc);
                relClose (P.rms.value, r.rms[(std::size_t) ci], 1e-9, "RMS agrees" + atc);
                test::approx (P.samplePeak.value, r.peak[(std::size_t) ci], 0.0, "the peak agrees exactly" + atc);
                test::approx (P.crestFactorDb.value, r.crest[(std::size_t) ci], 1e-6, "the crest factor agrees" + atc);
                relClose (P.infraLowFraction.value, r.infra[(std::size_t) ci], 1e-6,
                          "the infra-low fraction agrees" + atc);
            }
            relClose (o.R.rms.value, r.rmsAll, 1e-9, "the programme RMS agrees" + at);
            test::approx (o.R.crestFactorDb.value, r.crestAll, 1e-6, "the programme crest factor agrees" + at);
            relClose (o.R.programmeMeanSquare.value, r.progMeanSquare, 1e-9, "the mean square agrees" + at);
            test::ok (o.R.tailEnergyRatio.valid == r.tailValid, "the tail's validity agrees" + at);
            if (r.tailValid) relClose (o.R.tailEnergyRatio.value, r.tailRatio, 1e-6, "the tail ratio agrees" + at);
            // The window's DENOMINATOR, not just its sum: every frame of the window has every prepared
            // channel present and finite here, so the count is exactly the window length times the width.
            if (r.tailValid)
                test::ok (o.R.tailWindowFinite == (std::int64_t) std::lround (kFs * pm.tailWindowMs / 1000.0)
                                                  * (std::int64_t) c.ch,
                          "the tail window counts every sample in it" + at
                          + " (got " + std::to_string (o.R.tailWindowFinite) + ")");

            if (c.ch >= 2)
            {
                test::approx (o.R.stereoCorrelation.value, r.corr, 0.0, "the correlation agrees EXACTLY" + at);
                relClose (o.R.stereoBalanceDb.value, r.balanceDb, 1e-9, "the balance agrees" + at);
                relClose (o.R.stereoSideToMidRatio.value, r.sideToMid, 1e-9, "side/mid agrees" + at);
            }

            test::ok (o.R.shortTermObservations == r.observations, "the observation count agrees" + at);
            // THE SERIES ITSELF, not only its percentiles. The reference already builds the short-term
            // observations independently and summed them in the same (oldest-first) order, so these are
            // comparable directly — the cheapest second-program gate on the series available, and it was
            // sitting unused until a crew round pointed at it.
            {
                std::size_t k = 0;
                for (const Event& e : o.trace)
                    if (e.kind == analysis::ProgrammeTraceKind::ShortTerm && k < r.obsSeries.size())
                    {
                        relClose (e.a, r.obsSeries[k], 1e-9,
                                  "short-term observation " + std::to_string (k) + " agrees with the "
                                  "independent series" + at);
                        ++k;
                    }
                test::ok (k == r.obsSeries.size(), "…for every observation the reference produced" + at);
            }
            test::ok (o.R.shortTermGatedObservations == r.gated, "the gated count agrees" + at);
            if (r.statsValid && o.R.shortTermP50.valid)
            {
                // The reference sorts and returns an order statistic; this class returns a histogram bin's
                // lower bound, so they agree to within one 0.1 LU bin and no closer.
                test::approx (o.R.shortTermP10.value, r.p10, 0.1001, "P10 agrees to one bin" + at);
                test::approx (o.R.shortTermP50.value, r.p50, 0.1001, "P50 agrees to one bin" + at);
                test::approx (o.R.shortTermP95.value, r.p95, 0.1001, "P95 agrees to one bin" + at);
                if (o.R.lraLu.valid) test::approx (o.R.lraLu.value, r.lra, 0.2001, "LRA agrees to two bins" + at);
            }
        }
    }
}

static void testEdges()
{
    test::group ("lengths: 0, 1, W-1, W, W+H-1, W+H, and nothing in between is special");
    {
        const int H = 480, W = 300 * H;
        for (std::size_t T : { (std::size_t) 0, (std::size_t) 1, (std::size_t) (W - 1), (std::size_t) W,
                               (std::size_t) (W + H - 1), (std::size_t) (W + H) })
        {
            Planes p = T == 0 ? Planes (2, std::vector<float> {}) : programme (2, T, 99u);
            const auto o = run (p, kFs, 1024, { 333 });
            test::ok (o.prepared && o.accepted, "T = " + std::to_string (T) + " ran");
            test::ok (o.R.totalSamples == (std::int64_t) T, "T = " + std::to_string (T) + " counted every sample");
            const std::int64_t closed = (std::int64_t) T / H;
            test::ok (o.R.uncoveredSubHopSamples == (std::int64_t) T - closed * H,
                      "the sub-hop remainder is NAMED at T = " + std::to_string (T));
            if (T < (std::size_t) W)
                test::ok (! o.R.lraLu.valid, "shorter than the 3 s window: LRA refuses at T = " + std::to_string (T));
        }
    }

    test::group ("an empty and an all-silent programme say cannot-say, and never publish a comfortable zero");
    {
        const auto e = run (Planes (2, std::vector<float> {}), kFs, 256, { 256 });
        test::ok (e.prepared, "an empty stream still prepares");
        test::ok (e.R.totalSamples == 0, "nothing was processed");
        test::ok (! e.R.rms.valid && e.R.rms.reason == Reason::SilentProgramme, "the RMS refuses");
        test::ok (! e.R.tailEnergyRatio.valid, "the tail refuses");
        test::ok (! e.R.integratedLufs.valid
                  && e.R.integratedLufs.reason == Reason::ShorterThanLoudnessWindow,
                  "integrated loudness refuses instead of publishing -120 — and names the right cause: no "
                  "400 ms gating block EXISTED, which LoudnessMeter's -120 does not distinguish from "
                  "'every block was gated out' (gatingBlockCount() does)");
        test::ok (! e.R.truePeakDbtp.valid, "the true peak refuses instead of publishing -240 dBTP");
        test::ok (! e.R.plrDb.valid, "PLR refuses");
        test::ok (! e.R.lraLu.valid, "LRA refuses instead of publishing 0.0 LU");

        const auto s = run (silence (2, (std::size_t) (kFs * 5.0)), kFs, 1024, { 1024 });
        test::ok (s.R.leadingDigitalSilenceSamples == (std::int64_t) (kFs * 5.0),
                  "an all-silent programme is all leading digital silence");
        test::ok (s.R.programmeSpanSamples == 0, "…with an empty span");
        test::ok (! s.R.stereoCorrelation.valid && s.R.stereoCorrelation.reason == Reason::SilentProgramme,
                  "…and the correlation refuses rather than reading the primitive's +1.0");
        test::ok (! s.R.crestFactorDb.valid, "…and the crest factor refuses rather than dividing by zero");
        test::ok (! s.R.lraLu.valid, "…and LRA refuses, which a constant tone's honest 0.0 LU could not be told from");
    }

    test::group ("a constant tone's LRA is a MEASUREMENT of 0 LU, and is not the silent case");
    {
        const double amp = std::pow (10.0, -23.0 / 20.0) * std::sqrt (2.0);
        Params pm; pm.maxDurationSec = 60.0;
        const auto o = run (tone (2, (std::size_t) (kFs * 30.0), 1000.0, amp), kFs, 4096, { 4096 }, pm);
        test::ok (o.R.lraLu.valid, "a 30 s constant tone HAS an LRA");
        test::approx (o.R.lraLu.value, 0.0, 0.11, "…and it is 0 LU, to one histogram bin");
        test::ok (o.R.shortTermGatedObservations >= 20, "over many gated observations");
    }

    test::group ("non-finite input: excluded where exclusion is exact, fatal where state carries it");
    {
        auto p = programme (2, (std::size_t) (kFs * 4.0), 5u);
        p[0][60000] = std::numeric_limits<float>::quiet_NaN();
        p[1][60001] = std::numeric_limits<float>::infinity();
        const auto o = run (p, kFs, 1024, { 1024 });
        test::ok (o.R.nonFiniteInputSamples == 2, "both holes are counted");
        test::ok (o.R.channel[0].nonFiniteSamples == 1 && o.R.channel[1].nonFiniteSamples == 1,
                  "…per channel");
        test::ok (o.R.rms.valid && o.R.channel[0].dcOffset.valid,
                  "the SUMS stay valid — dropping one term of a sum is exact");
        test::ok (std::isfinite (o.R.rms.value) && std::isfinite (o.R.channel[0].dcOffset.value),
                  "…and finite");
        test::ok (! o.R.infraLowFraction.valid && o.R.infraLowFraction.reason == Reason::NonFiniteInput,
                  "the infra-low fraction refuses: the filter carried the substitution forward");
        test::ok (! o.R.integratedLufs.valid && ! o.R.lraLu.valid && ! o.R.truePeakDbtp.valid
                  && ! o.R.plrDb.valid && ! o.R.shortTermP95.valid,
                  "the whole recursive family refuses with a reason");
        test::ok (o.R.integratedLufs.reason == Reason::NonFiniteInput, "…and names the input as the cause");

        // Every field of the report stays finite: a NaN is never the absence marker.
        bool allFinite = true;
        o.R.visitValues ([&] (const char*, int, const analysis::ProgrammeValue& v)
                         { if (! std::isfinite (v.value)) allFinite = false; });
        test::ok (allFinite, "no field of the report is non-finite");
        for (const Event& e : o.trace)
            if (! (std::isfinite (e.a) && std::isfinite (e.b) && std::isfinite (e.c))) allFinite = false;
        test::ok (allFinite, "…and no element of the trace is either");
    }

    test::group ("nothing finite at all: a MAXIMUM has no value, and 0.0 would be the comfortable zero");
    {
        // A maximum is not exclusion-exact — the sample dropped for being non-finite may have been the
        // largest — so with no finite sample there is no peak. Every other fixture here has at most a
        // handful of holes among hundreds of thousands of good samples, which is why a crew round could
        // make samplePeak unconditionally valid and the suite stayed green.
        const auto nan = std::numeric_limits<float>::quiet_NaN();
        Planes p (2, std::vector<float> (24000, nan));
        const auto o = run (p, kFs, 1024, { 1024 });
        test::ok (o.R.channel[0].finiteSamples == 0 && o.R.channel[0].nonFiniteSamples == 24000,
                  "every sample of both channels is a hole");
        test::ok (! o.R.channel[0].samplePeak.valid
                  && o.R.channel[0].samplePeak.reason == Reason::NoFiniteSamples,
                  "the channel's sample peak is refused, not published as 0.0");
        test::ok (b64 (o.R.channel[0].samplePeak.value) == b64 (0.0), "…with the canonical +0.0 beside it");
        test::ok (! o.R.samplePeak.valid && o.R.samplePeak.reason == Reason::NoFiniteSamples,
                  "and so is the programme's");
        test::ok (! o.R.rms.valid && ! o.R.crestFactorDb.valid, "…as is everything derived from them");
        test::ok (o.R.nonFiniteInputSamples == 48000, "and the damage is counted in full");
    }

    test::group ("a peak under gainToDb's floor is refused, not published as -240 dBTP");
    {
        // `core::gainToDb` floors its argument at 1e-12 (Math.h:29), so a programme whose reference true
        // peak lands at or below that floor would publish exactly -240.0 dBTP — a floored constant reading
        // as a measurement. A crew round relaxed the gate from the floor to 0.0 and nothing noticed.
        Planes p (1, std::vector<float> (24000, 0.0f));
        for (std::size_t i = 0; i < 24000; ++i) p[0][i] = (i & 1u) ? 1.0e-14f : -1.0e-14f;
        const auto o = run (p, kFs, 1024, { 1024 });
        test::ok (! o.R.truePeakDbtp.valid && o.R.truePeakDbtp.reason == Reason::SilentProgramme,
                  "a true peak at 1e-14 is refused rather than floored to -240 dBTP");
        test::ok (! o.R.plrDb.valid, "…and PLR refuses with it");
        test::ok (o.R.channel[0].samplePeak.valid, "…while the SAMPLE peak, which has no floor, is reported");
        test::approx (o.R.channel[0].samplePeak.value, 1.0e-14, 1e-20, "…as the value it actually is");
    }

    test::group ("a hole inside a silence run makes that run's COUNT say cannot-say");
    {
        Planes p (1, std::vector<float> ((std::size_t) (kFs * 2.0), 0.0f));
        for (std::size_t i = 24000; i < 48000; ++i) p[0][i] = 0.3f;
        p[0][1000] = std::numeric_limits<float>::quiet_NaN();     // inside the LEADING silence
        const auto o = run (p, kFs, 512, { 512 });
        test::ok (! o.R.leadingSilenceSamples.valid
                  && o.R.leadingSilenceSamples.reason == Reason::NonFiniteInput,
                  "the leading run carries a hole, so its length is not a measurement");
        test::ok (o.R.trailingSilenceSamples.valid, "the trailing run is clean, and stays valid");
    }

    test::group ("finite input can still overflow a filter, and that is its own reason");
    {
        Planes p (1, std::vector<float> (4096, std::numeric_limits<float>::max()));
        const auto o = run (p, kFs, 1024, { 1024 });
        test::ok (o.R.nonFiniteInputSamples == 0, "every input sample was finite");
        test::ok (o.R.nonFiniteIntermediates > 0, "…and something downstream still overflowed ("
                                                  + std::to_string (o.R.nonFiniteIntermediates) + ")");
        test::ok (! o.R.infraLowFraction.valid
                  && o.R.infraLowFraction.reason == Reason::NonFiniteIntermediate,
                  "the infra-low fraction refuses with NonFiniteIntermediate, not with NonFiniteInput");
        bool allFinite = true;
        o.R.visitValues ([&] (const char*, int, const analysis::ProgrammeValue& v)
                         { if (! std::isfinite (v.value)) allFinite = false; });
        test::ok (allFinite, "and every published field is still finite");
    }

    test::group ("the tail window's COORDINATES, and the invariant behind them");
    {
        // `tailWritten_` is meant to equal `lastSignal_ + 1` at all times — that is what lets the ring hold
        // the last N frames ENDING AT the last signal frame while the trailing silence is held back. Nothing
        // asserted the coordinate, only the value computed from it, and a crew round showed the skip in
        // `tailWritten_ += pending - keep` could be dropped entirely with the suite staying green: both the
        // write and the read go through the same counter, so the last N entries are the same either way and
        // only the published coordinate moves. So assert the coordinate, against the span.
        const std::size_t head = (std::size_t) (kFs * 0.3), gap = (std::size_t) (kFs * 0.3);
        const std::size_t burst = (std::size_t) (kFs * 0.05), tailPad = (std::size_t) (kFs * 0.4);
        Planes p (1, std::vector<float> (head + gap + burst + tailPad, 0.0f));
        for (std::size_t i = 0; i < head; ++i)
            p[0][i] = (float) (0.4 * std::cos (kTwoPi * 180.0 * (double) i / kFs));
        for (std::size_t i = 0; i < burst; ++i)
            p[0][head + gap + i] = (float) (0.25 * std::cos (kTwoPi * 300.0 * (double) i / kFs));
        const auto o = run (p, kFs, 1024, { 1024 });
        const std::int64_t lastSignal = (std::int64_t) (head + gap + burst) - 1;
        const std::int64_t N = 4800;
        int seen = 0;
        for (const Event& e : o.trace)
            if (e.kind == analysis::ProgrammeTraceKind::TailWindow)
            {
                ++seen;
                test::ok (e.end == lastSignal + 1,
                          "the tail window ENDS one past the last signal frame (" + std::to_string (e.end)
                          + " vs " + std::to_string (lastSignal + 1) + ")");
                test::ok (e.begin == lastSignal + 1 - N, "…and begins exactly a window earlier ("
                                                         + std::to_string (e.begin) + ")");
                test::approx (e.b, (double) N, 0.0, "…and counts every frame in it (mono: N samples)");
            }
        test::ok (seen == 1, "finish() publishes the tail window exactly once");
    }

    test::group ("the histogram's BOTTOM edge: an observation under the gate is binned at 0, not at 1");
    {
        // `binOf` is called on EVERY observation for the trace, gated or not, and an observation below the
        // -70 LUFS gate has a LUFS under the histogram's floor, so the clamp at the bottom is live there
        // even though the percentile path never sees such a bin. A crew round moved that clamp from 0 to 1
        // and the suite stayed green, because nothing read the bin of an ungated observation.
        Params pm; pm.maxDurationSec = 60.0;
        const auto o = run (tone (2, (std::size_t) (kFs * 8.0), 1000.0, 1.0e-6), kFs, 4096, { 4096 }, pm);
        int ungated = 0;
        for (const Event& e : o.trace)
            if (e.kind == analysis::ProgrammeTraceKind::ShortTerm && e.c < 0.5)
            {
                ++ungated;
                test::approx (e.b, 0.0, 0.0, "an observation under the absolute gate is binned at 0");
            }
        test::ok (ungated >= 2, "the fixture produced ungated observations to bin (" + std::to_string (ungated) + ")");
        test::ok (! o.R.lraLu.valid, "…and the LRA is refused, there being nothing gated in");
    }

    test::group ("a span shorter than the tail window is refused, however long the FILE is");
    {
        // 0.5 s of leading digital silence, 20 ms of music, 0.1 s of trailing silence. The FILE position of
        // the last signal frame is 24 960, comfortably past the 4 800-frame window, but the PROGRAMME is
        // only 960 frames long — so a window would be 80 % leading silence, frames every other mean in this
        // report excludes. A crew round showed the file-position form of this test surviving the suite.
        const std::size_t lead = (std::size_t) (kFs * 0.5), body = (std::size_t) (kFs * 0.02);
        Planes p (1, std::vector<float> (lead + body + (std::size_t) (kFs * 0.1), 0.0f));
        for (std::size_t i = 0; i < body; ++i)
            p[0][lead + i] = (float) (0.4 * std::cos (kTwoPi * 500.0 * (double) i / kFs));
        const auto o = run (p, kFs, 512, { 512 });
        test::ok (o.R.programmeSpanSamples == (std::int64_t) body, "the span is the music, not the file");
        test::ok (o.R.totalSamples > 4800, "…while the file is far longer than the tail window");
        test::ok (! o.R.tailEnergyRatio.valid
                  && o.R.tailEnergyRatio.reason == Reason::ShorterThanTailWindow,
                  "the tail ratio is refused: the window does not fit inside the PROGRAMME");
    }

    test::group ("the trace's own overflow counter sees the event finish() adds");
    {
        // The tail event is emitted from finish(), AFTER the report's counters would naturally be read. Size
        // the buffer to hold every streaming event and not that one: the overflow must still be reported.
        auto p = programme (1, (std::size_t) (kFs * 1.0), 91u);
        const auto big = run (p, kFs, 1024, { 1024 }, Params {}, -1, -1, 1u << 16);
        const std::size_t total = big.trace.size();
        test::ok (big.traceOverflow == 0 && total > 2, "the big buffer held everything ("
                                                      + std::to_string (total) + " events)");
        const auto edge = run (p, kFs, 1024, { 1024 }, Params {}, -1, -1, total - 1u);
        test::ok (edge.trace.size() == total - 1u, "the tight buffer holds one event fewer");
        test::ok (edge.traceOverflow == 1,
                  "…and the report counts the ONE event it could not hold — the tail event finish() adds, "
                  "which a counter read at the top of buildReport() would have missed (got "
                  + std::to_string (edge.traceOverflow) + ")");
    }

    test::group ("both damage kinds at once: the INPUT is named as the cause, not the overflow it caused");
    {
        // A crew round reordered the reason ladder in buildDynamics() and every check stayed green, because
        // no fixture had both kinds of damage: a NaN input is sanitised before the filters, so it never
        // produces a non-finite intermediate, and FLT_MAX produces one with no bad input. With both
        // present, the order is observable — and the INPUT is the root cause a reader needs first.
        Planes p (1, std::vector<float> (24000, 0.0f));
        for (std::size_t i = 0; i < 24000; ++i)
            p[0][i] = (float) (0.3 * std::sin (kTwoPi * 440.0 * (double) i / kFs));
        for (std::size_t i = 8000; i < 9000; ++i) p[0][i] = std::numeric_limits<float>::max();
        p[0][12000] = std::numeric_limits<float>::quiet_NaN();
        const auto o = run (p, kFs, 1024, { 1024 });
        test::ok (o.R.nonFiniteInputSamples == 1, "one hole");
        test::ok (o.R.nonFiniteIntermediates > 0, "…and a filter overflow from the finite FLT_MAX run");
        test::ok (o.R.integratedLufs.reason == Reason::NonFiniteInput,
                  "with both present the reason names the INPUT, not the overflow");
        test::ok (o.R.infraLowFraction.reason == Reason::NonFiniteInput, "…and so does the infra-low fraction");
    }

    test::group ("a channel that disappears and comes back");
    {
        auto p = programme (2, (std::size_t) (kFs * 3.0), 31u);
        PR a;
        test::run (a.prepare (kFs, 4096, 2));
        const auto v = ptrs (p, 0);
        const float* one[1] { p[0].data() };
        test::run (a.process (v.data(), 2, 48000));
        const float* mid[1] { p[0].data() + 48000 };
        test::run (a.process (mid, 1, 24000));                    // channel 1 goes away
        const auto v2 = ptrs (p, 72000);
        test::run (a.process (v2.data(), 2, (int) (p[0].size() - 72000)));
        a.finish();
        const auto& R = a.report();
        test::ok (R.channel[1].absentSamples == 24000, "the absence is counted, sample for sample");
        test::ok (R.stereoFrames == (std::int64_t) p[0].size() - 24000,
                  "the stereo relation is measured only over the frames that had both channels");
        test::ok (R.channel[0].finiteSamples == (std::int64_t) p[0].size(), "channel 0 lost nothing");
        test::ok (R.integratedLufs.valid, "the loudness family is still a measurement — the gap was silence, not damage");
    }

    // TWO AVERAGES OF ONE PROGRAMME, and the header now says where they part — so the suite says it too.
    // `programmeMeanSquare` here is a SAMPLE total over the programme span, ungated;
    // `BandCrest::programmeMeanSquareDb()` is a mean of BLOCK mean-squares over blocks clearing -70 dBFS.
    // The consumer holds both and asked for the boundary in writing; a sentence in a header that nothing runs
    // is the thing that goes stale, so the two numbers the header cites are measured here.
    test::group ("programmeMeanSquare against BandCrest's: one on a tone, two different things on real shapes");
    {
        struct Case { const char* what; bool quietUnderGate; double wantDelta, tol; };
        // A stationary tone collapses the two: every block identical, every one over the gate, nothing trimmed.
        // A loud smooth envelope leaves only block-granularity against a sample total. States at -62 dBFS put
        // blocks UNDER the gate, which drops them from one population and keeps them in the other.
        const Case cases[] = { { "stationary tone",        false, 0.0,   1.0e-6 },
                               { "smooth loud envelope",   false, -0.039, 0.02  },
                               { "quiet states under -70", true,  3.048, 0.05   } };
        for (const Case& c : cases)
        {
            const int n = (int) (kFs * 20.0);
            Planes q (2, std::vector<float> ((std::size_t) n, 0.0f));
            for (int i = 0; i < n; ++i)
            {
                const double t = (double) i / kFs;
                double s;
                if (c.quietUnderGate)
                {
                    const double env = (((std::size_t) i / (std::size_t) (kFs * 0.7)) % 3 == 0) ? 1.0 : 0.0008;
                    s = env * (0.3 * std::sin (kTwoPi * 220.0 * t) + 0.1 * std::sin (kTwoPi * 1310.0 * t));
                }
                else if (std::string (c.what) == "stationary tone") s = 0.2 * std::sin (kTwoPi * 997.0 * t);
                else s = (0.5 + 0.45 * std::sin (kTwoPi * 0.31 * t)) * 0.3 * std::sin (kTwoPi * 220.0 * t);
                q[0][(std::size_t) i] = q[1][(std::size_t) i] = (float) s;
            }
            analysis::BandCrest bc;
            if (! test::run (bc.prepare (kFs, 2, (long long) n))) continue;
            for (int at = 0; at < n; )
            {
                const int k = std::min (4096, n - at);
                const float* qp[2] { q[0].data() + at, q[1].data() + at };
                if (! test::run (bc.process (qp, 2, k))) break;
                at += k;
            }
            bc.finish();
            Params mp; mp.maxDurationSec = 600.0;
            const auto o = run (q, kFs, 4096, { 4096 }, mp);
            const bool measured = o.R.programmeMeanSquare.valid;
            test::ok (measured, std::string (c.what) + ": PRECONDITION: the report measured this programme");
            if (! measured) continue;
            const double there = 10.0 * std::log10 (o.R.programmeMeanSquare.value);
            const double delta = bc.programmeMeanSquareDb() - there;
            test::approx (delta, c.wantDelta, c.tol,
                          std::string (c.what) + ": crest " + std::to_string (bc.programmeMeanSquareDb())
                          + " dB against report " + std::to_string (there) + " dB, apart by " + std::to_string (delta));
        }
    }

    test::group ("capacity: overflow is DATA — a partial answer is refused, the counters keep counting");
    {
        Params pm; pm.maxDurationSec = 1.0;                       // far shorter than the programme
        const double amp = std::pow (10.0, -23.0 / 20.0) * std::sqrt (2.0);
        const auto o = run (tone (2, (std::size_t) (kFs * 12.0), 1000.0, amp), kFs, 4096, { 4096 }, pm);
        test::ok (o.accepted, "no call was refused because a store filled up");
        test::ok (o.R.loudnessDroppedBlocks > 0, "the dropped gating blocks are counted ("
                                                 + std::to_string (o.R.loudnessDroppedBlocks) + ")");
        test::ok (! o.R.integratedLufs.valid
                  && o.R.integratedLufs.reason == Reason::LoudnessCapacityExceeded,
                  "integrated loudness refuses rather than describing the part that fitted");
        test::ok (! o.R.plrDb.valid, "PLR refuses with it");
        test::ok (o.R.shortTermDroppedObservations > 0 || o.R.lraLu.valid == false,
                  "the short-term series either fitted or says it did not");
        // WHERE IT ACTUALLY TURNS, which that disjunction does not say and nothing else pinned. The store holds
        // floor(hops) + 8 observations against a production of floor(hops) - 29, so a programme is refused once
        // it runs 38 hops past its declared length. THAT MOVED WITH v0.43: at one observation a second the same
        // arithmetic gave about a hundred hops of undeclared slack, so a file 5 s over its declared length used
        // to slip through and is now refused — measured against a 10 s declaration at 48 kHz, the turn was at
        // 20 s and is at 13.8. The contract ("past maxDurationSec it refuses") did not change; the slack behind
        // it shrank to match the words. In hops, not seconds, so the two rates give the same two numbers.
        for (const double rate : { 48000.0, 44100.0 })
        {
            const long long hop = 10LL * std::llround (0.01 * rate);
            const long long declHops = 100;
            for (const long long over : { 37LL, 38LL })
            {
                Params bp; bp.maxDurationSec = (double) (declHops * hop) / rate;
                const auto b2 = run (tone (2, (std::size_t) ((declHops + over) * hop), 1000.0, amp, rate),
                                     rate, 4096, { 4096 }, bp);
                const std::string where = std::to_string ((int) rate) + " Hz, " + std::to_string (over) + " hops over";
                test::ok ((b2.R.shortTermDroppedObservations > 0) == (over == 38),
                          where + ": dropped " + std::to_string (b2.R.shortTermDroppedObservations)
                                + " (37 over still fits, 38 does not)");
                test::ok (b2.R.lraLu.valid == (over == 37),
                          where + ": the LRA is " + (b2.R.lraLu.valid ? "a measurement" : "refused"));
            }
        }
        test::ok (o.R.totalSamples == (std::int64_t) (kFs * 12.0), "and every sample was still measured");
        test::ok (o.R.samplePeak.valid && o.R.rms.valid, "the sample-domain measurements are untouched");
    }

    test::group ("the trace buffer is the caller's, and a small one is counted rather than grown");
    {
        auto p = programme (1, (std::size_t) (kFs * 2.0), 77u);
        const auto big = run (p, kFs, 1024, { 1024 }, Params {}, -1, -1, 1u << 16);
        const auto tiny = run (p, kFs, 1024, { 1024 }, Params {}, -1, -1, 32u);
        test::ok (tiny.trace.size() == 32u, "the small buffer holds exactly its capacity");
        test::ok (tiny.traceOverflow > 0, "the overflow is counted");
        test::ok (big.traceOverflow == 0, "the large one overflows not at all");
        for (std::size_t i = 0; i < tiny.trace.size(); ++i)
        {
            const Event& x = tiny.trace[i]; const Event& y = big.trace[i];
            test::ok (x.kind == y.kind && x.index == y.index && b64 (x.a) == b64 (y.a),
                      "the prefix is the same trace");
            if (i > 4) break;
        }
        std::string why;
        // The report is the same either way, with nothing to zero out by hand first: the overflow count
        // lives on the accessor precisely so the REPORT is not a function of the caller's buffer size.
        test::ok (sameReport (big.R, tiny.R, why), "recording changes no measurement (" + why + ")");
    }

    test::group ("no allocation inside process() or finish()");
    {
        auto p = programme (2, (std::size_t) (kFs * 3.5), 4u);
        std::vector<Event> buf (1u << 16);
        PR a;
        test::run (a.prepare (kFs, 1024, 2));
        a.setTraceBuffer (buf.data(), buf.size());
        const auto v = ptrs (p, 0);
        const long long before = alloc::count.load();
        for (std::size_t at = 0; at < p[0].size(); )
        {
            const int n = (int) std::min<std::size_t> (997, p[0].size() - at);
            const auto vv = ptrs (p, at);
            test::run (a.process (vv.data(), 2, n));
            at += (std::size_t) n;
        }
        a.finish();
        const long long calls = alloc::count.load() - before;
        // `ptrs` itself allocates a vector per call, so the count is taken against a run that does the same
        // work without the analyzer: what is asserted is that the analyzer adds nothing.
        long baseline = 0;
        {
            const long long b0 = alloc::count.load();
            for (std::size_t at = 0; at < p[0].size(); )
            {
                const int n = (int) std::min<std::size_t> (997, p[0].size() - at);
                const auto vv = ptrs (p, at);
                (void) vv;
                at += (std::size_t) n;
            }
            baseline = alloc::count.load() - b0;
        }
        test::okNoAlloc (calls == baseline, "process() + finish() allocate nothing of their own (got "
                                            + std::to_string (calls) + " against " + std::to_string (baseline) + ")");
        (void) v;
    }

    test::group ("rates: the same samples, every supported rate, and the sub-hop follows the rate");
    {
        // 11 050 Hz IS THE POINT OF THIS LIST. Every other rate here has an integral 0.01*fs, so a crew
        // round turned the sub-hop's `lround` into `floor` and the whole suite stayed green. At 11 050 Hz
        // the sub-hop is lround(110.5) = 111 samples and floor gives 110 — a different measurement grid.
        for (double fs : { 8000.0, 11050.0, 22050.0, 44100.0, 48000.0, 88200.0, 96000.0, 192000.0 })
        {
            Params pm;
            const auto st = PR::storageFor (fs, 1024, 2, pm);
            test::ok (st.ok, "prepared at " + std::to_string ((int) fs) + " Hz");
            test::ok (st.subHopSamples == (std::int64_t) std::lround (0.01 * fs),
                      "the sub-hop is lround(0.01 fs) at " + std::to_string ((int) fs) + " Hz");
            if (fs > 11049.0 && fs < 11051.0)
                test::ok (st.subHopSamples == 111, "…and at 11 050 Hz that is 111 samples, not floor's 110");
            if (fs > 22049.0 && fs < 22051.0)
                test::ok (st.subHopSamples == 221, "…and at 22 050 Hz 221, not 220");
            const auto o = run (tone (2, (std::size_t) (fs * 1.5), 100.0, 0.4, fs), fs, 1024, { 333 }, pm);
            test::approx (o.R.crestFactorDb.value, 3.0103, 0.05,
                          "a sine still reads 3.01 dB of crest at " + std::to_string ((int) fs) + " Hz");
        }
    }

    test::group ("negative zero is digital silence, and it is value-equal but not bit-equal to +0.0");
    {
        // -0.0f has exactly zero amplitude, so treating it as digital silence is right, and IEEE equality
        // agrees (`exactlyEqual(-0.0f, 0.0f)` is true). What must NOT agree is the bit comparison: a file
        // whose left channel carries -0.0 and right +0.0 is value-identical and bit-different, and the
        // report publishes both counts precisely so that distinction survives.
        const std::size_t n = 24000;
        Planes p (2, std::vector<float> (n, 0.0f));
        for (std::size_t i = 0; i < n; ++i) { p[0][i] = -0.0f; p[1][i] = 0.0f; }
        const auto o = run (p, kFs, 1000, { 1000 });
        test::ok (o.R.programmeSpanSamples == 0, "an all -0.0 / +0.0 programme has an EMPTY span");
        test::ok (o.R.leadingDigitalSilenceSamples == (std::int64_t) n, "…it is digital silence end to end");
        test::ok (o.R.stereoValueEqualFrames == (std::int64_t) n, "every frame is VALUE-equal");
        test::ok (o.R.stereoBitIdenticalFrames == 0, "…and not one is BIT-equal");
        test::ok (! o.R.rms.valid && o.R.rms.reason == Reason::SilentProgramme, "and there is nothing to measure");
    }

    test::group ("one channel silent, the other loud");
    {
        // Missing until a crew round asked for it. The span is active (the left channel is not silent), so
        // the right channel's own fields must each say the right thing rather than inheriting the frame's.
        const std::size_t n = (std::size_t) (kFs * 1.0);
        Planes p (2, std::vector<float> (n, 0.0f));
        for (std::size_t i = 0; i < n; ++i)
            p[0][i] = (float) (0.4 * std::cos (kTwoPi * 220.0 * (double) i / kFs));
        const auto o = run (p, kFs, 997, { 997 });
        test::ok (o.R.programmeSpanSamples > 0, "the span is active — one channel carries signal");
        test::ok (o.R.channel[1].rms.valid && b64 (o.R.channel[1].rms.value) == b64 (0.0),
                  "the silent channel's RMS is a valid ZERO, not a refusal");
        test::ok (o.R.channel[1].dcOffset.valid && b64 (o.R.channel[1].dcOffset.value) == b64 (0.0),
                  "…and so is its DC");
        test::ok (! o.R.channel[1].crestFactorDb.valid
                  && o.R.channel[1].crestFactorDb.reason == Reason::SilentProgramme,
                  "…but its crest factor is refused rather than dividing by zero");
        test::ok (! o.R.channel[1].infraLowFraction.valid
                  && o.R.channel[1].infraLowFraction.reason == Reason::SilentProgramme,
                  "…and so is its energy fraction, which has no denominator");
        test::ok (o.R.channel[0].crestFactorDb.valid, "the loud channel is unaffected");
        test::ok (! o.R.stereoCorrelation.valid && o.R.stereoCorrelation.reason == Reason::SilentProgramme,
                  "the correlation is refused: one side contributes no energy at all");
        test::ok (! o.R.stereoBalanceDb.valid, "…and the balance would be an infinite ratio");
        // Side and mid are each a quarter of the live channel's energy, so their ratio is exactly 1.
        test::ok (o.R.stereoSideToMidRatio.valid, "side/mid IS defined here");
        test::approx (o.R.stereoSideToMidRatio.value, 1.0, 1e-12,
                      "…and a signal in one channel only is exactly half mid, half side");
    }

    test::group ("the two silence thresholds: the parameter's, and the depth-free one");
    {
        // A 16-bit dither tail: +-1 LSB = -90.31 dBFS. At the -96 dBFS default it is SIGNAL; a -90 dBFS
        // threshold — 0.31 dB above that LSB — would have called the whole tail silence.
        const float lsb = 1.0f / 32768.0f;
        Planes p (1, std::vector<float> ((std::size_t) (kFs * 1.0), 0.0f));
        for (std::size_t i = 24000; i < 48000; ++i) p[0][i] = (i & 1u) ? lsb : -lsb;
        const auto o = run (p, kFs, 512, { 512 });
        test::ok (o.R.trailingSilenceSamples.value == 0,
                  "at -96 dBFS a 16-bit dither tail is signal, not silence");
        Params trap; trap.silenceThresholdDb = -90.0;
        const auto t = run (p, kFs, 512, { 512 }, trap);
        test::ok (t.R.trailingSilenceSamples.value == 48000 && t.R.leadingSilenceSamples.value == 48000,
                  "at -90 dBFS not one sample of the file clears the threshold, so the WHOLE programme reads "
                  "silent — 1 LSB of 16-bit is -90.309 dBFS, 0.309 dB UNDER a -90 dBFS gate, and that is the "
                  "trap the -96 dBFS default exists to avoid");
        test::ok (o.R.trailingDigitalSilenceSamples == 0 && t.R.trailingDigitalSilenceSamples == 0,
                  "the digital-silence count has no threshold and does not move");
    }
}

int main()
{
    std::printf ("felitronics::analysis::ProgrammeReport tests\n");
    testStorageAndRefusals();
    testClosedForms();
    testNegatives();
    testInvariance();
    testAgainstReference();
    testEdges();
    return test::report();
}

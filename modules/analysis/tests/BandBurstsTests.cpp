// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// BandBursts self-tests.
//
// TWO KINDS OF TRUTH, BECAUSE ONE OF THEM PROVES NOTHING ON ITS OWN.
//   · AN INDEPENDENT ORACLE. The events, the hop energies and the baselines are recomputed by a
//     WHOLE-PROGRAMME reference written straight from the rule — hop energies into a plain array, the
//     baseline by `std::sort` over a copied slice (not `nth_element` over a ring), the event machine as a
//     flat scan — and the engine must null it field for field. This is the half that can catch a
//     deterministically WRONG schedule.
//   · RE-SLICING (law 8a). The same programme cut 14 ways plus five random splittings and four different
//     `maxBlock` values, compared BITWISE over every hop of the trace and every field of the report.
//     A re-slicing test alone certifies self-consistency and not correctness: one program under two
//     slicings is still one program, so a wrong-but-deterministic rule passes every slicing green. That
//     is why the oracle above exists, and why the mutant pass below includes a mutant that re-slicing
//     provably cannot see.
//
// MUTANT PASS — 18 of 19 mutants are caught (run by hand against this suite; each mutation is applied to
// the header, the object file is DELETED first, and a run whose build shows no `Building CXX` line is
// refused rather than scored — make's mtime granularity is one second and a same-second header edit is
// invisible to it, which is worth three false greens if you let it be):
//   m01 hop closed at the call boundary as well as on the absolute count   -> RED (40 checks)
//   m02 denormal flush moved to the end of process()                       -> RED  (poison+silence fixture)
//   m03 baseline made INCLUSIVE (ring written before the decision)         -> RED  (the oracle)
//   m04 eligibility k >= N  ->  k > N                                      -> RED
//   m05 eligibility k >= N  ->  k >= N-1                                   -> RED
//   m06 zero-baseline gate removed                                         -> RED  (silent lead-in: three
//       counter checks. NOT via a non-finite dB — `closeEvent`'s peakDen_ guard publishes 0.0 — so the
//       "no event carries a non-finite excess" check passes even on the mutant, and the counters are
//       what catch it.
//   m07 peak seed 0/0 instead of the opening hop                           -> RED
//   m08 peak tie-break > -> >=                                             -> RED  (the bitwise tie fixture)
//   m09 partial tail hop admitted to the ring                              -> RED
//   m10 partial tail hop allowed to open an event                          -> RED
//   m11 hysteresis removed (exit := enter)                                 -> RED
//   m12 onset counters stopped at maxEvents                                -> RED
//   m13 absent channel skipped instead of driven with 0.0f                 -> RED
//   m14 10*log10 -> 20*log10 on a power ratio                              -> RED
//   m15 armedness read from the flag alone (a moved-from detector stays armed) -> RED, by CRASHING: the
//       mutant indexes the empty channel list, which is the undefined behaviour the fix exists to stop.
//       A crash is a weaker witness than a failed check, and it is recorded as what it is.
//   m16 partial-tail damage not flagged on the event finish() closes       -> RED
//   m17 corners validated as doubles, not as the floats the filter is given -> RED (8 checks)
//   m18 the baselineMs upper bound removed                                 -> *** GREEN, NOT CAUGHT ***
//       and it is not catchable behaviourally. The bound exists so that `llround` is never handed an
//       argument outside int64 (baselineMs = 1e30 does); but whatever that undefined conversion returns,
//       the `nb < 1 || nb > kMaxBaselineHops` check downstream refuses it anyway, so the observable
//       answer is the same refusal either way. The only difference is whether the UB executes. It is
//       kept as a defensive bound and scored as uncaught rather than dropped to make the tally look
//       round.
//   m19 threshold finiteness/cap removed (enterDb = +inf silently disables detection) -> RED
// A gate that never reddened is not a gate; m02 and m13 are the two the brief mandated, and m02 needs the
// POISON fixture (a finite sign-alternating 3e38 that overflows the SVF) plus digital silence — on healthy
// material the flush cannot change a bit by construction (`Svf.h:172`).

#include <felitronics_test.h>
#include <felitronics/analysis/BandBursts.h>

#include <algorithm>
#include <atomic>
#include <bit>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <random>
#include <string>
#include <vector>

static std::atomic<long> g_allocs { 0 };
void* operator new      (std::size_t s) { g_allocs.fetch_add (1); return std::malloc (s ? s : 1); }
void* operator new[]    (std::size_t s) { g_allocs.fetch_add (1); return std::malloc (s ? s : 1); }
void  operator delete   (void* p) noexcept { std::free (p); }
void  operator delete[] (void* p) noexcept { std::free (p); }
void  operator delete   (void* p, std::size_t) noexcept { std::free (p); }
void  operator delete[] (void* p, std::size_t) noexcept { std::free (p); }

using namespace felitronics;
using BB    = analysis::BandBursts;
using Params = analysis::BandBurstsParams;
using Burst  = analysis::BandBurst;
using Trace  = analysis::BandBurstsHopTrace;
using Planes = std::vector<std::vector<float>>;

static constexpr double kTau = 6.283185307179586;
static std::uint64_t bits (double d) noexcept { return std::bit_cast<std::uint64_t> (d); }

// A small, fast configuration for the heavy invariance work: 10 ms hops, a 20-hop (200 ms) ring.
static Params fastParams()
{
    Params p;
    p.hopMs = 10.0; p.baselineMs = 200.0;
    p.enterDb = 6.0; p.exitDb = 3.0;
    p.maxEvents = 1 << 12;
    return p;
}

// ============================================================================================ material
static std::vector<const float*> ptrs (const Planes& p, std::size_t at)
{
    std::vector<const float*> v;
    for (const auto& ch : p) v.push_back (ch.data() + at);
    return v;
}

// Nonstationary material with bursts at KNOWN hop positions: a broadband noise bed plus a 7 kHz burst
// (windowed, so it has a real attack) at each requested hop. A flat tone would be a weak witness.
static Planes material (std::size_t frames, int nch, const std::vector<int>& burstHops, int hopSamples,
                        int burstHopLen, double bed, double burst, unsigned seed)
{
    std::mt19937 rng (seed);
    std::uniform_real_distribution<float> u (-1.0f, 1.0f);
    Planes p ((std::size_t) nch, std::vector<float> (frames, 0.0f));
    for (std::size_t i = 0; i < frames; ++i)
        for (int c = 0; c < nch; ++c)
            p[(std::size_t) c][i] = (float) (bed * (double) u (rng));

    for (int h : burstHops)
    {
        const std::size_t a = (std::size_t) h * (std::size_t) hopSamples;
        const std::size_t b = std::min (frames, a + (std::size_t) burstHopLen * (std::size_t) hopSamples);
        const double len = (double) (b - a);
        for (std::size_t i = a; i < b; ++i)
        {
            const double w = 0.5 - 0.5 * std::cos (kTau * (double) (i - a) / len);    // no click at the seam
            const double v = burst * w * std::sin (kTau * 7000.0 * (double) i / 48000.0);
            for (int c = 0; c < nch; ++c) p[(std::size_t) c][i] = (float) ((double) p[(std::size_t) c][i] + v);
        }
    }
    return p;
}

// ====================================================================================== the run harness
struct Report
{
    bool ok = false, finished = false, eventsValid = false, energyValid = false, eventsComplete = false;
    std::int64_t samples = 0, hops = 0, eligible = 0, zeroBase = 0, burstHops = 0, damaged = 0;
    std::int64_t tailSamples = 0, eventCount = 0, stored = 0, onsets = 0, intervals = 0, ioiOverflow = 0;
    std::int64_t overflowSamples = 0, firstNonFinite = 0;
    double tailEnergy = 0.0;
    int modal = 0; std::int64_t modalMass = 0;
    int hopSamples = 0, baselineHops = 0;
    std::vector<double> chEnergy;
    std::vector<std::int64_t> chNonFinite, chAbsent;
    std::vector<Burst> events;
    std::vector<std::int64_t> ioi, lag;
};

static Report snapshot (const BB& d)
{
    Report r;
    r.ok = true; r.finished = d.isFinished();
    r.eventsValid = d.eventsValid(); r.energyValid = d.programmeEnergyValid();
    r.eventsComplete = d.eventsComplete();
    r.samples = d.samplesProcessed(); r.hops = d.hopCount(); r.eligible = d.eligibleHops();
    r.zeroBase = d.zeroBaselineHops(); r.burstHops = d.burstHops(); r.damaged = d.damagedHops();
    r.tailSamples = d.tailPartialSamples(); r.tailEnergy = d.tailPartialEnergy();
    r.eventCount = d.eventCount(); r.stored = d.storedEventCount();
    r.onsets = d.onsetCount(); r.intervals = d.intervalCount(); r.ioiOverflow = d.intervalOverflow();
    r.overflowSamples = d.overflowSamples(); r.firstNonFinite = d.firstNonFiniteAt();
    r.modal = d.modalIntervalHops(); r.modalMass = d.modalIntervalMass();
    r.hopSamples = d.hopSamples(); r.baselineHops = d.baselineHops();
    for (int c = 0; c < d.channels(); ++c)
    {
        r.chEnergy.push_back (d.bandEnergy (c));
        r.chNonFinite.push_back (d.nonFiniteSamples (c));
        r.chAbsent.push_back (d.absentSamples (c));
    }
    for (std::int64_t i = 0; i < r.stored; ++i) r.events.push_back (d.event (i));
    for (int b = 1; b <= BB::kIoiBins; ++b) r.ioi.push_back (d.intervalBin (b));
    for (int b = 1; b <= BB::kMaxLag;  ++b) r.lag.push_back (d.lagBin (b));
    return r;
}

struct Cap { std::vector<Trace> v; };
static void observe (void* u, const Trace& t) { static_cast<Cap*> (u)->v.push_back (t); }

struct Run { Report rep; std::vector<Trace> trace; };

// Runs the detector over `p`, cutting the stream with the sizes in `split` (cycled). `nchCall` < 0 means
// "pass every prepared channel"; otherwise that many, so the absent-channel path is exercised.
static Run runIt (const Planes& p, const Params& par, double fs, int maxBlock, const std::vector<int>& split,
                  bool wantTrace = true, int nchCall = -1)
{
    Run out;
    Cap cap; if (wantTrace) cap.v.reserve (1u << 16);
    BB d;
    d.setParams (par);
    if (wantTrace) d.setObserver (&observe, &cap);
    const int nch = (int) p.size();
    if (! d.prepare (fs, maxBlock, nch)) return out;
    const std::size_t frames = p[0].size();
    std::size_t at = 0; std::size_t k = 0;
    while (at < frames)
    {
        const std::size_t n = std::min (frames - at, (std::size_t) std::max (1, split[k++ % split.size()]));
        const auto v = ptrs (p, at);
        if (! d.process (v.data(), nchCall < 0 ? nch : nchCall, (int) n)) return out;
        at += n;
    }
    d.finish();
    out.rep = snapshot (d);
    out.trace = cap.v;
    return out;
}

// ============================================================================== bitwise comparison
static bool sameBurst (const Burst& a, const Burst& b)
{
    return a.start == b.start && a.length == b.length && a.peakAt == b.peakAt && a.hops == b.hops
        && bits (a.peakPower) == bits (b.peakPower) && bits (a.peakBaseline) == bits (b.peakBaseline)
        && bits (a.peakExcessDb) == bits (b.peakExcessDb) && bits (a.peakWidePower) == bits (b.peakWidePower)
        && bits (a.energy) == bits (b.energy)
        && a.touchedNonFinite == b.touchedNonFinite
        && a.baselineTouchedNonFinite == b.baselineTouchedNonFinite
        && a.closedByFinish == b.closedByFinish;
}

// Field by field and BITWISE — not memcmp (padding) and not != (it cannot tell −0.0 from +0.0).
static bool sameTrace (const Trace& a, const Trace& b)
{
    return a.hopIndex == b.hopIndex && a.startSample == b.startSample && a.endSample == b.endSample
        && bits (a.energy) == bits (b.energy) && bits (a.wideEnergy) == bits (b.wideEnergy)
        && bits (a.baseline) == bits (b.baseline)
        && a.full == b.full && a.eligible == b.eligible && a.damaged == b.damaged
        && a.above == b.above && a.inEvent == b.inEvent && a.eventStart == b.eventStart;
}

static bool sameReport (const Report& a, const Report& b)
{
    if (a.ok != b.ok || a.finished != b.finished || a.eventsValid != b.eventsValid
        || a.energyValid != b.energyValid || a.eventsComplete != b.eventsComplete) return false;
    if (a.samples != b.samples || a.hops != b.hops || a.eligible != b.eligible || a.zeroBase != b.zeroBase
        || a.burstHops != b.burstHops || a.damaged != b.damaged || a.tailSamples != b.tailSamples
        || a.eventCount != b.eventCount || a.stored != b.stored || a.onsets != b.onsets
        || a.intervals != b.intervals || a.ioiOverflow != b.ioiOverflow
        || a.overflowSamples != b.overflowSamples || a.firstNonFinite != b.firstNonFinite
        || a.modal != b.modal || a.modalMass != b.modalMass
        || a.hopSamples != b.hopSamples || a.baselineHops != b.baselineHops) return false;
    if (bits (a.tailEnergy) != bits (b.tailEnergy)) return false;
    if (a.chEnergy.size() != b.chEnergy.size() || a.events.size() != b.events.size()) return false;
    for (std::size_t i = 0; i < a.chEnergy.size(); ++i)
        if (bits (a.chEnergy[i]) != bits (b.chEnergy[i]) || a.chNonFinite[i] != b.chNonFinite[i]
            || a.chAbsent[i] != b.chAbsent[i]) return false;
    for (std::size_t i = 0; i < a.events.size(); ++i) if (! sameBurst (a.events[i], b.events[i])) return false;
    for (std::size_t i = 0; i < a.ioi.size(); ++i) if (a.ioi[i] != b.ioi[i]) return false;
    for (std::size_t i = 0; i < a.lag.size(); ++i) if (a.lag[i] != b.lag[i]) return false;
    return true;
}

// ==================================================================================== THE ORACLE
// The same rule, written as a whole-programme program: no ring, no streaming state machine, energies in a
// plain array, the baseline by sorting a copied slice, the events by a flat scan. The filter itself is the
// primitive under composition, so the oracle drives the SAME `eq::Crossover2` pair — with its own explicit
// 64-sample flush counter — and what is independently constructed is everything above it.
struct OracleEvent { std::int64_t start, length, peakAt, hops; double peakPower, peakBaseline, excessDb, energy; };
struct Oracle
{
    std::vector<double> hopE;                 // full hops only
    std::vector<double> baseline;             // 0 where not eligible
    std::vector<char>   eligible;
    std::vector<OracleEvent> events;
    std::int64_t tailSamples = 0; double tailEnergy = 0.0;
};

static void finishOracle (OracleEvent& cur, double pn, double pd, int H, Oracle& o)
{
    cur.peakPower    = pn / (double) H;
    cur.peakBaseline = pd / (double) H;
    cur.excessDb     = pd > 0.0 ? 10.0 * std::log10 (pn / pd) : 0.0;
    o.events.push_back (cur);
}

static Oracle oracle (const Planes& p, const Params& par, double fs)
{
    Oracle o;
    const int nch = (int) p.size();
    const std::size_t frames = p[0].size();
    const int H = (int) std::llround (fs * par.hopMs / 1000.0);
    const int N = (int) std::llround (par.baselineMs / par.hopMs);
    const double enterR = std::pow (10.0, par.enterDb / 10.0);
    const double exitR  = std::pow (10.0, par.exitDb  / 10.0);

    eq::Crossover2 lo, hi;
    lo.prepare (fs, nch); hi.prepare (fs, nch);
    lo.setFrequency ((float) par.bandLowHz); hi.setFrequency ((float) par.bandHighHz);

    std::vector<double> all;                                     // every hop, the last one possibly partial
    double acc = 0.0; int inHop = 0; int gridPhase = 0;
    for (std::size_t i = 0; i < frames; ++i)
    {
        for (int c = 0; c < nch; ++c)
        {
            const float raw = p[(std::size_t) c][i];
            const float x = std::isfinite (raw) ? raw : 0.0f;
            float a = 0.0f, b = 0.0f, band = 0.0f, d2 = 0.0f;
            lo.processSample (c, x, a, b);
            hi.processSample (c, b, band, d2);
            if (std::isfinite (band)) { volatile double q; q = (double) band * (double) band; acc += q; }
        }
        if (++gridPhase == core::StateGrid::kPeriod) { gridPhase = 0; lo.flushDenormals(); hi.flushDenormals(); }
        if (++inHop == H) { all.push_back (acc); acc = 0.0; inHop = 0; }
    }
    o.hopE = all;
    if (inHop > 0) { o.tailSamples = inHop; o.tailEnergy = acc; }

    const std::size_t nHops = o.hopE.size();
    o.baseline.assign (nHops, 0.0); o.eligible.assign (nHops, 0);
    bool in = false; OracleEvent cur {}; double pn = 0.0, pd = 0.0;
    for (std::size_t k = 0; k < nHops; ++k)
    {
        double base = 0.0; bool elig = false;
        if (k >= (std::size_t) N)
        {
            std::vector<double> w (o.hopE.begin() + (std::ptrdiff_t) (k - (std::size_t) N),
                                   o.hopE.begin() + (std::ptrdiff_t) k);
            std::sort (w.begin(), w.end());                      // a SORT, not nth_element
            base = w[(std::size_t) (N / 2)];
            if (base > 0.0) { elig = true; }
        }
        o.baseline[k] = elig ? base : 0.0; o.eligible[k] = elig ? (char) 1 : (char) 0;
        const std::int64_t startS = (std::int64_t) k * (std::int64_t) H;
        if (! elig) { if (in) { cur.length = startS - cur.start; finishOracle (cur, pn, pd, H, o); in = false; } continue; }
        const bool above = o.hopE[k] > base * (in ? exitR : enterR);
        if (! in && above)
        {
            in = true; cur = OracleEvent {}; cur.start = startS; cur.hops = 0; cur.energy = 0.0;
            pn = o.hopE[k]; pd = base; cur.peakAt = startS;
        }
        if (in)
        {
            if (above)
            {
                ++cur.hops; cur.energy += o.hopE[k];
                volatile double l2, r2; l2 = o.hopE[k] * pd; r2 = pn * base;
                if (l2 > r2) { pn = o.hopE[k]; pd = base; cur.peakAt = startS; }
            }
            else { cur.length = startS - cur.start; finishOracle (cur, pn, pd, H, o); in = false; }
        }
    }
    if (in) { cur.length = (std::int64_t) frames - cur.start; finishOracle (cur, pn, pd, H, o); }
    return o;
}

// ======================================================================================= the groups
static void groupOracle()
{
    test::group ("the oracle — a whole-programme reference, independently constructed");
    const Params par = fastParams();
    const int H = 480, N = 20;
    // Nonstationary bed + bursts at hops the fixture CHOSE, well past the run-up.
    const std::vector<int> hops { 25, 40, 41, 42, 70, 95, 96, 130 };
    const Planes p = material (200u * (std::size_t) H, 2, hops, H, 2, 0.02, 0.35, 12345u);

    const Run r = runIt (p, par, 48000.0, 1024, { 1024 });
    const Oracle o = oracle (p, par, 48000.0);

    test::ok (r.rep.ok, "oracle fixture: the engine ran");
    test::ok ((std::size_t) r.rep.hops == o.hopE.size(), "hop COUNT matches the oracle");

    // Every full hop's RAW energy and RAW baseline, bitwise, before any log or threshold.
    std::size_t full = 0, eBad = 0, bBad = 0, elBad = 0;
    for (const auto& t : r.trace)
    {
        if (! t.full) continue;
        const std::size_t k = (std::size_t) t.hopIndex;
        if (k >= o.hopE.size()) break;
        if (bits (t.energy)   != bits (o.hopE[k]))     ++eBad;
        if (bits (t.baseline) != bits (o.baseline[k])) ++bBad;
        if ((t.eligible ? 1 : 0) != (int) o.eligible[k]) ++elBad;
        ++full;
    }
    test::ok (full == o.hopE.size(), "every full hop was traced");
    test::ok (eBad == 0, "hop energy is BITWISE the oracle's (" + std::to_string (eBad) + " differ)");
    test::ok (bBad == 0, "hop baseline is BITWISE the oracle's (" + std::to_string (bBad) + " differ)");
    test::ok (elBad == 0, "eligibility matches the oracle");

    test::ok (r.rep.eventCount == (std::int64_t) o.events.size(),
              "event COUNT matches the oracle (" + std::to_string (r.rep.eventCount) + " vs "
              + std::to_string (o.events.size()) + ")");
    test::ok (r.rep.eventCount > 0, "the fixture produced events at all");
    std::size_t evBad = 0;
    for (std::size_t i = 0; i < o.events.size() && i < r.rep.events.size(); ++i)
    {
        const Burst& a = r.rep.events[i]; const OracleEvent& b = o.events[i];
        if (a.start != b.start || a.length != b.length || a.peakAt != b.peakAt || a.hops != b.hops
            || bits (a.peakPower) != bits (b.peakPower) || bits (a.peakBaseline) != bits (b.peakBaseline)
            || bits (a.peakExcessDb) != bits (b.excessDb) || bits (a.energy) != bits (b.energy)) ++evBad;
    }
    test::ok (evBad == 0, "every event field is BITWISE the oracle's (" + std::to_string (evBad) + " differ)");
    test::ok (r.rep.tailSamples == o.tailSamples, "the uncovered tail matches the oracle");
    test::ok (bits (r.rep.tailEnergy) == bits (o.tailEnergy), "the tail's energy matches the oracle");

    // Every published double must be FINITE. A bound test alone is NaN-blind.
    bool fin = std::isfinite (r.rep.tailEnergy);
    for (const auto& e : r.rep.events)
        fin = fin && std::isfinite (e.peakPower) && std::isfinite (e.peakBaseline)
                  && std::isfinite (e.peakExcessDb) && std::isfinite (e.energy)
                  && std::isfinite (e.peakWidePower);
    for (double v : r.rep.chEnergy) fin = fin && std::isfinite (v);
    test::ok (fin, "every published double is finite");
}

static void groupSlicing()
{
    test::group ("law 8a — bit-identical under arbitrary re-slicing");
    const Params par = fastParams();
    const int H = 480, N = 20, W = H * N;
    const std::vector<int> hops { 25, 26, 55, 80, 81, 82, 110 };
    // 197 hops and a deliberate partial tail; plus digital silence and a POISON burst, so the flush and
    // the two non-finite gates are all inside the material the slicings have to agree on.
    Planes p = material (197u * (std::size_t) H + 133u, 2, hops, H, 2, 0.02, 0.35, 777u);
    for (std::size_t i = 60u * (std::size_t) H; i < 66u * (std::size_t) H; ++i)
        for (auto& ch : p) ch[i] = 0.0f;                                  // digital silence: the flush bites
    p[0][150u * (std::size_t) H + 7u] = std::numeric_limits<float>::quiet_NaN();
    p[1][150u * (std::size_t) H + 9u] = std::numeric_limits<float>::infinity();
    for (int j = 0; j < 8; ++j)                                           // FINITE, and it overflows the SVF
        p[0][170u * (std::size_t) H + (std::size_t) j] = (j % 2 == 0) ? 3.0e38f : -3.0e38f;

    const Run ref = runIt (p, par, 48000.0, 4096, { (int) p[0].size() });
    test::ok (ref.rep.ok, "reference run accepted");
    test::ok (ref.rep.damaged > 0, "the fixture really produced damaged hops");
    test::ok (ref.rep.overflowSamples > 0, "the fixture really overflowed the filter from a FINITE input");
    test::ok (ref.rep.tailSamples == 133, "the fixture really left a partial tail");

    // Boundaries put exactly on the transitions: a hop edge, the end of the ring, a StateGrid tick.
    const std::vector<std::vector<int>> splits {
        { (int) p[0].size() }, { 1 }, { 2 }, { 3 }, { 7 }, { 64 }, { 63 }, { 65 },
        { H - 1 }, { H }, { H + 1 }, { W - 1 }, { W }, { W + 1 }, { 4096 }, { 9000 },
        { 480, 1, 479 }, { 1, 63, 1, 4096 }, { H, 64, 1 }
    };
    std::size_t bad = 0;
    for (const auto& s : splits)
    {
        const Run r = runIt (p, par, 48000.0, 4096, s);
        if (! r.rep.ok || ! sameReport (ref.rep, r.rep) || r.trace.size() != ref.trace.size()) { ++bad; continue; }
        for (std::size_t i = 0; i < r.trace.size(); ++i)
            if (! sameTrace (ref.trace[i], r.trace[i])) { ++bad; break; }
    }
    test::ok (bad == 0, "19 fixed splittings agree bitwise, report AND trace (" + std::to_string (bad) + " differ)");

    std::mt19937 rng (4242u);
    std::uniform_int_distribution<int> ud (1, 3000);
    std::size_t rbad = 0;
    for (int trial = 0; trial < 5; ++trial)
    {
        std::vector<int> s; for (int j = 0; j < 40; ++j) s.push_back (ud (rng));
        const Run r = runIt (p, par, 48000.0, 4096, s);
        if (! r.rep.ok || ! sameReport (ref.rep, r.rep) || r.trace.size() != ref.trace.size()) { ++rbad; continue; }
        for (std::size_t i = 0; i < r.trace.size(); ++i)
            if (! sameTrace (ref.trace[i], r.trace[i])) { ++rbad; break; }
    }
    test::ok (rbad == 0, "5 ragged random splittings agree bitwise (" + std::to_string (rbad) + " differ)");

    // maxBlock sizes NOTHING: a different one must not move a bit, including one far under the calls used.
    std::size_t mbad = 0;
    for (int mb : { 1, 16, 512, 65536 })
    {
        const Run r = runIt (p, par, 48000.0, mb, { 1777 });
        if (! r.rep.ok || ! sameReport (ref.rep, r.rep) || r.trace.size() != ref.trace.size()) { ++mbad; continue; }
        for (std::size_t i = 0; i < r.trace.size(); ++i)
            if (! sameTrace (ref.trace[i], r.trace[i])) { ++mbad; break; }
    }
    test::ok (mbad == 0, "4 different maxBlock values agree bitwise (" + std::to_string (mbad) + " differ)");
}

static void groupNegatives()
{
    test::group ("the negative tests — level is not a burst");
    const Params par = fastParams();
    const int H = 480;

    // (1) THE MAIN ONE. Stationary and BRIGHT: a lot of band energy, evenly. No events.
    {
        std::mt19937 rng (99u);
        std::uniform_real_distribution<float> u (-0.7f, 0.7f);
        Planes p (2, std::vector<float> (600u * (std::size_t) H, 0.0f));
        for (auto& ch : p) for (auto& v : ch) v = u (rng);
        const Run r = runIt (p, par, 48000.0, 2048, { 2048 }, false);
        test::ok (r.rep.ok, "stationary-bright: accepted");
        test::ok (r.rep.eligible > 500, "stationary-bright: hops really were judged");
        test::ok (r.rep.eventCount == 0,
                  "STATIONARY BRIGHT gives no events (" + std::to_string (r.rep.eventCount) + ")");
    }
    // A steady 7 kHz tone in the middle of the band: the same answer, by a different route.
    {
        Planes p (1, std::vector<float> (400u * (std::size_t) H, 0.0f));
        for (std::size_t i = 0; i < p[0].size(); ++i)
            p[0][i] = (float) (0.5 * std::sin (kTau * 7000.0 * (double) i / 48000.0));
        const Run r = runIt (p, par, 48000.0, 2048, { 2048 }, false);
        test::ok (r.rep.eventCount == 0, "a steady in-band tone gives no events");
    }

    // (2) THE ONE THAT IS NOT IN THE BRIEF. The same programme x2 must not move the report. A power of
    // two is exact, so every energy scales by exactly 4 and every RATIO is bit-identical. The material
    // keeps a noise floor throughout: at exact digital silence the 1e-15 state flush is the one place
    // where x2 is NOT homogeneous (measured: it breaks on 27 of 256 grid alignments), and this test is
    // about the instrument, not about that.
    {
        const std::vector<int> hops { 25, 50, 51, 90, 120 };
        const Planes a = material (200u * (std::size_t) H, 2, hops, H, 2, 0.03, 0.3, 31337u);
        Planes b = a;
        for (auto& ch : b) for (auto& v : ch) v = v * 2.0f;
        const Run ra = runIt (a, par, 48000.0, 1024, { 1024 });
        const Run rb = runIt (b, par, 48000.0, 1024, { 1024 });
        test::ok (ra.rep.eventCount > 0, "x2 fixture: there are events to compare");
        test::ok (ra.rep.eventCount == rb.rep.eventCount, "x2 GAIN: the same number of events");
        std::size_t coordBad = 0, dbBad = 0, scaleBad = 0;
        for (std::size_t i = 0; i < ra.rep.events.size() && i < rb.rep.events.size(); ++i)
        {
            const Burst& x = ra.rep.events[i]; const Burst& y = rb.rep.events[i];
            if (x.start != y.start || x.length != y.length || x.peakAt != y.peakAt || x.hops != y.hops) ++coordBad;
            if (bits (x.peakExcessDb) != bits (y.peakExcessDb)) ++dbBad;
            if (bits (x.peakPower * 4.0) != bits (y.peakPower)
                || bits (x.peakBaseline * 4.0) != bits (y.peakBaseline)
                || bits (x.energy * 4.0) != bits (y.energy)) ++scaleBad;
        }
        test::ok (coordBad == 0, "x2 GAIN: every coordinate identical");
        test::ok (dbBad == 0, "x2 GAIN: every excess dB BITWISE identical");
        test::ok (scaleBad == 0, "x2 GAIN: every raw power is EXACTLY 4x — it is measuring a ratio");
        std::size_t ioiBad = 0;
        for (std::size_t i = 0; i < ra.rep.ioi.size(); ++i)
            if (ra.rep.ioi[i] != rb.rep.ioi[i] || ra.rep.lag[i] != rb.rep.lag[i]) ++ioiBad;
        test::ok (ioiBad == 0, "x2 GAIN: both periodicity histograms identical");
    }

    // (3) The run-up. A programme that is bright from its very first sample must not report the run-up
    // as a burst: before the ring is full nothing is judged at all.
    {
        std::mt19937 rng (5u);
        std::uniform_real_distribution<float> u (-0.6f, 0.6f);
        Planes p (1, std::vector<float> (80u * (std::size_t) H, 0.0f));
        for (auto& v : p[0]) v = u (rng);
        const Run r = runIt (p, par, 48000.0, 999, { 999 });
        test::ok (r.rep.eventCount == 0, "RUN-UP: a bright start is not a burst");
        test::ok (r.rep.eligible == r.rep.hops - 20, "RUN-UP: exactly baselineHops hops were not judged");
        std::size_t early = 0;
        for (const auto& t : r.trace) if (t.full && t.hopIndex < 20 && t.eligible) ++early;
        test::ok (early == 0, "RUN-UP: no hop before the ring was full was eligible");
    }

    // (4) A silent lead-in. More than half the ring is exact zero, so the baseline is zero and the hop is
    // NOT judged — without that gate the first non-zero sample opens an event at any level and publishes
    // +inf dB. This is the fixture for that gate.
    {
        Planes p (1, std::vector<float> (120u * (std::size_t) H, 0.0f));
        std::mt19937 rng (6u);
        std::uniform_real_distribution<float> u (-0.001f, 0.001f);          // -60 dBFS, well above subnormal
        for (std::size_t i = 60u * (std::size_t) H; i < p[0].size(); ++i) p[0][i] = u (rng);
        const Run r = runIt (p, par, 48000.0, 1024, { 1024 });
        test::ok (r.rep.zeroBase > 0, "SILENT LEAD-IN: zero-baseline hops are counted");
        bool fin = true;
        for (const auto& e : r.rep.events) fin = fin && std::isfinite (e.peakExcessDb);
        test::ok (fin, "SILENT LEAD-IN: no event carries a non-finite excess");
        std::size_t judgedInSilence = 0;
        for (const auto& t : r.trace) if (t.full && t.hopIndex >= 20 && t.hopIndex < 55 && t.eligible) ++judgedInSilence;
        test::ok (judgedInSilence == 0, "SILENT LEAD-IN: nothing inside the silence was judged");
    }
}

static void groupSemantics()
{
    test::group ("what the numbers mean");
    const Params par = fastParams();
    const int H = 480, N = 20;

    // A PINNED dB. A hop at exactly 4x its baseline is 6.0206 dB, not 12.04 — these are POWERS.
    // Built by construction: a constant-energy bed, then one hop whose amplitude is doubled, which is
    // exactly 4x the power. The ring holds 20 identical bed hops, so the median IS the bed.
    {
        Planes p (1, std::vector<float> (60u * (std::size_t) H, 0.0f));
        for (std::size_t i = 0; i < p[0].size(); ++i)
            p[0][i] = (float) (0.25 * std::sin (kTau * 7000.0 * (double) i / 48000.0));
        for (std::size_t i = 30u * (std::size_t) H; i < 33u * (std::size_t) H; ++i) p[0][i] *= 3.0f;
        const Run r = runIt (p, par, 48000.0, 1024, { 1024 });
        test::ok (r.rep.eventCount == 1, "a 3x-amplitude passage is one event ("
                  + std::to_string (r.rep.eventCount) + ")");
        if (r.rep.events.size() == 1)
        {
            const Burst& e = r.rep.events[0];
            // x3 amplitude is x9 POWER: 10*log10(9) = 9.54 dB. A 20*log10 mutant on a power ratio
            // reads 19.08, so this one number is the gate on the whole /10-versus-/20 question.
            test::approx (e.peakExcessDb, 9.542, 0.5, "a 9x POWER hop reads 9.54 dB (not 19.08)");
            test::ok (e.start == 30 * H, "the event starts on its hop boundary");
            test::ok (e.hops == 3, "the event is three hops long");
            test::ok (e.length == 3 * H, "its length is three hops of samples");
        }
    }

    // THE ADAPTATION CAP IS A COUNTING THEOREM: a step up that never comes back produces one event of
    // exactly baselineHops/2 hops, because at that point the median of the ring IS the burst level.
    {
        Planes p (1, std::vector<float> (120u * (std::size_t) H, 0.0f));
        for (std::size_t i = 0; i < p[0].size(); ++i)
            p[0][i] = (float) (0.05 * std::sin (kTau * 7000.0 * (double) i / 48000.0));
        for (std::size_t i = 40u * (std::size_t) H; i < p[0].size(); ++i) p[0][i] *= 10.0f;  // +20 dB
        const Run r = runIt (p, par, 48000.0, 1024, { 1024 }, false);
        test::ok (r.rep.eventCount == 1, "a permanent step up is ONE event ("
                  + std::to_string (r.rep.eventCount) + ")");
        if (r.rep.events.size() == 1)
            test::ok (r.rep.events[0].hops == N / 2,
                      "the step's event is exactly baselineHops/2 = " + std::to_string (N / 2)
                      + " hops (got " + std::to_string (r.rep.events[0].hops) + ") — the level-step signature");
    }

    // HYSTERESIS. One dipping hop inside one burst must not split it, because a split invents an
    // inter-onset interval that never happened.
    {
        Planes p (1, std::vector<float> (80u * (std::size_t) H, 0.0f));
        for (std::size_t i = 0; i < p[0].size(); ++i)
            p[0][i] = (float) (0.05 * std::sin (kTau * 7000.0 * (double) i / 48000.0));
        for (std::size_t k = 30; k < 36; ++k)
        {
            // 1.7x amplitude is 2.89x POWER = 4.6 dB: BELOW the 6 dB enter ratio (3.981) and ABOVE the
            // 3 dB exit ratio (1.995). That gap is the whole point of the dip — a value above enter
            // would be held by one threshold as well as two, and the fixture would prove nothing.
            const double g = (k == 33) ? 1.7 : 6.0;
            for (std::size_t i = k * (std::size_t) H; i < (k + 1) * (std::size_t) H; ++i) p[0][i] *= (float) g;
        }
        const Run r = runIt (p, par, 48000.0, 1024, { 1024 }, false);
        test::ok (r.rep.eventCount == 1, "HYSTERESIS: a dipping hop does not split the burst ("
                  + std::to_string (r.rep.eventCount) + " events)");
        test::ok (r.rep.intervals == 0, "HYSTERESIS: and so it invents no inter-onset interval");
    }

    // THE PEAK TIE GOES TO THE EARLIER HOP. Two hops of identical energy against an identical baseline.
    {
        Planes p (1, std::vector<float> (80u * (std::size_t) H, 0.0f));
        for (std::size_t i = 0; i < p[0].size(); ++i)
            p[0][i] = (float) (0.05 * std::sin (kTau * 6000.0 * (double) i / 48000.0));
        for (std::size_t k = 40; k < 43; ++k)
            for (std::size_t i = k * (std::size_t) H; i < (k + 1) * (std::size_t) H; ++i) p[0][i] *= 4.0f;
        const Run r = runIt (p, par, 48000.0, 1024, { 1024 });
        const Trace* t41 = nullptr; const Trace* t42 = nullptr;
        for (const auto& t : r.trace)
        {
            if (t.hopIndex == 41) t41 = &t;
            if (t.hopIndex == 42) t42 = &t;
        }
        const bool tied = t41 != nullptr && t42 != nullptr
                       && bits (t41->energy) == bits (t42->energy)
                       && bits (t41->baseline) == bits (t42->baseline);
        test::ok (tied, "tie fixture: hops 41 and 42 really are a BITWISE tie in energy and baseline");
        if (r.rep.events.size() >= 1)
        {
            const Burst& e = r.rep.events[0];
            test::ok (e.hops == 3, "tie fixture: the event spans all three boosted hops");
            test::ok (e.peakAt == 41 * H, "a peak TIE keeps the EARLIER hop (peakAt "
                      + std::to_string (e.peakAt) + ", want " + std::to_string (41 * H) + ")");
        }
        else test::ok (false, "tie fixture produced no event");
    }

    // THE WIDEBAND COORDINATE. The same band burst, once alone and once with a loud out-of-band bed
    // under it: the band numbers must be near-identical and peakWidePower must be far apart. That is the
    // question the ratio cannot answer, and the reason the field exists.
    {
        Planes a (1, std::vector<float> (80u * (std::size_t) H, 0.0f));
        for (std::size_t i = 0; i < a[0].size(); ++i)
            a[0][i] = (float) (0.02 * std::sin (kTau * 7000.0 * (double) i / 48000.0));
        for (std::size_t i = 40u * (std::size_t) H; i < 41u * (std::size_t) H; ++i) a[0][i] *= 8.0f;
        Planes b = a;
        for (std::size_t i = 0; i < b[0].size(); ++i)                     // 200 Hz, far below the band
            b[0][i] = (float) ((double) b[0][i] + 0.5 * std::sin (kTau * 200.0 * (double) i / 48000.0));
        const Run ra = runIt (a, par, 48000.0, 1024, { 1024 }, false);
        const Run rb = runIt (b, par, 48000.0, 1024, { 1024 }, false);
        test::ok (ra.rep.eventCount == 1 && rb.rep.eventCount == 1, "wideband fixture: one event each");
        if (ra.rep.events.size() == 1 && rb.rep.events.size() == 1)
        {
            test::approx (rb.rep.events[0].peakExcessDb, ra.rep.events[0].peakExcessDb, 0.5,
                          "out-of-band energy does not change the BAND excess");
            test::ok (rb.rep.events[0].peakWidePower > ra.rep.events[0].peakWidePower * 10.0,
                      "but peakWidePower sees it — the band's share is a second coordinate");
        }
    }
}

static void groupPeriodicity()
{
    test::group ("periodicity — counts, not a verdict");
    const Params par = fastParams();
    const int H = 480;

    // A HI-HAT: a burst every 12 hops. The adjacent histogram must have its mode at 12 and hold nearly
    // all of its mass there; the autocorrelation must comb at every multiple of 12.
    {
        std::vector<int> hops;
        for (int k = 25; k < 400; k += 12) hops.push_back (k);
        const Planes p = material (420u * (std::size_t) H, 1, hops, H, 1, 0.01, 0.30, 808u);
        const Run r = runIt (p, par, 48000.0, 2048, { 2048 }, false);
        test::ok (r.rep.eventCount > 20, "hat train: the hats were found ("
                  + std::to_string (r.rep.eventCount) + " events)");
        test::ok (r.rep.modal == 12, "hat train: the modal spacing is 12 hops (got "
                  + std::to_string (r.rep.modal) + ")");
        test::ok (r.rep.modalMass * 10 >= r.rep.intervals * 9,
                  "hat train: >=90% of the spacings sit at the mode +/-1");
        // The comb: lags 12, 24, 36 must all carry pairs, and 6 and 18 must not.
        test::ok (r.rep.lag[11] > 0 && r.rep.lag[23] > 0 && r.rep.lag[35] > 0,
                  "hat train: the autocorrelation combs at 12, 24, 36");
        test::ok (r.rep.lag[5] == 0 && r.rep.lag[17] == 0,
                  "hat train: and carries nothing at 6 or 18");
    }

    // A HAT WITH DROPOUTS — every 5th hit removed. The adjacent histogram smears into 12 and 24; the
    // autocorrelation keeps its peak at 12. This is the case that made the second histogram worth its
    // 4 KB, and it is why "periodicity" is not the adjacent histogram alone.
    {
        std::vector<int> hops; int n = 0;
        for (int k = 25; k < 400; k += 12) { if (++n % 5 != 0) hops.push_back (k); }
        const Planes p = material (420u * (std::size_t) H, 1, hops, H, 1, 0.01, 0.30, 809u);
        const Run r = runIt (p, par, 48000.0, 2048, { 2048 }, false);
        test::ok (r.rep.ioi[23] > 0, "dropout: the adjacent histogram picks up a spurious 24");
        test::ok (r.rep.lag[11] > r.rep.lag[23], "dropout: the autocorrelation still peaks at 12");
    }

    // TWO ISOLATED FLASHES, 100 hops apart. The mode carries 100% of the mass and means NOTHING — which
    // is why the core publishes the COUNT beside it and no score.
    {
        const Planes p = material (300u * (std::size_t) H, 1, { 30, 130 }, H, 1, 0.01, 0.30, 810u);
        const Run r = runIt (p, par, 48000.0, 2048, { 2048 }, false);
        test::ok (r.rep.intervals <= 1, "two flashes give at most one interval ("
                  + std::to_string (r.rep.intervals) + ")");
        test::ok (r.rep.modalMass == r.rep.intervals,
                  "and the mode holds all of it — the count is what says not to trust it");
    }

    // A spacing past the histogram's reach goes to the OVERFLOW count, not into the last bin.
    {
        const Planes p = material (1200u * (std::size_t) H, 1, { 30, 1000 }, H, 1, 0.01, 0.30, 811u);
        const Run r = runIt (p, par, 48000.0, 4096, { 4096 }, false);
        test::ok (r.rep.ioiOverflow == 1, "a 970-hop spacing lands in the overflow count");
        test::ok (r.rep.ioi[BB::kIoiBins - 1] == 0, "and NOT folded into the last bin");
    }
}

static void groupContract()
{
    test::group ("the law-11 contract, storage and lifecycle");
    Params par = fastParams();
    BB d;
    d.setParams (par);

    // Refusals. storageFor() must say no to exactly what prepare() says no to.
    struct Bad { const char* what; double fs; int nch; Params p; };
    std::vector<Bad> bad;
    { Params p = par; bad.push_back ({ "sampleRate 0",          0.0, 2, p }); }
    { Params p = par; bad.push_back ({ "sampleRate NaN",        std::numeric_limits<double>::quiet_NaN(), 2, p }); }
    { Params p = par; bad.push_back ({ "sampleRate too low",    999.0, 2, p }); }
    { Params p = par; bad.push_back ({ "channels 0",            48000.0, 0, p }); }
    { Params p = par; bad.push_back ({ "channels too many",     48000.0, core::kMaxChannels + 1, p }); }
    { Params p = par; p.bandHighHz = p.bandLowHz;      bad.push_back ({ "band hi == lo",      48000.0, 2, p }); }
    { Params p = par; p.bandLowHz  = 0.5;              bad.push_back ({ "band lo below 1 Hz", 48000.0, 2, p }); }
    { Params p = par; p.bandHighHz = 9000.0;           bad.push_back ({ "band above 0.49 fs", 16000.0, 2, p }); }
    { Params p = par; p.hopMs      = 0.4;              bad.push_back ({ "hopMs below the floor", 48000.0, 2, p }); }
    { Params p = par; p.hopMs      = 2000.0;           bad.push_back ({ "hopMs above the cap", 48000.0, 2, p }); }
    { Params p = par; p.hopMs      = std::numeric_limits<double>::quiet_NaN();
                                                       bad.push_back ({ "hopMs NaN",          48000.0, 2, p }); }
    { Params p = par; p.baselineMs = 5.0;              bad.push_back ({ "baselineMs < hopMs", 48000.0, 2, p }); }
    { Params p = par; p.baselineMs = std::numeric_limits<double>::quiet_NaN();
                                                       bad.push_back ({ "baselineMs NaN",     48000.0, 2, p }); }
    { Params p = par; p.enterDb    = 0.0;              bad.push_back ({ "enterDb 0",          48000.0, 2, p }); }
    { Params p = par; p.enterDb    = std::numeric_limits<double>::quiet_NaN();
                                                       bad.push_back ({ "enterDb NaN",        48000.0, 2, p }); }
    { Params p = par; p.exitDb     = 9.0;              bad.push_back ({ "exitDb > enterDb",   48000.0, 2, p }); }
    { Params p = par; p.maxEvents  = -1;               bad.push_back ({ "maxEvents negative", 48000.0, 2, p }); }
    { Params p = par; p.maxEvents  = BB::kMaxEventsLimit + 1;
                                                       bad.push_back ({ "maxEvents too many", 48000.0, 2, p }); }
    std::size_t leak = 0;
    for (const auto& b : bad)
    {
        const BB::Storage st = BB::storageFor (b.fs, b.nch, b.p);
        BB e; e.setParams (b.p);
        const bool accepted = e.prepare (b.fs, 512, b.nch);
        if (st.ok || accepted || st.bytes() != 0) { ++leak; test::ok (false, std::string ("must refuse: ") + b.what); }
    }
    test::ok (leak == 0, std::to_string (bad.size()) + " malformed configurations refused by BOTH storageFor and prepare");

    // The published storage is the allocated storage, and prepare() is the only allocator. Counted here
    // rather than derived: `okNoAlloc` is only strict on libc++, which is this row.
    {
        BB e; e.setParams (par);
        const BB::Storage st = BB::storageFor (48000.0, 2, par);
        test::ok (st.ok && st.bytes() > 0, "storageFor publishes a demand BEFORE prepare()");
        test::ok (st.hopSamples == 480 && st.baselineHops == 20, "and it publishes the geometry it derived");
        test::run (e.prepare (48000.0, 512, 2));
        Planes p (2, std::vector<float> (40000, 0.01f));
        const auto v = ptrs (p, 0);
        const long before = g_allocs.load();
        test::run (e.process (v.data(), 2, 40000));
        e.finish();
        const long after = g_allocs.load();
        test::okNoAlloc (after == before, "process() and finish() allocate NOTHING");
    }

    // Law 11 order: malformed -> unprepared -> finished -> too many channels -> n == 0 -> run.
    {
        BB e;
        Planes p (2, std::vector<float> (1000, 0.0f));
        const auto v = ptrs (p, 0);
        test::ok (! e.process (v.data(), 2, 100), "unprepared process() is refused");
        e.setParams (par);
        test::run (e.prepare (48000.0, 512, 2));
        test::ok (! e.process (v.data(), -1, 100), "a negative channel count is refused");
        test::ok (! e.process (v.data(), 2, -1), "a negative length is refused");
        test::ok (! e.process (v.data(), 3, 100), "more channels than prepared is refused");
        test::run (e.process (v.data(), 2, 0));
        test::ok (e.samplesProcessed() == 0, "n == 0 changes nothing and is accepted");
        test::run (e.process (v.data(), 2, 1000));
        // A legal call LONGER than maxBlock is consumed IN FULL.
        Planes q (2, std::vector<float> (5000, 0.0f));
        const auto w = ptrs (q, 0);
        test::run (e.process (w.data(), 2, 5000));
        test::ok (e.samplesProcessed() == 6000, "a call longer than maxBlock is consumed in full");
        e.finish();
        test::ok (! e.process (v.data(), 2, 100), "process() after finish() is refused");
        const Report a = snapshot (e);
        e.finish(); e.finish();
        test::ok (sameReport (a, snapshot (e)), "finish() is IDEMPOTENT — twice more changes no bit");
        e.reset();
        test::ok (e.samplesProcessed() == 0 && e.hopCount() == 0 && e.eventCount() == 0
                  && e.damagedHops() == 0 && e.onsetCount() == 0 && ! e.isFinished(),
                  "reset() re-anchors the clock, the ring, the events and the histograms");
        test::run (e.process (v.data(), 2, 1000));
        test::ok (e.samplesProcessed() == 1000, "and process() works again after reset()");
    }

    // reset() really is a fresh start: the same programme twice around a reset() is bit-identical.
    {
        const Planes p = material (100u * 480u, 2, { 25, 60 }, 480, 2, 0.02, 0.3, 4u);
        const Run one = runIt (p, par, 48000.0, 1024, { 1024 });
        BB e; e.setParams (par); test::run (e.prepare (48000.0, 1024, 2));
        const auto v = ptrs (p, 0);
        test::run (e.process (v.data(), 2, (int) p[0].size()));
        e.finish();
        e.reset();
        test::run (e.process (v.data(), 2, (int) p[0].size()));
        e.finish();
        test::ok (sameReport (one.rep, snapshot (e)), "a run AFTER reset() equals a run from prepare()");
    }
}

static void groupEdges()
{
    test::group ("edges — empty, short, the tail, holes, a missing channel, exhaustion");
    const Params par = fastParams();
    const int H = 480, N = 20;

    // T = 0, 1, W-1, W, W+H-1, W+H. Nothing may be judged before hop N closes.
    for (std::size_t T : { (std::size_t) 0, (std::size_t) 1, (std::size_t) (H * N - 1),
                           (std::size_t) (H * N), (std::size_t) (H * N + H - 1), (std::size_t) (H * N + H) })
    {
        // In-band content, not DC: a constant 0.1f has no energy in 5-9 kHz, so the baseline would be
        // zero and this sweep would measure the zero-baseline gate instead of the hop geometry.
        Planes p (1, std::vector<float> (T, 0.0f));
        for (std::size_t i = 0; i < T; ++i)
            p[0][i] = (float) (0.2 * std::sin (kTau * 7000.0 * (double) i / 48000.0));
        BB e; e.setParams (par);
        test::run (e.prepare (48000.0, 256, 1));
        if (T > 0) { const auto v = ptrs (p, 0); test::run (e.process (v.data(), 1, (int) T)); }
        e.finish();
        const std::int64_t wantHops = (std::int64_t) (T / (std::size_t) H);
        const bool wantJudged = wantHops > N;
        test::ok (e.hopCount() == wantHops && e.samplesProcessed() == (std::int64_t) T,
                  "T = " + std::to_string (T) + ": the hop count and the clock are exact");
        test::ok ((e.eligibleHops() > 0) == wantJudged,
                  "T = " + std::to_string (T) + ": judged only once hop N has closed");
        test::ok (e.eventsValid() == wantJudged,
                  "T = " + std::to_string (T) + ": eventsValid tracks that");
        if (! wantJudged)
            test::ok (e.eventsInvalidReason() == analysis::BandBurstsInvalid::NoEligibleHop,
                      "T = " + std::to_string (T) + ": and the reason is NoEligibleHop");
        test::ok (e.tailPartialSamples() == (std::int64_t) (T % (std::size_t) H),
                  "T = " + std::to_string (T) + ": the uncovered tail is named exactly");
    }

    // The tail is NAMED and does not enter the ring or open an event. A burst entirely inside the final
    // partial hop must therefore add no event, while its energy is still published.
    {
        Planes p (1, std::vector<float> (60u * (std::size_t) H + 200u, 0.0f));
        for (std::size_t i = 0; i < p[0].size(); ++i)
            p[0][i] = (float) (0.02 * std::sin (kTau * 7000.0 * (double) i / 48000.0));
        for (std::size_t i = 60u * (std::size_t) H; i < p[0].size(); ++i) p[0][i] *= 30.0f;
        const Run r = runIt (p, par, 48000.0, 1024, { 1024 }, false);
        test::ok (r.rep.tailSamples == 200, "the partial tail is 200 samples");
        test::ok (r.rep.tailEnergy > 0.0 && std::isfinite (r.rep.tailEnergy),
                  "its energy is published as a finite number");
        test::ok (r.rep.eventCount == 0, "a burst living only in the partial tail opens NO event");
        test::ok (r.rep.hops == 60, "and the tail is not counted as a hop");
    }

    // Digital silence throughout: no events, no damage, and the baseline is zero for every judged hop.
    {
        Planes p (2, std::vector<float> (60u * (std::size_t) H, 0.0f));
        const Run r = runIt (p, par, 48000.0, 1024, { 1024 }, false);
        test::ok (r.rep.eventCount == 0, "silence gives no events");
        test::ok (r.rep.damaged == 0, "silence damages nothing");
        test::ok (r.rep.zeroBase == r.rep.hops - N, "and every judged hop had a zero baseline");
        test::ok (bits (r.rep.chEnergy[0]) == bits (0.0), "the band energy of silence is exactly zero");
    }

    // Non-finite input: counted per channel, first coordinate recorded, hop marked damaged, and the
    // programme-energy flag goes invalid with a reason. The filter never sees it.
    {
        Planes p (2, std::vector<float> (60u * (std::size_t) H, 0.05f));
        p[0][30u * (std::size_t) H + 5u] = std::numeric_limits<float>::quiet_NaN();
        p[1][30u * (std::size_t) H + 6u] = -std::numeric_limits<float>::infinity();
        const Run r = runIt (p, par, 48000.0, 1024, { 1024 });
        test::ok (r.rep.chNonFinite[0] == 1 && r.rep.chNonFinite[1] == 1, "non-finite samples counted per channel");
        test::ok (r.rep.firstNonFinite == 30 * H + 5, "the first coordinate is recorded exactly");
        test::ok (r.rep.damaged == 1, "exactly one hop is damaged");
        test::ok (! r.rep.energyValid, "the programme energy is marked invalid");
        bool fin = true; for (double v : r.rep.chEnergy) fin = fin && std::isfinite (v);
        test::ok (fin, "and every per-channel total stays finite");
        std::size_t dmg = 0;
        for (const auto& t : r.trace) if (t.damaged) ++dmg;
        test::ok (dmg == 1, "the trace marks that one hop and no other");
    }

    // A FINITE input that overflows the filter. This is the case input sanitising cannot cover, and the
    // fixture the denormal-flush mutant needs.
    {
        Planes p (1, std::vector<float> (60u * (std::size_t) H, 0.05f));
        for (int j = 0; j < 16; ++j)
            p[0][30u * (std::size_t) H + (std::size_t) j] = (j % 2 == 0) ? 3.0e38f : -3.0e38f;
        const Run r = runIt (p, par, 48000.0, 1024, { 1024 }, false);
        test::ok (r.rep.chNonFinite[0] == 0, "every input sample was finite");
        test::ok (r.rep.overflowSamples > 0, "yet the FILTER overflowed and it is counted");
        test::ok (! r.rep.energyValid, "which invalidates the programme energy too");
        test::ok (std::isfinite (r.rep.chEnergy[0]), "and the accumulator was never poisoned");
    }

    // A channel that DISAPPEARS MID-STREAM is handed the silence it is receiving (law 11a), which must be
    // bit-identical to feeding it explicit zeros over the same samples. It has to leave mid-stream and
    // come back: a channel absent from sample 0 has nothing in its filters, so "driven with zeros" and
    // "skipped entirely" are the same state and the fixture would be blind. Here channel 1 is loud first,
    // so on leaving its four SVF columns must ring DOWN through the sum — and on returning it must not
    // replay a frozen tail.
    {
        const int T = 60 * H, A = 25 * H, B = 35 * H;
        Planes two = material ((std::size_t) T, 2, { 25, 40 }, H, 2, 0.02, 0.3, 21u);
        Planes zeroed = two;
        for (int i = A; i < B; ++i) zeroed[1][(std::size_t) i] = 0.0f;   // explicit zeros over the gap

        Cap ca; ca.v.reserve (1u << 12);
        BB a; a.setParams (par); a.setObserver (&observe, &ca);
        test::run (a.prepare (48000.0, 4096, 2));
        const auto p0 = ptrs (two, 0);
        const auto pA = ptrs (two, (std::size_t) A);
        const auto pB = ptrs (two, (std::size_t) B);
        test::run (a.process (p0.data(), 2, A));                  // both channels
        test::run (a.process (pA.data(), 1, B - A));              // channel 1 GONE
        test::run (a.process (pB.data(), 2, T - B));              // and back
        a.finish();
        const Report ra = snapshot (a);

        const Run b = runIt (zeroed, par, 48000.0, 4096, { 4096 });
        test::ok (ra.chAbsent[1] == B - A, "an absent channel's samples are counted exactly");
        test::ok (ra.chAbsent[0] == 0, "and the present channel counts none");
        // Everything but the absence counters themselves: those MUST differ, since recording the
        // absence is their whole job. Zeroed in a copy rather than skipped in the comparator, so the
        // one field this fixture is allowed to move is named right here.
        Report raX = ra, rbX = b.rep;
        for (auto& v : raX.chAbsent) v = 0;
        for (auto& v : rbX.chAbsent) v = 0;
        test::ok (sameReport (raX, rbX), "ABSENT == ZERO-FED, the whole report bitwise but for chAbsent");
        test::ok (ca.v.size() == b.trace.size(), "and the same number of hops closed");
        std::size_t tBad = 0;
        for (std::size_t i = 0; i < ca.v.size() && i < b.trace.size(); ++i)
            if (! sameTrace (ca.v[i], b.trace[i])) ++tBad;
        test::ok (tBad == 0, "every traced hop is bitwise the same (" + std::to_string (tBad) + " differ)");
    }

    // CAPACITY EXHAUSTION IS DATA. A tiny list keeps a prefix, says it is incomplete, and the count and
    // BOTH histograms go on — the previous onset is a scalar, so the spacings never notice.
    {
        Params small = par; small.maxEvents = 3;
        std::vector<int> hops;
        for (int k = 25; k < 300; k += 10) hops.push_back (k);
        const Planes p = material (320u * (std::size_t) H, 1, hops, H, 1, 0.01, 0.30, 55u);
        const Run big = runIt (p, par, 48000.0, 2048, { 2048 }, false);
        const Run r   = runIt (p, small, 48000.0, 2048, { 2048 }, false);
        test::ok (r.rep.eventCount == big.rep.eventCount,
                  "the event COUNT is unaffected by the capacity");
        test::ok (r.rep.stored == 3 && ! r.rep.eventsComplete, "a 3-entry list keeps 3 and says so");
        test::ok (big.rep.eventsComplete, "while the roomy run reports complete");
        test::ok (r.rep.intervals == big.rep.intervals, "the spacings kept counting past exhaustion");
        std::size_t hBad = 0;
        for (std::size_t i = 0; i < r.rep.ioi.size(); ++i)
            if (r.rep.ioi[i] != big.rep.ioi[i] || r.rep.lag[i] != big.rep.lag[i]) ++hBad;
        test::ok (hBad == 0, "and BOTH histograms are identical to the roomy run");
        std::size_t pBad = 0;
        for (std::size_t i = 0; i < r.rep.events.size(); ++i)
            if (! sameBurst (r.rep.events[i], big.rep.events[i])) ++pBad;
        test::ok (pBad == 0, "the kept prefix is bitwise the roomy run's prefix");
        // maxEvents == 0: nothing stored, nothing read, no crash.
        Params none = par; none.maxEvents = 0;
        const Run z = runIt (p, none, 48000.0, 2048, { 2048 }, false);
        test::ok (z.rep.eventCount == big.rep.eventCount && z.rep.stored == 0 && ! z.rep.eventsComplete,
                  "maxEvents == 0 counts everything and stores nothing");
        BB e; e.setParams (none); test::run (e.prepare (48000.0, 256, 1));
        test::ok (e.event (0).length == 0 && e.event (-1).length == 0,
                  "and every accessor stays total");
    }

    // The defaults must be usable at the rates that matter, and refused where the band cannot exist.
    {
        for (double fs : { 44100.0, 48000.0, 88200.0, 96000.0, 192000.0 })
        {
            BB e; e.setParams (Params {});
            test::ok (e.prepare (fs, 512, 2), "the DEFAULT band prepares at " + std::to_string ((int) fs) + " Hz");
        }
        BB e; e.setParams (Params {});
        test::ok (! e.prepare (16000.0, 512, 2), "and is refused at 16 kHz, where 9 kHz is past 0.49 fs");
        BB f; f.setParams (Params {});
        test::run (f.prepare (48000.0, 512, 2));
        test::approx (f.bandLowHz(), 5000.0, 0.01, "the EFFECTIVE low corner is published");
        test::approx (f.bandHighHz(), 9000.0, 0.01, "the EFFECTIVE high corner is published");
        test::ok (f.hopSamples() == 480 && f.baselineHops() == 200, "the default geometry is 480 / 200");
        test::ok (f.latencySamples() == 0, "a read-only sink has no latency");
    }
}


// ===================================================================== the post-build review round
// Everything in this group came from the two crew rounds run AFTER the module was first green. The
// first two checks are the important ones: an ANALYTIC band gain and an ANALYTIC median convention,
// both computed outside this program. The oracle above shares the filter code and the median-index
// convention with the engine, so a mistake in either would be invisible to it — these two numbers are
// the outside witnesses that close that hole.
static void groupOutsideWitnesses()
{
    test::group ("outside witnesses + the review round's findings");
    const Params par = fastParams();
    const int H = 480;

    // (1) THE BAND GAIN, ANALYTICALLY. A unit 7 kHz sine at 48 kHz is exactly 60 cycles per 480-sample
    // hop, so the hop's wideband power is exactly 480 * 1/2 = 240. The cascade HP4(5k)^2 * LP4(9k)^2,
    // evaluated from the BLT prewarp the SVF implements, has |H(7k)|^2 = 0.3966383155110891 (-4.0161 dB),
    // so a settled hop's BAND energy must be 240 * that = 95.1932. Neither number comes from this code.
    {
        Planes p (1, std::vector<float> (60u * (std::size_t) H, 0.0f));
        for (std::size_t i = 0; i < p[0].size(); ++i)
            p[0][i] = (float) std::sin (kTau * 7000.0 * (double) i / 48000.0);
        const Run r = runIt (p, par, 48000.0, 1024, { 1024 });
        const Trace* t = nullptr;
        for (const auto& e : r.trace) if (e.hopIndex == 40) t = &e;      // long settled
        test::ok (t != nullptr, "band-gain fixture: hop 40 was traced");
        if (t != nullptr)
        {
            test::approx (t->wideEnergy, 240.0, 0.05,
                          "a unit 7 kHz sine's WIDEBAND hop power is exactly 480/2 = 240");
            test::approx (t->energy, 95.1932, 0.5,
                          "and its BAND energy is 240*|H(7k)|^2 = 95.1932 — the analytic cascade gain");
        }
    }

    // (2) THE MEDIAN CONVENTION, ANALYTICALLY. Ten quiet hops and ten loud ones fill a 20-hop ring, so
    // the ordinary median would be the mean of ranks 9 and 10 while the UPPER median (sorted index
    // N/2 = 10) is the LOUD value. Alternating hops put exactly ten of each in the ring; the traced
    // baseline must equal a loud hop's energy, not a quiet one's and not anything in between.
    {
        Planes p (1, std::vector<float> (80u * (std::size_t) H, 0.0f));
        for (std::size_t k = 0; k < 80; ++k)
        {
            const double a = (k % 2 == 0) ? 0.02 : 0.20;                  // 100x in power
            for (std::size_t i = k * (std::size_t) H; i < (k + 1) * (std::size_t) H; ++i)
                p[0][i] = (float) (a * std::sin (kTau * 6000.0 * (double) i / 48000.0));
        }
        const Run r = runIt (p, par, 48000.0, 1024, { 1024 });
        const Trace* dec = nullptr; const Trace* loud = nullptr;
        for (const auto& e : r.trace)
        {
            if (e.hopIndex == 40) dec = &e;                               // ring holds hops 20..39
            if (e.hopIndex == 39) loud = &e;                              // odd hop => the loud level
        }
        test::ok (dec != nullptr && loud != nullptr, "median fixture: both hops traced");
        if (dec != nullptr && loud != nullptr)
            test::ok (bits (dec->baseline) == bits (loud->energy),
                      "the baseline is the UPPER median — a LOUD hop's energy, bit for bit");
    }

    // (3) A MOVED-FROM detector must not stay armed. The implicit move steals the vectors and copies the
    // flag, so without deriving armedness from the storage this indexes an empty channel list.
    {
        BB a; a.setParams (par);
        test::run (a.prepare (48000.0, 512, 2));
        BB b = std::move (a);
        Planes p (2, std::vector<float> (64, 0.1f));
        const auto v = ptrs (p, 0);
        test::ok (! a.process (v.data(), 2, 64), "a MOVED-FROM detector refuses process()");
        a.finish();                                                        // must not touch anything
        test::ok (a.samplesProcessed() == 0, "and finish() on it does nothing");
        test::run (b.process (v.data(), 2, 64));
        test::ok (b.samplesProcessed() == 64, "while the moved-TO detector works");
    }

    // (4) The refusals the review found were missing: an unbounded baselineMs reached llround out of
    // range, a non-finite enterDb passed `> 0` and disabled detection behind a valid-looking report, and
    // the corners were validated as doubles while the FILTER is handed floats.
    {
        const double inf = std::numeric_limits<double>::infinity();
        struct Bad { const char* what; double fs; Params p; };
        std::vector<Bad> bad;
        { Params p = par; p.baselineMs = 1e30;  bad.push_back ({ "baselineMs 1e30 (llround out of range)", 48000.0, p }); }
        { Params p = par; p.baselineMs = inf;   bad.push_back ({ "baselineMs +inf", 48000.0, p }); }
        { Params p = par; p.baselineMs = std::numeric_limits<double>::max();
                                               bad.push_back ({ "baselineMs DBL_MAX", 48000.0, p }); }
        { Params p = par; p.enterDb = inf;     bad.push_back ({ "enterDb +inf (detection silently off)", 48000.0, p }); }
        { Params p = par; p.enterDb = 1.0e6;   bad.push_back ({ "enterDb past the cap", 48000.0, p }); }
        { Params p = par; p.bandLowHz = 7000.0; p.bandHighHz = std::nextafter (7000.0, inf);
                          bad.push_back ({ "corners that collapse to one float", 48000.0, p }); }
        { Params p = par; p.bandHighHz = 5000.0001;
                          bad.push_back ({ "a band narrower than one float step", 48000.0, p }); }
        // 0.49*44101 = 21609.49, whose float is 21609.490234375 — just PAST the cap, so Svf would clamp
        // behind the accessor's back. Validated on the float now, so it is refused.
        { Params p = par; p.bandHighHz = 0.49 * 44101.0;
                          bad.push_back ({ "a corner whose FLOAT rounds past 0.49 fs", 44101.0, p }); }
        std::size_t leak = 0;
        for (const auto& b : bad)
        {
            const BB::Storage st = BB::storageFor (b.fs, 2, b.p);
            BB e; e.setParams (b.p);
            if (st.ok || e.prepare (b.fs, 512, 2)) { ++leak; test::ok (false, std::string ("must refuse: ") + b.what); }
        }
        test::ok (leak == 0, std::to_string (bad.size()) + " more malformed configurations refused");
        // And a band that IS a band, one float step wider than the collapsing one, still prepares.
        BB e; Params good = par; good.bandLowHz = 7000.0; good.bandHighHz = 7100.0;
        e.setParams (good);
        test::ok (e.prepare (48000.0, 512, 1), "a 100 Hz band at 7 kHz still prepares");
        test::ok (e.bandLowHz() < e.bandHighHz(), "and its published corners are strictly ordered");
    }

    // (5) Damage inside the PARTIAL TAIL must reach the flag of an event that finish() closes at T,
    // because the tail lies inside that event's published extent.
    {
        Planes p (1, std::vector<float> (60u * (std::size_t) H + 200u, 0.0f));
        for (std::size_t i = 0; i < p[0].size(); ++i)
            p[0][i] = (float) (0.02 * std::sin (kTau * 7000.0 * (double) i / 48000.0));
        for (std::size_t i = 55u * (std::size_t) H; i < p[0].size(); ++i) p[0][i] *= 8.0f;   // open a burst
        p[0][60u * (std::size_t) H + 100u] = std::numeric_limits<float>::quiet_NaN();        // in the tail
        const Run r = runIt (p, par, 48000.0, 1024, { 1024 }, false);
        test::ok (r.rep.eventCount >= 1, "tail-damage fixture: an event is open at the end");
        if (! r.rep.events.empty())
        {
            const Burst& e = r.rep.events[0];
            test::ok (e.closedByFinish, "it was closed by finish()");
            test::ok (e.length == r.rep.samples - e.start, "its extent reaches T, covering the tail");
            test::ok (e.touchedNonFinite, "and the tail's damage IS flagged on it");
        }
    }

    // (6) The step cap is ceil(N/2), not N/2 — measured at odd ring sizes, where they differ.
    {
        for (int n : { 1, 3, 5 })
        {
            Params p = par; p.baselineMs = par.hopMs * (double) n;
            Planes m (1, std::vector<float> (60u * (std::size_t) H, 0.0f));
            for (std::size_t i = 0; i < m[0].size(); ++i)
                m[0][i] = (float) (0.02 * std::sin (kTau * 7000.0 * (double) i / 48000.0));
            for (std::size_t i = 30u * (std::size_t) H; i < m[0].size(); ++i) m[0][i] *= 10.0f;
            const Run r = runIt (m, p, 48000.0, 1024, { 1024 }, false);
            const int want = (n + 1) / 2;                                 // ceil(n/2)
            test::ok (! r.rep.events.empty() && r.rep.events[0].hops == want,
                      "baselineHops = " + std::to_string (n) + ": a step lasts ceil(N/2) = "
                      + std::to_string (want) + " hops (got "
                      + std::to_string (r.rep.events.empty() ? -1 : r.rep.events[0].hops) + ")");
        }
    }

    // (7) A RAMP ESCAPES THE CAP. The median chases a rising signal from behind and never catches it, so
    // an event can run far past ceil(N/2). The header used to claim the cap held for every event.
    {
        Planes p (1, std::vector<float> (80u * (std::size_t) H, 0.0f));
        for (std::size_t i = 0; i < p[0].size(); ++i)
            p[0][i] = (float) (0.002 * std::sin (kTau * 7000.0 * (double) i / 48000.0));
        for (std::size_t k = 30; k < 75; ++k)                             // +15 % per hop, for 45 hops
        {
            const double g = std::pow (1.15, (double) (k - 30) + 1.0);
            for (std::size_t i = k * (std::size_t) H; i < (k + 1) * (std::size_t) H; ++i) p[0][i] *= (float) g;
        }
        const Run r = runIt (p, par, 48000.0, 1024, { 1024 }, false);
        test::ok (! r.rep.events.empty(), "ramp fixture: an event opened");
        if (! r.rep.events.empty())
            test::ok (r.rep.events[0].hops > 20 / 2,
                      "a RISING ramp holds an event open past ceil(N/2) = 10 hops (got "
                      + std::to_string (r.rep.events[0].hops) + ") — the cap is a STEP theorem only");
    }

    // (8) THE PROGRAMME-WIDE SUM DILUTES A SINGLE-CHANNEL BURST BY THE CHANNEL COUNT. One channel going
    // 9x in power is (n+8)/n of the total: 5.0 at two channels (over the 3.981 enter ratio) and 1.25 at
    // thirty-two (nowhere near). Named in the header; here it is measured.
    {
        for (int nch : { 2, 32 })
        {
            Planes p ((std::size_t) nch, std::vector<float> (60u * (std::size_t) H, 0.0f));
            for (int c = 0; c < nch; ++c)
                for (std::size_t i = 0; i < p[0].size(); ++i)
                    p[(std::size_t) c][i] = (float) (0.02 * std::sin (kTau * 7000.0 * (double) i / 48000.0));
            for (std::size_t i = 40u * (std::size_t) H; i < 43u * (std::size_t) H; ++i)
                p[0][i] *= 3.0f;                                          // channel 0 alone: 9x power
            const Run r = runIt (p, par, 48000.0, 2048, { 2048 }, false);
            if (nch == 2) test::ok (r.rep.eventCount == 1, "a 9x burst in 1 of 2 channels IS found");
            else          test::ok (r.rep.eventCount == 0,
                                    "the same burst in 1 of 32 channels is DILUTED away ("
                                    + std::to_string (r.rep.eventCount) + " events)");
        }
    }

    // (9) A SPECTRAL SHIFT AT CONSTANT LEVEL is a band burst, and the wideband coordinate is what says
    // so. 1 kHz -> 7 kHz at the same amplitude: the band rises enormously while the programme's level
    // does not move at all. The instrument is right to report it; peakWidePower is how a consumer tells
    // it apart from something getting louder.
    {
        Planes p (1, std::vector<float> (80u * (std::size_t) H, 0.0f));
        for (std::size_t i = 0; i < p[0].size(); ++i)
        {
            const double hz = (i < 40u * (std::size_t) H) ? 1000.0 : 7000.0;
            p[0][i] = (float) (0.1 * std::sin (kTau * hz * (double) i / 48000.0));
        }
        const Run r = runIt (p, par, 48000.0, 1024, { 1024 });
        test::ok (r.rep.eventCount >= 1, "a 1 kHz -> 7 kHz shift at constant level reports a band burst");
        if (! r.rep.events.empty())
        {
            const Burst& e = r.rep.events[0];
            test::ok (e.peakExcessDb > 20.0, "with a large BAND excess ("
                      + std::to_string ((int) e.peakExcessDb) + " dB)");
            // The wideband power of both halves is the same 0.1-amplitude sine: 480*0.005 = 2.4.
            test::approx (e.peakWidePower, 0.005, 0.0005,
                          "while the WIDEBAND mean square is unchanged — the level never moved");
        }
    }
}

int main()
{
    std::printf ("felitronics::analysis::BandBursts\n");
    groupOracle();
    groupSlicing();
    groupNegatives();
    groupSemantics();
    groupPeriodicity();
    groupContract();
    groupEdges();
    groupOutsideWitnesses();
    return test::report();
}

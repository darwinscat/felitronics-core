// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// The EXPOSURE of analysis::ClipDetector — fcore::ClipProbe, the `clips` text format and the wasm C ABI.
// The detector itself is tested by felitronics_clip_detector_tests and is not touched here; what is under
// test is the wiring P71 adds, which is where an instrument learns to lie.
//
// FOUR THINGS IT PROVES, IN ORDER OF WHAT WOULD HURT MOST IF IT WERE FALSE:
//
//  1. NO ROAD CAN CERTIFY A FILE IT DID NOT MEASURE. A report read before finish() is short by up to
//     decisionDelaySamples() of runs; a refused process() still lets the detector finish; and a refused
//     prepare() leaves the PREVIOUS file's answer behind isFinished() == true (measured on this tree: a good
//     run, finish(), then prepare(0.0, …) → isFinished() true, runCount() 1, samplePeak(0) 0.625). Each of
//     those three produces the same wrong output — `runs 0 / stored 0 / complete 1`, which reads as a
//     certificate of cleanliness. Every one of them is asserted refused, and the negative control is spelled
//     out: the mutation stand below removes each guard and the suite must go red.
//
//  2. LAW 8a ON THE SHIPPED ADAPTER, NOT ON THE DETECTOR. fcore::ClipProbe is the class both roads drive, and
//     it is the one that walks the planes and cuts the stream. So the slicing matrix runs THROUGH IT, and the
//     bug class it exists to exclude — `planar[c]` instead of `planar[c] + off` in the chunk loop — is one a
//     test that hands the detector correct plane pointers cannot see. The detector's own re-slicing is
//     already covered in modules/analysis; re-proving it here would prove nothing about this diff.
//
//  3. AN INVARIANCE TEST COMPARES CODE WITH ITSELF, SO AN ORACLE SITS BESIDE IT. Fourteen slicings agreeing
//     bit for bit is consistency, not truth: a deterministically wrong adapter agrees with itself everywhere.
//     Two independent witnesses are therefore pinned. (a) A bare analysis::ClipDetector, driven in one call
//     and read field by field in this file, never through readClips() — the wrapper must equal it. (b) A
//     LITERAL BIT PATTERN for a level: eleven samples at 40/64 and one at 39/64 give the mean 479/768, whose
//     double is 3fe3f55555555555, cross-checked here against the exact rational. That constant is what
//     catches the one mutant no invariance comparison can: rounding every level through float, which is
//     slicing-invariant, is identical on both tiers, and prints the same under any decimal format.
//
//  4. THE C ABI IS THE SAME NUMBERS AND REFUSES THE SAME THINGS. Compiled natively here (EMSCRIPTEN_KEEPALIVE
//     degrades to extern "C"), so its validation, its packing and its two capacities run under ctest, ASan
//     and UBSan — the arrangement felitronics_abi_tests uses, and for the same reason.
//
// MUTATION STAND — run on 2026-09-15 against an ISOLATED COPY of the sources (never the shared
// ClipDetector.h, which five branches are reading): **13 of 13 mutants turn this suite red.** The stand's
// output is in the P71 report; the mutants were
//   m1  a band closed on the call boundary          m6b the poison does not clear validity
//   m2  the pending queue drained per process()     m7  the shim feeds frames/2
//   m3  every run level rounded through float       m8  the report readable before finish()
//   m4  the chunk loop without its plane offset     m9  the peak copier fills the caller's capacity
//   m5  runsComplete() always true                  m10 the run copier returns doubles, not runs
//   m6  finish() ignores the refusal poison         m11 the three clocks no longer compared
//                                                   m12 the run copier reads `cap` as runs, not doubles
// m3 is the one that matters most: it is slicing-invariant, identical on both tiers and prints the same
// under any decimal format, so every slicing comparison here and the whole parity harness pass it — only the
// pinned bit pattern kills it. m4 and m11 are the two that a test driving the detector instead of the adapter
// could not see at all. m6b, m9 and m12 came from the crew's post-build rounds, and m12 was a live hazard
// rather than an invented mutation: the entry point took its capacity in RUNS while the other seven copiers
// of that ABI take ELEMENTS.

#include <felitronics_test.h>

#include "fcore_clips_format.h"

#include <algorithm>
#include <atomic>
#include <bit>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <random>
#include <string>
#include <vector>

// The allocation counter. Enforced only on libc++ (felitronics_test.h says why), which is this developer
// machine and the wasm row — so an allocation added to process() is caught HERE or not at all.
static std::atomic<long> g_allocs { 0 };
void* operator new      (std::size_t s) { g_allocs.fetch_add (1); return std::malloc (s ? s : 1); }
void* operator new[]    (std::size_t s) { g_allocs.fetch_add (1); return std::malloc (s ? s : 1); }
void  operator delete   (void* p) noexcept { std::free (p); }
void  operator delete[] (void* p) noexcept { std::free (p); }
void  operator delete   (void* p, std::size_t) noexcept { std::free (p); }
void  operator delete[] (void* p, std::size_t) noexcept { std::free (p); }

using namespace felitronics;
using CD = analysis::ClipDetector;
using Planes = std::vector<std::vector<float>>;

// The C ABI, declared here rather than included: fc_probe.cpp has no header, and writing the prototypes out
// is what makes this a test OF AN ABI instead of a test of a C++ class that happens to sit behind one. A
// signature that drifts fails to link, which is the diagnosis.
extern "C" {
    int           fc_probe_clips_run (const float* planar, std::uint32_t frames, std::uint32_t channels,
                                      double sampleRate, std::uint32_t maxRuns);
    double        fc_probe_clips_rate     (void);
    std::uint32_t fc_probe_clips_channels (void);
    std::uint32_t fc_probe_clips_samples  (void);
    std::uint32_t fc_probe_clips_max_runs (void);
    std::uint32_t fc_probe_clips_delay    (void);
    std::uint32_t fc_probe_clips_count    (void);
    std::uint32_t fc_probe_clips_stored   (void);
    int           fc_probe_clips_complete (void);
    std::uint32_t fc_probe_clips_stride   (void);
    std::uint32_t fc_probe_clips_peaks (double* out, std::uint32_t cap);
    std::uint32_t fc_probe_clips_runs  (double* out, std::uint32_t cap);
}

//================================================================================================ helpers

static std::uint64_t b64 (double d) noexcept { return std::bit_cast<std::uint64_t> (d); }

// BIT-FOR-BIT, FIELD BY FIELD. Not memcmp of the struct (ClipRun has padding, and padding is not a
// measurement) and not `!=` (which cannot tell -0.0 from +0.0 and calls two NaNs different). The pattern is
// WaveformShapeTests.cpp:163's; ClipDetectorTests.cpp:172 compares numerically and is not the model.
static bool sameRun (const analysis::ClipRun& a, const analysis::ClipRun& b)
{
    return a.start == b.start && a.length == b.length && b64 (a.level) == b64 (b.level)
        && a.channel == b.channel && a.sign == b.sign && a.evidence == b.evidence;
}

static bool sameReport (const fcore::ClipsReport& a, const fcore::ClipsReport& b)
{
    if (a.ok != b.ok || b64 (a.sampleRate) != b64 (b.sampleRate) || a.channels != b.channels
        || a.frames != b.frames || a.maxRuns != b.maxRuns || a.decisionDelay != b.decisionDelay
        || a.runCount != b.runCount || a.storedRunCount != b.storedRunCount || a.complete != b.complete
        || a.peak.size() != b.peak.size() || a.runs.size() != b.runs.size())
        return false;
    for (std::size_t c = 0; c < a.peak.size(); ++c) if (b64 (a.peak[c]) != b64 (b.peak[c])) return false;
    for (std::size_t i = 0; i < a.runs.size(); ++i) if (! sameRun (a.runs[i], b.runs[i])) return false;
    return true;
}

static std::vector<const float*> ptrsAt (const Planes& p, std::size_t off)
{
    std::vector<const float*> v (p.size());
    for (std::size_t c = 0; c < p.size(); ++c) v[c] = p[c].data() + off;
    return v;
}

// A slicing: the call lengths, repeated in order. `{0}` means one call for the whole stream; a negative
// length means "ragged, from this seed".
struct Slicing
{
    std::vector<long long> steps;
    int seed = 0;
    std::string name;
};

static std::vector<long long> sliceLengths (const Slicing& s, long long total)
{
    std::vector<long long> out;
    if (s.seed != 0)
    {
        std::mt19937 rng ((unsigned) s.seed);
        std::uniform_int_distribution<long long> U (0, std::max<long long> (1, total / 7));
        for (long long at = 0; at < total; )
        {
            const long long n = std::min (total - at, U (rng));
            out.push_back (n);                                            // a zero-length call is legal and must change nothing
            at += n;
        }
        return out;
    }
    if (s.steps.size() == 1 && s.steps[0] == 0) { out.push_back (total); return out; }
    for (long long at = 0, k = 0; at < total; ++k)
    {
        const long long n = std::min (total - at, s.steps[(std::size_t) k % s.steps.size()]);
        out.push_back (n);
        at += n;
    }
    return out;
}

// THE ADAPTER UNDER TEST. Everything the two roads share runs here: prepare, the chunked plane walk, finish,
// and the report reader that refuses an unfinished measurement.
static fcore::ClipsReport viaProbe (const Planes& p, double sr, std::int64_t maxRuns, const Slicing& sl,
                                    long long frames, bool* noAlloc = nullptr)
{
    fcore::ClipProbe probe;
    fcore::ClipsReport r;
    if (! probe.prepare (sr, (int) p.size(), maxRuns, frames)) return r;
    // EVERYTHING THE HARNESS ALLOCATES IS ALLOCATED BEFORE THE COUNTER STARTS — the call lengths and the
    // plane table both. The first draft of this function built them inside the loop and the no-alloc check
    // failed on every row: a counter that is also counting the fixture proves nothing about process().
    const std::vector<long long> lens = sliceLengths (sl, frames);
    const int nch = (int) p.size();
    const float* view[core::kMaxChannels] {};
    const long before = g_allocs.load();
    long long at = 0;
    for (long long n : lens)
    {
        for (int c = 0; c < nch; ++c) view[c] = p[(std::size_t) c].data() + at;
        if (! probe.process (view, nch, n)) return r;
        at += n;
    }
    if (! probe.finish()) return r;
    if (noAlloc != nullptr) *noAlloc = g_allocs.load() == before;
    readClips (probe, r);
    return r;
}

// THE ORACLE BESIDE IT: a bare analysis::ClipDetector, one call, read field by field here. It shares the
// report STRUCT with the adapter and not one line of the adapter's code, so "the wrapper equals this" is a
// claim about the wrapper rather than about a slicing.
static fcore::ClipsReport bareWhole (const Planes& p, double sr, std::int64_t maxRuns, long long frames)
{
    fcore::ClipsReport r;
    CD d;
    analysis::ClipDetectorParams prm;
    prm.maxRuns = (int) maxRuns;
    d.setParams (prm);
    if (! d.prepare (sr, 4096, (int) p.size())) return r;
    const auto v = ptrsAt (p, 0);
    if (frames > 0 && ! d.process (v.data(), (int) p.size(), (int) frames)) return r;
    d.finish();
    r.sampleRate = sr;
    r.channels = (int) p.size();
    r.frames = d.samplesProcessed();
    r.maxRuns = maxRuns;
    r.decisionDelay = d.decisionDelaySamples();
    r.runCount = d.runCount();
    r.storedRunCount = d.storedRunCount();
    r.complete = d.runsComplete();
    for (int c = 0; c < (int) p.size(); ++c) r.peak.push_back (d.samplePeak (c));
    for (std::int64_t i = 0; i < d.storedRunCount(); ++i) r.runs.push_back (d.run (i));
    r.ok = true;
    return r;
}

//================================================================================================ material

// THE PINNED WITNESS, on a 1/64 grid at 1 kHz (W = 20 samples, so every prefix of it is cheap to test).
// A linear approach of 10/64 per sample, a plateau of twelve, and a linear departure. With the plateau
// exactly flat the mean is 40/64 = 0.625 — a dyadic number, which survives a float round-trip unchanged and
// therefore witnesses nothing about the last bits. With ONE plateau sample at 39/64 the mean is 479/768,
// which is not dyadic; that is the constant this suite pins.
//
// The approach step is 10/64 on BOTH sides and load-bearing at that value. The chance bound asks that
// ((hi - lo + q)/sigma)^(L-1) <= 1e-6, with sigma the mean |step| over the flanks: at this slope sigma is
// about 8.2/64, the ratio 0.244, and 0.244^11 = 1.9e-7 — inside the bound. The sequence at
// ClipDetectorTests.cpp:642, whose approach reaches the plateau from 30/64 with a shorter ramp, gives sigma
// 6.9/64, ratio 0.288 and 0.288^11 = 1.1e-6 — OUTSIDE it, and the same 40→39 edit there finds NO run at all.
// Measured both ways before this constant was written down.
static Planes pinnedWitness (bool oneOdd)
{
    const double g = 1.0 / 64.0;
    std::vector<float> x;
    auto put = [&] (int code) { x.push_back ((float) ((double) code * g)); };
    for (int c : { 0, 1, 0, -1 }) put (c);                                // shows the stream a one-code step: q = 1/64
    for (int k = -2; k <= 2; ++k) put (k * 10);
    for (int i = 0; i < 12; ++i) put (i == 11 && oneOdd ? 39 : 40);       // the plateau
    for (int k = 2; k >= -2; --k) put (k * 10);
    for (int i = 0; i < 30; ++i) { put (0); put (-3); }                   // 20 ms of quiet, so the run is decided
    return { x };
}

// Non-stationary, seeded, CLAMPED — the material the slicing matrix runs on. A clamp is applied here, so
// "was a run found" is a comparison with a fact. Holes and a clipped tail are deliberate: the tail's runs
// exist only after finish(), and a hole splits a run so the second half is found by CEILING evidence.
static Planes clampedProgramme (long long n, double sr, int channels, unsigned seed, bool holes, bool tail)
{
    Planes p ((std::size_t) channels, std::vector<float> ((std::size_t) n, 0.0f));
    std::mt19937 rng (seed);
    std::uniform_int_distribution<int> N (-40, 40);
    const int ceilCode = 29196;
    for (int c = 0; c < channels; ++c)
    {
        const long long period = 96 + 16 * c;
        for (long long i = 0; i < n; ++i)
        {
            const long long k = (i + 7 * c) % (2 * period), t = k % period;
            const double sign = k < period ? 1.0 : -1.0;
            const double u = (double) (2 * t - period) / (double) period;
            double code = sign * 1.25 * (double) ceilCode * (1.0 - u * u) + (double) N (rng);
            if (! tail && i > n - (long long) (sr * 0.03)) code = 0.0;    // a quiet tail: nothing pending at finish()
            const double clamped = std::max (-(double) ceilCode, std::min ((double) ceilCode, std::round (code)));
            float v = (float) (clamped / 32768.0);
            if (holes && (i % 977) == 13 * (c + 1))
                v = (i % 3 == 0) ? std::numeric_limits<float>::quiet_NaN()
                                 : (i % 3 == 1) ? std::numeric_limits<float>::infinity()
                                                : -std::numeric_limits<float>::infinity();
            p[(std::size_t) c][(std::size_t) i] = v;
        }
    }
    return p;
}

//================================================================================================ the suite

static void theFormatIsBitwise()
{
    test::group ("the format is bitwise, so one flipped bit changes a character");
    const Planes w = pinnedWitness (true);
    const long long n = (long long) w[0].size();
    const fcore::ClipsReport r = viaProbe (w, 1000.0, 1 << 16, { { 0 }, 0, "whole" }, n);
    test::ok (r.ok && r.runs.size() == 1, "the pinned witness has exactly one run (" + std::to_string (r.runs.size()) + ")");
    if (r.runs.size() != 1) return;

    // The level is the mean of eleven 40/64 and one 39/64 = 479/768. Pinned as a literal AND derived from the
    // exact rational: two independent routes to the same double, so a typo in the constant cannot pass.
    const double want = 479.0 / 768.0;
    test::ok (b64 (want) == 0x3fe3f55555555555ull, "479/768 is 3fe3f55555555555");
    test::ok (b64 (r.runs[0].level) == 0x3fe3f55555555555ull,
              "the run's level is that bit pattern, not its float32 rounding (got "
              + std::to_string ((unsigned long long) b64 (r.runs[0].level)) + ")");
    test::ok (b64 ((double) (float) want) == 0x3fe3f55560000000ull,
              "and the float32 rounding of it IS a different pattern — so the check above can fail");
    test::ok (r.runs[0].start == 9 && r.runs[0].length == 12 && r.runs[0].sign == 1
              && r.runs[0].evidence == analysis::ClipEvidence::Ramp, "start 9, length 12, a positive Ramp");

    // The text, end to end. A level printed through %.6f would be 0.623698 for both patterns above.
    const std::string text = fcore::formatClips (r);
    test::ok (text.find ("run 9 12 3fe3f55555555555 0 1 1\n") != std::string::npos,
              "the run line carries the level as its bit pattern");
    test::ok (text.find ("# fcore clips v1 sr=408f400000000000 ch=1 frames=") == 0,
              "the header's sample rate is a bit pattern too (1000.0 = 408f400000000000)");
    test::ok (text.find ("delay=20\n") != std::string::npos, "and the decision delay it was measured with");

    // An exactly flat plateau has a dyadic mean, which is why the odd sample is there at all.
    const fcore::ClipsReport flat = viaProbe (pinnedWitness (false), 1000.0, 1 << 16, { { 0 }, 0, "whole" }, n);
    test::ok (flat.runs.size() == 1 && b64 (flat.runs[0].level) == b64 (0.625),
              "an exactly flat plateau's level is 0.625 exactly — dyadic, and no witness of rounding");
}

static void lawEightA()
{
    test::group ("law 8a: the report cannot see where the stream was cut — through the SHIPPED adapter");
    const std::vector<Slicing> slicings {
        { { 0 }, 0, "whole" },        { { 1 }, 0, "1" },            { { 2 }, 0, "2" },
        { { 3 }, 0, "3" },            { { 7 }, 0, "7 (prime)" },    { { 19 }, 0, "W-1" },
        { { 20 }, 0, "W" },           { { 21 }, 0, "W+1" },         { { 8191 }, 0, "kChunk-1" },
        { { 8192 }, 0, "kChunk" },    { { 8193 }, 0, "kChunk+1" },  { { 40000 }, 0, "> kChunk" },
        { { 1, 8192, 3, 20 }, 0, "mixed" },
        { {}, 12345, "ragged (seed 12345)" }, { {}, 999, "ragged (seed 999)" },
    };

    struct Case { long long n; double sr; int ch; unsigned seed; bool holes, tail; const char* what; };
    const std::vector<Case> cases {
        { 0,      1000.0, 1, 1, false, false, "T = 0 (an empty stream is a measurement)" },
        { 1,      1000.0, 1, 1, false, false, "T = 1" },
        { 19,     1000.0, 1, 1, false, false, "T = W-1" },
        { 20,     1000.0, 1, 1, false, false, "T = W" },
        { 21,     1000.0, 1, 1, false, false, "T = W+1" },
        { 231,    1000.0, 1, 1, false, true,  "the pinned witness's length" },
        { 24000,  8000.0, 2, 7, false, true,  "8 kHz stereo, clamped, clipped tail" },
        { 24000,  8000.0, 2, 7, true,  true,  "…with non-finite holes" },
        { 24000,  8000.0, 3, 9, true,  false, "three channels, quiet tail" },
        { 60000, 48000.0, 2, 3, true,  true,  "48 kHz stereo, past kChunk" },
    };

    int compared = 0;
    for (const Case& k : cases)
    {
        const Planes p = k.n == 231 ? pinnedWitness (true)
                                    : clampedProgramme (k.n, k.sr, k.ch, k.seed, k.holes, k.tail);
        const long long n = (long long) p[0].size();
        // The oracle first, and it is NOT the whole-in-one-call slicing of the adapter: a bare detector,
        // driven and read in this file.
        const fcore::ClipsReport oracle = bareWhole (p, k.sr, 1 << 16, n);
        bool noAlloc = true;
        const fcore::ClipsReport base = viaProbe (p, k.sr, 1 << 16, slicings[0], n, &noAlloc);
        test::ok (base.ok && oracle.ok, std::string ("both roads measured: ") + k.what);
        test::ok (sameReport (base, oracle),
                  std::string ("the adapter equals a bare detector read outside it: ") + k.what);
        test::okNoAlloc (noAlloc, std::string ("no allocation in process()/finish(): ") + k.what);
        for (std::size_t i = 1; i < slicings.size(); ++i)
        {
            bool na = true;
            const fcore::ClipsReport r = viaProbe (p, k.sr, 1 << 16, slicings[i], n, &na);
            test::ok (sameReport (base, r), std::string (k.what) + " — slicing " + slicings[i].name);
            test::okNoAlloc (na, std::string (k.what) + " — no alloc, slicing " + slicings[i].name);
            ++compared;
        }
        // maxBlock is not a parameter of this measurement — ClipDetector::prepare() names its second argument
        // `maxBlock: nothing is sized by it`. The line is here because law 8a asks for it explicitly, and the
        // cheapest way to keep a claim honest is to let it fail: the adapter fixes maxBlock at kChunk and the
        // oracle above passes 4096, so every row of this matrix already compares two different ones.
        test::ok (fcore::ClipProbe::kChunk != 4096, "the adapter's step and the oracle's differ, so that is tested");
    }
    std::printf ("    (%d slicing comparisons)\n", compared);
}

// A CHANNEL THAT DISAPPEARS. Law 11a: a call that does not carry a channel makes it a hole for those
// samples, which is the only way a stream whose width changes can be expressed. The presence timeline is part
// of the INPUT, not of the slicing — so every slicing here is additionally cut at the moment the channel
// goes away, and what is compared is the same measurement under different cuts of it.
static void aChannelThatDisappears()
{
    test::group ("a channel that disappears mid-stream is a hole, not an error (law 11a)");
    const Planes p = clampedProgramme (4800, 1000.0, 2, 17, false, true);
    const long long n = (long long) p[0].size(), gone = n / 2;

    auto drive = [&] (const std::vector<long long>& firstHalf, const std::vector<long long>& secondHalf)
    {
        fcore::ClipProbe probe;
        fcore::ClipsReport r;
        if (! probe.prepare (1000.0, 2, 1 << 16, n)) return r;
        const float* view[core::kMaxChannels] {};
        long long at = 0;
        for (int half = 0; half < 2; ++half)
            for (long long step : half == 0 ? firstHalf : secondHalf)
            {
                for (int c = 0; c < 2; ++c) view[c] = p[(std::size_t) c].data() + at;
                if (! probe.process (view, half == 0 ? 2 : 1, step)) return r;   // the second channel stops arriving
                at += step;
            }
        if (! probe.finish()) return r;
        readClips (probe, r);
        return r;
    };
    auto cut = [] (long long total, long long step)
    {
        std::vector<long long> v;
        for (long long a = 0; a < total; a += step) v.push_back (std::min (step, total - a));
        return v;
    };

    const fcore::ClipsReport base = drive (cut (gone, gone), cut (n - gone, n - gone));
    test::ok (base.ok && base.frames == n, "two calls, the second one channel narrower");
    test::ok (base.runCount > 0, "and it found runs (" + std::to_string (base.runCount) + ")");
    // Channel 1 stops at `gone`, so nothing of it can be reported past that coordinate.
    bool noneLate = true;
    for (const auto& u : base.runs) if (u.channel == 1 && u.start >= gone) noneLate = false;
    test::ok (noneLate, "no run of the vanished channel starts after it vanished");
    test::ok (base.peak.size() == 2 && base.peak[1] > 0.0, "its peak is the peak of the samples it did send");
    for (long long step : { (long long) 1, (long long) 3, (long long) 17, (long long) 20, (long long) 512, (long long) 8192 })
        test::ok (sameReport (base, drive (cut (gone, step), cut (n - gone, step))),
                  "the same report when each half is cut into " + std::to_string (step) + "-sample calls");

    // channels == 0 — every channel a hole — is legal, and then the plane table need not exist at all.
    {
        fcore::ClipProbe probe;
        test::run (probe.prepare (1000.0, 2, 16, 100));
        test::ok (probe.process (nullptr, 0, 100), "channels == 0 with a null plane table is accepted");
        test::run (probe.finish());
        test::ok (probe.valid() && probe.frames() == 100 && probe.runCount() == 0 && probe.peak (0) == 0.0,
                  "100 samples of pure hole: measured, no runs, no peak");
    }
    // A width the preparation does not have is still refused, and still poisons.
    {
        fcore::ClipProbe probe;
        const float* view[2] { p[0].data(), p[1].data() };
        test::run (probe.prepare (1000.0, 2, 16, n));
        test::ok (! probe.process (view, 3, 10), "a width wider than the preparation is refused");
        test::ok (! probe.finish(), "…and poisons the measurement, because the audio of that call is missing");
    }
}

// EVERY RUN'S LEVEL, RECOMPUTED FROM THE SAMPLES BY A DIFFERENT ROUTE. The pinned witness gives one bit
// pattern from a hand calculation; this gives one for all of them, from the plane data, with no streaming
// state involved. It is EXACT, not approximate: the fixture's samples are k/32768 with |k| <= 32768, so each
// is exact in a double, and a sum of at most a few thousand of them has a numerator below 2^28 — so the sum
// is exact and the single division that follows is the only rounding, which is precisely what the detector
// does (a compensated accumulation of exact values has a zero correction). A level that disagrees here means
// the run's coordinates are wrong, not that the arithmetic is.
static void everyLevelAgainstTheSamples()
{
    test::group ("every reported run's level, recomputed from the samples outside the detector");
    int checked = 0;
    for (int ch : { 1, 2, 3 })
    {
        const Planes p = clampedProgramme (24000, 8000.0, ch, 31u + (unsigned) ch, false, true);
        const long long n = (long long) p[0].size();
        const fcore::ClipsReport r = viaProbe (p, 8000.0, 1 << 16, { { 0 }, 0, "whole" }, n);
        test::ok (r.ok && r.runCount > 20, "a report with runs to check (" + std::to_string (r.runCount) + ")");
        bool levels = true, inside = true, flat = true;
        for (const auto& u : r.runs)
        {
            if (u.start < 0 || u.length < 2 || u.start + u.length > n || u.channel < 0 || u.channel >= ch
                || (u.sign != 1 && u.sign != -1)) { inside = false; continue; }
            const std::vector<float>& y = p[(std::size_t) u.channel];
            double sum = 0.0, lo = y[(std::size_t) u.start], hi = lo;
            for (std::int64_t i = u.start; i < u.start + u.length; ++i)
            {
                const double v = (double) y[(std::size_t) i];
                sum += v;
                lo = std::min (lo, v);
                hi = std::max (hi, v);
            }
            if (b64 (sum / (double) u.length) != b64 (u.level)) levels = false;
            // The run must actually be FLAT: its range cannot exceed tau = 2q, and q here is at most one
            // 16-bit code because the generator shows the stream a one-code step in its first samples.
            if (! (hi - lo <= 2.0 / 32768.0)) flat = false;
            ++checked;
        }
        test::ok (inside, "every run lies inside the stream, on a channel it has, with a polarity");
        test::ok (levels, "every level is the samples' mean, bit for bit, at " + std::to_string (ch) + " channels");
        test::ok (flat, "and every run is flat to within 2q");
    }
    std::printf ("    (%d runs recomputed from the samples)\n", checked);
}

// The rate goes into the report as a bit pattern, and the decision delay is derived from it. Both are places
// where a narrower type or a rounded constant would hide for a long time.
static void ratesAndTheirWindows()
{
    test::group ("the sample rate and the window it implies");
    const Planes p = clampedProgramme (4000, 48000.0, 1, 3, false, true);
    const long long n = (long long) p[0].size();

    // A FRACTIONAL RATE. 48000.1 is not representable in float32 without loss: the double is
    // 40e7700333333333 and its float32 rounding is 40e7700340000000. The report must carry the double.
    const fcore::ClipsReport frac = viaProbe (p, 48000.1, 1 << 16, { { 0 }, 0, "whole" }, n);
    test::ok (frac.ok && b64 (frac.sampleRate) == 0x40e7700333333333ull,
              "48000.1 is reported as the double it is, not as its float32 rounding");
    test::ok (b64 ((double) (float) 48000.1) == 0x40e7700340000000ull, "…and those two patterns do differ");
    test::ok (frac.decisionDelay == 960, "its window is floor(48000.1 * 20 / 1000) = 960");

    // THE WINDOW'S OWN STEP. floor(sr/50) changes at 1050, so 1049.9 and 1050 must not report the same delay
    // — a rate rounded anywhere on the way in would make them equal.
    test::ok (viaProbe (p, 1049.9, 16, { { 0 }, 0, "whole" }, n).decisionDelay == 20, "1049.9 -> a 20-sample window");
    test::ok (viaProbe (p, 1050.0, 16, { { 0 }, 0, "whole" }, n).decisionDelay == 21, "1050.0 -> 21");

    // THE ENDS OF THE RANGE, and one step past each.
    test::ok (viaProbe (p, 1000.0, 16, { { 0 }, 0, "whole" }, n).ok, "kMinSampleRate is accepted");
    test::ok (viaProbe (p, 768000.0, 16, { { 0 }, 0, "whole" }, n).decisionDelay == 15360,
              "kMaxSampleRate is accepted, with a 15360-sample window");
    fcore::ClipProbe probe;
    test::ok (! probe.prepare (std::nextafter (1000.0, 0.0), 1, 16, n), "one ulp below the minimum is refused");
    test::ok (! probe.prepare (std::nextafter (768000.0, 2e6), 1, 16, n), "one ulp above the maximum is refused");
    test::ok (! probe.prepare (std::numeric_limits<double>::quiet_NaN(), 1, 16, n), "NaN is refused");
    test::ok (! probe.prepare (std::numeric_limits<double>::infinity(), 1, 16, n), "infinity is refused");
}

// Streams that are nothing but one pathology each. Sparse holes in ordinary material do not reach these.
static void pathologicalStreams()
{
    test::group ("streams that are nothing but one pathology");
    const long long n = 4000;
    auto oneChannel = [&] (float v) { return Planes { std::vector<float> ((std::size_t) n, v) }; };

    // SILENCE WITH NO CAPACITY AT ALL. `complete` is `count <= capacity`, so a file with no runs is COMPLETE
    // even at capacity zero — the opposite of what "capacity 0 means incomplete" would suggest.
    {
        const fcore::ClipsReport r = viaProbe (oneChannel (0.0f), 8000.0, 0, { { 0 }, 0, "whole" }, n);
        test::ok (r.ok && r.runCount == 0 && r.storedRunCount == 0 && r.complete,
                  "pure silence at capacity 0 is COMPLETE: there was nothing to store");
        test::ok (b64 (r.peak[0]) == b64 (0.0), "and its peak is +0.0, bit for bit");
    }
    // ALL -0.0. The core treats signed zeros alike, and a peak is an absolute value, so +0.0 must come out.
    {
        const fcore::ClipsReport r = viaProbe (oneChannel (-0.0f), 8000.0, 16, { { 0 }, 0, "whole" }, n);
        test::ok (r.ok && r.runCount == 0, "a stream of -0.0 has no runs");
        test::ok (b64 (r.peak[0]) == b64 (0.0) && b64 (r.peak[0]) != b64 (-0.0),
                  "its peak is +0.0 and the comparison can tell the two zeros apart");
    }
    // SUBNORMALS. The smallest positive float32 is 2^-149, whose double is 36a0000000000000. If any tier
    // flushed subnormals to zero this peak would read 0 — which is why it is checked here and not only
    // asserted about the arithmetic.
    {
        Planes p = oneChannel (0.0f);
        p[0][10] = std::numeric_limits<float>::denorm_min();
        p[0][20] = -std::numeric_limits<float>::denorm_min();
        const fcore::ClipsReport r = viaProbe (p, 8000.0, 16, { { 0 }, 0, "whole" }, n);
        test::ok (r.ok && b64 (r.peak[0]) == 0x36a0000000000000ull,
                  "the smallest float32 subnormal survives into the peak as 2^-149");
    }
    // A CHANNEL THAT IS ENTIRELY NON-FINITE, beside one that is clipped. The dead channel must report
    // nothing at all and must not disturb its neighbour.
    {
        const Planes live = clampedProgramme (n, 8000.0, 1, 77, false, true);
        Planes p { std::vector<float> ((std::size_t) n, std::numeric_limits<float>::quiet_NaN()), live[0] };
        for (long long i = 0; i < n; i += 3) p[0][(std::size_t) i] = std::numeric_limits<float>::infinity();
        for (long long i = 1; i < n; i += 3) p[0][(std::size_t) i] = -std::numeric_limits<float>::infinity();
        const fcore::ClipsReport r = viaProbe (p, 8000.0, 1 << 16, { { 0 }, 0, "whole" }, n);
        test::ok (r.ok && b64 (r.peak[0]) == b64 (0.0), "an all-non-finite channel has a +0.0 peak");
        bool none = true;
        for (const auto& u : r.runs) if (u.channel == 0) none = false;
        test::ok (none, "…and no runs");
        const fcore::ClipsReport alone = viaProbe (live, 8000.0, 1 << 16, { { 0 }, 0, "whole" }, n);
        test::ok (alone.ok && r.runCount == alone.runCount, "and the live channel reports what it does alone");
        test::ok (sameReport (r, viaProbe (p, 8000.0, 1 << 16, { { 7 }, 0, "7" }, n)), "under any slicing");
    }
}

// A KNOWN RUN ON THE READER'S OWN BOUNDARY. ClipProbe cuts at kChunk internally and the CLI's file reader
// cuts at 8192 frames, so a run that starts at, ends at, or straddles that coordinate is the sharpest case
// for both. The witness's run is at start 9 of its own block, so the prefix does the placing.
static void runsOnTheChunkBoundary()
{
    test::group ("a known run placed on the adapter's own chunk boundary");
    const Planes w = pinnedWitness (true);
    const long long K = fcore::ClipProbe::kChunk;
    struct Place { long long pad; const char* what; };
    for (const Place& pl : std::vector<Place> { { K - 9, "starts exactly at kChunk" },
                                                { K - 21, "ends exactly at kChunk" },
                                                { K - 10, "straddles kChunk from one before" },
                                                { K - 15, "straddles kChunk in the middle" },
                                                { 2 * K - 9, "starts at 2*kChunk" } })
    {
        Planes p { std::vector<float> ((std::size_t) pl.pad, 0.0f) };
        p[0].insert (p[0].end(), w[0].begin(), w[0].end());
        const long long n = (long long) p[0].size();
        const fcore::ClipsReport base = viaProbe (p, 1000.0, 1 << 16, { { 0 }, 0, "whole" }, n);
        test::ok (base.ok && base.runs.size() == 1 && base.runs[0].start == pl.pad + 9
                  && base.runs[0].length == 12 && b64 (base.runs[0].level) == 0x3fe3f55555555555ull,
                  std::string ("the run is found where it was put, level unchanged — it ") + pl.what);
        for (const Slicing& sl : std::vector<Slicing> { { { 1 } }, { { 8191 } }, { { 8192 } }, { { 8193 } },
                                                        { { 20 } }, { {}, 4242, "ragged" } })
            test::ok (sameReport (base, viaProbe (p, 1000.0, 1 << 16, sl, n)),
                      std::string ("…and under every slicing when it ") + pl.what);
    }
}

static void theTraceOfDecisions()
{
    test::group ("the intermediate trace: every decision lands at the same sample under every slicing");
    // Comparing final reports is not enough — LAW8-KWEIGHTING.md:78 has this project's own example, where the
    // final LUFS matched bit for bit while 5 of 97 intermediate block energies had moved. ClipDetector has no
    // frames to hook, but it has something better: the run list IS the event trace, and samplesProcessed() is
    // its clock. So the reference is built by feeding ONE SAMPLE AT A TIME and snapshotting after each — the
    // exact coordinate at which every run becomes visible — and each slicing is then checked against that
    // reference AT EVERY ONE OF ITS OWN CALL BOUNDARIES. A run appended one call early or one call late is
    // caught there even when finish() would later hide it.
    const Planes p = clampedProgramme (2400, 1000.0, 2, 42, true, true);
    const long long n = (long long) p[0].size();

    std::vector<std::vector<analysis::ClipRun>> trace ((std::size_t) n + 1);
    std::vector<std::int64_t> counts ((std::size_t) n + 1, 0);
    {
        CD d;
        test::run (d.prepare (1000.0, 8192, 2));
        for (long long i = 0; i < n; ++i)
        {
            const auto v = ptrsAt (p, (std::size_t) i);
            test::run (d.process (v.data(), 2, 1));
            counts[(std::size_t) i + 1] = d.runCount();
            for (std::int64_t r = 0; r < d.storedRunCount(); ++r) trace[(std::size_t) i + 1].push_back (d.run (r));
        }
    }
    test::ok (counts[(std::size_t) n] > 0, "the trace saw runs decided while streaming ("
              + std::to_string (counts[(std::size_t) n]) + ")");
    // A trace that never moves would make every comparison below vacuous.
    int transitions = 0;
    for (long long i = 1; i <= n; ++i) if (counts[(std::size_t) i] != counts[(std::size_t) i - 1]) ++transitions;
    test::ok (transitions > 4, "and it moves at more than a handful of coordinates (" + std::to_string (transitions) + ")");

    for (const Slicing& sl : std::vector<Slicing> { { { 2 }, 0, "2" }, { { 3 }, 0, "3" }, { { 7 }, 0, "7" },
                                                    { { 19 }, 0, "19" }, { { 20 }, 0, "20" }, { { 21 }, 0, "21" },
                                                    { { 256 }, 0, "256" }, { {}, 777, "ragged" } })
    {
        CD d;
        test::run (d.prepare (1000.0, 8192, 2));
        long long at = 0;
        bool same = true;
        for (long long step : sliceLengths (sl, n))
        {
            const auto v = ptrsAt (p, (std::size_t) at);
            test::run (d.process (v.data(), 2, (int) step));
            at += step;
            if (d.samplesProcessed() != at) { same = false; break; }
            if (d.runCount() != counts[(std::size_t) at]) { same = false; break; }
            if ((std::size_t) d.storedRunCount() != trace[(std::size_t) at].size()) { same = false; break; }
            for (std::int64_t r = 0; r < d.storedRunCount(); ++r)
                if (! sameRun (d.run (r), trace[(std::size_t) at][(std::size_t) r])) { same = false; break; }
            if (! same) break;
        }
        test::ok (same, "every call boundary agrees with the per-sample trace — slicing " + sl.name);
    }
}

static void noCertificateWithoutAMeasurement()
{
    test::group ("no road can certify a file it did not measure");
    const Planes p = clampedProgramme (2400, 1000.0, 2, 5, false, true);
    const long long n = (long long) p[0].size();
    const auto v = ptrsAt (p, 0);

    // (a) BEFORE finish(): the run list is short by up to decisionDelaySamples() of runs, so there is no
    //     report at all — not a short one.
    {
        fcore::ClipProbe probe;
        test::run (probe.prepare (1000.0, 2, 1 << 16, n));
        test::run (probe.process (v.data(), 2, n));
        test::ok (! probe.valid(), "an unfinished measurement is not valid");
        test::ok (probe.runCount() == 0 && probe.storedRunCount() == 0 && ! probe.complete()
                  && probe.frames() == 0 && probe.peak (0) == 0.0 && probe.decisionDelay() == 0,
                  "and every accessor answers zero rather than a partial truth");
        // A report value that already holds a good answer must be CLEARED by the refusal, or the refusal is
        // worse than useless: the caller reads the previous file.
        fcore::ClipsReport reused = viaProbe (p, 1000.0, 1 << 16, { { 0 }, 0, "whole" }, n);
        test::ok (reused.ok && reused.runCount > 0, "a good report to reuse");
        test::ok (! readClips (probe, reused), "readClips refuses an unfinished probe");
        test::ok (! reused.ok && reused.runCount == 0 && reused.runs.empty(),
                  "and clears the report it was handed");
        test::ok (fcore::formatClips (reused).empty(), "a report that is not ok formats as nothing at all");
        test::run (probe.finish());
        test::ok (probe.valid() && probe.runCount() > 0, "after finish() there is a report, and it has runs");
        test::ok (probe.finish(), "finish() is idempotent — a second call answers the same verdict (law 12)");
        test::ok (probe.valid() && probe.runCount() > 0, "and does not throw the frozen report away");
    }

    // (b) A REFUSED process() POISONS THE MEASUREMENT. Without this, the detector still finishes and the
    //     report reads `runs 0 / stored 0 / complete 1` — a certificate of cleanliness for audio nothing
    //     looked at. Asserted here for the bare detector too, so the reason is on the record.
    {
        CD bare;
        test::run (bare.prepare (1000.0, 8192, 2));
        test::ok (! bare.process (v.data(), 9, 4), "the detector refuses a width it was not prepared for");
        bare.finish();
        test::ok (bare.isFinished() && bare.runCount() == 0 && bare.runsComplete() && bare.samplePeak (0) == 0.0,
                  "…and then finishes anyway, reporting a complete, empty, clean file — this is the trap");

        fcore::ClipProbe probe;
        test::run (probe.prepare (1000.0, 2, 1 << 16, n));
        test::ok (! probe.process (v.data(), 9, 4), "the adapter refuses it too");
        test::ok (! probe.finish(), "and then refuses to finish: a poisoned measurement has no result");
        test::ok (! probe.valid() && probe.runCount() == 0, "so nothing can be read out of it");
        test::ok (! probe.process (v.data(), 2, n), "a legal call after a refusal is refused as well");
    }

    // (b2) THE POISON'S EDGES. A refusal that carried samples must drop the report AT ONCE — not at the next
    //      finish(), because an accessor read in between would serve a report of less audio than was fed. A
    //      refusal that carried NOTHING must destroy nothing, or a drain loop that ends with a zero-length
    //      call, or a caller that probes a width with n == 0, loses a perfectly good measurement.
    {
        fcore::ClipProbe probe;
        test::run (probe.prepare (1000.0, 2, 1 << 16, n));
        test::run (probe.process (v.data(), 2, n));
        test::run (probe.finish());
        const std::int64_t was = probe.runCount();
        test::ok (was > 0, "a finished measurement");
        test::ok (! probe.process (v.data(), 2, 0), "a zero-length call after finish() is refused");
        test::ok (probe.valid() && probe.runCount() == was, "…and does NOT destroy the report");
        test::ok (! probe.process (v.data(), 9, 0), "a zero-length call of a bad width is refused");
        test::ok (probe.valid() && probe.runCount() == was, "…and does not destroy it either");
        test::ok (probe.finish() && probe.valid(), "finish() is still idempotent after those");
        test::ok (! probe.process (v.data(), 2, 4), "a call CARRYING samples after finish() is refused");
        test::ok (! probe.valid() && probe.runCount() == 0,
                  "…and drops the report immediately: those samples are not in it");
        test::ok (! probe.finish(), "and a later finish() cannot resurrect it");
    }
    // (b3) THE POLICY, STATED: a refusal poisons even when nothing was lost. The whole stream is fed, then one
    //      malformed call that consumes nothing, then finish() — and the report is refused although the three
    //      clocks agree. Without this assertion the flag in finish() would be indistinguishable from dead code.
    {
        fcore::ClipProbe probe;
        test::run (probe.prepare (1000.0, 2, 1 << 16, n));
        test::run (probe.process (v.data(), 2, n));
        test::ok (! probe.process (v.data(), 9, 4), "the whole stream, then one malformed call");
        test::ok (! probe.finish() && ! probe.valid(),
                  "…and there is no report, though every sample of the stream was measured");
    }

    // A CALL LONGER THAN THE STREAM IS REFUSED BEFORE THE FIRST SAMPLE MOVES. This is the check that needs
    // the declared length: without it `n` is simply believed and the adapter reads past the caller's buffer
    // (the first draft of this test did exactly that and died with SIGBUS, which is the honest demonstration
    // that no length-free adapter can refuse it).
    {
        fcore::ClipProbe probe;
        test::run (probe.prepare (1000.0, 2, 16, n));
        test::ok (! probe.process (v.data(), 2, n + 1), "one frame more than the stream holds is refused");
        test::ok (! probe.finish(), "…and poisons: it carried samples that were not consumed");
    }
    {
        fcore::ClipProbe probe;
        test::run (probe.prepare (1000.0, 2, 16, n));
        test::ok (! probe.process (v.data(), 2, std::numeric_limits<long long>::max()),
                  "and so is a length that would overflow the sample clock");
        test::ok (! probe.finish(), "…which is the same refusal, reached before any pointer arithmetic");
    }
    // A stream that stops short has no report either: the three clocks must agree.
    {
        fcore::ClipProbe probe;
        test::run (probe.prepare (1000.0, 2, 16, n));
        test::run (probe.process (v.data(), 2, n - 1));
        test::ok (! probe.finish(), "a file that delivered one frame less than it was sized for is refused");
        test::ok (! probe.valid(), "…and publishes nothing");
    }

    // (c) A REFUSED prepare() LEAVES THE PREVIOUS ANSWER BEHIND isFinished(). Measured, not assumed.
    {
        CD bare;
        test::run (bare.prepare (1000.0, 8192, 2));
        test::run (bare.process (v.data(), 2, n));
        bare.finish();
        const std::int64_t was = bare.runCount();
        test::ok (was > 0, "a good run on the bare detector");
        test::ok (! bare.prepare (0.0, 8192, 2), "prepare() refuses a rate below kMinSampleRate");
        test::ok (bare.isFinished() && bare.runCount() == was,
                  "and the refusal leaves isFinished() true with the old count — why the adapter has its own flag");

        fcore::ClipProbe probe;
        test::run (probe.prepare (1000.0, 2, 1 << 16, n));
        test::run (probe.process (v.data(), 2, n));
        test::run (probe.finish());
        test::ok (probe.valid() && probe.runCount() == was, "the adapter agrees on the count");
        test::ok (! probe.prepare (0.0, 2, 1 << 16, n), "a refused prepare on the adapter");
        test::ok (! probe.valid() && probe.runCount() == 0 && probe.frames() == 0 && probe.maxRuns() == 0,
                  "clears the result instead of serving the previous file");
        fcore::ClipsReport r;
        test::ok (! readClips (probe, r) && ! r.ok, "and the report path refuses");
    }

    // (d) prepare() is the re-anchor. A second measurement on the same object must be the first one's equal.
    {
        fcore::ClipProbe probe;
        for (int pass = 0; pass < 2; ++pass)
        {
            test::run (probe.prepare (1000.0, 2, 1 << 16, n));
            test::run (probe.process (v.data(), 2, n));
            test::run (probe.finish());
        }
        fcore::ClipsReport again;
        test::ok (readClips (probe, again), "a second measurement on the same adapter");
        test::ok (sameReport (again, bareWhole (p, 1000.0, 1 << 16, n)), "…is bit-identical to the first");
    }
}

static void capacityIsData()
{
    test::group ("capacity exhaustion is data: the count keeps counting and the list says it is short");
    const Planes p = clampedProgramme (24000, 8000.0, 2, 11, false, true);
    const long long n = (long long) p[0].size();
    const fcore::ClipsReport full = viaProbe (p, 8000.0, 1 << 16, { { 0 }, 0, "whole" }, n);
    test::ok (full.ok && full.complete && full.runCount == full.storedRunCount && full.runCount > 8,
              "the whole list first (" + std::to_string (full.runCount) + " runs)");

    for (std::int64_t cap : { (std::int64_t) 0, (std::int64_t) 1, (std::int64_t) 3,
                              full.runCount - 1, full.runCount, full.runCount + 1 })
    {
        const fcore::ClipsReport r = viaProbe (p, 8000.0, cap, { { 0 }, 0, "whole" }, n);
        const bool whole = cap >= full.runCount;
        test::ok (r.ok && r.runCount == full.runCount, "the COUNT is the same at capacity " + std::to_string (cap));
        test::ok (r.storedRunCount == std::min (cap, full.runCount), "stored is min(capacity, count)");
        test::ok (r.complete == whole, std::string ("complete is ") + (whole ? "true" : "false")
                  + " at capacity " + std::to_string (cap));
        test::ok ((std::int64_t) r.runs.size() == r.storedRunCount, "and the list is exactly that long");
        for (std::size_t i = 0; i < r.runs.size(); ++i)
            test::ok (sameRun (r.runs[i], full.runs[i]), "the stored runs are the PREFIX, unchanged, at capacity "
                      + std::to_string (cap));
        // The overflow happens at the same canonical run whatever the slicing: it is a property of the file.
        test::ok (sameReport (r, viaProbe (p, 8000.0, cap, { { 7 }, 0, "7" }, n)),
                  "and the same under a different slicing, at capacity " + std::to_string (cap));
    }
    // THE CAPACITY BOUNDARY WHERE EVERY CHANNEL DECIDES AT ONCE. The pinned witness copied into all sixteen
    // channels with identical timing gives sixteen runs at the same coordinate, so the capacity cuts across
    // the channel order rather than across time — and the prefix that survives must be channels 0..k, which
    // is finish()'s documented (end, channel) order.
    {
        const Planes one = pinnedWitness (true);
        Planes wide ((std::size_t) core::kMaxChannels, one[0]);
        const long long wn = (long long) one[0].size();
        const fcore::ClipsReport all = viaProbe (wide, 1000.0, 1 << 16, { { 0 }, 0, "whole" }, wn);
        test::ok (all.ok && all.runCount == core::kMaxChannels && all.complete,
                  "sixteen channels, sixteen simultaneous runs (" + std::to_string (all.runCount) + ")");
        bool ordered = true;
        for (std::size_t i = 0; i < all.runs.size(); ++i)
            if (all.runs[i].channel != (int) i || all.runs[i].start != 9
                || b64 (all.runs[i].level) != 0x3fe3f55555555555ull) ordered = false;
        test::ok (ordered, "in channel order, each at start 9 with the pinned level");
        const fcore::ClipsReport cut = viaProbe (wide, 1000.0, core::kMaxChannels - 1, { { 0 }, 0, "whole" }, wn);
        test::ok (cut.runCount == core::kMaxChannels && cut.storedRunCount == core::kMaxChannels - 1 && ! cut.complete,
                  "one short of sixteen: all counted, fifteen stored, incomplete");
        bool prefix = true;
        for (std::size_t i = 0; i < cut.runs.size(); ++i) if (! sameRun (cut.runs[i], all.runs[i])) prefix = false;
        test::ok (prefix, "and the fifteen are the first fifteen channels, unchanged");
        test::ok (sameReport (cut, viaProbe (wide, 1000.0, core::kMaxChannels - 1, { { 3 }, 0, "3" }, wn)),
                  "the same run is dropped under a different slicing");
    }
    test::ok (! fcore::ClipProbe {}.prepare (8000.0, 2, fcore::ClipProbe::kMaxRuns + 1, 1000),
              "a capacity past the ABI bound is refused (it would allocate hundreds of MB in a -fno-exceptions module)");
    test::ok (fcore::ClipProbe {}.prepare (8000.0, 2, fcore::ClipProbe::kMaxRuns, 1000),
              "the bound itself is accepted");
    test::ok (! fcore::ClipProbe {}.prepare (8000.0, 2, -1, 1000), "a negative capacity is refused");
}

static void theAbi()
{
    test::group ("the C ABI: the same numbers, the same refusals, and two capacities kept apart");
    const int ch = 2;
    const long long n = 12000;
    const Planes p = clampedProgramme (n, 8000.0, ch, 23, true, true);
    std::vector<float> planar ((std::size_t) n * (std::size_t) ch);
    for (int c = 0; c < ch; ++c)
        for (long long i = 0; i < n; ++i) planar[(std::size_t) (c * n + i)] = p[(std::size_t) c][(std::size_t) i];

    test::ok (fc_probe_clips_count() == 0 && fc_probe_clips_stored() == 0 && fc_probe_clips_complete() == 0
              && fc_probe_clips_samples() == 0 && fc_probe_clips_delay() == 0 && fc_probe_clips_rate() == 0.0,
              "before any run every getter is zero, complete included");
    test::ok (fc_probe_clips_stride() == 6, "the stride is 6 doubles per run");

    test::ok (fc_probe_clips_run (planar.data(), (std::uint32_t) n, (std::uint32_t) ch, 8000.0, 1 << 16) == 1,
              "a run through the ABI");
    const fcore::ClipsReport want = bareWhole (p, 8000.0, 1 << 16, n);
    test::ok (want.ok && want.runCount > 8, "the oracle has runs to compare (" + std::to_string (want.runCount) + ")");
    test::ok ((std::int64_t) fc_probe_clips_count() == want.runCount, "the ABI's count matches the oracle");
    test::ok ((std::int64_t) fc_probe_clips_stored() == want.storedRunCount, "…and its stored count");
    test::ok (fc_probe_clips_complete() == (want.complete ? 1 : 0), "…and its completeness");
    test::ok ((std::int64_t) fc_probe_clips_samples() == want.frames, "…and the samples it consumed");
    test::ok ((int) fc_probe_clips_channels() == want.channels, "…and the width");
    test::ok (b64 (fc_probe_clips_rate()) == b64 (8000.0), "…and the rate, bit for bit");
    test::ok ((std::int64_t) fc_probe_clips_delay() == want.decisionDelay, "…and the decision delay");
    test::ok (fc_probe_clips_max_runs() == (1u << 16), "…and the capacity it was prepared with");

    const std::uint32_t stored = fc_probe_clips_stored();
    {
        // AN OVERSIZED CAPACITY MUST NOT BE FILLED. The copier's source length is the CHANNEL COUNT and its
        // limit is the caller's capacity; swapping the two — returning `cap` and filling the surplus with
        // zeros, which is what an out-of-range peak() politely answers — is a real defect that every call
        // passing `cap == channels` would pass. So this one passes more, with canaries.
        std::vector<double> peaks ((std::size_t) ch + 2, -1.0);
        test::ok (fc_probe_clips_peaks (peaks.data(), (std::uint32_t) ch) == (std::uint32_t) ch, "the peaks copy whole");
        for (int c = 0; c < ch; ++c) test::ok (b64 (peaks[(std::size_t) c]) == b64 (want.peak[(std::size_t) c]),
                                               "peak " + std::to_string (c) + " bit for bit");
        std::fill (peaks.begin(), peaks.end(), -1.0);
        test::ok (fc_probe_clips_peaks (peaks.data(), (std::uint32_t) ch + 2u) == (std::uint32_t) ch,
                  "a capacity larger than the width returns the WIDTH, not the capacity");
        test::ok (b64 (peaks[(std::size_t) ch]) == b64 (-1.0) && b64 (peaks[(std::size_t) ch + 1]) == b64 (-1.0),
                  "…and writes nothing past the last channel");
        std::fill (peaks.begin(), peaks.end(), -1.0);           // the call above wrote two; refill before asking for one
        test::ok (fc_probe_clips_peaks (peaks.data(), 1u) == 1u, "a capacity of one returns one");
        test::ok (b64 (peaks[1]) == b64 (-1.0), "…and writes exactly one");
        test::ok (fc_probe_clips_peaks (peaks.data(), 0u) == 0u, "a capacity of zero returns zero");
        test::ok (fc_probe_clips_peaks (nullptr, (std::uint32_t) ch) == 0u, "a null peak buffer is refused");
    }
    {
        // A GUARD WORD PAST THE END. The claim is not only "it wrote the right numbers" but "it wrote exactly
        // as many as it said": an off-by-one here is a heap write in a page, and ASan catches it in this row.
        const std::uint32_t stride = fc_probe_clips_stride();
        std::vector<double> buf ((std::size_t) stored * stride + 1, -7.5);
        test::ok (fc_probe_clips_runs (buf.data(), stored * stride) == stored,
                  "the runs copy whole: the capacity is in DOUBLES, the answer is in RUNS");
        test::ok (b64 (buf.back()) == b64 (-7.5), "and not one double past them");
        bool all = true;
        for (std::uint32_t i = 0; i < stored; ++i)
        {
            const analysis::ClipRun& w = want.runs[i];
            const double* r = buf.data() + (std::size_t) i * stride;
            all = all && r[0] == (double) w.start && r[1] == (double) w.length && b64 (r[2]) == b64 (w.level)
                      && r[3] == (double) w.channel && r[4] == (double) w.sign
                      && r[5] == (double) (int) w.evidence;
        }
        test::ok (all, "every field of every run, bit for bit, in the detector's order");
    }
    {
        // The CALLER's capacity, which is not the file's. A short buffer is a short copy and nothing else:
        // `complete` must still describe the detector's list.
        const std::uint32_t stride = fc_probe_clips_stride();
        std::vector<double> buf ((std::size_t) stored * stride + 4, -7.5);
        test::ok (fc_probe_clips_runs (buf.data(), 3 * stride) == 3, "18 doubles of capacity return 3 RUNS");
        test::ok (b64 (buf[3 * stride]) == b64 (-7.5), "and write exactly 3 * stride doubles");
        test::ok (fc_probe_clips_complete() == 1, "a short copy does not make the FILE's list incomplete");
        // THE UNIT MISTAKE MUST BE SAFE. A caller who follows this ABI's other seven copiers and passes the
        // ELEMENT count it thinks it needs — `stored`, meaning runs — gets floor(stored/6) runs and nothing
        // written past its buffer. If `cap` were in runs this same call would overwrite six times the span.
        std::fill (buf.begin(), buf.end(), -7.5);
        test::ok (fc_probe_clips_runs (buf.data(), stored) == stored / stride,
                  "a capacity of `stored` DOUBLES yields floor(stored/6) runs, not stored runs");
        test::ok (b64 (buf[(std::size_t) (stored / stride) * stride]) == b64 (-7.5),
                  "…and writes nothing past the whole runs that fit");
        std::fill (buf.begin(), buf.end(), -7.5);              // the call above wrote most of the buffer
        test::ok (fc_probe_clips_runs (buf.data(), stride - 1) == 0,
                  "a capacity below one whole run writes nothing: a partial run is not a run");
        test::ok (b64 (buf[0]) == b64 (-7.5), "…literally nothing");
        test::ok (fc_probe_clips_runs (buf.data(), 0) == 0, "a capacity of zero writes nothing and returns 0");
        test::ok (fc_probe_clips_runs (nullptr, stored * stride) == 0, "a null buffer is refused");
        test::ok (fc_probe_clips_runs (reinterpret_cast<double*> (reinterpret_cast<char*> (buf.data()) + 1),
                                       stored * stride) == 0, "a misaligned buffer is refused, not rounded down");
        test::ok (fc_probe_clips_runs (buf.data(), 0xFFFFFFFFu) == stored,
                  "an absurd capacity clamps to stored rather than wrapping the stride multiply");
    }

    // A REFUSED RUN CLEARS THE PREVIOUS RESULT — the getters must not serve the last good file.
    struct Bad { const char* what; std::uint32_t frames, channels, maxRuns; double sr; int off; };
    for (const Bad& b : std::vector<Bad> {
            { "frames = 0",                0u, 2u, 1u << 16, 8000.0, 0 },
            { "channels = 0",   (std::uint32_t) n, 0u, 1u << 16, 8000.0, 0 },
            { "channels past the core's width", (std::uint32_t) n, 99u, 1u << 16, 8000.0, 0 },
            { "sampleRate = 0",  (std::uint32_t) n, 2u, 1u << 16,    0.0, 0 },
            { "sampleRate = NaN",(std::uint32_t) n, 2u, 1u << 16, std::numeric_limits<double>::quiet_NaN(), 0 },
            { "sampleRate = 1e300", (std::uint32_t) n, 2u, 1u << 16, 1e300, 0 },
            { "maxRuns past the ABI bound", (std::uint32_t) n, 2u, (std::uint32_t) fcore::ClipProbe::kMaxRuns + 1u, 8000.0, 0 },
            { "a misaligned planar pointer", (std::uint32_t) (n - 1), 2u, 1u << 16, 8000.0, 1 } })
    {
        test::ok (fc_probe_clips_run (planar.data(), (std::uint32_t) n, (std::uint32_t) ch, 8000.0, 1 << 16) == 1,
                  std::string ("a good run before: ") + b.what);
        const float* base = b.off == 0 ? planar.data()
                                       : reinterpret_cast<const float*> (reinterpret_cast<const char*> (planar.data()) + b.off);
        test::ok (fc_probe_clips_run (base, b.frames, b.channels, b.sr, b.maxRuns) == 0,
                  std::string ("refused: ") + b.what);
        std::vector<double> sink (16, -7.5);
        test::ok (fc_probe_clips_count() == 0 && fc_probe_clips_stored() == 0 && fc_probe_clips_complete() == 0
                  && fc_probe_clips_samples() == 0 && fc_probe_clips_rate() == 0.0
                  && fc_probe_clips_peaks (sink.data(), 16) == 0 && fc_probe_clips_runs (sink.data(), 16) == 0
                  && b64 (sink[0]) == b64 (-7.5),
                  std::string ("…and every getter is zero, with nothing written: ") + b.what);
    }

    // AND A GOOD RUN AFTER A REFUSED ONE WORKS. A refusal clears; it does not wedge.
    test::ok (fc_probe_clips_run (planar.data(), (std::uint32_t) n, (std::uint32_t) ch, 8000.0, 1 << 16) == 1,
              "a good run after all those refusals");
    test::ok ((std::int64_t) fc_probe_clips_count() == want.runCount, "…reports the same numbers as before");

    // maxRuns = 0 is legal: count, store nothing, say so.
    test::ok (fc_probe_clips_run (planar.data(), (std::uint32_t) n, (std::uint32_t) ch, 8000.0, 0) == 1,
              "maxRuns = 0 is a legal run");
    test::ok (fc_probe_clips_count() > 0 && fc_probe_clips_stored() == 0 && fc_probe_clips_complete() == 0,
              "…which counts, stores nothing and reports itself incomplete");

    // A ONE-CHANNEL, ONE-SAMPLE run: the smallest thing the ABI accepts.
    {
        const float one = 0.5f;
        test::ok (fc_probe_clips_run (&one, 1u, 1u, 8000.0, 4) == 1, "one sample, one channel");
        test::ok (fc_probe_clips_count() == 0 && fc_probe_clips_samples() == 1 && fc_probe_clips_complete() == 1,
                  "…no run in it, and the report is whole");
    }
}

int main()
{
    std::printf ("felitronics_clips_exposure_tests\n");
    theFormatIsBitwise();
    lawEightA();
    everyLevelAgainstTheSamples();
    ratesAndTheirWindows();
    pathologicalStreams();
    runsOnTheChunkBoundary();
    aChannelThatDisappears();
    theTraceOfDecisions();
    noCertificateWithoutAMeasurement();
    capacityIsData();
    theAbi();
    return test::report();
}

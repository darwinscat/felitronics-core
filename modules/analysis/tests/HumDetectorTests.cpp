// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// HumDetector self-tests — mains hum against the notes that look like it, and law 8a over the whole report.
//
//   · GEOMETRY pinned by HAND-COMPUTED numbers, not by the function the object sizes itself with: an auto
//     order of 17 at 48 kHz / 18 at 96 kHz / 13 at 4 kHz, and every derived bin extent at the test rate
//   · the MONOTONE-SKIRT oracle: a direct double DFT of a Hann-windowed 49.0 Hz tone, computed outside the
//     object, showing the sampled skirt has NO interior local maximum — which is what lets the detector
//     enumerate every local maximum without reading a note's sidelobe as a 50 Hz line
//   · POSITION accuracy against the frequency that was synthesised: 0.02 bins over a sub-bin sweep
//   · the NEGATIVE notes, G1 = 49.0 / A#1 = 58.3 / B1 = 61.7 Hz, sustained through two quiet stretches
//   · HUM UNDER A LOUDER NOTE — the case that fails if only the strongest peak of the window is examined
//   · 100/120 Hz with NO fundamental: base 2, `fundamentalObserved == false`, the fundamental marked derived
//   · a MOVING tone, both between stretches and inside them (two sweeps whose averages agree)
//   · the validity ladder: NotPrepared / NotFinished / InsufficientResolution / ShorterThanWindow /
//     AllFramesHoled / NoQuietStretch / SingleQuietStretch — a false "clean" is the one forbidden answer
//   · law 8a: the FULL frame trace, every accumulated band bin, every report/candidate/harmonic/stretch field,
//     compared BIT-EXACTLY by std::bit_cast across 14 slicings, seeded ragged partitions, and two maxBlocks
//   · capacity as DATA: maxStretches 1 against 64 classifies identically; the trace overflows the same way
//   · lifecycle: finish idempotent, process refused after it, reset replays bit-identically, no allocation
//
// MUTANT PASS — run, results recorded. 34 of 35 red; six of them were GREEN on an earlier pass and every
// one of those greens was a hole in the SUITE, which is what sections 10b-10k exist for. Each mutant was
// built with the object file DELETED first, and the stand REFUSES a verdict unless a `Building CXX` line
// appeared: a header edited in the same SECOND as the previous build is invisible to make's one-second
// mtime granularity and reads as a false GREEN.
//   ( 1) the frame consumed at the process() BOUNDARY, not inside tick()        -> RED, 102 of 304
//   ( 2) the active stretch closed at the end of every process()                -> RED,   3 of 315
//   ( 3) frameFinite() dropped from the quiet test                              -> RED,   3 of 315
//   ( 4) any bin inside the tolerance, not a strict interior local maximum      -> RED,   5 of 315
//   ( 5) the bin CENTRE returned instead of the interpolated position           -> RED,  25 of 315
//   ( 6) fundamentalObserved = baseFound (asserting an unseen h = 1)            -> RED,   2 of 315
//   ( 7) the candidate summaries frozen once the stretch list is full           -> RED,   6 of 315
//   ( 8) the interior-bin fold's factor 2 removed from the quiet measure        -> RED,   4 of 315  [was GREEN]
//   ( 9) prominence from tonePower over a one-bin floor (a 7.8 dB head start)   -> RED,   1 of 315
//   (10) the per-FRAME spread dropped from the stationarity gate                -> RED,   2 of 315
//   (11) the candidate bands left INSIDE the quiet gate                         -> RED,   4 of 315
//   (12) SingleQuietStretch downgraded to a verdict                             -> RED,   1 of 315
//   (13) the resolution gate removed                                            -> RED,   2 of 315
//   (14) the floor as a one-sided median, not paired geometric means            -> RED,   2 of 315  [was GREEN]
//   (15) normalise() skipped: published powers left as accumulated sums         -> RED,   3 of 315
//   (16) the search window narrowed to the acceptance tolerance                 -> RED,   9 of 315
//   (17) the exclusion list built WITHOUT merging, so it is unsorted            -> RED,   2 of 315  [was GREEN]
//   (18) the comb searched on the NOMINAL grid even when a base was found       -> RED,   3 of 315  [was GREEN]
//   (19) the stretch's own average replaced by the pooled band                  -> RED,   6 of 315
//   (20) the open stretch at the tail not closed by finish()                    -> RED,  81 of 309
//   (21) the STRETCH-spread half of the stationarity gate dropped               -> GREEN, and equivalent —
//        see the note below; the frame-spread gate is the binding one on every input reachable
//   (22) the prominence gate read as an amplitude ratio (10 dB becomes 5 dB)    -> RED,   1 of 315
//   (23) the ABSOLUTE level floor removed (a ratio alone certifies residue)     -> RED,   9 of 315
//   (24) the frames-per-observation rule removed                               -> RED,   2 of 315
//   (25) a digital-silence frame counted as quiet                              -> RED,   2 of 315
//   (26) a found-but-not-stationary candidate reported as clean                -> RED,   3 of 315
//   (27) the comb-without-base incompleteness removed                          -> RED,   1 of 315
//   (28) ONE harmonic at a multiple counted as a comb                          -> RED,   2 of 315
//   (29) the summary channel's comparison flipped                              -> RED,   1 of 315  [was GREEN]
//   (30) the winner tie-break disabled                                         -> RED,   2 of 315  [was GREEN]
//   (31) the trace's stretch index off by one                                  -> RED,   1 of 315  [was GREEN]
//   (32) the resolution bound back to a round 0.5 Hz (2.048 bins at 32 kHz)    -> RED,   4 of 315
//   (33) the lower median replaced by the upper (even populations only)        -> RED,   1 of 315  [was GREEN]
//   (34) tonePower loses the fold factor 2                                     -> RED,   1 of 315  [was GREEN]
//   (35) stretchesComplete off by one at the capacity boundary                 -> RED,   1 of 315  [was GREEN]
//
// WHY (21) IS AN EQUIVALENCE AND NOT A HOLE, measured rather than argued. Averaging frames raises a line and
// its background together, so the peak-to-floor RATIO barely improves with frame count — averaging buys the
// stability of the estimate and the accuracy of its position, not sensitivity (sensitivity comes from a
// longer WINDOW, which narrows the bin and so lowers the per-bin background while the line keeps its power).
// The consequence is that stretch observations and frame observations appear and disappear together: swept
// from amplitude 1.2e-3 down to 8e-5, at 3 frames per stretch and again at 16, the two counts crossed at the
// same step every time, and the frame spread was WIDER than the stretch spread in every row (0.58 vs 0.57 at
// 3 frames, 0.67 vs 0.57 at 16). So no input reached makes the stretch-spread gate the deciding one. It is
// kept because it is the gate that states the brief's requirement on the ROBUST estimate, and because the
// redundancy is an empirical property of these parameters rather than a theorem.
//
// A SEVENTH mutant came back green and was NOT a hole: deleting the gap check on consecutive selected frame
// indices changed nothing, because every closed frame reaches onFrame() and every path but "quiet and
// finite" has already closed the stretch. That check could not fire, so it was deleted from the header
// rather than kept as a gate that cannot redden.
// Mutant (2) is this instrument's stand-in for "the denormal flush moved to the call boundary": HumDetector
// holds no IIR state and no per-sample recurrence, so core::StateGrid has nothing to schedule here, and the
// state a call boundary could corrupt is the open stretch.

#include <felitronics_test.h>
#include <alloc_counter.h>   // installs the allocation counter: EVERY form of `new`, over-aligned included
#include <felitronics/analysis/HumDetector.h>

#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>
#include <random>
#include <string>
#include <vector>

using namespace felitronics;
using analysis::HumDetector;
using analysis::HumDetectorParams;
using analysis::HumMains;
using analysis::HumReason;
using analysis::HumReport;

namespace
{

// 3000 Hz is chosen so that the AUTO order (13, N = 8192) gives bin = 3000/8192 = 0.3662109375 Hz — bit for
// bit the bin of the production configuration (48000 / 2^17) and the same 2.731 s window, at a sixteenth of
// the FFT cost. The suite therefore runs at the resolution the instrument actually ships with, and the
// 49.0-vs-50.0 separation it is tested on is the real 2.73 bins rather than an easier number.
constexpr double kSr = 3000.0;
constexpr double kBin = kSr / 8192.0;

//------------------------------------------------------------------------------------------------------------
// THE MATERIAL. Planar float, built by adding tones and noise over sample ranges, so a "quiet stretch" is an
// explicit interval of the fixture rather than an emergent property nobody checked.
struct Sig
{
    std::vector<std::vector<float>> ch;
    Sig (int nch, std::size_t len) : ch ((std::size_t) nch, std::vector<float> (len, 0.0f)) {}
    std::size_t len() const { return ch[0].size(); }
    int channels() const { return (int) ch.size(); }

    void tone (std::size_t from, std::size_t to, double hz, double amp, double phase = 0.0)
    {
        for (std::size_t i = from; i < std::min (to, len()); ++i)
            for (auto& c : ch)
                c[i] += (float) (amp * std::sin (2.0 * core::kPi * hz * (double) i / kSr + phase));
    }
    // a linear sweep, for the "a note moves with the music" fixtures
    void sweep (std::size_t from, std::size_t to, double hzA, double hzB, double amp)
    {
        const double span = (double) (to - from);
        double ph = 0.0;
        for (std::size_t i = from; i < std::min (to, len()); ++i)
        {
            const double f = hzA + (hzB - hzA) * ((double) (i - from) / span);
            ph += 2.0 * core::kPi * f / kSr;
            for (auto& c : ch) c[i] += (float) (amp * std::sin (ph));
        }
    }
    void noise (std::size_t from, std::size_t to, double amp, unsigned seed)
    {
        std::mt19937 rng (seed);
        std::uniform_real_distribution<double> u (-1.0, 1.0);
        for (std::size_t i = from; i < std::min (to, len()); ++i)
            for (auto& c : ch) c[i] += (float) (amp * u (rng));
    }
    // Steeply tilted background: white noise through `poles` cascaded one-poles, scaled to an exact RMS,
    // so the level is a parameter and the SLOPE across the floor window is the thing being tested.
    void redNoise (std::size_t from, std::size_t to, double rms, double cutoffHz, int poles, unsigned seed)
    {
        const std::size_t hi2 = std::min (to, len());
        if (hi2 <= from) return;
        const std::size_t n = hi2 - from;
        std::mt19937 rng (seed);
        std::uniform_real_distribution<double> u (-1.0, 1.0);
        std::vector<double> t (n);
        for (auto& v : t) v = u (rng);
        const double a = std::exp (-2.0 * core::kPi * cutoffHz / kSr);
        for (int pp = 0; pp < poles; ++pp)
        {
            double z = 0.0;
            for (auto& v : t) { z = a * z + (1.0 - a) * v; v = z; }
        }
        double ss = 0.0;
        for (double v : t) ss += v * v;
        const double g = ss > 0.0 ? rms / std::sqrt (ss / (double) n) : 0.0;
        for (std::size_t i = 0; i < n; ++i) for (auto& c : ch) c[from + i] += (float) (g * t[i]);
    }
    void dc (std::size_t from, std::size_t to, double d)
    {
        for (std::size_t i = from; i < std::min (to, len()); ++i) for (auto& c : ch) c[i] += (float) d;
    }
    std::vector<const float*> ptrs (std::size_t at) const
    {
        std::vector<const float*> p;
        for (const auto& c : ch) p.push_back (c.data() + at);
        return p;
    }
};

// Geometry of the standard two-stretch fixture: three quiet frames, a loud gap, three more quiet frames.
constexpr std::size_t kN = 8192, kHop = 4096;
constexpr std::size_t kQuietA0 = 0,                  kQuietA1 = kN + 2 * kHop;          // 16384
constexpr std::size_t kLoud0   = kQuietA1,           kLoud1   = kQuietA1 + 4 * kHop;    // 32768
constexpr std::size_t kQuietB0 = kLoud1,             kQuietB1 = kLoud1 + kN + 2 * kHop; // 49152
constexpr std::size_t kLen     = kQuietB1;

// A fixture whose quiet stretches carry `hz` (plus an exact comb when `harmonics`), under a noise bed.
Sig twoStretch (int nch, double hz, double amp, bool harmonics, unsigned seed = 7)
{
    Sig s (nch, kLen);
    s.noise (kQuietA0, kQuietA1, 3.0e-4, seed);
    s.noise (kQuietB0, kQuietB1, 3.0e-4, seed + 1u);
    s.noise (kLoud0, kLoud1, 0.3, seed + 2u);                       // the programme: ~ -15 dBFS
    if (amp > 0.0)
    {
        s.tone (kQuietA0, kQuietA1, hz, amp);
        s.tone (kQuietB0, kQuietB1, hz, amp, 0.7);                  // a different phase: no coherent luck
        if (harmonics)
            for (int h = 2; h <= 4; ++h)
            {
                const double a = amp / (double) (h * h);
                s.tone (kQuietA0, kQuietA1, hz * (double) h, a);
                s.tone (kQuietB0, kQuietB1, hz * (double) h, a, 0.3);
            }
    }
    return s;
}

//------------------------------------------------------------------------------------------------------------
// Run a signal through the detector under a given slicing, and serialise EVERYTHING it published as bits.
// Not the final report alone: docs/LAW8-KWEIGHTING.md:78 is the project's own measured case of identical
// scalars over five moved intermediates, and section 3's re-slicing compares the object with itself, so the
// per-frame trace and every accumulated bin have to be in here too.
struct Out
{
    std::vector<std::uint64_t> bits;
    bool ok = false;
    void put (std::uint64_t v) { bits.push_back (v); }
    void put (double d) { bits.push_back (std::bit_cast<std::uint64_t> (d)); }
    void put (std::int64_t v) { bits.push_back (std::bit_cast<std::uint64_t> (v)); }
    void put (int v) { bits.push_back ((std::uint64_t) (std::int64_t) v); }
    void put (bool v) { bits.push_back (v ? 1u : 0u); }
    void put (const analysis::HumPeak& p)
    {
        put (p.found); put (p.prominent); put (p.accepted); put (p.bin);
        put (p.hz); put (p.tonePower); put (p.peakBinPower); put (p.floorPower); put (p.prominenceDb);
    }
    bool operator== (const Out& o) const { return ok == o.ok && bits == o.bits; }
    bool operator!= (const Out& o) const { return ! (*this == o); }
};

// Everything the detector published, as bits. Split out of runSliced so a reused object can be compared
// against a fresh one FIELD BY FIELD — the replay check used to compare a single frequency, and a mutant that
// dropped `traceCount_ = 0` from reset() survived all 198 checks because of it.
void serialise (const HumDetector& d, int nch, int maxHarmonic, Out& o)
{
    o.put ((std::int64_t) d.traceCount());
    o.put (d.traceComplete());
    for (std::int64_t i = 0; i < d.storedTraceCount(); ++i)
    {
        const auto t = d.traceAt (i);
        o.put (t.frameStart); o.put (t.frameEnd); o.put (t.stretchIndex);
        o.put (t.meanSquare); o.put (t.bandSum);
        for (int h = 0; h < 4; ++h) { o.put (t.hypHz[h]); o.put (t.hypTone[h]); o.put (t.hypFloor[h]); o.put ((int) t.hypBin[h]); }
        o.put ((std::uint64_t) t.hypAccepted); o.put ((int) t.channel);
        o.put ((std::uint64_t) t.finite); o.put ((std::uint64_t) t.quiet);
    }
    for (int c = 0; c < nch; ++c)
    {
        const HumReport r = d.report (c);
        o.put (r.valid); o.put ((int) r.reason); o.put ((int) r.mains); o.put (r.line);
        o.put (r.fundamentalHz); o.put (r.fundamentalObserved); o.put (r.fundamentalDerived);
        o.put (r.baseHarmonic); o.put (r.harmonicsObserved);
        o.put (r.frames); o.put (r.finiteFrames); o.put (r.holedFrames); o.put (r.quietFrames);
        o.put (r.quietFrames); o.put (r.silentFrames);
        o.put (r.quietStretches); o.put (r.eligibleStretches); o.put (r.storedStretches);
        o.put (r.stretchesComplete);
        o.put (r.totalSamples); o.put (r.tailUncoveredSamples); o.put (r.binHz); o.put (r.windowSamples);
        for (int cand = 0; cand < HumDetector::kCandidates; ++cand)
        {
            const auto k = d.candidate (c, cand);
            o.put (k.nominalHz); o.put (k.baseFound); o.put (k.baseHarmonic); o.put (k.fundamentalHz);
            o.put (k.fundamentalObserved); o.put (k.base); o.put (k.windowPeak);
            o.put (k.harmonicsObserved); o.put (k.lowestHarmonicObserved);
            o.put (k.combWithoutBase); o.put (k.combFundamentalHz);
            o.put (k.stretchObservations); o.put (k.stretchOffTolerance); o.put (k.frameObservations);
            o.put (k.stretchSpreadHz); o.put (k.frameSpreadHz); o.put (k.maxIntraStretchSpreadHz);
            o.put (k.stationary); o.put (k.passed);
            for (int h = 1; h <= maxHarmonic; ++h)
            {
                const auto hh = d.harmonic (c, cand, h);
                o.put (hh.index); o.put (hh.inBand); o.put (hh.peak);
            }
        }
        o.put (d.stretchCount (c)); o.put (d.stretchesComplete (c));
        for (std::int64_t i = 0; i < d.storedStretchCount (c); ++i)
        {
            const auto st = d.stretch (c, i);
            o.put (st.index); o.put (st.startSample); o.put (st.endSample); o.put (st.frames);
        }
        o.put (d.selectedFrames (c));
        const double* band = d.bandPowerSum (c);
        for (int k = 0; k < d.bandBins(); ++k) o.put (band[k]);     // EVERY accumulated bin, not a checksum
    }
}

Out runSliced (const Sig& s, const HumDetectorParams& p, const std::vector<int>& slices, int maxBlock)
{
    Out o;
    HumDetector d;
    d.setParams (p);
    if (! d.prepare (kSr, maxBlock, s.channels())) { o.put (std::uint64_t (0xDEAD)); return o; }
    const std::size_t total = s.len();
    std::size_t at = 0, si = 0;
    bool allOk = true;
    while (at < total)
    {
        const int want = slices[si % slices.size()]; ++si;
        const std::size_t take = want <= 0 ? 0u : std::min ((std::size_t) want, total - at);
        const auto v = s.ptrs (at);
        if (! d.process (v.data(), s.channels(), (int) take)) allOk = false;
        at += take;
    }
    d.finish();
    o.ok = allOk;
    serialise (d, s.channels(), p.maxHarmonic, o);
    return o;
}

// The whole signal in one call — the reference run, and the one the assertions read.
HumReport once (const Sig& s, const HumDetectorParams& p, HumDetector& d)
{
    d.setParams (p);                                           // the helper used to drop `p` on the floor
    if (! d.prepare (kSr, 4096, s.channels())) return HumReport {};
    const auto v = s.ptrs (0);
    if (! d.process (v.data(), s.channels(), (int) s.len())) return HumReport {};
    d.finish();
    return d.report (0);
}

std::string hz (double v) { char b[40]; std::snprintf (b, sizeof b, "%.4f", v); return b; }

} // namespace

int main()
{
    using felitronics::test::ok;
    using felitronics::test::approx;

    //---------------------------------------------------------------------------------------------------
    // 1. GEOMETRY, pinned by hand. A budget that sizes itself through the same function it is compared with
    //    proves nothing (the repository has that lesson written down), so these are numbers computed here.
    {
        HumDetectorParams p;                                        // defaults, fftOrder = 0 = AUTO
        const auto g48 = HumDetector::geometryFor (48000.0, p);
        ok (g48.ok && g48.order == 17, "geometry: AUTO picks order 17 at 48 kHz (got " + std::to_string (g48.order) + ")");
        approx (g48.binHz, 48000.0 / 131072.0, 1e-15, "geometry: bin is 0.3662 Hz at 48 kHz");
        ok (g48.binHz <= HumDetector::kMaxBinHz, "geometry: and that is inside the 0.5 Hz resolution bound");
        const auto g96 = HumDetector::geometryFor (96000.0, p);
        ok (g96.ok && g96.order == 18, "geometry: AUTO picks order 18 at 96 kHz (got " + std::to_string (g96.order) + ")");
        const auto g441 = HumDetector::geometryFor (44100.0, p);
        ok (g441.ok && g441.order == 17, "geometry: AUTO picks order 17 at 44.1 kHz");
        const auto g = HumDetector::geometryFor (kSr, p);
        ok (g.ok && g.order == 13 && g.n == 8192 && g.hop == 4096, "geometry: order 13 / N 8192 / hop 4096 at 4 kHz");
        ok (g.bins == 4097, "geometry: 4097 bins");
        // topHz = (60 + 0.5)*8 + 2*bin + 10 + bin = 495.0986..., / bin = 1351.9 -> 1352
        ok (g.bandLo == 2 && g.bandHi == 1352 && g.bandBins == 1351,
            "geometry: the accumulated band is bins 2..1015 (got " + std::to_string (g.bandLo) + ".."
            + std::to_string (g.bandHi) + ")");
        // margin = 3 + 10 + bin = 13.366; (50 - margin)/bin = 100.03 -> 100; (120 + margin)/bin = 364.2 -> 365
        ok (g.stretchLo == 100 && g.stretchHi == 365 && g.stretchBins == 266,
            "geometry: the active-stretch band is bins 74..274 (got " + std::to_string (g.stretchLo) + ".."
            + std::to_string (g.stretchHi) + ")");
        ok (g.floorSpanBins == 28 && g.floorExcludeBins == 5 && g.floorPairs == 23,
            "geometry: 28 background bins minus 5 excluded leaves 23 pairs");
        ok (g.searchBins == 9, "geometry: the search window is 9 bins each side");
        ok (g.spanHarmonics == 8, "geometry: the span covers all 8 reported harmonics");
        // the resolution bound is a SEPARATION, and this rate meets it with the production margin
        approx (HumDetector::kMinNoteSeparationHz / g.binHz, 2.731, 0.001,
                "geometry: 49.0 and 50.0 Hz are 2.73 bins apart here, as at 48 kHz / 2^17");
        HumDetectorParams coarse = p; coarse.fftOrder = 12;   // bin 0.732: 1.37 bins, a merged pair
        ok (! HumDetector::geometryFor (kSr, coarse).ok
            || HumDetector::geometryFor (kSr, coarse).binHz > HumDetector::kMaxBinHz,
            "geometry: one order coarser does not meet the separation");

        // the same numbers reach Storage, and Storage is published before anything is allocated
        const auto st = HumDetector::storageFor (kSr, 2, p);
        ok (st.ok && st.bandDoubles == 2702 && st.stretchDoubles == 532 && st.floorDoubles == 23,
            "storage: published from the same geometry (" + std::to_string (st.bandDoubles) + " band doubles)");
        ok (st.stretchEntries == 512 && st.harmonicEntries == 32 && st.channelEntries == 2,
            "storage: stretch, harmonic and channel rows are sized by the parameters");
        ok (st.bytes() >= st.frames.bytes()
                            + (std::uint64_t) sizeof (double) * (st.bandDoubles + st.stretchDoubles + st.floorDoubles),
            "storage: the total covers the producer's budget plus every double row it publishes");

        // refusals
        HumDetectorParams bad = p;
        ok (! HumDetector::geometryFor (0.0, p).ok, "refuse: sample rate zero");
        ok (! HumDetector::geometryFor (std::numeric_limits<double>::quiet_NaN(), p).ok, "refuse: rate NaN");
        ok (! HumDetector::geometryFor (999.0, p).ok, "refuse: rate below the bound");
        ok (! HumDetector::geometryFor (768001.0, p).ok, "refuse: rate above the bound");
        bad = p; bad.fftOrder = 3;            ok (! HumDetector::geometryFor (kSr, bad).ok, "refuse: order 3");
        bad = p; bad.fftOrder = 23;           ok (! HumDetector::geometryFor (kSr, bad).ok, "refuse: order 23");
        bad = p; bad.hop = -1;                ok (! HumDetector::geometryFor (kSr, bad).ok, "refuse: negative hop");
        bad = p; bad.hop = 8193;              ok (! HumDetector::geometryFor (kSr, bad).ok, "refuse: hop past the window");
        bad = p; bad.toleranceHz = 0.0;       ok (! HumDetector::geometryFor (kSr, bad).ok, "refuse: zero tolerance");
        bad = p; bad.searchHz = 0.1;          ok (! HumDetector::geometryFor (kSr, bad).ok, "refuse: search inside the tolerance");
        bad = p; bad.floorExcludeHz = 20.0;   ok (! HumDetector::geometryFor (kSr, bad).ok, "refuse: exclusion wider than the span");
        bad = p; bad.maxHarmonic = 0;         ok (! HumDetector::geometryFor (kSr, bad).ok, "refuse: no harmonic");
        bad = p; bad.maxHarmonic = 65;        ok (! HumDetector::geometryFor (kSr, bad).ok, "refuse: too many harmonics");
        bad = p; bad.maxStretches = -1;       ok (! HumDetector::geometryFor (kSr, bad).ok, "refuse: negative capacity");
        bad = p; bad.quietThresholdDb = std::numeric_limits<double>::infinity();
        ok (! HumDetector::geometryFor (kSr, bad).ok, "refuse: threshold not finite");
        HumDetector d;
        ok (! d.prepare (kSr, 64, 0), "refuse: zero channels");
        ok (! d.prepare (kSr, 64, core::kMaxChannels + 1), "refuse: too many channels");
        ok (d.report (0).reason == HumReason::NotPrepared, "refuse: the report of an unprepared object says so");
        ok (d.prepare (kSr, 64, 1), "refuse: and a good call is still accepted afterwards");
        ok (d.report (0).reason == HumReason::NotFinished, "lifecycle: before finish() the report is NotFinished");
        ok (d.report (7).reason == HumReason::NotPrepared, "lifecycle: a channel outside the width answers NotPrepared");
    }

    //---------------------------------------------------------------------------------------------------
    // 2. THE MONOTONE-SKIRT ORACLE, computed outside the object. The detector enumerates EVERY strict
    //    interior local maximum in its search window instead of taking the argmax, which is only safe if a
    //    lone Hann-windowed tone has no interior maximum on the bin grid. |W(m+f)| holds |sin(pi f)| fixed
    //    across bins and falls as 1/|m+f|^3, so it must not. Direct double DFT of the windowed 49.0 Hz tone.
    {
        const int n = 8192;
        std::vector<double> w ((std::size_t) n);
        for (int i = 0; i < n; ++i) w[(std::size_t) i] = 0.5 - 0.5 * std::cos (2.0 * core::kPi * i / (double) n);
        const int lo = (int) std::floor (47.0 / kBin), hi = (int) std::ceil (53.0 / kBin);
        std::vector<double> mag ((std::size_t) (hi - lo + 1), 0.0);
        for (int k = lo; k <= hi; ++k)
        {
            double re = 0.0, im = 0.0;
            for (int t = 0; t < n; ++t)
            {
                const double v = w[(std::size_t) t] * std::sin (2.0 * core::kPi * 49.0 * (double) t / kSr);
                const double a = -2.0 * core::kPi * (double) k * (double) t / (double) n;
                re += v * std::cos (a); im += v * std::sin (a);
            }
            mag[(std::size_t) (k - lo)] = re * re + im * im;
        }
        const int peak = (int) std::llround (49.0 / kBin);
        int interior = 0;
        for (int k = lo + 1; k < hi; ++k)
        {
            const std::size_t j = (std::size_t) (k - lo);
            if (mag[j] > mag[j - 1] && mag[j] >= mag[j + 1] && std::abs (k - peak) > 2) ++interior;
        }
        ok (interior == 0, "skirt oracle: a lone 49.0 Hz tone has " + std::to_string (interior)
                           + " interior maxima outside its main lobe (must be 0)");
        // ... the fall is MONOTONE from the main lobe outwards, which is the property being relied on
        int nonMonotone = 0;
        for (int k = peak + 2; k < hi; ++k)
        {
            const std::size_t j = (std::size_t) (k - lo);
            if (! (mag[j] > mag[j + 1])) ++nonMonotone;
        }
        ok (nonMonotone == 0, "skirt oracle: the skirt falls monotonically out to 53 Hz ("
                              + std::to_string (nonMonotone) + " steps do not)");
        // ... and it really does reach the 50 Hz tolerance band, so the claim above is not vacuous
        const int at50 = (int) std::llround (50.0 / kBin);
        const double share = mag[(std::size_t) (at50 - lo)] / mag[(std::size_t) (peak - lo)];
        ok (share > 1e-6 && share < 0.05, "skirt oracle: the 49 Hz tone reaches the bin nearest 50 Hz at "
                                          + hz (10.0 * std::log10 (share)) + " dB — real, and not a peak");
    }

    //---------------------------------------------------------------------------------------------------
    // 3. POSITION, against the frequency that was synthesised — a sub-bin sweep. The bin centre nearest
    //    50 Hz is 102 * 0.48828125 = 49.805 Hz, so "50.000" can only come from the interpolator.
    {
        double worst = 0.0;
        bool worstFinite = true;
        int detected = 0;
        for (double f : { 49.60, 49.81, 50.00, 50.14, 50.29, 50.45 })
        {
            HumDetectorParams p;
            HumDetector d;
            const auto s = twoStretch (1, f, 1.0e-3, true);
            const HumReport r = once (s, p, d);
            if (r.mains != HumMains::Hz50) continue;
            ++detected;
            const double e = std::fabs (r.line.hz - f);
            if (! std::isfinite (e)) worstFinite = false;
            if (e > worst) worst = e;
        }
        ok (detected == 6, "position: all six sub-bin offsets were detected (" + std::to_string (detected) + "/6)");
        ok (worstFinite, "position: every error is finite");
        ok (worst < 0.02 * kBin, "position: worst error " + hz (worst) + " Hz is under 0.02 bins ("
                                 + hz (0.02 * kBin) + " Hz)");
    }

    //---------------------------------------------------------------------------------------------------
    // 4. REAL HUM: 50.14 Hz with an exact comb, in two quiet stretches separated by programme.
    {
        HumDetectorParams p;
        HumDetector d;
        const auto s = twoStretch (1, 50.14, 1.0e-3, true);
        const HumReport r = once (s, p, d);
        ok (r.valid && r.reason == HumReason::Ok, "hum: the report is valid");
        ok (r.mains == HumMains::Hz50, "hum: mains reads 50 Hz");
        approx (r.line.hz, 50.14, 0.01, "hum: the peak is where it was synthesised");
        ok (r.fundamentalObserved && ! r.fundamentalDerived && r.baseHarmonic == 1,
            "hum: the fundamental itself was seen, so nothing is derived");
        ok (r.line.prominenceDb > 20.0, "hum: it stands " + hz (r.line.prominenceDb) + " dB over the local background");
        ok (r.harmonicsObserved >= 3, "hum: " + std::to_string (r.harmonicsObserved) + " harmonics of the comb were found");
        ok (r.quietStretches == 2, "hum: exactly two quiet stretches (" + std::to_string (r.quietStretches) + ")");
        ok (r.storedStretches == 2 && r.stretchesComplete, "hum: and both are in the list");
        // The coordinates are EXACT, computed from the fixture: frames start every 4096 and span 8192, the
        // loud gap is [16384, 32768), so the last quiet frame of A starts at 8192 and the first of B at 32768.
        const auto s0 = d.stretch (0, 0), s1 = d.stretch (0, 1);
        ok (s0.index == 0 && s0.startSample == 0 && s0.endSample == 16384 && s0.frames == 3,
            "hum: stretch 0 is exactly [0, 16384) over 3 frames (got [" + std::to_string (s0.startSample) + ", "
            + std::to_string (s0.endSample) + ") over " + std::to_string (s0.frames) + ")");
        ok (s1.index == 1 && s1.startSample == 32768 && s1.endSample == 49152 && s1.frames == 3,
            "hum: stretch 1 is exactly [32768, 49152) over 3 frames (got [" + std::to_string (s1.startSample)
            + ", " + std::to_string (s1.endSample) + ") over " + std::to_string (s1.frames) + ")");
        ok (r.quietFrames == 6 && r.eligibleStretches == 2,
            "hum: 6 quiet frames, both stretches eligible (" + std::to_string (r.quietFrames) + ", "
            + std::to_string (r.eligibleStretches) + ")");
        // the published level is a per-frame AVERAGE, not the accumulated sum: A^2/2 for a 1e-3 sine
        approx (r.line.tonePower, 5.0e-7, 7.0e-8, "hum: tonePower is the tone's mean square A^2/2 = 5e-7 ("
                                                  + hz (r.line.tonePower * 1.0e7) + "e-7)");
        const auto k50 = d.candidate (0, 0), k60 = d.candidate (0, 1);
        ok (k50.stretchObservations == 2, "hum: the base was accepted in both stretches ("
                                          + std::to_string (k50.stretchObservations) + ")");
        ok (k50.stationary && k50.passed, "hum: and it stood still");
        ok (k50.stretchSpreadHz <= p.toleranceHz && k50.frameSpreadHz <= p.toleranceHz,
            "hum: stretch spread " + hz (k50.stretchSpreadHz) + " Hz, frame spread " + hz (k50.frameSpreadHz) + " Hz");
        ok (! k60.passed, "hum: the 60 Hz hypothesis did not pass");
        ok (d.strongestChannel() == 0, "hum: the summary channel is the only channel");
        // the comb's harmonics sit on EXACT multiples of the measured fundamental
        int exact = 0;
        for (int h = 2; h <= 4; ++h)
        {
            const auto hh = d.harmonic (0, h);
            if (hh.peak.accepted && std::fabs (hh.peak.hz - 50.14 * (double) h) < 0.05) ++exact;
        }
        ok (exact == 3, "hum: harmonics 2..4 are at exact multiples (" + std::to_string (exact) + "/3)");
    }

    //---------------------------------------------------------------------------------------------------
    // 5. THE NEGATIVE NOTES. A sustained bass note in the quiet is not hum: G1 = 49.0 Hz is 1.0 Hz from 50,
    //    A#1 = 58.3 and B1 = 61.7 are 1.7 Hz from 60. Each is given a full harmonic series, which is what
    //    makes the comb useless as a discriminator and the POSITION decisive.
    {
        for (double note : { 49.0, 58.3, 61.7 })
        {
            HumDetectorParams p;
            HumDetector d;
            const auto s = twoStretch (1, note, 1.5e-3, true);
            const HumReport r = once (s, p, d);
            ok (r.valid, "note " + hz (note) + ": the report is valid — we could look");
            ok (r.mains == HumMains::None, "note " + hz (note) + ": mains reads none");
            ok (r.line.hz == 0.0 && r.harmonicsObserved == 0,
                "note " + hz (note) + ": and no conclusion field was filled in");
            // the evidence is still published: the window's strongest peak names the note
            const int cand = note < 55.0 ? 0 : 1;
            const auto k = d.candidate (0, cand);
            ok (k.windowPeak.found, "note " + hz (note) + ": the window peak is published as evidence");
            approx (k.windowPeak.hz, note, 0.05, "note " + hz (note) + ": and it is the note, not 50 or 60");
            ok (! k.baseFound, "note " + hz (note) + ": no comb base was accepted");
        }
    }

    //---------------------------------------------------------------------------------------------------
    // 6. HUM UNDER A LOUDER NOTE — the case that fails when only the strongest peak of the window is read.
    //    A 49.0 Hz note 5x the amplitude of a 50.15 Hz hum: the note owns the window's argmax and fails the
    //    tolerance; the hum is the SECOND local maximum, 2.36 bins away, and is what the report must name.
    {
        HumDetectorParams p;
        HumDetector d;
        Sig s (1, kLen);
        s.noise (kQuietA0, kQuietA1, 3.0e-4, 21);
        s.noise (kQuietB0, kQuietB1, 3.0e-4, 22);
        s.noise (kLoud0, kLoud1, 0.3, 23);
        s.tone (kQuietA0, kQuietA1, 49.0, 2.0e-3);  s.tone (kQuietB0, kQuietB1, 49.0, 2.0e-3, 1.1);
        s.tone (kQuietA0, kQuietA1, 98.0, 6.0e-4);  s.tone (kQuietB0, kQuietB1, 98.0, 6.0e-4, 0.4);
        s.tone (kQuietA0, kQuietA1, 50.15, 4.0e-4); s.tone (kQuietB0, kQuietB1, 50.15, 4.0e-4, 0.9);
        s.tone (kQuietA0, kQuietA1, 100.30, 1.5e-4); s.tone (kQuietB0, kQuietB1, 100.30, 1.5e-4, 0.2);
        const HumReport r = once (s, p, d);
        ok (r.valid && r.mains == HumMains::Hz50, "masked hum: found under a note 14 dB louder 1.15 Hz away");
        // the louder neighbour 2.36 bins away pulls the three-bin estimate by ~0.07 bins; that is the
        // interference, not the estimator, and it is an order under the 0.5 Hz tolerance that decides anything
        approx (r.line.hz, 50.15, 0.05, "masked hum: at the frequency it was synthesised (" + hz (r.line.hz) + ")");
        const auto k = d.candidate (0, 0);
        approx (k.windowPeak.hz, 49.0, 0.05, "masked hum: the window's STRONGEST peak is still the note");
        ok (k.windowPeak.tonePower > k.base.tonePower * 4.0,
            "masked hum: the note really is the louder of the two, so the argmax rule would have missed this");
    }

    //---------------------------------------------------------------------------------------------------
    // 7. 100/120 Hz WITH NO FUNDAMENTAL — full-wave rectification. The comb must be reported and the
    //    unobserved fundamental must NOT be asserted as seen.
    {
        HumDetectorParams p;
        HumDetector d;
        Sig s (1, kLen);
        s.noise (kQuietA0, kQuietA1, 3.0e-4, 31);
        s.noise (kQuietB0, kQuietB1, 3.0e-4, 32);
        s.noise (kLoud0, kLoud1, 0.3, 33);
        for (double f : { 120.0, 240.0, 360.0 })
        {
            const double a = 1.0e-3 * 120.0 / f;
            s.tone (kQuietA0, kQuietA1, f, a);
            s.tone (kQuietB0, kQuietB1, f, a, 0.5);
        }
        const HumReport r = once (s, p, d);
        ok (r.valid && r.mains == HumMains::Hz60, "rectifier: a 120/240/360 comb reads 60 Hz mains");
        ok (r.baseHarmonic == 2, "rectifier: the comb's base is h = 2 (got " + std::to_string (r.baseHarmonic) + ")");
        ok (! r.fundamentalObserved, "rectifier: the 60 Hz fundamental was NOT observed");
        ok (r.fundamentalDerived, "rectifier: and the fundamental it reports is marked derived");
        approx (r.fundamentalHz, 60.0, 0.02, "rectifier: derived as 120/2 = " + hz (r.fundamentalHz) + " Hz");
        ok (! d.harmonic (0, 1).peak.accepted, "rectifier: harmonic 1 is not accepted");
        ok (d.harmonic (0, 2).peak.accepted && d.harmonic (0, 4).peak.accepted,
            "rectifier: harmonics 2 and 4 are");
        approx (d.harmonic (0, 2).peak.hz, 120.0, 0.02, "rectifier: h = 2 sits at 120 Hz");
    }

    //---------------------------------------------------------------------------------------------------
    // 8. A NOTE MOVES. Between stretches (50.1 then 51.9) it gives ONE observation, not two. Inside them,
    //    two symmetric sweeps whose stretch averages BOTH land on 50.0 are what the per-frame spread is for.
    {
        {
            HumDetectorParams p;
            HumDetector d;
            Sig s (1, kLen);
            s.noise (kQuietA0, kQuietA1, 3.0e-4, 41);
            s.noise (kQuietB0, kQuietB1, 3.0e-4, 42);
            s.noise (kLoud0, kLoud1, 0.3, 43);
            s.tone (kQuietA0, kQuietA1, 50.10, 1.0e-3);
            s.tone (kQuietB0, kQuietB1, 51.90, 1.0e-3);              // a semitone-ish step: the note moved
            const HumReport r = once (s, p, d);
            // The pooled average DOES carry a 50 Hz-compatible line (stretch A's tone), so "no hum" would be
            // a false clean. The answer is the named incompleteness instead.
            ok (! r.valid && r.reason == HumReason::CandidateNotStationary,
                "moving note: a line was measured and did not stand still -> CandidateNotStationary");
            ok (r.mains == HumMains::None, "moving note: and no mains is claimed");
            const auto k = d.candidate (0, 0);
            ok (k.baseFound, "moving note: the candidate's evidence is still published");
            ok (k.stretchObservations == 1, "moving note: exactly one stretch showed a 50 Hz-compatible line ("
                                            + std::to_string (k.stretchObservations) + ")");
            ok (k.stretchOffTolerance == 1, "moving note: and the other stretch's off-tolerance peak is counted");
            ok (! k.stationary, "moving note: so it is not stationary");
        }
        {
            HumDetectorParams p;
            HumDetector d;
            Sig s (1, kLen);
            s.noise (kQuietA0, kQuietA1, 3.0e-4, 44);
            s.noise (kQuietB0, kQuietB1, 3.0e-4, 45);
            s.noise (kLoud0, kLoud1, 0.3, 46);
            s.sweep (kQuietA0, kQuietA1, 49.2, 50.8, 1.2e-3);        // both average to 50.0
            s.sweep (kQuietB0, kQuietB1, 50.8, 49.2, 1.2e-3);
            HumDetector dd;
            const HumReport r = once (s, p, dd);
            const auto k = dd.candidate (0, 0);
            ok (! r.valid && r.reason == HumReason::CandidateNotStationary,
                "sweeping note: a swept line is an incompleteness, not a clean report");
            ok (k.baseFound, "sweeping note: the pooled average DOES show a base — the sweeps average to 50.0");
            ok (k.stretchObservations >= 2, "sweeping note: and both stretch averages accept it ("
                                            + std::to_string (k.stretchObservations) + ")");
            ok (k.stretchSpreadHz <= p.toleranceHz, "sweeping note: the STRETCH spread is inside the tolerance ("
                                                    + hz (k.stretchSpreadHz) + " Hz) — it cannot see the sweep");
            ok (k.frameSpreadHz > p.toleranceHz, "sweeping note: the FRAME spread is " + hz (k.frameSpreadHz)
                                                 + " Hz, outside it — this is the check that refuses the sweep");
            ok (! k.stationary && r.mains == HumMains::None, "sweeping note: so mains reads none");
            ok (k.frameObservations >= 4, "sweeping note: over " + std::to_string (k.frameObservations)
                                          + " frame observations");
        }
    }

    //---------------------------------------------------------------------------------------------------
    // 9. THE VALIDITY LADDER. Every rung is a NAMED incompleteness; none of them may read as "no hum".
    {
        HumDetectorParams p;
        // no quiet stretch at all: loud from end to end
        {
            HumDetector d; Sig s (1, kLen);
            s.noise (0, kLen, 0.3, 51);
            s.tone (0, kLen, 50.1, 1.0e-3);                          // the hum IS there, and unreachable
            const HumReport r = once (s, p, d);
            ok (! r.valid && r.reason == HumReason::NoQuietStretch, "ladder: nothing quiet -> NoQuietStretch");
            ok (r.mains == HumMains::None && r.line.hz == 0.0, "ladder: and no conclusion is published");
            ok (r.frames > 0 && r.quietFrames == 0, "ladder: the counts still say why (" + std::to_string (r.frames)
                                                    + " frames, 0 quiet)");
        }
        // exactly one quiet stretch, with a textbook hum in it: still not an answer
        {
            HumDetector d; Sig s (1, kLen);
            s.noise (0, kQuietA1, 3.0e-4, 52);
            s.noise (kQuietA1, kLen, 0.3, 53);
            s.tone (0, kQuietA1, 50.1, 1.0e-3);
            const HumReport r = once (s, p, d);
            ok (! r.valid && r.reason == HumReason::SingleQuietStretch,
                "ladder: one quiet stretch -> SingleQuietStretch, not a verdict");
            ok (r.quietStretches == 1, "ladder: and it says there was one");
            const auto k = d.candidate (0, 0);
            ok (k.baseFound && k.stretchObservations == 1,
                "ladder: the provisional candidate is still published as evidence");
            ok (k.maxIntraStretchSpreadHz >= 0.0 && k.frameObservations >= 2,
                "ladder: with the intra-stretch stillness a consumer can read ("
                + hz (k.maxIntraStretchSpreadHz) + " Hz over " + std::to_string (k.frameObservations) + " frames)");
        }
        // shorter than one window
        {
            HumDetector d; Sig s (1, kN - 1);
            s.noise (0, kN - 1, 3.0e-4, 54);
            const HumReport r = once (s, p, d);
            ok (! r.valid && r.reason == HumReason::ShorterThanWindow, "ladder: shorter than the window");
            ok (r.frames == 0 && r.tailUncoveredSamples == (std::int64_t) kN - 1,
                "ladder: no frame, and the whole programme is named uncovered");
        }
        // the empty stream, and T = 0 / 1 / W-1 / W / W+H-1 / W+H
        {
            for (std::size_t t : { (std::size_t) 0, (std::size_t) 1, kN - 1, kN, kN + kHop - 1, kN + kHop })
            {
                HumDetector d; Sig s (1, std::max<std::size_t> (t, 1));
                s.noise (0, t, 3.0e-4, 55);
                if (! d.prepare (kSr, 512, 1)) { ok (false, "ladder: prepare"); continue; }
                const auto v = s.ptrs (0);
                ok (d.process (v.data(), 1, (int) t), "ladder T=" + std::to_string (t) + ": the call is accepted");
                d.finish();
                const HumReport r = d.report (0);
                const std::int64_t want = t < kN ? 0 : (std::int64_t) ((t - kN) / kHop + 1);
                ok (r.frames == want, "ladder T=" + std::to_string (t) + ": " + std::to_string (r.frames)
                                      + " frames, oracle says " + std::to_string (want));
                ok (r.totalSamples == (std::int64_t) t, "ladder T=" + std::to_string (t) + ": the clock counted every sample");
                if (want == 0) ok (r.reason == HumReason::ShorterThanWindow,
                                   "ladder T=" + std::to_string (t) + ": and says ShorterThanWindow");
            }
        }
        // a coarser window than the pair needs
        {
            HumDetectorParams cp = p; cp.fftOrder = 11;              // bin 1.95 Hz: 49.0 and 50.0 merge
            HumDetector d; const auto s = twoStretch (1, 50.1, 1.0e-3, true);
            d.setParams (cp);
            ok (d.prepare (kSr, 512, 1), "ladder: a coarse order still prepares");
            const auto v = s.ptrs (0);
            ok (d.process (v.data(), 1, (int) s.len()), "ladder: and still processes");
            d.finish();
            const HumReport r = d.report (0);
            ok (! r.valid && r.reason == HumReason::InsufficientResolution,
                "ladder: bin " + hz (r.binHz) + " Hz -> InsufficientResolution, never a guessed 50");
            ok (r.mains == HumMains::None, "ladder: and no mains is claimed at that resolution");
        }
        // every frame holed
        {
            HumDetector d; Sig s (1, kLen);
            for (auto& v : s.ch[0]) v = std::numeric_limits<float>::quiet_NaN();
            const HumReport r = once (s, p, d);
            ok (! r.valid && r.reason == HumReason::AllFramesHoled, "ladder: all frames holed");
            ok (r.holedFrames == r.frames && r.finiteFrames == 0, "ladder: and every frame is counted as holed");
            ok (r.quietFrames == 0, "ladder: a holed frame's zeroed spectrum is NOT read as quiet");
        }
        // a channel that disappears: law 11a holes for the absent one
        {
            HumDetector d; const auto s = twoStretch (2, 50.14, 1.0e-3, true);
            ok (d.prepare (kSr, 4096, 2), "ladder: prepare two channels");
            const auto v = s.ptrs (0);
            ok (d.process (v.data(), 1, (int) s.len()), "ladder: a one-channel call into a two-channel object");
            d.finish();
            ok (d.report (0).mains == HumMains::Hz50, "ladder: the fed channel is measured");
            const HumReport r1 = d.report (1);
            ok (! r1.valid && r1.reason == HumReason::AllFramesHoled,
                "ladder: the absent channel is holes, not silence");
        }
    }

    //---------------------------------------------------------------------------------------------------
    // 10. A DC OFFSET IS NOT LOUDNESS, and the candidate bands are out of the quiet gate. Both are what let
    //     the instrument look at material it would otherwise call "not quiet".
    {
        HumDetectorParams p;
        {
            HumDetector d; auto s = twoStretch (1, 50.14, 1.0e-3, true);
            s.dc (kQuietA0, kQuietA1, 0.02);                         // -34 dBFS of DC: 15 dB over the gate
            s.dc (kQuietB0, kQuietB1, -0.02);
            const HumReport r = once (s, p, d);
            ok (r.valid && r.mains == HumMains::Hz50, "DC: a large offset does not censor the quiet stretches");
        }
        {
            // a hum LOUDER than the quiet threshold: without the candidate-band exclusion every frame would
            // be "not quiet" and the answer would be NoQuietStretch on exactly the files that have the fault
            HumDetector d; Sig s (1, kLen);
            s.noise (kQuietA0, kQuietA1, 3.0e-4, 61);
            s.noise (kQuietB0, kQuietB1, 3.0e-4, 62);
            s.noise (kLoud0, kLoud1, 0.3, 63);
            s.tone (kQuietA0, kQuietA1, 50.14, 0.02);                // -37 dBFS, over the -50 dBFS gate
            s.tone (kQuietB0, kQuietB1, 50.14, 0.02, 0.8);
            const HumReport r = once (s, p, d);
            ok (r.valid && r.mains == HumMains::Hz50, "loud hum: a -37 dBFS line does not hide behind the gate");
            ok (r.line.prominenceDb > 40.0, "loud hum: " + hz (r.line.prominenceDb) + " dB over the background");
        }
    }

    //---------------------------------------------------------------------------------------------------
    // 10b. THE QUIET MEASURE, CALIBRATED AGAINST ARITHMETIC — and the threshold it is compared with.
    //      The mutation stand found this hole: with every fixture's quiet bed 25 dB under the gate and its
    //      programme 35 dB over, deleting the interior-bin fold's factor 2 (a 3 dB error in the measure)
    //      left the suite GREEN. A measure needs an ABSOLUTE witness, and one fixture standing on the gate.
    {
        // a 0.5-amplitude sine at 700 Hz — outside every candidate band (h*50 <= 400, h*60 <= 480), so
        // nothing is subtracted and the Hann-weighted mean square must be exactly A^2/2 = 0.125
        HumDetectorParams p; p.traceCapacity = 64;
        HumDetector d; d.setParams (p);
        Sig s (1, kN + kHop);
        s.tone (0, s.len(), 700.0, 0.5);
        ok (d.prepare (kSr, 4096, 1), "quiet measure: prepare");
        const auto v = s.ptrs (0);
        ok (d.process (v.data(), 1, (int) s.len()), "quiet measure: process");
        d.finish();
        ok (d.storedTraceCount() >= 1, "quiet measure: a frame was traced");
        approx (d.traceAt (0).meanSquare, 0.125, 1.0e-5,
                "quiet measure: a 0.5 sine reads A^2/2 = 0.125 exactly (" + hz (d.traceAt (0).meanSquare) + ")");
        ok (d.traceAt (0).quiet == 0, "quiet measure: and -9 dBFS is not quiet");

        // ... and a bed that STRADDLES the gate: ms = A^2/2 = 3.125e-5 against a -48 dBFS gate (1.585e-5).
        // True: loud. At half the measure (1.5625e-5) it would fall under the gate and be selected.
        HumDetectorParams q; q.traceCapacity = 64; q.quietThresholdDb = -48.0;
        HumDetector e; e.setParams (q);
        Sig t (1, kN + kHop);
        t.tone (0, t.len(), 700.0, 7.905694150420949e-3);          // A^2/2 = 3.125e-5, 1.97x the gate
        ok (e.prepare (kSr, 4096, 1), "gate straddle: prepare");
        const auto vt = t.ptrs (0);
        ok (e.process (vt.data(), 1, (int) t.len()), "gate straddle: process");
        e.finish();
        approx (e.traceAt (0).meanSquare, 3.125e-5, 1.0e-9, "gate straddle: the measure is 1.97x the gate");
        ok (e.report (0).quietFrames == 0, "gate straddle: so NO frame is selected ("
                                           + std::to_string (e.report (0).quietFrames) + ")");
        ok (e.report (0).reason == HumReason::NoQuietStretch, "gate straddle: and the reason says so");
    }

    //---------------------------------------------------------------------------------------------------
    // 10bb. THE CANDIDATE BANDS ARE EXCLUDED AT EVERY HARMONIC, not just the low ones. The mutation stand
    //       found this hole: building the exclusion list WITHOUT merging leaves it unsorted (50*7 = 350
    //       arrives after 60*6 = 360), the quiet gate's ascending gap walk then skips those ranges, and no
    //       fixture noticed because none carried energy at 350 / 400 / 450 Hz. This one is made of nothing
    //       else: every tone sits in a candidate band, so a correct gate measures ~0 and selects the frame.
    {
        HumDetectorParams p; p.traceCapacity = 64;
        HumDetector d; d.setParams (p);
        Sig s (1, kN + 3 * kHop);
        for (double f : { 50.0, 120.0, 250.0, 300.0, 350.0, 400.0, 420.0, 480.0 })
            s.tone (0, s.len(), f, 0.01);                          // 8 tones, each ~ -43 dBFS: far over the gate
        ok (d.prepare (kSr, 4096, 1), "band exclusion: prepare");
        const auto v = s.ptrs (0);
        ok (d.process (v.data(), 1, (int) s.len()), "band exclusion: process");
        d.finish();
        const HumReport r = d.report (0);
        ok (r.frames >= 3, "band exclusion: frames were produced (" + std::to_string (r.frames) + ")");
        ok (r.quietFrames == r.frames, "band exclusion: EVERY frame is quiet — all eight tones are in candidate "
                                       "bands (" + std::to_string (r.quietFrames) + "/" + std::to_string (r.frames) + ")");
        ok (d.traceAt (0).meanSquare < 1.0e-7, "band exclusion: the measure reads "
                                               + hz (d.traceAt (0).meanSquare) + ", i.e. essentially nothing");
    }

    //---------------------------------------------------------------------------------------------------
    // 10bc. THE COMB IS EXACT MULTIPLES OF THE MEASURED FUNDAMENTAL, not of the nominal. At 50.45 Hz — just
    //       inside the 0.5 Hz tolerance — the 8th harmonic sits 3.6 Hz off the nominal grid, nearly four
    //       times pass B's 2-bin window, so searching around 8*50 instead of 8*f0est loses it. (The mutation
    //       stand found this hole: the other fixtures are only 0.14 Hz off and reach h = 4, where the
    //       nominal window still happens to contain the harmonic.)
    {
        HumDetectorParams p;
        HumDetector d;
        Sig s (1, kLen);
        s.noise (kQuietA0, kQuietA1, 3.0e-4, 121);
        s.noise (kQuietB0, kQuietB1, 3.0e-4, 122);
        s.noise (kLoud0, kLoud1, 0.3, 123);
        for (int h = 1; h <= 8; ++h)
        {
            const double f = 50.45 * (double) h, a = 1.4e-3 / (double) h;
            s.tone (kQuietA0, kQuietA1, f, a, 0.13 * (double) h);
            s.tone (kQuietB0, kQuietB1, f, a, 0.41 * (double) h);
        }
        const HumReport r = once (s, p, d);
        ok (r.valid && r.mains == HumMains::Hz50, "exact comb: a 50.45 Hz line is inside the tolerance and found");
        approx (r.line.hz, 50.45, 0.01, "exact comb: measured at " + hz (r.line.hz) + " Hz");
        ok (r.harmonicsObserved == 8, "exact comb: all eight harmonics were found ("
                                      + std::to_string (r.harmonicsObserved) + "/8)");
        approx (d.harmonic (0, 8).peak.hz, 403.6, 0.05, "exact comb: h = 8 sits at 8*50.45 = 403.6 Hz, which is "
                                                        "3.6 Hz off the nominal grid");
        ok (d.harmonic (0, 7).peak.accepted && d.harmonic (0, 8).peak.accepted,
            "exact comb: harmonics 7 and 8 are accepted");
        ok (d.harmonic (0, 8).inBand, "exact comb: and h = 8 has a full background window inside the band");
    }

    //---------------------------------------------------------------------------------------------------
    // 10c. THE LOCAL BACKGROUND, against a SECOND implementation written here. The mutation stand found
    //      this hole too: on a flat noise bed the paired geometric mean and a one-sided bin agree, so
    //      replacing sqrt(P[kp-d])*sqrt(P[kp+d]) with the high side alone left the suite green. The floor
    //      is now recomputed from the PUBLISHED band sums by an independent routine and compared bit-exactly,
    //      and separately measured on a STEEPLY TILTED background, where the two definitions diverge.
    {
        HumDetectorParams p;
        HumDetector d;
        Sig s (1, kLen);
        s.noise (kQuietA0, kQuietA1, 3.0e-4, 111);
        s.noise (kQuietB0, kQuietB1, 3.0e-4, 112);
        s.noise (kLoud0, kLoud1, 0.3, 113);
        // a steep low-frequency tilt through the 50 Hz region: three cascaded poles at 12 Hz put the
        // background's power on a 1/f^6 slope there, so the two flanks of the floor window differ by ~11 dB
        s.redNoise (kQuietA0, kQuietA1, 1.5e-3, 12.0, 3, 114);
        s.redNoise (kQuietB0, kQuietB1, 1.5e-3, 12.0, 3, 115);
        s.tone (kQuietA0, kQuietA1, 50.14, 1.2e-3);
        s.tone (kQuietB0, kQuietB1, 50.14, 1.2e-3, 0.6);
        const HumReport r = once (s, p, d);
        ok (r.valid && r.mains == HumMains::Hz50, "tilted floor: the line is still found on a steep LF slope");
        // the second implementation of the floor, over the published band sums
        const auto g = d.geometry();
        const double* band = d.bandPowerSum (0);
        const int kp = r.line.bin;
        const bool usable = band != nullptr && kp - g.floorSpanBins >= g.bandLo && kp + g.floorSpanBins <= g.bandHi;
        ok (usable, "tilted floor: the peak bin " + std::to_string (kp) + " has a full background window");
        if (usable)
        {
            std::vector<double> pairs;
            for (int dd = g.floorExcludeBins + 1; dd <= g.floorSpanBins; ++dd)
            {
                const double a = band[kp - dd - g.bandLo], b = band[kp + dd - g.bandLo];
                if (! (a > 0.0) || ! (b > 0.0) || ! std::isfinite (a) || ! std::isfinite (b)) continue;
                pairs.push_back (std::sqrt (a) * std::sqrt (b));
            }
            ok (pairs.size() >= (std::size_t) HumDetector::kMinFloorPairs,
                "tilted floor: the oracle found " + std::to_string (pairs.size()) + " usable pairs");
            if (! pairs.empty())
            {
                std::sort (pairs.begin(), pairs.end());
                const double want = pairs[(pairs.size() - 1) / 2] * (1.0 / (double) d.selectedFrames (0));
                ok (std::bit_cast<std::uint64_t> (want) == std::bit_cast<std::uint64_t> (r.line.floorPower),
                    "tilted floor: an independent paired-median implementation agrees BIT-FOR-BIT");
            }
            // the tilt is real: the flanks of the floor window differ by a lot, which is what makes the
            // pairing matter rather than being an ornament
            const double left = band[kp - g.floorSpanBins - g.bandLo], right = band[kp + g.floorSpanBins - g.bandLo];
            ok (left > right * 4.0, "tilted floor: the flanks differ by " + hz (10.0 * std::log10 (left / right))
                                    + " dB, so a one-sided floor would be wrong");
        }
    }

    //---------------------------------------------------------------------------------------------------
    // 10d. ONE PERIODOGRAM IS NOT A LINE. A stretch holding a single frame has no averaging, so a noise
    //      excursion 10 dB over its neighbours is common — and if both stretches carry the SAME noise (a
    //      looped passage, a duplicated silent section) the excursion repeats and passes stationarity with
    //      it. Measured at 4 kHz with the rule off: 75 of 1200 noise-only files, and 68 of 1200 dither-only
    //      files, were reported as mains. The rule takes both to 0, and names its own price.
    {
        HumDetectorParams p;
        // two quiet stretches of exactly ONE frame each, carrying IDENTICAL noise
        const std::size_t q = kN, gap = 4 * kHop;
        Sig s (1, q + gap + q);
        s.noise (0, q, 3.0e-4, 131);
        s.noise (q, q + gap, 0.3, 132);
        s.noise (q + gap, s.len(), 3.0e-4, 131);               // the same seed: the same samples
        HumDetector d;
        const HumReport r = once (s, p, d);
        ok (r.quietStretches == 2, "one-frame: there ARE two quiet stretches (" + std::to_string (r.quietStretches) + ")");
        ok (r.eligibleStretches == 0, "one-frame: and neither is eligible (" + std::to_string (r.eligibleStretches) + ")");
        ok (! r.valid && r.reason == HumReason::StretchesTooShort,
            "one-frame: so the answer is StretchesTooShort, NOT a clean report");
        ok (r.mains == HumMains::None, "one-frame: and no mains is claimed");
        // the SAME material with the rule relaxed to one frame is where the false positives lived
        HumDetectorParams relaxed = p; relaxed.minFramesPerObservation = 1;
        HumDetector e;
        const HumReport re = once (s, relaxed, e);
        ok (re.valid && re.reason == HumReason::Ok,
            "one-frame: with the rule off the same file becomes a verdict instead of an incompleteness");
        ok (re.eligibleStretches == 2, "one-frame: both stretches are then eligible");
    }

    //---------------------------------------------------------------------------------------------------
    // 10e. NOISE ALONE IS NOT HUM — over many seeds, which is the fixture a single-seed suite cannot be.
    {
        HumDetectorParams p;
        int reported = 0, valid = 0;
        for (unsigned seed = 0; seed < 40; ++seed)
        {
            Sig s (1, kLen);
            s.noise (kQuietA0, kQuietA1, 3.0e-4, 900u + seed);
            s.noise (kQuietB0, kQuietB1, 3.0e-4, 1900u + seed);
            s.noise (kLoud0, kLoud1, 0.3, 2900u + seed);
            HumDetector d;
            const HumReport r = once (s, p, d);
            if (r.valid) ++valid;
            if (r.mains != HumMains::None) ++reported;
        }
        ok (valid == 40, "noise only: all 40 seeds produced a VALID report (" + std::to_string (valid) + "/40)");
        ok (reported == 0, "noise only: and none of them reported mains (" + std::to_string (reported) + "/40)");
        // amplitude-modulated noise has 50 Hz SIDEBANDS and no 50 Hz line; it must not read as one either
        int am = 0;
        for (unsigned seed = 0; seed < 12; ++seed)
        {
            Sig s (1, kLen);
            s.noise (kQuietA0, kQuietA1, 3.0e-4, 700u + seed);
            s.noise (kQuietB0, kQuietB1, 3.0e-4, 800u + seed);
            s.noise (kLoud0, kLoud1, 0.3, 850u + seed);
            for (std::size_t i = 0; i < s.len(); ++i)
                s.ch[0][i] = (float) ((double) s.ch[0][i] * (1.0 + std::cos (2.0 * core::kPi * 50.0 * (double) i / kSr)));
            HumDetector d;
            if (once (s, p, d).mains != HumMains::None) ++am;
        }
        ok (am == 0, "AM noise: 50 Hz modulation is not a 50 Hz line (" + std::to_string (am) + "/12 reported)");
    }

    //---------------------------------------------------------------------------------------------------
    // 10f. THE FLOOR'S TIE RULE, at an EVEN population. The default geometry gives 17 pairs, an odd count,
    //      where the lower median and the upper median are the same element — so `(used-1)/2` and `used/2`
    //      are indistinguishable. floorExcludeHz = 2.0 gives 16 pairs, where they are not.
    {
        HumDetectorParams p; p.floorExcludeHz = 2.0;
        const auto g = HumDetector::geometryFor (kSr, p);
        ok (g.ok && g.floorPairs == 22, "even floor: the population is " + std::to_string (g.floorPairs)
                                        + " pairs, an even count");
        HumDetector d;
        const auto s = twoStretch (1, 50.14, 1.2e-3, true, 141);
        const HumReport r = once (s, p, d);
        ok (r.valid && r.mains == HumMains::Hz50, "even floor: the line is found");
        const double* band = d.bandPowerSum (0);
        const int kp = r.line.bin;
        const bool room = band != nullptr && kp - g.floorSpanBins >= g.bandLo && kp + g.floorSpanBins <= g.bandHi;
        ok (room, "even floor: the background window fits");
        if (room)
        {
            std::vector<double> pairs;
            for (int dd = g.floorExcludeBins + 1; dd <= g.floorSpanBins; ++dd)
            {
                const double a = band[kp - dd - g.bandLo], b = band[kp + dd - g.bandLo];
                if (! (a > 0.0) || ! (b > 0.0)) continue;
                pairs.push_back (std::sqrt (a) * std::sqrt (b));
            }
            std::sort (pairs.begin(), pairs.end());
            ok (pairs.size() == 22, "even floor: the oracle sees 22 pairs too (" + std::to_string (pairs.size()) + ")");
            const double lower = pairs[(pairs.size() - 1) / 2] * (1.0 / (double) d.selectedFrames (0));
            const double upper = pairs[pairs.size() / 2] * (1.0 / (double) d.selectedFrames (0));
            ok (std::bit_cast<std::uint64_t> (lower) != std::bit_cast<std::uint64_t> (upper),
                "even floor: the two medians really are different elements here");
            ok (std::bit_cast<std::uint64_t> (lower) == std::bit_cast<std::uint64_t> (r.line.floorPower),
                "even floor: and the LOWER one is what the instrument published");
        }
    }

    //---------------------------------------------------------------------------------------------------
    // 10g. PARAMETERS ARE COMMITTED BY prepare(), NOT READ LIVE. Writing setParams() under a live object
    //      used to move a vector stride: maxHarmonic 1 -> 64 then read entry 127 of a two-entry vector
    //      (ASan-confirmed by the review round). And a geometry that cannot serve the measurement is now
    //      REFUSED rather than accepted into a report that reads clean.
    {
        HumDetectorParams p;
        HumDetector d; d.setParams (p);
        ok (d.prepare (kSr, 4096, 1), "committed: prepare with the defaults");
        HumDetectorParams wild = p; wild.maxHarmonic = 64; wild.maxStretches = 1 << 16; wild.toleranceHz = 40.0;
        d.setParams (wild);                                    // pending only
        ok (d.committedParams().maxHarmonic == p.maxHarmonic, "committed: the live object kept its own harmonics");
        ok (d.params().maxHarmonic == 64, "committed: while the pending set carries the new value");
        ok (d.harmonic (0, 0, 64).index == 0, "committed: a harmonic past the committed count answers empty");
        const auto s = twoStretch (1, 50.14, 1.0e-3, true, 151);
        const auto v = s.ptrs (0);
        ok (d.process (v.data(), 1, (int) s.len()), "committed: and the run is unaffected");
        d.finish();
        ok (d.report (0).mains == HumMains::Hz50, "committed: the measurement used the committed tolerance");

        // maxHarmonic = 1 must still be able to see the 100/120 Hz hypotheses: the band is sized for them
        // whatever the REPORT asks for. It used to cut the band at 72 Hz and answer valid/None on a 120 comb.
        HumDetectorParams one = p; one.maxHarmonic = 1;
        const auto g1 = HumDetector::geometryFor (kSr, one);
        ok (g1.ok && g1.spanHarmonics == 2, "committed: maxHarmonic 1 still spans 2 base harmonics");
        ok (g1.stretchHi > (int) (130.0 / kBin), "committed: and the active band reaches past 120 Hz (bin "
                                                 + std::to_string (g1.stretchHi) + ")");
        HumDetector d1;
        Sig s120 (1, kLen);
        s120.noise (kQuietA0, kQuietA1, 3.0e-4, 161);
        s120.noise (kQuietB0, kQuietB1, 3.0e-4, 162);
        s120.noise (kLoud0, kLoud1, 0.3, 163);
        s120.tone (kQuietA0, kQuietA1, 120.0, 1.0e-3);
        s120.tone (kQuietB0, kQuietB1, 120.0, 1.0e-3, 0.5);
        const HumReport r120 = once (s120, one, d1);
        ok (r120.valid && r120.mains == HumMains::Hz60,
            "committed: a 120 Hz comb is found with maxHarmonic = 1 (base " + std::to_string (r120.baseHarmonic) + ")");

        // a background window too small to serve kMinFloorPairs is refused, not accepted into a clean answer
        HumDetectorParams thin = p; thin.floorSpanHz = 2.0;
        ok (! HumDetector::geometryFor (kSr, thin).ok, "committed: floorSpanHz 2.0 leaves 1 pair -> refused");
        // and an absurd Hz parameter is refused rather than converted (the UB cast the review found)
        HumDetectorParams huge = p; huge.searchHz = 1.0e20;
        ok (! HumDetector::geometryFor (kSr, huge).ok, "committed: searchHz 1e20 is refused, not cast");
        huge = p; huge.floorSpanHz = 1.0e20;
        ok (! HumDetector::geometryFor (kSr, huge).ok, "committed: floorSpanHz 1e20 is refused");
        huge = p; huge.toleranceHz = 1.0e20;
        ok (! HumDetector::geometryFor (kSr, huge).ok, "committed: toleranceHz 1e20 is refused");
        huge = p; huge.minFramesPerObservation = 0;
        ok (! HumDetector::geometryFor (kSr, huge).ok, "committed: minFramesPerObservation 0 is refused");
        // pass B's enumeration follows its configured tolerance rather than searchHz
        HumDetectorParams wideH = p; wideH.harmonicToleranceBins = 12.0;
        const auto gw = HumDetector::geometryFor (kSr, wideH);
        ok (gw.ok, "committed: a wide harmonic tolerance is a legal configuration");
    }

    //---------------------------------------------------------------------------------------------------
    // 10h. THE RESOLUTION BOUND AT ANOTHER RATE. A bound of "bin <= 0.5 Hz" admits 0.48828 Hz at 4, 8, 16,
    //      32 and 64 kHz, where 1.0 Hz is 2.048 bins — the exact separation at which two equal Hann main
    //      lobes stop having a dip. Measured there: a 49.0 Hz note over a 50.0 Hz hum 14 dB below it was
    //      MISSED, and an equal-amplitude pair merged into one accepted peak at 49.66 Hz. The bound is a
    //      separation in bins now, so 32 kHz takes order 17 and the same case is found.
    {
        const double sr = 32000.0;
        HumDetectorParams p;
        const auto g = HumDetector::geometryFor (sr, p);
        ok (g.ok && g.order == 17, "32 kHz: AUTO takes order 17, not 16 (got " + std::to_string (g.order) + ")");
        ok (HumDetector::kMinNoteSeparationHz / g.binHz > HumDetector::kMinSeparationBins,
            "32 kHz: the pair is " + hz (HumDetector::kMinNoteSeparationHz / g.binHz) + " bins apart");
        // order 16 would be admitted by a 0.5 Hz bound and is refused by this one
        HumDetectorParams sixteen = p; sixteen.fftOrder = 16;
        const auto g16 = HumDetector::geometryFor (sr, sixteen);
        ok (g16.ok && g16.binHz > HumDetector::kMaxBinHz,
            "32 kHz: order 16 gives bin " + hz (g16.binHz) + " Hz, past the bound");

        // and the masked-hum case really is found at the chosen order
        const std::int64_t n = g.n, hop = g.hop;
        const std::size_t q = (std::size_t) (n + hop), gap = (std::size_t) (2 * hop);
        const std::size_t len = 2 * q + gap;
        std::vector<float> x (len, 0.0f);
        std::mt19937 rng (171);
        std::uniform_real_distribution<double> u (-1.0, 1.0);
        for (std::size_t i = 0; i < len; ++i)
        {
            const bool loud = i >= q && i < q + gap;
            double v = (loud ? 0.3 : 3.0e-4) * u (rng);
            if (! loud)
            {
                const double t = (double) i / sr;
                v += 2.0e-3 * std::sin (2.0 * core::kPi * 49.0 * t);       // the note, 14 dB louder
                v += 4.0e-4 * std::sin (2.0 * core::kPi * 50.0 * t + 0.4); // the hum, 1.0 Hz away
            }
            x[i] = (float) v;
        }
        HumDetector d; d.setParams (p);
        ok (d.prepare (sr, 8192, 1), "32 kHz: prepare");
        const float* pv[1] = { x.data() };
        ok (d.process (pv, 1, (int) len), "32 kHz: process");
        d.finish();
        const HumReport r = d.report (0);
        ok (r.valid && r.mains == HumMains::Hz50,
            "32 kHz: a 50.0 Hz hum under a 49.0 Hz note 14 dB louder is FOUND (reason "
            + std::to_string ((int) r.reason) + ")");
        approx (r.line.hz, 50.0, 0.06, "32 kHz: at " + hz (r.line.hz) + " Hz");
    }

    //---------------------------------------------------------------------------------------------------
    // 10i. THE GATES THE MUTATION STAND AND THE ADVERSARIAL ROUND FOUND UNPINNED.
    {
        HumDetectorParams p;
        // (a) DRIFT BEYOND THE TOLERANCE — the stretch-spread half of the stationarity gate. 49.72 Hz in one
        //     stretch, 50.28 in the other: both inside the tolerance, spread 0.56 Hz outside it.
        {
            HumDetector d;
            Sig s (1, kLen);
            s.noise (kQuietA0, kQuietA1, 3.0e-4, 181);
            s.noise (kQuietB0, kQuietB1, 3.0e-4, 182);
            s.noise (kLoud0, kLoud1, 0.3, 183);
            for (int h = 1; h <= 3; ++h)
            {
                s.tone (kQuietA0, kQuietA1, 49.72 * (double) h, 1.2e-3 / (double) h);
                s.tone (kQuietB0, kQuietB1, 50.28 * (double) h, 1.2e-3 / (double) h, 0.5);
            }
            const HumReport r = once (s, p, d);
            const auto k = d.candidate (0, 0);
            ok (k.baseFound && k.stretchObservations == 2,
                "drift: both stretches show a 50 Hz-compatible line (" + std::to_string (k.stretchObservations) + ")");
            ok (k.stretchSpreadHz > p.toleranceHz, "drift: the STRETCH spread is " + hz (k.stretchSpreadHz)
                                                   + " Hz, past the tolerance — this is the gate");
            ok (! k.stationary, "drift: so it is not stationary");
            ok (! r.valid && r.reason == HumReason::CandidateNotStationary,
                "drift: and hum that moved 0.56 Hz does NOT read as clean");
        }
        // (b) THE PROMINENCE GATE'S VALUE. A line 7 dB over its background must be refused: nothing else in
        //     the suite has a line between 5 and 10 dB, so halving the gate used to go unnoticed.
        {
            HumDetector d;
            Sig s (1, kLen);
            // scaled so the LEVEL is comfortably over minLevelDbfs: the prominence gate must be what
            // refuses this, not the absolute floor, or halving the gate would go unnoticed again
            s.noise (kQuietA0, kQuietA1, 3.0e-3, 191);
            s.noise (kQuietB0, kQuietB1, 3.0e-3, 192);
            s.noise (kLoud0, kLoud1, 0.3, 193);
            s.tone (kQuietA0, kQuietA1, 50.14, 1.05e-4);
            s.tone (kQuietB0, kQuietB1, 50.14, 1.05e-4, 0.6);
            const HumReport r = once (s, p, d);
            const auto k = d.candidate (0, 0);
            const double prom = k.windowPeak.prominenceDb;
            ok (prom > 3.0 && prom < p.minProminenceDb,
                "weak line: the window peak stands " + hz (prom) + " dB up, under the 10 dB gate");
            ok (k.windowPeak.tonePower > 1.0e-10,
                "weak line: and its LEVEL is over the absolute floor, so the prominence gate is what refused it");
            ok (! k.baseFound && r.mains == HumMains::None,
                "weak line: so no base is accepted and no mains is claimed");
        }
        // (c) TWO CANDIDATES BOTH PASSING — the winner rule. A 50 Hz comb and a stronger 60 Hz comb.
        {
            HumDetector d;
            Sig s (1, kLen);
            s.noise (kQuietA0, kQuietA1, 3.0e-4, 201);
            s.noise (kQuietB0, kQuietB1, 3.0e-4, 202);
            s.noise (kLoud0, kLoud1, 0.3, 203);
            for (int h = 1; h <= 3; ++h)
            {
                s.tone (kQuietA0, kQuietA1, 50.10 * h, 4.0e-4 / h);       s.tone (kQuietB0, kQuietB1, 50.10 * h, 4.0e-4 / h, 0.2);
                s.tone (kQuietA0, kQuietA1, 60.20 * h, 2.0e-3 / h, 0.7);  s.tone (kQuietB0, kQuietB1, 60.20 * h, 2.0e-3 / h, 0.9);
            }
            const HumReport r = once (s, p, d);
            ok (d.candidate (0, 0).passed && d.candidate (0, 1).passed, "dual comb: BOTH candidates pass");
            ok (r.valid && r.mains == HumMains::Hz60, "dual comb: the more prominent one wins (60 Hz)");
            approx (r.line.hz, 60.20, 0.02, "dual comb: reported at " + hz (r.line.hz) + " Hz");
            ok (d.candidate (0, 0).base.prominenceDb < d.candidate (0, 1).base.prominenceDb,
                "dual comb: and 50 Hz really was the weaker of the two");
        }
        // (d) TWO CHANNELS THAT DIFFER — the summary channel's rule. Hum on the right only.
        {
            HumDetector d;
            Sig s (2, kLen);
            s.noise (kQuietA0, kQuietA1, 3.0e-4, 211);
            s.noise (kQuietB0, kQuietB1, 3.0e-4, 212);
            s.noise (kLoud0, kLoud1, 0.3, 213);
            for (int h = 1; h <= 3; ++h)                               // added to channel 1 alone
                for (std::size_t i = 0; i < s.len(); ++i)
                    if (! (i >= kLoud0 && i < kLoud1))
                        s.ch[1][i] += (float) ((1.2e-3 / h) * std::sin (2.0 * core::kPi * 50.14 * h * (double) i / kSr));
            d.setParams (p);
            ok (d.prepare (kSr, 4096, 2), "two channels: prepare");
            const auto v = s.ptrs (0);
            ok (d.process (v.data(), 2, (int) s.len()), "two channels: process");
            d.finish();
            ok (d.report (1).mains == HumMains::Hz50, "two channels: the right channel carries the hum");
            ok (d.report (0).mains == HumMains::None, "two channels: the left does not");
            ok (d.strongestChannel() == 1, "two channels: and the summary speaks for the right one (got "
                                           + std::to_string (d.strongestChannel()) + ")");
            ok (d.report().mains == HumMains::Hz50, "two channels: so report() finds it");
        }
        // ... and with the hum on BOTH channels at different levels, which is what exercises the comparison
        // rather than the "only one candidate passed" shortcut
        {
            HumDetector d;
            Sig s (2, kLen);
            s.noise (kQuietA0, kQuietA1, 3.0e-4, 214);
            s.noise (kQuietB0, kQuietB1, 3.0e-4, 215);
            s.noise (kLoud0, kLoud1, 0.3, 216);
            for (int c = 0; c < 2; ++c)
            {
                const double a = c == 0 ? 3.0e-4 : 1.5e-3;         // the right channel is 14 dB louder
                for (std::size_t i = 0; i < s.len(); ++i)
                    if (! (i >= kLoud0 && i < kLoud1))
                        s.ch[(std::size_t) c][i] += (float) (a * std::sin (2.0 * core::kPi * 50.14 * (double) i / kSr));
            }
            d.setParams (p);
            ok (d.prepare (kSr, 4096, 2), "both channels: prepare");
            const auto v = s.ptrs (0);
            ok (d.process (v.data(), 2, (int) s.len()), "both channels: process");
            d.finish();
            ok (d.report (0).mains == HumMains::Hz50 && d.report (1).mains == HumMains::Hz50,
                "both channels: both carry the line");
            ok (d.report (1).line.prominenceDb > d.report (0).line.prominenceDb + 6.0,
                "both channels: the right one stands " + hz (d.report (1).line.prominenceDb - d.report (0).line.prominenceDb)
                + " dB further up");
            ok (d.strongestChannel() == 1, "both channels: so the summary picks the LOUDER one (got "
                                           + std::to_string (d.strongestChannel()) + ")");
        }
        // (e) DIGITAL BLACK IS NOT EVIDENCE OF ABSENCE. Exact zeros in both quiet passages, a hum in the
        //     programme: a selected all-zero frame would report "no hum" from frames that saw nothing.
        {
            HumDetector d;
            Sig s (1, kLen);
            s.noise (kLoud0, kLoud1, 0.3, 221);
            s.tone (kLoud0, kLoud1, 50.14, 1.0e-2);                    // the hum lives where it is masked
            const HumReport r = once (s, p, d);
            ok (r.silentFrames > 0, "digital black: " + std::to_string (r.silentFrames) + " frames are exactly silent");
            ok (r.quietFrames == 0, "digital black: and none of them is selected");
            ok (! r.valid && r.reason == HumReason::NoQuietStretch,
                "digital black: so the answer is NoQuietStretch, not a clean report");
            // dither, by contrast, IS a background a line could stand above
            HumDetector e;
            Sig t (1, kLen);
            t.noise (kQuietA0, kQuietA1, 3.0e-5, 231);
            t.noise (kQuietB0, kQuietB1, 3.0e-5, 232);
            t.noise (kLoud0, kLoud1, 0.3, 233);
            t.tone (kQuietA0, kQuietA1, 50.14, 2.0e-4);
            t.tone (kQuietB0, kQuietB1, 50.14, 2.0e-4, 0.4);
            const HumReport rt = once (t, p, e);
            ok (rt.silentFrames == 0 && rt.quietFrames == 6, "dither: every quiet frame is still selected");
            ok (rt.valid && rt.mains == HumMains::Hz50, "dither: and a -74 dBFS line is found over it");
        }
        // (f) THE TRACE'S STRETCH INDEX, absolutely. The 8a comparison is the code against itself, so an
        //     off-by-one in the published index survives it.
        {
            HumDetectorParams tp = p; tp.traceCapacity = 64;
            HumDetector d;
            const auto s = twoStretch (1, 50.14, 1.0e-3, true, 241);
            (void) once (s, tp, d);
            ok (d.traceCount() == 11, "trace: 11 frames, one record each (" + std::to_string (d.traceCount()) + ")");
            int first = -1, second = -1, unselected = 0;
            for (std::int64_t i = 0; i < d.storedTraceCount(); ++i)
            {
                const auto t = d.traceAt (i);
                if (t.quiet == 0) { ++unselected; ok (t.stretchIndex == -1, "trace: an unselected frame has no stretch"); }
                else if (first < 0) first = (int) t.stretchIndex;
                else if (t.stretchIndex != first && second < 0) second = (int) t.stretchIndex;
            }
            ok (first == 0 && second == 1, "trace: the two stretches are indexed 0 and 1 (got "
                                           + std::to_string (first) + ", " + std::to_string (second) + ")");
            ok (unselected == 5, "trace: and 5 frames straddle the programme (" + std::to_string (unselected) + ")");
            ok (d.traceAt (0).frameStart == 0 && d.traceAt (0).frameEnd == (std::int64_t) kN,
                "trace: the first record is the first window");
        }
        // (g) THE CAPACITY BOUNDARY: stretchCount == maxStretches is COMPLETE, one more is not.
        {
            HumDetectorParams two = p; two.maxStretches = 2;
            HumDetector d;
            const auto s = twoStretch (1, 50.14, 1.0e-3, true, 251);
            const HumReport r = once (s, two, d);
            ok (r.quietStretches == 2 && r.storedStretches == 2 && r.stretchesComplete,
                "capacity boundary: exactly at the capacity, the list is complete");
            ok (d.stretchesComplete (0), "capacity boundary: and the accessor agrees with the report");
            HumDetectorParams one = p; one.maxStretches = 1;
            HumDetector e;
            const HumReport re = once (s, one, e);
            ok (re.storedStretches == 1 && ! re.stretchesComplete, "capacity boundary: one less, and it is not");
            ok (! e.stretchesComplete (0), "capacity boundary: the accessor agrees there too");
            ok (e.stretch (0, 1).endSample == 0 && e.stretch (0, 1).frames == 0,
                "capacity boundary: past the stored prefix the accessor answers an EMPTY record");
        }
    }

    //---------------------------------------------------------------------------------------------------
    // 10j. A RATIO IS NOT EVIDENCE. A programme that is ONE pure tone has a spectrum of ~1e-23 away from
    //      that tone, and a local maximum of that residue stands 18-37 dB over the residue beside it — three
    //      such fixtures reported confident mains at -230 to -348 dBFS before the absolute floor existed.
    //      A pure tone's frames are identical, so the frames rule does not catch this; only the level does.
    {
        HumDetectorParams p;
        struct Case { const char* name; int kind; double v; };
        const Case cases[] = { { "a pure 440 Hz tone", 0, 440.0 }, { "a pure 700 Hz tone", 0, 700.0 },
                               { "a constant +1.0", 1, 1.0 }, { "a constant -0.5", 1, -0.5 } };
        for (const auto& cse : cases)
        {
            Sig s (1, kLen);
            if (cse.kind == 0) { s.tone (kQuietA0, kQuietA1, cse.v, 1.0e-3); s.tone (kQuietB0, kQuietB1, cse.v, 1.0e-3); }
            else               { s.dc (kQuietA0, kQuietA1, cse.v); s.dc (kQuietB0, kQuietB1, cse.v); }
            s.noise (kLoud0, kLoud1, 0.3, 261);
            HumDetector d;
            const HumReport r = once (s, p, d);
            ok (r.eligibleStretches == 2, std::string ("residue: ") + cse.name + " gives two eligible stretches");
            ok (r.mains == HumMains::None, std::string ("residue: ") + cse.name + " reports NO mains");
            ok (! d.candidate (0, 0).baseFound && ! d.candidate (0, 1).baseFound,
                std::string ("residue: ") + cse.name + " accepts no base at either nominal");
        }
        // a real line that is simply TOO QUIET to be one: 76 dB of prominence over a dither bed, at -123 dBFS
        {
            Sig s (1, kLen);
            s.noise (kQuietA0, kQuietA1, 1.0e-8, 271);
            s.noise (kQuietB0, kQuietB1, 1.0e-8, 272);
            s.noise (kLoud0, kLoud1, 0.3, 273);
            s.tone (kQuietA0, kQuietA1, 50.14, 1.0e-6);
            s.tone (kQuietB0, kQuietB1, 50.14, 1.0e-6, 0.4);
            HumDetector d;
            const HumReport r = once (s, p, d);
            const auto k = d.candidate (0, 0);
            ok (k.windowPeak.found && k.windowPeak.prominenceDb > 30.0,
                "level floor: the line stands " + hz (k.windowPeak.prominenceDb) + " dB over its background");
            ok (k.windowPeak.tonePower < 1.0e-10,
                "level floor: but its level is " + hz (10.0 * std::log10 (k.windowPeak.tonePower))
                + " dBFS, under the -100 dBFS floor");
            ok (! k.baseFound && r.mains == HumMains::None, "level floor: so it is not accepted");
            // and with the floor moved down, the very same file IS reported — the gate is what decided
            HumDetectorParams low = p; low.minLevelDbfs = -140.0;
            HumDetector e;
            ok (once (s, low, e).mains == HumMains::Hz50, "level floor: lowering the floor accepts it");
        }
    }

    //---------------------------------------------------------------------------------------------------
    // 10k. A COMB WHOSE BASE IS OUT OF SCOPE is not an absence. 150/300/450 Hz — a three-pulse rectifier's
    //      ripple — has no h = 1 or h = 2, so no `mains` can be claimed; two harmonics agreeing on one
    //      fundamental is still a comb, and it gets a named incompleteness instead of a clean report.
    {
        HumDetectorParams p;
        HumDetector d;
        Sig s (1, kLen);
        s.noise (kQuietA0, kQuietA1, 3.0e-4, 281);
        s.noise (kQuietB0, kQuietB1, 3.0e-4, 282);
        s.noise (kLoud0, kLoud1, 0.3, 283);
        for (int h : { 3, 6 })
        {
            const double a = 1.2e-3 * 3.0 / (double) h;
            s.tone (kQuietA0, kQuietA1, 50.0 * (double) h, a);
            s.tone (kQuietB0, kQuietB1, 50.0 * (double) h, a, 0.5);
        }
        const HumReport r = once (s, p, d);
        const auto k = d.candidate (0, 0);
        ok (! k.baseFound, "comb without base: no h = 1 or h = 2 was observed");
        ok (k.lowestHarmonicObserved == 3, "comb without base: the lowest one seen is h = "
                                           + std::to_string (k.lowestHarmonicObserved));
        ok (k.combWithoutBase, "comb without base: two harmonics agree on a fundamental");
        approx (k.combFundamentalHz, 50.0, 0.02, "comb without base: implied at " + hz (k.combFundamentalHz)
                                                 + " Hz — derived, never observed");
        ok (! r.valid && r.reason == HumReason::CombWithoutBase,
            "comb without base: so the answer is CombWithoutBase, not a clean report");
        ok (r.mains == HumMains::None, "comb without base: and no mains family is invented");
        // ONE line at a multiple is not a comb: 300 Hz alone must not trigger it
        HumDetector e;
        Sig t (1, kLen);
        t.noise (kQuietA0, kQuietA1, 3.0e-4, 291);
        t.noise (kQuietB0, kQuietB1, 3.0e-4, 292);
        t.noise (kLoud0, kLoud1, 0.3, 293);
        t.tone (kQuietA0, kQuietA1, 300.0, 1.2e-3);
        t.tone (kQuietB0, kQuietB1, 300.0, 1.2e-3, 0.5);
        const HumReport rt = once (t, p, e);
        ok (! e.candidate (0, 0).combWithoutBase && ! e.candidate (0, 1).combWithoutBase,
            "comb without base: a single line at 300 Hz is not a comb");
        ok (rt.valid && rt.mains == HumMains::None, "comb without base: and 300 Hz alone stays a valid None");
    }

    //---------------------------------------------------------------------------------------------------
    // 11. LAW 8a. The full trace, every accumulated bin and every report field, bit-exact by std::bit_cast
    //     across 14 slicings (n == 0 calls and calls larger than maxBlock included), seeded ragged
    //     partitions, and two prepared maxBlock values. Boundaries are placed exactly on the transitions:
    //     the hop, the window end, and the entry and exit of a quiet stretch.
    {
        HumDetectorParams p;
        p.traceCapacity = 4096;
        p.maxStretches = 64;
        const auto s = twoStretch (2, 50.14, 1.0e-3, true, 71);
        const std::vector<std::vector<int>> slicings = {
            { (int) kLen }, { 1 }, { 2 }, { 3 }, { 4093 }, { (int) kHop - 1 }, { (int) kHop }, { (int) kHop + 1 },
            { (int) kN - 1 }, { (int) kN }, { (int) kN + 1 }, { 512 }, { 65536 },
            { 0, 1, 0, 4095, 0, 700 },
            // exactly on the transitions: the end of stretch A, the end of the loud gap, the last hop
            { (int) kQuietA1, (int) (kLoud1 - kQuietA1), (int) (kLen - kLoud1) },
        };
        const Out ref = runSliced (s, p, slicings[0], 512);
        ok (ref.ok && ref.bits.size() > 2000, "8a: the trace is non-trivial (" + std::to_string (ref.bits.size()) + " words)");
        int differing = 0;
        for (std::size_t i = 1; i < slicings.size(); ++i)
            if (runSliced (s, p, slicings[i], 512) != ref) ++differing;
        ok (differing == 0, "8a: " + std::to_string (differing) + " of " + std::to_string (slicings.size() - 1)
                            + " slicings differ");
        int ragged = 0;
        for (unsigned seed = 0; seed < 6; ++seed)
        {
            std::mt19937 rng (seed + 101u);
            std::uniform_int_distribution<int> dist (1, 6000);
            std::vector<int> sl;
            for (int i = 0; i < 64; ++i) sl.push_back (dist (rng));
            if (runSliced (s, p, sl, 512) != ref) ++ragged;
        }
        ok (ragged == 0, "8a: " + std::to_string (ragged) + " of 6 seeded ragged partitions differ");
        ok (runSliced (s, p, { 512 }, 64) == ref && runSliced (s, p, { 512 }, 65536) == ref,
            "8a: maxBlock 64 and 65536 give the same report");
        // the fixture really does exercise selection: two stretches, six selected frames of eleven, both
        // channels — a constant tone selected end to end would compare bit-identically and prove much less
        {
            HumDetector probe; probe.setParams (p);
            ok (probe.prepare (kSr, 4096, 2), "8a: the witness prepare");
            const auto vp = s.ptrs (0);
            ok (probe.process (vp.data(), 2, (int) s.len()), "8a: the witness run");
            probe.finish();
            const HumReport pr = probe.report (0);
            ok (pr.frames == 11 && pr.quietFrames == 6 && pr.quietStretches == 2,
                "8a: the fixture closes 11 frames, selects 6 of them, in 2 stretches ("
                + std::to_string (pr.frames) + "/" + std::to_string (pr.quietFrames) + "/"
                + std::to_string (pr.quietStretches) + ")");
            ok (probe.traceCount() == 22, "8a: and traces every (frame, channel) pair ("
                                          + std::to_string (probe.traceCount()) + ")");
        }
    }

    //---------------------------------------------------------------------------------------------------
    // 12. CAPACITY IS DATA. A stretch list of 1 must classify exactly as a list of 64 — only the records
    //     and the completeness flag may differ. Same for the frame trace.
    {
        HumDetectorParams big; big.maxStretches = 64; big.traceCapacity = 4096;
        HumDetectorParams small = big; small.maxStretches = 1; small.traceCapacity = 3;
        const auto s = twoStretch (1, 50.14, 1.0e-3, true, 81);
        HumDetector a, b;
        a.setParams (big); b.setParams (small);
        const HumReport ra = once (s, big, a);
        const HumReport rb = once (s, small, b);
        ok (ra.valid && rb.valid, "capacity: both runs are valid");
        ok (ra.mains == rb.mains, "capacity: the same mains");
        ok (std::bit_cast<std::uint64_t> (ra.line.hz) == std::bit_cast<std::uint64_t> (rb.line.hz),
            "capacity: the peak position is BIT-identical");
        ok (std::bit_cast<std::uint64_t> (ra.line.prominenceDb) == std::bit_cast<std::uint64_t> (rb.line.prominenceDb),
            "capacity: so is the prominence");
        ok (ra.harmonicsObserved == rb.harmonicsObserved, "capacity: and the harmonic count");
        ok (ra.quietStretches == rb.quietStretches && rb.quietStretches == 2,
            "capacity: both counted every stretch (" + std::to_string (rb.quietStretches) + ")");
        ok (ra.storedStretches == 2 && rb.storedStretches == 1, "capacity: only the stored list is shorter");
        ok (ra.stretchesComplete && ! rb.stretchesComplete, "capacity: and the short one says it is incomplete");
        ok (b.stretch (0, 0).index == 0 && b.stretch (0, 1).startSample == 0,
            "capacity: the stored prefix is the first stretch, and past it the accessor is empty");
        ok (b.traceCount() > 3 && b.storedTraceCount() == 3 && ! b.traceComplete(),
            "capacity: the frame trace overflows the same way");
        const auto ka = a.candidate (0, 0), kb = b.candidate (0, 0);
        ok (ka.stretchObservations == kb.stretchObservations
            && std::bit_cast<std::uint64_t> (ka.frameSpreadHz) == std::bit_cast<std::uint64_t> (kb.frameSpreadHz),
            "capacity: the candidate summaries are identical — the 2nd stretch was measured, only not stored");
    }

    //---------------------------------------------------------------------------------------------------
    // 13. LIFECYCLE: finish idempotent, process refused after it, reset replays bit-identically, and
    //     process()/finish() touch no heap. The allocation counter is only strict on libc++.
    {
        HumDetectorParams p; p.traceCapacity = 256;
        const auto s = twoStretch (1, 50.14, 1.0e-3, true, 91);
        HumDetector d; d.setParams (p);
        const long long before = alloc::count.load();
        ok (d.prepare (kSr, 4096, 1), "lifecycle: prepare");
        ok (alloc::count.load() > before, "lifecycle: prepare is where the heap is touched");
        const auto v = s.ptrs (0);                                   // this vector allocates; do it before the mark
        // NOTHING that builds a message may run between the mark and the read: a std::string past libc++'s
        // small-buffer size allocates, and the counter would blame process() for the test's own prose.
        const long long mark = alloc::count.load();
        const bool okRun = d.process (v.data(), 1, (int) s.len());
        const bool okEmpty = d.process (v.data(), 1, 0);
        d.finish();
        const long long after = alloc::count.load();
        ok (okRun, "lifecycle: process");
        ok (okEmpty, "lifecycle: an n == 0 call is accepted and changes nothing");
        felitronics::test::okNoAlloc (after == mark, "lifecycle: process/finish allocate nothing");
        const HumReport r1 = d.report (0);
        d.finish();
        const HumReport r2 = d.report (0);
        ok (std::bit_cast<std::uint64_t> (r1.line.hz) == std::bit_cast<std::uint64_t> (r2.line.hz)
            && r1.quietStretches == r2.quietStretches && r1.harmonicsObserved == r2.harmonicsObserved,
            "lifecycle: finish() is idempotent");
        ok (! d.process (v.data(), 1, 16), "lifecycle: process after finish is refused");
        ok (! d.process (v.data(), 1, 0), "lifecycle: and so is an n == 0 call — law 11 checks finished first");
        ok (! d.process (v.data(), 2, 16), "lifecycle: more channels than prepared is refused");
        ok (! d.process (v.data(), 1, -1), "lifecycle: a negative length is refused");
        // ... and on an ACTIVE object, where the refusal is not just the finished check, and it moves nothing
        {
            HumDetector act; act.setParams (p);
            ok (act.prepare (kSr, 4096, 1), "lifecycle: an active object to refuse on");
            ok (act.process (v.data(), 1, 4096), "lifecycle: a good call first");
            const std::int64_t before2 = act.samplesProcessed();
            ok (! act.process (v.data(), 2, 16), "lifecycle: too many channels refused while ACTIVE");
            ok (! act.process (v.data(), 1, -1), "lifecycle: a negative length refused while ACTIVE");
            ok (act.samplesProcessed() == before2, "lifecycle: and a refused call consumed nothing");
        }
        // reset, then replay: bit-identical to a fresh object
        d.reset();
        ok (d.samplesProcessed() == 0 && d.stretchCount (0) == 0 && ! d.isFinished(), "lifecycle: reset re-anchors");
        const Out fresh = runSliced (s, p, { (int) kLen }, 4096);
        ok (fresh.bits.size() > 1000, "replay: the fresh reference is the whole report ("
                                      + std::to_string (fresh.bits.size()) + " words)");
        {
            // `d` has already run the same signal, been finished, and been reset above — replay on it and
            // compare EVERY field, not one frequency. A run that populated frames, stretches and the trace
            // first is the point: a reset() that forgot traceCount_ used to survive the whole suite.
            const auto v2 = s.ptrs (0);
            ok (d.process (v2.data(), 1, (int) s.len()), "replay: process after reset");
            d.finish();
            Out reused; reused.ok = fresh.ok;
            serialise (d, 1, p.maxHarmonic, reused);
            ok (reused == fresh, "replay: a reused object's WHOLE report is bit-identical to a fresh one");
        }
    }

    return felitronics::test::report();
}

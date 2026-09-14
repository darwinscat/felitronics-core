// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// SourceForensics self-tests — the spectral wall and the sample grid, and above all the fixtures on which a
// naive version of either LIES CONFIDENTLY. The positive cases are the cheap half; the negative ones are
// why this file is long.
//
//   · a REFERENCE NULL on the accumulated spectrum: the Welch mean against one computed independently here
//     (its own frame schedule k*hop, its own window, its own plain summation)
//   · a SECOND ORACLE of a different construction for the edge analysis: the whole cell / median-3 /
//     suffix-rank / quantile / scan pipeline re-implemented over the published mean spectrum, and its
//     coordinates compared with the report's. Two oracles because one is not enough — the analytic
//     fixtures say WHAT the answer is, the re-implementation says the published answer is the one the
//     documented rules produce, and only two constructions disagreeing exposes an under-read (P11).
//   · THE TRAPS, each a fixture: a deep notch WITH RECOVERY (a band-stop, a comb, a room null) is not a
//     wall; a 2 kHz-wide notch is not a wall; a monotone -15 dB/octave tilt is not a wall AND is not given
//     a frequency; flat full-band noise has no edge; a notch BELOW a real wall does not steal it; nested
//     edges report BOTH; a narrow sustained line above a real wall is forgiven up to a published rank and
//     no further, and what was forgiven is published
//   · the grid: 8/16/24-bit, 16-bit DITHERED into 24 (which must NOT read as 16), 32-bit PCM in float32
//     (which must read as "no <= 24-bit grid", never as 32) AND the documented limit of that — codes near
//     full scale ARE 24-bit once they are float32 — a non-dyadic gain, exactly +-1.0, one denormal in an
//     otherwise 16-bit programme (which the histogram catches where the max cannot), 256 distinct values
//     on a 16-bit grid, exact zeros, subnormals, and the distinct count against a std::set oracle
//     including its exact-capacity and overflow behaviour
//   · LAW 8a: the FULL trace — every event, every counter, every raw accumulator bin and its Neumaier
//     compensation, then every field of the final report — compared BIT-EXACTLY (std::bit_cast, field by
//     field, never memcmp of a padded struct and never `!=`) across 15 slicings, 4 prepared maxBlock
//     values, seeded ragged partitions, and call boundaries placed exactly on the transitions that matter
//   · lifecycle, refusals (including the ranges that would otherwise be undefined behaviour), the tail
//     contract, a disappearing channel, no allocation in process/finish
//
// FIXTURE NOTES, both learned the expensive way in this file's own measure-off:
//   · A DIGITAL LOW-PASS IS NOT A GENTLE-ROLLOFF FIXTURE. eq::Crossover2's SVF has a (1+z^-1)^2 numerator,
//     hence a double zero AT Nyquist, so its response really does fall off a cliff there and this
//     instrument really does find a sharp edge at 0.94 of Nyquist. The gentle case therefore has to be
//     built without a filter — an additive -15 dB/octave tilt — and the LR4 case is kept as a POSITIVE
//     test of the Nyquist zero instead.
//   · A FIXTURE'S OWN FLOAT32 FLOOR SETS THE DEEPEST MEASURABLE DROP. The nested-edge fixture first used a
//     -125 dB shelf, whose floor above the outer edge was the float32 quantisation of the fixture array at
//     -147 dB: a 22 dB drop, below minDropDb, so the outer edge was correctly refused and the test looked
//     like a code defect. The shelf is now -99 dB. Anything quieter than about -130 dB in a fixture is
//     measuring the fixture.
//
// MUTANT PASS — run once, recorded here rather than asserted. A gate that has never been red is not a
// gate. Each mutant was built with the object file DELETED first and refused a verdict unless a
// `Building CXX` line appeared: a header edited in the same SECOND as the previous build is invisible to
// make's 1-second mtime granularity, and that alone cost three false greens elsewhere today.
//   (1) the frame's power accumulated at the END of process() instead of at frame close  -> RED, 22 of 422
//   (2) a holed frame accumulated anyway (the frameFinite check dropped)                  -> RED,  4 of 422
//   (3) the aggregate dividing every channel by channel 0's frame count (pooling)         -> RED,  1 of 422
//   (4) the k <= 23 gate removed, so a finer grid reports a PCM word length anyway        -> RED,  5 of 422
//   (5) exact zeros no longer skipped when witnessing the grid                            -> RED, 13 of 422
//   (6) the distinct set marked incomplete on REACHING the limit, not on a new key it      -> RED,  1 of 422
//       could not store
//   (7) the plateau quantile taken as the MAXIMUM instead of the median                    -> RED, 17 of 422
//   (8) the suffix rank replaced by the cell just above (the pre-consilium construction)   -> RED, 13 of 422
//   (9) the top cell left as a partial cell of its own instead of absorbing the remainder  -> RED,  2 of 422
//  (10) the grid coordinates taken after the clock advanced instead of before              -> RED,  1 of 422
//  (11) the STRICT floor taken from the median-filtered cells instead of the raw ones      -> RED,  6 of 422
//  (12) the measurement reading the PENDING parameters instead of prepare()'s snapshot     -> RED,  1 of 422
//  (13) pcmCompatible ignoring the PCM range, so +1.0 gets a word length                   -> RED,  1 of 422
//  (14) a failed prepare() leaving the finished report armed                               -> RED,  3 of 422
//  (15) the emptiness minimum back to one bin, so "empty above Nyquist" can be claimed     -> RED,  1 of 422
// Fifteen of fifteen red, and the baseline rebuilt green afterwards. (11) through (15) exist because the
// CODE-REVIEW round found those five defects and the first ten mutants did not cover them; (11) in
// particular came back GREEN on the first stand, which is what proved the suite had a hole there rather
// than a gate. Mutant (8) is the one that matters most: it
// restores the construction the design consilium killed, and it is red on the notch fixtures alone.
// The mutant the brief names for every analyzer of this wave — moving the denormal flush to the end of
// process() — does not exist here: this analyzer has no feedback state, no IIR and no denormal cadence,
// so there is nothing to flush and nothing for a StateGrid to clock. Mutant (1) is its structural
// equivalent (state maintenance moved off an audio-time coordinate onto a call boundary).

#include <felitronics_test.h>
#include <felitronics/analysis/SourceForensics.h>
#include <felitronics/core/OfflineFft.h>
#include <felitronics/eq/Crossover2.h>

#include <algorithm>
#include <atomic>
#include <bit>
#include <cmath>
#include <complex>
#include <cstdint>
#include <limits>
#include <random>
#include <set>
#include <string>
#include <vector>

static std::atomic<long> g_allocs { 0 };
void* operator new (std::size_t s) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s); }
void  operator delete (void* p) noexcept { std::free (p); }
void  operator delete (void* p, std::size_t) noexcept { std::free (p); }

using namespace felitronics;
using analysis::ForensicsEvent;
using analysis::ForensicsReason;
using analysis::SampleGrid;
using analysis::SourceForensics;
using analysis::SourceForensicsParams;
using analysis::SpectralWall;

namespace
{

constexpr double kFs = 48000.0;
constexpr int    kOrder = 12;                      // N = 4096 at 48 kHz: 11.72 Hz bins, 46.875 Hz cells
constexpr std::size_t kLen = 16384;                // 7 frames at hop N/2

//==============================================================================
// fixtures
//==============================================================================

// Built in the FREQUENCY domain: a spectrum of length L with seeded random phases and an amplitude given
// by `amp(f)`, inverse-transformed once and scaled to a 0.5 peak. Exactly band-limited by construction,
// O(L log L) rather than O(tones * L), and every bin populated, so a cell is a smooth estimate rather
// than a comb of tones straddling cell edges.
//
// No taper: the rule that demands one (CLAUDE.md, learned in the F6 true-peak work) is about an
// UNWINDOWED measurement finding a genuine reconstruction overshoot at a file boundary. Every frame here
// is Hann-windowed, so each frame tapers itself and no analysis window can measure the programme's edge.
template <typename Amp>
std::vector<float> spectral (std::size_t L, double fs, unsigned seed, Amp amp)
{
    std::vector<std::complex<double>> X (L, std::complex<double> {});
    std::mt19937 rng (seed);
    std::uniform_real_distribution<double> ph (0.0, 2.0 * core::kPi);
    for (std::size_t k = 1; k < L / 2; ++k)
    {
        const double f = (double) k * fs / (double) L;
        X[k] = std::polar (amp (f), ph (rng));
        X[L - k] = std::conj (X[k]);
    }
    core::offline::fftInplace (X, +1);
    double peak = 0.0;
    for (std::size_t i = 0; i < L; ++i) peak = std::max (peak, std::fabs (X[i].real()));
    const double g = peak > 0.0 ? 0.5 / peak : 0.0;
    std::vector<float> x (L, 0.0f);
    for (std::size_t i = 0; i < L; ++i) x[(std::size_t) i] = (float) (X[i].real() * g);
    return x;
}

// A 1/f-tilted passband up to cutHz, then a flat floor `floorDb` below the level the passband had there.
std::vector<float> wallFixture (double cutHz, double floorDb, unsigned seed = 3, double fs = kFs,
                                std::size_t L = kLen)
{
    const double floorLin = (1000.0 / cutHz) * std::pow (10.0, floorDb / 20.0);
    return spectral (L, fs, seed, [=] (double f) { return f < cutHz ? 1000.0 / f : floorLin; });
}

// NON-stationary, structured, seeded — bursts against a bed, lead-in silence, a quiet tail. A steady tone
// is a weak witness of invariance, so this is what the law-8a trace runs on.
std::vector<std::vector<float>> ragged (std::size_t n, int nch, unsigned seed)
{
    std::mt19937 rng (seed);
    std::uniform_real_distribution<float> u (-1.0f, 1.0f);
    std::vector<std::vector<float>> x ((std::size_t) nch, std::vector<float> (n, 0.0f));
    for (std::size_t i = 0; i < n; ++i)
    {
        const double t = (double) i / kFs;
        const double e = i < n / 8 ? 0.0 : i > n * 7 / 8 ? 0.02 : (i % 967 < 41 ? 1.0 : 0.25);
        for (int c = 0; c < nch; ++c)
        {
            const double f = c == 0 ? 1237.0 : 613.0;
            x[(std::size_t) c][i] = (float) (e * (0.7 * std::sin (2.0 * core::kPi * f * t) + 0.3 * (double) u (rng)));
        }
    }
    return x;
}

// Onto the 2^-(bits-1) grid — what integer PCM of that depth looks like once it is float32.
void quantise (std::vector<float>& x, int bits)
{
    const double q = std::ldexp (1.0, bits - 1);
    for (auto& v : x) v = (float) (std::round ((double) v * q) / q);
}

void fillRange (std::vector<float>& x)                         // use the whole range before quantising
{
    double peak = 0.0;
    for (float v : x) peak = std::max (peak, std::fabs ((double) v));
    if (peak > 0.0) for (auto& v : x) v = (float) (0.9 * (double) v / peak);
}

//==============================================================================
// oracle 1 — the Welch mean spectrum, computed independently of the object
//==============================================================================

std::vector<double> welchMean (const std::vector<float>& x, int order, std::int64_t hop, std::int64_t& frames)
{
    const std::int64_t n = (std::int64_t) 1 << order;
    const int bins = (int) (n / 2 + 1);
    std::vector<double> w ((std::size_t) n);
    double sumW2 = 0.0;
    for (std::int64_t i = 0; i < n; ++i)
    {
        w[(std::size_t) i] = 0.5 - 0.5 * std::cos (2.0 * core::kPi * (double) i / (double) n);
        sumW2 += w[(std::size_t) i] * w[(std::size_t) i];
    }
    std::vector<double> acc ((std::size_t) bins, 0.0);
    std::vector<std::complex<double>> buf ((std::size_t) n);
    frames = 0;
    for (std::int64_t k = 0; k * hop + n <= (std::int64_t) x.size(); ++k)     // the schedule, from k*hop
    {
        const std::int64_t start = k * hop;
        for (std::int64_t i = 0; i < n; ++i)
            buf[(std::size_t) i] = std::complex<double> ((double) x[(std::size_t) (start + i)] * w[(std::size_t) i], 0.0);
        core::offline::fftInplace (buf, -1);
        const double norm = 1.0 / ((double) n * sumW2);
        for (int b = 0; b < bins; ++b)
        {
            const auto& z = buf[(std::size_t) b];
            acc[(std::size_t) b] += (z.real() * z.real() + z.imag() * z.imag()) * norm;
        }
        ++frames;
    }
    if (frames > 0) for (auto& v : acc) v /= (double) frames;
    return acc;
}

//==============================================================================
// oracle 2 — the edge analysis, re-implemented from the documented rules in a different shape
//==============================================================================

struct WallOracle
{
    bool found = false, clipped = false;
    double cutoffHz = 0.0, steepestHz = 0.0, endHz = 0.0, dropDb = 0.0, plateau = 0.0, maxAbove = 0.0;
    int cellCount = 0, binsPerCell = 0;
};

double median3 (double a, double b, double c)                    // by sorting, not by the header's min/max
{
    double v[3] { a, b, c };
    std::sort (std::begin (v), std::end (v));
    return v[1];
}

// NEAREST-RANK on a sorted copy: rank ceil(q*n), 1-based, clamped. Spelled with integers on purpose —
// ceil(0.9 * 40) in floating point is 37, not 36, because 0.9 is not representable.
double quantileOf (std::vector<double> v, int num, int den)
{
    if (v.empty()) return 0.0;
    std::sort (v.begin(), v.end());
    const std::size_t rank = (v.size() * (std::size_t) num + (std::size_t) den - 1u) / (std::size_t) den;
    return v[std::min (std::max<std::size_t> (rank, 1u), v.size()) - 1u];
}

WallOracle oracleWall (const std::vector<double>& mean, const SourceForensicsParams& p, double fs)
{
    WallOracle o;
    const int bins = (int) mean.size();
    const double binHz = fs / (double) ((std::int64_t) 1 << p.fftOrder);
    int bpc = std::max (1, (int) (p.cellWidthHz / binHz));
    bpc = std::min (bpc, bins);
    const double cellHz = (double) bpc * binHz;
    const int cellCount = std::max (1, bins / bpc);
    o.cellCount = cellCount;
    o.binsPerCell = bpc;
    std::vector<double> cells ((std::size_t) cellCount, 0.0);
    for (int j = 0; j < cellCount; ++j)
    {
        const int from = j * bpc, to = j == cellCount - 1 ? bins : from + bpc;
        double acc = 0.0;
        for (int b = from; b < to; ++b) acc += mean[(std::size_t) b];
        cells[(std::size_t) j] = acc / (double) (to - from);
    }
    std::vector<double> sm = cells;
    for (int j = 1; j < cellCount - 1; ++j)
        sm[(std::size_t) j] = median3 (cells[(std::size_t) (j - 1)], cells[(std::size_t) j], cells[(std::size_t) (j + 1)]);
    // the floor above each boundary: the (t+1)-th largest of the suffix, t <= exemptCells and <= a tenth
    std::vector<double> flo ((std::size_t) cellCount, 0.0);
    for (int j = 0; j < cellCount; ++j)
    {
        std::vector<double> above (sm.begin() + j, sm.end());
        std::sort (above.rbegin(), above.rend());
        const int t = std::min (p.exemptCells, ((int) above.size() - 1) / 10);
        flo[(std::size_t) j] = above[(std::size_t) t];
    }
    const int plateauCells = std::max (1, (int) (p.plateauSpanHz / cellHz));
    const int floorCells = std::max (1, (int) (p.floorSpanHz / cellHz));
    const int first = std::max (plateauCells, (int) std::ceil (p.searchFromHz / cellHz));
    auto spanOf = [&] (int from, int to)
    {
        std::vector<double> v;
        for (int j = std::max (0, from); j < std::min (to, cellCount); ++j) v.push_back (sm[(std::size_t) j]);
        return v;
    };
    int best = -1;
    double bp = 0.0, bm = 0.0;
    for (int j = first; j < cellCount; ++j)
    {
        const double pl = quantileOf (spanOf (j - plateauCells, j), 1, 2);
        const double mx = flo[(std::size_t) j];
        if (! (pl > 0.0) || ! (pl > mx)) continue;
        if (best < 0 || pl * bm > bp * mx) { best = j; bp = pl; bm = mx; }
    }
    if (best < 0) return o;
    o.steepestHz = (double) best * cellHz;
    o.plateau = bp;
    o.maxAbove = bm;
    const int floorTo = std::min (best + floorCells, cellCount);
    const double floorLocal = quantileOf (spanOf (best, floorTo), 1, 2);
    const double startRatio = std::pow (10.0, -p.transitionStartDb / 10.0);
    const double endRatio = std::pow (10.0, p.transitionEndDb / 10.0);
    int js = -1;
    for (int j = best; j > best - plateauCells; --j)
        if (sm[(std::size_t) (j - 1)] >= o.plateau * startRatio) { js = j; break; }
    if (js < 0) { o.clipped = true; js = best - plateauCells; }
    int je = -1;
    for (int j = js; j < floorTo; ++j)
        if (sm[(std::size_t) j] <= floorLocal * endRatio) { je = j; break; }
    if (je < 0) { o.clipped = true; je = floorTo; }
    o.cutoffHz = (double) js * cellHz;
    o.endHz = (double) je * cellHz;
    o.dropDb = 10.0 * std::log10 (o.plateau / o.maxAbove);
    o.found = o.dropDb >= p.minDropDb;
    return o;
}

//==============================================================================
// the law-8a trace
//==============================================================================

struct Trace { std::vector<std::uint64_t> w; };

void put (Trace& tr, std::int64_t v) { tr.w.push_back ((std::uint64_t) v); }
void putD (Trace& tr, double v) { tr.w.push_back (std::bit_cast<std::uint64_t> (v)); }
void putB (Trace& tr, bool v) { tr.w.push_back (v ? 1u : 0u); }

void putGrid (Trace& tr, const SampleGrid& g)
{
    putB (tr, g.valid); put (tr, (std::int64_t) g.reason); put (tr, g.gridExponent);
    putB (tr, g.pcmCompatible); putB (tr, g.outsidePcmRange); put (tr, g.minExactPcmBits);
    putD (tr, g.absPeak); putD (tr, g.sampleMin); putD (tr, g.sampleMax);
    put (tr, g.nonZeroSamples); put (tr, g.zeroSamples); put (tr, g.nonFiniteSamples);
    put (tr, g.absentSamples); put (tr, g.offGridSamples); put (tr, g.firstOffGridSample);
    put (tr, g.firstMaxGridSample); put (tr, g.distinctValues); putB (tr, g.distinctComplete);
}

void putWall (Trace& tr, const SpectralWall& w)
{
    putB (tr, w.valid); put (tr, (std::int64_t) w.reason); putB (tr, w.sharp); putB (tr, w.nearNyquist);
    putD (tr, w.cutoffHz); putD (tr, w.cutoffFractionOfNyquist); putD (tr, w.steepestHz);
    putD (tr, w.transitionEndHz); putD (tr, w.transitionHz);
    putB (tr, w.transitionClipped); putB (tr, w.truncatedAtNyquist);
    putD (tr, w.plateauPower); putD (tr, w.floorLocalPower); putD (tr, w.maxAbovePower);
    putD (tr, w.sufMaxPower); put (tr, w.exemptedCells);
    putD (tr, w.dropDb); putD (tr, w.strictDropDb); putD (tr, w.localDropDb); putD (tr, w.recoveryDb);
    putD (tr, w.plateauSpreadDb); putD (tr, w.steepnessDbPerOctave);
    putB (tr, w.secondValid); putB (tr, w.secondSharp);
    putD (tr, w.secondCutoffHz); putD (tr, w.secondDropDb); putD (tr, w.secondTransitionHz);
    putB (tr, w.emptyAboveValid); putD (tr, w.emptyAboveHz); putD (tr, w.emptyAboveFractionOfNyquist);
    putD (tr, w.emptyThresholdPower); putD (tr, w.peakCellPower);
    putD (tr, w.binHz); putD (tr, w.cellHz); putD (tr, w.nyquistHz);
    putD (tr, w.searchFromHz); putD (tr, w.searchToHz);
    put (tr, w.framesUsed); put (tr, w.framesHoled);
}

// Filled at the CANONICAL moment — the instant a frame closes, or the instant an off-grid sample or a
// saturating value is met — never at the exit of process(), where a re-slicing defect is invisible.
void observer (void* user, const SourceForensics& sf, ForensicsEvent ev, int channel) noexcept
{
    Trace& tr = *static_cast<Trace*> (user);
    put (tr, (std::int64_t) ev);
    put (tr, channel);
    put (tr, sf.samplesProcessed());
    put (tr, sf.frames().frameIndex());
    put (tr, sf.frames().frameStart());
    put (tr, sf.tailUncoveredSamples());
    for (int c = 0; c < sf.channels(); ++c)
    {
        putB (tr, sf.lastFrameUsed (c));
        put (tr, sf.framesUsed (c));
        put (tr, sf.framesHoled (c));
        putB (tr, sf.frames().frameFinite (c));
        putGrid (tr, sf.sampleGrid (c));
        const std::int64_t* h = sf.gridExponentHistogram (c);
        for (int k = 0; k < SourceForensics::gridExponentBuckets(); ++k) put (tr, h[(std::size_t) k]);
        const double* ps = sf.powerSum (c);
        const double* pc = sf.powerSumCompensation (c);
        for (int b = 0; b < sf.bins(); ++b) { putD (tr, ps[(std::size_t) b]); putD (tr, pc[(std::size_t) b]); }
    }
}

struct RunSpec
{
    std::vector<int> slices;
    int maxBlock = 512;
    int order = 8;
    int hop = 0;
    int distinct = 1 << 10;
    int channels = 2;
    std::int64_t narrowFrom = -1;   // from this SAMPLE on, the calls carry one channel (a disappearing channel)
};

// One full run: the trace of every canonical event, then every field of the final report.
std::vector<std::uint64_t> runTrace (const std::vector<std::vector<float>>& x, const RunSpec& spec)
{
    Trace tr;
    tr.w.reserve (1u << 20);
    SourceForensics sf;
    SourceForensicsParams p;
    p.fftOrder = spec.order;
    p.hop = spec.hop;
    p.maxDistinctValues = spec.distinct;
    sf.setParams (p);
    sf.setObserver (&observer, &tr);
    if (! sf.prepare (kFs, spec.maxBlock, spec.channels)) { tr.w.push_back (0xDEADull); return tr.w; }
    const std::size_t total = x[0].size();
    std::vector<const float*> planes ((std::size_t) spec.channels);
    std::size_t at = 0, s = 0;
    while (at < total)
    {
        const int want = spec.slices[s % spec.slices.size()];
        ++s;
        std::size_t take = want <= 0 ? 0u : std::min ((std::size_t) want, total - at);
        // The channel-presence mask is INPUT, not slicing: it is a function of the absolute sample index,
        // so a call that would straddle the transition is cut at it. Without this the "same input, other
        // slicing" premise is false and the test asserts a theorem that does not hold.
        int nch = spec.channels;
        if (spec.narrowFrom >= 0)
        {
            if ((std::int64_t) at >= spec.narrowFrom) nch = 1;
            else if ((std::int64_t) (at + take) > spec.narrowFrom) take = (std::size_t) (spec.narrowFrom - (std::int64_t) at);
        }
        for (int c = 0; c < nch; ++c) planes[(std::size_t) c] = x[(std::size_t) c].data() + at;
        if (! sf.process (planes.data(), nch, (int) take)) tr.w.push_back (0xBADDull);
        at += take;
    }
    sf.finish();
    put (tr, sf.samplesProcessed());
    put (tr, sf.tailUncoveredSamples());
    put (tr, sf.frames().frameCount());
    for (int c = 0; c < spec.channels; ++c) { putWall (tr, sf.wall (c)); putGrid (tr, sf.sampleGrid (c)); }
    putWall (tr, sf.wall());
    return tr.w;
}

// An independent grid-exponent oracle: multiply by two until the value is integral. Exact over the whole
// range, subnormals included (2^-149 needs 149 doublings, every one of them exact in binary64).
int oracleGridK (float x)
{
    double v = (double) x;
    for (int k = 0; k < 400; ++k)
    {
        if (core::exactlyEqual (std::floor (v), v)) return k;
        v *= 2.0;
    }
    return -1;
}

// Measure one mono programme with the given params.
struct Measured
{
    SpectralWall w;
    SampleGrid g;
    double cellHz = 0.0;
    int bins = 0, binsPerCell = 0, cellCount = 0;
    std::vector<double> mean;                 // the published per-bin Welch mean, for the oracles below
};

Measured measure (const std::vector<float>& x, SourceForensicsParams p, double fs = kFs, int maxBlock = 512)
{
    SourceForensics sf;
    sf.setParams (p);
    Measured m;
    if (! sf.prepare (fs, maxBlock, 1)) return m;
    const float* in[1] { x.data() };
    if (! sf.process (in, 1, (int) x.size())) return m;
    sf.finish();
    m.w = sf.wall (0);
    m.g = sf.sampleGrid (0);
    m.cellHz = sf.cellHz();
    m.bins = sf.bins();
    m.binsPerCell = sf.binsPerCell();
    m.cellCount = sf.cellCount();
    m.mean.resize ((std::size_t) m.bins);
    for (int b = 0; b < m.bins; ++b) m.mean[(std::size_t) b] = sf.meanPower (0, b);
    return m;
}

// The maximum RAW cell at or above the boundary `fromHz`, recomputed from the published mean spectrum.
// This is an ORACLE for `sufMaxPower`: the same quantity taken from the median-FILTERED cells is a
// DIFFERENT number, and the difference is the whole reason the strict floor is published at all.
double rawCellMaxFrom (const Measured& m, double fromHz)
{
    const int jFrom = (int) std::llround (fromHz / m.cellHz);
    double best = 0.0;
    for (int j = jFrom; j < m.cellCount; ++j)
    {
        const int from = j * m.binsPerCell, to = j == m.cellCount - 1 ? m.bins : from + m.binsPerCell;
        double acc = 0.0;
        for (int b = from; b < to; ++b) acc += m.mean[(std::size_t) b];
        best = std::max (best, acc / (double) (to - from));
    }
    return best;
}

SourceForensicsParams defaults()
{
    SourceForensicsParams p;
    p.fftOrder = kOrder;
    return p;
}

} // namespace

//==============================================================================

int main()
{
    using felitronics::test::approx;
    using felitronics::test::ok;
    using felitronics::test::okNoAlloc;
    using felitronics::test::run;

    // ---------- 1. the published geometry, and no coordinate above Nyquist ----------
    {
        SourceForensics sf;
        auto p = defaults();
        sf.setParams (p);
        ok (run (sf.prepare (kFs, 512, 2)), "geometry: prepare");
        ok (sf.bins() == (1 << (kOrder - 1)) + 1, "geometry: bins == N/2 + 1");
        approx (sf.binHz(), kFs / 4096.0, 1e-12, "geometry: binHz");
        approx (sf.cellHz(), (double) sf.binsPerCell() * sf.binHz(), 1e-12, "geometry: a cell is a whole number of bins");
        ok (sf.cellHz() <= p.cellWidthHz, "geometry: the effective cell is no wider than requested");
        ok (sf.cellCount() * sf.binsPerCell() <= sf.bins(), "geometry: the cells fit inside the bins");
        ok ((sf.cellCount() + 1) * sf.binsPerCell() > sf.bins(), "geometry: and the top cell absorbs the remainder");
        ok (sf.searchToHz() <= kFs * 0.5, "geometry: the top candidate is at or below Nyquist");
        ok (sf.searchFromHz() >= p.searchFromHz, "geometry: the effective search floor is published");
        approx (sf.searchFromHz(), (double) sf.plateauSpanCells() * sf.cellHz(), 1e-12,
                "geometry: and a full plateau span is what pushes it up");
        ok (sf.distinctTableSlots() >= (std::size_t) sf.distinctLimit() * 4u / 3u,
            "geometry: the table leaves the distinct limit at or below 3/4 load");
        ok (sf.exemptCells() == p.exemptCells, "geometry: the exemption rank is published");
        // A published emptiness coordinate must be STRICTLY below Nyquist at EVERY order, or the claim can
        // be "everything above Nyquist is empty" — vacuously true and read as a finding. At fftOrder 8 the
        // 200 Hz minimum floors to one bin, which is exactly how that happened (code-review round).
        for (int order : { 8, 9, 10, 12 })
        {
            SourceForensics probe;
            auto q = defaults();
            q.fftOrder = order;
            probe.setParams (q);
            ok (run (probe.prepare (kFs, 64, 1)), "empty-floor: prepare at order " + std::to_string (order));
            std::vector<float> x (4096, 0.0f);
            const double f = (double) (probe.bins() - 3) * probe.binHz();      // a tone in the top bins
            for (std::size_t i = 0; i < x.size(); ++i)
                x[(std::size_t) i] = (float) (0.5 * std::sin (2.0 * core::kPi * f * (double) i / kFs));
            const float* in[1] { x.data() };
            ok (run (probe.process (in, 1, (int) x.size())), "empty-floor: process at order " + std::to_string (order));
            probe.finish();
            const auto w = probe.wall (0);
            ok (! w.emptyAboveValid || w.emptyAboveHz < w.nyquistHz,
                "empty-floor: order " + std::to_string (order) + " never claims emptiness above Nyquist itself ("
                + std::to_string (w.emptyAboveHz) + " Hz)");
        }
        const auto st = SourceForensics::storageFor (kFs, 2, p);
        ok (st.ok && st.bytes() > 0, "storage: published before the allocation");
        ok (st.frames.ok && st.bytes() > st.frames.bytes(), "storage: the frame producer's budget is included");
    }

    // ---------- 2. REFERENCE NULL: the Welch mean against an independent one ----------
    {
        auto x = wallFixture (10000.0, -80.0);
        SourceForensics sf;
        sf.setParams (defaults());
        ok (run (sf.prepare (kFs, 997, 1)), "null: prepare");
        const float* in[1] { x.data() };
        ok (run (sf.process (in, 1, (int) x.size())), "null: process");
        sf.finish();
        std::int64_t oracleFrames = 0;
        const auto want = welchMean (x, kOrder, 1 << (kOrder - 1), oracleFrames);
        ok (oracleFrames > 4, "null: the oracle produced frames (" + std::to_string (oracleFrames) + ")");
        ok (sf.framesUsed (0) == oracleFrames, "null: the same frame count as the independent schedule");
        double worst = 0.0;
        bool finite = true;
        for (int b = 0; b < sf.bins(); ++b)
        {
            const double got = sf.meanPower (0, b), wnt = want[(std::size_t) b];
            const double d = std::fabs (got - wnt) / (wnt > 1e-30 ? wnt : 1e-30);
            if (! std::isfinite (d)) finite = false;
            if (d > worst) worst = d;
        }
        ok (finite, "null: every relative error is finite");
        ok (worst < 1e-9, "null: worst relative error " + std::to_string (worst) + " < 1e-9");
    }

    // ---------- 3. the sharp walls, analytically ----------
    {
        struct Case { double cut; double floorDb; double minDrop; double maxTrans; bool nearNyq; };
        const Case cases[] = {
            {  6000.0, -80.0, 55.0,  400.0, false },
            { 10000.0, -80.0, 55.0,  400.0, false },
            { 16000.0, -80.0, 55.0,  400.0, false },
            { 20500.0, -99.0, 70.0,  700.0, true  },
        };
        for (const auto& cs : cases)
        {
            const auto m = measure (wallFixture (cs.cut, cs.floorDb), defaults());
            const auto& w = m.w;
            const std::string tag = "wall " + std::to_string ((int) cs.cut) + ": ";
            ok (w.valid && w.reason == ForensicsReason::Ok, tag + "valid");
            ok (w.sharp, tag + "the edge is sharp");
            ok (std::fabs (w.cutoffHz - cs.cut) <= 3.0 * m.cellHz,
                tag + "cutoff " + std::to_string (w.cutoffHz) + " within 3 cells of the construction");
            ok (w.cutoffHz <= cs.cut + m.cellHz, tag + "and never ABOVE the constructed edge");
            ok (w.dropDb >= cs.minDrop, tag + "drop " + std::to_string (w.dropDb) + " >= " + std::to_string (cs.minDrop));
            ok (w.transitionHz <= cs.maxTrans, tag + "transition " + std::to_string (w.transitionHz) + " Hz is narrow");
            ok (! w.transitionClipped, tag + "the transition is measured, not censored");
            ok (w.nearNyquist == cs.nearNyq, tag + "nearNyquist agrees with the position");
            ok (w.cutoffHz <= w.nyquistHz && w.emptyAboveHz <= w.nyquistHz, tag + "no coordinate exceeds Nyquist");
            approx (w.cutoffFractionOfNyquist, w.cutoffHz / (kFs * 0.5), 1e-12, tag + "the fraction is the Hz over Nyquist");
            ok (w.emptyAboveValid && w.emptyAboveHz > cs.cut - 4.0 * m.cellHz,
                tag + "empty above " + std::to_string (w.emptyAboveHz));
            ok (w.plateauSpreadDb < 12.0, tag + "the plateau reference was flat (" + std::to_string (w.plateauSpreadDb) + " dB)");
            ok (w.steepnessDbPerOctave > 100.0, tag + "steepness " + std::to_string (w.steepnessDbPerOctave) + " dB/oct");
            ok (w.strictDropDb >= w.dropDb - 30.0,
                tag + "the raw suffix maximum is close to the forgiven one (strict "
                    + std::to_string (w.strictDropDb) + " vs " + std::to_string (w.dropDb) + " dB)");
            ok (core::exactlyEqual (w.sufMaxPower, rawCellMaxFrom (m, w.steepestHz)),
                tag + "and it IS the maximum raw cell above the boundary, recomputed from the mean spectrum");
            ok (w.strictDropDb > 0.0 && w.strictDropDb <= w.dropDb,
                tag + "the strict drop is published and is never the larger");
            ok (! w.secondValid, tag + "and there is no second edge above it");
        }
        // an edge placed EXACTLY on a cell boundary, where an off-by-one in the cell map would hide
        {
            SourceForensics probe;
            probe.setParams (defaults());
            ok (run (probe.prepare (kFs, 64, 1)), "on-grid: prepare");
            const double onGrid = 200.0 * probe.cellHz();               // 9375 Hz: a cell edge exactly
            const auto m = measure (wallFixture (onGrid, -80.0), defaults());
            ok (m.w.sharp && std::fabs (m.w.cutoffHz - onGrid) <= 2.0 * m.cellHz,
                "on-grid: an edge on a cell boundary is found at it (" + std::to_string (m.w.cutoffHz) + ")");
        }
        // the aggregate of a mono programme IS that channel, bit for bit
        {
            auto x = wallFixture (12000.0, -80.0);
            SourceForensics sf;
            sf.setParams (defaults());
            ok (run (sf.prepare (kFs, 512, 1)), "mono aggregate: prepare");
            const float* in[1] { x.data() };
            ok (run (sf.process (in, 1, (int) x.size())), "mono aggregate: process");
            sf.finish();
            ok (core::exactlyEqual (sf.wall().cutoffHz, sf.wall (0).cutoffHz)
                && core::exactlyEqual (sf.wall().dropDb, sf.wall (0).dropDb),
                "mono aggregate: the mean of one mean is that mean, bit for bit");
        }
    }

    // ---------- 4. THE TRAPS: what must NOT be an edge, and must not even get a frequency ----------
    {
        struct Neg { const char* name; std::vector<float> x; };
        std::vector<Neg> negs;
        // a deep notch WITH RECOVERY: full band to Nyquist, minus a gap
        for (double gapHi : { 10400.0, 12000.0 })
            negs.push_back ({ gapHi < 11000.0 ? "notch 400 Hz" : "notch 2 kHz",
                              spectral (kLen, kFs, 5, [=] (double f)
                              { return (f >= 10000.0 && f < gapHi) ? 1.0e-5 : 1000.0 / f; }) });
        // a monotone -15 dB/octave tilt, to Nyquist. Not a filter: a digital low-pass has a double zero AT
        // Nyquist, which IS an edge, so a filtered fixture cannot test "gentle" (see the header note).
        negs.push_back ({ "tilt -15 dB/oct", spectral (kLen, kFs, 7, [] (double f)
                          { return std::pow (1000.0 / f, 2.5); }) });
        negs.push_back ({ "flat noise", spectral (kLen, kFs, 11, [] (double) { return 1.0; }) });
        for (auto& ng : negs)
        {
            const auto m = measure (ng.x, defaults());
            const std::string tag = std::string ("neg ") + ng.name + ": ";
            ok (! m.w.sharp, tag + "NOT a sharp edge");
            ok (! m.w.valid && m.w.reason == ForensicsReason::ShallowerThanMinDrop,
                tag + "and NOT given a valid frequency — the best descent is only "
                    + std::to_string (m.w.dropDb) + " dB (reason " + std::to_string ((int) m.w.reason) + ")");
            ok (m.w.dropDb < 24.0, tag + "the conservative drop stays under minDropDb");
            ok (m.w.peakCellPower > 0.0, tag + "while the evidence is still published");
        }
        // A single narrow tone is NOT in that list: a programme whose only content is 20 Hz wide IS
        // band-limited, and the instrument finds that band limit — at the window's skirt rather than at the
        // tone, so `sharp` is withheld while `valid` is not. Measured, not assumed.
        {
            auto x = spectral (kLen, kFs, 13, [] (double f) { return f > 990.0 && f < 1010.0 ? 1.0 : 0.0; });
            const auto m = measure (x, defaults());
            ok (m.w.valid, "one sine: a 20 Hz-wide programme has a real band limit");
            ok (! m.w.sharp, "one sine: but the edge the WINDOW leaves is too wide to call sharp ("
                             + std::to_string (m.w.transitionHz) + " Hz)");
            ok (m.w.cutoffHz > 1000.0 && m.w.cutoffHz < 3000.0,
                "one sine: and it sits just above the tone (" + std::to_string (m.w.cutoffHz) + " Hz)");
            ok (m.w.emptyAboveValid && m.w.emptyAboveHz < 2000.0, "one sine: with everything above it empty");
        }
        // the notch's own local evidence: a large LOCAL drop, a small conservative one. That difference is
        // the whole reason a notch is not a wall.
        {
            auto x = spectral (kLen, kFs, 5, [] (double f)
                               { return (f >= 10000.0 && f < 10400.0) ? 1.0e-5 : 1000.0 / f; });
            auto p = defaults();
            p.searchFromHz = 5000.0;
            const auto m = measure (x, p);
            ok (! m.w.valid, "notch: still no edge with the search floor raised to 5 kHz");
        }
    }

    // ---------- 5. a notch BELOW a real wall must not steal it; nested edges report BOTH ----------
    {
        auto x = spectral (kLen, kFs, 17, [] (double f)
                           { return f >= 16000.0 ? 1.0e-5 : ((f >= 8000.0 && f < 8400.0) ? 1.0e-5 : 1000.0 / f); });
        const auto m = measure (x, defaults());
        ok (m.w.sharp, "notch+wall: the wall is still found");
        ok (std::fabs (m.w.cutoffHz - 16000.0) <= 3.0 * m.cellHz,
            "notch+wall: at 16 kHz, not at the notch (" + std::to_string (m.w.cutoffHz) + ")");
    }
    {
        // nested: full band to 16 kHz, a -99 dB shelf to 20.5 kHz, nothing above
        // The shelf is 1e-4 in amplitude (-80 dB): at 1e-5 the region above the OUTER edge would be the
        // fixture's own float32 quantisation at -147 dB, a 22 dB drop, and the outer edge would be refused
        // for a reason that has nothing to do with the code (see the header's fixture note).
        auto x = spectral (kLen, kFs, 19, [] (double f)
                           { return f < 16000.0 ? 1000.0 / f : (f < 20500.0 ? 1000.0 / f * 1.0e-4 : 1.0e-9); });
        const auto m = measure (x, defaults());
        ok (m.w.sharp && std::fabs (m.w.cutoffHz - 16000.0) <= 3.0 * m.cellHz,
            "nested: the INNER edge wins the primary, its plateau being the programme itself ("
            + std::to_string (m.w.cutoffHz) + ")");
        ok (m.w.dropDb > 55.0, "nested: with the shelf as its floor (" + std::to_string (m.w.dropDb) + " dB)");
        ok (m.w.secondValid, "nested: and the OUTER edge is published too");
        ok (std::fabs (m.w.secondCutoffHz - 20500.0) <= 3.0 * m.cellHz,
            "nested: at 20.5 kHz (" + std::to_string (m.w.secondCutoffHz) + ")");
        ok (m.w.secondDropDb > 24.0, "nested: with its own drop (" + std::to_string (m.w.secondDropDb) + " dB)");
    }

    // ---------- 6. a narrow line above a real wall: forgiven to the published rank, and no further ------
    {
        // The line is one cell wide at the default geometry and is the loudest thing above the edge. The
        // exemption must recover the wall while the forgiveness stays visible; a line too loud to forgive
        // must cost the report its validity, NOT move the frequency.
        // (a) a genuine one-frequency line, up to -30 dB of the midband: the wall SURVIVES — the 3-cell
        //     median already removes a line that occupies one cell, and the suffix rank adds 1-2 dB on top
        //     of that (measured). The line is still fully visible in recoveryDb and in emptyAboveHz.
        for (double lvl : { 3.0e-3, 1.0e-2, 3.0e-2 })
        {
            auto x = spectral (kLen, kFs, 23, [=] (double f)
                               { return f < 16000.0 ? 1000.0 / f
                                      : (f > 18996.0 && f < 19002.0 ? lvl : 1.0e-5); });
            const auto m = measure (x, defaults());
            const std::string tag = "1-bin line " + std::to_string (lvl) + ": ";
            ok (m.w.valid && m.w.sharp, tag + "the wall survives a narrow line above it (drop "
                                            + std::to_string (m.w.dropDb) + " dB)");
            ok (std::fabs (m.w.cutoffHz - 16000.0) <= 3.0 * m.cellHz,
                tag + "and is still located at 16 kHz (" + std::to_string (m.w.cutoffHz) + ")");
            ok (m.w.recoveryDb > 5.0, tag + "while recoveryDb publishes the line ("
                                          + std::to_string (m.w.recoveryDb) + " dB)");
            ok (m.w.dropDb > m.w.strictDropDb, tag + "and the strict drop shows what the rank forgave ("
                                                   + std::to_string (m.w.strictDropDb) + " dB)");
            // The line is what separates the raw floor from the filtered one: taken from the median-filtered
            // cells this number would be tens of dB lower and would forgive the very thing it exists to show.
            ok (core::exactlyEqual (m.w.sufMaxPower, rawCellMaxFrom (m, m.w.steepestHz)),
                tag + "the strict floor is the RAW maximum above the boundary, filter and rank included");
            ok (m.w.sufMaxPower > m.w.maxAbovePower,
                tag + "which stands above the forgiving floor the search used");
            ok (m.w.exemptedCells > 0, tag + "with the rank that was skipped");
            ok (m.w.emptyAboveValid && m.w.emptyAboveHz > 18000.0,
                tag + "the per-BIN emptiness test sees the line itself (" + std::to_string (m.w.emptyAboveHz) + " Hz)");
        }
        // (b) a 100 Hz-wide BAND at -20 dB is not a line and is not forgiven: 2-3 cells of real content
        //     above the edge means the band is NOT limited there, so the report refuses — and the refusal
        //     is what keeps a wrong FREQUENCY from being published, which is what happens without it.
        {
            auto x = spectral (kLen, kFs, 23, [] (double f)
                               { return f < 16000.0 ? 1000.0 / f
                                      : (f > 18950.0 && f < 19050.0 ? 1.0e-2 : 1.0e-5); });
            const auto m = measure (x, defaults());
            ok (! m.w.valid && m.w.reason == ForensicsReason::ShallowerThanMinDrop,
                "wide line: 100 Hz of content above the edge is refused, not forgiven (drop "
                + std::to_string (m.w.dropDb) + " dB)");
            ok (m.w.localDropDb > 60.0, "wide line: while the LOCAL drop is published and large ("
                                        + std::to_string (m.w.localDropDb) + " dB)");
            ok (m.w.recoveryDb > 40.0, "wide line: and so is the recovery ("
                                       + std::to_string (m.w.recoveryDb) + " dB)");
        }
        // (c) a band loud enough to move the argmax: the coordinate itself becomes meaningless, which is
        //     precisely why `valid` gates it rather than the number being published as fact.
        {
            auto x = spectral (kLen, kFs, 23, [] (double f)
                               { return f < 16000.0 ? 1000.0 / f
                                      : (f > 18950.0 && f < 19050.0 ? 1.0e-1 : 1.0e-5); });
            const auto m = measure (x, defaults());
            ok (! m.w.valid, "loud line: refused");
            ok (m.w.cutoffHz < 5000.0, "loud line: and the argmax has indeed walked away from the wall ("
                                       + std::to_string (m.w.cutoffHz) + " Hz) — the reason `valid` exists");
        }
        // exemptCells = 0 turns the RANK off, leaving the 3-cell median as the only forgiveness there is —
        // which is precisely why the strict number is taken from the RAW cells and not from the filtered
        // ones: "unforgiven" computed from a median-filtered spectrum would be neither.
        {
            auto x = spectral (kLen, kFs, 23, [] (double f)
                               { return f < 16000.0 ? 1000.0 / f
                                      : (f > 18996.0 && f < 19002.0 ? 3.0e-3 : 1.0e-5); });
            auto p0 = defaults();
            p0.exemptCells = 0;
            const auto m0 = measure (x, p0);
            const auto m2 = measure (x, defaults());
            ok (m0.w.dropDb <= m2.w.dropDb,
                "exemptCells = 0: the drop is no larger than with the rank on ("
                + std::to_string (m0.w.dropDb) + " vs " + std::to_string (m2.w.dropDb) + " dB)");
            ok (m0.w.dropDb >= m0.w.strictDropDb,
                "exemptCells = 0: and the median filter's own forgiveness is what remains between the "
                "forgiving and the raw floor (" + std::to_string (m0.w.dropDb) + " vs "
                + std::to_string (m0.w.strictDropDb) + " dB)");
            ok (core::exactlyEqual (m0.w.sufMaxPower, m2.w.sufMaxPower),
                "exemptCells = 0: while the raw suffix maximum does not depend on the rank at all");
        }
    }

    // ---------- 7. a pre-wall taper biases the edge DOWNWARD, never upward ----------
    {
        // cutoffHz is "the last cell within transitionStartDb of the MEDIAN of the span below", so a
        // passband that slopes across that span puts the contour inside it. The bias is real, systematic
        // and one-directional; this pins the direction and the bound rather than pretending it is zero.
        for (double taperDb : { 0.0, 6.0, 12.0 })
        {
            auto x = spectral (kLen, kFs, 29, [=] (double f)
            {
                if (f >= 16000.0) return 1.0e-5;
                const double t = f > 14000.0 ? (f - 14000.0) / 2000.0 : 0.0;
                return 1000.0 / f * std::pow (10.0, -taperDb * t / 20.0);
            });
            const auto m = measure (x, defaults());
            const std::string tag = "taper " + std::to_string ((int) taperDb) + " dB: ";
            ok (m.w.valid, tag + "an edge is still found");
            ok (m.w.cutoffHz <= 16000.0 + m.cellHz, tag + "the reported edge is never ABOVE the true one");
            ok (m.w.cutoffHz >= 16000.0 - 1000.0 - (taperDb / 6.0) * 500.0,
                tag + "and the downward bias stays inside the documented bound (" + std::to_string (m.w.cutoffHz) + ")");
        }
    }

    // ---------- 8. a digital low-pass has a double zero AT Nyquist — a fixture lesson, kept as a test ----
    {
        auto x = spectral (kLen, kFs, 31, [] (double) { return 1.0; });
        eq::Crossover2 xo;                                    // the repo's LR4, not a fifth private copy
        xo.prepare (kFs, 1);
        xo.setFrequency (4000.0f);
        for (std::size_t i = 0; i < x.size(); ++i)
        {
            float lo = 0.0f, hi = 0.0f;
            xo.processSample (0, x[i], lo, hi);
            x[i] = lo;
        }
        const auto m = measure (x, defaults());
        ok (m.w.valid, "LR4: an edge is found");
        ok (m.w.cutoffHz > 0.85 * m.w.nyquistHz,
            "LR4: and it is the SVF's Nyquist zero (" + std::to_string (m.w.cutoffHz) + " Hz), not the 4 kHz corner");
        ok (m.w.nearNyquist, "LR4: flagged as unattributable, being near Nyquist");
    }

    // ---------- 9. the sample grid: what the samples prove ----------
    {
        SourceForensics sf;
        auto p = defaults();
        p.fftOrder = 8;
        p.maxDistinctValues = 1 << 10;
        sf.setParams (p);
        ok (run (sf.prepare (kFs, 64, 1)), "grid: prepare");
        const float table[] = { 1.0f, 0.5f, 0.25f, -0.5f, 3.0f / 32768.0f, 1.0f / 32768.0f, 1.0f / 8388608.0f,
                                2.0f, 5.0f, 0.1f, 0.7f, 1.0f / 3.0f, 1.0e-30f,
                                std::numeric_limits<float>::denorm_min(),
                                3.0f * std::numeric_limits<float>::denorm_min(),
                                std::numeric_limits<float>::min() };
        int mismatches = 0;
        for (float v : table)
        {
            sf.reset();
            const float one[1] { v };
            const float* in[1] { one };
            if (! sf.process (in, 1, 1)) ++mismatches;
            if (sf.sampleGrid (0).gridExponent != oracleGridK (v)) ++mismatches;
        }
        ok (mismatches == 0, "grid: " + std::to_string (mismatches) + " of 16 known values disagree with the oracle");
        {
            std::mt19937 rng (99);
            std::uniform_int_distribution<std::uint32_t> ud (0u, 0x7F7FFFFFu);
            int bad = 0, seen = 0;
            for (int i = 0; i < 400; ++i)
            {
                const float v = std::bit_cast<float> (ud (rng));
                if (! std::isfinite (v) || core::exactlyEqual (v, 0.0f)) continue;
                sf.reset();
                const float one[1] { v };
                const float* in[1] { one };
                (void) sf.process (in, 1, 1);
                ++seen;
                if (sf.sampleGrid (0).gridExponent != oracleGridK (v)) ++bad;
            }
            ok (seen > 300, "grid: the random sweep actually ran (" + std::to_string (seen) + " values)");
            ok (bad == 0, "grid: " + std::to_string (bad) + " random floats disagree with the oracle");
        }
        {
            sf.reset();
            const float one[1] { std::numeric_limits<float>::denorm_min() };
            const float* in[1] { one };
            (void) sf.process (in, 1, 1);
            const auto g = sf.sampleGrid (0);
            ok (g.gridExponent == 149, "grid: the smallest subnormal answers k = 149, not a capped sentinel");
            ok (g.valid && ! g.pcmCompatible && g.reason == ForensicsReason::FinerThan24BitGrid,
                "grid: a measured grid, and NOT a PCM word length");
            ok (g.minExactPcmBits == 0, "grid: so no word length is published");
        }
    }
    {
        struct Depth { int bits; int wantK; };
        const Depth depths[] = { { 8, 7 }, { 16, 15 }, { 24, 23 } };
        for (const auto& d : depths)
        {
            auto x = wallFixture (12000.0, -80.0, 5, kFs, 4096);
            fillRange (x);
            quantise (x, d.bits);
            auto p = defaults();
            p.fftOrder = 8;
            p.maxDistinctValues = 1 << 20;
            const auto m = measure (x, p);
            const std::string tag = std::to_string (d.bits) + "-bit: ";
            ok (m.g.valid && m.g.pcmCompatible, tag + "a PCM grid holds it");
            ok (m.g.gridExponent == d.wantK, tag + "k = " + std::to_string (m.g.gridExponent)
                                           + ", expected " + std::to_string (d.wantK));
            ok (m.g.minExactPcmBits == d.bits, tag + "the shortest exact word is "
                                             + std::to_string (m.g.minExactPcmBits));
            ok (m.g.offGridSamples == 0 && m.g.firstOffGridSample < 0, tag + "no sample is off the grid");
            ok (m.g.reason == ForensicsReason::Ok, tag + "and the reading is published");
        }
    }
    {
        // THE TRAP: 16-bit dithered up into 24 must NOT read as 16 bits.
        auto x = wallFixture (12000.0, -80.0, 5, kFs, 4096);
        fillRange (x);
        quantise (x, 16);
        std::mt19937 rng (31);
        std::uniform_int_distribution<int> d (-1, 1);
        const double q24 = std::ldexp (1.0, 23);
        for (auto& v : x) v = (float) ((double) v + (double) d (rng) / q24);
        auto p = defaults();
        p.fftOrder = 8;
        p.maxDistinctValues = 1 << 20;
        const auto m = measure (x, p);
        ok (m.g.gridExponent == 23 && m.g.minExactPcmBits == 24,
            "dither: 16-bit dithered into 24 reads as <= 24 bits, never as 16");
    }
    {
        // THE OTHER TRAP: 32-bit PCM in float32 is NOT "32 bits" — it is "no <= 24-bit grid".
        auto x = wallFixture (12000.0, -80.0, 5, kFs, 4096);
        const double q32 = std::ldexp (1.0, 31);
        for (auto& v : x) v = (float) (std::round ((double) v * q32) / q32);
        auto p = defaults();
        p.fftOrder = 8;
        p.maxDistinctValues = 1 << 20;
        const auto m = measure (x, p);
        ok (m.g.valid && ! m.g.pcmCompatible, "pcm32: a 32-bit PCM stream has no <= 24-bit grid");
        ok (m.g.reason == ForensicsReason::FinerThan24BitGrid, "pcm32: and says exactly that");
        ok (m.g.minExactPcmBits == 0, "pcm32: no word length is invented — 32 least of all");
        ok (m.g.gridExponent > 23, "pcm32: the measured grid is still published (k = "
                                 + std::to_string (m.g.gridExponent) + ")");
        ok (m.g.offGridSamples > 0 && m.g.firstOffGridSample >= 0, "pcm32: with a count and a coordinate");
        // ...and the CONTAINER-versus-CONTENT boundary, measured rather than assumed. float32's top binade
        // [0.5, 1) is spaced 2^-24, so a 32-bit stream near full scale reads k = 24 — one bit finer than a
        // 24-bit grid, so still not a PCM word this names. What reads as 24-bit is a 32-bit CONTAINER
        // carrying 24-bit content, and as 16-bit one carrying 16-bit content.
        struct Codes { const char* name; int step; int wantK; bool wantPcm; int wantBits; };
        const Codes codes[] = { { "pcm32 near full scale", 1,     24, false,  0 },
                                { "pcm32 carrying 24-bit", 256,   23, true,  24 },
                                { "pcm32 carrying 16-bit", 65536, 15, true,  16 } };
        for (const auto& cd : codes)
        {
            std::vector<float> v (2048, 0.0f);
            std::mt19937 rng (7);
            std::uniform_int_distribution<int> off (0, 4095);
            for (std::size_t i = 0; i < v.size(); ++i)
                v[(std::size_t) i] = (float) ((std::ldexp (1.0, 30) + (double) (off (rng) * cd.step)) / q32);
            const auto mc = measure (v, p);
            const std::string tag = std::string (cd.name) + ": ";
            ok (mc.g.gridExponent == cd.wantK, tag + "k = " + std::to_string (mc.g.gridExponent)
                                             + ", expected " + std::to_string (cd.wantK));
            ok (mc.g.pcmCompatible == cd.wantPcm, tag + "and the PCM reading is withheld exactly when it must be");
            ok (mc.g.minExactPcmBits == cd.wantBits, tag + "word length " + std::to_string (mc.g.minExactPcmBits));
        }
    }
    {
        // one denormal in an otherwise 16-bit programme: the MAX is a one-sample statistic, and the
        // histogram is what keeps that honest
        auto x = wallFixture (12000.0, -80.0, 5, kFs, 4096);
        fillRange (x);
        quantise (x, 16);
        x[1234] = std::numeric_limits<float>::denorm_min();
        SourceForensics sf;
        auto p = defaults();
        p.fftOrder = 8;
        p.maxDistinctValues = 1 << 20;
        sf.setParams (p);
        ok (run (sf.prepare (kFs, 512, 1)), "denormal: prepare");
        const float* in[1] { x.data() };
        ok (run (sf.process (in, 1, (int) x.size())), "denormal: process");
        sf.finish();
        const auto g = sf.sampleGrid (0);
        ok (g.gridExponent == 149 && ! g.pcmCompatible,
            "denormal: one denormal really does put the stream off every PCM grid");
        ok (g.firstMaxGridSample == 1234, "denormal: and the coordinate where it did so is published");
        const std::int64_t* h = sf.gridExponentHistogram (0);
        ok (h != nullptr && h[SourceForensics::gridExponentBuckets() - 1] == 1,
            "denormal: the histogram holds exactly one sample in the finer-than-PCM bucket");
        ok (h[15] > 100, "denormal: while the programme itself is a spike at k = 15 ("
                       + std::to_string (h[15]) + " samples)");
        std::int64_t total = 0;
        for (int k = 0; k < SourceForensics::gridExponentBuckets(); ++k) total += h[(std::size_t) k];
        ok (total == g.nonZeroSamples, "denormal: and the histogram accounts for every grid witness");
    }
    {
        // a non-dyadic gain; exactly +-1.0; a sample past unity
        auto x = wallFixture (12000.0, -80.0, 5, kFs, 2048);
        fillRange (x);
        quantise (x, 16);
        for (auto& v : x) v = (float) (0.1 * (double) v);
        auto p = defaults();
        p.fftOrder = 8;
        const auto m = measure (x, p);
        ok (! m.g.pcmCompatible, "gain: a non-dyadic gain leaves no <= 24-bit PCM grid at all");

        // exactly +-1.0. A b-bit word carries i/2^(b-1) for i in [-2^(b-1), 2^(b-1) - 1], so -1.0 IS a PCM
        // sample and +1.0 is NOT, at any depth — however routinely a float master's peak sits there. The
        // grid is still measured; only the word length is withheld (code-review round: this used to
        // publish "1 bit", which no normalised word can hold).
        const float pm1[4] { 1.0f, -1.0f, 1.0f, -1.0f };
        SourceForensics one;
        one.setParams (p);
        ok (run (one.prepare (kFs, 64, 1)), "unity: prepare");
        const float* in1[1] { pm1 };
        ok (run (one.process (in1, 1, 4)), "unity: process");
        one.finish();
        const auto g1 = one.sampleGrid (0);
        ok (g1.valid && g1.outsidePcmRange && g1.reason == ForensicsReason::OutsidePcmRange,
            "unity: exactly +1.0 is outside every normalised PCM word's range");
        ok (g1.gridExponent == 0 && g1.minExactPcmBits == 0,
            "unity: the grid is measured, the word length withheld");
        approx (g1.absPeak, 1.0, 1e-12, "unity: with the peak published beside it");
        approx (g1.sampleMax, 1.0, 1e-12, "unity: and the signed range the claim rests on");
        approx (g1.sampleMin, -1.0, 1e-12, "unity: on both sides");
        // -1.0 alone, however, IS a PCM sample and keeps its word length
        const float minusOne[3] { -1.0f, -0.5f, 0.0f };
        SourceForensics neg;
        neg.setParams (p);
        ok (run (neg.prepare (kFs, 64, 1)), "minus one: prepare");
        const float* in3[1] { minusOne };
        ok (run (neg.process (in3, 1, 3)), "minus one: process");
        neg.finish();
        const auto g3 = neg.sampleGrid (0);
        ok (g3.pcmCompatible && ! g3.outsidePcmRange && g3.minExactPcmBits == 2,
            "minus one: -1.0 is inside the range and keeps its word length");

        const float over[4] { 0.5f, 1.5f, -0.25f, 0.0f };
        SourceForensics ovr;
        ovr.setParams (p);
        ok (run (ovr.prepare (kFs, 64, 1)), "over-unity: prepare");
        const float* in2[1] { over };
        ok (run (ovr.process (in2, 1, 4)), "over-unity: process");
        ovr.finish();
        const auto g2 = ovr.sampleGrid (0);
        ok (g2.valid && g2.outsidePcmRange && g2.reason == ForensicsReason::OutsidePcmRange,
            "over-unity: a sample past unity fits no normalised PCM word");
        ok (! g2.pcmCompatible, "over-unity: and pcmCompatible says so too, not only the reason");
        ok (g2.minExactPcmBits == 0 && g2.gridExponent == 2,     // -0.25 is the finest of the four
            "over-unity: the grid is still measured, the word length withheld");
        approx (g2.absPeak, 1.5, 1e-12, "over-unity: the peak is published, since the claim rests on it");
        ok (g2.zeroSamples == 1 && g2.nonZeroSamples == 3, "over-unity: zeros are counted apart from grid witnesses");
    }
    {
        // all-zero, and the -0.0 / +0.0 identity
        auto p = defaults();
        p.fftOrder = 8;
        const float z[6] { 0.0f, -0.0f, 0.0f, -0.0f, 0.0f, -0.0f };
        SourceForensics sf;
        sf.setParams (p);
        ok (run (sf.prepare (kFs, 64, 1)), "zeros: prepare");
        const float* in[1] { z };
        ok (run (sf.process (in, 1, 6)), "zeros: process");
        sf.finish();
        const auto g = sf.sampleGrid (0);
        ok (! g.valid && g.reason == ForensicsReason::NoNonZeroSample,
            "zeros: silence witnesses no grid, and says so rather than answering 1 bit");
        ok (g.zeroSamples == 6 && g.distinctValues == 1,
            "zeros: -0.0 is the same VALUE as +0.0 and is not counted twice");
        ok (g.minExactPcmBits == 0 && g.gridExponent == 0, "zeros: and no word length is published");
    }
    {
        // the distinct count against a std::set oracle, its exact capacity, and its overflow
        auto x = wallFixture (12000.0, -80.0, 5, kFs, 4096);
        fillRange (x);
        quantise (x, 12);
        std::set<std::uint32_t> want;
        for (float v : x)
        {
            const std::uint32_t b = std::bit_cast<std::uint32_t> (v);
            want.insert (b == 0x80000000u ? 0u : b);
        }
        auto p = defaults();
        p.fftOrder = 8;
        p.maxDistinctValues = 1 << 14;
        const auto m = measure (x, p);
        ok (m.g.distinctComplete, "distinct: the set fitted");
        ok (m.g.distinctValues == (std::int64_t) want.size(),
            "distinct: " + std::to_string (m.g.distinctValues) + " values, the std::set oracle says "
            + std::to_string (want.size()));
        {   // exactly at capacity: repeats past the limit must NOT mark it incomplete
            auto tp = p;
            tp.maxDistinctValues = (int) want.size();
            SourceForensics tight;
            tight.setParams (tp);
            ok (run (tight.prepare (kFs, 512, 1)), "distinct: tight prepare");
            const float* in[1] { x.data() };
            ok (run (tight.process (in, 1, (int) x.size())), "distinct: tight process");
            ok (run (tight.process (in, 1, (int) x.size())), "distinct: the same values again");
            tight.finish();
            ok (tight.sampleGrid (0).distinctComplete
                && tight.sampleGrid (0).distinctValues == (std::int64_t) want.size(),
                "distinct: a stream repeating exactly `limit` values stays exactly countable");
        }
        {   // one value too many: a lower bound, and every other counter keeps counting
            auto sp = p;
            sp.maxDistinctValues = 64;
            const auto sm = measure (x, sp);
            ok (! sm.g.distinctComplete && sm.g.distinctValues == 64,
                "distinct: exhaustion is DATA — 64 stored, and the set is marked incomplete");
            ok (sm.g.nonZeroSamples + sm.g.zeroSamples == (std::int64_t) x.size(),
                "distinct: and every sample is still counted past the overflow");
        }
        {   // 256 distinct values on a 16-bit grid — a mu-law-shaped source. The DEPTH reads 16 and the
            // DISTINCT COUNT is the actual answer; neither number is the other's substitute.
            std::vector<float> few (4096, 0.0f);
            std::mt19937 rng (43);
            std::uniform_int_distribution<int> pick (-128, 127);
            for (auto& v : few) v = (float) (std::ldexp ((double) (pick (rng) * 251), -15));
            auto fp = p;
            fp.maxDistinctValues = 1 << 12;
            const auto fm = measure (few, fp);
            ok (fm.g.pcmCompatible && fm.g.minExactPcmBits == 16,
                "mu-law-shaped: the grid reads 16 bits (" + std::to_string (fm.g.minExactPcmBits) + ")");
            ok (fm.g.distinctComplete && fm.g.distinctValues <= 256,
                "mu-law-shaped: while only " + std::to_string (fm.g.distinctValues) + " distinct values exist");
        }
    }

    // ---------- 10. LAW 8a: the full trace, bit-exact, across slicings ----------
    {
        const int order = 8, n = 1 << order;
        for (int hop : { n / 2, 97 })
        {
            auto x = ragged (4099, 2, 5);
            // a non-finite sample and an off-grid sample placed against the geometry: one inside a window
            // on a hop boundary, one just past a window end. Both are canonical events the trace must agree
            // on, whatever the call boundaries.
            x[0][(std::size_t) (2 * hop)] = std::numeric_limits<float>::quiet_NaN();
            x[1][(std::size_t) (n + 1)] = 0.1234567f;
            RunSpec base;
            base.order = order;
            base.hop = hop;
            base.slices = { 4099 };
            const auto ref = runTrace (x, base);
            ok (ref.size() > 4096, "8a hop=" + std::to_string (hop) + ": the trace is non-trivial ("
                                   + std::to_string (ref.size()) + " words)");
            const std::vector<std::vector<int>> slicings = {
                { 1 }, { 2 }, { 3 }, { 7 }, { 97 }, { hop - 1 }, { hop }, { hop + 1 },
                { n - 1 }, { n }, { n + 1 }, { 512 }, { 4097 }, { 8192 }, { 0, 1, 0, 700, 0 },
            };
            int differ = 0;
            for (const auto& sl : slicings)
            {
                RunSpec s = base;
                s.slices = sl;
                if (runTrace (x, s) != ref) ++differ;
            }
            ok (differ == 0, "8a hop=" + std::to_string (hop) + ": " + std::to_string (differ) + " of "
                             + std::to_string (slicings.size()) + " slicings differ");
            const std::vector<std::vector<int>> onEdges = {
                { 2 * hop, 1, 4099 },                                  // the non-finite sample alone
                { n, 1, 4099 },                                        // the first frame's last sample alone
                { n + 1, 1, 4099 },                                    // the off-grid sample alone
                { hop, hop, hop, hop, hop, hop, hop, hop, hop, hop, hop, hop, hop, hop, hop, hop },
            };
            int edgeDiffer = 0;
            for (const auto& sl : onEdges)
            {
                RunSpec s = base;
                s.slices = sl;
                if (runTrace (x, s) != ref) ++edgeDiffer;
            }
            ok (edgeDiffer == 0, "8a hop=" + std::to_string (hop) + ": " + std::to_string (edgeDiffer)
                                 + " of 4 boundary-on-transition partitions differ");
            int raggedDiffer = 0;
            for (unsigned seed = 0; seed < 6; ++seed)
            {
                std::mt19937 rng (seed + 71u);
                std::uniform_int_distribution<int> d (1, 1500);
                RunSpec s = base;
                s.slices.clear();
                for (int i = 0; i < 64; ++i) s.slices.push_back (d (rng));
                if (runTrace (x, s) != ref) ++raggedDiffer;
            }
            ok (raggedDiffer == 0, "8a hop=" + std::to_string (hop) + ": " + std::to_string (raggedDiffer)
                                   + " of 6 seeded ragged partitions differ");
            int mbDiffer = 0;
            for (int mb : { 1, 64, 4096, 65536 })
            {
                RunSpec s = base;
                s.maxBlock = mb;
                s.slices = { 333 };
                if (runTrace (x, s) != ref) ++mbDiffer;
            }
            ok (mbDiffer == 0, "8a hop=" + std::to_string (hop) + ": " + std::to_string (mbDiffer)
                               + " of 4 prepared maxBlock values differ");
            // a channel that disappears mid-stream. The presence mask is INPUT — a function of the absolute
            // sample index — so runTrace cuts any call that would straddle the transition.
            RunSpec dis = base;
            dis.narrowFrom = n + hop + 13;
            dis.slices = { 4099 };
            const auto disRef = runTrace (x, dis);
            int disDiffer = 0;
            for (const auto& sl : { std::vector<int> { 1 }, std::vector<int> { 3 }, std::vector<int> { hop },
                                    std::vector<int> { 997 }, std::vector<int> { 8192 } })
            {
                RunSpec s = dis;
                s.slices = sl;
                if (runTrace (x, s) != disRef) ++disDiffer;
            }
            ok (disDiffer == 0, "8a hop=" + std::to_string (hop) + ": a disappearing channel, "
                                + std::to_string (disDiffer) + " of 5 slicings differ");
            ok (disRef != ref, "8a hop=" + std::to_string (hop) + ": and its absence is VISIBLE in the report");
            // the distinct table's exhaustion is a canonical event too
            RunSpec tiny = base;
            tiny.distinct = 32;
            tiny.slices = { 4099 };
            const auto tinyRef = runTrace (x, tiny);
            int tinyDiffer = 0;
            for (const auto& sl : { std::vector<int> { 1 }, std::vector<int> { 5 }, std::vector<int> { 1024 } })
            {
                RunSpec s = tiny;
                s.slices = sl;
                if (runTrace (x, s) != tinyRef) ++tinyDiffer;
            }
            ok (tinyDiffer == 0, "8a hop=" + std::to_string (hop) + ": the table's exhaustion lands on the same "
                                 "sample in every slicing (" + std::to_string (tinyDiffer) + " of 3 differ)");
        }
    }

    // ---------- 11. the second oracle: the documented rules, re-implemented ----------
    {
        const double cuts[] = { 10000.0, 16000.0, 20500.0 };
        int coordDiffer = 0, foundDiffer = 0, geomDiffer = 0;
        double worstDrop = 0.0;
        for (double cut : cuts)
        {
            auto x = wallFixture (cut, cut > 20000.0 ? -99.0 : -80.0);
            SourceForensics sf;
            const auto p = defaults();
            sf.setParams (p);
            ok (run (sf.prepare (kFs, 512, 1)), "oracle2: prepare");
            const float* in[1] { x.data() };
            ok (run (sf.process (in, 1, (int) x.size())), "oracle2: process");
            sf.finish();
            std::vector<double> mean ((std::size_t) sf.bins());
            for (int b = 0; b < sf.bins(); ++b) mean[(std::size_t) b] = sf.meanPower (0, b);
            const auto o = oracleWall (mean, p, kFs);
            const auto w = sf.wall (0);
            if (o.cellCount != sf.cellCount() || o.binsPerCell != sf.binsPerCell()) ++geomDiffer;
            if (! core::exactlyEqual (o.cutoffHz, w.cutoffHz)) ++coordDiffer;
            if (! core::exactlyEqual (o.steepestHz, w.steepestHz)) ++coordDiffer;
            if (! core::exactlyEqual (o.endHz, w.transitionEndHz)) ++coordDiffer;
            if (o.found != w.valid) ++foundDiffer;
            worstDrop = std::max (worstDrop, std::fabs (o.dropDb - w.dropDb));
        }
        ok (geomDiffer == 0, "oracle2: the cell geometry agrees with the independent re-implementation");
        ok (coordDiffer == 0, "oracle2: " + std::to_string (coordDiffer) + " coordinates differ from it");
        ok (foundDiffer == 0, "oracle2: and so does the verdict on whether an edge exists");
        ok (worstDrop < 1e-9, "oracle2: the drop agrees to " + std::to_string (worstDrop) + " dB");
    }

    // ---------- 12. the tail, the short programme, holes, and the frame bookkeeping ----------
    {
        const int order = 8, n = 1 << order, hop = n / 2;
        auto p = defaults();
        p.fftOrder = order;
        p.hop = hop;
        const std::int64_t lengths[] { 0, 1, n - 1, n, n + hop - 1, n + hop, n + hop + 1 };
        for (std::int64_t total : lengths)
        {
            SourceForensics sf;
            sf.setParams (p);
            ok (run (sf.prepare (kFs, 64, 1)), "tail: prepare");
            std::vector<float> x ((std::size_t) std::max<std::int64_t> (total, 1), 0.25f);
            const float* in[1] { x.data() };
            if (total > 0) ok (run (sf.process (in, 1, (int) total)), "tail: process");
            sf.finish();
            const std::int64_t frames = total >= n ? (total - n) / hop + 1 : 0;
            const std::string tag = "tail T=" + std::to_string (total) + ": ";
            ok (sf.frames().frameCount() == frames, tag + std::to_string (sf.frames().frameCount())
                                                  + " frames, the oracle says " + std::to_string (frames));
            ok (sf.samplesProcessed() == total, tag + "every sample is on the clock");
            ok (sf.tailUncoveredSamples() == (frames > 0 ? total - ((frames - 1) * hop + n) : total),
                tag + "the uncovered tail is named exactly");
            const auto w = sf.wall (0);
            if (frames == 0)
                ok (! w.valid && w.reason == ForensicsReason::ShorterThanWindow,
                    tag + "no frame means ShorterThanWindow, not a zero that reads as 'nothing found'");
            else
                // A constant is not edgeless: all of its energy is at DC, so the window's own skirt IS a
                // band limit and the instrument finds it. What this section pins is the bookkeeping.
                ok (w.framesUsed == frames, tag + "every usable frame reached the accumulator");
            ok (sf.sampleGrid (0).nonZeroSamples == total, tag + "the grid family covers the uncovered tail too");
        }
    }
    {
        const int order = 8, n = 1 << order, hop = n / 2;
        auto p = defaults();
        p.fftOrder = order;
        p.hop = hop;
        SourceForensics sf;
        sf.setParams (p);
        ok (run (sf.prepare (kFs, 64, 1)), "holes: prepare");
        std::vector<float> x (1400, 0.5f);
        x[300] = std::numeric_limits<float>::infinity();
        const float* in[1] { x.data() };
        ok (run (sf.process (in, 1, 1400)), "holes: process");
        sf.finish();
        const std::int64_t frames = (1400 - n) / hop + 1;
        std::int64_t covering = 0;
        for (std::int64_t k = 0; k < frames; ++k)
            if (300 >= k * hop && 300 < k * hop + n) ++covering;
        ok (sf.framesHoled (0) == covering, "holes: " + std::to_string (sf.framesHoled (0)) + " frames refused, "
                                          + std::to_string (covering) + " windows contain the hole");
        ok (sf.framesUsed (0) + sf.framesHoled (0) == frames, "holes: and every frame is accounted for");
        const auto g = sf.sampleGrid (0);
        ok (g.nonFiniteSamples == 1 && g.nonZeroSamples == 1399,
            "holes: the non-finite sample is counted, never measured");
        ok (g.gridExponent == 1, "holes: and it does not touch the grid");
    }
    {
        auto p = defaults();
        p.fftOrder = 8;
        SourceForensics sf;
        sf.setParams (p);
        ok (run (sf.prepare (kFs, 64, 1)), "allholed: prepare");
        std::vector<float> x (1000, std::numeric_limits<float>::quiet_NaN());
        const float* in[1] { x.data() };
        ok (run (sf.process (in, 1, 1000)), "allholed: process");
        sf.finish();
        const auto w = sf.wall (0);
        ok (! w.valid && w.reason == ForensicsReason::AllFramesHoled, "allholed: says AllFramesHoled");
        ok (sf.framesUsed (0) == 0 && sf.framesHoled (0) > 0, "allholed: with the counts to prove it");
        ok (! sf.sampleGrid (0).valid && sf.sampleGrid (0).reason == ForensicsReason::NoNonZeroSample,
            "allholed: and the grid family has no witness either");
    }
    {
        // digital silence is FINITE, not a hole: frames exist and carry no power
        auto p = defaults();
        p.fftOrder = 8;
        SourceForensics sf;
        sf.setParams (p);
        ok (run (sf.prepare (kFs, 64, 1)), "silence: prepare");
        std::vector<float> x (1000, 0.0f);
        const float* in[1] { x.data() };
        ok (run (sf.process (in, 1, 1000)), "silence: process");
        sf.finish();
        const auto w = sf.wall (0);
        ok (sf.framesUsed (0) > 0, "silence: the frames are usable");
        ok (! w.valid && w.reason == ForensicsReason::NoSpectralEnergy,
            "silence: and the answer is NoSpectralEnergy, never a cutoff at the bottom of the band");
        ok (core::exactlyEqual (w.peakCellPower, 0.0) && ! w.emptyAboveValid,
            "silence: nothing is declared empty relative to a zero reference");
    }

    // ---------- 13. the aggregate is the mean of the per-channel MEANS ----------
    {
        // L is present throughout; R for one window only. Pooling the raw sums would weight L by its frame
        // count; the mean of means weights the two channels equally. The check is bit-exact, so the two
        // constructions cannot both pass it.
        const int order = 10, n = 1 << order, hop = n / 2;
        const std::size_t len = 8192;
        std::vector<std::vector<float>> x (2);
        x[0] = wallFixture (8000.0, -80.0, 3, kFs, len);
        x[1] = wallFixture (20000.0, -80.0, 4, kFs, len);
        auto p = defaults();
        p.fftOrder = order;
        p.hop = hop;
        SourceForensics sf;
        sf.setParams (p);
        ok (run (sf.prepare (kFs, 512, 2)), "aggregate: prepare");
        const float* both[2] { x[0].data(), x[1].data() };
        ok (run (sf.process (both, 2, n)), "aggregate: both channels for one window");
        const float* one[1] { x[0].data() + n };
        ok (run (sf.process (one, 1, (int) (len - (std::size_t) n))), "aggregate: then L alone");
        sf.finish();
        ok (sf.framesUsed (0) > sf.framesUsed (1) * 4, "aggregate: L has many more usable frames than R");
        ok (sf.framesUsed (1) >= 1, "aggregate: R still contributed one");
        ok (sf.wall (0).valid && sf.wall (1).valid, "aggregate: both channels measure an edge");
        ok (! core::exactlyEqual (sf.wall (0).cutoffHz, sf.wall (1).cutoffHz),
            "aggregate: and the two edges differ");
        // rebuild the aggregate independently from the published per-channel means and null its peak cell
        {
            const int bpc = sf.binsPerCell();
            double peak = 0.0;
            for (int j = 0; j < sf.cellCount(); ++j)
            {
                const int from = j * bpc, to = j == sf.cellCount() - 1 ? sf.bins() : from + bpc;
                double acc = 0.0;
                for (int b = from; b < to; ++b)
                {
                    double s = 0.0;
                    s += sf.meanPower (0, b);
                    s += sf.meanPower (1, b);
                    acc += s / 2.0;
                }
                peak = std::max (peak, acc / (double) (to - from));
            }
            ok (core::exactlyEqual (peak, sf.wall().peakCellPower),
                "aggregate: the mean of the per-channel means, bit for bit (pooled sums would differ)");
        }
        ok (sf.wall().framesUsed == sf.framesUsed (0) + sf.framesUsed (1),
            "aggregate: and it reports both channels' frames");
    }

    // ---------- 14. lifecycle, refusals, no allocation ----------
    {
        SourceForensics sf;
        auto p = defaults();
        p.fftOrder = 8;
        p.maxDistinctValues = 1 << 10;
        sf.setParams (p);
        const long before = g_allocs.load();
        ok (run (sf.prepare (kFs, 64, 2)), "life: prepare");
        ok (g_allocs.load() > before, "life: prepare is where the heap is touched");
        std::vector<float> a (1000, 0.3f), b (1000, -0.3f);
        const float* in[2] { a.data(), b.data() };
        const long b2 = g_allocs.load();
        ok (run (sf.process (in, 2, 1000)), "life: process");
        sf.finish();
        okNoAlloc (g_allocs.load() == b2, "life: process and finish allocate nothing");
        const long b3 = g_allocs.load();
        (void) sf.wall (0);
        (void) sf.wall();
        (void) sf.sampleGrid (0);
        okNoAlloc (g_allocs.load() == b3, "life: and reading the report allocates nothing");
        const auto after1 = sf.wall (0);
        sf.finish();
        ok (sf.isFinished(), "life: finish is idempotent");
        ok (core::exactlyEqual (after1.plateauPower, sf.wall (0).plateauPower)
            && core::exactlyEqual (after1.dropDb, sf.wall (0).dropDb),
            "life: and does not re-derive (a second normalisation would halve every power)");
        const std::int64_t frozen = sf.samplesProcessed();
        ok (! sf.process (in, 2, 10), "life: process refuses after finish");
        ok (sf.samplesProcessed() == frozen, "life: and consumes nothing");
        sf.reset();
        ok (sf.samplesProcessed() == 0 && ! sf.isFinished() && sf.framesUsed (0) == 0, "life: reset re-anchors");
        ok (! sf.sampleGrid (0).valid && sf.sampleGrid (0).distinctValues == 0, "life: including the grid family");
        ok (! sf.wall (0).valid && sf.wall (0).reason == ForensicsReason::NotFinished,
            "life: and the edge report is unarmed again");
        ok (run (sf.process (in, 2, 0)), "life: n == 0 is accepted");
        ok (sf.samplesProcessed() == 0, "life: and moves nothing");
        ok (run (sf.process (in, 2, 1000)), "life: a call 15x maxBlock is legal");
        ok (sf.samplesProcessed() == 1000, "life: and is consumed whole");
        ok (! sf.process (in, 3, 10), "refuse: more channels than prepared");
        ok (! sf.process (in, -1, 10), "refuse: a negative channel count");
        ok (! sf.process (in, 2, -1), "refuse: a negative length");
        ok (sf.samplesProcessed() == 1000, "refuse: and a refused call consumes nothing");
        ok (run (sf.process (in, 0, 5)), "life: a call carrying no channels is legal (all holes)");
        ok (sf.sampleGrid (0).absentSamples == 5, "life: and its samples are counted absent");
    }
    {
        // setParams() between prepare() and finish() must not move a published number: the measurement is
        // made with the snapshot prepare() took (code-review round — minDropDb was read at finish(), so
        // raising it afterwards turned a measured 79 dB edge from valid to invalid with no new prepare).
        auto x = wallFixture (12000.0, -80.0);
        SourceForensics sf;
        auto p = defaults();
        sf.setParams (p);
        ok (run (sf.prepare (kFs, 512, 1)), "snapshot: prepare");
        const float* in[1] { x.data() };
        ok (run (sf.process (in, 1, (int) x.size())), "snapshot: process");
        auto hostile = p;
        hostile.minDropDb = 400.0;
        hostile.nearNyquistFraction = 0.01;
        hostile.maxTransitionHz = 1.0;
        sf.setParams (hostile);                       // pending: must take effect only at the NEXT prepare
        sf.finish();
        const auto w = sf.wall (0);
        ok (w.valid && w.sharp, "snapshot: the edge is still the one the snapshot's thresholds accept");
        ok (! w.nearNyquist, "snapshot: and the pending nearNyquistFraction did not reach the report");
        approx (sf.activeParams().minDropDb, p.minDropDb, 1e-12, "snapshot: the active parameters are published");
        ok (core::exactlyEqual (sf.params().minDropDb, 400.0), "snapshot: while the pending ones are readable too");
        // ...and a FAILED prepare must disarm the report, not leave a valid one behind an unprepared object
        ok (! sf.prepare (0.0, 512, 1), "disarm: a bad prepare is refused");
        ok (! sf.isPrepared() && ! sf.isFinished(), "disarm: and the object is unarmed");
        ok (! sf.wall (0).valid && sf.wall (0).reason == ForensicsReason::NoChannel,
            "disarm: the channel report goes with it");
        ok (! sf.sampleGrid (0).valid && sf.sampleGrid (0).reason == ForensicsReason::NoChannel,
            "disarm: and so does the grid report");
    }
    {
        // WHICH nested edge becomes the primary is decided by the conservative drop and nothing else.
        // With a SHALLOW shelf the outer edge wins and the inner one is not reported at all — a second
        // search below the primary is structurally useless, since everything above such a candidate
        // includes the primary's own plateau. Documented, and pinned here (code-review round).
        auto x = spectral (kLen, kFs, 19, [] (double f)
                           { return f < 16000.0 ? 1.0 : (f < 20500.0 ? 1.0e-2 : 1.0e-6); });
        const auto m = measure (x, defaults());
        ok (m.w.valid, "shallow shelf: an edge is found");
        ok (m.w.cutoffHz > 20000.0, "shallow shelf: the OUTER edge wins on the conservative drop ("
                                    + std::to_string (m.w.cutoffHz) + " Hz)");
        ok (! m.w.secondValid, "shallow shelf: and nothing above it is reported as a second edge");
    }
    {
        // a replay after reset is bit-identical to a fresh object
        auto x = ragged (2000, 1, 19);
        RunSpec spec;
        spec.order = 8;
        spec.channels = 1;
        spec.slices = { 333 };
        const auto fresh = runTrace (x, spec);
        Trace tr;
        tr.w.reserve (1u << 18);
        SourceForensics sf;
        auto p = defaults();
        p.fftOrder = 8;
        p.maxDistinctValues = 1 << 10;
        sf.setParams (p);
        ok (run (sf.prepare (kFs, 512, 1)), "replay: prepare");
        const float* in[1] { x[0].data() };
        ok (run (sf.process (in, 1, 900)), "replay: a first, discarded pass");
        sf.finish();
        sf.reset();
        sf.setObserver (&observer, &tr);
        std::size_t at = 0;
        while (at < x[0].size())
        {
            const std::size_t take = std::min<std::size_t> (333, x[0].size() - at);
            const float* pl[1] { x[0].data() + at };
            ok (run (sf.process (pl, 1, (int) take)), "replay: process");
            at += take;
        }
        sf.finish();
        put (tr, sf.samplesProcessed());
        put (tr, sf.tailUncoveredSamples());
        put (tr, sf.frames().frameCount());
        putWall (tr, sf.wall (0));
        putGrid (tr, sf.sampleGrid (0));
        putWall (tr, sf.wall());
        ok (! fresh.empty() && tr.w == fresh, "replay: after reset the trace is bit-identical to a fresh object");
    }
    {
        // prepare() refusals. Every frequency parameter is bounded by the sample rate and every dB one by
        // kMaxDb, because a FINITE 1e300 narrowed to `int` is undefined behaviour, not a large cell.
        SourceForensics sf;
        auto p = defaults();
        p.fftOrder = SourceForensics::kMinFftOrder - 1;
        sf.setParams (p);
        ok (! sf.prepare (kFs, 64, 1), "refuse: fftOrder below the floor");
        p.fftOrder = SourceForensics::kMaxFftOrder + 1;
        sf.setParams (p);
        ok (! sf.prepare (kFs, 64, 1), "refuse: fftOrder above the ceiling");
        p = defaults();
        p.fftOrder = 8;
        p.hop = -1;
        sf.setParams (p);
        ok (! sf.prepare (kFs, 64, 1), "refuse: a negative hop");
        p.hop = 257;
        sf.setParams (p);
        ok (! sf.prepare (kFs, 64, 1), "refuse: a hop wider than the window");
        struct Bad { const char* what; SourceForensicsParams p; };
        std::vector<Bad> bads;
        auto mk = [] (auto fn) { SourceForensicsParams q; q.fftOrder = 8; fn (q); return q; };
        bads.push_back ({ "a distinct limit past the cap",
                          mk ([] (SourceForensicsParams& q) { q.maxDistinctValues = SourceForensics::kMaxDistinctLimit + 1; }) });
        bads.push_back ({ "a zero distinct limit", mk ([] (SourceForensicsParams& q) { q.maxDistinctValues = 0; }) });
        bads.push_back ({ "a zero cell width", mk ([] (SourceForensicsParams& q) { q.cellWidthHz = 0.0; }) });
        bads.push_back ({ "a NaN cell width",
                          mk ([] (SourceForensicsParams& q) { q.cellWidthHz = std::numeric_limits<double>::quiet_NaN(); }) });
        bads.push_back ({ "a finite but absurd 1e300 cell width",
                          mk ([] (SourceForensicsParams& q) { q.cellWidthHz = 1e300; }) });
        bads.push_back ({ "a plateau span past the sample rate",
                          mk ([] (SourceForensicsParams& q) { q.plateauSpanHz = 1e9; }) });
        bads.push_back ({ "a negative plateau span", mk ([] (SourceForensicsParams& q) { q.plateauSpanHz = -1.0; }) });
        bads.push_back ({ "a floor span of zero", mk ([] (SourceForensicsParams& q) { q.floorSpanHz = 0.0; }) });
        bads.push_back ({ "a search floor past the sample rate",
                          mk ([] (SourceForensicsParams& q) { q.searchFromHz = 1e300; }) });
        bads.push_back ({ "a near-Nyquist fraction above one",
                          mk ([] (SourceForensicsParams& q) { q.nearNyquistFraction = 1.5; }) });
        bads.push_back ({ "a near-Nyquist fraction of zero",
                          mk ([] (SourceForensicsParams& q) { q.nearNyquistFraction = 0.0; }) });
        bads.push_back ({ "a zero emptiness threshold", mk ([] (SourceForensicsParams& q) { q.emptyDb = 0.0; }) });
        bads.push_back ({ "an emptiness threshold past kMaxDb",
                          mk ([] (SourceForensicsParams& q) { q.emptyDb = 1e9; }) });
        bads.push_back ({ "a negative minimum drop", mk ([] (SourceForensicsParams& q) { q.minDropDb = -1.0; }) });
        bads.push_back ({ "an infinite minimum drop",
                          mk ([] (SourceForensicsParams& q) { q.minDropDb = std::numeric_limits<double>::infinity(); }) });
        bads.push_back ({ "a maximum transition of zero",
                          mk ([] (SourceForensicsParams& q) { q.maxTransitionHz = 0.0; }) });
        bads.push_back ({ "a negative transition threshold",
                          mk ([] (SourceForensicsParams& q) { q.transitionStartDb = -1.0; }) });
        bads.push_back ({ "an exemption rank past the cap",
                          mk ([] (SourceForensicsParams& q) { q.exemptCells = SourceForensics::kMaxExemptCells + 1; }) });
        bads.push_back ({ "a negative exemption rank", mk ([] (SourceForensicsParams& q) { q.exemptCells = -1; }) });
        bads.push_back ({ "a zero emptiness span", mk ([] (SourceForensicsParams& q) { q.emptyMinHz = 0.0; }) });
        int accepted = 0;
        for (const auto& bd : bads)
        {
            SourceForensics t;
            t.setParams (bd.p);
            if (t.prepare (kFs, 64, 1)) { ++accepted; ok (false, std::string ("refuse: ") + bd.what); }
            if (SourceForensics::storageFor (kFs, 1, bd.p).ok) ++accepted;
        }
        ok (accepted == 0, "refuse: all " + std::to_string (bads.size())
                           + " malformed parameter sets are refused by prepare AND by storageFor");
        p = defaults();
        p.fftOrder = 8;
        sf.setParams (p);
        ok (! sf.prepare (0.0, 64, 1), "refuse: a zero sample rate");
        ok (! sf.prepare (std::numeric_limits<double>::infinity(), 64, 1), "refuse: a non-finite sample rate");
        ok (! sf.prepare (999.0, 64, 1), "refuse: a sample rate below the floor");
        ok (! sf.prepare (kFs, 64, 0), "refuse: zero channels");
        ok (! sf.prepare (kFs, 64, core::kMaxChannels + 1), "refuse: too many channels");
        ok (! SourceForensics::storageFor (kFs, 0, p).ok, "refuse: storageFor agrees about the channel count");
        ok (SourceForensics::geometryFor (kFs, 1, p).ok, "refuse: and a good geometry is accepted");
        ok (run (sf.prepare (kFs, 64, 1)), "refuse: and a good call is accepted afterwards");
        SourceForensics raw;
        const float z[4] { 0.0f, 0.0f, 0.0f, 0.0f };
        const float* in[1] { z };
        ok (! raw.process (in, 1, 4), "refuse: an unprepared object refuses process");
        ok (! raw.wall (0).valid && ! raw.wall().valid, "refuse: and publishes nothing");
        raw.finish();
        ok (! raw.isFinished(), "refuse: finish on an unprepared object does nothing");
    }

    // ---------- 15. other sample rates, end to end ----------
    {
        for (double fs : { 8000.0, 96000.0 })
        {
            const double cut = fs * 0.25;                             // half of Nyquist: unambiguously below it
            auto x = wallFixture (cut, -80.0, 37, fs, 16384);
            auto p = defaults();
            p.fftOrder = 11;
            p.cellWidthHz = fs / 960.0;                               // the cells scale with the rate
            p.plateauSpanHz = fs / 24.0;
            p.floorSpanHz = fs / 24.0;
            p.searchFromHz = fs / 48.0;
            const auto m = measure (x, p, fs);
            const std::string tag = "fs=" + std::to_string ((int) fs) + ": ";
            ok (m.w.valid && m.w.sharp, tag + "the edge is found");
            ok (std::fabs (m.w.cutoffHz - cut) <= 4.0 * m.cellHz,
                tag + "at " + std::to_string (m.w.cutoffHz) + ", the construction says " + std::to_string (cut));
            ok (! m.w.nearNyquist, tag + "and half of Nyquist is not near it");
            ok (m.w.cutoffHz <= m.w.nyquistHz && m.w.emptyAboveHz <= m.w.nyquistHz,
                tag + "no coordinate exceeds Nyquist");
        }
    }

    return felitronics::test::report();
}

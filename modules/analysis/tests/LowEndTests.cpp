// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// analysis::LowEnd self-tests — the wide-bass and dominant-note instrument a lacquer asks for.
//
// WHAT IS NULLED AGAINST SOMETHING COMPUTED OUTSIDE THE OBJECT, which is the part that matters: a
// re-slicing test compares the implementation with ITSELF, so a schedule or a calibration that is
// deterministically wrong stays bit-identical across every slicing and passes green. This repository
// has already measured that (SpectrumFramesTests.cpp:31 — a frame trigger anchored to the wrong end
// survived all 14 slicings and took an outside witness to kill). So:
//   · THE LR4 ANALYTIC ORACLE. The settled side fraction of a "60 Hz mono bass + X Hz anti-phase tone"
//     fixture is predicted from the transfer function alone — |H_lp|^2 = 1/(1+r^4)^2 with the SVF's
//     PREWARPED r = tan(pi f/fs)/tan(pi fc/fs) — and compared with what the object measured, at
//     X = 120, 180, 240, 480 and 1000 Hz. The analogue-prototype r = f/fc is NOT this filter and is
//     visibly wrong here, so the oracle is a real second implementation, not a restatement.
//   · THE NEGATIVE TEST the task named: a programme WIDE above the crossover and MONO below must not
//     read as wide bass. Its bar is the number above, 4.8e-8, not "about zero".
//   · ABSOLUTE SPECTRAL CALIBRATION. A full-scale sine inside one semitone band must read A^2/2 — its
//     own mean square. That one number pins the one-sided FOLD, the |X|^2/(N*sum w^2) normalisation and
//     the half-bin cell offset at once; a constant calibration error survives every self-comparison.
//   · A DIRECT O(N^2) DFT BAND NULL. One frame's band energy against a windowed DFT and an
//     independently written fractional-overlap integration.
//   · THE CENTROID against a tone detuned by a known number of cents.
//
// LAW 8a — bit-identical under ARBITRARY re-slicing (stronger than law 11(a)'s same-boundaries promise).
// Compared BIT-EXACTLY with std::bit_cast, field by field, never memcmp (padding) and never != (which
// cannot tell -0.0 from +0.0) — the form of WaveformShapeTests.cpp:163:
//   · the whole TRACE, filled when a block CLOSES and when a frame is CONSUMED, not at the exit of
//     process(): coordinates, validity, hole counts and the RAW linear energies and per-frame band
//     powers, before any ratio, log or reduction;
//   · AND the whole final REPORT, field by field, because a trace does not cover the integrals, the
//     histogram, the extrema, the median or the centroid — and a defect in those cancels in a maximum,
//     a percentile or a rounding (docs/LAW8-KWEIGHTING.md:78: identical LUFS and dBTP over 5 changed
//     block energies out of 97);
//   · over 16 slicings — whole, 1, 2, 3, prime, B-1, B, B+1, H-1, H, H+1, W-1, W, W+1, larger than any
//     maxBlock, and seeded ragged partitions — with call boundaries landing deliberately on a hop, a
//     window end, a StateGrid tick (64) and a block edge; and across three prepared maxBlock values.
//
// MUTANT PASS — run once, not asserted, and recorded here because a gate that has never been red is not
// a gate. Each mutant was built with the object file DELETED first: a header edited in the same SECOND
// as the previous build is invisible to make's 1-second mtime granularity, and the stand refuses a
// verdict unless a `Building CXX` line appeared (this cost three false greens elsewhere today).
//   ( 1) frame consumed at the end of process() instead of on the sample that closed it -> RED, 24 of 322
//   ( 2) denormal flush moved to the end of process() instead of the StateGrid boundary  -> RED, 11 of 322
//   ( 3) the one-sided FOLD dropped (interior bins weighted 1 instead of 2)              -> RED,  7 of 322
//   ( 4) the half-bin cell offset dropped ([k,k+1) cells, not [k-1/2,k+1/2))             -> RED, 15 of 322
//   ( 5) the first moment taken at the BIN CENTRE instead of the overlap midpoint        -> RED,  6 of 322
//   ( 6) filter-then-encode instead of encode-then-filter                                -> RED,  1 of 322
//   ( 7) extrema ties broken by the LATEST block (> becomes >=)                          -> RED,  3 of 322
//   ( 8) the background median taken over ENERGIES instead of densities                  -> RED,  1 of 322
//   ( 9) the partial final block claims a FULL block length                              -> RED,  8 of 322
//   (10) Mid analysed alone, the Side axis dropped from the note spectrum                -> RED,  3 of 322
// THREE OF THESE SURVIVED THE FIRST PASS, and the holes they exposed are why three sections exist:
// (5) was green because the only cents assertions sat on a 13-bin-wide band where the edge cells carry
// almost nothing — the DFT null now nulls the first MOMENT too, at a band 1.2 bins wide; (6) was green
// because exact mono gives a bit-exact zero either way, so the discriminator had to be a NEAR-mono
// fixture (one sample, one ulp); (8) was green because the white-noise section proved the densities are
// flat without ever checking that backgroundDensity() is the median OF them, which an independent
// median now does. A gate that has never been red is not a gate, and three of these were not.
// Mutants (3), (4), (5), (8) and (10) are invisible to the re-slicing comparison BY CONSTRUCTION: they
// are deterministically wrong under every slicing, so all 16 partitions agree with each other and stay
// green. Every one of them is killed by an oracle computed outside the object, and by nothing else.

#include <felitronics_test.h>
#include <felitronics/analysis/LowEnd.h>

#include <atomic>
#include <bit>
#include <cmath>
#include <complex>
#include <cstdint>
#include <limits>
#include <algorithm>
#include <random>
#include <string>
#include <vector>

static std::atomic<long> g_allocs { 0 };
void* operator new (std::size_t s) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s); }
void  operator delete (void* p) noexcept { std::free (p); }
void  operator delete (void* p, std::size_t) noexcept { std::free (p); }

using namespace felitronics;
using analysis::LowEnd;
using analysis::LowEndParams;
using analysis::LowEndReason;
using analysis::LowEndTrace;

namespace
{

constexpr double kPi = core::kPi;

//==============================================================================
// --- the outside oracles ---

// The LR4 (two cascaded Butterworth Q=1/sqrt2 sections) low-pass POWER response of THIS filter: the
// Cytomic SVF prewarps its cutoff with tan(pi fc/fs) (Svf.h:61), so the ratio is a ratio of tangents.
// |H_lp|^2 = (1/(1+r^4))^2.  The analogue prototype r = f/fc gives visibly different numbers here.
double lr4LowPower (double f, double fc, double fs) noexcept
{
    const double r = std::tan (kPi * f / fs) / std::tan (kPi * fc / fs);
    const double g = 1.0 / (1.0 + r * r * r * r);
    return g * g;
}
double lr4HighPower (double f, double fc, double fs) noexcept
{
    const double r = std::tan (kPi * f / fs) / std::tan (kPi * fc / fs);
    const double ri = 1.0 / r;
    const double g = 1.0 / (1.0 + ri * ri * ri * ri);
    return g * g;
}

double noteHzOf (int midi, double tuning) noexcept { return tuning * std::exp2 ((double) (midi - 69) / 12.0); }

// A windowed direct DFT and an independently written fractional-overlap band integration: the second
// implementation of the fold, written from the definition rather than from the header.
void dftBand (const std::vector<float>& x, std::size_t from, int n, double fs,
              double centreHz, double semiUp, double& energy, double& centroidHz)
{
    std::vector<double> w ((std::size_t) n);
    double sumW2 = 0.0;
    for (int i = 0; i < n; ++i) { w[(std::size_t) i] = 0.5 - 0.5 * std::cos (2.0 * kPi * (double) i / (double) n); sumW2 += w[(std::size_t) i] * w[(std::size_t) i]; }
    const double binHz = fs / (double) n;
    const int bins = n / 2 + 1;
    const double lo = centreHz / semiUp, hi = centreHz * semiUp;
    const int ka = std::max (0, (int) std::floor (lo / binHz + 0.5));
    const int kb = std::min (bins - 1, (int) std::floor (hi / binHz + 0.5));
    double e = 0.0, mom = 0.0;
    for (int k = ka; k <= kb; ++k)
    {
        std::complex<double> acc { 0.0, 0.0 };
        for (int i = 0; i < n; ++i)
        {
            const double v = (double) x[from + (std::size_t) i] * w[(std::size_t) i];
            const double a = -2.0 * kPi * (double) k * (double) i / (double) n;
            acc += std::complex<double> (v * std::cos (a), v * std::sin (a));
        }
        const double p = (acc.real() * acc.real() + acc.imag() * acc.imag()) / ((double) n * sumW2);
        double cellLo = ((double) k - 0.5) * binHz, cellHi = ((double) k + 0.5) * binHz;
        if (cellLo < 0.0) cellLo = 0.0;
        if (cellHi > 0.5 * fs) cellHi = 0.5 * fs;
        const double ovLo = std::max (cellLo, lo), ovHi = std::min (cellHi, hi);
        const double ov = ovHi > ovLo ? ovHi - ovLo : 0.0;
        const double fold = (k == 0 || k == bins - 1) ? 1.0 : 2.0;
        const double contrib = p * fold * ov / (cellHi - cellLo);
        e += contrib;
        // the FIRST MOMENT of a partial cell sits at the midpoint of the OVERLAP. Weighting it by the
        // BIN CENTRE instead is the defect this null exists to catch: at a band only ~1.2 bins wide the
        // edge cells carry most of the weight and the centroid can leave its own band entirely.
        mom += contrib * 0.5 * (ovLo + ovHi);
    }
    energy = e;
    centroidHz = e > 0.0 ? mom / e : 0.0;
}

//==============================================================================
// --- fixtures ---

struct Stereo { std::vector<float> l, r; };

// Non-stationary and structured, with a fixed seed: a flat tone is a weak witness of invariance.
// Impulses are planted immediately before and after the hop, window, grid and block boundaries.
Stereo fixture (std::size_t n, unsigned seed, std::int64_t hop, std::int64_t win, std::int64_t block)
{
    std::mt19937 rng (seed);
    std::uniform_real_distribution<float> u (-1.0f, 1.0f);
    Stereo s;
    s.l.assign (n, 0.0f); s.r.assign (n, 0.0f);
    for (std::size_t i = 0; i < n; ++i)
    {
        const double t = (double) i / 6000.0;
        const double env = i < n / 11 ? 0.0                                 // lead silence
                         : i > n * 9 / 10 ? 0.01                            // a quiet tail
                         : (i % 4801 < 200 ? 1.0 : 0.3);                    // bursts against a bed
        const double bass = 0.6 * std::sin (2.0 * kPi * 82.41 * t) + 0.25 * std::sin (2.0 * kPi * 123.9 * t);
        const double mid  = 0.3 * std::sin (2.0 * kPi * 517.0 * t);
        s.l[i] = (float) (env * (bass + mid + 0.15 * (double) u (rng)));
        s.r[i] = (float) (env * (0.93 * bass - mid + 0.15 * (double) u (rng)));   // wide highs, near-mono bass
    }
    const std::int64_t marks[] = { hop, win, block, 64 };
    for (std::int64_t m : marks)
        for (std::int64_t d = -1; d <= 1; ++d)
        {
            const std::int64_t at = m + d;
            if (at >= 0 && at < (std::int64_t) n) { s.l[(std::size_t) at] += 0.8f; s.r[(std::size_t) at] -= 0.8f; }
        }
    return s;
}

//==============================================================================
// --- the law-8a witnesses: the trace, and the whole report ---

struct Collector
{
    std::vector<std::uint64_t> bits;
    static void fn (void* user, const LowEndTrace& t) { static_cast<Collector*> (user)->take (t); }
    void take (const LowEndTrace& t)
    {
        bits.push_back ((std::uint64_t) (t.kind == LowEndTrace::Kind::Frame ? 1 : 0));
        bits.push_back ((std::uint64_t) t.index);
        bits.push_back ((std::uint64_t) t.start);
        bits.push_back ((std::uint64_t) t.end);
        bits.push_back (t.valid ? 1u : 0u);
        bits.push_back ((std::uint64_t) t.finiteSamples);
        bits.push_back ((std::uint64_t) t.holes);
        bits.push_back (std::bit_cast<std::uint64_t> (t.midEnergy));
        bits.push_back (std::bit_cast<std::uint64_t> (t.sideEnergy));
        bits.push_back ((std::uint64_t) t.bandCount);
        for (int b = 0; b < t.bandCount; ++b)                       // the RAW per-frame band powers
        {
            bits.push_back (t.bandMid  != nullptr ? std::bit_cast<std::uint64_t> (t.bandMid[b])  : 0u);
            bits.push_back (t.bandSide != nullptr ? std::bit_cast<std::uint64_t> (t.bandSide[b]) : 0u);
        }
    }
};

void reportBits (const LowEnd& x, std::vector<std::uint64_t>& o)
{
    auto d = [&o] (double v) { o.push_back (std::bit_cast<std::uint64_t> (v)); };
    auto i = [&o] (std::int64_t v) { o.push_back ((std::uint64_t) v); };
    i ((std::int64_t) x.widthReason()); i ((std::int64_t) x.noteReason());
    i (x.samplesProcessed()); i (x.blockSamples()); i (x.analysedChannels());
    d (x.lowMidEnergy()); d (x.lowSideEnergy()); d (x.highMidEnergy()); d (x.highSideEnergy());
    d (x.rawMidEnergy()); d (x.rawSideEnergy()); d (x.lowBandEnergy());
    d (x.lowSideFraction()); d (x.highSideFraction()); d (x.rawSideFraction());
    i (x.finiteSamples()); i (x.holeSamples()); i (x.nonFiniteSamples()); i (x.absentSamples());
    i (x.filterNonFiniteSamples()); i (x.firstHoleSample()); i (x.lastHoleSample());
    i (x.blockCount()); i (x.storedBlockCount()); i (x.blocksComplete() ? 1 : 0);
    for (std::int64_t b = 0; b < x.storedBlockCount(); ++b)
    {
        const analysis::LowEndBlock r = x.block (b);
        i (r.index); i (r.start); i (r.samples); i (r.finiteSamples); i (r.holes);
        d (r.midEnergy); d (r.sideEnergy); i (r.valid ? 1 : 0); d (r.sideFraction()); d (r.energy());
    }
    for (int b = 0; b < LowEnd::kHistogramBins; ++b) i (x.histogram (b));
    i (x.histogramSamples());
    i (x.worstFractionBlock()); d (x.worstFraction()); d (x.worstFractionEnergy());
    i (x.peakEnergyBlock()); d (x.peakBlockEnergy()); d (x.peakEnergyBlockFraction());
    i (x.peakSideEnergyBlock()); d (x.peakBlockSideEnergy());
    d (x.peakLowSideAmplitude()); i (x.peakLowSideAmplitudeAt());
    i (x.usedFrames()); i (x.holedFrames()); i (x.tailUncoveredSamples());
    i (x.windowSamples()); i (x.hopSamples()); d (x.binHz());
    i (x.bandCount()); i (x.underResolvedBands());
    for (int b = 0; b < x.bandCount(); ++b)
    {
        const analysis::LowEndBand r = x.band (b);
        i (r.midi); d (r.centreHz); d (r.widthHz); d (r.binsPerBand);
        d (r.midEnergy); d (r.sideEnergy); d (r.energy); d (r.density);
        d (r.centroidHz); d (r.centsOffset); d (r.sideFraction());
    }
    i (x.peakBand()); i (x.peakDensityBand()); i (x.secondBand());
    d (x.backgroundDensity()); d (x.peakBandEnergy()); d (x.peakBandWidthHz());
    d (x.secondBandEnergy()); d (x.totalBandEnergy()); d (x.peakShare());
    i (x.peakMidi()); d (x.peakNoteHz()); d (x.peakCentroidHz()); d (x.peakCentsOffset());
    d (x.peakBandSideFraction());
}

// One run: feed the fixture in the given slicing, collect the trace AND the final report as bits.
std::vector<std::uint64_t> run (const Stereo& x, const LowEndParams& p, double fs, int nch,
                                const std::vector<int>& slices, int maxBlock)
{
    LowEnd le;
    le.setParams (p);
    Collector c;
    c.bits.reserve (1u << 16);
    std::vector<std::uint64_t> out;
    if (! le.prepare (fs, maxBlock, nch)) { out.push_back (0xDEADu); return out; }
    le.setTrace (&Collector::fn, &c);
    const float* planes[2] = { nullptr, nullptr };
    const std::size_t total = x.l.size();
    std::size_t at = 0, s = 0;
    while (at < total)
    {
        const int want = slices[s % slices.size()]; ++s;
        const std::size_t take = want <= 0 ? 0u : std::min ((std::size_t) want, total - at);
        planes[0] = x.l.data() + at;
        planes[1] = x.r.data() + at;
        if (! le.process (planes, nch, (int) take)) { out.push_back (0xBADu); return out; }
        at += take;
    }
    if (! le.finish()) { out.push_back (0xBAD2u); return out; }
    out = c.bits;
    reportBits (le, out);
    return out;
}

// Feeds a whole stereo buffer in one call and finishes. Used by every non-8a section.
bool feed (LowEnd& le, const Stereo& x, int nch)
{
    const float* planes[2] = { x.l.data(), x.r.data() };
    if (! le.process (planes, nch, (int) x.l.size())) return false;
    return le.finish();
}

// A settled two-tone fixture: `midHz` in phase (so it is pure Mid) and `sideHz` in ANTI-phase (pure
// Side), both at amplitude 1. This is the fixture the LR4 oracle predicts.
Stereo twoTone (std::size_t n, double fs, double midHz, double sideHz)
{
    Stereo s; s.l.assign (n, 0.0f); s.r.assign (n, 0.0f);
    for (std::size_t i = 0; i < n; ++i)
    {
        const double t = (double) i / fs;
        const double m = std::sin (2.0 * kPi * midHz * t);
        const double d = std::sin (2.0 * kPi * sideHz * t);
        s.l[i] = (float) (m + d);
        s.r[i] = (float) (m - d);
    }
    return s;
}

// THE CENTROID INVARIANT, checked wherever a band table exists: a band with energy must put its
// centroid strictly INSIDE its own [f*2^-1/24, f*2^+1/24), and a band without energy must report
// exactly 0. This is the property the overlap-midpoint first moment exists to guarantee — weighting a
// fractional edge cell by the BIN centre instead can push a centroid out of its own band — so it is
// asserted as an invariant rather than left to the one fixture that happens to notice.
int centroidsOutsideTheirBand (const LowEnd& le)
{
    int bad = 0;
    for (int b = 0; b < le.bandCount(); ++b)
    {
        const analysis::LowEndBand r = le.band (b);
        const double lo = r.centreHz * std::exp2 (-1.0 / 24.0), hi = r.centreHz * std::exp2 (1.0 / 24.0);
        if (! (r.energy > 0.0)) { if (r.centroidHz > 0.0) ++bad; continue; }   // an empty band claims nothing
        if (! (r.centroidHz >= lo && r.centroidHz <= hi)) ++bad;
        if (! std::isfinite (r.centroidHz) || ! std::isfinite (r.centsOffset)) ++bad;
        if (std::fabs (r.centsOffset) > 50.0 + 1e-6) ++bad;                    // a semitone is +-50 cents
    }
    return bad;
}

// Sums the LOW-band block energies over [fromBlock, toBlock), i.e. past the filter's startup.
void settledLow (const LowEnd& le, std::int64_t fromBlock, std::int64_t toBlock, double& mid, double& side)
{
    mid = 0.0; side = 0.0;
    for (std::int64_t b = fromBlock; b < toBlock && b < le.storedBlockCount(); ++b)
    {
        const analysis::LowEndBlock r = le.block (b);
        mid += r.midEnergy; side += r.sideEnergy;
    }
}

} // namespace

//==============================================================================
int main()
{
    using test::ok;
    using test::approx;

    // The test geometry: fs = 6000 with fftOrder 14 gives 0.3662 Hz bins — the SAME bin width, and so
    // the same semitone-vs-lobe geometry, as the default 2^17 at 48 kHz, for 1/8 of the samples.
    constexpr double kFs = 6000.0;
    LowEndParams base;
    base.fftOrder = 14;                 // N = 16384, 2.73 s at 6 kHz
    base.hop = 0;                       // N/2
    base.crossoverHz = 120.0;
    base.maxBlocks = 1 << 14;

    //==========================================================================
    test::group ("law 11d — storageFor is what prepare() allocates, and it refuses what prepare() refuses");
    {
        const LowEnd::Storage st = LowEnd::storageFor (kFs, 2, base);
        ok (st.ok, "storageFor accepts the default geometry");
        ok (st.bandCount == 40, "30..300 Hz at A4=440 is 40 semitone bands (MIDI 23..62), got " + std::to_string (st.bandCount));
        ok (st.blockSamples == 60, "10 ms at 6 kHz is lround(0.01*fs) = 60 samples, got " + std::to_string (st.blockSamples));
        ok (st.frames.ok && st.frames.bytes() > 0, "the nested SpectrumFrames ask is part of the budget");
        ok (st.bytes() > st.frames.bytes(), "the total exceeds the nested part");
        ok (st.accDoubles == 3u * (std::size_t) st.bandCount, "three per-band accumulators are budgeted");
        ok (st.binWeights > 0, "the precomputed fold table is budgeted");

        // the default parameters at 48 kHz: the nested ask alone is the number the header quotes
        const LowEnd::Storage big = LowEnd::storageFor (48000.0, 2, LowEndParams {});
        ok (big.ok && big.frames.bytes() == 5505040u,
            "the default 2^17 stereo frame store is 5505040 bytes, got " + std::to_string (big.frames.bytes()));

        // every refusal, and each one for its own reason
        LowEndParams p = base;
        ok (! LowEnd::storageFor (0.0, 2, p).ok, "sampleRate 0 refused");
        ok (! LowEnd::storageFor (std::numeric_limits<double>::quiet_NaN(), 2, p).ok, "NaN sampleRate refused");
        ok (! LowEnd::storageFor (kFs, 0, p).ok, "0 channels refused");
        ok (! LowEnd::storageFor (kFs, core::kMaxChannels + 1, p).ok, "too many channels refused");
        p = base; p.crossoverHz = std::numeric_limits<double>::quiet_NaN();
        ok (! LowEnd::storageFor (kFs, 2, p).ok, "a NaN crossover is refused — eq::Svf would take it into tan() and poison its coefficients permanently (Svf.h:44)");
        p = base; p.crossoverHz = 0.0;
        ok (! LowEnd::storageFor (kFs, 2, p).ok, "a crossover at 0 Hz refused");
        p = base; p.crossoverHz = 0.5 * kFs;
        ok (! LowEnd::storageFor (kFs, 2, p).ok, "a crossover above 0.49*fs refused rather than silently clamped");
        p = base; p.highNoteHz = 0.49 * kFs;
        ok (! LowEnd::storageFor (kFs, 2, p).ok, "a top band that would be clipped by Nyquist refused");
        p = base; p.lowNoteHz = 300.0; p.highNoteHz = 301.0;
        ok (! LowEnd::storageFor (kFs, 2, p).ok, "a range holding fewer than 2 notes refused");
        p = base; p.tuningHz = -440.0;
        ok (! LowEnd::storageFor (kFs, 2, p).ok, "a negative tuning refused");
        p = base; p.maxBlocks = -1;
        ok (! LowEnd::storageFor (kFs, 2, p).ok, "a negative capacity refused");
        p = base; p.fftOrder = 3;
        ok (! LowEnd::storageFor (kFs, 2, p).ok, "an fftOrder the frame producer refuses is refused here too");
        p = base; p.crossoverHz = 0.5;
        ok (! LowEnd::storageFor (kFs, 2, p).ok,
            "a crossover below 1 Hz is refused rather than accepted and then silently clamped to 1 Hz by eq::Svf");
        p = base; p.lowNoteHz = 1.0e-200; p.highNoteHz = 2.0e-200; p.tuningHz = 1.0e-200;
        ok (! LowEnd::storageFor (kFs, 2, p).ok,
            "a note range near zero is refused: its band-edge frequencies underflow the first moment to"
            " zero while the energies stay positive, and a centroid would then leave its own band");
        // the note-range rule is EXACT, not exact to within a rounding of its own logarithm
        {
            LowEndParams q = base; q.tuningHz = 440.0; q.lowNoteHz = 30.0;
            q.highNoteHz = std::nextafter (440.0, 0.0);
            const LowEnd::Storage st2 = LowEnd::storageFor (kFs, 2, q);
            ok (st2.ok, "a range ending one ulp below A4 is accepted");
            LowEnd le2; le2.setParams (q);
            ok (test::run (le2.prepare (kFs, 4096, 2)), "prepare");
            bool above = false;
            for (int b = 0; b < le2.bandCount(); ++b) if (le2.band (b).centreHz > q.highNoteHz) above = true;
            ok (! above, "and NO band centre lies above it — the integer bound is snapped against noteHz()"
                " itself, not left to the rounding of a log2");
        }
        // an OUTSIDE oracle on the weight table's SIZE: storageFor() and buildBands() reach the same
        // count through two different expressions, so the count is computed here a third way.
        {
            for (int order : { 4, 8, 12, 14 })
            {
                LowEndParams q = base; q.fftOrder = order;
                const LowEnd::Storage st3 = LowEnd::storageFor (kFs, 2, q);
                if (! st3.ok) continue;
                const std::int64_t nn = (std::int64_t) 1 << order;
                const double bh = kFs / (double) nn;
                const int bins = (int) (nn / 2 + 1);
                std::size_t want = 0;
                for (int b = 0; b < st3.bandCount; ++b)
                {
                    const double c = noteHzOf (23 + b, q.tuningHz);
                    const int ka = std::max (0, (int) std::floor (c * std::exp2 (-1.0 / 24.0) / bh + 0.5));
                    const int kb = std::min (bins - 1, (int) std::floor (c * std::exp2 (1.0 / 24.0) / bh + 0.5));
                    if (kb >= ka) want += (std::size_t) (kb - ka + 1);
                }
                ok (st3.binWeights == want, "the published weight-table size at order " + std::to_string (order)
                    + " matches an independently computed count (" + std::to_string (st3.binWeights)
                    + " vs " + std::to_string (want) + ")");
            }
        }

        // and prepare() refuses exactly the same arguments
        LowEnd le; le.setParams (base);
        ok (! le.prepare (0.0, 512, 2), "prepare refuses what storageFor refuses (rate)");
        ok (! le.isPrepared(), "a refused prepare leaves the object unprepared");
        const float z = 0.0f; const float* pl[2] = { &z, &z };
        ok (! le.process (pl, 2, 1), "process refuses before a successful prepare");
        ok (! le.finish(), "finish refuses before a successful prepare");
        ok (test::run (le.prepare (kFs, 512, 2)), "prepare accepts the default geometry");
    }

    //==========================================================================
    test::group ("the width matrix — silence, mono, hard-left, anti-phase");
    {
        const std::size_t n = 3000;                     // 0.5 s: 50 blocks, shorter than one window
        LowEndParams p = base;
        struct Case { const char* name; double lgain, rgain; double wantFrac; bool wantValid; };
        const Case cases[] = {
            { "mono (L = R): perfectly lateral",        1.0,  1.0, 0.0, true  },
            { "hard left (R = 0): half vertical",       1.0,  0.0, 0.5, true  },
            { "anti-phase (L = -R): pure vertical",     1.0, -1.0, 1.0, true  },
        };
        for (const Case& c : cases)
        {
            Stereo x; x.l.assign (n, 0.0f); x.r.assign (n, 0.0f);
            for (std::size_t i = 0; i < n; ++i)
            {
                const double v = std::sin (2.0 * kPi * 60.0 * (double) i / kFs);
                x.l[i] = (float) (c.lgain * v); x.r[i] = (float) (c.rgain * v);
            }
            LowEnd le; le.setParams (p);
            ok (test::run (le.prepare (kFs, 4096, 2)), "prepare");
            ok (test::run (feed (le, x, 2)), "feed");
            approx (le.lowSideFraction(), c.wantFrac, 2.0e-4, std::string ("side fraction — ") + c.name);
            ok (le.widthValid() == c.wantValid, std::string ("widthValid — ") + c.name);
        }
        // mono is EXACTLY zero, not nearly: 0.5f*(x-x) is exact, so no rounding can leak into Side
        {
            Stereo x; x.l.assign (n, 0.0f); x.r.assign (n, 0.0f);
            for (std::size_t i = 0; i < n; ++i)
            {
                const double v = std::sin (2.0 * kPi * 60.0 * (double) i / kFs);
                x.l[i] = (float) v; x.r[i] = (float) v;
            }
            LowEnd le; le.setParams (p);
            ok (test::run (le.prepare (kFs, 4096, 2)) && test::run (feed (le, x, 2)), "mono prepare+feed");
            ok (std::bit_cast<std::uint64_t> (le.lowSideEnergy()) == std::bit_cast<std::uint64_t> (0.0),
                "a mono programme's Side energy is bit-exactly +0.0 — not 'nearly zero'");
            ok (le.widthValid() && le.widthReason() == LowEndReason::Ok, "and it is VALID with zero width, not 'undefined'");
        }
        // digital silence is the 0/0, and the one place a zero would lie
        {
            Stereo x; x.l.assign (n, 0.0f); x.r.assign (n, 0.0f);
            LowEnd le; le.setParams (p);
            ok (test::run (le.prepare (kFs, 4096, 2)) && test::run (feed (le, x, 2)), "silence prepare+feed");
            ok (le.widthReason() == LowEndReason::NoEnergy, "digital silence reads NoEnergy, not a zero width");
            ok (! le.widthValid(), "and is not valid");
            ok (le.finiteSamples() == (std::int64_t) n, "…while every sample was still counted as finite");
        }
        // a MONO OBJECT: R is L by construction, so Side is identically zero and that is the right answer
        {
            std::vector<float> m (n, 0.0f);
            for (std::size_t i = 0; i < n; ++i) m[i] = (float) std::sin (2.0 * kPi * 60.0 * (double) i / kFs);
            LowEnd le; le.setParams (p);
            ok (test::run (le.prepare (kFs, 4096, 1)), "prepare for ONE channel");
            const float* pl[1] = { m.data() };
            ok (test::run (le.process (pl, 1, (int) n)) && test::run (le.finish()), "mono-object feed");
            ok (le.analysedChannels() == 1, "analysedChannels says 1");
            ok (std::bit_cast<std::uint64_t> (le.lowSideEnergy()) == std::bit_cast<std::uint64_t> (0.0),
                "a mono OBJECT's Side is bit-exactly zero");
            ok (le.widthValid(), "and valid");
        }
    }

    //==========================================================================
    test::group ("encode-then-filter — the precision the order was chosen for");
    {
        // Filtering L and R separately quantises two large correlated histories into float SVF state and
        // then subtracts them, MANUFACTURING side residue of order eps*|Mid| at every sample. Encoding
        // first forms the difference exactly (Sterbenz: two floats within a factor of 2 subtract exactly)
        // and only then filters it, so the Side axis keeps its RELATIVE precision.
        // The witness: a loud mono tone whose two channels differ at exactly ONE sample by ONE ulp. The
        // true Side signal is a single half-ulp impulse; anything far above that is the filter's own
        // rounding, i.e. the other order's noise floor. Exact mono cannot tell the two apart — both give
        // a bit-exact zero — which is why this fixture is near-mono and not mono.
        constexpr double fs = 48000.0;
        const std::size_t n = 96000;
        std::vector<float> L (n), R (n);
        for (std::size_t i = 0; i < n; ++i)
        {
            const float v = (float) std::sin (2.0 * kPi * 60.0 * (double) i / fs);
            L[i] = v; R[i] = v;
        }
        R[1234] = std::nextafter (R[1234], 2.0f);
        LowEndParams p = base; p.fftOrder = 12; p.crossoverHz = 120.0;
        LowEnd le; le.setParams (p);
        ok (test::run (le.prepare (fs, 8192, 2)), "prepare");
        const float* pl[2] = { L.data(), R.data() };
        ok (test::run (le.process (pl, 2, (int) n)) && test::run (le.finish()), "feed");
        const double ratio = le.lowSideEnergy() / le.lowMidEnergy();
        ok (le.lowSideEnergy() > 0.0, "the one-ulp difference IS measured, not lost");
        ok (ratio < 1.0e-18,
            "and the Side axis carries only that impulse: side/mid = " + std::to_string (ratio)
            + ", where per-sample filter rounding over 96000 samples would sit orders of magnitude higher");
    }

    //==========================================================================
    test::group ("the LR4 analytic oracle — the side fraction is PREDICTED, not just reproduced");
    {
        // 48 kHz here so the oracle's numbers are the ones the header quotes.
        constexpr double fs = 48000.0;
        LowEndParams p = base;
        p.fftOrder = 12;                         // the spectrum is not under test in this section
        p.crossoverHz = 120.0;
        const std::size_t n = 96000;             // 2 s = 200 blocks
        const double sideHz[] = { 120.0, 180.0, 240.0, 480.0, 1000.0 };
        for (double sh : sideHz)
        {
            const Stereo x = twoTone (n, fs, 60.0, sh);
            LowEnd le; le.setParams (p);
            ok (test::run (le.prepare (fs, 8192, 2)), "prepare");
            ok (test::run (feed (le, x, 2)), "feed");
            double mid = 0.0, side = 0.0;
            settledLow (le, 50, 200, mid, side);                        // from 0.5 s: past the LR4 transient
            const double got = side / (mid + side);
            const double want = lr4LowPower (sh, 120.0, fs)
                              / (lr4LowPower (60.0, 120.0, fs) + lr4LowPower (sh, 120.0, fs));
            // 1.5 % relative: the residual is the tones' non-integer periods over the window, not the filter
            approx (got / want, 1.0, 0.015,
                    "low side fraction at " + std::to_string ((int) sh) + " Hz vs the LR4 transfer function"
                    " (want " + std::to_string (want) + ", got " + std::to_string (got) + ")");
        }
        // AND THE ORACLE IS THE RIGHT FILTER. Low against Nyquist the prewarping barely shows (at
        // 240 Hz / 48 kHz the tangent form and the analogue prototype r = f/fc agree to 0.05 %), so the
        // table above would pass against either. Put the discriminator where prewarping bites — a high
        // crossover — and the two predictions are a factor of 5.7 apart: fc = 4 kHz with a 12 kHz side
        // tone is 2.650e-5 warped against 1.499e-4 for the prototype. Matching one REFUTES the other.
        {
            LowEndParams q = base; q.fftOrder = 12; q.crossoverHz = 4000.0;
            const std::size_t n = 96000;
            const Stereo x = twoTone (n, fs, 1000.0, 12000.0);
            LowEnd le; le.setParams (q);
            ok (test::run (le.prepare (fs, 8192, 2)) && test::run (feed (le, x, 2)), "prepare+feed at fc = 4 kHz");
            double mid = 0.0, side = 0.0;
            settledLow (le, 50, 200, mid, side);
            const double got = side / (mid + side);
            const double warped = lr4LowPower (12000.0, 4000.0, fs)
                                / (lr4LowPower (1000.0, 4000.0, fs) + lr4LowPower (12000.0, 4000.0, fs));
            auto protoPow = [] (double f, double fc) { const double g = 1.0 / (1.0 + std::pow (f / fc, 4.0)); return g * g; };
            const double prototype = protoPow (12000.0, 4000.0)
                                   / (protoPow (1000.0, 4000.0) + protoPow (12000.0, 4000.0));
            ok (warped / prototype < 0.25, "the two candidate oracles are a factor of "
                + std::to_string (prototype / warped) + " apart here, so this case decides between them");
            approx (got / warped, 1.0, 0.03, "the measurement follows the PREWARPED response (want "
                    + std::to_string (warped) + ", got " + std::to_string (got) + ")");
            ok (std::fabs (got / prototype - 1.0) > 0.5,
                "…and refutes the analogue prototype, so the oracle is this filter and not a restatement of it");
        }
    }

    //==========================================================================
    test::group ("the NEGATIVE test — wide above the crossover, mono below, must not read as wide bass");
    {
        constexpr double fs = 48000.0;
        LowEndParams p = base; p.fftOrder = 12; p.crossoverHz = 120.0;
        const std::size_t n = 96000;
        const Stereo x = twoTone (n, fs, 60.0, 1000.0);        // mono bass, anti-phase 1 kHz
        LowEnd le; le.setParams (p);
        ok (test::run (le.prepare (fs, 8192, 2)) && test::run (feed (le, x, 2)), "prepare+feed");
        double mid = 0.0, side = 0.0;
        settledLow (le, 50, 200, mid, side);
        const double got = side / (mid + side);
        const double want = lr4LowPower (1000.0, 120.0, fs)
                          / (lr4LowPower (60.0, 120.0, fs) + lr4LowPower (1000.0, 120.0, fs));
        ok (want < 1.0e-7, "the LR4 permits only " + std::to_string (want) + " here, so the bar is a NUMBER, not 'about zero'");
        ok (got < 1.0e-6, "the low band does NOT read as wide: " + std::to_string (got));
        ok (got > 0.0, "…and it is not exactly zero either — the leakage is real and measured, not suppressed");
        // the other half of the same claim: the crossover really did split
        const double wantHigh = lr4HighPower (1000.0, 120.0, fs)
                              / (lr4HighPower (60.0, 120.0, fs) + lr4HighPower (1000.0, 120.0, fs));
        ok (wantHigh > 0.99, "the LR4 high-pass oracle predicts a high-band side fraction of "
            + std::to_string (wantHigh));
        {
            LowEnd h; h.setParams (p);
            ok (test::run (h.prepare (fs, 8192, 2)) && test::run (feed (h, x, 2)), "a second pass for the high band");
            ok (h.highMidEnergy() > 0.0, "the high band's MID path is alive — without this the fraction"
                " reads exactly 1.0, which is within 0.35 % of the oracle and would pass a loose tolerance");
            const double wantRatio = lr4HighPower (1000.0, 120.0, fs) / lr4HighPower (60.0, 120.0, fs);
            approx ((h.highSideEnergy() / h.highMidEnergy()) / wantRatio, 1.0, 0.02,
                    "the HIGH band's side/MID RATIO follows its own oracle (want " + std::to_string (wantRatio)
                    + ", got " + std::to_string (h.highSideEnergy() / h.highMidEnergy())
                    + ") — a ratio a dead Mid path cannot fake, unlike the fraction");
            approx (h.highSideFraction() / wantHigh, 1.0, 0.002,
                    "and so does the fraction, at a tolerance tight enough to exclude 1.0");
        }
        ok (le.rawSideFraction() > 0.4 && le.rawSideFraction() < 0.6,
            "the UNFILTERED programme is half Side, as two equal-amplitude tones make it");
    }

    //==========================================================================
    test::group ("the crossover's own startup dominates a near-zero integral, and the series does not");
    {
        // A trap worth a gate rather than a comment: the LR4 starts at sample 0 with no history, and
        // while it charges it passes what the settled filter rejects. For a quantity as small as a
        // mono-bass side fraction that transient is nearly the whole integral.
        constexpr double fs = 48000.0;
        const std::size_t n = (std::size_t) (fs * 4.0);
        Stereo x; x.l.assign (n, 0.0f); x.r.assign (n, 0.0f);
        for (std::size_t i = 0; i < n; ++i)
        {
            const double t = (double) i / fs;
            const double b = 0.5 * std::sin (2.0 * kPi * 82.4069 * t);      // mono bass
            const double w = 0.3 * std::sin (2.0 * kPi * 900.0 * t);        // anti-phase highs
            x.l[i] = (float) (b + w); x.r[i] = (float) (b - w);
        }
        LowEndParams p = base; p.fftOrder = 12;
        LowEnd le; le.setParams (p);
        ok (test::run (le.prepare (fs, 8192, 2)) && test::run (feed (le, x, 2)), "prepare+feed");
        auto lp = [fs] (double f) {
            const double r = std::tan (kPi * f / fs) / std::tan (kPi * 120.0 / fs);
            const double g = 1.0 / (1.0 + r * r * r * r); return g * g; };
        const double pm = lp (82.4069) * 0.25 / 2.0, ps = lp (900.0) * 0.09 / 2.0;
        const double settledWant = ps / (pm + ps);
        // the SETTLED series matches the analytic prediction to four figures
        double mid = 0.0, side = 0.0;
        settledLow (le, 2, le.storedBlockCount(), mid, side);
        approx ((side / (mid + side)) / settledWant, 1.0, 2e-3,
                "from 20 ms on, the series matches the LR4 prediction (" + std::to_string (settledWant) + ")");
        // the WHOLE-FILE integral is an order of magnitude higher, and that is the startup, not a defect
        ok (le.lowSideFraction() > 20.0 * settledWant,
            "while the whole-file integral reads " + std::to_string (le.lowSideFraction())
            + " — over 20x the settled value, because the filter's charge-up is part of it");
        double s10 = 0.0;
        for (std::int64_t b = 0; b < 10; ++b) s10 += le.block (b).sideEnergy;
        ok (s10 / le.lowSideEnergy() > 0.9,
            "the first 100 ms holds " + std::to_string (100.0 * s10 / le.lowSideEnergy())
            + " % of the whole file's low SIDE energy — which is why the 10 ms series is published beside the integral");
    }

    //==========================================================================
    test::group ("absolute spectral calibration — a full-scale sine reads its own mean square");
    {
        // A tone exactly at the centre of the LOWEST band, where the semitone is only 4.87 bins wide.
        LowEndParams p = base;
        const double c0 = noteHzOf (23, 440.0);                    // B0, 30.8677 Hz
        const std::size_t n = 1u << 15;                            // two frames at hop N/2
        Stereo x; x.l.assign (n, 0.0f); x.r.assign (n, 0.0f);
        for (std::size_t i = 0; i < n; ++i)
        {
            const float v = (float) std::sin (2.0 * kPi * c0 * (double) i / kFs);
            x.l[i] = v; x.r[i] = v;                                 // mono, so Mid carries amplitude 1
        }
        LowEnd le; le.setParams (p);
        ok (test::run (le.prepare (kFs, 1 << 13, 2)) && test::run (feed (le, x, 2)), "prepare+feed");
        ok (le.noteValid(), "the note report is valid");
        const analysis::LowEndBand b0 = le.band (0);
        ok (b0.midi == 23, "band 0 is MIDI 23");
        approx (b0.centreHz, c0, 1e-9, "band 0's centre is the note");
        approx (b0.binsPerBand, 4.869, 0.01, "a semitone at 30.87 Hz is 4.87 bins at this geometry");
        // THE calibration number: a unit-amplitude sine's mean square is 0.5, and 99.96 % of it is in band
        approx (b0.midEnergy, 0.5, 0.002,
                "a full-scale sine inside one band reads A^2/2 = 0.5 (got " + std::to_string (b0.midEnergy) + ")");
        ok (b0.sideEnergy < 1e-30, "and nothing at all in the Side axis");
        ok (le.peakBand() == 0, "it is the peak band");
        approx (b0.centsOffset, 0.0, 0.5, "a tone at the centre reads 0 cents");
        ok (le.underResolvedBands() == 0, "no band is narrower than a Hann main lobe at this order");
        ok (centroidsOutsideTheirBand (le) == 0, "every band's centroid lies inside its own semitone");

        // the same tone one order LOWER does not resolve the bottom, and the report says so
        LowEndParams q = base; q.fftOrder = 12;
        LowEnd lo; lo.setParams (q);
        ok (test::run (lo.prepare (kFs, 1 << 13, 2)) && test::run (feed (lo, x, 2)), "prepare+feed at order 12");
        ok (lo.underResolvedBands() > 0,
            "at order 12 the low bands are narrower than the lobe, and underResolvedBands() names it: "
            + std::to_string (lo.underResolvedBands()));
    }

    //==========================================================================
    test::group ("the band fold against a direct DFT");
    {
        LowEndParams p = base;
        p.fftOrder = 12;                                    // N = 4096: one frame, and an O(N^2) DFT is affordable
        const std::size_t n = 1u << 12;                     // exactly one frame, [0, N)
        std::mt19937 rng (20260914u);
        std::uniform_real_distribution<float> u (-0.5f, 0.5f);
        Stereo x; x.l.assign (n, 0.0f); x.r.assign (n, 0.0f);
        for (std::size_t i = 0; i < n; ++i)
        {
            const double t = (double) i / kFs;
            const double v = 0.5 * std::sin (2.0 * kPi * 98.0 * t) + 0.2 * std::sin (2.0 * kPi * 171.3 * t) + 0.1 * (double) u (rng);
            x.l[i] = (float) v; x.r[i] = (float) v;          // mono: the Mid axis carries it all
        }
        LowEnd le; le.setParams (p);
        ok (test::run (le.prepare (kFs, 1 << 12, 2)) && test::run (feed (le, x, 2)), "prepare+feed");
        ok (le.usedFrames() == 1, "exactly one frame closed");
        const double semiUp = std::exp2 (1.0 / 24.0);
        int checked = 0;
        for (int b = 0; b < le.bandCount(); b += 7)          // a spread of bands: low, middle and top
        {
            const analysis::LowEndBand row = le.band (b);
            double want = 0.0, wantCentroid = 0.0;
            dftBand (x.l, 0, (int) n, kFs, row.centreHz, semiUp, want, wantCentroid);
            const double scale = std::max (want, 1e-18);
            approx (row.midEnergy / scale, 1.0, 2e-9,
                    "band " + std::to_string (b) + " (MIDI " + std::to_string (row.midi)
                    + ") against a direct windowed DFT and an independent fractional-overlap integration");
            // the CENTROID against the same outside oracle. At this order band 0 is only ~1.2 bins wide,
            // so its edge cells carry most of the weight and a moment taken at the bin centre instead of
            // the overlap midpoint is off by a large fraction of a bin.
            approx (row.centroidHz / std::max (wantCentroid, 1e-9), 1.0, 1e-9,
                    "band " + std::to_string (b) + "'s first MOMENT against the same oracle (want "
                    + std::to_string (wantCentroid) + " Hz, got " + std::to_string (row.centroidHz) + " Hz)");
            ++checked;
        }
        ok (centroidsOutsideTheirBand (le) == 0, "every band's centroid lies inside its own semitone");
        ok (checked >= 6, "the DFT null covered energy AND centroid over " + std::to_string (checked) + " bands");
    }

    //==========================================================================
    test::group ("the dominant note — placement, the runner-up, the centroid's cents");
    {
        LowEndParams p = base;
        const std::size_t n = 1u << 15;
        // E2 (MIDI 40) loud, A2 (MIDI 45) quieter: the peak and the runner-up are both known
        const double e2 = noteHzOf (40, 440.0), a2 = noteHzOf (45, 440.0);
        Stereo x; x.l.assign (n, 0.0f); x.r.assign (n, 0.0f);
        for (std::size_t i = 0; i < n; ++i)
        {
            const double t = (double) i / kFs;
            const float v = (float) (0.7 * std::sin (2.0 * kPi * e2 * t) + 0.2 * std::sin (2.0 * kPi * a2 * t));
            x.l[i] = v; x.r[i] = v;
        }
        LowEnd le; le.setParams (p);
        ok (test::run (le.prepare (kFs, 1 << 13, 2)) && test::run (feed (le, x, 2)), "prepare+feed");
        ok (le.peakMidi() == 40, "the peak band is E2 (MIDI 40), got MIDI " + std::to_string (le.peakMidi()));
        ok (std::string (LowEnd::pitchClassName (le.peakMidi())) == "E", "named E");
        ok (LowEnd::noteOctave (le.peakMidi()) == 2, "octave 2 (scientific pitch: MIDI 60 is C4)");
        ok (le.secondBand() >= 0 && le.band (le.secondBand()).midi == 45,
            "the runner-up is A2 (MIDI 45) — which is what says whether the 'dominant' note has a rival");
        approx (le.peakNoteHz(), e2, 1e-9, "the reported note frequency is the nominal note");
        ok (le.peakShare() > 0.8 && le.peakShare() <= 1.0, "the peak's share of the whole range is bounded and large");
        ok (le.peakBandSideFraction() < 1e-20, "the dominant note is entirely LATERAL — the cutting question");
        // dominance is derived from three published numbers, and none of them is a pole
        const double dom = le.peakBandEnergy() / (le.backgroundDensity() * le.peakBandWidthHz());
        ok (dom > 1e4, "the tone stands far above the background density (" + std::to_string (dom) + ")");

        // the CENTROID recovers a known detuning
        for (double cents : { -30.0, 20.0 })
        {
            Stereo y; y.l.assign (n, 0.0f); y.r.assign (n, 0.0f);
            const double f = e2 * std::exp2 (cents / 1200.0);
            for (std::size_t i = 0; i < n; ++i)
            {
                const float v = (float) std::sin (2.0 * kPi * f * (double) i / kFs);
                y.l[i] = v; y.r[i] = v;
            }
            LowEnd d; d.setParams (p);
            ok (test::run (d.prepare (kFs, 1 << 13, 2)) && test::run (feed (d, y, 2)), "prepare+feed");
            ok (d.peakMidi() == 40, "a tone " + std::to_string ((int) cents) + " cents off is still band E2");
            approx (d.peakCentsOffset(), cents, 2.0,
                    "the centroid recovers " + std::to_string ((int) cents) + " cents (got "
                    + std::to_string (d.peakCentsOffset()) + ")");
        }

        // AN ANTI-PHASE BASS NOTE MUST NOT VANISH. This is why the spectral axes are Mid AND Side: a
        // Mid-only note spectrum would report "no dominant note" for exactly the programme part 1 is
        // shouting about.
        {
            Stereo y; y.l.assign (n, 0.0f); y.r.assign (n, 0.0f);
            for (std::size_t i = 0; i < n; ++i)
            {
                const float v = (float) std::sin (2.0 * kPi * e2 * (double) i / kFs);
                y.l[i] = v; y.r[i] = -v;
            }
            LowEnd d; d.setParams (p);
            ok (test::run (d.prepare (kFs, 1 << 13, 2)) && test::run (feed (d, y, 2)), "prepare+feed");
            ok (d.noteValid(), "an anti-phase bass note is still a note");
            ok (d.peakMidi() == 40, "and it is found: MIDI " + std::to_string (d.peakMidi()));
            approx (d.peakBandSideFraction(), 1.0, 1e-6, "…reported as ENTIRELY vertical, which is the answer that matters");
            ok (d.band (d.peakBand()).midEnergy < 1e-30, "its Mid axis is empty");
            approx (d.lowSideFraction(), 1.0, 1e-6, "and part 1 agrees: the low end is pure Side");
        }

        // THE SEMITONE GRID'S TILT, measured: under a flat spectrum the top band holds 5.02 dB more
        // energy than the median band, with no note present. This is why the background is a median of
        // DENSITIES and why the density argmax is published beside the energy argmax.
        {
            std::mt19937 rng (7u);
            std::normal_distribution<float> g (0.0f, 0.25f);
            Stereo y; y.l.assign (n, 0.0f); y.r.assign (n, 0.0f);
            for (std::size_t i = 0; i < n; ++i) { const float v = g (rng); y.l[i] = v; y.r[i] = v; }
            LowEnd d; d.setParams (p);
            ok (test::run (d.prepare (kFs, 1 << 13, 2)) && test::run (feed (d, y, 2)), "prepare+feed white noise");
            const double top = d.band (d.bandCount() - 1).energy, bottom = d.band (0).energy;
            const double ratioDb = 10.0 * std::log10 (top / bottom);
            approx (ratioDb, 10.0, 2.0,
                    "under white noise the band ENERGIES tilt ~10 dB across a decade of the range (got "
                    + std::to_string (ratioDb) + " dB) — the grid's own tilt, not a note");
            const double dTop = d.band (d.bandCount() - 1).density, dBottom = d.band (0).density;
            approx (10.0 * std::log10 (dTop / dBottom), 0.0, 2.0, "while the DENSITIES are flat, which is why the background uses them");
            // and the background IS the median of those densities — computed here independently, over
            // the non-peak bands, with the same even-count convention
            std::vector<double> dens;
            for (int b = 0; b < d.bandCount(); ++b) if (b != d.peakBand()) dens.push_back (d.band (b).density);
            std::sort (dens.begin(), dens.end());
            const std::size_t k = dens.size();
            const double wantMedian = (k % 2) == 1 ? dens[k / 2] : 0.5 * (dens[k / 2 - 1] + dens[k / 2]);
            ok (std::bit_cast<std::uint64_t> (d.backgroundDensity()) == std::bit_cast<std::uint64_t> (wantMedian),
                "backgroundDensity() is bit-exactly the median DENSITY of the non-peak bands");
            std::vector<double> ener;
            for (int b = 0; b < d.bandCount(); ++b) if (b != d.peakBand()) ener.push_back (d.band (b).energy);
            std::sort (ener.begin(), ener.end());
            const double energyMedian = (k % 2) == 1 ? ener[k / 2] : 0.5 * (ener[k / 2 - 1] + ener[k / 2]);
            ok (std::fabs (energyMedian / wantMedian - 1.0) > 1.0,
                "…and that is a different number from the median ENERGY (" + std::to_string (energyMedian)
                + " against " + std::to_string (wantMedian) + "), so the distinction is under test");
        }
    }

    //==========================================================================
    test::group ("law 8a — the trace AND the whole report, bit-exact under 16 slicings and 3 maxBlocks");
    {
        LowEndParams p = base;
        p.fftOrder = 14;
        p.maxBlocks = 1 << 12;
        const std::int64_t W = 1 << 14, H = W / 2, B = 60;
        const std::size_t n = (std::size_t) (W + H + 4321);         // several frames, a ragged tail
        const Stereo x = fixture (n, 424242u, H, W, B);

        const std::vector<std::vector<int>> slicings = {
            { (int) n },                                            // the whole stream in one call
            { 1 }, { 2 }, { 3 }, { 7 }, { 97 },
            { (int) B - 1 }, { (int) B }, { (int) B + 1 },
            { (int) H - 1 }, { (int) H }, { (int) H + 1 },
            { (int) W - 1 }, { (int) W },
            { (int) W + 1 },                                        // larger than any maxBlock below
            { 64 },                                                 // every call boundary IS a StateGrid tick
            { 63, 1, 129, 5 },                                      // boundaries straddling the grid
            { 5, 1, 4000, 13, 1, 777, 2 },                          // ragged, fixed
            { (int) H, 1, (int) B, 3, (int) W, 11 },                // boundaries ON hop / block / window ends
        };
        const std::vector<std::uint64_t> want = run (x, p, kFs, 2, slicings[0], 1 << 13);
        ok (want.size() > 1000 && want[0] != 0xDEADu, "the reference run produced a trace of "
            + std::to_string (want.size()) + " words");
        int agreed = 0;
        for (std::size_t s = 1; s < slicings.size(); ++s)
        {
            const std::vector<std::uint64_t> got = run (x, p, kFs, 2, slicings[s], 1 << 13);
            bool same = got.size() == want.size();
            std::size_t firstDiff = 0;
            if (same)
                for (std::size_t i = 0; i < got.size(); ++i)
                    if (got[i] != want[i]) { same = false; firstDiff = i; break; }
            ok (same, "slicing " + std::to_string (s) + " is bit-identical"
                + (same ? "" : " (first difference at word " + std::to_string (firstDiff) + " of "
                                + std::to_string (want.size()) + ")"));
            if (same) ++agreed;
        }
        ok (agreed == (int) slicings.size() - 1, "all " + std::to_string (agreed) + " re-slicings agreed");

        // maxBlock sizes NOTHING: it selects no window, no hop, no reduction and no branch
        for (int mb : { 1, 64, 1 << 16 })
        {
            const std::vector<std::uint64_t> got = run (x, p, kFs, 2, slicings[17], mb);
            ok (got == want, "a run prepared with maxBlock " + std::to_string (mb) + " gives the same report");
        }
    }

    //==========================================================================
    test::group ("the tail, and the lengths around every boundary");
    {
        LowEndParams p = base; p.fftOrder = 12;                     // W = 4096, H = 2048
        const std::int64_t W = 1 << 12, H = W / 2, B = 60;
        for (std::int64_t T : { (std::int64_t) 0, (std::int64_t) 1, B - 1, B, B + 1,
                                W - 1, W, W + H - 1, W + H, W + H + 1 })
        {
            const Stereo x = fixture ((std::size_t) std::max<std::int64_t> (T, 1), 9u, H, W, B);
            LowEnd le; le.setParams (p);
            ok (test::run (le.prepare (kFs, 1 << 12, 2)), "prepare");
            if (T > 0)
            {
                const float* pl[2] = { x.l.data(), x.r.data() };
                ok (test::run (le.process (pl, 2, (int) T)), "process T=" + std::to_string (T));
            }
            ok (test::run (le.finish()), "finish T=" + std::to_string (T));
            const std::string tag = " (T=" + std::to_string (T) + ")";
            ok (le.samplesProcessed() == T, "every sample was consumed" + tag);
            // the TIME path has no tail: the blocks tile [0, T) exactly, the last one short
            const std::int64_t wantBlocks = (T + B - 1) / B;
            ok (le.blockCount() == wantBlocks, "blocks tile [0,T): want " + std::to_string (wantBlocks)
                + " got " + std::to_string (le.blockCount()) + tag);
            std::int64_t covered = 0;
            for (std::int64_t b = 0; b < le.storedBlockCount(); ++b) covered += le.block (b).samples;
            ok (covered == T, "and they cover exactly T samples, the last one short" + tag);
            ok (T % B != 0 || le.blockCount() * B == T, "no zero-length partial block when T is a whole number of blocks" + tag);
            // the SPECTRAL path has a tail, and names it rather than inventing a frame
            const std::int64_t wantFrames = T < W ? 0 : (T - W) / H + 1;
            ok (le.usedFrames() + le.holedFrames() == wantFrames,
                "frames are only the COMPLETE ones: want " + std::to_string (wantFrames) + tag);
            ok (le.tailUncoveredSamples() == (wantFrames == 0 ? T : T - ((wantFrames - 1) * H + W)),
                "the uncovered tail is NAMED" + tag);
            if (T < W)
            {
                ok (le.noteReason() == LowEndReason::ShorterThanWindow,
                    "a programme shorter than the window has no note report, with a reason" + tag);
                ok (le.peakBand() < 0, "and no peak band, rather than a zero that reads as 'nothing found'" + tag);
            }
            if (T == 0)
            {
                ok (le.widthReason() == LowEndReason::NoFiniteSamples, "an empty stream: NoFiniteSamples");
                ok (le.blockCount() == 0, "and no blocks at all");
            }
        }
    }

    //==========================================================================
    test::group ("holes — non-finite samples, a filter overflow, and a channel that disappears");
    {
        LowEndParams p = base; p.fftOrder = 12;
        const std::int64_t W = 1 << 12;
        const std::size_t n = (std::size_t) (2 * W);
        {
            Stereo x = fixture (n, 5u, W / 2, W, 60);
            x.l[100] = std::numeric_limits<float>::quiet_NaN();
            x.r[201] = std::numeric_limits<float>::infinity();
            x.l[3000] = -std::numeric_limits<float>::infinity();
            LowEnd le; le.setParams (p);
            ok (test::run (le.prepare (kFs, 1 << 13, 2)) && test::run (feed (le, x, 2)), "prepare+feed");
            ok (le.nonFiniteSamples() == 3, "three non-finite samples counted, got " + std::to_string (le.nonFiniteSamples()));
            ok (le.holeSamples() == 3, "and they are the three holes");
            ok (le.finiteSamples() == (std::int64_t) n - 3, "the rest reached the accumulators");
            ok (le.firstHoleSample() == 100 && le.lastHoleSample() == 3000, "the hole coordinates are published");
            ok (std::isfinite (le.lowSideEnergy()) && std::isfinite (le.lowMidEnergy()),
                "no non-finite value reached an accumulator");
            ok (le.widthValid(), "the width is still measured from what was usable");
            // the block holding a hole is marked, and its neighbours are not
            const analysis::LowEndBlock holed = le.block (100 / 60);
            ok (! holed.valid && holed.holes == 1, "the block holding the NaN is marked invalid with its hole count");
            ok (le.block (0).valid, "and block 0 is untouched");
            // a frame holding a hole is discarded whole, not guessed
            ok (le.holedFrames() > 0, "the frames containing a hole were discarded, not zero-filled: "
                + std::to_string (le.holedFrames()));
        }
        // a FILTER overflow from finite input: 3e38 is finite, its square is not, and neither is l+r
        {
            Stereo x; x.l.assign (2000, 0.0f); x.r.assign (2000, 0.0f);
            for (std::size_t i = 0; i < 2000; ++i) { x.l[i] = 0.1f; x.r[i] = 0.1f; }
            x.l[500] = 3.0e38f; x.r[500] = 3.0e38f;                 // l+r overflows float inside encode()
            LowEnd le; le.setParams (p);
            ok (test::run (le.prepare (kFs, 4096, 2)) && test::run (feed (le, x, 2)), "prepare+feed");
            ok (le.holeSamples() >= 1, "a finite input whose Mid/Side ENCODE overflows is a hole, counted");
            ok (std::isfinite (le.lowMidEnergy()) && std::isfinite (le.rawMidEnergy()),
                "and nothing non-finite reached any accumulator");
            ok (le.nonFiniteSamples() == 0, "the INPUT was finite, so it is not counted as a non-finite sample");
        }
        // PROMOTED BEFORE SQUARING, tested where it can actually fail. The fixture above never reaches
        // the squaring at all: l + r overflows inside encode, so the sample becomes a hole first. An
        // ANTI-PHASE 1.7e38 survives encode (m = 0 exactly, s = 0.5*(l-r) = 1.7e38, finite) and does get
        // squared — and in float (1.7e38)^2 is an infinity, so a square taken before the promotion to
        // double would publish one.
        {
            const std::size_t m = 600;
            Stereo x; x.l.assign (m, 0.0f); x.r.assign (m, 0.0f);
            for (std::size_t i = 200; i < 260; ++i) { x.l[i] = 1.7e38f; x.r[i] = -1.7e38f; }
            LowEnd le; le.setParams (p);
            ok (test::run (le.prepare (kFs, 4096, 2)) && test::run (feed (le, x, 2)), "prepare+feed");
            ok (le.holeSamples() == 0, "1.7e38 anti-phase survives the encode: not a hole");
            ok (le.lowSideEnergy() > 1.0e70, "the enormous side energy is measured, not clamped: "
                + std::to_string (le.lowSideEnergy()));
            ok (std::isfinite (le.lowSideEnergy()) && std::isfinite (le.rawSideEnergy())
                && std::isfinite (le.lowSideFraction()) && std::isfinite (le.peakLowSideAmplitude()),
                "and every published value is FINITE — the squares were taken in double, after the promotion");
            ok (le.widthValid(), "the measurement stands");
        }
        // a channel that disappears mid-stream is a hole for those samples, not a mono reading
        {
            Stereo x = fixture (1200, 6u, 2048, 4096, 60);
            LowEnd le; le.setParams (p);
            ok (test::run (le.prepare (kFs, 4096, 2)), "prepare");
            const float* both[2] = { x.l.data(), x.r.data() };
            const float* one[1] = { x.l.data() + 600 };
            ok (test::run (le.process (both, 2, 600)), "600 stereo samples");
            ok (test::run (le.process (one, 1, 600)), "then 600 with only the left channel");
            ok (test::run (le.finish()), "finish");
            ok (le.absentSamples() == 600, "the 600 samples missing their right channel are ABSENT, counted apart from non-finite: "
                + std::to_string (le.absentSamples()));
            ok (le.nonFiniteSamples() == 0, "and none of them is called a non-finite sample");
            ok (le.finiteSamples() == 600, "only the stereo half reached the accumulators — a missing R is NOT R = 0, which would read 50 % wide");
        }
    }

    //==========================================================================
    test::group ("capacity exhaustion is data — the prefix, the flag, and the counters that keep going");
    {
        LowEndParams p = base; p.fftOrder = 12; p.maxBlocks = 5;
        const std::size_t n = 1200;                                  // 20 blocks of 60
        const Stereo x = fixture (n, 11u, 2048, 4096, 60);
        LowEnd le; le.setParams (p);
        ok (test::run (le.prepare (kFs, 4096, 2)), "prepare with room for 5 blocks");
        const float* pl[2] = { x.l.data(), x.r.data() };
        ok (test::run (le.process (pl, 2, (int) n)), "a process() whose buffer fills MID-CALL is still accepted whole");
        ok (test::run (le.finish()), "finish");
        ok (le.blockCount() == 20, "every block was counted: " + std::to_string (le.blockCount()));
        ok (le.storedBlockCount() == 5, "five were stored");
        ok (! le.blocksComplete(), "and the report says the series is a prefix");
        ok (le.histogramSamples() > 5 * 60, "the histogram kept accumulating past the capacity: "
            + std::to_string (le.histogramSamples()));
        ok (le.worstFractionBlock() >= 0, "and so did the extrema");
        ok (le.widthValid(), "the measurement was not abandoned");
        // the report is IDENTICAL to a run with room to spare, apart from the series itself
        LowEndParams q = p; q.maxBlocks = 1 << 12;
        LowEnd big; big.setParams (q);
        ok (test::run (big.prepare (kFs, 4096, 2)) && test::run (feed (big, x, 2)), "the same stream with room");
        ok (std::bit_cast<std::uint64_t> (le.lowSideFraction()) == std::bit_cast<std::uint64_t> (big.lowSideFraction()),
            "the integral is bit-identical whether or not the series overflowed");
        ok (le.histogramSamples() == big.histogramSamples() && le.worstFractionBlock() == big.worstFractionBlock(),
            "and so are the histogram and the extrema — exhaustion changed the DATA KEPT, not the measurement");
    }

    //==========================================================================
    test::group ("the extrema are three different questions, and ties go to the earliest");
    {
        LowEndParams p = base; p.fftOrder = 12;
        // THREE REGIMES, 20 blocks each, because an LR4 RINGS ACROSS BLOCK BOUNDARIES: a loud block's
        // tail dominates the mid energy of the quiet block right after it, so a one-block regime does
        // not measure what it looks like it measures. The coordinates are read from the SETTLED middle.
        //   blocks  0..19  loud MONO          -> fraction 0,   energy 0.5  per sample
        //   blocks 20..39  quiet ANTI-PHASE   -> fraction 1,   energy 5e-5
        //   blocks 40..59  loud HALF-PANNED   -> fraction 0.5, energy 0.25, and the largest SIDE energy
        // so each of the three extrema is won by a different regime, by construction.
        const std::int64_t B = 60;
        const std::size_t n = (std::size_t) (60 * B);
        Stereo x; x.l.assign (n, 0.0f); x.r.assign (n, 0.0f);
        for (std::size_t i = 0; i < n; ++i)
        {
            const float v = (float) std::sin (2.0 * kPi * 40.0 * (double) i / kFs);
            const std::int64_t blk = (std::int64_t) i / B;
            if (blk < 20)      { x.l[i] = v;          x.r[i] = v; }
            else if (blk < 40) { x.l[i] = 0.01f * v;  x.r[i] = -0.01f * v; }
            else               { x.l[i] = v;          x.r[i] = 0.0f; }
        }
        LowEnd le; le.setParams (p);
        ok (test::run (le.prepare (kFs, 512, 2)) && test::run (feed (le, x, 2)), "prepare+feed");
        const std::int64_t wf = le.worstFractionBlock(), pe = le.peakEnergyBlock(), ps = le.peakSideEnergyBlock();
        ok (wf >= 22 && wf < 40, "the WORST FRACTION is in the quiet ANTI-PHASE regime (blocks 20..39), got "
            + std::to_string (wf));
        ok (le.worstFractionEnergy() < 0.01 * le.peakBlockEnergy(),
            "…and its energy is published beside it, so an accidental 1.0 can be weighed");
        ok (pe >= 2 && pe < 20, "the PEAK ENERGY block is in the loud MONO regime, got " + std::to_string (pe));
        ok (ps >= 42 && ps < 60, "the greatest VERTICAL modulation is in the loud HALF-PANNED regime, got "
            + std::to_string (ps));
        ok (wf != ps && pe != ps,
            "three different blocks — the cutting engineer's number has its own coordinate, which neither other extremum finds");
        ok (le.peakLowSideAmplitudeAt() >= 40 * B,
            "and the peak side AMPLITUDE lands in that regime too, at sample "
            + std::to_string (le.peakLowSideAmplitudeAt()));
        approx (le.block (30).sideFraction(), 1.0, 1e-6, "a settled anti-phase block reads fraction 1");
        approx (le.block (50).sideFraction(), 0.5, 1e-6, "a settled half-panned block reads fraction 0.5");
        approx (le.block (10).sideFraction(), 0.0, 1e-9, "a settled mono block reads fraction 0");

        // A TIE, EXACT BY CONSTRUCTION. Two identical blocks cannot be made by repeating a waveform —
        // the filter state differs at the two block starts and so do the energies to the last bit. What
        // CAN be made identical is a block that starts from an EXACTLY ZERO state: an anti-phase
        // impulse followed by enough silence for StateGrid's flush to zap the tail to exact zero
        // (Svf.h:161 flushes below 1e-15 on the grid boundary). The premise is ASSERTED before the rule
        // is, so this cannot become a blind fixture.
        {
            const std::size_t m = (std::size_t) (60 * B);
            Stereo y; y.l.assign (m, 0.0f); y.r.assign (m, 0.0f);
            for (int k = 0; k < 3; ++k)                                  // impulses 20 blocks apart
            {
                const std::size_t at = (std::size_t) (k * 20 * B);
                y.l[at] = 0.5f; y.r[at] = -0.5f;
            }
            LowEnd t; t.setParams (p);
            ok (test::run (t.prepare (kFs, 512, 2)) && test::run (feed (t, y, 2)), "prepare+feed the tie");
            const analysis::LowEndBlock b0 = t.block (0), b20 = t.block (20), b40 = t.block (40);
            const bool tied = std::bit_cast<std::uint64_t> (b0.sideEnergy) == std::bit_cast<std::uint64_t> (b20.sideEnergy)
                           && std::bit_cast<std::uint64_t> (b0.sideEnergy) == std::bit_cast<std::uint64_t> (b40.sideEnergy);
            ok (tied, "the premise: three blocks starting from an exactly zero filter state have BIT-IDENTICAL energies");
            ok (b0.sideEnergy > 0.0, "and they are not all zero");
            if (tied)
            {
                ok (t.peakSideEnergyBlock() == 0, "a tie goes to the EARLIEST block, got "
                    + std::to_string (t.peakSideEnergyBlock()));
                ok (t.worstFractionBlock() == 0, "and so does the worst-fraction tie");
                ok (t.peakEnergyBlock() == 0, "and the peak-energy tie");
            }
        }
    }

    //==========================================================================
    test::group ("the hostile corners of the parameter space that storageFor() still ACCEPTS");
    {
        const std::size_t n = 2400;
        const Stereo x = fixture (n, 13u, 2048, 4096, 60);
        // maxBlocks == 0: a legal capacity. Nothing is stored, everything is still counted, and no
        // accessor reads past an empty vector.
        {
            LowEndParams p = base; p.fftOrder = 12; p.maxBlocks = 0;
            ok (LowEnd::storageFor (kFs, 2, p).ok, "maxBlocks 0 is accepted");
            LowEnd le; le.setParams (p);
            ok (test::run (le.prepare (kFs, 4096, 2)) && test::run (feed (le, x, 2)), "prepare+feed");
            ok (le.blockCount() == 40 && le.storedBlockCount() == 0, "40 blocks counted, none stored");
            ok (! le.blocksComplete(), "and the series is reported as incomplete");
            ok (le.block (0).samples == 0 && le.block (-1).samples == 0 && le.block (1 << 20).samples == 0,
                "every block accessor is total — no read past an empty vector");
            ok (le.widthValid() && le.histogramSamples() > 0,
                "the measurement itself is unaffected: exhaustion changed the data KEPT, not the answer");
        }
        // the smallest window the frame producer allows, with hop == 1, so a frame closes on EVERY
        // sample past the first window and the 40 bands collapse onto one or two heavily clipped bins —
        // including bin 0, whose cell is a HALF cell and whose fold factor is 1, not 2.
        {
            LowEndParams p = base; p.fftOrder = 4; p.hop = 1;
            ok (LowEnd::storageFor (kFs, 2, p).ok, "fftOrder 4 with hop 1 is accepted");
            LowEnd le; le.setParams (p);
            ok (test::run (le.prepare (kFs, 4096, 2)) && test::run (feed (le, x, 2)), "prepare+feed");
            ok (le.usedFrames() + le.holedFrames() == (std::int64_t) n - 15,
                "a frame closed on every sample past the 16-sample window: "
                + std::to_string (le.usedFrames() + le.holedFrames()));
            ok (le.underResolvedBands() == le.bandCount(), "every band is under-resolved, and the report says so");
            int bad = 0;
            for (int b = 0; b < le.bandCount(); ++b)
            {
                const analysis::LowEndBand r = le.band (b);
                if (! (std::isfinite (r.energy) && r.energy >= 0.0 && std::isfinite (r.centroidHz)
                       && std::isfinite (r.centsOffset))) ++bad;
            }
            ok (bad == 0, "every one of the 40 band rows is finite and non-negative even here, got "
                + std::to_string (bad) + " that were not");
            ok (centroidsOutsideTheirBand (le) == 0,
                "and every centroid is still inside its own semitone, at a 16-point window where the bands"
                " collapse onto one or two heavily clipped bins — including bin 0's HALF cell");
            ok (le.noteValid(), "and the note report is produced rather than refused");
        }
        // hop == N (no overlap at all), fed through the widest channel count the core allows
        {
            LowEndParams p = base; p.fftOrder = 12; p.hop = 1 << 12;
            LowEnd le; le.setParams (p);
            ok (test::run (le.prepare (kFs, 4096, core::kMaxChannels)), "prepare for kMaxChannels");
            std::vector<std::vector<float>> ch ((std::size_t) core::kMaxChannels, std::vector<float> (n, 0.0f));
            for (std::size_t i = 0; i < n; ++i)
                for (int c = 0; c < core::kMaxChannels; ++c)
                    ch[(std::size_t) c][i] = (float) (0.3 * std::sin (2.0 * kPi * (60.0 + 7.0 * (double) c) * (double) i / kFs));
            std::vector<const float*> pp ((std::size_t) core::kMaxChannels);
            for (int c = 0; c < core::kMaxChannels; ++c) pp[(std::size_t) c] = ch[(std::size_t) c].data();
            ok (test::run (le.process (pp.data(), core::kMaxChannels, (int) n)) && test::run (le.finish()), "feed 16 channels");
            ok (le.analysedChannels() == 2, "only the L/R pair was analysed — channels 2.. take no part");
            ok (le.usedFrames() + le.holedFrames() == 0,
                "2400 samples cannot fill a 4096 window, so no frame closed at all");
            ok (le.noteReason() == LowEndReason::ShorterThanWindow, "…and the note report says exactly that");
            ok (le.widthValid(), "while the TIME path measured every sample, as it has no window");
        }
    }

    //==========================================================================
    test::group ("the centroid invariant across the geometries most likely to break it");
    {
        // Constant DC (all the energy at 0 Hz, i.e. OUTSIDE every band, pushing every band onto its own
        // lower edge), at five window orders and three sample rates — 1 kHz is the lowest the instrument
        // accepts, where 300 Hz is a third of Nyquist and the bands are enormous in bin terms.
        int bad = 0, cases = 0;
        for (double fs : { 1000.0, 8000.0, 48000.0 })
            for (int order : { 4, 6, 8, 12 })
            {
                LowEndParams p = base; p.fftOrder = order;
                if (! LowEnd::storageFor (fs, 2, p).ok) continue;
                LowEnd le; le.setParams (p);
                if (! le.prepare (fs, 4096, 2)) { ++bad; continue; }
                const std::size_t n = (std::size_t) (1u << order) * 3u;
                Stereo x; x.l.assign (n, 1.0f); x.r.assign (n, 1.0f);      // pure DC, mono
                if (! feed (le, x, 2)) { ++bad; continue; }
                bad += centroidsOutsideTheirBand (le);
                ++cases;
            }
        ok (cases >= 10, "the sweep covered " + std::to_string (cases) + " geometries");
        ok (bad == 0, "no band with energy puts its centroid outside its own semitone, and no empty band"
            " claims a frequency: " + std::to_string (bad) + " violations");
    }

    test::group ("finish() is idempotent, process() refuses after it, reset() replays identically");
    {
        LowEndParams p = base; p.fftOrder = 12;
        const std::size_t n = 1u << 13;
        const Stereo x = fixture (n, 31u, 2048, 4096, 60);
        LowEnd le; le.setParams (p);
        ok (test::run (le.prepare (kFs, 4096, 2)) && test::run (feed (le, x, 2)), "prepare+feed");
        std::vector<std::uint64_t> first;
        reportBits (le, first);
        ok (test::run (le.finish()), "a second finish() is accepted");
        std::vector<std::uint64_t> again;
        reportBits (le, again);
        ok (first == again, "…and changes not one bit of the report");
        const float* pl[2] = { x.l.data(), x.r.data() };
        ok (! le.process (pl, 2, 10), "process() refuses after finish()");
        ok (! le.process (pl, 2, 0), "…and so does an n == 0 call: after finish() EVERYTHING is refused until reset(), which is the law 11 order analysis::ClipDetector uses (the finished check precedes the n == 0 check)");
        le.reset();
        ok (! le.isFinished() && le.samplesProcessed() == 0, "reset() re-anchors the clock");
        ok (test::run (feed (le, x, 2)), "and the same stream can be replayed");
        std::vector<std::uint64_t> replay;
        reportBits (le, replay);
        ok (replay == first, "the replay is bit-identical — reset() re-anchored EVERYTHING: clock, ring, filters, grid, histogram, extrema, reasons");
    }

    //==========================================================================
    test::group ("no allocation in process() or finish()");
    {
        LowEndParams p = base; p.fftOrder = 12;
        const std::size_t n = 1u << 13;
        const Stereo x = fixture (n, 77u, 2048, 4096, 60);
        LowEnd le; le.setParams (p);
        ok (test::run (le.prepare (kFs, 1024, 2)), "prepare (this is the one call that allocates)");
        const float* pl[2] = { x.l.data(), x.r.data() };
        (void) pl;
        bool allAccepted = true;
        std::int64_t fed = 0;
        const long before = g_allocs.load();
        for (std::size_t at = 0; at < n; at += 997)
        {
            const float* q[2] = { x.l.data() + at, x.r.data() + at };
            const int take = (int) std::min<std::size_t> (997, n - at);
            allAccepted = le.process (q, 2, take) && allAccepted;
            fed += take;
        }
        allAccepted = le.finish() && allAccepted;
        const long after = g_allocs.load();
        // asserted AFTER the counter is read, so the assertion cannot allocate inside the measured region
        ok (allAccepted, "every call was ACCEPTED — without this a stage that had silently stopped"
            " processing would allocate nothing and pass this group");
        ok (le.samplesProcessed() == fed, "and it really consumed all " + std::to_string (fed) + " samples");
        ok (le.usedFrames() > 0, "…and really transformed frames while being measured");
        test::okNoAlloc (after == before, "process() and finish() allocated nothing ("
                         + std::to_string (after - before) + " allocations)");
    }

    return test::report();
}

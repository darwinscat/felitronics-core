// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// felitronics::analysis::MultiResSpectrumPaneFast — a paired reference NULL against the pane it is a
// sibling of. MultiResSpectrumPane is already property-tested against the physics
// (MultiResSpectrumTests.cpp, docs/ANALYZER-MULTIRES.md); this file does not re-derive any of it. It
// asserts that the fast pane, fed the SAME frames tick for tick, produces the same numbers — so all of
// that coverage transfers — and then pins the three places where the two are allowed to differ:
//
//   • the FILL is bit-identical, everywhere, always. Same accumulate, same float smoothing order, same
//     prefix sums, and buildColumns' column plan evaluates the band in the sibling's own operation
//     order. A single differing bit here is a bug, not a tolerance.
//   • the PEAK differs by the sibling's own dB rounding. It stores the hold as a float in dB and
//     reconstructs the power with an exp; the fast pane keeps power throughout. The divergence is the
//     OLD pane's quantisation (float ulp at −100 dB is 7.6e-6 dB there, 3.0e-7 dB in power), so the
//     budget is 1e-3 dB — anything looser would be hiding something.
//   • a NEGATIVE peakFallDb grows the hold without limit. In dB that merely got large; in power it
//     would reach an infinity, so the fast pane clamps at kMaxPower. Not a supported setting; pinned
//     so the divergence is deliberate rather than discovered.
//
// The scenarios are the ones the design consilium said have teeth — above all the LOUD NEIGHBOUR case,
// which is the one that fails if the prefix sums ever stop being double, and which a noise-only null
// would pass while shipping a broken pane.

#include <felitronics_test.h>
#include <felitronics/analysis/MultiResSpectrumPane.h>
#include <felitronics/analysis/MultiResSpectrumPaneFast.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

static std::atomic<long> g_allocs { 0 };
void* operator new      (std::size_t s) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s ? s : 1); }
void* operator new[]    (std::size_t s) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s ? s : 1); }
void  operator delete   (void* p) noexcept { std::free (p); }
void  operator delete[] (void* p) noexcept { std::free (p); }
void  operator delete   (void* p, std::size_t) noexcept { std::free (p); }
void  operator delete[] (void* p, std::size_t) noexcept { std::free (p); }

using namespace felitronics;
using felitronics::test::ok;
using felitronics::test::approx;
using felitronics::test::group;

using Ref  = analysis::MultiResSpectrumPane;
using Fast = analysis::MultiResSpectrumPaneFast;

namespace
{
constexpr double kPi = 3.14159265358979323846;
constexpr double fs  = 48000.0;
constexpr int    N   = Ref::kMaxSize;

struct Rng
{
    std::uint64_t s = 0x9E3779B97F4A7C15ull;
    double uni() noexcept { s ^= s >> 12; s ^= s << 25; s ^= s >> 27; return (double) ((s * 0x2545F4914F6CDD1Dull) >> 40) / 8388608.0 - 1.0; }
};

// A pair of panes that always see byte-identical frames, and the worst divergence seen so far.
struct Pair
{
    std::unique_ptr<Ref>  ref  { std::make_unique<Ref>() };
    std::unique_ptr<Fast> fast { std::make_unique<Fast>() };

    double worstFillDb = 0.0, worstPeakDb = 0.0, worstFillPow = 0.0;
    double deepestLifted = -1e9;                                          // worst peak read the fast pane shows where the sibling is below −150 dB
    double worstColFill = 0.0, worstColPeak = 0.0;
    bool   fillBitExact = true, peakAtLeastFill = true;

    void tune (float smooth, int cover) { ref->smoothCoeff = smooth; fast->smoothCoeff = smooth; ref->coverSamples = cover; fast->coverSamples = cover; }
    void reset() { ref->reset(); fast->reset(); }

    template <class Fill>
    void tick (Fill&& fill, int order = 14)
    {
        std::vector<float> buf ((std::size_t) (1 << order));
        fill (buf.data(), 1 << order);
        std::copy (buf.begin(), buf.end(), ref->frameInput());
        std::copy (buf.begin(), buf.end(), fast->frameInput());
        ref->ingest (order); fast->ingest (order);
    }
    void starveBoth() { ref->starve(); fast->starve(); }

    // The state check is cheap enough to run on EVERY tick; the dense read sweep behind it is not
    // (it is ~8 000 stitched reads through both panes), so callers run that one periodically.
    void compareBins()
    {
        for (int k = 0; k < ref->tierCount(); ++k)
            for (int i = 0; i < ref->tierBins (k); ++i)
            {
                if (ref->tierBinPower (k, i) != fast->tierBinPower (k, i)) fillBitExact = false;
                worstFillPow = std::max (worstFillPow, std::fabs ((double) ref->tierBinPower (k, i) - (double) fast->tierBinPower (k, i)));
            }
    }

    void compare()                                                        // state + the dense read sweep
    {
        compareBins();
        for (double f = 20.0; f < 24000.0; f *= 1.0116)                  // ~1/60 octave, through both seams
        {
            for (double tilt : { 0.0, 4.5, -6.0 })
            {
                const double a = ref->readDb (f, fs, tilt, 1000.0), b = fast->readDb (f, fs, tilt, 1000.0);
                const double c = ref->readPeakDb (f, fs, tilt, 1000.0), d = fast->readPeakDb (f, fs, tilt, 1000.0);
                worstFillDb = std::max (worstFillDb, std::fabs (a - b));
                // The peak budget applies where a reading can be DRAWN — the value here is already
                // tilted, and the deepest plot bottom in use is −120 dB. Below that the band integral
                // itself runs out of precision and neither pane is authoritative: a band is a
                // DIFFERENCE of two prefix sums whose magnitude is set by the loudest bin in the tier,
                // so a band 150 dB under a full-scale tone comes out around fifteen ulps of the
                // accumulation, where one ulp is worth ~0.27 dB. The group below measures that
                // directly and shows the SIBLING missing a cancellation-free oracle by more than the
                // budget, so the gate is the sibling's limit, not this pane's.
                // Above the line the two agree to ~1e-4 dB, four orders inside the budget.
                if (c > -120.0) worstPeakDb = std::max (worstPeakDb, std::fabs (c - d));
                else            deepestLifted = std::max (deepestLifted, d);
                if (d < b - 1e-9) peakAtLeastFill = false;                // the fast pane's own invariant
            }
        }
    }

    // buildColumns through several maps: the plan must reproduce the direct read exactly on the fill.
    void compareColumns()
    {
        for (auto [w, h, top, bot] : { std::tuple<float,float,double,double> { 900.0f, 300.0f,  6.0, -90.0 },
                                       { 360.0f, 120.0f,  0.0, -120.0 },
                                       { 256.0f,  96.0f, 12.0, -60.0 } })
            for (double tilt : { 0.0, 4.5, -6.0 })
            {
                analysis::PlotMap pm; pm.width = w; pm.height = h; pm.plotBottom = h; pm.specTop = top; pm.specBottom = bot;
                std::vector<float> ra, rp, fa, fp;
                ref ->buildColumns (pm, fs, tilt, 1000.0, [&] (int, float, float y, float yp) { ra.push_back (y); rp.push_back (yp); });
                fast->buildColumns (pm, fs, tilt, 1000.0, [&] (int, float, float y, float yp) { fa.push_back (y); fp.push_back (yp); });
                if (ra.size() != fa.size()) { worstColFill = 1e9; return; }
                for (std::size_t i = 0; i < ra.size(); ++i)
                {
                    worstColFill = std::max (worstColFill, (double) std::fabs (ra[i] - fa[i]));
                    worstColPeak = std::max (worstColPeak, (double) std::fabs (rp[i] - fp[i]));
                }
            }
    }
};

void fillSine (float* f, int n, double hz, double amp, double phase = 0.0)
{
    for (int i = 0; i < n; ++i) f[i] = (float) (amp * std::sin (2.0 * kPi * hz * (double) i / fs + phase));
}
}

int main()
{
    std::printf ("felitronics::analysis::MultiResSpectrumPaneFast — paired NULL against MultiResSpectrumPane\n");
    constexpr double kPeakTol = 1e-3;                                     // the sibling's dB quantisation, nothing more

    //==========================================================================================
    group ("rolling noise: three hops x three speeds, state compared every tick");
    {
        double worstFill = 0, worstPeak = 0, worstCol = 0, deepest = -1e9; bool bitExact = true, ordered = true;
        for (int hop : { 800, 1600, 3200 })
            for (float sm : { 0.05f, 0.25f, 0.5f })
            {
                Pair p; p.tune (sm, hop); Rng r;
                std::vector<float> ring ((std::size_t) N);
                for (auto& v : ring) v = (float) r.uni();
                for (int t = 0; t < 16; ++t)
                {
                    std::copy (ring.begin() + hop, ring.end(), ring.begin());
                    for (int i = N - hop; i < N; ++i) ring[(std::size_t) i] = (float) r.uni();
                    p.tick ([&] (float* f, int n) { std::copy (ring.begin(), ring.begin() + n, f); });
                    if (t % 5 == 4 || t == 15) p.compare(); else p.compareBins();
                }
                p.compareColumns();
                worstFill = std::max (worstFill, p.worstFillDb); worstPeak = std::max (worstPeak, p.worstPeakDb);
                worstCol  = std::max (worstCol, std::max (p.worstColFill, p.worstColPeak));
                deepest   = std::max (deepest, p.deepestLifted);
                bitExact = bitExact && p.fillBitExact; ordered = ordered && p.peakAtLeastFill;
            }
        ok (bitExact, "every bin of the fill state is bit-identical on every tick of all nine runs");
        ok (worstFill == 0.0, "the stitched FILL read is bit-identical too (worst " + std::to_string (worstFill) + " dB)");
        ok (worstPeak <= kPeakTol, "the peak read is within the sibling's own dB quantisation (worst " + std::to_string (worstPeak) + " dB)");
        ok (deepest < -120.0, "and where the sibling reads below −120 dB the fast pane does too — nothing invisible is lifted into view (worst " + std::to_string (deepest) + " dB)");
        ok (worstCol <= 0.02, "buildColumns agrees to well under a pixel (worst " + std::to_string (worstCol) + " px)");
        ok (ordered, "the fast pane's peak never reads below its fill");
    }

    //==========================================================================================
    // The test with teeth. A band's power is a DIFFERENCE of two large partial sums; if the prefix
    // arrays ever stop being double, a quiet band beside a loud one is erased by cancellation — and a
    // noise-only NULL would never notice, because noise has no loud partial sum to cancel against.
    group ("loud neighbour: a full-scale 50 Hz tone with a −100 dB probe at 8 kHz");
    {
        Pair p; p.tune (0.25f, 1600);
        for (int t = 0; t < 12; ++t)
            p.tick ([&] (float* f, int n)
            {
                for (int i = 0; i < n; ++i)
                    f[i] = (float) (std::sin (2.0 * kPi * 50.0 * (double) i / fs)
                                  + 1.0e-5 * std::sin (2.0 * kPi * 8000.0 * (double) i / fs));
            });
        p.compare(); p.compareColumns();
        const double probeRef = p.ref->readDb (8000.0, fs), probeFast = p.fast->readDb (8000.0, fs);
        ok (p.fillBitExact, "the fill state survives the 100 dB span bit-identically");
        ok (p.worstFillDb == 0.0, "and so does every stitched fill read");
        ok (p.worstPeakDb <= kPeakTol, "the peak read stays inside the budget beside a full-scale neighbour");
        ok (probeRef < -80.0 && probeRef > -120.0, "the probe is actually visible in the reference (" + std::to_string (probeRef) + " dB) — the test can fail");
        approx (probeFast, probeRef, 1e-9, "...and the fast pane reads the same quiet probe, not the floor");
    }

    //==========================================================================================
    group ("tones on every hard edge: bin centres, half bins, both seams, the blend ends, Nyquist");
    {
        Pair p; p.tune (1.0f, 1600);
        auto probe = [&] (double hz) { p.reset(); p.tick ([&] (float* f, int n) { fillSine (f, n, hz, 1.0); }); p.compare(); };
        const double s1 = p.ref->seamHz (1, fs), s2 = p.ref->seamHz (2, fs);
        for (double hz : { 20.0, 50.0, 58.59375, 60.0, 100.0, 440.0, 1000.0,
                           s1 * 0.999, s1, s1 * 1.001, s1 * 1.2599, s2 * 0.999, s2, s2 * 1.001, s2 * 1.2599,
                           5000.0, 11025.0, 19000.0, 23999.0 })
            probe (hz);
        probe (24000.0);                                                  // exactly Nyquist
        p.compareColumns();
        ok (p.fillBitExact, "bit-identical fill state on every edge tone");
        ok (p.worstFillDb == 0.0, "bit-identical fill reads (seams, blends, bin centres, Nyquist)");
        ok (p.worstPeakDb <= kPeakTol, "peak reads within budget (worst " + std::to_string (p.worstPeakDb) + " dB)");
        ok (p.deepestLifted < -120.0, "nothing below −120 dB in the sibling is lifted into view (worst " + std::to_string (p.deepestLifted) + " dB)");
        ok (p.worstColFill == 0.0, "and the column plan reproduces the direct read exactly");
    }

    //==========================================================================================
    group ("silence, starve, transients and poisoned input");
    {
        Pair p; p.tune (0.25f, 1600); Rng r;
        for (int t = 0; t < 6; ++t) p.tick ([&] (float* f, int n) { for (int i = 0; i < n; ++i) f[i] = (float) r.uni(); });
        for (int t = 0; t < 40; ++t) { p.starveBoth(); p.compare(); }
        ok (p.worstFillDb == 0.0 && p.worstPeakDb <= kPeakTol, "starve x40 tracks tick for tick");

        // A short decay first (both fills must track each other on the way down), then the exact-floor
        // check with the smoothing wide open, where one zero frame puts the state at 0 and the reading
        // must be the floor ITSELF, not merely near it. Waiting ~240 ticks for the 0.25 one-pole to get
        // there proves nothing extra and costs the sanitizer build minutes.
        p.reset();
        for (int t = 0; t < 40; ++t) { p.tick ([] (float* f, int n) { std::fill (f, f + n, 0.0f); }); p.compareBins(); }
        p.compare();
        p.reset(); p.tune (1.0f, 1600);
        p.tick ([] (float* f, int n) { std::fill (f, f + n, 0.0f); });
        p.compare();
        bool floorExact = true;
        for (double f = 25.0; f < 20000.0; f *= 1.07)
            if (p.fast->readDb (f, fs) != (double) Fast::kFloorDb || p.ref->readDb (f, fs) != (double) Ref::kFloorDb) floorExact = false;
        ok (floorExact, "digital silence reads exactly the floor in both, at every probe");

        // And a pane that has never seen a frame at all. The peak array rests on a float sentinel, so
        // "a hold on the floor counts as nothing" has to survive that float's rounding; if it did not,
        // every bin would contribute ~1e-20 and a fresh pane would read about −161 dB, not the floor.
        p.reset();
        bool freshFloor = true;
        for (double f = 25.0; f < 20000.0; f *= 1.07)
            if (p.fast->readPeakDb (f, fs) != (double) Fast::kFloorDb || p.fast->readDb (f, fs) != (double) Fast::kFloorDb) freshFloor = false;
        ok (freshFloor, "a pane that has never ingested reads exactly the floor on both traces");
        p.tune (0.25f, 1600);

        p.reset();                                                        // a burst confined to the last 512 samples
        p.tick ([] (float* f, int n) { std::fill (f, f + n, 0.0f); for (int i = n - 512; i < n; ++i) f[i] = (i % 2) ? -0.9f : 0.9f; });
        p.compare();
        ok (p.worstFillDb == 0.0, "a 512-sample burst: the short tier sees it identically");

        p.reset();                                                        // NaN / Inf samples must not diverge
        p.tick ([] (float* f, int n)
        {
            for (int i = 0; i < n; ++i) f[i] = 0.3f;
            f[10] = std::numeric_limits<float>::quiet_NaN();
            f[N - 3] = std::numeric_limits<float>::infinity();
        });
        p.compare();
        ok (p.worstFillDb == 0.0 && std::isfinite (p.fast->readDb (1000.0, fs)), "NaN/Inf samples: dropped identically, reads stay finite");
    }

    //==========================================================================================
    // A frame shorter than the longest tier: that tier HOLDS and stays invalid, so tierAt has to skip
    // it and the column plan's validity mask has to notice when it later becomes valid. This is the one
    // path where the plan can go stale without any of its numeric inputs moving.
    group ("a short frame: the long tier holds, and the plan follows it back to valid");
    {
        Pair p; p.tune (0.25f, 1600); Rng r;
        double worstFill = 0.0, worstPeak = 0.0; bool skipped = false;
        auto sweep = [&] ()
        {
            analysis::PlotMap pm; pm.width = 900.0f; pm.height = 300.0f; pm.plotBottom = 300.0f;
            std::vector<float> af, ap, bf, bp;
            p.ref ->buildColumns (pm, fs, 4.5, 1000.0, [&] (int, float, float y, float yp) { af.push_back (y); ap.push_back (yp); });
            p.fast->buildColumns (pm, fs, 4.5, 1000.0, [&] (int, float, float y, float yp) { bf.push_back (y); bp.push_back (yp); });
            for (std::size_t i = 0; i < af.size() && i < bf.size(); ++i)
            {
                worstFill = std::max (worstFill, (double) std::fabs (af[i] - bf[i]));
                worstPeak = std::max (worstPeak, (double) std::fabs (ap[i] - bp[i]));
            }
        };
        for (int t = 0; t < 5; ++t)                                       // 4096-sample frames: tier 0 never sees one
        {
            p.tick ([&] (float* f, int n) { for (int i = 0; i < n; ++i) f[i] = (float) r.uni(); }, 12);
            skipped = skipped || (p.ref->tierAt (30.0, fs) != 0);         // the read fell back off tier 0
            p.compare(); sweep();
        }
        ok (skipped, "with a short frame the reads really do skip the tier that never ran");
        for (int t = 0; t < 5; ++t)                                       // now the full frame arrives
        {
            p.tick ([&] (float* f, int n) { for (int i = 0; i < n; ++i) f[i] = (float) r.uni(); }, 14);
            p.compare(); sweep();
        }
        ok (p.ref->tierAt (30.0, fs) == 0, "and once the long tier has a frame the lows come from it again");
        ok (worstFill == 0.0, "the plan tracked the validity change with a bit-identical fill (worst " + std::to_string (worstFill) + " px)");
        ok (worstPeak < 0.01, "and the peak column agrees to far under a pixel (worst " + std::to_string (worstPeak) + " px)");
        ok (p.worstFillDb == 0.0 && p.worstPeakDb <= kPeakTol, "reads stayed in budget throughout");
    }

    //==========================================================================================
    // Where a power-domain hold could drift that its dB ancestor could not: the decay is now a repeated
    // MULTIPLY, so the question is whether 10^(-fall/10) applied thousands of times walks away from
    // subtracting `fall` thousands of times. starve() is that multiply with no FFT attached, so it tests
    // the thing directly and cheaply — 4000 of them cost less than a hundred ingests. The ingest loop
    // that follows covers the smoothing path over a long run.
    group ("endurance: the peak decay applied 4000 times, and 120 ticks of noise");
    {
        Pair p; p.tune (0.25f, 1600); Rng r;
        for (int t = 0; t < 6; ++t) p.tick ([&] (float* f, int n) { for (int i = 0; i < n; ++i) f[i] = (float) r.uni(); });
        for (int t = 0; t < 4000; ++t) { p.starveBoth(); if (t % 200 == 0 || t > 3990) p.compare(); }
        ok (p.worstPeakDb <= kPeakTol, "4000 decay steps: the multiply has not walked away from the subtraction (worst " + std::to_string (p.worstPeakDb) + " dB)");

        Pair q; q.tune (0.25f, 1600);
        for (int t = 0; t < 120; ++t)
        {
            q.tick ([&] (float* f, int n) { for (int i = 0; i < n; ++i) f[i] = (float) (0.3 * r.uni()); });
            if (t % 30 == 0 || t > 116) q.compare(); else q.compareBins();
        }
        ok (q.fillBitExact && q.worstFillDb == 0.0, "120 ticks of noise: the fill has not drifted a single bit");
        ok (q.worstPeakDb <= kPeakTol, "nor has the peak left its budget (worst " + std::to_string (q.worstPeakDb) + " dB)");
    }

    //==========================================================================================
    group ("the column plan rebuilds when, and only when, its inputs move");
    {
        Pair p; p.tune (0.25f, 1600); Rng r;
        for (int t = 0; t < 4; ++t) p.tick ([&] (float* f, int n) { for (int i = 0; i < n; ++i) f[i] = (float) r.uni(); });
        double worst = 0.0;
        auto sweepOnce = [&] (float w, double rate, double tilt)
        {
            analysis::PlotMap pm; pm.width = w; pm.height = 300.0f; pm.plotBottom = 300.0f;
            std::vector<float> a, b;
            p.ref ->buildColumns (pm, rate, tilt, 1000.0, [&] (int, float, float y, float) { a.push_back (y); });
            p.fast->buildColumns (pm, rate, tilt, 1000.0, [&] (int, float, float y, float) { b.push_back (y); });
            if (a.size() != b.size()) { worst = 1e9; return; }
            for (std::size_t i = 0; i < a.size(); ++i) worst = std::max (worst, (double) std::fabs (a[i] - b[i]));
        };
        for (float w : { 256.0f, 257.0f, 300.0f, 640.0f, 899.0f, 900.0f, 1600.0f, 0.0f })  // incl. the clamps
            sweepOnce (w, fs, 4.5);
        for (double rate : { 44100.0, 48000.0, 88200.0, 96000.0 }) sweepOnce (900.0f, rate, 1.5);
        for (double tilt : { -6.0, 0.0, 1.5, 4.5 })                        sweepOnce (900.0f, fs, tilt);
        p.fast->bandOctaves = 1.0 / 6.0;  p.ref->bandOctaves = 1.0 / 6.0;  sweepOnce (900.0f, fs, 4.5);
        p.fast->blendOctaves = 0.0;       p.ref->blendOctaves = 0.0;       sweepOnce (900.0f, fs, 4.5);
        p.fast->binsPerBand = 4.0;        p.ref->binsPerBand = 4.0;        sweepOnce (900.0f, fs, 4.5);
        p.fast->bandOctaves = 1.0 / 24.0; p.ref->bandOctaves = 1.0 / 24.0;
        p.fast->blendOctaves = 1.0 / 3.0; p.ref->blendOctaves = 1.0 / 3.0;
        p.fast->binsPerBand = 2.0;        p.ref->binsPerBand = 2.0;        sweepOnce (900.0f, fs, 4.5);
        const int two[2] = { 13, 10 };                                     // and a whole new ladder
        p.fast->setTiers (two, 2); p.ref->setTiers (two, 2);
        for (int t = 0; t < 3; ++t) p.tick ([&] (float* f, int n) { for (int i = 0; i < n; ++i) f[i] = (float) r.uni(); }, 13);
        sweepOnce (900.0f, fs, 4.5);
        ok (worst == 0.0, "every width, sample rate, tilt, tuning edit and setTiers still reads bit-identically (worst " + std::to_string (worst) + " px)");
    }

    //==========================================================================================
    // Why the peak budget stops at −120 dB, measured rather than asserted. Beside a full-scale 11 kHz
    // tone the 23.8 kHz band sits ~150 dB down while the prefix sums it is subtracted from carry the
    // tone's own power, so the result is a handful of double ulps. An independent whole-bin oracle,
    // built from each pane's OWN published per-bin peaks and summed without any cancellation, is the
    // arbiter: if the sibling misses it by more than the fast pane's budget, the sibling is not a
    // reference at that depth and no tolerance can pretend otherwise.
    group ("below the gate it is the band integral that runs out, not this pane");
    {
        Pair p; p.tune (1.0f, 1600);
        p.tick ([] (float* f, int n) { fillSine (f, n, 11025.0, 1.0); });
        const int    k = 2;
        const double f = 23792.82;
        const double binHz = fs / (double) (1 << p.ref->tierOrder (k));
        const double h = std::exp2 (0.5 / 24.0), uMax = (double) p.ref->tierBins (k) - 0.5;
        const double uLo = std::clamp ((f / h) / binHz + 0.5, 0.5, uMax);
        const double uHi = std::clamp ((f * h) / binHz + 0.5, 0.5, uMax);
        auto oracle = [&] (auto&& pane)                                   // whole-bin, no differencing
        {
            double acc = 0.0;
            for (int i = 0; i < pane->tierBins (k); ++i)
            {
                const double lo = std::max (uLo, (double) i), hi = std::min (uHi, (double) i + 1.0);
                const double db = (double) pane->tierBinPeak (k, i);
                if (hi > lo && db > -200.0) acc += (hi - lo) * std::pow (10.0, db / 10.0);
            }
            return 10.0 * std::log10 (acc);
        };
        const double oR = oracle (p.ref), oF = oracle (p.fast);
        const double bR = p.ref->tierBandPeakDb (k, f, fs), bF = p.fast->tierBandPeakDb (k, f, fs);
        approx (oF, oR, 1e-4, "both panes publish the same per-bin peaks here, so the oracle is one number");
        ok (bR < -140.0, "the band under test really is ~150 dB below the tone (" + std::to_string (bR) + " dB)");
        ok (std::fabs (bR - oR) > kPeakTol,
            "the SIBLING misses the cancellation-free oracle by more than the budget (" + std::to_string (std::fabs (bR - oR)) + " dB) — it is not a reference this deep");
        ok (std::fabs (bF - oF) < 1.0 && std::fabs (bR - oR) < 1.0,
            "both stay within a dB of the oracle, so this is lost precision, not a wrong answer");
    }

    //==========================================================================================
    group ("the two deliberate divergences, pinned");
    {
        auto f = std::make_unique<Fast>();
        f->smoothCoeff = 1.0f; f->peakFallDb = -3.0f;                      // unsupported, but must stay finite
        for (int t = 0; t < 400; ++t) { fillSine (f->frameInput(), N, 1000.0, 0.5); f->ingest (14); }
        // Per BIN the hold is clamped at kMaxPower (+120 dB); a BAND sums its bins, so its reading sits
        // above that by 10·log10(bins). The point of the clamp is that nothing reaches an infinity.
        const double held = f->readPeakDb (1000.0, fs);
        ok (std::isfinite (held) && held < 200.0,
            "a negative peakFallDb is clamped instead of reaching an infinity (" + std::to_string (held) + " dB)");
        f->peakFallDb = std::numeric_limits<float>::quiet_NaN();
        fillSine (f->frameInput(), N, 1000.0, 0.5); f->ingest (14);
        ok (std::isfinite (f->readPeakDb (1000.0, fs)), "a non-finite peakFallDb does not poison the hold");
    }

    //==========================================================================================
    group ("no allocation on the hot path, plan rebuild included");
    {
        auto f = std::make_unique<Fast>();
        analysis::PlotMap pm; pm.width = 900.0f; pm.height = 300.0f; pm.plotBottom = 300.0f;
        Rng r;
        for (int t = 0; t < 3; ++t) { for (int i = 0; i < N; ++i) f->frameInput()[i] = (float) r.uni(); f->ingest (14); }
        f->buildColumns (pm, fs, 4.5, 1000.0, [] (int, float, float, float) {});
        const long before = g_allocs.load();
        for (int i = 0; i < N; ++i) f->frameInput()[i] = (float) r.uni();
        f->ingest (14);
        f->starve();
        f->buildColumns (pm, fs, 4.5, 1000.0, [] (int, float, float, float) {});
        pm.width = 640.0f;                                                 // forces a full plan rebuild
        f->buildColumns (pm, fs, 4.5, 1000.0, [] (int, float, float, float) {});
        (void) f->readDb (1000.0, fs); (void) f->tierAt (1000.0, fs);
        test::okNoAlloc (g_allocs.load() == before, "ingest / starve / buildColumns / plan rebuild / reads allocate nothing");
    }

    return test::report();
}

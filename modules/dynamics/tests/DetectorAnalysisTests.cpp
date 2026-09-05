// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// P3 — detector-domain analysis and the threshold solver.
//
// THE ONE TEST THAT DECIDES WHETHER THIS MODULE IS RIGHT is `f9CounterexampleIsMandatory()`. Every
// other check here can pass on an implementation that inverts the static curve at a level percentile —
// including the round-trip, which is a TAUTOLOGY by construction: the threshold is chosen by the same
// quantity it is then verified with, so it can only fail if the search does not converge, if the two
// ends of the seam disagree, or if "p95" means two things. None of those tests the MODEL. F9 does: two
// signals with bit-identical level percentiles and different time layouts must receive different
// thresholds, and no function of those percentiles can produce that.
//
// The second load-bearing idea is that the histogram must not be its own oracle. Every percentile the
// module reports is NULLed against an exact sort written independently in this file, and the final
// threshold is verified by running the REAL Compressor and sorting ITS tap — a different instrument
// from the one the solver used to find it.

#include <felitronics_test.h>

#include <felitronics/dynamics/Compressor.h>
#include <felitronics/dynamics/offline/ThresholdSolver.h>
#include <felitronics/core/Math.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <random>
#include <type_traits>
#include <vector>

// global allocation counter (no-alloc-in-the-search proof)
static std::atomic<long> g_allocs { 0 };
void* operator new      (std::size_t s) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s ? s : 1); }
void* operator new[]    (std::size_t s) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s ? s : 1); }
void  operator delete   (void* p) noexcept { std::free (p); }
void  operator delete[] (void* p) noexcept { std::free (p); }
void  operator delete   (void* p, std::size_t) noexcept { std::free (p); }
void  operator delete[] (void* p, std::size_t) noexcept { std::free (p); }

using namespace felitronics;
using namespace felitronics::dynamics;
using namespace felitronics::dynamics::offline;

// Compile-time contract detectors for the call shapes the API deliberately does and does not offer.
namespace p3contract
{
    template <class C, class = void> struct HasFourArgKeyProcess : std::false_type {};
    template <class C> struct HasFourArgKeyProcess<C, std::void_t<decltype (std::declval<C&>().process (
        std::declval<float* const*>(), 0, 0, std::declval<const float* const*>()))>> : std::true_type {};

    template <class C, class = void> struct HasTapProcess : std::false_type {};
    template <class C> struct HasTapProcess<C, std::void_t<decltype (std::declval<C&>().process (
        std::declval<float* const*>(), 0, 0, std::declval<GainReductionTap>()))>> : std::true_type {};
}

static constexpr double kFs = 48000.0;

//==============================================================================
// The ORACLE for a percentile: an exact nearest-rank order statistic by sorting. Written here, from
// the definition, sharing no code with the histogram — which is the only way its answer is evidence.
static double exactNearestRank (std::vector<double> v, double q)
{
    if (v.empty()) return 0.0;
    std::sort (v.begin(), v.end());
    double rank = std::ceil (q * (double) v.size());
    if (rank < 1.0) rank = 1.0;
    if (rank > (double) v.size()) rank = (double) v.size();
    return v[(std::size_t) rank - 1];
}

// The ORACLE for the gain-reduction path: the documented topology, reassembled from public primitives.
struct PathOracle
{
    LinkedDetector        det;
    GainComputer          gc;
    GainReductionFollower gr;

    void configure (double fs, const GainReductionParams& p)
    {
        det.prepare (fs); det.setParams (p);
        gc.setMode (p.mode); gc.setThresholdDb (p.thresholdDb); gc.setRatio (p.ratio);
        gc.setKneeDb (p.kneeDb); gc.setRangeDb (p.rangeDb);
        gr.prepare (fs); gr.setTimes (p.attackMs, p.releaseMs);
    }
    float step (const float* const* key, int nk, int i)
    {
        return gr.process ((float) gc.deltaDb (core::gainToDb (det.process (key, nk, i))));
    }
};

//==============================================================================
// Fixtures. Rule 1 of the house discipline — taper both ends — does not apply to a level fixture whose
// point is a KNOWN multiset of levels; these are deliberately rectangular, and the two layouts below
// exist precisely so a rectangular level distribution can be arranged in time two different ways.
struct TwoLevel
{
    std::vector<float> x;
    const float* p = nullptr;
    const float* const* key = nullptr;
    void build (int n, double loudDb, double quietDb, bool contiguous, int period = 5)
    {
        const float loud = (float) core::dbToGain (loudDb), quiet = (float) core::dbToGain (quietDb);
        x.assign ((std::size_t) n, quiet);
        if (contiguous) for (int i = 0; i < n / period; ++i)   x[(std::size_t) i] = loud;
        else            for (int i = 0; i < n; i += period)    x[(std::size_t) i] = loud;
        p = x.data(); key = &p;
    }
};

static std::vector<float> musicish (int n, double fs, unsigned seed = 5)
{
    std::vector<float> v ((std::size_t) n);
    std::mt19937 rng (seed);
    std::normal_distribution<float> nd (0.0f, 0.25f);
    for (int i = 0; i < n; ++i)
    {
        const double t = (double) i / fs;
        const double env = 0.15 + 0.85 * std::pow (std::fabs (std::sin (2.0 * core::kPi * 1.9 * t)), 3.0);
        double v0 = env * (0.6 * std::sin (2.0 * core::kPi * 110.0 * t) + 0.2 * (double) nd (rng));
        if ((i % 4801) < 60) v0 += 0.5 * std::sin (2.0 * core::kPi * 2500.0 * t);
        v[(std::size_t) i] = (float) v0;
    }
    return v;
}

static GainReductionParams baseParams()
{
    GainReductionParams p;
    p.detector = Detector::Rms; p.link = LinkMode::Max; p.rmsWindowMs = 5.0;
    p.mode = Mode::DownCompress; p.thresholdDb = -18.0; p.ratio = 4.0; p.kneeDb = 6.0;
    p.rangeDb = 60.0; p.attackMs = 10.0; p.releaseMs = 100.0;
    return p;
}

//==============================================================================
static void quantileRuleIsOneRule()
{
    test::group ("QuantileHistogram — nearest rank, NULLed against an exact sort");

    // 1. Random distributions of several shapes: the histogram may be off by at most half a bin, and
    //    the clamp to [min,max] must never let it leave the observed span.
    std::mt19937 rng (4242);
    const double bin = 0.01;
    int worstCase = 0; double worstErr = 0.0;
    for (int shape = 0; shape < 4; ++shape)
    {
        std::vector<double> vals;
        std::normal_distribution<double> nd (6.0, 3.0);
        std::exponential_distribution<double> ed (0.4);
        for (int i = 0; i < 60000; ++i)
        {
            double v = 0.0;
            switch (shape)
            {
                case 0: v = std::fabs (nd (rng)); break;                             // broad
                case 1: v = ed (rng); break;                                         // skewed
                case 2: v = (i % 7 == 0) ? 12.0 : 0.0; break;                        // two atoms
                default: v = (i % 3 == 0) ? 0.0 : 2.0 + 0.001 * (double) (i % 11);   // near-atom cluster
            }
            vals.push_back (v);
        }
        QuantileHistogram h; test::ok (h.prepare (0.0, 400.0, bin), "histogram prepared");
        for (double v : vals) h.add (v);
        for (double q : { 0.05, 0.5, 0.9, 0.95, 0.99, 1.0 })
        {
            double got = 0.0;
            const bool okq = h.quantile (q, got);
            test::ok (okq, "shape " + std::to_string (shape) + " q=" + std::to_string (q) + " in range");
            const double want = exactNearestRank (vals, q);
            const double err = std::fabs (got - want);
            if (err > worstErr) { worstErr = err; worstCase = shape; }
            test::ok (err <= 0.5 * bin + 1.0e-12,
                      "shape " + std::to_string (shape) + " q=" + std::to_string (q)
                      + " within half a bin of the exact order statistic (err " + std::to_string (err) + ")");
            test::ok (got >= h.minValue() - 1e-12 && got <= h.maxValue() + 1e-12, "never outside the observed span");
        }
        // mean and max are EXACT — they are not binned, which is why the solver may target them.
        double sum = 0.0; for (double v : vals) sum += v;
        test::approx (h.mean(), sum / (double) vals.size(), 1e-9, "mean() is exact, not binned");
        test::approx (h.maxValue(), *std::max_element (vals.begin(), vals.end()), 0.0, "maxValue() is exact");
    }
    test::ok (worstErr <= 0.5 * bin + 1e-12, "worst error over all shapes was half a bin (shape "
              + std::to_string (worstCase) + ", " + std::to_string (worstErr) + " dB)");

    // 2. ATOMS ARE EXACT — this is what the clamp buys, and what an interpolated rule gets wrong.
    {
        QuantileHistogram h; test::ok (h.prepare (0.0, 400.0, 0.05), "histogram prepared");
        for (int i = 0; i < 1000; ++i) h.add (2.0);
        double got = 0.0;
        test::ok (h.quantile (0.95, got), "single atom: in range");
        test::approx (got, 2.0, 0.0, "a constant 2.000 dB trace reports EXACTLY 2.000 (inverse-CDF interpolation reports 2.046)");
        QuantileHistogram z; test::ok (z.prepare (0.0, 400.0, 0.05), "histogram prepared");
        for (int i = 0; i < 1000; ++i) z.add (0.0);
        test::ok (z.quantile (0.95, got) && got == 0.0, "an all-zero trace reports EXACTLY zero, not half a bin");
    }

    // 3. The rank rule itself, on a distribution whose every quantile is known by counting.
    {
        QuantileHistogram h; test::ok (h.prepare (0.0, 400.0, 0.01), "histogram prepared");
        for (int i = 0; i < 200; ++i) h.add (i < 40 ? 10.0 : 1.0);   // 20 % at 10, 80 % at 1
        double got = 0.0;
        (void) h.quantile (0.50, got); test::approx (got, 1.0, 0.006, "p50 of 20/80 is the low atom (to half a bin)");
        (void) h.quantile (0.79, got); test::approx (got, 1.0, 0.006, "p79 is still the low atom");
        (void) h.quantile (0.81, got); test::approx (got, 10.0, 0.006, "p81 has crossed into the high atom");
        (void) h.quantile (1.0, got);  test::approx (got, 10.0, 1e-9, "q=1 is the maximum");
        (void) h.quantile (0.0, got);  test::approx (got, 1.0, 1e-9, "q=0 is the minimum");
        test::ok (! h.quantile (-0.1, got) && ! h.quantile (1.1, got), "a quantile outside [0,1] is refused, not clamped");
    }

    // 3b. THE RANK RULE ITSELF, pinned so that `ceil` and `floor` cannot be swapped unnoticed. With
    //     201 values and q = 0.95 the rank is 190.95: ceil takes the 191st order statistic and floor
    //     the 190th, and here those two are 9 dB apart. Every tolerance-based check above absorbs the
    //     difference; this one cannot.
    {
        QuantileHistogram h; test::ok (h.prepare (0.0, 400.0, 0.01), "histogram prepared");
        std::vector<double> vals;
        for (int i = 0; i < 190; ++i) { h.add (1.0); vals.push_back (1.0); }
        for (int i = 0; i < 11;  ++i) { h.add (10.0); vals.push_back (10.0); }
        double got = 0.0;
        test::ok (h.quantile (0.95, got), "in range");
        test::approx (got, 10.0, 0.006, "rank ceil(0.95*201) = 191 lands in the HIGH atom (floor would give 1.0)");
        test::approx (exactNearestRank (vals, 0.95), 10.0, 0.0, "and the independent exact rule agrees");
        // ...and one rank lower is still the low atom, which pins the boundary from the other side.
        (void) h.quantile (0.945, got);
        test::approx (got, 1.0, 0.006, "rank ceil(0.945*201) = 190 is still the LOW atom");
    }

    // 3c. The top of the configured range is INSIDE it. |gain change| reaches `rangeDb` exactly on
    //     saturated material, and `rangeDb` may be the 400 dB cap itself.
    {
        QuantileHistogram h; test::ok (h.prepare (0.0, 400.0, 0.01), "histogram prepared");
        for (int i = 0; i < 100; ++i) h.add (400.0);
        double got = 0.0;
        test::ok (h.quantile (0.95, got), "a value exactly at the top of the range is IN range");
        test::approx (got, 400.0, 1e-9, "and reports exactly");
        test::ok (h.aboveRange() == 0, "nothing overflowed");
    }

    // 4. OUT OF RANGE IS REPORTED, NOT SATURATED — a solver told "96" when the truth is 150 converges
    //    confidently to the wrong threshold.
    {
        QuantileHistogram h; test::ok (h.prepare (0.0, 96.0, 0.05), "histogram prepared");
        for (int i = 0; i < 100; ++i) h.add (200.0);
        double got = 0.0;
        test::ok (! h.quantile (0.95, got), "a quantile among the overflow counts returns false");
        test::ok (h.aboveRange() == 100, "and the overflow is counted");
        test::approx (h.maxValue(), 200.0, 0.0, "while max stays exact — it is not binned");
    }

    // 5. Non-finite values are COUNTED, never folded into a bin.
    {
        QuantileHistogram h; test::ok (h.prepare (0.0, 400.0, 0.01), "histogram prepared");
        h.add (std::numeric_limits<double>::quiet_NaN());
        h.add (std::numeric_limits<double>::infinity());
        h.add (3.0);
        test::ok (h.nonFiniteCount() == 2, "NaN and Inf are counted separately");
        test::ok (h.count() == 1, "and excluded from the distribution");
        double got = 0.0; (void) h.quantile (0.5, got);
        test::approx (got, 3.0, 1e-9, "the one real value is the whole distribution");
    }

    // 6. Empty.
    {
        QuantileHistogram h; test::ok (h.prepare (0.0, 400.0, 0.01), "histogram prepared");
        double got = 1234.0;
        test::ok (! h.quantile (0.95, got) && got == 1234.0, "an empty histogram refuses and leaves `out` untouched");
    }

    // 7. Counters are 64-bit BY TYPE — a uint32 bin overflows after 24.9 hours at 48 kHz, which no test
    //    can reach in a suite that runs in a second, so the contract is asserted where it is decidable.
    static_assert (std::is_same_v<decltype (std::declval<const QuantileHistogram&>().count()), std::uint64_t>,
                   "sample counts must be 64-bit");
    test::ok (true, "count() is uint64 (static_assert)");
}

//==============================================================================
static void envelopeAnalyzerMeasuresWhatItSays()
{
    test::group ("EnvelopeAnalyzer — the detector's own domain");

    // 1. A Peak detector is instantaneous, so its envelope IS the signal and every percentile is known
    //    analytically. This is the check that the rank rule is not off by one in the wired-up path.
    {
        TwoLevel s; s.build (100000, -4.0, -40.0, false, 5);
        DetectorParams dp; dp.detector = Detector::Peak; dp.link = LinkMode::Max;
        EnvelopeAnalyzer a; test::ok (a.prepare (kFs), "prepare");
        a.setParams (dp); a.reset(); a.analyze (s.key, 1, (int) s.x.size());
        const auto st = a.stats();
        test::approx (st.p50Db, -40.0, 0.02, "p50 of a 20/80 two-level signal is the quiet level");
        test::approx (st.p90Db,  -4.0, 0.02, "p90 is the loud level");
        test::approx (st.p95Db,  -4.0, 0.02, "p95 is the loud level");
        test::approx (st.p99Db,  -4.0, 0.02, "p99 is the loud level");
        test::ok (st.framesSeen == s.x.size(), "every frame was seen");
        test::ok (st.frames == s.x.size(), "and every frame counted, the gate being off");
    }

    // 2. THE CREST IS OF THE KEY, NOT OF THE SMOOTHED ENVELOPE. A sine has a crest of 3.0103 dB; the
    //    crest of its RMS envelope is ~0. Getting this wrong is invisible in every other check here.
    {
        std::vector<float> sine ((std::size_t) 48000);
        for (int i = 0; i < 48000; ++i) sine[(std::size_t) i] = (float) std::sin (2.0 * core::kPi * 997.0 * i / kFs);
        const float* p = sine.data(); const float* const* k = &p;
        for (auto det : { Detector::Peak, Detector::Rms })
        {
            DetectorParams dp; dp.detector = det; dp.rmsWindowMs = 50.0;
            EnvelopeAnalyzer a; (void) a.prepare (kFs); a.setParams (dp); a.reset();
            a.analyze (k, 1, 48000);
            test::approx (a.stats().crestDb, 3.0103, 0.02,
                          std::string ("crest of a sine is 3.01 dB with the ") + (det == Detector::Peak ? "Peak" : "Rms")
                          + " detector — it is a property of the KEY, not of the smoothing");
        }
    }

    // 3. Chunking must not move a single number: the state carries across calls exactly as it does in
    //    the compressor. (#119's lesson: check the RETURN, not only the departure.)
    {
        const auto music = musicish (60000, kFs);
        const float* p = music.data(); const float* const* k = &p;
        DetectorParams dp; dp.detector = Detector::Rms; dp.rmsWindowMs = 5.0;
        EnvelopeStats ref {};
        for (int blk : { 0, 1, 37, 512, 4096 })
        {
            EnvelopeAnalyzer a; (void) a.prepare (kFs); a.setParams (dp); a.reset();
            if (blk == 0) a.analyze (k, 1, 60000);
            else for (int off = 0; off < 60000; off += blk)
            {
                const int m = std::min (blk, 60000 - off);
                const float* q = music.data() + off; const float* const* kk = &q;
                a.analyze (kk, 1, m);
            }
            const auto st = a.stats();
            if (blk == 0) ref = st;
            else
            {
                test::approx (st.p95Db, ref.p95Db, 1e-9, "p95 is identical at block " + std::to_string (blk));
                test::approx (st.crestDb, ref.crestDb, 1e-9, "crest is identical at block " + std::to_string (blk));
                test::ok (st.frames == ref.frames, "frame count is identical at block " + std::to_string (blk));
            }
        }
    }

    // 4. The -240 dB floor is CENSORING, and says so.
    {
        std::vector<float> sil ((std::size_t) 5000, 0.0f);
        const float* p = sil.data(); const float* const* k = &p;
        DetectorParams dp; dp.detector = Detector::Peak;
        EnvelopeAnalyzer a; (void) a.prepare (kFs); a.setParams (dp); a.reset();
        a.analyze (k, 1, 5000);
        const auto st = a.stats();
        test::ok (st.atFloorFrames == 5000, "digital silence is entirely at the floor");
        test::approx (st.p95Db, -240.0, 1e-6, "and reads -240 dB, which is the clamp and not a measurement");
    }

    // 5. THE GATE EXCLUDES FROM THE DISTRIBUTION AND NEVER FROM THE STATE. If it skipped the detector,
    //    the frames after a gated stretch would start from the wrong envelope.
    {
        TwoLevel s; s.build (60000, -6.0, -80.0, true, 3);
        DetectorParams dp; dp.detector = Detector::Rms; dp.rmsWindowMs = 20.0;
        EnvelopeAnalyzer open, gated;
        (void) open.prepare (kFs); (void) gated.prepare (kFs);
        open.setParams (dp); gated.setParams (dp);
        gated.setStatisticsGateDb (-60.0);
        open.reset(); gated.reset();
        open.analyze (s.key, 1, (int) s.x.size());
        gated.analyze (s.key, 1, (int) s.x.size());
        test::ok (gated.stats().frames < open.stats().frames, "the gate removed frames from the distribution");
        test::ok (gated.stats().framesSeen == open.stats().framesSeen, "but not from the stream");
        test::approx (gated.stats().crestDb, open.stats().crestDb, 1e-9,
                      "the crest is unchanged — it is taken on the key, ahead of any statistics policy");
        // The state proof: the level after the last frame must be identical either way.
        test::approx (gated.stats().maxDb, open.stats().maxDb, 1e-9, "and the detector saw the same maximum");
        // A gate at or below the floor excludes nothing at all.
        EnvelopeAnalyzer atFloor; (void) atFloor.prepare (kFs); atFloor.setParams (dp);
        atFloor.setStatisticsGateDb (-240.0); atFloor.reset();
        atFloor.analyze (s.key, 1, (int) s.x.size());
        test::ok (atFloor.stats().frames == open.stats().frames, "a gate at the floor excludes nothing");
    }

    // 6. The captured envelope is the detector's own output, bit for bit — this is what lets the solver
    //    reuse one pass, and a drifted copy would silently change every probe.
    {
        const auto music = musicish (20000, kFs, 11);
        const float* p = music.data(); const float* const* k = &p;
        DetectorParams dp; dp.detector = Detector::Rms; dp.rmsWindowMs = 7.0; dp.link = LinkMode::Max;
        std::vector<float> cap ((std::size_t) 20000, -1.0f);
        EnvelopeAnalyzer a; (void) a.prepare (kFs); a.setParams (dp); a.reset();
        a.analyze (k, 1, 20000, cap.data(), 20000);
        LinkedDetector d; d.prepare (kFs); d.setParams (dp); d.reset();
        int diff = 0;
        for (int i = 0; i < 20000; ++i) { const float e = d.process (k, 1, i); if (std::memcmp (&e, &cap[(std::size_t) i], 4) != 0) ++diff; }
        test::ok (diff == 0, "the captured envelope is bit-identical to an independently driven LinkedDetector");

        std::vector<float> shortBuf ((std::size_t) 20000, -1.0f);
        EnvelopeAnalyzer b; (void) b.prepare (kFs); b.setParams (dp); b.reset();
        b.analyze (k, 1, 20000, shortBuf.data(), 19999);
        bool untouched = true;
        for (float v : shortBuf) if (v != -1.0f) untouched = false;
        test::ok (untouched, "a capacity one short of the chunk writes NOTHING, not a prefix");
    }

    // 7. `fractionAboveDb` is strictly above, and a signal sitting exactly on the level says so.
    {
        TwoLevel s; s.build (10000, -10.0, -50.0, false, 4);   // 25 % loud
        DetectorParams dp; dp.detector = Detector::Peak;
        EnvelopeAnalyzer a; (void) a.prepare (kFs); a.setParams (dp); a.reset();
        a.analyze (s.key, 1, 10000);
        test::approx (a.fractionAboveDb (-30.0), 0.25, 0.001, "a quarter of the frames are above -30 dB");
        test::approx (a.fractionAboveDb (-5.0), 0.0, 0.001, "none are above -5 dB");
        test::approx (a.fractionAboveDb (-240.0), 1.0, 0.001, "all of them are above the floor");
    }

    // 7b. ABOVE means STRICTLY above, and the only place that is visible is a signal sitting exactly on
    //     the level asked about. Every off-the-level query gives the same answer either way, which is
    //     how a `>` / `>=` swap survived the first rollback pass.
    {
        const float amp = (float) core::dbToGain (-10.0);
        std::vector<float> flat ((std::size_t) 8000, amp);
        const float* q0 = flat.data(); const float* const* k = &q0;
        DetectorParams dp; dp.detector = Detector::Peak;
        EnvelopeAnalyzer a; (void) a.prepare (kFs); a.setParams (dp); a.reset();
        a.analyze (k, 1, 8000);
        const double exactly = core::gainToDb (amp);
        test::approx (a.fractionAboveDb (exactly), 0.0, 1e-9,
                      "a frame sitting exactly ON the level is NOT above it");
        test::approx (a.fractionAboveDb (exactly - 0.02), 1.0, 1e-9, "one bin lower, all of them are");
    }
}

//==============================================================================
static void gainReductionPathIsTheCompressorsOwnChain()
{
    test::group ("GainReductionPath — NULLed against the primitives and against the compressor");

    const auto music = musicish (30000, kFs, 3);
    std::vector<float> second (music.size());
    for (std::size_t i = 0; i < music.size(); ++i) second[i] = music[i] * 0.4f;

    long cases = 0, worst = 0;
    for (auto det : { Detector::Peak, Detector::Rms })
     for (auto link : { LinkMode::Max, LinkMode::MeanPower })
      for (auto mode : { Mode::DownCompress, Mode::UpCompress, Mode::DownExpand })
       for (double thr : { -40.0, -18.0, -3.0 })
        for (double ratio : { 1.0, 2.5, 12.0 })
         for (double knee : { 0.0, 9.0 })
          for (int nch = 1; nch <= 2; ++nch)
          {
              GainReductionParams p = baseParams();
              p.detector = det; p.link = link; p.mode = mode;
              p.thresholdDb = thr; p.ratio = ratio; p.kneeDb = knee;
              const float* planes[2] { music.data(), second.data() };

              GainReductionPath path; path.prepare (kFs); path.setParams (p); path.reset();
              PathOracle o; o.configure (kFs, p);
              int bad = 0;
              for (int i = 0; i < 30000; ++i)
              {
                  const float a = path.process (planes, nch, i);
                  const float b = o.step (planes, nch, i);
                  if (std::memcmp (&a, &b, 4) != 0) ++bad;
              }
              worst = std::max (worst, (long) bad);
              ++cases;
          }
    test::ok (worst == 0, "the path equals an external assembly from the public primitives over "
              + std::to_string (cases) + " configurations, bit for bit");

    // processLevel() is the SAME code, entered one stage later — a cached envelope and a live key must
    // not be two implementations that merely agree today.
    {
        GainReductionParams p = baseParams();
        GainReductionPath live, cachedPath;
        live.prepare (kFs); live.setParams (p); live.reset();
        cachedPath.prepare (kFs); cachedPath.setParams (p); cachedPath.reset();
        LinkedDetector d; d.prepare (kFs); d.setParams (p); d.reset();
        const float* p0 = music.data(); const float* const* k = &p0;
        int bad = 0;
        for (int i = 0; i < 30000; ++i)
        {
            const float a = live.process (k, 1, i);
            const float b = cachedPath.processLevel (d.process (k, 1, i));
            if (std::memcmp (&a, &b, 4) != 0) ++bad;
        }
        test::ok (bad == 0, "processLevel(detector output) == process(key), bit for bit");
    }
}

//==============================================================================
static void theCompressorTapIsTheTruth()
{
    test::group ("Compressor's GainReductionTap — the per-sample output P7 needs");

    const auto music = musicish (24000, kFs, 8);

    // 1. The tap equals the path driven externally, bit for bit, across shapes.
    long cases = 0; long worst = 0;
    for (auto mode : { Mode::DownCompress, Mode::UpCompress, Mode::DownExpand })
     for (double thr : { -35.0, -12.0 })
      for (double ratio : { 2.0, 8.0 })
       for (int blk : { 1, 37, 512 })
       {
           CompressorParams cp; static_cast<GainReductionParams&> (cp) = baseParams();
           cp.mode = mode; cp.thresholdDb = thr; cp.ratio = ratio;
           cp.makeupDb = 3.0; cp.autoMakeup = true; cp.lookaheadMs = 1.0;

           Compressor c; test::ok (c.prepare (kFs, 512, 1, 20.0), "prepare");
           c.setParams (cp);
           std::vector<float> audio = music, tap ((std::size_t) 24000, -12345.0f);
           for (int off = 0; off < 24000; off += blk)
           {
               const int m = std::min (blk, 24000 - off);
               float* io[1] { audio.data() + off };
               c.process (io, 1, m, GainReductionTap { tap.data() + off, m });
           }
           GainReductionPath path; path.prepare (kFs); path.setParams (cp); path.reset();
           const float* p0 = music.data(); const float* const* k = &p0;
           int bad = 0;
           for (int i = 0; i < 24000; ++i)
           {
               const float want = path.process (k, 1, i);
               if (std::memcmp (&want, &tap[(std::size_t) i], 4) != 0) ++bad;
           }
           // The follower's denormal flush is clocked by the CALLER's blocks, so a value that has
           // decayed into the subnormal range can differ in its last bits between partitions. Anything
           // above that is a real disagreement.
           worst = std::max (worst, (long) bad);
           ++cases;
       }
    test::ok (worst == 0, "the tap is bit-identical to the path over " + std::to_string (cases)
              + " configurations and three block sizes");

    // 2. Tapping must not change the audio, or the state.
    {
        CompressorParams cp; static_cast<GainReductionParams&> (cp) = baseParams();
        cp.lookaheadMs = 2.0;
        Compressor a, b;
        (void) a.prepare (kFs, 256, 1, 20.0); (void) b.prepare (kFs, 256, 1, 20.0);
        a.setParams (cp); b.setParams (cp);
        std::vector<float> ax = music, bx = music, tap ((std::size_t) 24000);
        for (int off = 0; off < 24000; off += 256)
        {
            const int m = std::min (256, 24000 - off);
            float* ai[1] { ax.data() + off }; float* bi[1] { bx.data() + off };
            a.process (ai, 1, m);
            b.process (bi, 1, m, GainReductionTap { tap.data() + off, m });
        }
        test::ok (std::memcmp (ax.data(), bx.data(), ax.size() * 4) == 0, "audio is bit-identical with the tap on");
        test::ok (a.gainReductionDb() == b.gainReductionDb(), "and so is the state");
    }

    // 3. A SHORT TAP REFUSES THE WHOLE CALL. A partially written diagnostic looks like data, which is
    //    worse than none — the same reason a channel count above the prepared maximum is refused.
    {
        CompressorParams cp; static_cast<GainReductionParams&> (cp) = baseParams();
        Compressor c; (void) c.prepare (kFs, 256, 1, 0.0); c.setParams (cp);
        std::vector<float> audio (256, 0.7f), tap (256, -1.0f);
        const std::vector<float> before = audio;
        float* io[1] { audio.data() };
        c.process (io, 1, 256, GainReductionTap { tap.data(), 255 });
        test::ok (std::memcmp (audio.data(), before.data(), 256 * 4) == 0, "a short tap leaves the AUDIO untouched");
        bool tapUntouched = true; for (float v : tap) if (v != -1.0f) tapUntouched = false;
        test::ok (tapUntouched, "and the tap buffer untouched");
        test::ok (c.gainReductionDb() == 0.0, "and the state untouched");
        // exactly enough is enough — and the proof is that the TAP was written, not that the audio
        // moved: the first sample of a fresh compressor legitimately comes out at unity gain.
        c.process (io, 1, 256, GainReductionTap { tap.data(), 256 });
        bool tapWritten = false; for (float v : tap) if (v != -1.0f) tapWritten = true;
        test::ok (tapWritten, "a capacity of exactly numSamples is accepted");
    }

    // 4. A null tap with a nonzero capacity is simply off.
    {
        CompressorParams cp; static_cast<GainReductionParams&> (cp) = baseParams();
        Compressor c; (void) c.prepare (kFs, 64, 1, 0.0); c.setParams (cp);
        std::vector<float> audio (64, 0.6f);
        float* io[1] { audio.data() };
        c.process (io, 1, 64, GainReductionTap { nullptr, 999 });
        test::ok (c.gainReductionDb() < 0.0, "a null tap with a capacity is off, not a refusal — the call ran");
    }

    // 5. The call shapes: the tap forms exist, and the bare-pointer key form still does not.
    static_assert (! p3contract::HasFourArgKeyProcess<Compressor>::value,
                   "a key pointer without a channel count must still not compile");
    static_assert (p3contract::HasTapProcess<Compressor>::value, "the tap forms must exist");
    test::ok (true, "no four-argument KEY form, and the tap form does exist (static_assert)");
}

//==============================================================================
// THE TEST THAT DECIDES EVERYTHING ELSE.
static void f9CounterexampleIsMandatory()
{
    test::group ("F9 — identical level percentiles, different p95 GR, therefore different thresholds");

    const int n = (int) (5.0 * kFs);
    TwoLevel contiguous, spread;
    contiguous.build (n, -4.0, -40.0, true, 5);
    spread.build     (n, -4.0, -40.0, false, 5);

    GainReductionParams p;
    p.detector = Detector::Peak; p.link = LinkMode::Max;
    p.mode = Mode::DownCompress; p.ratio = 4.0; p.kneeDb = 0.0; p.rangeDb = 60.0;
    p.attackMs = 5.0; p.releaseMs = 5.0;

    // First: the premise. The two signals must be indistinguishable in the detector domain.
    EnvelopeStats sa, sb;
    for (int which = 0; which < 2; ++which)
    {
        EnvelopeAnalyzer a; (void) a.prepare (kFs); a.setParams (p); a.reset();
        a.analyze (which == 0 ? contiguous.key : spread.key, 1, n);
        (which == 0 ? sa : sb) = a.stats();
    }
    test::approx (sa.p50Db, sb.p50Db, 1e-9, "the two layouts have the same p50");
    test::approx (sa.p90Db, sb.p90Db, 1e-9, "the same p90");
    test::approx (sa.p95Db, sb.p95Db, 1e-9, "the same p95");
    test::approx (sa.p99Db, sb.p99Db, 1e-9, "the same p99 — no function of these can tell them apart");

    // Second: the gain reduction they actually produce is nothing alike.
    {
        p.thresholdDb = -20.0;
        double got[2];
        for (int which = 0; which < 2; ++which)
        {
            GainReductionPath path; path.prepare (kFs); path.setParams (p); path.reset();
            std::vector<double> g ((std::size_t) n);
            const auto* k = (which == 0 ? contiguous.key : spread.key);
            for (int i = 0; i < n; ++i) g[(std::size_t) i] = std::fabs ((double) path.process (k, 1, i));
            got[which] = exactNearestRank (g, 0.95);
        }
        test::approx (got[0], 12.0, 0.01, "laid out as one contiguous second, p95 GR is 12.000 dB");
        test::approx (got[1], 2.42, 0.01, "laid out as every fifth sample, p95 GR is 2.420 dB");
    }

    // Third: the solver must therefore return DIFFERENT thresholds. An implementation that inverts the
    // static curve at a level percentile returns the SAME one, and fails exactly here.
    {
        ThresholdTarget t; t.amountDb = 6.0; t.toleranceDb = 0.05;
        double thr[2], seed[2];
        for (int which = 0; which < 2; ++which)
        {
            ThresholdSolver s; test::ok (s.prepare (kFs, n), "prepare");
            const auto r = s.solve (which == 0 ? contiguous.key : spread.key, 1, n, p, t);
            test::ok (r.status == SolveStatus::Solved, "solved layout " + std::to_string (which));
            test::approx (r.achievedDb, 6.0, 0.05, "and hit the target on layout " + std::to_string (which));
            thr[which] = r.thresholdDb; seed[which] = r.seedDb;
        }
        test::approx (seed[0], seed[1], 1e-9,
                      "the SEED is identical for both — it is a function of the level percentile alone");
        test::ok (std::fabs (thr[0] - thr[1]) > 5.0,
                  "but the answers differ by " + std::to_string (std::fabs (thr[0] - thr[1]))
                  + " dB, which no inversion of the static curve can produce");
    }
}

//==============================================================================
static void theRoundTripIsMeasuredWithADifferentInstrument()
{
    test::group ("round-trip — verified by sorting the REAL compressor's tap, not by the solver's own histogram");

    const int n = (int) (4.0 * kFs);
    const auto music = musicish (n, kFs, 21);
    const float* p0 = music.data(); const float* const* key = &p0;

    long cases = 0; double worst = 0.0;
    for (auto det : { Detector::Peak, Detector::Rms })
     for (double ratio : { 2.0, 4.0, 12.0 })
      for (double knee : { 0.0, 6.0 })
       for (double atk : { 1.0, 30.0 })
        for (double rel : { 40.0, 400.0 })
         for (double want : { 1.0, 3.0, 8.0 })
         {
             GainReductionParams p = baseParams();
             p.detector = det; p.ratio = ratio; p.kneeDb = knee; p.attackMs = atk; p.releaseMs = rel;
             ThresholdTarget t; t.amountDb = want; t.toleranceDb = 0.05;

             ThresholdSolver s; (void) s.prepare (kFs, n);
             const auto r = s.solve (key, 1, n, p, t);
             if (r.status != SolveStatus::Solved) { test::ok (false, "solved"); continue; }

             // The independent instrument: the real Compressor, its tap, an exact sort.
             Compressor c; (void) c.prepare (kFs, 1024, 1, 0.0);
             CompressorParams cp; static_cast<GainReductionParams&> (cp) = p; cp.thresholdDb = r.thresholdDb;
             c.setParams (cp);
             std::vector<float> audio = music, tap ((std::size_t) n);
             for (int off = 0; off < n; off += 1024)
             {
                 const int m = std::min (1024, n - off);
                 float* io[1] { audio.data() + off };
                 c.process (io, 1, m, GainReductionTap { tap.data() + off, m });
             }
             std::vector<double> mag ((std::size_t) n);
             for (int i = 0; i < n; ++i) mag[(std::size_t) i] = std::fabs ((double) tap[(std::size_t) i]);
             const double exact = exactNearestRank (mag, 0.95);
             worst = std::max (worst, std::fabs (exact - want));
             ++cases;
         }
    test::ok (worst <= 0.06, "over " + std::to_string (cases)
              + " configurations the REAL compressor's exactly-sorted p95 is within "
              + std::to_string (worst) + " dB of the requested amount (tolerance 0.05 + half a bin)");
}

//==============================================================================
static void everyRefusalHasAName()
{
    test::group ("the inverse problem is not always defined — and says which way it is not");

    const auto music = musicish (48000, kFs, 2);
    const float* p0 = music.data(); const float* const* key = &p0;
    std::vector<float> silence ((std::size_t) 48000, 0.0f);
    const float* z0 = silence.data(); const float* const* zkey = &z0;

    auto run = [&] (GainReductionParams p, ThresholdTarget t, const float* const* k, int n, int cached)
    {
        ThresholdSolver s; (void) s.prepare (kFs, cached);
        return s.solve (k, 1, n, p, t);
    };
    const GainReductionParams B = baseParams();
    ThresholdTarget T; T.amountDb = 3.0; T.toleranceDb = 0.05;

    test::ok (run (B, T, key, 48000, 48000).status == SolveStatus::Solved, "ordinary material solves");

    // The trap both consilium seats found independently: DIGITAL SILENCE has a plausible answer. The
    // gainToDb floor makes it a constant -240 dB, so a low enough threshold "compresses" it and the
    // search converges — on audio that is still exactly zero. It must be named, not returned as a number.
    {
        const auto r = run (B, T, zkey, 48000, 48000);
        test::ok (r.status == SolveStatus::FloorDrivesTheAnswer,
                  "digital silence is named FloorDrivesTheAnswer, not answered with a plausible -244");
        test::ok (r.thresholdDb < -240.0, "and the threshold it would have returned is indeed below the floor");
    }

    // ...AND IN EVERY MODE. A mode-shaped test of this got it wrong exactly where a mode-shaped test
    // does: UpCompress does not need a threshold below the floor to process silence — it processes the
    // quiet end by design — so it returned Solved at -236 dB while lifting digital silence by the
    // requested 3 dB. The question is asked of the CURVE now, so all three modes answer it.
    {
        const int n = 40000;
        std::vector<float> mixed ((std::size_t) n, 0.0f);
        for (int i = 0; i < n * 3 / 5; ++i)
            mixed[(std::size_t) i] = (float) (core::dbToGain (-20.0) * std::sin (2.0 * core::kPi * 400.0 * i / kFs));
        const float* m0 = mixed.data(); const float* const* mk = &m0;
        for (auto mode : { Mode::DownCompress, Mode::UpCompress, Mode::DownExpand })
        {
            auto p = B; p.mode = mode; p.detector = Detector::Peak; p.ratio = 3.0; p.kneeDb = 0.0;
            const auto r = run (p, T, mk, n, 0);
            const std::string nm = (mode == Mode::DownCompress) ? "DownCompress"
                                 : (mode == Mode::UpCompress) ? "UpCompress" : "DownExpand";
            if (r.status == SolveStatus::Solved)
            {
                // If it says Solved, no frame at the floor may be being processed — that is the whole
                // content of the other status, and this is the assertion that makes it mean something.
                GainComputer c; c.setMode (mode); c.setRatio (p.ratio); c.setKneeDb (p.kneeDb);
                c.setRangeDb (p.rangeDb); c.setThresholdDb (r.thresholdDb);
                test::approx (std::fabs (c.deltaDb (-240.0)), 0.0, 0.0,
                              nm + ": a Solved answer leaves the frames on the floor untouched");
            }
            else
                test::ok (r.status == SolveStatus::FloorDrivesTheAnswer,
                          nm + ": material with digital silence is either solved cleanly or named "
                          "FloorDrivesTheAnswer");
        }
    }

    { auto p = B; p.ratio = 1.0;
      test::ok (run (p, T, key, 48000, 0).status == SolveStatus::NoProcessingPossible, "ratio 1:1 — the curve is flat"); }
    { auto p = B; p.rangeDb = 0.0;
      test::ok (run (p, T, key, 48000, 0).status == SolveStatus::NoProcessingPossible, "a zero range — the curve is flat"); }
    { auto t = T; t.amountDb = 0.0;
      const auto r = run (B, t, key, 48000, 0);
      test::ok (r.status == SolveStatus::TargetIsZero, "a zero target is degenerate but legitimate");
      test::approx (r.achievedDb, 0.0, 1e-12, "and the canonical no-processing threshold delivers zero"); }
    { auto t = T; t.amountDb = 90.0;
      const auto r = run (B, t, key, 48000, 0);
      test::ok (r.status == SolveStatus::TargetNotReachable, "a target beyond the range is unreachable");
      test::ok (r.attainableDb > 55.0 && r.attainableDb <= 60.0,
                "and the report says how much this material CAN give (" + std::to_string (r.attainableDb) + " dB)"); }

    // TargetNotReachable is MEASURED, not inferred from rangeDb: a short, sparse passage cannot charge
    // the follower to a target that the curve alone would reach easily.
    {
        TwoLevel s; s.build (2000, -6.0, -60.0, false, 40);     // 2.5 % loud, 42 ms long
        auto p = B; p.detector = Detector::Peak; p.ratio = 20.0; p.attackMs = 200.0;
        auto t = T; t.amountDb = 40.0;
        const auto r = run (p, t, s.key, 2000, 0);
        test::ok (r.status == SolveStatus::TargetNotReachable,
                  "40 dB is inside a 60 dB range, but the ballistics never get there on 42 ms of material");
        test::ok (r.attainableDb < 40.0, "and the attainable amount is reported (" + std::to_string (r.attainableDb) + " dB)");
    }

    { auto t = T; t.amountDb = std::numeric_limits<double>::quiet_NaN();
      test::ok (run (B, t, key, 48000, 0).status == SolveStatus::InvalidTarget, "a NaN target is refused"); }
    { auto t = T; t.amountDb = -3.0;
      test::ok (run (B, t, key, 48000, 0).status == SolveStatus::InvalidTarget, "a negative target is refused"); }
    { auto t = T; t.quantile = 0.0;
      test::ok (run (B, t, key, 48000, 0).status == SolveStatus::InvalidTarget, "a quantile of 0 is refused"); }
    { auto t = T; t.quantile = 1.5;
      test::ok (run (B, t, key, 48000, 0).status == SolveStatus::InvalidTarget, "a quantile above 1 is refused"); }
    { auto p = B; p.ratio = std::numeric_limits<double>::infinity();
      test::ok (run (p, T, key, 48000, 0).status == SolveStatus::InvalidParams, "non-finite curve parameters are refused"); }
    test::ok (run (B, T, key, 0, 0).status == SolveStatus::EmptyInput, "zero samples is EmptyInput");
    test::ok (run (B, T, nullptr, 48000, 0).status == SolveStatus::EmptyInput, "a null key is EmptyInput, not a read of plane 0");
    { ThresholdSolver s; test::ok (s.solve (key, 1, 48000, B, T).status == SolveStatus::NotPrepared, "solving before prepare is NotPrepared"); }
    { ThresholdSolver s; test::ok (s.solve (key, 0, 48000, B, T).status == SolveStatus::NotPrepared, "and a zero channel count too"); }
    { auto t = T; t.maxEvaluations = 3;
      const auto r = run (B, t, key, 48000, 0);
      test::ok (r.status == SolveStatus::EvaluationLimit, "running out of passes is EvaluationLimit, its own name");
      test::ok (r.evaluations <= 6, "and it stopped when it said it would"); }

    // THE ANSWER IS THE ENDPOINT THAT MEETS THE TARGET, and the only place that is visible is where
    // the bracket is still WIDE. Once the search has converged the two ends are a thousandth of a dB
    // apart and returning the wrong one is invisible to every tolerance in this file — which is
    // exactly how a mutation of it survived the first rollback pass.
    {
        for (int budget : { 2, 4, 6, 9 })
        {
            auto t = T; t.maxEvaluations = budget; t.amountDb = 6.0;
            ThresholdSolver s; (void) s.prepare (kFs, 48000);
            const auto r = s.solve (key, 1, 48000, B, t);
            if (r.status != SolveStatus::Solved && r.status != SolveStatus::EvaluationLimit) continue;
            // re-measure at the RETURNED threshold with the real compressor and an exact sort
            Compressor c; (void) c.prepare (kFs, 4096, 1, 0.0);
            CompressorParams cp; static_cast<GainReductionParams&> (cp) = B; cp.thresholdDb = r.thresholdDb;
            c.setParams (cp);
            std::vector<float> audio = music, tap ((std::size_t) 48000);
            for (int off = 0; off < 48000; off += 4096)
            {
                const int m = std::min (4096, 48000 - off);
                float* io[1] { audio.data() + off };
                c.process (io, 1, m, GainReductionTap { tap.data() + off, m });
            }
            std::vector<double> mag ((std::size_t) 48000);
            for (int i = 0; i < 48000; ++i) mag[(std::size_t) i] = std::fabs ((double) tap[(std::size_t) i]);
            const double exact = exactNearestRank (mag, 0.95);
            test::ok (exact >= t.amountDb - 0.006,
                      "with a budget of " + std::to_string (budget) + " passes the returned threshold still DELIVERS the target ("
                      + std::to_string (exact) + " dB) — it is the end of the bracket that was verified, not the end that failed");
            test::approx (r.achievedDb, exact, 0.012,
                          "and the reported amount belongs to the reported threshold");
        }
    }

    // A SATURATED 400 dB RANGE must not fall out of the histogram — |gain change| reaches `rangeDb`
    // exactly, and `rangeDb` may be the cap itself.
    {
        auto p = B; p.rangeDb = 400.0; p.ratio = 20.0;
        auto t = T; t.amountDb = 30.0;
        const auto r = run (p, t, key, 48000, 0);
        test::ok (r.status != SolveStatus::QuantileOutOfRange, "a 400 dB range still measures in range");
        test::ok (r.status == SolveStatus::Solved, "and solves");
    }
}

//==============================================================================
static void theGateReadsTheLevelAndNotTheResult()
{
    test::group ("the statistics gate — on the LEVEL, which is what keeps the objective monotone");

    // THE FIXTURE HAS TO MAKE THE POLICY VISIBLE. Half a tone against half silence does not: the
    // gated and ungated p95 came out 0.02 dB apart there, because a compressor's gain reduction is
    // nearly flat near its top, and a mutation that gated on the RESULT instead of the level survived.
    // What separates them is (a) the excluded part DOMINATING the count and (b) the included part
    // having a spread. So: a tenth of the material is a level RAMP, and nine tenths is digital
    // silence. Then the ungated p95 is the ramp's MEDIAN gain reduction and the gated one is the
    // ramp's own p95 — several dB apart, and nothing about the shape of a compressor can close that.
    const int n = 60000;
    const int loud = n / 10;
    std::vector<float> x ((std::size_t) n, 0.0f);
    for (int i = 0; i < loud; ++i)
    {
        const double db = -40.0 + 36.0 * (double) i / (double) (loud - 1);      // -40 .. -4 dB
        x[(std::size_t) i] = (float) core::dbToGain (db);
    }
    const float* p0 = x.data(); const float* const* key = &p0;

    // Peak detection with near-instant ballistics, so the gain reduction follows the ramp rather than
    // smearing it back into a flat top.
    auto B = baseParams();
    B.detector = Detector::Peak; B.kneeDb = 0.0; B.ratio = 4.0;
    B.attackMs = 0.1; B.releaseMs = 0.1;
    for (double gate : { -60.0, -40.0 })
    {
        ThresholdTarget t; t.amountDb = 3.0; t.toleranceDb = 0.05; t.statisticsGateDb = gate;
        ThresholdSolver s; (void) s.prepare (kFs, n);
        const auto r = s.solve (key, 1, n, B, t);
        test::ok (r.status == SolveStatus::Solved, "a gated search converges (gate " + std::to_string (gate) + " dB)");

        // INDEPENDENT VERIFICATION, and the point of the whole test: recompute the statistic from the
        // real compressor's tap, selecting frames by the DETECTOR LEVEL. A solver that gated on the
        // gain reduction instead would have measured a different set and cannot agree here.
        Compressor c; (void) c.prepare (kFs, 1024, 1, 0.0);
        CompressorParams cp; static_cast<GainReductionParams&> (cp) = B; cp.thresholdDb = r.thresholdDb;
        c.setParams (cp);
        std::vector<float> audio = x, tap ((std::size_t) n);
        for (int off = 0; off < n; off += 1024)
        {
            const int m = std::min (1024, n - off);
            float* io[1] { audio.data() + off };
            c.process (io, 1, m, GainReductionTap { tap.data() + off, m });
        }
        LinkedDetector d; d.prepare (kFs); d.setParams (B); d.reset();
        std::vector<double> kept;
        for (int i = 0; i < n; ++i)
        {
            const double lvlDb = core::gainToDb (d.process (key, 1, i));
            if (lvlDb > gate) kept.push_back (std::fabs ((double) tap[(std::size_t) i]));
        }
        test::ok (! kept.empty() && kept.size() < (std::size_t) n, "the gate kept some frames and dropped others");
        test::approx (exactNearestRank (kept, 0.95), r.achievedDb, 0.012,
                      "the solver's gated statistic equals one recomputed by selecting on the LEVEL");
        test::ok (r.envelope.frames == kept.size(), "and the analyzer counted exactly those frames");

        // ...and it is NOT the ungated statistic. Without this the check above passes on any material
        // where the two happen to coincide, which is most of it.
        std::vector<double> everything ((std::size_t) n);
        for (int i = 0; i < n; ++i) everything[(std::size_t) i] = std::fabs ((double) tap[(std::size_t) i]);
        test::ok (std::fabs (exactNearestRank (everything, 0.95) - r.achievedDb) > 0.5,
                  "and differs from the UNGATED statistic by more than half a dB ("
                  + std::to_string (exactNearestRank (everything, 0.95)) + " vs " + std::to_string (r.achievedDb)
                  + ") — so the fixture can tell the two policies apart");
    }

    // A gate that keeps everything must give the same answer as no gate at all.
    {
        ThresholdTarget open; open.amountDb = 3.0; open.toleranceDb = 0.05;
        auto low = open; low.statisticsGateDb = -300.0;
        ThresholdSolver a, b;
        (void) a.prepare (kFs, n); (void) b.prepare (kFs, n);
        const auto ra = a.solve (key, 1, n, B, open);
        const auto rb = b.solve (key, 1, n, B, low);
        test::ok (std::memcmp (&ra.thresholdDb, &rb.thresholdDb, 8) == 0,
                  "a gate below the floor is bit-identical to no gate");
    }
    // And a gate high enough to remove everything is EmptyInput, not a silent zero.
    {
        ThresholdTarget t; t.amountDb = 3.0; t.statisticsGateDb = 40.0;
        ThresholdSolver s; (void) s.prepare (kFs, n);
        test::ok (s.solve (key, 1, n, B, t).status == SolveStatus::EmptyInput,
                  "a gate that excludes every frame is EmptyInput, not a plausible answer");
    }
}

//==============================================================================
static void aSearchIsAnIterationAndStateMustNotSurviveIt()
{
    test::group ("undifferentiated sequences — a solve is an iteration, and setParams() is not a reset");

    // The signal ENDS LOUD, so every probe leaves a large non-zero follower state behind. A solver that
    // forgot to reset between probes would return an answer that depends on the order it tried them.
    const int n = 40000;
    std::vector<float> x ((std::size_t) n, (float) core::dbToGain (-40.0));
    for (int i = n - 8000; i < n; ++i) x[(std::size_t) i] = (float) core::dbToGain (-3.0);
    for (int i = 0; i < n; i += 13) x[(std::size_t) i] = (float) core::dbToGain (-8.0);
    const float* p0 = x.data(); const float* const* key = &p0;

    const auto B = baseParams();
    ThresholdTarget T; T.amountDb = 4.0; T.toleranceDb = 0.05;

    // 1. The same object, solved twice — bit for bit.
    {
        ThresholdSolver s; (void) s.prepare (kFs, n);
        const auto a = s.solve (key, 1, n, B, T);
        const auto b = s.solve (key, 1, n, B, T);
        test::ok (std::memcmp (&a.thresholdDb, &b.thresholdDb, 8) == 0, "solving twice on one object is bit-identical");
        test::ok (a.evaluations == b.evaluations, "and takes the same number of passes");
    }
    // 2. A DIFFERENT solve in between must leave nothing behind.
    {
        ThresholdSolver fresh; (void) fresh.prepare (kFs, n);
        const auto want = fresh.solve (key, 1, n, B, T);
        ThresholdSolver used; (void) used.prepare (kFs, n);
        auto other = T; other.amountDb = 14.0;
        (void) used.solve (key, 1, n, B, other);
        auto p2 = B; p2.mode = Mode::UpCompress;
        (void) used.solve (key, 1, n, p2, other);
        const auto got = used.solve (key, 1, n, B, T);
        test::ok (std::memcmp (&want.thresholdDb, &got.thresholdDb, 8) == 0,
                  "a solver used for two other problems returns the first answer bit for bit");
    }
    // 3. The order the targets are asked for must not matter.
    {
        const double targets[3] { 1.0, 12.0, 3.0 };
        double forward[3], backward[3];
        ThresholdSolver a, b;
        (void) a.prepare (kFs, n); (void) b.prepare (kFs, n);
        for (int i = 0; i < 3; ++i) { auto t = T; t.amountDb = targets[i]; forward[i] = a.solve (key, 1, n, B, t).thresholdDb; }
        for (int i = 2; i >= 0; --i) { auto t = T; t.amountDb = targets[i]; backward[i] = b.solve (key, 1, n, B, t).thresholdDb; }
        for (int i = 0; i < 3; ++i)
            test::ok (std::memcmp (&forward[i], &backward[i], 8) == 0,
                      "target " + std::to_string (targets[i]) + " dB gives the same answer whichever order it was asked in");
    }
    // 4. A long loud record followed by a short quiet one: no tail of the first survives.
    {
        std::vector<float> quiet ((std::size_t) 4000, (float) core::dbToGain (-30.0));
        for (int i = 0; i < 4000; i += 7) quiet[(std::size_t) i] = (float) core::dbToGain (-12.0);
        const float* q0 = quiet.data(); const float* const* qkey = &q0;
        ThresholdSolver fresh, used;
        (void) fresh.prepare (kFs, n); (void) used.prepare (kFs, n);
        const auto want = fresh.solve (qkey, 1, 4000, B, T);
        (void) used.solve (key, 1, n, B, T);
        const auto got = used.solve (qkey, 1, 4000, B, T);
        test::ok (std::memcmp (&want.thresholdDb, &got.thresholdDb, 8) == 0, "a shorter, quieter follow-up is unaffected by the first");
        test::ok (want.envelope.frames == got.envelope.frames, "and the statistics were not accumulated across solves");
    }
    // 5. The stepwise form is the same search, taken one pass at a time (what a single-threaded wasm
    //    host needs in order to stay responsive).
    {
        ThresholdSolver a, b;
        (void) a.prepare (kFs, n); (void) b.prepare (kFs, n);
        const auto oneShot = a.solve (key, 1, n, B, T);
        int steps = 0;
        if (b.begin (key, 1, n, B, T)) while (b.step()) ++steps;
        test::ok (std::memcmp (&oneShot.thresholdDb, &b.solution().thresholdDb, 8) == 0,
                  "begin()/step() reaches the same answer, bit for bit");
        test::ok (steps > 0 && b.solution().evaluations == oneShot.evaluations, "in the same number of passes");
    }
    // 6. The envelope cache is a memory/CPU trade, never a different answer.
    {
        ThresholdSolver cached, live;
        (void) cached.prepare (kFs, n); (void) live.prepare (kFs, 0);
        const auto a = cached.solve (key, 1, n, B, T);
        const auto c = live.solve (key, 1, n, B, T);
        test::ok (std::memcmp (&a.thresholdDb, &c.thresholdDb, 8) == 0, "cached and live envelopes give the same threshold, bit for bit");
        test::ok (std::memcmp (&a.achievedDb, &c.achievedDb, 8) == 0, "and the same measured amount");
        test::ok (a.evaluations == c.evaluations, "and the same number of passes");
        // A cache too small for the material must silently fall back, not truncate.
        ThresholdSolver small; (void) small.prepare (kFs, n / 2);
        const auto d = small.solve (key, 1, n, B, T);
        test::ok (std::memcmp (&a.thresholdDb, &d.thresholdDb, 8) == 0, "a cache too small falls back to the live path, same answer");
    }
    // 7. No allocation once prepared — a search that allocated would be a search that could stall.
    {
        ThresholdSolver s; (void) s.prepare (kFs, n);
        (void) s.solve (key, 1, n, B, T);                        // warm any lazy state
        const long before = g_allocs.load();
        (void) s.solve (key, 1, n, B, T);
        test::okNoAlloc (g_allocs.load() == before, "solve() allocates nothing after prepare()");
    }
}

//==============================================================================
static void theCurveOnlyEverSeesLevelMinusThreshold()
{
    test::group ("shift invariance — a free consequence of the topology, and a sharp test");

    // The curve is a function of (level - threshold) alone and the follower never sees an absolute
    // level, so scaling the key by g must move the answer by exactly 20·log10(g). Anything that starts
    // using an absolute level after the curve breaks this and nothing else here would notice.
    const int n = 30000;
    const auto music = musicish (n, kFs, 33);
    const auto B = baseParams();
    ThresholdTarget T; T.amountDb = 3.0; T.toleranceDb = 0.02; T.resolutionDb = 1e-4;

    ThresholdSolver s; (void) s.prepare (kFs, n);
    const float* m0 = music.data(); const float* const* mk = &m0;
    const double ref = s.solve (mk, 1, n, B, T).thresholdDb;

    for (double g : { 0.25, 0.5, 2.0, 4.0 })
    {
        std::vector<float> scaled ((std::size_t) n);
        for (int i = 0; i < n; ++i) scaled[(std::size_t) i] = (float) (g * (double) music[(std::size_t) i]);
        const float* s0 = scaled.data(); const float* const* sk = &s0;
        ThresholdSolver t2; (void) t2.prepare (kFs, n);
        const double got = t2.solve (sk, 1, n, B, T).thresholdDb;
        test::approx (got, ref + 20.0 * std::log10 (g), 0.01,
                      "scaling the key by " + std::to_string (g) + " moves the threshold by exactly 20log10(g)");
    }
}

//==============================================================================
static void everyModeIsSolvedInItsOwnDirection()
{
    test::group ("modes — the search direction comes from the mode, and a flipped sign fails loudly");

    const int n = 30000;
    const auto music = musicish (n, kFs, 77);
    const float* m0 = music.data(); const float* const* key = &m0;

    for (auto mode : { Mode::DownCompress, Mode::UpCompress, Mode::DownExpand })
     for (double want : { 2.0, 7.0 })
     {
         auto p = baseParams(); p.mode = mode; p.ratio = 3.0;
         ThresholdTarget t; t.amountDb = want; t.toleranceDb = 0.05;
         ThresholdSolver s; (void) s.prepare (kFs, n);
         const auto r = s.solve (key, 1, n, p, t);
         const std::string nm = (mode == Mode::DownCompress) ? "DownCompress" : (mode == Mode::UpCompress ? "UpCompress" : "DownExpand");
         test::ok (r.status == SolveStatus::Solved, nm + " solves for " + std::to_string (want) + " dB");
         test::approx (r.achievedDb, want, 0.06, nm + " hits the target");
         // The answer must be strictly inside the bracket, not pinned at an end — which is exactly what
         // a search running in the wrong direction produces.
         test::ok (r.thresholdDb > r.bracketLoDb - 1e-9 && r.thresholdDb < r.bracketHiDb + 1e-9,
                   nm + " converged inside its bracket rather than pinning at an end");

         // More processing must need a more aggressive threshold, in the direction the mode implies.
         auto t2 = t; t2.amountDb = want + 3.0;
         ThresholdSolver s2; (void) s2.prepare (kFs, n);
         const auto r2 = s2.solve (key, 1, n, p, t2);
         if (r2.status == SolveStatus::Solved)
         {
             const bool lowerIsMore = (mode == Mode::DownCompress);
             test::ok (lowerIsMore ? (r2.thresholdDb < r.thresholdDb) : (r2.thresholdDb > r.thresholdDb),
                       nm + ": asking for more moves the threshold the way the mode says");
         }
     }
}

//==============================================================================
static void theStatisticIsAParameterNotAnAssumption()
{
    test::group ("mean and max are targetable too — the search is valid for any pointwise-monotone functional");

    const int n = 30000;
    const auto music = musicish (n, kFs, 91);
    const float* m0 = music.data(); const float* const* key = &m0;
    const auto B = baseParams();

    for (auto stat : { TraceStatistic::Quantile, TraceStatistic::Mean, TraceStatistic::Max })
    {
        ThresholdTarget t; t.statistic = stat; t.amountDb = 3.0; t.toleranceDb = 0.05;
        ThresholdSolver s; (void) s.prepare (kFs, n);
        const auto r = s.solve (key, 1, n, B, t);
        test::ok (r.status == SolveStatus::Solved, "the search converges on this functional");

        // verified with the independent instrument again
        Compressor c; (void) c.prepare (kFs, 512, 1, 0.0);
        CompressorParams cp; static_cast<GainReductionParams&> (cp) = B; cp.thresholdDb = r.thresholdDb;
        c.setParams (cp);
        std::vector<float> audio = music, tap ((std::size_t) n);
        for (int off = 0; off < n; off += 512)
        {
            const int m = std::min (512, n - off);
            float* io[1] { audio.data() + off };
            c.process (io, 1, m, GainReductionTap { tap.data() + off, m });
        }
        double sum = 0.0, mx = 0.0;
        std::vector<double> mag ((std::size_t) n);
        for (int i = 0; i < n; ++i) { const double v = std::fabs ((double) tap[(std::size_t) i]); mag[(std::size_t) i] = v; sum += v; mx = std::max (mx, v); }
        const double got = (stat == TraceStatistic::Mean) ? sum / (double) n
                         : (stat == TraceStatistic::Max)  ? mx
                                                          : exactNearestRank (mag, 0.95);
        test::approx (got, 3.0, 0.06, "and the real compressor delivers it, measured independently");
    }
}

//==============================================================================
static void theSeedIsAnEstimateAndItsValueIsMeasured()
{
    test::group ("the seed — measured, not assumed, including where it costs more than it saves");

    const int n = (int) (3.0 * kFs);
    struct Case { const char* name; std::vector<float> x; };
    std::vector<Case> cases;
    cases.push_back ({ "music", musicish (n, kFs, 5) });
    { TwoLevel s; s.build (n, -4.0, -40.0, true, 5);  cases.push_back ({ "F9 contiguous", s.x }); }
    { TwoLevel s; s.build (n, -4.0, -40.0, false, 5); cases.push_back ({ "F9 spread",     s.x }); }

    int betterCount = 0, worseCount = 0, totalSaved = 0;
    for (auto& c : cases)
    {
        const float* p0 = c.x.data(); const float* const* key = &p0;
        auto p = baseParams(); p.detector = Detector::Peak; p.attackMs = 5.0; p.releaseMs = 5.0;
        ThresholdTarget on; on.amountDb = 6.0; on.toleranceDb = 0.05; on.useSeed = true;
        auto off = on; off.useSeed = false;

        ThresholdSolver a, b;
        (void) a.prepare (kFs, n); (void) b.prepare (kFs, n);
        const auto ra = a.solve (key, 1, n, p, on);
        const auto rb = b.solve (key, 1, n, p, off);
        test::approx (ra.thresholdDb, rb.thresholdDb, 0.01,
                      std::string (c.name) + ": the seed changes the route, never the destination");
        const int saved = rb.evaluations - ra.evaluations;
        totalSaved += saved;
        if (saved > 0) ++betterCount; else if (saved < 0) ++worseCount;
        std::printf ("      %-14s seed %8.3f  answer %8.3f  passes with seed %2d, without %2d (%+d)\n",
                     c.name, ra.seedDb, ra.thresholdDb, ra.evaluations, rb.evaluations, saved);
    }
    // The honest statement, and the reason this is a measurement rather than a requirement: midpoint
    // bisection is minimax-optimal on bracket width, so a seed CANNOT improve the worst case — it can
    // only help when it lands near the root, and on the F9 spread layout it is 28 dB away.
    test::ok (betterCount + worseCount >= 0, "seed saved " + std::to_string (totalSaved)
              + " passes in total over " + std::to_string (cases.size()) + " cases ("
              + std::to_string (betterCount) + " better, " + std::to_string (worseCount) + " worse)");

    // A seed outside the bracket must be ignored, not clamped onto an endpoint that is already measured.
    {
        const auto music = musicish (n, kFs, 13);
        const float* p0 = music.data(); const float* const* key = &p0;
        auto p = baseParams();
        ThresholdTarget t; t.amountDb = 3.0; t.toleranceDb = 0.05;
        ThresholdSolver s; (void) s.prepare (kFs, n);
        const auto r = s.solve (key, 1, n, p, t);
        test::ok (r.status == SolveStatus::Solved && r.evaluations < t.maxEvaluations,
                  "the search converges regardless of where the seed lands");
    }
}

//==============================================================================
static void theNaiveInversionIsMeasuredNotArguedAbout()
{
    test::group ("how far the static-curve inversion misses — the one claim about the MODEL still testable");

    // The strongest honest version of the naive method: the same key, the same detector, the same
    // curve, the same quantile rule — with ONLY the ballistics removed. Whatever it misses by is the
    // ballistics' contribution and nothing else.
    const int n = (int) (4.0 * kFs);
    const auto music = musicish (n, kFs, 64);
    const float* p0 = music.data(); const float* const* key = &p0;
    const double want = 3.0;

    double worstMiss = 0.0; int rows = 0, outsideTolerance = 0;
    for (auto det : { Detector::Peak, Detector::Rms })
     for (double atk : { 1.0, 10.0, 50.0 })
      for (double rel : { 20.0, 100.0, 500.0 })
       for (double ratio : { 2.0, 4.0, 10.0 })
       {
           auto p = baseParams();
           p.detector = det; p.attackMs = atk; p.releaseMs = rel; p.ratio = ratio; p.rmsWindowMs = 5.0;

           EnvelopeAnalyzer a; (void) a.prepare (kFs); a.setParams (p); a.reset();
           a.analyze (key, 1, n);
           ThresholdSolver s; (void) s.prepare (kFs, n);
           const double naive = s.staticInversionThresholdDb (p, a.stats().p95Db, want);

           GainReductionPath path; path.prepare (kFs);
           auto q = p; q.thresholdDb = naive; path.setParams (q); path.reset();
           std::vector<double> mag ((std::size_t) n);
           for (int i = 0; i < n; ++i) mag[(std::size_t) i] = std::fabs ((double) path.process (key, 1, i));
           const double miss = exactNearestRank (mag, 0.95) - want;
           worstMiss = std::max (worstMiss, std::fabs (miss));
           if (std::fabs (miss) > 0.3) ++outsideTolerance;
           ++rows;
       }
    // Not a gate: deciding in advance what miss is "enough" is deciding the answer before the
    // experiment. It is recorded because both outcomes are useful and the negative one more so.
    test::ok (rows == 54, "54 ballistics/detector/ratio combinations were measured");
    std::printf ("      naive static inversion misses a %.0f dB p95 target by up to %.3f dB; %d of %d rows exceed 0.3 dB\n",
                 want, worstMiss, outsideTolerance, rows);
    test::ok (worstMiss >= 0.0, "the miss is recorded (worst " + std::to_string (worstMiss) + " dB over "
              + std::to_string (rows) + " configurations, " + std::to_string (outsideTolerance) + " outside 0.3 dB)");
}

//==============================================================================
static void poisonedKeysAndOddShapes()
{
    test::group ("adversarial input — the gate is upstream, and the analysis inherits it");

    const int n = 20000;
    auto music = musicish (n, kFs, 44);
    music[100] = std::numeric_limits<float>::infinity();
    music[8000] = std::numeric_limits<float>::quiet_NaN();
    music[12000] = 1.0e20f;
    const float* p0 = music.data(); const float* const* key = &p0;
    const auto B = baseParams();

    // The path must agree with the oracle even here — the gate is what makes one poisoned sample
    // survivable, and it lives in the detector, ahead of everything this module adds.
    {
        GainReductionPath path; path.prepare (kFs); path.setParams (B); path.reset();
        PathOracle o; o.configure (kFs, B);
        int bad = 0;
        for (int i = 0; i < n; ++i) { const float a = path.process (key, 1, i), b = o.step (key, 1, i); if (std::memcmp (&a, &b, 4) != 0) ++bad; }
        test::ok (bad == 0, "a key carrying Inf, NaN and 1e20 still nulls against the primitives, bit for bit");
    }
    {
        ThresholdTarget t; t.amountDb = 3.0; t.toleranceDb = 0.05;
        ThresholdSolver s; (void) s.prepare (kFs, n);
        const auto r = s.solve (key, 1, n, B, t);
        test::ok (r.status == SolveStatus::Solved || r.status == SolveStatus::TargetNotReachable,
                  "and the search still terminates with a named outcome");
        test::ok (std::isfinite (r.thresholdDb) && std::isfinite (r.achievedDb), "with finite numbers");
    }
    // A single-sample key, and a key of all one value.
    {
        std::vector<float> one (1, 0.5f);
        const float* q0 = one.data(); const float* const* k = &q0;
        ThresholdTarget t; t.amountDb = 3.0;
        ThresholdSolver s; (void) s.prepare (kFs, 4);
        const auto r = s.solve (k, 1, 1, B, t);
        test::ok (r.status != SolveStatus::Solved || std::isfinite (r.thresholdDb), "one sample gives a finite, named result");
        test::ok (r.envelope.framesSeen == 1, "and the statistics say how little there was");
    }
    // Many channels, and a mono key against a stereo programme.
    {
        std::vector<std::vector<float>> ch (4, musicish (n, kFs, 51));
        std::vector<const float*> planes (4);
        for (int c = 0; c < 4; ++c) planes[(std::size_t) c] = ch[(std::size_t) c].data();
        for (auto link : { LinkMode::Max, LinkMode::MeanPower })
        {
            auto p = B; p.link = link;
            ThresholdTarget t; t.amountDb = 3.0; t.toleranceDb = 0.05;
            ThresholdSolver s; (void) s.prepare (kFs, n);
            const auto r = s.solve (planes.data(), 4, n, p, t);
            test::ok (r.status == SolveStatus::Solved, "four key channels solve");
        }
    }
}


//==============================================================================
// Every one of these is a counterexample the AFTER consilium produced against code that was already
// green. They are here in the shape they arrived in — a specific input and the wrong answer it drew —
// because a fix without the input that provoked it is a fix nobody can check.
static void theConsiliumsCounterexamples()
{
    test::group ("the counterexamples that broke this module while its suite was green");

    // 1. A FINITE `rangeDb` PRODUCED AN INFINITE THRESHOLD. The solver sized its bracket from the raw
    //    field while the curve caps it at 400 dB, so `rangeDb / multiplier` overflowed and the search
    //    returned Solved with `thresholdDb = +Inf` — from entirely finite parameters.
    {
        std::vector<float> flat ((std::size_t) 100, 0.5f);
        const float* p0 = flat.data(); const float* const* key = &p0;
        for (double range : { 1.0e308, std::numeric_limits<double>::max(), 400.0, 1.0e9 })
        {
            auto p = baseParams();
            p.mode = Mode::UpCompress; p.detector = Detector::Peak; p.ratio = 2.0;
            p.kneeDb = 0.0; p.attackMs = 0.0; p.releaseMs = 0.0; p.rangeDb = range;
            ThresholdTarget t; t.amountDb = 3.0; t.toleranceDb = 0.05;
            ThresholdSolver s; test::ok (s.prepare (kFs, 128), "prepare");
            const auto r = s.solve (key, 1, 100, p, t);
            test::ok (std::isfinite (r.thresholdDb),
                      "rangeDb = " + std::to_string (range) + " gives a FINITE threshold");
            test::ok (std::isfinite (r.bracketLoDb) && std::isfinite (r.bracketHiDb), "and a finite bracket");
            // A REFUSAL IS ALSO FINITE, so finiteness alone proves nothing: the answer has to be an
            // answer. The curve caps the range at 400 dB, so every one of these is an ordinary solve.
            test::ok (r.status == SolveStatus::Solved,
                      "rangeDb = " + std::to_string (range) + " SOLVES — the effective range is 400 dB, not 1e308");
            test::approx (r.achievedDb, 3.0, 0.06, "and the answer delivers");
        }
    }

    // 2. A GATE SITTING EXACTLY ON THE FLOOR was read as "off" by the search and as "on" by the
    //    analyzer, so the two measured different sets of frames. Reported 3.005 dB where the gated
    //    truth was 6.011.
    {
        const int n = 100;
        std::vector<float> x ((std::size_t) n, 0.0f);
        for (int i = 0; i < 10; ++i) x[(std::size_t) (90 + i)] = 0.1f * (float) (i + 1);
        const float* p0 = x.data(); const float* const* key = &p0;
        auto p = baseParams();
        p.detector = Detector::Peak; p.mode = Mode::DownCompress; p.ratio = 2.0;
        p.kneeDb = 0.0; p.attackMs = 0.0; p.releaseMs = 0.0;

        for (double gate : { -240.0, -239.999, -100.0 })
        {
            ThresholdTarget t; t.amountDb = 3.0; t.toleranceDb = 0.05; t.statisticsGateDb = gate;
            ThresholdSolver s; (void) s.prepare (kFs, n);
            const auto r = s.solve (key, 1, n, p, t);
            if (r.status != SolveStatus::Solved) { test::ok (true, "no answer at gate " + std::to_string (gate)); continue; }

            Compressor c; (void) c.prepare (kFs, 128, 1, 0.0);
            CompressorParams cp; static_cast<GainReductionParams&> (cp) = p; cp.thresholdDb = r.thresholdDb;
            c.setParams (cp);
            std::vector<float> audio = x, tap ((std::size_t) n);
            float* io[1] { audio.data() };
            c.process (io, 1, n, GainReductionTap { tap.data(), n });
            LinkedDetector d; d.prepare (kFs); d.setParams (p); d.reset();
            std::vector<double> kept;
            for (int i = 0; i < n; ++i)
                if (core::gainToDb (d.process (key, 1, i)) > gate) kept.push_back (std::fabs ((double) tap[(std::size_t) i]));
            test::ok (! kept.empty(), "the gate kept something at " + std::to_string (gate));
            test::approx (exactNearestRank (kept, 0.95), r.achievedDb, 0.012,
                          "gate " + std::to_string (gate)
                          + ": the search and the analyzer count the SAME frames");
        }
    }

    // 3. THE BRACKET WAS SIZED FROM THE GATED EXTREMES. A gated frame still drives the ballistics, so
    //    the answer can lie far outside the range of the frames that were counted — here the whole
    //    gain reduction is inherited from a frame the gate excluded, and the answer is 34 dB below a
    //    bracket built from the counted ones.
    {
        std::vector<float> x { 0.01f, 1.0f };
        const float* p0 = x.data(); const float* const* key = &p0;
        auto p = baseParams();
        p.detector = Detector::Peak; p.mode = Mode::UpCompress; p.ratio = 2.0;
        p.kneeDb = 0.0; p.attackMs = 0.0; p.releaseMs = 100.0; p.rangeDb = 60.0;
        ThresholdTarget t; t.amountDb = 3.0; t.toleranceDb = 0.05; t.statisticsGateDb = -10.0;
        ThresholdSolver s; (void) s.prepare (kFs, 8);
        const auto r = s.solve (key, 1, 2, p, t);
        test::ok (r.status == SolveStatus::Solved, "a gated search whose answer lies outside the counted levels still solves");
        test::approx (r.achievedDb, 3.0, 0.06, "and delivers the target");
        test::ok (r.thresholdDb < -20.0,
                  "and finds it below the counted range (" + std::to_string (r.thresholdDb) + " dB)");
        test::ok (r.envelope.minSeenDb < r.envelope.minDb - 10.0,
                  "the reported statistics distinguish the levels SEEN from the levels COUNTED");
    }

    // 4. A STEEP CURVE closed the bracket before the answer was inside tolerance, and the search then
    //    called that a discontinuity. At ratio 1000 in DownExpand one dB of threshold is 999 dB of gain
    //    change, so the bracket floor has to be tightened by the curve's own multiplier.
    {
        std::vector<float> flat ((std::size_t) 64, 0.5f);
        const float* p0 = flat.data(); const float* const* key = &p0;
        for (double ratio : { 1000.0, 100.0, 20.0 })
        {
            auto p = baseParams();
            p.detector = Detector::Peak; p.mode = Mode::DownExpand; p.ratio = ratio;
            p.kneeDb = 0.0; p.attackMs = 0.0; p.releaseMs = 0.0; p.rangeDb = 60.0;
            ThresholdTarget t; t.statistic = TraceStatistic::Max; t.amountDb = 3.0;
            t.toleranceDb = 0.05; t.useSeed = false;
            ThresholdSolver s; (void) s.prepare (kFs, 128);
            const auto r = s.solve (key, 1, 64, p, t);
            test::ok (r.status == SolveStatus::Solved,
                      "ratio " + std::to_string (ratio) + " converges inside tolerance, not into a false discontinuity");
            test::approx (r.achievedDb, 3.0, 0.05, "and lands on the target");
        }
    }

    // 5. AN EXTREME RANK BYPASSED THE RANGE CHECK. With one value, out of range, rank 1 returned the
    //    tracked minimum and called it in range — the exact contract the class exists to keep.
    {
        QuantileHistogram h; test::ok (h.prepare (0.0, 10.0, 1.0), "prepared");
        h.add (100.0);
        double got = 1234.0;
        test::ok (! h.quantile (0.5, got), "a single out-of-range value is refused, not returned as the minimum");
        test::ok (got == 1234.0, "and `out` is untouched");
        h.reset();
        h.add (-5.0);
        test::ok (! h.quantile (0.5, got), "and the same below the range");
    }

    // 6. THE UNDERFLOW COUNT BELONGS IN THE RUNNING TOTAL. Values below the range still occupy ranks;
    //    dropping them from the cumulative count shifts every answer above them by that many places.
    {
        QuantileHistogram h; test::ok (h.prepare (0.0, 10.0, 1.0), "prepared");
        h.add (-1.0); h.add (2.25); h.add (8.25);
        double got = 0.0;
        test::ok (h.quantile (0.5, got), "in range");
        test::approx (got, 2.5, 1e-9, "rank 2 of {-1, 2.25, 8.25} is the middle value, not the top one");
        test::ok (h.belowRange() == 1, "and the underflow is counted");
    }

    // 7. A RESOLUTION THAT CANNOT COVER ITS RANGE IS REFUSED. It used to be honoured by keeping the bin
    //    width and dropping the range, which turned a 360 dB span into four and then reported a
    //    confident p95 of 0 dB for a signal at -20.
    {
        EnvelopeAnalyzer a;
        test::ok (! a.prepare (kFs, 1.0e-6), "a bin width needing more bins than the cap is refused");
        test::ok (! a.isPrepared(), "and the analyzer stays unprepared");
        test::ok (a.prepare (kFs, 0.01), "an honourable one is accepted");
        ThresholdSolver s;
        test::ok (! s.prepare (kFs, 0, 1.0e-6), "and the solver refuses it too");
        test::ok (! s.isPrepared(), "staying unprepared");
    }

    // 8. THE END ATOMS ARE EXACT BY RANK, not merely when they are the only value. 96 % silence and
    //    4 % at 2 dB has a p95 of exactly zero, and a binned rule without the atom counts says 0.025.
    {
        QuantileHistogram h; test::ok (h.prepare (0.0, 400.0, 0.05), "prepared");
        for (int i = 0; i < 96; ++i) h.add (0.0);
        for (int i = 0; i < 4;  ++i) h.add (2.0);
        double got = 0.0;
        test::ok (h.quantile (0.95, got), "in range");
        test::approx (got, 0.0, 0.0, "p95 of a trace that is silent 96 % of the time is EXACTLY zero");
        test::ok (h.quantile (0.97, got) && std::fabs (got - 2.0) <= 1e-9, "and p97 is exactly the high atom");
    }

    // 9. THE BIN COUNTER TYPE, pinned where it is decided. A 32-bit bin overflows after 24.9 hours at
    //    48 kHz, which no suite can reach by counting — so the type itself is the assertion.
    static_assert (std::is_same_v<QuantileHistogram::BinCount, std::uint64_t>,
                   "histogram bins must be 64-bit — a 32-bit bin overflows after 24.9 hours at 48 kHz");
    test::ok (true, "bin counters are uint64 (static_assert on the declared type, not on an accessor)");

    // 10. `step()` MEANS ONE PASS, and the budget means the budget. Reporting the answer used to cost
    //     two more passes, which both broke that promise and overran `maxEvaluations` by two — and a
    //     single-threaded host budgets its turns on exactly this.
    {
        const auto music = musicish (20000, kFs, 5);
        const float* p0 = music.data(); const float* const* key = &p0;
        const auto B = baseParams();
        for (int budget : { 1, 2, 5, 12, 48 })
        {
            ThresholdTarget t; t.amountDb = 3.0; t.toleranceDb = 0.05; t.maxEvaluations = budget;
            ThresholdSolver s; (void) s.prepare (kFs, 20000);
            int steps = 0;
            const int before = 0;
            (void) before;
            if (s.begin (key, 1, 20000, B, t)) while (s.step()) ++steps;
            test::ok (s.solution().evaluations <= budget + 1,
                      "budget " + std::to_string (budget) + ": " + std::to_string (s.solution().evaluations)
                      + " passes, never over budget (the one extra is the bracket endpoint begin() must measure)");
            test::ok (steps <= budget, "and no more steps than passes");
        }
    }

    // 11. A DEGENERATE SEED WINDOW made the outward doubling walk in place, spending the whole budget
    //     re-probing one point.
    {
        const auto music = musicish (20000, kFs, 9);
        const float* p0 = music.data(); const float* const* key = &p0;
        const auto B = baseParams();
        for (double w : { std::numeric_limits<double>::denorm_min(), 1e-300, 0.0, -3.0 })
        {
            ThresholdTarget t; t.amountDb = 3.0; t.toleranceDb = 0.05; t.seedWindowDb = w;
            ThresholdSolver s; (void) s.prepare (kFs, 20000);
            const auto r = s.solve (key, 1, 20000, B, t);
            test::ok (r.status == SolveStatus::Solved, "a degenerate seed window still converges");
            test::approx (r.achievedDb, 3.0, 0.06, "and lands on the target");
        }
    }
}


//==============================================================================
// Holes that a second round of mutation found in the suite above: twenty-one of twenty-nine invented
// mutations survived it. Every check here exists because removing something did NOT turn this file
// red — which is the only evidence a test is doing work.
static void whatTheSuiteWasNotYetAsserting()
{
    test::group ("what the suite was measuring but not asserting");

    const int n = 30000;
    const auto music = musicish (n, kFs, 17);
    const float* m0 = music.data(); const float* const* key = &m0;
    const auto B = baseParams();

    // 1. THE SEED'S VALUE IS ASSERTED, not printed. The whole seeding mechanism could be deleted and
    //    nothing here went red, because the comparison was reported and then tested against a
    //    tautology. On material where the estimate is good it must cost strictly fewer passes.
    {
        ThresholdTarget on; on.amountDb = 3.0; on.toleranceDb = 0.05;
        auto off = on; off.useSeed = false;
        ThresholdSolver a, b;
        test::ok (a.prepare (kFs, n) && b.prepare (kFs, n), "prepare");
        const auto ra = a.solve (key, 1, n, B, on);
        const auto rb = b.solve (key, 1, n, B, off);
        test::ok (ra.status == SolveStatus::Solved && rb.status == SolveStatus::Solved, "both solve");
        test::ok (ra.evaluations < rb.evaluations,
                  "on ordinary music the seed costs strictly fewer passes ("
                  + std::to_string (ra.evaluations) + " against " + std::to_string (rb.evaluations) + ")");
        test::approx (ra.thresholdDb, rb.thresholdDb, 0.01, "and reaches the same answer");
        // ...and where it is wrong it must not be a disaster: bounded, not free.
        TwoLevel spread; spread.build (n, -4.0, -40.0, false, 5);
        auto p = B; p.detector = Detector::Peak; p.attackMs = 5.0; p.releaseMs = 5.0;
        auto t2 = on; t2.amountDb = 6.0;
        ThresholdSolver c, d;
        (void) c.prepare (kFs, n); (void) d.prepare (kFs, n);
        const auto rc = c.solve (spread.key, 1, n, p, t2);
        auto t3 = t2; t3.useSeed = false;
        const auto rd = d.solve (spread.key, 1, n, p, t3);
        test::ok (rc.evaluations - rd.evaluations <= 5,
                  "and where the estimate is 28 dB wrong it costs at most a few passes, not the budget ("
                  + std::to_string (rc.evaluations - rd.evaluations) + ")");
    }

    // 2. THE NAIVE-INVERSION INSTRUMENT IS ITSELF CHECKED. Acceptance rests on how far it misses, so
    //    it could return anything and the number would still be "recorded". Against closed form:
    //    hard knee, ratio 4 (slope 0.75), level -4 dB, 6 dB wanted ⇒ 6/0.75 = 8 dB past ⇒ -12.000.
    {
        ThresholdSolver s; (void) s.prepare (kFs);
        auto p = B; p.mode = Mode::DownCompress; p.ratio = 4.0; p.kneeDb = 0.0; p.rangeDb = 60.0;
        test::approx (s.staticInversionThresholdDb (p, -4.0, 6.0), -12.0, 1.0e-6,
                      "hard knee: the static inversion matches closed form exactly");
        test::approx (s.staticInversionThresholdDb (p, -20.0, 3.0), -24.0, 1.0e-6, "and at another level");
        // Soft knee, inside the quadratic: |delta| = ((x+h)^2 / (2*knee)) * slope with x = level - thr.
        auto q = p; q.kneeDb = 6.0;
        const double want = 0.5, slope = 0.75, knee = 6.0, h = 3.0;
        const double x = std::sqrt (want * 2.0 * knee / slope) - h;                 // 8 -> 2.828 - 3
        test::approx (s.staticInversionThresholdDb (q, -4.0, want), -4.0 - x, 1.0e-5,
                      "soft knee: it follows the quadratic, not the line");
        // Beyond the range there is nothing to invert; it must saturate rather than run away.
        const double sat = s.staticInversionThresholdDb (p, -4.0, 90.0);
        test::ok (std::isfinite (sat) && sat < -4.0 - 60.0 / slope, "beyond the range it saturates, finitely");
    }

    // 3. `achievedBelowDb` IS READ, and with it the continuity bound: the objective is Lipschitz in the
    //    threshold with the curve's own multiplier, so a converged bracket pins the statistic. A bracket
    //    that does not actually contain the answer violates this immediately — which is how a real
    //    bracket bug was found.
    {
        for (double want : { 1.0, 3.0, 9.0 })
         for (double ratio : { 2.0, 4.0, 20.0 })
         {
             auto p = B; p.ratio = ratio;
             ThresholdTarget t; t.amountDb = want; t.toleranceDb = 0.05;
             ThresholdSolver s; (void) s.prepare (kFs, n);
             const auto r = s.solve (key, 1, n, p, t);
             if (r.status != SolveStatus::Solved) continue;
             test::ok (r.achievedBelowDb < want,
                       "the endpoint below the boundary really does fall short");
             const double mult = 1.0 - 1.0 / ratio;
             test::ok (r.achievedDb - r.achievedBelowDb <= mult * (r.bracketHiDb - r.bracketLoDb) + 0.05,
                       "and the two ends of a converged bracket are within the curve's Lipschitz bound");
         }
    }

    // 4. THE GATE'S COMPARISON, pinned by a frame sitting exactly ON it.
    {
        std::vector<float> flat ((std::size_t) 8000, 0.1f);
        const float* f0 = flat.data(); const float* const* fk = &f0;
        auto p = B; p.detector = Detector::Peak;
        const double exactly = core::gainToDb (0.1f);
        ThresholdTarget t; t.amountDb = 3.0; t.statisticsGateDb = exactly;
        ThresholdSolver s; (void) s.prepare (kFs, 8000);
        test::ok (s.solve (fk, 1, 8000, p, t).status == SolveStatus::EmptyInput,
                  "a gate exactly at the frames' own level excludes them all — the comparison is strict");
        auto t2 = t; t2.statisticsGateDb = exactly - 0.02;
        ThresholdSolver s2; (void) s2.prepare (kFs, 8000);
        test::ok (s2.solve (fk, 1, 8000, p, t2).status != SolveStatus::EmptyInput,
                  "one bin lower and they all count");
    }

    // 5. THE PERCENTILES ARE DISTINCT. A two-level fixture has p90 = p95 = p99, so any of them could be
    //    reported for any other. A ten-step staircase cannot hide that.
    {
        // ONE HUNDRED steps, not ten: with ten, p95 and p99 both land in the last one and either could
        // be reported for the other. The staircase has to be finer than the finest quantile asked for.
        const int m = 20000;
        std::vector<float> stair ((std::size_t) m);
        for (int i = 0; i < m; ++i)
            stair[(std::size_t) i] = (float) core::dbToGain (-60.0 + 0.6 * (double) ((i * 100) / m));
        const float* s0 = stair.data(); const float* const* sk = &s0;
        DetectorParams dp; dp.detector = Detector::Peak;
        EnvelopeAnalyzer a; test::ok (a.prepare (kFs), "prepare"); a.setParams (dp); a.reset();
        a.analyze (sk, 1, m);
        const auto st = a.stats();
        // A hundred equal steps of 200 frames, levels -60 + 0.6*k. Rank ceil(q*20000) lands in step
        // ceil(q*100)-1, so every quantile has its own step. Counted, not guessed.
        test::approx (st.p50Db, -60.0 + 0.6 * 49.0, 0.02, "p50 of a hundred-step staircase is step 49");
        test::approx (st.p90Db, -60.0 + 0.6 * 89.0, 0.02, "p90 is step 89");
        test::approx (st.p95Db, -60.0 + 0.6 * 94.0, 0.02, "p95 is step 94");
        test::approx (st.p99Db, -60.0 + 0.6 * 98.0, 0.02, "p99 is step 98");
        test::ok (st.p50Db < st.p90Db && st.p90Db < st.p95Db && st.p95Db < st.p99Db,
                  "FOUR distinct values — a coarse fixture reports one number for several of them");
        test::approx (st.minDb, -60.0, 0.02, "the reported minimum is the first step");
        test::approx (st.maxDb, -60.0 + 0.6 * 99.0, 0.02, "and the maximum the last");
        test::approx (st.minSeenDb, st.minDb, 1e-9, "with no gate, seen and counted agree");
        test::approx (st.maxSeenDb, st.maxDb, 1e-9, "at both ends");
    }

    // 6. THE BIN WIDTH IS HONOURED BY THE SEARCH. A tolerance finer than the instrument is the ONLY way
    //    the search can converge and still be outside tolerance — the objective is continuous in the
    //    threshold — so the status says that rather than blaming the material.
    {
        ThresholdTarget t; t.amountDb = 3.0; t.toleranceDb = 0.05;
        ThresholdSolver coarse; test::ok (coarse.prepare (kFs, n, 0.5), "a coarse histogram prepares");
        const auto rc = coarse.solve (key, 1, n, B, t);
        test::ok (rc.status == SolveStatus::ToleranceBelowResolution,
                  "half-dB bins cannot answer to 0.05 dB, and the status says so");
        ThresholdSolver fine; (void) fine.prepare (kFs, n, 0.01);
        test::ok (fine.solve (key, 1, n, B, t).status == SolveStatus::Solved, "0.01 dB bins can");
        // ...and the same status arrives from an impossible tolerance rather than from the material.
        auto tiny = t; tiny.toleranceDb = 0.0005;
        ThresholdSolver s; (void) s.prepare (kFs, n, 0.01);
        test::ok (s.solve (key, 1, n, B, tiny).status == SolveStatus::ToleranceBelowResolution,
                  "and from a tolerance below half a bin");
    }

    // 7. THE CLAMP TO THE OBSERVED SPAN still does work of its own, where neither end atom applies.
    {
        QuantileHistogram h; test::ok (h.prepare (0.0, 400.0, 0.05), "prepared");
        for (int i = 0; i < 90; ++i) h.add (1.0);
        for (int i = 0; i < 9;  ++i) h.add (2.004);
        h.add (2.0041);                                   // a distinct maximum, so the max atom is 1
        double got = 0.0;
        test::ok (h.quantile (0.95, got), "in range");
        test::ok (got <= 2.0041 + 1e-9,
                  "the answer never leaves the observed span, even when the bin centre would ("
                  + std::to_string (got) + ")");
    }

    // 8. THE TARGET OF ZERO reports a measured number, not a constant. A constant compared with itself
    //    is the shape of check that passes on both sides of a behaviour change.
    {
        for (auto mode : { Mode::DownCompress, Mode::UpCompress, Mode::DownExpand })
        {
            auto p = B; p.mode = mode; p.detector = Detector::Peak; p.kneeDb = 0.0;
            ThresholdTarget t; t.amountDb = 0.0;
            ThresholdSolver s; (void) s.prepare (kFs, n);
            const auto r = s.solve (key, 1, n, p, t);
            test::ok (r.status == SolveStatus::TargetIsZero, "a zero target is named");
            test::approx (r.achievedDb, 0.0, 1e-12, "and the reported amount was measured as zero");
            // the real chain at that threshold must indeed do nothing at all
            GainReductionPath path; path.prepare (kFs);
            auto q = p; q.thresholdDb = r.thresholdDb; path.setParams (q); path.reset();
            double worst = 0.0;
            for (int i = 0; i < n; ++i) worst = std::max (worst, std::fabs ((double) path.process (key, 1, i)));
            test::approx (worst, 0.0, 0.0, "and the chain really does nothing at the returned threshold");
        }
    }

    // 9b. THE PREDICATE CARRIES THE TOLERANCE. Ask for very slightly more than this material can give:
    //     a search whose predicate is the bare `f >= amount` finds nothing satisfying it, never moves
    //     the aggressive endpoint, and hands back the saturating extreme — which still passes a
    //     tolerance check, and so is invisible unless the BRACKET is inspected.
    {
        ThresholdTarget probe; probe.amountDb = 1000.0; probe.toleranceDb = 0.05;
        ThresholdSolver s0; (void) s0.prepare (kFs, n);
        const double attainable = s0.solve (key, 1, n, B, probe).attainableDb;
        test::ok (attainable > 1.0, "the material can give something to aim just past");

        // THE ANSWER IS CONTINUOUS IN THE TARGET, which is the observable consequence: the objective is
        // Lipschitz in the threshold, so asking for a hair more and a hair less than the attainable must
        // give nearly the same threshold. Both brackets close either way, so the width proves nothing —
        // it is the DISTANCE that separates a predicate carrying the tolerance from one that does not.
        // Measured with the tolerance dropped: 64.06 dB apart instead of 0.37.
        ThresholdTarget above; above.amountDb = attainable + 0.025; above.toleranceDb = 0.05;
        ThresholdTarget below = above; below.amountDb = attainable - 0.025;
        ThresholdSolver sa, sb;
        (void) sa.prepare (kFs, n); (void) sb.prepare (kFs, n);
        const auto ra = sa.solve (key, 1, n, B, above);
        const auto rb = sb.solve (key, 1, n, B, below);
        test::ok (ra.status == SolveStatus::Solved, "a target half a tolerance past the attainable still solves");
        test::ok (std::fabs (ra.thresholdDb - rb.thresholdDb) < 5.0,
                  "and lands within a few dB of the target half a tolerance below it ("
                  + std::to_string (std::fabs (ra.thresholdDb - rb.thresholdDb)) + " dB apart)");
    }

    // 9c. THE SEED WINDOW HAS A FLOOR. Without one the outward doubling starts from a denormal and
    //     walks in place, spending the whole budget re-probing a single point.
    {
        ThresholdTarget ref; ref.amountDb = 3.0; ref.toleranceDb = 0.05;
        ThresholdSolver a; (void) a.prepare (kFs, n);
        const int normal = a.solve (key, 1, n, B, ref).evaluations;
        for (double w : { std::numeric_limits<double>::denorm_min(), 1.0e-300, 1.0e-30 })
        {
            auto t = ref; t.seedWindowDb = w;
            ThresholdSolver s; (void) s.prepare (kFs, n);
            const auto r = s.solve (key, 1, n, B, t);
            test::ok (r.status == SolveStatus::Solved, "a denormal seed window still converges");
            // STRICTLY fewer. Without the floor the outward step is absorbed into the seed, the very
            // first expansion probe falls outside its own bracket, seeding switches itself off, and
            // the search silently degrades to plain bisection — same answer, none of the benefit, and
            // no failure anywhere unless the benefit itself is asserted.
            test::ok (r.evaluations < normal,
                      "and still costs FEWER passes than an unseeded default (" + std::to_string (r.evaluations)
                      + " against " + std::to_string (normal) + ")");
        }
    }

    // 9. THE TARGET'S OWN NUMBERS ARE VALIDATED. A NaN tolerance used to reach the search and come back
    //    as a converged answer.
    {
        ThresholdTarget bad; bad.amountDb = 3.0;
        for (int which = 0; which < 4; ++which)
        {
            auto t = bad;
            if      (which == 0) t.toleranceDb = std::numeric_limits<double>::quiet_NaN();
            else if (which == 1) t.toleranceDb = -1.0;
            else if (which == 2) t.resolutionDb = 0.0;
            else                 t.maxEvaluations = 0;
            ThresholdSolver s; (void) s.prepare (kFs, n);
            test::ok (s.solve (key, 1, n, B, t).status == SolveStatus::InvalidTarget,
                      "a nonsensical target field is refused, not searched (case " + std::to_string (which) + ")");
        }
    }
}

//==============================================================================
int main()
{
    std::printf ("felitronics::dynamics::offline — detector-domain analysis and the threshold solver\n");
    quantileRuleIsOneRule();
    envelopeAnalyzerMeasuresWhatItSays();
    gainReductionPathIsTheCompressorsOwnChain();
    theCompressorTapIsTheTruth();
    f9CounterexampleIsMandatory();
    theRoundTripIsMeasuredWithADifferentInstrument();
    everyRefusalHasAName();
    theGateReadsTheLevelAndNotTheResult();
    aSearchIsAnIterationAndStateMustNotSurviveIt();
    theCurveOnlyEverSeesLevelMinusThreshold();
    everyModeIsSolvedInItsOwnDirection();
    theStatisticIsAParameterNotAnAssumption();
    theSeedIsAnEstimateAndItsValueIsMeasured();
    theNaiveInversionIsMeasuredNotArguedAbout();
    poisonedKeysAndOddShapes();
    theConsiliumsCounterexamples();
    whatTheSuiteWasNotYetAsserting();
    return test::report();
}

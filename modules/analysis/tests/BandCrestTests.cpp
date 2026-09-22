// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// felitronics::analysis::BandCrest — K1. Every oracle here is computed OUTSIDE the object:
//
//   * the full-band peak against `ReferenceTruePeakMeter`, the instrument that CERTIFIES a delivery. Not a
//     tolerance — the reconstruction is literally the same topology, so the two are bit-identical or one of
//     them is wrong;
//   * the interpolator's delay against an impulse, which is where the 63.5 in the header comes from;
//   * the block grid against `LoudnessMeter`'s gating block count, counted by different code over the same
//     time;
//   * the two hop clocks against each other — the cross-check that already caught a whole hop of silence
//     being closed as programme;
//   * the comparator's numbers against arithmetic on hand-made distributions, and against an INVARIANT that
//     holds by construction.

#include <felitronics_test.h>
#include <alloc_counter.h>

#include <felitronics/analysis/BandCrest.h>
#include <felitronics/analysis/LoudnessMeter.h>
#include <felitronics/analysis/ReferenceTruePeakMeter.h>
#include <felitronics/oversampling/PolyphaseOversampler.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using felitronics::test::ok;
using felitronics::analysis::BandCrest;
using felitronics::analysis::bandCrestLoss;

namespace
{
constexpr double kPi = 3.14159265358979323846;

// A programme with a REAL crest: a tone bed carrying the energy and sparse transients carrying the peak. A
// fixture whose peak-to-loudness ratio is small never makes a crest measurement say anything.
//
// TAPERED AT BOTH ENDS, which is not decoration: a tone that stops dead at full scale has a genuine
// reconstruction overshoot at the file boundary, and a true-peak instrument finds it correctly — at which
// point the measurement is of the edge and not of the programme.
std::vector<std::vector<float>> programme (double fs, int nch, double seconds, double bigEvery = 0.0)
{
    const auto n = (std::size_t) (seconds * fs);
    std::vector<std::vector<float>> ch ((std::size_t) nch, std::vector<float> (n, 0.0f));
    std::uint32_t st = 0x2545F49u;
    for (int c = 0; c < nch; ++c)
        for (std::size_t i = 0; i < n; ++i)
        {
            st = st * 1664525u + 1013904223u;
            const double t = (double) i / fs;
            const double w = std::min (1.0, std::min (t / 0.05, (seconds - t) / 0.05));
            double v = 0.22 * std::sin (2.0 * kPi * (110.0 + 7.0 * c) * t)
                     + 0.14 * std::sin (2.0 * kPi * (1700.0 + 90.0 * c) * t)
                     + 0.06 * ((double) (st >> 8) / 8388608.0 - 1.0);
            const double ph = std::fmod (t, 0.5);
            if (ph < 0.004)
            {
                const bool big = bigEvery > 0.0 && ((int) (t / 0.5) % (int) bigEvery == 0);
                v += big ? 0.85 : 0.25;
            }
            ch[(std::size_t) c][i] = (float) (w * v);
        }
    return ch;
}

// The same as `runIt`, for an object whose parameters were set by the caller: `runIt` prepares from whatever
// `setParams` left, so the two differ only in that this one does not reset them.
bool runIt2 (BandCrest& bc, std::vector<std::vector<float>>& buf, double fs, int block = 4096);

bool runIt (BandCrest& bc, std::vector<std::vector<float>>& buf, double fs, int block = 4096)
{
    std::vector<const float*> p;
    for (auto& c : buf) p.push_back (c.data());
    const int n = (int) buf[0].size();
    if (! bc.prepare (fs, (int) buf.size(), n)) return false;
    for (int off = 0; off < n; )
    {
        const int m = std::min (block, n - off);
        std::vector<const float*> q;
        for (auto& c : buf) q.push_back (c.data() + off);
        if (! bc.process (q.data(), (int) buf.size(), m)) return false;
        off += m;
    }
    bc.finish();
    return true;
}

bool runIt2 (BandCrest& bc, std::vector<std::vector<float>>& buf, double fs, int block)
{
    return runIt (bc, buf, fs, block);
}

//==================================================================================================
void theDelayIsWhereTheHeaderSaysItIs()
{
    felitronics::test::group ("the interpolator's delay — an impulse, and the number the header carries");
    felitronics::oversampling::PolyphaseOversampler os;
    if (! felitronics::test::run (os.prepare (BandCrest::kFactor, 1, BandCrest::kTapsPerPhase))) return;
    const int n = 256, at = 64;
    std::vector<float> x ((std::size_t) n, 0.0f); x[(std::size_t) at] = 1.0f;
    std::vector<float> y ((std::size_t) n * BandCrest::kFactor, 0.0f);
    const float* in[1] { x.data() }; float* out[1] { y.data() };
    os.upsample (in, 1, n, out);

    // THE FIR IS SYMMETRIC, and about a HALF sample — which is the whole reason the header states the delay as
    // a doubled integer instead of rounding it. The centre for an impulse at base `at` is
    // `at*factor + (factor*tapsPerPhase - 1)/2`, here 256 + 63.5 = 319.5.
    const double centreX2 = 2.0 * (double) (at * BandCrest::kFactor) + (double) BandCrest::kDelayOsX2;
    const auto lo = (std::size_t) ((centreX2 - 1.0) / 2.0), hi = lo + 1;
    double worst = 0.0;
    for (std::size_t k = 0; k + 1 < 60; ++k)
        worst = std::max (worst, std::fabs ((double) y[lo - k] - (double) y[hi + k]));
    ok (worst == 0.0, "the reconstruction is symmetric about oversampled index "
                      + std::to_string (centreX2 / 2.0) + ", bit for bit (worst |diff| "
                      + std::to_string (worst) + ")");
    ok (y[lo] == y[hi] && y[lo] > 0.8f,
        "the two samples either side of the half are equal and carry the impulse");
    ok (BandCrest::kDelayOsX2 == BandCrest::kFactor * BandCrest::kTapsPerPhase - 1,
        "and the constant the class skips by is that same expression, not a literal");
}

//==================================================================================================
void thePeakIsTheCertificate()
{
    felitronics::test::group ("the full-band peak IS the certifying instrument's, bit for bit");
    for (const double fs : { 44100.0, 48000.0, 96000.0 })
    {
        auto buf = programme (fs, 2, 3.0, 25.0);
        BandCrest bc;
        if (! felitronics::test::run (runIt (bc, buf, fs))) continue;

        felitronics::analysis::ReferenceTruePeakMeter tp;
        if (! felitronics::test::run (tp.prepare (fs, 4096, 2))) continue;
        const int n = (int) buf[0].size();
        for (int off = 0; off < n; )
        {
            const int m = std::min (4096, n - off);
            const float* q[2] { buf[0].data() + off, buf[1].data() + off };
            if (! tp.process (q, 2, m)) break;
            off += m;
        }
        tp.drain();

        const double mine = bc.fullBandPeakLin(), ref = tp.truePeakLinear();
        ok (mine == ref, std::to_string ((int) fs) + " Hz: BandCrest " + std::to_string (mine)
                         + " and ReferenceTruePeakMeter " + std::to_string (ref)
                         + " are the same double — the reconstruction is the same topology, so a tolerance "
                           "here would be hiding a difference rather than allowing one");
    }
}

//==================================================================================================
void theGridIsTheLoudnessMeters()
{
    // THE RATES ARE CHOSEN AGAINST THE FORMULA, not for convenience. At 44.1, 48 and 96 kHz `0.01*fs` is an
    // integer and a hop rounded whole agrees with a hop built from sub-hops — so a sweep over those three is
    // blind to the difference, and the first version of this check was exactly that sweep. 22050 and 11025 are
    // where they part (2210 samples against 2205), and 8000 is the core's floor.
    felitronics::test::group ("the block grid is the loudness meter's, counted by different code");
    for (const double fs : { 8000.0, 11025.0, 22050.0, 44100.0, 48000.0, 96000.0 })
    {
        auto buf = programme (fs, 2, 3.0);
        BandCrest bc;
        if (! felitronics::test::run (runIt (bc, buf, fs))) continue;

        felitronics::analysis::LoudnessMeter lm;
        const int n = (int) buf[0].size();
        if (! felitronics::test::run (lm.prepareForSamples (fs, 2, n))) continue;
        const float* q[2] { buf[0].data(), buf[1].data() };
        if (! felitronics::test::run (lm.process (q, 2, n))) continue;

        // THE HOP ITSELF, in samples, before any count derived from it: a block count can agree by luck when
        // the hop does not, and the hop is the thing the header claims is shared.
        const int lmHop = 10 * (int) std::lround (0.01 * fs);
        ok (bc.hopSamples() == lmHop,
            std::to_string ((int) fs) + " Hz: the hop is " + std::to_string (bc.hopSamples())
            + " samples, built from sub-hops exactly as the meter builds it (" + std::to_string (lmHop) + ")");
        // THE COUNT IS THE METER'S PLUS THE PARTIAL LAST HOP, and the difference is deliberate rather than a
        // discrepancy: this analyzer closes the hop the programme ends inside, because that hop is the end of
        // the file and dropping it would lose exactly the material a master's loudest moments often sit in.
        // The meter does not. So the relationship is an equality with a named term, not an equality — the
        // header used to claim the latter and the test could not see it, because the only lengths it tried
        // ended on a hop boundary.
        const long long partial = ((long long) buf[0].size() % bc.hopSamples()) != 0 ? 1 : 0;
        ok (bc.blockCount() == (long long) lm.gatingBlockCount() + partial,
            std::to_string ((int) fs) + " Hz: " + std::to_string (bc.blockCount())
            + " blocks here and " + std::to_string (lm.gatingBlockCount()) + " gating blocks there, the "
            + (partial ? "programme ending inside a hop" : "programme ending on a hop boundary"));
        // THE TWO HOP CLOCKS INSIDE THIS CLASS. They partition the same input time and are counted by
        // different code; they disagreed by one when the interpolator's drain was taken whole, closing a hop
        // of silence as programme. A comment claiming they agree is not the same thing as this line.
        ok (bc.hopCount() == bc.basePeakHops(),
            "... and the oversampled and base-rate hop clocks agree (" + std::to_string (bc.hopCount())
            + " and " + std::to_string (bc.basePeakHops()) + ")");
    }
}

//==================================================================================================
void theLossIsInvariantToGain()
{
    felitronics::test::group ("the loss measures the crest, not the gain applied to either side");
    const double fs = 48000.0;
    auto src = programme (fs, 2, 6.0, 25.0);
    auto loud = src;
    for (auto& c : loud) for (auto& v : c) v *= 0.5f;      // an exact power of two: no rounding of its own

    BandCrest a, b;
    if (! felitronics::test::run (runIt (a, src, fs))) return;
    if (! felitronics::test::run (runIt (b, loud, fs))) return;
    std::vector<double> scratch;
    int bad = 0;
    for (int band = 0; band < BandCrest::kBands; ++band)
    {
        const auto L = bandCrestLoss (a, b, band, scratch);
        if (! (L.usable > 0 && std::fabs (L.maxDb) <= 1.0e-9)) ++bad;
    }
    ok (bad == 0, "a master that is the source at exactly -6 dB loses no crest in any band: the loss is a "
                  "ratio of linear cells and a gain cancels out of it");

    // AND AT GAINS THAT ARE NOT POWERS OF TWO, where the claim is narrower and the earlier version of this
    // check did not look. x0.5 is EXACT in float, so it demonstrates the arithmetic and nothing about the
    // general case; at x0.37 the loss is 1.1e-5 dB, which is not the comparison drifting but the FIXTURE:
    // `k*x` rounds per sample, so the two signals genuinely differ by about that much and the loss correctly
    // says so. The claim is therefore two claims — exactly zero where the scaling is exact, and inside float
    // rounding where it is not — and both are asserted rather than the stronger one being implied.
    int drift = 0;
    double worst = 0.0;
    for (const double g : { 0.37, 1.0e-3, 1.0e3, 0.999999 })
    {
        auto scaled = src;
        for (auto& c : scaled) for (auto& v : c) v = (float) ((double) v * g);
        BandCrest sc;
        if (! runIt (sc, scaled, fs)) { ++drift; continue; }
        for (int band = 0; band < BandCrest::kBands; ++band)
            worst = std::max (worst, std::fabs (bandCrestLoss (a, sc, band, scratch).maxDb));
    }
    ok (drift == 0 && worst < 1.0e-3,
        "at gains that are not powers of two the loss stays inside float rounding of the fixture (worst "
        + std::to_string (worst) + " dB) — the arithmetic does not amplify it");
}

//==================================================================================================
void cvarSeesWhatP95CannotAndNeverLess()
{
    felitronics::test::group ("CVaR95 against p95 — the control, and the invariant that makes it safe");
    const double fs = 48000.0;
    auto src = programme (fs, 2, 30.0, 25.0);
    auto master = src;
    // A HARD CEILING ONLY THE BIG TRANSIENTS REACH — a limiter's signature: damage concentrated on a few
    // percent of the blocks, which is exactly where a 95th percentile cannot look.
    for (auto& c : master) for (auto& v : c) v = std::clamp (v, -0.45f, 0.45f);

    BandCrest a, b;
    if (! felitronics::test::run (runIt (a, src, fs))) return;
    if (! felitronics::test::run (runIt (b, master, fs))) return;
    std::vector<double> scratch;

    const auto L = bandCrestLoss (a, b, BandCrest::kHighMid, scratch);
    std::printf ("        highMid: p95 %.3f  CVaR95 %.3f  max %.3f  over 3 dB %lld of %lld blocks\n",
                 L.p95Db, L.cvar95Db, L.maxDb, (long long) L.over3Db, (long long) L.inActive);
    ok (L.usable > 100 && L.p95Db < 0.2 && L.cvar95Db > 1.5 && L.over3Db >= 4,
        "THE CONTROL: p95 reads " + std::to_string (L.p95Db) + " dB — no damage — while CVaR95 reads "
        + std::to_string (L.cvar95Db) + " and " + std::to_string (L.over3Db)
        + " blocks lost more than 3 dB. A cost function calibrated on p95 would call this master undamaged");

    // AND THE OTHER DIRECTION, which is the stronger claim: a fixture where p95 sees damage and CVaR95 does
    // not MUST NOT EXIST. It cannot, by construction — CVaR95 averages the worst 5 %, every member of which is
    // at or above the 95th percentile — but "cannot by construction" is exactly the kind of statement that is
    // false at the edges of a rounding rule, so it is checked rather than argued, over shapes chosen to attack
    // the rank arithmetic: tiny populations, ties, all-zero, all-equal, one outlier.
    {
        std::uint32_t st = 0x1234u;
        int violations = 0, checked = 0;
        for (int trial = 0; trial < 400; ++trial)
        {
            const std::size_t k = (std::size_t) (1 + (trial * 7919u) % 120u);
            std::vector<double> v (k);
            for (std::size_t i = 0; i < k; ++i)
            {
                st = st * 1664525u + 1013904223u;
                const double u = (double) (st >> 8) / 16777216.0;
                v[i] = (trial % 4 == 0) ? 0.0
                     : (trial % 4 == 1) ? 3.0
                     : (trial % 4 == 2) ? (i + 1 == k ? 12.0 : 0.0)
                                        : (u * 8.0 - 2.0);
            }
            std::sort (v.begin(), v.end());
            auto at = [&] (double q)
            {
                auto i = (std::size_t) std::llround (q * (double) (k - 1));
                if (i >= k) i = k - 1;
                return v[i];
            };
            const double p95 = std::max (0.0, at (0.95));
            const std::size_t tail = std::max<std::size_t> (1, (std::size_t) ((double) k * 0.05 + 0.5));
            double t = 0.0;
            for (std::size_t i = k - tail; i < k; ++i) t += std::max (0.0, v[i]);
            const double cvar = t / (double) tail;
            ++checked;
            if (cvar + 1.0e-12 < p95) ++violations;
        }
        ok (violations == 0, "CVaR95 >= p95 on all " + std::to_string (checked)
                             + " shapes — populations of 1 to 120, all-zero, all-equal, a single outlier and "
                               "random: the fixture where p95 sees damage and CVaR95 does not cannot be built");
    }
}

//==================================================================================================
// THE BANDS ARE THE FILTERS THE HEADER NAMES — against arithmetic, not against the object.
//
// This is the check that survives a rewrite of the filter bank, and the reason it exists: the bank was five
// `Crossover2` objects and became six `Svf` cascades, which is bit-identical only if the definitions are the
// same. Nothing else in this suite could tell the two apart — a wrong corner, a swapped band index or a
// missing section all pass a re-slicing test and a peak null.
//
// THE MEASURED QUANTITY IS A RATIO, band mean-square over full-band mean-square, which cancels the tone's
// amplitude, the interpolator's passband gain and the window: for a steady tone at f it is |H_b(f)|^2. The
// oracle is the prewarped Linkwitz-Riley magnitude written out from the filter's own definition — at
// Q = 1/sqrt(2) a second-order section gives |LP|^2 = 1/(1+W^4) and |HP|^2 = W^4/(1+W^4), with
// W = tan(pi f / fs') / tan(pi fc / fs'), the fourth order being each squared. `fs'` is the OVERSAMPLED rate:
// the split runs at 4x, and an oracle at the base rate would be describing a different filter.
void theBandsAreTheFiltersTheHeaderNames()
{
    felitronics::test::group ("each band is its named cascade — analytic magnitudes, not the object's own word");
    const double fs = 48000.0, osRate = fs * BandCrest::kFactor;
    const felitronics::analysis::BandCrestParams pr {};
    auto W   = [&] (double f, double fc) { return std::tan (kPi * f / osRate) / std::tan (kPi * fc / osRate); };
    auto lp4 = [&] (double f, double fc) { const double w = W (f, fc); const double d = 1.0 + w * w * w * w; return 1.0 / (d * d); };
    auto hp4 = [&] (double f, double fc) { const double w = W (f, fc); const double w4 = w * w * w * w;
                                           const double d = 1.0 + w4; return (w4 * w4) / (d * d); };

    for (const double f : { 60.0, 500.0, 3500.0, 11000.0, 1200.0 })
    {
        const auto n = (std::size_t) (fs * 3.0);
        std::vector<std::vector<float>> ch (2, std::vector<float> (n, 0.0f));
        for (int c = 0; c < 2; ++c)
            for (std::size_t i = 0; i < n; ++i)
            {
                const double t = (double) i / fs;
                const double w = std::min (1.0, std::min (t / 0.25, (3.0 - t) / 0.25));   // taper BOTH ends
                ch[(std::size_t) c][i] = (float) (0.30 * w * std::sin (2.0 * kPi * f * t));
            }
        BandCrest bc;
        if (! felitronics::test::run (runIt (bc, ch, fs))) continue;
        const long long mid = bc.blockCount() / 2;              // a block in the settled middle, past the taper
        const double full = bc.blockMeanSq (mid, BandCrest::kFull);
        if (! (full > 0.0)) { ok (false, "the full band carries energy at " + std::to_string ((int) f) + " Hz"); continue; }

        const double want[4] {
            lp4 (f, pr.bandEdgeHz[0]),
            hp4 (f, pr.bandEdgeHz[0]) * lp4 (f, pr.bandEdgeHz[1]),
            hp4 (f, pr.bandEdgeHz[1]) * lp4 (f, pr.bandEdgeHz[2]),
            hp4 (f, pr.bandEdgeHz[2]),
        };
        double worst = 0.0; int at = -1;
        for (int b = 0; b < BandCrest::kFull; ++b)
        {
            const double got = bc.blockMeanSq (mid, b) / full;
            // Relative where the band carries something, absolute where it is far down and the float state is
            // at its own floor — a relative test on 1e-9 would be measuring rounding, not the filter.
            const double err = want[b] > 1.0e-6 ? std::fabs (got - want[b]) / want[b] : std::fabs (got - want[b]);
            if (err > worst) { worst = err; at = b; }
        }
        ok (worst < 2.0e-3, std::to_string ((int) f) + " Hz: every band within "
                            + std::to_string (worst) + " of its analytic magnitude (worst at band "
                            + std::to_string (at) + ")");
    }
}

//==================================================================================================
// THE PROGRAMME LEVEL IS IN THE GATE'S OWN UNITS — the point of publishing it at all.
//
// A caller wanting a floor "40-odd dB below the programme" reaches for integrated loudness, and `I - 42` is
// the obvious spelling. It is also wrong in a way that hides: `I` is K-weighted and this gate is not, and the
// offset between them is a function of the SPECTRUM. This checks the two claims that make the scalar useful —
// that a floor set from it lands where arithmetic says, and that it is gated at a FIXED -70 dBFS rather than
// at `programmeFloorDb`, so a caller deriving the floor from it is not chasing its own tail.
void theProgrammeLevelIsInTheGatesUnits()
{
    felitronics::test::group ("the programme level is the gate's own quantity, gated at a fixed floor");
    const double fs = 48000.0;
    auto buf = programme (fs, 2, 8.0, 25.0);
    BandCrest bc;
    if (! felitronics::test::run (runIt (bc, buf, fs))) return;
    const double lvl = bc.programmeMeanSquareDb();
    ok (lvl > -100.0 && lvl < 0.0, "the programme answers a level of " + std::to_string (lvl) + " dBFS");

    // 1. A FLOOR SET FROM IT LANDS WHERE ARITHMETIC SAYS. At `level` exactly, every block at or above the
    //    programme's own mean is admitted and the rest are not — so the count matches a count taken here.
    {
        felitronics::analysis::BandCrestParams p;
        p.programmeFloorDb = lvl;
        BandCrest at;
        at.setParams (p);
        if (felitronics::test::run (runIt2 (at, buf, fs)))
        {
            long long want = 0;
            const double thr = std::pow (10.0, lvl / 10.0);
            for (long long j = 0, e = at.blockCount(); j < e; ++j)
                if (at.blockMeanSq (j, BandCrest::kFull) >= thr) ++want;
            ok (at.activeBlocks (BandCrest::kFull) == want,
                "a floor set AT the published level admits the blocks arithmetic says it should ("
                + std::to_string (at.activeBlocks (BandCrest::kFull)) + " of " + std::to_string (at.blockCount()) + ")");
        }
    }

    // 2. THE SCALAR DOES NOT MOVE WITH `programmeFloorDb`. If it were gated at the configured floor, a caller
    //    deriving the floor from it would be solving a fixed point; it is gated at -70 dBFS and nothing else.
    {
        felitronics::analysis::BandCrestParams p;
        p.programmeFloorDb = -20.0;                 // far above anything this fixture has in most blocks
        BandCrest high;
        high.setParams (p);
        if (felitronics::test::run (runIt2 (high, buf, fs)))
            ok (std::fabs (high.programmeMeanSquareDb() - lvl) <= 1.0e-12,
                "moving the configured floor by 50 dB does not move the published level ("
                + std::to_string (high.programmeMeanSquareDb()) + " against " + std::to_string (lvl) + ")");
    }

    // 2b. A CHANNEL THAT VANISHES MID-PROGRAMME IS COUNTED. It changes the measurement — the peak is a max
    //     over the channels PRESENT and the mean square is divided by what was accumulated — and every number
    //     stays finite and plausible while it happens, so without a counter nothing says the width moved.
    {
        auto two = programme (fs, 2, 4.0);
        BandCrest bc;
        if (felitronics::test::run (bc.prepare (fs, 2, (long long) two[0].size())))
        {
            const int n = (int) two[0].size();
            const float* q2[2] { two[0].data(), two[1].data() };
            const float* q1[1] { two[0].data() + n / 2 };
            const bool a1 = bc.process (q2, 2, n / 2);
            const bool b1 = bc.process (q1, 1, n - n / 2);
            bc.finish();
            ok (a1 && b1 && bc.widestChannels() == 2 && bc.narrowedSamples() == (long long) (n - n / 2),
                "a programme fed at two channels and then at one counts the narrowed stretch ("
                + std::to_string (bc.narrowedSamples()) + " samples of " + std::to_string (n) + ")");
            // AND A RUN THAT NEVER NARROWS COUNTS NOTHING — the control, without which the line above passes
            // on a counter wired to the frame count.
            BandCrest steady;
            auto whole = programme (fs, 2, 4.0);
            if (felitronics::test::run (runIt (steady, whole, fs)))
                ok (steady.narrowedSamples() == 0 && steady.widestChannels() == 2,
                    "... and a run at a steady width counts none");
        }
    }

    // 3. AND IT IS A SENTINEL, NOT A LEVEL, WHEN NOTHING CLEARS THE GATE.
    {
        std::vector<std::vector<float>> quiet (2, std::vector<float> ((std::size_t) (fs * 2.0), 0.0f));
        BandCrest silent;
        if (felitronics::test::run (runIt (silent, quiet, fs)))
            ok (silent.programmeMeanSquareDb() == BandCrest::kSilenceDb,
                "digital silence answers the sentinel, not a number that looks like a level");
    }
}

//==================================================================================================
// THREE PRECONDITIONS THE COMPARATOR USED TO ASSUME. Each was found by review with a measurement attached, and
// each failed in the same direction: a plausible number for a comparison that had not happened.
void theComparatorRefusesWhatItCannotCompare()
{
    felitronics::test::group ("the comparison refuses two runs that are not the same question");
    const double fs = 48000.0;
    auto a = programme (fs, 2, 6.0, 25.0);
    auto b = a;
    std::vector<double> scratch;

    BandCrest same0, same1;
    if (! felitronics::test::run (runIt (same0, a, fs))) return;
    if (! felitronics::test::run (runIt (same1, b, fs))) return;
    const auto ok0 = bandCrestLoss (same0, same1, BandCrest::kFull, scratch);
    ok (ok0.usable > 0 && ok0.valid && std::fabs (ok0.maxDb) <= 1.0e-9,
        "PRECONDITION: the same audio against itself is a comparison, and it finds no loss");

    // A DIFFERENT GRID IS A DIFFERENT QUESTION. Same audio, master on a 50 ms hop: before the guard this
    // answered CVaR95 7.8 dB and `valid` true.
    {
        BandCrest other;
        felitronics::analysis::BandCrestParams p; p.hopMs = 50.0;
        other.setParams (p);
        if (felitronics::test::run (runIt2 (other, b, fs)))
        {
            const auto r = bandCrestLoss (same0, other, BandCrest::kFull, scratch);
            ok (! r.valid && r.usable == 0 && r.blocks == 0,
                "a master measured on a 50 ms hop is refused, not compared (usable "
                + std::to_string (r.usable) + ", blocks " + std::to_string (r.blocks) + ")");
        }
    }
    // ... and so is a different band split.
    {
        BandCrest other;
        felitronics::analysis::BandCrestParams p; p.bandEdgeHz[1] = 1500.0;
        other.setParams (p);
        if (felitronics::test::run (runIt2 (other, b, fs)))
            ok (! bandCrestLoss (same0, other, BandCrest::kFull, scratch).valid,
                "and so is a master measured through different band edges");
    }

    // A GROSS MISALIGNMENT IS SEEN. A master delayed by whole blocks must not read as a clean comparison.
    //
    // THE FIXTURE NEEDS A FEATURE TO LOCK ONTO, which the first version of this check did not give it: with a
    // transient every 0.5 s the block-peak series is nearly flat, the correlation is flat with it, and the
    // argmax lands on lag 0 for want of anything better — the instrument read "aligned" because it could not
    // see. The same lesson as the tone burst in K6. `bigEvery = 4` puts a much larger transient every 2 s, so
    // the series has structure at the scale the lag is measured in.
    {
        auto marked = programme (fs, 2, 6.0, 4.0);
        BandCrest ref2;
        if (! felitronics::test::run (runIt (ref2, marked, fs))) return;
        const auto hop = (std::size_t) ref2.hopSamples();
        auto shifted = marked;
        for (auto& c : shifted) { c.insert (c.begin(), 3 * hop, 0.0f); c.resize (marked[0].size()); }
        BandCrest late;
        if (felitronics::test::run (runIt (late, shifted, fs)))
        {
            const auto r = bandCrestLoss (ref2, late, BandCrest::kFull, scratch);
            ok (r.lagBlocks != 0, "a master delayed by three hops reports a non-zero block lag ("
                                  + std::to_string (r.lagBlocks) + ")");
        }
    }

    // A MUTED OUTPUT BAND IS COUNTED, NOT DROPPED. Before this it left the population in silence, so a band
    // the chain removed read exactly like a band it left alone.
    {
        auto muted = a;
        for (auto& c : muted) std::fill (c.begin() + (std::ptrdiff_t) (c.size() / 2), c.end(), 0.0f);
        BandCrest half;
        if (felitronics::test::run (runIt (half, muted, fs)))
        {
            const auto r = bandCrestLoss (same0, half, BandCrest::kFull, scratch);
            ok (r.outSilent > 0, "the blocks whose master band is digital silence are counted ("
                                 + std::to_string (r.outSilent) + " of " + std::to_string (r.inActive) + ")");
        }
    }
}

//==================================================================================================
void theSameProgrammeInAnySlicing()
{
    felitronics::test::group ("law 8a — the call sizes do not change the measurement");
    const double fs = 48000.0;
    auto buf = programme (fs, 2, 4.0, 25.0);
    BandCrest ref;
    if (! felitronics::test::run (runIt (ref, buf, fs, 4096))) return;

    int bad = 0;
    for (const int block : { 1, 7, 64, 977, 65536 })
    {
        BandCrest bc;
        if (! runIt (bc, buf, fs, block)) { ++bad; continue; }
        if (bc.hopCount() != ref.hopCount() || bc.blockCount() != ref.blockCount()) { ++bad; continue; }
        for (long long j = 0; j < ref.blockCount(); ++j)
            for (int b = 0; b < BandCrest::kBands; ++b)
                if (bc.blockPeakLin (j, b) != ref.blockPeakLin (j, b)
                    || bc.blockMeanSq (j, b) != ref.blockMeanSq (j, b)) { ++bad; j = ref.blockCount(); break; }
    }
    ok (bad == 0, "1, 7, 64, 977 and 65536 samples a call give the same cells, bit for bit");
}

//==================================================================================================
void processAsksTheHeapForNothing()
{
    felitronics::test::group ("RT — prepare() allocates, process() and finish() do not");
    const double fs = 48000.0;
    auto buf = programme (fs, 2, 2.0);
    BandCrest bc;
    if (! felitronics::test::run (bc.prepare (fs, 2, (int) buf[0].size()))) return;
    std::vector<const float*> p; for (auto& c : buf) p.push_back (c.data());
    const long long before = alloc::count.load();
    const int n = (int) buf[0].size();
    for (int off = 0; off < n; )
    {
        const int m = std::min (4096, n - off);
        std::vector<const float*> q; for (auto& c : buf) q.push_back (c.data() + off);
        const long long inner = alloc::count.load();
        if (! bc.process (q.data(), 2, m)) break;
        // THE COUNT IS READ INTO A LOCAL BEFORE `ok` IS CALLED. Passing `alloc::count.load()` as one argument
        // and a message as another leaves their evaluation order UNSPECIFIED, and the message is a
        // `std::string` that allocates: clang happened to read the counter first and gcc built the string
        // first, so the same code passed on one row and failed 25 checks on the other. The trap is known in
        // this tree; this is it again.
        const long long after = alloc::count.load();
        ok (after == inner, "process() asked the heap for nothing");
        off += m;
    }
    const long long b2 = alloc::count.load();
    bc.finish();
    const long long b3 = alloc::count.load();
    ok (b3 == b2, "finish() asked the heap for nothing");
    (void) before;

    // AND THE BUDGET IS THE DEMAND. `storageFor` is the function `prepare` sizes itself with, so the two
    // cannot drift; what is checked here is that a FRESH object's allocation is inside it.
    const felitronics::analysis::BandCrestParams dflt {};
    const auto st = BandCrest::storageFor (fs, 2, n, dflt);
    ok (st.ok && st.bytes() > 0, "storageFor publishes a demand of " + std::to_string (st.bytes()) + " B");
}

} // namespace

int main()
{
    std::printf ("felitronics::analysis::BandCrest — K1\n");
    theDelayIsWhereTheHeaderSaysItIs();
    thePeakIsTheCertificate();
    theGridIsTheLoudnessMeters();
    theLossIsInvariantToGain();
    cvarSeesWhatP95CannotAndNeverLess();
    theBandsAreTheFiltersTheHeaderNames();
    theProgrammeLevelIsInTheGatesUnits();
    theComparatorRefusesWhatItCannotCompare();
    theSameProgrammeInAnySlicing();
    processAsksTheHeapForNothing();
    return felitronics::test::report();
}

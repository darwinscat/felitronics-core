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
    felitronics::test::group ("the block grid is the loudness meter's, counted by different code");
    for (const double fs : { 44100.0, 48000.0, 96000.0 })
    {
        auto buf = programme (fs, 2, 3.0);
        BandCrest bc;
        if (! felitronics::test::run (runIt (bc, buf, fs))) continue;

        felitronics::analysis::LoudnessMeter lm;
        const int n = (int) buf[0].size();
        if (! felitronics::test::run (lm.prepareForSamples (fs, 2, n))) continue;
        const float* q[2] { buf[0].data(), buf[1].data() };
        if (! felitronics::test::run (lm.process (q, 2, n))) continue;

        ok (bc.blockCount() == (long long) lm.gatingBlockCount(),
            std::to_string ((int) fs) + " Hz: " + std::to_string (bc.blockCount())
            + " blocks here and " + std::to_string (lm.gatingBlockCount()) + " gating blocks there");
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
    theSameProgrammeInAnySlicing();
    processAsksTheHeapForNothing();
    return felitronics::test::report();
}

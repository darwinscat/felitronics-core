// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.
//
// K14 — the Side-only high shelf that rides inside MonoBass's M/S island.
//
// THE FIXTURES ARE CHOSEN SO THAT ONE WRONG IMPLEMENTATION DIES IN EACH, not so that a correct one
// passes. Side-only content (L = x, R = -x) makes Mid exactly zero, so the measured gain IS |H| and a
// filter put on Mid, or on L and R, reads 0 dB there. Hard-panned content makes Mid and Side equal, so
// the OTHER channel shows the bleed a Side boost necessarily causes and an L/R shelf never does. And
// the numbers are computed from the tree's own designer, not quoted from anywhere.

#include <felitronics_test.h>

#include <felitronics/eq/MatchedBiquad.h>
#include <felitronics/stereo/MonoBass.h>

#include <cmath>
#include <complex>
#include <cstdio>
#include <random>
#include <vector>

namespace
{
using felitronics::stereo::MonoBass;
using felitronics::stereo::MonoBassParams;
using felitronics::stereo::StereoAirParams;
namespace test = felitronics::test;
using test::ok;
using test::approx;

constexpr double kPi = 3.14159265358979323846;
constexpr double kFs = 48000.0;

// The amplitude of `hz` in x, by projection — the same trick the stereo suite uses, and it ignores
// everything that is not that tone.
double toneAmp (const std::vector<float>& x, double hz, double fs, std::size_t from)
{
    double re = 0.0, im = 0.0; std::size_t n = 0;
    for (std::size_t i = from; i < x.size(); ++i, ++n)
    {
        const double w = 2.0 * kPi * hz * (double) i / fs;
        re += (double) x[i] * std::cos (w);
        im += (double) x[i] * std::sin (w);
    }
    return n > 0 ? 2.0 * std::sqrt (re * re + im * im) / (double) n : 0.0;
}

// |H| of the tree's own shelf, so the expectation and the code share one definition of the filter and
// differ only in the path the samples took.
double shelfDb (double f0, double gainDb, double f, double fs)
{
    const auto c = felitronics::eq::matched::highShelfDb (f0, fs, gainDb);
    const std::complex<double> z = std::exp (std::complex<double> (0.0, -2.0 * kPi * f / fs));
    return 20.0 * std::log10 (std::abs ((c.b0 + c.b1 * z + c.b2 * z * z) / (1.0 + c.a1 * z + c.a2 * z * z)));
}

struct Rig
{
    MonoBass mb;
    bool prepare (const MonoBassParams& bp, const StereoAirParams& ap, double fs = kFs)
    {
        if (! mb.prepare (fs, 4096, 2)) return false;
        mb.setParams (bp);
        mb.setAir (ap);
        return true;
    }
    void run (std::vector<float>& l, std::vector<float>& r, int block = 512)
    {
        float* io[2] { l.data(), r.data() };
        for (std::size_t i = 0; i < l.size(); i += (std::size_t) block)
        {
            const int k = (int) std::min ((std::size_t) block, l.size() - i);
            float* p[2] { l.data() + i, r.data() + i };
            (void) mb.process (p, 2, k);
            (void) io;
        }
    }
};

MonoBassParams bassOff() { MonoBassParams p; p.enabled = false; return p; }
} // namespace

int main()
{
    std::printf ("felitronics::stereo K14 — the Side-only air shelf\n");
    const std::size_t n = 24000;

    test::group ("Side-only content: the measured gain IS the shelf, and the bass is not in the way");
    {
        // L = x, R = -x  ->  Mid exactly 0, Side = x. A shelf on Mid, or on L and R, reads 0 dB here.
        std::vector<float> l (n), r (n);
        for (std::size_t i = 0; i < n; ++i)
        {
            const double t = (double) i / kFs;
            const float v = (float) (0.3 * std::sin (2.0 * kPi * 1000.0 * t)
                                   + 0.3 * std::sin (2.0 * kPi * 10000.0 * t));
            l[i] = v; r[i] = -v;
        }
        const double in1 = toneAmp (l, 1000.0, kFs, 0), in10 = toneAmp (l, 10000.0, kFs, 0);
        StereoAirParams ap; ap.enabled = true; ap.frequencyHz = 6000.0f; ap.gainDb = 3.0f;
        Rig rig;
        if (test::run (rig.prepare (bassOff(), ap)))
        {
            rig.run (l, r);
            const std::size_t skip = 4800;              // past the 20 ms ramp
            const double g1  = 20.0 * std::log10 (toneAmp (l, 1000.0,  kFs, skip) / in1);
            const double g10 = 20.0 * std::log10 (toneAmp (l, 10000.0, kFs, skip) / in10);
            approx (g1,  shelfDb (6000.0, 3.0, 1000.0,  kFs), 0.02,
                    "1 kHz is essentially untouched (" + std::to_string (g1) + " dB)");
            approx (g10, shelfDb (6000.0, 3.0, 10000.0, kFs), 0.02,
                    "10 kHz gets the shelf (" + std::to_string (g10) + " dB)");
            ok (g10 > 2.5 && g10 < 2.8,
                "…and that is +2.65, NOT the +3 the parameter says: the plateau is above the corner");
        }
    }

    test::group ("a hard-panned top bleeds into the other channel, inverted — no L/R shelf does this");
    {
        // L = x, R = 0  ->  Mid = Side = x/2, so R' = m - H*s and the other channel stops being silent.
        std::vector<float> l (n, 0.0f), r (n, 0.0f);
        for (std::size_t i = 0; i < n; ++i)
            l[i] = (float) (0.4 * std::sin (2.0 * kPi * 10000.0 * (double) i / kFs));
        StereoAirParams ap; ap.enabled = true; ap.frequencyHz = 6000.0f; ap.gainDb = 3.0f;
        Rig rig;
        if (test::run (rig.prepare (bassOff(), ap)))
        {
            rig.run (l, r);
            const std::size_t skip = 4800;
            const double aL = toneAmp (l, 10000.0, kFs, skip), aR = toneAmp (r, 10000.0, kFs, skip);
            ok (aR > 0.0 && aL > 0.0, "both channels carry the tone now");
            const double bleedDb = 20.0 * std::log10 (aR / aL);
            approx (bleedDb, -15.39, 0.3,
                    "the silent channel comes up to -15.4 dB of the loud one (" + std::to_string (bleedDb) + ")");
        }
    }

    test::group ("the mono fold: within one ulp, and bit-exact when there is no Side at all");
    {
        std::mt19937 rng (7); std::uniform_real_distribution<float> d (-0.7f, 0.7f);
        std::vector<float> l0 (n), r0 (n);
        for (std::size_t i = 0; i < n; ++i) { l0[i] = d (rng); r0[i] = d (rng); }
        std::vector<float> l = l0, r = r0;
        StereoAirParams ap; ap.enabled = true; ap.frequencyHz = 6000.0f; ap.gainDb = 6.0f;
        Rig rig;
        if (test::run (rig.prepare (bassOff(), ap)))
        {
            rig.run (l, r);
            double worst = 0.0;
            for (std::size_t i = 0; i < n; ++i)
                worst = std::max (worst, std::fabs (((double) l[i] + r[i]) - ((double) l0[i] + r0[i])));
            // ONE ULP, not zero, and not a percentage. (M+S)+(M-S) rounds each term separately, so the
            // sum cannot be bit-exact; the BOUND is what has a meaning at any level.
            ok (worst <= 1.2e-07,
                "the mono sum moves by at most one ulp (" + std::to_string (worst) + ")");

            // …and where Side is exactly zero the shelf of a zero stream is exactly zero, so the
            // reconstruction is the input, bit for bit.
            std::vector<float> lm (n), rm (n);
            for (std::size_t i = 0; i < n; ++i) lm[i] = rm[i] = l0[i];
            const std::vector<float> keep = lm;
            Rig rig2;
            if (test::run (rig2.prepare (bassOff(), ap)))
            {
                rig2.run (lm, rm);
                std::size_t differ = 0;
                for (std::size_t i = 0; i < n; ++i) if (! (lm[i] == keep[i]) || ! (rm[i] == keep[i])) ++differ;
                ok (differ == 0, "a mono programme comes back bit-identical (" + std::to_string (differ) + " differ)");
            }
        }
    }

    test::group ("zero is a BRANCH: it must be bit-identical to the tool being off");
    {
        std::mt19937 rng (11); std::uniform_real_distribution<float> d (-0.7f, 0.7f);
        std::vector<float> l0 (n), r0 (n);
        for (std::size_t i = 0; i < n; ++i) { l0[i] = d (rng); r0[i] = d (rng); }
        MonoBassParams bp; bp.enabled = true; bp.frequencyHz = 120.0f; bp.lowWidth = 0.0f;

        auto render = [&] (const StereoAirParams& ap)
        {
            std::vector<float> l = l0, r = r0;
            Rig rig; if (! rig.prepare (bp, ap)) return std::vector<float> {};
            rig.run (l, r);
            std::vector<float> out = l; out.insert (out.end(), r.begin(), r.end());
            return out;
        };
        StereoAirParams off;                                   // enabled = false
        StereoAirParams zero; zero.enabled = true; zero.gainDb = 0.0f;
        StereoAirParams neg;  neg.enabled = true;  neg.gainDb = -2.0f;   // clamps to 0
        const auto a = render (off), b = render (zero), c = render (neg);
        std::size_t d1 = 0, d2 = 0;
        for (std::size_t i = 0; i < a.size() && i < b.size(); ++i) if (! (a[i] == b[i])) ++d1;
        for (std::size_t i = 0; i < a.size() && i < c.size(); ++i) if (! (a[i] == c[i])) ++d2;
        ok (! a.empty() && d1 == 0, "0 dB with the flag on is bit-identical to the flag off (" + std::to_string (d1) + ")");
        ok (! a.empty() && d2 == 0, "…and so is a negative request, which clamps to 0 (" + std::to_string (d2) + ")");
    }

    test::group ("the island belongs to both tools, and each keeps its own latch");
    {
        // The bass settled full-wide is the case that used to leave the loop for the rest of the block.
        std::vector<float> l (n), r (n);
        for (std::size_t i = 0; i < n; ++i)
        {
            const float v = (float) (0.3 * std::sin (2.0 * kPi * 10000.0 * (double) i / kFs));
            l[i] = v; r[i] = -v;
        }
        const double in10 = toneAmp (l, 10000.0, kFs, 0);
        MonoBassParams bp; bp.enabled = true; bp.frequencyHz = 120.0f; bp.lowWidth = 1.0f;   // settles into bypass
        StereoAirParams ap; ap.enabled = true; ap.frequencyHz = 6000.0f; ap.gainDb = 3.0f;
        const std::vector<float> l0 = l, r0 = r;

        // A PROJECTION OVER THE WHOLE TAIL CANNOT SEE THIS, and finding that out is why the planted
        // control is run first. The old mid-loop `return` drops the shelf only for the REMAINDER OF THE
        // BLOCK in which the bass settles — about a hundred samples out of twenty-four thousand — and a
        // tone projection averages them away completely. What it cannot hide from is RE-SLICING: the
        // number of dropped samples depends on where the caller cut, so the same programme through
        // different block sizes stops matching. That is law 8a, and it is the shape this defect has.
        auto renderAt = [&] (int block)
        {
            std::vector<float> a = l0, b2 = r0;
            Rig rg; if (! rg.prepare (bp, ap)) return std::vector<float> {};
            rg.run (a, b2, block);
            return a;
        };
        const auto small = renderAt (128), whole = renderAt ((int) n);
        double worst = 0.0;
        for (std::size_t i = 0; i < small.size() && i < whole.size(); ++i)
            worst = std::max (worst, std::fabs ((double) small[i] - (double) whole[i]));
        ok (! small.empty() && worst < 1.0e-5,
            "the same programme through 128-sample and whole-file calls agrees (worst "
                + std::to_string (worst) + ") — a shelf dropped for the rest of a block would not");

        Rig rig;
        if (test::run (rig.prepare (bp, ap)))
        {
            rig.run (l, r, 128);                     // small blocks: the settle lands mid-loop
            const double g10 = 20.0 * std::log10 (toneAmp (l, 10000.0, kFs, 8000) / in10);
            approx (g10, shelfDb (6000.0, 3.0, 10000.0, kFs), 0.05,
                    "and the shelf is still there after the bass went full-wide (" + std::to_string (g10) + " dB)");
        }
    }

    test::group ("the corner is clamped against the rate, and the clamp is readable");
    {
        MonoBass mb;
        if (test::run (mb.prepare (8000.0, 512, 2)))
        {
            StereoAirParams ap; ap.enabled = true; ap.frequencyHz = 12000.0f; ap.gainDb = 3.0f;
            mb.setAir (ap);
            approx ((double) mb.air().frequencyHz, 3600.0, 1e-4,
                    "at 8 kHz the corner clamps to 0.45 fs = 3600 Hz");
            ap.frequencyHz = 100.0f; mb.setAir (ap);
            approx ((double) mb.air().frequencyHz, 3000.0, 1e-4, "…and up to the 3 kHz floor from below");
            ap.gainDb = 99.0f; mb.setAir (ap);
            approx ((double) mb.air().gainDb, (double) MonoBass::kMaxAirDb, 1e-6,
                    "the plateau clamps to the documented maximum, and reads back clamped");
            ap.gainDb = std::nan (""); mb.setAir (ap);
            approx ((double) mb.air().gainDb, (double) MonoBass::kMaxAirDb, 1e-6,
                    "a non-finite request keeps the last good value, as every parameter here does");
        }
    }

    test::group ("the width report: three energies, because the fraction is blind where it matters");
    {
        StereoAirParams ap; ap.enabled = true; ap.frequencyHz = 6000.0f; ap.gainDb = 3.0f;

        // (a) A MONO TOP. Side is exactly zero, so the shelf of a zero stream is exactly zero and both
        //     Side energies are exactly 0 — not "small". Width is 0, which is a READING, not a refusal.
        {
            std::vector<float> l (n), r (n);
            std::mt19937 rng (3); std::uniform_real_distribution<float> d (-0.5f, 0.5f);
            for (std::size_t i = 0; i < n; ++i) { l[i] = d (rng); r[i] = l[i]; }
            Rig rig;
            if (test::run (rig.prepare (bassOff(), ap)))
            {
                rig.run (l, r);
                ok (rig.mb.airSideEnergyBefore() == 0.0 && rig.mb.airSideEnergyAfter() == 0.0,
                    "a mono top has exactly zero Side energy, before and after");
                ok (rig.mb.airMidEnergy() > 0.0, "…while the Mid band is not empty");
                approx (rig.mb.airWidthBefore(), 0.0, 1e-12, "width reads 0, which is a reading");
                approx (rig.mb.airWidthAfter(),  0.0, 1e-12, "…and still 0 afterwards");
            }
        }

        // (b) AN ANTI-PHASE TOP — the case the whole three-energy decision exists for. Mid is exactly 0,
        //     so BOTH widths read 1.000 and neither moves, while the Side energy grows by the band's own
        //     integral of the shelf. A page publishing only the fraction would report "nothing happened".
        {
            std::vector<float> l (n), r (n);
            std::mt19937 rng (5); std::uniform_real_distribution<float> d (-0.5f, 0.5f);
            for (std::size_t i = 0; i < n; ++i) { l[i] = d (rng); r[i] = -l[i]; }
            Rig rig;
            if (test::run (rig.prepare (bassOff(), ap)))
            {
                // MEASURED AFTER THE RAMP, and the difference is not small enough to wave away: the gain
                // glides over 20 ms, so 960 of 24000 samples carry less than the plateau and the ratio
                // comes out 1.1 % low — 1.892 against 1.9135, which is exactly what the ramp's own
                // average predicts. `reset()` snaps the smoother to its target and starts the width
                // interval over, which is what a caller wanting a settled reading does too.
                std::vector<float> warm (2048, 0.05f), warmR (2048, -0.05f);
                rig.run (warm, warmR);
                rig.mb.reset();
                rig.run (l, r);
                approx (rig.mb.airWidthBefore(), 1.0, 1e-9, "anti-phase reads width 1.000 before");
                approx (rig.mb.airWidthAfter(),  1.0, 1e-9, "…and exactly 1.000 after: the fraction is BLIND here");
                const double ratio = rig.mb.airSideEnergyAfter() / rig.mb.airSideEnergyBefore();
                // White noise through an LR4 high-pass at the corner, weighted by the shelf: computed
                // from the tree's own designer, and NOT the plateau (1.9953) — the band includes the skirt.
                approx (ratio, 1.9135, 0.02,
                        "…while the energies say it plainly: x" + std::to_string (ratio));
                ok (ratio < 1.96, "and that is BELOW the plateau — the band is a weighting, not a wall");
            }
        }

        // (c) SIDE CONTENT BELOW THE BAND. An LR4 high-pass at 6 kHz rejects 500 Hz by tens of dB, so a
        //     wide bass with a mono top must not read as a wide top.
        {
            std::vector<float> l (n), r (n);
            for (std::size_t i = 0; i < n; ++i)
            {
                const double t = (double) i / kFs;
                const float side = (float) (0.4 * std::sin (2.0 * kPi * 500.0 * t));
                const float mid  = (float) (0.4 * std::sin (2.0 * kPi * 10000.0 * t));
                l[i] = mid + side; r[i] = mid - side;
            }
            Rig rig;
            if (test::run (rig.prepare (bassOff(), ap)))
            {
                rig.run (l, r);
                ok (rig.mb.airWidthBefore() >= 0.0 && rig.mb.airWidthBefore() < 1e-3,
                    "a 500 Hz Side under a 10 kHz Mid reads a NARROW top ("
                        + std::to_string (rig.mb.airWidthBefore()) + ")");
            }
        }

        // (d) WITH THE TOOL OFF nothing is judged — which is not the same as judging zero.
        {
            std::vector<float> l (n, 0.1f), r (n, -0.1f);
            StereoAirParams off;
            Rig rig;
            if (test::run (rig.prepare (bassOff(), off)))
            {
                rig.run (l, r);
                ok (rig.mb.airJudgedSamples() == 0, "no samples were judged");
                ok (rig.mb.airWidthBefore() < 0.0, "and the width refuses rather than answering 0.0");
            }
        }
    }

    return test::report();
}

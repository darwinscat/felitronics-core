// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// codecgrid's curve and its reading. The level accumulator is checked against a logarithm taken per coefficient
// (it takes one per sixteen frames); the reading is checked on curves built to have exactly one property each —
// a dip in every cell, a dip in one cell, a dip repeated at many offsets, a dip at the very end of the range.

#include <felitronics/codecgrid/GridCurve.h>

#include "felitronics_test.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>

using namespace felitronics::codecgrid;
using felitronics::test::approx;
using felitronics::test::group;
using felitronics::test::ok;

namespace
{
struct Lcg
{
    std::uint32_t s;
    double next() noexcept { s = s * 1664525u + 1013904223u; return (double) (s >> 8) / 16777216.0; }   // [0, 1)
    double gauss() noexcept { double a = 0.0; for (int i = 0; i < 12; ++i) a += next(); return a - 6.0; }
};

struct Reader
{
    std::vector<double> doubles;
    std::vector<float> floats;
    ReadScratch scratch;
    explicit Reader (int offsets) : doubles (ReadScratch::doublesFor (offsets)), floats (ReadScratch::floatsFor (offsets))
    {
        scratch.local = doubles.data();
        scratch.sort = doubles.data() + offsets;
        scratch.dip = floats.data();
    }
};

// A curve of `offsets` x 32 cells: a level per cell plus small noise.
std::vector<float> flatCurve (int offsets, std::uint32_t seed, double noiseDb)
{
    std::vector<float> c ((std::size_t) offsets * kCells);
    Lcg rng { seed };
    for (int o = 0; o < offsets; ++o)
        for (int k = 0; k < kCells; ++k) c[(std::size_t) o * kCells + (std::size_t) k] = (float) (-20.0 - k + noiseDb * rng.gauss());
    return c;
}
}

int main()
{
    std::printf ("felitronics::codecgrid curve tests\n");

    group ("the helpers: median, mirrored neighbour, phases on a circle");
    {
        double a[5] { 5.0, 1.0, 4.0, 2.0, 3.0 }, b[4] { 4.0, 1.0, 3.0, 2.0 }, one[1] { 7.0 };
        approx (detail::medianOf (a, 5), 3.0, 0.0, "odd count: the middle value");
        approx (detail::medianOf (b, 4), 2.5, 0.0, "even count: the mean of the two middle values");
        approx (detail::medianOf (one, 1), 7.0, 0.0, "one value: itself");
        approx (detail::medianOf (one, 0), 0.0, 0.0, "no values: zero, not a read past the array");
        ok (detail::mirrored (-1, 10) == 1 && detail::mirrored (-8, 10) == 8, "below the range: reflected about offset 0, the end not repeated");
        ok (detail::mirrored (10, 10) == 8 && detail::mirrored (17, 10) == 1, "above the range: reflected about the last offset");
        ok (detail::mirrored (0, 10) == 0 && detail::mirrored (9, 10) == 9 && detail::mirrored (4, 10) == 4, "inside the range: itself");
        ok (detail::mirrored (5, 1) == 0, "a range of one offset has one neighbour");
        ok (phasesAgree (0, 1023, 1024, 1) && phasesAgree (1023, 0, 1024, 1), "0 and hop - 1 are neighbours on the circle");
        ok (! phasesAgree (0, 2, 1024, 1) && phasesAgree (0, 2, 1024, 2), "the tolerance is inclusive and no wider");
        ok (phasesAgree (575, 575 + 576, 576, 0) && phasesAgree (-1, 575, 576, 0), "phases are compared modulo the hop, negatives included");
    }

    group ("the accumulator gives the mean of logarithms, with a logarithm per cell and per 64 frames");
    {
        const int bins = 2048, frames = 150;                     // 256 coefficients a cell: more than one mantissa batch; 150 frames: two full log batches and a short one, not a multiple of four
        std::vector<float> floats (LevelAccumulator::floatsFor (bins));
        std::vector<double> doubles (LevelAccumulator::doublesFor (bins));
        LevelAccumulator acc;
        acc.attach (floats.data(), doubles.data(), bins);

        Lcg rng { 99u };
        std::vector<float> l ((std::size_t) frames * bins), r ((std::size_t) frames * bins);
        for (auto& v : l) v = (float) (rng.gauss() * 0.1);
        for (auto& v : r) v = (float) (rng.gauss() * 0.05);
        // a few true zeros and a few huge values: the floor and the product must both hold
        for (int f = 0; f < frames; f += 3) { l[(std::size_t) f * bins + 5] = 0.0f; r[(std::size_t) f * bins + 5] = 0.0f; }
        l[7] = 1.0e6f;

        acc.beginReference();
        for (int f = 0; f < frames; ++f) acc.addReferenceFrame (l.data() + (std::size_t) f * bins, r.data() + (std::size_t) f * bins);
        ok (acc.finishReference(), "a reference over 37 frames is accepted");
        acc.beginOffset();
        for (int f = 0; f < frames; ++f) acc.addFrame (l.data() + (std::size_t) f * bins, r.data() + (std::size_t) f * bins);
        float cells[kCells];
        acc.finishOffset (cells);
        ok (acc.frames() == frames, "every frame was counted");

        // the same thing with a logarithm per coefficient, written independently
        double worst = 0.0;
        for (int s = 0; s < kSignals; ++s)
            for (int g = 0; g < kGroups; ++g)
            {
                double sum = 0.0;
                for (int k = g * (bins / kGroups); k < (g + 1) * (bins / kGroups); ++k)
                {
                    double power = 0.0;
                    std::vector<double> v ((std::size_t) frames);
                    for (int f = 0; f < frames; ++f)
                    {
                        const float a = l[(std::size_t) f * bins + (std::size_t) k], b = r[(std::size_t) f * bins + (std::size_t) k];
                        v[(std::size_t) f] = s == 0 ? a : s == 1 ? b : s == 2 ? 0.5f * (a + b) : 0.5f * (a - b);       // mid and side are single precision
                        power += v[(std::size_t) f] * v[(std::size_t) f];
                    }
                    const double rms = std::sqrt (power / frames) + 1.0e-30;
                    for (int f = 0; f < frames; ++f) sum += 20.0 * std::log10 (std::fmin (std::fmax (std::fabs (v[(std::size_t) f]) / rms, kLevelFloor), kLevelCeiling));
                }
                const double want = sum / ((bins / kGroups) * frames);
                worst = std::fmax (worst, std::fabs ((double) cells[s * kGroups + g] - want));
            }
        std::printf ("    worst difference from a logarithm per coefficient: %.3e dB\n", worst);
        ok (std::isfinite (worst) && worst < 1.0e-4, "the batched product is the same mean (within single precision of the ratios)");
        bool finite = true;
        for (float c : cells) finite = finite && std::isfinite (c);
        ok (finite, "true zeros and a huge value leave every cell finite");

        // the floor: a frame of exact zeros reads -80 dB in every cell, not -infinity
        acc.beginOffset();
        std::vector<float> z ((std::size_t) bins, 0.0f);
        acc.addFrame (z.data(), z.data());
        acc.finishOffset (cells);
        bool floor = true;
        for (float c : cells) floor = floor && std::fabs ((double) c + 80.0) < 1.0e-9;
        ok (floor, "an all-zero frame reads exactly the floor, -80 dB");

        // a NaN is the floor, never kept
        acc.beginOffset();
        std::vector<float> bad ((std::size_t) bins, std::nanf (""));
        acc.addFrame (bad.data(), bad.data());
        acc.finishOffset (cells);
        bool nanFloor = true;
        for (float c : cells) nanFloor = nanFloor && std::isfinite (c) && std::fabs ((double) c + 80.0) < 1.0e-9;
        ok (nanFloor, "a frame of NaN reads the floor in every cell, and nothing non-finite comes out");

        // no frame: zeros, and a reference without frames is refused
        acc.beginOffset();
        acc.finishOffset (cells);
        bool zero = true;
        for (float c : cells) zero = zero && std::fpclassify (c) == FP_ZERO;
        ok (zero, "an offset without frames reads zero in every cell");
        acc.beginReference();
        ok (! acc.finishReference(), "a reference without frames is refused");

        // the zero share: the coefficient zeroed in every third frame is counted there and nowhere else
        acc.beginReference();
        for (int f = 0; f < frames; ++f) acc.addReferenceFrame (l.data() + (std::size_t) f * bins, r.data() + (std::size_t) f * bins);
        ok (acc.finishReference(), "the reference again");
        std::vector<std::uint32_t> counts ((std::size_t) kSignals * kZeroBands, 0u);
        for (int f = 0; f < frames; ++f) acc.addZeroFrame (l.data() + (std::size_t) f * bins, r.data() + (std::size_t) f * bins, counts.data());
        const int zeroed = (frames + 2) / 3, band = 5 / (bins / kZeroBands);
        ok ((int) counts[(std::size_t) band] >= zeroed && (int) counts[(std::size_t) (2 * kZeroBands + band)] >= zeroed,
            "the zeroed coefficient is counted in its band for left and for mid");
    }

    group ("the reading: one dip in every cell is found where it is, and stands out in every cell");
    {
        for (int offsets : { 576, 960, 1024 })
            for (int at : { 0, 1, 300, offsets - 1 })
            {
                auto c = flatCurve (offsets, 1000u + (std::uint32_t) at, 0.05);
                for (int k = 0; k < kCells; ++k) c[(std::size_t) at * kCells + (std::size_t) k] -= 3.0f;
                Reader rd (offsets);
                const GridReading g = readCurve (c.data(), offsets, rd.scratch);
                const std::string tag = "offsets " + std::to_string (offsets) + ", dip at " + std::to_string (at) + ": ";
                ok (g.offset == at, tag + "the offset is found, the ends of the range included");
                ok (g.cells == kCells, tag + "it stands out in all 32 cells");
                approx (g.localDipSum, 3.0 * kCells, 1.5, tag + "the summed dip is 32 x 3 dB");
                ok (g.score > 50.0 && g.second < 0.2 * g.score, tag + "a large score and no runner-up");
                ok (GridRule {}.found (g) && ! GridRule {}.several (g), tag + "the rule finds it");
            }
    }

    group ("the reading: what must NOT be found");
    {
        const int offsets = 1024;
        Reader rd (offsets);
        // noise alone
        int worstCells = 0;
        double worstScore = 0.0;
        for (std::uint32_t seed = 1; seed <= 40; ++seed)
        {
            const auto c = flatCurve (offsets, seed, 0.05);
            const GridReading g = readCurve (c.data(), offsets, rd.scratch);
            worstScore = std::fmax (worstScore, g.score);
            worstCells = std::max (worstCells, g.cells);
            ok (! GridRule {}.found (g), "noise seed " + std::to_string (seed) + " is not a grid");
        }
        std::printf ("    40 noise curves: largest score %.2f, most cells %d\n", worstScore, worstCells);

        // a dip in ONE cell, however deep: the score may be large, the breadth is not
        auto one = flatCurve (offsets, 7u, 0.05);
        one[(std::size_t) 500 * kCells + 3] -= 40.0f;
        const GridReading g1 = readCurve (one.data(), offsets, rd.scratch);
        ok (g1.offset == 500 && g1.cells == 1, "a deep dip in one cell is seen at its offset, in one cell");
        ok (g1.score > 9.2 && g1.localDipSum > 5.0 && ! GridRule {}.found (g1), "score and depth pass, the breadth condition refuses it");

        // the same dip at every 64th offset: broad, deep — and not unique
        auto many = flatCurve (offsets, 8u, 0.05);
        for (int o = 10; o < offsets; o += 64)
            for (int k = 0; k < kCells; ++k) many[(std::size_t) o * kCells + (std::size_t) k] -= 3.0f;
        const GridReading gm = readCurve (many.data(), offsets, rd.scratch);
        ok (gm.second > 0.9 * gm.score, "a repeated dip has a runner-up as strong as the winner");
        ok (! GridRule {}.found (gm) && GridRule {}.several (gm), "the uniqueness condition refuses it, and it is named `several`");

        // a dip three offsets wide is one dip, not a winner and a runner-up: the neighbours within two samples
        // are not candidates for `second`
        auto wide = flatCurve (offsets, 9u, 0.05);
        for (int o = 699; o <= 701; ++o)
            for (int k = 0; k < kCells; ++k) wide[(std::size_t) o * kCells + (std::size_t) k] -= 3.0f;
        const GridReading gw = readCurve (wide.data(), offsets, rd.scratch);
        ok (gw.offset >= 699 && gw.offset <= 701 && GridRule {}.found (gw), "a dip three offsets wide is found as one");

        // a tiny dip, exactly repeated in every cell: many sigmas of nothing
        std::vector<float> still ((std::size_t) offsets * kCells, -30.0f);
        for (int k = 0; k < kCells; ++k) still[(std::size_t) 200 * kCells + (std::size_t) k] -= 0.01f;
        const GridReading gs = readCurve (still.data(), offsets, rd.scratch);
        ok (gs.offset == 200 && gs.score > 1.0e3, "on a noiseless curve a 0.01 dB dip is an enormous score");
        ok (! GridRule {}.found (gs), "...and the absolute depth condition refuses it");

        // too few offsets to have a neighbourhood: an empty reading, not a read past the curve
        const GridReading gz = readCurve (still.data(), 2 * kLocalHalf, rd.scratch);
        ok (gz.offset == 0 && std::fpclassify (gz.score) == FP_ZERO && gz.cells == 0, "fewer offsets than one neighbourhood: nothing is read");
    }

    group ("the reading is a function of the curve alone: the same curve twice, with a dirty scratch between");
    {
        const int offsets = 576;
        auto c = flatCurve (offsets, 4242u, 0.1);
        for (int k = 0; k < kCells; k += 2) c[(std::size_t) 123 * kCells + (std::size_t) k] -= 2.0f;
        Reader rd (offsets);
        const GridReading a = readCurve (c.data(), offsets, rd.scratch);
        for (auto& v : rd.doubles) v = 1.0e300;
        for (auto& v : rd.floats) v = -1.0e30f;
        const GridReading b = readCurve (c.data(), offsets, rd.scratch);
        ok (a.offset == b.offset && a.cells == b.cells && std::memcmp (&a.score, &b.score, sizeof (double)) == 0
                && std::memcmp (&a.second, &b.second, sizeof (double)) == 0 && std::memcmp (&a.localDipSum, &b.localDipSum, sizeof (double)) == 0,
            "bit-identical readings: nothing is carried in the scratch");
        ok (a.offset == 123 && a.cells == kCells / 2, "and the dip in every other cell is counted in exactly 16 cells");
    }

    return felitronics::test::report();
}

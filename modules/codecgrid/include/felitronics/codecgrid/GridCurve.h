// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/core/DetMath.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>

//==============================================================================
// felitronics::codecgrid — the EVIDENCE of a codec's frame grid, and how it is read.
//
// THE CURVE. A stretch of stereo PCM is analysed with one codec's transform at every frame offset it could have
// been coded at. For one offset, one signal (left, right, mid, side) and one band group, the number kept is
//
//     the mean over frames and coefficients of 20 log10 max (|X| / rms, 1e-4)
//
// where rms is that coefficient's own level over the stretch, taken once, at offset 0. A coefficient is read no
// deeper than 80 dB below its own level: a quantiser's zero is "far below where this coefficient usually is", and
// how far below stops mattering — the decoder's own rounding, the file's word length and the arithmetic here all
// live down there. At the encoder's offset the zeros pull the mean down; one sample off they are gone. The curve
// is those numbers for every offset: [offsets][4 signals][8 groups], 32 CELLS per offset.
//
// THE READING. A ratio is not evidence, so four numbers are read off a curve and all four are needed:
//   score        how far the deepest offset stands out — the dip summed over the cells, each cell measured against
//                the LOCAL median of the 16 neighbouring offsets, in robust sigmas of that sum over all offsets;
//   localDipSum  the same dip in decibels, not in sigmas: a pure tone gives a large score on a dip of nothing;
//   cells        in how many of the 32 cells the dip stands out on its own (over 4 of that cell's own sigmas) —
//                a real grid shows in most signals and bands, an accident in one;
//   second       the score of the best offset more than two samples away: a periodic signal has many such
//                offsets, a coded one has one.
// Every one of them was a plausible single number before it was shown not to be enough (a global median instead
// of the local one read 8 where the local one reads 186; a median wrapped around the ends of the offset range
// put every false maximum at offset 0 — so the neighbourhood is MIRRORED at the ends, never wrapped).
//
// THE PRODUCT INSTEAD OF A LOGARITHM PER COEFFICIENT. A mean of logarithms is the logarithm of a product, so the
// ratios are multiplied, sixteen frames at a time, and one logarithm is taken per sixteen — exactly the same sum,
// a sixteenth of the logarithms. Sixteen ratios in [1e-4, about 1e3] cannot leave a double.
//
// No allocation anywhere here: every function works in storage the caller owns.
//==============================================================================
namespace felitronics::codecgrid
{

inline constexpr int kSignals = 4;                 // left, right, mid, side
inline constexpr int kGroups  = 8;                 // band groups per signal
inline constexpr int kCells   = kSignals * kGroups;
inline constexpr int kZeroBands = 32;              // the zero share is published in 32 bands
inline constexpr double kLevelFloor = 1.0e-4;      // a coefficient is read no deeper than 80 dB below its level
inline constexpr double kZeroLevel  = 1.0e-2;      // "a zero": 40 dB below the coefficient's own level
inline constexpr int kLocalHalf = 8;               // the local median takes this many offsets each side
inline constexpr int kLogBatch  = 16;              // frames per logarithm

//==============================================================================
// LevelAccumulator — turns frames of coefficients (left and right, `bins` each) into the 32 cells of one offset.
// The caller owns the storage: inverseRef, product, logSum — kSignals * bins doubles each.
//==============================================================================
class LevelAccumulator
{
public:
    static constexpr std::size_t doublesFor (int bins) noexcept { return (std::size_t) kSignals * (std::size_t) bins; }

    void attach (double* inverseRef, double* product, double* logSum, int bins) noexcept
    {
        inverseRef_ = inverseRef; product_ = product; logSum_ = logSum; bins_ = bins;
    }

    int bins() const noexcept { return bins_; }

    //--- the reference level: the rms of every (signal, coefficient) over the frames of ONE offset
    void beginReference() noexcept
    {
        std::fill (logSum_, logSum_ + doublesFor (bins_), 0.0);
        frames_ = 0;
    }

    void addReferenceFrame (const float* left, const float* right) noexcept
    {
        for (int k = 0; k < bins_; ++k)
        {
            const double l = (double) left[k], r = (double) right[k];
            const double m = 0.5 * (l + r), s = 0.5 * (l - r);
            logSum_[(std::size_t) k]             += l * l;
            logSum_[(std::size_t) (bins_ + k)]     += r * r;
            logSum_[(std::size_t) (2 * bins_ + k)] += m * m;
            logSum_[(std::size_t) (3 * bins_ + k)] += s * s;
        }
        ++frames_;
    }

    // False when there was no frame: no reference, and nothing may be accumulated against it.
    bool finishReference() noexcept
    {
        if (frames_ <= 0) return false;
        const std::size_t n = doublesFor (bins_);
        for (std::size_t i = 0; i < n; ++i)
            inverseRef_[i] = 1.0 / (std::sqrt (logSum_[i] / (double) frames_) + 1.0e-30);
        return true;
    }

    //--- one offset
    void beginOffset() noexcept
    {
        const std::size_t n = doublesFor (bins_);
        std::fill (product_, product_ + n, 1.0);
        std::fill (logSum_, logSum_ + n, 0.0);
        frames_ = 0;
        pending_ = 0;
    }

    void addFrame (const float* left, const float* right) noexcept
    {
        const double* ir = inverseRef_;
        for (int k = 0; k < bins_; ++k)
        {
            const double l = (double) left[k], r = (double) right[k];
            const double v[kSignals] { l, r, 0.5 * (l + r), 0.5 * (l - r) };
            for (int s = 0; s < kSignals; ++s)
            {
                const std::size_t i = (std::size_t) (s * bins_ + k);
                double ratio = std::fabs (v[s]) * ir[i];
                if (! (ratio >= kLevelFloor)) ratio = kLevelFloor;     // a NaN is read as the floor, never kept
                if (ratio > 1.0e12) ratio = 1.0e12;                    // sixteen of these stay inside a double
                product_[i] *= ratio;
            }
        }
        ++frames_;
        if (++pending_ == kLogBatch) flush();
    }

    // cells: kCells floats, [signal][group]. The mean log level of the offset, in decibels. Zero frames: zeros.
    void finishOffset (float* cells) noexcept
    {
        if (pending_ > 0) flush();
        const int perGroup = bins_ / kGroups;
        for (int s = 0; s < kSignals; ++s)
            for (int g = 0; g < kGroups; ++g)
            {
                double acc = 0.0;
                const double* p = logSum_ + (std::size_t) (s * bins_ + g * perGroup);
                for (int k = 0; k < perGroup; ++k) acc += p[k];
                cells[s * kGroups + g] = frames_ > 0 ? (float) (20.0 * acc / ((double) perGroup * (double) frames_)) : 0.0f;
            }
    }

    int frames() const noexcept { return frames_; }

    //--- the zero share at one offset: per signal, how many coefficients of each of 32 bands sit 40 dB under
    // their own level. counts: kSignals * kZeroBands integers, zeroed by the caller before the first frame.
    void addZeroFrame (const float* left, const float* right, std::uint32_t* counts) const noexcept
    {
        const int perBand = bins_ / kZeroBands;
        for (int k = 0; k < perBand * kZeroBands; ++k)
        {
            const double l = (double) left[k], r = (double) right[k];
            const double v[kSignals] { l, r, 0.5 * (l + r), 0.5 * (l - r) };
            const int band = k / perBand;
            for (int s = 0; s < kSignals; ++s)
                if (std::fabs (v[s]) * inverseRef_[(std::size_t) (s * bins_ + k)] < kZeroLevel)
                    ++counts[s * kZeroBands + band];
        }
    }

private:
    void flush() noexcept
    {
        const std::size_t n = doublesFor (bins_);
        for (std::size_t i = 0; i < n; ++i)
        {
            logSum_[i] += core::det::log10 (product_[i]);
            product_[i] = 1.0;
        }
        pending_ = 0;
    }

    double* inverseRef_ = nullptr;
    double* product_ = nullptr;
    double* logSum_ = nullptr;
    int bins_ = 0, frames_ = 0, pending_ = 0;
};

//==============================================================================
// The reading of one curve.
//==============================================================================
struct GridReading
{
    int    offset = 0;            // the offset that stands out most, 0 .. offsets - 1
    double score = 0.0;           // in robust sigmas of the summed local dip
    double second = 0.0;          // the best offset more than two samples away, same units
    double localDipSum = 0.0;     // the dip at `offset`, dB, summed over the cells
    int    cells = 0;             // cells in which the dip stands out on its own
};

// The rule the measurement campaign fixed before its held-out set was opened: a grid is FOUND in this stretch.
struct GridRule
{
    double minScore = 9.2;        // 1.5 x the largest score of 168 lossless stretches
    double minDipDb = 5.0;
    int    minCells = 6;
    double maxSecondShare = 0.5;  // the runner-up scores at most this share of the winner

    constexpr bool broad (const GridReading& r) const noexcept
    {
        return r.score > minScore && r.localDipSum > minDipDb && r.cells >= minCells;
    }
    constexpr bool found (const GridReading& r) const noexcept { return broad (r) && r.second <= maxSecondShare * r.score; }
    // broad but not unique: several offsets stand out alike — spliced coded material, or a periodic signal
    constexpr bool several (const GridReading& r) const noexcept { return broad (r) && ! (r.second <= maxSecondShare * r.score); }
};

namespace detail
{
    // The median of n doubles; the mean of the two middle ones when n is even. Reorders `v`.
    inline double medianOf (double* v, int n) noexcept
    {
        if (n <= 0) return 0.0;
        std::sort (v, v + n);
        return (n & 1) != 0 ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
    }

    // The neighbour of offset i at distance d, mirrored at the ends without repeating the end sample.
    constexpr int mirrored (int i, int n) noexcept
    {
        if (n <= 1) return 0;
        const int period = 2 * (n - 1);
        int m = i % period;
        if (m < 0) m += period;
        return m < n ? m : period - m;
    }
}

// What readCurve() needs beside the curve: `offsets` doubles for the summed dip, `offsets` for a scratch, and
// offsets * kCells floats for the dip of every cell.
struct ReadScratch
{
    double* local = nullptr;
    double* sort = nullptr;
    float*  dip = nullptr;
    static constexpr std::size_t doublesFor (int offsets) noexcept { return 2u * (std::size_t) offsets; }
    static constexpr std::size_t floatsFor (int offsets) noexcept { return (std::size_t) offsets * kCells; }
};

// curve: offsets * kCells floats, [offset][cell]. offsets >= 2 * kLocalHalf + 1.
inline GridReading readCurve (const float* curve, int offsets, const ReadScratch& s) noexcept
{
    GridReading out;
    if (offsets < 2 * kLocalHalf + 1) return out;

    for (int c = 0; c < kCells; ++c)
        for (int o = 0; o < offsets; ++o)
        {
            float w[2 * kLocalHalf + 1];
            for (int d = -kLocalHalf; d <= kLocalHalf; ++d)
                w[d + kLocalHalf] = curve[(std::size_t) detail::mirrored (o + d, offsets) * kCells + (std::size_t) c];
            std::sort (w, w + 2 * kLocalHalf + 1);
            s.dip[(std::size_t) o * kCells + (std::size_t) c] = w[kLocalHalf] - curve[(std::size_t) o * kCells + (std::size_t) c];
        }

    for (int o = 0; o < offsets; ++o)
    {
        double acc = 0.0;
        for (int c = 0; c < kCells; ++c) acc += (double) s.dip[(std::size_t) o * kCells + (std::size_t) c];
        s.local[o] = acc;
    }

    std::copy (s.local, s.local + offsets, s.sort);
    const double med = detail::medianOf (s.sort, offsets);
    for (int o = 0; o < offsets; ++o) s.sort[o] = std::fabs (s.local[o] - med);
    const double sigma = 1.4826 * detail::medianOf (s.sort, offsets) + 1.0e-9;

    int best = 0;
    for (int o = 1; o < offsets; ++o) if (s.local[o] > s.local[best]) best = o;

    bool any = false;
    double second = 0.0;
    for (int o = 0; o < offsets; ++o)
    {
        int d = o - best;
        if (d < 0) d = -d;
        if (offsets - d < d) d = offsets - d;                  // the offsets are a circle: 0 follows the last one
        if (d <= 2) continue;
        if (! any || s.local[o] > second) { second = s.local[o]; any = true; }
    }

    int cells = 0;
    for (int c = 0; c < kCells; ++c)
    {
        for (int o = 0; o < offsets; ++o) s.sort[o] = (double) s.dip[(std::size_t) o * kCells + (std::size_t) c];
        const double cellMed = detail::medianOf (s.sort, offsets);
        for (int o = 0; o < offsets; ++o) s.sort[o] = std::fabs ((double) s.dip[(std::size_t) o * kCells + (std::size_t) c] - cellMed);
        const double cellSigma = 1.4826 * detail::medianOf (s.sort, offsets) + 1.0e-6;
        if ((double) s.dip[(std::size_t) best * kCells + (std::size_t) c] / cellSigma > 4.0) ++cells;
    }

    out.offset = best;
    out.score = s.local[best] / sigma;
    out.second = any ? second / sigma : 0.0;
    out.localDipSum = s.local[best];
    out.cells = cells;
    return out;
}

// Two grid phases agree when they are within `tolerance` samples of each other on the circle of `hop`.
constexpr bool phasesAgree (int a, int b, int hop, int tolerance) noexcept
{
    int d = (a - b) % hop;
    if (d < 0) d += hop;
    if (hop - d < d) d = hop - d;
    return d <= tolerance;
}

} // namespace felitronics::codecgrid

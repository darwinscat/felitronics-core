// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <span>
#include <vector>

namespace felitronics::analysis
{

//==============================================================================
// The five sums one stretch of stereo reduces to, and the three numbers read out of them. Shared by the column
// band and the playhead needle below, which are the SAME formula over different stretches — a second copy of it
// is how the two would start to disagree.
//
// THE CORRELATION IS NOT PEARSON'S, whatever the spec's comment calls it: nothing is centred. It is the
// normalised dot product ΣLR / √(ΣLL·ΣRR) — the classic phase-correlation meter's number. On L = (1, 2),
// R = (2, 1) it reads 0.8, where Pearson's coefficient reads −1. That is the quantity the page draws, so it is
// the quantity here; an oracle written from the word "Pearson" would reject a correct port.
//
// BINARY64, AND THE FIVE PRODUCTS ARE STORED BEFORE THEY ARE ADDED. The spec is JavaScript, which never fuses
// `s += a*b` into one FMA; a C++ build may (clang's default `-ffp-contract=on` fuses the single expression on
// arm64, gcc's `fast` fuses across statements too), and a fused sum is a different number — one rounding where
// the spec has two — so a column and a verdict threshold drift by an ulp. Each product goes through a `volatile`
// store, which no contraction mode can see through (law 10; `analysis::ClipDetector` pins its variance the same
// way). Everything is promoted to double before it is touched, exactly as the spec's Float32Array reads are.
struct StereoSums
{
    double ll = 0.0, rr = 0.0, lr = 0.0, mid = 0.0, side = 0.0;

    inline void add (float lf, float rf) noexcept
    {
        const double l = (double) lf, r = (double) rf;
        volatile double p;
        p = l * l; ll += p;
        p = r * r; rr += p;
        p = l * r; lr += p;
        const double m = (l + r) * 0.5, s = (l - r) * 0.5;
        p = m * m; mid  += p;
        p = s * s; side += p;
    }

    // corr = ΣLR / √(ΣLL·ΣRR), clamped to [−1, 1]; a denominator not above 1e-12 reads +1 (silence is neutral).
    // The clamp is JavaScript's Math.max(−1, Math.min(1, x)), which lets a NaN through — std::clamp would too,
    // but std::min / std::max would not, so it is spelled out.
    double correlation() const noexcept
    {
        const double denom = std::sqrt (ll * rr);
        if (! (denom > 1e-12)) return 1.0;
        double x = lr / denom;
        if (x > 1.0)  x = 1.0;
        if (x < -1.0) x = -1.0;
        return x;
    }

    // width = side_rms / (mid_rms + side_rms), 0 when the sum is not above 1e-9. `n` is the stretch's frame count
    // floored at 1, as the spec floors it.
    double width (std::uint64_t n) const noexcept
    {
        const double dn = (double) std::max<std::uint64_t> (1, n);
        const double midRms = std::sqrt (mid / dn), sideRms = std::sqrt (side / dn);
        return (midRms + sideRms) > 1e-9 ? sideRms / (midRms + sideRms) : 0.0;
    }

    // rms = √((ΣLL + ΣRR) / (2n)) — the stretch's RMS over both channels. The spec calls it `loud`; it is NOT a
    // loudness in the BS.1770 sense (no K-weighting, no gate, no LUFS), and a name that sat next to `fc_probe_lufs`
    // in the same ABI as "loud" would be two numbers under one label — so here it is called what it is.
    double rms (std::uint64_t n) const noexcept
    {
        const double dn = (double) std::max<std::uint64_t> (1, n);
        return std::sqrt ((ll + rr) / (2.0 * dn));
    }
};

//==============================================================================
// felitronics::analysis::StereoColumns — the stereo band: a whole file reduced to `columns` stretches, each with
// its width, its phase correlation and its RMS (the spec's `loud`). Plus `needle()`, the same three over any stretch — the
// playhead readout.
//
// THIS IS A PORT, AND THE SPEC IS EXECUTABLE: the site's `computeStereoColumns`, `correlationOf` and `widthOf`
// (stereo-meter.js). Bit for bit, on the same float32 PCM, and a cross-language NULL says so. The VERDICT the page
// draws from these numbers (`evaluateStereo`, `widthCode`, `verdictForCorr` and their thresholds) is product
// policy and stays on the page; this header is the measurement only.
//
// THE DEFINITION (the spec's names):
//   cols  = max(1, min(columns, frames))        a short file gets FEWER columns, not zero-filled ones. `columns < 1`
//                                               is REFUSED here, where the spec clamps it to 1: a column count of
//                                               zero or less is a caller's mistake, not a picture
//   per   = frames / cols                       a double
//   column i covers [floor(i · per), min(frames, floor((i + 1) · per)))
//   L = channel 0; R = channel 1, or channel 0 again on a mono file (width 0, correlation +1); channels past 1
//   are not read.
// OUTPUT TYPES FOLLOW THE SPEC, and they are not all one type. `width`, `corr` and `loud` (here `rms`) are
// Float32Array there, so they are float32 here. `maxLoud` (here `maxRms`) is the maximum of the UNROUNDED double
// RMS — every threshold of the page's verdict compares a column against `maxLoud · 0.15`, so recomputing it from the
// rounded array moves it by an ulp and flips a column sitting on the threshold. And the needle returns doubles, not
// the column's float32.
//
// THE SPEC'S KNOWN PROPERTIES ARE KEPT, each with a pinned witness:
//   * A COLUMN BOUNDARY IS floor(i · (frames / cols)), not the exact floor(i · frames / cols): at frames = 1206,
//     cols = 1200 column 200 starts at frame 200, where the exact partition says 201.
//   * Frames past the last column's end are read by no column (floor(cols · per) can land below `frames`).
//   * A NaN sample makes its column's sums NaN: width 0 (the comparison fails), correlation +1, RMS NaN, and it
//     cannot become `maxRms`. An infinite one makes the RMS +Inf and the rest depends on its partner: L = +Inf against
//     R = 1 reads correlation NaN and width NaN (the clamp lets NaN through); against R = 0, ΣLL·ΣRR is Inf·0 = NaN and
//     the correlation guard reads +1; L = R = +Inf reads width 0 (the side sum is Inf − Inf).
//
// WHAT PARITY DOES NOT COVER, because it is runtime state and not code (law 10): FTZ/DAZ. The spec measures a
// subnormal pair at an RMS of 1e-40; a host thread with denormals flushed reads 0, and wasm cannot flush at
// all. This class never sets either — run it with them off (an offline tool is) when the bits must match.
//
// STREAMING, like WaveformPeaks: the total length is a `prepare()` argument, any split of `process()` gives the
// bits one call gives, a call past the prepared length is refused whole. The width is at least 1 and EXACT.
// RT-safe: prepare() is the only allocation.
class StereoColumns
{
public:
    static constexpr int           kDefaultColumns = 1200;              // the spec's default; the page asks for 1100
    static constexpr int           kMaxColumns     = 1 << 20;            // capacity guard, not part of the picture
    static constexpr std::uint64_t kMaxFrames      = std::uint64_t (1) << 53;

    // The three numbers over one stretch [from, to) of two planes — `correlationOf` and `widthOf`, plus the
    // stretch's RMS. `right` may be `left` (a mono file). REFUSES a stretch that is not inside `frames` —
    // and that is a difference from the spec on purpose: JavaScript reads past an array as `undefined`, the sums
    // turn NaN, and the guards answer correlation +1 / width 0, i.e. "mono, in phase", for a stretch that does
    // not exist. A clamp to the planes would be worse still (on [0.5, 0.5] against [-0.5, -0.5] it reads -1 / 1
    // where the spec reads 1 / 0). A caller asked a question with no answer, and is told so.
    struct Needle { double correlation = 1.0, width = 0.0, rms = 0.0; };
    [[nodiscard]] static bool needle (const float* left, const float* right, std::uint64_t frames,
                                      std::uint64_t from, std::uint64_t to, Needle& out) noexcept
    {
        if (left == nullptr || right == nullptr || from > to || to > frames) return false;
        StereoSums s;
        for (std::uint64_t j = from; j < to; ++j) s.add (left[j], right[j]);
        const std::uint64_t n = to - from;
        out.correlation = s.correlation();
        out.width       = s.width (n);
        out.rms    = s.rms (n);
        return true;
    }

    [[nodiscard]] bool prepare (int numChannels, std::uint64_t totalFrames, int columns = kDefaultColumns)
    {
        prepared_ = false;
        if (numChannels < 1) return false;
        if (totalFrames < 1 || totalFrames > kMaxFrames) return false;
        if (columns < 1 || columns > kMaxColumns) return false;

        nch_    = numChannels;
        frames_ = totalFrames;
        cols_   = (int) std::min<std::uint64_t> ((std::uint64_t) columns, totalFrames);   // >= 1: both are
        per_    = (double) totalFrames / (double) cols_;
        width_.assign ((std::size_t) cols_, 0.0f);
        corr_.assign  ((std::size_t) cols_, 0.0f);
        rms_.assign  ((std::size_t) cols_, 0.0f);
        prepared_ = true;
        reset();
        return true;
    }

    void reset() noexcept
    {
        std::fill (width_.begin(), width_.end(), 0.0f);
        std::fill (corr_.begin(),  corr_.end(),  0.0f);
        std::fill (rms_.begin(),  rms_.end(),  0.0f);
        maxRms_ = 0.0; seen_ = 0; sums_ = {};
        col_ = 0; enterColumn();
    }

    [[nodiscard]] bool process (const float* const* channels, int numChannels, int n) noexcept
    {
        if (! prepared_ || channels == nullptr || numChannels != nch_ || n < 0) return false;
        if ((std::uint64_t) n > frames_ - seen_) return false;
        if (n > 0 && (channels[0] == nullptr || (nch_ > 1 && channels[1] == nullptr))) return false;

        const float* L = channels[0];
        const float* R = nch_ > 1 ? channels[1] : channels[0];
        for (int k = 0; k < n; ++k)
        {
            const std::uint64_t j = seen_ + (std::uint64_t) k;
            while (col_ < cols_ && j >= end_) closeColumn();
            if (col_ < cols_ && j >= start_) sums_.add (L[k], R[k]);
        }
        seen_ += (std::uint64_t) n;
        if (seen_ == frames_) while (col_ < cols_) closeColumn();
        return true;
    }

    [[nodiscard]] bool prepared() const noexcept { return prepared_; }
    [[nodiscard]] bool complete() const noexcept { return prepared_ && seen_ == frames_; }

    int    columns() const noexcept { return cols_; }
    bool   isMono()  const noexcept { return nch_ < 2; }
    double maxRms() const noexcept { return maxRms_; }
    std::span<const float> width()       const noexcept { return { width_.data(), width_.size() }; }
    std::span<const float> correlation() const noexcept { return { corr_.data(),  corr_.size() }; }
    std::span<const float> rms()    const noexcept { return { rms_.data(),  rms_.size() }; }
    int           numChannels() const noexcept { return nch_; }
    std::uint64_t totalFrames() const noexcept { return frames_; }
    std::uint64_t framesSeen()  const noexcept { return seen_; }

    // Where column i starts and ends — the spec's own expressions, exposed so a caller can map a column to time.
    std::uint64_t columnStart (int i) const noexcept { return (std::uint64_t) std::floor ((double) i * per_); }
    std::uint64_t columnEnd   (int i) const noexcept
    {
        return std::min (frames_, (std::uint64_t) std::floor ((double) (i + 1) * per_));
    }

private:
    void enterColumn() noexcept
    {
        if (col_ >= cols_) return;
        start_ = columnStart (col_);
        end_   = columnEnd (col_);
        sums_  = {};
    }

    void closeColumn() noexcept
    {
        // `end - start` floored at 1, as the spec's `Math.max(1, end - start)`; an empty column is legal.
        const std::uint64_t n = end_ > start_ ? end_ - start_ : 0;
        const auto i = (std::size_t) col_;
        width_[i] = (float) sums_.width (n);
        corr_[i]  = (float) sums_.correlation();
        const double rms = sums_.rms (n);
        rms_[i]  = (float) rms;
        if (rms > maxRms_) maxRms_ = rms;
        ++col_;
        enterColumn();
    }

    std::vector<float> width_, corr_, rms_;
    StereoSums    sums_;
    double        per_ = 0.0, maxRms_ = 0.0;
    std::uint64_t frames_ = 0, seen_ = 0, start_ = 0, end_ = 0;
    int           nch_ = 0, cols_ = 0, col_ = 0;
    bool          prepared_ = false;
};

} // namespace felitronics::analysis

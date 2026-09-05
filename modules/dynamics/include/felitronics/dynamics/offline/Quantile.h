// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <cmath>
#include <cstdint>
#include <vector>

namespace felitronics::dynamics::offline
{

//==============================================================================
// felitronics::dynamics::offline::QuantileHistogram — ONE definition of "the p95", written once.
//
// It exists because a percentile defined in two places is a way to pass a round-trip test that
// measures nothing: the solver picks a threshold by a quantile and the acceptance then checks it with
// a quantile, so if the two agree only with each other, a wrong shared rule is invisible. Everything
// in P3 that says "p95" comes through here, and the tests NULL it against an independently written
// exact sort.
//
// THE RULE, named exactly, because "the 95th percentile" is not one thing: NEAREST-RANK. With N values
// added, the answer is the order statistic of rank `ceil(q·N)` (1-based, clamped to [1, N]) — a value
// that actually occurred, not a point invented between two that did. The histogram reports the CENTRE
// of the bin that rank falls in, so the error is at most half a bin and is unbiased, and the result is
// then CLAMPED to the observed [min, max].
//
// It is deliberately NOT an interpolated rule. Interpolating across a bin assumes the values are
// spread through it, which is a second definition needing its own justification and its own tests, and
// which invents a spread exactly where there is none: an ATOM. Measured, with 0.05 dB bins, on a trace
// that is a constant 2.000 dB: inverse-CDF interpolation reports 2.0475, and mixing in zeros moves it
// only between 2.0450 and 2.0475 rather than away — the error is the bin, not the mixture. Atoms are
// not an exotic case here: a steady tone produces a constant gain reduction, which is what a synthetic
// acceptance fixture is made of.
//
// WHAT IS EXACT, precisely — because the useful version of this promise is the narrow one. The
// observed minimum and maximum are TRACKED rather than binned, and so are the counts of values equal
// to each. So the answer is returned exactly whenever the rank lands inside the atom at either end:
// an all-zero trace, a constant trace, and — the case a synthetic fixture actually produces — a trace
// that is silent 96 % of the time and 2 dB the rest, whose p95 is exactly zero. A binned rule without
// those counters reports 0.025 there, which is the sort of number that makes an analytic fixture look
// approximately right instead of wrong. Everything strictly between the two atoms carries the
// half-bin bound and nothing better.
//
// IT IS ALSO NOT THE RULE `analysis::LoudnessMeter` USES FOR LRA. That one takes the bin's LOWER bound
// with a rank of `(N-1)·p + 0.5`, which is faithful to libebur128 and is correct there because LRA is a
// DIFFERENCE of two percentiles and the half-bin offset cancels. Nothing cancels here, so the offset
// would be a bias. Two rules, both named, neither pretending to be the other.
//
// COUNTERS ARE 64-BIT. A uint32 bin overflows after 2^32 samples: 24.9 hours at 48 kHz, 6.2 at 192.
// That is not reachable in one offline pass today and is entirely reachable by a meter left running.
//
// OUT OF RANGE IS REPORTED, NEVER SILENTLY SATURATED. Values above the configured top land in an
// overflow counter, and a quantile that falls among them returns false rather than the top of the
// range — a compressor may legitimately be asked for a 400 dB range (GainComputer's own cap), and a
// solver told "the p95 is 96" when it is 150 would converge, confidently, to the wrong threshold.
//
// OFFLINE / MESSAGE THREAD: prepare() allocates. add() does not, and is one branch and one increment.
class QuantileHistogram
{
public:
    // `binDb` is the resolution the quantile is reported to: the error of the rule above is bounded by
    // one bin, and is zero wherever the clamp applies. Silently reasonable defaults are avoided here —
    // the caller states the range it can produce, because only the caller knows it.
    // REFUSES rather than silently narrowing. A bin width fine enough to need more bins than the cap
    // used to be honoured by keeping the WIDTH and dropping the RANGE, so a caller asking for 1e-6 dB
    // resolution over 360 dB got a histogram covering four of them — and then a valid-looking p95 of
    // zero for a signal at -20 dB. A measurement that cannot be made must fail where it is configured.
    [[nodiscard]] bool prepare (double loDb, double hiDb, double binDb)
    {
        bins_.clear();
        if (! std::isfinite (loDb) || ! std::isfinite (hiDb) || ! (binDb > 0.0) || ! std::isfinite (binDb)) return false;
        if (! (hiDb > loDb)) return false;
        lo_  = loDb;
        bin_ = binDb;
        // `floor(span/bin) + 1` bins, so `hiDb` itself lands in the LAST bin rather than one past the
        // end, and nothing above it does. It is not a rounding nicety: |gain change| is bounded BY
        // `rangeDb` and reaches it exactly on saturated material, so a half-open top sent the saturated
        // case — the one a solver meets when asked for more than the material can give — into overflow.
        const double nb = std::floor ((hiDb - lo_) / bin_) + 1.0;
        if (! (nb >= 1.0) || nb > 4.0e6) return false;    // 32 MB ceiling; absurd resolutions are refused
        bins_.assign ((std::size_t) nb, 0);
        reset();
        return true;
    }

    void reset() noexcept
    {
        for (auto& c : bins_) c = 0;
        n_ = 0; over_ = 0; under_ = 0; nonFinite_ = 0; sum_ = 0.0;
        min_ = 0.0; max_ = 0.0; minCount_ = 0; maxCount_ = 0; seen_ = false;
    }

    // One value. Non-finite values are COUNTED, not silently dropped and not folded into a bin: a NaN
    // in a gain-reduction trace means the caller's chain is poisoned, and turning that into a plausible
    // number is how a broken measurement looks healthy.
    inline void add (double v) noexcept
    {
        if (! std::isfinite (v)) { ++nonFinite_; return; }
        if (! seen_) { min_ = max_ = v; minCount_ = maxCount_ = 1; seen_ = true; }
        else
        {
            // The equality is deliberate — it counts an ATOM, and an atom is exactly a set of values
            // that compare equal. Spelled through `<` so that `-Wfloat-equal` keeps meaning something
            // everywhere else, where a float equality usually IS a mistake.
            if      (v < min_)      { min_ = v; minCount_ = 1; }
            else if (! (min_ < v))  ++minCount_;
            if      (v > max_)      { max_ = v; maxCount_ = 1; }
            else if (! (v < max_))  ++maxCount_;
        }
        ++n_;
        sum_ += v;
        const double idx = std::floor ((v - lo_) / bin_);
        if (idx < 0.0)                       { ++under_; return; }
        if (idx >= (double) bins_.size())    { ++over_;  return; }
        ++bins_[(std::size_t) idx];
    }

    std::uint64_t count()         const noexcept { return n_; }
    std::uint64_t aboveRange()    const noexcept { return over_; }
    std::uint64_t belowRange()    const noexcept { return under_; }
    std::uint64_t nonFiniteCount()const noexcept { return nonFinite_; }
    double        minValue()      const noexcept { return seen_ ? min_ : 0.0; }
    double        maxValue()      const noexcept { return seen_ ? max_ : 0.0; }
    double        binWidth()      const noexcept { return bin_; }

    // The q-quantile by the rule above. Returns false — leaving `out` untouched — when there is
    // nothing to report (no values) or when the answer lies outside the configured range, which is a
    // different thing from an extreme answer and must not be confused with one.
    [[nodiscard]] bool quantile (double q, double& out) const noexcept
    {
        if (n_ == 0 || ! (q >= 0.0) || ! (q <= 1.0)) return false;
        // Nearest rank, 1-based. `ceil` and not rounding: the q-quantile is the smallest value the
        // fraction q of the sample does not exceed, and for q = 1 that has to be the maximum.
        double rank = std::ceil (q * (double) n_);
        if (rank < 1.0) rank = 1.0;
        if (rank > (double) n_) rank = (double) n_;

        // THE TWO END ATOMS ARE TRACKED, not binned, so any rank landing inside either is exact — and
        // that is the case a synthetic fixture keeps producing. Guarded by the range counters, because
        // a minimum that fell BELOW the configured range is still out of range and must be refused
        // rather than reported: `out of range returns false` is the contract, and an extreme rank is
        // exactly where it would otherwise be quietly bypassed.
        if (rank <= (double) minCount_ && under_ == 0 && inRange (min_))      { out = min_; return true; }
        if (rank > (double) (n_ - maxCount_) && over_ == 0 && inRange (max_)) { out = max_; return true; }

        double cum = (double) under_;
        if (rank <= cum) return false;                     // the answer is below the configured range
        for (std::size_t b = 0; b < bins_.size(); ++b)
        {
            const double c = (double) bins_[b];
            if (c <= 0.0) continue;
            cum += c;
            if (cum >= rank)
            {
                double v = lo_ + ((double) b + 0.5) * bin_;      // the bin's CENTRE
                if (v < min_) v = min_;
                if (v > max_) v = max_;
                out = v;
                return true;
            }
        }
        return false;                                      // it fell among the overflow counts
    }

    // The arithmetic mean and the maximum, both exact — they need no binning, and the solver can
    // target them for the same reason it can target a quantile: every one of them is monotone under
    // pointwise domination of the trace, which is what makes the search valid.
    double mean() const noexcept { return n_ > 0 ? sum_ / (double) n_ : 0.0; }

    // The fraction of values strictly above `vDb`, to within one bin. Reported as a fraction of the
    // values ADDED, so the out-of-range counts are included in the denominator and in the numerator
    // where they belong — an overflowing value is above any in-range threshold.
    double fractionAbove (double vDb) const noexcept
    {
        if (n_ == 0) return 0.0;
        double idx = std::floor ((vDb - lo_) / bin_);
        double above = (double) over_;
        for (std::size_t b = 0; b < bins_.size(); ++b)
            if ((double) b > idx) above += (double) bins_[b];
        return above / (double) n_;
    }

private:
public:
    // THE COUNTER TYPE, named once so it can be pinned once. A 32-bit bin overflows after 2^32 samples
    // — 24.9 hours at 48 kHz, 6.2 at 192 — which no suite can reach by counting, so the contract is
    // asserted where it is decidable rather than left to a comment nobody can fail.
    using BinCount = std::uint64_t;

private:
    // Whether a value landed in a bin at all. The extreme-atom shortcuts above need it: with every
    // value ABOVE the configured top, `under_` is still zero and the minimum is still tracked, so a
    // guard on the underflow counter alone would have handed back an out-of-range minimum as an
    // in-range answer — which is the contract this class exists to keep.
    bool inRange (double v) const noexcept
    {
        const double idx = std::floor ((v - lo_) / bin_);
        return idx >= 0.0 && idx < (double) bins_.size();
    }

    std::vector<BinCount> bins_;
    std::uint64_t n_ = 0, over_ = 0, under_ = 0, nonFinite_ = 0, minCount_ = 0, maxCount_ = 0;
    double lo_ = 0.0, bin_ = 0.01, min_ = 0.0, max_ = 0.0, sum_ = 0.0;
    bool   seen_ = false;
};

} // namespace felitronics::dynamics::offline

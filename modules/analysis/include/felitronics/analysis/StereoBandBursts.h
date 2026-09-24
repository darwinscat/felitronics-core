// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.
#pragma once

#include <felitronics/analysis/BandBursts.h>
#include <felitronics/core/DetMath.h>
#include <felitronics/stereo/MidSide.h>

#include <cmath>
#include <cstdint>
#include <vector>

namespace felitronics::analysis
{

//======================================================================================================
// THE SAME BAND-BURST DETECTOR ON TWO AXES — Mid and Side — so that an event can be placed as CENTRED or
// WIDE. Two `BandBursts` engines prepared mono, fed in per-sample lockstep from one stereo input.
//
// WHAT THIS ANSWERS, AND WHAT IT DOES NOT. It separates centred from wide. It does NOT separate a hi-hat
// from a sibilant, which is the question the request opened with: a centred hat, shaker or snare top is
// centred exactly as an "s" is, and a doubled or reverbed vocal is wide although it is a voice. Position
// is real evidence and worth having, but the hat-versus-sibilant discriminator is PERIODICITY — a hat sits
// on the grid and an "s" does not — and each axis publishes the inter-onset histograms that carry it. A
// rule shaped "centred AND aperiodic AND a few hops long" uses both; one shaped "centred AND loud in band"
// uses neither.
//
// WHY THE "NO BURST ON THE OTHER AXIS" TEST CANNOT WORK, stated because it is the obvious rule to reach
// for. If one channel is a scaled copy of the other — any channel imbalance at all — then `R = k·L` makes
// `M = (1+k)/2·L` and `S = (1-k)/2·L`, so BOTH axes are scaled copies of the same signal. This detector
// fires on a RELATIVE excess over a trailing baseline, and a scaled copy has a proportionally scaled
// baseline, so the relative excess is IDENTICAL on the two axes. Measured: a 0.18 dB imbalance puts a Side
// event under every Mid event with the same 15.33 dB excess while Side sits 39.9 dB down. So a caller must
// compare the axes' POWERS, not their event lists — which is what `crossPower()` below is for.
//
// ONE PARAMETER SET, TWO ENGINES, ONE CLOCK. The hop schedule is `llround(fs·hopMs/1000)` and depends on
// nothing but the rate and the parameters, so both engines partition the same time into the same hops;
// feeding them per sample makes that an identity rather than an argument, and it is what lets a peak on
// one axis be joined to the other axis's hop at that instant.
//
// ABSENCE IS STRUCTURAL, NOT A THRESHOLD. Side is absent when the programme has ONE channel — and only
// then. Bit-identical channels give an exactly zero Side, which the engine already reports as
// `eventsValid() == false` with `NoEligibleHop` and every baseline zero: that is "explicitly absent"
// without anyone inventing a dual-mono predicate. A near-mono programme has a real Side carrying a real
// residue, and `exactZeroHops()` beside the band energies lets a caller tell the two apart and set its own
// threshold. An exact-zero test here would be a cliff and would also be redundant.
class StereoBandBursts
{
public:
    static constexpr int kMid = 0, kSide = 1, kAxes = 2;

    using Params = BandBurstsParams;

    // The other axis's hop at THIS event's peak. `eligible` is load-bearing: without it, "the other axis
    // was quiet" and "the other axis was not judged yet" are the same reading, and the second is what a
    // programme's first two seconds always look like.
    struct Cross
    {
        // MEAN SQUARES, DIVIDED BY THE HOP — the same units `BandBurst::peakPower` and `peakBaseline` are
        // in, and that identity is the whole point of these fields: a caller is meant to divide one by the
        // other. The hop trace publishes SUMS (`peakPower = peakNum_ / hop` at the event, `energy` raw at
        // the trace), so storing the trace's number straight would have put the two sides of that division
        // in different units and the ratio would have been wrong by `hopSamples` — 480 at 48 kHz, and a
        // DIFFERENT number at another rate, which is the shape of error that survives a fixture.
        //
        // Caught by arithmetic, not by a test: for `R = 0.98·L` the axes are scaled copies, so Mid must sit
        // (0.99/0.01)^2 = +39.9 dB over Side, and the wrapper reported +13.1 — short by exactly 10log10(480).
        double       power    = 0.0;    // the other axis's band mean square in the peak hop
        double       baseline = 0.0;    // …and its baseline there, 0 when it was not eligible
        bool         eligible = false;
        std::int64_t hop      = -1;     // the hop, as a start sample; -1 when nothing was captured
    };

    struct Storage
    {
        bool        ok      = false;
        std::size_t crosses = 0;                  // maxEvents per axis
        BandBursts::Storage axis {};              // ONE engine's demand; the total carries two

        std::uint64_t bytes() const noexcept
        {
            // ELEMENTS ONLY, never a container's own footprint. The two engines are MEMBERS of this
            // object, not heap allocations, so `sizeof (BandBursts)` is already inside
            // `sizeof (StereoBandBursts)` and counting it here both double-counts and — because a
            // std::vector holds POINTERS — publishes a different number on wasm32 than on a 64-bit
            // host. The cross-tier storage gate caught exactly that: 3 697 472 against 3 697 328,
            // 144 bytes, which is two engines' worth of three pointers each. A demand a page sizes a
            // heap from cannot depend on which tier answered.
            return 2u * axis.bytes()
                 + (std::uint64_t) crosses * kAxes * sizeof (Cross);
        }
    };

    static Storage storageFor (double sampleRate, const Params& p) noexcept
    {
        Storage st;
        // ONE CHANNEL EACH, always: the engines are fed Mid and Side, never L and R.
        const BandBursts::Storage one = BandBursts::storageFor (sampleRate, 1, p);
        if (! one.ok) return st;
        st.ok      = true;
        st.axis    = one;
        st.crosses = (std::size_t) (p.maxEvents > 0 ? p.maxEvents : 0);
        return st;
    }

    void setParams (const Params& p) noexcept { params_ = p; }
    const Params& params() const noexcept { return params_; }

    // Law 11b: disarmed first, so a refused prepare cannot be read as a measurement.
    [[nodiscard]] bool prepare (double sampleRate, int maxChannels) noexcept
    {
        prepared_ = false;
        finished_ = false;
        if (maxChannels < 1 || maxChannels > core::kMaxChannels) return false;
        const Storage st = storageFor (sampleRate, params_);
        if (! st.ok) return false;

        sampleRate_  = sampleRate;
        maxChannels_ = maxChannels;
        for (int a = 0; a < kAxes; ++a)
        {
            eng_[a].setParams (params_);
            if (! eng_[a].prepare (sampleRate, 1, 1)) return false;   // maxBlock sizes nothing here; one channel each
        }
        for (int a = 0; a < kAxes; ++a) cross_[a].assign (st.crosses, Cross {});

        // THE DOME, IN CLOSED FORM. The band is an LR4 pair, so for a pure tone the band-over-wideband
        // ratio peaks strictly below 1 — and where it peaks is not the band's geometric mean but the
        // geometric mean of the PREWARPED corners. Writing a = tan(pi·low/fs), b = tan(pi·high/fs), the
        // maximum of the product response sits at tan(pi·f/fs) = sqrt(a·b) and equals b^8/(a^2+b^2)^4.
        // Derived, then checked against a 200 001-point scan of the response at four rates: it agrees to
        // four decimals at every one. -3.9885 dB at 48 kHz, -3.8559 at 44.1, -4.5072 at 96.
        //
        // IT IS PUBLISHED BECAUSE THE SHARE ALONE IS NOT COMPARABLE. A caller thresholding "the band holds
        // most of this hop" against a textbook number will be wrong by the dome and wrong again at another
        // rate; against the dome the ratio is a fraction of what the instrument can even report.
        {
            const double a = core::det::tan (core::kPi * params_.bandLowHz  / sampleRate);
            const double b = core::det::tan (core::kPi * params_.bandHighHz / sampleRate);
            const double d = a * a + b * b;
            domeShare_ = (a > 0.0 && b > 0.0 && d > 0.0) ? (b * b * b * b * b * b * b * b) / (d * d * d * d) : 0.0;
            const double t = (a > 0.0 && b > 0.0) ? std::sqrt (a * b) : 0.0;
            // INVERTED BY BISECTION ON det::tan, NOT BY atan. The obvious spelling is
            // `sampleRate * atan(t) / pi`, and it is the one this shipped with until the libm audit
            // refused it: atan is a system transcendental and is not the same function on every row, so a
            // published frequency would differ between macOS, glibc and wasm — and this number crosses the
            // ABI. `sqrt` is fine and stays: IEEE-754 requires it correctly rounded, so it IS the same
            // everywhere. det::tan is strictly increasing on (0, pi/2) and t lies between a and b by
            // construction, so the band's own corners bracket the answer; 60 halvings take the interval
            // below 1e-14 Hz at any rate this library accepts. Comparisons and halves only — nothing here
            // reaches libm. Same idiom as LowEnd::settlingBlocks.
            if (t > 0.0)
            {
                double lo = params_.bandLowHz, hi = params_.bandHighHz;
                for (int it = 0; it < 60; ++it)
                {
                    const double mid = 0.5 * (lo + hi);
                    if (core::det::tan (core::kPi * mid / sampleRate) < t) lo = mid; else hi = mid;
                }
                domeHz_ = 0.5 * (lo + hi);
            }
            else domeHz_ = 0.0;
            if (! std::isfinite (domeShare_)) domeShare_ = 0.0;
            if (! std::isfinite (domeHz_))    domeHz_    = 0.0;
        }

        prepared_ = true;
        reset();
        return true;
    }

    void reset() noexcept
    {
        finished_ = false;
        ranNc_    = 0;
        samples_  = 0;
        for (int a = 0; a < kAxes; ++a)
        {
            eng_[a].reset();
            eng_[a].setObserver (&trampoline, this);
            std::fill (cross_[a].begin(), cross_[a].end(), Cross {});
            last_[a] = Cross {};
            exactZero_[a] = 0;
            pendingHop_[a] = -1;
            pending_[a] = Cross {};
            // AND THE COMMIT CURSOR, which this loop forgot and which is the whole of the defect it had.
            // `committed_` counts how many events' cross readings have been written; after a run it equals
            // the event count, and `while (committed_ < stored)` then fires for NONE of the next run's
            // first events. Their `cross_` entries stay at Cross {} — power 0, eligible false, hop -1 —
            // which is exactly what this class's own header calls "reads as 'the other axis was silent
            // there'". A consumer running the same instance twice, as a render's output witness does,
            // got a correct first answer and a confidently wrong second one. Shipped in v0.48.0; the
            // fixture below runs the same instance twice and demands the same bits.
            committed_[a] = 0;
        }
    }

    // Read-only, allocation-free. Law 11a: a narrower call than the first is REFUSED — Mid and Side are a
    // statement about a PAIR, and a pair that loses a channel mid-programme is not the same measurement.
    [[nodiscard]] bool process (const float* const* in, int numChannels, int n) noexcept
    {
        if (! prepared_ || finished_ || in == nullptr) return false;
        if (numChannels < 1 || numChannels > maxChannels_ || n < 0) return false;
        if (ranNc_ == 0) ranNc_ = numChannels;
        else if (numChannels != ranNc_) return false;
        for (int c = 0; c < numChannels; ++c) if (in[c] == nullptr) return false;
        if (n == 0) return true;

        for (int i = 0; i < n; ++i)
        {
            float m = in[0][i], s = 0.0f;
            if (numChannels >= 2)
            {
                // ENCODED BEFORE THE FILTERS, as LowEnd does: subtracting two filtered histories would
                // lose Side's relative precision on near-mono material, where Side IS the small difference
                // of two large numbers. A non-finite either side propagates into both axes, which is
                // correct — a hole in R is a hole in both — and each engine counts and locates it itself.
                stereo::MidSide::encode (in[0][i], in[1][i], m, s);
            }
            // …and for one channel, x goes to Mid UNTOUCHED rather than through encode(x, x): the identity
            // holds mathematically but `x + x` overflows above 1.7e38 and would manufacture a hole.
            const float* pm[1] { &m };
            const float* ps[1] { &s };
            axis_ = kMid;  if (! eng_[kMid].process (pm, 1, 1)) return false;
            axis_ = kSide; if (! eng_[kSide].process (ps, 1, 1)) return false;
        }
        samples_ += n;
        return true;
    }

    [[nodiscard]] bool finish() noexcept
    {
        if (! prepared_) return false;
        if (finished_) return true;
        for (int a = 0; a < kAxes; ++a)
        {
            axis_ = a;
            eng_[a].finish();                 // void: the engine has nothing left to refuse at this point
        }
        finished_ = true;
        return true;
    }

    //--------------------------------------------------------------------------------------------------
    bool isPrepared() const noexcept { return prepared_ && eng_[kMid].isPrepared() && eng_[kSide].isPrepared(); }
    bool isFinished() const noexcept { return finished_; }
    double sampleRate() const noexcept { return sampleRate_; }
    int    channels()   const noexcept { return ranNc_; }
    std::int64_t samplesProcessed() const noexcept { return samples_; }

    // SIDE IS ABSENT ONLY WHEN THERE IS ONE CHANNEL. Everything else is a real axis reporting for itself.
    bool sideAbsent() const noexcept { return ranNc_ == 1; }

    const BandBursts& axis (int a) const noexcept { return eng_[a == kSide ? kSide : kMid]; }
    const BandBursts& mid()  const noexcept { return eng_[kMid]; }
    const BandBursts& side() const noexcept { return eng_[kSide]; }

    // Hops whose band energy was EXACTLY zero on this axis — the count beside the threshold, which is this
    // repo's rule for an exact-zero predicate. Bit-identical channels put every judged hop here on Side.
    std::int64_t exactZeroHops (int a) const noexcept
    {
        return a >= 0 && a < kAxes ? exactZero_[a] : 0;
    }

    // The OTHER axis at this event's peak hop. `axis` names the event's own axis.
    Cross crossAt (int a, std::int64_t event) const noexcept
    {
        if (a < 0 || a >= kAxes) return {};
        const auto& v = cross_[a];
        return event >= 0 && event < (std::int64_t) v.size() && event < eng_[a].storedEventCount()
             ? v[(std::size_t) event] : Cross {};
    }

    // The LR4 pair's own ceiling on `peakPower / peakWidePower`, and where it sits. A caller comparing an
    // event's share against a textbook band share is comparing against a number this instrument cannot
    // reach; against the dome it is reading a fraction of what is reachable, which is rate-comparable.
    double domeShare() const noexcept { return domeShare_; }
    double domeHz()    const noexcept { return domeHz_; }

private:
    static void trampoline (void* user, const BandBurstsHopTrace& t) noexcept
    {
        static_cast<StereoBandBursts*> (user)->onHop (t);
    }

    // THE JOIN, AND WHY IT HAS TO WORK IN BOTH DIRECTIONS. Both engines are fed one sample at a time, Mid
    // then Side, so within one input sample the Mid trace for hop k fires and then the Side trace for the
    // same hop k. A request that is only filled by a LATER trace therefore works for Mid and is dead for
    // Side: when Side asks about hop k, Mid's hop k has already gone by, and Mid's next trace is hop k+1,
    // whose start sample does not match. That is exactly what shipped in the first draft — every Side
    // event carried power 0, baseline 0, eligible false, hop -1, which reads as "Mid was silent there"
    // and is the opposite of the truth. The suite was green because it asserted only the Mid direction.
    //
    // So each axis also REMEMBERS its latest hop, and a request is answered from that memory when the
    // other axis has already passed. Forwards and backwards, and a peak that moves later overwrites. The
    // peak rule itself is never re-derived here — the detector publishes where its peak is.
    void onHop (const BandBurstsHopTrace& t) noexcept
    {
        const int a = axis_, other = a == kMid ? kSide : kMid;
        // `fabs(x) <= 0.0` and not `x == 0.0`: exact equality is what is meant — an exactly zero Side is
        // a STRUCTURE, not a small number — but gcc's -Wfloat-equal is an error in this tree and clang
        // says nothing, so the literal spelling passes here and fails the matrix. Same idiom as
        // core::StreamResampler and dynamiceq::LaneDynamics.
        if (std::fabs (t.energy) <= 0.0) ++exactZero_[a];

        // MEAN SQUARES, divided by THIS axis's hop — the units `BandBurst::peakPower` is in, so that a
        // caller can divide one by the other. (Both engines run the same geometry, so the two hop lengths
        // agree; taking it from the axis the energy actually came from is the spelling that stays true.)
        const double hop = (double) eng_[a].hopSamples();
        last_[a].power    = hop > 0.0 ? t.energy / hop : 0.0;
        last_[a].baseline = hop > 0.0 ? t.baseline / hop : 0.0;
        last_[a].eligible = t.eligible;
        last_[a].hop      = t.startSample;

        // Forwards: the other axis asked about this hop before it arrived.
        if (pendingHop_[other] == t.startSample) pending_[other] = last_[a];

        // Ask for this hop on the other axis when it is where our own peak now sits — and if that axis has
        // ALREADY passed this hop, take its answer from memory instead of waiting for a trace that is
        // never coming.
        if (t.inEvent && t.peakAt == t.startSample)
        {
            pendingHop_[a] = t.startSample;
            if (last_[other].hop == t.startSample) pending_[a] = last_[other];
        }

        // Commit when this axis's event count has grown past what we have stored.
        const std::int64_t stored = eng_[a].storedEventCount();
        while (committed_[a] < stored)
        {
            if ((std::size_t) committed_[a] < cross_[a].size()) cross_[a][(std::size_t) committed_[a]] = pending_[a];
            ++committed_[a];
            pending_[a] = Cross {};
            pendingHop_[a] = -1;
        }
    }

    Params params_ {};
    bool   prepared_ = false, finished_ = false;
    double sampleRate_ = 48000.0, domeShare_ = 0.0, domeHz_ = 0.0;
    int    maxChannels_ = 0, ranNc_ = 0, axis_ = 0;
    std::int64_t samples_ = 0;

    BandBursts eng_[kAxes];
    std::vector<Cross> cross_[kAxes];
    Cross        pending_[kAxes] {};
    Cross        last_[kAxes] {};          // each axis's most recent hop, so a passed hop can still answer
    std::int64_t pendingHop_[kAxes] { -1, -1 };
    std::int64_t committed_[kAxes] { 0, 0 };
    std::int64_t exactZero_[kAxes] { 0, 0 };
};

} // namespace felitronics::analysis

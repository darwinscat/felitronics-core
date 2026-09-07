// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/eq/EqTypes.h>
#include <felitronics/eq/MatchedBiquad.h>
#include <felitronics/core/Smoother.h>
#include <felitronics/core/StateGrid.h>
#include <felitronics/eq/Svf.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <cstddef>

namespace felitronics::eq
{

//==============================================================================
// Stateless band design + response — used by EqBand for audio AND by the GUI for a race-free
// curve (the GUI computes from a BandParams snapshot it owns, never reading engine internals).
struct BandDesign
{
    static constexpr int kMaxSections = 8;   // up to 96 dB/oct (16-pole HP/LP = 8 biquads)
    BiquadCoeffs sec[kMaxSections];
    int n = 1;
};

// Design one lane's filter. The design fields (freq/Q/gainDb/slope) are read from the PRIMARY slot,
// lanes[Stereo] — laneView() packs whichever lane is being evaluated into that slot, so this one
// routine serves every lane. `type`/`swept` are the point's shared fields.
inline BandDesign designBand (const BandParams& in, double fs) noexcept
{
    auto finiteOr = [] (double x, double fb) noexcept { return std::isfinite (x) ? x : fb; };
    const LaneParams& src = in.lane (Lane::Stereo);
    const double     freq  = std::clamp (finiteOr (src.freq, 1000.0), 10.0, 0.49 * fs);
    const double     Q     = std::clamp (finiteOr (src.Q, 1.0), 0.05, 40.0);
    const double     gainDb = std::clamp (finiteOr (src.gainDb, 0.0), -30.0, 30.0);
    const int        slope = src.slope;
    const FilterType type  = in.type;
    const bool       swept = in.swept;

    BandDesign d;
    const bool isCut = (type == FilterType::HighPass || type == FilterType::LowPass);

    if (isCut && ! swept)
    {
        // Butterworth cascade: order = slope/6 poles {1..16}; sections = order/2 (+ a first-order
        // section for odd order, e.g. 6/18/30 dB/oct). Per-section Qs are the standard Butterworth
        // pole Qs — a pole pair at s = −sin θ ± j·cos θ has Q = 1/(2 sin θ) — realised by
        // Nyquist-matched biquads. The SIN enumeration is correct for BOTH parities; the previous
        // cos form was even-order-only (a 3rd-order LP read −7.8 dB at fc instead of −3 dB, its
        // corner ~42% low). Pairs are emitted low-Q-first, which for even order reproduces the
        // previous cascade section-for-section (same Q per slot), so shipped even slopes don't move.
        const int  order = std::clamp (slope / 6, 1, 16);
        const bool isHP  = (type == FilterType::HighPass);
        int idx = 0;
        if (order % 2 == 1)
            d.sec[idx++] = isHP ? matched::highpass1 (freq, fs) : matched::lowpass1 (freq, fs);
        const int nSec = order / 2;
        for (int k = 1; k <= nSec && idx < BandDesign::kMaxSections; ++k)
        {
            const double Qk = 1.0 / (2.0 * std::sin ((2.0 * (nSec - k) + 1.0) * kPi / (2.0 * order)));
            d.sec[idx++] = isHP ? matched::highpass (freq, fs, Qk) : matched::lowpass (freq, fs, Qk);
        }
        d.n = idx;
    }
    else if (type == FilterType::Notch && ! swept)
    {
        // Variable-steepness matched band-stop. `order` mirrors HP/LP (slope/6); a notch biquad is
        // inherently 2-sided, so the Butterworth prototype order — and the biquad count — is
        // ceil(order/2), capped at kMaxSections. order∈{1,2} (slope 6/12) → 1 section == today's
        // single matched notch BIT-FOR-BIT (legacy sessions don't drift); 24→2, 48→4, 96→8. Q stays
        // the overall −3 dB bandwidth, independent of order.
        const int order = std::clamp (slope / 6, 1, 16);
        const int m     = std::clamp ((order + 1) / 2, 1, BandDesign::kMaxSections);
        d.n = matched::notchCascade (freq, fs, Q, m, d.sec);
    }
    else if (type == FilterType::BandPass && ! swept)
    {
        // Variable-steepness band-pass — the exact band-PASS mirror of the Notch branch above. `order`
        // mirrors HP/LP (slope/6); a band-pass biquad is inherently 2-sided, so the Butterworth
        // prototype order — and the biquad count — is ceil(order/2) = (order+1)/2, capped at
        // kMaxSections. order∈{1,2} (slope 6/12) → 1 section == today's single matched band-pass
        // BIT-FOR-BIT (legacy sessions don't drift); 24→2, 48→4, 96→8. Q stays the overall −3 dB
        // bandwidth, independent of order; the swept search band keeps the single section / SVF below.
        const int order = std::clamp (slope / 6, 1, 16);
        const int m     = std::clamp ((order + 1) / 2, 1, BandDesign::kMaxSections);
        d.n = matched::bandpassCascade (freq, fs, Q, m, d.sec);
    }
    else if (type == FilterType::Tilt)
    {
        d.sec[0] = matched::lowShelfDb  (freq, fs, -gainDb);   // lows down
        d.sec[1] = matched::highShelfDb (freq, fs,  gainDb);   // highs up -> spectral tilt about f0
        d.n = 2;
    }
    else
    {
        switch (type)
        {
            case FilterType::Bell:      d.sec[0] = matched::peakingDb    (freq, fs, Q, gainDb);     break;
            case FilterType::LowShelf:  d.sec[0] = matched::lowShelfQDb  (freq, fs, gainDb, Q);     break;  // resonant (Q) shelf
            case FilterType::HighShelf: d.sec[0] = matched::highShelfQDb (freq, fs, gainDb, Q);     break;
            case FilterType::HighPass:  d.sec[0] = matched::highpass     (freq, fs, Q);         break;  // swept fallback (single)
            case FilterType::LowPass:   d.sec[0] = matched::lowpass      (freq, fs, Q);         break;  // swept fallback (single)
            case FilterType::BandPass:  d.sec[0] = matched::bandpass     (freq, fs, Q);         break;
            case FilterType::Notch:     d.sec[0] = matched::notch        (freq, fs, Q);         break;
            case FilterType::AllPass:   d.sec[0] = matched::allpass      (freq, fs, Q);         break;
            case FilterType::Tilt:      break;   // handled above (two shelves)
        }
        d.n = 1;
    }
    return d;
}

inline std::complex<double> evalCoeffs (const BiquadCoeffs& c, double w) noexcept
{
    const std::complex<double> z1 = std::polar (1.0, -w), z2 = z1 * z1;
    return (c.b0 + c.b1 * z1 + c.b2 * z2) / (1.0 + c.a1 * z1 + c.a2 * z2);
}

// Race-free GUI response at digital w — computed purely from a caller-owned BandParams snapshot. Reads
// the PRIMARY slot (lanes[Stereo]); feed it a laneView() to evaluate a specific lane.
inline std::complex<double> bandResponse (const BandParams& p, double fs, double w) noexcept
{
    if (! p.on || p.bypass) return { 1.0, 0.0 };
    const BandDesign d = designBand (p, fs);
    std::complex<double> h { 1.0, 0.0 };
    for (int s = 0; s < d.n; ++s) h *= evalCoeffs (d.sec[s], w);
    return h;
}

// The point is "unsplit" (single-ST) when no other lane is ENABLED — the only configuration in which
// the swept search band bites (the search→treat workflow operates on unsplit points). Shared by the
// runtime (EqBand::sweptActive) and the stateless analytics (laneView), so the displayed/FIR'd curve
// and the audio can never disagree about the swept gating.
inline bool onlyStereoEnabled (const BandParams& p) noexcept
{
    return ! (p.lane (Lane::Left).on || p.lane (Lane::Right).on
           || p.lane (Lane::Mid).on  || p.lane (Lane::Side).on);
}

// A single-lane view of a band: `which`'s design fields folded into the primary design slot
// (lanes[Stereo]), so the same designBand()/bandResponse() machinery evaluates that lane. `on` folds
// the point's on with the lane's on; `bypass` folds the point's whole-point bypass with the lane's.
// `swept` survives only for the ST lane of an UNSPLIT point — the same gate the runtime applies — so
// a split point's ST lane evaluates as the matched cascade the engine actually runs. Generalises the
// old sideView: laneView(p, Lane::Side) is the former Side view.
inline BandParams laneView (const BandParams& p, Lane which) noexcept
{
    const LaneParams src = p.lane (which);
    BandParams v = p;
    v.lane (Lane::Stereo) = src;                            // design fields come from `which`
    v.swept  = (which == Lane::Stereo) && p.swept && onlyStereoEnabled (p);
    v.on     = p.on && src.on;                              // fold the point's + lane's enable
    v.bypass = p.bypass || src.bypass;                      // point bypass OR lane bypass mutes it
    return v;
}

// A stereo display / analysis axis (decision #7 display math; matrixResponse is the exact version).
// Axis::Stereo is the ST-lanes-only composite — what a non-stereo bus actually runs (the FIR path's
// mono IR is built from it); the four domain axes each fold the ST lanes in.
enum class Axis { Stereo, Left, Right, Mid, Side };

inline constexpr Lane axisLane (Axis a) noexcept
{
    switch (a)
    {
        case Axis::Stereo: return Lane::Stereo;
        case Axis::Left:   return Lane::Left;
        case Axis::Right:  return Lane::Right;
        case Axis::Mid:    return Lane::Mid;
        case Axis::Side:   return Lane::Side;
    }
    return Lane::Mid;
}

// Composite complex response of the whole bank on ONE stereo axis at digital w, from a caller-owned
// snapshot (race-free). Industry-standard display math (decision #7): axis a = ∏ over bands of
// H_ST · H_a — the ST lane (folds into every axis) times that axis's own-domain lane. An idle lane
// contributes exact identity (laneView folds it off → bandResponse early-returns 1 → skipped, exactly
// as an inactive band is today). Cross-domain coupling is deliberately ignored HERE; matrixResponse is
// the exact 2×2 the FIR path and the tests use.
inline std::complex<double> compositeResponse (const BandParams* bands, int numBands,
                                               double fs, double w, Axis a) noexcept
{
    const Lane axl = axisLane (a);
    std::complex<double> h { 1.0, 0.0 };
    for (int i = 0; i < numBands; ++i)
    {
        h *= bandResponse (laneView (bands[i], Lane::Stereo), fs, w);   // ST folds into every axis
        if (a != Axis::Stereo)
            h *= bandResponse (laneView (bands[i], axl), fs, w);        // the axis's own-domain lane
    }
    return h;
}

// The exact 2×2 complex transfer of the whole bank in the L/R basis at digital w. This is the FIR
// builder's and the tests' source of truth (the GUI uses compositeResponse). Entry hXY = output X per
// unit input Y.
struct ResponseMatrix { std::complex<double> hLL, hLR, hRL, hRR; };

// Fold ONE band into the running 2×2 accumulator, in the NORMATIVE processing order (matches
// EqBand::process). The matrix operators do NOT commute — diag(H_L,H_R) and H_MS only commute when
// H_L == H_R — so this order is load-bearing and identical in the engine, the analytic forms and the
// FIR builder:
//   H_band = H_MS · diag(H_L, H_R) · (H_ST · I)          (H_ST is scalar; its position is free)
//   H_MS   = [[ (H_M+H_S)/2, (H_M−H_S)/2 ],
//             [ (H_M−H_S)/2, (H_M+H_S)/2 ]]              (identity when both M/S lanes idle)
//   acc   <- H_band · acc                                (later bands multiply on the LEFT = series flow)
// This is the SINGLE source of truth for the composition; matrixResponse (complex) and
// matrixResponseZeroPhase (linear-phase, real signed) differ only in how the 5 scalar lane responses are
// evaluated — never in how they combine (so polarity/order can never diverge between the two paths).
inline void accumulateBand (ResponseMatrix& acc, std::complex<double> hST, std::complex<double> hL,
                            std::complex<double> hR, std::complex<double> hM, std::complex<double> hS) noexcept
{
    const std::complex<double> dLL = hST * hL;              // diag(H_L,H_R) · (H_ST·I)
    const std::complex<double> dRR = hST * hR;
    const std::complex<double> mp  = 0.5 * (hM + hS);       // H_MS diagonal
    const std::complex<double> mm  = 0.5 * (hM - hS);       // H_MS off-diagonal
    // H_band = H_MS · diag(dLL, dRR) = [[mp·dLL, mm·dRR], [mm·dLL, mp·dRR]]
    const ResponseMatrix hb { mp * dLL, mm * dRR, mm * dLL, mp * dRR };
    acc = ResponseMatrix {
        hb.hLL * acc.hLL + hb.hLR * acc.hRL,  hb.hLL * acc.hLR + hb.hLR * acc.hRR,
        hb.hRL * acc.hLL + hb.hRR * acc.hRL,  hb.hRL * acc.hLR + hb.hRR * acc.hRR
    };
}

inline ResponseMatrix matrixResponse (const BandParams* bands, int numBands, double fs, double w) noexcept
{
    // The bank product is taken in process order: band 0 first ⇒ rightmost factor.
    ResponseMatrix acc { { 1.0, 0.0 }, { 0.0, 0.0 }, { 0.0, 0.0 }, { 1.0, 0.0 } };   // identity
    for (int i = 0; i < numBands; ++i)
    {
        const BandParams& b = bands[i];
        accumulateBand (acc,
                        bandResponse (laneView (b, Lane::Stereo), fs, w),
                        bandResponse (laneView (b, Lane::Left),   fs, w),
                        bandResponse (laneView (b, Lane::Right),  fs, w),
                        bandResponse (laneView (b, Lane::Mid),    fs, w),
                        bandResponse (laneView (b, Lane::Side),   fs, w));
    }
    return acc;
}

// The ZERO-PHASE 2×2 transfer: every SCALAR LANE response replaced by its magnitude |H_lane| (real,
// non-negative) BEFORE composing. The matrix entries stay REAL but SIGNED — e.g. (|H_M|−|H_S|)/2 keeps
// its sign; taking |·| of the entries instead would break M/S reconstruction. This is the linear-phase
// FIR path's source of truth: a real (signed) spectrum inverse-transforms to a SYMMETRIC IR ⇒ exact
// linear phase, while the sign is preserved. Same accumulateBand composition as matrixResponse.
inline ResponseMatrix matrixResponseZeroPhase (const BandParams* bands, int numBands, double fs, double w) noexcept
{
    auto mag = [&] (const BandParams& b, Lane l) noexcept
    { return std::complex<double> (std::abs (bandResponse (laneView (b, l), fs, w)), 0.0); };

    ResponseMatrix acc { { 1.0, 0.0 }, { 0.0, 0.0 }, { 0.0, 0.0 }, { 1.0, 0.0 } };   // identity
    for (int i = 0; i < numBands; ++i)
    {
        const BandParams& b = bands[i];
        accumulateBand (acc, mag (b, Lane::Stereo), mag (b, Lane::Left), mag (b, Lane::Right),
                        mag (b, Lane::Mid), mag (b, Lane::Side));
    }
    return acc;
}

// Does lane `l` of band `b` ACTUALLY run? (point on && !point bypass && lane on && !lane bypass) — the
// same predicate the engine uses (EqBand::laneOnIn). Basis detection and the FIR builder read this to
// pick a convolver topology.
inline bool laneActive (const BandParams& b, Lane l) noexcept
{
    const LaneParams& lp = b.lane (l);
    return b.on && ! b.bypass && lp.on && ! lp.bypass;
}

// The convolver topology a bank needs (matches convolution::MatrixConvolver::Topology, kept here so the
// eq module owns no convolution dependency): no active L/R lane anywhere ⇒ the whole matrix is
// [[a,b],[b,a]] (diagonal in M/S) ⇒ MSDiag (the plain single-ST case lands here too — today's 2-IR cost);
// no active M/S lane ⇒ diag(H_LL,H_RR) ⇒ LRDiag; otherwise the general 4-entry Full matrix.
enum class MatrixBasis { MSDiag, LRDiag, Full };

inline MatrixBasis detectBasis (const BandParams* bands, int numBands) noexcept
{
    bool anyLR = false, anyMS = false;
    for (int i = 0; i < numBands; ++i)
    {
        const BandParams& b = bands[i];
        anyLR = anyLR || laneActive (b, Lane::Left) || laneActive (b, Lane::Right);
        anyMS = anyMS || laneActive (b, Lane::Mid)  || laneActive (b, Lane::Side);
    }
    if (! anyLR) return MatrixBasis::MSDiag;
    if (! anyMS) return MatrixBasis::LRDiag;
    return MatrixBasis::Full;
}

//==============================================================================
// EqBand — one EQ point, split across up to five placement lanes (ST / L / R / M / S). Owns each
// lane's parameter smoothers, picks the right design per lane, and processes a block in place. Two
// engines under the ST lane:
//   * static treatment band  -> matched biquad(s) (Nyquist-accurate; 24 dB/oct = 2 sections)
//   * swept / search band    -> zero-delay SVF    (clean under a fast fc sweep; single-ST config only)
// Smoothers advance in closed form per block, so coefficients are recomputed once per block while
// freq/Q/gain still move in real time — no zipper, no per-sample biquad redesign. Idle lanes cost zero
// (never designed; their smoothers snap; excluded from the moving/settled check). Lane on/off and lane
// bypass are HARD operator steps (topology reset + snap), exactly the shipped band-enable semantics —
// not amplitude crossfades. `response()` always reports the matched (display) curve so the GUI stays
// honest near Nyquist.
class EqBand
{
public:
    static constexpr int kMaxChannels = Svf::kMaxChannels;
    static constexpr int kMaxSections = BandDesign::kMaxSections;

    // REFUSES a sample rate it cannot honour, and says so, rather than half-honouring it. `fs` reaches
    // every coefficient through tan(pi*f/fs) and its friends, so a zero or non-finite rate does not
    // degrade the filter, it poisons it: measured, prepare(0, 2) and prepare(NaN, 2) each put 63 of 64
    // output samples non-finite on an ordinary 0.25 input, and prepare(1.0, 2) does the same on a
    // HighPass or a Notch. Spelled POSITIVELY so NaN fails — `NaN <= 0.0` is false, which is exactly how
    // the same hole survived in Saturator until P6 F12. A refused band stays unprepared and processBlock
    // does nothing at all, so a caller that ignores the verdict gets silence rather than poison; the
    // return is there for one that wants to know. (EqEngine, the module's front door, marks its own
    // prepare [[nodiscard]]; here that would cost a hundred casts in call sites the gate already covers.)
    bool prepare (double sampleRate, int numChannels, double smoothMs = 30.0) noexcept
    {
        prepared_ = false;                                    // any early return below leaves it unprepared
        // The rate must leave the DESIGN DOMAIN NON-EMPTY, which is a stronger demand than "positive" and
        // is read off the design rather than chosen: every band clamps its frequency to [10 Hz, 0.49*fs],
        // and `std::clamp` with lo > hi is a violated precondition — undefined behaviour, not a clamp.
        // So 0.49*fs must reach 10 Hz, i.e. fs >= 20.41 Hz. Accepting less was measured to produce finite
        // nonsense rather than an obvious failure: at fs = 20 the output looks plausible at 0.377, and at
        // fs = 1 it is 7.2e+28 — which is exactly why "the output is finite" is not a test for this.
        // Spelled positively so NaN fails, as in Saturator (P6 F12) and TruePeakLimiter.
        if (! (std::isfinite (sampleRate) && 0.49 * sampleRate >= 10.0 && sampleRate <= 3.0e6)) return false;
        fs = sampleRate;
        if (numChannels < 1 || numChannels > kMaxChannels) return false;   // law 11(b): BINDING
        ch = numChannels;
        stFreqS_.prepare (fs, smoothMs); stQS_.prepare (fs, smoothMs); stGainS_.prepare (fs, smoothMs);
        for (const Lane l : kMonoLanes)
        {
            LaneRt& rt = laneRt_[(std::size_t) l];
            rt.freqS.prepare (fs, smoothMs); rt.qS.prepare (fs, smoothMs); rt.gainS.prepare (fs, smoothMs);
        }
        svf_.prepare (fs, ch);
        deltaST_.prepare (fs, ch);
        for (const Lane l : kMonoLanes) laneRt_[(std::size_t) l].delta.prepare (fs, 1);
        lastType_ = p.type;
        reset();                                              // prepare's tail IS a stream restart — see reset()
        prepared_ = true;                                     // fully built — processBlock may now run
        return true;
    }

    // Clear everything the PREVIOUS AUDIO left behind, and nothing else: filter memory, the participation
    // ledgers, the dynamic seams and the design-key caches. What a STOP needs — a stage that is skipped for
    // a while and comes back must not re-emit audio from before the gap, but it is still the same stream and
    // the same parameter trajectory, so the smoothers, the targets and the grid phase are none of its
    // business. This is the operation `MasteringChain` and OrbitCab reach for at a bypass edge, and it is
    // the granularity P18 settled on inside this file: a stop clears signal memory, only an explicit
    // reset() restarts the parameter epoch.
    void clearAudioState() noexcept
    {
        resetST();
        for (const Lane l : kMonoLanes) resetLane (laneRt_[(std::size_t) l]);
        // Nothing has run since a stream restart, so nothing can be STOPPING: forget what ran before, or the
        // next block would compute a falling edge against a stream that no longer exists.
        ranMatchedST_ = ranSweptST_ = ranDeltaST_ = 0;
        for (int i = 0; i < kNumLanes; ++i) ranLane_[i] = ranLaneDelta_[i] = false;
        // No live delta may survive: the seam is state like any other, and a restart that kept it would
        // apply the previous stream's gain reduction to the first block of the new one.
        deltaST_.reset(); deltaSTDb_ = 0.0;
        deltaAppliedST_ = std::numeric_limits<double>::quiet_NaN();
        for (const Lane l : kMonoLanes)
        {
            LaneRt& rt = laneRt_[(std::size_t) l];
            rt.delta.reset(); rt.deltaDb = 0.0;
            // Invalidate the design keys too — prepare() may have changed the sample rate, and Svf
            // keeps coefficients baked with g = tan(pi*f/fs) across a prepare().
            rt.deltaApplied = rt.freqApplied = rt.qApplied = std::numeric_limits<double>::quiet_NaN();
        }
        stFreqApplied_ = stQApplied_ = std::numeric_limits<double>::quiet_NaN();
    }

    // A STREAM RESTART, and now it means what the house signature says: the band is left in exactly the
    // state prepare() leaves it in, so `prepare()` is this plus the rate/width setup and cannot drift from
    // it. It used to clear the filters and deliberately leave the freq/Q/gain smoothers wherever the
    // previous stream had pushed them, which made a second render of the same programme start from a
    // different design — measured before the workaround that hid it: a band edited from 0 to +9 dB, 10 ms
    // of the 30 ms ramp rendered, reset, rendered again, differed from a freshly prepared chain by 0.51
    // FULL SCALE. `MasteringChain` carried that as `forceSnap_`, writing every band with its lanes off and
    // then writing them back so the real write would snap; the workaround is gone with the defect.
    // `initialized = false` is the half that is easy to miss and is REQUIRED: without it a write arriving
    // AFTER the reset ramps from the old value instead of snapping, which is the same defect one call
    // later — and it is the right contract on its own, since "the first write after a restart snaps" is
    // exactly what a first write after prepare() already does.
    void reset() noexcept
    {
        clearAudioState();
        grid_.reset();               // a stream restart re-anchors the audio-time maintenance grid
        snapAll();
        updateCoeffs();
        recomputePending = false;
        settlePending_   = false;
        initialized      = false;
    }

    //==========================================================================
    // DYNAMIC GAIN DELTA — the seam a dynamic EQ composes through.
    //
    // The band applies a per-lane gain offset (dB) with a Cytomic SVF bell sitting INSIDE that lane,
    // immediately after its matched static sections and BEFORE the M/S delta-fold. That placement is
    // the whole point and cannot be reproduced from outside: for the M/S lanes the fold stays
    // dM = filt(m) - m with filt = static o svfDelta, so a lane whose delta is 0 still leaves its
    // axis bit-exact.
    //
    // WHO computes the delta is none of this module's business — a detector, a gain computer and
    // ballistics live in whatever composes them (see felitronics::dynamiceq). `eq` deliberately takes
    // a number, so it gains no dependency on `dynamics` and stays readable.
    //
    // The static curve remains a MATCHED biquad; only the moving part is an SVF, because its gain is
    // cheap to modulate without a redesign. Feeding the whole gain through one SVF instead — as a
    // simpler dynamic band does — would forfeit the Nyquist-accurate static response.
    //
    // Call from the processBlock thread. Cheap: it stores a number; coefficients follow at control
    // rate inside process. Deltas are IGNORED unless params().dyn.on, so a non-dynamic band is
    // bit-identical to one that never heard of dynamics.
    void setLaneDeltaDb (Lane l, double db) noexcept
    {
        const double v = std::isfinite (db) ? std::clamp (db, -60.0, 60.0) : 0.0;
        if (l == Lane::Stereo) deltaSTDb_ = v;
        else                   laneRt_[(std::size_t) l].deltaDb = v;
    }

    double laneDeltaDb (Lane l) const noexcept
    {
        return l == Lane::Stereo ? deltaSTDb_ : laneRt_[(std::size_t) l].deltaDb;
    }

    // Call from the SAME thread as processBlock() (typically the audio thread, where the host adapter
    // reads its atomic parameters), or synchronise externally — the engine takes no internal lock. The
    // FIRST call SNAPS every lane (a freshly loaded plugin starts at its settings, no ramp); later
    // calls smooth freq/Q/gain of ACTIVE lanes and SNAP idle lanes (cost-zero) and lanes that just
    // toggled on/off (a hard operator step). Identical params are a no-op. type/slope/swept/on apply on
    // the next block.
    void setParams (const BandParams& npIn) noexcept
    {
        const BandParams np = clamped (npIn);

        if (! initialized)
        {
            p = np;
            snapAll();
            initialized = true;
            recomputePending = true;
            return;
        }
        if (np == p) return;

        // Snapshot each lane's active-ness BEFORE overwriting p: a lane that stays active ramps (smooth
        // edit); one that toggles on/off — or is idle — snaps (idle = cost-zero, toggle = hard step).
        // Classify the change: a write confined to IDLE lanes' freq/Q/gain snaps those smoothers but
        // needs NO recompute — parking automation on a disabled lane burns nothing.
        bool material = (np.type != p.type) || (np.swept != p.swept) || (np.on != p.on) || (np.bypass != p.bypass);
        bool wasOn[kNumLanes];
        for (int i = 0; i < kNumLanes; ++i)
        {
            wasOn[i] = laneOnIn (p, i);
            const LaneParams& a = p.lanes[(std::size_t) i];
            const LaneParams& b = np.lanes[(std::size_t) i];
            const bool topo = (a.on != b.on) || (a.bypass != b.bypass) || (a.slope != b.slope);   // routing / section count
            material = material || topo || ((wasOn[i] || laneOnIn (np, i)) && ! (a == b));         // an active lane's design
        }

        p = np;

        applyLaneTargets (Lane::Stereo, stFreqS_, stQS_, stGainS_, wasOn[(std::size_t) Lane::Stereo]);
        for (const Lane l : kMonoLanes)
        {
            LaneRt& rt = laneRt_[(std::size_t) l];
            applyLaneTargets (l, rt.freqS, rt.qS, rt.gainS, wasOn[(std::size_t) l]);
        }
        if (material) recomputePending = true;
    }

    const BandParams& params() const noexcept { return p; }

    // Audio thread. In-place. RT-safe (no alloc / lock / IO). Normative order: ST per channel, then
    // L on ch0 / R on ch1, then the M/S delta-fold (2-channel only; mono/surround run the ST lane only).
    // Law 11 (DSP-ARCHITECTURE.md §2). Geometry first, before anything moves: a malformed or too-wide
    // call is refused whole and leaves the band exactly as it was. `numSamples` has no limit here — the
    // band sizes no scratch by it.
    [[nodiscard]] bool processBlock (float* const* channels, int numChannels, int numSamples) noexcept
    {
        if (numChannels < 0 || numSamples < 0) return false;   // malformed — never clamped into an index
        if (! prepared_) return false;                         // a refused rate has no coefficients to run
        if (numChannels > ch) return false;                    // width is a LIMIT — law 11(b)
        const int nc = numChannels;

        const bool stRun = laneOn (Lane::Stereo);              // ST runs on every channel (any nc)
        const bool lRun  = nc == 2 && laneOn (Lane::Left);
        const bool rRun  = nc == 2 && laneOn (Lane::Right);
        const bool mRun  = nc == 2 && laneOn (Lane::Mid);
        const bool sRun  = nc == 2 && laneOn (Lane::Side);
        const bool anyRun = stRun || lRun || rRun || mRun || sRun;

        // A cell that STOPS executing loses its sample memory here, on the falling edge. Frozen state is
        // not decayed state: a cell that resumes with it replays a signal from before the gap — measured
        // +8.12 dBFS out of DIGITAL SILENCE when a Side lane sat out a mono stretch. Every gate that can
        // stop a real recursion counts, not just the channel count: a lane switched off, the swept/matched
        // branch, and dyn.on stop one just as completely. This generalises what the band already did for
        // itself when EVERY lane went idle; that case is now simply the one where no cell is left running.
        // A call carrying no samples ran nothing, so it stopped nothing — it must not move an edge.
        // LAW 11a: the edge is clocked by `numSamples`, NOT by the width. This used to also require
        // `nc > 0`, on the reading that "a call that processes no channel ran nothing". That reading is
        // wrong, and it is the same defect P18 closed: a stretch of zero-width calls is a real gap in the
        // stream, and every cell sat it out. Measured on the sibling case (Compressor, 5 ms of lookahead,
        // 4800 zero-width samples, then stereo DIGITAL SILENCE): 0.280315 out of the silence, -11.05 dBFS,
        // the last non-zero sample at index 239 — exactly the frozen 240-sample lookahead line replaying.
        // The other job of the old `nc > 0` — keeping a negative width out of the half-open ranges below,
        // which would index bqST_[s][-1] inside the object where a sanitizer cannot see it — is now done
        // earlier and better by law 11's malformed-call refusal at the top.
        if (numSamples > 0) dropStoppedCells (nc, stRun, p.dyn.on);

        if (numSamples == 0) return true;   // no samples, no time: nothing advances and no edge moves

        // A MATERIAL change — type, slope, swept, on/bypass, or an active lane's design — takes effect
        // HERE, at the call boundary, exactly as it always did. It cannot wait for the next grid tick:
        // the run flags and dropStoppedCells above already read the NEW params, so up to kPeriod-1
        // samples would run the new topology on the old coefficients (a swept band re-enabled after a
        // static stretch would drive `svf_` with whatever its last swept episode left there). What waits
        // for the grid is the RAMP, which is a continuous approximation and has no such coupling.
        // `settlePending_` is ARMED here, not just cleared: this design is made from the smoother values
        // as they stand BEFORE the first tick advances them, so the tick owes one redesign whatever it
        // then finds. Without the arming, a glide that LANDS on its first tick — which is every glide
        // when `smoothMs` is 0, and any glide that starts within `settled()`'s epsilon — read
        // `moving == false` and was never designed at its target at all: measured, a 500 -> 4000 Hz
        // write at smoothMs 0 left the band answering +0.919 dB at 4 kHz where +12 was asked, for ever.
        if (recomputePending) { updateCoeffs(); recomputePending = false; settlePending_ = true; }

        // Dynamics is opt-in per point: with dyn.on false nothing below touches the signal, so a
        // static band is bit-identical to one built before dynamics existed. The delta stays on the
        // CALLER's clock deliberately — it is a value pushed in by setLaneDeltaDb() at whatever cadence
        // its producer runs (felitronics::dynamiceq drives it every 16 samples), so quantising it to
        // kPeriod would round a 1 ms attack up to 1.33 ms. It is an arrival, not maintenance.
        const bool dyn = p.dyn.on;
        if (dyn) updateDeltaCoeffs();

        // THE SEGMENT LOOP. Everything periodic — the parameter ramp, the redesign it earns, and the
        // law-8 flush — happens at `core::StateGrid` boundaries, which are counted in AUDIO samples and
        // therefore fall on the same absolute indices however the caller sliced the stream. The audio
        // itself is untouched by the split: each lane's recursion is per-sample and the M/S fold reads
        // L and R after they were written for the SAME sample, so processing [a,b) lane by lane inside a
        // segment is the identical sequence of operations.
        for (int off = 0; off < numSamples; )
        {
            if (grid_.phase() == 0) tick (anyRun);
            const int seg = grid_.segment (numSamples - off);
            if (anyRun) runAudio (channels, nc, off, seg, stRun, lRun, rRun, mRun, sRun, dyn);
            grid_.advance (seg);
            off += seg;
        }

        // The POISON half, per call and on top of the grid. See eq::Biquad::healPoison(): a NaN's only
        // quality is how soon it goes, and a 16-sample host must not start waiting 64. On a stream whose
        // state stays finite this line cannot change a bit, so the invariance claim above survives it.
        healState();
        return true;
    }

private:
    // One grid tick: advance the parameter ramp by exactly kPeriod samples, redesign if it moved, then
    // run law 8. ALWAYS exactly kPeriod — a partial segment must not advance control, or the ramp would
    // be back on the caller's clock (advance(17) then advance(47) is not advance(64) in binary64).
    void tick (bool anyRun) noexcept
    {
        stFreqS_.advance (core::StateGrid::kPeriod); stQS_.advance (core::StateGrid::kPeriod); stGainS_.advance (core::StateGrid::kPeriod);
        for (const Lane l : kMonoLanes)
        {
            LaneRt& rt = laneRt_[(std::size_t) l];
            rt.freqS.advance (core::StateGrid::kPeriod); rt.qS.advance (core::StateGrid::kPeriod); rt.gainS.advance (core::StateGrid::kPeriod);
        }

        // Recompute only when an ENABLED lane actually moves — a static, settled band skips the trig.
        bool moving = false;
        if (laneOn (Lane::Stereo)) moving = ! (stFreqS_.settled() && stQS_.settled() && stGainS_.settled());
        // ENABLED, not running: a lane that is switched on keeps earning redesigns even while the channel
        // count parks it, because its smoothers advance regardless (above) and updateCoeffs() is what turns
        // a smoothed value into coefficients. Gating this on xRun let a lane settle unseen during a mono
        // stretch and come back filtering at the value of the FIRST mono block — measured 798 Hz for a
        // 500->4000 Hz edit, +0.64 dB where +12 was asked. Design eligibility is nc-agnostic, exactly as
        // updateCoeffs() is, so the mono design sequence now matches the stereo one block for block and the
        // coefficients after a return are BIT-EQUAL to the run where the channel never left. The trig costs
        // what it costs in stereo, and only while the lane actually ramps.
        if (laneOn (Lane::Left))  moving = moving || laneMoving (Lane::Left);
        if (laneOn (Lane::Right)) moving = moving || laneMoving (Lane::Right);
        if (laneOn (Lane::Mid))   moving = moving || laneMoving (Lane::Mid);
        if (laneOn (Lane::Side))  moving = moving || laneMoving (Lane::Side);
        // `settlePending_` is the ramp's last step: the tick on which a smoother finally lands reads
        // `settled()`, so without carrying one more redesign the settled value would never be designed.
        if (settlePending_ || moving) { updateCoeffs(); settlePending_ = moving; }

        // AND THE MOVING BELL FOLLOWS THE STATIC ONE. updateDeltaCoeffs() reads the SMOOTHED freq/Q —
        // "the moving part travels with the static curve during a ramp instead of jumping ahead of it" —
        // so it has to be redesigned wherever those move, which is here. Doing it only at the call
        // boundary (where the arriving delta VALUE is consumed) left the bell designed at the frequency
        // the smoothers held BEFORE this tick, and on a whole-stream call left it there for the whole
        // render: measured on a 500 -> 4000 Hz ramp with a constant -12 dB delta, the static band was at
        // 652.1 Hz while the delta bell sat at 500 Hz, and a 1024-sample call differed from 16-sample
        // calls on 1008 of 1024 samples, worst 0.27. Free when nothing moved — every branch inside is
        // keyed on the applied freq/Q/delta bits.
        if (p.dyn.on) updateDeltaCoeffs();

        // ONLY THE AUDIO STOPS. A parameter ramp runs on the caller's clock, and the smoothers above
        // advance for lanes that are not running — the fully-idle band was the one case that fell out of
        // that rule, because the early return used to sit above them. It made the SAME edit arrive at two
        // different times depending on whether some unrelated lane happened to be on: with a flat 0 dB
        // companion keeping the band alive the design tracked through the gap, and without one the ramp
        // froze and finished ~200 ms AFTER the stream came back (measured: a 500 -> 4000 Hz edit parked at
        // 651.5 Hz for a one-second mono stretch, then 1548 Hz at 10 ms, 3536 at 50 ms, 3984 at 200 ms).
        // Nothing below this line touches a stopped cell: with no lane running there is no state to flush
        // — dropStoppedCells cleared every cell on the falling edge — so the flush is skipped, not owed.
        if (anyRun) flushState();
    }

    // `numSamples` samples of one grid segment, starting at `off`. Byte-for-byte the loop that used to
    // run over the whole call.
    void runAudio (float* const* channels, int nc, int off, int numSamples,
                   bool stRun, bool lRun, bool rRun, bool mRun, bool sRun, bool dyn) noexcept
    {
        // (1) ST lane — per channel. The swept SVF path only runs in the single-ST config; matched
        //     biquads otherwise. (The swept engine legitimately runs on mono and surround too.)
        if (stRun)
        {
            if (sweptActive())
                for (int c = 0; c < nc; ++c)
                {
                    float* d = channels[c] + off;
                    for (int n = 0; n < numSamples; ++n) d[n] = svf_.processSample (c, d[n]);
                }
            else
                for (int c = 0; c < nc; ++c)
                {
                    float* d = channels[c] + off;
                    for (int n = 0; n < numSamples; ++n)
                    {
                        float x = d[n];
                        for (int s = 0; s < designNST_; ++s) x = bqST_[s][c].processSample (x);
                        if (dyn) x = deltaST_.processSample (c, x);   // moving part, after the matched static
                        d[n] = x;
                    }
                }
        }

        // (2) L lane on ch0, R lane on ch1, then (3) the M/S delta-fold — 2-channel only.
        if (nc == 2)
        {
            float* L = channels[0] + off;
            float* R = channels[1] + off;

            if (lRun)
            {
                LaneRt& rt = laneRt_[(std::size_t) Lane::Left];
                for (int n = 0; n < numSamples; ++n)
                {
                    float x = L[n];
                    for (int s = 0; s < rt.designN; ++s) x = rt.bq[s].processSample (x);
                    if (dyn) x = rt.delta.processSample (0, x);
                    L[n] = x;
                }
            }
            if (rRun)
            {
                LaneRt& rt = laneRt_[(std::size_t) Lane::Right];
                for (int n = 0; n < numSamples; ++n)
                {
                    float x = R[n];
                    for (int s = 0; s < rt.designN; ++s) x = rt.bq[s].processSample (x);
                    if (dyn) x = rt.delta.processSample (0, x);
                    R[n] = x;
                }
            }
            if (mRun || sRun)
            {
                LaneRt& M = laneRt_[(std::size_t) Lane::Mid];
                LaneRt& S = laneRt_[(std::size_t) Lane::Side];
                for (int n = 0; n < numSamples; ++n)
                {
                    const float m = 0.5f * (L[n] + R[n]);
                    const float s = 0.5f * (L[n] - R[n]);
                    float dM = 0.0f, dS = 0.0f;
                    // filt = static o svfDelta, so the fold still nulls to zero on an idle lane.
                    if (mRun) { float x = m; for (int k = 0; k < M.designN; ++k) x = M.bq[k].processSample (x); if (dyn) x = M.delta.processSample (0, x); dM = x - m; }
                    if (sRun) { float y = s; for (int k = 0; k < S.designN; ++k) y = S.bq[k].processSample (y); if (dyn) y = S.delta.processSample (0, y); dS = y - s; }
                    L[n] += dM + dS;   // L=M+S, R=M-S: fold deltas back. An idle lane (d=0) leaves its axis bit-exact.
                    R[n] += dM - dS;
                }
            }
        }
    }

public:
    // Complex frequency response of ONE axis at digital w (rad/sample) from the band's current smoothed
    // coefficients: H_ST · H_a (an idle lane's column is designN==0 → contributes identity). Best-effort
    // LIVE readout — for a guaranteed race-free GUI curve prefer the free compositeResponse().
    std::complex<double> response (double w, Axis a = Axis::Mid) const noexcept
    {
        if (! p.on || p.bypass) return { 1.0, 0.0 };
        std::complex<double> h { 1.0, 0.0 };
        for (int s = 0; s < designNST_; ++s) h *= evalCoeffs (coeffsST_[s], w);       // ST folds into every axis
        if (a != Axis::Stereo)
        {
            const LaneRt& rt = laneRt_[(std::size_t) axisLane (a)];                    // laneRt_[Stereo] would be
            for (int s = 0; s < rt.designN; ++s) h *= evalCoeffs (rt.coeffs[s], w);    // empty anyway — explicit > implicit
        }
        return h;
    }

private:
    // A single-signal lane (L / R / M / S): one biquad column (they are single-signal filters by
    // construction). The ST lane is separate — it keeps per-channel columns + the swept SVF.
    struct LaneRt
    {
        Smoother     freqS, qS, gainS;
        Biquad       bq[kMaxSections];
        BiquadCoeffs coeffs[kMaxSections];
        int          designN = 0;
        bool         active  = false;   // designed on the last updateCoeffs()? (drives the topology reset)
        Svf          delta;             // dynamic gain-delta bell, INSIDE the lane (see setLaneDeltaDb)
        double       deltaDb = 0.0;
        double       deltaApplied = std::numeric_limits<double>::quiet_NaN();   // last coeffs pushed,
        double       freqApplied  = std::numeric_limits<double>::quiet_NaN();   // ... and the freq/Q
        double       qApplied     = std::numeric_limits<double>::quiet_NaN();   // they were designed at
    };

    // Same bit-pattern discipline the params use, so a redesign-skip cannot be tripped by -Wfloat-equal
    // and NaN (the "never applied" sentinel) always compares unequal.
    static bool sameBits (double a, double b) noexcept
    {
        return std::bit_cast<std::uint64_t> (a) == std::bit_cast<std::uint64_t> (b);
    }

    // Control-rate: push any changed lane delta into its SVF. The bell tracks the lane's own freq/Q,
    // so the moving part sits exactly where the static curve does.
    void updateDeltaCoeffs() noexcept
    {
        // Keyed on freq/Q as well as the delta itself: holding a steady delta while the user drags the
        // node would otherwise leave the moving bell parked at the frequency it was designed for. The
        // SMOOTHED values are used, so the moving part travels with the static curve during a ramp
        // instead of jumping to the target ahead of it.
        {
            const double f = stFreqS_.value(), q = stQS_.value();
            if (! sameBits (deltaSTDb_, deltaAppliedST_) || ! sameBits (f, stFreqApplied_)
                || ! sameBits (q, stQApplied_))
            {
                deltaST_.setParams (FilterType::Bell, f, q, deltaSTDb_);
                deltaAppliedST_ = deltaSTDb_;
                stFreqApplied_ = f; stQApplied_ = q;
            }
        }
        for (const Lane l : kMonoLanes)
        {
            LaneRt& rt = laneRt_[(std::size_t) l];
            const double f = rt.freqS.value(), q = rt.qS.value();
            if (! sameBits (rt.deltaDb, rt.deltaApplied) || ! sameBits (f, rt.freqApplied)
                || ! sameBits (q, rt.qApplied))
            {
                rt.delta.setParams (FilterType::Bell, f, q, rt.deltaDb);
                rt.deltaApplied = rt.deltaDb;
                rt.freqApplied = f; rt.qApplied = q;
            }
        }
    }

    static constexpr Lane kMonoLanes[4] { Lane::Left, Lane::Right, Lane::Mid, Lane::Side };

    bool pointActive() const noexcept { return p.on && ! p.bypass; }

    // Does lane `l` want to run? (point active AND lane on AND lane not bypassed) — nc-agnostic; the
    // process() step additionally gates L/R/M/S on nc==2.
    bool laneOn (Lane l) const noexcept
    {
        const LaneParams& lp = p.lane (l);
        return pointActive() && lp.on && ! lp.bypass;
    }

    static bool laneOnIn (const BandParams& q, int i) noexcept
    {
        const LaneParams& lp = q.lanes[(std::size_t) i];
        return q.on && ! q.bypass && lp.on && ! lp.bypass;
    }

    bool laneMoving (Lane l) const noexcept
    {
        const LaneRt& rt = laneRt_[(std::size_t) l];
        return ! (rt.freqS.settled() && rt.qS.settled() && rt.gainS.settled());
    }

    // Unsplit (single-ST) check — delegates to the free onlyStereoEnabled() so the runtime gate and
    // the stateless analytics (laneView) can never disagree about swept.
    bool onlyStereo() const noexcept { return onlyStereoEnabled (p); }

    // The swept (SVF) engine only runs for types the single SVF stage can realise, AND only in the
    // single-ST configuration. Tilt has no one-SVF realisation with a unity pivot beyond ~6 dB, so a
    // swept Tilt runs the matched two-shelf design instead.
    bool sweptActive() const noexcept
    {
        return p.swept && p.type != FilterType::Tilt && onlyStereo();
    }

    void applyLaneTargets (Lane l, Smoother& fS, Smoother& qS2, Smoother& gS, bool wasOn) noexcept
    {
        const LaneParams& lp = p.lane (l);
        if (laneOn (l) && wasOn) { fS.setTarget (lp.freq); qS2.setTarget (lp.Q); gS.setTarget (lp.gainDb); }  // smooth edit
        else                     { fS.snap (lp.freq);      qS2.snap (lp.Q);      gS.snap (lp.gainDb); }        // idle / hard step
    }

    void snapAll() noexcept
    {
        const LaneParams& st = p.lane (Lane::Stereo);
        stFreqS_.snap (st.freq); stQS_.snap (st.Q); stGainS_.snap (st.gainDb);
        for (const Lane l : kMonoLanes)
        {
            LaneRt& rt = laneRt_[(std::size_t) l];
            const LaneParams& lp = p.lane (l);
            rt.freqS.snap (lp.freq); rt.qS.snap (lp.Q); rt.gainS.snap (lp.gainDb);
        }
    }

    BandParams clamped (const BandParams& in) const noexcept
    {
        auto finiteOr = [] (double x, double fb) noexcept { return std::isfinite (x) ? x : fb; };
        BandParams np = in;
        for (int i = 0; i < kNumLanes; ++i)
        {
            LaneParams& lp = np.lanes[(std::size_t) i];
            lp.freq   = std::clamp (finiteOr (lp.freq, 1000.0), 10.0, 0.49 * fs);
            lp.Q      = std::clamp (finiteOr (lp.Q, 1.0), 0.05, 40.0);
            lp.gainDb = std::clamp (finiteOr (lp.gainDb, 0.0), -30.0, 30.0);
        }
        // Dynamics rails. Without them an absurd threshold makes digital silence earn the full range
        // (the detector level is floored, so a threshold below that floor is permanently exceeded),
        // and a huge range overflows float on its way into the gain follower, poisoning the meter
        // with NaN until the next reset. +24 dBFS on the threshold because float hosts legitimately
        // run hot inside a chain.
        np.dyn.rangeDb = std::clamp (finiteOr (np.dyn.rangeDb, 0.0), -30.0, 30.0);
        np.dyn.thrDb   = std::clamp (finiteOr (np.dyn.thrDb, -24.0), -120.0, 24.0);
        np.dyn.atk     = std::clamp (finiteOr (np.dyn.atk, 0.5), 0.0, 1.0);
        np.dyn.rel     = std::clamp (finiteOr (np.dyn.rel, 0.5), 0.0, 1.0);
        return np;
    }

    void updateCoeffs() noexcept
    {
        const bool typeChanged = (p.type != lastType_);   // shared type change resets EVERY lane's column
        lastType_ = p.type;

        // ---- ST lane (per-channel columns; the only lane that can sweep) ----
        {
            const bool active = laneOn (Lane::Stereo);
            const bool sw     = sweptActive();
            BandDesign d; d.n = 0;
            if (active)
            {
                BandParams sp = p;
                LaneParams& s = sp.lane (Lane::Stereo);
                s.freq = stFreqS_.value(); s.Q = stQS_.value(); s.gainDb = stGainS_.value();
                sp.swept = sw;   // single-stage swept design ONLY in the single-ST config; else the matched cascade
                d = designBand (sp, fs);   // sp.type is the point's shared field
            }
            // Topology switch (section count, swept<->static, activation, OR the shared type) → clear
            // the ST columns so a re-activated section never resumes a stale tail. (A toggle set AND
            // unset between two process blocks intentionally does NOT reset: no block ran in the
            // interim, so the state is one continuous stream — resetting would be the artifact.)
            if (d.n != designNST_ || sw != lastSwept_ || active != stActive_ || typeChanged) resetST();
            designNST_ = d.n;
            lastSwept_ = sw;
            stActive_  = active;
            for (int s = 0; s < d.n; ++s) coeffsST_[s] = d.sec[s];
            for (int c = 0; c < ch; ++c)
                for (int s = 0; s < d.n; ++s) bqST_[s][c].setCoeffs (coeffsST_[s]);
            if (active && sw) svf_.setParams (p.type, stFreqS_.value(), stQS_.value(), stGainS_.value());  // single SVF stage
        }

        // ---- L / R / M / S lanes (one column each) ----
        for (const Lane l : kMonoLanes)
        {
            LaneRt& rt = laneRt_[(std::size_t) l];
            const bool active = laneOn (l);
            BandDesign d; d.n = 0;
            if (active)
            {
                BandParams sp = p;
                LaneParams& s = sp.lane (Lane::Stereo);
                s = p.lane (l);                       // this lane's design fields (incl. slope) into the slot
                s.freq = rt.freqS.value(); s.Q = rt.qS.value(); s.gainDb = rt.gainS.value();
                sp.swept = false;                     // only the ST lane sweeps
                d = designBand (sp, fs);
            }
            if (d.n != rt.designN || active != rt.active || typeChanged) resetLane (rt);
            rt.designN = d.n;
            rt.active  = active;
            for (int s = 0; s < d.n; ++s) { rt.coeffs[s] = d.sec[s]; rt.bq[s].setCoeffs (rt.coeffs[s]); }
        }
    }

    // Law 8, once per `core::StateGrid` period. DESIGNED SECTIONS ONLY, which is not an optimisation
    // dressed as one: a section past `designN` was zeroed by resetST()/resetLane() at the topology change
    // that shrank the count and is never written by the audio loop, so visiting it was always a
    // guaranteed no-op — 48 Biquads scanned where a stereo one-section bell has 2. That headroom is what
    // pays for running this eight times more often than a 512-sample host used to.
    void flushState() noexcept
    {
        if (p.dyn.on)   // Law 8 for the moving part too — these are feedback kernels like any other
        {
            deltaST_.flushDenormals();
            for (const Lane l : kMonoLanes) laneRt_[(std::size_t) l].delta.flushDenormals();
        }
        if (sweptActive()) svf_.flushDenormals();
        else for (int c = 0; c < ch; ++c)
                 for (int s = 0; s < designNST_; ++s) bqST_[s][c].flushDenormals();
        for (const Lane l : kMonoLanes)
        {
            LaneRt& rt = laneRt_[(std::size_t) l];
            for (int s = 0; s < rt.designN; ++s) rt.bq[s].flushDenormals();
        }
    }

    // The poison half of the same sweep, run at the END OF EVERY CALL rather than on the grid — see
    // eq::Biquad::healPoison(). Cannot change a bit while the state is finite, so it is invisible to the
    // slicing-invariance claim; it exists so a host with a block SHORTER than a grid period keeps the
    // recovery it has today instead of waiting for the next boundary.
    void healState() noexcept
    {
        if (p.dyn.on)
        {
            deltaST_.healPoison();
            for (const Lane l : kMonoLanes) laneRt_[(std::size_t) l].delta.healPoison();
        }
        if (sweptActive()) svf_.healPoison();
        else for (int c = 0; c < ch; ++c)
                 for (int s = 0; s < designNST_; ++s) bqST_[s][c].healPoison();
        for (const Lane l : kMonoLanes)
        {
            LaneRt& rt = laneRt_[(std::size_t) l];
            for (int s = 0; s < rt.designN; ++s) rt.bq[s].healPoison();
        }
    }

    void resetST() noexcept
    {
        for (int c = 0; c < kMaxChannels; ++c)
            for (int s = 0; s < kMaxSections; ++s) bqST_[s][c].reset();
        svf_.reset();
    }

    static void resetLane (LaneRt& rt) noexcept
    {
        for (int s = 0; s < kMaxSections; ++s) rt.bq[s].reset();
    }

    // Clear the sample memory of every cell that ran on the previous sample-bearing call and does not run
    // on this one. SIGNAL HISTORY ONLY: seams (deltaSTDb_, LaneRt::deltaDb), smoothers and the applied-value
    // caches are current CONTROL, not memory of past audio — they say what the band should do now, and a
    // lane coming back must obey the latest command, not the one in force when it left. (setLaneDeltaDb's
    // producer decides when a delta has expired; LaneDynamics already zeroes its own seam on the same edge.)
    // The running set of channels is always [0, nc), so recording the COUNT that ran says which columns did,
    // and the falling edge is the half-open range [now, then). Per channel, never wholesale: a channel that
    // never left owes nothing to one that did, which is why Svf grew resetChannel().
    void dropStoppedCells (int nc, bool stRun, bool dynOn) noexcept
    {
        const bool swept = sweptActive();
        const int nowMatched = (stRun && ! swept)          ? nc : 0;
        const int nowSwept   = (stRun &&   swept)          ? nc : 0;
        const int nowDelta   = (stRun && ! swept && dynOn) ? nc : 0;

        for (int c = nowMatched; c < ranMatchedST_; ++c)
            for (int s = 0; s < kMaxSections; ++s) bqST_[s][c].reset();
        for (int c = nowSwept; c < ranSweptST_; ++c) svf_.resetChannel (c);
        for (int c = nowDelta;  c < ranDeltaST_;  ++c) deltaST_.resetChannel (c);
        ranMatchedST_ = nowMatched; ranSweptST_ = nowSwept; ranDeltaST_ = nowDelta;

        for (const Lane l : kMonoLanes)
        {
            const std::size_t i = (std::size_t) l;
            LaneRt&    rt   = laneRt_[i];
            const bool runs = (nc == 2) && laneOn (l);          // L/R/M/S are stereo-only, as everywhere here
            if (ranLane_[i]      && ! runs)            for (int s = 0; s < kMaxSections; ++s) rt.bq[s].reset();
            if (ranLaneDelta_[i] && ! (runs && dynOn)) rt.delta.resetChannel (0);   // one column by construction
            ranLane_[i]      = runs;
            ranLaneDelta_[i] = runs && dynOn;
        }
    }

    BandParams p;
    double fs = 44100.0;
    int    ch = 2;
    bool   prepared_ = false;                 // true only after a fully-successful prepare()
    bool   recomputePending = true;   // a MATERIAL change is pending: applied at the next call boundary
    bool   settlePending_   = false;  // the ramp's last redesign, owed to the tick AFTER a smoother lands
    bool   initialized = false;
    core::StateGrid grid_;            // audio-time maintenance + control clock (see core/StateGrid.h)
    FilterType lastType_ = FilterType::Bell;

    // What actually ADVANCED state on the previous sample-bearing call — the channel count that ran, 0 for
    // "this cell did not run at all". Deliberately not stActive_/LaneRt::active: those describe the DESIGN,
    // are nc-agnostic, and are refreshed only when updateCoeffs() happens to run, so they lag execution.
    int    ranMatchedST_ = 0, ranSweptST_ = 0, ranDeltaST_ = 0;
    bool   ranLane_[kNumLanes] {}, ranLaneDelta_[kNumLanes] {};

    // ST lane: per-channel biquad columns (mono→surround→ambisonics) + the swept SVF.
    Smoother     stFreqS_, stQS_, stGainS_;
    Biquad       bqST_[kMaxSections][kMaxChannels];
    BiquadCoeffs coeffsST_[kMaxSections];
    int          designNST_ = 1;
    bool         stActive_  = false;
    bool         lastSwept_ = false;
    Svf          svf_;
    Svf          deltaST_;             // dynamic gain-delta bell for the ST lane (per channel)
    double       deltaSTDb_ = 0.0;
    double       deltaAppliedST_ = std::numeric_limits<double>::quiet_NaN();   // last coeffs pushed,
    double       stFreqApplied_  = std::numeric_limits<double>::quiet_NaN();   // ... and the freq/Q
    double       stQApplied_     = std::numeric_limits<double>::quiet_NaN();   // they were designed at

    // L / R / M / S lanes. Indexed by Lane; the [Lane::Stereo] entry is unused (the ST lane keeps the
    // per-channel columns above) — a trivial, deliberate slot for index-by-enum clarity.
    LaneRt laneRt_[kNumLanes];
};

} // namespace felitronics::eq

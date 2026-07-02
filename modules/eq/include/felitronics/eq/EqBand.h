// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/eq/EqTypes.h>
#include <felitronics/eq/MatchedBiquad.h>
#include <felitronics/core/Smoother.h>
#include <felitronics/eq/Svf.h>

#include <algorithm>
#include <complex>
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

inline ResponseMatrix matrixResponse (const BandParams* bands, int numBands, double fs, double w) noexcept
{
    // Per band, composed in the NORMATIVE processing order (matches EqBand::process). The matrix
    // operators do NOT commute — diag(H_L,H_R) and H_MS only commute when H_L == H_R — so this order
    // is load-bearing and identical in the engine, this analytic form and the FIR builder:
    //   H_band = H_MS · diag(H_L, H_R) · (H_ST · I)          (H_ST is scalar; its position is free)
    //   H_MS   = [[ (H_M+H_S)/2, (H_M−H_S)/2 ],
    //             [ (H_M−H_S)/2, (H_M+H_S)/2 ]]              (identity when both M/S lanes idle)
    // The bank product is taken in process order: band 0 first ⇒ rightmost factor.
    ResponseMatrix acc { { 1.0, 0.0 }, { 0.0, 0.0 }, { 0.0, 0.0 }, { 1.0, 0.0 } };   // identity
    for (int i = 0; i < numBands; ++i)
    {
        const BandParams& b = bands[i];
        const std::complex<double> hST = bandResponse (laneView (b, Lane::Stereo), fs, w);
        const std::complex<double> hL  = bandResponse (laneView (b, Lane::Left),   fs, w);
        const std::complex<double> hR  = bandResponse (laneView (b, Lane::Right),  fs, w);
        const std::complex<double> hM  = bandResponse (laneView (b, Lane::Mid),    fs, w);
        const std::complex<double> hS  = bandResponse (laneView (b, Lane::Side),   fs, w);

        const std::complex<double> dLL = hST * hL;              // diag(H_L,H_R) · (H_ST·I)
        const std::complex<double> dRR = hST * hR;
        const std::complex<double> mp  = 0.5 * (hM + hS);       // H_MS diagonal
        const std::complex<double> mm  = 0.5 * (hM - hS);       // H_MS off-diagonal
        // H_band = H_MS · diag(dLL, dRR) = [[mp·dLL, mm·dRR], [mm·dLL, mp·dRR]]
        const ResponseMatrix hb { mp * dLL, mm * dRR, mm * dLL, mp * dRR };
        // acc <- H_band · acc  (later bands multiply on the LEFT, matching the series signal flow)
        acc = ResponseMatrix {
            hb.hLL * acc.hLL + hb.hLR * acc.hRL,  hb.hLL * acc.hLR + hb.hLR * acc.hRR,
            hb.hRL * acc.hLL + hb.hRR * acc.hRL,  hb.hRL * acc.hLR + hb.hRR * acc.hRR
        };
    }
    return acc;
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

    void prepare (double sampleRate, int numChannels, double smoothMs = 30.0) noexcept
    {
        fs = sampleRate;
        ch = numChannels < 1 ? 1 : (numChannels > kMaxChannels ? kMaxChannels : numChannels);
        stFreqS_.prepare (fs, smoothMs); stQS_.prepare (fs, smoothMs); stGainS_.prepare (fs, smoothMs);
        for (const Lane l : kMonoLanes)
        {
            LaneRt& rt = laneRt_[(std::size_t) l];
            rt.freqS.prepare (fs, smoothMs); rt.qS.prepare (fs, smoothMs); rt.gainS.prepare (fs, smoothMs);
        }
        snapAll();
        svf_.prepare (fs, ch);
        lastType_ = p.type;
        updateCoeffs();
        recomputePending = false;
        initialized = false;
        bandWasActive = false;
        reset();
    }

    void reset() noexcept
    {
        resetST();
        for (const Lane l : kMonoLanes) resetLane (laneRt_[(std::size_t) l]);
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
    void processBlock (float* const* channels, int numChannels, int numSamples) noexcept
    {
        const int nc = numChannels < ch ? numChannels : ch;

        const bool stRun = laneOn (Lane::Stereo);              // ST runs on every channel (any nc)
        const bool lRun  = nc == 2 && laneOn (Lane::Left);
        const bool rRun  = nc == 2 && laneOn (Lane::Right);
        const bool mRun  = nc == 2 && laneOn (Lane::Mid);
        const bool sRun  = nc == 2 && laneOn (Lane::Side);
        const bool anyRun = stRun || lRun || rRun || mRun || sRun;

        if (! anyRun || numSamples <= 0)
        {
            if (! anyRun && bandWasActive) { reset(); bandWasActive = false; }   // clear tails so re-enabling doesn't pop
            return;
        }
        bandWasActive = true;

        // Advance every lane's smoothers (closed form; an idle/snapped smoother is a settled no-op).
        stFreqS_.advance (numSamples); stQS_.advance (numSamples); stGainS_.advance (numSamples);
        for (const Lane l : kMonoLanes)
        {
            LaneRt& rt = laneRt_[(std::size_t) l];
            rt.freqS.advance (numSamples); rt.qS.advance (numSamples); rt.gainS.advance (numSamples);
        }

        // Recompute only when a RUNNING lane actually moves — a static, settled band skips the trig.
        bool moving = false;
        if (stRun) moving = moving || ! (stFreqS_.settled() && stQS_.settled() && stGainS_.settled());
        if (lRun)  moving = moving || laneMoving (Lane::Left);
        if (rRun)  moving = moving || laneMoving (Lane::Right);
        if (mRun)  moving = moving || laneMoving (Lane::Mid);
        if (sRun)  moving = moving || laneMoving (Lane::Side);
        if (recomputePending || moving) { updateCoeffs(); recomputePending = moving; }

        // (1) ST lane — per channel. The swept SVF path only runs in the single-ST config; matched
        //     biquads otherwise. (The swept engine legitimately runs on mono and surround too.)
        if (stRun)
        {
            if (sweptActive())
                for (int c = 0; c < nc; ++c)
                {
                    float* d = channels[c];
                    for (int n = 0; n < numSamples; ++n) d[n] = svf_.processSample (c, d[n]);
                }
            else
                for (int c = 0; c < nc; ++c)
                {
                    float* d = channels[c];
                    for (int n = 0; n < numSamples; ++n)
                    {
                        float x = d[n];
                        for (int s = 0; s < designNST_; ++s) x = bqST_[s][c].processSample (x);
                        d[n] = x;
                    }
                }
        }

        // (2) L lane on ch0, R lane on ch1, then (3) the M/S delta-fold — 2-channel only.
        if (nc == 2)
        {
            float* L = channels[0];
            float* R = channels[1];

            if (lRun)
            {
                LaneRt& rt = laneRt_[(std::size_t) Lane::Left];
                for (int n = 0; n < numSamples; ++n)
                {
                    float x = L[n];
                    for (int s = 0; s < rt.designN; ++s) x = rt.bq[s].processSample (x);
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
                    if (mRun) { float x = m; for (int k = 0; k < M.designN; ++k) x = M.bq[k].processSample (x); dM = x - m; }
                    if (sRun) { float y = s; for (int k = 0; k < S.designN; ++k) y = S.bq[k].processSample (y); dS = y - s; }
                    L[n] += dM + dS;   // L=M+S, R=M-S: fold deltas back. An idle lane (d=0) leaves its axis bit-exact.
                    R[n] += dM - dS;
                }
            }
        }

        flushState();   // per-block denormal guard
    }

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
    };

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

    void flushState() noexcept
    {
        if (sweptActive()) svf_.flushDenormals();
        else for (int c = 0; c < ch; ++c)
                 for (int s = 0; s < kMaxSections; ++s) bqST_[s][c].flushDenormals();
        for (const Lane l : kMonoLanes)
        {
            LaneRt& rt = laneRt_[(std::size_t) l];
            for (int s = 0; s < kMaxSections; ++s) rt.bq[s].flushDenormals();
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

    BandParams p;
    double fs = 44100.0;
    int    ch = 2;
    bool   recomputePending = true;
    bool   initialized = false;
    bool   bandWasActive = false;
    FilterType lastType_ = FilterType::Bell;

    // ST lane: per-channel biquad columns (mono→surround→ambisonics) + the swept SVF.
    Smoother     stFreqS_, stQS_, stGainS_;
    Biquad       bqST_[kMaxSections][kMaxChannels];
    BiquadCoeffs coeffsST_[kMaxSections];
    int          designNST_ = 1;
    bool         stActive_  = false;
    bool         lastSwept_ = false;
    Svf          svf_;

    // L / R / M / S lanes. Indexed by Lane; the [Lane::Stereo] entry is unused (the ST lane keeps the
    // per-channel columns above) — a trivial, deliberate slot for index-by-enum clarity.
    LaneRt laneRt_[kNumLanes];
};

} // namespace felitronics::eq

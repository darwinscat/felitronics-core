// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/eq/EqTypes.h>
#include <felitronics/eq/EqTypes.h>
#include <felitronics/eq/MatchedBiquad.h>
#include <felitronics/core/Smoother.h>
#include <felitronics/eq/Svf.h>

#include <algorithm>
#include <complex>

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

inline BandDesign designBand (const BandParams& in, double fs) noexcept
{
    auto finiteOr = [] (double x, double fb) noexcept { return std::isfinite (x) ? x : fb; };
    BandParams p = in;
    p.freq   = std::clamp (finiteOr (p.freq, 1000.0), 10.0, 0.49 * fs);
    p.Q      = std::clamp (finiteOr (p.Q, 1.0), 0.05, 40.0);
    p.gainDb = std::clamp (finiteOr (p.gainDb, 0.0), -30.0, 30.0);

    BandDesign d;
    const bool isCut = (p.type == FilterType::HighPass || p.type == FilterType::LowPass);

    if (isCut && ! p.swept)
    {
        // Butterworth cascade: order = slope/6 poles {1..16}; sections = order/2 (+ a first-order
        // section for odd order, e.g. 6/18/30 dB/oct). Per-section Qs are the standard Butterworth
        // pole Qs — a pole pair at s = −sin θ ± j·cos θ has Q = 1/(2 sin θ) — realised by
        // Nyquist-matched biquads. The SIN enumeration is correct for BOTH parities; the previous
        // cos form was even-order-only (a 3rd-order LP read −7.8 dB at fc instead of −3 dB, its
        // corner ~42% low). Pairs are emitted low-Q-first, which for even order reproduces the
        // previous cascade section-for-section (same Q per slot), so shipped even slopes don't move.
        const int  order = std::clamp (p.slope / 6, 1, 16);
        const bool isHP  = (p.type == FilterType::HighPass);
        int idx = 0;
        if (order % 2 == 1)
            d.sec[idx++] = isHP ? matched::highpass1 (p.freq, fs) : matched::lowpass1 (p.freq, fs);
        const int nSec = order / 2;
        for (int k = 1; k <= nSec && idx < BandDesign::kMaxSections; ++k)
        {
            const double Qk = 1.0 / (2.0 * std::sin ((2.0 * (nSec - k) + 1.0) * kPi / (2.0 * order)));
            d.sec[idx++] = isHP ? matched::highpass (p.freq, fs, Qk) : matched::lowpass (p.freq, fs, Qk);
        }
        d.n = idx;
    }
    else if (p.type == FilterType::Notch && ! p.swept)
    {
        // Variable-steepness matched band-stop. `order` mirrors HP/LP (slope/6); a notch biquad is
        // inherently 2-sided, so the Butterworth prototype order — and the biquad count — is
        // ceil(order/2), capped at kMaxSections. order∈{1,2} (slope 6/12) → 1 section == today's
        // single matched notch BIT-FOR-BIT (legacy sessions don't drift); 24→2, 48→4, 96→8. Q stays
        // the overall −3 dB bandwidth, independent of order.
        const int order = std::clamp (p.slope / 6, 1, 16);
        const int m     = std::clamp ((order + 1) / 2, 1, BandDesign::kMaxSections);
        d.n = matched::notchCascade (p.freq, fs, p.Q, m, d.sec);
    }
    else if (p.type == FilterType::Tilt)
    {
        d.sec[0] = matched::lowShelfDb  (p.freq, fs, -p.gainDb);   // lows down
        d.sec[1] = matched::highShelfDb (p.freq, fs,  p.gainDb);   // highs up -> spectral tilt about f0
        d.n = 2;
    }
    else
    {
        switch (p.type)
        {
            case FilterType::Bell:      d.sec[0] = matched::peakingDb    (p.freq, fs, p.Q, p.gainDb);     break;
            case FilterType::LowShelf:  d.sec[0] = matched::lowShelfQDb  (p.freq, fs, p.gainDb, p.Q);     break;  // resonant (Q) shelf
            case FilterType::HighShelf: d.sec[0] = matched::highShelfQDb (p.freq, fs, p.gainDb, p.Q);     break;
            case FilterType::HighPass:  d.sec[0] = matched::highpass     (p.freq, fs, p.Q);         break;  // swept fallback (single)
            case FilterType::LowPass:   d.sec[0] = matched::lowpass      (p.freq, fs, p.Q);         break;  // swept fallback (single)
            case FilterType::BandPass:  d.sec[0] = matched::bandpass     (p.freq, fs, p.Q);         break;
            case FilterType::Notch:     d.sec[0] = matched::notch        (p.freq, fs, p.Q);         break;
            case FilterType::AllPass:   d.sec[0] = matched::allpass      (p.freq, fs, p.Q);         break;
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

// Race-free GUI response at digital w — computed purely from a caller-owned BandParams snapshot.
inline std::complex<double> bandResponse (const BandParams& p, double fs, double w) noexcept
{
    if (! p.on || p.bypass) return { 1.0, 0.0 };
    const BandDesign d = designBand (p, fs);
    std::complex<double> h { 1.0, 0.0 };
    for (int s = 0; s < d.n; ++s) h *= evalCoeffs (d.sec[s], w);
    return h;
}

// A "Side lane" view of a band: its s* fields mapped into the primary design slots, so the same
// designBand()/bandResponse() machinery evaluates the Side lane. Only meaningful when p.ms (the
// engine splits M/S only on a 2-channel signal). on/bypass fold in the Side lane's enable state.
inline BandParams sideView (const BandParams& p) noexcept
{
    BandParams s = p;
    s.on     = p.on && p.sOn;     // Side runs only if the band is on AND its Side lane is enabled
    s.bypass = p.sBypass;
    s.type   = p.sType;
    s.freq   = p.sFreq;
    s.Q      = p.sQ;
    s.gainDb = p.sGainDb;
    s.slope  = p.sSlope;
    s.swept  = false;             // M/S lanes are matched-only (no SVF sweep)
    return s;
}

// Composite complex response of the whole bank on ONE stereo axis at digital w, from a caller-owned
// snapshot (race-free). side=false → the Mid axis (also the plain L=R response when no band is M/S);
// side=true → the Side axis. A non-M/S band contributes its single (stereo) response to both axes; an
// M/S band contributes its Mid lane to the Mid axis and its Side lane to the Side axis.
inline std::complex<double> compositeResponse (const BandParams* bands, int numBands,
                                               double fs, double w, bool side) noexcept
{
    std::complex<double> h { 1.0, 0.0 };
    for (int i = 0; i < numBands; ++i)
        h *= (side && bands[i].ms) ? bandResponse (sideView (bands[i]), fs, w)
                                   : bandResponse (bands[i], fs, w);
    return h;
}

//==============================================================================
// teq::EqBand — one EQ band. Owns its parameter smoothers, picks the right design, and
// processes a block in place (mono or stereo). Two engines under one band:
//   * static treatment band  -> matched biquad(s) (Nyquist-accurate; 24 dB/oct = 2 sections)
//   * swept / search band    -> zero-delay SVF    (clean under a fast fc sweep)
// Smoothers advance in closed form per block, so coefficients are recomputed once per block
// while freq/Q/gain still move in real time — no zipper, no per-sample biquad redesign.
// `response()` always reports the matched (display) curve so the GUI stays honest near Nyquist.
class EqBand
{
public:
    static constexpr int kMaxChannels = Svf::kMaxChannels;
    static constexpr int kMaxSections = BandDesign::kMaxSections;

    void prepare (double sampleRate, int numChannels, double smoothMs = 30.0) noexcept
    {
        fs = sampleRate;
        ch = numChannels < 1 ? 1 : (numChannels > kMaxChannels ? kMaxChannels : numChannels);
        freqS.prepare (fs, smoothMs);
        qS.prepare    (fs, smoothMs);
        gainS.prepare (fs, smoothMs);
        sFreqS.prepare (fs, smoothMs);
        sQS.prepare    (fs, smoothMs);
        sGainS.prepare (fs, smoothMs);
        freqS.snap (p.freq); qS.snap (p.Q); gainS.snap (p.gainDb);
        sFreqS.snap (p.sFreq); sQS.snap (p.sQ); sGainS.snap (p.sGainDb);
        svf.prepare (fs, ch);
        updateCoeffs();
        recomputePending = false;
        initialized = false;
        wasActive = false;
        reset();
    }

    void reset() noexcept
    {
        for (int c = 0; c < kMaxChannels; ++c)
            for (int s = 0; s < kMaxSections; ++s) bq[s][c].reset();
        for (int s = 0; s < kMaxSections; ++s) bqSide[s].reset();
        svf.reset();
    }

    // Call from the SAME thread as processBlock() (typically the audio thread, where the host
    // adapter reads its atomic parameters), or synchronise externally — the engine takes no
    // internal lock. The FIRST call SNAPS (a freshly loaded plugin starts at its settings, no
    // ramp); later calls smooth freq/Q/gain. Identical params are a no-op (recompute stays skipped
    // even if called every block). type/slope/swept/on apply on the next block.
    void setParams (const BandParams& npIn) noexcept
    {
        auto finiteOr = [] (double x, double fb) noexcept { return std::isfinite (x) ? x : fb; };
        BandParams np = npIn;
        np.freq    = std::clamp (finiteOr (np.freq, 1000.0), 10.0, 0.49 * fs);
        np.Q       = std::clamp (finiteOr (np.Q, 1.0), 0.05, 40.0);
        np.gainDb  = std::clamp (finiteOr (np.gainDb, 0.0), -30.0, 30.0);
        np.sFreq   = std::clamp (finiteOr (np.sFreq, 1000.0), 10.0, 0.49 * fs);
        np.sQ      = std::clamp (finiteOr (np.sQ, 1.0), 0.05, 40.0);
        np.sGainDb = std::clamp (finiteOr (np.sGainDb, 0.0), -30.0, 30.0);

        if (! initialized)
        {
            p = np;
            freqS.snap (p.freq);   qS.snap (p.Q);   gainS.snap (p.gainDb);
            sFreqS.snap (p.sFreq); sQS.snap (p.sQ); sGainS.snap (p.sGainDb);
            initialized = true;
            recomputePending = true;
        }
        else if (! (np == p))
        {
            p = np;
            freqS.setTarget (p.freq);   qS.setTarget (p.Q);   gainS.setTarget (p.gainDb);
            sFreqS.setTarget (p.sFreq); sQS.setTarget (p.sQ); sGainS.setTarget (p.sGainDb);
            recomputePending = true;
        }
    }

    const BandParams& params() const noexcept { return p; }

    // Audio thread. In-place. RT-safe (no alloc / lock / IO).
    void processBlock (float* const* channels, int numChannels, int numSamples) noexcept
    {
        const bool midActive  = p.on && ! p.bypass;                       // Mid/main lane
        const bool sideActive = p.on && p.ms && p.sOn && ! p.sBypass;     // Side lane (M/S only)
        const bool active     = midActive || sideActive;                  // a bypassed Mid still lets Side run
        if (! active || numSamples <= 0)
        {
            if (! active && wasActive) { reset(); wasActive = false; }   // clear the tail so re-enabling doesn't pop
            return;
        }
        wasActive = true;

        freqS.advance (numSamples);  qS.advance (numSamples);  gainS.advance (numSamples);
        sFreqS.advance (numSamples); sQS.advance (numSamples); sGainS.advance (numSamples);

        // Recompute only when something actually moves — a static, settled band skips the
        // per-block tan/exp/sin/sqrt entirely.
        const bool moving = ! (freqS.settled() && qS.settled() && gainS.settled()
                            && sFreqS.settled() && sQS.settled() && sGainS.settled());
        if (recomputePending || moving) { updateCoeffs(); recomputePending = moving; }

        const int nc = numChannels < ch ? numChannels : ch;

        // M/S dual-design (2-ch only): exact local encode -> filter Mid (col 0) / Side (col 1) -> decode.
        // Matched filters only (no swept in M/S). Composes in series — the next band sees L/R again.
        if (p.ms && nc == 2)
        {
            float* L = channels[0];
            float* R = channels[1];
            for (int n = 0; n < numSamples; ++n)
            {
                const float m = 0.5f * (L[n] + R[n]);
                const float s = 0.5f * (L[n] - R[n]);
                float dM = 0.0f, dS = 0.0f;
                if (midActive)  { float x = m; for (int k = 0; k < designN;     ++k) x = bq[k][0].processSample (x); dM = x - m; }
                if (sideActive) { float y = s; for (int k = 0; k < designNside; ++k) y = bqSide[k].processSample (y); dS = y - s; }
                L[n] += dM + dS;   // L=M+S, R=M-S: fold deltas back. An idle lane (d=0) leaves its axis bit-exact.
                R[n] += dM - dS;
            }
            flushState();
            return;
        }

        // Non-M/S path runs the Mid/main lane only. If the Mid lane is off (e.g. M/S band on mono with
        // Mid bypassed) there is nothing to do here.
        if (! midActive) { flushState(); return; }

        // One sample through this band's filter chain on a given per-channel state column.
        auto filt = [this] (int col, float x) noexcept
        {
            if (sweptActive()) return svf.processSample (col, x);
            for (int s = 0; s < designN; ++s) x = bq[s][col].processSample (x);
            return x;
        };

        // Stereo / mono / surround: the Mid/main lane runs per channel (M/S is handled above).
        for (int c = 0; c < nc; ++c)
        {
            float* d = channels[c];
            for (int n = 0; n < numSamples; ++n) d[n] = filt (c, d[n]);
        }

        flushState();   // per-block denormal guard
    }

    // Complex frequency response at digital w (rad/sample) from the band's current smoothed state.
    // For a race-free GUI curve prefer the free bandResponse() with your own BandParams snapshot.
    std::complex<double> response (double w) const noexcept
    {
        if (! p.on || p.bypass) return { 1.0, 0.0 };
        std::complex<double> h { 1.0, 0.0 };
        for (int s = 0; s < designN; ++s) h *= evalCoeffs (coeffs[s], w);
        return h;
    }

private:
    void updateCoeffs() noexcept
    {
        BandParams sp = p;
        sp.freq = freqS.value(); sp.Q = qS.value(); sp.gainDb = gainS.value();

        const BandDesign d = designBand (sp, fs);

        // Side lane (M/S only): design from a side-view of the params (s* fields into the design slots).
        BandDesign dS; dS.n = 0;
        if (p.ms)
        {
            BandParams ssp = p;
            ssp.type = p.sType; ssp.slope = p.sSlope; ssp.swept = false;
            ssp.freq = sFreqS.value(); ssp.Q = sQS.value(); ssp.gainDb = sGainS.value();
            dS = designBand (ssp, fs);
        }

        // A topology switch (section count changed, swept<->static, route, OR Stereo<->M/S) changes
        // WHICH filters run — clear all state so a re-activated section never resumes a stale tail.
        const bool sw = sweptActive();
        if (d.n != designN || dS.n != designNside || sw != lastSwept || p.ms != lastMs) reset();

        designN     = d.n;
        designNside = dS.n;
        lastSwept = sw;
        lastMs    = p.ms;

        for (int s = 0; s < d.n; ++s) coeffs[s] = d.sec[s];
        for (int c = 0; c < ch; ++c)
            for (int s = 0; s < d.n; ++s) bq[s][c].setCoeffs (coeffs[s]);

        for (int s = 0; s < dS.n; ++s) { coeffsSide[s] = dS.sec[s]; bqSide[s].setCoeffs (coeffsSide[s]); }   // Side -> its OWN state

        if (sw) svf.setParams (p.type, sp.freq, sp.Q, sp.gainDb);   // swept = single SVF stage (12 dB/oct)
    }

    void flushState() noexcept
    {
        if (sweptActive()) { svf.flushDenormals(); return; }
        for (int c = 0; c < ch; ++c)                                   // Mid/main lane, per channel
            for (int s = 0; s < kMaxSections; ++s) bq[s][c].flushDenormals();
        for (int s = 0; s < kMaxSections; ++s) bqSide[s].flushDenormals();   // Side lane (M/S)
    }

    // The swept (SVF) engine only runs for types the single SVF stage can actually realise. Tilt has
    // no one-SVF realisation with a unity pivot beyond ~6 dB of tilt, so a swept Tilt band runs the
    // matched two-shelf design instead (the audio finally matches the displayed curve; previously a
    // swept Tilt was a silent pass-through while the GUI showed the full ±gain tilt).
    bool sweptActive() const noexcept { return p.swept && p.type != FilterType::Tilt; }

    BandParams p;
    double fs = 44100.0;
    int    ch = 2;
    int    designN = 1;
    int    designNside = 0;
    bool   recomputePending = true;
    bool   initialized = false;
    bool   wasActive = false;
    bool   lastSwept = false;
    bool   lastMs = false;

    Smoother freqS, qS, gainS, sFreqS, sQS, sGainS;
    Biquad   bq[kMaxSections][kMaxChannels];
    Biquad   bqSide[kMaxSections];   // Side lane's OWN state (M/S runs it on 2-ch only). It must never
                                     // alias a Mid channel column: writing Side coeffs into bq[s][1]
                                     // corrupted channel 1 on surround (nc>=3), where the plain
                                     // per-channel loop — not the M/S branch — consumes column 1.
    Svf      svf;
    BiquadCoeffs coeffs[kMaxSections], coeffsSide[kMaxSections];
};

} // namespace felitronics::eq

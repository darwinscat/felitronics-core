// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/core/Math.h>   // felitronics::core::kPi

#include <algorithm>
#include <cmath>

namespace felitronics::dynamics
{

//==============================================================================
// felitronics::dynamics::BandBallistics — attack/release times DERIVED FROM THE BAND ITSELF, so a
// dynamic EQ band can be "auto" (the user sets a range and nothing else) and still time itself
// correctly from a de-esser at 7 kHz to a boom tamer at 120 Hz.
//
// The law is the band's RINGING time, not its period. A 2nd-order bandpass has poles at real part
// -alpha, alpha = w0/(2Q); driven at fc its envelope rises as (1 - e^-alpha*t), so the 10-90% rise is
// ln(9)/alpha = ln9*Q/(pi*fc) — LINEAR IN Q/fc, i.e. the inverse bandwidth. An earlier design used
// "three periods, log-widened for narrow bands" and was 5.6x too fast at Q=40: the follower would
// modulate the gain faster than the filter can ring, which is intermodulation on every narrow band.
//
// NEAR NYQUIST the analog form is wrong. The TPT/SVF probe is bilinear-prewarped in fc only
// (eq::Svf, g = tan(pi*f/fs)), so its -3 dB edges compress toward Nyquist and the EFFECTIVE Q rises:
// at fc=16k/Q=40 the analog form under-predicts by more than 2x. We therefore use the DIGITAL
// bandwidth, Bw = fs*sin(2*pi*fc/fs)/(2*pi*Q), which degenerates to the analog form as fc/fs -> 0:
//
//     ringSec = 2*ln9*Q / (fs * sin (2*pi*fc/fs))
//
// Two corrections on top of the raw ring time:
//   • PERIOD FLOOR. Below Q = pi/ln9 ~ 1.43 the ring is shorter than one cycle, and a sub-period
//     "attack" on an envelope is undefined — it modulates gain INSIDE the cycle (the classic
//     low-frequency compressor distortion). Never go below a few periods.
//   • ATTACK FRACTION. The sidechain bandpass ALREADY imposes that same rise on the probe. Setting
//     the gain follower equal to it stacks two matched lags (~1.4x slower than either), so to react
//     as fast as the band physically allows the follower takes a FRACTION of the ring time.
//
// Release is a small multiple of attack (the ring-DOWN is the same time constant — symmetric poles),
// NOT a large one: at 12x, release reaches seconds and collides with the program-level estimator's
// own time constant, and the two servos then chase each other and the gain reduction breathes. Keep
// the release well clear of a programme-level estimator's own constant (RelativeLevel defaults to a
// 2 s symmetric average against a 500 ms release ceiling). The topology is feed-forward — the
// estimator watches the unprocessed probe — so a collision costs breathing, not stability.
//
// UNITS — the trap worth naming. ringMs() reports a 10-90% RISE TIME, because that is what can be
// measured against a real filter (and the tests do exactly that). GainReductionFollower, like every
// follower here, reads its milliseconds as the exponential TIME CONSTANT tau. compute() converts
// (tau = rise / ln9) before applying attackFraction; skipping that conversion silently makes the
// follower 2.197x slower than the fraction claims.
//
// Q RAILS: clamped to [0.05, 100], which contains what eq::EqBand admits ([0.05, 40]) but is WIDER
// than eq::Svf's own guard (it only rejects Q < 1e-3). A consumer that drives an Svf outside this
// range directly would get ballistics describing a different filter than the one running, so
// consumers should sanitise fc/Q once, centrally, and hand the same values to both.
//
// Pure and stateless — a calculator, not a processor. Non-finite / out-of-range inputs fall back to
// sane values rather than propagating; the result is always finite and inside the rails.
class BandBallistics
{
public:
    struct Config
    {
        double periodFloor    = 2.5;    // attack never shorter than this many periods of fc
        double attackFraction = 0.35;   // fraction of the probe's own ring time (see above)
        double releaseFactor  = 3.0;    // release = releaseFactor * attack
        // AUTO rails — where the automatic value is allowed to land.
        double attackMinMs    = 1.0,  attackMaxMs  = 300.0;
        double releaseMinMs   = 12.0, releaseMaxMs = 500.0;
        // ABSOLUTE rails — where the value may land AFTER the user's deviation knob. Wider on
        // purpose: clamping auto and knob to the same rails made the knob a placebo wherever auto
        // already sat on a rail (every treble band pinned to the same release), which is precisely
        // the de-esser case this module exists to serve.
        double attackAbsMinMs  = 0.2,  attackAbsMaxMs  = 1000.0;
        double releaseAbsMinMs = 3.0,  releaseAbsMaxMs = 1000.0;
        double knobRange      = 4.0;    // deviation knob spans 2^(+-knobRange/2) => 1/4x .. 4x
    };

    struct Times { double attackMs = 10.0, releaseMs = 100.0; };

    // The probe's own 10-90% envelope rise, in ms — digital-bandwidth form (exact near Nyquist).
    static double ringMs (double fs, double fc, double Q) noexcept
    {
        fs = finite (fs, 48000.0); if (fs <= 0.0) fs = 48000.0;
        fc = finite (fc, 1000.0);
        Q  = finite (Q, 1.0);
        // Rails wide enough to contain what the EQ actually admits (eq::EqBand clamps Q to [0.05, 40];
        // eq::Svf itself only guards Q >= 1e-3), so the ballistics never describe a band the filter
        // could not be running.
        fc = std::clamp (fc, 1.0, 0.49 * fs);
        Q  = std::clamp (Q, 0.05, 100.0);

        const double s = std::sin (2.0 * core::kPi * fc / fs);
        if (! (s > 1.0e-12)) return 1.0e3;                       // degenerate: treat as very slow
        return 1.0e3 * 2.0 * kLn9 * Q / (fs * s);
    }

    // knob in [0,1]: 0.5 == auto, 0 == 1/4x (faster), 1 == 4x (slower). Out-of-range/NaN -> auto.
    static double knobMultiplier (double knob, double knobRange = 4.0) noexcept
    {
        knob = finite (knob, 0.5);
        knob = std::clamp (knob, 0.0, 1.0);
        return std::pow (2.0, (knob - 0.5) * finite (knobRange, 4.0));
    }

    // The deviation multipliers are applied BEFORE the final clamp, so "faster" cannot push a treble
    // band below the floor. That dead zone is deliberate and documented, not a bug.
    static Times compute (double fs, double fc, double Q,
                          double attackKnob, double releaseKnob) noexcept
    {
        return compute (fs, fc, Q, attackKnob, releaseKnob, Config{});
    }

    static Times compute (double fs, double fc, double Q) noexcept
    {
        return compute (fs, fc, Q, 0.5, 0.5, Config{});
    }

    static Times compute (double fs, double fc, double Q,
                          double attackKnob, double releaseKnob,
                          const Config& cfgRaw) noexcept
    {
        // The rails come from the CALLER, so they get the same finite/ordering discipline as every
        // other input: std::clamp with lo > hi, or with a NaN bound, is undefined behaviour.
        const Config cfg = sanitize (cfgRaw);

        fs = finite (fs, 48000.0); if (fs <= 0.0) fs = 48000.0;
        const double fcHi      = std::max (1.0, 0.49 * fs);            // fs < ~2 Hz would invert these
        const double fcClamped = std::clamp (finite (fc, 1000.0), 1.0, fcHi);

        // UNITS. ringMs() is a 10-90% rise; GainReductionFollower reads its milliseconds as the
        // exponential TIME CONSTANT (tau). Converting here is not cosmetic: feeding a 10-90 figure
        // to a tau-based follower makes it ln9 = 2.197x slower than the fraction claims.
        const double tauBandMs = ringMs (fs, fc, Q) / kLn9;
        const double period    = cfg.periodFloor * 1.0e3 / fcClamped;

        const double atkAuto = std::clamp (std::max (cfg.attackFraction * tauBandMs, period),
                                           cfg.attackMinMs, cfg.attackMaxMs);
        const double relAuto = std::clamp (atkAuto * cfg.releaseFactor,
                                           cfg.releaseMinMs, cfg.releaseMaxMs);

        double atk = std::clamp (atkAuto * knobMultiplier (attackKnob, cfg.knobRange),
                                 cfg.attackAbsMinMs, cfg.attackAbsMaxMs);
        double rel = std::clamp (relAuto * knobMultiplier (releaseKnob, cfg.knobRange),
                                 cfg.releaseAbsMinMs, cfg.releaseAbsMaxMs);
        rel = std::max (rel, atk);   // release is never shorter than attack, knobs notwithstanding

        return { atk, rel };
    }

private:
    static constexpr double kLn9 = 2.1972245773362196;   // ln(9) — the 10->90% factor

    static double finite (double v, double fallback) noexcept { return std::isfinite (v) ? v : fallback; }

    static Config sanitize (const Config& in) noexcept
    {
        Config c;
        c.periodFloor    = std::max (0.0, finite (in.periodFloor, 2.5));
        c.attackFraction = std::max (0.0, finite (in.attackFraction, 0.35));
        c.releaseFactor  = std::max (1.0, finite (in.releaseFactor, 3.0));
        c.attackMinMs      = std::max (1.0e-3, finite (in.attackMinMs, 1.0));
        c.attackMaxMs      = std::max (c.attackMinMs, finite (in.attackMaxMs, 300.0));
        c.releaseMinMs     = std::max (1.0e-3, finite (in.releaseMinMs, 12.0));
        c.releaseMaxMs     = std::max (c.releaseMinMs, finite (in.releaseMaxMs, 500.0));
        c.attackAbsMinMs   = std::min (c.attackMinMs,  std::max (1.0e-3, finite (in.attackAbsMinMs, 0.2)));
        c.attackAbsMaxMs   = std::max (c.attackMaxMs,  finite (in.attackAbsMaxMs, 1000.0));
        c.releaseAbsMinMs  = std::min (c.releaseMinMs, std::max (1.0e-3, finite (in.releaseAbsMinMs, 3.0)));
        c.releaseAbsMaxMs  = std::max (c.releaseMaxMs, finite (in.releaseAbsMaxMs, 1000.0));
        c.knobRange        = std::max (0.0, finite (in.knobRange, 4.0));
        return c;
    }
};

} // namespace felitronics::dynamics

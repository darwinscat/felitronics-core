// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/core/FlushToZero.h>
#include <felitronics/dynamics/ChannelLinker.h>
#include <felitronics/dynamics/EnvelopeFollower.h>

#include <algorithm>
#include <cmath>
#include <vector>

namespace felitronics::dynamics
{

//==============================================================================
// felitronics::dynamics::NoiseGate — a "dual-detection" (ISP Decimator "G-String" style) noise gate.
//
// Architecturally DISTINCT from the module's continuous Compressor / GainComputer(DownExpand): a gate is
// a bistable Schmitt trigger — hysteresis (open ≠ close) + a HOLD timer bridge the rectified ripple of low
// notes, and the gain slews with a LINEAR-fast OPEN (preserves the pick attack — a dB-linear open from a
// deep floor sits near-silent for most of the ramp and swallows the transient) and an EXP-slow CLOSE.
//
// It COMPOSES the module kit: the linked key uses ChannelLinker (max|ch| over ≤2 lanes) on a per-channel
// sidechain HPF, and the detector envelope is a Peak EnvelopeFollower. On top sit the gate-specific parts
// the compressor primitives don't have: the Schmitt state machine + hold, the asymmetric VCA ramp, the
// closed floor, and an ENABLE crossfade so toggling the gate on/off never pops.
//
// TWO-PHASE, to support a KEYED gate where extra stages run between detect and attenuate (the OrbitCab
// case: detector on the clean guitar, VCA after the amp EQ, with the preamp + EQ in between):
//   • analyse() reads the KEY and fills a per-sample gain curve;
//   • applyGain() multiplies a possibly-DIFFERENT downstream buffer by that curve.
// A latency gap between the two acts as free lookahead (uncompensated). process() is the self-keyed
// convenience (analyse + applyGain on the same buffer).
//
// 🔴 RT-safe: analyse()/applyGain()/process() never allocate, lock, do IO, or throw — all storage is sized
// once in prepare(). NaN/Inf-safe and denormal-safe in software (Law 8), no host FTZ assumed. Zero latency.
//==============================================================================
class NoiseGate
{
public:
    // Fixed voicing. Only the THRESHOLD is dynamic (an analyse() argument) — the rest shapes the "feel"
    // and rarely needs touching. Defaults = the values dialled in and shipped by OrbitCab's in-amp gate.
    struct Config
    {
        float hysteresisDb  = 6.0f;    // open threshold − close threshold (Schmitt)
        float envAttackMs   = 0.5f;    // detector attack — catch the pick
        float envReleaseMs  = 20.0f;   // detector release — with the hold, bridges low-note rectified ripple
        float holdMs        = 40.0f;   // stay-open time after the level drops below close (the real chatter bridge)
        float openMs        = 0.3f;    // VCA LINEAR-in-amplitude open (transient-safe)
        float closeMs       = 100.0f;  // VCA exponential close
        float floorDb       = -90.0f;  // closed attenuation (deep, but never exact zero → no denormal / conv-tail chill)
        float enableMs      = 25.0f;   // on/off crossfade (feature toggle / host bypass) — no pop
        float sidechainHpHz = 75.0f;   // sidechain rumble / DC reject
        float keyClampAbs   = 16.0f;   // detector-input bound (+24 dBFS): a pathological ±FLT_MAX can't ring the HPF/env forever
        LinkMode link       = LinkMode::Max;   // image-preserving multichannel link (reacts to the loudest lane)
    };

    void prepare (double sampleRate, int maxBlock, int /*maxChannels*/)
    {
        sr = sampleRate > 1000.0 ? sampleRate : 48000.0;
        gainCurve.assign ((std::size_t) std::max (1, maxBlock), 1.0f);
        setConfig (cfg);
        // enable is seeded by the first analyse() block (on ? 1 : 0 via the ramp); start unity.
        enable = 0.0f;
        reset();
    }

    // Recompute coefficients from `c`. Not the audio thread (or call it before the stream starts).
    void setConfig (const Config& c)
    {
        cfg = c;
        env.prepare (sr);
        env.setDetector (Detector::Peak);
        env.setTimes (cfg.envAttackMs, cfg.envReleaseMs);
        holdSamples = (int) std::lround ((double) cfg.holdMs * 0.001 * sr);
        hpR        = (float) std::exp (-6.283185307179586 * (double) std::max (0.0f, cfg.sidechainHpHz) / sr);
        floorGain  = (float) std::pow (10.0, (double) cfg.floorDb / 20.0);
        openStep   = (float) ((1.0 - (double) floorGain) / std::max (1.0, (double) cfg.openMs  * 0.001 * sr));
        closeRatio = (float) std::pow ((double) floorGain, 1.0 / std::max (1.0, (double) cfg.closeMs * 0.001 * sr));
        enableStep = (float) (1.0 / std::max (1.0, (double) cfg.enableMs * 0.001 * sr));
    }

    // Stream restart: clear the detector + close the gate. `enable` is NOT cleared here — it's the on/off
    // crossfade position (a parameter state), which a stream restart must not jump.
    void reset()
    {
        env.reset();
        xPrev[0] = xPrev[1] = 0.0f;
        hpPrev[0] = hpPrev[1] = 0.0f;
        coreGain = floorGain;   // start CLOSED
        lastCoreGain = floorGain;
        lastGain = 1.0f + enable * (floorGain - 1.0f);
        open = false;
        hold = 0;
    }

    // Seed the on/off crossfade to a known state (a JUMP, no ramp): call after prepare() when the host knows
    // the restored on-state, so a session saved gate-ON starts already gated instead of a 25 ms fade-in from
    // unity (which would leak a block of ungated signal on load).
    void seedEnabled (bool on) noexcept
    {
        enable = on ? 1.0f : 0.0f;
        lastGain = 1.0f + enable * (lastCoreGain - 1.0f);
    }

    // PHASE A — read the LINKED key from `key` (≤2 lanes are linked) and fill the per-sample gain curve.
    // `on` = the gate feature toggle; `thresholdDb` = the open threshold (dBFS, vs the key level).
    void analyse (const float* const* key, int numChannels, int n, bool on, float thresholdDb) noexcept
    {
        n = std::min (n, (int) gainCurve.size());
        const int keyCh = std::min (numChannels, 2);
        if (keyCh < 1 || n <= 0) return;                                  // no lane to key off — never deref key[0] blindly
        if (! std::isfinite (thresholdDb)) thresholdDb = -50.0f;          // a NaN threshold would freeze both compares
        const float openLin  = dbToGain (thresholdDb);
        const float closeLin = dbToGain (thresholdDb - cfg.hysteresisDb);
        const float enTarget = on ? 1.0f : 0.0f;

        // per-sample HPF outputs, wrapped as 1-sample "channels" so ChannelLinker::linkAmplitude applies.
        float hp[2] = { 0.0f, 0.0f };
        const float* hpPlanes[2] = { &hp[0], &hp[1] };

        for (int i = 0; i < n; ++i)
        {
            for (int ch = 0; ch < keyCh; ++ch)
            {
                float x = key[ch][i];
                if (! std::isfinite (x)) x = 0.0f;                        // heal a NaN/±Inf key sample
                x = std::clamp (x, -cfg.keyClampAbs, cfg.keyClampAbs);    // bound the detector (audio path untouched)
                float h = x - xPrev[ch] + hpR * hpPrev[ch];              // one-pole sidechain high-pass
                if (! std::isfinite (h)) h = 0.0f;
                xPrev[ch] = x; hpPrev[ch] = h;
                core::flushDenormal (xPrev[ch]);
                core::flushDenormal (hpPrev[ch]);
                hp[ch] = h;
            }
            const float r = linkAmplitude (cfg.link, hpPlanes, keyCh, 0);  // max|hp| (image-preserving link)
            const float e = env.process (r);                              // Peak follower: fast attack, slow release

            // Schmitt state machine: hysteresis (open ≠ close) + hold, in the LINEAR domain.
            if (open) { if (e < closeLin) { if (--hold <= 0) open = false; } else hold = holdSamples; }
            else      { if (e > openLin)  { open = true; hold = holdSamples; } }
            const float target = open ? 1.0f : floorGain;

            // VCA: LINEAR-fast OPEN (transient-safe), EXP-slow CLOSE.
            if      (coreGain < target) coreGain = std::min (coreGain + openStep, target);
            else if (coreGain > target) coreGain = std::max (coreGain * closeRatio, target);

            // ENABLE crossfade: gate off ⇒ unity passthrough, gate on ⇒ gated; ramped so a toggle never pops.
            enable += std::clamp (enTarget - enable, -enableStep, enableStep);
            gainCurve[(std::size_t) i] = 1.0f + enable * (coreGain - 1.0f);
        }
        env.flushDenormals();
        lastCoreGain = coreGain;
        lastGain     = gainCurve[(std::size_t) (n - 1)];
    }

    // PHASE B — apply the stored curve to `io` (the SAME curve on every lane = a linked gate).
    void applyGain (float* const* io, int numChannels, int n) const noexcept
    {
        n = std::min (n, (int) gainCurve.size());
        for (int ch = 0; ch < numChannels; ++ch)
            for (int i = 0; i < n; ++i)
                io[ch][i] *= gainCurve[(std::size_t) i];
    }

    // Self-keyed convenience: detect + attenuate the same buffer in one call.
    void process (float* const* io, int numChannels, int n, bool on, float thresholdDb) noexcept
    {
        analyse (io, numChannels, n, on, thresholdDb);
        applyGain (io, numChannels, n);
    }

    // EFFECTIVE applied gain at the last block end (incl. the enable crossfade) — the honest GR-meter reading:
    // 1 = no reduction (open OR disabled), floorGain = fully closed & enabled.
    float currentGain()     const noexcept { return lastGain; }
    // Raw VCA gain (open/close only, before the enable crossfade) — for state-machine tests / diagnostics.
    float currentCoreGain() const noexcept { return lastCoreGain; }
    float floorGainLinear() const noexcept { return floorGain; }
    int   latencySamples()  const noexcept { return 0; }

private:
    static float dbToGain (float db) noexcept { return std::pow (10.0f, db * 0.05f); }

    Config cfg {};
    double sr = 48000.0;
    float  hpR = 0.0f, floorGain = 0.0f, openStep = 0.0f, closeRatio = 0.0f, enableStep = 0.0f;
    int    holdSamples = 0;

    EnvelopeFollower env;                 // reused: the Peak detector ballistics + Law-8 flush
    float  xPrev[2]  { 0.0f, 0.0f };      // per-channel one-pole sidechain HPF state
    float  hpPrev[2] { 0.0f, 0.0f };
    float  coreGain = 0.0f;               // gate VCA gain (open/close), before the enable crossfade
    float  lastCoreGain = 0.0f;
    float  lastGain = 1.0f;               // effective (enable-weighted) gain at last block end
    float  enable = 0.0f;                 // 0 = off (unity), 1 = on (gated) — crossfaded, survives reset()
    bool   open = false;
    int    hold = 0;

    std::vector<float> gainCurve;         // per-sample multiplier (incl. enable), sized maxBlock in prepare()
};

} // namespace felitronics::dynamics

// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/core/Config.h>
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

    // Law 11: every argument prepare() takes is BINDING. `maxChannels` used to be ignored here, which
    // made the gate's declared width a fiction — it would happily attenuate any number of lanes, and the
    // fused process() below needs a bound to build its chunk pointers from. Refuses rather than clamps,
    // for the reason Compressor states: a silently reduced width refuses every later call instead.
    [[nodiscard]] bool prepare (double sampleRate, int maxBlock, int maxChannels)
    {
        prepared_ = false;                                   // any early return below leaves it unprepared
        if (maxBlock < 1) return false;
        if (maxChannels < 1 || maxChannels > core::kMaxChannels) return false;
        sr = sampleRate > 1000.0 ? sampleRate : 48000.0;
        maxCh_ = maxChannels;
        gainCurve.assign ((std::size_t) maxBlock, 1.0f);
        setConfig (cfg);
        // enable is seeded by the first analyse() block (on ? 1 : 0 via the ramp); start unity.
        enable = 0.0f;
        reset();
        prepared_ = true;
        return true;
    }

    bool isPrepared() const noexcept { return prepared_; }
    int  maxBlock()   const noexcept { return prepared_ ? (int) gainCurve.size() : 0; }
    int  maxChannels() const noexcept { return prepared_ ? maxCh_ : 0; }

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
        analysedN_ = 0;                   // a restarted stream carries no analysed curve
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
    // LAW 11(a), the named exception: this phase's RESULT is the curve, and the curve is exactly
    // maxBlock long, so a longer call cannot be chunked — the second phase would have nothing to apply.
    // It therefore REFUSES instead of clamping. Clamping is what it used to do, and it let the samples
    // past maxBlock out UNGATED: measured 3840 of 4096 at +89.99 dB over the gated ones, which is
    // 100% of the construction ceiling (-floorDb = 90 dB). Use process() for a call of any length.
    [[nodiscard]] bool analyse (const float* const* key, int numChannels, int n, bool on, float thresholdDb) noexcept
    {
        if (numChannels < 0 || n < 0) return false;                       // malformed
        if (! prepared_) return false;
        if (numChannels > maxCh_) return false;                           // width is a LIMIT
        if (n > (int) gainCurve.size()) return false;                     // capacity-bearing: see above
        if (n == 0) return true;                                          // law 11(d): n == 0 is the ONE true
                                                                          // no-op — it must not invalidate a
                                                                          // curve a previous call produced
        const int keyCh = std::min (numChannels, 2);
        analysedN_ = 0;                                                   // no valid curve until this call writes one
        if (keyCh < 1) return true;                                       // no lane to key off — never deref key[0] blindly
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
        analysedN_   = n;                                                 // this is how far the curve is valid
        return true;
    }

    // How many samples of the curve the last accepted analyse() produced. applyGain() may not reach past it.
    int analysedSamples() const noexcept { return analysedN_; }

    // PHASE B — apply the stored curve to `io` (the SAME curve on every lane = a linked gate).
    [[nodiscard]] bool applyGain (float* const* io, int numChannels, int n) const noexcept
    {
        if (numChannels < 0 || n < 0) return false;
        if (! prepared_) return false;
        if (numChannels > maxCh_) return false;
        // THE CURVE HAS A LENGTH, AND ONLY THE CURVE KNOWS IT. Bounding phase B by the buffer's CAPACITY
        // instead of by what phase A actually wrote let a short analysis be followed by a long apply, and
        // the tail was then multiplied by the PREVIOUS call's curve: measured, analyse(quiet, 512) settled
        // closed, then analyse(LOUD, 256), then applyGain(512) attenuated samples [256,512) by 90.0 dB —
        // the floor — on material that was never analysed at all. Same class as the width EqEngine's
        // captureSectionInput closed (a capture narrower than the read), one axis over.
        if (numChannels == 0) return true;   // no lanes to attenuate — the same geometry analyse()
                                             // accepts, and law 11's order answers it before the curve
        if (n > analysedN_) return false;
        for (int ch = 0; ch < numChannels; ++ch)
            for (int i = 0; i < n; ++i)
                io[ch][i] *= gainCurve[(std::size_t) i];
        return true;
    }

    // Self-keyed convenience: detect + attenuate the same buffer in one call. LAW 11(a) — the length is
    // a CAPACITY: this form chunks to maxBlock, so any n is gated IN FULL and the result is bit-identical
    // to the caller having chunked it itself (every piece of state — env, hold, open, coreGain, enable —
    // lives in members and carries across the seam). The two-phase form above cannot do this, which is
    // exactly why this one exists.
    [[nodiscard]] bool process (float* const* io, int numChannels, int n, bool on, float thresholdDb) noexcept
    {
        if (numChannels < 0 || n < 0) return false;
        if (! prepared_) return false;
        if (numChannels > maxCh_) return false;
        if (n == 0) return true;                        // no samples: the ONE true no-op
        if (numChannels == 0)
        {
            // An accepted call with no lanes produced no curve, and the fused form has to say so or it
            // disagrees with its own two-phase halves: analyse(key, 0, n) sets the length to 0, so a
            // later applyGain() is refused — while this path used to leave the OLD curve standing and
            // the same applyGain() then attenuated by the floor, about -90 dB.
            analysedN_ = 0;
            return true;
        }

        const int mb = (int) gainCurve.size();
        float* sub[core::kMaxChannels] {};
        for (int off = 0; off < n; )
        {
            const int m = std::min (n - off, mb);
            for (int c = 0; c < numChannels; ++c) sub[(std::size_t) c] = io[c] + off;
            // NEITHER OF THESE CAN FIRE, and that is load-bearing rather than lucky: every condition the
            // two phases refuse on — malformed extents, unprepared, a width past maxCh_, a length past the
            // curve — is settled above this loop, and `m <= mb` by construction. It matters because a
            // mid-loop `return false` would break law 11's own invariant: the chunks BEFORE it have
            // already run, so the call would be both refused and half-done. If a future guard makes one
            // of these reachable, the loop has to become all-or-nothing FIRST.
            if (! analyse (sub, numChannels, m, on, thresholdDb)) return false;
            if (! applyGain (sub, numChannels, m)) return false;
            off += m;                                   // `off += mb` could step past INT_MAX
        }
        return true;
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
    int    maxCh_ = 0;                    // declared width — law 11(b), a LIMIT, not a clamp
    int    analysedN_ = 0;                // samples of gainCurve the last analyse() actually wrote
    bool   prepared_ = false;             // true only after a prepare() that succeeded
};

} // namespace felitronics::dynamics

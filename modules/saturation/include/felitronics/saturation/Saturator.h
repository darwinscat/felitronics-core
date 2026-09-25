// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/saturation/WaveShaper.h>
#include <felitronics/oversampling/Oversampler.h>
#include <felitronics/oversampling/PolyphaseOversampler.h>
#include <felitronics/core/Config.h>
#include <felitronics/core/DelayLine.h>
#include <felitronics/core/Math.h>
#include <felitronics/core/StateGrid.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace felitronics::saturation
{

//==============================================================================
// felitronics::saturation::Saturator — the production soft-saturation stage (the mastering "glue/color"
// before the make-loud gain). Wraps a stateless WaveShaper curve in OVERSAMPLING (default 4×) so the new
// harmonics fold above the base Nyquist instead of aliasing back, runs a DC blocker inside the oversampled
// region (the asymmetric curve's even harmonics shift the mean), then returns to base rate and applies
// DRIVE-COMPENSATION + a linear dry/wet + output trim.
//
// Gain-staging (a reference tool reverted its saturator twice over this): the curve
// is peak-normalised (|x|≤1 → |y|≤1), `autoComp` undoes the small-signal-gain bump so loudness doesn't jump
// with drive, and the dry/wet blend is LINEAR (convex combination of bounded signals → peak-safe). Place it
// BEFORE the loudness gain + final limiter.
//
// RT-safe: prepare() does all allocation; process() is alloc/lock/throw-free, in place. n may exceed the
// maxBlock passed to prepare() — process() chunks internally, state carries across chunks. With oversampling
// the stage reports a round-trip latency.
//
// PARAMETERS GLIDE — and the first write of a stream SNAPS. See `kGlideMs` for the rule, the clock and the
// number it was chosen by. In one line: a write that moves drive, bias, auto-compensation, mix or output trim
// starts a linear parameter-space ramp on the 64-sample `core::StateGrid`, and the curve's coefficients are
// interpolated per sample between two designed sets, so no knob steps the output; a write before the first
// sample after prepare()/reset() applies at once, exactly as setParams() always did. `shape` and `dcBlockHz`
// are not continuous and land at once, as they always did (a shape change snaps every parameter with it).
//
// THE OVERSAMPLER IS A CHOICE (P31), made at prepare(): `Topology::Kaiser` (the default, and the only one
// for a factor that is not a power of two) is PolyphaseOversampler — cutoff fixed at 0.45 fs, round trip
// tapsPerPhase - 1, and at 44.1 kHz -1.80 dB at 19 kHz and -15.55 at 20 kHz. `Topology::Cascade` is
// CascadeOversampler — just as strict, flat to 20 kHz at every rate, lengths from the sample rate, round
// trip 131 samples at 44.1 kHz and 76 at 48 kHz (4x). The dry path is delayed by whichever round trip was
// built, so a mix < 1 stays comb-free under both. Under Cascade `tapsPerPhase` is range-checked and
// otherwise unused, and the RATE is binding: the cascade designs itself from it and refuses one outside
// [1 kHz, 3 MHz] — a window the Kaiser stage, which never looked at the rate, does not have.
//
// Poison-hardened: non-finite params fall back to defaults (applyParams), each input sample is
// sanitized at the gate (NaN/Inf → 0, huge finite → ±1e6 clamp) so one bad sample can't lodge in the
// oversampler/DC/dry-delay state, and the DC-blocker state flushes non-finite values per block. All
// three guards are bit-transparent on finite, in-range signals (the NULL test proves it).
class Saturator
{
public:
    struct Params
    {
        WaveShaper::Shape shape = WaveShaper::Shape::Tanh;
        float driveDb   = 3.0f;    // 0..36 (mastering 1..6); k = dbToGain(driveDb) − 1
        float bias      = 0.0f;    // −0.5..0.5 — Asym only, sets the even-harmonic content
        float mix       = 1.0f;    // 0..1 dry/wet
        float outputDb  = 0.0f;    // −24..24 post trim
        float autoComp  = 0.5f;    // 0..1 — drive-compensation amount (1 = unity small-signal gain)
        float dcBlockHz = 10.0f;   // DC blocker corner (in the oversampled domain)
    };

    // WHAT prepare() ASKS THE HEAP FOR (law 11d) — the one function it sizes itself with, so a caller
    // budgeting memory reads the counts the six buffers and the dry bank are actually built from. FALSE,
    // with `out` untouched, exactly where prepare() refuses the same arguments (it IS prepare()'s gate),
    // and a refused prepare() allocates nothing. Asked of a FRESH saturator.
    struct Storage
    {
        std::size_t osBuf = 0, wetBuf = 0;      // floats
        std::size_t ptrs  = 0;                  // float* in EACH of the two pointer tables
        std::size_t dc    = 0;                  // floats in EACH of the two DC-blocker state vectors
        std::size_t dryLines = 0;               // core::DelayLine per channel
        int         dryDelaySamples = 0;        // the oversampler round trip each of them holds
        oversampling::Oversampler::Storage os {};            // empty when the factor is 1
        std::uint64_t bytes() const noexcept
        {
            return (std::uint64_t) sizeof (float)  * ((std::uint64_t) osBuf + wetBuf + 2u * (std::uint64_t) dc)
                 + (std::uint64_t) sizeof (float*) * (2u * (std::uint64_t) ptrs)
                 + core::delayBankBytes (dryLines, dryDelaySamples)
                 + os.bytes();
        }
    };

    [[nodiscard]] static bool storageFor (double sampleRate, int maxBlock, int maxChannels,
                                          int oversampleFactor, int tapsPerPhase, Storage& out,
                                          oversampling::Topology topology = oversampling::Topology::Kaiser) noexcept
    {
        // Spelled positively so a NaN FAILS — see prepare() below for what a NaN rate used to produce.
        if (! (sampleRate > 0.0) || ! std::isfinite (sampleRate) || maxBlock < 1) return false;
        if (maxChannels < 1 || maxChannels > core::kMaxChannels) return false;   // law 11(b): BINDING
        const int os = (oversampleFactor >= 2) ? oversampleFactor : 1;
        Storage st;
        if (os > 1 && ! oversampling::Oversampler::storageFor (topology, sampleRate, os, maxChannels, tapsPerPhase, st.os))
            return false;
        const std::size_t ch = (std::size_t) maxChannels;
        st.osBuf    = ch * (std::size_t) maxBlock * (std::size_t) os;
        st.wetBuf   = ch * (std::size_t) maxBlock;
        st.ptrs     = ch;
        st.dc       = ch;
        st.dryLines = ch;
        st.dryDelaySamples = os > 1 ? oversampling::Oversampler::latencyFor (topology, sampleRate, os, tapsPerPhase) : 0;
        out = st;
        return true;
    }

    // THE LATENCY a prepared saturator will report for this geometry, without preparing one — the dry
    // aligner of a composite is sized by it, and a second derivation of it is the drift law 11d's budgets
    // exist to make impossible. 0 where prepare() refuses the same arguments.
    //
    // THAT IS NOT THE SAME AS WHAT `latencySamples()` READS AFTER A REFUSAL, and the difference is this
    // class's, not this function's: `latencySamples()` is not gated by `prepared_` (unlike the compressor's
    // and the limiter's), so after a refused preparation it still reports the topology of the last one that
    // succeeded — 63 where this answers 0. A composite sizes itself from THIS, which is the number the
    // preparation it is about to make will produce.
    static int latencyFor (double sampleRate, int maxBlock, int maxChannels,
                           int oversampleFactor, int tapsPerPhase,
                           oversampling::Topology topology = oversampling::Topology::Kaiser) noexcept
    {
        Storage st;
        return storageFor (sampleRate, maxBlock, maxChannels, oversampleFactor, tapsPerPhase, st, topology)
             ? st.dryDelaySamples : 0;
    }

    bool prepare (double sampleRate, int maxBlock, int maxChannels, int oversampleFactor = 4,
                  int tapsPerPhase = oversampling::PolyphaseOversampler::kDefaultTapsPerPhase,
                  oversampling::Topology topology = oversampling::Topology::Kaiser)
    {
        prepared_ = false;                                             // any early return below leaves it unprepared
        // Spelled positively so a NaN FAILS. `NaN <= 0.0` is false, so the previous form accepted a
        // NaN rate, returned true, and went on to emit NaN: fsOs is NaN, so the DC blocker's
        // exp(-2π·fc/NaN) is NaN and the Asym path carries it into the output (measured on
        // prepare(NaN, 64, 1, 4, 32) — returns TRUE, output sample 40 is nan). Same shape and same
        // reason as the guard in TruePeakLimiter::prepare. The gate is storageFor()'s now, so the budget
        // and the preparation refuse the same arguments by construction rather than by agreement.
        //
        // ONE OBSERVABLE CHANGE COMES WITH THAT: `os_` used to be written BEFORE the oversampler was asked,
        // so a preparation the oversampler refused (`tapsPerPhase` past 1024) left the new factor standing
        // beside the OLD oversampler, and `latencySamples()` — which is not gated by `prepared_` — reported
        // the round trip of a topology that was never built. Measured, `prepare(4, 64)` then `prepare(1, 64)`
        // then a refused `prepare(4, 2000)`: 63 before, 0 now. The new answer is law 11(b)'s (a refused call
        // touches nothing); it is stated here because this header is shared with the plug-ins.
        Storage st;
        if (! storageFor (sampleRate, maxBlock, maxChannels, oversampleFactor, tapsPerPhase, st, topology)) return false;
        fs_       = sampleRate;
        maxBlock_ = maxBlock;
        channels_ = maxChannels;
        os_       = (oversampleFactor >= 2) ? oversampleFactor : 1;
        if (os_ > 1 && ! ovs_.prepare (topology, sampleRate, os_, channels_, tapsPerPhase)) return false;
        // A factor of 1 builds no oversampler, so it holds none: whatever an earlier preparation built (the
        // Kaiser buffers, or the heap-held cascade) is released here rather than kept beside a budget that
        // says the oversampler half is empty. Assigning a fresh switch only frees.
        if (os_ == 1) ovs_ = oversampling::Oversampler {};

        osBuf_.assign  (st.osBuf,  0.0f);
        wetBuf_.assign (st.wetBuf, 0.0f);
        osPtrs_.assign (st.ptrs, nullptr);
        wetPtrs_.assign(st.ptrs, nullptr);
        dcX1_.assign   (st.dc, 0.0f);
        dcY1_.assign   (st.dc, 0.0f);
        const int lat = st.dryDelaySamples;                      // align the dry to the wet's round-trip
        core::prepareDelayBank (dryDelay_, st.dryLines, lat);
        for (auto& d : dryDelay_) d.setDelay (lat);
        applyParams();
        // The ledger must die with the buffers it indexes. prepare() REALLOCATES every per-channel vector
        // above, so a stale count from a wider previous life would send the next drop past the end of the
        // new ones — measured as an AddressSanitizer container-overflow on dryDelay_ after
        // prepare(2) -> process -> prepare(1) -> process. prepare() does not call reset(), so this cannot
        // be left to reset() to do.
        ranNc_ = ranDcNc_ = 0;
        // The glide restarts with the stream, in the rate's own units: kGlideMs in 64-sample grid periods.
        glideTicks_ = ticksFor (sampleRate, glideMs_);
        snapGlide();
        grid_.reset();
        fresh_ = true;
        prepared_ = true;                                              // fully built — process() may now run
        return true;
    }

    void reset() noexcept
    {
        if (os_ > 1) ovs_.reset();
        std::fill (dcX1_.begin(), dcX1_.end(), 0.0f);
        std::fill (dcY1_.begin(), dcY1_.end(), 0.0f);
        for (auto& d : dryDelay_) d.reset();
        ranNc_ = ranDcNc_ = 0;   // nothing has run, so nothing can be stopping (see dropStoppedCells)
        // A STREAM RESTART LANDS EVERY GLIDE on its target and re-anchors the grid, and the next write snaps:
        // a restart that resumed a ramp would make a second render of the same programme start somewhere else.
        if (glideActive()) applyParams();
        snapGlide();
        grid_.reset();
        fresh_ = true;
    }

    int  latencySamples() const noexcept { return os_ > 1 ? ovs_.latencySamples() : 0; }

    // The round trip a KAISER oversampler of this shape costs, in baseband samples — `PolyphaseOversampler`'s
    // own `tpp - 1`. Kaiser only: the topology-aware answer is latencyFor(), which the budget now reads.
    static int latencyForFactor (int oversampleFactor, int tapsPerPhase) noexcept
    {
        return oversampleFactor > 1 && tapsPerPhase > 0 ? tapsPerPhase - 1 : 0;
    }

    // THE GLIDE. A write that changes a CONTINUOUS parameter — driveDb, bias, autoComp, mix, outputDb — does not
    // apply at once: it becomes a target, and at the next `core::StateGrid` boundary (64 samples of audio time) a
    // linear ramp in PARAMETER space starts toward it, `glideTicks()` boundaries long. At each boundary the curve
    // is DESIGNED at the ramp's next point (the same WaveShaper + drive-compensation arithmetic as a settled
    // stage), and inside the period every coefficient — k, bias, tanh(k·b) and the peak normaliser in the
    // oversampled loop; drive-compensation, mix and trim at base rate — is interpolated per sample between the
    // two designed sets. So the output moves continuously, and the transcendental cost of a glide is one design
    // per 64 samples, not one per sample.
    //
    // WHY, measured on a -12 dBFS 227 Hz sine through the chain (K = 128), max|Δ²y| where the change reaches
    // the output against -67.2 dBFS for the steady saturated tone: driveDb 3 -> 9 was -18.3 dBFS as a step and
    // 3 -> 3.5, a small knob move, -44.6. The length was chosen against those two and the other continuous
    // parameters — see the item's note in CHANGELOG.
    //
    // CLOCKED BY AUDIO SAMPLES, NEVER BY CALLS (law 8a): the ramp steps on the grid, and a coefficient inside a
    // period is a function of the sample's index in it, so a stream cut any other way — whole, per sample, around
    // maxBlock — renders the same bits for the same event timeline. A clock-only call (nch == 0) spends the same
    // audio time on the glide. It changes no latency and allocates nothing.
    //
    // SNAPPED, NOT RAMPED, before the first sample of a stream: every write between prepare()/reset() and the
    // first accepted call with samples applies at once — so a stage whose parameters were set before its first
    // sample renders exactly what it rendered before the glide existed, which is the whole offline contract of
    // the mastering chain above it. A write that changes nothing continuous costs a comparison (it used to cost
    // a full re-design), so a caller that re-sends its parameters every block does not pay for the glide.
    // `shape` is a topology switch and SNAPS every parameter with it; `dcBlockHz` is a filter coefficient and
    // lands at once, as both always did.
    static constexpr double kGlideMs = 30.0;

    // THE GLIDE LENGTH IS THE CALLER'S TO CHANGE, and 0 turns it off: every write then lands at once, which is
    // the stage exactly as it was before the glide existed, bit for bit (pinned against the frozen pre-change
    // engine in the suite). A product with its own parameter smoothing, or one that wants a hard step, says 0.
    // Takes effect for the next write; a glide in progress when it is set to 0 lands at once. Non-finite or
    // negative is 0. RT-safe; call it from the thread that calls setParams().
    void setGlideMs (double ms) noexcept
    {
        glideMs_ = (std::isfinite (ms) && ms > 0.0) ? ms : 0.0;
        glideTicks_ = ticksFor (fs_, glideMs_);
        if (glideTicks_ == 0 && glideActive()) { applyParams(); snapGlide(); }
    }
    double glideMs() const noexcept { return glideMs_; }

    void setParams (const Params& p) noexcept
    {
        const bool shapeMoved = p.shape != params_.shape;
        params_ = p;
        if (fresh_ || shapeMoved || glideTicks_ == 0) { applyParams(); snapGlide(); return; }
        applyDc();                                                       // a coefficient: lands at once
        float t[kNumP];
        targetOf (params_, t);
        if (std::equal (t, t + kNumP, pt_)) return;                      // nothing continuous moved
        std::copy (t, t + kNumP, pt_);
        pending_ = true;                                                 // the ramp starts at the next boundary
    }

    // How many grid periods a glide takes at the prepared rate: glideMs() / 64 samples, rounded, at least 1 —
    // or 0 when the glide is off (or before prepare()).
    int  glideTicks() const noexcept { return glideTicks_; }
    // True while a written parameter has not landed: a target waiting for the grid, or a ramp in progress.
    bool isGliding() const noexcept { return glideActive(); }

    // In place, planar. RT-safe. n may exceed maxBlock — chunked internally, fully processed.
    // Law 11 (DSP-ARCHITECTURE.md §2). The width used to be clamped and the surplus left BIT-IDENTICAL to
    // its input — measured at drive +24 dB, prepared 2 / called 4: max |out-in| was 0.94633 on a prepared
    // plane and exactly 0 on a surplus one, i.e. no saturation at all. Refused whole now.
    [[nodiscard]] bool process (float* const* io, int numChannels, int n) noexcept
    {
        if (numChannels < 0 || n < 0) return false;
        if (! prepared_) return false;                                   // unprepared / failed-prepare → no OOB
        if (numChannels > channels_) return false;                       // width is a LIMIT — law 11(b)
        if (n == 0) return true;                                         // no samples: no time, no edge
        fresh_ = false;                                                  // the stream has started: writes glide now
        const int nc = numChannels;
        dropStoppedCells (nc);                                           // law 11(d): the edge is clocked by n
        if (nc == 0) { advanceClock (n); return true; }                  // a pause is audio time for the glide too
        // Chunk to maxBlock so a caller passing n > maxBlock is FULLY processed instead of silently
        // dropped. State carries across chunks via the members → bit-identical to one big call.
        // maxBlock_ ≥ 1 whenever prepared_ (prepare() rejects less), so the loop always advances.
        float* sub[core::kMaxChannels];
        for (int off = 0; off < n; )
        {
            const int m = std::min (n - off, maxBlock_);
            for (int c = 0; c < nc; ++c) sub[c] = io[c] + off;
            runChunk (sub, nc, m);
            off += m;                                                    // `off += maxBlock_` could step past INT_MAX
        }
        return true;
    }

private:
    //==========================================================================================================
    // THE GLIDE MACHINERY — see kGlideMs for the contract.
    enum { kDrive, kBias, kAutoComp, kMix, kOutDb, kNumP };

    // The continuous parameters as the stage uses them: non-finite -> the struct default (the house rule),
    // autoComp and mix clamped to [0, 1]. The ramp runs in THIS space, so it never spends time past a clamp.
    static void targetOf (const Params& p, float* t) noexcept
    {
        t[kDrive]    = finite (p.driveDb,  3.0f);
        t[kBias]     = finite (p.bias,     0.0f);
        t[kAutoComp] = std::clamp (finite (p.autoComp, 0.5f), 0.0f, 1.0f);
        t[kMix]      = std::clamp (finite (p.mix,      1.0f), 0.0f, 1.0f);
        t[kOutDb]    = finite (p.outputDb, 0.0f);
    }

    // Everything the audio loops read, at one point of the parameter space.
    struct Consts { WaveShaper::Coeffs sh {}; float comp = 1.0f, mix = 1.0f, out = 1.0f; };

    // THE DESIGN — the one arithmetic both the settled stage and every glide point use, so a glide that lands
    // hands over to the settled path on the very floats it arrived with. A WaveShaper's coefficients are a pure
    // function of (shape, bias, drive): the three setters below leave it where applyParams() leaves `shaper_`.
    Consts design (const float* v, WaveShaper& w) const noexcept
    {
        w.setShape (params_.shape);
        w.setBias  (v[kBias]);
        w.setDrive ((float) (core::dbToGain (v[kDrive]) - 1.0));    // driveDb 0 → k≈0 (linear)
        Consts c;
        c.sh   = w.coeffs();
        c.comp = (float) std::pow ((double) std::max (1.0e-6f, w.slopeAtZero()), (double) -v[kAutoComp]);
        c.mix  = v[kMix];
        c.out  = (float) core::dbToGain (v[kOutDb]);
        return c;
    }

    Consts settled() const noexcept { Consts c; c.sh = shaper_.coeffs(); c.comp = comp_; c.mix = mix_; c.out = outGain_; return c; }

    bool glideActive() const noexcept { return pending_ || ticksLeft_ > 0 || interp_; }

    static int ticksFor (double fs, double ms) noexcept
    {
        if (! (ms > 0.0) || ! (fs > 0.0)) return 0;
        return std::max (1, (int) std::lround (ms * 1.0e-3 * fs / (double) core::StateGrid::kPeriod));
    }

    // The stream's parameters ARE the written ones: no target waits, no ramp runs, the vector sits on it.
    void snapGlide() noexcept
    {
        targetOf (params_, pt_);
        for (int i = 0; i < kNumP; ++i) pv_[i] = (double) pt_[i];
        pending_ = false; ticksLeft_ = 0; interp_ = false;
    }

    // ONE GRID BOUNDARY of a glide. The period that starts here interpolates from where the last one ended
    // (`ce_`, or the settled design when no ramp was running) to the design at the ramp's next point. A target
    // written since the last boundary restarts the ramp from the point it stands on, so a retarget never jumps.
    // When the ramp has landed, the boundary after its last period hands over to the settled path — on the same
    // floats, because `design()` is applyParams()'s own arithmetic.
    void tick() noexcept
    {
        const Consts start = interp_ ? ce_ : settled();
        if (pending_)
        {
            pending_   = false;
            ticksLeft_ = glideTicks_;
            for (int i = 0; i < kNumP; ++i) pd_[i] = ((double) pt_[i] - pv_[i]) / (double) glideTicks_;
        }
        if (ticksLeft_ > 0)
        {
            --ticksLeft_;
            // IN DOUBLE, and landing on the target itself: the ramp ACCUMULATES, and a float accumulator drifts
            // by ~ticks·ulp/2 — at a few MHz enough to overshoot a small target before the last step lands it.
            float v[kNumP];
            for (int i = 0; i < kNumP; ++i)
            {
                pv_[i] = ticksLeft_ > 0 ? pv_[i] + pd_[i] : (double) pt_[i];
                v[i]   = (float) pv_[i];
            }
            WaveShaper w;
            cs_ = start;
            ce_ = design (v, w);
            // Per-sample increments over one period: os samples for the curve, base samples for the rest.
            const float nOs = (float) (core::StateGrid::kPeriod * os_), nB = (float) core::StateGrid::kPeriod;
            dK_  = (ce_.sh.drive    - cs_.sh.drive)    / nOs;
            dB_  = (ce_.sh.bias     - cs_.sh.bias)     / nOs;
            dBt_ = (ce_.sh.biasTanh - cs_.sh.biasTanh) / nOs;
            dN_  = (ce_.sh.norm     - cs_.sh.norm)     / nOs;
            dC_  = (ce_.comp - cs_.comp) / nB;
            dM_  = (ce_.mix  - cs_.mix)  / nB;
            dO_  = (ce_.out  - cs_.out)  / nB;
            interp_ = true;
            return;
        }
        interp_ = false;
        applyParams();                                                   // the target, on ce_'s own floats
    }

    // A whole ≤ maxBlock slice. Settled, it is one pass exactly as before the glide existed; gliding, it is cut at
    // the grid so each piece lies inside one period and knows its place in it.
    void runChunk (float* const* io, int nc, int m) noexcept
    {
        if (! glideActive()) { processChunk (io, nc, m); grid_.skip (m); return; }
        float* sub[core::kMaxChannels];
        for (int o = 0; o < m; )
        {
            if (grid_.phase() == 0 && glideActive()) tick();
            const int seg = grid_.segment (m - o);
            for (int c = 0; c < nc; ++c) sub[c] = io[c] + o;
            if (interp_) processChunkGlide (sub, nc, seg, grid_.phase());
            else         processChunk (sub, nc, seg);
            grid_.advance (seg);
            o += seg;
        }
    }

    // A clock-only call: the same boundaries a call with audio would have met.
    void advanceClock (int n) noexcept
    {
        if (! glideActive()) { grid_.skip (n); return; }
        for (int o = 0; o < n; )
        {
            if (grid_.phase() == 0 && glideActive()) tick();
            const int seg = grid_.segment (n - o);
            grid_.advance (seg);
            o += seg;
        }
    }

    // The interpolated curve over one os-rate run starting at os index `r0` of the period: every coefficient is
    // `start + d * index`, the product and the sum in SEPARATE statements so `-ffp-contract=on` fuses neither
    // (law 10) — the value is a function of the index alone, never of an accumulator a cut could reset.
    template <WaveShaper::Shape S>
    void shapeGlide (float* b, int osN, int r0, int c) noexcept
    {
        const bool dc = dcEnabled_;
        float x1 = dc ? dcX1_[(std::size_t) c] : 0.0f, y1 = dc ? dcY1_[(std::size_t) c] : 0.0f;
        for (int i = 0; i < osN; ++i)
        {
            const float r = (float) (r0 + i);
            const float ek = dK_ * r, eb = dB_ * r, et = dBt_ * r, en = dN_ * r;
            WaveShaper::Coeffs k;
            k.drive = cs_.sh.drive + ek; k.bias = cs_.sh.bias + eb; k.biasTanh = cs_.sh.biasTanh + et; k.norm = cs_.sh.norm + en;
            const float w = WaveShaper::shapeAt<S> (k, b[i]);
            if (dc)
            {
                const float d = w - x1 + dcR_ * y1;                      // the settled loop's DC blocker, verbatim
                x1 = w;
                y1 = (std::fabs (d) < 1e-30f) ? 0.0f : d;
                b[i] = y1;
            }
            else b[i] = w;
        }
        if (dc)
        {
            dcX1_[(std::size_t) c] = std::isfinite (x1) ? x1 : 0.0f;
            dcY1_[(std::size_t) c] = std::isfinite (y1) ? y1 : 0.0f;
        }
    }

    // Clear the sample memory of every cell that ran on the previous accepted call and does not run on this
    // one. A channel that leaves and RETURNS is the case: its oversampler FIR, its DC blocker and its dry
    // delay are frozen, not decayed, and it replays them into a stream that has moved on — measured 0.9337
    // out of DIGITAL SILENCE, 29 samples after a return. The DC blocker has a second gate of its own,
    // dcEnabled_, which the Asym→symmetric→Asym path opens and closes at a constant channel count and which
    // freezes x1/y1 exactly the same way. Per channel, never wholesale: a channel that never left keeps its
    // history bit-exact, which is why the oversampler grew resetChannel(). A call that carries no samples
    // (or none this stage accepts) ran nothing, so it stops nothing and never reaches here.
    void dropStoppedCells (int nc) noexcept
    {
        const int nowDc = dcEnabled_ ? nc : 0;
        for (int c = nc; c < ranNc_; ++c)
        {
            if (os_ > 1) ovs_.resetChannel (c);
            dryDelay_[(std::size_t) c].reset();
        }
        for (int c = nowDc; c < ranDcNc_; ++c)
        {
            dcX1_[(std::size_t) c] = 0.0f;
            dcY1_[(std::size_t) c] = 0.0f;
        }
        ranNc_ = nc; ranDcNc_ = nowDc;
    }

    // One ≤ maxBlock slice; nc already clamped by process(). Everything stateful streams across calls.
    void processChunk (float* const* io, int nc, int n) noexcept
    {
        if (n <= 0) return;                                              // belt-and-suspenders (process() guarantees n ≥ 1)
        const int osN = n * os_;
        for (int c = 0; c < nc; ++c)
        {
            osPtrs_[(std::size_t) c]  = &osBuf_[(std::size_t) c * (std::size_t) (maxBlock_ * os_)];
            wetPtrs_[(std::size_t) c] = &wetBuf_[(std::size_t) c * (std::size_t) maxBlock_];
        }

        // 0) sanitize AT THE GATE, in place — so the wet path AND the dry tap below both read the
        //    sanitized signal: a NaN/Inf sample would lodge in the oversampler FIR / DC-blocker /
        //    dry-delay state (NaN never decays out of a recursion); the clamp also catches an extreme
        //    finite input that would overflow downstream. Finite in-range samples pass bit-identically
        //    (std::clamp returns the value untouched inside the bounds).
        for (int c = 0; c < nc; ++c)
            for (int i = 0; i < n; ++i)
            {
                const float v = io[c][i];
                io[c][i] = std::clamp (std::isfinite (v) ? v : 0.0f, -1.0e6f, 1.0e6f);
            }

        // 1) to the oversampled domain (io is only READ here, so it still holds the dry signal below)
        if (os_ > 1) ovs_.upsample (io, nc, n, osPtrs_.data());
        else for (int c = 0; c < nc; ++c) std::copy (io[c], io[c] + n, osPtrs_[(std::size_t) c]);

        // 2) waveshape per channel in the oversampled domain; DC-block only the asymmetric curve (the
        //    even-harmonic offset). Symmetric curves are zero-mean → no blocker → no needless phase shift.
        for (int c = 0; c < nc; ++c)
        {
            float* b = osPtrs_[(std::size_t) c];
            if (dcEnabled_)
            {
                float x1 = dcX1_[(std::size_t) c], y1 = dcY1_[(std::size_t) c];
                for (int i = 0; i < osN; ++i)
                {
                    const float w  = shaper_.processSample (b[i]);
                    const float dc = w - x1 + dcR_ * y1;     // one-pole DC blocker  H(z)=(1-z⁻¹)/(1-R·z⁻¹)
                    x1 = w;
                    y1 = (std::fabs (dc) < 1e-30f) ? 0.0f : dc;   // law 8 — see the note below the loop
                    b[i] = y1;
                }
                // LAW 8, and the rationale that used to stand here is retracted rather than edited. It
                // read: "no denormal threshold: zapping a tiny finite tail would break the bit-identical
                // NULL for quiet legitimate signals." Right about the risk, wrong about the cost of doing
                // nothing, and — this is the part worth recording — the NULL it appealed to never tested
                // the claim either way: scenario 4's 1e-20 block is absorbed to exact zero by the Asym
                // bias (-0.3f + 1e-20f == -0.3f), so the fixture never puts a tiny FINITE value in this
                // state at all. Four threshold variants were built against it and all 194 checks pass.
                //
                // What doing nothing actually cost. On silence the Asym curve gives w == 0 exactly (it is
                // normalised so y(0) == 0), so x1 goes to zero and y1 becomes a pure R^n decay — which
                // STICKS: with R = 0.99967 (10 Hz at 4x 48 kHz) every subnormal k*u with k <= 0.5/(1-R)
                // ~= 1528 maps to itself under round-to-nearest, so y1 froze at 2.14e-42 after 1.5 s and
                // never underflowed. The stuck value then fed the 128-tap downsampling FIR, every operand
                // subnormal, and LEFT the stage: -3.40249e-41 on every output sample, bit-identical on
                // arm64 and x86, i.e. contagion into whatever comes next. Measured on x86-64 (i9-13900H,
                // gcc 14.2 -O2, 2 s mono, best of 5): music 18.97 ms, stalled silence 452.72 ms — 23.9x,
                // 22.6 %RT per channel. Forcing FTZ/DAZ on the same binary: 18.33 ms, no penalty at all,
                // which is what proves it is purely the denormal path. With the flush: 18.25 ms.
                //
                // The flush is PER SAMPLE (in the loop above), not here, because a threshold applied at
                // the end of a call would put a numerical event wherever the CALLER cut the stream — the
                // chunk-invariance argument from docs/LAW8-KWEIGHTING.md, and this module asserts
                // bit-identical chunking explicitly. 1e-30f rather than core::flushDenormal's 1e-15f is a
                // margin choice, not a NULL-forced one: it is 8 decades above the float subnormal floor
                // and ~600 dB below anything audible, so it cannot be reached by a signal while still
                // being unreachable from below by the stall. This per-call line keeps its ORIGINAL job,
                // the poison guard: a NaN in recursive state never decays out on its own.
                dcX1_[(std::size_t) c] = std::isfinite (x1) ? x1 : 0.0f;
                dcY1_[(std::size_t) c] = std::isfinite (y1) ? y1 : 0.0f;
            }
            else
            {
                for (int i = 0; i < osN; ++i) b[i] = shaper_.processSample (b[i]);
            }
        }

        // 3) back to base rate (wet)
        if (os_ > 1) ovs_.downsample (osPtrs_.data(), nc, n, wetPtrs_.data());
        else for (int c = 0; c < nc; ++c) std::copy (osPtrs_[(std::size_t) c], osPtrs_[(std::size_t) c] + n, wetPtrs_[(std::size_t) c]);

        // 4) drive-compensate the wet, then a linear (peak-safe) dry/wet, then output trim. The dry runs
        //    through a DelayLine matching the oversampler round-trip — an undelayed dry combs the wet at
        //    mix < 1 (63 samples at the default ≈ −5.1 dB of comb ripple at 1 kHz / 48 k; it was 31 and
        //    −7 dB before the taps default rose, and the notch spacing halved with it).
        for (int c = 0; c < nc; ++c)
        {
            core::DelayLine& dl = dryDelay_[(std::size_t) c];
            for (int i = 0; i < n; ++i)
            {
                const float dry = dl.process (io[c][i]);
                const float wet = comp_ * wetPtrs_[(std::size_t) c][i];
                io[c][i] = outGain_ * ((1.0f - mix_) * dry + mix_ * wet);
            }
        }
    }

    // processChunk() for a slice INSIDE ONE GLIDING GRID PERIOD, starting at base index `ph0` of it: the same
    // gate, the same oversampler, the same DC blocker and the same dry path, with the curve's coefficients and
    // comp/mix/trim interpolated per sample (see shapeGlide). A separate function on purpose: the settled slice
    // above is left exactly as it was, code and all, so a stage that is not gliding pays nothing for this one.
    void processChunkGlide (float* const* io, int nc, int n, int ph0) noexcept
    {
        if (n <= 0) return;
        const int osN = n * os_;
        for (int c = 0; c < nc; ++c)
        {
            osPtrs_[(std::size_t) c]  = &osBuf_[(std::size_t) c * (std::size_t) (maxBlock_ * os_)];
            wetPtrs_[(std::size_t) c] = &wetBuf_[(std::size_t) c * (std::size_t) maxBlock_];
        }
        for (int c = 0; c < nc; ++c)
            for (int i = 0; i < n; ++i)
            {
                const float v = io[c][i];
                io[c][i] = std::clamp (std::isfinite (v) ? v : 0.0f, -1.0e6f, 1.0e6f);
            }
        if (os_ > 1) ovs_.upsample (io, nc, n, osPtrs_.data());
        else for (int c = 0; c < nc; ++c) std::copy (io[c], io[c] + n, osPtrs_[(std::size_t) c]);
        for (int c = 0; c < nc; ++c)
        {
            float* b = osPtrs_[(std::size_t) c];
            switch (params_.shape)
            {
                case WaveShaper::Shape::Tanh:  shapeGlide<WaveShaper::Shape::Tanh>  (b, osN, ph0 * os_, c); break;
                case WaveShaper::Shape::Atan:  shapeGlide<WaveShaper::Shape::Atan>  (b, osN, ph0 * os_, c); break;
                case WaveShaper::Shape::Cubic: shapeGlide<WaveShaper::Shape::Cubic> (b, osN, ph0 * os_, c); break;
                case WaveShaper::Shape::Asym:  shapeGlide<WaveShaper::Shape::Asym>  (b, osN, ph0 * os_, c); break;
            }
        }
        if (os_ > 1) ovs_.downsample (osPtrs_.data(), nc, n, wetPtrs_.data());
        else for (int c = 0; c < nc; ++c) std::copy (osPtrs_[(std::size_t) c], osPtrs_[(std::size_t) c] + n, wetPtrs_[(std::size_t) c]);
        for (int c = 0; c < nc; ++c)
        {
            core::DelayLine& dl = dryDelay_[(std::size_t) c];
            for (int i = 0; i < n; ++i)
            {
                const float j = (float) (ph0 + i);
                const float ec = dC_ * j, em = dM_ * j, eo = dO_ * j;         // separate statements: law 10
                const float comp = cs_.comp + ec, mix = cs_.mix + em, out = cs_.out + eo;
                const float dry = dl.process (io[c][i]);
                const float wet = comp * wetPtrs_[(std::size_t) c][i];
                io[c][i] = out * ((1.0f - mix) * dry + mix * wet);
            }
        }
    }

    static float finite (float v, float fallback) noexcept { return std::isfinite (v) ? v : fallback; }

    // The SETTLED design, at the written parameters. Non-finite params fall back to the struct defaults (house
    // rule) — std::clamp passes NaN through. `design()` is the arithmetic, shared with every glide point.
    void applyParams() noexcept
    {
        float t[kNumP];
        targetOf (params_, t);
        const Consts c = design (t, shaper_);
        comp_    = c.comp;
        mix_     = c.mix;
        outGain_ = c.out;
        applyDc();
    }

    void applyDc() noexcept
    {
        const float dcHz = finite (params_.dcBlockHz, 10.0f);
        const double fsOs = fs_ * (double) os_;
        const double fc   = std::clamp ((double) dcHz, 0.0, 0.49 * fsOs);
        dcR_ = (fc <= 0.0) ? 0.0f : (float) std::exp (-2.0 * core::kPi * fc / fsOs);
        dcEnabled_ = (dcHz > 0.0f) && (params_.shape == WaveShaper::Shape::Asym);
    }

    Params      params_ {};
    WaveShaper  shaper_ {};
    oversampling::Oversampler ovs_;

    double fs_ = 0.0;
    int    maxBlock_ = 0, channels_ = 0, os_ = 1;
    float  comp_ = 1.0f, dcR_ = 0.0f, mix_ = 1.0f, outGain_ = 1.0f;
    bool   dcEnabled_ = false;
    bool   prepared_  = false;                             // true only after a fully-successful prepare()

    // The glide (see kGlideMs). `pv_` is the parameter vector at the last grid boundary, `pt_` the target, `pd_`
    // the ramp's per-boundary step; `cs_`/`ce_` the designs the current period interpolates between, and the
    // d*_ members their per-sample increments (os rate for the curve, base rate for comp/mix/trim).
    core::StateGrid grid_;
    double glideMs_ = kGlideMs;
    int    glideTicks_ = 0, ticksLeft_ = 0;
    bool   fresh_ = true, pending_ = false, interp_ = false;
    double pv_[kNumP] {}, pd_[kNumP] {};
    float  pt_[kNumP] {};
    Consts cs_ {}, ce_ {};
    float  dK_ = 0.0f, dB_ = 0.0f, dBt_ = 0.0f, dN_ = 0.0f, dC_ = 0.0f, dM_ = 0.0f, dO_ = 0.0f;
    int    ranNc_ = 0, ranDcNc_ = 0;                       // what advanced state on the previous accepted call

    std::vector<float>  osBuf_, wetBuf_;
    std::vector<float*> osPtrs_, wetPtrs_;
    std::vector<float>  dcX1_, dcY1_;
    std::vector<core::DelayLine> dryDelay_;    // per channel, os round-trip — keeps the dry/wet mix comb-free
};

} // namespace felitronics::saturation

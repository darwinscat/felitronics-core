// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/core/Config.h>
#include <felitronics/core/Math.h>
#include <felitronics/core/DelayLine.h>
#include <felitronics/core/FlushToZero.h>
#include <felitronics/oversampling/PolyphaseOversampler.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

namespace felitronics::limiter
{

namespace detail
{
    // O(1)-amortized sliding-window maximum (monotonic-decreasing deque). Front = max over the last W.
    // prepare() fixes the CAPACITY; setWindow() picks the effective window ≤ capacity. The lookahead is
    // prepare-time now, so the window is set once per prepare rather than retuned. Expiry runs BEFORE the
    // insert — insert-first can overwrite the head (the current max) once a strictly-decreasing run fills
    // the deque, under-reading the max exactly when a decaying peak is still inside the lookahead. Indices
    // are int64 (a `long` is 32-bit on Windows → wraps after ~3 h of oversampled pushes).
    class SlidingMax
    {
    public:
        void prepare (int maxWindow)
        {
            cap = maxWindow < 1 ? 1 : maxWindow;
            v.assign ((std::size_t) cap, 0.0f); ix.assign ((std::size_t) cap, 0);
            W = cap;
            reset();
        }

        void reset() noexcept { head = tail = count = 0; n = 0; }

        void setWindow (int w) noexcept { W = w < 1 ? 1 : (w > cap ? cap : w); }   // stale entries expire on the next pushes

        float push (float x) noexcept
        {
            while (count > 0 && ix[(std::size_t) head] <= n - (std::int64_t) W) { head = (head + 1) % cap; --count; }
            while (count > 0) { const int b = (tail - 1 + cap) % cap; if (v[(std::size_t) b] <= x) { tail = b; --count; } else break; }
            v[(std::size_t) tail] = x; ix[(std::size_t) tail] = n; tail = (tail + 1) % cap; ++count;
            ++n;
            return v[(std::size_t) head];
        }

    private:
        int cap = 1, W = 1; std::vector<float> v; std::vector<std::int64_t> ix; int head = 0, tail = 0, count = 0; std::int64_t n = 0;
    };
}

//==============================================================================
// TOPOLOGY, fixed for the life of a prepared stream. These three decide the FIR design, the buffer
// sizes and the reported latency, so none of them can be a per-block parameter: changing the factor
// has no state-preserving mapping between a live 2×/4×/8× history, and changing the lookahead moves
// latencySamples(), which in any host is a resynchronisation event rather than automation. They live
// here, and not in TruePeakLimiterParams, so that a caller cannot express the change at all — an
// earlier revision accepted both silently and applied neither, and raising `lookaheadMs` mid-stream
// re-emitted already-delivered audio at +28.8 dB over the ceiling.
struct TruePeakLimiterConfig
{
    double lookaheadMs      = 1.0;
    int    oversampleFactor = 4;     // ≥ 2; a requested 1 becomes 2 — there is no 1× path
    int    tapsPerPhase     = 32;    // ≥ 4
};

// PER-BLOCK parameters — safe to change at any time, from the audio thread, mid-stream.
struct TruePeakLimiterParams
{
    double ceilingDbTp = -1.0;       // output true-peak ceiling (dBTP) — see the derate note below
    double releaseMs   = 100.0;
};

//==============================================================================
// felitronics::limiter::TruePeakLimiter — a lookahead limiter that BOUNDS EVERY SAMPLE ON ITS OWN
// F×fs GRID. It oversamples, limits at the oversampled rate (so inter-sample peaks are real samples),
// then downsamples — Option B, the only structure that can bound inter-sample peaks at all, since a
// fast-moving gain applied at baseband would itself create them.
//
// Gain law: a sliding-window MAXIMUM of the channel-linked oversampled peak over the lookahead window
// → required gain = ceiling / slidingMax (so EVERY sample in the window, including the one being output
// `lookahead` samples behind, is ≤ ceiling) → instant attack + rate-limited (release) recovery. The
// emitted sample is necessarily inside its own detector window and the release branch can only keep
// MORE reduction than required, so on that grid the bound is algebra, not a heuristic — which is the
// whole difference between this and a lookahead gain ramp.
//
// ====================================================================================
// WHAT IT DOES NOT PROMISE — read before setting a ceiling from a platform's number.
// ====================================================================================
// The grid bound is not a bound on the reconstructed output a meter or a listener sees. Two terms open
// a gap between them, and they behave nothing alike (measured in ctest, felitronics_limiter_ceiling_tests):
//
//   * GRID GEOMETRY — closed form, and oversampling shrinks it. A crest can hide between detector
//     samples; for a tone at fs·p/q the excess is exactly −20·log10(cos(π/M')), where M' is the number
//     of distinct magnitude phases it visits on the F×fs grid. The worst tone the round trip passes FLAT
//     is fs/3: +1.250 / +0.302 / +0.076 dB at 2× / 4× / 8×. (Geometry alone is worse higher up — 2fs/5
//     gives 0.436 dB at 4× — but the prototype's own droop removes that content before delivery, so it
//     cannot reach the output. Widen the pass band and this term has to be re-derived.)
//   * GAIN MODULATION — and neither oversampling nor a slower release removes it. The attack is
//     instantaneous, so the limited product is not band-limited and the downsampler overshoots
//     re-band-limiting it. It SATURATES near 0.92 dB (4×) / 0.87 (8×) above ~5 dB of reduction, and
//     does so at a 100 ms release as readily as at 1 ms. What it answers to is spectral tilt: at 4× and
//     a 50 ms release, flat noise costs 0.92 dB, roughly music-like material (one-pole, a = 0.9) costs
//     0.10 dB, and dark material 0.04 dB.
//
// Worst measured total over the witness matrix INSIDE the domain below: +1.25 / +0.99 / +0.93 dB at
// 2× / 4× / 8× (the fs/3 tone at 2×, dense material at a 1 ms release at 4× and 8×).
// NEITHER term grows with crest factor or with makeup gain (verified at +64, +76 and +88 dB) — which
// is what makes the gap characterisable, where a gain-ramp limiter misses by more the harder the
// material gets.
//
// READ THE STATUS OF THOSE NUMBERS CORRECTLY. They are a CHARACTERISATION over the witness matrix in
// felitronics_limiter_ceiling_tests at a lookahead of ≥ 1 ms and a release of ≥ 1 ms — not a theorem,
// and the floors below are NOT inside that domain. Two measured examples of leaving it: alternating
// ±500000 (i.e. +114 dBFS, inside the gate) at the release floor delivers +2.67 dB over, where the same
// shape at a 1 ms release delivers +0.44 and at 50 ms +0.004; and with BOTH parameters on their floors
// the dense witness reaches +1.59 (4×) / +1.56 (8×), past the figures above. Bringing even that corner
// inside would need floors at 24 baseband samples (0.5 ms), which is a musical setting and would change
// the sound of a legitimate one — so the corner is documented instead of clamped away. The only thing
// proven for every input is the on-grid bound.
//
// SO: a product ceiling needs a DERATE, or a loop closed on measured true peak. Setting ceilingDbTp
// to C does not deliver C dBTP. Over the characterised domain, budget ~1.2 dB for a bright, dense,
// hard-limited master and ~0.3 dB for ordinary material, at 4× or 8× alike — and if you need a
// GUARANTEE rather than a budget, measure the delivered peak and close the loop on it.
//
// ALWAYS IN THE PATH, even when nothing is being limited: the 0.90 × Nyquist prototype is a round trip,
// so the top of the band is attenuated — 0.00 dB at 0.36 fs, −0.40 at 0.40, −4.0 at 0.44, −6.0 at 0.45,
// identically at every factor (the transition width is constant in base-rate units). That is a mastering
// decision, not a rounding error, and it is pinned in the suite so a wider pass band has to be chosen
// deliberately.
//
// WHAT IS CLAMPED, all of it visible rather than silent (oversampleFactor(), lookaheadSamples(),
// effectiveReleaseMs(), effectiveCeilingDbTp()): the oversampling factor into [2, 16]; the lookahead
// into [2 baseband samples, 20 ms]; the release time constant to at least 8 baseband samples and to
// strictly less than an infinite hold; the ceiling into [-200, +60] dBTP. The two floors are measured —
// the smallest values that keep the witness matrix inside the figures above — and both sit far below any
// musical setting (42 µs and 167 µs at 48 kHz). They prevent the DEGENERACY (a lookahead of zero is an
// OS-rate clipper; a release of zero lets the gain jump every oversampled sample); they do not turn the
// characterisation below into a bound.
//
// CHANNEL COUNT is part of the topology too. Only the channels passed to process() advance their own
// history, while the detector and the gain are shared, so CHANGING the count clears the state (a
// discontinuity, which is what such a change already is) — without that, a channel that sat out a few
// blocks came back re-emitting audio at +19.76 dB over the ceiling. Passing MORE channels than
// prepared limits the first maxChannels of them and leaves the rest untouched: prepare for the most you
// will ever pass.
//
// RT-safe: prepare() allocates; process() does no alloc/lock/throw and accepts any block length. One
// linked gain for all channels (no image shift). Non-finite and absurd input samples are sanitised at
// the gate — bit-transparently for finite audio within ±1e6.
class TruePeakLimiter
{
public:
    // Topology is chosen HERE and nowhere else. Returns false and leaves the limiter unprepared on an
    // unusable stream configuration; a false return is the only way to learn that, so check it.
    [[nodiscard]] bool prepare (double sampleRate, int maxBlock, int maxChannels, const TruePeakLimiterConfig& config = {})
    {
        prepared_ = false;                                     // any early return below leaves it unprepared
        // Spelled positively so NaN fails: `sampleRate <= 0.0` is FALSE for NaN, which let a NaN rate
        // through to std::lround(NaN) — undefined behaviour, and a signed overflow right after it.
        if (! (sampleRate > 0.0) || ! std::isfinite (sampleRate) || sampleRate > kMaxSampleRate) return false;
        if (maxBlock < 1) return false;
        // Both ends. Only the lower bound was checked, so tapsPerPhase = INT_MAX reached N = L*tpp in
        // the oversampler and overflowed a signed int before allocating (UBSan-confirmed, then a
        // length_error — a terminate under the wasm tier's -fno-exceptions).
        if (config.tapsPerPhase < 4 || config.tapsPerPhase > kMaxTapsPerPhase) return false;
        if (! std::isfinite (config.lookaheadMs) || config.lookaheadMs < 0.0) return false;

        fs    = sampleRate;
        // Refused rather than clamped: a silently reduced channel count would leave the surplus
        // channels passing through this module UNLIMITED, and a silently reduced factor would run a
        // topology the caller did not ask for. Both are the class of defect this task closes.
        if (maxChannels < 1 || maxChannels > core::kMaxChannels) return false;
        if (config.oversampleFactor > kMaxFactor) return false;
        maxCh = maxChannels;
        // CLAMPED, not refused. Since process() chunks, maxBlock is a scratch-buffer size rather than a
        // limit on what a caller may pass — and refusing a big one would fail prepare(), after which
        // process() returns the buffer UNTOUCHED, i.e. exactly the unlimited-passthrough defect this
        // task exists to close. An offline caller sizing maxBlock to a whole file is a normal thing to do.
        maxBlock_ = std::min (maxBlock, kMaxBlock);
        tpp   = config.tapsPerPhase;
        F     = config.oversampleFactor < 2 ? 2 : config.oversampleFactor;   // a requested 1 becomes 2; above kMaxFactor was refused
        if (! os.prepare (F, maxCh, tpp)) return false;        // oversampler rejected → stay unprepared

        osBuf.assign ((std::size_t) maxCh, std::vector<float> ((std::size_t) maxBlock * (std::size_t) F, 0.0f));
        osPtrs.assign ((std::size_t) maxCh, nullptr);

        // A rate so low that 20 ms cannot hold the minimum lookahead would make the clamp below
        // std::clamp(x, 2, 1) — lo > hi is undefined behaviour, and the reported lookahead would then
        // disagree with the delay the DelayLine actually took. Refuse instead.
        const int maxLookBb  = (int) std::ceil (kMaxLookaheadMs * 0.001 * fs);
        if (maxLookBb < kMinLookaheadSamples) return false;
        const int maxLookOS  = maxLookBb * F;
        osDelays.assign ((std::size_t) maxCh, core::DelayLine {});
        for (auto& d : osDelays) d.prepare (maxLookOS);
        slide.prepare (maxLookOS + 1);

        // The lookahead floor is not taste. At zero the structure degenerates into a clipper at the
        // oversampled rate: +2.05 dB over the ceiling, measured, and unlike the terms above that one has
        // no closed form to derate against. Two baseband samples is where the witness matrix comes back
        // inside the envelope with the release at a musical value, at both conformance rates and every
        // factor. (One sample is already inside at a 50 ms release; two is the value that also holds the
        // click at the release floor, which is why the pair is stated together.) At any musical lookahead
        // this clamp is invisible: it is 42 µs at 48 kHz.
        // Clamp in the DOUBLE domain BEFORE lround: `std::lround(1e300)` is out of range for a long,
        // and the result of that then clamped upward landed on the MINIMUM lookahead — a caller asking
        // for an absurd value silently got the smallest one instead of the largest.
        const double lookMs = std::clamp (config.lookaheadMs, 0.0, kMaxLookaheadMs);
        lookBaseband = (int) std::lround (lookMs * 0.001 * fs);
        lookBaseband = std::clamp (lookBaseband, kMinLookaheadSamples, maxLookBb);
        const int lookOS = lookBaseband * F;
        for (auto& d : osDelays) d.setDelay (lookOS);
        slide.setWindow (lookOS + 1);                          // the window must equal the ACTUAL lookahead + the emitted sample

        apply (params);
        reset();
        prepared_ = true;                                      // fully built — process() may now run
        return true;
    }

    void reset() noexcept
    {
        os.reset();                                            // the up/down FIR histories are state too: without
                                                               // this, silence after a reset came out at the ceiling
                                                               // for ~110 samples of the previous stream
        for (auto& b : osBuf) std::fill (b.begin(), b.end(), 0.0f);
        for (auto& d : osDelays) d.reset();
        slide.reset();
        grDb = 0.0f;
        lastNc_ = 0;
    }

    void setParams (const TruePeakLimiterParams& p) noexcept { params = p; apply (p); }

    // All of these read as "nothing prepared" after a failed prepare(), rather than reporting the
    // topology of whatever was prepared before it — a stale latency is worse than an obvious zero.
    int    latencySamples()  const noexcept { return prepared_ ? os.latencySamples() + lookBaseband : 0; }
    double gainReductionDb() const noexcept { return grDb; }
    bool   isPrepared()      const noexcept { return prepared_; }

    // The EFFECTIVE topology, after the clamps above — a caller that asked for something outside the
    // supported range can see what it actually got instead of guessing.
    int    oversampleFactor()  const noexcept { return prepared_ ? F : 0; }
    int    lookaheadSamples()  const noexcept { return prepared_ ? lookBaseband : 0; }
    double effectiveReleaseMs()   const noexcept { return relMsEffective; }
    double effectiveCeilingDbTp() const noexcept { return ceilingDb; }

    // Audio thread, in place (baseband). RT-safe. `numSamples` may exceed the maxBlock passed to
    // prepare() — chunked internally, state carries across chunks, so the result is bit-identical to
    // the caller having made maxBlock-sized calls. It used to return the buffer UNTOUCHED instead,
    // i.e. unlimited and silently, which for a module whose only promise is a ceiling is the worst
    // possible failure.
    void process (float* const* channels, int numChannels, int numSamples) noexcept
    {
        const int nc = numChannels < maxCh ? numChannels : maxCh;
        if (! prepared_ || nc <= 0 || numSamples <= 0) return;

        // A CHANGE OF CHANNEL COUNT IS A TOPOLOGY CHANGE, and it has to clear the state. Only the
        // supplied channels advance their oversampler history and lookahead delay, while the detector
        // deque and the gain are SHARED and advance regardless — so a channel that sits out a few blocks
        // comes back holding audio whose detector entries expired globally while it was away. Measured
        // before this guard: stereo, then 40 samples with a +20 dBFS impulse on the right, then mono
        // silence, then stereo again → the right channel re-emitted at +19.76 dB OVER the ceiling, at
        // every factor. It falsifies the proof above directly: the emitted sample was no longer inside
        // its own detector window. Resetting costs a discontinuity, which is what a channel-count change
        // already is, and it restores the bound from the first sample after it.
        if (lastNc_ != 0 && nc != lastNc_) reset();
        lastNc_ = nc;

        float* sub[core::kMaxChannels] {};
        for (int off = 0; off < numSamples; )
        {
            const int n = std::min (numSamples - off, maxBlock_);
            for (int c = 0; c < nc; ++c) sub[(std::size_t) c] = channels[c] + off;
            processChunk (sub, nc, n);
            off += n;                                          // `off += maxBlock_` could step past INT_MAX
        }
    }

private:
    static constexpr double kMaxLookaheadMs      = 20.0;
    static constexpr double kMaxSampleRate       = 3.0e6;   // far above any audio rate; keeps the derived
                                                            // sizes below INT_MAX with the 20 ms capacity
    // maxBlock is a WORKING block size, not a file length — process() chunks anything larger — so the
    // cap is set by what the oversampled scratch buffer may cost per channel (maxBlock x F floats), not
    // by any use case. At the extremes allowed here that is 64 MB/channel; without a cap a hostile or
    // mistaken prepare() could ask for gigabytes.
    static constexpr int    kMaxBlock            = 1 << 20;
    static constexpr int    kMaxFactor           = 16;
    static constexpr int    kMaxTapsPerPhase     = 1024;    // N = factor*taps must stay well inside int
    static constexpr int    kMinLookaheadSamples = 2;        // measured floor — see prepare()
    static constexpr int    kMinReleaseSamples   = 8;        // measured floor — see apply()
    static constexpr double kMinCeilingDb        = -200.0;   // below the 24-bit floor; -1e308 overflows the
    static constexpr double kMaxCeilingDb        =   60.0;   // float cast and kills the gain permanently

    void processChunk (float* const* channels, int nc, int numSamples) noexcept
    {
        const int osN = numSamples * F;

        // GATE. One bad sample used to be fatal: gainToDb(inf) is inf, so rawRedDb went to -inf, and
        // -inf * relCoef stays -inf — the rest of the stream became digital silence. Sanitising here,
        // BEFORE upsample(), is what makes that unreachable: the oversampler copies raw input into its
        // history, where poison would then survive for tapsPerPhase samples and no downstream guard
        // could remove it. Identical shape to saturation::Saturator's gate, and bit-transparent by
        // construction — a finite sample within ±1e6 takes neither branch.
        for (int c = 0; c < nc; ++c)
            for (int i = 0; i < numSamples; ++i)
            {
                const float v = channels[c][i];
                channels[c][i] = std::clamp (std::isfinite (v) ? v : 0.0f, -1.0e6f, 1.0e6f);
            }

        for (int c = 0; c < nc; ++c) osPtrs[(std::size_t) c] = osBuf[(std::size_t) c].data();
        os.upsample (channels, nc, numSamples, osPtrs.data());

        for (int i = 0; i < osN; ++i)
        {
            float linkedPeak = 0.0f;
            for (int c = 0; c < nc; ++c) { const float a = std::fabs (osBuf[(std::size_t) c][(std::size_t) i]); if (a > linkedPeak) linkedPeak = a; }

            const float  smax    = slide.push (linkedPeak);
            const double smaxDb  = core::gainToDb (smax);
            double rawRedDb = ceilingDb - smaxDb;
            if (rawRedDb > 0.0) rawRedDb = 0.0;

            grDb = std::min ((float) rawRedDb, grDb * relCoef);   // instant attack, exponential release toward 0 dB
            const float gain = (float) core::dbToGain ((double) grDb);

            for (int c = 0; c < nc; ++c)
            {
                const float x = osBuf[(std::size_t) c][(std::size_t) i];
                osBuf[(std::size_t) c][(std::size_t) i] = osDelays[(std::size_t) c].process (x) * gain;
            }
        }

        os.downsample ((const float* const*) osPtrs.data(), nc, numSamples, channels);
        // The flush ITSELF is pinned — removing it leaves the gain state at a denormal instead of exact
        // zero after limiting plus silence, and the suite fails. What is NOT distinguishable is this
        // variant from flushDenormal, and that is correct rather than a gap: given the gate above and the
        // ceiling clamp in apply(), rawRedDb is bounded and grDb cannot go non-finite by any path, so the
        // two behave identically on every reachable state. Kept as the house-standard form for recursive
        // state; the thing that actually closed the poison defect is the gate.
        core::flushPoison (grDb);
    }

    void apply (const TruePeakLimiterParams& p) noexcept
    {
        // Non-finite params fall back to the defaults (house rule) — a NaN ceiling would poison grDb.
        // And a FINITE but absurd one is just as fatal, which is less obvious: at -1e308 the float cast
        // of (ceiling - smaxDb) overflows to -inf, grDb sticks at -inf and the stream goes to digital
        // silence; with the poison flush below that instead becomes a full release to unity, i.e. the
        // material ships UNLIMITED — worse, for a module whose only promise is a ceiling. So the range
        // is clamped to what a true-peak ceiling can mean at all. Nothing a caller would ever set moves.
        ceilingDb = std::clamp (std::isfinite (p.ceilingDbTp) ? p.ceilingDbTp : -1.0, kMinCeilingDb, kMaxCeilingDb);

        // The release floor, like the lookahead floor, is measured rather than chosen: at zero the
        // coefficient is zero and the gain may jump on every oversampled sample. Measured with the floor
        // removed, the click train — the worst witness for this, not the plateau — reaches +2.7 dB over
        // the ceiling at 8×, against +0.63 at this floor. Eight baseband samples is the smallest floor at
        // which each witness lands back inside the envelope with the OTHER parameter at a musical value,
        // at both conformance rates and every factor; with both parameters on their floors at once the
        // corner is outside it, and that is documented in the header rather than clamped away. 167 µs at
        // 48 kHz, i.e. far below any musical release.
        const double relMs = std::isfinite (p.releaseMs) ? p.releaseMs : 100.0;
        const double tMin  = (double) kMinReleaseSamples * (double) F;         // in OVERSAMPLED samples
        const double t     = std::max (relMs * 0.001 * fs * (double) F, tMin);
        // relCoef must stay STRICTLY below 1, or the release stops existing: measured, a release of
        // 175 s at 48 kHz/4x rounds exp(-1/t) to exactly 1.0f and the gain then never recovers at all,
        // while effectiveReleaseMs() still reports the finite value the caller asked for. Backing off to
        // the largest float below 1 keeps recovery monotone without touching any normal setting (the
        // coefficient is left in float on purpose — a double would change the arithmetic and with it
        // every existing sample).
        relCoef = (float) std::exp (-1.0 / t);
        if (! (relCoef < 1.0f)) relCoef = std::nextafterf (1.0f, 0.0f);
        const double tEff = (relCoef > 0.0f) ? -1.0 / std::log ((double) relCoef) : 0.0;
        relMsEffective = tEff / (0.001 * fs * (double) F);
    }

    double fs = 48000.0;
    int maxCh = 0, tpp = 32, F = 4, maxBlock_ = 0;
    bool prepared_ = false;                     // true only after a fully-successful prepare()
    TruePeakLimiterParams params;

    oversampling::PolyphaseOversampler os;
    std::vector<std::vector<float>>    osBuf;
    std::vector<float*>                osPtrs;
    std::vector<core::DelayLine>       osDelays;
    detail::SlidingMax                 slide;

    double ceilingDb = -1.0, relMsEffective = 100.0;
    int    lookBaseband = 0;
    int    lastNc_ = 0;                         // channel count of the previous process() call
    float  relCoef = 0.0f, grDb = 0.0f;
};

} // namespace felitronics::limiter

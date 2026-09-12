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
#include <cstddef>
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
        // What prepare() asks the heap for (law 11d): two vectors of `cap` entries. The floor is
        // modelled here because prepare() clamps rather than refuses — the budget of a call is
        // storageFor() with the SAME arguments.
        struct Storage
        {
            std::size_t entries = 1;
            std::uint64_t bytes() const noexcept
            {
                return ((std::uint64_t) sizeof (float) + (std::uint64_t) sizeof (std::int64_t))
                     * (std::uint64_t) entries;
            }
        };
        static Storage storageFor (int maxWindow) noexcept
        {
            Storage s; s.entries = (std::size_t) (maxWindow < 1 ? 1 : maxWindow); return s;
        }

        void prepare (int maxWindow)
        {
            const Storage st = storageFor (maxWindow);
            cap = (int) st.entries;
            v.assign (st.entries, 0.0f); ix.assign (st.entries, 0);
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
    int    tapsPerPhase     = oversampling::PolyphaseOversampler::kDefaultTapsPerPhase;   // ≥ 4
};

// PER-BLOCK parameters — safe to change at any time, from the audio thread, mid-stream.
struct TruePeakLimiterParams
{
    double ceilingDbTp = -1.0;       // output true-peak ceiling (dBTP) — see the derate note below
    double releaseMs   = 100.0;
};

//==============================================================================
// THE TRACE, for a caller that needs the SHAPE of what happened rather than one instantaneous number.
// `gainReductionDb()` is a display value: polling it at block ends samples the gain at the host's block
// rate, which is not a property of the signal at all — the same objection `Compressor` records for its
// own meter, one module over, and the same answer (`dynamics::GainReductionTap`).
//
// IT IS ON THE OVERSAMPLED GRID, and that is a decision rather than an implementation leak. The gain is
// decided and applied once per F×fs sample (see processChunk), so F×fs is where the trace exists;
// folding it to baseband needs a RULE, and the rule belongs to whoever reads the statistic.
//
// BE HONEST ABOUT THE SIZE OF THAT, because the first draft of this comment carried an invented number
// and the measurement is an order of magnitude smaller. Taking the MINIMUM over each group of F
// preserves the OVERALL MAXIMUM exactly — that one is a theorem, since the group maximum of the
// magnitude is kept. NOTHING ELSE IS: a fold changes the distribution, so p95 is not preserved in
// general and the second draft of this comment claiming "every upper quantile" was wrong. On the
// measured fixture the two happened to agree (p95 9.525000, max 9.538984 either way), which is one
// fixture and not a proof. What the fold biases measurably is the mean and the active fraction, and
// only upward. Swept over the whole release range on a fixture built to move the gain
// INSIDE a baseband sample (transient pairs, instant attack, 4×): the mean bias is **+0.0029 to
// +0.0080 dB** and the active-fraction bias **+0.0000 to +0.0009**, at release 0.05 ms through 50 ms.
// Small — but it is a choice, it is not zero, and it is not the caller's to discover. A caller that
// wants a baseband curve for a display folds it itself and says which rule it used.
//
// `linkedPeakLin` is THE RECONSTRUCTED PEAK THE LIMITER ACTUALLY SAW — max over channels of |x| on the
// oversampled grid, taken BEFORE the sliding window and before any gain is applied. It is the quantity
// a mastering report means by "the highest inter-sample peak inside the chain", and it cannot be
// recovered from the delivered file: the limiter's whole job is to remove it.
//
// CAPACITY IS BINDING. A non-null buffer whose capacity is short of what the call will produce REFUSES
// the whole call, leaving audio, state and the buffers untouched — the same rule as the compressor's
// tap, and for the same reason: a partial trace looks like data.
struct TruePeakLimiterTap
{
    float* gainReductionDb = nullptr;   // signed dB (<= 0), one per OVERSAMPLED sample
    float* linkedPeakLin   = nullptr;   // linear, one per oversampled sample, BEFORE the gain
    int    capacity        = 0;         // in OVERSAMPLED samples; < numSamples * F refuses the call
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
//     is 2fs/5 at 4× and 8× (+0.436 and +0.108); at 2× fs/3 still costs more (+1.250 against 2fs/5's
//     +0.436), so the allowance is a MAXIMUM over the two rather than a swap of one for the other.
//     RE-DERIVED WHEN tapsPerPhase ROSE TO 64, exactly as the sentence that used to stand here
//     demanded: at 32 taps the prototype's own droop took 0.775 dB off a 0.40 fs tone and so removed
//     2fs/5 before delivery, which left fs/3 (+1.250 / +0.302 / +0.076) worst by default. The pass band
//     is flat there now, so the tone reaches the output and the 4× and 8× figures rise by 0.135 and
//     0.032 dB. Read that correctly: the old numbers were not a property of the limiter, they were a
//     property of its lowpass, and hiding the term is not the same as not having it — `mastering` has
//     run at 64 taps, and therefore at +0.436, since it was written. A sweep of every fs·p/q with
//     q ≤ 64 below 0.46 fs, net of the round-trip droop, confirms 2fs/5 is the maximum and that it
//     holds at every taps count anything here uses (64, 80, 96, 128, 256 all give +0.436 at 4×), but NOT
//     without limit: at 512 taps and 8× the winner becomes 4fs/9 (+0.132 net against 2fs/5's +0.108),
//     because the pass band finally reaches 0.4444 fs. Re-derive when the cutoff moves, and re-derive
//     when the taps grow far enough to deliver a smaller-M tone higher up.
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
// the dense witness reaches +1.87 (4×) / +1.82 (8×), past the figures above. Those two rose with the
// taps default — they were +1.59 and +1.56 at 32 taps — because at both floors the limiter is an OS-rate
// clipper and a longer, sharper FIR rings more re-band-limiting a signal of steps. Both are pinned per
// taps count in the suite. Bringing even that corner
// inside would need floors at 24 baseband samples (0.5 ms), which is a musical setting and would change
// the sound of a legitimate one — so the corner is documented instead of clamped away. The only thing
// proven for every input is the on-grid bound.
//
// SO: a product ceiling needs a DERATE, or a loop closed on measured true peak. Setting ceilingDbTp
// to C does not deliver C dBTP. Over the characterised domain, budget ~1.2 dB for a bright, dense,
// hard-limited master and ~0.3 dB for ordinary material, at 4× or 8× alike — and if you need a
// GUARANTEE rather than a budget, measure the delivered peak and close the loop on it.
//
// ALWAYS IN THE PATH, even when nothing is being limited: the 0.90 × Nyquist prototype is applied TWICE,
// once interpolating and once decimating, so the top of the band is attenuated by |H|² — the dB figures
// of one pass, DOUBLED. Measured on the round trip (and derived independently from designFilter()'s
// coefficients, agreeing to three decimals) at the SHIPPED default of 64 taps, against the 32 it used
// to be:
//
//        f/fs      2× (64)   4× (64)   8× (64)  |  2× (32)   4× (32)   8× (32)
//        0.36     +0.000    +0.000    +0.000   |  −0.001    −0.001    −0.000
//        0.40     +0.000    +0.000    +0.000   |  −0.800    −0.775    −0.762
//        0.42     −0.313    −0.305    −0.301   |  −3.057    −3.017    −2.996
//        0.44     −5.111    −5.091    −5.081   |  −8.092    −8.065    −8.052
//        0.45    −12.041   −12.041   −12.041   | −12.041   −12.041   −12.041
//
// THIS COMMENT USED TO READ −0.40 / −4.0 / −6.0 AND CALL THEM THE ROUND TRIP. Those are one pass. The
// correction matters because the number is a mastering decision, not a rounding error, and a product
// reading it was sizing that decision on half the truth — doubly so with a saturator in front, since a
// second oversampled stage doubles it again: at 44.1 kHz, clipper + limiter cost −1.549 dB at 17.6 kHz
// and −6.033 at 18.5 WHILE THE DEFAULT WAS 32, and +0.000 / −0.610 now. The default moved for a bigger
// reason than the droop, though — at 32 taps the prototype delivered 27 dB of stopband where its own
// Kaiser design declares 90, which is an ALIASING defect; see PolyphaseOversampler.h for the derivation.
// NB it is only APPROXIMATELY factor-independent: exact at 0.45, spread over 0.038 dB at 0.40 (32 taps).
// The suite's own budget check stops at 5fs/14 ≈ 0.357 fs, which is why the doc defect stood
// uncorrected; the frequencies the figures name are measured two-sidedly in the ceiling suite now, and
// the whole tpp × factor surface is pinned by felitronics_oversampling_tests.
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
// prepared is REFUSED — `process()` returns false having touched nothing, because width is a LIMIT
// under law 11(b) and a partially limited block is worse than a refused one. This paragraph used to
// promise the opposite ("limits the first maxChannels and leaves the rest untouched"); law 11 replaced
// that behaviour and the sentence outlived it. Prepare for the most you will ever pass.
//
// RT-safe: prepare() allocates; process() does no alloc/lock/throw and accepts any block length. One
// linked gain for all channels (no image shift). Non-finite and absurd input samples are sanitised at
// the gate — bit-transparently for finite audio within ±1e6.
class TruePeakLimiter
{
public:
    // WHAT prepare() ASKS THE HEAP FOR (law 11d) — the one function it sizes itself with, so a caller
    // budgeting memory reads the counts the scratch, the delay bank and the sliding window are actually
    // built from and the two cannot drift. FALSE, with `out` untouched, exactly where prepare() refuses
    // the same arguments (it IS prepare()'s gate), and a refused prepare() allocates nothing. Asked of a
    // FRESH limiter. `maxBlock` is CLAMPED rather than refused, as prepare() clamps it, and the budget
    // models the clamp: sizing the scratch from the uncapped argument is the defect the cap exists to
    // prevent, and a budget that carried the uncapped number would re-introduce it on paper.
    struct Storage
    {
        std::size_t channels    = 0;       // one oversampled scratch buffer and one delay line per channel
        std::size_t osBufSamples = 0;      // floats in EACH scratch buffer: the capped block x the factor
        int         osDelaySamples = 0;    // the 20 ms lookahead capacity, in OVERSAMPLED samples
        oversampling::PolyphaseOversampler::Storage os {};
        detail::SlidingMax::Storage slide {};
        std::uint64_t bytes() const noexcept
        {
            return (std::uint64_t) channels * ((std::uint64_t) sizeof (std::vector<float>)
                                               + (std::uint64_t) sizeof (float) * (std::uint64_t) osBufSamples
                                               + (std::uint64_t) sizeof (float*))
                 + core::delayBankBytes (channels, osDelaySamples)
                 + os.bytes() + slide.bytes();
        }
    };

    [[nodiscard]] static bool storageFor (double sampleRate, int maxBlock, int maxChannels,
                                          const TruePeakLimiterConfig& config, Storage& out) noexcept
    {
        // Spelled positively so NaN fails: `sampleRate <= 0.0` is FALSE for NaN, which let a NaN rate
        // through to std::lround(NaN) — undefined behaviour, and a signed overflow right after it.
        if (! (sampleRate > 0.0) || ! std::isfinite (sampleRate) || sampleRate > kMaxSampleRate) return false;
        if (maxBlock < 1) return false;
        // Both ends. Only the lower bound was checked, so tapsPerPhase = INT_MAX reached N = L*tpp in
        // the oversampler and overflowed a signed int before allocating.
        if (config.tapsPerPhase < 4 || config.tapsPerPhase > kMaxTapsPerPhase) return false;
        if (! std::isfinite (config.lookaheadMs) || config.lookaheadMs < 0.0) return false;
        // Refused rather than clamped: a silently reduced channel count would leave the surplus channels
        // passing through this module UNLIMITED, and a silently reduced factor would run a topology the
        // caller did not ask for.
        if (maxChannels < 1 || maxChannels > core::kMaxChannels) return false;
        if (config.oversampleFactor > kMaxFactor) return false;
        Storage st;
        const int f = oversampleFactorFor (config);
        if (! oversampling::PolyphaseOversampler::storageFor (f, maxChannels, config.tapsPerPhase, st.os))
            return false;
        // A rate so low that 20 ms cannot hold the minimum lookahead would make prepare()'s clamp
        // std::clamp(x, 2, 1) — lo > hi is undefined behaviour. Refuse instead, here and there.
        const int maxLookBb = maxLookaheadSamplesFor (sampleRate);
        if (maxLookBb < kMinLookaheadSamples) return false;
        st.channels       = (std::size_t) maxChannels;
        st.osBufSamples   = (std::size_t) blockFor (maxBlock) * (std::size_t) f;
        st.osDelaySamples = maxLookBb * f;
        st.slide          = detail::SlidingMax::storageFor (st.osDelaySamples + 1);
        out = st;
        return true;
    }

    // THE EFFECTIVE FACTOR, in one place: a requested 1 becomes 2, and anything above kMaxFactor was
    // refused before this is asked. `MasteringChain` reads this to size a tap buffer without a limiter in
    // hand, and reading it back off a prepared limiter has to give the same answer (pinned).
    static int oversampleFactorFor (const TruePeakLimiterConfig& config) noexcept
    {
        return config.oversampleFactor < 2 ? 2 : config.oversampleFactor;
    }

    // The scratch block prepare() will actually use — the cap, in one place, so the budget cannot size
    // itself from the uncapped argument the preparation ignores.
    static int blockFor (int maxBlock) noexcept { return std::min (maxBlock, kMaxBlock); }

    // The 20 ms capacity in BASEBAND samples — what bounds the lookahead below and sizes the rings.
    static int maxLookaheadSamplesFor (double sampleRate) noexcept
    {
        return (int) std::ceil (kMaxLookaheadMs * 0.001 * sampleRate);
    }

    // THE LOOKAHEAD a config actually takes at this rate, after the 20 ms ceiling and the two-sample
    // floor — prepare()'s own clamp chain, in ONE place.
    static int lookaheadSamplesFor (double sampleRate, double lookaheadMs, int maxLookBaseband) noexcept
    {
        // Clamp in the DOUBLE domain BEFORE lround: `std::lround(1e300)` is out of range for a long, and
        // the result of that then clamped upward landed on the MINIMUM lookahead — a caller asking for an
        // absurd value silently got the smallest one instead of the largest.
        const double lookMs = std::clamp (lookaheadMs, 0.0, kMaxLookaheadMs);
        const int    n      = (int) std::lround (lookMs * 0.001 * sampleRate);
        return std::clamp (n, kMinLookaheadSamples, maxLookBaseband);
    }

    // THE LATENCY a prepared limiter will report for this geometry, without preparing one: the
    // oversampler's round trip plus the lookahead. `MasteringChain` sizes its dry aligner from the
    // latency it READS BACK off the stage, and this is how a budget reaches the same number without a
    // second derivation of it. 0 where prepare() refuses the same arguments.
    static int latencyFor (double sampleRate, int maxBlock, int maxChannels,
                           const TruePeakLimiterConfig& config) noexcept
    {
        Storage st;
        if (! storageFor (sampleRate, maxBlock, maxChannels, config, st)) return 0;
        return (config.tapsPerPhase - 1)
             + lookaheadSamplesFor (sampleRate, config.lookaheadMs, maxLookaheadSamplesFor (sampleRate));
    }

    // Topology is chosen HERE and nowhere else. Returns false and leaves the limiter unprepared on an
    // unusable stream configuration; a false return is the only way to learn that, so check it.
    [[nodiscard]] bool prepare (double sampleRate, int maxBlock, int maxChannels, const TruePeakLimiterConfig& config = {})
    {
        prepared_ = false;                                     // any early return below leaves it unprepared
        // THE GATE IS storageFor()'s, so the budget and the preparation refuse the same arguments by
        // construction rather than by agreement — every check that used to stand here is there, with its
        // reasons, and nothing has been dropped: the positive NaN spelling, both ends of tapsPerPhase, the
        // channel and factor refusals, and the 20 ms floor that would otherwise make the lookahead clamp
        // std::clamp(x, 2, 1).
        //
        // ONE OBSERVABLE CHANGE COMES WITH THAT: `fs` used to be written between the rate check and the
        // width check, so a preparation refused on its WIDTH had already replaced the rate, and
        // `effectiveReleaseMs()` — which is not gated by `prepared_`, and is the ONLY readout here that is
        // a function of the rate — then answered in terms of it. Measured, `prepare(96000, 64, 2)` then a
        // refused `prepare(48000, 64, 33)` then `setParams(releaseMs = 0.01)`: 0.166667 ms before, 0.083333
        // now, the latter being the release the limiter is actually running. (`effectiveCeilingDbTp()` is
        // NOT affected — the ceiling is a clamp on the parameter and never sees the rate.) The new answer is
        // law 11(b)'s; it is stated here because this header is shared with the plug-ins.
        Storage st;
        if (! storageFor (sampleRate, maxBlock, maxChannels, config, st)) return false;

        fs    = sampleRate;
        maxCh = maxChannels;
        // CLAMPED, not refused. Since process() chunks, maxBlock is a scratch-buffer size rather than a
        // limit on what a caller may pass — and refusing a big one would fail prepare(), after which
        // process() returns the buffer UNTOUCHED, i.e. exactly the unlimited-passthrough defect this
        // task exists to close. An offline caller sizing maxBlock to a whole file is a normal thing to do.
        maxBlock_ = blockFor (maxBlock);
        tpp   = config.tapsPerPhase;
        F     = oversampleFactorFor (config);                  // a requested 1 becomes 2; above kMaxFactor was refused
        if (! os.prepare (F, maxCh, tpp)) return false;        // oversampler rejected → stay unprepared

        // st.osBufSamples is maxBlock_ x F, NOT maxBlock x F: the cap above exists to bound exactly this
        // allocation ("without a cap a hostile or mistaken prepare() could ask for gigabytes"), and sizing
        // the buffer from the UNCAPPED argument left it doing nothing. It is reachable by ordinary use,
        // not only by a hostile one — this header tells an offline caller that sizing maxBlock to a whole
        // file is a normal thing to do, and a 10-minute stereo file at 48 kHz asked for 460 MB per channel
        // instead of the 16 MB the cap allows. processChunk() already works in maxBlock_-sized pieces.
        //
        // NO TEMPORARY: `assign(maxCh, std::vector<float>(n))` built one buffer to copy maxCh times, and
        // that temporary alone was 524 288 B of peak nobody could see in the end state on the biggest
        // topology the MASTERING CHAIN builds (its quantum stops at 8192). Asked of this class directly the
        // same temporary reached 64 MiB, since the block cap above is 1 Mi samples. Resizing the outer vector and filling each buffer in place is the same end state, one
        // allocation class fewer, and — the point — a re-preparation at the same geometry now asks for
        // nothing at all.
        // `reserve` BEFORE `resize`, and that is not decoration: a `resize` past the current capacity grows
        // GEOMETRICALLY (libc++ doubles), so widening a bank of 3 buffers to 4 asked for 6 and the published
        // number — which pays for 4 — was no longer an upper bound on what the call held at once. `reserve`
        // asks for exactly what is needed. (The code-review round; measured 48 B over, 2 x sizeof(vector).)
        osBuf.reserve (st.channels);
        osBuf.resize  (st.channels);
        for (auto& b : osBuf) b.assign (st.osBufSamples, 0.0f);
        osPtrs.assign (st.channels, nullptr);

        const int maxLookBb  = maxLookaheadSamplesFor (fs);
        const int maxLookOS  = st.osDelaySamples;
        core::prepareDelayBank (osDelays, st.channels, maxLookOS);
        slide.prepare (maxLookOS + 1);

        // The lookahead floor is not taste. At zero the structure degenerates into a clipper at the
        // oversampled rate: +2.05 dB over the ceiling, measured, and unlike the terms above that one has
        // no closed form to derate against. Two baseband samples is where the witness matrix comes back
        // inside the envelope with the release at a musical value, at both conformance rates and every
        // factor. (One sample is already inside at a 50 ms release; two is the value that also holds the
        // click at the release floor, which is why the pair is stated together.) At any musical lookahead
        // this clamp is invisible: it is 42 µs at 48 kHz.
        // THROUGH the static, not beside it: `latencyFor()` publishes this number to a caller that has no
        // limiter yet, and two spellings of one clamp chain is the drift law 11d's budgets exist to make
        // impossible.
        lookBaseband = lookaheadSamplesFor (fs, config.lookaheadMs, maxLookBb);
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
        linkedPeakLin_ = 0.0f;                                 // a free-running maximum SINCE RESET, like
                                                               // TruePeakMeter's — so it means "this stream"
        lastNc_ = 0;
    }

    void setParams (const TruePeakLimiterParams& p) noexcept { params = p; apply (p); }

    // All of these read as "nothing prepared" after a failed prepare(), rather than reporting the
    // topology of whatever was prepared before it — a stale latency is worse than an obvious zero.
    int    latencySamples()  const noexcept { return prepared_ ? os.latencySamples() + lookBaseband : 0; }
    double gainReductionDb() const noexcept { return grDb; }
    bool   isPrepared()      const noexcept { return prepared_; }

    // The highest RECONSTRUCTED peak this limiter has seen since reset() — the channel-linked maximum of
    // |x| on the F x fs grid, taken before the sliding window and before any gain. Free-running, exactly
    // like `analysis::TruePeakMeter::truePeakDb()`, and for the same reason: a maximum with a ballistic
    // on it is a display, not a measurement.
    //
    // IT IS NOT "THE" TRUE PEAK OF THE INPUT, and the difference is the point of reporting it separately.
    // This is what the limiter's OWN reconstruction saw — a 0.90 x Nyquist Kaiser prototype at this
    // instance's factor and tapsPerPhase — and it is the number that EXPLAINS the gain reduction this
    // instance applied. A meter of a different design reads something else on the same signal, and the
    // gap is the material's, not a defect: against `analysis::TruePeakMeter` (the spec's 12-tap filter)
    // at 4x, measured, a 1 kHz burst train agrees to **-0.0012 dB** and a pair of adjacent full-scale
    // impulses — maximally broadband, i.e. the worst case for two different low-passes — disagrees by
    // **-1.4183 dB**. Report it as the limiter's reconstruction, never as the file's true peak, and use
    // `analysis::TruePeakMeter` for the latter. It cannot be recovered from the delivered file at all:
    // removing it is the limiter's whole job.
    double maxReconstructedPeakDb() const noexcept { return core::gainToDb ((double) linkedPeakLin_); }

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
    // Law 11 (DSP-ARCHITECTURE.md §2). The width used to be CLAMPED, and the header above documented the
    // surplus channels leaving UNLIMITED as if it were a contract. Measured: prepared for 2, called with
    // 4, ceiling -1 dBFS, input +6.02 dBFS — the surplus came out at +6.02, i.e. +7.02 dB over the ceiling
    // it was told to hold (unbounded in general: it is the caller's own input, untouched). Worse than the
    // level, the two processed planes carry latencySamples() = 111 of oversampler + lookahead that the two
    // untouched ones do not, so a fold-down combs at fs/(2*111) = 216 Hz — a fault that sounds like a
    // timbre rather than like a fault. Refused whole, before anything moves.
    [[nodiscard]] bool process (float* const* channels, int numChannels, int numSamples) noexcept
    {
        return process (channels, numChannels, numSamples, TruePeakLimiterTap {});
    }

    // The full form: the same call, with the oversampled gain-reduction and reconstructed-peak traces
    // written out. `tap.gainReductionDb == nullptr && tap.linkedPeakLin == nullptr` is off and costs
    // nothing; a non-null tap with `tap.capacity < numSamples * oversampleFactor()` REFUSES the whole
    // call, before anything moves. RT-safe.
    [[nodiscard]] bool process (float* const* channels, int numChannels, int numSamples,
                                TruePeakLimiterTap tap) noexcept
    {
        if (numChannels < 0 || numSamples < 0) return false;
        if (! prepared_) return false;
        if (numChannels > maxCh) return false;                 // width is a LIMIT — law 11(b)
        // Checked BEFORE anything moves, so a refused call is indistinguishable from one never made.
        // The multiplication is in `long long` on purpose: `numSamples * F` overflows a signed int at
        // 537 million samples per channel at 4x, which a whole-file offline call can reach, and the
        // overflow would make a SHORT buffer compare as large enough.
        if ((tap.gainReductionDb != nullptr || tap.linkedPeakLin != nullptr)
            && (long long) tap.capacity < (long long) numSamples * (long long) F) return false;
        if (numSamples == 0) return true;
        const int nc = numChannels;

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
        if (nc == 0) return true;                              // law 11(d): the reset above IS the edge

        float* sub[core::kMaxChannels] {};
        for (int off = 0; off < numSamples; )
        {
            const int n = std::min (numSamples - off, maxBlock_);
            for (int c = 0; c < nc; ++c) sub[(std::size_t) c] = channels[c] + off;
            // The tap advances on the OVERSAMPLED clock, so the chunk loop has to step it by n*F and
            // not by n. Spelled as a separate `TruePeakLimiterTap` rather than by mutating the caller's
            // copy, so the capacity check above stays the statement about the WHOLE call that it is.
            TruePeakLimiterTap sTap;
            const std::size_t osOff = (std::size_t) off * (std::size_t) F;
            if (tap.gainReductionDb != nullptr) sTap.gainReductionDb = tap.gainReductionDb + osOff;
            if (tap.linkedPeakLin   != nullptr) sTap.linkedPeakLin   = tap.linkedPeakLin   + osOff;
            processChunk (sub, nc, n, sTap);
            off += n;                                          // `off += maxBlock_` could step past INT_MAX
        }
        return true;
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

    void processChunk (float* const* channels, int nc, int numSamples, TruePeakLimiterTap tap) noexcept
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

            // BEFORE the sliding window and before any gain: this is the reconstructed peak the limiter
            // saw, which is what a mastering report means and what the delivered file no longer holds.
            if (linkedPeakLin_ < linkedPeak) linkedPeakLin_ = linkedPeak;
            if (tap.linkedPeakLin != nullptr) tap.linkedPeakLin[(std::size_t) i] = linkedPeak;

            const float  smax    = slide.push (linkedPeak);
            const double smaxDb  = core::gainToDb (smax);
            double rawRedDb = ceilingDb - smaxDb;
            if (rawRedDb > 0.0) rawRedDb = 0.0;

            grDb = std::min ((float) rawRedDb, grDb * relCoef);   // instant attack, exponential release toward 0 dB
            if (tap.gainReductionDb != nullptr) tap.gainReductionDb[(std::size_t) i] = grDb;
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
    // Read only after a successful prepare() overwrites it, so this is a seed rather than a policy —
    // but a seed that disagrees with the shipped default is the next reader's wrong answer about what
    // the default IS, which is the whole subject of this class of defect.
    int maxCh = 0, tpp = oversampling::PolyphaseOversampler::kDefaultTapsPerPhase, F = 4, maxBlock_ = 0;
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
    float  linkedPeakLin_ = 0.0f;               // free-running max of the reconstructed linked peak
};

} // namespace felitronics::limiter

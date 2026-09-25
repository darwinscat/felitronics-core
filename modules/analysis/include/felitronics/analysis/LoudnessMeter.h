// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/core/Config.h>
#include <felitronics/core/DetMath.h>
#include <felitronics/core/Math.h>
#include <felitronics/analysis/KWeightingFilter.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <span>
#include <vector>

namespace felitronics::analysis
{

//==============================================================================
// felitronics::analysis::LoudnessMeter — ITU-R BS.1770-4 / EBU R128 loudness: MOMENTARY (400 ms),
// SHORT-TERM (3 s) and INTEGRATED (gated). LUFS = -0.691 + 10·log10(Σ weightₖ·meanSquareₖ) of the
// K-weighted signal.
//
// Energy is accumulated in 10 ms SUB-HOPS: momentary is the mean of the last 40, short-term of the last 300,
// so either window lands within 10 ms of any event. EBU Tech 3341's file-based cases 10 and 13 slide a 3 s /
// 400 ms tone in 150 ms / 20 ms steps and expect the maximum reading to be the tone ±0.1 LU at EVERY offset —
// a window that only moves in 100 ms steps cannot do that (a 400 ms burst 40 ms off its grid reads 0.45 LU
// low). The integrated measure keeps its 100 ms HOP: every tenth sub-hop closes a 400 ms gating block
// (75 % overlap), and integrated = the gated mean over those blocks — absolute gate at -70 LUFS, then a
// -10 LU relative gate (a two-pass over all absolute-gated blocks — the threshold moves as more program
// arrives, so it must NOT be a one-pass running sum).
//
// RT-safe: prepare() — or prepareForSamples() — allocates the sub-hop ring + the integrated-block buffer;
// process() only indexes them (no alloc/lock/throw). The gated measure keeps every block's energy to the end
// (the relative gate is a two-pass), so the block store is sized for the prepared capacity — maxDurationSec at
// the prepared rate, or a count of samples; blocks past it are counted in droppedBlocks() and not kept — size
// for the longest program and that reads 0. Channel
// weights default to 1.0 (correct for mono/stereo); set them per the BS.1770 roles (Ls/Rs = 1.41, LFE
// excluded) for surround — the host-layout→role mapping is product glue.
template <class Math = core::SystemMath>
class BasicLoudnessMeter
{
public:
    // THE LOWEST RATE THIS METER MEASURES AT — the core's floor (P103; the entries that call it have had it since
    // P51). Below twice the K-weighting shelf the filter is aliased and, in most of that range, unstable: a 0 dBFS
    // 400 Hz sine read +3043 LUFS at 3300 Hz. A rate <= 0 or NaN is not "wrong", it is "not given", and still reads
    // as 48 kHz — the published default; a rate in (0, kMinSampleRate) is refused. storageFor() decides it.
    static constexpr double kMinSampleRate = core::kMinSampleRate;

    // SECONDS are the convenience; the store is counted in SAMPLES — see prepareForSamples(). A NaN or negative
    // duration reads as 0 s, exactly as it always has (std::max keeps its first argument on NaN).
    [[nodiscard]] bool prepare (double sampleRate, int numChannels, double maxDurationSec = 3600.0)
    {
        const double rate = sampleRate > 0.0 ? sampleRate : 48000.0;
        return prepareForSamples (sampleRate, numChannels, std::max (0.0, maxDurationSec) * rate);
    }

    // THE SAME PREPARATION, SIZED IN SAMPLES — the unit the store is actually counted in. A whole-programme
    // caller knows its length in frames. (A trip through seconds used to be able to lose it: `frames / fs` was
    // +inf at a finite rate this meter accepted before P103 — 1e-305 Hz, 2000 frames — and the store was sized
    // from `(std::size_t) inf`, undefined behaviour that kept 3 blocks on arm64 and wasm32 and 4 on x86-64 gcc
    // for the same call. At the floor `frames / fs` is at most INT_MAX / 8000 s.)
    //
    // A REFUSED CALL LEAVES NOTHING OF THE PREVIOUS PROGRAMME READABLE (law 11b), whichever check refused it:
    // `reset()` runs FIRST, before any of them, so every reading below answers what a never-prepared meter
    // answers. (Before P103 a refusal only cleared the flag, and momentaryLufs(), integratedLufs(),
    // gatingBlockEnergies() and droppedBlocks() went on serving the previous programme.)
    [[nodiscard]] bool prepareForSamples (double sampleRate, int numChannels, double maxSamples)
    {
        prepared_ = false;
        reset();
        const double rate = sampleRate > 0.0 ? sampleRate : 48000.0;    // fs<=0 → subSamples 0 → /0 in finishSubHop
        if (numChannels < 1 || numChannels > kMaxChannels) return false;   // law 11(b): BINDING
        Storage st;
        if (! storageFor (sampleRate, maxSamples, st)) return false;    // law 11(b): validate, THEN write
        fs = rate;
        ch = numChannels;
        kw.prepare (fs, ch);
        subSamples = st.subSamples;
        for (int c = 0; c < kMaxChannels; ++c) w[c] = 1.0;
        subRing.assign (kSubRing, 0.0);
        blockE.assign (st.blocks, 0.0);
        stE.assign (st.shortTerm, 0.0);
        reset();
        prepared_ = true;
        return true;
    }

    // WHAT prepareForSamples() ALLOCATES — the one function it sizes itself with, so a caller budgeting memory
    // reads the numbers the store is built from and the two cannot drift. Every byte of the meter's heap is here:
    // the three vectors below and nothing else (the K-weighting state is inline) — asked of a FRESH meter; a prepared
    // one keeps whatever storage still fits.
    struct Storage
    {
        int         subSamples = 0;     // one 10 ms sub-hop, lround (0.01·fs) and at least 1
        std::size_t blocks     = 0;     // 400 ms gating blocks kept for the integrated measure
        std::size_t shortTerm  = 0;     // 3 s short-term samples kept for LRA, one per 100 ms hop (10 Hz)
        // 64 bits because the product is the point: on wasm32 a `size_t` byte count wraps long before the element
        // counts above do.
        std::uint64_t bytes() const noexcept
        {
            return (std::uint64_t) sizeof (double) * ((std::uint64_t) kSubRing + blocks + shortTerm);
        }
    };

    // Sized by HOPS at this rate, not by seconds: a hop is 10 × lround (0.01·fs) samples, which is 100 ms only
    // where fs is a multiple of 100 — elsewhere a per-second count drifts from the blocks that actually arrive.
    // The +4 blocks / +8 short-term samples are MARGIN, not need: the first gating block is born on the 4th hop
    // and the first short-term sample on the 30th, so a store sized for exactly n samples already holds every
    // block n samples produce — LoudnessConformanceTests pins the margin by the count an OVERFEED drops, and pins
    // the bytes by a table.
    //
    // FALSE, with `out` untouched, when the capacity is not REPRESENTABLE. "Finite" is not the property the
    // arithmetic needs. The hop, ten sub-hops, is an `int`, so the sub-hop lround (0.01·fs) is bounded by
    // kMaxSubHop — and the test is made on the COMPUTED 0.01·fs, not on fs, because near the edge that product
    // rounds up to .5 and lround() goes with it. And every block index is an `int`, so the count is bounded by
    // kMaxBlocks. A rate <= 0 or NaN is read as 48 kHz, exactly as prepare() reads it: the budget of a call is
    // storageFor() with the SAME arguments. A rate in (0, kMinSampleRate) is refused, and it is decided HERE and
    // only here, BEFORE that substitution — prepare() and prepareForSamples() refuse through this function.
    [[nodiscard]] static bool storageFor (double sampleRate, double maxSamples, Storage& out) noexcept
    {
        if (sampleRate > 0.0 && sampleRate < kMinSampleRate) return false;   // NaN, <= 0 and +inf are not in it
        const double rate   = sampleRate > 0.0 ? sampleRate : 48000.0;   // prepare()'s own substitution
        const double subHop = 0.01 * rate;                      // a sub-hop is lround (subHop) samples
        if (! (subHop < (double) kMaxSubHop + 0.5)) return false; // lround <= kMaxSubHop, so 10 of them fit an int
        if (! (maxSamples >= 0.0)) return false;              // NaN and -inf fail here; +inf fails the block bound below
        const int s = std::max (1, (int) std::lround (subHop));
        const double hops = std::ceil (maxSamples) / (double) (s * kSubHopsPerHop);
        if (! (hops <= (double) kMaxBlocks)) return false;
        out.subSamples = s;
        out.blocks     = (std::size_t) hops + 4;
        // ONE PER HOP since the cadence became 10 Hz, where it used to be one per ten. Ten times the entries
        // and ten times this term of the demand: for a ten-minute programme it is 48 KB where it was 4.8.
        out.shortTerm  = (std::size_t) std::ceil (hops) + 8;
        return true;
    }

    void reset() noexcept
    {
        kw.reset();
        for (int c = 0; c < kMaxChannels; ++c) subSumSq[c] = 0.0;
        ranNc_ = 0;                             // nothing has run, so nothing can be stopping
        subCount = 0; subWrite = 0; subFilled = 0; subInHop = 0; blockCount = 0; droppedBlocks_ = 0;
        nonFiniteSubHops_ = 0;
        stCount = 0; droppedShortTerm_ = 0;
        std::fill (subRing.begin(), subRing.end(), 0.0);
    }

    // Non-finite is refused, not stored: `w[c] * energy` would make every sub-hop non-finite from a
    // CONFIGURATION mistake, which no amount of healing downstream can undo. House rule for params.
    void setChannelWeight (int c, double weight) noexcept
    {
        if (c >= 0 && c < kMaxChannels && std::isfinite (weight)) w[c] = weight;
    }

    // Law 11 (DSP-ARCHITECTURE.md §2). NB the sub-hop counter runs on the SAMPLES, not on the channels,
    // so a zero-width call still spends measurement time — that was already true and is law 11(d).
    [[nodiscard]] bool process (const float* const* channels, int numChannels, int n) noexcept
    {
        if (numChannels < 0 || n < 0) return false;
        if (! prepared_) return false;                                 // unprepared — subRing/blockE/stE empty
        if (numChannels > ch) return false;                            // width is a LIMIT — law 11(b)
        if (n == 0) return true;                                      // law 11(d): the ONE true no-op —
                                                                      // and the falling edge below is
                                                                      // clocked by `n`, exactly as 11a says
        const int nc = numChannels;
        // THE SIBLING OF TruePeakMeter'S HISTORY. `KWeightingFilter` holds a per-channel TDF-II state,
        // and a channel that stops and comes back is weighted against audio from before the gap: the
        // RLB shelf's tail spills out on return. Measured — a 60 Hz tone, one second of zero-width
        // calls, then 400 ms of DIGITAL SILENCE — momentary read -29.19 LUFS where real silence reads
        // -120.00. Same class, same answer, one file over.
        //
        // ONLY THE FILTER. `subSumSq` is not stale memory, it is MEASUREMENT ALREADY TAKEN in the
        // sub-hop now in progress, and clearing it here threw that away: one zero-width call in the
        // middle of a tone moved the reading by -3.02 LU, a single zero-width sample turned -9.11 LUFS
        // into -120, and the evidence of a non-finite input (`nonFiniteSubHops()`) went from 1 to 0.
        // A falling edge drops what a cell REMEMBERS, never what a meter has already counted.
        for (int c = nc; c < ranNc_; ++c) kw.resetChannel (c);
        ranNc_ = nc;
        for (int i = 0; i < n; ++i)
        {
            for (int c = 0; c < nc; ++c) { const double y = kw.process (c, (double) channels[c][i]); subSumSq[c] += y * y; }
            if (++subCount >= subSamples) finishSubHop (nc);
        }
        return true;
    }

    double momentaryLufs()  const noexcept { return lufsOf (meanLastSubHops (kMomentarySubHops)); }  // 400 ms
    double shortTermLufs()  const noexcept { return lufsOf (meanLastSubHops (kSubRing)); }           // 3 s
    double integratedLufs() const noexcept { return integrated(); }                                  // gated
    double loudnessRangeLu() const noexcept { return lra(); }                                        // EBU Tech 3342 (P95−P10)
    // How many 3 s short-term samples the range was read from. Published because the CADENCE is part of what
    // Tech 3342 3.1 specifies — "a minimum block overlap of 2.9 s … i.e. >=10 Hz sampling" — and a count is
    // the only way a caller, or a test, can see which cadence a build runs.
    int shortTermCount() const noexcept { return stCount; }

    // Gating blocks that arrived past the prepared capacity and were not kept. Non-zero means integratedLufs()
    // describes only the part of the program that fitted — a caller that must not lose a block sizes prepare() /
    // prepareForSamples() for its longest program and checks this reads 0.
    //
    // THE SENTENCE USED TO NAME loudnessRangeLu() HERE TOO, and that was an overstatement worth spelling out:
    // the range is read from a SECOND store, the two overflow at different moments, and a non-zero count here
    // does NOT mean the range was truncated. Blocks: capacity floor(hops) + 4 against a production of
    // floor(hops) - 3, so an overfeed of 8 hops loses the first one. Short-term: capacity ceil(hops) + 8
    // against floor(hops) - 29, so it takes 38. Thirty hops — three seconds — where this counter climbs and
    // the range is still whole. Integers, not rates: LoudnessConformanceTests feeds exact hop multiples at
    // 48 and 44.1 kHz and gets the same two thresholds.
    //
    // So the block store always goes first, and a caller that sizes for its longest programme and checks THIS
    // reads 0 does have both answers whole — the half of the old sentence that was true. What it could not do
    // before droppedShortTermSamples() existed was tell the two apart once this one was non-zero.
    int droppedBlocks() const noexcept { return droppedBlocks_; }

    // Short-term samples that arrived past the prepared capacity and were not kept — the same statement for
    // loudnessRangeLu() that droppedBlocks() makes for integratedLufs(), and until K12 nothing made it. A
    // truncated range is not a smaller range, it is a different number: the percentiles are taken over whatever
    // survived. The margin this counter watches NARROWED with the cadence — at one sample a second the store
    // held floor(hops/10) + 8 (a truncating cast, not a ceiling) against a production of about hops/10 - 2,
    // which took an overfeed of about a hundred hops; at one a hop it is ceil(hops) + 8 against
    // floor(hops) - 29, and 38 hops are enough. Both were silent; this is the counter that is not.
    int droppedShortTermSamples() const noexcept { return droppedShortTerm_; }

    // Completed 10 ms sub-hops whose channel-weighted mean square came out NON-FINITE — the damage counter,
    // and the reason the readings above may be believed or may not. STICKY until reset()/prepare(): a live
    // display recovers within a few hundred ms (the K-weighting state is healed at each sub-hop boundary),
    // but the whole-programme integrated and LRA answers are compromised for good, so anything non-zero
    // means "this number is best effort, not a measurement".
    //
    // Why the SUB-HOP is the unit and not the 400 ms gating block: the sub-hop is the first energy boundary
    // that does not depend on how the caller chunks its calls, and one poisoned sub-hop fans out into
    // exactly FOUR overlapping gating blocks (the 400 ms window advances every 100 ms), so a block counter
    // would report one event four times. It counts non-finite ENERGY, not bad input samples — a non-finite
    // channel weight would land here too, and healing the filter cannot heal a poisoned configuration.
    //
    // uint64 and saturating: at 100 sub-hops a second a plain int overflows — UB — in about 249 days of
    // continuous poison, and the one invariant this counter must keep is that non-zero never becomes zero
    // again before reset().
    std::uint64_t nonFiniteSubHops() const noexcept { return nonFiniteSubHops_; }

    // The 400 ms gating blocks' K-weighted mean-square energies, in arrival order, exactly as finishHop()
    // recorded them — BEFORE either gate. Always FINITE: a sub-hop whose K-weighted energy came out
    // non-finite is recorded as silence and counted in nonFiniteSubHops(), so this vector stays a portable
    // bit-comparison surface (a computed NaN would not be — its sign bit differs between arm64 and x86-64). This is what a cross-toolchain bit-exactness check compares:
    // integratedLufs() is DISCONTINUOUS in these energies (a block landing within ~1e-12 of a gate flips its
    // inclusion and moves the reading by ~0.01 dB), so it cannot carry a bit-identity claim, while the vector
    // itself is continuous in the input samples and can. Raw energy, not dB, deliberately — a dB accessor
    // would route the comparison back through log10, whose 1-ulp disagreement between libms is precisely what
    // such a check must not inherit.
    //
    // The span is valid until the next prepare() or prepareForSamples() (the only places the storage is
    // reallocated); process() only appends, and reset() zeroes the count rather than the storage. Empty before
    // either.
    //
    // WHAT THIS PINS. It commits the meter to keeping every block's energy as a contiguous array of double in
    // arrival order. Retention itself is already forced by the two-pass relative gate, so the new constraint
    // is only the storage shape — but it does stand in the way of the one plausible future change here, a
    // libebur128-style HISTOGRAM mode for unbounded streams (which is the real answer to droppedBlocks()).
    // Should that arrive, this accessor returns an empty span in histogram mode rather than pretending.
    std::span<const double> gatingBlockEnergies() const noexcept
    {
        return { blockE.data(), (std::size_t) blockCount };
    }
    int gatingBlockCount() const noexcept { return blockCount; }

private:
    static constexpr int kMaxChannels      = core::kMaxChannels;
    static constexpr int kSubHopsPerHop    = 10;                                    // 10 × 10 ms = the 100 ms gating hop
    static constexpr int kMomentarySubHops = 40;                                    // 400 ms
    static constexpr int kSubRing          = 300;                                   // 3 s — the short-term window, and the ring
    // The bounds storageFor() refuses beyond, each named for the property it protects.
    static constexpr int kMaxSubHop        = std::numeric_limits<int>::max() / 10;  // ten sub-hops make a hop, and a hop is an int
    // THE BOUND IS THE TIGHTER OF THE TWO STORES, and until K12 it was written for the looser one. It read
    // INT_MAX - 4, sized for `blocks = floor(hops) + 4`, because the short-term store was a tenth of the hop
    // count and could not be the binding constraint. At one sample a hop it is `ceil(hops) + 8`, so INT_MAX - 4
    // ACCEPTS a capacity of INT_MAX + 4 — and `stCount < (int) stE.size()` then compares against a narrowed
    // negative, rejecting every short-term sample from the first one on, with `droppedBlocks()` still reading 0.
    // Unreachable in practice (the store would be 17 GB) but reachable on paper, which is where a bound lives.
    static constexpr int kMaxBlocks        = std::numeric_limits<int>::max() - 8;   // both stores' indices are ints

    void finishSubHop (int nc) noexcept
    {
        // LAW 8 lives HERE, not at the end of process(), and the difference is not cosmetic. This boundary
        // is a deterministic 10 ms of AUDIO (lround(0.01*fs) samples); the end of process() is wherever the
        // caller happened to cut the stream. Flushing there would make the numbers depend on the host's
        // block size, and the family claims — and tests — that they do not (here within 1e-9 by
        // LoudnessConformanceTests' "chunk invariance"; bit-exact across call sizes 1 … 100 003 by
        // felitronics-mastering-core's ProbeTests, which drives this meter through its measurement probe). It also failed outright at low rates: this meter took any
        // rate before P103 (the probe took 1 kHz before P51), and at 1 kHz a single 8192-sample call spans 8.2 s. And the
        // interval to beat is not the RLB's 2.85 s but the SHELF's 90 ms (4 324 samples at 48 kHz) — 8192
        // samples is already 170 ms, so a per-host-block flush would let the shelf sit subnormal for half of
        // every silent block. Here the arrears can never exceed 10 ms of audio at any rate, and the cost is
        // ~100 flushes a second.
        kw.flushDenormals();

        // Non-finite energy is caught HERE — the first place it exists — per CHANNEL, and the poisoned
        // channel's 10 ms is recorded as silence rather than propagated.
        //
        // Per channel, not on the sum, because BS.1770 gives LFE w = 0 and `0 * NaN` is NaN: one excluded
        // channel would otherwise poison a sub-hop whose audible channels were perfectly fine.
        //
        // Recorded as 0.0 and NOT kept raw, which is the part that looks wrong and is not. Keeping the
        // poison would read as the honest choice — gatingBlockEnergies() is the forensic surface — but that
        // surface is ALSO the tier's cross-toolchain bit-exactness comparison, and a computed NaN is not
        // portable: `inf - inf` (which is what +inf input produces inside the biquad) is 0x7ff8000000000000
        // on arm64 and 0xfff8000000000000 on x86-64 — measured on this Mac against Debian/gcc 14.2, and the
        // sign bit alone makes every block line differ between machines. A NaN that ARRIVED as input carries
        // its payload identically on both, which is exactly why the old pinned test never showed this.
        // Zero is canonical, and the event is not lost: it is in nonFiniteSubHops(). The cost is bounded and
        // derivable — a zeroed sub-hop is 1/40 of a gating block, so 10*log10(39/40) = -0.11 dB on that
        // block, and nothing on the rest of the programme.
        double subMS = 0.0; bool poisoned = false;
        for (int c = 0; c < nc; ++c)
        {
            const double e = subSumSq[c];              // NaN*NaN or inf*inf — any bad sample lands here
            // AN EXCLUDED CHANNEL CANNOT POISON A MEASUREMENT IT IS NOT IN. `nonFiniteSubHops()` is
            // documented as counting sub-hops whose CHANNEL-WEIGHTED mean square came out non-finite,
            // and the code counted them unweighted: a NaN in a `w = 0` channel — which is exactly what
            // BS.1770 gives LFE — flagged a sub-hop whose weighted energy was perfectly fine. The
            // energy path next door already reasons this way, and reasoned it first; the counter had
            // simply not been made to agree with it. That mattered the moment the counter became
            // load-bearing: `TargetLoudnessSolver::measureInputLoudnessRange()` (felitronics-mastering-core) REFUSES on it, so
            // an unweighted count would refuse a 5.1 programme over a channel the standard excludes.
            //
            // EXCLUDED MEANS EXACTLY ZERO, and the comparison is exact ON PURPOSE — `core::exactlyEqual`
            // is this repository's own name for that, written so core headers stay warning-clean under
            // strict downstream flags (`-Wfloat-equal` is one, and it is an ERROR on the GCC rows). An
            // epsilon would be WRONG here, not merely loose: it would silently drop a channel somebody
            // weighted at 1e-12, and the weight is a caller's statement about the layout, not a
            // measurement. `w[c] > 0.0` would be wrong too, and that is not obvious — the SOLVER's
            // setter refuses a negative weight, this one refuses only a NON-FINITE one
            // (`setChannelWeight` above), so a negative weight reaches here and it DOES contribute to
            // the sum below. "Not zero" is the predicate; "positive" is a different one.
            if (! std::isfinite (e))
            {
                const bool excluded = core::exactlyEqual (w[c], 0.0);
                if (! excluded) poisoned = true;
                continue;
            }
            subMS += w[c] * (e / (double) subSamples);
        }
        if (poisoned && nonFiniteSubHops_ != ~std::uint64_t {}) ++nonFiniteSubHops_;
        if (! std::isfinite (subMS)) { subMS = 0.0; if (! poisoned && nonFiniteSubHops_ != ~std::uint64_t {}) ++nonFiniteSubHops_; }
        subRing[(std::size_t) subWrite] = subMS;
        subWrite = (subWrite + 1) % kSubRing;
        if (subFilled < kSubRing) ++subFilled;
        if (++subInHop >= kSubHopsPerHop) { subInHop = 0; finishHop(); }
        // Clear EVERY channel's accumulator, not just c < nc: a channel that vanishes mid-hop (host drops
        // the channel count) must not park its partial energy and leak it into a later hop when it returns.
        for (int c = 0; c < kMaxChannels; ++c) subSumSq[c] = 0.0;
        subCount = 0;
    }

    // Every 100 ms: a 400 ms gating block for the integrated measure and, once 3 s are in, a short-term
    // sample for LRA — one PER HOP since K12, which is the 10 Hz Tech 3342 3.1 asks for.
    void finishHop() noexcept
    {
        if (subFilled >= kMomentarySubHops)                                         // a 400 ms block every 100 ms
        {
            if (blockCount < (int) blockE.size()) blockE[(std::size_t) blockCount++] = meanLastSubHops (kMomentarySubHops);
            else ++droppedBlocks_;                                                  // past the capacity: counted, not kept
        }
        // A 3 s SHORT-TERM SAMPLE EVERY HOP — 10 Hz, which is what EBU Tech 3342 requires and what this meter
        // did not do. §3.1, verbatim: "using a sliding analysis-window of length 3 seconds for integration …
        // A minimum block overlap of 2.9 s between consecutive analysis windows (i.e. >=10 Hz sampling of the
        // loudness level) is required". With a 3 s window a 2.9 s overlap IS a 100 ms step, so the two
        // phrasings are one requirement. The requirement entered in V3 (January 2016); this meter kept
        // libebur128's 1 Hz cadence, which is a 2 s overlap.
        //
        // "OUR HOP IS EXACTLY IT" ONLY WHERE 0.01·fs IS A WHOLE NUMBER. A sub-hop is lround (0.01·fs) samples,
        // so at 48 and 44.1 kHz a hop is 100 ms to the sample and the cadence is 10 Hz exactly. At 22050 Hz the
        // sub-hop rounds 220.5 up to 221, a hop is 2210 samples = 100.227 ms, and the cadence is 9.9774 Hz —
        // UNDER the standard's minimum, by 0.23 %. At 8050 Hz it is 9.9383 Hz. The window rounds with it: 3.0068 s
        // at 22050, 3.0186 s at 8050, against a specified 3 s. It cuts the other way too — 8049 Hz rounds 80.49
        // DOWN to 80, giving 10.0613 Hz and a window of 2.9817 s — fast enough, and too short. That is the
        // meter's whole sub-hop grid, not something K12 introduced, and moving it would move every number this
        // class produces at those rates — so it is named here rather than quietly fixed, and the compliance
        // claim above is made for rates whose hundredth is whole.
        //
        // WHAT IT CHANGES — and the answer that stood here was WRONG. It read "on programmes of 30 s and
        // longer, nothing", which was true of the three fixtures it was measured on and false as a claim.
        // A square envelope whose states last exactly the window's 3 s reads 20.0 LU at 1 Hz and 9.5 LU at
        // 10 Hz, and that gap does NOT close with length: 1 Hz answers a flat 20.0 at 12, 20, 30, 60 and 120 s,
        // while 10 Hz answers 9.5 at all of them but 20 s, which is 9.3 — the same figures at 44.1 kHz and at
        // 48 kHz. (The sentence here used to call all five "the same two numbers"; the 20 s row was in the
        // sweep that produced it, reading 9.3, and the summary rounded a measurement away.)
        //
        // The mechanism is the PERCENTILES, not the sample count. At 1 Hz such a programme offers the window
        // six phases; a sixth of them sit on a pure loud block and a sixth on a pure quiet one, so P95 and P10
        // both land on an extreme and the range comes out the full 20 dB the envelope swings. At 10 Hz there
        // are sixty phases, each aligned one is 1.7 % — under the 5 % P95 reaches for — and the fifty-odd
        // mixed windows that 1 Hz never looked at fill the distribution in between. So the move is away from
        // a number that depended on where the grid happened to land, which is what the overlap is for: the
        // standard asks for it "to prevent loss of precision in the measurement of shorter programmes".
        // ON THIS FIXTURE the two cadences agree exactly away from the window's own timescale — 0 LU at 1.5 s
        // states, 20 LU at 5 s and at 7.5 s — and that is a property of the FIXTURE, not a regime. Its two
        // states are 20 dB apart, so both clear the -20 LU relative gate and the distribution keeps them both.
        // Drop the quiet state to -30 dB and it does not: at 60 s, 5 s states read 4.7 LU at 1 Hz against 7.7 at
        // 10 Hz, 7.5 s states 4.7 against 6.9, 10 s states 4.7 against 5.7 — slow envelopes, disagreeing. Even
        // "far faster" is not flat: 1 s states read 3.0 against 2.6 on both shapes. The 1 Hz figures throughout
        // were read off the pre-change meter and are not pinned by anything; the 10 Hz ones are —
        // LoudnessConformanceTests measures them.
        if (subFilled >= kSubRing)
        {
            if (stCount < (int) stE.size()) stE[(std::size_t) stCount++] = meanLastSubHops (kSubRing);
            else ++droppedShortTerm_;                                               // past the capacity: counted, not kept
        }
    }

    double meanLastSubHops (int k) const noexcept
    {
        const int kk = k < subFilled ? k : subFilled;
        if (kk <= 0) return 0.0;
        double s = 0.0;
        for (int j = 0; j < kk; ++j) { const int idx = (subWrite - 1 - j + kSubRing) % kSubRing; s += subRing[(std::size_t) idx]; }
        return s / kk;
    }

    static double lufsOf (double meanSquare) noexcept { return meanSquare > 1e-12 ? -0.691 + 10.0 * Math::log10 (meanSquare) : -120.0; }

    double integrated() const noexcept
    {
        if (blockCount <= 0) return -120.0;
        const double absT = Math::pow10 ((-70.0 + 0.691) / 10.0);               // energy for -70 LUFS
        double sum = 0.0; int cnt = 0;
        // isfinite FIRST, at every gate. `NaN > absT` is already false, but `+inf > absT` is TRUE, and one
        // +inf energy then poisons `sum`, makes relT infinite, and the second gate admits nothing — the
        // −120 "reads as silence" lie. Healing the filter makes this MORE reachable, not less: with the
        // state repaired at the next sub-hop the surrounding energies are finite, so a single infinity can
        // sit alone in a window instead of being swamped by NaNs.
        for (int j = 0; j < blockCount; ++j) { const double z = blockE[(std::size_t) j]; if (std::isfinite (z) && z > absT) { sum += z; ++cnt; } }
        if (cnt == 0) return -120.0;
        const double relT = 0.1 * (sum / cnt);                                     // -10 LU relative to the abs-gated mean
        double s2 = 0.0; int c2 = 0;
        for (int j = 0; j < blockCount; ++j) { const double z = blockE[(std::size_t) j]; if (std::isfinite (z) && z > absT && z > relT) { s2 += z; ++c2; } }
        return c2 > 0 ? lufsOf (s2 / c2) : -120.0;
    }

    // LRA (EBU Tech 3342) = P95 − P10 of the gated 3 s short-term loudness distribution. Two-pass gate over
    // the 10 Hz-cadence stE[] energies: absolute −70 LUFS, then −20 LU below the energy-mean of the abs-gated
    // set; percentiles via a fixed 0.1 LU histogram (−70..+30 LUFS), libebur128-faithful (no sort, no alloc).
    double lra() const noexcept
    {
        if (stCount <= 0) return 0.0;
        const double absT = Math::pow10 ((-70.0 + 0.691) / 10.0);                // energy for −70 LUFS
        double sum = 0.0; int cnt = 0;
        // isfinite first here too, and here it is not merely a wrong answer but UNDEFINED BEHAVIOUR: an
        // infinite short-term energy passes `inf >= absT`, makes relT infinite, then passes `inf >= inf`,
        // and the histogram index below computes (int) inf. Reproduced on main with UBSan — "inf is outside
        // the range of representable values of type 'int'" — from ONE +inf input sample, aligned so its
        // sub-hop is the last in a 3 s window. The clamp two lines down is too late: the conversion is the UB.
        for (int j = 0; j < stCount; ++j) { const double e = stE[(std::size_t) j]; if (std::isfinite (e) && e >= absT) { sum += e; ++cnt; } }
        if (cnt == 0) return 0.0;
        const double relT = 0.01 * (sum / cnt);                                     // −20 LU relative (energy ×0.01)

        constexpr int kBins = 1000;                                                 // −70..+30 LUFS, 0.1 LU bins
        int hist[kBins] = { 0 }; int total = 0;
        for (int j = 0; j < stCount; ++j)
        {
            const double e = stE[(std::size_t) j];
            if (std::isfinite (e) && e >= absT && e >= relT)                        // ≥ gates (libebur128-faithful)
            {
                int b = (int) ((lufsOf (e) + 70.0) * 10.0);                         // 0.1 LU bins from −70 LUFS
                b = b < 0 ? 0 : (b >= kBins ? kBins - 1 : b);
                ++hist[b]; ++total;
            }
        }
        if (total < 2) return 0.0;
        auto pct = [&] (double p) {                                                 // libebur128 rank rule, bin lower bound
            const int rank = (int) ((double) (total - 1) * p + 0.5);
            int cum = 0, b = 0;
            for (; b < kBins; ++b) { cum += hist[b]; if (cum > rank) break; }
            return -70.0 + (double) (b < kBins ? b : kBins - 1) * 0.1;
        };
        return pct (0.95) - pct (0.10);                                             // the 0.05-LU bin offset cancels
    }

    double fs = 48000.0; int ch = 2, subSamples = 480;
    int ranNc_ = 0;                             // channels that advanced K-weighting on the previous call
    bool prepared_ = false;                     // true only after prepare() (subRing/blockE/stE allocated)
    BasicKWeightingFilter<Math> kw;
    double w[kMaxChannels] {};
    double subSumSq[kMaxChannels] {};
    int subCount = 0;
    std::vector<double> subRing;                                                    // 300 × 10 ms sub-hop energies
    int subWrite = 0, subFilled = 0, subInHop = 0;
    std::vector<double> blockE;
    int blockCount = 0, droppedBlocks_ = 0, droppedShortTerm_ = 0;
    std::uint64_t nonFiniteSubHops_ = 0;   // sticky until reset()/prepare() — see the accessor
    std::vector<double> stE;                                                        // 3 s short-term energies @10 Hz (LRA)
    int stCount = 0;
};


using LoudnessMeter = BasicLoudnessMeter<core::SystemMath>;
using DeterministicLoudnessMeter = BasicLoudnessMeter<core::DetMath>;

} // namespace felitronics::analysis

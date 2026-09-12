// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/core/Config.h>
#include <felitronics/core/Math.h>

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <vector>

namespace felitronics::core
{

//==============================================================================
// felitronics::core::DeliveryResampler — the RATIONAL L:M sample-rate converter for DELIVERY: the
// first stage of a render whose file leaves at 44.1 / 48 / 88.2 / 96 / 176.4 / 192 kHz, where the
// conversion has to disappear under everything that follows it.
//
// THE BAR, and it is the reason this class exists rather than a call into a neighbour: stopband
// >= 140 dB, passband ripple <= 0.001 dB to 20 kHz. For scale, TPDF-dithered N-bit quantisation
// leaves 1/2^N rms re digital full scale — -96.3 dBFS at 16 bits, -144.5 at 24 — as a BROADBAND
// floor; an alias is a LINE, so the two are not directly comparable, and this file does not pretend
// they are.
//
// ------------------------------------------------------------------------------------------------
// THE FAMILY, now three. The other two are NOT interchangeable with this one and neither was widened
// to cover it:
//   * `core::StreamResampler` — the live rate-match feeding a NAM model locked to 48 kHz. It
//     interpolates LINEARLY between two rows of a 512-row phase table, and that interpolation floors
//     its error near -105 dB NO MATTER HOW MANY TAPS the kernel has: the floor is the phase grid, not
//     the kernel. A rational converter has no phase grid to interpolate — at L:M there are exactly L
//     phases and every one of them is exact.
//   * `convolution::resampleIr` — the offline one-shot that rate-converts an impulse response on load.
//     Same kernel family (Kaiser windowed sinc, designed in double), different budget: it evaluates the
//     kernel per output sample for a few thousand samples once, the wrong shape for a whole render.
//
// ------------------------------------------------------------------------------------------------
// THE DECOMPOSITION, which is the design and not an implementation detail.
//
// A single-stage prototype for 192 -> 44.1 kHz (147:640) is designed at the intermediate rate
// L*fs = 28.224 MHz, where the 20 kHz -> 22.05 kHz transition band is 7.3e-5 of the sample rate:
// 140,827 taps, 141,120 padded coefficients (1,102 KB) and 42.34 MMAC/s. The route this file picks is
// 192 -> 88.2 -> 44.1 — 12,792 coefficients (100 KB) and 26.99 MMAC/s, i.e. 11.0x less memory and
// 1.57x less work. It is NOT the textbook 192 -> 96 -> 48 -> 44.1 (29.96 MMAC/s): going to 88.2 first
// does the RATIONAL hop where its transition band is wide (20 -> 44.1 kHz, 84 taps per phase) and
// leaves a single 2:1 stage to pay for the narrow 20 -> 22.05 kHz one. Nobody wrote that route down;
// it is what costing the candidates returns.
//
// THE COST LAW the routes are chosen by. A stage's MAC rate is
//
//     K * f_in * f_out / df ,        K = (A - 8) / (2.285 * 2*pi) ,   df = its transition width
//
// with NO L in it: the polyphase prototype grows with L but only N/L of it is touched per output, so L
// is a MEMORY term and never a compute term. The lever a cascade has is putting the narrow transition
// where f_in*f_out is smallest — and on 176.4 <-> 192 there is no lever at all: both rates are at the
// top, the direct 160:147 is already optimal, and "put the rational stage at the lowest rate" applied
// literally would route 176.4 -> 88.2 -> 44.1 -> 48 -> 96 -> 192, a 22.05 kHz brickwall on a 192 kHz
// file. So the route is not a rule: `plan()` costs every candidate and keeps the cheapest.
//
// 🔴 THE CANDIDATES ARE THE ENDPOINTS' POWER-OF-TWO RELATIVES, AND THAT IS A CHOICE WITH A MEASURED
// PRICE. It is not "the whole useful space" — an earlier version of this comment said so and was
// wrong. Searching every rate a*p/q and b*p/q with p, q <= 32 finds cheaper routes on 20 of the 30
// pairs, by 3.4 % to 10.3 %, 4.6 % over the whole matrix — and it never improves the WORST pair, which
// stays at 43.04 MMAC/s, because that pair is direct either way. What it buys that with is
// intermediates like 56.448 and 63.7 kHz, and memory: 192 -> 44.1 through 63.7 kHz is 24.81 MMAC/s
// against 26.99 for 8.6x the coefficients (109,896 against 12,792), and the plain 1:2 pairs pay up to
// 20x (444 -> 8,892) to shave 3.7 %. Not taken; the render-time bound is the worst pair's, and the
// worst pair does not move.
//
// ------------------------------------------------------------------------------------------------
// THE PASSBAND EDGE IS PROPORTIONAL, and this is a product decision, not a numeric one.
//
//     passband edge = passbandFraction * min(inRate, outRate) / 2 ,   default 20000/22050 = 0.90703
//
// so it lands on exactly 20.000 kHz whenever 44.1 kHz is one of the two rates — which is where the
// stated bar applies — and rises with the rates otherwise: 40 kHz for 88.2 <-> 96, 80 kHz for
// 176.4 <-> 192. The fixed 20 kHz of the specification is the band the ACCEPTANCE is measured in, not
// the edge of the filter. Reading it as the edge is cheap and wrong: 176.4 -> 192 with a 20 kHz
// brickwall costs 4,481 taps instead of 35,201 precisely BECAUSE it throws away the
// 20-80 kHz a 192 kHz deliverable is supposed to carry. It is the one number here that changes the
// output's CONTENT rather than its cost, so it is a parameter with a stated default, not a constant.
//
// ------------------------------------------------------------------------------------------------
// 🔴 THE DESIGN TARGET IS SET BY MEASUREMENT, NOT BY THE FORMULA. Kaiser's attenuation fit is stated
// for A <= 120 and under-delivers past it: asking for 140 dB and measuring the prototype back gives
// -139.30 dB here and -138.81 dB through the tree's own `designFilter()`. The formula sizes; it does
// not certify. So the target was chosen by sweeping the whole 30-pair matrix through the converter —
// worst spur of the matrix, and the work of its worst pair:
//
//     edge margin   A = 143     A = 146     A = 148      worst pair at A = 146
//        0          -140.41     -142.01     -139.00 !      40.70 MMAC/s
//        0.03       -138.40 !   -141.44     -142.26        41.63
//        0.06       -142.32     -146.58     -148.43        43.04      <- shipped
//        0.10       -143.13     -145.72     -147.42        45.16
//
// WITHOUT AN EDGE MARGIN THE MATRIX IS NOT MONOTONE IN A — asking for 2 dB MORE attenuation took it
// from -142.0 to -139.0, through the bar. The worst line in that row always sits EXACTLY on a Nyquist:
// an input tone at the output Nyquist, or the literal (+1, -1) input Nyquist on an interpolating pair.
// With the transition ending exactly there, that point lives on the end of the transition SLOPE, and
// which lobe the rounded tap count puts on it decides +-10 dB. So the design stop edge is moved 6 % of
// the transition width inside the Nyquist (`kStopEdgeMargin`); the passband edge does not move. Then
// the matrix is monotone and the Nyquist sits in developed stopband: 6.6 dB under the bar for +5.7 %
// work on the worst pair. (A first version of this block quoted -144.46 dB for A = 146 with no margin.
// That was a broken instrument — rectangular window, sparse stopband probes, no probe at the Nyquist —
// and the true figure for that design was -142.01; see DeliveryResamplerBarTests.cpp.)
// `DeliveryResamplerBarTests.cpp` re-measures the shipped coefficients and gates the DELIVERED figure,
// not only the -140 dB requirement — a gate at the requirement alone stayed green when a mutation
// stand dropped the target to 140, which is the regression this default exists to prevent.
//
// 🔴 THE RIPPLE BAR IS NOT THE BINDING ONE. A Kaiser's passband deviation is about the size of its
// stopband ripple, so at this attenuation it is ~1e-6 dB — three orders of magnitude under the
// 0.001 dB asked for. Ripple comes free with the stopband and needs no cascade budget. Worst measured
// over all 30 pairs: 0.0000013 dB; worst absolute gain error 0.0000014 dB.
//
// ------------------------------------------------------------------------------------------------
// 🔴 DOUBLE IN THE SAMPLE LOOP — a law-3 carve-out, and float32 does not merely cost margin here, it
// FAILS THE BAR. Law 3 (DSP-ARCHITECTURE.md §2) is "float in the hot path, double only in offline
// coefficient design", with per-module carve-outs named at the point of use. This is one. Change ONLY
// the four accumulators in `dot()` to float32 — coefficients, history and products still double — and
// the bar suite goes red on 13 of the 30 pairs, worst -132.19 dB on 44.1 -> 88.2. On a periodic input
// the rounding error is a deterministic function of sample and phase, so it lands as spectral LINES
// rather than a floor, and the dyadic stages are hit hardest: with one or two phases the pattern
// repeats every sample or two and concentrates into very few, very high lines. Storing the
// coefficients or the history in float32 instead costs the worst pair a few tenths of a dB; it is the
// accumulation that cannot be narrowed.
//
// So `core::firDot` is deliberately NOT reused. That is not a criticism of it: its float32 arithmetic
// and its bit-identical-across-five-rows contract are the same decision, and P57 explicitly does not
// ask for bit agreement between rows — the bar is a tolerance, not a bit pattern. The design also
// goes through libm's `sin` (the Bessel series and `sqrt` are exact-rounding-safe), whose last places
// differ between rows; that is the other half of why bits are not promised.
//
// ------------------------------------------------------------------------------------------------
// LATENCY IS EXACT, BY CONSTRUCTION AND NOT BY MEASUREMENT.
//
// Every stage is built with N = 2*L*k + 1 taps, symmetric about index L*k, so output n of that stage
// is centred on input position n*M/L - k and the stage's delay is EXACTLY k of ITS OWN input samples —
// an integer, with no half-grid-sample residue. (An even N puts the centre between two design-grid
// points; a one-way converter has no return leg to cancel that against. The +1 forbids it.) The
// cascade's total is the sum of the stages' delays carried to the output rate, and it is a RATIONAL
// number: 48 -> 44.1 -> 48 is 120 + 110*(48000/44100) = 239.727891 samples at 48 kHz, and 192 -> 44.1
// is 41*(44100/192000) + 220*(44100/88200) = 119.4171875. `latencyOutputSamples()` returns it exactly;
// do NOT round it and then null against the result — this file's own first probe did, and read -10 dB.
// The suite measures the two-stage figure back out of the carrier phase, because recomputing the
// formula from the plan's own fields is a tautology that cannot see a mis-composed cascade.
//
// ------------------------------------------------------------------------------------------------
// WHAT THE ROUND TRIP MEASURES, AND WHAT IT CANNOT. 48 -> 44.1 -> 48 against an analytic oracle nulls
// at -143.3 dBFS peak / -160.3 dBFS RMS, well past the -120 dB asked of it. What sets that floor is
// deliberately NOT claimed here. Three float32 quantisations sit in the chain (float in, float between
// the two converters, float out) and the cascade's passband deviation is of the same order, and every
// one of those terms scales with the signal — so how the residual moves with level cannot say which
// dominates. Two earlier versions of this paragraph each read the scaling as proof of a mechanism
// ("passband ripple", then "three float32 binades"), and the next design change broke both readings.
// What the round trip CAN do is bound the time-domain error; it cannot see a -140 dB stopband, which is
// why the bar is gated on the spectrum instead.
//
// ------------------------------------------------------------------------------------------------
// WHY THE HOUSE `process (io, nch, n)` SIGNATURE IS NOT USED. A rational converter emits a VARIABLE
// number of samples per call — it depends on the phase state, compounds through the cascade's
// per-stage rounding, and is not in general `n` — so the shape is push/pull, like `StreamResampler`'s,
// plus the thing that one does not have and a file render cannot do without: `flush()`, which drains
// the group delay. Without it every rendered file is short by the cascade's delay — 110.25 samples at
// 48 -> 44.1, 520.3 at 44.1 -> 192, the longest in the matrix (2.71 ms).
//
// LAW 11, clause by clause, because the first version of this class claimed it and did not honour it:
//   * (a) THE LENGTH IS A CAPACITY. `maxInputBlock` sizes scratch; `process()` CHUNKS a longer call
//     and processes it in full. A short OUTPUT buffer is still refused — that is a request that cannot
//     be honoured, not a long one.
//   * (b) THE WIDTH IS A LIMIT. `nch > maxChannels` is refused before anything moves.
//   * 11a, THE FALLING EDGE. A channel whose index is >= `nch` on a call with n > 0 has STOPPED; its
//     history is dropped, so it cannot replay audio when it returns (it did: a two-channel render that
//     went mono for three blocks brought the stopped channel back at +0.4 dBFS out of digital
//     silence). `nch == 0` with n > 0 stops every channel AND spends the audio time on the clock, so a
//     gap is a gap and not a pause — resuming after it is bit-identical to having fed the gap as
//     silence, once the silence is longer than the filter.
//   * `flush()` ends the programme. `process()` is refused after it until `reset()`: accepting it
//     used to append audio whose tail a second, idempotent `flush()` then silently withheld.
//
// 🔴 NON-INTEGER RATES ARE REFUSED, not approximated, and so is any plan past a hard size. The ratio is
// reduced with an exact gcd on integer rates; a reduced L or M past `kMaxRatioTerm`, or a stage past
// `kMaxStageCoefficients`, is refused BEFORE a single narrowing conversion. That bound is not
// decoration: `passbandFraction = nextafter(1.0, 0.0)` is finite, positive and under 1, and it used to
// produce an out-of-range double->int conversion and a signed overflow inside `plan()` and a plan that
// reported `ok` with zero coefficients. A finiteness check is not a range check.
//==============================================================================
struct DeliveryResampler
{
    static constexpr int       kMaxStages            = 3;          // the search's longest route
    static constexpr int       kMaxRatioTerm         = 8192;       // reduced L or M per hop
    static constexpr long long kMaxStageCoefficients = 1LL << 22;  // taps AND padded storage, per stage
    static constexpr double    kDesignStopbandDb     = 146.0;      // see "set by measurement" above
    static constexpr double    kStopEdgeMargin       = 0.06;       // of the transition width; see the same block
    static constexpr double    kDefaultPassbandFraction = 20000.0 / 22050.0;   // exactly 20 kHz at 44.1

    struct Params
    {
        double inRate  = 48000.0;
        double outRate = 48000.0;
        // Passband edge as a fraction of the LOWER Nyquist. The default puts it on 20.000 kHz whenever
        // 44.1 kHz is involved; see the header block for why it is not a fixed 20 kHz.
        double passbandFraction = kDefaultPassbandFraction;
        // Kaiser sizing target in dB — not the delivered attenuation; see the header block.
        double stopbandDb = kDesignStopbandDb;
    };

    //==========================================================================
    // THE PLAN. A pure function of the parameters: no allocation, no state, nothing to prepare, so the
    // suite walks every delivery pair through it without constructing a resampler.
    struct StagePlan
    {
        long long inRate = 0, outRate = 0;   // this stage's own rates
        int L = 1, M = 1;                    // reduced ratio, outRate:inRate
        int halfLen = 0;                     // k — the stage's delay in ITS OWN input samples
        int tapsPerPhase = 0;                // 2k+1, padded to a multiple of 4
        double passbandHz = 0.0;             // flat to here
        double stopbandHz = 0.0;             // the DESIGN stop edge: inside the Nyquist by kStopEdgeMargin
        double nyquistHz = 0.0;              // min(inRate, outRate) / 2 — where aliasing and imaging begin
    };

    struct Plan
    {
        bool ok = false;
        bool identity = false;               // exactly equal rates: a copy, no filter, no delay
        int  count = 0;
        StagePlan stage[kMaxStages];
        double macsPerSecond = 0.0;          // the cost the route was chosen by
        long long coefficients = 0;          // padded doubles across every stage
    };

    // ---- the cost law -------------------------------------------------------
    static constexpr double kKaiserSlope = 2.285;

    // Half-length k for a Kaiser prototype, N ~ (A - 8) / (2.285 * dw), dw = 2*pi*df/Fp, returned in
    // the 2*L*k+1 form. REFUSES — rather than narrowing — anything past the stage bound: the double is
    // range-checked before it is ever converted, because a finite double is not a representable int.
    static bool planHalfLen (double stopbandDb, double transitionHz, double protoRateHz, int L,
                             int& kOut) noexcept
    {
        const double dw = 2.0 * kPi * transitionHz / protoRateHz;
        const double n  = (stopbandDb - 8.0) / (kKaiserSlope * dw);
        const double cap = (double) kMaxStageCoefficients;
        // `n <= cap` refuses early; the `k` test below and `makeStage`'s tap test would refuse the same
        // plans without it (a mutant removing it left 144,360 plans bit-identical, 288 of them refusals).
        // It stays because it is the one that runs before any further arithmetic on an absurd `n`.
        if (! std::isfinite (n) || ! (n > 0.0) || ! (n <= cap)) return false;
        const double k = std::max (1.0, std::ceil (n / (2.0 * (double) L)));
        if (! (k <= cap)) return false;
        kOut = (int) k;
        return true;
    }

    static double kaiserBeta (double stopbandDb) noexcept
    {
        if (stopbandDb > 50.0) return 0.1102 * (stopbandDb - 8.7);
        if (stopbandDb >= 21.0) return 0.5842 * std::pow (stopbandDb - 21.0, 0.4) + 0.07886 * (stopbandDb - 21.0);
        return 0.0;
    }

    // One stage of a route, or false when this hop cannot carry the wanted band or would exceed a bound.
    static bool makeStage (long long inRate, long long outRate, double passbandHz,
                           double stopbandDb, StagePlan& s) noexcept
    {
        const long long g = std::gcd (inRate, outRate);
        const long long L = outRate / g, M = inRate / g;
        if (L > kMaxRatioTerm || M > kMaxRatioTerm) return false;

        const double nyq = 0.5 * (double) std::min (inRate, outRate);
        // THE guard against a route that loses band — the planner's candidate filter defers to this.
        if (! (passbandHz < nyq)) return false;
        // The DESIGN stop edge sits a fixed fraction of the transition INSIDE the Nyquist, so that the
        // Nyquist itself — the first point of the stopband — lies in developed attenuation rather than on
        // the end of the transition slope. See "set by measurement" in the header block for what that
        // slope did to the matrix without it.
        const double stop = nyq - kStopEdgeMargin * (nyq - passbandHz);

        int k = 0;
        if (! planHalfLen (stopbandDb, stop - passbandHz, (double) L * (double) inRate, (int) L, k)) return false;
        const long long taps = 2LL * L * (long long) k + 1;
        const long long tpp  = ((2LL * (long long) k + 1) + 3) & ~3LL;
        if (taps > kMaxStageCoefficients || L * tpp > kMaxStageCoefficients) return false;

        s.inRate = inRate; s.outRate = outRate;
        s.L = (int) L; s.M = (int) M;
        s.passbandHz = passbandHz; s.stopbandHz = stop; s.nyquistHz = nyq;
        s.halfLen = k;
        s.tapsPerPhase = (int) tpp;
        return true;
    }

    static double stageMacsPerSecond (const StagePlan& s) noexcept
    {
        return (double) s.tapsPerPhase * (double) s.outRate;
    }

    static Plan plan (const Params& p) noexcept
    {
        Plan best;
        if (! (std::isfinite (p.inRate) && std::isfinite (p.outRate))) return best;
        if (! (p.inRate > 0.0) || ! (p.outRate > 0.0)) return best;
        if (! (std::isfinite (p.passbandFraction) && p.passbandFraction > 0.0 && p.passbandFraction < 1.0)) return best;
        if (! (std::isfinite (p.stopbandDb) && p.stopbandDb >= 20.0 && p.stopbandDb <= 300.0)) return best;

        // Integer rates only — the ratio has to be exact. The range test comes FIRST, so `llround`
        // never sees a value it cannot represent; then an equality test rejects 44100.5.
        if (p.inRate < 1000.0 || p.outRate < 1000.0 || p.inRate > 3.0e6 || p.outRate > 3.0e6) return best;
        const long long a = (long long) std::llround (p.inRate);
        const long long b = (long long) std::llround (p.outRate);
        if (! exactlyEqual ((double) a, p.inRate) || ! exactlyEqual ((double) b, p.outRate)) return best;

        if (a == b) { best.ok = true; best.identity = true; best.count = 0; return best; }

        const double bandHz = p.passbandFraction * 0.5 * (double) std::min (a, b);

        // Candidates: the endpoints' power-of-two relatives — a deliberate restriction with a measured
        // price, see the header block.
        long long cand[16];
        int nc = 0;
        for (long long base : { a, b })
            for (int e = -3; e <= 3; ++e)
            {
                const long long r = (e >= 0) ? (base << e) : (base >> (-e));
                if (r == a || r == b) continue;
                if ((e < 0) && ((r << (-e)) != base)) continue;             // not divisible — skip
                // A PRUNING, not the guard. `makeStage` refuses any hop that cannot carry the band, and
                // it is the only thing that has to: removing this line was run as a mutant and left all
                // 624 plans of a 13-rate x 4-passband sweep bit-identical. It only spares the search
                // from costing routes that were never going to be legal.
                if (0.5 * (double) r <= bandHz) continue;
                bool dup = false;
                for (int i = 0; i < nc; ++i) if (cand[i] == r) dup = true;
                if (! dup && nc < 16) cand[nc++] = r;
            }

        auto consider = [&] (const long long* path, int hops) noexcept
        {
            Plan cur;
            cur.count = hops;
            double macs = 0.0;
            long long coefs = 0;
            for (int i = 0; i < hops; ++i)
            {
                if (! makeStage (path[i], path[i + 1], bandHz, p.stopbandDb, cur.stage[i])) return;
                macs  += stageMacsPerSecond (cur.stage[i]);
                coefs += (long long) cur.stage[i].L * (long long) cur.stage[i].tapsPerPhase;
            }
            cur.ok = true;
            cur.macsPerSecond = macs;
            cur.coefficients  = coefs;
            if (! best.ok || macs < best.macsPerSecond) best = cur;
        };

        { const long long path[2] = { a, b }; consider (path, 1); }
        for (int i = 0; i < nc; ++i)
        {
            const long long path[3] = { a, cand[i], b };
            consider (path, 2);
            for (int j = 0; j < nc; ++j)
            {
                if (j == i) continue;
                const long long p3[4] = { a, cand[i], cand[j], b };
                consider (p3, 3);
            }
        }
        return best;
    }

    //==========================================================================
    // Total delay in OUTPUT samples — exact, and generally not an integer. Each stage delays by
    // `halfLen` of its own input samples; carried to the final output rate that is
    // halfLen * outRate_final / inRate_stage.
    static double latencyOutputSamples (const Plan& pl) noexcept
    {
        if (! pl.ok || pl.count == 0) return 0.0;
        const double fsOut = (double) pl.stage[pl.count - 1].outRate;
        double d = 0.0;
        for (int i = 0; i < pl.count; ++i)
            d += (double) pl.stage[i].halfLen * fsOut / (double) pl.stage[i].inRate;
        return d;
    }

    //==========================================================================
    // WHAT prepare() ASKS THE HEAP FOR (law 11d), and the ONE function it sizes itself with — so an owner
    // budgeting before it builds reads the numbers the buffers are actually built from, and a clamp or a
    // refusal cannot live in one of the two and not the other. `ok == false` is exactly the set of
    // arguments prepare() refuses, and then it asks for nothing.
    //
    // REQUESTED bytes, counted per `operator new`: the coefficient tables, the histories, the stage and
    // input scratch, the transient full-length prototype `designPhases` builds and frees once per stage,
    // and the arrays of vector headers that hold them. Not a promise that a heap can serve them —
    // allocator headers are not counted — and a SUM, not a peak: each prototype is freed before the next
    // stage's is made, so the peak is lower by all but the largest of them.
    struct Storage
    {
        bool ok = false;
        Plan plan {};
        int flushInputs = 0, sizeBlock = 0;
        int stageBound[kMaxStages] {};
        std::uint64_t doubles = 0, ints = 0, vectorHeaders = 0;
        std::uint64_t bytes() const noexcept
        {
            return ok ? doubles * sizeof (double) + ints * sizeof (int)
                          + vectorHeaders * sizeof (std::vector<double>)
                      : 0u;
        }
    };

    static Storage storageFor (const Params& p, int maxChannels, int maxInputBlock) noexcept
    {
        Storage st;
        if (maxChannels < 1 || maxChannels > kMaxChannels) return st;
        if (maxInputBlock < 1 || maxInputBlock > (1 << 22)) return st;
        const Plan pl = plan (p);
        if (! pl.ok) return st;

        // flush() pushes the whole group delay in as zeros, counted in stage-0 input samples. It is
        // computed BEFORE the scratch is sized, because the scratch has to hold a flush as well as a
        // full block: sizing from `maxInputBlock` alone refused the flush of any caller whose blocks are
        // smaller than the delay, which is most of them.
        double dIn = 0.0;
        if (pl.count > 0)
        {
            const double fsIn = (double) pl.stage[0].inRate;
            for (int i = 0; i < pl.count; ++i)
                dIn += (double) pl.stage[i].halfLen * fsIn / (double) pl.stage[i].inRate;
        }
        st.flushInputs = pl.identity ? 0 : ((int) std::ceil (dIn) + 1);
        st.sizeBlock   = std::max (maxInputBlock, st.flushInputs);

        const std::uint64_t ch = (std::uint64_t) maxChannels;
        long long bound = st.sizeBlock;
        for (int s = 0; s < pl.count; ++s)
        {
            const StagePlan& sp = pl.stage[s];
            bound = (bound * (long long) sp.L) / (long long) sp.M + 1;     // worst case over the phase state
            if (bound > (1 << 24)) return st;                               // absurd; refuse, do not allocate it
            st.stageBound[s] = (int) bound;
            st.doubles += (std::uint64_t) (2LL * sp.L * sp.halfLen + 1)     // the transient prototype
                        + (std::uint64_t) sp.L * (std::uint64_t) sp.tapsPerPhase   // its phase-major table
                        + ch * (std::uint64_t) (2 * sp.tapsPerPhase)        // histories
                        + ch * (std::uint64_t) bound;                       // stage scratch
        }
        st.doubles += ch * (std::uint64_t) st.sizeBlock;                    // input scratch
        st.ints = (std::uint64_t) pl.count;                                 // stageBound
        // The arrays of vector headers: the coefficient tables' (built aside, then swapped in), the
        // histories', the stage scratch's and the input scratch's. For a FRESH object — a re-prepare can
        // reuse header arrays of the same length, and then asks for less than this.
        st.vectorHeaders = (std::uint64_t) pl.count + 2u * (std::uint64_t) pl.count * ch + ch;
        st.plan = pl;
        st.ok = true;
        return st;
    }

    [[nodiscard]] static std::uint64_t prepareBytes (const Params& p, int maxChannels, int maxInputBlock) noexcept
    {
        return storageFor (p, maxChannels, maxInputBlock).bytes();
    }

    //==========================================================================
    [[nodiscard]] bool prepare (const Params& p, int maxChannels, int maxInputBlock)
    {
        // Law 11b: DISARM, VALIDATE, WRITE. Nothing below is adopted until every argument is honoured —
        // and the validation IS the budget: storageFor() refuses exactly what this refuses.
        prepared = false;
        const Storage st = storageFor (p, maxChannels, maxInputBlock);
        if (! st.ok) return false;
        const Plan& pl = st.plan;

        std::vector<std::vector<double>> newCoef ((std::size_t) pl.count);
        for (int s = 0; s < pl.count; ++s)
            newCoef[(std::size_t) s] = designPhases (pl.stage[s], p.stopbandDb);

        // ---- adopt
        thePlan = pl;
        channels = maxChannels;
        maxBlock = maxInputBlock;
        flushInputs = st.flushInputs;
        coef.swap (newCoef);
        stageBound.assign (st.stageBound, st.stageBound + pl.count);

        hist.assign ((std::size_t) (pl.count * maxChannels), {});
        for (int s = 0; s < pl.count; ++s)
            for (int c = 0; c < maxChannels; ++c)
                hist[(std::size_t) (s * maxChannels + c)].assign ((std::size_t) (2 * pl.stage[s].tapsPerPhase), 0.0);

        scratch.assign ((std::size_t) (pl.count * maxChannels), {});
        for (int s = 0; s < pl.count; ++s)
            for (int c = 0; c < maxChannels; ++c)
                scratch[(std::size_t) (s * maxChannels + c)].assign ((std::size_t) st.stageBound[s], 0.0);

        inScratch.assign ((std::size_t) maxChannels, {});
        for (int c = 0; c < maxChannels; ++c) inScratch[(std::size_t) c].assign ((std::size_t) st.sizeBlock, 0.0);

        prepared = true;
        reset();
        return true;
    }

    void reset() noexcept
    {
        for (auto& h : hist) std::fill (h.begin(), h.end(), 0.0);
        for (int s = 0; s < kMaxStages; ++s)
        {
            phaseQ[s]    = 0;
            ringPos[s]   = 0;
            pendingIn[s] = 1;        // one push is owed before output 0, whose base index is 0
        }
        activeWidth = 0;             // every ring is silent, so no channel has anything to drop
        flushed = false;
    }

    // The most outputs a call of `n` inputs can write — independent of how the call is chunked, because
    // which outputs exist depends only on how many inputs have been fed. SATURATES at INT_MAX; a call
    // whose bound does not fit in `nOut` cannot be honoured and is refused.
    [[nodiscard]] int maxOutputFor (int n) const noexcept
    {
        const long long b = outputBound (n);
        return b > (long long) INT_MAX ? INT_MAX : (int) b;
    }

    // What `flush()` needs. A caller cannot derive it — the drain length is the cascade's own group
    // delay, not the caller's block — and getting it wrong is a refused flush, i.e. a truncated file.
    [[nodiscard]] int maxFlushOutput() const noexcept { return prepared ? maxOutputFor (flushInputs) : 0; }
    [[nodiscard]] int flushInputSamples() const noexcept { return flushInputs; }

    double latencyOutputSamples() const noexcept { return latencyOutputSamples (thePlan); }
    const Plan& currentPlan() const noexcept { return thePlan; }
    bool isPrepared() const noexcept { return prepared; }

    //==========================================================================
    // Law 11's order: malformed -> unprepared -> ended -> too wide -> n == 0 -> nch == 0 -> run.
    [[nodiscard]] bool process (const float* const* in, int nch, int n,
                                float* const* out, int outCapacity, int& nOut) noexcept
    {
        nOut = 0;
        if (n < 0 || nch < 0 || outCapacity < 0) return false;
        if (! prepared) return false;
        if (flushed) return false;                                   // the programme ended; reset() first
        if (nch > channels) return false;
        if (n == 0) return true;
        if (nch == 0)
        {
            // 11a: every channel stopped, and n samples of audio time passed. Planes may be null here.
            stopChannelsFrom (0);
            advanceClock (n);
            return true;
        }
        if (in == nullptr || out == nullptr) return false;
        for (int c = 0; c < nch; ++c)
            if (in[c] == nullptr || out[c] == nullptr) return false;
        if ((long long) outCapacity < outputBound (n)) return false;

        stopChannelsFrom (nch);

        if (thePlan.identity)
        {
            // Exactly equal rates. A 0.907-Nyquist low-pass is a real filter, and a caller asking for
            // no rate change must not silently get one — so this path is a copy, bit for bit.
            for (int c = 0; c < nch; ++c) std::copy (in[c], in[c] + n, out[c]);
            nOut = n;
            return true;
        }

        // 11(a): the length is a capacity. Chunks of `maxBlock`; the outputs land contiguously.
        int written = 0;
        for (int off = 0; off < n; off += maxBlock)
        {
            const int chunk = std::min (maxBlock, n - off);
            for (int c = 0; c < nch; ++c)
            {
                double* dst = inScratch[(std::size_t) c].data();
                for (int i = 0; i < chunk; ++i) dst[i] = (double) in[c][off + i];
            }
            written += runChain (nch, chunk, out, written);
        }
        nOut = written;
        return true;
    }

    // Drains the group delay and ENDS THE PROGRAMME: `process()` is refused afterwards until `reset()`.
    // A second `flush()` is accepted and writes nothing. A render that skips it loses its tail.
    [[nodiscard]] bool flush (int nch, float* const* out, int outCapacity, int& nOut) noexcept
    {
        nOut = 0;
        if (nch < 0 || outCapacity < 0) return false;
        if (! prepared) return false;
        if (nch > channels) return false;
        if (flushed) return true;
        if (nch == 0) { stopChannelsFrom (0); flushed = true; return true; }
        if (out == nullptr) return false;
        for (int c = 0; c < nch; ++c)
            if (out[c] == nullptr) return false;
        if (thePlan.identity) { flushed = true; return true; }      // a copy has no tail to drain
        if ((long long) outCapacity < outputBound (flushInputs)) return false;

        stopChannelsFrom (nch);
        for (int c = 0; c < nch; ++c)
            std::fill_n (inScratch[(std::size_t) c].data(), flushInputs, 0.0);
        nOut = runChain (nch, flushInputs, out, 0);
        flushed = true;
        return true;
    }

private:
    //==========================================================================
    // The prototype, phase-major. Designed and stored entirely in double — see the law-3 carve-out in
    // the header block. N = 2*L*k + 1 taps, symmetric about index L*k, normalised GLOBALLY to a DC gain
    // of L — which is NOT "every phase sums to 1": the phase sums spread by the stopband ripple
    // (+-4e-9 on the default 48 -> 44.1 stage), and that spread is an LPTV image of DC at that level.
    //
    // 🔴 DO NOT "FIX" THAT WITH PER-PHASE NORMALISATION, the way `StreamResampler` does per row. It was
    // done here and measured: the worst spur of the matrix went from -146.58 dB to -143.74 dB, 2.8 dB
    // worse, through the delivered-figure gate. Dividing phase q by its own sum (1 + e_q) turns the image of a
    // tone at f from a(f) — the stopband response at that image — into a(f) - e_q, i.e. it subtracts
    // the DC image from EVERY frequency: exact at DC, and up to double the error wherever a(f) sits in
    // a null. `StreamResampler` needs it because a 64-tap kernel's phase sums drift by 6 ppm, far above
    // its own stopband; here they already sit inside a -146 dB stopband, so there is nothing to gain
    // and a lobe to lose.
    static std::vector<double> designPhases (const StagePlan& s, double stopbandDb)
    {
        const long long N = 2LL * s.L * (long long) s.halfLen + 1;
        const double beta = kaiserBeta (stopbandDb);
        const double i0b  = besselI0 (beta);
        const double mid  = (double) ((long long) s.L * (long long) s.halfLen);   // (N-1)/2, an integer
        const double fcn  = 0.5 * (s.passbandHz + s.stopbandHz) / ((double) s.L * (double) s.inRate);

        std::vector<double> h ((std::size_t) N);
        for (long long i = 0; i < N; ++i)
        {
            const double x = (double) i - mid;
            const double sinc = (std::fabs (x) < 1.0e-13) ? (2.0 * fcn)
                                                          : std::sin (2.0 * kPi * fcn * x) / (kPi * x);
            const double u = x / (mid + 1.0);
            const double w = besselI0 (beta * std::sqrt (std::max (0.0, 1.0 - u * u))) / i0b;
            h[(std::size_t) i] = sinc * w;
        }

        double dc = 0.0;
        for (double v : h) dc += v;
        const double g = exactlyEqual (dc, 0.0) ? 1.0 : ((double) s.L / dc);
        for (double& v : h) v *= g;

        // Phase q holds h[q], h[q+L], h[q+2L], … — the taps that multiply x[base], x[base-1], … in
        // that order, so the history window is read newest-first and contiguously.
        const std::size_t T = (std::size_t) s.tapsPerPhase;
        std::vector<double> out ((std::size_t) s.L * T, 0.0);
        for (int q = 0; q < s.L; ++q)
            for (std::size_t i = 0; i < T; ++i)
            {
                const long long idx = (long long) q + (long long) i * (long long) s.L;
                if (idx < N) out[(std::size_t) q * T + i] = h[(std::size_t) idx];
            }
        return out;
    }

    static double besselI0 (double x) noexcept
    {
        double sum = 1.0, term = 1.0;
        const double y = x * x * 0.25;
        for (int k = 1; k < 200; ++k)
        {
            term *= y / ((double) k * (double) k);
            sum  += term;
            if (term < 1.0e-18 * sum) break;
        }
        return sum;
    }

    // Four double accumulators — here for the shorter dependency chain, NOT for bit agreement across
    // rows, which this class does not promise. Narrowing THESE to float32 is what fails the bar.
    static double dot (const double* a, const double* b, int len) noexcept
    {
        double s0 = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0;
        for (int i = 0; i < len; i += 4)
        {
            s0 += a[i + 0] * b[i + 0];
            s1 += a[i + 1] * b[i + 1];
            s2 += a[i + 2] * b[i + 2];
            s3 += a[i + 3] * b[i + 3];
        }
        return (s0 + s1) + (s2 + s3);
    }

    // Push one sample into a channel's double-length newest-first ring.
    static void push (double* ring, int T, int& pos, double x) noexcept
    {
        pos = (pos == 0) ? (T - 1) : (pos - 1);
        ring[pos] = x;
        ring[pos + T] = x;
    }

    // Worst-case output count for `n` more inputs, independent of the phase state.
    long long outputBound (long long n) const noexcept
    {
        if (! prepared || n < 0) return 0;
        if (thePlan.identity) return n;
        long long b = n;
        for (int s = 0; s < thePlan.count; ++s)
            b = (b * (long long) thePlan.stage[s].L) / (long long) thePlan.stage[s].M + 1;
        return b;
    }

    // 11a: channels [w, activeWidth) have stopped — drop their memory. Channels that start again later
    // find silent rings, which is exactly "silence while they were away".
    void stopChannelsFrom (int w) noexcept
    {
        if (w < activeWidth)
            for (int s = 0; s < thePlan.count; ++s)
                for (int c = w; c < activeWidth; ++c)
                {
                    auto& h = hist[(std::size_t) (s * channels + c)];
                    std::fill (h.begin(), h.end(), 0.0);
                }
        activeWidth = w;
    }

    //--------------------------------------------------------------------------
    // 🔴 THE PHASE AUTOMATON. Output n reads phase q = (n*M) mod L, centred on input index
    // base = (n*M) div L, so from one output to the next `base` advances by (q + M) div L — 0, 1, or
    // more. `pend` counts that advance DOWN in pushes still owed: decremented once per input sample,
    // and every output it reaches zero for is emitted. The easy way to write this wrongly is "one
    // output per input sample", and it is wrong in BOTH directions: at L:M = 2:1 each input owes two
    // outputs and it emits one; at 1:2 each input owes half of one and it emits one — twice as many.
    // Zero-priming makes the first output's base 0, hence `pend = 1` out of reset().
    //
    // The state is identical for every active channel of a stage, so the channels run from one entry
    // state and the exit state is adopted once.
    struct StageState { int q, pend, pos; };

    int runStage (int s, int c, const double* src, int nIn, double* dst, StageState& st) noexcept
    {
        const StagePlan& sp = thePlan.stage[s];
        const int T = sp.tapsPerPhase;
        double* ring = hist[(std::size_t) (s * channels + c)].data();
        const double* cf = coef[(std::size_t) s].data();
        int pos = st.pos, q = st.q, pend = st.pend, k = 0;
        for (int i = 0; i < nIn; ++i)
        {
            push (ring, T, pos, src[i]);
            --pend;
            while (pend == 0)
            {
                dst[k++] = dot (cf + (std::size_t) q * (std::size_t) T, ring + pos, T);
                const int t = q + sp.M;
                pend = t / sp.L;
                q    = t % sp.L;
            }
        }
        st.pos = pos; st.q = q; st.pend = pend;
        return k;
    }

    // The same automaton, CLOCK ONLY, in closed form — for a gap (`nch == 0`), whose length the caller
    // chooses and can be two billion samples. Output j after the current state is due at push
    // pend + floor((q + j*M)/L), so of the next nIn pushes it emits every j with
    // floor((q + j*M)/L) <= nIn - pend, and resumes owing pend + floor((q + count*M)/L) - nIn. The
    // gap-equivalence test checks this against the loop above, sample for sample.
    long long advanceStage (int s, long long nIn) noexcept
    {
        const StagePlan& sp = thePlan.stage[s];
        const long long L = sp.L, M = sp.M;
        const long long q = phaseQ[s], pend = pendingIn[s];
        const long long T = sp.tapsPerPhase;
        ringPos[s] = (int) (((ringPos[s] - nIn) % T + T) % T);
        if (nIn < pend) { pendingIn[s] = (int) (pend - nIn); return 0; }
        const long long count = ((nIn - pend + 1) * L - 1 - q) / M + 1;
        const long long t = q + count * M;
        pendingIn[s] = (int) (pend + t / L - nIn);
        phaseQ[s]    = (int) (t % L);
        return count;
    }

    void advanceClock (long long n) noexcept
    {
        for (int s = 0; s < thePlan.count; ++s) n = advanceStage (s, n);
    }

    int runChain (int nch, int nIn, float* const* out, int outOffset) noexcept
    {
        int n = nIn;
        for (int s = 0; s < thePlan.count; ++s)
        {
            const StageState entry { phaseQ[s], pendingIn[s], ringPos[s] };
            StageState exitState = entry;
            int produced = 0;
            for (int c = 0; c < nch; ++c)
            {
                StageState st = entry;
                const double* src = (s == 0) ? inScratch[(std::size_t) c].data()
                                             : scratch[(std::size_t) ((s - 1) * channels + c)].data();
                produced = runStage (s, c, src, n, scratch[(std::size_t) (s * channels + c)].data(), st);
                exitState = st;
            }
            phaseQ[s] = exitState.q; pendingIn[s] = exitState.pend; ringPos[s] = exitState.pos;
            n = produced;
        }
        const int last = thePlan.count - 1;
        for (int c = 0; c < nch; ++c)
        {
            const double* src = scratch[(std::size_t) (last * channels + c)].data();
            for (int i = 0; i < n; ++i) out[c][outOffset + i] = (float) src[i];
        }
        return n;
    }

    Plan   thePlan {};
    bool   prepared = false;
    bool   flushed  = false;
    int    channels = 0, maxBlock = 0, flushInputs = 0, activeWidth = 0;
    int    phaseQ[kMaxStages] {};
    int    pendingIn[kMaxStages] {};
    int    ringPos[kMaxStages] {};
    std::vector<std::vector<double>> coef;        // [stage]
    std::vector<std::vector<double>> hist;        // [stage * channels + ch]
    std::vector<std::vector<double>> scratch;     // [stage * channels + ch]
    std::vector<std::vector<double>> inScratch;   // [ch]
    std::vector<int> stageBound;
};

} // namespace felitronics::core

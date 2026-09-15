// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/core/Math.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <limits>
#include <vector>

namespace felitronics::convolution
{

namespace detail
{
    // Modified Bessel I0 (series) — for the Kaiser window. Converges fast for moderate beta.
    inline double besselI0 (double x) noexcept
    {
        double sum = 1.0, term = 1.0;
        const double y = x * x * 0.25;
        for (int k = 1; k < 64; ++k)
        {
            term *= y / ((double) k * (double) k);
            sum  += term;
            if (term < 1e-17 * sum) break;
        }
        return sum;
    }

    // THE PHASE MEMO — why one output sample's window IS another's, to the bit.
    //
    // A tap's weight depends on the output position only through xx = t - k, and the window's 2R indices
    // are c + j for the same offsets j about c = floor(t). The whole argument rests on ONE step being
    // exact: frac = t - c. When it is, t IS c + frac as a real number, so the exact real behind every
    // xx is (c + frac) - (c + j) = frac - j — the SAME real for two outputs that share a frac, whatever
    // their c — and IEEE subtraction is correctly rounded, so both compute the same double from it. It
    // does NOT matter whether xx itself is representable; it matters that both roundings round the same
    // number. Then the 2R weights agree, and so does their sum, accumulated in the same order. Reusing
    // them moves no bit.
    //
    // `frac = t - c` is exact for every c >= 0: at c = 0 it is t itself, and above it Sterbenz's lemma
    // applies (c <= t < c + 1 <= 2c). BELOW zero it is not, and that is not hypothetical — at c = -1,
    // frac = t + 1 has to carry t's last bits at four times t's magnitude and loses them. Measured: at
    // 8 -> 48 kHz, outputs 0 and 6 have the same frac to the bit, yet their window arguments differ by
    // 4.4e-16 at tap 28 of 64. So c >= 0 is the gate, and the outputs whose CENTRE falls before input
    // sample 0 take the long road — one of them at 48 -> 44.1 kHz, nine at 11.025 -> 192.
    //
    // The argument has a second precondition, quieter than the first: every tap index must itself be an
    // exact double, or the real behind xx stops being frac - j — past 2^53 the integer k would round.
    // The position guard in `resampleIr` keeps c + R inside int, three decades below that, so it holds;
    // but it holds by that guard and not by luck, which is one reason the guard may not be relaxed.
    //
    // It pays because audio rates are small rationals. The phase of t repeats with period
    // outSr / gcd(inSr, outSr), and rounding spreads each of those onto a handful of neighbouring
    // doubles, so a one-second 48 -> 44.1 kHz resample has 913 distinct phases for 44100 outputs and
    // reuses 98% of its windows; 48 -> 96 kHz has two. Measured over every standard rate pair at a
    // four-second IR, the mean reuse is 99%. One standard pair does NOT fit: 11.025 -> 192 kHz has
    // 16686 distinct phases against a budget of 16131 rows at R = 32, so its store fills and the rest
    // of its outputs miss. It keeps what it stored (the hits it does get keep `live` true), and it is
    // the pair nobody loads — an 11 kHz IR on a 192 kHz host.
    //
    // A ratio with no period — a corrupt file rate, 48000 -> 44101 — reuses NOTHING, and the memo may
    // not make that case slower: it gives up after kGiveUpPhases MISSES without a single hit. That
    // threshold is above the worst standard period (2560, for 11.025 -> 192 kHz), so no real pair
    // trips it, and the futile case pays one fill of at most 2 MiB (0.7% of its own runtime) and then
    // stops probing altogether.
    //
    // THE MEMO IS A CACHE AND NOTHING ELSE. Every value in it was computed by the loop below, and a miss
    // computes the same value the same way, so the result cannot depend on the memo's size, its hash or
    // its give-up. What the suite checks is the consequence: the same resample run against a frozen copy
    // of the unmemoized kernel, output float by output float, over cases that reach the hit path, the
    // give-up and the full-store path alike.
    struct PhaseMemo
    {
        // 8 MiB of STORED weights — what a full store holds, not what the object costs: `store` is a vector,
        // so its capacity overshoots (9.85 MB measured at R = 300), and the two hash arrays are extra.
        static constexpr int kBudgetDoubles = 1 << 20;
        static constexpr int kGiveUpPhases  = 4096;      // misses tolerated before a ratio is declared aperiodic

        PhaseMemo (int strideIn, int outLen) : stride (strideIn)
        {
            maxPhases = (int) std::min<long long> ((long long) outLen, (long long) kBudgetDoubles / stride);
            if (maxPhases < 1) maxPhases = 1;
            int b = 1;
            while (b < 2 * maxPhases) b <<= 1;                            // TWICE the rows, and that is load-bearing:
            // at most maxPhases keys ever land in b buckets, so the table is never more than half full and the
            // linear probe below always meets an empty bucket. Size it at `maxPhases` and the probe spins for
            // ever on a full table — a hang, not a wrong answer, which only a timeout will catch.
            key.assign ((std::size_t) b, kEmpty);
            at.assign ((std::size_t) b, 0);
            mask = (unsigned) (b - 1);
        }

        // The bit pattern of a phase is its key; an all-ones double (a NaN) marks an empty bucket, and a
        // phase — which lies in [0, 1) — can never spell one.
        static constexpr unsigned long long kEmpty = ~0ULL;

        struct Probe { const double* row = nullptr; unsigned slot = 0; unsigned long long k = 0; bool storable = false; };

        Probe probe (double frac) noexcept
        {
            Probe p;
            if (! live) return p;
            std::memcpy (&p.k, &frac, sizeof (double));
            unsigned long long h = p.k + 0x9e3779b97f4a7c15ULL;          // splitmix64's finalizer
            h = (h ^ (h >> 30)) * 0xbf58476d1ce4e5b9ULL;
            h = (h ^ (h >> 27)) * 0x94d049bb133111ebULL;
            p.slot = (unsigned) (h ^ (h >> 31)) & mask;
            while (key[p.slot] != kEmpty && key[p.slot] != p.k) p.slot = (p.slot + 1) & mask;
            if (key[p.slot] == p.k) { p.row = &store[(std::size_t) at[p.slot] * (std::size_t) stride]; ++hits; }
            else
            {
                p.storable = (int) (store.size() / (std::size_t) stride) < maxPhases;
                // COUNTED HERE, NOT IN keep(): a miss that cannot be stored is still a miss, and once the
                // store is full there are no more insertions to count. Counting insertions made the give-up
                // unreachable for every halfTaps >= 128, where the budget stops the store at 4080 rows —
                // below the threshold — so an aperiodic ratio went on hashing and probing to the last output.
                if (hits == 0 && ++misses >= kGiveUpPhases) { live = false; p.storable = false; }
            }
            return p;
        }

        void keep (const Probe& p, const double* row)
        {
            at[p.slot] = (int) (store.size() / (std::size_t) stride);
            key[p.slot] = p.k;
            store.insert (store.end(), row, row + stride);
        }

        std::vector<unsigned long long> key;
        std::vector<int> at;
        std::vector<double> store;
        long long hits = 0, misses = 0;
        int stride = 0, maxPhases = 0;
        unsigned mask = 0;
        bool live = true;
    };
}

// EVERY FIELD HAS A RANGE. A value that is not a REQUEST — a NaN, a negative beta, a radius below one —
// is repaired to a documented default; a value that IS a request but lies past a CEILING is refused, and
// `resampleIr` returns empty (see "THE CONFIG" in the function). The ceilings are not taste:
//
//  * `kMaxBeta` is where the Bessel series in `detail::besselI0` STOPS MEETING ITS OWN convergence test.
//    That series runs at most 64 terms and quits early when a term falls below 1e-17 of the sum; the
//    largest beta that still quits by k = 63 is 53.038057, measured. Past it the loop simply runs out
//    of terms and returns a number that is silently wrong — not an infinity, not a NaN, just wrong:
//    1.20e-12 relative at beta 64, 9.3e-5 at 90, and the whole window shape gone by 200. (The series
//    does eventually overflow to +inf and make every window inf/inf = NaN, but not until beta ≈ 1.36e4
//    — far past where its answers stopped being answers, which is why the gate is a RANGE and not an
//    `isfinite`. An isfinite gate is what was here, and beta = 1e300 was the value it let through.)
//    A Kaiser beta of 53 is a ~490 dB stopband (beta/0.1102 + 8.7); float32 carries 144. Nothing
//    real is refused.
//  * `kMaxHalfTaps` bounds ONE output's tap loop. Without it `halfTaps` = 1e9 is two billion iterations
//    per output sample on the message thread. The ceiling's job is to stop an unbounded loop, not to
//    pick a design point, and 4096 is far past where the choice stops mattering: a Kaiser transition is
//    about (A-8)/(2.285*2R) rad/sample, so 8192 taps at beta 8 transition in ~29 Hz at 48 kHz against
//    the default 64 taps' ~3.8 kHz — 130 times narrower than the kernel every load actually uses.
//  * `halfTaps` < 1 is not a narrow kernel, it is NO kernel (an empty tap loop, an all-zero "IR"), so it
//    is repaired UP to 1 — the falsification test that pins that is ResamplerTests.cpp's "sanitizes a
//    degenerate config". The CEILINGS are not the same rule at the other end: a radius of 8192 or a beta of
//    90 is a meaningful request, and answering it with a smaller kernel would be a different filter than the
//    caller asked for, so `resampleIr` REFUSES those instead of repairing them. Repair what is not a
//    request; refuse a request we will not serve.
//
// A WORK BUDGET WAS CONSIDERED AND DECLINED. `kMaxHalfTaps` bounds one window, not the job: 2R*outLen at
// the ceilings is still 1.4e11 tap evaluations. A limit on the PRODUCT would bound the job — and would also
// refuse a legitimately long IR at an ordinary radius, which is exactly what `kMaxResampleSamples` is there
// to allow. Two bounds that contradict each other are worse than one that is honest about what it bounds.
// THE ALLOCATION BOUND (see `resampleIr`): 16.7M samples is 64 MiB of float per channel, 87 seconds at
// 192 kHz, 349 at 48 kHz. It is a SEPARATE constant from `MatrixConvolverNupc::kMaxIrSamples` — that one
// bounds a schedule's address arithmetic in a convolution backend this file knows nothing about — but their
// equality is deliberate and is an INVARIANT, not a coincidence: the resampler can produce anything the
// convolver can hold. Set it lower and there would be IR lengths the backend accepts and the loader cannot
// deliver; set it higher and the extra is memory nothing downstream can use.
inline constexpr int kMaxResampleSamples = 1 << 24;

struct IrResampleConfig
{
    static constexpr int    kMaxHalfTaps = 4096;   // window radius ceiling — see above
    static constexpr double kMaxBeta     = 53.0;   // where the Bessel series stops converging — see above

    int    halfTaps    = 32;     // window radius in input samples (64-tap filter) — more = sharper transition.
                                 // Below 1 → repaired to 1; above kMaxHalfTaps → the call is REFUSED
    double beta        = 8.0;    // Kaiser beta (~80 dB stopband); 5.65 ≈ 60 dB. NaN or negative → repaired to
                                 // 8.0; above kMaxBeta (+inf included) → the call is REFUSED
    double cutoffScale = 0.95;   // fraction of the (lower) Nyquist used as the passband edge. Outside (0,1] → 0.95
};

//==============================================================================
// A RESAMPLED IR IS NOT A RESAMPLED SIGNAL, and this is the factor that tells the two apart.
//
// `resampleIr` below normalizes every output tap to unity DC, which preserves the WAVEFORM's
// AMPLITUDE: tap for tap, the output carries the values the input would have had at the new rate.
// That is what a signal wants. A convolution's gain is not an amplitude — it is a SUM over taps, so
// it is proportional to how many of them fit into a second. The same response resampled from `inSr`
// to `outSr` therefore convolves `outSr / inSr` times as loud: +6.02 dB for a 48 kHz IR on a 96 kHz
// host, +5.28 at 88.2, -0.74 at 44.1, +12.04 at 192. Measured on a real 48 kHz spring IR the reading
// matches that ratio to four decimals at every standard rate, and multiplying by the factor below
// puts all five within 0.0002 dB of each other (P68).
//
// SO WHOEVER HANDS THE RESULT TO A CONVOLVER MULTIPLIES BY THIS. `CabConvolver` does it for the IRs
// it resamples itself; a caller that resamples an IR on its own must do it too, or its wet path
// changes level with the host's clock. It is exactly 1 when the rates match, so a load that is not
// resampled never sees it, and 1 is also what an unusable pair of rates gets: a factor that is not a
// positive finite number is not a compensation, and silencing or blasting an IR is a worse answer
// than leaving it at the level it came with (the rule P67 settled for broken metadata).
//
// It does NOT make the resample gain-exact by itself. The density term is all it takes out; the
// kernel's pre-ringing that would fall BEFORE output sample 0 is still dropped (see resampleIr's
// note below), and how much that costs depends on how abruptly the IR starts and on how far its
// first energy sits from sample 0 — not on this factor.
//
// AND IT IS A RIPPLE, NOT A LOSS — the sign alternates, because the kernel's nearest pre-ring lobes
// are negative and cutting them ADDS level. A lone impulse (the sharpest onset there is) at input
// index `lead`, DC gain after the factor, in dB:
//
//     lead:            0      1      2      3      4      6      8     12     16     32
//     96 -> 48 kHz  -2.343 +0.827 +0.781 -0.340 -0.515 +0.337 -0.247 -0.119 -0.049 -0.000
//     96 -> 44.1    -2.574 +0.653 +0.949 -0.049 -0.606 +0.327 -0.161 +0.021 +0.050 +0.000
//     48 -> 44.1    -0.484 +0.512 -0.494 +0.389 -0.303 -0.101 +0.035 +0.063 -0.010 +0.000
//     48 -> 192     -0.697 +0.241 -0.181 +0.153 -0.138 -0.111 -0.084 -0.038 -0.009 -0.000
//
// So an onset-trimmed IR can come out nearly a dB HOT as easily as quiet, and a measurement that
// trims the front must not charge the difference here IN EITHER DIRECTION. It dies at `halfTaps`
// input samples of lead — the window's own backward reach — and only there. A smoother onset pays
// far less: a one-pole `exp(-n/tau)` starting at full scale with no lead loses 0.21 dB at 44.1 kHz
// and 0.39 at 192 for tau = 1 INPUT SAMPLE (8.7 dB per sample, an impulse with a smear), but only
// 0.004 and 0.011 dB at tau = 1 ms, which is what a real decay looks like. Pad the front if it
// matters; the numbers above are why a measurement of this wants a lead of at least `halfTaps`.
// TWO GATES, AND BOTH EARN THEIR KEEP. Guarding only the RESULT is not the same rule and was wrong here:
// -48 kHz against -96 kHz divides to a perfectly finite 0.5, so a pair of NEGATIVE rates — which
// `resampleIr` refuses outright — would have been compensated for a resample that never ran. Guarding only
// the inputs is not enough either: two finite rates 600 decades apart still divide to an infinity or a
// zero. Positivity covers NaN as well (every comparison against it is false), and it is the result gate
// that answers for the infinities: inf/48000 is inf, 48000/inf is 0, inf/inf is a NaN, and all three leave
// through the same door, which is why there is no `isfinite` on the inputs — it would be a third spelling
// of a rule already stated twice.
//
// It is `resampleIr`'s RATE gate, not the loader's. `CabConvolver` decides "no resample" by its own
// `kRateMatchTolerance` — two rates a part per million apart are ONE rate there and this is never
// consulted — and `resampleIr` has geometry gates further down (a length or a tap position `int` cannot
// address) that this does not repeat, so it can still answer 1e12 for a pair whose resample would be
// refused. Nothing is scaled in that case because nothing is staged: the loader gives up on the empty
// result before it reaches a gain.
[[nodiscard]] inline double convolutionRateGain (double inSr, double outSr) noexcept
{
    if (! (inSr > 0.0) || ! (outSr > 0.0)) return 1.0;                  // NaN, zero and negatives, both sides
    const double g = inSr / outSr;
    return (std::isfinite (g) && g > 0.0) ? g : 1.0;                    // ...and every infinity, plus over/underflow
}

//==============================================================================
// Offline windowed-sinc (Kaiser) IR resampler. MESSAGE-THREAD ONLY (double math, allocates) — for
// rate-converting an impulse response to the host SR on load. DC gain is normalized to 1.
//
// OUTSIDE THE INPUT IS SILENCE. Every output sample divides by the weight of its WHOLE window, and a tap
// that falls past either end of the input adds its weight but no signal: the IR is resampled as if it were
// padded with zeros, which is what it is — a cabinet trimmed at its onset had silence before it, and the
// convolver plays silence after its last tap. Those taps used to be skipped BEFORE their weight was added,
// so an edge sample was divided by only the part of its window that landed on the input, and the samples
// that were not there counted as the weighted mean of the ones that were. A cabinet's onset is its loudest
// edge; that invented a broadband floor over the top octave, which is exactly where a cabinet is quietest.
// What zeros cannot give back is the kernel's pre-ringing that would fall BEFORE output sample 0: an IR that
// starts at sample 0 keeps it only by adding delay, and it costs more the more abruptly the IR starts.
// `felitronics_convolution_resampler_tests` owns the numbers — shift invariance, and a cabinet-like IR's
// band response against its own response and against the untruncated resample.
//
// A LOAD IS NEVER RESAMPLED TO NOTHING. The length is inLen*ratio rounded but at least one sample, as JUCE's
// resampleImpulseResponse had it: a one-tap IR at 96 -> 44.1 kHz rounded to zero taps, and the loader had
// nothing to publish. The result is empty only for no input (a null pointer or a non-positive length), a rate
// that is not a positive finite number, an output `int` cannot address — longer than INT_MAX samples, or a
// ratio so small that output sample 0 alone sits past INT_MAX (the floor of one sample is what makes that
// reachable: `(int) floor(t)` would be undefined) — an output longer than `kMaxResampleSamples`, or a config
// asking for a kernel past `IrResampleConfig`'s ceilings.
//
// AND IT IS NEVER RESAMPLED TO GIGABYTES. The length gate above only kept the arithmetic addressable: a
// length of INT_MAX is 8.6 GB asked of the heap in ONE call, on the message thread, and a file rate does
// not have to be absurd to get there — 1.2 Hz against a 48 kHz host is a ratio of 40000, and a one-second
// IR then asks for 7.7 GB. `kMaxResampleSamples` is the bound that was missing. An output past it is
// REFUSED, not truncated: that is the rule P67 already ratified for a rate whose result cannot be
// addressed, and the caller sees the same empty vector for the same reason. It costs nothing real — at
// 192 kHz the ceiling is 87 seconds of taps, and the convolver downstream caps a cab at four.
//
// FAMILY SPLIT vs core::StreamResampler — restated, because the other half of it changed under this
// comment. That one used to be a Catmull-Rom cubic and "too low-SNR for IRs" was the whole argument.
// Since P34 it is a 64-tap polyphase windowed sinc, i.e. the SAME family as this one, so the split is
// no longer about quality. It is about BUDGET and THREAD: this is an offline one-shot that may
// allocate, work in double and size its kernel to the job; that is a streaming rate-match on the audio
// thread with a fixed table built in reset() and a per-block cost that has to stay inside a couple of
// percent of a neural stage. Still not interchangeable — for the opposite reason to the one that used
// to be written here.
inline std::vector<float> resampleIr (const float* in, int inLen, double inSr, double outSr,
                                      IrResampleConfig cfg = {})
{
    std::vector<float> out;
    if (in == nullptr || inLen <= 0 || ! (inSr > 0.0) || ! (outSr > 0.0) || ! std::isfinite (inSr)
        || ! std::isfinite (outSr)) return out;

    const double ratio = outSr / inSr;
    if (! (ratio > 0.0) || ! std::isfinite (ratio)) return out;            // finite rates can still under/overflow here
    const double want = (double) inLen * ratio;
    if (! (want < (double) std::numeric_limits<int>::max())) return out;   // before the cast: a wrapped length is garbage
    const int outLen = std::max (1, (int) std::llround (want));            // never zero taps (see above)
    if (outLen > kMaxResampleSamples) return out;                          // the allocation bound, on the ROUNDED
                                                                           // length: `want` of M + 0.25 still fits
    // THE CONFIG: A VALUE THAT IS NOT A REQUEST IS REPAIRED, A REQUEST WE WILL NOT SERVE IS REFUSED.
    // A NaN, a negative beta, a window radius below one — none of those is a kernel anybody asked for (a
    // radius of zero is an EMPTY tap loop, i.e. an all-zero "IR"), so they become the documented default,
    // which is what this struct has always done and what ResamplerTests' falsification case pins. But a
    // radius of 8192, or a beta of 90, IS a request — a perfectly meaningful one that this implementation
    // will not serve — and quietly answering it with a 4096-tap kernel or with beta 8 would hand back a
    // different filter from the one asked for. A kernel is what a resampler sounds like, so that is refused
    // instead. One rule, two kinds of fact.
    // (+inf and a NaN part company here, and the comparison is exactly why: an infinity is ORDERED — it is
    // above the ceiling, the same answer every other value above the ceiling gets — while a NaN is nowhere
    // on the line at all, and a value that is not anywhere is not a request.)
    if (cfg.halfTaps > IrResampleConfig::kMaxHalfTaps) return out;
    if (cfg.beta     > IrResampleConfig::kMaxBeta)     return out;
    const int    R      = cfg.halfTaps < 1 ? 1 : cfg.halfTaps;
    const double beta   = cfg.beta >= 0.0 ? cfg.beta : 8.0;             // NaN and negatives → the default
    const double cScale = (std::isfinite (cfg.cutoffScale) && cfg.cutoffScale > 0.0 && cfg.cutoffScale <= 1.0)
                        ? cfg.cutoffScale : 0.95;
    const double fc     = 0.5 * std::min (1.0, ratio) * cScale;            // cycles per INPUT sample
    const double i0beta = detail::besselI0 (beta);

    // Every tap index is an int: the last output's input position plus the window radius has to fit.
    const double tLast = ((double) outLen - 0.5) / ratio - 0.5;
    if (! (tLast + (double) R + 2.0 < (double) std::numeric_limits<int>::max())) return out;
    out.assign ((std::size_t) outLen, 0.0f);

    const int taps = 2 * R;
    detail::PhaseMemo memo (taps + 1, outLen);
    std::vector<double> row ((std::size_t) (taps + 1));

    for (int n = 0; n < outLen; ++n)
    {
        const double t = ((double) n + 0.5) / ratio - 0.5;                 // output n → input position (sample-centred)
        const int    c = (int) std::floor (t);
        double acc = 0.0, wsum = 0.0;
        // c >= 0 is what makes `t - c` exact, and that is what makes a phase reproducible — see
        // detail::PhaseMemo. c < 0 means the output's CENTRE sits before input sample 0 (the window itself
        // reaches back further than that at every output, and always has); there are as many such outputs
        // as the ratio puts before the first input sample — one for 48 -> 44.1 kHz, nine for 11.025 -> 192
        // — and they take the long road, as does every output when a ratio has no period at all.
        const detail::PhaseMemo::Probe p = c >= 0 ? memo.probe (t - (double) c)
                                                  : detail::PhaseMemo::Probe {};
        if (p.row != nullptr)
        {
            wsum = p.row[(std::size_t) taps];                              // the same terms, summed in the same order
            for (int j = 0; j < taps; ++j)
            {
                const int k = c - R + 1 + j;
                if (k >= 0 && k < inLen) acc += (double) in[k] * p.row[(std::size_t) j];
            }
        }
        else
        {
            // THE MEMO RESTS ON THIS LOOP BODY BEING A FUNCTION OF `xx` AND THE PER-CALL CONSTANTS, AND
            // OF NOTHING ELSE. Reach for `t`, `k`, `n` or `c` directly in here — in a faster sinc, say —
            // and two outputs that share a phase stop sharing an answer, with no compile error anywhere
            // to say so. The differential test in the resampler suite is the only thing that would notice —
            // and only if the damage is bigger than a last place or two: a perturbation of about 1e-16
            // reaches a float32 result roughly once in a billion samples, which is why the gate above has a
            // hand-found witness rather than a sweep.
            for (int j = 0; j < taps; ++j)
            {
                const int    k    = c - R + 1 + j;
                const double xx   = t - (double) k;
                const double sinc = (std::fabs (xx) < 1e-12) ? (2.0 * fc)
                                                             : std::sin (2.0 * core::kPi * fc * xx) / (core::kPi * xx);
                const double r    = xx / (double) R;                       // window argument in [-1, 1]
                const double win  = (r <= -1.0 || r >= 1.0) ? 0.0
                                  : detail::besselI0 (beta * std::sqrt (1.0 - r * r)) / i0beta;
                const double w    = sinc * win;
                wsum += w;                                                 // the WHOLE window: OUTSIDE THE INPUT IS SILENCE
                if (k >= 0 && k < inLen) acc += (double) in[k] * w;
                if (p.storable) row[(std::size_t) j] = w;
            }
            if (p.storable) { row[(std::size_t) taps] = wsum; memo.keep (p, row.data()); }
        }
        out[(std::size_t) n] = (float) (! core::exactlyEqual (wsum, 0.0) ? acc / wsum : 0.0);    // normalize → unity DC (intentional exact ==)
    }
    return out;
}

// The `int` length is checked BEFORE it is narrowed: a vector longer than INT_MAX would wrap into a small
// positive length and resample a plausible-looking prefix of something the caller never asked about.
inline std::vector<float> resampleIr (const std::vector<float>& in, double inSr, double outSr,
                                      IrResampleConfig cfg = {})
{
    if (in.size() > (std::size_t) std::numeric_limits<int>::max()) return {};
    return resampleIr (in.data(), (int) in.size(), inSr, outSr, cfg);
}

} // namespace felitronics::convolution

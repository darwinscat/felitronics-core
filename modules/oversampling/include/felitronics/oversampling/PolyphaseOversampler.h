// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/core/Config.h>
#include <felitronics/core/Math.h>
#include <felitronics/core/PolyphaseFir.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace felitronics::oversampling
{

namespace detail
{
    inline double besselI0 (double x) noexcept
    {
        double sum = 1.0, term = 1.0;
        const double y = x * x * 0.25;
        for (int k = 1; k < 64; ++k) { term *= y / ((double) k * (double) k); sum += term; if (term < 1e-17 * sum) break; }
        return sum;
    }
}

//==============================================================================
// felitronics::oversampling::PolyphaseOversampler — integer up/down sampling by a polyphase
// windowed-sinc (Kaiser) FIR. Linear phase → a constant group delay (reported). For inter-sample-peak
// detection (the true-peak limiter) and, with up+down around a process, alias-free nonlinear processing.
// RT-safe: prepare() allocates the FIR + per-channel histories; up/downsample() do no alloc/lock/throw.
// Both inner loops are `core::firDot` — ONE kernel with a nailed-down summation order, scalar / SSE2 /
// NEON / wasm-SIMD128, the same bits on every row (P56). The coefficients are laid out for it in
// prepare(): contiguous for downsample, phase-major for upsample, both padded with +0.0f to a multiple
// of four; the histories became double-length backwards rings so the window is contiguous.
//
// WHY THE DEFAULT IS 64, AND WHAT SETS IT. The cutoff is FIXED at 0.90 x baseband Nyquist
// (designFilter() below), so the transition band has to fit between 0.45 fs and the fold at 0.50 fs,
// and tapsPerPhase is the only thing that decides whether it does. The design DECLARES its own target
// one line down — `beta = 9.0` is a ~90 dB Kaiser stopband — and everything above 0.50 fs folds
// straight back into the audio band, so that number is the filter's job rather than a nicety.
// Worst |H| over the whole fold region [0.5, factor/2] fs — the WORST of factors 2, 4 and 8:
//
//        tpp      32       40       50       56       57      58      59      60      64      72
//        dB    -26.9    -37.0    -55.7    -76.2    -82.6   -90.7   -89.8   -90.7   -90.5   -91.9
//
// Read that carefully, because it is NOT a step. Below ~58 the transition is genuinely unfinished and
// the rejection is monotonically poor. At ~58 it reaches the Kaiser window's floor and then RIPPLES
// there, +-1 dB, as the last sidelobe slides relative to the fold: 59 (-89.85) and 61 (-89.89) sit a
// tenth of a dB the wrong side of the declared 90, while 63 and 65 are near -91.3. So "the first taps
// count meeting -90.0" is 58, not 60, and that integer is a knife-edge artefact of testing a rippling
// quantity against a hard bar — what the sweep establishes is where the transition FINISHES. THE
// DEFAULT USED TO BE 32, i.e. the prototype delivered 27 dB where it promised 90; 64 sits at -90.46,
// past the knee and on the floor. Above the knee further taps buy pass-band width, not rejection.
//
// End to end on a Saturator (0.17 fs tone, 4x), TOTAL non-harmonic energy — every bin that is not a
// harmonic BELOW Nyquist, so a harmonic that folded to get there counts as the aliasing it is:
//
//        drive     +0       +6      +12      +18      +24      +30      +36  dB
//        32     -135.3    -60.1    -48.2    -44.3    -31.8    -25.3    -22.9  dBc
//        64     -130.8   -132.3   -101.0    -50.3    -31.3    -24.8    -22.4  dBc
//
// At the drive this stage is actually used at — its own Params doc calls 1..6 dB the mastering range —
// the taps remove 50 to 72 dB of aliasing. They stop helping above about +24 dB, where the harmonic
// series reaches past the OS Nyquist and folds INSIDE the oversampled domain, which no decimation
// filter can reach: that is the FACTOR's axis, and there the total is 0.5 dB worse, because a flatter
// pass band also delivers what had already folded. The component the taps own at every drive is the
// transition-band leakage — the 3rd harmonic of a 0.17 fs tone sits at 0.51 fs and folds to 0.49 fs,
// where 32 taps left it at -45.5 dBc and 64 puts it at -148.0, i.e. into the float noise.
// (That last number read -139.4 until P56 replaced the inner loop with `core::firDot`: four partial
// accumulators take a quarter of the sequential sum's rounding steps each, so the arithmetic noise
// floor this measurement was reading dropped ~8.6 dB. The FILTER did not change — the measurement was
// measuring the summation, which is what a residual 90 dB under the design's own stopband does.)
// Pass-band flatness is a COROLLARY of the same width, not a second requirement: two oversampled
// stages in series (the real assembly — a clipper in front of a limiter) then droop under 0.04 dB to
// 0.41 fs, against -3.25 dB at 32. The pass-band edge follows 0.45 fs - c/tpp with c = 2.385 (0.1 dB),
// 2.505 (0.05 dB), 2.723 (0.01 dB), so a WIDER flat band is bought hyperbolically: 0.42 fs costs 80
// taps, 0.43 costs 128, 0.44 costs 239. Do not buy it that way — 20 kHz at 44.1 kHz is 0.4535 fs, ABOVE
// the fixed cutoff, so no tapsPerPhase reaches it and more taps make it WORSE (-27.4 dB at 32 against
// -31.1 at 64). Widening the audible band is the CUTOFF's axis, not this one.
// Cost, both linear in tapsPerPhase: the FIR (measured 2.09x for the clipper+limiter pair, 2.13 %RT ->
// 4.44 %RT at 48 kHz stereo 4x) and the reported latency (tapsPerPhase - 1 baseband samples).
class PolyphaseOversampler
{
public:
    // The shipped topology decision, named so a consumer can state "the core's default" rather than
    // re-spell the number — and so a test can assert the DEFAULT moved, not just that some number did.
    static constexpr int kDefaultTapsPerPhase = 64;
    // BOTH factors of N = factor*tapsPerPhase are bounded, HERE rather than at each caller, because the
    // product is what overflows and either argument alone can do it. tapsPerPhase = INT_MAX overflowed the
    // multiply below and then threw a length_error out of assign() — a terminate under the wasm tier's
    // -fno-exceptions. TruePeakLimiter::prepare refused both at its own gate (its contract, kept), but
    // Saturator passes an unbounded oversampleFactor straight through, and `prepare(INT_MAX, 1)` on the
    // default taps is UBSan-confirmed signed overflow followed by that same length_error. The ceiling is
    // generous rather than tight: the largest factor anything in the tree asks for is 32 (PowerAmpStage's
    // alias-free reference), and 64 x 1024 taps is still a trivially small allocation.
    static constexpr int kMaxTapsPerPhase = 1024;
    static constexpr int kMaxFactor       = 64;

    // WHAT prepare() ASKS THE HEAP FOR (law 11d) — the one function it sizes itself with, so a caller
    // budgeting memory reads the counts the five buffers are actually built from. FALSE, with `out`
    // untouched, exactly where prepare() refuses the same arguments: the budget of a call is storageFor()
    // with the SAME arguments, and a refused prepare() allocates nothing.
    //
    // ⚠ IT MODELS THE CHANNEL CLAMP, because prepare() clamps rather than refuses — `prepare(4, 17, 64)`
    // answers true and sizes itself for `core::kMaxChannels`. That is a law-11(b) violation older than
    // this budget (a width it clamped is a width it lied about) and it is **P55**, not this function's to
    // fix: a budget that refused where the call succeeds would be wrong about the call in front of it.
    // ⚠ THE COUNTS MOVED WITH P56, and the reason is `core::firDot`, not a bigger filter:
    // the histories are DOUBLE-LENGTH rings (every sample stored twice, so the window is contiguous)
    // and the coefficients exist twice, once contiguous for `downsample` and once phase-major for
    // `upsample`. Both lengths are rounded up to a multiple of four, which is what lets the kernel run
    // without a tail. Nothing here is a heuristic: it is exactly what prepare() below assigns.
    struct Storage
    {
        std::size_t proto = 0, protoPhase = 0, upHist = 0, downHist = 0;   // floats
        std::size_t upPos = 0, downPos = 0;                                // ints
        std::uint64_t bytes() const noexcept
        {
            return (std::uint64_t) sizeof (float) * ((std::uint64_t) proto + protoPhase + upHist + downHist)
                 + (std::uint64_t) sizeof (int)   * ((std::uint64_t) upPos + downPos);
        }
        // Does an oversampler already holding `other` have room for this without asking the heap again?
        // It lives HERE rather than being spelled out field by field at each caller, because a caller's
        // copy of the list is a copy that can go stale: the one in `MasteringChain::Storage::fitsWithin`
        // silently omitted `downPos`, and would have omitted `protoPhase` the moment P56 added it.
        bool fitsWithin (const Storage& other) const noexcept
        {
            return proto <= other.proto && protoPhase <= other.protoPhase
                && upHist <= other.upHist && downHist <= other.downHist
                && upPos <= other.upPos && downPos <= other.downPos;
        }
    };

    [[nodiscard]] static bool storageFor (int factor, int maxChannels, int tapsPerPhase, Storage& out) noexcept
    {
        if (factor < 2 || factor > kMaxFactor) return false;
        if (tapsPerPhase < 4 || tapsPerPhase > kMaxTapsPerPhase) return false;
        const std::size_t ch   = (std::size_t) channelsFor (maxChannels);
        const std::size_t tp   = (std::size_t) core::firPadLen (tapsPerPhase);
        const std::size_t np   = (std::size_t) core::firPadLen (factor * tapsPerPhase);
        out.proto      = np;
        out.protoPhase = (std::size_t) factor * tp;
        out.upHist     = ch * 2u * tp;
        out.downHist   = ch * 2u * np;
        out.upPos      = ch;
        out.downPos    = ch;
        return true;
    }

    // The width prepare() will actually run at — the clamp above, in ONE place so the budget and the
    // preparation cannot disagree about it. See P55.
    static int channelsFor (int maxChannels) noexcept
    {
        return maxChannels < 1 ? 1 : (maxChannels > core::kMaxChannels ? core::kMaxChannels : maxChannels);
    }

    // factor 2/4/8; tapsPerPhase = FIR taps per polyphase branch (filter length = factor*tapsPerPhase).
    // The default is the topology decision above; pass it explicitly to pin one against this default.
    bool prepare (int factor, int maxChannels, int tapsPerPhase = kDefaultTapsPerPhase)
    {
        Storage st;
        if (! storageFor (factor, maxChannels, tapsPerPhase, st)) return false;
        L = factor; tpp = tapsPerPhase; N = L * tpp;
        tppPad = core::firPadLen (tpp); nPad = core::firPadLen (N);
        channels_ = channelsFor (maxChannels);
        designFilter();
        upHist.assign   (st.upHist,   0.0f);
        downHist.assign (st.downHist, 0.0f);
        upPos.assign    (st.upPos,    0);
        downPos.assign  (st.downPos,  0);
        return true;
    }

    void reset() noexcept
    {
        std::fill (upHist.begin(),   upHist.end(),   0.0f);
        std::fill (downHist.begin(), downHist.end(), 0.0f);
        std::fill (upPos.begin(),    upPos.end(),    0);
        std::fill (downPos.begin(),  downPos.end(),  0);
    }

    // Clear ONE channel's histories, leaving every other channel BIT-EXACT — for an owner whose
    // channel stopped being fed and will be fed again. The ring POSITIONS go with the samples: they
    // are state like any other, and a channel whose ring is zeroed while its cursor stays put sits at
    // a phase reset() can never produce, so "this column now equals a fresh one" would be false about
    // a column that looks clean. Unprepared (channels_ == 0) is a no-op, not a bad index.
    void resetChannel (int c) noexcept
    {
        if (c < 0 || c >= channels_) return;
        std::fill_n (upHist.begin()   + (std::ptrdiff_t) c * (std::ptrdiff_t) (2 * tppPad), 2 * tppPad, 0.0f);
        std::fill_n (downHist.begin() + (std::ptrdiff_t) c * (std::ptrdiff_t) (2 * nPad),   2 * nPad,   0.0f);
        upPos[(std::size_t) c]   = 0;
        downPos[(std::size_t) c] = 0;
    }

    int factor() const noexcept { return L; }
    // Group delay of ONE filter pass, in OVERSAMPLED samples ((N-1)/2 for a linear-phase FIR).
    double filterLatencyOversampled() const noexcept { return (double) (N - 1) * 0.5; }
    // Round-trip (up THEN down) latency in BASEBAND samples. EXACT integer: the two (N-1)/2 OS group
    // delays plus the decimation phase (L-1) combine to (N - L)/L = tpp - 1 baseband samples.
    // Unprepared (tpp == 0) reports 0, not -1 — hosts query latency before prepare().
    int latencySamples() const noexcept { return tpp > 0 ? tpp - 1 : 0; }

    // n baseband samples per channel → n*L oversampled samples. out[ch] holds n*L.
    //
    // The ring runs BACKWARDS: `pos` is where the NEWEST sample lives, so reading forward from it
    // already gives the oldest-coefficient-first order the FIR wants, and every sample is stored a
    // second time `tppPad` slots up so the window [pos, pos+tppPad) never wraps. That is what leaves
    // the inner loop free of a wrap test — and free of everything else: see core::firDot.
    void upsample (const float* const* in, int channels, int n, float* const* out) noexcept
    {
        const int nc = channels < channels_ ? channels : channels_;
        for (int c = 0; c < nc; ++c)
        {
            int   pos  = upPos[(std::size_t) c];
            float* h   = &upHist[(std::size_t) c * (std::size_t) (2 * tppPad)];
            for (int m = 0; m < n; ++m)
            {
                if (pos == 0) pos = tppPad;
                --pos;                                            // pos = newest x[m]
                h[pos] = in[c][m];
                h[pos + tppPad] = in[c][m];
                for (int p = 0; p < L; ++p)
                {
                    const float acc = core::firDot (&protoPhase[(std::size_t) p * (std::size_t) tppPad],
                                                    &h[pos], tppPad);
                    out[c][m * L + p] = (float) L * acc;          // *L restores the zero-stuff gain
                }
            }
            upPos[(std::size_t) c] = pos;
        }
    }

    // n*L oversampled samples per channel → n baseband samples (lowpass + decimate). out[ch] holds n.
    // Same backwards double-length ring as upsample; the coefficients need no repacking here, they
    // were already contiguous — only the zero padding to a multiple of four is new.
    void downsample (const float* const* in, int channels, int n, float* const* out) noexcept
    {
        const int nc = channels < channels_ ? channels : channels_;
        for (int c = 0; c < nc; ++c)
        {
            int   pos = downPos[(std::size_t) c];
            float* h  = &downHist[(std::size_t) c * (std::size_t) (2 * nPad)];
            for (int m = 0; m < n; ++m)
            {
                for (int p = 0; p < L; ++p)
                {
                    if (pos == 0) pos = nPad;
                    --pos;
                    h[pos] = in[c][m * L + p];
                    h[pos + nPad] = in[c][m * L + p];
                }
                out[c][m] = core::firDot (proto.data(), &h[pos], nPad);   // proto sums to 1 → unity passband
            }
            downPos[(std::size_t) c] = pos;
        }
    }

private:
    void designFilter()
    {
        proto.assign ((std::size_t) nPad, 0.0f);    // == Storage::proto: N taps then the +0.0f padding
        const double fc   = 0.5 / (double) L * 0.90;              // cutoff (cycles/OS-sample), guard below baseband Nyquist
        const double cen  = (double) (N - 1) * 0.5;
        const double beta = 9.0;                                  // Kaiser ~ -90 dB stopband
        const double i0b  = detail::besselI0 (beta);
        double sum = 0.0;
        for (int i = 0; i < N; ++i)
        {
            const double x    = (double) i - cen;
            const double sinc = (std::fabs (x) < 1e-9) ? (2.0 * fc)
                                                       : std::sin (2.0 * core::kPi * fc * x) / (core::kPi * x);
            const double r    = (double) (2 * i - (N - 1)) / (double) (N - 1);   // ∈ [-1,1]
            const double win  = detail::besselI0 (beta * std::sqrt (std::max (0.0, 1.0 - r * r))) / i0b;
            const double v    = sinc * win;
            proto[(std::size_t) i] = (float) v;
            sum += v;
        }
        const float inv = (float) (1.0 / sum);                   // normalize Σ → 1 (unity DC)
        for (auto& v : proto) v *= inv;                          // the padding is 0 and stays 0

        // The phase-major copy upsample() reads: protoPhase[p][k] == proto[k*L + p], contiguous in k
        // and zero-padded to tppPad. The strided gather that used to sit inside the inner loop — once
        // per output sample, per phase, per tap — now happens exactly once, here, at prepare() time.
        protoPhase.assign ((std::size_t) L * (std::size_t) tppPad, 0.0f);
        for (int p = 0; p < L; ++p)
            for (int k = 0; k < tpp; ++k)
                protoPhase[(std::size_t) p * (std::size_t) tppPad + (std::size_t) k]
                    = proto[(std::size_t) (k * L + p)];
    }

    int L = 0, tpp = 0, N = 0, channels_ = 0;
    int tppPad = 0, nPad = 0;                  // tpp and N rounded up to a multiple of 4 (core::firDot)
    std::vector<float> proto;                  // N taps, Σ = 1, then +0.0f padding to nPad
    std::vector<float> protoPhase;             // the same taps phase-major: L blocks of tppPad
    std::vector<float> upHist, downHist;       // per-channel DOUBLE-LENGTH rings, newest-first
    std::vector<int>   upPos, downPos;         // index of the NEWEST sample in its ring
};

} // namespace felitronics::oversampling

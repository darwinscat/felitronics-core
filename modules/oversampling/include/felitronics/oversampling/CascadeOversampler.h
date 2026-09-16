// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/core/Config.h>
#include <felitronics/core/DetMath.h>
#include <felitronics/core/Math.h>
#include <felitronics/core/PolyphaseFir.h>
#include <felitronics/oversampling/PolyphaseOversampler.h>   // detail::besselI0 — one copy of the window

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace felitronics::oversampling
{

//==============================================================================
// felitronics::oversampling::CascadeOversampler — a STRICT, linear-phase up/down sampler that is flat to
// 20 kHz at 44.1 kHz, for the nonlinear stages (P31). Same job, same call shape and same kernel
// (`core::firDot`) as `PolyphaseOversampler`; a different filter, and a different price.
//
// WHY IT EXISTS. `PolyphaseOversampler` centres a fixed cutoff at 0.45 fs so that its transition band
// finishes below the fold — the guard band its header derives. At 44.1 kHz that cutoff sits BELOW the top
// of the audio band: one round trip reads -1.80 dB at 19 kHz and -15.55 at 20 kHz, and no taps count fixes
// it. This class keeps the guard (STRICT: the image of content anywhere below fs/2 is rejected, both
// ways) and moves the band edge instead, by paying for the one thing that buys it: a transition only
// fs/2 - 20 kHz wide.
//
// THE STRUCTURE, and why it is not the cheaper one. Factor 2^S as S stages of 2x:
//   * stage 1 (fs -> 2fs) is an ordinary Kaiser windowed sinc — the only stage whose transition is narrow,
//     so the only expensive one, and at 2fs rather than at F*fs it costs half of what one 4x stage would;
//   * stages 2..S are HALFBANDS — every other tap exactly zero, one inner product per output sample of
//     half the length — because their job is easy: content is already below fs/2, so their images sit
//     at >= 1.5 fs and the transition may be wide.
// Stage 1 is NOT a halfband, and cannot be: a halfband's transition is centred ON its quarter rate (for
// stage 1 that is fs/2 itself), since H(f) + H(Fs/2 - f) = 1 makes its image rejection at fs - a equal to
// its pass-band deviation at a. A pair of halfbands (79 + 23 taps) that loses 0.41 dB at 20 kHz leaves
// that tone's image at -32.5 dB (both pinned in OversamplingTests, "the cutoff axis"). What that costs
// through a waveshaper was measured by the P31 design stand, not by the suite: -34 dBc of aliasing at
// +6 dB of drive, against -141 for the fixed-cutoff design.
//
// THE DESIGN RULE, computed in designFor() from the sample rate and nothing else:
//   * band edge fp = 20 kHz, or 20/44.1 of fs below 44.1 kHz (so every rate up to 44.1 kHz gets the
//     44.1 kHz geometry, and a 32 kHz stream is not asked for a 20 kHz band it does not have);
//   * stage 1 (Kaiser, beta 9.5) takes the smallest tapsPerPhase t >= 16 with
//         (kPassA0 + kPassA1/t + kStopB0 + kStopB1/t) * fs / t  <=  fs/2 - fp
//     and puts its cutoff at fs/2 - (kStopB0 + kStopB1/t) * fs / t. The constants are measured, in units of
//     fs/t: one pass is within 0.005 dB at (kPassA0 + kPassA1/t) below the cutoff, and |H| stays under
//     -90 dB from (kStopB0 + kStopB1/t) above it, for every t from 12 to 240 (the fitted bounds sit above
//     the measured distances at every point). So the stop edge lands just under fs/2 and 20 kHz is flat.
//     The constants come from the P31 design stand; what the suite pins is their OUTCOME, at every one of
//     the 110 taps counts the rule can produce.
//   * beta is 9.5, not the 9 of PolyphaseOversampler, and the reason is the RULE: a beta-9 window's floor
//     ripples around -90 dB itself (that header's "59 taps is worse than 58"), so the -90 point jumps
//     between sidelobes as t moves (the stand read it wandering between about 2.9 and 3.6 fs/t) and no
//     smooth rule can place it. At 9.5 the floor is ~-94.6 and the -90 point sits on the skirt (3.01 to
//     3.13 fs/t on the same stand), at a cost the stand put at ~4 % of taps.
//   * halfband lengths by stage: 27, 23, 23, 15, 15 (beta 9.5) — the shortest that keep the cascade at the
//     stage-1 floor (the design stand read a 23-tap stage 2 at -80.5; the suite's strictness bar fails it).
//   * decimation phases are chosen so that the round trip is an INTEGER number of base samples and the
//     composite response is symmetric about it; latency is reported exactly.
//
// WHAT IT DELIVERS, pinned in CascadeOversamplerTests (both composites recovered through the public API):
// over 14 rates from 8 kHz to 768 kHz and factors 2..64, images and aliases of anything below fs/2 at
// -91.0 dB or lower and one pass within 0.0042 dB of flat up to the band edge; over EVERY taps count the
// rule can produce (110 of them, 16..125), -90.9 dB and 0.0049 dB — the rule's real margin; and a
// palindromic composite.
//
//        rate     stage-1 t   4x latency    firDot MAC per base sample, 4x (PolyphaseOversampler: 512)
//        44.1 k      125         131             572
//        48   k       70          76             348
//        88.2 k       22          28             156
//        96   k       21          27             156
//
// THE PRICE IS LATENCY, not CPU: 131 base samples at 44.1 kHz against 63. That is the width of the
// transition, not the structure — strictness with a 2.05 kHz transition costs about 250 taps at 2 fs
// whatever the topology (Kaiser's length estimate; the rule builds 250). Above 44.1 kHz the transition
// widens and the price falls, but its two halves cross at different rates: the MULTIPLY count drops below
// the fixed-cutoff design's from about 46 kHz (348 at 48 kHz against 512), the LATENCY only from about
// 50 kHz (4x: 76 at 48 kHz against 63; 28 at 88.2 kHz).
//
// WHAT IT IS NOT FOR. The certified true-peak reference (`analysis::ReferenceTruePeakMeter`) is
// PolyphaseOversampler at 4x/32 by contract, and stays so. Factors that are not powers of two are
// refused here and remain PolyphaseOversampler's.
//
// LAW 11(b): unlike PolyphaseOversampler (whose channel CLAMP is P55), prepare() REFUSES a channel count
// outside [1, core::kMaxChannels], a factor that is not a power of two in [2, kMaxFactor], and a rate
// outside [kMinSampleRate, kMaxSampleRate] ([8 kHz, 3 MHz]: the core's floor, P51, and the real-time stages'
// ceiling) — and a refused call touches nothing. A stage that offers this topology inherits the rate window:
// PolyphaseOversampler never looked at the rate, so under Cascade a Saturator or a limiter refuses rates below
// 8 kHz that it accepts under Kaiser. (No tap depends on where the floor is: every rate up to 44.1 kHz gets
// the same geometry.)
// LAW 2 (P56): every inner loop is `core::firDot`. The one sum outside it — the halfband decimator's
// centre tap — adds a value that was halved and STORED a sample earlier, so no contraction can fuse the
// multiply into the add (the tree builds -ffp-contract=on; gcc's `fast` fuses across statements). Halving is
// exact above the subnormal range, so a fusion would move bits only there — the store makes it moot.
// Coefficients are designed with core::det::sin and core::det::mul, so the doubles are the same on every
// row before they are narrowed.
// RT-safe: prepare() allocates; upsample()/downsample() do no alloc/lock/throw and take any n (they walk
// it in kChunk pieces through two scratch buffers sized at prepare()).
class CascadeOversampler
{
public:
    static constexpr int    kMaxFactor      = 64;
    static constexpr int    kMaxStages      = 6;           // log2 (kMaxFactor)
    static constexpr double kBandEdgeHz     = 20000.0;
    static constexpr double kEdgeRate       = 44100.0;     // below this rate the edge is 20/44.1 of fs
    static constexpr double kMinSampleRate  = core::kMinSampleRate;   // P51: the core's one floor (8 kHz)
    static constexpr double kMaxSampleRate  = 3.0e6;
    static constexpr double kBeta           = 9.5;
    static constexpr double kPassA0 = 2.745, kPassA1 = 1.4;    // one pass within 0.005 dB, in fs/t below the cutoff
    static constexpr double kStopB0 = 3.01,  kStopB1 = 1.6;    // |H| <= -90 dB, in fs/t above the cutoff
    static constexpr int    kMinFirstTaps   = 16;
    static constexpr int    kMaxFirstTaps   = 1024;        // a loop bound; the rule never reaches it (125 at most)
    static constexpr int    kChunk          = 64;          // base samples per internal pass
    static constexpr int    kHalfbandM[kMaxStages] = { 0, 6, 5, 5, 3, 3 };   // stage k>=2 is 4m+3 taps

    //==============================================================================
    // THE DESIGN a rate and a factor get — a pure function of the two, so a caller can price a topology
    // (latency above all) without preparing one. FALSE, `out` untouched, exactly where prepare() refuses.
    struct Design
    {
        int    stages = 0;                          // log2 (factor)
        int    factor = 0;
        double bandEdgeHz = 0.0;                    // one pass within 0.005 dB up to here
        int    firstTapsPerPhase = 0;               // stage 1 is 2 * this taps
        double firstCutoffHz = 0.0;                 // stage 1's -6 dB point
        double firstCutoff = 0.0;                   // the same, in cycles per sample of the 2 fs stream
        int    taps[kMaxStages] = {};               // every stage's length (stage 1: 2t; the rest: 4m+3)
        int    decimationPhase[kMaxStages] = {};    // 0 or 1: which of the two input samples a decimator emits at
        int    latencySamples = 0;                  // round trip, base samples, exact
        int    upLegTwice = 0;                      // the UP leg alone, in top-rate samples, times two (it may be a half)
    };

    [[nodiscard]] static bool designFor (double sampleRate, int factor, Design& out) noexcept
    {
        if (! (sampleRate >= kMinSampleRate) || ! (sampleRate <= kMaxSampleRate)) return false;   // NaN fails both
        if (factor < 2 || factor > kMaxFactor || (factor & (factor - 1)) != 0) return false;
        Design d;
        d.factor = factor;
        while ((1 << d.stages) < factor) ++d.stages;
        d.bandEdgeHz = std::min (kBandEdgeHz, core::det::mul (kBandEdgeHz / kEdgeRate, sampleRate));
        const double room = 0.5 * sampleRate - d.bandEdgeHz;
        int t = kMinFirstTaps;
        for (; t < kMaxFirstTaps; ++t)
        {
            const double need = (kPassA0 + kStopB0 + (kPassA1 + kStopB1) / (double) t) * sampleRate / (double) t;
            if (need <= room) break;
        }
        d.firstTapsPerPhase = t;
        d.firstCutoffHz = 0.5 * sampleRate - (kStopB0 + kStopB1 / (double) t) * sampleRate / (double) t;
        d.firstCutoff   = d.firstCutoffHz / (2.0 * sampleRate);
        d.taps[0] = 2 * t;
        for (int k = 1; k < d.stages; ++k) d.taps[k] = 4 * kHalfbandM[k] + 3;

        // Integer latency. The round trip, in top-rate samples, is sum_k R_k * ((N_k - 1) - d_k), where R_k
        // is how many top-rate samples one output sample of stage k spans; the phases are searched in a
        // fixed order and the first that makes it a multiple of the factor wins — deterministic, and the
        // same condition makes the composite response symmetric about that integer (a non-multiple choice
        // can still put a peak on an integer, asymmetrically; the tests check the symmetry, not the peak).
        const int combos = 1 << d.stages;
        bool found = false;
        for (int mask = 0; mask < combos && ! found; ++mask)
        {
            long total = 0;
            int  span  = factor;
            for (int k = 0; k < d.stages; ++k)
            {
                span /= 2;
                total += (long) span * ((long) (d.taps[k] - 1) - (long) ((mask >> k) & 1));
            }
            if (total % factor == 0)
            {
                found = true;
                d.latencySamples = (int) (total / factor);
                long up = 0; span = factor;
                for (int k = 0; k < d.stages; ++k) { span /= 2; up += (long) span * (long) (d.taps[k] - 1); }
                d.upLegTwice = (int) up;
                for (int k = 0; k < d.stages; ++k) d.decimationPhase[k] = (mask >> k) & 1;
            }
        }
        if (! found) return false;                  // not reachable for the lengths above (pinned for every factor)
        out = d;
        return true;
    }

    static int latencyFor (double sampleRate, int factor) noexcept
    {
        Design d;
        return designFor (sampleRate, factor, d) ? d.latencySamples : 0;
    }

    //==============================================================================
    // WHAT prepare() ASKS THE HEAP FOR (law 11d). FALSE, `out` untouched, exactly where prepare() refuses.
    struct Storage
    {
        std::size_t coeffs = 0, rings = 0, scratch = 0;   // floats
        std::size_t cursors = 0;                          // ints
        std::uint64_t bytes() const noexcept
        {
            return (std::uint64_t) sizeof (float) * ((std::uint64_t) coeffs + rings + scratch)
                 + (std::uint64_t) sizeof (int) * (std::uint64_t) cursors;
        }
        bool fitsWithin (const Storage& other) const noexcept
        {
            return coeffs <= other.coeffs && rings <= other.rings && scratch <= other.scratch && cursors <= other.cursors;
        }
    };

    [[nodiscard]] static bool storageFor (double sampleRate, int factor, int maxChannels, Storage& out) noexcept
    {
        if (maxChannels < 1 || maxChannels > core::kMaxChannels) return false;
        Design d;
        if (! designFor (sampleRate, factor, d)) return false;
        Layout lay;
        layoutFor (d, lay);
        Storage st;
        st.coeffs  = lay.coeffs;
        st.rings   = (std::size_t) maxChannels * lay.ringsPerChannel;
        st.scratch = 2u * (std::size_t) kChunk * (std::size_t) (factor / 2);
        st.cursors = (std::size_t) maxChannels * (std::size_t) d.stages * 3u;
        out = st;
        return true;
    }

    //==============================================================================
    bool prepare (double sampleRate, int factor, int maxChannels)
    {
        Design d;
        Storage st;
        if (! designFor (sampleRate, factor, d) || ! storageFor (sampleRate, factor, maxChannels, st)) return false;
        design_ = d;
        layoutFor (d, lay_);
        channels_ = maxChannels;
        coeffs_.assign (st.coeffs, 0.0f);
        designCoefficients();
        rings_.assign (st.rings, 0.0f);
        cursors_.assign (st.cursors, 0);
        scratch_.assign (st.scratch, 0.0f);
        return true;
    }

    void reset() noexcept
    {
        std::fill (rings_.begin(), rings_.end(), 0.0f);
        std::fill (cursors_.begin(), cursors_.end(), 0);
    }

    // One channel's histories and cursors, every other channel bit-exact (the same contract as
    // PolyphaseOversampler::resetChannel). Unprepared or out of range: a no-op.
    void resetChannel (int c) noexcept
    {
        if (c < 0 || c >= channels_) return;
        std::fill_n (rings_.begin() + (std::ptrdiff_t) c * (std::ptrdiff_t) lay_.ringsPerChannel,
                     (std::ptrdiff_t) lay_.ringsPerChannel, 0.0f);
        std::fill_n (cursors_.begin() + (std::ptrdiff_t) c * (std::ptrdiff_t) design_.stages * 3,
                     (std::ptrdiff_t) design_.stages * 3, 0);
    }

    int factor() const noexcept          { return channels_ > 0 ? design_.factor : 0; }
    int latencySamples() const noexcept  { return channels_ > 0 ? design_.latencySamples : 0; }
    const Design& design() const noexcept { return design_; }

    // n base samples per channel -> n*factor. Channels past the prepared count are not touched (as in
    // PolyphaseOversampler; the consumers refuse such a call before it gets here).
    void upsample (const float* const* in, int channels, int n, float* const* out) noexcept
    {
        const int nc = channels < channels_ ? channels : channels_;
        const int S = design_.stages, F = design_.factor;
        for (int c = 0; c < nc; ++c)
            for (int off = 0; off < n; )
            {
                const int m = std::min (kChunk, n - off);
                const float* src = in[c] + off;
                int len = m;
                for (int s = 0; s < S; ++s)
                {
                    float* dst = (s == S - 1) ? out[c] + (std::ptrdiff_t) off * F : scratchAt (s & 1);
                    if (s == 0) kaiserUp (c, src, len, dst); else halfbandUp (s, c, src, len, dst);
                    src = dst;
                    len *= 2;
                }
                off += m;
            }
    }

    // n*factor samples per channel -> n base samples.
    void downsample (const float* const* in, int channels, int n, float* const* out) noexcept
    {
        const int nc = channels < channels_ ? channels : channels_;
        const int S = design_.stages, F = design_.factor;
        for (int c = 0; c < nc; ++c)
            for (int off = 0; off < n; )
            {
                const int m = std::min (kChunk, n - off);
                const float* src = in[c] + (std::ptrdiff_t) off * F;
                int len = m * F;
                for (int s = S - 1; s >= 0; --s)
                {
                    float* dst = (s == 0) ? out[c] + off : scratchAt ((S - 1 - s) & 1);
                    if (s == 0) kaiserDown (c, src, len, dst); else halfbandDown (s, c, src, len, dst);
                    src = dst;
                    len /= 2;
                }
                off += m;
            }
    }

private:
    // Where everything lives. Per stage: the tapped length the inner product runs over (padded to a
    // multiple of four), and the offsets of its coefficients and of its three rings inside one channel's
    // slab. Every ring is double-length and newest-first, as in PolyphaseOversampler.
    struct Layout
    {
        int pad[kMaxStages] = {};                 // stage 1: tppPad (up) — its down ring is nPad, below
        int nPad = 0;                             // stage 1's full length, padded
        std::size_t coefUp[kMaxStages] = {}, coefDown[kMaxStages] = {};
        std::size_t ringUp[kMaxStages] = {}, ringDown[kMaxStages] = {}, ringCtr[kMaxStages] = {};
        std::size_t coeffs = 0, ringsPerChannel = 0;
    };

    static void layoutFor (const Design& d, Layout& lay) noexcept
    {
        lay = Layout {};
        std::size_t co = 0, ro = 0;
        const int t = d.firstTapsPerPhase;
        lay.pad[0] = core::firPadLen (t);
        lay.nPad   = core::firPadLen (2 * t);
        lay.coefUp[0] = co;   co += 2u * (std::size_t) lay.pad[0];         // phase-major, two phases
        lay.coefDown[0] = co; co += (std::size_t) lay.nPad;               // contiguous
        lay.ringUp[0] = ro;   ro += 2u * (std::size_t) lay.pad[0];
        lay.ringDown[0] = ro; ro += 2u * (std::size_t) lay.nPad;
        for (int k = 1; k < d.stages; ++k)
        {
            const int m = kHalfbandM[k];
            lay.pad[k] = core::firPadLen (2 * m + 2);                     // the non-centre taps
            lay.coefUp[k] = co;   co += (std::size_t) lay.pad[k];
            lay.coefDown[k] = co; co += (std::size_t) lay.pad[k];
            lay.ringUp[k] = ro;   ro += 2u * (std::size_t) lay.pad[k];
            lay.ringDown[k] = ro; ro += 2u * (std::size_t) lay.pad[k];
            lay.ringCtr[k] = ro;  ro += 2u * (std::size_t) (m + 1);
        }
        lay.coeffs = co;
        lay.ringsPerChannel = ro;
    }

    float* scratchAt (int which) noexcept
    {
        return scratch_.data() + (std::ptrdiff_t) which * (std::ptrdiff_t) kChunk * (std::ptrdiff_t) (design_.factor / 2);
    }
    float* ring (int c, std::size_t off) noexcept
    {
        return rings_.data() + (std::ptrdiff_t) c * (std::ptrdiff_t) lay_.ringsPerChannel + (std::ptrdiff_t) off;
    }
    int& cursor (int c, int s, int which) noexcept
    {
        return cursors_[(std::size_t) ((c * design_.stages + s) * 3 + which)];
    }

    // Push one sample into a newest-first double-length ring of `len`; returns the window start.
    static float* push (float* h, int& pos, int len, float x) noexcept
    {
        if (pos == 0) pos = len;
        --pos;
        h[pos] = x;
        h[pos + len] = x;
        return h + pos;
    }

    void kaiserUp (int c, const float* x, int n, float* y) noexcept
    {
        const int len = lay_.pad[0];
        float* h = ring (c, lay_.ringUp[0]);
        int& pos = cursor (c, 0, 0);
        const float* p0 = coeffs_.data() + lay_.coefUp[0];
        const float* p1 = p0 + len;
        for (int i = 0; i < n; ++i)
        {
            const float* w = push (h, pos, len, x[i]);
            y[2 * i]     = 2.0f * core::firDot (p0, w, len);        // *2 restores the zero-stuff gain
            y[2 * i + 1] = 2.0f * core::firDot (p1, w, len);
        }
    }

    void kaiserDown (int c, const float* u, int n2, float* y) noexcept
    {
        const int len = lay_.nPad, d = design_.decimationPhase[0];
        float* h = ring (c, lay_.ringDown[0]);
        int& pos = cursor (c, 0, 1);
        const float* proto = coeffs_.data() + lay_.coefDown[0];
        for (int i = 0; 2 * i < n2; ++i)
        {
            const float* w = push (h, pos, len, u[2 * i]);
            if (d == 0) y[i] = core::firDot (proto, w, len);
            w = push (h, pos, len, u[2 * i + 1]);
            if (d == 1) y[i] = core::firDot (proto, w, len);
        }
    }

    // Halfband 2x up: y[2n] = 2 * sum h[2i] x[n-i] (the tapped phase), y[2n+1] = x[n-m] (centre 1/2, times 2).
    void halfbandUp (int s, int c, const float* x, int n, float* y) noexcept
    {
        const int len = lay_.pad[s], m = kHalfbandM[s];
        float* h = ring (c, lay_.ringUp[s]);
        int& pos = cursor (c, s, 0);
        const float* a = coeffs_.data() + lay_.coefUp[s];
        for (int i = 0; i < n; ++i)
        {
            const float* w = push (h, pos, len, x[i]);
            y[2 * i]     = core::firDot (a, w, len);
            y[2 * i + 1] = w[m];
        }
    }

    // Halfband 2x down, decimation phase d: samples of parity d feed the tapped ring; the other parity feeds
    // the centre line, which STORES 0.5*x so the add below never meets a multiply.
    void halfbandDown (int s, int c, const float* u, int n2, float* y) noexcept
    {
        const int len = lay_.pad[s], m = kHalfbandM[s], d = design_.decimationPhase[s];
        float* h = ring (c, lay_.ringDown[s]);
        float* ctr = ring (c, lay_.ringCtr[s]);
        int& pos = cursor (c, s, 1);
        int& cpos = cursor (c, s, 2);
        const float* b = coeffs_.data() + lay_.coefDown[s];
        for (int i = 0; 2 * i < n2; ++i)
        {
            if (d == 1)
            {
                const float* cw = push (ctr, cpos, m + 1, 0.5f * u[2 * i]);
                const float* w  = push (h, pos, len, u[2 * i + 1]);
                const float acc = core::firDot (b, w, len);
                y[i] = acc + cw[m];
            }
            else
            {
                const float* w  = push (h, pos, len, u[2 * i]);
                const float acc = core::firDot (b, w, len);
                y[i] = acc + ctr[cpos + m];
                push (ctr, cpos, m + 1, 0.5f * u[2 * i + 1]);
            }
        }
    }

    void designCoefficients() noexcept
    {
        std::fill (coeffs_.begin(), coeffs_.end(), 0.0f);
        const double i0b = detail::besselI0 (kBeta);
        auto window = [i0b] (int i, int n)
        {
            const double r = (double) (2 * i - (n - 1)) / (double) (n - 1);
            return detail::besselI0 (kBeta * std::sqrt (std::max (0.0, 1.0 - core::det::mul (r, r)))) / i0b;
        };

        // Stage 1: a Kaiser windowed sinc at 2 fs, the same construction as PolyphaseOversampler's, with the
        // cutoff the rule chose. Normalised to unity DC in float, as there.
        {
            const int n = 2 * design_.firstTapsPerPhase, t = design_.firstTapsPerPhase;
            const double fc  = design_.firstCutoff;
            const double cen = (double) (n - 1) * 0.5;
            float* proto = coeffs_.data() + lay_.coefDown[0];
            double sum = 0.0;
            for (int i = 0; i < n; ++i)
            {
                const double x = (double) i - cen;                     // never 0: n is even
                const double v = core::det::sin (2.0 * core::kPi * fc * x) / (core::kPi * x) * window (i, n);
                proto[i] = (float) v;
                sum += v;
            }
            const float inv = (float) (1.0 / sum);
            for (int i = 0; i < lay_.nPad; ++i) proto[i] *= inv;
            float* phase = coeffs_.data() + lay_.coefUp[0];
            for (int p = 0; p < 2; ++p)
                for (int k = 0; k < t; ++k)
                    phase[p * lay_.pad[0] + k] = proto[k * 2 + p];
        }

        // Stages 2..S: halfbands of 4m+3 taps. Centre 1/2 exactly; even offsets structurally zero; odd
        // offsets a windowed sinc scaled to sum to 1/2. sin(pi t / 2) is +-1 for odd t, so no sine at all.
        for (int s = 1; s < design_.stages; ++s)
        {
            const int m = kHalfbandM[s], n = 4 * m + 3, cidx = 2 * m + 1, nz = 2 * m + 2;
            double v[32] = {};                                        // nz <= 2 * 6 + 2
            double odd = 0.0;
            for (int i = 0; i < nz; ++i)
            {
                const int j = 2 * i, at = j > cidx ? j - cidx : cidx - j;   // odd distance from the centre
                const double sgn = ((at - 1) / 2) % 2 == 0 ? 1.0 : -1.0;
                v[i] = sgn / (core::kPi * (double) at) * window (j, n);
                odd += v[i];
            }
            const double scale = 0.5 / odd;
            float* up = coeffs_.data() + lay_.coefUp[s];
            float* dn = coeffs_.data() + lay_.coefDown[s];
            for (int i = 0; i < nz; ++i)
            {
                const float f = (float) (v[i] * scale);
                dn[i] = f;
                up[i] = 2.0f * f;                                     // exact: a power of two
            }
        }
    }

    Design design_ {};
    Layout lay_ {};
    int    channels_ = 0;
    std::vector<float> coeffs_, rings_, scratch_;
    std::vector<int>   cursors_;
};

} // namespace felitronics::oversampling

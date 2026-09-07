// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// Theory-first FALSIFICATION suite for felitronics::core::StreamResampler.
//
// Distinct from StreamResamplerTests.cpp (which pins OrbitCab's golden values): every expectation here
// is DERIVED from the documented contract + first-principles DSP theory, then asserted — never fitted
// to the code.
//
// 🔴 P34 REPLACED THE KERNEL, AND THREE OF THESE DERIVATIONS DIED WITH IT. They are written out below
// rather than deleted, because "this used to be provable and now is not" is the finding, and a reader
// who does not know it will try to restore them.
//
//   * GONE — Catmull-Rom order (Keys a = -1/2). The cubic reproduced constant / linear / QUADRATIC
//     exactly but not the cubic term, so its pointwise error was O(ω^3) and its SNR fell ~18 dB per
//     octave of frequency (SNR ≈ 7560/ω^6, ~111 dB at f = 0.01). A 64-tap Kaiser-windowed sinc has no
//     such law: inside the passband its floor is set by the phase table and float arithmetic, so the
//     SNR is HIGH AND FLAT and only the designed band edge bends it. The suite now asserts flatness —
//     and asserts it in a form the old kernel would FAIL, which is what keeps it a test.
//   * GONE — bit-exact DC. catmull(c,c,c,c,t) collapsed to exactly c in float. The new rows are a
//     partition of unity by construction, normalised in double and stored as float, so a 64-term float
//     dot product lands within a derived 4.6e-6 (measured 3.6e-7) instead of on the nose.
//   * GONE — exact polynomial reproduction. The cubic was an INTERPOLATING kernel (it passed through
//     the input samples); a windowed sinc is an APPROXIMATING one and does not. What survives is
//     agreement to float precision, which the falsify suite pins with derived bounds.
//   * KEPT — Linearity: feed + a fixed weight set are linear in the samples, so
//     resample(a·x + b·y) == a·resample(x) + b·resample(y) up to float rounding only.
//   * KEPT AND UNCHANGED — Count: an output is emitted while floor(kHalf + k·r) + kHalf < N + kTaps,
//     i.e. k < N/r, so K = ⌈N/r⌉ exactly as before: the priming and the availability test moved by the
//     same amount. The produced-count drift stays O(1), independent of N.
//   * MOVED — Latency: at r = 1 the identity path copies, so out[k] = x[k - kHalf]; a general impulse
//     onset sits at k* ≈ (p + kHalf)/inPerOut. Every one of those numbers is read from the class
//     (StreamResampler::kHalf / delayInputSamples()), never restated here — restating it is exactly how
//     the shipped latency formula stayed 2.16 samples wrong through a whole release cycle.
//
// Last group falsifies two KNOWN latent hazards on CONTRACT-VIOLATION inputs (feeds larger than the
// buffer; negative/absurd capacity). On the pre-hardening header these tripped ASan (feed() negative-size
// memmove; a backstop that drove pos below its own history depth → an out-of-bounds tap read in produce).
// The clamp value is the KERNEL's history depth, so it moved from 1 to kBehind = 31 with the taps: a
// hard-coded 1 there would now be a live out-of-bounds read rather than a stale comment.

#include <felitronics/core/StreamResampler.h>

#include "felitronics_test.h"

#include <algorithm>
#include <cmath>
#include <vector>

using felitronics::core::StreamResampler;
using felitronics::test::group;
using felitronics::test::ok;

namespace
{
    constexpr double kPi = 3.14159265358979323846;

    bool anyBad (const std::vector<float>& v)
    { for (float x : v) if (std::isnan (x) || std::isinf (x)) return true; return false; }

    // Feed the whole input as one block (capacity sized to hold it), then produce everything in one call.
    // With a single feed + single produce, output k reads at buf-position pos_k = 1 + k·inPerOut with no
    // intervening compaction, so the analytic mapping below is exact.
    std::vector<float> resampleOneShot (double inRate, double outRate, const std::vector<float>& in)
    {
        StreamResampler r;
        r.reset (inRate, outRate, (int) in.size() + 16);
        r.feed (in.data(), (int) in.size());
        std::vector<float> out ((std::size_t) ((double) in.size() * outRate / inRate) + 64);
        const int k = r.produceAvailable (out.data(), (int) out.size());
        out.resize ((std::size_t) k);
        return out;
    }

    // Stream `in` through in `block`-sized feeds, draining after each; return the whole output.
    std::vector<float> resampleStreaming (double inRate, double outRate, const std::vector<float>& in, int block)
    {
        StreamResampler r;
        r.reset (inRate, outRate, block * 2 + 16);
        std::vector<float> out, tmp ((std::size_t) (block * 4 + 16));
        for (std::size_t i = 0; i < in.size(); i += (std::size_t) block)
        {
            const int n = (int) std::min ((std::size_t) block, in.size() - i);
            r.feed (in.data() + i, n);
            const int got = r.produceAvailable (tmp.data(), (int) tmp.size());
            out.insert (out.end(), tmp.begin(), tmp.begin() + got);
        }
        return out;
    }

    // Pure sine at input rate: f is cycles/sample (normalized frequency).
    std::vector<float> sineNorm (int n, double f, double amp)
    {
        std::vector<float> v ((std::size_t) n);
        for (int i = 0; i < n; ++i) v[(std::size_t) i] = (float) (amp * std::sin (2.0 * kPi * f * i));
        return v;
    }

    // SNR (dB) of the one-shot resample of a normalized-frequency sine vs the analytic interpolation target.
    double sineResampleSnrDb (double inRate, double outRate, double f, int n, double amp)
    {
        const auto in  = sineNorm (n, f, amp);
        const auto out = resampleOneShot (inRate, outRate, in);
        const double ipo = inRate / outRate;
        const int K = (int) out.size();
        double sig = 0.0, noise = 0.0; int c = 0;
        for (int k = 64; k + 64 < K; ++k)                            // skip startup / tail edges
        {
            // out[k] targets input position k·ipo - kHalf. Asked of the class, not restated: this is
            // the number that moved when the kernel did, and hard-coding it is how it went stale before.
            const double truev = amp * std::sin (2.0 * kPi * f * (k * ipo - StreamResampler::delayInputSamples()));
            const double e = (double) out[(std::size_t) k] - truev;
            noise += e * e; sig += truev * truev; ++c;
        }
        return (c > 0 && noise > 0.0) ? 10.0 * std::log10 (sig / noise) : 999.0;
    }
}

int main()
{
    std::printf ("felitronics::core stream-resampler THEORY (falsification) tests\n");

    // ------------------------------------------------------------------------------------------------
    group ("DC — a settled constant resamples to the same constant, to a DERIVED float bound");
    {
        // 🔴 THIS USED TO BE BIT-EXACT AND CANNOT BE ANY MORE. That is a real contract change, not a
        // slackened test. The cubic's weights summed to exactly 1 in float arithmetic at every phase
        // (catmull(c,c,c,c,t) = 0.5·2c), so DC came back bit-identical. This kernel's rows are a
        // partition of unity BY CONSTRUCTION — each row is divided by its own sum inside reset() — but
        // that normalisation happens in double and is stored as float, and the run-time sum is a
        // 64-term float dot product. Bound: 64 terms each rounded at 2^-24, against a coefficient
        // vector whose ABSOLUTE sum is ~1.2 (a windowed sinc has negative lobes) → 64·2^-24·1.2 ≈ 4.6e-6.
        // Measured worst deviation across four ratios × four constants: 3.6e-7, with 46 % of samples
        // still landing bit-exact. Asserted at 2e-6 — inside the derivation, an order above the
        // measurement, and far below any real defect (a mis-normalised row is a per-mille effect).
        const double ratios[][2] = { {48000, 48000}, {44100, 48000}, {96000, 48000}, {48000, 44100} };
        const float  consts[]    = { 1.0f, 0.5f, 0.25f, -0.75f };
        int checked = 0; double worst = 0.0;
        for (auto& rr : ratios)
            for (float C : consts)
            {
                std::vector<float> in (6000, C);
                const auto out = resampleOneShot (rr[0], rr[1], in);
                for (std::size_t k = 512; k + 128 < out.size(); ++k)   // skip the startup ramp and the tail
                {
                    ++checked;
                    worst = std::max (worst, std::fabs ((double) out[k] - (double) C));
                }
            }
        std::printf ("      worst settled-DC deviation over %d samples: %.3e\n", checked, worst);
        ok (checked > 5000 && worst < 2.0e-6,
            "every settled DC sample equals the input constant to the derived 64-tap float bound");
    }

    // ------------------------------------------------------------------------------------------------
    group ("linearity — resample(a·x + b·y) == a·resample(x) + b·resample(y) (FP-tight)");
    {
        const double A = 0.7, B = -0.4;
        const int N = 8000;
        auto x = sineNorm (N, 0.013, 0.5);
        auto y = sineNorm (N, 0.041, 0.5);          // (phase-shifting one is unnecessary; linearity is exact in the samples)
        std::vector<float> z ((std::size_t) N);
        for (int i = 0; i < N; ++i) z[(std::size_t) i] = (float) (A * x[(std::size_t) i] + B * y[(std::size_t) i]);

        const int block = 333;                       // odd block: same produce schedule for all three
        const auto ox = resampleStreaming (44100, 48000, x, block);
        const auto oy = resampleStreaming (44100, 48000, y, block);
        const auto oz = resampleStreaming (44100, 48000, z, block);

        const std::size_t n = std::min ({ ox.size(), oy.size(), oz.size() });
        ok (n > 7000, "all three resamples produced a full-length output");
        float worst = 0.0f;
        for (std::size_t k = 0; k < n; ++k)
            worst = std::max (worst, std::fabs (oz[k] - (float) (A * ox[k] + B * oy[k])));
        ok (worst < 5.0e-6f, "superposition holds within float rounding");   // theory: 0 in reals; measured ~1.8e-7
    }

    // ------------------------------------------------------------------------------------------------
    group ("SNR vs frequency — the cubic's ω^6 LAW IS GONE, and flatness is the new claim");
    {
        // 🔴 THE LAW THIS GROUP USED TO ASSERT NO LONGER HOLDS, and that is the acceptance rather than a
        // regression. Catmull-Rom's error is O(f^3), so its SNR fell ~18 dB per octave of frequency and
        // the suite pinned that slope, plus a "design pin" that the top of the band is low-SNR BY
        // DESIGN. A windowed sinc has no such slope: inside the passband its error floor is set by the
        // phase table and float arithmetic, not by the frequency, so the right assertion is the
        // opposite one — SNR must be HIGH and FLAT right up to the band edge.
        //
        // Measured (44.1 -> 48, against an ideally delayed sine, delay = kHalf input samples):
        //   f = 0.005 .. 0.40 -> 104.4 / 101.2 / 100.9 / 100.7 / 102.4 / 100.6 / 96.8 / 95.8 dB
        //   f = 0.45 (19.8 kHz, into the transition band) -> 85.6 dB
        // The cubic at f = 0.30 read about 55 dB and fell away steeply; here the whole passband sits
        // within 9 dB of itself.
        const double fs[] = { 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4 };
        double lo = 1e9, hi = -1e9;
        for (double f : fs)
        {
            const double s = sineResampleSnrDb (44100, 48000, f, 40000, 0.5);
            std::printf ("      SNR(f=%.3f, %5.0f Hz) = %6.1f dB\n", f, f * 44100.0, s);
            lo = std::min (lo, s); hi = std::max (hi, s);
            ok (s >= 90.0, "SNR at f=" + std::to_string (f) + " clears 90 dB — the cubic could not do "
                           "this above f=0.05, and did not claim to");
        }
        ok (hi - lo < 12.0,
            "…and the SNR spread across SIX octaves is " + std::to_string (hi - lo) + " dB. The cubic's "
            "own law was 18 dB per octave, i.e. ~100 dB across this span: a kernel that still obeyed it "
            "would fail this line, which is exactly what makes the line worth asserting");

        // The band edge is a DESIGNED feature and must still be visible: the cutoff sits at 0.99 of the
        // lower Nyquist, so f = 0.45 is inside the transition and must read measurably worse than the
        // passband without collapsing. A kernel with no band edge at all would pass the flatness line
        // above and fail this one.
        const double edge = sineResampleSnrDb (44100, 48000, 0.45, 40000, 0.5);
        std::printf ("      SNR(f=0.450, 19845 Hz, into the transition) = %.1f dB\n", edge);
        ok (edge > 70.0 && edge < lo, "the designed band edge is visible at f=0.45 and does not collapse");
    }

    // ------------------------------------------------------------------------------------------------
    group ("count conservation — produced count tracks n_in/ratio with O(1) drift (no accumulator drift)");
    {
        struct RB { double in, out; };
        const RB ratios[] = { {44100,48000}, {48000,44100}, {96000,48000}, {48000,96000}, {44100,96000} };
        int worstDrift = 0; bool allBounded = true;
        for (auto rb : ratios)
            for (int N : { 10000, 500000 })                          // 50× length: drift bound must NOT grow
            {
                std::vector<float> in ((std::size_t) N, 0.3f);
                const auto out = resampleStreaming (rb.in, rb.out, in, 512);
                const double expected = (double) N * rb.out / rb.in;
                const int drift = (int) std::llround (std::fabs ((double) out.size() - expected));
                worstDrift = std::max (worstDrift, drift);
                if (drift > 2) allBounded = false;                   // derived <1; assert <=2 with margin
            }
        std::printf ("      worst produced-count drift across all ratios/lengths: %d samples\n", worstDrift);
        ok (allBounded, "|K - n_in·outRate/inRate| stays O(1) (<=2) at every ratio, unchanged at 50× length");
    }

    // ------------------------------------------------------------------------------------------------
    group ("latency — exactly kHalf samples at ratio 1; onset shift tracks (p+kHalf)/inPerOut");
    {
        // The whole point of asking the class rather than restating it: D comes from the header, so a
        // change to kTaps moves this test with the code instead of leaving it pinning a stale integer.
        // That is precisely how the old `ceil(3·sr/48000)+3` formula survived being 2.16 samples wrong.
        const int D = StreamResampler::kHalf;

        // Ratio 1: unit impulse at input 0 → the identity path copies, so out[D] is bit-exactly 1.
        {
            std::vector<float> in (256, 0.0f); in[0] = 1.0f;
            const auto out = resampleOneShot (48000, 48000, in);
            bool leadingSilent = out.size() > (std::size_t) D;
            for (int k = 0; k < D && leadingSilent; ++k) leadingSilent = (out[(std::size_t) k] == 0.0f);
            ok (leadingSilent && out[(std::size_t) D] == 1.0f,
                "identity ratio delays the impulse by exactly kHalf samples, bit-exact");
        }
        // Other ratios: an impulse at input p peaks at output k* ≈ (p+D)/inPerOut. At the two integer
        // ratios the peak lands EXACTLY there (measured delta 0.00); at 44.1↔48 the true maximum falls
        // between output samples, so the nearest one is up to half a sample off (measured 0.33 / −0.27).
        {
            struct RB { double in, out; };
            const RB ratios[] = { {48000,96000}, {96000,48000}, {44100,48000}, {48000,44100} };
            const int p = 100;
            bool allAligned = true;
            for (auto rb : ratios)
            {
                std::vector<float> in (900, 0.0f); in[(std::size_t) p] = 1.0f;
                const auto out = resampleOneShot (rb.in, rb.out, in);
                int pk = 0;
                for (std::size_t k = 1; k < out.size(); ++k)
                    if (std::fabs (out[k]) > std::fabs (out[(std::size_t) pk])) pk = (int) k;
                const double ipo = rb.in / rb.out;
                if (std::fabs ((double) pk - ((double) p + D) / ipo) > 1.0) allAligned = false;
            }
            ok (allAligned, "impulse onset sits at (p+kHalf)/inPerOut (±1) — latency scales with the ratio");
        }
    }

    // ------------------------------------------------------------------------------------------------
    group ("backstop hardening — contract-violation feeds stay in-bounds (was UB; ASan/UBSan-verified)");
    {
        // (a) A single feed block LARGER than the whole buffer. Pre-fix: drop = len+n-cap > len, so the
        //     memmove size (size_t)(len-drop) underflowed → ASan negative-size-param. Now: keep the newest
        //     `cap` samples, drop history; everything stays in-bounds.
        {
            StreamResampler r; r.reset (48000, 48000, 64);
            std::vector<float> in (4096, 0.2f);
            r.feed (in.data(), 4096);
            ok (r.len >= 0 && r.len <= (int) r.buf.size() && r.pos >= (double) StreamResampler::kBehind,
                "oversized feed leaves len/pos in-bounds");
            std::vector<float> out (256);
            const int k = r.produceAvailable (out.data(), 256);
            out.resize ((std::size_t) std::max (0, k));
            ok (! anyBad (out), "produce after an oversized feed emits no NaN/Inf and never reads OOB");
        }
        // (a2) A backstop that would drive the read head past its own history. Pre-fix: pos -= drop pushed
        //      pos below 1 (here to -3) → produce read buf[i-1] with i<0 (UBSan unsigned-offset-overflow).
        //      Now: pos is clamped to >= 1, so buf[i-1] never dips below buf[0].
        {
            StreamResampler r; r.reset (48000, 48000, 16);
            std::vector<float> a (20, 0.5f), b (5, -0.5f);
            r.feed (a.data(), 20);                               // fills near capacity, pos still 1
            r.feed (b.data(), 5);                                // triggers the backstop
            ok (r.pos >= (double) StreamResampler::kBehind,
                "backstop never drops pos below kBehind (the lowest tap stays >= buf[0])");
            std::vector<float> out (64);
            const int k = r.produceAvailable (out.data(), 64);
            out.resize ((std::size_t) std::max (0, k));
            ok (! anyBad (out), "produce after a pos-clamping backstop stays in-bounds (no OOB read)");
        }
        // (b) Absurd capacity. Pre-fix: capacity=-1 wrapped (size_t)(-1)+8 → a 7-element buffer, and a normal
        //     feed then underflowed the memmove. Now: negative clamps to 0 → an 8-element buffer, and the
        //     oversized-block backstop keeps it in-bounds. A large negative must not attempt a wrapped huge alloc.
        {
            StreamResampler r; r.reset (48000, 48000, -1);
            ok (r.buf.size() >= (std::size_t) StreamResampler::kTaps,
                "negative capacity is clamped, not wrapped into a tiny/huge buffer — and the floor is now "
                "the kernel's own window, because a buffer shorter than kTaps could never produce at all");
            std::vector<float> in (512, 0.3f);
            r.feed (in.data(), 512);                             // n > cap → newest-cap backstop, in-bounds
            ok (r.len >= 0 && r.len <= (int) r.buf.size() && r.pos >= (double) StreamResampler::kBehind,
                "feed after clamped-negative reset stays in-bounds");

            StreamResampler r2; r2.reset (48000, 48000, -1000000);   // would wrap to a huge alloc pre-fix
            ok (r2.buf.size() >= (std::size_t) StreamResampler::kTaps && r2.buf.size() < 128,
                "large-negative capacity clamps to a small buffer (kTaps + 8 = 72), no wrapped alloc");
        }
        // NOTE: the INT_MAX-side truncation guard (capacity near 2^31 → (int) buf.size() wrap) is likewise fixed
        // in reset(), but exercising it would demand an ~8 GB allocation, so it is not driven at runtime here.
    }

    return felitronics::test::report();
}

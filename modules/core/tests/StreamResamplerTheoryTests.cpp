// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// Theory-first FALSIFICATION suite for felitronics::core::StreamResampler.
//
// Distinct from StreamResamplerTests.cpp (which pins OrbitCab's golden values): every expectation here
// is DERIVED from the documented contract + first-principles DSP theory, then asserted — never fitted to
// the code. The derivations:
//
//   * Catmull-Rom order (Keys a=-1/2). Expanding the four weights against e^{iωj} shows the scheme
//     reproduces constant / linear / QUADRATIC exactly (S0=1, S1=t, S2=t^2) but S3 = t-3t^2+3t^3 != t^3,
//     so the leading pointwise error is O(ω^3) (NOT O(f^4)) with shape P(t)=t(1-t)(1-2t). ∫_0^1 P^2 = 1/210.
//     For a unit-amplitude sine of angular freq ω=2πf per input sample, error power ≈ (ω^3/6)^2·(1/210)·(1/2),
//     signal power 1/2, so  SNR ≈ 7560 / ω^6 . At f=0.01 that is ~111 dB; the SNR falls ~18 dB per octave
//     of f (the ω^6 law). We assert a LOWER bound with margin at low f, and only PIN the high-f value
//     (the header documents Catmull-Rom as low-SNR-by-design — a high-f case can characterize, not fail).
//   * DC: catmull(c,c,c,c,t) collapses to 0.5·(2c)=c; for c whose 5·c is exact (powers-of-two-ish) it is
//     bit-exact, so a settled constant must resample to the EXACT same constant at any ratio.
//   * Linearity: feed + catmull are linear in the samples, so resample(a·x+b·y) == a·resample(x)+b·resample(y)
//     up to float rounding only (measured worst ~1.8e-7).
//   * Count: an output is emitted while floor(1+k·r) <= N (r=inPerOut), i.e. k < N/r, so K=⌈N/r⌉ and the
//     produced-count drift |K - N·outRate/inRate| < 1 — O(1), independent of N (the phase accumulator does
//     not drift). We prove the bound holds unchanged at 50× the length.
//   * Latency: at r=1, out[k]=buf[k+1]=x[k-2] → exactly 2 samples; a general impulse onset sits at
//     k* ≈ (p+2)/inPerOut.
//
// Last group falsifies two KNOWN latent hazards on CONTRACT-VIOLATION inputs (feeds larger than the buffer;
// negative/absurd capacity). On the pre-hardening header these tripped ASan (feed() negative-size-param
// memmove; a backstop that drove pos below 1 → an out-of-bounds buf[i-1] read in produce). They were fixed
// ADDITIVELY (guards on the backstop path only; the valid path stays bit-identical — the golden suite proves
// it). These cases run clean under the ASan/UBSan gate.

#include <felitronics/core/StreamResampler.h>

#include "felitronics_test.h"

#include <algorithm>
#include <cmath>
#include <vector>

using felitronics::core::StreamResampler;
using felitronics::test::approx;
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
            const double truev = amp * std::sin (2.0 * kPi * f * (k * ipo - 2.0));   // out[k] targets input-pos k·ipo - 2
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
    group ("DC exactness — a settled constant resamples to the EXACT same constant, any ratio");
    {
        // catmull(c,c,c,c,t) = 0.5·(2c) = c bit-exactly for constants whose 5·c is exact.
        const double ratios[][2] = { {48000, 48000}, {44100, 48000}, {96000, 48000}, {48000, 44100} };
        const float  consts[]    = { 1.0f, 0.5f, 0.25f, -0.75f };
        int checked = 0, exact = 0;
        for (auto& rr : ratios)
            for (float C : consts)
            {
                std::vector<float> in (6000, C);
                const auto out = resampleOneShot (rr[0], rr[1], in);
                for (std::size_t k = 512; k < out.size(); ++k)       // skip the startup ramp through the 3 leading zeros
                {
                    ++checked;
                    if (out[k] == C) ++exact;                        // BIT-exact, no tolerance
                }
            }
        ok (checked > 5000 && exact == checked, "every settled DC sample equals the input constant exactly");
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
    group ("cubic SNR — low f meets the O(f^3)-error (SNR≈7560/ω^6) bound; high f is a pin");
    {
        const double snr_lo   = sineResampleSnrDb (44100, 48000, 0.01,  20000, 0.5);   // theory ~111 dB
        const double snr_2f   = sineResampleSnrDb (44100, 48000, 0.02,  20000, 0.5);   // theory ~93 dB (-18 dB/oct)
        const double snr_half = sineResampleSnrDb (44100, 48000, 0.005, 40000, 0.5);   // theory ~129 dB (+18 dB/oct)
        std::printf ("      SNR(f=0.01)=%.1f dB  SNR(f=0.02)=%.1f  SNR(f=0.005)=%.1f\n", snr_lo, snr_2f, snr_half);

        ok (snr_lo   >= 100.0, "low-f (f=0.01) resample SNR meets the derived >100 dB floor");
        ok (snr_half >= 118.0, "halving f gains ~18 dB (the ω^6 SNR law), stays above 118 dB");
        // The ω^6 law predicts a ~18 dB step per octave of f; accept a generous band around it.
        const double dropPerOct = snr_lo - snr_2f;
        approx (dropPerOct, 18.0, 8.0, "SNR degrades ~18 dB/octave of f, as O(f^3) error → O(f^6) power predicts");

        // High-f PIN (characterization, not a spec): Catmull-Rom is low-SNR-by-design up top. It must stay
        // finite/positive and drop far below the low-f figure — proving the images are NOT brickwalled.
        const double snr_hi = sineResampleSnrDb (44100, 48000, 0.3, 20000, 0.5);
        std::printf ("      PIN SNR(f=0.30)=%.1f dB (low-SNR-by-design)\n", snr_hi);
        ok (snr_hi > 0.0 && snr_hi < snr_lo - 40.0, "high-f SNR is finite/positive yet far below low-f (design pin)");
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
    group ("latency — exactly 2 samples at ratio 1; onset shift tracks (p+2)/inPerOut");
    {
        // Ratio 1: unit impulse at input 0 → out[0]=out[1]=0, out[2]=1 (bit-exact; catmull @ t=0 = center tap).
        {
            std::vector<float> in (64, 0.0f); in[0] = 1.0f;
            const auto out = resampleOneShot (48000, 48000, in);
            ok (out.size() >= 3 && out[0] == 0.0f && out[1] == 0.0f && out[2] == 1.0f,
                "identity ratio delays the impulse by exactly 2 samples, bit-exact");
        }
        // Other ratios: an impulse at input p peaks at output k* ≈ (p+2)/inPerOut.
        {
            struct RB { double in, out; };
            const RB ratios[] = { {48000,96000}, {96000,48000}, {44100,48000}, {48000,44100} };
            const int p = 50;
            bool allAligned = true;
            for (auto rb : ratios)
            {
                std::vector<float> in (600, 0.0f); in[(std::size_t) p] = 1.0f;
                const auto out = resampleOneShot (rb.in, rb.out, in);
                int pk = 0;
                for (std::size_t k = 1; k < out.size(); ++k)
                    if (std::fabs (out[k]) > std::fabs (out[(std::size_t) pk])) pk = (int) k;
                const double ipo = rb.in / rb.out;
                if (std::fabs ((double) pk - (p + 2) / ipo) > 1.0) allAligned = false;
            }
            ok (allAligned, "impulse onset sits at (p+2)/inPerOut (±1) — latency scales with the ratio");
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
            ok (r.len >= 0 && r.len <= (int) r.buf.size() && r.pos >= 1.0, "oversized feed leaves len/pos in-bounds");
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
            ok (r.pos >= 1.0, "backstop never drops pos below 1 (read head stays >= buf[0])");
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
            ok (r.buf.size() >= 8, "negative capacity is clamped, not wrapped into a tiny/huge buffer");
            std::vector<float> in (512, 0.3f);
            r.feed (in.data(), 512);                             // n > cap → newest-cap backstop, in-bounds
            ok (r.len >= 0 && r.len <= (int) r.buf.size() && r.pos >= 1.0, "feed after clamped-negative reset stays in-bounds");

            StreamResampler r2; r2.reset (48000, 48000, -1000000);   // would wrap to a huge alloc pre-fix
            ok (r2.buf.size() >= 8 && r2.buf.size() < 64, "large-negative capacity clamps to a small buffer, no wrapped alloc");
        }
        // NOTE: the INT_MAX-side truncation guard (capacity near 2^31 → (int) buf.size() wrap) is likewise fixed
        // in reset(), but exercising it would demand an ~8 GB allocation, so it is not driven at runtime here.
    }

    return felitronics::test::report();
}

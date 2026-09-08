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
#include <cstring>
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
        // 🔴 THE SKIP IS DERIVED FROM THE KERNEL, not a round number. Output k is centred on input
        // position k·ipo − D and the aperture reaches kHalf either side of that centre, so the window
        // is not yet clear of reset()'s leading zeros until k·ipo − D ≥ kHalf, i.e. k ≥ (D+kHalf)/ipo.
        // A fixed 64 was fine for the cubic (D = 2, aperture ±2) and is NOT for this kernel: at
        // ipo = 0.5 the first analysed output still had half its window in the priming zeros, and the
        // SNR it reported was the startup transient, not the passband. That cost a round.
        const int edge = (int) std::ceil ((StreamResampler::delayInputSamples() + StreamResampler::kHalf) / ipo) + 8;
        for (int k = edge; k + edge < K; ++k)
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
        // 64-term float dot product.
        //
        // 🔴 THE DERIVATION PUBLISHED HERE FIRST WAS WRONG AND A CREW ROUND CAUGHT IT. It used an
        // absolute coefficient sum of "~1.2", which is what a naive look at a unit-sum kernel suggests.
        // MEASURED on the shipped table: max Σ|w| = 2.77365 on the interpolating leg and 2.26297 on the
        // decimating one — a windowed sinc's negative lobes carry far more weight than its unit sum
        // admits. So the worst-case bound is
        //     kTaps · 2^-24 · max Σ|w| = 64 · 5.96e-8 · 2.774 = 1.06e-5,
        // not 4.6e-6, and the 2e-6 that used to be asserted here was BELOW its own stated derivation —
        // i.e. fitted to the measurement while claiming to be derived. Asserted at 1.2e-5 now: just
        // above the adversarial bound, so no correct implementation can fail it on any toolchain, and
        // still 100x below a real defect (a mis-normalised row is a per-mille effect). The measurement
        // is PRINTED so a regression is visible in the log even while it passes: 4.768e-7 today, and a
        // crew round confirmed the same value bit-for-bit on Apple clang, gcc 14 and MSVC.
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
        ok (checked > 5000 && worst < 1.2e-5,
            "every settled DC sample equals the input constant to the derived 64-tap float bound "
            "(kTaps·2^-24·max sum|w| = 1.06e-5; measured " + std::to_string (worst) + ")");
    }

    // ------------------------------------------------------------------------------------------------
    group ("THE TABLE AS AN OBJECT — invariants no dB measurement in this repo can see");
    {
        // 🔴 WHY THIS GROUP EXISTS, and it is the most expensive lesson of the P34 crew rounds. A
        // diverse-testing round mutated reset() to early-return when `tab` was already populated — a
        // STALE TABLE, i.e. a resampler that keeps designing for the previous ratio forever. That
        // mutant passed **345 of 345 checks** across every suite in this repository. Not one number
        // moved, because every fixture in every suite constructs a FRESH instance, and a fresh instance
        // has an empty table and therefore builds the right one. The production path nobody was
        // exercising is the ordinary one: a DAW changes sample rate, NamBackend::prepare re-enters
        // configureRates, and the SAME down/up objects are reset() to a new ratio.
        //
        // The lesson generalises past this defect: a kernel is a data structure as well as a filter,
        // and its structural invariants have to be asserted directly. Everything below is a property
        // of the TABLE, checkable in microseconds, and each one is a defect class that costs 90 dB or
        // more of transparency while every passband and stopband pin in the tree stays green.

        // (a) reset() to a different ratio must REBUILD, not reuse.
        {
            StreamResampler r;
            r.reset (44100.0, 48000.0, 512);
            const std::vector<float> up = r.tab;
            r.reset (48000.0, 44100.0, 512);
            // Precondition first, so this cannot pass by both tables being empty or both identical for
            // an unrelated reason: the two ratios design genuinely different kernels (cutoff 0.99 of
            // the input Nyquist one way, 0.99·44100/48000 = 0.909 the other).
            ok (! up.empty() && up.size() == r.tab.size(),
                "precondition: both ratios build a table of the same shape");
            bool differs = false;
            for (std::size_t i = 0; i < up.size() && ! differs; ++i) differs = (up[i] != r.tab[i]);
            ok (differs, "reset() to a DIFFERENT ratio rebuilds the phase table — a stale table passes "
                         "every dB assertion in this repository, so it has to be caught structurally");
        }

        // (b) …and the rendered consequence, which is the shape the DAW actually takes: a stage reset
        // to ratio A, RUN, then reset to ratio B must produce exactly what a virgin B produces.
        {
            auto render = [] (bool viaOtherRatio)
            {
                StreamResampler r;
                if (viaOtherRatio)
                {
                    r.reset (48000.0, 44100.0, 2048);
                    std::vector<float> warm (1024, 0.3f), sink (4096);
                    r.feed (warm.data(), 1024);
                    (void) r.produceAvailable (sink.data(), 4096);   // make it carry real state
                }
                r.reset (44100.0, 48000.0, 2048);
                std::vector<float> in (1024), out (4096);
                for (int i = 0; i < 1024; ++i) in[(std::size_t) i] = 0.5f * (float) std::sin (0.21 * i);
                r.feed (in.data(), 1024);
                const int k = r.produceAvailable (out.data(), 4096);
                out.resize ((std::size_t) k);
                return out;
            };
            const auto viaOther = render (true), virginRun = render (false);
            bool identical = (viaOther.size() == virginRun.size());
            for (std::size_t i = 0; identical && i < viaOther.size(); ++i)
                identical = (std::memcmp (&viaOther[i], &virginRun[i], sizeof (float)) == 0);
            ok (identical && viaOther.size() > 900,
                "a resampler reset to one ratio, RUN, and then reset to another renders bit-identically "
                "to one that only ever saw the second — reset() is a full re-design, not a top-up");
        }

        // (c) THE WINDOW EDGE. kernelAt is defined to return exactly zero at |x| >= kHalf, and the seam
        // row depends on it: a `>` instead of `>=` there (the header names it "the classic off-by-one")
        // leaves a ~1.3e-5 tap alive at the edge, which survives the whole suite except two float
        // tolerances that catch it by luck. Asserted directly, it cannot hide.
        {
            const double fc = StreamResampler::kCutoff;
            ok (StreamResampler::kernelAt ((double) StreamResampler::kHalf, fc) == 0.0
                && StreamResampler::kernelAt (-(double) StreamResampler::kHalf, fc) == 0.0,
                "the window is EXACTLY zero at |x| = kHalf, both signs");
            ok (StreamResampler::kernelAt ((double) StreamResampler::kHalf + 1.0, fc) == 0.0,
                "…and beyond it");
            ok (std::fabs (StreamResampler::kernelAt (0.0, fc) - fc) < 1e-12,
                "…and the centre tap is fc, i.e. the window is 1 there and the sinc is unwindowed at 0");
        }

        // (d) THE SEAM ROW. Row kPhases must be row 0 shifted by exactly one tap; that is the whole
        // reason a (kPhases+1)-row table exists instead of a wrapping index. If it is not, one output
        // in every 147 at the shipped ratio uses a kernel misaligned by a full input sample.
        {
            StreamResampler r;
            r.reset (44100.0, 48000.0, 512);
            const float* row0   = r.tab.data();
            const float* rowEnd = r.tab.data() + (std::size_t) StreamResampler::kPhases * StreamResampler::kTaps;
            bool shifted = true;
            for (int j = 1; j < StreamResampler::kTaps; ++j) shifted = shifted && (rowEnd[j] == row0[j - 1]);
            ok (shifted, "row kPhases IS row 0 shifted one tap — the seam is consistent, so linear "
                         "interpolation never has to wrap");
            ok (rowEnd[0] == 0.0f && row0[StreamResampler::kTaps - 1] == 0.0f,
                "…and the taps that enter and leave at the seam are exactly zero, which is what makes "
                "the shift exact rather than approximate");
        }

        // (e) EVERY row is a partition of unity, checked as data rather than through a DC signal.
        {
            StreamResampler r;
            r.reset (48000.0, 44100.0, 512);
            double worst = 0.0;
            for (int p = 0; p <= StreamResampler::kPhases; ++p)
            {
                double sum = 0.0;
                for (int j = 0; j < StreamResampler::kTaps; ++j)
                    sum += (double) r.tab[(std::size_t) p * StreamResampler::kTaps + (std::size_t) j];
                worst = std::max (worst, std::fabs (sum - 1.0));
            }
            std::printf ("      worst |sum(row) - 1| over all %d rows: %.3e\n", StreamResampler::kPhases + 1, worst);
            ok (worst < 1.0e-6, "every phase row sums to 1 to float storage precision");
        }

        // (f) IDENTITY builds no table at all — a bit-copy has nothing to read, and 128 KiB per
        // instance is not free when a NamStage holds four of them and retires up to 64 backends.
        {
            StreamResampler r;
            r.reset (48000.0, 48000.0, 512);
            ok (r.identity && r.tab.empty(), "an identity ratio designs no kernel and allocates no table");
            r.reset (48000.0, 44100.0, 512);
            ok (! r.identity && ! r.tab.empty(), "…and moving off identity builds one");
        }
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
        // 🔴 TWO RATIOS, and the second one is not decoration. A crew mutation that corrupted the taps
        // with a ZERO-SUM symmetric perturbation — invisible to DC, invisible to the 44.1↔48 rows —
        // cost 4.4 dB at 3 kHz through a 2x INTERPOLATION and survived every suite, because nothing
        // measured an integer up-ratio's passband at all. 48 -> 96 is that case: the phase alternates
        // between exactly two rows, so a defect living in one row has nowhere to average out.
        struct RB { double in, out; const char* name; };
        const double fs[] = { 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4 };
        double lo = 1e9, hi = -1e9;
        for (const RB rb : { RB {44100, 48000, "44.1->48"}, RB {48000, 96000, "48->96 (integer 2x)"} })
            for (double f : fs)
            {
                const double s = sineResampleSnrDb (rb.in, rb.out, f, 40000, 0.5);
                std::printf ("      %-20s SNR(f=%.3f, %5.0f Hz) = %6.1f dB\n", rb.name, f, f * rb.in, s);
                lo = std::min (lo, s); hi = std::max (hi, s);
                // 🔴 THE FLOOR IS DERIVED, and the first version of this line was not — it was read off
                // the 44.1->48 column alone (min 95.8) and set at 92, which the 48->96 column then
                // failed at three frequencies. The failure was correct and the threshold was wrong.
                // What bounds this SNR is the window's own PASSBAND RIPPLE: an ideal-sine oracle counts
                // a systematic gain error of delta as "noise", so SNR <= -20*log10(delta), and for a
                // Kaiser beta = 8.6, delta = 10^(-86.7/20) = 4.6e-5 -> 86.7 dB. That is a property of
                // the DESIGN, identical at every ratio; the 44.1->48 column simply happens to sample
                // the ripple pattern at kinder points. Asserted at 84 dB: below the derived floor, so
                // no correct kernel can fail it anywhere, and 40 dB above what the cubic reached in
                // this band (its SNR at f=0.3 was about 55 dB and falling 18 dB per octave).
                ok (s >= 84.0, std::string (rb.name) + " SNR at f=" + std::to_string (f)
                               + " clears the derived 86.7 dB Kaiser ripple floor (asserted at 84)");
            }
        // The spread is bounded by the same ripple: the ceiling is float/table precision (~105 dB) and
        // the floor is the 86.7 dB ripple, so ~18 dB is the designed range and 25 dB is the assertion.
        // A kernel still obeying the cubic's 18 dB-per-octave law would spread ~100 dB across this span
        // and fail by a factor of four, which is what makes a loose-looking bound still worth having.
        ok (hi - lo < 25.0,
            "…and the SNR spread across SIX octaves and two ratios is " + std::to_string (hi - lo)
            + " dB, inside the range the window's own ripple allows — the cubic's law would give ~100");

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

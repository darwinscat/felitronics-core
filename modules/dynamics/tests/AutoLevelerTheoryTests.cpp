// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// THEORY-FIRST FALSIFICATION suite for felitronics::dynamics::AutoLeveler.
//
// This battery is derived from the header's DOCUMENTED contract + energy-matching theory FIRST,
// then run against the code — every tolerance below is derived from the header's own constants
// (τ=150 ms, the 9 / 40 / 20 dB/s ceilings, the −18 dBFS seed, the −60 dBFS floor, the ±clamp),
// NOT copied from the sibling ported tests. It exists to BREAK the gain dynamics, not to confirm
// them. Never loosen a check to go green: a fingerprint-pinned surprise gets pinned + reported,
// a contract-path defect gets an additive fix — the tolerances here are load-bearing.
//
// Attack surface (one group each):
//   1. steady-state convergence to sqrt(dryMS/mixMS) across ±40 dB, both directions + clamp;
//   2. the per-sample dB/s ceiling as a HARD invariant under thousands of hostile segments
//      (random level jumps, enable flips, transition windows, block sizes 1..2048, 44.1..192 kHz);
//   3. block-size invariance of the realized slew RATE (the frozen-leveler bug this design replaced);
//   4. sample-rate invariance of that rate;
//   5. seed() semantics: instant applied-gain seed, clamp, the (−18 dBFS, ·/g²) fixed point,
//      and the −18 dBFS "prompt handover" (a 0 dBFS seed would out-mass a real signal);
//   6. pump-proofing: a tremolo'd wet must NOT drive gain ripple past what the 150 ms follower admits;
//   7. silence handling: sub-floor energy HOLDS (never slams to the +36 dB clamp);
//   8. determinism + reset() bit-identity.

#include <felitronics_test.h>
#include <felitronics/dynamics/AutoLeveler.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <random>
#include <string>
#include <vector>

using felitronics::dynamics::AutoLeveler;
namespace test = felitronics::test;

namespace
{
    //--------------------------------------------------------------- documented contract constants
    // Mirrored from AutoLeveler.h's private section — the CONTRACT the suite derives every bound from.
    constexpr double kTau        = 0.15;    // s  — follower time constant
    constexpr double kSeedMeanSq = 0.0158;  //    — −18 dBFS² seed energy
    constexpr double kSlew       = 9.0;     // dB/s — normal makeup ceiling
    constexpr double kSnapSlew   = 40.0;    // dB/s — deterministic-retarget / fast ceiling
    constexpr double kSnapWindow = 0.35;    // s  — fast-rate budget per snap
    constexpr double kFloorMeanSq= 1.0e-6;  //    — ~−60 dBFS: below this the target HOLDS
    constexpr float  kMinGain    = 0.0631f; // −24 dB
    constexpr float  kMaxGain    = 63.10f;  // +36 dB
    constexpr float  kSilentDb   = -120.0f;

    // juce::Decibels::gainToDecibels (g, -120.0f), replicated exactly (== what the class returns).
    float gainDb (float g) { return g > 0.0f ? std::max (kSilentDb, std::log10 (g) * 20.0f) : kSilentDb; }

    // dB delta between two positive gains, evaluated in DOUBLE — avoids compounding the float-log10
    // noise of differencing two float-dB readings when we assert the per-sample ceiling.
    double dbDelta (float g1, float g0) { return 20.0 * std::log10 ((double) g1 / (double) g0); }

    float clampGain (float g) { return std::clamp (g, kMinGain, kMaxGain); }

    // The applied glide is a GEOMETRIC progression by a per-sample ratio stored as a FLOAT
    // (prepare(): ratioUpNormal = (float) pow (10, (kSlew/sr)/20)). Quantising that ratio to float
    // makes the realised rate a deterministic function of sr — 9 dB/s ± ~0.6% (e.g. +0.51% @ 96 kHz,
    // −0.60% @ 192 kHz). This reproduces the SAME float ratio the code stores, so the rate tests can
    // pin the exact shipped number instead of a nominal 9 that the float ramp never quite equals.
    double predictedRate (double sr)
    {
        const float ratio = (float) std::pow (10.0, (kSlew / sr) / 20.0);
        return 20.0 * std::log10 ((double) ratio) * sr;
    }

    // Drive a STATIONARY (dry, mix) scenario from a fresh prepare() for `seconds` and return the
    // final applied gain in dB. From reset the follower ratio is preserved immediately (both climb
    // from 0), so convergence is limited only by the slew ramp |Δdb|/kSlew plus a few τ.
    float settleDb (double dryMS, double mixMS, int blk, double sr, double seconds)
    {
        AutoLeveler a; a.prepare (sr);
        const long blocks = std::lround (seconds * sr / blk);
        for (long i = 0; i < blocks; ++i)
        {
            a.processBlock (dryMS, mixMS, true, blk);
            for (int s = 0; s < blk; ++s) a.getNextGain();
        }
        return a.currentGainDb();
    }

    // Per-sample applied-gain trajectory (dB) for a STATIONARY (dry, mix) scenario, chopped by a
    // repeating block-size pattern. `nSamples` samples captured. Used to break block-size invariance.
    std::vector<float> perSampleDb (double dryMS, double mixMS,
                                    const std::vector<int>& blkPattern, double sr, long nSamples)
    {
        AutoLeveler a; a.prepare (sr);
        std::vector<float> out; out.reserve ((size_t) nSamples);
        size_t p = 0;
        while ((long) out.size() < nSamples)
        {
            const int blk = blkPattern[p++ % blkPattern.size()];
            a.processBlock (dryMS, mixMS, true, blk);
            for (int s = 0; s < blk && (long) out.size() < nSamples; ++s)
            {
                a.getNextGain();
                out.push_back (a.currentGainDb());
            }
        }
        return out;
    }
}

int main()
{
    std::printf ("felitronics::dynamics AutoLeveler THEORY-first falsification suite\n");

    const std::vector<double> rates { 44100.0, 48000.0, 88200.0, 96000.0, 176400.0, 192000.0 };

    //================================================================ 1. steady-state convergence
    // Stationary dry at D, wet at M => applied gain must converge to sqrt(D/M). Test both directions
    // across the full non-clamped span, plus saturation beyond the ±clamp. Tolerance: run 6 s
    // (>> |Δdb|/9 + 5τ for every case below), the follower residual is ~0 (ratio preserved from
    // reset), so only the float dB round-trip remains => 0.01 dB. Beyond the clamp: land on the rail.
    test::group ("steady-state: applied gain -> sqrt(dryMS/mixMS), both directions, ±40 dB");
    {
        const double D = 1.0;
        for (double wantDb : { -40.0, -36.0, -24.0, -18.0, -12.0, -6.0, -3.0,
                                3.0, 6.0, 12.0, 18.0, 24.0, 36.0, 40.0 })
        {
            const double g = std::pow (10.0, wantDb / 20.0);   // desired makeup
            const double M = D / (g * g);                       // wet energy that demands it
            const double expectDb = std::clamp ((double) gainDb ((float) g),
                                                (double) gainDb (kMinGain), (double) gainDb (kMaxGain));
            const bool clamped = (wantDb < -24.0 || wantDb > 36.0);
            test::approx (settleDb (D, M, 512, 48000.0, 6.0), expectDb, clamped ? 0.02 : 0.01,
                "did not converge to sqrt(D/M) for " + std::to_string (wantDb) + " dB");
        }
    }

    //================================================================ 2. per-sample dB/s HARD ceiling
    // The centrepiece. Across thousands of adversarial segments the per-sample applied-gain step in
    // dB must NEVER exceed the ceiling of the ACTIVE regime: fast (kSnapSlew) while snapGliding(),
    // else normal (kSlew). snapGliding() read BEFORE getNextGain() == the `fast` the step uses.
    // NO seed()/snapRatioTo() here (those legitimately teleport the applied gain) — every sample-to-
    // sample move must be a getNextGain glide, so any breach is a real ceiling violation.
    // Bound = step + max(5%·step, 1e-5 dB): the 5% covers the float rounding of the per-sample ratio
    // (~1..3% of the step at 192 kHz / large gains); it is still 88x below the fast/normal gap, so a
    // regime confusion or a genuinely too-fast glide is caught outright.
    test::group ("HARD ceiling: per-sample applied dB/s never exceeds the active regime (property)");
    {
        std::mt19937_64 rng (0xA1E7ULL);
        auto U = [&] (double lo, double hi) { return std::uniform_real_distribution<double> (lo, hi) (rng); };
        auto pick = [&] (int lo, int hi) { return std::uniform_int_distribution<int> (lo, hi) (rng); };

        double worstOverBy = -1.0;   // how far the worst sample exceeded its bound (dB); <0 => never
        bool   allFinite = true;
        long   totalSamples = 0;
        const int kSegments = 1500;

        for (int seg = 0; seg < kSegments; ++seg)
        {
            const double sr = rates[(size_t) pick (0, (int) rates.size() - 1)];
            AutoLeveler a; a.prepare (sr);
            bool enabled = true;
            const int blocks = pick (3, 18);
            for (int b = 0; b < blocks; ++b)
            {
                const int blk = pick (1, 2048);
                // hostile energies: log-uniform over ~ −80..+20 dBFS², sometimes sub-floor / huge
                const double dry = std::pow (10.0, U (-8.0, 1.0));
                const double mix = std::pow (10.0, U (-9.0, 2.0));
                if (U (0.0, 1.0) < 0.20) enabled = ! enabled;          // random enable/disable flips
                if (U (0.0, 1.0) < 0.15) a.openTransitionWindow();     // random transition windows
                a.processBlock (dry, mix, enabled, blk);
                for (int s = 0; s < blk; ++s)
                {
                    const bool   fast = a.snapGliding();               // == the regime getNextGain will use
                    const float  g0   = a.currentGain();
                    const float  g1   = a.getNextGain();
                    allFinite = allFinite && std::isfinite (g1) && g1 > 0.0f;
                    const double step  = (fast ? kSnapSlew : kSlew) / sr;
                    const double bound = step + std::max (0.05 * step, 1.0e-5);
                    worstOverBy = std::max (worstOverBy, std::fabs (dbDelta (g1, g0)) - bound);
                    ++totalSamples;
                }
            }
        }
        std::printf ("      (%ld samples; worst overshoot vs ceiling = %.3g dB)\n", totalSamples, worstOverBy);
        test::ok (worstOverBy <= 0.0, "a per-sample step exceeded the active dB/s ceiling by "
                                      + std::to_string (worstOverBy) + " dB");
        test::ok (allFinite, "applied gain went non-finite / non-positive under the hostile schedule");
    }

    //================================================================ 3. block-size invariance
    // The frozen-leveler bug this design replaced was block-size-DEPENDENT crawl (0.4 dB/s @ 64-sample
    // blocks vs the designed 9 dB/s). Pin it hard: drive a big +30 dB target from reset (rawDb ≈ +30
    // from block 1, both followers climb together) and measure the realized applied-gain slope over a
    // window strictly INSIDE the ramp (0.5 s..2.8 s, both < 30/9 = 3.33 s). It must read kSlew for
    // primes, size-1, and alternating 1/2048 alike; and the per-sample trajectories must coincide.
    test::group ("block-size invariance: realized slew rate INDEPENDENT of block size (primes, size-1, alt 1/2048)");
    {
        const double sr = 48000.0, D = 1.0, M = 1.0e-3;   // sqrt(D/M) = 31.62 => +30 dB
        const long n1 = std::lround (0.5 * sr), n2 = std::lround (2.8 * sr), nCap = n2 + 4096;
        const auto ref = perSampleDb (D, M, { 512 }, sr, nCap);
        const double refRate = (ref[(size_t) n2] - ref[(size_t) n1]) / ((double) (n2 - n1) / sr);
        // The reference chop is a correct float-quantised 9 dB/s ramp (the frozen-leveler bug crawled
        // at ~0.4 dB/s here). Pin it to the predicted float rate, not a nominal 9 it never equals.
        test::approx (refRate, predictedRate (sr), 0.02, "reference (blk=512) slew rate off the 9 dB/s ramp");

        const std::vector<std::vector<int>> patterns {
            { 1 }, { 2 }, { 3 }, { 7 }, { 64 }, { 1021 }, { 2048 }, { 1, 2048 }, { 1, 2, 3, 5, 7 } };
        for (const auto& pat : patterns)
        {
            const auto t = perSampleDb (D, M, pat, sr, nCap);
            const double rate = (t[(size_t) n2] - t[(size_t) n1]) / ((double) (n2 - n1) / sr);
            std::string tag = "blk={"; for (int b : pat) tag += std::to_string (b) + ","; tag += "}";
            // THE anti-freeze contract: every chop realises the SAME rate as blk=512 (invariance),
            // not a block-size-dependent crawl. 0.01 dB/s tol => any real chop dependence is caught.
            test::approx (rate, refRate, 0.01, tag + " slew rate depends on block size (frozen-leveler bug)");

            // trajectory coincidence through the ramp: the algorithm is block-size invariant, but the
            // per-block float summation of the target ramp (more block boundaries => more clamp
            // corrections) diverges the traces by O(1e-3) dB — well inside the RATE-invariance claim.
            float sup = 0.0f;
            for (long i = 0; i <= n2; i += 97) sup = std::max (sup, std::fabs (t[(size_t) i] - ref[(size_t) i]));
            test::ok (sup < 0.01f, tag + " trajectory diverges from blk=512 by "
                                   + std::to_string (sup) + " dB across the ramp");
        }
    }

    //================================================================ 4. sample-rate behaviour of the rate
    // The rate is defined in dB/SECOND, so it is *time*-invariant by construction — up to float. Two
    // independent quantisations govern the realised APPLIED ramp, and the slower one wins:
    //   • the per-sample GEOMETRIC ratio (a float member) → 9 dB/s ± ~0.6% (e.g. +0.51% @ 96k, −0.60% @ 192k);
    //   • the block-rate TARGET slew stepDb=(float)(9·blk/sr) → ≈ 9.000 dB/s (coarser step, tighter).
    // So realised = min(geometric, target). FINDING (reported): the float geometric ratio alone would
    // exceed the nominal 9 dB/s ceiling by ~0.5% at 88.2/96 kHz — but here the target-slew clamp holds
    // the applied ramp to 9.000, so the audible ramp stays on-ceiling; the raw ~0.5% overshoot only
    // surfaces per-sample when the target is far ahead (see the HARD-ceiling property's 5% float band).
    // Benign, sub-1%, deterministic; kept bit-identical to shipped OrbitCab (pinned, not "corrected").
    test::group ("sample-rate: realized rate == min(float geometric ramp, target slew) (44.1..192 kHz)");
    {
        const double D = 1.0, M = 1.0e-3;   // +30 dB target
        const int blk = 256;
        for (double sr : rates)
        {
            const long n1 = std::lround (0.5 * sr), n2 = std::lround (2.8 * sr), nCap = n2 + 4096;
            const auto t = perSampleDb (D, M, { blk }, sr, nCap);
            const double rate = (t[(size_t) n2] - t[(size_t) n1]) / ((double) (n2 - n1) / sr);
            const double targetRate = (double) (float) (kSlew * (double) blk / sr) * (sr / (double) blk);
            const double predicted  = std::min (predictedRate (sr), targetRate);
            test::approx (rate, predicted, 0.02,
                "sr=" + std::to_string ((int) sr) + " realized rate off min(geometric,target)");
            test::ok (std::fabs (rate - kSlew) < 0.01 * kSlew,   // stays inside the sub-1% float band of 9
                "sr=" + std::to_string ((int) sr) + " realized rate strayed > 1% from 9 dB/s ("
                + std::to_string (rate) + ")");
        }
    }

    //================================================================ 5. seed() semantics
    test::group ("seed(): instant applied-gain seed + clamp");
    {
        for (float g : { 0.001f, kMinGain, 0.5f, 1.0f, 2.0f, 10.0f, kMaxGain, 1000.0f })
        {
            AutoLeveler a; a.prepare (48000.0);
            a.seed (g);
            // "the initial applied gain matches the seed math before any measurement has settled"
            test::approx (a.currentGain(), clampGain (g), 1.0e-6,
                "seed did not place applied gain at clamp(g) instantly, g=" + std::to_string (g));
            test::approx (a.currentGainDb(), gainDb (clampGain (g)), 1.0e-4,
                "seed dB mismatch for g=" + std::to_string (g));
        }
    }

    test::group ("seed(): (−18 dBFS, ·/g²) is a follower fixed point — feeding it holds g");
    {
        // seed installs dryMeanSq = kSeedMeanSq, mixMeanSq = kSeedMeanSq/g². Feeding EXACTLY that pair
        // must hold the makeup dead still (rawTarget = sqrt(g²) = g) — pins the −18 dBFS seed dry energy
        // AND the mixMeanSq = kSeedMeanSq/g² formula.
        for (float g : { 0.25f, 1.0f, 4.0f })
        {
            AutoLeveler a; a.prepare (48000.0);
            a.seed (g);
            const double mixSeed = kSeedMeanSq / ((double) g * g);
            float lo = 1.0e9f, hi = -1.0e9f;
            for (int i = 0; i < 4000; ++i)   // ~5 s
            {
                a.processBlock (kSeedMeanSq, mixSeed, true, 64);
                for (int s = 0; s < 64; ++s) a.getNextGain();
                lo = std::min (lo, a.currentGainDb()); hi = std::max (hi, a.currentGainDb());
            }
            test::ok (hi - lo < 0.01f, "seed pair is not a fixed point for g=" + std::to_string (g)
                                       + " (drifted " + std::to_string (hi - lo) + " dB)");
            test::approx (a.currentGainDb(), gainDb (g), 0.02, "seed fixed point settled off g");
        }
    }

    test::group ("seed(): −18 dBFS gives PROMPT handover (a 0 dBFS seed would out-mass real signal)");
    {
        // seed(1.0) => dry=mix=kSeedMeanSq. Engage a +6 dB program at realistic −18 dBFS energy
        // (dry = kSeedMeanSq pinned, mix = kSeedMeanSq/4). With the −18 dBFS seed the dry follower is
        // already at program energy, so the makeup hands over at the follower τ + 9 dB/s slew:
        //   • by t=0.3 s it must already have climbed > +2 dB (slew-limited: 0.3·9 = 2.7 dB reachable);
        //   • by t=1.0 s (> 6/9 slew + ~5τ) it must have converged to +6 dB.
        // A 0 dBFS (energy 1.0) seed out-masses the −18 dBFS program ~63:1 → the makeup is pinned near
        // 0 dB for ~5τ and misses BOTH bounds. So these two checks encode the −18 dBFS choice itself.
        AutoLeveler a; a.prepare (48000.0);
        a.seed (1.0f);
        const double mix = kSeedMeanSq / 4.0;   // +6 dB program
        auto runTo = [&] (double sec) {
            const long blocks = std::lround (sec * 48000.0 / 64.0);
            for (long i = 0; i < blocks; ++i)
            { a.processBlock (kSeedMeanSq, mix, true, 64); for (int s = 0; s < 64; ++s) a.getNextGain(); }
        };
        runTo (0.3);
        test::ok (a.currentGainDb() > 2.0f, "handover stalled at t=0.3 s (currentGainDb="
                                            + std::to_string (a.currentGainDb()) + ") — seed out-massed the signal");
        runTo (1.0 - 0.3);
        test::approx (a.currentGainDb(), 6.0, 0.15, "did not converge to +6 dB by t=1.0 s");
    }

    //================================================================ 6. pump-proofing
    // A tremolo'd wet over a stationary dry must NOT make the makeup track the tremolo (= pumping, the
    // design's documented anti-goal). The 150 ms follower attenuates a wet-energy sinusoid at f by
    // G = 1/sqrt(1+(2πfτ)²); the makeup peak-to-peak is therefore bounded by 10·log10((1+βG)/(1−βG))
    // — and the 9 dB/s slew clamps it further. A leveler that TRACKED the tremolo would show the full
    // 10·log10((1+β)/(1−β)). Assert the applied ripple sits under the follower bound (with a small
    // slew/bias margin) and is a fraction of the tracking ripple.
    test::group ("pump-proof: tremolo'd wet -> gain ripple bounded by the 150 ms follower (no pump)");
    {
        const double sr = 48000.0, D = 0.25, M0 = 0.25;   // makeup ~0 dB
        const double f = 5.0, beta = 0.5, twoPi = 2.0 * std::acos (-1.0);
        const int blk = 64;
        const double G = 1.0 / std::sqrt (1.0 + std::pow (twoPi * f * kTau, 2.0));
        const double followerPP = 10.0 * std::log10 ((1.0 + beta * G) / (1.0 - beta * G));   // ~0.9 dB
        const double trackingPP = 10.0 * std::log10 ((1.0 + beta) / (1.0 - beta));           // ~4.77 dB

        AutoLeveler a; a.prepare (sr);
        const long blocks = std::lround (5.0 * sr / blk);
        float lo = 1.0e9f, hi = -1.0e9f;
        for (long i = 0; i < blocks; ++i)
        {
            const double t = (double) (i * blk) / sr;
            const double mix = M0 * (1.0 + beta * std::sin (twoPi * f * t));
            a.processBlock (D, mix, true, blk);
            for (int s = 0; s < blk; ++s) a.getNextGain();
            if (t > 2.5)   // measure after the followers reach steady oscillation
                { lo = std::min (lo, a.currentGainDb()); hi = std::max (hi, a.currentGainDb()); }
        }
        const double appliedPP = hi - lo;
        std::printf ("      (applied pp=%.3f dB ; follower bound=%.3f ; tracking=%.3f)\n",
                     appliedPP, followerPP, trackingPP);
        test::ok (appliedPP <= followerPP + 0.10, "gain ripple " + std::to_string (appliedPP)
                  + " dB exceeds the follower bound " + std::to_string (followerPP) + " dB (pumping)");
        test::ok (appliedPP < 0.5 * trackingPP, "gain tracked the tremolo (pumping)");
    }

    //================================================================ 7. silence handling
    // Below the documented floor the target HOLDS — it must never slam to the +36 dB clamp. Real
    // silence decays dry AND wet together, so the ratio is preserved right up to the gate, and once
    // dryMeanSq crosses kFloorMeanSq the makeup is frozen. Converge to +6 dB, then feed dry=mix=0.
    test::group ("silence: sub-floor energy holds the makeup (never slams to +36 dB)");
    {
        const double sr = 48000.0, D = 1.0, g = 2.0, M = D / (g * g);   // +6 dB
        AutoLeveler a; a.prepare (sr);
        for (int i = 0; i < 4000; ++i)   // converge to +6 dB
        { a.processBlock (D, M, true, 64); for (int s = 0; s < 64; ++s) a.getNextGain(); }
        test::approx (a.currentGainDb(), 6.0, 0.05, "did not converge to +6 dB before silence");

        float hiGain = -1.0e9f;
        for (int i = 0; i < 8000; ++i)   // ~10 s of true silence
        {
            a.processBlock (0.0, 0.0, true, 64);
            for (int s = 0; s < 64; ++s) a.getNextGain();
            hiGain = std::max (hiGain, a.currentGain());
        }
        test::approx (a.currentGainDb(), 6.0, 0.05, "makeup drifted during silence");
        test::ok (hiGain < (float) g * 1.02f, "silence slammed the makeup upward toward the clamp (peaked at "
                  + std::to_string (gainDb (hiGain)) + " dB)");
        test::ok (hiGain < 0.25f * kMaxGain, "makeup approached the +36 dB clamp under silence");
        // and the floor is on DRY energy: dry sub-floor holds even if wet collapses hard.
        for (int i = 0; i < 4000; ++i)
        { a.processBlock (1.0e-9, 1.0e-15, true, 64); for (int s = 0; s < 64; ++s) a.getNextGain(); }
        test::ok (std::isfinite (a.currentGain()) && a.currentGain() < 0.25f * kMaxGain,
                  "sub-floor dry + collapsed wet slammed the makeup to the clamp");
    }

    //================================================================ 8. determinism + reset()
    // Identical schedule from a fresh prepare, and from reset() on a used instance, must produce a
    // BIT-IDENTICAL applied-gain sequence.
    test::group ("determinism: identical schedule -> bit-identical gains (fresh + after reset)");
    {
        std::mt19937_64 rng (0xD37EULL);
        struct Step { double dry, mix; bool en; int blk; bool trans; };
        std::vector<Step> sched;
        for (int i = 0; i < 400; ++i)
        {
            std::uniform_real_distribution<double> U (0.0, 1.0);
            sched.push_back ({ std::pow (10.0, -6.0 + 7.0 * U (rng)),
                               std::pow (10.0, -7.0 + 8.0 * U (rng)),
                               U (rng) > 0.15,
                               std::uniform_int_distribution<int> (1, 1024) (rng),
                               U (rng) < 0.1 });
        }
        auto play = [&] (AutoLeveler& a) {
            std::vector<float> g;
            for (const auto& st : sched)
            {
                if (st.trans) a.openTransitionWindow();
                a.processBlock (st.dry, st.mix, st.en, st.blk);
                for (int s = 0; s < st.blk; ++s) g.push_back (a.getNextGain());
            }
            return g;
        };
        AutoLeveler a1; a1.prepare (96000.0);
        AutoLeveler a2; a2.prepare (96000.0);
        const auto g1 = play (a1);
        const auto g2 = play (a2);
        test::ok (g1 == g2, "two fresh instances diverged on an identical schedule");

        a1.reset();
        const auto g3 = play (a1);
        test::ok (g1 == g3, "reset() did not restore a bit-identical starting state");
    }

    return test::report();
}

// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// THEORY-FIRST FALSIFICATION suite for measurement::detectLeadingSilence. Distinct from
// IrPostLeadingSilenceTests.cpp (which pins hand-derived GOLDEN indices + a differential oracle):
// every expectation here is DERIVED from the documented contract BEFORE trusting the code, then the
// implementation is challenged to meet it. The fingerprint is OrbitCab's SHIPPED behavior — if a
// derivation and the code ever disagree, the code wins and the deviation gets pinned loudly, never
// papered over.
//
// Attack surfaces (all numbers derived, not tuned to pass):
//   1. Pre-roll exactness   — the trim lands EXACTLY (int)(0.0002·sr) before the first supra-threshold
//      sample, clamped at 0. The cast TRUNCATES (44100 → 8.82 → 8, not round()'s 9 — observable at the
//      gate boundary and pinned here); the source doc says "~0.2 ms", so truncation is in-contract.
//   2. Minimum-lead gate    — strict '>': a lead of exactly (int)(0.0005·sr) is rejected ("nothing
//      meaningful (>~0.5 ms) precedes"), one sample more is returned. Both sides, several rates.
//   3. Threshold boundary   — strict '>': a sample exactly equal to 0.001f·peak (the float product,
//      peaks 1.0 and 0.7 — the latter's product rounds) is silence; +1 ULP is signal. ± sign.
//   4. Scale invariance     — k·IR → the SAME index. EXACT for powers of two even on boundary samples
//      (pow2 scaling is exact and the rounding grid scales uniformly while everything stays normal,
//      so 0.001f·(k·peak) == k·(0.001f·peak) bit-for-bit — verified constants). For non-pow2 k it
//      holds away from the threshold boundary (rounding may flip an exact-boundary sample).
//   5. Channel fold         — the threshold is 0.001·GLOBAL peak (a quiet channel above its OWN 0.1 %
//      but below the global one must NOT fire); the onset is the MIN over channels of the first
//      crossing (the shrinking-bound scan == min-fold, so channel order is irrelevant).
//   6. Rate equivalence     — sr enters ONLY via the two truncated constants: 48000 and 48100 share
//      (9, 24) → identical results on any input; 47999 has gate 23 → a lead of 24 flips from rejected
//      to accepted across that boundary. Constants verified as exact IEEE-754 doubles.
//   7. Idempotence          — detect(trim(x, detect(x))) == 0 ALWAYS: a returned trim r = onset−preRoll
//      relocates the first crossing to exactly preRoll (samples before onset are sub-threshold on every
//      channel, and the peak survives, so the threshold is unchanged), giving lead' = 0 ≤ gate.
//   8. Contract fuzz        — a fully independent model (double-fold peak → float threshold, full-scan
//      min-fold onset, table constants) must predict the result on every random input; plus the derived
//      lemma that the peak sample itself always crosses (|peak| > 0.001f·peak whenever peak > 0), so an
//      onset always exists. Runs under ASan+UBSan in the FELITRONICS_ENABLE_SANITIZERS build.

#include <felitronics_test.h>
#include <felitronics/measurement/IrPost.h>

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstddef>
#include <random>
#include <string>
#include <utility>
#include <vector>

using namespace felitronics;
using felitronics::test::ok;
using felitronics::test::group;

using Buf = std::vector<std::vector<float>>;

namespace
{
    int detect (const Buf& buf, int numSamples, double sr)
    {
        std::vector<const float*> ptrs;
        ptrs.reserve (buf.size());
        for (const auto& ch : buf) ptrs.push_back (ch.data());
        return measurement::detectLeadingSilence (ptrs.data(), (int) buf.size(), numSamples, sr);
    }

    Buf deltaAt (int len, int idx, float amp = 1.0f)
    {
        Buf b (1, std::vector<float> ((std::size_t) len, 0.0f));
        b[0][(std::size_t) idx] = amp;
        return b;
    }

    // Derived constants: preRoll = trunc(0.0002·sr), gate = trunc(0.0005·sr). Every double product
    // below was verified exactly (IEEE-754, round-to-nearest — platform-stable):
    //   8000 → 1.6/4.0 · 22050 → 4.41/11.025 · 44100 → 8.82/22.05 · 48000 → 9.6/24 (exact) ·
    //   88200 → 17.64/44.1 · 96000 → 19.2/48 (exact) · 192000 → 38.4/96 (exact)
    struct RateConsts { double sr; int preRoll, gate; };
    constexpr RateConsts kRates[] = {
        {   8000.0,  1,  4 }, {  22050.0,  4, 11 }, {  44100.0,  8, 22 }, {  48000.0,  9, 24 },
        {  88200.0, 17, 44 }, {  96000.0, 19, 48 }, { 192000.0, 38, 96 },
    };
} // namespace

//==============================================================================
int main()
{
    // --- 1. pre-roll exactness: delta at 400 → 400 − preRoll, at every rate ---
    group ("pre-roll exactness (truncating cast, derived per rate)");
    for (const auto& rc : kRates)
    {
        const int got = detect (deltaAt (500, 400), 500, rc.sr);
        ok (got == 400 - rc.preRoll, "delta at 400 @ " + std::to_string ((int) rc.sr) + " Hz → 400-"
                                         + std::to_string (rc.preRoll) + " (got " + std::to_string (got) + ")");
    }
    {
        // Clamp at 0: onset earlier than the pre-roll cannot go negative (and is gated anyway).
        ok (detect (deltaAt (500, 4), 500, 48000.0) == 0, "onset before the pre-roll clamps at 0");
    }

    // --- 2. minimum-lead gate: strict '>' on both sides, at every rate ---
    group ("minimum-lead gate (strict, both sides of trunc(0.0005*sr))");
    for (const auto& rc : kRates)
    {
        const int reject = rc.preRoll + rc.gate;        // lead == gate → rejected
        ok (detect (deltaAt (600, reject), 600, rc.sr) == 0,
            "lead == gate rejected @ " + std::to_string ((int) rc.sr) + " Hz");
        ok (detect (deltaAt (600, reject + 1), 600, rc.sr) == rc.gate + 1,
            "lead == gate+1 returned @ " + std::to_string ((int) rc.sr) + " Hz");
    }

    // --- 3. threshold boundary: exactly 0.001f·peak is silence, +1 ULP is signal ---
    group ("threshold boundary (strict '>', +-1 ULP around 0.001f*peak)");
    for (float peak : { 1.0f, 0.7f })                   // 0.001f·0.7f rounds — probe the real product
    {
        const float thr = 0.001f * peak;                // the contract's threshold, evaluated in float
        auto base = deltaAt (2000, 1900, peak);

        auto atThr = base; atThr[0][100] = thr;
        ok (detect (atThr, 2000, 48000.0) == 1891,
            "sample == threshold ignored (peak " + std::to_string (peak) + ") → onset stays at 1900");

        auto above = base; above[0][100] = std::nextafter (thr, 1.0f);
        ok (detect (above, 2000, 48000.0) == 91,
            "one ULP above fires (peak " + std::to_string (peak) + ") → onset 100");

        auto neg = base; neg[0][100] = -std::nextafter (thr, 1.0f);
        ok (detect (neg, 2000, 48000.0) == 91, "negative one-ULP-above fires the same (|x|)");
    }

    // --- 4. scale invariance ---
    group ("scale invariance (pow2 exact incl. boundary; non-pow2 away from it)");
    {
        // Boundary-adversarial base: peak 0.7f at 1900, a sample EXACTLY at threshold at 100 (must
        // stay silent), the real onset at 300 (one ULP above). Expected 300-9 = 291 at 48 kHz.
        const float p = 0.7f, thr = 0.001f * p;
        auto base = deltaAt (2000, 1900, p);
        base[0][100] = thr;
        base[0][300] = std::nextafter (thr, 1.0f);
        ok (detect (base, 2000, 48000.0) == 291, "boundary-adversarial base case → 291");

        for (float k : { 0.5f, 2.0f, 1024.0f, 1.0f / 1024.0f, 6.103515625e-05f })  // powers of two
        {
            auto scaled = base;
            for (auto& v : scaled[0]) v *= k;
            ok (detect (scaled, 2000, 48000.0) == 291,
                "pow2 scale k=" + std::to_string (k) + " → identical index (boundary preserved)");
        }
        for (float k : { 3.0f, 0.37f })                 // non-pow2: margin-safe signals only
        {
            auto margin = deltaAt (2000, 1900, 1.0f);
            margin[0][300] = 0.5f;                      // 500x the threshold — no rounding hazard
            for (auto& v : margin[0]) v *= k;
            ok (detect (margin, 2000, 48000.0) == 291,
                "non-pow2 scale k=" + std::to_string (k) + " → identical index away from the boundary");
        }
    }

    // --- 5. channel fold semantics ---
    group ("channel fold (global peak, min over channels, order-free)");
    {
        Buf quiet (2, std::vector<float> (1200, 0.0f));
        quiet[0][1000] = 1.0f;
        quiet[1][100]  = 0.0005f;                       // 2500x its OWN channel's 0.1 %, half the global
        ok (detect (quiet, 1200, 48000.0) == 991, "quiet channel below the GLOBAL threshold cannot fire");

        Buf stag (3, std::vector<float> (1200, 0.0f));
        stag[0][800] = 1.0f; stag[1][300] = 0.01f; stag[2][500] = 0.4f;
        const int want = 291;                           // min(800, 300, 500) - 9
        ok (detect (stag, 1200, 48000.0) == want, "onset = min over channels of the first crossing");
        std::swap (stag[0], stag[1]);
        ok (detect (stag, 1200, 48000.0) == want, "channel permutation (swap 0,1) → same result");
        std::swap (stag[1], stag[2]);
        ok (detect (stag, 1200, 48000.0) == want, "channel permutation (rotate) → same result");
    }

    // --- 6. rate equivalence classes ---
    group ("rate equivalence (sr enters only via the two truncated constants)");
    {
        // 48000 and 48100 both derive (preRoll 9, gate 24) → identical everywhere.
        for (int idx : { 33, 34, 100, 9, 0 })
            ok (detect (deltaAt (400, idx), 400, 48000.0) == detect (deltaAt (400, idx), 400, 48100.0),
                "48000 vs 48100 identical (delta at " + std::to_string (idx) + ")");
        // 47999 derives gate 23 (0.0005*47999 = 23.9995) → a lead of exactly 24 flips to ACCEPTED.
        ok (detect (deltaAt (400, 33), 400, 48000.0) == 0,  "lead 24 rejected at 48000 (gate 24)");
        ok (detect (deltaAt (400, 33), 400, 47999.0) == 24, "lead 24 accepted at 47999 (gate 23)");
    }

    // --- 7 + 8. contract fuzz: full independent model + idempotence + never-eats-energy ---
    group ("contract fuzz — independent model, idempotence, never-eats-energy (2500 cases)");
    {
        std::mt19937 rng (0x1ead511eu);
        auto u01 = [&] { return std::uniform_real_distribution<float> (0.0f, 1.0f)(rng); };
        int modelMisses = 0, notIdempotent = 0, ateEnergy = 0, noOnsetLemma = 0;
        for (int it = 0; it < 2500; ++it)
        {
            const auto& rc  = kRates[rng() % (sizeof (kRates) / sizeof (kRates[0]))];
            const int nch   = 1 + (int) (rng() % 4);
            const int total = (int) (rng() % 1025);
            Buf b ((std::size_t) nch, std::vector<float> ((std::size_t) std::max (total, 1), 0.0f));
            const int scenario = (int) (rng() % 5);
            for (auto& chan : b)
            {
                if (scenario == 0) continue;                                        // silence
                if (scenario == 1 || scenario == 2)                                 // deltas (+ floor)
                {
                    if (scenario == 2)
                    {
                        const float floorAmp = std::pow (10.0f, -6.0f + 4.0f * u01());
                        for (auto& v : chan) v = floorAmp * (2.0f * u01() - 1.0f);
                    }
                    for (int k = 0, nk = 1 + (int) (rng() % 4); k < nk && total > 0; ++k)
                        chan[rng() % (std::size_t) total] =
                            (u01() < 0.5f ? -1.0f : 1.0f) * std::pow (10.0f, -6.0f + 6.3f * u01());
                }
                else if (scenario == 3)                                             // dense noise
                    for (auto& v : chan) v = 2.0f * u01() - 1.0f;
                else                                                                // denormal sprinkle
                    for (int k = 0; k < 16 && total > 0; ++k)
                        chan[rng() % (std::size_t) total] = 1.0e-41f * (float) (1 + rng() % 1000);
            }

            // Independent model. Peak via an exact double-fold (|float| and max are exact, and the
            // result is a widened member of the input set, so the cast back to float is lossless —
            // bit-equal to any in-float fold); threshold per the contract, 0.001f·peak in float;
            // onset via a full per-channel scan + min (== the shrinking-bound fold by derivation).
            double peakD = 0.0;
            for (const auto& chan : b)
                for (int i = 0; i < total; ++i)
                    peakD = std::max (peakD, (double) std::fabs (chan[(std::size_t) i]));
            const float peak = (float) peakD;
            int expected = 0;
            int firstCross = INT_MAX;
            if (total > 0 && peak > 0.0f)
            {
                const float thr = 0.001f * peak;
                for (const auto& chan : b)
                    for (int i = 0; i < total; ++i)
                        if (std::abs (chan[(std::size_t) i]) > thr) { firstCross = std::min (firstCross, i); break; }
                if (firstCross == INT_MAX) ++noOnsetLemma;      // derived: impossible while peak > 0
                else
                {
                    const int lead = std::max (0, firstCross - rc.preRoll);
                    expected = lead > rc.gate ? lead : 0;
                }
            }

            const int got = detect (b, total, rc.sr);
            if (got != expected) ++modelMisses;

            // Never-eats-energy: everything before the trim is <= 0.001f·peak on EVERY channel.
            if (got > 0)
            {
                const float thr = 0.001f * peak;
                for (const auto& chan : b)
                    for (int i = 0; i < got; ++i)
                        if (std::abs (chan[(std::size_t) i]) > thr) { ++ateEnergy; break; }
            }

            // Idempotence: apply the returned trim, re-detect → must be 0.
            {
                Buf t (b.size());
                for (std::size_t ch = 0; ch < b.size(); ++ch)
                    t[ch].assign (b[ch].begin() + got, b[ch].begin() + std::max (total, got));
                if (detect (t, total - got, rc.sr) != 0) ++notIdempotent;
            }
        }
        ok (modelMisses == 0,   "independent contract model predicts every result ("
                                    + std::to_string (modelMisses) + " misses)");
        ok (noOnsetLemma == 0,  "peak sample always crosses its own threshold (onset exists)");
        ok (ateEnergy == 0,     "never eats early energy: all trimmed samples <= 0.001*peak on every channel");
        ok (notIdempotent == 0, "idempotent: trim once, re-detect → 0 (" + std::to_string (notIdempotent) + ")");
    }

    return felitronics::test::report();
}

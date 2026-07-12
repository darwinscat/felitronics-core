// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// PARITY suite for measurement::detectLeadingSilence — the promoted OrbitCab HEAD-trim detector
// (orbitcab src/core/IRMath.h, cab::detectLeadingSilence). The fingerprint is OrbitCab's SHIPPED
// behavior: trim indices must stay BIT-IDENTICAL so the consumer's class-1 shim can't move a single
// shipped HEAD-trim point. Two independent oracles, both must agree:
//   1. Hand-derived golden indices — every expectation below is derived by executing the documented
//      algorithm on paper (constants verified as exact IEEE-754 doubles, platform-stable):
//        sr      preRoll=(int)(0.0002·sr)   gate=(int)(0.0005·sr)   (returned lead must be > gate)
//        44100   8.8200000000000003 → 8     22.050000000000001 → 22
//        48000   9.5999999999999996 → 9     24 (exact)         → 24
//        96000   19.199999999999999 → 19    48 (exact)         → 48
//   2. A transcription of the OrbitCab source (below), faithful for FINITE inputs — differential
//      fuzz (finite-only by construction) proves the promoted function equals it bit-for-bit across
//      random channel counts / lengths / rates / content, including threshold-adversarial samples
//      exactly at 0.001f·peak and denormals. On NaN the transcription mirrors JUCE's SCALAR fold
//      (seed-first) while shipped OrbitCab may take a SIMD path with different NaN semantics — so
//      the oracle is never consulted outside the finite domain.
//   3. Outside the golden domain core DEFINES the behavior (see the header's DOMAIN note): the NaN
//      fold, the Inf-threshold degeneration, and the input guards are pinned directly against
//      hand-derived values, without the transcription oracle.

#include <felitronics_test.h>
#include <felitronics/measurement/IrPost.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <random>
#include <string>
#include <utility>
#include <vector>

using namespace felitronics;
using felitronics::test::ok;
using felitronics::test::group;

using Buf = std::vector<std::vector<float>>;

//==============================================================================
// Transcription of OrbitCab's cab::detectLeadingSilence (src/core/IRMath.h), faithful for FINITE
// inputs — the only domain where the tests consult it — with juce::AudioBuffer<float> replaced by
// raw channel vectors. JUCE equivalences, exact on that domain:
//   * juce::jmax (a, b) == (a < b ? b : a)  — transcribed literally.
//   * buf.getMagnitude (ch, 0, n) == findMinAndMax over the range, then the max of {min, -min,
//     max, -max}. Compare/negate only — every intermediate is a member of the input set or its
//     negation, no rounding anywhere — so for finite inputs it equals the scalar max-|x| fold
//     bit-for-bit (on NaN it seeds from the first element like JUCE's scalar path, while shipped
//     JUCE may take a SIMD path with different NaN semantics — hence the finite-domain scope).
namespace orbitcab
{
    inline float jmax (float a, float b) { return a < b ? b : a; }
    inline int   jmax (int a, int b)     { return a < b ? b : a; }

    inline float getMagnitude (const float* d, int n)
    {
        if (n <= 0) return 0.0f;                                     // JUCE: empty range → Range() → 0
        float mn = d[0], mx = d[0];                                  // FloatVectorOperations::findMinAndMax
        for (int i = 1; i < n; ++i) { mn = d[i] < mn ? d[i] : mn; mx = d[i] > mx ? d[i] : mx; }
        return jmax (jmax (mn, -mn), jmax (mx, -mx));                // max of {start, -start, end, -end}
    }

    inline int detectLeadingSilence (const Buf& buf, int numSamples, double sampleRate)
    {
        const int total = numSamples;
        const int nch   = (int) buf.size();
        if (total <= 0 || nch <= 0)
            return 0;

        float peak = 0.0f;
        for (int ch = 0; ch < nch; ++ch)
            peak = jmax (peak, getMagnitude (buf[(std::size_t) ch].data(), total));
        if (peak <= 0.0f)
            return 0;

        const float thresh = 0.001f * peak;
        int onset = total;
        for (int ch = 0; ch < nch && onset > 0; ++ch)
        {
            const float* d = buf[(std::size_t) ch].data();
            for (int i = 0; i < onset; ++i)
                if (std::abs (d[i]) > thresh) { onset = i; break; }
        }

        const int preRoll = (int) (0.0002 * sampleRate);    // ~0.2 ms
        const int lead    = jmax (0, onset - preRoll);
        return lead > (int) (0.0005 * sampleRate) ? lead : 0;
    }
} // namespace orbitcab

//==============================================================================
namespace
{
    int promoted (const Buf& buf, int numSamples, double sr)
    {
        std::vector<const float*> ptrs;
        ptrs.reserve (buf.size());
        for (const auto& ch : buf) ptrs.push_back (ch.data());
        return measurement::detectLeadingSilence (ptrs.data(), (int) buf.size(), numSamples, sr);
    }

    // Both oracles at once (finite inputs only): promoted == OrbitCab transcription == hand-derived.
    void golden (const Buf& buf, int n, double sr, int want, const std::string& what)
    {
        const int got = promoted (buf, n, sr);
        ok (got == want, what + " — hand-derived index (got " + std::to_string (got)
                              + ", want " + std::to_string (want) + ")");
        ok (got == orbitcab::detectLeadingSilence (buf, n, sr), what + " — matches the OrbitCab source");
    }

    Buf mono (int len) { return Buf (1, std::vector<float> ((std::size_t) len, 0.0f)); }

    Buf deltaAt (int len, int idx, float amp = 1.0f)
    {
        auto b = mono (len);
        b[0][(std::size_t) idx] = amp;
        return b;
    }
} // namespace

//==============================================================================
int main()
{
    // --- deltas across rates: pre-roll subtraction + minimum-lead gate, hand-derived ---
    group ("48 kHz deltas (preRoll 9, gate 24)");
    {
        golden (deltaAt (400, 0),   400, 48000.0, 0,   "delta at 0 → nothing precedes");
        golden (deltaAt (400, 100), 400, 48000.0, 91,  "delta at 100 → 100-9");
        golden (deltaAt (400, 33),  400, 48000.0, 0,   "delta at 33 → lead 24 == gate, strict '>' rejects");
        golden (deltaAt (400, 34),  400, 48000.0, 25,  "delta at 34 → lead 25, first accepted");
        golden (deltaAt (400, 9),   400, 48000.0, 0,   "delta at 9 → lead 0");
        golden (deltaAt (400, 5),   400, 48000.0, 0,   "delta at 5 → onset-preRoll clamps at 0");
        golden (deltaAt (400, 100, -0.8f), 400, 48000.0, 91, "negative delta → |x| detected the same");
    }

    group ("44.1 kHz deltas (preRoll 8 — the cast truncates 8.82, gate 22)");
    {
        golden (deltaAt (400, 30),  400, 44100.0, 0,   "delta at 30 → lead 22 == gate → 0");
        golden (deltaAt (400, 31),  400, 44100.0, 23,  "delta at 31 → lead 23 (round-to-9 would give 22 → 0)");
        golden (deltaAt (400, 400 - 1), 400, 44100.0, 391, "delta at 399 → 399-8");
    }

    group ("96 kHz deltas (preRoll 19, gate 48)");
    {
        golden (deltaAt (400, 67),  400, 96000.0, 0,   "delta at 67 → lead 48 == gate → 0");
        golden (deltaAt (400, 68),  400, 96000.0, 49,  "delta at 68 → lead 49");
    }

    // --- ramped onset: quadratic ramp x[500+j] = j²·1e-6f, j = 0..1000 (peak ≈ 1 at 1500) ---
    // thresh ≈ 1e-3·peak; first j with j²·1e-6 > 1e-3 ⇒ j² > 1000 ⇒ j = 32 (31² = 961 is 3.9 % below,
    // 32² = 1024 is 2.4 % above — no ULP hazard). Onset 532 → lead 532-9 = 523.
    group ("ramped onset (48 kHz)");
    {
        auto b = mono (2000);
        for (int j = 0; j <= 1000; ++j)
            b[0][(std::size_t) (500 + j)] = (float) (j * j) * 1.0e-6f;
        golden (b, 2000, 48000.0, 523, "quadratic ramp crosses 0.001·peak at j=32");
    }

    // --- pre-onset noise floors relative to a unit delta at 1000 (48 kHz; thresh = 0.001) ---
    // -50 dB = 3.16e-3 > thresh → the forward scan fires at sample 0 → under-trims to 0 (the safe
    // direction). -70 dB = 3.16e-4 and -90 dB = 3.16e-5 < thresh → floor ignored → 1000-9 = 991.
    group ("pre-onset noise floors");
    {
        for (double dB : { -50.0, -70.0, -90.0 })
        {
            auto b = deltaAt (1200, 1000);
            const float amp = (float) std::pow (10.0, dB / 20.0);
            for (int i = 0; i < 1000; ++i)
                b[0][(std::size_t) i] = (i & 1) ? -amp : amp;
            golden (b, 1200, 48000.0, dB == -50.0 ? 0 : 991,
                    "floor at " + std::to_string ((int) dB) + " dB");
        }
    }

    // --- multi-channel: global peak, min-over-channels onset ---
    group ("multi-channel fold");
    {
        Buf b (2, std::vector<float> (1200, 0.0f));
        b[0][800] = 1.0f;                       // the peak channel
        b[1][300] = 0.5f;                       // earlier, quieter onset — still > 0.001·globalPeak
        golden (b, 1200, 48000.0, 291, "staggered onsets → earliest channel wins (300-9)");
        std::swap (b[0], b[1]);
        golden (b, 1200, 48000.0, 291, "channel order swapped → same trim");

        Buf q (2, std::vector<float> (1200, 0.0f));
        q[0][1000] = 1.0f;
        q[1][100]  = 0.0005f;                   // above its own channel's 0.001·chPeak, BELOW the global one
        golden (q, 1200, 48000.0, 991, "threshold is 0.001·GLOBAL peak — quiet channel can't fire");
    }

    // --- degenerate inputs ---
    group ("degenerate");
    {
        golden (mono (0), 0, 48000.0, 0, "empty buffer");
        golden (mono (100), 100, 48000.0, 0, "all-zero");
        Buf negZero (1, std::vector<float> (100, -0.0f));
        golden (negZero, 100, 48000.0, 0, "all minus-zero → peak <= 0 bail");
        golden (deltaAt (1, 0), 1, 48000.0, 0, "single nonzero sample");
        ok (promoted (Buf {}, 100, 48000.0) == 0, "zero channels → 0");
        golden (deltaAt (400, 100), -3, 48000.0, 0, "negative numSamples → 0");
    }

    // --- denormal-range signal: thresh itself lands in the denormals, scan still exact ---
    group ("denormal range");
    {
        auto b = deltaAt (600, 500, 1.0e-37f);              // normal peak; thresh 1e-40 is denormal
        for (int i = 0; i < 500; ++i)
            b[0][(std::size_t) i] = 1.4e-45f;               // min denormal — below thresh
        golden (b, 600, 48000.0, 491, "denormal floor below a denormal threshold");
    }

    // --- defined domain: non-finite samples + input guards (core-pinned, no transcription oracle —
    //     see the header's DOMAIN note; the shipped JUCE NaN fold is platform-dependent SIMD) ---
    group ("defined domain: NaN/Inf pins + input guards");
    {
        const float nan = std::numeric_limits<float>::quiet_NaN();
        const float inf = std::numeric_limits<float>::infinity();

        // NaN is always the SECOND argument of the fold's std::max → deterministically ignored,
        // and an IEEE compare with NaN in the scan is false → never fires.
        auto nanLead = deltaAt (400, 100);
        nanLead[0][0] = nan;
        ok (promoted (nanLead, 400, 48000.0) == 91, "NaN-leading channel: NaN ignored → 100-9 = 91");

        Buf allNaN (1, std::vector<float> (100, nan));
        ok (promoted (allNaN, 100, 48000.0) == 0, "all-NaN → peak stays 0 → don't trim");

        Buf nanCh (2, std::vector<float> (400, 0.0f));
        for (auto& v : nanCh[0]) v = nan;
        nanCh[1][100] = 1.0f;
        ok (promoted (nanCh, 400, 48000.0) == 91, "all-NaN channel beside a clean one → clean decides");

        // Inf is ordered → propagates into the peak exactly as in OrbitCab: thresh = Inf, nothing
        // crosses, onset degenerates to numSamples → the gated tail (numSamples - preRoll).
        ok (promoted (deltaAt (400, 200, inf), 400, 48000.0) == 391,
            "Inf peak → thresh Inf, no crossing → gated tail 400-9 = 391 (matches OrbitCab: Inf ordered)");
        auto infPlus = deltaAt (1200, 800);
        infPlus[0][500] = inf;
        ok (promoted (infPlus, 1200, 48000.0) == 1191,
            "Inf hijacks the threshold over finite signal → gated tail 1200-9 = 1191");

        // Input guards: every invalid input → 0, noexcept-safe, no UB (the sampleRate guard also
        // kills the (int)(0.0002*NaN/Inf) casts; the clamp kills the absurd-finite-rate overflow).
        const auto good = deltaAt (400, 100);
        ok (promoted (good, 400, std::numeric_limits<double>::quiet_NaN()) == 0, "NaN sampleRate → 0");
        ok (promoted (good, 400, std::numeric_limits<double>::infinity()) == 0,  "Inf sampleRate → 0");
        ok (promoted (good, 400, 0.0) == 0,      "zero sampleRate → 0");
        ok (promoted (good, 400, -48000.0) == 0, "negative sampleRate → 0");
        ok (promoted (good, 400, 1.0e18) == 0,   "absurd finite rate → clamped (int) casts, gated to 0");

        const float* oneNull[2] = { good[0].data(), nullptr };
        ok (measurement::detectLeadingSilence (oneNull, 2, 400, 48000.0) == 0, "one null channel → 0");
        ok (measurement::detectLeadingSilence (nullptr, 2, 400, 48000.0) == 0, "null channel array → 0");
    }

    // --- differential fuzz: promoted vs the OrbitCab transcription — finite inputs only ---
    group ("differential fuzz vs the OrbitCab source (4000 finite cases)");
    {
        std::mt19937 rng (0x0c4b1eadu);
        auto u01 = [&] { return std::uniform_real_distribution<float> (0.0f, 1.0f)(rng); };
        const double rates[] = { 8000.0, 11025.0, 22050.0, 44100.0, 48000.0, 88200.0,
                                 96000.0, 176400.0, 192000.0, 44100.5, 12345.6 };
        int mismatches = 0, invariantBreaks = 0;
        for (int it = 0; it < 4000; ++it)
        {
            const int nch   = (int) (rng() % 5);            // 0..4 — includes the no-channel edge
            const int total = (int) (rng() % 1025);         // 0..1024
            const double sr = (rng() % 4 == 0) ? 1000.0 + 399000.0 * (double) u01()
                                               : rates[rng() % (sizeof (rates) / sizeof (rates[0]))];
            Buf b ((std::size_t) nch, std::vector<float> ((std::size_t) std::max (total, 1), 0.0f));
            const int scenario = (int) (rng() % 6);
            for (auto& chan : b)
            {
                if (scenario == 0) continue;                                    // all zero
                if (scenario == 1 || scenario == 2)                             // sparse deltas (+ floor)
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
                else if (scenario == 3)                                         // dense noise
                    for (auto& v : chan) v = 2.0f * u01() - 1.0f;
                else if (scenario == 4 && total > 0)                            // threshold-adversarial:
                {                                                               // exact 0.001f·peak samples
                    chan[(std::size_t) (total - 1)] = 1.0f;
                    const float thr = 0.001f * 1.0f;
                    for (int k = 0; k < 8 && total > 1; ++k)
                    {
                        const std::size_t i = rng() % (std::size_t) (total - 1);
                        const int pick = (int) (rng() % 3);
                        chan[i] = pick == 0 ? thr
                                : pick == 1 ? std::nextafter (thr, 1.0f)
                                            : -std::nextafter (thr, 0.0f);
                    }
                }
                else if (scenario == 5)                                         // denormal sprinkle
                    for (int k = 0; k < 16 && total > 0; ++k)
                        chan[rng() % (std::size_t) total] = 1.0e-41f * (float) (1 + rng() % 1000);
            }
            const int got  = promoted (b, total, sr);
            const int want = orbitcab::detectLeadingSilence (b, total, sr);
            if (got != want) ++mismatches;
            const int gate = (int) (0.0005 * sr);
            if (! (got == 0 || (got > gate && got < total))) ++invariantBreaks;
        }
        ok (mismatches == 0, "bit-identical to the OrbitCab source on every fuzz case ("
                                 + std::to_string (mismatches) + " mismatches)");
        ok (invariantBreaks == 0, "result is 0 or in (gate, numSamples) on every fuzz case");
    }

    return felitronics::test::report();
}

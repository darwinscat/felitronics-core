// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// Theory-first falsification suite for saturation::Saturator. Every expectation was derived from the
// DOCUMENTED contract (Saturator.h + WaveShaper.h header docs) plus FP/DSP theory BEFORE running, and
// none may be weakened to make a run pass:
//
//  1) ODD SYMMETRY. Tanh/Atan are odd curves; Cubic is an odd polynomial; Asym with bias=0 degenerates
//     to the odd tanh (biasTanh=tanh(0)=0, |r(1)|=|r(-1)|). The documented topology for symmetric
//     shapes has NO DC blocker ("DC-block only the asymmetric curve"), and every remaining stage —
//     gate clamp, FIR up/down-sampling (sums of c*x), constant gains (drive, autoComp, trim), the
//     linear dry/wet, the dry delay — is EXACTLY negation-equivariant under IEEE-754 round-to-nearest
//     (rounding is sign-symmetric: fl(-a op -b) == -fl(a op b); libm tanh/atan compute on |x| and
//     transfer the sign). Asym's DC blocker is LTI from zero state, so it is negation-equivariant too.
//     => two independent, identically-prepared instances fed x and -x must produce EXACT negations
//     (FP ==, which is bit-identity except the +0/-0 pair). Any asymmetry beyond that is a finding.
//
//  2) MONOTONICITY. The Tanh static curve y=tanh(kx)/tanh(k) is strictly monotone; multiplication by
//     a positive constant is monotone under RN rounding. At os=1 there is no FIR and no DC blocker
//     (symmetric shape) => a monotone ramp must map to an EXACTLY non-decreasing output, zero
//     tolerance. At os=4 the polyphase FIRs are LTI with tiny passband ripple; a full-scale ramp of
//     16384 samples is spectrally near-DC, so ripple can bend adjacent samples by at most
//     ~ripple x per-sample-step (~1e-4 x 1e-3) plus float ulp noise => bound any dip by 1e-6.
//
//  3) POISON CONTAINMENT. The gate sanitizes IN PLACE BEFORE any state ("each input sample is
//     sanitized at the gate (NaN/Inf -> 0, huge finite -> +/-1e6 clamp) so one bad sample can't lodge
//     in the oversampler/DC/dry-delay state") => a poisoned stream must be BIT-IDENTICAL to a clean
//     twin fed the documented substitution, everywhere, for every shape. Against the TRUE clean twin,
//     outputs may differ only causally after the hit and must reconverge:
//       - shapes without the DC blocker have FINITE memory: up FIR (tapsPerPhase base samples) +
//         down FIR (~taps) + dry delay (latency) => bit-exact reconvergence within a horizon
//         H = 4*latency + 128 (generous cover of that sum);
//       - Asym adds the one-pole DC blocker: pole a = exp(-2*pi*10Hz/(4*48kHz)) per oversampled
//         sample => over a 2000-base-sample window the residual shrinks by a^8000 ~ 0.073; assert a
//         conservative x0.5 per-window geometric decay (floored at 2e-7 float-ulp noise) and a final
//         residual <= 1e-6. Sustained Inf bursts recover identically (gate feeds zeros during burst).
//
//  4) CHUNK ASSOCIATIVITY. "n may exceed maxBlock - chunked internally, state carries across chunks
//     -> bit-identical to one big call" => ANY call-size decomposition of the same stream is
//     bit-identical: thousands of random chunkings x random finite params, incl. zero-length calls
//     (a no-op by theory: no samples => no state advance).
//
//  5) +/-1e6 GATE BOUNDARY. Clamp is inclusive: exactly +/-1e6 passes AS IS, 1 ulp above maps to
//     exactly the bound (bit-identical to a twin fed the bound), 1 ulp below passes untouched;
//     outputs stay finite ("the clamp also catches an extreme finite input that would overflow
//     downstream").
//
//  6) PARAM SANITY. "Non-finite params fall back to the struct defaults" — per field, other fields
//     unaffected => poisoning one field must be bit-identical to setting that field to its default.

#include <felitronics_test.h>
#include <felitronics/saturation/WaveShaper.h>
#include <felitronics/saturation/Saturator.h>

#include <bit>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <random>
#include <string>
#include <vector>

using namespace felitronics;
using Shape  = saturation::WaveShaper::Shape;
using Sat    = saturation::Saturator;
using Params = saturation::Saturator::Params;

static const char* shapeName (Shape s)
{
    switch (s) { case Shape::Tanh: return "Tanh"; case Shape::Atan: return "Atan";
                 case Shape::Cubic: return "Cubic"; default: return "Asym"; }
}

//==============================================================================
// Planar test buffers.
struct Planar
{
    Planar (int channels, int n) : ch ((size_t) channels, std::vector<float> ((size_t) n)), ptrs ((size_t) channels)
    { for (size_t i = 0; i < ch.size(); ++i) ptrs[i] = ch[i].data(); }

    float* const* data() { return ptrs.data(); }
    std::vector<std::vector<float>> ch;
    std::vector<float*> ptrs;
};

static bool bitEq (float a, float b) { return std::bit_cast<std::uint32_t> (a) == std::bit_cast<std::uint32_t> (b); }

static void fillSignal (Planar& b, unsigned seed, float amp = 0.95f)
{
    std::mt19937 rng (seed);
    std::uniform_real_distribution<float> uni (-0.35f, 0.35f);
    for (size_t c = 0; c < b.ch.size(); ++c)
        for (size_t n = 0; n < b.ch[c].size(); ++n)
        {
            const double ph = 0.0371 + (double) c * 0.313;
            float v = amp * (0.6f * (float) std::sin (2.0 * 3.14159265358979 * (0.007 + 0.002 * (double) c) * (double) n + ph)) + uni (rng);
            if (v == 0.0f) v = 1e-3f;   // keep the odd-symmetry comparison free of +0/-0 bookkeeping
            b.ch[c][n] = v;
        }
}

static bool prep (Sat& s, const Params& p, double fs = 48000.0, int maxBlock = 512, int chans = 1, int os = 4)
{
    const bool okPrep = s.prepare (fs, maxBlock, chans, os);
    if (okPrep) s.setParams (p);
    return okPrep;
}

// Documented gate mapping (oracle written from the doc, not the body): NaN/Inf -> 0, clamp to +/-1e6 inclusive.
static float gateOracle (float x)
{
    if (! std::isfinite (x)) return 0.0f;
    if (x >  1.0e6f) return  1.0e6f;
    if (x < -1.0e6f) return -1.0e6f;
    return x;
}

static const float kQNaN = std::bit_cast<float> (0x7FC00000u);
static const float kSNaN = std::bit_cast<float> (0x7FA00000u);
static const float kInf  = std::numeric_limits<float>::infinity();

//==============================================================================
int main()
{
    std::printf ("SaturatorTheoryTests (theory-first falsification: contract derived before running)\n");

    //==========================================================================
    // 1) Odd symmetry — steady-state paired instances: process(x) and process(-x) must be exact negations.
    test::group ("odd symmetry: f(-x) == -f(x) exactly (paired instances; symmetric shapes bypass the DC blocker)");
    for (Shape shape : { Shape::Tanh, Shape::Atan, Shape::Cubic, Shape::Asym })   // Asym with bias=0 is odd by its formula
        for (float driveDb : { 0.0f, 3.0f, 36.0f })
            for (int os : { 1, 4 })
            {
                Params p; p.shape = shape; p.driveDb = driveDb; p.bias = 0.0f;
                p.mix = 0.7f; p.outputDb = 1.5f; p.autoComp = 0.5f;

                const int N = 4096;
                Planar a (1, N), b (1, N);
                fillSignal (a, 42u + (unsigned) os);
                for (int n = 0; n < N; ++n) b.ch[0][(size_t) n] = -a.ch[0][(size_t) n];

                Sat sa, sb;
                test::ok (prep (sa, p, 48000.0, 512, 1, os) && prep (sb, p, 48000.0, 512, 1, os), "prepare");
                sa.process (a.data(), 1, N);
                sb.process (b.data(), 1, N);

                int bad = 0;
                for (int n = 0; n < N; ++n)
                    if (! (b.ch[0][(size_t) n] == -a.ch[0][(size_t) n])) ++bad;   // FP ==: bit-identity except +/-0
                test::ok (bad == 0, std::string ("odd symmetry ") + shapeName (shape)
                                    + " drive=" + std::to_string (driveDb) + " os=" + std::to_string (os)
                                    + " (mismatches: " + std::to_string (bad) + "/" + std::to_string (N) + ")");
            }

    //==========================================================================
    // 2) Monotonicity — Tanh ramp.
    test::group ("monotonicity: Tanh ramp non-decreasing (os=1 exact; os=4 dip <= 1e-6 from FIR passband ripple)");
    for (float driveDb : { 0.0f, 36.0f })
        for (int os : { 1, 4 })
        {
            Params p; p.shape = Shape::Tanh; p.driveDb = driveDb; p.mix = 1.0f; p.outputDb = 0.0f; p.autoComp = 0.5f;

            const int N = 16384;
            Planar r (1, N);
            for (int n = 0; n < N; ++n) r.ch[0][(size_t) n] = -1.0f + 2.0f * (float) n / (float) (N - 1);

            Sat s;
            test::ok (prep (s, p, 48000.0, 512, 1, os), "prepare");
            const int skip = s.latencySamples() + 256;   // settle: FIR group delay + full kernel span
            s.process (r.data(), 1, N);

            float maxDip = 0.0f;
            for (int n = skip + 1; n < N; ++n)
            {
                const float dip = r.ch[0][(size_t) n - 1] - r.ch[0][(size_t) n];
                if (dip > maxDip) maxDip = dip;
            }
            const float tol = (os == 1) ? 0.0f : 1e-6f;   // os=1: memoryless monotone chain, zero tolerance
            test::ok (maxDip <= tol, "Tanh ramp monotone drive=" + std::to_string (driveDb) + " os=" + std::to_string (os)
                                     + " (max dip " + std::to_string (maxDip) + ", tol " + std::to_string (tol) + ")");
        }

    //==========================================================================
    // 3a) Gate substitution equivalence: poisoned stream == clean twin fed the documented substitution, bitwise.
    test::group ("poison gate: poisoned stream bit-identical to the doc-mapped substitution twin (all shapes)");
    for (Shape shape : { Shape::Tanh, Shape::Atan, Shape::Cubic, Shape::Asym })
    {
        Params p; p.shape = shape; p.driveDb = 6.0f; p.bias = 0.25f; p.mix = 0.7f;

        const int N = 2048, CH = 2;
        Planar poisoned (CH, N);
        fillSignal (poisoned, 777u);

        std::mt19937 rng (901u);
        const float poisonVals[] = { kQNaN, kSNaN, kInf, -kInf, 3.7e10f, -2.9e38f, -7.7e33f, 8.8e6f };
        for (int c = 0; c < CH; ++c)
        {
            for (int k = 0; k < 40; ++k)
                poisoned.ch[(size_t) c][rng() % (unsigned) N] = poisonVals[k % 8];
            for (int n = 700;  n < 764;  ++n) poisoned.ch[(size_t) c][(size_t) n] = (c == 0 ? kInf : -kInf);   // Inf burst
            for (int n = 1200; n < 1264; ++n) poisoned.ch[(size_t) c][(size_t) n] = kQNaN;                     // NaN burst
        }

        Planar twin (CH, N);
        for (int c = 0; c < CH; ++c)
            for (int n = 0; n < N; ++n)
                twin.ch[(size_t) c][(size_t) n] = gateOracle (poisoned.ch[(size_t) c][(size_t) n]);

        Sat sp, st;
        test::ok (prep (sp, p, 48000.0, 512, CH) && prep (st, p, 48000.0, 512, CH), "prepare");
        sp.process (poisoned.data(), CH, N);
        st.process (twin.data(), CH, N);

        int bad = 0, nonFinite = 0;
        for (int c = 0; c < CH; ++c)
            for (int n = 0; n < N; ++n)
            {
                if (! bitEq (poisoned.ch[(size_t) c][(size_t) n], twin.ch[(size_t) c][(size_t) n])) ++bad;
                if (! std::isfinite (poisoned.ch[(size_t) c][(size_t) n])) ++nonFinite;
            }
        test::ok (bad == 0,       std::string ("gate substitution bit-identical, ") + shapeName (shape) + " (" + std::to_string (bad) + " diffs)");
        test::ok (nonFinite == 0, std::string ("output all-finite under poison, ")  + shapeName (shape));
    }

    //==========================================================================
    // 3b) Containment horizon vs the TRUE clean twin.
    test::group ("poison horizon: single NaN — finite-memory shapes reconverge bit-exactly within H = 4*latency + 128");
    for (Shape shape : { Shape::Tanh, Shape::Atan, Shape::Cubic })
    {
        Params p; p.shape = shape; p.driveDb = 6.0f; p.mix = 0.7f;

        const int N = 4096, j = 1500;
        Planar clean (1, N), pois (1, N);
        fillSignal (clean, 4242u);
        for (int n = 0; n < N; ++n) pois.ch[0][(size_t) n] = clean.ch[0][(size_t) n];
        pois.ch[0][j] = kQNaN;

        Sat sc, sp;
        test::ok (prep (sc, p) && prep (sp, p), "prepare");
        const int H = 4 * sc.latencySamples() + 128;   // >= up-FIR + down-FIR + dry-delay memory, generous
        sc.process (clean.data(), 1, N);
        sp.process (pois.data(), 1, N);

        int preDiff = 0, tailDiff = 0, nonFinite = 0;
        for (int n = 0; n < N; ++n)
        {
            if (! std::isfinite (pois.ch[0][(size_t) n])) ++nonFinite;
            if (n < j     && ! bitEq (pois.ch[0][(size_t) n], clean.ch[0][(size_t) n])) ++preDiff;    // causality
            if (n >= j + H && ! bitEq (pois.ch[0][(size_t) n], clean.ch[0][(size_t) n])) ++tailDiff;  // finite memory
        }
        test::ok (nonFinite == 0, std::string ("all-finite, ") + shapeName (shape));
        test::ok (preDiff == 0,   std::string ("causal: prefix before the hit bit-identical, ") + shapeName (shape));
        test::ok (tailDiff == 0,  std::string ("bit-exact reconvergence after H, ") + shapeName (shape)
                                  + " (H=" + std::to_string (H) + ", diffs " + std::to_string (tailDiff) + ")");
    }

    test::group ("poison horizon: Asym one-pole DC tail — geometric decay (<=0.5x per 2000-sample window; theory 0.073x) and final residual <= 1e-6");
    for (int burstLen : { 1, 256 })   // single NaN hit, then a sustained Inf burst
    {
        Params p; p.shape = Shape::Asym; p.driveDb = 6.0f; p.bias = 0.25f; p.mix = 0.7f;

        const int j = 1500, tail = 26000, N = j + burstLen + tail;
        Planar clean (1, N), pois (1, N);
        fillSignal (clean, 515u);
        for (int n = 0; n < N; ++n) pois.ch[0][(size_t) n] = clean.ch[0][(size_t) n];
        for (int n = j; n < j + burstLen; ++n) pois.ch[0][(size_t) n] = (burstLen == 1 ? kQNaN : kInf);

        Sat sc, sp;
        test::ok (prep (sc, p) && prep (sp, p), "prepare");
        sc.process (clean.data(), 1, N);
        sp.process (pois.data(), 1, N);

        int nonFinite = 0;
        for (int n = 0; n < N; ++n) if (! std::isfinite (pois.ch[0][(size_t) n])) ++nonFinite;
        test::ok (nonFinite == 0, "Asym all-finite under poison, burst=" + std::to_string (burstLen));

        const int W = 2000, kWins = tail / W;   // 13 windows from the end of the burst
        std::vector<float> wmax ((size_t) kWins, 0.0f);
        for (int k = 0; k < kWins; ++k)
            for (int n = j + burstLen + k * W; n < j + burstLen + (k + 1) * W; ++n)
            {
                const float d = std::fabs (pois.ch[0][(size_t) n] - clean.ch[0][(size_t) n]);
                if (d > wmax[(size_t) k]) wmax[(size_t) k] = d;
            }

        bool decays = true;
        for (int k = 1; k < kWins; ++k)
            if (wmax[(size_t) k] > std::max (0.5f * wmax[(size_t) k - 1], 2e-7f)) decays = false;   // 2e-7 = float-ulp noise floor
        test::ok (decays, "Asym windowed residual decays geometrically, burst=" + std::to_string (burstLen));
        test::ok (wmax[(size_t) kWins - 1] <= 1e-6f, "Asym final residual <= 1e-6, burst=" + std::to_string (burstLen)
                                                     + " (got " + std::to_string (wmax[(size_t) kWins - 1]) + ")");
    }

    //==========================================================================
    // 4) Chunk associativity — thousands of random chunkings, bit-identical to the one-call reference.
    test::group ("chunk associativity: 3000 random chunkings x random finite params — bit-identical");
    {
        std::mt19937 rng (0xACDCu);
        std::uniform_real_distribution<float> uSig (-1.2f, 1.2f);
        std::uniform_real_distribution<float> uDrive (0.0f, 36.0f), uBias (-0.5f, 0.5f), uMix (0.0f, 1.0f),
                                              uOut (-6.0f, 6.0f), uComp (0.0f, 1.0f), uDc (5.0f, 40.0f);
        const int maxBlock = 192;
        int badTrials = 0, prepFails = 0;

        for (int trial = 0; trial < 3000; ++trial)
        {
            Params p;
            p.shape    = static_cast<Shape> (trial % 4);
            p.driveDb  = uDrive (rng); p.bias = uBias (rng); p.mix = uMix (rng);
            p.outputDb = uOut (rng);   p.autoComp = uComp (rng); p.dcBlockHz = uDc (rng);
            const int os = (trial % 2 == 0) ? 4 : 1;
            const int ch = (trial % 3 == 0) ? 2 : 1;
            const int N  = 200 + (int) (rng() % 600u);

            Planar ref (ch, N), split (ch, N);
            for (int c = 0; c < ch; ++c)
                for (int n = 0; n < N; ++n)
                    split.ch[(size_t) c][(size_t) n] = ref.ch[(size_t) c][(size_t) n] = uSig (rng);

            Sat sr, ss;
            if (! prep (sr, p, 48000.0, maxBlock, ch, os) || ! prep (ss, p, 48000.0, maxBlock, ch, os)) { ++prepFails; continue; }

            sr.process (ref.data(), ch, N);   // one call; N > maxBlock exercises the internal chunker

            int pos = 0;
            float* sub[core::kMaxChannels] {};
            while (pos < N)
            {
                if (rng() % 10u == 0u)                       // zero-length call: a no-op by theory
                {
                    for (int c = 0; c < ch; ++c) sub[c] = split.ch[(size_t) c].data() + pos;
                    ss.process (sub, ch, 0);
                }
                int m = 1 + (int) (rng() % (unsigned) (2 * maxBlock + 41));
                if (m > N - pos) m = N - pos;
                for (int c = 0; c < ch; ++c) sub[c] = split.ch[(size_t) c].data() + pos;
                ss.process (sub, ch, m);
                pos += m;
            }

            bool same = true;
            for (int c = 0; c < ch && same; ++c)
                for (int n = 0; n < N; ++n)
                    if (! bitEq (ref.ch[(size_t) c][(size_t) n], split.ch[(size_t) c][(size_t) n])) { same = false; break; }
            if (! same) ++badTrials;
        }
        test::ok (prepFails == 0, "all trial prepares succeeded (" + std::to_string (prepFails) + " failures)");
        test::ok (badTrials == 0, "3000 random chunkings bit-identical to one-call reference ("
                                  + std::to_string (badTrials) + " diverging trials)");
    }

    //==========================================================================
    // 5) +/-1e6 gate boundary — inclusive clamp, exact at the bound, 1 ulp either side.
    test::group ("gate boundary: +/-1e6 inclusive — at-bound passes, +1 ulp clamps to the bound bitwise, -1 ulp untouched");
    for (Shape shape : { Shape::Tanh, Shape::Cubic })
    {
        Params p; p.shape = shape; p.driveDb = (shape == Shape::Cubic ? 36.0f : 3.0f); p.mix = 0.7f;

        const float bound = 1.0e6f;
        const float specials[] = {
            bound,  std::nextafterf (bound,  kInf), std::nextafterf (bound,  0.0f),
            -bound, std::nextafterf (-bound, -kInf), std::nextafterf (-bound, 0.0f),
            std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(), 8.5e7f, -3.3e6f,
        };

        const int N = 256;
        Planar in (1, N), twin (1, N);
        fillSignal (in, 33u, 0.4f);
        for (int n = 0, k = 0; n < N; n += 17, ++k) in.ch[0][(size_t) n] = specials[(size_t) (k % 10)];
        for (int n = 0; n < N; ++n) twin.ch[0][(size_t) n] = gateOracle (in.ch[0][(size_t) n]);

        Sat si, st;
        test::ok (prep (si, p) && prep (st, p), "prepare");
        si.process (in.data(), 1, N);
        st.process (twin.data(), 1, N);

        int bad = 0, nonFinite = 0;
        for (int n = 0; n < N; ++n)
        {
            if (! bitEq (in.ch[0][(size_t) n], twin.ch[0][(size_t) n])) ++bad;
            if (! std::isfinite (in.ch[0][(size_t) n])) ++nonFinite;
        }
        test::ok (bad == 0,       std::string ("boundary mapping bit-identical to doc oracle, ") + shapeName (shape));
        test::ok (nonFinite == 0, std::string ("finite output at extreme finite inputs, ")       + shapeName (shape));
    }

    //==========================================================================
    // 6) Non-finite params fall back to that field's default; other fields keep their set values.
    test::group ("param sanity: each non-finite field falls back to its struct default (bit-identical to the explicit twin)");
    {
        const Params defaults;                    // struct defaults per the header
        Params base; base.shape = Shape::Asym; base.driveDb = 9.0f; base.bias = 0.3f; base.mix = 0.6f;
        base.outputDb = -3.0f; base.autoComp = 0.8f; base.dcBlockHz = 25.0f;

        struct FieldRef { const char* name; float Params::* member; };
        const FieldRef fields[] = { { "driveDb", &Params::driveDb }, { "bias", &Params::bias },
                                    { "mix", &Params::mix },         { "outputDb", &Params::outputDb },
                                    { "autoComp", &Params::autoComp }, { "dcBlockHz", &Params::dcBlockHz } };

        for (const auto& f : fields)
            for (float poison : { kQNaN, kInf, -kInf })
            {
                Params pp = base; pp.*(f.member) = poison;
                Params pt = base; pt.*(f.member) = defaults.*(f.member);

                const int N = 1024;
                Planar a (1, N), b (1, N);
                fillSignal (a, 606u);
                for (int n = 0; n < N; ++n) b.ch[0][(size_t) n] = a.ch[0][(size_t) n];

                Sat sa, sb;
                test::ok (prep (sa, pp) && prep (sb, pt), "prepare");
                sa.process (a.data(), 1, N);
                sb.process (b.data(), 1, N);

                int bad = 0;
                for (int n = 0; n < N; ++n)
                    if (! bitEq (a.ch[0][(size_t) n], b.ch[0][(size_t) n])) ++bad;
                test::ok (bad == 0, std::string ("non-finite ") + f.name + " == default-field twin ("
                                    + std::to_string (bad) + " diffs)");
            }
    }

    //==========================================================================
    // Totality storm: heavy random poison mixture, every shape — output always finite.
    test::group ("totality: 30% poison storm (NaN/sNaN/Inf/huge/random-bits) — output always finite, every shape");
    for (Shape shape : { Shape::Tanh, Shape::Atan, Shape::Cubic, Shape::Asym })
    {
        Params p; p.shape = shape; p.driveDb = 12.0f; p.bias = 0.4f; p.mix = 0.5f;

        const int N = 8192, CH = 2;
        Planar buf (CH, N);
        std::mt19937 rng (1313u + (unsigned) shape);
        std::uniform_real_distribution<float> uni (-2.0f, 2.0f);
        for (int c = 0; c < CH; ++c)
            for (int n = 0; n < N; ++n)
            {
                if (rng() % 10u < 3u)
                {
                    switch (rng() % 6u)
                    {
                        case 0:  buf.ch[(size_t) c][(size_t) n] = kQNaN; break;
                        case 1:  buf.ch[(size_t) c][(size_t) n] = kSNaN; break;
                        case 2:  buf.ch[(size_t) c][(size_t) n] = kInf;  break;
                        case 3:  buf.ch[(size_t) c][(size_t) n] = -kInf; break;
                        case 4:  buf.ch[(size_t) c][(size_t) n] = std::numeric_limits<float>::max(); break;
                        default: buf.ch[(size_t) c][(size_t) n] = std::bit_cast<float> ((std::uint32_t) rng()); break;
                    }
                }
                else
                    buf.ch[(size_t) c][(size_t) n] = uni (rng);
            }

        Sat s;
        test::ok (prep (s, p, 48000.0, 512, CH), "prepare");
        s.process (buf.data(), CH, N);

        int nonFinite = 0;
        for (int c = 0; c < CH; ++c)
            for (int n = 0; n < N; ++n)
                if (! std::isfinite (buf.ch[(size_t) c][(size_t) n])) ++nonFinite;
        test::ok (nonFinite == 0, std::string ("poison storm all-finite, ") + shapeName (shape)
                                  + " (" + std::to_string (nonFinite) + " bad)");
    }

    return test::report();
}

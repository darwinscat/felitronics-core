// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// felitronics::mastering — the checks that must NOT trust the chain's own arithmetic.
//
// WHY THIS IS A SEPARATE BINARY. The obvious way to test the latency contract is "bypass everything,
// null the output against the input delayed by latencySamples()". That test is unsound on its own: the
// declared latency and the delay the bypass aligners hold come from the SAME number, so a wrong number
// satisfies both and the null passes. The falsifier is exact — `Compressor::prepare(.., maxLookaheadMs
// = 50)` caps its lookahead at 2400 samples, so a chain asking for 60 ms and computing lround(2880)
// itself would hold latency AND aligner at 2880 and pass, while the active chain sat 480 samples out.
//
// So everything here re-derives the answer from outside the chain:
//   * the delay is FOUND by searching for the shift that nulls, not read from the chain;
//   * the active chain's group delay is measured from the phase of a DFT bin, which knows nothing about
//     any of it (and is measured at a frequency whose period exceeds twice the delay, so the answer is
//     not ambiguous modulo a period — the trap an impulse-peak or a 1 kHz phase reading falls into);
//   * the composition is nulled against the same stages driven by hand, so a wrong ORDER, a missing
//     stage or a mis-forwarded parameter fails even though the chain is perfectly self-consistent;
//   * the passband is measured and pinned, because the stages' own suite stops at 0.357*fs and did not
//     notice that its header's droop figures were a factor of two out.

#include <felitronics/mastering/OfflineRenderer.h>
#include <felitronics_test.h>

#include <algorithm>
#include <bit>
#include <cmath>
#include <complex>
#include <cstdint>
#include <string>
#include <vector>

using namespace felitronics;
using felitronics::test::ok;
using felitronics::test::approx;
using felitronics::test::group;

namespace
{
constexpr double kPi = 3.14159265358979323846;
using Buf = std::vector<std::vector<float>>;

std::uint32_t bits (float f) noexcept { return std::bit_cast<std::uint32_t> (f); }

std::vector<float*> planes (Buf& b, int off = 0)
{
    std::vector<float*> p;
    p.reserve (b.size());
    for (auto& v : b) p.push_back (v.data() + off);
    return p;
}

Buf silence (int nch, int n) { return Buf ((std::size_t) nch, std::vector<float> ((std::size_t) n, 0.0f)); }

Buf tone (int nch, int n, double hz, float amp, double fs = 48000.0)
{
    Buf x = silence (nch, n);
    for (int i = 0; i < n; ++i)
        for (int c = 0; c < nch; ++c)
            x[(std::size_t) c][(std::size_t) i] = amp * (float) std::sin (2.0 * kPi * hz * i / fs);
    return x;
}

// One DFT bin, exactly on a bin centre (the caller picks n so that hz*n/fs is an integer).
std::complex<double> binAt (const std::vector<float>& v, int from, int n, double hz, double fs)
{
    std::complex<double> acc { 0.0, 0.0 };
    const double w = 2.0 * kPi * hz / fs;
    for (int i = 0; i < n; ++i)
    {
        const double t = (double) (from + i);
        acc += (double) v[(std::size_t) (from + i)] * std::complex<double> (std::cos (w * t), -std::sin (w * t));
    }
    return acc;
}

mastering::MasteringChainConfig cfgAll (int tpp = 64)
{
    mastering::MasteringChainConfig c;
    c.internalBlock = 128;
    c.eq = c.compressor = c.limiter = c.dither = true;
    c.monoBass = false;                                    // stereo-only; kept out so mono cases work too
    c.clipper  = true;
    c.tapsPerPhase = tpp;
    c.compressorLookaheadMs = 1.0;
    c.limiterLookaheadMs    = 1.0;
    return c;
}

// A chain that is fully ACTIVE — every stage is called — but linear at the frequencies used below, so
// its group delay is the honest algorithmic one and not something a gain reduction bent.
mastering::MasteringChainParams transparentButActive()
{
    mastering::MasteringChainParams p;
    p.compressor.ratio       = 1.0;                        // called, delays, applies exactly unity
    p.compressor.makeupDb    = 0.0;
    p.clipper.driveDb        = 0.0f;                       // called, oversampled round trip, ~linear
    p.clipper.mix            = 1.0f;
    p.clipper.autoComp       = 0.0f;
    p.limiter.ceilingDbTp    = 40.0;                       // called, round trip, nothing limited
    p.dither.bits            = 32;                         // called, early-returns (float export)
    return p;
}
} // namespace

//==============================================================================
// ORACLE 1 — find the delay by SEARCHING for it, then compare with what the chain claims.
static void testLatencyFoundNotAsked()
{
    group ("latency — found by search, not read from the chain");

    const int nch = 2, n = 30000;
    for (int tpp : { 32, 64 })
        for (double compLookMs : { 0.0, 1.0, 5.0 })
        {
            auto cfg = cfgAll (tpp);
            cfg.compressorLookaheadMs = compLookMs;

            mastering::MasteringChain chain;
            if (! chain.prepare (48000.0, nch, cfg)) { ok (false, "prepare"); continue; }

            mastering::MasteringChainParams p;             // everything present, everything bypassed
            p.bypassEq = p.bypassMonoBass = p.bypassCompressor = true;
            p.bypassClipper = p.bypassLimiter = p.bypassDither = true;
            chain.setParams (p);

            // A signal with no periodicity, so exactly one shift can null it.
            Buf x = silence (nch, n);
            std::uint32_t s = 2463534242u;
            for (int i = 0; i < n; ++i)
                for (int c = 0; c < nch; ++c)
                {
                    s ^= s << 13; s ^= s >> 17; s ^= s << 5;
                    x[(std::size_t) c][(std::size_t) i] = (float) ((double) s / 4294967296.0 - 0.5);
                }

            Buf y = x;
            { auto pl = planes (y); chain.process (pl.data(), nch, n); }

            int found = -1, matches = 0;
            for (int d = 0; d <= 4096; ++d)
            {
                bool all = true;
                for (int c = 0; c < nch && all; ++c)
                    for (int i = d; i < n && all; ++i)
                        if (bits (y[(std::size_t) c][(std::size_t) i]) != bits (x[(std::size_t) c][(std::size_t) (i - d)])) all = false;
                if (all) { ++matches; if (found < 0) found = d; }
            }
            const std::string tag = "tpp=" + std::to_string (tpp) + " compLookahead=" + std::to_string ((int) compLookMs) + "ms";
            ok (matches == 1, "exactly one shift nulls the bypassed chain — " + tag);
            ok (found == chain.latencySamples(),
                "the FOUND delay (" + std::to_string (found) + ") is the declared one ("
                + std::to_string (chain.latencySamples()) + ") — " + tag);
        }
}

//==============================================================================
// ORACLE 2 — the ACTIVE chain's group delay, measured from a DFT bin's phase.
static void testActiveGroupDelayByPhase()
{
    group ("latency — the ACTIVE chain, measured from DFT phase");

    const int nch = 1, n = 96000;
    const auto cfg = cfgAll (64);
    mastering::MasteringChain chain;
    ok (chain.prepare (48000.0, nch, cfg), "prepare");
    chain.setParams (transparentButActive());
    const int D = chain.latencySamples();

    // 25 Hz: its period is 1920 samples, comfortably more than twice any delay this chain can have, so
    // the phase determines the delay outright rather than modulo a period. Every FIR in the path is
    // flat and linear-phase down here, and the compressor's ring is an integer delay.
    const double f = 25.0, fs = 48000.0;
    Buf x = tone (nch, n, f, 0.25f);
    Buf y = x;
    { auto pl = planes (y); chain.process (pl.data(), nch, n); }

    const int from = 24000, len = 48000;                   // an integer number of periods, clear of both edges
    const auto bx = binAt (x[0], from, len, f, fs);
    const auto by = binAt (y[0], from, len, f, fs);
    double dphase = std::arg (by) - std::arg (bx);
    while (dphase >  kPi) dphase -= 2.0 * kPi;
    while (dphase < -kPi) dphase += 2.0 * kPi;
    const double measured = -dphase / (2.0 * kPi * f / fs);

    approx (measured, (double) D, 0.05,
            "the active chain's group delay equals its declared latency (" + std::to_string (D) + " samples)");

    // MUTATION WITNESS: this is the check that survives when latencySamples() and the aligners are
    // wrong TOGETHER. A one-sample error here is 0.0033 rad at 25 Hz and lands well outside the
    // tolerance above.
}

//==============================================================================
// ORACLE 3 — the composition, nulled against the same stages driven by hand.
static void testCompositionAgainstAHandBuiltChain()
{
    group ("composition — nulled against the stages driven by hand");

    const int nch = 2, n = 24000, K = 128;
    auto cfg = cfgAll (64);
    cfg.internalBlock = K;
    cfg.monoBass = true;
    cfg.sidechainHpfHz = 90.0;

    mastering::MasteringChainParams p;
    p.inputGainDb = 2.0;
    p.preLimiterGainDb = 4.0;
    p.eqBands[0].on = true;
    p.eqBands[0].type = eq::FilterType::Bell;
    p.eqBands[0].lane (eq::Lane::Stereo).on = true;
    p.eqBands[0].lane (eq::Lane::Stereo).freq = 2500.0;
    p.eqBands[0].lane (eq::Lane::Stereo).Q = 1.1;
    p.eqBands[0].lane (eq::Lane::Stereo).gainDb = -3.0;
    p.monoBass = { true, 130.0f, 0.0f };
    p.compressor.thresholdDb = -20.0;
    p.compressor.ratio = 3.0;
    p.compressor.attackMs = 8.0;
    p.compressor.releaseMs = 140.0;
    p.compressor.makeupDb = 2.0;
    p.clipper.driveDb = 5.0f;
    p.limiter.ceilingDbTp = -1.0;
    p.limiter.releaseMs = 50.0;
    p.dither.bits = 24;

    Buf x = silence (nch, n);
    for (int i = 0; i < n; ++i)
    {
        const double t = (double) i / 48000.0;
        const float e = (float) (0.5 * std::sin (2.0 * kPi * 70.0 * t) + 0.35 * std::sin (2.0 * kPi * 1800.0 * t));
        x[0][(std::size_t) i] = e;
        x[1][(std::size_t) i] = -0.7f * e;
    }

    // --- the chain ---
    mastering::MasteringChain chain;
    ok (chain.prepare (48000.0, nch, cfg), "prepare the chain");
    chain.setParams (p);
    Buf viaChain = x;
    { auto pl = planes (viaChain); chain.process (pl.data(), nch, n); }

    // --- the same stages, by hand, at the same quantum ---
    eq::EqEngine          eqE;
    stereo::MonoBass      mb;
    dynamics::Compressor  comp;
    saturation::Saturator sat;
    limiter::TruePeakLimiter lim;
    dither::Dither        dit;
    eq::Biquad            hpf[2];

    eqE.prepare (48000.0, K, nch);
    for (int i = 0; i < eq::EqEngine::kMaxBands; ++i) eqE.setBand (i, p.eqBands[i]);
    mb.prepare (48000.0, K, nch);
    mb.setParams (p.monoBass);
    ok (comp.prepare (48000.0, K, nch, std::max (cfg.compressorLookaheadMs, 1.0)), "hand: compressor prepare");
    { auto cp = p.compressor; cp.lookaheadMs = cfg.compressorLookaheadMs; comp.setParams (cp); }
    ok (sat.prepare (48000.0, K, nch, cfg.oversampleFactor, cfg.tapsPerPhase), "hand: saturator prepare");
    sat.setParams (p.clipper);
    limiter::TruePeakLimiterConfig lc;
    lc.lookaheadMs = cfg.limiterLookaheadMs; lc.oversampleFactor = cfg.oversampleFactor; lc.tapsPerPhase = cfg.tapsPerPhase;
    ok (lim.prepare (48000.0, K, nch, lc), "hand: limiter prepare");
    lim.setParams (p.limiter);
    dit.prepare (48000.0, K, nch);
    dit.setParams (p.dither);
    const eq::BiquadCoeffs hc = eq::matched::highpass (cfg.sidechainHpfHz, 48000.0, 0.70710678118654752);
    for (int c = 0; c < nch; ++c) hpf[c].setCoeffs (hc);

    const float gIn  = (float) core::dbToGain (p.inputGainDb);
    const float gPre = (float) core::dbToGain (p.preLimiterGainDb);

    Buf byHand = x;
    std::vector<float> key ((std::size_t) K * (std::size_t) nch, 0.0f);
    for (int off = 0; off + K <= n; off += K)
    {
        float* ch[2] { byHand[0].data() + off, byHand[1].data() + off };
        for (int c = 0; c < nch; ++c)
            for (int i = 0; i < K; ++i)
            {
                const float v = ch[c][i];
                ch[c][i] = std::clamp (std::isfinite (v) ? v : 0.0f, -1.0e6f, 1.0e6f);
            }
        for (int c = 0; c < nch; ++c) for (int i = 0; i < K; ++i) ch[c][i] *= gIn;
        eqE.process (ch, nch, K);
        mb.process (ch, nch, K);
        const float* kp[2] {};
        for (int c = 0; c < nch; ++c)
        {
            float* k = key.data() + (std::size_t) c * (std::size_t) K;
            for (int i = 0; i < K; ++i) k[i] = hpf[c].processSample (ch[c][i]);
            hpf[c].flushDenormals();
            kp[c] = k;
        }
        comp.process (ch, nch, K, kp, nch);
        sat.process (ch, nch, K);
        for (int c = 0; c < nch; ++c) for (int i = 0; i < K; ++i) ch[c][i] *= gPre;
        lim.process (ch, nch, K);
        dit.process (ch, nch, K);
    }

    // The chain's output trails the hand-built one by exactly the internal quantum.
    long long bad = 0;
    int firstBad = -1;
    for (int c = 0; c < nch; ++c)
        for (int i = 0; i + K < n; ++i)
            if (bits (viaChain[(std::size_t) c][(std::size_t) (i + K)]) != bits (byHand[(std::size_t) c][(std::size_t) i]))
            { ++bad; if (firstBad < 0) firstBad = i; }
    ok (bad == 0, "the chain is bit-identical to the same stages driven by hand ("
                  + std::to_string (bad) + " differ, first at " + std::to_string (firstBad) + ")");

    // MUTATION WITNESSES this one catches and the self-consistent tests do not: swapping two stages,
    // dropping the gate, applying inputGain after the EQ, forwarding the wrong params object, taking
    // the key from after the compressor instead of before it, or dropping the sidechain filter's flush.
}

//==============================================================================
// ORACLE 4 — the passband, measured and PINNED.
static void testPassbandPinned()
{
    group ("passband — the oversampled stages' cost, measured and pinned");

    // The stages' own suite stops at 0.357*fs, which is why nobody noticed that the limiter header's
    // "-0.40 dB at 0.40 fs, -4.0 at 0.44" are ONE FILTER PASS while the path is a round trip. Clipper
    // and limiter in series double it again. These numbers are measured on this tree; a change to the
    // prototype has to edit them on purpose.
    struct Row { double r; double tpp32; double tpp64; };
    const Row rows[] = {
        { 0.36,  -0.001,  0.000 },
        { 0.40,  -1.549,  0.000 },
        { 0.42,  -6.033, -0.610 },
        { 0.44, -16.131, -10.182 },
    };

    for (const Row& row : rows)
        for (int which = 0; which < 2; ++which)
        {
            const int tpp = which == 0 ? 32 : 64;
            auto cfg = cfgAll (tpp);
            cfg.internalBlock = 256;
            cfg.eq = false;
            cfg.dither = false;
            cfg.compressor = false;
            cfg.clipper = true;
            cfg.limiter = true;

            const double fs = 44100.0;                     // the CD target the chain is voiced for
            const int n = 44100;
            mastering::MasteringChain chain;
            if (! chain.prepare (fs, 1, cfg)) { ok (false, "prepare"); continue; }
            auto p = transparentButActive();
            chain.setParams (p);

            Buf x = tone (1, n, row.r * fs, 0.25f, fs);
            Buf y = x;
            { auto pl = planes (y); chain.process (pl.data(), 1, n); }

            double pin = 0.0, pout = 0.0;
            for (int i = n / 2; i < n; ++i)
            {
                pin  += (double) x[0][(std::size_t) i] * x[0][(std::size_t) i];
                pout += (double) y[0][(std::size_t) i] * y[0][(std::size_t) i];
            }
            const double db = 10.0 * std::log10 (pout / pin);
            const double want = which == 0 ? row.tpp32 : row.tpp64;
            approx (db, want, 0.06, "clipper+limiter at " + std::to_string (row.r) + "*fs, tpp="
                                    + std::to_string (tpp) + " (" + std::to_string ((int) (row.r * fs)) + " Hz)");
        }

    // And the reason the chain defaults to 64 rather than the stages' own 32, stated as a check.
    ok (mastering::MasteringChainConfig {}.tapsPerPhase == 64,
        "the chain's default tapsPerPhase is 64 — 32 costs -1.55 dB at 17.6 kHz with the clipper in");
}

//==============================================================================
// ORACLE 5 — the two gain nodes do different jobs, and it is measurable which.
static void testTwoGainNodes()
{
    group ("the two gain nodes — one drives the compressor, the other does not");

    const int nch = 1, n = 24000;
    auto cfg = cfgAll (64);
    cfg.clipper = false;
    cfg.dither  = false;
    cfg.eq      = false;

    mastering::MasteringChainParams base;
    base.compressor.thresholdDb = -24.0;
    base.compressor.ratio       = 6.0;
    base.compressor.attackMs    = 5.0;
    base.compressor.releaseMs   = 100.0;
    base.bypassLimiter          = true;                    // keep the limiter out of the arithmetic
    base.limiter.ceilingDbTp    = 0.0;

    Buf x = tone (nch, n, 220.0, 0.30f);

    auto run = [&] (double inDb, double preDb) {
        mastering::MasteringChain c;
        ok (c.prepare (48000.0, nch, cfg), "gain-node probe: prepare");
        auto p = base;
        p.inputGainDb = inDb;
        p.preLimiterGainDb = preDb;
        c.setParams (p);
        Buf y = x;
        auto pl = planes (y);
        c.process (pl.data(), nch, n);
        double acc = 0.0;
        for (int i = n / 2; i < n; ++i) acc += (double) y[0][(std::size_t) i] * y[0][(std::size_t) i];
        return std::sqrt (acc / (double) (n - n / 2));
    };

    const double a = run (0.0, 0.0);
    const double b = run (0.0, 6.0);
    const double c = run (6.0, 0.0);

    ok (a > 1.0e-4, "the reference render is not silence");

    // The two nodes are the same 6 dB in two places, and the compressor is what tells them apart.
    // Downstream of every dynamic stage the node is exactly a scale: 6 dB in, 6 dB out. Upstream of a
    // 6:1 compressor whose threshold the tone is already over, 6 dB in comes back as 6/ratio = 1 dB.
    // That ratio IS the reason P7 gets its own node: it needs a control variable that moves loudness
    // without re-running the compression it already solved for.
    const double dbPre = 20.0 * std::log10 (b / a);
    const double dbIn  = 20.0 * std::log10 (c / a);
    approx (dbPre, 6.0, 0.05, "preLimiterGainDb is a pure scale: +6 dB in gives +6 dB out");
    approx (dbIn,  1.0, 0.30, "inputGainDb DRIVES the compressor: +6 dB in gives back 6/ratio = 1 dB");
    ok (dbPre - dbIn > 4.0, "the same 6 dB at the two nodes differ by "
                            + std::to_string (dbPre - dbIn) + " dB — they are not interchangeable");
    // MUTATION WITNESS: swap the two nodes and both approx() checks take each other's value.
}

//==============================================================================
// ORACLE 6 — the sidechain high-pass is real, and 0 Hz is exactly the self-keyed path.
static void testSidechainKey()
{
    group ("sidechain — the internal key filter");

    const int nch = 1, n = 24000;
    auto cfg = cfgAll (64);
    cfg.clipper = false;
    cfg.dither  = false;
    cfg.eq      = false;

    mastering::MasteringChainParams p;
    p.compressor.thresholdDb = -26.0;
    p.compressor.ratio       = 8.0;
    p.compressor.attackMs    = 5.0;
    p.compressor.releaseMs   = 120.0;
    p.bypassLimiter          = true;

    auto rms = [&] (double hpfHz, double toneHz) {
        auto c2 = cfg;
        c2.sidechainHpfHz = hpfHz;
        mastering::MasteringChain c;
        ok (c.prepare (48000.0, nch, c2), "sidechain probe: prepare");
        c.setParams (p);
        Buf y = tone (nch, n, toneHz, 0.35f);
        auto pl = planes (y);
        c.process (pl.data(), nch, n);
        double acc = 0.0;
        for (int i = n / 2; i < n; ++i) acc += (double) y[0][(std::size_t) i] * y[0][(std::size_t) i];
        return std::sqrt (acc / (double) (n - n / 2));
    };

    // A loud 50 Hz tone: with a 500 Hz key filter the detector barely sees it, so it is barely
    // compressed; self-keyed it is squashed.
    const double keyed  = rms (500.0, 50.0);
    const double selfK  = rms (0.0,   50.0);
    ok (keyed > selfK * 1.5, "a 500 Hz key filter stops a 50 Hz tone from ducking the programme ("
                             + std::to_string (20.0 * std::log10 (keyed / std::max (1e-12, selfK))) + " dB louder)");

    // ...and it must NOT stop a tone above the corner from compressing.
    const double keyedHi = rms (500.0, 3000.0);
    const double selfHi  = rms (0.0,   3000.0);
    approx (keyedHi / selfHi, 1.0, 0.05, "above the corner the key filter changes almost nothing");

    // 0 Hz is not "a filter set flat" — it is the self-keyed three-argument path, and the chain must
    // take it. Nulled against a hand-built self-keyed run at the same quantum.
    {
        const int K = cfg.internalBlock;
        auto c2 = cfg;
        c2.sidechainHpfHz = 0.0;
        // The limiter has to come OUT of the topology here, not just be bypassed: a bypassed
        // latency-bearing stage still holds its PDC through the aligner, so leaving it in would shift
        // the chain against the hand-built reference by its 111 samples and this null would be
        // measuring the alignment rather than the key path. (It failed exactly that way when written.)
        c2.limiter = false;
        mastering::MasteringChain c;
        ok (c.prepare (48000.0, nch, c2), "prepare with the key filter off");
        c.setParams (p);
        Buf x = tone (nch, n, 180.0, 0.4f);
        Buf viaChain = x;
        { auto pl = planes (viaChain); c.process (pl.data(), nch, n); }

        dynamics::Compressor comp;
        ok (comp.prepare (48000.0, K, nch, std::max (cfg.compressorLookaheadMs, 1.0)), "hand: compressor prepare");
        { auto cp = p.compressor; cp.lookaheadMs = cfg.compressorLookaheadMs; comp.setParams (cp); }
        Buf byHand = x;
        for (int off = 0; off + K <= n; off += K)
        {
            float* ch[1] { byHand[0].data() + off };
            for (int i = 0; i < K; ++i)
            {
                const float v = ch[0][i];
                ch[0][i] = std::clamp (std::isfinite (v) ? v : 0.0f, -1.0e6f, 1.0e6f);
            }
            comp.process (ch, nch, K);                     // the THREE-argument, self-keyed form
        }
        long long bad = 0;
        for (int i = 0; i + K < n; ++i)
            if (bits (viaChain[0][(std::size_t) (i + K)]) != bits (byHand[0][(std::size_t) i])) ++bad;
        ok (bad == 0, "sidechainHpfHz = 0 is bit-identical to the self-keyed compressor (" + std::to_string (bad) + " differ)");
    }
}

//==============================================================================
int main()
{
    std::printf ("felitronics::mastering — oracles (nothing here trusts the chain's own arithmetic)\n");
    testLatencyFoundNotAsked();
    testActiveGroupDelayByPhase();
    testCompositionAgainstAHandBuiltChain();
    testPassbandPinned();
    testTwoGainNodes();
    testSidechainKey();
    return felitronics::test::report();
}

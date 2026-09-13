// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// felitronics_truepeak_instrument_gap_tests — HOW FAR APART THE CORE'S TWO TRUE-PEAK METERS READ (P62).
//
// `analysis::ReferenceTruePeakMeter` (4x, 32 taps per phase, 0.90 x Nyquist, at every rate — what certifies a
// delivered file) against `analysis::TruePeakMeter` (the spec's 12 taps per phase, full Nyquist, 4x / 2x / 1x
// by rate — what TargetLoudnessSolver used to aim with). A number about that gap used to travel between
// sessions with no measurement behind it; this suite is its owner now, and a header that mentions the gap
// points here instead of quoting it.
//
// WHAT IS MEASURED, and why it is not a single song. The quantity that matters is the gap on a DELIVERED
// master, so every material goes the way a delivery goes: made at 48 kHz, converted FIRST by
// `mastering::DeliveryConverter` to each of the six delivery rates, rendered through the default
// `MasteringChain` with 6 dB of drive into the limiter at a -1 dBTP ceiling, and then read by both meters as
// their callers read them (the reference drained by its own drain(), the cheap one with the 64 zeros the solver
// fed it). The corpus is synthetic because this repository carries no audio, and it is chosen by mechanism
// rather than by taste — the crew round named each mechanism before the fixture was built:
//   * music         tones plus beat-gated noise plus sparse transients: the calm case;
//   * drums         kick, snare, and hats made of FIRST-DIFFERENCED noise — bright, dense transients;
//   * bright-noise  flat white noise, hard limited;
//   * hf-burst      one 16 kHz burst, the worst of eight phases: the interpolation-GRID mechanism, which is
//                   what the cheap meter's 2x and 1x regimes at 88.2-192 kHz lose;
//   * click         one click flat to 0.45 fs, the worst of eight sub-sample offsets: the FILTER mechanism.
// The same comparison was run outside the tree on fourteen real programmes delivered at all six rates; at every
// rate its largest gap stayed below this table's largest. Its numbers belong to that run and are not copied here.
//
// THE TOLERANCE is 5e-4 dB per cell. The readings are not bit-exact across rows — the chain and the generators
// go through libm — but rounding moves them by orders of magnitude less than that.

#include <felitronics/analysis/ReferenceTruePeakMeter.h>
#include <felitronics/analysis/TruePeakMeter.h>
#include <felitronics/mastering/DeliveryConverter.h>
#include <felitronics/mastering/MasteringChain.h>
#include <felitronics/mastering/OfflineRenderer.h>
#include <felitronics_test.h>
#include <truepeak_oracle.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

namespace
{
namespace fm = felitronics::mastering;
namespace fa = felitronics::analysis;
using felitronics::test::ok;
using felitronics::test::group;

constexpr double kPi       = 3.14159265358979323846;
constexpr double kSourceFs = 48000.0;
constexpr int    kNch      = 2;
constexpr double kRates[]  = { 44100.0, 48000.0, 88200.0, 96000.0, 176400.0, 192000.0 };

using Planes = std::vector<std::vector<float>>;

struct Lcg { std::uint32_t s; double next() { s = s * 1664525u + 1013904223u; return (double) (std::int32_t) s * 4.656612873077393e-10; } };

double taper (int i, int n, int fade)
{
    const int k = std::min (i, n - 1 - i);
    return k >= fade ? 1.0 : 0.5 * (1.0 - std::cos (kPi * (double) k / (double) fade));
}

void normalise (Planes& p, double peakLin)
{
    double m = 0.0;
    for (auto& c : p) for (float v : c) m = std::max (m, (double) std::fabs (v));
    if (m > 0.0) for (auto& c : p) for (float& v : c) v = (float) ((double) v * peakLin / m);
}

// ---- the materials, all at the 48 kHz source ----
Planes music (double seconds)
{
    const int n = (int) (seconds * kSourceFs);
    Planes p (kNch, std::vector<float> ((std::size_t) n, 0.0f));
    Lcg r { 12345u };
    double lp[2] = { 0.0, 0.0 };
    for (int i = 0; i < n; ++i)
    {
        const double t = (double) i / kSourceFs;
        const double beat = std::exp (-8.0 * std::fmod (t, 0.5));
        for (int c = 0; c < kNch; ++c)
        {
            lp[c] = 0.90 * lp[c] + 0.10 * r.next();
            const double tone = std::sin (2.0 * kPi * (110.0 + 3.0 * c) * t) + 0.55 * std::sin (2.0 * kPi * (735.0 + 7.0 * c) * t)
                              + 0.30 * std::sin (2.0 * kPi * (2810.0 + 11.0 * c) * t);
            double v = 0.45 * tone + 0.9 * beat * lp[c];
            if (i % (int) (0.12 * kSourceFs) == 0) v += 1.6 * beat;
            p[(std::size_t) c][(std::size_t) i] = (float) (v * taper (i, n, 480));
        }
    }
    normalise (p, 0.9);
    return p;
}

Planes drums (double seconds)
{
    const int n = (int) (seconds * kSourceFs);
    Planes p (kNch, std::vector<float> ((std::size_t) n, 0.0f));
    Lcg r { 777u };
    double prev[2] = { 0.0, 0.0 };
    for (int i = 0; i < n; ++i)
    {
        const double t = (double) i / kSourceFs;
        const double tk = std::fmod (t, 0.5), ts = std::fmod (t + 0.25, 0.5), th = std::fmod (t, 0.125);
        const double kick = std::exp (-tk / 0.08) * std::sin (2.0 * kPi * (50.0 * tk + 40.0 * (1.0 - std::exp (-tk / 0.02)) * 0.02));
        for (int c = 0; c < kNch; ++c)
        {
            const double w = r.next();
            const double snare = std::exp (-ts / 0.06) * (0.7 * w + 0.3 * std::sin (2.0 * kPi * 190.0 * ts));
            const double hat   = std::exp (-th / 0.015) * (w - prev[c]);             // first difference: bright
            prev[c] = w;
            p[(std::size_t) c][(std::size_t) i] = (float) ((kick + 0.8 * snare + 0.6 * hat) * taper (i, n, 480));
        }
    }
    normalise (p, 0.9);
    return p;
}

Planes brightNoise (double seconds)
{
    const int n = (int) (seconds * kSourceFs);
    Planes p (kNch, std::vector<float> ((std::size_t) n, 0.0f));
    Lcg r { 4242u };
    for (int i = 0; i < n; ++i)
        for (int c = 0; c < kNch; ++c)
            p[(std::size_t) c][(std::size_t) i] = (float) (r.next() * taper (i, n, 480));
    normalise (p, 0.9);
    return p;
}

// One 16 kHz burst, 60 ms with 15 ms raised-cosine ends, both channels, starting phase `phase`.
Planes hfBurst (double phase)
{
    const int n = (int) (0.25 * kSourceFs);
    Planes p (kNch, std::vector<float> ((std::size_t) n, 0.0f));
    const int b0 = (int) (0.08 * kSourceFs), bl = (int) (0.06 * kSourceFs);
    for (int i = 0; i < bl; ++i)
    {
        const double v = std::sin (2.0 * kPi * 16000.0 * (double) i / kSourceFs + phase) * taper (i, bl, 720);
        for (int c = 0; c < kNch; ++c) p[(std::size_t) c][(std::size_t) (b0 + i)] = (float) (0.9 * v);
    }
    return p;
}

// One band-limited click flat to 0.45 fs of the source, Blackman over +-256 samples, centred `delta` samples
// after a grid point.
Planes click (double delta)
{
    const int n = (int) (0.25 * kSourceFs);
    Planes p (kNch, std::vector<float> ((std::size_t) n, 0.0f));
    const int centre = n / 2;
    for (int k = -256; k <= 256; ++k)
    {
        const double t = (double) k - delta;
        const double z = 0.9 * t;
        const double s = std::fabs (z) < 1e-12 ? 1.0 : std::sin (kPi * z) / (kPi * z);
        const double w = std::fabs (t) >= 256.0 ? 0.0 : 0.42 + 0.5 * std::cos (kPi * t / 256.0) + 0.08 * std::cos (2.0 * kPi * t / 256.0);
        for (int c = 0; c < kNch; ++c) p[(std::size_t) c][(std::size_t) (centre + k)] = (float) (0.9 * s * w);
    }
    return p;
}

// ---- delivery: SRC first, then the chain at the delivery rate with the limiter working ----
Planes deliver (const Planes& src, double rate)
{
    const long long inF = (long long) src[0].size();
    const long long outF = fm::DeliveryConverter::deliveredFrames (kSourceFs, rate, inF);
    Planes conv (kNch, std::vector<float> ((std::size_t) outF, 0.0f));
    fm::DeliveryConverter dc;
    felitronics::test::run (dc.prepare (kSourceFs, rate, kNch, 4096));
    const float* ip[kNch] { src[0].data(), src[1].data() };
    float*       cp[kNch] { conv[0].data(), conv[1].data() };
    felitronics::test::run (dc.convert (ip, kNch, inF, cp, outF));

    fm::MasteringChain chain;
    fm::OfflineRenderer renderer;
    felitronics::test::run (chain.prepare (rate, kNch, fm::MasteringChainConfig {}));
    felitronics::test::run (renderer.prepare (kNch, 4096));
    fm::MasteringChainParams params;
    params.compressor.thresholdDb = -20.0; params.compressor.ratio = 2.0; params.compressor.kneeDb = 6.0;
    params.compressor.attackMs = 15.0; params.compressor.releaseMs = 180.0;
    params.preLimiterGainDb = 6.0;
    params.limiter.ceilingDbTp = -1.0; params.limiter.releaseMs = 100.0;
    chain.setParams (params);
    Planes out (kNch, std::vector<float> ((std::size_t) outF, 0.0f));
    const float* ci[kNch] { conv[0].data(), conv[1].data() };
    float*       co[kNch] { out[0].data(), out[1].data() };
    felitronics::test::run (renderer.render (chain, ci, co, kNch, (int) outF));
    return out;
}

struct Reading { double r = 0.0, c = 0.0, cLin = 0.0, samplePeakLin = 0.0; int cFactor = 0; };

// Both instruments, each driven the way its certifying / aiming caller drives it: the reference drained by its
// own drain(), the cheap meter drained with the 64 zeros TargetLoudnessSolver used to feed it.
Reading measure (const Planes& x, double rate)
{
    const int n = (int) x[0].size();
    const float* p[kNch] { x[0].data(), x[1].data() };
    Reading out;
    fa::ReferenceTruePeakMeter R;
    felitronics::test::run (R.prepare (rate, n, kNch));
    felitronics::test::run (R.process (p, kNch, n));
    R.drain();
    fa::TruePeakMeter C;
    felitronics::test::run (C.prepare (rate, n, kNch));
    felitronics::test::run (C.process (p, kNch, n));
    std::vector<float> z (64, 0.0f);
    const float* zp[kNch] { z.data(), z.data() };
    felitronics::test::run (C.process (zp, kNch, 64));
    out.r = R.truePeakDb(); out.c = C.truePeakDb(); out.cLin = C.truePeakLinear(); out.samplePeakLin = R.samplePeakLinear();
    out.cFactor = C.oversampleFactor();
    return out;
}

// THE PINS. Every cell is R - C in dB (reference minus cheap) on the named material delivered at the named rate;
// the witness rows are the WORST over their grid offsets. Written from this suite's own printout, so a filter,
// a factor rule, the chain or the delivery converter that moves any reading fails here first — which is the
// point: this table is the only place the number lives.
constexpr int kMaterials = 5;
constexpr const char* kNames[kMaterials] = { "music", "drums", "bright-noise", "hf-burst", "click" };
constexpr double kPinned[kMaterials][6] = {
    //   44.1         48          88.2        96          176.4       192
    { +0.000725,  -0.006831,  +0.000031,  +0.000000,  +0.003574,  +0.000000 },   // music
    { +0.299474,  +0.136137,  -0.063778,  -0.006664,  +0.000295,  +0.000000 },   // drums
    { +0.000000,  +0.012345,  -0.006144,  +0.000000,  +0.000000,  +0.000000 },   // bright-noise
    { +0.002190,  +0.039063,  +0.009563,  +0.226378,  +0.010992,  +0.282293 },   // hf-burst (worst of 8 phases)
    { +0.172170,  +0.189417,  +0.159602,  +0.136960,  +0.179681,  +0.170299 },   // click    (worst of 8 offsets)
};
constexpr double kTolDb = 5.0e-4;
constexpr int    kWorstMaterial = 1, kWorstRate = 0;                  // drums @ 44.1 kHz
constexpr double kWorstTruthMinusReference = 0.025064;   // -0.637080 truth, -0.662144 reference
constexpr double kWorstTruthMinusCheap     = 0.324539;   //                 -0.961619 cheap
constexpr int kOffsets = 8;

Planes makeMaterial (int m, int offset)
{
    switch (m)
    {
        case 0:  return music (1.5);
        case 1:  return drums (1.5);
        case 2:  return brightNoise (1.5);
        case 3:  return hfBurst (2.0 * kPi * (double) offset / (double) kOffsets);
        default: return click ((double) offset / (double) kOffsets);
    }
}
bool searched (int m) { return m >= 3; }

} // namespace

int main()
{
    std::printf ("felitronics_truepeak_instrument_gap_tests\n");

    double gap[kMaterials][6] {};
    Reading at[kMaterials][6] {};
    std::printf ("  R - C, dB          ");
    for (double r : kRates) std::printf ("%11.1f", r / 1000.0);
    std::printf ("\n");
    for (int m = 0; m < kMaterials; ++m)
    {
        std::printf ("  %-18s ", kNames[m]);
        for (int k = 0; k < 6; ++k)
        {
            double worst = -1e9;
            for (int o = 0; o < (searched (m) ? kOffsets : 1); ++o)
            {
                const Reading rd = measure (deliver (makeMaterial (m, o), kRates[k]), kRates[k]);
                if (rd.r - rd.c > worst) { worst = rd.r - rd.c; at[m][k] = rd; }
            }
            gap[m][k] = worst;
            std::printf ("%+11.6f", worst);
        }
        std::printf ("\n");
    }

    group ("every cell of the table is the pinned number");
    for (int m = 0; m < kMaterials; ++m)
        for (int k = 0; k < 6; ++k)
            felitronics::test::approx (gap[m][k], kPinned[m][k], kTolDb,
                                       std::string (kNames[m]) + " @ " + std::to_string ((int) kRates[k]));

    group ("the worst cell is named: the material and the rate, not only the number");
    {
        int wm = 0, wk = 0;
        for (int m = 0; m < kMaterials; ++m)
            for (int k = 0; k < 6; ++k)
                if (gap[m][k] > gap[wm][wk]) { wm = m; wk = k; }
        std::printf ("    worst: %s @ %.1f kHz, R - C = %+.6f dB\n", kNames[wm], kRates[wk] / 1000.0, gap[wm][wk]);
        ok (wm == kWorstMaterial && wk == kWorstRate, "the worst reading gap is drums delivered at 44.1 kHz");
    }

    // The two rows a closed form owns, so they do not rest on the pins alone. A tone at fs/q with its crest
    // half-way between two points of an evaluation grid is under-read by cos(2*pi*s/q), s the crest's distance to
    // the nearest point in samples. The reference evaluates the odd eighths of a sample at every rate (a 128-tap
    // prototype's half-integer group delay), so its miss is at most 1/8. The cheap meter at 2x (88.2-176.4 kHz)
    // evaluates the quarters and the grid itself — a miss of 1/4 — and at 1x (>= 176.4 kHz) only the grid, a miss
    // of 1/2. 16 kHz is fs/6 at 96 kHz and fs/12 at 192 kHz; the search over eight burst phases lands near the
    // worst placement without being able to exceed it.
    group ("the grid-geometry rows are the closed form, from below");
    {
        const auto miss = [] (double s, double q) { return -20.0 * std::log10 (std::cos (2.0 * kPi * s / q)); };
        const double bound96  = miss (0.25, 6.0)  - miss (0.125, 6.0);
        const double bound192 = miss (0.5, 12.0)  - miss (0.125, 12.0);
        std::printf ("    hf-burst @ 96 kHz  %+.6f against %+.6f;  @ 192 kHz  %+.6f against %+.6f\n",
                     gap[3][3], bound96, gap[3][5], bound192);
        ok (gap[3][3] <= bound96 + 1.0e-3 && gap[3][3] >= bound96 - 2.0e-3, "96 kHz: the cheap meter's quarter-sample grid");
        ok (gap[3][5] <= bound192 + 1.0e-3 && gap[3][5] >= bound192 - 2.0e-3, "192 kHz: the cheap meter is the sample grid alone");
    }

    group ("at 176.4 kHz and above the cheap meter IS a sample-peak meter, on every material");
    for (int m = 0; m < kMaterials; ++m)
        for (int k = 4; k < 6; ++k)
            ok (at[m][k].cFactor == 1 && at[m][k].cLin == at[m][k].samplePeakLin,
                std::string (kNames[m]) + " @ " + std::to_string ((int) kRates[k]) + ": cheap reading == sample peak");

    // WHY THE CHEAP METER CANNOT AIM A DELIVERED PROMISE, and why the reference is not the last word either. The
    // worst cell against the band-limited truth, by spectral zero-padding (test_support/truepeak_oracle.h) and by
    // its sinc cross-check, which share no construction: the drums at 44.1 kHz leave the limiter ABOVE its own
    // -1 dBTP ceiling (its downsampler's step overshoot on bright, dense material — truepeak_witnesses.h), the
    // cheap meter under-reads that by more than TargetLoudnessSolver's whole 0.05 dB aim, and the reference
    // under-reads it too, by less.
    group ("the worst cell against the truth");
    {
        const Planes d = deliver (makeMaterial (1, 0), 44100.0);
        const Reading rd = measure (d, 44100.0);
        const double fft  = std::max (felitronics::test::tp::truePeakDbFft (d[0], 16), felitronics::test::tp::truePeakDbFft (d[1], 16));
        const double sinc = std::max (felitronics::test::tp::truePeakDbSinc (d[0], 256, 16, 2048), felitronics::test::tp::truePeakDbSinc (d[1], 256, 16, 2048));
        std::printf ("    drums @ 44.1 kHz: truth %+.6f (sinc %+.6f)  reference %+.6f  cheap %+.6f dBTP\n", fft, sinc, rd.r, rd.c);
        felitronics::test::approx (sinc, fft, 5.0e-3, "the two oracles agree");
        felitronics::test::approx (fft - rd.r, kWorstTruthMinusReference, kTolDb, "the reference under-reads the truth by the pinned amount");
        felitronics::test::approx (fft - rd.c, kWorstTruthMinusCheap, kTolDb, "the cheap meter under-reads it by the pinned amount");
        ok (fft - rd.c > 0.05, "which is more than the solver's 0.05 dB aim");
    }

    return felitronics::test::report();
}

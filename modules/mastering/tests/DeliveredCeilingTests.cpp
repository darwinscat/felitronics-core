// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// felitronics_delivered_ceiling_tests — A DELIVERED FILE MEETS ITS PROMISE AS THE CERTIFICATE READS IT, AT EVERY
// PAIR OF RATES (P62).
//
// The promise `maxTruePeakDbTp` is about the file `DeliveredMastering::solve` hands back, and the file is certified
// by `analysis::ReferenceTruePeakMeter` — the arithmetic `fcore_measure` runs (tools/fcore_probe.h measures through
// that class; felitronics_probe_tests nulls it against the hand-rolled loop). So the check is literal: for every
// source rate and every delivery rate of the six, solve a programme, read the delivered samples with the certifying
// class, and require
//   * the certificate at or under the promise;
//   * the solver's own report EQUAL to the certificate, bit for bit — the report is the certificate, not an estimate
//     of it (before P62 the solver read with the other meter and the two were different numbers). In dB, and so for
//     a delivered peak above the dB floors: the solver spells anything under 1e-10 as -200, the meter's truePeakDb()
//     floors at gainToDb's 1e-12 and fcore_measure's at 1e-9 — three spellings of silence, none of which is a file a
//     loudness target delivers;
//   * the loudness inside the request's tolerance.
//
// WHY THIS FIXTURE CAN FAIL. The material is the bright-transient programme felitronics_truepeak_instrument_gap_tests
// names as the worst reading gap between the two meters. The suite prints what `analysis::TruePeakMeter` — the
// meter the solver aimed with before P62 — reads on each delivered file, and requires at least one pair where it
// reads under the certificate by more than the request's whole aim margin: on that pair a solver aiming with it
// would have called a file feasible that the certificate puts over the promise. (Checked the direct way too, outside
// the suite: with the solver's measurement put back on the cheap meter, this file fails — pairs deliver over the
// promise, and on most of the others the report is no longer the certificate.)
//
// THE FIXTURE'S LENGTH AND PASS BUDGET ARE FOR THE SEARCH, not for the ceiling: 1.2 s at four passes left this bright,
// dense material in `PassLimit` on most pairs, a fact about gating blocks on a short programme and nothing about the
// meter. At 3 s and eight passes every pair solves, so every row is a delivered file and not a best effort.

#include <felitronics/analysis/ReferenceTruePeakMeter.h>
#include <felitronics/analysis/TruePeakMeter.h>
#include <felitronics/mastering/DeliveredMastering.h>
#include <felitronics_test.h>

#include <algorithm>
#include <bit>
#include <cmath>
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

constexpr double kPi      = 3.14159265358979323846;
constexpr int    kNch     = 2;
constexpr double kRates[] = { 44100.0, 48000.0, 88200.0, 96000.0, 176400.0, 192000.0 };
constexpr int    kBlock   = 4096;

using Planes = std::vector<std::vector<float>>;
struct Lcg { std::uint32_t s; double next() { s = s * 1664525u + 1013904223u; return (double) (std::int32_t) s * 4.656612873077393e-10; } };

// The gap suite's "drums", generated at `fs`: kick, snare, and hats of first-differenced noise, tapered ends.
Planes drums (double fs, double seconds)
{
    const int n = (int) (seconds * fs), fade = (int) (0.01 * fs);
    Planes p (kNch, std::vector<float> ((std::size_t) n, 0.0f));
    Lcg r { 777u };
    double prev[2] = { 0.0, 0.0 }, peak = 0.0;
    for (int i = 0; i < n; ++i)
    {
        const double t = (double) i / fs;
        const double tk = std::fmod (t, 0.5), ts = std::fmod (t + 0.25, 0.5), th = std::fmod (t, 0.125);
        const double kick = std::exp (-tk / 0.08) * std::sin (2.0 * kPi * (50.0 * tk + 40.0 * (1.0 - std::exp (-tk / 0.02)) * 0.02));
        const int k = std::min (i, n - 1 - i);
        const double env = k >= fade ? 1.0 : 0.5 * (1.0 - std::cos (kPi * (double) k / (double) fade));
        for (int c = 0; c < kNch; ++c)
        {
            const double w = r.next();
            const double snare = std::exp (-ts / 0.06) * (0.7 * w + 0.3 * std::sin (2.0 * kPi * 190.0 * ts));
            const double hat   = std::exp (-th / 0.015) * (w - prev[c]);
            prev[c] = w;
            const double v = (kick + 0.8 * snare + 0.6 * hat) * env;
            p[(std::size_t) c][(std::size_t) i] = (float) v;
            peak = std::max (peak, std::fabs (v));
        }
    }
    for (auto& c : p) for (float& v : c) v = (float) ((double) v * 0.5 / peak);
    return p;
}

struct Row
{
    fm::MasteringSolveStatus status = fm::MasteringSolveStatus::NotPrepared;
    double integrated = 0.0, reported = 0.0, certificate = 0.0, cheap = 0.0;
    int passes = 0;
};

Row deliver (double source, double delivery, const fm::LoudnessRequest& req)
{
    Row row;
    const Planes in = drums (source, 3.0);
    const long long inF  = (long long) in[0].size();
    const long long outF = fm::DeliveredMastering::deliveredFrames (source, delivery, inF);
    Planes out (kNch, std::vector<float> ((std::size_t) outF, 0.0f));

    fm::MasteringChain chain;
    fm::OfflineRenderer renderer;
    fm::TargetLoudnessSolver solver;
    fm::DeliveredMastering dm;
    if (! felitronics::test::run (chain.prepare (delivery, kNch, fm::MasteringChainConfig {}))) return row;
    if (! felitronics::test::run (renderer.prepare (kNch, kBlock))) return row;
    if (! felitronics::test::run (solver.prepare (delivery, kNch, kBlock, chain.internalBlock(), chain.tapOversampleFactor()))) return row;
    if (! felitronics::test::run (dm.prepare (source, delivery, kNch, kBlock))) return row;

    fm::MasteringChainParams params;
    params.compressor.thresholdDb = -20.0; params.compressor.ratio = 2.0; params.compressor.kneeDb = 6.0;
    params.compressor.attackMs = 15.0; params.compressor.releaseMs = 180.0;
    params.limiter.ceilingDbTp = -1.0; params.limiter.releaseMs = 100.0;

    const float* ip[kNch] { in[0].data(), in[1].data() };
    float*       op[kNch] { out[0].data(), out[1].data() };
    const fm::LoudnessSolution sol = dm.solve (solver, chain, renderer, params, ip, kNch, inF, op, outF, req);
    row.status = sol.status; row.passes = sol.passes;
    row.integrated = sol.measured.integratedLufs; row.reported = sol.measured.truePeakDbTp;

    const float* dp[kNch] { out[0].data(), out[1].data() };
    fa::ReferenceTruePeakMeter cert;
    if (! felitronics::test::run (cert.prepare (delivery, (int) outF, kNch))) return row;
    if (! felitronics::test::run (cert.process (dp, kNch, (int) outF))) return row;
    cert.drain();
    row.certificate = cert.truePeakDb();

    fa::TruePeakMeter cheap;                                           // what the solver used to aim with, as it fed it
    if (! felitronics::test::run (cheap.prepare (delivery, (int) outF, kNch))) return row;
    if (! felitronics::test::run (cheap.process (dp, kNch, (int) outF))) return row;
    const std::vector<float> z (64, 0.0f);
    const float* zp[kNch] { z.data(), z.data() };
    if (! felitronics::test::run (cheap.process (zp, kNch, 64))) return row;
    row.cheap = cheap.truePeakDb();
    return row;
}

} // namespace

int main()
{
    std::printf ("felitronics_delivered_ceiling_tests\n");
    fm::LoudnessRequest req;
    req.targetLufs = -12.0; req.maxTruePeakDbTp = -1.0; req.toleranceLu = 0.1; req.maxPasses = 8;

    group ("every source rate to every delivery rate: under the promise, as the certificate reads it");
    int breaches = 0, pairs = 0;
    double worstCheapGap = -1e9;
    std::string worstPair;
    std::printf ("      source -> delivery   status  passes   I (LUFS)   certificate  reported     cheap\n");
    for (double s : kRates)
        for (double d : kRates)
        {
            const Row r = deliver (s, d, req);
            ++pairs;
            const std::string pair = std::to_string ((int) s) + " -> " + std::to_string ((int) d);
            std::printf ("      %6d -> %6d   %6d  %6d   %8.4f   %+11.6f  %+.6f  %+.6f\n", (int) s, (int) d, (int) r.status,
                         r.passes, r.integrated, r.certificate, r.reported, r.cheap);
            ok (r.status == fm::MasteringSolveStatus::Solved, pair + ": solved");
            ok (r.certificate <= req.maxTruePeakDbTp, pair + ": the certificate is at or under the promise");
            ok (std::bit_cast<std::uint64_t> (r.reported) == std::bit_cast<std::uint64_t> (r.certificate),
                pair + ": the solver's report IS the certificate, bit for bit");
            ok (std::fabs (r.integrated - req.targetLufs) <= req.toleranceLu, pair + ": the loudness is inside the tolerance");
            if (r.certificate > req.maxTruePeakDbTp) ++breaches;
            if (r.certificate - r.cheap > worstCheapGap) { worstCheapGap = r.certificate - r.cheap; worstPair = pair; }
        }
    std::printf ("      %d pairs, %d over the promise; the cheap meter reads under the certificate by up to %+.6f dB (%s)\n",
                 pairs, breaches, worstCheapGap, worstPair.c_str());

    group ("the fixture can fail: aiming with the cheap meter would have misjudged at least one pair");
    ok (worstCheapGap > req.truePeakAimDb,
        "on " + worstPair + " the cheap meter reads under the certificate by more than the whole aim margin");

    return felitronics::test::report();
}

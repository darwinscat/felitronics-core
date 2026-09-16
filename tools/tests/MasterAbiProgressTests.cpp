// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#include <felitronics_test.h>

#include "fc_master_abi.h"

#include <cmath>
#include <cstdint>
#include <cstring>
#include <random>
#include <string>
#include <vector>

using felitronics::test::group;
using felitronics::test::ok;

namespace
{
constexpr double   kFs = 48000.0;
constexpr int      kNch = 2;
constexpr uint32_t kFrames = 4u * 48000u;

std::vector<float> programme (double fs, uint32_t frames)
{
    std::vector<float> v ((std::size_t) frames * kNch);
    std::mt19937 rng (99u);
    std::uniform_real_distribution<float> u (-1.0f, 1.0f);
    for (int c = 0; c < kNch; ++c)
        for (uint32_t i = 0; i < frames; ++i)
        {
            const double t = (double) i / fs;
            const double env = 0.4 + 0.6 * (0.5 + 0.5 * std::sin (2.0 * 3.141592653589793 * 0.7 * t));
            double x = env * (0.5 * std::sin (2.0 * 3.141592653589793 * (220.0 + 5.0 * c) * t) + 0.2 * u (rng));
            if (i % 5760u == 0u) x += 0.6;
            v[(std::size_t) c * frames + i] = (float) (0.3 * x);
        }
    return v;
}

fc_master makeHandle (double deliveryRate = 0.0)
{
    fc_master_config c {};
    FC_INIT (c);
    if (fc_master_config_defaults (&c) != FC_OK) return 0;
    c.sampleRate = kFs; c.channels = kNch; c.deliveryRate = deliveryRate;
    fc_master h = 0;
    return fc_master_create (&c, &h) == FC_OK ? h : 0;
}

struct Inputs
{
    fc_master_params    p {};
    fc_loudness_request r {};
    bool ok = false;
    Inputs()
    {
        FC_INIT (p); FC_INIT (r);
        ok = fc_master_params_defaults (&p) == FC_OK && fc_loudness_request_defaults (&r) == FC_OK;
        r.targetLufs = -14.0; r.maxTruePeakDbTp = -1.0;
    }
};

struct Recorder
{
    std::vector<fc_progress> seen;
    long long stopAt = -1;
    static int32_t fn (void* ctx, const fc_progress* e)
    {
        auto& r = *static_cast<Recorder*> (ctx);
        r.seen.push_back (*e);
        return (long long) r.seen.size() - 1 == r.stopAt ? 0 : 1;
    }
};

bool sameBytes (const void* a, const void* b, std::size_t n) { return std::memcmp (a, b, n) == 0; }
bool samePass (const fc_solve_pass& a, const fc_solve_pass& b)
{
    return sameBytes (&a.gainDb, &b.gainDb, 8) && sameBytes (&a.ceilingDb, &b.ceilingDb, 8)
        && sameBytes (&a.integratedLufs, &b.integratedLufs, 8) && sameBytes (&a.truePeakDbTp, &b.truePeakDbTp, 8)
        && sameBytes (&a.plrDb, &b.plrDb, 8) && sameBytes (&a.limiterMaxGrDb, &b.limiterMaxGrDb, 8)
        && sameBytes (&a.loudnessRangeLu, &b.loudnessRangeLu, 8) && a.violated == b.violated;
}

struct Solved
{
    fc_status st = FC_ERR_HANDLE;
    fc_solution sol = 0xDEADBEEFu;
    std::vector<float> out;
    fc_solution_summary sum {};
    std::vector<fc_solve_pass> log;
};

Solved solve (fc_master h, const Inputs& in, const std::vector<float>& audio)
{
    Solved s;
    s.out.assign (audio.size(), 0.25f);
    s.st = fc_master_solve (h, &in.p, &in.r, audio.data(), s.out.data(), kFrames, &s.sol);
    if (s.st != FC_OK) return s;
    FC_INIT (s.sum);
    (void) fc_solution_summary_get (s.sol, &s.sum);
    s.log.resize (32);
    uint32_t w = 0;
    (void) fc_solution_log (s.sol, s.log.data(), 32u, &w);
    s.log.resize (w);
    (void) fc_solution_destroy (s.sol);
    return s;
}

bool sameSolve (const Solved& a, const Solved& b)
{
    bool same = a.st == FC_OK && b.st == FC_OK && a.out.size() == b.out.size()
             && sameBytes (a.out.data(), b.out.data(), a.out.size() * sizeof (float))
             && a.sum.status == b.sum.status && a.sum.passes == b.sum.passes && a.sum.logCount == b.sum.logCount
             && sameBytes (&a.sum.preLimiterGainDb, &b.sum.preLimiterGainDb, 8)
             && sameBytes (&a.sum.ceilingDbTp, &b.sum.ceilingDbTp, 8) && a.log.size() == b.log.size();
    for (std::size_t i = 0; same && i < a.log.size(); ++i) same = samePass (a.log[i], b.log[i]);
    return same;
}
}   // namespace

static void testSetProgress()
{
    group ("fc_master_set_progress: the handle, and removal");
    Recorder r;
    ok (fc_master_set_progress (0u, &Recorder::fn, &r) == FC_ERR_HANDLE, "handle 0: FC_ERR_HANDLE");
    const fc_master h = makeHandle();
    ok (h != 0u && fc_master_set_progress (h, &Recorder::fn, &r) == FC_OK, "a live handle: FC_OK");
    ok (fc_master_set_progress (h, nullptr, nullptr) == FC_OK, "fn NULL: FC_OK");
    const Inputs in;
    const auto audio = programme (kFs, kFrames);
    const Solved s = solve (h, in, audio);
    ok (in.ok && s.st == FC_OK && r.seen.empty(), "a removed callback is not called");
    (void) fc_master_destroy (h);
    ok (fc_master_set_progress (h, &Recorder::fn, &r) == FC_ERR_HANDLE, "a destroyed handle: FC_ERR_HANDLE");
}

static void testSolve()
{
    group ("fc_master_solve with a callback: the same bits, the stages, the records");
    const Inputs in;
    const auto audio = programme (kFs, kFrames);
    const fc_master a = makeHandle(), b = makeHandle(), other = makeHandle();
    Recorder r, bystander;
    ok (in.ok && fc_master_set_progress (b, &Recorder::fn, &r) == FC_OK
        && fc_master_set_progress (other, &Recorder::fn, &bystander) == FC_OK, "PRECONDITION: set");
    const Solved plain = solve (a, in, audio);
    const Solved heard = solve (b, in, audio);
    ok (plain.st == FC_OK && plain.sum.passes >= 2, "PRECONDITION: solved in " + std::to_string (plain.sum.passes) + " renders");
    ok (sameSolve (plain, heard), "audio, summary and log: the same bits with a callback");
    ok (bystander.seen.empty(), "another handle's callback is not called");

    int latency = 0;
    (void) fc_master_latency (b, &latency);
    const double units = 2.0 * kFrames + latency, bound = kFrames / 100.0;
    int stages = 0, shape = 0, gaps = 0, records = 0;
    for (std::size_t i = 0; i < r.seen.size(); ++i)
    {
        const fc_progress& e = r.seen[i];
        if (e.fraction == 0.0)
        {
            ++stages;
            if (e.stage != FC_PROGRESS_PASS || e.pass != stages || e.maxPasses != in.r.maxPasses) ++shape;
            if (i > 0 && r.seen[i - 1].fraction != 1.0) ++shape;
        }
        else
        {
            const fc_progress& p = r.seen[i - 1];
            if (e.stage != p.stage || e.pass != p.pass || ! (e.fraction > p.fraction)) ++shape;
            if ((e.fraction - p.fraction) * units > bound * (1.0 + 1e-9)) ++gaps;
        }
        if (e.hasRecord != (e.fraction == 1.0 ? 1 : 0)) ++shape;
        if (e.hasRecord == 1)
        {
            if (e.pass < 1 || (std::size_t) e.pass > heard.log.size() || ! samePass (e.record, heard.log[(std::size_t) e.pass - 1])) ++records;
        }
        else
        {
            const fc_solve_pass zero {};
            if (! samePass (e.record, zero)) ++records;
        }
    }
    ok (stages == heard.sum.passes && shape == 0 && ! r.seen.empty() && r.seen.back().fraction == 1.0,
        "one PASS stage per render, numbered from 1 against maxPasses, each from exactly 0 to exactly 1 ("
        + std::to_string (r.seen.size()) + " events)");
    ok (gaps == 0, "no two events more than 1 % of the programme apart");
    ok (records == 0, "each closing event carries its log entry, every other event a zero record");
    (void) fc_master_destroy (a); (void) fc_master_destroy (b); (void) fc_master_destroy (other);
}

static void testStop()
{
    group ("a stop: FC_ERR_CANCELLED, nothing after it, the handle's state, then the same bits");
    const Inputs in;
    const auto audio = programme (kFs, kFrames);
    const fc_master ref = makeHandle();
    const Solved plain = solve (ref, in, audio);
    Recorder all;
    const fc_master probe = makeHandle();
    (void) fc_master_set_progress (probe, &Recorder::fn, &all);
    (void) solve (probe, in, audio);
    ok (plain.st == FC_OK && all.seen.size() > 10, "PRECONDITION: a reference and its events");

    for (const long long at : { 0LL, (long long) all.seen.size() / 2, (long long) all.seen.size() - 1 })
    {
        const fc_master h = makeHandle();
        Recorder r; r.stopAt = at;
        (void) fc_master_set_progress (h, &Recorder::fn, &r);
        Solved s = solve (h, in, audio);
        const std::string where = "stop at event " + std::to_string (at);
        ok (s.st == FC_ERR_CANCELLED && s.sol == 0xDEADBEEFu && (long long) r.seen.size() == at + 1,
            where + ": FC_ERR_CANCELLED, *out_solution untouched, nothing after");
        std::vector<float> one (kNch, 0.0f), got (kNch, 0.0f);
        const fc_status proc = fc_master_process (h, one.data(), got.data(), 1u);
        ok (at == 0 ? proc == FC_OK : proc == FC_ERR_STATE,
            where + (at == 0 ? ": no render begun, the handle streams as before" : ": a render began, process is FC_ERR_STATE"));
        if (at == 0) (void) fc_master_reset (h);
        (void) fc_master_set_progress (h, nullptr, nullptr);
        ok (sameSolve (solve (h, in, audio), plain), where + ": the same handle then solves to the same bits");
        (void) fc_master_destroy (h);
    }
    (void) fc_master_destroy (ref); (void) fc_master_destroy (probe);
}

static void testRange()
{
    group ("fc_master_measure_lra: one LRA stage, the same range, a stop");
    const auto audio = programme (kFs, kFrames);
    const fc_master h = makeHandle();
    double plain = 0.0;
    ok (fc_master_measure_lra (h, audio.data(), kFrames, &plain) == FC_OK, "PRECONDITION: measured");
    Recorder r;
    (void) fc_master_set_progress (h, &Recorder::fn, &r);
    double heard = -1.0;
    ok (fc_master_measure_lra (h, audio.data(), kFrames, &heard) == FC_OK && sameBytes (&plain, &heard, 8),
        "the same range with a callback");
    int shape = 0;
    for (std::size_t i = 0; i < r.seen.size(); ++i)
        if (r.seen[i].stage != FC_PROGRESS_LRA || r.seen[i].pass != 0 || r.seen[i].maxPasses != 0 || r.seen[i].hasRecord != 0
            || (i > 0 && ! (r.seen[i].fraction > r.seen[i - 1].fraction))) ++shape;
    ok (r.seen.size() > 2 && shape == 0 && r.seen.front().fraction == 0.0 && r.seen.back().fraction == 1.0,
        "LRA, pass 0 of 0, no record, from exactly 0 to exactly 1");
    Recorder s; s.stopAt = (long long) r.seen.size() / 2;
    (void) fc_master_set_progress (h, &Recorder::fn, &s);
    double v = 123.0;
    ok (fc_master_measure_lra (h, audio.data(), kFrames, &v) == FC_ERR_CANCELLED && v == 123.0
        && (long long) s.seen.size() == s.stopAt + 1, "a stop: FC_ERR_CANCELLED, *out untouched, nothing after");
    double again = 0.0;
    ok (fc_master_measure_lra (h, audio.data(), kFrames, &again) == FC_OK && sameBytes (&plain, &again, 8), "then the same range");
    (void) fc_master_destroy (h);
}

static void testDelivered()
{
    group ("fc_master_solve_delivered 48 -> 44.1 kHz: CONVERT first, a stop there moves nothing");
    const Inputs in;
    const auto audio = programme (kFs, kFrames);
    const fc_master h = makeHandle (44100.0);
    uint32_t d = 0;
    ok (h != 0u && fc_master_delivered_frames (h, kFrames, &d) == FC_OK, "PRECONDITION: a delivering handle");
    Recorder r;
    (void) fc_master_set_progress (h, &Recorder::fn, &r);
    std::vector<float> out ((std::size_t) d * kNch);
    fc_solution sol = 0;
    ok (fc_master_solve_delivered (h, &in.p, &in.r, audio.data(), kFrames, out.data(), d, &sol) == FC_OK
        && ! r.seen.empty() && r.seen.front().stage == FC_PROGRESS_CONVERT && r.seen.back().stage == FC_PROGRESS_PASS,
        "CONVERT, then PASS");
    (void) fc_solution_destroy (sol);

    const fc_master g = makeHandle (44100.0);
    std::size_t convertEvents = 0;
    while (convertEvents < r.seen.size() && r.seen[convertEvents].stage == FC_PROGRESS_CONVERT) ++convertEvents;
    Recorder s; s.stopAt = (long long) convertEvents / 2;
    (void) fc_master_set_progress (g, &Recorder::fn, &s);
    fc_solution none = 0xDEADBEEFu;
    ok (fc_master_solve_delivered (g, &in.p, &in.r, audio.data(), kFrames, out.data(), d, &none) == FC_ERR_CANCELLED
        && none == 0xDEADBEEFu && s.seen.back().stage == FC_PROGRESS_CONVERT, "a stop in CONVERT: FC_ERR_CANCELLED");
    ok (fc_master_render_delivered (g, audio.data(), kFrames, out.data(), d) == FC_OK,
        "and the handle still renders without a configure");
    (void) fc_master_destroy (h); (void) fc_master_destroy (g);
}

int main()
{
    std::printf ("fc_master — progress and cancellation through the C ABI\n");
    testSetProgress();
    testSolve();
    testStop();
    testRange();
    testDelivered();
    return felitronics::test::report();
}

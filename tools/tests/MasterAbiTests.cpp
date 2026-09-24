// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// The mastering C ABI, exercised natively. `tools/wasm/fc_master.cpp` compiles anywhere (its
// EMSCRIPTEN_KEEPALIVE degrades to a plain extern "C"), so the whole validation, ownership and mapping
// layer runs under ctest, ASan and UBSan like anything else.
//
// WHAT THIS SUITE IS FOR, said plainly because it decides what belongs in it: every REFUSAL the ABI
// makes has a test here that fails if the refusal disappears. The facade's job is to say no to what a
// page can hand it, and a refusal with no test is a refusal that will be deleted by someone tidying up.
//
// What is NOT here: the bit-exactness of the render, which needs a whole programme and lives in
// `fcore_master selftest`; and the detached-view discipline, which needs a wasm heap and lives in the
// browser harness.

#include <felitronics_test.h>
#include <alloc_counter.h>   // installs the allocation counter: EVERY form of `new`, over-aligned included

#include "fc_master_abi.h"

#include <felitronics/eq/EqEngine.h>
#include <felitronics/mastering/DeliveryConverter.h>
#include <felitronics/mastering/LoudnessSolver.h>
#include <felitronics/mastering/MasteringChain.h>
#include <felitronics/mastering/OfflineRenderer.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <memory>
#include <new>
#include <string>
#include <tuple>
#include <vector>

// The allocation counter is the shared one (test_support/alloc_counter.h, included above): every form of
// `new`, the over-aligned included. This file used to carry its own copy of it — the reference copy, in
// fact, since it was already complete — and the counter's fault switch, its byte accounting and the
// plain-object exemption all moved into that header verbatim, so the other sixty suites get them too.

// One vector's allocation, as the counter sees it — written through `volatile`, so the optimizer cannot remove it.
static long long vectorRequest (std::size_t n)
{
    const long long before = alloc::bytes.load();
    {
        std::vector<char> v;
        v.assign (n, 0);
        volatile char* sink = v.data();
        sink[0] = 1;
    }
    return alloc::bytes.load() - before;
}

using felitronics::test::ok;
using felitronics::test::approx;
using felitronics::test::group;
using felitronics::test::okNoAlloc;

namespace
{

constexpr double kFs  = 48000.0;
constexpr int    kNch = 2;

// The two quantile histograms a solve leaves in its solution (`GainReductionSummariser`) — a ONE-TIME allocation of the
// first render, like the traces, and 0.01 dB over 400 dB on every row. Spelled through the histogram's own sizing
// function so a range or a bin width that moves moves this with it.
// THREE since K11: the compressor's distribution, the limiter's, and the limiter's gated one. The count is part
// of the oracle — adding the third turned every budget check in three suites red at once, which is the
// allocation oracle doing its job rather than an inconvenience to be edited away.
const std::uint64_t kGrWindowBytes =
    3u * felitronics::dynamics::offline::QuantileHistogram::storageBytes (
             0.0, felitronics::mastering::TargetLoudnessSolver::kGrRangeDb, 0.01);

// What a `LoudnessSolution` costs BESIDES the parts spelled out at its check: the passes' log, the measurement,
// the verdict fields. Named so that the parts that DO move — the traces, the histograms, K11's gated summary —
// are added one by one and a change to any of them names itself instead of being absorbed into a total.
// K13 (v11) took it from 2368 to 2416: the measurement gained three doubles and three int64s for what
// the peak clipper did. K14 (v12) took it to 2464 for the air band's five doubles and its sample count.
// Exactly 48 bytes each time, and the pin is here so that a growth nobody meant names itself.
constexpr std::uint64_t kSolutionRecordRest = 2464u;

// The topology axis of the P41 create/configure matrix — see the switch that reads it.
constexpr int kTopologies = 9;

fc_master_config goodConfig()
{
    fc_master_config c {};
    fc_master_config_default (&c);
    c.sampleRate = kFs;
    c.channels   = kNch;
    c.monoBass   = 1;         // stereo-only stage, present so its mapping is exercised
    c.clipper    = 1;
    return c;
}

fc_master_params goodParams()
{
    fc_master_params p {};
    fc_master_params_default (&p);
    return p;
}

// A handle that is definitely alive, or 0 if creation failed.
fc_master make()
{
    const fc_master_config c = goodConfig();
    fc_master h = 0;
    (void) fc_master_create (&c, &h);
    return h;
}

// EVERY CHANNEL DIFFERENT, and that is not decoration: with identical planes no fault in the planar
// addressing — a stride of zero, every plane pointed at channel 0, the two channels crossed — can
// change a single sample, so the whole suite would be blind to a whole class while looking busy.
std::vector<float> tone (std::size_t frames, int nch)
{
    std::vector<float> v (frames * (std::size_t) nch, 0.0f);
    for (int c = 0; c < nch; ++c)
        for (std::size_t i = 0; i < frames; ++i)
            v[(std::size_t) c * frames + i] =
                (float) ((0.4 - 0.11 * c) * std::sin (2.0 * 3.14159265358979 * (440.0 + 137.0 * c)
                                                      * (double) i / kFs));
    return v;
}

//==============================================================================
// v7 — THE EQ CURVE AND ITS ORACLE.
//
// The curve is held against AUDIO: the same bands are driven through `eq::EqEngine` and the steady-state
// gain of a sine is measured. The two sides share the numbers in `CurveBand` and nothing else — the ABI's
// C struct and the core's own `BandParams` are written from them separately, so the mapping the curve goes
// through is not the mapping the oracle goes through.
struct CurveBand
{
    int    lane;                    // the index into `fc_eq_band::lanes`, which is eq::Lane's order
    int    type;                    // fc_filter_type
    double freqHz, q, gainDb;
    int    slope;
};

void placeBand (fc_eq_band& fb, felitronics::eq::BandParams& cb, const CurveBand& s)
{
    fb.on = 1; fb.type = s.type; fb.swept = 0; fb.bypass = 0;
    for (int l = 0; l < FC_MAX_EQ_LANES; ++l) { fb.lanes[l].on = 0; fb.lanes[l].bypass = 0; }
    fb.lanes[s.lane].on     = 1;
    fb.lanes[s.lane].freq   = s.freqHz;
    fb.lanes[s.lane].q      = s.q;
    fb.lanes[s.lane].gainDb = s.gainDb;
    fb.lanes[s.lane].slope  = s.slope;

    cb = felitronics::eq::BandParams {};
    cb.on   = true;
    cb.type = (felitronics::eq::FilterType) s.type;
    for (int l = 0; l < felitronics::eq::kNumLanes; ++l) cb.lanes[l].on = false;
    cb.lanes[s.lane].on     = true;
    cb.lanes[s.lane].freq   = s.freqHz;
    cb.lanes[s.lane].Q      = s.q;
    cb.lanes[s.lane].gainDb = s.gainDb;
    cb.lanes[s.lane].slope  = s.slope;
}

// Which axis the measurement reads, and how the two channels are driven to isolate it. Mono runs a
// ONE-channel engine, where the band runs its Stereo lane and nothing else; Left/Right run two
// independent channels; Mid and Side drive L = R and L = −R, so the other domain is exactly zero.
enum class Drive { Mono, Left, Right, Mid, Side };

double domain (Drive d, double l, double r)
{
    switch (d)
    {
        case Drive::Mono:  case Drive::Left: return l;
        case Drive::Right: return r;
        case Drive::Mid:   return 0.5 * (l + r);
        case Drive::Side:  return 0.5 * (l - r);
    }
    return l;
}

// The measured magnitude, in dB, at each frequency of `freqs`. Every frequency is a multiple of fs/N, so
// N samples hold a whole number of periods and the ratio of the two sums is the steady-state gain.
std::vector<double> measuredCurve (const felitronics::eq::BandParams* bands, int n, double fs,
                                   const std::vector<double>& freqs, Drive d)
{
    namespace E = felitronics::eq;
    const int nch = (d == Drive::Mono) ? 1 : 2;
    const int N   = 4800;
    std::vector<double> out (freqs.size(), std::numeric_limits<double>::quiet_NaN());

    auto eng = std::make_unique<E::EqEngine>();
    if (! eng->prepare (fs, N, nch)) return out;
    for (int i = 0; i < E::EqEngine::kMaxBands; ++i)
        eng->setBand (i, i < n ? bands[i] : E::BandParams {});

    std::vector<float> L ((std::size_t) N), R ((std::size_t) N);
    for (std::size_t fi = 0; fi < freqs.size(); ++fi)
    {
        const double dp = 2.0 * felitronics::core::kPi * freqs[fi] / fs;
        auto fill = [&]
        {
            for (int k = 0; k < N; ++k)
            {
                const double v = 0.25 * std::sin (dp * (double) k);
                L[(std::size_t) k] = (float) v;
                R[(std::size_t) k] = (float) (d == Drive::Side ? -v : d == Drive::Mid ? v : 0.7 * v);
            }
        };
        float* ch[2] = { L.data(), R.data() };
        eng->reset();
        for (int b = 0; b < 5; ++b) { fill(); if (! eng->process (ch, nch, N)) return out; }

        fill();
        std::vector<double> before ((std::size_t) N);
        for (int k = 0; k < N; ++k) before[(std::size_t) k] = domain (d, L[(std::size_t) k], R[(std::size_t) k]);
        if (! eng->process (ch, nch, N)) return out;

        double si = 0.0, so = 0.0;
        for (int k = 0; k < N; ++k)
        {
            const double o = domain (d, L[(std::size_t) k], R[(std::size_t) k]);
            si += before[(std::size_t) k] * before[(std::size_t) k];
            so += o * o;
        }
        out[fi] = 10.0 * std::log10 (so / si);
    }
    return out;
}

}   // namespace

int main()
{
    std::printf ("fc_master — the mastering C ABI, exercised natively\n");

    //==========================================================================
    group ("build identity");
    {
        ok (fc_master_abi_version() == FC_MASTER_ABI_VERSION, "the version this build speaks");
        ok (fc_master_max_eq_bands() == FC_MAX_EQ_BANDS, "the band count the ABI mirrors");
        // FROZEN at v1, so a v1 loader fails on the version rather than on a size it would blame on its own layout.
        ok (fc_master_sizeof_params() == 6560u, "sizeof(params) at v1 — frozen");
        ok (fc_master_sizeof_config() == 80u, "sizeof(config) at v1 — frozen, though this build's is 88");
        ok (fc_master_sizeof (FC_STRUCT_PARAMS, FC_MASTER_ABI_VERSION) == sizeof (fc_master_params)
            && fc_master_sizeof (FC_STRUCT_CONFIG, FC_MASTER_ABI_VERSION) == sizeof (fc_master_config),
            "and the size table's current row is this build's sizeof");
        ok (fc_master_max_channels() == (std::uint32_t) felitronics::core::kMaxChannels,
            "and the channel ceiling is the CORE's, not a number this file chose");
    }

    //==========================================================================
    // THE DEFAULTS ARE THE CORE'S. A zeroed struct is not a valid parameter set and is not close to one:
    // measured, a memset(0) set against `MasteringChainParams{}` differs in 287 998 of 288 000 samples.
    // These checks pin the three fields that make that true, so a `fc_master_params_default` that
    // regressed to a memset would fail HERE rather than in someone's master.
    group ("fc_master_params_default carries the CORE's defaults, not zero");
    {
        fc_master_params p {};
        std::memset (&p, 0xAB, sizeof p);      // poison, so "written" is distinguishable from "left"
        fc_master_params_default (&p);
        ok (p.header.abiVersion == 1u && p.header.structSize == 6560u,
            "the header is stamped — at v1, where the frozen writer lives — so the result is usable as an argument");
        ok (p.dither.bits == 24, "dither defaults to 24 bits, not 0");
        ok (p.dither.shaping == FC_SHAPING_WEIGHTED, "and to weighted shaping, not None");
        ok ((p.dither.seedLo | p.dither.seedHi) != 0u, "and to a NON-ZERO seed");
        ok (p.eqBands[0].lanes[0].on == 1, "the Stereo lane is ON by default (a zeroed struct has it off)");
        ok (p.eqBands[0].lanes[1].on == 0, "and the other lanes are not");
        ok (p.monoBass.enabled == 1, "MonoBassParams::enabled defaults true — unlike the CONFIG's monoBass");
        ok (p.compressor.ratio > 1.0, "the compressor ratio is a real default, not 0");

        fc_master_config c {};
        std::memset (&c, 0xAB, sizeof c);
        fc_master_config_default (&c);
        ok (c.monoBass == 0, "the CONFIG's monoBass defaults OFF — the two defaults really do disagree");
        // The core has no default rate and no default width, so neither does this. Writing one would be
        // the facade choosing a geometry for every caller who forgot to state one.
        ok (c.sampleRate == 0.0 && c.channels == 0, "and it invents NO sample rate and NO channel count");
        fc_master hFromDefault = 0;
        ok (fc_master_create (&c, &hFromDefault) == FC_ERR_REFUSED_BY_CORE,
            "so a create straight from the defaults is refused rather than silently given 48 kHz stereo");
        ok (c.internalBlock == 256, "the internal quantum's default");
        ok (c.tapsPerPhase == 64, "and the oversampler's tap count");

        fc_loudness_request r {};
        std::memset (&r, 0xAB, sizeof r);
        fc_loudness_request_default (&r);
        ok (std::isnan (r.targetLufs) && std::isnan (r.maxTruePeakDbTp),
            "the two REQUIRED request fields stay NaN: the core ships no delivery policy and neither does this");
        ok (std::isinf (r.limiterGr.limitDb) && r.limiterGr.limitDb > 0.0,
            "an OFF gain-reduction limit is +infinity and nothing else");
    }

    //==========================================================================
    group ("create refuses what a page can hand it");
    {
        fc_master h = 1234;
        ok (fc_master_create (nullptr, &h) == FC_ERR_NULL, "a null config");
        const fc_master_config c = goodConfig();
        ok (fc_master_create (&c, nullptr) == FC_ERR_NULL, "a null out-handle");

        fc_master_config bad = c; bad.header.abiVersion = 999u;
        h = 0xDEADBEEFu;
        ok (fc_master_create (&bad, &h) == FC_ERR_ABI_VERSION, "a version this build does not know");
        // NOT cleared. A caller reusing a variable that still holds a LIVE handle would otherwise have
        // it wiped by a call that failed on the version field, leaving that object unreachable and
        // undestroyable. A refused call is indistinguishable from one never made, arguments included.
        ok (h == 0xDEADBEEFu, "and the out-handle is left exactly as it was");

        bad = c; bad.header.structSize = (std::uint32_t) sizeof (bad) + 8u;
        ok (fc_master_create (&bad, &h) == FC_ERR_STRUCT_SIZE, "a struct size that disagrees with this build");
        bad = c; bad.header.structSize = (std::uint32_t) sizeof (bad) - 4u;
        ok (fc_master_create (&bad, &h) == FC_ERR_STRUCT_SIZE, "in either direction — a SHORTER one too");

        bad = c; bad.sampleRate = std::numeric_limits<double>::quiet_NaN();
        ok (fc_master_create (&bad, &h) == FC_ERR_NON_FINITE, "a NaN sample rate");
        bad = c; bad.sampleRate = std::numeric_limits<double>::infinity();
        ok (fc_master_create (&bad, &h) == FC_ERR_NON_FINITE, "an infinite one");
        bad = c; bad.compressorLookaheadMs = std::numeric_limits<double>::quiet_NaN();
        ok (fc_master_create (&bad, &h) == FC_ERR_NON_FINITE, "a NaN lookahead");

        bad = c; bad.channels = 0;
        ok (fc_master_create (&bad, &h) == FC_ERR_REFUSED_BY_CORE, "zero channels");
        bad = c; bad.channels = 99;
        ok (fc_master_create (&bad, &h) == FC_ERR_REFUSED_BY_CORE, "more channels than the core has");
        bad = c; bad.sampleRate = -48000.0;
        ok (fc_master_create (&bad, &h) == FC_ERR_REFUSED_BY_CORE, "a negative rate — the CORE's refusal, forwarded");
        bad = c; bad.internalBlock = 0;
        ok (fc_master_create (&bad, &h) == FC_ERR_REFUSED_BY_CORE, "an impossible internal quantum");
        bad = c; bad.monoBass = 1; bad.channels = 1;
        ok (fc_master_create (&bad, &h) == FC_ERR_REFUSED_BY_CORE,
            "mono-bass on a mono chain — a stage that would silently do nothing");
    }

    //==========================================================================
    // The handle is an index and a generation, not a pointer. Under emscripten, destroy-then-create
    // returned the SAME address 19 times out of 19 (emmalloc and dlmalloc both), so a pointer handle
    // would silently address the NEXT object and a linear memory has no unmapped page to trap on.
    group ("handles: a stale one is refused, never aliased");
    {
        fc_master a = make();
        ok (a != 0, "a live handle is non-zero");
        ok (fc_master_destroy (a) == FC_OK, "destroy accepts it");
        ok (fc_master_destroy (a) == FC_ERR_HANDLE, "and refuses the second destroy");

        std::int32_t lat = 0;
        ok (fc_master_latency (a, &lat) == FC_ERR_HANDLE, "every other entry point refuses it too");
        fc_master_stats st {}; FC_INIT (st);
        ok (fc_master_get_stats (a, &st) == FC_ERR_HANDLE, "including the getters");
        ok (fc_master_reset (a) == FC_ERR_HANDLE, "and reset");

        fc_master b = make();
        ok (b != 0 && b != a, "the next handle is a DIFFERENT value, whatever the allocator reused");
        ok (fc_master_latency (b, &lat) == FC_OK && lat > 0, "and the new one works");
        ok (fc_master_latency (a, &lat) == FC_ERR_HANDLE, "while the old one still does not");
        fc_master_destroy (b);

        ok (fc_master_latency (0, &lat) == FC_ERR_HANDLE, "handle 0 is never valid — a zeroed JS variable");
        ok (fc_master_reset (0xFFFFFFFFu) == FC_ERR_HANDLE, "nor is a made-up one");
    }

    group ("handles: a master is not a solution and vice versa");
    {
        fc_master h = make();
        fc_solution_summary sum {}; FC_INIT (sum);
        ok (fc_solution_summary_get ((fc_solution) h, &sum) == FC_ERR_HANDLE,
            "a master handle passed to a solution getter is refused, not reinterpreted");
        // The reverse, with a REAL solution handle. The earlier form of this check passed handle 0, which
        // is refused before any kind test runs — a check that looked like it covered the confusion and
        // could not have caught it.
        fc_master_params p = goodParams();
        fc_loudness_request req {}; fc_loudness_request_default (&req);
        const std::size_t frames = (std::size_t) (kFs * 4.0);
        auto in = tone (frames, kNch);
        std::vector<float> out (in.size(), 0.0f);
        fc_solution sol = 0;
        ok (fc_master_solve (h, &p, &req, in.data(), out.data(), (std::uint32_t) frames, &sol) == FC_OK
            && sol != 0, "PRECONDITION: a real solution handle exists");
        std::int32_t lat = 0;
        ok (fc_master_latency ((fc_master) sol, &lat) == FC_ERR_HANDLE,
            "a SOLUTION handle at a master entry point is refused, not reinterpreted");
        ok (fc_master_destroy ((fc_master) sol) == FC_ERR_HANDLE, "and cannot be destroyed as a master");
        ok (fc_solution_destroy (sol) == FC_OK, "while its own destroy accepts it");
        fc_master_destroy (h);
    }

    group ("the handle table is finite and says so");
    {
        std::vector<fc_master> live;
        fc_master h = 0;
        const fc_master_config c = goodConfig();
        while (fc_master_create (&c, &h) == FC_OK) live.push_back (h);
        ok (! live.empty(), "some handles were issued");
        ok (fc_master_create (&c, &h) == FC_ERR_EXHAUSTED, "and the table refuses rather than overwriting");
        for (fc_master x : live) fc_master_destroy (x);
        ok (fc_master_create (&c, &h) == FC_OK, "and recovers once they are destroyed");
        fc_master_destroy (h);
    }

    //==========================================================================
    group ("configure: refusals, and a refused call moves nothing");
    {
        fc_master h = make();
        fc_master_params p = goodParams();
        fc_master_resolved r {}; FC_INIT (r);
        ok (fc_master_configure (h, &p, &r) == FC_OK, "a good set is accepted");
        const double ceilingBefore = r.limiterCeilingDbTp;

        ok (fc_master_configure (h, nullptr, &r) == FC_ERR_NULL, "a null parameter block");
        ok (fc_master_configure (h, &p, nullptr) == FC_ERR_NULL, "a null resolved-out block");

        fc_master_resolved rBad {};                       // NOT stamped
        ok (fc_master_configure (h, &p, &rBad) == FC_ERR_ABI_VERSION,
            "an OUT struct with no header — the caller says which layout it has room for");

        fc_master_params bad = p; bad.header.structSize = 4u;
        ok (fc_master_configure (h, &bad, &r) == FC_ERR_STRUCT_SIZE, "a mis-sized parameter block");

        bad = p; bad.eqBands[3].type = 4242;
        ok (fc_master_configure (h, &bad, &r) == FC_ERR_ENUM, "a filter-type code that names nothing");
        bad = p; bad.eqBands[3].type = -1;
        ok (fc_master_configure (h, &bad, &r) == FC_ERR_ENUM, "including a negative one");
        bad = p; bad.compressor.mode = 3;
        ok (fc_master_configure (h, &bad, &r) == FC_ERR_ENUM, "one past the last compressor mode");
        bad = p; bad.compressor.link = 2;
        ok (fc_master_configure (h, &bad, &r) == FC_ERR_ENUM, "one past the last link mode");
        bad = p; bad.clipper.shape = 4;
        ok (fc_master_configure (h, &bad, &r) == FC_ERR_ENUM, "one past the last clipper shape");
        bad = p; bad.dither.shaping = 3;
        ok (fc_master_configure (h, &bad, &r) == FC_ERR_ENUM, "one past the last noise shaping");

        bad = p; bad.inputGainDb = std::numeric_limits<double>::quiet_NaN();
        ok (fc_master_configure (h, &bad, &r) == FC_ERR_NON_FINITE,
            "a NaN gain — the ONE class the core has no verdict for (it maps NaN to 0 dB silently)");
        bad = p; bad.eqBands[7].lanes[2].freq = std::numeric_limits<double>::infinity();
        ok (fc_master_configure (h, &bad, &r) == FC_ERR_NON_FINITE,
            "and it reaches into every lane of every band, not just the first");
        bad = p; bad.limiter.ceilingDbTp = std::numeric_limits<double>::quiet_NaN();
        ok (fc_master_configure (h, &bad, &r) == FC_ERR_NON_FINITE, "a NaN ceiling");
        bad = p; bad.clipper.driveDb = std::numeric_limits<float>::quiet_NaN();
        ok (fc_master_configure (h, &bad, &r) == FC_ERR_NON_FINITE, "a NaN in a FLOAT field too");

        // A refused configure must leave the chain exactly as it was. Without this the ABI would have a
        // refusal that still changed the render, which is the disease law 11 exists to cure.
        fc_master_resolved after {}; FC_INIT (after);
        ok (fc_master_resolved_get (h, &after) == FC_OK, "the chain still answers");
        ok (after.limiterCeilingDbTp == ceilingBefore,
            "and its resolved values are the ones from the last ACCEPTED configure");
        fc_master_destroy (h);
    }

    //==========================================================================
    // Ranges are the core's business, not this file's. A facade with its own range table is a second
    // copy of every stage's limits, and the first stage to move one leaves the two disagreeing.
    group ("configure does NOT range-check: an out-of-range value is CLAMPED by the core and reported");
    {
        fc_master h = make();
        fc_master_params p = goodParams();
        p.limiter.ceilingDbTp = 1.0e308;                 // finite, absurd
        fc_master_resolved r {}; FC_INIT (r);
        ok (fc_master_configure (h, &p, &r) == FC_OK, "accepted, not refused");
        ok (r.limiterCeilingDbTp <= 60.0 && r.limiterCeilingDbTp > 0.0,
            "and comes back CLAMPED to the limiter's own ceiling range");
        approx (r.limiterCeilingDbTp, 60.0, 1e-9, "which is +60 dBTP");
        fc_master_destroy (h);
    }

    //==========================================================================
    // The resolved values are the whole reason `configure` exists in this shape. Measured before it was
    // written: reading `resolved()` straight after a deferred `setParams` reports the PREVIOUS set —
    // 5.0000 dB out on the ceiling, 149.968 ms on the release — and on a fresh chain it reports the
    // core's defaults rather than anything the caller asked for.
    group ("configure hands back what the chain will ACTUALLY run, not the previous set");
    {
        fc_master h = make();
        fc_master_params p = goodParams();
        p.limiter.ceilingDbTp = -6.0;
        p.limiter.releaseMs   = 300.0;
        p.monoBass.enabled = 1; p.monoBass.frequencyHz = 200.0f; p.monoBass.lowWidth = 0.5f;
        fc_master_resolved r {}; FC_INIT (r);
        ok (fc_master_configure (h, &p, &r) == FC_OK, "configured");
        approx (r.limiterCeilingDbTp, -6.0, 1e-9, "the ceiling is the one just asked for");
        // NOT -6.0-exactly-round on the release: the limiter's coefficient makes a round trip through a
        // float, so 300 ms comes back as 300.277. That is the TRUTH about what will run, and reporting
        // the request instead would be this file inventing a number.
        ok (std::fabs (r.limiterReleaseMs - 300.0) < 1.0 && r.limiterReleaseMs != 300.0,
            "the release is the EFFECTIVE one, not the requested one");
        ok (r.monoBass.frequencyHz == 200.0f && r.monoBass.lowWidth == 0.5f, "and the mono-bass fold");
        ok (r.latencySamples > 0 && r.internalBlock == 256, "the geometry comes with it");
        ok (r.tapOversampleFactor >= 1, "and the tap stride, which is not oversampleFactor");
        fc_master_destroy (h);
    }

    group ("configure is REFUSED once audio has been handed over");
    {
        fc_master h = make();
        fc_master_params p = goodParams();
        fc_master_resolved r {}; FC_INIT (r);
        ok (fc_master_configure (h, &p, &r) == FC_OK, "before any audio");

        auto buf = tone (512, kNch);
        ok (fc_master_process (h, buf.data(), buf.data(), 512) == FC_OK, "a block goes through");
        ok (fc_master_configure (h, &p, &r) == FC_ERR_STATE,
            "and now configure refuses — re-preparing would silently discard the stream");
        ok (fc_master_reset (h) == FC_OK, "reset is the way back to the head of a stream");
        ok (fc_master_configure (h, &p, &r) == FC_OK, "after which configure is legal again");
        fc_master_destroy (h);
    }

    //==========================================================================
    group ("process: spans, alignment and overlap");
    {
        fc_master h = make();
        fc_master_params p = goodParams();
        fc_master_resolved r {}; FC_INIT (r);
        (void) fc_master_configure (h, &p, &r);

        auto a = tone (256, kNch);
        auto b = tone (256, kNch);
        ok (fc_master_process (h, nullptr, b.data(), 256) == FC_ERR_NULL, "a null input");
        ok (fc_master_process (h, a.data(), nullptr, 256) == FC_ERR_NULL, "a null output");
        ok (fc_master_process (h, a.data(), b.data(), 0) == FC_OK,
            "zero frames is the one true no-op, not an error");

        // 4-byte misalignment reads garbage in a release wasm build and only traps under -sSAFE_HEAP.
        const char* raw = reinterpret_cast<const char*> (a.data());
        const float* skew = reinterpret_cast<const float*> (raw + 1);
        ok (fc_master_process (h, skew, b.data(), 16) == FC_ERR_ALIGNMENT, "a misaligned input pointer");
        ok (fc_master_process (h, a.data(), const_cast<float*> (skew), 16) == FC_ERR_ALIGNMENT,
            "a misaligned output pointer");

        ok (fc_master_process (h, a.data(), a.data(), 256) == FC_OK, "in == out is legal and skips the copy");
        // A PARTIAL overlap is refused: memmove would keep the audio and eat part of a buffer the caller
        // declared const, and a contract that silently consumes its own input is worse than one that says no.
        ok (fc_master_process (h, a.data(), a.data() + 1, 128) == FC_ERR_SPAN, "a partial overlap forward");
        ok (fc_master_process (h, a.data() + 1, a.data(), 128) == FC_ERR_SPAN, "and backward");

        ok (fc_master_process (h, a.data(), b.data(), 0x80000000u) == FC_ERR_RANGE,
            "a frame count the core's own `int` cannot receive");
        fc_master_destroy (h);
    }

    //==========================================================================
    group ("process: RT-safety — the audio path allocates nothing");
    {
        fc_master h = make();
        fc_master_params p = goodParams();
        fc_master_resolved r {}; FC_INIT (r);
        (void) fc_master_configure (h, &p, &r);

        auto buf = tone (4096, kNch);
        (void) fc_master_process (h, buf.data(), buf.data(), 4096);   // warm every lazy path first
        const long long before = alloc::count.load();
        for (int i = 0; i < 8; ++i)
            (void) fc_master_process (h, buf.data(), buf.data(), 4096);
        std::int32_t lat = 0; (void) fc_master_latency (h, &lat);
        fc_master_stats st {}; FC_INIT (st);
        (void) fc_master_get_stats (h, &st);
        const long long after = alloc::count.load();
        okNoAlloc (after == before, "eight process() calls and two getters allocate nothing");
        ok (st.framesIn == 4096u * 9u, "and the frame counter agrees with what was handed over");
        fc_master_destroy (h);
    }

    //==========================================================================
    group ("process: a non-finite sample is COUNTED, not refused");
    {
        fc_master h = make();
        fc_master_params p = goodParams();
        fc_master_resolved r {}; FC_INIT (r);
        (void) fc_master_configure (h, &p, &r);

        auto buf = tone (1024, kNch);
        buf[100] = std::numeric_limits<float>::quiet_NaN();
        buf[200] = std::numeric_limits<float>::infinity();
        ok (fc_master_process (h, buf.data(), buf.data(), 1024) == FC_OK,
            "the call is ACCEPTED — refusing would throw away every good sample travelling with it");
        fc_master_stats st {}; FC_INIT (st);
        ok (fc_master_get_stats (h, &st) == FC_OK, "stats readable");
        ok (st.nonFiniteIn == 2u, "and both substitutions are counted, so the damage is visible");
        bool finite = true;
        for (float x : buf) if (! std::isfinite (x)) finite = false;
        ok (finite, "nothing non-finite reached the output");
        fc_master_destroy (h);
    }

    //==========================================================================
    group ("flush: capacity is binding, and a short one is refused rather than half-served");
    {
        fc_master h = make();
        fc_master_params p = goodParams();
        fc_master_resolved r {}; FC_INIT (r);
        (void) fc_master_configure (h, &p, &r);
        std::int32_t lat = 0;
        ok (fc_master_latency (h, &lat) == FC_OK && lat > 0, "the latency is answerable");

        std::vector<float> out ((std::size_t) lat * kNch, 0.0f);
        std::uint32_t got = 12345;
        ok (fc_master_flush (h, out.data(), 0, &got) == FC_ERR_CAPACITY, "a zero capacity");
        // NOT cleared: `written` is written only once every refusal is behind us. Clearing on entry
        // meant a call refused for ALIASING had already put a zero into the caller's audio — a refusal
        // that moved something, which is the one thing a refusal may not do.
        ok (got == 12345, "and a refused call leaves `written` exactly as it was");
        ok (fc_master_flush (h, out.data(), (std::uint32_t) lat - 1, &got) == FC_ERR_CAPACITY,
            "one frame short of the latency — the core keeps no arrears, so a partial drain is a trap");
        ok (fc_master_flush (h, nullptr, (std::uint32_t) lat, &got) == FC_ERR_NULL, "a null buffer");
        ok (fc_master_flush (h, out.data(), (std::uint32_t) lat, nullptr) == FC_ERR_NULL, "a null count");
        ok (fc_master_flush (h, out.data(), (std::uint32_t) lat, &got) == FC_OK, "the exact latency works");
        ok (got == (std::uint32_t) lat, "and drains all of it");
        fc_master_destroy (h);
    }

    //==========================================================================
    group ("measure_lra: a measurement DOES refuse a poisoned programme");
    {
        fc_master h = make();
        const std::size_t frames = (std::size_t) (kFs * 30.0);
        std::vector<float> in (frames * kNch, 0.0f);
        for (int c = 0; c < kNch; ++c)
            for (std::size_t i = 0; i < frames; ++i)
            {
                const int sec = (int) (i / 48000);
                const double amp = (sec % 6 < 3) ? 0.5 : 0.03;
                in[(std::size_t) c * frames + i] =
                    (float) (amp * std::sin (2.0 * 3.14159265358979 * (220.0 + 57.0 * c) * (double) i / kFs));
            }
        double lra = 0.0;
        ok (fc_master_measure_lra (h, in.data(), (std::uint32_t) frames, &lra) == FC_OK, "a clean programme");
        ok (lra > 1.0, "PRECONDITION: the fixture has a real loudness range to lose");
        const double clean = lra;

        for (int sec = 0; sec < 30; ++sec)
            if (sec % 6 < 3)
                for (std::size_t i = (std::size_t) sec * 48000; i < (std::size_t) (sec + 1) * 48000; ++i)
                    in[i] = std::numeric_limits<float>::quiet_NaN();
        double poisoned = -1.0;
        ok (fc_master_measure_lra (h, in.data(), (std::uint32_t) frames, &poisoned) == FC_ERR_REFUSED_BY_CORE,
            "and a poisoned one is REFUSED — unlike process(), a measurement has no gate to hide behind");
        ok (poisoned == -1.0, "the out-parameter is untouched by the refusal");

        // The number the refusal is protecting the caller from. Without the core's fix this read 21.4 LU
        // against a true 4.8 and came back as a success.
        ok (clean > 1.0 && clean < 12.0, "PRECONDITION: the clean answer is a plausible range");
        ok (fc_master_measure_lra (h, in.data(), 0, &lra) == FC_ERR_REFUSED_BY_CORE, "zero frames");
        ok (fc_master_measure_lra (h, nullptr, 1024, &lra) == FC_ERR_NULL, "a null programme");
        ok (fc_master_measure_lra (h, in.data(), 1024, nullptr) == FC_ERR_NULL, "a null out");
        fc_master_destroy (h);
    }

    //==========================================================================
    group ("solve: the verdict is DATA, the status is the ABI's");
    {
        fc_master h = make();
        fc_master_params p = goodParams();
        fc_loudness_request req {}; fc_loudness_request_default (&req);
        const std::size_t frames = (std::size_t) (kFs * 6.0);
        auto in = tone (frames, kNch);
        std::vector<float> out (in.size(), 0.0f);
        fc_solution sol = 0;

        ok (fc_master_solve (h, &p, &req, in.data(), out.data(), (std::uint32_t) frames, &sol) == FC_OK,
            "a request with NO target still returns a VERDICT rather than an ABI error");
        fc_solution_summary sum {}; FC_INIT (sum);
        ok (fc_solution_summary_get (sol, &sum) == FC_OK, "the summary is readable");
        ok (sum.status == FC_SOLVE_INVALID_REQUEST,
            "and it says InvalidRequest — which is the SOLVER's answer, not a marshalling failure");
        ok (fc_solution_destroy (sol) == FC_OK, "the solution is destroyed separately from the chain");
        ok (fc_solution_destroy (sol) == FC_ERR_HANDLE, "and its handle then goes stale like any other");

        req.targetLufs = -16.0;
        req.maxTruePeakDbTp = -1.0;
        req.maxPasses = 2;
        sol = 0;
        ok (fc_master_solve (h, &p, &req, in.data(), out.data(), (std::uint32_t) frames, &sol) == FC_OK,
            "a real request");
        FC_INIT (sum);
        ok (fc_solution_summary_get (sol, &sum) == FC_OK && sum.passes > 0, "which spent passes");
        fc_measurement meas {}; FC_INIT (meas);
        ok (fc_solution_measurement (sol, &meas) == FC_OK, "and carries a measurement of the DELIVERED render");
        ok (meas.gatingBlocks > 0, "with gating blocks in it");

        std::vector<fc_solve_pass> log (sum.logCount + 4);
        std::uint32_t written = 999;
        ok (fc_solution_log (sol, log.data(), (std::uint32_t) log.size(), &written) == FC_OK, "the log copies out");
        ok (written == (std::uint32_t) sum.logCount, "exactly logCount records");
        std::vector<fc_solve_pass> small (1);
        ok (fc_solution_log (sol, small.data(), 1, &written) == FC_OK && written == 1,
            "a short buffer is FILLED to its capacity, not overrun");
        ok (fc_solution_log (sol, nullptr, 4, &written) == FC_ERR_NULL, "a null buffer with a non-zero capacity");
        ok (fc_solution_log (sol, log.data(), 0, &written) == FC_OK && written == 0,
            "asking for nothing is not an error");

        // P65 — `written` MAY NOT POINT INTO THE RECORDS. It used to be cleared on entry and set after the copy, so a
        // `written` inside the `cap` records took the count over one of them on FC_OK, and a zero on a refusal. At the
        // boundary, both sides: the first and the last four bytes the span covers are refused, the four right after it
        // and right before it are not.
        {
            const std::uint32_t cap = (std::uint32_t) sum.logCount;
            ok (cap >= 2, "PRECONDITION: a log of two records or more, so the first and the last are different ones");
            std::vector<fc_solve_pass> buf ((std::size_t) cap + 2);
            const std::size_t total = buf.size() * sizeof (fc_solve_pass), span = (std::size_t) cap * sizeof (fc_solve_pass);
            std::memset (buf.data(), 0xA5, total);
            const std::vector<fc_solve_pass> pristine = buf;
            const auto unmoved = [&] { return std::memcmp (buf.data(), pristine.data(), total) == 0; };
            fc_solve_pass* const rec = buf.data() + 1;                 // one record of room on each side
            const auto at = [&] (std::ptrdiff_t byte)
            { return reinterpret_cast<std::uint32_t*> (reinterpret_cast<unsigned char*> (rec) + byte); };

            ok (fc_solution_log (sol, rec, cap, at (0)) == FC_ERR_SPAN && unmoved(),
                "`written` on the first record's first four bytes is refused, and nothing in the buffer moves");
            ok (fc_solution_log (sol, rec, cap, at ((std::ptrdiff_t) span - 4)) == FC_ERR_SPAN && unmoved(),
                "and on the last record's last four bytes");
            // THE WHOLE CAPACITY, not the records the log fills: one record more than logCount, `written` on its
            // `violated`. A check over the `n` records actually copied passed every line above and wrote the count here
            // (the code-review round, by mutation) — a refusal that would depend on how many renders the search took.
            ok (fc_solution_log (sol, rec, cap + 1, at ((std::ptrdiff_t) span + (std::ptrdiff_t) offsetof (fc_solve_pass, violated)))
                    == FC_ERR_SPAN && unmoved(),
                "and inside capacity the log does not fill");

            // The same fill, because the copy is field by field: a record's four bytes of tail padding are not written.
            std::vector<fc_solve_pass> want (cap);
            std::memset (want.data(), 0xA5, span);
            std::uint32_t wn = 0;
            ok (fc_solution_log (sol, want.data(), cap, &wn) == FC_OK && wn == cap, "the log, read into a buffer of its own");
            std::memcpy (buf.data(), pristine.data(), total);   // each check on a clean buffer, whatever the last one did
            ok (fc_solution_log (sol, rec, cap, at ((std::ptrdiff_t) span)) == FC_OK && *at ((std::ptrdiff_t) span) == cap
                && std::memcmp (rec, want.data(), span) == 0,
                "`written` right AFTER the records is legal: FC_OK, the count, and every record whole");
            std::memcpy (buf.data(), pristine.data(), total);
            ok (fc_solution_log (sol, rec, cap, at (-4)) == FC_OK && *at (-4) == cap
                && std::memcmp (rec, want.data(), span) == 0,
                "and right BEFORE them");

            // No span is named by cap == 0, so there is nothing for `written` to alias.
            std::memcpy (buf.data(), pristine.data(), total);
            ok (fc_solution_log (sol, rec, 0, at (0)) == FC_OK && *at (0) == 0u,
                "cap == 0 names no span: `written` anywhere is FC_OK and a zero");

            // A refused call leaves `written` as it was — it is cleared only once every check is behind the call.
            std::uint32_t kept = 777u;
            ok (fc_solution_log (sol, nullptr, 4, &kept) == FC_ERR_NULL && kept == 777u,
                "a refusal for a null buffer leaves `written` untouched");
            auto* const misaligned = reinterpret_cast<fc_solve_pass*> (reinterpret_cast<unsigned char*> (rec) + 4);
            ok (fc_solution_log (sol, misaligned, 1, &kept) == FC_ERR_ALIGNMENT && kept == 777u,
                "and so does one for alignment");
        }
        fc_solution_destroy (sol);

        // A search reads its input again on every pass, so rendering in place would make every pass after
        // the first read the previous pass's master.
        ok (fc_master_solve (h, &p, &req, in.data(), in.data(), (std::uint32_t) frames, &sol) == FC_ERR_SPAN,
            "solving IN PLACE is refused");
        ok (fc_master_solve (h, &p, &req, in.data(), in.data() + 4, (std::uint32_t) 64, &sol) == FC_ERR_SPAN,
            "and so is a partial overlap");

        // P64 — THE CORE REFUSES WHAT THE FACADE REFUSES. Output channel 0 laid over input channel 1 is, planar, `out =
        // in + frames`: the facade refuses it on its span check. The core used to refuse only a channel against itself,
        // and answered this call with a master read from its own previous pass — a direct C++ caller had a road the
        // facade exists to close. The same buffer through both, then the layout one plane further on through both.
        {
            using namespace felitronics::mastering;
            const std::size_t F = (std::size_t) kFs;
            const auto t = tone (F, kNch);
            std::vector<float> pool (4 * F, 0.25f);
            std::copy (t.begin(), t.end(), pool.begin());
            const std::vector<float> pristine = pool;
            const auto unmoved = [&] { return std::memcmp (pool.data(), pristine.data(), pool.size() * sizeof (float)) == 0; };

            fc_solution held = 0xABCDu;
            ok (fc_master_solve (h, &p, &req, pool.data(), pool.data() + F, (std::uint32_t) F, &held) == FC_ERR_SPAN
                && held == 0xABCDu && unmoved(),
                "the facade refuses output channel 0 over input channel 1, writing nothing");

            MasteringChain chain;
            OfflineRenderer rend;
            TargetLoudnessSolver solver;
            const bool prepared = rend.prepare (kNch, 4096) && chain.prepare (kFs, kNch, MasteringChainConfig {})
                               && solver.prepare (kFs, kNch, rend.blockSize(), chain.internalBlock(),
                                                  chain.tapOversampleFactor());
            ok (prepared, "the direct C++ search is prepared");
            if (prepared)
            {
                LoudnessRequest lr {};
                lr.targetLufs = req.targetLufs; lr.maxTruePeakDbTp = req.maxTruePeakDbTp; lr.maxPasses = req.maxPasses;
                const float* ip[2] { pool.data(), pool.data() + F };
                float*       op[2] { pool.data() + F, pool.data() + 2 * F };
                const LoudnessSolution direct = solver.solve (chain, rend, MasteringChainParams {}, ip, op, kNch, (int) F, lr);
                ok (direct.status == MasteringSolveStatus::InvalidRequest && direct.passes == 0 && unmoved(),
                    "and the core refuses the SAME call, before a render, writing nothing");

                // One plane further on the two spans are edge to edge: legal through both, and neither calls it a refusal.
                fc_solution sol2 = 0;
                ok (fc_master_solve (h, &p, &req, pool.data(), pool.data() + 2 * F, (std::uint32_t) F, &sol2) == FC_OK,
                    "the facade accepts output right after input");
                fc_solution_summary s2 {}; FC_INIT (s2);
                ok (fc_solution_summary_get (sol2, &s2) == FC_OK && s2.status != FC_SOLVE_INVALID_REQUEST,
                    "with a verdict that is not a refusal");
                fc_solution_destroy (sol2);
                std::memcpy (pool.data(), pristine.data(), pool.size() * sizeof (float));
                const float* ip2[2] { pool.data(), pool.data() + F };
                float*       op2[2] { pool.data() + 2 * F, pool.data() + 3 * F };
                const LoudnessSolution direct2 = solver.solve (chain, rend, MasteringChainParams {}, ip2, op2, kNch, (int) F, lr);
                ok (direct2.status != MasteringSolveStatus::InvalidRequest && direct2.passes > 0,
                    "and so does the core");
            }
        }
        fc_master_destroy (h);
    }


    //==========================================================================
    // Every check below was written because a MUTATION SURVIVED without it, or because a crew round
    // named the hole. Recorded that way rather than tidied in: a test whose reason is invisible is the
    // first one somebody deletes.
    group ("holes the mutation stand and the review round found");
    {
        // 1. `reset()` must clear the non-finite counter, or it describes two streams at once.
        {
            fc_master h = make();
            fc_master_params p = goodParams();
            fc_master_resolved r {}; FC_INIT (r);
            (void) fc_master_configure (h, &p, &r);
            auto buf = tone (1024, kNch);
            buf[10] = std::numeric_limits<float>::quiet_NaN();
            (void) fc_master_process (h, buf.data(), buf.data(), 1024);
            fc_master_stats st {}; FC_INIT (st);
            (void) fc_master_get_stats (h, &st);
            ok (st.nonFiniteIn == 1u, "PRECONDITION: the counter moved, so the check below is live");
            ok (fc_master_reset (h) == FC_OK, "reset");
            FC_INIT (st);
            (void) fc_master_get_stats (h, &st);
            ok (st.nonFiniteIn == 0u, "and the counter describes THIS stream, not the previous one");
            ok (st.framesIn == 0u && st.framesFlushed == 0u, "as do the frame counters");
            fc_master_destroy (h);
        }

        // 2. The frames x channels product must be refused where it leaves 32 bits — and the test has to
        //    reach THAT check rather than the INT_MAX one, which fires first for a larger count.
        //    0x40000000 frames is under INT_MAX; times 2 channels times 4 bytes it is 2^33.
        {
            fc_master h = make();
            fc_master_params p = goodParams();
            fc_master_resolved r {}; FC_INIT (r);
            (void) fc_master_configure (h, &p, &r);
            auto buf = tone (64, kNch);
            ok (fc_master_process (h, buf.data(), buf.data(), 0x40000000u) == FC_ERR_SPAN,
                "a byte span past 32 bits is SPAN, not RANGE — the two checks are different and both live");
            ok (fc_master_process (h, buf.data(), buf.data(), 0x80000000u) == FC_ERR_RANGE,
                "while a frame count past INT_MAX is RANGE");
            fc_master_destroy (h);
        }

        // 3. A handle is only the value `packHandle` produced. Without the high-bit test each live
        //    handle has 65 536 accepted spellings and a made-up one can destroy somebody's object.
        {
            fc_master h = make();
            std::int32_t lat = 0;
            ok (fc_master_latency (h, &lat) == FC_OK, "the real handle works");
            ok (fc_master_latency (h | 0x12340000u, &lat) == FC_ERR_HANDLE,
                "a value with bits this ABI never sets is refused, not accepted as the same handle");
            ok (fc_master_destroy (h | 0xFFFF0000u) == FC_ERR_HANDLE, "and cannot destroy it either");
            ok (fc_master_latency (h, &lat) == FC_OK, "the object is still there");
            fc_master_destroy (h);
        }

        // 4. A refused create must not spend a generation: the slot never issued a handle, so there is
        //    nothing for a stale one to alias. Without this, refused creates retire the whole table.
        {
            fc_master_config bad = goodConfig();
            bad.sampleRate = -1.0;
            fc_master h = 0;
            int refused = 0;
            for (int i = 0; i < 600; ++i)
                if (fc_master_create (&bad, &h) == FC_ERR_REFUSED_BY_CORE) ++refused;
            ok (refused == 600, "600 creates, all refused by the core");
            const fc_master_config good = goodConfig();
            fc_master live = 0;
            ok (fc_master_create (&good, &live) == FC_OK,
                "600 refused creates later the table still issues handles");
            fc_master_destroy (live);
        }

        // 5. A solve the SOLVER refused before rendering must leave the stream alone. Measured before
        //    the fix: process(64), then a targetless solve, and the configure that had correctly
        //    answered FC_ERR_STATE was accepted and destroyed the 64 frames sitting in the FIFO.
        {
            fc_master h = make();
            fc_master_params p = goodParams();
            fc_master_resolved r {}; FC_INIT (r);
            (void) fc_master_configure (h, &p, &r);
            auto buf = tone (64, kNch);
            ok (fc_master_process (h, buf.data(), buf.data(), 64) == FC_OK, "a partial quantum is in flight");
            ok (fc_master_configure (h, &p, &r) == FC_ERR_STATE, "PRECONDITION: the guard is armed");

            fc_loudness_request req {}; fc_loudness_request_default (&req);   // no target: InvalidRequest
            const std::size_t frames = (std::size_t) (kFs * 4.0);
            auto in = tone (frames, kNch);
            std::vector<float> out (in.size(), 0.0f);
            fc_solution sol = 0xABCD;
            // A SEARCH IS NOT EXEMPT FROM THE STREAM GUARD. It used to reset the chain out from under a
            // stream in progress and answer FC_OK, so a solve destroyed exactly what `configure` had
            // just refused to destroy. One policy for both, and `reset` is the way out.
            ok (fc_master_solve (h, &p, &req, in.data(), out.data(), (std::uint32_t) frames, &sol)
                    == FC_ERR_STATE, "solve is refused while a stream is in progress");
            ok (sol == 0xABCD, "and the out-handle is untouched, so a live solution cannot be wiped by it");
            ok (fc_master_configure (h, &p, &r) == FC_ERR_STATE, "the guard is still armed");
            ok (fc_master_reset (h) == FC_OK, "reset ends the stream");
            ok (fc_master_solve (h, &p, &req, in.data(), out.data(), (std::uint32_t) frames, &sol) == FC_OK,
                "and now the solve returns a verdict");
            fc_solution_summary sum {}; FC_INIT (sum);
            (void) fc_solution_summary_get (sol, &sum);
            ok (sum.status == FC_SOLVE_INVALID_REQUEST, "which is the solver's own InvalidRequest");
            // PROCESS FIRST, CONFIGURE AFTER, and the order is the whole check. `configure` CLEARS the
            // solved mark, so asking it first and `process` second tested nothing at all — a mutant
            // that set the mark unconditionally sailed through, because by the time `process` ran the
            // mark had already been wiped by the check in front of it. A dead check that looks alive.
            auto blk2 = tone (256, kNch);
            ok (fc_master_process (h, blk2.data(), blk2.data(), 256) == FC_OK,
                "a solve that never rendered did not mark the handle as solved");
            ok (fc_master_reset (h) == FC_OK, "reset");
            ok (fc_master_configure (h, &p, &r) == FC_OK,
                "and the chain is configurable, as a chain that never rendered should be");
            fc_solution_destroy (sol);
            fc_master_destroy (h);
        }

        // 6. A solve that DID render has reset the chain, so the guard must drop with it. Both halves
        //    are asserted, because only the pair says the condition is on the verdict and not on luck.
        {
            fc_master h = make();
            fc_master_params p = goodParams();
            fc_master_resolved r {}; FC_INIT (r);
            (void) fc_master_configure (h, &p, &r);
            fc_loudness_request req {}; fc_loudness_request_default (&req);
            req.targetLufs = -16.0; req.maxTruePeakDbTp = -1.0; req.maxPasses = 1;
            const std::size_t frames = (std::size_t) (kFs * 4.0);
            auto in = tone (frames, kNch);
            std::vector<float> out (in.size(), 0.0f);
            fc_solution sol = 0;
            ok (fc_master_solve (h, &p, &req, in.data(), out.data(), (std::uint32_t) frames, &sol) == FC_OK,
                "a real solve");
            ok (fc_master_configure (h, &p, &r) == FC_OK,
                "and the guard is DOWN, because the search reset the chain underneath it");
            fc_solution_destroy (sol);
            fc_master_destroy (h);
        }

        // 7. BS.1770 weights, without which a surround search cannot be expressed through this ABI at
        //    all — the standard weights Ls/Rs at 1.41 and excludes LFE.
        {
            fc_master h = make();
            ok (fc_master_set_channel_weight (h, 0, 1.0) == FC_OK, "a legal weight");
            ok (fc_master_set_channel_weight (h, 1, 1.41) == FC_OK, "the surround weight");
            ok (fc_master_set_channel_weight (h, 2, 0.0) == FC_OK, "and zero, which is what LFE gets");
            ok (fc_master_set_channel_weight (h, -1, 1.0) == FC_ERR_RANGE, "a negative channel");
            ok (fc_master_set_channel_weight (h, 999, 1.0) == FC_ERR_RANGE, "one past the last channel");
            ok (fc_master_set_channel_weight (h, 0, -1.0) == FC_ERR_RANGE,
            "a negative weight is out of RANGE — calling a number you can see 'non-finite' is a lie");
            ok (fc_master_set_channel_weight (h, 0, std::numeric_limits<double>::quiet_NaN())
                    == FC_ERR_NON_FINITE, "a NaN weight");
            ok (fc_master_set_channel_weight (0, 0, 1.0) == FC_ERR_HANDLE, "and it checks the handle first");
            fc_master_destroy (h);
        }

        // 8. Null scalar out-parameters, at every entry point that has one. The heap-bounds half of the
        //    same guard cannot be exercised natively — `inHeap` has nothing to bound outside a linear
        //    memory — and saying so is better than a test that passes for that reason.
        {
            fc_master h = make();
            std::int32_t lat = 0;
            ok (fc_master_latency (h, nullptr) == FC_ERR_NULL, "latency");
            double lra = 0.0;
            ok (fc_master_measure_lra (h, nullptr, 16, &lra) == FC_ERR_NULL, "measure_lra input");
            fc_master_config c = goodConfig();
            ok (fc_master_create (&c, nullptr) == FC_ERR_NULL, "create out-handle");
            fc_master_params p = goodParams();
            fc_loudness_request req {}; fc_loudness_request_default (&req);
            auto in = tone (256, kNch);
            std::vector<float> out (in.size(), 0.0f);
            ok (fc_master_solve (h, &p, &req, in.data(), out.data(), 256, nullptr) == FC_ERR_NULL,
                "solve out-solution");
            ok (fc_master_latency (h, &lat) == FC_OK, "and the handle survived all of it");
            fc_master_destroy (h);
        }
    }


    //==========================================================================
    // The second crew round measured that none of this was covered: `solve`'s DELIVERED audio was never
    // compared to anything, seven of the thirteen resolved fields were never read, and a handle that had
    // solved could go straight back to `process` and render a third thing nobody asked for.
    group ("solve: the delivered audio is the CORE's, bit for bit");
    {
        using namespace felitronics::mastering;
        fc_master h = make();
        fc_master_params p = goodParams();
        p.limiter.ceilingDbTp = -1.3;
        p.compressor.thresholdDb = -17.3;
        fc_loudness_request req {}; fc_loudness_request_default (&req);
        req.targetLufs = -16.0; req.maxTruePeakDbTp = -1.0; req.maxPasses = 2;
        // NOT the default 0.1: a field left at its default cannot catch a mapping that drops it, and a
        // mutant that dropped `toleranceLu` survived for exactly that reason.
        req.toleranceLu = 0.002;
        req.truePeakAimDb = 0.073;
        // EVERY request field off its default. A crew round dropped ten of them one at a time and every
        // one survived, because every fixture set target/tp/maxPasses and left the rest where the
        // default writer had put them — the same "a field at its default cannot catch a missing
        // mapping" that the parameter fixture had already been repaired for.
        req.minPlrDb = 3.7;
        req.maxLraLossLu = 7.3;
        req.inputLoudnessRangeLu = 5.3;
        req.activityThresholdDb = 0.37;
        req.initialGainDb = 1.7;
        req.limiterGr.limitDb = 17.3;    req.limiterGr.statistic = FC_GR_P95;
        req.compressorGr.limitDb = 23.7; req.compressorGr.statistic = FC_GR_MEAN;
        // A TOLERANCE FINER THAN THE SEARCH'S OWN STEP, so dropping it (back to the default 0.1) really
        // does change where the search stops. At 0.037 a mutant that dropped it still converged inside
        // both tolerances and survived — a field is only pinned where its value is the binding one.

        const std::size_t frames = (std::size_t) (kFs * 6.0);
        auto in = tone (frames, kNch);
        std::vector<float> viaAbi (in.size(), 0.0f);
        fc_solution sol = 0;
        ok (fc_master_solve (h, &p, &req, in.data(), viaAbi.data(), (std::uint32_t) frames, &sol) == FC_OK,
            "the ABI solve ran");
        fc_solution_summary sum {}; FC_INIT (sum);
        (void) fc_solution_summary_get (sol, &sum);

        // The same search through the C++ API, with the same lifecycle order `fc_master_configure` takes.
        MasteringChainConfig cc {};
        cc.internalBlock = 256; cc.monoBass = true; cc.clipper = true;
        MasteringChainParams cp {};
        cp.limiter.ceilingDbTp = -1.3; cp.compressor.thresholdDb = -17.3;
        MasteringChain chain;
        OfflineRenderer rend;
        TargetLoudnessSolver solver;
        std::vector<float> viaCpp (in.size(), 0.0f);
        const bool prepared = rend.prepare (kNch, 4096)
                           && chain.prepare (kFs, kNch, cc)
                           && solver.prepare (kFs, kNch, rend.blockSize(), chain.internalBlock(),
                                              chain.tapOversampleFactor());
        ok (prepared, "the direct C++ search is prepared the same way");
        if (prepared)
        {
            LoudnessRequest lr {};
            lr.targetLufs = -16.0; lr.maxTruePeakDbTp = -1.0; lr.maxPasses = 2;
            lr.toleranceLu = 0.002; lr.truePeakAimDb = 0.073;
            lr.minPlrDb = 3.7; lr.maxLraLossLu = 7.3; lr.inputLoudnessRangeLu = 5.3;
            lr.activityThresholdDb = 0.37; lr.initialGainDb = 1.7;
            lr.limiterGr.limitDb = 17.3;    lr.limiterGr.statistic = GrStatistic::P95;
            lr.compressorGr.limitDb = 23.7; lr.compressorGr.statistic = GrStatistic::Mean;
            const float* ip[2] { in.data(), in.data() + frames };
            float*       op[2] { viaCpp.data(), viaCpp.data() + frames };
            const LoudnessSolution direct = solver.solve (chain, rend, cp, ip, op, kNch, (int) frames, lr);
            ok ((int) direct.status == sum.status, "the two verdicts agree");
            ok (direct.preLimiterGainDb == sum.preLimiterGainDb, "and the gain, bit for bit");
            ok (direct.ceilingDbTp == sum.ceilingDbTp, "and the ceiling");
            std::size_t diff = 0;
            for (std::size_t i = 0; i < viaAbi.size(); ++i)
                if (std::memcmp (&viaAbi[i], &viaCpp[i], sizeof (float)) != 0) ++diff;
            ok (diff == 0, "and the DELIVERED audio is bit-identical through the ABI");
            double amp = 0.0, moved = 0.0;
            for (std::size_t i = 0; i < viaAbi.size(); ++i)
            {
                amp = std::max (amp, (double) std::fabs (viaAbi[i]));
                moved = std::max (moved, (double) std::fabs (viaAbi[i] - in[i]));
            }
            ok (amp > 0.05 && moved > 0.01, "PRECONDITION: the search actually rendered something");

            fc_measurement meas {}; FC_INIT (meas);
            (void) fc_solution_measurement (sol, &meas);
            ok (meas.integratedLufs == direct.measured.integratedLufs, "the measurement crosses unchanged");
            ok (meas.truePeakDbTp == direct.measured.truePeakDbTp, "true peak too");
            ok (meas.loudnessRangeLu == direct.measured.loudnessRangeLu, "and the loudness range");
            // `valid` is a FLAG, and a mutant that hard-coded it to 1 survived because nothing read it.
            ok (meas.compressor.valid == (direct.measured.compressor.valid ? 1 : 0)
                && meas.limiter.valid == (direct.measured.limiter.valid ? 1 : 0),
                "and both gain-reduction summaries carry the core's own validity, not a constant");
            ok (meas.limiter.maxDb == direct.measured.limiter.maxDb
                && meas.compressor.p95Db == direct.measured.compressor.p95Db,
                "and their numbers cross unchanged");
            // THE PASS BUDGET REACHED THE SOLVER. Dropping `maxPasses` from the mapping left the default
            // (4) in its place, and a search that converged inside two passes anyway did not notice.
            ok (sum.passes <= 2 + 1,
                "the search spent no more than the budget the request carried (+1 for delivering)");

            // EVERY field of the summary, the measurement, both gain-reduction summaries and the log —
            // not a representative sample. A crew round crossed, zeroed or constant-folded 32 of them
            // one at a time and every mutant survived, because the comparison read six numbers.
            ok (sum.binding == (std::int32_t) direct.binding, "the binding constraint");
            ok (sum.alsoViolated == direct.alsoViolated, "the violation mask");
            ok (sum.logCount == direct.logCount, "the log length");
            ok (sum.activityThresholdDb == direct.activityThresholdDb, "the echoed activity threshold");
            ok (sum.achievedBelowLufs == direct.achievedBelowLufs
                && sum.achievedAboveLufs == direct.achievedAboveLufs, "both bracketing loudnesses");
            ok (sum.gainBelowDb == direct.gainBelowDb && sum.gainAboveDb == direct.gainAboveDb,
                "and both bracketing gains, in the right order");

            ok (meas.samplePeakDb == direct.measured.samplePeakDb, "the sample peak");
            ok (meas.plrDb == direct.measured.plrDb, "the peak-to-loudness ratio");
            ok (meas.limiterMaxReconstructedPeakDb == direct.measured.limiterMaxReconstructedPeakDb,
                "the reconstructed peak the limiter saw");
            ok (meas.latencySamples == direct.measured.latencySamples, "the latency");
            ok (meas.gatingBlocks == direct.measured.gatingBlocks
                && meas.droppedBlocks == direct.measured.droppedBlocks
                && meas.nonFiniteSubHops == direct.measured.nonFiniteSubHops, "the three block counters");
            ok (meas.loudnessValid == (direct.measured.loudnessValid ? 1 : 0)
                && meas.lraValid == (direct.measured.lraValid ? 1 : 0), "and both validity flags");

            auto sameGr = [] (const fc_gr_stats& a, const GainReductionStats& b)
            {
                return a.meanDb == b.meanDb && a.p95Db == b.p95Db && a.maxDb == b.maxDb
                    && a.activeFraction == b.activeFraction && a.frames == b.frames
                    && a.nonFinite == b.nonFinite && a.aboveRange == b.aboveRange
                    && a.valid == (b.valid ? 1 : 0);
            };
            ok (sameGr (meas.compressor, direct.measured.compressor)
                && sameGr (meas.limiter, direct.measured.limiter),
                "every field of both gain-reduction summaries, and NOT crossed between the two");
            ok (meas.compressor.maxDb != meas.limiter.maxDb || direct.measured.compressor.maxDb
                                                            == direct.measured.limiter.maxDb,
                "PRECONDITION: the two summaries differ, so a swap between them would be visible");

            std::vector<fc_solve_pass> log ((std::size_t) sum.logCount + 1);
            std::uint32_t written = 0;
            (void) fc_solution_log (sol, log.data(), (std::uint32_t) log.size(), &written);
            ok (written == (std::uint32_t) direct.logCount, "the log has the core's length");
            bool logSame = (written > 0);
            for (std::uint32_t i = 0; i < written; ++i)
            {
                const SolvePassRecord& d = direct.log[i];
                if (log[i].gainDb != d.gainDb || log[i].ceilingDb != d.ceilingDb
                    || log[i].integratedLufs != d.integratedLufs || log[i].truePeakDbTp != d.truePeakDbTp
                    || log[i].plrDb != d.plrDb || log[i].limiterMaxGrDb != d.limiterMaxGrDb
                    || log[i].loudnessRangeLu != d.loudnessRangeLu
                    || log[i].violated != d.violated) logSame = false;
            }
            ok (logSame, "and EVERY field of every pass record crosses unchanged — all eight");

            // v4 — both traces cross unchanged, every field of every bucket, and are not crossed between the stages.
            fc_measurement m4 {}; FC_INIT (m4);
            (void) fc_solution_measurement (sol, &m4);
            ok (m4.compressorGrTraceBuckets == direct.compressorTrace.buckets && m4.limiterGrTraceBuckets == direct.limiterTrace.buckets
                && m4.compressorGrTraceValid == (direct.compressorTrace.valid ? 1 : 0)
                && m4.limiterGrTraceValid == (direct.limiterTrace.valid ? 1 : 0), "v4: the traces' counts and validity");
            auto sameTrace = [&] (std::int32_t stage, const GainReductionTrace& d)
            {
                std::vector<fc_gr_trace_bucket> tb ((std::size_t) std::max (d.buckets, 0) + 1u);   // one past the trace: `w` is not the capacity
                std::uint32_t w = 0;
                if (fc_solution_gr_trace (sol, stage, tb.data(), (std::uint32_t) tb.size(), &w) != FC_OK || (int) w != d.buckets || w == 0) return false;
                for (std::uint32_t i = 0; i < w; ++i)
                    if (tb[i].maxDb != d.bucket[i].maxDb || tb[i].meanDb != d.bucket[i].meanDb
                        || tb[i].samples != d.bucket[i].samples || tb[i].nonFinite != d.bucket[i].nonFinite) return false;
                return true;
            };
            ok (sameTrace (FC_GR_STAGE_COMPRESSOR, direct.compressorTrace) && sameTrace (FC_GR_STAGE_LIMITER, direct.limiterTrace),
                "v4: EVERY field of every bucket of both traces crosses unchanged, compressor as compressor, limiter as limiter");
            ok (direct.compressorTrace.bucket[0].samples != direct.limiterTrace.bucket[0].samples,
                "PRECONDITION: the two traces differ (the limiter counts oversampled sub-samples), so a swap would show");
            {
                double lm = 0.0, cm = 0.0;
                for (int i = 0; i < direct.limiterTrace.buckets; ++i) lm = std::max (lm, direct.limiterTrace.bucket[i].maxDb);
                for (int i = 0; i < direct.compressorTrace.buckets; ++i) cm = std::max (cm, direct.compressorTrace.bucket[i].maxDb);
                // The limiter is idle on this fixture (its constraints keep the gain down); the limiter's live case is
                // compared in the v4 group, on a request loud enough to engage it.
                ok (cm > 0.1, "PRECONDITION: the compressor trace is live (" + std::to_string (cm) + " dB; limiter "
                    + std::to_string (lm) + " dB) — a comparison of zeros would pass a getter that lost its values");
            }
        }
        fc_solution_destroy (sol);
        fc_master_destroy (h);
    }

    group ("solve leaves the handle holding the SOLVER's parameters, and says so");
    {
        fc_master h = make();
        fc_master_params p = goodParams();
        fc_master_resolved r {}; FC_INIT (r);
        (void) fc_master_configure (h, &p, &r);
        fc_loudness_request req {}; fc_loudness_request_default (&req);
        req.targetLufs = -16.0; req.maxTruePeakDbTp = -1.0; req.maxPasses = 1;
        const std::size_t frames = (std::size_t) (kFs * 4.0);
        auto in = tone (frames, kNch);
        std::vector<float> out (in.size(), 0.0f);
        fc_solution sol = 0;
        ok (fc_master_solve (h, &p, &req, in.data(), out.data(), (std::uint32_t) frames, &sol) == FC_OK,
            "a solve that ran");

        auto blk = tone (256, kNch);
        ok (fc_master_process (h, blk.data(), blk.data(), 256) == FC_ERR_STATE,
            "process is REFUSED — the chain holds the solver's gain and ceiling, not the caller's");
        std::vector<float> drain (2048 * kNch, 0.0f);
        std::uint32_t w = 0;
        ok (fc_master_flush (h, drain.data(), 2048, &w) == FC_ERR_STATE, "and so is flush");
        ok (fc_master_reset (h) == FC_OK, "reset is accepted");
        ok (fc_master_process (h, blk.data(), blk.data(), 256) == FC_ERR_STATE,
            "and does NOT lift it: reset clears audio state and keeps the solver's parameters");
        ok (fc_master_configure (h, &p, &r) == FC_OK, "configure is the way back");
        ok (fc_master_process (h, blk.data(), blk.data(), 256) == FC_OK, "and process works again");
        fc_solution_destroy (sol);
        fc_master_destroy (h);
    }

    group ("resolved: every field is read out of the chain, not invented");
    {
        // Seven of these were never asserted anywhere, so a crossed or zeroed field was invisible.
        // `tapOversampleFactor` is the sharp one: the header spends a paragraph on its NOT being
        // `oversampleFactor`, and with a clipper but no limiter the two really do differ.
        fc_master_config c = goodConfig();
        c.clipper = 1; c.limiter = 0; c.oversampleFactor = 4;
        fc_master h = 0;
        ok (fc_master_create (&c, &h) == FC_OK, "a chain with a clipper and NO limiter");
        fc_master_params p = goodParams();
        fc_master_resolved r {}; FC_INIT (r);
        ok (fc_master_configure (h, &p, &r) == FC_OK, "configured");
        ok (r.oversampleFactor == 4, "oversampleFactor reports the CLIPPER's factor when there is no limiter");
        ok (r.tapOversampleFactor == 1, "while the TAP stride is 1 — the two are different questions");
        ok (r.limiterLatency == 0 && r.limiterLookahead == 0, "and the absent limiter costs no latency");
        ok (r.clipperLatency > 0, "while the clipper's is real");
        ok (r.latencySamples == r.internalBlock + r.compressorLookahead + r.clipperLatency + r.limiterLatency,
            "and the total is the sum of the present stages plus the quantum");
        ok (r.compressorTapOffset == 0, "the compressor tap sits at the chain input");
        ok (r.limiterTapOffset == r.compressorLookahead + r.clipperLatency,
            "and the limiter tap behind the stages in front of it");
        fc_master_destroy (h);

        fc_master_config c2 = goodConfig();
        c2.limiter = 1; c2.clipper = 1;
        fc_master h2 = 0;
        ok (fc_master_create (&c2, &h2) == FC_OK, "and with a limiter present");
        fc_master_resolved r2 {}; FC_INIT (r2);
        ok (fc_master_configure (h2, &p, &r2) == FC_OK, "configured");
        ok (r2.tapOversampleFactor == r2.oversampleFactor,
            "the tap stride becomes the limiter's factor — PRECONDITION that the check above is live");
        ok (r2.limiterLookahead > 0 && r2.limiterLatency > r2.limiterLookahead,
            "the limiter's lookahead and its total latency are different numbers, and neither is zero");
        fc_master_destroy (h2);
    }

    group ("a scalar out-parameter may not alias the audio it reports on");
    {
        fc_master h = make();
        fc_master_params p = goodParams();
        fc_master_resolved r {}; FC_INIT (r);
        (void) fc_master_configure (h, &p, &r);
        std::int32_t lat = 0; (void) fc_master_latency (h, &lat);
        std::vector<float> out ((std::size_t) lat * kNch + 8, 0.0f);
        // `written` pointed INTO the drain buffer: flush would write the audio and then overwrite its
        // first sample with the frame count. The in/out overlap rule does not see this class at all.
        auto* aliased = reinterpret_cast<std::uint32_t*> (out.data());
        ok (fc_master_flush (h, out.data(), (std::uint32_t) lat, aliased) == FC_ERR_SPAN,
            "flush refuses a `written` that points inside its own output");
        std::uint32_t w = 0;
        ok (fc_master_flush (h, out.data(), (std::uint32_t) lat, &w) == FC_OK && w == (std::uint32_t) lat,
            "PRECONDITION: the same call with a separate `written` works");
        fc_master_stats st {}; FC_INIT (st);
        (void) fc_master_get_stats (h, &st);
        ok (st.framesFlushed == (std::uint64_t) lat, "and framesFlushed counts what was drained");
        fc_master_destroy (h);
    }

    group ("solve forwards the CORE's verdict for a degenerate length rather than answering for it");
    {
        fc_master h = make();
        fc_master_params p = goodParams();
        fc_loudness_request req {}; fc_loudness_request_default (&req);
        req.targetLufs = -16.0; req.maxTruePeakDbTp = -1.0;
        auto in = tone (16, kNch);
        std::vector<float> out (in.size(), 0.0f);
        fc_solution sol = 0;
        // `frames == 0` used to be an ABI refusal. The solver has its own answer for it and the contract
        // says FC_OK means a verdict was obtained — answering for the core on one input and forwarding
        // it on every other is two policies for one question.
        ok (fc_master_solve (h, &p, &req, in.data(), out.data(), 0, &sol) == FC_OK,
            "zero frames still produces a VERDICT");
        fc_solution_summary sum {}; FC_INIT (sum);
        ok (fc_solution_summary_get (sol, &sum) == FC_OK && sum.status == FC_SOLVE_INVALID_REQUEST,
            "and the verdict is the solver's InvalidRequest");
        fc_solution_destroy (sol);
        fc_master_destroy (h);
    }


    //==========================================================================
    // The pre-merge diff pass ran its own mutations and eight survived. These are the checks that would
    // have killed them — written from the mutants, not from the code, so each one names what it kills.
    group ("the diff pass's survivors");
    {
        using namespace felitronics::mastering;

        // 1. `measure_lra` must be the CORE's number, not a number this file touched. A mutant adding
        //    0.5 LU to the result survived: the only check on it was a wide interval.
        {
            fc_master h = make();
            const std::size_t frames = (std::size_t) (kFs * 8.0);
            auto in = tone (frames, kNch);
            // THE RANGE LIVES IN CHANNEL 1 ONLY, and channel 0 is a steady tone. With the modulation in
            // both planes, a mutant that read channel 0 for every plane measured the same range and
            // survived — the fixture could not tell "both channels" from "channel 0 twice".
            //
            // FOUR SECONDS A STATE, WHERE IT WAS ONE, and the reason is the instrument, not this check. The
            // short-term window is 3 s: against 1 s states it never sat inside one, so it averaged them and the
            // fixture's range came out 0.5–0.6 LU — a fifth of the 1.9 LU the programme actually swings, and a
            // precondition of "> 0.5" that passed by a single 0.1 LU histogram bin. K12 spent that bin (0.6 → 0.5)
            // and the check went red, which is the margin being one bin, not the cadence being wrong. At 4 s the
            // window fits inside a state, so the measured range IS the swing — 1.9 LU at BOTH cadences, at 8, 12
            // and 20 s — and a mutant that shifts the result by half a LU has somewhere to be seen.
            for (std::size_t i = 0; i < frames; ++i)
                if ((i / 192000) % 2 == 0) in[frames + i] *= 0.05f;
            double viaAbi = 0.0;
            ok (fc_master_measure_lra (h, in.data(), (std::uint32_t) frames, &viaAbi) == FC_OK, "measured");

            TargetLoudnessSolver solver;
            ok (solver.prepare (kFs, kNch, 4096, 256, 4), "a direct solver for the same measurement");
            const float* ip[2] { in.data(), in.data() + frames };
            double direct = 0.0;
            ok (solver.measureInputLoudnessRange (ip, kNch, (int) frames, direct), "the core measured too");
            ok (viaAbi == direct, "and the two are the SAME NUMBER, bit for bit");
            // WHAT THE RANGE IS FOR, stated precisely: the equality above already catches an added offset even
            // on a range of zero, because it compares two measurements of the SAME audio. What needs a range is
            // the OTHER mutant this block exists for — the one that reads channel 0 for every plane. Channel 0 is
            // the steady tone, so that mutant measures 0.0 LU where the truth is 1.9, and the precondition is
            // what guarantees those two are far apart. (The sentence here used to say the range was needed so
            // "half a LU added would show", which is the check that does not need it.)
            ok (direct > 1.5, "PRECONDITION: the fixture swings " + std::to_string (direct)
                              + " LU, so reading channel 0 for both planes would read 0.0 instead");
            fc_master_destroy (h);
        }

        // 2. `framesFlushed` ACCUMULATES. A mutant that made it the size of the LAST flush survived,
        //    because nothing ever flushed twice.
        {
            fc_master h = make();
            fc_master_params p = goodParams();
            fc_master_resolved r {}; FC_INIT (r);
            (void) fc_master_configure (h, &p, &r);
            std::int32_t lat = 0; (void) fc_master_latency (h, &lat);
            std::vector<float> out ((std::size_t) lat * kNch, 0.0f);
            std::uint32_t w = 0;
            ok (fc_master_flush (h, out.data(), (std::uint32_t) lat, &w) == FC_OK, "one flush");
            ok (fc_master_flush (h, out.data(), (std::uint32_t) lat, &w) == FC_OK, "and another");
            fc_master_stats st {}; FC_INIT (st);
            (void) fc_master_get_stats (h, &st);
            ok (st.framesFlushed == (std::uint64_t) lat * 2u,
                "framesFlushed is the RUNNING TOTAL, not the size of the last call");
            fc_master_destroy (h);
        }

        // 3. A solution OUTLIVES the chain that produced it. The header promises it; a mutant that
        //    destroyed every solution along with its master survived, because nothing read one after.
        {
            fc_master h = make();
            fc_master_params p = goodParams();
            fc_loudness_request req {}; fc_loudness_request_default (&req);
            req.targetLufs = -16.0; req.maxTruePeakDbTp = -1.0; req.maxPasses = 1;
            const std::size_t frames = (std::size_t) (kFs * 4.0);
            auto in = tone (frames, kNch);
            std::vector<float> out (in.size(), 0.0f);
            fc_solution sol = 0;
            ok (fc_master_solve (h, &p, &req, in.data(), out.data(), (std::uint32_t) frames, &sol) == FC_OK,
                "a solve");
            fc_solution_summary before {}; FC_INIT (before);
            ok (fc_solution_summary_get (sol, &before) == FC_OK, "readable while the master lives");
            ok (fc_master_destroy (h) == FC_OK, "the master is destroyed");
            fc_solution_summary after {}; FC_INIT (after);
            ok (fc_solution_summary_get (sol, &after) == FC_OK, "and the SOLUTION is still readable");
            ok (after.preLimiterGainDb == before.preLimiterGainDb && after.status == before.status,
                "with the same verdict it had");
            ok (fc_solution_destroy (sol) == FC_OK, "and it is destroyed separately, as documented");
        }

        // 4. The channel weight must REACH the solver. A mutant that dropped the forwarding survived,
        //    because the weight tests asserted statuses and never an effect.
        {
            fc_master_config c = goodConfig();
            c.channels = 3; c.monoBass = 0;                     // mono-bass is stereo-only
            fc_master h = 0;
            ok (fc_master_create (&c, &h) == FC_OK, "a three-channel chain");
            const std::size_t frames = (std::size_t) (kFs * 8.0);
            std::vector<float> in (frames * 3, 0.0f);
            for (std::size_t i = 0; i < frames; ++i)            // content ONLY in channel 2
                in[2 * frames + i] = (float) (((i / 48000) % 2 ? 0.5 : 0.05)
                                              * std::sin (2.0 * 3.14159265358979 * 300.0 * (double) i / kFs));
            double withIt = 0.0;
            ok (fc_master_measure_lra (h, in.data(), (std::uint32_t) frames, &withIt) == FC_OK,
                "PRECONDITION: at weight 1 the programme is measurable");
            ok (withIt > 0.5, "PRECONDITION: and it has a real range to lose");
            ok (fc_master_set_channel_weight (h, 2, 0.0) == FC_OK, "now exclude the only channel with audio");
            double without = -1.0;
            ok (fc_master_measure_lra (h, in.data(), (std::uint32_t) frames, &without) == FC_OK,
                "the call still succeeds");
            ok (without == 0.0 && without != withIt,
                "but the RANGE collapses to zero — so the weight really reached the meter");
            // ⚠️ NB the value it collapses to is `0.0`, which is also what "no dynamic range at all"
            // reports. `measureInputLoudnessRange` cannot tell those apart — the same objection
            // `MasterMeasurement::lraValid` exists for one level up. That is a CORE gap, recorded in
            // the findings rather than papered over here, and this check pins today's behaviour so a
            // future fix to it is a deliberate change rather than a surprise.
            fc_master_destroy (h);
        }

        // 4b. AN EXCLUDED CHANNEL CANNOT POISON A MEASUREMENT IT IS NOT IN. Making the counter
        //     load-bearing for a refusal exposed that `LoudnessMeter` counted poisoned sub-hops
        //     UNWEIGHTED, against its own documentation — so a NaN in a `w = 0` channel, which is what
        //     BS.1770 gives LFE, refused a programme whose weighted energy was perfectly fine.
        {
            fc_master_config c = goodConfig();
            c.channels = 3; c.monoBass = 0;
            fc_master h = 0;
            ok (fc_master_create (&c, &h) == FC_OK, "a three-channel chain");
            const std::size_t frames = (std::size_t) (kFs * 8.0);
            std::vector<float> in (frames * 3, 0.0f);
            for (int ch = 0; ch < 2; ++ch)
                for (std::size_t i = 0; i < frames; ++i)
                    in[(std::size_t) ch * frames + i] =
                        (float) (((i / 48000) % 2 ? 0.5 : 0.05)
                                 * std::sin (2.0 * 3.14159265358979 * (300.0 + 90.0 * ch) * (double) i / kFs));
            in[2 * frames + 12345] = std::numeric_limits<float>::quiet_NaN();   // poison in channel 2 ONLY

            double weighted = 0.0;
            ok (fc_master_measure_lra (h, in.data(), (std::uint32_t) frames, &weighted) == FC_ERR_REFUSED_BY_CORE,
                "PRECONDITION: at weight 1 the poisoned channel refuses the measurement");
            ok (fc_master_set_channel_weight (h, 2, 0.0) == FC_OK, "exclude it, as BS.1770 does for LFE");
            double excluded = 0.0;
            ok (fc_master_measure_lra (h, in.data(), (std::uint32_t) frames, &excluded) == FC_OK,
                "and the measurement is accepted — the poison is in a channel that is not in it");
            ok (excluded > 0.5, "with a real range from the two channels that ARE in it");
            fc_master_destroy (h);
        }

        // 5. The declared check order. A call that is illegal for the HANDLE says so before it says
        //    anything about the structs it carries.
        {
            fc_master h = make();
            fc_master_params p = goodParams();
            fc_master_resolved r {}; FC_INIT (r);
            (void) fc_master_configure (h, &p, &r);
            auto blk = tone (256, kNch);
            (void) fc_master_process (h, blk.data(), blk.data(), 256);
            fc_master_params bad = p; bad.header.abiVersion = 999u;
            ok (fc_master_configure (h, &bad, &r) == FC_ERR_STATE,
                "mid-stream AND mis-versioned answers STATE — the handle's state comes first");
            fc_master_destroy (h);
        }
    }


    //==========================================================================
    group ("the alignment and aliasing guards the pre-merge round found unpinned");
    {
        fc_master h = make();
        fc_master_params p = goodParams();
        fc_master_resolved r {}; FC_INIT (r);
        (void) fc_master_configure (h, &p, &r);

        // A misaligned STRUCT pointer. Only the audio pointers were pinned; the struct and scalar
        // guards could both be deleted without a red test.
        auto raw = std::vector<char> (sizeof (fc_master_params) + 8);
        auto* skewParams = reinterpret_cast<fc_master_params*> (raw.data() + 1);
        ok (fc_master_configure (h, skewParams, &r) == FC_ERR_ALIGNMENT, "a misaligned parameter struct");
        auto rawR = std::vector<char> (sizeof (fc_master_resolved) + 8);
        auto* skewRes = reinterpret_cast<fc_master_resolved*> (rawR.data() + 1);
        ok (fc_master_configure (h, &p, skewRes) == FC_ERR_ALIGNMENT, "a misaligned resolved struct");

        // A misaligned SCALAR out-parameter.
        auto rawL = std::vector<char> (sizeof (std::int32_t) + 8);
        auto* skewLat = reinterpret_cast<std::int32_t*> (rawL.data() + 1);
        ok (fc_master_latency (h, skewLat) == FC_ERR_ALIGNMENT, "a misaligned scalar out-parameter");

        // `flush` with MORE capacity than the latency: the planar stride is `capacity`, not the
        // latency, and nothing exercised the two being different.
        std::int32_t lat = 0; (void) fc_master_latency (h, &lat);
        const std::uint32_t cap = (std::uint32_t) lat + 100u;
        std::vector<float> wide ((std::size_t) cap * kNch, -7.0f);
        std::uint32_t w = 0;
        ok (fc_master_flush (h, wide.data(), cap, &w) == FC_OK && w == (std::uint32_t) lat,
            "a capacity above the latency drains exactly the latency");
        ok (wide[(std::size_t) cap - 1] == -7.0f && wide[(std::size_t) cap + 0] != -7.0f,
            "and the planes are laid out at stride CAPACITY: the gap after plane 0 is untouched "
            "while plane 1 starts at cap");
        fc_master_destroy (h);
    }

    group ("the solve out-handle may not point into anything the call reads");
    {
        fc_master h = make();
        fc_master_params p = goodParams();
        fc_loudness_request req {}; fc_loudness_request_default (&req);
        req.targetLufs = -16.0; req.maxTruePeakDbTp = -1.0; req.maxPasses = 1;
        const std::size_t frames = (std::size_t) (kFs * 2.0);
        auto in = tone (frames, kNch);
        std::vector<float> out (in.size(), 0.0f);

        ok (fc_master_solve (h, &p, &req, in.data(), out.data(), (std::uint32_t) frames,
                             reinterpret_cast<fc_solution*> (out.data())) == FC_ERR_SPAN,
            "an out-handle inside the OUTPUT");
        ok (fc_master_solve (h, &p, &req, in.data(), out.data(), (std::uint32_t) frames,
                             reinterpret_cast<fc_solution*> (const_cast<float*> (in.data() + 8))) == FC_ERR_SPAN,
            "an out-handle inside the INPUT — which used to zero a sample before the search read it");
        ok (in[8] != 0.0f, "and that input sample is still what it was");
        ok (fc_master_solve (h, &p, &req, in.data(), out.data(), (std::uint32_t) frames,
                             reinterpret_cast<fc_solution*> (&p.bypassLimiter)) == FC_ERR_SPAN,
            "and one inside the PARAMETER STRUCT, which used to be zeroed before the mapping read it");
        ok (p.bypassLimiter == goodParams().bypassLimiter, "that field is still what it was too");
        fc_master_destroy (h);
    }

    //==========================================================================
    // P41 — THE BUDGET IS THE ALLOCATION. Every delta is read into a local BEFORE its check: a call's arguments are
    // evaluated in an unspecified order and a message string allocates (gcc builds it first — measured).
    group ("fc_master_need: each budget is exactly what the call allocates");
    {
        // The counter's own rule first (see LoudnessConformanceTests.cpp): around the STL's big-block threshold and far
        // above it, a vector counts as exactly the bytes it asked for. Literal sizes, not the counter's own constants.
        for (const std::size_t nb : { std::size_t { 4095 }, std::size_t { 4096 }, std::size_t { 4097 }, std::size_t { 1 } << 20 })
        {
            const long long got = vectorRequest (nb);
            ok (got == (long long) nb, "the byte counter counts a " + std::to_string (nb) + "-byte vector as "
                                       + std::to_string (nb) + " bytes (read " + std::to_string (got) + ")");
        }
        fc_master h = make();
        fc_master_params p = goodParams();
        fc_master_resolved r {}; FC_INIT (r);
        ok (fc_master_configure (h, &p, &r) == FC_OK, "PRECONDITION: a configured handle");
        const std::uint32_t n = 192000;                  // 4 s at 48 kHz: measure_lra builds a meter only from 3 s
        fc_need solve {}; FC_INIT (solve);
        fc_need lra {};   FC_INIT (lra);
        fc_need bad {};   FC_INIT (bad);
        ok (fc_master_need (h, FC_NEED_SOLVE, n, &solve) == FC_OK, "the solve budget is answered");
        ok (fc_master_need (h, FC_NEED_MEASURE_LRA, n, &lra) == FC_OK, "the measure_lra budget is answered");
        ok (fc_master_need (h, 7, n, &bad) == FC_ERR_ENUM, "an op this ABI does not define is refused");
        // Narrowing before the op, as the check order says: a count past INT_MAX is RANGE whatever the op, and INT_MAX
        // itself is a count the core takes.
        fc_need big {}; FC_INIT (big);
        ok (fc_master_need (h, FC_NEED_SOLVE, 0x80000000u, &big) == FC_ERR_RANGE, "frames past INT_MAX: FC_ERR_RANGE");
        ok (fc_master_need (h, 7, 0x80000000u, &big) == FC_ERR_RANGE, "and RANGE before ENUM: narrowing comes first");
        ok (fc_master_need (h, FC_NEED_SOLVE, 0x7FFFFFFFu, &big) == FC_OK && big.callBytes > 0, "INT_MAX frames are budgeted");
        ok (solve.solverPrepared == 0, "PRECONDITION: the solver is not prepared yet");

        // THE ORACLE, literal on purpose (the one place a restatement is mandatory). 48 kHz stereo, 4 s: the loudness
        // meter is sized for 192000 + 48000 samples = 50 hops of 4800 → 8·(300 + 54 + 58) = 3296 B. The last term
        // is what K12 moved: one short-term sample per hop, the cadence EBU Tech 3342 §3.1 asks for, where it used
        // to be one per ten — so 5 + 8 became 50 + 8, the margin of 8 unchanged. The REFERENCE
        // true-peak meter (P62) is, per channel, one 4x / 32-tap PolyphaseOversampler — prototype 128, phase-major copy
        // 4·32, up ring 2·32, down ring 2·128 floats = 2304 B, plus two int cursors = 2312 B — and one shared scratch of
        // 1024·4 floats = 16 384 B: 2·2312 + 16 384 = 21 008 B. It drains from a fixed array, so there is no drain term
        // (the pre-P62 solve carried a 392 B TruePeakMeter and a 512 B drain buffer).
        // And two traces of 1000 x 32 B: 64 000 B.
        ok (lra.callBytes == 3296u, "the measure_lra budget is the hand-derived 3296 B");
        // And, since the quantile read-back (`fc_solution_gr_quantile`), the window histograms the solution keeps —
        // three of them since K11. THE FIGURE IS PRINTED rather than spelled: the literal that used to stand here
        // (727 960 B) described the two-histogram build and would have gone on reading as a measurement.
        ok (solve.callBytes == 3296u + 21008u + 64000u + kGrWindowBytes,
            "the solve budget is meter + reference true-peak meter + two traces + three window histograms = "
            + std::to_string (solve.callBytes) + " B");

        long long before = alloc::bytes.load();
        const fc_status w = fc_master_set_channel_weight (h, 0, 1.0);
        const long long prepared = alloc::bytes.load() - before;
        ok (w == FC_OK, "the first weight prepares the solver");
        ok (prepared == (long long) solve.solverPrepareBytes, "and allocates exactly `solverPrepareBytes`");
        fc_need after {}; FC_INIT (after);
        ok (fc_master_need (h, FC_NEED_SOLVE, n, &after) == FC_OK && after.solverPrepared == 1,
            "and from then on the budget says the preparation is spent");

        auto in = tone ((int) n, kNch);
        double v = 0.0;
        before = alloc::bytes.load();
        (void) fc_master_measure_lra (h, in.data(), n, &v);
        const long long lraBytes = alloc::bytes.load() - before;
        ok (lraBytes == (long long) lra.callBytes, "measure_lra allocates exactly its budget");

        // Under 3 s there is no range: the call refuses BEFORE it builds a meter, and its budget says so (the
        // code-review round found it promising a meter for exactly this call).
        fc_need shortLra {}; FC_INIT (shortLra);
        ok (fc_master_need (h, FC_NEED_MEASURE_LRA, 48000, &shortLra) == FC_OK && shortLra.callBytes == 0,
            "a 1 s programme: the measure_lra budget is 0");
        before = alloc::bytes.load();
        const fc_status sr = fc_master_measure_lra (h, in.data(), 48000, &v);
        const long long shortBytes = alloc::bytes.load() - before;
        ok (sr == FC_ERR_REFUSED_BY_CORE && shortBytes == 0, "and the refused call allocates nothing");

        // The range rule's own edge, at 48 kHz: exactly 3 s is measurable, one frame less is not — in the call AND in
        // its budget, which read the same `rangeMeasurable`.
        fc_need at3 {}, under3 {}; FC_INIT (at3); FC_INIT (under3);
        ok (fc_master_need (h, FC_NEED_MEASURE_LRA, 144000, &at3) == FC_OK && at3.callBytes > 0,
            "exactly 3 s: a meter is budgeted");
        ok (fc_master_need (h, FC_NEED_MEASURE_LRA, 143999, &under3) == FC_OK && under3.callBytes == 0,
            "one frame under 3 s: nothing is");
        // A solve builds its meters whatever the length — the range rule is NOT the solve's. 1 s still costs
        // meter + reference true-peak meter + two traces: 8·(300 + 24 + 28) + 21 008 + 64 000 = 87 824 B, the
        // short-term term being the one K12 moved (one sample per hop, not one per ten). A length the solve
        // refuses costs 0. THE TOTAL IS PRINTED rather than spelled: the literal that used to close this line
        // (727 696 B) was parameterised on one side and written out on the other, so it went on reading as a
        // measurement after the histograms moved it.
        fc_need s1 {}, s0 {}; FC_INIT (s1); FC_INIT (s0);
        ok (fc_master_need (h, FC_NEED_SOLVE, 48000, &s1) == FC_OK && s1.callBytes == 87824u + kGrWindowBytes,
            "a 1 s solve is budgeted in full: " + std::to_string (s1.callBytes) + " B");
        ok (fc_master_need (h, FC_NEED_SOLVE, 0, &s0) == FC_OK && s0.callBytes == 0,
            "a 0-frame solve, which the core refuses before any pass, costs 0");

        std::vector<float> out (in.size(), 0.0f);
        fc_loudness_request req {}; fc_loudness_request_default (&req);
        req.targetLufs = -14.0; req.maxTruePeakDbTp = -1.0;
        fc_solution sol = 0;
        // THE SOLUTION RECORD IS A PLAIN OBJECT, told to the counter so that MSVC's container padding is never taken off
        // it (as for the instance record above).
        alloc::plainObjectSize.store ((std::size_t) solve.facadeBytes, std::memory_order_relaxed);
        before = alloc::bytes.load();
        const fc_status sv = fc_master_solve (h, &p, &req, in.data(), out.data(), n, &sol);
        const long long solveBytes = alloc::bytes.load() - before;
        alloc::plainObjectSize.store (0, std::memory_order_relaxed);
        fc_solution_summary sum {}; FC_INIT (sum);
        ok (sv == FC_OK && fc_solution_summary_get (sol, &sum) == FC_OK && sum.passes > 0, "PRECONDITION: the search rendered");
        const long long traces = 2LL * 1000LL * 32LL;
        const long long perPass = (long long) solve.callBytes - traces - (long long) kGrWindowBytes;
        ok (solveBytes == (long long) sum.passes * perPass + traces + (long long) kGrWindowBytes
                          + (long long) solve.facadeBytes,
            "a solve allocates passes × (both meters) + its two traces + its window histograms + the facade's record: "
            "its budget's parts");
        (void) fc_solution_destroy (sol);
        (void) fc_master_destroy (h);
    }


    //==========================================================================
    // P41 part 3 — THE CHAIN'S OWN SIDE OF THE FORMULA. `fc_master_need_create` is a DRY RUN of the create
    // (every refusal it can reach before its first allocation, with the same status), and FC_NEED_CONFIGURE
    // is the re-preparation's. Both are held against the counter byte for byte.
    group ("fc_master_need_create: a dry run of the create, and its budget to the byte");
    {
        // THE MATRIX: four rates x three widths x four topologies — the default, the clipper, the key
        // filter, and the legal maximum.
        const double rates[]  = { 44100.0, 48000.0, 96000.0, 192000.0 };
        const std::int32_t widths[] = { 1, 2, 16 };
        // THE COUNTER IS TOLD WHICH SIZE IS A PLAIN OBJECT — see containerBytes(). The number is the core's,
        // read out of a budget rather than written here, so it cannot fall out of step with the facade.
        {
            fc_master_config probe = goodConfig();
            fc_need pn {}; FC_INIT (pn);
            ok (fc_master_need_create (&probe, &pn) == FC_OK && pn.facadeBytes > 0u,
                "PRECONDITION: the facade publishes the size of its own instance record");
            alloc::plainObjectSize.store ((std::size_t) pn.facadeBytes, std::memory_order_relaxed);
            // AND THE RULE IS CALIBRATED, on this row's own STL, rather than trusted: a plain `new` of that
            // size counts as exactly that size, and a VECTOR of it counts as exactly its own bytes. The
            // second is the rule's known collision — a container that happens to be exactly as long as the
            // facade's record would be left uncorrected — and naming it here is what keeps it from being
            // discovered as a byte-for-byte failure with no explanation.
            const std::size_t n = (std::size_t) pn.facadeBytes;
            const long long before = alloc::bytes.load();
            {
                void* raw = ::operator new (n);
                volatile char* sink = static_cast<char*> (raw);
                sink[0] = 1;
                ::operator delete (raw, n);
            }
            const long long got = alloc::bytes.load() - before;
            ok (got == (long long) n, "the counter reads a plain `new` of " + std::to_string (n)
                                      + " B as " + std::to_string (n) + " B (read " + std::to_string (got) + ")");
        }
        int rows = 0, statusOff = 0, bytesOff = 0, zeroBudget = 0, solverOff = 0, cfgOff = 0, headerOff = 0;
        long long worst = 0;
        for (const double fs : rates)
            for (const std::int32_t nch : widths)
                for (int topo = 0; topo < kTopologies; ++topo)
                {
                    fc_master_config c {};
                    fc_master_config_default (&c);
                    c.sampleRate = fs; c.channels = nch;
                    // EVERY OPTIONAL STAGE IS ABSENT ON SOME ROW. With only the clipper moving, a budget
                    // that added the EQ engine's 331 KiB unconditionally — for a chain that never builds
                    // one — was green on every row (the diverse-testing round found it as a surviving
                    // mutation, in the core's matrix; this one carries the same axis for the facade).
                    switch (topo)
                    {
                        case 0: break;                                      // the default
                        case 1: c.clipper = 1; break;
                        case 2: c.sidechainHpfHz = 80.0; break;
                        case 3: c.eq = 0; break;                            // no engine at all
                        case 4: c.compressor = 0; break;
                        case 5: c.limiter = 0; break;
                        case 6: c.dither = 0; break;
                        case 7: c.monoBass = 1; break;                      // stereo only: refused at 1 and 16
                        default:
                            c.clipper = 1; c.sidechainHpfHz = 80.0; c.internalBlock = 8192;
                            c.oversampleFactor = 16; c.tapsPerPhase = 1024;
                            c.compressorLookaheadMs = 250.0; c.limiterLookaheadMs = 20.0;
                            break;
                    }
                    ++rows;
                    fc_need nd {}; FC_INIT (nd);
                    const fc_status ns = fc_master_need_create (&c, &nd);
                    fc_master h = 0;
                    const long long before = alloc::bytes.load();
                    const fc_status cs = fc_master_create (&c, &h);
                    const long long got = alloc::bytes.load() - before;
                    if (ns != cs) ++statusOff;
                    if (cs != FC_OK) continue;
                    if (nd.callBytes == 0u) ++zeroBudget;
                    if (nd.solverPrepareBytes != 0u || nd.solverPrepared != 0) ++solverOff;
                    // THE OUT-STRUCT'S HEADER IS PART OF THE ANSWER. A mutation that returned FC_OK with
                    // `abiVersion = 0` left the whole suite green: every field was checked except the two
                    // that tell a caller which ABI wrote them. (The diverse-testing round.)
                    if (nd.header.abiVersion != FC_MASTER_ABI_VERSION
                        || nd.header.structSize != (std::uint32_t) sizeof (fc_need)) ++headerOff;
                    const long long budget = (long long) (nd.callBytes + nd.facadeBytes);
                    if (got != budget) { ++bytesOff; if (got > worst) worst = got; }

                    // THE RE-PREPARATION, on the same handle: published 0 and measured 0.
                    fc_master_params p = goodParams();
                    fc_master_resolved r {}; FC_INIT (r);
                    fc_need cn {}; FC_INIT (cn);
                    const fc_status cns = fc_master_need (h, FC_NEED_CONFIGURE, 0, &cn);
                    const long long b2 = alloc::bytes.load();
                    const fc_status ccs = fc_master_configure (h, &p, &r);
                    const long long got2 = alloc::bytes.load() - b2;
                    if (cns != FC_OK || ccs != FC_OK || cn.callBytes != 0u || got2 != 0
                        || cn.facadeBytes != 0u || cn.solverPrepareBytes != 0u || cn.solverPrepared != 0) ++cfgOff;
                    (void) fc_master_destroy (h);
                }
        ok (rows == 4 * 3 * kTopologies, "PRECONDITION: 4 rates x 3 widths x " + std::to_string (kTopologies)
            + " topologies (" + std::to_string (rows) + " rows)");
        ok (statusOff == 0, "need_create answers exactly what create answers, on every row ("
                            + std::to_string (statusOff) + " off)");
        ok (bytesOff == 0, "and its budget plus the facade's record is what the create allocates, byte for byte ("
                           + std::to_string (bytesOff) + " rows off, worst " + std::to_string (worst) + " B)");
        ok (zeroBudget == 0, "FC_OK never carries a budget of 0: the number means ONE thing on this op");
        ok (solverOff == 0, "the solver's fields come back neutral for a create");
        ok (headerOff == 0, "and the budget it wrote carries THIS build's header (" + std::to_string (headerOff) + " off)");
        ok (cfgOff == 0, "a configure is published as 0 and allocates 0, on every row ("
                         + std::to_string (cfgOff) + " off)");

        // THE REFUSALS, each with the status the create gives it — and NOTHING allocated on the way to it.
        // Every one of these used to cost between 51 288 and 394 456 bytes before it said no at this
        // geometry, and up to 1 668 312 at sixteen channels and an 8192-sample quantum.
        struct Refusal { const char* what; fc_status want; void (*edit) (fc_master_config&); };
        const Refusal refusals[] = {
            { "a 20 Hz rate",                      FC_ERR_REFUSED_BY_CORE, [] (fc_master_config& c) { c.sampleRate = 20.0; } },
            { "a rate of zero",                    FC_ERR_REFUSED_BY_CORE, [] (fc_master_config& c) { c.sampleRate = 0.0; } },
            { "a NaN rate",                        FC_ERR_NON_FINITE,      [] (fc_master_config& c) { c.sampleRate = std::numeric_limits<double>::quiet_NaN(); } },
            { "no channels",                       FC_ERR_REFUSED_BY_CORE, [] (fc_master_config& c) { c.channels = 0; } },
            { "more channels than the core has",   FC_ERR_REFUSED_BY_CORE, [] (fc_master_config& c) { c.channels = 999; } },
            { "a 300 ms compressor lookahead",     FC_ERR_REFUSED_BY_CORE, [] (fc_master_config& c) { c.compressorLookaheadMs = 300.0; } },
            { "2000 taps per phase",               FC_ERR_REFUSED_BY_CORE, [] (fc_master_config& c) { c.tapsPerPhase = 2000; } },
            { "a factor past the limiter's 16",    FC_ERR_REFUSED_BY_CORE, [] (fc_master_config& c) { c.oversampleFactor = 32; } },
            { "mono-bass on a mono chain",         FC_ERR_REFUSED_BY_CORE, [] (fc_master_config& c) { c.channels = 1; c.monoBass = 1; } },
            { "a quantum under the floor",         FC_ERR_REFUSED_BY_CORE, [] (fc_master_config& c) { c.internalBlock = 4; } },
            // The boundaries the diverse-testing round found missing: every one of them is a refusal the
            // core states somewhere, and none of them had a row here proving the two entry points agree on
            // it or that nothing was allocated reaching it.
            { "a rate past the core's ceiling",     FC_ERR_REFUSED_BY_CORE, [] (fc_master_config& c) { c.sampleRate = 3.0e6 + 1.0; } },
            { "a quantum past the ceiling",         FC_ERR_REFUSED_BY_CORE, [] (fc_master_config& c) { c.internalBlock = 8193; } },
            { "an oversampling factor of 1",        FC_ERR_REFUSED_BY_CORE, [] (fc_master_config& c) { c.oversampleFactor = 1; } },
            { "three taps per phase",               FC_ERR_REFUSED_BY_CORE, [] (fc_master_config& c) { c.tapsPerPhase = 3; } },
            { "a negative compressor lookahead",    FC_ERR_REFUSED_BY_CORE, [] (fc_master_config& c) { c.compressorLookaheadMs = -1.0; } },
            { "a negative limiter lookahead",       FC_ERR_REFUSED_BY_CORE, [] (fc_master_config& c) { c.limiterLookaheadMs = -1.0; } },
            { "a negative key filter",              FC_ERR_REFUSED_BY_CORE, [] (fc_master_config& c) { c.sidechainHpfHz = -1.0; } },
            { "a key filter above Nyquist",         FC_ERR_REFUSED_BY_CORE, [] (fc_master_config& c) { c.sidechainHpfHz = 0.5 * kFs; } },
            { "a non-finite key filter",            FC_ERR_NON_FINITE,      [] (fc_master_config& c) { c.sidechainHpfHz = std::numeric_limits<double>::infinity(); } },
            // P51 — THE RATE FLOOR, 8000 Hz, and the one row that sees it alone: with the limiter and the EQ off no
            // stage refuses a low rate, so a create that took one ulp under the floor would have nothing else to
            // stop it. The rest are the mistakes the floor is for; every one of them was CREATED on origin/main.
            { "one ulp under the rate floor, no stage that refuses", FC_ERR_REFUSED_BY_CORE,
              [] (fc_master_config& c) { c.sampleRate = std::nextafter (8000.0, 0.0); c.eq = 0; c.limiter = 0; } },
            { "7999 Hz",                            FC_ERR_REFUSED_BY_CORE, [] (fc_master_config& c) { c.sampleRate = 7999.0; } },
            { "88.2 — kilohertz passed as hertz",   FC_ERR_REFUSED_BY_CORE, [] (fc_master_config& c) { c.sampleRate = 88.2; } },
            { "44.1 with the limiter off",          FC_ERR_REFUSED_BY_CORE, [] (fc_master_config& c) { c.sampleRate = 44.1; c.limiter = 0; } },
            { "3300 Hz, where the K-weighting is past Nyquist", FC_ERR_REFUSED_BY_CORE, [] (fc_master_config& c) { c.sampleRate = 3300.0; } },
            // A DELIVERING handle runs its chain at the DELIVERY rate, so the chain's floor never sees the source —
            // the resampler's plan is what refuses it, on both sides. THE STAMP IS RE-WRITTEN FIRST: goodConfig() is a
            // v1 struct, and `deliveryRate` past a v1 stamp is not read (the trap pinned in the v2 group below) — these
            // rows were plain handles, refused or not for the wrong reason, until they said so.
            { "delivering from 7999 Hz",            FC_ERR_REFUSED_BY_CORE, [] (fc_master_config& c) { FC_INIT (c); c.sampleRate = 7999.0; c.deliveryRate = 48000.0; } },
            { "delivering from 4000 Hz",            FC_ERR_REFUSED_BY_CORE, [] (fc_master_config& c) { FC_INIT (c); c.sampleRate = 4000.0; c.deliveryRate = 48000.0; } },
            { "delivering to 7999 Hz",              FC_ERR_REFUSED_BY_CORE, [] (fc_master_config& c) { FC_INIT (c); c.deliveryRate = 7999.0; } },
        };
        int refOff = 0, refLeak = 0, refTouched = 0;
        for (const Refusal& rf : refusals)
        {
            fc_master_config c = goodConfig();
            rf.edit (c);
            // A SENTINEL, NOT ZEROS. "The budget is left untouched" is not proved by reading 0 out of a
            // struct that was 0 going in — an implementation that wrote zeros, or re-stamped the header,
            // would pass. The whole struct is filled with a pattern, given a valid input header, copied,
            // and compared byte for byte afterwards. (The diverse-testing round.)
            fc_need nd;
            std::memset (&nd, 0xA5, sizeof (nd));
            FC_INIT (nd);
            fc_need before_nd = nd;
            const fc_status ns = fc_master_need_create (&c, &nd);
            fc_master h = 0xDEADBEEFu;
            const fc_master before_h = h;
            const long long before = alloc::bytes.load();
            const fc_status cs = fc_master_create (&c, &h);
            const long long got = alloc::bytes.load() - before;
            if (ns != rf.want || cs != rf.want) ++refOff;
            if (std::memcmp (&nd, &before_nd, sizeof (nd)) != 0 || h != before_h) ++refTouched;
            if (got != 0) ++refLeak;
            if (cs == FC_OK) (void) fc_master_destroy (h);
        }
        ok (refOff == 0, "every geometry the core refuses is refused by BOTH, with the same status ("
                         + std::to_string (refOff) + " off, over " + std::to_string (sizeof (refusals) / sizeof (refusals[0]))
                         + " refusals)");
        ok (refTouched == 0, "and neither out-argument is touched, byte for byte, by either call ("
                             + std::to_string (refTouched) + " touched)");
        ok (refLeak == 0, "and a refused create allocates NOTHING — it used to ask for 394 456 bytes on its "
                          "way to `false` here, and 1 668 312 at sixteen channels ("
                          + std::to_string (refLeak) + " leaked)");

        // AND THE FLOOR ITSELF IS A RATE, on every kind of handle — `>=`, not `>`. Both calls, and both say OK.
        {
            struct Accept { const char* what; double sr, dr; };
            const Accept accepts[] = { { "a plain handle at 8000 Hz", 8000.0, 0.0 },
                                       { "delivering 8000 -> 48000", 8000.0, 48000.0 },
                                       { "delivering 48000 -> 8000", 48000.0, 8000.0 },
                                       { "delivering 8000 -> 8000", 8000.0, 8000.0 } };
            for (const Accept& a : accepts)
            {
                fc_master_config c = goodConfig();
                FC_INIT (c);                                        // a stamp under which `deliveryRate` is read
                c.sampleRate = a.sr; c.deliveryRate = a.dr;
                fc_need an {}; FC_INIT (an);
                fc_master h = 0;
                const bool okNeed = fc_master_need_create (&c, &an) == FC_OK && an.callBytes > 0u;
                const bool okMade = fc_master_create (&c, &h) == FC_OK;
                ok (okNeed && okMade, std::string ("accepted at the floor: ") + a.what);
                // and it is the KIND of handle the row names: a delivering one answers the delivered length — one
                // second in, one second out
                std::uint32_t d = 0;
                const fc_status ds = okMade ? fc_master_delivered_frames (h, (std::uint32_t) a.sr, &d) : FC_ERR_HANDLE;
                ok (a.dr == 0.0 ? ds == FC_ERR_STATE : (ds == FC_OK && d == (std::uint32_t) a.dr),
                    std::string ("and it is the handle the row names: ") + a.what);
                if (okMade) (void) fc_master_destroy (h);
            }
        }

        // THE ARGUMENT CHECKS ARE THIS CALL'S OWN, and they come before the core's — the same order the
        // create takes.
        fc_master_config gc = goodConfig();
        fc_need nd {}; FC_INIT (nd);
        ok (fc_master_need_create (&gc, nullptr) == FC_ERR_NULL, "a null budget: FC_ERR_NULL");
        ok (fc_master_need_create (nullptr, &nd) == FC_ERR_NULL, "a null config: FC_ERR_NULL");
        { fc_need bad {}; FC_INIT (bad); bad.header.abiVersion = FC_MASTER_ABI_VERSION + 1u;
          ok (fc_master_need_create (&gc, &bad) == FC_ERR_ABI_VERSION, "a budget struct from another ABI: refused"); }
        { fc_master_config bc = gc; bc.header.structSize = 3u;
          ok (fc_master_need_create (&bc, &nd) == FC_ERR_STRUCT_SIZE, "a config of the wrong size: refused"); }
        { fc_master_config bc = gc; bc.sidechainHpfHz = std::numeric_limits<double>::quiet_NaN();
          ok (fc_master_need_create (&bc, &nd) == FC_ERR_NON_FINITE,
              "a field the mapping rejects: its own status, ahead of any geometry"); }

        // THE TABLE IS PART OF THE ANSWER. With every slot taken the create answers FC_ERR_EXHAUSTED, and a
        // budget published for it would be a number for a call that cannot be made.
        fc_master live[8] {};
        int made = 0;
        for (int i = 0; i < 8; ++i) if (fc_master_create (&gc, &live[i]) == FC_OK) ++made;
        ok (made == 8, "PRECONDITION: the handle table is full (" + std::to_string (made) + " live)");
        fc_need full {}; FC_INIT (full);
        fc_master extra = 0;
        ok (fc_master_need_create (&gc, &full) == FC_ERR_EXHAUSTED && full.callBytes == 0u,
            "with no free slot the budget is refused, exactly as the create is");
        ok (fc_master_create (&gc, &extra) == FC_ERR_EXHAUSTED, "PRECONDITION: and the create really is refused");
        for (int i = 0; i < made; ++i) (void) fc_master_destroy (live[i]);

        // FC_NEED_CONFIGURE's OWN ARGUMENT: a re-preparation has no programme length, so a non-zero count is
        // refused rather than ignored — narrowing before field values, as everywhere here.
        fc_master h = make();
        fc_master_params p = goodParams();
        fc_master_resolved r {}; FC_INIT (r);
        ok (fc_master_configure (h, &p, &r) == FC_OK, "PRECONDITION: a configured handle");
        fc_need cn {}; FC_INIT (cn);
        ok (fc_master_need (h, FC_NEED_CONFIGURE, 1, &cn) == FC_ERR_RANGE, "a frame count with a configure: FC_ERR_RANGE");
        ok (fc_master_need (h, FC_NEED_CONFIGURE, 0, &cn) == FC_OK && cn.callBytes == 0u, "and 0 is the count it takes");

        // THE DEMAND ANSWERS A NUMBER, NEVER A PERMISSION (the decision recorded with this work). A stream in
        // progress makes the configure itself FC_ERR_STATE; its COST is unchanged by the moment, and a page
        // deciding whether to reset and re-configure needs the number exactly then. `fc_master_need` already
        // behaves this way for a solve, and this is the same rule, not a second one.
        auto audio = tone (256, kNch);
        std::vector<float> outBuf (audio.size(), 0.0f);
        ok (fc_master_process (h, audio.data(), outBuf.data(), 256) == FC_OK, "PRECONDITION: a stream is in progress");
        ok (fc_master_configure (h, &p, &r) == FC_ERR_STATE, "PRECONDITION: the configure itself is refused now");
        fc_need mid {}; FC_INIT (mid);
        ok (fc_master_need (h, FC_NEED_CONFIGURE, 0, &mid) == FC_OK && mid.callBytes == 0u,
            "and its budget is still answered: the demand is a number, not a permission");
        fc_need sd {}; FC_INIT (sd);
        ok (fc_master_need (h, FC_NEED_SOLVE, 48000, &sd) == FC_OK && sd.callBytes > 0u,
            "PRECONDITION: which is what a solve's budget already does mid-stream");
        (void) fc_master_destroy (h);
        alloc::plainObjectSize.store (0, std::memory_order_relaxed);   // the exemption is this group's only
    }

    //==========================================================================
    // P57b — ABI v2, THE COMPATIBILITY RULE, held against bytes rather than read off the header: which versions are
    // read, how many bytes are read and written, that a caller's header comes back as it went in, and that the
    // frozen writers stay frozen. Every buffer below is a struct followed by CANARY bytes, because "wrote nothing
    // past the caller's size" is only a claim until something past it is looked at.
    group ("the version rule — what is read, what is written, and nothing past the caller's size");
    {
        const std::uint32_t kCur = FC_MASTER_ABI_VERSION;
        ok (kCur == 12u, "PRECONDITION: this group is written for v12 (v2: deliveryRate; v3: compressorMix; v4: the GR trace; "
                         "v5: fc_master_set_progress, no struct grew; v6: M2 — params, resolved and request grew; "
                         "v7: fc_master_eq_curve, no struct grew; v8: the request's two quantiles and "
                         "fc_solution_gr_quantile; v9: fc_master_eq_dyn_times, no struct grew; v10: K11 — the request's "
                         "limiter input gate, fc_gr_active_stats and fc_solution_gr_active_stats; v11: K13 — the peak "
                         "clipper inside the limiter, params/resolved/measurement grew; v12: K14 — the Side air "
                         "shelf: config, params, resolved and measurement grew)");

        // THE TABLE (rule 5), every (struct, version) pair of today.
        ok (fc_master_sizeof (FC_STRUCT_CONFIG, 1) == 80u && fc_master_sizeof (FC_STRUCT_CONFIG, 2) == 88u
            && fc_master_sizeof (FC_STRUCT_CONFIG, 3) == 88u && fc_master_sizeof (FC_STRUCT_CONFIG, 4) == 88u
            && fc_master_sizeof (FC_STRUCT_CONFIG, 5) == 88u && fc_master_sizeof (FC_STRUCT_CONFIG, 6) == 88u
            && fc_master_sizeof (FC_STRUCT_CONFIG, 7) == 88u && fc_master_sizeof (FC_STRUCT_CONFIG, 8) == 88u
            && fc_master_sizeof (FC_STRUCT_CONFIG, 9) == 88u && fc_master_sizeof (FC_STRUCT_CONFIG, 10) == 88u,
            "config: 80 at v1, 88 from v2");
        ok (fc_master_sizeof (FC_STRUCT_PARAMS, 1) == 6560u && fc_master_sizeof (FC_STRUCT_PARAMS, 2) == 6560u
            && fc_master_sizeof (FC_STRUCT_PARAMS, 3) == 6568u && fc_master_sizeof (FC_STRUCT_PARAMS, 4) == 6568u
            && fc_master_sizeof (FC_STRUCT_PARAMS, 5) == 6568u && fc_master_sizeof (FC_STRUCT_PARAMS, 6) == 6584u
            && fc_master_sizeof (FC_STRUCT_PARAMS, 7) == 6584u && fc_master_sizeof (FC_STRUCT_PARAMS, 8) == 6584u
            && fc_master_sizeof (FC_STRUCT_PARAMS, 9) == 6584u && fc_master_sizeof (FC_STRUCT_PARAMS, 10) == 6584u,
            "params: 6560 at v1 and v2, 6568 from v3, 6584 from v6");
        ok (fc_master_sizeof (FC_STRUCT_RESOLVED, 1) == 80u && fc_master_sizeof (FC_STRUCT_RESOLVED, 2) == 80u
            && fc_master_sizeof (FC_STRUCT_RESOLVED, 3) == 88u && fc_master_sizeof (FC_STRUCT_RESOLVED, 4) == 88u
            && fc_master_sizeof (FC_STRUCT_RESOLVED, 5) == 88u && fc_master_sizeof (FC_STRUCT_RESOLVED, 6) == 96u
            && fc_master_sizeof (FC_STRUCT_RESOLVED, 7) == 96u && fc_master_sizeof (FC_STRUCT_RESOLVED, 8) == 96u
            && fc_master_sizeof (FC_STRUCT_RESOLVED, 9) == 96u && fc_master_sizeof (FC_STRUCT_RESOLVED, 10) == 96u,
            "resolved: 80 at v1 and v2, 88 from v3, 96 from v6");
        ok (fc_master_sizeof (FC_STRUCT_MEASUREMENT, 1) == 208u && fc_master_sizeof (FC_STRUCT_MEASUREMENT, 3) == 208u
            && fc_master_sizeof (FC_STRUCT_MEASUREMENT, 4) == 224u && fc_master_sizeof (FC_STRUCT_MEASUREMENT, 5) == 224u
            && fc_master_sizeof (FC_STRUCT_MEASUREMENT, 6) == 224u && fc_master_sizeof (FC_STRUCT_MEASUREMENT, 7) == 224u
            && fc_master_sizeof (FC_STRUCT_MEASUREMENT, 8) == 224u && fc_master_sizeof (FC_STRUCT_MEASUREMENT, 9) == 224u && fc_master_sizeof (FC_STRUCT_MEASUREMENT, 10) == 224u,
            "measurement: 208 to v3, 224 from v4 — the trace's four fields");
        ok (fc_master_sizeof (FC_STRUCT_REQUEST, 1) == 120u && fc_master_sizeof (FC_STRUCT_REQUEST, 5) == 120u
            && fc_master_sizeof (FC_STRUCT_REQUEST, 6) == 128u && fc_master_sizeof (FC_STRUCT_REQUEST, 7) == 128u
            && fc_master_sizeof (FC_STRUCT_REQUEST, 8) == 144u && fc_master_sizeof (FC_STRUCT_REQUEST, 9) == 144u && fc_master_sizeof (FC_STRUCT_REQUEST, 10) == 152u,
            "request: 120 to v5, 128 from v6 — `grTraceBuckets` and its named padding — 144 from v8, the two quantiles");
        int inherit = 0;
        for (int id = FC_STRUCT_STATS; id <= FC_STRUCT_SUMMARY; ++id)
        {
            if (id == FC_STRUCT_MEASUREMENT || id == FC_STRUCT_REQUEST) continue;
            for (std::uint32_t v = 2u; v <= kCur; ++v)
                if (fc_master_sizeof (id, 1) == 0u || fc_master_sizeof (id, v) != fc_master_sizeof (id, 1)) ++inherit;
        }
        ok (inherit == 0, "every struct that did not grow inherits its v1 row, at every version");
        ok (fc_master_sizeof (FC_STRUCT_CONFIG, 0) == 0u && fc_master_sizeof (FC_STRUCT_CONFIG, kCur + 1u) == 0u
            && fc_master_sizeof (99, 1) == 0u && fc_master_sizeof (-1, 1) == 0u,
            "and 0 for a version or an id it does not have");

        // IN (rule 6): the version x size matrix, through both calls that read a config.
        struct Row { std::uint32_t v, s; fc_status want; const char* what; };
        const Row rows[] = {
            { 1u, 80u, FC_OK,              "v1 at 80 bytes" },
            { 2u, 88u, FC_OK,              "v2 at 88 bytes" },
            { 3u, 88u, FC_OK,              "v3 at 88 bytes — the config did not grow at v3" },
            { 4u, 88u, FC_OK,              "v4 at 88 bytes — nor at v4" },
            { 5u, 88u, FC_OK,              "v5 at 88 bytes — nor at v5" },
            { 6u, 88u, FC_OK,              "v6 at 88 bytes — nor at v6" },
            { 7u, 88u, FC_OK,              "v7 at 88 bytes — nor at v7" },
            { 8u, 88u, FC_OK,              "v8 at 88 bytes — nor at v8" },
            { 9u, 88u, FC_OK,              "v9 at 88 bytes — nor at v9" },
            { 10u, 88u, FC_OK,             "v10 at 88 bytes — nor at v10" },
            { 1u, 88u, FC_ERR_STRUCT_SIZE, "v1 claiming v2's size" },
            { 2u, 80u, FC_ERR_STRUCT_SIZE, "v2 claiming v1's size" },
            { 0u, 80u, FC_ERR_ABI_VERSION, "version 0" },
            { kCur + 1u, 96u, FC_ERR_ABI_VERSION, "a version newer than this build — refused by decision" },
        };
        for (const Row& rw : rows)
        {
            fc_master_config c = goodConfig();
            c.header.abiVersion = rw.v; c.header.structSize = rw.s;
            fc_need nd {}; FC_INIT (nd);
            fc_master h = 0;
            const fc_status ns = fc_master_need_create (&c, &nd);
            const fc_status cs = fc_master_create (&c, &h);
            ok (ns == rw.want && cs == rw.want, std::string ("config ") + rw.what);
            if (cs == FC_OK) (void) fc_master_destroy (h);
        }

        // ONLY THE CALLER'S BYTES ARE READ. A v1 config with a delivery rate sitting in the bytes after its 80 — which
        // is exactly what a caller that filled the struct from a FROZEN writer and then set the v2 field has — makes
        // a handle that does NOT deliver. The same bytes under a v2 stamp make one that does.
        {
            fc_master_config c = goodConfig();                              // the frozen writer: a v1 stamp
            ok (c.header.abiVersion == 1u && c.header.structSize == 80u, "PRECONDITION: goodConfig() is a v1 struct");
            c.deliveryRate = 96000.0;
            fc_master h = 0; std::uint32_t d = 0;
            ok (fc_master_create (&c, &h) == FC_OK && fc_master_delivered_frames (h, 1000u, &d) == FC_ERR_STATE,
                "THE TRAP, pinned: `_default` then `deliveryRate` — the field lies past the stamp and is not read");
            (void) fc_master_destroy (h);
            FC_INIT (c);
            ok (fc_master_create (&c, &h) == FC_OK && fc_master_delivered_frames (h, 48000u, &d) == FC_OK && d == 96000u,
                "and under a v2 stamp the same bytes make a delivering handle");
            (void) fc_master_destroy (h);
        }

        // THE SAME, FOR PARAMETERS, and here the new field has a value that shows: a v2 parameter set carrying a mix
        // past its 6560 bytes renders at mix 1 — the read-back says so — and the same bytes under v3 apply 0.25.
        {
            fc_master h = make();
            fc_master_resolved r {}; FC_INIT (r);
            fc_master_params p2 = goodParams();
            p2.header.abiVersion = 2u; p2.header.structSize = 6560u;
            p2.compressorMix = 0.25;
            ok (fc_master_configure (h, &p2, &r) == FC_OK && r.compressorMix == 1.0,
                "a v2 parameter set: the mix past its 6560 bytes is not read — the chain applies 1");
            fc_master_params p3 = p2; FC_INIT (p3);
            ok (fc_master_configure (h, &p3, &r) == FC_OK && r.compressorMix == 0.25,
                "and under v3 the same bytes apply 0.25");
            (void) fc_master_destroy (h);
        }

        // THE CALLER'S SIZE BOUNDS EVERY CHECK, not only the copy: an out-handle placed right after a v2 parameter set
        // is outside it, and a v3 caller whose struct is 8 bytes longer has it inside.
        {
            fc_master h = make();
            std::vector<std::uint64_t> buf (6568u / 8u + 4u, 0u);
            fc_master_params base = goodParams();
            std::memcpy (buf.data(), &base, 6560u);
            auto* pv = reinterpret_cast<fc_master_params*> (static_cast<void*> (buf.data()));
            pv->header.abiVersion = 2u; pv->header.structSize = 6560u;
            auto* after = reinterpret_cast<fc_solution*> (static_cast<void*> ((char*) buf.data() + 6560));
            fc_loudness_request req {}; fc_loudness_request_default (&req);
            req.targetLufs = -14.0; req.maxTruePeakDbTp = -1.0;
            auto in = tone (4800, kNch);
            std::vector<float> out (in.size(), 0.0f);
            const fc_status s2 = fc_master_solve (h, pv, &req, in.data(), out.data(), 4800u, after);
            ok (s2 == FC_OK, "an out-handle right after a v2 parameter set is NOT inside it: FC_OK");
            if (s2 == FC_OK) (void) fc_solution_destroy (*after);
            (void) fc_master_destroy (h);
            h = make();
            pv->header.abiVersion = 3u; pv->header.structSize = 6568u;
            ok (fc_master_solve (h, pv, &req, in.data(), out.data(), 4800u, after) == FC_ERR_SPAN,
                "the same address inside a v3 parameter set: SPAN");
            (void) fc_master_destroy (h);
        }

        // OUT (rule 7): the header is an echo, and nothing past the caller's size is written — at two sizes where a
        // struct has two (`fc_master_resolved`: 80 bytes to v2, 88 from v3) and at the one size of the others.
        auto canaried = [] (std::size_t bytes) { std::vector<unsigned char> b (bytes + 64u, 0xC3); return b; };
        auto intact   = [] (const std::vector<unsigned char>& b, std::size_t from)
        { for (std::size_t i = from; i < b.size(); ++i) if (b[i] != 0xC3) return false; return true; };
        {
            fc_master h = make();
            fc_master_params p = goodParams();
            auto rb = canaried (88u);
            auto* r = reinterpret_cast<fc_master_resolved*> (rb.data());
            r->header.abiVersion = 1u; r->header.structSize = 80u;
            ok (fc_master_configure (h, &p, r) == FC_OK && r->header.abiVersion == 1u && r->header.structSize == 80u
                && intact (rb, 80u) && r->latencySamples > 0,
                "configure: a v1 resolved comes back stamped v1 — not this build's version — with nothing written past "
                "its 80, where this build's `compressorMix` would go");
            auto r3b = canaried (88u);
            auto* r3 = reinterpret_cast<fc_master_resolved*> (r3b.data());
            r3->header.abiVersion = 3u; r3->header.structSize = 88u;
            ok (fc_master_configure (h, &p, r3) == FC_OK && r3->header.abiVersion == 3u && r3->compressorMix == 1.0
                && intact (r3b, 88u), "and a v3 one gets its 88 — `compressorMix` included — and not a byte more");
            auto rg = canaried (80u);
            auto* r2 = reinterpret_cast<fc_master_resolved*> (rg.data());
            r2->header.abiVersion = 1u; r2->header.structSize = 80u;
            ok (fc_master_resolved_get (h, r2) == FC_OK && r2->header.abiVersion == 1u && intact (rg, 80u),
                "resolved_get: the same");
            auto sb = canaried (32u);
            auto* st = reinterpret_cast<fc_master_stats*> (sb.data());
            st->header.abiVersion = 1u; st->header.structSize = 32u;
            ok (fc_master_get_stats (h, st) == FC_OK && st->header.abiVersion == 1u && intact (sb, 32u), "get_stats: the same");
            auto nb = canaried (40u);
            auto* nd = reinterpret_cast<fc_need*> (nb.data());
            nd->header.abiVersion = 1u; nd->header.structSize = 40u;
            ok (fc_master_need (h, FC_NEED_SOLVE, 48000u, nd) == FC_OK && nd->header.abiVersion == 1u
                && nd->_pad0 == 0 && nd->callBytes > 0u && intact (nb, 40u),
                "need: the same, and the named padding is written 0");
            auto in = tone (48000, kNch);
            std::vector<float> out (in.size(), 0.0f);
            fc_loudness_request req {}; fc_loudness_request_default (&req);
            req.targetLufs = -14.0; req.maxTruePeakDbTp = -1.0;
            fc_solution sol = 0;
            ok (fc_master_solve (h, &p, &req, in.data(), out.data(), 48000u, &sol) == FC_OK, "PRECONDITION: a solution");
            auto mb = canaried (208u);
            auto* ms = reinterpret_cast<fc_measurement*> (mb.data());
            ms->header.abiVersion = 1u; ms->header.structSize = 208u;
            ok (fc_solution_measurement (sol, ms) == FC_OK && ms->header.abiVersion == 1u && intact (mb, 208u),
                "solution_measurement: the same");
            auto mb3 = canaried (224u);
            auto* m3 = reinterpret_cast<fc_measurement*> (mb3.data());
            m3->header.abiVersion = 3u; m3->header.structSize = 208u;
            ok (fc_solution_measurement (sol, m3) == FC_OK && m3->header.abiVersion == 3u && intact (mb3, 208u),
                "a v3 measurement: 208 bytes, and the 16 where v4 keeps the trace's counts untouched");
            auto mb4 = canaried (224u);
            auto* m4 = reinterpret_cast<fc_measurement*> (mb4.data());
            m4->header.abiVersion = 4u; m4->header.structSize = 224u;
            ok (fc_solution_measurement (sol, m4) == FC_OK && m4->header.abiVersion == 4u && intact (mb4, 224u)
                && m4->limiterGrTraceBuckets == 1000 && m4->compressorGrTraceBuckets == 1000
                && m4->limiterGrTraceValid == 1 && m4->compressorGrTraceValid == 1,
                "and a v4 one gets its 224 — 1000 buckets per stage for 48000 frames, both valid — and not a byte more");
            auto ub = canaried (88u);
            auto* su = reinterpret_cast<fc_solution_summary*> (ub.data());
            su->header.abiVersion = 1u; su->header.structSize = 88u;
            ok (fc_solution_summary_get (sol, su) == FC_OK && su->header.abiVersion == 1u && su->passes > 0 && intact (ub, 88u),
                "solution_summary_get: the same");
            (void) fc_solution_destroy (sol);
            (void) fc_master_destroy (h);
        }

        // THE FROZEN WRITERS (rule 8) write exactly v1's bytes and a v1 stamp — COUNTED, not inferred from a formula.
        {
            auto cb = canaried (88u);
            fc_master_config_default (reinterpret_cast<fc_master_config*> (cb.data()));
            auto* c = reinterpret_cast<fc_master_config*> (cb.data());
            std::size_t firstUntouched = cb.size();
            for (std::size_t i = cb.size(); i-- > 0;) if (cb[i] != 0xC3) { firstUntouched = i + 1; break; }
            ok (c->header.abiVersion == 1u && c->header.structSize == 80u && firstUntouched <= 80u && intact (cb, 80u),
                "fc_master_config_default writes a v1 stamp and no byte at or past 80 — though this build's config is 88 ("
                + std::to_string (firstUntouched) + ")");
            auto pb = canaried (6568u);
            fc_master_params_default (reinterpret_cast<fc_master_params*> (pb.data()));
            ok (reinterpret_cast<fc_master_params*> (pb.data())->header.abiVersion == 1u && intact (pb, 6560u),
                "fc_master_params_default: v1, 6560 bytes and none past — though this build's params are 6584");
            auto qb = canaried (120u);
            fc_loudness_request_default (reinterpret_cast<fc_loudness_request*> (qb.data()));
            ok (reinterpret_cast<fc_loudness_request*> (qb.data())->header.abiVersion == 1u && intact (qb, 120u),
                "fc_loudness_request_default: v1, 120 bytes and none past");
        }

        // THE VERSIONED WRITERS write the caller's version and nothing else.
        {
            auto cb = canaried (88u);
            auto* c = reinterpret_cast<fc_master_config*> (cb.data());
            std::memset (cb.data(), 0, 8);
            ok (fc_master_config_defaults (c) == FC_ERR_ABI_VERSION, "an unstamped struct is refused: the writer cannot know its size");
            c->header.abiVersion = 1u; c->header.structSize = 80u;
            ok (fc_master_config_defaults (c) == FC_OK && c->header.abiVersion == 1u && c->header.structSize == 80u
                && intact (cb, 80u), "a v1 stamp: 80 bytes written, the header echoed, the canaries after 80 intact");
            fc_master_config frozen {};
            fc_master_config_default (&frozen);
            ok (std::memcmp ((const char*) c + 8, (const char*) &frozen + 8, 72u) == 0,
                "and those 80 bytes are the frozen writer's 80");
            auto vb = canaried (88u);
            auto* v2 = reinterpret_cast<fc_master_config*> (vb.data());
            v2->header.abiVersion = 2u; v2->header.structSize = 88u;
            ok (fc_master_config_defaults (v2) == FC_OK && v2->deliveryRate == 0.0 && intact (vb, 88u)
                && std::memcmp ((const char*) v2 + 8, (const char*) &frozen + 8, 72u) == 0,
                "a v2 stamp: 88 bytes, `deliveryRate` = 0 — the value under which v2 is v1 — and nothing past");
            v2->header.structSize = 80u;
            ok (fc_master_config_defaults (v2) == FC_ERR_STRUCT_SIZE, "a v2 stamp with v1's size: refused");
            ok (fc_master_config_defaults (nullptr) == FC_ERR_NULL, "null: refused");
            fc_master_params pv {}; FC_INIT (pv);
            fc_loudness_request qv {}; FC_INIT (qv);
            ok (fc_master_params_defaults (&pv) == FC_OK && fc_loudness_request_defaults (&qv) == FC_OK
                && std::isnan (qv.targetLufs) && pv.dither.bits == goodParams().dither.bits && pv.compressorMix == 1.0,
                "the params and request writers carry the core's defaults too — `compressorMix` 1 among them");
            auto pb2 = canaried (6568u);
            auto* p2 = reinterpret_cast<fc_master_params*> (pb2.data());
            p2->header.abiVersion = 2u; p2->header.structSize = 6560u;
            ok (fc_master_params_defaults (p2) == FC_OK && p2->header.abiVersion == 2u && intact (pb2, 6560u),
                "a v2 parameter set: 6560 bytes written, and the 8 where v3 keeps its mix untouched");
        }
    }

    //==========================================================================
    // v4 — the gain-reduction trace. Its checks in the header's order, its capacity rule, and its OWNERSHIP: the trace
    // is the solution's, so neither a second solve on the same chain handle nor destroying that handle moves it.
    group ("v4: fc_solution_gr_trace — the header's order, the capacity rule, and a trace that outlives its chain");
    {
        fc_master h = make();
        fc_master_params p = goodParams();
        fc_loudness_request req {}; fc_loudness_request_default (&req);
        req.targetLufs = -4.0; req.maxTruePeakDbTp = -1.0; req.maxPasses = 3;    // loud enough that the limiter works
        const std::size_t frames = (std::size_t) (kFs * 3.0);
        auto in = tone (frames, kNch);
        std::vector<float> out (in.size(), 0.0f);
        fc_solution sol = 0;
        ok (fc_master_solve (h, &p, &req, in.data(), out.data(), (std::uint32_t) frames, &sol) == FC_OK, "PRECONDITION: a solution");

        std::vector<fc_gr_trace_bucket> first (1001u);
        std::uint32_t w = 777;
        ok (fc_solution_gr_trace (sol, FC_GR_STAGE_LIMITER, first.data(), 1001u, &w) == FC_OK && w == 1000u,
            "1000 buckets for 3 s, written into a buffer of 1001");
        double firstMax = 0.0;
        for (std::uint32_t i = 0; i < w; ++i) firstMax = std::max (firstMax, first[i].maxDb);
        ok (firstMax > 0.5, "PRECONDITION: the limiter worked on this solve (" + std::to_string (firstMax) + " dB peak)");
        std::vector<fc_gr_trace_bucket> small (10u);
        ok (fc_solution_gr_trace (sol, FC_GR_STAGE_LIMITER, small.data(), 10u, &w) == FC_OK && w == 10u
            && std::memcmp (small.data(), first.data(), 10u * sizeof (fc_gr_trace_bucket)) == 0,
            "a buffer of 10 is FILLED to its capacity with the first 10, not overrun");
        ok (fc_solution_gr_trace (sol, FC_GR_STAGE_LIMITER, nullptr, 0u, &w) == FC_OK && w == 0u,
            "asking for nothing is not an error, and `written` says 0");
        ok (fc_solution_gr_trace (sol, FC_GR_STAGE_LIMITER, nullptr, 4u, &w) == FC_ERR_NULL, "a null buffer with a capacity");
        auto* odd = reinterpret_cast<fc_gr_trace_bucket*> (static_cast<void*> ((char*) first.data() + 4));
        ok (fc_solution_gr_trace (sol, FC_GR_STAGE_LIMITER, odd, 4u, &w) == FC_ERR_ALIGNMENT, "a buffer off the 8-byte grid");
        ok (fc_solution_gr_trace (sol, FC_GR_STAGE_LIMITER, first.data(), 4u, nullptr) == FC_ERR_NULL, "a null `written`");
        {
            // `written` INSIDE the buckets: refused with SPAN before anything is written — not the count over bucket 0.
            std::vector<fc_gr_trace_bucket> alias (4u);
            alias[0].samples = 4242u; alias[1].nonFinite = 4343u;
            auto* inside  = &alias[0].samples;
            auto* inside2 = &alias[1].nonFinite;
            ok (fc_solution_gr_trace (sol, FC_GR_STAGE_LIMITER, alias.data(), 4u, inside) == FC_ERR_SPAN && alias[0].samples == 4242u
                && fc_solution_gr_trace (sol, FC_GR_STAGE_LIMITER, alias.data(), 4u, inside2) == FC_ERR_SPAN && alias[1].nonFinite == 4343u,
                "`written` pointing into the buckets is SPAN, and the refusal writes nothing");
            std::uint32_t w9 = 55u;
            ok (fc_solution_gr_trace (sol, 9, first.data(), 4u, &w9) == FC_ERR_ENUM && w9 == 55u,
                "and a refusal for the stage code leaves `written` as it was");
        }
        ok (fc_solution_gr_trace (sol, 2, first.data(), 4u, &w) == FC_ERR_ENUM && fc_solution_gr_trace (sol, -1, first.data(), 4u, &w) == FC_ERR_ENUM,
            "a stage code that names nothing is FC_ERR_ENUM");
        ok (fc_solution_gr_trace (sol, 7, nullptr, 0u, &w) == FC_ERR_ENUM, "... with a capacity of 0 too — never FC_OK");
        ok (fc_solution_gr_trace (sol, 7, nullptr, 4u, &w) == FC_ERR_NULL,
            "and the out-parameter comes before the field value: null buffer AND bad stage answer NULL (the header's order)");
        ok (fc_solution_gr_trace (0u, FC_GR_STAGE_LIMITER, first.data(), 4u, &w) == FC_ERR_HANDLE
            && fc_solution_gr_trace (h, FC_GR_STAGE_LIMITER, first.data(), 4u, &w) == FC_ERR_HANDLE,
            "a null handle and a CHAIN handle are not solutions");

        // THE LIVE LIMITER, bit for bit against the core. The same facade and C++ geometry as the group above, a request
        // loud enough that BOTH stages work — without it, a getter that zeroed the limiter's values would compare equal.
        {
            using namespace felitronics::mastering;
            fc_master hl = make();
            fc_master_params pl = goodParams();
            pl.limiter.ceilingDbTp = -1.3; pl.compressor.thresholdDb = -17.3;
            fc_loudness_request rl {}; fc_loudness_request_default (&rl);
            rl.targetLufs = -4.0; rl.maxTruePeakDbTp = -1.0; rl.maxPasses = 3;
            std::vector<float> viaAbi (in.size(), 0.0f), viaCpp (in.size(), 0.0f);
            fc_solution sl = 0;
            ok (fc_master_solve (hl, &pl, &rl, in.data(), viaAbi.data(), (std::uint32_t) frames, &sl) == FC_OK, "PRECONDITION: the loud ABI solve ran");
            MasteringChainConfig cc {};
            cc.internalBlock = 256; cc.monoBass = true; cc.clipper = true;
            MasteringChainParams cp {};
            cp.limiter.ceilingDbTp = -1.3; cp.compressor.thresholdDb = -17.3;
            MasteringChain chain; OfflineRenderer rend; TargetLoudnessSolver solver;
            const bool prep = rend.prepare (kNch, 4096) && chain.prepare (kFs, kNch, cc)
                           && solver.prepare (kFs, kNch, rend.blockSize(), chain.internalBlock(), chain.tapOversampleFactor());
            LoudnessRequest lr {}; lr.targetLufs = -4.0; lr.maxTruePeakDbTp = -1.0; lr.maxPasses = 3;
            const float* ip[2] { in.data(), in.data() + frames };
            float*       op[2] { viaCpp.data(), viaCpp.data() + frames };
            LoudnessSolution direct;
            if (prep) direct = solver.solve (chain, rend, cp, ip, op, kNch, (int) frames, lr);
            ok (prep && std::memcmp (viaAbi.data(), viaCpp.data(), viaAbi.size() * sizeof (float)) == 0,
                "PRECONDITION: the direct search delivers the same audio");
            for (const auto& [code, name, tr] : { std::tuple<std::int32_t, const char*, const GainReductionTrace*> { FC_GR_STAGE_LIMITER, "limiter", &direct.limiterTrace },
                                                  std::tuple<std::int32_t, const char*, const GainReductionTrace*> { FC_GR_STAGE_COMPRESSOR, "compressor", &direct.compressorTrace } })
            {
                std::vector<fc_gr_trace_bucket> tb (1000u);
                std::uint32_t wl = 0;
                double live = 0.0; std::size_t bad = 0;
                const bool got = fc_solution_gr_trace (sl, code, tb.data(), 1000u, &wl) == FC_OK && (int) wl == tr->buckets;
                for (std::uint32_t i = 0; got && i < wl; ++i)
                {
                    live = std::max (live, tr->bucket[i].maxDb);
                    if (std::memcmp (&tb[i].maxDb, &tr->bucket[i].maxDb, 8) != 0 || std::memcmp (&tb[i].meanDb, &tr->bucket[i].meanDb, 8) != 0
                        || tb[i].samples != tr->bucket[i].samples || tb[i].nonFinite != tr->bucket[i].nonFinite) ++bad;
                }
                ok (live > 0.5, std::string ("PRECONDITION: the ") + name + " trace is live (" + std::to_string (live) + " dB)");
                ok (got && bad == 0, std::string ("the live ") + name + " trace through the ABI is the core's, every bucket bit for bit ("
                    + std::to_string (bad) + " differ)");
            }
            (void) fc_solution_destroy (sl);
            (void) fc_master_destroy (hl);
        }

        // OWNERSHIP. A second solve on the same chain handle at another target, then the chain destroyed: the first
        // solution's trace is the same bytes.
        fc_master_params pc = goodParams();
        fc_master_resolved rr {}; FC_INIT (rr);
        ok (fc_master_configure (h, &pc, &rr) == FC_OK, "PRECONDITION: the handle configured again after its solve");
        req.targetLufs = -20.0;
        fc_solution sol2 = 0;
        ok (fc_master_solve (h, &p, &req, in.data(), out.data(), (std::uint32_t) frames, &sol2) == FC_OK, "PRECONDITION: a second solve");
        std::vector<fc_gr_trace_bucket> second (1000u), again (1000u);
        (void) fc_solution_gr_trace (sol2, FC_GR_STAGE_LIMITER, second.data(), 1000u, &w);
        ok (std::memcmp (second.data(), first.data(), 1000u * sizeof (fc_gr_trace_bucket)) != 0,
            "PRECONDITION: the second solve's trace differs from the first's");
        ok (fc_master_destroy (h) == FC_OK, "the chain handle destroyed");
        ok (fc_solution_gr_trace (sol, FC_GR_STAGE_LIMITER, again.data(), 1000u, &w) == FC_OK && w == 1000u
            && std::memcmp (again.data(), first.data(), 1000u * sizeof (fc_gr_trace_bucket)) == 0,
            "the FIRST solution's trace is unchanged by the second solve and by its chain's destruction");
        (void) fc_solution_destroy (sol2);
        (void) fc_solution_destroy (sol);
        ok (fc_solution_gr_trace (sol, FC_GR_STAGE_LIMITER, again.data(), 1000u, &w) == FC_ERR_HANDLE, "a destroyed solution is stale");

        // THE PRICE, PINNED: a solution record is 2368 B, two `GainReductionTrace`s and two `QuantileHistogram`s,
        // neither one's store included. (2336 before the percentile work: each of the measurement's two stage
        // summaries gained the quantile it was read at and the fraction it was read for, 16 B apiece.)
        {
            using felitronics::mastering::GainReductionTrace;
            fc_master hb = make();
            fc_need nd {}; FC_INIT (nd);
            using felitronics::dynamics::offline::QuantileHistogram;
            ok (fc_master_need (hb, FC_NEED_SOLVE, 48000u, &nd) == FC_OK
                && nd.facadeBytes == kSolutionRecordRest + 2u * sizeof (GainReductionTrace)
                                      + 3u * sizeof (QuantileHistogram)
                                      + sizeof (felitronics::mastering::ActiveGainReductionStats)
                && sizeof (GainReductionTrace) == (24u + sizeof (GainReductionTrace::bucket) + 7u) / 8u * 8u,
                "a solution record costs " + std::to_string (kSolutionRecordRest) + " B plus two traces of "
                + std::to_string (sizeof (GainReductionTrace)) + " B, three window histograms of "
                + std::to_string (sizeof (QuantileHistogram)) + " B and K11's gated summary of "
                + std::to_string (sizeof (felitronics::mastering::ActiveGainReductionStats))
                + " B, no store included (" + std::to_string (nd.facadeBytes) + ")");
            (void) fc_master_destroy (hb);
        }

        // A VERDICT BEFORE ANY RENDER carries no trace: zero buckets, not valid, and zero written.
        fc_master h2 = make();
        fc_loudness_request bad {}; fc_loudness_request_default (&bad);            // no target: InvalidRequest
        fc_solution sol3 = 0;
        ok (fc_master_solve (h2, &p, &bad, in.data(), out.data(), (std::uint32_t) frames, &sol3) == FC_OK, "PRECONDITION: a verdict");
        fc_measurement m3 {}; FC_INIT (m3);
        ok (fc_solution_measurement (sol3, &m3) == FC_OK && m3.limiterGrTraceBuckets == 0 && m3.compressorGrTraceBuckets == 0
            && m3.limiterGrTraceValid == 0 && m3.compressorGrTraceValid == 0, "InvalidRequest: no buckets, not valid");
        ok (fc_solution_gr_trace (sol3, FC_GR_STAGE_COMPRESSOR, first.data(), 1000u, &w) == FC_OK && w == 0u, "and nothing written");
        (void) fc_solution_destroy (sol3);
        (void) fc_master_destroy (h2);
    }

    //==========================================================================
    // v6
    group ("v6: the dual release, grTraceBuckets, fc_solution_gr_trace64 and fc_master_need_solve");
    {
        using namespace felitronics::mastering;
        auto canaried = [] (std::size_t bytes) { std::vector<unsigned char> b (bytes + 64u, 0xC3); return b; };
        auto intact   = [] (const std::vector<unsigned char>& b, std::size_t from)
        { for (std::size_t i = from; i < b.size(); ++i) if (b[i] != 0xC3) return false; return true; };

        // THE DUAL RELEASE: off by default, mapped when on, floored by the core and read back out of it, a NaN refused.
        {
            fc_master h = make();
            fc_master_resolved r {}; FC_INIT (r);
            fc_master_params p {}; FC_INIT (p);
            ok (fc_master_params_defaults (&p) == FC_OK && p.limiterDualRelease == 0 && p._pad0 == 0
                && p.limiterSlowReleaseMs == felitronics::limiter::TruePeakLimiterParams {}.slowReleaseMs,
                "the defaults writer: the dual release off, the core's slow release, the padding 0");
            ok (fc_master_configure (h, &p, &r) == FC_OK && r.limiterSlowReleaseMs == 0.0, "off, the resolved slow release is 0");

            // The same chain through C++, for the number the core reads back.
            MasteringChainConfig cc {}; cc.internalBlock = 256; cc.monoBass = true; cc.clipper = true;
            auto directSlow = [&] (double slowMs)
            {
                MasteringChain chain;
                MasteringChainParams cp {};
                cp.limiter.dualRelease = true; cp.limiter.slowReleaseMs = slowMs;
                chain.setParams (cp);
                return chain.prepare (kFs, kNch, cc) ? chain.resolved().limiterSlowReleaseMs : -1.0;
            };
            p.limiterDualRelease = 1; p.limiterSlowReleaseMs = 180.0;
            const double want180 = directSlow (180.0);
            // The effective release, read back from a float coefficient.
            ok (fc_master_configure (h, &p, &r) == FC_OK && r.limiterSlowReleaseMs == want180 && std::fabs (want180 - 180.0) < 0.5,
                "on, the slow release the core runs, read out of it (180 ms asked, " + std::to_string (r.limiterSlowReleaseMs) + " read)");
            p.limiterDualRelease = -7;
            ok (fc_master_configure (h, &p, &r) == FC_OK && r.limiterSlowReleaseMs == want180, "any non-zero flag is on");
            p.limiterDualRelease = 1; p.limiterSlowReleaseMs = 0.0;
            const double wantFloor = directSlow (0.0);
            ok (fc_master_configure (h, &p, &r) == FC_OK && r.limiterSlowReleaseMs == wantFloor && wantFloor > 0.16 && wantFloor < 0.17,
                "a slow release of 0 is floored by the core at 8 samples and read back (" + std::to_string (r.limiterSlowReleaseMs) + " ms)");
            p.limiterSlowReleaseMs = std::numeric_limits<double>::quiet_NaN();
            ok (fc_master_configure (h, &p, &r) == FC_ERR_NON_FINITE, "a NaN slow release: NON_FINITE");
            p.limiterDualRelease = 0; p.limiterSlowReleaseMs = std::numeric_limits<double>::infinity();
            ok (fc_master_configure (h, &p, &r) == FC_ERR_NON_FINITE, "and an infinite one with the dual release off, too");

            // A v5 parameter set: the two fields past its 6568 bytes are not read. The same bytes under v6 are.
            fc_master_params p5 = p;
            p5.limiterDualRelease = 1; p5.limiterSlowReleaseMs = 180.0;
            p5.header.abiVersion = 5u; p5.header.structSize = 6568u;
            ok (fc_master_configure (h, &p5, &r) == FC_OK && r.limiterSlowReleaseMs == 0.0,
                "a v5 parameter set: the dual release past its 6568 bytes is not read — the limiter runs its single release");
            FC_INIT (p5);
            ok (fc_master_configure (h, &p5, &r) == FC_OK && r.limiterSlowReleaseMs == want180, "and under v6 the same bytes switch it on");

            // A v5 resolved gets its 88 and not a byte more; a v6 one gets its 96.
            auto rb = canaried (96u);
            auto* r5 = reinterpret_cast<fc_master_resolved*> (rb.data());
            r5->header.abiVersion = 5u; r5->header.structSize = 88u;
            ok (fc_master_configure (h, &p5, r5) == FC_OK && r5->header.structSize == 88u && intact (rb, 88u),
                "configure: a v5 resolved, nothing written past its 88 where v6 keeps the slow release");
            auto rb6 = canaried (96u);
            auto* r6 = reinterpret_cast<fc_master_resolved*> (rb6.data());
            r6->header.abiVersion = 6u; r6->header.structSize = 96u;
            ok (fc_master_configure (h, &p5, r6) == FC_OK && r6->limiterSlowReleaseMs == want180 && intact (rb6, 96u),
                "and a v6 one gets its 96, the slow release included");
            (void) fc_master_destroy (h);
        }

        // 4096 buckets through both copiers against the core, and the 64-bit copier's checks.
        {
            const std::size_t frames = (std::size_t) (kFs * 3.0);
            auto in = tone (frames, kNch);
            fc_master h = make();
            fc_master_params p {}; FC_INIT (p); (void) fc_master_params_defaults (&p);
            p.limiter.ceilingDbTp = -1.3; p.compressor.thresholdDb = -17.3;
            fc_loudness_request req {}; FC_INIT (req); (void) fc_loudness_request_defaults (&req);
            ok (req.grTraceBuckets == 1000 && req._pad0 == 0, "the request's defaults writer: 1000 buckets, the padding 0");
            req.targetLufs = -4.0; req.maxTruePeakDbTp = -1.0; req.maxPasses = 3; req.grTraceBuckets = 4096;
            std::vector<float> viaAbi (in.size(), 0.0f), viaCpp (in.size(), 0.0f);
            fc_solution sol = 0;
            ok (fc_master_solve (h, &p, &req, in.data(), viaAbi.data(), (std::uint32_t) frames, &sol) == FC_OK, "PRECONDITION: a solution");

            MasteringChainConfig cc {}; cc.internalBlock = 256; cc.monoBass = true; cc.clipper = true;
            MasteringChainParams cp {}; cp.limiter.ceilingDbTp = -1.3; cp.compressor.thresholdDb = -17.3;
            MasteringChain chain; OfflineRenderer rend; TargetLoudnessSolver solver;
            const bool prep = rend.prepare (kNch, 4096) && chain.prepare (kFs, kNch, cc)
                           && solver.prepare (kFs, kNch, rend.blockSize(), chain.internalBlock(), chain.tapOversampleFactor());
            LoudnessRequest lr {}; lr.targetLufs = -4.0; lr.maxTruePeakDbTp = -1.0; lr.maxPasses = 3; lr.grTraceBuckets = 4096;
            const float* ip[2] { in.data(), in.data() + frames };
            float*       op[2] { viaCpp.data(), viaCpp.data() + frames };
            LoudnessSolution direct;
            if (prep) direct = solver.solve (chain, rend, cp, ip, op, kNch, (int) frames, lr);
            ok (prep && std::memcmp (viaAbi.data(), viaCpp.data(), viaAbi.size() * sizeof (float)) == 0,
                "PRECONDITION: the direct search delivers the same audio");
            fc_measurement m {}; FC_INIT (m);
            ok (fc_solution_measurement (sol, &m) == FC_OK && m.limiterGrTraceBuckets == 4096 && m.compressorGrTraceBuckets == 4096
                && m.limiterGrTraceValid == 1 && m.compressorGrTraceValid == 1, "the measurement: 4096 buckets per stage, both valid");
            for (const auto& [code, name, tr] : { std::tuple<std::int32_t, const char*, const GainReductionTrace*> { FC_GR_STAGE_LIMITER, "limiter", &direct.limiterTrace },
                                                  std::tuple<std::int32_t, const char*, const GainReductionTrace*> { FC_GR_STAGE_COMPRESSOR, "compressor", &direct.compressorTrace } })
            {
                std::vector<fc_gr_trace_bucket64> t64 (4097u);
                std::vector<fc_gr_trace_bucket>   t32 (4096u);
                std::uint32_t w64 = 0, w32 = 0;
                const bool got = fc_solution_gr_trace64 (sol, code, t64.data(), 4097u, &w64) == FC_OK && w64 == 4096u
                              && fc_solution_gr_trace (sol, code, t32.data(), 4096u, &w32) == FC_OK && w32 == 4096u
                              && tr->buckets == 4096 && tr->bucket.size() == 4096u;
                double live = 0.0; std::size_t bad = 0;
                for (std::uint32_t i = 0; got && i < 4096u; ++i)
                {
                    const GainReductionTraceBucket& b = tr->bucket[i];
                    live = std::max (live, b.maxDb);
                    if (std::memcmp (&t64[i].maxDb, &b.maxDb, 8) != 0 || std::memcmp (&t64[i].meanDb, &b.meanDb, 8) != 0
                        || t64[i].samples != b.samples || t64[i].nonFinite != b.nonFinite
                        || std::memcmp (&t32[i].maxDb, &b.maxDb, 8) != 0 || std::memcmp (&t32[i].meanDb, &b.meanDb, 8) != 0
                        || (std::uint64_t) t32[i].samples != b.samples || (std::uint64_t) t32[i].nonFinite != b.nonFinite) ++bad;
                }
                ok (live > 0.5, std::string ("PRECONDITION: the ") + name + " trace is live (" + std::to_string (live) + " dB)");
                ok (got && bad == 0, std::string ("the ") + name + " trace at 4096 buckets, through both copiers, is the core's bit for bit ("
                    + std::to_string (bad) + " differ)");
            }

            std::vector<fc_gr_trace_bucket64> b64 (16u);
            std::uint32_t w = 555u;
            ok (fc_solution_gr_trace64 (sol, FC_GR_STAGE_LIMITER, b64.data(), 10u, &w) == FC_OK && w == 10u
                && b64[0].samples == direct.limiterTrace.bucket[0].samples && b64[10].samples == 0u,
                "the capacity is in BUCKETS: 10 asked, 10 written, the 11th untouched");
            ok (fc_solution_gr_trace64 (sol, FC_GR_STAGE_LIMITER, nullptr, 0u, &w) == FC_OK && w == 0u, "asking for nothing: FC_OK, 0 written");
            w = 555u;
            ok (fc_solution_gr_trace64 (sol, FC_GR_STAGE_LIMITER, nullptr, 4u, &w) == FC_ERR_NULL && w == 555u, "a null buffer with a capacity: NULL");
            auto* odd = reinterpret_cast<fc_gr_trace_bucket64*> (static_cast<void*> ((char*) b64.data() + 4));
            ok (fc_solution_gr_trace64 (sol, FC_GR_STAGE_LIMITER, odd, 4u, &w) == FC_ERR_ALIGNMENT && w == 555u, "off the 8-byte grid: ALIGNMENT");
            ok (fc_solution_gr_trace64 (sol, FC_GR_STAGE_LIMITER, b64.data(), 4u, nullptr) == FC_ERR_NULL, "a null `written`: NULL");
            b64[1].samples = 4242u;
            auto* inside = reinterpret_cast<std::uint32_t*> (static_cast<void*> (&b64[1].samples));
            ok (fc_solution_gr_trace64 (sol, FC_GR_STAGE_LIMITER, b64.data(), 4u, inside) == FC_ERR_SPAN && b64[1].samples == 4242u,
                "`written` inside the buckets: SPAN, and nothing written");
            ok (fc_solution_gr_trace64 (sol, 2, b64.data(), 4u, &w) == FC_ERR_ENUM && fc_solution_gr_trace64 (sol, 7, nullptr, 0u, &w) == FC_ERR_ENUM
                && w == 555u, "a stage that names nothing: ENUM, with a capacity of 0 too, `written` untouched");
            ok (fc_solution_gr_trace64 (sol, 7, nullptr, 4u, &w) == FC_ERR_NULL, "the out-parameter before the field value: NULL");
            ok (fc_solution_gr_trace64 (h, FC_GR_STAGE_LIMITER, b64.data(), 4u, &w) == FC_ERR_HANDLE, "a chain handle is not a solution");
            (void) fc_solution_destroy (sol);

            // A v5 request: the count past its 120 bytes is not read.
            fc_master_params pc {}; FC_INIT (pc); (void) fc_master_params_defaults (&pc);
            fc_master_resolved rr {}; FC_INIT (rr);
            fc_loudness_request q5 = req; q5.header.abiVersion = 5u; q5.header.structSize = 120u;
            std::vector<float> out (in.size(), 0.0f);
            ok (fc_master_configure (h, &pc, &rr) == FC_OK && fc_master_solve (h, &p, &q5, in.data(), out.data(), (std::uint32_t) frames, &sol) == FC_OK
                && fc_solution_measurement (sol, &m) == FC_OK && m.limiterGrTraceBuckets == 1000,
                "a v5 request: `grTraceBuckets` past its 120 bytes is not read — 1000 buckets");
            (void) fc_solution_destroy (sol);

            // The count's range is the core's verdict, delivered with FC_OK.
            struct Row { std::int32_t buckets; bool admitted; int want; };
            for (const Row rw : { Row { 0, false, 0 }, Row { -1, false, 0 }, Row { std::numeric_limits<std::int32_t>::min(), false, 0 },
                                  Row { 65537, false, 0 }, Row { 1, true, 1 }, Row { 65536, true, 65536 } })
            {
                fc_loudness_request q = req; q.grTraceBuckets = rw.buckets; q.maxPasses = 1;
                fc_solution_summary su {}; FC_INIT (su);
                const bool solved = fc_master_configure (h, &pc, &rr) == FC_OK
                                 && fc_master_solve (h, &p, &q, in.data(), out.data(), (std::uint32_t) frames, &sol) == FC_OK
                                 && fc_solution_summary_get (sol, &su) == FC_OK && fc_solution_measurement (sol, &m) == FC_OK;
                ok (solved && (su.status == FC_SOLVE_INVALID_REQUEST) == ! rw.admitted && m.limiterGrTraceBuckets == rw.want,
                    "grTraceBuckets " + std::to_string (rw.buckets) + ": " + (rw.admitted ? "a render" : "InvalidRequest, no trace")
                    + " (" + std::to_string (m.limiterGrTraceBuckets) + " buckets)");
                (void) fc_solution_destroy (sol);
            }
            (void) fc_master_destroy (h);
        }

        // fc_master_need_solve: FC_NEED_SOLVE's at 1000, the traces' difference at 65536, 0 for a refused count, the checks,
        // and a 65 536-bucket solve's allocation.
        {
            fc_master h = make();
            const std::uint32_t n = 192000u;
            fc_loudness_request req {}; FC_INIT (req); (void) fc_loudness_request_defaults (&req);
            fc_need base {}, same {}, big {}, zero {}; FC_INIT (base); FC_INIT (same); FC_INIT (big); FC_INIT (zero);
            ok (fc_master_need (h, FC_NEED_SOLVE, n, &base) == FC_OK && fc_master_need_solve (h, &req, n, &same) == FC_OK
                && same.callBytes == base.callBytes && same.facadeBytes == base.facadeBytes && same.solverPrepareBytes == base.solverPrepareBytes
                && same.solverPrepared == base.solverPrepared && same._pad0 == 0,
                "at the default 1000 buckets it is FC_NEED_SOLVE's budget, field for field");
            req.grTraceBuckets = 65536;
            ok (fc_master_need_solve (h, &req, n, &big) == FC_OK && big.callBytes == base.callBytes + 2u * (65536u - 1000u) * 32u,
                "65 536 buckets: two traces of 65 536 x 32 B in place of 1000 x 32 B");
            req.grTraceBuckets = 0;
            ok (fc_master_need_solve (h, &req, n, &zero) == FC_OK && zero.callBytes == 0u, "a count the core refuses: FC_OK, callBytes 0");
            req.grTraceBuckets = 4096;
            fc_need s1 {}, d1 {}; FC_INIT (s1); FC_INIT (d1);
            ok (fc_master_need_solve (h, &req, 1000u, &s1) == FC_OK && fc_master_need (h, FC_NEED_SOLVE, 1000u, &d1) == FC_OK
                && s1.callBytes == d1.callBytes && s1.callBytes > 0u,
                "4096 buckets over 1000 frames hold one bucket per frame — the budget of 1000 buckets");
            fc_loudness_request q5 = req; q5.header.abiVersion = 5u; q5.header.structSize = 120u;
            fc_need n5 {}; FC_INIT (n5);
            ok (fc_master_need_solve (h, &q5, n, &n5) == FC_OK && n5.callBytes == base.callBytes,
                "a v5 request is budgeted at its 1000 buckets, whatever lies past its 120 bytes");
            fc_loudness_request bad = req; bad.header.abiVersion = 0u;
            ok (fc_master_need_solve (h, &bad, n, &big) == FC_ERR_ABI_VERSION, "a request of version 0: ABI_VERSION");
            ok (fc_master_need_solve (h, nullptr, n, &big) == FC_ERR_NULL && fc_master_need_solve (h, &req, n, nullptr) == FC_ERR_NULL,
                "a null request or out: NULL");
            fc_need badOut {}; FC_INIT (badOut); badOut.header.structSize = 32u;
            ok (fc_master_need_solve (h, &bad, n, &badOut) == FC_ERR_STRUCT_SIZE, "`out` before `req`: a bad out and a bad request answer STRUCT_SIZE");
            ok (fc_master_need_solve (h, &req, 0x80000000u, &big) == FC_ERR_RANGE, "frames past INT_MAX: RANGE");
            ok (fc_master_need_solve (0u, &req, n, &big) == FC_ERR_HANDLE, "a null handle: HANDLE");

            ok (fc_master_set_channel_weight (h, 0, 1.0) == FC_OK, "PRECONDITION: the solver prepared");
            req.grTraceBuckets = 65536;
            fc_need nb {}; FC_INIT (nb);
            ok (fc_master_need_solve (h, &req, n, &nb) == FC_OK, "PRECONDITION: budgeted");
            req.targetLufs = -14.0; req.maxTruePeakDbTp = -1.0;
            fc_master_params p {}; FC_INIT (p); (void) fc_master_params_defaults (&p);
            auto in = tone (n, kNch);
            std::vector<float> out (in.size(), 0.0f);
            fc_solution sol = 0;
            alloc::plainObjectSize.store ((std::size_t) nb.facadeBytes, std::memory_order_relaxed);
            const long long before = alloc::bytes.load();
            const fc_status sv = fc_master_solve (h, &p, &req, in.data(), out.data(), n, &sol);
            const long long got = alloc::bytes.load() - before;
            alloc::plainObjectSize.store (0, std::memory_order_relaxed);
            fc_solution_summary sum {}; FC_INIT (sum);
            fc_measurement m {}; FC_INIT (m);
            ok (sv == FC_OK && fc_solution_summary_get (sol, &sum) == FC_OK && sum.passes > 0 && fc_solution_measurement (sol, &m) == FC_OK
                && m.limiterGrTraceBuckets == 65536, "PRECONDITION: the search rendered, 65 536 buckets");
            const long long traces = 2LL * 65536LL * 32LL;
            const long long perPass = (long long) nb.callBytes - traces - (long long) kGrWindowBytes;
            ok (perPass > 0 && got == (long long) sum.passes * perPass + traces + (long long) kGrWindowBytes
                               + (long long) nb.facadeBytes,
                "a 65 536-bucket solve allocates passes x its meters + two traces of 2 097 152 B + its window histograms "
                "+ the record (" + std::to_string (got) + ")");
            (void) fc_solution_destroy (sol);
            (void) fc_master_destroy (h);
        }

        // A solve stopped on its first pass's record allocates exactly its budget, on a plain and on a delivering handle.
        {
            struct StopOnRecord { static std::int32_t fn (void*, const fc_progress* e) { return e->hasRecord != 0 ? 0 : 1; } };
            for (const double dr : { 0.0, 96000.0 })
            {
                fc_master_config c {}; FC_INIT (c);
                (void) fc_master_config_defaults (&c);
                c.sampleRate = kFs; c.channels = kNch; c.deliveryRate = dr;
                fc_master h = 0;
                const std::uint32_t n = 65536u;
                const std::uint32_t outN = dr == 0.0 ? n : 2u * n;
                fc_master_params p {}; FC_INIT (p); (void) fc_master_params_defaults (&p);
                fc_loudness_request req {}; FC_INIT (req); (void) fc_loudness_request_defaults (&req);
                req.targetLufs = -14.0; req.maxTruePeakDbTp = -1.0; req.grTraceBuckets = 65536;
                fc_need nb {}; FC_INIT (nb);
                const bool ready = fc_master_create (&c, &h) == FC_OK && fc_master_set_channel_weight (h, 0, 1.0) == FC_OK
                                && fc_master_set_progress (h, &StopOnRecord::fn, nullptr) == FC_OK
                                && fc_master_need_solve (h, &req, n, &nb) == FC_OK && nb.solverPrepared == 1;
                auto in = tone (n, kNch);
                std::vector<float> out ((std::size_t) outN * kNch, 0.0f);
                fc_solution sol = 0;
                alloc::plainObjectSize.store ((std::size_t) nb.facadeBytes, std::memory_order_relaxed);
                const long long before = alloc::bytes.load();
                const fc_status sv = ! ready ? FC_ERR_STATE
                                   : dr == 0.0 ? fc_master_solve (h, &p, &req, in.data(), out.data(), n, &sol)
                                               : fc_master_solve_delivered (h, &p, &req, in.data(), n, out.data(), outN, &sol);
                const long long got = alloc::bytes.load() - before;
                alloc::plainObjectSize.store (0, std::memory_order_relaxed);
                ok (ready && sv == FC_ERR_CANCELLED && got == (long long) (nb.callBytes + nb.facadeBytes),
                    std::string (dr == 0.0 ? "plain" : "48 -> 96 kHz") + ": stopped on the first record, " + std::to_string (got)
                    + " B allocated against a budget of " + std::to_string (nb.callBytes + nb.facadeBytes) + " B");
                (void) fc_master_destroy (h);
            }
        }
    }

    //==========================================================================
    // v3 — `compressorMix` crosses like every other value: finite is this file's check, the range is the core's.
    group ("v3: compressorMix — mapped, clamped by the core and read back, a NaN refused");
    {
        fc_master h = make();
        fc_master_resolved r {}; FC_INIT (r);
        fc_master_params p = goodParams(); FC_INIT (p);
        (void) fc_master_params_defaults (&p);
        p.compressorMix = 0.37;
        ok (fc_master_configure (h, &p, &r) == FC_OK && r.compressorMix == (double) (float) 0.37,
            "0.37 is applied, and the read-back is the value the stage holds (a float, as the core keeps it)");
        p.compressorMix = 2.0;
        ok (fc_master_configure (h, &p, &r) == FC_OK && r.compressorMix == 1.0, "2 is CLAMPED by the core to 1, and reported");
        p.compressorMix = -1.0;
        ok (fc_master_configure (h, &p, &r) == FC_OK && r.compressorMix == 0.0, "and -1 to 0");
        fc_master_resolved before = r;
        p.compressorMix = std::numeric_limits<double>::quiet_NaN();
        ok (fc_master_configure (h, &p, &r) == FC_ERR_NON_FINITE, "a NaN is refused — the core would have made it 1 silently");
        ok (std::memcmp (&before, &r, sizeof r) == 0, "and the refused configure moved nothing");
        (void) fc_master_destroy (h);

        fc_master_config nc = goodConfig(); nc.compressor = 0;
        fc_master hn = 0;
        p.compressorMix = 0.37;
        ok (fc_master_create (&nc, &hn) == FC_OK && fc_master_configure (hn, &p, &r) == FC_OK && r.compressorMix == 0.0,
            "with no compressor in the topology there is nothing to mix: the read-back is 0");
        (void) fc_master_destroy (hn);
    }

    //==========================================================================
    // P57b — THE DELIVERING HANDLE: states, lengths, spans and budgets. The bit-exactness of its renders against the
    // core lives in `fcore_master selftest`; the budgets are here, held against the counter.
    group ("v2: the delivering handle — refusals, the delivered length, and budgets to the byte");
    {
        auto deliveringConfig = [] (double fs, double dr)
        {
            fc_master_config c {}; FC_INIT (c);
            (void) fc_master_config_defaults (&c);
            c.sampleRate = fs; c.channels = kNch; c.deliveryRate = dr;
            return c;
        };

        // CREATE: the delivering budget is the create's allocation, and a pair the resampler cannot plan is refused
        // by both calls with nothing allocated.
        {
            fc_master_config c = deliveringConfig (48000.0, 96000.0);
            fc_need nd {}; FC_INIT (nd);
            ok (fc_master_need_create (&c, &nd) == FC_OK && nd.callBytes > 0u, "a delivering create is budgeted");
            fc_master_config legacy = c; legacy.deliveryRate = 0.0;
            fc_need nl {}; FC_INIT (nl);
            ok (fc_master_need_create (&legacy, &nl) == FC_OK && nd.callBytes > nl.callBytes,
                "and costs more than the same geometry without the converter");
            alloc::plainObjectSize.store ((std::size_t) nd.facadeBytes, std::memory_order_relaxed);
            fc_master h = 0;
            const long long before = alloc::bytes.load();
            const fc_status cs = fc_master_create (&c, &h);
            const long long got = alloc::bytes.load() - before;
            const long long budget = (long long) (nd.callBytes + nd.facadeBytes);
            ok (cs == FC_OK && got == budget, "the delivering create allocates its budget plus the facade's record ("
                + std::to_string (got) + " against " + std::to_string (budget) + ")");
            (void) fc_master_destroy (h);
            alloc::plainObjectSize.store (0, std::memory_order_relaxed);

            struct Refusal { const char* what; double dr; fc_status want; };
            const Refusal rs[] = {
                { "a non-integer delivery rate", 44100.5, FC_ERR_REFUSED_BY_CORE },
                { "a negative delivery rate",   -48000.0, FC_ERR_REFUSED_BY_CORE },
                { "a NaN delivery rate", std::numeric_limits<double>::quiet_NaN(), FC_ERR_NON_FINITE },
                { "an infinite delivery rate", std::numeric_limits<double>::infinity(), FC_ERR_NON_FINITE },
            };
            for (const Refusal& r : rs)
            {
                fc_master_config bc = deliveringConfig (48000.0, r.dr);
                fc_need bn {}; FC_INIT (bn);
                fc_master bh = 0;
                const fc_status ns = fc_master_need_create (&bc, &bn);
                const long long b0 = alloc::bytes.load();
                const fc_status bs = fc_master_create (&bc, &bh);
                const long long leaked = alloc::bytes.load() - b0;
                ok (ns == r.want && bs == r.want && leaked == 0, std::string (r.what) + ": refused by both, nothing allocated");
                if (bs == FC_OK) (void) fc_master_destroy (bh);
            }
        }

        // STATE: each kind of handle refuses the other kind's calls.
        {
            fc_master_config c = deliveringConfig (48000.0, 44100.0);
            fc_master dh = 0, lh = make();
            ok (fc_master_create (&c, &dh) == FC_OK && lh != 0, "PRECONDITION: one delivering handle, one plain");
            fc_master_params p = goodParams();
            fc_loudness_request req {}; fc_loudness_request_default (&req);
            req.targetLufs = -14.0; req.maxTruePeakDbTp = -1.0;
            auto in = tone (4800, kNch);
            std::vector<float> out (in.size(), 0.0f);
            std::uint32_t written = 0, d = 0;
            fc_solution sol = 0;
            ok (fc_master_process (dh, in.data(), out.data(), 4800u) == FC_ERR_STATE
                && fc_master_process (dh, in.data(), out.data(), 0u) == FC_ERR_STATE, "process on a delivering handle: STATE, even for 0 frames");
            ok (fc_master_flush (dh, out.data(), 4800u, &written) == FC_ERR_STATE, "flush: STATE");
            ok (fc_master_solve (dh, &p, &req, in.data(), out.data(), 4800u, &sol) == FC_ERR_STATE, "solve: STATE");
            ok (fc_master_delivered_frames (lh, 4800u, &d) == FC_ERR_STATE, "delivered_frames on a plain handle: STATE");
            ok (fc_master_render_delivered (lh, in.data(), 4800u, out.data(), 4800u) == FC_ERR_STATE, "render_delivered: STATE");
            ok (fc_master_solve_delivered (lh, &p, &req, in.data(), 4800u, out.data(), 4800u, &sol) == FC_ERR_STATE,
                "solve_delivered: STATE");
            (void) fc_master_destroy (dh);
            (void) fc_master_destroy (lh);
        }

        // THE DELIVERED LENGTH is the converter's, on every pair — forwarded, not recomputed — with the traps the
        // converter's own suite names, and the narrowing on both sides of INT_MAX.
        {
            constexpr double kRates[] = { 44100.0, 48000.0, 88200.0, 96000.0, 176400.0, 192000.0 };
            constexpr std::uint32_t kLens[] = { 0u, 1u, 147u, 1000u, 176401u, 10000019u };
            int pairs = 0, off = 0;
            for (double a : kRates)
                for (double b : kRates)
                {
                    fc_master_config c = deliveringConfig (a, b);
                    fc_master h = 0;
                    if (fc_master_create (&c, &h) != FC_OK) { ++off; continue; }
                    ++pairs;
                    for (std::uint32_t n : kLens)
                    {
                        std::uint32_t d = 0xFFFFFFFFu;
                        const fc_status st = fc_master_delivered_frames (h, n, &d);
                        if (st != FC_OK || (long long) d != felitronics::mastering::DeliveryConverter::deliveredFrames (a, b, n)) ++off;
                    }
                    (void) fc_master_destroy (h);
                }
            ok (pairs == 36 && off == 0, "every pair of the six rates (equal rates included) answers the converter's length ("
                + std::to_string (pairs) + " pairs, " + std::to_string (off) + " off)");

            fc_master_config c48 = deliveringConfig (44100.0, 48000.0);
            fc_master h = 0;
            std::uint32_t d = 0;
            ok (fc_master_create (&c48, &h) == FC_OK && fc_master_delivered_frames (h, 147u, &d) == FC_OK && d == 160u,
                "147 frames at 44.1 -> 48 kHz are 160 — the double formula says 161");
            (void) fc_master_destroy (h);
            fc_master_config c441 = deliveringConfig (176400.0, 44100.0);
            ok (fc_master_create (&c441, &h) == FC_OK && fc_master_delivered_frames (h, 176401u, &d) == FC_OK && d == 44101u,
                "176,401 frames at 176.4 -> 44.1 kHz are 44,101");
            (void) fc_master_destroy (h);

            fc_master_config up = deliveringConfig (44100.0, 192000.0);
            ok (fc_master_create (&up, &h) == FC_OK, "PRECONDITION: a 44.1 -> 192 kHz handle");
            std::uint32_t untouched = 0xABCDu;
            ok (fc_master_delivered_frames (h, 0x7FFFFFFFu, &untouched) == FC_ERR_RANGE && untouched == 0xABCDu,
                "INT_MAX input frames deliver past INT_MAX: FC_ERR_RANGE, and the out-parameter is untouched");
            ok (fc_master_delivered_frames (h, 0x80000000u, &untouched) == FC_ERR_RANGE, "an input count past INT_MAX: FC_ERR_RANGE");
            const std::uint32_t maxIn = (std::uint32_t) ((0x7FFFFFFFLL * 147LL) / 640LL);   // 44.1 -> 192 is 147:640 up
            ok (fc_master_delivered_frames (h, maxIn, &d) == FC_OK && d <= 0x7FFFFFFFu
                && fc_master_delivered_frames (h, maxIn + 1u, &d) == FC_ERR_RANGE,
                "and the edge is exactly where the delivered length stops fitting an int");
            ok (fc_master_delivered_frames (h, 1u, nullptr) == FC_ERR_NULL, "a null out-parameter: FC_ERR_NULL");
            fc_need big {}; FC_INIT (big);
            ok (fc_master_need (h, FC_NEED_SOLVE, 0x7FFFFFFFu, &big) == FC_ERR_RANGE,
                "and a budget for a delivered length the call would refuse is refused the same way");
            // THE COMPUTED LENGTH IS NARROWED BEFORE IT IS COMPARED. Both counts fit an int here; the delivered one
            // does not, so the answer is RANGE — CAPACITY would send the caller to resize a buffer that cannot exist.
            // (The mutation stand: swapping the two checks in the shared helper left every other check green.)
            {
                float dummy[4] {};
                fc_master_params pp = goodParams();
                fc_loudness_request rq {}; fc_loudness_request_default (&rq);
                rq.targetLufs = -14.0; rq.maxTruePeakDbTp = -1.0;
                fc_solution so = 0;
                ok (fc_master_render_delivered (h, dummy, 0x7FFFFFFFu, dummy + 2, 1000u) == FC_ERR_RANGE,
                    "render_delivered: a delivered length past INT_MAX is RANGE, not CAPACITY, whatever outFrames says");
                ok (fc_master_solve_delivered (h, &pp, &rq, dummy, 0x7FFFFFFFu, dummy + 2, 1000u, &so) == FC_ERR_RANGE,
                    "solve_delivered: the same");
            }
            (void) fc_master_destroy (h);
        }

        // RENDER: the order of the two lengths, the spans, and an audio path that asks the heap for nothing.
        {
            fc_master_config c = deliveringConfig (48000.0, 96000.0);
            fc_master h = 0;
            ok (fc_master_create (&c, &h) == FC_OK, "PRECONDITION: a 48 -> 96 kHz handle");
            fc_master_params p = goodParams();
            fc_master_resolved r {}; FC_INIT (r);
            ok (fc_master_configure (h, &p, &r) == FC_OK, "PRECONDITION: configured");
            const std::uint32_t n = 4800u, dn = 9600u;
            auto in = tone (n, kNch);
            std::vector<float> out ((std::size_t) dn * kNch, 0.0f);
            ok (fc_master_render_delivered (h, in.data(), n, out.data(), 0x80000000u) == FC_ERR_RANGE,
                "an outFrames past INT_MAX is RANGE — narrowing comes before the CAPACITY it would also fail");
            ok (fc_master_render_delivered (h, in.data(), 0x80000000u, out.data(), dn) == FC_ERR_RANGE,
                "and so is an inFrames past INT_MAX");
            ok (fc_master_render_delivered (h, in.data(), n, out.data(), dn - 1u) == FC_ERR_CAPACITY
                && fc_master_render_delivered (h, in.data(), n, out.data(), dn + 1u) == FC_ERR_CAPACITY,
                "a delivered length one short or one long: CAPACITY — the stride IS the length");
            ok (fc_master_render_delivered (h, nullptr, n, out.data(), dn) == FC_ERR_NULL
                && fc_master_render_delivered (h, in.data(), n, nullptr, dn) == FC_ERR_NULL, "null spans: NULL");
            std::vector<float> both ((std::size_t) dn * kNch + (std::size_t) n * kNch, 0.0f);
            ok (fc_master_render_delivered (h, both.data(), n, both.data() + 4, dn) == FC_ERR_SPAN,
                "an output that overlaps the input: SPAN");
            ok (fc_master_render_delivered (h, both.data(), n, both.data(), dn) == FC_ERR_SPAN,
                "and so does the SAME buffer, which a conversion cannot share");
            const long long allocs = alloc::count.load();
            const fc_status rs = fc_master_render_delivered (h, in.data(), n, out.data(), dn);
            const long long asked = alloc::count.load() - allocs;
            ok (rs == FC_OK && asked == 0, "the delivered render asks the heap for nothing (" + std::to_string (asked) + " requests)");
            double peak = 0.0;
            for (float v : out) peak = std::max (peak, (double) std::fabs (v));
            ok (peak > 0.05, "PRECONDITION: and it rendered something");

            // After a search the chain holds the SOLVER's parameters: refused until a configure puts a set back.
            fc_loudness_request req {}; fc_loudness_request_default (&req);
            req.targetLufs = -14.0; req.maxTruePeakDbTp = -1.0;
            auto longIn = tone (48000, kNch);
            std::vector<float> sOut ((std::size_t) 96000 * kNch, 0.0f);
            fc_solution sol = 0;
            ok (fc_master_solve_delivered (h, &p, &req, longIn.data(), 48000u, sOut.data(), 96000u, &sol) == FC_OK,
                "PRECONDITION: a delivered solve ran");
            ok (fc_master_render_delivered (h, in.data(), n, out.data(), dn) == FC_ERR_STATE, "then render_delivered: STATE");
            ok (fc_master_configure (h, &p, &r) == FC_OK && fc_master_render_delivered (h, in.data(), n, out.data(), dn) == FC_OK,
                "and a configure lifts it");
            ok (fc_master_solve_delivered (h, &p, &req, longIn.data(), 48000u, longIn.data(), 96000u, &sol) == FC_ERR_SPAN,
                "a delivered search over its own input: SPAN");
            ok (fc_master_solve_delivered (h, &p, &req, longIn.data(), 48000u, sOut.data(), 96000u,
                                           (fc_solution*) (void*) sOut.data()) == FC_ERR_SPAN,
                "and an out-handle inside the output: SPAN");
            (void) fc_solution_destroy (sol);
            (void) fc_master_destroy (h);
        }

        // BUDGETS: a delivered search and a delivered range, against the counter. The literal terms are the
        // converted programme (2 ch x delivered frames x 4 B); since P62 there is no drain buffer to add.
        {
            fc_master_config c = deliveringConfig (48000.0, 96000.0);
            fc_master h = 0;
            ok (fc_master_create (&c, &h) == FC_OK, "PRECONDITION: a 48 -> 96 kHz handle for the budgets");
            ok (fc_master_set_channel_weight (h, 0, 1.0) == FC_OK, "PRECONDITION: the solver prepared");
            const std::uint32_t n = 192000u;                     // 4 s in, 384 000 delivered
            const long long programme = 2LL * 384000LL * 4LL;
            fc_need solve {}; FC_INIT (solve);
            ok (fc_master_need (h, FC_NEED_SOLVE, n, &solve) == FC_OK && solve.solverPrepared == 1,
                "the delivered solve is budgeted, for the INPUT count");
            fc_master_params p = goodParams();
            fc_loudness_request req {}; fc_loudness_request_default (&req);
            req.targetLufs = -14.0; req.maxTruePeakDbTp = -1.0;
            auto in = tone (n, kNch);
            std::vector<float> out ((std::size_t) 384000 * kNch, 0.0f);
            fc_solution sol = 0;
            alloc::plainObjectSize.store ((std::size_t) solve.facadeBytes, std::memory_order_relaxed);   // the record: see above
            long long before = alloc::bytes.load();
            const fc_status sv = fc_master_solve_delivered (h, &p, &req, in.data(), n, out.data(), 384000u, &sol);
            const long long solveBytes = alloc::bytes.load() - before;
            alloc::plainObjectSize.store (0, std::memory_order_relaxed);
            fc_solution_summary sum {}; FC_INIT (sum);
            ok (sv == FC_OK && fc_solution_summary_get (sol, &sum) == FC_OK && sum.passes > 0, "PRECONDITION: the delivered search rendered");
            const long long traces = 2LL * 1000LL * 32LL;
            const long long perPass = (long long) solve.callBytes - programme - traces - (long long) kGrWindowBytes;
            ok (perPass > 0 && solveBytes == (long long) sum.passes * perPass + traces + (long long) kGrWindowBytes
                               + programme + (long long) solve.facadeBytes,
                "a delivered solve allocates passes x (both meters at 96 kHz) + two traces + its window histograms + the "
                "converted programme + the record (" + std::to_string (solveBytes) + ")");
            (void) fc_solution_destroy (sol);
            (void) fc_master_destroy (h);

            // EQUAL RATES read the input in place: no programme term at all.
            fc_master_config same = deliveringConfig (48000.0, 48000.0);
            fc_master hs = 0, hl = make();
            ok (fc_master_create (&same, &hs) == FC_OK, "PRECONDITION: an equal-rate delivering handle");
            fc_need ns {}, nl {}; FC_INIT (ns); FC_INIT (nl);
            ok (fc_master_need (hs, FC_NEED_SOLVE, n, &ns) == FC_OK && fc_master_need (hl, FC_NEED_SOLVE, n, &nl) == FC_OK
                && ns.callBytes == nl.callBytes, "at equal rates the delivered solve costs what a plain one does");
            (void) fc_master_destroy (hs);
            (void) fc_master_destroy (hl);

            // THE RANGE, JUDGED ON THE DELIVERED LENGTH. 132 300 frames at 44.1 kHz are 288 000 at 96 kHz — exactly 3 s,
            // measurable — while the same count read as 96 kHz frames is 1.38 s and would have budgeted nothing.
            fc_master_config lc = deliveringConfig (44100.0, 96000.0);
            ok (fc_master_create (&lc, &h) == FC_OK && fc_master_set_channel_weight (h, 0, 1.0) == FC_OK,
                "PRECONDITION: a 44.1 -> 96 kHz handle, solver prepared");
            ok (felitronics::mastering::TargetLoudnessSolver::measureRangeBytes (96000.0, 132300) == 0u,
                "PRECONDITION: the input count read at the delivery rate budgets 0 — the defect this op must not have");
            auto lin = tone (132300u, kNch);
            fc_need lra {}; FC_INIT (lra);
            ok (fc_master_need (h, FC_NEED_MEASURE_LRA, 132300u, &lra) == FC_OK && lra.callBytes > (std::uint64_t) (2 * 288000 * 4),
                "the delivered range is budgeted: the converted programme plus a meter");
            double v = 0.0;
            before = alloc::bytes.load();
            const fc_status ls = fc_master_measure_lra (h, lin.data(), 132300u, &v);
            const long long lraBytes = alloc::bytes.load() - before;
            ok (ls == FC_OK && lraBytes == (long long) lra.callBytes,
                "and allocates exactly that (" + std::to_string (lraBytes) + " against " + std::to_string (lra.callBytes) + ")");
            fc_need shortLra {}; FC_INIT (shortLra);
            ok (fc_master_need (h, FC_NEED_MEASURE_LRA, 132299u, &shortLra) == FC_OK && shortLra.callBytes == 0u,
                "one input frame fewer delivers under 3 s: budgeted 0");
            before = alloc::bytes.load();
            const fc_status ss = fc_master_measure_lra (h, lin.data(), 132299u, &v);
            const long long shortBytes = alloc::bytes.load() - before;
            ok (ss == FC_ERR_REFUSED_BY_CORE && shortBytes == 0, "and refused having converted — and allocated — nothing");
            (void) fc_master_destroy (h);
        }
    }

    //==========================================================================
    // P57b — THE DIVERSE-TESTING ROUND'S GAPS on the delivering handle, each one a change to the facade that left the
    // suite green before it was here.
    group ("v2/v3: what the diverse-testing round found untested");
    {
        auto deliveringConfig = [] (double fs, double dr)
        {
            fc_master_config c {}; FC_INIT (c);
            (void) fc_master_config_defaults (&c);
            c.sampleRate = fs; c.channels = kNch; c.deliveryRate = dr;
            return c;
        };
        fc_loudness_request req {}; FC_INIT (req);
        (void) fc_loudness_request_defaults (&req);
        req.targetLufs = -14.0; req.maxTruePeakDbTp = -1.0;

        // THE LENGTH EDGE, for every call that takes one: at 44.1 -> 192 kHz, 493 250 150 input frames deliver
        // exactly INT_MAX and 493 250 151 one more. The call's own arguments fit an int both times.
        {
            fc_master_config c = deliveringConfig (44100.0, 192000.0);
            fc_master h = 0;
            ok (fc_master_create (&c, &h) == FC_OK, "PRECONDITION: a 44.1 -> 192 kHz handle");
            std::uint32_t d = 0;
            ok (fc_master_delivered_frames (h, 493250150u, &d) == FC_OK && d == 0x7FFFFFFFu,
                "493 250 150 frames deliver exactly INT_MAX");
            fc_need atEdge {}; FC_INIT (atEdge);
            ok (fc_master_need (h, FC_NEED_SOLVE, 493250150u, &atEdge) == FC_OK
                && atEdge.callBytes >= 2ull * 0x7FFFFFFFull * 4ull,
                "and its budget is answered, with the converted programme in it — not 0 at the edge");
            float dummy[4] {};
            fc_master_params pp = goodParams();
            fc_solution so = 0;
            double lra = 0.0;
            ok (fc_master_render_delivered (h, dummy, 493250150u, dummy + 2, 1000u) == FC_ERR_CAPACITY,
                "render_delivered at the edge: the length is narrowable, so the answer is CAPACITY, not RANGE");
            ok (fc_master_solve_delivered (h, &pp, &req, dummy, 493250150u, dummy + 2, 1000u, &so) == FC_ERR_CAPACITY,
                "solve_delivered at the edge: CAPACITY");
            fc_need past {}; FC_INIT (past);
            ok (fc_master_delivered_frames (h, 493250151u, &d) == FC_ERR_RANGE
                && fc_master_need (h, FC_NEED_SOLVE, 493250151u, &past) == FC_ERR_RANGE
                && fc_master_render_delivered (h, dummy, 493250151u, dummy + 2, 0x7FFFFFFFu) == FC_ERR_RANGE
                && fc_master_solve_delivered (h, &pp, &req, dummy, 493250151u, dummy + 2, 0x7FFFFFFFu, &so) == FC_ERR_RANGE
                && fc_master_measure_lra (h, dummy, 493250151u, &lra) == FC_ERR_RANGE,
                "one frame past the edge: RANGE from every call, before a span is looked at");
            (void) fc_master_destroy (h);

            fc_master_config same = deliveringConfig (48000.0, 48000.0);
            ok (fc_master_create (&same, &h) == FC_OK, "PRECONDITION: an equal-rate handle");
            ok (fc_master_render_delivered (h, dummy, 0x7FFFFFFFu, dummy + 2, 0u) == FC_ERR_CAPACITY,
                "equal rates, INT_MAX frames: a count the core takes — CAPACITY for the wrong outFrames, not RANGE");
            (void) fc_master_destroy (h);
        }

        // THE SOLVE'S ALIASES, on the delivering entry point itself — it is a separate implementation of the rule.
        {
            fc_master_config c = deliveringConfig (48000.0, 96000.0);
            fc_master h = 0;
            ok (fc_master_create (&c, &h) == FC_OK, "PRECONDITION: a 48 -> 96 kHz handle");
            auto in = tone (4800, kNch);
            std::vector<float> out ((std::size_t) 9600 * kNch, 0.0f);
            std::vector<std::uint64_t> buf (6568u / 8u + 4u, 0u);
            fc_master_params base = goodParams();
            std::memcpy (buf.data(), &base, 6560u);
            auto* pv = reinterpret_cast<fc_master_params*> (static_cast<void*> (buf.data()));
            pv->header.abiVersion = 2u; pv->header.structSize = 6560u;
            auto* after = reinterpret_cast<fc_solution*> (static_cast<void*> ((char*) buf.data() + 6560));
            const fc_status s2 = fc_master_solve_delivered (h, pv, &req, in.data(), 4800u, out.data(), 9600u, after);
            ok (s2 == FC_OK, "an out-handle right after a v2 parameter set is outside it — FC_OK");
            if (s2 == FC_OK) (void) fc_solution_destroy (*after);
            (void) fc_master_destroy (h);
            ok (fc_master_create (&c, &h) == FC_OK, "PRECONDITION: a fresh one");
            pv->header.abiVersion = 3u; pv->header.structSize = 6568u;
            ok (fc_master_solve_delivered (h, pv, &req, in.data(), 4800u, out.data(), 9600u, after) == FC_ERR_SPAN,
                "and inside a v3 one: SPAN");
            fc_master_params p = goodParams();
            fc_loudness_request rq = req;
            ok (fc_master_solve_delivered (h, &p, &rq, in.data(), 4800u, out.data(), 9600u,
                                           (fc_solution*) (void*) in.data()) == FC_ERR_SPAN,
                "an out-handle inside the INPUT: SPAN");
            ok (fc_master_solve_delivered (h, &p, &rq, in.data(), 4800u, out.data(), 9600u,
                                           (fc_solution*) (void*) &p.eqBands[3].on) == FC_ERR_SPAN,
                "inside the parameters: SPAN");
            ok (fc_master_solve_delivered (h, &p, &rq, in.data(), 4800u, out.data(), 9600u,
                                           (fc_solution*) (void*) &rq.maxPasses) == FC_ERR_SPAN,
                "inside the request: SPAN");
            std::vector<float> both ((std::size_t) 9600 * kNch + (std::size_t) 4800 * kNch, 0.0f);
            fc_solution so = 0;
            ok (fc_master_solve_delivered (h, &p, &rq, both.data(), 4800u, both.data() + 4, 9600u, &so) == FC_ERR_SPAN,
                "an output that partly overlaps the input: SPAN");
            fc_loudness_request bad = rq; bad.header.abiVersion = FC_MASTER_ABI_VERSION + 1u;
            ok (fc_master_solve_delivered (h, &p, &bad, in.data(), 4800u, out.data(), 9600u, &so) == FC_ERR_ABI_VERSION,
                "a request from a newer ABI: refused");
            bad = rq; bad.header.structSize = 112u;
            ok (fc_master_solve_delivered (h, &p, &bad, in.data(), 4800u, out.data(), 9600u, &so) == FC_ERR_STRUCT_SIZE,
                "a request of the wrong size: refused");

            // A SEARCH THE SOLVER REFUSED MOVED NOTHING: no target is InvalidRequest, and render_delivered stays open.
            fc_loudness_request none {}; FC_INIT (none);
            (void) fc_loudness_request_defaults (&none);
            fc_master_resolved r {}; FC_INIT (r);
            ok (fc_master_configure (h, &p, &r) == FC_OK, "PRECONDITION: configured");
            fc_solution sv = 0;
            fc_solution_summary sum {}; FC_INIT (sum);
            ok (fc_master_solve_delivered (h, &p, &none, in.data(), 4800u, out.data(), 9600u, &sv) == FC_OK
                && fc_solution_summary_get (sv, &sum) == FC_OK && sum.status == FC_SOLVE_INVALID_REQUEST,
                "a delivered solve with no target: FC_OK carrying InvalidRequest");
            ok (fc_master_render_delivered (h, in.data(), 4800u, out.data(), 9600u) == FC_OK,
                "and the handle still renders: a refused search did not mark the chain as the solver's");
            (void) fc_solution_destroy (sv);

            // THE RE-PREPARATION of a delivering handle costs nothing, as a plain one's does — at the delivery rate.
            fc_need cn {}; FC_INIT (cn);
            const long long b0 = alloc::bytes.load();
            const fc_status cs = fc_master_configure (h, &p, &r);
            const long long got = alloc::bytes.load() - b0;
            ok (fc_master_need (h, FC_NEED_CONFIGURE, 0u, &cn) == FC_OK && cn.callBytes == 0u && cs == FC_OK && got == 0,
                "a delivering configure is budgeted 0 and allocates 0");
            (void) fc_master_destroy (h);
        }

        // THE AUDIO OUTPUT OVER THE CALLER'S STRUCTS, the equal-rate INT_MAX edge, a measure / solve / measure
        // sequence, and the version matrix on every OUT entry point that takes a resolved.
        {
            fc_master_config c = deliveringConfig (48000.0, 48000.0);
            fc_master h = 0;
            ok (fc_master_create (&c, &h) == FC_OK, "PRECONDITION: an equal-rate handle");
            std::uint32_t d = 0;
            ok (fc_master_delivered_frames (h, 0x7FFFFFFFu, &d) == FC_OK && d == 0x7FFFFFFFu,
                "equal rates: INT_MAX frames deliver INT_MAX — the edge is inclusive");
            // A buffer holding the parameters and, right after them, the start of the "output".
            std::vector<std::uint64_t> buf (6568u / 8u + 4096u, 0u);
            fc_master_params base = goodParams(); FC_INIT (base);
            (void) fc_master_params_defaults (&base);
            std::memcpy (buf.data(), &base, sizeof base);
            auto* pv = reinterpret_cast<fc_master_params*> (static_cast<void*> (buf.data()));
            auto in = tone (4800, kNch);
            auto* outOverParams = reinterpret_cast<float*> (static_cast<void*> ((char*) buf.data() + 6560));
            fc_solution so = 0;
            ok (fc_master_solve_delivered (h, pv, &req, in.data(), 4800u, outOverParams, 4800u, &so) == FC_ERR_SPAN,
                "an output that starts inside a v3 parameter set: SPAN — the search would write audio over it");
            fc_loudness_request rq = req;
            std::vector<float> rqOut (4800u * kNch, 0.0f);
            ok (fc_master_solve_delivered (h, &base, &rq, in.data(), 4800u, (float*) (void*) &rq.minPlrDb, 4800u, &so) == FC_ERR_SPAN,
                "and one that starts inside the request: SPAN");
            (void) fc_master_destroy (h);

            fc_master_config lc = deliveringConfig (48000.0, 96000.0);
            ok (fc_master_create (&lc, &h) == FC_OK, "PRECONDITION: a 48 -> 96 kHz handle");
            auto prog = tone (192000u, kNch);
            double first = -1.0, second = -2.0;
            std::vector<float> out ((std::size_t) 384000 * kNch, 0.0f);
            fc_solution sv = 0;
            fc_master_params p = goodParams();
            ok (fc_master_measure_lra (h, prog.data(), 192000u, &first) == FC_OK
                && fc_master_solve_delivered (h, &p, &req, prog.data(), 192000u, out.data(), 384000u, &sv) == FC_OK
                && fc_master_measure_lra (h, prog.data(), 192000u, &second) == FC_OK && first == second,
                "measure, solve, measure: the range is the same number before and after a search");
            (void) fc_solution_destroy (sv);
            for (const std::uint32_t ver : { 1u, 2u, 3u })
            {
                std::vector<unsigned char> rb (96u, 0xC3);
                auto* r = reinterpret_cast<fc_master_resolved*> (static_cast<void*> (rb.data()));
                const std::uint32_t size = ver < 3u ? 80u : 88u;
                r->header.abiVersion = ver; r->header.structSize = size;
                bool tail = true;
                const fc_status st1 = fc_master_resolved_get (h, r);
                for (std::size_t i = size; i < rb.size(); ++i) tail = tail && rb[i] == 0xC3;
                const bool stamp = r->header.abiVersion == ver && r->header.structSize == size;
                r->header.abiVersion = ver; r->header.structSize = size;
                const fc_status st2 = fc_master_configure (h, &p, r);
                for (std::size_t i = size; i < rb.size(); ++i) tail = tail && rb[i] == 0xC3;
                ok (st1 == FC_ERR_STATE || st1 == FC_OK, "PRECONDITION: resolved_get answered");
                ok (st1 == FC_OK && st2 == FC_OK && stamp && tail,
                    "a v" + std::to_string (ver) + " resolved through resolved_get and configure: accepted, stamp echoed, nothing past its "
                    + std::to_string (size) + " bytes");
            }
            (void) fc_master_destroy (h);
        }

        // A FULL TABLE refuses a delivered search before the solver prepares anything.
        {
            fc_master_config c = deliveringConfig (48000.0, 44100.0);
            fc_master h = 0;
            ok (fc_master_create (&c, &h) == FC_OK, "PRECONDITION: a delivering handle, solver not prepared");
            fc_master live[7] {};
            int made = 0;
            const fc_master_config gc = goodConfig();
            for (auto& l : live) if (fc_master_create (&gc, &l) == FC_OK) ++made;
            ok (made == 7, "PRECONDITION: the table is full");
            fc_master_params p = goodParams();
            auto in = tone (4800, kNch);
            std::vector<float> out ((std::size_t) 4410 * kNch, 0.0f);
            fc_solution so = 0;
            fc_need before {}; FC_INIT (before);
            (void) fc_master_need (h, FC_NEED_SOLVE, 4800u, &before);
            const long long b0 = alloc::bytes.load();
            const fc_status st = fc_master_solve_delivered (h, &p, &req, in.data(), 4800u, out.data(), 4410u, &so);
            const long long got = alloc::bytes.load() - b0;
            fc_need afterN {}; FC_INIT (afterN);
            (void) fc_master_need (h, FC_NEED_SOLVE, 4800u, &afterN);
            ok (st == FC_ERR_EXHAUSTED && got == 0 && before.solverPrepared == 0 && afterN.solverPrepared == 0,
                "EXHAUSTED, nothing allocated, and the solver still unprepared (" + std::to_string (got) + " B)");
            for (int i = 0; i < made; ++i) (void) fc_master_destroy (live[i]);
            (void) fc_master_destroy (h);
        }

        // need_create ECHOES its OUT header too, whatever the config's version.
        {
            const fc_master_config c = deliveringConfig (48000.0, 96000.0);
            fc_need nd {}; nd.header.abiVersion = 1u; nd.header.structSize = 40u;
            ok (fc_master_need_create (&c, &nd) == FC_OK && nd.header.abiVersion == 1u && nd.header.structSize == 40u,
                "need_create: a v1-stamped budget comes back stamped v1");
        }

        // compressorMix: the infinities are refused like a NaN, not clamped to the core's 1.
        {
            fc_master h = make();
            fc_master_resolved r {}; FC_INIT (r);
            fc_master_params p {}; FC_INIT (p);
            (void) fc_master_params_defaults (&p);
            p.compressorMix = std::numeric_limits<double>::infinity();
            ok (fc_master_configure (h, &p, &r) == FC_ERR_NON_FINITE, "compressorMix +inf: NON_FINITE");
            p.compressorMix = -std::numeric_limits<double>::infinity();
            ok (fc_master_configure (h, &p, &r) == FC_ERR_NON_FINITE, "and -inf");
            (void) fc_master_destroy (h);
        }

        // THE CREATE BUDGET ON A DOWNSAMPLE AND AT EQUAL RATES — the chain at the delivery rate, not at the higher one.
        for (const auto& pr : { std::pair<double, double> { 192000.0, 44100.0 }, std::pair<double, double> { 48000.0, 48000.0 } })
        {
            fc_master_config c = deliveringConfig (pr.first, pr.second);
            fc_need nd {}; FC_INIT (nd);
            ok (fc_master_need_create (&c, &nd) == FC_OK, "PRECONDITION: budgeted");
            alloc::plainObjectSize.store ((std::size_t) nd.facadeBytes, std::memory_order_relaxed);
            fc_master h = 0;
            const long long b0 = alloc::bytes.load();
            const fc_status cs = fc_master_create (&c, &h);
            const long long got = alloc::bytes.load() - b0;
            ok (cs == FC_OK && got == (long long) (nd.callBytes + nd.facadeBytes),
                std::to_string ((int) pr.first) + " -> " + std::to_string ((int) pr.second)
                + ": the create allocates its budget, to the byte (" + std::to_string (got) + ")");
            (void) fc_master_destroy (h);
            alloc::plainObjectSize.store (0, std::memory_order_relaxed);
        }

        // A NON-FINITE INPUT SAMPLE on a delivering handle costs what it costs a plain one: converting a NaN is
        // converting a zero in its place, bit for bit, and it is counted once — not smeared over a kernel.
        {
            fc_master_config c = deliveringConfig (44100.0, 48000.0);
            fc_master h = 0;
            ok (fc_master_create (&c, &h) == FC_OK, "PRECONDITION: a 44.1 -> 48 kHz handle");
            fc_master_params p = goodParams();
            fc_master_resolved r {}; FC_INIT (r);
            ok (fc_master_configure (h, &p, &r) == FC_OK, "PRECONDITION: configured");
            auto clean = tone (44100, kNch);
            clean[1000] = 0.0f; clean[44100 + 20000] = 1.0e6f;
            auto dirty = clean;
            dirty[1000] = std::numeric_limits<float>::quiet_NaN();
            dirty[44100 + 20000] = std::numeric_limits<float>::infinity();   // not finite: replaced by 0, like a NaN
            clean[44100 + 20000] = 0.0f;
            dirty[44100 + 30000] = 3.0e30f;                                  // finite, far outside: clamps to 1e6
            clean[44100 + 30000] = 1.0e6f;
            std::vector<float> oc ((std::size_t) 48000 * kNch, 0.0f), od (oc.size(), 0.0f);
            ok (fc_master_render_delivered (h, clean.data(), 44100u, oc.data(), 48000u) == FC_OK
                && fc_master_render_delivered (h, dirty.data(), 44100u, od.data(), 48000u) == FC_OK, "both rendered");
            std::size_t differ = 0, nonFinite = 0;
            for (std::size_t i = 0; i < oc.size(); ++i)
            {
                if (std::memcmp (&oc[i], &od[i], sizeof (float)) != 0) ++differ;
                if (! std::isfinite (od[i])) ++nonFinite;
            }
            ok (differ == 0 && nonFinite == 0, "NaN, inf and 3e30 render exactly as 0, 0 and 1e6 do (" + std::to_string (differ)
                + " samples differ)");
            fc_master_stats st {}; FC_INIT (st);
            ok (fc_master_get_stats (h, &st) == FC_OK && st.nonFiniteIn == 2u,
                "and the two non-finite samples are counted as two (" + std::to_string (st.nonFiniteIn) + ")");

            // THE COUNT IS THE LAST CALL'S THAT REACHED IT — the core's contract, read through the facade. Refused before
            // its count, by this facade or by the core, a call leaves the number where it was; the next call that
            // reaches its count replaces it.
            const auto count = [&] { fc_master_stats s {}; FC_INIT (s); (void) fc_master_get_stats (h, &s); return s.nonFiniteIn; };
            ok (fc_master_render_delivered (h, dirty.data(), 44100u, dirty.data(), 48000u) == FC_ERR_SPAN && count() == 2u,
                "a render refused by the facade on its spans leaves the count at 2");
            ok (fc_master_render_delivered (h, clean.data(), 44100u, oc.data(), 47999u) == FC_ERR_CAPACITY && count() == 2u,
                "and one refused on its length leaves it too");
            auto poisoned = tone (44100 * 4, kNch);                           // long enough for a range to be asked for
            poisoned[777] = std::numeric_limits<float>::quiet_NaN();
            double lra = -1.0;
            ok (fc_master_measure_lra (h, poisoned.data(), 44100u * 4u, &lra) == FC_OK && count() == 1u,
                "a range measurement that reaches its count replaces it: 1 (" + std::to_string (count()) + ")");
            ok (fc_master_render_delivered (h, clean.data(), 44100u, oc.data(), 48000u) == FC_OK && count() == 0u,
                "and so does the next render: 0");
            (void) fc_master_destroy (h);

            // AT EQUAL RATES the range is measured on the caller's samples in place, so a poisoned programme is REFUSED —
            // after its count, which the call keeps.
            fc_master_config ce = deliveringConfig (48000.0, 48000.0);
            fc_master he = 0;
            ok (fc_master_create (&ce, &he) == FC_OK && fc_master_configure (he, &p, &r) == FC_OK,
                "PRECONDITION: an equal-rate delivering handle");
            auto poisoned48 = tone (48000 * 4, kNch);
            poisoned48[777] = std::numeric_limits<float>::quiet_NaN();
            fc_master_stats se {}; FC_INIT (se);
            ok (fc_master_measure_lra (he, poisoned48.data(), 48000u * 4u, &lra) == FC_ERR_REFUSED_BY_CORE
                && fc_master_get_stats (he, &se) == FC_OK && se.nonFiniteIn == 1u,
                "at equal rates a range measurement refused AFTER its count keeps its own: 1");
            (void) fc_master_destroy (he);
        }
    }

    //==========================================================================
    // v8 — THE GAIN-REDUCTION PERCENTILE: the `q` each limit binds, and `fc_solution_gr_quantile`. Two things
    // are pinned here and nowhere else — that the two new fields REACH the core, a mapping dropped on either of
    // them leaving a valid default behind, and the entry point's own check order.
    group ("v8: the gain-reduction quantiles — the request's two fields, and fc_solution_gr_quantile");
    {
        using namespace felitronics::mastering;

        // THE DEFAULTS WRITER carries the core's own 0.95, and a request older than v8 is read with it.
        {
            fc_loudness_request d {}; FC_INIT (d);
            const LoudnessRequest cd;
            ok (fc_loudness_request_defaults (&d) == FC_OK
                && d.limiterGrQuantile == cd.limiterGr.quantile && d.compressorGrQuantile == cd.compressorGr.quantile
                && d.limiterGrQuantile == 0.95,
                "the versioned defaults writer: both quantiles are the core's 0.95");
            fc_loudness_request v1 {}; fc_loudness_request_default (&v1);
            ok (v1.header.abiVersion == 1u && v1.header.structSize == 120u,
                "PRECONDITION: the frozen writer still stamps v1 at 120 bytes");
        }

        // THE FIELDS REACH THE CORE, proven by a refusal only the core makes: `q` outside (0, 1] is
        // `InvalidRequest`, which crosses back as FC_OK plus a verdict. A dropped mapping would leave 0.95 in the
        // core's request and the solve would run.
        {
            const std::uint32_t n = 48000u;
            auto in = tone (n, kNch);
            std::vector<float> out (in.size(), 0.0f);
            fc_master_params p {}; FC_INIT (p); (void) fc_master_params_defaults (&p);
            const double bad[] = { 0.0, -0.25, 1.0000001, std::numeric_limits<double>::quiet_NaN(),
                                   std::numeric_limits<double>::infinity() };
            for (const double q : bad)
                for (int which = 0; which < 2; ++which)
                {
                    fc_master h = make();
                    fc_loudness_request req {}; FC_INIT (req); (void) fc_loudness_request_defaults (&req);
                    req.targetLufs = -14.0; req.maxTruePeakDbTp = -1.0;
                    (which == 0 ? req.limiterGrQuantile : req.compressorGrQuantile) = q;
                    fc_solution sol = 0;
                    fc_solution_summary sum {}; FC_INIT (sum);
                    const fc_status st = fc_master_solve (h, &p, &req, in.data(), out.data(), n, &sol);
                    const bool refused = st == FC_OK && fc_solution_summary_get (sol, &sum) == FC_OK
                                      && sum.status == FC_SOLVE_INVALID_REQUEST && sum.passes == 0;
                    ok (refused, std::string (which == 0 ? "limiterGrQuantile" : "compressorGrQuantile") + " = "
                                 + std::to_string (q) + ": the core refuses the request (a dropped mapping would solve)");
                    if (st == FC_OK) (void) fc_solution_destroy (sol);
                    (void) fc_master_destroy (h);
                }
            // A v6 request carries them past its 128 bytes, where they are not read — the trap rule 8 names.
            {
                fc_master h = make();
                fc_loudness_request req {}; FC_INIT (req); (void) fc_loudness_request_defaults (&req);
                req.targetLufs = -14.0; req.maxTruePeakDbTp = -1.0; req.limiterGrQuantile = -1.0;
                req.header.abiVersion = 6u; req.header.structSize = 128u;
                fc_solution sol = 0;
                fc_solution_summary sum {}; FC_INIT (sum);
                ok (fc_master_solve (h, &p, &req, in.data(), out.data(), n, &sol) == FC_OK
                    && fc_solution_summary_get (sol, &sum) == FC_OK && sum.status != FC_SOLVE_INVALID_REQUEST,
                    "a v6 request: the quantile past its 128 bytes is not read, and 0.95 is what the core gets");
                (void) fc_solution_destroy (sol);
                (void) fc_master_destroy (h);
            }
        }

        // THE READ-BACK IS THE CORE'S OWN NUMBER, at every fraction and on both stages, and it is the number an
        // FC_GR_PERCENTILE limit was judged by. The direct solve is the same request through C++.
        {
            const std::uint32_t n = 4u * 48000u;
            auto in = tone (n, kNch);
            std::vector<float> out (in.size(), 0.0f), outC (in.size(), 0.0f);
            fc_master_params p {}; FC_INIT (p); (void) fc_master_params_defaults (&p);
            fc_loudness_request req {}; FC_INIT (req); (void) fc_loudness_request_defaults (&req);
            req.targetLufs = -9.0; req.maxTruePeakDbTp = -1.0; req.maxPasses = 4;
            req.limiterGr.limitDb = 2.0; req.limiterGr.statistic = FC_GR_PERCENTILE;
            req.limiterGrQuantile = 0.62; req.compressorGrQuantile = 0.31;

            fc_master h = make();
            fc_solution sol = 0;
            const fc_status st = fc_master_solve (h, &p, &req, in.data(), out.data(), n, &sol);
            fc_solution_summary sum {}; FC_INIT (sum);
            fc_measurement meas {}; FC_INIT (meas);
            const bool solved = st == FC_OK && fc_solution_summary_get (sol, &sum) == FC_OK
                             && fc_solution_measurement (sol, &meas) == FC_OK;
            ok (solved && sum.passes > 0 && meas.limiter.valid, "PRECONDITION: the search rendered and measured");

            // The same thing directly, through the C++ the facade forwards to.
            MasteringChain chain; OfflineRenderer r; TargetLoudnessSolver solver;
            MasteringChainConfig cc {}; MasteringChainParams cp {};
            const fc_master_config gc = goodConfig();
            cc.internalBlock = gc.internalBlock; cc.eq = gc.eq != 0; cc.monoBass = gc.monoBass != 0;
            cc.compressor = gc.compressor != 0; cc.clipper = gc.clipper != 0; cc.limiter = gc.limiter != 0;
            cc.dither = gc.dither != 0;
            cc.compressorLookaheadMs = gc.compressorLookaheadMs; cc.limiterLookaheadMs = gc.limiterLookaheadMs;
            cc.oversampleFactor = gc.oversampleFactor; cc.tapsPerPhase = gc.tapsPerPhase;
            cc.sidechainHpfHz = gc.sidechainHpfHz;
            LoudnessRequest lr;
            lr.targetLufs = -9.0; lr.maxTruePeakDbTp = -1.0; lr.maxPasses = 4;
            lr.limiterGr.limitDb = 2.0; lr.limiterGr.statistic = GrStatistic::Percentile;
            lr.limiterGr.quantile = 0.62; lr.compressorGr.quantile = 0.31;
            const bool built = r.prepare (kNch, 4096) && chain.prepare (gc.sampleRate, kNch, cc)
                            && solver.prepare (gc.sampleRate, kNch, r.blockSize(), chain.internalBlock(),
                                               chain.tapOversampleFactor());
            const float* ip[16] {}; float* op[16] {};
            for (int c = 0; c < kNch; ++c) { ip[c] = in.data() + (std::size_t) c * n; op[c] = outC.data() + (std::size_t) c * n; }
            LoudnessSolution direct;
            if (built) direct = solver.solve (chain, r, cp, ip, op, kNch, (int) n, lr);
            ok (built && direct.passes == sum.passes && direct.status == (MasteringSolveStatus) sum.status,
                "PRECONDITION: the direct C++ solve is the same search");

            int diff = 0, refusedDiff = 0, answered = 0;
            double biggest = 0.0;
            for (const int stage : { FC_GR_STAGE_COMPRESSOR, FC_GR_STAGE_LIMITER })
                for (const double q : { 0.05, 0.31, 0.5, 0.62, 0.95, 1.0 })
                {
                    double want = 0.0, got = 0.0;
                    const bool coreOk = direct.grQuantile (stage == FC_GR_STAGE_LIMITER ? GrStage::Limiter
                                                                                        : GrStage::Compressor, q, want);
                    const fc_status qs = fc_solution_gr_quantile (sol, stage, q, &got);
                    if (coreOk != (qs == FC_OK)) ++refusedDiff;
                    else if (coreOk && std::memcmp (&want, &got, sizeof (double)) != 0) ++diff;
                    if (qs == FC_OK) { ++answered; if (got > biggest) biggest = got; }
                }
            // ANSWERED, and not only agreed on: twelve refusals agree with twelve refusals.
            ok (answered == 12 && biggest > 0.0,
                "PRECONDITION: all twelve readings were ANSWERED and the trace is live (largest "
                + std::to_string (biggest) + " dB)");
            ok (solved && diff == 0 && refusedDiff == 0,
                "twelve readings through the ABI are the core's, bit for bit (" + std::to_string (diff)
                + " differ, " + std::to_string (refusedDiff) + " disagree about answerability)");
            double judged = 0.0;
            ok (solved && fc_solution_gr_quantile (sol, FC_GR_STAGE_LIMITER, req.limiterGrQuantile, &judged) == FC_OK
                && judged >= 0.0 && judged <= meas.limiter.maxDb,
                "and the reading at the limit's own q sits inside the render it describes (" + std::to_string (judged)
                + " dB, maximum " + std::to_string (meas.limiter.maxDb) + ")");

            // THE CHECK ORDER, and `outDb` untouched by every refusal — the rule every out-parameter here follows.
            double sink = 12345.0;
            ok (fc_solution_gr_quantile (sol, FC_GR_STAGE_LIMITER, 0.5, nullptr) == FC_ERR_NULL
                && fc_solution_gr_quantile (sol, 7, 0.5, &sink) == FC_ERR_ENUM
                && fc_solution_gr_quantile (sol, -1, 0.5, &sink) == FC_ERR_ENUM
                && fc_solution_gr_quantile (sol, FC_GR_STAGE_LIMITER, std::numeric_limits<double>::quiet_NaN(), &sink) == FC_ERR_NON_FINITE
                && fc_solution_gr_quantile (sol, FC_GR_STAGE_LIMITER, std::numeric_limits<double>::infinity(), &sink) == FC_ERR_NON_FINITE
                && fc_solution_gr_quantile (sol, FC_GR_STAGE_LIMITER, 0.0, &sink) == FC_ERR_RANGE
                && fc_solution_gr_quantile (sol, FC_GR_STAGE_LIMITER, 1.5, &sink) == FC_ERR_RANGE
                && sink == 12345.0,
                "null, enum, non-finite and range are named apart, and none of them writes `outDb`");
            // THE ENUM BEFORE THE VALUE: a call that is wrong in both ways is answered for the stage.
            ok (fc_solution_gr_quantile (sol, 7, std::numeric_limits<double>::quiet_NaN(), &sink) == FC_ERR_ENUM
                && sink == 12345.0, "a call wrong in both: the stage is the answer, as the header's order says");
            // A misaligned `outDb`, and the handle before everything.
            {
                std::vector<unsigned char> raw (sizeof (double) + 8u, 0u);
                auto* mis = reinterpret_cast<double*> (static_cast<void*> (raw.data() + 1));
                ok (fc_solution_gr_quantile (sol, FC_GR_STAGE_LIMITER, 0.5, mis) == FC_ERR_ALIGNMENT,
                    "a misaligned `outDb`: ALIGNMENT");
            }
            (void) fc_solution_destroy (sol);
            ok (fc_solution_gr_quantile (sol, FC_GR_STAGE_LIMITER, 0.5, &sink) == FC_ERR_HANDLE && sink == 12345.0,
                "a destroyed solution is stale, and still writes nothing");
            (void) fc_master_destroy (h);
        }

        // A SOLUTION THAT RENDERED NOTHING has no distribution and says so, rather than answering 0.
        {
            fc_master h = make();
            fc_loudness_request bad {}; fc_loudness_request_default (&bad);       // no target: InvalidRequest
            auto in = tone (4800, kNch);
            std::vector<float> out (in.size(), 0.0f);
            fc_master_params p {}; FC_INIT (p); (void) fc_master_params_defaults (&p);
            fc_solution sol = 0;
            double sink = -5.0;
            ok (fc_master_solve (h, &p, &bad, in.data(), out.data(), 4800u, &sol) == FC_OK
                && fc_solution_gr_quantile (sol, FC_GR_STAGE_LIMITER, 0.5, &sink) == FC_ERR_REFUSED_BY_CORE
                && fc_solution_gr_quantile (sol, FC_GR_STAGE_COMPRESSOR, 0.5, &sink) == FC_ERR_REFUSED_BY_CORE
                && sink == -5.0,
                "a verdict before any render: both stages refuse and write nothing");
            (void) fc_solution_destroy (sol);
            (void) fc_master_destroy (h);
        }
    }

    //==========================================================================
    // v7 — THE EQ CURVE. Held against audio: `measuredCurve` drives the same bands through `eq::EqEngine` and
    // measures the steady-state gain of a sine at each frequency of the grid.
    group ("v7: fc_master_eq_curve — the curve is the response the EQ runs");
    {
        namespace E = felitronics::eq;
        constexpr double kEqFs = 48000.0;
        // Every frequency is a multiple of 48000/4800 = 10 Hz, which is what makes the oracle's window hold a
        // whole number of periods.
        const std::vector<double> grid { 120.0, 400.0, 1000.0, 3000.0, 9000.0 };
        const std::uint32_t nf = (std::uint32_t) grid.size();

        std::vector<double> got (grid.size(), 0.0);
        std::uint32_t written = 0xFFFFFFFFu;
        auto ask = [&] (const fc_master_params& pp, std::int32_t lane, std::int32_t band)
        {
            got.assign (grid.size(), 0.0);
            written = 0xFFFFFFFFu;
            return fc_master_eq_curve (&pp, kEqFs, lane, band, grid.data(), nf,
                                       got.data(), (std::uint32_t) got.size(), &written);
        };
        auto worstAgainst = [&] (const std::vector<double>& want)
        {
            double w = 0.0;
            for (std::size_t i = 0; i < want.size(); ++i)
                w = std::max (w, std::fabs (got[i] - want[i]));
            return w;
        };

        // ---- fixture LR: two Stereo bands, one Left band, one Right band. No M/S lane, so the Left and Right
        // axes ARE what channel 0 and channel 1 run, and a ONE-channel engine runs the Stereo axis alone.
        fc_master_params lr {}; FC_INIT (lr);
        ok (fc_master_params_defaults (&lr) == FC_OK, "PRECONDITION: a parameter set at this build's version");
        E::BandParams lrCore[4];
        const CurveBand lrSrc[4] = {
            { 0, FC_FILTER_BELL,       1000.0, 1.4,   6.0, 12 },
            { 0, FC_FILTER_HIGH_PASS,   200.0, 0.707, 0.0, 12 },
            { 1, FC_FILTER_BELL,       3000.0, 3.0,   8.0, 12 },
            { 2, FC_FILTER_HIGH_SHELF, 4000.0, 0.7,  -5.0, 12 },
        };
        for (int i = 0; i < 4; ++i) placeBand (lr.eqBands[i], lrCore[i], lrSrc[i]);

        // ---- fixture MS: a Stereo band, a Mid band and a Side band. No L/R lane, so a pure-Mid and a
        // pure-Side drive isolate those two axes exactly.
        fc_master_params ms {}; FC_INIT (ms);
        ok (fc_master_params_defaults (&ms) == FC_OK, "PRECONDITION: the second parameter set");
        E::BandParams msCore[3];
        const CurveBand msSrc[3] = {
            { 0, FC_FILTER_BELL,      9000.0, 1.0, -4.0, 12 },
            { 3, FC_FILTER_LOW_SHELF,  250.0, 0.8,  5.0, 12 },
            { 4, FC_FILTER_BELL,      1000.0, 2.5, -7.0, 12 },
        };
        for (int i = 0; i < 3; ++i) placeBand (ms.eqBands[i], msCore[i], msSrc[i]);

        const double kTol = 0.005;     // dB, against a measurement in steady state
        double worst = 0.0;
        struct Row { const fc_master_params* p; const E::BandParams* bands; int n; Drive drive; std::int32_t axis;
                     const char* what; };
        const Row rows[] = {
            { &lr, lrCore, 4, Drive::Mono,  FC_EQ_AXIS_STEREO, "Stereo: the Stereo lane alone, on a mono bus" },
            { &lr, lrCore, 4, Drive::Left,  FC_EQ_AXIS_LEFT,   "Left: the Stereo lane times the Left lane" },
            { &lr, lrCore, 4, Drive::Right, FC_EQ_AXIS_RIGHT,  "Right: the Stereo lane times the Right lane" },
            { &ms, msCore, 3, Drive::Mid,   FC_EQ_AXIS_MID,    "Mid: a pure-Mid drive" },
            { &ms, msCore, 3, Drive::Side,  FC_EQ_AXIS_SIDE,   "Side: a pure-Side drive" },
        };
        std::vector<std::vector<double>> axisCurve;
        for (const Row& rw : rows)
        {
            const std::vector<double> meas = measuredCurve (rw.bands, rw.n, kEqFs, grid, rw.drive);
            const fc_status st = ask (*rw.p, rw.axis, -1);
            axisCurve.push_back (got);
            bool finite = true;
            for (double v : meas) finite = finite && std::isfinite (v);
            const double w = worstAgainst (meas);
            worst = std::max (worst, w);
            ok (st == FC_OK && written == nf && finite && w <= kTol,
                std::string ("measured null — ") + rw.what + " (worst " + std::to_string (w) + " dB)");
        }
        std::printf ("      worst curve-vs-measurement error over %zu points: %.6g dB\n",
                     grid.size() * (sizeof rows / sizeof rows[0]), worst);

        // THE FIVE AXES ARE FIVE CURVES. Fixture LR places its bands on the Stereo, Left and Right lanes, so
        // the three differ where they are placed; fixture MS does the same for Mid and Side.
        {
            const std::size_t at120 = 0, at1k = 2, at3k = 3, at9k = 4;
            ok (std::fabs (axisCurve[1][at3k] - axisCurve[0][at3k]) > 1.0
                && std::fabs (axisCurve[2][at9k] - axisCurve[0][at9k]) > 1.0
                && std::fabs (axisCurve[1][at3k] - axisCurve[2][at3k]) > 1.0,
                "Stereo, Left and Right are three different curves where their bands sit");
            ok (std::fabs (axisCurve[3][at1k] - axisCurve[4][at1k]) > 1.0,
                "and Mid and Side at 1 kHz");
            const fc_status st = ask (ms, FC_EQ_AXIS_STEREO, -1);
            ok (st == FC_OK && std::fabs (got[at120] - axisCurve[3][at120]) > 1.0
                && std::fabs (got[at1k] - axisCurve[4][at1k]) > 1.0,
                "the Stereo axis of a set with Mid and Side lanes is neither of them");
        }

        // ONE BAND AT A TIME. Each per-band curve is the response of that band alone — measured for two of
        // them — and the whole bank is their product, which in dB is their sum.
        {
            const std::vector<double> only0 = measuredCurve (&lrCore[0], 1, kEqFs, grid, Drive::Mono);
            ok (ask (lr, FC_EQ_AXIS_STEREO, 0) == FC_OK && worstAgainst (only0) <= kTol,
                "band 0 alone on the Stereo axis matches an engine carrying only that band");
            const std::vector<double> only2 = measuredCurve (&lrCore[2], 1, kEqFs, grid, Drive::Left);
            ok (ask (lr, FC_EQ_AXIS_LEFT, 2) == FC_OK && worstAgainst (only2) <= kTol,
                "band 2 alone on the Left axis, likewise");

            std::vector<double> sum (grid.size(), 0.0);
            for (std::int32_t b = 0; b < FC_MAX_EQ_BANDS; ++b)
            {
                ok (ask (lr, FC_EQ_AXIS_LEFT, b) == FC_OK && written == nf,
                    std::string ("band ") + std::to_string (b) + " answers on its own");
                for (std::size_t i = 0; i < sum.size(); ++i) sum[i] += got[i];
            }
            ok (ask (lr, FC_EQ_AXIS_LEFT, -1) == FC_OK && worstAgainst (sum) <= 1e-9,
                "band = -1 is the sum in dB of all 24 bands — the ones that are off contributing exactly 0");
        }

        // A parameter set with nothing on is exactly unity at every frequency.
        {
            fc_master_params flat {}; FC_INIT (flat);
            ok (fc_master_params_defaults (&flat) == FC_OK, "PRECONDITION: the core's own defaults");
            double m = 0.0;
            for (std::int32_t axis = FC_EQ_AXIS_STEREO; axis <= FC_EQ_AXIS_SIDE; ++axis)
            {
                ok (ask (flat, axis, -1) == FC_OK, "a default parameter set answers on every axis");
                for (double v : got) m = std::max (m, std::fabs (v));
            }
            ok (m <= 1e-12, "and every band being off is 0 dB everywhere, on all five");
        }

        // `bypassEq` and the dynamics are not read — the contract says so, so it is pinned.
        {
            fc_master_params q = lr;
            q.bypassEq = 1;
            for (int b = 0; b < FC_MAX_EQ_BANDS; ++b)
            { q.eqBands[b].dyn.on = 1; q.eqBands[b].dyn.rangeDb = -12.0; }
            ok (ask (q, FC_EQ_AXIS_LEFT, -1) == FC_OK && worstAgainst (axisCurve[1]) <= 0.0,
                "bypassEq and dyn move nothing: the curve is the same bits");

            // AND A `dyn` THE CALLER NEVER FILLED. A block this call does not read may not refuse it either, on
            // or off — otherwise a page drawing a curve has to supply the one struct it was told to ignore.
            for (const int on : { 0, 1 })
            {
                fc_master_params z = lr;
                for (int b = 0; b < FC_MAX_EQ_BANDS; ++b)
                {
                    z.eqBands[b].dyn.on      = on;
                    z.eqBands[b].dyn.rangeDb = std::numeric_limits<double>::quiet_NaN();
                    z.eqBands[b].dyn.thrDb   = std::numeric_limits<double>::infinity();
                    z.eqBands[b].dyn.atk     = -std::numeric_limits<double>::infinity();
                    z.eqBands[b].dyn.rel     = std::numeric_limits<double>::quiet_NaN();
                }
                ok (ask (z, FC_EQ_AXIS_LEFT, -1) == FC_OK && written == nf && worstAgainst (axisCurve[1]) <= 0.0,
                    std::string ("a non-finite `dyn` (on = ") + std::to_string (on)
                    + ") neither refuses the call nor moves the curve");
            }
        }

        // A BAND SWITCHED OFF, OR BYPASSED, CONTRIBUTES UNITY — the same curve as one that was never placed.
        {
            fc_master_params q = lr;
            q.eqBands[2].on = 0;                      // the Left band, off
            const fc_status s1 = ask (q, FC_EQ_AXIS_LEFT, -1);
            const std::vector<double> off = got;
            q.eqBands[2].on = 1; q.eqBands[2].bypass = 1;
            const fc_status s2 = ask (q, FC_EQ_AXIS_LEFT, -1);
            double d = 0.0;
            for (std::size_t i = 0; i < got.size(); ++i) d = std::max (d, std::fabs (got[i] - off[i]));
            q.eqBands[2].bypass = 0; q.eqBands[2].lanes[1].bypass = 1;
            const fc_status s3 = ask (q, FC_EQ_AXIS_LEFT, -1);
            double dl = 0.0;
            for (std::size_t i = 0; i < got.size(); ++i) dl = std::max (dl, std::fabs (got[i] - off[i]));
            ok (s1 == FC_OK && s2 == FC_OK && s3 == FC_OK && d <= 0.0 && dl <= 0.0,
                "off, point-bypassed and lane-bypassed are one answer: unity");
            ok (std::fabs (off[3] - axisCurve[1][3]) > 1.0, "PRECONDITION: and that answer is not the placed one");
        }

        // IT ASKS THE HEAP FOR NOTHING.
        {
            (void) ask (lr, FC_EQ_AXIS_LEFT, -1);                  // warm any lazy path first
            const long long before = alloc::count.load();
            for (int i = 0; i < 8; ++i) (void) ask (lr, FC_EQ_AXIS_LEFT, -1);
            okNoAlloc (alloc::count.load() == before, "eight curves allocate nothing");
        }

        // ---- THE REFUSALS. `written` carries a sentinel into every one of them and must still carry it out,
        // and the buffer carries a canary that must survive.
        {
            std::vector<double> buf (grid.size(), -12345.0);
            const std::vector<double> canary = buf;
            std::uint32_t w = 0x5A5A5A5Au;
            auto refused = [&] (fc_status want, fc_status st, const char* what)
            {
                bool intact = true;
                for (std::size_t i = 0; i < buf.size(); ++i) intact = intact && std::fabs (buf[i] - canary[i]) <= 0.0;
                ok (st == want && w == 0x5A5A5A5Au && intact,
                    std::string (what) + ": refused, `written` untouched, nothing written");
            };

            refused (FC_ERR_CAPACITY,
                     fc_master_eq_curve (&lr, kEqFs, FC_EQ_AXIS_MID, -1, grid.data(), nf, buf.data(), nf - 1u, &w),
                     "a cap below the count");
            refused (FC_ERR_RANGE,
                     fc_master_eq_curve (&lr, kEqFs, FC_EQ_AXIS_MID, -1, grid.data(), 0u, buf.data(), nf, &w),
                     "a grid of no frequencies");
            refused (FC_ERR_ENUM,
                     fc_master_eq_curve (&lr, kEqFs, -1, -1, grid.data(), nf, buf.data(), nf, &w),
                     "an axis code below the first");
            refused (FC_ERR_ENUM,
                     fc_master_eq_curve (&lr, kEqFs, FC_EQ_AXIS_SIDE + 1, -1, grid.data(), nf, buf.data(), nf, &w),
                     "an axis code past the last");
            refused (FC_ERR_RANGE,
                     fc_master_eq_curve (&lr, kEqFs, FC_EQ_AXIS_MID, -2, grid.data(), nf, buf.data(), nf, &w),
                     "a band index below -1");
            refused (FC_ERR_RANGE,
                     fc_master_eq_curve (&lr, kEqFs, FC_EQ_AXIS_MID, FC_MAX_EQ_BANDS, grid.data(), nf, buf.data(), nf, &w),
                     "a band index past the last");

            for (const double bad : { std::numeric_limits<double>::quiet_NaN(),
                                      std::numeric_limits<double>::infinity(),
                                      -std::numeric_limits<double>::infinity() })
            {
                refused (FC_ERR_NON_FINITE,
                         fc_master_eq_curve (&lr, bad, FC_EQ_AXIS_MID, -1, grid.data(), nf, buf.data(), nf, &w),
                         "a non-finite sample rate");
                // THE LAST frequency of the grid, so a call that wrote as it went would have written every
                // earlier point before meeting it.
                std::vector<double> g = grid;
                g.back() = bad;
                refused (FC_ERR_NON_FINITE,
                         fc_master_eq_curve (&lr, kEqFs, FC_EQ_AXIS_MID, -1, g.data(), nf, buf.data(), nf, &w),
                         "a non-finite frequency, last in the grid");
            }

            // The `eq` module's own domain: it clamps a band's frequency into [10, 0.49*rate], so a rate under
            // 20.408163... Hz has no such interval, and 3 MHz is its ceiling.
            for (const double bad : { 0.0, -48000.0, 20.0, 3.0e6 + 1.0 })
                refused (FC_ERR_REFUSED_BY_CORE,
                         fc_master_eq_curve (&lr, bad, FC_EQ_AXIS_MID, -1, grid.data(), nf, buf.data(), nf, &w),
                         "a sample rate the eq module will not honour");
            {
                std::vector<double> eight (grid.size(), 0.0);
                std::uint32_t w8 = 0u;
                ok (fc_master_eq_curve (&lr, 8000.0, FC_EQ_AXIS_MID, -1, grid.data(), nf, eight.data(), nf, &w8) == FC_OK
                    && w8 == nf && std::isfinite (eight[0]),
                    "and 8 kHz, which is inside it, answers");
            }

            refused (FC_ERR_NULL,
                     fc_master_eq_curve (nullptr, kEqFs, FC_EQ_AXIS_MID, -1, grid.data(), nf, buf.data(), nf, &w),
                     "a null parameter set");
            refused (FC_ERR_NULL,
                     fc_master_eq_curve (&lr, kEqFs, FC_EQ_AXIS_MID, -1, nullptr, nf, buf.data(), nf, &w),
                     "a null grid");
            refused (FC_ERR_NULL,
                     fc_master_eq_curve (&lr, kEqFs, FC_EQ_AXIS_MID, -1, grid.data(), nf, nullptr, nf, &w),
                     "a null buffer");
            {
                fc_master_params bad = lr; bad.header.abiVersion = FC_MASTER_ABI_VERSION + 1u;
                refused (FC_ERR_ABI_VERSION,
                         fc_master_eq_curve (&bad, kEqFs, FC_EQ_AXIS_MID, -1, grid.data(), nf, buf.data(), nf, &w),
                         "a parameter set from a newer ABI");
                bad = lr; bad.header.structSize = (std::uint32_t) sizeof (bad) - 8u;
                refused (FC_ERR_STRUCT_SIZE,
                         fc_master_eq_curve (&bad, kEqFs, FC_EQ_AXIS_MID, -1, grid.data(), nf, buf.data(), nf, &w),
                         "a parameter set claiming another size");
                bad = lr; bad.eqBands[0].type = FC_FILTER_TILT + 1;
                refused (FC_ERR_ENUM,
                         fc_master_eq_curve (&bad, kEqFs, FC_EQ_AXIS_MID, -1, grid.data(), nf, buf.data(), nf, &w),
                         "a filter type that names nothing");
                bad = lr; bad.eqBands[0].lanes[0].q = std::numeric_limits<double>::quiet_NaN();
                refused (FC_ERR_NON_FINITE,
                         fc_master_eq_curve (&bad, kEqFs, FC_EQ_AXIS_MID, -1, grid.data(), nf, buf.data(), nf, &w),
                         "a non-finite field of the parameter set");
            }

            // NO TWO OF THE THREE MAY TOUCH, and `written` INSIDE THE GRID is the pair that is not a matter of
            // taste: `freqHz` is the caller's `const`, and a `*written` landing in it before the curve is
            // evaluated answers FC_OK for a frequency nobody asked about. Three placements — the low half of the
            // first frequency, its high half (the one that moves the value), and the second frequency.
            {
                std::vector<double> g (2, 1000.0);
                const std::vector<double> gWas = g;
                std::vector<double> outv (2, -1.0);
                for (const std::size_t at : { (std::size_t) 0, (std::size_t) 4, (std::size_t) 8 })
                {
                    auto* w2 = reinterpret_cast<std::uint32_t*> (reinterpret_cast<unsigned char*> (g.data()) + at);
                    const fc_status st = fc_master_eq_curve (&lr, kEqFs, FC_EQ_AXIS_STEREO, 0,
                                                             g.data(), 2u, outv.data(), 2u, w2);
                    bool gridIntact = true, outIntact = true;
                    for (std::size_t i = 0; i < g.size(); ++i)
                        gridIntact = gridIntact && std::fabs (g[i] - gWas[i]) <= 0.0;
                    for (double v : outv) outIntact = outIntact && std::fabs (v + 1.0) <= 0.0;
                    ok (st == FC_ERR_SPAN && gridIntact && outIntact,
                        std::string ("`written` at byte ") + std::to_string (at)
                        + " of the grid: refused, and neither the grid nor the buffer is written");
                }
            }

            // ALIGNMENT AND ALIASING. The grid and the buffer may not touch, and `written` may not point into
            // the buffer — the rule `fc_master_flush` and `fc_solution_log` follow.
            {
                std::vector<double> big (grid.size() + 2u, 0.0);
                refused (FC_ERR_SPAN,
                         fc_master_eq_curve (&lr, kEqFs, FC_EQ_AXIS_MID, -1, big.data(), nf, big.data() + 1, nf, &w),
                         "a grid overlapping the buffer");
                ok (fc_master_eq_curve (&lr, kEqFs, FC_EQ_AXIS_MID, -1, grid.data(), nf,
                                        big.data(), nf, (std::uint32_t*) big.data()) == FC_ERR_SPAN,
                    "`written` inside the buffer: refused");
                auto* mis = reinterpret_cast<double*> (reinterpret_cast<unsigned char*> (big.data()) + 4);
                refused (FC_ERR_ALIGNMENT,
                         fc_master_eq_curve (&lr, kEqFs, FC_EQ_AXIS_MID, -1, grid.data(), nf - 1u, mis, nf - 1u, &w),
                         "a buffer that is not 8-byte aligned");
                refused (FC_ERR_ALIGNMENT,
                         fc_master_eq_curve (&lr, kEqFs, FC_EQ_AXIS_MID, -1, mis, nf - 1u, big.data(), nf, &w),
                         "a grid that is not 8-byte aligned");
                bool untouched = true;
                for (double v : big) untouched = untouched && std::fabs (v) <= 0.0;
                ok (untouched, "and none of those four wrote a byte of the buffer they were handed");
            }
            refused (FC_ERR_NULL,
                     fc_master_eq_curve (&lr, kEqFs, FC_EQ_AXIS_MID, -1, grid.data(), nf, buf.data(), nf, nullptr),
                     "a null `written`");
        }

        // A CALLER STILL AT v1 reads the same curve: the eq bands are v1's own fields.
        {
            fc_master_params old = lr;
            old.header.abiVersion = 1u; old.header.structSize = 6560u;
            ok (ask (old, FC_EQ_AXIS_LEFT, -1) == FC_OK && written == nf && worstAgainst (axisCurve[1]) <= 0.0,
                "a v1-stamped parameter set answers the same bits");
        }
    }

    //==========================================================================
    // v9 — fc_master_eq_dyn_times. The knobs are DEVIATIONS around a value nothing publishes, so this call is the
    // only way to read what the follower is set to. What makes the group worth writing is the ORACLE: `times()`
    // below is the arithmetic of BandBallistics' HEADER PROSE, re-derived, not a second call into the class — and
    // it simplifies where the class does not (the band's tau is `ringMs/ln9`, and the ln9 cancels), so agreement
    // is two constructions meeting rather than one object agreeing with itself.
    group ("v9: fc_master_eq_dyn_times — the milliseconds behind a deviation knob");
    {
        constexpr double kLn9 = 2.1972245773362196;

        // THE ORACLE, from BandBallistics.h's prose: the digital-bandwidth ring, the period floor, the auto rails,
        // the 2^((knob-0.5)*4) deviation, the absolute rails, and release never shorter than attack. The lane rails
        // in front of it are LaneDynamics' ([10, 0.49*fs] and [0.05, 40]), which is the pair the chain applies.
        auto times = [&] (double fs, double laneF, double laneQ, double atkKnob, double relKnob)
        {
            auto fin  = [] (double v, double fb) { return std::isfinite (v) ? v : fb; };
            const double fc = std::clamp (fin (laneF, 1000.0), 10.0, 0.49 * fs);
            const double Q  = std::clamp (fin (laneQ, 1.0), 0.05, 40.0);
            const double tauBandMs = 2.0e3 * Q / (fs * std::sin (2.0 * felitronics::core::kPi * fc / fs));
            const double period    = 2.5 * 1.0e3 / fc;
            const double atkAuto   = std::clamp (std::max (0.35 * tauBandMs, period), 1.0, 300.0);
            const double relAuto   = std::clamp (3.0 * atkAuto, 12.0, 500.0);
            auto mult = [&] (double k) { return std::pow (2.0, (std::clamp (fin (k, 0.5), 0.0, 1.0) - 0.5) * 4.0); };
            const double atk = std::clamp (atkAuto * mult (atkKnob), 0.2, 1000.0);
            double       rel = std::clamp (relAuto * mult (relKnob), 3.0, 1000.0);
            rel = std::max (rel, atk);
            return std::pair<double, double> { atk, rel };
        };

        // A parameter set with ONE point whose lanes sit apart on purpose: lane 0 low, lane 1 high, so "per lane,
        // not per band" is a measurable claim rather than a sentence in the header.
        auto dynParams = [&] (double f0, double q0, double f1, double q1, double atk, double rel)
        {
            fc_master_params p = goodParams();
            p.eqBands[3].on = 1;
            p.eqBands[3].lanes[0].on = 1; p.eqBands[3].lanes[0].freq = f0; p.eqBands[3].lanes[0].q = q0;
            p.eqBands[3].lanes[1].on = 1; p.eqBands[3].lanes[1].freq = f1; p.eqBands[3].lanes[1].q = q1;
            p.eqBands[3].dyn.on = 1; p.eqBands[3].dyn.rangeDb = -6.0;
            p.eqBands[3].dyn.atk = atk; p.eqBands[3].dyn.rel = rel;
            return p;
        };

        double gotA = -1.0, gotR = -1.0;
        auto ask = [&] (const fc_master_params& p, double fs, std::int32_t band, std::int32_t lane)
        {
            gotA = -1.0; gotR = -1.0;
            return fc_master_eq_dyn_times (&p, fs, band, lane, &gotA, &gotR);
        };
        auto agrees = [&] (double wantA, double wantR)
        {
            const double ea = std::fabs (gotA - wantA) / std::max (1.0, std::fabs (wantA));
            const double er = std::fabs (gotR - wantR) / std::max (1.0, std::fabs (wantR));
            return ea <= 1.0e-12 && er <= 1.0e-12;
        };

        // ---- THE VALUE, against the oracle, over the rates and the shapes the product actually asks for ----
        {
            struct Row { double fs, f, q, atk, rel; const char* what; };
            const Row rows[] = {
                { 48000.0,  7000.0,  4.0, 0.5, 0.5, "a de-esser at 7 kHz, both knobs on auto" },
                { 48000.0,   120.0,  1.0, 0.5, 0.5, "a boom tamer at 120 Hz" },
                { 48000.0,   250.0,  0.7, 0.0, 1.0, "fastest attack, slowest release — the knobs at their ends" },
                { 44100.0, 16000.0, 40.0, 0.5, 0.5, "near Nyquist at the highest Q the EQ admits" },
                { 96000.0,    30.0,  0.3, 0.2, 0.8, "a wide low point at 96 kHz" },
                {  8000.0,  3000.0,  2.0, 0.5, 0.5, "the lowest rate a chain will run at" },
                { 192000.0,  9000.0, 12.0, 0.9, 0.1, "a high rate, and a release knob FASTER than its attack" },
            };
            int bad = 0;
            for (const Row& r : rows)
            {
                const fc_master_params p = dynParams (r.f, r.q, 1000.0, 1.0, r.atk, r.rel);
                const auto want = times (r.fs, r.f, r.q, r.atk, r.rel);
                if (! (ask (p, r.fs, 3, 0) == FC_OK && agrees (want.first, want.second))) { ++bad; }
                ok (ask (p, r.fs, 3, 0) == FC_OK && agrees (want.first, want.second),
                    std::string (r.what) + ": " + std::to_string (gotA) + " / " + std::to_string (gotR) + " ms");
            }
            ok (bad == 0, "every row met the independently derived oracle");
        }

        // ---- RELEASE IS NEVER SHORTER THAN ATTACK, knobs notwithstanding — the last line of compute() ----
        {
            const fc_master_params p = dynParams (9000.0, 12.0, 1000.0, 1.0, 1.0, 0.0);
            ok (ask (p, 48000.0, 3, 0) == FC_OK && gotR >= gotA,
                "attack knob at its slowest and release at its fastest still leaves release >= attack");
        }

        // ---- PER LANE, NOT PER BAND. Same `dyn`, two lanes, two answers — and that is the contract's surprise ----
        {
            const fc_master_params p = dynParams (120.0, 1.0, 9000.0, 8.0, 0.5, 0.5);
            const auto w0 = times (48000.0, 120.0, 1.0, 0.5, 0.5);
            const auto w1 = times (48000.0, 9000.0, 8.0, 0.5, 0.5);
            const bool l0 = ask (p, 48000.0, 3, 0) == FC_OK && agrees (w0.first, w0.second);
            const double a0 = gotA;
            const bool l1 = ask (p, 48000.0, 3, 1) == FC_OK && agrees (w1.first, w1.second);
            ok (l0 && l1 && std::fabs (a0 - gotA) > 1.0e-9,
                "one point's two lanes answer their own times — the shared `dyn` does not make them one number");
        }

        // ---- A LANE THAT IS OFF, AND A BAND THAT IS OFF OR BYPASSED, ARE ANSWERED ANYWAY. The core computes a
        // lane's ballistics whatever its switch says (LaneDynamics::setParams loops over every lane), so a refusal
        // here would be this facade inventing a rule the audio does not have.
        {
            fc_master_params p = dynParams (3000.0, 2.0, 500.0, 1.0, 0.5, 0.5);
            const auto want = times (48000.0, 500.0, 1.0, 0.5, 0.5);
            p.eqBands[3].lanes[1].on = 0;
            const bool offLane = ask (p, 48000.0, 3, 1) == FC_OK && agrees (want.first, want.second);
            p.eqBands[3].on = 0; p.eqBands[3].bypass = 1;
            const bool offBand = ask (p, 48000.0, 3, 1) == FC_OK && agrees (want.first, want.second);
            ok (offLane && offBand, "an off lane, and an off and bypassed band, answer the times they would run at");
        }

        // ---- THE RAILS ARE VISIBLE IN THE ANSWER, which is the whole point of a readback. A caller that wrote
        // MILLISECONDS into a knob reads back the slow rail — the failure the site had, made legible.
        {
            const fc_master_params ms   = dynParams (4000.0, 3.0, 1000.0, 1.0, 50.0, 250.0);   // ms into a knob
            const fc_master_params rail = dynParams (4000.0, 3.0, 1000.0, 1.0,  1.0,   1.0);
            const bool a = ask (ms, 48000.0, 3, 0) == FC_OK;
            const double msA = gotA, msR = gotR;
            const bool b = ask (rail, 48000.0, 3, 0) == FC_OK;
            ok (a && b && std::fabs (msA - gotA) <= 0.0 && std::fabs (msR - gotR) <= 0.0,
                "50 and 250 written into the knobs read back as the SLOW rail, bit for bit — not as 50 and 250 ms");

            const fc_master_params below = dynParams (4000.0, 3.0, 1000.0, 1.0, -3.0, -3.0);
            const fc_master_params zero  = dynParams (4000.0, 3.0, 1000.0, 1.0,  0.0,  0.0);
            const bool c = ask (below, 48000.0, 3, 0) == FC_OK;
            const double loA = gotA, loR = gotR;
            const bool d = ask (zero, 48000.0, 3, 0) == FC_OK;
            ok (c && d && std::fabs (loA - gotA) <= 0.0 && std::fabs (loR - gotR) <= 0.0,
                "a knob below 0 reads back as 0 — the fast rail");
        }

        // ---- THE LANE'S OWN RAILS, and they are LaneDynamics', not BandBallistics' wider pair. A reader that
        // asked BandBallistics with the RAW freq would answer for 5 Hz where the probe sits at 10.
        {
            const fc_master_params lo  = dynParams (5.0,  3.0, 1000.0, 1.0, 0.5, 0.5);
            const fc_master_params ten = dynParams (10.0, 3.0, 1000.0, 1.0, 0.5, 0.5);
            const bool a = ask (lo, 48000.0, 3, 0) == FC_OK; const double loA = gotA, loR = gotR;
            const bool b = ask (ten, 48000.0, 3, 0) == FC_OK;
            ok (a && b && std::fabs (loA - gotA) <= 0.0 && std::fabs (loR - gotR) <= 0.0,
                "a lane at 5 Hz answers for 10 Hz — the probe's floor, bit for bit");

            const fc_master_params hiQ = dynParams (2000.0, 1.0e6, 1000.0, 1.0, 0.5, 0.5);
            const fc_master_params q40 = dynParams (2000.0,   40.0, 1000.0, 1.0, 0.5, 0.5);
            const bool c = ask (hiQ, 48000.0, 3, 0) == FC_OK; const double hiA = gotA, hiR = gotR;
            const bool d = ask (q40, 48000.0, 3, 0) == FC_OK;
            ok (c && d && std::fabs (hiA - gotA) <= 0.0 && std::fabs (hiR - gotR) <= 0.0,
                "Q at 1e6 answers for Q 40 — the EQ's own ceiling, bit for bit");

            // And the ceiling moves WITH the rate: 0.49*fs, not a constant. AT Q 40, deliberately — at the Q 3 this
            // check first used, BOTH rates land on the 1 ms auto floor and answer the same number, so the fixture
            // proved nothing while reading as if it had. The rail is what has to be cleared for the claim to be
            // measurable, and the ring time is linear in Q.
            const fc_master_params up = dynParams (1.0e9, 40.0, 1000.0, 1.0, 0.5, 0.5);
            const auto w44 = times (44100.0, 1.0e9, 40.0, 0.5, 0.5);
            const auto w96 = times (96000.0, 1.0e9, 40.0, 0.5, 0.5);
            const bool e = ask (up, 44100.0, 3, 0) == FC_OK && agrees (w44.first, w44.second);
            const double at44 = gotA;
            const bool f = ask (up, 96000.0, 3, 0) == FC_OK && agrees (w96.first, w96.second);
            ok (e && f && std::fabs (at44 - gotA) > 1.0e-9,
                "the frequency ceiling is 0.49*sampleRate and the answer moves with the rate");
        }

        // ---- THE REFUSALS, in the header's order, and NEITHER OUTPUT IS TOUCHED BY ANY OF THEM ----
        {
            const fc_master_params p = dynParams (3000.0, 2.0, 1000.0, 1.0, 0.5, 0.5);
            double a = -7.0, r = -9.0;
            auto refused = [&] (fc_status want, fc_status st, const char* what)
            {
                ok (st == want && std::fabs (a + 7.0) <= 0.0 && std::fabs (r + 9.0) <= 0.0,
                    std::string (what) + " — refused, and neither output is written");
            };

            refused (FC_ERR_NULL,  fc_master_eq_dyn_times (&p, 48000.0, 3, 0, nullptr, &r), "a null attack output");
            refused (FC_ERR_NULL,  fc_master_eq_dyn_times (&p, 48000.0, 3, 0, &a, nullptr), "a null release output");
            refused (FC_ERR_NULL,  fc_master_eq_dyn_times (nullptr, 48000.0, 3, 0, &a, &r), "a null parameter set");
            refused (FC_ERR_SPAN,  fc_master_eq_dyn_times (&p, 48000.0, 3, 0, &a, &a), "both outputs at one address");
            refused (FC_ERR_RANGE, fc_master_eq_dyn_times (&p, 48000.0, -1, 0, &a, &r), "band -1");
            refused (FC_ERR_RANGE, fc_master_eq_dyn_times (&p, 48000.0, FC_MAX_EQ_BANDS, 0, &a, &r), "band past the bank");
            refused (FC_ERR_RANGE, fc_master_eq_dyn_times (&p, 48000.0, 3, -1, &a, &r), "lane -1");
            refused (FC_ERR_RANGE, fc_master_eq_dyn_times (&p, 48000.0, 3, FC_MAX_EQ_LANES, &a, &r), "lane past the point");

            const double nan = std::numeric_limits<double>::quiet_NaN();
            const double inf = std::numeric_limits<double>::infinity();
            refused (FC_ERR_NON_FINITE, fc_master_eq_dyn_times (&p, nan, 3, 0, &a, &r), "a NaN rate");
            refused (FC_ERR_NON_FINITE, fc_master_eq_dyn_times (&p, inf, 3, 0, &a, &r), "an infinite rate");
            // REFUSED, NOT SUBSTITUTED — the header's rule, and the boundary is `>=`, so 8000 itself is audio.
            refused (FC_ERR_REFUSED_BY_CORE, fc_master_eq_dyn_times (&p, 0.0, 3, 0, &a, &r), "a rate of zero");
            refused (FC_ERR_REFUSED_BY_CORE, fc_master_eq_dyn_times (&p, -48000.0, 3, 0, &a, &r), "a negative rate");
            refused (FC_ERR_REFUSED_BY_CORE, fc_master_eq_dyn_times (&p, 7999.0, 3, 0, &a, &r),
                     "a rate below the core's floor");
            ok (fc_master_eq_dyn_times (&p, 8000.0, 3, 0, &a, &r) == FC_OK, "and 8000 itself is answered");
            a = -7.0; r = -9.0;

            // THE HEADER, and the refusal set of `fc_master_configure` reached through the mapping.
            { fc_master_params bad = p; bad.header.abiVersion = FC_MASTER_ABI_VERSION + 1u;
              refused (FC_ERR_ABI_VERSION, fc_master_eq_dyn_times (&bad, 48000.0, 3, 0, &a, &r), "a version this build does not know"); }
            { fc_master_params bad = p; bad.header.structSize = (std::uint32_t) sizeof (bad) + 8u;
              refused (FC_ERR_STRUCT_SIZE, fc_master_eq_dyn_times (&bad, 48000.0, 3, 0, &a, &r), "a size that is not its version's row"); }
            { fc_master_params bad = p; bad.eqBands[3].dyn.atk = nan;
              refused (FC_ERR_NON_FINITE, fc_master_eq_dyn_times (&bad, 48000.0, 3, 0, &a, &r),
                       "a non-finite knob — refused where configure refuses it, so it never reaches the [0,1] rail"); }
            { fc_master_params bad = p; bad.eqBands[3].lanes[0].freq = nan;
              refused (FC_ERR_NON_FINITE, fc_master_eq_dyn_times (&bad, 48000.0, 3, 0, &a, &r), "a non-finite lane frequency"); }

            // ORDER: a call wrong in two ways answers the EARLIER check. `band` is read before `sampleRate`, and a
            // null output before either — one malformed call, one answer, and the header says which.
            ok (fc_master_eq_dyn_times (&p, nan, FC_MAX_EQ_BANDS, 0, &a, &r) == FC_ERR_RANGE,
                "a bad band AND a NaN rate answers RANGE — the band is checked first");
            ok (fc_master_eq_dyn_times (&p, nan, FC_MAX_EQ_BANDS, 0, nullptr, &r) == FC_ERR_NULL,
                "... and a null output beats both");
        }

        // ---- AN OUTPUT INSIDE `params` IS REFUSED, and this is the pair that is not a matter of taste: the store
        // that writes the attack lands in the caller's `const` parameter set BEFORE the release is read out of it.
        {
            fc_master_params p = dynParams (3000.0, 2.0, 1000.0, 1.0, 0.5, 0.5);
            const fc_master_params was = p;
            double r = -9.0;
            auto* inside = reinterpret_cast<double*> (reinterpret_cast<unsigned char*> (&p) + 8);
            const fc_status st = fc_master_eq_dyn_times (&p, 48000.0, 3, 0, inside, &r);
            ok (st == FC_ERR_SPAN && std::memcmp (&p, &was, sizeof (p)) == 0 && std::fabs (r + 9.0) <= 0.0,
                "the attack output pointing into `params`: refused, and the parameter set is untouched");

            double aa = -7.0;
            auto* inside2 = reinterpret_cast<double*> (reinterpret_cast<unsigned char*> (&p) + 16);
            const fc_status st2 = fc_master_eq_dyn_times (&p, 48000.0, 3, 0, &aa, inside2);
            ok (st2 == FC_ERR_SPAN && std::memcmp (&p, &was, sizeof (p)) == 0 && std::fabs (aa + 7.0) <= 0.0,
                "and the release output pointing into `params` likewise");
        }

        // ---- A CALLER STILL AT v1 reads the same pair: the eq bands and their `dyn` are v1's own fields ----
        {
            fc_master_params p = dynParams (7000.0, 4.0, 1000.0, 1.0, 0.25, 0.75);
            const auto want = times (48000.0, 7000.0, 4.0, 0.25, 0.75);
            p.header.abiVersion = 1u; p.header.structSize = 6560u;
            ok (ask (p, 48000.0, 3, 0) == FC_OK && agrees (want.first, want.second),
                "a v1-stamped parameter set answers the same milliseconds");
        }
    }

    //==========================================================================
    // v10 / K11 — `fc_solution_gr_active_stats`: the limiter's statistics over the windows its input reached the
    // gate. The behaviour is pinned in LoudnessSolverTests against a second programme; what is pinned HERE is
    // the SURFACE — which stage is answerable, which refusal each wrong call gets, and that the gate crosses.
    group ("v10: fc_solution_gr_active_stats — the limiter's active-window statistics across the ABI");
    {
        fc_master h = make();
        fc_master_params p = goodParams();
        fc_master_resolved rr {}; FC_INIT (rr);
        (void) fc_master_configure (h, &p, &rr);
        // THE VERSIONED WRITER, and the choice is the point. `fc_loudness_request_default` is FROZEN at v1: it
        // stamps 120 bytes, and a v10 field written after it lands past the stamp where nothing reads it. The
        // first version of this group used it and read the gate back as 0 — which is not a defect, it is rule 6
        // working, and it is pinned as such at the end of this group.
        fc_loudness_request req {}; FC_INIT (req);
        ok (fc_loudness_request_defaults (&req) == FC_OK, "PRECONDITION: a request at THIS version's layout");
        req.targetLufs = -7.0; req.maxTruePeakDbTp = -1.0;
        req.maxPasses = 1; req.initialGainDb = 12.0;
        req.limiterActiveInputDb = -47.5;                       // off its default, so a lost field shows
        const std::size_t frames = (std::size_t) (kFs * 3.0);
        const std::vector<float> in = tone (frames, kNch);
        std::vector<float> out (in.size(), 0.0f);
        fc_solution sol = 0;
        const bool solved = fc_master_solve (h, &p, &req, in.data(), out.data(), (std::uint32_t) frames, &sol) == FC_OK;
        ok (solved && sol != 0, "PRECONDITION: a one-pass solve at a fixed drive");

        fc_gr_active_stats a {}; FC_INIT (a);
        const fc_status st = fc_solution_gr_active_stats (sol, FC_GR_STAGE_LIMITER, &a);
        ok (st == FC_OK && a.windows > 0, "the limiter answers, over " + std::to_string (a.windows) + " windows");
        // THE GATE CROSSED THE ABI. A field that is marshalled and a field that is forgotten are the same
        // number when the value equals its default, which is why the request above moved it.
        ok (st == FC_OK && a.thresholdDb == -47.5,
            "and it echoes the gate it was read at, so a lost field cannot pass for a default");
        ok (st == FC_OK && a.activeWindows <= a.windows,
            "the active count is a subset of the windows (" + std::to_string (a.activeWindows) + " of "
            + std::to_string (a.windows) + ")");

        // THE COMPRESSOR IS A REFUSAL, NOT ZEROES — the distinction this ABI exists to keep. FC_ERR_STATE says
        // "a stage this ABI knows, with no such measurement here"; zeroes would say "the compressor never
        // worked", which is a claim about the audio that nothing measured.
        fc_gr_active_stats c {}; FC_INIT (c);
        c.windows = 0xABCDu;
        ok (fc_solution_gr_active_stats (sol, FC_GR_STAGE_COMPRESSOR, &c) == FC_ERR_STATE && c.windows == 0xABCDu,
            "the compressor is FC_ERR_STATE and the out-struct is untouched — its input is not tapped, so there "
            "is no gated distribution to answer with, and zeroes would be a measurement nobody made");
        ok (fc_solution_gr_active_stats (sol, 7, &a) == FC_ERR_ENUM, "a code this ABI does not define: ENUM");
        ok (fc_solution_gr_active_stats (sol, -1, &a) == FC_ERR_ENUM, "... and a negative one");
        ok (fc_solution_gr_active_stats (sol, FC_GR_STAGE_LIMITER, nullptr) == FC_ERR_NULL, "a null out-struct");
        ok (fc_solution_gr_active_stats (0, FC_GR_STAGE_LIMITER, &a) == FC_ERR_HANDLE, "a handle that is not one");
        { fc_gr_active_stats bad {}; FC_INIT (bad); bad.header.abiVersion = FC_MASTER_ABI_VERSION + 1u;
          ok (fc_solution_gr_active_stats (sol, FC_GR_STAGE_LIMITER, &bad) == FC_ERR_ABI_VERSION,
              "a version this build does not know"); }
        { fc_gr_active_stats bad {}; FC_INIT (bad); bad.header.structSize = (std::uint32_t) sizeof (bad) + 8u;
          ok (fc_solution_gr_active_stats (sol, FC_GR_STAGE_LIMITER, &bad) == FC_ERR_STRUCT_SIZE,
              "a size that is not its version's row"); }
        // THE HANDLE IS CHECKED BEFORE THE STAGE, the header's order: a bad handle AND a bad code is HANDLE.
        ok (fc_solution_gr_active_stats (0, 7, &a) == FC_ERR_HANDLE, "and the handle is checked before the code");

        // RULE 6, ON THIS FIELD: a caller still at v1 gets the DEFAULT gate, whatever it wrote past its stamp.
        // That is the compatibility rule's whole promise — only the caller's bytes are read — and the number it
        // lands on is -60, the documented default, not the -47.5 sitting in memory at offset 144.
        {
            fc_loudness_request old {}; FC_INIT (old);
            (void) fc_loudness_request_defaults (&old);
            old.targetLufs = -7.0; old.maxTruePeakDbTp = -1.0;
            old.maxPasses = 1; old.initialGainDb = 12.0;
            old.limiterActiveInputDb = -47.5;
            old.header.abiVersion = 1u; old.header.structSize = 120u;      // a v1 caller, with v10 bytes after it
            fc_master h1 = make();
            fc_master_params p1 = goodParams();
            fc_master_resolved r1 {}; FC_INIT (r1);
            (void) fc_master_configure (h1, &p1, &r1);
            std::vector<float> o1 (in.size(), 0.0f);
            fc_solution s1 = 0;
            const bool ran = fc_master_solve (h1, &p1, &old, in.data(), o1.data(), (std::uint32_t) frames, &s1) == FC_OK;
            fc_gr_active_stats v1 {}; FC_INIT (v1);
            const fc_status vs = ran ? fc_solution_gr_active_stats (s1, FC_GR_STAGE_LIMITER, &v1) : FC_ERR_HANDLE;
            ok (vs == FC_OK && v1.thresholdDb == -60.0,
                "a v1-stamped request is read with the DEFAULT gate (-60), not the v10 bytes past its size ("
                + std::to_string (v1.thresholdDb) + ")");
            if (s1 != 0) fc_solution_destroy (s1);
            fc_master_destroy (h1);
        }

        if (sol != 0) fc_solution_destroy (sol);
        fc_master_destroy (h);
    }

    //==========================================================================
    // Q2 — `maxPasses = 1` WITH NO DRIVE-BOUND LIMIT IS EXACTLY ONE RENDER. An external search wants an oracle:
    // render once at a chosen gain and hand back the full measurement, without a search of the solver's own. It
    // does not need a new entry point, and this group is the evidence for that rather than a claim about it.
    //
    // WHY IT HOLDS, from the solver: the search is `for (pass = 0; pass < maxPasses; ++pass)`; the delivery
    // re-render is skipped when the reported candidate is the last one rendered, which a single-pass search
    // always is; and the bracket rescue needs `bound.cap.have`, which is set ONLY by a render that broke a
    // `kDriveBound` limit — `limiterGr`, `minPlrDb`, `maxLraLossLu`. With those three off, nothing can break
    // one, so the rescue is unreachable. `passes == maxPasses + 2` is real, and this is the corner where it
    // cannot happen.
    //
    // THE COUNT IS TAKEN FROM OUTSIDE. `summary.passes` is the solver's own counter, and a test that read only
    // that would be the object agreeing with itself; the renders are counted here through the PROGRESS
    // callback, which is a different road into the same fact. And the instrument gets a control: the same
    // material at `maxPasses = 4` must make it count more than one, or "exactly one" is a statement about a
    // counter that never moves.
    group ("Q2: maxPasses = 1 with no drive-bound limit is exactly one render");
    {
        struct Counter { std::int32_t lastStage = -1, lastPass = -1; int passStarts = 0, finalStarts = 0; };
        const fc_progress_fn fn = [] (void* ctx, const fc_progress* e) -> std::int32_t
        {
            auto& c = *static_cast<Counter*> (ctx);
            if (e->stage != c.lastStage || e->pass != c.lastPass)
            {
                c.lastStage = e->stage; c.lastPass = e->pass;
                if (e->stage == FC_PROGRESS_PASS)  ++c.passStarts;
                if (e->stage == FC_PROGRESS_FINAL) ++c.finalStarts;
            }
            return 1;
        };
        // The three drive-bound limits, OFF — and each off-value is the request's own contract, not a guess:
        // `limitDb == +inf` is "no limit", `minPlrDb == -inf` is "no limit", a NaN input range switches the
        // range constraint off.
        auto oracleRequest = [] (std::int32_t passes)
        {
            fc_loudness_request r {}; fc_loudness_request_default (&r);
            r.targetLufs = -16.0; r.maxTruePeakDbTp = -1.0;
            r.limiterGr.limitDb       = std::numeric_limits<double>::infinity();
            r.compressorGr.limitDb    = std::numeric_limits<double>::infinity();
            r.minPlrDb                = -std::numeric_limits<double>::infinity();
            r.maxLraLossLu            = std::numeric_limits<double>::infinity();
            r.inputLoudnessRangeLu    = std::numeric_limits<double>::quiet_NaN();
            r.maxPasses               = passes;
            r.initialGainDb           = 6.0;      // a chosen drive, which is the whole point of the oracle
            return r;
        };

        const std::size_t frames = (std::size_t) (kFs * 4.0);
        const std::vector<float> in = tone (frames, kNch);
        auto runIt = [&] (std::int32_t passes, Counter& c, fc_solution_summary& sum)
        {
            fc_master h = make();
            fc_master_params p = goodParams();
            fc_master_resolved rr {}; FC_INIT (rr);
            (void) fc_master_configure (h, &p, &rr);
            (void) fc_master_set_progress (h, fn, &c);
            fc_loudness_request req = oracleRequest (passes);
            std::vector<float> out (in.size(), 0.0f);
            fc_solution sol = 0;
            const fc_status st = fc_master_solve (h, &p, &req, in.data(), out.data(), (std::uint32_t) frames, &sol);
            FC_INIT (sum);
            if (st == FC_OK) (void) fc_solution_summary_get (sol, &sum);
            if (sol != 0) fc_solution_destroy (sol);
            fc_master_destroy (h);
            return st;
        };

        Counter one {}; fc_solution_summary s1 {};
        const fc_status st1 = runIt (1, one, s1);
        ok (st1 == FC_OK && one.passStarts == 1 && one.finalStarts == 0,
            "one search pass and no final re-render (" + std::to_string (one.passStarts) + " + "
            + std::to_string (one.finalStarts) + " renders seen from the progress road)");
        ok (st1 == FC_OK && s1.passes == 1,
            "... and the solver's own counter agrees: passes == " + std::to_string (s1.passes));

        // THE CONTROL, and it is the half that makes the check above mean anything: the same material and the
        // same request at four passes has to move the counter. Without this, an instrument wired to a constant
        // 1 would read as a proof.
        Counter four {}; fc_solution_summary s4 {};
        const fc_status st4 = runIt (4, four, s4);
        ok (st4 == FC_OK && four.passStarts > 1,
            "CONTROL: the same material at maxPasses = 4 spends " + std::to_string (four.passStarts)
            + " search passes, so the counter is live");
    }

    //==========================================================================
    // K6 — THE DELIVERED RENDER IS ALIGNED WITH ITS INPUT. A consumer comparing input against output block by
    // block (a crest or spectrum loss per 400 ms window) has to know whether output sample n is input sample n.
    // `OfflineRenderer`'s contract says it is — `out[n] = y[n + D]` — and `MasteringChainTests` nulls that for the
    // renderer. What is pinned HERE is the same claim through the DELIVERED path and the ABI, where a second axis
    // exists that is not latency: at `deliveryRate != sampleRate` the two sides have different LENGTHS and
    // different sample grids, so the correspondence is in TIME and the index correspondence is simply gone.
    group ("K6: the delivered render is aligned with its input — no latency offset, and the rest is time, not index");
    {
        auto deliveringConfig = [] (double fs, double dr)
        {
            fc_master_config c {}; FC_INIT (c);
            (void) fc_master_config_defaults (&c);
            c.sampleRate = fs; c.channels = kNch; c.deliveryRate = dr;
            c.clipper = 1;                  // ON so the chain CARRIES latency: an alignment claim on a chain with
            c.limiter = 1;                  // D == 0 would hold for a renderer that never compensated anything
            return c;
        };
        auto allBypassed = [] ()
        {
            fc_master_params p = goodParams();
            p.inputGainDb = 0.0; p.preLimiterGainDb = 0.0;
            p.bypassEq = 1; p.bypassMonoBass = 1; p.bypassCompressor = 1;
            p.bypassClipper = 1; p.bypassLimiter = 1; p.bypassDither = 1;
            return p;
        };
        // A burst with SILENCE either side, so "where is the feature" is a question with an answer, and the energy
        // centroid is the instrument rather than a first-sample-above-a-threshold (the oversamplers' FIRs ring
        // symmetrically about an edge, so an onset index moves with the filter while the centroid does not).
        auto burst = [] (std::size_t frames, int nch, double fs, double atSec, double lenSec)
        {
            std::vector<float> v (frames * (std::size_t) nch, 0.0f);
            const auto from = (std::size_t) (atSec * fs);
            const auto to   = std::min (frames, from + (std::size_t) (lenSec * fs));
            for (int c = 0; c < nch; ++c)
                for (std::size_t i = from; i < to; ++i)
                {
                    const double w = 0.5 - 0.5 * std::cos (2.0 * felitronics::core::kPi
                                                           * (double) (i - from) / (double) (to - from));
                    v[(std::size_t) c * frames + i] =
                        (float) (0.5 * w * std::sin (2.0 * felitronics::core::kPi * (700.0 + 300.0 * c)
                                                     * (double) i / fs));
                }
            return v;
        };
        auto centroidSec = [] (const std::vector<float>& v, std::size_t frames, int nch, double fs)
        {
            double num = 0.0, den = 0.0;
            for (int c = 0; c < nch; ++c)
                for (std::size_t i = 0; i < frames; ++i)
                {
                    const double e = (double) v[(std::size_t) c * frames + i] * (double) v[(std::size_t) c * frames + i];
                    num += e * (double) i; den += e;
                }
            return den > 0.0 ? (num / den) / fs : -1.0;
        };

        // ---- (a) EQUAL RATES, EVERY STAGE BYPASSED: the delivered render is the input, BIT FOR BIT, at the same
        // index. The chain still carries its PDC — a bypassed clipper and limiter hold theirs in a `core::DryAligner`
        // — so this is the strongest available form of "the latency is compensated": not close, identical.
        {
            fc_master_config c = deliveringConfig (kFs, kFs);
            fc_master h = 0;
            const bool made = fc_master_create (&c, &h) == FC_OK;
            fc_master_params p = allBypassed();
            fc_master_resolved rr {}; FC_INIT (rr);
            const bool cfg = made && fc_master_configure (h, &p, &rr) == FC_OK;
            ok (cfg && rr.latencySamples > 0, "PRECONDITION: a delivering handle whose chain carries latency ("
                + std::to_string (rr.latencySamples) + " samples)");

            const std::size_t n = 48000;
            const std::vector<float> in = burst (n, kNch, kFs, 0.25, 0.10);
            std::vector<float> out (in.size(), -1.0f);
            const bool ran = cfg && fc_master_render_delivered (h, in.data(), (std::uint32_t) n,
                                                                out.data(), (std::uint32_t) n) == FC_OK;
            std::size_t differing = 0;
            for (std::size_t i = 0; i < in.size(); ++i) if (! (out[i] == in[i])) ++differing;
            ok (ran && differing == 0,
                "every stage bypassed: the delivered render IS the input, bit for bit, at the same index ("
                + std::to_string (differing) + " samples differ)");
            if (made) (void) fc_master_destroy (h);
        }

        // A BROADBAND burst, because the lag instrument below needs an unambiguous autocorrelation. A TONE burst
        // does not have one: at 700 Hz the correlation peaks every 68.6 samples and the first version of this check
        // reported a best lag of -206, which is -3 periods to within a tenth of a sample. The Hann envelope over
        // 4800 samples cannot outvote three periods of phase, so the instrument was reading the carrier, not the
        // alignment. Deterministic noise (an LCG, a different seed per channel) has one peak and no second-best.
        auto noiseBurst = [] (std::size_t frames, int nch, double fs, double atSec, double lenSec)
        {
            std::vector<float> v (frames * (std::size_t) nch, 0.0f);
            const auto from = (std::size_t) (atSec * fs);
            const auto to   = std::min (frames, from + (std::size_t) (lenSec * fs));
            for (int c = 0; c < nch; ++c)
            {
                std::uint32_t s = 0x9E3779B9u + 0x7F4A7C15u * (std::uint32_t) c;
                for (std::size_t i = from; i < to; ++i)
                {
                    s = s * 1664525u + 1013904223u;
                    const double w = 0.5 - 0.5 * std::cos (2.0 * felitronics::core::kPi
                                                           * (double) (i - from) / (double) (to - from));
                    v[(std::size_t) c * frames + i] = (float) (0.4 * w * ((double) (s >> 8) / 8388608.0 - 1.0));
                }
            }
            return v;
        };

        // ---- (b) EQUAL RATES, THE STAGES RUNNING. A bit-compare is gone once the clipper and the limiter act, so
        // the instrument is the lag that maximises the cross-correlation. It must be 0 — not "small". A SECOND
        // instrument of a different construction stands beside it (the energy centroid), because a single one
        // agreeing with itself is what the tone-burst version above did while it was wrong.
        {
            fc_master_config c = deliveringConfig (kFs, kFs);
            fc_master h = 0;
            fc_master_params p = goodParams();
            p.bypassDither = 1;                        // the RNG would put noise in the silence and blunt the lag
            fc_master_resolved rr {}; FC_INIT (rr);
            const bool cfg = fc_master_create (&c, &h) == FC_OK && fc_master_configure (h, &p, &rr) == FC_OK;

            const std::size_t n = 48000;
            const std::vector<float> in = noiseBurst (n, kNch, kFs, 0.25, 0.10);
            std::vector<float> out (in.size(), 0.0f);
            const bool ran = cfg && fc_master_render_delivered (h, in.data(), (std::uint32_t) n,
                                                                out.data(), (std::uint32_t) n) == FC_OK;
            long bestLag = 0x7FFFFFFF; double best = -1.0;
            const long span = 4L * (long) std::max (1, rr.latencySamples);
            for (long lag = -span; lag <= span; ++lag)
            {
                double acc = 0.0;
                for (std::size_t i = 0; i < n; ++i)
                {
                    const long j = (long) i + lag;
                    if (j < 0 || j >= (long) n) continue;
                    acc += (double) in[i] * (double) out[(std::size_t) j];
                }
                if (acc > best) { best = acc; bestLag = lag; }
            }
            ok (ran && bestLag == 0, "the stages running: the lag of best correlation is 0 over +-"
                + std::to_string (span) + " samples (found " + std::to_string (bestLag) + ")");
            // THE SECOND INSTRUMENT, of a different construction: WHERE THE PEAK OF A LONE IMPULSE LANDS. It reads
            // one index rather than a whole distribution, so nothing about the envelope can move it.
            //
            // NOT the energy centroid, which this check used first and which failed by 3.2 ms while the lag search
            // said 0. That disagreement was the limiter telling the truth, not a misalignment: it attenuates ahead
            // of the peak through its lookahead and stays down through its release, so more of the burst's TAIL is
            // pulled down than its head and the energy really does sit earlier. A centroid measures the envelope a
            // dynamics stage is there to reshape, so it can only answer an alignment question with every stage
            // bypassed — which is where (a) and (c) use it and where it is exact.
            std::vector<float> imp ((std::size_t) n * (std::size_t) kNch, 0.0f);
            const std::size_t at = 20000;
            for (int ch = 0; ch < kNch; ++ch) imp[(std::size_t) ch * n + at] = 0.3f;    // modest: the alignment is
            std::vector<float> impOut (imp.size(), 0.0f);                               // the claim, not the ceiling
            const bool ran2 = cfg && fc_master_render_delivered (h, imp.data(), (std::uint32_t) n,
                                                                 impOut.data(), (std::uint32_t) n) == FC_OK;
            std::size_t peakAt = 0; double peak = -1.0;
            for (std::size_t i = 0; i < n; ++i)
            {
                const double a2 = std::fabs ((double) impOut[i]);
                if (a2 > peak) { peak = a2; peakAt = i; }
            }
            ok (ran2 && peakAt == at, "... and a lone impulse comes back at its own index (" + std::to_string (at)
                + " in, " + std::to_string (peakAt) + " out)");
            (void) fc_master_destroy (h);
        }

        // ---- (c) DIFFERENT RATES: the index correspondence is GONE and the TIME correspondence holds. This is the
        // half a consumer gets wrong — the output is not the input's length, and sample n of one is not sample n of
        // the other, while 0.25 s in is still 0.25 s in.
        {
            constexpr double kPairs[][2] = { { 48000.0, 44100.0 }, { 44100.0, 48000.0 }, { 48000.0, 96000.0 } };
            for (const auto& pr : kPairs)
            {
                const double sr = pr[0], dr = pr[1];
                fc_master_config c = deliveringConfig (sr, dr);
                fc_master h = 0;
                fc_master_params p = allBypassed();
                fc_master_resolved rr {}; FC_INIT (rr);
                const bool cfg = fc_master_create (&c, &h) == FC_OK && fc_master_configure (h, &p, &rr) == FC_OK;

                const std::size_t n = (std::size_t) sr;                     // one second
                std::uint32_t d = 0;
                const bool len = cfg && fc_master_delivered_frames (h, (std::uint32_t) n, &d) == FC_OK;
                const std::vector<float> in = burst (n, kNch, sr, 0.25, 0.10);
                std::vector<float> out ((std::size_t) d * (std::size_t) kNch, 0.0f);
                const bool ran = len && fc_master_render_delivered (h, in.data(), (std::uint32_t) n,
                                                                    out.data(), d) == FC_OK;
                const double tIn  = centroidSec (in, n, kNch, sr);
                const double tOut = centroidSec (out, (std::size_t) d, kNch, dr);
                const bool lengthIsTheConverters =
                    len && (long long) d == felitronics::mastering::DeliveryConverter::deliveredFrames (sr, dr, (long long) n);
                ok (ran && lengthIsTheConverters && std::fabs (tIn - tOut) < 1.0e-3,
                    std::to_string ((int) sr) + " -> " + std::to_string ((int) dr)
                    + " Hz: the delivered length is the converter's (" + std::to_string (d)
                    + " frames, not " + std::to_string (n) + ") and the feature stays at the same TIME ("
                    + std::to_string (tIn) + " s in, " + std::to_string (tOut) + " s out)");
                (void) fc_master_destroy (h);
            }
        }
    }

    //==========================================================================
    // P41 — A CALL THAT NEVER RETURNED POISONS THE INSTANCE. LAST in this file on purpose: the poison belongs to the
    // whole module and is permanent — that is the contract — so nothing may run after it in this binary.
#if defined(__cpp_exceptions)
    group ("a call that never returned: every later status call answers POISONED and touches nothing");
    {
        fc_master h = make();
        fc_master_params p = goodParams();
        fc_master_resolved r {}; FC_INIT (r);
        ok (fc_master_configure (h, &p, &r) == FC_OK, "PRECONDITION: a configured handle");
        ok (fc_master_set_channel_weight (h, 0, 1.0) == FC_OK, "PRECONDITION: the solver is prepared, so the failing "
                                                                "allocation is INSIDE the search, after its first render");
        const std::uint32_t n = 48000;
        auto in = tone ((int) n, kNch);
        std::vector<float> out (in.size(), 0.0f);
        fc_loudness_request req {}; fc_loudness_request_default (&req);
        req.targetLufs = -14.0; req.maxTruePeakDbTp = -1.0;
        req.initialGainDb = 12.0;          // pass 1 renders at a gain the caller never configured — the wasm replay
        fc_solution earlier = 0;
        ok (fc_master_solve (h, &p, &req, in.data(), out.data(), n, &earlier) == FC_OK && earlier != 0,
            "PRECONDITION: a solution handed out BEFORE the poison");
        fc_solution sol = 0;
        bool escaped = false;
        alloc::failNext = true;
        try { (void) fc_master_solve (h, &p, &req, in.data(), out.data(), n, &sol); }
        catch (const std::bad_alloc&) { escaped = true; }
        alloc::failNext = false;
        ok (escaped, "PRECONDITION: the allocation failure escaped the entry point, as an abort would");
        ok (sol == 0, "and no solution handle was handed out");

        std::vector<float> a ((std::size_t) 64 * kNch, 0.25f), b ((std::size_t) 64 * kNch, 0.0f);
        ok (fc_master_process (h, a.data(), b.data(), 64) == FC_ERR_POISONED,
            "process on the same handle: POISONED — it used to render at the search's +12 dB under FC_OK");
        bool untouched = true; for (float v : b) untouched = untouched && v == 0.0f;
        ok (untouched, "and the output buffer is untouched");
        std::int32_t lat = -7;
        ok (fc_master_latency (h, &lat) == FC_ERR_POISONED && lat == -7, "latency: POISONED, out-parameter untouched");
        ok (fc_master_configure (h, &p, &r) == FC_ERR_POISONED, "configure cannot rescue it");
        ok (fc_master_reset (h) == FC_ERR_POISONED, "reset cannot rescue it");
        fc_master_config cfg = goodConfig();
        fc_master h2 = 12345u;
        ok (fc_master_create (&cfg, &h2) == FC_ERR_POISONED && h2 == 12345u,
            "a NEW handle is refused too — the poison is the module's, not the handle's");
        ok (fc_master_solve (h, &p, &req, in.data(), out.data(), n, &sol) == FC_ERR_POISONED, "solve: POISONED");
        ok (fc_master_destroy (h) == FC_ERR_POISONED, "even destroy: the page throws the whole instance away");
        fc_solution_summary sb {}; FC_INIT (sb);
        ok (fc_solution_summary_get (earlier, &sb) == FC_ERR_POISONED, "and a solution handed out before the poison too");
        // EVERY guarded entry point, not a sample (the code-review round removed the guard from fc_solution_log and
        // moved it behind the handle check in process — both passed a sample of seven), and POISON BEFORE HANDLE: an
        // invalid handle after the poison answers 14, not FC_ERR_HANDLE.
        fc_master_resolved rr {}; FC_INIT (rr);
        fc_master_stats sst {};   FC_INIT (sst);
        fc_need nd {};            FC_INIT (nd);
        fc_measurement ms {};     FC_INIT (ms);
        fc_solve_pass lg[4] {};
        std::uint32_t wrote = 7u;
        double lraOut = -1.0;
        ok (fc_master_resolved_get (h, &rr) == FC_ERR_POISONED, "resolved_get: POISONED");
        ok (fc_master_flush (h, b.data(), 64, &wrote) == FC_ERR_POISONED && wrote == 7u, "flush: POISONED, count untouched");
        ok (fc_master_get_stats (h, &sst) == FC_ERR_POISONED, "get_stats: POISONED");
        ok (fc_master_measure_lra (h, in.data(), n, &lraOut) == FC_ERR_POISONED && lraOut == -1.0, "measure_lra: POISONED");
        ok (fc_master_set_channel_weight (h, 0, 1.0) == FC_ERR_POISONED, "set_channel_weight: POISONED");
        ok (fc_master_need (h, FC_NEED_SOLVE, n, &nd) == FC_ERR_POISONED, "need: POISONED");
        {
            // It takes no handle at all, and is still guarded: the poison belongs to the MODULE, and a
            // budget answered after an abandoned call would be arithmetic over objects nobody can vouch for.
            const fc_master_config pc = goodConfig();
            fc_need pn {}; FC_INIT (pn);
            ok (fc_master_need_create (&pc, &pn) == FC_ERR_POISONED && pn.callBytes == 0u,
                "need_create: POISONED, and the budget is untouched");
        }
        ok (fc_solution_measurement (earlier, &ms) == FC_ERR_POISONED, "solution_measurement: POISONED");
        ok (fc_solution_log (earlier, lg, 4, &wrote) == FC_ERR_POISONED && wrote == 7u, "solution_log: POISONED");
        fc_gr_trace_bucket tb[4] {};
        tb[0].samples = 99u;
        ok (fc_solution_gr_trace (earlier, FC_GR_STAGE_LIMITER, tb, 4, &wrote) == FC_ERR_POISONED && wrote == 7u && tb[0].samples == 99u,
            "solution_gr_trace (v4): POISONED, count and buckets untouched");
        ok (fc_solution_destroy (earlier) == FC_ERR_POISONED, "solution_destroy: POISONED");
        // POISON BEFORE HANDLE, entry point by entry point: an INVALID handle after the poison answers 14, never
        // FC_ERR_HANDLE — the code-review round moved the guard behind the handle check in `process` and a single
        // entry point's check could not see it.
        std::int32_t lat0 = -7;
        ok (fc_master_latency (0, &lat0) == FC_ERR_POISONED && lat0 == -7, "latency(0): 14, out-parameter untouched");
        fc_solution sx = 0;
        const bool all = fc_master_process (0, a.data(), b.data(), 64)             == FC_ERR_POISONED
                      && fc_master_flush (0, b.data(), 64, &wrote)                  == FC_ERR_POISONED
                      && fc_master_configure (0, &p, &rr)                           == FC_ERR_POISONED
                      && fc_master_resolved_get (0, &rr)                            == FC_ERR_POISONED
                      && fc_master_get_stats (0, &sst)                              == FC_ERR_POISONED
                      && fc_master_reset (0)                                        == FC_ERR_POISONED
                      && fc_master_destroy (0)                                      == FC_ERR_POISONED
                      && fc_master_measure_lra (0, in.data(), n, &lraOut)           == FC_ERR_POISONED
                      && fc_master_set_channel_weight (0, 0, 1.0)                   == FC_ERR_POISONED
                      && fc_master_solve (0, &p, &req, in.data(), out.data(), n, &sx) == FC_ERR_POISONED
                      && fc_master_need (0, FC_NEED_SOLVE, n, &nd)                  == FC_ERR_POISONED
                      && fc_master_need_create (nullptr, nullptr)                     == FC_ERR_POISONED
                      && fc_solution_summary_get (0, &sb)                           == FC_ERR_POISONED
                      && fc_solution_measurement (0, &ms)                           == FC_ERR_POISONED
                      && fc_solution_log (0, lg, 4, &wrote)                         == FC_ERR_POISONED
                      && fc_solution_destroy (0)                                    == FC_ERR_POISONED
                      && fc_solution_gr_trace (0, FC_GR_STAGE_LIMITER, tb, 4, &wrote) == FC_ERR_POISONED
                      && fc_solution_gr_quantile (0, FC_GR_STAGE_LIMITER, 0.5, &lraOut) == FC_ERR_POISONED   // v8
                      // v2/v3 — the list is EVERY entry point, so a new one joins it (the diverse-testing round
                      // found the five below missing: removing their guard left the suite green)
                      && fc_master_delivered_frames (0, 64u, &wrote)                == FC_ERR_POISONED
                      && fc_master_render_delivered (0, in.data(), 64u, out.data(), 64u) == FC_ERR_POISONED
                      && fc_master_solve_delivered (0, &p, &req, in.data(), n, out.data(), n, &sx) == FC_ERR_POISONED
                      && fc_master_params_defaults (&p)                             == FC_ERR_POISONED
                      && fc_master_config_defaults (nullptr)                        == FC_ERR_POISONED
                      && fc_loudness_request_defaults (&req)                        == FC_ERR_POISONED;
        ok (all, "every status entry point with an INVALID handle answers 14 after the poison, not FC_ERR_HANDLE");
        // The entry points without a status read no instance state and stay callable.
        ok (fc_master_abi_version() == FC_MASTER_ABI_VERSION, "build identity still answers");
        fc_master_params d {}; fc_master_params_default (&d);
        ok (d.header.abiVersion == 1u, "and the frozen defaults writer still writes");
        fc_master_params dv {}; FC_INIT (dv);
        ok (fc_master_params_defaults (&dv) == FC_ERR_POISONED,
            "while the VERSIONED one returns a status, and is guarded like every such entry point");
    }
#endif

    return felitronics::test::report();
}

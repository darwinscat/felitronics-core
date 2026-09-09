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

#include "fc_master_abi.h"

#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <vector>

// The allocation counter, same idiom as the module suites: a global operator new so that "process()
// does not allocate" is COUNTED rather than read off the source.
static std::atomic<long long> g_allocs { 0 };
void* operator new      (std::size_t s) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s ? s : 1); }
void* operator new[]    (std::size_t s) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s ? s : 1); }
void  operator delete   (void* p) noexcept { std::free (p); }
void  operator delete[] (void* p) noexcept { std::free (p); }
void  operator delete   (void* p, std::size_t) noexcept { std::free (p); }
void  operator delete[] (void* p, std::size_t) noexcept { std::free (p); }

using felitronics::test::ok;
using felitronics::test::approx;
using felitronics::test::group;
using felitronics::test::okNoAlloc;

namespace
{

constexpr double kFs  = 48000.0;
constexpr int    kNch = 2;

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

std::vector<float> tone (std::size_t frames, int nch)
{
    std::vector<float> v (frames * (std::size_t) nch, 0.0f);
    for (int c = 0; c < nch; ++c)
        for (std::size_t i = 0; i < frames; ++i)
            v[(std::size_t) c * frames + i] =
                (float) (0.4 * std::sin (2.0 * 3.14159265358979 * 440.0 * (double) i / kFs));
    return v;
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
        ok (fc_master_sizeof_params() == sizeof (fc_master_params), "sizeof(params) as this build sees it");
        ok (fc_master_sizeof_config() == sizeof (fc_master_config), "sizeof(config)");
        ok (fc_master_max_channels() >= 2, "kMaxChannels is at least stereo");
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
        ok (p.header.abiVersion == FC_MASTER_ABI_VERSION && p.header.structSize == sizeof p,
            "the header is stamped, so the result is usable as an argument");
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
        std::int32_t lat = 0;
        ok (fc_master_latency ((fc_master) 0, &lat) == FC_ERR_HANDLE, "and the reverse cannot even be formed");
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
        const long long before = g_allocs.load();
        for (int i = 0; i < 8; ++i)
            (void) fc_master_process (h, buf.data(), buf.data(), 4096);
        std::int32_t lat = 0; (void) fc_master_latency (h, &lat);
        fc_master_stats st {}; FC_INIT (st);
        (void) fc_master_get_stats (h, &st);
        const long long after = g_allocs.load();
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
        ok (got == 0, "and `written` is cleared even on a refusal");
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
        fc_solution_destroy (sol);

        // A search reads its input again on every pass, so rendering in place would make every pass after
        // the first read the previous pass's master.
        ok (fc_master_solve (h, &p, &req, in.data(), in.data(), (std::uint32_t) frames, &sol) == FC_ERR_SPAN,
            "solving IN PLACE is refused");
        ok (fc_master_solve (h, &p, &req, in.data(), in.data() + 4, (std::uint32_t) 64, &sol) == FC_ERR_SPAN,
            "and so is a partial overlap");
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
            fc_solution sol = 0;
            ok (fc_master_solve (h, &p, &req, in.data(), out.data(), (std::uint32_t) frames, &sol) == FC_OK,
                "the solve returns a verdict");
            fc_solution_summary sum {}; FC_INIT (sum);
            (void) fc_solution_summary_get (sol, &sum);
            ok (sum.status == FC_SOLVE_INVALID_REQUEST, "PRECONDITION: it refused BEFORE rendering");
            ok (fc_master_configure (h, &p, &r) == FC_ERR_STATE,
                "and the guard is STILL armed — a refused solve moved nothing");
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
            auto warm = tone (64, kNch);
            (void) fc_master_process (h, warm.data(), warm.data(), 64);

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
            ok (fc_master_set_channel_weight (h, 0, -1.0) == FC_ERR_NON_FINITE, "a negative weight");
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

    return felitronics::test::report();
}

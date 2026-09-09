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

#include <felitronics/mastering/LoudnessSolver.h>
#include <felitronics/mastering/MasteringChain.h>
#include <felitronics/mastering/OfflineRenderer.h>

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
        req.toleranceLu = 0.037;
        req.truePeakAimDb = 0.073;

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
            lr.toleranceLu = 0.037; lr.truePeakAimDb = 0.073;
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

            std::vector<fc_solve_pass> log ((std::size_t) sum.logCount + 1);
            std::uint32_t written = 0;
            (void) fc_solution_log (sol, log.data(), (std::uint32_t) log.size(), &written);
            ok (written == (std::uint32_t) direct.logCount, "the log has the core's length");
            bool logSame = (written > 0);
            for (std::uint32_t i = 0; i < written; ++i)
                if (log[i].gainDb != direct.log[i].gainDb || log[i].ceilingDb != direct.log[i].ceilingDb
                    || log[i].integratedLufs != direct.log[i].integratedLufs
                    || log[i].violated != direct.log[i].violated) logSame = false;
            ok (logSame, "and every pass record crosses unchanged (gain, ceiling, loudness, mask)");
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
            for (std::size_t i = 0; i < frames; ++i)
                if ((i / 48000) % 2 == 0) in[frames + i] *= 0.05f;
            double viaAbi = 0.0;
            ok (fc_master_measure_lra (h, in.data(), (std::uint32_t) frames, &viaAbi) == FC_OK, "measured");

            TargetLoudnessSolver solver;
            ok (solver.prepare (kFs, kNch, 4096, 256, 4), "a direct solver for the same measurement");
            const float* ip[2] { in.data(), in.data() + frames };
            double direct = 0.0;
            ok (solver.measureInputLoudnessRange (ip, kNch, (int) frames, direct), "the core measured too");
            ok (viaAbi == direct, "and the two are the SAME NUMBER, bit for bit");
            ok (direct > 0.5, "PRECONDITION: the fixture has a range, so an added offset would show");
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

    return felitronics::test::report();
}

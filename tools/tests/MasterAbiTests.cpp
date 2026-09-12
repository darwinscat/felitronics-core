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
#include <new>
#include <string>
#include <vector>

// The allocation counter, same idiom as the module suites: a global operator new so that "process()
// does not allocate" is COUNTED rather than read off the source.
static std::atomic<long long> g_allocs { 0 };
// P41: a switch that makes the NEXT allocation fail — natively an exception, the exact analogue of the wasm abort
// (nothing after the failing statement runs, the facade's own bookkeeping included). Exceptions only: the wasm tier
// builds this file with -fno-exceptions, where the analogue is the abort itself and the suite cannot outlive it.
static std::atomic<bool> g_failNextAlloc { false };
// BYTES too: fc_master_need's budgets are held against them, so the counter counts what the CONTAINER asked for. MSVC's
// STL on x86/x64 asks operator new for sizeof(void*) + 31 bytes more (one word more under _DEBUG) on every block of 4096
// or more — its own alignment, which a budget leaves to the caller — and the counter takes that back off. The rule and
// its reasons are written out in LoudnessConformanceTests.cpp; here, as there, a check fails if it is not this STL's —
// and under iterator debugging, which neither file models (a proxy allocation per container), it fails by name.
#if defined(_MSVC_STL_VERSION) && (defined(_M_IX86) || defined(_M_X64))
#  if defined(_DEBUG)
static constexpr std::size_t kStlBigPad = 2 * sizeof (void*) + 31;
#  else
static constexpr std::size_t kStlBigPad = sizeof (void*) + 31;
#  endif
#else
static constexpr std::size_t kStlBigPad = 0;
#endif
static constexpr std::size_t kStlBigBlock = 4096;
static std::atomic<long long> g_bytes { 0 };
// A PLAIN `new T` IS NOT A CONTAINER, and on this row the difference is 39 bytes. The correction above
// undoes what MSVC's STL adds to a CONTAINER's request; a direct `new` asks for exactly `sizeof(T)` and
// taking the correction off it subtracts bytes nobody added. A create makes two such allocations: the EQ
// engine, which is over-aligned and therefore counted raw on its own path, and THIS FILE's instance record,
// which is an ordinary `new` of 18 496 bytes — over the big-block threshold, and indistinguishable from a
// container by size alone. So the test tells the counter that size (it is `fc_need::facadeBytes`, which the
// core publishes) and the counter leaves requests of exactly it alone. The `win` row found this: every one
// of the 48 create rows read 39 bytes short while macOS and Linux, whose STL adds nothing, were exact.
static std::atomic<std::size_t> g_plainObjectSize { 0 };
static long long containerBytes (std::size_t s) noexcept
{
    if (kStlBigPad == 0) return (long long) s;
    if (s == g_plainObjectSize.load (std::memory_order_relaxed)) return (long long) s;
    return (long long) (s >= kStlBigBlock + kStlBigPad ? s - kStlBigPad : s);
}
static void* countedNew (std::size_t s)
{
    g_allocs.fetch_add (1, std::memory_order_relaxed);
    g_bytes.fetch_add (containerBytes (s), std::memory_order_relaxed);
#if defined(__cpp_exceptions)
    if (g_failNextAlloc.exchange (false)) throw std::bad_alloc();
#endif
    return std::malloc (s ? s : 1);
}
// One vector's allocation, as the counter sees it — written through `volatile`, so the optimizer cannot remove it.
static long long vectorRequest (std::size_t n)
{
    const long long before = g_bytes.load();
    {
        std::vector<char> v;
        v.assign (n, 0);
        volatile char* sink = v.data();
        sink[0] = 1;
    }
    return g_bytes.load() - before;
}
// EVERY FORM, not two. `eq::EqEngine` has `alignof` 64, so `fc_master_create` builds it through the
// OVER-ALIGNED `operator new(size_t, align_val_t)` — 331 KiB, the single largest request a create makes,
// invisible to a counter that overrides only the two sized forms. That blindness is P52, and a suite that
// holds a published budget against a counter cannot have it: the create budgets below would have passed
// while silently comparing two numbers that both left the engine out.
static void* countedAlignedNew (std::size_t s, std::size_t a)
{
    g_allocs.fetch_add (1, std::memory_order_relaxed);
    // RAW, not `containerBytes`. That correction undoes what MSVC's STL adds ON TOP of a container's own
    // request; an over-aligned `new` here is an OBJECT (`eq::EqEngine`, alignof 64) whose size is exactly
    // `sizeof`, and taking the correction off it subtracts bytes nobody added. The `win` row found it.
    g_bytes.fetch_add ((long long) s, std::memory_order_relaxed);
#if defined(__cpp_exceptions)
    if (g_failNextAlloc.exchange (false)) throw std::bad_alloc();
#endif
#if defined(_MSC_VER)
    return _aligned_malloc (s ? s : 1, a);
#else
    const std::size_t al = a < sizeof (void*) ? sizeof (void*) : a;
    void* p = nullptr;
    return posix_memalign (&p, al, s ? s : 1) == 0 ? p : nullptr;
#endif
}
// The nothrow forms take the SAME counter and the same fault switch, but answer a failure the way their
// contract does — a null pointer, never an exception out of a noexcept function.
static void* countedNothrowNew (std::size_t s, std::size_t a) noexcept
{
    g_allocs.fetch_add (1, std::memory_order_relaxed);
    g_bytes.fetch_add (a != 0 ? (long long) s : containerBytes (s), std::memory_order_relaxed);   // raw when over-aligned — see above
    if (g_failNextAlloc.exchange (false)) return nullptr;
#if defined(_MSC_VER)
    return a != 0 ? _aligned_malloc (s ? s : 1, a) : std::malloc (s ? s : 1);
#else
    if (a == 0) return std::malloc (s ? s : 1);
    const std::size_t al = a < sizeof (void*) ? sizeof (void*) : a;
    void* p = nullptr;
    return posix_memalign (&p, al, s ? s : 1) == 0 ? p : nullptr;
#endif
}
static void alignedFree (void* p) noexcept
{
#if defined(_MSC_VER)
    _aligned_free (p);
#else
    std::free (p);
#endif
}
void* operator new      (std::size_t s) { return countedNew (s); }
void* operator new[]    (std::size_t s) { return countedNew (s); }
void* operator new      (std::size_t s, std::align_val_t a) { return countedAlignedNew (s, (std::size_t) a); }
void* operator new[]    (std::size_t s, std::align_val_t a) { return countedAlignedNew (s, (std::size_t) a); }
void* operator new      (std::size_t s, const std::nothrow_t&) noexcept { return countedNothrowNew (s, 0); }
void* operator new[]    (std::size_t s, const std::nothrow_t&) noexcept { return countedNothrowNew (s, 0); }
void* operator new      (std::size_t s, std::align_val_t a, const std::nothrow_t&) noexcept { return countedNothrowNew (s, (std::size_t) a); }
void* operator new[]    (std::size_t s, std::align_val_t a, const std::nothrow_t&) noexcept { return countedNothrowNew (s, (std::size_t) a); }
void  operator delete   (void* p) noexcept { std::free (p); }
void  operator delete[] (void* p) noexcept { std::free (p); }
void  operator delete   (void* p, std::size_t) noexcept { std::free (p); }
void  operator delete[] (void* p, std::size_t) noexcept { std::free (p); }
void  operator delete   (void* p, std::align_val_t) noexcept { alignedFree (p); }
void  operator delete[] (void* p, std::align_val_t) noexcept { alignedFree (p); }
void  operator delete   (void* p, std::size_t, std::align_val_t) noexcept { alignedFree (p); }
void  operator delete[] (void* p, std::size_t, std::align_val_t) noexcept { alignedFree (p); }

using felitronics::test::ok;
using felitronics::test::approx;
using felitronics::test::group;
using felitronics::test::okNoAlloc;

namespace
{

constexpr double kFs  = 48000.0;
constexpr int    kNch = 2;

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
        // meter is sized for 192000 + 48000 samples = 50 hops of 4800 → 8·(300 + 54 + 13) = 2936 B; the true-peak
        // meter 4·48 + 4·(2 ch · 2 · 12) + 4·2 = 392 B (the history ring became DOUBLE-LENGTH with P56, so
        // core::firDot reads a contiguous window — P56; it was 296 B); the drain 2·64·4 = 512 B.
        ok (lra.callBytes == 2936u, "the measure_lra budget is the hand-derived 2936 B");
        ok (solve.callBytes == 2936u + 392u + 512u, "the solve budget is meter + true-peak meter + drain = 3840 B");

        long long before = g_bytes.load();
        const fc_status w = fc_master_set_channel_weight (h, 0, 1.0);
        const long long prepared = g_bytes.load() - before;
        ok (w == FC_OK, "the first weight prepares the solver");
        ok (prepared == (long long) solve.solverPrepareBytes, "and allocates exactly `solverPrepareBytes`");
        fc_need after {}; FC_INIT (after);
        ok (fc_master_need (h, FC_NEED_SOLVE, n, &after) == FC_OK && after.solverPrepared == 1,
            "and from then on the budget says the preparation is spent");

        auto in = tone ((int) n, kNch);
        double v = 0.0;
        before = g_bytes.load();
        (void) fc_master_measure_lra (h, in.data(), n, &v);
        const long long lraBytes = g_bytes.load() - before;
        ok (lraBytes == (long long) lra.callBytes, "measure_lra allocates exactly its budget");

        // Under 3 s there is no range: the call refuses BEFORE it builds a meter, and its budget says so (the
        // code-review round found it promising a meter for exactly this call).
        fc_need shortLra {}; FC_INIT (shortLra);
        ok (fc_master_need (h, FC_NEED_MEASURE_LRA, 48000, &shortLra) == FC_OK && shortLra.callBytes == 0,
            "a 1 s programme: the measure_lra budget is 0");
        before = g_bytes.load();
        const fc_status sr = fc_master_measure_lra (h, in.data(), 48000, &v);
        const long long shortBytes = g_bytes.load() - before;
        ok (sr == FC_ERR_REFUSED_BY_CORE && shortBytes == 0, "and the refused call allocates nothing");

        // The range rule's own edge, at 48 kHz: exactly 3 s is measurable, one frame less is not — in the call AND in
        // its budget, which read the same `rangeMeasurable`.
        fc_need at3 {}, under3 {}; FC_INIT (at3); FC_INIT (under3);
        ok (fc_master_need (h, FC_NEED_MEASURE_LRA, 144000, &at3) == FC_OK && at3.callBytes > 0,
            "exactly 3 s: a meter is budgeted");
        ok (fc_master_need (h, FC_NEED_MEASURE_LRA, 143999, &under3) == FC_OK && under3.callBytes == 0,
            "one frame under 3 s: nothing is");
        // A solve builds its meters whatever the length — the range rule is NOT the solve's. 1 s still costs
        // meter + true-peak meter + drain: 8·(300 + 24 + 10) + 392 + 512 = 3576 B. A length the solve refuses costs 0.
        fc_need s1 {}, s0 {}; FC_INIT (s1); FC_INIT (s0);
        ok (fc_master_need (h, FC_NEED_SOLVE, 48000, &s1) == FC_OK && s1.callBytes == 3576u,
            "a 1 s solve is budgeted in full: 3576 B");
        ok (fc_master_need (h, FC_NEED_SOLVE, 0, &s0) == FC_OK && s0.callBytes == 0,
            "a 0-frame solve, which the core refuses before any pass, costs 0");

        std::vector<float> out (in.size(), 0.0f);
        fc_loudness_request req {}; fc_loudness_request_default (&req);
        req.targetLufs = -14.0; req.maxTruePeakDbTp = -1.0;
        fc_solution sol = 0;
        before = g_bytes.load();
        const fc_status sv = fc_master_solve (h, &p, &req, in.data(), out.data(), n, &sol);
        const long long solveBytes = g_bytes.load() - before;
        fc_solution_summary sum {}; FC_INIT (sum);
        ok (sv == FC_OK && fc_solution_summary_get (sol, &sum) == FC_OK && sum.passes > 0, "PRECONDITION: the search rendered");
        const long long perPass = (long long) solve.callBytes - 512;
        ok (solveBytes == (long long) sum.passes * perPass + 512 + (long long) solve.facadeBytes,
            "a solve allocates passes × (both meters) + the drain + the facade's record: its budget's parts");
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
            g_plainObjectSize.store ((std::size_t) pn.facadeBytes, std::memory_order_relaxed);
            // AND THE RULE IS CALIBRATED, on this row's own STL, rather than trusted: a plain `new` of that
            // size counts as exactly that size, and a VECTOR of it counts as exactly its own bytes. The
            // second is the rule's known collision — a container that happens to be exactly as long as the
            // facade's record would be left uncorrected — and naming it here is what keeps it from being
            // discovered as a byte-for-byte failure with no explanation.
            const std::size_t n = (std::size_t) pn.facadeBytes;
            const long long before = g_bytes.load();
            {
                void* raw = ::operator new (n);
                volatile char* sink = static_cast<char*> (raw);
                sink[0] = 1;
                ::operator delete (raw, n);
            }
            const long long got = g_bytes.load() - before;
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
                    const long long before = g_bytes.load();
                    const fc_status cs = fc_master_create (&c, &h);
                    const long long got = g_bytes.load() - before;
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
                    const long long b2 = g_bytes.load();
                    const fc_status ccs = fc_master_configure (h, &p, &r);
                    const long long got2 = g_bytes.load() - b2;
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
            const long long before = g_bytes.load();
            const fc_status cs = fc_master_create (&c, &h);
            const long long got = g_bytes.load() - before;
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
        g_plainObjectSize.store (0, std::memory_order_relaxed);   // the exemption is this group's only
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
        g_failNextAlloc = true;
        try { (void) fc_master_solve (h, &p, &req, in.data(), out.data(), n, &sol); }
        catch (const std::bad_alloc&) { escaped = true; }
        g_failNextAlloc = false;
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
                      && fc_solution_destroy (0)                                    == FC_ERR_POISONED;
        ok (all, "every status entry point with an INVALID handle answers 14 after the poison, not FC_ERR_HANDLE");
        // The entry points without a status read no instance state and stay callable.
        ok (fc_master_abi_version() == FC_MASTER_ABI_VERSION, "build identity still answers");
        fc_master_params d {}; fc_master_params_default (&d);
        ok (d.header.abiVersion == FC_MASTER_ABI_VERSION, "and the defaults writer still writes");
    }
#endif

    return felitronics::test::report();
}

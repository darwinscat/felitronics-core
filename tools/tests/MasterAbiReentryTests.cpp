// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// fc_master — RE-ENTRY POISONS THE MODULE, AND FOR GOOD (P41). A binary of its own because the poison is permanent by
// contract: felitronics_master_abi_tests spends its one poisoning on a call that never returned, and a second scenario
// cannot share a process with it. This one needs no exceptions, so the wasm tier runs it too.
//
// The scenario is the one fc_master_abi.h names under "THE MODULE IS NOT RE-ENTRANT": an entry point called while
// another is still running — here from inside an allocation the module made, a native new_handler's position. The inner
// call cannot be told apart from the first call after an abandoned one, so it answers POISONED; the outer call returns
// normally; and from then on every status call answers POISONED, although every call returned. "For good" is the latch
// in FC_GUARD: without it the module would answer again once the outer call returned, and the 14 the inner call handed
// out would have been a lie.

#include <felitronics_test.h>

#include "fc_master_abi.h"

#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <new>

static std::atomic<bool> g_reenterOnNextAlloc { false };
static fc_master    g_target       = 0;
static fc_status    g_inner        = FC_OK;
static int          g_innerCalls   = 0;
static std::int32_t g_innerLatency = -7;

static void* hookedNew (std::size_t s)
{
    if (g_reenterOnNextAlloc.exchange (false))
    {
        g_inner = fc_master_latency (g_target, &g_innerLatency);
        ++g_innerCalls;
    }
    return std::malloc (s ? s : 1);
}
void* operator new      (std::size_t s) { return hookedNew (s); }
void* operator new[]    (std::size_t s) { return hookedNew (s); }
void  operator delete   (void* p) noexcept { std::free (p); }
void  operator delete[] (void* p) noexcept { std::free (p); }
void  operator delete   (void* p, std::size_t) noexcept { std::free (p); }
void  operator delete[] (void* p, std::size_t) noexcept { std::free (p); }

int main()
{
    using felitronics::test::ok;
    felitronics::test::group ("an entry point re-entered from inside an allocation poisons the module, for good");

    fc_master_config cfg {};
    fc_master_config_default (&cfg);
    cfg.sampleRate = 48000.0;
    cfg.channels   = 2;
    fc_master h = 0;
    ok (fc_master_create (&cfg, &h) == FC_OK && h != 0, "PRECONDITION: a handle");
    fc_master_params p {};
    fc_master_params_default (&p);
    fc_master_resolved r {}; FC_INIT (r);
    ok (fc_master_configure (h, &p, &r) == FC_OK, "PRECONDITION: configured");
    std::int32_t lat = -7;
    ok (fc_master_latency (h, &lat) == FC_OK && lat >= 0, "PRECONDITION: it answers");

    g_target = h;
    fc_master h2 = 0;
    g_reenterOnNextAlloc = true;                        // the next allocation is inside fc_master_create
    const fc_status outer = fc_master_create (&cfg, &h2);
    ok (g_innerCalls == 1, "PRECONDITION: an entry point really was called from inside another");
    ok (g_inner == FC_ERR_POISONED && g_innerLatency == -7, "the inner call answers POISONED and writes nothing");
    ok (outer == FC_OK && h2 != 0, "the outer call returns normally: it cannot know");

    lat = -7;
    ok (fc_master_latency (h, &lat) == FC_ERR_POISONED && lat == -7,
        "and from then on every status call answers POISONED — for good, although every call returned");
    fc_master h3 = 12345u;
    ok (fc_master_create (&cfg, &h3) == FC_ERR_POISONED && h3 == 12345u, "a new handle is refused too");
    ok (fc_master_destroy (h2) == FC_ERR_POISONED, "even destroy");
    ok (fc_master_abi_version() == FC_MASTER_ABI_VERSION, "the build identity still answers");
    return felitronics::test::report();
}

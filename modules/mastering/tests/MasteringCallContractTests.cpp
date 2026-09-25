// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

//==================================================================================================
// LAW 11b (docs/DSP-ARCHITECTURE.md §2) for mastering::OfflineRenderer and
// mastering::MasteringChain — the named case the law-11 census (modules/laws/tests/CallContractTests)
// carried for them, verbatim, through the same harness (test_support/law11_call_contract.h).
// The chain's other law-11 claims (its EXACT width, block invariance) are MasteringChainTests' own.
//==================================================================================================

#include <felitronics_test.h>
#include <law11_call_contract.h>

#include <felitronics/mastering/MasteringChain.h>
#include <felitronics/mastering/OfflineRenderer.h>

#include <cstdio>
#include <vector>

namespace {

using felitronics::test::ok;
using felitronics::test::group;
using namespace felitronics::test::law11;

void refusedPrepareCases()
{
    group ("law 11b — a refused prepare() adopts nothing and leaves the object UNUSABLE");
    {
        // The two halves are only compatible in one order — disarm, validate, write — and the
        // observable contract is these three lines, not "the old fields survive": storing the width
        // before validating the block size was a heap-buffer-overflow in render().
        felitronics::mastering::OfflineRenderer r;
        ok (r.prepare (1, 256), "precondition: a 1-channel, 256-sample renderer");
        ok (! r.prepare (2, 0), "prepare(2, 0) is refused — the block size is impossible");
        ok (r.maxChannels() != 2, "...and the REFUSED call's own width was not adopted");
        felitronics::mastering::MasteringChain chain;
        ok (chain.prepare (kFs, 1), "precondition: a 1-channel chain to render through");
        std::vector<float> buf ((std::size_t) 256, 0.1f);
        float* p[1] { buf.data() };
        ok (! r.render (chain, (const float* const*) p, p, 1, 256),
            "...and the renderer is UNUSABLE until a prepare() succeeds (it used to write past its scratch)");
        ok (r.prepare (1, 256) && r.render (chain, (const float* const*) p, p, 1, 256),
            "...and a successful prepare() brings it back");
    }
}

} // namespace

int main()
{
    std::printf ("felitronics law 11 — the caller's contract: mastering\n");

    refusedPrepareCases();

    return felitronics::test::report();
}

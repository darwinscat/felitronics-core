// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// P88 (law 11e). The one hazard of the new reset() that no single-threaded gate can reach: an
// UNCONDITIONAL Idle store would wipe a publication that lands between reset()'s load of the state and its
// store — the same loss the fix removes, reached through a race instead of a sequence. The stores are made
// only when a swap was in flight, and this is the gate that says so.
//
// Two real threads. The audio thread spins reset(); once it is really spinning the loader publishes ONCE; then
// two more restarts, which must adopt whatever was accepted. A two-sample DC probe, made while both threads
// are parked, reads the live head tap, and it must be the head of what was just published — a literal from
// outside the class, a different one every round.
// Measured: the bodies before P88 lose every publication. The fixed bodies with the Idle store made
// unconditional lose 12–95 % of them per route, in every run measured, on macOS/arm64, Debian/x86-64 (gcc)
// and Windows/x64 (MSVC) — two review rounds taught this gate to use a fresh head each round (with two
// alternating ones a loss after a loss read "right") and to publish only once the restarts are running
// (without that some routes lost nothing on Debian). The fixed bodies lose none, on any of the three, and
// 2000 publications a route keeps the suite under a second on Windows.
//
// Outside the wasm-audio tier by construction, like core's RtStreams suite: that tier is single-threaded by
// contract and links no pthreads.

#include <felitronics_test.h>
#include <felitronics/convolution/ConvolutionEngine.h>
#include <felitronics/convolution/MatrixConvolver.h>
#include <felitronics/core/Math.h>

#include <atomic>
#include <cstdio>
#include <string>
#include <thread>
#include <type_traits>

using namespace felitronics;

namespace
{
enum class Via { setIr, setOperator, stageThenPublish };

// One publication of a one-tap operator whose only tap is `h`, by the route under test.
template <class C> bool publish (C& c, Via via, float h)
{
    if constexpr (std::is_same_v<C, convolution::ConvolutionEngine<core::fft::DefaultRealFft, 1>>)
    {
        (void) via;
        return c.setIr (&h, 1);
    }
    else
    {
        const float* b[1] { &h };
        if (via == Via::setOperator) return c.setOperator (C::Topology::LRDiag, b, 1, 1);
        if (! c.stageOperator (C::Topology::LRDiag, b, 1, 1)) return false;
        c.publishStaged();
        return true;
    }
}

template <class C>
long lostPublications (const std::string& name, Via via, long iterations)
{
    C c;
    felitronics::test::ok (c.prepare (4, 4, 1, 1), name + ": prepare");
    // A DIFFERENT head for every publication: with two alternating values, a loss that follows a loss would
    // read the "right" head by chance and the count could never pass half (a review round caught that).
    const auto headOf = [] (long i) { return 0.25f + 1.0e-4f * (float) i; };
    std::atomic<int>  phase { 0 };            // 0 parked · 1 racing · 2 audio thread done for this round
    std::atomic<bool> loaderDone { false }, quit { false }, spinning { false };
    std::atomic<long> which { 0 };
    std::atomic<long> refused { 0 };

    std::thread audio ([&] {
        while (! quit.load())
        {
            if (phase.load() != 1) { std::this_thread::yield(); continue; }
            c.reset();
            spinning.store (true);             // the loader publishes only once restarts are really running
            while (! loaderDone.load()) c.reset();
            c.reset(); c.reset();              // an accepted publication must survive, and be adopted by, these
            phase.store (2);
        }
    });
    std::thread loader ([&] {
        unsigned spins = 1u;
        while (! quit.load())
        {
            if (phase.load() != 1 || loaderDone.load() || ! spinning.load()) { std::this_thread::yield(); continue; }
            spins = spins * 1103515245u + 12345u;                     // land the publish at varying points
            for (int k = (int) ((spins >> 16) % 51u); k > 0; --k) (void) c.isBusy();
            const float h = headOf (which.load());
            const bool accepted = publish (c, via, h);
            if (! accepted) refused.fetch_add (1);
            loaderDone.store (true);
        }
    });

    long lost = 0;
    for (long i = 0; i < iterations; ++i)
    {
        which.store (i);
        loaderDone.store (false);
        spinning.store (false);
        phase.store (1);
        while (phase.load() != 2) std::this_thread::yield();
        float in[2] { 1.0f, 1.0f }, out[2] { 0.0f, 0.0f };          // both threads parked: the probe is synchronised
        const float* ins[1] { in };
        float*       outs[1] { out };
        felitronics::test::run (c.process (ins, outs, 1, 2));
        if (! core::sameBits (out[1], headOf (i))) ++lost;
        phase.store (0);
    }
    quit.store (true);
    audio.join();
    loader.join();
    std::printf ("      %s: %ld publications, %ld lost, %ld refused\n", name.c_str(), iterations, lost, refused.load());
    return lost + refused.load();
}

} // namespace

int main()
{
    std::printf ("felitronics::convolution restart race (P88)\n");
    const long n = 2000;                         // see the header: overwhelming for the mutant, cheap on Windows
    felitronics::test::group ("reset() spinning beside a publishing loader never loses an accepted publication");
    {
        using Eng = convolution::ConvolutionEngine<core::fft::DefaultRealFft, 1>;
        using Mat = convolution::MatrixConvolver<>;
        const long e  = lostPublications<Eng> ("ConvolutionEngine, setIr", Via::setIr, n);
        const long m1 = lostPublications<Mat> ("MatrixConvolver, setOperator", Via::setOperator, n);
        const long m2 = lostPublications<Mat> ("MatrixConvolver, stage + publish", Via::stageThenPublish, n);
        felitronics::test::ok (e == 0,  "ConvolutionEngine: none lost, none refused (" + std::to_string (e) + ")");
        felitronics::test::ok (m1 == 0, "MatrixConvolver via setOperator: none lost, none refused (" + std::to_string (m1) + ")");
        felitronics::test::ok (m2 == 0, "MatrixConvolver via stage + publish: none lost, none refused (" + std::to_string (m2) + ")");
    }
    return felitronics::test::report();
}

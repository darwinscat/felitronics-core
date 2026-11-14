// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// MatrixConvolverNupc with REAL threads, through the spread B_max stage: a loader publishes operators of random
// length and topology — growing a slot's banks, shrinking, switching topology — while the audio thread processes
// random blocks and now and then clears its history or resets. The hazards are the ones a single thread cannot reach: the
// loader growing (reallocating) the INACTIVE slot's spectra while the audio thread reads the live one; a spread
// operator's preparation reading the incoming slot across many calls; the fade's end releasing a slot the loader
// may regrow at once. The arithmetic is checked elsewhere; here every call must be accepted, every sample finite,
// and the loader must get its publications through. Under the ThreadSanitizer job (FELITRONICS_ENABLE_TSAN) a data
// race fails it outright.
//
// The threshold is lowered (the class's template parameter) so short operators spread, and the cap is small, so a
// run is quick and passes through many swaps. Outside the wasm-audio tier by construction: that tier is
// single-threaded by contract.

#include <felitronics_test.h>
#include <felitronics/convolution/MatrixConvolverNupc.h>

#include <atomic>
#include <cmath>
#include <cstdio>
#include <string>
#include <thread>
#include <vector>

using namespace felitronics;

namespace
{
struct Lcg
{
    unsigned long long s;
    unsigned next() { s = s * 6364136223846793005ULL + 1442695040888963407ULL; return (unsigned) (s >> 33); }
    float bipolar() { return (float) (next() & 0xffffu) / 32768.0f - 1.0f; }
};

using Spread = convolution::MatrixConvolverNupc<core::fft::DefaultRealFft, 1>;

struct Outcome { long published = 0, refused = 0, calls = 0, rejectedCalls = 0, clears = 0, resets = 0; bool finite = true; };

Outcome run (int nch, int publications, unsigned long long seed)
{
    constexpr int cap = 8 * 2048 + 2048, xfade = 300;
    Spread mc;
    Outcome o;
    if (! mc.prepare (128, cap, xfade, nch)) { o.rejectedCalls = -1; return o; }

    // The operators the loader picks from: lengths from one head to past the cap, every topology.
    const int lengths[] { 100, 3000, 4097, 9000, cap / 2, cap - 1, cap, cap + 700 };
    std::vector<float> taps[4];
    Lcg g { seed };
    for (auto& t : taps) { t.resize ((std::size_t) cap + 700); for (float& v : t) v = 0.05f * g.bipolar(); }
    const float* banks[4] { taps[0].data(), taps[1].data(), taps[2].data(), taps[3].data() };

    std::atomic<bool> done { false };
    std::thread loader ([&]
    {
        Lcg r { seed * 7ULL + 1ULL };
        const Spread::Topology topos[] { Spread::Topology::LRDiag, Spread::Topology::Full, Spread::Topology::MSDiag };
        while (o.published < publications)
        {
            const int len = lengths[r.next() % 8u];
            const Spread::Topology t = topos[r.next() % 3u];
            if (mc.setOperator (t, banks, 4, len)) ++o.published;   // grows the inactive slot when it reaches further
            else { ++o.refused; std::this_thread::yield(); }        // a swap in flight: try again
        }
        done.store (true, std::memory_order_release);
    });

    // The audio thread: random blocks, a clear now and then.
    Lcg r { seed * 13ULL + 5ULL };
    std::vector<float> x[2], y[2];
    for (int c = 0; c < 2; ++c) { x[c].assign (700, 0.0f); y[c].assign (700, 0.0f); }
    while (! done.load (std::memory_order_acquire))
    {
        const int n = 1 + (int) (r.next() % 600u);
        for (int c = 0; c < nch; ++c) for (int i = 0; i < n; ++i) x[c][(std::size_t) i] = 0.3f * r.bipolar();
        const float* in[2] { x[0].data(), x[1].data() };
        float* out[2] { y[0].data(), y[1].data() };
        if (! mc.process (in, out, nch, n)) ++o.rejectedCalls;
        ++o.calls;
        for (int c = 0; c < nch; ++c) for (int i = 0; i < n; ++i) if (! std::isfinite (y[c][(std::size_t) i])) o.finite = false;
        const unsigned roll = r.next() % 97u;
        if (roll == 0u) { mc.clearAudioState(); ++o.clears; }   // the audio thread's own restarts: a clear keeps a
        if (roll == 1u) { mc.reset(); ++o.resets; }             // swap running; a reset adopts one (and reads its banks)
    }
    loader.join();
    return o;
}
} // namespace

int main()
{
    std::printf ("felitronics::convolution MatrixConvolverNupc — a loader beside the audio thread, through the spread stage\n");
    for (int nch = 1; nch <= 2; ++nch)
    {
        const Outcome o = run (nch, 400, 0x5eedULL + (unsigned long long) nch);
        const std::string tag = std::to_string (nch) + " ch: ";
        std::printf ("      %s%ld published, %ld refused while busy, %ld calls, %ld clears, %ld resets\n", tag.c_str(), o.published, o.refused, o.calls, o.clears, o.resets);
        test::ok (o.published == 400, tag + "every publication got through (" + std::to_string (o.published) + ")");
        test::ok (o.rejectedCalls == 0, tag + "every process() call accepted (" + std::to_string (o.rejectedCalls) + " refused)");
        test::ok (o.finite, tag + "every sample finite");
        test::ok (o.calls > 400, tag + "the audio thread ran beside the loader (" + std::to_string (o.calls) + " calls)");
    }
    return test::report();
}

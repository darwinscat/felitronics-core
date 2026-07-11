// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// JUCE-free self-tests for core::DryAligner — the PDC-alignment gate ported from OrbitCab's
// PowerAmpRouterAlignTests. The aligned dry must be delayed by EXACTLY the requested tap: bit-exact
// by construction (a pure copy through the ring — no arithmetic, no tolerance). Designed to BREAK on
// any off-by-one in the ring, a bad ring-wrap at ragged block sizes, a stale tap after a per-block
// tap change, or channel cross-talk in the shared-cursor multi-channel ring.

#include <felitronics_test.h>
#include <felitronics/core/DryAligner.h>

#include <atomic>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

static std::atomic<long> g_allocs { 0 };
void* operator new      (std::size_t s) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s ? s : 1); }
void* operator new[]    (std::size_t s) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s ? s : 1); }
void  operator delete   (void* p) noexcept { std::free (p); }
void  operator delete[] (void* p) noexcept { std::free (p); }
void  operator delete   (void* p, std::size_t) noexcept { std::free (p); }
void  operator delete[] (void* p, std::size_t) noexcept { std::free (p); }

using namespace felitronics;

namespace
{
    // Distinct-per-sample signal so a delay is detectable sample-exact (no value repeats within the
    // window). Pure LCG → reproducible, seedable per channel.
    std::vector<float> distinctSignal (int n, std::uint32_t seed = 0x1234567u)
    {
        std::vector<float> v ((std::size_t) n);
        std::uint32_t s = seed;
        for (int i = 0; i < n; ++i)
        {
            s = s * 1664525u + 1013904223u;
            v[(std::size_t) i] = (float) ((double) s / 4294967296.0 - 0.5);
        }
        return v;
    }

    // Run `in` through a DryAligner (1 ch) at fixed delay D in `block`-sized chunks; return the output.
    std::vector<float> runAligner (core::DryAligner& a, const std::vector<float>& in, int D, int block)
    {
        std::vector<float> out (in.size(), 0.0f), buf ((std::size_t) block);
        for (std::size_t off = 0; off < in.size(); off += (std::size_t) block)
        {
            const int n = (int) std::min ((std::size_t) block, in.size() - off);
            for (int i = 0; i < n; ++i) buf[(std::size_t) i] = in[off + (std::size_t) i];
            const float* io[1] = { buf.data() };
            a.advance (io, 1, n, D);
            for (int i = 0; i < n; ++i) out[off + (std::size_t) i] = a.delayed (0)[i];
        }
        return out;
    }

    // True iff `out[n] == in[n - d]` for every n >= d (bit-exact; the warm-up region n < d is ignored).
    bool isDelayedBy (const std::vector<float>& out, const std::vector<float>& in, int d)
    {
        for (std::size_t n = (std::size_t) d; n < in.size(); ++n)
            if (out[n] != in[n - (std::size_t) d]) return false;
        return true;
    }
}

int main()
{
    std::printf ("felitronics::core DryAligner tests\n");

    // --- exact delay at arbitrary taps across ragged block sizes (bit-exact) ---
    test::group ("exact delay at arbitrary taps across ragged block sizes (bit-exact)");
    {
        const auto in = distinctSignal (7000);
        for (int D : { 0, 1, 5, 31, 100, 255 })
            for (int blk : { 1, 7, 64, 128, 300, 512 })   // incl. blocks LARGER than the 256 ring capacity
            {
                core::DryAligner a; a.prepare (1, blk, 256);
                const auto out = runAligner (a, in, D, blk);
                test::ok (isDelayedBy (out, in, D),
                          "exact delay D=" + std::to_string (D) + " at block size " + std::to_string (blk));
            }
    }

    // --- zero-tap passthrough: D == 0 is an exact identity copy, sample 0 included (no warm-up) ---
    test::group ("zero-tap passthrough is a bit-identical copy");
    {
        const auto in = distinctSignal (3000, 0xBEEFu);
        core::DryAligner a; a.prepare (1, 128, 256);
        const auto out = runAligner (a, in, 0, 128);
        bool identical = out.size() == in.size();
        for (std::size_t n = 0; identical && n < in.size(); ++n) identical = (out[n] == in[n]);
        test::ok (identical, "D=0 output == input bit-for-bit from sample 0");
    }

    // --- a tap beyond capacity clamps (no overrun / crash) ---
    test::group ("a tap beyond capacity clamps to capacity-1");
    {
        const auto in = distinctSignal (3000);
        core::DryAligner a; a.prepare (1, 512, 256);
        const auto out = runAligner (a, in, 1000, 64);      // capacity 256 → clamps to 255
        test::ok (isDelayedBy (out, in, 255), "over-capacity tap clamps to capacity-1");
        test::ok (a.capacity() == 256, "capacity() reports the prepared ring size");
        const auto neg = runAligner (a, in, -7, 64);        // negative tap clamps to 0 (ring is warm by now,
        bool id = true;                                     //   D=0 reads back the freshest write regardless)
        for (std::size_t n = 0; id && n < in.size(); ++n) id = (neg[n] == in[n]);
        test::ok (id, "negative tap clamps to 0 (exact identity)");
    }

    // --- variable tap per block: the ring is fed EVERY block, so the tap may move block-to-block ---
    test::group ("variable tap per block stays sample-exact (warm ring, moving tap)");
    {
        // With a continuously fed ring, a block using tap D must emit in[n - D] for every global index
        // n >= D in that block — regardless of what tap earlier blocks used. This is the tube↔capture /
        // model-load-moves-latency case: no stale-tap, no re-warm glitch.
        const int blk = 64, N = 6400;
        const auto in = distinctSignal (N, 0xCAFEu);
        const int taps[] = { 0, 17, 96, 255, 3, 128, 1, 200 };   // moves up AND down, incl. the extremes
        core::DryAligner a; a.prepare (1, blk, 256);
        std::vector<float> buf ((std::size_t) blk);
        bool exact = true;
        for (int off = 0; off < N; off += blk)
        {
            const int D = taps[(std::size_t) (off / blk) % (sizeof (taps) / sizeof (taps[0]))];
            for (int i = 0; i < blk; ++i) buf[(std::size_t) i] = in[(std::size_t) (off + i)];
            const float* io[1] = { buf.data() };
            a.advance (io, 1, blk, D);
            for (int i = 0; i < blk && exact; ++i)
            {
                const int n = off + i;
                if (n >= D) exact = (a.delayed (0)[i] == in[(std::size_t) (n - D)]);
            }
            if (! exact)
            {
                test::ok (false, "block at offset " + std::to_string (off) + " (tap " + std::to_string (D) + ") broke alignment");
                break;
            }
        }
        test::ok (exact, "every block bit-exact under a per-block moving tap");
    }

    // --- multi-channel independence: one shared cursor, per-channel history ---
    test::group ("multi-channel: channels delay independently (no cross-talk)");
    {
        const int N = 5000, blk = 128, D = 77;
        const auto in0 = distinctSignal (N, 0x1111u);
        const auto in1 = distinctSignal (N, 0x2222u);
        const auto in2 = distinctSignal (N, 0x3333u);
        core::DryAligner a; a.prepare (3, blk, 256);
        std::vector<float> b0 ((std::size_t) blk), b1 ((std::size_t) blk), b2 ((std::size_t) blk);
        std::vector<float> o0 ((std::size_t) N), o1 ((std::size_t) N), o2 ((std::size_t) N);
        for (int off = 0; off < N; off += blk)
        {
            const int n = std::min (blk, N - off);
            for (int i = 0; i < n; ++i)
            {
                b0[(std::size_t) i] = in0[(std::size_t) (off + i)];
                b1[(std::size_t) i] = in1[(std::size_t) (off + i)];
                b2[(std::size_t) i] = in2[(std::size_t) (off + i)];
            }
            const float* io[3] = { b0.data(), b1.data(), b2.data() };
            a.advance (io, 3, n, D);
            for (int i = 0; i < n; ++i)
            {
                o0[(std::size_t) (off + i)] = a.delayed (0)[i];
                o1[(std::size_t) (off + i)] = a.delayed (1)[i];
                o2[(std::size_t) (off + i)] = a.delayed (2)[i];
            }
        }
        test::ok (isDelayedBy (o0, in0, D), "channel 0 delayed by exactly D");
        test::ok (isDelayedBy (o1, in1, D), "channel 1 delayed by exactly D");
        test::ok (isDelayedBy (o2, in2, D), "channel 2 delayed by exactly D");
    }

    // --- the staged copy survives the in-place stage overwriting the live buffer ---
    test::group ("delayed() is a real copy — survives the caller's buffer being overwritten");
    {
        const int blk = 32;
        const auto in = distinctSignal (blk, 0x7777u);
        core::DryAligner a; a.prepare (1, blk, 64);
        std::vector<float> live (in);
        const float* io[1] = { live.data() };
        a.advance (io, 1, blk, 0);
        std::fill (live.begin(), live.end(), -99.0f);       // the in-place stage trashes the live buffer
        bool intact = true;
        for (int i = 0; i < blk && intact; ++i) intact = (a.delayed (0)[i] == in[(std::size_t) i]);
        test::ok (intact, "scratch holds the pre-stage input after the live buffer changed");
    }

    // --- reset(): history cleared, warm-up region reads silence ---
    test::group ("reset clears the ring (warm-up reads zeros)");
    {
        const int blk = 64, D = 50;
        const auto in = distinctSignal (1000, 0x4444u);
        core::DryAligner a; a.prepare (1, blk, 256);
        (void) runAligner (a, in, D, blk);                  // warm the ring
        a.reset();
        const auto out = runAligner (a, in, D, blk);
        bool warmupSilent = true;
        for (int i = 0; i < D && warmupSilent; ++i) warmupSilent = (out[(std::size_t) i] == 0.0f);
        test::ok (warmupSilent, "first D samples after reset are exact zeros (no stale history)");
        test::ok (isDelayedBy (out, in, D), "post-reset stream realigns bit-exactly");
    }

    // --- 🔴 RT: advance() performs zero heap allocations ---
    test::group ("advance() no-alloc");
    {
        core::DryAligner a; a.prepare (2, 512, 4096);
        std::vector<float> l (512, 0.25f), r (512, -0.5f);
        const float* io[2] = { l.data(), r.data() };
        a.advance (io, 2, 512, 100);                        // warm call
        const long before = g_allocs.load();
        a.advance (io, 2, 512, 100);
        a.advance (io, 2, 512, 3000);                       // tap jump — still no realloc
        a.reset();
        a.advance (io, 2, 512, 0);
        test::okNoAlloc (g_allocs.load() == before, "advance()/reset() did not allocate");
    }

    return test::report();
}

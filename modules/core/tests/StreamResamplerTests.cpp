// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// Streaming Catmull-Rom rate-matcher (promoted from OrbitCab, where it rate-locks the NAM stage to
// 48 kHz on any host rate). The first four groups are OrbitCab's AmpStageTests resampler cases,
// ported with the SAME expected values — they are the golden this extraction must not move:
// identity ratio = a clean 2-sample delay, up/down conversion keeps level and output count, DC
// shows no drift over many blocks. On top: block-size invariance (same count, samples within one
// float ulp — the phase accumulator runs at different magnitudes per chunking, measured 5.96e-8
// worst-case; bit-exact at ratio 1), startup underflow padding, and the no-alloc-in-process proof.

#include <felitronics/core/StreamResampler.h>

#include "felitronics_test.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdlib>
#include <vector>

// global allocation counter (no-alloc-in-process proof)
static std::atomic<long> g_allocs { 0 };
void* operator new      (std::size_t s) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s ? s : 1); }
void* operator new[]    (std::size_t s) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s ? s : 1); }
void  operator delete   (void* p) noexcept { std::free (p); }
void  operator delete[] (void* p) noexcept { std::free (p); }
void  operator delete   (void* p, std::size_t) noexcept { std::free (p); }
void  operator delete[] (void* p, std::size_t) noexcept { std::free (p); }

using felitronics::core::StreamResampler;
using felitronics::test::approx;
using felitronics::test::group;
using felitronics::test::ok;

namespace
{
    constexpr double kPi = 3.14159265358979323846;

    // Run a signal through a resampler in `block`-sized chunks; return the concatenated output.
    std::vector<float> runResampler (double inRate, double outRate, const std::vector<float>& in, int block)
    {
        StreamResampler r;
        r.reset (inRate, outRate, block * 2 + 16);
        std::vector<float> out, tmp ((std::size_t) (block * 4 + 16));
        out.reserve ((std::size_t) ((double) in.size() * outRate / inRate) + 32);
        for (std::size_t i = 0; i < in.size(); i += (std::size_t) block)
        {
            const int n = (int) std::min ((std::size_t) block, in.size() - i);
            r.feed (in.data() + i, n);
            const int got = r.produceAvailable (tmp.data(), (int) tmp.size());
            out.insert (out.end(), tmp.begin(), tmp.begin() + got);
        }
        return out;
    }

    bool anyBad (const std::vector<float>& v)
    { for (float x : v) if (std::isnan (x) || std::isinf (x)) return true; return false; }

    float peak (const std::vector<float>& v)
    { float p = 0; for (float x : v) p = std::max (p, std::fabs (x)); return p; }

    std::vector<float> sine (int n, double rate, double f, float amp = 0.5f)
    {
        std::vector<float> v ((std::size_t) n);
        for (int i = 0; i < n; ++i) v[(std::size_t) i] = amp * (float) std::sin (2.0 * kPi * f * i / rate);
        return v;
    }
}

int main()
{
    std::printf ("felitronics::core stream-resampler tests\n");

    group ("identity at ratio 1 — clean 2-sample delay (catmull @ t=0)");
    {
        auto in  = sine (4000, 48000.0, 600.0);
        auto out = runResampler (48000.0, 48000.0, in, 512);
        ok (! anyBad (out), "no NaN/Inf at identity ratio");
        int matched = 0, checked = 0;
        for (std::size_t k = 2; k + 2 < out.size() && k < in.size() && k < 2000; ++k, ++checked)
            if (std::fabs (out[k] - in[k - 2]) < 1.0e-4f) ++matched;
        ok (checked > 1000 && matched > checked - 4, "out[k] == in[k-2] over the checked run");
    }

    group ("upsample 44100 -> 48000 (no NaN, level kept, count ~ ratio)");
    {
        auto in  = sine (44100, 44100.0, 100.0);
        auto out = runResampler (44100.0, 48000.0, in, 512);
        ok (! anyBad (out), "no NaN/Inf on upsample");
        approx (peak (out), 0.5, 0.04, "sine level survives the upsample");
        approx ((double) out.size() / (double) in.size(), 48000.0 / 44100.0, 0.01, "output count tracks the ratio");
    }

    group ("downsample 96000 -> 48000 (no NaN, level kept, count ~ half)");
    {
        auto in  = sine (96000, 96000.0, 100.0);
        auto out = runResampler (96000.0, 48000.0, in, 512);
        ok (! anyBad (out), "no NaN/Inf on downsample");
        approx (peak (out), 0.5, 0.04, "sine level survives the downsample");
        approx ((double) out.size() / (double) in.size(), 0.5, 0.01, "output count is ~ half");
    }

    group ("no drift on DC over many blocks");
    {
        std::vector<float> in (48000, 1.0f);
        auto out = runResampler (44100.0, 48000.0, in, 512);
        ok (! anyBad (out), "no NaN/Inf on DC");
        double tail = 0; int c = 0;
        for (std::size_t k = (out.size() > 256 ? out.size() - 256 : 0); k < out.size(); ++k) { tail += out[k]; ++c; }
        approx (c ? tail / c : 0.0, 1.0, 0.02, "DC tail still sits at 1 after a second of resampling");
    }

    group ("block-size invariance — chunking never changes the stream");
    {
        auto in = sine (9000, 44100.0, 440.0);
        const auto a = runResampler (44100.0, 48000.0, in, 512);
        const auto b = runResampler (44100.0, 48000.0, in, 64);
        const auto c = runResampler (44100.0, 48000.0, in, 48);   // not a divisor of the input length
        ok (a.size() == b.size() && a.size() == c.size(), "output count is independent of the feed block size");
        float worst = 0.0f;
        const std::size_t n = std::min ({ a.size(), b.size(), c.size() });
        for (std::size_t k = 0; k < n; ++k)
            worst = std::max ({ worst, std::fabs (a[k] - b[k]), std::fabs (a[k] - c[k]) });
        ok (worst <= 1.5e-7f, "samples match within ~1 float ulp across chunkings");   // measured 5.96e-8 worst-case

        const auto i1 = runResampler (48000.0, 48000.0, in, 512);
        const auto i2 = runResampler (48000.0, 48000.0, in, 48);
        ok (i1 == i2, "identity ratio is BIT-exact across chunkings (integer phase)");
    }

    group ("produceExact pads startup underflow with silence");
    {
        StreamResampler r;
        r.reset (48000.0, 48000.0, 64);
        float out[8] = { 9, 9, 9, 9, 9, 9, 9, 9 };
        r.produceExact (out, 8);                                 // nothing fed yet — only the 3 history zeros
        bool allZero = true;
        for (float x : out) allZero = allZero && x == 0.0f;
        ok (allZero, "unfed resampler emits exact silence, never junk");
    }

    group ("RT safety — feed/produce never allocate");
    {
        StreamResampler r;
        r.reset (44100.0, 48000.0, 512 * 2 + 16);                // ALLOC lives here, message thread
        auto in = sine (512 * 8, 44100.0, 220.0);
        std::vector<float> tmp ((std::size_t) (512 * 4 + 16));
        const long before = g_allocs.load();
        for (int blk = 0; blk < 8; ++blk)
        {
            r.feed (in.data() + (std::size_t) (blk * 512), 512);
            (void) r.produceAvailable (tmp.data(), (int) tmp.size());
        }
        r.produceExact (tmp.data(), 32);
        felitronics::test::okNoAlloc (g_allocs.load() == before, "feed/produceAvailable/produceExact allocated nothing");
    }

    return felitronics::test::report();
}

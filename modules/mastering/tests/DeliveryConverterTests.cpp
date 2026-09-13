// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// P57 step 2 — mastering::DeliveryConverter: SRC first, as one whole-programme operation. Every oracle
// here is independent of the arithmetic it checks: the length against a second integer expression and
// the double trap it replaces, the trim and drain against the same programme followed by real silence,
// the content against the analytic signal at the stated sub-sample offset.

#include <felitronics/mastering/DeliveryConverter.h>
#include <felitronics_test.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <new>
#include <string>
#include <vector>

static std::atomic<long long> g_allocs { 0 }, g_allocBytes { 0 };
void* operator new      (std::size_t s) { g_allocs.fetch_add (1); g_allocBytes.fetch_add ((long long) s); return std::malloc (s ? s : 1); }
void* operator new[]    (std::size_t s) { g_allocs.fetch_add (1); g_allocBytes.fetch_add ((long long) s); return std::malloc (s ? s : 1); }
void  operator delete   (void* p) noexcept { std::free (p); }
void  operator delete[] (void* p) noexcept { std::free (p); }
void  operator delete   (void* p, std::size_t) noexcept { std::free (p); }
void  operator delete[] (void* p, std::size_t) noexcept { std::free (p); }

using felitronics::mastering::DeliveryConverter;
using felitronics::test::group;
using felitronics::test::ok;
using felitronics::test::okNoAlloc;

namespace
{
    constexpr double kPi = 3.14159265358979323846;
    constexpr double kRates[] = { 44100.0, 48000.0, 88200.0, 96000.0, 176400.0, 192000.0 };

    struct Planar
    {
        std::vector<std::vector<float>> ch;
        std::vector<float*> ptr;
        Planar (int n, long long frames) : ch ((std::size_t) n, std::vector<float> ((std::size_t) frames, 0.0f)), ptr ((std::size_t) n)
        { for (int c = 0; c < n; ++c) ptr[(std::size_t) c] = ch[(std::size_t) c].data(); }
    };

    bool convertWith (double a, double b, const Planar& in, long long inFrames, Planar& out, long long outFrames, int block)
    {
        DeliveryConverter dc;
        if (! dc.prepare (a, b, (int) in.ch.size(), block)) return false;
        std::vector<const float*> ip (in.ch.size());
        for (std::size_t c = 0; c < in.ch.size(); ++c) ip[c] = in.ch[c].data();
        return dc.convert (ip.data(), (int) in.ch.size(), inFrames, out.ptr.data(), outFrames);
    }
}

static void testLength()
{
    group ("the delivered length is exact integer arithmetic, not a ceil of a double");
    ok (DeliveryConverter::deliveredFrames (44100.0, 48000.0, 147) == 160,
        "147 frames at 44.1 -> 48 is 160 — std::ceil of the double quotient says 161");
    ok (DeliveryConverter::deliveredFrames (176400.0, 44100.0, 176401) == 44101, "176,401 at 176.4 -> 44.1 is 44,101");
    ok (DeliveryConverter::deliveredFrames (48000.0, 48000.0, 12345) == 12345, "equal rates deliver every frame");
    ok (DeliveryConverter::deliveredFrames (48000.0, 44100.0, 0) == 0, "an empty programme delivers nothing");
    ok (DeliveryConverter::deliveredFrames (44100.5, 48000.0, 100) == -1, "a rate the resampler refuses has no length");
    bool match = true;
    long long doubleWrong = 0;
    for (double a : kRates)
        for (double b : kRates)
        {
            const long long g = std::gcd ((long long) a, (long long) b), L = (long long) b / g, M = (long long) a / g;
            for (long long n = 1; n < 3000; n += (n < 400 ? 1 : 37))
            {
                const long long want = (n * L - 1) / M + 1;            // a second expression of ceil(nL/M)
                if (DeliveryConverter::deliveredFrames (a, b, n) != want) match = false;
                // "frames x ratio", the way it is naturally written — the quotient first. (`n * b / a`
                // happens to be exact on these rates, which is why the trap is easy to miss by testing.)
                if ((long long) std::ceil ((double) n * (b / a)) != want) ++doubleWrong;
            }
        }
    ok (match, "deliveredFrames == floor((n*L - 1)/M) + 1 on every pair and length tried");
    ok (doubleWrong > 0, "…and the double expression it replaces really is wrong on some of them");
    std::printf ("      the double ceil disagrees on %lld of the lengths tried\n", doubleWrong);
}

static void testConvert()
{
    group ("the trim and the drain: convert(P) is the head of convert(P followed by real silence)");
    for (double a : kRates)
        for (double b : kRates)
        {
            if (a == b) continue;
            const long long N = 3001, K = 4096;
            Planar p (2, N + K);
            for (long long i = 0; i < N; ++i)                  // audio up to the LAST frame, then silence
            {
                p.ch[0][(std::size_t) i] = (float) (0.7 * std::sin (2.0 * kPi * 997.0 * (double) i / a));
                p.ch[1][(std::size_t) i] = (float) (0.5 * std::sin (2.0 * kPi * 3001.0 * (double) i / a + 1.0));
            }
            const long long outN = DeliveryConverter::deliveredFrames (a, b, N);
            const long long outNK = DeliveryConverter::deliveredFrames (a, b, N + K);
            Planar o1 (2, outN), o2 (2, outNK);
            ok (convertWith (a, b, p, N, o1, outN, 512) && convertWith (a, b, p, N + K, o2, outNK, 512), "both conversions accepted");
            bool head = true;
            for (int c = 0; c < 2; ++c)
                if (std::memcmp (o1.ch[(std::size_t) c].data(), o2.ch[(std::size_t) c].data(), (std::size_t) outN * sizeof (float)) != 0) head = false;
            ok (head, "the drain delivers exactly what real silence after the programme would have");
        }

    group ("the content: the delivered tone is the analytic tone at the stated sub-sample offset");
    for (const auto& pr : { std::pair<double, double> { 48000.0, 44100.0 }, { 192000.0, 44100.0 }, { 44100.0, 192000.0 }, { 176400.0, 48000.0 } })
    {
        const double a = pr.first, b = pr.second, f = 1000.0;
        const long long N = (long long) a * 2;
        Planar p (1, N);
        for (long long i = 0; i < N; ++i) p.ch[0][(std::size_t) i] = (float) (0.8 * std::sin (2.0 * kPi * f * (double) i / a));
        const long long outN = DeliveryConverter::deliveredFrames (a, b, N);
        Planar o (1, outN);
        DeliveryConverter dc;
        ok (dc.prepare (a, b, 1, 1024), "prepare");
        const float* ip[1] = { p.ch[0].data() };
        ok (dc.convert (ip, 1, N, o.ptr.data(), outN), "convert");
        const double d = dc.latencyOutputSamples() - (double) dc.trimSamples();
        ok (std::fabs (d) <= 0.5, "the residual offset is within half a sample");
        double worst = 0.0;
        for (long long n = (long long) b / 4; n < outN - (long long) b / 4; ++n)
            worst = std::max (worst, std::fabs ((double) o.ch[0][(std::size_t) n] - 0.8 * std::sin (2.0 * kPi * f * ((double) n - d) / b)));
        ok (20.0 * std::log10 (std::max (1e-300, worst)) <= -120.0, "the delivered samples null against x((n - d)/fs) to -120 dB");
    }

    group ("block size is free, equal rates copy the bits, and a wrong length is refused");
    {
        const long long N = 5000;
        Planar p (2, N);
        for (long long i = 0; i < N; ++i) { p.ch[0][(std::size_t) i] = (float) std::sin (0.01 * (double) i); p.ch[1][(std::size_t) i] = (float) std::cos (0.02 * (double) i); }
        const long long outN = DeliveryConverter::deliveredFrames (96000.0, 44100.0, N);
        Planar x1 (2, outN), x2 (2, outN), x3 (2, outN);
        ok (convertWith (96000.0, 44100.0, p, N, x1, outN, 1) && convertWith (96000.0, 44100.0, p, N, x2, outN, 777)
            && convertWith (96000.0, 44100.0, p, N, x3, outN, 65536), "three block sizes accepted");
        ok (x1.ch == x2.ch && x1.ch == x3.ch, "…and bit-identical");

        Planar y (2, N);
        ok (convertWith (48000.0, 48000.0, p, N, y, N, 256) && y.ch == p.ch, "equal rates: the caller's bits come back");
        Planar z (2, outN + 1);
        ok (! convertWith (96000.0, 44100.0, p, N, z, outN + 1, 256), "a length the formula did not give is refused");
    }

    group ("law 11d and RT: the budget is what prepare() asks for, and convert() asks for nothing");
    for (const auto& pr : { std::pair<double, double> { 48000.0, 44100.0 }, { 192000.0, 44100.0 }, { 44100.0, 176400.0 }, { 96000.0, 96000.0 } })
        for (int ch : { 1, 2 })
        {
            const std::uint64_t want = DeliveryConverter::prepareBytes (pr.first, pr.second, ch, 1024);
            DeliveryConverter dc;
            const long long before = g_allocBytes.load();
            const bool prepared = dc.prepare (pr.first, pr.second, ch, 1024);
            const long long asked = g_allocBytes.load() - before;
            ok (prepared && want > 0u, "prepared, with a budget");
            okNoAlloc ((std::uint64_t) asked == want, "prepareBytes == the bytes prepare() requested");

            const long long N = 20000, outN = DeliveryConverter::deliveredFrames (pr.first, pr.second, N);
            Planar p (ch, N), o (ch, outN);
            std::vector<const float*> ip ((std::size_t) ch);
            for (int c = 0; c < ch; ++c) ip[(std::size_t) c] = p.ch[(std::size_t) c].data();
            const long long allocs = g_allocs.load();
            const bool converted = dc.convert (ip.data(), ch, N, o.ptr.data(), outN);
            const long long after = g_allocs.load();
            ok (converted, "convert accepted");
            okNoAlloc (after == allocs, "convert() allocated nothing");
        }
}

// P64. `convert` checked its planes for null and nothing else, and it writes `out` at the delivery stride while still
// reading `in` at the source one. So an output plane over an input plane was accepted, and upsampling — where the
// writes run AHEAD of the reads — overwrote programme not yet read. Every layout here against the same conversion
// into buffers of its own, bit for bit; each of the two lengths used for its own side, in both directions of the
// ratio, because a rule judged at one length misses an overlap that lies only in the longer span.
static void testPlanes()
{
    group ("planes: an output over any input or another output is refused, each side at its own length");
    const long long N = 48000;
    for (const auto& pr : { std::pair<double, double> { 44100.0, 48000.0 }, { 48000.0, 44100.0 } })
    {
        const bool up = pr.second > pr.first;
        const long long D = DeliveryConverter::deliveredFrames (pr.first, pr.second, N);
        const auto Nz = (std::size_t) N, Dz = (std::size_t) D, M = std::max (Nz, Dz);
        Planar src (2, N);
        for (long long i = 0; i < N; ++i)
        {
            src.ch[0][(std::size_t) i] = (float) (0.3 * std::sin (0.01 * (double) i));
            src.ch[1][(std::size_t) i] = (float) (0.3 * std::sin (0.037 * (double) i));
        }
        Planar honest (2, D);
        ok (convertWith (pr.first, pr.second, src, N, honest, D, 4096), up ? "up: the honest conversion" : "down: the honest conversion");

        DeliveryConverter dc;
        if (! dc.prepare (pr.first, pr.second, 2, 4096)) { ok (false, "prepared"); continue; }
        const auto same = [] (const float* a, const float* b, std::size_t n) { return std::memcmp (a, b, n * sizeof (float)) == 0; };

        // (1) THE WITNESS: out[0] IS in[1], its buffer long enough for either length.
        {
            std::vector<float> A = src.ch[0], B (M, 0.25f), C (Dz, 0.25f);
            std::copy (src.ch[1].begin(), src.ch[1].end(), B.begin());
            const std::vector<float> B0 = B, C0 = C;
            const float* in[2] = { A.data(), B.data() };
            float* out[2] = { B.data(), C.data() };
            const bool accepted = dc.convert (in, 2, N, out, D);
            if (accepted)
            {
                long long wrong = 0;
                for (std::size_t i = 0; i < Dz; ++i) wrong += (C[i] != honest.ch[1][i]) ? 1 : 0;
                std::printf ("      %s: out[0] = in[1] was ACCEPTED; channel 1 differs from the honest conversion in %lld of %lld frames\n",
                             up ? "up" : "down", wrong, D);
            }
            ok (! accepted && B == B0 && C == C0, "out[0] = in[1] is refused, and nothing is written");
        }

        // (2) THE LONGER SIDE ONLY. Upsampling, the output is the longer plane: it starts N frames before an input
        //     plane and reaches into it. Downsampling, the input is: an output plane starts D frames into it. At the
        //     SHORTER length both would look disjoint.
        {
            std::vector<float> pool (Nz + Dz + M, 0.25f), other (M, 0.25f), spare (Dz, 0.25f);
            const float* in[2] {};
            float* out[2] {};
            if (up)
            {
                std::copy (src.ch[1].begin(), src.ch[1].end(), pool.begin() + (std::ptrdiff_t) Nz);
                std::copy (src.ch[0].begin(), src.ch[0].end(), other.begin());
                in[0] = other.data(); in[1] = pool.data() + Nz;
                out[0] = pool.data(); out[1] = spare.data();                 // [0, D) reaches into [N, 2N)
            }
            else
            {
                std::copy (src.ch[0].begin(), src.ch[0].end(), pool.begin());
                std::copy (src.ch[1].begin(), src.ch[1].end(), other.begin());
                in[0] = pool.data(); in[1] = other.data();
                out[0] = pool.data() + Nz + Dz; out[1] = pool.data() + Dz;   // [D, 2D) starts inside [0, N)
            }
            const std::vector<float> pool0 = pool;
            ok (! dc.convert (in, 2, N, out, D) && pool == pool0,
                up ? "up: an output plane reaching an input only past N frames is refused"
                   : "down: an output plane starting inside an input only past D frames is refused");
        }

        // (3) LEGAL: all four planes edge to edge in one allocation, each at its own length — which a rule judged at
        //     the LONGER length for both would refuse.
        {
            std::vector<float> pool (2 * Nz + 2 * Dz, 0.25f);
            std::copy (src.ch[0].begin(), src.ch[0].end(), pool.begin());
            std::copy (src.ch[1].begin(), src.ch[1].end(), pool.begin() + (std::ptrdiff_t) (Nz + Dz));
            const float* in[2] = { pool.data(), pool.data() + Nz + Dz };
            float* out[2] = { pool.data() + Nz, pool.data() + 2 * Nz + Dz };
            ok (dc.convert (in, 2, N, out, D) && same (out[0], honest.ch[0].data(), Dz) && same (out[1], honest.ch[1].data(), Dz),
                "planes edge to edge at their own lengths convert, bit-identical to buffers of their own");
        }

        // (3b) LEGAL: the two OUTPUT planes edge to edge, at the output length. Downsampling, the input length is the
        //      longer one, and outputs judged at it would overlap where they do not — a case the separate buffers
        //      above leave to wherever the allocator put them (the mutation stand).
        {
            std::vector<float> pool (2 * Dz, 0.25f);
            const float* in[2] = { src.ch[0].data(), src.ch[1].data() };
            float* out[2] = { pool.data(), pool.data() + Dz };
            ok (dc.convert (in, 2, N, out, D) && same (out[0], honest.ch[0].data(), Dz) && same (out[1], honest.ch[1].data(), Dz),
                "two output planes edge to edge convert, bit-identical to buffers of their own");
        }

        // (4) Two OUTPUT planes one frame apart: every frame the conversion writes lands on a sample of the other plane.
        {
            std::vector<float> pool (Dz + 1, 0.25f);
            const std::vector<float> pool0 = pool;
            const float* in[2] = { src.ch[0].data(), src.ch[1].data() };
            float* out[2] = { pool.data(), pool.data() + 1 };
            ok (! dc.convert (in, 2, N, out, D) && pool == pool0, "two output planes one frame apart are refused");
        }

        // (4b) Upsampling, two OUTPUT planes N frames apart: they share the last D − N frames, and only the output
        //      length sees it — outputs judged at the shorter of the two lengths passed every line above (the review
        //      round, by mutation).
        if (up)
        {
            std::vector<float> pool (Nz + Dz, 0.25f);
            const std::vector<float> pool0 = pool;
            const float* in[2] = { src.ch[0].data(), src.ch[1].data() };
            float* out[2] = { pool.data(), pool.data() + Nz };
            ok (! dc.convert (in, 2, N, out, D) && pool == pool0,
                "up: two output planes N frames apart, sharing only the frames past N, are refused");
        }
    }
}

// The count of non-finite input samples is THE LAST convert()'s THAT REACHED IT. A call refused before reading its
// input leaves the previous count exactly where it was — the rule `fc_solution_log` keeps for `written` — and the next
// call that reads its programme replaces it, an empty one with 0.
static void testCountIsTheLastCallThatReachedIt()
{
    group ("nonFiniteInputSamples: the last convert() that reached its count");
    const long long N = 4800;
    for (const auto& pr : { std::pair<double, double> { 44100.0, 48000.0 }, { 48000.0, 48000.0 } })
    {
        const long long D = DeliveryConverter::deliveredFrames (pr.first, pr.second, N);
        Planar dirty (2, N), clean (2, N), out (2, D);
        for (long long i = 0; i < N; ++i)
            for (int c = 0; c < 2; ++c)
                dirty.ch[(std::size_t) c][(std::size_t) i] = clean.ch[(std::size_t) c][(std::size_t) i] = (float) std::sin (0.01 * (double) (i + c));
        dirty.ch[0][10] = std::numeric_limits<float>::quiet_NaN();
        dirty.ch[1][N - 1] = std::numeric_limits<float>::infinity();
        dirty.ch[1][20] = std::numeric_limits<float>::quiet_NaN();
        DeliveryConverter dc;
        if (! dc.prepare (pr.first, pr.second, 2, 1024)) { ok (false, "prepared"); continue; }
        const std::string tag = pr.first == pr.second ? "equal rates: " : "44.1 -> 48: ";
        const float* di[2] = { dirty.ch[0].data(), dirty.ch[1].data() };
        const float* ci[2] = { clean.ch[0].data(), clean.ch[1].data() };
        ok (dc.convert (di, 2, N, out.ptr.data(), D) && dc.nonFiniteInputSamples() == 3u, tag + "a call that reaches its count: 3");
        ok (! dc.convert (ci, 2, N, out.ptr.data(), D + 1) && dc.nonFiniteInputSamples() == 3u, tag + "refused on its length: still 3");
        float* over[2] = { dirty.ch[1].data(), out.ptr[1] };                  // an output over an input
        ok (! dc.convert (di, 2, N, over, D) && dc.nonFiniteInputSamples() == 3u, tag + "refused on its planes: still 3");
        ok (dc.convert (ci, 2, N, out.ptr.data(), D) && dc.nonFiniteInputSamples() == 0u, tag + "the next that reaches it: 0");
        ok (dc.convert (di, 2, N, out.ptr.data(), D) && dc.nonFiniteInputSamples() == 3u
            && dc.convert (ci, 2, 0, out.ptr.data(), 0) && dc.nonFiniteInputSamples() == 0u,
            tag + "and an empty programme reaches it with nothing to read: 0");
    }
}

int main()
{
    std::printf ("felitronics::mastering::DeliveryConverter — SRC first, whole programme\n");
    testLength();
    testConvert();
    testPlanes();
    testCountIsTheLastCallThatReachedIt();
    return felitronics::test::report();
}

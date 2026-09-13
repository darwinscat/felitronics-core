// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// P57 step 2 — mastering::DeliveryConverter: SRC first, as one whole-programme operation. Every oracle
// here is independent of the arithmetic it checks: the length against a second integer expression and
// the double trap it replaces, the trim and drain against the same programme followed by real silence,
// the content against the analytic signal at the stated sub-sample offset.

#include <felitronics/mastering/DeliveryConverter.h>
#include <felitronics_test.h>

#include <atomic>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <new>
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

int main()
{
    std::printf ("felitronics::mastering::DeliveryConverter — SRC first, whole programme\n");
    testLength();
    testConvert();
    return felitronics::test::report();
}

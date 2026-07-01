// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// JUCE-free self-tests for the zero-latency partitioned convolver: equality with direct (time-domain)
// convolution, including under HOSTILE variable block splits + impulses at partition boundaries (the
// off-by-one alignment trap), the zero-latency assertion, and a no-allocation-in-process() proof.

#include <felitronics_test.h>
#include <felitronics/convolution/PartitionedConvolver.h>

#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

// --- global allocation counter (proves process() does not allocate) ---
static std::atomic<long> g_allocs { 0 };
void* operator new      (std::size_t s) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s ? s : 1); }
void* operator new[]    (std::size_t s) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s ? s : 1); }
void  operator delete   (void* p) noexcept { std::free (p); }
void  operator delete[] (void* p) noexcept { std::free (p); }
void  operator delete   (void* p, std::size_t) noexcept { std::free (p); }
void  operator delete[] (void* p, std::size_t) noexcept { std::free (p); }

using namespace felitronics;

// A tiny deterministic LCG so the test is reproducible (no <random> alloc surprises).
struct Lcg { unsigned long long s; float next() { s = s * 6364136223846793005ULL + 1442695040888963407ULL; return (float) ((s >> 40) & 0xffff) / 32768.0f - 1.0f; } };

static std::vector<float> directConv (const std::vector<float>& x, const std::vector<float>& h)
{
    std::vector<float> y (x.size(), 0.0f);
    for (int t = 0; t < (int) x.size(); ++t)
    {
        double acc = 0.0;
        for (int i = 0; i < (int) h.size() && i <= t; ++i) acc += (double) h[i] * x[t - i];
        y[(std::size_t) t] = (float) acc;
    }
    return y;
}

static double maxAbsDiff (const std::vector<float>& a, const std::vector<float>& b)
{
    double m = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) m = std::max (m, (double) std::fabs (a[i] - b[i]));
    return m;
}

int main()
{
    std::printf ("felitronics::convolution tests\n");

    // --- correctness vs direct convolution, processed as ONE block ---
    test::group ("Conv == direct (single block)");
    {
        const int P = 128, irLen = 500, n = 1000;
        Lcg r { 12345 };
        std::vector<float> ir (irLen), x (n);
        for (auto& v : ir) v = 0.2f * r.next();
        for (auto& v : x)  v = 0.5f * r.next();

        convolution::PartitionedConvolver<> conv;
        test::ok (conv.prepare (P, irLen), "prepare");
        conv.setIr (ir.data(), irLen);

        std::vector<float> y (n, 0.0f);
        conv.process (x.data(), y.data(), n);
        test::ok (maxAbsDiff (y, directConv (x, ir)) < 2e-3, "single-block output == direct conv");
    }

    // --- the alignment trap: HOSTILE variable block splits must still equal direct conv ---
    test::group ("Conv == direct (hostile variable block splits)");
    {
        const int P = 64, irLen = 300, n = 900;
        Lcg r { 99 };
        std::vector<float> ir (irLen), x (n), y (n, 0.0f);
        for (auto& v : ir) v = 0.3f * r.next();
        for (auto& v : x)  v = 0.4f * r.next();

        convolution::PartitionedConvolver<> conv;
        conv.prepare (P, irLen);
        conv.setIr (ir.data(), irLen);

        const int splits[] = { 1, 7, 31, 2, 64, 3, 127, 5, 16, 1, 100, 65, 8, 4, 200 };
        int off = 0, si = 0;
        while (off < n)
        {
            int blk = splits[si++ % (int) (sizeof (splits) / sizeof (splits[0]))];
            if (off + blk > n) blk = n - off;
            conv.process (x.data() + off, y.data() + off, blk);
            off += blk;
        }
        test::ok (maxAbsDiff (y, directConv (x, ir)) < 2e-3, "variable-block output == direct conv (no off-by-one)");
    }

    // --- zero latency + exact delay: IR = delta at d (incl. partition boundaries), input = delta at 0 ---
    test::group ("Conv zero-latency + exact delay at partition boundaries");
    {
        const int P = 64, n = 400;
        const int delays[] = { 0, 1, 63, 64, 65, 127, 128, 129, 200 };
        for (int d : delays)
        {
            std::vector<float> ir ((std::size_t) d + 1, 0.0f); ir[(std::size_t) d] = 1.0f;
            std::vector<float> x (n, 0.0f); x[0] = 1.0f;
            std::vector<float> y (n, 0.0f);
            convolution::PartitionedConvolver<> conv;
            conv.prepare (P, d + 1);
            conv.setIr (ir.data(), d + 1);
            conv.process (x.data(), y.data(), n);
            test::approx (y[(std::size_t) d], 1.0, 2e-4, "impulse appears at exactly delay d");
            double leak = 0.0; for (int i = 0; i < n; ++i) if (i != d) leak = std::max (leak, (double) std::fabs (y[(std::size_t) i]));
            test::ok (leak < 2e-4, "no energy off the delay");
        }
        // d == 0 specifically pins ZERO added latency: out[0] is produced in the same call.
        std::vector<float> ir { 1.0f }, x (n, 0.0f), y (n, 0.0f); x[0] = 1.0f;
        convolution::PartitionedConvolver<> conv; conv.prepare (P, 1); conv.setIr (ir.data(), 1);
        conv.process (x.data(), y.data(), 1);   // a SINGLE-sample call
        test::approx (y[0], 1.0, 1e-6, "latency 0: out[0] ready in the same 1-sample call");
        test::ok (convolution::PartitionedConvolver<>::latencySamples() == 0, "reports 0 latency");
    }

    // --- IR shorter than P (head only, no FFT tail) ---
    test::group ("Conv short IR (head only)");
    {
        const int P = 128, irLen = 40, n = 300;
        Lcg r { 7 };
        std::vector<float> ir (irLen), x (n), y (n, 0.0f);
        for (auto& v : ir) v = 0.5f * r.next();
        for (auto& v : x)  v = 0.5f * r.next();
        convolution::PartitionedConvolver<> conv; conv.prepare (P, irLen); conv.setIr (ir.data(), irLen);
        conv.process (x.data(), y.data(), n);
        test::ok (maxAbsDiff (y, directConv (x, ir)) < 2e-3, "head-only IR == direct conv");
    }

    // --- no allocation during process() ---
    test::group ("Conv no-alloc in process()");
    {
        const int P = 128, irLen = 600, n = 512;
        std::vector<float> ir (irLen, 0.01f), x (n, 0.1f), y (n, 0.0f);
        convolution::PartitionedConvolver<> conv;
        conv.prepare (P, irLen);
        conv.setIr (ir.data(), irLen);          // allocations allowed up to here
        const long before = g_allocs.load();
        conv.process (x.data(), y.data(), n);
        conv.process (x.data(), y.data(), n);   // cross a chunk boundary a few times
        const long after = g_allocs.load();
        test::okNoAlloc (after == before, "process() performed zero heap allocations");
    }

    // --- lifecycle/misuse: setIr()/process() before / after a failed prepare must not /0 or touch empty buffers ---
    // Unprepared P_==0 → setIr divided by P_ (the same fault as the ConvolutionEngine ctor bug); process wrote
    // frame_[0] into an empty buffer. A failed prepare (non-pow2 P) must stay unprepared. Guarded now.
    test::group ("PartitionedConvolver: reject setIr/process before / after failed prepare");
    {
        convolution::PartitionedConvolver<> conv;                    // NOT prepared (P_ == 0)
        std::vector<float> ir { 1.0f, 0.5f, 0.25f }, x (8, 0.0f), y (8, 0.0f); x[0] = 1.0f;
        conv.setIr (ir.data(), (int) ir.size());                     // must NOT divide by P_==0
        conv.process (x.data(), y.data(), 8);                        // must NOT write into empty frame_
        test::ok (! conv.prepare (48, 256), "prepare(non-pow2 P=48) fails");
        conv.setIr (ir.data(), (int) ir.size());                     // still unprepared → no-op
        conv.process (x.data(), y.data(), 8);
        test::ok (conv.prepare (64, 256), "prepare(P=64) works");
        conv.setIr (ir.data(), (int) ir.size());
        conv.process (x.data(), y.data(), 8);
        test::approx (y[0], 1.0, 1e-6, "convolves correctly once prepared");
    }

    return test::report();
}

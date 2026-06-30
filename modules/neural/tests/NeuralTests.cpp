// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// JUCE-free self-tests for the neural seam: a live model swap changes the output; the epoch GC frees a
// retired model only AFTER the audio thread stepped past its block (the `>` boundary — never on the
// audio thread); process() does no alloc/delete; latency follows the model; the retire queue is bounded.

#include <felitronics_test.h>
#include <felitronics/neural/NeuralStage.h>

#include <atomic>
#include <cstdio>
#include <cstdlib>
#include <memory>

static std::atomic<long> g_allocs { 0 };
void* operator new      (std::size_t s) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s ? s : 1); }
void* operator new[]    (std::size_t s) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s ? s : 1); }
void  operator delete   (void* p) noexcept { std::free (p); }
void  operator delete[] (void* p) noexcept { std::free (p); }
void  operator delete   (void* p, std::size_t) noexcept { std::free (p); }
void  operator delete[] (void* p, std::size_t) noexcept { std::free (p); }

using namespace felitronics;

// Reference backend: a gain "inference" with a static destructor counter (so GC timing is observable).
struct GainInference
{
    static inline std::atomic<int> dtors { 0 };
    float gain = 1.0f; int latency = 0;
    explicit GainInference (float g = 1.0f, int lat = 0) noexcept : gain (g), latency (lat) {}
    ~GainInference() { dtors.fetch_add (1, std::memory_order_relaxed); }
    void prepare (double, int, int) noexcept {}
    void reset() noexcept {}
    void process (float* const* io, int nc, int n) noexcept { for (int c = 0; c < nc; ++c) for (int i = 0; i < n; ++i) io[c][i] *= gain; }
    int  latencySamples() const noexcept { return latency; }
};
static_assert (neural::Inference<GainInference>, "GainInference must satisfy the seam");

int main()
{
    std::printf ("felitronics::neural tests\n");

    // --- a live swap changes the model the audio runs ---
    test::group ("NeuralStage swap changes the live model");
    {
        neural::NeuralStage<GainInference> stage; stage.prepare ({ 48000.0, 512, 2 });
        test::ok (! stage.hasModel(), "starts with no model");
        float a[4] { 1, 1, 1, 1 }, b[4] { 1, 1, 1, 1 }; float* io[2] { a, b };
        stage.process (io, 2, 4);
        test::approx (a[0], 1.0, 1e-6, "no model → passthrough");

        stage.swapPrepared (std::make_unique<GainInference> (2.0f));
        test::ok (stage.hasModel(), "has model after swap");
        float c[4] { 1, 1, 1, 1 }, d[4] { 1, 1, 1, 1 }; float* io2[2] { c, d };
        stage.process (io2, 2, 4);
        test::approx (c[0], 2.0, 1e-6, "model applied (×2)");

        stage.swapPrepared (std::make_unique<GainInference> (0.5f));
        float e[4] { 1, 1, 1, 1 }, f[4] { 1, 1, 1, 1 }; float* io3[2] { e, f };
        stage.process (io3, 2, 4);
        test::approx (e[0], 0.5, 1e-6, "swapped model applied (×0.5)");
    }

    // --- epoch GC: a retired model is freed only AFTER the audio steps past its retire block ---
    test::group ("NeuralStage epoch GC (no audio-thread delete, > boundary)");
    {
        GainInference::dtors.store (0);
        neural::NeuralStage<GainInference> stage; stage.prepare ({ 48000.0, 64, 1 });
        stage.swapPrepared (std::make_unique<GainInference> (1.0f));        // model A live
        float x[4] { 1, 1, 1, 1 }; float* io[1] { x };
        stage.process (io, 1, 4);                                           // audioBlock → 1 (A "in use")
        stage.swapPrepared (std::make_unique<GainInference> (2.0f));        // B live; A retired at block 1
        stage.collectGarbage();                                            // now=1, 1>1 false → A kept
        test::ok (GainInference::dtors.load() == 0, "retired model NOT freed at the same block");
        stage.process (io, 1, 4);                                           // audioBlock → 2
        stage.collectGarbage();                                            // now=2 > 1 → A freed
        test::ok (GainInference::dtors.load() == 1, "retired model freed after audio stepped past");
    }

    // --- no alloc + no delete during process() ---
    test::group ("NeuralStage no-alloc + no-delete in process()");
    {
        neural::NeuralStage<GainInference> stage; stage.prepare ({ 48000.0, 512, 2 });
        stage.swapPrepared (std::make_unique<GainInference> (1.5f));
        float a[512], b[512]; for (int i = 0; i < 512; ++i) { a[i] = 0.2f; b[i] = 0.2f; }
        float* io[2] { a, b };
        GainInference::dtors.store (0);
        const long beforeNew = g_allocs.load();
        stage.process (io, 2, 512); stage.process (io, 2, 512);
        test::okNoAlloc (g_allocs.load() == beforeNew, "process() did not allocate");
        test::ok (GainInference::dtors.load() == 0, "process() did not delete");
    }

    // --- latency follows the live model ---
    test::group ("NeuralStage latency");
    {
        neural::NeuralStage<GainInference> stage; stage.prepare ({ 48000.0, 512, 2 });
        test::ok (stage.latencySamples() == 0, "no model → 0 latency");
        stage.swapPrepared (std::make_unique<GainInference> (1.0f, 37));
        test::ok (stage.latencySamples() == 37, "latency reported from the model");
    }

    // --- retire queue is bounded (RT-safe): swap returns false when full (no unbounded growth) ---
    test::group ("NeuralStage retire-queue overflow");
    {
        neural::NeuralStage<GainInference, 4> stage; stage.prepare ({ 48000.0, 64, 1 });
        stage.swapPrepared (std::make_unique<GainInference> (1.0f));        // live (retire null = no-op)
        bool last = true;
        for (int i = 0; i < 6; ++i) last = stage.swapPrepared (std::make_unique<GainInference> ((float) i));  // no GC/process → fills
        test::ok (! last, "swap returns false when the retire queue is full");
    }

    return test::report();
}

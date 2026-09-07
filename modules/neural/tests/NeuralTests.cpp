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
    static inline std::atomic<int> alive { 0 }, dtors { 0 };   // alive = net live (never reset); dtors is store(0)'d by GC-timing tests
    float gain = 1.0f; int latency = 0;
    explicit GainInference (float g = 1.0f, int lat = 0) noexcept : gain (g), latency (lat) { alive.fetch_add (1, std::memory_order_relaxed); }
    ~GainInference() { alive.fetch_sub (1, std::memory_order_relaxed); dtors.fetch_add (1, std::memory_order_relaxed); }
    void prepare (double, int, int) noexcept {}
    void reset() noexcept {}
    bool refuseEverything = false;   // a backend-SPECIFIC refusal, so the stage has a verdict to lose
    [[nodiscard]] bool process (float* const* io, int nc, int n) noexcept   // law 11: the verdict is returned
    { if (nc < 0 || n < 0 || refuseEverything) return false;
      for (int c = 0; c < nc; ++c) for (int i = 0; i < n; ++i) io[c][i] *= gain; return true; }
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
        felitronics::test::run (stage.process (io, 2, 4));
        test::approx (a[0], 1.0, 1e-6, "no model → passthrough");

        stage.swapPrepared (std::make_unique<GainInference> (2.0f));
        test::ok (stage.hasModel(), "has model after swap");
        float c[4] { 1, 1, 1, 1 }, d[4] { 1, 1, 1, 1 }; float* io2[2] { c, d };
        felitronics::test::run (stage.process (io2, 2, 4));
        test::approx (c[0], 2.0, 1e-6, "model applied (×2)");

        stage.swapPrepared (std::make_unique<GainInference> (0.5f));
        float e[4] { 1, 1, 1, 1 }, f[4] { 1, 1, 1, 1 }; float* io3[2] { e, f };
        felitronics::test::run (stage.process (io3, 2, 4));
        test::approx (e[0], 0.5, 1e-6, "swapped model applied (×0.5)");
    }

    // --- epoch GC: a retired model is freed only AFTER the audio steps past its retire block ---
    test::group ("NeuralStage epoch GC (no audio-thread delete, > boundary)");
    {
        GainInference::dtors.store (0);
        neural::NeuralStage<GainInference> stage; stage.prepare ({ 48000.0, 64, 1 });
        stage.swapPrepared (std::make_unique<GainInference> (1.0f));        // model A live
        float x[4] { 1, 1, 1, 1 }; float* io[1] { x };
        felitronics::test::run (stage.process (io, 1, 4));                                           // audioBlock → 1 (A "in use")
        stage.swapPrepared (std::make_unique<GainInference> (2.0f));        // B live; A retired at block 1
        stage.collectGarbage();                                            // now=1, 1>1 false → A kept
        test::ok (GainInference::dtors.load() == 0, "retired model NOT freed at the same block");
        felitronics::test::run (stage.process (io, 1, 4));                                           // audioBlock → 2
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
        felitronics::test::run (stage.process (io, 2, 512)); felitronics::test::run (stage.process (io, 2, 512));
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

    // Leak check (local proxy for the CI's LeakSanitizer): once every NeuralStage scope above has closed,
    // LAW 11 — the stage returns the BACKEND's verdict, it does not manufacture one. A stage with no
    // model installed is a documented clean passthrough and therefore an ACCEPTED call; a stage whose
    // backend refuses must say so, or the whole point of returning a verdict is lost one layer up.
    test::group ("law 11: NeuralStage returns the backend's verdict");
    {
        neural::NeuralStage<GainInference> stage; stage.prepare ({ 48000.0, 512, 2 });
        float a[4] { 1, 1, 1, 1 }, b[4] { 1, 1, 1, 1 }; float* io[2] { a, b };
        test::ok (stage.process (io, 2, 4), "no model installed is a passthrough, and an ACCEPTED call");
        stage.swapPrepared (std::make_unique<GainInference> (2.0f));
        test::ok (stage.process (io, 2, 4), "...and so is a normal call with a model");
        test::ok (! stage.process (io, -1, 4), "a malformed call is REFUSED");
        test::ok (! stage.process (io, 2, -4), "...on either extent");
        test::ok (! stage.process (io, 3, 4),  "...and a width past the prepared spec");
        {
            // The geometry has to be judged the SAME WAY with and without a model: it used to be
            // delegated, so an EMPTY stage answered "accepted" to a malformed call and a loaded one
            // answered "refused". And the refusal must not step the retire counter the GC reads.
            neural::NeuralStage<GainInference> empty; empty.prepare ({ 48000.0, 512, 2 });
            test::ok (! empty.process (io, -1, 4), "a stage with NO model refuses a malformed call too");
            test::ok (! empty.process (io, 3, 4),  "...and a too-wide one");
        }
        {
            // ...and a refusal that is the BACKEND's own must reach the caller, not be replaced by
            // the stage's optimism.
            neural::NeuralStage<GainInference> s2; s2.prepare ({ 48000.0, 512, 2 });
            auto b = std::make_unique<GainInference> (2.0f); b->refuseEverything = true;
            s2.swapPrepared (std::move (b));
            test::ok (! s2.process (io, 2, 4), "a BACKEND-specific refusal is returned, not swallowed");
            s2.collectGarbage();
        }
        stage.collectGarbage();
    }

    // every GainInference ever constructed must be destroyed — including each stage's LIVE model, which is a
    // raw pointer freed by ~NeuralStage. Before that destructor existed, the live models leaked (LSan caught it).
    test::ok (GainInference::alive.load() == 0, "no leak: every model destroyed (net alive == 0)");

    return test::report();
}

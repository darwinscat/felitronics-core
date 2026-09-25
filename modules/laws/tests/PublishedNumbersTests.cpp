// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

//==================================================================================================
// LAW 11d (DSP-ARCHITECTURE.md §2) — THE NUMBERS A STAGE PUBLISHES ARE THE NUMBERS ITS prepare() RUNS.
// A stage's static `storageFor()` is its budget and `latencyFor()` its latency, published so a composite can
// size itself before its stages exist; the law is that neither can drift from what `prepare()` then does.
//
// These checks are core's own properties of core's stages, and until the mastering chain moved to
// felitronics-mastering-core they had exactly one gate: the chain's suite (MasteringChainTests there, which
// keeps its copy — it sizes the chain from these numbers). A satellite suite does not run on a pull request
// to this repository, so the checks are here as well, verbatim, next to the census of law 11. Only the
// stages are touched; nothing here builds a chain.
//==================================================================================================

#include <felitronics/core/DelayLine.h>
#include <felitronics/core/DryAligner.h>
#include <felitronics/core/Math.h>
#include <felitronics/dynamics/Compressor.h>
#include <felitronics/eq/EqEngine.h>
#include <felitronics/limiter/TruePeakLimiter.h>
#include <felitronics/oversampling/PolyphaseOversampler.h>
#include <felitronics/saturation/Saturator.h>
#include <felitronics/stereo/MonoBass.h>
#include <felitronics_test.h>
#include <alloc_counter.h>   // installs the allocation counter: EVERY form of `new`, over-aligned included

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

using namespace felitronics;
using felitronics::test::ok;
using felitronics::test::approx;
using felitronics::test::group;

namespace
{
constexpr double kPi = 3.14159265358979323846;
std::uint32_t bits (float f) noexcept { return std::bit_cast<std::uint32_t> (f); }

void theStagesOwnBudgets()
{
    group ("law 11d — each stage's own budget is what preparing it directly allocates, at its own clamps");

    // THE STAGES' OWN BUDGETS, ASKED DIRECTLY, at arguments no chain can reach. Three mutations survived
    // the chain's own budget matrix (MasteringChainTests) for this reason alone: the chain caps its quantum at 8192 (so the limiter's own
    // 1 Mi-sample block cap is never exercised), refuses a width past `core::kMaxChannels` (so the
    // oversampler's channel CLAMP is never exercised), and never asks a dry aligner for a capacity under 2.
    // A budget is only pinned where its own clamps can be reached.
    {
        int off = 0;
        auto measure = [&off] (const char* /*what*/, std::uint64_t budget, auto&& build)
        {
            const long long before = alloc::bytes.load();
            const bool okPrep = build();
            const long long got = alloc::bytes.load() - before;
            if (! okPrep || got != (long long) budget) ++off;
        };
        {   // the limiter's block cap: a whole-file maxBlock is a normal thing for an offline caller to pass
            limiter::TruePeakLimiterConfig lc;
            limiter::TruePeakLimiter::Storage st;
            const bool okSt = limiter::TruePeakLimiter::storageFor (48000.0, (1 << 20) + 7, 2, lc, st);
            if (! okSt) ++off;
            auto l = std::make_unique<limiter::TruePeakLimiter>();
            measure ("limiter, block past the cap", st.bytes(),
                     [&] { return l->prepare (48000.0, (1 << 20) + 7, 2, lc); });
        }
        {   // the oversampler's channel clamp — law 11(b)'s unfinished application, a known open item. The
            // budget must model what prepare() DOES, not what it ought to do.
            oversampling::PolyphaseOversampler::Storage st;
            const bool okSt = oversampling::PolyphaseOversampler::storageFor (4, core::kMaxChannels + 1, 64, st);
            if (! okSt) ++off;
            auto o = std::make_unique<oversampling::PolyphaseOversampler>();
            measure ("oversampler, width past the maximum", st.bytes(),
                     [&] { return o->prepare (4, core::kMaxChannels + 1, 64); });
        }
        {   // the aligner's floors: a capacity under 2 and a block under 1 are RAISED, not refused, and the
            // budget has to raise them too. Two channels, so the raised counts still exceed the seed the
            // aligner's constructor took (a one-channel aligner at the floor is already the seed, and then
            // preparing it asks for nothing — which would make this check pass whatever the floors did).
            const core::DryAligner::Storage st = core::DryAligner::storageFor (2, 1, 0);
            auto a = std::make_unique<core::DryAligner>();
            const long long before = alloc::bytes.load();
            a->prepare (2, 1, 0);
            const long long got = alloc::bytes.load() - before;
            if (got != (long long) st.bytes() || got != (long long) st.freshBytes()) ++off;
        }
        {   // ...and INSIDE the seed, where `bytes()` and the request part company: one channel, a ring at its
            // 2-slot floor. `freshBytes()` is the request — the scratch alone — and `bytes()` over-states it by
            // the 8 B ring it does not ask for. This is the aligner a mono compressor with no lookahead gets.
            const core::DryAligner::Storage st = core::DryAligner::storageFor (1, 256, 2);
            auto a = std::make_unique<core::DryAligner>();
            const long long before = alloc::bytes.load();
            a->prepare (1, 256, 2);
            const long long got = alloc::bytes.load() - before;
            if (got != (long long) st.freshBytes() || st.freshBytes() != 256u * sizeof (float)
                || st.bytes() != st.freshBytes() + 2u * sizeof (float)) ++off;
            // Fully inside it: nothing at all.
            const core::DryAligner::Storage seed = core::DryAligner::storageFor (1, 1, 2);
            auto b = std::make_unique<core::DryAligner>();
            const long long b0 = alloc::bytes.load();
            b->prepare (1, 1, 2);
            if (alloc::bytes.load() - b0 != 0 || seed.freshBytes() != 0u) ++off;
        }
        ok (off == 0, "each stage's own budget is what preparing it directly allocates, at the arguments its "
                      "OWN clamps live at (" + std::to_string (off) + " off)");
    }

    // AND THE CLAMPS ARE PINNED AGAINST A PROPERTY, NOT AGAINST THE ALLOCATION — which is the one thing the
    // byte-for-byte checks above CANNOT do. `prepare()` now sizes itself THROUGH `storageFor`, so a mutation
    // inside that function moves the budget and the allocation together and the comparison stays green: the
    // stand proved it, with three clamps removed one at a time and the whole suite still passing. What a
    // clamp needs is an independent statement about itself, and idempotence past the bound is the natural
    // one: beyond the clamp the answer must stop moving.
    {
        limiter::TruePeakLimiterConfig lc;
        limiter::TruePeakLimiter::Storage cap1 {}, cap2 {};
        const bool a1 = limiter::TruePeakLimiter::storageFor (48000.0, 1 << 20, 2, lc, cap1);
        const bool a2 = limiter::TruePeakLimiter::storageFor (48000.0, 1 << 24, 2, lc, cap2);
        ok (a1 && a2 && cap1.bytes() == cap2.bytes() && cap1.osBufSamples == cap2.osBufSamples,
            "the limiter's scratch stops growing at its block cap: a 16x bigger block is the same budget ("
            + std::to_string (cap1.bytes()) + " B vs " + std::to_string (cap2.bytes()) + " B)");

        oversampling::PolyphaseOversampler::Storage w1 {}, w2 {};
        const bool b1 = oversampling::PolyphaseOversampler::storageFor (4, core::kMaxChannels, 64, w1);
        const bool b2 = oversampling::PolyphaseOversampler::storageFor (4, core::kMaxChannels + 8, 64, w2);
        ok (b1 && b2 && w1.bytes() == w2.bytes(),
            "the oversampler's width stops at kMaxChannels, which is what its prepare() CLAMPS to (that "
            "clamp should be a refusal, and until it is, the budget has to describe the clamp)");

        // The two smallest published clamps, which no chain reaches and which therefore had no gate at all
        // until the fix round asked for them: a negative delay is a ring of one slot, and a window under one
        // is a window of one. Both are what their `prepare()` does, and both are stated in `storageFor`.
        ok (core::DelayLine::storageFor (-1).samples == 1 && core::DelayLine::storageFor (-1000).samples == 1
            && core::DelayLine::storageFor (0).samples == 1,
            "a negative delay budgets the one slot prepare() gives it");
        ok (limiter::detail::SlidingMax::storageFor (0).entries == 1
            && limiter::detail::SlidingMax::storageFor (-5).entries == 1,
            "a window under one budgets the one entry prepare() gives it");

        const core::DryAligner::Storage f0 = core::DryAligner::storageFor (2, 0, 0);
        const core::DryAligner::Storage f1 = core::DryAligner::storageFor (2, 1, 2);
        ok (f0.bytes() == f1.bytes() && f0.ring == 4 && f0.scratch == 2,
            "the aligner's floors are floors: a capacity of 0 budgets the 2 slots it will be given, and a "
            "block of 0 the 1 sample (" + std::to_string (f0.bytes()) + " B)");
    }

    // THE DELAY BANK LEAVES WHAT THE `assign` IT REPLACED LEFT. The old form built every line fresh, so a
    // re-prepared bank had no tap; the helper re-uses the lines, and a tap that outlived a SHRINKING
    // preparation used to index before the start of the ring (the code-review round found it, and the
    // invariant it breaks predates this work: `prepare(8); setDelay(8); prepare(2)` on a bare DelayLine).
    {
        std::vector<core::DelayLine> bank;
        core::prepareDelayBank (bank, 2, 8);
        bank[0].setDelay (8);
        core::prepareDelayBank (bank, 1, 2);
        ok (bank.size() == 1 && bank[0].capacity() == 2 && bank[0].delay() == 0,
            "a shrunk bank has the capacity asked for and NO tap, as a freshly built one does");
        core::DelayLine d;
        d.prepare (8);
        d.setDelay (8);
        d.prepare (2);
        ok (d.delay() <= d.capacity(), "and a bare line's tap is re-clamped into the capacity it was re-prepared at");
        float last = 0.0f;
        for (int i = 0; i < 16; ++i) last = d.process ((float) (i + 1));
        ok (std::isfinite (last), "PRECONDITION: reading it afterwards stays inside the ring");
    }

    // THE STAGES' OWN STATICS, against prepared stages. This is what makes the chain's aligner sizing and
    // its budget the same arithmetic the stage runs rather than a copy of it.
    {
        int off = 0;
        for (const double fs : { 44100.0, 48000.0, 192000.0 })
            for (const double look : { 0.0, 1.0, 10.0, 100.0, 250.0 })
            {
                dynamics::Compressor c;
                const double maxLook = std::max (look, 1.0);
                if (! c.prepare (fs, 256, 2, maxLook)) { ++off; continue; }
                dynamics::CompressorParams cp; cp.lookaheadMs = look; c.setParams (cp);
                if (c.latencySamples() != dynamics::Compressor::latencyFor (fs, 256, 2, maxLook, look)) ++off;

                limiter::TruePeakLimiterConfig lc; lc.lookaheadMs = look > 20.0 ? 20.0 : look;
                limiter::TruePeakLimiter l;
                if (! l.prepare (fs, 256, 2, lc)) { ++off; continue; }
                if (l.latencySamples() != limiter::TruePeakLimiter::latencyFor (fs, 256, 2, lc)) ++off;
                if (l.oversampleFactor() != limiter::TruePeakLimiter::oversampleFactorFor (lc)) ++off;

                // EVERY FACTOR, the ones BELOW two included: the saturator's own domain has a second half
                // where the oversampler is not built at all and the latency is 0 whatever the taps say. The
                // chain never asks for it (its `admits` requires a factor of 2 or more), so nothing else in
                // this suite reaches it — measured on the stand, a `latencyFor` that ignored the factor
                // entirely and answered `tpp - 1` survived the whole suite.
                for (const int f : { -1, 0, 1, 2, 4, 16 })
                {
                    saturation::Saturator sat;
                    if (! sat.prepare (fs, 256, 2, f, 64)) { ++off; continue; }
                    if (sat.latencySamples() != saturation::Saturator::latencyFor (fs, 256, 2, f, 64)) ++off;
                }
            }
        ok (off == 0, "every stage's latencyFor() is the latency the prepared stage reports ("
                      + std::to_string (off) + " off)");
    }
}

void theEqEngineObject()
{
    group ("eq::EqEngine — objectBytes() is what building one asks the heap for");

    // AND THAT IT SEES AN OVER-ALIGNED `new` AT ALL — the form `eq::EqEngine` is built through. A counter
    // blind to it reads the default geometry's largest single request as zero and every budget below passes,
    // which is how this tree's counter once failed, and why this check is a PRECONDITION rather than a nicety.
    {
        const long long before = alloc::bytes.load();
        {
            auto e = std::make_unique<eq::EqEngine>();
            volatile const void* sink = e.get();
            (void) sink;
        }
        const long long got = alloc::bytes.load() - before;
        ok (got == (long long) eq::EqEngine::objectBytes(),
            "the counter sees the over-aligned `new` the EQ engine is built through: "
            + std::to_string (got) + " B");
    }
}

void theMonoBassParams()
{
    group ("stereo::MonoBassParams — the added type is exactly the three setters");

    // The stage grew a parameter type in this change; the point of the type is that it CANNOT drift
    // from the setters, so that is what gets checked: same object, bit-identical audio.
    const int n = 4000;
    std::vector<float> l1 ((std::size_t) n), r1 ((std::size_t) n), l2, r2;
    for (int i = 0; i < n; ++i)
    {
        l1[(std::size_t) i] = 0.6f * (float) std::sin (2.0 * kPi * 80.0 * i / 48000.0)
                            + 0.3f * (float) std::sin (2.0 * kPi * 900.0 * i / 48000.0);
        r1[(std::size_t) i] = -0.55f * l1[(std::size_t) i];
    }
    l2 = l1; r2 = r1;

    stereo::MonoBass a, b;
    felitronics::test::run (a.prepare (48000.0, n, 2));
    felitronics::test::run (b.prepare (48000.0, n, 2));
    a.setEnabled (true); a.setFrequency (137.0f); a.setLowWidth (0.25f);
    b.setParams ({ true, 137.0f, 0.25f });
    float* pa[2] { l1.data(), r1.data() };
    float* pb[2] { l2.data(), r2.data() };
    felitronics::test::run (a.process (pa, 2, n));
    felitronics::test::run (b.process (pb, 2, n));
    long long bad = 0;
    for (int i = 0; i < n; ++i)
        if (bits (l1[(std::size_t) i]) != bits (l2[(std::size_t) i]) || bits (r1[(std::size_t) i]) != bits (r2[(std::size_t) i])) ++bad;
    ok (bad == 0, "setParams() is bit-identical to the three setters");

    const auto rp = b.params();
    ok (rp.enabled && core::exactlyEqual (rp.frequencyHz, 137.0f) && core::exactlyEqual (rp.lowWidth, 0.25f),
        "params() reads back what was applied");
    b.setFrequency (2.0f);                                 // below kMinFreq
    approx ((double) b.params().frequencyHz, 20.0, 1e-4, "params() reports the CLAMPED frequency, not the request");
    b.setLowWidth (std::nan (""));
    approx ((double) b.params().lowWidth, 0.25, 1e-9, "a non-finite width is rejected, the last good value stands");
}
} // namespace

int main()
{
    std::printf ("felitronics law 11d — the numbers a stage publishes are the numbers its prepare() runs\n");

    theStagesOwnBudgets();
    theEqEngineObject();
    theMonoBassParams();

    return felitronics::test::report();
}

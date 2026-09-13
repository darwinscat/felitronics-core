// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/core/Config.h>
#include <felitronics/mastering/DeliveryConverter.h>
#include <felitronics/mastering/LoudnessSolver.h>
#include <felitronics/mastering/MasteringChain.h>
#include <felitronics/mastering/OfflineRenderer.h>

#include <climits>
#include <cmath>
#include <cstdint>
#include <vector>

namespace felitronics::mastering
{

//==============================================================================
// felitronics::mastering::DeliveredMastering — the three whole-programme operations of a delivery, SRC FIRST:
// a render at fixed parameters, the loudness search, and the input's loudness range. Each one converts the
// programme to the delivery rate with `DeliveryConverter` and then hands it, unchanged, to the class that
// already does the job at one rate — `OfflineRenderer`, `TargetLoudnessSolver` — on a chain prepared at the
// DELIVERY rate.
//
// WHY THIS IS A CLASS AND NOT THREE LINES IN EACH CALLER. Two paths deliver a programme: the C ABI facade and
// the direct C++ call `fcore_master selftest` compares it with, bit for bit. What sits between the converter
// and the solver — which buffer the converted programme lives in, when it is allocated, that a render
// converts straight into its output and renders in place, that a range too short to measure is refused
// before a sample is converted — is exactly the kind of composition two hand-written copies get wrong in the
// same way. Written here once, both paths call it, and the bit-compare then compares an ABI against a core
// rather than a copy against a copy.
//
// WHAT IT OWNS: the converter, prepared for one (source rate, delivery rate, width, block). WHAT IT DOES NOT:
// the chain, the renderer and the solver stay the caller's, because they already have a lifecycle of their
// own (configure, channel weights, lazy preparation) and a second owner would be a second place to get it
// wrong. Every operation checks that the chain AND the solver it is handed run at the delivery rate: either one
// prepared at the source rate would render or meter a converted programme at the wrong speed, every number
// plausible.
//
// A REFUSAL WRITES NOTHING AND ALLOCATES NOTHING (law 11b). Every verdict an operation can reach is reached before
// the converter writes a sample and before a programme buffer is asked for: the renderer's width and preparation,
// overlapping planes, and — for the search — everything `TargetLoudnessSolver::admits` would refuse before its first
// pass. Converting first and letting the renderer or the solver say no afterwards left the caller's output full of
// unmastered, resampled audio, or spent a whole programme's memory on a request the solver was always going to
// refuse.
//
// MEMORY (law 11d). `render` asks the heap for NOTHING: it converts into the caller's output and renders there
// in place, which `OfflineRenderer` supports. `solve` and `measureInputLoudnessRange` cannot do that — a search
// reads its input again on every pass and may not render over it — so each holds the converted programme for
// the length of the call, and `solveBytes` / `measureRangeBytes` say how much, through the very functions the
// solver's own budgets are made of. At EQUAL RATES there is nothing to convert and the caller's input is read
// directly, so no programme buffer exists and the budget says so.
//==============================================================================
class DeliveredMastering
{
public:
    [[nodiscard]] static long long deliveredFrames (double sourceRate, double deliveryRate, long long inFrames) noexcept
    {
        return DeliveryConverter::deliveredFrames (sourceRate, deliveryRate, inFrames);
    }

    // WHAT A DELIVERING `create` ASKS THE HEAP FOR, AND WHETHER IT CAN BE BUILT AT ALL — one expression for both,
    // as `createBytes` is for a chain at one rate: the chain and the renderer at the DELIVERY rate plus the
    // converter. 0 for exactly the geometries either half refuses, so a facade that decides by this number and
    // publishes it cannot disagree with itself.
    [[nodiscard]] static std::uint64_t createBytes (double sourceRate, double deliveryRate, int numChannels,
                                                    const MasteringChainConfig& config, int rendererBlock) noexcept
    {
        const std::uint64_t chain = mastering::createBytes (deliveryRate, numChannels, config, rendererBlock);
        const std::uint64_t conv  = DeliveryConverter::prepareBytes (sourceRate, deliveryRate, numChannels, rendererBlock);
        return (chain == 0u || conv == 0u) ? 0u : chain + conv;
    }

    // solve(): the converted programme, held for the whole call, plus the search's own PEAK at the delivered
    // length (one pass — see `TargetLoudnessSolver::solveBytes`). 0 where the call converts nothing: an empty or
    // unrepresentable programme, which the solver refuses before any pass.
    [[nodiscard]] static std::uint64_t solveBytes (double sourceRate, double deliveryRate, int numChannels,
                                                   long long inFrames) noexcept
    {
        const long long d = deliveredFrames (sourceRate, deliveryRate, inFrames);
        if (d <= 0 || d > INT_MAX || numChannels < 1 || numChannels > core::kMaxChannels) return 0u;
        return programmeBytes (sourceRate, deliveryRate, numChannels, d)
             + TargetLoudnessSolver::solveBytes (deliveryRate, numChannels, (int) d);
    }

    // measureInputLoudnessRange(): the converted programme and one meter — and NOTHING for a delivered length too
    // short to have a range, which the call refuses before converting a sample. THE SHORTNESS IS JUDGED ON THE
    // DELIVERED LENGTH, which is what the meter sees: 132,300 frames at 44.1 kHz are 1.38 s of input frames read
    // at 96 kHz but 3.0 s of programme, and a budget computed from the input count said 0 for a call that builds
    // a meter.
    [[nodiscard]] static std::uint64_t measureRangeBytes (double sourceRate, double deliveryRate, int numChannels,
                                                          long long inFrames) noexcept
    {
        const long long d = deliveredFrames (sourceRate, deliveryRate, inFrames);
        if (d <= 0 || d > INT_MAX || numChannels < 1 || numChannels > core::kMaxChannels) return 0u;
        const std::uint64_t meter = TargetLoudnessSolver::measureRangeBytes (deliveryRate, (int) d);
        return meter == 0u ? 0u : programmeBytes (sourceRate, deliveryRate, numChannels, d) + meter;
    }

    [[nodiscard]] static std::uint64_t prepareBytes (double sourceRate, double deliveryRate, int numChannels, int block) noexcept
    {
        return DeliveryConverter::prepareBytes (sourceRate, deliveryRate, numChannels, block);
    }

    [[nodiscard]] bool prepare (double sourceRate, double deliveryRate, int numChannels, int block)
    {
        prepared_ = false;
        if (! conv_.prepare (sourceRate, deliveryRate, numChannels, block)) return false;
        sourceRate_ = sourceRate; deliveryRate_ = deliveryRate; nch_ = numChannels;
        identity_ = conv_.plan().identity;
        prepared_ = true;
        return true;
    }

    bool   isPrepared()   const noexcept { return prepared_; }
    double sourceRate()   const noexcept { return sourceRate_; }
    double deliveryRate() const noexcept { return deliveryRate_; }
    const DeliveryConverter& converter() const noexcept { return conv_; }

    // Non-finite input samples in the programme the last render, solve or range measurement was handed — counted by
    // the converter's gate, or at equal rates (where the input is read in place and the chain's gate replaces them)
    // by the same test over the input. A number about the CALLER's programme, at the source rate.
    std::uint64_t nonFiniteInputSamples() const noexcept { return nonFinite_; }

    // A render at the parameters the chain already holds. `out` must be exactly `deliveredFrames(inFrames)`
    // frames per channel and `outFrames` that number; the planes must be `planesUsable` (Planes.h) — the conversion
    // writes `out` while it still reads `in`, at another stride, and the render then runs in place over `out`. No
    // allocation.
    [[nodiscard]] bool render (MasteringChain& chain, OfflineRenderer& renderer,
                               const float* const* in, int numChannels, long long inFrames,
                               float* const* out, long long outFrames) noexcept
    {
        if (! admits (chain, numChannels, inFrames, outFrames)) return false;
        if (renderer.blockSize() < 1 || numChannels > renderer.maxChannels()) return false;
        // THE RULE AT EVERY LENGTH, AN EMPTY PROGRAMME INCLUDED. It used to be skipped at `outFrames == 0`, and the
        // path then ran past it into `ro[c] = out[c]` below: a null table crashed a call whose programme is legal.
        // At two zero lengths the rule judges nothing but null — no span has a byte to overlap with — so an empty
        // programme in real tables is rendered as before, and one without tables is refused.
        if (! planesUsable (in, out, numChannels, inFrames, outFrames)) return false;
        if (! conv_.convert (in, numChannels, inFrames, out, outFrames)) return false;
        nonFinite_ = conv_.nonFiniteInputSamples();
        const float* ro[core::kMaxChannels] {};
        for (int c = 0; c < numChannels; ++c) ro[c] = out[c];
        return renderer.render (chain, ro, out, numChannels, (int) outFrames);
    }

    // The loudness search over the delivered programme. Refusals of this class's own are the solver's verdict
    // type and the solver's words for them: a converter that is not ready is `NotPrepared`, a length, width,
    // rate or buffer that does not fit is `InvalidRequest` — the same two a caller already has to read.
    LoudnessSolution solve (TargetLoudnessSolver& solver, MasteringChain& chain, OfflineRenderer& renderer,
                            const MasteringChainParams& params,
                            const float* const* in, int numChannels, long long inFrames,
                            float* const* out, long long outFrames, const LoudnessRequest& req)
    {
        LoudnessSolution refused;
        refused.activityThresholdDb = req.activityThresholdDb;
        if (! prepared_) { refused.status = MasteringSolveStatus::NotPrepared; return refused; }
        if (! admits (chain, numChannels, inFrames, outFrames))
            { refused.status = MasteringSolveStatus::InvalidRequest; return refused; }
        // AN EMPTY PROGRAMME IS THE SOLVER'S TO ANSWER, not this class's: it answers `InvalidRequest` before any
        // pass, and answering for it here would be a second policy for one question.
        if (outFrames == 0)
            return solver.solve (chain, renderer, params, in, out, numChannels, 0, req);
        // THE SOLVER'S OWN VERDICT, ASKED BEFORE A BYTE IS SPENT. Its words, its order, one definition.
        if (! solver.admits (chain, renderer, numChannels, (int) outFrames, req, refused.status)) return refused;
        // The planes, at the CALLER's two lengths — not left to the solver, which sees the converted programme and
        // never the caller's input. (At equal rates the search reads `in` directly on every pass, so an overlap there
        // would read its own master.)
        if (! planesUsable (in, out, numChannels, inFrames, outFrames))
            { refused.status = MasteringSolveStatus::InvalidRequest; return refused; }

        const float* src[core::kMaxChannels] {};
        std::vector<float> programme;
        if (! converted (in, numChannels, inFrames, outFrames, programme, src))
            { refused.status = MasteringSolveStatus::InvalidRequest; return refused; }
        return solver.solve (chain, renderer, params, src, out, numChannels, (int) outFrames, req);
    }

    // The input's loudness range, measured on the DELIVERED programme — the one the search will meter, so the
    // range constraint compares a programme with itself. False, having converted nothing, where the solver
    // would refuse a range for that length.
    [[nodiscard]] bool measureInputLoudnessRange (const TargetLoudnessSolver& solver, const float* const* in,
                                                  int numChannels, long long inFrames, double& out)
    {
        if (! prepared_ || ! solver.isPrepared()) return false;
        if (! (std::fabs (solver.sampleRate() - deliveryRate_) < 1.0e-9)) return false;   // meters at the delivery rate
        if (in == nullptr || numChannels != nch_ || inFrames <= 0) return false;
        for (int c = 0; c < numChannels; ++c) if (in[c] == nullptr) return false;
        const long long d = deliveredFrames (sourceRate_, deliveryRate_, inFrames);
        if (d <= 0 || d > INT_MAX) return false;
        if (TargetLoudnessSolver::measureRangeBytes (deliveryRate_, (int) d) == 0u) return false;
        const float* src[core::kMaxChannels] {};
        std::vector<float> programme;
        if (! converted (in, numChannels, inFrames, d, programme, src)) return false;
        return solver.measureInputLoudnessRange (src, numChannels, (int) d, out);
    }

private:
    static std::uint64_t programmeBytes (double sourceRate, double deliveryRate, int numChannels, long long d) noexcept
    {
        core::DeliveryResampler::Params p;
        p.inRate = sourceRate; p.outRate = deliveryRate;
        if (core::DeliveryResampler::plan (p).identity) return 0u;       // the caller's input is read in place
        return (std::uint64_t) sizeof (float) * (std::uint64_t) numChannels * (std::uint64_t) d;
    }

    bool admits (const MasteringChain& chain, int numChannels, long long inFrames, long long outFrames) const noexcept
    {
        if (! prepared_ || numChannels != nch_ || inFrames < 0 || outFrames < 0 || outFrames > INT_MAX) return false;
        if (outFrames != deliveredFrames (sourceRate_, deliveryRate_, inFrames)) return false;
        // THE RATE IS CHECKED, as the solver checks its own: a chain at the source rate would play a converted
        // programme at the wrong speed and every number it reported would still look like a measurement.
        return chain.isPrepared() && chain.numChannels() == numChannels
            && std::fabs (chain.sampleRate() - deliveryRate_) < 1.0e-9;
    }

    // The programme at the delivery rate, in `src`. At equal rates that is the caller's own input, read in place.
    bool converted (const float* const* in, int numChannels, long long inFrames, long long outFrames,
                    std::vector<float>& programme, const float** src)
    {
        if (identity_)
        {
            // Read in place; the chain's own gate replaces a bad sample one for one, so only the count is taken here.
            std::uint64_t bad = 0;
            for (int c = 0; c < numChannels; ++c)
            {
                src[c] = in[c];
                for (long long i = 0; i < inFrames; ++i) bad += std::isfinite (in[c][i]) ? 0u : 1u;
            }
            nonFinite_ = bad;
            return true;
        }
        programme.assign ((std::size_t) numChannels * (std::size_t) outFrames, 0.0f);
        float* dst[core::kMaxChannels] {};
        for (int c = 0; c < numChannels; ++c)
        {
            dst[c] = programme.data() + (std::size_t) c * (std::size_t) outFrames;
            src[c] = dst[c];
        }
        const bool ok = conv_.convert (in, numChannels, inFrames, dst, outFrames);
        nonFinite_ = conv_.nonFiniteInputSamples();
        return ok;
    }

    DeliveryConverter conv_;
    double sourceRate_ = 0.0, deliveryRate_ = 0.0;
    int nch_ = 0;
    std::uint64_t nonFinite_ = 0;
    bool identity_ = false, prepared_ = false;
};

} // namespace felitronics::mastering

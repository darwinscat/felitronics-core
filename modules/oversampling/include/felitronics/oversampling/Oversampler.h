// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/oversampling/CascadeOversampler.h>
#include <felitronics/oversampling/PolyphaseOversampler.h>

#include <cstddef>
#include <cstdint>
#include <vector>

namespace felitronics::oversampling
{

//==============================================================================
// Which oversampler a nonlinear stage wraps (P31). A topology choice, fixed for a prepared stream.
//
//   Kaiser   — PolyphaseOversampler: one stage, cutoff fixed at 0.45 fs, `tapsPerPhase` taps per phase.
//              Round trip tapsPerPhase - 1. At 44.1 kHz -1.80 dB at 19 kHz and -15.55 at 20 kHz.
//              THE DEFAULT everywhere, and the only one for a factor that is not a power of two.
//   Cascade  — CascadeOversampler: strict, flat to 20 kHz at every rate, lengths derived from the sample
//              rate. Round trip 131 base samples at 44.1 kHz, 76 at 48 kHz, 28 at 88.2 kHz (4x).
//              Powers of two only.
enum class Topology { Kaiser, Cascade };

// One member that is either, for a stage that offers both. Under Kaiser it IS PolyphaseOversampler — the
// same calls with the same arguments, the same refusals, the same bits — so a stage that adopts this type
// changes nothing for a caller that does not ask for the cascade.
//
// `tapsPerPhase` configures Kaiser only. Under Cascade it is still RANGE-CHECKED ([4, kMaxTapsPerPhase]),
// so the set of arguments a stage refuses does not shrink when the topology changes; it is otherwise not
// used — the cascade takes its lengths from the rate.
//
// THE CASCADE LIVES ON THE HEAP, and only when it is chosen. Held by value it added about half a kilobyte
// (sizeof (CascadeOversampler)) to every Saturator and TruePeakLimiter whether or not anyone asked for it,
// and twice that to a MasteringChain — a change to the default. A vector of zero or one keeps the switch
// copyable and costs the Kaiser path 32 bytes (sizeof (Oversampler) 200 against 168; Saturator 432 -> 464,
// TruePeakLimiter 408 -> 440, MasteringChain 18016 -> 18080 on arm64); the object it allocates is counted
// in Storage::heapObjects, so a stage's published budget still equals what its prepare() asks for.
class Oversampler
{
public:
    struct Storage
    {
        Topology topology = Topology::Kaiser;
        PolyphaseOversampler::Storage kaiser {};     // empty under Cascade
        CascadeOversampler::Storage   cascade {};    // empty under Kaiser
        std::size_t heapObjects = 0;                 // CascadeOversampler instances the switch allocates (0 or 1)
        std::uint64_t bytes() const noexcept
        {
            return kaiser.bytes() + cascade.bytes() + (std::uint64_t) heapObjects * sizeof (CascadeOversampler);
        }
        bool fitsWithin (const Storage& other) const noexcept
        {
            return kaiser.fitsWithin (other.kaiser) && cascade.fitsWithin (other.cascade) && heapObjects <= other.heapObjects;
        }
    };

    [[nodiscard]] static bool storageFor (Topology topology, double sampleRate, int factor, int maxChannels,
                                          int tapsPerPhase, Storage& out) noexcept
    {
        Storage st;
        st.topology = topology;
        if (topology == Topology::Kaiser)
        {
            if (! PolyphaseOversampler::storageFor (factor, maxChannels, tapsPerPhase, st.kaiser)) return false;
        }
        else
        {
            if (tapsPerPhase < 4 || tapsPerPhase > PolyphaseOversampler::kMaxTapsPerPhase) return false;
            if (! CascadeOversampler::storageFor (sampleRate, factor, maxChannels, st.cascade)) return false;
            st.heapObjects = 1;
        }
        out = st;
        return true;
    }

    // The round trip a preparation with these arguments will report, without preparing one. Under Kaiser
    // this is the formula every stage already used (tapsPerPhase - 1, independent of the factor and the
    // rate); 0 where the cascade refuses.
    static int latencyFor (Topology topology, double sampleRate, int factor, int tapsPerPhase) noexcept
    {
        if (topology == Topology::Kaiser) return tapsPerPhase > 0 ? tapsPerPhase - 1 : 0;
        return CascadeOversampler::latencyFor (sampleRate, factor);
    }

    // Refused arguments touch nothing. A successful preparation of one topology releases the other's memory.
    bool prepare (Topology topology, double sampleRate, int factor, int maxChannels, int tapsPerPhase)
    {
        if (topology == Topology::Kaiser)
        {
            if (! kaiser_.prepare (factor, maxChannels, tapsPerPhase)) return false;
            std::vector<CascadeOversampler>().swap (cascade_);        // frees; allocates nothing
        }
        else
        {
            Storage st;
            if (! storageFor (topology, sampleRate, factor, maxChannels, tapsPerPhase, st)) return false;
            if (cascade_.empty()) cascade_.emplace_back();
            if (! cascade_.front().prepare (sampleRate, factor, maxChannels)) return false;
            kaiser_ = PolyphaseOversampler {};
        }
        topology_ = topology;
        return true;
    }

    // Under Cascade every call goes through `cascade()`, which is null when the vector is empty — and it IS
    // empty in a MOVED-FROM switch, whose topology_ was copied while the vector was stolen. A moved-from
    // switch therefore reads as unprepared (latency 0, factor 0, calls are no-ops) instead of dereferencing
    // an empty vector: measured as a SIGSEGV on `latencySamples()` of a moved-from Saturator before this.
    void reset() noexcept
    {
        if (topology_ == Topology::Kaiser) kaiser_.reset();
        else if (auto* c = cascade()) c->reset();
    }
    void resetChannel (int ch) noexcept
    {
        if (topology_ == Topology::Kaiser) kaiser_.resetChannel (ch);
        else if (auto* c = cascade()) c->resetChannel (ch);
    }
    int factor() const noexcept
    {
        if (topology_ == Topology::Kaiser) return kaiser_.factor();
        const auto* c = cascade();
        return c != nullptr ? c->factor() : 0;
    }
    int latencySamples() const noexcept
    {
        if (topology_ == Topology::Kaiser) return kaiser_.latencySamples();
        const auto* c = cascade();
        return c != nullptr ? c->latencySamples() : 0;
    }
    Topology topology() const noexcept { return topology_; }

    void upsample (const float* const* in, int channels, int n, float* const* out) noexcept
    {
        if (topology_ == Topology::Kaiser) kaiser_.upsample (in, channels, n, out);
        else if (auto* c = cascade()) c->upsample (in, channels, n, out);
    }
    void downsample (const float* const* in, int channels, int n, float* const* out) noexcept
    {
        if (topology_ == Topology::Kaiser) kaiser_.downsample (in, channels, n, out);
        else if (auto* c = cascade()) c->downsample (in, channels, n, out);
    }

private:
    CascadeOversampler*       cascade() noexcept       { return cascade_.empty() ? nullptr : &cascade_.front(); }
    const CascadeOversampler* cascade() const noexcept { return cascade_.empty() ? nullptr : &cascade_.front(); }

    Topology             topology_ = Topology::Kaiser;
    PolyphaseOversampler kaiser_;
    std::vector<CascadeOversampler> cascade_;    // empty, or the one prepared cascade
};

} // namespace felitronics::oversampling

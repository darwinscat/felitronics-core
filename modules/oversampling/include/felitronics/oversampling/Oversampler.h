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
// THE CASCADE LIVES ON THE HEAP, and only when it is chosen. Held by value it would have added ~500 bytes
// to every Saturator and TruePeakLimiter whether or not anyone asked for it (measured: 432 -> 936 and
// 408 -> 912, MasteringChain 18016 -> 19024), which is a change to the default. A vector of zero or one
// keeps the switch copyable and costs the Kaiser path its 24 bytes; the object it allocates is counted
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

    // Invariant: topology_ == Cascade implies cascade_ holds exactly one prepared object.
    void reset() noexcept                { if (topology_ == Topology::Kaiser) kaiser_.reset(); else cascade_.front().reset(); }
    void resetChannel (int c) noexcept   { if (topology_ == Topology::Kaiser) kaiser_.resetChannel (c); else cascade_.front().resetChannel (c); }
    int  factor() const noexcept         { return topology_ == Topology::Kaiser ? kaiser_.factor() : cascade_.front().factor(); }
    int  latencySamples() const noexcept { return topology_ == Topology::Kaiser ? kaiser_.latencySamples() : cascade_.front().latencySamples(); }
    Topology topology() const noexcept   { return topology_; }

    void upsample (const float* const* in, int channels, int n, float* const* out) noexcept
    {
        if (topology_ == Topology::Kaiser) kaiser_.upsample (in, channels, n, out);
        else                               cascade_.front().upsample (in, channels, n, out);
    }
    void downsample (const float* const* in, int channels, int n, float* const* out) noexcept
    {
        if (topology_ == Topology::Kaiser) kaiser_.downsample (in, channels, n, out);
        else                               cascade_.front().downsample (in, channels, n, out);
    }

private:
    Topology             topology_ = Topology::Kaiser;
    PolyphaseOversampler kaiser_;
    std::vector<CascadeOversampler> cascade_;    // empty, or the one prepared cascade
};

} // namespace felitronics::oversampling

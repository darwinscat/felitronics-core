// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <cstdint>

namespace felitronics::mastering
{

//==============================================================================
// felitronics::mastering::planesUsable — the planes a whole-programme operation may read and write, judged once for
// the module. Every operation that takes `in` and `out` tables here reads its input while it writes its output, and
// none of them may do that over memory it shares: a search reads its input again on every pass, a conversion writes
// `out` at one stride while still reading `in` at another — so whether an overlap survives depends on the ratio and
// the offset, which is not a contract — and a render of two channels into one buffer keeps only the second channel's
// samples. So a call is usable when
//
//   * both tables, and every plane in them, are non-null;
//   * no INPUT plane's bytes touch any OUTPUT plane's — every pair, half-open, EACH SIDE AT ITS OWN LENGTH:
//     `[in[c], in[c] + inFrames)` against `[out[k], out[k] + outFrames)`. A conversion's two lengths differ, and
//     judging both at one of them either misses an overlap that lies only in the longer span or refuses planes that
//     never meet;
//   * no two OUTPUT planes' bytes touch — `[out[k], out[k] + outFrames)` against `[out[j], ...)`, k != j.
//
// What stays legal is what never overlaps: one buffer feeding two INPUT channels (nothing writes it), planes edge to
// edge in one allocation, a call shorter than its buffers (judged on the frames it uses). THE LENGTHS ARE NON-NEGATIVE,
// AND ZERO ONLY TOGETHER: at two zero lengths no span has a byte, so the rule judges null tables and planes and nothing
// else — `DeliveredMastering::render` asks it so for an empty programme. One zero beside a non-zero length is not
// "touches nothing" (an empty span whose address lies strictly inside the other still counts as inside it), and a
// negative length wraps into a meaningless span; no caller asks either.
//
// ONE DEFINITION, because two hand-written copies had already drifted: the loudness search tested only
// `in[c] == out[c]` and no null plane, and accepted `out[0] = in[1]` — a master read back as the next pass's input —
// while `DeliveredMastering` refused it. `TargetLoudnessSolver::solve` asks with one length for both sides,
// `DeliveryConverter::convert` and `DeliveredMastering` with a conversion's two. It trusts `numChannels`: every caller
// has checked the width, and the tables are that long.
[[nodiscard]] inline bool planesUsable (const float* const* in, float* const* out, int numChannels,
                                        long long inFrames, long long outFrames) noexcept
{
    if (in == nullptr || out == nullptr) return false;
    const auto bytesIn  = (std::uint64_t) inFrames  * sizeof (float);
    const auto bytesOut = (std::uint64_t) outFrames * sizeof (float);
    const auto at = [] (const float* p) noexcept { return (std::uint64_t) reinterpret_cast<std::uintptr_t> (p); };
    for (int c = 0; c < numChannels; ++c)
        if (in[c] == nullptr || out[c] == nullptr) return false;
    for (int k = 0; k < numChannels; ++k)
    {
        const std::uint64_t o = at (out[k]);
        for (int c = 0; c < numChannels; ++c)
        {
            const std::uint64_t i = at (in[c]);
            if (i < o + bytesOut && o < i + bytesIn) return false;
        }
        for (int j = k + 1; j < numChannels; ++j)
        {
            const std::uint64_t p = at (out[j]);
            if (p < o + bytesOut && o < p + bytesOut) return false;
        }
    }
    return true;
}

} // namespace felitronics::mastering

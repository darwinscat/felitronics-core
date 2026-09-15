// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

// The `clips` stream, formatted — the C++ half of a format whose JavaScript half is tools/wasm/clips-format.mjs.
//
// SEPARATE FROM fcore_clips.h ON PURPOSE. That header is the measurement wiring and is compiled INTO THE WASM
// MODULE; this one is the CLI's text and is not. The wasm side formats in JavaScript, as it does for `blocks`
// and for the shapes, so linking snprintf and a std::string builder into the module would be payload for a
// function it never calls — dead-stripped at -O3, but the checked artifact is built -O1 -g and the module's
// size is an acceptance criterion of its own (tools/wasm/build.sh).

#include "fcore_clips.h"

#include <bit>
#include <cstdint>
#include <cstdio>
#include <string>

namespace fcore
{

// THE FORMAT. Every floating-point number goes out as a raw IEEE-754 bit pattern — 16 lowercase hex digits —
// and never as decimal or %a: the other side of this comparison is JavaScript, which has no hex-float
// printing and whose decimal formatting is not C's, so a 16-digit pattern is the one representation both
// sides produce identically and `diff` IS the parity test. The same rule and the same reason as
// tools/wasm/blocks-format.mjs. It also makes the comparison sharp: a single flipped bit in a level or a
// peak changes a character here, where "%.17g" would still print it and "%.6f" would hide it entirely.
//
// The sample rate is a bit pattern for the same reason: C's %g and JavaScript's Number-to-string do not
// agree on a fractional rate. The integers are plain decimal — every one of them is bounded by the frame
// count, so it is below 2^32 and prints identically through printf("%lld") and JavaScript's String(Number).
//
// NOTHING IS IN dB. samplePeakDb() routes through core::gainToDb → log10, and libm is not bit-identical
// across toolchains; a dB column would break the diff while every measured bit matched.
// A report that is not `ok` formats as NOTHING, not as a plausible empty report. readClips() already returns
// false and clears its output, but a caller can forget a bool; it cannot forget an empty stdout, and the
// alternative — printing `runs 0 / complete 1` for a measurement that never happened — is the certificate of
// cleanliness this whole path exists to make impossible.
inline std::string formatClips (const ClipsReport& r)
{
    if (! r.ok) return {};
    auto bits = [] (double d) noexcept { return std::bit_cast<std::uint64_t> (d); };
    // 256, against a worst case of 137 for the header (five int64 at 20 digits plus a 16-digit pattern) and
    // 100 for a run line. snprintf would truncate rather than overflow, but a truncated line is a format
    // that silently stopped matching the JavaScript half — the diff would then report a measurement failure
    // for what is really a buffer.
    char line[256];
    std::string out;
    std::snprintf (line, sizeof line, "# fcore clips v1 sr=%016llx ch=%d frames=%lld maxruns=%lld delay=%lld\n",
                   (unsigned long long) bits (r.sampleRate), r.channels, (long long) r.frames,
                   (long long) r.maxRuns, (long long) r.decisionDelay);
    out += line;
    std::snprintf (line, sizeof line, "runs %lld stored %lld complete %d\n",
                   (long long) r.runCount, (long long) r.storedRunCount, r.complete ? 1 : 0);
    out += line;
    for (std::size_t c = 0; c < r.peak.size(); ++c)          // indexed by the vector, never by r.channels:
    {                                                        // a hand-built report cannot walk off the end
        std::snprintf (line, sizeof line, "peak %d %016llx\n", (int) c, (unsigned long long) bits (r.peak[c]));
        out += line;
    }
    for (const auto& u : r.runs)
    {
        std::snprintf (line, sizeof line, "run %lld %lld %016llx %d %d %d\n", (long long) u.start,
                       (long long) u.length, (unsigned long long) bits (u.level), u.channel, u.sign,
                       (int) u.evidence);
        out += line;
    }
    return out;
}

} // namespace fcore

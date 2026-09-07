// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <bit>
#include <cmath>
#include <cstdint>

namespace felitronics::core
{

constexpr double kPi = 3.14159265358979323846;

// dB <-> linear amplitude (20·log10). `double` on purpose: these are offline / coefficient-design /
// GUI helpers, never the per-sample loop (Law 3 carve-out). The floor keeps log10 finite.
inline double dbToGain (double dB)   noexcept { return std::pow (10.0, dB / 20.0); }

// THE FLOOR `gainToDb` PUTS UNDER ITS ARGUMENT IS LOAD-BEARING, so it is named and then USED rather than
// repeated. Law 11c ("a pause is silence") collapses a stretch of digital silence into two cheap
// recurrences the moment the detector level reaches this floor, because from there the dB conversion —
// and the whole static curve behind it — returns the SAME BITS for every smaller level. A constant
// sitting BESIDE the function instead of inside it would pin nothing: the two could drift apart in
// silence, which is exactly what a diff round caught here. The guard against lowering it is a RUNTIME
// test (`felitronics_core_tests`, "the gainToDb floor"), not a static_assert: `std::log10` is not
// portably constexpr in C++20 — it is not on the MSVC row — so a compile-time assertion could only pin
// the constant and never that the function still uses it.
inline constexpr double kGainToDbFloor = 1.0e-12;
inline double gainToDb (double gain) noexcept { return 20.0 * std::log10 (gain > kGainToDbFloor ? gain : kGainToDbFloor); }

// Fast 20*log10 for DETECTOR paths — accurate to ~0.001 dB and several times cheaper than
// std::log10, which is enough for deciding how hard to compress and nowhere near enough for
// measurement. Use gainToDb() for anything a user reads as a number.
//
// A float is m * 2^e with m in [1,2), both free from its bit pattern, so ln(x) = e*ln2 + ln(m) and
// only ln(m) needs work. On that interval the atanh series in t = (m-1)/(m+1) converges fast — t is
// at most 1/3, so four terms are already past float precision. No magic minimax constants, nothing
// to mistype: the series is its own proof.
inline float fastGainToDb (float gain) noexcept
{
    if (! (gain > 1.0e-9f)) return -180.0f;          // also rejects NaN, denormals and zero
    const std::uint32_t bits = std::bit_cast<std::uint32_t> (gain);
    const int   e = (int) ((bits >> 23) & 0xFFu) - 127;
    const float m = std::bit_cast<float> ((bits & 0x007FFFFFu) | 0x3F800000u);   // mantissa in [1,2)
    const float t  = (m - 1.0f) / (m + 1.0f);
    const float t2 = t * t;
    const float lnM = 2.0f * t * (1.0f + t2 * (0.333333333f + t2 * (0.2f + t2 * 0.142857143f)));
    return 8.685889638f * ((float) e * 0.693147181f + lnM);   // 20/ln(10) * ln(x)
}

// An INTENTIONAL exact floating-point `==`: documents the intent at the call site and silences
// -Wfloat-equal in exactly this one place, so core headers stay warning-clean under a downstream
// consumer's strict flags (JUCE's recommended set fires it). JUCE-free analogue of juce::exactlyEqual.
template <typename T>
constexpr bool exactlyEqual (T a, T b) noexcept
{
#if defined(__clang__)
    #pragma clang diagnostic push
    #pragma clang diagnostic ignored "-Wfloat-equal"
#elif defined(__GNUC__)
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wfloat-equal"
#endif
    return a == b;
#if defined(__clang__)
    #pragma clang diagnostic pop
#elif defined(__GNUC__)
    #pragma GCC diagnostic pop
#endif
}

// DID AN AUTONOMOUS STEP LEAVE THE STATE UNCHANGED, BIT FOR BIT? The silence advancement of law 11c is a
// deterministic, memoryless map `S -> F(S)`: if one step maps a state to ITSELF, every later step does
// too, so the rest of the pause is a no-op and can be skipped — exactly, not approximately.
// `==` will not do as the test. It calls +0.0 and -0.0 equal, and a follower whose state is -0.0 is one
// step away from the +0.0 the sum `0.0f + (-0.0f)` produces; and it calls a NaN unequal to itself, which
// would spin the loop for the whole pause on a state that cannot move. Bits answer both correctly.
inline bool sameBits (float a, float b) noexcept
{
    return std::bit_cast<std::uint32_t> (a) == std::bit_cast<std::uint32_t> (b);
}

} // namespace felitronics::core

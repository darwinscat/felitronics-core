// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/core/DetMath.h>   // felitronics::core::det — the deterministic dB spelling below

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
// test — `felitronics_pause_is_silence_tests`, the group "law 11c premise — gainToDb's floor is where
// the collapse thinks it is" — and not a static_assert: `std::log10` is not portably constexpr in C++20
// (it is not on the MSVC row), so a compile-time assertion could only pin the constant and never that
// the function still uses it.
inline constexpr double kGainToDbFloor = 1.0e-12;

// AND THE CLAMP ITSELF IS NOW WRITTEN ONCE, for the reason the paragraph above gives about the constant.
// There are two spellings of this conversion below — the system one, whose bits are a shipped product's
// sound, and the deterministic one, whose bits are compared across rows — and a floor RETYPED into the
// second could drift from the first in silence, which is the exact failure the paragraph above describes
// one level up. Law 11c's collapse is pinned to this predicate by a runtime test; both spellings inherit
// that proof rather than each needing its own.
namespace detail
{
    inline double gainToDbFloor (double gain) noexcept { return gain > kGainToDbFloor ? gain : kGainToDbFloor; }
}

inline double gainToDb (double gain) noexcept { return 20.0 * std::log10 (detail::gainToDbFloor (gain)); }

// THE SAME CONVERSION FOR A NUMBER THAT LEAVES THIS MACHINE. `gainToDb` above is `std::log10`, and that is
// not one function across rows: measured over the ranges this library uses, Apple's and glibc's log10
// disagree at 4325 of 200000 points and glibc's and musl's at 4278. A reported dB that moves by an ulp is
// a printed digit that moves, and the analyzers' outputs are diffed BYTE FOR BYTE between the native CLI
// and the wasm module — so a value destined for a report goes through `det::log10`, one implementation
// compiled into every build, and a value destined for the audio path does not.
//
// WHY THIS IS A SECOND FUNCTION AND NOT A POLICY ON THE FIRST. `gainToDb` is called ONCE PER SAMPLE on
// four real paths — dynamics/GainReductionPath.h:173, deesser/DeEsser.h:184, dynamiceq/DynamicEqBand.h:169
// and, worst, limiter/TruePeakLimiter.h:594, which is inside the OVERSAMPLED loop of the module that is
// most of a render's cost. `det::log10` is ~2.9x a system call, so routing the shared function through a
// policy would tax every one of those to make a handful of reported values reproducible. The split is by
// CONSUMER, not by function, and the lint (tools/lint/check-det-math.mjs) is what keeps it that way.
inline double gainToDbDet (double gain) noexcept { return 20.0 * det::log10 (detail::gainToDbFloor (gain)); }

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

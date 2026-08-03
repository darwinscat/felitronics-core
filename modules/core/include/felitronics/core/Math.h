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
inline double gainToDb (double gain) noexcept { return 20.0 * std::log10 (gain > 1.0e-12 ? gain : 1.0e-12); }

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

} // namespace felitronics::core

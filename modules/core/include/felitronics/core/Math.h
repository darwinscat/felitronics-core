// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <cmath>

namespace felitronics::core
{

constexpr double kPi = 3.14159265358979323846;

// dB <-> linear amplitude (20·log10). `double` on purpose: these are offline / coefficient-design /
// GUI helpers, never the per-sample loop (Law 3 carve-out). The floor keeps log10 finite.
inline double dbToGain (double dB)   noexcept { return std::pow (10.0, dB / 20.0); }
inline double gainToDb (double gain) noexcept { return 20.0 * std::log10 (gain > 1.0e-12 ? gain : 1.0e-12); }

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

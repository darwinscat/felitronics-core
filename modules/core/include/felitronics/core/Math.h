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

} // namespace felitronics::core

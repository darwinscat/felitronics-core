// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa. Part of felitronics-core — see LICENSE.

#pragma once

// Transitional compat shim — re-exports the eq types + the kMaxChannels SSOT (now in felitronics::core)
// into teq::.
#include <felitronics/eq/EqTypes.h>

namespace teq
{
    using felitronics::eq::FilterType;
    using felitronics::eq::BandParams;
    using felitronics::core::kMaxChannels;
}

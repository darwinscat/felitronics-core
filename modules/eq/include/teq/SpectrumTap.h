// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

// Transitional compat shim — re-exports felitronics::analysis::SpectrumTap (+ the FFT-size contract)
// as teq::SpectrumTap. teq::SpectrumTap is now the default-order alias of the templated tap.
#include <felitronics/analysis/SpectrumTap.h>

namespace teq
{
    using felitronics::analysis::SpectrumTap;
    using felitronics::analysis::kSpectrumFftOrder;
    using felitronics::analysis::kSpectrumFftSize;
}

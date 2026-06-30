// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

// Transitional compat shim — the entry point. `#include <teq/EqEngine.h>` brings the whole teq:: API
// (EqEngine + BandParams/EqBand/FilterType/Svf/Smoother/SpectrumTap/matched/…) via the sibling shims,
// re-exported from felitronics::eq / felitronics::core / felitronics::analysis.
#include <felitronics/eq/EqEngine.h>
#include <teq/EqBand.h>
#include <teq/SpectrumTap.h>

namespace teq { using felitronics::eq::EqEngine; }

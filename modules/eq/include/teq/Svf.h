// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

// Transitional compat shim — re-exports the Cytomic SVF as teq::Svf.
#include <felitronics/eq/Svf.h>
#include <teq/EqTypes.h>   // teq::FilterType, teq::kMaxChannels

namespace teq { using felitronics::eq::Svf; }

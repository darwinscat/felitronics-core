// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

// Transitional compat shim — re-exports felitronics::core::Smoother (exp) + LinearSmoother as teq::.
#include <felitronics/core/Smoother.h>

namespace teq { using felitronics::core::Smoother; using felitronics::core::LinearSmoother; }

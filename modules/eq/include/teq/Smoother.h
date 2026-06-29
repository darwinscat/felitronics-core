// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa. Part of felitronics-core — see LICENSE.

#pragma once

// Transitional compat shim — re-exports felitronics::core::Smoother as teq::Smoother.
#include <felitronics/core/Smoother.h>

namespace teq { using felitronics::core::Smoother; }

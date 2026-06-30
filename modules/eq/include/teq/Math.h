// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

// Transitional compat shim: re-exports felitronics::core into the legacy `teq::` namespace so existing
// consumers (`#include <teq/Math.h>`; `teq::kPi`) compile unchanged while they repoint FetchContent.
#include <felitronics/core/Math.h>

namespace teq { using felitronics::core::kPi; }

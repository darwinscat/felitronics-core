// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

//==============================================================================
// felitronics::blend — umbrella header for the mic-blend engine: parameter model + kernels + the engine
// itself. OFFLINE, JUCE-free, headless-testable (double coefficient design, float sample buffers).
// Promoted from OrbitCapture's app-local mixer (model/MixModel.h + core/BlendKernels.h + core/BlendEngine.h).
//==============================================================================

#include <felitronics/blend/BlendParams.h>
#include <felitronics/blend/BlendKernels.h>
#include <felitronics/blend/IrBlend.h>

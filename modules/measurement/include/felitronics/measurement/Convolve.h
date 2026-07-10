// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

//==============================================================================
// felitronics::measurement — OFFLINE, MESSAGE-THREAD-ONLY double-precision convolution for IR
// measurement (allocates; NOT RT-safe by design). This is deliberately NOT the RT FFT: real-time
// consumers use the float `felitronics::core::fft` SEAM.
//
// As of v0.7.0 the double FFT floor (nextPow2 / detail::fftInplace / convolve / magSpectrum) lives in
// `felitronics::core::offline` (a SECOND offline-double-FFT consumer — the analysis display curves —
// appeared, which was the documented PROMOTION TRIGGER). This header now just RE-EXPORTS them into
// `felitronics::measurement` so the whole measurement API + all its tests are unchanged. The precision
// rationale (why double, why the -208 dB floor is enough) now lives in OfflineFft.h.
//==============================================================================

#include <felitronics/core/OfflineFft.h>

namespace felitronics::measurement
{
// Re-export the promoted offline FFT floor — measurement's public spelling is unchanged.
using core::offline::nextPow2;
using core::offline::convolve;
using core::offline::magSpectrum;
} // namespace felitronics::measurement

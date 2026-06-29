// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa. Part of felitronics-core — see LICENSE.

#pragma once

// Transitional compat shim — re-exports the matched-biquad design + runtime into teq:: (including the
// `matched`, `rbj`, `detail` sub-namespaces the tests reference).
#include <felitronics/eq/MatchedBiquad.h>
#include <teq/Math.h>   // teq::kPi

namespace teq
{
    using felitronics::eq::BiquadCoeffs;
    using felitronics::eq::Biquad;
    namespace detail  { using namespace felitronics::eq::detail; }
    namespace matched { using namespace felitronics::eq::matched; }
    namespace rbj     { using namespace felitronics::eq::rbj; }
}

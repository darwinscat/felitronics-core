// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

// Transitional compat shim — re-exports the band engine + free response helpers into teq:: (and pulls
// the teq:: shims for everything <teq/EqBand.h> transitively provided in the original).
#include <felitronics/eq/EqBand.h>
#include <teq/EqTypes.h>
#include <teq/Math.h>
#include <teq/MatchedBiquad.h>
#include <teq/Smoother.h>
#include <teq/Svf.h>

namespace teq
{
    using felitronics::eq::EqBand;
    using felitronics::eq::BandDesign;
    using felitronics::eq::designBand;
    using felitronics::eq::evalCoeffs;
    using felitronics::eq::bandResponse;
    using felitronics::eq::laneView;
    using felitronics::eq::Axis;
    using felitronics::eq::axisLane;
    using felitronics::eq::compositeResponse;
    using felitronics::eq::ResponseMatrix;
    using felitronics::eq::matrixResponse;
}

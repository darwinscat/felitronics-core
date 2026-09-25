// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

//==================================================================================================
// LAW 11c ("a pause is silence") — the shared fixture helpers of PauseIsSilenceTests: the gap lengths,
// the bit comparisons, the drain to exact rest and the charging tone. Moved out of
// modules/mastering/tests/PauseIsSilenceTests.cpp verbatim, so felitronics-guitar-core tests its
// address of the law (poweramp::PowerAmpStage) with the SAME fixtures instead of a copy that could drift.
//
// Include in ONE translation unit per executable, then `using namespace felitronics::test::law11c;`.
//==================================================================================================
#pragma once

#include <felitronics/core/Math.h>

#include <cmath>
#include <cstring>
#include <vector>

namespace felitronics::test::law11c {

constexpr double kFs = 48000.0;

inline bool bitsEqual (float a, float b) noexcept { return core::sameBits (a, b); }
inline bool bitsEqual (double a, double b) noexcept { return std::memcmp (&a, &b, sizeof (double)) == 0; }

// The gap lengths. Deliberately NOT multiples of 16 (DynamicEqBand's control period), of 64 (the block
// sizes the older gap tests use) or of the lookahead: a length aligned to the mechanism's own period is
// the third recorded form of a blind fixture, and the counter shortcut in DynamicEqBand is exactly what
// such a length would hide.
const int kGaps[] = { 1, 2, 7, 15, 16, 17, 31, 63, 64, 65, 127, 480, 1000, 4801, 48000 };

// The collapse's own list adds the interval BETWEEN the two events that bound it: with a 5 ms Rms window
// the detector level crosses `kGainToDbFloor` at about 13 200 samples and the recurrence parks at about
// 23 609, and phase 2 is the stretch in between. `kGaps` jumps 4 801 -> 48 000 straight over it, so a
// silent loop that stopped one sample early was invisible: in Peak mode the level reaches zero in a
// single step, and by 48 000 both runs are long since parked on the same word.
const int kGapsCollapse[] = { 1, 2, 7, 15, 16, 17, 31, 63, 64, 65, 127, 480, 1000, 4801,
                              14000, 16000, 20000, 23000, 48000 };

// BRING A STAGE'S PER-CHANNEL PATH TO EXACT REST by feeding it real digital silence until its output is
// bit-zero on every sample of a whole block. It has to be MEASURED, not counted: the fixture's first
// version drained a fixed four blocks, which is enough for a band that CUTS and not for the same band
// BOOSTING — a +18 dB bell rings down about x80 per block, so it needs nine, and at four the residue was
// 1.1e-5 and produced a 1.4e-7 difference on the return that looked exactly like a defect in the code.
// Returns false if rest was not reached inside `maxBlocks`, and every caller asserts that.
template <class Stage>
bool drainToRest (Stage& a, Stage& b, int blockSize, int maxBlocks = 64)
{
    for (int k = 0; k < maxBlocks; ++k)
    {
        std::vector<float> z1 ((std::size_t) blockSize, 0.0f), z2 ((std::size_t) blockSize, 0.0f);
        float* io[2] = { z1.data(), z2.data() };
        std::vector<float> w1 ((std::size_t) blockSize, 0.0f), w2 ((std::size_t) blockSize, 0.0f);
        float* jo[2] = { w1.data(), w2.data() };
        if (! (a.process (io, 2, blockSize) && b.process (jo, 2, blockSize))) return false;
        bool rest = true;
        for (int i = 0; i < blockSize; ++i)
            rest = rest && bitsEqual (z1[(std::size_t) i], 0.0f) && bitsEqual (z2[(std::size_t) i], 0.0f)
                        && bitsEqual (w1[(std::size_t) i], 0.0f) && bitsEqual (w2[(std::size_t) i], 0.0f);
        if (rest && k > 0) return true;          // k > 0: one more block AFTER the first silent one
    }
    return false;
}

// A tone that charges a detector: full scale, so the shared state has somewhere to travel FROM.
inline void fillTone (std::vector<float>& v, int startSample, double hz, float amp)
{
    for (std::size_t i = 0; i < v.size(); ++i)
        v[i] = amp * (float) std::sin (2.0 * core::kPi * hz * (double) (startSample + (int) i) / kFs);
}

} // namespace felitronics::test::law11c

// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa. Part of felitronics-core — see LICENSE.

#pragma once

namespace felitronics::stereo
{

//==============================================================================
// felitronics::stereo::MidSide — the L/R ⟷ Mid/Side matrix. M = ½(L+R) is the mono sum (what a mono
// fold or a vinyl lateral cut hears); S = ½(L-R) is the stereo difference (the part that goes to zero as
// the image narrows). The ½ convention makes encode∘decode an identity. Processing the two axes
// independently is the basis of width / mono-bass / elliptical tools.
struct MidSide
{
    static inline void encode (float l, float r, float& m, float& s) noexcept { m = 0.5f * (l + r); s = 0.5f * (l - r); }
    static inline void decode (float m, float s, float& l, float& r) noexcept { l = m + s; r = m - s; }
};

} // namespace felitronics::stereo

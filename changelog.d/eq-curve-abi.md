<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

### tools — the EQ's curve through the C ABI; v7

**`fc_master_eq_curve (params, sampleRate, lane, band, freqHz, count, outDb, cap, written)`.** The magnitude
response, in dB, of a parameter set's EQ at frequencies the caller names — `eq::EqEngine::magnitudeDbFor`, read
out. No handle and no render: the parameters travel with the call, so the curve is answerable while a knob is
moving, and the call asks the heap for nothing.

`lane` is an `fc_eq_axis` (`FC_EQ_AXIS_STEREO`, `_LEFT`, `_RIGHT`, `_MID`, `_SIDE` — `eq::Axis`, where the four
domain axes each fold the Stereo lane in and `STEREO` is that lane alone); a code that names no axis is
`FC_ERR_ENUM`. `band` is −1 for the whole bank, or 0..`FC_MAX_EQ_BANDS`−1 for one band, anything else
`FC_ERR_RANGE`. `sampleRate` is held to the `eq` module's own domain — `FC_ERR_NON_FINITE` for a non-finite
rate, `FC_ERR_REFUSED_BY_CORE` for one outside it. The buffer is the caller's and `cap` is binding: the call
writes all `count` values or none, so `cap` below `count` is `FC_ERR_CAPACITY` and `count == 0` is
`FC_ERR_RANGE`; a non-finite frequency anywhere in the grid is `FC_ERR_NON_FINITE` before anything is written.
`freqHz` and `outDb` may not touch and `written` may not point into `outDb` (`FC_ERR_SPAN`); `written` is left
as it was by every refusal. `params.bypassEq` and each band's `dyn` are not read — the curve is the static
response of the bands.

**C ABI v7.** One entry point and no struct, so no row moves in the size table and every struct keeps its v6
size. `fc-master-layout.mjs` is v7.

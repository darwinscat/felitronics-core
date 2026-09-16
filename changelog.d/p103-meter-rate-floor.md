<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

### analysis — the loudness meter refuses a rate below 8000 Hz itself, and a refused prepare forgets the last programme

**What `analysis::LoudnessMeter` now promises about the rate:** a rate of 8000 Hz and up
(`LoudnessMeter::kMinSampleRate`, the core's floor) is measured; a rate in (0, 8000) is refused by
`prepare()`, `prepareForSamples()` and `storageFor()` alike; a rate that is **not given** — zero, negative,
NaN — still reads as 48 kHz, exactly as published. +inf is refused, as it was. Since P51 every entry that
hands the meter a rate from outside already refused these rates; now a direct C++ caller gets the same answer
instead of a K-weighting filter past Nyquist (a 0 dBFS 400 Hz sine read +3043 LUFS at 3300 Hz).
`KWeightingFilter::prepare()` cannot refuse and still takes any rate; the meter is its only owner.

**A refused prepare now disarms the meter completely** — a change on the OLD refusals too (a channel count
outside 1…16, a capacity that cannot be represented), not only on the new one. Before, a refusal cleared only
the prepared flag, and `momentaryLufs()`, `shortTermLufs()`, `integratedLufs()`, `loudnessRangeLu()`,
`droppedBlocks()`, `nonFiniteSubHops()` and `gatingBlockEnergies()` went on answering for the previous
programme. Now every one of them answers what a never-prepared meter answers.

Nothing changes for a successful prepare: a verdict-and-readings oracle over 3663 rates is bit-identical
outside (0, 8000) Hz, including zero, negative, NaN and +inf. The conformance table that sized the store at
150, 149 and 100 Hz now does it at 8050, 8049 and 8000 Hz with the same block counts.

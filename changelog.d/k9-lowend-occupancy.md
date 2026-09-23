### analysis · tools — LowEnd reaches 20 Hz, says where it stops resolving, and counts how OFTEN a band is there

The band table started at 30 Hz and answered one question: how loud is each semitone on average. A sub that
plays on an eighth of a programme is 10·log10(1/8) = 9.03 dB down in that average against one that plays
throughout, so a threshold set for the second misses the first entirely. K9 adds the other question.

**The range goes down to E0** — `lowNoteHz` is 20.0, 47 bands, MIDI 16..62. Raising `fftOrder` to 18 to
resolve the new bottom was costed and **refused**: it doubles the window to 5.46 s, which halves the frame
count a duty rests on (130 frames in three minutes becomes 64) and smears exactly the intermittency
occupancy exists to find. One request item would have eaten another. Memory would have risen 78 %; CPU per
second of programme is not the argument, at +6 %.

So the limit is published rather than hidden: `firstResolvedBand()`, `resolvedAboveHz()` — the closed form
`kLobeBins·binHz/(kSemiUp−kSemiDown)`, 25.356 Hz at 48 kHz/order 17 — and `lobeBins()` so no caller
re-types the 4. **This was already needed:** at 96 kHz the old 30 Hz table has nine under-resolved bands up
to 49.0 Hz, and nothing said where the limit was.

**Occupancy.** `dutyThresholdDb` (20 dB, a parameter rather than a constant — it is a hypothesis a consumer
will calibrate), counted once per accepted frame. Published: `dutyFrames()` and `dutyCount(b)` raw, so a
consumer never defines the population itself; `duty(b)`; `levelWhenOnDb(b)` — how loud the band is against
its frame's loudest **when it is on**, which duty does not carry and a filter choice needs; and
`lowestOccupiedBand(dutyMin)` as an accessor, because that threshold gets swept during calibration.

**Duty is a count of FRAMES, not a fraction of time**, and it is documented as such at the accessor. The
frames overlap, so one event can mark two, and at 48 kHz/order 17 the hop is 1.365 s: an event on "an
eighth of the time" does not read 0.125.

**The gate is the design, and the control is built from the failure it prevents.** "Within the threshold of
the frame's loudest band" is `E >= max·q`, which on digital silence reads `0 >= 0` — true of every band.
Planted back, it puts 15 silent frames into the denominator and marks all 47 bands present in each: 705
occupancy claims on nothing. A frame now enters only if its band total clears the note floor, the same
share of `frameEnergy()` the whole-programme report already uses. `count > 0` in `lowestOccupiedBand` is
load-bearing for the same reason: at `dutyMin <= 0` a band present in no frame would come back as
"occupied 0 % of the time", a number that reads as a finding.

**A crossover sweep, free, from the raw table.** The spectral axes are fed Mid and Side *before* the
crossover, so any LR4 low-pass can be applied to the band table afterwards as a weight:
`sideFractionBelow(fc)`. Against a real filter it agrees to **4.7e-5** across seven crossovers from 60 to
300 Hz — the design round had bounded it at 1e-2 — and it is bit-identical whichever crossover the build
installed. That accuracy is why no bank of three secondary filters was built: the sweep answers "where
should it sit" for the whole curve, and what still needs a real filter is the one frequency chosen, which
`fc_probe_lowend_run_with` now provides at any value instead of three fixed ones.

**And the sweep's limit is asserted as a failure, not hedged in a comment.** The table covers
[lowNoteHz, highNoteHz] and nothing else, so anti-phase energy outside it is invisible to the sweep while
the real low-pass passes it. Add anti-phase 15 Hz and 700 Hz and at fc = 60 Hz the sweep reads 0.0036 where
the installed filter reads 0.488 — the opposite answer to a question asked at 0.25, and the 15 Hz term is
precisely the sub-20 Hz vertical hazard this work is for. Rule: sweep to choose, install and read
`lowSideFraction()` for the number.

**`skipBlocks`**, default 0, moves the **histogram and nothing else** — the series, the integrals and the
extrema stay complete, because `worstFractionBlock()` is a coordinate into a programme and a coordinate
over a silently shortened programme is a lie. Pinned: with `skipBlocks = 5` the histogram is exactly the
first five blocks' samples lighter, every block of the series is bit-identical, the integral is
bit-identical, and the worst-block index is unchanged. `skippedBlocks()` is published, because a population
whose size a caller cannot see is one it cannot divide by.

**`settlingBlocks(sampleRate, crossoverHz, dB)`**, static and pure, so that skip comes from the filter
instead of being written down. LR4's double pole gives an envelope `t·exp(−t/τ)` with `τ = √2/(2π·fc)`;
normalised to its own peak that is `u·exp(1−u)`, and `u` is 10.233 at −60 dB. At 48 kHz: 2 blocks at
120 Hz, 3 at 100, 12 at 20. The "first 10 blocks" this replaces appears nowhere in the header it was read
from, and is five times too many at 120 Hz and too few at 20.

**Surface.** `fc_probe_lowend_run_with` (crossover, note range, fftOrder, duty threshold, skipBlocks) with
`storage_bytes_with` beside it, because law 11d's budget must answer about the geometry the call will use.
`side_fraction_below(hz)` and `lowest_occupied(dutyMin)`. Scalars 60 → 67, band row 11 → 13 (`dutyCount`,
`levelWhenOnDb`), and **both format roads go v1 → v2 together** — the JS stride check is the version gate,
and leaving it behind is how a parity contract drifts in silence.

**Callers must address bands by `midi` or `centreHz`, never by index.** Band 0 was MIDI 23 and is now
MIDI 16: every position moved by seven semitones. Four oracles in this repo's own suite spelled MIDI 23 as
"band 0" and all four went red; they now derive the first note from the range. Two more in the ABI suite
had stopped being oracles by spelling a width — a band stride of 11 that sized a 46-row buffer for a
13-wide row, and a 64-double buffer whose 60-double ask read lowend's correct refusal as "this mode has no
result".

**A measurement that corrected itself,** recorded because the two readings have opposite consequences. The
narrowest band (3.25 bins) reads about 1 dB light in density under noise — 24 realisations, mean −0.985 dB,
sd 1.343. On a *deterministic* flat spectrum the same ratio is exact to 1e-6, so the geometry and the
fractional edge rule are clean: what shifts is the log of a noisy estimate over few degrees of freedom, and
it shrinks as frames accumulate. Read `lowestOccupiedBand` on the longest excerpt available, and do not
tighten `dutyMin` on short clips to compensate for something that is not a bias.

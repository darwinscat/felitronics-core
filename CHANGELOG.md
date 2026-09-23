<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

# Changelog

Notable changes to felitronics-core. Releases are git tags (`vX.Y.Z`); the project VERSION lives in
`CMakeLists.txt`.

## v0.44.0 — 2026-09-23

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

## v0.43.0 — 2026-09-23

### analysis · mastering — every short-term series is sampled at 10 Hz, as EBU Tech 3342 requires

`LoudnessMeter` took one 3 s short-term sample a second. Tech 3342 §3.1 asks for more: *"A minimum block
overlap of 2.9 s between consecutive analysis windows (i.e. ≥10 Hz sampling of the loudness level) is
required."* With a 3 s window a 2.9 s overlap **is** a 100 ms step, so the two phrasings are one requirement,
and 100 ms is this meter's hop (exactly, at every rate whose hundredth is a whole number — see the last
paragraph). The requirement entered in V3 (January 2016); the meter had kept libebur128's 1 Hz cadence, a
2 s overlap. It now takes one sample per hop.

**This is compliance with the clause's text, not the repair of a failed conformance test.** Tech 3342's own
Table 1 passes at either cadence — 10.000 / 5.000 / 20.000 / 15.000 LU, well inside the ±1 LU tolerance —
and that table is now in the suite as `tech3342MinimumRequirements()`, cases 1–4, with 5–6 named as not
covered rather than left to be inferred.

**What it changes.** Not "nothing on long programmes", which is what a first measurement on three fixtures
suggested and is wrong. On a square envelope whose states last exactly the window's 3 s, the range reads
20.0 LU at 1 Hz and 9.5 LU at 10 Hz, and the gap does not close with length: 1 Hz answers a flat 20.0 at 12,
20, 30, 60 and 120 s, and 10 Hz answers 9.5 at all of them but 20 s, which is 9.3 — the same figures at 44.1
and 48 kHz. The mechanism is the percentiles. At 1 Hz such a programme offers the window six phases, a sixth
of them sitting wholly inside a loud state and a sixth wholly inside a quiet one,
so P95 and P10 both land on an extreme and the answer is the envelope's full 20 dB swing. At 10 Hz there are
sixty phases, each aligned one is 1.7 % — under the 5 % P95 reaches for — and the fifty-odd mixed windows
that 1 Hz never looked at fill the distribution between them. On that fixture the two cadences agree exactly
away from the window's timescale — 0 LU at 1.5 s states, 20 LU at 5 s and 7.5 s — but that is the fixture, not
a regime: its states are 20 dB apart, so both clear the −20 LU relative gate. Drop the quiet state to −30 dB
and 5 s states read 4.7 LU against 7.7, 7.5 s read 4.7 against 6.9, 10 s read 4.7 against 5.7.

**Cost.** Ten times the short-term entries: for a ten-minute programme the store is 48 KB where it was 4.8.
Every allocation oracle that spells a meter moved with it and says why — the solver's per-pass figure, the
`measure_lra` and `solve` budgets in the C ABI, and the storage table in the conformance suite.

**New:** `LoudnessMeter::shortTermCount()` — the number of short-term samples a range was read from. The
cadence is part of what the standard specifies, and a count is the only way a caller, or a test, can see
which cadence a build runs.

**Fixed, a defect this change would have introduced.** The accepted programme length was bounded at
`INT_MAX - 4`, sized for the block store's `floor(hops) + 4`; the short-term store's `ceil(hops) + 8` then
reaches `INT_MAX + 4`, and the guard `stCount < (int) stE.size()` compares against a narrowed negative, so
every short-term sample is refused from the first one on — a range of 0.0 LU presented as a measurement,
with `droppedBlocks()` reading 0. The bound is now `INT_MAX - 8`, the tighter of the two stores, pinned by
an allocation-free check at the exact boundary. (Unreachable in practice — the store would be 17 GB — and
reachable on paper, which is where a bound lives.)

**New, and a defect that predates this change:** `LoudnessMeter::droppedShortTermSamples()`. The range is
read from its own store, which overflows on its own schedule, and nothing counted the overflow — a programme
fed past its declared length got a range computed over whatever survived, with no way to know. The margin
also narrowed with the cadence: at one sample a second it took an overfeed of about a hundred hops; at one a
hop, 38 are enough. `droppedBlocks()` documented itself as covering the range too, which is half right in a
way worth stating: the block store always fills FIRST, 30 hops earlier, so a caller that checks it reads 0
does have both answers whole — but a non-zero count does **not** mean the range was truncated, and for three
seconds of overfeed it climbs while the range is still perfect. The ordering is now pinned as an ordering —
blocks always first, at both rates, at whole hops AND at a fractional one, and stated in hops rather than in
seconds, because the first draft of that test asked for "3.8 s" and got a length one hop short of what it
meant.

Two things the change made visible rather than caused:

- The C ABI's `measure_lra` fixture swung 1.9 LU on paper but measured 0.5–0.6: its states were 1 s long, so
  the 3 s window never sat inside one and always averaged them. The precondition asking for more than 0.5 LU
  was passing by a single 0.1 LU histogram bin, and the cadence spent that bin. The states are now 4 s, the
  measured range is the programme's real 1.9 LU at **both** cadences. What needs that range is not the check
  it was written under — a bit-equality between the ABI's answer and the core's catches an added offset even
  on a range of zero — but the other defect the block guards against, one that reads channel 0 for every
  plane: channel 0 is the steady tone, so that reads 0.0 LU against a truth of 1.9.
- The LRA-constrained solver search closes this fixture in seven renders where it closed in six. The bracket
  is judged against the measured input range, and that range moved. But the SEARCH still takes six: the
  seventh render is the delivery re-render, bit-identical to the fourth. Before the change the sixth probe
  happened to keep the limit and was delivered as it stood; now it breaks it by a 0.05 LU quantum — the range
  reads on a 0.1 LU staircase and the 2.05 LU limit sits between the two readings that probe can land on. The
  bound is therefore not loosened to seven: the re-render is identified and subtracted, so the six the test
  used to carry (regula falsi, not bisection) stays a search count, and the stale-end halving that makes it
  regula falsi is now asserted directly from the log.

**`ProgrammeReport` moves with it — the same clause governs it.** That class computes EBU Tech 3342 from its
own short-term series rather than reading `LoudnessMeter::loudnessRangeLu()`, so that the gated observation
count is known exactly, and it sampled that series once a second. For the length of this work the core
therefore held **two numbers under one standard's name that disagreed by 10.5 LU** — on the 3 s-state
envelope, at 44.1 and 48 kHz, over 12, 30 and 120 s programmes, `lraLu` read 20.0 where the meter read 9.5.
`kObservationHops` is now 10, both series sample at 10 Hz, and they agree to 0.3 LU again — the tolerance
being summation order (this class oldest-first, the meter newest-first) and nothing else.

Six published fields move with that cadence: `lraLu`, `shortTermP10`, `shortTermP50`, `shortTermP95`,
`shortTermSpreadLu` and `shortTermObservations` — the last by a factor of ten, since it counts observations. A consumer calibrated
against them will see the change; that is the point, not a side effect.

**The guard that should have caught it, and why it did not.** A check asserting the two classes agree to
0.3 LU already existed and passed the whole time. Its fixture is long steady steps, where every window that
is not straddling a step reads the same value and the sampling cadence cannot matter — a guard naming the
right pair and blind by fixture. It is now held on a square envelope with 3 s states as well, which is the
shape that can tell two cadences apart, and the independent oracle beside it writes out BOTH numbers (300
sub-hops for the window, 10 for the step) rather than importing either from the header.

**One boundary moved with it, worth naming.** `ProgrammeReport` refuses with `LoudnessCapacityExceeded`
past `maxDurationSec` — but the store carries a fixed margin behind that line, and the margin is counted in
OBSERVATIONS, so ten times as many observations is a tenth as much slack in seconds. Measured against a 10 s
declaration at 48 kHz: the refusal used to turn at 20 s of programme and now turns at 13.8. The contract did
not change; the undeclared slack behind it shrank to match its own words, from about a hundred hops to 38. A
caller that fed a few seconds over its declaration and got away with it will now be refused. Pinned in hops
at both rates, which nothing did before — the check that stood there was a disjunction ("either fitted or
says it did not") and asserted nothing either way.

**Budgets, for anything that reads them.** `fc_master_need(FC_NEED_MEASURE_LRA)` at the ABI suite's geometry
is **3296 B** where it was 2936; `fc_master_need(FC_NEED_SOLVE)` for 1 s of stereo at the default 1000
buckets is **1 047 848 B** where it was 1 047 704. Struct layouts and the ABI version are untouched — these
are computed, not frozen. `ProgrammeReport`'s short-term store at the default one-hour ceiling is 288 KB
where it was 28.9.

**A compliance caveat the work surfaced, deliberately left alone.** A sub-hop is `lround(0.01·fs)` samples,
so the hop is 100 ms only where a hundredth of the rate is whole. At 48 and 44.1 kHz the cadence is 10 Hz
exactly; at 22050 Hz it is 9.9774 Hz — *under* §3.1's minimum, by 0.23 % — and at 8050 Hz 9.9383 Hz, with the window
rounding to 3.0068 s and 3.0186 s against a specified 3 s. The rounding cuts the other way too: 8049 Hz takes
80.49 down to 80, giving 10.0613 Hz and a 2.9817 s window — fast enough, and too short. That is the meter's
whole sub-hop grid, older than this change, and moving it would move every number the class produces at
those rates. Recorded as its own item rather than touched here: the rates that ship are 44.1, 48 and 96 kHz,
where the cadence is 10 Hz exactly.

## v0.42.0 — 2026-09-22

### analysis · probe — BandCrest: the crest a band lost between a source and its master

The micro-dynamics damage counter an automatic mastering mode needs. Macro damage is LRA loss and already
exists; timbre is a band energy shift and does not yet; this is the third. Run it on the source, run it on the
delivered master, and the damage is the crest a band lost, block by block.

**What a band crest can and cannot mean** is in the header, because the honest limits are part of the
instrument: it is the peak-to-RMS of a signal seen through a NAMED filter, never a statement about "the content
between 2 and 6 kHz". The band is a dome and an LR4 skirt is only 24 dB/octave down, so a band's peak can be a
neighbour's transient leaking in. In the DIFFERENCE of two runs through the same filter that leakage largely
cancels; in one run's absolute number it does not.

#### Four decisions, each where the obvious choice is wrong

**Reconstructed peak, not sample peak — and the reason is not the sample rate.** Even with both runs on the
same grid a sample peak is biased in ONE direction: the input carries inter-sample peaks and the output has
been through a true-peak limiter whose whole job is to remove them. A sample-peak crest would under-read the
input's crest and therefore UNDER-REPORT the loss, which is the one number this exists to produce. The
reconstruction is the certifying instrument's own 4x/32, so the full-band peak can be nulled against the
certificate rather than believed.

**One reconstruction, then the split at 4x** — cheaper than five reconstructions, and it makes the band
magnitudes nearly rate-independent as a side effect, since an LR4 corner is prewarped by `tan(pi*f/fs)` and at
the base rate the skirt moves with the sample rate.

**The bands are independent filters of the input, not a reconstructing tree.** A tree is cheaper and sums back
to the input, but pays for that with allpass compensation on every band — and an allpass changes the waveform,
which is to say the peak, which is to say the very quantity being measured. Nothing here needs the bands to
sum, so nothing here pays for summing. They do not reconstruct the input and are not meant to.

**The activity mask is built from the SOURCE only.** This is the decision that can quietly invalidate
everything else: a mask taken from the output, or from each side separately, can drop exactly the blocks the
chain damaged most — a loud block the chain crushed can fall below a gate its input passed, and the damage
leaves the population with it. The gate is a conjunction, because each half fails in its own direction: an
absolute programme floor (a relative loudness gate would throw away the quiet passages where a compressor
works) AND a band share floor (without it a band with no content contributes the crest of its own noise floor,
and a source that is digital black in a band against a master carrying dither reads as several dB of "damage").

#### Four numbers, not one — the obvious one is blind to the instrument doing the damage

A limiter works on the loudest few percent of blocks. The suite's own fixture has **p95 reading 0.000 dB while
CVaR95 reads 2.214 and 8 of 297 blocks lost more than 3 dB**: a cost function calibrated on p95 would call that
master undamaged. `cvar95Db` — the mean of the worst 5 % — is the number to read; `p95Db` stays beside it
because its disagreement with CVaR95 is itself the signal that damage is concentrated; the exceedance counts
say the same as a count; and `p5Db` is the negative side, where a limiter's pumping shows (it lowers a band's
RMS around a kick while the peak elsewhere is untouched, so the crest RISES and the positive part hides it).

`peakShiftDb` and `levelShiftDb` ride along so a crest loss with no peak loss — a floor RISE, a clipper putting
a steady harmonic into a band that held only a click — can be told apart from a flattening.

The requester asked for the check in the other direction too: a fixture where p95 sees damage and CVaR95 does
not must not exist. It cannot, by construction — CVaR95 averages the worst 5 %, every member of which is at or
above the 95th percentile — but "cannot by construction" is exactly the kind of statement that is false at the
edges of a rounding rule, so it is checked over 400 shapes including populations of one, all-zero, all-equal
and a single outlier.

#### Every oracle is outside the object, and two of them caught defects before any test was written

The full-band peak is **bit-identical** to `ReferenceTruePeakMeter` at 44.1, 48 and 96 kHz — the same topology,
so a tolerance would hide a difference rather than allow one. The block grid equals `LoudnessMeter`'s gating
block count at all three. The interpolator's delay is measured from an impulse: symmetric about oversampled
index 319.5 with a worst difference of exactly 0, which is where the 63.5 in the header comes from — a HALF
sample, named rather than rounded away.

Two independent hop clocks, counted by different code over the same time, **disagreed 31 against 30**: the
interpolator's drain was being taken whole, so a hop of silence past the end of the programme closed as
programme. And the RT check passed on clang and failed 25 times on gcc — `alloc::count.load()` was one argument
to `ok` and an allocating `std::string` was another, and argument evaluation order is unspecified.

#### The surface

`fc_probe_crest_run[_with](slot, ...)`, `_scalars`, `_blocks`, `_loss`, `_storage_bytes[_with]`, plus
`fcore_measure crest [--against other.f32]`. **Two slots, not one**, because the loss is arithmetic that must
not be reimplemented across the ABI: formed from linear cells it is one logarithm and invariant to the gain a
master has over its source, and reassembled in JavaScript from two dB columns it would be neither. The
per-band accepted-block counts are published rather than left to be counted from the block table — a count a
consumer derives itself is a second definition of "active" waiting to drift.

An unrecognised `--` token is a refusal on both roads from the start. Native and wasm are byte-identical on the
default road, on a parameterised one and on a run carrying the loss rows; their refusal sets agree on twelve
argument lists, ten of them refusals.

#### What two review rounds found, and what that says about the tests

A code review found seven things and a falsification round — whose mandate was to make this analyzer report a
wrong number rather than to read it — found nine more. Both are worth listing, because six of the sixteen were
cases where a TEST was narrower than the claim it stood under.

**The full band's crest was a hybrid and wrong near Nyquist.** Its peak is the reconstruction; its mean square
came from the oversampled stream, which the 32-tap interpolator droops above about 0.375·fs. A pure sine has a
crest of exactly 3.0103 dB — an oracle that owes this code nothing — and it read **9.03 dB at 0.45·fs and
23.8 dB at 0.49·fs**. The per-band crests were never affected, because their peak and their mean square come
from the same stream; only that one number was wrong, and it fed the absolute gate and the published programme
level with it. The full band's mean square is the base-rate one now and the share gate keeps the oversampled
denominator, because a ratio across two domains is a share of nothing.

**The refusal set differed by platform.** Six sibling analyzers carry a 768 kHz ceiling; this one had none, so
`llround` decided — and it saturates on arm64 while it wraps on x86-64 glibc. Measured on both rows:
`fs = 1e308` was refused on one and accepted on the other with a 10-sample hop, answering a measurement that
looked entirely valid. The bounds are comparisons made before any conversion now, and the two rows agree.

**The filters were `SystemMath`**, where every sibling is deterministic and the math-policy suite asserts it of
each — `std::tan` disagrees between libms and moved 12 block rows of 3752 on a 6.3-minute mix. **The hop was
`llround(fs·hopMs/1000)`** where the loudness meter builds it from sub-hops, so 22050 Hz gave 2205 samples
against 2210 — and the test that claimed the identity ran at 44.1, 48 and 96 kHz, the three rates where the two
agree. **A refused `prepare()` answered the previous run.** **The comparator compared anything**: the same audio
with the master on a 50 ms hop gave CVaR95 7.8 dB with `valid` true. **No evidence was reported as a lag**: a
muted master read "misaligned by 8 blocks". **The non-finite count meant interpolator outputs**, 128 for one
NaN, and its location depended on the call sizes. **A master scaled by 2⁻⁵⁰** — float subnormals, strictly
positive — was compared as real. **A channel that stopped was not cleared**, so a returning one spliced a
32-sample-old ring onto new audio: 40 dB of crest on a block. **And the budget was a tautology**: `bytes() > 0`
passed while the demand omitted the interpolators themselves, 168 B a channel at every length.

**The bed placed every transient on the hop grid.** At t = k·0.5 s, which is a hop boundary at every rate this
suite uses — the feature sat exactly on the grid the measurement is made of, in the fixture the alignment tests
are built on. That is this repository's own rule 2 inverted. Moved off it by a third of a hop, which turned the
CVaR95 control red, because that control's threshold was a number fitted to the fixture of the day rather than
a statement of the mechanism. It states the mechanism now.

Also measured and then gated: gain invariance is **exact only where the scaling is exact** — ×0.5 is, ×0.37 is
not, and the 1.1e-5 dB that appears there is the fixture rounding rather than the comparison drifting, so the
claim is now the two claims it always was. And the filter bank runs the twelve evaluations it reads instead of
the twenty it computed: **2.84–3.12 s → 2.04–2.05 s** on a 4-minute stereo programme, bit-identical — but only
after a diff said it was not, because the first version passed `kQ` cast to float where `Crossover2` passes a
double.

## v0.41.0 — 2026-09-22

### mastering · ABI — the limiter's gain reduction over the windows it had something to work on

`fc_gr_stats` answers what a LIMIT is judged on, and it has to: every 4 ms window of the programme is an entry,
the silent ones included, so a caller cannot buy headroom with silence. That makes it the wrong answer to the
other question with the same units — how hard the stage works WHERE it works. On the same distribution, 20 % of
silence turns the music's p95 into its p93.75, and enough of it drags the quantile into the silent mass and
zeroes the statistic outright. Two questions, two distributions. The first one does not move: the solver's own
constraint still reads it, and K11 changes no number any limit is judged by.

**`fc_solution_gr_active_stats(solution, stage, out)`** (ABI v10) answers the second, for the limiter.
`fc_loudness_request::limiterActiveInputDb` is the gate, in dBFS at the limiter's node, default −60.

**The gate is on the stage's INPUT, never on its gain reduction.** Gating on |GR| would define "where the stage
works" as "where the stage worked" and report a statistic of a set the statistic itself chose. For the limiter
the input is already measured at the right place and on the right clock: `limiterPeakLin` is the reconstructed
peak the limiter SAW, one per oversampled sample, at the same index as `limiterGrDb`. Nothing new is tapped and
nothing new is aligned — the estimate that this would need a new tap and half a sprint was wrong by inspection.

**A whole window is accepted or dropped**, never part of one: the window is the quantile's unit and half a
window is an entry whose mean is over a denominator nobody stated. The decision is the window's PEAK input, not
its mean — a limiter reacts to peaks, so a window holding one transient in an otherwise quiet stretch is a
window it worked in. And every field is over the accepted windows, the samples' mean, max and active fraction
too, not only the quantiles: a mean over the whole programme beside a quantile over part of it would be two
bases under one name.

`windows` and `activeWindows` are published beside the numbers, so a caller can see how much of the programme
the answer rests on. A gate that accepts nothing answers `valid = 0` with the counts, not a zero: "no windows
qualified" and "the stage never worked" are different statements and only one of them was measured. The
compressor is **FC_ERR_STATE** for the same reason — its input is not tapped, and answering it with zeroes
would be a measurement nobody made.

**What the test had to be fixed to show.** The first fixture targeted −14 LUFS, where this material's limiter
touches only the peaks: 95 % of windows held no reduction and BOTH figures read 0.000 dB — a fixture that could
not tell the two distributions apart while looking like a measurement. The second version made the limiter work
but let the solver re-converge on each fixture, so a moved statistic could not be told from a moved denominator.
The third renders both at ONE fixed drive (`maxPasses = 1`, every drive-bound limit off — Q2's corner), so the
music is rendered identically in both and the only difference is the silence appended to one. The gated p95 is
then bit-identical between them and the ungated one falls, which is the whole claim.

The control's threshold was also wrong on its own terms: it demanded a 25 % collapse and failed where the
mechanism had worked perfectly. How far an ungated p95 falls is a fact about the distribution's shape near its
top — with half the programme silent it becomes the music's own p90, and a dense distribution moves little. The
check now claims the direction, which is what is true.

The suite PRINTS the four figures, because a passing check prints nothing and a number quoted from a green run
would then be one nobody could reproduce — which had already happened once, an earlier fixture's numbers
reaching a report describing this one. As it stands: ungated 6.0850 → 5.4950 dB, gated 6.0850 → 6.0850 dB.

What is invariant is the GATED QUANTILE and the maximum, not every field: the seam window between the music and
the silence passes the gate (the limiter's release carries into it), so the accepted-window count and the mean
move by that one window. "The same music at the same drive" is a statement about the population the quantile is
read from, and the test asserts exactly that and nothing wider.

**Two things a review round found afterwards, and one of them is a hole of a class already closed once.**
`fc_gr_active_stats` was missing from the size oracle on BOTH sides — the native `headered[]` table and the
JavaScript `STRUCT_IDS` — so `layout-check` passed while neither half knew the struct existed. Exactly the shape
of the flag hole fixed a commit earlier: when both halves are silent, the diff between them says nothing. Both
now carry it, so the row has an oracle and the page's `Struct` can stamp it.

The other was a **contradiction between the code and three doc sites about a NaN gate**. The code accepts
NOTHING for a NaN (it falls through both tests and lands on +inf); the header, the facade comment and the domain
row all claimed the opposite. The code is what was kept, and deliberately: a NaN is a caller's mistake, and
taking it as the widest gate hands back a full set of statistics that look like a measurement at a gate nobody
chose, while taking it as the narrowest hands back zero active windows, `valid = 0` and the NaN itself echoed in
`thresholdDb` — three signals a reader cannot miss. The three doc sites now say that, and a test pins the
MEANING rather than the admissibility: the domains gate can only report that a non-finite value crossed without
a refusal, not which end of the range it landed on.

**The budget grew by one histogram**, and every allocation oracle in three suites went red at once — which is
the oracle doing its job. Two of them were also carrying numbers their own run no longer produced (a "727 696 B"
in a sentence beside an assertion parameterised by the histogram count); those now print what they measured. One
precondition pinned three literals, two of which were functions of things that test is not about; it now pins
the per-pass cost and spells the rest through the expressions the budget is made of.

ABI v10 — one new struct with its own id and row, one entry point, one field at the end of `fc_loudness_request`
(120 → 152 at v10). `fc_gr_stats` is nested by value and frozen with it; nothing existing moved.

### analysis · dynamiceq · ABI — the measurement surface the automatic mastering mode asks for: a band that is an argument, the milliseconds behind a deviation knob, and the alignment stated as a contract

Three items requested by darwinscat.com's `mastering-v2`, which decides the shape of a chain from measurements
of the input and then has to compare its output against that input. Each one closes a question the surface
could not answer, and two of them close a way it could answer WRONGLY.

**K3 — `fc_probe_bursts_run_with`: the same detector, pointed somewhere else.** `analysis::BandBursts` at its
documented 5–9 kHz band answers "sibilance". The same machinery at 80 Hz – 8 kHz with a shorter baseline
answers "how dense are the transients", which decides whether a clipper goes before a limiter or after — and
that is a different question, not a different analyzer. The new entry point takes `bandLowHz`, `bandHighHz`,
`hopMs`, `baselineMs`, `enterDb`, `exitDb`; `maxEvents` is deliberately not exposed, because the counters keep
counting past the event list and `eventsComplete` already says whether the list is whole, so a caller that
could shrink it could only make the list lie about itself. `fc_probe_bursts_storage_bytes_with` prices it: the
baseline ring is `round(baselineMs / hopMs)` hops of `hopMs` each, so both of those move the allocation and a
page sizing itself by the default figure would be short exactly where it asked for a longer memory.

`fcore_measure bursts` takes the same six values as flags, with the same names, defaults and refusals. That is
not a convenience: a road the parity harness cannot drive is a road with no gate, and `fc_probe_bursts_run_with`
would otherwise have shipped with its only evidence being that it compiles. Native and wasm are byte-identical
on the parameterised road, and their refusal sets agree on nine malformed or out-of-domain argument lists.

**A number that was true and was about to become a lie.** The scalar block reported `enterDb` and `exitDb` by
reading the compile-time default, which was correct while one road existed and would have described 6 / 3 dB
while the events came from the caller's numbers the moment a second one did. The module now remembers what the
last successful run installed and reports that.

**`BandBursts::onsetsPerSecond()`** — the whole of the question is the DENOMINATOR, which is why it lives in the
analyzer rather than in each consumer's arithmetic. It is the JUDGED programme, `eligibleHops()` long, not the
file: the first `baselineHops` have no surroundings to be measured against and are not judged, so dividing by
the file's length reports a density over a stretch where no onset could have been found, and reads low on
exactly the short programmes where it matters most.

**K7 — `fc_master_eq_dyn_times`: the milliseconds behind a deviation knob.** `fc_eq_dyn::atk` and `::rel` are
not times. They are deviations in `[0, 1]` around an automatic value the core derives from the band's own fc/Q
(0.5 is auto, 0 four times faster, 1 four times slower), and that automatic value is a field of no struct — so
a caller could not tell "50 ms" from "the slow rail" by looking at what it had written. The call answers what
the follower is actually set to.

It is **per lane, not per band**, which is the core's shape rather than a choice made here: the times come from
the sidechain probe, the probe sits on the lane, and a point's up-to-five lanes carry their own freq/Q while
`dyn` is shared. And it runs `dynamiceq::LaneDynamics::ballisticsFor`, the one expression the chain runs. That
function is new, and it exists because the rails are part of the answer: a lane's freq and Q reach the producers
RAW (`MasteringChain` hands them the caller's parameters, not the band's clamped copy) and are railed inside
`LaneDynamics` to `[10 Hz, 0.49·fs]` and `[0.05, 40]`. A readback that had asked `BandBallistics` with the raw
pair would have answered for 5 Hz where the probe sits at 10. Both audio call sites now go through it too, so
there is one place where the rails live.

A rate below `core::kMinSampleRate` is refused rather than substituted, unlike the audio path, which reads an
unusable rate as 48000 because a `prepare` has already refused one by then. Nothing has refused anything here,
and an answer computed at a rate the caller did not ask about is indistinguishable from one it did.

**K6 — the delivered render is aligned with its input, and the rest is time, not index.** No code changed: what
was missing was the claim, stated where a consumer reads it and pinned where it can fail. `OfflineRenderer`'s
contract is arithmetic — `out[n] = y[n + D]` — so the latency is already off, output sample n is input sample n
processed, and the chain's tail is in the output rather than cut off; a caller that subtracts `latencySamples`
is introducing an error rather than removing one. The other half is not latency at all: on a delivering handle
the two sides have different rates, different lengths and different sample grids, so windows belong in TIME and
an index-for-index comparison measures the resampling ratio and calls it drift.

Three instruments pin it, deliberately of different construction, because one instrument agreeing with itself is
what the first version of this check did while it was wrong: bit-for-bit identity with every stage bypassed,
correlation lag with the stages running, and where a lone impulse lands. The energy centroid was tried as the
second instrument and rejected on evidence — it disagreed by 3.2 ms with the stages running, and that was the
limiter telling the truth (it attenuates ahead of the peak through its lookahead and stays down through its
release, so more of a burst's tail is pulled down than its head). A centroid measures the envelope a dynamics
stage exists to reshape, so it can only answer an alignment question with every stage bypassed, which is where
it is used and where it is exact.

**ABI v9** — one entry point, no struct grew, and therefore no row in the size table, exactly as v7.

## v0.40.0 — 2026-09-21

<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

### mastering · tools — the EQ points of the mastering chain are DYNAMIC: `eqBands[].dyn` renders

**No ABI change.** `FC_MASTER_ABI_VERSION` stays 8 and no struct moves. The six `fc_eq_dyn` fields were
already carried, clamped and handed to the core; what they were not was rendered. `eq::EqBand` applies a
gain delta a PRODUCER pushes into its per-lane seam, and the chain's EQ stage had none — so a point with
`dyn.on = 1` was, measurably, the same samples as one with it off.

**`mastering::MasteringChain` now drives one `dynamiceq::LaneDynamics` per EQ band.** The detector key is the
EQ section's INPUT, captured through `eq::EqEngine::captureSectionInput` before any band moves the signal, so
one point's moving delta cannot modulate a later point's detector at an overlapping frequency. Each point's
five placement lanes get their own probe, programme estimate and delta from one shared setting, which is the
composition `dynamiceq` was written for; nothing new is computed here.

**The default render did not move.** With every point unarmed — `dyn.on = 0`, which is the default — the EQ
stage is `eq::EqEngine::process` over the same quantum, sample for sample, and a producer whose point is not
armed runs its band over the WHOLE quantum in one call rather than in control-rate chunks. Checked against
v0.39.0 through the C ABI on twenty seconds of real programme with three points, mono-bass, compressor,
clipper, limiter and 24-bit dither: byte for byte identical, including a run with every `dyn` field set to a
far end while `dyn.on` stayed 0. `dyn.rangeDb = 0` is no dynamics whatever `dyn.on` says, and renders
identically too.

**It costs no latency.** The producer applies the delta derived from the PREVIOUS control chunk, so nothing is
read before it is written: `latencySamples()`, `compressorTapOffset` and `limiterTapOffset` are unchanged, and
an impulse leaves at the same index armed or not. Block invariance survives with a point armed — the
producer's 16-sample control grid restarts at every call, and the chain's call is always the internal quantum,
so the restart lands on the same absolute sample at every caller block size (measured bit-identical at 1, 7,
63, 256, 1021, 4096 and 65536).

**And the loudness search is undisturbed**, which is the property the EQ's place in the chain buys: the stage
sits ahead of the node the search moves, so the PRE-LIMITER tap is bit-identical at every drive with a point
armed, the limiter's worst gain reduction stays monotone in that drive, and `y(g, c) = 10^(c/20)·y(d, 0)`
holds to 4.9e-07 at a peak of 0.891 on a 3x3 grid — the same float-rounding order `LoudnessSolver.h` states
it at for a chain with no EQ at all.

**Memory is declared before it is asked for**, as every other stage's is: `Storage::dynBands` counts the
producers, `bytes()` carries them, and `prepareBytes()` is still what a fresh `prepare()` allocates byte for
byte on every row of the rate x width x topology matrix. One array of 24, held on the heap beside the engine
for the same reason — about 60 KiB on a chain that has an EQ, nothing on one that does not. `process()`
allocates nothing with points armed, arming and disarming included.

**The four `dyn` domain rows are pinned by behaviour now.** `FC_DOMAINS` carried `rangeDb`, `thrDb`, `atk` and
`rel` as clamps with no read-back — they were the last four rows of the table that nothing could pin, and the
reason was that the fields were inert. They now read back through a render, and `MasterDomainsTests` asserts
the opposite of what it used to: an armed point is NOT the same samples as `dyn.on = 0`. Pinning them needed
armings that make each bound visible — a threshold every sample is over for the `rangeDb` cap, an input gain
that moves the programme to meet each end of `thrDb`, and a threshold inside the programme's own band level
for the ballistics, since a reduction pinned at its cap never releases and a release knob then changes
nothing. Zero of 89 rows are now unpinned.

**`fcore_master` can drive it**: `band<i>.dyn.on`, `.dyn.range`, `.dyn.thr`, `.dyn.thrAuto`, `.dyn.atk`,
`.dyn.rel`.

**Not yet accepted by ear, and shipping that way on purpose.** Oleh listened to the A/B on 2026-09-21 — it
works, but a dynamic point is hard to hear on programme material, so this feature is lightly tested where
only ears can test it. That listening continues in TabbyEQ, where the same dynamics are still experimental;
until it concludes, treat the feature as one the measurements vouch for and the ear does not yet.

## v0.39.0 — 2026-09-21

<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

### tools — every enum code and every field domain of the mastering ABI, published once and held against real calls

**No ABI change.** `FC_MASTER_ABI_VERSION` stays 8 and no struct moves: what follows is metadata a page reads
out of `tools/wasm/fc-master-layout.mjs`, beside the layouts that already live there.

**Ten enum lists, all of them.** `FC_FILTER_TYPE`, `FC_DETECTOR`, `FC_LINK_MODE`, `FC_COMP_MODE`, `FC_SHAPE`,
`FC_NOISE_SHAPING`, `FC_GR_STATISTIC` (with `Percentile`), `FC_GR_STAGE` and `FC_PROGRESS_STAGE` join
`FC_EQ_AXIS`, each an array whose INDEX is the code, and `FC_ENUMS` says which header enum each one mirrors.
`layout-check.mjs` reads the enumerators out of `tools/fc_master_abi.h` and compares every list entry by entry
— declaration order, value and letters — so a code added, renumbered or permuted in the header fails there
instead of silently renaming a menu. Consumers held five of these by hand and had just drifted on the sixth.

**`FC_DOMAINS` — 89 rows, one per input field** of `fc_master_config`, `fc_master_params` (the EQ bands, the
mono bass, the compressor, the clipper, the limiter and the dither) and `fc_loudness_request`. Each row carries
the unit, the admitted interval, what a value outside it does — `refuse` with the FC_ERR_* that names it,
`clamp`, `verdict` (`fc_master_solve` answers FC_OK and the SUMMARY carries `FC_SOLVE_INVALID_REQUEST`), `free`
where the code bounds nothing and only finiteness is checked, or `any` — what a non-finite value does, where
the applied value can be read back, and what the domain depends on. Bounds that move with the sample rate are
written as `0.49*sr` / `8000/sr` and evaluated by `domainBound (bound, sampleRate)`; the four that depend on
another FIELD say so rather than averaging (the mono-bass stage narrows the width to exactly 2; the oversample
factor stops at 16 with the limiter, 64 with the clipper alone and nowhere with neither; the compressor's
250 ms lookahead ceiling is that stage's own). **Silent clamps are marked**, because the difference decides an
interface: a refusal arrives on the call that made it, a clamp arrives as nothing at all. No field of this ABI
has a list-valued domain — `oversampleFactor` admits every integer in its interval and `dither.bits` refuses
nothing — and the table says so instead of inventing one. Bounds that move with the oversample factor as well
as the rate are written `0.49*sr*os`, and `domainBound (bound, sampleRate, oversampleFactor)` evaluates all
three forms. `layout-check.mjs` holds the table against the struct layouts in both directions: a row must name
a value field that exists, and EVERY input field of the three structs must have a row, so a field added in a
later version cannot arrive without a domain.

**`FC_CONSTRAINT_BITS` and `constraintsOf (mask)`** — `fc_solution_summary.alsoViolated` is a bitmask whose bit
`i` is `FC_CONSTRAINT[i + 1]`, with `binding` included in it. The list is derived from `FC_CONSTRAINT` rather
than written out, and `layout-check.mjs` holds the core's `constraintBit()` expression against what it assumes;
`fc_constraint` and `fc_solve_status` join the enums pinned to the header, and a reordering of the C++
`MasteringConstraint` behind them is already a build error in the facade's static asserts.

**Two rows the first draft got wrong, both found by review and both reproduced through the ABI before they were
changed.** `initialGainDb` said "no bound"; the search clamps the starting gain to the chain's own ±60 dB
before its first render, and at `maxPasses = 1` a request of 100 comes back as 60 in the summary — so the row
is a clamp now, read back through `summary.preLimiterGainDb`. `clipper.dcBlockHz` said it had no ceiling; the
corner is clamped to `0.49*sr*os`, measured as 100000 Hz applying as 94080 at 48 kHz and 4×.

**`felitronics_master_domains_tests`** holds the table against the running ABI with a real C-ABI call per
bound, at two sample rates for every rate-dependent row. The table is the test's ARGUMENT and every probe is
derived from its own numbers, so a bound moved there without the code moving with it lands on the wrong side of
the real boundary. A clamp is pinned by READING BACK what was applied — `fc_master_resolved` where it publishes
one, `fc_master_eq_curve` for the EQ lane fields, the rendered audio for the rest — and asserting that the value
at the bound and beyond it are the same number while one step inside it is a different one; acceptance alone
would pass against a clamp anywhere at all. A refusal is probed twice — one step out and WELL out — because a
bound moved INWARDS still refuses one step past itself, and that near probe alone could not tell a delivery
floor of 8000 Hz from one of 22050.

**The step is the RESOLUTION OF WHAT ANSWERS, not a fixed fraction.** A whole number steps by 1; a status is
sharp, so a refusal steps by a billionth; a read-back is not, and each field states its own (`stepRel` /
`stepAbs`) — which is how tightly that bound is pinned. Where a bound is unreachable at ordinary settings the
field also states the parameter set that reaches it: the compressor's 400 dB range cap needs a threshold far
under the programme, the time constants need a step above the point where the ballistics stop being
distinguishable from instant, and the dither has to be off wherever the quantiser step is coarser than the
probe.

**The four `eqBands[].dyn` rows cannot be pinned here, and the reason is now a checked fact.** `eq::EqBand`
applies a delta that a producer pushes in through `setLaneDeltaDb`, and the engine the mastering chain drives
has none — so those fields are carried and clamped and change nothing this ABI renders. A named check renders
with the band dynamics armed and with `dyn.on` cleared and requires the two to be identical; if a producer is
ever wired in, that check goes red and the four rows become pinnable like any other clamp. They are the whole
of the list the run prints of clamps nothing can pin.

1280 checks; a planted-mutation round over the table killed 56 of 59, and all three survivors are a row DEMOTED
to a weaker claim, which the suite's header names as what it does not catch. Sixteen more planted violations —
an enum permuted in the header, a list shortened, a row naming a field that no longer exists, an input field
with no row, `constraintBit()` rewritten — are all caught by `layout-check.mjs`.

<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

### mastering · tools — the gain-reduction percentile, a general quantile read-back, and a bound that is a guarantee; C ABI v8

**A gain-reduction quantile is a WINDOW statistic now — a behaviour change, and `p95Db` moves with it.**
`GainReductionStats::p95Db` used to be the 0.95 quantile over tap SAMPLES; it is now the 0.95 quantile over the
programme cut into fixed 4 ms windows, each window contributing the MEAN `|GR|` over it (`GainReductionSummariser`,
`kGrQuantileWindowSeconds`). A single click is averaged down inside its window and no longer moves a percentile; a
sustained reduction does. The last window is short and is averaged over its own length, so it is an entry like any
other, and EVERY window of the programme is an entry — the silent ones included, because the denominator of "the p95 of
the gain reduction" is the programme and not the part of it the stage worked in. The window is fixed in the core and is
not a request field: it is part of what the number means, and two callers with two windows would be comparing different
quantities under one name. `GainReductionStats::aboveRange` therefore counts WINDOWS past the histogram's top.
`meanDb`, `maxDb` and `activeFraction` are unchanged sample statistics, and `GrStatistic::Max` stays sample-wise.

**`GrStatistic::Percentile`.** A fourth statistic beside `Mean`, `P95` and `Max`, with its fraction in
`GainReductionLimit::quantile` — in (0, 1] and finite, else `InvalidRequest` before any pass, whatever the statistic.
The default is 0.95, so a `Percentile` limit at its default is the `P95` limit, bit for bit. `grStatisticValue` is the
one function that reads a statistic off a summary. `GainReductionStats::quantileDb` is the number the limit was
judged by — NaN where the distribution cannot answer, which is "an unanswerable statistic is not a violation" in
arithmetic rather than in a second branch.

**`LoudnessSolution::grQuantile (stage, q, outDb)`.** The q-quantile of a stage's `|GR|`, read off the very
distribution the limits were judged on, so a reading and its limit are the same number and a reading can say how close
a render came. Each solution owns its two distributions (`compressorGrWindows`, `limiterGrWindows`), written by every
render like the traces, so a later solve cannot move a number already reported. `solveBytes (…, binDb = 0.01)` and
`DeliveredMastering::solveBytes (…, binDb)` count them: 640 016 B at the default bin width, making a 1 s stereo solve
727 696 B against 87 680 B before.

**`TargetUnreachable` now delivers a render that HOLDS the constraint it names.** `DriveBound` can only guarantee a
limit that grows with drive once it has seen a render holding one — it brackets `ok`, the loudest render that broke
nothing, under `cap`, the quietest that broke something. A search that started past the boundary and ended there never
got an `ok`, and what came back was the gentlest BROKEN render carrying the broken limit's name. The search now spends
ONE render at the drive the limiter idles at (`kIdleDriveMarginDb` under the engagement point, where the reduction is
zero and zero holds any reduction limit), and then walks the bracket that render completes with whatever budget the
search did not spend — `DriveBound::probe`, the same one the ordinary search uses — so what comes back is the LOUDEST
render that holds the limit rather than the quietest. Measured on this tree, 0.40 to 0.81 LU above the idle point; on
the three mixes this was built for, 0.4 to 1.8. It is all taken after the search, so a solve that finds its own
holding render is unchanged render for render and bit for bit — the idle render included, which is why its ceiling
moves only where the margin cannot absorb the meters' disagreement. With the budget exhausted the idle render is
still what comes back: the one case where a quiet render is delivered on purpose, because the guarantee outranks the
loudness.
- `pairFor` is the one place a drive becomes a `(gain, ceiling)` pair: `d = g - c` and the two are clamped to ±60 dB
  one number at a time, so a ceiling picked for the true-peak aim alone put the gain past its clamp and the drive
  RENDERED was not the drive chosen. The ceiling is chosen for the drive, inside the window that keeps the gain in
  range and at or under the promise; where that window is empty no pair expresses the drive and the rescue is not
  taken. Where the window forces the ceiling up, the delivered peak can pass the promise — that is reported rather
  than hidden (below).
- The certifying meter does not read the peak the limiter aims at; `DriveBound::overshootAt` measures the
  difference on the renders the search made (0.02 dB on this tree's music, 0.25 dB on a 15 kHz tone). Every PROBE
  subtracts it from its ceiling, as the ordinary search does. The IDLE render subtracts it only where it exceeds
  `truePeakAimDb`, the margin that exists to absorb it: below that size the ceiling is the aim exactly, because a
  ceiling moved by a fraction of a margin that already covers it moves the render — and a `Solved` that render
  reached becomes a `TargetUnreachable` a hair outside the tolerance.
- A probe that breaks a limit the bracket is about is a FAILED probe: it moves `cap` and is not offered as the
  render to deliver. `Best` ranks infeasible candidates by the worst excess across all constraints, so a reduction
  broken by a millionth of a decibel ranked gentler than a peak broken by a tenth, and the render handed back broke
  the very limit the verdict named.
- `alsoViolated` now carries what the DELIVERED render breaks as well as what stopped the search. Where nothing
  holds every constraint the render handed back breaks something of its own, and a caller told only why the search
  stopped was not told its file is above the promise.
- "The loudest render that holds the limit" is the contract only where the loudness is a measurement: renders under
  the absolute gate all read the meter's -120 sentinel, tie, and the first one offered is kept.
- A render taken aside need not have a measurable loudness: the reduction comes off the tap and the peak off the peak
  meter, so a render under the absolute gate still holds the limit and is still the answer, though never `Solved`.
  `violatedMask` and `worstExcess` no longer judge the peak-to-loudness ratio without one.
- None of these renders enters `Best::nearest*`, so `binding` is read from `cap` — what was broken at the smallest
  drive that broke anything, which the walk has just tightened — and the rest go into `alsoViolated`.
The `TargetUnreachable` line in the header says all this instead of "the best FEASIBLE render".

**C ABI v8.** `fc_loudness_request` gains `limiterGrQuantile` and `compressorGrQuantile` (144 B), defaulting to 0.95;
at their defaults a call is v7's. `FC_GR_PERCENTILE = 3` joins `fc_gr_statistic`. `fc_solution_gr_quantile (s, stage,
q, outDb)` answers the same distribution at any `q` — one entry point rather than a field per fraction — refusing with
`FC_ERR_ENUM` for a stage code that names nothing, `FC_ERR_NON_FINITE` for a non-finite `q`, `FC_ERR_RANGE` for a
finite `q` outside (0, 1] and `FC_ERR_REFUSED_BY_CORE` where the distribution cannot answer, and writing `*outDb` only
on `FC_OK`. `fc-master-layout.mjs` is v8; `fcore_master` takes `limGrQ=`, `compGrQ=` and `percentile` for
`limGrStat=` / `compGrStat=`.

## v0.38.0 — 2026-09-18

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
No two of `freqHz`, `outDb` and `written` may touch — all three pairs `FC_ERR_SPAN`, `written` inside the grid
included, since the grid is the caller's `const`; `written` is left as it was by every refusal.
`params.bypassEq` and each band's `dyn` are not read, and neither can refuse the call: the curve is the static
response of the bands, so a `dyn` the caller never filled is not this call's business.

**C ABI v7.** One entry point and no struct, so no row moves in the size table and every struct keeps its v6
size. `fc-master-layout.mjs` is v7 and exports `FC_EQ_AXIS` — the lane names in the order `fc_eq_axis` declares
them, so the index is the code `lane` takes and a consumer needs no hand-written copy of that order.
`layout-check.mjs` reads the `FC_EQ_AXIS_*` codes out of `tools/fc_master_abi.h` and holds the array against
them entry by entry — declaration order, value and spelling — so a permutation in the header fails the check.

## v0.37.0 — 2026-09-18

<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

### mastering — the loudness search stops at the limit it may not break

**`mastering::TargetLoudnessSolver`**, for the three limits that grow with drive — the limiter's gain reduction, the
peak-to-loudness ratio and the loudness-range loss:

- **A target past such a limit is answered at the limit.** When the next step asks for a drive at or past the
  smallest one a render broke a limit at, the search renders inside the bracket between that render and the loudest
  one that kept every limit instead — regula falsi on the limit's own statistic, held inside the bracket, the
  midpoint where a statistic is not measured — and delivers the loudest render that keeps them as
  `TargetUnreachable`, `binding` what the smallest breaking drive broke. It stops once the breaking render, at the ceiling its true peak asks for, is under the
  target's tolerance and within `toleranceLu` of that render, or at `maxPasses`. It used to step on loudness alone and
  deliver the feasible render it happened to have: after a first render at the starting gain, the source at its own
  level. Nothing is assumed about how a limit moves with drive beyond where to look: every render the bound chooses is
  measured, and a render that keeps a limit above one that broke it ends the bracket.
- **Whenever the search without the limit solves on a render that keeps it, the search with the limit is the same
  search**, render for render and bit for bit. Without a render that kept every limit there is no bracket, and the
  search is the one before.
- **The first render whose limiter works, after an idle one, is aimed 0.15 dB under the aim**, or by the largest
  overshoot measured on a working render at or under its drive. It used to be aimed at the promise itself, land over
  it by the between-sample overshoot and spend a render on the correction: on the suite's programme a -12 LUFS target
  now takes two renders where it took three.

No field and no ABI change.

<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

### limiter · mastering · tools — the limiter's dual release, the caller's gain-reduction trace size; C ABI v6

**`limiter::TruePeakLimiter` — dual release.** `TruePeakLimiterParams::dualRelease` (off) and `slowReleaseMs` (200 ms).
On, a second envelope with instant attack and release `slowReleaseMs` takes as input the least, over the last 70 ms
(`kSlowWindowMs`), of the largest reduction required within 20 ms (`kSlowBridgeMs`); the applied reduction is the larger
of the two envelopes', and `releaseMs` is the fast one's. Both envelopes read only the required reduction, so the
loudness search's scale law and the reduction's monotonicity in drive hold. Off, the output is the single release bit
for bit, whatever `slowReleaseMs` holds. Switching on starts the slow envelope and its windows from 0 dB; switching off
continues the fast envelope from the applied reduction. `effectiveSlowReleaseMs()` (0 while off) and
`MasteringChainResolved::limiterSlowReleaseMs`. Both windows are prepared whether or not it is on
(`Storage::slowWindow`, `slowBridge`).

**`mastering::TargetLoudnessSolver` — the trace's size is the caller's.** `LoudnessRequest::grTraceBuckets`: 1000 by
default, 1..65536, otherwise `InvalidRequest` before any pass. Each trace holds min(grTraceBuckets, frames) buckets on
the heap (`GainReductionTrace::bucket`, `kDefaultBuckets`, `kMaxBuckets`); a bucket's `samples` and `nonFinite` are
64-bit. `solveBytes (fs, nch, frames, grTraceBuckets)` and `DeliveredMastering::solveBytes (…, grTraceBuckets)` count
both traces. At 1000 buckets the trace is the one before, bit for bit.

**C ABI v6.** `fc_master_params` gains `limiterDualRelease`, `_pad0`, `limiterSlowReleaseMs` (6584 B);
`fc_master_resolved` gains `limiterSlowReleaseMs` (96 B); `fc_loudness_request` gains `grTraceBuckets`, `_pad0`
(128 B). `fc_solution_gr_trace64` copies `fc_gr_trace_bucket64`; `fc_solution_gr_trace` answers `FC_ERR_RANGE` for a
count past 32 bits. `fc_master_need_solve (h, req, frames, out)` budgets a solve for a request; `FC_NEED_SOLVE` budgets
1000 buckets. `fc-master-layout.mjs` is v6; `fcore_master` takes `lim.dual=`, `lim.slowRelease=`, `grTraceBuckets=`.

<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

### mastering · tools — progress and cancellation of the whole-programme calls; C ABI v5

**`felitronics/mastering/Progress.h`.** `ProgressCallback { bool (*fn) (void* context, const ProgressEvent&); void* context; }`,
taken by new overloads of `TargetLoudnessSolver::solve` / `measureInputLoudnessRange` and `DeliveredMastering::solve` /
`measureInputLoudnessRange`. An event carries the stage (`Convert`, `LoudnessRange`, `SearchPass`, `FinalRender`), the
render's number and its bound, the fraction of the stage — exactly 0 first, exactly 1 last, no two events more than 1 %
of the programme's frames apart — and, on a render's last event, its log record. `true` continues, `false` stops:
`solve` answers the new `MasteringSolveStatus::Cancelled`, the range measurements `false`, and the same call made again
gives the same result. The overloads without a callback are unchanged, and the results are the same bits either way.

**C ABI v5.** `fc_master_set_progress (h, fc_progress_fn fn, void* context)` sets a handle's callback; it receives an
`fc_progress` (stage, pass, maxPasses, fraction, hasRecord, record) and returns 0 to stop, and the call then answers
`FC_ERR_CANCELLED` (15). On the wasm tier a handle without one calls `Module.onProgress(msg)` when that is a function;
`false` or an exception stops. No struct with a header grew. Checked on a real 5:21 stereo programme in node by
`tools/wasm/progress-check.mjs`.

**The render, too.** `DeliveredMastering::render` gained the same `progress` overload as its siblings, reporting
`Convert` then a new `Render` stage — no log record, since there is no search. `fc_master_render_delivered` takes the
handle's callback the same way `fc_master_solve_delivered` and `fc_master_measure_lra` already did (a new code,
`FC_PROGRESS_RENDER`; no struct grew) and cancels the same way. Checked in `tools/tests/MasterAbiProgressTests.cpp`.

## v0.36.0 — 2026-09-16

<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

### tools · wasm — `fc_probe` ships as an ES module too, for a module worker

**`tools/wasm/build.sh` now builds `fcprobe.web.mjs`**: the `fcprobe.web.js` line plus `-sEXPORT_ES6=1`, as
`fcmaster.web.mjs` has been built all along. A module worker cannot load the classic glue; it can `import` this
one, whose default export is `createFcProbe`. `fcprobe.web.js` stays, because `probe.html` loads it with a plain
`<script>`.

**One wasm under both glues.** Both web builds write `fcprobe.web.wasm`, so the node module is built first and
the script compares the web module's hash with `fcprobe.node.wasm` after EACH web build, not once after both,
since a single check at the end would only see the second. The `.mjs` glue goes through `check-no-threads.mjs`
and appears in the size table. Nothing in `fc_probe.cpp` changed, and the `fcprobe.web.js` line is the one it
was.

## v0.35.0 — 2026-09-16

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

<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

### convolution — an IR whose header claims a rate under 8000 Hz loads as is, like any other broken rate

**What `convolution::CabConvolver::loadIR` now promises about the IR's rate:** a KNOWN rate is a finite rate of
8000 Hz (`core::kMinSampleRate`) or more, and only a known rate off the host's is resampled. Everything else —
NaN, zero, negative, ±inf, and now any rate under 8000 Hz — is UNKNOWN, and an unknown rate loads the taps as
they are, with no rate factor (P67's rule: the samples are fine, only the metadata is broken, and refusing
would play silence). Above, nothing changed.

**Why:** a header that says 44.1 is kilohertz written as hertz, and it was trusted: a 4096-tap cabinet on a
48 kHz host was resampled x1088 into 4 458 231 taps — measured 3.76 s of `loadIR`, now 0.2 ms — and on the
verbatim (reverb) path scaled by the rate factor as well. Now it is the NaN load, bit for bit — the staged taps, `irNormalizationGain()`,
`irNormalizationGainDb()` and what the convolver plays, on both paths.

**What that moves, and what it does not:** the rate factor on the verbatim path now starts at 8000 / 3e6 =
2.67e-3 (−51.48 dB, an 8 kHz IR on a 3 MHz host) instead of about 4.7e-10. A 4096-tap cabinet at any known
rate costs at most x375 of itself now (a 3 MHz host); on a 48 kHz host it is x6 and 25 ms at 8001 Hz. What a
known rate can still ask for is unchanged and is bounded by the resampler's 2^24-tap output, not by the rate:
an 8001 Hz IR of 2.8 million taps (a 350-second file) took 17.4 s of `loadIR` on a 48 kHz host, measured.

<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

### tools · wasm — loudness and clipped runs of a stream, read between the pieces

**`fc_probe` gains a streaming surface** for two instruments that were already streaming:
`analysis::DeterministicLoudnessMeter` and `analysis::ClipDetector`, behind one class, `fcore::StreamProbe`
(`tools/fcore_stream.h`). A page opens a stream with `fc_stream_create (rate, channels)`, feeds planar PCM
with `fc_stream_process`, and between pieces reads `fc_stream_loudness` (momentary, short-term and integrated
LUFS, samples consumed, gating blocks dropped past the one-hour store), `fc_stream_clips_count` and
`fc_stream_clips (h, from, out, cap)` — the runs decided so far, polled by index. `fc_stream_finish` decides the
last runs; `fc_stream_destroy` frees the stream. Nothing in either instrument changed.

**Handles, not a singleton**, so two streams run at once — at most 16. A handle is a serial looked up in a table
and never reused, so a stale, forged or failed handle is refused rather than dereferenced. Buffers are checked
as the rest of `fc_probe` checks them, `cap` is in doubles and the return in runs as in
`fc_probe_clips_runs`, and a refused piece that carried samples poisons the stream: every reader answers 0.

**The deterministic meter, not the system one** that `fc_probe_run` uses: the system meter's block energies
differ between native and wasm at 8000, 88200 and 192000 Hz. `fcore_measure stream` prints the same bytes as
`tools/wasm/stream-parity.mjs`, and CI diffs them at 44100, 48000 and 88200 Hz on the release and checked
modules. `felitronics_stream_abi_tests` feeds pieces of 1, 4096 and random sizes through three handles at once
and requires the same integrated loudness, bit for bit, and the same runs as the same instruments given the
whole buffer in one call.

## v0.34.0 — 2026-09-16

<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

### oversampling · saturation · limiter · poweramp — a strict oversampler that is flat to 20 kHz, as an option

`PolyphaseOversampler` centres a fixed cutoff at 0.45 fs, and at 44.1 kHz that sits below the top of the
audio band: **one round trip costs −1.80 dB at 19 kHz and −15.55 dB at 20 kHz**, and more taps make
20 kHz worse (−13.71 at 32, −19.17 at 120). The cutoff is not an oversight, and its derivation is now in
the header: the class wraps a nonlinearity, so an image of content between 20 kHz and fs/2 that the
interpolator only partly rejects is multiplied with that content and lands in the audio band. The
design is therefore **strict** — the transition must finish below fs/2.

- **`oversampling::CascadeOversampler`** keeps the guard and moves the band edge: a Kaiser 2x stage whose
  length and cutoff follow from the sample rate by a measured rule, then halfband 2x stages. Images and
  aliases of anything below fs/2 at **−91 dB or lower**, one pass within **0.0043 dB** up to 20 kHz, over
  14 rates from 8 kHz to 768 kHz and factors 2–64 — and **−90.9 dB / 0.0049 dB over every one of the 110
  first-stage lengths the rule can produce**. Powers of two only, rates 8 kHz (the core's floor, P51) to 3 MHz; refuses (law 11b)
  what it cannot build.
- **The price is latency, not CPU**: **131 base samples at 44.1 kHz 4x** (63 for the Kaiser stage) for
  about the same multiply count (572 against 512). Above 44.1 kHz it gets cheaper, but the two halves cross
  at different rates: in multiplies almost at once (508 at 44.7 kHz, 348 at 48 kHz), in latency only from
  50.5 kHz (**76 at 48 kHz, 28 at 88.2 kHz**).
- **A halfband first stage — the cheaper-looking cascade — cannot do this**, at any length:
  H(f) + H(Fs/2 − f) = 1 makes its image rejection at fs − a equal its pass-band deviation at a. A 79 + 23
  pair that loses 0.41 dB at 20 kHz leaves that tone's image at −32.5 dB; one flat to 20 kHz still passes
  don't-care content's images at −58 / −35 / −22 dB (r = 0.46 / 0.47 / 0.48) — both pinned in
  `felitronics_oversampling_tests`. Through a single-ended tube stage at +12 dB that is −48 … −56 dBFS of
  intermodulation in the audio band (its 2nd-order line at 0.9–1.8 kHz) from −20 dBFS of content at
  0.48–0.49 fs, against under −150 for the strict designs (P31 stand, not a test).
- **`oversampling::Topology`** (`Kaiser` | `Cascade`) and `oversampling::Oversampler`, which is either.
  `Saturator::prepare`, `TruePeakLimiterConfig::topology` and `PowerAmpStage::prepare` take it; **the
  default is `Kaiser` everywhere, and under it nothing changes** — same refusals, same latency, same bits
  (the whole suite passes unchanged). Under `Cascade`: `tapsPerPhase` is still range-checked and
  otherwise unused; the rate becomes binding (a Saturator or limiter at 500 Hz is refused under the
  cascade, accepted under Kaiser); the Saturator's dry path is delayed by the cascade's round trip;
  `PowerAmpStage` rounds the factor down to a power of two and clamps the design rate, as it clamps
  everything, and under the cascade its 4x and 32x are no longer sample-aligned (76 and 80 at 48 kHz).
- **The switch costs the default 32 bytes per stage** (the cascade is heap-held, only when chosen).
- **One source-level change for a product that reads budgets field by field**: the oversampler half of
  `Saturator::Storage` and `TruePeakLimiter::Storage` is now `oversampling::Oversampler::Storage` (the Kaiser
  fields live under `.kaiser`); `bytes()` and `fitsWithin()` are unchanged, and nothing in the tree reads
  deeper. `PowerAmpStage.h` now includes `<bit>` (C++20, which every module already requires).
- **The limiter's ceiling under the cascade**, measured with the ceiling suite's own witnesses: at
  44.1 kHz 4fs/9 is now delivered flat, so the grid allowance at 8x rises from 0.108 to **0.133 dB**
  (16x: 0.027 → 0.033); the 1.15 dB modulation envelope holds (worst 1.344 dB at 2x, inside 2.399), and at
  44.1 kHz, 4x and 8x, the dense excess is lower than Kaiser's (0.49 / 0.48 against 0.81 / 0.75 dB).
- **The shipped wasm module carries it without being able to use it**: `fcmaster.web.wasm` grows by about
  12 KB (222 389 → 234 728 B against the main this branch sits on), because the chain's stages now
  contain the switch, while `fc_master_config` cannot select the cascade. The chain's output does not move:
  `fcore_master` render (with and without the clipper) and solve are bit-identical to main's at 44.1 and
  48 kHz. `fcprobe.web.wasm` is unchanged.
- **Unchanged on purpose**: `ReferenceTruePeakMeter` stays `PolyphaseOversampler` at 4x/32 — it is the
  certified unit — and no stage's default moved. `MasteringChain` and the wasm ABI do not offer the
  topology.

<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

### nam · neural · core — `reset()` restarts the stream, and the price of it is published

`felitronics::nam::NamStage::reset()` was EMPTY, and the comment above it recorded what that cost with a
number: a dense 2001-tap capture that had played a tone answered **digital silence with 0.224604502320**,
and through `reset()` the same. The house verb means a stream RESTART — `eq::EqBand` separated it from
`clearAudioState()` for exactly this reason — so a stale network window speaking into the first receptive
field of the next stream was the defect and not the design. P24 closed the half that is a falling edge (a
lane the host stops handing over is fed the silence it is receiving); this is the other half, the lane that
is PRESENT and gets the caller's own samples.

- **What it now promises is INDEPENDENCE:** nothing the caller fed before the restart can be heard after
  it. Every lane that carried audio is fed the digital silence it still owes, in the call, in full, until
  its state is provably the state of a lane that was silent all along. That promise is **exact** — two
  stages fed different audio before the restart answer the next programme with the same bits, on real
  captures and synthetic ones, at every rate and both widths.
  What it does **not** promise is silence out: a capture answers digital zero with whatever its own biases
  make of it, fresh and restarted alike — measured on NAM's shipped examples at 48 kHz, **0.001195220510**
  for a real Standard and **9.266554832458** for the A2-max feature set. Exact zero is a property of the
  bias-free fixtures, which is what lets them witness the promise.
  Nor, in general, bit-identity with a stage prepared a moment ago: NAM's answer depends on how the stream
  is **cut into calls**, and a restart's chunking is its own. Measured against a stage prepared a moment
  ago: exactly 0 for every fixture in the suite and for a real `slimmable_wavenet`, and **1.037e-06** for a
  real Standard at blocks 64…512, where the restart's last chunk is short — with independence still
  exactly 0 for that same capture, so it is the arithmetic and not the state.
- **The rate-matchers are re-primed too**, because a restart re-anchors the audio-time clocks the way
  `eq::EqBand::reset()` re-anchors its `StateGrid`. Leave the two `core::StreamResampler` legs where the
  previous stream left them and the next programme runs at its sub-sample phase: **1.039e-06 over 5091 of
  5120 samples at 44.1 kHz**, with everything else fixed.
- **NEW API, `core::StreamResampler::clearAudioState()`** — the leg's state alone (zero history, `len`,
  `pos`), keeping the rates, the capacity and the 513 × 64 coefficients. `reset (rates, capacity)` is that
  class's `prepare()`: it reassigns both vectors, `shrink_to_fit()`s on the identity path and re-derives
  every coefficient through a windowed sinc with a Bessel evaluation per tap — **1.77 ms** for one lane's
  two filtering legs (the price is the kernel, not the memory: at an unchanged capacity the allocations are
  reused), and deliberately not `noexcept`. A live stream cannot restart through that; `clearAudioState()`
  is 10 ns. It cannot
  resurrect a refused configuration either: `reset()` leaves `len = 0` when its second allocation throws,
  and a restart that took its length from the buffer's size alone would put that object back to work
  through a coefficient table that is not there, so the class now carries an explicit validity bit.
- **NEW API, `NamStage::clearedSamples()`** — the restart's own odometer, kept apart from
  `drainedSamples()`. The audio cannot witness "the full length for a lane that was playing, the remainder
  for one mid-drain, and NOTHING for a lane that never played": past the debt the output is zero either
  way. A mutation that spends one sample less than the debt survives every audio gate in the suite and is
  caught only here.
- **`NamStage::reset()` is now `noexcept`, and it is an AUDIO-THREAD call whose cost is not the block's:**
  a whole drain length of inference per dirty lane — the field, plus the partitioned-FFT ring any `Linear`
  is charged, plus each rate-matcher leg's tap window, so it is bigger than `prewarmSamples()` and that
  getter is not an estimate of it. On an M-series core, per lane, a real Standard WaveNet is
  **3.77 ms at a 64-sample block — 282 % of that callback** — 3.46 at 256, 3.43 at 512; a real LSTM 1.3 ms;
  a dense 2001-tap `Linear` 0.13 ms. It allocates, locks and throws exactly where `process()` does (nowhere
  for `Linear`/WaveNet; upstream's per-sample Eigen temporaries for LSTM/ConvNet), and it is IDEMPOTENT —
  the debt is re-armed only by audio actually being fed, so a second restart with nothing in between costs
  nothing and a mono host pays for one lane (for a RECURRENT capture it is deliberately not idempotent —
  see below). It also grows faster than linearly as the block shrinks, because NAM's per-call overhead is
  paid `debt / maxBlock` times: a real Standard is 3.61 ms per lane at block 256 and **19.47 ms at block
  1**. There is no cheaper exact mechanism to substitute: NAM's own `Reset` with the prewarm off zeroes
  the Conv1D rings in 0.014 ms and still misses the prepared state by **4089 samples, worst 0.324**, and on
  a `Linear` with the FFT engine it allocates 46 times.
  The allocation carve-out is `process()`'s — the same code path, inherited and not added — and it is wider
  than the header's architecture names suggested: `wavenet_a2_max.nam`, a WaveNet in NAM's own example set,
  allocates **4 times per sample** in `process()` (1024 for one stereo 256-block, measured) and therefore
  inside a restart too. The LSTM/ConvNet row of that carve-out went the other way: at this pin a stereo
  LSTM restart of 48 000 samples allocated **nothing** through either gate, the house operator-new counter
  or Eigen's own.
- **A restart that arrives while the backend is UNPREPARED writes nothing, and is not dropped.** `prepare()` writes the
  new `maxBlock` before it can refuse, so an unprepared backend can carry a block of a billion beside a
  256-sample scratch — a restart that touched it is a heap-buffer-overflow, which is why it touches
  nothing there, ledgers included. Re-arming the debt at the next prepare does not cover the lane that is
  PLAYING (its debt is overwritten on every chunk it is fed), so without the parked intent "play, a
  refused prepare, `reset()`, a prepare that succeeds, play" hands back the old stream: **242 samples** of
  a delay(514) capture, measured. The request was honoured at the end of the prepare that can honour it
  — and `p85-p86-restart-reachable.md` in this directory has since made that bit unnecessary: every
  successful prepare restarts, asked or not, so the parked intent is gone and the sequence is closed
  by the stronger rule.
- **A recurrent capture stays the NAMED exception, in the mechanism and not only in the comment.** An LSTM
  lane that has already spent its drain reads a debt of zero and is still not empty — the repository's
  slow-cell fixture leaves 0.419413 there — so a recurrent lane that ever played is charged the whole
  half-second heuristic again at every restart (which is what NAM's own `Reset` does) and is never marked
  clean. What that leaves is measured, not promised away: 300 samples differing from a fresh instance,
  worst 1.49e-07, on a real capture.
- **Scope, by byte comparison against v0.33.0:** 6 082 560 samples over six captures × eight rates × three
  block sizes, with loads, clears, mid-stream re-prepares, refused calls, zero-length calls and width
  changes — **byte-identical** where `reset()` is not called. Where it is, 8.394 % of those samples move, by
  up to 0.586079 full scale, and that movement is the fix. On v0.33.0 the same run with `reset()` called
  three times per cell is byte-identical to the run without it, which is the defect stated as a measurement.
- **What the restart does NOT reach, each named with its number** rather than promised away: a recurrent
  cell; NAM's partitioned-FFT clock (swept over nine block sizes x eight rates it peaks at **1.788139e-07**
  against a stage prepared a moment ago, and is EXACTLY ZERO against one clocked to the same point — the
  engine's own arithmetic, not state this stage kept); and a capture
  whose conditioner is a model of its own (`config.condition_dsp`), whose memory neither NAM nor
  `detail::receptiveFieldFromConfig` counted at the time. That last one is a hole in the LEDGER and not in this verb:
  the identical number comes back through the untouched drain — **0.905147969723** after a full drain
  against 0.905148267746 after a restart — and a test now pins both halves of it.
- **Not fixed here, registered:** the conditioner ledger above (now closed — see the entry on the ledger
  answering for the whole model); `prepare()`, which has the same stale
  window (it is where the 0.224604502320 was first measured); and `rigplayer::RigPlayer`, which has no
  restart verb at all, so a consumer reaching this stage through the player cannot yet call the fix. The
  last two are what a consumer actually hits, and both are closed by
  `p85-p86-restart-reachable.md` in this directory.

<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

### core · analysis · mastering · tools — one sample-rate floor, 8000 Hz, and a rate in kilohertz is refused

**What each entry now promises about the rate** — one number, `felitronics::core::kMinSampleRate = 8000.0`,
compared as `sampleRate >= kMinSampleRate` (8000 is a rate; NaN is not):

- **`mastering::TargetLoudnessSolver`** measures at 8000 Hz and up. `prepare()` refuses a lower rate and
  disarms; `solveBytes()` and `measureRangeBytes()` answer 0 exactly there. It took any finite rate > 0.
- **`mastering::MasteringChain`** is built at 8000 Hz up to its own 3 MHz, on every topology — the floor is the
  chain's, not a stage's. Before, what refused a low rate was whichever stage happened to be on (the limiter
  above 50 Hz, the EQ above 20.4 Hz), so a chain without both took 1e-305 Hz. `admits`, `prepareBytes`,
  `reprepareBytes` and `mastering::createBytes` say the same.
- **`core::DeliveryResampler::plan`** plans integer rates from 8000 Hz (its literal floor was 1000), so
  `mastering::DeliveryConverter` and `DeliveredMastering` refuse a SOURCE or a DELIVERY rate below 8000 — in
  `prepare()`, in every budget, and in `deliveredFrames()`, which answers -1 there. This is the line that covers a delivering handle's source rate: its chain
  runs at the delivery rate and never sees the source.
- **The mastering C ABI:** `fc_master_create` and `fc_master_need_create` answer `FC_ERR_REFUSED_BY_CORE` for
  `sampleRate` (or a non-zero `deliveryRate`) below 8000, with nothing allocated. No ABI version change: the
  surface did not grow, and the refusal comes from the core the facade already asks.
- **The probe ABI and `fcore_measure`:** `fcore::Probe`, `fcore::ShapeProbe` and the six analyzers behind
  `fc_probe_*` — `ClipDetector`, `ProgrammeReport`, `HumDetector`, `LowEnd`, `SourceForensics`,
  `BandBursts` — take 8000 Hz and up (to 768 kHz). The six analyzers and the probe had a floor of 1000 Hz, each
  in its own copy, and the shapes had none; each class keeps its `kMinSampleRate` name, and every one is
  now the core's constant. `fc_probe_<mode>_storage_bytes` answers 0 below the floor. The ceilings did not
  move here; the two that were missing — the shapes' and the search's — are P104's note. `fcore_measure`'s
  `correlation` and `needle` modes do not read the rate and still take any finite positive one; every other mode
  refuses below 8000 and names the range.
  A refused `fcore::Probe` or `ShapeProbe` now reads like a fresh one — its meters and its peak and stereo parts
  are replaced by new ones, and their memory is released — instead of serving the previous file (the rates that
  used to re-prepare them are refusals now), and
  `DeliveredMastering::sourceRate()` / `deliveryRate()` read 0 after a refused prepare, as the chain's and the
  search's rates already did.

**What starts to be refused, and why that is a finding and not a loss.** Every call below 8000 Hz that was
accepted before — the probe and the analyzers between 1000 and 7999 Hz; the chain, the search and
`fc_master_create` anywhere below 8000 their stages let through (60, 88.2, 96, 192 Hz on the default
topology; 44.1 with the limiter off); delivering handles from or to 1000…7999 Hz; the shapes at any positive
rate their decimation could hold. Measured on the tree before this change:

- the BS.1770 K-weighting shelf is designed at 1681.97 Hz, so below 3364 Hz it is past Nyquist — aliased
  everywhere there, and **unstable** wherever the bilinear tan() comes out negative (1682–3364 Hz, 841–1121 Hz,
  …): `fcore_measure lufs 3300 2 fixture.f32` (the CI fixture) printed **+3048.86 LUFS**; a 0 dBFS 400 Hz sine
  reads −3.72 LUFS at 48 kHz and +3043 LUFS at 3300 Hz;
- a search at 3363 Hz answered TargetUnreachable **at its −60 dB gain rail**, chasing +2448 LUFS;
- a search at **88.2 Hz** — 88.2 kHz spelled in kilohertz — answered **Solved at −14.0 LUFS**, exit 0, with a
  whole mastering chain running a thousand times too slow; 384 and 768 did the same; `fcore_master render 96`
  wrote its file with no status at all.

8000 Hz is the lowest standard audio rate, so no real programme sits below it, and it is clear of the shelf's
edge (its pole radius is 0.43 at 8 kHz, 0.99997 at 3364 Hz). A `static_assert` in `KWeightingFilter.h` keeps
the floor above twice the shelf. **Not changed:** `analysis::LoudnessMeter` and `KWeightingFilter` still take any
rate, and the meter still reads a rate ≤ 0 as 48 kHz; every caller that takes a rate from outside is floored
above them, and the meter's own contract is a separate decision. A page that builds a Web Audio context below
8 kHz on purpose is refused too.

**Tests that ran below the floor moved above it, with their fixtures kept:** `HumDetectorTests` from 3 kHz to
12 kHz (order 15 — the same 0.366 Hz bin; the noise is the 3 kHz draws, interpolated, and every measured
prominence moved by at most 2e-4 dB), `LowEndTests` from 6 kHz to 12 kHz (order +1, every sample count ×2, the
noise of its band-statistics rows interpolated, the kick's click scaled to its per-bin power; the fixtures that
are only compared with an oracle computed from the same samples draw fresh noise), `ClipDetectorTests`'
short-stream null and pending-queue rows and `ClipsExposureTests` from 1 kHz (W = 20) to 8 kHz (W = 160), their
quiet gaps kept at the same number of windows. P41 F1's solve at 1e-305 Hz is unreachable
now; its property — a meter sized in samples — is pinned on `LoudnessMeter::prepareForSamples`, which still
takes such a rate. **CI:** the clips gate refuses 7999, one ulp under 8000, 44.1 and 1000, and measures 8000;
the price table's rates are 999, 1000, 7999, 8000, 11025, 12000, 16000 … and its floors are 475 rows / 284
priced (425 / 280 before — the 1000 Hz rows are refusals now, and the 2000 Hz rows left the grid).

<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

### test_support · tools — one allocation counter, and it proves itself before the suite runs

"`process()` does not allocate" is the repo's oldest RT claim and it was being measured by sixty-one
private instruments, fifty of which could not see the allocations that matter most.

A suite proves that claim by replacing the global allocation functions and reading a delta. The set to
replace is **eight**, not two: C++17 routes any object whose `alignof` exceeds
`__STDCPP_DEFAULT_NEW_ALIGNMENT__` (16 on every desktop row, **8** on wasm32) through
`operator new(std::size_t, std::align_val_t)`, a different function. On `v0.33.0` the idiom had been copied
into 61 translation units; **11 installed the over-aligned form and 50 did not**, and both kinds said "no
heap allocation in the audio path" in the same words. `core::SeamAllocator<64>` — and through it every
convolver buffer — goes that way, and so does `eq::EqEngine`, whose `alignof` is 64 and whose own header
already said a two-form counter would not see its 331 KiB.

The hole was **latent**, and that is part of the finding rather than a softening of it: re-measuring with a
fixed instrument found nothing that had fallen through (see the last bullet). `modules/eq/tests/` is the
sharpest case — the suite the blindness would have hurt most makes zero over-aligned allocations of its
own, because it puts the engine on the stack.

- **`test_support/alloc_counter.h`** is now the only counter in the tree: all eight `new`s and all twelve
  `delete`s, the byte accounting (including MSVC's container padding and the plain-object exemption the
  `win` row found), a one-shot fault switch and a one-shot re-entry hook. The 61 copies are gone —
  **1390 lines deleted from the 61 suites against 348 added**, and one 363-line header in their place. Including the header IS installing it; there is no macro,
  because a macro can be forgotten and a forgotten one is silent, whereas including it twice in one
  executable is a duplicate-symbol link error on every linker in the matrix.
- **The instrument is proven in the run whose conclusions depend on it.** Before `main()`, the header asks
  the counter for one allocation through **each of the eight forms** and requires the number to move; any
  that goes unseen is named and the run aborts. Not a check the suite could fail to reach, and not in
  `report()` — two suites here total in a harness of their own and never call it. Eight probes rather than
  two because the distinction is not the obvious one: *deleting* a form is harmless (the standard's
  defaults forward `new[]`→`new`, nothrow→throwing, array-aligned→scalar-aligned, so the allocation is
  still counted — measured on libc++, libstdc++ and the UCRT alike), while a form that is **present and
  stops counting** is exactly the original defect in a new shape. One mutation per form: eight of eight
  red, each naming itself.
- **`tools/lint/check-alloc-counter.mjs`** keeps it that way: a private `operator new` anywhere under
  `modules/`, `tools/` or `test_support/` fails CI by name. Lexed rather than grepped — the words
  "operator new" are ordinary prose in this tree — with a 32-case self-test and two negative controls in
  the workflow, one for the lint's reach and one that deletes the two over-aligned lines from the shared
  counter and requires the binary to refuse to run. It anchors on the operator rather than on the return
  type and decides scope by brace depth, because a return-type matcher lets through every replacement whose
  signature is spelled differently — a newline after `void*`, a `[[nodiscard]]`, a trailing return type —
  and falsely flags a class's own allocator.
- **What the fixed instrument then found, stated as a number rather than a reassurance.** Every
  over-aligned allocation in the tree was traced: 1741 distinct call stacks across all 128 test binaries.
  They come from `prepare()`, `setIr()` and `Bank::build()` — `AlignedVector::assign` under a convolver's
  `prepare` accounts for 1305 of them — and **not one stack names `process`, `analyse` or `applyGain`**.
  A second stand gave each over-aligned allocation a weight of 10⁶, so one landing inside a measured
  region could not be absorbed by a tolerance; exactly two suites moved, and both were among the eleven
  whose counters already saw them. So the blindness was real and nothing had slipped through it.

**What the counter still cannot see is named in the header** rather than left to be discovered: it replaces
the C++ allocation functions, not `malloc`, and `pffft.c` calls one directly. Measured rather than assumed —
with `--wrap=malloc,calloc,realloc` on the gcc row, 200 `process()` calls through
`MatrixConvolver<PffftRealFft>` make **zero** C allocations (four in the whole run, all inside `prepare`).

No audio moved: the full verbose output of all 132 tests is identical to `v0.33.0` apart from wall-clock
timings, on arm64 macOS and on MSVC alike, and all 128 `N checks, M failures` lines match line for line.

<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

### convolution — 🔴 SOUND CHANGE: an un-normalized IR keeps its level across sample rates

`CabConvolver` prepared with `normalize=false` convolved **louder in proportion to the host's sample
rate**. A 48 kHz IR played **+6.02 dB** hotter in a 96 kHz session than in a 48 kHz one, **+5.28** at
88.2, **+12.04** at 192, and **−0.74** at 44.1 — the same file, the same knob, a different clock. It now
plays at the same level at every rate, which **moves shipped sound** by exactly those figures. **Nothing
changes for a load that is not resampled** — the IR already at the host's rate (within the 1 ppm match
tolerance), or one whose rate is unknown — where the multiplier is exactly 1 and the staged taps are
byte-identical to what v0.33.0 published. A 48 kHz host is NOT by itself such a case: a 96 kHz IR played
6.02 dB quiet there and now does not.

- **Why it happened.** `resampleIr` normalizes each output tap to unity DC, which preserves a WAVEFORM's
  amplitude — right for a signal, and what the function is for. A convolution's gain is not an amplitude:
  it is a sum over taps, so it is proportional to how many of them fit into a second. Resampling 48 → 96 kHz
  produces twice as many taps of the same amplitude, hence twice the gain. JUCE's convolver applies exactly
  this factor on exactly this branch — `juce_Convolution.cpp:790-792` in 8.0.14: `normalise == yes` measures
  the resampled taps, `else` calls `resampled.applyGain (originalSampleRate / processSpec.sampleRate)`. The
  interim JUCE backend this replaced therefore had it; the replacement carried over the normalized half of
  that `if` and not the other one.
- **And this repo's own spec said so.** `docs/migration/convolution-core-api.md:30` — the parity table whose
  file header calls it "the authoritative spec … without changing the sound" — has carried the row
  **`Normalise::no still applies a gain | makeEngine, applyGain(irSr/hostSr) (:792) | after any resample,
  multiply IR by irSr/hostSr`** since the migration was written. The requirement was not missed for want of
  being written down: nothing ever checked a row of that table, and the same file's pseudocode ten lines
  later contradicted it ("every KNOWN rate" where the table says "after any resample"). The pseudocode is
  corrected here; the table was right all along. And the golden-test tier that spec defines for a
  non-host-rate IR is `Normalise::yes` — the one tier that structurally cannot show this.
- **Why only this path.** With `normalize=true` the reference-unity gain is measured from the FINAL,
  resampled taps, so it absorbs the factor whole. The defect was visible only on the one path nobody
  measured that way.
- **The fix** is one factor, `irSr/hostSr`, applied where the loader already applies its one gain. It is
  named — `convolution::convolutionRateGain(inSr, outSr)`, beside the resampler — because a second copy of
  this arithmetic was already written by hand in a shipped product, which is the class of defect where two
  copies must agree forever. A caller that resamples an IR itself should use it instead of spelling it again.
  Rates it cannot make a factor from (a NaN, a zero, a negative, an infinity, or a finite pair whose quotient
  overflows) return 1 — the same rates `resampleIr` refuses get no compensation.
- **The factor is NOT clamped, and its sibling is.** The normalization gain is measured from the IR's own
  content, so ±30 dB stops a pathological IR blasting; the rate factor is arithmetic on two numbers the
  caller supplied, and clamping it would quietly deliver a different filter than the caller asked for. It
  spans 4.7e−10 to 4.3e9 (−186 to **+192.6 dB**), which is what `resampleIr`'s length and position gates
  leave reachable. A WAV's rate is a `uint32` and nothing in this family validates it, so a garbage-but-
  **finite** header is the one broken metadata a real file can carry — and it now plays LOUD where it used
  to play quiet: a file claiming 352800 Hz on a 48 kHz host is +17.3 dB, one claiming 5e6 Hz is +40.4 dB.
  Both are correct by this contract and both are garbage. A consumer that loads untrusted files should
  bound the rate before this, the way orbit-amp's loader already refuses anything outside 8 kHz…768 kHz.
- **What it does NOT fix.** The density term is all it removes. A resample still drops the kernel's
  pre-ringing that would fall before output sample 0, and what that costs depends on the onset, not on this
  factor — and it is a RIPPLE, not a loss: the kernel's nearest pre-ring lobes are negative, so cutting
  them ADDS level. Measured on a lone impulse at input index `lead`, 96 → 44.1 kHz: **−2.57 dB** at lead 0,
  **+0.65** at 1, **+0.95** at 2, −0.05 at 3, −0.61 at 4, and so on down to exactly nothing from `halfTaps`
  (32 input samples) of lead onward, which is the window's own backward reach. A smoother onset pays far
  less: a one-pole `exp(−n/τ)` at full scale with no lead loses 0.21 dB at 44.1 kHz for τ = 1 input sample
  and 0.004 dB for τ = 1 ms, which is what a real decay looks like. It is older than this change and
  untouched by it, but it is the one case where an **onset-trimmed** IR's level still moves with the rate —
  so a measurement that trims the front must not charge it here **in either direction**.
- **A diagnostic that had gone stale.** `irNormalizationGainDb()` floored its argument at 1e−6, which sat
  below everything the old code could apply (the normalization clamps at −30 dB) and above what this one
  can: an IR file claiming 0.024 Hz on a 48 kHz host applies 5.0e−7 and the reading said −120.0000 dB for a
  gain that is −126.0206. The floor is now below anything the loader can produce — **and a NaN gain now
  reports a NaN**, which is the one reading that did move on the normalized path: `std::max(a, b)` is
  `(a < b) ? b : a` and every comparison against a NaN is false, so `max(floor, NaN)` returned the FLOOR,
  a plausible −120.0000 dB for a gain that is not a number (a NaN tap in a float32 WAV is enough). The
  linear accessor was right throughout.
- **`irNormalizationGain()` now reports whichever gain was applied** — reference-unity on the normalized
  path, the rate factor on the verbatim one — and is still exactly `1.0f` when nothing was applied. The name
  is historical; the contract is "what was multiplied in", which is what a reference nulling against the
  engine needs.
- **One deliberate departure from the JUCE branch cited above:** JUCE applies the factor even when its
  resampler early-returned unchanged (equal rates), so it multiplies by a number that is only approximately
  1. This applies it only after a real resample, which is what keeps the byte-verbatim guarantee at matched
  rates. The two differ by at most `kRateMatchTolerance` relative — **8.7e−6 dB** — and the spec table's own
  wording ("after any resample") is the one this follows.
- **Pinned by** a property over eight rate pairs × five physical frequencies (44.1 / 48 / 88.2 / 96 / 192 kHz
  hosts, 44.1 / 48 / 96 kHz sources), read by a DTFT that shares no code with the loader, at a tolerance of
  0.01 dB — **13×** the 64-tap Kaiser kernel's own passband ripple (β = 8 → δ = 8.6e−5 → 0.00075 dB), **18×**
  the worst reading over the grid (5.65e−4 dB), and **74× below** the 0.736 dB the weakest cell moves without
  the fix. The kernel's column sums — the quantity the factor actually inverts — were measured across 64
  phases of every standard rate pair and stray at most 1.6e−4 dB from the ratio, which is the part of that
  argument that is a measurement rather than an engineering formula. Fifteen mutations of the change are each
  caught, among them "remove it", "invert it", "apply it only when upsampling", "apply it only to mono",
  "apply it on the normalized path too", "pre-scale the taps before the normalization is measured" and
  "delete the resample and keep the gain" — that last one is why an existing P67 assertion was rewritten:
  "every tap moved" stopped witnessing that the resampler ran the moment a gain could move every tap by
  itself, so the witness is now "no single scalar maps the input onto these taps".

### convolution — the bounds the IR resampler did not have, and the fifteenfold it did not need to cost

**`convolution::resampleIr` is 12-22x faster, and not one output bit moved.** `besselI0` ran once per tap:
37.7 ms for a one-second IR at 44.1 kHz, 164.6 ms at 192 kHz, on the message thread, every cabinet change.
Measured split of that loop: the Bessel series 78% of it, `std::sin` 13%. What replaces it is not a table
and not an approximation — it is an EXACT memo on the fractional phase, so the kernel arithmetic is
untouched and the answer is the same double it always was:

* a tap's weight depends on the output position only through `xx = t - k`, and the window's indices are
  `c + j` for the same offsets `j` about `c = floor(t)`;
* `frac = t - c` is EXACT for every `c >= 0` (at zero it is `t`; above it, Sterbenz), so `t` IS `c + frac`
  and the exact real behind every `xx` is `frac - j` — the SAME real for two outputs that share a `frac`,
  whatever their `c`. IEEE subtraction is correctly rounded, so both produce the same double from it. `xx`
  itself need not be representable; only the sameness of the number being rounded matters.
* `c < 0` is excluded because `frac` is not exact there: at 8 -> 48 kHz outputs 0 and 6 carry the same
  `frac` to the bit and five of their 64 weights still differ in the last place. A fixture in the suite
  holds an input where that difference reaches the float32 result.

It pays because audio rates are small rationals: a one-second 48 -> 44.1 kHz resample has 913 distinct
phases for 44100 outputs and reuses 98% of its windows, and 48 -> 96 kHz has two. A rate with no period at
all (a corrupt file rate, 48000 -> 44101) reuses nothing, and the memo gives up after 4096 misses without a
hit rather than make that case slower — measured 1.0x. Measured after: 2.7 ms at 44.1 kHz, 8.0 ms at
192 kHz. The suite runs every case through a frozen copy of the un-memoized kernel and compares every
output float as bits.

**Bounds that were missing.**

* **`kMaxResampleSamples` (16.7M samples, 64 MiB of float per channel).** The old length gate only kept the
  arithmetic addressable, so a result one sample under `INT_MAX` was 8.6 GB asked of the heap in a single
  call, on the message thread — and a file rate does not have to be absurd to ask: 1.2 Hz against a 48 kHz
  host is a ratio of 40000, and a one-second IR then wants 7.7 GB. An output past the bound is REFUSED, the
  rule P67 already ratified for a result that cannot be addressed. The figure equals
  `MatrixConvolverNupc::kMaxIrSamples` on purpose: the resampler can produce anything the convolver can hold.
* **`IrResampleConfig` has ranges.** `beta` was gated by `isfinite` where a RANGE was meant, so `beta = 1e300`
  passed and the 64-term series returned a silently wrong number (it only overflows to `+inf`, and the window
  to NaN, at about 1.36e4 — long past where its answers stopped being answers). `kMaxBeta = 53.0` is where
  that series stops meeting its OWN convergence test, measured at 53.038057. `kMaxHalfTaps = 4096` bounds one
  output's tap loop, which `halfTaps = 1e9` did not. A value that is not a REQUEST (a NaN, a negative beta, a
  radius below one) is repaired as before; a request past a ceiling is refused, because answering it with a
  smaller kernel would hand back a filter the caller did not ask for.
* The vector overload checks a length before narrowing it to `int`.

**`CabConvolver::prepare` refuses what it cannot honour (law 11(b)), where it used to guess.** A rate that
was not given was answered with the factory 48 kHz and a convolver was then sized from the answer; an
infinite rate made two float-to-int conversions undefined (on arm64 that produced `LLONG_MAX` and a
one-sample "crossfade", i.e. a hard switch). The rate must now be in `(0, kMaxSampleRate = 3e6]` — the house
figure — and `maxIrSeconds` must be a non-negative number; the IR budget saturates in double BEFORE the cast,
through the public `maxIrSamplesFor`, which is total for every argument. A refused prepare now also clears
the pending-retry geometry it can no longer publish: that survived a refusal, and `isBusy()` then answered
true for the life of the object while `flushPending()` could never publish, because it returns on
`! prepared_`. The clear happens before every refusal, so a refused WIDTH (`prepare(..., 4)`) drops a
pending load too, which it did not before — three behaviour changes in this call, all three named.

**An out-of-bounds write in the reference-gain path.** `normalizationGain` took the first second of the IR
as its analysis window but sized its transform at `N`, which stops doubling at 1<<21; above a megasample of
window the two parted company and the copy ran off the end of the buffer — measured under AddressSanitizer
as a 12 MB heap-buffer-overflow WRITE for a 3 000 000-sample IR at a 4 MHz host, reachable through the public
`loadIR`. The window is now the first second OR the transform, whichever is shorter.

### core · analysis · oversampling · tools — the three libm derivations, closed, and the lint that keeps them closed

v0.33.0 named three places where a number that crosses a platform boundary was still derived through the
system libm, and left them open. All three are closed here, none of them moved a bit of any shipped audio,
and the audit that found them is a gate rather than a reading: `tools/lint/check-det-math.mjs` plus
`tools/lint/det-math-manifest.txt` enumerate every governed transcendental call in the tree — 331 of them
across 61 files — classify each one, and go red when a call is added, removed, or swapped for another.

**The oversampler did not need the frozen tap table the plan assumed.** `PolyphaseOversampler::designFilter`
moves to `det::sin` outright: the taps are narrowed to float, and that narrowing discards 29 of the bits the
libms can disagree about. The system and deterministic designs differ at 18 of 128 taps in double and at
zero of 128 in float, and the final float arrays are byte-identical across Apple clang/arm64, gcc 14/glibc,
emcc/musl and MSVC/UCRT. The margin is not luck: perturbing every `sin()` by a deliberate k ulp leaves all
128 taps unmoved up to k = 2^24, against a real spread of one to three. The tests pin all 128 taps against
the table the OLD path produced on those four rows.

**The FFT needed its seeds converted and nothing else.** `core::offline::fftInplace` takes its stage seeds
from `det::cos`/`det::sin`; the butterfly is left alone, because it was measured not to contract under
`-ffp-contract=on` on any row, with both controls in place (an FMA canary that fuses in the same build, and
a 1-ulp twiddle perturbation that moves the checksum). Pinning O(N log N) products against a flag no shipped
road uses would have been a real cost for nothing. Note the scope of what was OBSERVED, because it is not
the scope of what now holds: the before-and-after hash equality was taken on a 2^16 transform, whose 16 seed
angles are the ones that measurement exercised, while the consumers run at order 17 by default and
`HumDetector`'s AUTO order reaches 21 at 768 kHz. Cross-row agreement at those larger sizes is not observed,
it is by construction — `det` is one implementation compiled into every row — and that is the whole reason
for converting the seeds rather than measuring them again. `core::fft::ScalarRadix2Real` keeps its system seeds by
decision, listed in the manifest as `retain-rt`: that one is the partitioned convolver's, i.e. shipped audio,
and a deterministic route there is an additive class beside it, not a replacement.

**`gainToDb` was worse than described, and the shared function does not move at all.** It is called once per
sample on four paths, worst of them inside `TruePeakLimiter`'s oversampled loop, in the module that is most
of a render, and `det::log10` is 7.4x a system call. So the CONSUMERS move instead: `core::gainToDbDet` sits
beside `core::gainToDb`, sharing the one floor constant, and the programme report, the clip detector, the
reference true-peak meter and the loudness solver call it. `gainToDb` itself is bit-identical and the same
speed before and after.

**Three rounds against the work, and each found something the reading could not.** The gate did not catch
the regression it was written after — reverting `LoudnessSolver::peakDb` to `core::gainToDb` left the lint
green, because carrier calls were only checked inside the deterministic zone; carrier calls now join the
manifest's multiset. The matcher had four evasions and the manifest four phantom entries. The `det::` pin
was six points, which is not a pin: dropping one `volatile` from `exp2Frac` changed `det::pow10` across
100 000 arguments while all six pinned values still matched, so the pin is now a dense per-function checksum
captured on three rows. A NaN clamp written as `std::max (gain, floor)` returns NaN instead of the floor
with the whole suite green. A solver gate pinned only through digital silence let its threshold be raised
tenfold while a peak of 2e-10 made the report read -200 dB and the certificate -193.98.

### tools · analysis — the price of a measurement, asked before it is paid

`fc_probe_report_storage_bytes`, `_bursts_`, `_hum_`, `_forensics_` and `_lowend_storage_bytes` publish,
through the wasm ABI, the number each offline analyzer already computes for itself before it allocates a
byte (law 11d). Until now nothing carried it out of the module: JavaScript called `_run`, which went
straight into an allocating `prepare()`, and there was no way to ask the price first. The prices are not
small — `HumDetector` at 768 kHz and 16 channels asks for 352 688 184 bytes, the five together for
377 423 824, against a linear memory that stops at 2 GiB and also has to hold the caller's decoded input.

Each query takes the GEOMETRY only — `(channels, sampleRate)` — because the whole use is to ask before the
input buffer exists; `_run` keeps its own buffer checks, and a price above zero is not a promise about a
particular call. A refused geometry quotes a canonical `+0.0`. The answer is a `double`, not the
`std::uint64_t` the budget is accounted in: an i64 does cross this boundary on the pinned toolchain, as a
BigInt, and a BigInt is the wrong shape for a number a page has to add to its own input size and put in a
report — `bigint + number` throws and `JSON.stringify` refuses it, so every use site would need a
conversion first.

What the module does today when the allocation fails was measured rather than assumed, and the assumption
was wrong in both directions. It is worse than a bad return value: `-fno-exceptions` turns the failed
`operator new` into `abort()`, so the C entry point never returns at all — the caller reading
`_fc_probe_hum_run(...) === 1` is handed a thrown `WebAssembly.RuntimeError` instead of a status. It is
also not the death of the page: that error is catchable, and after catching it the module went on
answering, with the aborted mode's getters reading 0 as the `have*` discipline promises. No ceiling on
geometry was added — what a page can afford is the page's number, not this file's.

Also in this change, found while building the gate for it:

- **The export whitelist in `tools/wasm/build.sh` had four ways past it**, none of which its guard could
  see: it recognised only a return type of `int|double|std::uint32_t`, it matched the declaration at
  column zero, and its guard was a hand-written floor (`-ge 39` against 86 actual entry points). So an
  indented declaration was invisible to every scanner at once while `EMSCRIPTEN_KEEPALIVE` exported it
  anyway; two declarations on one line dropped the second from the list with the counts still agreeing;
  a one-line body calling another `fc_*` function exported the wrong symbol; and a comma declarator
  (`FC_EXPORT int fc_a (void), fc_b (void);`) is one token, one line and one name while the second is
  simply absent. All four were reproduced against the gate before they were fixed. Both ABIs' lists are now read as the identifier
  before the first `(`, one declaration AND one declarator per line are enforced, and the count is
  checked for EQUALITY against the declarations; a second gate refuses a return type the boundary does not carry rather than
  dropping it. `fc_master`'s copy had the type defect too, plus a name pattern without digits that would
  have TRUNCATED `fc_render_v2` into the list rather than omitting it.
- **Two rationales in the tree were false** and are corrected with the measurement that settled them:
  `-sWASM_BIGINT` is ON by default in emscripten 6.0.9, and `EMSCRIPTEN_KEEPALIVE` does keep a function
  that is missing from `-sEXPORTED_FUNCTIONS`.
- **`planarSpan()`'s channel bound was pinned by nothing in the repository.** Deleting it left all 123
  tests green, because every other entry point refuses the width a second time in `prepare()` — every one
  but `fc_probe_needle`, which has no `prepare()`. `felitronics_abi_tests` now pins it there.
- The five modes' default parameters are one `constexpr` constant each, read by the run, by the price and
  by the result getters that used to construct their own.

The gate: `felitronics_analysis_abi_tests` goes from 389 checks to 575 — the allocation a first run
actually asks for against the published demand, measured at the WIDEST geometry the ABI has and equal to
the byte; the price against the core's own `storageFor()` at every width and rate, each row asserted
positive so it cannot pass as `0 == 0`; a canonical `+0.0` on 75 refused geometries; a grid of five modes
by twenty-six rates by nine widths by three programme lengths, where the price must be positive exactly
where `_run` is accepted and every refused row must leave its getters silent; and the asymmetries that
are deliberate, written down so nobody "fixes" them. `tools/wasm/storage-probe.mjs` makes the assertions
only the wasm tier can make — starting with the one no native test can, that all five names reached the
artifact — and prints the demand table that `felitronics_analysis_abi_tests --storage-table` prints
natively: 425 rows, byte-identical across the two tiers, on the release and the checked module alike. (P51, in
the same release, moved the rate floor under this gate: 596 checks, 90 refused geometries and 475 rows since.)

Twenty-six mutants were run against it and twenty-three died. **The crew's testing round found the hole
the first twenty missed**: with `_run` ignoring what `prepare()` returned, the whole suite stayed green
while an empty programme at a refused geometry was ACCEPTED and its getters served the previous
programme's numbers — `_run` skips `process()` when there are no frames, so an unprepared analyzer still
reached `finish()`. Both grids had a single non-zero programme length. They now walk an empty one, a
short one and a full one, in both tiers, and the mutant of that line dies in four of the five modes —
in lowend it is equivalent, because `planarSpan` refuses an empty programme before `prepare()` is reached. Two more
survivors from the same round — an allocation that only happens above two channels, and one that only
happens away from 48 kHz — are what moved the allocation oracle to the widest geometry. The three that
survive are equivalent and named where they live: removing `setParams (kReportParams)` from report's run,
where the constant holds the instance's own defaults; hard-coding the width in the lowend query, whose
demand does not depend on it; and masking the low three bits off every price, since every demand in the
accepted domain is a multiple of eight. Six mutants of the build gate were killed too — the indented
declaration, the two on one line, the comma declarator, the empty extraction, the unsupported return
type, and the name taken from a body call.

<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

### nam · rigplayer · convolution — a prepared stage holds no audio, and a player can ask for a restart

P47 made `felitronics::nam::NamStage::reset()` an exact stream restart, and no product could call it.
Two doors stood in the way, and each was half the same defect. This opens both.

- **`NamStage::prepare()` now performs that restart, ALWAYS.** It carried the identical stale window —
  `::nam::DSP::Reset` calls `SetMaxBufferSize` and then a prewarm that is ZERO samples for a `Linear`,
  so `Buffer`'s per-channel input window survived the call. Measured over eight host rates x two block
  sizes x three capture shapes x {re-prepare at the same rate, re-prepare at a different one}, both
  lanes PRESENT: **72 of 96 cells answered digital silence with something, worst 0.567861497402**, and
  0 of 96 do now. A prepared stage and a just-constructed one name one state, so the two verbs of the
  class say one thing.
  There is deliberately **no predicate on what changed**: a re-prepare at the SAME rate and block is the
  common case — a host's buffer-size slider moves more often than its rate one, and a driver stops the
  stream for either — and it is where the leak was loudest (0.468718945980 at 48 kHz, against
  0.469410002232 across a rate change). A "fix it only when something moved" predicate leaves 8 of those
  96 cells leaking, which the mutation stand shows as red.
  The mechanism is P24's ledger and P47's drain, unchanged: `configureRates` charges every lane that has
  EVER been fed a full debt at the new rates, and the tail of `prepare()` spends it. So a first prepare
  after a load costs nothing, and a model change costs nothing either: it prepares a never-fed backend in
  `prepareModel()`, and again in `install()` when the host's numbers moved between the two halves. The `restartOwed_` bit P47 needed
  is **deleted**: once every successful prepare restarts, a parked request has no second question to
  answer, and the sequence it was written for is closed by the stronger rule.
- **THE PRICE, published rather than hidden.** For an architecture whose own `Reset` already prewarms
  (every WaveNet), the drain is a SECOND pass over the field and roughly doubles the call: on an
  M-series core, a stereo real Standard WaveNet `prepare()` on a dirty stage goes **6.5 ms → 13.1 ms**
  at 48 kHz and 7.6 → 14.2 at a 64-sample block; a real LSTM 2.5 → 5.8. That is the message thread, with
  no callback to miss. A first `prepare()` after a load does not move (7.16 → 7.25 ms), which is what
  `everFed_` buys and what the suite gates. Skipping the drain where NAM's own prewarm provably covers
  the ledger is a real optimisation and is registered as one rather than taken here — it is a predicate.
  One consequence is stated because it is a number that moved: against a stage prepared a moment ago, a
  real Standard WaveNet used to be **exactly 0** and is now **1.1e-06**, because the extra drain re-chunks
  a network NAM had already flushed and NAM's answer depends on how the stream is cut into calls. That is
  the same residue `reset()` publishes. INDEPENDENCE — the promise — stays exactly 0 on every capture.
- **`rigplayer::RigPlayer::reset()` — the verb a product can actually call.** The class had none, and
  orbit-amp reaches a `NamStage` only through it (`releaseResources()` is empty there and no host reset
  is overridden). It restarts both model slots (including one the blend law has put to SLEEP), the three
  convolvers **bypassed or not** — a bypassed one is skipped, so its history freezes and is replayed when
  the curve comes back — the dry path's alignment ring, the per-slot whole-sample alignment tails, the
  band filters, the audio-time grid and the scratch; and it SNAPS the gain ramps, as `eq::EqBand::reset()`
  does and for the reason that file records with a number. `prepare()` now calls the same two private
  bodies, so the two verbs cannot drift.
  It deliberately does **not** touch the blend law's state — a restart is not a device change. Re-arming
  the law's warm-up ledger would not deliver its invariant anyway: with `fed = 0` on both slots the law
  ramps its gain down over four blocks, so a by-hypothesis wrong-sounding network is audible regardless,
  and what it buys is a hole — 18 blocks, 192 ms, at every restart. Re-arming only the sounding slot is
  worse: the law rails the goal to the neighbour and plays a full spurious crossfade to the other capture
  and back. The price of the verb is FOUR networks, not one: two stages, each up to two lanes.
- **`convolution::MatrixConvolverNupc::clearAudioState()` (new, additive).** `reset()` there zeroed the
  history AND cancelled a swap in flight (`xfadePos_ = 0; state_ = 0`, keeping `cur_`), so a filter
  published a block ago and still crossfading in was dropped — and `CabConvolver::pendingRetry_` is
  already false after a successful publish, so nothing ever re-staged it: the knob move was lost until the
  next knob move. `clearAudioState()` is the history alone, touching only buffers `process()` writes, so
  a composite can restart on the audio thread without losing a filter or racing the message thread.
  `CabConvolver` forwards it. (The verb itself is fixed in the same release — see the P88 entry: `reset()`
  now ADOPTS the publication instead of dropping it, in all three convolvers. That removes the reason
  `RigPlayer::reset()` reaches for `clearAudioState()`, which leaves a fade running and so does not give a
  restart's independence mid-fade; switching the player is registered as its own task.)

Green: macOS/clang `ctest` 131/131 with `-DFELITRONICS_WITH_NAM=ON`, deb/gcc-14.2 131/131.

<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

### nam — the ledger answers for the WHOLE model, conditioner and heads included

`felitronics::nam::detail` is one registry with three questions — how far a capture's memory reaches,
whether anything in it is recurrent, and whether anything in it is charged NAM's partitioned-FFT ring —
and **two readers spend its answers**: law 11a's drain (a lane the host stops handing over is fed silence
for that long) and P47's stream restart. It walked `config.submodels` and stopped there, so a capture
whose **CONDITIONER is a whole model of its own** (`config.condition_dsp`, which NAM builds with `get_dsp`
like any other model) hid that model's memory from both. The two were short by exactly the same amount,
which is what said it was one defect in one ledger rather than two: **0.905147969723** left after a full
drain against **0.905148267746** after a restart, on the same capture.

- **The composition is a SUM, and that is the number rather than the wording.** The conditioner's output
  is the network's conditioning input, so the two memories are in SERIES — which is also how NAM's own
  arithmetic composes them. Measured from outside the registry, on a loaded two-layer stack with a
  2500-sample Linear conditioner: the impulse's last non-zero sample is **5000**, where the same stack
  without a conditioner reaches 2501 and a worst-of would have answered 2502.
  What the registry promises is an **upper bound**, and the slack is named: the condition enters a layer
  AFTER that layer's own convolution, so the true reach is `Mₒ + H + max(0, M_c − L₀)` and the answer is
  `Mₒ + H + M_c + 1` — over by `L₀ + 1`, never short.
- **All three questions walk the branch, not one of them.** A `Linear` conditioner is charged the
  2048-sample ring (its instance owns the same engine, and its ring feeds the layer arrays through the
  mixin); an `LSTM` conditioner makes the whole capture recurrent, because the cell's state enters every
  layer through a memoryless 1×1 and no finite silence empties it. A mutation stand confirmed each of the
  four one-function-only variants goes red, as do "walk it but discard the field" and "worst-of instead
  of the sum".
- **And two more branches of the same config were not being read**, both of which the sum needs in order
  to be an upper bound at all. A layer array ends in a causal head rechannel whose kernel is the layer's
  own `head.kernel_size`, worth `kernel − 1` samples — `example_models/A2.nam` spells it 16 — and a
  post-stack `config.head` is worth `Σ(kᵢ − 1)`. A `ConvNet` keeps its whole stack in a TOP-LEVEL
  `dilations` array and answered **zero**. Each is normally hidden by NAM's own answer, which `NamStage`
  raises this number with; none of them is hidden behind a **slimmable** WaveNet, which answers zero for
  everything. Measured on ones that load: a slimmable WaveNet with a 16-tap head reads 2 against an
  impulse reaching 16, and one with a ConvNet conditioner reads 2 against an impulse reaching 15.
- **The area of the change, bit-for-bit.** Rendered through both readers at 44.1 / 48 / 96 kHz, NAM's nine
  shipped example captures are **byte-identical** except the two that carry a conditioner, whose drain
  grows by **one sample** (`wavenet_a2_max` 31 → 32, `wavenet_condition_dsp` 45 → 46).
- **The answer is monotone, and it took a review round to make that true.** The first cut of this change
  added the ConvNet reader as a link in a first-non-zero CHAIN, ahead of the declared-field reader — and
  a chain can take a number away. A `Linear` capture carrying a stray `dilations` array loads (its parser
  reads neither key) and the chain answered **2** for a 4999-sample impulse response that answers 4999
  without the stray key: a lane draining 2050 where it needs 4999, i.e. a regression introduced by the
  fix. The three sources are now the stack, then the WORST of the other two. In the same round the
  container branch stopped RETURNING: NAM dispatches on the `architecture` string and never on shape, so
  a `"WaveNet"` carrying a stray `submodels` array loads with its whole stack behind that return
  (measured: 1 answered for a model reaching 4200), and the return also made `isRecurrent` charge a
  conditioner the field was ignoring — the very divergence this change exists to close. Both are gated.
  With those two, every term is added and none is replaced, so no capture can drain shorter than before.
  The one named exception is a number that is not representable at all: a field spelled `1e300` is now
  REFUSED rather than cast, because the cast is undefined behaviour.
- **A dead key must never silence a live reader, and a diverse-testing round found two more of those.**
  NAM's ConvNet parser does not read `layers`, so a ConvNet carrying `"layers": []` had its whole stack
  suppressed — behind a slimmable capture's conditioner, where NAM answers zero, that drained 102 for a
  model reaching 131 and handed back 29 samples. And a boolean is a number to `get<int>()` but not to
  `is_number()`: NAM loads `"dilations":[true,true]` as `[1,1]` and answers 3 for it, where this file
  answered nothing. Both gated.
- **A value that does not fit an `int` is now REFUSED rather than clamped**, because the model NAM built
  does not contain it: `"dilations":[4294967396]` builds a network whose dilation is 100, and clamping to
  INT_MAX made `reset()` spend 2 147 483 646 samples — a measured **23.7 s of synchronous audio-thread
  work** — for a capture that remembers a hundred. Refused reads as absent, and NAM's own answer covers
  what it did build.
- **The gates were checked by mutation, twice, in an isolated copy of the tree:** 21 variants of this
  file, 20 red. The one survivor is an `is_object()` guard proved equivalent (nlohmann's `contains()` is
  false for every non-object, without throwing). A sweep of 3200 generated loadable configs found no
  answer short of the model's measured impulse reach outside law 11a's named recurrent exception.
- **Not fixed here, registered:** the hybrid slimmable-wrapper shape, where NAM dispatches on a top-level
  `config.layers[i].slimmable` marker and the real config then hides under `config.model`. It loads, NAM
  answers 0 and the registry answers 0 however long its stack is (measured 314 samples of leak on one).
  Reading it means restating NAM's dispatch heuristic, which is its own decision and its own number.
  *(Closed by P92 without that restatement — see `p92-unplaced-shape-ceiling.md`.)*
- **And one that is not this module's at all, registered with its reproducer:** a `condition_dsp` that is
  a `SlimmableContainer` with a WaveNet submodel, or a slimmable WaveNet, **crashes the host on load** —
  a null write inside NAM's own `WaveNet::_set_condition_array` during `DSP::prewarm()`, because neither
  class overrides `SetMaxBufferSize` and the conditioner path is the one place that is the only call they
  get. Reproduces on pure `nam::get_dsp` + `DSP::Reset` with none of this code in the path, and
  `prepareModel`'s catch-all cannot catch a SIGSEGV.

<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

### convolution · lineareq — a restart keeps the operator the caller published, in all three swap-safe convolvers

A swap-safe convolver publishes an operator from the message thread and the audio thread adopts it at the
END of its crossfade, so between those two moments the live slot still names the PREVIOUS operator. All
three `reset()` bodies wiped the publication flag and kept that slot: the operator the caller had already
replaced came back, and nothing was left to re-stage the new one — `setIr()` had already returned true, and
a consumer's retry flag is clear after a successful publish. On those bodies each class's new independence
gate reads the OLD operator in every cell {Pending, Crossfading} x {mono, stereo} x both slot parities
(a worst sample **7.494e-01** `ConvolutionEngine`, **1.017e+00** `MatrixConvolver`, **1.927e+00**
`MatrixConvolverNupc`, from a restart made after the fade settled). The promise those bodies carried was
already the right one — "flush the tail, not revert the EQ" — and reverting to a superseded operator is
exactly reverting the EQ. Law 11e in `docs/DSP-ARCHITECTURE.md` states the rule.

- **Through the product classes** — `lineareq::LinearPhaseEq` and `NaturalPhaseEq` at 48 kHz, a +12 dB
  bell at 1 kHz published over a settled flat curve, the host restarting before a block picked it up — the
  gain at 1 kHz read **+0.00 dB** at 1, 2 and 4 channels (above two channels each channel has its own
  convolver, staged then published) and now reads +12 dB (gated to 0.2 dB). When the curve was the FIRST
  one ever published, the EQ answered **exactly zero** (the fixture's -600 dB floor) until the next band
  move, because the live slot still held the zero taps `prepare()` leaves; it now plays the bell. A restart
  with the move in flight now answers the next programme with the same bits as a restart after it settled
  — a level alone cannot tell `reset()` from an EQ that forwards to `clearAudioState()`, which reaches the
  same curve 20 ms later. `LinearPhaseEq` is also gated through an M/S move (the MSDiag operator, different
  L and R), by level and by independence. The shipped products did not reach this — each re-publishes
  after a re-prepare, and orbit-amp drains silence rather than calling `reset()`, naming this defect as one
  of its two reasons (`src/core/CabinetIr.h:220-222`; the other, the concurrency contract, stays) — so it
  was latent in the products and live in the API.
- **`reset()` now ENDS a swap in flight in favour of the new operator** in
  `convolution::{ConvolutionEngine, MatrixConvolver, MatrixConvolverNupc}`: `cur_` moves to the published
  slot and the state returns to Idle, so `isBusy()` is false at once and the consumer may publish again.
  An operator merely STAGED and not yet published (`stageOperator()` without `publishStaged()`) is
  untouched — it has not been accepted, and the later publish still finds it. The rule is scoped to the
  restart; `prepare()` still discards everything, by contract.
- **WHAT MAKES ADOPTION RIGHT IS LAW 11a INDEPENDENCE, NOT CLICK-FREEDOM.** A half-finished fade is a
  dependency on what came before: two convolvers holding the same published operator, one mid-fade and one
  settled, would answer the next programme differently for up to the length of the fade. Through
  `CabConvolver`, a restart one block into a 50 ms fade and one after it settled now differ in **no** sample.
  **The price is at the seam and it is published** (the engine's suite prints every number here): flushing
  the history is itself a cut — the previous stream's tail stops mid-decay — so the first sample after
  `reset()` steps whichever slot is live, **1.8203e-01** on DC 0.5 into a settled 700-tap operator.
  Adoption puts the new head tap where the old one was at that sample, so it moves the step by at most the
  head-tap difference times the input, in EITHER direction: **7.3203e-01** for a pair whose head tap flips
  +0.70 -> -0.40, **1.3933e-01** for a +1 dB broadband move and **2.2018e-01** for a -1 dB one, and exactly
  the flush for a change that leaves the head tap alone. That difference IS the change the caller asked for.
- **`clearAudioState()` is now on all three** (it was on `MatrixConvolverNupc` alone). It is the history
  and nothing else: it decides nothing about the operator and leaves a crossfade running from where it was
  — which is exactly why it does NOT give independence mid-fade: through `CabConvolver`, a clear one block
  into a 50 ms fade and a clear after it settled differ in exactly the fade's remainder, **2143 of 5120**
  samples a channel (worst 4.186e-01), and the engine's own suite shows 383 of 1600 between two clears
  inside one fade. `rigplayer::RigPlayer::reset()` calls this verb, a choice made when `reset()` still
  dropped the filter; switching it is registered as its own task, not taken here. `reset()` is built from
  `clearAudioState()` in every class, so the two verbs cannot drift.
- **`ConvolutionEngine`'s per-channel `buildIr()` no longer zeroes the staged slot's cached tail.** That
  was a MESSAGE-thread write into a history buffer the audio thread also writes — the same `0.0f`, so
  nothing could be lost, but a data race all the same (put back, ThreadSanitizer reports it), and the one
  thing that made `clearAudioState()`'s promise untrue for this class; its two siblings never write a tail
  off the audio thread. It is redundant: `primeTail()` overwrites that buffer whole at fade start for every
  channel being processed, and a channel NOT being processed was zeroed when it dropped out and is cold,
  for which zero is the right value. Byte-identical output across the change, measured against the tree
  before it over 24 cells of {no reset, reset at Idle} x width {1, 2, 2->1->2, 2->0->2} x block size.
- **The memory order is the one the new write needs.** `reset()` writes `cur_`, which the loader reads
  after acquiring `state_`, so the Idle store is a RELEASE store; and it is made only when a swap was
  actually in flight, because an unconditional `store(0)` would wipe a publication that landed between the
  load and the store. A new real-thread gate, `ConvolutionRestartRaceTests` (outside the wasm-audio tier,
  like core's RtStreams suite), loses **none** of its publications, where the bodies before lose **all** of
  them and the fixed bodies with the store made unconditional lose **12–95 %** per route in every run
  measured (macOS/arm64, Debian/x86-64 with gcc, Windows/x64 with MSVC). With that and the tail write gone,
  `reset()` no longer races a single-producer loader at all: ThreadSanitizer, a loader publishing in a loop
  beside process / reset / clearAudioState on all three classes, reports **0** races against **54–62** on
  the bodies before; the release store put back to relaxed draws **63–71**, and the tail write put back
  **2** (Apple clang needs `-fno-builtin` for the last — it does not instrument a `std::fill` otherwise,
  and a first run here was blind to exactly that; gcc 14 sees it as is). The documented "must not run
  concurrently" contract is kept until a sanitizer row carries that.
- Gates, as a property and not a point: `{Pending, Crossfading} x {mono, stereo} x both slot parities`
  and, for the matrix siblings, every topology with four distinct banks and different L/R inputs, with
  **law 11a independence** as the oracle (two convolvers fed different audio of different lengths, one
  restarted with the publication in flight and one after it settled, answer the next programme — longer
  than the whole FDL span — bit for bit), plus the published operator certified from OUTSIDE the class
  against a reference convolution; `clearAudioState()` mid-fade (still busy, ends on the sample it always
  would have, history gone, lands on the new operator, and not independent of where it was called); for
  the engine, the seam bound and no allocation in either verb; the product classes above; `CabConvolver`'s
  two verbs side by side, the clear by exactly the fade it left running; a fuzz of random publish /
  process / reset / clearAudioState interleavings on all three classes against a settled shadow of the
  same class (with a random width and a reference convolution for the engine), red on the first seed of
  every class before the fix; and the real-thread race gate above. A mutation stand over **17** mutants and
  9 test binaries (the engine, matrix, Nupc, gate-state, `CabConvolver`, fuzz, race, `LinearPhaseEq` and
  `NaturalPhaseEq` suites) is red on every one — including "adopt but leave the state in flight" (the
  operator rolls back one fade later), "`cur_ = 0` instead of `1 - cur_`" (caught only at the second slot
  parity), "leave the FDL", "leave the frame", "forget one slot's cached tail", "the adopted slot's
  topology slips", "the clear restarts the fade", "an EQ forwards `reset()` to `clearAudioState()`", "the
  Idle store made unconditional" (caught by the race gate alone) and each class's pre-P88 body. Four more
  are green on the suites by design: a dropped `ranNc_ = 0` and a dropped `xfadePos_ = 0` are truly
  equivalent; the release store put back to relaxed and the tail write put back are races that only
  ThreadSanitizer sees (above).

<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

### rigplayer · nam — a prepare restates what it counts, and a drained lane is not billed twice

`RigPlayer::prepare()` rebuilt everything the host rate DESIGNS — filters, rings, both stages, the
cold threshold kept in seconds — and nothing the host rate COUNTS. This closes that half, and one
double charge in `NamStage`.

- **A woken slot now warms for the field the NEW rate owes.** The blend law's warm-up debt was latched
  from `warmFor()` when a model landed and read for ever after, including by the wake of a sleeping
  slot. Measured through the player on a 6x6 host-rate grid, two block sizes and three capture shapes:
  the warm-up after a rate change did not depend on the new rate at all. **48 → 96 kHz warmed a slot
  for 2048 host samples where 4096 are owed — half the field; 48 → 192 kHz a quarter; 192 → 48 kHz held
  a slot silent four times too long.** Of 180 cells, 85 under-warmed, 85 over-warmed and 10 read right
  — eight of those 44.1 ↔ 48 kHz, so the pair a fixture reaches for first could not see it. **180 of 180
  now match a player prepared at the new rate from the start.**
- **The ledger is restated by the rule each count's algebra allows**, in a new law verb,
  `nam::blendRestated()` — not `blendLanded()`, which would also clear a load in flight, a sleep and a
  refusal. The debt is RECOMPUTED: it is not homogeneous in the rate (the rate-matcher's latency is 0 at
  48 kHz, 61 at 44.1, 96 at 96, and the block term does not scale), so rescaling it is a second copy of
  `warmFor()` with a different answer. The progress is mapped by its PREDICATE: an audible slot stays
  audible, a warming slot starts over. Scaling it by the rate ratio is the trap — 96 → 48 kHz turns a
  just-audible slot inaudible for any block over 96 samples and plays a spurious full crossfade. The rest
  count IS rescaled: it is pure elapsed time, so the ratio is exact, and is 1 where nothing moved.
- **A restart inside a warm-up re-arms it, at any rate — `reset()` included.** Both verbs flush every
  network, and the law went on crediting a warming slot the field it heard before the flush. Restarting
  1, 2 and 3 blocks into a 2001-tap field left 14, 13 and 12 blocks of warm-up against the 15 a restart
  at the landing costs. A warming slot is at weight zero, so re-arming it is silent; an AUDIBLE slot is
  still not re-armed, for the reasons P85 gave.
  *This supersedes the P85 note above that `RigPlayer::reset()` does not touch the blend law's state:
  it does not, except for this.*
- **The rest a slot has served survives a restart.** It was left in the old rate's samples against a
  threshold recomputed in the new: 48 → 96 kHz slept 0.73 s late, 96 → 48 kHz after one block. A first
  draft zeroed it instead — a behaviour change at the UNCHANGED rate, postponing every parked dial's sleep
  by two seconds at each restart — and was caught before it shipped.
- **The per-slot alignment delays are restated too** — they are host samples, written only by a knob, a
  table or a landing — and **a landing delivered before a rate change and taken after it** warms and
  aligns at the new rate. A `prepare()` between `load()` and the next block leaves the previous pack's
  model ids alone: they index a model list that now belongs to the new pack.
- **Slot trims switched OFF stay off across a restart.** The snap ignored `setInputTrims(false)`, so
  a restart ducked the slot to the pack's trim and ramped back over ~43 ms (−4.29 dB on the first block
  of a −6 dB entry). The shared body fixes `prepare()` and `reset()` at once.
- **`NamStage`: a lane whose falling-edge drain ran to the end is not billed again by the next
  `prepare()`.** Only a restart cleared its "may be holding audio" flag, so the next prepare charged a
  whole drain for a clean lane: 4093 samples for `wavenet_a1_standard`, 6347 for `A2`, and the
  re-prepare measured 14.7 ms where 12.7 is owed. A PARTLY drained lane is still charged in full, and a
  recurrent one always is. *This supersedes the P85 note above that `configureRates` charges every lane
  that has EVER been fed.*
- **`rigplayer` has a permanent RT allocation gate** (`RigPlayerRtAllocTests`), on the shared counter
  that sees all eight forms of `new`: `process()` on five shapes, across a turn, through a slot falling
  asleep and draining, a host dropping to mono and to a gap, a call four times the block; `reset()`; the
  per-block read-outs. Its control plants an allocation and must read one more than the same window
  without it.

**Registered, not taken:** skipping `prepare()`'s second warm-up pass for architectures whose `Reset`
already clears their state (P98). With the drain removed outright, every capture NAM ships keeps
independence at exactly 0 — and this tree's own `Buffer`-based fixtures leak on 72 of 96 cells, worst
0.567861. NAM's example set has no such capture in any tree, so a check against real captures alone
would have approved a blanket skip. Found on the way and registered: the alignment ceiling
`kBlendMaxDelay` is 128 samples, so a 64-sample lag at 48 kHz is clamped to half at 192 kHz (P99).

Independence (P47/P85) is unchanged: nothing the caller fed is audible after a restart — exactly 0
across 384 rows of NAM's own captures through the player. One number did move, and it is not
independence: since a lane drained to the end is no longer drained again by `prepare()`, two stages fed
the same audio but a different channel-width history now differ after a prepare by the drain's own
chunking residue — up to 2.4e-06 on `slimmable_wavenet` at a 17-sample block, 0 on `A2` — where they
used to be identical.

Found by the adversarial round and registered, not this branch's: a lane that comes back after its drain
has run out resumes at a frozen sub-sample phase, so at a non-integer rate ratio it differs from a lane
fed silence throughout by 0.114 on `wavenet_a1_standard` at 44.1 kHz (P101); and `RigPlayer::prepare()`
accepts positive rates far below any audio rate (1e-3 Hz) that its stages cannot honour (P102).

<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

### nam — what the ledger cannot place costs ONE allowance: never zero, never the face value

`felitronics::nam::detail`'s receptive-field registry promises an **upper bound** on a capture's memory,
and two readers spend its answer: law 11a's drain and P47's stream restart feed a lane that much digital
silence. Where it could not read a shape it answered **zero**, which is the tail of the previous sound
coming out of silence. It now answers what it read, plus one named allowance.

- **The measured case is NAM's own shipped capture.** `example_models/slimmable_wavenet.nam` with its
  config moved under `config.model`, and a top-level `layers` carrying nothing but the slimmable marker,
  **loads** (NAM's WaveNet parser delegates on the marker and the delegate reads `config.model`), makes
  the same sound, and its impulse reaches sample 2046 either way. The registry answered **2047** for the
  flat file and **0** for the wrapped one, and NAM's own answer for the architecture is `return 0`. All
  three questions were blind together: an LSTM conditioner inside the wrapper read as not recurrent, and a
  dense `Linear` one lost its 2048-sample ring.
- **Three events mean "cannot place", all read off the config and none off NAM's dispatch:** an object
  carrying a model's vocabulary under a key the registry does not read (at a config, a layer entry or a
  submodel entry, singly or in an array — under any key name); a value that is there and cannot be read
  (a slimmable dilation of `4294967396`, which NAM builds as 100, answered 0); and a reading the registry
  sets aside (a declared field beside a stack). Each adds `kUnreadShapeCeiling` = **48 000 samples**,
  **once per tree**, to what the registry did read, and charges the ring. It fires on **none** of the
  1229 distinct captures on the author's machine; the flat shipped slimmable still reads 2047.
- **Why once, why added, why not the face value** — each is a measured failure of the alternative. An
  allowance per node turned a 1.6 MB file of 30-byte dead siblings into an **INT_MAX** drain; a max at the
  root let a known 100 001-sample stack swallow the allowance of an unreadable stage in series with it;
  and trusting the face value of a dead `"receptive_field": 2147483647` beside a real Standard's stack —
  a file NAM loads unchanged — made `reset()` run for **about half an hour per lane** (2³¹ samples at 0.79–0.89 µs each).
- **The price, measured through `NamStage::reset()`** on real captures rewrapped, per lane at 48 kHz:
  **39.7 ms** on `wavenet_a1_standard` at a 256 block (44.5 ms at 64), **24.1 ms** on A2's submodel,
  **4.5 ms** on the slimmable itself. On shapes that ship: identity.
- **What moved, on purpose, and only in synthetic rows:** a stale `receptive_field` beside a stack
  (9 → 48 009), a stray top-level `dilations` beside `layers` (2 → 48 002), refused values (1e300,
  `4294967396`, a string) from "absent" to the allowance. **Closed on the way:** a 5000-tap `Linear`
  carrying a readable stray `layers` array answered 2 (reach 4999; now 48 002), and a `ConvNet` carrying a
  non-empty dead `layers` array answered 0 (NAM: 256; now 256).
- **What it does not close, stated rather than hidden — two doors, and shutting either opens the other.**
  A LIVE memory the registry cannot place and that is LONGER than the allowance drains short by the
  difference (a wrapped model with an inner field of 60 000; a 60 001-tap `Linear` whose declared field is
  set aside beside a dead `layers` array) — the longest real field is 6347. And a DEAD number the registry
  PLACES is still trusted at face value, exactly as before: a lower reading with no stack beside it (a dead
  `dilations:[2e9]` beside a `Linear`'s declared field), the wrapped form's own decoy stack, a `layers`
  array on an architecture that never reads it. Which door stays open is registered as a policy question.
  A differential fuzz of 2 000 loadable configs against NAM found no short answer outside these two and
  the recurrent exception. Separately, NAM's own recursive copy of the config takes the host down on a
  file about 2 000 levels deep before the registry runs. The only new recursion — looking for an
  architecture through unplaced nodes — is bounded at 32 unplaced hops; the nesting the registry already
  walked is read exactly as before, and a guard on it was measured to drain a 65-deep conditioner chain
  short.

<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

### mastering · tools — the search and the shapes get the ceilings their neighbours already had

**What each entry now promises about the rate, from above:**

- **`mastering::TargetLoudnessSolver`** measures up to 3 MHz — `kMaxSampleRate`, which is
  `MasteringChain::kMaxSampleRate`: the search measures what a chain renders, at the chain's rate. `prepare()`
  refuses a higher rate and disarms, and `solveBytes()` / `measureRangeBytes()` answer 0 there. It took any
  finite rate over the floor: at 1e300 Hz `prepare()` said yes while `solveBytes()` said 0, and one ulp over
  3 MHz it prepared for a chain that cannot exist.
- **`fcore::ShapeProbe` and `fc_probe_shapes_run`** take 8000 Hz to 768 kHz — `fcore::Probe`'s range. They
  kept no ceiling and drew a 1 MHz file that every other run entry of the probe ABI refuses. `fcore_measure
  waveform` / `stereo` refuse above 768 kHz too, and name the range.

Nothing changes at or under either ceiling: a verdict oracle over 3817 rates and 39 entries finds differences
only in the search's and the shapes' entries, only above their ceilings, and only from an acceptance to a
refusal.

## v0.33.0 — 2026-09-15

<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

### analysis · tools · wasm — the clip detector has a way out

`analysis::ClipDetector` has found sample clipping since v0.28.0 and nothing outside the core could ask it.
It can now, through the two roads every other measurement of this repo takes, and the detector itself is
unchanged.

- **`fcore_measure clips <sampleRate> <channels> <raw.f32le> [--max-runs N] [--chunk N]`** prints how many
  runs were found, whether the list is whole, the sample peak of every channel, and each stored run — start,
  length, level, channel, polarity, evidence. Levels, peaks and the sample rate go out as **raw IEEE-754 bit
  patterns**, as `blocks` does, so the comparison can tell a double from its own float32 rounding: a run of
  eleven samples at 40/64 and one at 39/64 has level 479/768 = `3fe3f55555555555`, whose float32 is
  `3fe3f55560000000`, and `%.6f` prints 0.623698 for both. Nothing is in dB — `samplePeakDb()` routes through
  `log10`, and libm is not bit-identical across toolchains.
- **`fc_probe_clips_run` + eleven getters** (twelve exports) in the wasm ABI, shaped like
  `fc_probe_shapes_run`: one run, then getters, with the run buffer owned by the caller and its capacity
  mandatory. Each run is six doubles (`start, length, level, channel, sign, evidence`); positions are below
  2^30 by the ABI's own span bound and so exact in a binary64. The export list is generated from the source,
  so the twelve came along by existing. Both copiers take their capacity in **elements**, like the other
  seven copiers of that ABI, and `fc_probe_clips_runs` answers in **runs** — a caller who passes the element
  count it would pass to any neighbour gets too few runs, which it can see, rather than a six-fold heap
  overwrite, which neither `outSpan` nor `SAFE_HEAP` could see.
- **Roads:** `fcore::ClipProbe` (`tools/fcore_clips.h`), shared by both, exactly as `fcore::Probe` and
  `fcore::ShapeProbe` are. It shares the **lifecycle**, not merely the report reader, because that is where
  every trap of this mode lives. It is the one class here that asks for something the detector does not need —
  the stream's length — because without it an adapter has to believe the `n` it is handed and reads past the
  caller's buffer; and because the two roads knew the length anyway, so two duplicated "did the file deliver
  what it was sized for" checks became one.
- **Four ways an exposure could have certified a file it never measured**, all closed, each with a negative
  control in the suite. (1) A report read before `finish()` is short by up to `decisionDelaySamples()` —
  20 ms — of runs. (2) A refused `process()` does not stop `finish()`, so the detector reports a complete,
  empty, clean file; a refusal that carried samples now poisons the measurement and clears validity on the
  spot. (3) `isFinished()` is **not** a validity flag: `prepare()` disarms and then returns early on a refused
  argument without clearing `finished_`, the run list or the counters — measured on this tree, a good run
  followed by `prepare(0.0, …)` leaves `isFinished()` true, `runCount()` 1 and `samplePeak(0)` 0.625, the
  previous file's answer behind a flag that says the measurement is done. (4) A short read: the CLI sizes the
  file first and the adapter refuses unless the length the file was sized for, the frames the reader handed
  over and the samples the detector consumed are all one number.
- **Two capacities, kept apart.** `complete` is the detector's list overflowing `maxRuns` — a property of the
  file, reported as data with the count still counting. A short copy into the caller's buffer is the caller's
  business, and a small buffer never makes a file incomplete. Note that capacity 0 does **not** make a report
  incomplete by itself: completeness is `count <= capacity`, so a file with no runs is complete at capacity 0.
  `maxRuns` is bounded at 2^20 on both roads, not at `ClipDetector::kMaxRunsLimit`: 2^24 runs is 665 MiB at
  16 channels and 768 kHz, and the module is built `-fno-exceptions`, where a failed allocation aborts the page.
- **Law 8a is tested on the shipped adapter, not on the detector.** 140 slicing comparisons through
  `fcore::ClipProbe` — whole, 1, 2, 3, 7, W±1, kChunk±1, past kChunk, mixed and ragged — over ten programmes
  (T = 0, 1, W−1, W, W+1, clamped material at 8/48 kHz, one to three channels, non-finite holes, a clipped
  tail), each also compared with a bare detector read outside the wrapper, and each with the allocation
  counter on. Beside them an **intermediate trace**: the reference is built one sample at a time, recording
  the exact coordinate at which every run becomes visible, and each slicing is checked at every one of its own
  call boundaries — a comparison of final reports is not enough (`docs/LAW8-KWEIGHTING.md:78`). And beside
  that an **outside oracle**: every one of 1368 reported runs has its level recomputed from the plane data as
  an exact rational mean and compared bit for bit, because fourteen slicings agreeing with each other is
  consistency and not truth.
- **Parity:** 64 byte-identical native-vs-wasm comparisons over 69 658 runs locally (four rates from 22.05 to
  96 kHz, one/two/six/sixteen channels, eight capacities, thirteen `--chunk` values), and a CI step that pins
  38 of them — 26 successful diffs and 12 REFUSALS, where both roads must exit non-zero and print nothing,
  because a byte diff of two successful commands says nothing about the command lines both sides are supposed
  to reject. Release and checked (`SAFE_HEAP` + `ASSERTIONS=2`) artifacts both. `--chunk` is honoured natively
  and ignored by the module, which makes each of those rows a **cross-tier law-8a test**: without it both sides
  cut the stream at multiples of `kChunk` and the diff would say nothing about re-slicing. The comparison's
  sharpness has its own control: one flipped bit in a level, injected into a scratch copy of the formatter, is
  caught.
- **A clamped fixture generator** (`tools/wasm/make-clip-fixture.mjs`), transcendental-free like its sibling —
  the existing one is clean audio as far as this detector is concerned, so a parity run on it only proved the
  two sides agree nothing is there. Its passages are fractions of the requested length, so even a
  one-second fixture carries all seven, holes and clipped tail included, and its NaN is written as an explicit
  bit pattern because `writeFloatLE(NaN)` stores an encoding ECMAScript lets the engine choose.
- **Mutation stand, 13 of 13 red**, on isolated copies of the sources: a band closed at a call boundary, the
  pending queue drained per call, a level rounded through float, the chunk loop without its plane offset,
  `runsComplete()` always true, `finish()` ignoring the poison, the poison not clearing validity, the shim fed
  half the buffer, the report published without `finish()`, the peak copier filling the caller's capacity, the
  run copier answering in doubles, the three clocks no longer compared, and the run copier reading its capacity
  as runs. The float-rounded level is the one every slicing comparison and the whole parity harness pass —
  only the pinned bit pattern sees it.
- **Green on four toolchains**: Apple clang arm64 (114/114), gcc 14.2 x86-64 on Debian (114/114), MSVC 14.44
  (111/111) and the `wasm-audio` tier under node (113/113), where the allocation counter is enforced rather
  than informational.

<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

### analysis — `ProgrammeReport`: one offline pass over a delivered programme, and what it refuses to say

A new offline analyzer, `felitronics::analysis::ProgrammeReport`, shaped like `analysis::ClipDetector`
(`setParams` / `prepare` / `process` / `finish` / `reset`, `Storage` + `storageFor()` published before a byte
is allocated). One pass produces: DC offset per channel; leading and trailing silence, both at a named
threshold and as threshold-free exact digital silence; the tail — the energy of the last `tailWindowMs` of
programme relative to the programme's own mean square, plus the last sample and the last sample before the
trailing silence; the infra-low energy fraction through `eq::Crossover2`; the stereo relations (L/R balance,
bit- and value-identity counts, the integral correlation via `analysis::StereoSums`, side/mid energy); and
the dynamics — integrated loudness, reference true peak (`analysis::ReferenceTruePeakMeter`, 4× / 32 taps,
drained), PLR, EBU Tech 3342 LRA, the short-term loudness percentiles and the crest factor. Also
`fcore_measure report`, which prints every scalar as a raw IEEE-754 bit pattern with its validity and
reason, through the report's own field visitor.

**Every mean is taken over the PROGRAMME SPAN** — the first to the last frame whose loudest finite present
channel exceeds `silenceThresholdDb` — and the tail window ends at the last such frame rather than at the
last sample of the file. The property that buys: padding a master does not move its measurements. The DC,
sample peak, RMS, crest factor, infra-low fraction, programme mean square, tail ratio, stereo relations and
reference true peak of the same music padded with five seconds of silence at each end come back
bit-identical, pinned as an allow-list so a new padding-sensitive field cannot widen the claim quietly.
Measured over the whole file instead, every one of them moves for a reason that has nothing to do with the
music: a 60 s programme padded to 70 s reads 0.67 dB lower, and the tail window lands entirely inside the
padding and reads 0, so a truncated fade reports as a perfect one.

The span is thresholded rather than keyed on exact zero, which was the first design and is a cliff exactly
where masters live: a 24-bit dithered file — the dominant delivery format — has no exact zero anywhere, so
its five "silent" seconds are ±2e-7 of dither, an exact-zero span swallows them, the RMS is diluted by
5.4 dB and a programme cut mid-fade reads a tail ratio of 9.1e-12 instead of 1.0. What padding does still
move is named in the header: `lastSample`, which is by definition the last sample of the file, and the
loudness family, because BS.1770 anchors its gating blocks and 3 s windows at the start of the stream —
anchoring the sub-meters at the first active frame would make them invariant and put this report's loudness
at odds with `fcore_measure lufs` and every other EBU tool on the same file, which is the worse trade.

**"Cannot say" is a result.** A mono file does not get a correlation of +1.0; a programme shorter than the
tail window does not get a ratio of exactly 1.0; a programme with fewer than two gated short-term
observations does not get an LRA of 0.0, which is also the honest answer for a constant tone. Each is
`valid = false` with a reason and a canonical `+0.0`, never a NaN. Two shipped primitives are wrapped
because their own defaults are verdicts rather than measurements: `StereoSums::correlation()` answers +1.0
below its 1e-12 denominator ("silence is neutral" — right for a meter, wrong for a report), and
`LoudnessMeter::integratedLufs()` answers −120.0 when nothing passed the gates. LRA is computed from this
class's own short-term series rather than from `loudnessRangeLu()`, whose 0.0 cannot be told from outside
the meter to mean "no range" or "fewer than two observations survived the gates".

**The infra-low fraction is a filter, not a band, and the header says so in the arithmetic.**
`eq::Crossover2` is LR4 — two cascaded Butterworth sections — so |LP|² = 1/(1+(f/fc)⁴)², which is −6 dB at
the crossover, **not** the 4th-order Butterworth 1/(1+(f/fc)⁸) the request was written against. With
fc = 30 Hz a pure 30 Hz tone reads 25 % rather than 50 %, a 40 Hz tone 5.78 % rather than 9 %, and flat
noise 0.104 % at 48 kHz. A reader told the wrong formula would see a few per cent "below 30 Hz" on a
bass-heavy master and conclude there was infrasound, where it is the fundamental at 40 Hz through the
filter's skirt.

**Law 8a, gated on the intermediates and not only on the result.** Every field, and every event of a
test-only frame trace, is bit-identical under any re-split of the stream into `process()` calls — checked
over 19 slicings, 8 `maxBlock` values and anchored random slicings, field by field through `std::bit_cast`.
The trace exists because the finished report is not a sufficient gate, and that is measured rather than
argued: of a 23-mutant stand, 22 die and **not one of them produces a single report-comparison failure**.
A mutant that closes the 10 ms sub-hop at the end of `process()` leaves the report bit-identical at all 19
slicings and all 8 `maxBlock` values; so does one that moves `eq::Crossover2::flushDenormals()` off
`core::StateGrid`. Beside the invariance sit closed-form oracles computed outside this repository (a sine's
crest factor is 20·log10√2 = 3.0103 dB; the LR4 magnitudes above; EBU Tech 3341's −23 dBFS tone) and an
independent whole-file reference program that takes its percentiles by sorting rather than from a histogram
and writes out the 3 s / 1 s cadence as literals rather than reading the header's constants — an oracle
that imports the subject's parameters moves with a mutation of them instead of opposing it. Two real
defects in this work were caught by that reference alone and by no invariance check.

Non-finite input is a hole: a canonical zero into every filter, counted per channel, and then excluded where
exclusion is exact (the sums stay valid) and fatal where state carries it forward (the infra-low fraction
and the whole loudness family refuse, with the reason). Nor is a maximum exclusion-exact — the sample
dropped for being non-finite may have been the largest — so with nothing finite the sample peak is refused
rather than published as a comfortable zero. Finite input is not enough either: `eq::Svf` updates its
integrators as `(float)(2·v − ic)`, and that intermediate can leave the float range even where the state is
not growing, so a stream of finite `FLT_MAX` drives the 30 Hz crossover non-finite at sample 849 exactly,
grid flush running (a 1.7e38 DC input does not overflow, so "the state doubles every sample" is the wrong
mechanism to quote). The crossover's overflow and the K-weighting's are counted separately, so a refusal
never names a filter the field does not go through. Nothing non-finite reaches an accumulator, a published
field or the trace. Capacity overflow never refuses a call or stops a counter; it invalidates exactly the
fields it damages.

475 checks, green under Apple clang 21 / libc++ / arm64 and gcc 14.2 / libstdc++ / x86-64, including the
strict header-hygiene gate (`-Wconversion -Wfloat-equal -Werror`) on both. `fcore_measure report` on a
stereo fixture is bit-identical across those two rows — every count and all 31 scalars.

(When this was written the header added that cross-platform bit-identity was NOT promised, because `log10`
and `pow` are not bit-portable. P79, later in this same release, narrowed that caveat rather than removing
it: those derivations run through `core::det` now, and the report is byte-identical across Apple
clang/arm64, gcc/glibc x86-64 AND wasm32/musl — **on the fixtures this release ships, for a build that
does not contract**, which both shipped roads are. Two limits are named rather than glossed: rebuilt with
the library's default `-ffp-contract=on`, 4 lines move again, because the arithmetic AROUND the
deterministic calls is still contractible; and three derivations still reach the system libm — the FFT's
twiddle seeds, `core::gainToDb`, and the reference true-peak meter's tap design — so an input outside the
measured fixtures could still diverge. Those three are the next task's scope. What `core::det` removed is
the share no build flag could reach.)

### Added
- **`analysis::SourceForensics` — what a delivered file actually was, as far as the samples can prove it.**
  Two families of evidence and no attribution: a long-term (Welch) **spectral wall** — where the upper band
  limit is, in Hz *and* as a fraction of this rate's Nyquist, how deep the drop is, how wide the transition
  is, how far the spectrum comes back above it, and how far up it is empty (the upsampling tell) — and the
  **sample grid**: the coarsest dyadic grid `2^-k` every finite sample lies on, the shortest normalised PCM
  word that holds them all, the distribution of `k` behind that maximum, and the count of distinct sample
  values. Composed over the `analysis::SpectrumFrames` producer; `setParams` / `prepare` / `process` /
  `finish` / `reset`, `Storage` + `storageFor()` published before the allocation (law 11d), `maxBlock` sizes
  nothing, and every published field is bit-identical under arbitrary re-slicing of the stream into
  `process()` calls (law 8a).
- **`fcore_measure forensics`** prints the whole report, every float as a raw IEEE-754 bit pattern.

### Notes
- **A band limit is a suffix property, and its FLOOR is read past the transition.** The strict suffix
  maximum is anchored at the transition's end, not at the winning boundary: the winner sits on the cell
  that contains the band limit — that is what makes it the winner — so anchored there it reports a clean
  wall with 80 dB of "recovery" above it and a strict drop near zero (measured: 805 of 2350 constructed
  cutoffs, worst 80.4 dB, and one NEGATIVE drop).
- **Is there a floor to reach at all?** The floor reference is the median of the span above the edge, so on
  a stopband that is still descending it sits halfway down the descent and the transition "ends" half a
  span past the edge: an identical 16 kHz brickwall measured 49.8 Hz of transition over a flat floor and
  1594 Hz over one decaying at 15 dB/kHz, with `sharp` flipping false at 10 dB/kHz and true again at 25.
  The span's two halves are now compared; disagreement means there is no single floor, the width is a lower
  bound, and `transitionClipped` says so — a flag that was otherwise provably unreachable, so the promise
  attached to it had never been kept.
- **The grid's exact reading is a maximum, so a robust one is published beside it.** On float-rendered
  material — a 16-bit programme with a float fade-out, which is every real render — the fade's samples take
  `gridExponent` to 85 or 149 and the exact word length goes dark. `robustPcmBits` is the shortest word
  holding all but `gridOutlierFraction` (5 % by default, set by the length of a real fade) of the non-zero
  samples, computed from the histogram already accumulated; a tolerance of 0 makes it the exact reading.
- **A geometry that cannot finish is refused, not discovered at finish().** The edge search sorts one
  plateau span per candidate, so `fftOrder 22` with 0.01 Hz cells and a 12 kHz plateau span is 2.1 million
  candidates sorting a million doubles each. `storageFor()` and `prepare()` refuse it; the default geometry
  is 19240.
- **A band limit is a SUFFIX property.** The obvious construction — the steepest local descent, then the
  level just above it — reports a deep notch (a band-stop, a comb, a room null) as a 100 dB wall with a
  one-cell transition while full-power spectrum resumes 400 Hz higher, and it can select that notch and
  thereby discard a genuine wall above it. So the floor of the drop is the loudest cell *anywhere* above the
  candidate, and the search maximises that conservative drop. A monotone roll-off then fails on the drop
  itself (~10 dB) rather than on a transition-width technicality, and both negative cases are handled by one
  mechanism.
- **An edge that is not there is not given a frequency.** An argmax always returns something — the lowest
  candidate on a falling spectrum, the highest on a flat one — so `minDropDb` gates whether an edge exists
  at all: below it the report is `valid = false, ShallowerThanMinDrop`, with every number still published as
  evidence.
- **A narrow line above the edge is forgiven to a published rank and no further.** A 19 kHz whine or a
  leaked pilot tone is the loudest thing above a real wall; the search therefore runs on a 3-cell
  median-filtered copy of the cells and takes the `(t+1)`-th largest of the suffix, `t` capped by
  `exemptCells` and by a tenth of the suffix, so the exemption self-disables near Nyquist. What was forgiven
  is published (`sufMaxPower`, `strictDropDb`, `exemptedCells`, `recoveryDb`), and the emptiness test runs
  on per-bin maxima, which see every line.
- **The position near Nyquist is not attributed.** A converter's anti-alias filter, a 320 kbit/s codec and a
  genuinely band-limited master all live at 0.9-0.95 of Nyquist; `nearNyquist` marks that band and the
  instrument stops there. `nearNyquist == false` is not a claim that anything compressed the file.
- **Two edges, and which one is the report.** The primary is whichever candidate has the stronger
  conservative drop, so a 16 kHz codec wall inside a 20.5 kHz export filter reports the inner edge as the
  primary and the outer as the second — but only when the shelf between them is deep. With a shallow shelf
  the outer edge wins instead and the inner one is not reported at all: a second search BELOW the primary
  is structurally useless, since everything above such a candidate includes the primary's own plateau.
- **The PCM range is part of the claim.** A b-bit word carries `i/2^(b-1)` for `i` in
  `[-2^(b-1), 2^(b-1) - 1]`, so `-1.0` is a PCM sample and `+1.0` is not, at any depth. A sample outside
  `[-1, +1)` therefore withholds `minExactPcmBits` (reason `OutsidePcmRange`) while `gridExponent` stays
  valid, and the signed minimum and maximum the claim rests on are published.
- **What the grid proves, and in which direction.** `minExactPcmBits = k + 1` is the *shortest* normalised
  PCM word that holds the observed samples exactly; it does not bound the source's word length from above (a
  24-bit file carrying a padded 16-bit master is indistinguishable from a 16-bit one), and the *absence* of
  zero low bits proves nothing at all — 16-bit dithered up into 24 reports `<= 24`, never 16. A grid finer
  than `2^-23` is refused as a PCM word rather than reported as 32: float32's top binade is spaced `2^-24`,
  so a 32-bit stream reads `k = 24` near full scale and `k = 31` lower down, while a 32-bit *container*
  carrying 24-bit content reads 24 and one carrying 16-bit content reads 16.

### analysis — `HumDetector`: mains hum, and the notes that look like it

`analysis::HumDetector` measures mains hum: a narrow, stationary line at 50 or 60 Hz with a comb of exact
multiples, looked for only inside the QUIET stretches of a programme, because music masks it. It reports the
line's interpolated position in Hz, its level over the local background in dB, which harmonics were found,
and which stretches were used — in sample coordinates. It runs on `analysis::SpectrumFrames`, is bit-identical
under arbitrary re-slicing of the stream into `process()` calls (law 8a), and publishes its whole budget
through `storageFor()` before allocating (law 11d). `fcore_measure hum` prints the report, every float as a
raw bit pattern.

**No answer reads as "clean" when the instrument could not look.** Seven named incompletenesses, not one of
them a zero that could be mistaken for absence: no quiet stretch at all; exactly one (because "stands still
between stretches" cannot be tested inside one); fewer than two stretches long enough to show the line twice;
a window too short to separate 49.0 Hz from 50.0 Hz; every frame holed; a mains-compatible line that was
measured and did NOT stand still; and a comb of harmonics whose base lies outside this scope. The evidence —
positions, levels, prominences, counts, the strongest peak of each search window whether accepted or not —
is published in every case.

What separates hum from a bass note (G1 is 1.0 Hz from the mains; A#1 and B1 are 1.7 Hz from 60) is the
position measured INSIDE the bin by a three-bin parabola, a 0.5 Hz tolerance that a drifting grid fits and a
note cannot, and the requirement that the same line stand still across at least two quiet stretches and
across every frame of them. The comb is reported, not required: a bass guitar's partials are near-exact
multiples too, so gating on a harmonic count would cost a false clean and buy little.

Four things in it are not obvious and are the reason it works:

* **Every local maximum of the search window is examined, not the strongest one.** A bass note 1 Hz away and
  14 dB louder owns the window's argmax; reading only that would lose the hum underneath it, with the
  resolution paid for and unused. A lone Hann-windowed tone's sampled skirt is monotone, so enumerating
  costs nothing in false positives.
* **Resolution is a duration, stated in bins.** Separating two Hann main lobes 1.0 Hz apart needs 2.7 bins
  between them — 2.0 bins is where the dip disappears — so the bin must be 0.37 Hz or finer: N = 2^17 at
  48 kHz, a 2.73 s window. `fftOrder = 0` picks the shortest window that achieves it at the file's rate.
* **A line must be prominent AND loud.** A ratio alone certifies arithmetic: a programme that is one pure
  tone has a spectrum of 1e-23 elsewhere, and a maximum of that residue stands 18–37 dB over the residue
  beside it. `minLevelDbfs` (−100 dBFS) is the other half of the test. A frame of exact digital silence is
  excluded for the same reason — it cannot have seen anything, so it is not evidence of absence.
* **A quiet stretch must show the line twice.** One frame is one periodogram, whose tail crosses a 10 dB gate
  often enough that a looped or duplicated quiet passage repeated the excursion and passed stationarity with
  it — measured, 75 of 1200 noise-only files and 68 of 1200 dither-only files. Two frames per stretch takes
  both to zero and costs no sensitivity.

The quiet gate measures the programme with the candidate bands removed, so a hum loud enough to fail the gate
cannot censor its own detection.

`mains` means a mains-COMPATIBLE stationary line was measured, not a causal claim: a synthesised 50.000 Hz
pedal with exact harmonics can be sample-for-sample what an interference pickup leaves. Four limits no
threshold can remove are named in the header, each with the field a consumer reads instead.

### analysis — `LowEnd`: wide bass, and which note owns the bottom

New offline instrument `analysis::LowEnd`, for a master headed to a lacquer. A cutter head writes the
mono sum laterally and the difference vertically, so out-of-phase low end is physics rather than taste.
Two halves, one clock:

- **Wide bass.** An LR4 split at `crossoverHz` (default 120 Hz, `eq::Crossover2`), Mid/Side, and the
  **side energy fraction** `S/(M+S)` of the low band — over absolute 10 ms blocks (the same
  `lround(0.01·fs)` grid `analysis::LoudnessMeter` uses, so the two instruments name the same
  intervals), as a duration-weighted 100-bin distribution, and integrated. A fraction and not the
  requested `S/M` ratio because `L = −R` is a real master that makes Mid *exactly* zero; both raw
  energies are published, so `S/M = f/(1−f)` is one line away. Mono reads exactly 0 and is **valid**;
  digital silence is the 0/0 and reads `NoEnergy` rather than a zero that would look mono.
  Three extremum coordinates, because they are three different questions: the worst fraction (with its
  energy beside it, so an accidental 1.0 can be weighed), the loudest block, and the greatest **vertical
  modulation** — the one a cutting engineer asks for first and which neither of the others identifies.
  The high band's own Mid/Side pair and the unfiltered pair are published too; low + high is an LR4
  allpass, not the input, and the header says so instead of implying additivity.
- **The dominant low note.** 30–300 Hz folded to semitone bands (40 of them at A4 = 440), from the
  shared `analysis::SpectrumFrames`, with the repository's fractional-edge power integration — but
  computed from a precomputed per-band weight table rather than prefix sums, because a fractional edge
  cell's power sits at the midpoint of the *overlap* and weighting it by the bin centre can place a
  centroid outside its own band. One-sided bins are folded, so a full-scale sine inside a band reads its
  own mean square. Both axes are analysed: an anti-phase bass note would vanish from a Mid-only
  spectrum, and that is exactly the programme the first half is shouting about — so the peak band
  reports its own side fraction, i.e. whether the dominant note is cuttable.
  The background is the median **density** of the non-peak bands: semitone bands widen with frequency,
  so under a flat spectrum the top band of the range already holds 5.02 dB more energy than the median
  with no note present. The argmax is published twice (energy and density) and the dominance ratio is
  deliberately not a stored field — its denominator is exactly zero for a tone in digital silence.

Law 8a throughout (bit-identical under arbitrary re-slicing, stronger than law 11(a)): one integer
clock, samples outside and channels inside, `maxBlock` sizing nothing, absolute frame and block grids,
the denormal flush on `core::StateGrid` rather than at the end of `process()`, and no transformed tail —
`finish()` invents no frame and names `tailUncoveredSamples()` instead. Capacity exhaustion is data: the
block series keeps a prefix and says so while every integral, histogram and extremum keeps counting.
`fcore_measure lowend` prints the report as raw IEEE-754 bit patterns.

The suite nulls against oracles computed outside the object — the LR4 prewarped transfer function
predicts the settled side fraction at 120/180/240/480/1000 Hz, a direct O(N²) DFT with an independently
written band integration nulls the fold, and a full-scale sine pins the absolute calibration — because a
re-slicing test compares the implementation with itself and cannot see a deterministically wrong
schedule or a constant calibration error.

<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

### analysis — `BandBursts`: bursts in a band, measured against the band's own past

A new offline analyzer, `felitronics::analysis::BandBursts`, finds where one frequency band (5–9 kHz by
default, the corners are parameters) rises above ITS OWN SURROUNDINGS: the coordinates of each burst, how
far above it rose, and how regularly the bursts recur. It is the honest way to talk about sibilance and
harshness without a target curve — an excess over a programme's own past needs no reference the core does
not have, and it survives the two things that destroy a balance measurement, a different genre and a
different overall level.

**No FFT.** A band is a filter, so the band is `high` of an `eq::Crossover2` at the low corner fed into
`low` of a second at the high corner — two LR4 skirts, zero latency, no transform per hop and no second
source of cross-platform divergence. Group delay at the band edges (≈ 0.14 ms, 1.4 % of a 10 ms hop) is
named in the header and deliberately not compensated.

**The baseline is the MEDIAN of a trailing ring of hops, not the mean**, and that is the decision the
module turns on. A mean baseline switches the detector off exactly where bursts are densest: a steady
train of duty cycle `d` over a silent floor reads `1/d`, so it is invisible once `d ≥ 1/10^(enterDb/10)`
— 25.1 % at the default 6 dB, which sixteenth-note hats are past. Three further failures were measured
during the design round: a 50 %-duty pattern reads 3.01 dB and never fires; moving one onset 10 ms
earlier flips it from found to missed; and 151 zeroed hops make the programme resuming at its ordinary
level read as a 6.1 dB burst. A median is immune to all four, raises the duty bound to 50 %, makes that
bound threshold-independent, and — because it SELECTS one of the observed hop energies rather than
summing them — removes the summation-order and drift questions instead of answering them.

**Periodicity is reported, never judged.** Onsets feed two integer histograms: the spacing of adjacent
events, and the bounded-lag autocorrelation of the onset train. The second is not redundant — a hi-hat
with one hit in five missing turns a clean period `P` into `P, 2P` in the first while the second still
peaks at `P`. Both survive event-list exhaustion. The core publishes the counts and the modal spacing; it
does not publish "this is a hi-hat" or "this is a problem".

**Law 8a**, bit-identical under arbitrary re-slicing (same binary, same channel-presence timeline): one
integer clock, the sample loop outside and the channel loop inside, hop boundaries from a counted
`nextHopEnd_` rather than a modulo, thresholds turned into ratios once in `prepare()` so no logarithm
decides anything, and denormal maintenance clocked by `core::StateGrid` instead of by the `process()`
boundary. `maxBlock` sizes nothing. Law 11d: `storageFor()` publishes the whole demand before anything is
allocated, and `process()`/`finish()` allocate nothing.

Two gates that input sanitising cannot cover are named and closed: `eq::Svf` narrows its state to float,
so a **finite** input can overflow it (successive ±3e38 make the 5 kHz high-pass emit −inf on the second
sample), and a hop whose baseline is **exactly zero** — a silent lead-in — would otherwise open an event
at any level and publish +∞ dB. The filter output is gated as well as the input, and a zero-baseline hop
is not judged and is counted. Events also carry whether damage sat in the BASELINE they were measured
against, because a burst can be manufactured entirely by a hole before it.

Absolute powers in the report (`peakPower`, `bandEnergy`, …) are uncalibrated filter-output power: the LR4
pair's passband peak is −3.99 dB at 6791 Hz at 48 kHz and moves with the sample rate (−1.11 dB at
22.05 kHz, −4.66 at 384 kHz), so they are comparable within one prepared instrument and not across rates.
Every ratio is immune to a GAIN — exactly, over the whole range a delivered programme occupies — and only
partly immune to a change of SPECTRUM: when a burst and its baseline have different shapes each side is
weighted by a different point of that rate-dependent dome, which is worth about half a dB between 44.1 and
96 kHz. Both limits are measured and named in the header rather than claimed away.

`fcore_measure bursts <rate> <channels> <raw.f32le>` prints the whole report as raw IEEE-754 bit patterns,
the way `blocks` does, so a future wasm comparison catches a flipped bit that decimal printing would
round away.

### Added

- **The five offline analyzers are callable from JavaScript, and the wasm module's answer is the native
  tool's answer BIT FOR BIT.** `ProgrammeReport`, `SourceForensics`, `HumDetector`, `LowEnd` and
  `BandBursts` reach `tools/wasm/fc_probe.cpp` through one `_run` entry point and caller-owned row
  buffers with mandatory capacities, and each has a parity harness that reproduces
  `fcore_measure <mode>` exactly — a `diff` of the two IS the test. Measured on every mode at 48 and
  44.1 kHz, against the release module and the checked debug module: **zero differing bytes** — 1299 lines
  of raw IEEE-754 bit patterns at 48 kHz and 1302 at 44.1, where the three extra are burst events the
  lower rate happens to find.

  This is the acceptance the task was written with. It had been weakened to a wasm-only comparison after
  a crew seat measured that byte-exact native-vs-wasm parity was **unattainable** for `hum`, `lowend` and
  `forensics` on the old numerics; `core::det` made the original criterion reachable, so it is the one
  being met.

- **A CI step that enforces it**, five analyzers × two rates × two modules, refusing to count a
  comparison whose output is empty.

### Fixed

- `tools/wasm/build.sh` had no include root for `modules/eq` or `modules/stereo`, which three of the five
  analyzers include from — the module could not compile at all. The same omission on the native side left
  `felitronics_abi_tests` and `felitronics_clips_exposure_tests` linking `felitronics::analysis` when they
  needed `felitronics::analysis_offline`.
- `tools/CMakeLists.txt` carried **four** `target_link_libraries(fcore_measure …)` lines, a residue of
  merging six branches that each added the one they needed. Collapsed to one.

### Notes

- `report` is the only mode whose text is NAMED rather than positional, and its ~100 field names cross
  the ABI **from the module**, produced by the same visitor walk that produced the rows. The visitor is
  deliberately the single enumeration of the report's fields; a name list rebuilt in C++ and again in
  JavaScript would be the second and third copies of it, and the first field added would stop being
  covered without anything failing. The same rule puts `lowend`'s note NAME on the module's side of the
  boundary rather than rebuilding a pitch-class table in JavaScript.

### Added

- **`felitronics::core::det` — a transcendental floor that computes the same bits on every row.**
  `cos`, `sin`, `tan`, `log2`, `log10`, `exp2` and `pow10`, each within 2 ulp of the correctly rounded
  value (measured against mpmath at 60 digits, not against a libm), each pinned against FP contraction so
  the answer does not depend on whether the row has an FMA instruction.

### Changed

- **The five offline analyzers and `analysis::SpectrumFrames` now measure on `core::det`, not on the
  system libm.** They were never reproducible across rows: over the ranges they actually use, `pow10`
  differed in 41 % of results between Apple's libm and both Linux ones, `tan` in 35 %, `log10` in 2 %,
  and the Hann window differed in 4032 of 131072 coefficients at order 17 — which every power bin is
  multiplied by. The odd row was not wasm but **Apple**: glibc and musl agree almost everywhere, so the
  developer's Mac computed something neither CI nor the browser did.

  `fcore_measure report|forensics|hum|lowend|bursts` on the same 10 s programme is now **byte-identical
  across Apple clang/arm64, gcc 14/glibc x86-64 and wasm32/musl** — 1299 lines, zero differences, where
  `lowend` alone differed in 28 lines before.

- **`SpectrumFrames` builds one quadrant of its Hann window and reflects the rest by index.** The window
  was never symmetric — `w[i] != w[N-i]` in 10314 of 16383 places, because `2*pi*i/N` and `2*pi*(N-i)/N`
  are different doubles — on the system libm and on `det` alike, which is why no comparison between them
  could ever have shown it. Reflecting makes the symmetry exact and costs a quarter of the calls.

### Notes

- The deterministic path is for COEFFICIENTS, never for a per-sample loop: a window built once in
  `prepare()`, log2(N) twiddle seeds per transform, a handful of thresholds per `setParams`, one `log10`
  per reported value. Measured cost is +0.5 ms on one order-17 window against a 20-30 ms analyzer run.
  The RT modules (`eq::Svf`, `analysis::KWeightingFilter`) are untouched by this entry.

- **The coefficient math of a filter is now a TYPE, not a flag.** `eq::Svf`, `eq::Crossover2`,
  `analysis::KWeightingFilter` and `analysis::LoudnessMeter` are templates on a math policy, with the
  shipped names bound to `core::SystemMath` exactly as before — TabbyEQ's and OrbitCab's coefficients do
  not move — and `eq::DeterministicSvf` and friends bound to `core::DetMath` for the offline analyzers,
  which is what they now own. `static_assert` pins both directions, so changing an alias is a deliberate
  act that must also edit a test.

  Measured and stated rather than assumed: the two policies give a different `tan` at 7709 of 119880
  filter arguments, and in **zero** of them does that difference reach the filter's float output over
  8192 samples — `Svf` carries its state in float and a one-ulp difference in a double coefficient does
  not survive the rounding. So this routing is DEFENSIVE; the cross-row divergence P79 actually removed
  travelled the double paths (window, `log10`, `pow10`, `log2`). The suite asserts the zero, so the day a
  change makes that path reachable it is a finding rather than a silent regression.

## v0.32.0 — 2026-09-14

### `analysis` — `ClipDetector`: clipping is a flat top, not a loud sample

A new offline analyser that reports whether delivered audio is sample-clipped, and where: every clipped run with its
channel, start, length, sign and level, plus each channel's sample peak and DC offset. Written from the physics of a
clamp rather than ported: a clamp turns the top of every excursion past its ceiling into a run of equal samples, so the
evidence is FLATNESS — consecutive loud samples are not evidence (a full-scale sine and a master normalised to 0 dBFS are
clean), and full scale plays no part (clipped and then turned down is still clipped).

- **The rule.** Bands of samples within 2q (q: the file's quantum as the stream shows it, never below the true one), at a
  local extreme that nothing within ±20 ms passes. A band is a clipped run on RAMP evidence — both sides known; entered
  and left faster than a parabolic crest held inside the band could be (V(L) = (tau+q)·((L+1)/(L−1))² + q, derived, not
  fitted); a coincidence chance ((r+q)/sigma)^(L−1) <= 1e-6 against the flank activity of the UNCLIPPED signal; not a
  step from another flat level; and, for a band that is not exactly flat, met by a line rather than a turn — or on
  CEILING evidence: at the level a ramp run of the same channel and polarity already proved. The step test is the whole
  answer to square waves, pulse trains and sample-and-hold staircases: they jump onto their flat tops from another one.
- **Measured on a bench clamped by construction** (every clamped run known exactly), by run length: 2 samples 99.4 %,
  3–5 99.6 %, 6–20 99.9 %, over 20 99.8 %; every file of the "must report" bucket (84/84) — float, 16 and 24 bit,
  TPDF-dithered, turned down by up to 40 dB, 44.1–96 kHz, one channel clamped, asymmetric clamps, a clamped passage
  inside a louder clean file, non-finite samples. On real recordings clamped by construction across 0.5–20 dB of drive,
  0 to −40 dB of gain afterwards and 44.1–192 kHz: 340 of 342 files, runs of 2 samples 99.7 %, longer 100.0 %. A single
  clamped sample has no flat top and is not found (0.7 %); neither are the few short runs of a polarity that come
  before its first ramp, since a ceiling is only proven forward.
- **Zero false positives** on every clean signal tried: the bench's "must stay silent" bucket (104/104); 890 further
  clean signals (full-scale sines at every sampling phase including exactly symmetric crests, 220 squares and PWMs,
  16/24/8-bit sub-bass with long identical codes, limited and normalised masters, impulses, gates, DC, sample-and-hold);
  52 minutes of synthetic programme; and — the set that set the chance bound — 18 full-length real tracks and stems as
  delivered, normalised to the 16- and 24-bit rails, at 0 / −12 / −24 dB in 16 bit and at −6 dB in 24 bit (126 files).
  Synthetic programme was clean already at a chance of 1e-3; the 72 level-shifted real files still gave 135 false runs
  there, 19 at 1e-4, one at 1e-5 and none from 1e-6. Real crests are often flatter than a parabola (bass through an amp, a compressor); the
  turn test is what separates those from a clamp when the band is not exactly flat.
- **What it does not see, as numbers.** A clamp later passed through AAC (256/128 kbit/s), MP3 (320/128) or a resampler
  (to 44.1 or 96 kHz) is no longer flat: 0 files of 18 per chain, 0 % of runs. A 20 Hz DC blocker after the clamp: 10 of
  18 files, ~5 % of runs. Noise added after a 16-bit clamp: runs of 3+ found 96 % at ±1 LSB, 50 % at ±2, none from ±4.
  A pure tone locked to the sample grid can hide its clamp (the stream never shows it a slow step, so q stays loose).
- **Against the detector it replaces**, on the same bench: that one reports 48 of the 104 clean files and misses 3 of
  the 84 clamped ones; of the 94 files where the two verdicts differ, 51 have a decided answer and this engine is right
  on all 51 (the other 43 are the bench's open questions). On the real recordings normalised to the rail it reports 14
  of 36 files; this engine none.
- `process()` is read-only, allocation-free and bit-identical under any slicing; `finish()` decides the last 20 ms; no
  libm call decides a verdict (the chance is a product, not a logarithm), so every tier agrees. A whole-file reference
  written straight from the rule nulls the streaming engine run for run (randomised material and 4000 tiny streams in
  the suite, 2263 corpus files out of it); a 41-mutant stand kills 40, the survivor equivalent. ASan found a ring
  overrun no other check could see — a long flat band trimmed the window deque only at queries — now expired on every
  push. A code-review round with one mandate, "find an input that reads or writes past a buffer or carries state
  across calls", found none: an instrumented copy asserting every ring index (the in-vector overrun ASan cannot see)
  over 3024 lifecycle cases, and 15 million frames of stress input, reached both ring capacities exactly and never
  past them (deque W+1 of W+2, pending ceil(W/2) of W/2+2). It did find the report accessors unchecked — a channel
  or run index outside the report, or any call before prepare(), read past a vector; they now answer empty.
  ~13 ns/sample on arm64 (a heavily clipped 5-minute stereo file in 0.38 s).

### `core` — a delivery resampler: rational L:M, 146.58 dB of stopband measured through it on all 30 pairs

A new primitive, **`core::DeliveryResampler`**: the rational polyphase windowed-sinc converter for the rate a
file is DELIVERED at — 44.1 / 48 / 88.2 / 96 / 176.4 / 192 kHz, every pair, both directions. Neither
neighbour was touched or widened: `StreamResampler` interpolates between phase-table rows and floors near
-105 dB whatever its kernel, and `resampleIr` is an IR-load one-shot. Non-integer rates are refused, not
approximated.

- **The bar, measured through the converter, not read off a formula.** Stopband >= 140 dB and ripple
  <= 0.001 dB to 20 kHz were asked. Delivered, over all 30 directed pairs: worst line **-146.58 dB**, worst
  ripple **0.0000013 dB**, worst absolute gain error **0.0000014 dB**; a 48 -> 44.1 -> 48 round trip nulls
  at **-143.3 dBFS peak / -160.3 RMS** against an analytic oracle (asked: -120). The Kaiser formula does
  not certify at this depth (asking 140 measures -139.3), so the design target is 146 dB and the transition
  ends 6 % of its width INSIDE the Nyquist: without that margin the matrix was not even monotone in the
  target — asking 2 dB more took the worst line from -142.0 to -139.0, through the bar.
- **Routes are costed, not written down.** A stage's work is `K * f_in * f_out / df` with no L in it, so a
  cascade's only lever is where the narrow transition runs. 192 -> 44.1 goes **192 -> 88.2 -> 44.1**
  (26.99 MMAC/s, 12,792 coefficients) against 42.34 MMAC/s and 141,120 coefficients single-stage — and not
  the textbook 192 -> 96 -> 48 -> 44.1 (29.96). 176.4 <-> 192 is where decomposition has nothing to give:
  it goes direct, because every lower rung discards band that pair must carry. At most two stages anywhere.
- **The passband edge is proportional** — 20.000 kHz whenever 44.1 is involved, 80 kHz on 176.4 <-> 192 —
  and a parameter, because it is the one number that changes what the output contains.
- **Double in the sample loop, as a named law-3 carve-out.** Narrowing only the four accumulators to float32
  fails the bar on 13 of the 30 pairs (worst -132.19 dB), so `firDot` is not reused; bit agreement across
  rows is not promised.
- **Latency is an exact rational**, each stage an integer number of its own input samples (N = 2Lk + 1), and
  it is measured back out of the carrier phase on single stages and cascades alike. `flush()` drains it;
  law 11 holds clause by clause (long calls chunk, a stopped channel drops its memory, a gap spends the
  clock, `process()` after `flush()` is refused).
- **Two suites, and the instrument is itself under test.** The spectral bar suite carries four negative
  controls — including an off-bin -138 dB image that the first, rectangular-window version of the suite read
  as -141.3 and passed. A mutation stand of 33 mutants catches all 32 non-equivalent ones; the survivor, a
  redundant range test, and an earlier pruning mutant are proven equivalent (144,360 and 624 plans
  bit-identical with and without them).

### `analysis` · `tools` — the waveform bars and the stereo band: one definition, two roads

The site draws two pictures from audio, and each had several definitions: the waveform "peaks" came from a Java
sidecar generator (ffmpeg to 8 kHz **s16**, integer buckets), a python script, the page's own JavaScript, and a
client-side costume that re-scaled demo sidecars by their true peak; the stereo band's correlation had the page's,
`fcore_measure`'s `long double` one and `CorrelationMeter`'s (a different quantity). A demo file and an uploaded one
could be drawn by different definitions in the same interface, and both pictures looked plausible. This is the core
half — the one definition and its two roads; moving the site's generators and page onto it is the site's work.

- **`analysis::WaveformPeaks`** — a port of `computePeaksFromBuffer` / `peaksFromWav` (audio-peaks.js), bit for bit:
  box-average decimation to ~8 kHz (`Math.round`, ties up), then max-abs per bucket; mix modes `avr` · `L` · `R`
  (the LAST channel) · `max`; the double output and the float32 form, because the two JS functions differ exactly
  there. **Not a metering peak**: not above the sample peak except by the rounding of a box mean, and neither of the core's
  two true peaks.
- **`analysis::StereoColumns` / `StereoSums`** — a port of `computeStereoColumns` / `correlationOf` / `widthOf`
  (stereo-meter.js): per column width, **uncentred** phase correlation (not Pearson's, whatever the JS comment
  says — (1,2)/(2,1) reads 0.8) and RMS (the JS `loud`; not a loudness), `maxRms` as the unrounded double the page's
  verdict thresholds use, and the playhead needle over any stretch. The verdict stays on the page.
- **The spec's known properties are kept, each pinned by a witness the site's own JS computed in node**: the last
  waveform bucket is never emitted at decLen 2007 / 1000 buckets (the boundary is 2007.0000000000002); a column
  boundary is floor(i·(len/cols)), so frame 200 of 1206/1200 is in column 200 and the last frame can be in none; a
  partial box is dropped; short files zero-fill buckets but shrink the column count.
- **Contraction:** the stereo products go through `volatile` stores (law 10 now names the pin): `mid += m*m` fused on
  arm64 reads a width of …0010 where the JS reads …0012. The pin was measured on Apple clang 21 arm64, gcc 14.2 x86-64 and gcc 14 arm64 under `on`/`fast`; the removed pin fuses where each compiler fuses, and the witnesses see it there.
- **Roads:** `fcore::ShapeProbe`, shared by `fcore_measure waveform|stereo|needle` and the new `fc_probe_shapes_run`,
  `fc_probe_waveform_*`, `fc_probe_stereo_*`, `fc_probe_needle` exports (names that cannot be read as
  `fc_probe_sample_peak`, `fc_probe_tp_linear` or `fc_probe_lufs`); existing `fc_probe_*` unchanged. The wasm export
  list is now generated from the source, as fc_master's is. A CI step diffs native against wasm (release and checked)
  over every mix mode, bucket and column count and a mono / stereo / six-channel fixture.
- **`fcore_measure correlation`** is the stereo band's binary64 formula over the whole file; law 9's one sanctioned
  `long double` exception is gone, and the artifact gate lost its one named exclusion.
- **Measured (out of tree, `.private/harness/p59a-shapes/`):** NULL site-JS vs C++ (`-ffp-contract=fast`) vs wasm on
  560 synthetic, 35 real and 79 decoded gate items — 0 differences, 0 split mismatches. Demo road vs upload road:
  lossless WAV/FLAC at the native rate bit-identical after PCM equality (46 of 53 files; the 7 others are float64 WAV,
  which the site's WAV reader reads as zeros — a site finding); 21 lossy decoder pairs all inside the bound derived
  from their own per-sample difference, including fixed-point decoders clamping a +20 dB master at 1.0. Mutation
  stand: 31 mutants, 28 killed (the two-statement removal of the product pin only by the out-of-tree NULL); the three
  survivors are two equivalent guards and a null check whose removal is undefined behaviour only UBSan sees.

### `mastering` — parallel compression: `MasteringChainParams::compressorMix` (P60)

The compressor stage's output is now `(1 - mix) * dry + mix * compressed`, where `dry` is the stage's own
input delayed by exactly the compressor's lookahead through a third `core::DryAligner`. It is the chain's
field and not the compressor's: `dynamics::Compressor` keeps dry/wet out on purpose.

- **mix = 1, the default, is the chain before the field, bit for bit.** Pinned in the suite against a bare
  compressor and the whole hand-built chain, and measured once against a render built at `bddb303`: seven
  scenarios, every output sample and every tap, identical. `mix = 0` is the input delayed by the lookahead,
  bit for bit, found from the input.
- **The blend is computed in double.** The float spelling is a multiply-add the compiler may fuse, and the
  desktop tier builds with `-ffp-contract=on` while wasm builds with `off`: measured, the float form gives
  different bits fused and unfused at ten of eleven mixes. In double both products are exact for
  `mix >= 2^-6`, so every row agrees.
- **A step at the quantum boundary, not a ramp** — the same contract as the two gain nodes. It keeps
  blending while the compressor is bypassed: skipping it there dropped the level by the whole gain
  reduction on the quantum the bypass engaged (20.2 dB at mix 0), and a steady bypass is bit-identical
  either way. `compressorGrDb` is still the compressed path's gain reduction at any mix.
- **Memory:** `4 * channels * (K + lookahead + 2)` bytes more per chain with a compressor — 2 448 B for
  the stereo 48 kHz default. `DryAligner::Storage::freshBytes()` is new: a one-channel aligner whose ring
  fits the constructor's seed asks for nothing, which `bytes()` over-stated by 8 B.
- **The C ABI does not carry it yet.** It waits for the ABI's version rule (P57); every params set through
  the ABI renders at mix 1. `tools/fc_master_abi.h` records the omission as a debt.

### `tools` · `mastering` — the C ABI's compatibility rule (v2, v3): delivery at another rate, and parallel compression

`tools/fc_master_abi.h` states the rule every later ABI change follows, so that a version is a row in a table
rather than an event, and makes the first two bumps by it.

- **The rule** (VERSIONING, rules 1-8). One version for the whole ABI, moved by a struct that grows or an entry
  point that is added. A bump appends fields to the END of structs that begin with a header, each IN field with a
  default under which the call is the previous version bit for bit; nested structs and header-less array elements
  never grow. Tail padding is a named field. The size of every (struct, version) is ONE row of
  `FC_MASTER_STRUCT_SIZES`, published by `fc_master_sizeof(id, version)`. An IN struct of any known version is read
  over this build's defaults at the CALLER's size, and every bound and alias check uses that size; a version newer
  than the build is refused by decision. An OUT struct keeps the caller's header and gets exactly its size.
- **The frozen writers.** `fc_*_default(out)` cannot know the caller's size, so they are frozen at v1 (stamp v1,
  write v1's bytes); `fc_*_defaults(out)` write the caller's stamped version. A field set after a frozen writer
  lies past the stamp and is not read — the JS accessor refuses to write or read such a struct at all, and the
  CLI, the parity harness and the suites moved to the versioned writers. `fc_master_sizeof_params/config` are
  frozen at v1 too, so a v1 page fails on the version, not on a size.
- **v2 — `fc_master_config::deliveryRate`.** Non-zero makes a DELIVERING handle: SRC first, chain and solver at the
  delivery rate, equal rates included. New entry points `fc_master_delivered_frames`, `fc_master_render_delivered`
  and `fc_master_solve_delivered`; `fc_master_measure_lra` converts on such a handle; `process`/`flush`/`solve`
  refuse. The composition is one core class, `mastering::DeliveredMastering`, which the facade and the selftest's
  direct path both call. Budgets: `FC_NEED_SOLVE` / `FC_NEED_MEASURE_LRA` on a delivering handle include the
  converted programme and judge a range on the DELIVERED length; `render_delivered` allocates nothing.
- **v3 — `compressorMix`** (P60's field) at the end of `fc_master_params` and `fc_master_resolved`. Clamped by the
  core and read back; a non-finite mix is refused.
- **A non-finite input sample costs a delivered render what it costs a plain one.** `DeliveryConverter` now gates
  every input sample with the chain's own rule (NaN/inf -> 0, clamp +-1e6) before the conversion and counts it:
  a windowed sinc used to spread one NaN over its kernel, and the chain's gate then zeroed every delivered sample
  it reached. Converting a NaN is bit-identical to converting a zero in its place.
- **Refusals move nothing.** `DeliveredMastering` reaches every verdict — renderer, width, overlapping planes, and
  everything the solver would refuse before a pass (`TargetLoudnessSolver::admits`, factored out of `solve()` with
  no change of verdict) — before the converter writes a sample or a programme buffer is allocated.
- **Proof.** A render, a solve and a range through a v1 build (`179e1c5`) and this one are byte-identical on the
  same fixture; v1- and v2-stamped structs render the v3 caller's bits at the previous values; the delivered render
  and search are bit-identical to the core on down, up and equal rates; `need` equals the counted allocation;
  `tools/wasm/layout-check.mjs` holds every JS field offset against the compiler (wasm32 on that tier) and
  `layout-gate-test.mjs` the JS gates; `master-parity.mjs` compares the delivered paths wasm vs native at the
  delivered length. Product check, 30 s of music: 96 -> 44.1 kHz -14.003 LUFS / -1.044 dBTP, 192 -> 48 kHz
  -14.003 / -1.046, measured by `fcore_measure` at the delivery rate.
- **Site transition.** The rule cannot reach a page already shipped against v1: its loader requires `version === 1`.
  The site moves to v3 as one release of the worker, its layout copy and the module.

### `mastering` · `tools` — the gain-reduction trace: WHERE the compressor and the limiter worked (ABI v4)

A loudness search used to report its gain reduction as three numbers for the whole programme — mean, p95, max. The
solution now also carries each stage's TRACE: the delivered programme cut into uniform buckets, and per bucket the
largest and the mean |GR|.

- **`mastering::GainReductionTrace`** in `LoudnessSolution` (`compressorTrace`, `limiterTrace`), built by
  `GainReductionTraceBuilder` in the solver's tap sink — the same tap values, windows and units the statistics use
  (frames for the compressor, frames × `tapOversampleFactor` sub-samples for the limiter), so where the statistics
  are valid the maximum over the buckets IS their `maxDb`, bit for bit. A maximum per bucket, not a point sample.
  `min(1000, frames)` buckets over the programme's frames, boundaries `floor(k·F/B)`, never an empty bucket; per bucket
  a sample count and a non-finite count, so a poisoned or empty stretch cannot read as "idle".
- **It describes the audio in `out`.** Reset and rewritten on every render, including the delivery re-render of a
  search that did not end on its best point; `valid` only when the render ran to its end, the window saw a sample and
  none was non-finite; no buckets for a verdict reached before any render. It lives in the solution, so a later solve
  does not touch it. The price is the solution's size: 2336 → 50 384 B per solution record, pinned with its formula.
- **ABI v4** (the rule's third bump): the new entry point `fc_solution_gr_trace(solution, stage, out, cap, written)` on
  the pattern of `fc_solution_log`, the header-less `fc_gr_trace_bucket` (frozen from v4), the `fc_gr_stage` codes,
  and the traces' bucket counts and validity at the end of `fc_measurement` — one row of the size table (224 B from v4).
  `fc-master-layout.mjs` is at v4; the site's copy of it must follow before a page reads the trace.
- **Verified:** the trace nulls bit for bit against the chain driven by hand and bucketed another way (three lengths,
  including one shorter than 1000 frames); the ABI's traces are the core's bit for bit — on a plain handle with both
  stages live (MasterAbiTests) and on the delivered-rate path (`fcore_master selftest`, whose fixture leaves the
  limiter idle) — and native agrees with wasm through `master-parity.mjs`; the re-render branch, the no-re-render
  branch, an early exit, the refusals and a failed first render each carry the trace their render left; an impulse's
  GR sits in its own bucket and ends where the limiter's hold plus release puts it. Built and run with MSVC as well.
  Mutation stand: 22 mutants, 18 killed; the four survivors are equivalent on every reachable input (the two stages'
  bucket counts and validity are equal by construction, the one reachable `RenderFailed` refuses before any tap, and a
  defensive branch for a stream that goes backwards).

### `analysis` · `mastering` · `tools` — one true-peak instrument aims and certifies a delivered ceiling

`TargetLoudnessSolver` promised a delivered file `<= maxTruePeakDbTp`, judged every render with
`analysis::TruePeakMeter`, and the file was certified by `fcore_measure`'s reference filter — two different
instruments. On bright material, and at the high delivery rates where the short filter stops interpolating, the
first read under the second by more than the solver's whole 0.05 dB aim, so a render the solver called feasible
was delivered over its promise. Fourteen real programmes delivered at six rates, outside the tree: 12 of 84
deliveries certified above -1 dBTP before, 0 after.

- **`analysis::ReferenceTruePeakMeter`** (new). The reference true peak as a module class: `PolyphaseOversampler`
  at 4x / 32 taps per phase at every rate, the maximum of |x| over the 4x stream floored at the sample peak, and
  `drain()` for the FIR's tail. Law 11 in the house order; a channel that stops is DRAINED at that moment rather than
  dropped — its pending peak was submitted and belongs to the reading, and it returns from silence — so a zero-width
  call is a pause that drains every channel. `prepare()` refuses a non-positive block as well as a bad rate or width;
  `storageFor (rate, maxBlock, channels)` equals what it allocates; nothing is allocated in `process()`/`drain()`.
  `analysis` now links `oversampling` (which depends only on `core`).
- **`fcore::Probe` measures through it.** `fcore_measure` and the browser's `fc_probe` report the same numbers,
  byte for byte (`--precise`, 40 runs over a real corpus against the previous binary). One sequence behaves
  differently, and neither shipped caller makes it: a channel left out of a narrower call and then given audio again
  no longer replays its pre-gap history. A narrowing stream that simply ends reads what it read before, to the bit.
- **The solver reads every render with the reference.** The ceiling is aimed, feasibility judged and
  `measured.truePeakDbTp` reported on the certificate's arithmetic: the reported number is now the certificate of
  the delivered samples, bit for bit (above the dB floors — silence is still spelled -200 dB, from the same float
  threshold as before).
  - **Behaviour:** a render that only the old meter called feasible now costs one more pass (10 extra passes over
    the 84 real deliveries), and the loudness it reaches is unchanged (worst move 0.0004 LU); a pass costs 10-15 %
    more (the reference filter on the whole programme, 60 s stereo, 44.1 -> 192 kHz).
  - **Budgets:** `solveBytes()` — and so `FC_NEED_SOLVE` — is the loudness meter plus 21 008 B for a stereo
    reference meter at any rate, where it was the loudness meter plus 392-296-248 B of `TruePeakMeter` plus a
    512 B drain buffer. `TargetLoudnessSolver::kDrainFrames` is gone: nothing is drained from a buffer any more.
- **`felitronics_truepeak_instrument_gap_tests`** (new) owns the number "how far apart the two meters read": five
  materials (music, drums, bright noise, a 16 kHz burst at its worst phase, a click at its worst offset) delivered
  at 44.1-192 kHz through the chain and pinned per cell. The worst is the drums at 44.1 kHz, where the cheap meter
  reads 0.2995 dB under the reference and 0.3267 dB under the band-limited truth (read at a continuous time, not on
  a zero-padding grid); the burst rows are the grid's closed form at 96 and 192 kHz, and at 176.4 kHz and above the
  cheap meter is shown to be a sample-peak meter. The reference is not the truth either, and that is pinned too:
  0.0272 dB under it on the drums, 0.3270 dB under it on a click flat to 0.45 fs at 48 kHz.
- **The promise names its instrument.** `LoudnessRequest`'s documentation now says what "-1 dBTP" means: the
  ceiling as `analysis::ReferenceTruePeakMeter` (and so `fcore_measure`) reads it. Like any BS.1770-class meter the
  reference under-reads the band-limited peak — up to 0.33 dB on a full-band click, pinned in the gap suite — so a
  third-party meter may read a delivered file above the promise. Kept on purpose: changing the certifying instrument
  would move every certificate already issued.
- **`felitronics_delivered_ceiling_tests`** (new): every source rate to every delivery rate of the six solves,
  certifies at or under the promise, and reports the certificate exactly.

### `mastering` — `LoudnessSolver.h` and `DeliveredMastering.h` enter the strict header gate

Both headers were outside `felitronics_header_hygiene`: the solver did not compile under the downstream flag set,
and `DeliveredMastering.h` includes it, so every mastering header built on top would have stayed out too. gcc 14.2
reported seven diagnostics, all in `LoudnessSolver.h`; `DeliveredMastering.h` added none of its own.

- **Six `-Wfloat-equal`** — the "no limit" sentinels (`GainReductionLimit::off()`, the `minPlrDb = -inf` and
  `maxLraLossLu = +inf` tests in `worstExcess()` and `violatedMask()`) and the infeasible tie-break on equal excess
  — now go through `core::exactlyEqual`, like the rest of core. `!=` became `! exactlyEqual`, which is the same
  predicate for every IEEE value, NaN included.
- **One `-Wconversion`, `uint64_t -> int` for `MasterMeasurement::nonFiniteSubHops`** — explicit cast, not a
  defect. The meter in `measure()` is a local, zeroed by its own prepare, fed exactly one `process()` of `frames`
  samples, and its counter moves by at most one per sub-hop of at least one sample, so it cannot exceed `frames`,
  an `int`. The comment at the cast says so, and names what would break it.
- **Every public header is now in the gate BY NAME.** `BlendKernels.h`, `BlendParams.h`, `IrBlend.h`,
  `MatrixConvolverNupc.h` and `StateGrid.h` were reached only through other headers' includes, so dropping one
  intermediate `#include` would have taken them out of the gate silently. None warned when named.
- **Proof.** The gate builds warning-free on gcc 14.2 (Debian), Apple clang and emscripten 6.0.9; the full suite
  passes. Behaviour is unchanged by construction: a probe that renders, solves and measures LRA (plain, constrained,
  infeasible, delivered 44.1 -> 48 kHz) prints byte-identical results against `main`, and its `.text` section
  compiled against `main` and against this branch (gcc 14.2 `-O3`) is byte-identical.

### `mastering` · `tools` — one rule for the planes a whole-programme operation reads and writes; `fc_solution_log` refuses a `written` inside its records

- **One rule, `mastering::planesUsable` (new `Planes.h`).** Tables and planes non-null; no input plane's bytes
  touching any output plane's, every pair, half-open, EACH SIDE AT ITS OWN LENGTH (`inFrames` for the input,
  `outFrames` for the output — a conversion's two differ, and judging both at one of them misses an overlap that lies
  only in the longer span or refuses planes that never meet); and no two output planes touching. The loudness search
  asks it with one length for both sides, `DeliveryConverter::convert` and `DeliveredMastering` with a conversion's
  two. Still legal: one buffer feeding two input channels, planes edge to edge in one allocation, a call shorter than
  its buffers, buffers reused across calls — each pinned bit for bit against buffers of their own.
- **The loudness search refused only `in[c] == out[c]`.** `out[0] = in[1]` was accepted: the render wrote channel 0's
  master where the next pass reads channel 1, and the call returned an ordinary verdict at a plausible gain over a
  master that is not the programme's — in the suite's witness the aliased call answered `TargetUnreachable` where the
  honest solve is `Solved`, at 12.2532 dB against 12.3175, reporting −10.072 LUFS against −12.093, with a delivered
  master different in every one of 144 000 frames. `out[0] = out[1]` was accepted too: `Solved` at 12.3041 dB against
  12.3175, with only the second channel's render left in the buffer where two channels were asked for. Both are
  `InvalidRequest` now, before a render. A search that happens to end after ONE render over cross-aliased buffers was correct, and is refused too:
  whether it is correct would depend on how many passes it took, the reason `in == out` was already refused. The
  direct C++ call now refuses what the facade refuses on the same memory (`MasterAbiTests`). **Also:** a single null
  plane beside good ones used to reach the renderer and crash; it is `InvalidRequest` now, from the same rule.
- **`DeliveryConverter::convert` checked its planes for null only.** It writes `out` at the delivery stride while still
  reading `in` at the source one, so an output plane over an input plane overwrote programme not yet read: at
  44.1 → 48 kHz `out[0] = in[1]` returned true with channel 1 wrong in 25 990 of 52 245 frames of the suite's
  witness; at 48 → 44.1 kHz it came out right only because there the writes lag the reads. Refused now, in both
  directions, before a sample is written — so is an overlap that exists only in the longer of the two spans, and two
  output planes on one buffer. **A behaviour change for a direct caller:** an identity conversion (equal rates) with
  `in[c] == out[c]` copied the bits correctly in place and is refused too — whether an overlap is safe would otherwise
  depend on the ratio; `DeliveredMastering` and the C ABI already refused it, and nothing in the tree calls it so.
- **`DeliveredMastering::render` of an empty programme skipped the plane rule** and then read the output table:
  `render (…, nullptr, 0, nullptr, 0)` was a SIGSEGV on a call whose programme is legal. The rule now applies at
  every length; at two zero lengths it judges only null, so an empty programme in real tables still renders, and
  one with null tables — or with null planes in them, which used to be accepted — is refused.
- **The non-finite input count is THE LAST CALL'S THAT REACHED THE COUNT** — stated in the same words for
  `DeliveryConverter::nonFiniteInputSamples`, `DeliveredMastering::nonFiniteInputSamples` and
  `fc_master_stats::nonFiniteIn`, and pinned by tests at all three. A call refused before its count (by the facade or
  the core) leaves the previous number, the rule `fc_solution_log` keeps for `written`; one refused after it keeps its
  own (at equal rates a range measurement refuses a poisoned programme having counted it). To make the words exact:
  the converter publishes its count once every input sample is read, and a conversion that does not complete no
  longer leaves a partial count behind. Not a version.
- **`fc_solution_log`: `written` may not point into the `cap` records** (FC_ERR_SPAN), and is cleared only once
  every check is behind the call — the order `fc_master_flush` takes. It used to be cleared on entry, so a
  `written` inside the buffer took the count over a copied record on FC_OK and a zero into the buffer on a
  refusal. **A behaviour change for a caller that read `written` after a refused call:** a refusal on `out` (null,
  alignment, span) now leaves it as it was instead of zeroing it; `cap == 0` is FC_OK with a zero, as before. The
  refused span is the whole capacity, not the records a given log fills. Not a version: this change does not move
  `FC_MASTER_ABI_VERSION` (VERSIONING rule 1 — no struct grew and no entry point was added).
- Comments that said what the code does not: `renderTapped` claimed the chain's drain produces no gain reduction
  (it carries the release, and an expanding or upward mode acts on its silence; the windows exclude it, which is
  why the numbers were right); `fc_master.cpp` claimed a field retyped or inserted mid-struct is a build error (only
  where the change moves an offset or a size the pins read: an `int32_t` dropped into padding, `int32_t` to
  `uint32_t`, or a `double` narrowed to `float` before another `double` all build outside v4's type-pinned fields,
  and layout-check compares offsets, not types). And `fc_master_abi.h` said a count out-parameter is cleared FIRST,
  which no entry point with one does any more — including the sentence under `fc_solution_gr_trace` that contrasted
  it with `fc_solution_log`.

### `convolution` — the IR resampler stops inventing a cabinet's top octave; a one-tap IR and a nearly-equal rate load

`convolution::resampleIr` treated what lies past either end of an impulse response as the weighted mean of the
samples it had, not as silence: a tap outside the input was skipped before its weight was added, so every edge
sample was divided by only the part of its window that landed on the input. A cabinet IR starts at its onset, and
`CabConvolver` resamples a cabinet whose file is at another rate than the host — the ordinary case — so a 48 kHz
cabinet in a 44.1 kHz session played a broadband floor over its own top octave. Twenty-one factory cabinets,
measured outside the tree, worst 1/6-octave band against each cabinet's own response: at 44.1 kHz +7.2 dB at 16 kHz
and +25.1 dB at 18 kHz, now +0.08 and -0.45 dB; at 96 kHz +18.5 dB at 18 kHz, now -0.12 dB.

- **Zeros outside the input.** Every output sample divides by the weight of its whole window; a tap past either end
  adds weight, not signal. Resampling `[zeros, x, zeros]` now equals resampling `x` shifted by whole samples. What
  zeros cannot give back is the kernel's pre-ringing that would fall before output sample 0 — an IR that starts at
  sample 0 keeps it only by adding delay. It is a cut of at most 0.45 dB at 18 kHz on the factory set and grows with
  how abruptly an IR starts: the test's cabinet-like fixture reads +1.23 dB at 20 kHz at 96 kHz.
- **At least one sample.** The length is `inLen * ratio` rounded but never zero, as JUCE's `resampleImpulseResponse`
  had it. A one-tap IR at 96 -> 44.1 kHz rounded to nothing, and `CabConvolver` then published nothing without a
  word: the previous cabinet kept playing. The result is empty now only for no input (a null pointer or a
  non-positive length), a rate that is not a positive finite number, a length past `INT_MAX` — refused before the
  cast, where 2^32 + 1 used to wrap to ONE sample — or a ratio so small that output sample 0's position is past
  `INT_MAX` (the one-sample floor is what made that `(int) floor(t)` reachable; it is refused rather than
  undefined).
- **`CabConvolver::kRateMatchTolerance`** (new, relative `1e-6`, the same as orbit-amp's `sameRate`). Two rates this
  close are one rate and the IR loads verbatim; the exact comparison sent a host reporting 48000.0000001 through the
  resampler, which band-limits at 0.95 of Nyquist and moves every tap. Measured outside the tree on two factory
  cabinets, against a 2048-tap resample: played verbatim 1 ppm off its rate, one reads -87.5 dB over its first 100
  ms and the other -72 dB over its first second, where running the 64-tap resampler costs -73 and -81 dB; at 10 ppm
  verbatim is the worse of the two.
- **A load that stages nothing is ignored whole.** `buildAndStage` overwrote the retained taps before it knew the
  load was empty, and returned after. With a load still pending from mid-crossfade, the retry then used the old
  length on an emptied or narrowed store: a read past the end of a vector (a zero-length load, or a mono load over a
  pending stereo one — ASan container-overflow, libc++ hardening abort), or, with a null data pointer, a retry
  refused forever, `isBusy()` stuck true and neither the pending IR nor the new one ever reaching the convolver. The
  taps, their gain and the pending geometry are now staged in locals and committed together; a zero or negative
  length, a null channel array or plane, or a known rate so far off that `resampleIr` cannot address the result (its
  length or its last position past `INT_MAX`) leaves the playing IR, `stagedTaps()`, the normalization gain and any
  pending retry exactly as they were.
  - **Behaviour:** after such a load `stagedTaps()` still holds the taps of the last load that staged any (it used
    to be emptied) — what the convolver plays once a pending retry has published. A null plane with a positive
    length is ignored instead of dereferenced.
- **An unknown IR rate loads the taps as is — one rule.** A rate that is not a positive finite number (NaN, zero,
  negative, ±inf) is unknown: the samples in the file are fine, only the metadata is broken, and refusing would drop
  the cabinet and play silence where as-is at worst plays an impulse of the wrong length. NaN, zero and negative
  already loaded as is on `main`; +inf went to the resampler, came back empty and the load was dropped without a
  word.
- **Tests.** `felitronics_convolution_resampler_tests`: shift invariance at both edges over ten rate and radius
  cases and five inputs (in front to 1e-6, behind to the bit — on `main` 1437 and 1053 misses, 3.25e-2 worst on the
  cabinet fixture), the edge taps against an independent long-double recomputation of the specification, a
  hand-worked 105/104, a cabinet-like IR's bands against its own response and against the untruncated resample, the
  rounding and the one-sample floor at seven short IRs, and the refusals (zero, negative and non-finite rates among
  them). `felitronics_convolution_cabconvolver_tests`: the tolerance witness on literal rates, so the one part per
  million is pinned from both sides (verbatim to the bit at 48000.0000001, at 0.9 ppm and at 0.5 ppm of 192 kHz,
  both ways; resampled at 1.1 ppm, both ways, every tap moved), an unknown rate (NaN, zero, negative, ±inf;
  normalized and not) loading as is and playing, a one-tap IR off-rate staging, publishing and playing at four rate
  pairs, and nine loads that stage nothing over a pending one — taps, gain and what plays, normalized and not — plus
  a one-tap load that wins as the latest.

### `tools` — release notes move to `changelog.d/`, one file per task

Every branch that changes behaviour now writes its note as a NEW FILE under `changelog.d/` instead of appending to
`## Unreleased`. Git cannot conflict on two branches adding two different files; it conflicted on that one section
five times in a single day, and each time cost a full rebase-build-push-wait round on work that was already green.

`node tools/changelog-collect.mjs --release vX.Y.Z` folds the fragments into `CHANGELOG.md` under the new heading and
deletes them — one commit on the release branch, where there is nobody to conflict with. `--preview` prints what the
next release would say, which is how the accumulated notes stay readable now that they live apart.

## v0.31.0 — 2026-09-12

### `core` · `oversampling` · `analysis` — one polyphase FIR kernel for the whole tree, and five rows that agree on its bits

The same inner loop was written three times — `PolyphaseOversampler::upsample`, the same class's `downsample`,
and `TruePeakMeter::process` — and between them they were **~83 % of a mastering render** (measured on a
600 s programme: `downsample` 50 %, `upsample` 33 %). They are now one function, **`core::firDot`**, with
four hand-written kernels (scalar / SSE2 / NEON / wasm-SIMD128) that compute the identical summation:
four partial accumulators, `(s0+s1) + (s2+s3)`, every multiply and add rounded separately.

- **The work was the REPACKING, not the kernel.** The three loops were not the same loop: two gathered
  coefficients with a stride of `L`, one read them contiguously, and all three walked a modulo ring backwards
  with a wrap test inside the inner loop. `prepare()` now stores the coefficients **phase-major** where the
  gather used to be, keeps the sample history in a **double-length backwards ring** (every sample written
  twice, so the window is always contiguous and already in the coefficients' order), and pads the run to a
  multiple of four with `+0.0f` so there is no tail. What is left is a dot product of two contiguous spans.
- **Speed, on the same programme.** A 60 s stereo render: **1.51 s → 0.48 s (3.15×)** on arm64 macOS
  (NEON), **2.75 s → 1.09 s (2.52×)** on x86-64 Linux/gcc 14.2 (SSE2).
- **Bit-identity across five rows, for the first time and now gated.** The kernel returns one constant on
  `win` (MSVC 19.44), `deb` (gcc 14.2), `mac` (Apple clang 14.0.3), `docker --platform linux/arm64`
  (gcc 14.4) and `wasm` (emsdk 6.0.9 in node, scalar kernel and `-msimd128` alike) — and so, measured, do
  `PolyphaseOversampler` and `TruePeakMeter` end to end, libm-designed Kaiser coefficients included.
  `felitronics_core_polyphasefir_tests` pins that constant, so a row that loses the property goes red.
- **It overrides law 10 locally, and that is the point** (`docs/DSP-ARCHITECTURE.md`): `acc += a*b` is the
  contractible form, arm64 fuses it, baseline x86-64 cannot. Before this change the two rows agreed only by
  accident — `std::vector::operator[]` happened to block gcc's contraction under the tree's stated
  `-ffp-contract=on`; under gcc's OWN default (`fast`, what any consumer TU of these INTERFACE targets gets)
  arm64 and x86-64 **did** diverge, measured on all three hashes. They no longer do.
- **The three pragmas were measured, not looked up.** gcc ignores `#pragma STDC FP_CONTRACT` in C++ and
  needs `#pragma GCC optimize("fp-contract=off")` (which costs 1.4 % by blocking inlining); clang needs
  `#pragma clang fp contract(off)` and must NOT be given `#pragma float_control(precise, on)`, which turns
  contraction back **on**; MSVC takes `#pragma fp_contract(off)`. A clang build with `-ffast-math` defeats
  all of them and nothing in a header can stop that — hence the gate.
- **What the claim does not cover, named rather than hidden:** FTZ/DAZ (wasm cannot flush at all, and a
  normal × normal product can land subnormal, so a flushing host splits native from wasm on NORMAL inputs),
  the rounding mode, and NaN sign/payload — all runtime state, none reachable by writing the loop differently.
- ⚠ **`PolyphaseOversampler`'s and `TruePeakMeter`'s output moves bit-for-bit.** The summation order changed;
  that is the content of the change, not a regression. Downstream: `limiter::TruePeakLimiter` and
  `poweramp::PowerAmpStage`, hence OrbitCab and orbit-amp. Three storage budgets moved with it (law 11d):
  `TruePeakMeter::storageFor` 296 → 392 B at 48 kHz stereo, and the solve budgets that carry it.

### `wasm-audio` — the tier gets the SIMD ISA it was already tested against

The `wasm-audio` preset now compiles with `-msimd128`, so `core::firDot` reaches its `wasm_f32x4_*`
kernel instead of falling back to scalar. This is an ISA switch, not a numeric one: the summation order
is stated in the source, `-ffast-math` stays off, and the tier's suite was measured to move not one bit
with the flag on. Still NOT `-mrelaxed-simd` — `relaxed_madd` is implementation-defined and breaks
determinism between machines rather than merely between tiers.

### `core` · `eq` · `dynamics` · `saturation` · `limiter` · `oversampling` · `mastering` · `tools` — the chain says what building it costs, and re-preparing it costs nothing

Law 11d's remaining half (`docs/DSP-ARCHITECTURE.md`): its budgets covered the solver's calls and said in as
many words that `create` and `configure` — the chain's own storage, and much the larger number — were **not
budgeted yet**. They are now, and closing that turned up two things worth more than the budgets.

- **A refused `fc_master_create` used to ask the heap for as much as 1 668 312 bytes on its way to saying no.** It
  built the instance, prepared the renderer, and let the chain allocate down to the first stage that refused —
  on the DEFAULT geometry a 20 Hz rate cost 392 408 B and a 300 ms compressor lookahead 394 456 B, even a
  refusal on the chain's own front door cost 51 288 B, and the same lookahead refusal at 48 kHz with 16
  channels and an 8192-sample quantum cost **1 668 312 B in 9 allocations** — the numbers scale with the
  geometry, so the small ones are examples and not a maximum. On the wasm tier an allocation that cannot be served is not a refusal at all but
  the end of the module, which is the whole subject of law 11d. **`MasteringChain::admits (fs, nch, cfg)`**
  now reaches the entire verdict — every stage's own gate included — **without a single allocation**, and
  `create` calls it before it builds anything. A refused create now allocates nothing, pinned over ten
  refusals.
- **`fc_master_configure` used to ask for 345 224 bytes to change nothing.** It built a SECOND 331 KiB
  `eq::EqEngine` before releasing the first — a transient peak larger than everything else the call asks for
  put together — plus three temporaries from copying `assign`s. The engine is now **re-used** and the
  temporaries are gone: a re-preparation at the handle's own geometry asks the heap for **0 bytes**, measured
  on every row of the matrix. The re-use is bit-identical, proven by null rather than by argument: a chain
  prepared three times with parameters written in between renders the same programme sample for sample as one
  prepared once, on three topologies — one of them the EQ alone, where no compressor, limiter or dither
  stands between the engine and the comparison — and again across a chain MOVED to another rate and quantum.
- **Nine headers publish what their `prepare()` allocates, through the function `prepare()` itself sizes by** —
  `core::DelayLine`, `core::DryAligner`, `oversampling::PolyphaseOversampler`,
  `limiter::TruePeakLimiter` (and its sliding window), `dynamics::Compressor`, `saturation::Saturator`,
  `eq::EqEngine` (object and scratch), `mastering::OfflineRenderer`, `mastering::MasteringChain`. Each
  `storageFor(...)` answers FALSE on exactly the arguments its `prepare()` refuses, so the budget and the
  preparation cannot drift; `MasteringChain::prepareBytes` and `mastering::createBytes` aggregate them.
- **`TruePeakLimiter`, `Saturator` and `Compressor` publish `latencyFor(...)`**, and their own `prepare()`
  runs through it — so a composite sizing a dry aligner, and a budget sizing the same aligner, read the
  number the prepared stage will report rather than deriving it a second time.
- **Temporaries removed from three copying `assign`s** (`orbitcab` and `orbit-amp` share these modules):
  `osBuf.assign (n, std::vector<float> (K·F))` built one buffer to copy `n` times — 524 288 B of peak on the
  biggest topology the mastering chain builds, and 64 MiB asked of `TruePeakLimiter` directly at its own
  block cap — and `assign (n, DelayLine {})` cost a temporary, `n` copies and then `n` reallocations.
  `core::prepareDelayBank` builds each line at its final size instead. Same audio, bit for bit; a fresh
  `create` is 4 120 B and 7 allocations lighter on the default topology, and a re-preparation asks for nothing.
- **C ABI (version unchanged at 1 — no struct moves):** `fc_master_need_create (cfg, fc_need*)` and
  `FC_NEED_CONFIGURE` for `fc_master_need`. `need_create` is a **dry run**: every refusal the create can reach
  before its first allocation comes back with the same status, so `FC_OK` always carries a non-zero budget and
  a caller can validate a configuration without paying for the attempt. `fc_master_need` keeps answering a
  NUMBER rather than a permission — a configure is budgeted mid-stream, where the call itself is
  `FC_ERR_STATE` — because the cost of a call does not depend on the moment it is made.
- **Three repairs the review round found in the diff itself**, each with the input that found it:
  `core::DelayLine::prepare` now re-clamps its tap into the new capacity (`prepare(8); setDelay(8);
  prepare(2)` left a tap outside the ring and `process()` then read before its start — an invariant this
  call has broken since it was written, which the delay bank made reachable in one more shape), and
  `core::prepareDelayBank` clears the tap so a re-used bank ends where the `assign` it replaces ended;
  `MasteringChain::reprepareBytes` asks the chain's OWN containers rather than trusting the geometry it
  remembers (a MOVED-FROM chain keeps `prepared_` and its scalars and has given its buffers away — the
  budget answered 0 for a re-preparation that really asked for 512 B); and the limiter reserves before it
  resizes its per-channel scratch bank, which removes a growth step the chain was introducing itself
  (a `resize` past the capacity grows geometrically, so widening 3 buffers to 4 asked for 6 — 48 B over).
  What a container does on its OWN growth stays the caller's margin, exactly as law 11d says: the `win` row
  measured MSVC's `assign` asking for 72 floats where 68 were wanted, which is why the bound on a GROWING
  re-preparation is a bound and only the FRESH budget is exact.
- **A re-preparation now KEEPS storage the old form gave back — on ONE edge, and being exact about which one
  took a correction.** The diff pass found the retention; the fix round found that the obvious explanation for
  it was wrong. `assign (n, DelayLine {})` did NOT free a narrower bank's rings in general — it copy-assigns
  into the lines that survive, and a `vector` copy-assigned from a shorter one keeps its capacity, so the old
  form retained too. The two diverge when the bank GROWS PAST ITS CAPACITY: `assign` reallocated and built
  fresh lines from the temporary, destroying the old ones; `reserve` + `emplace_back` MOVES them and their
  rings travel along. Measured on that edge: `Compressor::prepare (48000, 64, 1, 50 ms)` then
  `(48000, 64, 2, 1 ms)` holds **9 880 B** where it held 472; `TruePeakLimiter::prepare (48000, 65536, 1)`
  then `(48000, 256, 2)` holds **1 133 236 B** where it held 88 756. Separately, re-using the EQ engine means
  a chain moved from an 8192-sample quantum to a 256-sample one asks for **0** instead of 341 120 B and holds
  63 488 B more. This is what law 11d's budgets are stated over ("one already prepared keeps storage that
  still fits"); a consumer that must give a large geometry back destroys the stage rather than re-preparing it.
- **A REFUSED re-preparation no longer leaves its own arguments visible through the ungated readouts.** Moving
  every check ahead of the first write changed three answers, all of them on calls that return `false`, and
  all of them toward what law 11(b) asks for — a refused call touches nothing. Stated rather than left to be
  met, because these headers are shared with the plug-ins: `eq::EqEngine::sampleRate()` reported the rate of
  a preparation refused on its WIDTH (48000 → a refused 96000 answered 96000, now 48000);
  `saturation::Saturator::latencySamples()` reported the round trip of a topology the oversampler had just
  refused (63, now 0); and `limiter::TruePeakLimiter::effectiveReleaseMs()` answered in terms of a rate the
  refused call brought (0.166667 ms, now 0.083333) — that one only, since the ceiling is a clamp on the
  parameter and never sees the rate. Found by the diff pass and the code-review round, not by a test — no
  suite reads those three after a refusal.
- **Two products that were undefined are now merely large:** `saturation::Saturator`'s oversampled scratch
  and `eq::EqEngine`'s sidechain scratch were sized by an `int` product (`maxBlock * os`, `maxBlock * ch`)
  and are now computed in `size_t`. Past `INT_MAX` the old form was signed overflow — in practice a
  wrapped, far too small buffer — and the new one is an honest request the heap will refuse. Only a direct
  consumer can reach it: the mastering chain caps its quantum at 8192.
- **The allocation counters THIS BUDGET IS PINNED WITH install every form of `operator new`,** the
  over-aligned one included — `MasteringChainTests` and `MasterAbiTests`, the two that carry the numbers
  above. Without it `eq::EqEngine`'s 331 KiB — the largest single request a create makes on the default
  geometry; at 16 channels and an 8192-sample quantum the saturator's flat scratch is 8 MiB — is invisible,
  and a budget check would have compared two numbers that both left it out (the blindness P52 names). Pinned
  over 4 rates × 3 widths × 4 topologies, byte for byte, for `create` and for `configure`.
  ⚠️ **This sentence used to say "the suites' allocation counters", and that was false of the suites at
  large.** Counted on the tree: **11 test TUs of 61** install the over-aligned form. The other 50 — among
  them `LimiterTests`, `NamStageTests`, `LinearPhaseEqTests` and `EqEngineTests` — assert "no heap
  allocation" with a counter that CANNOT SEE an over-aligned request, so what they prove is narrower than
  what they say: no call to two of the eight replaceable forms, not no call to any of them.
  **What is NOT claimed here: that anything was actually slipping past.** Measured while closing P52 —
  three stands, including a backtrace over all 128 binaries — no over-aligned allocation reaches any
  `process()`, and 45 suites make none of their own at all. `EqEngineTests` is one of them: it builds
  `EqEngine eng;` on the STACK, so the 331 KiB above is a cost a `make_unique` consumer pays, not one that
  slipped through this suite's blind spot. The hole was real and latent, which is the honest description
  and the reason it survived. Closing it is **P52**.

### `analysis` · `mastering` · `tools` — a meter's store is counted in samples, a call publishes what it will allocate, and an instance that aborted refuses

- **`analysis::LoudnessMeter` sizes its gating-block store through ONE function, in SAMPLES:**
  `storageFor (sampleRate, maxSamples, Storage&)`, which the new `prepareForSamples()` sizes itself with;
  `prepare (seconds)` is now a thin wrapper and still reads a NaN or negative duration as 0 s. **A capacity
  that cannot be represented is REFUSED** — +inf samples, a block count past the `int` index, a rate whose hop
  overflows an `int` — where it used to reach an out-of-range float-to-`size_t` conversion. On every
  representable input `prepare (seconds)` sizes exactly the store it always did (pinned byte for byte by a
  hand-derived table in `LoudnessConformanceTests`).
- **`mastering::TargetLoudnessSolver` sizes both of its meters in samples** — `frames + ceil(fs)` — and no longer
  through `frames / fs + 1` seconds, which is +inf at a finite rate the chain accepts (without EQ and limiter),
  where the same solve used to keep a different number of gating blocks on different platforms. Pinned:
  `LoudnessSolverTests` keeps all 197 blocks of that programme, 0 dropped. At ordinary rates the store can now
  differ from before by ONE block, where the trip through seconds rounded across an integer; those blocks are
  margin, not need — the store still holds every block a programme produces — and no measurement changes.
- **Budgets — what a call will ask the heap for, from the very functions it sizes itself with:**
  `TargetLoudnessSolver::prepareBytes / solveBytes / measureRangeBytes`, `TruePeakMeter::storageFor`,
  `QuantileHistogram::storageBytes`, `LoudnessMeter::Storage::bytes()`. REQUESTED bytes — allocator headers,
  alignment and fragmentation are the caller's margin, and none of this promises that a heap can serve them.
  Each is exact for a FRESH object — a prepared one keeps storage that still fits — and 0 (`storageFor`:
  false) wherever the matching call refuses, a channel count included.
- **`dynamics::offline::QuantileHistogram::prepare`, refusing, writes nothing past its disarm** (law 11(b)), now
  that it sizes itself through `binsFor`. A refusal drops the bins, as it always did; one by the 4e6-bin ceiling
  used to go on and leave the refused range and width behind — `binWidth()` read them, and `add()` split
  below/above range by them — and now the last successful preparation's stay, as after every other refusal.
- **`tools` (C ABI) — `FC_ERR_POISONED` (14).** An exhausted heap aborts a wasm module inside the core, and the
  module is not stopped by it: the next call used to be answered by objects the abort had left half-changed.
  Every status-returning entry point now refuses, touching nothing, once any call has failed to return (an
  abort, a trap, natively an escaped exception); the page discards the instance. The `*_default` writers and
  the identity queries stay callable. **The module is not re-entrant, and now says so:** a status call made
  while another is still running — from a native `new_handler`, a signal handler — cannot be told apart from
  the first call after an abandoned one and is answered 14 too, so a native host that used to call back in from
  a `new_handler` now poisons the module.
- **`tools` (C ABI) — `fc_master_need (h, op, frames, fc_need*)`**: the budgets of `solve` and `measure_lra`,
  forwarded field by field and never summed (`callBytes`, `solverPrepareBytes`, `facadeBytes`,
  `solverPrepared`) — `facadeBytes`, the facade's own `sizeof` of the solution record, is the one number this
  ABI hands back that the core did not compute. A `frames` past INT_MAX is `FC_ERR_RANGE`, then an unknown `op`
  `FC_ERR_ENUM`. `create` and `configure` are not budgeted yet. Both additions are additive: no struct
  moved, `FC_MASTER_ABI_VERSION` stays 1, and the header now states that rule for new codes.
- **`docs` — law 11d: memory that cannot be had is not a refusal.** Exhaustion is fatal on every row — never a
  `false` — and the explicit exception to 11b. In its place the core publishes a DEMAND (a bound on what an
  allocating call holds at once, from the functions its `prepare()` sizes itself with) and the C ABI poisons a
  module whose call never returned. The chain's own storage (`create`, `configure`) is not budgeted yet.
- **`tools` — `fcore::Probe::prepare` answers with its meter.** It ignored the meter's return value, which
  could only fail on a channel count the probe had already checked; the meter now also refuses a store it
  cannot represent (3e8 s), and a probe that ignored that would report prepared and measure nothing — where it
  used to ask the heap for the impossible store.

## v0.30.0 — a pause is silence, a refusal has a name, and the chain answers through a C ABI (`core`, `nam`, `rigplayer`, `mastering`, `tools`)

### `nam` · `rigplayer` — a lane that stops being fed is DRAINED, not frozen (law 11a)

- **`nam::NamStage` HAS a falling edge now**, and the entry that stood here — "no falling edge …
  recorded, not silently claimed" — undersold it twice. A lane the host stops handing over used to be
  skipped whole, so BOTH its network window and its two `core::StreamResampler`s froze and were
  replayed on the return; and a slot the blend law puts to SLEEP was not handed to the stage at all,
  which leaks the same way for a different reason. Measured through `rigplayer::RigPlayer`
  against the commit this fix was made on, worst |out| out of DIGITAL SILENCE **over every leaving
  phase of the probe tone** / tail in host samples — a single leaving point is a lower bound and not a
  size, because what comes back is whatever the frozen state was holding: a memoryless capture
  **0.525665 / 125 at 44.1 kHz** · 0.528352 / 183 at 88.2 · **0.524339 / 193 at 96** · 0.521520 / 300
  at 176.4 · 0.522065 / 319 at 192 — i.e. **≈0.52 wherever a rate-matcher is installed at all**, and
  the spread an earlier single-phase table showed was the probe's own phase rather than the defect's
  shape; a 2001-tap capture **0.499533 / 2003 at 48 kHz**, where no rate-matcher exists at all; and a
  sleeping slot **0.500000** for a whole receptive field. ⚠️ **Those TAILS are this release's own
  doing.** The frozen state has a closed form — a stopped pair holds `kTaps` host samples in the down leg
  and `kTaps` model samples in the up leg — so the tail is **`kTaps·(1 + hostSR/modelRunSR)` host
  samples, which is exactly TWICE the reported latency** (2 × 61.4 = 122.8 against a measured 125;
  2 × 96 = 192 against 193; 2 × 160 = 320 against 319, the slack being the probe's block floor). The
  64-tap sinc announced below therefore made the tail **sixteen times longer by construction** (64 taps
  against the cubic's 4); on `v0.29.0` the same probe reads 10 · 13 · 12 · 20 · 20 host samples, whose
  measured ratios of 12.5×…16.1× fall short of the 16 only because a few-sample floor is a larger share
  of 10 than of 125. ⚠️ **The AMPLITUDE, though, did not wait for this release:** a reader upgrading from
  `v0.29.0` already had **≈0.52 out of digital silence over 10–20 samples** from the frozen rate-matcher
  (0.518397 at 44.1 kHz, 0.533270 at 88.2, 0.515475 at 96), on top of the delay line's 0.249992 in the
  first 3 samples — which is all there was at 48 kHz, where no rate-matcher exists. What this release
  lengthened is the TAIL, not the leak. An absent
  lane is now fed the digital silence it is actually receiving, for as long as its state can still be
  heard, and a sleeping slot gets a width-zero call — **exact zero** at 8 · 22.05 · 44.1 · 48 · 88.2 ·
  96 · 176.4 · 192 kHz, on a memoryless capture and on a real 6332-sample WaveNet, at gap widths 0 and
  1 and on both slots.
  - **The ceiling in `RigPlayer::process` was derived for the delay line and spent on both halves.**
    It is 3.8x larger for the models: both slots freeze at once, so their weights SUM rather than pick,
    and a frozen rate-matcher is not a replay — its phase rows are normalised by their SUM, so their
    MODULUS exceeds one. The return path is linear for a Linear capture, so the ceiling
    `A · max_n ‖h_n‖₁` is ATTAINABLE and was attained: **0.949383 (-0.45 dBFS) at 44.1 kHz, 100.00 %**
    by the sign pattern of the worst row, where a 220 Hz sine swept over every leaving phase reaches
    55 %. The published 0.25 / -12.0 dBFS / "in the first three samples" were the delay line's.
  - **BREAKING (a number, not an API): `NamStage::prewarmSamples()` answers for a Linear capture.**
    It reported **0** for an impulse response of any length — the config declares the field as a plain
    number and nothing read it, while NAM's own answer for that architecture is zero. It now reports
    the MEMORY, which is one less than the declared taps: a one-tap capture is a gain and still
    answers 0. A consumer that fades a model in by this number (`RigPlayer::warmFor`) waits a real
    interval for a real IR now, where it used to wait none.
  - The field is no longer capped at `1<<20`: the number is SPENT now, and a cap on a spent number is
    a silent under-drain. A `Linear` capture is also charged its partitioned-FFT ring (NAM runs one
    past 256 taps by default, and it holds input spectra past the field — measured, 1.909e-08 on 46
    samples of the return with a dense kernel), and a recurrent architecture is floored at half a
    second of the RUN rate rather than trusting `GetPrewarmSamples()`, which is half a second of the
    model's TAG and answers **1** when there is no tag.
  - **The dilated stack's LEGACY spelling is read.** NAM takes either `kernel_sizes` (an array) or a
    single `kernel_size` for every layer; only the array was read, so a legacy capture reported a
    field of zero — and on a `SlimmableWavenet`, whose own answer is also zero, nothing knew it at
    all: measured **0.462117** out of digital silence on a loaded model. A CONTAINER is now asked
    through for the FFT ring and the recurrent floor as well, not only for the field (measured on
    loaded models: 1.48e-08, and 0.499275 against 0.419115).
  - **`nam::NamStage::drainedSamples()`** is new: how many samples of silence the stage has fed to
    lanes the caller stopped handing over. It exists to be tested rather than acted on — past the debt
    an absent lane's output is zero whether it is still clocked or not, so "it drains, and then it
    STOPS" has no witness in the audio, and three mutations of the drain's LENGTH survived a suite of
    960 checks before it existed.
  - **What it COSTS, measured on a real 6332-sample capture:** a departing lane goes on costing the
    stereo price for its drain and NOT more — the drain runs one lane-block per block, exactly what that
    lane cost while it played, so there is no spike: **12.9 ms of CPU per stereo sleep at 48 kHz and 9.3
    at 44.1**, and **0.0000 ms per block** once the debt is spent. A host that alternates widths every
    block re-arms the debt each time and therefore pays for the absent lane indefinitely — that is what
    "a pause is silence" means when the pause is one block long, not a defect.
  - **An RT fix that came with this and outlives it:** NAM grows a `Buffer` capture's window on demand
    INSIDE `process()`, and `Reset` pre-grows it only through a prewarm that is zero samples long for
    a Linear capture — so the first audio call after a prepare allocated. It was dormant while
    instance 1 was never touched on a mono host; the drain touches it. Both instances are now walked
    once on the message thread, and the allocation counter starts at the FIRST call. It costs a real
    capture NOTHING — a `SlimmableContainer` of plain WaveNets answers a prewarm of its own, so the walk
    is skipped: 9.7 / 10.7 / 12.6 / 14.4 ms to prepare one, with the walk and without, identical. Where
    it does run (a Linear capture, a true `SlimmableWavenet`) it is +1.1 ms on an 8193-tap FFT capture at
    maxBlock 8192 and ~0 at 512.
  - Two things this does NOT close, said plainly: an LSTM's cell has no flush length, so its drain is
    a bound on NAM's own heuristic and not on the memory (0.419 against 0.023 for a lane clocked
    through the whole gap); and `prepare()`/`reset()` still do not clear a network's window, so a lane
    that is PRESENT can be handed silence and reply with 0.2246 — a stream-restart question, recorded
    with its number rather than claimed.


- **`nam`: THE RATE CONTRACT IS A FIXED WINDOW, AND THE CONSUMER MITIGATIONS IT FORCED ARE GONE
  (`nam::NamStage`, `rigplayer::RigPlayer`).** `install()` used to admit a model whose tag was within
  half a hertz of the rate the stage was RUNNING, and `prepare()` then adopted the accepted tag — a
  fuzz on equality spelled as a moving reference, so every accepted load moved the goalposts for the
  next one. Measured on the base commit, half-hertz steps at a 48 kHz host: `{load}` 1 step,
  `{load, process}` 1, `{load, prepare}` 66 (stopped by the retire queue, not by rates),
  `{load, process, prepare}` **5000 with no refusal at all**, run rate walked to 45500.0. The walk was
  clocked by `prepare()`; audio only drained the retire queue, which corrects the two numbers the
  header used to carry. And its cost was the opposite of what it looked like: a stage walked to 47900
  **refuses an ordinary 48000 capture** — the window moved, it never widened.
  - **The reference is the constant now.** `NamStage::acceptsModelRate(modelSR)` is a pure predicate on
    the model's own tag — "this stage takes a model exactly when a factory-rate host would not resample
    it" — so the accepted window is `[47999.5, 48000.5]` for the life of the process, whatever has been
    loaded before. It reads no stage state, so it is settled in `prepareModel()` and a doomed model
    never pays for its prewarm. New public surface: `kModelRateTolerance`, `acceptsModelRate`,
    `maxLatencySamples(hostSR)`.
  - **`install()`'s reconfiguration test is EXACT.** Two half-hertz tolerances of the same size do not
    compose: a backend prepared for host 48000.4 and installed into a stage at 48000.6 kept a
    rate-match computed for the wrong host and reported **0 samples of latency where the policy charges
    64**; the mirror reported **64 where the policy charges 0**. Both now equal the policy.
  - **The model scratch is sized for the path that RUNS**, from the real host rate. It was one ratio
    serving two paths and wrong at both ends: with no resampler the model is clocked by the HOST, so a
    whole `maxBlock` chunk is staged in a buffer the ratio had sized just UNDER `maxBlock` (reachable
    at `maxBlock >= 34*hostSR`, i.e. 34 seconds in one call — which law 11(a) explicitly invites);
    and `max(8000, hostSR)` substituted an assumed host rate, losing frames in silence below 8 kHz —
    **-1.22 dB at a 6 kHz host, -3.00 at 4 kHz, -6.05 at 2 kHz, -9.09 at 1 kHz, 0.00 at 8 kHz exactly**.
  - **Behaviour change where a conversion cannot be sized** — a host rate that is not positive, or so
    slow that one block exceeds the arithmetic (about 0.023 Hz at a 512-sample block), or a `maxBlock`
    past `(INT_MAX-16)/2`, which is doubled one line later and used to overflow there: the load is now
    REFUSED, visibly, where it used to succeed and then produce garbage or convert an out-of-range
    double to an `int`. Nothing shipped reaches them; `RigPlayer` maps every host outside `(0, 3e6]` to
    the factory rate first.
  - **What does NOT move, and it is stated because it nearly did:** `maxModelFrames` is also the block
    handed to NAM's `Reset`, and NAM prewarms in WHOLE blocks — so it decides how many samples of
    silence a stateful capture is warmed with, and therefore its state at the first real sample. The
    direct branch keeps the same `+16` the old expression had, preserving equal model/host rates
    (measured at blocks 64/256/512/1024/4096). Fractional hosts can change even an integer-tagged
    model: host 47999.75, tag 48000, block 512 changes Reset from 529 to 528; the decaying-cell LSTM's
    first output changes from 0.1597609967 to 0.1600718498. Larger blocks can differ by more than one.
  - **Exact install reconfiguration also adds stateful prewarming where latency does not change.**
    A handle prepared at 48000.1 and installed at 48000.2 (block 512) is now Reset again. NAM's LSTM
    Reset advances the existing cell rather than restoring it: first output 0.0548291542 instead of
    0.1600718498 for the decaying-cell fixture, with zero latency in both cases. Matching split loads
    and fused loads still prewarm once. The behaviour is recorded by tests rather than redesigned here:
    separating the scratch capacity from the Reset/prewarm schedule is its own task (plan P39a).
  - **`RigPlayer::dryAlignerCapacity` lost its `max(256, …)` floor and its `+2`** — both existed only to
    mitigate the walk, and their price was measured (first silent clamp after 82256 half-hertz steps at
    a 48 kHz host, 68511 at 96 k, 41021 at 192 k, 560 at 384 k, 72 at 3 MHz). It is now
    `NamStage::maxLatencySamples(usableSampleRate(fs)) + 1`: the bound asks the accepted window's LOW
    edge, because a lower model rate is a longer round trip and the nominal rate reads up to 0.0208
    samples short — enough to land on the wrong side of a rounding boundary at 2999249 Hz. The `+1` has
    exactly one reason, `DryAligner`'s usable range being `capacity-1`.

- **`tools`: A C-ABI FACADE OVER `mastering`, AND A CLI THAT PROVES IT ADDS NOTHING (`fc_master_*`,
  `fcore_master`).** `felitronics::mastering` is now callable from a browser worker; the desktop
  application links the same module AS C++ past this surface entirely, and that second half is the
  constraint that shapes the first. What the facade duplicates is exactly two things — the enum codes
  and the field mapping — and no arithmetic at all.
  - **`fcore_master selftest` is the acceptance, not a smoke test:** the same programme rendered
    through the C entry points and through a direct C++ call, in ONE binary on ONE machine, compared
    bit for bit — **0 of 384 000 samples differ**, with a precondition asserting the chain actually did
    something (measured on the shipped fixture, and printed by the test rather than pinned here). Block-independence
    survives the boundary at call sizes **1, 337, 4096 and whole-file, 0 differing samples each**,
    because nothing here re-blocks anything: the call reaches `MasteringChain::process()` in one piece
    and the chain's own fixed quantum stays the only clock. It runs in the `wasm-audio` tier too.
  - **Four disciplines the comparison needs, each measured rather than reasoned out.** The target does
    NOT inherit `tools/`' `-ffp-contract=off` — the two flavours differ in **203 269 of 288 000
    samples**, and two TUs of one binary built with different flags measured as AGREEING because the
    linker merged the header-only instantiations, which is agreement by link order. `setParams` then
    `prepare` is not the same render as `prepare` then `setParams` (**59 259 of 80 000 samples, 0.0715
    full scale**, all of it `stereo::MonoBass`, whose `reset()` snaps the width where its setter ramps
    it over 20 ms). A zeroed parameter struct is not `MasteringChainParams{}` (**287 998 of 288 000**),
    so `fc_master_params_default()` exists. And a second programme through one handle without
    `fc_master_reset` is a different render — the selftest prints the count for the fixture it
    ran, rather than carrying a figure that a later fixture change would quietly falsify.
  - **`configure` re-prepares, and is refused once audio has been handed over.** Reading `resolved()`
    straight back from a deferred `setParams` reports the PREVIOUS parameter set — **5.0000 dB** on the
    limiter ceiling, **150.0 ms** on its release — those are the ERRORS, not values the chain reports:
    a first set at −1 dBTP / 50 ms read back after a second asked for −6 / 200 — and the core's own
    defaults on a fresh chain.
    Applying early instead is worse: `Dither::setParams` reseeds on a seed change and
    `EqBand::setParams` snaps while uninitialised and glides after, so N configure calls with no audio
    between them would stop equalling one call with the last set — a render that depended on how many
    times a knob moved before the button was pressed.
  - **Handles are an index and a generation, not pointers.** Under emscripten, destroy-then-create
    returned the same address **19 times out of 19** with emmalloc and with dlmalloc, against **0 of
    20** for native malloc: a pointer handle is a use-after-free the developer's machine never
    reproduces and the shipping tier reproduces always, in a linear memory with no page to trap on.
  - **Ranges are the core's business and the facade does not check them** — the core clamps by design
    and reports what it clamped to, and a facade with its own range table is a second copy of every
    stage's limits. What it does refuse is a non-finite parameter (where the core has no verdict: a NaN
    gain becomes 0 dB silently), a bad memory span, an unknown enum code, and a struct whose version or
    size this build does not know.

- **`mastering::TargetLoudnessSolver::measureInputLoudnessRange()` now refuses a POISONED programme.**
  It checked `gatingBlockCount` and `droppedBlocks` and not `nonFiniteSubHops`, so it returned SUCCESS
  on a programme its own meter had already flagged — and the meter is a local, so the caller could not
  check for itself. Measured on 30 s alternating 3 s loud / 3 s quiet with every loud second poisoned:
  **4.8000 LU clean against 21.4000 LU poisoned, both `true`.** That number is the far end of the
  `maxLraLossLu` DELTA, so the constraint was being judged against a range the programme does not have.
  NB the counter's granularity is a COMPLETED 10 ms sub-hop, so a non-finite sample inside the final
  partial sub-hop is still not refused; closing that would need a second scan of the audio.

- **`mastering::MasteringChain` COUNTS the samples its input gate substitutes
  (`nonFiniteInputSamples()`).** The count is what turns a plausible render of a programme nobody
  submitted into a visible one. Counting rather than refusing is deliberate and measured: a NaN through
  the chain is bit-identical to a sanitised sample, while refusing the call would throw away every good
  sample travelling with it — and the size of that loss would depend on the caller's block size, which
  is the one thing the internal quantum exists to make irrelevant.
  - **THE SHAPE OF THE COUNTER IS NOT FREE, and the first version of it was not.** Writing it as
    `if (! isfinite(v)) { ...; ++member; } else ...` reads as the same code and stops the compiler
    vectorising a pass that runs over every input sample of every quantum. Measured, arm64 Release,
    best of seven interleaved runs over 20 000 quanta of the gate loop: **1.357 ms** branchless without
    a counter, **3.536 ms** with the branch (**x2.61**), **1.369 ms** for a branchless select plus a
    LOCAL counter (**x1.01**) — which is what ships. Whole chain, same source against `main`'s headers
    and against this branch's with every DSP stage off: **0.531 ms against 0.407 ms**, identical
    checksums. The shipped form keeps the original expression `std::clamp(isfinite(v) ? v : 0.0f, ...)`
    verbatim, so the arithmetic is unarguably unchanged.

- **`mastering`: TARGET-LOUDNESS SOLVER, NAMED CONSTRAINTS AND THE STATISTICS BEHIND THEM
  (`mastering::TargetLoudnessSolver`).** Hits a target integrated loudness under a stated true-peak
  ceiling in a bounded number of renders, and refuses BY NAME instead of crushing the programme when
  the target cannot be had.
  - **It is one scalar search, not two loops.** `preLimiterGainDb` (g) and `limiter.ceilingDbTp` (c)
    look like two knobs and are not: the limiter's reduction is `min(0, c - smaxDb)` and `smaxDb` is
    taken after the gain node, so everything the limiter does depends on `d = g - c` alone and `c` is a
    pure output scale — `y(g,c) = 10^(c/20) y(d,0)`, measured over a 3x3 grid at
    `9.1e-07 .. 1.7e-06` on a programme peaking at 0.89. So the ceiling programmed into the limiter is
    an OUTPUT of the solve, and the P11 inter-sample derate costs whatever the material's actually is
    (**+0.0005 to +0.1028 dB measured in situ** on five real mixes) instead of a flat 1.2 dB off
    every track.
  - **New: `limiter::TruePeakLimiterTap`** — the limiter's gain-reduction trace and the reconstructed
    peak it saw, on its OWN `F*fs` grid, plus `maxReconstructedPeakDb()`. Folding the trace to baseband
    needs a rule and the rule belongs to whoever reads the statistic: a min-fold biases the mean by
    **+0.0029 to +0.0080 dB** and the active fraction by up to **+0.0009** across the whole release
    range, while `max` and the upper quantiles are untouched. A tap too short REFUSES the whole call.
  - **New: `mastering::MasteringChainTaps`** — the chain forwards the limiter's traces, the
    compressor's own tap and the signal at the pre-limiter node, with each stage's offset from the
    chain's INPUT stated rather than implied, so a statistic is cropped to the window that carries
    programme. `OfflineRenderer::render` grows one templated tap sink; there is still exactly one copy
    of the `out[n] = y[n + D]` arithmetic.
  - **BREAKING (behaviour), `mastering::MasteringChain`: a parameter set written BEFORE `prepare()` is
    now KEPT.** `prepare()` ended with `pendingParams_ = params_`, which replaced the caller's pending
    write with the last APPLIED set — defaults, on a fresh object. Measured:
    `setParams(inputGainDb = 12)` then `prepare()` then `process()` delivered the input **unchanged**,
    12 dB that simply did not happen, with no refusal and no way to find out. Every stage already
    honours that order (`Compressor`, `TruePeakLimiter` and `Dither` all re-apply their stored
    parameters inside `prepare()`); the composite was the only place that did not. The parameters are
    now applied inside `prepare()`, so `params()` and `resolved()` describe the prepared chain rather
    than the previous one — mid-stream they still lag a `setParams()` by up to one internal quantum,
    which is the documented design.
  - **A limit already broken at the least drive the search will use is UPSTREAM, and says so.** The
    compressor's gain reduction cannot move at all (its node is before the gain), and the loudness
    range, the peak-to-loudness ratio and the limiter's own gain reduction all get WORSE with drive —
    so naming the loudness target as the reason points the user at the wrong number. Measured on the
    corpus: with a 0.4 LU range allowance, an ordinary compressor setting spends **0.80 to 3.70 LU**
    before the solver applies a single dB.
  - **Digital silence gets `MeasurementInvalid`, not a plausible answer.** `analysis::LoudnessMeter`
    returns the literal -120.0 when nothing passes its absolute gate, and a ceiling derived from a
    -200 dBTP peak reading clamps to +60 dBTP — i.e. it would switch the limiter OFF and report
    success. A programme that is merely too QUIET to measure at the starting gain is a different case
    and is bootstrapped from the peak instead: a flat tone at -71.7 LUFS is unmeasurable at 0 dB and
    ordinary 56 dB up, and refusing it would be a verdict about the starting gain.
  - **The answer does not depend on where the search started, and it used to.** A step that moved the
    gain and the ceiling together is exact — it is the scale law above — but exact at a FROZEN DRIVE.
    Measured on one programme and one request (-14 LUFS, -1 dBTP): from a 0 dB start the answer was
    8.5 dB of drive with no limiting, PLR 11.9 and LRA 4.1; from a 55 dB start it was **47.6 dB of
    drive, 38.05 dB of limiter gain reduction, PLR 4.6 and LRA 0.10 — also reported `Solved`.** Adding
    `minPlrDb = 8` then made the second one `TargetUnreachable` while the first stayed Solved, so the
    VERDICT depended on the start too. The ceiling now tracks its aim on every step, in both
    directions, and a warm start unwinds instead of freezing. Pinned over starts of 0 to 55 dB.
  - **`in == out` is refused.** One render in place is well defined and `OfflineRenderer` still
    supports it; a SEARCH is not, because every pass after the first reads the previous pass's master.
    Measured: a solve reported -22.996 LUFS and the gain it returned, applied to the untouched source,
    gives -29.000 — the answer missed its own programme by 6.0 LU.
  - **The request carries no delivery policy and no hidden state.** `targetLufs` and
    `maxTruePeakDbTp` have no defaults (NaN, refused) — "-14 LUFS, -1 dBTP" is a product's decision,
    not a core's. The input's loudness range moved out of the solver and into the request for the same
    reason: held as solver state it outlived the programme it described, and track B was judged against
    track A's range.
  - **A target past the +-60 dB gain node is `TargetUnreachable` with `GainRange` named**, not a
    pass limit: the movement test compared the UNCLAMPED step, so a saturated actuator re-rendered the
    same point until the budget ran out — measured, 29 identical renders of a 33-render budget.
  - Also fixed in the same pass, each found by a review round and reproduced before being acted on:
    the tap capacity arithmetic overflowed in `int` before its cast to `long long`; the tap counters
    advanced even when no tap was requested; `OfflineRenderer` checked the tap capacity per BLOCK, so a
    short tap failed half way through a render with output already written; the limiter's statistics
    window ignored the limiter's own interpolator latency, which cost the whole reaction of a peak in
    the last 32 samples; an infinity of the wrong sign disabled a constraint the caller meant to be
    unsatisfiable; and the solver did not check that its sample rate was the chain's.

- **BREAKING (behaviour), `core`, `nam`: `StreamResampler`'s interpolation kernel is now a 64-tap
  polyphase windowed sinc, not a Catmull-Rom cubic. Every model at 44.1 kHz sounds different — brighter
  in the top octave, and without a bass artefact it should never have had.** (Read "brighter" as scoped:
  it holds for hosts up to 96 kHz. At 176.4 and 192 kHz the fixed 64-tap window under a ~4:1 decimation
  costs up to 0.8 dB at 20 kHz — still far better than a kernel with no anti-aliasing at all, but not
  flat; the rows are pinned in both directions and §6.7 of the doc says so.) The cubic had no
  anti-aliasing of any kind, and P32 measured what that cost; this is the fix, and it is not neutral.
  - **Passband.** One round trip 44.1 ↔ 48 kHz, coherent carrier / worst phase, at 17.64 kHz:
    **−4.17 / −9.27 dB → +0.0002 / +0.0000 dB**. At 20 kHz **−5.48 / −14.79 → −0.0133 / −0.0135**. The
    carrier droop and the phase modulation were the same mechanism in two coordinates and both are gone:
    modulation depth is at most **0.0004 dB** where it reached 8.7, inside the Kaiser window's own
    derived passband ripple.
  - **The decimating leg had no stopband at all** — a flat −3 dB rms and a **0.0 dB sample peak** above
    the output Nyquist, because at phase t = 0 the cubic's weights were (0,1,0,0), a bare sample pick.
    Now **−9.08 / −15.41 / −27.45 / −47.56 / −88.77 dB** at 22.1 / 22.5 / 23 / 23.5 / 23.9 kHz, and the
    sample peak tracks the rms. A 23 kHz tone used to come back as TWO components 1.7 dB apart
    (21.1 kHz −5.33, 19.1 kHz −7.04); they are now −27.45 and **−91.58**.
  - **🔴 The product half: a driven nonlinearity does not MASK those images, it DEMODULATES them into
    the bass.** Through a real high-gain capture, a 20 kHz tone at −18 dBFS came back with a 100 Hz line
    at **−17.65 dBFS — 14.5 dB LOUDER than its own carrier**. It is now **−96.54 dBFS, 60.7 dB under the
    carrier: the line dropped 78.9 dB.** On a clean capture, −39.07 → −117.36. Measured through the real
    `NamStage` by one probe run against both trees.
  - **What it costs.** Latency at 44.1 kHz goes from **3.84 to 61.40 host samples (0.087 → 1.392 ms)**,
    derived from the kernel geometry and measured back from the carrier phase to four decimals; a host
    summing two rate-matching stages reports **122** samples. CPU rises from 0.0129 to **0.1447 %RT**
    for one mono round trip on arm64 and 0.0525 → 0.3239 on x86-64 gcc — **+2.3 % and +2.4 % of what the
    whole `NamStage` already spends on its model**, i.e. the same fraction on both toolchains.
  - **`nam::NamStage::latencySamples()` reports 61 at 44.1 kHz** (**6 at `v0.29.0`**; a 4 that only ever
    existed mid-sprint, see the reported-latency entry below), 96 at 96, 91 at 88.2, 53 at 32,
    and it no longer restates the geometry: it asks `StreamResampler::pairDelayHostSamples (hostSR,
    modelRunSR)` (the nullary `delayInputSamples()` this line used to name is GONE — see the entry
    below; the surviving form takes the two rates). Restating it
    is how the previous formula stayed 2.16 samples wrong for a release cycle. **Consumers that delay a
    dry/bypass path by this number must be re-checked** — **10× larger than the number they were given**
    (6 → 61 at 44.1 kHz) and 16× the delay that was physically there (3.8375 → 61.4000). `orbit-amp` was
    built against this branch and its whole suite passes (seven targets, zero failures), but its
    `src/core/BypassWire.h:37` caps the delay it can carry at 64 samples on a comment quoting a formula
    two generations old: 44.1 kHz fits by three samples and every host above 48 kHz silently
    under-delays its bypass path (96 samples short at 192 kHz). Fix it with the geometry, not a bigger
    constant.
  - **Two contract changes that are not tuning.** An exactly-equal in/out rate now short-circuits to a
    **bit-exact copy delayed by 32 samples** (the 0.99 cutoff is a real low-pass, and a caller asking for
    no rate change must not silently get one) — it was a 2-sample delay. And the kernel is
    **approximating, not interpolating**: it no longer passes through its input samples, so bit-exact DC
    and exact polynomial reproduction are gone (DC within a derived 4.6e−6, measured 3.6e−7), and an
    integer ratio like 96 → 48 kHz is a filtered decimation instead of sample-picking.
  - **At a 48 kHz host nothing changed, and that is gated**: the resampler is not installed at all when
    `|hostSR − modelRunSR| ≤ 0.5`, and the render is bit-identical to the pre-change tree (verified
    cross-tree, FNV-1a over a full render).
  - Full measurement, both oracles and the protocol: [`docs/STREAM-RESAMPLER-COST.md`](docs/STREAM-RESAMPLER-COST.md).

- **BREAKING (behaviour), `dynamics`, `deesser`, `dynamiceq`, `poweramp`, `multiband`, `core`: LAW 11c —
  A PAUSE IS SILENCE.** A call with `nch == 0, n > 0` now advances a stage's SHARED, one-per-instance
  ballistics exactly as `n` samples of digital silence at a live width would, instead of freezing them.
  Law 11(d) already called such a call a GAP IN THE STREAM rather than a no-op — audio time PASSED — and
  a detector standing still through passing time contradicts the sentence that made the grid advance.
  Per-channel memory keeps being dropped exactly as before (11a/11d are untouched), and `reset()` is not
  the answer either: it claims a stream RESTART where the caller said a gap.
  - **What the freeze cost, measured.** `dynamics::Compressor` held **-25.311 dB of gain reduction
    through a full second of gap** where the same second of silence releases to -2.076, and dipped the
    return by **-22.11 dB** on material below its threshold (-24.09 dB through ten seconds);
    `dynamiceq::DynamicEqBand` held -21.774 against -2.985; `deesser::DeEsser` -8.000 against -1.104;
    `dynamics::TransientShaper` moved the return by 1.76-2.60 dB. The loudest is `dynamics::NoiseGate`,
    and it is not a shifted envelope but a state machine that never fired: a gate that has to CLOSE
    through a pause stayed wide open, and a -54 dBFS tone on the return came out at -54 where silence
    gates it to -144 — **89.99 dB, 100 % of the construction ceiling.** All of those are now 0.00 dB.
  - **Seven addresses, chosen by MECHANISM.** The five above, plus `dynamiceq::LaneDynamics` (which
    deliberately DISENGAGED — the "dynamics were switched off" verb, not the "time passed" one; at width
    zero its STEREO lane now runs the control loop on silence, while L/R/M/S take the same "this lane
    stopped" branch they take at width one, because `laneRuns()` gates them on `nc == 2`) and
    `poweramp::PowerAmpStage`, whose one shared sag supply
    and thirteen block-rate glides stopped dead on a gap. `multiband::MultibandProcessor` forwards a gap
    to every band exactly once — it used to forward it TWICE to a bypassed band, invisible under freeze
    and a double clock under this law.
  - **`limiter::TruePeakLimiter` is deliberately NOT included.** It does a full `reset()` on a gap, which
    is wrong under any answer (gain reduction -5.08 dB to 0.00, a 48-sample hole on the return), but its
    ballistics are a sliding lookahead window rather than an exponential, and its per-channel state is
    audio that has not been emitted yet. That is its own item (P29) with its own rule.
  - **Bit-exactness, with its area.** Every call with `nch > 0` is bit-identical to the base — verified by
    hash over all seven stages, every mode x detector x link x control period, with width sweeps,
    self-keyed and externally keyed, and the gain-reduction tap. Two exceptions, each with its own test:
    `NoiseGate` now drops a stopped lane's sidechain high-pass (law 11a, which it never had — worth
    **89.99 dB** on a narrowing as well as on a gap), and `MultibandProcessor` no longer hands a bypassed
    band a row of NULL planes on a narrowing call, **which was a segfault, not a wrong number** — and the
    same fix has a second, non-crashing half: a bypassed band whose planes were already valid used to be
    clocked on the PREVIOUS chunk's split audio at a narrowing edge and is now clocked on digital silence
    (3.12 dB on its meter, 1.6 dB out of silence on un-bypass). The segfault itself needs the band
    bypassed BEFORE its first call at that width — once it has run live, the planes are filled.
  - **Two contract changes at width zero**, both consequences of "a pause is the same call carrying
    silence" and both tested: a `GainReductionTap` handed to `Compressor::process` is now FILLED on a
    zero-width call where it used to be left untouched, and an external key handed to the same call is now
    DEREFERENCED where it previously read nothing — so a key passed at width zero must be valid for `n`
    samples, exactly as at any other width.
  - **Cost.** The silent recurrence is autonomous, so it reaches a bitwise fixed point and the rest of the
    pause is free; and once the detector level reaches `core::kGainToDbFloor` the curve's output is a
    constant, so the per-sample work collapses to one multiply-add. Three of the seven have no such floor
    and run their full body until they park: one zero-width call covering a MINUTE costs `Compressor`
    1.49 ms, `NoiseGate` 0.20, `TransientShaper` 3.65, `DeEsser` 4.08 and `DynamicEqBand` 6.43 — against a
    2.67 ms callback budget. The same minute delivered as 128-sample calls costs **0.00304 ms** in its
    worst call, so an RT caller is unaffected and an offline caller that hands a whole transport jump as
    one call pays it once. `pow(c, n)` is deliberately not used —
    it is a different number from `n` rounded multiplications. The horizon is a property of the TIME
    CONSTANT, not of the pause: 23 609 samples at a 5 ms release, 457 808 at 100 ms, 4 461 677 at 1 s, and
    not reached in 200 000 000 at the coefficient cap.
  - **New API.** `core::kGainToDbFloor` and `core::sameBits`; `EnvelopeFollower::stateWord`/
    `advanceSilence`, `GainReductionFollower::advanceConstant`, `GainReductionPath::advanceSilence`,
    `LinkedDetector::stateWord`/`advanceSilence`, `poweramp::PowerAmpStage::sagDroop`. Additive only.
  - `dynamiceq::LaneDynamics`' park counter is now a saturating `long long`: as a `long` it is 32 bits on
    the MSVC row, where `parked += n` at 48 kHz is signed overflow after 12.4 hours — unreachable while a
    gap disengaged a lane in one step, reachable the moment a pause became something a lane spends.
- **BREAKING (reported latency), `nam`:** **`NamStage::latencySamples()` was 2.16 samples too long at
  44.1 kHz and up to 3.3 at 88.2 kHz.** It reported `ceil(3·hostSR/modelRunSR) + 3` — a guess at "~3
  samples of lookahead per stage" — where the geometry is exact and one line away in the same
  repository: `StreamResampler.h` says the identity ratio "passes the signal with a clean 2-sample
  delay", and that is true at EVERY ratio, because `reset()` leaves 3 leading history zeros with
  `pos = 1.0`, so output *k* reads input position *k·inPerOut − 2*. The round trip is therefore
  `2 + 2·hostSR/modelRunSR` host samples for the cubic kernel that was in the path when this was
  written — **3.8375 at 44.1 kHz, 6.0000 at 96, 5.6750 at 88.2, 2.9188 at 22.05**, which is where the
  2.16 and the 3.3 come from. ⚠️ **That formula does NOT survive this release, and the numbers a host
  will see are the ones in the 64-tap kernel entry above, not the ones this paragraph fixed.** The
  kernel became a 64-tap polyphase sinc in the same release, so the geometry is now `D·(1 + hostSR/
  modelRunSR)` with `D = kHalf = 32` and `latencySamples()` reports **61 at 44.1 kHz, 96 at 96, 91 at
  88.2, 47 at 22.05** — **ten times what the old API reported** (6 → 61) and **sixteen times the delay
  that was physically there** (3.8375 → 61.4000) — not the CUT this paragraph announced, which was not
  one fraction either (6→4 and 9→6 are two thirds, 5→3 is three fifths); unchanged (0) at the
  model's own rate, where the resampler is not in the path at all. Hosts using the reported number for
  delay compensation move by the DIFFERENCE — **55 samples at 44.1 kHz** (6 → 61) and **87 at 96**
  (9 → 96) — and **twice that** in the two shipped hosts, which sum a preamp and a poweramp stage
  (**110 and 174**), landing them at 122 and 192 samples of total reported latency.
  **Slot ALIGNMENT does not move**: it runs on `AlignmentTable::delayOf()` → `blendDelay()`/`lagTail_`,
  and none of those reads `latencySamples()` at all. ⚠️ **But "nothing inside `rigplayer` moves" is no
  longer true, and it stopped being true inside this release:** `RigPlayer::process` holds its own DRY
  leg back by exactly this number (`RigPlayer.h:933`, `dryLatency_.advance (…, latencySamples())`) and
  `warmFor()` spends it too (`RigPlayer.h:1219`), so the player's internal dry path moves with it —
  from 6 samples to 61 at 44.1 kHz. That delay line was added in this same release, after the sentence
  it falsifies was written. **Downstream, audio moves into alignment as well:** OrbitCab delays its dry/bypass path by this same number
  (`src/poweramp/PowerAmpRouter.cpp`, `src/core/CabEngine.cpp`) and orbit-amp does the same at the dry
  end of its crossfade, so the wet path sat at the true 3.84 samples while the dry was held at the
  reported 6 — a 2.16-sample mismatch whose first comb notch fell at ~10.2 kHz during an on↔off
  crossfade (3.0 samples and ~16 kHz at 96 kHz). Against the SHIPPED kernel the residual a rounded
  report leaves is **0.40 samples at 44.1 kHz and 0.00 at 96** (61.4000 against 61, 96.0000 against 96),
  which puts the first notch at 55 kHz — **above Nyquist, so there is none in the band at all**. A host
  summing two stages carries twice the residual (**0.80 samples**, 122.8 against 122 reported) and its
  first notch is 27.6 kHz, still above the band.
  ⚠️ **Three OrbitCab tests pin a number that is now known to be WRONG**
  (`tests/PowerAmpRouterAlignTests.cpp`, three `expectEquals(L, ceil(3·sr/48000) + 3)`). They will fail
  on the next core bump — and so would a fix that pinned `2 + 2·hostSR/modelRunSR`, which is the
  intermediate formula this paragraph installed and the same release then replaced. **In the PRODUCT,
  pin nothing: ask `nam::NamStage::rateMatch (hostSR, modelSR).latencySamples`** (or the stage's own
  `latencySamples()`). ⚠️ **Not the bare geometry:** `StreamResampler::pairDelayHostSamples` answers
  what a pair WOULD cost and returns 64 at a 48 kHz host on a 48 kHz model, where the policy installs
  no resampler at all and the answer is 0 — geometry and policy are two owners, and the second one is
  the one a delay line wants. In the TEST, keep measuring the delay independently, from the carrier
  phase of the shipped round trip, as core's own suite does: an expectation computed by the code under
  test is not an expectation. Restating a geometry by hand is exactly how the old formula stayed 2.16
  samples wrong for a release cycle, and doing it again with a newer constant only resets the clock.
- **`core`, docs:** **`StreamResampler`'s header claimed transparency it does not have, and now carries
  the measurement instead.** ⚠️ **Read this whole entry in the PAST tense: every cost in it is the
  CUBIC's, and it is the case FOR the kernel change announced above, not a description of what ships.**
  The shipped header (`StreamResampler.h:34-44`) already says it that way; only this note did not. The old justification — *"the driven nonlinear stage masks the
  interpolation images"* — had no number behind it, and the quantity that had since been measured was a
  different one. Measured (`docs/STREAM-RESAMPLER-COST.md`, new): a phase-dependent kernel is a linear
  periodically time-varying filter whose per-phase gain has Fourier coefficients `H(Ω+2πk)`, so the
  "amplitude modulation" and the "interpolation images" are **one mechanism**, not two. The NAM round
  trip at 44.1 ↔ 48 kHz costs **−4.17 dB coherent and −9.27 dB worst-phase at 17.64 kHz** (−2.59/−5.14
  at 15 kHz, −5.48/−14.79 at 20 kHz), from a composite period of exactly 147 output samples. **That is
  ONE round trip; OrbitCab runs two `NamStage`s IN SERIES** (preamp → EQ → poweramp, and its own
  `updateLatency()` comment says so), i.e. four of these stages — measured **−9.03/−13.16 at 17.64 kHz
  and −12.09/−17.85 at 20 kHz**, where doubling the decibels would say −8.35/−18.53. On that chain the
  BEST phase falls from −0.61 dB to −5.20, so the top octave is down at every phase rather than only at
  some. `rigplayer` runs its two stages in parallel and stays on the one-round-trip row. That set
  is complete for the shipped priming but is a LINE through the two stages' phase torus, not the full
  product: over all 160 integer alignments the worst phase barely moves (−9.29 against −9.27) while the
  coherent carrier spans −3.59…−6.83, because it is an interference term between the stages. **The decimating direction has no stopband at
  all**: at phase *t = 0* the weights are `(0,1,0,0)`, a bare sample pick, so a tone above the output
  Nyquist survives at −3 dB rms / 0.0 dB SAMPLE peak (a time-domain fact with a kernel reason: at
  `t = 0` the weights are `(0,1,0,0)` so |M(0)| = 1 at every frequency, and the phase grid has points
  within 1/147 of zero where |M| is still −0.002 dB; across the measured rows no spectral line exceeds
  −4.67 dB) and folds back as TWO strong components — `44100 − g` at
  about −5 dB and `g − 3900` at about −7 dB — i.e. across **18.15–22.05 kHz**, not one top slice. Against the model's own
  aliasing floor the OUTPUT leg alone sits 3–24 dB below it on a high-gain capture but **up to +9.8 dB
  above it on a clean one**, at every level from 17.5 kHz up (the whole rate-match: above in 22 of 30 tone × level cells, up to
  +12.3 dB), and **driving harder does not help** — across a 42 dB sweep the error-to-signal ratio is
  flat in 16–22 kHz where the artifacts live and grows +10…+14 dB in 0–4 kHz where they do not, because
  the nonlinearity **demodulates** the input leg's images into the audible range: a 20 kHz tone at
  −18 dBFS into a high-gain capture returns a **100 Hz line at −17.7 dBFS, 14.5 dB louder than its own
  carrier**, against −174.6 dBFS through an ideal round trip. **No kernel change here**: the header now states the cost, the
  new `felitronics_core_streamresampler_lptv_tests` (121 checks) pins the table, the period-147 closure,
  the 0 dB decimation peak and the criterion itself — the round trip adds **−8.84 dBc at 17.5 kHz**, and
  a `tanh` has to be driven to **`tanh(6.2x)`** (bisected) before its own folding reaches that, so below
  a near-square-wave drive the rate-match is the LOUDER artifact. The candidate comparison is in the
  document for the product decision — including the part the two tone axes get wrong on their own: a
  32-tap sinc flattens both axes and is still **5 dB worse in the bass** on real DI through a driven
  capture, because its band edge feeds the same demodulation from a different cause. A 64-tap one at a
  0.99 cutoff is better in every band at every drive, for **+2.5 %** of what the stage already spends on
  the model and **61.4 host samples** of delay against the cubic's 3.84 — that comparison is what the
  DECISION was made against; 61.4 is what this release ships.

- **BREAKING (behaviour + latency), `oversampling`, `saturation`, `limiter`, `poweramp`:** **the shipped
  `tapsPerPhase` default rises from 32 to 64, and the reason is aliasing, not the pass band.**
  `PolyphaseOversampler`'s cutoff is FIXED at 0.90 × baseband Nyquist, so its transition band has to fit
  between 0.45 fs and the fold at 0.50 fs, and `tapsPerPhase` is the only thing that decides whether it
  does. The design DECLARES its own target one line from the taps — `beta = 9.0`, a ~90 dB Kaiser
  stopband — and **at 32 taps it delivered 27 dB.** Everything above 0.50 fs folds straight back into the
  audio band, so that is not a nicety. Measured end to end on a `Saturator` at 0.17 fs and 4×, TOTAL
  non-harmonic energy — every bin that is not a harmonic BELOW Nyquist, so a harmonic that folded to get
  there counts as the aliasing it is:

    drive     +0       +6      +12      +18      +24      +30      +36  dB
    32     −135.3    −60.1    −48.2    −44.3    −31.8    −25.3    −22.9  dBc
    64     −130.8   −132.3   −101.0    −50.3    −31.3    −24.8    −22.4  dBc

  **At the drive this stage is used at — `Saturator::Params` calls 1..6 dB the mastering range — the
  taps remove 50 to 72 dB of aliasing.** They stop helping above about +24 dB, where the harmonic series
  reaches past the OS Nyquist and folds INSIDE the oversampled domain, which no decimation filter can
  reach: that is the oversampling FACTOR's axis, and there the total is 0.5 dB worse, because a flatter
  pass band also delivers what had already folded. The component the taps own at every drive is the
  transition-band leakage — the 3rd harmonic of a 0.17 fs tone lands at 0.51 fs and folds to 0.49 fs,
  where 32 taps left it at −45.5 dBc and 64 puts it at −139.4, into the float noise. (That last pair is
  one tone landing in one of the kernel's nulls; the honest figure for what 64 taps give ANYWHERE in the
  fold region is the worst case below, −90.5 dB, not the null.) Worst rejection over
  the whole fold region, worst of factors 2/4/8: 32 → −26.9 dB, 48 → −51.1, 56 → −76.2, 57 → −82.6,
  **58 → −90.7**, 59 → −89.8, 60 → −90.7, **64 → −90.5**, 96 → −94.3 — and the region is CLOSED at the
  fold, which takes two instruments to read. Below the knee the transition is unfinished and |H| is
  monotone into 0.5 fs, so the supremum sits exactly AT the fold, where a sine projection cannot go
  (it reads its own sampling phase there, ±3.01 dB at 2× and −8.17 at 8×) and the prototype's DTFT
  can: −26.93 / −51.05 / −76.17 / −82.63. At and above the knee the worst is a SIDELOBE inside the
  band, at 0.502–0.505 fs, where the swept projection is exact and the DTFT at the fold reads better
  (64 taps: −99.5 at the fold against −90.46 at the sidelobe). Read as a knee, not a step:
  below ~58 the transition is genuinely unfinished, at ~58 it reaches the Kaiser window's floor and then
  RIPPLES there by about a dB, so 59 and 61 fall a tenth of a dB short while 58, 60 and 64 clear it.
  "The first taps count meeting −90.0" is therefore 58, and that integer is an artefact of a hard bar on
  a rippling quantity. 64 is the round number past the knee; above it further taps buy pass-band width
  rather than rejection — which is why the answer is not 96.
  - **The pass band is the COROLLARY, and it is the half that was already written down.** Two oversampled
    stages in series — a clipper in front of a limiter, the real assembly — cost **−1.549 dB at 17.6 kHz
    and −6.033 at 18.5 kHz at 44.1 kHz** on the old default, i.e. every consumer that built that chain got
    a ~19 kHz lowpass silently. They now cost **+0.000 and −0.610**. `felitronics::mastering` has passed 64
    explicitly since it was written and is **bit-identical** across this change (verified over 942 912
    float32 values, 18 configurations; the same stand shows the DEFAULT paths differing, so it is not a
    blind null). A caller passing `tapsPerPhase` explicitly is likewise bit-identical — with one
    exception, which is the new guard below rather than the default: a topology that used to be accepted
    and is now refused (`factor` above 64, `tapsPerPhase` above 1024) does not "produce the same bits",
    it produces a `false` from `prepare()`.
  - **LATENCY MOVES, and it is host-visible.** Every affected stage reports `tapsPerPhase − 1`, so
    **31 → 63** samples; `TruePeakLimiter` at its default 1 ms lookahead and 48 kHz goes **79 → 111**.
    Read it from `latencySamples()`, which is what `mastering::MasteringChain` already does.
  - **CPU ROUGHLY DOUBLES in the FIR** — "no CPU cost" would be false here. Measured at 48 kHz, stereo,
    4×, block 512: `Saturator` 0.97 → 2.19 %RT, `TruePeakLimiter` 1.10 → 2.22 %RT, the pair in series
    **2.13 → 4.44 %RT** (2.09×).
  - **`poweramp::PowerAmpStage` gets the knob it never had** — `prepare (sampleRate, maxBlock,
    oversampleFactor, tapsPerPhase)` — and takes the same default. It used to hardcode 32 with no way for
    a caller to say otherwise, and its own aliasing gate certified that as adequate **because the gate's
    analysis window stopped at 10 kHz**, below where the transition-band leakage lands. Over the whole
    band the same gate reads **−68.7 dBc at 32 taps, which FAILS its own −70 dBc bar**, against −77.4 at
    64; a 3 kHz fundamental at +12 dB of drive goes −56.2 → −75.9 dBc. The window is widened to Nyquist
    and the suite prints the whole map at BOTH taps counts. The stage's CPU roughly doubles with everyone
    else's — **1.24 → 2.43 %RT** — and its round trip goes 31 → 63 samples (+0.67 ms at 48 kHz). What the
    taps do NOT fix there: the map's hot cells (a 3 kHz fundamental at +24 dB reads −40.9 dBc at either
    taps count) are the tube's harmonics folding inside the oversampled domain, which is the FACTOR's
    axis. OrbitCab keeps its own copy of the stage at an explicit 32 and is unaffected.
  - **`TruePeakLimiter`'s delivered-excess characterisation is RE-DERIVED**, as its own header demanded
    ("widen the pass band and this term has to be re-derived"). The worst tone the round trip passes flat
    is now **2fs/5 rather than fs/3**, so the grid-geometry term goes **+0.302 → +0.436 dB at 4×** and
    **+0.075 → +0.108 at 8×** (unchanged at 2×, where fs/3 still dominates at +1.250). The old figures
    were a property of the lowpass, not of the limiter — `mastering` has been running at +0.436 all along.
    2fs/5 is now a witness in the ceiling matrix, and the CHARACTERISATION battery — everything reached
    through `Setup`/`renderAt`, including the reference oracle it nulls against — runs on the SHIPPED
    topology instead of a hardcoded 32. Roughly two dozen guard and plumbing checks still spell 32 out on
    purpose; those are about refusal and resumption, not about the filter. **The one place the sharper filter costs** is the degenerate corner with
    lookahead AND release both at their floors, where re-band-limiting a step rings more: **+0.36 / +0.28 /
    +0.26 dB at 2× / 4× / 8×**, pinned. The ON-GRID bound, the only thing actually promised, is unchanged.
  - **Source-level API break beyond the defaults:** `PowerAmpStage::prepare` gains a fourth parameter,
    so its *type* changes. Ordinary calls still compile (the parameter is defaulted), but anything that
    names the function's type — `void (PowerAmpStage::*)(double, int, int)`, a `std::function` built
    from it, an explicit `&PowerAmpStage::prepare` cast — does not.
  - **`PolyphaseOversampler::prepare` now refuses `tapsPerPhase` above `kMaxTapsPerPhase` (1024) and
    `factor` above `kMaxFactor` (64)** —
    both factors of `N = factor * tapsPerPhase` had to be bounded because the PRODUCT is what overflows a
    signed int before allocating, and either argument alone can do it. `TruePeakLimiter` guarded both at
    its own gate already; `Saturator` passes an unbounded `oversampleFactor` straight through, and
    `prepare(INT_MAX, 1)` on the default taps is UBSan-confirmed overflow followed by a `length_error` —
    a terminate under the wasm tier's `-fno-exceptions`.
  - **`oversampling`'s own suite measured none of this**: it ran entirely at a hardcoded 32 on tones of
    500 Hz and 2 kHz at 48 kHz (0.010 and 0.042 fs), owning the default and never looking where the
    default decides anything. It now pins the stopband, the aliasing and the two-stage droop surface,
    two-sidedly, with its oracle's liveness asserted.

- **BREAKING (API + behaviour), every module with a block-level `process()`:** **the core had four
  different answers to "the caller passed something other than what was prepared", and every one of them
  was SILENT.** A census over all 24 modules found three answers for the width, four for the length, and a
  case that fits neither. The single rule is now **law 11** in `docs/DSP-ARCHITECTURE.md`, and every
  block-level entry point returns its verdict: **`[[nodiscard]] bool process(...)`**, `true` = "accepted
  and honoured in full". That is the API break — a call site that ignored the result now warns, and warns
  in exactly the place where the contract is decided.
  - **The LENGTH is a capacity, not a limit.** `maxBlock` sizes scratch; `process()` chunks, so any
    `n >= 0` is processed IN FULL and the chunked pass is bit-identical to the caller having chunked it
    itself at the same boundaries. Three stages truncated silently before: **`dynamics::NoiseGate`** let
    everything past `maxBlock` out UNGATED — 3840 of 4096 samples at **+89.99 dB** over the gated ones,
    which is 100 % of the construction ceiling (`-floorDb` = 90 dB); **`nam::NamStage`** let it bypass the
    amp model **bit-identical to its input** — 448 of 512 samples on a 64-sample prepare; and
    **`multiband::MultibandProcessor`** dropped the ENTIRE call, not one sample touched. (Dropping
    `NamStage`'s clamp without chunking would have been far worse than the defect: an unclamped `n` is a
    heap overflow in the backend's scratch and a resize — an allocation — inside NAM, on the audio thread.)
  - **The exception, and it is named in the law:** a two-phase API whose first phase RETURNS a buffer of
    `maxBlock` cannot chunk, because the result has to outlive the call. `NoiseGate::analyse`/`applyGain`
    and `EqEngine::captureSectionInput` therefore REFUSE an over-long call rather than truncate it.
    `NoiseGate` also now records how far the curve is valid: phase B used to be bounded by the buffer's
    CAPACITY rather than by what phase A wrote, so a short analysis followed by a long apply multiplied
    the tail by the PREVIOUS call's curve — measured, **90.0 dB of attenuation on material that was never
    analysed**. `NoiseGate::analysedSamples()` is the new accessor; `NoiseGate::prepare()` now returns
    `bool` and honours `maxChannels`, which it used to ignore.
  - **The WIDTH is a limit: `nch > maxChannels` refuses the WHOLE call, before anything moves.** It used
    to process a prefix in eleven modules. A prefix is not the safer half-measure it looks like — the
    surplus channels are unprocessed either way, and processing some of them only hides the fault while
    the processed ones acquire a latency and a gain the others do not: `limiter::TruePeakLimiter` prepared
    for 2 and called with 4 emitted the surplus **+7.02 dB over its ceiling** AND left the two it did
    process **79 samples late** relative to them, which combs at 304 Hz on any fold-down.
    `saturation::Saturator`'s surplus came out **bit-identical to the input**, with no saturation at all.
    **Affected:** `TruePeakLimiter`, `Saturator`, `EqBand`, `EqEngine`, `Dither` (which was bounded by
    `core::kMaxChannels` rather than by its own prepared width), `DeEsser`, `DynamicEqBand`,
    `TransientShaper` (which ignored `maxChannels` entirely), `PowerAmpStage`, `ConvolutionEngine`,
    `MultibandProcessor`, `LinearPhaseEq`, `NaturalPhaseEq`, `LoudnessMeter`, `TruePeakMeter`, `MonoBass`,
    `StereoWidth`, `RigPlayer`, `NamStage`.
  - **A NARROWER call stays legal** (it is the falling-edge mode P18 defined) **except where it is
    meaningless**, and there it is refused OBSERVABLY: `convolution::MatrixConvolver` and
    `MatrixConvolverNupc` need exactly `channels_` planes — a 2x2 matrix needs both inputs to compute
    either output. They used to drop such a call and write NOTHING, so a caller that pre-zeroed its output
    got digital silence and no way to find out (P18 F35). Their width is EXACT now in both directions: a
    WIDER call used to be accepted with the extra planes left dry. This propagates to `CabConvolver`,
    `LinearPhaseEq` and `NaturalPhaseEq`, whose "mono path" on a stereo-prepared engine had in fact been
    doing nothing at all.
  - **BREAKING (behaviour): a call with `n > 0` and NO channels is a GAP in the stream, not a no-op.** It
    spends audio time (the law-8a grid advances) AND every channel at index >= `nch` counts as stopped, so
    P18's falling edge fires — at `nch == 0`, for all of them. Saying "time passed" and "nobody stopped"
    in one breath reopens exactly the defect P18 closed: measured on untouched `main`,
    `dynamics::Compressor` with 5 ms of lookahead, a tone, 4800 samples of zero-width calls, then stereo
    DIGITAL SILENCE emitted **0.280315 out of the silence (-11.05 dBFS)**, last non-zero at sample 239 —
    the whole 240-sample lookahead line, note for note. `Compressor` also gains the NARROWING half of that
    edge, which it never had. `eq::EqEngine` used to drop a zero-width call outright while the same
    `EqBand` driven directly advanced its grid: on one 500 -> 4000 Hz glide and 10240 samples of such
    calls the two had diverged by **11.08 dB**. Negative `n` or `nch` is malformed and refused, never
    clamped into an index.
  - **`prepare()` is binding too, and refuses what it cannot honour** — an observable refusal in
    `process()` is worth nothing if `prepare()` already lied. `convolution::CabConvolver::prepare()`
    silently clamped `numChannels` to 2, after which `process(io, 4, n)` was a well-formed call that left
    planes 2-3 DRY; it returns `bool` now, as do `NoiseGate::prepare`, `Dither::prepare`,
    `DeEsser::prepare`, `DynamicEqBand::prepare` and `TransientShaper::prepare`.
    `MultibandProcessor` gained a `prepared_` flag (its `prepare()` could already fail and `process()`
    ran anyway), and `ConvolutionEngine::prepare()` now clears the participation ledger it left stale.
  - **`neural::Inference` changed:** the concept requires `process(io, nc, n) -> bool`. Any custom
    inference backend must return its verdict.
  - Two integer-overflow fixes that came with the "any `n`" promise: the chunk loops in
    `saturation::Saturator` and `poweramp::PowerAmpStage` advanced by `maxBlock` rather than by the length
    actually taken (signed overflow near `INT_MAX`, which `TruePeakLimiter` had already fixed), and
    `dither::Dither`'s blank counter added `n` to an `int` before clamping it.
  - **A refused `prepare()` writes NOTHING and leaves the object UNPREPARED.** Both halves cost a defect
    while this law was being applied and are now part of it: storing one argument before validating the
    next left a new WIDTH beside an old buffer (`mastering::OfflineRenderer` — a heap-buffer-overflow
    ASan caught), and refusing before reaching an inner `prepare()` left the object ARMED on its previous
    build (`multiband::MultibandCompressor` — `prepare(2)`, a refused `prepare(0)`, and `process()` still
    ran). **Also `[[nodiscard]] bool` now:** `rigplayer::RigPlayer::prepare` (it CLAMPED the width — the
    literal defect this law describes), `multiband::MultibandWidth::prepare` (its ceiling is 2, not
    `kMaxChannels`, because its band is a fixed stereo stage), `analysis::TruePeakMeter::prepare`,
    `analysis::LoudnessMeter::prepare`, `dynamiceq::LaneDynamics::prepare`,
    `mastering::OfflineRenderer::prepare`.
  - **`stereo::MonoBass` and `stereo::StereoWidth` now refuse `process()` before `prepare()`.** Their
    member defaults looked like a valid configuration and are not: MonoBass's crossover has no
    coefficients until `prepare()` runs, so a default-constructed object passes the side band it exists
    to fold at **-6.02 dB where a prepared one kills it to -54.22** — 48.2 dB at 30 Hz.
  - **The gap reaches per-channel state one layer further down than the first pass saw.** A composite
    must PASS a zero-width call to what it wraps, not answer for it: `multiband::MultibandProcessor` now
    also drops the parallel DRY delay (the whole 240-sample line replayed — 0.25 out of digital silence,
    -12.0 dBFS, invisible at mix 1 where the sum cancels exactly) and notifies BYPASSED bands (-11.83
    dBFS); `dynamiceq::LaneDynamics` releases its own lanes (a lane held -12.000 dB through a second of
    gap and pushed the return down by -5.25 dB); `analysis::LoudnessMeter` clears the K-weighting state
    of a stopped channel, like `TruePeakMeter`'s history (momentary read **-29.19 LUFS** into digital
    silence, against -120.00). New: `analysis::KWeightingFilter::resetChannel(c)`.
  - **`eq::EqEngine::captureSectionInput` is the law's other two-phase API and owes the same:** a width
    above the prepared one is refused rather than narrowed, and a refused OR empty call is inert — an
    `n == 0` capture used to wipe a valid 32-sample one.
  - **Also newly REFUSED where these used to clamp or run:** `prepare(..., maxChannels = 0)` on
    `eq::EqBand`, `eq::EqEngine`, `saturation::Saturator`, `lineareq::LinearPhaseEq`,
    `lineareq::NaturalPhaseEq`, `multiband::MultibandProcessor` (hence `MultibandCompressor`) — a JUCE
    host with a disabled bus passes exactly that; `prepare(..., maxBlock/blockSize <= 0)` on
    `dynamics::NoiseGate`, `multiband::MultibandProcessor`, `mastering::OfflineRenderer` and
    `rigplayer::RigPlayer`; and `process()` BEFORE `prepare()` on `dither::Dither`,
    `deesser::DeEsser`, `dynamics::TransientShaper`, `dynamiceq::DynamicEqBand`, `dynamiceq::LaneDynamics`
    and `dynamics::NoiseGate`, which used to run on their member defaults.
  - **`limiter::TruePeakLimiter` now treats a zero-width call as the channel-count change it is** and
    resets — measured, gain reduction from -5.08 dB to 0.00 and a 48-sample hole on the return. That is
    the discontinuity its header already accepts for a width change; it is called out here because it is
    a THIRD answer to law 11c's open question, arrived at by the width rule rather than chosen.
  - **`convolution::PartitionedConvolver::process` and `NonUniformConvolver::process` also return
    `bool`.** They are mono and take no channel count, so "every block-level entry point" would otherwise
    have been a claim with two exceptions.
  - **`stereo::MonoBass::prepare` and `StereoWidth::prepare` return `[[nodiscard]] bool` and honour the
    `maxChannels` they take** — they used to accept it and ignore it, which is the thing law 11(b) calls
    lying about a contract. `mastering::MasteringChain::prepare` propagates the MonoBass verdict.
  - **Migration:** check the return value. **Two-phase callers first:** `NoiseGate::analyse` now REFUSES
    a block longer than `maxBlock` where it used to truncate, so a consumer that hands it a raw host
    block and ignores the verdict gets the WHOLE block ungated instead of only the tail. Use the fused
    `process()`, which chunks, or chunk at the call site. `if (! stage.process (io, nch, n)) { /* your geometry is wrong */ }`
    A consumer that today passes more channels than it prepared for, or a block longer than a
    capacity-bearing `analyse()`, was already getting broken audio — it was just not being told.

- **fix(core, eq, stereo):** **law 8's denormal flush ran once per `process()` call, so the caller's
  block size decided where a numerical event landed.** New `core::StateGrid` — a phase counter over
  AUDIO samples, period 64, re-anchored by `reset()` — and `eq::EqBand` (hence `EqEngine`) and
  `stereo::MonoBass` now do their periodic maintenance there. The rule is written up as **law 8a** in
  `docs/DSP-ARCHITECTURE.md`; two modules had already reached it independently (`analysis::LoudnessMeter`
  per 10 ms sub-hop, `saturation::Saturator` per sample).
  - **BREAKING (behaviour), `eq` and `stereo`:** output moves. On settled parameters the change is
    confined to the sub-audible: 8.3 % of samples over a 27 M-sample sweep, worst difference 1.95e-13,
    loudest differing sample −127 dBFS. With a parameter ramp in flight it is large and intended
    (worst 1.25 full scale) — the smoothers now advance by exactly one grid period per tick instead of
    by the whole call, which is what makes a ramping band slicing-invariant at all.
  - What it buys, measured: a whole-file render of a bell into silence went from **38 522 of 40 000 tail
    samples differing** from a one-sample-at-a-time render to **0**; from **37 678 subnormal tail
    samples** (the 10–100× stall on any CPU without hardware FTZ) to **0**; and one `+Inf` input sample,
    which used to poison **479 900 of the next 480 000 samples** of a whole-file render, is now healed
    within one grid period (28 samples). A host with a block SHORTER than a period keeps its immediate
    recovery: the `isfinite` half of the flush still runs at the end of every call
    (`eq::Biquad::healPoison()`, `eq::Svf::healPoison()`, `eq::Crossover2::healPoison()` — additive).
  - **BREAKING (behaviour), `dynamiceq::LaneDynamics`:** it drives `eq::EqBand` in 16-sample control
    chunks, so the band's STATIC freq/Q/gain glide used to advance every 16 samples and now advances
    every 64. Measured on a 500 → 5000 Hz +12 dB edit with 30 ms smoothing: designs per 100 ms
    300 → 75, largest single step 2.15 → 3.22 dB — still finer than any real host block gave before
    (a 512-sample host took 10 steps of 11.8 dB), but it moved, and a dynamic band is where edit
    smoothness is most visible. The gain DELTA itself is unaffected: it is an arrival, still consumed
    at the producer's 16-sample cadence.
  - **BREAKING (behaviour), `stereo::MonoBass`:** settling into the full-wide bypass is now decided per
    sample instead of at the top of the next call. The samples in between used to take the M/S round
    trip, which is not the identity in float — 1 LSB of 24 bit (5.96e-08) on 22 953 samples of the
    re-slicing sweep, and it made the output depend on where the caller cut.
- **BREAKING (behaviour), `eq::EqBand::reset()`:** it is now a real STREAM RESTART — it re-snaps the
  freq/Q/gain smoothers, redesigns, and clears `initialized` so the first parameter write after it snaps,
  exactly as the first write after `prepare()` does. It used to clear filter state and leave the
  smoothers mid-glide, so a second render of the same programme started from a different design: measured
  **0.51 full scale** against a freshly prepared chain. `mastering::MasteringChain` carried that as a
  `forceSnap_` workaround (writing every band with its lanes off and then writing them back); the
  workaround is deleted, and the chain is bit-identical without it across 20 scenarios × 1 152 000 samples.
- **NEW API, `eq::EqBand::clearAudioState()` / `eq::EqEngine::clearAudioState()`:** a STOP — clear what
  the previous audio left behind (filter memory, participation ledgers, dynamic seams, design-key caches)
  and leave the parameter epoch and the grid phase alone. This is what a consumer skipping the engine for
  a while actually wants at a bypass edge; `reset()` there would now also snap every ramp in flight and
  turn the next parameter write into a hard step. **`mastering::MasteringChain` was updated; an external
  consumer that calls `EqEngine::reset()` on an EQ-off edge (OrbitCab's `AmpEq::process`) should switch
  to `clearAudioState()` to keep today's behaviour.**
- **fix(core):** `core::Smoother::advance(n)` memoises `pow(coeff, n)` on `n`, dropped whenever `coeff`
  changes. Bit-identical by construction (`pow` is a pure function) — pinned by a test, and by the
  bit-exact null of every existing consumer — and it is what keeps a transcendental per parameter per
  grid period off the audio thread.

- **fix(eq, saturation, dynamiceq, deesser, poweramp, convolution, multiband, dither):** **a filter that
  is not called does not decay — it FREEZES**, and replays a signal from before the gap when it is called
  again. Eleven instances of one shape: a cell of per-channel or per-lane recursion behind a gate the
  surrounding object keeps running through. All eleven were present in every release up to this one.
  Measured out of DIGITAL SILENCE on the input: `ConvolutionEngine` 7.44e-01 (**−2.6 dBFS**),
  `MultibandProcessor` 1.55e-01 (−16.2), `DeEsser` in SplitBand 6.25e-02 (−24.1), `PowerAmpStage`
  1.67e-02 (−35.5) on its own knob gates and 0.649 (−3.8) on a channel change, `EqBand` **+8.12 dBFS**,
  `LaneDynamics` 11.97 dB of unearned gain reduction, `DynamicEqBand` 0.388 eleven samples in,
  `Saturator` 0.9337, `EqBand`'s swept branch 0.690, `Dither` six LSB of 24 bit.
  - The rule is **the falling edge of PARTICIPATION**, per cell: a cell that ran on the previous
    sample-bearing call and does not run on this one loses its sample memory, whatever stopped it. NOT
    "the channel count grew" — `EqBand`'s single-signal lanes gate on `nc == 2` exactly, so `2→3→2`
    retires them on the INCREASE and revives them on the DECREASE, and the compressor's `nc > lastNc_`
    rule never fires there. Only sample memory is dropped: gain seams, smoothers and applied-value caches
    are current CONTROL and are preserved.
- **BREAKING (behaviour):** the gates above stop and restart at a CONSTANT channel count too, so fixing
  them changes what those transitions sound like. Named, because each is now pinned by a test:
  - a whole-band idle in `eq::EqBand` clears signal history only — it used to call the full `reset()`, so
    an inert companion lane on 0 dB decided whether another lane's commanded gain survived. The full
    reset, seams included, now belongs to explicit `reset()` alone.
  - switching a lane off, flipping `swept`, or turning `dyn.on` off now clears that path's delta filter
    state; `Saturator`'s Asym↔symmetric toggle clears its DC blocker; `PowerAmpStage`'s presence, depth,
    load and iron gates clear their filters and its sag envelope.
  - a parameter ramp in `eq::EqBand` now advances while the band is idle. It always did for lanes that
    were not running; the fully-idle band was the one case that fell out of the rule, so the same edit
    landed at two different times depending on whether an unrelated lane happened to be on.
  - a parked lane's PROGRAMME ESTIMATE in `dynamiceq::LaneDynamics` is now duration-aware: kept below a
    quarter of the estimator's own averaging constant, fast-adapted below four times it, discarded beyond.
    Keeping it unconditionally cost 17.87 dB for 3.36 s against a lane that heard the change — in boost
    mode, 18 dB of unearned BOOST. The trade is named rather than hidden: after a long SILENT park the
    first loud material is no longer ducked.
- **BREAKING (eq, source):** `EqEngine::prepare` is now `[[nodiscard]] bool` and REFUSES a sample rate it
  cannot honour; `EqBand::prepare` returns `bool` and a refused band stays inert. It validated nothing
  before, and `fs` reaches every coefficient through `tan(pi*f/fs)`: measured, `prepare(0)` and
  `prepare(NaN)` each put 63 of 64 output samples non-finite on a 0.25 input, and `prepare(1.0)` does the
  same to a HighPass or a Notch. Only the mastering chain was protected, because it validates at its own
  entrance. This is the second half of the fix `Saturator::prepare` received in v0.26.0.
- **BREAKING (eq, behaviour):** `EqEngine::captureSectionInput` records the WIDTH it captured, not only
  the length, and hands back `nullptr` for the columns outside it — the same refusal it already made for
  an over-long block. They used to point at the PREVIOUS block's audio, so a consumer reading one column
  too far detected on a signal that was not that block's, silently and on plausible data.
  `sectionInputChannels()` reports the width.
- **feat(eq, oversampling):** `resetChannel(int)` on `Svf`, `Crossover2`, `MultibandSplitter` and
  `PolyphaseOversampler` — purely additive. Clearing one channel's state without restarting the channels
  that never left is what the isolation half of the fix above needs, and none of them could do it.
  `PolyphaseOversampler::resetChannel` clears the ring POSITIONS with the samples.
- **perf(eq):** `Svf::flushDenormals` visits the PREPARED channels instead of `kMaxChannels`. The columns
  past `ch` are never written and `reset()` zeroes them, so visiting them was a no-op that cost the same
  as real work — since v0.26.0 the body is two `isfinite` tests and two stores per column rather than a
  flush, so a mono `Svf` paid for sixteen, once per block, in the primitive every band, lane, crossover
  and probe calls.

- **🔴 BREAKING (`core`): `StreamResampler::delayInputSamples()` now TAKES THE TWO RATES** —
  `delayInputSamples (inRate, outRate)`. It ignores them today and returns `kHalf` either way; the
  signature exists so that making the kernel length depend on the ratio (the open `kTaps` item) edits
  one body instead of every consumer's arithmetic. The nullary form is GONE, and downstream callers
  must pass the rates. **`felitronics-core` has no nullary caller left; `orbit-amp` has four.**
- **`core` + `nam`: ONE OWNER FOR THE RATE-MATCH DELAY, and the restatements are gone.** The composition
  is `core::StreamResampler::pairDelayHostSamples (hostSR, modelRunSR)` — pure geometry, "what a
  down+up pair costs in host samples", with no opinion about whether one is installed. The POLICY is
  `nam::NamStage::rateMatch (hostSR, modelSR) -> {modelRunSR, resampling, latencySamples}`: normalise an
  unknown model rate, then gate at half a hertz, then round to nearest — and the ORDER of those three is
  part of the contract. `nam::NamStage::kModelSampleRate` (48 kHz) is public, because a consumer sizing a
  delay line before any model exists could not previously even name the number and one invented an
  8 kHz stand-in and sized itself wrong. Seven restatements across three modules now ask instead.
- **`rigplayer`: a host rate is judged in ONE place, and every buffer derived from one is bounded.** New
  `RigPlayer::usableSampleRate (hostSR)` and `RigPlayer::kMaxSampleRate` (3.0e6 — the same ceiling
  `dynamics::Compressor` and `limiter::TruePeakLimiter` already use, for the same stated reason). A rate
  outside `(0, kMaxSampleRate]` falls back to 48 kHz exactly as a non-positive one always did. This is a
  FIX, not a tidy-up: the dry-aligner capacity became `(int) ceil(<a function of the rate>)` earlier on
  this branch, and a guard spelled `isfinite` does not make that conversion safe — measured through
  `prepare()` with UBSan, `prepare(1e300, 64, 2)` returned **true** while converting out of range and
  then overflowing `INT_MAX + 2`, and asked the heap for **16 GiB** in that one call. It now asks for
  none. The capacity itself is `RigPlayer::dryAlignerCapacity (hostSR)`, pure and public so that it can
  be pinned at rates this repository does not run: the floor at the shipped 256 had otherwise no test in
  the tree that could see it at all.

## v0.29.0 — the loudness tag answers for the slot that is sounding, and the chain gives the same bits however you cut it (`rigplayer`, `mastering`)

- **feat(rigplayer):** `soundingLoudness()` — the tag of the model actually carrying the sound, with
  `tagged`, `blended` and `slot` beside the number. `modelLoudness()` / `modelHasLoudness()` always
  read slot 0, and slots are handed out by knot parity: on an odd capture the whole sound comes from
  slot 1 while slot 0 holds a neighbour at zero weight, so a host drawing "the tag of the model
  sounding" drew the silent one's — or called a tagged model untagged. Three reviewers found it
  independently.
  - The slot is chosen by the APPLIED weight, not the requested one: while a model is loading or
    warming, the old slot is what is audible, and a face has to name what is heard.
  - `blended` is set only when the weight is strictly between the ends AND the two slots hold
    different files — at rest both slots hold the same capture, and one number is then the whole
    truth. During a real crossfade there are two tags in the sound, and an API that returns one
    without saying so invites the same quiet wrongness this readout exists to end.
- **BREAKING:** `modelLoudness()` and `modelHasLoudness()` are REMOVED rather than kept as aliases.
  Their documented meaning is "slot 0" — the defect itself — so an alias would go on answering wrongly
  for anyone who did not recompile.

<!-- `mastering` shipped in this tag (PR #132, merged before the v0.29.0 release PR #134) but its
     entry sat under Unreleased until it was moved here, so this file claimed the module was not
     released for three versions. The text below is that entry, unchanged. -->

- **feat(mastering):** a new module, and the first one that is a **composition** rather than a stage:
  `MasteringChain` — `gain → EQ → [M/S mono-bass] → compressor (optional internal sidechain HPF) →
  [soft clipper] → gain → true-peak limiter → dither` — plus `OfflineRenderer` over it. No new DSP; every
  stage already shipped. What is new is the part composition kept getting wrong.
  - **A FIXED INTERNAL QUANTUM, not a promise.** The chain never hands a caller's block boundary to a
    stage: it buffers and calls every stage with exactly `internalBlock` samples. That is not tidiness,
    it is the only thing that makes "same input, same output, whatever the block size" TRUE here.
    Measured on this tree: `eq::EqBand` and `stereo::MonoBass` flush filter state once per *call*, so
    renders at block 4096 and block 1 differ in 4721 and 8908 samples; a band parked at 0 dB moves that
    divergence out of the tail and into the programme (2083 samples of a 50000-sample tone, from sample
    24); `dither::Dither`'s auto-blank then amplifies 1e-15 into **three LSB of a 24-bit master** by
    deciding whether the tail is exactly zero; and the compressor is not exempt either — a key of
    `[1.3e-15, 1.5e-12]` with the RMS window set so the follower coefficient is exactly 0.5 gives
    0.478396237 in one 2-sample call against 0.478396297 in two 1-sample calls. Costs `internalBlock`
    samples of declared latency and, measured, **no CPU at all**: 2.26–2.30 %RT flat for every quantum
    from 16 to 4096, so the size is chosen for latency and never for speed. It also closes a Law 8 hole
    nobody had looked at: a whole-file call left 37678 of a 38000-sample tail SUBNORMAL, because
    `FlushToZero`'s "a block is too short to re-traverse the gap" is false for a big block.
  - **Latency is READ BACK from the stages, never computed.** A chain that derives the number itself can
    agree with its own bypass aligners while disagreeing with the stage — `Compressor::prepare(..,
    maxLookaheadMs = 50)` caps the lookahead at 2400 samples, so a chain asking for 60 ms and computing
    `lround(2880)` would hold both at 2880 and pass a bypass null while sitting 480 samples out. The
    suite finds the delay by SEARCHING for the shift that nulls, and measures the active chain's group
    delay from a DFT bin's phase at 25 Hz — a frequency whose period exceeds twice the delay, so the
    answer is not ambiguous modulo a period.
  - **Bypass is decided per stage by measurement, not by uniformity.** The compressor gets an exactly
    transparent WARM bypass out of its own curve (`ratio = 1` → slope exactly 0 → gain exactly `1.0f`,
    sign of zero included, detector still tracking). The clipper and the limiter are skipped with a
    `core::DryAligner` holding their PDC — the clipper's own `mix = 0` is bit-exact for ordinary audio
    but normalises `-0.0f`, and the limiter has no bypass at all: at a ceiling of +60 dBTP it still costs
    the 0.90·Nyquist round trip. Latency never moves when a bypass is toggled.
  - **`tapsPerPhase` defaults to 64**, against the stages' own 32. The prototype loss is a round trip, so
    it doubles in dB, and clipper + limiter in series double it again: at 44.1 kHz that is **−1.549 dB at
    17.6 kHz**, −6.033 at 18.5 and −16.131 at 19.4. At 64 taps: +0.000 / −0.610 / −10.182, for +64
    samples. The table is pinned in the suite — the limiter's own tests stop at 0.357·fs and did not
    notice that its header's droop figures are one filter pass labelled as the round trip.
  - **The channel count is EXACT.** `process()` refuses any count but the prepared one, touching neither
    buffer nor state, so `[A, refused, B]` is bit-identical to `[A, B]`. The stages disagree among
    themselves — the compressor refuses a wider call, the limiter, saturator and EQ process a PREFIX,
    mono-bass ignores anything that is not exactly two — and reconciling that from outside is the
    field-mapping that falls out of step.
  - **One input gate**, in the shape `Saturator` and `TruePeakLimiter` already use, so a poisoned sample
    is bit-identical to the sanitised one it stands for, whichever stages are on. And the chain validates
    its own sample rate positively: `Saturator::prepare(NaN)` returns **true** and emits NaN, and
    `EqEngine::prepare` does not validate at all.
  - `OfflineRenderer` is defined by a formula rather than a description: `out[n] = y[n + D]`, where `y`
    is the chain's output for the input followed by `D` zeros. Same length as the input, aligned, and the
    last `D` frames present — the ffmpeg chain this replaces never emitted them, so a click 4 ms before
    the end of a file disappeared.
  - **Coverage: 193 checks across two binaries, and 27 of 27 mutations of the module are caught.** The
    second binary exists because the obvious latency test is unsound on its own; everything in it
    re-derives the answer from outside the chain. Two of the mutations found real defects while the suite
    was green — `reset()` resumed an interrupted EQ parameter ramp (0.51 of difference, full scale,
    against a fresh chain) and a bypass toggle read a cold aligner.
- **feat(stereo):** `MonoBassParams` + `MonoBass::setParams()` / `params()`. Additive; `setParams` calls
  the three existing setters, so every clamp and rejection is unchanged and the audio is bit-identical
  (pinned). `params()` reads back the RESOLVED values — a corner below 20 Hz comes back clamped.
- **docs(ADR §4):** the "no cross-module glue in the core" rule is amended to what the tree has actually
  been doing since `dynamiceq`: a composite belongs here when IT is the unit under test; the voicing
  stays in the product.

## v0.28.0 — the host states the whole number, and the player does the subtracting (`rigplayer`)

- **feat(rigplayer):** `setHostInputDb()` / `setHostOutputDb()` — a host with a fader of its own now
  states the WHOLE level it wants a device played at, and the player uses it in place of the pack's
  `chain[].input_db` / `chain[].output_db`. `std::nullopt` (the default) hands the level back to the
  pack, so a plugin that only plays packs is unaffected and needs no call.
  - v0.27.0 left the arithmetic with the host: it sent `hand − what the pack says`, applied outside
    the player. That cannot be made correct. The subtraction has to know which pack is loaded at the
    moment it happens, and a host reads that from its own document — which diverges from the loaded
    pack whenever somebody edits it, a rebuild is in flight, or a device is switched mid-build. Each
    of those left the level silently wrong, in the same class as the double application the two keys
    were added to end. The player is the only place that cannot disagree with itself about which pack
    it holds. `hostInputDb()` / `hostOutputDb()` read the hand back.
  - The hand survives a `load()`: it belongs to the bench, not to the pack.
- **fix(rigplayer):** a pack swap no longer publishes unity levels in the middle of itself. `load()`
  called `unload()`, which resets both gains, and republished the real ones only at the end — a
  callback landing in between heard neither pack's level. The stage's levels are now published as
  soon as the stage is known, and one function decides what is applied.
- **test(rigplayer):** a fixture that is not a pure scalar. Every model in these tests was a linear
  gain, and through a pure gain the two levels are indistinguishable — moving where the player applies
  them would have failed nothing. A Linear NAM WITH BIAS adds an offset that input scaling leaves
  alone and output scaling takes down, so the release's central claim (input is drive, output is
  volume) is now something a swap would break.

## v0.27.0 — a pack carries its own two levels, and normalizing is the default (`rigplayer`, `nam`)

- **feat(dynamics):** the gain-reduction path became **one object**: `GainReductionPath` carries
  detector → curve → ballistics together, `CompressorParams : GainReductionParams : DetectorParams`
  completes that layering, and `GainReductionTap` publishes the per-sample gain reduction. Alongside them
  `dynamics::offline::{QuantileHistogram, EnvelopeAnalyzer, ThresholdSolver}` — the detector-domain
  analysis that finds the threshold delivering a wanted gain reduction by MEASURING rather than by
  inverting a curve. None of this was declared when it shipped; the refactor itself was bit-identical to
  the previous compressor across 23 328 configurations × 210 million samples.
- **BREAKING (dynamics, source):** completing that layering moved `thresholdDb` and its neighbours into
  `GainReductionParams`, so **`decltype (&CompressorParams::thresholdDb)` is now
  `double GainReductionParams::*`, not `double CompressorParams::*`**. Exact-type reflection,
  serialisation tables and any `T C::*` template argument spelled with the derived class stop compiling.
  The header says so at `Compressor.h`; this file did not. See also the v0.26.0 note above, which is the
  first half of the same change.
- **feat(rigplayer):** a pack states **how hard it is fed and how loud it leaves** — namz 4.1.0's
  `chain[].input_db` and `chain[].output_db`, applied by the player and by nothing else. Packs are not
  balanced against each other (a Big Muff leaves some 12 dB louder than a clean preamp, and a boost is
  hotter still than the preamp it feeds) while inside a pack the models sit within ±0.4 dB of one
  another, so the level that differs belongs to the DEVICE. Until now the only place to put it was
  `files[].input_db`, spread across every entry, which made one key mean two things at once.
  - The input level goes **first, ahead of the dry copy**, and deliberately NOT into `chainGain_` where
    `extendDb` lives: a blend knob mixes one guitar with itself, and feeding the models less while the
    dry side is fed as ever turns the mix into two instruments at two volumes. `extendDb` stays where it
    is — it is a trick played inside the pack past the top capture, which the dry path leaving at the
    input jack never saw.
  - The output level goes **last, after the mix**. One number for the whole stage, and applied any
    earlier it would ride the wet side alone and move the blend the pack states for that position.
  - Neither is behind `setInputTrims()`, which keeps gating exactly what it always gated:
    `files[].input_db`, the trim of one alias against its neighbour. `stageInputDb()` /
    `stageOutputDb()` read back what the player is applying — **a host with a fader of its own must
    send only its deviation from these**, or the level lands twice.
- **change(rigplayer):** `setNormalize()` now defaults to **true**. A model's `metadata.loudness` tag is
  a contract, not a listener's option: with it off, every capture plays at whatever level the hardware
  happened to give, and no two packs can be compared at all. A host that wants the old behaviour must
  now ask for it.
- **build(nam):** namz pinned at **v4.1.0**, which is where the two keys are.

## v0.26.0 — a tier that stops being aspirational, an external key for the compressor, and six kernels that finally arrive at zero (`build`, `dynamics`, `analysis`, `saturation`, `tools`)

- **feat(build):** the **`wasm-audio` tier is a gate**, not a paragraph. DSP-ARCHITECTURE.md §2 had
  named it since it was written and always marked it aspirational. There is now a `wasm-audio` CMake
  preset and a CI job that builds every default module for wasm32 with exceptions and RTTI off and no
  pthreads, runs the whole suite in node, and audits every emitted artifact — green under a checked
  `SAFE_HEAP` configuration too. What it does and does not prove is written down in
  `docs/WASM-AUDIO-TIER.md`: a thread is not a build error, and law 8 is not covered by it.
- **fix(core):** the one `throw` that closed the exception-free tier. `core/Fft.h`'s
  `SeamAllocator::allocate` held the only `throw` in the shipped core, and `-fno-exceptions` makes a
  `throw` a hard PARSE error — so the tier the ADR promises could not be built at all. Guarded the way
  libc++ guards `__throw_bad_array_new_length`: throw where exceptions exist, `std::abort()` where they
  do not.
- **BREAKING (dynamics, source):** `CompressorParams` stopped being a flat struct and became
  `CompressorParams : GainReductionParams : DetectorParams`, so its own fields moved into base classes.
  Three things break, and this release did not say so:
  - **Designated initializers stop compiling.** `CompressorParams{ .ratio = 4.0 }` is now
    *"field designator 'ratio' does not refer to any field in type 'CompressorParams'"* — a designator
    may only name a DIRECT member, and `ratio` lives in `GainReductionParams` now.
  - **Positional brace-initialization stops compiling too**, loudly rather than silently: brace elision
    fills the base first, so the pre-existing `CompressorParams{ -12.0, 4.0 }` now tries to initialize
    `DetectorParams::Detector` and `DetectorParams::LinkMode` from doubles. Loud is the good outcome —
    the alternative would have been the same spelling quietly assigning to different fields.
  - **`std::is_standard_layout_v<CompressorParams>` is now false** (it was true), which matters to
    anything doing `offsetof`, C interop, or exact-layout serialisation. It remains an aggregate and
    trivially copyable.
  - Assigning to the fields by name (`p.ratio = 4.0;`) is unaffected, which is how most callers write it
    and why this went unnoticed.
- **feat(dynamics):** an **external key for the compressor**. On a mastering bus the kick and the bass
  decide the gain reduction of the whole mix — defect number one of the ffmpeg chain this core
  replaces, and it would have been reproduced exactly, because the ADR rightly forbids putting an EQ
  inside the compressor. So the EQ does not come in: the input goes out. Eleven ways the module could
  not be trusted with a key are closed with it.
- **fix(dynamics):** the limiter's promise matches its measurement — **twenty** contract defects, against
  the six the task named. The worst was on no list: only the channels passed to `process()` advanced
  their state, so a channel-count change mid-stream went **+19.8 dB** over the ceiling. Topology moved
  to `prepare()`, input gated, oversized blocks chunked, floors stated.
- **fix:** **law 8, finished** — six kernels that never actually arrived at zero, not the three the
  earlier pass left open, and the entry recorded as harmless was the most expensive of them. Stated
  properly the mechanism is not "denormals": `x <- t + r*(x-t)` never ARRIVES, because a residual of
  `k` ulps maps to itself for every `k <= 0.5/(1-r)`, so the decay stops dead while the state is still
  finite and keeps feeding whatever is downstream.
- **fix(analysis,multiband):** silence stops costing **54x** more than music. `KWeightingFilter` was the
  last kernel in core+analysis+oversampling with no software denormal flush; on zero input its two
  TDF-II biquads reach a subnormal pair that maps to itself exactly and freeze there — measured out to
  600 s. On a CPU with no hardware FTZ, which is exactly what a browser gives, that is what silence
  then costs.
- **fix(analysis):** one bad sample stops making the loudness meter **lie quietly**. K-weighting is an
  IIR, so a single NaN made its state NaN forever, every later 400 ms block energy NaN, and
  `NaN > absT` false — so the absolute gate silently dropped all of them and the meter averaged only
  the part that predated the NaN, reporting a healthy, plausible number for a programme it had stopped
  measuring.
- **feat(analysis):** the **pre-gate block energies**, for a comparison that cannot jump. Integrated
  LUFS is discontinuous in its own inputs — BS.1770's gates are strict comparisons, so a block within
  ~1e-12 of a threshold flips inclusion between two builds and moves the reading by ~0.01 dB. No
  scalar tolerance on a gated quantity is a robust equivalence measure; the energies before either
  gate are continuous in the input samples, and cross-toolchain work needs that surface.
- **test(analysis,tools):** true-peak gets its **first external judge**. Both paths had been tested only
  against their own documented behaviour; EBU Tech 3341 §2.6 defines an acceptance envelope nobody here
  wrote, and it is now applied — tests 15-23 of the official EBU Loudness Test Set v05 pass on both
  implementations.
- **fix(build):** **law 10** — say what FP contraction we want instead of inheriting three answers. The
  arm64 Linux CI row went red on its first run on untouched `main`, which is what it was added for:
  `a*b + c` may fuse into one FMA with a single rounding, and the toolchains disagree about when.
  `-ffp-contract=on` is stated; it costs the shipping tier nothing, because it is already clang's
  default.
- **test(saturation):** a noise floor the build **measures for itself**, and the DC blocker's pole and
  numerator, neither of which was tested anywhere. The poison-containment floor had been a constant
  fitted twice to whichever toolchain last went red; below it the assertion has no defined answer,
  because the state residual is a random walk with an absorbing zero. It is now a control run — every
  input sample moved one ulp — because a real build's own noise floor reads 20 ulp = 1.19e-06, above
  the constant the line used to carry. Mutation testing then showed the corner frequency could be
  hardcoded and the numerator's zero moved without a single one of 197 checks noticing; both are
  asserted directly now, two-sided.
- **feat(tools):** `fcore::Probe` — one measurement body for the native reference and for wasm — plus
  the C ABI, build recipe, parity harness and page behind it, and an artifact audit that proves
  no-threads from the binary and keeps the two JS sides from drifting.

## v0.25.0 — a knob that clicks states its filter, and NAM is pinned to a release (`rigplayer`, `nam`)

- **feat(rigplayer):** a switch shipping its bands **per position**. namz 4.0.0 gave the format
  `positions[].sections` — each position states the filter it IS, with no travel law, because a
  switch's positions are words with an order and no angle for a gain to travel on. The player had no
  path to it: such a knob went down the curve path, met an empty grid, matched `0 == 0` on the length
  test and bypassed the filter. Flat, with every check along the way passing. `sectionsAtValue` builds
  that position's bands by name and `bandsPerPosition` sends the knob to the band path.
- **fix(rigplayer):** `setSwitch` refuses a value no position declares. It used to store any word,
  report success and read it back on the knob while the sound was the reference — a typo of one letter
  was accepted and heard as nothing. A dial still takes any degree of its travel, swept or not.
- **build(nam):** NeuralAmpModelerCore is pinned to **v0.5.4**, a release rather than a commit.
  `b5a68c3` was the tip of upstream main when OrbitCab pinned it and stopped being so forty-five
  minutes later; it was inherited here with the NAM path and never moved. A real model through a
  second of deterministic signal is **bit-identical** across the two pins, so nothing already captured
  sounds new.
- **build:** namz is fetched at **v4.0.0**, the schema the player now speaks.

## v0.24.0 — the constant-Q analyzer for half the CPU, and silence stops costing more than sound (`analysis`)

- **feat(analysis):** `MultiResSpectrumPaneFast` — the same pane as `MultiResSpectrumPane`, computed
  differently. A `perf` profile said **43.98 % of a tick was inside libm**, and the band integration
  everyone assumed was the cost was 7.6 % while the FFT itself was 4.6 %. Three exact changes remove
  35 349 of the tick's 37 111 transcendental calls: the peak trace is kept in **power** instead of dB
  (which had cost a `log10` per bin to make the value the peak law compares against and then an `exp`
  per bin to undo it — a transform and its own inverse); each column's geometry (owning tier, seam,
  blend, fractional bin edges, display tilt) is derived **once into a fixed-size plan and shared by
  the fill and the peak**, which previously each computed all of it; and DC and Nyquist are peeled out
  of the magnitude loop. **1.82× on an M5 Pro and 2.33× on an i9-13900H with pffft** (1.20× / 1.34× on
  the scalar FFT — the same absolute saving, since none of it depends on the transform). On x86 the
  fast pane now costs less than a single classic 16384 pane. It is a **sibling**: `MultiResSpectrumPane`
  is untouched, and the two are held together by a paired NULL in which **the fill is bit-identical**
  and the peak is inside the sibling's own float-dB quantisation. Two divergences are pinned rather
  than hidden — a negative `peakFallDb` is clamped instead of reaching an infinity, and below −120 dB
  the band integral's own conditioning (a difference of two prefix sums scaled by the loudest bin in
  the tier) stops the sibling being a reference at all, which a test measures rather than asserts.
  New: `docs/PERF-ANALYZER-MULTIRES.md`, with the profile, the tables, the machines they were taken
  on, and the lossy ideas that were measured and rejected.
- **fix(analysis):** `MultiResSpectrumPane` no longer leaves its state stuck in the subnormals on
  digital silence. The fill smoother is `pw ← (1−c)·pw`, and in float32 that descent does not end at
  zero: at the smallest subnormal `c·pw` rounds away and `pw` stops moving, leaving all ~10 755 bins
  doing subnormal arithmetic on the message thread with nothing setting FTZ/DAZ. On x86 the pane
  therefore got **slower the longer the transport stayed stopped** — 505.4 µs on settled silence
  against 167.6 µs with FTZ forced, a **3.0× penalty**, while on ARM it was free and so invisible.
  The smoothed power is now flushed to a true zero below `1e-30` — deliberately **not**
  `core::flushDenormal`, whose `1e-15` is an amplitude threshold and would erase the −150…−200 dB bins
  this pane deliberately sums. No reading moves. New `tierBinPower` accessor: `tierBinDb` floors at
  −200 and could not tell a bin that reached zero from one stuck at `1e-45`. The classic `SpectrumPane`
  never had this — it smooths dB, which is already floored.
- **docs(analysis):** `ANALYZER-MULTIRES.md` drops two claims the code stopped honouring in v0.22.2 —
  a shelf repeating the last band above Nyquist (it reads the floor), and a −120 dB floor (it is −200).

## v0.23.0 — the loudness meter meets Tech 3341 (`analysis`)

- **fix(analysis):** `LoudnessMeter` momentary and short-term now accumulate on a 10 ms sub-hop
  (M = the last 40, S = the last 300) instead of the 100 ms gating hop, so either window lands
  within 10 ms of any event. EBU Tech 3341's file-based cases 10 and 13 slide a 3 s / 400 ms tone
  in 150 ms / 20 ms steps and expect the maximum to read the tone ±0.1 LU at every offset; a 400 ms
  burst 40 ms off the hop grid read 0.45 LU low, and case 13 failed at 16 of its 20 offsets. **M and
  S readings on transients change** by up to that much; the integrated measure is untouched in
  substance — every tenth sub-hop closes the same 400 ms block at the same 100 ms hop, LRA keeps its
  1 s cadence — and reads as before to the suite's tolerances (not bit-for-bit: 40 partial sums
  where there were 4). The sub-hop ring is 300 doubles; `process()` still allocates nothing.
- **fix(analysis):** `LoudnessMeter` sizes its block store by hops at the prepared rate rather than
  by seconds (a hop is 10 × lround (0.01·fs) samples — 100 ms only where fs is a multiple of 100),
  and blocks that arrive past `maxDurationSec` are counted in the new `droppedBlocks()` accessor
  instead of vanishing under a straight-faced reading. Additive API; nothing throws on the audio
  thread. A caller that must not lose a block sizes `prepare()` for its longest program and checks
  it reads 0.
- **test(analysis):** `felitronics_loudness_conformance_tests` — EBU Tech 3341 (2023) Table 1
  synthesized from the spec's text: cases 1–5 at 48 and 44.1 kHz (I, M and S on cases 1–2), case 6
  (the 5.0 channel weights), cases 9 and 12 (S and M settle on a periodic program), cases 10 and 13
  (S and M at every offset); then what Table 1 lets a wrong meter get away with, each pinned by a
  signal only the property under test can move — the absolute gate isolated from the relative one
  and at its boundary, the relative gate at −10 LU, every channel as power (in phase or not), the K
  shape at both rates against the published 48 kHz coefficients, a burst that tells 75 % block
  overlap from none, chunk invariance, and `process()` allocating nothing. 59 checks. Grown from
  the Looper Cat suite that gated the product's move onto this meter; crew-reviewed (Codex,
  DeepSeek, Antigravity, a fresh Opus mutating the meter).

## v0.22.2 — fill and peak keep the same bins (`analysis`)

- **fix(analysis):** `MultiResSpectrumPane` shares the classic pane's −200 dB floor. Its fill
  counted every positive bin while its peak-hold was floored at −120 dB, so a band of quiet bins
  read higher on the fill than on the peak trace and the two crossed near the plot's right edge. A
  −130 dB tone now reads −130 untilted — below any plot bottom, not flattened — and the peak trace
  never sits below the fill.

## v0.22.1 — the floor is silence, not a shelf (`analysis`)

- **fix(analysis):** both spectrum panes clamped a reading at −120 dB and then added the display
  tilt, so silence came out as a straight line rising at the tilt's slope — with +6 dB/oct it stood
  at −96 dB at 16 kHz, in plain view on a 120 dB range. `MultiResSpectrumPane` now applies the tilt
  to the power and floors afterwards (`readDb (f, fs, tilt, pivot)`); `SpectrumPane` drops its
  internal floor to −200 dB, deep below any plot bottom, so a tilted floor can never surface and a
  bin with real energy at −130 dB keeps it. Both panes freeze the tilt past Nyquist, where the
  columns repeat the last band that fits. Consumers that pinned the classic pane's −120 floor in
  their own tests move to `SpectrumPane::kFloorDb`.

## v0.22.0 — the analyzer reads constant-Q, and rides SIMD (`analysis`, `fftpffft`)

- **feat(analysis):** `MultiResSpectrumPane` — a constant-Q analyzer from several FFT lengths at
  once. One frame from the rolling tap feeds 16384 / 4096 / 1024-point Hann suffixes that share the
  frame's end (the short tier reports a transient first); a reading is the power in a 1/24-octave
  band, integrated over the tier's bins with fractional edges from double prefix sums and normalised
  by the window's measured ENBW — the one quantity two FFT lengths agree on for both a sine and
  noise, which is what makes a seam invisible. Tiers are used where their bin is at least two per
  band (seams 811 / 3246 Hz at 48 kHz), crossfaded in power over a third of an octave; below the
  longest tier's bin the lows are bin-limited and say so. A tier shorter than the frame hop
  Welch-averages as many half-overlapped windows as reach back over it, so a click between two
  frames is never missed; the per-bin smoothing runs on power, not dB, because a log-domain average
  carries a bias that depends on how the bin was fed. Design and physics: `docs/ANALYZER-MULTIRES.md`.
  Crew-reviewed twice (codex, deepseek, Fable — the latter with simulations); property tests, not
  golden files (`felitronics_multires_spectrum_tests`, 166 checks).
- **feat(analysis):** `RollingSpectrumTap::tryPull (dst, order, hop)` reports the samples that
  entered the ring since the previous publish — the hop that happened, which a block boundary or a
  missed UI tick stretches past the one requested. The two-argument form stays.
- **feat(analysis):** `SpectrumPane` becomes `SpectrumPaneT<Fft>` (constrained to the packed-Hermitian
  layout its bin loop reads) with `SpectrumPane` the scalar alias, so every consumer reads as before;
  it gains `reset()` (the next ingest seeds) for a consumer that returns to it after drawing another
  pane. `MultiResSpectrumPaneT<MaxOrder, MaxTiers, Fft>` likewise.
- **feat(fftpffft):** `PffftOrderedRealFft` — pffft's real transform in canonical order, which is
  exactly the packed-Hermitian layout the scalar reference writes (F(0) and F(N/2) in the first slot,
  then interleaved Re/Im, the e^{−jωn} sign), so it advertises `kPackedHermitianSpectrum` and is
  admissible wherever bins are read. Beside the z-order `PffftRealFft`, not instead: the convolvers
  keep their vectorised MAC; the analyzers get their SIMD. Nulled float for float against the scalar
  (DC / Nyquist exact, Re / Im ≤ 2.4e-7 of full scale) with basis vectors pinning every slot, and both
  panes on both backends within 0.02 dB. Per tick on an M-series Mac: a classic 16384 pane 267 → 73 µs,
  the multi-res pane 381 → 126 µs.
- **test:** `felitronics_spectrum_pane_perf_tests` — the panes' cost per UI tick on the scalar FFT,
  printed for the record with loose ceilings; the pffft suite gains the scalar-vs-SIMD comparison.

## v0.21.2 — the pack's input trims become a switch

- **feat(rigplayer):** the per-file `input_db` trims (a linked setting plays its neighbour softer)
  can be switched off — `setInputTrims(bool)`, ON by default. OFF feeds every capture at unity:
  for a library shot at one honest level the stated attenuations only push a capture's drive
  around, and into a nonlinear model a few dB less in is a lot less out. Read on the audio side,
  riding the existing slot ramps — the toggle lands on the next block, click-free.

## v0.21.1 — a stale landing dies with its pack; a failed load is refused

- **fix(rigplayer):** `unload()` right after `deliver()`, with no audio block between: the published
  landing no longer lands a stale model in the wiped law over an emptied stage — it dies before the
  forget is posted.
- **fix(nam):** a load that failed is refused (`BlendState::refused`), not asked for again on every
  block: a file that fails identically every time cost a fetch and a parse per service tick, for
  ever. The refusal lifts when the request names something else for the slot, or on any landing; the
  slot keeps its old capture and counts as at rest, so it may sleep.

## v0.21.0 — a slot at rest goes cold (`rigplayer`, `nam`)

- **feat(rigplayer):** a slot that has stood silent — weight exactly zero, request unchanged — for
  `RigPlayer::kColdAfterSeconds` (2 s) goes cold: its model is not run, not mixed, and stays loaded;
  its delay line is cleared as it falls asleep. The next change of request wakes it warm-up first, by
  the landing path, with the model's own field — no load, nothing unfed heard. A slot with any weight
  never sleeps: between two captures both models run, on one they do not. `setColdAfterSeconds()`
  (zero or less = never), `slotCold()`, `coldBlocks()`. Measured in OrbitAmp's block with the dial
  at rest: two passes for one sound, 4.5 % of a P-core and 14.3 % of an E-core — now one pass.
- **feat(nam):** BlendLaw owns the flag — `BlendPolicy::coldAfterSamples` (default never),
  `BlendState::cold` and `still`, `blendSameRequest()`. A cold slot counts as unfed; the wake is
  `blendLanded()` with the model and the need already held; the rest is counted as played; a slot
  asked for nothing is at rest too.
- **fix(rigplayer):** the models run on the planes they are given, not on the width the player was
  prepared for: a host prepared for stereo that plays one plane no longer pushes the silent second
  plane through both networks (8.8 % → 4.5 % of real time in OrbitAmp's block).

## v0.20.0 — the pack player comes home (`rigplayer`)

- **feat(rigplayer):** a new header-only module, `felitronics::rigplayer` — one device of a
  `.orbitrig` pack, playing: which captures sound for a panel (`namz::rig`'s policy joined to
  `pickBlend` along the gain dial), the crossfade between them (`nam`'s BlendLaw, once per block on
  the audio thread), the tone knobs as the pack describes them (sections → biquads, a curve → one
  minimum-phase FIR per side), a blend knob's dry path, the models' alignment from the pack's
  `lag_samples`. JUCE-free, no thread of its own: a load is a job the host runs anywhere but the audio
  thread (`takeLoadJob` / `run` / `deliver`). Moved, not copied, from OrbitCapture NAM with its tests,
  so the capture app and a plugin play through ONE player; `modules/rigplayer/README.md` is the
  host's contract. Self-gates on `FELITRONICS_WITH_NAM`.
- **build(nam):** the module's own namz pin rises from v1.1.1 to v3.1.0 (`namz_rig.h` with `tone`
  and `lag_samples`).
- **fix(nam):** BlendLaw's exact-zero swap test is spelled `<= 0.0` on a weight clamped to [0, 1] —
  the same zero, and GCC's `-Wfloat-equal` in the header-hygiene gate is satisfied.

## v0.19.0 — a model is built apart from its stage (`nam`)

- **feat(nam):** `NamStage::prepareModel(bytes, sampleRate, maxBlock)` is the heavy half of a load —
  the bytes parsed, both instances built, the rate-match and the prewarm prepared — done on any
  thread but the audio one and touching no stage; `install(prepared)` is the light half, a pointer
  swap on the message thread through the same pending machinery a load uses. A pack player crossing
  a capture paid up to 20 ms of its drawing thread per landing (measured in OrbitCapture NAM,
  2026-08-30); now it hands that work to a worker and installs the result. `loadModelFromMemory` is
  the two calls in a row — one path, two entry points — and behaves as before. A model prepared for
  other numbers than the stage runs at is prepared again inside `install`, at the old cost; the rate
  contract is judged at `install`, since only a stage can. `NamBackend` binds the normalize flag when
  it meets its stage (`bindNormalize`) instead of at construction.

## v0.18.0 — a WAV can be handed over as bytes (`io`)

- **feat(io):** `writeWavMemory` returns the encoded WAV image instead of writing it to a path. A
  library that keeps its audio in a database needs the image itself, and a temp file on the way there
  is a file to lose. This is the encoder that was already inside `writeWav`, lifted out unchanged:
  `writeWav` is now that call plus one `fwrite`, so every rule holds for both by construction — an
  empty result means the request was refused rather than a truncated file written, ragged channels
  are guarded against reading past an end, and a header RIFF cannot represent (more than 65535
  channels, a sample rate that is not finite, data past 4 GB) is refused outright rather than
  silently wrapped into a lying field. A test asserts the file and the memory image are byte for byte
  identical, so the two cannot drift apart later.

## v0.17.0 — one law owns which capture is audible (`nam`)

- **feat(nam):** `BlendLaw.h` — two model slots play the same input, a knob between two captures asks
  for a fraction of each, and this header is the single writer of that number. It comes from the
  OrbitCapture NAM player, where three places wrote it — the blend, the code following the nearest
  capture, and a warm-up gate bolted on later. Each was defensible alone; together they replaced
  models under a live gain and swung the weight across its whole range inside single 10 ms blocks
  (measured on one sweep of a nine-capture device: 62 loads, 401 parked retries, a full 1.000 swing).
  Seven attempts to fix that by ear in one day, three of which made it worse — hence a law with
  properties provable by construction: the weight moves at most one step per block whatever happens,
  a model is replaced only in a slot whose weight is exactly zero, and progress needs no timeout
  because some slot's goal is always a rail. What gives is instantaneous accuracy — the incoming
  model arrives ~132 ms late and the outgoing one carries alone until it does, which reads as a knob
  rather than as a fault.
- **feat(nam):** the mixing arithmetic lives with the law rather than in the host: `blendMix` and
  `blendDelay`, with the properties that make a slip arithmetic instead of a listening exercise —
  weight 0 is the first slot alone, 1 the second, halfway their linear average (two coherent halves
  stay at unity, not +3 dB), a ramp lands exactly on the weight its block ends with, and a delay line
  carries across a block boundary with nothing lost. The host had lost one of its two model calls, so
  the mix was raw DI against a model; that is audible only as loudness rippling with the knob, and a
  law about weights cannot see it.
- **feat(nam):** `NamStage` reports the receptive field a model must be fed before it means anything,
  and the law counts it BEFORE the block rather than after — one block of conservatism instead of a
  few hundred samples of a network that has not heard the last 132 ms.
- **fix(nam):** an unfed slot stays silent even when that means silence. The old "something must
  sound" rule named slot 0, and a `NamStage` with no model is a PASSTHROUGH — so a device change put
  the raw DI on the output, some ten decibels above any normalised capture. The step now carries its
  own gain, ramped by the same clamp as the weight.
- **refactor(nam):** `blendDelay` takes a pointer and a capacity instead of a template over the
  tail's array size, which no caller with a `std::array` or a runtime-sized buffer could satisfy
  without copying the function — which is exactly what the host had done, leaving 45 assertions
  guarding code that never ran.

## v0.16.0 — the plot's own maths comes home, and 6 dB/oct stops diving (`analysis`, `eq`)

- **feat(analysis):** `PlotMap.h` + `SpectrumPane.h` graduate from TabbyEQ's `eqview` incubator.
  `PlotMap` is the log-frequency / dB coordinate map of an EQ plot (freq↔x, dB↔y, both ways);
  `SpectrumPane` is the analyzer pipeline — Hann → real FFT → smoothed dB with peak-hold, resolved
  into liquid log-frequency columns. Header-only, JUCE-free, no plugin includes. A second consumer
  is what unlocked the move: OrbitAmp's EQ links need the same two, and TabbyEQ now reads them from
  here like everyone else.
- **fix(eq):** `matched::lowpass1` no longer dives at Nyquist. Bilinear puts a zero at z = −1, so
  every 6 dB/oct low pass fell to −inf in the top octaves instead of rolling off gently — and since
  the first-order section is the odd member of every variable-slope cascade, orders 6 / 18 / 30 / 42
  inherited it. The replacement is a magnitude-matched one-pole: exact at DC, at f0 (−3.01 dB) and
  at Nyquist (the analog value), its pole taken from the quadratic whose roots are reciprocal, so
  the small root is always the stable one. `highpass1` keeps bilinear — its zero belongs at DC,
  which is where the analog filter has one.
- **test(eq):** the first-order sections are pinned against the analog prototype
  `|H| = 1/sqrt(1 + (f/f0)²)` over 44.1 / 48 / 96 / 192 kHz × 20 Hz … 20 kHz: unity at DC, −3.01 dB
  at f0 (4e−08), the analog magnitude **at Nyquist** (7e−15 — a bilinear section reads −inf and
  fails outright), monotone descent with no ripple between the matched points, and ≤ 0.97 dB from
  the prototype in between.

## v0.15.0 — a magnitude table becomes a filter (`lineareq`)

- **feat(lineareq):** `MagnitudeCurve.h` — dB against a log-frequency grid turned into a runnable
  minimum-phase FIR. A measurement produces a table, a pack ships a table, and every reader of one
  needs exactly this; it lived in OrbitCapture NAM until a second consumer appeared, which is one
  copy earlier than the house rule allows. `logFreqGrid`, `curveDbAt`, `heldOutsideBand`,
  `magnitudeCurveToFir`, `curveAtPosition`.
- Minimum phase via `MixedPhaseFir` at k = 1: the magnitude is exact, there is no pre-ring and no
  bulk delay, so switching a control does not shift the audio in time. It is not a fitted 1-pole —
  measured on real hardware the best fit missed a Big Muff tone control by 15 to 65 dB across its
  30 dB of tilt, so the curve is the only honest description.
- `curveAtPosition` interpolates against the producer's stated `norm`, never the array index: a
  control measured at 0/30/150/300 on a 300-degree dial has norms 0, 0.1, 0.5, 1.0, and treating
  those as evenly spaced puts a knob halfway four to six decibels wrong.
- **No product policy travels with it.** Trust thresholds and what a pack may claim stay with the
  format that owns them; the trusted band arrives here as two numbers.

## v0.14.0 — the auto-first dynamic EQ across placement lanes (`dynamics`, `eq`, `dynamiceq`)

- **feat(dynamics):** `RelativeLevel` — a slow programme-level estimator with an offset, so a
  threshold can be expressed relative to what a band's own region normally does rather than in
  absolute dBFS (meaningless when the same setting is shared by lanes sitting 20 dB apart).
  `BandBallistics` derives attack/release from a band's own `fc` and `Q`: a band's envelope rises
  with its inverse bandwidth `Q/fc`, not with its period, so one 0..1 "deviation" pair means the
  same thing on a 60 Hz band and a 7 kHz one. Both JUCE-free, software-flushed, control-rate safe.
- **feat(eq):** `BandParams` gains point-level `DynParams { on, rangeDb, thrDb, thrAuto, atk, rel }`
  and `EqBand` gains the **delta seam** `setLaneDeltaDb(Lane, dB)` — the band accepts a NUMBER and
  applies it with its own `Svf` inside the lane, after the matched static biquad and before the M/S
  fold. `eq` therefore takes no dependency on `dynamics`. With `dyn.on` false a band is
  bit-identical to one built before dynamics existed.
- **feat(eq):** `EqEngine::captureSectionInput()` + `bandAt()`/`bandCount()`. A dynamics layer must
  detect on the section's own input: in a series chain a band's input is the previous bands' OUTPUT,
  so their moving deltas would modulate later detectors at overlapping frequencies and the chain
  would pump. The engine deliberately does not hard-code the interleaved loop — it hands out the
  bands so a composition layer can run "delta, then band" in chain order.
- **feat(dynamiceq):** `LaneDynamics` — the lanes-aware composition layer: per-lane sidechain probe,
  envelope follower, `RelativeLevel`, `GainComputer` (ratio and knee fixed internally), GR
  ballistics and the seam write, at a K=16 control rate. Reuses the primitives rather than
  `DynamicEqBand`, whose single-SVF topology would forfeit the matched static curve.
- **fix:** the hardening rounds behind the above — a per-sample control path (interval-correct
  coefficients; a block-size-invariance test pins it), gate parity, no live state carried across
  `reset()`, per-lane state drop on disengage, no seed-from-warmup duck at transport start, a
  hot-signal clamp, and `core::fastGainToDb` on the detector paths (max error 0.0001 dB) to keep the
  per-sample path affordable.
- **fix(dynamiceq):** `laneSignal`'s `switch` names the Stereo lane instead of leaning on `default:`
  — `-Wswitch-enum` is part of JUCE's recommended warning set, so the old form fired in every
  consumer that included the header.

## v0.13.1 — nlohmann rides with the exported namz headers (`felitronics::nam`)

- **fix(nam):** namz's rig headers (`namz_rig*.h`) publicly include `<nlohmann/json.hpp>`;
  consumers reaching them through `felitronics::nam` (OrbitCab's v0.13.0 migration) could not
  resolve the include. The module now exports namz + nlohmann include dirs as a usage
  requirement of its shipped headers — NAM/Eigen stay PRIVATE.

## v0.13.0 — CabConvolver + opt-in NAM backend (`felitronics::convolution`, `felitronics::nam`)

- **feat(convolution):** `CabConvolver` — the product-level cab IR wrapper moves in from OrbitCab,
  DSP byte-identical: reference-unity RMS normalization (2 kHz-shaped reference, ±30 dB clamp,
  −60 dB near-silence floor), fixed NUPC schedule (head 128, 50 ms click-free crossfade),
  mono-broadcast / LRDiag true-stereo publishing, latest-wins retry of a load rejected
  mid-crossfade, staged-taps accessor for offline blend analysis. The
  `FELITRONICS_WITH_PFFFT` backend selection now propagates uniformly through the module's
  INTERFACE (fftpffft registers before convolution; hygiene/test targets deduped).
- **feat(nam):** `felitronics::nam` — NEW opt-in compiled module (`FELITRONICS_WITH_NAM`, default
  OFF, CMake ≥ 3.24). `NamStage` is OrbitCab's `cab::AmpStage` ported verbatim onto
  `felitronics::neural::NeuralStage` + `felitronics::core::StreamResampler`: dual-instance true
  stereo, −18 dB loudness makeup with per-model trim, model-rate contract, bounded retire queue
  with a lossless last-wins pending intent. Loads raw `.nam` and packed `.namz` (namz used
  directly on std buffers; NAMZ_IMPLEMENTATION is compiled exactly ONCE here — consumers include
  `<namz.h>` without the define and link this archive). NeuralAmpModelerCore pinned by SHA and
  linked WHOLE_ARCHIVE (architecture self-registration) with PRIVATE usage requirements — no
  Eigen/nlohmann/NAM headers leak into consumer TUs; namz pinned by immutable commit; both
  overridable by fail-loud local-source cache vars. ctest: analytic Linear FIR, namz round-trip +
  v1 wire compat, loudness/trim makeup, rate-contract refusal against a live model, 70+
  frozen-audio swap stress (mirror-publication discipline, deferred clear), true-stereo
  independence, exact 96 kHz latency pin, resampled/truncated IR loads, and a dedicated
  `EIGEN_RUNTIME_NO_MALLOC` gate (Linear + WaveNet allocation-free; LSTM/ConvNet documented as
  upstream-allocating at the pin, load-compatible by design). Consumed by OrbitCab (cab::AmpStage
  / cab::Convolver become aliases) and OrbitCapture NAM's Queue audition player.
- **test(support):** `approx()` NaN-proofed (unordered comparisons now FAIL); CI adds NAM rows on
  all three OSes (first MSVC nam_core build) + the sanitizer job, and keeps default OFF/OFF rows.

## v0.12.0 — decoupled-hop rolling analyzer tap (`felitronics::analysis`)

- **feat(analysis):** `RollingSpectrumTap` — a lock-free SPSC analyzer tap whose snapshot cadence
  (hop) is decoupled from its analysis window size. One rolling ring (max order 14 = 16384) serves any
  FFT size ≤ MaxOrder from a single buffer: `publishIfDue(order, hop)` copies the most-recent
  `N = 1<<order` samples — chronologically, across the ring wrap — into a single-slot immutable mailbox
  and force-publishes when the order changes, so a consumer can offer a **selectable analyzer FFT size
  with a click-free live switch** at a steady UI-rate cadence (overlapping windows for large N, gapped
  for small N) — no per-order ring duplication, one write per sample. `tryPull(dst, outOrder)` reports
  the order the frame was captured at so a wrong-size frame is discarded across a switch. `reset()`
  restarts the producer only and never revokes a mid-pull reader — closing a torn-frame race that a
  reset-clears-ready design would have when `prepareToPlay` runs against a live GUI reader. Tear-free by
  the same acquire/release ownership handoff `SpectrumTap` uses, plus per-frame order/size metadata.
  Header-only, JUCE-free. ctest: variable-N snapshot, warmup gate, hop cadence, forced publish on order
  change, chronological wrap copy, race-free reset. Consumed by tabby-eq's analyzer-resolution feature.

## v0.11.0 — UTF-8 → ASCII romanization for filename slugs (`felitronics::text`)

- **feat(text):** `felitronics::text` — JUCE-free UTF-8 → ASCII romanization for filename slugs.
  `decodeUtf8` (malformed bytes skipped, never guessed), `romanize` (one code point → ASCII: German
  umlauts / ß / Œ, common Latin-1 accents, and full Russian Cyrillic incl. the multi-letter cases
  Ж→Zh, Х→Kh, Щ→Shch, Ю→Yu; hard/soft signs vanish), and `transliterate` (the whole-string fold).
  Extracted from OrbitCapture NAM's `Matrix.h` so every Darwin's Cat product shares ONE romanization
  table — the RAW text stays in metadata; this is for FILENAMES, where cross-platform sync (macOS NFD
  ↔ Windows/Linux NFC) and legacy ASCII parsers demand plain ASCII and a never-empty result.
  Header-only, zero deps, header-hygiene clean under the strict downstream warning set. 41
  falsification checks (every mapping asserted; intentional collisions like е/э→"e" pinned so callers
  detect them).

## v0.10.0 — file I/O module + more OrbitCapture/OrbitCab promotions

Each promotion behind the extraction bar (theory-first falsification tests + adversarial crew
codex/deepseek + NULL where possible), landed as OrbitCapture NAM (the second capture product) and
OrbitCab consume them.

- **feat(io):** new zero-dependency `felitronics::io` module — minimal self-contained WAV
  read/write (`readWav`/`readWavMemory`/`writeWav`/`writeWavMonoF32`), moved from OrbitCapture's
  `oc/wav.hpp`. Crew-hardened: corrupt/truncated chunks are rejected loudly (never clamped),
  `WAVE_FORMAT_EXTENSIBLE` requires its full body + SubFormat GUID, checked chunk advance (no
  32-bit wrap), writer refuses headers it cannot represent (u16/u32 overflow, non-finite rates),
  data must be frame-aligned. WavTests pin every case.
- **feat(measurement):** `PeakClip.h` — standalone `scanPeakClip` (peak dBFS + flat-top clip run)
  extracted from `gateRecording` for sweepless consumers (NAM reamp takes); a NaN breaks a run,
  `clipRunSamples` clamped to >= 1. `detectLeadingSilence` — OrbitCab's forward-scan head-trim onset.
- **feat(core,saturation):** NaN/Inf poison guards + chunking hardening (promoted from OrbitCab),
  with theory-first falsification suites.

## v0.9.0 — RT stream types (`core::RtStreams`) + model guess + the mix-view overlay facade

Three promotions from OrbitCapture, each behind the extraction bar (tests + adversarial crew
codex/deepseek + NULL where possible) — landed as a second capture product arrived to consume them.

- **feat(core):** `RtStreams` — the RT buffer-swap discipline AS TYPES (`AuditionStream`,
  `ConvStream`, `RecStream`): flag DOWN → mutate within fixed/reserved storage → flag UP; std +
  atomics only. Crew-hardened with four pinned fixes: a zero-length publish no longer arms an
  empty looping clip (reader OOB), `RecStream::push` release-publishes `len` (an ARM harvest could
  read an unsynchronized sample), the audition buffer is fixed-size (no reliance on
  assign-within-capacity), `reserve()` guards int lengths. Reader rules (ACQUIRE gates, per-block
  flag reload) and the known one-block quiescence gap (+ the deferred reader-ack epoch design) are
  documented in the header. Suite includes a discipline-respecting concurrent smoke; ASan/UBSan
  clean.
- **feat(measurement):** `ModelGuess` — gear-model detection in a free-form file name against a
  catalog. Conservative contract: exact fingerprint (the entry's last token) beats the 3+-digit
  bare-number fallback; ANY ambiguity → no guess (a wrong guess poisons imported metadata
  forever). Crew fix: locale-FREE tokenizer (`std::isalnum` could admit high-bit bytes under a
  single-byte locale).
- **feat(blend):** `Overlay`/`makeOverlay` — the one-call mix-view facade: per-mic curves (post
  filters + gain), the blend curve (post Master), and the interference column (PRE-Master basis so
  a Master rolloff never reads as phase cancellation, then faded by the Master's |H| — a display
  weight by deliberate product decision, documented). NULLed to the hand-composed primitives at
  machine epsilon; analytical pins (+3.01 dB in-phase twins, deep 180° notch, |H| only
  attenuates). Crew fixes: bit-parity gain expression, no partial overlays on degenerate input,
  non-finite params heal to defaults (the offline convention).

## v0.8.0 — mic-blend engine (`felitronics::blend`) + fine IR alignment (`measurement::XcorrAlign`)

The MIX side of the IR-capture family joins the capture side (v0.7.0) in core — both extracted from
OrbitCapture (written portable by design), so OrbitCab and other consumers share ONE numerical
fingerprint for blending and aligning multi-mic IRs.

- **feat(blend):** a new `felitronics::blend` module — the **offline multi-mic IR blend engine**:
  `StripParams`/`MasterParams` (per-mic gain / phase / fractional time shift / HPF / LPF + master,
  solo/mute audibility rules), per-section-Q Butterworth HPF/LPF, Hilbert-based phase rotation,
  windowed-sinc fractional shift (positive `shiftMs` = delay), and `blendIrs` (weighted sum +
  master chain). Canonical defaults live HERE (80 Hz/24 dB · 8 kHz/12 dB) — consumers must not
  re-declare them (a drifted default silently re-voices saved mixes). Extraction was gated by a
  **byte-NULL** against the app's previous in-tree engine; OrbitCapture's `ocap::` blend names are
  now `using`-shims over this module.
- **feat(measurement):** `XcorrAlign` — **fine time/polarity alignment** of an IR against a
  reference by normalized cross-correlation (`xcorrAlign` / `xcorrAlignSet`, ±maxLag samples,
  fractional result): per-lag normalization with a Cauchy–Schwarz corr ≤ 1 bound, an onset-delta
  guard that REFUSES (corr = 0) when the two onsets sit further apart than the search range,
  polarity from the best normalized lag, subnormal-safe denominators. corr = 0 is the "no
  confident suggestion" contract — callers must leave such a channel untouched.
- **robustness (crew-hardened):** the adversarial consilium (codex + deepseek) hit XcorrAlign
  before merge and found real bugs: a fixed-denominator normalization that could rank a wrong lag
  above the true one AND report corr > 1; a confident wrong lag + false invert when the true lag
  lies beyond `maxLag` (repro'd at corr 0.91); `inf` via `√(eR·eC)` underflow on subnormal
  energies; a window off-by-one (+ the analysis-window floor raised 16 → 64). Each fix carries a
  pinned counterexample test, and the whole search is NULLed against a brute-force oracle over
  randomized signals.

## v0.7.0 — offline IR measurement + display curves (`felitronics::measurement`, `felitronics::analysis::offline`)

The **offline** (message-thread, allocating, double-precision) half of an IR-capture pipeline, extracted
from OrbitCapture (written portable by design) so the capture math lives in core instead of the app. This is
deliberately **NOT** the RT path — real-time consumers still use the float `felitronics::core::fft` seam;
these transform whole ~6 s captures where the numerical floor must sit far below the analog chain's, so they
run in `double`. Every function clamps/heals non-finite params; correctness is oracle- + **numpy-cross-NULL**-
anchored (an independent `numpy`/direct-time-domain recompute nulls the C++ to machine epsilon) and
ASan/UBSan-clean.

- **feat(measurement):** a new `felitronics::measurement` module — an **exponential sine sweep (ESS / Farina)**
  generator + matched inverse (`Sweep`), **Farina deconvolution** with a latency-absorbing onset search past
  the harmonic region (`Deconvolve`), IR post (onset / trim / peak-normalize, `IrPost`), a pre-deconv
  **capture-quality gate** (clip / non-finite / sweep-presence / SNR, `CaptureGate`), and **multi-mic
  common-onset alignment** that preserves the inter-mic comb (`MicSetAlign`). NULL-verified: sweep⊛inverse≈δ,
  known-answer in-band magnitude, and the direct-time-domain convolution vs `numpy.convolve`.
- **feat(analysis):** `felitronics::analysis::offline` display curves — `logMagnitudeCurve` (a 1/N-octave
  RMS-power-smoothed magnitude on a log-f grid, in dB; energy-preserving `10·log10(mean|X|²) = 20·log10(rms)`)
  and `interferenceDb` (where a multi-mic blend cancels or reinforces vs an incoherent power sum). Namespaced
  `::offline` to keep the RT-metering contract of the rest of `felitronics::analysis` intact.
- **refactor(core):** the shared offline double FFT (`nextPow2` / `detail::fftInplace` / `convolve` /
  `magSpectrum`) is **promoted** `measurement` → `felitronics::core::offline` (`core/OfflineFft.h`) now that a
  second offline-FFT consumer (the display curves) exists. `measurement/Convolve.h` re-exports it — the
  measurement API and all its tests are unchanged.
- **robustness (crew-hardened):** an adversarial "break-it" consilium (deepseek + antigravity + Fable) found
  9 edge-case bugs in `measurement` and 4 in the display curves (a tiny/inf smoothing band → `(int)ceil(inf)`
  UB; a non-pow2 `minNfft` → binHz↔bins skew; a huge-finite sample → `ΣM²` overflow → NaN; a top-bin band
  inversion) — each fixed **with a regression test**. Verified false alarms were rejected against the code.

## v0.6.0 — noise gate (`felitronics::dynamics::NoiseGate`)

A new `felitronics::dynamics` primitive: a dual-detection (ISP Decimator "G-String" style) **noise gate** —
architecturally distinct from the continuous `Compressor` (a bistable **Schmitt trigger** with hysteresis +
hold, a **LINEAR-fast open / EXP-slow close** VCA, a closed floor, and an on/off **enable crossfade**). It
**composes the module kit** (`ChannelLinker` linked key + a Peak `EnvelopeFollower`) and adds the gate-specific
state machine on top. Extracted from OrbitCab's in-amp gate (written portable by design) so plugins don't
reinvent it.

- **feat(dynamics):** `NoiseGate` — a two-phase **keyed** API (`analyse()` fills a per-sample gain curve from
  the clean KEY; `applyGain()` attenuates a possibly-different downstream buffer, so any latency between the two
  is free lookahead) plus a self-keyed `process()` convenience. A `Config` voicing struct (defaults = OrbitCab's
  shipped tuning), `seedEnabled()` for a restored on-state, `currentGain()`/`currentCoreGain()` for a GR meter;
  zero latency.
- **RT / robustness:** NaN/Inf-safe (detector-input clamp + non-finite heal) and denormal-safe in software
  (Law 8); no-allocation-in-`process()` proven; block-split NULL + sample-rate-invariance + adversarial
  fail-open / transient / low-note-chatter tests.

## v0.5.0 — the non-uniform (Gardner) convolver: block-independent, cheaper at small buffers

The v0.4.0 follow-up is delivered. A **non-uniform partitioned (Gardner 1995) convolver** — a time-domain head +
geometrically growing overlap-save FFT stages — replaces the fixed-`P=128` path. It is **block-INDEPENDENT**,
**true sample-zero-latency**, and **cheaper than `juce::dsp::Convolution` across the small, low-latency buffers a
live rig runs** (JUCE only wins the mean at large power-of-two blocks — a theorem of true zero-latency, not a
shortfall). Flat mean **~0.6 %RT on an Apple M5 Pro / ~1.08 % on an Intel i9-13900H** at every DAW buffer from 16
to 4096, from one `prepare()`.

- **feat(convolution):** `NonUniformConvolver<Fft>` — the mono zero-latency NUPC primitive (head `P0` + capped
  octave-doubling to `B_max`, its own frequency-domain delay line per stage).
- **feat(convolution):** `MatrixConvolverNupc<Fft>` — the shipping 2×2 matrix convolver on one raw-L/R history:
  all four routings (mono / LRDiag / MSDiag / Full) + a **click-free 2-slot smoothstep crossfade** for live IR
  swaps. A NULL-verified drop-in for `MatrixConvolver`.
- **feat(lineareq):** the linear- & mixed-phase EQ now convolves on `MatrixConvolverNupc` (A/B transparent — the
  change is CPU + zero latency only).
- **fix(convolution):** removed the long cold-prime crossfade — a cold FDL already yields the exact causal
  convolution, so the ~2.7 s first-activation fade only attenuated correct output and never completed for renders
  shorter than it (offline/short renders came out ~10 dB down). Every swap now uses the short anti-click fade
  (in `NonUniformConvolver`, `MatrixConvolver`, `ConvolutionEngine`).
- **fix(convolution):** `setIr()` broadcasts a mono IR to both channels on a stereo instance (was rejected); a
  `static_assert` pins `state_` lock-free; `maxIrSamples` capped against a stage-offset overflow.
- **docs/tools:** `PERF-NUPC-VS-JUCE.md` (a two-machine 4→8192 fine log-ladder sweep — Apple M5 Pro + Intel i9-13900H,
  adaptive 3–10 warmed reps — rendered to an in-repo SVG chart), `PERF-CONVOLVER-JUCE-GAP.md` (the design ADR);
  `fftbench` head-to-head + `FCORE_FINE_SWEEP`; `tools/plot-convolver-sweep.py`.

Verified by an architecture consilium (per phase), a Popper falsification campaign (differential fuzzer vs
`PartitionedConvolver` + double precision, ASan/UBSan/TSan), and a release-candidate crew review. `ctest` green.

## v0.4.0 — the #1 performance debt resolved

The scalar-FFT / `O(P)` direct-head convolution bottleneck — long-convolution cost that **exploded at host
block 2048+** — is fixed. Cost is now **block-INDEPENDENT**, zero-latency, and JUCE-free.

- **feat(fftpffft):** a new optional, compiled SIMD FFT backend `felitronics::fftpffft::PffftRealFft`
  (`-DFELITRONICS_WITH_PFFFT=ON`, default OFF) — vendored 2-file pffft, hidden-visibility, cross-backend NULL
  parity tested on x86-64/SSE + arm64/NEON incl. ASan/UBSan. (#25)
- **feat(lineareq):** the audio FFT backend is now a template parameter; the design-time FFTs stay pinned to
  the scalar packed-Hermitian layout and the split is **compile-enforced** — a SIMD backend cannot silently
  corrupt a designed FIR. (#24)
- **feat(lineareq):** the convolver partition is **decoupled from the host block** (fixed internal `P=128`).
  A 131072-tap linear-phase EQ is ~2.0 %RT with pffft at every host block (was ~39 %RT @ block 8192 on the
  old scalar path). (#26)
- **feat(convolution):** SIMD-aligned every FFT-seam buffer (`core::fft::SeamAllocator<64>`). (#23)
- **docs:** `PERF-SCALAR-FFT-BOTTLENECK.md` marked RESOLVED.

Known follow-up (tracked separately): the fixed-`P=128` convolver is block-independent but 2–13× more CPU than
`juce::dsp::Convolution` at host blocks ≥ 512; a non-uniform (Gardner) partitioned convolver is planned to
beat it. See `docs/PERF-CONVOLVER-JUCE-GAP.md`.

Prior versions: see the `v0.1.x` – `v0.3.0` git tags.

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

The control's threshold was also wrong on its own terms: it demanded a 25 % collapse and failed at 8.635 →
8.225 dB, where the mechanism had worked perfectly. How far an ungated p95 falls is a fact about the
distribution's shape near its top — with half the programme silent it becomes the music's own p90, and a dense
distribution moves little. The check now claims the direction, which is what is true.

**The budget grew by one histogram**, and every allocation oracle in three suites went red at once — which is
the oracle doing its job. Two of them were also carrying numbers their own run no longer produced (a "727 696 B"
in a sentence beside an assertion parameterised by the histogram count); those now print what they measured. One
precondition pinned three literals, two of which were functions of things that test is not about; it now pins
the per-pass cost and spells the rest through the expressions the budget is made of.

ABI v10 — one new struct with its own id and row, one entry point, one field at the end of `fc_loudness_request`
(120 → 152 at v10). `fc_gr_stats` is nested by value and frozen with it; nothing existing moved.

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
holding render is unchanged render for render and bit for bit, and with the budget exhausted the idle render is still
what comes back — the one case where a quiet render is delivered on purpose, because the guarantee outranks the
loudness.
- `pairFor` is the one place a drive becomes a `(gain, ceiling)` pair: `d = g - c` and the two are clamped to ±60 dB
  one number at a time, so a ceiling picked for the true-peak aim alone put the gain past its clamp and the drive
  RENDERED was not the drive chosen. The ceiling is chosen for the drive, inside the window that keeps the gain in
  range and at or under the promise, less the between-grid overshoot already measured (a quarter of a decibel on a
  15 kHz tone, five times the aim's own margin); where that window is empty no pair expresses the drive and the
  rescue is not taken.
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

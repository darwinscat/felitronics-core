### analysis · tools — how far, how often and for how long a render goes over its ceiling

`analysis::PeakExcursions`: the measurement behind the choice between a clipper and a limiter. The caller
supplies a render made **without** the limiter and the ceiling it means to deliver at; this reports the
excursions that would have to be removed.

**On the reconstruction, not the sample grid**, and that is the first thing the design round settled. A dBTP
ceiling is a statement about the reconstructed signal. Two samples of 0.85 read −1.41 dBFS on the grid and
**+0.4224 dBTP** reconstructed — a grid finder reports *nothing* for an excursion of 1.4 dB. The under-read
of a crest at *f* is `20·log10(cos(π·f/fs))`: 0.000 dB at 60 Hz, −0.019 at 1 kHz, −2.011 at 10 kHz, −5.105
at 15 kHz. Bass crests survive the grid intact and transients do not, which is exactly backwards for this
decision. The 4×/32 topology is `ReferenceTruePeakMeter`'s, spelled and `static_assert`ed equal: the
oversampler's own default is 64 taps per phase, and taking it would have built an instrument that disagrees
with the delivery certificate with nothing in a passing suite to say so.

**The certificate's sample-peak floor is published, not faked.** `truePeakLinear()` is
`max(reconstructed, samplePeak)`, so a certificate can exceed the ceiling because of a raw sample while the
reconstruction never does. There is no interval over the ceiling in that case, so no run is invented — both
peaks are published separately and a caller can see which one carried the reading.

**What moves with the merge window and what does not.** `runCount`, the duration histogram and the
percentile are functions of a parameter meant to be swept; `occupancy` and the total dose are accumulated
over raw above-samples and are **merge-free**, so a sweep keeps two fixed points. Dose is per-sample linear
excess, which is the only version additive across programmes: "peak × duration once per run" over-reads a
parabolic crest by 3/2 and counts merged gaps as excess, and dB excess is not additive at all.

**The percentile never reads the run list.** It comes from a fixed-time duration histogram — 1/32 ms bins,
because oversampled bins would cover 2.7 ms at 768 kHz, under the longest class edge — so it, the class
counts and doses, the crest histogram and both ceiling densities stay exact after the list fills. Law 11's
"exhaustion is data", applied to statistics and not only to coordinates.

**A rule of the form "p90 ≤ 2 ms" is satisfied trivially by a programme with no excursions at all**, so
`runCount()` is published beside the percentile and must be tested with it. "Nothing went over" and "things
went over but briefly" are different findings and read differently.

#### The field the request asked for could not be measured

`lowShare` — "the fraction of run energy below 120 Hz" — is not computable as specified, and its natural
implementation inverts the answer. A run is 0.5–5 ms and one cycle at 120 Hz is 8.3 ms, but the killing
defect is **phase**, not resolution: a causal LR4 low-pass at 120 Hz lags 60 Hz by 86.63°, so at the input's
crest the filtered signal sits near its own zero. Measured on a **pure 60 Hz sine crest** the field reads
**1.4 %** — and a rule reading "clipper first when lowShare ≤ 25 %" would then reach for a clipper on bass,
which is the one thing a clipper must not touch. The quantity is monotone in phase, not in bass content.

In its place, **the crest's own frequency, from its shape**: for `A·cos(ωt)` over a ceiling by `e` for a
duration `D`, `e = A·ω²·D²/8`, so `ω = sqrt(8e/(A·D²))`. No filter, no phase, no window, gain-invariant. A
60 Hz sine reads 59.3 Hz, a 40 Hz sine 39.6, a two-sample click 10.8 kHz. It is published per run and as a
third-octave histogram of **dose**, so a consumer sweeps "is this bass" on its own corpus rather than
receiving a boundary baked in at one frequency.

It reads **low by a known amount**, because a cosine crest is deeper than the parabola estimating it:
exactly `sqrt(2e/A)/acos(T/A)` — 0.991 at one decibel over the ceiling, 0.955 at six. A 40 Hz tone one
decibel over reads 39.63 Hz and measures 39.6. Not corrected, because the correction needs the true crest
shape, which is the thing being estimated; documented with its table and pinned against that closed form
rather than against the tone with a loose tolerance.

**`ceilingDensity()` and `ceilingDensityAbove(minusDb)`** — the fraction of the reconstruction's local
maxima sitting within a window of the programme's own maximum, from a 0.05 dB histogram of maximum levels.
The tell for a source that was true-peak *limited* rather than clipped: its peaks cluster at the ceiling
with no flat tops for a clipping test to find. Two denominators, because the first is dominated by the
musical waveform and reads low on a dynamic programme for reasons that are not about peak control.

#### A defect found because one number did not match physics

`samples_` was incremented **after** `consumeOs`, whose bound is `samples_ · kFactor` — so the first chunk
was bounded against zero and stopped after 128 oversampled samples, swallowing input samples 32 to 1024 of
every programme. The aggregates looked healthy: 117 of 118 runs read the right crest frequency and the
damaged one carried 0.34 % of the dose. The tell was that one run reading **135.9 Hz for a 60 Hz tone** —
because a truncated crest is a faster crest. The suite now pins the count at exactly `2·f·T`, which makes
the same defect loud: 120 runs in one second of 60 Hz, not 118.

#### Also in this version — three follow-ups the consumer found on v0.44.0

Two of them were parameters **accepted and then ignored**, which is worse than a refusal because
`params()` echoed the caller's own number back while the measurement ran on something else:

- `skipBlocks = -1` compared as `blockIndex >= -1`, true of every block, so nothing was skipped.
- `dutyThresholdDb = -5` makes the linear gate `10^0.5 = 3.16`, which no band can reach against its own
  frame's maximum — every duty read 0 and the measurement was silently empty.

Both refused now. And `sideFractionBelow()` answered `0.0` for a frequency the crossover would refuse,
which is indistinguishable from a legitimate reading of a perfectly mono low end; it is **−1.0** now, and a
fraction cannot be negative. `LowEnd::settlingBlocks` is exported on the probe ABI, because the consumer was
about to re-implement `u = 10.233` on its side and a constant copied across a boundary is a second
definition waiting to drift.

**Surface:** `fc_probe_excursions_run`, `_run_with`, `_scalars` (22), `_runs` (6 per row), `_classes_out`,
`_crest`, `_ceiling_density`, `_ceiling_maxima`, `_storage_bytes`, and `fc_probe_lowend_settling_blocks`.
Both roads print `excursions v1`, and the JS formatter's width checks are the version gate.

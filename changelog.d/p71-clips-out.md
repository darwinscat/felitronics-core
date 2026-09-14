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

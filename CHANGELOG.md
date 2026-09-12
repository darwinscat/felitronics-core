<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

# Changelog

Notable changes to felitronics-core. Releases are git tags (`vX.Y.Z`); the project VERSION lives in
`CMakeLists.txt`.

## Unreleased

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
- **The suites' allocation counters install EVERY form of `operator new`,** the over-aligned one included.
  Without it `eq::EqEngine`'s 331 KiB — the largest single request a create makes on the default geometry;
  at 16 channels and an 8192-sample quantum the saturator's flat scratch is 8 MiB — is invisible, and a budget
  check would have compared two numbers that both left it out (the blindness P52 names). Pinned over 4 rates ×
  3 widths × 4 topologies, byte for byte, for `create` and for `configure`.

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

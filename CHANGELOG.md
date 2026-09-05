<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

# Changelog

Notable changes to felitronics-core. Releases are git tags (`vX.Y.Z`); the project VERSION lives in
`CMakeLists.txt`.

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

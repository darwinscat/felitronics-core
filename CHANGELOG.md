<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

# Changelog

Notable changes to felitronics-core. Releases are git tags (`vX.Y.Z`); the project VERSION lives in
`CMakeLists.txt`.

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

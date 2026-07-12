<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

# Changelog

Notable changes to felitronics-core. Releases are git tags (`vX.Y.Z`); the project VERSION lives in
`CMakeLists.txt`.

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

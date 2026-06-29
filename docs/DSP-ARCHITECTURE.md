<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->
<!-- Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa. -->

# felitronics-core — DSP architecture (living ADR)

**Status:** draft / in progress · **Owner:** Darwin's Cat (Felitronics line) · **Started:** 2026-06 ·
**Reviewed:** 2026-06-29 (panel — see §9).

`felitronics-core` is a **shared, framework-agnostic, JUCE-free DSP core** for the whole product
family. One set of battle-tested, real-time-safe DSP primitives that every product builds on:
plugins (TabbyEQ, a future compressor / true-peak limiter), the **guitar amp** plugin (already
reusing the EQ core), a future **browser / WASM** build, and eventually **hardware with a processor
inside** (SoC / DSP pedal).

The model already works: TabbyEQ's `teq::` EQ core is JUCE-free and is **already reused by the guitar
amp plugin** (OrbitCab) for its tone controls. This document promotes that one-off reuse into a
deliberate, modular foundation — and adds the constraints that WASM + embedded targets impose, so we
bake them in now (while the core is small) instead of retrofitting.

---

## 1. Why a shared core (and why JUCE-free)

- **Reuse pays off, proven.** `teq::` already serves two products. Dynamics, convolution, limiting,
  analysis will each be reused the same way → write once, harden once, ship everywhere.
- **JUCE doesn't go where we want to go.** A plain C++ DSP core compiles to **WASM** (Emscripten)
  and to **embedded** toolchains (arm-none-eabi / RTOS / bare-metal); the JUCE framework does not.
  Keeping the core JUCE-free is the entire point — JUCE (or any host framework) lives only in thin
  per-platform *adapters*.
- **Consistency + one test surface.** Same convolver, same detector, same analyzer across products,
  with one self-test suite per module (the `teq` discipline: measured audio == analytic curve).

**Platform priority (product-driven, decides what we fund):**
1. **Desktop** plugins + apps — **primary** (where the revenue is). This is what ships and what CI gates.
2. **WASM** (browser demos / lightweight web) — **secondary, aspirational**.
3. **Embedded / hardware** — **far-future, low priority.**

We keep the core JUCE-free *so tiers 2–3 stay open*, but we **don't pay their tax early** (fixed-point,
no-heap, embedded CI) — a tier earns its CI gate and its constraints only when a product funds it.

---

## 2. Portability laws (non-negotiable — keep WASM + embedded reachable)

These are what make the core actually portable. Violating one silently breaks a target.

**Target tiers (each module declares which it supports).** Generic "WASM + embedded" is too coarse —
name the profiles. Per the **platform priority (§1)**, **CI gates `desktop` (primary) today, with
`wasm-audio` aspirational; both embedded tiers stay documented-but-not-gated** until a hardware product
exists (3rd review: don't pay embedded's CI tax early).
- `desktop` — heap in `prepare`, SIMD, JUCE only in adapters.
- `wasm-audio` — AudioWorklet-safe; no blocking / filesystem in the audio path; bounded memory growth;
  **no FP-control register → relies on software denormal flushing (law 8)**.
- `embedded-fpu` — float; fixed max channels / block / IR / model; no heap after init. Reachable from
  the same code via templating.
- `bare-mcu` — **fixed-point: a SEPARATE codebase, not a `using Sample` flag-flip** (different filter
  topology / overflow scaling / coeff quantization). Documented, not promised from one source. Excludes
  NAM/Eigen and long-IR convolution.

1. **No threads in the core.** Threading is the *adapter's* job. The core is called synchronously.
   *(Why: WASM threads need SharedArrayBuffer + COOP/COEP; bare-metal may have none. Precedent:
   OrbitCab builds/loads NAM models + IRs on a background/message thread and atomic-swaps them into the
   live pointer — the DSP core (`cab::AmpStage::process`, `cab::Convolver`) never spawns or blocks.)*
2. **No allocation / lock / IO / syscall / throw in `process()`** (the existing RT rule). Preallocate
   in a `prepare()` step. Ideally no heap at all in hot classes (fixed-size state) so bare-MCU works.
3. **Float in the hot path; `double` only in offline coefficient design.** No FPU → `double` is
   expensive. (teq already designs coeffs in `double`, processes samples in `float`.) Introduce a
   `using Sample = float` alias and keep raw `float*` out of public signatures **now**, so full
   sample-type templating later is a flag-flip, not a fork-rewrite. *(Carve-out: coefficient recompute
   on the audio thread, and meter/LUFS/true-peak accumulators, may legitimately use `double` — specify
   per module; the ban is on gratuitous `double` in the sample loop.)*
4. **Dependencies behind a seam, never hard-wired.** Heavy primitives (FFT, neural inference) are
   reached through a thin interface so each platform plugs its own impl: JUCE adapter →
   `juce::dsp::FFT`; WASM / embedded → pffft / kissfft / CMSIS-DSP. (Precedent: both plugins run the
   spectrum-analyzer `juce::dsp::FFT` *in the adapter/UI* — TabbyEQ `ui/EqCurveDisplay.h:111`, OrbitCab
   `ui/SpectrumAnalyser.h` — while the `teq::`/`cab::` cores stay FFT-free; OrbitCab also hides
   `juce::dsp::Convolution` behind `cab::Convolver`'s raw-float seam. NB: TabbyEQ has **no** linear-phase
   convolution yet — that is roadmap, not shipped.)
5. **Configurable sizes.** Not just `kMaxChannels` (16 desktop → mono/stereo on hardware) but also
   `kMaxBlockSize` / `kMaxBands` / `kMaxPartitions` / `kMaxIrSamples`. Set via **template param /
   policy + a per-tier CMake preset**, not a fragile pre-include `#define`. Each module reports its
   memory footprint.
6. **No exceptions / RTTI / OS / filesystem / locale** in core code. No global mutable state (watch
   header-only function-local statics).
7. **C++ subset that all toolchains accept** (C++20 desktop; keep an eye on what Emscripten and the
   embedded toolchain support).
8. **Denormals: every feedback kernel flushes in SOFTWARE; hardware FTZ is a desktop optimization,
   never a correctness crutch.** *(Fixed by the 3rd review — agy/Gemini found the hole: WASM and many
   embedded ARMs expose no FP-control register, so an "adapter sets FTZ, core assumes it" rule silently
   fails on the `wasm-audio` tier → feedback filters decay into subnormals on silence → 10–100× CPU
   spike.)* `teq` already does the right thing: `Biquad/Svf::flushDenormals()` zaps `|state| < 1e-15f`
   to exact zero every block, so the state never reaches the subnormal range — `eq` is denormal-safe
   **without** hardware FTZ. **The new `dynamics` envelope/release followers (and any future feedback
   module) MUST adopt the same per-block software flush.** The core never sets a global FTZ/DAZ mode; an
   adapter MAY set FTZ on desktop as a bonus. **Do NOT use `-ffast-math`** (it breaks the NaN/inf
   semantics the tests assert). Every kernel also ships a scalar fallback + a scalar↔SIMD parity test.

**These laws are CI-enforced for the funded tiers, not aspirational.** From day one: a
no-allocation-in-`process()` test, `-fno-exceptions` / `-fno-rtti` build configs, and an **Emscripten**
(`wasm-audio`) build — so a thread / alloc / dep / denormal creep fails the build immediately. An
**arm-none-eabi** job is documented but **not gated** until an embedded product funds it (3rd review:
avoid embedded-grade CI scope creep). Third-party deps (Eigen, kissfft, pffft) must be verified to
build under `-fno-exceptions` (some use throwing asserts).

---

## 3. Module layout

One umbrella, **independent modules**, each separately includable so a consumer pulls only what it
needs (an EQ plugin must NOT drag in the neural runtime). **Hybrid build model (panel decision):** the
*light* modules (`core`, `eq`, `dynamics`, simple analysis) are header-only `INTERFACE` targets; the
*heavy / platform-specific* ones (`convolution`, `neural`, the FFT backends, SIMD kernels) are
**compiled `STATIC`/`OBJECT` targets with narrow public headers**, so intrinsics, denormal handling,
size config, and per-platform flags stay out of public headers (no `#ifdef`-soup, ODR hazards, or
binary bloat). Every module is `felitronics::<module>`, JUCE-free, declares its target tiers, and has
its own JUCE-free self-tests.

| Module          | What                                                                 | Build / deps / notes |
|-----------------|---------------------------------------------------------------------|--------------|
| `core`          | Math, `Smoother`, `ScopedFlushToZero`, `kMaxChannels` (SSOT) + size config, fixed-size ring / lock-free SPSC FIFO, `Sample` alias, the **FFT seam** | header-only, zero deps; the shared base |
| `eq`            | matched biquads (Vicanek) + Cytomic SVF + `EqBand` + `EqEngine`      | header-only = today's `teq::` (becomes `eq`) |
| `dynamics`      | `EnvelopeFollower` (peak/RMS, attack/release) + `GainComputer` (threshold/ratio/knee/range, downward+upward) | header-only, zero deps; **first NEW module** |
| `analysis`      | spectrum tap (the existing `SpectrumTap`), correlation, LUFS/loudness | header-only; FFT **via the seam** |
| `convolution`   | partitioned (uniform / **zero-latency**) FFT convolution            | **compiled**; FFT **via the seam** |
| `oversampling`  | polyphase up/down-sampling                                           | compiled if SIMD; for true-peak + nonlinear |
| `limiter`       | true-peak limiter                                                    | compiled; uses `oversampling` |
| `neural`        | a thin **inference-object seam** (process-only); model *loading* lives in the adapter | **compiled, heavy, isolated**; backend chosen at build time per tier (NAM/Eigen desktop+wasm). NOT a runtime model-swap |

**The FFT seam** (`core/Fft`) — the architecture's keystone, **designed FIRST** (see §6). NOT a bare
`virtual fft(float*)`: a **plan object created in `prepare`** with explicit scratch ownership,
supported sizes, transform direction, normalization convention, packed real/complex bin layout, and
alignment — plus scalar↔backend parity tests. `convolution` owns its *own* algorithm seam on top
(zero-latency partitioning, latency reporting, IR replacement/crossfade, tail flush), not just raw FFT
calls. Adapters plug an engine per tier (`juce::dsp::FFT` desktop · pffft/kissfft wasm · CMSIS
embedded — CMSIS uses a different / fixed-point layout, so the seam must not collapse to a
lowest-common-denominator that kills desktop performance).

---

## 4. What stays OUT of the core (lives in adapters / consumers)

- Threading, scheduling, background rebuilds.
- Parameter systems (JUCE APVTS, etc.), state save/load, GUI.
- Host glue (AudioProcessor, CLAP, AU), file/IR loading, **neural model loading**, the concrete FFT
  engine, setting the denormal CPU mode.
- **Cross-module "glue."** Example: TabbyEQ's *dynamic EQ* = `eq` band + `dynamics` detector fed by
  the band's band-pass, applying the computed gain as an SVF gain-delta. That composition is
  TabbyEQ-specific and lives in TabbyEQ, **not** in `eq` or `dynamics`. A standalone compressor uses
  the same `dynamics` broadband. The core gives primitives; products compose them.

---

## 5. Build & repo strategy

- **Hybrid targets** (see §3): light modules header-only `INTERFACE`, heavy modules compiled
  `STATIC`/`OBJECT`; the umbrella `CMakeLists.txt` exposes `felitronics::core`, `felitronics::eq`, …
  Consumers link only what they use.
- **Consumption:** pinned tag via CMake `FetchContent` (exactly how the products already pin JUCE).
  OrbitCab already consumes `teq` this way (`SOURCE_SUBDIR teq`, pinned tag) — so `teq → felitronics::eq`
  is just repointing the URL/tag, not un-vendoring a folder copy.
- **Versioning is a contract, not just tags.** SemVer + `CHANGELOG.md` + a deprecation policy (compat
  aliases that warn). **DSP output is versioned behaviour** — a changed filter curve, limiter release
  shape, or convolver latency breaks presets/sessions as surely as an API rename → treat it as a
  breaking change. A **consumer-matrix CI** builds TabbyEQ + OrbitCab against each release; golden-audio
  vectors guard behaviour.
- **Portability CI for the funded tiers, from day one** (not "later"): the no-alloc /
  `-fno-exceptions` / `-fno-rtti` configs, scalar↔SIMD parity tests, and an **Emscripten**
  (`wasm-audio`) build of the light modules. An **arm-none-eabi** job is documented but **not gated**
  until an embedded product funds it (see §2).
- **License:** AGPL-3.0-or-later, SPDX header on every file. Each heavy module records its third-party
  deps + AGPL-compatibility in `THIRD_PARTY_NOTICES.md` (NAM MIT, Eigen MPL-2.0, nlohmann/json MIT,
  pffft/kissfft BSD = OK; watch ONNX/others).

---

## 6. Migration plan (keep every product green throughout)

1. **Establish conventions + the `core` API contract** (this doc): `prepare()` / `process()`, block
   views, the size constants, error returns, `noexcept` audio APIs, and the no-alloc test harness.
2. **FFT-seam spike FIRST** *(panel decision — the keystone)*: define `core/Fft` (the plan object) and
   prove a **JUCE-free zero-latency partitioned convolution** against synthetic IRs + one real backend,
   *before* building the core around the seam. This de-risks `convolution` / `analysis` / `limiter` and
   OrbitCab's de-JUCE-ing — the single highest architectural risk.
3. **Land `dynamics`** (the next TabbyEQ feature) within the established conventions. See TabbyEQ's
   upcoming `docs/DYNAMICS.md`.
4. **Move `teq/` → `felitronics-core/eq`**, keeping a thin `teq::` compat alias during transition;
   TabbyEQ + OrbitCab repoint their `FetchContent` pin one at a time.
5. **Extract heavier modules** (`convolution`, `limiter`, `neural`) as they prove out — `convolution`
   built on the spike's seam.

---

## 7. Consumers (each product fills in its own section)

### TabbyEQ (premium track EQ) — desktop plugin
The **reference consumer**: TabbyEQ's `teq::` core IS today's `eq` module (header-only, JUCE-free,
CMake `INTERFACE` target `teq::core`). The JUCE adapter (`src/`) maps APVTS → `teq::BandParams` and
calls the engine. Formats VST3/AU/Standalone/CLAP; future target a **WASM** demo. The de-JUCE state
here is the *goal* state — only the adapter touches JUCE.

**1. What `eq` actually is (verified in `teq/`).**
- `teq::EqEngine` — a bank of **24 bands** (`EqEngine.h:33`) over **≤16 channels** (`EqTypes.h:18`, the
  SSOT every per-channel state array derives from), in series, with built-in **pre/post `SpectrumTap`s**
  (`EqEngine.h:103-112`, fed channel 0 — `:70,72`) and a race-free GUI magnitude readout
  (`magnitudeDbFor`, `:93`).
- `teq::EqBand` runs **two engines under one band** (`EqBand.h:99-106`): a Nyquist-accurate **matched
  biquad** cascade (Vicanek, up to 96 dB/oct = 8 sections, `MatchedBiquad.h`) for static treatment, and
  a **Cytomic SVF** (`teq::Svf`, applied at `:283`) for the swept/search band (clean under fast `fc`
  sweeps). Coeffs recompute once per block, skipped when settled (`:192-194`).
- **M/S dual-lane**: in `ms` mode each band runs an independent Mid (col 0) + Side (col 1) design
  (`EqBand.h:200-216`; Side params `EqTypes.h:46-55`); mono/surround = Mid lane only.
- **Software denormal flush** in both kernels (`Biquad::flushDenormals` `MatchedBiquad.h:426`,
  `Svf::flushDenormals` `Svf.h:81`, per block via `EqBand::flushState` `:286`) → `eq` is already
  `wasm-audio`-safe without hardware FTZ (Law 8 — the reference pattern). JUCE-free **self-tests**
  (`teq/tests/`, ~227 checks: measured audio == analytic curve, per-channel independence, NaN/denormal).

**2. The adapter (`src/`) — the only JUCE surface.**
- `processBlock(juce::AudioBuffer<float>&)` + `juce::ScopedNoDenormals` (`PluginProcessor.cpp:69-71`);
  reads APVTS atomics → `teq::BandParams` (`readBand`, `:165-188`); feeds `setBand` then `process` on
  the audio thread (the engine's same-thread contract, `:139-142`).
- Output trim = `juce::LinearSmoothedValue` (`PluginProcessor.h:95`) — a **de-JUCE candidate** (→ `core`
  `Smoother`). The drag-audition / band-solo path uses `teq::Svf` directly in the adapter
  (`PluginProcessor.cpp:101-135`) — already JUCE-free DSP.
- GUI spectrum FFT = `juce::dsp::FFT` (`ui/EqCurveDisplay.h:111`) — the real "FFT lives in the adapter,
  core stays clean" precedent (mirrors OrbitCab's analyzer).

**3. Future core modules TabbyEQ will use.**
- **`dynamics`** (NEW — the next feature; full design in `tabby-eq/docs/DYNAMICS.md`): core primitives
  `EnvelopeFollower` + `GainComputer` (JUCE-free, software-flush per Law 8), developed inside `teq/`
  first and migrated with the `eq` extraction. The **dynamic-EQ composition** (detector = a Cytomic-SVF
  band-pass probe on the band's region; gain applied as a **matched-static × SVF gain-delta**, per M/S
  lane; GR-metering atomics for the host/UI) stays **product glue in TabbyEQ** (the §4 rule). De-esser
  = a preset on the same engine.
- **`analysis`** ← the pre/post `SpectrumTap`s (already JUCE-free here) + a future correlation meter.
- **`convolution` (linear-phase) — NOT IMPLEMENTED YET.** *Correction to the prior claim:* linear-phase
  is **deferred / roadmap**, not "already implemented in the adapter" — there is **no convolution / FIR
  / overlap-add in TabbyEQ's code** today (only the analyzer `juce::dsp::FFT`); docs list it as later
  (`docs/ROADMAP.md:10`, `docs/MS-DUAL-MODE.md:26`). When built it rides `felitronics::convolution` on
  the FFT seam. The dynamics×Linear interaction (dynamics is time-varying → bypassed in Linear mode,
  `docs/DYNAMICS.md:112-117`) is a constraint for *future* features, not current behaviour.
- possibly **`limiter`** on the output (future).

### Guitar amp plugin (OrbitCab) — desktop plugin
Signal chain: `input → preamp (NAM) → tone EQ (teq) → poweramp (NAM) → cab (IR convolution) → output`.
The headless DSP lives in `src/core/` (namespace `cab::`); the JUCE adapter (APVTS / processor / editor)
is `src/`. Formats VST3/AU/CLAP/Standalone.

**1. Reuse of `teq::` — already on the FetchContent model, not a copy.**
- **Correction to §1/§5's premise:** OrbitCab does **not** copy-vendor `teq/`. It pulls it via CMake
  `FetchContent` from `github.com/darwinscat/tabby-eq`, **pinned tag `v0.1.0`**, `SOURCE_SUBDIR teq`
  (only the header-only `teq/` core, never the TabbyEQ plugin or its JUCE), linked as `teq::core`.
  So this product is *already* on the pinned-fetch model the doc targets — there is no folder-copy to
  migrate; the `teq → felitronics::eq` move here is just repointing the `FetchContent` URL/tag.
- **Used for the amp tone stack only:** `cab::AmpEq` (`src/core/AmpEq.h`) wraps `teq::EqEngine`, using
  **6 of teq's 24 bands** — HPF, Bass (low shelf), Mid (bell), Treble (high shelf), Presence (high
  shelf), LPF — at fixed (generic) frequencies. Recorded in `THIRD_PARTY_NOTICES.md`.
- **Footprint gotcha (relevant to the `eq` module):** `teq::EqEngine` is ~**200 KB** (a fixed 24-band
  bank), so `cab::AmpEq` holds it on the **heap** (`unique_ptr`, built in `prepare()`), never by value
  — a by-value member overflowed MSVC's **1 MB** main-thread stack (the integration test stack-allocates
  the processor; macOS's 8 MB hid it). → keep `eq` heap-placeable and make `kMaxBands` shrinkable (law 5).
- Note: the cab's *own* per-slot HPF/LPF are **not** teq — they are `juce::dsp::StateVariableTPTFilter`
  (see §6). `teq` is the amp EQ only.

**2. Amp-specific DSP we have today** (all `cab::`, `src/core/`):
- **Neural amp (preamp + poweramp):** `cab::AmpStage` (`AmpStage.{h,cpp}`) over **NeuralAmpModelerCore**
  (sdatkinson, pinned commit `b5a68c3…`), inference on **Eigen** (MPL-2.0) + nlohmann/json (MIT),
  `NAM_SAMPLE_FLOAT` (float hot path). **Not RTNeural, not ONNX.** WaveNet/LSTM/ConvNet architectures
  self-register via static initializers → must be linked `WHOLE_ARCHIVE`. Two instances = preamp +
  poweramp. Effectively **JUCE-free** (NAM + std + StreamResampler). Model load is off-thread + atomic
  swap (threading in the adapter, law 1).
- **`cab::StreamResampler` (`StreamResampler.h`):** **JUCE-free** Catmull-Rom resampler that rate-matches
  a model's native SR to the host SR; the only source of reported latency (0 when SRs match).
- **IR cabinet convolution:** `cab::Convolver` (`Convolver.h`) + `cab::IRSlot` (`IRSlot.{h,cpp}`) over
  **`juce::dsp::Convolution`** in **zero-latency** mode. **JUCE-dependent** (see §6 conflicts).
- **Per-slot cab HPF/LPF:** `juce::dsp::StateVariableTPTFilter` (12 dB/oct Butterworth), pre-convolution,
  per A/B slot. **JUCE-dependent.**
- **Auto-leveler:** `cab::AutoLeveler` (`AutoLeveler.h`) — wet/dry RMS followers + silence gate →
  makeup `sqrt(dryMS/mixMS)`. JUCE-free math.
- **Spectrum tap:** `cab::SpectrumTap` (`SpectrumTap.h`) — the same struct as the doc's `analysis`
  SpectrumTap, but a **JUCE-coupled diverged copy** (`juce::FloatVectorOperations`, `SpectrumTap.h:6,37`);
  the JUCE-free `teq::SpectrumTap` twin is the consolidation base (see §7.3 + the shared-DSP review).
- **NOT present today:** dedicated noise gate, drive / boost / clipper / waveshaper, oversampling.
  ("boost" is a NAM *capture* variant, not DSP; "gate" in the code = bypass / mute / auto-level gating.)

**3. Core vs product-specific — recommended boundaries.**
- **Lift to the shared core:**
  - `eq` ← already `teq`. ✓
  - `analysis` ← `SpectrumTap` — **but `cab::SpectrumTap` is NOT JUCE-free** (`SpectrumTap.h:6,37` use
    `juce::FloatVectorOperations`); the **`teq::SpectrumTap` copy IS** (`std::copy` + `reset()` +
    `tryPull()`). They are **diverged copies** of one struct → consolidate to the `teq::` shape as
    `felitronics::analysis::SpectrumTap`; OrbitCab swaps `FloatVectorOperations`→`std::copy` and
    `CabEngine::pullSpectrum` (`CabEngine.cpp:323`) → `tap.tryPull()`. (See the shared-DSP review.)
  - `neural` ← the NAM/Eigen runtime — **isolated behind the inference seam** so EQ/comp products never
    drag it in. (Runtime is NAM, not RTNeural; model loading stays in the adapter, see §8.)
  - `convolution` ← the cab IR conv — but only **after** a JUCE-free FFT impl exists behind the seam
    (today it is `juce::dsp::Convolution`); must preserve **zero-latency** partitioning.
  - A small **resampler** util (StreamResampler) — generic rate-match, candidate for `core`.
- **Stays product-specific (OrbitCab):** the chain glue (`cab::CabEngine`: stage order, dry/wet, A↔B
  mix, phase, trim, auto-level), the tone-stack **voicing** (fixed freqs / future per-model measured
  stacks — composition, not a primitive), the NAM **library/selector** (`PreampLibrary` /
  `PowerampLibrary` filename → channel/gain/PP·SE/hours), IR loading, UI, APVTS. `AutoLeveler` is a thin
  matcher — leave product-specific unless a second product wants it.

**4. Target platforms.**
- **Shipping: desktop only** — macOS universal / Windows x64+arm64 / Linux x64+arm64.
- **WASM:** aspirational, not built. *Partly* ready — `cab::Params` is deliberately JUCE-free ("compiles
  under Emscripten / embedded"); `teq` + `cab::AmpStage` are JUCE-free. But the **cab path
  (convolution + SVF + `juce::AudioBuffer`/`SmoothedValue`/`FloatVectorOperations`) is JUCE-coupled** →
  not WASM-ready as-is.
- **Embedded / hardware:** **no concrete SoC / fixed-point / no-heap target defined for this plugin yet.**
  Honest blocker: the neural amp (**Eigen + WaveNet/LSTM**) is desktop/WASM-grade, not bare-MCU-grade — a
  real hardware amp would need a *different, lighter* neural runtime (RTNeural fixed-size, or a tiny
  model), not NAM-on-Eigen. So "embedded amp" changes the `neural` impl, not just shrinks constants.

**5. Constraints the core must respect for the amp.**
- **Zero-latency convolution** — the `convolution` module must offer a zero-latency partitioned mode
  (not only uniform-partitioned with PDC); OrbitCab advertises ~0-sample latency.
- **Separable `prepare()`/`process()` + swap-safe hot path** — NAM models and IRs are large; build/load
  on the message thread, atomic-swap into the live pointer (law 1, no alloc in `process()`).
- **Float hot path** — NAM (`NAM_SAMPLE_FLOAT`) and teq both process float (law 3 ✓).
- **Large fixed-size state** — `teq::EqEngine` ~200 KB; `eq` must stay heap-placeable, `kMaxBands`
  shrinkable for small RAM.
- **Licenses (all AGPL-compatible):** NAM = MIT, Eigen = MPL-2.0, nlohmann/json = MIT, teq = AGPL,
  JUCE = AGPL option. A JUCE-free `convolution` FFT must use a permissive lib (pffft / kissfft BSD = OK).

**6. Constraints / conflicts with the portability laws.**
- **Law 4 (deps behind a seam) is VIOLATED by the cab path today.** `cab::Convolver` hard-wires
  `juce::dsp::Convolution` and `cab::IRSlot` the JUCE SVF — so OrbitCab's "core" is **not fully
  JUCE-free**: only `teq` (eq) and `cab::AmpStage` (neural) are. Lifting the cab path into
  `felitronics::convolution` requires real de-JUCE-ing: a JUCE-free FFT behind the seam, replacing the
  JUCE SVF (teq's Cytomic SVF is the natural swap), and dropping `juce::AudioBuffer`/`SmoothedValue`/
  `FloatVectorOperations` for `core` primitives. **This is the guitar plugin's main migration cost — size
  it explicitly.** (This is exactly why the panel put the FFT-seam spike first, §6.)
- **Neural-on-embedded conflicts with "bare-MCU friendly."** NAM-on-Eigen will not fit a small MCU →
  the inference seam must allow a *different backend per build* (not a runtime swap of one model), §8.

### Future products
- Standalone **compressor** (Felitronics) — `dynamics` broadband + optional sidechain `eq`.
- **True-peak limiter** — `oversampling` + `limiter`.
- **Browser / WASM** demos of any of the above.
- **Hardware** pedal/box with a processor inside.

---

## 8. Open decisions

- **RESOLVED (panel): build model = hybrid** — light modules header-only `INTERFACE`, heavy modules
  compiled `STATIC`/`OBJECT`. §3, §5.
- **RESOLVED (panel): the FFT-seam + zero-latency convolution spike goes FIRST**, before dynamics. §6.
  **DONE (2026-06-29) — landed + green:** `felitronics::core::fft` (compile-time seam + a scalar radix-2
  reference backend), `felitronics::convolution::PartitionedConvolver` (zero-latency direct-head +
  uniform-partitioned overlap-save tail; verified == direct convolution under hostile variable block
  splits + boundary impulses; no-alloc-in-`process`), and an offline Kaiser `resampleIr` (≥55 dB). 120
  checks across the FFT / convolution / resampler suites. The seam is a **compile-time backend with an
  OPAQUE spectrum + a backend `spectralMultiplyAdd`** (so pffft/vDSP never pay an O(N) repack — the
  Gemini fix); real backends (pffft/kissfft/juce::dsp::FFT) plug in later. Production still needs a
  crossfade on live IR swap + Gardner non-uniform partitions for long IRs (both noted, not spike scope).
- **Sample type:** add a `using Sample = float` alias + non-`float*`-locked signatures **now**; full
  templating reaches `embedded-fpu` (float). **`bare-mcu` (fixed-point) is a SEPARATE codebase** — the
  alias does NOT flip a float SVF / NAM into fixed-point (different topology / scaling / quantization);
  don't promise it from one source (3rd review).
- **RESOLVED (3rd review): the FFT seam is compile-time** (template / C++20 concept), **not `virtual` in
  the convolver hot path** (a vtable kills inlining / vectorization). The plan-object contract (sizes /
  scratch / layout / normalization) is finalized in the spike. *(Nuance: the FFT runs ~once per partition
  block, not per sample, so a vtable wouldn't be catastrophic — but compile-time stays the default.)*
- **Double carve-outs:** specify, per module, where `double` is allowed on the audio thread
  (coefficient recompute) and for accumulators (RMS / LUFS / true-peak).
- Single repo (this) vs multi-repo per module. (Lean: single repo, multiple targets — possibly layered:
  `core` / `dsp` (eq+dynamics+filters+meters) / `fft_backends` / `convolution` / `neural`.)
- Naming: `felitronics::` vs keep `teq::`/`dyn::` per-module namespaces under the umbrella.
- **`convolution` is the first real de-JUCE-ing job** (raised by the guitar amp, §7.6): OrbitCab's cab
  path uses `juce::dsp::Convolution` + the JUCE SVF, so it's *not* JUCE-free today. Choosing the FFT-seam
  engine (pffft / kissfft / CMSIS) and porting the partitioned (zero-latency) convolver off `juce::dsp`
  is what unblocks both `felitronics::convolution` and the guitar plugin's WASM/embedded path → the spike.
- **`neural` seam = inference-object level, build-time backend choice** (refined by the panel): the core
  defines a process-only inference interface; **model loading lives in the adapter** (the core accepts a
  pre-built inference object). "Per-platform" = the adapter links a different backend (NAM/Eigen on
  desktop+WASM; a lighter runtime / tiny model on hardware) — NOT a runtime swap of one model. Don't list
  RTNeural as a dependency until something actually uses it.

- **RESOLVED (shared-DSP review, 2026-06-29 — `dsp-shared-dsp-review.md`): `analysis::SpectrumTap` = the
  `teq::` shape.** `cab::SpectrumTap` and `teq::SpectrumTap` are diverged copies of one SPSC struct;
  the `teq::` one is the strict superset (JUCE-free `std::copy` + `reset()` + `tryPull()`). OrbitCab
  adopts it; the mono-sum vs channel-0 *feeding* policy stays per-product.
- **RESOLVED (shared-DSP review): SVF + Smoother dedup = REPLACEMENT, not code-merge.** OrbitCab's
  per-slot `juce::dsp::StateVariableTPTFilter` → `teq::Svf` (same Zavalishin TPT topology; HP/LP @
  12 dB/oct, Butterworth at Q≈0.707) and its `juce::SmoothedValue`s → a `core` `Smoother`. Both change
  the impl, not the math family → a **versioned-behaviour** swap (guard with a magnitude/phase golden test).

---

## 9. Panel review — 2026-06-29

Reviewed by an independent panel — **Codex**, **DeepSeek**, and (added later) **agy / Gemini 3.1 Pro**
run from the OrbitCab session; full notes in the companion `dsp-architecture-third-review.md`. All three
**independently converged**:

> **Verdict: proceed-with-changes.** The JUCE-free, seam-based, modular core is sound; the module
> decomposition is sensible — but several seams were so underspecified that building heavy modules
> around them would force a costly re-architecture.

> **#1 recommendation (all three, independently):** design the **FFT / convolution seam + a zero-latency
> partitioned-convolution prototype FIRST**, before expanding the core around an unproven seam — it's
> the keystone the heaviest modules (`convolution`, `analysis`, `limiter`) and OrbitCab's de-JUCE-ing
> all depend on. "FFT seam ≠ `virtual fft(float*)`."

Accepted-changes checklist (folded into §2–§8 above):

- [x] Build model → **hybrid** (light header-only, heavy compiled). §3, §5.
- [x] **FFT-seam spike goes first**, before dynamics. §6.
- [x] Explicit **target tiers** (desktop / wasm-audio / embedded-fpu / bare-mcu); modules declare support. §2.
- [x] FFT seam = **plan object** (scratch / sizes / layout / normalization); `convolution` owns its algorithm seam. §3.
- [x] **Neural seam = inference-object**; model loading in adapter; build-time backend, not runtime swap. §3, §8.
- [x] **Software denormal flush in every feedback kernel** (teq's `<1e-15 → 0` per block; FTZ a desktop
      bonus only; no `-ffast-math`) + scalar↔SIMD parity tests. §2. *(corrected by the 3rd review)*
- [x] **FFT seam = compile-time** (template/concept), not `virtual` in the hot path. §8.
- [x] Size config beyond `kMaxChannels` (block / bands / partitions / IR) via template/policy + CMake preset; footprint reporting. §2.
- [x] **Versioning is a contract**: SemVer + CHANGELOG + deprecation; *DSP output is versioned behaviour*; consumer-matrix + golden-audio CI. §5.
- [x] `using Sample = float` alias + non-`float*`-locked signatures **now**. §2, §8.
- [x] **CI-enforced laws from day one (funded tiers)**: no-alloc-in-`process`, `-fno-exceptions` / `-fno-rtti`, Emscripten (`wasm-audio`); arm-none-eabi documented-not-gated. §2.
- [x] `double` carve-outs (coeff recompute on audio thread; meter/LUFS/true-peak accumulators). §2, §8.

### Third reviewer — agy / Gemini 3.1 Pro (2026-06-29, via the OrbitCab session)
Converged with Codex + DeepSeek on the sound calls (FFT-first, neural = build-time backend,
versioning-as-behaviour, hybrid build) and found **one real correctness hole + two scope refinements**,
all **verified against the actual code** (see `dsp-architecture-third-review.md`):
- **Law 8 was wrong for WASM (now fixed, §2):** WASM / embedded ARM expose no FP-control register, so
  "adapter sets FTZ, core assumes it" silently fails. `teq` already SOFTWARE-flushes (`<1e-15 → 0` per
  block) so `eq` is WASM-safe; the rule is now "every feedback kernel software-flushes" — the new
  `dynamics` followers must too. `-ffast-math` rejected (breaks the NaN/inf the tests assert).
- **`bare-mcu` fixed-point = a separate codebase**, not a `using Sample` flag-flip (§2, §8).
- **Tiers CI-gated only when a product funds them** (drop `bare-mcu` / arm-none-eabi from CI for now) (§2).
- **FFT seam resolved to compile-time** (template / concept), not `virtual` in the hot path (§8).

Still to do: resolve the remaining §8 open decisions during the FFT-seam spike.

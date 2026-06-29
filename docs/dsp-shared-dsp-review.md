<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->
<!-- Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa. -->

# Cross-plugin shared-DSP review + de-JUCE plan

**By:** the felitronics-core lead session · **Date:** 2026-06-29 · companion to `DSP-ARCHITECTURE.md`
and `dsp-architecture-third-review.md`. Code-grounded (every claim cites `file:line`, verified against
the actual sources, not taken on a reviewer's word). Panel: **Codex + DeepSeek + Gemini** (see §6).

**Priority frame (product-set, governs every call here):** desktop plugins/apps are **#1** (ship +
CI gate); `wasm-audio` is secondary "only if it doesn't tax desktop"; embedded is below the floor. The
core stays JUCE-free because that is what makes it *shareable and desktop-testable* — not for the lower
tiers' sake.

---

## 0. Status — what is already built (this session)

`felitronics-core` is no longer docs-only. Four light modules are migrated/written, JUCE-free, and
**green from a clean build** (`cmake -S . -B build -DFELITRONICS_BUILD_TESTS=ON && ctest`, AppleClang,
C++20):

| Module | Namespace | Contents | Self-tests |
|---|---|---|---|
| `core` | `felitronics::core` | `Config`(`Sample`,`kMaxChannels` SSOT,`kMaxBlockSize`), `Math`(`kPi`,dB↔gain), `Smoother`, `FlushToZero`(`flushDenormal`+`ScopedFlushToZero`) | pass |
| `analysis` | `felitronics::analysis` | `SpectrumTap` (consolidated, see §2) | pass |
| `dynamics` | `felitronics::dynamics` | `EnvelopeFollower` + `GainComputer` (the first NEW module) | **15 checks** pass |
| `eq` | `felitronics::eq` (+ `teq::` compat) | matched biquads + Cytomic SVF + `EqBand` + `EqEngine`, migrated from `teq/` | **596 checks** pass (the full teq suite, run verbatim through the `<teq/*.h>` shims) |

The 596-check pass is the proof that the `eq` re-home is **behaviour-identical** to TabbyEQ's proven
`teq::` core. Heavy modules (`convolution`, `neural`, `oversampling`, `limiter`) remain future work
gated on the FFT-seam spike (§4).

---

## 1. What is genuinely shared between OrbitCab and TabbyEQ

| Concern | OrbitCab (`cab::`) | TabbyEQ (`teq::`) | Verdict |
|---|---|---|---|
| **EQ engine** | `cab::AmpEq` wraps `teq::EqEngine`, 6 of 24 bands at fixed amp-tone freqs (`orbitcab/src/core/AmpEq.h:39-64`) | `teq::EqEngine` full 24-band/16-ch/M·S (`tabby-eq/teq/include/teq/EqEngine.h:33-34`) | **Already shared** via FetchContent (`orbitcab/CMakeLists.txt:106-114`). → `felitronics::eq`. |
| **Spectrum tap** | `cab::SpectrumTap` (`SpectrumTap.h`) | `teq::SpectrumTap` (`SpectrumTap.h`) | **Diverged copies** → consolidated (§2). |
| **Param smoothing** | `juce::SmoothedValue` (`CabEngine.h:123-127`, `AutoLeveler.h:80`) | `teq::Smoother` + `juce::LinearSmoothedValue` in the adapter (`src/PluginProcessor.h:95`) | **Replacement** → `felitronics::core::Smoother` (§2). |
| **2-pole TPT SVF** | `juce::dsp::StateVariableTPTFilter` (`IRSlot.h:56`) | `teq::Svf` Cytomic (`Svf.h`) | **Replacement**, not a merge (§2). |
| **Convolution** | `cab::Convolver` over `juce::dsp::Convolution` (`Convolver.h:70`) | none (linear-phase is roadmap, §5) | OrbitCab-only; → future `convolution` after the spike. |
| **Neural** | `cab::AmpStage` over NAM/Eigen (`AmpStage.{h,cpp}`) | none | OrbitCab-only; → future `neural` seam. |
| **Resampler** | `cab::StreamResampler` (Catmull-Rom, JUCE-free) | none | OrbitCab-only utility; candidate for `core` (but NOT for IR resampling — §4). |

Everything else is **product glue**: OrbitCab's chain order / dry-wet / A·B / auto-level / phase
(`cab::CabEngine`), the tone-stack voicing + NAM library/selector; TabbyEQ's M·S composition, the
search→treat workflow, the de-esser/dynamic-EQ composition. The §4 rule (glue stays in the product)
holds for both.

---

## 2. Dedup verdicts (the explicit dup-candidate adjudication)

### SpectrumTap — **diverged copies → one canonical primitive** (consolidation #1)
Same SPSC struct (2048-point window, `fifo`/`data`/`idx`/`atomic<bool> ready`), but:
- `cab::SpectrumTap` is **NOT JUCE-free** — `#include <juce_audio_basics>` + `juce::FloatVectorOperations::copy` (`orbitcab/src/core/SpectrumTap.h:6,37`); has only `push()`, and its consumer re-arms the tap **externally** (`cab::CabEngine::pullSpectrum`, `CabEngine.cpp:323`).
- `teq::SpectrumTap` is **JUCE-free** (`std::copy`) and a strict superset — adds `reset()` + `tryPull()` (`tabby-eq/teq/include/teq/SpectrumTap.h:31-41`).

→ Canonical = the `teq::` shape, now `felitronics::analysis::SpectrumTap`. **Two desktop-win
optimizations folded in** (from Gemini, verified): (1) the per-sample-written producer cursor
(`fifo`/`idx`) is `alignas(64)`-separated from the GUI-polled handshake (`ready`/`data`) → kills the
false sharing the audio thread would otherwise suffer; (2) templated on FFT order (`SpectrumTapT<N>`,
default 2048) so a different analysis resolution needs no fork. OrbitCab adopts it by swapping
`FloatVectorOperations`→`std::copy` and `pullSpectrum`→`tap.tryPull()`. The mono-sum (OrbitCab) vs
channel-0 (TabbyEQ) *feeding* policy stays per-product (it's caller code, not the tap).

### SVF — **replacement, not a code-merge**
OrbitCab's per-slot cab HPF/LPF are `juce::dsp::StateVariableTPTFilter<float>` at resonance 0.707
(`orbitcab/src/core/IRSlot.cpp:22-27`); `teq::Svf` (now `felitronics::eq::Svf`) is the same Zavalishin
TPT topology family with HP/LP modes (`Svf.h` switch cases). "Dedup SVF" = OrbitCab **moves onto
`felitronics::eq::Svf`** (Butterworth at Q≈0.707, 12 dB/oct) — not a merge of identical code. It is a
**versioned-behaviour swap** (the two impls are not guaranteed bit-identical): guard the de-JUCE with a
magnitude/phase golden test. Bonus: `felitronics::eq::Svf` clamps `f≤0.49·fs` and software-flushes
denormals — both safer than the bare JUCE object.

### Smoother — **replacement** (same class as SVF)
OrbitCab's `juce::SmoothedValue`s and TabbyEQ's `juce::LinearSmoothedValue` → `felitronics::core::Smoother`.
Note this is **exponential** one-pole, where the JUCE ones are linear → a versioned-behaviour change
(inaudible for a 30 ms glide, but guard the OrbitCab smoothers' migration with a golden test).

---

## 3. Core vs product boundary, per module (the lift targets)

- **`eq`** ← `teq` (done). TabbyEQ uses the full engine; OrbitCab keeps `cab::AmpEq` (6-band voicing) as product glue *on top of* `eq`.
- **`analysis`** ← `SpectrumTap` (done, §2). Correlation/LUFS later; the display FFT + window + bin policy stay in each product's UI.
- **`dynamics`** ← `EnvelopeFollower` + `GainComputer` (done). The dynamic-EQ / de-esser composition (detector probe = a `eq::Svf` band-pass on the band region; gain applied as a matched-static × SVF-delta, per M·S lane; GR atomics) stays **product glue in TabbyEQ**.
- **`convolution`** ← `cab::Convolver` — only **after** the FFT-seam spike (§4). Must preserve zero-latency.
- **`neural`** ← `cab::AmpStage` — behind the inference-object seam (§5 mismatch D).
- **`resampler`** — `cab::StreamResampler` is a fine generic rate-match utility for `core`, but is the **wrong tool for IR resampling** (§4).

---

## 4. The de-JUCE-ing cost (the migration's main risk) — desktop-first

**TabbyEQ:** trivial. The `teq::` core is already JUCE-free; the only JUCE is the adapter — `processBlock` + `juce::ScopedNoDenormals` (`src/PluginProcessor.cpp:69-71`), the output `juce::LinearSmoothedValue` (`PluginProcessor.h:95`, → `core::Smoother`), and the analyzer `juce::dsp::FFT` (`ui/EqCurveDisplay.h:111`). Repointing `teq → felitronics::eq` is a CMake change (§7), not a port.

**OrbitCab:** this is the real cost. The `cab::` "core" is **not** fully JUCE-free. JUCE-bound surface:

| Piece | JUCE dependency | file:line |
|---|---|---|
| `cab::Convolver` | `juce::dsp::Convolution` + `AudioBlock` + `ProcessSpec` | `Convolver.h:6,26,45-70` |
| `cab::IRSlot` | `juce::dsp::StateVariableTPTFilter` (→ `eq::Svf`) + `AudioBlock` | `IRSlot.h:56`, `IRSlot.cpp:22-27,120-123` |
| `cab::CabEngine` | `juce::AudioBuffer` scratch (wet/dry) + 6× `juce::SmoothedValue` + `FloatVectorOperations` | `CabEngine.h:118-127`, `CabEngine.cpp:111,323` |
| `cab::AutoLeveler` | `juce::SmoothedValue` | `AutoLeveler.h:80` |
| `cab::SpectrumTap` | `juce::FloatVectorOperations` (→ `analysis::SpectrumTap`) | `SpectrumTap.h:6,37` |

Already JUCE-free (no porting): `cab::Params`, `cab::AmpStage` (NAM behind a pImpl — `AmpStage.cpp`
includes NAM/std, no JUCE), `cab::StreamResampler`.

**Mapping onto the FFT-seam-first spike (the ADR keystone):** the mechanical part (swap
`AudioBuffer`→raw `float*`/`core` block views, `SmoothedValue`→`core::Smoother`, the SVF + tap dedups)
is low-risk and can land incrementally. The **architectural** risk is the convolver, and the spike must
de-risk **more than the raw FFT** — Gemini's refinement, validated: a JUCE-free `juce::dsp::Convolution`
replacement also needs (1) a **high-quality offline IR resampler** (windowed-sinc, ≳60 dB SNR — the
Catmull-Rom `StreamResampler` is too low-SNR for IRs, so it is **not** reusable here), (2) **variable
host block size** handling at zero latency (JUCE does this for free), (3) **RT-safe retired-partition
GC** on IR swap, and (4) the off-thread partition build living in the adapter (Law 1) without leaking
partition layout. Until the spike lands, OrbitCab stays desktop-only on `juce::dsp` and consumes only
the already-JUCE-free pieces (`eq`, then `analysis`, then `neural`).

---

## 5. ADR ↔ real-code mismatches (found + corrected in the ADR)

- **(A) "`cab::SpectrumTap` already JUCE-free" — FALSE.** It uses `juce::FloatVectorOperations`
  (`SpectrumTap.h:6,37`); only the `teq::` twin is JUCE-free. *Fixed:* ADR §7.3 + §7.2.
- **(B) "TabbyEQ linear-phase convolution already implemented in the JUCE adapter" — FALSE.** There is
  **no** convolution / FIR / overlap-add anywhere in `tabby-eq/src` or `tabby-eq/teq` — only the
  analyzer `juce::dsp::FFT` (`ui/EqCurveDisplay.h:111`); linear-phase is explicitly **deferred**
  (`docs/ROADMAP.md:10`, `docs/MS-DUAL-MODE.md:26`). The ADR's §2 law-1/law-4 *precedents* even cited a
  non-existent "TabbyEQ linear-phase FIR builder". *Fixed:* §7 corrected; §2 precedents re-pointed to
  the real one (both plugins run the analyzer FFT in the adapter while the cores stay FFT-free).
- **(C) `cab::AutoLeveler` lacks a denormal flush** on its one-pole RMS followers (`AutoLeveler.h:55-65`).
  *Severity: LOW* — these are `double` accumulators updated **once per block**, not a per-sample feedback
  kernel, so Law 8's "10–100× spike" does **not** apply (re-checked the arithmetic; Gemini's "MEDIUM for
  a 1-sample WASM block" is the only edge, and desktop is #1). Fold a `core::flushDenormal` in when
  `AutoLeveler` migrates to `core` — cheap Law-8 conformance, not a scheduling item.
- **(D) Neural boundary gap.** The ADR's `neural` seam wants the core to take a **pre-built inference
  object**, but `cab::AmpStage::loadModelFromMemory` parses JSON + builds the NAM DSP **inside** the core
  (`AmpStage.cpp:237-256`, with `try/catch`). Today it's isolated (pImpl, message thread) so it works,
  but for the `neural` extraction the parse+build must move to the adapter/factory (a static lib still
  carries the JSON/exception code even behind a pImpl). *Refinement, not a today-bug.*

---

## 6. Panel — convergence + who I corrected

All three ran at max reasoning on (A) an adversarial review of the new `core`/`dynamics`/`analysis`
code and (B) the `eq` re-home design. **Convergence:** the re-home = canonical `felitronics::*` +
`<teq/*.h>` re-export shims (using-declarations), `kMaxChannels` SSOT in `core`, a **real `teq_core`
target** — exactly what shipped (§0, §7).

- **Codex** (strongest): caught the **real `teq_core` target** need (an alias alone breaks OrbitCab's
  `get_target_property(... teq_core ...)`, `orbitcab/CMakeLists.txt:117`) — applied; the **GainComputer
  NaN-param hole** — fixed (finite guards); the RMS-TC nuance (RMS reaches `sqrt(1-1/e)`, not `1-1/e`)
  — documented; confirmed the knee is C1 and the SpectrumTap false-sharing is real.
- **Gemini** (adversarial third): caught the **EnvelopeFollower subnormal-coefficient** spike at an
  absurdly-short attack — **fixed** (flush the coeff); confirmed linear-phase is absent (B). *Corrected:*
  its "Smoother `std::pow` every block = MEDIUM CPU" is **negligible** (~0.05 %) and I will **not** change
  `core::Smoother` (bit-identity keeps the 596 eq checks green); its "GainComputer non-atomic params = bug"
  is the **documented same-thread contract**, not a bug.
- **DeepSeek** (weakest on the code pass — missed the NaN hole and the asm barrier, and wrongly declared
  the RMS follower "no bug"): raised one novel point — the soft-knee "bleeds into the inactive half" for
  up-compress/down-expand. *Corrected by math:* my knee is the **standard symmetric soft knee centred on
  the threshold** (±knee/2) — exactly as down-compress eases in *below* threshold, up-compress eases out
  *above* it; this is textbook, not a bug. On (B) DeepSeek agreed with the shim design but proposed a
  single `teq_compat.h` (we chose per-header shims to keep include granularity) and missed the
  `teq_core` target.

Fixes applied to the new code from this pass: subnormal-coeff flush (`EnvelopeFollower`), finite/NaN
guards (`GainComputer`), an aarch64 `"memory"` clobber on `ScopedFlushToZero`, plus the `alignas`
false-sharing fix already in `analysis::SpectrumTap`. The migrated `eq::Svf` was kept **verbatim**
(Gemini's "flush only active channels" micro-opt deferred — it would perturb the proven, test-locked
code for a negligible gain).

---

## 7. Verification table (re-checked against code / build)

| Claim | Status |
|---|---|
| `eq` re-home is behaviour-identical to `teq` | **VERIFIED** — the full teq suite (596 checks) passes verbatim through the `<teq/*.h>` shims. |
| `cab::SpectrumTap` is JUCE-coupled; `teq::` one is not | **VERIFIED** — `SpectrumTap.h:6,37` (juce) vs `std::copy`+`tryPull`. |
| TabbyEQ has no linear-phase convolution | **VERIFIED** — grep of `tabby-eq/src`+`teq` finds only the analyzer `juce::dsp::FFT`; docs say deferred. |
| `cab` cab-path hard-wires `juce::dsp::Convolution` + JUCE SVF | **VERIFIED** — `Convolver.h:70`, `IRSlot.h:56`. |
| OrbitCab reads the `teq_core` *target* (alias insufficient) | **VERIFIED** — `orbitcab/CMakeLists.txt:117`; the new `teq_core` is a real INTERFACE target. |
| `cab::AmpStage` builds NAM inside the core | **VERIFIED** — `AmpStage.cpp:237-256` (`nlohmann::json::parse` + `nam::get_dsp`). |
| `AutoLeveler` denormal severity is LOW | **VERIFIED (arithmetic)** — per-block `double` accumulator; not a per-sample feedback kernel. |
| New `dynamics`/`core`/`analysis` compile + pass | **VERIFIED** — clean ctest, 0 failures. |

---

## 8. Migration status + next

- **Done:** `core`, `analysis`, `dynamics`, `eq` built + green; ADR §7/§2/§8 corrected; dedup verdicts locked.
- **Next (plugin-side, agent tasks — see `docs/migration/`):** repoint TabbyEQ then OrbitCab's
  `FetchContent` from `teq` → `felitronics-core` (link `teq::core` compat or `felitronics::eq`); then the
  OrbitCab dedups (`SpectrumTap`→`analysis`, SVF→`eq::Svf`, smoothers→`core::Smoother`) each behind a
  golden test.
- **DONE (core-side, keystone):** the **FFT-seam + zero-latency partitioned-convolution spike** landed
  green — `core::fft` (compile-time seam + scalar backend, opaque spectrum + backend `spectralMultiplyAdd`),
  `convolution::PartitionedConvolver` (zero-latency, hostile-variable-block-tested, no-alloc), offline
  Kaiser `resampleIr` (≥55 dB). 120 checks. This unblocks `convolution`, `analysis` FFT, `limiter`, and
  OrbitCab's de-JUCE-ing. Remaining for production: a crossfade on live IR swap + Gardner non-uniform
  tail for long IRs (uniform is adequate for short cab IRs on desktop).
- **Deferred polish:** a test pinning the symmetric soft-knee as intentional; `eq::Svf` flush-active-channels micro-opt (only if it ever matters).

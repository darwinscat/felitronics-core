<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->
<!-- Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa. -->

# felitronics-core — DSP architecture (living ADR)

**Status:** draft / in progress · **Owner:** Darwin's Cat (Felitronics line) · **Started:** 2026-06.

`felitronics-core` is a **shared, framework-agnostic, JUCE-free DSP core** for the whole product
family. One set of battle-tested, real-time-safe DSP primitives that every product builds on:
plugins (TabbyEQ, a future compressor / true-peak limiter), the **guitar amp** plugin (already
reusing the EQ core), a future **browser / WASM** build, and eventually **hardware with a processor
inside** (SoC / DSP pedal).

The model already works: TabbyEQ's `teq::` EQ core is JUCE-free and is **already vendored into the
guitar amp plugin** for its tone controls. This document promotes that one-off reuse into a
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

---

## 2. Portability laws (non-negotiable — WASM + embedded are first-class targets)

These are what make the core actually portable. Violating one silently breaks a target.

1. **No threads in the core.** Threading is the *adapter's* job. The core is called synchronously.
   *(Why: WASM threads need SharedArrayBuffer + COOP/COEP; bare-metal may have none. Precedent:
   TabbyEQ's background linear-phase FIR builder lives in the JUCE adapter, not in `teq::`.)*
2. **No allocation / lock / IO / syscall / throw in `process()`** (the existing RT rule). Preallocate
   in a `prepare()` step. Ideally no heap at all in hot classes (fixed-size state) so bare-MCU works.
3. **Float in the hot path; `double` only in offline coefficient design.** No FPU → `double` is
   expensive. (teq already designs coeffs in `double`, processes samples in `float`.) Plan to
   *template the sample type* once a real embedded target lands.
4. **Dependencies behind a seam, never hard-wired.** Heavy primitives (FFT, neural inference) are
   reached through a thin interface so each platform plugs its own impl: JUCE adapter →
   `juce::dsp::FFT`; WASM / embedded → pffft / kissfft / CMSIS-DSP. (Precedent: TabbyEQ's
   linear-phase convolution uses `juce::dsp::FFT` *in the adapter*; the `teq::` core stays clean.)
5. **Configurable sizes.** Caps like `kMaxChannels` (16 on desktop) shrink to mono/stereo on
   hardware via a single constant. Optionally a fixed max block size.
6. **No exceptions / RTTI / OS / filesystem / locale** in core code. No global mutable state.
7. **C++ subset that all toolchains accept** (C++20 desktop; keep an eye on what Emscripten and the
   embedded toolchain support — avoid exotic stdlib if an embedded target needs `-fno-exceptions`).

---

## 3. Module layout

One umbrella, **independent modules**, each separately includable so a consumer pulls only what it
needs (an EQ plugin must NOT drag in the neural runtime). Each module = a CMake `INTERFACE` target
(`felitronics::<module>`), header-only + JUCE-free by default, with its own JUCE-free self-tests.

| Module          | What                                                                 | Deps / notes |
|-----------------|---------------------------------------------------------------------|--------------|
| `core`          | Math, `Smoother`, denormal flush, `kMaxChannels` (SSOT), fixed-size ring / lock-free FIFO, base types | zero deps, the shared base everyone uses |
| `eq`            | matched biquads (Vicanek) + Cytomic SVF + `EqBand` + `EqEngine`      | = today's `teq::` (becomes `eq`) |
| `dynamics`      | `EnvelopeFollower` (peak/RMS, attack/release) + `GainComputer` (threshold/ratio/knee/range, downward+upward) | zero deps; **first NEW module** |
| `analysis`      | spectrum tap (the existing `SpectrumTap`), correlation, LUFS/loudness | FFT **behind a seam** |
| `convolution`   | partitioned (uniform/non-uniform) FFT convolution                   | FFT **behind a seam** |
| `oversampling`  | polyphase up/down-sampling                                           | for true-peak + nonlinear stages |
| `limiter`       | true-peak limiter                                                    | uses `oversampling` |
| `neural`        | amp / NAM-style inference                                            | **heavy, isolated**: RTNeural / Eigen; the guitar amp's domain |

**The FFT seam** (`core/Fft`): an abstract interface (or concept) so `analysis` / `convolution` /
`limiter` never hard-code an engine. Adapters provide the impl.

---

## 4. What stays OUT of the core (lives in adapters / consumers)

- Threading, scheduling, background rebuilds.
- Parameter systems (JUCE APVTS, etc.), state save/load, GUI.
- Host glue (AudioProcessor, CLAP, AU), file/IR loading, the concrete FFT engine.
- **Cross-module "glue."** Example: TabbyEQ's *dynamic EQ* = `eq` band + `dynamics` detector fed by
  the band's band-pass, applying the computed gain as an SVF gain-delta. That composition is
  TabbyEQ-specific and lives in TabbyEQ, **not** in `eq` or `dynamics`. A standalone compressor uses
  the same `dynamics` broadband. The core gives primitives; products compose them.

---

## 5. Build & repo strategy

- Each module is a header-only `INTERFACE` library; the umbrella `CMakeLists.txt` exposes
  `felitronics::core`, `felitronics::eq`, … Consumers link only what they use.
- **Consumption:** pinned tag via CMake `FetchContent` (exactly how the products already pin JUCE),
  replacing the current *copy-the-folder* vendoring of `teq/`.
- **Portability CI (later):** build `core` + each light module under Emscripten and under an embedded
  toolchain (arm-none-eabi) so a thread/alloc/dep creep is caught immediately.
- **License:** AGPL-3.0-or-later, SPDX header on every file. Each heavy module records its
  third-party deps + AGPL-compatibility in `THIRD_PARTY_NOTICES.md` (RTNeural MIT/BSD = OK; pffft
  BSD = OK; watch ONNX/others).

---

## 6. Migration plan (keep every product green throughout)

1. **Establish conventions here** (this doc) before adding code.
2. **Land `dynamics` as the first module of the new structure** (it's the next TabbyEQ feature) —
   developed against these rules from day one. See TabbyEQ's upcoming `docs/DYNAMICS.md`.
3. **Move `teq/` → `felitronics-core/eq`**, keeping a thin `teq::` compat alias during transition;
   TabbyEQ + the guitar amp switch from vendored-copy to `FetchContent` pin one at a time.
4. **Extract heavier modules** (`convolution`, `limiter`, `neural`) as they prove out in products.

---

## 7. Consumers (each product fills in its own section)

### TabbyEQ (premium track EQ) — desktop plugin
- Uses **`eq`** (full: 24 bands, matched + SVF, M/S dual-mode).
- Will use **`dynamics`** (per-band dynamic EQ + de-esser, per Mid/Side lane), **`analysis`**
  (spectrum + correlation), **`convolution`** (linear-phase mode — already implemented in the JUCE
  adapter), possibly **`limiter`** on the output.
- Already split JUCE-free core (`teq::`) + thin JUCE adapter (`src/`). Formats VST3/AU/Standalone/CLAP.
- Future target: a **WASM** demo.
- **Note:** dynamics is time-varying → incompatible with the static linear-phase FIR. In Linear mode,
  dynamics is disabled / forced to Natural. (Bake into `docs/DYNAMICS.md`.)

### Guitar amp plugin — *TO BE FILLED BY THE GUITAR-AMP SESSION*
> The session working in the guitar amp repo knows this product best. Please fill in:
> - Which parts of `teq::` you already vendor, and how (copy? which files? what for — tone stack?).
> - The amp-specific DSP you have: NAM / neural amp model (which runtime — RTNeural/Eigen/other?),
>   IR cabinet convolution (which FFT?), tone stack, noise gate, boost/drive, oversampling, etc.
> - What you'd want from the shared core (`dynamics`? `convolution`? `oversampling`? `neural`?), and
>   what amp-specific code should stay in the product vs become a reusable module.
> - Target platforms (desktop only? WASM? **hardware / embedded** — which SoC/DSP, fixed-point?,
>   no-heap?, mono/stereo only?). These drive the portability laws above.
> - Any constraints/gotchas the core design must respect for the amp.

### Future products
- Standalone **compressor** (Felitronics) — `dynamics` broadband + optional sidechain `eq`.
- **True-peak limiter** — `oversampling` + `limiter`.
- **Browser / WASM** demos of any of the above.
- **Hardware** pedal/box with a processor inside.

---

## 8. Open decisions (to confirm)

- Sample-type templating: do it now, or when the first embedded target is real? (Lean: when real.)
- FFT seam API shape (interface vs C++20 concept vs a tiny vtable).
- Single repo (this) vs multi-repo per module. (Lean: single repo, multiple targets.)
- Naming: `felitronics::` vs keep `teq::`/`dyn::` per-module namespaces under the umbrella.

<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->
<!-- Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa. -->

# DSP architecture — third review + conclusions

**By:** the guitar-amp (OrbitCab) session · **Date:** 2026-06-29 · companion to `DSP-ARCHITECTURE.md`
(this is review notes, not the ADR — fold the accepted items into §2/§8/§9 of the ADR when ready).

## Verdict
**Proceed-with-changes.** The keystone call (FFT-seam + zero-latency convolution FIRST) is right; the
module decomposition and the new rigor (versioning-as-behaviour, neural-as-build-time-backend, hybrid
build) are sound. Two things to fix: **one real correctness hole (WASM denormals)** and **scope creep**
(building embedded-grade CI before an embedded product pays for it).

## Third reviewer — agy / Gemini 3.1 Pro (independent; did NOT ratify §9)
Converged with the first panel (Codex, DeepSeek) on the sound calls — FFT-first, neural=build-time
backend, versioning-as-behaviour. Found one new hole + sharpened two:

1. **WASM denormals / Law 8 is broken (NEW, correct).** Law 8 says "the adapter sets `ScopedFlushToZero`;
   the core assumes FTZ but never toggles it." But **WASM sandboxes have no access to the FP control
   register (MXCSR); many embedded ARMs likewise**. So on the `wasm-audio` tier the core's FTZ
   assumption fails → feedback filters (Cytomic SVF, envelope followers) decay into denormals on silence
   → 10–100× CPU spike → AudioWorklet glitch. "Adapter sets FTZ, core assumes it" doesn't hold where the
   product actually wants to go.
2. **Fixed-point `bare-mcu` is a fantasy under a type alias.** A `using Sample = float` → `q31_t` is not a
   flag-flip: a float Cytomic SVF / NAM-Eigen needs different topology, overflow scaling, headroom, coeff
   quantization. `bare-mcu` (fixed-point) is a **separate codebase**, not the same code with an alias.
3. **`WHOLE_ARCHIVE` (NAM self-registration) disables LTO/dead-code-elim** → blows the size budget of the
   small tiers. (Another reason `neural` must be a per-tier *build-time backend*, not one runtime.)

On decisions: FFT-first **sound but must be compile-time** (template/concept, *not* `virtual` in the
convolver hot path — vtable kills inlining/vectorization); 4 tiers + day-one CI **over-engineered**
(active scope = desktop / wasm-audio / embedded-fpu; drop `bare-mcu` from CI until a real product);
neural-build-time-backend / versioning-as-behaviour / portable-denormals-intent **sound**.

## My assessment (guitar-amp session)
- **Law 8 (denormals): agy is right — VERIFIED in code, with a severity refinement.** `teq` ALREADY
  flushes denormals **in software** every block: `Biquad::flushDenormals()` / `Svf::flushDenormals()`
  zap state to exact zero when `|state| < 1e-15f` (MatchedBiquad.h:426, Svf.h:81), called per block by
  `EqBand::flushState()` (EqBand.h:268,317). `1e-15` is **far above** the subnormal range (~1.2e-38) and
  a block is too short to re-traverse ~23 orders of magnitude → the state never reaches subnormal across
  blocks, so **`eq` is denormal-safe WITHOUT hardware FTZ** (teq's own comment: "host should *also* set
  FTZ/DAZ" = belt-and-suspenders, not a requirement). So Law 8 ("core … assumes FTZ … never toggles")
  **contradicts the actual `eq` implementation**, and agy's "10–100× *sustained* spike" is **overstated
  for teq** (already guarded) — but **valid for UNGUARDED feedback filters**, i.e. the new `dynamics`
  followers if they don't adopt the same flush. Correct law: *every feedback module flushes denormal
  state in software (teq's `<1e-15 → 0` per-block pattern); hardware FTZ is a desktop optimization,
  never a correctness crutch.*
- **Pushback on agy's fix:** he suggested `-ffast-math` — **reject it.** It breaks IEEE NaN/inf semantics
  (OrbitCab's tests explicitly assert finiteness). The safe fix is the **software flush / anti-denormal
  bias in the kernel**, which `eq` already does and the new modules (`dynamics` followers) must adopt.
- **Fixed-point + WHOLE_ARCHIVE:** both correct; both reinforce existing positions (tier over-scope; neural
  per-tier backend). Nothing to re-litigate, just record honestly.
- **FFT-seam compile-time:** agree; panel left it "TBD in the spike", agy says decide it — pick
  template/concept. (Minor: FFT is called once per partition block, not per sample, so a vtable is less
  catastrophic than stated — but compile-time is still the right default for the inner kernels.)

## Verification (re-checked against code — not taken on the reviewers' word)
| Claim | Status |
|---|---|
| teq flushes denormals in software | **VERIFIED** — `Biquad/Svf::flushDenormals()` zap `\|state\| < 1e-15f → 0` per block (MatchedBiquad.h:426 / Svf.h:81 / EqBand.h:268,317). **Refined:** this *prevents* subnormals without FTZ → `eq` is already WASM-safe; the spike severity applies to *unguarded* feedback filters, not teq. |
| WASM has no FTZ/DAZ (FP-control access) | **Spec fact** — WASM is strict IEEE-754 with no FP-environment control; not runnable here but well-established. The risk is real for any feedback kernel that does NOT software-flush. |
| `WHOLE_ARCHIVE` "disables LTO" | **Refined** — it prevents *dead-code elimination* of NAM's self-registering objects (binary bloat), not LTO as such; still conflicts with small-tier size budgets. Requirement verified in OrbitCab `CMakeLists.txt`. |
| `juce::dsp::Convolution` zero-latency in OrbitCab | **VERIFIED** — `Convolver.h`; `updateLatency()` sums only the NAM resampler, convolution contributes 0. |
| OrbitCab tests assert finiteness (⇒ reject `-ffast-math`) | **VERIFIED** — `tests/AmpEqTests.cpp` / `AmpStageTests.cpp` `anyBad()` NaN/inf checks. |
| fixed-point `bare-mcu` ≠ float `using Sample` alias | **Conceptual** — sound (fixed-point DSP = different topology/scaling/quantization), not code-verifiable; no embedded target exists to test against. |
| FFT seam must be compile-time, not `virtual` in the hot loop | **Mostly sound** — but FFT is called ~once per partition block (not per sample), so a vtable is less catastrophic than stated; compile-time still the right default for the inner kernels. |

## Recommended ADR edits (when applied)
1. **§2 Law 8 — rewrite:** core does **software** denormal flushing in feedback kernels (teq already does);
   FTZ is a desktop optimization, not assumed for correctness; `wasm-audio` relies on the software path.
2. **§2 / §8 — caveat:** `embedded-fpu` (float) is reachable via templating; **`bare-mcu` (fixed-point) is a
   separate codebase**, not a `using Sample` flag-flip — don't promise it from the same code.
3. **§8 — FFT seam resolves to compile-time** (template/concept), not runtime `virtual` in the hot path.
4. **§2 — tiers gated vs aspirational:** CI-enforce the tiers actually built (desktop, wasm-audio,
   embedded-fpu); keep `bare-mcu` / `arm-none-eabi` documented-but-not-gated until a product funds it.
5. **§9 — record agy (Gemini 3.1 Pro) as the third reviewer** + its dissent on Law 8 (the denormal hole).

## OrbitCab migration takeaways
- **De-JUCE-ing the cab path (`juce::dsp::Convolution` + JUCE SVF) is the guitar plugin's main cost** —
  exactly what the FFT-seam-first spike de-risks. Until then OrbitCab stays desktop-only on `juce::dsp`
  and consumes only the already-JUCE-free pieces (`eq` = teq; later `neural`, `analysis`).
- **teq's per-block software denormal flush is the reference** for how the `eq` module must behave on the
  `wasm-audio` tier (it already does the right thing — the ADR just needs to say so).

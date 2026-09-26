<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->
<!-- Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. -->

# felitronics-core — DSP architecture (living ADR)

**Status:** draft / in progress · **Owner:** Darwin's Cat (Felitronics line) · **Started:** 2026-06.

> **Convolution performance — two different debts, and neither is the one this box used to describe.**
>
> **✅ CLOSED in v0.4.0 (PRs #23–#26) — the scalar-FFT / block-explosion debt.** It was never "the default
> backend is scalar"; it was that scalar was the *only* backend and that the convolver's partition was tied
> to the host block. Both are gone: a SIMD **pffft** `RealFftBackend` ships as the optional compiled module
> `felitronics::fftpffft` (`-DFELITRONICS_WITH_PFFFT=ON`) behind the existing template seam, and the
> partition is decoupled from the host block. Long-convolution cost is now **block-INDEPENDENT** — a
> 131072-tap linear-phase EQ costs ~2.0 %RT at *every* host block, against ~39 %RT at block 8192 before.
> (`DefaultRealFft` is still `ScalarRadix2Real`, `core/Fft.h:233` — that is the seam's zero-dependency
> default, which is the point of a seam.) [`PERF-SCALAR-FFT-BOTTLENECK.md`](PERF-SCALAR-FFT-BOTTLENECK.md).
>
> **🟢 ACCEPTANCE MET — convolver CPU vs `juce::dsp::Convolution`.** The successor item, and now delivered
> in full: `MatrixConvolverNupc` (non-uniform / Gardner — a 128-sample time-domain head plus geometrically
> growing overlap-save FFT tail stages) is **complete**, with every topology (mono / LRDiag / MSDiag / Full)
> and the click-free 2-slot smoothstep warm crossfade, and `lineareq`'s linear- and natural-phase EQs
> convolve on it. Flat ~0.8 %RT at every block, **3–9× cheaper than JUCE at the 64–128 blocks live rigs
> run**, 2.4× cheaper than v0.4.0, true sample-zero-latency, NULL-verified.
>
> **🟡 What is actually still open** is one named thing, not the item as a whole: the **worst-buffer spike**
> — all stages' FFTs coincide every `lcm = B_max` samples (7.8 % @ block 64 for `B_max=4096`; `B_max=2048`
> is the default because it halves the spike at the same mean). Time-distributing the large FFTs over their
> deadline is the RT-hardening that closes it. Separately and permanently NOT a goal: JUCE stays cheaper on
> the *mean* at large oracle-tuned blocks — that is a theorem, the price of near-field zero-latency
> coverage, not a gap to chase. [`PERF-CONVOLVER-JUCE-GAP.md`](PERF-CONVOLVER-JUCE-GAP.md) ·
> [`PERF-NUPC-VS-JUCE.md`](PERF-NUPC-VS-JUCE.md).

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
2. **WASM** (browser demos / lightweight web) — **secondary, and GATED in CI** since the `wasm-audio`
   job landed. Read §2 for what that gate does and does not prove.
3. **Embedded / hardware** — **far-future, low priority.**

We keep the core JUCE-free *so tiers 2–3 stay open*, but we **don't pay their tax early** (fixed-point,
no-heap, embedded CI) — a tier earns its CI gate and its constraints only when a product funds it.

---

## 2. Portability laws (non-negotiable — keep WASM + embedded reachable)

These are what make the core actually portable. Violating one silently breaks a target.

**Target tiers (each module declares which it supports).** Generic "WASM + embedded" is too coarse —
name the profiles. Per the **platform priority (§1)**, **CI gates `desktop` (primary) and, since the
`wasm-audio` job, that tier too; both embedded tiers stay documented-but-not-gated** until a hardware
product exists (no point paying embedded's CI tax early).

**What the `wasm-audio` gate proves, exactly.** It builds every default module for wasm32 with
`-fno-exceptions -fno-rtti`, runs the whole self-test suite in node, and audits every emitted `.wasm`
for shared memory or thread-shaped imports. So: an exception or RTTI use is a **compile** error; an
incompatible dependency is a compile/link error; an allocation in `process()` is a **runtime** failure
(emscripten is libc++, so the alloc counter — `test_support/alloc_counter.h`, asserted through
`felitronics_test.h`'s `okNoAlloc` — actually enforces there, unlike the libstdc++ rows). **A thread on a live code path is a LINK error** — but only because the tier links
`--wrap=pthread_create`: on its own, emscripten without `-pthread` links pthread *stubs* returning
`ENOTSUP`, so `std::thread` compiles, links and merely aborts at runtime (measured, emsdk 6.0.9). What
the wrap still cannot see: a thread in code never emitted (an unused inline, an uninstantiated template),
and `std::mutex`, whose lock lives inside `libc++-noexcept.a` and whose stub silently succeeds — a lint
over module headers covers the latter. **Law 8 is not covered at all**: a denormal stall is a property of
the CPU at runtime, invisible to any build. Full write-up:
[`WASM-AUDIO-TIER.md`](WASM-AUDIO-TIER.md).
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
   a plugin builds/loads models + IRs on a background/message thread and atomic-swaps them into the
   live pointer — the DSP core never spawns or blocks.)*
2. **No allocation / lock / IO / syscall / throw in `process()`** (the existing RT rule). Preallocate
   in a `prepare()` step. Ideally no heap at all in hot classes (fixed-size state) so bare-MCU works.
3. **Float in the hot path; `double` only in offline coefficient design.** No FPU → `double` is
   expensive. (teq already designs coeffs in `double`, processes samples in `float`.) Introduce a
   `using Sample = float` alias and keep raw `float*` out of public signatures **now**, so full
   sample-type templating later is a flag-flip, not a fork-rewrite. *(Carve-out: coefficient recompute
   on the audio thread, and meter/LUFS/true-peak accumulators, may legitimately use `double` — specify
   per module; the ban is on gratuitous `double` in the sample loop.)* **`double` is the ceiling — `long
   double` is not a third option anywhere near a tier boundary, see law 9.**
4. **Dependencies behind a seam, never hard-wired.** Heavy primitives (FFT, neural inference) are
   reached through a thin interface so each platform plugs its own impl: JUCE adapter →
   `juce::dsp::FFT`; WASM / embedded → pffft / kissfft / CMSIS-DSP. (Precedent: both plugins run the
   spectrum-analyzer `juce::dsp::FFT` *in the adapter/UI* — TabbyEQ `ui/EqCurveDisplay.h:111`, OrbitCab
   `ui/SpectrumAnalyser.h` — while the `teq::`/`cab::` cores stay FFT-free; OrbitCab also hides
   `juce::dsp::Convolution` behind `cab::Convolver`'s raw-float seam. NB: TabbyEQ now **ships** linear- and
   mixed-phase EQ via `felitronics::lineareq` over `felitronics::convolution` — the `teq::` core itself stays FFT-seam-based.)
5. **Configurable sizes.** Not just `kMaxChannels` (16 desktop → mono/stereo on hardware) but also
   `kMaxBlockSize` / `kMaxBands` / `kMaxPartitions` / `kMaxIrSamples`. Set via **template param /
   policy + a per-tier CMake preset**, not a fragile pre-include `#define`. Each module reports its
   memory footprint.
6. **No exceptions / RTTI / OS / filesystem / locale** in core code. No global mutable state (watch
   header-only function-local statics).
7. **C++ subset that all toolchains accept** (C++20 desktop; keep an eye on what Emscripten and the
   embedded toolchain support).
8. **Denormals: every feedback kernel flushes in SOFTWARE; hardware FTZ is a desktop optimization,
   never a correctness crutch.** Every recursive kernel in the repo now has a written verdict —
   including the clean ones and why — in [`LAW8-AUDIT.md`](LAW8-AUDIT.md); the mechanism is a
   **fixed-point band**, not "it eventually underflows". *(The subtle hole: WASM and many
   embedded ARMs expose no FP-control register, so an "adapter sets FTZ, core assumes it" rule silently
   fails on the `wasm-audio` tier → feedback filters decay into subnormals on silence → 10–100× CPU
   spike.)* `teq` already does the right thing: `Biquad/Svf::flushDenormals()` zaps `|state| < 1e-15f`
   to exact zero every block, so the state never reaches the subnormal range — `eq` is denormal-safe
   **without** hardware FTZ. **The new `dynamics` envelope/release followers (and any future feedback
   module) MUST adopt the same per-block software flush.** The core never sets a global FTZ/DAZ mode; an
   adapter MAY set FTZ on desktop as a bonus. **Do NOT use `-ffast-math`** (it breaks the NaN/inf
   semantics the tests assert). Every kernel also ships a scalar fallback + a scalar↔SIMD parity test.

   **8a. THE FLUSH IS CLOCKED BY AUDIO TIME, NEVER BY THE `process()` CALL — `core::StateGrid`.** The
   sentence above says *"every block"*, and for a decade of host callbacks that read as a cadence. It is
   not one: the end of `process()` is wherever the CALLER chose to cut the stream, so the flush is an
   event whose position the core does not control. That fails in both directions at once. **The output
   becomes a function of the slicing** — measured on `eq::EqEngine`, 38522 of 40000 tail samples differ
   between a whole-file call and one-sample calls of the same programme, and on `stereo::MonoBass`
   37180 of 40000. **And on a call that spans a whole file the flush never fires inside it at all**, so
   the stall this law exists to prevent happens in full (37678 of 38000 tail samples subnormal), and a
   single non-finite input sample is never healed either — one +Inf poisoned 479900 of the next 480000
   samples of a whole-file render against 1 sample when the same stream was fed one sample at a time.
   A whole-file call is not exotic: `Compressor`, `EqEngine`, `Dither` and `TruePeakLimiter` all invite
   one in their headers.

   So: **a kernel drives its periodic maintenance from `core::StateGrid`** — a phase counter over AUDIO
   samples, period `kPeriod = 64`, re-anchored by `reset()` — and splits its per-sample work at those
   boundaries. Two kernels reached this conclusion on their own before it was a rule and are the
   reference implementations: `analysis::LoudnessMeter` (its 10 ms sub-hop, reasoned out in
   [`LAW8-KWEIGHTING.md`](LAW8-KWEIGHTING.md)) and `saturation::Saturator` (per sample). The claim a
   gridded kernel may then make is **bit-identical output under arbitrary re-slicing from the same
   `reset()`, for an identical event timeline and while the state stays FINITE** — nothing about events
   that ARRIVE per call (a parameter write, a bypass toggle, a channel-count change), which are the
   caller's own timeline, and nothing once a filter has overflowed: poison recovery deliberately keeps
   the call's clock (below), so from the first non-finite state onward the output is call-dependent
   again. That last clause is not hypothetical — a swept band-pass at 20 kHz Q 2 fed `{3e38f, 0}` emits
   a finite first sample and overflows its integrator, after which two one-sample calls give 0 on the
   second sample and one two-sample call gives NaN.

   Three things this rule does NOT say. **(a) The poison half keeps its own clock.** A NaN's only quality
   is how soon it goes, so the `isfinite` half of a flush runs at the END OF EVERY CALL as well — it
   cannot change a bit while the state is finite, so it costs the invariance claim nothing, and a host
   with a block shorter than a period keeps the immediate recovery it had (`eq::Biquad::healPoison()`).
   **(b) `kPeriod` is chosen, not derived.** Any period makes a stall unreachable, because a state under
   the 1e-15 threshold is 23 decades above the subnormal floor and cannot outlive one period; 64 is the
   smallest block a live rig runs and a power of two, so every host block that is a multiple of it is a
   whole number of periods and no segment is ever split. What the period does set is the one-time subnormal exposure per silence event (measured worst
   run for a real pole from 1e-15: 8 samples at 32, 16 at 64, 180 at 256, 8116 at 8192) and the poison
   recovery bound. Do not argue a period from pole radius: a DF2T biquad with `a1 = −65/128`,
   `a2 = 9/128` — poles at radius 0.265 — walks from 1e-15 to the exact nonzero fixed point
   `(z1, z2) = (2⁻¹⁴⁹, −0)` in 53 silent updates and stays there for ever. **(c) An oversampled kernel
   counts its OWN clock**: 64 base-rate frames is 512 recursive updates at 8×.
9. **No `long double` in anything that crosses a tier boundary.** It is not a type with a defined
   precision — it is whatever the target's ABI happens to say, and the answers disagree in *size* as well
   as in precision. Measured, one source file, three toolchains:

   | tier | `sizeof(long double)` | mantissa bits | what it actually is |
   |---|---:|---:|---|
   | arm64 macOS (Apple clang 21.0.0) | 8 | 53 | **plain `double`** — no extra precision at all |
   | x86-64 Linux (gcc 14.2.0) | 16 | 64 | 80-bit x87 extended, padded to 16 |
   | wasm32 (emsdk 6.0.9) | 16 | 113 | IEEE binary128 quad, software-emulated |

   Two traps live in that table. **The two 16-byte answers are different formats** — `sizeof` agreeing
   proves nothing, and a size assertion would pass while the arithmetic diverged. And **the dev machine is
   the least precise of the three**: on arm64 macOS `1.0L + 2⁻⁶⁰ == 1.0L`, so an accumulator written and
   tested *here* because it "needed the headroom" silently has none, while the same line carries 64 bits on
   the Linux CI row and 113 in the browser. (arm64 Linux is a fourth answer — AAPCS64 mandates binary128 —
   documented, not measured here.) None of this is a toolchain bug; it is what the type is.
   **Use `double`, or compensated (Kahan / Neumaier) summation**, which is deterministic on every tier
   because it never asks for anything but IEEE `double` arithmetic. **There is no exception any more.**
   The law used to sanction one, the `correlation` accumulator in `tools/fcore_measure.cpp` (now
   felitronics-mastering-core's); the waveform-and-stereo port found it
   was a third definition of a number the page and the stereo band already defined, and it became their
   binary64 formula (`analysis::StereoSums`, felitronics-mastering-core). The artifact gate below lost its one named exclusion with it. **This law IS enforced**, in the two halves the
   `wasm-audio` job pairs everywhere: a text lint over `modules/*/include` and `modules/*/src` that lexes
   before it matches (the words appear in prose constantly, including in this paragraph) and also catches
   an L-suffixed float literal, which is a `long double` that never names itself; and an **artifact** gate,
   because on wasm32 the type is software binary128, so every operation on it is a compiler-rt libcall
   standing in the object file as an undefined symbol. The text half sees code that is never emitted; the
   artifact half cannot be fooled by a comment, a macro or an alias. Both carry negative controls.
10. **FP contraction is STATED, never inherited: the build sets `-ffp-contract=on`.** `a*b + c` may be
   fused into one FMA with a single rounding instead of two — a *different* number, slightly more
   accurate, and the toolchains disagree about when they are allowed to do it: **clang defaults to `on`**
   (fusion within one expression — the C/C++ standard's own rule), **GCC defaults to `fast`** (fusion
   *across statements*, beyond what the standard permits), and **MSVC's `/fp:precise` does not contract
   at all**. Three defaults, three sets of numbers, from one source file.
   *(Not academic. It broke two exactness claims the day an arm64 Linux row first reached CI, on
   untouched `main`: a tube stage's "silence in ⇒ output identically 0" read **1.26e-08**, and the
   multi-res fast pane missed its budget by **0.005841 dB**. Both are cancellations — TT1 wants
   `g(0+vb) − g(−0+vb)` to vanish exactly — and a cancellation stops cancelling the moment one half is
   computed with one rounding and the other with two. Nothing could see it before: baseline x86-64 has
   no FMA instruction to contract with at all, and the only arm64 rows were also the only Apple rows,
   where clang's default already behaved.)*
   **`on`, not `off`.** It fixes both, and it costs the shipping tier nothing because it is already
   clang's default — macOS builds exactly as before, so no perf change and no golden/NULL churn. `off`
   fixes them too but bans FMA outright, measured at **+12 % on the `Svf` sample loop** (16.25 → 18.24 ms,
   arm64 macOS) for no benefit anyone has demonstrated. On the shipping convolver both are within noise
   (~0.58 %RT either way — its hot path is pffft's intrinsics, which no contraction flag reaches).
   NB **GCC ≤ 13 implements `on` as `off`** (no `fmadd` emitted at all); GCC 14 emits the standard
   behaviour — one `fmadd` for the single-expression form, none for the cross-statement one. Verified on
   both. **This law is enforced, and it is enforced differently from the two before it**: not by a lint
   and not by a scan, but because every non-MSVC row compiles with the flag — remove it and the arm64
   Linux row goes red on its own. (It said "unlike 8 and 9" when it was written, hours before law 9 got
   its own two gates. Law 8 remains the odd one out: a denormal stall is a property of the CPU at
   runtime and no build can see it.)
   **THE LAW HAS EXACTLY ONE LOCAL OVERRIDE, and it is a header that turns contraction OFF for itself:**
   `core::firDot` (`modules/core/include/felitronics/core/PolyphaseFir.h`, P56), the single polyphase-FIR
   inner product that `oversampling::PolyphaseOversampler`, `oversampling::CascadeOversampler` (P31) and
   `analysis::TruePeakMeter` all run. (The cascade has the one sum outside it — its halfband decimator's
   centre tap — and adds a value halved and STORED a sample earlier, so no contraction can join the two.) The law
   is right about the tree and wrong about that one loop: `acc += a*b` is the contractible form, arm64 fuses
   it and baseline x86-64 cannot, so contraction there does not make one number better — it makes five rows
   disagree. With the order nailed down and contraction off locally, they do not: one measured constant on
   `win` (MSVC 19.44), `deb` (gcc 14.2), `mac` (Apple clang 14), `docker --platform linux/arm64` (gcc 14.4)
   and `wasm` (emsdk 6.0.9). A header cannot use a flag, so it takes three pragmas, all of them measured on
   the toolchain rather than read out of a manual: gcc ignores `#pragma STDC FP_CONTRACT` in C++ and needs
   `#pragma GCC optimize("fp-contract=off")` (which also blocks inlining — **1.4 %** of a mastering render,
   measured); clang needs `#pragma clang fp contract(off)` and must NOT be given
   `#pragma float_control(precise, on)`, which turns contraction back ON; MSVC takes `#pragma fp_contract(off)`,
   and its `float_control(push)`/`(pop)` really does save it. The override reaches only as far as the language
   does — **a clang build with `-ffp-contract=fast` or `-ffast-math` defeats every pragma**, so the property is
   GATED (`felitronics_core_polyphasefir_tests` pins the cross-row bits and fails on a row that fused) rather
   than asserted. What stays outside it, because it is runtime state and not code: FTZ/DAZ (wasm cannot flush
   at all, so a host that flushes splits it from every native row — with NORMAL inputs, since a normal times a
   normal can land subnormal), the rounding mode, and NaN sign/payload.
   **AND A PIN THAT IS NOT AN OVERRIDE.** Where a number is gated against a reference that never contracts
   (JavaScript), the right answer is the unfused one, and a `volatile` store of each product before it is added
   holds under every contraction mode, gcc's `fast` included, with no pragma. The worked case — the stereo band's
   `analysis::StereoSums`, the width it read fused on arm64 and every cell the pin was measured in — moved with
   the offline analyzers to felitronics-mastering-core (`docs/LAW10-STEREOSUMS.md` there).

11. **THE CALL IS A REQUEST AGAINST A PREPARED CAPACITY, AND A REQUEST THAT CANNOT BE HONOURED IN FULL
   IS REFUSED AS A WHOLE — `[[nodiscard]] bool process(...)`.** `prepare(sampleRate, maxBlock,
   maxChannels)` states what the object was built for; every argument it takes is binding, and a module
   that ignores one is lying about its contract. What a caller may then hand `process()` split into four
   questions, and the core used to answer each of them four different ways — the census that opened this
   law found **three answers for the width, four for the length, and a case that fits none of them**:

   **(a) THE LENGTH `n` IS A CAPACITY, NOT A LIMIT.** `maxBlock` sizes scratch. `process()` chunks
   internally, so any `n >= 0` is processed IN FULL, and the chunked pass is **bit-identical to the
   caller having chunked it itself AT THE SAME BOUNDARIES** — that is the invariant, and it is what makes
   the rule testable. (It is deliberately not "bit-identical under arbitrary re-slicing": that is law 8a's
   claim, it is owned by the grid, and a rate-matched stage cannot make it at all.)
   Never truncate the tail, never refuse a long call. *(An offline caller sizing `maxBlock` to a whole
   file is a normal thing to do — `Compressor`, `EqEngine`, `Dither` and `TruePeakLimiter` all invite one
   in their headers — and refusing it would return the buffer UNTOUCHED, i.e. exactly the
   unlimited-passthrough defect the refusal was meant to prevent.)* Measured before this law:
   `dynamics::NoiseGate` clamped `n` to its curve and let **3840 of 4096** samples out **+89.99 dB**
   louder than the gated ones — 100 % of the construction ceiling, which is `-floorDb` = 90 dB;
   `multiband::MultibandProcessor` dropped the **entire** call.
   **The one exception is a two-phase API whose first phase RETURNS a buffer of `maxBlock`**
   (`NoiseGate::analyse`/`applyGain`, `EqEngine::captureSectionInput`): it cannot chunk, because the
   result must outlive the call. Those refuse — observably — and their fused convenience form chunks.
   **PHASE A OWNS THE CLOCK, PHASE B OWNS NOTHING.** The two phases describe ONE block of audio time, so
   only the first may spend it; a phase B that also advanced would count the same samples twice. And
   phase A owes a LENGTH as well as a buffer — what it actually produced — because phase B bounded by
   the buffer's capacity instead applies a curve from a previous call
   (`NoiseGate::analysedSamples()`, `EqEngine::sectionInputSamples()`).

   **(b) THE WIDTH `nch` IS A LIMIT.** State for a channel that was never prepared cannot be invented,
   so `nch > maxChannels` is refused **as a whole call, before anything moves** — the refused call is
   indistinguishable from one never made. Never a prefix. A prefix is not the safer half-measure it
   looks like: the surplus channels are unprocessed either way, and processing some of them only hides
   the fault while the processed ones acquire a latency and a gain the others do not — a level and comb
   mismatch on fold-down. Measured: `limiter::TruePeakLimiter` prepared for 2 and called with 4 emitted
   the surplus at **+7.02 dB over the ceiling it was told to hold** (unbounded in general — it is the
   caller's own input, untouched); `saturation::Saturator`'s surplus came out **bit-identical to the
   input**, with no saturation at all.

   **(c) A NARROWER CALL IS LEGAL.** `nch < maxChannels` is a supported mode, not an error: it means
   those channels are playing and the rest stopped, and law 11a below says what happens to the state of
   the ones that stopped. A module for which a narrow call is *meaningless* — a matrix convolution needs
   both input planes to compute either output; a stereo stage needs two — refuses it, says so in its
   header, and the refusal is **returned**. Silence is not a refusal: `convolution::MatrixConvolverNupc`
   used to drop a narrow call and write nothing at all, so a caller that pre-zeroed its output buffer
   got digital silence and no way to find out. **A module whose width is EXACT refuses `nch == 0` too,
   and its gap is `reset()`, not a zero-width call** — the clock-only form of 11(d) is not expressible
   for an operator that needs its planes. A composite therefore does not forward a gap to such a
   module: `lineareq::LinearPhaseEq`, `NaturalPhaseEq` and both matrix convolvers are the list here, and
   felitronics-mastering-core's `mastering::MasteringChain` is one downstream.

   **(d) FOR AN ACCEPTED CALL, `n > 0` IS THE TRIGGER FOR BOTH CLOCKS: AUDIO TIME *AND* THE FALLING
   EDGE.** "Accepted" is not decoration: a REFUSED call moves nothing at all, the clock included, so the
   checks of (b) and (c) come first and this clause never applies to a call that failed them. Where a
   width is EXACT (the matrix convolvers) a zero-width call is a refusal, not a clock-only call, and the
   width check settles it before `n == 0` is even looked at. A call carrying
   samples spent that many samples of audio time whether or not any channel was processed, so the law-8a
   grid ADVANCES — `nch == 0, n > 0` is a **clock-only** call, exactly what `stereo::MonoBass` already
   does in bypass (`grid_.skip(n)`). Measured before this law, on the same band and the same stream:
   `eq::EqBand::processBlock(io, 0, N)` advanced its grid and `eq::EqEngine::process(io, 0, N)` dropped
   the call, and after 10240 samples of zero-width calls their glides had diverged by **11.08 dB**.
   **And by the same token every channel at index >= `nch` STOPPED for those samples** — at `nch == 0`,
   all of them — so law 11a's falling edge fires too. These two halves are not independent: saying "time
   passed" and "nobody stopped" in one breath is a contradiction, and it reopens the exact defect P18
   closed. Measured on untouched `main`: `Compressor` prepared for 2 with 5 ms of lookahead, a tone, then
   4800 samples of zero-width calls, then stereo DIGITAL SILENCE — **0.280315 out of the silence
   (-11.05 dBFS), last non-zero at sample 239**, i.e. exactly the 240-sample lookahead line replaying
   audio from before the gap. `n == 0` is the only true no-op: no time, no edge, nothing.
   Negative `n` or `nch` is a malformed call: refused, never clamped into an index.

   **The verdict is RETURNED, because a `void` refusal is the disease this law exists to cure.** Every
   block-level `process()` is `[[nodiscard]] bool`; `true` means "this call was accepted and honoured in
   full". A compile-time diagnostic at the call site is the strongest signal available inside a
   `noexcept` function that may not allocate, lock or throw, it costs nothing at runtime, and it is
   already this repository's idiom for a refusal — `prepare()` has been `[[nodiscard]] bool` since P2.
   `bool` rather than a status enum on purpose: the reason is always visible at the call site (the caller
   knows what it passed), a second refusal idiom for the same concept is how a core ends up with four
   policies again, and `bool` composes — a composite ANDs its stages' verdicts the way `prepare()`
   already does. **Order of the checks is part of the law**, so that one malformed call has one answer:
   malformed (`n < 0 || nch < 0`) → unprepared → `nch > maxChannels` → `n == 0` → `nch == 0` → run.

   **11a. THE FALLING EDGE IS CLOCKED BY `n`, NOT BY `nch`.** P18 gave a channel that stops being fed
   the rule: drop its sample memory, or it replays it on return. What it left ambiguous is when a channel
   counts as having stopped, and `eq::EqBand` answered "only on a call that ran at least one channel"
   (`numSamples > 0 && nc > 0`). That half is wrong, and the number above is what it costs: a stretch of
   zero-width calls is a real gap in the stream, and a stage that treats it as "nothing happened" hands
   the gap's far side a frozen delay line. **A channel stopped for this call iff `n > 0` and its index is
   >= `nch`.** The `nc > 0` guard survives only in its other job — keeping a negative width out of a
   half-open range — and law 11's malformed-call refusal now does that job earlier and better.

   **AND WHERE THE MEMORY CANNOT BE DROPPED, IT IS DRAINED.** "Drop its sample memory" assumes the
   memory is ours to clear. For a stage that owns a black box — a third-party network whose window it
   cannot reach, a resampler leg beside it — the answer is to hand the stopped channel the DIGITAL
   SILENCE it is receiving, through the same code path into the stage's own scratch (`io` need not carry
   that plane, and at `nch == 0` may be null), until its state is provably that of a channel that was
   silent all along, and then to STOP: bounded by the memory the stage can count, each part in its own
   rate, so a permanently mono host pays for one lane, not two. "It drains, and then it stops" has no
   witness in the audio — past the debt the output is zero either way — so such a stage publishes an
   odometer a test can read. A recurrent cell has no flush length; its bound is a heuristic, and says so.

   **AND THE OTHER HALF IS `reset()`, WHICH IS A DIFFERENT OPERATION AND NOT A LONGER ONE.** A lane the
   caller STOPPED handing over is still the same stream; a lane that is PRESENT gets the caller's own
   samples, and a stale window speaking into them is a stream RESTART, which is `reset()`'s job. A drain
   SIMULATES silence; a restart puts the stage back where a freshly prepared one is, and for a
   finite-memory stage the two states coincide. What a restart buys is INDEPENDENCE, and it is exact: two
   instances fed different audio before it answer the next programme with the same bits. It re-anchors
   the audio-time clocks as well, exactly as `eq::EqBand::reset()` re-anchors its `StateGrid`. Where its
   cost is not the block's — a whole drain length — the price is published and the operation made
   IDEMPOTENT: the debt is re-armed only by audio actually being fed. And it flushes what the stage's
   LEDGER of its own memory can see, so the ledger answers for the whole model; what it cannot place is
   charged ONE bounded allowance — never zero, never the face value.

   **AND `prepare()` PERFORMS THAT RESTART TOO, ALWAYS — THE TWO VERBS NAME ONE STATE.** A prepared stage
   holds no audio the caller fed, whether or not the rate or the block size changed. There is
   deliberately no predicate on what changed: a re-prepare at the SAME numbers is the common case.
   **A prepare restates what it COUNTS, not only what it designs:** state kept as a count of host samples
   is recomputed, mapped by its predicate, or rescaled by the rate ratio — whichever its own algebra
   allows. **And a composite owes its consumer the same verb:** its `reset()` reaches every child that
   holds audio, bypassed or not.

   The NAM backend and the `.orbitrig` pack player are the worked case of all of this, with the
   measurements: felitronics-guitar-core `docs/LAW11-NAM.md`.

   **11b. `prepare()` IS BINDING, AND REFUSES WHAT IT CANNOT HONOUR.** An observable refusal in
   `process()` is worth nothing if `prepare()` already lied about the width: `convolution::CabConvolver`
   silently clamped `prepare(..., 4)` to 2, after which `process(io, 4, n)` was a perfectly legal call
   that left planes 2-3 DRY, and fifteen more modules clamped the same way. Every argument `prepare()`
   takes is binding; a value it cannot honour is refused there, the way `Compressor` and
   `TruePeakLimiter` already do — which is what makes law 11(b) reachable at all. And "the defaults are
   a valid configuration" is a claim to CHECK, not to assume: `stereo::MonoBass` looked like one and is
   not — its crossover has no coefficients until `prepare()` runs, so a default-constructed object
   passes the side band it is supposed to fold at **-6.02 dB where a prepared one kills it to -54.22**.

   **A REFUSED `prepare()` ADOPTS NOTHING AND LEAVES THE OBJECT UNUSABLE — IN THAT ORDER: DISARM,
   VALIDATE, WRITE.** Read literally, "writes nothing" and "leaves it unprepared" contradict each other,
   because disarming is itself a write; the order is what reconciles them, and getting it wrong costs a
   defect in either direction. Validate-then-disarm leaves the PREVIOUS preparation standing and still
   answering `process()`; write-then-validate leaves an UNVALIDATED argument in place. So: clear the
   readiness flag on entry, validate every argument, and only then store any of them. Both halves are
   load-bearing, and both were got wrong three times each while this law was being applied. Validating one
   argument, storing it, and then refusing on the next leaves a new WIDTH standing beside an old buffer:
   in `mastering::OfflineRenderer` (felitronics-mastering-core) that was a **heap-buffer-overflow**, a write past
   the scratch region, which ASan caught only because a reviewer went looking for it. And a refusal that returns
   before reaching the inner `prepare()` leaves the object ARMED on its previous build — measured on
   `multiband::MultibandCompressor`: `prepare(2)`, then a refused `prepare(0)`, then `process(io, 2, 64)`
   still returned true and still processed.

   **11c. A PAUSE IS SILENCE.** A clock-only call advances a stage's SHARED, one-per-instance ballistics
   **exactly as `n` samples of digital silence at a live width would**. Not frozen, not reset: the same
   arithmetic the audio path runs, with the detector's input at digital zero. This follows from 11(d)
   rather than adding to it — that clause already says the call is a GAP IN THE STREAM and not a no-op,
   so audio time PASSED, and a detector that stands still through passing time contradicts the same
   sentence that made the grid advance. `reset()` is not the answer either, for the reason P19 separated
   the verbs: `reset()` is a stream RESTART and claims more than the caller said, while
   `clearAudioState()` is a stop. The core used to give **three** answers here — freeze (`Compressor`,
   `DeEsser`, `DynamicEqBand`, `TransientShaper`, `NoiseGate`), reset (`limiter::TruePeakLimiter`), and
   disengage (`dynamiceq::LaneDynamics`) — and the freeze is what it cost: `Compressor` held **-25.311 dB
   of gain reduction through a full second of zero-width calls** where the same second of silence
   releases to -2.076, and dipped the return by **-22.11 dB** on material below its threshold (-24.09 dB
   through a ten-second gap); `DynamicEqBand` held -21.774 against -2.985; `DeEsser` -8.000 against
   -1.104. The loudest was `NoiseGate`, and it is not a shifted envelope but a state machine that never
   fired: a gate that has to CLOSE through a pause stayed wide open, and a -54 dBFS tone on the return
   came out at -54 where silence gates it to -144 — **89.99 dB, 100 % of the construction ceiling**
   (`-floorDb` = 90 dB).

   **THE ADDRESSES ARE THE MECHANISM, NOT A LIST**: shared, one-per-instance ballistics that a zero-width
   call leaves without a clock. Six in this repository — the five above plus `dynamiceq::LaneDynamics`
   (whose lanes now run the control loop at width zero, where the Stereo lane's linked probe over zero
   columns is exactly `+0.0f` and L/R/M/S take the same "this lane stopped" branch they take at width
   one) — and a seventh, `poweramp::PowerAmpStage`, in felitronics-guitar-core. **`LaneDynamics` is width-dependent BY DESIGN and its entry here is narrower than it looks:**
   `laneRuns()` gates L/R/M/S on `nc == 2`, so at width zero only the Stereo lane runs on silence and the
   other four take the same "this lane stopped" branch they take at width ONE — a hard drop of their
   detector, exactly as before. Its answer to a pause is therefore "what this stage does at that width",
   which is the honest reading of the rule for a stage whose topology is a function of the width, and NOT
   "the trajectory does not depend on the width" — that sentence below is about the four stages whose
   detector input is the linked frame, and it is false for these four lanes. A composite forwards the gap rather than swallowing it (`multiband::MultibandProcessor` does,
   per band, exactly once — it used to do it twice for a bypassed band, which was invisible under freeze
   and a double clock under this law). `analysis::LoudnessMeter` had already answered this way on its own:
   at `nch == 0` its sub-hop windows keep sliding as zero-energy hops.

   **WHY IT IS WELL POSED AT WIDTH ZERO.** The linked level of an all-zero frame is exactly `+0.0f` at
   EVERY width, width zero included (`linkAmplitudeImpl` returns 0 for `numChannels <= 0`), so the silent
   trajectory does not depend on how wide the caller's silence was — which is what lets a pause be
   compared against a silence at all. The comparison is at the SAME CALL BOUNDARIES: every stage here
   flushes once per call, so a gap cut into three pieces equals a silence cut into the same three. And it
   is stated after the law-11a falling edge has fired: per-channel memory is dropped exactly as before,
   so the detector meets zeros rather than a ring-down. **An external key is still consumed** — the
   PROGRAMME is what stopped, not the key, and `process(zeros, nch, n, key, nk)` runs the detector on the
   key, so a pause that ignored it would differ from the silence it is defined to equal. **Two contract
   changes follow at width zero and are stated rather than left to be discovered:** a caller that passes a
   key now has it DEREFERENCED on a call that previously read nothing, so the key must be valid for `n`
   samples exactly as at any other width; and a caller that passes a `GainReductionTap` now has it FILLED,
   sample by sample, where a zero-width call used to leave it untouched. Both are what "this call equals
   the same call carrying silence" means, and both have their own test.

   **THE COST IS BOUNDED BY THE BALLISTICS AND BY THE CALL, NOT BY THE PAUSE — and "free past the fixed
   point" is true of the collapsing path only.** Say the whole of it, because the short version is wrong
   for three of the seven: `TransientShaper`, `NoiseGate` and `LaneDynamics` have no dB floor to collapse
   against and run their full per-sample body until they park. Measured at 48 kHz, ONE zero-width call
   covering a full minute (2 880 000 samples): `Compressor` **1.49 ms** (it collapses), `NoiseGate`
   **0.20 ms**, `TransientShaper` **3.65 ms**, `DeEsser` **4.08 ms**, `DynamicEqBand` **6.43 ms** — three
   of which are past a 128-sample callback's 2.67 ms budget — and `LaneDynamics`, the most expensive,
   is further past it again. **Read those as orders of magnitude, not as figures:**
   they are wall-clock timings and they moved by 2x between runs of the same binary on the same machine
   depending on what else was building, which is exactly why the arithmetic claims above are stated in
   samples and these are not. What does not move is the SHAPE: that is a statement about ONE CALL carrying
   a minute, which is an offline pattern. The same minute delivered the way a host delivers it, 128
   samples at a time, costs **a few microseconds** in its worst single call — measured between 0.0009 and
   0.003 ms per stage on an idle machine and 0.014 ms on a loaded one, i.e. two to three orders of
   magnitude inside the budget — because the fixed point is reached in the first calls and every later one
   exits on its first step. An RT caller is safe; an offline caller that hands a whole transport jump as
   one call pays the numbers above, once.
   The silent recurrence is AUTONOMOUS, so it reaches a bitwise fixed point and everything past that point
   is free; and once the detector level
   reaches `core::kGainToDbFloor` the dB conversion returns the same bits for every smaller level, so the
   curve's output is a constant and the per-sample work collapses to one multiply-add. `pow(c, n)` is
   deliberately NOT used: it is a different number from `n` rounded multiplications, and it would buy
   speed with the bit-exactness this law is stated in. Be precise about the horizon, because the obvious
   claim is false: the recurrence does NOT reach zero — it parks on a SUBNORMAL, and `flushDenormals()` is
   what turns that into a real zero — and the settling length is a property of the TIME CONSTANT, measured
   at 48 kHz as **23 609 samples for a 5 ms release, 457 808 for 100 ms, 4 461 677 for 1 s, and not
   reached in 200 000 000 at the coefficient cap.** "A ten-hour gap costs what a ten-second one costs" is
   true only past that horizon.

   **THE LIMITER IS DELIBERATELY OUT.** `limiter::TruePeakLimiter` reaches `nc == 0` through its width
   guard and does a full `reset()` — wrong under any answer (measured: gain reduction -5.08 dB to 0.00 and
   a 48-sample hole on the return), but wrong in its own way: its ballistics are not exponential but a
   sliding lookahead window with a running maximum, where "a pause of `n`" means "`n` zeros entered the
   window and it shifted", and its per-channel state is not detector history but AUDIO THAT HAS NOT BEEN
   EMITTED YET, which cannot be dropped without a hole by definition. That is **P29**, and it will state
   its own rule as a window shift, with the oversampler's phase to prove as well.

   **11d. MEMORY THAT CANNOT BE HAD IS NOT A REFUSAL: EXHAUSTION IS FATAL, AND THE CORE PUBLISHES ITS DEMAND
   INSTEAD.** Law 11b refuses an ARGUMENT that cannot be honoured. It does not refuse memory that cannot be had —
   and on the wasm tier it could not: under `-fno-exceptions` a throwing `new` that fails aborts inside the call
   (emsdk 6.0.9: `bad_alloc` → `abort()` → a JavaScript `RuntimeError`). Natively `bad_alloc` escapes the call — or
   ends the process where it meets a `noexcept` boundary, and in a composite the two can be one line apart: a
   composite that allocates a stage itself (an escape), then calls that stage's `noexcept` `prepare()`, which
   allocates its scratch (a `terminate`). So exhaustion is **outside the refusal contract on every row**: it ends the module's
   usefulness and is never answered with `false` — save where a third-party backend throws (the NAM backend,
   in felitronics-guitar-core, catches what its preparation throws and stays unprepared), a refusal this law
   neither asks of the other modules nor forbids there. This is the explicit exception to 11b, and it was chosen over nothrow
   storage plus a status in every allocating module on measured grounds: what the core holds is a constant of its
   CONFIGURATION plus a small fraction of the programme, a caller can read it before committing (below) and stay
   clear of exhaustion by arithmetic, and a later move to the nothrow form is ADDITIVE — a new status code, no
   struct moves — while making it now would cost every module the desktop products share. Two obligations replace
   the refusal:

   * **THE DEMAND.** An allocating call on the worker path can state, before it is made, a bound on how much of the
     heap its OWN requests will occupy at once — not what the object already holds — computed by the very functions
     its `prepare()` sizes itself with, so the bound cannot drift from the allocation; a C ABI over it forwards
     it. In 64 bits. What each number bounds is written where it is defined, because "at once" is
     not one formula: a call that keeps what it asks for is bounded by the sum of its requests, exact on a FRESH
     object (one already prepared keeps storage that still fits); a call that builds and frees per pass (a
     search that renders once per pass) by one pass. REQUESTED bytes, not a promise that a heap can serve them:
     allocator headers, the standard library's own alignment (MSVC's STL, in a release build, asks for
     `sizeof(void*) + 31` more on a block of 4096 bytes or more), fragmentation and the runtime's growth step are
     the caller's margin.

     **THE FUNCTION THAT SIZES IS THE FUNCTION THAT VALIDATES, and that is what makes "cannot drift" a
     construction rather than a promise.** Each allocating `prepare()` in this tree now has a static
     `storageFor(...)` beside it that answers FALSE on exactly the arguments the preparation refuses and
     otherwise fills in the element counts the buffers are built from — and `prepare()` calls it as its own
     gate. A budget is that function's `bytes()`. A stage cannot change what it refuses, or what it allocates,
     without changing both at once. Where a composite needs a stage's LATENCY to size something of its own, the
     stage publishes a static `latencyFor(...)` on the same terms and its own `prepare()` runs through it, so the
     number a budget reads and the number the prepared object reports are one expression.

     **A CALL THAT REFUSES MAY HAVE ASKED FOR PART OF ITS BOUND ON THE WAY — except where it can be decided for
     nothing, and then it must be.** A composite that can reach the whole verdict, every stage's included,
     without a single allocation publishes it as its own `admits()`, and a boundary over it then refuses an
     impossible geometry having touched no heap at all — on the tier where an allocation that cannot be served is
     not a refusal but the end of the module, which is what this law is about.

     **THE DEMAND IS A NUMBER, NOT A PERMISSION.** A boundary's demand query answers what a call would REQUEST and
     does not consult the handle's state: the cost of a call does not depend on when it is made, and a caller
     deciding whether to reset a stream and re-configure needs the number precisely when the other calls would be
     refused. The one exception is the call that has no handle to ask: the create's own query is a DRY RUN,
     returning every status the create would return before its first allocation, because the configuration it is
     handed has never been admitted anywhere — and because a budget of 0 must keep meaning one thing.
   * **A MODULE WHOSE CALL NEVER RETURNED ANSWERS EVERY STATUS CALL WITH "DISCARD ME".** The runtime does not stop a
     module that aborted; it answers the next call with objects wherever the abort left them (measured: a stream
     answering `OK` at the gain of a search that never returned, and a handle table left full for good). So a
     boundary that cannot outlive an abort marks every call in progress and, finding the mark on entry, answers
     "poisoned" for good and touches nothing — ahead of every other check, for every handle, and for a re-entrant call too, which it cannot
     tell apart. Entry points that read no instance state (build identity, the defaults writers) stay callable.

   Not promised: that a demand will be admitted, that anything survives exhaustion, or that a native host which
   catches `bad_alloc` holds a usable object. RT law 2 is unchanged — `process()` allocates nothing — so none of
   this reaches the audio path. Gated: each allocating stage's suite pins its published budget against the bytes its
   `prepare()` requests, with the ONE counter in `test_support/alloc_counter.h`, which installs EVERY form of
   `operator new`, the over-aligned one included. The worked case of this whole law — the mastering chain, its C ABI
   and the measurements behind every clause above (the escape and the `terminate` one line apart, the bytes a
   refusal used to ask for, the poisoned replay), with the gates that pin them over a matrix of rates, widths and
   topologies — moved with the chain to felitronics-mastering-core, `docs/LAW11D-MASTERING.md`.

   **11e. A RESTART NEVER LOSES AN ACCEPTED PUBLICATION — IT ADOPTS IT.** A swap-safe convolver publishes an
   operator from the message thread (`setIr()`/`setOperator()` return true) and the audio thread adopts it
   at the END of a crossfade, so between those two moments the live slot still names the PREVIOUS operator.
   All three of `convolution::{ConvolutionEngine, MatrixConvolver, MatrixConvolverNupc}::reset()` used to
   wipe the publication flag and keep that slot, which silently reverted to the operator the caller had
   already replaced, with nothing left to re-stage it — `setIr()` had already said true, and a consumer's
   retry flag is clear after a successful publish. On those bodies each class's independence group reads the
   OLD operator in every cell {Pending, Crossfading} x {mono, stereo} x both slot parities (a worst sample
   7.494e-01, 1.017e+00 and 1.927e+00 from the settled restart), and through `lineareq::LinearPhaseEq` a +12
   dB bell published over a flat curve read +0.00 dB after the restart — and exactly zero output (the
   fixture's -600 dB floor) when it was the first curve ever published. `reset()` now ENDS a swap in flight
   in favour of the new operator: `cur_` moves to the published slot and the state returns to Idle, so the
   consumer may publish again at once. The rule is scoped to the RESTART: `prepare()` is a re-initialisation
   and discards everything by contract, which is why every consumer re-publishes after it. What makes
   adoption right is 11a INDEPENDENCE, not click-freedom: a half-finished fade is a dependency on what came
   before, and two engines holding the same published operator — one mid-fade, one settled — would otherwise
   answer the next programme differently for up to the length of the fade. The price is at the seam and it
   is published: flushing the history is itself a cut (1.8203e-01 on DC 0.5 into a settled 700-tap
   operator), and adoption puts the new head tap where the old one was, which moves that step by at most the
   head-tap difference times the input — in EITHER direction (7.3203e-01 for a pair whose head tap flips
   +0.70 -> -0.40; 1.3933e-01 for a +1 dB broadband move and 2.2018e-01 for a -1 dB one; exactly the flush
   for a change that leaves the head tap alone). An operator merely STAGED and not yet published
   (`stageOperator()` without `publishStaged()`) is untouched: it has not been accepted, so the later
   publish still finds it. The Idle store is made ONLY when a swap was in flight, and it is a release store
   because the loader reads `cur_` after acquiring it; with no cached tail written off the audio thread
   either, `reset()` no longer races a single-producer loader at all — ThreadSanitizer, a loader publishing
   in a loop beside process / reset / clearAudioState on all three classes, reports 0 races against 54–62 on
   the bodies before (Apple clang needs `-fno-builtin`, or a `std::fill` write goes uninstrumented; gcc 14
   sees it as is). The documented "must not run concurrently" contract stays until a sanitizer row carries
   that; the store's condition is gated already — with it removed, a real-thread test loses 12–95 % of
   publications in every run. `clearAudioState()` is the history alone and leaves a fade running, which is
   exactly why it does NOT give independence mid-fade. The house precedent is `eq::EqBand::reset()`, which
   SNAPS a pending design onto the target rather than dropping it or playing out its ramp.


**These laws are CI-enforced for the funded tiers, not aspirational** — but not all of them, and the
difference is worth reading rather than assuming. Today: a
no-allocation test on the paths that install an allocation counter — one counter, in
`test_support/alloc_counter.h`, which replaces every form of `new` including the over-aligned one, probes
one allocation through each of the eight before `main()` and aborts naming any that went uncounted, and is
kept the only one by a lint (P52; before it, 50 of 61 private counters could not see an over-aligned
allocation at all) — a compile-only `-fno-exceptions` /
`-fno-rtti` probe over every public header on the Clang/GCC rows (MSVC spells the flags differently and is
not a gate for it), an **Emscripten** (`wasm-audio`) build + node run, a **`long double` scan** over
`modules/*/include` and `modules/*/src` paired with an **artifact** gate over the tier's objects (law 9),
and **`-ffp-contract=on`** on every non-MSVC row, whose removal reds the arm64 Linux row (law 10). That
catches an exception, a
dep, an alloc on an instrumented path, (via `--wrap=pthread_create`) a thread that is actually reached, and
a `long double` in core code whether it is emitted or not;
it does **not** catch a denormal at all, nor a thread in code that is never emitted — see the box in the tier list above rather
than trusting this sentence's older, broader wording. **Law 8 is now the only law here with no gate of any
kind**, for the reason it states: a denormal stall lives in the CPU at runtime, not in the build. An
**arm-none-eabi** job is documented but **not gated** until an embedded product funds it (avoid
embedded-grade CI scope creep). Third-party deps (Eigen, kissfft, pffft) must be verified to
build under `-fno-exceptions` (some use throwing asserts).

---

## 3. Module layout

One umbrella, **independent modules**, each separately includable so a consumer pulls only what it
needs (an EQ plugin must NOT drag in the neural runtime). **Hybrid build model:** the
*light* modules (`core`, `eq`, `dynamics`, simple analysis) are header-only `INTERFACE` targets; the
*heavy / platform-specific* ones (`convolution`, `neural`, the FFT backends, SIMD kernels) are
**compiled `STATIC`/`OBJECT` targets with narrow public headers**, so intrinsics, denormal handling,
size config, and per-platform flags stay out of public headers (no `#ifdef`-soup, ODR hazards, or
binary bloat). Every module is `felitronics::<module>`, JUCE-free, declares its target tiers, and has
its own JUCE-free self-tests.

| Module          | What                                                                 | Build / deps / notes |
|-----------------|---------------------------------------------------------------------|--------------|
| `core`          | Math, `Smoother`, `ScopedFlushToZero`, `kMaxChannels` (SSOT) + size config, fixed-size ring / lock-free SPSC FIFO, `Sample` alias, `DelayLine` / `DryAligner` / `StreamResampler`, the **FFT seam** | header-only, zero deps; the shared base. `StreamResampler` is a 64-tap polyphase windowed sinc since P34 (it was a Catmull-Rom cubic with no anti-aliasing through v0.26.0); what the swap fixed and what it costs in latency and CPU is [`STREAM-RESAMPLER-COST.md`](STREAM-RESAMPLER-COST.md) |
| `eq`            | matched biquads (Vicanek) + Cytomic SVF + `EqBand` + `EqEngine`      | header-only = today's `teq::` (becomes `eq`) |
| `dynamics`      | `EnvelopeFollower` (peak/RMS, attack/release) + `GainComputer` (threshold/ratio/knee/range, downward+upward) | header-only, zero deps; **first NEW module** |
| `analysis`      | spectrum tap (the existing `SpectrumTap`), correlation, LUFS/loudness | header-only; FFT **via the seam** |
| `convolution`   | partitioned (uniform / **zero-latency**) FFT convolution            | **compiled**; FFT **via the seam** |
| `oversampling`  | polyphase up/down-sampling                                           | compiled if SIMD; for true-peak + nonlinear |
| `limiter`       | true-peak limiter                                                    | compiled; uses `oversampling` |
| `neural`        | a thin **inference-object seam** (process-only) + the swap-safe model holder; model *loading* lives in the adapter | header-only, backend-free; a backend (the NAM one is in felitronics-guitar-core) is chosen at build time. NOT a runtime model-swap of one backend for another |

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
- **PRODUCT-SPECIFIC cross-module "glue."** Example: TabbyEQ's *dynamic EQ* = `eq` band + `dynamics`
  detector fed by the band's band-pass, applying the computed gain as an SVF gain-delta. That
  composition is TabbyEQ-specific and lives in TabbyEQ, **not** in `eq` or `dynamics`. A standalone
  compressor uses the same `dynamics` broadband. The core gives primitives; products compose them.

  **AMENDED — where the line actually is.** As written this read as "no composite belongs in the core",
  and the tree had already outgrown that in five places (`dynamiceq`, `deesser`, `multiband`, and two
  guitar composites since moved to felitronics-guitar-core) before `mastering` arrived. The rule those
  follow, stated properly — and "the core" in it is the family's shared layer, this repository or a
  satellite repository beside it:

  > A composite belongs in the core when **it is the unit under test** and its behaviour is shared.
  > It gets its OWN module, which may depend on many others; the primitive modules stay independent of
  > each other. What stays in the product is the **VOICING** — preset tables, target curves, defaults,
  > the taste.

  The discriminator is not "how many modules does it touch" but "can it be wrong on its own". A
  mastering chain can: its defects are latency arithmetic, block dependence, stale state across a
  bypass and a lost tail — none of which live in any stage, and all of which are only reachable by
  testing the composition. `felitronics::mastering` therefore ships in that shared layer — in
  felitronics-mastering-core, beside this repository, since this one keeps the base functions — while
  `mastering-config.json` (JAZZ/METAL, LIGHT/HEAVY) stays in the product, and that split is the rule, not an
  exception to it.

  **AMENDED AGAIN — the mastering brain moves into the shared layer too.** For mastering the product is now
  two shells, a web page and a desktop app, that must decide identically on the same file. The decisions —
  which devices to place and how, the voicing they read (targets, device defaults), the texts of the report
  and the project file itself — therefore live in `felitronics::session`, a module of
  felitronics-mastering-core, and the product keeps only what a shell is: the screen, the playback, the
  storage, the threads. `session` is offline code, not a real-time stage; the laws it keeps and the ones it
  does not are stated in that repository (`docs/SESSION.md`). The line above still holds for everything that
  has one product and no second shell: there the voicing stays in the product.

---

## 5. Build & repo strategy

- **Hybrid targets** (see §3): light modules header-only `INTERFACE`, heavy modules compiled
  `STATIC`/`OBJECT`; the umbrella `CMakeLists.txt` exposes `felitronics::core`, `felitronics::eq`, …
  Consumers link only what they use.
- **Consumption:** pinned tag via CMake `FetchContent` (exactly how the products already pin JUCE).
  OrbitCab already consumes `teq` this way (`SOURCE_SUBDIR teq`, pinned tag) — so `teq → felitronics::eq`
  is just repointing the URL/tag, not un-vendoring a folder copy.
- **Versioning is a contract, not just tags.** SemVer + `CHANGELOG.md`, and **no compatibility
  aliases**: a name whose documented meaning IS the defect is removed, and the consumers stop
  compiling until they move — that refusal is the notice. An alias that goes on answering the old way
  is a wrong answer with a warning attached, and the warning is the half nobody reads;
  `modelLoudness()` always meant "slot 0", which was the bug, and v0.29.0 deleted it rather than
  keep it politely lying. **DSP output is versioned behaviour** — a changed filter curve, limiter
  release shape, or convolver latency breaks presets/sessions as surely as an API rename → treat it
  as a breaking change. A **consumer-matrix CI** builds TabbyEQ + OrbitCab against each release; golden-audio
  vectors guard behaviour.
- **Portability CI for the funded tiers, from day one** (not "later"): the no-alloc /
  `-fno-exceptions` / `-fno-rtti` configs, scalar↔SIMD parity tests, and an **Emscripten**
  (`wasm-audio`) build of the light modules — all three now real, plus a native↔wasm bit-exactness NULL
  test on a generated fixture. An **arm-none-eabi** job is documented but **not gated**
  until an embedded product funds it (see §2).
- **License:** AGPL-3.0-or-later, SPDX header on every file. Each heavy module records its third-party
  deps + AGPL-compatibility in `THIRD_PARTY_NOTICES.md` (pffft/kissfft BSD = OK; watch ONNX/others).

---

## 6. Migration plan (keep every product green throughout)

1. **Establish conventions + the `core` API contract** (this doc): `prepare()` / `process()`, block
   views, the size constants, error returns, `noexcept` audio APIs, and the no-alloc test harness.
2. **FFT-seam spike FIRST** *(the keystone)*: define `core/Fft` (the plan object) and
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
- **Variable-order Notch** (`matched::notchCascade`, `designBand` Notch branch): the band-stop mirrors
  HP/LP — `slope` → order (`clamp(slope/6,1,16)`), realised as a Butterworth LP→BS cascade of
  `ceil(order/2)` matched sections (capped at 8) with the **zeros pinned at f0** (an order-fold infinite
  null that never drifts) and only the **poles staggered**; `Q` stays the −3 dB **width**, independent
  of order. Mapping: 6/12 dB/oct → 1 section (**== the legacy single matched notch, bit-for-bit** →
  old sessions don't drift), 24→2, 48→4, 96→8. The GUI curve picks it up for free (it iterates `d.sec[]`).
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

**3. Core modules TabbyEQ uses (now shipped) + future.**
- **`dynamics`** (TabbyEQ Phase 2; full design in `tabby-eq/docs/DYNAMICS.md`): core primitives
  `EnvelopeFollower` + `GainComputer` + `Compressor` (JUCE-free, software-flush per Law 8) now live in
  **`felitronics::dynamics`** (migrated with the `eq` extraction). The **dynamic-EQ composition** (detector = a Cytomic-SVF
  band-pass probe on the band's region; gain applied as a **matched-static × SVF gain-delta**, per M/S
  lane; GR-metering atomics for the host/UI) stays **product glue in TabbyEQ** (the §4 rule). De-esser
  = a preset on the same engine.
- **`analysis`** ← the pre/post `SpectrumTap`s + the `CorrelationMeter` (both now in `felitronics::analysis`).
- **`convolution` + `lineareq` — SHIPPED.** TabbyEQ now offers **Linear Phase** (linear-phase FIR, 5
  quality steps) and **Natural Phase** (mixed-phase FIR, φ=k·φ_min) via **`felitronics::lineareq`**
  (`LinearPhaseEq` / `NaturalPhaseEq` / `MixedPhaseFir`) over **`felitronics::convolution`**'s
  partitioned, click-free IR-swap engine (`ConvolutionEngine`) on the FFT seam — so dragging a band
  re-swaps the IR **artifact-free**. The dynamics×FIR interaction (dynamics is time-varying → bypassed
  in the FIR phase modes, `docs/DYNAMICS.md:112-117`) is now live behaviour, not a future constraint.
- possibly **`limiter`** on the output (future).

### Guitar products (OrbitCab, orbit-amp, the capture apps)
Their DSP beyond this core — the NAM backend, the `.orbitrig` pack player, the tube power-amp trunk —
lives in felitronics-guitar-core. The pre-migration inventory of OrbitCab's DSP that this section used to
carry left core with those modules.

### Future products
- Standalone **compressor** (Felitronics) — `dynamics` broadband + optional sidechain `eq`.
- **True-peak limiter** — `oversampling` + `limiter`.
- **Browser / WASM** demos of any of the above.
- **Hardware** pedal/box with a processor inside.

---

## 8. Open decisions

- **RESOLVED: build model = hybrid** — light modules header-only `INTERFACE`, heavy modules
  compiled `STATIC`/`OBJECT`. §3, §5.
- **RESOLVED: the FFT-seam + zero-latency convolution spike goes FIRST**, before dynamics. §6.
  **DONE (2026-06-29) — landed + green:** `felitronics::core::fft` (compile-time seam + a scalar radix-2
  reference backend), `felitronics::convolution::PartitionedConvolver` (zero-latency direct-head +
  uniform-partitioned overlap-save tail; verified == direct convolution under hostile variable block
  splits + boundary impulses; no-alloc-in-`process`), and an offline Kaiser `resampleIr` (≥55 dB). 120
  checks across the FFT / convolution / resampler suites. The seam is a **compile-time backend with an
  OPAQUE spectrum + a backend `spectralMultiplyAdd`** (so pffft/vDSP never pay an O(N) repack);
  real backends (pffft/kissfft/juce::dsp::FFT) plug in later. Production still needs a
  crossfade on live IR swap + Gardner non-uniform partitions for long IRs (both noted, not spike scope).
- **Sample type:** add a `using Sample = float` alias + non-`float*`-locked signatures **now**; full
  templating reaches `embedded-fpu` (float). **`bare-mcu` (fixed-point) is a SEPARATE codebase** — the
  alias does NOT flip a float SVF / NAM into fixed-point (different topology / scaling / quantization);
  don't promise it from one source.
- **RESOLVED: the FFT seam is compile-time** (template / C++20 concept), **not `virtual` in
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
- **`neural` seam = inference-object level, build-time backend choice**: the core
  defines a process-only inference interface; **model loading lives in the adapter** (the core accepts a
  pre-built inference object). "Per-platform" = the adapter links a different backend (NAM/Eigen on
  desktop+WASM; a lighter runtime / tiny model on hardware) — NOT a runtime swap of one model. Don't list
  RTNeural as a dependency until something actually uses it.

- **RESOLVED: `analysis::SpectrumTap` = the
  `teq::` shape.** `cab::SpectrumTap` and `teq::SpectrumTap` are diverged copies of one SPSC struct;
  the `teq::` one is the strict superset (JUCE-free `std::copy` + `reset()` + `tryPull()`). OrbitCab
  adopts it; the mono-sum vs channel-0 *feeding* policy stays per-product.
- **RESOLVED: SVF + Smoother dedup = REPLACEMENT, not code-merge.** OrbitCab's
  per-slot `juce::dsp::StateVariableTPTFilter` → `teq::Svf` (same Zavalishin TPT topology; HP/LP @
  12 dB/oct, Butterworth at Q≈0.707) and its `juce::SmoothedValue`s → a `core` `Smoother`. Both change
  the impl, not the math family → a **versioned-behaviour** swap (guard with a magnitude/phase golden test).

<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->
<!-- Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. -->

# #1 performance debt — the scalar-FFT / direct-head convolution bottleneck

**Status:** OPEN · highest-priority core optimisation · **Owner:** Darwin's Cat (Felitronics line)
**Impact:** latent for shipped OrbitCab · **blocker** for any FFT/convolution product (linear-phase EQ,
matrix / true-stereo convolution).

---

## TL;DR

`felitronics-core` ships **exactly one** real-FFT backend — `ScalarRadix2Real`
(`modules/core/include/felitronics/core/Fft.h:145`, `using DefaultRealFft = ScalarRadix2Real`) — a
correctness-first scalar Cooley–Tukey with `std::complex<double>` twiddles. **No SIMD backend
(pffft / vDSP / juce) exists yet.** Every FFT- and partitioned-convolution module rides it.

Compounding it, the partitioned convolvers also run a **scalar time-domain "direct head" FIR** — an
`O(P)` multiply-accumulate **per output sample** (`P` = partition / head length) — that is slow
*independently* of the FFT.

Net effect: long-convolution cost is fine at small host blocks but **grows pathologically at block
2048 and larger** — the exact failure OrbitCab's IR cab hit before it was moved off the core convolver
onto `juce::dsp::Convolution`.

There are therefore **two independent "тормозяки"**, and a real fix needs **both** addressed.

---

## How we found it

OrbitCab's IR cab was built on `convolution::ConvolutionEngine`. On an M5 Pro:

| Measurement | Scalar core path | `juce::dsp::Convolution` |
|---|---|---|
| Plugin cab, small block | ~6% RT | ~1% RT |
| **Plugin cab, block 2048** | **~15% RT** | **~1% RT (flat)** |
| Standalone convolver bench, block 512 | ≈2.24% RT | ≈0.42% RT (~5×) |
| Raw real-FFT bench (2048-pt) | baseline | ~7–20× faster (juce / pffft) |

Reverting the cab to `juce::dsp::Convolution` dropped it to **~1% flat across all block sizes**. That
sidestepped the problem *for OrbitCab* but left it fully intact in the core for every other convolution
consumer. This document is the core-level record + fix plan.

> The "**~15% @ 2048 → ~1% flat**" line is the whole story: the scalar path is acceptable at ≤128 and
> explodes as the block (and thus the partition size `P`) grows. Any fix must be judged **at 2048+**.

---

## The two bottlenecks

### 1. Scalar FFT backend — `modules/core/include/felitronics/core/Fft.h`
- `transform()` (`:109–136`) — scalar radix-2 Cooley–Tukey with per-butterfly `std::complex<double>`
  twiddles. Correctness-first, no SIMD, large constant factor.
- `spectralMultiplyAdd()` (`:94–104`) — scalar per-bin complex MAC (`spectrumFloats(N) = N`, `:51`).
- `using DefaultRealFft = ScalarRadix2Real` (`:145`) — the **only** concrete backend in-repo.

Used for every partitioned-tail forward/inverse FFT and every spectral multiply on the audio thread.

### 2. Scalar direct-head FIR — the partitioned convolvers
A zero-latency partitioned convolver runs a **time-domain "head"** in parallel with the FFT "tail". The
head is a length-`P` dot product executed for **every output sample**:

```
for (int i = 0; i < P; ++i) head += h0[i] * fr[-i];   // O(P) scalar MAC / output sample
```

- `modules/convolution/include/felitronics/convolution/PartitionedConvolver.h:111`
- `modules/convolution/include/felitronics/convolution/ConvolutionEngine.h:210` (doubles during a
  click-free IR crossfade — `:307–308`)
- `modules/convolution/include/felitronics/convolution/MatrixConvolver.h:259 / 334 / 372` (×bank: 2× for
  LR/MS-diag, 4× for Full; doubles during a swap — `:442`)

Because `P` is tied to the partition / host-block size, **this is what blows up at 2048+** — and it is
*not* fixed by a faster FFT alone.

---

## Scope — who rides it (most-severe first)

| Module | FFT / conv? | Backend | Hot-path cost | RT thread? | Bottleneck? | Sev |
|---|---|---|---|:--:|:--:|:--:|
| `convolution/ConvolutionEngine` + `PartitionedConvolver` | yes | partitioned on `ScalarRadix2Real` | `O(P)` scalar MAC / **output sample** (direct head) + per-chunk scalar FFT + `O(numParts)` spectral MAC; ×2 while crossfading | yes | **yes** | 🔴 high |
| `convolution/MatrixConvolver` | yes | own 2×2 partitioned on `ScalarRadix2Real` | `O(P)` / sample **× bank** (2–4×) + tail FFT / spectral MAC | yes | **yes** | 🔴 high |
| `lineareq/LinearPhaseEq` + `NaturalPhaseEq` | yes | `MatrixConvolver<DefaultRealFft>` | `O(P)` direct head + amortised `O(firLen/P)` tail; FIR up to **131072** taps | yes | **yes** (partial¹) | 🔴 high |
| `oversampling/PolyphaseOversampler` | FIR, no FFT | scalar polyphase Kaiser-sinc | `O(taps)` / sample, **short** filters (tpp=32, N≤256) | yes | no (inherent) | 🟢 low |
| `lineareq/MixedPhaseFir` | FFT, design-time | `ScalarRadix2Real` in `build()` | 3× `O(D log D)` FFT on the **message** thread; 0 on audio | no | no | ⚪ none |
| `neural/Inference` + `NeuralStage` | no | plumbing (concept + swap holder) | `O(1)` / block dispatch | no | no | ⚪ none |
| cheap-sweep² | no | biquad/SVF · envelope · waveshaper · short FIR | `O(1)` / sample / band | yes | no | ⚪ none |

¹ *partial* = the FFTs are correctly **amortised** (one chunk per `P` samples via overlap-save), not a
naive whole-IR FFT per block — but every kernel (butterflies, `spectralMultiplyAdd` over ~255
partitions, the `O(P)` head) is scalar with zero SIMD, i.e. the same cost *class* the cab shed.

² cheap-sweep = `eq · dynamics · saturation · stereo · dither · limiter · multiband · deesser ·
dynamiceq · analysis`. All confirmed clean: pure per-sample scalar DSP; the only `fft`/`convol` grep
hits are **comments / compat shims** — `analysis/SpectrumTap.h:63` is copy-only on the audio thread (its
FFT runs on the GUI thread via `tryPull()`), and `eq/EqBand.h:274`'s `MatrixBasis` enum only *describes*
a downstream topology (the `eq` module includes no convolution headers; it is `Svf` / `MatchedBiquad`).

---

## Root cause (one paragraph)

There is exactly one concrete real-FFT backend, `ScalarRadix2Real`, aliased as `DefaultRealFft`
(`Fft.h:145`). It is scalar Cooley–Tukey with `std::complex<double>` twiddles (`transform()`
`Fft.h:109–136`) plus a scalar per-bin complex MAC (`spectralMultiplyAdd` `Fft.h:94–104`); no
SIMD backend (pffft / vDSP / juce) exists in-repo. Every partitioned convolution rides it —
`PartitionedConvolver` and `ConvolutionEngine` (`convolution/*`), `MatrixConvolver`
(`MatrixConvolver.h:50`) and, through it, `LinearPhaseEq` / `NaturalPhaseEq` (`LinearPhaseEq.h:189`,
`NaturalPhaseEq.h:198`). On top of the scalar FFT there is a **second, FFT-independent bottleneck** —
the scalar direct-head FIR `for (i<P) h += h0[i]*fr[-i]` run **every output sample**
(`PartitionedConvolver.h:111`, `ConvolutionEngine.h:210`, `MatrixConvolver.h:259/334`), doubling during
a crossfade (`ConvolutionEngine.h:307–308`, `MatrixConvolver.h:442`). Precision note: the per-output-sample
head cost is `O(P)` where **`P` is the partition SIZE (head length), not `numParts`**; `numParts` drives
only the per-chunk spectral MAC. Host block size does not change the per-sample cost directly (the
internal `2P` frame decouples it) — but `P` scales *with* the block/partition, which is why the cost
climbs with larger host buffers.

---

## Impact

- **OrbitCab (shipped): latent — not burning CPU today.** Its only long-convolution path (the IR cab)
  already delegates to `juce::dsp::Convolution` (SIMD, off `ScalarRadix2Real`); its `teq` tone EQ is
  biquads; `analysis::SpectrumTap` is copy-only on the audio thread. Everything OrbitCab still runs on
  the scalar path is short filters / biquads.
- **Future & other products: blocker.** The moment `LinearPhaseEq` / `NaturalPhaseEq` (FIR up to 131072
  taps — `LinearPhaseEq.h:52`) or `MatrixConvolver` / true-stereo convolution ship in a real plugin, they
  hit the audio thread on the scalar backend — the same cost class the cab already had to abandon.

---

## Fix plan (in leverage order)

### (a) One fast `RealFftBackend` behind the existing concept seam — **maximum leverage**
Every convolver is **already templated** on the `RealFftBackend` concept:
`MatrixConvolver.h:50` (`template <class Fft = core::fft::DefaultRealFft>`), instantiated at
`LinearPhaseEq.h:189` / `NaturalPhaseEq.h:198`; `PartitionedConvolver` / `ConvolutionEngine` likewise.
Provide a SIMD real-FFT backend satisfying the same concept and switch `DefaultRealFft` (`Fft.h:145`):

- **pffft** — BSD, desktop + WASM (NEON needs `PFFFT_ENABLE_NEON=1`); native `zconvolve_accumulate` maps
  onto our `spectralMultiplyAdd` seam.
- **vDSP / Accelerate** — Apple desktop.
- **`juce::dsp::FFT`** — in the JUCE adapter layer.
- **CMSIS** — far-future embedded.

This accelerates **every** amortised forward/inverse FFT (`PartitionedConvolver.h:124/140`,
`MatrixConvolver.h:389/271`) and every `spectralMultiplyAdd` (`Fft.h:94–104`; call sites
`PartitionedConvolver.h:132`, `MatrixConvolver.h:276/287`) across **all** modules **at zero API change**.
Bonus: design-time `MixedPhaseFir` (`MixedPhaseFir.h:63/105/120`) speeds up too.

> **Prototype already exists.** Validated `JuceRealFft` + `PffftRealFft` adapters (both satisfy the
> `RealFftBackend` concept) plus a 3-way correctness-null + benchmark harness are stashed in the OrbitCab
> repo as dev tooling (tracked under felitronics-core issue #21). The bench is how the numbers above were
> produced — reuse it as the acceptance gate.

### (b) SIMD direct-head FIR — the per-sample floor
Even after (a), the scalar `O(P)` head (`PartitionedConvolver.h:111`, `ConvolutionEngine.h:210`,
`MatrixConvolver.h:259/334/372`) stays scalar and **dominates at large `P` / short IRs** (a guitar-cab IR
has few tail partitions, so the head *is* the cost). Vectorise the dot product (SIMD MAC / partition-0
unroll). This is the part a faster FFT does **not** fix, and the reason acceptance is measured at 2048+.

### (c) Polyphase oversampler SIMD — low priority, inherent
`PolyphaseOversampler` (`:84/:105`) is a correct short-filter scalar FIR (tpp=32, N≤256). Not a
bottleneck — FFT/partitioning would not help (filters are far too short). Optional ~2–4× SIMD MAC only
(documented future backend, `:32`).

---

## Acceptance criteria

A fix is accepted **only when BOTH hold** — output-match alone is **not** sufficient:

1. **Correctness (necessary, not sufficient).** The existing reference-NULL property tests stay green
   (ULP-bounded vs a trusted convolver); bit-for-bit behaviour unchanged. *A null match ("совпадение")
   proves nothing about speed.*

2. **Comparable performance at scale — the real bar.** Measured RT cost at host blocks **{2048, 4096,
   8192}** (with representative IR / FIR lengths) must be **within a small constant factor of a reference
   SIMD convolver** (`juce::dsp::Convolution` / pffft) **and must NOT grow pathologically with block
   size**. "Faster than the scalar baseline" is *not* enough. Target: reproduce, for the core convolvers,
   the flat-low profile the cab got from juce (**~1% at 2048**, vs the scalar **~15%**). **Benchmark at
   2048 and larger is the gate** — that is precisely where the scalar path breaks and where the win must
   show. Small-block parity is necessary but does not close the item.

Both `(a)` and `(b)` are effectively required to pass criterion 2: at block 2048 the `O(P)` scalar head
alone is ~2048 MACs/sample and will keep the numbers red even with a perfect FFT.

---

## Reference — the hot lines

- **Backend:** `core/Fft.h:145` (alias), `:109–136` (scalar transform), `:94–104` (scalar spectral MAC), `:51` (`spectrumFloats`).
- **Direct head (scalar FIR):** `PartitionedConvolver.h:111`, `ConvolutionEngine.h:210` (+ `:307–308` fade), `MatrixConvolver.h:259 / 334 / 372` (+ `:442` fade).
- **Tail FFT / spectral MAC:** `PartitionedConvolver.h:124 / 132 / 140`, `ConvolutionEngine.h:219`, `MatrixConvolver.h:271 / 276 / 287 / 389`.
- **Template seam (where a SIMD backend plugs in):** `MatrixConvolver.h:50`, `LinearPhaseEq.h:189`, `NaturalPhaseEq.h:198`.
- **FIR sizes:** `LinearPhaseEq.h:52` (up to 131072), `NaturalPhaseEq.h:55` (up to 16384).
- **Clean (for the record):** `MixedPhaseFir.h:34` (RT-unsafe / design-time), `NeuralStage.h:67` (`O(1)` dispatch), `PolyphaseOversampler.h:32` (short scalar FIR), `analysis/SpectrumTap.h:63` (copy-only tap).

<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->
<!-- Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. -->

# WASM spike — does the core actually live under Emscripten?

**Verdict: yes, and bit-exactly.** The `wasm-audio` tier that [`DSP-ARCHITECTURE.md`](DSP-ARCHITECTURE.md) §2
declares is real, not aspirational: `analysis` + `oversampling` compile under Emscripten with
`-fno-exceptions -fno-rtti` and no pthreads, run in a browser with no COOP/COEP headers, and produce
**bit-identical doubles** to the native reference — provided one compiler flag is set, which is the whole
substance of this document.

The spike exists because *"it compiled and linked"* is not evidence. The precedent that motivated it: in a
sibling project libsoxr compiled and linked cheaply under Emscripten and then **crashed at runtime**,
discovered only after the build.

**Toolchain, recorded explicitly.** emsdk **6.0.9**, emscripten commit
`4e4223852a0835923411059a3929907d7df1232e`, bundled node 24.19.0. Native reference: Apple clang 21.0.0
(arm64 macOS) and gcc 14.2.0 (x86-64 Debian 13). Reproduce with `tools/wasm/build.sh`.

---

## What was built

| | |
|---|---|
| `tools/fcore_probe.h` | `fcore::Probe` — the measurement body, shared **verbatim** by the native CLI and the wasm shim. Chunks internally at 8192 frames, so a whole-file buffer costs the same bounded scratch as a stream. |
| `tools/fcore_measure.cpp` | The native reference, refactored onto `Probe`. New `blocks` mode; `--precise`. |
| `tools/wasm/fc_probe.cpp` | The C ABI: planar-contiguous, one pointer, channel `c` at `planar + c*frames`. |
| `tools/wasm/build.sh` | Three artifacts — `web` (the specified one), `node` (harness), `debug` (SAFE_HEAP + assertions). |
| `tools/wasm/parity.mjs` | The wasm side of the parity check; emits exactly what the CLI emits, so `diff` **is** the test. |
| `tools/wasm/probe.html` | The page. Deliberately ugly. |

---

## What is compared, and why it is not the obvious thing

The naive acceptance test — "the LUFS from wasm matches the native LUFS to 1e-9" — **cannot be sound**, and
it took building the thing to see why.

Integrated LUFS is **discontinuous in its own inputs**. BS.1770 gating is a pair of strict comparisons
(`LoudnessMeter.h`, the absolute gate at −70 LUFS then the −10 LU relative gate). A 400 ms block sitting
within ~1e-12 of a threshold flips its inclusion between two builds and moves the answer by ~**0.01 dB** —
seven orders above any tolerance worth writing down. No scalar tolerance on a gated quantity is a robust
equivalence measure; it merely happens to pass.

So the comparison surface is the **pre-gate 400 ms block-energy vector**, newly exposed as
`analysis::LoudnessMeter::gatingBlockEnergies()`, plus the true-peak **linear** maximum. Both are continuous
in the input samples and can carry a bit-exactness claim. The scalar LUFS and dBTP remain, for humans.

This turned out to be far more than a principled nicety — see the negative control below. **The vector is
roughly three thousand times more sensitive than the scalar it replaced.**

---

## Acceptance

### 1. Bit-exactness — native vs wasm

| file | rate / ch | block energies | result |
|---|---|---|---|
| `click120` | 44.1 kHz mono | 137 | **bit-identical** |
| `gr-clip` | 48 kHz stereo | 247 | **bit-identical** |
| `cold-gaze` | 48 kHz stereo | 3213 | **bit-identical** |

3597 doubles plus three true-peak maxima, every bit equal, native arm64 vs wasm32.

### 2. Size

| file | raw | gzip | brotli |
|---|---|---|---|
| `fcprobe.web.wasm` | 30 673 | 16 311 | **14 660** |
| `fcprobe.web.js` | 8 719 | 3 092 | 2 788 |

The two release `.wasm` (web glue and node glue) are **byte-identical** — `-sENVIRONMENT` shapes JS, not
code — which is what lets a parity result proven under node transfer to the browser.

### 3. No COOP/COEP, and no threads

Verified in a real browser against `python3 -m http.server`, which sends no isolation headers:
`crossOriginIsolated === false` and `SharedArrayBuffer === undefined` — threads are not merely unused, they
are **unavailable**, and the module works anyway. Its printed output was bit-identical to the native CLI.

Stronger and cheaper than the page test, and now part of `build.sh`: the artifact itself carries **zero**
occurrences of `SharedArrayBuffer`, `pthread`, `new Worker` or `Atomics`, and its memory is not shared.

### 4. `-fno-exceptions` clean

Clean **for the spike**, which does not include `core/Fft.h`. It will not be clean for a CI gate over the
whole `analysis` module — see the findings.

### 5. Rate sweep — attached as evidence, deliberately not a gate

Bit-exactness was measured at 44.1 / 48 / 88.2 / 96 / 176.4 / 192 kHz. Two of those rates do **not** match.
That is a property of the libms involved rather than of this build, so it feeds the determinism question
below instead of failing the spike.

---

## The load-bearing flag: `-ffp-contract=off`

Baseline WASM has **no scalar FMA instruction**, so Emscripten cannot contract `a*b+c`. A native build that
does contract computes different doubles. This is not a tolerance question; it is a different computation.

**Negative control.** The same native tool, same source, built with contraction (the compiler default):

| file | block lines differing from wasm |
|---|---|
| `click120` | **137 of 137** |
| `gr-clip` | **248 of 248** |
| `cold-gaze` | **3188 of 3213** |

Nearly every block. The scalar LUFS hides this entirely — it averages the difference down to ~7e-15 dB — which
is precisely why the criterion moved onto the vector.

**The flag is architecture-sensitive, and that is a trap.** On arm64 `fmadd` is in the baseline ISA, so clang
contracts by default and the parity check fails immediately. On baseline x86-64 there is no FMA instruction,
so gcc *cannot* contract even at its `-ffp-contract=fast` default — parity passes **by accident**. Add
`-march=native` on the same x86 machine and gcc emits 20 FMA instructions and every line differs again.

> A CI job on default-flags x86-64 would therefore **never catch** the flag going missing, while the
> developer's own arm64 Mac fails instantly. The flag is set explicitly in `tools/CMakeLists.txt`, and
> `felitronics_probe_tests` asserts it actually took effect: dropping it turns a silent parity failure into a
> red test.

---

## Cross-tier determinism: the outlier is Apple's libm, not wasm

Bit-exactness holds at the rates the corpus uses — and **not** everywhere:

| rate | native arm64 (clang) vs wasm | native x86-64 (gcc) vs wasm |
|---|---|---|
| 44.1 kHz | identical | identical |
| 48 kHz | identical | identical |
| **88.2 kHz** | **395 of 397 blocks differ** | **identical** |
| 96 kHz | identical | — |
| 176.4 kHz | identical | — |
| **192 kHz** | **396 of 397 blocks differ** | — |

The cause is `analysis::KWeightingFilter::computeCoeffs`, which designs the BS.1770 coefficients with
`std::tan` and `std::pow` in `double`. Emscripten's musl `tan` and Apple's `tan` disagree by 1 ulp at several
rates; at 88.2 and 192 kHz that disagreement survives into the coefficients instead of being absorbed by
rounding, and a changed coefficient changes the filter for *every* sample — hence nearly every block, not one.

The three-way comparison is the useful part: **x86-64 Linux/gcc agrees with wasm at 88.2 kHz while arm64
macOS does not.** glibc and musl agree; Apple's libm is the odd one out. So the divergence lives on the
development machine, and Linux CI would never see it.

Two consequences worth stating plainly:

- 44.1, 48 and 96 kHz pass **by rounding absorption, not by construction**. Any libm update on either side
  could end that. 96 kHz is heavily used, so this is not an exotic corner.
- The core has no stated cross-platform determinism contract. It needs one before any facade can promise
  "bit-identical output" across tiers. The cheap route is to remove libm from the compared path entirely —
  the coefficients are designed from two transcendentals at a handful of standard rates.

---

## Emscripten gotchas found

1. **Exports are opt-in.** In 6.x, `Module._malloc`, `Module._free` and the `HEAP*` views are `undefined`
   unless named in `-sEXPORTED_FUNCTIONS` / `-sEXPORTED_RUNTIME_METHODS`. The page dies on first use.
2. **Memory grows *inside* the C call.** With `ALLOW_MEMORY_GROWTH`, `memory.grow` replaces the underlying
   `ArrayBuffer` and **detaches every existing view**: writes are silently dropped, reads come back
   undefined. It is not enough to re-read the heap after `_malloc` — the meter's own ~1.27 MB `prepare()`
   can cross a growth boundary during the call. Never hold a `HEAP*` view or a `subarray` across any exported
   call; return results by value.
3. **Never `-mrelaxed-simd` on anything checked by a null test.** It lets LLVM emit `f64x2.relaxed_madd`,
   which is implementation-defined — fused on hosts with FMA, unfused elsewhere. That breaks determinism
   between *machines*, not merely between tiers. Plain `-msimd128` lowers to `f64x2.mul/add` and is safe.
4. **`long double` is 16 bytes on wasm32** (IEEE quad, software) against 8 on arm64 macOS and 80-bit x87 on
   x86-64 Linux — three tiers, three answers.
5. **`>>> 2`, never `>> 2`,** when converting a pointer to a `HEAPF32` index: above 2 GB a signed shift goes
   negative.
6. **A checked build is not optional.** `-O1 -g -sASSERTIONS=2 -sSAFE_HEAP=1 -sSTACK_OVERFLOW_CHECK=2` runs
   clean and bit-identical to release. That, not the release build working, is what rules out the libsoxr
   failure mode.

Non-problems, checked and dismissed: `LoudnessMeter::prepare(fs, nc, 4 h)` allocates **1 269 696 bytes**, not
a heap threat; 64-byte `align_val_t` allocation works under emmalloc; the 64 KB default stack is ample.

---

## Fixed here, because the spike's reference tool cannot be knowingly wrong

The post-build review consilium found that **the reference tool was under-reporting true peak on any file
that ends on a transient**, and had been doing so since it was written. The polyphase oversampler is causal
with a group delay of 63.5 oversampled samples, so the reconstruction of the final ~16 baseband samples never
left the filter while input was still arriving. Measured, before the fix:

| signal | reported | truth |
|---|---|---|
| unit impulse in the last sample | 0.000071 | 0.881 |
| last ten samples at 0.95, abrupt end | 0.0152 (**−36 dBTP**) | 1.0625 (**+0.53 dBTP** — over full scale) |

A true-peak tool that reports −36 dBTP for a clipping ending is worse than no tool, and both tiers would have
agreed on it bit-exactly. Two fixes, both in `fcore::Probe`:

- **`finish()`** drains the FIR with `kOsTapsPerPhase` silent samples (the oversampler only — feeding the
  loudness meter would append program that was never submitted). Idempotent; ends the stream.
- **the sample peak is a hard floor** under the true peak. The reconstruction passes through the samples by
  construction, so `truePeak >= samplePeak` always. On its own the floor recovers most of the error above
  (0.0152 → 0.95); `finish()` supplies the inter-sample part the floor cannot know about.

Also hardened while there, all of it facing input a page computes: the sample rate is bounded to
1 kHz … 768 kHz (positive-and-finite is not enough — `lround(0.01*fs)` on 1e300 is out of range, and a page
can pass `Number.MIN_VALUE` as easily as 48000), `maxDurationSec` is validated, a non-positive channel count
is refused rather than silently recording invented silent blocks, and the C ABI checks that the span it was
handed lies inside the wasm heap instead of trapping on a pointer near the top of linear memory.

---

## Findings for the core (not fixed here — this is a spike)

- **`analysis::KWeightingFilter` violates portability law 8.** It is the only feedback kernel in the core with
  no software denormal flush. After 2.96 s of digital silence its state reaches the subnormal range and then
  **sticks permanently** — measured out to 600 s, frozen at `-1.240105e-321`. On x86 without hardware FTZ,
  which is what a browser gives, **processing silence then costs 51× more than processing music**
  (22.62 ms vs 0.44 ms per 2 s, i9-13900H). Invisible on Apple Silicon, which handles subnormals at full
  speed — which is why it survived this long.
- **One `throw` blocks the `-fno-exceptions` CI gate.** `core/Fft.h`'s `SeamAllocator::allocate` is the only
  one in `core` + `analysis` + `oversampling`. The spike does not include that header; a gate over the whole
  `analysis` module will, through the three spectrum panes.
- **The core holds two different true-peak filters.** This tool's path (`PolyphaseOversampler` at 4×, 32
  taps/phase, cutoff 0.90× Nyquist, Kaiser β9) and `analysis::TruePeakMeter` (48 taps, full base Nyquist, β8,
  factor by rate). They agree to 0.0009–0.0028 dB on music but diverge by **~1.1 dB** on single-sample
  impulses, where the 0.90 guard discards the top of the band and this path reads **low**. Which is
  spec-correct is unresolved: ffmpeg's `ebur128` prints one decimal and cannot arbitrate.
- **A single non-finite sample silently freezes the loudness.** The two paths fail differently, and only one
  of them recovers. The true peak is blind for 32 samples while the NaN sits in the polyphase ring
  (`std::max(x, NaN)` returns `x`) and then works again — feed it louder audio afterwards and it duly rises.
  K-weighting is an IIR, so its state is NaN *forever*: every later block energy is NaN, `NaN > absT` is
  false, and the absolute gate silently drops all of them. The meter then goes on reporting a healthy number
  computed over the fraction of the program that predates the NaN, and `droppedBlocks()` says nothing.
  Pinned by `felitronics_probe_tests`.

---

## Reproducing

```sh
source ~/emsdk/emsdk_env.sh
./tools/wasm/build.sh                                        # three artifacts + size + thread audit
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DFELITRONICS_BUILD_TESTS=ON && cmake --build build -j

ffmpeg -i song.flac -f f32le -acodec pcm_f32le song.f32
./build/tools/fcore_measure blocks 48000 2 song.f32              > native.txt
node tools/wasm/parity.mjs tools/wasm/build/fcprobe.node.js 48000 2 song.f32 > wasm.txt
diff native.txt wasm.txt          # empty output IS the acceptance criterion
```

`-DCMAKE_BUILD_TYPE=Release` is not optional: the analyzer perf tests fail an unoptimised build.

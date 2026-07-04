<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

# Convolver CPU gap vs juce::dsp::Convolution (the next perf item)

**Status: OPEN.** The original #1 debt (scalar-FFT + O(P) direct-head cost that EXPLODED with the host block)
is RESOLVED by #23–#26 (SIMD pffft backend + fixed `P=128` → block-independent). This document records a
*separate, newly-measured* finding: our partitioned convolver, while now block-independent and JUCE-free, is
**2–13× more CPU than `juce::dsp::Convolution` at host blocks ≥ 512** (both zero-latency). Written down so we
don't forget the problem or the plan.

## The measurement (tools/fftbench.cpp, `-DFELITRONICS_BENCH_JUCE=ON`)

131072-tap stereo IR, 48 kHz, steady-state, FTZ on, %RT = CPU / audio time. Both report `getLatency()==0`.

| host block | `juce::dsp::Convolution` | our `MatrixConvolver<pffft>` |
|---:|---:|---:|
| 256  | 1.80% | 1.99% |
| 512  | 0.87% | 2.00% |
| 1024 | 0.46% | 2.04% |
| 2048 | 0.26% | 1.99% |
| 4096 | 0.17% | 1.99% |
| 8192 | 0.15% | 2.00% |

JUCE gets ~12× cheaper from block 256→8192; ours is flat ~2%. (i9-13900H/SSE: same shape, ~2× the M5 %RT in
absolute terms.) At the common 512–2048 blocks JUCE is **2–8× cheaper**. We are comparable only at block 256.

This was verified NOT to be a benchmark artifact:
- **Correctness probe** — a unit impulse through JUCE reproduces the full 131072-tap IR: `maxErr 6.5e-9`
  (rel `3.3e-7`), `energy out/IR = 1.0000`. JUCE really convolves all 131072 taps.
- **codex** read JUCE's source: default `Latency{0}`, `fftSize = 2·maximumBlockSize`, `numSegments = L/P + 1`
  — genuinely zero-latency, uniform-partitioned, not truncating.
- **deepseek** derived the per-sample cost; the ×12 matches (see below).

## Why (the physics — independently derived + crew-confirmed)

Per-output-sample cost of a **uniform-partitioned** overlap-save FFT convolver with partition size `P` and IR
length `L`:

```
cost/sample  ≈  c1·log2(2P)        (forward+inverse FFT, grows slowly)
             +  c2·(L/P)           (spectral MAC over L/P partitions, shrinks with P)
```

- **JUCE** sets `P = maximumBlockSize` (from `prepare`). Bigger block ⇒ far fewer partitions
  (`L/P` = 512 at P=256 → 16 at P=8192) ⇒ the MAC term collapses ⇒ ~12× cheaper. Zero-latency because
  overlap-save needs only the current + past blocks, never the future.
- **Ours** (`MatrixConvolver`) is a *fixed-`P=128` head+tail*: a time-domain direct head of `P` taps/sample
  (scalar) + an FFT tail of `L/P = 1023` partitions. Chosen in #26 to make cost **block-INDEPENDENT** (the old
  `P≥maxBlock` grew the O(P) scalar head until it exploded — 39%RT@8192). But 1023 partitions + the scalar head
  ⇒ a flat ~2% that does not benefit from large blocks. deepseek: ~98% of our 2% is the 1023-partition MAC.

The key insight: the old code's bug was **the scalar O(P) head**, not `P=block` per se. JUCE uses `P=maxBlock`
with **no scalar head** (pure uniform FFT), so it never explodes AND wins at large blocks. #26 fixed the
explosion by shrinking `P`; it did not adopt JUCE's block-adaptive partitioning.

## The solution — a uniform-FFT convolver (P = maxBlock, no scalar head)

Add a **uniform-partitioned FFT convolution** path: partition size = `maxBlock` (fixed at `prepare`), overlap-
save, **no time-domain head**, zero-latency for block-aligned processing (host block ≤ P). This is what JUCE
does and what our OLD code almost did — minus the scalar head that made the old code explode.

Open design questions (for the architecture consilium):
1. **A new convolver class** (`UniformConvolver`) vs **a mode of MatrixConvolver**. The head+tail machinery
   (MSDiag/Full topologies, warm FDL, click-free 2-slot swap) must be preserved either way — the lineareq
   consumers depend on it.
2. **Block-size handling.** Uniform `P=maxBlock` is zero-latency only for host blocks ≤ P. Our current head+
   tail handles ANY block (per-sample). Decide: require block ≤ maxBlock (the plugin norm), or keep both paths.
3. **Latency semantics.** Uniform `P=maxBlock` is "block-granular zero latency" (getLatency()=0 for block-
   aligned hosts). Our head+tail is "sample zero latency". Both report 0; document the distinction.
4. **The C1 design/audio split (from #24) still holds** — design FFTs stay scalar packed-Hermitian; only the
   audio uniform-FFT rides the SIMD backend.
5. **Reuse the FDL + spectral MAC** already SIMD-aligned (#23) + the pffft backend (#25).

**Acceptance:** match JUCE's %RT within ~2× across blocks 512–8192, `getLatency()==0`, and NULL-correct output
vs the current `MatrixConvolver` across all topologies. Bench: `tools/fftbench.cpp` (add a uniform column).

## Reproduce

```
cmake -S . -B build-juce -DCMAKE_BUILD_TYPE=Release -DFELITRONICS_WITH_PFFFT=ON -DFELITRONICS_BENCH_JUCE=ON
cmake --build build-juce --target fcore_fftbench -j
./build-juce/tools/fcore_fftbench          # prints the correctness probe, the OLD/NEW/pffft sweep, and the JUCE head-to-head
```

<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

# Convolver CPU gap vs juce::dsp::Convolution (the next perf item)

**Status: Phase-1 LANDED** (the mono/LRDiag `NonUniformConvolver` primitive — see the Phase-1 result below;
matrix routing + click-free swap are later phases). The original #1 debt (scalar-FFT + O(P) direct-head cost
that EXPLODED with the host block) is RESOLVED by #23–#26 (SIMD pffft backend + fixed `P=128` → block-
independent). This document records a
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

## The solution — a Non-Uniform Partitioned Convolver (Gardner 1995)

The consilium (me + codex + deepseek + Fable-5) chose a **Non-Uniform Partitioned Convolver** over a uniform
`P=maxBlock` path. A uniform-`P=maxBlock` would only match JUCE by *becoming* JUCE (its cost tied to the block);
NUPC instead keeps **true sample-zero-latency** AND a **block-INDEPENDENT** cost.

**The schedule.** The IR is tiled by a time-domain **head** `[0, P0)` (direct FIR, immediate) + a list of
overlap-save **stages** of *growing* block size: stage `s` has block `B_s`, count `C_s`, covers
`[offset_s, offset_s + C_s·B_s)`, FFT size `2·B_s`, its own frequency-domain delay line of depth `C_s`.
**Zero-latency invariant:** `offset_s == B_s == cumulative coverage of head + earlier stages`, i.e. the
recurrence `B_{s+1} = (C_s + 1)·B_s`, `B_0 = P0`. Each stage's inherent `B_s`-sample overlap-save latency is
exactly absorbed by the earlier coverage (Gardner's condition) → total added latency 0. The schedule is a
**tunable `{blockSize, count}` list** (Fable: pure octave-doubling is ~26% suboptimal for a short near-field —
a few DISTINCT block sizes each REPEATED, e.g. head 128 + 128×3 + 512×7 + 4096×31, is cheaper); the default is
capped-doubling to `B_max` then one uniform `B_max` stage. **Cost is `10·log2(2B) + 6·C` flops/sample per
stage** (the prior "~20 MAC" estimate ignored FFT cost, which dominates — verified by all three advisors).

**Why we do NOT beat JUCE everywhere (a theorem, not a shortfall).** True sample-zero-latency FORCES small
partitions in the near field `[0, B_max)` — an FFT-overhead tax JUCE never pays because it uses one large
partition for the whole IR. So JUCE is cheaper on the *mean* at large blocks. But NUPC wins outright exactly
where it matters: the small host blocks (64–256) that low-latency/live rigs run.

**C1 + reuse.** The design/audio FFT split (#24) still holds — design FFTs stay scalar packed-Hermitian; only
the audio path rides the SIMD backend. Reuses the SeamAllocator-aligned FDL + spectral MAC (#23) + pffft (#25).

## Phase-1 result (mono/LRDiag primitive — `NonUniformConvolver<Fft>`, LANDED)

Phase 1 = the bare mono primitive (LRDiag = two instances), NULL-verified, benched. Later phases add the
matrix-routing facade (MSDiag/LRDiag/Full) + the click-free warm swap on top.

**Measured (M5 Pro, 131072-tap stereo LRDiag, 48 kHz, pffft, `getLatency()==0`):**

| host block | `juce::dsp::Convolution` | `NonUniformConvolver<pffft>` mean (worst buffer) | ratio |
|---:|---:|---:|---:|
| 64   | 7.26% | **0.80%** (max 7.82%) | **9.1× us** |
| 128  | 2.77% | **0.81%** (max 4.58%) | **3.4× us** |
| 256  | 1.79% | **0.78%** (max 2.37%) | **2.3× us** |
| 512  | 0.90% | 0.79% (max 1.59%) | ~tie |
| 1024 | 0.44% | 0.78% (max 1.16%) | 1.8× JUCE |
| 8192 | 0.15% | 0.79% (max 0.82%) | 5.3× JUCE |

- **NUPC mean is FLAT ~0.8% at every block** (block-INDEPENDENT); JUCE swings ~48× (0.15%→7.26%).
- **NUPC is 3–9× cheaper at the 64–128 blocks a live guitar amp / monitor runs** — its target regime.
- **2.4× cheaper than our own v0.4.0 `MatrixConvolver<pffft>`** (flat ~2.1% → ~0.8%).
- **Both are zero-latency.** (Correcting an earlier framing: JUCE's `Convolution` reports `getLatency()==0`
  even when prepared for maxBlock=8192 and fed a smaller block — it is NOT merely "block-granular". Its *cost*,
  not its latency, tracks the actual block: prepared for 8192 yet fed 256, it stays 0-latency but pays 2.82%.)
- **The worst-buffer spike** (all stages' FFTs coincide every `lcm=B_max` samples) is the RT trade-off: bounded,
  worst at tiny blocks (7.8% @block 64 for `B_max=4096`). A `B_max` sweep shows the *mean* is flat ~0.8% for any
  `B_max`, but the spike grows with `B_max` (block 64: `B_max=2048`→6.4%, `4096`→7.9%, `8192`→12%) — so
  **`B_max=2048` is the small-block default** (same mean, smaller spike). Distributing the large FFTs over their
  deadline (time-distributed scheduling) is a later-phase RT-hardening.

**Correctness:** NULL-verified sample-by-sample vs a direct time-domain convolution AND the proven
`PartitionedConvolver` (rel ~2e-7); true zero latency proven by an impulse (in → IR out, no shift); a pffft-vs-
scalar cross-null on the shipping backend; no-alloc in `process()`; schedule/coverage validation guards.

**Acceptance (bar "A-sharpened", Oleh):** beat JUCE at small blocks + block-INDEPENDENT + true sample-zero-
latency + ≥2× better than v0.4.0 — MET. Beating JUCE at large oracle-tuned blocks is out of scope (the theorem
above). Consumers keep true zero latency vs the block-granular latency a large-`P` uniform path would impose.

## Reproduce

```
cmake -S . -B build-juce -DCMAKE_BUILD_TYPE=Release -DFELITRONICS_WITH_PFFFT=ON -DFELITRONICS_BENCH_JUCE=ON
cmake --build build-juce --target fcore_fftbench -j
./build-juce/tools/fcore_fftbench          # prints the correctness probe, the OLD/NEW/pffft + NUPC sweep, and the JUCE head-to-head
```

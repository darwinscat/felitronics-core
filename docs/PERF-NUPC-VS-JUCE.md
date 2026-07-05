<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

# Convolver benchmark — MatrixConvolverNupc vs juce::dsp::Convolution

The `lineareq` EQ convolves on **`MatrixConvolverNupc`** — a JUCE-free, non-uniform (Gardner 1995) partitioned
convolver: a 128-sample time-domain head + geometrically growing overlap-save FFT tail stages, on the SIMD
`pffft` backend. It is **block-INDEPENDENT** (one `prepare()`, the same cost at every host block), **true
sample-zero-latency**, and a NULL-verified drop-in for the old fixed-`P=128` `MatrixConvolver`.

**📊 [Interactive chart →](https://claude.ai/code/artifact/1a118004-7d2f-4dc5-8522-bd95d743e4d9)**

## Head-to-head (131072-tap stereo IR, LRDiag, 48 kHz, Apple M5 Pro, pffft · JUCE 8.0.4)

All three report `getLatency() == 0`. `%RT = CPU / audio time`; lower is better.

| host block | `juce::dsp::Convolution` | `MatrixConvolver<pffft>` (old) | **`MatrixConvolverNupc<pffft>`** | ratio |
|---:|---:|---:|---:|---|
| 64   | 10.50% | 2.17% | **0.62%** | **17× vs JUCE · 3.5× vs old** |
| 128  | 3.40%  | 2.23% | **0.62%** | 5.5× vs JUCE · 3.6× vs old |
| 256  | 2.00%  | 2.30% | **0.67%** | 3.0× vs JUCE · 3.4× vs old |
| 512  | 0.98%  | 2.21% | **0.64%** | 1.5× vs JUCE · 3.5× vs old |
| 1024 | 0.50%  | 2.19% | **0.62%** | JUCE 1.2× · 3.5× vs old |
| 2048 | 0.29%  | 2.24% | **0.63%** | JUCE 2.2× · 3.6× vs old |
| 4096 | 0.19%  | 2.29% | **0.65%** | JUCE 3.4× · 3.5× vs old |
| 8192 | 0.15%  | 2.22% | **0.64%** | JUCE 4.3× · 3.5× vs old |

**The result.** `MatrixConvolverNupc` is a flat **~0.62 % RT at every host block** — **≈3.5× cheaper than the old
fixed-`P=128` matrix convolver everywhere**, and **3–17× cheaper than JUCE at the 64–128-sample buffers a guitar
amp or live monitor actually runs**. JUCE's cost swings ~70× with the block (10.50 % → 0.15 %): it is expensive
at the small blocks and only overtakes us past ~512 samples, where its few large partitions cost less.

## Every real DAW buffer — the nose-to-nose sweep (all 65 sizes)

The head-to-head above walks powers of two. But hosts offer **every 32-sample buffer** (16, 32, 64, 96, 128,
160, …, 2048) — and most of the ones users pick are **not powers of two**. Measured at all 65, LRDiag, JUCE 8.0.4
oracle-tuned (`maxBlock` = the actual buffer, its best case). NUPC is **dead-flat 0.58–0.67 % (mean ~0.6 %)** at
every one; JUCE is a sawtooth that peaks at each power of two. Representative rows (`*` = non-power-of-two; JUCE's
tiny-block cost carries a few-% run-to-run variance, so this sweep and the powers-of-two table above are independent
runs that differ slightly on JUCE — e.g. 9.74 % vs 10.50 % at block 64):

| host buffer | ms | `juce::dsp::Convolution` 8.0.4 | **NUPC mean** | NUPC worst buffer | cheaper |
|---:|---:|---:|---:|---:|---|
| 16    | 0.33 | **102.6 %** ⛔ (over real-time) | **0.59 %** | 52.5 % | **174× us** |
| 32    | 0.67 | 32.1 % | **0.59 %** | 15.2 % | **54× us** |
| 64    | 1.33 | 9.74 % | **0.60 %** | 10.6 % | **16× us** |
| 96\*  | 2.00 | 3.05 % | **0.59 %** | 4.97 % | 5.2× us |
| 128   | 2.67 | 2.84 % | **0.60 %** | 3.64 % | 4.7× us |
| 192\* | 4.00 | 1.99 % | **0.62 %** | 2.11 % | 3.2× us |
| 256   | 5.33 | 1.94 % | **0.59 %** | 1.83 % | 3.3× us |
| 512   | 10.67 | 0.95 % | **0.61 %** | 1.51 % | 1.6× us |
| 992\* | 20.67 | 0.57 % | **0.59 %** | 0.90 % | ≈ even |
| 1024  | 21.33 | 0.49 % | **0.59 %** | 1.48 % | JUCE 1.2× |
| 2048  | 42.67 | 0.32 % | **0.62 %** | 0.70 % | JUCE 1.9× |

**Reading it.** JUCE exceeds real-time only at the extreme 16-sample buffer (102.6 %), but stays **2–170× more
expensive than NUPC up to ~512 samples** — the whole low-latency / live-rig range. The two lines cross near **544
samples**: below it NUPC's flat mean wins, above it JUCE's few large partitions win the mean (see *Honest bounds*).
Both are zero-latency throughout. The `NUPC worst buffer` column is the single most expensive buffer in the window
(the once-per-`B_max` coincident-FFT): it spikes at tiny blocks (52 % at 16) but stays **under 100 % — no xrun**,
while JUCE's *mean* there already isn't.

> **Measurement note (reproducibility).** The sweep's warm-up **must outlast NUPC's cold-prime crossfade**
> (`coldXfade_ = max_s(C_s·B_s)` ≈ 2.5 s for `B_max = 2048` over a 131072-tap IR). A shorter warm leaves the
> measure window still double-convolving the warm + cold slots and **inflates the mean by ~0.35 %** — the sweep
> uses a 3.0 s warm so it settles to the same steady state as the head-to-head.

## Per topology (all flat, block-independent)

`MatrixConvolverNupc<pffft>`, mean %RT (worst single-buffer %RT in parentheses):

| host block | LRDiag | MSDiag | Full |
|---:|---:|---:|---:|
| 64   | 0.63% (24.8%) | 0.76% (8.8%) | 0.81% (11.8%) |
| 256  | 0.61% (1.7%)  | 0.75% (3.0%) | 0.81% (2.9%)  |
| 1024 | 0.62% (0.8%)  | 0.74% (1.1%) | 0.82% (1.2%)  |
| 8192 | 0.62% (0.6%)  | 0.75% (0.8%) | 0.83% (0.9%)  |

MSDiag / Full add the ½(X_L±X_R) view / the 4-bank cross sums — a modest, flat surcharge over LRDiag.

## Honest bounds

- **JUCE wins the mean at large blocks (≥1024).** That is a theorem, not a shortfall: true *sample*-zero-latency
  forces small partitions in the near field (FFT overhead) that JUCE never pays. The win is small-blocks +
  block-independence; both are zero-latency. Chasing JUCE at 8192 would mean *becoming* JUCE (block-granular,
  `P=maxBlock`) — surrendering the thing it structurally cannot do.
- **The worst-buffer spike.** Every stage's FFT coincides once per `lcm = B_max` samples; at tiny blocks that one
  buffer is expensive (LRDiag ~25 % of a single 64-sample buffer, dropping to ~1 % by block 512). The *mean* is
  flat ~0.62 %; the spike is bounded and < 100 % (no xrun on its own). Distributing the large FFTs over their
  deadline (time-distributed scheduling) is a later-phase RT-hardening; a smaller `B_max` already trims it.
- **Variable buffers — the real-world edge.** JUCE's cost tracks the *actual* block: prepared for maxBlock=8192
  but fed 256, it stays zero-latency yet pays 2.67 % (worse than its own 256-tuned 2.00 %). NUPC is flat ~0.62 %
  from one prepare — cheaper and steadier for a host that delivers small or variable buffers.

## Reproduce

```sh
cmake -S . -B build-juce -DCMAKE_BUILD_TYPE=Release -DFELITRONICS_WITH_PFFFT=ON -DFELITRONICS_BENCH_JUCE=ON
cmake --build build-juce --target fcore_fftbench -j
./build-juce/tools/fcore_fftbench                       # correctness probe + OLD/NEW/pffft + NUPC + matrix-NUPC + JUCE head-to-head
FCORE_SWEEP_ONLY=1 ./build-juce/tools/fcore_fftbench    # only the 65-buffer nose-to-nose sweep, as CSV (buffer,ms,juce,nupc,nupc_max)
```

Correctness: `MatrixConvolverNupc`'s output is NULL-verified sample-by-sample against a direct time-domain
convolution and the proven `PartitionedConvolver`, across all topologies + the click-free swap, on both the
scalar and pffft backends (a randomized differential fuzzer + ASan/UBSan/TSan). See
[`PERF-CONVOLVER-JUCE-GAP.md`](PERF-CONVOLVER-JUCE-GAP.md) for the design (the non-uniform schedule, the
zero-latency invariant, the phased build).

<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

# Convolver benchmark — MatrixConvolverNupc vs juce::dsp::Convolution

The `lineareq` EQ convolves on **`MatrixConvolverNupc`** — a JUCE-free, non-uniform (Gardner 1995) partitioned
convolver: a 128-sample time-domain head + geometrically growing overlap-save FFT tail stages, on the SIMD
`pffft` backend. It is **block-INDEPENDENT** (one `prepare()`, the same cost at every host block), **true
sample-zero-latency**, and a NULL-verified drop-in for the old fixed-`P=128` `MatrixConvolver`.

**📊 [Interactive chart →](https://claude.ai/code/artifact/1a118004-7d2f-4dc5-8522-bd95d743e4d9)**

## Head-to-head (131072-tap stereo IR, LRDiag, 48 kHz, Apple M5 Pro, pffft)

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
./build-juce/tools/fcore_fftbench   # the correctness probe, the OLD/NEW/pffft + NUPC + matrix-NUPC sweep, and the JUCE head-to-head
```

Correctness: `MatrixConvolverNupc`'s output is NULL-verified sample-by-sample against a direct time-domain
convolution and the proven `PartitionedConvolver`, across all topologies + the click-free swap, on both the
scalar and pffft backends (a randomized differential fuzzer + ASan/UBSan/TSan). See
[`PERF-CONVOLVER-JUCE-GAP.md`](PERF-CONVOLVER-JUCE-GAP.md) for the design (the non-uniform schedule, the
zero-latency invariant, the phased build).
